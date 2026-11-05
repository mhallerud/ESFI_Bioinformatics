#!/usr/bin/env Rscript
##
## dada2_denoise.R -- stage 4 of the DADA2 branch: the denoising itself.
##
##   Rscript dada2_denoise.R --in-dir <filtered dir> --out-dir <batch dir> \
##                           --mode single|paired [--err <err.rds>] \
##                           [--fasta-dir <run>/dada2] [--pattern '_LOCUS_trnL\.']
##
## One invocation = one batch (library x locus by default). If --err is not
## given the error model is learned inside the batch, which keeps every task
## self-contained and lets the scheduler run all batches at once with no
## barrier. Pass --err when you want a model shared across libraries.
##
## Outputs, all inside --out-dir:
##   seqtab.rds        ASV table before chimera removal  (samples x ASVs)
##   seqtab_nochim.rds ASV table after chimera removal
##   asv_table.tsv     the same, as text
##   asvs.fasta        the ASV sequences, >ASV_n
##   track.tsv         reads surviving each step, per sample
##   err.rds / err.pdf the in-batch error model, when one was learned here
##
## And, unless --no-fasta, per-sample drop-ins for the existing downstream code:
##   <fasta-dir>/LIB_<lib>/<SampleID>.clust.fasta   in fastx_collapser's
##                                                  ">rank-count" format
##
## What this replaces: `fastx_collapser` at the end of demult.sh. That step
## dereplicates; this one dereplicates *and* infers which unique sequences are
## real and which are sequencing error, then folds the errors into their parent.

suppressPackageStartupMessages({
  args0 <- commandArgs(trailingOnly = FALSE)
  here <- dirname(sub("^--file=", "", args0[grep("^--file=", args0)]))
  if (!length(here) || !nzchar(here)) here <- "."
  source(file.path(here, "../dada2/dada2_common.R"))
})

args <- parse_args(list(
  in_dir = NULL, out_dir = NULL, mode = "single", pattern = "",
  err = "", fasta_dir = "", no_fasta = "FALSE",
  pool = "pseudo", omega_a = "1e-40", threads = "8",
  minOverlap = "10", maxMismatch = "0",
  chimera = "consensus", errfun = "binned", nbases = "1e8",
  filter_stats = "", min_reads = "1"
))
if (is.null(args$in_dir) || is.null(args$out_dir))
  stop("--in-dir and --out-dir are required")

mode    <- match.arg(args$mode, c("single", "paired"))
threads <- arg_threads(args$threads)
pool    <- arg_pool(args$pool)
omegaA  <- arg_num(args, "omega_a", 1e-40)
min_reads <- as.integer(arg_num(args, "min_reads", 1))
dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

getN <- function(x) sum(getUniques(x))

## ------------------------------------------------------------ input files
if (mode == "single") {
  filtFs <- list_reads(args$in_dir, args$pattern, "single")
  filtRs <- character(0)
  if (!length(filtFs)) stop("no filtered FASTQ found under ", args$in_dir)
} else {
  filtFs <- list_reads(args$in_dir, args$pattern, "R1")
  ok <- file.exists(sub("\\.R1\\.", ".R2.", filtFs))
  filtFs <- filtFs[ok]
  filtRs <- sub("\\.R1\\.", ".R2.", filtFs)
  if (!length(filtFs)) stop("no filtered R1/R2 pairs found under ", args$in_dir)
}

## Samples with almost no reads make dada() unhappy and contribute nothing;
## drop them here and record them in track.tsv as zero.
nreads <- vapply(filtFs, function(f) {
  tryCatch(as.integer(ShortRead::countFastq(f)$records), error = function(e) 0L)
}, integer(1))
names(nreads) <- filtFs
dropped <- filtFs[nreads < min_reads]
keep <- nreads >= min_reads
filtFs <- filtFs[keep]
if (mode == "paired") filtRs <- filtRs[keep]
if (!length(filtFs)) stop("every sample in this batch is below --min-reads")

names(filtFs) <- sample_id_of(filtFs)
if (mode == "paired") names(filtRs) <- sample_id_of(filtRs)
msg("batch: ", length(filtFs), " samples, ", sum(nreads[keep]), " filtered reads",
    if (length(dropped)) paste0(" (", length(dropped), " empty samples skipped)") else "")

## ------------------------------------------------------------ error model
if (nzchar(args$err)) {
  e <- readRDS(args$err)
  errF <- e$errF; errR <- e$errR
  msg("using error model from ", args$err)
} else {
  efun <- pick_errfun(args$errfun)
  nb   <- arg_num(args, "nbases", 1e8)
  msg("learning in-batch error model (", args$errfun, ")")
  errF <- learnErrors(filtFs, nbases = nb, errorEstimationFunction = efun,
                      multithread = threads, randomize = TRUE, verbose = 1)
  errR <- if (mode == "paired")
    learnErrors(filtRs, nbases = nb, errorEstimationFunction = efun,
                multithread = threads, randomize = TRUE, verbose = 1) else NULL
  saveRDS(list(mode = mode, errF = errF, errR = errR, errfun = args$errfun),
          file.path(args$out_dir, "err.rds"))
  pdf(file.path(args$out_dir, "err.pdf"), width = 8, height = 8)
  print(plotErrors(errF, nominalQ = TRUE))
  if (!is.null(errR)) print(plotErrors(errR, nominalQ = TRUE))
  invisible(dev.off())
}

## ------------------------------------------------------------- denoising
msg("dada() on forward/single reads, pool=", format(pool))
ddF <- dada(filtFs, err = errF, pool = pool, multithread = threads,
            OMEGA_A = omegaA, verbose = FALSE)
if (length(filtFs) == 1L) ddF <- list(ddF); names(ddF) <- names(filtFs)

if (mode == "paired") {
  msg("dada() on reverse reads")
  ddR <- dada(filtRs, err = errR, pool = pool, multithread = threads,
              OMEGA_A = omegaA, verbose = FALSE)
  if (length(filtRs) == 1L) ddR <- list(ddR); names(ddR) <- names(filtRs)

  msg("mergePairs(minOverlap=", args$minOverlap, ", maxMismatch=",
      args$maxMismatch, ")")
  mergers <- mergePairs(ddF, filtFs, ddR, filtRs,
                        minOverlap = arg_num(args, "minOverlap", 10),
                        maxMismatch = arg_num(args, "maxMismatch", 0),
                        verbose = FALSE)
  if (!is.list(mergers) || is.data.frame(mergers)) {
    mergers <- list(mergers); names(mergers) <- names(filtFs)
  }
  seqtab <- makeSequenceTable(mergers)
} else {
  seqtab <- makeSequenceTable(ddF)
}
msg("sequence table: ", nrow(seqtab), " samples x ", ncol(seqtab), " ASVs")

## An empty table is a configuration problem, not a crash: it means nothing
## survived filtering/merging. Say so, write the empty outputs so the sample
## list downstream still matches, and stop cleanly.
if (ncol(seqtab) == 0L) {
  msg("WARNING: no ASVs in this batch. Check track.tsv, the primer sequences, ",
      "and --minOverlap.")
  saveRDS(seqtab, file.path(args$out_dir, "seqtab.rds"))
  saveRDS(seqtab, file.path(args$out_dir, "seqtab_nochim.rds"))
  if (nzchar(args$fasta_dir)) {
    for (sid in unique(c(sample_id_of(filtFs), sample_id_of(dropped)))) {
      lb <- lib_of(sid)
      dir.create(file.path(args$fasta_dir, paste0("LIB_", lb)),
                 recursive = TRUE, showWarnings = FALSE)
      file.create(file.path(args$fasta_dir, paste0("LIB_", lb),
                            paste0(sid, ".clust.fasta")))
    }
  }
  quit(save = "no", status = 0)
}

## ------------------------------------------------------------- chimeras
if (identical(args$chimera, "none")) {
  seqtab_nc <- seqtab
  msg("chimera removal skipped (--chimera none)")
} else {
  seqtab_nc <- removeBimeraDenovo(seqtab, method = args$chimera,
                                  multithread = threads, verbose = FALSE)
  msg("after ", args$chimera, " chimera removal: ", ncol(seqtab_nc), " ASVs (",
      round(100 * sum(seqtab_nc) / max(1, sum(seqtab)), 1), "% of reads kept)")
}

## ---------------------------------------------------------------- outputs
saveRDS(seqtab,    file.path(args$out_dir, "seqtab.rds"))
saveRDS(seqtab_nc, file.path(args$out_dir, "seqtab_nochim.rds"))

asv_names <- paste0("ASV_", seq_len(ncol(seqtab_nc)))
writeLines(as.vector(rbind(paste0(">", asv_names), colnames(seqtab_nc))),
           file.path(args$out_dir, "asvs.fasta"))
tab <- data.frame(SampleID = rownames(seqtab_nc), seqtab_nc,
                  check.names = FALSE, stringsAsFactors = FALSE)
colnames(tab)[-1] <- asv_names
write.table(tab, file.path(args$out_dir, "asv_table.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## Reads surviving each step. Always read this before trusting the output --
## a big drop at mergePairs means minOverlap is above the shortest amplicon,
## a big drop at chimera removal on a short locus usually means real variants
## are being called chimeric.
ids <- rownames(seqtab_nc)
track <- data.frame(SampleID = ids, stringsAsFactors = FALSE)
if (nzchar(args$filter_stats) && file.exists(args$filter_stats)) {
  fs <- read.delim(args$filter_stats, stringsAsFactors = FALSE)
  track <- merge(track, fs, by = "SampleID", all.x = TRUE, sort = FALSE)
  track <- track[match(ids, track$SampleID), ]
}
track$denoisedF <- vapply(ddF[ids], getN, numeric(1))
if (mode == "paired") {
  track$denoisedR <- vapply(ddR[ids], getN, numeric(1))
  track$merged    <- vapply(mergers[ids], getN, numeric(1))
}
track$nonchim <- rowSums(seqtab_nc)[ids]
if (length(dropped)) {
  empty <- data.frame(SampleID = sample_id_of(dropped))
  for (cn in setdiff(colnames(track), "SampleID")) empty[[cn]] <- 0
  track <- rbind(track, empty[, colnames(track), drop = FALSE])
}
write.table(track, file.path(args$out_dir, "track.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

if (nzchar(args$fasta_dir)) {
  if (!arg_lgl(args, "no_fasta")) {
    n <- write_seqtab_as_clust(seqtab_nc, args$fasta_dir)
    msg("wrote ", length(n), " .clust.fasta files under ", args$fasta_dir)
  }
  ## Every sample the batch started with needs a file, even if nothing survived
  ## -- demult.sh always leaves one behind, and LocalTaxaTool reports zero-size
  ## files rather than silently shrinking the sample list.
  expected <- unique(c(sample_id_of(filtFs), sample_id_of(dropped),
                       if (nzchar(args$filter_stats) && file.exists(args$filter_stats))
                         read.delim(args$filter_stats, stringsAsFactors = FALSE)$SampleID))
  missing <- setdiff(expected, rownames(seqtab_nc))
  for (sid in missing) {
    lb <- lib_of(sid)
    dir.create(file.path(args$fasta_dir, paste0("LIB_", lb)),
               recursive = TRUE, showWarnings = FALSE)
    f <- file.path(args$fasta_dir, paste0("LIB_", lb), paste0(sid, ".clust.fasta"))
    if (!file.exists(f)) file.create(f)
  }
  if (length(missing))
    msg(length(missing), " sample(s) produced no ASVs; empty .clust.fasta written")
}
msg("batch done -> ", args$out_dir)
