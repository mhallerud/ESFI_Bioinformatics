#!/usr/bin/env Rscript
##
## dada2_filter.R -- stage 2 of the DADA2 branch: quality filtering.
##
##   Rscript dada2_filter.R --in-dir <demux dir> --out-dir <filtered dir> \
##                          --mode single|paired [--pattern '_LOCUS_trnL\.'] ...
##
## Runs one batch (a library x locus slice, or a whole locus) in a single R
## session. filterAndTrim is itself multithreaded, so one array task per batch
## keeps the scheduler busy without paying R's startup cost 2000 times.
##
## Truncated to 200 bp for BF3-BR2 since trimmed amplicons should be 418 bp, 
## plus handles quality issues with R2
##
## Writes <out-dir>/ESFI_<lib>/<SampleID>[.R1|.R2].fastq.gz and a filter_stats.tsv.

suppressPackageStartupMessages({
  args0 <- commandArgs(trailingOnly = FALSE)
  here <- dirname(sub("^--file=", "", args0[grep("^--file=", args0)]))
  if (!length(here) || !nzchar(here)) here <- "."
  source(file.path(here, "./dada2_common.R"))
})

## IMPORTANT: Filtering Settings
args <- parse_args(list(
  in_dir = NULL, out_dir = NULL, mode = "paired", pattern = "",
  maxEE = "c(2,5)", truncQ = "2", maxN = "50", minLen = "100", maxLen = "Inf",
  trimLeft = "0", rm_phix = "TRUE", threads = "10", stats = ""
))

if (is.null(args$in_dir) || is.null(args$out_dir))
  stop("--in-dir and --out-dir are required")

mode    <- match.arg(args$mode, c("single", "paired"))
threads <- arg_threads(args$threads)
maxEE   <- as.numeric(strsplit(args$maxEE, ",")[[1]])
if (mode == "paired" && length(maxEE) == 1L) maxEE <- c(maxEE, maxEE)

dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)
stats_path <- if (nzchar(args$stats)) args$stats else
  file.path(args$out_dir, "filter_stats.tsv")
dir.create(dirname(stats_path), recursive = TRUE, showWarnings = FALSE)

## Mirror the ESFI_<lib>/ layout of the input under the output directory.
out_path <- function(in_paths) {
  rel <- sub(paste0("^", normalizePath(args$in_dir, mustWork = FALSE), "/?"), "",
             normalizePath(in_paths, mustWork = FALSE))
  p <- file.path(args$out_dir, rel)
  p <- sub("\\.fastq$|\\.fq$", ".fastq.gz", p)
  p <- sub("\\.fastq\\.gz$|\\.fq\\.gz$", ".fastq.gz", p)
  for (d in unique(dirname(p))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  p
}

if (mode == "single") {
  fns <- list_reads(args$in_dir, args$pattern, "single")
  if (!length(fns)) stop("no input FASTQ found under ", args$in_dir)
  msg("filtering ", length(fns), " single-end files")
  filts <- out_path(fns)
  out <- filterAndTrim(fns, filts,
                       maxEE = maxEE[1], truncQ = arg_num(args, "truncQ"),
                       maxN = arg_num(args, "maxN"),
                       minLen = arg_num(args, "minLen"),
                       maxLen = arg_num(args, "maxLen"),
                       trimLeft = arg_num(args, "trimLeft"),
                       rm.phix = arg_lgl(args, "rm_phix"),
                       compress = TRUE, verbose = FALSE,
                       multithread = threads)
  df <- data.frame(SampleID = sample_id_of(fns),
                   reads_in = out[, 1], reads_filtered = out[, 2],
                   stringsAsFactors = FALSE)
} else {
  files <- list.files(args$in_dir, full.names=TRUE)
  fnFs <- files[grepl("R1.fastq", files)] #list_reads(args$in_dir, args$pattern, "R1")
  fnRs <- files[grepl("R2.fastq", files)] #list_reads(args$in_dir, args$pattern, "R2")
  if (!length(fnFs)) stop("no R1 FASTQ found under ", args$in_dir)
  ## Keep only samples that have both mates.
  #common <- intersect(sample_id_of(fnFs), sample_id_of(fnRs))
  #fnFs <- fnFs[sample_id_of(fnFs) %in% common]
  #fnRs <- fnRs[sample_id_of(fnRs) %in% common]
  #stopifnot(identical(sample_id_of(fnFs), sample_id_of(fnRs)))
  msg("filtering ", length(fnFs), " read pairs")
  print(fnFs)
  print(out_path(fnFs))
  out <- filterAndTrim(fnFs, out_path(fnFs), fnRs, out_path(fnRs),
                       maxEE = c(2,5), truncQ = arg_num(args, "truncQ"),
                       maxN = arg_num(args, "maxN"),
                       minLen = arg_num(args, "minLen"),
                       maxLen = arg_num(args, "maxLen"),
                       trimLeft = arg_num(args, "trimLeft"),
                       rm.phix = arg_lgl(args, "rm_phix"),
                       matchIDs = TRUE,
                       compress = TRUE, verbose = FALSE,
                       multithread = threads)
  df <- data.frame(SampleID = sample_id_of(fnFs),
                   reads_in = out[, 1], reads_filtered = out[, 2],
                   stringsAsFactors = FALSE)
}

write.table(df, stats_path, sep = "\t", quote = FALSE, row.names = FALSE)
msg("kept ", sum(df$reads_filtered), " / ", sum(df$reads_in), " reads (",
    round(100 * sum(df$reads_filtered) / max(1, sum(df$reads_in)), 1), "%)")
msg("stats -> ", stats_path)
