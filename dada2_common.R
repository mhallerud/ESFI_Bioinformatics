## dada2_common.R -- shared helpers for the optional DADA2 branch.
## source()d by dada2_filter.R / dada2_learn_errors.R / dada2_denoise.R /
## dada2_combine.R. Not meant to be run on its own.

suppressPackageStartupMessages(library(dada2))

## ---------------------------------------------------------------- arguments
## Accepts --key=value and --key value. Returns a named character list with
## `defaults` filled in for anything not supplied.
parse_args <- function(defaults = list()) {
  raw <- commandArgs(trailingOnly = TRUE)
  out <- defaults
  i <- 1L
  while (i <= length(raw)) {
    a <- raw[[i]]
    if (!grepl("^--", a)) stop("unexpected argument: ", a)
    a <- sub("^--", "", a)
    if (grepl("=", a, fixed = TRUE)) {
      k <- sub("=.*$", "", a); v <- sub("^[^=]*=", "", a)
      i <- i + 1L
    } else {
      k <- a
      if (i + 1L <= length(raw) && !grepl("^--", raw[[i + 1L]])) {
        v <- raw[[i + 1L]]; i <- i + 2L
      } else { v <- "TRUE"; i <- i + 1L }
    }
    out[[gsub("-", "_", k)]] <- v
  }
  out
}

arg_num <- function(args, key, default = NA_real_) {
  v <- args[[key]]
  if (is.null(v) || is.na(v) || identical(v, "")) return(default)
  if (identical(toupper(v), "INF")) return(Inf)
  as.numeric(v)
}

arg_lgl <- function(args, key, default = FALSE) {
  v <- args[[key]]
  if (is.null(v)) return(default)
  toupper(v) %in% c("TRUE", "T", "YES", "1")
}

## dada(pool=) takes FALSE, TRUE or "pseudo".
arg_pool <- function(v) {
  if (is.null(v)) return(FALSE)
  if (toupper(v) == "FALSE") return(FALSE)
  if (toupper(v) == "TRUE")  return(TRUE)
  "pseudo"
}

## multithread= takes TRUE or an integer core count.
arg_threads <- function(v) {
  n <- suppressWarnings(as.integer(v))
  if (is.na(n) || n < 1L) return(TRUE)
  n
}

msg <- function(...) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ..., "\n", sep = "")
  flush.console()
}

## ------------------------------------------------------- file / name helpers
## Demultiplexed file names are "<SampleID>.fastq[.gz]" (merged mode) or
## "<SampleID>.R1.fastq[.gz]" / ".R2." (paired mode). SampleID keeps the exact
## grammar demult.sh produces, e.g.
##   TNFM096_c_LIB_A_FBAR_ACAGTGAC_RBAR_ACAGTGAC_LOCUS_trnL
## LocalTaxaTool splits that on "_" and demands exactly 10 fields, so it must
## survive untouched all the way to the .clust.fasta name.
sample_id_of <- function(paths) {
  b <- basename(paths)
  b <- sub("\\.gz$", "", b)
  b <- sub("\\.fastq$", "", b)
  b <- sub("\\.fq$", "", b)
  sub("\\.R[12]$", "", b)
}

lib_of <- function(sample_ids) {
  bad <- !grepl("_LIB_[^_]+_", sample_ids)
  if (any(bad)) stop("could not parse _LIB_ from sample id: ",
                     paste(head(sample_ids[bad], 3), collapse = ", "))
  sub("^.*_LIB_([^_]+)_.*$", "\\1", sample_ids)
}

locus_of <- function(sample_ids) sub("^.*_LOCUS_", "", sample_ids)

## List the R1 (or single-end) files of a batch directory, sorted for
## reproducibility.
list_reads <- function(dir, pattern = NULL, mate = c("single", "R1", "R2")) {
  mate <- match.arg(mate)
  pat <- switch(mate,
                single = "\\.fastq(\\.gz)?$|\\.fq(\\.gz)?$",
                R1     = "\\.R1\\.fastq(\\.gz)?$|\\.R1\\.fq(\\.gz)?$",
                R2     = "\\.R2\\.fastq(\\.gz)?$|\\.R2\\.fq(\\.gz)?$")
  f <- list.files(dir, pattern = pat, full.names = TRUE, recursive = TRUE)
  if (mate == "single") f <- f[!grepl("\\.R[12]\\.f(ast)?q(\\.gz)?$", f)]
  if (!is.null(pattern) && nzchar(pattern)) f <- f[grepl(pattern, basename(f))]
  sort(f)
}

## ------------------------------------------------------------- error models
## Binned-quality loess fit.
##
## This flowcell reports 4 quality values (2, 12, 24, 40) instead of a
## continuous scale, which is a known failure mode for dada2's stock
## loessErrfun: it fits a degree-2 loess with span 0.75 through what is
## effectively 4 points, and the result is frequently non-monotonic -- error
## rates that *rise* with quality -- which makes dada() split real variants or
## merge distinct ones depending on which way the fit bends.
##
## The fix, same one used for NovaSeq/NextSeq data elsewhere: weight by
## log10(total observations), fit degree 1 with a wide span so the 4 points
## define a line, then force the fitted rate to be non-increasing in Q. The
## monotonic step raises low-Q rates to the maximum seen at or above them, which
## is the conservative direction -- it denoises less, not more.
loessErrfun_binned <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow = 0, ncol = length(qq))
  for (nti in c("A", "C", "G", "T")) {
    for (ntj in c("A", "C", "G", "T")) {
      if (nti != ntj) {
        errs <- trans[paste0(nti, "2", ntj), ]
        tot  <- colSums(trans[paste0(nti, "2", c("A", "C", "G", "T")), ])
        rlogp <- log10((errs + 1) / tot)
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q = qq, errs = errs, tot = tot, rlogp = rlogp)
        mod.lo <- suppressWarnings(
          loess(rlogp ~ q, df, weights = log10(tot + 1), span = 2, degree = 1)
        )
        pred <- suppressWarnings(predict(mod.lo, qq))
        if (all(is.na(pred))) {
          ## Degenerate transition (no observations at all): fall back to a
          ## flat, deliberately pessimistic rate rather than dropping the row.
          pred <- rep(log10(1e-4), length(qq))
        }
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred) > maxrli] <- pred[[maxrli]]
        pred[seq_along(pred) < minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      }
    }
  }

  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est > MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est < MIN_ERROR_RATE] <- MIN_ERROR_RATE

  ## Force each transition's error rate to be non-increasing in Q.
  est <- t(apply(est, 1, function(x) rev(cummax(rev(x)))))
  est[est > MAX_ERROR_RATE] <- MAX_ERROR_RATE

  err <- rbind(1 - colSums(est[1:3, ]), est[1:3, ],
               est[4, ], 1 - colSums(est[4:6, ]), est[5:6, ],
               est[7:8, ], 1 - colSums(est[7:9, ]), est[9, ],
               est[10:12, ], 1 - colSums(est[10:12, ]))
  rownames(err) <- paste0(rep(c("A", "C", "G", "T"), each = 4), "2",
                          c("A", "C", "G", "T"))
  colnames(err) <- colnames(trans)
  err
}

pick_errfun <- function(name) {
  if (identical(name, "binned")) loessErrfun_binned else dada2:::loessErrfun
}

## --------------------------------------------------------------- fasta out
## Write one sample's ASVs in exactly the format fastx_collapser produces, so
## LocalTaxaTool (and blast.sh, and every downstream script) needs no change:
##   >rank-count
##   SEQUENCE
write_clust_fasta <- function(seqs, counts, path) {
  keep <- counts > 0
  seqs <- seqs[keep]; counts <- counts[keep]
  o <- order(counts, decreasing = TRUE)
  seqs <- seqs[o]; counts <- counts[o]
  if (!length(seqs)) { file.create(path); return(invisible(0L)) }
  lines <- character(2L * length(seqs))
  lines[c(TRUE, FALSE)] <- paste0(">", seq_along(seqs), "-", counts)
  lines[c(FALSE, TRUE)] <- seqs
  writeLines(lines, path)
  invisible(length(seqs))
}

## Write a whole sequence table out as per-sample .clust.fasta files under
## <out_dir>/LIB_<lib>/<SampleID>.clust.fasta
write_seqtab_as_clust <- function(seqtab, out_dir) {
  seqs <- colnames(seqtab)
  ids  <- rownames(seqtab)
  libs <- lib_of(ids)
  for (l in unique(libs)) dir.create(file.path(out_dir, paste0("LIB_", l)),
                                     recursive = TRUE, showWarnings = FALSE)
  n <- integer(length(ids))
  for (i in seq_along(ids)) {
    p <- file.path(out_dir, paste0("LIB_", libs[[i]]),
                   paste0(ids[[i]], ".clust.fasta"))
    n[[i]] <- write_clust_fasta(seqs, as.integer(seqtab[i, ]), p)
  }
  invisible(setNames(n, ids))
}
