# set up environment
.libPaths(.libPaths()[-1])
library(dada2)
library(ShortRead)
library(Biostrings)

# set input variables
## ARGUMENTS: 1 = path to library (e.g., LIB_ITSP1); 2 = prefix for outputs (e.g., ITSP1_cutadapt) 
run_dada2 <- function(lib.path, prefix){
#args <- commandArgs(trailingOnly = TRUE)
#lib.path=args[1]
#prefix=args[2]
path.cut=file.path(lib.path, "cutadapt")

currdir=getwd()
#setwd(file.path(currdir, path.cut))

# set up output directory
outdir=file.path(lib.path, "dada2")
if (!dir.exists(outdir)){
	dir.create(outdir)
}#if

# make list of input reads
cutFs <- list.files(path.cut, "R1_trim.fastq")
cutRs <- list.files(path.cut, "R2_trim.fastq")
cutFs <- file.path(path.cut, cutFs)
cutRs <- file.path(path.cut, cutRs)

# plot initial quality profiles
pdf(file.path(outdir, paste0(prefix,"_InputQualityProfile_FWD.pdf")))
plotQualityProfile(cutFs[1:10])
dev.off()

pdf(file.path(outdir, paste0(prefix,"_InputQualityProfile_REV.pdf")))
plotQualityProfile(cutRs[1:10])
dev.off()

# make list for output filenames
filtFs <- file.path(outdir, "filtered", basename(cutFs))
filtRs <- file.path(outdir, "filtered", basename(cutRs))

# quality filtering of reads- we’re going to use some looser quality filtering than the maxEE=c(2,2) in the tutorial, esp. on the R reads
# maxEE set to (4,4) for 3 allowed errors per sequence. Given that we're clustering to 97% similarity, this would be 3% error for even short sequences 
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN=0, maxEE = c(3, 5), rm.phix = TRUE, compress = TRUE, multithread = TRUE)
#write.csv(out, file.path(outdir, 'filtered', paste0(prefix,'_filterstats.csv')))

# list output filtered files (this might not be the same as the original, for example if no reads in blank samples didn’t clear filtering)
filtFs <- list.files(file.path(outdir,"filtered"), "_R1_trim.fastq")
filtRs <- list.files(file.path(outdir,"filtered"), "_R2_trim.fastq")

# learn error rates separately for FW and REV reads
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# plot error rate info
pdf(file.path(outdir, paste0(prefix, "_errorRates_FWD.pdf")))
plotErrors(errF, nominalQ=TRUE)
dev.off()

pdf(file.path(outdir, paste0(prefix, "_errorRates_REV.pdf")))
plotErrors(errR, nominalQ=TRUE)
dev.off()

# run dada2 denoising algorithm
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

# merge read pairs
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

## this step no longer necessary with improved demultiplexing...
# now, let’s reverse the reversed sequences…
#rev <- which(endsWith(names(mergers), "rev1.1.fastq"))
#for (sample in rev){
#	mergers[[sample]]$sequence <- sapply(mergers[[sample]]$sequence, rc)
#}#for

# make table of sequences
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# save sequence table- we’ll do chimera removal for all runs at the same time!
saveRDS(seqtab, file.path(outdir, paste0(prefix, "_seqtab.RDS")))

# save data
save.image(file.path(outdir, paste0(prefix, "_dada2.RData")))

# for some reason, R can't find dada2 at this point, so we'll check our libpaths...
if(length(.libPaths()==3)) .libPaths(.libPaths()[-1])

# output read counts from each step
getN <- function(x) sum(getUniques(x))
countNonZeros <- function(x) { sum(x!=0) }

remove <- which(!unique(rownames(out)) %in% names(dadaFs))

ifelse (length(remove)==0,
{track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), apply(seqtab, 1, countNonZeros))},
{track <- cbind(out[-remove,], sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), apply(seqtab, 1, countNonZeros))})

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "#ASVs")
rownames(track) <- names(dadaFs)
#head(track)

write.csv(track, file.path(outdir, paste0(prefix, "_dada2ReadCounts.csv")))

}#end function
