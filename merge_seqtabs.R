library(dada2)
library(LocaTT)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
RUN=args[1]
setwd(file.path("/nfs4/FW_HMSC/Levi_Lab/Hallerud/ESRF/Arthropods", RUN, "3_clustering"))

# grab seqtab filepaths
allfiles <- list.files(".")
files <- allfiles[grepl("_ESVs_mapped.txt", allfiles)]
#files <- sapply(files, function(x) file.path("3_clustering", x))

# function to read in all seqtabs & reformat
reformatSeqtab <- function(seqtab, seqfa){
	seqtab <- read.csv(seqtab, sep="\t")
	rownames(seqtab) <- seqtab[,1]
	seqtab <- seqtab[,-1]	

	seq_fa <- LocaTT::read.fasta(seqfa)
	seq_fa$Name <- stringr::str_split(seq_fa$Name, ";", 2, TRUE)[,1]
	rownames(seqtab) <- seq_fa$Sequence[match(rownames(seqtab), seq_fa$Name)]

	return(t(seqtab))
}#reformatSeqtab

# function to export sequences
exportSeqs <- function(seqtab, prefix, outfile){
	outseqs <- colnames(seqtab)
	outnames <- paste0(prefix, 1:length(outseqs))
	outsizes <- paste0("size=", colSums(seqtab))
	outnamesfull <- paste(outnames, outsizes, sep=";")
	LocaTT::write.fasta(outnamesfull, outseqs, outfile)
}#exportSeqs

# set up initial seqtab
fa1 <- gsub("_mapped.txt", ".fasta", files[1])
full_seqtab <- reformatSeqtab(files[1], fa1)

# merge in each additional seqtab
if (length(files)>1){
	for (seqtab in files[2:length(files)]){
		# grab corresponding fa file
		fa <- gsub("_mapped.txt", ".fasta", seqtab)
		# read in reformatted seqtab
		st <- reformatSeqtab(seqtab, fa)
		# merge
		full_seqtab <- dada2::mergeSequenceTables(full_seqtab, st)
	}#for
}#if

# check headers.
#colnames(pf01)
#sum(rownames(pf01)!=rownames(pf02))

# save
write.csv(full_seqtab, file=paste(RUN,"ESV_seqtab.csv", sep="_"))
exportSeqs(full_seqtab, "ESV", paste0(RUN, "_ESVs.fasta"))

# run chimera filtering on full seqtab
seqtab.nochim <- removeBimeraDenovo(full_seqtab, method = "consensus", multithread = TRUE, verbose = FALSE)

# export sequences
write.csv(seqtab.nochim, paste0(RUN,"_ASVs_mapped.csv"))
exportSeqs(seqtab.nochim, "ASV", paste0(RUN, "_ASVs.fasta"))


# reporting:
print(paste("Samples in ESV table:", dim(full_seqtab)[1]))
print(paste("Unique sequences in ESV table:", dim(full_seqtab)[2]))
print(paste("Samples in ASV table:", dim(seqtab.nochim)[1]))
print(paste("Unique sequences in ASV table:", dim(seqtab.nochim)[2]))
