library(dada2)
library(LocaTT)

args <- commandArgs(trailingOnly = TRUE)
RUN=args[1]
setwd(file.path("/nfs4/FW_HMSC/Levi_Lab/Hallerud/ESRF/Arthropods", RUN))

# grab seqtab filepaths
allfiles <- list.files("3_clustering/")
files <- allfiles[grepl("_ESVs_mapped.txt", allfiles)]
files <- sapply(files, function(x) file.path("3_clustering", x))

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

# set up initial seqtab
fa1 <- gsub("_mapped.txt", ".fasta", files[1])
full_seqtab <- reformatSeqtab(files[1], fa1)

# merge in each additional seqtab
for (seqtab in files[2:length(files)]){
	# grab corresponding fa file
	fa <- gsub("_mapped.txt", ".fasta", seqtab)
	# read in reformatted seqtab
	st <- reformatSeqtab(seqtab, fa)
	# merge
	full_seqtab <- dada2::mergeSequenceTables(full_seqtab, st)
}#for

# check headers.
#colnames(pf01)
#sum(rownames(pf01)!=rownames(pf02))

# save
write.csv(full_seqtab, file=paste(RUN,"ESV_seqtab.csv", sep="_"))
