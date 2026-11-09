######################################################################################
#### Removes chimeras from a library using dada2's algorithm                      ####
#### usage: Rscript dada2_removeChimeras.R <seqtab> <fasta> <outprefix> <threads> ####
######################################################################################
library(dada2)
library(LocaTT)
library(stringr)


# read in arguments
args <- commandArgs(trailingOnly = TRUE)
SEQTAB=args[1]
SEQS=args[2]
OUT=args[3]
THREADS=args[4]

# read in sequences
seqs <- LocaTT::read.fasta(SEQS)
seqs$Name <- stringr::str_split(seqs$Name, ";", 2, TRUE)[,1]

# read in sequence table from vsearch & reformat
seqtab <- read.csv(SEQTAB, sep="\t")
rownames(seqtab) <- seqtab[,"X.OTU.ID"]
rownames(seqtab) <- seqs$Sequence[match(rownames(seqtab), seqs$Name)]
seqtab <- seqtab[,-1]
seqtab <- t(seqtab)

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = THREADS, verbose = FALSE)
write.csv(seqtab.nochim, paste0(OUT,"_ASVs.csv"))

# export sequences
outseqs <- colnames(seqtab.nochim)
outnames <- seqs$Name[match(outseqs, seqs$Sequence)]
outsizes <- paste0("size=", colSums(seqtab))
outnamesfull <- paste(paste(OUT, outnames, sep="_"), outsizes, sep=";")
LocaTT::write.fasta(outnamesfull, outseqs, paste0(OUT, "_ASVs.fasta"))
