library(dada2)
library(here)
script_path=here::here()


# set input variables
## USAGE: Rscript 6_run_dada2_assign.R infile prefix [type]
args <- commandArgs(trailingOnly = TRUE)
INFILE=args[1]
OUT=args[2]

ifelse(length(args)==3, type<-args[3], type<-"FASTA")

currdir=getwd()
#setwd(file.path(currdir, path.cut))

# set up sequence data
print("Reading in data...........")
if(type=="SEQTAB"){
	SEQTAB=INFILE
	seqs <- readRDS(SEQTAB)
}else if(type=="FASTA"){
	library(LocaTT)
	seqs <- LocaTT::read.fasta(INFILE)
	seqs <- seqs[,2]
}#IFELSE


# taxonomic assignment
print("Starting taxonomic assignment....")
set.seed(100)
taxa <- assignTaxonomy(seqs, file.path(script_path, "../BOLD/BOLD_Public.04-Apr-2025.dada2.fasta"), 
                       minBoot=50, outputBootstraps=TRUE, multithread=TRUE)
write.csv(taxa, paste0(PREFIX, "_taxa.csv"))
