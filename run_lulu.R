### USAGE: Rscript run_lulu.R otutab alignfile outfile
# assumes otutab is in dada2 output format


library(data.table)
library(lulu)
setwd(".")

main <- function(){
	# read in command line arguments
	args <- commandArgs(trailingOnly=TRUE)
	otufile = args[1]
	alignfile = args[2]
	outfile = args[3]
	print(paste("OTU table:", otufile))
	print(paste("Self-aligns:", alignfile))
	print(paste("Output:", outfile))
	
	# read in OTU table, reformat
	seqtab <- read.csv(otufile)
	#otutab <- t(seqtab[,-1])
	rownames(otutab) <- otutab[,1] #paste0("otu",1:nrow(otutab))
	otutab <- otutab[,-1]
	#colnames(otutab) <- seqtab$V1

	#otutab <- otutab[,-1]
	#otutab <- otutab[,-ncol(otutab)]
	#dim(otutab)
	
	# read in self-alignment file
	aligns <- fread("Runs1to3_OTUs_selfblast.txt", sep="\t")
	
	# run LULU	
	curated <- lulu::lulu(as.data.frame(otutab), as.data.frame(aligns), 
                             minimum_ratio_type='min', #error needs to have lower abundance in all samples 'min' or just on avg 'avg'
			     minimum_ratio=1, #minimum abundance ratio between potential error and parent
                             minimum_match=0.84,#minimum threshold of sequence similarity between parent-error
                             minimum_relative_cooccurence=0.95) #minimum cooccurrence rate of error and parent
	
	# save output
	write.csv(curated$curated_table, outfile)
}#main

main()

