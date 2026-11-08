#!/bin/R

library(data.table)
library(stringr)
library(tidyr)


main <- function(){
	# read in command line arguments
	print("INPUTS:")
	args <- commandArgs(trailingOnly=TRUE)
	seqtab = args[1]
	min_reps = args[2]
	min_reads = args[3]
	outfile = args[4]
	prop_reads = FALSE
	print(paste("Seqtab:", seqtab))
	print(paste("Min Reads:", min_reads))
	print(paste("Min Reps:", min_reps))
	print(paste("Outfile:", outfile))
	print("")
	
	# read in OTUtab
	otutab <- data.table::fread(seqtab)
	colnames(otutab)[2:ncol(otutab)] <- paste0("otu", 1:(ncol(otutab)-1))
	
	# if prop_reads, convert to proportion of reads per sample
	if(prop_reads){
		apply(otutab, 2, FUN=function(x) x / sum(x))
	}#prop_reads
	
	# convert otutab to long format
	print("Converting table to long format....")
	otulong <- tidyr::pivot_longer(otutab, cols=names(otutab)[2:ncol(otutab)])
	names(otulong) <- c("PCR.ID","OTU","reads")
	otulong[,c("Plate","SampleID")] <- stringr::str_split(otulong$PCR.ID, "\\.", n=2, simplify=TRUE)
	otulong[,c("SampleID","Rep")] <- stringr::str_split(otulong$SampleID, "_", n=2, simplify=TRUE)
	
	# filter by reads first
	print("Filtering by read counts....")
	otulong$reads[otulong$reads < min_reads] = 0
	
	# filter by reps next
	print("Filtering by replicates....")
	# first make table of reps per OTU and sample combination
	otulong2 <- otulong[-which(otulong$reads==0),]
	samplexreps <- as.data.frame(table(otulong2$OTU, otulong2$SampleID))
	samplexreps <- samplexreps[which(samplexreps$Freq!=0),]
	# now make list of OTU and sample combinations where minimum rep count isn't met
	exclude <- samplexreps[samplexreps$Freq<min_reps,]
	names(exclude) <- c("OTU","Sample","Reps")
	# set these to 0 in main DF
	for (i in 1:nrow(exclude)){
		rm <- exclude[i,]
		otulong$reads[which(otulong$OTU==rm$OTU & otulong$SampleID==rm$Sample)] <- 0
	}#for i
	#otulong <- otulong[-which(otulong$reads==0),]
	
	# convert back to OTU table with OTUs in rows, samples in columns
	print("Converting back to wide format....")
	out_long <- otulong[,c("OTU","PCR.ID","reads")]
	out <- tidyr::pivot_wider(out_long, names_from="PCR.ID", values_from="reads")
	replace <- as.list(rep(0, ncol(out)))
	names(replace) <- names(out)
	out <- tidyr::replace_na(out, replace)
	# save
	print("Exporting filtered table....")
	saveRDS(out, paste0(outfile,".RData"))
	write.csv(as.data.frame(out), paste0(outfile,".csv"), row.names=FALSE)
	#return(out)
}#main



main()
