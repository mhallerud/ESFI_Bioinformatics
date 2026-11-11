## usage : Rscript merge_otutabs.R <RUN>
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
RUN=args[1]
setwd(file.path("/nfs4/FW_HMSC/Levi_Lab/Hallerud/ESRF/Arthropods", RUN, "3_clustering"))

# grab otutab filepaths
allfiles <- list.files(".")
files <- allfiles[grepl("_OTUs_mapped.txt", allfiles)]
files <- files[!startsWith(files, RUN)]

# function to read & reformat OTUtab
readOTUtab <- function(file){
	otutab <- read.csv(file, sep="\t")
	rownames(otutab) <- otutab$X.OTU.ID
	otutab <- otutab[,-which(names(otutab)%in%c("X.OTU.ID"))]
	print(dim(otutab))
	return(otutab)
}#readOTUtab

full_otutab <- readOTUtab(files[1])
if(length(files)>1){
	# merge in each additional seqtab
	for (otutab in files[2:length(files)]){
		print(paste("Merging in", otutab))
		ot <- readOTUtab(otutab)
		full_otutab <- cbind(full_otutab, ot)
	}#for
}#if

# save
outpath = paste0(RUN, "_OTUs_mapped.csv")
write.csv(full_otutab, outpath)
print(paste("Successfully merged at :", outpath))
