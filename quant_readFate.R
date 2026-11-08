library(stringr)

args <- commandArgs(trailingOnly = TRUE)
TSV <- args[1]

# read in & reformat read count info
reads <- read.csv(TSV, header=FALSE)
reads$V1 <- trimws(reads$V1)
reads <- reads[!grepl("total", reads$V1),]
reads[,c("Reads","File")] <- stringr::str_split(reads$V1, " ", 2, TRUE)
reads[,c("Sample1","Sample2")] <- stringr::str_split(gsub("_R1.fastq","",reads$File),"-",2,TRUE)
reads$Reads <- as.numeric(reads$Reads)

# classify read type
reads$Type <- NA
reads$Type[reads$Sample1==reads$Sample2] <- "diag"
reads$Type[which(reads$Sample1!=reads$Sample2 & reads$Sample1!="unknown" & reads$Sample2!="unknown")] <- "tagjump"
reads$Type[reads$Sample1=="unknown" & reads$Sample2!="unknown"] <- "FWD_missing"
reads$Type[reads$Sample1!="unknown" & reads$Sample2=="unknown"] <- "REV_missing"
reads$Type[reads$Sample1=="unknown" & reads$Sample2=="unknown"] <- "BOTH_missing"

# tabulate
tbl <- aggregate(Reads~Type, reads, FUN=sum)
tbl$Prop <- tbl$Reads/sum(reads$Reads)
tbl
