run_summary_plots <- function(run){
  # set up env
  library(ggplot2)
  library(viridis)
  library(LocaTT)
  setwd(file.path("/nfs4/FW_HMSC/Levi_Lab/Hallerud/ESRF/Arthropods/",run,"1_demult"))
  dir.create("plots")
  
  # read in & process read counts
  dirs <- list.dirs(".", recursive=FALSE)
  dirs <- dirs[grepl("ESFI",dirs)]
  all_reads <- data.frame()
  for (dir in dirs){
    name <- gsub("_","",basename(dir))
    # read in TSV
    tsv <- read.table(file.path(dir, paste0(name,"_reads.tsv")))
    names(tsv) <- c("Lines","File")
    tsv$Library <- name
    # assign sample names
    tsv$File <- basename(tsv$File)
    tsv$File <- gsub("_ACAGTG", "", tsv$File)
    tsv[c("Sample1","Sample2")] <- stringr::str_split(tsv$File, "-", 2, TRUE)
    tsv$Sample2 <- gsub("_R1.fastq","", tsv$Sample2)
    tsv <- tsv[!grepl("total",tsv$File),]
    # match sampleIDs to tags (to track if any particular tags are problematic)
    barcodes <- LocaTT::read.fasta(paste0("../barcodes/i5_",name,".fa"))
    barcodes$Sequence <- gsub("CCHGAYATRGCHTTYCCHCG","", barcodes$Sequence)
    barcodes$Name <- gsub(paste0(name,"_"),"", barcodes$Name)
    barcodes$Name <- sub("_[ACGT]{6}$", "", sub(paste0("^", name, "_"), "", barcodes$Name))
    tsv$Barcode1 <- barcodes$Sequence[match(tsv$Sample1, barcodes$Name)]
    tsv$Barcode2 <- barcodes$Sequence[match(tsv$Sample2, barcodes$Name)]
    print(table(is.na(tsv$Barcode1)))
    tsv$Barcode1[is.na(tsv$Barcode1)] <- "Missing"
    tsv$Barcode2[is.na(tsv$Barcode2)] <- "Missing"
    # merge
    all_reads <- rbind(all_reads, tsv)
  }#for

  lv <- sort(unique(c(as.character(all_reads$Barcode1), as.character(all_reads$Barcode2))))
  all_reads$Barcode1 <- factor(all_reads$Barcode1, levels = lv)
  all_reads$Barcode2 <- factor(all_reads$Barcode2, levels = lv)
  #all_reads$Lines <- as.numeric(trimws(all_reads$Lines))

  # check how many reads are missing a tag
  total <- sum(all_reads$Lines)
  missing_both <- sum(all_reads$Lines[all_reads$Sample1=="unknown" & all_reads$Sample2=="unknown"], na.rm=TRUE) / total
  missing_fwd <- sum(all_reads$Lines[all_reads$Sample1=="unknown" & all_reads$Sample2!="unknown"], na.rm=TRUE) / total
  missing_rev <- sum(all_reads$Lines[all_reads$Sample1!="unknown" & all_reads$Sample2=="unknown"], na.rm=TRUE) / total
  missing_tot <- missing_fwd+missing_rev+missing_both
  print(" ")
  print(paste("Reads missing 1+ tags: ", scales::percent(missing_tot)))
  print(paste("     Missing both tags: ", scales::percent(missing_both)))
  print(paste("     Missing FWD tag: ", scales::percent(missing_fwd)))
  print(paste("     Missing REV tag: ", scales::percent(missing_rev)))
  
  # check proportion of reads demultiplexed to tag jumps (12%)
  tagjumps <- sum(all_reads$Lines[all_reads$Sample1!=all_reads$Sample2 & all_reads$Sample1!="unknown" & all_reads$Sample2!="unknown"], na.rm=TRUE) / total
  print(paste("Demultiplexed but mismatched tags (i.e., jumps): ", scales::percent(tagjumps)))

  # check proportion of 'on-diagonal' reads- only 34%!
  print(paste("Proportion of successfully demultiplexed reads: ",
              scales::percent(sum(all_reads$Lines[all_reads$Sample1==all_reads$Sample2 & all_reads$Sample1!="unknown"]) / total, na.rm=TRUE)))

  # classify read fates
  all_reads$Type <- NA
  all_reads$Type[all_reads$Sample1==all_reads$Sample2] <- "proper"
  all_reads$Type[all_reads$Sample1!=all_reads$Sample2 & all_reads$Sample1!="unknown" & all_reads$Sample2!="unknown"] <- "tagjump"
  all_reads$Type[all_reads$Sample1!=all_reads$Sample2 & all_reads$Sample1=="unknown"] <- "missingFWD"
  all_reads$Type[all_reads$Sample1!=all_reads$Sample2 & all_reads$Sample2=="unknown"] <- "missingREV"
  all_reads$Type[all_reads$Sample2=="unknown" & all_reads$Sample1=="unknown"] <- "missingBOTH"
  all_reads$Type <- factor(all_reads$Type, levels=c("proper","tagjump","missingFWD","missingREV","missingBOTH"))
  
  # calculate on-target, missing, off-target rates per sample
  tags <- aggregate(Lines~Barcode1+Type+Library, all_reads, FUN=function(x) sum(x,na.rm=TRUE))
  tags <- tidyr::pivot_wider(tags, id_cols=c("Library","Barcode1"), names_from="Type", values_from="Lines", values_fill=0)
  tags$onesided <- tags$missingFWD+tags$missingREV
  tags$frac_diag <- tags$proper / (tags$proper + tags$onesided)
  head(tags)
  
  ## library level plots
  #all_reads$Lines[all_reads$Lines==0] <- NA
  for (p in unique(all_reads$Library)){
    print(paste("Making plots for",p))
    # plot...
    plot <- ggplot(data=all_reads[which(all_reads$Library==p),])+
      geom_tile(aes(x=Barcode1, y=Barcode2, fill=Lines+1))+
      scale_fill_viridis(option="turbo")+
      theme(axis.text=element_text(size=6),
            axis.text.x=element_text(angle=90), 
            axis.ticks=element_blank())+
      labs(fill="Reads", sub=p)+
      ggtitle("Proper Combinations vs. Tag-Jumping (includes missing)")
    ggsave(file.path("plots", paste(p,"TagJumpingPlot.png",sep="_")),
           plot, height=10, width=12)
  
    # on-target rates per tag
    plot <- ggplot(all_reads[which(all_reads$Barcode1!="Missing" & all_reads$Barcode2!="Missing" & all_reads$Library==p),])+
      geom_tile(aes(x=Barcode1, y=Barcode2, fill=Lines+1))+
      scale_fill_viridis(trans="log10")+
      theme(axis.text=element_text(size=6),
            axis.text.x=element_text(angle=90), 
            axis.ticks=element_blank())+
      labs(fill="Reads", sub=p)+
      ggtitle("Proper Combinations vs. Tag-Jumping")
    ggsave(file.path("plots", paste(p,"OnTargetsPerTag.png",sep="_")), 
           plot, width=12, height=12)
  
    # plot on-target vs. off-target
    plot <- ggplot(tags[which(tags$Library==p),])+ 
      geom_abline(slope=1, linetype=2)+
      geom_point(aes(onesided+1, proper+1), alpha=0.7, size=3)+
      geom_text(data = tags[which(tags$Library==p & tags$frac_diag<0.1),],
                aes(x=onesided+1, y=proper+1, label=Barcode1), size=2, col="red")+
      theme_classic()+
      scale_x_log10()+ 
      scale_y_log10()+
      xlab("Missing a barcode")+
      ylab("On-target")+
      ggtitle("Barcode Missingness Check")+
      labs(subtitle=p)
    ggsave(file.path("plots", paste(p,"TagsVsMissing.png",sep="_")), 
           plot, width=12, height=12)

    # plot reads per sample
    plot <- ggplot(all_reads[which(all_reads$Type!="missingBOTH" & all_reads$Library==p),])+ 
      geom_bar(aes(x=Sample1, y=Lines+1, fill=Type), stat="identity")+
      theme_classic()+
      xlab("Sample")+
      ylab("Demultiplexed Reads")+
      ggtitle("Read Fates")+
      labs(fill=element_blank(), sub=p)+
      facet_wrap(~Library)+
      coord_flip()
    ggsave(file.path("plots", paste(p,"ReadsPerRep.png",sep="_")),
                     plot, width=12, height=12)
  }#for
  
  ## run-level plots
  # tag jumping
  plot <- ggplot(all_reads)+
    geom_tile(aes(x=Barcode1, y=Barcode2, fill=Lines))+
    scale_fill_viridis(option="turbo", trans="log10")+
    theme(axis.text=element_text(size=6),
          axis.text.x=element_text(angle=90), 
          axis.ticks=element_blank())+
    labs(fill="Reads", sub=run)+
    ggtitle("Proper Combinations vs. Tag-Jumping (includes missing)")
  ggsave(file.path("plots", paste(run,"TagJumpingPlot.png",sep="_")),
         plot, height=10, width=12)
  
  # on-target rates per tag
  plot <- ggplot(all_reads[all_reads$Barcode1!="Missing" & all_reads$Barcode2!="Missing",])+
    geom_tile(aes(x=Barcode1, y=Barcode2, fill=Lines+1))+
    scale_fill_viridis(trans="log10")+
    theme(axis.text=element_text(size=6),
          axis.text.x=element_text(angle=90), 
          axis.ticks=element_blank())+
    labs(fill="Reads", sub=run)+
    ggtitle("Proper Combinations vs. Tag-Jumping")
  ggsave(file.path("plots", paste(run,"OnTargetsPerTag.png",sep="_")), 
         plot, width=12, height=12)
  
  # plot on-target vs. off-target
  plot <- ggplot(tags)+ 
    geom_abline(slope=1, linetype=2)+
    geom_point(aes(onesided+1, proper+1, col=Library), alpha=0.7, size=3)+
    theme_classic()+
    scale_x_log10()+ 
    scale_y_log10()+
    labs(subtitle=run)+
    xlab("Missing a barcode")+
    ylab("On-target")+
    ggtitle("Barcode Missingness Check")
  ggsave(file.path("plots", paste(run,"TagsVsMissing.png",sep="_")), 
         plot, width=12, height=12)
  
  # plot tags x reads
  plot <- ggplot(all_reads[which(all_reads$Barcode1!="Missing"),])+ 
    geom_bar(aes(x=Barcode1, y=Lines, fill=Library), stat="identity")+
    theme_classic()+
    xlab("Tag 1")+
    ylab("All Reads")+
    labs(subtitle=run)+
    ggtitle(paste(run, "Tag Fates"))+
    labs(fill=element_blank())+
    facet_wrap(~Type)+
    coord_flip()
  ggsave(file.path("plots", paste(run,"ReadsPerBarcode.png",sep="_")), 
         plot, width=12, height=12)
  
  # read allocation per lib
  plot <- ggplot(all_reads)+
    geom_bar(aes(x=Library, y=Lines, fill=Type), stat="identity")+
    xlab("Library")+
    ylab("Demultiplexed Reads")+
    labs(subtitle=run)+
    coord_flip()+
    theme_classic()+
    ggtitle(run)
  ggsave(file.path("plots",paste(run,"ReadsPerLibrary.png",sep="_")),
         plot, width=7,height=7)
}#function


args <- commandArgs(trailingOnly = TRUE)
run <- args[1]

run_summary_plots(run)

# zip folder for easier download
system("zip -r plots.zip plots")
