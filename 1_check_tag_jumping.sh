#!/bin/bash
######################################################################################################################################
### *must be run from within run directory, which must contain:                                                                    ###
###         0_raw : directory containing raw fastq.gz files                                                                        ###
###         barcodes : barcodes directory containing TSV files in 11-field format                                                  ###
### USAGE: 1_check_tag_jumping.sh <type> <lib> <threads>                                                                           ###
###          type : MF or PF                                                                                                       ###
###          lib : 01-44                                                                                                           ###
###          threads : number of threads                                                                                           ###
### type & lib are used to locate the library directory& barcodes  (specific to ESFI project)                                      ###
### OUTPUTS: ESFI_PF_P01/trimmed_R1_ESVs.fasta = exact sequence variants for R1 files                                              ###
###          demultiplexing_read_category_ESFI${type}P${lib}.tsv = assignment of each read to a demultiplexing category            ###
###          ESFI${type}P${lib}_reads.tsv =  lines per demultiplexed file (based on barcodes only)                                 ###
###          ESFI${type}P${lib}_R1_cutadapt.fastq.gz = concatenated reads into single file,                                        ###
###                                                         with demultiplexed fate added to read header                           ###
###          ESFI${type}P${lib}_R1_ASVs_mapped.txt = mapping of raw R1 reads against most abundant 100 ASVs                        ###
######################################################################################################################################
type=$1
lib=$2
threads=$3

OUTDIR=1_demult/ESFI_${type}_P${lib}
if [[ ! -d ${OUTDIR} ]]; then mkdir ${OUTDIR}; fi
cd ${OUTDIR}

LIB_INDEX="$(awk '{print $4}' ../../barcodes/ESFI${type}P${lib}_barcodes.tsv | head -n1)" 

# set up barcode only files, account for anchoring of all NNN${barcode} options 
awk '{print ">"$2"\n^NNN"$5}' ../../barcodes/ESFI${type}P${lib}_barcodes.tsv > i5_barcodes.fa

# run cutadapt for anchored unique barcodes with 0 errors, min length 110 (primers=40 bp, insert=70 bp)
if [[ ! -d cutadapt ]]; then mkdir cutadapt; fi
cutadapt -e 0 --no-indels --overlap 8 --minimum-length 70 -g file:i5_barcodes.fa -G file:i5_barcodes.fa -o cutadapt/{name1}-{name2}_R1.fastq -p cutadapt/{name1}-{name2}_R2.fastq \
../../0_raw/*-${LIB_INDEX}_*_R1_001.fastq.gz ../../0_raw/*-${LIB_INDEX}_*_R2_001.fastq.gz -j ${threads}

# grab word counts to check tag jumping (check in R)
wc -l cutadapt/*R1.fastq > ESFI${type}P${lib}_reads.tsv


# label reads with demultiplex outcome based on these outputs
cd cutadapt
for f in *_R1.fastq; do 
	b=${f%_R1.fastq}
	a=${b%%-*}
	c=${b#*-} 
	if [[ $a == unknown && $c == unknown ]]; then k=both_unknown; 
	elif [[ $a == unknown && $c != unknown ]]; then k=one_sided_FWD; 
	elif [[ $a != unknown && $c == unknown ]]; then k=one_sided_REV; 
	elif [[ $a == $c ]]; then k=diag; 
	else k=offdiag; 
	fi 
	awk -v k=$k 'NR%4==1{sub(/^@/,"",$1)
	print $1"\t"k}' $f 
done > ../demultiplexing_read_category_ESFI${type}P${lib}.tsv
cd ..

# concatenate R1s into single file with samples added to header filename
for f in `ls -1 cutadapt/*R1.fastq`; do
	sed -e "s|\(^@VH.*\) .*$|\1;sample=${f%.*};|g" $f | sed 's/.assembled//g'| sed 's|2_paired/||g' >> ESFI${type}P${lib}_R1_cutadapt.fastq
done

# map against top 100 most abundant R1 ASVs in library
head -n 200 trimmed_long_ESVs.fasta > trimmed_R1_ESVs_top100.fasta
vsearch --usearch_global ESFI${type}P${lib}_R1_cutadapt.fastq --db trimmed_R1_ESVs_top100.fasta --id 0.97 --otutabout ESFI${type}P${lib}_R1_ASVs_mapped.txt
#makeblastdb -in trimmed_R1_ESVs_top100.fasta -dbtype nucl
#seqtk seq -A ../0_raw/*${LIB_INDEX}*_R1_001.fastq.gz > raw_R1.fa
#hqsub -c "blastn -query raw_R1.fa -db mock_refs -outfmt "6 qseqid sseqid pident length qcovs" -max_target_seqs 1 -perc_identity 95 -num_threads 10 > raw_map_morphos.tsv" -p 10 -q levi_lab -r map_mcraw_morphos_log


# clean up
rm i5_barcodes.fa
rm raw_R1.fa
rm cutadapt/*fastq -f
rmdir cutadapt
gzip ESFI${type}P${lib}_R1_cutadapt.fastq
