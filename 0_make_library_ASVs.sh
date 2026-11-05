#!/bin/bash
######################################################################################################
### *must be run from within 1_demult directory                                                    ###
### USAGE: 0_make_library_ASVs.sh <type> <lib> <threads>                                           ###
###         type : PF or MT                                                                        ###
###         lib : 01-44                                                                            ###
###         threads : number of processors to use in multi-threading                               ###
### type & lib are used to locate the library directory& barcodes  (specific to ESFI project)      ###
### OUTPUTS: ESFI_PF_P01/trimmed_R1_ESVs.fasta = exact sequence variants for R1 files              ###
######################################################################################################

type=$1
lib=$2
threads=$3

LIB_INDEX="$(awk '{print $3}' ../barcodes/ESFI${type}P${lib}_barcodes.tsv | head -n1)"

cd ESFI_${type}_P${lib}

# trim primers (min length 200)
cutadapt -e 0.15 --no-indels --minimum-length 100 --revcomp -g CCHGAYATRGCHTTYCCHCG -G TCDGGRTGNCCRAARAAYCA -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -o primers_trimmed_R1.fastq.gz -p primers_trimmed_R2.fastq.gz \
../../0_raw/*${LIB_INDEX}*R1_001.fastq.gz ../../0_raw/*${LIB_INDEX}*R2_001.fastq.gz -j ${threads}

# dereplicate
seqtk seq primers_trimmed_R1.fastq.gz -A > primers_trimmed_R1.fasta
vsearch --derep_fulllength primers_trimmed_R1.fasta --output trimmed_R1_ESVs.fasta --sizeout --relabel esv --minuniquesize 100 --threads ${threads}

# clean up
rm primers_trimmed*fasta -f
rm primers_trimmed*fastq.gz -f

