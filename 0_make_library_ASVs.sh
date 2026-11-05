#!/bin/bash
######################################################################################################
### *must be run from within run directory, which must contain:                                    ###
###         0_raw : directory containing raw fastq.gz files                                        ###
###         barcodes : barcodes directory containing TSV files in 11-field format                  ###
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

LIB_INDEX="$(awk '{print $3}' barcodes/ESFI${type}P${lib}_barcodes.tsv | head -n1)"

OUTDIR=1_demult/ESFI_${type}_P${lib}
if [[ ! -d ${OUTDIR} ]]; then mkdir ${OUTDIR} -p; fi
cd ${OUTDIR}

# trim primers (min length 200)
cutadapt -e 0.15 --overlap 20 --no-indels --revcomp -g CCHGAYATRGCHTTYCCHCG -G TCDGGRTGNCCRAARAAYCA -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -o primers_trimmed_R1.fastq.gz -p primers_trimmed_R2.fastq.gz \
../../0_raw/*${LIB_INDEX}*R1_001.fastq.gz ../../0_raw/*${LIB_INDEX}*R2_001.fastq.gz -j ${threads}

# dereplicate- check these long reads that are passing
seqtk seq primers_trimmed_R1.fastq.gz -A > primers_trimmed_R1.fasta
#vsearch --fastx_filter primers_trimmed_R1.fasta --minseqlength 200 --output primers_trimmed_R1_long.fasta
vsearch --derep_fulllength primers_trimmed_R1_long.fasta --minseqlength 200 --output trimmed_long_ESVs.fasta --sizeout --relabel esv --minuniquesize 100 --threads ${threads}

# now do the same with first 1 million reads to check short reads
#vsearch --fastx_filter primers_trimmed_R1.fasta --maxseqlength 200 --output primers_trimmed_R1_short.fasta
vsearch --derep_fulllength primers_trimmed_R1_short.fasta --maxseqlength 200 --output trimmed_short_ESVs.fasta --sizeout --relabel esv --minuniquesize 100 --threads ${threads}


# clean up
rm primers_trimmed*fasta -f
#rm primers_trimmed*fastq.gz -f

