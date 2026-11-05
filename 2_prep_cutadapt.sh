#!/bin/bash
type=$1
lib=$2
R1=$3
R2=$4

while read well sample lib_index LIB FWDbarcode REVbarcode FWDprimer REVprimer locus primer region project; 
do
OUTDIR=1_demult/ESFI_${type}_P${lib}/demux
if [[ ! -d ${OUTDIR} ]]; then mkdir ${OUTDIR} -p; fi
echo "cutadapt -e 1 --no-indels --minimum-length 70 --overlap 28 --revcomp --discard-untrimmed -g ${FWDbarcode}${FWDprimer} -G ${REVbarcode}${REVprimer} -o ${OUTDIR}/LIB${type}${lib}_${sample}_${FWDbarcode}_R1.fastq -p ${OUTDIR}/LIB${type}${lib}_${sample}_${FWDbarcode}_R2.fastq 0_raw/lane1*-${LIB}_*R1_001.fastq.gz 0_raw/lane1*-${LIB}_*R2_001.fastq.gz"
done  < barcodes/ESFI${type}P${lib}_barcodes.tsv


