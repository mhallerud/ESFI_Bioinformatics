#!/bin/bash
DIR=$1

# set up names, dirs
OUT="$(echo $DIR | sed 's/_//g')"
OUTDIR=1_demult/${DIR}/demux
if [[ ! -d ${OUTDIR} ]]; then mkdir ${OUTDIR} -p; fi

# set up barcodes files
BARCODES=barcodes/${OUT}_barcodes.tsv
awk -v lib="${OUT}" '{print ">"lib"_"$2"_"$3"\n"$5$7}' ${BARCODES} > barcodes/i5_${OUT}.fa
awk -v lib="${OUT}" '{print ">"lib"_"$2"_"$3"\n"$6$8}' ${BARCODES} > barcodes/i7_${OUT}.fa
LIB="$(awk '{print $4}' $BARCODES | head -n1)"

# run cutadapt
# excluding adapter trimming b/c seqs with both barcodes & adapters get shunted to 1-2.fastq instead of {name1}-{name2}.fastq
cutadapt -e 0 --no-indels --minimum-length 100 --overlap 28 --revcomp --discard-untrimmed -g file:barcodes/i5_${OUT}.fa -G file:barcodes/i7_${OUT}.fa -o ${OUTDIR}/{name1}-{name2}_R1.fastq -p ${OUTDIR}/{name1}-{name2}_R2.fastq 0_raw/lane1*-${LIB}_*R1_001.fastq.gz 0_raw/lane1*-${LIB}_*R2_001.fastq.gz -j 10

# clean up tag-jumping files
for f in `ls -1 ${OUTDIR}`; do
	b=$(basename "$f" | sed 's/.fastq//')
	n1=${b%%-*}; n2=${b#*-}
	n2b="$(echo $n2 | sed 's/_R[12]//')"
	if [[ $n1 != $n2b ]]; then rm ${OUTDIR}/$f; 
	else mv ${OUTDIR}/${f} ${OUTDIR}/${n2}.fastq; fi
done
