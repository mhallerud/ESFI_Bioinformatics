#!/bin/bash
type=$1
lib=$2

OUTDIR=2_paired/ESFI_${type}_P${lib}
if [[ ! -d ${OUTDIR} ]]; then mkdir ${OUTDIR} -p; fi

# -v : minimum overlap for pairing
# -q : minimum quality score- if two consecutive bases are below this score, the read is trimmed
# -u : maximum proportion of uncalled bases per read

for file in `ls -1 1_demult/ESFI_${type}_P${lib}/demux/*R1.fastq | cut -d "." -f1 | sed 's/_R1//g'`; do
	OUT="$(basename $file)"
	echo "pear -f ${file}_R1.fastq -r ${file}_R2.fastq -o ${OUTDIR}/${OUT} -q 20 -v 20 -u 0.05"
done
