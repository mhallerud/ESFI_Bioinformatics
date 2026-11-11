#!/bin/bash
prefix=$1
threads=$2

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

blastn -db ${SCRIPT_DIR}/../BOLD/BOLD_Public.04-Apr-2025 -query ${prefix}_OTUs.fasta -max_target_seqs 10 -evalue 1e-5 -out ${prefix}_OTUs_BOLD_BLAST.tsv \
-outfmt "6 qseqid sseqid sscinames scomnames staxids pident qcovs evalue bitscore gaps gapopen qseq sseq"  -num_threads ${threads}


