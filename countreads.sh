#!/bin/bash

p=$1

echo "COUNTING READS FOR PLATE $p"

echo "Demultiplexed"
for f in 1_demult/ESFI_*_P${p}/demux/*R1.fastq; do echo "$f $(( $(wc -l < $f) / 4 ))"; done | awk '{sum+=$2}END{print sum}'

echo "Paired"
for f in 2_paired/ESFI_*_P${p}/*assembled.fastq; do echo "$f $(( $(wc -l < $f) / 4 ))"; done | awk '{sum+=$2}END{print sum}'

echo "Denoised"
grep size 3_clustering/ESFI*P${p}_ESVs.fasta | sed 's/size=/\t/g' | awk '{sum+=$2}END{print sum}'

echo "Length filtered"
vsearch --fastx_filter 3_clustering/ESFI*P${p}_ESVs.fasta --fastq_maxlen 436 --fastq_minlen 400 --fastaout 3_clustering/ESFI*P${p}_ESVs_400to436bp.fasta
grep size 3_clustering/ESFI*P${p}_ESVs_400to436bp.fasta | sed 's/size=/\t/g' | awk '{sum+=$2}END{print sum}'

