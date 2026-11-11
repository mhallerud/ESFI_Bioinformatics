#!/bin/bash
## Must be run from within Run directory

LIB=$1
THREADS=$2

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )


# set filepaths / names
OUT="$(echo ${LIB} | sed 's/_//g')"
if [[ ! -d 3_clustering ]]; then mkdir 3_clustering; fi

# clean up extra pear outputs
peardir=2_paired/${LIB}/pear
if [[ ! -d ${peardir} ]]; then mkdir ${peardir} -p; fi
mv 2_paired/${LIB}/*unassembled*fastq ${peardir}
mv 2_paired/${LIB}/*discarded*fastq ${peardir}

# concatenate reads into single FASTQ per library- sampleID is appended to each read header
echo "MAKING FASTQ FILES FOR " ${LIB} "............"
cd 2_paired/${LIB}
rm ../../3_clustering/${OUT}_full.fastq -f
for f in `ls -1 *assembled.fastq`; do
        sed -e "s|^\(@VH[^ ]*\) .*$|\1;sample=${f%.*};|" $f | sed 's/.assembled//g'| sed 's|2_paired/||g' >> ../../3_clustering/${OUT}_full.fastq
done

cd ../../3_clustering
# Quality filtering
# maxee_rate=0.01 is expected error of 1 per 100 bp, e.g. where EE within 100 bp ~ Pois(lambda=1)
# see vsearch docs: https://torognes.github.io/vsearch/misc/vsearch-expected_error.7.html
echo "QUALITY FILTERING FOR " ${LIB} "............"
vsearch --fastx_filter ${OUT}_full.fastq --fastq_maxee_rate 0.01 --fastaout ${OUT}_filt.fasta --sizeout --threads ${THREADS}

# Dereplicate & add abundance to headers
vsearch --derep_fulllength ${OUT}_filt.fasta --output ${OUT}_uniques.fasta --sizeout --minuniquesize 5 --threads 10

# Denoising & frequency filtering: create exact sequence variants
# paper describing algorithm: https://www.biorxiv.org/content/10.1101/081257v1.full
# minimum reads per ESV: 10
# defaults: minsize 8, unoise_alpha 2
# sizein: read size info of input; centroid_sizeout: output size of centroids
# clusterout_sort : sort outputs by abundance
# centroids: writes sequence centroids (i.e., the first seq in the )
# msaout : multiple sequence alignments & consensus sequences
# consout : output of consensus sequences
# alnout : output of (human-readable) alignments
echo "DENOISING " ${LIB} "............"
rm ${OUT}_ESVs* -f
vsearch --cluster_unoise ${OUT}_uniques.fasta --minsize 10 --unoise_alpha 2 --sizein --sizeout --clusterout_sort --relabel esv \
--centroids ${OUT}_ESVs.fasta --msaout ${OUT}_ESVs.msa --consout ${OUT}_ESVs_consensus.fa --alnout ${OUT}_ESVs.aln --threads ${THREADS}

# Length filtering
# keep 6-codon difference from expected amplicon
vsearch --fastx_filter ${OUT}_ESVs.fasta --fastq_minlen 400 --fastq_maxlen 436 --fastaout ${OUT}_ESVs_400to436bp.fasta

# map these back to samples to create current seqtab
vsearch --search_exact ${OUT}_filt.fasta --db ${OUT}_ESVs_400to436bp.fasta --otutabout ${OUT}_ESVs_mapped.txt --threads ${THREADS}

# Chimera filtering with UCHIME2
# key parameters:
# dn : pseudo-count prior on 'no votes' in chimera scoring fn-- increase to increase likelihood of identifying chimeras (default: 1.4)
# xn : weight of 'no votes' in chimera scoring fn-- increase to reduce likelihood of identifying chimeras (default: 8)
# abskew : abundance skew between (default: 16)
#vsearch --uchime3_denovo ${OUT}_ESVs.fasta --sizein --fasta_score --sizeout --relabel ASV_${OUT} --nonchimeras ${OUT}_ASVs.fasta \
#--chimeras ${OUT}_chimeras.fasta --uchimealns ${OUT}_chimeras.aln --uchimeout ${OUT}_uchime.tsv

# We'll do chimera filtering at the run level with dada2 instead
