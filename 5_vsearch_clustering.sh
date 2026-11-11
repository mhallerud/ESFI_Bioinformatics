#!/bin/bash
RUN=$1
THREADS=$2

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

cd /nfs4/FW_HMSC/Levi_Lab/Hallerud/ESRF/Arthropods/${RUN}/3_clustering


# setup for previous version that uses dada2 to create quality filtering...
# create re-replicated reads based on dada2 denoised output (per PCR replicate)
# this is necessary so that read counts can be dereplicated at the library level
#for dir in `ls -d 2_denoise/ESFI*`; do
#	dir="$(basename $dir)"
#	[[ ! -d 3_clustering/${dir} ]]; then mkdir 3_clustering/${dir} -p; fi
#	for f in `ls -1 2_denoise/${dir}/*clust.fasta`; do
#		vsearch --rereplicate 2_denoise/${dir}/${f} > 3_clustering/${dir}/${f}
#	;done
#;done


# merge sequence tables for this run & filter chimeras 
# dada2 is better than vsearch because it accounts for 'consensus' chimera identity,
# which benefits from pulling information across the entire run
echo "MERGING SEQUENCE TABLES & FILTERING CHIMERAS........."
/local/cqls/R/R-4.4.2/Rscript ${SCRIPT_DIR}/merge_seqtabs.R ${RUN}
echo "\n\n"


# Cluster ASVs into 97% OTUs
echo "CLUSTERING INTO 97% OTUs............"
vsearch --cluster_size ${RUN}_ASVs.fasta --id 0.97 --centroids ${RUN}_OTUs.fasta --sizein --relabel ${RUN}otu --sizeout --clusterout_sort \
--consout ${RUN}_OTUs_consensus.fasta --msaout ${RUN}_OTUs.msa --uc ${RUN}_OTUs.uc --alnout ${RUN}_OTUs.aln

# clean up uc file
grep -v "^C" ${RUN}_OTUs.uc > asvgroups.temp
cut -f2 asvgroups.temp | while read l;  do echo "otu$(($l + 1))"; done > otus.temp
cut -f3-9 asvgroups.temp > rest.temp
paste otus.temp rest.temp > ${RUN}_OTUs.clean.uc
rm *temp -f

# alternative: swarm OTUs
#/nfs4/FW_HMSC/Levi_Lab/Hallerud/ESRF/Arthropods/swarm/bin/swarm -z ${RUN}_ASVs.fasta -d 9 -w ${RUN}_swarmOTUs.fasta -o ${RUN}_swarmOTUs.txt

# map OTUs within samples
echo "MAPPING TO LIBRARIES..............."
# run PER LIBRARY (otherwise will take forever)
for fasta in `ls -1 ESFI*filt.fasta`; do
	echo "            " $fasta
	lib=${fasta%_filt.fasta}
	vsearch --usearch_global ${fasta} --db ${RUN}_OTUs.fasta --id 0.97 --otutabout ${lib}_OTUs_mapped.txt
done

# merge across runs
/local/cqls/R/R-4.4.2/Rscript ${SCRIPT_DIR}/merge_otutabs.R ${RUN}
