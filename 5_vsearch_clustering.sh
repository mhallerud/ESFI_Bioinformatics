#!/bin/bash
RUN=$1
cd /nfs4/FW_HMSC/Levi_Lab/Hallerud/ESRF/Arthropods/${RUN}

# create re-replicated reads based on dada2 denoised output (per PCR replicate)
# this is necessary so that read counts can be dereplicated at the library level
for dir in `ls -d 2_denoise/ESFI*`; do
	dir="$(basename $dir)"
	[[ ! -d 3_clustering/${dir} ]]; then mkdir 3_clustering/${dir} -p; fi
	for f in `ls -1 2_denoise/${dir}/*clust.fasta`; do
		vsearch --rereplicate 2_denoise/${dir}/${f} > 3_clustering/${dir}/${f}
	;done
;done

# concatenate reads for each library - sampleIDs appended to read headers
for dir in `ls -d 3_clustering/ESFI*`; do
	for f in `ls -1 3_clustering/${dir}/*.fasta`; do
		dir="$(basename $dir)"
		sed -e "s|\(^@VH.*\) .*$|\1;sample=${f%.*};|g" $f | sed 's/.assembled//g'| sed 's|2_paired/||g' >> 3_clustering/${dir}_full.fasta
	;done
;done

# create amplicon sequence variants per library
if [[ ! -d 3_clustering/ASVs ]]; then mkdir -p 3_clustering/ASVs; fi
for dir in `ls -d 2_denoise`; do
	dir="$(basename $dir)"
	vsearch --derep_fulllength 3_clustering/${dir}_full.fasta --sizeout --relabel ${dir} --output 3_clustering/ASVs/${dir}_ASVs.fasta --uc 3_clustering/ESVs/${dir}_ASVs.uc --threads 10
done

# cluster OTUs
cd 3_clustering
cat ASVs/*ASVs.fasta > ${RUN}_ASVs.fasta
vsearch --cluster_size ${RUN}_ASVs.fasta --id 0.97 --centroids ${RUN}_OTUs.fasta --sizein --relabel otu --uc ${RUN}_OTUs.uc

# clean up uc file
grep -v "^C" ${RUN}_OTUs.uc > asvgroups.temp
cut -f2 asvgroups.temp | while read l;  do echo "otu$(($l + 1))"; done > otus.temp
cut -f3-9 asvgroups.temp > rest.temp
paste otus.temp rest.temp > ${RUN}_OTUs.clean.uc
rm *temp -f

# alternative: swarm OTUs
/nfs4/FW_HMSC/Levi_Lab/Hallerud/ESRF/Arthropods/swarm/bin/swarm -z ${RUN}_ASVs.fasta -d 9 -w ${RUN}_swarmOTUs.fasta -o ${RUN}_swarmOTUs.txt

# Mapping reads to each OTU- this is using a strict option (100% match) for ASVs and less strict option for OTUs
# Note that this method for generating an OTU table only works for OTUs defined via greedy clustering (i.e., as in Usearch and Vsearch), making an OTU table for other methods is outlined in the tutorial
#vsearch --search_exact ${LIB}.fasta --db #f	${LIB}_ASVs.fasta --otutabout ${LIB}_ASVs_mapped.txt

# set up for multithreading with vsearch_mapping.txt
# run PER LIBRARY (otherwise will take forever)
for f in `ls -1 ASVs/*ASVs.fasta`; do
	vsearch --usearch_global ${f} --db ${RUN}_OTUs.fasta --id 0.97 --otutabout ${RUN}_OTUs_mapped.txt
done

 
