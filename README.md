# Bioinformatics Pipeline for Arthropod Metabarcoding
## Input Data
Bulk arthropod samples were collected using malaise traps, pitfall traps, and blue vane traps with a 1-week sampling period. Samples were transferred into 95% ethanol on the day of collection and transferred to the lab. For each sample, arthropods were manually sorted from debris and placed in 10-ml or 50-ml tubes with ethanol for storage. Prior to extraction, samples were dried in a fume hood and weighed for biomass. Half of the pitfall traps were pooled per site, then DNA was extracted from each sample (or pool) following a non-destructive lysis protocol. Metabarcoding for each DNA extract proceeded using the primers BF3 (5'-CCHGAYATRGCHTTYCCHCG-3') + BR2 (5'-TCDGGRTGNCCRAARAAYCA-3') to amplify ~458 bp of the COI Folmer region, with 3 PCR replicates per extract. Metabarcoding libraries were constructed using twin-tagged primers with unique 8-bp barcodes, and each library was pooled equimolar after PCR. Libraries were indexed and Illumina adapters were added using the NEBNext Ultra II library prep kits, then quantified and pooled for sequencing. Libraries were sequenced on an Illumina NextSeq 2000 P2 lane for 2x300 bp, with sequencing allocation doubled for libraries containing high diversity samples (i.e., pooled pitfalls & malaise traps) relative to libraries containing lower diversity samples (i.e., single pitfall and blue vane trap samples).

## Basic Bioinformatic Steps
0. Initial assessment of library quality:
   A) Assemble counts of amplicon sequence variants of all sequences containing primers within each library.
   B) Quantify tag-jumping and mis-assignment rates for reads within each tag combination, per library.
1. Demultiplex samples, requiring both barcodes and primers to be present, with no errors allowed.
2. Merge paired-end sequences.
3. Run quality filtering, denoising, length filtering, and chimera removal to create amplicon sequence variants (ASVs) per library.
4. Merge ASVs across libraries and cluster into project-level molecular operational taxonomic units (mOTUs) at a 97% similarity threshold.
5. Assign taxonomy to each OTU using the Barcode of Life (BOLD) & NCBI reference databases.
6. Post-filtering of sequences to remove chimeras, NUMTs, off-target, and poorly assigned sequences.

## Pipeline setup
The pipeline depends on the run directory being set up as follows:
- 0_raw: contains raw gzipped FASTQ files from sequencing facility. Files must follow the convention *`LIB_INDEX`*_R`1|2`_001.fastq.gz
- barcodes: contains barcode information for each library in 12-column TSV format

## Running the Pipeline
0. A) Assemble top amplicon sequence variants (with primers) per library.
```
#!/bin/bash
# usage: 0_make_library_ASVs.sh <type (PF | MT)> <lib#> <threads>
0_make_library_ASVs.sh PF 01 10
```  
0. B) Quantify tag-jumping and mis-assignment rates per library. Mapping to the top ASVs requires the previous step.
```
#!/bin/bash
# usage: 1_check_tag_jumping.sh <type> <lib#> <threads>
1_check tag_jumping.sh PF 01 10
```
1. Demultiplex all samples in run using cutadapt.
```
#!/bin/bash
# usage: 2_prep_cutadapt.sh <type> <lib#>
2_prep_cutadapt.sh PF 01
```
2. Merging read pairs using PEAR.
```
#!/bin/bash
# usage: 3_prep_pear <type> <lib>
3_prep_pear.sh PF 01
```
3. Denoising and filtering using vsearch, and chimera removal using dada2.
```
#!/bin/bash (requires R path to be set within script)
# usage: 4_denoise_library.sh <lib> <threads>
4_denoise_library.sh ESFI_PF_P01 10
```
4. OTU clustering across the run using vsearch.
Make sure to increase memory allocation for this step!
```
#!/bin/bash
# usage: 5_run_vsearch.sh <run> <threads>
5_run_vsearch.sh NextSeq676 10
# mapping is slow- this will output library-level commands under ${RUN}/3_clustering/run_mapping.txt
# that can be submitted for parallel processing
# once mapping is complete, merge OTUtabs across libraries:
/local/cqls/R/R-4.4.2/Rscript ${SCRIPT_DIR}/merge_otutabs.R ${RUN}
# this will create a run-level OTU table with PCR reps in columns, 97% OTUs in rows
```
5. Merging OTUs across runs.
```
# merge across runs (not the most elegant merge method but it works)
cat NextSeq*/3_clustering/ESFI*_OTUs.fasta > MERGED_OTUs.fasta
# re-cluster 
vsearch --cluster_size MERGED_OTUs.fasta --id 0.97 --centroids MERGED_97_OTUs.fasta --sizein --sizeout --uc MERGED_97_OTUs.uc --aln MERGED_97_OTUs.aln --msa MERGED_97_OTUs.msa
```
6. Taxonomic assignment with blast & dada2.
NOTE: dada2 requires a large memory allocation to run because it is not parallelized!
```
# run BLAST assignment
6_BLAST_assignment.sh MERGED_97 <threads>
# run dada2 assignment- exports confidence scores
/local/cqls/R/R-4.4.2/Rscript 6_dada2_taxonAssign.R MERGED_97_OTUs.fasta MERGED_97_OTUs_BOLD_dada2Assign.csv
```
7. Additional chimera filtering with LULU.
....

Note: For each of these, commands can be written to a file and submitted as a job to enable multi-threading, e.g.:
```
#!/bin/bash
for p in {01..18}; do
	2_prep_cutadapt.txt PF ${p}
done > cutadapt_commands.txt

cat cutadapt_commands.txt | hqsub -t array -q fw -r cutadapt_demult_log
```
