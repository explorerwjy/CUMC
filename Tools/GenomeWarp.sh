#!/bin/bash

genomewarp=/home/local/users/jw/software_packages/genomewarp/target/verilylifesciences-genomewarp-1.0.0-runnable.jar

# INPUTS
queryvcf=
targetvcf=
targetbed=

# QUERY TO TARGET
chain=/home/local/users/jw/resources/reference_genomes/LiftOver/hg38ToHg19.over.chain
querybed=/home/local/users/jw/resources/reference_genomes/GRCh38/HG38.WGS.bed
queryfasta=/home/local/users/jw/resources/reference_genomes/GRCh38/GRCh38.no_alt_analysis_set.fa
targetfasta=/home/local/users/jw/resources/reference_genomes/hg19/hg19.fasta
workdir=`pwd`

java -jar $genomewarp \
	--lift_over_chain_path "${chain}" \
	--raw_query_vcf "${queryvcf}" \
	--raw_query_bed "${querybed}" \
	--ref_query_fasta "${queryfasta}" \
	--ref_target_fasta "${targetfasta}" \
	--work_dir "${workdir}" \
	--output_variants_file "${targetvcf}" \
	--output_regions_file "${targetbed}"
