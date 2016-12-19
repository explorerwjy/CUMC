## Resource Directories
export EXOMPPLN="/home/local/ARCS/hz2408/pipeline/Exome_pipeline_scripts_GATKv3" # Directory containing pipeline shell scripts
export EXOMRES="/home/local/ARCS/hq2130/Exome_Seq/resources" # Directory containing resources/references for pipeline

#bwa,samtools,picard,vcftools,bcftools
export PATH="/home/local/ARCS/hz2408/bin/R-3.2.3/bin:${PATH}"
export PATH="/home/local/ARCS/hq2130/src/bwa/:${PATH}"
export PATH="/home/local/ARCS/hq2130/src/samtools/:${PATH}"
export PATH="/home/local/ARCS/hq2130/src/bcftools/:${PATH}"
export PATH="/home/local/ARCS/hq2130/src/htslib/:${PATH}"
export PATH="/home/local/ARCS/hq2130/src/picard-tools-1.139/:${PATH}"
export PATH="/home/local/ARCS/hz2408/bin/vcftools_0.1.13/bin/:${PATH}" 
export PERL5LIB="/home/local/ARCS/hz2408/bin/vcftools_0.1.13/perl/:${PATH}"
export PATH="/home/local/ARCS/hz2408/bin/Annovar:${PATH}"
export PATH="/home/local/ARCS/hq2130/src/ExpeDat/:${PATH}"

#jar files  and directories for software
GATKJAR="/home/local/ARCS/hq2130/src/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar" #Current GATK jar file
PICARD="/home/local/ARCS/hq2130/src/picard-tools-1.139/" #directory containing Picard jar files
SNPEFF="/home/local/ARCS/hq2130/src/snpEff/snpEff.jar" # Current snpEff jar file


## References
export BUILD="b37" # shorthand for build
export REF='$EXOMRES/hg19.fasta' #it is actually b37
export HAPMAP="$EXOMRES/hapmap_3.3.b37.vcf" # hapmap vcf from GATK
export INDEL="$EXOMRES/Mills_and_1000G_gold_standard.indels.b37.sites.vcf" # Gold standard INDEL reference from GATK
export TGVCF="$EXOMRES/1000G_omni2.5.b37.sites.vcf " 
export INDEL1KG="$EXOMRES/1000G_phase1.indels.b37.vcf" # INDEL reference from 1000 genomes
export DBSNP="$EXOMRES/dbsnp_138.b37.vcf" # dbSNP vcf from GATK
export ONEKG="$EXOMRES/1000G_phase1.snps.high_confidence.b37.vcf" # 1000 genome SNPs vcf
export ANNOVAR='/home/local/ARCS/hq2130/src/annovar' 
export ANNHDB='/home/local/ARCS/hz2408/resources/humandb/' #Location of annovar databases 
export HUMANREF="$EXOMRES/human_g1k_v37.fasta" # human 1000 genome assembly from GATK


# no use below 
export STHSH="$EXOMRES/b37/stampy_b37" # hash file for Stampy - omit ".sthash" extension for compatibility with Stampy
export STIDX="$EXOMRES/b37/stampy_b37" # genome index file for Stampy - omit ".stidx" extension for compatibility with Stampy


#GATK no-phone-home key
export ETKEY="$EXOMRES/ads2202_c2b2.columbia.edu.key"

#Capture Kit Target Files
export AgtV2="$EXOMRES/SureSelect_All_Exon_V2.b37.ordered.bed"
export AgtV4="$EXOMRES/SureSelect_All_Exon_V4_b37.ordered.bed"
export AgtV5="$EXOMRES/SureSelect_Human_All_Exon_V5_Covered.ordered.bed"
export AgtV5UTR="$EXOMRES/SureSelect_Human_All_Exon_V5_UTRs_Covered.ordered.bed"
export NbgV2="$EXOMRES/SeqCap_EZ_Exome_v2.hg19.targets.bed"
export NbgV3="$EXOMRES/SeqCap_EZ_Exome_v3.b37.targets.bed"
export IllTS="$EXOMRES/truseq_exome_targeted_regions.b37.ordered.bed"
export BigTgt="$EXOMRES/custom_intervals/FullIntervals.bed"
export TGTCODES="AgtV2:AgtV4:AgtV5:AgtV5UTR:NbgV2:NbgV3:IllTS:BigTgt"

