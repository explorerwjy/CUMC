## Resource Directories
export EXOMPPLN="/home/yufengshen/CUMC/Exome-pipeline-Jiayao"
export EXOMRES="/home/yufengshen/resources" # Directory containing resources/references for pipeline

#jar files  and directories for software
GATKJAR="//home/yufengshen/bin/GenomeAnalysisTK-3.2.jar" #Current GATK jar file
PICARD="/home/yufengshen/bin/picard.jar" #directory containing Picard jar files
SNPEFF="/home/yufengshen/bin/snpEff.jar" # Current snpEff jar file

## References
export BUILD="b37" # shorthand for build
export REF="$EXOMRES/reference_genomes/hg19/hg19.fasta" # human 1000 genome assembly from GATK
export HAPMAP="$EXOMRES/references/b37/hapmap_3.3.b37.vcf" # hapmap vcf from GATK
export INDEL="$EXOMRES/references/b37/Mills_and_1000G_gold_standard.indels.b37.sites.vcf" # Gold standard INDEL reference from GATK
export TGVCF="$EXOMRES/references/b37/1000G_omni2.5.b37.sites.vcf" 
export INDEL1KG="$EXOMRES/references/b37/1000G_phase1.indels.b37.vcf" # INDEL reference from 1000 genomes
export DBSNP="$EXOMRES/references/b37/dbsnp_138.b37.vcf" # dbSNP vcf from GATK
export ONEKG="$EXOMRES/references/b37/1000G_phase1.snps.high_confidence.b37.vcf" # 1000 genome SNPs vcf
export ANNOVAR='/home/local/ARCS/hq2130/src/annovar'
export ANNHDB='/home/local/ARCS/hz2408/resources/humandb/' #Location of annovar databases
export HUMANREF="$EXOMRES/human_g1k_v37.fasta" # human 1000 genome assembly from GATK

export STHSH="$EXOMRES/references/b37/stampy_b37" # hash file for Stampy - omit ".sthash" extension for compatibility with Stampy
export STIDX="$EXOMRES/references/b37/stampy_b37" # genome index file for Stampy - omit ".stidx" extension for compatibility with Stampy
#GATK no-phone-home key
export ETKEY="$EXOMRES/ads2202_c2b2.columbia.edu.key"

#Capture Kit Target Files
export xgen="/home/yufengshen/resources/references/b37/CaptureKitBeds/xgen-exome-research-panel-targets.bed"
