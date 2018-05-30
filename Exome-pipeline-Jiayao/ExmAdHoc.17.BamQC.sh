#!/bin/bash
#$ -S /bin/bash
#$ -j y
#$ -N BAMQC 
#$ -l h_rt=12:00:00
#$ -l h_vmem=16G
#$ -cwd

# This script do some QC on bam files. Including flagstats, idxstats, duplications, insertion sizes and gc bias.

#set default arguments
usage="
-t 1-{number of fastq files] ExmAln.1a.Align_Fastq_to_Bam_with_BWAmem.sh -i <InputFile> -r <reference_file> -t <target intervals file> -l <logfile> -PH

     -i (required) - Table containing the path to the fastq file and the RG read header
     -r (required) - shell file containing variables with locations of reference files and resource directories (WES_Pipeline_References.b37.sh)
     -l (optional) - Log file
     -t (optional) - Exome capture kit targets or other genomic intervals bed file (must end .bed for GATK compatability); this file is required if calling the pipeline but otherwise can be omitted
     -P (flag) - Initiate exome analysis pipeline after completion of script
     -F (flag) - Fix mis-encoded base quality scores - Rewrite the qulaity scores from Illumina 1.5+ to Illumina 1.8+
     -H (flag) - echo this message and exit
"

PipeLine="false"

#get arguments
while getopts i:r:a:l:PFH opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        r) RefFil="$OPTARG";; 
        a) ArrNum="$OPTARG";; 
		l) LogFil="$OPTARG";;
        P) PipeLine="true";;
        F) FixMisencoded="true";;
        H) echo "$usage"; exit;;
    esac
done

#check all required paramaters present
if [[ ! -e "$InpFil" ]] || [[ ! -e "$RefFil" ]]; then echo "Missing/Incorrect required arguments"; echo "$usage"; exit; fi

#Call the RefFil to load variables
RefFil=`readlink -f $RefFil`
source $RefFil 

#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func"

#set local variables
if [[ -z "${ArrNum}" ]]
then
    ArrNum=$SGE_TASK_ID
fi

InpFil=`readlink -f $InpFil`  # resolve input file path
BAM=`readlink -f $(tail -n+$ArrNum $InpFil | head -n 1 | cut -f1)`
BamNam=$(basename $BAM | sed s/.bam// ) # a name for the output files - basically the original file name
FlgStat=$BamNam.flagstat #output file for bam flag stats
IdxStat=$BamNam.idxstats #output file for bam index stats
if [[ -z "$LogFil" ]]; then LogFil=$BamNam.FqB.log; fi # a name for the log file
TmpLog=$BamNam.FqB.temp.log #temporary log file
TmpDir=$BamNam.FqB.tempdir; mkdir -p $TmpDir #temporary directory

#start log
ProcessName="Collect BAMQC Metrics with Samtools&Picard"
#funcWriteStartLog
echo " Build of reference files: "$BUILD >> $TmpLog
echo "----------------------------------------------------------------" >> $TmpLog

#Get flagstat
StepName="Output flag stats using Samtools"
StepCmd="samtools flagstat $BAM > $FlgStat"
#funcRunStep

#get index stats
StepName="Output idx stats using Samtools"
StepCmd="samtools idxstats $BAM > $IdxStat"
#funcRunStep

StepName="LibraryComplexity"
StepCmd="java -Xmx4G -XX:ParallelGCThreads=1 -Djava.io.tmpdir=$TmpDir -jar $PICARD EstimateLibraryComplexity 
 INPUT=$BAM
 OUTPUT=$BamNam.est_lib_complex_metrics.txt
 2>>$TmpLog"
#funcRunStep

StepName="CollectInsertSize"
StepCmd="java -Xmx4G -XX:ParallelGCThresdfsdfsads=1 -Djava.io.tmpdir=$TmpDir -jar $PICARD CollectInsertSizeMetrics 
 INPUT=$BAM
 OUTPUT=$BamNam.insert_size_metrics.txt
 H=$BamNam.insert_size_histogram.pdf
 M=0.05
 DEVIATIONS=10.0
 INCLUDE_DUPLICATES=falsegc_bias_metrics
 ASSUME_SORTED=true
 2>>$TmpLog"
#funcRunStep

#======================================================================================================
# Notice This Step Require Exactly The Same Reference Genome as Bam Files. You May Want to Change
# $REF Here.
#======================================================================================================
StepName="GC Bias"
StepCmd="java -Xmx4G -XX:ParallelGCThreads=1 -Djava.io.tmpdir=$TmpDir -jar $PICARD CollectGcBiasMetrics 
 INPUT=$BAM
 OUTPUT=$BamNam.gc_bias_metrics.txt
 CHART=$BamNam.gc_bias_metrics.pdf
 S=$BamNam.summary_metrics.txt
 R=$REF
 2>>$TmpLog"
funcRunStep

#R=/home/local/users/jw/Genetics_Projects/SPARK/Regeneron_NovoSeq_Freeze_One/resouces/genome.fa
#End Log
funcWriteEndLog
echo "Done BAM QC "$BAM
