#!/bin/bash
ProjectHome=`pwd`
RefFil=$HOME/CUMC/Exome-pipeline-Jiayao/WES_Pipeline_References.b37.sh
FastQC=
BWA=
MergeBAM=$HOME/CUMC/Exome-pipeline-Jiayao/ExmVC.2.MergeVCF.sh
DoC=
HapCaller=/home/local/users/jw/CUMC/Exome-pipeline-Jiayao/ExmAln.2.HaplotypeCaller_GVCFmode.sh
JointGT=/home/local/users/jw/CUMC/Exome-pipeline-Jiayao/ExmVC.1hc.GenotypeGVCFs.sh
VCFMerge=/home/local/users/jw/CUMC/Exome-pipeline-Jiayao/ExmVC.2.MergeVCF.sh
VQSR=/home/local/users/jw/CUMC/Exome-pipeline-Jiayao/ExmVC.4.RecalibrateVariantQuality.sh
Annovar=/home/local/users/jw/CUMC/Exome-pipeline-Jiayao/Test.AnnotateVCF_direct.sh

BedFil=/home/local/users/jw/resources/references/b37/CaptureKitBeds/CCDS_hg19.bed
FastQList=
RawBamList=
BamList=/home/local/users/jw/tmp/Jhe/bam.list
GVCFList=
SplitedDir=
RawVCF=
VQSRVCF=
AnnovarVCF=


#HaploytypeCaller
#NJob=`wc -l $BamList|cut -f 1 -d ' '`
#echo $NJob
#seq $NJob | parallel -j 20 --eta $HapCaller -i $BamList -r $RefFil -t $BedFil -a {}

#Joint Genotyping
#AC=/home/local/users/jw/tmp/Jhe/AC.gvcf.list
#BD=/home/local/users/jw/tmp/Jhe/BD.gvcf.list
#CUAC=/home/local/users/jw/tmp/Jhe/CUAC.gvcf.list
#seq 20 | parallel -j 20 --eta sh $JointGT -i $AC -r $RefFil -a {} -j 20 -t $BedFil
#seq 20 | parallel -j 20 --eta sh $JointGT -i $BD -r $RefFil -a {} -j 20 -t $BedFil
#seq 20 | parallel -j 20 --eta sh $JointGT -i $CUAC -r $RefFil -a {} -j 20 -t $BedFil

#MergeVCF
#SplitedDir=*.gvcf.splitfiles
#echo $SplitedDir
#parallel -j 3 --eta $VCFMerge -i {} -r $RefFil ::: $SplitedDir

#VQSR
#RawVcf=*.rawvariants.vcf.gz 
#RawVcf=CUAC.rawvariants.vcf.gz
#echo $RawVcf
#parallel -j 1 --eta $VQSR -i {} -r $RefFil ::: $RawVcf

#Annovar
InpFils=*recalibrated.vcf
parallel -j 3 --eta $Annovar -i {} -r $RefFil ::: $InpFils
