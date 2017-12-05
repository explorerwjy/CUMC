#!/bin/bash
RefFil=$HOME/CUMC/Exome-pipeline-Jiayao/WES_Pipeline_References.b37.sh
FastQC=
BWA=
MergeBAM=$HOME/CUMC/Exome-pipeline-Jiayao/ExmVC.2.MergeVCF.sh
DoC=$HOME/CUMC/Exome-pipeline-Jiayao/ExmAln.8a.DepthofCoverage.sh
HapCaller=$HOME/CUMC/Exome-pipeline-Jiayao/ExmAln.2.HaplotypeCaller_GVCFmode.sh
JointGT=$HOME/CUMC/Exome-pipeline-Jiayao/ExmVC.1hc.GenotypeGVCFs.sh
VCFMerge=$HOME/CUMC/Exome-pipeline-Jiayao/ExmVC.2.MergeVCF.sh
VQSR=$HOME/CUMC/Exome-pipeline-Jiayao/ExmVC.4.RecalibrateVariantQuality.sh
Annovar=$HOME/CUMC/Exome-pipeline-Jiayao/Test.AnnotateVCF_direct.sh

ProjectHome=
ProjectName=
BedFil=
FastQList=
RawBamList=
BamList=
GVCFList=
SplitedDir=
RawVCF=
VQSRVCF=
AnnovarVCF=


#==========================================================================================================
#HaploytypeCaller
GVCF=${ProjectHome}/GVCF
mkdir -p $GVCF
cd $GVCF
NJob=`wc -l $BamList|cut -f 1 -d ' '`
echo $NJob
seq $NJob | parallel -j 20 --eta $HapCaller -i $BamList -r $RefFil -t $BedFil -a {}
find `pwd` -name '*.g.vcf.gz' > $GVCFList
#==========================================================================================================

mkdir -p ${ProjectHome}/JointGenotyping
##==========================================================================================================
##Joint Genotyping
cd ${ProjectHome}/JointGenotyping
NUM_JOB=30
seq $NUM_JOB | parallel -j $NUM_JOB --eta sh $JointGT -i $GVCFList -r $RefFil -a {} -j $NUM_JOB -t $BedFil -n $ProjectName
##==========================================================================================================


#==========================================================================================================
#DoC
DoC_DIR=${ProjectHome}/DoC
mkdir -p $DoC_DIR
cd $DoC_DIR
NJob=`wc -l $BamList|cut -f 1 -d ' '`
echo $NJob
nohup seq $NJob | parallel -j 20 --eta $DoC -i $BamList -r $RefFil -t $BedFil -a {} &
#==========================================================================================================

#==========================================================================================================
#MergeVCF
cd ${ProjectHome}/JointGenotyping
echo $SplitedDir
$VCFMerge -i $SplitedDir -r $RefFil 
#==========================================================================================================

#==========================================================================================================
#VQSR
cd ${ProjectHome}/JointGenotyping
echo $RawVcf
$VQSR -i $RawVCF -r $RefFil 
#==========================================================================================================

#==========================================================================================================
#Annovar
cd ${ProjectHome}/JointGenotyping
if [ -e $VQSRVCF ]
then
	echo $VQSRVCF
	echo $Annovar -i $VQSRVCF -r $RefFil 
	$Annovar -i $VQSRVCF -r $RefFil 
else
	$Annovar -i $RawVCF -r $RefFil 
fi
#==========================================================================================================
