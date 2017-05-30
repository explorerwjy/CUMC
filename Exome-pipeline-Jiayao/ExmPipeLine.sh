#!/bin/bash
ProjectHome=
RefFil=$HOME/CUMC/Exome-pipeline-Jiayao/WES_Pipeline_References.b37.biocluster.sh
FastQC=
BWA=
MergeBAM=$HOME/CUMC/Exome-pipeline-Jiayao/ExmVC.2.MergeVCF.sh
DoC=
HapCaller=
JointGT=
VCFMerge=
VQSR=
Annovar=

BedFil=
FastQList=
RawBamList=
BamList=
GVCFList=
SplitedDir=
RawVCF=
VQSRVCF=
AnnovarVCF=


#HaploytypeCaller
NJob=`wc -l $BamList|cut -f 1 -d ' '`
seq $NJob | parallel -j 20 --eta $HapCaller -i $BamList -r $RefFil -t $BedFil -a {}

#Joint Genotyping
seq 20 | parallel -j 20 --eta sh $JointGT -i $GVCFList -r $RefFil -a {} -t $BedFil

#MergeVCF
bash $VCFMerge -i $SplitedDir -r $RefFil

#VQSR
bash $VQSR -i $RawVCF -r $RefFil

#Annovar
bash $Annovar -i $VQSRVCF -r $RefFil
