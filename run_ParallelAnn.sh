#!/bin/bash

Ref=/home/local/users/jw/CUMC/Exome-pipeline-Jiayao/WES_Pipeline_References.b37.sh
AnnotateVCF=/home/local/users/jw/CUMC/Exome-pipeline-Jiayao/ExmVC.3.AnnotateVCF.sh
MergeVCF=/home/local/users/jw/CUMC/Exome-pipeline-Jiayao/ExmVC.2.MergeVCF.sh
source $Ref
InpFil=`readlink -f $1`
if [[ $InpFil != *.vcf.gz ]];then
	echo "UnCompressed $InpFil, Let's Compress it"
	bgzip $InpFil
	tabix -p vcf $InpFil.gz
	InpFil=$InpFil.gz
else
	echo "Compressed $InpFil"
fi
if [[ ! -f $InpFil.tbi ]];then
	tabix -p vcf $InpFil
fi
VcfNam=$(basename $InpFil|sed s/.gz$//|sed s/.vcf$//)
mkdir -p $VcfNam
cd $VcfNam
Contigs=`tabix -l $InpFil`
for contig in $Contigs
do
	#echo $contig
	#Modify the small vcf
	Header=Tmp."$contig"."$VcfNam".Header
	ContigFil=Tmp."$contig"."$VcfNam".vcf
	tabix -H $InpFil >$Header
	tabix $InpFil $contig |cat $Header - >$ContigFil
	rm $Header
	echo $ContigFil
	#nohup $AnnotateVCF -i $ContigFil -r $Ref -t 1&
done
parallel $AnnotateVCF -i Tmp.{}."$VcfNam".vcf -r $Ref ::: $Contigs
cd ../
echo "All Jobs finished."
VcfFil="$VcfNam.Annotated.vcf"
echo $VcfFil
SrtDir="Tmp.sort.$VcfNam"
mkdir -p $SrtDir
vcf-concat -p $VcfNam/*.$BUILD.vcf.gz | vcf-sort -t $SrtDir -c > $VcfFil
bgzip -f $VcfFil; tabix -f -p vcf $VcfFil.gz
rm -rf $VcfNam $SrtDir
