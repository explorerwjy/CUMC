#!/bin/bash
#$ -cwd -l mem=1G,time=1:: -N VCFGetFreq

usage="xVCFGetFrequencies.sh -v <vcf file> -s <Sample List File> -o <Output Name> -H <this message>"

#get arguments
while getopts v:s:o:H opt; do
    case "$opt" in
        v) VcfFil="$OPTARG";;
        s) SamFil="$OPTARG";;
        o) OutNam="$OPTARG";;
        H) echo "$usage"; exit;;
    esac
done

VcfFil=`readlink -f $VcfFil`
SamFil=`readlink -f $SamFil`
LogFil=$OutNam.log

echo "`date`" > $LogFil
echo "Variant frequency acquistion:" >> $LogFil
echo -e "\tVCF file: $VcfFil" >> $LogFil
if [[ $SamFil ]]; then 
    echo -e "\tSample List:" >> $LogFil
    cat $SamFil >> $LogFil


CMD="vcftools --vcf $VcfFil --freq --out $OutNam"
if [[ $SamFil ]]; then 
    CMD=$CMD" --keep $SamFil"
fi


