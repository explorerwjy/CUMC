#Pipeline for single family with PSAP
OriginalVCF=$1
FamList=$2
OutName=`basename $FamList|sed s/.list//g`
echo $OutName

step1="vcftools --gzvcf $OriginalVCF --recode --recode-INFO-all --out $OutName --keep $FamList"
echo $step1
$step1

step2="python /home/local/users/jw/CUMC/Exome-Filters-Jiayao/Clean.py -v $OutName.recode.vcf -o $OutName.vcf"
echo $step2
$step2
rm $OutName.recode.vcf

step3="python /home/local/users/jw/CUMC/Exome-Filters-Jiayao/ExmFilter.py -v $OutName.vcf -p $OutName.ped -o $OutName.RareCoding.vcf"
echo $step3
$step3
rm $OutName.vcf

step4="/home/local/users/jw/CUMC/psap/family_psap_pipeline.sh $OutName.RareCoding.vcf $OutName $OutName.ped"
echo $step4
$step4


step5="PSAP_Report.py -r ./annotated/$OutName.report.txt -p $OutName.ped -v $OutName.RareCoding.vcf"
echo $step5
$step5

