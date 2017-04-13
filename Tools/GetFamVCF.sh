
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

step3="python /home/local/users/jw/CUMC/Exome-Filters-Jiayao/Filters.py -v $OutName.vcf -o $OutName.RareCoding.vcf"
echo $step3
$step3
