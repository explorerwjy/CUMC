#!/bin/bash
#$ -l mem=4G,time=1:: -cwd -N Primus
#set default arguments
usage="
 -t 1-<X> xExmCheckKinship_FINAL.sh -i <InputFile> -g <GenFil> -r <RltFil> -s <SexFil> -l <logfile> -H

     -i (required) - Input Table in correct format:
        Family <TAB> Indv1,Indv2,...IndvN (comma separate indiviudal IDs for all members of the family) ...[any further columns are ignored]
     -g (required) - output of plink IBD analysis (--genome)
     -s (required) - output of plink sex check analysis (--impute-sex)
     -H (flag) - echo this message and exit
"

KinFil=/ifs/scratch/c2b2/af_lab/ads2202/Exome_Seq/resources/kintab.txt
#get arguments
while getopts i:g:s:FBH opt; do
    case "$opt" in
        i) InpTab="$OPTARG";;
        g) GenFil="$OPTARG";; 
        s) SexFil="$OPTARG";;
        H) echo "$usage"; exit;;
    esac
done

KinFil=/ifs/scratch/c2b2/af_lab/ads2202/Exome_Seq/resources/kintab.txt
KinTab=${InpTab/.txt/}"_Checked.txt"


InpTab=`readlink -f $InpTab`
GenFil=`readlink -f $GenFil`
SexFil=`readlink -f $SexFil`


##use PRIMUS to generate pedigree images
FAM=`cut -f 1 $InpTab | head -n $SGE_TASK_ID | tail -n 1`
echo $FAM
if [[ ${FAM:0:1} == [0-9] ]]; then
    FAM=Family$FAM
fi
echo $FAM
mkdir PRIMUS.$FAM
cd PRIMUS.$FAM

SamList=`head -n $SGE_TASK_ID $InpTab | tail -n 1 | cut -f 2 | sed s/,/" "/g`
for i in $SamList; do echo $i >> SamList.$FAM ; done

awk 'NR==FNR{a[$0];next}NF<3||($1 in a)&&($3 in a)' SamList.$FAM $GenFil > $FAM.genome 
~/scratch/src/PRIMUS_v1.8.0/bin/run_PRIMUS.pl -p $FAM.genome --sexes FILE=$SexFil SEX=4
for j in $(find | grep ps$); do
     NetNum=${j##*_}
     NetNum=${NetNum/.ps/}
     OutFil=Fam$FAM.PrimusPedigree.$NetNum
     LEN1=`grep -n "Page: 1" $j | cut -f 1 -d ":"`
     LEN2=`grep -n "Page: 2" $j | cut -f 1 -d ":"`
     LEN2=$(( LEN2 + 1 ))
     head -n $LEN1 $j > Temp.ps
     tail -n +$LEN2 $j | sed s/%%Pages" "2/%%Pages" "1/ | sed s/PRIMUS-Network/$FAM-Network/ >> Temp.ps
     ps2pdf Temp.ps Temp.pdf
     convert Temp.pdf $OutFil.png
done
rm -rf Temp*
