
PipePath="/home/local/ARCS/nz2274/Pipeline/NA_script/"
source $PipePath/mappability_generate.sh
function annotatedwithMappability {
   len=$1
   fMAP=$2
   fvcf=$3
   fbname=$(basename $fvcf)
   fbname=${fbname%.vcf.gz}
   fbname=${fbname%.vcf}
   fout="${fbname}.mappability.vcf.gz"
   echo $fMAP

  if [ ! -f $fMAP ]; then
         Redo_mappability $len
   fi
   keywd="Mappability"
   des="Mappability is obtained in file $fMAP"
   fhead="header_add.txt"

   echo "##INFO=<ID=$keywd,Number=1,Type=String,Description=\"$des\">" > $fhead
   bcftools annotate -a $fMAP  -c CHROM,FROM,TO,Mappability  -h  $fhead  $fvcf |bgzip -c > $fout
}


echo "input format annotatedwithMappability len *bed.gz  *.vcf.gz"
echo " bash $PipePath/mappability.annotation.sh 152 /home/local/ARCS/nz2274/Resources/mappability/hg19_152bp_mappability.bed.gz  /home/local/ARCS/nz2274/PAH/PAH_2017/GVCF/Rare.coding.PAH_new_old_PCGC_4650.vcf.gz"
annotatedwithMappability   $1 $2 $3
