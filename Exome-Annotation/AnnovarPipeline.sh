#peline for single family with PSAP

usage="AnnovarPipeline.sh -i INPUT_VCF.vcf -b TARGET_INTERVEL.bed -N 100"

while getopts i:v:N:H opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        b) BedFil="$OPTARG";;
	    N) NumJobs=100;;
        H) echo "$usage"; exit;;
	esac
done

if [[ ! -e "$InpFil" ]] || [[ ! -e "$BedFil" ]] ; then echo "Missing/Incorrect required arguments"; echo "$usage"; exit; fi

#SplitBed and VCF
TarLen=$(cat $BedFil | wc -l)
RemTar=$(( TarLen % NumJobs )) # get remainder of target file length and number of jobs
QuoTar=$(( TarLen / NumJobs )) # get quotient of target file length and number of jobs
SttLn=1
DivLen=0
echo $RemTar
echo $SttLn
for ((i=1; i <= $ArrNum; i++)); do
    SttLn=$(( SttLn + DivLen ))
    if [[ $i -le $RemTar ]]; then
        DivLen=$(( QuoTar + 1 ))
        else
        DivLen=$QuoTar
    fi
done

if [[ -z "$VcfNam" ]];then VcfNam=`basename $InpFil`; VcfNam=${VcfNam/.vcf.gz/}; fi # a name for the output files
if [[ -z $LogFil ]]; then LogFil=$VcfNam.GgVCF.log; fi # a name for the log file
VcfDir=$VcfNam.splitfiles; mkdir -p $VcfDir # Directory to output slices to
TmpDir=$VcfNam.GgVCF.tempdir; mkdir -p $TmpDir #temporary directory
TgtFil=$TmpDir/Range.$VcfNam.bed #exome capture range
tail -n+$SttLn $TgtBed | head -n $DivLen > $TgtFil #get exome capture range
echo "Target file line range: $SttLn - $(( $SttLn + $DivLen - 1 ))" >> $TmpLog





