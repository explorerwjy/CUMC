BamLocation=$1
Dest=$2
Dest=$(readlink -f $Dest)
echo $Dest
if [[ -z $Dest ]]
then
	echo "No Destination"
	exit
fi

for line in `less $BamLocation`
do
	echo $line
	DirNam=$(dirname $line)
	BamNam=$(basename $line|sed s/.bam//)
	bai="${DirNam}/${BamNam}"
	OutNam=$(samtools view -H $line|grep -P  "RG" |head -n1|cut -d " " -f1|grep -P -o "SM:.+$" |uniq|sed s/SM://)
	echo $OutNam
	BAMLK=${OutNam}.bam
	BAILK=${OutNam}.bai
	#echo $BAMLK $BAILK
	ln -s $line ${Dest}/$BAMLK
	if [ -e $bai.bai ] 
	then
		bai=$bai.bai	
	elif [ -e $bai.bam.bai ] 
	then
		bai=$bai.bam.bai
	else
		echo "Cant Decide BAI file"
	fi
	ln -s $bai ${Dest}/$BAILK
done
