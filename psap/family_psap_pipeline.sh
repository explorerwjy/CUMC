#!/usr/bin/env bash

# $1 = vcf file, $2 = output file, $3 = ped file

PSAP_PATH=/home/local/users/jw/CUMC/
ANNOVAR_PATH=/home/local/users/jw/software_packages/annovar/ #INSERT PATH TO ANNOVAR DIRECTORY HERE eg. /scratch/dclab/annovar/
curpath=$(pwd)
echo $PWD
echo "PSAP path is "${PSAP_PATH}"psap/"
echo "ANNOVAR PATH is "${ANNOVAR_PATH}
#module load R
echo "Note: ped file requires header"
echo "Example: popScore_analysis_pipeline.sh FAM001.vcf FAM.OUTFILE FAM001.ped"
echo "Example: popScore_analysis_pipeline.sh FAM001.vcf FAM.OUTFILE FAM001.ped"
if [ $# -gt 0 ] && [ $1 == "-h" ]
then
	echo "arg1 =  VCF file"
	echo "arg2 = output file name"
	echo "arg3 = family pedigree file"
	echo "Example: popScore_analysis_pipeline.sh FAM001.vcf FAM.OUTFILE FAM001.ped"
	exit
fi

if [ $# == 3 ]
then
# Check that all required ANNOVAR annotation files are in the humandb directory
	MISSING=0
	for FILE in "hg19_ALL.sites.2014_09.txt" "hg19_cadd.txt" "hg19_esp6500si_all.txt" "hg19_snp137.txt" "hg19_wgEncodeGencodeBasicV19Mrna.fa" "hg19_wgEncodeGencodeBasicV19.txt" "hg19_mac63kFreq_ALL.txt"
	do
		if [ ! -f ${ANNOVAR_PATH}humandb/$FILE ]
		then
			MISSING=$(( $MISSING+1 ))
		fi
	done
# If any required ANNOVAR annotations are missing print an error message and then exit
	if [ $MISSING -gt 0 ]
	then
		echo "ERROR: Missing required ANNOVAR annotations.  Please run get_annovar_annos.sh prior to running this script."
		exit
	fi

# Extract and move to VCF file directory
	FILE_LOC=${1%/*.vcf} # Extract location of VCF file
	if [[ $FILE_LOC != *".vcf"* ]]
	then
		cd $FILE_LOC # Use location of  VCF file as working directory, this is where all output will be written
	fi
	echo $PWD
	VCF=${1##/*/} # Extract VCF file name
	echo $VCF
	OUTFILE=$2 # Name of output file (no directory should be included here)
	PED_FILE=$(readlink -f $3) # Name of pedigree file (directory should be included here)
    
	
# Convert vcf file to annovar file
  echo "$curpath/annotated/${OUTFILE}.avinput.hg19_multianno.txt"
#exit;  
if [ ! -e $curpath/annotated/${OUTFILE}.avinput.hg19_multianno.txt  ]; then    
     echo "PROGRESS: Converting VCF file to annovar input"
   echo "perl ${ANNOVAR_PATH}convert2annovar.pl -format vcf4old $VCF -outfile ${OUTFILE}.avinput -includeinfo"
      perl ${ANNOVAR_PATH}convert2annovar.pl -format vcf4old $VCF -outfile $curpath/${OUTFILE}.avinput -includeinfo
 # Write column names from VCF file to header file (will be used later)
        

# If there is no annotated directory create annotated directory
        if [ $(ls -d $PWD/*/ | grep -c -w "annotated") == 0 ]
        then
                echo "PROGRESS: Creating directory annotated/ to store annotation and analysis output"
                mkdir $curpath/annotated
        fi
# Annotate with ANNOVAR
	echo "PROGRESS: Annotating data with ANNOVAR"
	perl ${ANNOVAR_PATH}table_annovar.pl $curpath/${OUTFILE}.avinput -remove -outfile $curpath/annotated/${OUTFILE}.avinput ${ANNOVAR_PATH}humandb/ -buildver hg19 -protocol wgEncodeGencodeBasicV19,mac63kFreq_ALL,esp6500si_all,1000g2014sep_all,snp137,cadd -operation g,f,f,f,f,f -nastring NA -otherinfo -argument -separate,,,,,-otherinfo
else 
    echo "$curpath/annotated/${OUTFILE}.avinput.hg19_multianno.txt existed!"
fi
less $VCF|head -n 1500| grep '#' | tail -n 1 > $curpath/annotated/${OUTFILE}.avinput.header # Extract all VCF header lines.  Last line of VCF header contains column names.  Write last line of VCF header to .avinput.header file
# Annotate with PSAP (requires ped file)

# Annotate with PSAP (requires ped file)
	echo "bash ${PSAP_PATH}psap/annotate_PSAP.sh ${OUTFILE}.avinput $PED_FILE"
	
	bash ${PSAP_PATH}psap/annotate_PSAP.sh $curpath/annotated/${OUTFILE}.avinput $PED_FILE

else
	echo "ERROR: Incorrect number of arguments." $# "arguments provided"
	echo "Please provide the VCF file, output file name and family pedigree file.  Use the -h argument for help."
fi
