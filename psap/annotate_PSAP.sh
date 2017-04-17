# $1 = annovar input file name (eg. FILE.avinput)
# $2 = PED file name (eg. FAMILY.ped)

#module load R
#PSAP_PATH=/home/local/ARCS/nz2274/Application/ #INSERT PATH TO PSAP DIRECTORY HERE eg. "/scratch/dclab/"
PSAP_PATH=/home/local/users/jw/CUMC/

if [ $# == 2 ]
then  
# annotate with popstat (requires ped file)
	echo "Annotating data with population model"
	echo "Rscript ${PSAP_PATH}psap/RScripts/generic_apply_popStat.R ${1}.avinput ${2} $PSAP_PATH"

	if [ ! -s ${1}_popStat.txt ]; then
		Rscript ${PSAP_PATH}psap/RScripts/generic_apply_popStat.R ${1}.avinput ${2} $PSAP_PATH	
	fi   
# perform family analysis and validate
	echo "Performing family analysis and candidate validation"
	echo "Rscript ${PSAP_PATH}psap/RScripts/generic_candidate_analysis.R ${1}.avinput ${2} $PSAP_PATH"
	if [ ! -s ${1}.report.txt ]; then
		Rscript ${PSAP_PATH}psap/RScripts/generic_candidate_analysis.R ${1}.avinput ${2} $PSAP_PATH
	fi
	echo "Done.  All resutls can be found in the path ${PWD}/annotated/"
else
	echo $#
	echo "Incorrect number of arguments"
	echo "Please provide ANNOVAR input file name and family pedigree file"
fi
