#!/bin/bash
# Download ANNOVAR Database for hg19

dbToDownload=$1

AnnovarDirectory='/home/local/users/jw/software_packages/annovar'



# Annotation database
refGene=refGene #FASTA sequences for all annotated transcripts in RefSeq Gene. Date: 20151211
MultiScore=dbnsfp30a #whole-exome SIFT, PolyPhen2 HDIV, PolyPhen2 HVAR, LRT, MutationTaster, MutationAssessor, FATHMM, MetaSVM, MetaLR, VEST, CADD, GERP++, DANN, fitCons, PhyloP and SiPhy scores from dbNSFP version 3.0a. Date: 20151015
dbnsfp=dbnsfp31a_interpro #protein domain for variants. Date: 20151219
dbscSNV=dbscsnv11 #dbscSNV version 1.1 for splice site prediction by AdaBoost and Random Forest. Date: 20151218
COSMIC=cosmic70 #COSMIC database version 70. Date: 20140911
ExAC=exac03 #ExAC 65000 exome allele frequency data for ALL, AFR (African), AMR (Admixed American), EAS (East Asian), FIN (Finnish), NFE (Non-finnish European), OTH (other), SAS (South Asian)). version 0.3. Left normalization done. Date: 20151129
KAVIAR=kaviar_20150923 #170 million Known VARiants from 13K genomes and 64K exomes in 34 projects. Date: 20151203
HRCRL=hrcr1 #40 million variants from 32K samples in haplotype reference consortium. Date: 20151203
OneKG=1000g2015aug #The 1000G team fixed a bug in chrX frequency calculation. Based on 201508 collection v5b (based on 201305 alignment). Date: 20150824
MCAP=mcap #M-CAP scores for non-synonymous variants. Date: 20161104
dbSNP=avsnp147 #dbSNP147 with allelic splitting and left-normalization. Date: 20160606
ICGC=icgc21 #International Cancer Genome Consortium version 21. Date: 20160622
CLINVAR=clinvar_20160302 #Clinvar with separate columns (CLINSIG CLNDBN CLNACC CLNDSDB CLNDSDBID). Date: 20160303
MITIMPACT=mitimpact24 #pathogenicity predictions of human mitochondrial missense variants. Date: 20160123
CADD=cadd13 #CADD version 1.3. Date: 20160607
CADD10=cadd13gt10 #CADD version 1.3 score>10. Date: 20160621
CADD20=cadd13gt20 #CADD version 1.3 score>20. Date: 20160621
FATHMM=fathmm #whole-genome FATHMM_coding and FATHMM_noncoding scores (noncoding and coding scores in the 2015 version was reversed). Date: 20160315
GWAVA=gwava #whole genome GWAVA_region_score GWAVA_tss_score GWAVA_unmatched_score. Date: 20150623
EIGEN=eigen #whole-genome Eigen scores. Date: 20160330
ESPall=esp6500siv2_all #alternative allele frequency in All subjects in the NHLBI-ESP project with 6500 exomes, including the indel calls and the chrY calls. This is lifted over from hg19 by myself. Date: 20141222
ESPaa=esp6500siv2_aa #alternative allele frequency in African American subjects in the NHLBI-ESP project with 6500 exomes, including the indel calls and the chrY calls. This is lifted over from hg19 by myself. Date: 20141222
ESPea=esp6500siv2_ea #alternative allele frequency in European American subjects in the NHLBI-ESP project with 6500 exomes, including the indel calls and the chrY calls. This is lifted over from hg19 by myself. Date: 20141222


cd $AnnovarDirectory

case $dbToDownload in
    1)
        nam="Download: refseq hg19 gene reference"
        cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar $refGene humandb/"
        ;;
    2)
		nam="Download: Multi-Scores"
        cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar $MultiScore humandb/"
        ;;
    3)
		nam="Download: dbSFP reference"
        cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar $dbsfp humandb/"
        ;;
    4)
		nam="Download: dbscSNV"
        cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar $dbscSNV humandb/"
        ;;
    5)
		nam="Download: ESP alternative allele frequency - all"
        cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar $ESPall humandb/"
        ;;
    6)
		nam="Download: ESP alternative allele frequency - African Americans"
        cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar $ESPaa humandb/"
        ;;
    7)
		nam="Download: SP alternative allele frequency - European Americans"
        cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar $ESPea humandb/"
        ;;
    8)
		nam="Download: superdups "
        cmd="perl annotate_variation.pl -downdb -buildver hg19 genomicSuperDups humandb/"
        ;;
    9)
		nam="Download: Cadd top 10% "
        cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar $CADD10 humandb/"
        ;;
    10)
		nam="Download: Cadd top 20% "
        cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar $CADD20 humandb/"
        ;;
    11)
		nam="Download: exac03"
        cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar $ExAC humandb/"
        ;;
    12)
		nam="Download: Cadd full "
        cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar $CADD humandb/"
        ;;
    13)
		nam="Download: CLINVAR database with Variant Clinical Significance  "
        cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar $CLINVAR humandb/"
        ;;
    14)
		nam="Download: COSMIC "
        cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar $COSMIC humandb/"
        ;;
	15)
		nam="Download: MCAP "
		cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar $MCAP humandb/"
		;;
	16)
		nam="Download: KAVIAR "
		cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar $KAVIAR humandb/"
		;;
	17)
		nam="Download: HRCRL "
		cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar $HRCRL humandb/"
		;;
	18)
		nam="Download: 1KG "
		cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar $OneKG humandb/"
		;;
	19)
		nam="Download: dbSNP "
		cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar $dbSNP humandb/"
		;;
	20)
		nam="Download: ICGC "
		cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar $ICGC humandb/"
		;;
	21)
		nam="Download: MITIMPACT "
		cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar $IMITIMPACT humandb/"
		;;
	22)
		nam="Download: FATHMM "
		cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar $FATHMM humandb/"
		;;
	23)
		nam="Download: GWAVA "
		cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar $GWAVA humandb/"
		;;
	24)
		nam="Download: EIGEN "
		cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar $EIGEN humandb/"
		;;
	25) 
		nam="Download: gnomad_exome"
		cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar gnomad_exome humandb/"
		;;
	26)	
		nam="Download: gnomad_exome"
		cmd="perl annotate_variation.pl -downdb -buildver hg38 -webfrom annovar gnomad_exome humandb/"
		;;
	27)
		nam="Download: gnomad_genome"
		cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar gnomad_genome humandb/"
		;;
	28)
		nam="Download: gnomad_genome"
		cmd="perl annotate_variation.pl -downdb -buildver hg38 -webfrom annovar gnomad_genome humandb/"
		;;
esac

echo "Start $nam - `date`"
echo $cmd
eval $cmd
if [[ $? == 0 ]]; then
    echo "Finish $nam - `date`"
else
    echo "Error during "${nam/load:/loading}
fi
