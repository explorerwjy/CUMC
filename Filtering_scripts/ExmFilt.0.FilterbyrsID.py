#!/usr/bin/env python
#$ -cwd -l mem=2G,time=:15: -N AAFFilt

# The purpose of this script is to filter a VCF file by minor/alternate allele frequency as provided by 1KG and GO-ESP
#    -v/--vcf     <required>    The script requires an input vcf
#    -o/--out     <required>    The use should specify a base name for the output files.The script outputs the filtered results as a vcf file. The script also outputs a log file. 
#    -m/--maf     <optional>    Minor allele frequency for filtering. Default is 0.01
#    -G/--greater <flag>        Filter for variants with maf greater than or equal to the filter level. Default is less than or equal to.
#    -W/--within  <flag>        Filter for allele frequency within the cohort. Default is just 1KG and ESP frequencies.


from optparse import OptionParser
import gzip
import os

parser = OptionParser()
# Basic Input Files, required
parser.add_option("-v", "--vcf", dest="VCFfile",help="input VCF file", metavar="VCFfile")
parser.add_option("-o", "--output", dest="OutputFileName",help="user specified name of file for output to be written", metavar="OutputFileName")

(options, args) = parser.parse_args()

FilTyp=str(options.VCFfile)
FilTyp=FilTyp.split(".")
FilTyp=str(FilTyp[-1])
if FilTyp == 'gz':
    VCF=gzip.open(options.VCFfile,'r')
elif FilTyp == 'vcf':
    VCF=open(options.VCFfile,'r')
else:
    print "Incorrect vcf file type"
    sys.exit(1)

BaseName=str(options.OutputFileName)
#open input and output files
VcfOutputFilename=BaseName+'.filter.aaf.vcf'
LogOutputFilename=BaseName+'.filter.aaf.log'
Outvcf=open(VcfOutputFilename,'w')
Outlog=open(LogOutputFilename,'w')

with open('/home/local/ARCS/hq2130/Exome_Seq/scripts/Filtering_scripts/All_HapMap.snp') as f:
    all_hapmap = set()
    for line in f:
        all_hapmap.add(line.strip())


#start log file
import datetime
TimeNow=str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
Outlog.write("Allele Frequency Filtering log: "+TimeNow+"\n")
Outlog.write("Input VCF: "+os.path.abspath(options.VCFfile)+"\n")
Outlog.write("filtered by rsid in hapmap\n")

OrigCount=0
FiltCount=0
ChrCount=0
print "Start filtering..."
for line in VCF:
    OrigCount=OrigCount+1
    # Output vcf header to new vcf
    if '#' in line:
        Outvcf.write(line)
    # Start filtering variants
    if '#' not in line:
        linelist=line.split("\t")
        ChrPresent=linelist[0]
        if ChrCount != ChrPresent:
            print "Chromosome "+str(ChrPresent)
            ChrCount=ChrPresent
        if linelist[2] not in all_hapmap:
            continue
        else:
            all_hapmap.remove(linelist[2])
        # filter variants with rsid
        Outvcf.write(line)
        FiltCount=FiltCount+1

Outlog.write("  Number of variants in original VCF: "+str(OrigCount)+"\n")
Outlog.write("  Number of filtered variants: "+str(FiltCount)+"\n")
TimeNow=str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
Outlog.write("  Filtering Finished: "+TimeNow+"\n")
print "Done"

