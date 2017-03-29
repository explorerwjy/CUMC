#!/usr/bin/env python

# The purpose of this script is to filter out unwanted variants from a VCF file by ExAC frequency, and variant function
#    -v/--vcf     <required>    The script requires an input vcf
#    -o/--out     <required>    The use should specify a base name for the output files.The script outputs the filtered results as a vcf file. The script also outputs a log file. 
#    -m/--maf     <optional>    Minor allele frequency for filtering. Default is 0.1


from optparse import OptionParser
import gzip
import os

parser = OptionParser()
# Basic Input Files, required
parser.add_option("-v", "--vcf", dest="VCFfile",help="input VCF file", metavar="VCFfile")
parser.add_option("-o", "--output", dest="OutputFileName",help="user specified name of file for output to be written", metavar="OutputFileName")
parser.set_defaults(MAFFilter="0.1")
parser.add_option("-m", "--maf", dest="MAFFilter",help="user specified allele frequency", metavar="MAFFilter")


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
MafCutOff=float(options.MAFFilter)

## Variant classes
MissenseClass=['nonsynonymousSNV','unknown']
NonFrameshiftClass=['nonframeshiftdeletion','nonframeshiftinsertion','nonframeshiftsubstitution']
FrameshiftClass=['frameshiftdeletion','frameshiftinsertion','frameshiftsubstitution']
NonsenseClass=['stopgain','stoploss']
InDelClass=NonFrameshiftClass+FrameshiftClass
SplicingClass=['splicing','exonic,splicing']

#open input and output files
VcfOutputFilename=BaseName+'.filter.aaf.vcf'
LogOutputFilename=BaseName+'.filter.aaf.log'
Outvcf=open(VcfOutputFilename,'w')
Outlog=open(LogOutputFilename,'w')

#start log file
import datetime
TimeNow=str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
Outlog.write("Allele Frequency Filtering log: "+TimeNow+"\n")
Outlog.write("Input VCF: "+os.path.abspath(options.VCFfile)+"\n")
Outlog.write("  ExAC allele frequency maximum: "+str(MafCutOff)+"\n")
    
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
        # Get INFO
        INFOstring=linelist[7]
        INFOcolumnList=INFOstring.split(";")
        INFOdict={}
        for element in INFOcolumnList:
            if '=' in element:
                FieldName,FieldValue=element.split('=',1)
                INFOdict[FieldName]=FieldValue
        
        #Variant details
        GeneName=INFOdict.get('GeneName','.')
        VariantFunction=INFOdict.get('VarFunc','none')
        VariantClassList=INFOdict.get('VarClass','none').split(',')
        
        #ExAC frequency:
        ExACscore=str(INFOdict.get('ExACfreq',0))
        ExACscore=ExACscore.split(",")
        ExACscore=ExACscore[0]
        if ExACscore == ".":
            ExACscore=float(0)
        else:
            ExACscore=float(ExACscore)
        
        # Check if KG passes threshold
        PassMAF=False
        if float(ExACscore) <= MafCutOff:
            PassMAF=True
        
         ## Check if Variant class passes
        PassFunction=False
        if VariantFunction=='exonic' or VariantFunction in SplicingClass or VariantFunction=='none':
            PassFunction=True
            
        PassClass=False
        if any( str(i) != 'synonymousSNV' for i in VariantClassList):
            PassClass=True
        
        if PassMAF and PassFunction:
            Outvcf.write(line)
            FiltCount=FiltCount+1

Outlog.write("  Number of variants in original VCF: "+str(OrigCount)+"\n")
Outlog.write("  Number of variants in output file: "+str(FiltCount)+"\n")
TimeNow=str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
Outlog.write("  Filtering Finished: "+TimeNow+"\n")
print "Done"

