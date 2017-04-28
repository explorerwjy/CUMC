#!/home/local/users/jw/anaconda2/bin/python
# Author: jywang	explorerwjy@gmail.com

#=========================================================================
# Filters on VCF file
# Load an yml file and perform certein filters.
#=========================================================================

import time
import argparse
import yaml
import csv
import re
import os
import pprint

PGLINE="##ExmFilter.py"
VQSR = re.compile('([\d.]+)')
pp = pprint.PrettyPrinter(indent=4)

class ExmFilter:
    def __init__(self, VCFName, YMLFilters, Ped=None, OutName=None, Debug=False):
        self.VCFName = VCFName
        if OutName == None:
            self.OutName = re.sub('\.gz$', '', self.VCFName)
            self.OutName = re.sub('\.vcf$', '', self.OutName)
            self.OutName = self.OutName + '.RareCoding.vcf'
        self.Filter = YML_Filter(YMLFilters)
        self.Pedigree = PEDIGREE(Ped)
        self.Debug = Debug
    def run(self):
        fin = GetVCF(self.VCFName)
        fout = open(self.OutName, 'wb')
        if self.Debug:
            ferr = open('FilterdOut.'+self.OutName, 'wb')
        for l in fin:
            if l.startswith('##'):
                fout.write(l)
            elif l.startswith('#'):
                fout.write(PGLINE)
                fout.write(l)
                self.Header = l.strip().split('\t')
            else:
                var = Variant(l)
                if var.FilterByCriteria(self.Filter, self.Pedigree):
                    fout.write(l)
                elif self.Debug:
                    ferr.write(var.out())
        fin.close()
        fout.close()
        if self.Debug:
            ferr.close()

class YML_Filter():
    def __init__(self, YML_FILTER_FILE):
        global PGLINE
        PGLINE = PGLINE + " YMLFilters: " +YML_FILTER_FILE + "\n"
        with open(YML_FILTER_FILE, 'rb') as ymlfile:
            self.cfg = yaml.load(ymlfile)
        self.INFO = self.cfg['INFO']
        self.READS = self.cfg['READS']
        self.Filter = self.cfg['FILTER']
        self.SNP = self.cfg['SNP']
        self.INDEL = self.cfg['INDEL']

    def show(self):
        pp.pprint(self.INFO)
        pp.pprint(self.READS)
        pp.pprint(self.FILTER)
        pp.pprint(self.SNP)
        pp.pprint(self.INDEL)

class Sample():
    # GT:AD:DP:GQ:PL    0/0:7,0:7:18:0,18,270
    def __init__(self, Format, Genotype):
        self.Format = Format.split(':')
        self.Genotype = Genotype.split(':')
        if self.Genotype[0] == './.':
            self.isCalled = False
            return
        self.isCalled = True
        self.GT = map(int, re.findall('[\d.]', self.Genotype[0]))
        self.Dict = {}
        for i, item in enumerate(self.Format):
            self.Dict[item] = self.Genotype[i]

class Variant():
    def __init__(self, record, headers):
        record = record.strip().split('\t')
        self.Chrom, self.Pos, self.Id, self.Ref, self.Alt, self.Qual, self.Filter, self.Info_str, self.Format = record[0:9]
        self.Genotypes = record[9:]
        self.Alts = self.Alt.split(',')
        self.Alleles = [self.Ref] + self.Alts
        self.GetInfo()

    def GetInfo(self):
        self.Info = {}
        tmp = self.Info_str.split(';')
        for kv in tmp:
            try:
                k, v = kv.split('=')
                if k not in self.Info:
                    self.Info[k] = v.split(',')
                else:
                    self.Info[k].extend(v.split(','))
            except:
                pass

    def CheckGT(self, idxAllele, Samples, Filters):
        AVGAD = [0, 0]
        AD1, AD2 = 0,0
        AVGPL = 0
        Total = 0
        for sample in Samples:
            if sample.info[0] == './.':
                continue
            if sample.AD[0] != 0:
                AVGAD[0] += float(sample.AD[0])
                AD1 += 1
            if sample.AD[1] != 0:
                AVGAD[1] += float(sample.AD[1])
                AD2 += 1
            AVGPL += float(sample.fmt['GQ'])
            Total += 1
        try:
            AVGAD[0] = AVGAD[0] / AD1 
        except:
            AVGAD[0] = 10
        try:
            AVGAD[1] = AVGAD[1] / AD2 
        except:
            AVGAD[1] = 10
        AVGPL = AVGPL / Total
        if Filters.READS['min_proband_AD'] != None:
            if (AVGAD[0] < Filters.READS['min_proband_AD'] and AVGAD[1] < Filters.READS['min_proband_AD']):
                return False
        if Filters.READS['min_proband_PL'] != None:
            if AVGPL < Filters.READS['min_proband_PL']:
                return False
        return True

    def CheckInfo(self, idx, Filters):
        # =============================================================================
        # Filter on Exon
        if Filters.INFO.get(['exon'], None)  == True:
            if Filters.INFO['exon_flag'] != None:
                VarFunc = self.Info.get('Func.refGene',[None]*len(self.Alts))[idx]
                if VarFunc not in Filters.INFO['exon_flag']:
                    return 'exon_flag'
        # =============================================================================
        # =============================================================================
        # Filter on AC
        if Filters.INFO.get('max_AC', None) != None:
            if int(self.Info['AC'][idx]) > int(Filters.INFO['max_AC']):
                return 'max_AC'
        # =============================================================================
        # =============================================================================
        # Filter on Population MAF        
        if Filters.INFO['max_ExAC'] != None:
            if self.GetAF(self.Info['ExAC_ALL'][idx]) > Filters.INFO['max_ExAC']:
                return 'max_ExAC'
        if Filters.INFO['max_gnomAD'] != None:
            if GetAF(self.Info['gnomAD_genome_ALL'][idx]) > Filters.INFO['max_gnomAD']:
                return 'max_gnomAD'
        if Filters.INFO['max_1KG'] != None:
            if GetAF(self.Info['1000g2015aug_all'][idx]) > Filters.INFO['max_1KG']:
                return 'max_1KG'
        # =============================================================================
        # =============================================================================
        # Filter on Exclude Gene/Chrom 
        if Filters.INFO.get('excluded_gene', None) != None:
            for gene in Filters.INFO['excluded_gene']:
                if gene in self.Info['Gene.refGene'][idx]:
                    return 'excluded_gene'
        if Filters.INFO.get('excluded_chrom', None) != None:
            if self.Chrom in Filters.INFO['excluded_chrom']:
                return 'excluded_chrom'
        # =============================================================================
        # Filter on GenomeRegion/Mappability        
        if Filters.INFO['max_seqdup'] != None:
            segdupScore = self.Info.get('genomicSuperDups',[0])[0]
            if segdupScore != '.':
                segdupScore = re.search(
                        'Score:(\d+?\.?\d+)', segdupScore).group(1)
                if float(segdupScore) >= Filters.INFO['max_seqdup']:
                    return 'max_seqdup'
        if Filters.INFO.get('Mappability', None) != None:
            if float(self.Info.get('Mappability', [1])[0]) < Filters.INFO['Mappability']:
                return 'Mappability'
        return True

    def isSNP(self, idxAllele):
        if len(self.Ref) == 1 and len(self.Alleles[idxAllele + 1]) == 1:
            return True
        else:
            return False

    def isINDEL(self, idxAllele):
        return not self.isSNP(idxAllele)

    def CheckSNP(self, Filters):
        if Filters.SNP.get('max_FS', None) != None:
            if float(self.Info.get('FS', [0])[0]) > Filters.SNP['max_FS']:
                return 'max_FS'
        if Filters.SNP.get('min_QD',None) != None:
            if float(self.Info.get('QD', [100])[0]) < Filters.SNP['min_QD']:
                return 'min_QD'
        if Filters.SNP.get('min_ReadPosRankSum', None) != None:
            if float(self.Info.get('ReadPosRankSum',[100])[0]) < Filters.SNP['min_ReadPosRankSum']:
                return 'min_ReadPosRankSum'
        return True

    def CheckINDEL(self, Filters):
        if Filters.SNP.get('max_FS', None) != None:
            if float(self.Info.get('FS', [0])[0]) > Filters.SNP['max_FS']:
                return 'max_FS'
        if Filters.SNP.get('min_QD',None) != None:
            if float(self.Info.get('QD', [100])[0]) < Filters.SNP['min_QD']:
                return 'min_QD'
        if Filters.SNP.get('min_ReadPosRankSum', None) != None:
            if float(self.Info.get('ReadPosRankSum',[100])[0]) < Filters.SNP['min_ReadPosRankSum']:
                return 'min_ReadPosRankSum'
        return True

    def CheckFilter(self, Filters):
        if self.Filter == 'PASS' or self.Filter == '.':
            return True
        v1, v2 = VQSR.findall(self.Filter)
        if float(v2) > Filters.Filter['VQSRSNP']:
            return 'VQSR'
        elif float(v2) > Filters.Filter['VQSRINDEL']:
            return 'VQSR'
        else:
            return True

    def FindAllele(self):
        Samples = []
        AllelePool = []
        for item in self.List[9:]:
            indi = Sample(self.Format, item)
            Samples.append(indi)
            AllelePool.extend(indi.GT)
        AllelePool = set(AllelePool)
        if len(AllelePool) != 2:
            return None, Samples
        return sorted(list(AllelePool))[1], Samples

    def FilterByCriteria(self, Filters):
        FLAG_INFO = self.CheckInfo(idxAllele, Filters)
        if FLAG_INFO != True:
            return FLAG_INFO
        FLAG_GT = self.CheckGT(idxAllele, Samples, Filters)
        if FLAG_GT != True:
            return FLAG_GT
        FLAG_FILTER = self.CheckFilter(Filters)
        if FLAG_FILTER != True:
            return FLAG_FILTER
        FLAG_ = self.CheckFilter(Filters)
        
        if self.isINDEL(idxAllele):
            FLAG_INDEL = self.CheckINDEL(idxAllele, Filters)
            if FLAG_INDEL != True:
                return FLAG_INDEL
        else:
            FLAG_SNP = self.CheckSNP(idxAllele, Filters)
            if FLAG_SNP != True:
                return FLAG_SNP            
        return True

    def GetAF(self, Freq):
        if Freq == '.':
            return 0
        else:
            return float(Freq)

def GetVCF(VCFin):
    if VCFin.endswith('.gz'):
        return gzip.open(VCFin, 'rb')
    else:
        return open(VCFin, 'rb')

def GetOptions():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', require=True, type=str, help='VCF record')
    parser.add_argument('-o', '--out', type=str, help='Output Prifix')
    parser.add_argument('-p', '--ped', type=str, help='Pedigree file of the VCF samples')
    parser.add_argument('-f', '--filter', default='/home/local/users/jw/CUMC/Exome-Filters-Jiayao/ALL_FILTER.yml' ,type=str, help='YML file for filter criteria')
    parser.add_argument('--debug', default=False, type=bool, help='Whether Output Filtered Variants for Debug Purpose')
    args = parser.parse_args()
    return args.vcf, args.out, args.filter, args.ped, args.debug

def main():
    VCFin, VCFout, Filters, PedFil, Debug = GetOptions()
    YMLFilters = YML_Filter(Filters)
    Ped = Pedigree(PedFil)
    instance = ExmFilter(VCFin, YMLFilters, OutName=None, Ped=Ped, Debug=Debug)
    instance.run()

if __name__ == '__main__':
    main()
