#!/home/yufengshen/anaconda2/bin/python
# Given set of vcf files, get the avg of N denovo and N inherited variants
import argparse
import re
import os
import utils
import time
import yaml
import csv
import pprint

CSV_HEADER = ['Sample', 'Phenotype', 'Chrom', 'Pos', 'Ref', 'Alt', 'AC', 'Gene', 'GeneFunc', 'ExonicFunc', 'AAchange', 'ExAC', 'gnomAD', 'MCAP', 'MetaSVM', 'CADD', 'PP2', 'ESP', '1KG', 'mis_z', 'lof_z','pLI', 'pRec','HeartRank', 'LungRank', 'BrainRank','Filter', 'QUAL', 'ProbandGT', 'FatherGT', 'MotherGT', 'Relateness']
VQSR = re.compile('([\d.]+)')
pp = pprint.PrettyPrinter(indent=4)

EXAC_GENE_SCORE = '/home/yufengshen/resources/ExAC/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt'
DIAPHRAGM_EXP = '/home/yufengshen/resources/GeneRanks/diaphragm_rank.csv'
LUNG_EXP = '/home/yufengshen/resources/GeneRanks/Lung_rank_RNAseq_Asselin-Labat-GPL13112_human.csv'
MOUSEBRAIN_EXP = '/home/yufengshen/resources/GeneRanks/mousebrain.csv'

class GeneScore():
    def __init__(self, ExAC=True, Diaphragm=False, Lung=False, MouseBrain=False):
        self.ExAC, self.Diaphram, self.Lung, self.MouseBrain = None, None, None, None
        if ExAC == True:
            self.Load_ExAC()
        if Diaphragm == True:
            self.Load_Diaphragm()
        if Lung == True:
            self.Load_Lung()
        if MouseBrain == True:
            self.Load_MouseBrain()
    def Load_ExAC(self):
        stime = time.time()
        fin = open(EXAC_GENE_SCORE,'rb')
        header = fin.readline().strip().split('\t')
        self.ExAC = {}
        for l in fin:
            llist = l.strip().split('\t')
            tmp = {}
            for k,v in zip(header, llist):
                tmp[k] = v
            self.ExAC[tmp['gene']] = tmp
        print 'Finished Reading Gene Score ExAC %.3f'%-(stime - time.time())
    def Load_Diaphragm(self):
        stime = time.time()
        Heart = open(DIAPHRAGM_EXP,'rb')
        Reader = csv.reader(Heart)
        Heart_head = Reader.next()
        heart = Heart_head.index('rank')
        self.Diaphragm = {}
        for row in Reader:
            self.Diaphragm[row[0]]=str(100 - float(row[heart]))
        print 'Finished Reading Gene Score Diaphragm %.3f'%(time.time() - stime) 
    def Load_Lung(self):
        stime = time.time()
        fin = open(LUNG_EXP,'rb')
        Reader = csv.reader(fin)
        head = Reader.next()
        lung = head.index('Control-Stroma rank')
        self.Lung = {}
        for row in Reader:
            self.Lung[row[1]]= row[lung]
        print 'Finished Reading Gene Score Lung %.3f'%(time.time() - stime) 
    def Load_MouseBrain(self):
        stime = time.time()
        Brain = open(MOUSEBRAIN_EXP,'rb')
        Reader = csv.reader(Brain)
        Header = Reader.next()
        Brain = Header.index('brain_rank')
        self.MouseBrain = {}
        for row in Reader:
            self.MouseBrain[row[0]] = row[Brain] 
        print 'Finished Reading Gene Score MouseBrain %.3f'%(time.time() - stime) 

class YML_Filter():
    def __init__(self, yml_dict):
        print yml_dict
        self.INFO = yml_dict['INFO']
        self.READS = yml_dict['READS']
        self.FILTER = yml_dict['FILTER']
        self.SNP = yml_dict['SNP']
        self.INDEL = yml_dict['INDEL']
    def show(self):
        pp.pprint(self.INFO)
        pp.pprint(self.READS)
        pp.pprint(self.FILTER)
        pp.pprint(self.SNP)
        pp.pprint(self.INDEL)

class Sample():
    # GT:AD:DP:GQ:PL	0/0:7,0:7:18:0,18,270
    def __init__(self, name, Format, tmp):
        self.name = name
        self.Format = tmp
        self.formats = Format.split(':')
        self.info = tmp.split(':')
        self.GT = map(int, re.findall('[\d.]', self.info[0]))
        self.fmt = {}
        for i, item in enumerate(self.formats):
            self.fmt[item] = self.info[i]
        if 'AD' in self.fmt:
            tmp = self.fmt['AD'].split(',')
            self.AD = [0]*2
            self.AD[0] = float(tmp[self.GT[0]])
            self.AD[1] = float(tmp[self.GT[1]])
    def show(self):
        return '{}: {}:{}:{}\t'.format(self.name,self.fmt['GT'],self.fmt['AD'],self.fmt['GQ'])
   
class Variant():
    def __init__(self, record, headers):
        record = record.strip().split('\t')
        self.List = record
        CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT = record[0:9]
        self.Chrom = CHROM
        self.Pos = POS
        self.Id = ID
        self.Ref = REF
        self.Alt = ALT
        self.Alts = ALT.split(',')
        self.Alleles = [self.Ref] + self.Alts
        self.Qual = QUAL
        self.Filter = FILTER
        self.Info_str = INFO
        self.GetInfo(INFO)
        self.Format = FORMAT

    def GetInfo(self, INFO):
        res = {}
        tmp = INFO.split(';')
        for kv in tmp:
            try:
                k,v = kv.split('=')
                if k not in res:
                    res[k] = v.split(',')
                else:
                    res[k].extend(v.split(','))
            except:
                continue
        self.Info = res
    
    def CheckGT(self, Proband, Father, Mother, Filters):
        #print Filters.READS
        #exit()
        Pool = set([Father.GT[0], Father.GT[1], Mother.GT[0], Mother.GT[1]])
        if Proband.GT[0] not in Pool or Proband.GT[1] not in Pool:
            #print Proband.GT, Father.GT, Mother.GT
            if Filters.READS['min_proband_AD'] != None: 
                if (Proband.AD[0] < Filters.READS['min_proband_AD'] or Proband.AD[1] < Filters.READS['min_proband_AD']) :
                    return False
                #print Proband.AD
            if Filters.READS['min_proband_PL'] != None :
                if (Proband.fmt['GQ']) < Filters.READS['min_proband_PL']:
                    return False
            if Filters.READS['min_proband_alt_freq'] != None :
                if  (Proband.AD[0]/float(Proband.fmt['DP']) < Filters.READS['min_proband_alt_freq'] or Proband.AD[1]/float(Proband.fmt['DP']) < Filters.READS['min_proband_alt_freq']):
                    return False
            if Filters.READS['min_parents_DP'] != None :
                try:
                    if  (int(Father.fmt['DP']) < Filters.READS['min_parents_DP'] or int(Mother.fmt['DP'])< Filters.READS['min_parents_DP']):
                        return False
                except ValueError:
                    return False
            if Filters.READS['min_parents_GQ'] != None :
                if  (Proband.fmt['GQ']) < Filters.READS['min_parents_GQ']:
                    return False
            if Filters.READS['min_parents_ref_freq'] != None :
                if (float(Father.AD[0])/float(Father.fmt['DP']) < Filters.READS['min_parents_ref_freq'] and Mother.AD[0] / float(Mother.fmt['DP']) < Filters.READS['min_parents_ref_freq'] ):
                    return False
            #print Proband.show(), Father.show(), Mother.show()
            return True
        else:
            return False
    def CheckInfo(self, Proband, Father, Mother, Filters):
        idx = Proband.GT[1]-1
        if Filters.INFO['exon'] != None and Filters.INFO['exon'] == True:
            if Filters.INFO['exon_flag'] != None:
                VarFunc = self.Info['Func.refGene'][idx]
                #VarFunc = self.Info.get('Func.refGene',[None]*len(self.Alts))[idx]
                #print VarFunc
                if VarFunc not in Filters.INFO['exon_flag']:
                    return False
            #print VarFunc, Proband.show(), Father.show(), Mother.show()
        if Filters.INFO['max_AC'] != None:
            if  int(self.Info['AC'][idx]) > Filters.INFO['max_AC']:
                return False
        #print VarFunc, self.Info['AC'][idx], Proband.show(), Father.show(), Mother.show()
        if Filters.INFO['max_ExAC'] != None:
            if AF(self.Info['ExAC_ALL'][idx]) > Filters.INFO['max_ExAC']:
                return False
        if Filters.INFO['max_gnomAD'] != None:
            if AF(self.Info['gnomAD_genome_ALL'][idx]) > Filters.INFO['max_gnomAD']:
                return False
        if Filters.INFO['max_1KG'] != None:
            if AF(self.Info['1000g2015aug_all'][idx]) > Filters.INFO['max_1KG']:
                return False
       # print AF(self.Info['ExAC_ALL'][idx]), VarFunc, self.Info['AC'][idx], Proband.show(), Father.show(), Mother.show()
        if Filters.INFO['excluded_gene'] != None:
            for gene in Filters.INFO['excluded_gene']:
                if gene in self.Info['Gene.refGene'][idx]:
                    return False
        if Filters.INFO['excluded_chrom'] != None:
            if self.Chrom in Filters.INFO['excluded_chrom']:
                return False
        #print self.Info['Gene.refGene'][idx], AF(self.Info['ExAC_ALL'][idx]), VarFunc, self.Info['AC'][idx], Proband.show(), Father.show(), Mother.show()
        if Filters.INFO['max_seqdup'] != None:
            segdupScore = self.Info['genomicSuperDups'][idx]
            if segdupScore != '.':
                #print self.Info['genomicSuperDups'][idx]
                segdupScore = re.search('Score:(\d+?\.?\d+)',segdupScore).group(1)
                #print segdupScore
                if float(segdupScore) >= Filters.INFO['max_seqdup']:
                    return False
        print self.Info['Gene.refGene'][idx], AF(self.Info['ExAC_ALL'][idx]), VarFunc, self.Info['AC'][idx], Proband.show(), Father.show(), Mother.show(), segdupScore
        return True

    def isSNP(self, Proband):
        if len(self.Alleles[Proband.GT[0]]) == 1 and len(self.Alleles[Proband.GT[1]]) == 1:
            return True
        else:
            return False
    
    def CheckSNP(self, Proband, Filters):
        idx = Proband.GT[1] - 1
        if Filters.SNP['max_FS'] != None:
            #print self.Info['FS']
            if float(self.Info['FS'][0]) > Filters.SNP['max_FS']:
                return False
        if Filters.SNP['min_QD'] != None:
            if float(self.Info['QD'][0]) < Filters.SNP['min_QD']:
                return False
        if Filters.SNP.get('min_ReadPosRankSum',None) != None:
            if float(self.Info['ReadPosRankSum']) > Filters.SNP['min_ReadPosRankSum']:
                return False
        return True

    def CheckINDEL(self, Proband, Filters):
        idx = Proband.GT[1] - 1
        #print self.Alleles, self.Info['FS'], self.Info['QD'], self.Info['ReadPosRankSum']
        if Filters.INDEL['max_FS'] != None:
            if float(self.Info['FS'][0]) > Filters.INDEL['max_FS']:
                return False
        if Filters.INDEL['min_QD'] != None:
            if float(self.Info['QD'][0]) < Filters.INDEL['min_QD']:
                return False
        if Filters.INDEL['min_ReadPosRankSum'] != None:
            try:
                if float(self.Info['ReadPosRankSum'][0]) > Filters.INDEL['min_ReadPosRankSum']:
                    return False
            except KeyError:
                pass
        return True

    def isINDEL(self, Proband):
        return not self.isSNP(Proband)

    def CheckFilter(self, Proband, Filter):
        if self.Filter == 'PASS' or Filter == '.':
            return True
        v1, v2 = VQSR.findall(Filter)
        if float(v2) > Filters.Filter['VQSRSNP']:
            return False
        else:
            return True

    def CheckDeNovo(self, headers, Pedigree, Filters):
        Proband = self.List[headers.index(Pedigree.Proband.Sample)]
        Father = self.List[headers.index(Pedigree.Father.Sample)]
        Mother = self.List[headers.index(Pedigree.Mother.Sample)]
        try:    
            Proband = Sample(Pedigree.Proband.Sample, self.Format, Proband)
            Father = Sample(Pedigree.Father.Sample, self.Format, Father)
            Mother = Sample(Pedigree.Mother.Sample, self.Format, Mother)
        except ValueError:
            #print self.List[headers.index(Pedigree.Proband.Sample)], self.List[headers.index(Pedigree.Father.Sample)], self.List[headers.index(Pedigree.Mother.Sample)] 
            return False, None, None, None 
        if not self.CheckGT(Proband, Father, Mother, Filters):
            return False, None, None, None
        else:
            #print 'pass GT'
            pass 
        if not self.CheckInfo(Proband, Father, Mother, Filters):
            return False, None, None, None
        else:
            #print 'pass info'
            pass
        if not ( (self.isINDEL(Proband) and self.CheckINDEL(Proband, Filters)) or (self.isSNP(Proband) and self.CheckSNP(Proband, Filters)) ):
            return False, None, None, None
        else:
            print 'pass snp/indel'
        return True, Proband, Father, Mother

    def show(self):
        print '\t'.join(self.List)

    def OutAsCSV(self, Proband, Father, Mother, Pedigree, genescore):
        Gene = self.Info['Gene.refGene'][Proband.GT[1]-1]
        try:    
            mis_z = genescore.ExAC[Gene]['mis_z']
            lof_z = genescore.ExAC[Gene]['lof_z']
            pLI = genescore.ExAC[Gene]['pLI']
            pRec = genescore.ExAC[Gene]['pRec']
        except:
            print 'None ExAC find'
            mis_z ,lof_z ,pLI ,pRec = ['NA'] * 4
        try:
            LungRank = genescore.Lung[Gene]
        except:
            print 'None LungRank find'
            LungRank = 'NA'
        try:
            HeartRank = genescore.Diaphragm[Gene]
        except:
            print 'None HeartRank find'
            HeartRank = 'NA'
        try:
            BrainRank = genescore.MouseBrain[Gene]
        except:
            print 'None BrainRank find'
            BrainRank = 'NA'
        #print Proband
        Indi = Pedigree.GetIndi(Proband.name)
        Phenotype = Indi.PhenotypeDetail
        Relateness = Indi.Relateness
        return [Proband.name, Phenotype, self.Chrom, self.Pos, self.Ref, self.Alt, ','.join(str(x) for x in self.Info['AC']),','.join( self.Info['Gene.refGene']),','.join( self.Info['Func.refGene']), ','.join(self.Info['ExonicFunc.refGene']), ','.join(self.Info['AAChange.refGene']),','.join(str(AF(x)) for x in self.Info['ExAC_ALL']), ','.join(str(AF(x)) for x in self.Info['gnomAD_genome_ALL']), ','.join(self.Info['MCAP']), ','.join(self.Info['MetaSVM_pred']),','.join(self.Info['CADD_phred']),','.join(self.Info['Polyphen2_HDIV_pred']) , ','.join(str(AF(x)) for x in self.Info['esp6500siv2_all']), ','.join(str(AF(x)) for x in self.Info['1000g2015aug_all']), mis_z, lof_z, pLI, pRec, HeartRank, LungRank, BrainRank, self.Filter, self.Qual, Proband.Format, Father.Format, Mother.Format ,Relateness   ]


class Individual():
    def __init__(self, List, Header):
        self.Fam, self.Sample, self.Father, self.Mother, self.Gender, self.Pheno = List[:6]
        self.PhenotypeDetail = List[Header.index('PhenotypeDetail')]
        self.Relateness = List[Header.index('Relateness')]
        self.SexCheck = List[Header.index('SexCheck')]
    def show(self):
        print self.Fam, self.Sample, self.Father, self.Mother

class Pedigree():
    def __init__(self, PedFil):
        fin = open(PedFil, 'r')
        self.individuals = []
        for l in fin:
            if l.startswith('#'):
                Header = l.strip().split('\t') 
            indi = l.strip().split('\t')
            indi = Individual(indi, Header)
            indi.show()
            self.individuals.append(indi)
        self.Proband, self.Father, self.Mother = None, None, None
        for ind in self.individuals:
            if ind.Fam == ind.Sample:
                self.Proband = ind
        for ind in self.individuals:
            if self.Proband.Father == ind.Sample:
                self.Father = ind
            if self.Proband.Mother == ind.Sample:
                self.Mother = ind
    def show(self):
        print 'Proband:%s\tFather:%s\tMother:%s' % (self.Proband.Sample, self.Father.Sample, self.Mother.Sample)
    def GetIndi(self, Sample):
        for Indi in self.individuals:
            if Indi.Sample == Sample:
                return Indi
        return None

def GetOptions():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dir", type=str, required=True, help="<Required> Dir contains family Ped and VCF files")
    parser.add_argument("-f", "--filter", type=str, required=True, help="<Required> ymal file contains filter criteria")
    args = parser.parse_args()
    return args.dir, args.filter

def MatchVcfPed(vcfs, peds):
    res = []
    for vcf in vcfs:
        for ped in peds:
            if utils.GetBaseName(vcf) == utils.GetBaseName(ped):
                print utils.GetBaseName(vcf), utils.GetBaseName(ped)
                res.append((vcf, ped))
                break
    print "Finded vcf-ped pairs:"
    for v, p in res:
        print v, p
        print '\n'
    return res

def AF(CHR):
    if CHR == '.':
        return 0
    else:
        return float(CHR)

def For_one_vcf(VCF, Ped, Writer, Filters, CSV_HEADER, gene_score):
    Ped = Pedigree(Ped)
    fin = open(VCF, 'rb')
    for line in fin:
        if line.startswith("##"):
            # Meta info
            continue 
        elif line.startswith('#'):
            # header
            header = line.strip().split('\t')
        else:
            # variants line
            var = Variant(line, header)
            isDenovo, Proband, Father, Mother = var.CheckDeNovo(header, Ped, Filters)
            if isDenovo:
                #var.show()
                Writer.writerow(var.OutAsCSV(Proband, Father, Mother, Ped, gene_score))

def Call_DeNovo(InpDir, Filters):
    vcfs = utils.get_files(InpDir, '.vcf')
    peds = utils.get_files(InpDir, '.ped')
    vcf_peds = MatchVcfPed(vcfs, peds)
    OutFil = open('Denovo.csv','wb')
    Writer = csv.writer(OutFil, delimiter=',')
    #CSV_HEADER[0] = '#'+CSV_HEADER[0]
    Writer.writerow(CSV_HEADER)
    gene_score = GeneScore(ExAC=True, Diaphragm=True, Lung=True, MouseBrain=True)
    print gene_score
    for vcf, ped in vcf_peds:
        print "Processing vcf: %s\tpedigree: %s" % (vcf, ped)
        s_time = time.time()
        #if not vcf.endswith('CARE19.vcf'):
        #    continue
        For_one_vcf(vcf, ped, Writer, Filters, CSV_HEADER, gene_score)
        print "%s Finished, used %.3fs"%(vcf, time.time()-s_time)

def Parse_YAML(yaml_fname):
    with open(yaml_fname, 'rb') as ymlfile:
        cfg = yaml.load(ymlfile)
        Filters = YML_Filter(cfg)
    Filters.show()
    return Filters

def main():
    dirc, yaml_fname = GetOptions()
    Filters = Parse_YAML(yaml_fname)
    dirc = os.path.abspath(dirc)
    if not dirc.endswith('/'):
        dirc += '/'
    print "Vcf list is", dirc
    Call_DeNovo(dirc, Filters)


if __name__ == '__main__':
    main()
