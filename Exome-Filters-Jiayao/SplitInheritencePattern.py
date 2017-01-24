#!/home/local/users/jw/anaconda2/bin/python
#Author: jywang	explorerwjy@gmail.com

#========================================================================================================
# Split vcf file by inheritence pattern given by genmod
#========================================================================================================

from optparse import OptionParser
import utils
import re

def GetOptions():
	parser = OptionParser()
	parser.add_option('-d','--dir',dest = 'dir', metavar = 'dir', help = 'InputDir contains VCF file and pedigree files')
	(options,args) = parser.parse_args()
	
	return options.dir
def MatchVcfPed(vcfs,peds):
	res = []
	for vcf in vcfs:
		for ped in peds:
			if utils.GetBaseName(vcf) == utils.GetBaseName(ped): 
				print utils.GetBaseName(vcf), utils.GetBaseName(ped)
				res.append((vcf,ped))
				break
	print "Finded vcf-ped pairs:"
	for v,p in res:
		print v,p
	print '\n'
	return res

def getGT(GTstring):
	return re.findall('[\d.]',GTstring.split(':')[0])

def CompoundTrios(CH_buffer,fout_CH,Proband_idx,Father_idx,Mother_idx):
	Flag1 = False
	Flag2 = False
	res = []
	for l in CH_buffer:
		llist = l.strip().split('\t')
		GT_Proband = getGT(llist[Proband_idx])
		if GT_Proband[0] == GT_Proband[1]:
			continue
		GT_Father = getGT(llist[Father_idx])
		GT_Mother = getGT(llist[Mother_idx])
		if GT_Proband[1] in GT_Father and GT_Proband[1] not in GT_Mother:
			Flag1 = True
			res.append(l)
		if GT_Proband[1] in GT_Mother and GT_Proband[1] not in GT_Father:
			Flag2 = True
			res.append(l)
	if Flag1 and Flag2:
		for l in res:
			fout_CH.write(l)
def CompoundDual(CH_buffer,fout_CH,Proband_idx,Parent_idx):
	Flag1 = False
	Flag2 = False
	res = []
	for l in CH_buffer:
		llist = l.strip().split('\t')
		GT_Proband = getGT(llist[Proband_idx])
		if GT_Proband[0] == GT_Proband[1]:
			continue
		GT_Parent = getGT(llist[Parent_idx])
		if GT_Proband[1] in GT_Parent:
			Flag1 = True
			res.append(l)
		if GT_Proband[1] not in GT_Parent:
			Flag2 = True
			res.append(l)
	if Flag1 and Flag2:
		for l in res:
			fout_CH.write(l)
def CompoundQuart(CH_buffer,fout_CH,Proand1_idx,Proband2_idx,i):
	return

def LookCompoundHeter(CH_buffer,fout_CH,RecordLen,Format_idx):
	if RecordLen < 2:
		return
	if RecordLen - Format_idx - 1 == 1: # Singleton
		return
	if RecordLen - Format_idx - 1 == 2: # Dual
		Proband_idx, Parent_idx = [Format_idx+1, Format_idx+2]
		CompoundDual(CH_buffer,fout_CH,Proband_idx,Parent_idx)
	if RecordLen - Format_idx - 1 == 3: # Trios
		Proband_idx,Father_idx,Mother_idx = [Format_idx+i+1 for i in range(3)]
		CompoundTrios(CH_buffer,fout_CH,Proband_idx,Father_idx,Mother_idx)
	if RecordLen - Format_idx - 1 == 4: # quart	
		CompoundQuart()
		
	return


def SplitInheritencePattern(InpDir):
	vcfs = utils.get_files(InpDir,'.tsv')
	peds = utils.get_files(InpDir,'.ped')
	vcf_peds = MatchVcfPed(vcfs,peds)
	for vcf,ped in vcf_peds:
		print "Processing vcf: %s\tpedigree: %s" % (vcf,ped)
		fout_AR = open(utils.GetBaseName(vcf)+'_AR.tsv','wb')
		fout_AR.write('Autsomal Recessive Variants\n')
		fout_AD = open(utils.GetBaseName(vcf)+'_AD.tsv','wb')
		fout_AD.write('Autsomal Dominant Variants\n')
		fout_XL = open(utils.GetBaseName(vcf)+'_XL.tsv','wb')
		fout_XL.write('X-linked Variants\n')
		fout_CH = open(utils.GetBaseName(vcf)+'_CH.tsv','wb')
		fout_CH.write('Compound Heterozygote\n')

		fin = open(vcf,'rb')
		header = fin.readline()
		fout_AR.write(header)
		fout_AD.write(header)
		fout_XL.write(header)
		fout_CH.write(header)
		headerList = header.strip().split('\t')
		RecordLen = len(headerList)
		Gene_idx = headerList.index('Annotation')
		Model_idx = headerList.index('GeneticModels')
		Format_idx = headerList.index('FORMAT')	
		CH_buffer = []
		LastGene = None
		for l in fin:
			llist = l.strip().split('\t')
			Gene = llist[Gene_idx]
			Model = llist[Model_idx].split(':')[-1]
			if 'AR' in Model:
				fout_AR.write(l)
			if 'AD' in Model:
				fout_AD.write(l)
			if 'X' in Model:
				fout_XL.write(l)

			if Gene != LastGene:
				if LastGene != None:
					LookCompoundHeter(CH_buffer,fout_CH,RecordLen,Format_idx)
				CH_buffer = [l]
				LastGene = Gene
			else:
				CH_buffer.append(l)
		LookCompoundHeter(CH_buffer,fout_CH,RecordLen,Format_idx)
		fout_AR.close()
		fout_AD.close()
		fout_XL.close()
		fout_CH.close()

def main():
	InpDir = GetOptions()
	SplitInheritencePattern(InpDir)
	return

if __name__=='__main__':
	main()
