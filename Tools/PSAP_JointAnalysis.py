#!/home/local/users/jw/anaconda2/bin/python
#Author: jywang	explorerwjy@gmail.com

#========================================================================================================
# PSAP_JointAnalysis.py
#========================================================================================================

import argparse
import csv

cutoff = 0.05

class Ped:
	def __init__(self, fname):
		self.fname = fname
		self.hand = open(fname, 'rb')
		self.header = self.hand.readline().strip().split('\t')
	def GetProband(self):
		self.hand.seek(1)
		for l in self.hand:
			llist = l.strip().split('\t')
			fam, _id, phenotpye, phenotpyeDetail = llist[0], llist[1], llist[5], llist[7]
			if fam == _id and phenotpye == '2':
				father, mother = llist[2], llist[3]
				return _id, father, mother, phenotpyeDetail
		print "Error when parsing pedigree, No proband find."

Header = ['Proband', 'phenotpyeDetail', 'Chrom', 'Pos', 'Ref', 'Alt', 'CADD_Phred', 'Dz.Model.Proband', 'popScore.Proband', 'vid','validation','Flag', 'Gene','GeneName', 'Allele Count','GeneFunc','ExonicFunc','AAchange',	'ExAC_ALL',	'gnomAD_genome_ALL', '1KG',	'VarType', 'MetaSVM', 'CADD13', 'PP2', 'MCAP', 'mis_z', 'lof_z', 'pLI', 'pRec', 'HeartRank', 'LungRank', 'BrainRank', 'Filter', 'Proband', 'Father', 'Mother']
HeaderToKeep = ['Chrom', 'Pos', 'Ref', 'Alt', 'CADD_Phred', 'vid','validation','Flag', 'Gene','GeneName', 'Allele Count','GeneFunc','ExonicFunc','AAchange', 'ExAC_ALL', 'gnomAD_genome_ALL', '1KG', 'VarType', 'MetaSVM', 'CADD13', 'PP2', 'MCAP', 'mis_z', 'lof_z', 'pLI', 'pRec', 'HeartRank', 'LungRank', 'BrainRank', 'Filter']
class OnePSAPReport:
	def __init__(self, fname, ped):
		self.fname = fname
		self.hand = open(fname, 'rb')
		self.reader = csv.reader(self.hand)
		self.header = self.reader.next()
		self.ped = ped
	def FilterHeader(self):
		self.proband, self.father, self.mother, self.phenotpyeDetail = self.ped.GetProband()
		print self.proband, self.father, self.mother, self.phenotpyeDetail
		self.idxproband = self.header.index(self.proband+'.GT')
		self.header2keep = HeaderToKeep + ['Dz.Model.'+proband, 'popScore.'proband, proband+'.GT', father+'.GT', mother+'.GT']
		self.idx2keep = []
		for i in range(len(self.header)):
			if self.header[i] in self.header2keep:
				self.idx2keep.append(i)
		self.idx.extend(range(len(self.header)-3:len(self.header)))
		self.idx_VarType = self.header.index('VarType')
		self.idx_Pvalue = self.header.index('popScore.'proband)
	def GetRecords(self):
		for row in self.reader:
			if self.Pass(row):
				yield self.trimmRecord(row)
	def Pass(self, row):
		idx_proband_allele = self.GetProbandAlleleIdx(row)
		# Not Include slient variants
		VarType = row(self.idx_VarType).split(',')
		if VarType[idx_proband_allele] == 'slient':
			return False
		if self.idx_Pvalue < cutoff:
			return False
		return True
	def GetProbandAlleleIdx(self):
		GT = row[self.idxproband].split(':')[0]
		GT = GT.split('/')
		if len(GT) == 1:
			GT = GT.split('|')
		return int(GT[1]) - 1


def GetOptions():
	parser = argparse.ArgumentParser()
	parser.add_argument('-d','--dir', type=str, help = 'dir contains psap result from PSAP_Report.py')
	parser.add_argument('-o','--out', type=str, help = 'OutPut Name')
	parser.add_argument('-p','--pvalue', type=str, help = 'PSAP Pvalue to cutoff')
	args = parser.parse_args()
	if args.pvalue == None:
		args.pvalue = 0.05
	if args.out == None:
		args.out = args.pvalue + '.csv'
	return args.dir, args.out, args.pvalue

def main():
	Dir, Out, Pvalue = GetOptions()
	
	return

if __name__=='__main__':
	main()
