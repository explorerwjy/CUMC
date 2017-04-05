#!/home/local/users/jw/anaconda2/bin/python
# Author: jywang	explorerwjy@gmail.com

#=========================================================================
# Annotate the variant list
# Variant list according to PSAP report.
#=========================================================================

import time
import argparse
import yaml
import csv
import os
import pprint
from MyVCF import *

#CSV_HEADER = ['Chrom', 'Pos', 'Ref', 'Alt', 'AC', 'Gene', 'GeneName', 'GeneFunc', 'ExonicFunc', 'AAchange', 'ExAC', 'gnomAD_genome', 'VariantType', 'MetaSVM', 'CADD',
		'PP2', 'MCAP', '1KG', 'mis_z', 'lof_z', 'pLI', 'pRec', 'HeartRank', 'LungRank', 'BrainRank', 'Filter', 'QUAL']

#header = ['CHROM', 'POS', 'REF', 'ALT', 'Gene', 'GeneName', 'GeneFunc', 'ExonicFunc', 'AAchange', 'gnomAD_genome', 'ExACfreq', '1KGfreq', 'VariantType', 'MetaSVM', 'CADD', 'Polyphen2', 'Gene pLI',
			'Gene mis-Z', 'Gene lof-Z', 'Gene pRec', 'HeartExpressionRank', 'LungExpressionRank', 'MouseBrianRank', Proband + '(Proband)', Father + '(Father)', Mother + '(Mother)', 'Dz.Model', 'PSAP p-value', 'Flag']

HEADERS_PSAP_TO_BE_POP = ['Gene.wgEncodeGencodeBasicV19', 'Func.wgEncodeGencodeBasicV19', 'ExonicFunc.wgEncodeGencodeBasicV19','AAChange.wgEncodeGencodeBasicV19','mac63kFreq_ALL','1000g2014sep_all','esp6500si_all']
CSV_HEADER = ['Gene', 'GeneName', 'Allele Count', 'GeneFunc', 'ExonicFunc', 'AAchange', 'ExAC_ALL', 'gnomAD_genome_ALL', '_1KG', 'VarType', 'MetaSVM', 'CADD13', 'PP2', 'MCAP', 'mis_z', 'lof_z', 'pLI', 'pRec','HeartRank', 'LungRank', 'BrainRank']

def GetOptions():
	parser = argparse.ArgumentParser()
	parser.add_argument('-r', '--report', type=str, help='PSAP Report')
	parser.add_argument('-v', '--vcf', type=str, help='VCF record')
	parser.add_argument('-p', '--ped', type=str, help='Pedigree')
	parser.add_argument('-o', '--out', type=str, help='Output Prifix')
	args = parser.parse_args()
	if args.out == None:
		args.out = args.ped.rstrip('ped')
	return args.report, args.vcf, args.ped, args.out

class PSAP_REPORT(object):
	"""docstring for PSAP_REPORT"""
	def __init__(self, PsapReportTxT, VarDict, Genes):
		#super(PSAP_REPORT, self).__init__()
		self.VarDict = VarDict
		self.Genes = Genes.Genes
		self.fin = open(PsapReportTxT, 'rb')
		self.FullHeader = self.fin.readline.strip().split('\t')
		TmpHeader = self.FullHeader
		pop_idxes = [TmpHeader.index(item) for item in HEADERS_PSAP_TO_BE_POP]
		pop_idxes.sort()
		self.pop_idxes = pop_idxes
		self.Pop_non_display(TmpHeader)
		self.Header = TmpHeader
		

	def RunReport(self, Writer):
		for l in self.fin:
			row = l.strip().split('\t')
			self.Pop_non_display(row)
			record = Record()
			record.read_PsapReport(row, self.Header)
			record.fetch_VCF(self.VarDict)
			others = record.FormRecord(self.Genes)
			res = row + others
			Writer.writerow(res)

	def Pop_non_display(self, row):
		for idx in self.pop_idxes:
			row.pop(idx)

	def read_VCF(self, vcf_dict):
		if self.ALT != '-':
			self.VCF = vcf_dict.get(self.CHROM + ':' + self.POS)
		elif self.ALT == '-':
			self.VCF = vcf_dict.get(self.CHROM + ':' + str(int(self.POS) - 1))

class Record():
	def __init__(self):
		pass

	def read_PsapReport(self, l, header, Ped):
		self.CHROM = llist[header.index('Chr')]
		self.POS = llist[header.index('Start')]
		self.REF = llist[header.index('Ref')]
		self.ALT = llist[header.index('Alt')]


	def fetch_VCF(self, vcf_dict):
		try:
			if self.ALT != '-':
				self.VCF = vcf_dict[self.CHROM + ':' + self.POS]
			elif self.ALT == '-':
				self.VCF = vcf_dict[self.CHROM + ':' + str(int(self.POS) - 1)]
		except KeyError:
			print "KeyError while fetch record in VCF, {}:{}\t{}-{} ".format(self.CHROM, self.POS, self.REF, self.ALT)
			exit()

	def FormRecord(self, Ped, genescore):
		Gene = self.Info['Gene.refGene'][0]
		GeneName = genescore[Gene].Name
		AC = ','.join(self.Info['AC'])
		GeneFunc = ','.join(self.Info['Func.refGene'])
		ExonicFunc = ','.join(self.Info['ExonicFunc.refGene'])
		AAchange = ','.join(self.Info['AAchange'])
		ExAC_ALL = ','.join(str(AF(self.Info['ExAC_ALL'])))
		gnomAD_genome_ALL = ','.join(str(AF(self.Info['gnomAD_genome_ALL'])))
		MCAP = ','.join(self.Info['MCAP'])
		MetaSVM = ','.join(self.Info['MetaSVM_pred'])
		CADD = ','.join(self.Info['CADD_phred'])
		PP2 = ','.join(self.Info['Polyphen2_HDIV _pred'])
		VarType = self.VCF.GetVarType(GeneFunc, ExonicFunc, MetaSVM, CADD, PP2)
		_1KG = ','.join(str(AF(self.Info['1000g2015aug_all'])))
		mis_z = genesocre[Gene].Mis_z
		lof_z = genesocre[Gene].Lof_z
		pLI = genesocre[Gene].pLI
		pRec = genesocre[Gene].pRec
		HeartRank = genesocre[Gene].DiaphragmRank
		LungRank = genesocre[Gene].LungRank
		BrainRank = genesocre[Gene].MouseBrainRank

		#return '\t'.join([self.CHROM, self.POS, self.REF, self.ALT, self.Gene, self.GeneName, self.GeneFunc, self.ExonicFunc, self.AAchange, self.ExACfreq, self.gnomAD_genome, self._1KGfreq, self.MetaSVM, self.CADD, self.PolyPhen2, self.pLI, self.misZ, self.lofZ, self.pREC, self.HRank, self.LRank, self.BRank, self.Proband, self.Father, self.Mother, self.DzModel, self.Pval, self.Flag])
		INFO = [Gene, GeneName, AC, GeneFunc, ExonicFunc, AAchange, ExAC_ALL, gnomAD_genome_ALL, _1KG, VarType, MetaSVM, CADD, PP2, MCAP, mis_z, lof_z, pLI, pRec, HeartRank, LungRank, BrainRank]
		GTs = self.VCF.getSampleGenotypes()
		return INFO + GTs


def GenerateReport(psap_report, vcf, ped, OUT):
	f_psap = open(psap_report, 'rb')
	f_vcf = open(vcf, 'rb')
	Ped = Pedigree(ped)
	Proband = Ped.Proband.Sample
	Father = Ped.Father.Sample
	Mother = Ped.Mother.Sample
	gene_score = GENE_ANNOTATION()
	vcf_dict = GetVCF(f_vcf, Ped, gene_score, OUT)
	f_vcf.close()

	fout = open(OUT + 'PSAP.csv', 'wb')
	fout.write('\t'.join(header) + '\n')
	psap_header = f_psap.readline().strip().split('\t')
	for l in f_psap:
		record = Record()
		record.read_PsapReport(l, psap_header)
		record.read_VCF(vcf_dict)
		# record.read_MCAP(MCAP_dict)
		record.read_Gene(Gene_dict)
		record.read_GeneRank(HRank, LRank, PPT)
		line = record.FormRecord()
		fout.write(line + '\n')
	f_psap.close()
	fout.close()


def GetVCF(f_vcf, Ped, gene_score, OUT):
	print 'Reading VCF'
	stime = time.time()
	res = {}
	INDEL_OUT = open(OUT + 'INDEL.csv', 'wb')
	Writer = csv.writer(INDEL_OUT, delimiter=',')
	Writer.writerow(CSV_HEADER)
	Filters = Parse_YAML(ALL_FILTER)
	for l in f_vcf:
		if l.startswith('##'):
			continue
		elif l.startswith('#'):
			vcf_header = l.strip().split('\t')
		else:
			llist = l.strip().split('\t')
			var = Variant(l, vcf_header)
			res[llist[0] + ':' + llist[1]] = var
			_pass, Proband, Father, Mother= var.ProbandisINDEL(Ped, Filters)
			if _pass:
				Writer.writerow(var.OutAsCSV(Proband, Father, Mother, Ped, gene_score))

	print 'Finished Reading VCF %.3f' % -(stime - time.time())
	return res



def main():
	PSAP, VCF, PED, OUT = GetOptions()
	GenerateReport(PSAP, VCF, PED, OUT)
	return


if __name__ == '__main__':
	main()
