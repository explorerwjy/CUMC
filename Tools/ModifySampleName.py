#!/home/local/users/jw/anaconda2/bin/python
#Author: jywang	explorerwjy@gmail.com

#========================================================================================================
# ModifySampleName.py
#========================================================================================================

import argparse
import re
import gzip

def GetOptions():
	parser = argparse.ArgumentParser()
	parser.add_argument('-v','--vcf',required=True, type=str, help = 'VCF file to be modify sample names')
	parser.add_argument('-o','--out',required=True, type=str, help = 'Output name of modify sample names VCF')
	parser.add_argument('-m','--mapping', type=str, help = 'ID Mapping File')
	args = parser.parse_args()
	
	return args.vcf, args.out, args.mapping

def ModifySampleID(samples, IDMapping):
	for i in range(len(samples)):
		#P = N0.search(samples[i]).group(1)
		#samples[i] = P
		samples[i] = IDMapping[samples[i]]

def Load_ID_Mapping(IDMappingFil):
	fin = open(IDMappingFil, 'rb')
	res = {}
	for l in fin:
		k,v = l.strip().split(',')
		res[k] = v
	return res

def main():
	VCF, OUT, IDMappingFil = GetOptions()
	if VCF.endswith(".gz"):
		fin = gzip.open(VCF, 'rb')
	else:
		fin = open(VCF, 'rb')
	fout = open(OUT+'.vcf', 'wb')
	IDMapping = Load_ID_Mapping(IDMappingFil)
	for l in fin:
		if l.startswith('##'):
			fout.write(l)
		elif l.startswith('#'):
			#fout.write('##ModifySampleNameWithPythonScript ModifySampleName.py\n')
			headers = l.strip().split('\t')
			Part1 = headers[:9]
			samples = headers[9:]
			#print Part1
			ModifySampleID(samples, IDMapping)
			#print samples
			tmp = Part1 + samples
			fout.write('\t'.join(tmp)+'\n')
		else:
			fout.write(l)
	return

if __name__=='__main__':
	main()
