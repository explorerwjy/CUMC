#!/home/local/users/jw/bin/python2.7
#read a vcf file and write as tsv file
from utils import *
import re
class Variant():
	def __init__(self,record,headers):
		self.record = record.strip()
		record=record.strip().split('\t')
		CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT = record[0:9]
		self.Chrom = CHROM  
		self.Pos = POS 
		self.Id = ID  
		self.Ref = REF 
		self.Alt = ALT 
		self.Qual = QUAL    
		self.Filter = FILTER  
		self.GetInfo(INFO)    
		self.Format = FORMAT
		self.GetSample(record,headers)
	
	def GetInfo(self,INFO):
		info = {}
		tmp = INFO.split(';')
		for item in tmp:
			item = item.split('=')
			if len(item) == 2:
				if item[0] in info:
					info[item[0]] = info[item[0]]+','+item[1]
				else:
					info[item[0]] = item[1]
			elif len(item) == 1:
				info[item[0]] = None
			else:
				print "Error Parsing INFO:", item
		self.Info = info
	
	def GetSample(self,record,headers):
		self.Sample = {}
		for i,sample in enumerate(headers[9:]):
			self.Sample[sample] = record[i+9]

	def show_samples(self):
		for k,v in self.Sample.items():
			print '%s:%v'%(k,v)


def get_header(line):
	return line.strip().split('\t')


def variant2tsv(variant,header,info_num):
	tmp = [variant.Chrom,variant.Pos,variant.Id,variant.Ref,variant.Alt,variant.Qual,variant.Filter]
	for field in header[7:info_num+7]:
		try:
			tmp.append(variant.Info[field])
		except KeyError:
			#print field
			tmp.append('.')
	tmp.append(variant.Format)
	for field in header[info_num+8:]:
		#print field
		try:
			tmp.append(variant.Sample[field])
		except KeyError:
			print field
			tmp.append('.')
	#print '\t'.join(tmp)+'\n'
	for i,item in enumerate(tmp):
		if item == None:
			tmp[i] = '.'
	return '\t'.join(tmp)+'\n'

def format_info(INFOs):
	res = []
	info = re.compile('##INFO=<ID=([^,]+),')
	for item in INFOs:
		res.append(info.match(item).group(1))
	#print res
	return res

Custome_info = ['Func.refGene','Gene.refGene','GeneDetail.refGene','ExonicFunc.refGene','AAChange.refGene','esp6500siv2_all','1000g2015aug_all','ExAC_ALL','MetaSVM_pred','CADD13_PHRED','MCAP','avsnp147','Annotation','GeneticModels','ModelScore','Compounds',]

def format_header(header,Infos):
	res = header[:7]
	res.extend(Infos)
	res.extend(header[8:])
	return res,len(Infos)

def GetOptions():
	parser=OptionParser()
	parser.add_option('-v','--vcf',dest='VCF',metavar='VCF',help='Input VCF file name')
	parser.add_option('-o','--outname',dest='OutName',metavar='OutName',help='Name of Output table')
	
	(options,args) = parser.parse_args()
	if options.OutName == None:
		options.OutName = GetBaseName(options.VCF)
	return options.VCF,options.OutName
def VCF2TABLE(vcf,out):
	fin = open(vcf,'rb')
	INFOs = []
	header = []
	for line in fin:
		if line.startswith('##'): #Mate Info
			if line.startswith('##INFO'): #Info
				INFOs.append(line)
		elif line.startswith('#C'): #header line, next line is variant.
			header = get_header(line)
			#Infos = format_info(INFOs)
			Infos = Custome_info 
			new_header,info_num = format_header(header,Infos)
			fout = open(out+'.tsv','wb')
			fout.write('\t'.join(new_header)+'\n')
		else:
			variant = Variant(line,header)
			fout.write(variant2tsv(variant,new_header,info_num))

def main():
	VCF,OUT=GetOptions()
	VCF2TABLE(VCF,OUT)
if __name__ == '__main__':
	main()
