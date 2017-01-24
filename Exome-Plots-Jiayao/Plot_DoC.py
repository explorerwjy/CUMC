#!/home/local/users/jw/bin/python2.7
#Author: jywang	explorerwjy@gmail.com

#========================================================================================================
# Plot The Depth of Cov
#========================================================================================================

from optparse import OptionParser
from utils import *
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

def GetOptions():
	parser = OptionParser()
	parser.add_option('-d','--dir',dest = 'DIR', metavar = 'DIR', help = 'Dirctory contains sample summary')
	parser.add_option('-m','--mode',dest='MODE',metavar='MODE',help = 'Mode to excute, 1, prepare input table; 2, plot D10-D15/Mean plot')
	parser.add_option('-t','--table',dest='Table',metavar='Table',help = 'DP Table to draw plot from')
	parser.add_option('-o','--out',dest='Outname',metavar='Outname',help = 'Outname')
	(options,args) = parser.parse_args()
	
	return options.MODE, options.DIR, options.Table, options.Outname
def Draw_Table(DIR,outname):
	sample_file = get_files(DIR,'DoC.sample_summary')
	fout = open(outname,'wb')
	headers = ['sample_id','total','mean','granular_third_quartile','granular_median','granular_first_quartile','%_bases_above_1','%_bases_above_5','%_bases_above_10','%_bases_above_15','%_bases_above_20']
	fout.write('\t'.join(headers)+'\n')
	for sample in sample_file:
		with open(sample) as fin:
			data = fin.readlines()[1]
			fout.write(data)
	fout.close()
def Plot_D15(Table):
	fin = open(Table,'rb')
	fin.readline()
	Means,D10,D15 = [],[],[]
	for l in fin:
		tmps = l.strip().split('\t')
		if float(tmps[8]) > 80:
			Means.append(float(tmps[2]))
			D10.append(float(tmps[8]))
			D15.append(float(tmps[9]))
	fin.close()
	Plot_scatter(Means,D10,D15)
def Plot_scatter(Means,D10,D15):
	#print Means
	#print Means,D10,D15
	#D10:
	#print "Plot"
	plt.scatter(Means,D10,color='r',marker='o',label='D10',alpha=1,s=2)
	plt.scatter(Means,D15,color='b',marker='o',label='D15',alpha=1,s=2)
	plt.xlabel('Mean Coverage')
	plt.ylabel('% Targeted Bases')
	plt.title('Mean-DP15')
	#plt.axis([0, 201, 0, 0.05])ppppp
	#plt.xticks(np.arange(0, 501, 10),rotation='vertical')
	plt.legend(loc='upper right')
	plt.grid(True)
	plt.show()

def make_key(chrom,start,end):
	if int(end) - int(start) > 1:
		return chrom+':'+str(int(start)+1)+'-'+end
	else:
		return chrom+':'+str(int(start)+1)

def Get_Genes(BedFile):
	fin = open(BedFile,'rb')
	last_gene = None
	res = []
	res_gene_dict = {}
	res_gene_list = []
	for l in fin:
		llist = l.strip().split('\t')
		chrom,start,end,gene = llist
		res_gene_dict[make_key(chrom,start,end)] = gene
		if gene != last_gene:
			if last_gene != None:
				res.append(One_gene)
			One_gene = [make_key(chrom,start,end)]
			res_gene_list.append(gene)
			last_gene = gene
		else:
			One_gene.append(make_key(chrom,start,end))
	res.append(One_gene)
	return res, res_gene_dict, res_gene_list

def slice_interv(intv_string):
	chrom,start_end = intv_string.split(":")
	start_end = start_end.split('-')
	if len(start_end) == 1:
		return int(start_end[0]),int(start_end[0])
	else:
		return int(start_end[0]),int(start_end[1])
# Intvs: list of Target:D15% pair
def Gene_Cov(Intvs):
	Lengths = []
	D15s = []
	for Intv in Intvs:
		start,end = slice_interv(Intv[0])

		length = end - start + 1
		Lengths.append(length)
		D15s.append(float(Intv[1]))

	TotalLen = sum(Lengths)
	TotalD15 = 0
	for length,D15 in zip(Lengths,D15s):
		TotalD15 += length*D15

	return float(TotalD15/TotalLen)

def Get_Sample_IntervStat(InpDir,Intv_Gene_Dict):
	sample_file = get_files(InpDir,'.DoC.sample_interval_summary')
	res = []
	sample_names = []
	for sample in sample_file:
		sample_names.append(sample.split('/')[-1].split('.')[0])
		One_sample_intervals = []
		One_sample_genes = []
		with open(sample) as fin:
			fin.readline() #Skip Header
			for l in fin:
				llist = l.strip().split('\t')
				#[0]:TargetKey [3]:TotalCov [4]:meanCov [11]:D15%
				One_sample_intervals.append([llist[0],llist[11]])
		#Reduce Interv by genes
		last_gene = None
		tmp = []
		for i in xrange(len(One_sample_intervals)):
			Intv,D15 = One_sample_intervals[i]
			gene = Intv_Gene_Dict[Intv]
			if gene != last_gene :
				if last_gene != None:
					#yield one gene result
					One_sample_genes.append(Gene_Cov(tmp))
					
				last_gene = gene
				tmp = [(Intv,D15)]
			else:
				tmp.append((Intv,D15))
		One_sample_genes.append(Gene_Cov(tmp))


		res.append(One_sample_genes)

	return sample_names,res

def Write(sample_names,Samples,res_gene_list):
	fout = open('DoC81_by_gene.tsv','wb')
	fout.write('ID'+'\t'+'\t'.join(res_gene_list)+'\n')
	for i in range(len(Samples)):
		fout.write(sample_names[i]+'\t'+'\t'.join(map(str,Samples[i]))+'\n')
	fout.close()

def ReadData():
	InpDir = './DoC81'
	BedFile = 'Modified_bed.bed'

	Gene_List,res_gene_dict, res_gene_list = Get_Genes(BedFile)
	#print Gene_List[0]
	sample_names,Samples = Get_Sample_IntervStat(InpDir,res_gene_dict)
	#print Samples[0][0]
	Write(sample_names,Samples,res_gene_list)

def PlotData():
	fin = open("DoC81_by_gene.tsv",'rb')
	Data = []
	genes = fin.readline().strip().split('\t')[1:]
	IDs = []
	for l in fin:
		llist = l.strip().split('\t')
		IDs.append(llist[0])
		Data.append(llist[1:])
	"""
	for i in xrange(len(Data)):
		for j in xrange(len(Data[0])):
	"""
	print "GeneName\t%D15>90\t%D15>80\t%D15>70\t%D15>60\t%D15>50"
	for i in xrange(len(Data[0])):
		Gene = []
		for j in xrange(len(Data)):
			Gene.append(float(Data[j][i]))
		#print Gene
		Stat_Gene(Gene,genes[i])
		#Plot_histo(Gene,genes[i])

def Stat_Gene(Gene,gene_name):
	bins = [0,0,0,0,0]
	for sample in Gene:
		if sample > 90:
			bins[0] += 1
		if sample > 80:
			bins[1] += 1
		if sample > 70:
			bins[2] += 1
		if sample > 60:
			bins[3] += 1
		if sample > 50:
			bins[4] += 1
	print gene_name+'\t'+"\t".join(map(str,bins))


def Plot_histo(Gene,gene_name):
	#n, bins, patches = plt.hist(Gene, 50, normed=1, facecolor='green', alpha=0.75)
	#l = plt.plot(bins, y, 'r--', linewidth=1)
	
	plt.clf()

	bins=[0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
	plt.hist(Gene, bins=bins)

	plt.xlabel(r'%D15')
	plt.ylabel('Num of Samples')
	plt.title('D15 of Gene: {}'.format(gene_name))
	plt.xticks(np.arange(0, 100, 10),rotation='vertical')
	plt.grid(True)

	#plt.show()
	plt.savefig('{}.png'.format(gene_name))

def main():
	#ReadData()
	PlotData()

if __name__=='__main__':
	main()
