#!/home/local/users/jw/bin/python2.7
#Author: jywang	explorerwjy@gmail.com

#========================================================================================================
# Make FastqTable According to list of fastq names
#========================================================================================================

from optparse import OptionParser

def get_name(path):
	return path.split('/')[-1]
def search_mate(fq_name,fq_list):
	fq_name = fq_name.split('_')
	for j in xrange(1,len(fq_list)):
		mate_name = get_name(fq_list[j])
		mate_name = mate_name.split('_')
		flag = True
		for k in [0,1,2,4]:
			if fq_name[k] != mate_name[k]:
				flag = False
				break
		if flag == True: # Find Mate
			if fq_name[3] == 'R1' and mate_name[3] == 'R2':
				return fq_list[0], fq_list[j] , j
			elif fq_name[3] == 'R2' and mate_name[3] == 'R1':
				return fq_list[j], fq_list[0] , j 
			else:
				print 'Error with File Name','_'.join(fq_name), '_'.join(mate_name)
				exit()
	#Didn't find a mate:
	return '_'.join(fq_name), None, None
def get_RG(fq_name):
	fq_name = get_name(fq_name).split('.')[0]
	return '{}\\t{}\\t{}\\t{}\\t{}\\t{}'.format('@RG','ID:'+fq_name,'SM:'+fq_name,'LB:'+fq_name,'PL:'+fq_name,'CN:'+fq_name)
def MakePair(fq_list):
	while len(fq_list) > 0:	
		fq_path = fq_list[0]
		fq_name = get_name(fq_path)
		R1, R2, j = search_mate(fq_name, fq_list)
		RG = get_RG(fq_list[0])
		if R2 != None:	
			yield R1 + '\t' + RG + '\t' + R2 + '\n'
			fq_list.pop(j)
			fq_list.pop(0)
		else:
			yield R1 + '\t' + RG + '\t' + '\n'

def MakeFastqTable(Input, Output):
	fq_list = [ tmp.strip() for tmp in open(Input, 'rb').readlines() ]
	fout = open(Output, 'wb')
	for record in MakePair(fq_list):
		fout.write(record)
	fout.close()


def GetOptions():
	parser = OptionParser()
	parser.add_option('-i','--input',dest = 'input', metavar = 'input', help = 'Input fastq file listSRA')
	parser.add_option('-o','--output',dest = 'output', metavar = 'output', help = 'Output fastq table name')

	(options,args) = parser.parse_args()
	if options.output == None:
		options.output = 'Fastq_table.txt'
	return options.input, options.output

def main():
	Input, Output = GetOptions()
	MakeFastqTable(Input, Output)

if __name__=='__main__':
	main()
