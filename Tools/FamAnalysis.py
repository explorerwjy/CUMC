#!/home/local/users/jw/anaconda2/bin/python
#Author: jywang	explorerwjy@gmail.com

#========================================================================================================
# FamAnalysis.py
# 
#========================================================================================================

import argparse
import os
import subprocess

PrepareVCF = '/home/local/users/jw/CUMC/Tools/PSAP_Family.sh'

def GetOptions():
	parser = argparse.ArgumentParser()
	parser.add_argument('-v','--vcf', type=str, help = 'VCF File')
	parser.add_argument('-p','--ped', type=str, help = 'Pedigree File')
	args = parser.parse_args()
	
	return args.vcf, args.ped

def CleanPed():
	fin = open(TXT,'rb')
	fout = open(PED, 'wb')
	for l in fin:
		l = l.strip() + '\n'
		fout.write(l)
	fin.close()
	fout.close()

def GetIndividualList(Path):
	with open(Path) as fin:
		for l in fin:
			if l.startswith('##'):
				continue
			elif l.startswith('#'):
				return l.strip().split('\t')[9:]

def ParsePed(Ped, Path):
	Samples = GetIndividualList(Path)
	WorkDir = os.getcwd()
	fin = open(Ped, 'rb')
	NewFam = 'None'
	for l in fin:
		if l.startswith('#'):
			header = l
			continue
		llist = l.strip().split()
		FamID = llist[0]
		if FamID != NewFam:
			if NewFam != 'None':
				#print NewFam
				#showFam(Fam)
				ProcessOneFam(WorkDir, Samples, NewFam, Fam, Path)
				NewFam = FamID
				Fam = [header]
				Fam.append(l)
			else:
				NewFam = FamID
				Fam = [header]
				Fam.append(l)
		else:
			Fam.append(l)
	ProcessOneFam(WorkDir, Samples, NewFam, Fam, Path)
	#CheckAndSeperate('{}.list'.format(NewFam), Mapper)		
	fin.close()

def ProcessOneFam(WorkDir, SampleList, NewFam, Fam, Path):
	for item in printSample(Fam):
		if item not in SampleList:
			os.chdir(WorkDir)
			return
	print NewFam
	Path = WorkDir + '/' + Path
	subdir = WorkDir + '/' + NewFam
	if not os.path.isdir(subdir):
		os.mkdir(subdir)
	os.chdir(subdir)
	fout = open('{}.ped'.format(NewFam), 'wb')	
	for item in Fam:
		fout.write(item)
	fout.close()
	fout = open('{}.list'.format(NewFam), 'wb')	
	for item in printSample(Fam):
		fout.write(item+'\n')
	fout.close()
	SampleList = '{}.list'.format(NewFam)
	if Path != None:
		cmd = 'nohup {} {} {} &'.format(PrepareVCF, Path, SampleList)
		print cmd
		process = subprocess.Popen(cmd, shell=True)
	os.chdir(WorkDir)

def showFam(Fam):
	for l in Fam:
		print l.split('\t')[0:2]
def printSample(Fam):
	for l in Fam:
		if l.startswith('#'):
			continue
		else:
			yield l.strip().split('\t')[1]
			#res.append(l.strip().split('\t')[1])

def CheckAndSeperate(LIST, Mapper):
	# Check
	try:	
		List = [ x.strip() for x in open(LIST,'rb').readlines() ]
		Path = set([])
		for item in List:	
			Path.add(Mapper[item])
		if len(Path) != 1:
			print 'Error with',LIST
			return None
		return Mapper[item] 
	except:
		print 'No {} Sample Here {}'.format(item, LIST)
			

def SeperatePedMkdir(PED, VCF):
	Peds = ParsePed(PED, VCF)

def main():
	#CleanPed()
	VCF, PED = GetOptions()
	SeperatePedMkdir(PED, VCF)
	return

if __name__=='__main__':
	main()
