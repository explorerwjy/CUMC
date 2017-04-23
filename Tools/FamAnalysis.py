#!/home/local/users/jw/anaconda2/bin/python
#Author: jywang	explorerwjy@gmail.com

#========================================================================================================
# FamAnalysis.py
# 
#========================================================================================================

import argparse
import os
import subprocess

PrepareVCF = '/home/local/users/jw/CUMC/Tools/GetFamVCF.sh'
PED = '/home/local/users/jw/CARE/hg19/CARE.hg19.Annotated.ped'
IDMapper = ''
VCF = '/home/local/users/jw/CARE/hg19/CARE.hg19.RareCoding.vcf'

def GetOptions():
	parser = argparse.ArgumentParser()
	parser.add_argument('-','--', type=str, help = '')
	args = parser.parse_args()
	
	return

def CleanPed():
	fin = open(TXT,'rb')
	fout = open(PED, 'wb')
	for l in fin:
		l = l.strip() + '\n'
		fout.write(l)
	fin.close()
	fout.close()

def ParsePed(Ped):
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
				#Path = CheckAndSeperate('{}.list'.format(NewFam), Mapper)		
				Path = VCF
				SampleList = '{}.list'.format(NewFam)
				if Path != None:
					cmd = 'nohup {} {} {} &'.format(PrepareVCF, Path, SampleList)
					process = subprocess.Popen(cmd, shell=True)
				os.chdir(WorkDir)

				NewFam = FamID
				Fam = [header]
				Fam.append(l)
			else:
				NewFam = FamID
				Fam = [header]
				Fam.append(l)
		else:
			Fam.append(l)
	
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
	Path = VCF
	SampleList = '{}.list'.format(NewFam)
	if Path != None:
		cmd = 'nohup {} {} {} &'.format(PrepareVCF, Path, SampleList)
		process = subprocess.Popen(cmd, shell=True)
	os.chdir(WorkDir)
	#CheckAndSeperate('{}.list'.format(NewFam), Mapper)		
	fout.close()
	fin.close()

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
			

def SeperatePedMkdir():
	Peds = ParsePed(PED)

def main():
	#CleanPed()
	SeperatePedMkdir()
	return

if __name__=='__main__':
	main()
