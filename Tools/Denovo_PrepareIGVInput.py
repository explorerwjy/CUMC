#!/home/yufengshen/anaconda2/bin/python
# Author: jywang	explorerwjy@gmail.com

#=========================================================================
# DNM_PrepareIGVInput.py
# Prepare bamout script and variant table for IGV shot
# Input: Dir contain subdirs of PSAP results, (PSAP.csv, .ped file is need)
#        TXT file contains (sampleID,BAM) tuple
# Output: List of (sampleID Interval1 Interval2 ... Inverbaln)
#         Scripts to run BAMOUT.
#         Scripts to run IGV
#=========================================================================

import argparse
import csv
import re
import os


ProjectName = re.compile('\w+')
SampleName = re.compile('(COL-CHUNG_[a-zA-z0-9]+_Diabetes_[\d]+)*.bam')
SampleName = re.compile('(CARE[a-zA-Z0-9-]+).bam')
ProjectName = re.compile('\w+')
EXM_REF = '$HOME/CUMC/Exome-pipeline-Jiayao/WES_Pipeline_References.b37.sh'
BAMOUT_CMD = '$HOME/CUMC/Exome-pipeline-Jiayao/ExmAdHoc.15b.BamOut.sh'
IGV_GENERATE_CMD = '$HOME/CUMC/Exome-IGV-Jiayao/Generate_IGV_plots.sh'
RUN = 'nohup'

def GetOptions():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dir', type=str, required=True,
                        help='Directory contains PSAP by fam results, each subdir should have .csv and .ped')
    parser.add_argument('-b', '--bam', type=str, required=True,
                        help='File contains bam locations. ')
    parser.add_argument('-l1', '--list1', type=str,
                        help='list of sample and corresponding intervals, For Bamout')
    parser.add_argument('-l2', '--list2', type=str,
                        help='list of variants and corresponding bam, For IGV')
    parser.add_argument('-s1', '--script1', type=str,
                        help='Script to run GATK Bamout. ')
    parser.add_argument('-s2', '--script2', type=str,
                        help='Script to run IGV. ')
    parser.add_argument('-p', '--parallel', type=int,
                        default=20, help='Number of Parallel to go. ')
    args = parser.parse_args()
    args.dir = os.path.abspath(args.dir)
    print args.dir
    if args.list1 == None:
        args.list1 = ProjectName.search(
            args.dir.split('/')[-1]).group(0) + '.samplelist'
    if args.list2 == None:
        args.list2 = ProjectName.search(
            args.dir.split('/')[-1]).group(0) + '.varlist'
    if args.script1 == None:
        args.script1 = 'run_bamout.sh'
    if args.script2 == None:
        args.script2 = 'run_IGV.sh'
    return args

class DNM_PrepareIGVInput:
	#def __init__(self, Dir, BAMList, VarList1, VarList2, Script1, Script2):
	def __init(self, args):
		self.Dir = os.path.abspath(args.dir)
		self.bampath = BamLocation(args.bam)
		self.VarList1 = args.l1  
		self.VarList1_hand = open(self.VarList1, 'wb')
		self.VarList2 = args.l2 
		self.VarList2_hand = open(self.VarList2, 'wb')
		self.Script1 = open(args.s1 , 'wb')
		self.Script2 = open(args.s2 , 'wb')

	def run(self):
		pedigree = Pedigree(ped)
		Bams = BamLocation(bam)
		with open(_csv, 'rb') as csvfile:
			reader = csv.reader(csvfile)
			header = reader.next()
			for row in reader:
				var = Variant(header, row)
				var.addPed(pedigree)
				var.addBam(Bams)
				res.append(var)
				fout.write(var.OutVar())
		fout.close()
		return res
		self.FlushBamout()
		self.FlushIGV()

	def LoadCSV(self, _csv, ped, bam, output):
		pedigree = Pedigree(ped)
		Bams = BamLocation(bam)
		res = []
		fout = open(output,'wb')
		with open(_csv, 'rb') as csvfile:
			reader = csv.reader(csvfile)
			header = reader.next()
			for row in reader:
				var = Variant(header, row)
				var.addPed(pedigree)
				var.addBam(Bams)
				res.append(var)
				fout.write(var.OutVar())
		fout.close()
		return res

	def ReduceVarian2Sample(self):
		res = {}
		for variant in variants:
			if variant.Sample not in res:
				res[variant.Sample] = Sample(variant)
			else:
				res[variant.Sample].AddVar(variant)
		return res.values()

	def FlushBamout(self):
        fout = self.Script1
        fout.write('REF={}\n'.format(EXM_REF))
        fout.write('CMD={}\n'.format(BAMOUT_CMD))
        fout.write('InpFil={}\n'.format(self.VarList1))
        fout.write("NumJob=`wc -l $InpFil|cut -f1 -d ' '`\n\n")
        fout.write("mkdir -p Bamouts;cd Bamouts\n")
        fout.write(
            "seq $NumJob| parallel -j {} --eta $CMD -i $InpFil -a {} -r $REF\n\n".format(self.parallel, "{}"))
        fout.write("cd ../;find `pwd`/Bamouts/ -name \"*.bamout.bam\" > Bamout.bam.list\n")
        fout.write("cat {} {} > ALL.bam.list\n".format(
            self.bampath.bamlocationfile, "Bamout.bam.list"))
        fout.close()

	def FlushIGV(self):
        fout = self.Script2
        fout.write(
            '{} -b ALL.bam.list -v {}'.format(IGV_GENERATE_CMD, self.VarList2))
        fout.close()

class BAM:
	def __init__(self, fpath):
		self.FullPath = fpath.strip()
		self.BamName = self.FullPath.split('/')[-1]
		self.BamOutName = '.'.join(self.BamName.split('.')[:-1]) + ".bamout.bam"
		self.Sample = SampleName.search(self.BamName).group(1)
		#self.Sample = self.BamName
 
class BamLocation:
	def __init__(self, bamlocationfile):
		fin = open(bamlocationfile, 'rb')
		self.Bams = {}
		for l in fin:
			if 'bamout' in l:
				continue
			bam = BAM(l)
			self.Bams[bam.BamName] = bam
	def search(self,sample):
		for k,v in self.Bams.items():
			if sample in v.FullPath:
				return v
		return None

class Variant:
	def __init__(self, header, row):
		self.Chrom = row[header.index('Chrom')]
		self.Pos = int(row[header.index('Pos')])
		self.Sample = row[header.index('Sample')]
		self.start = max(0, self.Pos - 150)
		self.end = self.Pos + 150

	def addPed(self, pedigree):
		self.Father = pedigree.Probands[self.Sample].Father
		self.Mother = pedigree.Probands[self.Sample].Mother

	def addBam(self, bamlocation):
		print self.Sample
		for bam in bamlocation.Bams.values():
			if self.Sample == bam.Sample:
				self.SampleBam = bam
				self.SampleBamout = bam.BamName.rstrip('.bam')+'.bamout.bam'
			if self.Father == bam.Sample:
				self.FatherBam = bam
				self.FatherBamout = bam.BamName.rstrip('.bam')+'.bamout.bam'
			if self.Mother == bam.Sample:
				self.MotherBam = bam
				self.MotherBamout = bam.BamName.rstrip('.bam')+'.bamout.bam'
		try:
			for var in [self.SampleBam, self.SampleBamout, self.FatherBam, self.FatherBamout, self.MotherBam, self.MotherBamout ]:
				if var not in locals():
					print 'Cant find {} in locals.'.format(var)
					exit()
		except:
			print self.Sample
	def OutVar(self):
		return '{}\t{}\t{}\n'.format(self.Chrom, self.Pos, ','.join([self.SampleBam.BamName, self.SampleBamout, self.FatherBam.BamName, self.FatherBamout, self.MotherBam.BamName, self.MotherBamout]))
	def Out2L(self):
		return "-L {}:{}-{}".format(self.Chrom, self.start, self.end)
	def Out2I(self):
		return "{}\t{}".format(self.Chrom, self.Pos)

def main():
	args = GetOptions()
	instance = DNM_PrepareIGVInput(args)
	instance.run()
	return


if __name__ == '__main__':
	main()
