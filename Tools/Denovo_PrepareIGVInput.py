#!/home/yufengshen/anaconda2/bin/python
# Author: jywang	explorerwjy@gmail.com

#=========================================================================
# PSAP_PrepareInput.py
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
EXM_REF = '$HOME/CUMC/Exome-pipeline-Jiayao/WES_Pipeline_References.b37.sh'
BAMOUT_CMD = '$HOME/CUMC/Exome-pipeline-Jiayao/ExmAdHoc.15b.BamOut.sh'
RUN = 'nohup'

def GetOptions():
	parser = argparse.ArgumentParser()
	parser.add_argument('-d', '--dir', type=str, required=True, help='Directory contains PSAP by fam results, each subdir should have .csv and .ped')
	parser.add_argument('-b', '--bam', type=str, required=True, help='File contains bam locations. ')
	parser.add_argument('-l1', '--list1', type=str, help='list of sample and corresponding intervals, For Bamout')
	parser.add_argument('-l2', '--list2', type=str, help='list of variants and corresponding bam, For IGV')
	parser.add_argument('-s1', '--script1', type=str, help='Script to run GATK Bamout. ')
	parser.add_argument('-s2', '--script2', type=str, help='Script to run IGV. ')
	args = parser.parse_args()
	args.dir = os.path.abspath(args.dir)
	print args.dir
	if args.list1 == None:
		args.list1 = ProjectName.search(args.dir.split('/')[-1]).group(0) + '.samplelist'
	if args.list2 == None:
		args.list2 = ProjectName.search(args.dir.split('/')[-1]).group(0) + '.varlist'
	if args.script1 == None:
		args.script1 = 'run_bamout.sh'
	if args.script2 == None:
		args.script2 = 'run_IGV.sh'
	return args.dir, args.bam, args.list1, args.list2, args.script1, args.script2

class Denovo_Prepare_IGV_Input:
	def __init__(self, Dir, BAMList, VarList1, VarList2, Script1, Script2):
		self.Dir = os.path.abspath(Dir)
		self.bampath = BamLocation(BAMList)
		self.VarList1 = VarList1 
		self.VarList1_hand = open(self.VarList1, 'wb')
		self.VarList2 = VarList2
		self.VarList2_hand = open(self.VarList2, 'wb')
		self.Script1 = open(Script1, 'wb')
		self.Script2 = open(Script2, 'wb')

	def run(self):
		dirList = [self.Dir + '/' + x for x in os.listdir(self.Dir)]
		for subdir in dirList:
			if not os.path.isdir(subdir):
				continue
			subList = [ subdir+'/' + x for x in os.listdir(subdir)]
			self.OnePed(subList)
		self.FlushBamout()
		#self.FlushIGV()

	def OnePed(self, subList):
		PedFil = None
		PSAPFil = None
		for _file in subList:
			if _file.endswith('.ped'):
				PedFil = _file
			if _file.endswith('.csv'):
				PSAPFil = _file
		if PedFil == None or PSAPFil == None:
			return None
		print PedFil.split('/')[-1], PSAPFil.split('/')[-1]
		self.ReadCSV(PSAPFil, PedFil)

	def ReadCSV(self, PSAPFil, ped):
		pedigree = Pedigree(ped)
		Proband = pedigree.GuessProband()
		Bam = self.bampath.search(Proband)
		res = []
		with open(PSAPFil, 'rb') as csvfile:
			reader = csv.reader(csvfile)
			header = reader.next()
			idx_Chr = header.index("Chrom")
			idx_Pos = header.index("Pos")
			idx_Ref = header.index("Ref")
			idx_Alt = header.index("Alt")
			idx_Proband = header.index(Proband)
			idx_DModel = header.index("Dz.Model."+Proband)
			idx_PopScore = header.index("popScore."+Proband)
			for row in reader:
				var = Variant(row, idx_Chr, idx_Pos, idx_Ref, idx_Alt, idx_Proband, idx_DModel, idx_PopScore)
				if var.needIGV():
					res.append(var.Out2L())
					self.VarList2_hand.write("{}\t{},{}\n".format(var.Out2I(), Bam.BamName, Bam.BamOutName))
			self.VarList1_hand.write("{}\t{}\n".format(Bam.FullPath, '\t'.join(res)))
		return res

	def FlushBamout(self):
		fout = self.Script1
		fout.write('REF={}\n'.format(EXM_REF))
		fout.write('CMD={}\n'.format(BAMOUT_CMD)) 
		fout.write('InpFil={}\n'.format(self.VarList1)) 
		fout.write("NumJob=`wc -l $InpFil|cut -f1 -d ' '`\n\n")
		fout.write("seq $NumJob| parallel -j 30 --eta $CMD -i $InpFil -a {} -r $REF\n\n")
		fout.write("find `pwd` -name \"*.bamout.bam\" > Bamout.bam.list")
		fout.close()
	def FlushIGV(self):
		fout = self.Script2
		fout.write(''.format())

class Indv:
	def __init__(self,FamID, SampleID, FatherID, MotherID, Sex, Affected):
		self.FamID = FamID
		self.SampleID = SampleID
		self.FatherID = FatherID
		self.MotherID = MotherID
		self.Sex = Sex
		self.Affected = Affected
	def show(self):
		print '\t'.join([self.FamID,self.SampleID, self.FatherID, self.MotherID, self.Sex, self.Affected])

#Type: Singleton, Trios and complex Fam.
class Pedigree:
	def __init__(self, pedfile):
		fin = open(pedfile, 'rb')
		self.Indvs = []
		for l in fin:
			if l.startswith('#'):
				continue
			llist = l.strip().split('\t')
			FamID, SampleID, FatherID, MotherID, Sex, Affected = llist[:6]
			indv = Indv(FamID, SampleID, FatherID, MotherID, Sex, Affected)
			self.Indvs.append(indv)
		self.FamID = FamID
	def GuessProband(self):
		for Indv in self.Indvs:
			Indv.show()
			if Indv.Affected == '2':
				return Indv.SampleID
		print self.FamID, "Can't find a Proband"
		exit()
		return self.Indvs[0].SampleID

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
	def __init__(self, row, idx_Chr, idx_Pos, idx_Ref, idx_Alt, idx_Proband, idx_DModel, idx_PopScore):
		self.Chrom = row[idx_Chr]
		self.Pos = int(row[idx_Pos])
		self.Ref = row[idx_Ref]
		self.Alt = row[idx_Alt].split(',')
		self.Sample = row[idx_Proband]
		self.DModel = row[idx_DModel]
		self.PopScore = row[idx_PopScore]
		self.start = max(0, self.Pos - 150)
		self.end = self.Pos + 150
	# INDELs and Chet with top PopScore needs IGV comfirm.
    def addPed(self, pedigree):
        self.Father = pedigree.Probands[self.Sample].Father
        self.Mother = pedigree.Probands[self.Sample].Mother
    def addBam(self, bamlocation ):
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
	Dir, Bam, List1, List2, script1, script2 = GetOptions()
	instance = Denovo_Prepare_IGV_Input(Dir, Bam, List1, List2, script1, script2)
	instance.run()
	return


if __name__ == '__main__':
	main()
