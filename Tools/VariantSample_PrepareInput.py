#!/home/yufengshen/anaconda2/bin/python
# Author: jywang	explorerwjy@gmail.com

#=========================================================================
# VariantSample_PrepareIGVInput.py
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

SampleName = re.compile('([A-Za-z0-9-]+).bam')
EXM_REF = '$HOME/CUMC/Exome-pipeline-Jiayao/WES_Pipeline_References.b37.biocluster.sh'
BAMOUT_CMD = '$HOME/CUMC/Exome-pipeline-Jiayao/ExmAdHoc.15b.BamOut.sh'
IGV_GENERATE_CMD = '$HOME/CUMC/Exome-IGV-Jiayao/Generate_IGV_plots.sh'
RUN = 'nohup'


def GetOptions():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--var', type=str, required=True,
                        help='Input Variant Table, "chr\tpos\tSampleID" format.')
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
    parser.add_argument('--parallel', type=int,
                        default=20, help='Number of Parallel to go. ')
    args = parser.parse_args()
    if args.list1 == None:
        args.list1 = 'Variants.bamout.input'
    if args.list2 == None:
        args.list2 = 'Variants.IGV.varlist'
    if args.script1 == None:
        args.script1 = 'run_bamout.sh'
    if args.script2 == None:
        args.script2 = 'run_IGV.sh'
    return args


class Variant_PrepareIGVInput:
    def __init__(self, args):
        self.VARFil = os.path.abspath(args.var)
        self.Bams = BamLocation(args.bam)
        self.VarList1 = os.path.abspath(args.list1)
        self.VarList1_hand = open(self.VarList1, 'wb')
        self.VarList2 = os.path.abspath(args.list2)
        self.VarList2_hand = open(self.VarList2, 'wb')
        self.Script1 = open(args.script1, 'wb')
        self.Script2 = open(args.script2, 'wb')
        self.parallel = args.parallel

    def run(self):
        self.variants = []
        with open(self.VARFil, 'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter="\t")
            header = reader.next()
            for row in reader:
                var = Variant(header, row)
                var.addBam(self.Bams)
                self.variants.append(var)
                self.VarList2_hand.write(var.OutVar())
        self.FlushBamout()
        self.FlushIGV()

    def ReduceVarian2Sample(self):
        res = {}
        for variant in self.variants:
            for sample in variant.Samples:
                if sample not in res:
                    res[sample] = Sample(variant)
                else:
                    res[sample].AddVar(variant)
        #print res
        #print res.values()
        return res.values()

    def FlushBamout(self):
        # Write Bamout List
        Samples = self.ReduceVarian2Sample()
        for sample in Samples:
            #print sample.
            for line in sample.FlushBamout():
                self.VarList1_hand.write(line)
        # Write Bamout Script
        fout = self.Script1
        fout.write('REF={}\n'.format(EXM_REF))
        fout.write('CMD={}\n'.format(BAMOUT_CMD))
        fout.write('InpFil={}\n'.format(self.VarList1))
        fout.write("NumJob=`wc -l $InpFil|cut -f1 -d ' '`\n\n")
        fout.write("mkdir -p Bamouts;cd Bamouts\n")
        fout.write(
            "seq $NumJob| parallel -j {} --eta $CMD -i $InpFil -a {} -r $REF\n\n".format(self.parallel, "{}"))
        fout.write(
            "cd ../;find `pwd`/Bamouts/ -name \"*.bamout.bam\" > Bamout.bam.list\n")
        fout.write("cat {} {} > ALL.bam.list\n".format(
            self.Bams.bamlocationfile, "Bamout.bam.list"))
        fout.close()

    def FlushIGV(self):
        fout = self.Script2
        fout.write(
            '{} -b ALL.bam.list -v {}'.format(IGV_GENERATE_CMD, self.VarList2))
        fout.close()


class Sample:
    def __init__(self, variant):
        self.Sample = variant.Samples
        self.SampleBam = variant.SampleBam
        self.SampleBamout = variant.SampleBamout
        # self.variants = []
        chrom, start, end = variant.Chrom, variant.start, variant.end
        self.Intervals = ['{}:{}-{}'.format(str(chrom), str(start), str(end))]

    def AddVar(self, variant):
        # self.vairants.append(variant)
        chrom, start, end = variant.Chrom, variant.start, variant.end
        self.Intervals.append(
            '{}:{}-{}'.format(str(chrom), str(start), str(end)))

    def FlushBamout(self):
        # Format is
        # BAMPATH	-L	VAR1	-L VAR2	 .... -L VARN
        # SampleCMD =  ''.format(self.SampleBam.FullPath, '\t'.join(['-L '+L for L in self.Intervals]))
        SampleCMD = []
        for bam in self.SampleBam:
            print bam.FullPath
            SampleCMD.append('{}\t{} \n'.format(bam.FullPath, '\t'.join([
                                       '-L ' + L for L in self.Intervals])))
        return SampleCMD


class BAM:
    def __init__(self, fpath):
        self.FullPath = fpath.strip()
        self.BamName = self.FullPath.split('/')[-1]
        self.BamOutName = '.'.join(
            self.BamName.split('.')[:-1]) + ".bamout.bam"
        #self.Sample = SampleName.search(self.BamName).group(1)
        self.Sample = SampleName.search(self.BamName).group(1)


class BamLocation:
    def __init__(self, bamlocationfile):
        self.bamlocationfile = os.path.abspath(bamlocationfile)
        fin = open(bamlocationfile, 'rb')
        self.Bams = {}
        for l in fin:
            if 'bamout' in l:
                continue
            bam = BAM(l)
            self.Bams[bam.BamName] = bam

    def search(self, sample):
        for k, v in self.Bams.items():
            if sample in v.FullPath:
                return v
        return None


class Variant:
    def __init__(self, header, row):
        #print header
        self.Chrom = row[header.index('Chrom')]
        self.Pos = int(row[header.index('Pos')])
        self.Samples = row[header.index('Sample')].split(",")
        self.start = max(0, self.Pos - 150)
        self.end = self.Pos + 150
        self.SampleBam = []
        self.SampleBamout = []
    def addBam(self, bamlocation):
        # print self.Sample
        for bam in bamlocation.Bams.values():
            for sample in self.Samples:
                if sample == bam.Sample:
                    self.SampleBam.append(bam)
                    self.SampleBamout.append(bam.BamName.rstrip('.bam') + '.bamout.bam')
        #try:
        #for var1,var2 in zip(self.SampleBam, self.SampleBamout):
        #    if var1 not in locals() or var2 not in locals():
        #        print 'Cant find {} in locals.'.format(var1)
        #        print 'Cant find {} in locals.'.format(var2)
        #        exit()
        #except:
        #    print self.Samples

    def OutVar(self):
        #return '{}\t{}\t{}\n'.format(self.Chrom, self.Pos, ','.join([self.SampleBam.BamName, self.SampleBamout]))
        pairs = []
        for a,b in zip(self.SampleBam, self.SampleBamout):
            pairs.append(a.BamName + "," + b)
        return '{}\t{}\t{}\n'.format(self.Chrom, self.Pos, ','.join(pairs))

    def Out2L(self):  # To run_Bamout.sh
        return "-L {}:{}-{}".format(self.Chrom, self.start, self.end)

    def Out2I(self):  # To run_IGV.sh
        return "{}\t{}".format(self.Chrom, self.Pos)

def main():
    args = GetOptions()
    instance = Variant_PrepareIGVInput(args)
    instance.run()
    return


if __name__ == '__main__':
    main()
