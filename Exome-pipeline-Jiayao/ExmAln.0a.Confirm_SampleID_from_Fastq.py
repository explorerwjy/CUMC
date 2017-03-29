#!/home/yufengshen/anaconda2/bin/python
# Author: jywang	explorerwjy@gmail.com

#=========================================================================
# ExmAln.0a.Confirm_SampleID_from_Fastq.py
#=========================================================================

import argparse
import re

name = re.compile('(OMG\d+-\d+-[A-Za-z0-9-]+)') #PIPseq ID

class Sample:
    def __init__(self, ID, F_path):
        self.SampleID = ID
        self.F_paths = [F_path]
    def addPath(self, F_path):
        self.F_paths.append(F_path)

def GetOptions():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fastq', type=str, help='Fastq List')
    args = parser.parse_args()

    return args.fastq

def ExtractSample(FastqList):
    res = {}
    fin = open(FastqList, 'rb')
    for l in fin:
        fpath = l.strip()
        fname = fpath.split('/')[-1]
        sample = name.search(fname).group(1)
        #print sample
        if sample in res:
            res[sample].addPath(fpath)
        else:
            res[sample] = Sample(sample, fpath)
    fin.close()
    for sample in sorted(res.values(),key=lambda x:x.SampleID):
        print sample.SampleID

def main():
    FastqList = GetOptions()
    ExtractSample(FastqList)
    return


if __name__ == '__main__':
    main()
