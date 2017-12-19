#!~/anaconda2/bin/python
# Author: jywang	explorerwjy@gmail.com

#=========================================================================
# Filters on VCF file
# Select Biallelic High Quaily Common Variants For PLINK IBD and KING
#=========================================================================

import argparse
import gzip
import os
import re
import sys

def GetOptions():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', required=True, type=str, help="Input VCF file")
    parser.add_argument('-o', '--out', type=str, help="Output Filtered VCF file")
    parser.add_argument('-f', '--freq', type=float, default=5e-3, help="WGS AF greater than freq (default = 5e-3)")
    args = parser.parse_args()
    return args

class WGS_Common:
    def __init__(self, args):
        self.vcf = args.vcf
        self.out = args.out
        self.freq = args.freq
        if self.out == None:
            self.out = self.vcf.split("/")[-1].rstrip(".gz").rstrip(".vcf") + ".wgscommon.vcf"
        elif self.out == 'stdout':
            self.out = sys.stdout
    def run(self):
        if self.vcf.endswith('.gz'):
            hand = gzip.open(self.vcf, 'rb')
        else:
            hand = open(self.vcf, 'rb')
        fout = open(self.out, 'wb')
        Count_All = 0
        Count_Pass = 0
        for l in hand:
            if l.startswith('##'):
                fout.write(l)
                continue
            elif l.startswith('#'):
                fout.write(l)
                continue
            llist = l.strip().split('\t')
            Count_All += 1
            if self.F_WGSAF(llist):
                fout.write(l)
                Count_Pass +=1

            if Count_All % 10000 == 0:
                sys.stderr.write( "Read %d Variants\r" % Count_All)
        fout.close()
        sys.stderr.write("Finish Reading All %d Variants. %d Variants Pass the Filter\n" % (Count_All, Count_Pass))

    def F_WGSAF(self, llist):
        if len(llist[4].split(",")) > 1:
            return False
        if llist[6] != "PASS":
            return False
        INFO = llist[7]
        infolist = INFO.split(';')
        infodict = {}
        for kv in infolist:
            kv = kv.split('=')
            if len(kv) == 2:
                k, v = kv
                if k in infodict:
                    infodict[k] = infodict[k] + ',' + v
                else:
                    infodict[k] = v
        AF = infodict['gnomAD_genome_ALL'].split(',')[-1]
        if AF == '.':
            return False
        if float(AF) > self.freq:
            return True
        else:
            return False


def main():
    args = GetOptions()
    ins = WGS_Common(args) 
    ins.run()
    #Filter(VCFin, VCFout, lower, upper)


if __name__ == '__main__':
    main()
