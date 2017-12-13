#!/home/local/users/jw/anaconda2/bin/python
# Author: jywang	explorerwjy@gmail.com

#=========================================================================
# CalculateAC.py
#=========================================================================

import argparse
import gzip


def GetFil(Fname):
    if Fname.endswith(".gz"):
        return gzip.open(Fname, 'rb')
    else:
        return open(Fname, 'rb')


def ParseINFO(INFO_string):
    INFO_dict = {}
    for item in INFO_string.split(";"):
        kv = item.split("=")
        if len(kv) == 2:
            k, v = kv
            if k not in INFO_dict:
                INFO_dict[k] = v.split(",")
            else:
                INFO_dict[k].extend(v.split(","))
    return INFO_dict


def GetAC(Alts, Samples):
    Alleles = {}
    for i,alt in enumerate(Alts):
        Alleles[i+1] = 0
    NotCalled = 0
    for sample in Samples:
        GTs = sample.split(":")[0]
        try:
            GTs = map(int, GTs.split("/"))
        except:
            NotCalled += 1
            continue
        for gt in GTs:
            if gt == 0:
                continue
            #if gt not in Alleles:
            #    Alleles[gt] = 1
            #else:
            #    Alleles[gt] += 1
            Alleles[gt] += 1
    res = []
    for gt, count in sorted(Alleles.items(), key=lambda x: x[0]):
        res.append(count)
    return map(str, res), NotCalled


def GetAF(ACs, NumSamples, NotCalled):
    Num = NumSamples - NotCalled
    res = []
    for AC in ACs:
        res.append(round(float(AC) / (Num * 2), 7))
    return map(str, res)


class CalculateAC:
    def __init__(self, args):
        self.InpFil = args.vcf
        if args.out != None:
            self.OutFil = args.out
        else:
            self.OutFil = self.InpFil.rstrip(
                ".gz").rstrip(".vcf") + "ModifyAC.vcf"
        self.Calculate()


    def Calculate(self):
        print "Calculate AC"
        fin = GetFil(self.InpFil)
        fout = open(self.OutFil, 'wb')
        for l in fin:
            if l.startswith("##"):
                fout.write(l)
                continue
            elif l.startswith("#"):
                fout.write(l)
                NumSamples = len(l.strip().split("\t")) - 9
                print NumSamples
                continue
            llist = l.strip().split("\t")
            INFO_string = llist[7]
            INFO_dict = ParseINFO(INFO_string)
            ACs = INFO_dict["AC"]
            AFs = INFO_dict["MLEAF"]
            AN = int(INFO_dict["AN"][0])
            Samples = llist[9:]
            MyACs, NotCalled = GetAC(Samples)
            MyAFs = GetAF(MyACs, NumSamples, NotCalled)
            MyAN = (NumSamples - NotCalled) * 2
            INFO_list = INFO_string.split(";")
            FLAG = 0
            for i in xrange(len(INFO_list)):
                kv = INFO_list[i]
                try:
                    k, v = kv.split("=")
                    if k == "AC":
                        INFO_list[i] = "AC={}".format(",".join(MyACs))
                        FLAG += 1
                    if k == "AF":
                        INFO_list[i] = "AF={}".format(",".join(MyAFs))
                        FLAG += 1
                    if k == "AC":
                        INFO_list[i] = "AC={}".format(",".join(MyACs))
                        FLAG += 1
                except:
                    continue
                if FLAG == 3:
                    break
            NewINFO = ";".join(INFO_list)
            llist[7] = NewINFO
            NewLine = "\t".join(llist) + "\n"
            fout.write(NewLine)
        fin.close()
        fout.close()


def AFCompare(AFs, MyAFs):
    AFs = [round(float(x), 3) for x in AFs]
    MyAFs = [round(float(x), 3) for x in MyAFs]
    if AFs != MyAFs:
        return True
    else:
        return False


def GetOptions():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', type=str, help='Input VCF file')
    parser.add_argument('-o', '--out', type=str, help='Output VCF file')
    args = parser.parse_args()
    return args


def main():
    args = GetOptions()
    ins = CalculateAC(args)
    return


if __name__ == '__main__':
    main()
