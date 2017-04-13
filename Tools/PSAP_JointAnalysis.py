#!/home/local/users/jw/anaconda2/bin/python
#Author: jywang	explorerwjy@gmail.com

#========================================================================================================
# PSAP_JointAnalysis.py
#========================================================================================================

import argparse
import csv


def GetOptions():
	parser = argparse.ArgumentParser()
	parser.add_argument('-d','--dir', type=str, help = 'dir contains psap result from PSAP_Report.py')
	parser.add_argument('-o','--out', type=str, help = 'OutPut Name')
	parser.add_argument('-p','--pvalue', type=str, help = 'PSAP Pvalue to cutoff')
	args = parser.parse_args()
	
	return args.

def main():

	return

if __name__=='__main__':
	main()
