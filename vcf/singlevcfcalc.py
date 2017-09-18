#!/bin/env python3.5

'''
singlevcfcalc.py - calculate ld stats between two regions within a single vcf file.

usage: python singlevcfcalc.py -v [vcf (.gz)] -r [region 1] [region 2] -l [ld stats] -f [filter] -w [windowsize] > [outfile]

LD stats can be any combination of d, dprime, or r2 -
separate by a slash (i.e. r2/dprime) for more than one.

if looking to do an intra-region comparison, simply enter the same region name
for both region 1 and region 2 (i.e. -v myvcf.gz -r chromosome_2 chromosome_2 -l r2/dprime > chr2.ld)

AH - 06/2017
'''

import vcf
import random
import argparse
from ldcalc import *
from popgen import *

parser = argparse.ArgumentParser(description = 'Calculate LD stats between two regions in a VCF file.', 
                                usage = 'singlevcfcalc.py [options]')

parser.add_argument('-v', '--vcfinput', required = True,
                   type = str, help = 'Input VCF')
parser.add_argument('-r', '--regions', required = True,
                   type = str, nargs = '+', help = 'Two regions to compare (as they appear in the vcf)')
parser.add_argument('-l', '--ldstats', required = True,
                   type = str, help = 'LD stats to calculate, separated by forward slashes (i.e. d/dprime)')
parser.add_argument('-f', '--filter', required = False,
                   type = float, help = 'Proportion of records to perform calculations on. Optional.')
parser.add_argument('-w', '--windowsize', required = False,
                   type = int, help = 'Window size to calculate LD between. Only use for intra-region calculations. Optional.')
parser.add_argument('-p', '--haps', required = False,
                   action = 'store_true', help = 'Whether to output the number of observed haplotypes for each comparison. Optional.')
parser.add_argument('-q', '--sequential', required = False,
                   action = 'store_true', help = 'Calculate in sequential pairs instead of all pairs (Lewontin 1995). Optional.')


args = parser.parse_args()
vcfin = args.vcfinput
regions = args.regions
ldstats = args.ldstats
filt = args.filter
haps = args.haps
sequential = args.sequential
windowsize = args.windowsize

if sequential == True:
    sequentialvcfcalc(vcfin, ref = regions[0], target = regions[1], stat = ldstats, windowsize = windowsize, haps = haps)  
else:
    singlevcfcalc(vcfin, ref = regions[0], target = regions[1], stat = ldstats, filter = filt, windowsize = windowsize, haps = haps)
