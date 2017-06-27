#!/bin/env python3.5

'''
singlevcfcalc.py - calculate ld stats between two regions within a single vcf file.

usage: python singlevcfcalc.py -v [vcf (.gz)] -r [region 1] [region 2] -l [ld stats] > [outfile]

LD stats can be any combination of d, dprime, or r2 -
separate by a slash (i.e. r2/dprime) for more than one.

if looking to do an intra-region comparison, simply enter the same region name
for both region 1 and region 2 (i.e. -v myvcf.gz -r chromosome_2 chromosome_2 -l r2/dprime > chr2.ld)

AH - 06/2017
'''

import vcf
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

args = parser.parse_args()
vcfin = args.vcfinput
regions = args.regions
ldstats = args.ldstats

singlevcfcalc(vcfin, ref = regions[0], target = regions[1], stat = ldstats)
