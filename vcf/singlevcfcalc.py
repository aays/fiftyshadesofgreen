'''
singlevcfcalc.py - calculate ld stats between two regions
within a single vcf

usage: python singlevcfcalc.py [vcf (.gz)] [region 1] [region 2] [ld stats] > [outfile]

ld stats can be any combination of d, dprime, or r2 -
separate by a slash (i.e. r2/dprime) for more than one

AH - 06/2017
'''

import sys
import vcf
from ldcalc import *
from popgen import *

vcfin = sys.argv[1]
region1 = sys.argv[2]
region2 = sys.argv[3]
ldstats = sys.argv[4]

singlevcfcalc(vcfin, ref = region1, target = region2, stat = ldstats)
