'''
allpairscalc.py - calculates LD stats between _all_ SNPs in a given vcf

usage: python allpairscalc.py -v [vcf (.gz)] -l [ld stats] > [outfile]

AH - 09/2017
'''

import vcf
import argparse
from tqdm import tqdm
from popgen import * # snppuller, header
from ldcalc import *

parser = argparse.ArgumentParser(description = 'Calculate LD stats between all SNP pairs in a VCF file.', 
                                usage = 'allpairscalc.py [options]')

parser.add_argument('-v', '--vcfinput', required = True,
                   type = str, help = 'Input VCF')
parser.add_argument('-l', '--ldstats', required = True,
                   type = str, help = 'LD stats to calculate, separated by forward slashes (i.e. d/dprime)')
parser.add_argument('-p', '--haps', required = False,
                   action = 'store_true', help = 'Whether to output the number of observed haplotypes for each comparison. Optional.')


args = parser.parse_args()
vcfin = args.vcfinput
ldstats = args.ldstats
haps = args.haps

def allpairscalc(vcf_file, stat, haps = False):
    '''In a single VCF, calculates linkage stats between all possible pairs, regardless of region.
    The stat parameter can take any of 'd', 'dprime', or 'r2' as input. 
    Multiple parameters can be provided if separated by a forward slash (ie 'd/dprime' or 'r2/d').
    Output is printed in a space separated format.
    
    Setting haps to True will also print out how many of the four possible haplotypes
    were actually present in the final comparison.
    '''
    
    def metadata(record1, record2):
        out = record1.CHROM + ' ' + str(record1.POS) + ' ' + record2.CHROM + ' ' + str(record2.POS)
        return out
      
    def ldgetter(record1, record2):
        if haps == False: # proceed w/o haps
            if len(stat) == 1:
                if 'd' in stat:
                    print(metadata(record1, record2), dcalc(record1, record2))
                elif 'dprime' in stat:
                    print(metadata(record1, record2), dprimecalc(record1, record2))
                elif 'r2' in stat:
                    print(metadata(record1, record2), r2calc(record1, record2))
            elif len(stat) == 2:
                if 'd' in stat and 'dprime' in stat:
                    print(metadata(record1, record2), dcalc(record1, record2), dprimecalc(record1, record2))
                elif 'd' in stat and 'r2' in stat:
                    print(metadata(record1, record2), dcalc(record1, record2), r2calc(record1, record2))
                elif 'dprime' in stat and 'r2' in stat:
                    print(metadata(record1, record2), dcalc(record1, record2), r2calc(record1, record2))
            elif len(stat) == 3:
                print(metadata(record1, record2), dcalc(record1, record2), dprimecalc(record1, record2), r2calc(record1, record2))

        elif haps == True: # show haps/4 for each comparison
            observed_haps = list(freqsgetter(record1, record2)[2].values()) # get hap frequencies
            hapcount = 4 - observed_haps.count(0)
            if len(stat) == 1:
                if 'd' in stat:
                    print(metadata(record1, record2), dcalc(record1, record2), hapcount)
                elif 'dprime' in stat:
                    print(metadata(record1, record2), dprimecalc(record1, record2), hapcount)
                elif 'r2' in stat:
                    print(metadata(record1, record2), r2calc(record1, record2), hapcount)
            elif len(stat) == 2:
                if 'd' in stat and 'dprime' in stat:
                    print(metadata(record1, record2), dcalc(record1, record2), dprimecalc(record1, record2), hapcount)
                elif 'd' in stat and 'r2' in stat:
                    print(metadata(record1, record2), dcalc(record1, record2), r2calc(record1, record2), hapcount)
                elif 'dprime' in stat and 'r2' in stat:
                    print(metadata(record1, record2), dcalc(record1, record2), r2calc(record1, record2), hapcount)
            elif len(stat) == 3:
                print(metadata(record1, record2), dcalc(record1, record2), dprimecalc(record1, record2), r2calc(record1, record2), hapcount) 

    left = snppuller(vcf_file) # create ref vcf record generator
    stat = stat.split('/') # get stat
    header(stat, haps) # print header
    
    for record1 in tqdm(left):
        right = snppuller(vcf_file)
        if len(record1.ALT) > 1:
            continue
        for record2 in right:
            if record1.CHROM == record2.CHROM and record1.POS == record2.POS:
                continue
            else:
                ldgetter(record1, record2)

allpairscalc(vcf_file = vcfin, stat = ldstats, haps = haps)
