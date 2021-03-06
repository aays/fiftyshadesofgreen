'''
ldcalc.py - functions to calculate linkage disequilibrium stats from vcf files

Uses LD functions defined in popgen.py.

AH - 06/2017
'''

import vcf
import random
import itertools
from tqdm import tqdm
from popgen import *

def snppuller(vcf_file, chrom = None, start = None, end = None, return_none = False):
    '''Returns a generator object for a specified VCF snippet that returns only
    SNPs, while filtering out singletons.
    '''
    
    vcfin = vcf.Reader(filename = vcf_file, compressed = True)
    
    # filters
    def hardsnpcheck(record): # ensure biallelic SNP
        if len(record.REF) == 1 and len(record.ALT) == 1 and len(record.ALT[0]) == 1 and len(record.alleles) == 2:
            return True
        elif len(record.REF) != 1 or len(record.ALT) != 1 or len(record.ALT[0]) != 1 or len(record.alleles) != 2:
            return False
        
    def issingleton(record): # ensure not singleton
        if type(record.INFO['AN']) == list:
            count = record.INFO['AN'][0] - record.INFO['AC'][0]
        else:
            count = record.INFO['AN'] - record.INFO['AC'][0]
        if count == 1:
            return True
        if record.INFO['AC'][0] == 1:
            return True
        else:
            return False
        
    def isinvariant(record):
        if record.INFO['AC'][0] == 0:
            return True
        if record.INFO['AF'][0] == 1.0:
            return True
        else:
            return False
        
    # fetch
    acquired = False
    if chrom and start and end:
        for record in vcfin.fetch(chrom = chrom, start = start, end = end):
            if hardsnpcheck(record) and not issingleton(record) and not isinvariant(record):
                acquired = True
                yield record
            else:
                pass
    elif chrom and not start or not end:
        for record in vcfin.fetch(chrom = chrom):
            if hardsnpcheck(record) and not issingleton(record) and not isinvariant(record):
                acquired = True
                yield record
            else:
                pass
    else:
        for record in vcfin:
            if hardsnpcheck(record) and not issingleton(record) and not isinvariant(record):
                acquired = True
                yield record
            else:
                pass
    if return_none and not acquired:
        return None
        
def header(stat, haps = False):
    '''Helper function that determines output headers in singlevcfcalc.
    '''
    
    if haps == False:
        if len(stat) == 1:
            if 'd' in stat:
                print('chrom1', 'pos1', 'chrom2', 'pos2', 'd')
            elif 'dprime' in stat:
                print('chrom1', 'pos1', 'chrom2', 'pos2', 'dprime')
            elif 'r2' in stat:
                print('chrom1', 'pos1', 'chrom2', 'pos2', 'r2')
        elif len(stat) == 2:
            if 'd' in stat and 'dprime' in stat:
                print('chrom1', 'pos1', 'chrom2', 'pos2', 'd', 'dprime')
            elif 'd' in stat and 'r2' in stat:
                print('chrom1', 'pos1', 'chrom2', 'pos2', 'd', 'r2')        
            elif 'dprime' in stat and 'r2' in stat:
                print('chrom1', 'pos1', 'chrom2', 'pos2', 'dprime', 'r2')
        elif len(stat) == 3:
            print('chrom1', 'pos1', 'chrom2', 'pos2', 'd', 'dprime', 'r2')
    elif haps == True:
        if len(stat) == 1:
            if 'd' in stat:
                print('chrom1', 'pos1', 'chrom2', 'pos2', 'd', 'hapcount', 'haplist')
            elif 'dprime' in stat:
                print('chrom1', 'pos1', 'chrom2', 'pos2', 'dprime', 'hapcount', 'haplist')
            elif 'r2' in stat:
                print('chrom1', 'pos1', 'chrom2', 'pos2', 'r2', 'hapcount', 'haplist')
        elif len(stat) == 2:
            if 'd' in stat and 'dprime' in stat:
                print('chrom1', 'pos1', 'chrom2', 'pos2', 'd', 'dprime', 'hapcount', 'haplist')
            elif 'd' in stat and 'r2' in stat:
                print('chrom1', 'pos1', 'chrom2', 'pos2', 'd', 'r2', 'hapcount', 'haplist')
            elif 'dprime' in stat and 'r2' in stat:
                print('chrom1', 'pos1', 'chrom2', 'pos2', 'dprime', 'r2', 'hapcount', 'haplist')
        elif len(stat) == 3:
            print('chrom1', 'pos1', 'chrom2', 'pos2', 'd', 'dprime', 'r2', 'hapcount', 'haplist')
        

def gethaps(record1, record2, missing = True):
    '''Returns haps observed between two records as a list. Helper function for
    singlevcfcalc. Based on doublegtcounts() from popgen.
    '''
    # check strains b/w compared records are identical
    strainlist = [record1.samples[i].sample for i in range(len(record1.samples))] 
    assert strainlist == [record2.samples[i].sample for i in range(len(record2.samples))]
    # parse through VCF calls
    out = []
    for strain in strainlist:
        gt1 = record1.genotype(strain)['GT']
        gt2 = record2.genotype(strain)['GT']
        if gt1 == '.' and gt2 == '.' and missing == True:
            out.append('--')
            continue
        elif gt1 == '.' and missing == True:
            if gt2 == '0':
                out.append('-' + record2.REF)
            elif gt2 == '1':
                out.append('-' + str(record2.ALT[0]))
            else:
                continue
        elif gt1 == '0':
            if gt2 == '.' and missing == True:
                out.append(record1.REF + '-')
            elif gt2 == '0':
                out.append(record1.REF + str(record2.REF))
            elif gt2 == '1':
                out.append(record1.REF + str(record2.ALT[0]))
        elif gt1 == '1':
            if gt2 == '.' and missing == True:
                out.append(str(record1.ALT[0]) + '-')
            elif gt2 == '0':
                out.append(str(record1.ALT[0]) + record2.REF)
            elif gt2 == '1':
                out.append(str(record1.ALT[0]) + str(record2.ALT[0]))
    out = ','.join([str(hap) for hap in out]) # print nicer
    return out
        
def singlevcfcalc(vcf_file, ref, target, stat, filter = None, windowsize = None, haps = False):
    '''
    In a single VCF, calculates linkage stats between two entire regions.
    The stat parameter can take any of 'd', 'dprime', or 'r2' as input. 
    Multiple parameters can be provided if separated by a forward slash (ie 'd/dprime' or 'r2/d').
    Output is printed in a space separated format.
    
    usage - singlevcfcalc('myfile.vcf.gz', 'chromosome_6', 'chromosome_7', 'd/r2')
    will calculate both d and r2 between sites on chr6 and chr7.
    
    For intrachromosomal stats, simply input the same region twice:
    singlevcfcalc('myfile.vcf.gz', 'chromosome_6', 'chromosome_6', 'd/dprime')
    
    A filter can also be provided - setting filter = 0.8 will drop records in both
    *just target* roughly 20% of the time.
    
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
            haplist = gethaps(record1, record2)
            if len(stat) == 1:
                if 'd' in stat:
                    print(metadata(record1, record2), dcalc(record1, record2), hapcount, haplist)
                elif 'dprime' in stat:
                    print(metadata(record1, record2), dprimecalc(record1, record2), hapcount, haplist)
                elif 'r2' in stat:
                    print(metadata(record1, record2), r2calc(record1, record2), hapcount, haplist)
            elif len(stat) == 2:
                if 'd' in stat and 'dprime' in stat:
                    print(metadata(record1, record2), dcalc(record1, record2), dprimecalc(record1, record2), hapcount, haplist)
                elif 'd' in stat and 'r2' in stat:
                    print(metadata(record1, record2), dcalc(record1, record2), r2calc(record1, record2), hapcount, haplist)
                elif 'dprime' in stat and 'r2' in stat:
                    print(metadata(record1, record2), dcalc(record1, record2), r2calc(record1, record2), hapcount, haplist)
            elif len(stat) == 3:
                print(metadata(record1, record2), dcalc(record1, record2), dprimecalc(record1, record2), r2calc(record1, record2), hapcount, haplist)
            
    reflocus = snppuller(vcf_file, chrom = ref) # create ref vcf record generator
    stat = stat.split('/') # get stat
    header(stat, haps) # print header
    
    if not filter:
        for record1 in tqdm(reflocus):
            targetlocus = snppuller(vcf_file, chrom = target)
            if len(record1.ALT) > 1:
                continue
            for record2 in targetlocus:
                if not windowsize:
                    ldgetter(record1, record2)
                elif windowsize: # if a windowsize is provided
                    if abs(record2.POS - record1.POS) <= windowsize:
                        ldgetter(record1, record2)
                    elif abs(record2.POS - record1.POS) > windowsize:
                        continue
                        
    elif filter:
        for record1 in tqdm(reflocus):
            targetlocus = snppuller(vcf_file, chrom = target)
            if len(record1.ALT) > 1:
                continue
            for record2 in targetlocus:
                if random.random() <= filter:
                    if not windowsize:
                        ldgetter(record1, record2)
                    elif windowsize:
                        if abs(record2.POS - record1.POS) <= windowsize:
                            ldgetter(record1, record2)
                        elif abs(record2.POS - record1.POS) > windowsize:
                            continue
                elif random.random() > filter:
                    continue

def sequentialvcfcalc(vcf_file, ref, target, stat, windowsize = None, haps = False):
    '''
    Calculates LD for 'sequential' pairs as described in Lewontin 1995.
    
    Consider three loci A, B, and C - singlevcfcalc would calculate LD for the AB, BC, and AC pairs.
    However, if AB and BC are both in LD, it follows that AC are in LD. This means certain pairs are
    technically nonindependent, which would bias tests for 'LD hotspots' and inflate the number of
    pairwise comparisons in apparent complete linkage.
    
    sequentialvcfcalc computes LD between two regions in a 'forward' and then 'reverse' fashion.
    If loci A, B, and C are on region 1, while region 2 features a, b, and c, the function will first
    calculate LD for Aa, Bb, and Cc in the 'forward' direction, and then compute aB, bC, and so on in the
    'reverse' direction. The net effect is acquiring LD values for Aa, aB, Bb, bC, and so forth.
    This preserves independence of pairwise tests.
    
    Removed filter option from singlevcfcalc - was unused. For filtering, filter input vcf using 
    vcf_subset prior to calculation.
    
    Caution: sequentialvcfcalc behaves erratically when presented with two regions
    of different lengths, either repeating earlier SNPs for new comparisons or ending its iterations early. 
    As such, it is best used for intra-region calculations. There are currently no plans to introduce 
    mathematically robust inter-region compatibility.
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
            
    stat = stat.split('/') # get stat
    header(stat, haps) # print header
    # create ref vcf record list. this cannot be a generator bc itertools.cycle()
    reflocus = [record for record in reclist(vcf_file, chrom = ref, snpsonly = True)]
    targetlocus = [record for record in reclist(vcf_file, chrom = target, snpsonly = True)]

    # forward
    for record1, record2 in tqdm(zip(itertools.cycle(reflocus), targetlocus)):
        if not windowsize:
            ldgetter(record1, record2)
        elif windowsize:
            if abs(record2.POS - record1.POS) <= windowsize:
                ldgetter(record1, record2)
            elif abs(record2.POS - record1.POS) > windowsize:
                continue
        continue

    # reverse
    # load in records again
    reflocus_rev = [record for record in reclist(vcf_file, chrom = ref, snpsonly = True)][1:] # 'waste' first record to make offset
    targetlocus_rev = [record for record in reclist(vcf_file, chrom = target, snpsonly = True)]

    for record2, record1 in tqdm(zip(itertools.cycle(targetlocus_rev), reflocus_rev)): # keep ordering of rec1, rec2 consistent with previous for loop
        if len(record1.ALT) > 1:
            continue
        if not windowsize:
            ldgetter(record1, record2)
        elif windowsize: # if a windowsize is provided
            if abs(record2.POS - record1.POS) <= windowsize:
                ldgetter(record1, record2)
            elif abs(record2.POS - record1.POS) > windowsize:
                continue

  
