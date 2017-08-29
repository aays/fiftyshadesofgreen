'''
ldcalc.py - functions to calculate linkage disequilibrium stats from vcf files

Uses LD functions defined in popgen.py.

AH - 06/2017
'''

import vcf
import random
from tqdm import tqdm
from popgen import *

def snppuller(vcf_file, chrom = None, pos = None):
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
    if chrom is not None and pos is not None:
        pos = pos.split('-')
        for record in vcfin.fetch(chrom = chrom, start = pos[0], end = pos[1]):
            if hardsnpcheck(record) == True and issingleton(record) == False and isinvariant(record) == False:
                yield record
            else:
                pass
    elif chrom is not None and pos is None:
        for record in vcfin.fetch(chrom = chrom):
            if hardsnpcheck(record) == True and issingleton(record) == False and isinvariant(record) == False:
                yield record
            else:
                pass
    else:
        for record in vcfin:
            if hardsnpcheck(record) == True and issingleton(record) == False and isinvariant(record) == False:
                yield record
            else:
                pass
        
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
                print('chrom1', 'pos1', 'chrom2', 'pos2', 'd', 'hapcount')
            elif 'dprime' in stat:
                print('chrom1', 'pos1', 'chrom2', 'pos2', 'dprime', 'hapcount')
            elif 'r2' in stat:
                print('chrom1', 'pos1', 'chrom2', 'pos2', 'r2', 'hapcount')
        elif len(stat) == 2:
            if 'd' in stat and 'dprime' in stat:
                print('chrom1', 'pos1', 'chrom2', 'pos2', 'd', 'dprime', 'hapcount')
            elif 'd' in stat and 'r2' in stat:
                print('chrom1', 'pos1', 'chrom2', 'pos2', 'd', 'r2', 'hapcount')
            elif 'dprime' in stat and 'r2' in stat:
                print('chrom1', 'pos1', 'chrom2', 'pos2', 'dprime', 'r2', 'hapcount')
        elif len(stat) == 3:
            print('chrom1', 'pos1', 'chrom2', 'pos2', 'd', 'dprime', 'r2', 'hapcount')
        
        
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
