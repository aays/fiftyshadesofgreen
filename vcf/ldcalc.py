'''
ldcalc.py - functions to calculate linkage disequilibrium stats from vcf files

Uses LD functions defined in popgen.py.

AH - 06/2017
'''

import vcf
from popgen import *

def snppuller(vcf_file, chrom = None, pos = None):
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
        
def header(stat):
    '''Helper function that determines output headers in singlevcfcalc.'''
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

def singlevcfcalc(vcf_file, ref, target, stat):
    '''
    In a single VCF, calculates linkage stats between two entire regions. (Positional options coming soon)
    The stat parameter can take any of 'd', 'dprime', or 'r2' as input. 
    Multiple parameters can be provided if separated by a forward slash (ie 'd/dprime' or 'r2/d').
    Output is printed in a space separated format.
    
    usage - singlevcfcalc('myfile.vcf.gz', 'chromosome_6', 'chromosome_7', 'd/r2')
    will calculate both d and r2 between sites on chr6 and chr7.
    
    For intrachromosomal stats, simply input the same region twice:
    singlevcfcalc('myfile.vcf.gz', 'chromosome_6', 'chromosome_6', 'd/dprime')
    '''
    def metadata(record1, record2):
        out = record1.CHROM + ' ' + str(record1.POS) + ' ' + record2.CHROM + ' ' + str(record2.POS)
        return out
    reflocus = snppuller(vcf_file, chrom = ref)
    stat = stat.split('/')
    header(stat)
    for record1 in reflocus:
        targetlocus = snppuller(vcf_file, chrom = target)
        if len(record1.ALT) > 1:
            continue
        for record2 in targetlocus:
            if len(record2.ALT) > 1:
                continue
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
