#!/usr/bin/python3.5

'''
usage:
python3.5 snpsample2vcf.py [vcf (.gz)] [windowsize] [chrom] [outfile]

e.g.
python3.5 snpsample2vcf.py chromosome10.vcf.gz 100000 chromosome_10 out
will create a vcf called outsampled.vcf.

Number of SNPs picked out from each window are hardcoded in the snpcounts dictionary,
based on calculations from snpcounter.py. These were calculated such that an equivalent
amount of SNPs are obtained from each chromosome so as not to bias downstream inter-chromosome 
r2 comparisons, and correspond to a window size of 100kb.

In order to reuse this script for a different window size/genome, one would have to
run snpcounter.py for the desired window size and with the chromosome lengths modified as
necessary, after inputting the respective SNP counts (as well as chromosome lengths) into this script.

reference - notebook 8.6c

AH - 04/2017
'''

import sys
import vcf
import random

random.seed(42) 

snpcounts = {'chromosome_1': 25,
'chromosome_10': 30,
'chromosome_11': 51,
'chromosome_12': 20,
'chromosome_13': 38,
'chromosome_14': 48,
'chromosome_15': 100,
'chromosome_16': 26,
'chromosome_17': 28,
'chromosome_2': 22,
'chromosome_3': 22,
'chromosome_4': 49,
'chromosome_5': 56,
'chromosome_6': 22,
'chromosome_7': 31,
'chromosome_8': 39,
'chromosome_9': 25} 

lengths = {'chromosome_1': 8033585,
'chromosome_2': 9223677,
'chromosome_3': 9219486,
'chromosome_4': 4091191,
'chromosome_5': 3500558,
'chromosome_6': 9023763,
'chromosome_7': 6421821,
'chromosome_8': 5033832,
'chromosome_9': 7956127,
'chromosome_10': 6576019,
'chromosome_11': 3826814,
'chromosome_12': 9730733,
'chromosome_13': 5206065,
'chromosome_14': 4157777,
'chromosome_15': 1922860,
'chromosome_16': 7783580,
'chromosome_17': 7188315} 

def chromgetter(vcfinput):
    '''Fetches chrom name from vcf. Helper function for randvariantgetter.'''
    thisvcf = vcf.Reader(filename = vcfinput, compressed = True)
    chromname = next(thisvcf).CHROM
    return chromname

def windowranger(windowsize, totalsize):
    '''Returns a list of mins and maxes given a window size and
    max length. These lists will be in a tuple that needs to be
    unpacked.'''
    windowmins = list(range(0, totalsize, windowsize))
    windowmaxes = windowmins[1:len(windowmins)]
    windowmaxes.append(totalsize)
    windowmins = [num + 1 for num in windowmins]
    return windowmins, windowmaxes

def randvariantgetter(vcfinput, chrom, windowmin, windowmax, snplist):
    '''Returns randomly selected SNP positions within current window.
    Number of SNPs picked out for given chrom is determined from snpcounts 
    dict. Feed empty list if first iteration.'''
    vcfin = vcf.Reader(filename = vcfinput, compressed = True) 
    snppositions = [] # reset list
    snippet = vcfin.fetch(chrom = chrom, start = windowmin, end = windowmax - 1)
    try:
        snppositions = [record.POS for record in snippet]
    except IndexError:
        pass
    try:
        snpsample = sorted(random.sample(snppositions, snpcounts[chrom]))
    except ValueError:
        snpsample = snppositions
    snpsample = snplist + snpsample
    return snpsample

def variantwriter(vcfinput, inchrom, outfile, snplist):
    '''Given a list of a SNPs and an input, writes a vcf file
    that only contains those records.'''
    vcfin = vcf.Reader(filename = vcfinput, compressed = True) 
    snippet = vcfin.fetch(chrom = inchrom)
    with open(outfile, 'w') as outvcf:
        writer = vcf.Writer(outvcf, snippet)
        try:
            for record in snippet:
                if record.POS in snplist:
                    writer.write_record(record)
                else:
                    pass
        except IndexError:
            pass 

# analysis
vcfin = str(sys.argv[1])
windowsize = int(sys.argv[2])
chrom = str(sys.argv[3])
outfilename = str(sys.argv[4])

ranges = windowranger(windowsize, lengths[chrom])

snplist = []

for start, end in zip(ranges[0], ranges[1]):
    snplist = randvariantgetter(vcfin, chrom, start, end, snplist)
    
outname = outfilename + 'sampled.vcf'
    
variantwriter(vcfin, chrom, outname, snplist)    

