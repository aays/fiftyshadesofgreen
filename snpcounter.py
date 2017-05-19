#!/usr/bin/python3.5

'''
Given a zipped vcf file and a window size, calculates the number of SNPs in each given window.
Chromosome lengths are currently hardcoded for the Chlamydomonas reinhardtii genome.

usage:
python3.5 snpsampler.py [vcf (.gz)] [windowsize] > [outfile]

reference - notebook 8.6

AH - 04/2017
'''

import sys
import vcf

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
    '''Fetches chrom name from vcf. Helper function for variantgetter.'''
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

def variantgetter(vcfinput, windowmin, windowmax):
    '''Returns number of SNPs within a given window'''
    vcfin = vcf.Reader(filename = vcfinput, compressed = True) 
    chromname = chromgetter(vcfinput)
    counter = 0
    for record in vcfin.fetch(chrom = chromname, start = windowmin, end = windowmax - 1):
        if record.POS >= windowmin and record.POS <= windowmax and record.is_snp == True:
            counter = counter + 1
        else:
            pass
    return counter


# analysis

infile = sys.argv[1]
windowsize = int(sys.argv[2])

chromname = chromgetter(infile)
ranges = windowranger(windowsize, lengths[chromname]) # get ranges

print('start end snpcount') # initialize headers

for start, end in zip(ranges[0], ranges[1]):
    print(start, end, variantgetter(infile, start, end)) # add data
