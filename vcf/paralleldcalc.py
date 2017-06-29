'''
parallelldcalc.py - an attempt to parallelize ldcalc's operations

Uses LD functions defined in popgen.py and ldcalc.py

probably doesn't work yet so be wary

AH - 06/2017
'''

import vcf
import itertools
import sys
from multiprocessing import Pool
from popgen import *
from ldcalc import *

vcfin = sys.argv[1]
region1 = sys.argv[2]
region2 = sys.argv[3]
stats = sys.argv[4]
processes = sys.argv[5]

def metadata(record1, record2):
        out = record1.CHROM + ' ' + str(record1.POS) + ' ' + record2.CHROM + ' ' + str(record2.POS)
        return out

def ldgetter(record1, record2, stat):
    '''Helper function for parallelvcfcalc.
    '''
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

def parallelvcfcalc(vcf_file, ref, target, stat, num_processes = 1):
    '''
    In a single VCF, calculates linkage stats between two entire regions.
    The stat parameter can take any of 'd', 'dprime', or 'r2' as input. 
    Multiple parameters can be provided if separated by a forward slash (ie 'd/dprime' or 'r2/d').
    Output is printed in a space separated format.
    
    usage - parallelvcfcalc('myfile.vcf.gz', 'chromosome_6', 'chromosome_7', 'd/r2', 4)
    will calculate both d and r2 between sites on chr6 and chr7 using 4 parallel processes.
    
    For intrachromosomal stats, simply input the same region twice:
    parallelvcfcalc('myfile.vcf.gz', 'chromosome_6', 'chromosome_6', 'd/dprime', 4)
    
    num_processes will set how many processes Python will use in completing
    the operation.
    '''
    reflocus = snppuller(vcf_file, chrom = ref)
    stat = stat.split('/')
    header(stat) # print header
    # parallel process begins here
    for record1 in reflocus:
        targetlocus = snppuller(vcf_file, chrom = target)
        chunk = [i for i in itertools.islice(targetlocus, num_processes)]
        while True:
            pool = Pool(processes = num_processes)
            pool.starmap(ldgetter, zip(itertools.repeat(record1), chunk, itertools.repeat(stat)))
            pool.close()
            pool.join()
            chunk = [i for i in itertools.islice(targetlocus, num_processes)] # get next chunk
            if len(chunk) == 0:
                break

parallelvcfcalc(vcfin, region1, region2, stat = stats, num_processes = processes)



