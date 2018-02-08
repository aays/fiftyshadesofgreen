'''
pick out variants at random from the VCF and calculate LD between them

can parallelize this by suppressing col header printing, running in parallel,
and then concatenating the outfiles

usage:
python3.5 ld_decay.py --vcf [vcf.gz] --count [number of comparisons] --filter [optional distance limit] > [outfile]
'''

import vcf
import random
import argparse
from popgen import *
from ldcalc import *

parser = argparse.ArgumentParser(description = 'Calculate LD stats between randomly selected sites in a VCF file.',
                                usage = 'ld_decay.py [options] > [outfile]')

parser.add_argument('-v', '--vcf', required = True, type = str, help = 'Input VCF(.gz)')
parser.add_argument('-c', '--count', required = True, type = int, help = 'Number of pairwise draws to consider in calculations')
parser.add_argument('-f', '--filter', required = False, type = int, help = 'Maximum distance between SNPs [optional]')
parser.add_argument('-t', '--title', required = False, action = 'store_true', help = 'Print column headers? [optional]')

args = parser.parse_args()
filename = str(args.vcf)
total_count = int(args.count)
filt = args.filter
incl_title = args.title

vcfin = vcf.Reader(filename = filename, compressed = True)
chrom = next(vcfin).CHROM
chrom_length = vcfin.contigs[chrom][1] # get chromosome length

def get_random_record(filename, chrom_length, offset, limit_center = None, distance_limit = None):
    '''
    filename = vcf file
    chrom_length = length of chromosome
    offset = when creating generator objects, how far ahead of drawn position to look?
    limit_center = [default None] an integer value, representing a site position around which to look distance_limit bp
    on either side
    distance_limit = [default None] an integer value, describing how far around the limit center to look

    returns a random record in the vcf given these parameters
    '''
    found = False
    while not found:
        try:
            if not distance_limit:
                pos = random.randint(1, chrom_length)
                record = next(snppuller(filename, chrom, pos, pos + offset))
            elif distance_limit:
                in_range = False
                while not in_range:
                    pos = random.randint(limit_center - distance_limit, limit_center + distance_limit)
                    if pos > chrom_length or pos < 0: # out of range - fetch another record
                        continue
                    else:
                        record = next(snppuller(filename, chrom, pos, pos + offset))
                        in_range = True
            found = True
        except StopIteration:
            continue
    return record

if incl_title:
    print('chrom1 pos1 chrom2 pos2 r2')
    
for i in tqdm(range(total_count)):
    if not filt:
        record1 = get_random_record(filename, chrom_length, 100)
        record2 = get_random_record(filename, chrom_length, 100)
    elif filt:
        record1 = get_random_record(filename, chrom_length, 100)
        record2 = get_random_record(filename, chrom_length, 100, \
            limit_center = record1.POS, distance_limit = int(filt))
    r2 = r2calc(record1, record2, snpcheck = False)
    print(record1.CHROM, record1.POS, record2.CHROM, record2.POS, r2)
