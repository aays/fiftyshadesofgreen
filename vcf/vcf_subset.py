#!usr/bin/env python3.5

'''
Given a tabix-indexed vcf, start/end coordinates, and an outfile name, writes a new vcf 
containing that snippet. Can also be provided with a filter value in the form of a 
proportion (float, 0 < x < 1) or a number of sites to keep. Filtering uses the random
module and is less exact for smaller numbers of records. 

usage:
python3.5 vcf_subset.py -v <input vcf> -c <chromosome> -p <start-end> -o <outname>

e.g.
python3.5 vcf_subset.py -v genome.vcf -c chromosome_6 -p 280000-290000 -f 0.5 -o snippet.vcf

will write 0.5 of the records from position 280000 to 290000 to a new file called snippet.vcf.

AH - 05/2017

'''

import sys
import vcf
import random
import argparse
from tqdm import tqdm

parser = argparse.ArgumentParser(description = 'Subset and/or filter a vcf file.', 
                                usage = 'vcf_subset.py [options]')

parser.add_argument('-v', '--vcfinput', required = True,
                   type = str, help = 'Input VCF')
parser.add_argument('-c', '--chrom', required = False,
                   type = str, nargs = '+', help = 'Chromosome name(s) (as they appears in the vcf). Optional')
parser.add_argument('-p', '--positions', required = False,
                   type = str, help = 'Positions, in the format "start-end" (ie 100-200). Optional.')
parser.add_argument('-f', '--filter_fraction', required = False,
                   type = float, help = 'Fraction of sites to keep. Optional. Cannot be used with -n')
parser.add_argument('-n', '--filter_num', required = False,
                   type = int, help = 'Number of sites to keep. Optional. Cannot be used with -f. Only works if chrom(s) provided.')
parser.add_argument('-o', '--out', required = True,
                   type = str, help = 'Name of desired output file')

# define args
args = parser.parse_args()
vcfin = args.vcfinput
chrom = args.chrom
pos = args.positions
filt_frac = args.filter_fraction
filt_num = args.filter_num
outfile = args.out
    
# load in file    
file = vcf.Reader(filename = vcfin, compressed = True)

# assert both not present
if filt_frac and filt_num:
    print('Error - both filter fraction and filter-to-number args provided.')
    print('Please provide only one of the two args.')
    sys.exit(1)

# initial subset
if pos:
    pos = pos.split('-')
    start = int(pos[0])
    end = int(pos[1])
    snippet = file.fetch(chrom = chrom, start = start, end = end) 
elif chrom and not pos and len(chrom) == 1: # single chrom provided
    snippet = file.fetch(chrom = chrom[0])
elif chrom and not pos and len(chrom) > 1: # multiple chroms provided
    snippet = []
    for c in chrom:
        file = vcf.Reader(filename = vcfin, compressed = True) # since the generator resets each time
        snippet.append(file.fetch(chrom = c))
elif not chrom and not pos: # no chrom or pos provided - just filtering
    snippet = file

# filter and write
### given a fraction
if filt_frac: # given fraction between 0 and 1
    assert not filt_num
    with open(outfile, 'a') as out:
        if len(chrom) == 1:
            writer = vcf.Writer(out, snippet)
            for record in tqdm(snippet):
                if random.random() <= filt_frac: # filtering
                    writer.write_record(record)
        elif len(chrom) > 1:
            for region in snippet: # write and filter each region independently
                writer = vcf.Writer(out, region)
                for record in tqdm(region):
                    if random.random() <= filt_frac:
                        writer.write_record(record)
        elif not chrom:
            writer = vcf.Writer(out, snippet)
            for record in tqdm(snippet):
                if random.random() <= filt_frac:
                    writer.write_record(record)

### given a number                
elif filt_num and len(chrom) == 1: # given number of sites to keep - one region
    assert not filt_frac
    if pos:
        totalrecs = end - start # length of specified snippet
    elif not pos:
        totalrecs = snippet.contigs[chrom][1] # get length of entire chrom from metadata
    frac_from_num = filt_num/totalrecs # get filter fraction
    with open(outfile, 'w') as out:
        writer = vcf.Writer(out, snippet)
        for record in tqdm(snippet):
            if random.random() <= frac_from_num:
                writer.write_record(record)

elif filt_num and len(chrom) > 1: # multiple regions
    assert not filt_frac
    for region in snippet:
        totalrecs = snippet.contigs[region][1] # get length of current region
        frac_from_num  = filt_num/totalrecs # filter fraction
        with open(outfile, 'a') as out:
            writer = vcf.Writer(out, region)
            for record in tqdm(region):
                if random.random() <= frac_from_num:
                    writer.write_record(record)


# no filtering - just write                
elif not filt_frac and not filt_num:
    if len(chrom) == 1:
        with open(outfile, 'w') as out:
            writer = vcf.Writer(out, snippet)
            for record in tqdm(snippet):
                writer.write_record(record) 
    elif len(chrom) > 1:
        for region in snippet:
            with open(outfile, 'a') as out:
                writer = vcf.Writer(out, region)
                for record in tqdm(region):
                    writer.write_record(record)

print('VCF written to', outfile)         
