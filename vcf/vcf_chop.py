#!usr/bin/env python3.5

'''
Given a tabix-indexed vcf, start/end coordinates, and an outfile name, writes a new vcf 
containing that snippet.

usage:
python3.5 vcf_chop.py -v <input vcf> -c <chromosome> -p <start-end> -o <outname>

e.g.
python3.5 vcf_chop.py -v genome.vcf -c chromosome_6 -p 280000-290000 -o snippet.vcf

will write records from position 280000 to 290000 to a new file called snippet.vcf.

AH - 05/2017

'''

import vcf
import argparse

parser = argparse.ArgumentParser(description = 'Subset a vcf file.', 
                                usage = 'vcf_chop.py [options]')

parser.add_argument('-v', '--vcfinput', required = True,
                   type = str, help = 'Input VCF')
parser.add_argument('-c', '--chrom', required = True,
                   type = str, help = 'Chromosome name (as it appears in the vcf)')
parser.add_argument('-p', '--positions', required = False,
                   type = str, help = 'Positions, in the format "start-end" (ie 100-200)')
parser.add_argument('-o', '--out', required = True,
                   type = str, help = 'Name of desired output file')

args = parser.parse_args()
vcfin = args.vcfinput
chrom = args.chrom
pos = args.positions
outfile = args.out
    
file = vcf.Reader(filename = vcfin, compressed = True)

if pos:
    pos = pos.split('-')
    start = int(pos[0])
    end = int(pos[1])
    snippet = file.fetch(chrom = chrom, start = start, end = end) 
else:
    snippet = file.fetch(chrom = chrom) 

with open(outfile, 'w') as out:
    writer = vcf.Writer(out, snippet)
    for record in snippet:
        writer.write_record(record) 

print('VCF written to', outfile)         
