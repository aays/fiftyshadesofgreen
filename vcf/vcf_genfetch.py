#!usr/bin/env python3.5

'''
Given a tabix-indexed vcf, start/end coordinates, and an outfile name, writes a new vcf 
containing that snippet.

usage:
python3.5 vcf_genfetch.py <input vcf> <chromosome> <start> <end> <outname>

e.g.
python3.5 vcf_genfetch.py genome.vcf chromosome_6 280000 290000 snippet.vcf

will write records from position 280000 to 290000 to a new file called snippet.vcf.

I might eventually update this with a failsafe that uses subprocess to create a tabix file
should there not be one.

AH - 05/2017

'''

import vcf
import sys

vcfin = str(sys.argv[1])
chrom = str(sys.argv[2])
start = int(sys.argv[3])
end = int(sys.argv[4])
outfile = str(sys.argv[5]) 

file = vcf.Reader(filename = vcfin, compressed = True)
snippet = file.fetch(chrom = chrom, start = start, end = end) 

with open(outfile, 'w') as out:
    writer = vcf.Writer(out, snippet)
    for record in snippet:
        writer.write_record(record) 

print('VCF written to', outfile)         
