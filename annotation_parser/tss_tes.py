'''
returns recombination rate as a function of distance from
TSS + TES, where TSS is the start of UTR5 and 
TES is the start of UTR3 for a given gene in a GFF

usage:
python3.5 tss_tes.py --gff [.gff] --table [.txt.gz] --distance [int] --chromosome [matching gff] > out.txt
'''

import argparse
import sys
from tqdm import tqdm

try:
    import antr
except ImportError:
    sys.path.append('/scratch/research/projects/chlamydomonas/\
        genomewide_recombination/analysis/fiftyshadesofgreen/annotation_parser/')

    # args
parser = argparse.ArgumentParser(description = 'Calculate \
        rho as a fxn of dist from TSS/TES',
        usage = 'tss_tes.py [options]')

parser.add_argument('-t', '--table', required = True,
    type = str, help = 'Annotation table (.txt.gz)')
parser.add_argument('-g', '--gff', required = True,
    type = str, help = 'GFF file (.gff/gff3)')
parser.add_argument('-d', '--distance', required = True,
    type = int, help = 'Distance from TES/TSS to consider')
parser.add_argument('-c', '--chromosome', required = True,
    type = str, help = 'Chromosome to do calculations for')

args = parser.parse_args()

table = args.table
gff = args.gff
distance = int(args.distance)
chromosome = args.chromosome

print('type distance rho') # colnames

utr_lines = []
with open(gff, 'r') as f:
    for line in tqdm(f):
        if chromosome in line:
            p = antr.Reader(table)

            if '_prime_UTR' in line:

                sp = line.split('\t')[0:7] # type start end . strand
                chrom, origin, utr, start, end, junk, strand = sp
                start = int(start)
                end = int(end)

                if 'five' in utr: # 5' UTR
                    if strand == '+':
                        region_outside = (start - distance, start)
                        region_inside = (start, start + distance)

                        for record in p.fetch(chrom, region_outside[0], region_outside[1]): # 'left side' of TSS
                            rho = record.ld_rho
                            dist = record.pos - start
                            assert dist <= 0
                            out = ' '.join(['TSS', dist, rho])
                            print(out)

                        for record in p.fetch(chrom, region_inside[0], region_inside[1]): # 'right side' of TSS
                            rho = record.ld_rho
                            dist = record.pos - start
                            assert dist >= 0
                            out = ' '.join(['TSS', dist, rho])
                            print(out)

                    if strand == '-':
                        region_outside = (end, end + distance)
                        region_inside = (end - distance, end)

                        for record in p.fetch(chrom, region_outside[0], region_outside[1]):
                            rho = record.ld_rho
                            dist = end - record.pos
                            assert dist <= 0
                            out = ' '.join(['TSS', dist, rho])
                            print(out)

                        for record in p.fetch(chrom, region_inside[0], region_inside[1]):
                            rho = record.ld_rho
                            dist = end - record.pos
                            assert dist >= 0
                            out = ' '.join(['TSS', dist, rho])

                elif 'three' in utr: # 3' UTR - outside gene is POSITIVE and inside gene is NEGATIVE
                    if strand == '+':
                        region_outside = (end, end + distance)
                        region_inside = (end - distance, end)

                        for record in p.fetch(chrom, region_outside[0], region_outside[1]):
                            rho = record.ld_rho
                            dist = record.pos - end
                            assert dist >= 0
                            out = ' '.join(['TES', dist, rho])
                            print(out)

                        for record in p.fetch(chrom, region_inside[0], region_inside[1]):
                            rho = record.ld_rho
                            dist = record.pos - end
                            assert dist <= 0
                            out = ' '.join(['TES', dist, rho])
                            print(out)

                    if strand == '-':
                        region_outside = (start - distance, start)
                        region_inside = (start, start + distance)

                        for record in p.fetch(chrom, region_outside[0], region_outside[1]):
                            rho = record.ld_rho
                            dist = start - record.pos
                            assert dist >= 0
                            out = ' '.join(['TES', dist, rho])
                            print(out)

                        for record in p.fetch(chrom, region_inside[0], region_inside[1]):
                            rho = record.ld_rho
                            dist = start - record.pos
                            assert dist <= 0
                            out = ' '.join(['TES', dist, rho])
