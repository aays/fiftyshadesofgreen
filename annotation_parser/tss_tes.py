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
parser = argparse.ArgumentParser(description = 'Calculate rho as a fxn of dist from TSS/TES',
        usage = 'tss_tes.py [options]')

parser.add_argument('-t', '--table', required = True,
    type = str, help = 'Annotation table (.txt.gz)')
parser.add_argument('-g', '--gff', required = True,
    type = str, help = 'GFF file (.gff/gff3)')
parser.add_argument('-d', '--distance', required = True,
    type = int, help = 'Distance from TES/TSS to consider')
parser.add_argument('-c', '--chromosome', required = True,
    type = str, help = 'Chromosome to do calculations for')
parser.add_argument('-w', '--windowsize', required = False,
    type = int, help = 'Windowsize, if wanting to do a windowed comparison (optional)')

args = parser.parse_args()

table = args.table
gff = args.gff
distance = int(args.distance)
chromosome = args.chromosome
windowsize = args.windowsize

if windowsize:
    
    print('type windowleft windowright rho rho_total rho_count') # colnames

    def windowcalc(feature_type, strand, distance, windowsize, table, chromosome, region, start, end):
        '''Where region is a tuple of size 2, indicating the start and the end of the region
        i.e. windowcalc('TES', '+', 20, 'table.txt.gz', 'chromosome_2', (900, 1000), 16200, 17300)'''
        
        windowlist = list(range(region[0], region[1] + 1, windowsize))
        
        p = antr.Reader(table)
        
        for i in range(len(windowlist) - 1):
            windowleft, windowright = windowlist[i], windowlist[i + 1]
            rho_cumulative = 0.0
            count = 0
                          
            for record in p.fetch(chromosome, windowleft, windowright):
                rho_cumulative += record.ld_rho
                count += 1
            
            try:
                rho_out = rho_cumulative / count
            except ZeroDivisionError:
                assert count == 0
                rho_out = 0
                
            if feature_type == 'TSS' and strand == '+':
                windowleft_out = windowleft - start
                windowright_out = windowright - start
            elif feature_type == 'TSS' and strand == '-':
                windowleft_out = end - windowleft
                windowright_out = end - windowright
            elif feature_type == 'TES' and strand == '+':
                windowleft_out = windowleft - end
                windowright_out = windowright - end
            elif feature_type == 'TES' and strand == '-':
                windowleft_out = start - windowleft
                windowright_out = start - windowright


            windowout = ' '.join([str(feature_type), str(windowleft_out), str(windowright_out), \
                                  str(rho_out), str(rho_cumulative), str(count)])
            print(windowout)
                
    
else: # no windowsize
    print('type distance rho')

    def singlecalc(feature_type, strand, distance, table, chromosome, region, start, end):
        p = antr.Reader(table)

        for record in p.fetch(chromosome, region[0], region[1]):
            rho = record.ld_rho

            if feature_type == 'TSS' and strand == '+':
                dist = record.pos - start
            elif feature_type == 'TSS' and strand == '-':
                dist = end - record.pos
            elif feature_type == 'TES' and strand == '+':
                dist = record.pos - end
            elif feature_type == 'TES' and strand == '-':
                dist = start - record.pos

            out = ' '.join([feature_type, str(dist), str(rho)])
            print(out)


with open(gff, 'r') as f:
    for line in tqdm(f):
        if chromosome in line:
            p = antr.Reader(table)

            if '_prime_UTR' in line:

                sp = line.split('\t')[0:7]
                chrom, origin, utr, start, end, junk, strand = sp
                start = int(start)
                end = int(end)

                if 'five' in utr:
                    if strand == '+':
                        if start > distance:
                            region_outside = (start - distance, start)
                            region_inside = (start, start + distance)
                        else:
                            region_outside = (0, start)
                            region_inside = (start, start + distance)

                        if not windowsize:
                            singlecalc('TSS', strand, distance, table, chromosome, region_outside, start, end)
                            singlecalc('TSS', strand, distance, table, chromosome, region_inside, start, end)

                        elif windowsize:
                            windowcalc('TSS', strand, distance, windowsize, table, chromosome, region_outside, start, end)
                            windowcalc('TSS', strand, distance, windowsize, table, chromosome, region_inside, start, end)

                    elif strand == '-':
                        region_outside = (end, end + distance)
                        region_inside = (end - distance, end)
                        
                        if not windowsize:
                            singlecalc('TSS', strand, distance, table, chromosome, region_outside, start, end)
                            singlecalc('TSS', strand, distance, table, chromosome, region_inside, start, end)

                        elif windowsize:
                            windowcalc('TSS', strand, distance, windowsize, table, chromosome, region_outside, start, end)
                            windowcalc('TSS', strand, distance, windowsize, table, chromosome, region_inside, start, end)

                elif 'three' in utr:
                    if strand == '+':
                        region_outside = (end, end + distance)
                        region_inside = (end - distance, end)

                        if not windowsize:
                            singlecalc('TES', strand, distance, table, chromosome, region_outside, start, end)
                            singlecalc('TES', strand, distance, table, chromosome, region_inside, start, end)

                        elif windowsize:
                            windowcalc('TES', strand, distance, windowsize, table, chromosome, region_outside, start, end)
                            windowcalc('TES', strand, distance, windowsize, table, chromosome, region_inside, start, end)

                    elif strand == '-':
                        region_outside = (start - distance, start)
                        region_inside = (start, start + distance)
                        
                        if not windowsize:
                            singlecalc('TES', strand, distance, table, chromosome, region_outside, start, end)
                            singlecalc('TES', strand, distance, table, chromosome, region_inside, start, end)

                        elif windowsize:
                            windowcalc('TES', strand, distance, windowsize, table, chromosome, region_outside, start, end)
                            windowcalc('TES', strand, distance, windowsize, table, chromosome, region_inside, start, end)
