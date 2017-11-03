'''
test_distal.py - correlate ldhelmet results to distance from closest chromosome end

usage:
python3.5 test_genic_utr.py -l [ldhelmetfile (.txt)] -c chromosome_5 -a annotation_table.txt.gz > chromosome_5.txt

AH - 11/2017
'''

import ant
import argparse
from tqdm import tqdm

parser = argparse.ArgumentParser(description = 'Returns flat file with distance to closest end of chr + recombination rate'
                                usage = 'test_distal.py [options]')

parser.add_argument('-l', '--ldhelmetfile', required = True,
                   type = str, help = 'LDhelmet file')
parser.add_argument('-c', '--chrom', required = True,
                   type = str, help = 'Chromosome (formatted according to the annotation table).')
parser.add_argument('-a', '--annotation', required = True,
                   type = str, help = 'Annotation table. Must have corresponding .tbi file.')
parser.add_argument('-d', '--header', required = False,
                   action = 'store_true', help = 'Add header to top of file. Optional - default False.')

args = parser.parse_args()
ldhelmetfile = str(args.ldhelmetfile)
chrom = str(args.chrom)
annotation = str(args.annotation)
header = args.header

lengths = {'chromosome_1': 8033585, # hardcoded for chlamy
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

chromname = ldhelmetfile[:-7] # file name is chromosome_x_50.txt 
# though hardcoding is probably less than ideal...

length = lengths[chromname]
midpoint = int(length / 2) # rounds down

if header:
    print('distance rho')

with open(ldhelmetfile, 'r') as f:
    for line in tqdm(f):
        if not line.split(' ')[0].startswith(('#', 'ver')): # filter LDhelmet preamble
            line_content = line.split(' ')
            line_start, line_end, rho = int(line_content[0]), int(line_content[1]), float(line_content[2])
            p = ant.Reader(annotation).fetch(chrom = chrom, start = line_start, end = line_end)
            for record in p:
                if record.pos < midpoint: # closer to start
                    print(record.pos, rho)
                elif record.pos > midpoint: # closer to end
                    recdist = length - record.pos
                    print(recdist, rho)
