'''
test_genic_utr.py - correlate ldhelmet results to whether sites are genic or not

usage:
python3.5 test_genic_utr.py -l [ldhelmetfile (.txt)] -c chromosome_5 -a annotation_table.txt.gz > chromosome_5.txt

AH - 10/2017
'''

import ant
import argparse
from tqdm import tqdm

parser = argparse.ArgumentParser(description = 'Returns rho for genic (UTR3 and UTR5) and intergenic regions', 
                                usage = 'test_genic_utr.py [options]')

parser.add_argument('-l', '--ldhelmetfile', required = True,
                   type = str, help = 'LDhelmet file')
parser.add_argument('-c', '--chrom', required = True,
                   type = str, help = 'Chromosome (formatted according to the annotation table).')
parser.add_argument('-a', '--annotation', required = True,
                   type = str, help = 'Annotation table. Must have corresponding .tbi file.')

args = parser.parse_args()
ldhelmetfile = str(args.ldhelmetfile)
chrom = str(args.chrom)
annotation = str(args.annotation)

with open(ldhelmetfile, 'r') as f:
    for line in tqdm(f):
        if line.split(' ')[0].startswith(('#', 'ver')):
            pass
        else:
            line_content = line.split(' ')
            line_start, line_end, rho = int(line_content[0]), int(line_content[1]), float(line_content[2])
            p = ant.Reader(annotation).fetch(chrom = chrom, start = line_start, end = line_end)
            for record in p:
                if record.is_genic:
                    if record.is_utr3:
                        print(chrom, 'genic_utr3', rho)
                    elif record.is_utr5:
                        print(chrom, 'genic_utr5', rho)
                    elif not record.is_utr3 and not record.is_utr5:
                        print(chrom, 'genic', rho)
                elif record.is_intergenic:
                    print(chrom, 'intergenic', rho)
                else:
                    continue
            
