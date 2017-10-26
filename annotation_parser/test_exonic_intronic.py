'''
test_exonic_intronic.py - correlate ldhelmet results to whether sites are exonic or not

usage:
python3.5 test_exonic_intronic.py -l [ldhelmetfile (.txt)] -c chromosome_5 -a annotation_table.txt.gz > chromosome_5_exonic.txt

AH - 10/2017
'''

import ant
import argparse
from tqdm import tqdm

parser = argparse.ArgumentParser(description = 'Returns rho for exonic and intronic regions', 
                                usage = 'test_exonic_intronic.py [options]')

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
    for line in f:
        if line.split(' ')[0].startswith(('#', 'ver')):
            pass
        else:
            line_content = line.split(' ')
            line_start, line_end, rho = int(line_content[0]), int(line_content[1]), float(line_content[2])
            p = ant.Reader(annotation).fetch(chrom = chrom, start = line_start, end = line_end)
            for record in p:
                if record.is_exonic:
                    print(chrom, 'exonic', rho)
                elif record.is_intronic:
                    print(chrom, 'intronic', rho)
                else:
                    continue
            
