'''
where are hotspots located?

script should take in a flat file, select columns of interest, and loop through the annotation table using the boolean attributes

python3.5 hotspot_annotation.py -t [table (.txt.gz)] -f [hotspot file (.txt/csv)] > [out (.txt)]

AH - 01/2018
'''

import sys
import argparse
from tqdm import tqdm

try:
    import antr
except ImportError:
    sys.path.append('/scratch/research/projects/chlamydomonas/genomewide_recombination/analysis/fiftyshadesofgreen/annotation_parser/')

# args
parser = argparse.ArgumentParser(description = 'Uses annotation table to examine where hotspots are preferentially located in the genome.',
                                usage = 'hotspot_annotation.py [options]')

parser.add_argument('-t', '--table', required = True,
                   type = str, help = 'Annotation table file (.txt.gz)')
parser.add_argument('-f', '--filename', required = False,
                   type = str, help = 'csv file, exported from Singhal hotspot detection script.')

args = parser.parse_args()
table = args.table
filename = args.filename

def check_gene_proximity(record, dist, direction):
    '''(rec, int, str) -> bool
    dir = 'u' for upstream, 'd' for downstream
    '''
    p = antr.Reader(table)

    try: # in case site hits end of chrom
        if direction == 'u':
            region = p.fetch(record.chrom, record.pos, record.pos + dist)
        elif direction == 'd':
            region = p.fetch(record.chrom, record.pos - dist, record.pos)
        else:
            print('Invalid argument provided to direction.')
            print('Valid arguments are: u (upstream) and d (downstream)')
        for record in region:
            if record.is_genic:
                out = True
                break
            else:
                continue
        if out:
            return True
        elif not out:
            return False
    except:
        return False

correlates = dict.fromkeys(['is_genic', 'is_intergenic', 'is_exonic', 'is_intronic',
                            'is_utr5', 'is_utr3', 'is_in_CDS', 'upstream', 'downstream', 'both'], 0)

with open(filename, 'r') as f:
    for line in f:
        if 'chr,block_start,block_end' in line: # skip header
            continue
        else:
            sp = line.rstrip().split(',')
            chrom, start, end = sp[0], float(sp[1]), float(sp[2])

            p = antr.Reader(table)

            for record in tqdm(p.fetch(chrom, start, end)):
                for key in correlates:
                    if key in ['upstream', 'downstream', 'both']:
                        continue
                    if key != 'is_intergenic' and getattr(record, key):
                        correlates[key] += 1
                    elif key == 'is_intergenic' and getattr(record, key):
                        upstream = False
                        downstream = False
                        neither = True
                        if check_gene_proximity(record, 2000, 'u'):
                            upstream = True
                            neither = False
                        if check_gene_proximity(record, 2000, 'd'):
                            downstream = True
                            neither = False
                        if upstream and downstream:
                            correlates['both'] += 1
                        elif upstream and not downstream:
                            correlates['upstream'] += 1
                        elif downstream and not upstream:
                            correlates['downstream'] += 1
                        elif neither and not upstream and not downstream:
                            correlates['is_intergenic'] += 1

print('correlate count')
for key in correlates:
    print(key, correlates[key])
