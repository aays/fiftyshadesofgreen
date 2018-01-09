'''where are hotspots located?

script should take in a flat file, select columns of interest, and loop through the annotation table using the boolean attributes

python3.5 hotspot_annotation.py

'''

import sys
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

correlates = dict.fromkeys(['is_genic', 'is_intergenic', 'is_exonic', 'is_intronic', 'is_utr5', 'is_utr3', 'is_in_CDS'], 0)

with open(filename, 'r') as f:
	for line in tqdm(f):
		if 'chr,block_start,block_end' in line: # skip header
			continue
		else:
			sp = line.rstrip().split(',')
			chrom, start, end = sp[1:4]

			p = antr.Reader(table)

			for record in p.fetch(chrom, start, end):
				for key in correlates:
					if getattr(record, key):
						correlates[key] += 1

print('correlate count')
for key in correlates:
	print(key, correlates[key])
    
