'''
ant_correlation.py - correlate ldhelmet results to genomic features

usage:
python3.5 ant_correlation.py -l [ldhelmetfile (.txt)] -c chromosome_5 -a annotation_table.txt.gz > chromosome_5_ant.txt

warning: horrendously, horrifically slow. would avoid using, in all honesty

AH - 10/2017
'''


import ant
import argparse
from tqdm import tqdm

parser = argparse.ArgumentParser(description = 'Correlate stuff b/w LDhelmet output + annotation table.', 
                                usage = 'ant_correlation.py [options]')

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
    lines = f.readlines()[3:] # skip preamble
    
def linetodict(line):
    line = line.split(' ')
    start, end = int(line[0]), int(line[1])
    return range(start, end), float(line[2]) # range + mean rho

rholist = [linetodict(line) for line in tqdm(lines)]
ldh_dict = dict(rholist)
del(rholist) # clear up

# introduce annotation table

def checkrho(record, ld_dict):
    currentrho = 0.0
    for key in ld_dict.keys():
        if record.pos in key:
            currentrho = ld_dict[key]
            continue
        elif record.pos not in key:
            pass
    return currentrho

print('chrom type avgrho numrecords')

p = ant.Reader(annotation).fetch(chrom)
exonic_rho = [checkrho(r, ldh_dict) for r in p if r.is_exonic]
were_exonic = len(exonic_rho)
exonic_rho = sum(exonic_rho) / were_exonic # overwrites massive list
print(chrom, 'exonic', exonic_rho, were_exonic)

p = ant.Reader(annotation).fetch(chrom)
intronic_rho = [checkrho(r, ldh_dict) for r in p if r.is_intronic]
were_intronic = len(intronic_rho)
intronic_rho = sum(intronic_rho) / were_intronic
print(chrom, 'intronic', intronic_rho, were_intronic)

p = ant.Reader(annotation).fetch(chrom)
genic_rho = [checkrho(r, ldh_dict) for r in p if r.is_genic]
were_genic = len(genic_rho)
genic_rho = sum(genic_rho) / were_genic
print(chrom, 'genic', genic_rho, were_genic)

p = ant.Reader(annotation).fetch(chrom)
intergenic_rho = [checkrho(r, ldh_dict) for r in p if r.is_intergenic]
were_intergenic = len(intergenic_rho)
intergenic_rho = sum(intergenic_rho) / were_intergenic
print(chrom, 'intergenic', intergenic_rho, were_intergenic)

p = ant.Reader(annotation).fetch(chrom)
CDS_rho = [checkrho(r, ldh_dict) for r in p if r.is_in_CDS]
were_CDS = len(CDS_rho)
CDS_rho = sum(CDS_rho) / were_CDS
print(chrom, 'CDS', CDS_rho, were_CDS)

p = ant.Reader(annotation).fetch(chrom)
fold0_rho = [checkrho(r, ldh_dict) for r in p if r.is_fold0]
were_fold0 = len(fold0_rho)
fold0_rho = sum(fold0_rho) / were_fold0
print(chrom, 'fold0', fold0_rho, were_fold0)

p = ant.Reader(annotation).fetch(chrom)
fold2_rho = [checkrho(r, ldh_dict) for r in p if r.is_fold2]
were_fold2 = len(fold2_rho)
fold2_rho = sum(fold2_rho) / were_fold2
print(chrom, 'fold2', fold2_rho, were_fold2)

p = ant.Reader(annotation).fetch(chrom)
fold3_rho = [checkrho(r, ldh_dict) for r in p if r.is_fold3]
were_fold3 = len(fold3_rho)
fold3_rho = sum(fold3_rho) / were_fold3
print(chrom, 'fold3', fold3_rho, were_fold3)

p = ant.Reader(annotation).fetch(chrom)
fold4_rho = [checkrho(r, ldh_dict) for r in p if r.is_fold4]
were_fold4 = len(fold4_rho)
fold4_rho = sum(fold4_rho) / were_fold4
print(chrom, 'fold4', fold4_rho, were_fold4)

p = ant.Reader(annotation).fetch(chrom)
utr3_rho = [checkrho(r, ldh_dict) for r in p if r.is_utr3]
were_utr3 = len(utr3_rho)
utr3_rho = sum(utr3_rho) / were_utr3
print(chrom, 'utr3', utr3_rho, were_utr3)

p = ant.Reader(annotation).fetch(chrom)
utr5_rho = [checkrho(r, ldh_dict) for r in p if r.is_utr5]
were_utr5 = len(utr5_rho)
utr5_rho = sum(utr5_rho) / were_utr5
print(chrom, 'utr5', utr5_rho, were_utr5)
