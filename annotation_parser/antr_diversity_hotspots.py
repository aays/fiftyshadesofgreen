'''
antr_diversity_hotspots.py - calculate nucleotide diversity (theta_pi) in and out of recombination hotspots

requires a specially made csv file that modifies find_hotspots.py's output to include whether windows are hotspots or not

usage:
python antr_diversity_hotspots.py --table [table (.txt.gz)] --dist [dist file (.csv/.txt)] > [outfile]

AH - 04/2018
'''

from tqdm import tqdm
from ness_vcf import SFS
import antr
import argparse

parser = argparse.ArgumentParser(description = 'Get diversity (theta pi) in and out of hotspots.',
                                 usage = 'antr_diversity_hotspots.py -t [table] -d [dist table] > [outfile]')
parser.add_argument('-t', '--table', required = True,
                    type = str, help = 'Annotation table file (.txt.gz)')
parser.add_argument('-d', '--dist', required = True,
                    type = str, help = 'File containing rho LD distribution with hotspot column')
parser.add_argument('-m', '--min_alleles', required = False,
                    type = int, help = 'Filter sites with less than n alleles called. [Optional]')
parser.add_argument('-n', '--neutral_only', required = False,
                    action = 'store_true', help = 'Only consider neutral (intergenic, intronic, 4-fold degenerate) sites? [Optional]')

args = parser.parse_args()

table = str(args.table)
dist = str(args.dist)

if args.min_alleles:
    min_alleles = args.min_alleles
else:
    min_alleles = None
    
if args.neutral_only:
    neutral_only = True

def MAF_from_allele_count(allele_counts, min_alleles = None):
    minor_allele_count = sorted(allele_counts)[-2] # second most common allele count
    total_alleles_called = sum(allele_counts)
    if min_alleles and total_alleles_called <= min_alleles:
        return None
    try:
        MAF = minor_allele_count / float(total_alleles_called)
        return (MAF, total_alleles_called)
    except ZeroDivisionError:
        return None

def SFS_from_antr(table, chromosome, start, end, min_alleles = None, neutral_only = False, counter = False):
    SFSs = {}
    p = antr.Reader(table)
    if counter:
        record_count = 0
    for record in p.fetch(chromosome, start, end):
        # diversity calc
        allele_counts = record.quebec_alleles
        if neutral_only and True not in [record.is_intergenic, record.is_intronic, record.is_fold4]:
            continue
        try:
            MAF, total_alleles_called = MAF_from_allele_count(allele_counts, min_alleles = min_alleles)
        except TypeError:
            continue
        if min_alleles and total_alleles_called < min_alleles: # filter sites that don't have enough called alleles
            continue
        if total_alleles_called not in SFSs:
            SFSs[total_alleles_called] = SFS([0]*(total_alleles_called + 1))
        SFSs[total_alleles_called].add(MAF, total_alleles_called)
        if counter:
            record_count += 1
    diversity = sum([sfs.theta_pi() * sfs.sites() for sfs in SFSs.values()]) / sum([sfs.sites() for sfs in SFSs.values()])
    if not counter:
        return diversity
    elif counter:
        return diversity, record_count

# column headers
print('chromosome start end diversity record_count is_hotspot')

with open(dist, 'r') as f:
    for line in tqdm(f):
        if line.startswith('chr,block_start'): # header
            continue
        else:
            sp = [l.rstrip('\n') for l in line.split(',')]
            chrom, start, end, rho, chrom_mean, actual_ratio, is_hotspot = sp
            try:
                curr_theta, c = SFS_from_antr(table, chrom, int(start), int(end), 
                                              min_alleles = min_alleles, neutral_only = neutral_only, counter = True)
                c = int(c)
            except ZeroDivisionError: # nothing in window
                curr_theta, c = 0, 0
            print(chrom, start, end, curr_theta, c, is_hotspot)
