'''
takes in a table + windowsize, and returns mean rho values / window + nucleotide diversity / window

usage:
python3.5 antr_correlate_diversity.py -t [table (.txt.gz)] -w [windowsize] > output.txt
python3.5 antr_correlate_diversity.py -t table.txt.gz -w 100000 > output.txt

AH - 12/2017
'''

import argparse
import sys
from tqdm import tqdm
from collections import OrderedDict
from ness_vcf import SFS

try:
    import antr
except ImportError:
    sys.path.append('/scratch/research/projects/chlamydomonas/genomewide_recombination/analysis/fiftyshadesofgreen/annotation_parser/')

# args
parser = argparse.ArgumentParser(description = 'General purpose calculation of rho and other correlates in defined windows.',
                                usage = 'antr_correlate_diversity.py [options]')

parser.add_argument('-t', '--table', required = True,
                   type = str, help = 'Annotation table file (.txt.gz)')
parser.add_argument('-w', '--windowsize', required = True,
                   type = int, help = 'Window size')
parser.add_argument('-m', '--min_alleles', required = False,
                    type = int, help = 'Filter sites with less than n alleles called. [Optional]')
parser.add_argument('-n', '--neutral_only', required = False,
                    action = 'store_true', help = 'Only consider neutral (intergenic, intronic, 4-fold degenerate) sites? [Optional]')

args = parser.parse_args()

table = str(args.table)
windowsize = int(args.windowsize)
min_alleles = args.min_alleles
neutral_only = args.neutral_only

# chromosome lengths - hardcoded for chlamy
lengths = {'chromosome_1': 8033585,
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

def attr_fetch(rec, attribute):
    '''(rec, str) -> bool/float
    Used for fetching desired attributes from a record.'''
    rec_attr = [item for item in dir(rec) if '__' not in item and attribute in item]
    try:
        assert len(rec_attr) == 1
    except:
        raise AssertionError('{} is not a valid attribute. {} matches found - {}'.format(attribute, len(rec_attr), rec_attr))
    rec_attr = rec_attr[0] # extract item from list
    out = getattr(rec, rec_attr)
    return out

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

def SFS_from_antr(table, chromosome, start, end, min_alleles = None, neutral_only = False):
    SFSs = {}
    p = antr.Reader(table)
    for record in tqdm(p.fetch(chromosome, start, end)):
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
    diversity = sum([sfs.theta_pi() * sfs.sites() for sfs in SFSs.values()]) / sum([sfs.sites() for sfs in SFSs.values()])
    return diversity

# print column headers
print('chromosome', 'start', 'end', 'diversity', 'rho', 'rho_total', 'rho_count', 'iter_count')

for chrom in range(1, 18):
    current_chrom = 'chromosome_{}'.format(str(chrom))
    windows = list(range(0, lengths[current_chrom], windowsize)) + [lengths[current_chrom]]

    p = antr.Reader(table)

    for i in range(len(windows) - 1):
        window = (windows[i], windows[i + 1])
        rho = 0.0
        count = 0
        record_counter = 0

        # iterate through records in window
        for record in tqdm(p.fetch(current_chrom, window[0], window[1])):
            if record.ld_rho != 'NA':
                rho += record.ld_rho
                count += 1
            else:
                continue
            record_counter += 1

        try:
            rho_out = rho / count
            curr_div = SFS_from_antr(table, current_chrom, window[0], window[1], min_alleles = min_alleles, neutral_only = neutral_only)
        except ZeroDivisionError: # nothing in window
            rho_out = 0
            curr_div = 0

        print(current_chrom, window[0], window[1], curr_div, rho_out, rho, count, record_counter)
