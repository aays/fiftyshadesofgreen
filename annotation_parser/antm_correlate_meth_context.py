'''
modification of antm_correlate_meth.py to account for context

annotation table must have both ld_rho values + methylation beta values/sequence contexts appended as additional columns.
this uses the _antm_ parser variant, not antr

usage:
python3.5 antm_correlate_meth_context.py -t [table (.txt.gz)] -w [windowsize] -c [chromosome] > output.txt
python3.5 antm_correlate_meth_context.py -t table.txt.gz -w 100000 -c chromosome_2 > output.txt

AH - 01/2018
'''

import argparse
import sys
from tqdm import tqdm
from collections import OrderedDict

try:
    import antm
except ImportError:
    sys.path.append('/scratch/research/projects/chlamydomonas/genomewide_recombination/analysis/fiftyshadesofgreen/annotation_parser/')

# args
parser = argparse.ArgumentParser(description = 'General purpose calculation of rho and other correlates in defined windows.',
                                usage = 'antm_correlate_meth.py [options]')

parser.add_argument('-t', '--table', required = True,
                   type = str, help = 'Annotation table file (.txt.gz)')
parser.add_argument('-w', '--windowsize', required = True,
                   type = int, help = 'Window size')
parser.add_argument('-c', '--chromosome', required = True,
                   type = str, help = 'Chromosome name as it appears in the table.')
parser.add_argument('-a', '--all_sites', required = False,
                   action = 'store_true', help = 'Include if considering _all_ sites, not just methylated ones')

args = parser.parse_args()

table = args.table
windowsize = int(args.windowsize)
current_chrom = args.chromosome
all_sites = args.all_sites

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

# print column headers

print('chromosome', 'start', 'end', 'rho_values', 'rho_count', 'methylation_values', 'methylation_count',
      'CG_values', 'CG_count', 'CHG_values', 'CHG_count', 'CHH_values', 'CHH_counts',
      'CN_values', 'CN_count', 'record_count')

windows = list(range(0, lengths[current_chrom], windowsize)) + [lengths[current_chrom]]

p = antm.Reader(table)

for i in tqdm(range(len(windows) - 1)):
    window = (windows[i], windows[i + 1])
    rho = 0.0
    meth = 0.0
    rho_count = 0
    meth_count = 0
    record_counter = 0
    meth_dict = dict.fromkeys(['CG', 'CHG', 'CHH', 'CN'], 0.0)
    meth_context_counts = dict.fromkeys(['CG', 'CHG', 'CHH', 'CN'], 0)

    # iterate through records in window
    for record in p.fetch(current_chrom, window[0], window[1]):
        if record.ld_rho != 'NA' and all_sites:
            if record.methylation:
                beta_val, context = record.methylation
                meth += beta_val # all methylation
                rho += record.ld_rho
                meth_dict[context] += beta_val # per-context methylation
                meth_context_counts[context] += 1
                rho_count += 1
                meth_count += 1
            elif not record.methylation:
                rho += record.ld_rho
                rho_count += 1
            record_counter += 1
        elif record.ld_rho != 'NA' and not all_sites:
            if record.methylation:
                beta_val, context = record.methylation
                meth += beta_val
                rho += record.ld_rho
                meth_dict[context] += beta_val
                meth_context_counts[context] += 1
                rho_count += 1
                meth_count += 1
            record_counter += 1

    print(current_chrom, window[0], window[1], rho, rho_count, meth, meth_count, 
          meth_dict['CG'], meth_context_counts['CG'], meth_dict['CHG'], meth_context_counts['CHG'],
          meth_dict['CHH'], meth_context_counts['CHH'], meth_dict['CN'], meth_context_counts['CN'], record_counter)

