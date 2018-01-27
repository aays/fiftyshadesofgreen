'''
modification of antr_correlate_num.py specifically for methylation data

annotation table must have both ld_rho values + methylation beta values appended as additional columns.
this uses the _antm_ parser variant, not antr

usage:
python3.5 antm_correlate_meth.py -t [table (.txt.gz)] -w [windowsize] > output.txt
python3.5 antm_correlate_meth.py -t table.txt.gz -w 100000 > output.txt

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

args = parser.parse_args()

table = args.table
windowsize = int(args.windowsize)

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

print('chromosome', 'start', 'end', 'rho', 'rho_values', 'methylation', 'methylation_values', 'iter_count', 'record_count')

for chrom in range(1, 18):
    current_chrom = 'chromosome_{}'.format(str(chrom))
    windows = list(range(0, lengths[current_chrom], windowsize)) + [lengths[current_chrom]]

    p = antm.Reader(table)

    for i in range(len(windows) - 1):
        window = (windows[i], windows[i + 1])
        rho = 0.0
        meth = 0.0
        count = 0
        record_counter = 0

        # iterate through records in window
        for record in tqdm(p.fetch(current_chrom, window[0], window[1])):
            if record.ld_rho != 'NA' and record.methylation != 'NA':
                meth += record.methylation
                rho += record.ld_rho
                count += 1
            else:
                continue
            record_counter += 1

        try:
            rho_out = rho / count
            meth_out = meth / count
        except ZeroDivisionError: # nothing in window
            continue
            #rho_out = 0
            #meth_out = 0

        print(current_chrom, window[0], window[1], rho_out, rho, meth_out, meth, count, record_counter)


