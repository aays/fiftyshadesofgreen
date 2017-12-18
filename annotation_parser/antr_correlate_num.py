'''
takes in a table + windowsize + list of correlates, and returns mean rho values / window + correlate values / window
for downstream correlations. distinct from correlate_dict in that this accounts for numerical data, not just T/F situations.

usage:
python3.5 antr_correlate_num.py -t [table (.txt.gz)] -w [windowsize] -c [correlates] > output.txt
python3.5 antr_correlate_num.py -t table.txt.gz -w 100000 -c map_rho > output.txt

AH - 12/2017
'''

import argparse
import sys
from tqdm import tqdm
from collections import OrderedDict

try:
    import antr
except ImportError:
    sys.path.append('/scratch/research/projects/chlamydomonas/genomewide_recombination/analysis/fiftyshadesofgreen/annotation_parser/')

# args
parser = argparse.ArgumentParser(description = 'General purpose calculation of rho and other correlates in defined windows.',
                                usage = 'antr_correlate_num.py [options]')

parser.add_argument('-t', '--table', required = True,
                   type = str, help = 'Annotation table file (.txt.gz)')
parser.add_argument('-w', '--windowsize', required = True,
                   type = int, help = 'Window size')
parser.add_argument('-c', '--correlates', required = False,
                   type = str, nargs = '+', help = 'Space separated list of correlates. Be as exact with naming as possible.')

args = parser.parse_args()

table = args.table
windowsize = int(args.windowsize)
correlates = args.correlates

if not correlates:
    print('Please provide one or more correlates.')
    sys.exit(1)

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

title1 = ' '.join([item + '_values' for item in correlates])
title2 = ' '.join([item + '_values_totals' for item in correlates])
title3 = ' '.join([item + '_rho' for item in correlates])
title4 = ' '.join([item + '_rho_totals' for item in correlates])
title5 = ' '.join([item + '_count' for item in correlates])
print('chromosome', 'start', 'end', title1, title2, title3, \
    title4, title5, 'iter_count', 'record_count')

for chrom in range(1, 18):
    current_chrom = 'chromosome_{}'.format(str(chrom))
    windows = list(range(0, lengths[current_chrom], windowsize)) + [lengths[current_chrom]]

    p = antr.Reader(table)

    for i in range(len(windows) - 1):
        window = (windows[i], windows[i + 1])

        if correlates:
            rho = OrderedDict.fromkeys(correlates, 0.0)
            corr = OrderedDict.fromkeys(correlates, 0.0)
            count = OrderedDict.fromkeys(correlates, 0)
            total_counter = 0
            record_counter = 0

            # iterate through records in window
            for record in tqdm(p.fetch(current_chrom, window[0], window[1])):
                for key in correlates:
                    if attr_fetch(record, key) != 'NA' and record.ld_rho != 'NA':
                        rho[key] += record.ld_rho
                        corr[key] += attr_fetch(record, key)
                        count[key] += 1
                        total_counter += 1
                    else:
                        continue
                record_counter += 1

            corrvals = list(corr.values())
            rhovals = list(rho.values())
            countvals = list(count.values())

            try:
                corr_out = ' '.join([str(corrvals[i] / countvals[i]) for i in range(len(corrvals))])
                corr_totals = ' '.join([str(v) for v in corrvals])
                rho_out = ' '.join([str(rhovals[i] / countvals[i]) for i in range(len(rhovals))])
                rho_totals = ' '.join([str(v) for v in rhovals])
                counts = ' '.join([str(v) for v in countvals])
            except ZeroDivisionError: # nothing in window
                corr_out = ' '.join([str(0) for i in range(len(corrvals))])
                corr_totals = ' '.join([str(0) for v in corrvals])
                rho_out = ' '.join([str(0) for i in range(len(rhovals))])
                rho_totals = ' '.join([str(0) for v in rhovals])
                counts = ' '.join([str(0) for v in countvals])

            print(current_chrom, window[0], window[1], \
                corr_out, corr_totals, rho_out, rho_totals, counts, total_counter, record_counter)

