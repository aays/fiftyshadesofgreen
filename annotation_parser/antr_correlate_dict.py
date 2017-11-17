'''
takes in a table + windowsize + list of correlates, and returns mean rho values / window for records
satisfying a given correlate.

usage:
python3.5 antr_correlate.py -t [table (.txt.gz)] -w [windowsize] -c [correlates] > output.txt
python3.5 antr_correlate.py -t table.txt.gz -w 100000 -c exonic intronic > output.txt

AH - 11/2017
'''


import argparse
import sys
from tqdm import tqdm
from collections import OrderedDict

try:
    import antr
except ImportError:
    sys.path.append('/scratch/research/projects/chlamydomonas/genomewide_recombination/analysis/fiftyshadesofgreen/annotation_parser/')


parser = argparse.ArgumentParser(description = 'General purpose calculation of rho w/ attributes in defined windows.',
                                usage = 'antr_correlation.py [options]')

parser.add_argument('-t', '--table', required = True,
                   type = str, help = 'Annotation table file (.txt.gz)')
parser.add_argument('-w', '--windowsize', required = True,
                   type = int, help = 'Window size')
parser.add_argument('-c', '--correlates', required = True,
                   type = str, nargs = '+', help = 'Space separated list of correlates. Be as exact with naming as possible.')

args = parser.parse_args()

table = args.table
windowsize = int(args.windowsize)
correlates = args.correlates

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
    rec_attr = [item for item in dir(rec) if '__' not in item]
    rec_attr = [item for item in rec_attr if attribute in item]
    try:
        assert len(rec_attr) == 1
    except:
        raise AssertionError('{} is not a valid attribute.'.format(attribute))
    rec_attr = rec_attr[0] # extract item from list
    out = getattr(rec, rec_attr)
    return out

title1 = ' '.join(correlates)
title2 = ' '.join([item + '_total' for item in correlates])
print('chromosome', 'start', 'end', title1, title2, 'count')

for chrom in range(1, 18):
    current_chrom = 'chromosome_{}'.format(str(chrom))
    windows = list(range(0, lengths[current_chrom], windowsize)) + [lengths[current_chrom]]

    p = antr.Reader(table) 

    for i in range(len(windows) - 1):
        rho = OrderedDict.fromkeys(correlates, 0.0)
        count = OrderedDict.fromkeys(correlates, 0)
        total_counter = 0

        for record in tqdm(p.fetch(current_chrom, window[0], window[1])):
            for key in rho.keys():
                if attr_fetch(record, key) and not record.ld_rho == 'NA':
                    rho[key] += record.ld_rho
                    count[key] += 1
                    total_counter += 1
                else:
                    continue
        
        rhovals = list(rho.values())
        countvals = list(count.values())

        allvals = ' '.join([str(rhovals[i] / countvals[i]) for i in range(len(rhovals))])
        totals = ' '.join([str(v) for v in rhovals])

        print(current, window[0], window[1], allvals, totals, total_counter)




