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
    rec_attr = [item for item in dir(rec) if '__' not in item and attribute in item]
    try:
        assert len(rec_attr) == 1
    except:
        raise AssertionError('{} is not a valid attribute.'.format(attribute))
    rec_attr = rec_attr[0] # extract item from list
    out = getattr(rec, rec_attr)
    return out

## single correlate
if len(correlates) == 1:
    
    corr = correlates[0]
    print('chromosome', 'start', 'end', corr, corr + '_total')

    for chrom in range(1, 18):
        current_chrom = 'chromosome_{}'.format(str(chrom))
        windows = list(range(0, lengths[current_chrom], windowsize)) + [lengths[current_chrom]]

        p = antr.Reader(table)

        for i in range(len(windows) - 1):
            window = (windows[i], windows[i + 1])
            corr_rho_out = 0.0
            corr_count = 0 # no need for a 'total counter' - only one correlate

            for record in tqdm(p.fetch(current_chrom, window[0], window[1])):
                if not record.ld_rho == 'NA' and attr_fetch(record, corr):
                    corr_rho_out += record.ld_rho
                    corr_count += 1
                else:
                    continue
            print(current_chrom, window[0], window[1], corr_rho_out / corr_count, corr_rho_out)

## two correlates
if len(correlates) == 2:
    
    corr1, corr2 = correlates # unpack list
    print('chromosome', 'start', 'end', corr1, corr2, corr1 + '_total', corr2 + '_total', 'count')

    for chrom in range(1, 18):
        current_chrom = 'chromosome_{}'.format(str(chrom))
        windows = list(range(0, lengths[current_chrom], windowsize)) + [lengths[current_chrom]]

        p = antr.Reader(table) 
        
        for i in range(len(windows) - 1):
            window = (windows[i], windows[i + 1])
            corr1_rho_out = 0.0
            corr2_rho_out = 0.0
            corr1_count, corr2_count, total_counter = 0, 0, 0

            for record in tqdm(p.fetch(current_chrom, window[0], window[1])):
                if not record.ld_rho == 'NA' and attr_fetch(record, corr1):
                    corr1_rho_out += record.ld_rho
                    corr1_count += 1
                    total_counter += 1
                elif not record.ld_rho == 'NA' and attr_fetch(record, corr2):
                    corr2_rho_out += record.ld_rho
                    corr2_count += 1
                    total_counter += 1
                else:
                    continue
            print(current_chrom, window[0], window[1], corr1_rho_out / corr1_count, 
                  corr2_rho_out / corr2_count, corr1_rho_out, corr2_rho_out, total_counter)

elif len(correlates) == 3:

    corr1, corr2, corr3 = correlates 
    print('chromosome', 'start', 'end', corr1, corr2, corr3,
    corr1 + '_total', corr2 + '_total', corr3 + '_total', 'count')

    for chrom in range(1, 18):
        current_chrom = 'chromosome_{}'.format(str(chrom))
        windows = list(range(0, lengths[current_chrom], windowsize)) + [lengths[current_chrom]]

        p = antr.Reader(table) 
        
        for i in range(len(windows) - 1):
            window = (windows[i], windows[i + 1])
            corr1_rho_out = 0.0
            corr2_rho_out = 0.0
            corr3_rho_out = 0.0
            corr1_count, corr2_count, corr3_count, total_counter = 0, 0, 0, 0

            for record in tqdm(p.fetch(current_chrom, window[0], window[1])):
                if not record.ld_rho == 'NA' and attr_fetch(record, corr1):
                    corr1_rho_out += record.ld_rho
                    corr1_count += 1
                    total_counter += 1
                elif not record.ld_rho == 'NA' and attr_fetch(record, corr2):
                    corr2_rho_out += record.ld_rho
                    corr2_count += 1
                    total_counter += 1
                elif not record.ld_rho == 'NA' and attr_fetch(record, corr3):
                    corr3_rho_out += record.ld_rho
                    corr3_count += 1
                    total_count += 1
                else:
                    continue
            print(current_chrom, window[0], window[1], corr1_rho_out / corr1_count, 
                  corr2_rho_out / corr2_count, corr3_rho_out / corr3_count, corr1_rho_out, 
                  corr2_rho_out, corr3_rho_out, total_counter)

elif len(correlates) == 4:

    corr1, corr2, corr3, corr4 = correlates 
    print('chromosome', 'start', 'end', corr1, corr2, corr3, corr4,
    corr1 + '_total', corr2 + '_total', corr3 + '_total', corr4 + '_total', 'count')

    for chrom in range(1, 18):
        current_chrom = 'chromosome_{}'.format(str(chrom))
        windows = list(range(0, lengths[current_chrom], windowsize)) + [lengths[current_chrom]]

        p = antr.Reader(table) 
        
        for i in range(len(windows) - 1):
            window = (windows[i], windows[i + 1])
            corr1_rho_out = 0.0
            corr2_rho_out = 0.0
            corr3_rho_out = 0.0
            corr4_rho_out = 0.0
            corr1_count, corr2_count, corr3_count, corr4_count, total_counter = 0, 0, 0, 0, 0

            for record in tqdm(p.fetch(current_chrom, window[0], window[1])):
                if not record.ld_rho == 'NA' and attr_fetch(record, corr1):
                    corr1_rho_out += record.ld_rho
                    corr1_count += 1
                    total_counter += 1
                elif not record.ld_rho == 'NA' and attr_fetch(record, corr2):
                    corr2_rho_out += record.ld_rho
                    corr2_count += 1
                    total_counter += 1
                elif not record.ld_rho == 'NA' and attr_fetch(record, corr3):
                    corr3_rho_out += record.ld_rho
                    corr3_count += 1
                    total_count += 1
                elif not record.ld_rho == 'NA' and attr_fetch(record, corr3):
                    corr4_rho_out += record.ld_rho
                    corr4_rho_out += 1
                    total_count += 1
                else:
                    continue
            print(current_chrom, window[0], window[1], corr1_rho_out / corr1_count, 
                  corr2_rho_out / corr2_count, corr3_rho_out / corr3_count, corr4_rho_out / corr4_count, 
                  corr1_rho_out, corr2_rho_out, corr3_rho_out, corr4_rho_out, total_counter)

elif len(correlates) > 4:
    raise InputError('Too many correlates. Please provide 1-4 correlates as input.')




