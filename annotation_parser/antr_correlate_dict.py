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

# args
parser = argparse.ArgumentParser(description = 'General purpose calculation of rho w/ attributes in defined windows.',
                                usage = 'antr_correlation.py [options]')

parser.add_argument('-t', '--table', required = True,
                   type = str, help = 'Annotation table file (.txt.gz)')
parser.add_argument('-w', '--windowsize', required = True,
                   type = int, help = 'Window size')
parser.add_argument('-c', '--correlates', required = False,
                   type = str, nargs = '+', help = 'Space separated list of correlates. Be as exact with naming as possible.')
parser.add_argument('-g', '--gc_content', required = False,
                   action = 'store_true', help = 'Include GC content in output file? (Optional)')

args = parser.parse_args()

table = args.table
windowsize = int(args.windowsize)
correlates = args.correlates
gc = args.gc_content

if not correlates and not gc:
    print('Please provide one or more correlates.')
    print('This could also be just GC content (use --gc_content without providing anything to -c)')
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

def gc_calc(chromosome, window, table):
    '''Returns GC content in given window as proportion.'''
    seq = ''.join([record.ref for record in antr.Reader(table).fetch(chromosome, window[0], window[1])])
    total = len(seq)
    GC = seq.count('G') + seq.count('C')
    GC_content = GC / total
    return GC_content

# print column headers
if correlates:
    title1 = ' '.join(correlates)
    title2 = ' '.join([item + '_total' for item in correlates])
    title3 = ' '.join([item + '_count' for item in correlates])
    if gc:
        print('chromosome', 'start', 'end', title1, title2, title3, 'GC%', 'count')
    elif not gc:
        print('chromosome', 'start', 'end', title1, title2, title3, 'count')    
elif gc and not correlates:
    print('chromosome', 'start', 'end', 'GC%', 'rho', 'rho_total', 'count')

# iterate through chromosomes
for chrom in range(1, 18):
    current_chrom = 'chromosome_{}'.format(str(chrom))
    windows = list(range(0, lengths[current_chrom], windowsize)) + [lengths[current_chrom]]

    p = antr.Reader(table)

    for i in range(len(windows) - 1):
        window = (windows[i], windows[i + 1])
        
        if correlates:
            rho = OrderedDict.fromkeys(correlates, 0.0)
            count = OrderedDict.fromkeys(correlates, 0)
            total_counter = 0
            
            # iterate through records in window
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

            try:
                allvals = ' '.join([str(rhovals[i] / countvals[i]) for i in range(len(rhovals))])
                totals = ' '.join([str(v) for v in rhovals])
                counts = ' '.join([str(v) for v in countvals])
            except ZeroDivisionError: # nothing in window
                allvals = ' '.join([str(0) for i in range(len(rhovals))])
                totals = ' '.join([str(0) for v in rhovals])
                counts = ' '.join([str(0) for v in countvals])
            
        if gc: # gc content option selected
            gc_rho = 0.0
            gc_counter = 0
            
            for record in tqdm(p.fetch(current_chrom, window[0], window[1])):
                gc_rho += record.ld_rho
                gc_counter += 1
            
            gc_window = gc_calc(current_chrom, window, table)
            gc_rho_perbp = gc_rho / gc_counter
            
            if correlates:
                print(current_chrom, window[0], window[1], allvals, totals, counts, gc_window, total_counter)
            elif not correlates:
                print(current_chrom, window[0], window[1], gc_window, gc_rho_perbp, gc_rho, gc_counter)
        elif not gc:
            print(current_chrom, window[0], window[1], allvals, totals, counts, total_counter)
