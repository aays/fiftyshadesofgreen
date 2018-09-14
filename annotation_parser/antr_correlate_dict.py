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
def args():
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
    parser.add_argument('-x', '--gene_context', required = False,
                       type = int, help = 'Report rho for regions of x kb downstream and upstream of all genes. (Optional)')

    args = parser.parse_args()

    if not args.correlates and not args.gc_content:
        print('Please provide one or more correlates.')
        print('This could also be just GC content (use --gc_content without providing anything to -c)')
        sys.exit(1)
    
    return [args.table, int(args.windowsize), args.correlates, 
            args.gc_content, args.gene_context]

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
    if attribute == 'upstream' or attribute == 'downstream':
        return None
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

def check_gene_proximity(record, dist, direction):
    '''(rec, int, str) -> bool
    dir = 'u' for upstream, 'd' for downstream
    '''
    p = antr.Reader(table)

    try: # in case site hits end of chrom
        if direction == 'u':
            region = p.fetch(record.chrom, record.pos, record.pos + dist)
        elif direction == 'd':
            region = p.fetch(record.chrom, record.pos - dist, record.pos)
        else:
            print('Invalid argument provided to direction.')
            print('Valid arguments are: u (upstream) and d (downstream)')
        for record in region:
            if record.is_genic:
                out = True
                break
            else:
                continue
        if out:
            return True
        elif not out:
            return False
    except:
        return False

def main(table, windowsize, correlates, gc_content, gene_context):
    # print column headers
    if correlates:
        if context_size: # ie upstream/downstream of genes
            correlates.extend(['upstream', 'downstream', 'both'])
        title1 = ' '.join([item + '_total' for item in correlates])
        title2 = ' '.join([item + '_count' for item in correlates])
        if gc:
            print('chromosome', 'start', 'end', title1, title2, 'GC', 'count')
        elif not gc:
            print('chromosome', 'start', 'end', title1, title2, 'count')    
    elif gc and not correlates:
        print('chromosome', 'start', 'end', 'GC', 'rho', 'rho_total', 'count')

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
                        if key in ['upstream', 'downstream', 'both']:
                            continue

                        elif attr_fetch(record, key) and not record.ld_rho == 'NA':
                            if key == 'intergenic' and attr_fetch(record, 'intergenic') and context_size:
                                neither = True
                                upstream = False
                                downstream = False

                                if check_gene_proximity(record, context_size, 'u'):
                                    neither = False
                                    upstream = True
                                if check_gene_proximity(record, context_size, 'd'):
                                    neither = False
                                    downstream = True

                                if upstream and downstream:
                                    rho['both'] += record.ld_rho
                                    count['both'] += 1
                                    total_counter += 1
                                    continue # don't class as intergenic
                                elif upstream and not downstream:
                                    rho['upstream'] += record.ld_rho
                                    count['upstream'] += 1
                                    total_counter += 1
                                    continue
                                elif downstream and not upstream:
                                    rho['downstream'] += record.ld_rho
                                    count['downstream'] += 1
                                    total_counter += 1
                                    continue
                                elif neither: # continue to code below
                                    pass

                            rho[key] += record.ld_rho
                            count[key] += 1
                            total_counter += 1
                        else:
                            continue

                rhovals = list(rho.values())
                countvals = list(count.values())

                totals = ' '.join([str(v) for v in rhovals])
                counts = ' '.join([str(v) for v in countvals])

            if gc: # gc content option selected
                gc_rho = 0.0
                gc_counter = 0

                for record in tqdm(p.fetch(current_chrom, window[0], window[1])):
                    gc_rho += record.ld_rho
                    gc_counter += 1

                gc_window = gc_calc(current_chrom, window, table)
                gc_rho_perbp = gc_rho / gc_counter

                if correlates:
                    print(current_chrom, window[0], window[1], totals, counts, gc_window, total_counter)
                elif not correlates:
                    print(current_chrom, window[0], window[1], gc_window, gc_rho_perbp, gc_rho, gc_counter)
            elif not gc:
                print(current_chrom, window[0], window[1], totals, counts, total_counter)

if __name__ == 'main':
    arguments = args()
    main(*arguments)
