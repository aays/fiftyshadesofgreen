'''
modification of antr_correlate_dict.py for intergenic regions only
that allows for examination of regions upstream/downstream of genes

AH - 02/2018
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
parser.add_argument('-x', '--gene_context', required = False,
                   type = int, help = 'Report rho for regions of x kb downstream and upstream of all genes. (Optional)')

args = parser.parse_args()

table = args.table
windowsize = int(args.windowsize)
context_size = args.gene_context

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

def check_gene_proximity(record, dist, direction):
    '''(rec, int, str) -> bool
    direction = 'u' for upstream, 'd' for downstream
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


# print column headers
print('chromosome', 'start', 'end', 'upstream_total', 'downstream_total', \
    'intergenic_total', 'both_total', 'upstream_count', 'downstream_total', \
    'intergenic_total', 'both_total', 'total_counter')

for chrom in range(1, 18):
    current_chrom = 'chromosome_{}'.format(str(chrom))
    windows = list(range(0, lengths[current_chrom], windowsize)) + [lengths[current_chrom]]

    for i in range(len(windows) - 1):
        window = (windows[i], windows[i + 1])

        rho = OrderedDict.fromkeys(['upstream', 'downstream', 'intergenic', 'both'], 0.0)
        count = OrderedDict.fromkeys(['upstream', 'downstream', 'intergenic', 'both'], 0)
        total_counter = 0

        for record in tqdm(p.fetch(current_chrom), window[0], window[1]):
            if record.is_intergenic:
                total_counter += 1
                if check_gene_proximity(record, context_size, 'u'):
                    rho['upstream'] += record.ld_rho
                    count['upstream'] += 1
                if check_gene_proximity(record, context_size, 'd'):
                    rho['downstream'] += record.ld_rho
                    count['downstream'] += 1
                if check_gene_proximity(record, context_size, 'u') and check_gene_proximity(record, context_size, 'd'):
                    rho['both'] += record.ld_rho
                    count['both'] += 1
                elif not check_gene_proximity(record, context_size, 'u') and not check_gene_proximity(record, context_size, 'd'):
                    rho['intergenic'] += record.ld_rho
                    count['intergenic'] += 1
            else:
                continue

        rhovals = list(rho.values())
        countvals = list(count.values())

        totals = ' '.join([str(v) for v in rhovals])
        counts = ' '.join([str(v) for v in countvals])

        print(current_chrom, window[0], window[1], totals, counts, total_counter)
