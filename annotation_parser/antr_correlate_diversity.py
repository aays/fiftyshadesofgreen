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

def args():
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
    parser.add_argument('-d', '--gene_density', required = False,
                        action = 'store_true', help = 'Compute gene density per window (CDS + intron + UTR sites) [Optional]')
    parser.add_argument('-s', '--measure', required = False,
                        type = str, help = "Diversity measure to use ([theta_pi, theta_w, both]) - default is theta pi")
    parser.add_argument('-r', '--neutral_regions', required = False,
                        type = str, nargs = '+', help = '[if --neutral_only] What regions? [intergenic, intronic, fold4]. Defaults to all')

    args = parser.parse_args()
    if not args.measure:
        measure = 'theta_pi'
    else:
        measure = args.measure

    try:
        assert measure in ['theta_pi', 'theta_w', 'both']
    except:
        raise AssertionError('Invalid measure provided. Valid options include [theta_pi, theta_w, both]')
                        
    if args.neutral_regions:
        print(args.neutral_regions)
        try:
            assert 1 <= len(args.neutral_regions) <= 3
        except:
            raise AssertionError('Too many regions specified in --neutral_regions! Valid options are [intronic, intergenic, 4fold]')
    else:
        neutral_regions = None

    return [str(args.table), int(args.windowsize), args.min_alleles,
            args.neutral_only, args.gene_density, measure, args.neutral_regions]
    
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

def test_neutral(record, neutral_regions):
    '''Helper function for SFS_from_antr if custom neutral regions defined'''
    report = [getattr(record, 'is_' + region) for region in neutral_regions]
    if True not in report:
        return False
    elif True in report:
        return True

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

def SFS_from_antr(table, chromosome, start, end, neutral_regions, min_alleles = None, neutral_only = False, measure = 'theta_pi'):
    SFSs = {}
    p = antr.Reader(table)
    for record in tqdm(p.fetch(chromosome, start, end)):
        # diversity calc
        allele_counts = record.quebec_alleles
        if neutral_only and not neutral_regions: # default - all three of intergenic, intronic, and fold 4
            if True not in [record.is_intronic, record.is_intergenic, record.is_fold4]:
                continue
        elif neutral_only and len(neutral_regions) >= 1:
            if test_neutral(record, neutral_regions): # whatever the specified neutral regions were
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
    if measure == 'theta_pi':
        diversity = sum([sfs.theta_pi() * sfs.sites() for sfs in SFSs.values()]) / sum([sfs.sites() for sfs in SFSs.values()])
        return diversity
    elif measure == 'theta_w':
        diversity = sum([sfs.theta_w() * sfs.sites() for sfs in SFSs.values()]) / sum([sfs.sites() for sfs in SFSs.values()])
        return diversity
    elif measure == 'both':
        try:
            diversity_tajima = sum([sfs.theta_pi() * sfs.sites() for sfs in SFSs.values()]) / sum([sfs.sites() for sfs in SFSs.values()])
        except ZeroDivisionError:
            diversity_tajima = 0
        try:
            diversity_watterson = sum([sfs.theta_w() * sfs.sites() for sfs in SFSs.values()]) / sum([sfs.sites() for sfs in SFSs.values()])
        except ZeroDivisionError:
            diversity_watterson = 0
        return diversity_tajima, diversity_watterson

def calculate_gene_density(table, chromosome, start, end):
    p = antr.Reader(table)
    genic_annotations = ['CDS', 'intronic', 'utr5', 'utr3']
    counts = OrderedDict.fromkeys(genic_annotations, 0)
    total_count = 0
    for record in p.fetch(chromosome, start, end):
        for annotation in genic_annotations:
            if attr_fetch(record, annotation):
                counts[annotation] += 1
                total_count += 1
    return counts, total_count

def eprint(*args, **kwargs): # throws out progress info
    print(*args, file=sys.stderr, **kwargs)

def main(table, windowsize, min_alleles, neutral_only, gene_density, measure, neutral_regions):
    if neutral_regions:
        eprint('Regions selected - {}'.format(neutral_regions))
    if measure == 'both':
        div_colname = 'theta_pi theta_w'
    else:
        div_colname = measure
    if gene_density:
        print('chromosome', 'start', 'end', '{}'.format(div_colname), 'rho', 'rho_total', 'rho_count', 'iter_count',
              'CDS_count', 'intronic_count', 'utr5_count', 'utr3_count', 'total_gene_count')
    else:
        print('chromosome', 'start', 'end', '{}'.format(div_colname), 'rho', 'rho_total', 'rho_count', 'iter_count')

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
                if not measure == 'both':
                    curr_div = SFS_from_antr(table, current_chrom, window[0], window[1], 
                                             min_alleles = min_alleles, neutral_only = neutral_only, 
                                             measure = measure, neutral_regions = neutral_regions)
                elif measure == 'both':
                    theta_pi, theta_w = SFS_from_antr(table, current_chrom, window[0], window[1],
                                                      min_alleles = min_alleles, neutral_only = neutral_only, 
                                                      measure = measure, neutral_regions = neutral_regions)
                    curr_div = ' '.join([str(theta) for theta in [theta_pi, theta_w]])
                if gene_density:
                    gene_counts, total_gene_count = calculate_gene_density(table, current_chrom, window[0], window[1])
            except ZeroDivisionError: # nothing in window
                rho_out = 0
                if not measure == 'both':
                    curr_div = 0
                elif measure == 'both':
                    curr_div = ' '.join(['0', '0'])

            if gene_density:
                print(current_chrom, window[0], window[1], curr_div, rho_out, rho, count, record_counter,
                      gene_counts['CDS'], gene_counts['intronic'], gene_counts['utr5'], gene_counts['utr3'], total_gene_count)
            elif not gene_density:
                print(current_chrom, window[0], window[1], curr_div, rho_out, rho, count, record_counter)
                
if __name__ == '__main__':
    arguments = args()
    main(*arguments)
