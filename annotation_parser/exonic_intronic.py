'''
computes rho in first/last/other exons and introns

usage:
python3.5 exonic_intronic.py --gff [.gff] --table [.txt.gz] > [output (.txt.)]

AH - 12/2017
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

parser.add_argument('-g', '--gff', required = True,
                   type = str, help = 'GFF file (.GFF/.GFF3)')
parser.add_argument('-t', '--table', required = True,
                   type = str, help = 'Annotation table file (.txt.gz)')

args = parser.parse_args()
table = args.table
gff_file = args.gff

def eprint(*args, **kwargs): # throws out progress info
    print(*args, file=sys.stderr, **kwargs)

# analysis

with open(gff_file) as f:
    gff = f.readlines()

gff = [line.split('\t') for line in gff]
gff = [line for line in gff if line[0] not in ['cpDNA', 'mtMinus', 'mtDNA'] and 'chromosome' in line[0]]

gff_filt = [line for line in gff if line[2] in ['gene', 'exon', 'intron']] # remove everything else

genes = dict.fromkeys(['chromosome_' + str(i) for i in range(1,18)], [])

for key in genes.keys():
    snippet = [line for line in gff_filt if line[0] == key and 'gene' in line[2]]
    genes[key] = [[int(line[3]), int(line[4])] for line in snippet]

# create dicts to store coords in
first_exons = dict.fromkeys(['chromosome_' + str(i) for i in range(1,18)], [])
other_exons = dict.fromkeys(['chromosome_' + str(i) for i in range(1,18)], [])
last_exons = dict.fromkeys(['chromosome_' + str(i) for i in range(1,18)], [])
first_introns = dict.fromkeys(['chromosome_' + str(i) for i in range(1,18)], [])
other_introns = dict.fromkeys(['chromosome_' + str(i) for i in range(1,18)], [])
last_introns = dict.fromkeys(['chromosome_' + str(i) for i in range(1,18)], [])

for key in genes.keys():
    for gene in tqdm(genes[key]):

        # pull all exons in current gene
        all_exons = [line for line in gff_filt if line[0] == key and 'exon' in line[2]]
        current_exons = [[int(line[3]), int(line[4])] for line in all_exons \
                        if int(line[3]) >= gene[0] and int(line[4]) <= gene[1]]

        # pull all introns in current gene
        all_introns = [line for line in gff_filt if line[0] == key and 'intron' in line[2]]
        current_introns = [[int(line[3]), int(line[4])] for line in all_introns \
                        if int(line[3]) >= gene[0] and int(line[4]) <= gene[1]]

        try:
            first_exons[key] = first_exons[key] + [[current_exons[0][0], current_exons[0][1]]]
            other_exons[key] = other_exons[key] + current_exons[1:-1]
            last_exons[key] = last_exons[key] + [[current_exons[-1][0], current_exons[-1][1]]]
        except IndexError:
            eprint('exon not found at', key, gene)
            
        try:
            first_introns[key] = first_introns[key] + [[current_introns[0][0], current_introns[0][1]]]
            other_introns[key] = other_introns[key] + current_introns[1:-1]
            last_introns[key] = last_introns[key] + [[current_introns[-1][0], current_introns[-1][1]]]
        except IndexError: # some genes don't have intronic regions?
            continue


# col headers for file
print('chromosome', 'start', 'end', 'length', 'type', 'order', 'rho', 'total_rho', 'count')

# iterate through chromosomes
for chrom in range(1, 18):
    eprint('starting', chrom)
    current_chrom = 'chromosome_{}'.format(str(chrom))
    eprint('current_chrom = ', current_chrom)
    p = antr.Reader(table)

    # exons
    exon_order = 0 # keep track of which dict for the 'order' column
    # 0 - first, 1 - other, 2 - last
    for coord_dict in [first_exons, other_exons, last_exons]:
        for exon in tqdm(coord_dict[current_chrom]):
            
            p = antr.Reader(table)

            exon = [int(v) for v in exon]

            exon_start, exon_end = exon # unpack
            exon_length = exon_end - exon_start
            exon_total_rho = 0.0
            count = 0

            # iterate through exon records
            for record in p.fetch(current_chrom, exon_start, exon_end):
                try:
                    assert record.is_exonic
                except AssertionError:
                    eprint(current_chrom, exon_start, exon_end, counter)
                    break
                exon_total_rho += record.ld_rho
                count += 1

            try:
                rho = exon_total_rho / count
            except ZeroDivisionError:
                rho = 0

            print(current_chrom, exon_start, exon_end, exon_length, 'e', exon_order, rho, exon_total_rho, count)

        exon_order += 1

    # introns
    intron_order = 0 # keep track of which dict for the 'order' column
    # 0 - first, 1 - other, 2 - last
    for coord_dict in [first_introns, other_introns, last_introns]:
        for intron in tqdm(coord_dict[current_chrom]):
            
            p = antr.Reader(table)

            intron = [int(v) for v in exon]

            intron_start, intron_end = exon # unpack
            intron_length = intron_end - intron_start
            intron_total_rho = 0.0
            count = 0

            # iterate through exon records
            for record in p.fetch(current_chrom, intron_start, intron_end):
                try:
                    assert record.is_exonic
                except AssertionError:
                    eprint(current_chrom, intron_start, intron_end, counter)
                    break
                intron_total_rho += record.ld_rho
                count += 1

            try:
                rho = intron_total_rho / count
            except ZeroDivisionError:
                rho = 0

            print(current_chrom, intron_start, intron_end, intron_length, 'i', intron_order, rho, intron_total_rho, count)

        intron_order += 1
