'''
add_ld_rho.py but it takes into account methylation

the methylation bed file has been split by chromosome in data/methylation/bed_split - should make it quicker to parse them

python3.5 add_rho_meth.py > [output table]

AH - 02/2018
'''

import ant
from tqdm import tqdm

p = ant.Reader('data/annotation_table.txt.gz')
for item in p.header:
    print(item.strip())

print('##ld_rho=rate of recombination measured by LDhelmet on Quebec dataset. p/bp [FLOAT]')
print('##methylation=beta value at CpG sites in three clones of CC2937 [FLOAT]')

# add column names
p.cols.append('ld_rho')
p.cols.append('methylation')
p.cols = [item.strip() for item in p.cols]
print('#', '\t'.join(p.cols), sep = '') # add new cols

def check_methylation(chrom, pos):
    filename = 'data/methylation/bed_split/{}.bed'.format(chrom)
    with open(filename, 'r') as m:
        found = False
        for line in m:
            split = [i.rstrip() for i in line.split('\t')]
            chrom, c_pos, beta = str(split[0]), int(split[1]), float(split[3])
            if pos == c_pos:
                found = True
                return beta
            else:
                continue
        if not found:
            return None

def makelookup(chrom):
    filename = 'data/methylation/bed_split/{}.bed'.format(chrom)
    lookup = {}
    with open(filename, 'r') as m:
        for line in tqdm(m):
            split = [i.rstrip() for i in line.split('\t')]
            c_pos, beta = int(split[1]), float(split[3])
            lookup[c_pos] = beta
    return lookup

# write records
for i in range(1, 18):
    methylation_lookup = makelookup('chromosome_{}'.format(i))
    with open('analysis/ldh_test/out{}.txt'.format(i)) as f: # hardcoded for proj dir
        for line in tqdm(f):
            if line.startswith(('#', 'ver')):
                continue
            else:
                split = line.split(' ')
                start, end, rho = int(split[0]), int(split[1]), float(split[2])
                p = ant.Reader('data/annotation_table.txt.gz').fetch('chromosome_{}'.format(i), 
                                                                        start - 1, end - 1, raw = True)
                for record in p:
                    pos = int(record.split('\t')[1])
                    if pos in methylation_lookup.keys():
                        record = record + '\t' + str(rho) + '\t' + str(methylation_lookup[pos])
                    else:
                        record = record + '\t' + str(rho) + '\t' + '.'
                    print(record)
