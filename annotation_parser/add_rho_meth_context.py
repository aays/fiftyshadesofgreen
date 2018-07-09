'''
add_ld_rho.py but it takes into account methylation

the methylation bed file has been split by chromosome in data/methylation/bed_split - should make it quicker to parse them

python3.5 add_rho_meth.py > [output table]

updated - this takes into account methylation context

methylation is therefore encoded as a list of size 2
this is because we'd otherwise have to add an extra column and modify the antm parser a whole lot more

AH - 07/2018
'''

import ant
from tqdm import tqdm

p = ant.Reader('data/annotation_table.txt.gz')
for item in p.header:
    print(item.strip())

print('##ld_rho=rate of recombination measured by LDhelmet on Quebec dataset. p/bp [FLOAT]')
print('##methylation=beta value and context in three clones of CC2937 [LIST]')

# add column names
p.cols.append('ld_rho')
p.cols.append('methylation')
p.cols = [item.strip() for item in p.cols]
print('#', '\t'.join(p.cols), sep = '') # add new cols

def makelookup(chrom):
    filename = 'data/methylation/bed_split_full/{}.bed'.format(chrom)
    lookup = {}
    with open(filename, 'r') as m:
        for line in tqdm(m):
            split = [i.rstrip() for i in line.split('\t')]
            c_pos, beta, context = int(split[1]), float(split[3]), str(split[5])
            lookup[c_pos] = [beta, context]
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
                        record += '\t' + str(rho) + '\t' + str(methylation_lookup[pos][0]) + str(methylation_lookup[pos][1])
                    else:
                        record = record + '\t' + str(rho) + '\t' + '[]'
                    print(record)
