'''
quick and dirty hardcoded-paths script to create a new annotation with LDhelmet rho values appended

AH - 11/2017
'''

import ant
from tqdm import tqdm

p = ant.Reader('data/annotation_table.txt.gz')
for item in p.header:
    print(item.strip())

print('##ld_rho=rate of recombination measured by LDhelmet on Quebec dataset. p/bp [FLOAT]')

# add column names
p.cols.append('ld_rho')
p.cols = [item.strip() for item in p.cols]
print('#', '\t'.join(p.cols), sep = '') # add ld_rho to end of colnames

# write records
for i in range(1, 18):
    with open('data/ldhelmet_files/finals/genome50/chromosome_{}_50.txt'.format(i)) as f: # hardcoded for proj dir
        for line in tqdm(f):
            if line.startswith(('#', 'ver')):
                continue
            else:
                split = line.split(' ')
                start, end, rho = int(split[0]), int(split[1]), float(split[2])
                p = ant.Reader('data/annotation_table.txt.gz').fetch('chromosome_{}'.format(i), 
                                                                        start - 1, end - 1, raw = True)
                for record in p:
                    record = record + '\t' + str(rho)
                    print(record)
