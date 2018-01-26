import antr
from tqdm import tqdm
import sys

table = sys.argv[-1]

p = antr.Reader(table)
for item in p.header:
    print(item.strip())

print('##methylation=beta value at CpG sites in three clones of CC2937. [FLOAT]')

# add column names
p.cols.append('methylation')
p.cols = [item.strip() for item in p.cols]
print('#', '\t'.join(p.cols), sep = '') # add ld_rho to end of colnames

# write records
with open('data/methylation/beta_vals_no_context.bed', 'f') as f:
    for line in tqdm(f):
        split = [i.rstrip() for i in line.split('\t')]
        chrom, c_pos, beta = str(split[0]), int([split[1]]), float(split[3])

        p = antr.Reader(table).fetch(chrom, c_pos, c_pos + 1, raw = True)

        for record in p:
            record = record + '\t' + str(beta)
            print(record)
