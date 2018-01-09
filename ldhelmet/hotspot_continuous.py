''' python3.5 hotspot_continuous.py [hotspot csv output (.txt/csv)] > [out (csv)]

labels continuous hotspots for downstream dplyr groupings (ie hotspots larger than 1 x window size)
'''

import sys

filename = sys.argv[-1]

with open(filename, 'r') as f:
    all_lines = f.readlines() # have to read into an object since iterating through two lines at a time

header = all_lines[0].rstrip() + ',hotspot_group'
print(header)
all_lines = [line.rstrip().split(',') for line in all_lines[1:]]

counter = 1
for i in range(len(all_lines) - 1):
    left = all_lines[i]
    right = all_lines[i + 1]
    if left[2] != right[1]:
        out = ','.join(all_lines[i]) + ',' + str(counter)
        print(out)
        counter += 1
    elif left[2] == right[1]: # continuous hotspot
        out = ','.join(all_lines[i]) + ',' + str(counter)
        print(out)

final_line = all_lines[-1]
print(','.join(final_line) + ',' + str(counter))
