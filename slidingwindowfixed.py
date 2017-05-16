#!usr/bin/env python3.5

"""
Takes in an ldhelmet .txt output file as well as
a window size. Returns means in that windowsize as
well as the mean of those means. Prints every tenth
window value to terminal for the sake of progress
monitoring.

usage:
python3.5 slidingwindowfixed.py [ldhelmet txt infile] [windowsize] > [outfile]

can be done across blocks with a nested for loop:

for w in 500 1000 5000 10000 20000 30000; do
    for b in 10 50 100; do
        python3.5 slidingwindowfixed.py chromosome_12_$b\.txt $w\ >> slide$w\.txt ;
        done;
    done


reference - notebook 7.1
"""

import pandas as pd
import sys

inputsize = int(sys.argv[2])

# functions
def colfixer(df):
    df.drop(df.columns[[3, 4, 5]], axis = 1, inplace = True)
    df.columns = ['left_snp', 'right_snp', 'mean']

def meaner(df, i):
    value = df['right_snp'][i] - df['left_snp'][i] # get range length
    value = value * df['mean'][i]
    return value

def slider(df, windowsize, block):
    dfmeans = []
    windowmeans = []
    counter = 0
    windowlist = list(range(int(df.iloc[0,0]), int(df.iloc[-1,1]), windowsize))
    for window in windowlist:
        counter = counter + 1
        winstart = window
        winend = window + windowsize
        subdf = df[(df['left_snp'] >= winstart) & (df['right_snp'] <= winend)].reset_index()
        currentwindowmean = [meaner(df, i) for i in range(subdf.shape[0])]
        if sum(currentwindowmean) == 0:
            print(winstart, 0, block)
            continue
        else:
            currentwindowmean = sum(currentwindowmean)/len(currentwindowmean)
            print(winstart, currentwindowmean, block)
            windowmeans.append(currentwindowmean)
#    return windowmeans

# analysis
df = pd.read_csv(sys.argv[1], sep = ' ',
                 skiprows = 2, header = 'infer')

colfixer(df)
# print('fixed columns')

outname = sys.argv[1]
outname = outname[outname.find('chromosome'):-4]
block = outname[outname.rfind('_') + 1:]

# creates machine-readable space-separated structure
slider(df, inputsize, block)
# print('done sliding')



