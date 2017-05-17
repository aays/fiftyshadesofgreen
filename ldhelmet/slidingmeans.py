#!usr/bin/env python3.5

"""
Takes in an LDHelmet (Chan et al., 2012) .txt output file as well as a window size. Returns mean
values of r within each window.

usage:
python3.5 slidingmeans.py [ldhelmet txt infile] [windowsize] > [outfile]

This can be done across blocks with a nested for loop, if the user is looking to compare the effects of block
penalty changes at different window sizes:

for w in 500 1000 5000 10000 20000 30000; do
    for b in 10 50 100; do
        python3.5 slidingmeans.py chromosome_12_$b\.txt $w\ >> chr12slide$w\.txt ;
        done;
    done


AH - 03/2017
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

# analysis
df = pd.read_csv(sys.argv[1], sep = ' ',
                 skiprows = 2, header = 'infer')

colfixer(df)

outname = str(sys.argv[1])
outname = outname[outname.find('chromosome'):-4]
block = outname[outname.rfind('_') + 1:]

# creates machine-readable space-separated structure
slider(df, inputsize, block)



