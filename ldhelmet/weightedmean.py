"""
Returns a weighted mean r value for the entirety of a given LDHelmet (Chan et al. 2012) output file.
Each r value is weighted by the range it was calculated over.

Usage:

python weightedmean.py <infile> >> <outfile>

"""

import sys
import pandas as pd

def colfixer(data):
    data.drop(df.columns[[3, 4, 5]], axis = 1, inplace = True) 
    data.columns = ['left_snp', 'right_snp', 'mean'] 

def valuegetter(df, i):
    value = df.iloc[i,1] - df.iloc[i,0] # get range length
    value = value * df.iloc[i,2]
    return value
  
  
# analysis
infile = str(sys.argv[-1]) 
df = pd.read_csv(infile, sep = ' ', 
                 skiprows = 2, header = 'infer')   

colfixer(df)  
filelength = df.iloc[-1,1] - df.iloc[0,0] # get range of file
dfmeans = [valuegetter(df, row) for row in list(range(df.shape[0]))] 

infile = infile[infile.find('chromosome'):] # isolate chromosome name
finalvalue = sum(dfmeans)/filelength

print(infile, finalvalue) 
