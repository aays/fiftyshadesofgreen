"""
This script takes an ldhelmet .txt outfile as input
and returns a mean value for the whole thing. If an 
outfile name is specified, it will append the value 
to said outfile.

Usage:

python ldhelmetmeans.py <infile> >> <outfile>

"""

import sys
import pandas as pd

# arg setup
infile = str(sys.argv[-1]) 

# initializing df and means object
df = pd.read_csv(infile, sep = ' ', 
                 skiprows = 2, header = 'infer') 

dfmeans = []

# data wrangling
def colfixer(data):
    data.drop(df.columns[[3, 4, 5]], axis = 1, inplace = True) 
    data.columns = ['left_snp', 'right_snp', 'mean'] 
    
colfixer(df)

# calculating weighted means
for i in list(range(df.shape[0])):
    value = df.iloc[i,1] - df.iloc[i,0] # get range length
    value = value * df.iloc[i,2]
    dfmeans.append(value)     
    
# get range of file
filelength = df.iloc[-1,1] - df.iloc[0,0]

# final calc    
infile = infile[infile.find('chromosome'):] # isolate chromosome name

finalvalue = sum(dfmeans)/filelength

print(infile, finalvalue) 
