'''
GC_content.py - returns windowed GC content values from fasta files (overlapping windows)

usage: python3.5 GC_content.py -i [fasta] -w [windowsize] > [outfile]

useful for correlating GC content with LDhelmet rho outputs

AH - 10/2017
'''

import argparse
from Bio import SeqIO
from tqdm import tqdm

parser = argparse.ArgumentParser(description = 'Get GC content in specified windows from a fasta file.', 
                                usage = 'gc_content.py [options]')
parser.add_argument('-i', '--input', required = True,
                   type = str, help = 'Input fasta. (Required)')
parser.add_argument('-w', '--windowsize', required = True,
                   type = int, help = 'Window size. (Required)')
args = parser.parse_args()

# function(s)
def GC_content(chrname, seq, windowsize):
    windows = list(range(0, len(seq), windowsize))
    windows.append(len(seq))
    halfwindow = int(windowsize / 2)
    
    for i in tqdm(range(len(windows) - 1)):
        subseq = seq[windows[i]:windows[i + 1]]
        GC = subseq.count('G') + subseq.count('C')
        out1 = float(GC) / (windows[i + 1] - windows[i]) # divide GC count by total bases
        print(chrname, windows[i], windows[i + 1], out1)
        
        if i < len(windows) - 2:
            coord1 = windows[i] + halfwindow
            coord2 = windows[i + 1] + halfwindow
            subseq_ahead = seq[coord1:coord2]
            GC_ahead = subseq_ahead.count('G') + subseq.count('C')
            out2 = float(GC_ahead) / ((windows[i + 1] + halfwindow) - (windows[i] + halfwindow))
            print(chrname, windows[i] + halfwindow, windows[i + 1] + halfwindow, out2)

# analysis
infile = SeqIO.parse(args.input, 'fasta')

print('chr block_start block_end block_GC')
for rec in infile:
    GC_content(rec.id, rec.seq, args.windowsize)
