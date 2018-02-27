import sys
import subprocess

"""
This quick hardcoded script exists to be run in parallel in order to generate variant files
from a .vcf file containing Chlamydomonas variants. The bash command would be

parallel -j 17 -i sh -c "python /scratch/research/projects/chlamydomonas/genomewide_recombination/analysis/vcfconvert.py chromosome_{}" -- {1..17}

and it would generate one fasta per chromosome containing variants. These can
later be concatenated with reference chromosome fastas using bash operations.

"""

lengths = {'chromosome_1': 8033585,
'chromosome_2': 9223677,
'chromosome_3': 9219486,
'chromosome_4': 4091191,
'chromosome_5': 3500558,
'chromosome_6': 9023763,
'chromosome_7': 6421821,
'chromosome_8': 5033832,
'chromosome_9': 7956127,
'chromosome_10': 6576019,
'chromosome_11': 3826814,
'chromosome_12': 9730733,
'chromosome_13': 5206065,
'chromosome_14': 4157777,
'chromosome_15': 1922860,
'chromosome_16': 7783580,
'chromosome_17': 7188315}

chromosome = str(sys.argv[-1])

script = 'time /scratch/research/repos/vcf2fasta/vcf2fasta.py '
vcf = '-v /scratch/research/projects/chlamydomonas/genomewide_recombination/data/vcfs/all_quebec.HC.vcf.gz -r /scratch/research/references/chlamydomonas/5.3_chlamy_w_organelles_mt_minus/chlamy.5.3.w_organelles_mtMinus.fasta '
region = '-i {chrname}:'.format(chrname = chromosome)
bases = '1-{length} '.format(length = lengths[chromosome])
gq = '-g 30 '
samples = '-s CC3071 CC2937 CC2936 CC3076 CC3086 CC3064 CC3068 CC3065 CC3060 CC3059 CC3062 CC2935 CC2938 CC3079 CC3075 CC3084 CC3073 CC3061 CC3063 '
outfile = '> /scratch/research/projects/chlamydomonas/genomewide_recombination/data/fastas/no_clones/{chrname}.fasta'.format(chrname = chromosome)

cmd = script + vcf + region + bases + gq + samples + outfile
print(cmd)
subprocess.call(cmd, shell = True)
