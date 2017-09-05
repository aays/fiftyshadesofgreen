# attempt at all inter-chr pairs pipeline

java -jar analysis/GenomeAnalysisTK.jar -T GenotypeGVCFs -R /scratch/research/references/chlamydomonas/5.3_chlamy_w_organelles_mt_minus/chlamy.5.3.w_organelles_mtMinus.fasta -XL exclusions.intervals \
--variant data/gvcfs/CC3071.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC2937.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC2936.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3076.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3086.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3064.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3068.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3065.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3060.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC2935.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC2938.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3079.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3075.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3084.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3062.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3073.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3059.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3061.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3063.haplotypeCalled.g.vcf.gz \
-o allvariants.vcf

bgzip allvariants.vcf
tabix allvariants.vcf.gz
mv -v allvariants.vcf* data/vcfs/organellesredux/no_clones

date
echo 'called initial vcf'

python3.5 analysis/fiftyshadesofgreen/vcf/vcf_subset.py -v data/vcfs/organellesredux/no_clones/allvariants.vcf.gz -c chromosome_1 chromosome_2 chromosome_3 chromosome_4 chromosome_5 chromosome_6 chromosome_7 chromosome_8 chromosome_9 chromosome_10 chromosome_11 chromosome_12 chromosome_13 chromosome_14 chromosome_15 chromosome_16 chromosome_17 -f 0.0002 -o data/vcfs/organellesredux/no_clones/allvariantsfiltered.vcf

date
echo 'subset complete'

grep '#CHROM' data/vcfs/organellesredux/no_clones/allvariantsfiltered.vcf | head -n 1 > header
grep -v '#' data/vcfs/organellesredux/no_clones/allvariantsfiltered.vcf | sponge data/vcfs/organellesredux/no_clones/allvariantsfiltered.vcf
cat header data/vcfs/organellesredux/no_clones/allvariantsfiltered.vcf | sponge data/vcfs/organellesredux/no_clones/allvariantsfiltered.vcf
rm header

date
echo 'formatting corrected'

bgzip data/vcfs/organellesredux/no_clones/allvariantsfiltered.vcf
tabix data/vcfs/organellesredux/no_clones/allvariantsfiltered.vcf.gz

date

time python3.5 analysis/fiftyshadesofgreen/vcf/allpairscalc.py -v data/vcfs/organellesredux/no_clones/allvariantsfiltered.vcf.gz -l d/dprime/r2 --haps > data/organellelinkage/no_clones/allpairs.txt
