# provide chromosome name as positional argument - ie chromosome_5
# usage: 
# chmod +x decay_pipeline.sh
# ./decay_pipeline.sh chromosome_15

java -jar analysis/GenomeAnalysisTK.jar -T GenotypeGVCFs \
-R /scratch/research/references/chlamydomonas/5.3_chlamy_w_organelles_mt_minus/chlamy.5.3.w_organelles_mtMinus.fasta \
-L $1 \
--variant data/gvcfs/CC2935.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC2938.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3079.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3075.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3084.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3073.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3061.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3063.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3071.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC2937.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC2936.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3076.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3086.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3064.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3068.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3065.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3060.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3059.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3062.haplotypeCalled.g.vcf.gz \
-o data/vcfs/no_clones/$1.vcf

bgzip data/vcfs/no_clones/$1.vcf
tabix data/vcfs/no_clones/$1.vcf.gz

# make sampled version in /no_clones/sampled
python3.5 analysis/fiftyshadesofgreen/vcf/snpsample2vcf.py data/vcfs/no_clones/$1.vcf.gz 100000 $1 data/vcfs/no_clones/sampled/$1

bgzip data/vcfs/no_clones/sampled/$1\sampled.vcf
tabix data/vcfs/no_clones/sampled/$1\sampled.vcf.gz

./analysis/plink \
--vcf data/vcfs/no_clones/sampled/$1\sampled.vcf.gz \
--make-bed \
--out data/vcfs/no_clones/temp_plink_files/$1\sampled \
--allow-extra-chr 

./analysis/plink \
--bfile data/vcfs/no_clones/temp_plink_files/$1\sampled \
--r2 \
--out data/vcfs/no_clones/plink_out/$1 \
--ld-window-r2 0 \
--ld-window-kb 100 \
--allow-extra-chr

Rscript analysis/fiftyshadesofgreen/plink/r2corrector.R data/vcfs/no_clones/plink_out/$1\.ld

python3.5 analysis/fiftyshadesofgreen/plink/zerocorrector.py \
data/vcfs/no_clones/plink_out/$1\.ldz > data/vcfs/no_clones/plink_out/$1\.txt

rm data/vcfs/no_clones/plink_out/$1\.ldz
