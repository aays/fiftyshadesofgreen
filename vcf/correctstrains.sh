# full organelle linkage pipeline

touch minus.intervals
for line in mtDNA mtMinus:286502-287123 mtMinus:308874-311600;
    do echo $line >> minus.intervals;
done

touch plus.intervals
for line in chromosome_6:480316-484847 chromosome_6:559654-561075 cpDNA;
    do echo $line >> plus.intervals;
done

java -jar analysis/GenomeAnalysisTK.jar -T GenotypeGVCFs \
-R /scratch/research/references/chlamydomonas/5.3_chlamy_w_organelles_mt_minus/chlamy.5.3.w_organelles_mtMinus.fasta \
-L minus.intervals \
--variant data/gvcfs/CC2935.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC2938.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3079.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3075.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3084.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3073.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3061.haplotypeCalled.g.vcf.gz \
--variant data/gvcfs/CC3063.haplotypeCalled.g.vcf.gz \
-o mtmtd1midminus.vcf

java -jar analysis/GenomeAnalysisTK.jar -T GenotypeGVCFs \
-R /scratch/research/references/chlamydomonas/5.3_chlamy_w_organelles_mt_minus/chlamy.5.3.w_organelles_mtMinus.fasta \
-L plus.intervals \
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
-o cpmtafusplus.vcf

bgzip mtmtd1midminus.vcf
tabix mtmtd1midminus.vcf.gz

bgzip cpmtafusplus.vcf
tabix cpmtafusplus.vcf.gz    

rm plus.intervals
rm minus.intervals

mv -v mtmtd1midminus.vcf* data/vcfs/organellesredux/no_clones/correctstrains/
mv -v cpmtafusplus.vcf* data/vcfs/organellesredux/no_clones/correctstrains/

python3.5 analysis/fiftyshadesofgreen/vcf/singlevcfcalc.py \
-v data/vcfs/organellesredux/no_clones/correctstrains/mtmtd1midminus.vcf.gz \
-r mtMinus mtDNA \
-l d/dprime/r2 > data/organellelinkage/no_clones/correctstrains/mtmtd1mid.txt

python3.5 analysis/fiftyshadesofgreen/vcf/singlevcfcalc.py \
-v data/vcfs/organellesredux/no_clones/correctstrains/cpmtafusplus.vcf.gz \
-r chromosome_6 cpDNA \
-l d/dprime/r2 > data/organellelinkage/no_clones/correctstrains/cpmtafus.txt
