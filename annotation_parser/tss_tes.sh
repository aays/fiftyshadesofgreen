for i in {1..17};
    do echo chrom $i ;
    python3.5 analysis/fiftyshadesofgreen/annotation_parser/tss_tes.py \
    --gff data/final.strict.GFF3 \
    --table data/annotation_table_rho.txt.gz \
    --distance 50000 \
    --windowsize 2000 \
    --chromosome chromosome_$i > data/correlates/overall/tss_tes/take2/windowed/d50k_w2k/tss_d50k_w2k_chromosome_$i\.txt;
done

cd data/correlates/overall/tss_tes/take2/windowed/d50k_w2k/
touch genomewide_d50k_w2k.txt
head -n +1 tss_d50k_w2k_chromosome_1.txt >> genomewide_d50k_w2k.txt

for i in {1..17};
    do tail -n +2 tss_d50k_w2k_chromosome_$i\.txt >> genomewide_d50k_w2k.txt;
done

Rscript ../../graphmaker_windowed.R genomewide_d50k_w2k.txt

prol
