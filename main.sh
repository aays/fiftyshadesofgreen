for i in {1..17}; do
    ../ldhelmet find_confs -w 50 --num_threads 30 -i out$i\.conf ../../data/fastas/no_clones/chromosome_$i\.fasta
    sleep 1
    ../ldhelmet table_gen --num_threads 30 -o out$i\.lk -t 0.01 -r 0.0 0.1 10.0 1.0 100.0 -c out$i\.conf
    sleep 1
    ../ldhelmet pade --num_threads 30 -c out$i\.conf -t 0.01 -x 11 -i out$i\.pade
    sleep 1
    ../ldhelmet rjmcmc --num_threads 30 -w 50 -l out$i\.lk -p out$i\.pade -b 50 -s \
    ../../data/fastas/no_clones/chromosome_$i\.fasta --burn_in 100000 -n 1000000 -o out$i\.post
    sleep 1
    rm out$i\.lk
    ../ldhelmet post_to_text --mean -p 0.025 -p 0.975 -o out$i\.txt out$i\.post
    echo done $i;
done

cd ../../ # back to project dir

python3.5 analysis/fiftyshadesofgreen/annotation_parser/add_ld_rho.py > data/annotation_table_rho.txt # create new annotation table
echo 'made table'
bgzip data/annotation_table_rho.txt
tabix -p vcf data/annotation_table_rho.txt.gz
echo 'complete'

# GC content
python3.5 analysis/fiftyshadesofgreen/annotation_parser/antr_correlate_dict.py \
-t data/annotation_table_rho.txt.gz \
-w 50000 \
--gc_content > gc_content_50k.txt

# hotspots - using find_hotspots.py from Singhal 2015
for i in {1..17}; do 
    python3.5 find_hotspots.py \
    --input analysis/ldh_test/out$i\.txt \
    --out data/ldhelmet_files/hotspots/b2000_f1000000/chromosome_$i\.txt \
    --chr chromosome_$i \
    --block 2000 \
    --flank 1000000; 
done

# 20 kb windows
for i in {1..17}; do 
    python3.5 find_hotspots.py \
    --input analysis/ldh_test/out$i\.txt \
    --out data/ldhelmet_files/hotspots/b20k_f10mil/chromosome_$i\.txt \
    --chr chromosome_$i \
    --block 20000 \
    --flank 1000000; 
done

# 50 kb windows
for i in {1..17}; do 
    python3.5 find_hotspots.py \
    --input analysis/ldh_test/out$i\.txt \
    --out data/ldhelmet_files/hotspots/b50k_f10mil/chromosome_$i\.txt \
    --chr chromosome_$i \
    --block 50000 \
    --flank 1000000; 
done

# split annotation table into chromosomes
cd data
cp -v annotation_table_rho.txt.gz annotation_table_temp.txt.gz

bgzip -d annotation_table_temp.txt.gz

grep '#' annotation_table_temp.txt > annotation_header.txt

touch chroms.txt
for i in {1..17}; echo "chromosome_$i[^0-9]" >> chroms.txt; done

while read -r line; do
    grep "$line" annotation_table_temp.txt >> $line\.txt
done < chroms.txt
# for some reason, this was stuck on chr1 after copying in all lines - probably iterating through all lines
# I deleted the chromosome_1 regex pattern from chroms.txt and then it was fine again

for i in {1..17}; do
    touch chromosome_$i\.txt
    cat annotation_header.txt >> chromosome_$i\.txt
    cat chromosome_$i[^0-9].txt >> chromosome_$i\.txt
done

rm *0-9].txt
rm annotation_table_temp.txt
rm annotation_header.txt
rm chroms.txt

for i in {1..17}; do
    echo "zipping + tabixing chromosome $i..."
    bgzip chromosome_$i\.txt
    tabix -p vcf chromosome_$i\.txt.gz
done

cd ../

# correlates - run in parallel over each chromosome
parallel -j 17 -i sh -c 'python3.5 analysis/fiftyshadesofgreen/annotation_parser/antr_correlate_dict_chrom.py \
--table data/annotation_tables/chromosome_{}.txt.gz \
--windowsize 100000 \
--correlates is_genic intronic exonic utr3 utr5 intergenic CDS \
--gene_context 2000 \
--chromosome chromosome_{} > correlates_context_chromosome_{}.txt' -- {1..17}

mv -v correlates_context_chromosome_* data/correlates/gen_context/
cd !$

# concatenate
touch gen_correlates_context_100k.txt
cat correlates_context_chromosome_1.txt > gen_correlates_context_100k.txt

for i in {2..17}; do
    tail -n +2 correlates_context_chromosome_$i\.txt >> gen_correlates_context_100k.txt;
done
