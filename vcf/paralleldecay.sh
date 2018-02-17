parallel -j 5 -i sh -c "python3.5 analysis/fiftyshadesofgreen/vcf/ld_decay.py --vcf data/vcfs/no_clones/$1\.vcf.gz --count 100000 > data/new_lddecay/$1\_{}.txt" -- {1..5}

touch data/new_lddecay/$1\_500k.txt

cd data/new_lddecay

for i in {1..5};
    do cat $1\_$i\.txt >> $1\_500k.txt;
    rm $1\_$i\.txt;
done

cd -

### 100k

parallel -j 5 -i sh -c "python3.5 analysis/fiftyshadesofgreen/vcf/ld_decay.py --vcf data/vcfs/no_clones/$1\.vcf.gz --count 100000 --filter 100000 > data/new_lddecay/$1\_{}_filt.txt" -- {1..5}

touch data/new_lddecay/$1\_500k_filt.txt

cd data/new_lddecay

for i in {1..5}; do
    cat $1\_$i\_filt.txt >> $1\_500k_filt.txt;
    rm $1\_$i\_filt.txt;
done

cd -
