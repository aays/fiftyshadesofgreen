mkdir annotatin_tables
cp -v annotation_table_rho.txt.gz annotation_table_temp.txt.gz

bgzip -d annotation_table_temp.txt.gz

mv -v annotation_table_temp.txt annotation_tables/

cd annotation_tables/

grep '#' annotation_table_temp.txt > annotation_header.txt

touch chroms.txt # create regex patterns
for i in {1..17}; echo "chromosome_$i[^0-9]" >> chroms.txt; done

while read -r line; do
    grep "$line" annotation_table_temp.txt >> $line\.txt
done < chroms.txt

for i in {1..17}; do
    touch chromosome_$i\.txt
    cat annotation_header.txt >> chromosome_$i\.txt
    cat chromosome_$i[^0-9].txt >> chromosome_$i\.txt
done

# cleanup
rm *0-9].txt
rm annotation_table_temp.txt
rm annotation_header.txt
rm chroms.txt

for i in {1..17}; do
    echo "zipping + tabixing chromosome $i..."
    bgzip chromosome_$i\.txt
    tabix -p vcf chromosome_$i\.txt.gz
done
