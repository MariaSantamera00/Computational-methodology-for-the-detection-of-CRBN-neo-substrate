#!/bin/bash

## This script helps to statistically process the results obtained with the scripts "stat_proteome_filter.sh"

# Activate conda environment
conda activate maria

# Root directory
root="/home/l061003/TFM_MariaSantamera/"
filter_results_dir="" #ARGUMENT

# Create a single file with the names of filter proteins
cd ${filter_results_dir}
cat aa_AF-*.txt > aa_results.txt
sed -i '/^[[:space:]]*$/d' aa_results.txt

# Get the average length of the loops: 
awk '{ sum += length($0) } END { if (NR > 0) print sum/NR }' aa_results.txt

# Get the frequence of each aminoacid: 
car=$(cat aa_results.txt | wc -c)
awk -v var1="$car" '{ split($0, chars, ""); for (i = 1; i <= length($0); i++)
char_count[chars[i]]++ } END { for (char in char_count) print char, (char_count[char]/var1) }' aa_results.txt | sort


# Position of Gly in each loop
rm pos_gly.txt
while IFS= read -r line
do

echo "$line" | awk -v string="G" 'BEGIN{ SLENGTH = length(string) }
{
    skipped = 0
    starts = ""
    while ( SSTART = index($0,string) ) {
        starts = starts (starts?" ":"") (skipped + SSTART)
        $0 = substr($0,SSTART + SLENGTH)
        skipped += (SSTART + SLENGTH - 1)
    }
}
starts {print starts }' >> pos_gly.txt

done < aa_results.txt

# Put same length loops in the same file

count=0
while IFS= read -r line
do
	count=$((count+1))
	len="${#line}"	

		cat pos_gly.txt | sed -n "${count}p" >> pos_gly_${len}.txt
done < aa_results.txt

for j in $(seq 5 7)
do
	echo "Loop length: $j"
	for i in $(seq 4 6)
	do
		num=$(cat pos_gly_${j}.txt | grep "$i" | wc -l)
		tot=$(cat pos_gly_${j}.txt | wc -l)
		por=$(printf "%.4f\n" $((10**4 * num/tot))e-4)
		echo "Gly in pos $i : $num ($por)" 
	done
done


