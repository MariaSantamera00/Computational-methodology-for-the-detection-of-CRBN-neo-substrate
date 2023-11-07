#!/bin/bash

## This script searches among all the given proteins (the whole human proteome) for the motif: "Beta-hairpin loop, 4-8 aa, Gly in any position of the loop". This script creates one file for each protein passing the filter.

# Script designed to be sent to queue



# Activate conda environment
conda activate maria

# Root directory
root="/home/l061003/TFM_MariaSantamera/"
pdb_files="/lrlhps/users/l001803/TMP/pdb_structure/" #ARGUMENT


# Files and directoriesi
names_pdb_files="${root}validation_set/test_AF.pdbs"

# Results directory
if [ ! -d "${root}filter1_AF/" ]
then
        mkdir "${root}filter1_AF/"
fi

results="${root}filter1_AF/"


name=$(sed -n "${SGE_TASK_ID}p" $names_pdb_files)
i=$(basename "$name" | cut -f 1 -d '.')

# Get dssp file
mkdssp -i ${pdb_files}${i}.pdb -o ${pdb_files}${i}_dssp.txt
archivo_dssp="${pdb_files}${i}_dssp.txt"

# Write on a single line SS characters
cat ${archivo_dssp} | awk '/^  #/{p=1;next}p' | awk '{print substr($0, 17, 1)}' | awk '{printf "%s", $0}' > temp_1_${i}.txt

# Look for the degrons loop pattern in the previous file
cat temp_1_${i}.txt | grep -o -E "E+[TS ]{4,8}E+" > temp_2_${i}.txt

count=0
while IFS= read -r patron
do
	#rm "${results}proteins_filter.txt"
	# Look for the position of the motifs found in the sequence.
	init=$(grep -bo "$patron" temp_1_${i}.txt | awk -F ":" '{print $1 +1}')
	end=$(grep -bo "$patron" temp_1_${i}.txt | awk -F ":" '{print $1 + length($2)-1+1}')
		
	# Look for the Gly
	pattern=$(cat ${archivo_dssp} | awk '/^  #/{p=1;next}p' | sed -n "${init},${end}p" | awk 'substr($0, 14, 1) == "G" && substr($0, 17, 1) ~ /[TS ]/{print $0}' | wc -l) 

	if [ "$pattern" -ne 0 ]
	then 
		count=$((count+1))
		cat ${archivo_dssp} | awk '/^  #/{p=1;next}p' | sed -n "${init},${end}p" | awk 'substr($0, 17, 1) ~ /[TS ]/{print substr($0, 14, 1)}' | awk '{printf "%s", $0}'  > "${results}aa_${i}_${count}.txt"
		echo -e "\n" >> "${results}aa_${i}_${count}.txt"
	fi
done < temp_2_${i}.txt

	if [ "$count" -ne 0 ]
        then
                echo -e "$i" > "${results}proteins_filter_${i}.txt"
        fi


rm temp_${i}.txt
rm temp_1_${i}.txt
rm temp_2_${i}.txt


