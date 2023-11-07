#!/bin/bash

## This script searches for the binding motif of the target to the CRBN pocket (beta-hairpin loop) and generates a ligand rest file that can be used in other experiments (such as docking with Lightdock).

# Activate conda environment
conda activate maria

# Root directory
root="/home/l061003/TFM_MariaSantamera/"

# Files and directories
pdb_files="${root}validation_set/input_pdb/"

# Results directory

if [ ! -d "${root}ss_motifs/" ]
then
        mkdir "${root}ss_motifs/"
fi

results="${root}ss_motifs"


list=("5fqd" "5hxb" "6h0f" "6h0g" "6uml" "7lps")


for i in ${list[@]}
do
# Get dssp file
mkdssp -i ${pdb_files}${i}_target_a.pdb -o ${pdb_files}${i}_dssp.txt
archivo_dssp="${pdb_files}${i}_dssp.txt"

# Write on a single line SS characters
cat ${archivo_dssp} | awk '/^  #/{p=1;next}p' | awk '{print substr($0, 17, 1)}' | awk '{printf "%s", $0}' > temp_1.txt

#Look for the degrons loop pattern in the previous file
cat temp_1.txt | grep -o -E "E+[TS ]{4,8}E+" > temp_2.txt

count=0

while IFS= read -r patron
do
	count=$((count+1))
	
	#Look for the position of the motifs found in the sequence.
	init=$(grep -bo "$patron" temp_1.txt | awk -F ":" '{print $1 +1}')
	end=$(grep -bo "$patron" temp_1.txt | awk -F ":" '{print $1 + length($2)-1+1}')
	
	#Create a ligand restrictions file if there is a G in the motif
	cat ${archivo_dssp} | awk '/^  #/{p=1;next}p' | sed -n "${init},${end}p" > temp.txt
	if grep -q -E "^.{13}G" temp.txt && grep -q -E "^.{16}[TS ]" temp.txt
	then 
		awk '{print "L L." substr($0, 14, 1) "." substr($0, 8, 3)}' temp.txt > "${results}/rest_lig_${count}_${i}.txt"
	fi

	#Replace G with Gly (all aminoacids) in the restrictions file
	while IFS= read -r line
	do
		aa1=$(echo $line | cut -f2 -d ".")
                aa3=$(python ${root}bin/transl_python.py "$aa1")
                echo $line | awk -F "." -v OFS='.' -v var1="$aa1" -v var2="$aa3" '{$2 = ( $2 ==var1 ? var2 : $2 ); print $0}'

	done < "${results}/rest_lig_${count}_${i}.txt"
		
done < temp_2.txt


rm temp.txt
rm temp_1.txt
rm temp_2.txt

done

