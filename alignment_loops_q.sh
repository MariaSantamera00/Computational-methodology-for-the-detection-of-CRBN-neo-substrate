#!/bin/bash

## This script searches for the beta-hairpin loop motif in the alphafold models of all proteins in the human proteome and aligns it with all protein loops in the validation set. It calculates a score of that alignment and aligns the whole protein with respect to the loop with the highest score.  

# Script designed to be sent to queue


# Activate conda environment
conda activate maria

# Root directory
root="/home/l061003/TFM_MariaSantamera/"
alphafold_files="/lrlhps/users/l001803/TMP/pdb_structure/" #ARGUMENT

# Files and directoriesi
pdb_files="${root}validation_set/input_pdb/"
names_pdb_files="${root}validation_set/test.pdbs"
names_alphafold_files="${root}validation_set/test_AF.pdbs"

# Results directory
if [ ! -d "${root}alignment_AF/" ]
then
        mkdir "${root}alignment_AF/"
fi

results="${root}alignment_AF/"


#Ligand (target) and receptor (CRBN)
pdb=("5fqd" "5hxb" "6h0f" "6h0g" "6uml" "7lps")

ap=$(sed -n "${SGE_TASK_ID}p" $names_alphafold_files)
iap=$(basename "$ap" | cut -f 1 -d '.')

alphafold="${alphafold_files}${iap}.pdb"


## STEP 1.  Create a library with validation set motifs
awk '$6 >= 35 && $6 <=43 {print $0}' ${pdb_files}5fqd_target_a_l.pdb > ${pdb_files}loops/5fqd_target_loop.pdb
awk '$6 >= 569 && $6 <= 578 {print $0}' ${pdb_files}5hxb_target_a_l.pdb > ${pdb_files}loops/5hxb_target_loop.pdb
awk '$6 >= 145 && $6 <= 154 {print $0}' ${pdb_files}6h0f_target_a_l.pdb > ${pdb_files}loops/6h0f_target_loop.pdb
awk '$6 >= 417 && $6 <= 426 {print $0}' ${pdb_files}6h0g_target_a_l.pdb > ${pdb_files}loops/6h0g_target_loop.pdb
awk '$6 >= 410 && $6 <= 419 {print $0}' ${pdb_files}6uml_target_a_l.pdb > ${pdb_files}loops/6uml_target_loop.pdb
awk '$6 >= 140 && $6 <= 149 {print $0}' ${pdb_files}7lps_target_a_l.pdb > ${pdb_files}loops/7lps_target_loop.pdb


## STEP 2. GET THE LOOP  

#Get dssp file
mkdssp -i ${alphafold_files}${iap}.pdb -o ${alphafold_files}${iap}_dssp.txt
archivo_dssp="${alphafold_files}${iap}_dssp.txt"

#Write on a single line SS characters
cat ${archivo_dssp} | awk '/^  #/{p=1;next}p' | awk '{print substr($0, 17, 1)}' | awk '{printf "%s", $0}' > temp_1_${iap}.txt

#Look for the degrons loop pattern in the previous file
cat temp_1_${iap}.txt | grep -o -E "E+[TS ]{5,7}E+" > temp_2_${iap}.txt


count=0
while IFS= read -r patron
do
        #Look for the position of the motifs found in the sequence
        init=$(grep -bo "$patron" temp_1_${iap}.txt | awk -F ":" '{print $1 +1}')
        end=$(grep -bo "$patron" temp_1_${iap}.txt | awk -F ":" '{print $1 + length($2)-1+1}')

	#Look for the glycine
        pattern=$(cat ${archivo_dssp} | awk '/^  #/{p=1;next}p' | sed -n "${init},${end}p" | awk 'substr($0, 17, 1) ~ /[TS ]/{print $0}' | awk 'NR>=4 && NR<=6 && substr($0, 14, 1) == "G" {print $0}' | wc -l)
	
	if [ "$pattern" -ne 0 ]
        then
		count=$((count+1))
		st=$(cat ${archivo_dssp} | awk '/^  #/{p=1;next}p' | sed -n "${init},${end}p" | awk -v found=0 '{if ($5 ~ /[TS ]/ && found==0) {found=1; ini=NR-3}}END{print ini}')
		en=$(cat ${archivo_dssp} | awk '/^  #/{p=1;next}p' | sed -n "${init},${end}p" | awk '{if ($5 ~ /[TS ]/) {end=NR+4}}END{print end}')
		res1=$(cat ${archivo_dssp} | awk '/^  #/{p=1;next}p' | sed -n "${init},${end}p" | sed -n "${st},${en}p" | awk '{print substr($0, 8, 3)}' | sed 's/\s\+\([0-9]\)/\1/g' | head -n1)
		res2=$(cat ${archivo_dssp} | awk '/^  #/{p=1;next}p' | sed -n "${init},${end}p" | sed -n "${st},${en}p" | awk '{print substr($0, 8, 3)}' | sed 's/\s\+\([0-9]\)/\1/g' | tail -n1)
		if [ ! -d "${results}loops/" ]
		then
        		mkdir "${results}loops/"
		fi
		cat ${alphafold_files}${iap}.pdb | awk -v res1=$res1 -v res2=$res2 '($1=="ATOM" && $6>=res1 && $6<=res2) {print $0}' >  "${results}loops/${iap}_${count}_loop.pdb"



## STEP 3. ALIGNMENT SCORE between the found motif and the loops of the validation set

elec=-1
for i in ${pdb[@]}
do
	score=$(python ${root}bin/score_alignment.py ${pdb_files}loops/${i}_target_loop.pdb ${results}loops/${iap}_${count}_loop.pdb)
	if (( $(echo "$score > $elec" |bc -l) ))
	then 
		elec=$score
		loop=$i
	fi
done


## STEP 4. ALIGNMENT BETWEEN MOTIF AND BEST VALIDATION SET LOOP

if [ ! -d "${results}loops_al/" ]
then
	mkdir "${results}loops_al/"
fi

echo "load ${results}loops/${iap}_${count}_loop.pdb; load ${pdb_files}loops/${loop}_target_loop.pdb; super ${iap}_${count}_loop, ${loop}_target_loop ; save ${results}loops_al/${iap}_${count}_loop_al.pdb, ${iap}_${count}_loop" > pymol_${iap}_${count}.pml

pymol -c pymol_${iap}_${count}.pml >PYMOL.log 2>PYMOL_error.out
rm pymol_${iap}_${count}.pml

## STEP 5. ALIGNMENT BETWEEN AF PROTEIN AND ITS LOOP
if [ ! -d "${results}prot_al/" ]
then
        mkdir "${results}prot_al/"
fi

	./TMalign ${alphafold_files}${iap}.pdb ${results}loops_al/${iap}_${count}_loop_al.pdb -o ${results}prot_al/${iap}_${count}_al


fi

done < temp_2_${iap}.txt

rm temp_1_${iap}.txt
rm temp_2_${iap}.txt





