#!/bin/bash

## This script divides AlphaFold proteins into domains based on the DPAM database. Then, it searches for the beta-hairpin loop motif in these domains and aligns it with all protein loops in the validation set. It calculates a score of that alignment and aligns the whole protein domain with respect to the loop with the highest score.  

# Script designed to be sent to queue


# Activate conda environment
#conda activate maria

# Function to print help in terminal
print_help(){
	echo "alignment_loop_domains_q.sh"
	echo "María Santamera Lastras 2023"
	echo -e "\nusage:alignment_loop_domains_q.sh <root_folder> <AF_files_folder> <DPAM_files_folder> <output_folder>\n"
	echo -e "\troot_folder : root folder (should contain the GitHub validation set)"
	echo -e "\tAF_files_folder : folder in which AlphaFold files are stored"
	echo -e "\tDPAM_files_folder : folder in which DPAM files are stored"
	echo -e "\toutput_folder : folder in which data and results will be stored"
}



# Argument assignation 
root=$1 #root="/home/l061003/TFM_MariaSantamera/"
alphafold_files=$2 #alphafold_files="/lrlhps/users/l001803/TMP/pdb_structure/"
dpam_files=$3 #dpam_files="/home/l061003/Documents/DPAM/HomSa/" #ARGUMENT
results=$4 #results="${root}alignment_AF_domains/"

#Controls help message (print "print_help" function if user writes "-h" or "-help" in terminal)
if [ "$root" == "-h" ] || [ "$root" == "-help" ] ;then
        print_help
        exit
fi


#Control of arguments (show an error message if user doesn't write enough arguments)
if [ "$#" -ne 4 ]; then
    echo -e "ERROR: too few arguments\n"
    print_help
    exit
fi



# Files and directoriesi
pdb_files="${root}validation_set/input_pdb/"
names_pdb_files="${root}validation_set/test.pdbs"
names_alphafold_files="${root}validation_set/test_AF.pdbs"

# Results directory
if [ -d "${results}" ]
then
        echo "The folder exists"
	echo "Do you want to delete the previous results and recreate the directory? y/n"
	read answer
	
	# The directory will only be deleted and recreated if the user says "yes"
  	case $answer in
    		y)
      			rm -r ${results}
      			mkdir ${results}
    		;;
    		n)
      			exit
    		;;
    		*)
      			echo "You must choose y/n"
    		;;
  	esac
else
	mkdir ${results}
fi


#Ligand (target) and receptor (CRBN)
pdb=("5fqd" "5hxb" "6h0f" "6h0g" "6uml" "7lps")

ap=$(sed -n "${SGE_TASK_ID}p" $names_alphafold_files)
iap=$(basename "$ap" | cut -f 1 -d '.')

alphafold="${alphafold_files}${iap}.pdb"


## STEP 0.  Create a library with validation set motifs
awk '$6 >= 35 && $6 <=43 {print $0}' ${pdb_files}5fqd_target_a_l.pdb > ${pdb_files}loops/5fqd_target_loop.pdb
awk '$6 >= 569 && $6 <= 578 {print $0}' ${pdb_files}5hxb_target_a_l.pdb > ${pdb_files}loops/5hxb_target_loop.pdb
awk '$6 >= 145 && $6 <= 154 {print $0}' ${pdb_files}6h0f_target_a_l.pdb > ${pdb_files}loops/6h0f_target_loop.pdb
awk '$6 >= 417 && $6 <= 426 {print $0}' ${pdb_files}6h0g_target_a_l.pdb > ${pdb_files}loops/6h0g_target_loop.pdb
awk '$6 >= 410 && $6 <= 419 {print $0}' ${pdb_files}6uml_target_a_l.pdb > ${pdb_files}loops/6uml_target_loop.pdb
awk '$6 >= 140 && $6 <= 149 {print $0}' ${pdb_files}7lps_target_a_l.pdb > ${pdb_files}loops/7lps_target_loop.pdb


## STEP 1: DOMAIN DETECTION IN ALPHAFOLD STRUCTURES

uniprot_id=$(echo "${iap}" | cut -f2 -d "-")
dpam_file="${dpam_files}${uniprot_id}.domains"

if [ ! -d "${results}domains/" ]
then
        mkdir "${results}domains/"
fi


num_dom=0
while IFS= read -r line
do
        num_dom=$((num_dom+1))

        if [ -f "${results}domains/${iap}_d${num_dom}.pdb" ]
        then
                rm "${results}domains/${iap}_d${num_dom}.pdb"
        fi

        pos=$(echo "$line" | cut -f2)
        nfragments=$(echo "$pos" | tr "," "\n" | wc -l)

        for (( i=1; i<=$nfragments; i++ ))
        do
                fragment=$(echo "$pos" | tr "," "\n" | sed -n "${i}p")
                st=$(echo "$fragment" | cut -f1 -d "-")
                end=$(echo "$fragment" | cut -f2 -d "-")
                cat $alphafold | awk -v st=$st -v end=$end '$1 == "ATOM" && $6 >= st && $6 <= end {print $0}' >> "${results}domains/${iap}_d${num_dom}.pdb"
        done

done < "${dpam_file}"

dom_files="${results}domains/"



## STEP 2. GET THE LOOP  

for (( dom=1; dom<=$num_dom; dom++ ))
do


#Get dssp file
mkdssp -i ${dom_files}${iap}_d${dom}.pdb -o ${dom_files}${iap}_d${dom}_dssp.txt
archivo_dssp="${dom_files}${iap}_d${dom}_dssp.txt"

#Write on a single line SS characters
cat ${archivo_dssp} | awk '/^  #/{p=1;next}p' | awk '{print substr($0, 17, 1)}' | awk '{printf "%s", $0}' > temp_1_${iap}_${dom}.txt

#Look for the degrons loop pattern in the previous file
cat temp_1_${iap}_${dom}.txt | grep -o -E "E+[TS ]{5,7}E+" > temp_2_${iap}_${dom}.txt


count=0
while IFS= read -r patron
do
        #Look for the position of the motifs found in the sequence
        init=$(grep -bo "$patron" temp_1_${iap}_${dom}.txt | awk -F ":" '{print $1 +1}')
        end=$(grep -bo "$patron" temp_1_${iap}_${dom}.txt | awk -F ":" '{print $1 + length($2)-1+1}')

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

		cat ${dom_files}${iap}_d${dom}.pdb | awk -v res1=$res1 -v res2=$res2 '($1=="ATOM" && $6>=res1 && $6<=res2) {print $0}' > "${results}loops/${iap}_d${dom}_${count}_loop.pdb"


## STEP 3. ALIGNMENT SCORE between the found motif and the loops of the validation set

elec=-1
loop=""
for i in ${pdb[@]}
do
	score=$(python ${root}bin/score_alignment.py ${pdb_files}loops/${i}_target_loop.pdb ${results}loops/${iap}_d${dom}_${count}_loop.pdb)
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

echo "load ${results}loops/${iap}_d${dom}_${count}_loop.pdb; load ${pdb_files}loops/${loop}_target_loop.pdb; super ${iap}_d${dom}_${count}_loop, ${loop}_target_loop ; save ${results}loops_al/${iap}_d${dom}_${count}_loop_al.pdb, ${iap}_d${dom}_${count}_loop" > pymol_${iap}_d${dom}_${count}.pml

pymol -c pymol_${iap}_d${dom}_${count}.pml >PYMOL.log 2>PYMOL_error.out

rm pymol_${iap}_d${dom}_${count}.pml



## STEP 5. ALIGNMENT BETWEEN AF PROTEIN AND ITS LOOP
if [ ! -d "${results}prot_al/" ]
then
        mkdir "${results}prot_al/"
fi

	./TMalign ${dom_files}${iap}_d${dom}.pdb ${results}loops_al/${iap}_d${dom}_${count}_loop_al.pdb -o ${results}prot_al/${iap}_d${dom}_${count}_al


fi

done < temp_2_${iap}_d${dom}.txt

rm temp_1_${iap}_d${dom}.txt
rm temp_2_${iap}_d${dom}.txt


done


# Filtering (check if any protein domain has at least one loop)

f=$(ls ${results}prot_al/ | grep "${iap}" | grep -c "al.pdb")

if [ "$f" -ne 0 ]
then
	if [ ! -d "${results}filter/" ]
	then
        	mkdir "${results}filter/"
	fi

        echo "${iap}" > "${results}filter/sum_${iap}.txt"
fi





