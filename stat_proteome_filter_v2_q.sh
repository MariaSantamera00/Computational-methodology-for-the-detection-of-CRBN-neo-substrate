#!/bin/bash

## This script searches among all the given proteins (the whole human proteome) for the motif: "Beta-hairpin loop, 5-7 aa, Gly in position 4, 5 or 6 of the loop". This script creates one file for each protein passing the filter.

# Script designed to be sent to queue


# Activate conda environment
#conda activate maria

# Function to print help in terminal
print_help(){
        echo "stat_proteome_filter_v2_q.sh"
        echo "María Santamera Lastras 2023"
        echo -e "\nusage:stat_proteome_filter_v2_q.sh <root_folder> <AF_files_folder> <names_AF_files> <output_folder>\n"
        echo -e "\troot_folder : root folder (should contain the GitHub validation set)"
        echo -e "\tAF_files_folder : folder in which AlphaFold files are stored"
        echo -e "\tnames_AF_files: file with AlphaFold files names"
        echo -e "\toutput_folder : folder in which data and results will be stored"
}

# Argument assignation
root=$1 
pdb_files=$2
names_pdb_files=$3 #names_pdb_files="${root}validation_set/test_AF.pdbs"
results=$4 #results="${root}filter1_AF/"


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


name=$(sed -n "${SGE_TASK_ID}p" $names_pdb_files)
i=$(basename "$name" | cut -f 1 -d '.')

#Get dssp file
mkdssp -i ${pdb_files}${i}.pdb -o ${pdb_files}${i}_dssp.txt
archivo_dssp="${pdb_files}${i}_dssp.txt"

#Write on a single line SS characters
cat ${archivo_dssp} | awk '/^  #/{p=1;next}p' | awk '{print substr($0, 17, 1)}' | awk '{printf "%s", $0}' > temp_1_${i}.txt

#Look for the degrons loop pattern in the previous file
cat temp_1_${i}.txt | grep -o -E "E+[TS ]{5,7}E+" > temp_2_${i}.txt

count=0
while IFS= read -r patron
do
	#rm "${results}proteins_filter.txt"
	#Look for the position of the motifs found in the sequence.
	init=$(grep -bo "$patron" temp_1_${i}.txt | awk -F ":" '{print $1 +1}')
	end=$(grep -bo "$patron" temp_1_${i}.txt | awk -F ":" '{print $1 + length($2)-1+1}')
		
	#Look for the glycine
	pattern=$(cat ${archivo_dssp} | awk '/^  #/{p=1;next}p' | sed -n "${init},${end}p" | awk 'substr($0, 17, 1) ~ /[TS ]/{print $0}' | awk 'NR>=4 && NR<=6 && substr($0, 14, 1) == "G" {print $0}' | wc -l) 

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


