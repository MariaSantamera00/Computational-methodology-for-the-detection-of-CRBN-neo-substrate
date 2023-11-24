#!/bin/bash

## This script allows to generate the input file of the dpocket program for all given files.

# Function to print help in terminal
print_help(){
        echo "dpocket_proteome.sh"
        echo "María Santamera Lastras 2023"
        echo -e "\nusage:dpocket_proteome.sh <AF_files_folder> <Input_AF_files_names> <output_folder>\n"
        echo -e "\tAF_files_folder : folder in which AlphaFold files are stored"
        echo -e "\tInput_AF_files_names : file with the file names to be used to generate the output file"
        echo -e "\toutput_folder : folder in which data and results will be stored"
}


# Argument assignation
alphafold_files=$1 #alphafold_files="/lrlhps/users/l001803/TMP/crbn_eval/results/prime_output/"
names_pdb_files=$2 #names_pdb_files="/lrlhps/users/l001803/TMP/crbn_eval/results/temp_protein_filter_prime.txt"
results=$3 #results="/lrlhps/users/l001803/TMP/crbn_eval/results/"


# Controls help message (print "print_help" function if user writes "-h" or "-help" in terminal)
if [ "$alphafold_files" == "-h" ] || [ "$alphafold_files" == "-help" ] ;then
        print_help
        exit
fi


#Control of arguments (show an error message if user doesn't write enough arguments)
if [ "$#" -ne 3 ]; then
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



comb=$(cat ${names_pdb_files} | wc -l)
for num in $(seq 1 $comb)
do

name=$(sed -n "${num}p" $names_pdb_files)
i=$(echo "ps_${name}_al_5fqd_crbn_r_thalidomide-out")

echo -e "${alphafold_files}${i}.pdb\tEF2"


done > "${results}gpocket_input.txt"


