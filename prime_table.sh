#!/bin/bash

## This script generates a table with the energy obtained after the minimization performed with prime (obtained from the log file derived from the minimization).


# Activate conda environment
#conda activate maria

# Function to print help in terminal
print_help(){
        echo "prime_table.sh"
        echo "María Santamera Lastras 2023"
        echo -e "\nusage:prime_table.sh <root_folder> <prime_log_folder>\n"
        echo -e "\troot_folder : root folder (should contain the GitHub validation set)"
        echo -e "\tprime_log_folder : folder in which prime log files (with energy values) are stored. This is also the folder in which results will be stored"
}


# Argument assignation
root=$1 #root="/home/l061003/TFM_MariaSantamera/"
output_dir=$2 #output_dir="/lrlhps/users/l001803/TMP/crbn_eval/results_val/"

#Controls help message (print "print_help" function if user writes "-h" or "-help" in terminal)
if [ "$root" == "-h" ] || [ "$root" == "-help" ] ;then
        print_help
        exit
fi


#Control of arguments (show an error message if user doesn't write enough arguments)
if [ "$#" -ne 2 ]; then
    echo -e "ERROR: too few arguments\n"
    print_help
    exit
fi


# Files and directories
input_dir="${root}validation_set/input_pdb/"
names_pdb_files="${root}validation_set/test_cross.pdbs"


comb=$(cat ${names_pdb_files} | wc -l)
for num in $(seq 1 $comb)
do

	ap=$(sed -n "${num}p" $names_pdb_files)
	iap=$(echo "${ap}" | cut -f2 -d "-" | cut -f1 -d '.')
	crbn=$(echo "${ap}" | cut -f1 -d "-" | cut -f1 -d '.')

	prime_log="${output_dir}prime_output/ps_${iap}_${crbn}.log"
	energy=$(cat ${prime_log} | tail | grep "TOTALE"| awk '{print $2}')

	echo -e "${iap}-${crbn}\t${energy}"

done > "${output_dir}total_energy_cross_ligres.txt"

