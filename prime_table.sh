#!/bin/bash

## This script generates a table with the energy obtained after the minimization performed with prime (obtained from the log file derived from the minimization).


# Activate conda environment
conda activate maria

# Root directory
root="/home/l061003/TFM_MariaSantamera/"

# Files and directoriesi
input_dir="${root}validation_set/input_pdb/"
names_pdb_files="${root}validation_set/test_cross.pdbs"
output_dir="/lrlhps/users/l001803/TMP/crbn_eval/results_val/" #ARGUMENT (prime results)


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

