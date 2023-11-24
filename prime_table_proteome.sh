#!/bin/bash

#This script allows you to sort the results of interest from the "prime_sidechain_3a_proteome_q.sh" script into a table with protein IDs and total binding energy values. 

# Script designed to be sent to queue

# Change this field depending on the folder where the MMISMSA and structural alignment results are stored.
output_dir="/lrlhps/users/l001803/TMP/crbn_eval/results/"
input_dir="/lrlhps/users/l001803/TMP/prot_al/"
names_pdb_files="/lrlhps/users/l001803/TMP/test_filter_go.txt"


comb=$(cat ${names_pdb_files} | wc -l)

for num in $(seq 1 $comb)
do
	
	ap=$(sed -n "${num}p" $names_pdb_files)
	iap=$(echo "${ap}" | cut -f2 -d "-" | cut -f1 -d '.')
	crbn=$(echo "${ap}" | cut -f1 -d "-" | cut -f1 -d '.')

	energy=$(cat ${prime_log} | tail | grep "TOTALE"| awk '{print $2}')
	prime_log="${output_dir}prime_output/ps_${iap}_${crbn}.log"
	energy=$(cat ${prime_log} | tail | grep "TOTALE"| awk '{print $2}')

	echo -e "${iap}\t${energy}"

done > "${output_dir}total_energy_cross_ligres.txt"

