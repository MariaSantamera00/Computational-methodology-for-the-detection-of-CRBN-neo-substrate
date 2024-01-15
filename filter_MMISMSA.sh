#!/bin/bash

# This script filters the total energy results of MMISMSA and allows to obtain a list of the filtered complexes. 

# Change this field depending on the folder where the MMISMSA and prime results are stored.
# dir=""



count=0
while IFS= read -r line
do

cd "${dir}phcal_3/${line}"
score=$(cat "mm_result_energy" | sed -n '2p' | cut -f6 -d ";")

if (( $(echo "$score < -30" |bc -l) && $(echo "$score > -35" |bc -l) ))
then
	count=$((count+1))
	echo -e "${dir}prime_output/${line}\tEF2" | sed 's/phcal_3//g' | sed 's/-out/-out.pdb/g'
fi

done < "${dir}ph_proteome_6756.txt" > "${dir}dpocket_input_6756_30-35.txt"

echo $count
