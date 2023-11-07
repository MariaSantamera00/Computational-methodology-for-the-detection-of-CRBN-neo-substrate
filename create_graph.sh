#!/bin/bash

## This script uses the "distances.py" script to select the residues whose alpha carbons are less than 8 A and creates a graph with these connections (using the graph.py script). It then calculates the connectivity. 

# Activate conda environment
conda activate maria

# Root directory
root="/home/l061003/TFM_MariaSantamera/"

# Files and directories
pdb_dir="${root}validation_set/AF_proteins_examples/"
pdb_file="${pdb_dir}AF-Q14191-F1-model_v4.pdb"

# Results directory

if [ ! -d "${root}graph_results/" ]
then
        mkdir "${root}graph_results/"
fi

results="${root}graph_results/"

name=$(basename "$pdb_file" | cut -f 1 -d '.')

# Get pdb file only with CA
awk '($1=="ATOM" && $3=="CA"){print $0}' ${pdb_file} > "${results}${name}_ca.pdb"

# Get the input for python script 
# Coordinates
xyz=$(echo "[")
while IFS= read -r line
do
	xyz=$(echo "${xyz}[$(echo "$line" | awk '{print $7}'),$(echo "$line" | awk '{print $8}'),$(echo $line | awk '{print $9}')]")

done < "${results}${name}_ca.pdb"

xyz=$(echo "$(echo "$xyz" | sed 's/]\[/],[/g')]")


# Compute matrix distance. Get nodes and edges file. 
python ${root}bin/distances.py "${xyz}" > "${results}${name}_edges.txt"

# Create the graph
python ${root}bin/graph.py "${results}${name}_edges.txt" 


