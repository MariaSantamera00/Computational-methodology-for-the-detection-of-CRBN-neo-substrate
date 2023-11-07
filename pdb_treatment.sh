#!/bin/bash

## This script allows to preprocess pdb files to obtain some of the files in the directory "validation_set". It is neccessary to prepare the files for subsequent steps.

# Root directory
root="/home/l061003/TFM_MariaSantamera/"


pdb_files="${root}validation_set/input_pdb"


# Remove the molecular glue
files_glues=$(find $pdb_files -type f | grep -E "*lig*")

for file in $files_glues
do
	echo "Removing molecular glue in $file..."
	name=$(basename "$file" | cut -f 1 -d '.')
	new_name=$(echo "$name" | sed 's/_lig_union/_crbn/')
	grep '^ATOM' $file > "$pdb_files/${new_name}.pdb"
done


# Remove non-protein lines (in target)
target=$(find $pdb_files -type f | grep -E "*target.pdb")

for file in $target
do
	echo "Removing HETATM in $file..."
	name=$(basename "$file" | cut -f 1 -d '.')
	new_name=$(echo "$name" | sed 's/_target/_target_a/')
	grep '^ATOM' $file > "$pdb_files/${new_name}.pdb"
done

# Rename CRBN chains (to R)
files_crbn=$(find $pdb_files -type f | grep -E "*crbn.pdb")
for file in $files_crbn
do
	echo "Replacing chain in $file..."
	name=$(basename "$file" | cut -f 1 -d '.')
	awk '($1=="ATOM" || $1=="ANISOU") {$0 = substr($0, 1, 21) "R" substr($0,23)} 1' $file > "$pdb_files/${name}_r.pdb"
done


# Rename target chains (to L)
files_target=$(find $pdb_files -type f | grep -E "*target_a.pdb")
for file in $files_target
do
        echo "Replacing chain in $file..."
        name=$(basename "$file" | cut -f 1 -d '.')
        awk '($1=="ATOM" || $1=="ANISOU") {$0 = substr($0, 1, 21) "L" substr($0,23)} 1' $file > "$pdb_files/${name}_l.pdb"
done


# Renumber non-human CRBN files
#Chicken
awk '($1=="ATOM" || $1=="ANISOU") {sub($6, $6 - 2);sub($5, "R")} 1' "$pdb_files/4ci1_crbn.pdb" > "$pdb_files/4ci1_crbn_r.pdb"

#Mouse
awk '($1=="ATOM" || $1=="ANISOU") {sub($6, $6 - 3);sub($5, "R")} 1' "$pdb_files/5yiz_crbn.pdb" > "$pdb_files/5yiz_crbn_r.pdb"


# Join the ligand and receptor files in a single one (crbn_r+target_a)
lista=("5fqd" "5hxb" "6bn7" "6h0f" "6h0g" "6uml" "7lps")
for i in ${lista[@]}
do
	cat $pdb_files/${i}_crbn_r.pdb $pdb_files/${i}_target_a.pdb > $pdb_files/${i}_crbn_target.pdb	

done

