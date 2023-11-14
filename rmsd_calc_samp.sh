#!/bin/bash

## Script function: Automatic RMSD calculation for the complexes obtained with LightDock

# Activate conda environment
#conda activate maria

# Root directory
root="/home/l061003/TFM_MariaSantamera/"

# Files and directoriesi
pdb_files="${root}validation_set/input_pdb/"
docking_files="/home/l061003/Documents/docking_200/" # ARGUMENT
swarm_number=284 #ARGUMENT

list=("5fqd" "5hxb" "6h0f" "6h0g" "6uml" "7lps")

for i in ${list[@]}
do
    for swarm in {0..${swarm_number}}
    do	
	for glow in {0..199}
	do
           original_file="${pdb_files}${i}_crbn_target.pdb"
           original="${i}_crbn_target"
           dock_file="${docking_files}${i}/swarm_${swarm}/lightdock_${glow}.pdb"
           dock="lightdock_${glow}"
           echo "load ${original_file}, quiet = 1; load ${dock_file}, quiet = 1; sele chain R and ${original}, quiet = 1; align ${dock}, sele, quiet = 1; rms_cur ${dock}, ${original}" > rmsd_test_pymol_${i}.pml
	   pymol -c rmsd_test_pymol_${i}.pml >>PYMOL_${i}.log 2>>PYMOL_error_${i}.out
	done
    done
    #grep "Executive: RMSD =" PYMOL_${i}.log | cut -f7 -d " " > RMSD_${i}
    #max=$(sort -n RMSD_${i} | head -n 1)
    #echo "Max RMSD = $max" > temporal.txt
    #cat "RMSD_${i}" >> temporal.txt
    #mv temporal.txt RMSD_${i}
done

