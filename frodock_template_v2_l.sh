#!/bin/bash

## This script uses Frodock for docking of the proteins of the validation set using following parameters:
        # Self-docking
        # CRBN restrictions (mid-point between the most important residues in the interaction)



# Activate conda environment
conda activate maria

# Root directory
root="/home/l061003/TFM_MariaSantamera/"

# Files and directories
pdb_files="${root}validation_set/input_pdb/"


# Results directory
if [ ! -d "${root}frodockv2_loop_MSL_results/" ]
then
        mkdir "${root}frodockv2_loop_MSL_results/"
fi

docking_files="${root}frodockv2_loop_MSL_results"


list=("5fqd" "5hxb" "6bn7" "6h0f" "6h0g" "6uml" "7lps")

for pdb in ${list[@]}
do

# Ligand (target) and receptor (CRBN)
receptor="${pdb_files}${pdb}_target_a_l.pdb"
ligand="${pdb_files}${pdb}_crbn_r.pdb"



if [ ! -d "${docking_files}/${pdb}" ]
then
        mkdir "${docking_files}/${pdb}"
fi

cd "${docking_files}/${pdb}"

cp ~/.conda/envs/maria/bin/frodock3_linux64/bin/soap.bin soap.bin



##1. CHECK PDB INPUT FILES
pdb2pqr --ff=CHARMM ${receptor} ${pdb}_crbn_r_pqr.pdb
pdb2pqr --ff=CHARMM ${ligand} ${pdb}_target_a_pqr.pdb

receptor=${pdb}_crbn_r_pqr.pdb
ligand=${pdb}_target_a_pqr.pdb


##2. POTENTIAL MAP GENERATION
	#Creation of receptor vdw potential map:
	frodock3_linux64/bin/frodockgrid ${receptor} -o ${pdb}_crbn_W.ccp4

	#Creation of the receptor electrostatic potential map:
	frodock3_linux64/bin/frodockgrid ${receptor} -o ${pdb}_crbn_E.ccp4 -m 1

	#Creation of the receptor desolvation potential map:
	frodock3_linux64/bin/frodockgrid ${receptor} -o ${pdb}_crbn_DS.ccp4 -m 3

	#Creation of the ligand desolvation potential map:
	frodock3_linux64/bin/frodockgrid ${ligand} -o ${pdb}_target_DS.ccp4 -m 3


##3. PERFORM THE EXHAUSTIVE DOCKING SEARCH
frodock3_linux64/bin/frodock ${pdb}_crbn_r_pqr_ASA.pdb ${pdb}_target_a_pqr_ASA.pdb -w ${pdb}_crbn_W.ccp4 -e ${pdb}_crbn_E.ccp4 --th 10 -d ${pdb}_crbn_DS.ccp4,${pdb}_target_DS.ccp4 -o dock.dat --around -59.0,52.0,118.0


##4. CLUSTERING AND VISUALIZATION OF PREDICTIONS
	#Sort de solutions by the RMSD value
	frodock3_linux64/bin/frodockcluster dock.dat ${ligand} --nc 100 -o clust_dock.dat
		# The RMSD default value used in the clusterization is 5Å.

	#To visualize the 10th first solutions:
	frodock3_linux64/bin/frodockview clust_dock.dat -r 1-10

	#Coordinate generation of the predicted solutions:
	frodock3_linux64/bin/frodockview clust_dock.dat -r 1-5 -p ${ligand}

done




