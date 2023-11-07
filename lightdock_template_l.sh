#!/bin/bash

## This script uses Lightdock for docking of the proteins of the validation set using following parameters:
        # Self-docking
        # No flexibility
        # Dfire
        # CRBN restrictions (rest_crbn file)
        # NO ligand restrictions
        # 100 simulations



# Activate conda environment
conda activate maria

# Root directory
root="/home/l061003/TFM_MariaSantamera/"


# Files and directories
pdb_files="${root}validation_set/input_pdb/"
rest_crbn="${root}validation_set/rest_crbn"
names_pdb_files="${root}validation_set/test.pdbs"

# Results directory

if [ ! -d "${root}lightdock_v1_loop_MSL_results/" ]
then
        mkdir "${root}lightdock_v1_loop_MSL_results/"
fi

docking_files="${root}lightdock_v1_loop_MSL_results"


list=("5fqd" "5hxb" "6bn7" "6h0f" "6h0g" "6uml" "7lps")

for pdb in ${list[@]}
do

# Ligand (taget) and receptor (CRBN)
ligand="${pdb_files}${pdb}_target_a_l.pdb"
receptor="${pdb_files}${pdb}_crbn_r_phe.pdb"



if [ ! -d "${docking_files}/${pdb}" ]
then
	mkdir "${docking_files}/${pdb}"
fi

cd "${docking_files}/${pdb}"

# First proof of concept. 
	# Self-docking without target restrictions
	# No flexibility
	# Dfire
	# CRBN restrictions 
	# 100 simulations
#Step 1. Setup
lightdock3_setup.py --noxt --noh --now -r $rest_crbn ${receptor} ${ligand}

#Step 2. Simulation
lightdock3.py -s dfire setup.json 100

#Step 3. Generating structures
#Step 4. Clustering structures intra-swarm
#Step 5. Generating rank and filtering

### Calculate the number of swarms ###
CORES=8
s=`ls -d ./swarm_* | wc -l`
swarms=$((s-1))

### Create files for Ant-Thony ###
for i in $(seq 0 $swarms)
  do
    echo "cd swarm_${i}; lgd_generate_conformations.py ${receptor} ${ligand}  gso_100.out 200 > /dev/null 2> /dev/null;" >> generate_lightdock.list;
    echo "cd swarm_${i}; lgd_cluster_bsas.py gso_100.out > /dev/null 2> /dev/null;" >> cluster_lightdock.list;

  done

### Generate LightDock models ###
ant_thony.py -c ${CORES} generate_lightdock.list;
### Clustering BSAS (rmsd) within swarm ###
ant_thony.py -c ${CORES} cluster_lightdock.list;
### Generate ranking files for filtering ###
lgd_rank.py $s 100;
### Filtering models by >40% of satisfied 		restraints ### 
lgd_filter_restraints.py --cutoff 5.0 --fnat 0.8 rank_by_scoring.list $rest_crbn A B > /dev/null 2> /dev/null;

done
