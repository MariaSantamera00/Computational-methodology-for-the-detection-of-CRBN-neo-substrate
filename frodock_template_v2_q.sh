#!/bin/bash

## This script uses Frodock for docking of the proteins of the validation set using following parameters:
        # Self-docking
        # CRBN restrictions (mid-point between the most important residues in the interaction)
# Script designed to be sent to queue


# Activate conda environment
#conda activate maria

# Function to print help in terminal
print_help(){
        echo "frodock_template_v2_q.sh"
        echo "María Santamera Lastras 2023"
        echo -e "\nusage:frodock_template_v2_q.sh <root_folder> <output_folder>\n"
        echo -e "\troot_folder : root folder (should contain the GitHub files)"
        echo -e "\toutput_folder : folder in which data and results will be stored"
}



# Argument assignation
root=$1 
docking_files=$2 

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
pdb_files="${root}validation_set/input_pdb/"
names_pdb_files="${root}validation_set/test.pdbs"

# Results directory
if [ -d "${docking_files}" ]
then
        echo "The folder exists"
        echo "Do you want to delete the previous results and recreate the directory? y/n"
        read answer

        # The directory will only be deleted and recreated if the user says "yes"
        case $answer in
                y)
                        rm -r ${docking_files}
                        mkdir ${docking_files}
                ;;
                n)
                        exit
                ;;
                *)
                        echo "You must choose y/n"
                ;;
        esac
else
        mkdir ${docking_files}
fi



# Ligand (target) and receptor (CRBN)
pdb=$(sed -n "${SGE_TASK_ID}p" $names_pdb_files)
receptor="${pdb_files}${pdb}_target_a_l.pdb"
ligand="${pdb_files}${pdb}_crbn_r.pdb"



if [ ! -d "${docking_files}${pdb}" ]
then
        mkdir "${docking_files}${pdb}"
fi

cd "${docking_files}${pdb}"

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






