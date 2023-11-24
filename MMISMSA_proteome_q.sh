#!/bin/bash

## This script uses GROMACS to perform a binary complex minimization and then calculate the binding energy of that system with MMISMSA.


# Script designed to be sent to queue

# Load GROMACS
#module load gromacs/2022.3

# Activate conda environment
#conda activate maria

# Function to print help in terminal
print_help(){
        echo "MMISMSA_proteome_q.sh.sh"
        echo "María Santamera Lastras 2023"
        echo -e "\nusage:MMISMSA_proteome_q.sh <root_folder> <AF_files_folder> <Input_AF_files_names> <output_folder> <loops_folder>\n"
        echo -e "\troot_folder : root folder (should contain the GitHub validation set)"
        echo -e "\tAF_files_folder : folder in which AlphaFold files are stored"
        echo -e "\tInput_AF_files_names : file with the file names to be used to generate the output file"
        echo -e "\toutput_folder : folder in which data and results will be stored"
}


# Arguments assignation
root=$1 #root="/home/l061003/"
dir=$2 # dir="/lrlhps/users/l001803/TMP/crbn_eval/results/prime_output/"
names_files=$3 #names_files="/lrlhps/users/l001803/TMP/crbn_eval/results/temp_6756_protein_filter_lys.txt"
output_dir=$4 #output_dir="/lrlhps/users/l001803/TMP/crbn_eval/results/phcal_3/"


#Controls help message (print "print_help" function if user writes "-h" or "-help" in terminal)
if [ "$root" == "-h" ] || [ "$root" == "-help" ] ;then
        print_help
        exit
fi


#Control of arguments (show an error message if user doesn't write enough arguments)
if [ "$#" -ne 4 ]; then
    echo -e "ERROR: too few arguments\n"
    print_help
    exit
fi


# Files and directories
min_mdp="${root}/validation_set/gromacs_input_files/min.mdp"
ion_mdp="${root}/validation_set/gromacs_input_files/ion.mdp"


# Results directory
if [ -d "${output_dir}" ]
then
        echo "The folder exists"
        echo "Do you want to delete the previous results and recreate the directory? y/n"
        read answer

        # The directory will only be deleted and recreated if the user says "yes"
        case $answer in
                y)
                        rm -r ${output_dir}
                        mkdir ${output_dir}
                ;;
                n)
                        exit
                ;;
                *)
                        echo "You must choose y/n"
                ;;
        esac
else
        mkdir ${output_dir}
fi


mkdir "${output_dir}pdb_files"



pdb=$(sed -n "${SGE_TASK_ID}p" $names_files)
pdb_file="${dir}${pdb}"
min=$(basename "$pdb_file" | cut -f 1 -d '.')
dir=$(pwd)



# Remove molecular glue
cat $pdb_file | grep -v "HETATM" | grep -v "ACE" | grep -v "NMA" > "${output_dir}pdb_files/${min}_nohet.pdb"
pdb_nohet="${output_dir}pdb_files/${min}_nohet.pdb"

# Minimización GROMACs
mkdir ${output_dir}${min}
cd ${output_dir}${min}
echo "6" | gmx pdb2gmx -f ${pdb_nohet} -ignh -water tip3p

gmx editconf -f conf.gro -o box.gro -d 1.2

gmx solvate -cp box.gro -cs spc216 -o solv.pdb -p topol.top

gmx grompp -f ${ion_mdp} -c solv.pdb -o ion -p topol.top > charge_output.txt 2>&1

charge=$(cat charge_output.txt | grep "System has non-zero total charge:" | awk '{print $NF}' | cut -f1 -d ".")

if (( $(echo "$charge > 0" |bc -l) ))
then
	echo "13" | gmx genion -s ion -o neutral.pdb -p topol.top -nn ${charge}
elif (( $(echo "$charge < 0" |bc -l) ))
then
	charge_abs=$(echo $charge | sed 's/^\-//')
	echo "13" | gmx genion -s ion -o neutral.pdb -p topol.top -np ${charge_abs}
else
	cp solv.pdb neutral.pdb
fi


gmx grompp -f ${min_mdp} -o min -c neutral.pdb -p topol.top

gmx mdrun -v -deffnm min -nt 16


# Convert gromac files to amber python -m install parmed
python ${root}parmed_ph.py "topol.top" "min.gro" "min_amber.top" "min_amber.crd" 


# Search for hydrogen bridges
cMMISMSA/src/MM --topology min_amber.top --rst min_amber.crd --output mm_result --mask 108-50000


cd $dir

