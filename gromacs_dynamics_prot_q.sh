#!/bin/bash

## This script performs a classical molecular dynamics with gromacs (minimization, balancing and simulation) and allows to analyze the obtained results of energy, RMSD distances. 

# Script designed to be sent to queue

# Load GROMACS
#module load gromacs/2022.3

# Activate conda environment
#conda activate maria

# Function to print help in terminal
print_help(){
        echo "gromacs_dynamics_prot_q.sh"
        echo "María Santamera Lastras 2023"
        echo -e "\nusage:gromacs_dynamics_prot_q.sh <root_folder> <AF_files_folder> <Input_AF_files_names> <output_folder> <loops_folder>\n"
        echo -e "\troot_folder : root folder (should contain the GitHub validation set)"
        echo -e "\tAF_files_folder : folder in which AlphaFold files are stored"
        echo -e "\tInput_AF_files_names : file with the file names to be used to generate the output file"
	echo -e "\toutput_folder : folder in which data and results will be stored"
	echo -e "\tloops_folder : folder in which AlphaFold files loops ares stored"
}


# Arguments assignation
root=$1 
dir=$2 
names_files=$3 
output_dir=$4 
loop_dir=$5 

#Controls help message (print "print_help" function if user writes "-h" or "-help" in terminal)
if [ "$root" == "-h" ] || [ "$root" == "-help" ] ;then
        print_help
        exit
fi


#Control of arguments (show an error message if user doesn't write enough arguments)
if [ "$#" -ne 5 ]; then
    echo -e "ERROR: too few arguments\n"
    print_help
    exit
fi


# Files and directories
min_mdp="${root}/validation_set/gromacs_input_files/min.mdp" 
ion_mdp="${root}/validation_set/gromacs_input_files/ion.mdp"
npt_mdp="${root}/validation_set/gromacs_input_files/npt.mdp"
nvt_mdp="${root}/validation_set/gromacs_input_files/nvt.mdp"
md_mdp="${root}/validation_set/gromacs_input_files/md.mdp"


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
pdb_file="${dir}${pdb}.pdb"
min=$(basename "$pdb_file" | cut -f 1 -d '.')
dir=$(pwd)



# Remove molecular glue
cat $pdb_file | grep -v "HETATM" | grep -v "ACE" | grep -v "NMA" > "${output_dir}pdb_files/${min}_nohet.pdb"
pdb_nohet="${output_dir}pdb_files/${min}_nohet.pdb"


if [ ! -d "${output_dir}phcal_3${min}" ]
then
	mkdir "${output_dir}phcal_3${min}"
fi

cd "${output_dir}phcal_3${min}"


# Minimization GROMACs

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


# NVT Equilibration 
gmx grompp -f ${nvt_mdp} -c min -r min -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt #-nb gpu

# NPT Equilibration
gmx grompp -f ${npt_mdp} -c nvt -r nvt -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt #-nb gpu


# Production MD
gmx grompp -f ${md_mdp} -c npt -t npt -p topol.top -o md_0_1.tpr
gmx mdrun -deffnm md_0_1 #-nb gpu


# Results analysis
#source /lrlhps/users/l001803/TMP/gmx_old/bin/GMXRC.bash

# Calculate energy

echo "1" | gmx trjconv -f md_0_1.xtc -s md_0_1.tpr -o md_pbc1.xtc -pbc whole

echo "1" | gmx trjconv -f md_pbc1.xtc -s md_0_1.tpr -o md_pbc2.xtc -pbc nojump

echo "q" | gmx make_ndx -f min.gro
echo "[ crbn ]" > temp.ndx
for i in `seq 1 1663`; do echo "$i"; done >> temp.ndx
cat "temp.ndx" "index.ndx" > "index_temp.ndx"
rm "index.ndx"
mv "index_temp.ndx" "index.ndx"
rm temp.ndx

echo "0 2" | gmx trjconv -f md_pbc2.xtc -s md_0_1.tpr -o md_pbc3.xtc -fit rot+trans -n index.ndx

cMMISMSA/src/MM --topology min_amber.top --xtc md_pbc3.xtc --output md_result --mask 108-50000

# Calculate diference between energy before and after de simulation
init=$(cat "mm_result_energy" | tail -n1 | cut -f6 -d ";")
fin=$(cat "md_result_energy" | tail -n1 | cut -f6 -d ";")
av=$(cat "md_result_energy_averages" | tail -n1 | cut -f6 -d ";")

d_if=$(echo "$init - $fin" | bc)
d_if_ab=$(echo "if ($d_if < 0) - $d_if else $d_if" | bc)
d_av=$(echo "$init - $av" | bc)
d_av_ab=$(echo "if ($d_av < 0) - $d_av else $d_av" | bc)

echo -e "before-after;before-average\n${d_if_ab};${d_av_ab}" > dm_analysis.txt

# DISTANCES
# Create index file

# Calculate distances
pdb_id=$(echo "${min}" | sed 's/ps_//g' | sed 's/_al_5fqd_crbn_r_thalidomide-out//g')
loop="${loop_dir}${pdb_id}_loop_al.pdb"

cat ${loop} | awk '{print $6}' | sort | uniq > "loop_res_temp"

cat min.gro | tail -n +1666 | sed 's/.//6; s/.//6; s/.//6' > "target_temp"

echo "[ loop ]" > index_loop.ndx

awk 'NR==FNR{val[$1]; next} $1 in val' "loop_res_temp" "target_temp" | awk '{print $3}' >> index_loop.ndx

echo -e "\n[ crbn ]" >> index_loop.ndx
for i in `seq 1 1663`; do echo "$i"; done >> index_loop.ndx


rm loop_res
rm target_temp

# Calculate distances
echo -e "0\n1" | gmx pairdist -s md_0_1.tpr -f md_0_1.xtc -n index_loop.ndx -o distances

cat distances.xvg | grep -v "@" | grep -v "#" | awk '{print $2}' > distances_fil.txt


# RMSD
echo "1 0" | gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -center
echo "4 4" | gmx rms -s md_0_1.tpr -f md_0_1_noPBC.xtc -o rmsd.xvg -tu ns
echo "4 4" | gmx rms -s min.tpr -f md_0_1_noPBC.xtc -o rmsd_xtal.xvg -tu ns



cd $dir

