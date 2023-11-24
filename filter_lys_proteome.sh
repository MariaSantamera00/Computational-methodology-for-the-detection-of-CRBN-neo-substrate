#!/bin/bash

# This script checks the distance of the lysines from the targets (the entire proteome or previously filtered proteins) to the geometric center of the E2 ligase in the complex that forms CRBN together with DDB1, Cullin-4A, RBX1, NEDD8 and the UBC12 ligase.  

# Script designed to be sent to queue


# Activate conda environment
#conda activate maria

# Function to print help in terminal
print_help(){
        echo "filter_lys_proteome.sh"
        echo "María Santamera Lastras 2023"
        echo -e "\nusage:filter_lys_proteome.sh <AF_files_folder> <Input_AF_files_names> <root_folder> <output_folder>\n"
        echo -e "\tAF_files_folder : folder in which AlphaFold files are stored"
        echo -e "\tInput_AF_files_names : file with the file names to be used to generate the output file"
	echo -e "\troot_folder : root folder (should contain the GitHub repo folder with validation set and utilities)"
        echo -e "\toutput_folder : folder in which data and results will be stored"
}


# Argument assignation
alphafold_files=$1 #alphafold_files="/lrlhps/users/l001803/TMP/crbn_eval/results/prime_output/"
names_pdb_files=$2 #names_pdb_files="/lrlhps/users/l001803/TMP/crbn_eval/results/temp_6756_protein_filter_prime.txt"
root=$3 #root="/home/l061003/"
results=$4 #results="/lrlhps/users/l001803/TMP/crbn_eval/results/"
e2_complex="${root}validation_set/ubiq_complex.pdb"

# Controls help message (print "print_help" function if user writes "-h" or "-help" in terminal)
if [ "$alphafold_files" == "-h" ] || [ "$alphafold_files" == "-help" ] ;then
        print_help
        exit
fi


#Control of arguments (show an error message if user doesn't write enough arguments)
if [ "$#" -ne 4 ]; then
    echo -e "ERROR: too few arguments\n"
    print_help
    exit
fi


# Results directory
if [ -d "${results}" ]
then
        echo "The folder exists"
        echo "Do you want to delete the previous results and recreate the directory? y/n"
        read answer

        # The directory will only be deleted and recreated if the user says "yes"
        case $answer in
                y)
                        rm -r ${results}
                        mkdir ${results}
                ;;
                n)
                        exit
                ;;
                *)
                        echo "You must choose y/n"
                ;;
        esac
else
        mkdir ${results}
fi



e2_i=$(basename "$e2_complex"| cut -f 1 -d '.')

#name=$(sed -n "${SGE_TASK_ID}p" $names_pdb_files)
#i=$(echo "ps_${name}_al_5fqd_crbn_r_thalidomide-out")

comb=$(cat ${names_pdb_files} | wc -l)
for num in $(seq 1 $comb)
do

name=$(sed -n "${num}p" $names_pdb_files)
i=$(echo "ps_${name}_al_5fqd_crbn_r_thalidomide-out")



# Align the target protein with the ubiquitination complex

if [ ! -d "${results}align_files" ]
then
        mkdir "${results}align_files"
fi


echo "load ${alphafold_files}${i}.maegz, format=mae; load ${e2_complex}; select crbn, ${i} and chain R; align ${e2_i}, crbn; save ${results}align_files/${e2_i}_${i}.pdb,${e2_i}; save ${results}align_files/${i}.pdb,${i}" > pymol_${i}.pml

pymol -c pymol_${i}.pml >PYMOL.log 2>PYMOL_error.out

rm pymol_${i}.pml
rm PYMOL.log
rm PYMOL_error.out

# Distance of all lysines of the target (its alpha carbon) to the geometric center of UBC12

if [ ! -d "${results}distances" ]
then
        mkdir "${results}distances"
fi


# Target lysines: #OJO: CAMBIAR LA L DE LA CADENA POR A PARA ALPHAFOLD
xyz_lys=$(echo "[")

cat "${results}align_files/${i}.pdb" | awk '($5 == "A" && $4 == "LYS" && $3 == "CA"){print $0}' > "${results}distances/${i}_lys.pdb" 

while IFS= read -r lys_line
do

	xyz_lys=$(echo "${xyz_lys}$(echo "$lys_line" | awk '($1=="ATOM" || $1=="HETATM"){print "["substr($0, 32, 8)","substr($0, 39, 8)","substr($0, 47, 8)"]"}'| sed 's/ //g')")

done < "${results}distances/${i}_lys.pdb"

xyz_lys=$(echo "$(echo "${xyz_lys}" | sed 's/]\[/],[/g')]")
echo "$xyz_lys" > temp_lys.txt


# Geometrical center E2 CA:
#x=$(cat "${results}align_files/${e2_i}_${i}.pdb" | awk '($5 == "G" && $3 == "CA"){print substr($0, 32, 8)}' | awk '{sum += $1 } END {print sum/ NR}')
#y=$(cat "${results}align_files/${e2_i}_${i}.pdb" | awk '($5 == "G" && $3 == "CA"){print substr($0, 39, 8)}' | awk '{sum += $1 } END {print sum/ NR}')
#z=$(cat "${results}align_files/${e2_i}_${i}.pdb" | awk '($5 == "G" && $3 == "CA"){print substr($0, 47, 8)}' | awk '{sum += $1 } END {print sum/ NR}')

#xyz_e2=$(echo "[[${x},${y},${z}]]")
#echo "$xyz_e2" > temp_e2.txt



# E2 CAs: 
xyz_e2ca=$(echo "[")

cat "${results}align_files/${e2_i}_${i}.pdb" | awk '($5 == "G" && $3 == "CA"){print $0}' > "${results}distances/${e2_i}_${i}_lys.pdb"
while IFS= read -r e2_line
do

	xyz_e2ca=$(echo "${xyz_e2ca}$(echo "$e2_line" | awk '($1=="ATOM" || $1=="HETATM"){print "["substr($0, 32, 8)","substr($0, 39, 8)","substr($0, 47, 8)"]"}'| sed 's/ //g')")

done < "${results}distances/${e2_i}_${i}_lys.pdb"

xyz_e2ca=$(echo "$(echo "${xyz_e2ca}" | sed 's/]\[/],[/g')]")
echo "$xyz_e2ca" > temp_e2ca.txt




# Calculate distances
#~/.conda/envs/maria/bin/python /home/l061003/distances_calc.py temp_lys.txt temp_e2.txt > "${results}distances/dist_lys_${i}.txt"
python ${root}bin/distances_calc.py temp_lys.txt temp_e2ca.txt > "${results}distances/e2ca_dist_lys_${i}.txt"

rm temp_lys.txt
#rm temp_e2.txt
rm temp_e2ca.txt

rm "${results}distances/${i}_lys.pdb"
rm "${results}distances/${e2_i}_${i}_lys.pdb"
#done


#cat "${results}distances"/dist_lys_*.txt > "${results}total_dist_lys.txt"
#cat "${results}distances"/e2ca_dist_lys_*.txt > "${results}e2ca_total_dist_lys.txt"


# Filter
if [ ! -d "${results}filter" ]
then
        mkdir "${results}filter"
fi


treshold=10
hits=$(cat "${results}distances/e2ca_dist_lys_${i}.txt" | awk '($1 < 10){print $0}' | wc -l)
if [ "$hits" -ne 0 ]
then
	echo "${name}" > "${results}filter/${name}.txt"
fi


done





