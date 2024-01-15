#!/bin/bash

## This script generates the necessary input to run Schrodinger's Prime and executes the program. It also allows to get the input for the dfe program (the configuration file).

# Script designed to be sent to queue



# Load Schrodinger
#module load schrodinger/2022.3

# Activate conda environment
#conda activate maria

# Function to print help in terminal
print_help(){
        echo "prime_sidechain_3a_proteome_q.sh"
        echo "María Santamera Lastras 2023"
        echo -e "\nusage:prime_sidechain_3a_proteome_q.sh <root_folder> <AF_files_folder> <Input_AF_files_names> <output_folder> <loops_folder>\n"
        echo -e "\troot_folder : root folder (should contain the GitHub validation set)"
        echo -e "\tAF_files_folder : folder in which AlphaFold files are stored"
        echo -e "\tInput_AF_files_names : file with the file names to be used to generate the output file"
        echo -e "\toutput_folder : folder in which data and results will be stored"
        echo -e "\tloops_folder : folder in which AlphaFold files loops ares stored"
}



# Arguments assignation
root=$1 
input_dir=$2 
names_pdb_files=$3 
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



dir=$(pwd)

comb=$(cat ${names_pdb_files} | wc -l)
for num in $(seq 1 $comb)
do



ap=$(sed -n "${num}p" $names_pdb_files)
iap=$(basename "$ap" | cut -f 1 -d '.')

crbn_dir="${root}validation_set/input_pdb/"
crbn="5fqd_crbn_r_thalidomide"


## PREPROCESSING PDB FILES WITH PREPWIZARD --> MAE FILES

if [ ! -d "${output_dir}prep_mae/" ]
then
        mkdir "${output_dir}prep_mae/"
fi

cd "${output_dir}prep_mae/"

schrodinger/schrodinger-2023.3/utilities/prepwizard -captermini -fillloops -fillsidechains -disulfides -propka_pH 7.0 -rmsd 0.6 -WAIT "${input_dir}${iap}.pdb" "${output_dir}prep_mae/${iap}.mae" 

schrodinger/schrodinger-2023.3/utilities/prepwizard -captermini -fillloops -fillsidechains -disulfides -propka_pH 7.0 -rmsd 0.3 -WAIT  "${crbn_dir}${crbn}.pdb" "${output_dir}prep_mae/${crbn}.mae"

cd "${dir}"

## MERGE MAE FILES --> MAE FILE

if [ ! -d "${output_dir}merge/" ]
then
	mkdir "${output_dir}merge/"
fi


run merge.py "${output_dir}prep_mae/${crbn}.mae" "${output_dir}prep_mae/${iap}.mae" "${output_dir}merge/${iap}_${crbn}.mae"


## CREATE prime_sidechain.inp FILE
# Create a directory
if [ ! -d "${output_dir}prime_sidechain/" ]
then
        mkdir "${output_dir}prime_sidechain/"
fi

# Write the header of the file

echo -e "STRUCT_FILE\t${output_dir}merge/${iap}_${crbn}.mae\nJOB_TYPE\tREFINE\nPRIME_TYPE\tSIDE_PRED\nSELECT\t pick" > "${output_dir}prime_sidechain/header_ps_${iap}${crbn}.inp"


# Write the body of the file (receptor lines)
res_crbn=(351 352 353 355 357 375 377 378 390 397)
chain_r=$(cat ${crbn_dir}${crbn}.pdb | awk '{print $5}' | head -n1)

count=-1
for i in ${res_crbn[@]}
do
	count=$((count+1))
	echo -e "RESIDUE_${count}\t${chain_r}:${i}" 
done > "${output_dir}prime_sidechain/bodyr_ps_${iap}${crbn}.inp"

# Write the body of the file (ligand lines)
# Convert .mae files to .pdb files

run schrodinger/schrodinger-2023.3/utilities/pdbconvert -imae "${output_dir}prep_mae/${iap}.mae" -opdb "${output_dir}prep_mae/${iap}_mae.pdb"
run schrodinger/schrodinger-2023.3/utilities/pdbconvert -imae "${output_dir}prep_mae/${crbn}.mae" -opdb "${output_dir}prep_mae/${crbn}_mae.pdb"

# Residues < 3A
xyz_lig=$(echo "[")

cat ${output_dir}prep_mae/${iap}_mae.pdb | awk '$1=="ATOM"{print $0}' | awk '{print substr($0, 24)}' > ${output_dir}prep_mae/${iap}_mae_atom.pdb
while IFS= read -r ligand_line
do

	xyz_lig=$(echo "${xyz_lig}$(echo "$ligand_line" | awk '($1=="ATOM" || $1=="HETATM"){print "["substr($0, 32, 8)","substr($0, 39, 8)","substr($0, 47, 8)"]"}' | sed 's/ //g')")

done < "${output_dir}prep_mae/${iap}_mae.pdb"

xyz_lig=$(echo "$(echo "${xyz_lig}" | sed 's/]\[/],[/g')]")
echo "$xyz_lig" > prueba_lig.txt


xyz_cer=$(echo "[")

cat ${output_dir}prep_mae/${crbn}_mae.pdb | awk '$1=="ATOM"{print $0}' | awk '{print substr($0, 24)}' > ${output_dir}prep_mae/${crbn}_mae_atom.pdb

while IFS= read -r crbn_line
do
	xyz_cer=$(echo "${xyz_cer}$(echo "$crbn_line" | awk '($1=="ATOM" || $1=="HETATM"){print "["substr($0, 32, 8)","substr($0, 39, 8)","substr($0, 47, 8)"]"}'| sed 's/ //g')")

done < "${output_dir}prep_mae/${crbn}_mae.pdb"

xyz_cer=$(echo "$(echo "${xyz_cer}" | sed 's/]\[/],[/g')]")
echo "$xyz_cer" > prueba_cer.txt


python ${root}bin/distances_file.py prueba_lig.txt prueba_cer.txt > "temp_${iap}${crbn}.pdb"

rm prueba_lig.txt
rm prueba_cer.txt

awk 'NR==FNR { lines[$0]; next } FNR in lines' "temp_${iap}${crbn}.pdb" "${output_dir}prep_mae/${iap}_mae_atom.pdb" | awk '{print $1}'| sort | uniq > "temp_3a_${iap}${crbn}.pdb"



#Loop residues

pdb_id=$(echo "$iap" | cut -f1,2,3,4 -d "_")
loop="${loop_dir}${pdb_id}_loop_al.pdb"
cat ${loop} | awk '{print substr($0, 24, 3)}'| sed 's/ //g' | sort | uniq > "temp_loop_${iap}${crbn}.pdb"



lig_crbn=$(cat "temp_3a_${iap}${crbn}.pdb" "temp_loop_${iap}${crbn}.pdb" | sort | uniq | sort -n | tr "\n" " " | sed 's/$//')

chain_l=$(cat ${loop} | awk '{print $5}' | head -n1)


for i in ${lig_crbn[@]}
do
        count=$((count+1))
	echo -e "RESIDUE_${count}\t${chain_l}:${i}"
done > "${output_dir}prime_sidechain/bodyl_ps_${iap}${crbn}.inp"

rm "temp_${iap}${crbn}.pdb"
rm "temp_loop_${iap}${crbn}.pdb"
rm "temp_3a_${iap}${crbn}.pdb"

# Write the final lines of the file

echo -e "NUM_SC_OUTPUT_STRUCT\t1\nUSE_CRYSTAL_SYMMETRY\tno\nUSE_RANDOM_SEED\tno\nSEED\t0\nOPLS_VERSION\tS-OPLS\nEXT_DIEL\t80.00\nUSE_MEMBRANE\tno" > "${output_dir}prime_sidechain/end_ps_${iap}${crbn}.inp"


# Join all the lines in a single file

cat "${output_dir}prime_sidechain/header_ps_${iap}${crbn}.inp" "${output_dir}prime_sidechain/bodyr_ps_${iap}${crbn}.inp" "${output_dir}prime_sidechain/bodyl_ps_${iap}${crbn}.inp" "${output_dir}prime_sidechain/end_ps_${iap}${crbn}.inp" > "${output_dir}prime_sidechain/ps_${iap}_${crbn}.inp"


# Remove unnecessary files
rm "${output_dir}prime_sidechain/header_ps_${iap}${crbn}.inp"
rm "${output_dir}prime_sidechain/bodyr_ps_${iap}${crbn}.inp"
rm "${output_dir}prime_sidechain/bodyl_ps_${iap}${crbn}.inp"
rm "${output_dir}prime_sidechain/end_ps_${iap}${crbn}.inp"






## MINIMIZATION
if [ ! -d "${output_dir}prime_output/" ]
then
        mkdir "${output_dir}prime_output/"
fi

cd "${output_dir}prime_output/"

schrodinger/schrodinger-2022.3/prime "${output_dir}prime_sidechain/ps_${iap}_${crbn}.inp" -HOST cluster_mpi:50 -NJOBS 8


cd "${dir}"

## CREATE A FILE WITH CRBN RESIDUES (METADYNAMICS)
crbn_line="(res.n 380 and chain.name R) or (res.n 386 and chain.name R) or (res.n 400 and chain.name R)"


## CREATE A FILE WITH LIGAND LOOP RESIDUES (METADYNAMICS)
pdb_id=$(echo "$iap" | cut -f1,2,3,4 -d "_")
loop="${loop_dir}${pdb_id}_loop_al.pdb"
cat ${loop} | awk '{print $6}' | sort | uniq > "temp_loop_${iap}${crbn}.pdb"

chain_l=$(cat ${loop} | awk '{print $5}' | head -n1)
loop_res=$(cat "temp_loop_${iap}${crbn}.pdb" | sort | uniq | sort -n | tr "\n" " " | sed 's/$//')

lig_line=""
for i in ${loop_res[@]}
do
        lig_line=$(echo "${lig_line}$(echo "(res.n ${i} and chain.name ${chain_l}) or ")")

done



rm "temp_loop_${iap}${crbn}.pdb"

## CRETATE CONFIG FILE (METADYNAMICS)
if [ ! -d "${output_dir}config_run/" ]
then
        mkdir "${output_dir}config_run/"
fi


echo -e "{\"DEFAULT\":\n  {\"pocket_asl\":  \"$(echo "${crbn_line}")\",\n   \"ligand_asl\": \"$(echo "${lig_line::-4}")\",\n   \"wall\": 80.0,\n   \"width\": 0.08,\n   \"height\": 0.02,\n   \"interval\": 0.09,\n   \"time\": 15000.0,\n   \"args\":\n     {\"fl\": 2.9,\n      \"fr\": 4.0,\n      \"fd\": 3.2,\n      \"rd\": null,\n      \"cs\": 1000.0,\n      \"bs\": 0.2,\n      \"im\": \"rectangle\",\n      \"fm\": \"mask_mean\",\n      \"cwp\": false,\n      \"dcf\": 0.02,\n      \"ngt\": 0.05,\n      \"bt\": true},\n    \"correl\":\n      {\"a\": -0.09,\n       \"b\": 0.51}\n  }\n}" > "${output_dir}config_run/config_run_${iap}_${crbn}.json"



done
