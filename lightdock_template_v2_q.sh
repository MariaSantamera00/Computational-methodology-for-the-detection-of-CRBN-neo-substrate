#!/bin/bash

## This script uses Lightdock for docking of the proteins of the validation set using following parameters: 
	# Self-docking
        # No flexibility
        # Dfire
        # CRBN restrictions (rest_crbn file)
        # Ligand restrictions (beta-hairpin loop with Gly)
        # 100 simulations
# Script designed to be sent to queue


# Activate conda environment 
#conda activate maria

# Function to print help in terminal
print_help(){
        echo "ligthdock_template_v2_q.sh"
        echo "María Santamera Lastras 2023"
        echo -e "\nusage:ligthdock_template_v2_q.sh <root_folder> <output_folder>\n"
        echo -e "\troot_folder : root folder (should contain the GitHub files)"
        echo -e "\toutput_folder : folder in which data and results will be stored"
}



# Argument assignation
root=$1 
docking_files=$2 #results="${root}ligthdock_loop_MSL_results/"

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
rest_crbn="${root}validation_set/rest_crbn"
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


#Ligand (target) and receptor (CRBN)
pdb=$(sed -n "${SGE_TASK_ID}p" $names_pdb_files)
ipdb=$(basename "$pdb" | cut -f 1 -d '.')

ligand="${pdb_files}${ipdb}_target_a_l.pdb"
receptor="${pdb_files}${ipdb}_crbn_r_phe.pdb"


#GET RESTRICTION FILES (LIGAND)
#Get dssp file
mkdssp -i ${ligand} -o ${pdb_files}${ipdb}_target_a_l_dssp.txt
archivo_dssp="${pdb_files}${ipdb}_target_a_l_dssp.txt"


#Write on a single line SS characters
cat ${archivo_dssp} | awk '/^  #/{p=1;next}p' | awk '{print substr($0, 17, 1)}' | awk '{printf "%s", $0}' > temp_1_${ipdb}.txt

#Look for the degrons loop pattern in the previous file
cat temp_1_${ipdb}.txt | grep -o -E "E+[TS ]{5,7}E+" > temp_2_${ipdb}.txt

count=0
while IFS= read -r patron
do

	#Look for the position of the motifs found in the sequence. 
        init=$(grep -bo "$patron" temp_1_${ipdb}.txt | awk -F ":" '{print $1 +1}')
        end=$(grep -bo "$patron" temp_1_${ipdb}.txt | awk -F ":" '{print $1 + length($2)-1+1}')
	#Look for the glycine
        pattern=$(cat ${archivo_dssp} | awk '/^  #/{p=1;next}p' | sed -n "${init},${end}p" | awk 'substr($0, 17, 1) ~ /[TS ]/{print $0}' | awk 'NR>=4 && NR<=6 && substr($0, 14, 1) == "G" {print $0}' | wc -l)
	
	if [ "$pattern" -ne 0 ]
        then
                count=$((count+1))
		cat ${archivo_dssp} | awk '/^  #/{p=1;next}p' | sed -n "${init},${end}p" | awk 'substr($0, 17, 1) ~ /[TS ]/{print $0}' | awk '{print "L L." substr($0, 14, 1) "." substr($0, 8, 3)}' > temp_${ipdb}.txt
	
	#Replace G with Gly (all aminoacids) in the restrictions file
        if [ -f "${docking_files}rest_lig_${ipdb}_${count}.txt" ]
	then 
		rm "${docking_files}rest_lig_${ipdb}_${count}.txt"
	fi
		
	while IFS= read -r line
        do
                aa1=$(echo $line | cut -f2 -d ".")
		aa3=$(python ${root}bin/transl_python.py "$aa1")
		echo $line | awk -F "." -v OFS='.' -v var1="$aa1" -v var2="$aa3" '{$2 = ( $2 ==var1 ? var2 : $2 ); print $0}'
        done < temp_${ipdb}.txt > "${docking_files}rest_lig_${ipdb}_${count}.txt"

	fi

done < temp_2_${ipdb}.txt

rm temp_${ipdb}.txt
rm temp_1_${ipdb}.txt
rm temp_2_${ipdb}.txt


if [ "$count" -ne 0 ]
then	
	for (( i=1; i<=count; i++))
	do
		# Directory for docking results

		if [ ! -d "${docking_files}${ipdb}_${i}" ]
		then
        		mkdir "${docking_files}${ipdb}_${i}"
		fi
		cd "${docking_files}${ipdb}_${i}"
		
		# Copy ligand and receptor files
		cp "${pdb_files}${ipdb}_target_a_l.pdb" "${ipdb}_target_a_l.pdb"
		cp "${pdb_files}${ipdb}_crbn_r_phe.pdb" "${ipdb}_crbn_r_phe.pdb"
		
		rec="${ipdb}_crbn_r_phe.pdb"
		lig="${ipdb}_target_a_l.pdb"
		
		# Second proof of concept.
		        # Self-docking 
        		# No flexibility
        		# Dfire
        		# CRBN restrictions
        		# Ligand restrictions
        		# 100 simulations
				
		#Step 0. Join rest files
		sed -i 's/\.\s*/\./g' "${docking_files}rest_lig_${ipdb}_${i}.txt"
		cat $rest_crbn "${docking_files}rest_lig_${ipdb}_${i}.txt" > rest_file
	
		#Step 1. Setup
		lightdock3_setup.py --noxt --noh --now -r rest_file ${rec} ${lig}

		#Step 2. Simulation
		lightdock3.py -s dfire -min setup.json 100

		#Step 3. Generating structures
		#Step 4. Clustering structures intra-swarm
		#Step 5. Generating rank and filtering

		### Calculate the number of swarms ###
		CORES=8
		s=`ls -d ./swarm_* | wc -l`
		swarms=$((s-1))

		### Create files for Ant-Thony ###
		for j in $(seq 0 $swarms)
		  do
		    echo "cd swarm_${j}; lgd_generate_conformations.py ${docking_files}${ipdb}_${i}/${rec} ${docking_files}${ipdb}_${i}/${lig}  gso_100.out 200 > /dev/null 2> /dev/null;" >> generate_lightdock.list;
		    echo "cd swarm_${j}; lgd_cluster_bsas.py gso_100.out > /dev/null 2> /dev/null;" >> cluster_lightdock.list;
		
		  done
		
		### Generate LightDock models ###
		ant_thony.py -c ${CORES} generate_lightdock.list;
		### Clustering BSAS (rmsd) within swarm ###
		ant_thony.py -c ${CORES} cluster_lightdock.list;
		### Generate ranking files for filtering ###
		lgd_rank.py $s 100;
		### Filtering models by >40% of satisfied restraints ###
		lgd_filter_restraints.py --cutoff 5.0 --fnat 0.8 rank_by_scoring.list rest_file A B > /dev/null 2> /dev/null;
	
	done
fi

