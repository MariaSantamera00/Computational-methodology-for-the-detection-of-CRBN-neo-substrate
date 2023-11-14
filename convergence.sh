#!/bin/bash

## To study the convergence of the docking performed with LightDock, this script calculates the median, the maximum value and the mean of all the values in the "score" column of all the gso_*.out files (0, 10, 20...). A file is obtained for each complex generated by docking and this file allows to create a graph later for analysis. 


# Activate conda environment
#conda activate maria

# Function to print help in terminal
print_help(){
        echo "convergence.sh"
        echo "María Santamera Lastras 2023"
        echo -e "\nusage:convergence.sh <root_folder> <docking_folder>\n"
	echo -e "\troot_folder : root folder (should contain the GitHub validation set)"
        echo -e "\tdocking_folder : folder with docking results obtained with Lightdock"
}

# Argument assignation
root=$1 #root="/home/l061003/TFM_MariaSantamera/"
directory=$2 #directory="${root}lightdock_loop_MSL_results/"


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


list=("5fqd" "5hxb" "6bn7" "6h0f" "6h0g" "6uml" "7lps")


for i in ${list[@]}
do
	docking_directory="${directory}${i}" 
	mkdir ${docking_directory}/convergence
	convergence_directory="${docking_directory}/convergence"



	if [ -f "${convergence_directory}/scores_gso_total" ]
	then
		rm "${convergence_directory}/scores_gso_total"
	fi


	for step in {0..100..10} # Change it depending on the number of steps
	do
		if [ -f "${convergence_directory}/scores_gso_${step}" ]
		then
        		rm "${convergence_directory}/scores_gso_${step}"
		fi

		find ${docking_directory} -type f -name "gso_${step}.out" | while IFS= read -r line
		do
			awk '{print $13}' ${line} >> ${convergence_directory}/scores_gso_${step}
		done
		
		# Mean:
		cat ${convergence_directory}/scores_gso_${step}| awk '{sum += $1 } END {print sum/ NR}' >> ${convergence_directory}/scores_gso_total_mean
		# Max value:
		cat ${convergence_directory}/scores_gso_${step}| awk 'BEGIN{a= 0}{if ($1>0+a) a=$1} END{print a}' >> ${convergence_directory}/scores_gso_total_max 
		# Median:
		cat scores_gso_${step} | awk '{ count [NR] = $ 1 ; } END { if (NR % 2) { print count [ (NR + 1) / 2 ]; } else { print (count [ (NR / 2)] + count [ (NR / 2) + 1 ]) / 2.0 ; } }' >> scores_gso_total_median
	
	done

done

