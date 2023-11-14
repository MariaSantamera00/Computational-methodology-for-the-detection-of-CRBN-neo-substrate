#!/bin/bash

## This script uses the "distances.py" script to select the residues whose alpha carbons are less than 8 A and creates a graph with these connections (using the graph.py script). It then calculates the connectivity. 

# Activate conda environment
#conda activate maria

# Function to print help in terminal
print_help(){
        echo "create_graph.sh"
        echo "María Santamera Lastras 2023"
        echo -e "\nusage:create_graph.sh <root_folder> <AF_file> <output_folder>\n"
        echo -e "\troot_folder : root folder (should contain the GitHub files)"
        echo -e "\tAF_file : AlphaFold file"
        echo -e "\toutput_folder : folder in which data and results will be stored"
}



# Argument assignation
root=$1 #root="/home/l061003/TFM_MariaSantamera/"
pdb_file=$2 #pdb_file="${root}validation_set/AF_proteins_examples/AF-Q14191-F1-model_v4.pdb"
results=$3 #results="${root}graph_results/"

#Controls help message (print "print_help" function if user writes "-h" or "-help" in terminal)
if [ "$root" == "-h" ] || [ "$root" == "-help" ] ;then
        print_help
        exit
fi


#Control of arguments (show an error message if user doesn't write enough arguments)
if [ "$#" -ne 3 ]; then
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



name=$(basename "$pdb_file" | cut -f 1 -d '.')

# Get pdb file only with CA
awk '($1=="ATOM" && $3=="CA"){print $0}' ${pdb_file} > "${results}${name}_ca.pdb"

# Get the input for python script 
# Coordinates
xyz=$(echo "[")
while IFS= read -r line
do
	xyz=$(echo "${xyz}[$(echo "$line" | awk '{print $7}'),$(echo "$line" | awk '{print $8}'),$(echo $line | awk '{print $9}')]")

done < "${results}${name}_ca.pdb"

xyz=$(echo "$(echo "$xyz" | sed 's/]\[/],[/g')]")


# Compute matrix distance. Get nodes and edges file. 
python ${root}bin/distances.py "${xyz}" > "${results}${name}_edges.txt"

# Create the graph
python ${root}bin/graph.py "${results}${name}_edges.txt" 


