#!/bin/bash

## This script parses the output of dpocket to choose only the fields of interest for the study. The results of this filtering are written to two new files in the folder with the input. 


# Root directory 
root="/lrlhps/users/l001803/TMP/crbn_eval/"

# Files and directories
dir="${root}results/"
file="${dir}dpout_fpocketp.txt"

cat $file | awk -v OFS="\t" '{print $1,$11,$20,$21,$22,$24,$26,$31}' > "${dir}results_dpocket.txt"

#Drug score distribution
cat $file | awk '{print $31}' > "${dir}results_dpocket_drugscore.txt"

