#!/bin/bash

## This script allows to obtain molecular glues of the aligned pdbs and overlap Phe with the molecular glue of the pdb file. 

## OVERLAP PHE AND MOLECULAR GLUE
obabel -ipdb Documents/phe.pdb -omol2 -O Documents/phe.mol2 -h
obabel -ipdb 5fqd_mglue.pdb -omol2 -O 5fqd_mglue.mol2 -h
/lrlhps/users/l001803/TMP/CROCK.exe -r /home/l061003/Documents/pdb_files/pdb_alineado/5fqd_mglue.mol2 -d /home/l061003/Documents/phe.mol2 -o phe_5fqd
obabel -imol2 /home/l061003/phe/phe_5fqd.mol2 -opdb -O /home/l061003/phe/phe_5fqd.pdb

cat phe_5fqd.pdb | grep "ATOM" | awk '($1=="ATOM") {$0 = substr($0, 1, 21) "R" substr($0,23)} 1' > phe_5fqd_r.pdb


# OBTAIN MOLECULAR GLUES
cat 6uml_al.pdb | grep "Y70 C" > 6uml_mglue.pdb
cat 7lps_al.pdb | grep "RN9 B" > 7lps_mglue.pdb
cat 6bn7_al.pdb | grep "RN3 B" > 6bn7_mglue.pdb
cat 6h0f_al.pdb | grep "Y70 B " > 6h0f_mglue.pdb
cat 6h0g_al.pdb | grep "Y70 B " > 6h0g_mglue.pdb
cat 5hxb_al.pdb | grep "85C C " > 5hxb_mglue.pdb
