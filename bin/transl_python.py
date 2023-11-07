#!/usr/bin/python

# This script translate 1 letter aminoacids into 3 letters aminoacids

import sys

def traduction(aa):

	protein_letters_1to3 = {
    	"A": "ALA",
    	"C": "CYS",
    	"D": "ASP",
    	"E": "GLU",
    	"F": "PHE",
    	"G": "GLY",
    	"H": "HIS",
    	"I": "ILE",
    	"K": "LYS",
    	"L": "LEU",
    	"M": "MET",
    	"N": "ASN",
    	"P": "PRO",
    	"Q": "GLN",
    	"R": "ARG",
    	"S": "SER",
    	"T": "THR",
    	"V": "VAL",
    	"W": "TRP",
    	"Y": "TYR",
	}
	return protein_letters_1to3.get(aa)

if __name__ == "__main__":
	aa_bash = sys.argv[1]
	result = traduction(aa_bash)
	print(result)

