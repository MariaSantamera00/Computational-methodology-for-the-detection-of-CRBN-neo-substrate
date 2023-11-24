#!/usr/bin/python

## This script calculates the distances between the elements of the list (each element with all the others)
        ## Input: List of coordinates
        ## Returns: Distance values


from scipy.spatial.distance import cdist
import pandas as pd
import sys


# Load data
file_lig = sys.argv[1]
file_cer = sys.argv[2]

with open(file_lig, 'r') as file1:
	xyz_lig = eval(file1.read())

with open(file_cer, 'r') as file2:
	xyz_cer = eval(file2.read())



# Calculate distances between all CA
dist_matrix = cdist(xyz_lig, xyz_cer)


# Write distances

for i in range(len(dist_matrix)):
	for j in range (len(dist_matrix[i])):
		print(f"{dist_matrix[i][j]}")

