#!/usr/bin/python

## This script calculates the distances between the elements of the lists (each element of xyz_lig with each element of xyz_cer), filtering those where the distance is less than 3
        ## Input: Two list of coordinates
        ## Returns: Elements whose distance is less than 3 


from scipy.spatial.distance import cdist
import pandas as pd
import sys


# Load data
xyz_lig = eval(sys.argv[1])
xyz_cer = eval(sys.argv[2])

# Calculate distances between all CA
dist_matrix = cdist(xyz_lig, xyz_cer)

# Atoms with less than 8A of distance
lim = 3
less_lim = []


for i in range(len(dist_matrix)):
	for j in range (len(dist_matrix[i])):
		rules = [dist_matrix[i][j] < lim]
		if all(rules) :
			less_lim.append((i,j,dist_matrix[i][j]))

for i in less_lim:
	print(f"{i[0]}")


