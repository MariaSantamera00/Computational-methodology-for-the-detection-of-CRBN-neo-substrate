#!/usr/bin/python

## This script calculates the distances between the elements of the list (each element with all the others), filtering those where the distance is less than 8
	## Input: List of coordinates 
	## Returns: Elements whose distance is less than 8 (and the distance value)

from scipy.spatial.distance import cdist
import pandas as pd
import sys


# Load data
xyz = eval(sys.argv[1])

# Calculate distances between all CA
dist_matrix = cdist(xyz, xyz)

# Atoms with less than 8A of distance
lim = 8
less_lim = []

for i in range(len(dist_matrix)):
	for j in range (len(dist_matrix[i])):
		dif=abs(i-j)
		rules = [dist_matrix[i][j] < lim , dif >= 3]
		if all(rules) :
			less_lim.append((i,j,dist_matrix[i][j]))

for i in less_lim:
	print(f"{i[0]}-{i[2]}")


