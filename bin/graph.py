#!/usr/bin/python

## This script creates a Graph using networkx and calculates connected components
	## Input: List of linked nodes (node1-node2...)
	## Returns: Sorted connected components of the created graph

import sys
import networkx as nx

# Load data
edges = sys.argv[1]

# Create an empty graph
G = nx.Graph()

with open(sys.argv[1], 'r') as file:
	for line in file:
		node1, node2 = line.strip().split('-')
		G.add_edge(node1, node2)

#print(list(G.edges)

#Conected components
ccomp=list(nx.connected_components(G))
for i in range(len(ccomp)):
	print(f"Componente {i}:")
	print(sorted(ccomp[i],key=lambda x: int(x)))

