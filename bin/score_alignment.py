## This script that calculate a score to assess protein structure similarity following the guidelines set forth in the article: "Guyon, F., & Tufféry, P. (2014). Fast protein fragment similarity scoring using a Binet-Cauchy kernel. Bioinformatics (Oxford, England), 30(6), 784–791. https://doi.org/10.1093/bioinformatics/btt618"



import sys
import numpy as np
import math


from Bio.PDB.PDBParser import PDBParser
 

def bc_score(st1, st2):

    c1 = []
    c2 = []
    for at1 in st1.get_atoms():
        at_name = at1.get_name()
        if at_name == 'CA':
            c1.append(at1.get_coord())

    for at2 in st2.get_atoms():
        at_name = at2.get_name()
        if at_name == 'CA':
            c2.append(at2.get_coord())


    X = np.array(c1)
    Xt = np.transpose(X)

    max_res = X.shape[0]

    Y = np.array(c2)

    Y = Y[0:max_res,:]
    Yt = np.transpose(Y)


    up = np.linalg.det(np.matmul(Xt,Y))
    down1 = np.linalg.det(np.matmul(Xt,X))
    down2 = np.linalg.det(np.matmul(Yt,Y))

    down = math.sqrt(down1*down2)


    return up/down

pp = PDBParser()

 

st1 = pp.get_structure('prot1',sys.argv[1])
st2 = pp.get_structure('prot1',sys.argv[2])
 

print(bc_score(st1,st2))
