# -*- coding: utf-8 -*-

import numpy as np
from scipy.linalg import inv

def Chain_Green2(Go, Self, N_matrix):
    

#    G = np.zeros([N_matrix**2,N_matrix**2,4,4])
#    G=np.array(G,np.clongdouble)
    Id = np.identity(4*N_matrix**2)

    
    matrx_inv = inv(Id -np.dot(Go,Self))
    gg = np.dot(matrx_inv , Go)
        

#    for i_matrix in range(N_matrix**2):
#        for j_matrix in range(N_matrix**2):
#            for i in range (4):
#                for j in range (4):
#                    G[i_matrix,j_matrix,i,j] = gg [(i_matrix)*4+i,(j_matrix)*4+j]



    return(gg)
