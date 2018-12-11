#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 14:14:50 2018

@author: cristina
"""
import numpy as np

N_atoms = 3
N_matrix = 4
lamda = 3

Self = np.zeros([N_matrix ** 2, N_matrix ** 2], dtype= int)


for i_matrix in range(N_matrix - 1):
        for j_matrix in range(N_matrix):
            
            g_i = (i_matrix)*N_matrix + j_matrix
            g_j = (i_matrix)*N_matrix +( j_matrix + 1 )
            
            Self [g_i, g_j]= lamda
            Self [g_j, g_i]= -lamda
            
i_matrix = 3
for j_matrix in range(N_matrix - 1):
            
    g_i = (i_matrix)*N_matrix + j_matrix
    g_j = (i_matrix)*N_matrix +( j_matrix + 1 )
            
    Self [g_i, g_j]= lamda
    Self [g_j, g_i]= -lamda


for i_matrix in range(N_matrix - 1):
        for j_matrix in range(N_matrix):
            
            g_i = (i_matrix)*N_matrix + j_matrix
            g_j = (i_matrix + 1)*N_matrix + j_matrix
            
            Self [g_i, g_j]= 2*lamda
            Self [g_j, g_i]= -2*lamda
            
#i_matrix = 3
#for j_matrix in range(N_matrix - 1):
#            
#    g_i = (i_matrix)*N_matrix + j_matrix
#    g_j = (i_matrix)*N_matrix +( j_matrix + 1 )
#            
#    Self [g_i, g_j]= lamda
#    Self [g_j, g_i]= -lamda
            
            
            

#for i_matrix in range(N_matrix - 1):
#    j_matrix = 3
#    g_i = (i_matrix-1)*N_matrix + j_matrix
#    g_j = (i_matrix-1)*N_matrix + (j_matrix+1)
#            
#    Self [g_i, g_j]= 3 * lamda
#    Self [g_j, g_i]= - 3 * lamda


#for j_matrix in range(N_matrix - 1):
#            
#    g_i = (i_matrix-1)*N_matrix + j_matrix
#    g_j = (i_matrix-1)*N_matrix + (j_matrix+1)
#            
#    Self [g_i, g_j]= lamda
#    Self [g_j, g_i]= -lamda

#i = np.arange(N_matrix)
#j = np.arange(N_matrix)
#I, J = np.meshgrid(i, j)
#ii = np.reshape(I, ((N_matrix) * (N_matrix), 1))
#jj = np.reshape(J, ((N_matrix) * (N_matrix), 1))
#
#g_i = (ii-1)*N_matrix + jj
#g_j = (ii-1)*N_matrix + (jj + 1)
#
#
#Self [g_i, g_j]= lamda
#Self [g_j, g_i]= -lamda
#
#n = np.arange(N_matrix)
#m = np.arange(N_matrix)
#N, M = np.meshgrid(n, m)
#nn = np.reshape(N, ((N_matrix) * (N_matrix), 1))
#mm = np.reshape(M, ((N_matrix) * (N_matrix), 1))
#
#g_i = (nn-1)*N_matrix + mm
#g_j = (nn)*N_matrix + mm
#    
#
#Self [g_i, g_j]= 2 * lamda
#Self [g_j, g_i]= -2 * lamda
#    
#

print(Self)