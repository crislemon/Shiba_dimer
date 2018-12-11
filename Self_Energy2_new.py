#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
sin = np.sin
cos = np.cos
exp=np.exp
ui = complex(0.0, 1.0)

#Local self-energy for potential plus exchange and SOC

def Self_Energy2(J,S,thetaS,phi,U,N_atoms,N_matrix,lamda):


    Self = np.zeros([N_matrix ** 2 * 4, N_matrix ** 2 * 4], dtype=complex)
    Self2 = np.zeros([N_matrix ** 2 * 4, N_matrix ** 2 * 4], dtype=complex)
    
    i = np.arange(N_atoms)
    I = np.meshgrid(i)
    ii = np.reshape(I, (N_atoms, 1))
    
    i_matrix = int(N_matrix/2)
    
    g_i = (i_matrix) * N_matrix + (ii+1)
    
    
    theta_i=thetaS[ii]
    phi_i=phi[ii]

    "diagonal in the atom space"

    Self [g_i * 4 + 0, g_i * 4 + 0]=J*S*cos(theta_i)-U
    Self [g_i * 4 + 1, g_i * 4 + 1]=-J*S*cos(theta_i)-U
    Self [g_i * 4 + 0, g_i * 4 + 1]=J*S*sin(theta_i)*exp(-ui*phi_i)
    Self [g_i * 4 + 1, g_i * 4 + 0]=J*S*sin(theta_i)*exp(ui*phi_i)
    Self [g_i * 4 + 2, g_i * 4 + 2]=-J*S*cos(theta_i)+U
    Self [g_i * 4 + 3, g_i * 4 + 3]=J*S*cos(theta_i)+U
    Self [g_i * 4 + 2, g_i * 4 + 3]=-J*S*sin(theta_i)*exp(ui*phi_i)
    Self [g_i * 4 + 3, g_i * 4 + 2]=-J*S*sin(theta_i)*exp(-ui*phi_i)


    "Non - diagonal in the atom space"
    "Spin orbit interaction"

    #the coupling along x direction (sigma_y)
    
    se = np.zeros([N_matrix ** 2, N_matrix ** 2, 4, 4], dtype=complex)
     
    for i_matrix in range(N_matrix - 1):
        for j_matrix in range(N_matrix):
            
            g_i = (i_matrix)*N_matrix + j_matrix
            g_j = (i_matrix)*N_matrix + (j_matrix + 1)
            
            se[g_i, g_j, 0, 1]= - lamda
            se[g_i, g_j, 1, 0]= lamda
            se[g_i, g_j, 2, 3]= lamda
            se[g_i, g_j, 3, 2]= - lamda
            
            
            se[g_j, g_i, 0, 1]= lamda
            se[g_j, g_i, 1, 0]= - lamda
            se[g_j, g_i, 2, 3]= - lamda
            se[g_j, g_i, 3, 2]= lamda

    i_matrix = N_matrix - 1
    for j_matrix in range(N_matrix - 1):
            
        g_i = (i_matrix)*N_matrix + j_matrix
        g_j = (i_matrix)*N_matrix +( j_matrix + 1 )
            
        se[g_i, g_j, 0, 1]= - lamda
        se[g_i, g_j, 1, 0]= lamda
        se[g_i, g_j, 2, 3]= lamda
        se[g_i, g_j, 3, 2]= - lamda
            
            
        se[g_j, g_i, 0, 1]= lamda
        se[g_j, g_i, 1, 0]= - lamda
        se[g_j, g_i, 2, 3]= - lamda
        se[g_j, g_i, 3, 2]= lamda
    
    #the coupling along y direction (sigma_x)
    
    for i_matrix in range(N_matrix - 1):
        for j_matrix in range(N_matrix):
            
            g_i = (i_matrix)*N_matrix + j_matrix
            g_j = (i_matrix + 1)*N_matrix + j_matrix
            
            se [g_i, g_j, 0, 1]= ui * lamda
            se [g_i, g_j, 1, 0]= ui * lamda
            se [g_i, g_j, 2, 3]= -ui * lamda
            se [g_i, g_j, 3, 2]= -ui * lamda

            se [g_j, g_i, 0, 1]= -ui * lamda
            se [g_j, g_i, 1, 0]= -ui * lamda
            se [g_j, g_i, 2, 3]= ui * lamda
            se [g_j, g_i, 3, 2]= ui * lamda
    
    
    
    
    for i_matrix in range(N_matrix**2):
            for j_matrix in range(N_matrix**2):
                for i in range (4):
                    for t in range (4):
                        Self2[i_matrix*4+i,j_matrix*4+t] = se[i_matrix,j_matrix,i,t]
    
    
    return(Self + Self2)







