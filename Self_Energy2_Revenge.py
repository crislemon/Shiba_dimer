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
     
    i = np.arange(N_atoms)
    I, J = np.meshgrid(i, i)
    ii = np.reshape(I, (N_atoms ** 2, 1))
    jj = np.reshape(J, (N_atoms ** 2, 1))

    g_i = (ii-1)*N_matrix + jj
    g_j = (ii-1)*N_matrix + (jj + 1)

#    Self [g_i * 4 + 0, g_j * 4 + 2]= lamda
#    Self [g_i * 4 + 1, g_j * 4 + 3]= -lamda
#    Self [g_i * 4 + 2, g_j * 4 + 0]= -lamda
#    Self [g_i * 4 + 3, g_j * 4 + 1]= lamda
#
#    Self [g_j * 4 + 0, g_i * 4 + 2]= -lamda
#    Self [g_j * 4 + 1, g_i * 4 + 3]= lamda
#    Self [g_j * 4 + 2, g_i * 4 + 0]= lamda
#    Self [g_j * 4 + 3, g_i * 4 + 1]= -lamda

    Self [g_i * 4 + 0, g_j * 4 + 1]= -lamda
    Self [g_i * 4 + 1, g_j * 4 + 0]= lamda
    Self [g_i * 4 + 2, g_j * 4 + 3]= lamda
    Self [g_i * 4 + 3, g_j * 4 + 2]= -lamda

    Self [g_i * 4 + 0, g_j * 4 + 1]= lamda
    Self [g_i * 4 + 1, g_j * 4 + 0]= -lamda
    Self [g_i * 4 + 2, g_j * 4 + 3]= -lamda
    Self [g_i * 4 + 3, g_j * 4 + 2]= lamda


    #the coupling along y direction (sigma_x)
    
    g_i = (ii-1)*N_matrix + jj
    g_j = (ii)*N_matrix + jj 
    
    Self [g_i * 4 + 0, g_j * 4 + 1]= ui * lamda
    Self [g_i * 4 + 1, g_j * 4 + 0]= ui * lamda
    Self [g_i * 4 + 2, g_j * 4 + 3]= -ui * lamda
    Self [g_i * 4 + 3, g_j * 4 + 2]= -ui * lamda

    Self [g_j * 4 + 0, g_i * 4 + 1]= - ui * lamda
    Self [g_j * 4 + 1, g_i * 4 + 0]= - ui * lamda
    Self [g_j * 4 + 2, g_i * 4 + 3]= ui * lamda
    Self [g_j * 4 + 3, g_i * 4 + 2]= ui * lamda
    
    
    return(Self)







