#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 12:37:27 2018

@author: cristina
"""
#everything in atomic units

def Shiba_Chain2(nstep, N_atoms, N_omega, spin1, spin2, alpha, borde, ancho, k_F, U, j, DOS, s, delta):

    import numpy as np
    pi=np.pi
    
    
    N_x = N_atoms + 2*borde
    N_y = ancho


    "Magnetic impurities parameters"
    S=5.0/2.0 #Cr
    #S=1.0
    
    phi = np.zeros(N_atoms)
    thetaS = [spin1, spin2]
    thetaS = np.array(thetaS)


    "Material data Bi2Pd"
    Damping=0.02/27211.6 #Dynes damping
    Delta=0.75/27211.6 #SC gap
    DOS_o=1 #Normal phase DOS
    Fermi_k= k_F
    mass_eff=1 #SC Band effective mass

    a_interatomic=nstep*3.36/0.529
    

    "spin-orbit coupling"
    lamda = (alpha/(2 * 3.36))/27.2116
    

    "we define the omega vector"
    N_delta = 3
    Romega = np.zeros([N_omega])
    Romega=np.array(Romega, np.longdouble)

    step_omega=N_delta*Delta/(N_omega-1)

    for i_omega in range(N_omega):
        Romega[i_omega] = (-N_delta/2.*Delta+(i_omega)*step_omega)
        
    vv=Romega*27211.6

    "Kondo hamiltonian"
    J = j
    
    "We calculate the Green's functions and solve Dyson eq"
    
    #impurity Hamiltonian
    import Self_Energy2D as SE2
    Self = SE2.Self_Energy(J, S, thetaS, phi, U, N_atoms, N_x, N_y, borde, lamda)
    
    
    GG = np.zeros([4 * N_y * N_x , 4 * N_y * N_x, N_omega], dtype=complex)
    
    
    for i_omega in range(N_omega):
        
        lomega = Romega[i_omega]
    
    
        #BCS Green's function
        import Free_Green_new as FG2
        Go = FG2.Free_Green(N_x, N_y, lomega, Damping, Fermi_k, mass_eff, DOS_o, Delta, a_interatomic)
    
        
    
        import Dyson as Dy
        gg = Dy.Dyson_eq(Go , Self , N_x, N_y)
        
        GG[:,:, i_omega] = gg
        
    

    
    spectro = np.zeros([N_x * N_y, N_x * N_y, N_omega])
    #spectro = np.array(spectro,np.longdouble)

    for i_atom in range(N_y):
        for j_atom in range(N_x):
            I = i_atom*N_x + j_atom

            for i_omega in range(N_omega):
             
                tr = GG[I*4 + 0, I*4 + 0, i_omega] + GG[I*4 + 1, I*4 + 1, i_omega]
                #tr = GG[I*4 + 0, I*4 + 0, i_omega] + GG[I*4 + 3, I*4 + 3, N_omega - (i_omega+1)]
                spectro[i_atom , j_atom, i_omega]= - (tr.imag)/(2*pi)
             

    
    return(vv, spectro)

