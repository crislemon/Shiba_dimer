#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 14:22:20 2018

@author: cristina
"""

import numpy as np
import matplotlib.pyplot as plt
import Shiba_Chain2 as sc
import detect_peaks as dp


alpha_ini = 0.0
alpha_final = 5.0
alpha = np.linspace(alpha_ini, alpha_final, 10)
pi=np.pi
borde = 1
ancho = 3
d = 1.0
N_omega = 2001
U = 5500./27211.6
U = 0
k_f = 0.5
DOS = 1.0
s = 5.0/2.0 #spin
delta = 0.75/27211.6 #SC gap
j = 1800./27211.6 #coupling

row = int(ancho/2)


Vpeak_AF_plus =np.zeros(len(alpha))
Vpeak_AF_minus =np.zeros(len(alpha))



for n_i in range(len(alpha)):
    
    #AF case
    spin1=0
    spin2= pi
    (vv, spectro) = sc.Shiba_Chain2(d, 2, N_omega, spin1, spin2, alpha[n_i], borde, ancho, k_f, U, j, DOS, s, delta)
    
    plt.figure(1)
    plt.style.use('seaborn-bright')
    plt.plot(vv, spectro[1,row,:],linewidth=0.8, label = '%i' % n_i)
    ndexes = dp.detect_peaks(spectro[2,row,:])
    peaks = vv[ndexes]
    plt.plot(peaks,spectro[1,row,ndexes],'y*')
    plt.xlabel('mV')
    plt.ylabel('PDOS')
    plt.legend()
    plt.title('AF case')
    
    
    Shiba=vv[ndexes]
    
    Shiba_plus=[]
    Shiba_minus=[]
    
    for i in ndexes:
        if vv[i]<=0:
            Shiba_plus.append(vv[i])
        else:
            Shiba_minus.append(vv[i])
    
    Vpeak_AF_plus[n_i]= Shiba_plus[-1]     
    Vpeak_AF_minus[n_i]= Shiba_minus[0]
    
 

"Ferromagnetic spin1 = 0 spin2 = 0"


Vpeak_FM1_plus =np.zeros(len(alpha))
Vpeak_FM1_minus =np.zeros(len(alpha))
Vpeak_FM2_plus =np.zeros(len(alpha))
Vpeak_FM2_minus =np.zeros(len(alpha))

for n_i in range(len(alpha)):
    
    #FM case
    spin1=0
    spin2=0

    (vv, spectro) = sc.Shiba_Chain2(d, 2, N_omega, spin1, spin2, alpha[n_i], borde, ancho, k_f, U, j, DOS, s, delta)
    
    plt.figure(2)
    plt.style.use('seaborn-bright')
    plt.plot(vv, spectro[1,row,:],linewidth=0.8, label = '%i' % n_i)
    ndexes = dp.detect_peaks(spectro[1,row,:])
    peaks = vv[ndexes]
    plt.plot(peaks,spectro[1,row,ndexes],'y*')
    plt.xlabel('mV')
    plt.ylabel('PDOS')
    plt.legend()
    
    
    Shiba=vv[ndexes]
    
    Shiba_plus=[]
    Shiba_minus=[]
    
    for i in ndexes:
        if vv[i]<=0:
            Shiba_plus.append(vv[i])
        else:
            Shiba_minus.append(vv[i])

    Vpeak_FM1_plus[n_i]=Shiba_plus[-1]
    Vpeak_FM1_minus[n_i]=Shiba_minus[0]
    Vpeak_FM2_plus[n_i]=Shiba_plus[0]
    Vpeak_FM2_minus[n_i]=Shiba_minus[-1]
  

#Delta= 1.0#meV
    
plt.figure(3)
plt.style.use('seaborn-pastel')
plt.plot(alpha,Vpeak_AF_plus,'b.-',alpha,Vpeak_AF_minus,'b.-',label = 'AF')
plt.plot(alpha,Vpeak_FM1_plus,'r.-',alpha,Vpeak_FM1_minus,'r.-', label = 'FM')
#plt.plot(alpha,Vpeak_FM2_plus,'r.-',alpha,Vpeak_FM2_minus,'r.-', label = 'FM')
plt.show()

plt.legend()
plt.xlabel('alpha')
plt.ylabel('Shiba peak (meV)')
plt.title('SOC')


plt.figure(4)
plt.plot(alpha,Vpeak_FM1_minus/delta,'r.-')
plt.xlabel('alpha')
plt.ylabel('E/Delta')
plt.ylim(0, 0.3)
