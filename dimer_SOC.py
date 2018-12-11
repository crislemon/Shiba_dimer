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

pi=np.pi

n=1.0

alpha=np.linspace(0,2,15)

Vpeak_AF_plus =np.zeros(len(alpha))
Vpeak_AF_minus =np.zeros(len(alpha))



for n_i in range(len(alpha)):
    
    #AF case
    (gg , N_matrix , N_omega , vv, spectro)=sc.Shiba_Chain2(n, 3, 'AF', alpha[n_i])
    
    plt.figure(1)
    plt.style.use('seaborn-bright')
    row = int(N_matrix/2)
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
    (gg , N_matrix , N_omega , vv, spectro)=sc.Shiba_Chain2(n, 3, 'FM', alpha[n_i])
    
    plt.figure(2)
    plt.style.use('seaborn-bright')
    row = int(N_matrix/2)
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
  

Delta= 1.0#meV
    
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
plt.plot(alpha,Vpeak_FM1_minus/Delta,'r.-')
plt.xlabel('alpha')
plt.ylabel('E/Delta')
plt.ylim(0, 0.3)
