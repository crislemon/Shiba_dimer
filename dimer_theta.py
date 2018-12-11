#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 16:01:34 2018

@author: cristina
"""

import numpy as np
import matplotlib.pyplot as plt
import Shiba_Chain2 as sc
import detect_peaks as dp

pi=np.pi

n= 0.7

spin2 = np.linspace(0, 2*pi, 30)

Vpeak_AF_plus =np.zeros(len(spin2))
Vpeak_AF_minus =np.zeros(len(spin2))

spin1 = 0
alpha = 0

N_matrix = 4
row = int(N_matrix/2.0)

#nstep, N_matrix, spin1, spin2, alpha


Vpeak_FM_plus =np.zeros(len(spin2))
Vpeak_FM_minus =np.zeros(len(spin2))
    

(vv, spectro)=sc.Shiba_Chain2(n, N_matrix, spin1, spin2[0], 0)
spectro = spectro[1,row,:]
indexes = dp.detect_peaks(spectro)
Shiba=vv[indexes]
    
Shiba_plus=[]
Shiba_minus=[]
    
for i in indexes:
    if vv[i]>=0:
        Shiba_plus.append(vv[i])
    else:
        Shiba_minus.append(vv[i])
        
Vpeak_FM_plus[0]= Shiba_plus[-1]     
Vpeak_FM_minus[0]= Shiba_minus[0]    

for n_i in range(1,len(spin2)):
    
    (vv, spectro)=sc.Shiba_Chain2(n, N_matrix, spin1, spin2[n_i], 0)
    spectro = spectro[1,row,:]
    #plt.plot(vv, spectro, label='AF')
    indexes = dp.detect_peaks(spectro)
    #plt.plot(vv[indexes],spectro[indexes], '*')
    Shiba=vv[indexes]
    
    Shiba_plus=[]
    Shiba_minus=[]
    
    for i in indexes:
        if vv[i]>=0:
            Shiba_plus.append(vv[i])
        else:
            Shiba_minus.append(vv[i])
    
    #we search for the peaks atminimum distance from respect the previous one
    dif = abs(Vpeak_FM_plus[n_i-1]-Shiba_plus)
    dif = np.array(dif)
    idxplus = np.argmin(dif)
    Vpeak_FM_plus[n_i] = Shiba_plus[idxplus]
    
    dif = abs(Vpeak_FM_minus[n_i-1]-Shiba_minus)
    dif = np.array(dif)
    idxminus = np.argmin(dif)
    Vpeak_FM_minus[n_i] = Shiba_minus[idxminus]
    
    
Vpeak_FM2_plus =np.zeros(len(spin2))
Vpeak_FM2_minus =np.zeros(len(spin2))
    

(vv, spectro)=sc.Shiba_Chain2(n, N_matrix, spin1, spin2[0], 0)
spectro = spectro[1,row,:]
indexes = dp.detect_peaks(spectro)
Shiba=vv[indexes]
    
Shiba_plus=[]
Shiba_minus=[]
    
for i in indexes:
    if vv[i]>=0:
        Shiba_plus.append(vv[i])
    else:
        Shiba_minus.append(vv[i])
        
Vpeak_FM2_plus[0]= Shiba_plus[-1]     
Vpeak_FM2_minus[0]= Shiba_minus[0]    

for n_i in range(1,len(spin2)):
    
    (vv, spectro)=sc.Shiba_Chain2(n, N_matrix, spin1, spin2[n_i], 0)
    spectro = spectro[1,row,:]
    indexes = dp.detect_peaks(spectro)
    Shiba=vv[indexes]
    
    Shiba_plus=[]
    Shiba_minus=[]
    
    for i in indexes:
        if vv[i]>=0:
            Shiba_plus.append(vv[i])
        else:
            Shiba_minus.append(vv[i])
    
    #we search for the peaks atminimum distance from respect the previous one
    dif = abs(Vpeak_FM2_plus[n_i-1]-Shiba_plus)
    dif = np.array(dif)
    idxplus = np.argmin(dif)
    Vpeak_FM2_plus[n_i] = Shiba_plus[idxplus]
    
    if abs(Vpeak_FM2_plus[n_i] - Vpeak_FM_plus[n_i])<0.01 and len(Shiba_plus)==2:
        Vpeak_FM2_plus[n_i] = Shiba_plus[idxplus-1]
  
    
    dif = abs(Vpeak_FM2_minus[n_i-1]-Shiba_minus)
    dif = np.array(dif)
    idxminus = np.argmin(dif)
    Vpeak_FM2_minus[n_i] = Shiba_minus[idxminus]
    
#    if Vpeak_FM2_minus[n_i] == Vpeak_FM_minus[n_i] and len(Shiba_minus)==2:
#        Vpeak_FM2_minus[n_i] = Shiba_minus[idxminus-1]
    
    
    
plt.figure(3)
plt.style.use('seaborn-pastel')
plt.plot(spin2,Vpeak_FM_plus,'.',label = 'peak 1')#,spin2,Vpeak_FM_minus,'r.', label = 'FM')
plt.plot(spin2,Vpeak_FM2_plus,'.', label = 'peak 2')#,spin2,Vpeak_FM2_minus,'r.')


plt.legend()
plt.xlabel('theta')
plt.ylabel('Shiba peak (meV)')
plt.title('Shiba peak vs theta')

plt.show()



