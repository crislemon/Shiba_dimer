#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 17:12:41 2018

@author: cristina
"""

import numpy as np
import matplotlib.pyplot as plt
import Shiba_Chain2 as sc
import detect_peaks as dp


d_ini = 0.5
d_final = 2.5
N = np.linspace(d_final, d_ini, 31)
pi=np.pi
borde = 1
ancho = 3
alpha = 0.0
N_omega = 2001
U = 5500./27211.6
U = 0
k_f = 0.5
DOS = 1.0
s = 5.0/2.0 #spin
delta = 0.75/27211.6 #SC gap
j = 1800./27211.6 #coupling

row = int(ancho/2)


""""""""""""""""""""""""
"Anti-ferromagnetic spin1 = 0 spin2 = pi"
""""""""""""""""""""""""
spin1=0
spin2= pi


Vpeak_AF_plus =np.zeros(len(N))
Vpeak_AF_minus =np.zeros(len(N))

(vv, spectro) = sc.Shiba_Chain2(N[0], 2, N_omega, spin1, spin2, alpha, borde, ancho, k_f, U, j, DOS, s, delta)
    
spectro = spectro[1,row,:]
indexes = dp.detect_peaks(spectro)
    
Shiba=vv[indexes]
    
Shiba_plus=[]
Shiba_minus=[]
    
for i in indexes:
    if vv[i] > 0.0:
        Shiba_plus.append(vv[i])
    elif vv[i] < 0.0:
        Shiba_minus.append(vv[i])
    elif vv[i] == 0.0:
        Shiba_plus.append(vv[i])
        Shiba_minus.append(vv[i])
    
    
Vpeak_AF_plus[0]= Shiba_plus[0]     
Vpeak_AF_minus[0]= Shiba_minus[-1]



for n_i in range(1, len(N)):
    
    (vv, spectro) = sc.Shiba_Chain2(N[n_i], 2, N_omega, spin1, spin2, alpha, borde, ancho, k_f, U, j, DOS, s, delta)
    
    spectro = spectro[1,row,:]
    indexes = dp.detect_peaks(spectro)
    
    Shiba=vv[indexes]
    
    Shiba_plus=[]
    Shiba_minus=[]
    
    for i in indexes:
        if vv[i] > 0.0:
            Shiba_plus.append(vv[i])
        elif vv[i] < 0.0:
            Shiba_minus.append(vv[i])
        elif vv[i] == 0.0:
            Shiba_plus.append(vv[i])
            Shiba_minus.append(vv[i])
    
    
    Vpeak_AF_plus[n_i]= Shiba_plus[-1]     
    Vpeak_AF_minus[n_i]= Shiba_minus[0]
    
    #we search for the peaks atminimum distance from respect the previous one
    dif = abs(Vpeak_AF_plus[n_i-1]-Shiba_plus)
    dif = np.array(dif)
    idxplus = np.argmin(dif)
    Vpeak_AF_plus[n_i] = Shiba_plus[idxplus]
    
    dif = abs(Vpeak_AF_minus[n_i-1]-Shiba_minus)
    dif = np.array(dif)
    idxminus = np.argmin(dif)
    Vpeak_AF_minus[n_i] = Shiba_minus[idxminus]
    



#################
""""""""""""""""""
"Ferromagnetic spin1 = 0 spin2 = 0"
""""""""""""""""""

Vpeak_FM_plus =np.zeros(len(N))
Vpeak_FM_minus =np.zeros(len(N))
    
spin1=0
spin2=0

(vv, spectro) = sc.Shiba_Chain2(N[0], 2, N_omega, spin1, spin2, alpha, borde, ancho, k_f, U, j, DOS, s, delta)
spectro = spectro[1,row,:]
indexes = dp.detect_peaks(spectro)
Shiba=vv[indexes]
    
Shiba_plus=[]
Shiba_minus=[]
    
for i in indexes:
    
        if vv[i] > 0.0:
            Shiba_plus.append(vv[i])
        elif vv[i] < 0.0:
            Shiba_minus.append(vv[i])
        elif vv[i] == 0.0:
            Shiba_plus.append(vv[i])
            Shiba_minus.append(vv[i])
        
Vpeak_FM_plus[0]= Shiba_plus[0]     
Vpeak_FM_minus[0]= Shiba_minus[-1]    

for n_i in range(1,len(N)):
    
    (vv, spectro) = sc.Shiba_Chain2(N[n_i], 2, N_omega, spin1, spin2, alpha, borde, ancho, k_f, U, j, DOS, s, delta)
    spectro = spectro[1,row,:]
    #plt.plot(vv, spectro, label='AF')
    indexes = dp.detect_peaks(spectro)
    #plt.plot(vv[indexes],spectro[indexes], '*')
    Shiba=vv[indexes]
    
    Shiba_plus=[]
    Shiba_minus=[]
    
    for i in indexes:
        if vv[i] > 0.0:
            Shiba_plus.append(vv[i])
        elif vv[i] < 0.0:
            Shiba_minus.append(vv[i])
        elif vv[i] == 0.0:
            Shiba_plus.append(vv[i])
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
    
    

#################
""""""""""""""""""
"Ferromagnetic spin1 = pi spin2 = pi"
""""""""""""""""""

Vpeak_FM2_plus =np.zeros(len(N))
Vpeak_FM2_minus =np.zeros(len(N))
    
spin1=pi
spin2=pi

(vv, spectro) = sc.Shiba_Chain2(N[0], 2, N_omega, spin1, spin2, alpha, borde, ancho, k_f, U, j, DOS, s, delta)
spectro = spectro[1,row,:]
indexes = dp.detect_peaks(spectro)
Shiba=vv[indexes]
    
Shiba_plus=[]
Shiba_minus=[]
    
for i in indexes:
    
    if vv[i] > 0.0:
            Shiba_plus.append(vv[i])
    elif vv[i] < 0.0:
            Shiba_minus.append(vv[i])
    elif vv[i] == 0.0:
            Shiba_plus.append(vv[i])
            Shiba_minus.append(vv[i])
        
Vpeak_FM2_plus[0]= Shiba_plus[0]     
Vpeak_FM2_minus[0]= Shiba_minus[-1]    

for n_i in range(1,len(N)):
    
    (vv, spectro) = sc.Shiba_Chain2(N[n_i], 2, N_omega, spin1, spin2, alpha, borde, ancho, k_f, U, j, DOS, s, delta)
    spectro = spectro[1,row,:]
    indexes = dp.detect_peaks(spectro)
    Shiba=vv[indexes]
    
    Shiba_plus=[]
    Shiba_minus=[]
    
    for i in indexes:
        if vv[i] > 0.0:
            Shiba_plus.append(vv[i])
        elif vv[i] < 0.0:
            Shiba_minus.append(vv[i])
        elif vv[i] == 0.0:
            Shiba_plus.append(vv[i])
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
    
    if Vpeak_FM2_minus[n_i] == Vpeak_FM_minus[n_i] and len(Shiba_minus)==2:
        Vpeak_FM2_minus[n_i] = Shiba_minus[idxminus-1]
    
    
    
    
plt.figure(1)
plt.style.use('seaborn-pastel')
plt.plot(N,Vpeak_AF_plus,'b',N,Vpeak_AF_minus,'b',label = 'AF', linewidth=1.0)
plt.plot(N,Vpeak_FM_plus,'ro',N,Vpeak_FM_minus,'ro', label = 'FM', linewidth=0.8)
plt.plot(N,Vpeak_FM2_plus,'ro',N,Vpeak_FM2_minus,'ro', linewidth=0.8)
plt.show()

plt.legend()
plt.xlabel('distance (a)')
plt.ylabel('Shiba peak (meV)')
plt.title('U = 0 2D matrix')
plt.savefig('results/AF_FMoscillations.pdf')

  
