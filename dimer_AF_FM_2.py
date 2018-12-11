#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 16:23:26 2018

@author: cristina
"""

import numpy as np
import matplotlib.pyplot as plt
import Shiba_Chain2 as sc
#import detect_peaks as dp

#distance vector
d_ini = 0.4
d_final = 2.5
N = np.linspace(d_final, d_ini, 31)
pi=np.pi
borde = 1
ancho = 3
alpha = 3.5
N_omega = 2001
U = 5500./27211.6
#U = 0
k_f = 0.5
DOS = 1.0
s = 5.0/2.0 #spin
delta = 0.75/27211.6 #SC gap
j = 1800./27211.6 #coupling

row = int(ancho/2)

spectro_AF = np.zeros([len(N), N_omega ])
spectro_FM = np.zeros([len(N), N_omega ])
spectro_FM2 = np.zeros([len(N), N_omega ])


#######################
""""""""""""""""""""""""
"Anti-ferromagnetic spin1 = 0 spin2 = pi"
""""""""""""""""""""""""
spin1=0
spin2= pi


for n_i in range(len(N)):
    
    (vv, spectro) = sc.Shiba_Chain2(N[n_i], 2, N_omega, spin1, spin2, alpha, borde, ancho, k_f, U, j, DOS, s, delta)
    spectro_AF[n_i, :] = spectro[1,row,:]
    
#################
""""""""""""""""""
"Ferromagnetic spin1 = 0 spin2 = 0"
""""""""""""""""""
spin1 = 0
spin2 = 0

for n_i in range(len(N)):
    
    (vv, spectro) = sc.Shiba_Chain2(N[n_i], 2, N_omega, spin1, spin2, alpha, borde, ancho, k_f, U, j, DOS, s, delta)
    spectro_FM[n_i, : ] = spectro[1,row,:]
    

#################
""""""""""""""""""
"Ferromagnetic spin1 = pi spin2 = pi"
""""""""""""""""""

spin1=pi
spin2=pi

for n_i in range(len(N)):
    
    (vv, spectro) = sc.Shiba_Chain2(N[n_i], 2, N_omega, spin1, spin2, alpha, borde, ancho, k_f, U, j, DOS, s, delta)
    spectro_FM2[n_i, :] = spectro[1,row,:]
    
    
   
plt.figure(1)
plt.imshow(spectro_AF, aspect='auto', cmap = plt.cm.gnuplot)

ticks = np.linspace(0, N_omega - 1, 3, dtype = 'int')
ticklabels = vv[ticks]
ticklabels = np.around(ticklabels, decimals=2)
plt.xticks(ticks, ticklabels)

ticks2 = np.linspace(0, len(N) - 1, 3, dtype = 'int')
ticklabels2 = N[ticks2]
ticklabels2 = np.around(ticklabels2, decimals=2)
plt.yticks(ticks2, ticklabels2)

plt.xlabel('Energy (meV)')
plt.ylabel('d (a)')
plt.title('AF dimer SOC 2.5 eV')
plt.colorbar()
plt.savefig('results/AF_dimer.pdf')

#FM up
plt.figure(2)
plt.imshow(spectro_FM, aspect='auto', cmap = plt.cm.gnuplot)

ticks = np.linspace(0,N_omega - 1, 3, dtype = 'int')
ticklabels = vv[ticks]
ticklabels = np.around(ticklabels, decimals=2)
plt.xticks(ticks, ticklabels)

ticks2 = np.linspace(0, len(N) - 1, 3, dtype = 'int')
ticklabels2 = N[ticks2]
ticklabels2 = np.around(ticklabels2, decimals=2)
plt.yticks(ticks2, ticklabels2)

plt.xlabel('Energy (meV)')
plt.ylabel('d (a)')
plt.title('FM up dimer SOC 2.5 eV')
plt.colorbar()
plt.savefig('results/FM_dimer.pdf')

#FM down
plt.figure(3)
plt.imshow(spectro_FM2, aspect='auto', cmap = plt.cm.gnuplot)

ticks = np.linspace(0, N_omega - 1, 3, dtype = 'int')
ticklabels = vv[ticks]
ticklabels = np.around(ticklabels, decimals=2)
plt.xticks(ticks, ticklabels)

ticks2 = np.linspace(0, len(N) - 1, 3, dtype = 'int')
ticklabels2 = N[ticks2]
ticklabels2 = np.around(ticklabels2, decimals=2)
plt.yticks(ticks2, ticklabels2)

plt.xlabel('Energy (meV)')
plt.ylabel('d (a)')
plt.title('FM down dimer SOC 2.5 eV')
plt.colorbar()
plt.savefig('results/FM2_dimer.pdf')

#save data
np.savetxt('results/spectro_AF.txt', spectro_AF)
np.savetxt('results/spectro_FM.txt', spectro_FM)
np.savetxt('resuts/spectro_FM2.txt', spectro_FM2)


  
