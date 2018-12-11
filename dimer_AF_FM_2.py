#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 16:23:26 2018

@author: cristina
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 17:12:41 2018

@author: cristina
"""

#dimer

import numpy as np
import matplotlib.pyplot as plt
import Shiba_Chain2 as sc
import detect_peaks as dp

#distance vector
d_ini = 0.4
d_final = 2.5
N = np.linspace(d_final, d_ini, 21)
pi=np.pi
borde = 1
ancho = 3
alpha = 0.0
N_omega = 2001
U = 5500./27211.6
#U = 0
k_F = 0.183

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
    
    (vv, spectro) = sc.Shiba_Chain2(N[n_i], 2, borde, ancho, spin1, spin2, alpha, N_omega, U, k_F)
    spectro_AF[n_i, :] = spectro[1,row,:]
    
#################
""""""""""""""""""
"Ferromagnetic spin1 = 0 spin2 = 0"
""""""""""""""""""
spin1 = 0
spin2 = 0

for n_i in range(len(N)):
    
    (vv, spectro) = sc.Shiba_Chain2(N[n_i], 2, borde, ancho, spin1, spin2, alpha, N_omega, U, k_F)
    spectro_FM[n_i, : ] = spectro[1,row,:]
    

#################
""""""""""""""""""
"Ferromagnetic spin1 = pi spin2 = pi"
""""""""""""""""""

spin1=pi
spin2=pi

for n_i in range(len(N)):
    
    (vv, spectro) = sc.Shiba_Chain2(N[n_i], 2, borde, ancho, spin1, spin2, alpha, N_omega, U, k_F)
    spectro_FM2[n_i, :] = spectro[1,row,:]
    
    
#    if Vpeak_FM2_minus[n_i] == Vpeak_FM_minus[n_i] and len(Shiba_minus)==2:
#        Vpeak_FM2_minus[n_i] = Shiba_minus[idxminus-1]
    
    
    
    
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
plt.savefig('AF_dimer.pdf')

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
plt.savefig('FM_dimer.pdf')

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
plt.savefig('FM2_dimer.pdf')

#save data
np.savetxt('spectro_AF.txt', spectro_AF)
np.savetxt('spectro_FM.txt', spectro_FM)
np.savetxt('spectro_FM2.txt', spectro_FM2)


  
