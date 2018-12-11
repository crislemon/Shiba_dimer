#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 17:53:46 2018

@author: cristina
"""

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import Shiba_Chain2 as sc2
import detect_peaks as dp
import time

pi=np.pi

d = 1.0
N_matrix = 4

#nstep, N_matrix, spin1, spin2, alpha

spin1 = 0
spin2 = 0

t1=time.time()

(vv, spectro, N_omega, gg) = sc2.Shiba_Chain2(d, N_matrix, spin1, spin2, 0)

t2 = time.time()

print('The program is finished after', t2 - t1)



spectro = np.zeros([N_matrix,N_matrix,N_omega])
spectro=np.array(spectro,np.longdouble)

for i_atom in range(N_matrix):
    for j_atom in range(N_matrix):
         I = i_atom + (j_atom)*N_matrix

         for i_omega in range(N_omega):
             
             tr2 = gg[I*4 + 0, I*4 + 0, i_omega] + gg[I*4 + 1, I*4 + 1, i_omega]
             spectro[i_atom,j_atom,i_omega]= - (tr2.imag)/(2*pi)


plt.figure(5)
plt.style.use('seaborn-bright')
row = int(N_matrix/2)
plt.plot(vv, spectro[1,row,:],label='%s atom' % i_atom,linewidth=0.8)
ndexes = dp.detect_peaks(spectro[1,row,:])
peaks = vv[ndexes]

minpeak = min(abs(peaks))
peaks = peaks.tolist()
i=peaks.index(minpeak)


plt.plot(peaks,spectro[1,row,ndexes],'y*')
plt.xlabel('mV')
plt.ylabel('PDOS')
plt.title('We use peak # %i ' %i)
plt.savefig('spectro.pdf')

z = np.zeros([N_matrix,N_matrix])
z = np.array(z,np.float)



#spectra = spectro[:,:,ndexes[i+1]]
#titulo = vv[ndexes[i+1]]

spectra = spectro[:,:,ndexes[i]]
titulo = vv[ndexes[i]]


for i_atom in range(N_matrix):
    for j_atom in range(N_matrix):

        z[i_atom,j_atom] = spectra[i_atom,j_atom]

#save data
data_spectro = np.array([vv,spectro[1,row,:]])
data_3D = z
np.savetxt('data_spectro.txt', data_spectro)
np.savetxt('data_3D.txt', data_3D)


# Plot the surface.
X = list(range(N_matrix))
Y = list(range (N_matrix))
X, Y = np.meshgrid(X, Y)
Z = z

#3D plot
fig2 = plt.figure(6)
ax = fig2.add_subplot((111), projection='3d')
ax.plot_wireframe(X, Y, z)
plt.title('FM omega = %f mV' %titulo)
ax.set_zlabel('PDOS')
plt.savefig('3D.pdf')

#2D plot
plt.figure(7)
plt.imshow(z)
plt.colorbar()
plt.title('FM omega = %f mV' %titulo)
plt.savefig('2D.pdf')

