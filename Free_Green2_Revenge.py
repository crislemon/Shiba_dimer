#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

# Functions
pi = np.pi
sin = np.sin
cos = np.cos
sqrt = np.sqrt
exp = np.exp


def Free_Green2(N_matrix, lomega, Damping, Fermi_k, mass_eff, DOS_o, Delta, a_interatomic):

    G = np.zeros([N_matrix ** 2 * 4, N_matrix ** 2 * 4], dtype=complex)
    
    omega = lomega + 1j * Damping

    # Non diagonal in atom
    i = np.arange(N_matrix)
    #o = np.arange(N_omega)
    I, J, N, M = np.meshgrid(i, i, i, i)
    ii = np.reshape(I, (N_matrix ** 4, 1))
    jj = np.reshape(J, (N_matrix ** 4, 1))
    mm = np.reshape(N, (N_matrix ** 4, 1))
    nn = np.reshape(M, (N_matrix ** 4, 1))
    #oo = np.reshape(O, (N_matrix ** 4, 1))

    rr = sqrt((ii - nn) ** 2 + (jj - mm) ** 2) * a_interatomic
    rr[np.where(rr == 0)] = 1   # avoid 1 / 0 errors !!

    #omega = Romega[oo] + 1j * Damping
    SS = sqrt(Delta**2 - omega**2)
    xi = Fermi_k / (mass_eff * SS)
    factor = - pi * DOS_o * exp(-rr / xi) / (SS * Fermi_k * rr)

    g_i = (ii - 1) * N_matrix + jj
    g_j = (nn - 1) * N_matrix + mm

    G[g_i * 4 + 0, g_j * 4 + 0] = omega * sin(Fermi_k * rr) +\
     SS * cos(Fermi_k * rr) * factor
    G[g_i * 4 + 1, g_j * 4 + 1] = omega * sin(Fermi_k * rr) -\
     SS * cos(Fermi_k * rr) * factor
    G[g_i * 4 + 2, g_j * 4 + 2] = omega * sin(Fermi_k * rr) +\
     SS * cos(Fermi_k * rr) * factor
    G[g_i * 4 + 3, g_j * 4 + 3] = omega * sin(Fermi_k * rr) -\
     SS * cos(Fermi_k * rr) * factor

    G[g_i * 4 + 0, g_j * 4 + 3] = - Delta * sin(Fermi_k * rr) * factor
    G[g_i * 4 + 1, g_j * 4 + 2] = Delta * sin(Fermi_k * rr) * factor
    G[g_i * 4 + 2, g_j * 4 + 1] = Delta * sin(Fermi_k * rr) * factor
    G[g_i * 4 + 3, g_j * 4 + 0] = - Delta * sin(Fermi_k * rr) * factor

    # Diagonal in atom
    factor = - pi * DOS_o / SS

    G[g_i * 4 + 0, g_i * 4 + 0] = omega * factor
    G[g_i * 4 + 1, g_i * 4 + 1] = omega * factor
    G[g_i * 4 + 2, g_i * 4 + 2] = omega * factor
    G[g_i * 4 + 3, g_i * 4 + 3] = omega * factor
    G[g_i * 4 + 0, g_i * 4 + 3] = -Delta * factor
    G[g_i * 4 + 1, g_i * 4 + 2] = Delta * factor
    G[g_i * 4 + 2, g_i * 4 + 1] = Delta * factor
    G[g_i * 4 + 3, g_i * 4 + 0] = -Delta * factor

    return G

