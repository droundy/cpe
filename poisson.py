#!/usr/bin/python

from __future__ import division

import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt

x = np.arange(0, 1000.0, 0.001)
dx = x[1] - x[0]

epsilon = 80.0

beta = 1.0
n0 = 0.1
Lambda_cation = 1.0
Lambda_anion = 0.0
K_cation = 0.0 # Kohlrausch's constant
K_anion = 0.0 # Kohlrausch's constant

K = np.sqrt(2*beta*n0/epsilon)

def psi_from_psi0(psi0):
    y0 = beta*psi0
    return  2/beta*np.log((np.exp(y0/2)+1 + (np.exp(y0/2)-1)*np.exp(-K*x))
                          /(np.exp(y0/2)+1 - (np.exp(y0/2)-1)*np.exp(-K*x)))

def ncation_from_psi0(psi0):
    psi = psi_from_psi0(psi0)
    return n0*np.exp(-beta*psi)

def nanion_from_psi0(psi0):
    psi = psi_from_psi0(psi0)
    return n0*np.exp(beta*psi)

# the following is the charge density
def rho_from_psi0(psi0):
    return ncation_from_psi0(psi0) - nanion_from_psi0(psi0)

def resistivity_from_psi0(psi0):
    nanion = nanion_from_psi0(psi0)
    anion_conductivity = Lambda_anion*nanion - K_anion*nanion*np.sqrt(nanion)
    ncation = ncation_from_psi0(psi0)
    cation_conductivity = Lambda_cation*ncation - K_cation*ncation*np.sqrt(ncation)
    return 1.0/(anion_conductivity + cation_conductivity)

def sigma_from_psi0(psi0):
    return sum(rho_from_psi0(psi0))*dx

psi0 = 10.0
psi = psi_from_psi0(psi0)

plt.subplot(211)

plt.plot(x, psi, 'g-')
plt.ylabel(r'$\psi$')

plt.subplot(212)

plt.semilogy(x, ncation_from_psi0(psi0), 'r-')
plt.semilogy(x, nanion_from_psi0(psi0), 'b-')
plt.semilogy(x, resistivity_from_psi0(psi0), 'k:')

plt.xlabel('$x$')
plt.ylabel('$n$')

plt.savefig('poisson.pdf')

plt.figure()

resistivity = resistivity_from_psi0(psi0)

plt.loglog(x, resistivity, 'k-')
plt.ylabel(r'resistivity $\rho$')
plt.xlabel('$x$')

plt.savefig('resistivity.pdf')

plt.figure()

psi0s = np.arange(0, 10, 0.1)
sigmas = np.zeros_like(psi0s)
for i in range(len(psi0s)):
    sigmas[i] = sigma_from_psi0(psi0s[i])

plt.plot(psi0s, sigmas)
plt.xlabel(r'$\Delta\Phi$')
plt.ylabel(r'$\sigma$')

plt.savefig('sigma-vs-V.pdf')


# plt.figure()

# omegas = 10**np.arange(0.0, 5, 0.1)
# Z = np.zeros(len(omegas), dtype=complex)

# make resistivity a complex array?

# for i in range(len(omegas)):
#     Z[i] = dx*sum(1.0/(resistivity + 1.0j*omegas[i]*epsilon))

# plt.loglog(omegas, Z.real, 'r-')
# plt.loglog(omegas, Z.imag, 'b-')

# plt.savefig('impedance.pdf')
