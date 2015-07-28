#!/usr/bin/python2
from __future__ import division

import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import sys
import scipy.optimize as optimize

import compute

e_w = 80 #Permiativity in water
e_0 = 8.85e-12 #F/m Vacuum Permiativity

nm = 1e-9
resistivity = 10.0**np.arange(200, 0.0, -0.5)
dxs = np.zeros_like(resistivity) + 1.0*nm

experiment = np.loadtxt('RawDataGraphene.txt')
omega = experiment[:,0]/2/np.pi
Vsd = 25e-3 # experimental applied voltage = 25 mV
experimental_impedance_real = Vsd/experiment[:,5]*15800 # scale by area
experimental_impedance_imag = Vsd/experiment[:,6]*15800 # scale by area
experimental_impedance = np.concatenate((experimental_impedance_real,
                                         experimental_impedance_imag))

mymatrix = np.zeros((len(experimental_impedance), len(dxs)))

for i in range(len(resistivity)):
    for j in range(len(omega)):
        o = omega[j]
        Zr = compute.impedance_real(1.0, [resistivity[i]], o, e_0*e_w)
        Zi = compute.impedance_imag(1.0, [resistivity[i]], o, e_0*e_w)
        # FIXME THE FOLLOWING ARE SWAPPED.  DO I HAVE EXPERIMENT BACKWARDS?
        mymatrix[j,i] = Zi/experimental_impedance_real[j]
        mymatrix[j + len(omega),i] = Zr/experimental_impedance_imag[j]

def chisq(dxs):
    Zs = np.dot(mymatrix, dxs)
    return np.sum((Zs - 1)**2)

v = optimize.minimize(chisq, dxs,
                      bounds = [(0, None) for dx in dxs],
                      options = {'maxiter': 2000,
                                 'maxfun': 10000000,
                                 'maxfev': 10000000,
                                 'disp': True})

print('Done optimizing!')
# print(v)

def x_from_dxs(dxs):
    x = np.zeros_like(dxs)
    x[0] = dxs[0]/2
    for i in range(1,len(x)):
        x[i] = x[i-1] + dxs[i-1]/2 + dxs[i]/2
    return x

dxs = v.x
x = x_from_dxs(dxs)

np.savetxt('rho.csv', np.transpose(np.array([x, dxs, resistivity])))

plt.figure()

plt.loglog(x, resistivity)

# plt.ylim(-1, 1e15)
# plt.xlim(1e-3, 1e4)
plt.xlabel(r'$x$ (m) ?!')
plt.ylabel(r'$\rho$ ($\Omega m$)')
# plt.legend(loc='best')

plt.savefig('optimal-resistivity.pdf')

plt.figure()

Zi = compute.impedance_imag(dxs, resistivity, omega, e_0*e_w)
Zr = compute.impedance_real(dxs, resistivity, omega, e_0*e_w)

plt.loglog(omega, Zr, 'b-', label='real')
plt.loglog(omega, Zi, 'r-', label='imag')

plt.loglog(omega, experimental_impedance_real, 'b+', label='exp real')
plt.loglog(omega, experimental_impedance_imag, 'r+', label='exp imag')

plt.xlabel(r'$\omega$ (rad/sec)')
plt.ylabel(r'$Z\times A$ ($\Omega \mu m^2$)')

plt.legend(loc='best')

plt.savefig('optimal-Z-fit.pdf')

plt.show()
