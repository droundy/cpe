#!/usr/bin/python2
from __future__ import division

from matplotlib import pyplot as plt
import numpy as np
from scipy import optimize
import sys

import compute

e_w = 80 #Permiativity in water
e_0 = 8.85e-12 #F/m Vacuum Permiativity

nm = 1e-9
dx = 0.001*nm
xmax = 10*nm
x = np.arange(dx/2,xmax,dx) #range of distances

def resistivity_from_log_differences(log_differences):
    rho = np.zeros_like(log_differences)
    for i in xrange(len(rho)-2,-1,-1):
        if np.isfinite(np.exp(log_differences[i])):
            rho[i] = rho[i+1] + np.exp(log_differences[i])
        else:
            rho[i] = rho[i+1]
    return rho

try:
    oldrho = np.loadtxt('resistivity.csv')
    oldrho[-1] = 0 # since that is enforced elsewhere...
    log_differences = np.zeros_like(x)
    #log_differences[:-1] = np.log(oldrho[:-1] - oldrho[1:])
    for i in xrange(len(x)-1):
        log_differences[i] = np.log(oldrho[i] - oldrho[i+1])
        if not np.isfinite(log_differences[i]):
            log_differences[i] = 1e-300
except:
    print sys.exc_info()[0]
    print 'failed to read in resistivity!'
    log_differences = -10-np.log(x)

experiment = np.loadtxt('RawDataGraphene.txt')
omega = experiment[:,0]/2/np.pi
Vsd = 25e-3 # experimental applied voltage = 25 mV
experimental_impedance_real = Vsd/experiment[:,5]*15800 # scale by area
experimental_impedance_imag = Vsd/experiment[:,6]*15800 # scale by area
print experimental_impedance_real
print experimental_impedance_imag
experimental_impedance = np.concatenate((experimental_impedance_real,
                                         experimental_impedance_imag))

def impedances(log_differences):
    rho = resistivity_from_log_differences(log_differences)
    Z_real = compute.impedance_real(x[1]-x[0], rho, omega, e_0*e_w)
    Z_imag = compute.impedance_imag(x[1]-x[0], rho, omega, e_0*e_w)
    # the factors of 1e12 below convert from Ohm*meter^2 to
    # Ohm*micron^2, which is what we are plotting.
    return np.concatenate((Z_real*1e12, Z_imag*1e12))

mycount = 0
def errors(log_differences):
    global mycount
    err = sum((np.log(impedances(log_differences))
               - np.log(experimental_impedance))**2)
    if mycount % 100 == 0:
        np.savetxt('resistivity.csv',
                   resistivity_from_log_differences(log_differences))
        print 'error', mycount, err
    mycount += 1
    return err

def callme(x):
    print
    print 'hello'
    np.savetxt('resistivity.csv', resistivity_from_log_differences(x))
    print

datas = optimize.minimize(errors,
                          log_differences,
                          method='Powell', # 'Nelder-Mead'
                          callback=callme,
                          options={'disp': True,
                                   'maxiter': 100})
print datas

log_differences = datas.x
np.savetxt('resistivity.csv', resistivity_from_log_differences(log_differences))

Zs = impedances(log_differences)

plt.loglog(omega, experimental_impedance[0:len(omega)], 'ro')
plt.loglog(omega, experimental_impedance[len(omega):], 'rx')

plt.loglog(omega, Zs[0:len(omega)], 'b-')
plt.loglog(omega, Zs[len(omega):], 'b--')

rho = resistivity_from_log_differences(log_differences)

plt.ylabel(r'$Z$')
plt.xlabel(r'$\omega$ ')

plt.figure()
rhomax = 1e100
plt.loglog(x[rho < rhomax]/nm, rho[rho < rhomax], '-')
plt.xlabel(r'x ($nm$)')
plt.ylabel(r'Resistivity ($\Omega m$)')
plt.title('')


plt.show()
