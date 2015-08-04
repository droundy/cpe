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
n_0 = 0.01 #molarity of ions (10 mM phoshpate buffer from Crosser et
           #al., with no sodium chloride or potasium.)

nm = 1e-9
micron = 1e-6
NA = 6.0221413e+23 # Avogadro's number
Molar = NA/1e3 # moles/liter in SI units of 1/m**3


resistivity = 10.0**np.arange(200, 1.0, -0.5)
rootdxs = np.zeros_like(resistivity)

experiment = np.loadtxt('RawDataGraphene.txt')
omega = experiment[:,0]/2/np.pi
Vsd = 25e-3 # experimental applied voltage = 25 mV

experimental_impedance_real = Vsd/experiment[:,6]*15800 # scale by area
experimental_impedance_imag = Vsd/experiment[:,5]*15800 # scale by area
experimental_impedance = np.concatenate((experimental_impedance_real,
                                         experimental_impedance_imag))

mymatrix = np.zeros((len(experimental_impedance), len(rootdxs)))

for i in range(len(resistivity)):
    for j in range(len(omega)):
        o = omega[j]
        Zr = compute.impedance_real(1.0, [resistivity[i]], o, e_0*e_w)/micron**2
        Zi = compute.impedance_imag(1.0, [resistivity[i]], o, e_0*e_w)/micron**2
        mymatrix[j,i] = Zi/experimental_impedance_imag[j]
        mymatrix[j + len(omega),i] = Zr/experimental_impedance_real[j]

def chisq(rootdxs):
    Zs = np.dot(mymatrix, rootdxs**2)
    return np.sum((Zs - 1)**2)

chisq_old = 1e100
chisq_new = chisq(rootdxs)
while chisq_new < chisq_old:
    print('chisq(%g) == %g' % (rootdxs[0], chisq_new))
    chisq_old = chisq_new
    rootdxs /= 2
    chisq_new = chisq(rootdxs)
    sys.stdout.flush()

class MyBounds(object):
    def __init__(self, xmin, xmax ):
        self.xmax = xmax
        self.xmin = xmin
    def __call__(self, **kwargs):
        x = kwargs["x_new"]
        tmax = bool(np.all(x <= self.xmax))
        tmin = bool(np.all(x >= self.xmin))
        return tmax and tmin

v = optimize.minimize(chisq, rootdxs,
                      bounds = [(0, None) for dx in rootdxs],
                      options = {'maxiter': 2000,
                                 'ftol': 1e-14,
                                 'gtol': 0,
                                 'maxfun': 10000000,
                                 'maxfev': 10000000,
                                 'disp': True})

print('Done optimizing!')
# print(v)
mychisq = chisq(v.x)
print('\nfinal chisq:  %g\n' % mychisq)

def x_from_rootdxs(rootdxs):
    x = np.zeros_like(rootdxs)
    x[0] = rootdxs[0]**2/2
    for i in range(1,len(x)):
        x[i] = x[i-1] + rootdxs[i-1]**2/2 + rootdxs[i]**2/2
    return x

rootdxs = v.x
x = x_from_rootdxs(rootdxs)

np.savetxt('rho.csv', np.transpose(np.array([x, rootdxs, resistivity])))

plt.figure()

Zi = compute.impedance_imag(rootdxs**2, resistivity, omega, e_0*e_w)/micron**2
Zr = compute.impedance_real(rootdxs**2, resistivity, omega, e_0*e_w)/micron**2

plt.loglog(omega*2*np.pi, Zr, 'b-', label='real')
plt.loglog(omega*2*np.pi, Zi, 'r-', label='imag')

plt.loglog(omega*2*np.pi, experimental_impedance_real, 'b+', label='exp real')
plt.loglog(omega*2*np.pi, experimental_impedance_imag, 'r+', label='exp imag')

plt.xlabel(r'$f$ (Hz)')
plt.ylabel(r'$Z\times A$ ($\Omega m^2$)')

plt.title('$\chi^2 = %g$ (in silly units)' % mychisq)

plt.legend(loc='best')

plt.savefig('optimal-Z-fit.pdf')

plt.figure(figsize=(8,9))

plt.subplot(2,1,1)

plt.loglog(x/nm, resistivity)

# plt.ylim(-1, 1e15)
# plt.xlim(1e-3, 1e4)
plt.ylabel(r'resistivity ($\Omega m$)')
# plt.legend(loc='best')

plt.subplot(2,1,2)

v_an = 0 #number anion ARTIFICIALLY SET TO ZERO TO IGNORE PHOSPHATE
v_cat = 1 #number cation
l_an = 5.011 #mS m^2/mol, limiting molar conductivity of each ion.
l_cat = 7.634 #mS m^2/mol
lim_mol_cond = v_an*l_an + v_cat*l_cat #total limiting molar conductivity

molar_concentration = 1.0/(lim_mol_cond*resistivity)

plt.loglog(x/nm, molar_concentration)

# plt.ylim(-1, 1e15)
# plt.xlim(1e-3, 1e4)
plt.xlabel(r'$x$ (nm)')
plt.ylabel(r'$n$ (mol/L)')
# plt.legend(loc='best')

plt.tight_layout()

plt.savefig('optimal-resistivity-and-concentration.pdf')

plt.figure()

k = 1.38e-23 #m^2kg/s^2K Boltzman Constant
T = 300 #room temp. Kelven
kT = k*T
eV = 1.6e-19 # J
V = compute.potential_from_concentration(kT, n_0, molar_concentration)
plt.plot(x/nm, V/eV, '-')

q = -1.6e-19 #C electron charge
e_g = 4 #Permiativity in glass
def image_charge_potential(x): #Joules
    return q**2/(16*np.pi*e_w*e_0)*(e_w-e_g)/(e_w+e_g)*1/x

# plt.plot(x/nm, image_charge_potential(x)/eV, 'r-')

Na_solvation_free_energy = 414*1e3/NA # see Horinek et al.
plt.axhline(Na_solvation_free_energy/eV, color='r', ls='--')

plt.ylim(0, 12)
plt.xlim(0, 10)
plt.xlabel(r'$x$ (nm)')
plt.ylabel(r'$V$ (eV)')
# plt.legend(loc='best')

ax2 = plt.axes([.35, .25, .5, .6])
ax2.semilogx(x/nm, V/eV, '-')

# plt.semilogx(x/nm, image_charge_potential(x)/eV, 'r-')

plt.xlabel(r'$x$ (nm)')
plt.ylabel(r'$V$ (eV)')
plt.ylim(0,12)
plt.xlim(xmin=0.01)
plt.axhline(Na_solvation_free_energy/eV, color='r', ls='--')

plt.savefig('optimal-potential.pdf')


plt.show()
