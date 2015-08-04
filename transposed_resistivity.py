#!/usr/bin/python2
from __future__ import division

import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import sys
import scipy.optimize as optimize

import compute, experiment

e_w = 80 #Permiativity in water
e_0 = 8.85e-12 #F/m Vacuum Permiativity
n_0 = 0.01 #molarity of ions (10 mM phoshpate buffer from Crosser et
           #al., with no sodium chloride or potasium.)

nm = 1e-9
micron = 1e-6
NA = 6.0221413e+23 # Avogadro's number
Molar = NA/1e3 # moles/liter in SI units of 1/m**3


resistivity = 10.0**np.arange(12, -1.0, -0.1)

omega = experiment.omega()

symbols = ['<', 'o', 's', 'v', '^', '>']
colors  = ['r', 'b', 'g', 'c', 'm', 'k']
e_impedance_real = []
e_impedance_imag = []

def chisq(rootdxs):
    Zs = np.dot(mymatrix, rootdxs**2)
    return np.sum((Zs - 1)**2)

class MyBounds(object):
    def __init__(self, xmin, xmax ):
        self.xmax = xmax
        self.xmin = xmin
    def __call__(self, **kwargs):
        x = kwargs["x_new"]
        tmax = bool(np.all(x <= self.xmax))
        tmin = bool(np.all(x >= self.xmin))
        return tmax and tmin

def x_from_rootdxs(rootdxs):
    x = np.zeros_like(rootdxs)
    x[0] = rootdxs[0]**2/2
    for i in range(1,len(x)):
        x[i] = x[i-1] + rootdxs[i-1]**2/2 + rootdxs[i]**2/2
    return x

for w in range(len(symbols)):
    area = experiment.area_sqmicron(w)
    e_impedance_real.append(experiment.impedance_real(w)*area)
    e_impedance_imag.append(experiment.impedance_imag(w)*area)
    experimental_impedance = np.concatenate((e_impedance_real[w],
                                             e_impedance_imag[w]))

    rootdxs = np.zeros_like(resistivity)
    mymatrix = np.zeros((len(experimental_impedance), len(rootdxs)))

    for i in range(len(resistivity)):
        for j in range(len(omega)):
            o = omega[j]
            Zr = compute.impedance_real(1.0, [resistivity[i]], o, e_0*e_w)/micron**2
            Zi = compute.impedance_imag(1.0, [resistivity[i]], o, e_0*e_w)/micron**2
            mymatrix[j,i] = Zi/e_impedance_imag[w][j]
            mymatrix[j + len(omega),i] = Zr/e_impedance_real[w][j]

    chisq_old = 1e100
    chisq_new = chisq(rootdxs)
    while chisq_new < chisq_old:
        print('chisq(%g) == %g' % (rootdxs[0], chisq_new))
        chisq_old = chisq_new
        rootdxs /= 2
        chisq_new = chisq(rootdxs)
        sys.stdout.flush()

    v = optimize.minimize(chisq, rootdxs,
                          bounds = [(0, None) for dx in rootdxs],
                          options = {'maxiter': 2000,
                                     'ftol': 1e-14,
                                     'gtol': 0,
                                     'maxfun': 10000000,
                                     'maxfev': 10000000,
                                     'disp': True})

    print('Done optimizing!', w)
    # print(v)
    mychisq = chisq(v.x)
    print('\nfinal chisq(%d):  %g\n' % (w, mychisq))

    rootdxs = v.x
    x = x_from_rootdxs(rootdxs)

    np.savetxt('rho.csv', np.transpose(np.array([x, rootdxs, resistivity])))

    plt.figure('impedance')

    Zi = compute.impedance_imag(rootdxs**2, resistivity, omega, e_0*e_w)/micron**2
    Zr = compute.impedance_real(rootdxs**2, resistivity, omega, e_0*e_w)/micron**2

    if w == 0:
        plt.loglog(omega*2*np.pi, Zr, colors[w]+'--', label='real')
        plt.loglog(omega*2*np.pi, Zi, colors[w]+'-', label='imag')
    else:
        plt.loglog(omega*2*np.pi, Zr, colors[w]+'--')
        plt.loglog(omega*2*np.pi, Zi, colors[w]+'-')

    plt.loglog(omega*2*np.pi, e_impedance_real[w], 'b'+symbols[w],
               markerfacecolor='none',
               markeredgecolor=colors[w])
    plt.loglog(omega*2*np.pi, e_impedance_imag[w], 'r'+symbols[w],
               markerfacecolor=colors[w],
               markeredgecolor='none',
               label='$A = %g$, $\chi^2 = %.2g$' % (area, mychisq))

    plt.xlabel(r'$f$ (Hz)')
    plt.ylabel(r'$Z\times A$ ($\Omega m^2$)')

    #plt.title('$\chi^2 = %g$ (in silly units)' % mychisq)

    plt.figure('resistivity', figsize=(8,9))

    plt.subplot(2,1,1)
    plt.loglog(x/nm, resistivity)

    plt.ylabel(r'resistivity ($\Omega m$)')
    plt.xlim(xmin=1e-2)

    plt.subplot(2,1,2)

    v_an = 0 #number anion ARTIFICIALLY SET TO ZERO TO IGNORE PHOSPHATE
    v_cat = 1 #number cation
    l_an = 5.011 #mS m^2/mol, limiting molar conductivity of each ion.
    l_cat = 7.634 #mS m^2/mol
    lim_mol_cond = v_an*l_an + v_cat*l_cat #total limiting molar conductivity

    molar_concentration = 1.0/(lim_mol_cond*resistivity)

    plt.loglog(x/nm, molar_concentration)
    plt.xlim(xmin=1e-2)

    plt.xlabel(r'$x$ (nm)')
    plt.ylabel(r'$n$ (mol/L)')

    plt.figure('potential')

    k = 1.38e-23 #m^2kg/s^2K Boltzman Constant
    T = 300 #room temp. Kelven
    kT = k*T
    eV = 1.6e-19 # J
    V = compute.potential_from_concentration(kT, n_0, molar_concentration)

    plt.plot(x/nm, V/eV, '-', label='$A=%g$' % area)

    Na_solvation_free_energy = 414*1e3/NA # see Horinek et al.
    # plt.axhline(Na_solvation_free_energy/eV, color='r', ls='--')

    plt.xlim(0, 10)
    plt.xlabel(r'$x$ (nm)')
    plt.ylabel(r'$V$ (eV)')

    # ax2 = plt.axes([.35, .25, .5, .6])
    # ax2.semilogx(x/nm, V/eV, '-')

    # plt.xlabel(r'$x$ (nm)')
    # plt.ylabel(r'$V$ (eV)')
    # plt.xlim(xmin=0.01)
    # plt.axhline(Na_solvation_free_energy/eV, color='r', ls='--')

plt.figure('impedance')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('optimal-Z-fit.pdf')

plt.figure('resistivity')
plt.tight_layout()
plt.savefig('optimal-resistivity-and-concentration.pdf')


plt.figure('potential')
plt.legend(loc='best')
plt.savefig('optimal-potential.pdf')

plt.show()
