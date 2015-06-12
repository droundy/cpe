from __future__ import division

import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt

import compute

#Concentration of ions as function of distance
#
#Constants
k = 1.38e-23 #m^2kg/s^2K Boltzman Constant
T = 300 #room temp. Kelven
beta = 1/(k*T)
q = -1.6e-19 #C electron charge
e_w = 80 #Permiativity in water
e_g = 4 #Permiativity in glass
e_0 = 8.85e-12 #F/m Vacuum Permiativity
n_0 = 1 #molar number or ions
nm = 1e-5 #meters

xmin = 1e-12
x = np.arange(xmin,1e-9,1e-13) #range of distances

#Potential energy
def potential(x): #Joules
    #return q**2/(16*pi*e_w*e_0)*(e_w-e_g)/(e_w+e_g)*(xmin/(x-xmin))**10
    return q**2/(16*np.pi*e_w*e_0)*(e_w-e_g)/(e_w+e_g)*1/x

plt.figure('Potential')
plt.loglog(x,potential(x),'-')
plt.xlabel(r'Distance from graphene($nm$)')
plt.ylabel(r'Potential($J$)')
plt.title('Potential Energy')

plt.savefig('potential-energy-vs-distance.pdf')

plt.figure('Concentration')
plt.semilogx(x,compute.concentration_no_bias(n_0, k*T, potential, x), '-')
plt.xlabel(r'Distance from graphene($nm$)')
plt.ylabel(r'Concentration of ions ($\frac{mol}{m^3}$)')
plt.title('Concentration as function of distance form graphene')

plt.savefig('concentration-vs-distance.pdf')


#Conductivity

def conductivity(x): #1e-3/ohm/m
    v_an = 1 #number anion
    v_cat = 1 #number cation
    l_an = 5.011 #mS m^2/mol, limiting molar conductivity of each ion.
    l_cat = 7.634 #mS m^2/mol
    lim_mol_cond = v_an*l_an + v_cat*l_cat #total limiting molar conductivity
    return lim_mol_cond*compute.concentration_no_bias(n_0, k*T, potential, x)



#Resistivity
def resistivity(x): #1e3*ohm*m
    return 1/conductivity(x)


plt.figure('Resistivity')
plt.loglog(x/nm, resistivity(x), '-')
plt.xlabel(r'x ($nm$)')
plt.ylabel(r'Resistivity ($\Omega m$)')
plt.title('')

plt.savefig('resistivity-vs-distance.pdf')


#Impedance
omega = np.arange(1, 1e3, .5)/2/np.pi

Z_real = compute.impedance_real(x, resistivity(x), omega, e_0*e_w)

print 'We are half-way done!'

Z_imag = compute.impedance_imag(x, resistivity(x), omega, e_0*e_w)


alpha = 0.9 #fits finding slope of Impedance curves

plt.figure('Impedance')
plt.loglog(omega*2*np.pi,Z_real*1e12, 'r-', label='Real Component')
plt.loglog(omega*2*np.pi,Z_imag*1e12, 'b-', label='Imaginary Component')
plt.loglog(omega*2*np.pi,omega**(-alpha)*6e-4*1e12, 'r:', label='Real Fit')
plt.loglog(omega*2*np.pi,omega**(-alpha)*8e-3*1e12, 'b:', label='Imaginary Fit')
plt.xlabel( r'Frequency ($Hz$)')
plt.ylabel(r'Impedance ($\Omega \mu m^2$)')
plt.legend()

plt.savefig('impedance-vs-frequency.pdf')

plt.figure('Experimental impedance')

f = np.arange(1, 1000, 0.1)

experimental_data = [(9e7, 3e7, .8, .8, 23400, '-'),
                     (2e8, 9e7, .8, .8, 19200, ':'),
                     (2e8, 7e7, .7, .7, 15800, '--'),
                     (2.3e8, 1e8, .75, .75, 7500, ':'),
                     (2e8, 1.5e8, .75, .75, 6100, ':'),
                     (3e8, 1e8, .75, .75, 5000, '--'),
                     ]

for im1, re1, alphare, alphaim, area, symbol in experimental_data:
    plt.loglog(f, area*re1/f**alphare, 'r'+symbol, label=r'$A = %g \mu m^2$' % area)
    plt.loglog(f, area*im1/f**alphaim, 'b'+symbol)

plt.loglog(omega*2*np.pi,Z_real*1e12, 'r-', label='Real Component')
plt.loglog(omega*2*np.pi,Z_imag*1e12, 'b-', label='Imaginary Component')

plt.xlabel( r'Frequency ($Hz$)')
plt.ylabel(r'Impedance ($\Omega \mu m^2$)')
plt.legend().get_frame().set_alpha(0.25)

plt.savefig('experimental-impedance-vs-frequency.pdf')



plt.show()
