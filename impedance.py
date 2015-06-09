from __future__ import division 
from pylab import *
import matplotlib.pyplot as plt

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

xmin = 1e-10
x = arange(xmin,1e-9,1e-12) #range of distances 

#Potential energy
def potential(x): #Joules
    return q**2/(16*pi*e_w*e_0)*(e_w-e_g)/(e_w+e_g)*(xmin/(x-xmin))**10
    return q**2/(16*pi*e_w*e_0)*(e_w-e_g)/(e_w+e_g)*1/x

#Concentration
def concentration(x): #molar
    return n_0*exp(-beta*potential(x))


plt.figure('Potential')
plt.loglog(x,potential(x),'-')
plt.xlabel(r'Distance from graphene($nm$)')
plt.ylabel(r'Potential($J$)')
plt.title('Potential Energy')

plt.figure('Concentration')
plt.semilogx(x,concentration(x), '-')
plt.xlabel(r'Distance from graphene($nm$)')
plt.ylabel(r'Concentration of ions ($\frac{mol}{m^3}$)')
plt.title('Concetration as function of distance form graphene')



#Conductivity

def conductivity(x): #1e-3/ohm/m
    v_an = 1 #number anion
    v_cat = 1 #number cation
    l_an = 5.011 #mS m^2/mol, limiting molar conductivity of each ion. 
    l_cat = 7.634 #mS m^2/mol
    lim_mol_cond = v_an*l_an + v_cat*l_cat #total limiting molar conductivity
    return lim_mol_cond*concentration(x) 



#Resistivity
def resistivity(x): #1e3*ohm*m
    return 1/conductivity(x)


plt.figure('Resistivity')
plt.loglog(x/nm, resistivity(x), '-')
plt.xlim(xmax=10)
plt.xlabel(r'x ($nm$)')
plt.ylabel(r'Resistivity ($\Omega m$)')
plt.title('')
plt.show()


#Impedance
omega = arange(1, 1e3, .5)/2/pi 
def initial_r(x,omega): #ohms*m
    return resistivity(x)/(1+(omega*e_0*e_w*resistivity(x))**2)

Z_real = 0    
for xval in x: 
    dx = .1*nm
    Z_real += initial_r(xval,omega)*dx #ohms*m**2
    
print 'hello'

def initial_i(x,omega): #ohms*m
    return (omega*e_0*e_w*resistivity(x)**2)/(1+(omega*e_0*e_w*resistivity(x))**2)
    
Z_imag = xmin/(omega*e_0*e_w)
for xval in x:
    dx = .1*nm 
    Z_imag += initial_i(xval,omega)*dx #ohms*m**2

alpha = 0.9 #fits finding slope of Impedance curves

plt.figure('Impedance')
plt.loglog(omega*2*pi,Z_real, label='Real Component')
plt.loglog(omega*2*pi,Z_imag, label='Imaginary Component')
plt.loglog(omega*2*pi,omega**(-alpha)*6e4, ':', label='Real Fit')
plt.loglog(omega*2*pi,omega**(-alpha)*8e5, ':', label='Imaginary Fit')
plt.xlabel( r'Frequency ($Hz$)')
plt.ylabel(r'Impedance ($\Omega m^2$)')
plt.legend()
plt.show()
