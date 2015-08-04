#!/usr/bin/python3
from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt

qe = 1.6e-19 # Coulomb
eV = qe*1 # J
nm = 1e-9 # meters

e_w = 80 #Permiativity in water
e_0 = 8.85e-12 #F/m Vacuum Permiativity

k = 1.38e-23 #m^2kg/s^2K Boltzman Constant
T = 300 #room temp. Kelven
kT = k*T
beta = 1/kT

NA = 6.0221413e+23 # Avogadro's number
Molar = NA/1e3 # moles/liter in SI units of 1/m**3

class Species:
    def __init__(self, name, *, q, molarity,
                 molar_conductivity, external_potential):
        self.name = name
        self.molar_conductivity = molar_conductivity
        self.PE = external_potential
        self.molarity = molarity
        self.q = q

def simple_pe(x):
    """A very crude approximation to the potential energy

    FIXME This comes from staring at the output of the
    transposed_resistivity.py program, and needs to be replaced.

    """
    return 12*eV*(x-nm)/nm

Na = Species('Na',
             q=qe,
             molarity = 0.137,
             molar_conductivity=7.634,
             external_potential=simple_pe)
Cl = Species('Cl',
             q=-qe,
             molarity = Na.molarity,
             molar_conductivity=5.011,
             external_potential=simple_pe)

def charge_density_from_Phi(species, x, Phi):
    rho = np.zeros_like(x)
    for s in species:
        n_i = s.molarity*np.exp(-beta*(s.PE(x) + s.q*Phi))
        rho += s.q*n_i*Molar
    return rho

