from __future__ import division

import numpy as np

def concentration_no_bias(n_0, kT, potential, x): #molar
    """
       Find the concentration of ions given a potential, assuming no
       applied voltage.

       n_0 is the bulk concentration of ions.

       kT is Boltzmann's constant time the temperature.

       potential is a function of position that gives the potential felt by ions.

       x is the position (which may be a numpy array).
    """
    return n_0*np.exp(-potential(x)/kT)

def impedance_real(xs, resistivity, omega, epsilon): # ohm*m**2
    """
       Find the real part of the impedance.  This function assumes
       that the x values stored in xs are equally spaced.

       epsilon is the dielectric constant of the material.

       resistivity is an array with the same size as xs.
    """
    dx = xs[1] - xs[0]
    Z_real = 0
    for i in xrange(len(xs)):
        dZ = dx*resistivity[i]/(1+(omega*epsilon*resistivity[i])**2)
        if not np.isnan(dZ).any():
            Z_real += dZ
    return Z_real


def impedance_imag(xs, resistivity, omega, epsilon): # ohm*m**2
    """
       Find the imaginary part of the impedance.  This function
       assumes that the x values stored in xs are equally spaced.

       epsilon is the dielectric constant of the material.

       resistivity is an array with the same size as xs.
    """
    dx = xs[1] - xs[0]
    Z_imag = 0
    for i in xrange(len(xs)):
        dZ = dx*omega*epsilon*resistivity[i]**2/(1+(omega*epsilon*resistivity[i])**2)
        if not np.isnan(dZ).any():
            Z_imag += dZ
        else:
            Z_imag += dx/(omega*epsilon)
    return Z_imag
