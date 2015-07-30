from __future__ import division

import numpy as np
import collections

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

def potential_from_concentration(kT, n_0, n): #molar
    """Given the concentration of ions find the potential felt
       by the ions.

       kT is Boltzmann's constant time the temperature.

       n_0 is the bulk concentration of ions, i.e. the concentration
       where the potential is zero.

       n is the actual concentration of ions.

    """
    return -kT*np.log(n/n_0)

def impedance_real(dxs, resistivity, omega, epsilon): # ohm*m**2
    """Find the real part of the impedance.  The dxs argument is either a
       scalar value for dx, or an array holding the dx value for each
       resistivity.

       epsilon is the dielectric constant of the material.

       resistivity is an array with the same size as xs.

    """
    try:
        dxs[0]
    except:
        dx = dxs
    Z_real = 0
    for i in range(len(resistivity)):
        try:
            dx = dxs[i]
        except:
            pass
        dZ = dx*resistivity[i]/(1+(omega*epsilon*resistivity[i])**2)
        if np.isnan(dZ).any() or np.isinf(dZ).any():
            # The resistivity is presumably infinite... which would
            # give us a correct contribution of zero.
            pass
        else:
            Z_real += dZ
    return Z_real


def impedance_imag(dxs, resistivity, omega, epsilon): # ohm*m**2
    """Find the imaginary part of the impedance.  The dxs argument is either a
       scalar value for dx, or an array holding the dx value for each
       resistivity.

       epsilon is the dielectric constant of the material.

       resistivity is an array with the same size as xs.
    """
    try:
        dxs[0]
    except:
        dx = dxs
    Z_imag = 0
    for i in range(len(resistivity)):
        try:
            dx = dxs[i]
        except:
            pass
        dZ = dx*omega*epsilon*resistivity[i]**2/(1+(omega*epsilon*resistivity[i])**2)
        if np.isnan(dZ).any() or np.isinf(dZ).any():
            Z_imag += dx/(omega*epsilon)
        else:
            Z_imag += dZ
    return Z_imag
