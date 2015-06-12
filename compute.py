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
