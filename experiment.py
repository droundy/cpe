import numpy as np

def __experiment():
    return np.loadtxt('RawDataGraphene.txt')

def frequency():
    return __experiment()[:,0]

def omega():
    return frequency()/2/np.pi

Vsd = 25e-3 # experimental applied voltage = 25 mV

def area_sqmicron(i):
    areas_sqmicron = [23400, 19200, 15800, 7500, 6100, 5000]
    return areas_sqmicron[i]

def area_sqmeter(i):
    return area_sqmicron(i)*1e-12

def impedance_real(i):
    return Vsd/__experiment()[:,2*i+2]

def impedance_imag(i):
    return Vsd/__experiment()[:,2*i+1]
