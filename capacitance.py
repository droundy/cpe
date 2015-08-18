from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

# initial conditions 

E_0 = 0 #N/C
Phi_0 = 1 #J/C
n_0 = 1
k = 1.38e-23 #m^2kg/s^2K Boltzman Constant
T = 300 #room temp. Kelven
beta = 1/(k*T)
e_w = 80 #Permiativity in water
e_g = 4 #Permiativity in glass
e_0 = 8.85e-12 #F/m Vacuum Permiativity
eV = 1.6*1e-19 #C

nm = 1e-9
dx = 0.001*nm
xmax = 10*nm
x = np.arange(dx/2,xmax,dx) #range of distances

nop = 2 #number of particles
q = [1,-1] #charges for particles 


V_graph = (-0.1*x/nm+0.6)*eV
V_graph[V_graph < 0] = 0

n = np.zeros_like(x)
rho = np.zeros_like(x)
E = np.zeros_like(x)
Phi = np.zeros_like(x)
Phi[len(x)-2] = .001 


for j in range(len(x)-2,0,-1):
    rho[j] = 0
    for i in range(nop):
        #print n_0*np.exp(-beta*(V_graph[j] + q[i]*Phi[j]))
        n[j] = n_0*np.exp(-beta*(V_graph[j] + q[i]*Phi[j])) #Concentration 
        rho[j] += q[i]*n[j] #Charge Density 
    E[j] = rho[j]/(e_w*e_0)*dx + E[j+1]
    Phi[j-1] = E[j]*dx + Phi[j]
    
plt.figure()
plt.plot(x,V_graph)

plt.figure()
plt.plot(x,Phi)
plt.show()



