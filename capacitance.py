from __future__ import division

import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt

# initial conditions

E_0 = 0 #N/C
Phi_0 = 1 #J/C
k = 1.38e-23 #m^2kg/s^2K Boltzman Constant
T = 300 #room temp. Kelven
beta = 1/(k*T)
e_w = 80 #Permiativity in water
e_g = 4 #Permiativity in glass
e_0 = 8.85e-12 #F/m Vacuum Permiativity
eV = 1.6*1e-19 #J
N_avogadro = 6.02214129e23
n_0 = 1*N_avogadro/1e3

nm = 1e-9
dx = 0.1*nm
xmax = 2000*nm
x = np.arange(dx/2,xmax,dx) #range of distances

nop = 2 #number of particles
q = [eV,-eV] #charges for particles


V_graph = (-0.1*x/nm+0.6)*eV
V_graph[V_graph < 0] = 0

n = np.zeros_like(x)
E = np.zeros_like(x)

def find_Phi_rho_from_Phi0(x, V_graph, Phi0):
    rho = np.zeros_like(x)
    Phi = np.zeros_like(x)
    Phi[len(x)-2] = Phi0

    for j in range(len(x)-2,0,-1):
        rho[j] = 0
        for i in range(nop):
            #print n_0*np.exp(-beta*(V_graph[j] + q[i]*Phi[j]))
            n[j] = n_0*np.exp(-beta*(V_graph[j] + q[i]*Phi[j])) #Concentration
            rho[j] += q[i]*n[j] #Charge Density
        E[j] = -rho[j]/(e_w*e_0)*dx + E[j+1]
        Phi[j-1] = E[j]*dx + Phi[j]
    return Phi, rho

def find_V_from_Phi(x, Phi):
    return Phi[0]

def find_Q_from_rho(x, rho):
    return -rho.sum()*(x[1]-x[0])

# plt.figure()
# plt.plot(x/nm,V_graph/eV)
# plt.xlabel('$x$ (nm)')
# plt.ylabel('$V$ (eV)')

# plt.figure()
# plt.plot(x,find_Phi_rho_from_Phi0(x, V_graph, 0.001)[0])

Vmax = 2.0
dVgoal = 0.01

phi0max = 1e-4
dphi0 = 1e-5 # phi0max/20
phi0s = np.arange(0, phi0max + dphi0/2, dphi0)
Vs = np.zeros_like(phi0s)
Qs = np.zeros_like(phi0s)
for i in range(len(phi0s)):
    Phi,rho = find_Phi_rho_from_Phi0(x, V_graph, phi0s[i])
    Qs[i] = find_Q_from_rho(x, rho)
    Vs[i] = find_V_from_Phi(x, Phi)
    # print Vs[i], Qs[i]
    # plt.figure()
    # plt.plot(x/nm,Phi)
    # plt.xlabel('$x$ (nm)')
    # plt.ylabel('$\Phi$ (V)')
while Vs[-1] < Vmax:
    phi0 = phi0s[-1] + dphi0
    Phi,rho = find_Phi_rho_from_Phi0(x, V_graph, phi0)
    Q = find_Q_from_rho(x, rho)
    V = find_V_from_Phi(x, Phi)
    print("V = %g\tQ = %g" % (V, Q))
    if V - Vs[-1] > dVgoal:
        dphi0 /= 2
        print("Decreasing dphi0 to %g due to dV = %g" % (dphi0, V-Vs[-1]))
    else:
        Vs = np.append(Vs, [V])
        Qs = np.append(Qs, [Q])
        phi0s = np.append(phi0s, [phi0])

plt.figure()
plt.plot(Vs, Qs, 'k-')
plt.plot(-Vs, -Qs, 'k-')
plt.xlabel('$V$ (Volts)')
plt.ylabel('$Q$ (Coulombs/meter$^2$)')

plt.xlim(-Vmax, Vmax)

plt.savefig('Q-vs-V.pdf')

plt.show()



