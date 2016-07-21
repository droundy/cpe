from __future__ import division

import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt

concentrations = [1,10,100] # in mM

B = 0.5 # Tesla
e = 1.60217662e-19 # coulombs
I = 5e-6 # Amps
hbar = 1.0545718e-34 # J s
v_F = 1e6 # m/s

colors = {
    1: 'b',
    10: 'g',
    100: 'r',
}

for c in concentrations:
    data = np.loadtxt('hall-measurements/%d_mM.csv' % c)
    Vg = data[:,0]
    V_xx = data[:,1]
    V_Hall = data[:,2]

    plt.figure('Hall voltage versus gate voltage')
    plt.plot(Vg, V_Hall, label="%d mM" % c)
    plt.xlabel('$V_g$')
    plt.ylabel('$V_H$')
    plt.legend(loc='best')
    plt.savefig('hall-measurements/vh-vs-vg.pdf')

    plt.figure('xx voltage versus gate voltage')
    plt.plot(Vg, V_xx, label="%d mM" % c)
    plt.xlabel('$V_g$')
    plt.ylabel('$V_{xx}$')
    plt.legend(loc='best')
    plt.savefig('hall-measurements/vxx-vs-vg.pdf')

    ns = B*I/V_Hall/e

    plt.figure('simple carrier density versus gate voltage')
    plt.plot(Vg, ns, label="%d mM" % c)
    plt.xlabel('$V_g$')
    plt.ylabel('$n_s$')
    plt.legend(loc='best')
    plt.savefig('hall-measurements/ns-vs-vg.pdf')

    mu = V_Hall/V_xx/B

    plt.figure('mobility versus gate voltage')
    plt.plot(Vg, mu, label="%d mM" % c)
    plt.xlabel('$V_g$')
    plt.ylabel(r'$\mu$ (m$^2$/Vs)')
    plt.legend(loc='best')
    plt.savefig('hall-measurements/mu-vs-vg.pdf')

    Vqc = hbar*v_F*np.sqrt(np.pi*abs(ns))/e*np.sign(ns)

    plt.figure('quantum capacitance voltage versus gate voltage')
    plt.plot(Vg, Vqc, label="%d mM" % c)
    plt.xlabel('$V_g$')
    plt.ylabel('$V_{qc}$')
    plt.legend(loc='best')
    plt.savefig('hall-measurements/vqc-vs-vg.pdf')

    Vdl = Vg - Vqc

    plt.figure('simple carrier density versus double-layer voltage')
    plt.plot(Vdl, ns, '.', label="%d mM" % c)
    plt.xlabel('$V_{dl}$')
    plt.ylabel('$n_s$')
    plt.legend(loc='best')
    plt.savefig('hall-measurements/ns-vs-vdl.pdf')

    # mu = 0.6 # read by eye off of the simple analysis plot

    # xi = V_Hall/(mu*V_xx*B)
    for i in range(len(V_Hall)-1):
        if V_Hall[i+1]*V_Hall[i] < 0:
            VHi = abs(V_Hall[i])
            VHip = abs(V_Hall[i+1])
            V_dirac = (Vg[i]*VHip + Vg[i+1]*VHi)/(VHi + VHip)
            print i, V_Hall[i], V_Hall[i+1], Vg[i], Vg[i+1], V_dirac
    xi = np.tanh((Vg - V_dirac)/.13) # try a tanh

    plt.figure('carrier homogeneity versus gate voltage')
    plt.plot(Vg-V_dirac, xi, '-', label="%d mM" % c)
    plt.xlabel('$V_{g} - V_D$')
    plt.ylabel(r'$\xi$')
    plt.legend(loc='best')
    plt.savefig('hall-measurements/xi-vs-vg.pdf')

    mu = V_Hall/(xi*V_xx*B)

    plt.figure('new mobility versus gate voltage')
    plt.plot(Vg-V_dirac, mu, '-', label="%d mM" % c)
    plt.xlabel('$V_{g}-V_D$')
    plt.ylabel(r'$\mu$')
    plt.legend(loc='best')
    plt.savefig('hall-measurements/new-mu-vs-vg.pdf')

    nt = xi*I*B/V_Hall/e
    n_plus = (xi + 1)*nt/2
    n_minus = nt - n_plus

    plt.figure('total carrier density versus double-layer voltage')
    plt.plot(Vg-V_dirac, nt, colors[c]+'-', label="%d mM" % c)
    plt.plot(Vg-V_dirac, n_plus, colors[c]+'--')
    plt.plot(Vg-V_dirac, n_minus, colors[c]+':')
    plt.xlabel('$V_{g}-V_D$')
    plt.ylabel(r'$n_t$')
    plt.ylim(0, 0.2e17)
    plt.legend(loc='best')
    plt.savefig('hall-measurements/nt-vs-vg.pdf')

plt.show()
