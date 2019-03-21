from scipy import *
import matplotlib.pyplot as plt

import numpy as np
import numpy.linalg as la

def constructV(N, k, a, V_1, V_2 = 0.):
    dim = 2 * N + 1
    hbar = 6.582e-16 ## hbar in (eV s)
    m = 0.5109e6 ## mass of electron in (eV/c^2)
    c = 2.998e8 
    n = linspace(-N, N, 2 * N + 1) 
    k_n = k + (2 * pi * n / a)
    eps = (hbar ** 2 * (k_n * c) ** 2) / (2 * m)

    V = (V_1 * np.eye(dim, k = -1)) + (V_1 * np.eye(dim, k = 1)) + np.diag(eps) + (V_2* np.eye(dim, k = -2)) + (V_2 * np.eye(dim, k = 2))

    return V

def bandgap_plot(ax, N = 10, V_1 = 0.2, V_2 = 0.4):
    a = 5e-10
    kstep = 100
    k = linspace(-1, 1, kstep)
    energy = zeros((kstep, 3))

    for i in range(kstep):
        V = constructV(N, k[i] * pi / a, a, V_1, V_2)
        energy[i,:] = la.eigvalsh(V)[0:3]

    ax.plot(k, energy[:, 0])
    ax.plot(k, energy[:, 1])
    ax.plot(k, energy[:, 2])

fig = plt.figure(1)
fig.suptitle('Band gap diagrams')
ax1 = fig.add_subplot(1,3,1)
ax2 = fig.add_subplot(1,3,2)
ax3 = fig.add_subplot(1,3,3)

bandgap_plot(ax1, V_2 = 0.)
ax1.set_xlabel(r'$k$ / ($\frac{\pi}{a}$)') 
ax1.set_ylabel(r'Energy (eV)')
ax1.set_title(r'$V_1 = 0.2 eV$')

bandgap_plot(ax2)    
ax2.set_xlabel(r'$k$ / ($\frac{\pi}{a}$)') 
ax2.set_ylabel(r'Energy (eV)')
ax2.set_title(r'$V_1 = 0.2 eV, V_2 = 0.4 eV$')

bandgap_plot(ax3, V_1 = 4.0, V_2 = 0.)
ax3.set_xlabel(r'$k$ / ($\frac{\pi}{a}$)') 
ax3.set_ylabel(r'Energy (eV)')
ax3.set_title(r'$V_1 = 4.0 eV, V_2 = 0.0 eV$')

fig.tight_layout()
plt.show()