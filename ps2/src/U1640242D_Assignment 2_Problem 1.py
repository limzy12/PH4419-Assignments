from scipy import *
import matplotlib.pyplot as plt

import numpy as np
import numpy.linalg as la

######################
#      Part (a)      # 
######################

## Define a function to construct the V matrix
def constructV(N, k, a, V_1, V_2 = 0.):
    ## Set up relevant parameters
    dim = 2 * N + 1
    hbar = 6.582e-16 ## hbar in (eV s)
    m = 0.5109e6 ## mass of electron in (eV/c^2)
    c = 2.998e8 
    n = linspace(-N, N, 2 * N + 1) 
    k_n = k + (2 * pi * n / a)
    eps = (hbar ** 2 * (k_n * c) ** 2) / (2 * m)

    ## Construct the matrix
    V = (V_1 * np.eye(dim, k = -1)) + (V_1 * np.eye(dim, k = 1)) + np.diag(eps) + (V_2* np.eye(dim, k = -2)) + (V_2 * np.eye(dim, k = 2))

    return V

## Define a function to plot the band gap diagram for given N, V_1 and V_2
def bandgap_plot(ax, N = 10, V_1 = 0.2, V_2 = 0.4):
    ## Define dome relevant parameters 
    a = 5e-10
    kstep = 100
    k = linspace(-1, 1, kstep)
    energy = zeros((kstep, 3))

    ## For each value of k, construct V and solve for 3 smallers eigenvalues
    for i in range(kstep):
        V = constructV(N, k[i] * pi / a, a, V_1, V_2)
        energy[i,:] = la.eigvalsh(V)[0:3]

    ## Plot the bandgap diagram
    ax.plot(k, energy[:, 0])
    ax.plot(k, energy[:, 1])
    ax.plot(k, energy[:, 2])

### Here we plot the band gap diagrams of the three scenarios given
fig = plt.figure(1)
fig.suptitle('Band gap diagrams')
ax1 = fig.add_subplot(1,3,1)
ax2 = fig.add_subplot(1,3,2)
ax3 = fig.add_subplot(1,3,3)

## V_1 = 0.2 eV, V_2 = 0.0 eV
bandgap_plot(ax1, V_2 = 0.)
ax1.set_xlabel(r'$k$ / ($\frac{\pi}{a}$)') 
ax1.set_ylabel(r'Energy (eV)')
ax1.set_title(r'$V_1 = 0.2 eV, V_2 = 0.0 eV$')
'''
Analytically, the bandgap has the expression 2V_1 = 0.4 eV in this case. From the plot, we can see that indeed the bandgap is approximately 0.4 eV
'''

## V_1 = 0.2 eV, V_2 = 0.4 eV
bandgap_plot(ax2)    
ax2.set_xlabel(r'$k$ / ($\frac{\pi}{a}$)') 
ax2.set_ylabel(r'Energy (eV)')
ax2.set_title(r'$V_1 = 0.2 eV, V_2 = 0.4 eV$')
'''
Here, we expect to see two bandgaps of 2V_1 and 2V_2 respectively. Indeed, this is observed as the first bandgap is approximately 0.4 eV (like before) and the second bandgap is approximately 0.8 eV.
'''

bandgap_plot(ax3, V_1 = 4.0, V_2 = 0.)
ax3.set_xlabel(r'$k$ / ($\frac{\pi}{a}$)') 
ax3.set_ylabel(r'Energy (eV)')
ax3.set_title(r'$V_1 = 4.0 eV, V_2 = 0.0 eV$')
'''
Now we see that the energy bands do not really exist, but rather they are discrete energy levels. The electrons can no longer travel freely though the lattice. Furthermore, we see that the lowest enery level is in fact negative. This corresponds to an electron bounded to the large potential. 
'''

fig.tight_layout()


######################
#      Part (b)      # 
######################

## Define a function to return the energy given k_x, k_y, a and gamma.
def E_k(kx, ky, gamma, a = 1.42e-10):
    E1 = 4 * cos(sqrt(3) * a * ky / 2)
    E2 = cos(3 * a * kx / 2) + cos(sqrt(3) * a * ky / 2)
    return gamma * sqrt(1 + E1 * E2) 

## Define a function to plot the required plots
def graphene_plot(a, gamma = 3.0):

    ## Setting up of plots
    fig = plt.figure(2)
    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)

    ## We shall normalise with respect to enegry at the Gamma point
    gammaEnergy = E_k(0., 0., gamma)

    ## Plotting E_k vs k_x with k_y = 0
    steps = 1000
    kx = linspace(0, 2 * pi / (3 * a), steps) ## Plot Gamma to M
    # Calculate and normalise E_k
    energy = E_k(kx, 0., gamma) / gammaEnergy
    ax1.plot(linspace(0, 1, steps), energy)
    ax1.set_title(r'Plot of $E(k)$ (normalised w.r.t $E(\Gamma)$) vs $k_x$, with $k_y = 0$')
    ax1.set_xlabel(r'k / ($\frac{2 \pi}{3a}$)')
    ax1.set_ylabel('Normalised Energy')
    ax1.set_ylim([0, 1])

    ## Plotting E_k vs k_y with k_x = 0
    ky = linspace(0, 2 * pi / (sqrt(3) * a), steps) ## Plot Gamma to K on ky
    # Calculate and normalise E_k
    energy = E_k(0., ky, gamma) / gammaEnergy
    ax2.plot(linspace(0, 1, steps), energy)
    ax2.set_title(r'Plot of $E(k)$ (normalised w.r.t $E(\Gamma)$) vs $k_y$, with $k_x = 0$')
    ax2.set_xlabel(r'k / ($\frac{2 \pi}{\sqrt{3}a}$)')
    ax2.set_ylabel('Normalised Energy')
    ax2.set_ylim([0, 1])
    fig.tight_layout()

    ### Here we plot the energy as we move between the high symmetry points
    ## Setting up plots
    fig = plt.figure(3)
    fig.suptitle(r'Band energy across $\Gamma \rightarrow M \rightarrow K \rightarrow \Gamma$ normalised to energy at $\Gamma$')
    ax1 = fig.add_subplot(1,3,1)
    ax2 = fig.add_subplot(1,3,2,sharey=ax1)
    ax3 = fig.add_subplot(1,3,3,sharey=ax1)

    ## Plot for Gamma to M
    energy = E_k(kx, 0., gamma) / gammaEnergy
    ax1.plot(linspace(0, 1, steps), energy)
    ax1.set_ylim([0, max(energy)])
    ax1.set_xlim([0, 1])
    ax1.set_ylabel('Normalised energy')
    ax1.set_xlabel(r'$\frac{2\pi}{3a}(k, 0)$')
    ax1.set_title(r'$\Gamma \rightarrow M$')

    ## Plot for M to K
    k_KM = linspace(0, pi / (sqrt(3) * a), steps)
    energy = E_k(kx[-1], k_KM, gamma) / gammaEnergy
    ax2.plot(linspace(0, 0.5, steps), energy)
    plt.setp(ax2.get_yticklabels(), visible = False)
    ax2.set_xlim([0, 0.5])
    ax2.set_xlabel(r'$\frac{2\pi}{\sqrt{3}a}(\frac{1}{\sqrt{3}}, k)$')
    ax2.set_title(r'$M \rightarrow K$')
    
    ## Plot for K to Gamma
    k_KGamma = linspace(pi /(sqrt(3) * a), 0, steps)
    energy = E_k(k_KGamma * 2 / sqrt(3), k_KGamma, gamma) / gammaEnergy
    ax3.plot(linspace(1, 0, steps), energy)
    ax3.set_xlim([0, 1])
    ax3.invert_xaxis()
    ax3.set_title(r'$K \rightarrow \Gamma$')
    ax3.set_xlabel(r'$\frac{\pi}{\sqrt{3}a}(\frac{2}{\sqrt{3}} k, k)$')
    plt.setp(ax3.get_yticklabels(), visible = False)

graphene_plot(1.42e-10)

plt.show()