from scipy import *
import matplotlib.pyplot as plt

import numpy as np

def set_boundaries(psi, xsteps, ysteps, Lx, Ly, xCutoff, yCutoff, u_flow):
    psi[-1, :] = u_flow * Ly
    psi[0, :] = 0
    psi[0:yCutoff, xCutoff:] = 0
    for i in range(ysteps):
        psi[i, 0] = i * (u_flow * Ly) / (ysteps - 1)
        if i > yCutoff: 
            idx = i - yCutoff 
            psi[i, -1] = idx * (u_flow * Ly) / (ysteps - yCutoff - 1)
    return psi

def laplace_solver(epsilon = 1.e-3, nmax = 1.e4, u_flow = 5.):
    ## Setting up some other required constants
    xsteps = 100
    ysteps = 100
    Lx = 1.
    Ly = 1.
    xCutoff = int(0.6 * xsteps) - 1
    yCutoff = int(0.4 * ysteps) + 1
    error = 1
    itercount = 0

    ## Initialise the matrix
    psi = zeros((ysteps, xsteps)) ## row index -> y; column index -> x; 
    next_psi = zeros((ysteps, xsteps))
    set_boundaries(psi, xsteps, ysteps, Lx, Ly, xCutoff, yCutoff, u_flow)

    while (error > epsilon) and (itercount < nmax):
        next_psi[1:-1, 1:-1] = 0.25 * (psi[1:-1, 2:] + psi[1:-1, :-2] + psi[2:,1:-1] + psi[:-2, 1:-1])
        error = abs(mean(next_psi[1:-1, 1:-1] - psi[1:-1, 1:-1]))

        itercount += 1
        psi = copy(next_psi)
        set_boundaries(psi, xsteps, ysteps, Lx, Ly, xCutoff, yCutoff, u_flow)
        
    print("Number of iterations: " + str(itercount))
    print("Error: " + str(error))

    X = linspace(0, Lx, xsteps)
    Y = linspace(0, Ly, ysteps)
    streamX, streamY = np.gradient(psi)

    plt.figure(1)
    plt.streamplot(X, Y, streamX, -streamY, 2, color = np.hypot(streamX, streamY), cmap = 'jet')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    plt.title('Incompressible fluid flow in a channel') 
    plt.colorbar().set_label('Velocity of fluid')
    plt.xlim(0, Lx)
    plt.ylim(0, Ly)
    plt.vlines(0.6 * Lx, 0, 0.4 * Ly)
    plt.hlines(0.4 * Ly, 0.6 * Lx, Lx)
    plt.show()

laplace_solver()