from scipy import *
import matplotlib.pyplot as plt

import numpy as np

## We define a function to set the required boundary conditions
def set_boundaries(psi, xsteps, ysteps, Lx, Ly, xCutoff, yCutoff, u_flow):

    # Top boundary
    psi[-1, :] = u_flow * Ly

    # Bottom boundary
    psi[0, :] = 0

    # Cutout in corner
    psi[0:yCutoff, xCutoff:] = 0

    # Left and right boundaries (linear)
    for i in range(ysteps):
        psi[i, 0] = i * (u_flow * Ly) / (ysteps - 1)
        if i > yCutoff: 
            idx = i - yCutoff 
            psi[i, -1] = idx * (u_flow * Ly) / (ysteps - yCutoff - 1)
    return psi

## We define a function to solve the Laplace equation
def laplace_solver(epsilon = 1.e-3, nmax = 1.e4, u_flow = 5.):
    ## Setting up some other required constants
    xsteps = 150
    ysteps = 150
    Lx = 1.
    Ly = 1.
    xCutoff = int(0.6 * xsteps) - 1
    yCutoff = int(0.4 * ysteps) + 1
    error = 1
    itercount = 0

    ## Initialise the matrix
    psi = zeros((ysteps, xsteps)) ## row index -> y; column index -> x; 
    next_psi = zeros((ysteps, xsteps)) 
    # Set the boundary conditions
    set_boundaries(psi, xsteps, ysteps, Lx, Ly, xCutoff, yCutoff, u_flow)

    while (error > epsilon) and (itercount < nmax):
        # Compute the values at the next time step
        next_psi[1:-1, 1:-1] = 0.25 * (psi[1:-1, 2:] + psi[1:-1, :-2] + psi[2:,1:-1] + psi[:-2, 1:-1])
        # Calculate the error
        error = abs(mean(next_psi[1:-1, 1:-1] - psi[1:-1, 1:-1]))   

        # Count the iterations
        itercount += 1
        psi = copy(next_psi)
        # Reset boundary conditions
        set_boundaries(psi, xsteps, ysteps, Lx, Ly, xCutoff, yCutoff, u_flow)
        
    print("Number of iterations: " + str(itercount))
    print("Error: " + str(error))

    # Set up plot variables
    X = linspace(0, Lx, xsteps)
    Y = linspace(0, Ly, ysteps)
    streamX, streamY = np.gradient(psi)

    # Plotting of streamlines
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

'''
In order to make less iterations, we could increase number of steps. i.e. smaller but more finite elements. However, this should come with a corresponding decrease in epsilon to ensure that enough iterations have been run to properly evaluate all points in the domain. It should be noted that although we take less iterations, each iteration may take a longer time to run. 
'''