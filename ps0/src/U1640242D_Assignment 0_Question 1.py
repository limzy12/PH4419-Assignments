from scipy import *
import matplotlib.pyplot as plt

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

### Define a function to calculate the potential due to a single point charge
def calcPotential(r, Q, X, Y):
    # r = position vector of the point charge
    # Q = charge magnitude of the point charge
    # X, Y = meshgrid coordinates

    # Sanity checks to prevent errors
    assert X.ndim == Y.ndim == 2
    assert shape(X) == shape(Y)

    rx = r[0]
    ry = r[1]

    # Calculating required coefficients
    dist = sqrt((X - rx) ** 2 + (Y - ry) ** 2)
    coeff = 1 / (4 * pi) # Since we take permittivity = 1

    # Setting up the result matrix
    V = zeros(shape(X))

    # Calculating the potential
    V = coeff * Q / dist

    return V

### Define a function to draw the required plots for part a)
def plotEquipotentialAndField():
    ## Define parameters for the charges
    r = [[-1., -1.], [1., 1.]]  # Positions of the charges
    q = 1                       # Strength of the charges

    ## Define plot parameters
    xMin = yMin = -2.
    xMax = yMax = 2.
    numPtsX = 50
    numPtsY = 50

    ## Seting up the meshgrid for plotting
    x = linspace(xMin, xMax, numPtsX)
    y = linspace(yMin, yMax, numPtsY)
    X, Y = meshgrid(x, y)

    ## Define some variables to store results
    V = zeros(shape(X))
    EX = zeros(shape(X))
    EY = zeros(shape(X))

    ## Calculating the potential
    for i in range(len(r)):
        V +=  calcPotential(r[i], q, X, Y)

    ## Calculating the Electric field
    Efield = np.gradient(V)
    norm = hypot(Efield[0], Efield[1])
    EY = -1 * Efield[0] / norm
    EX = -1 * Efield[1] / norm

    ## Setting up the plot
    plt.figure(1)
    plt.contour(X, Y, V, 100)
    plt.quiver(X, Y, EX, EY)
    plt.title("Quiver plot of the electric field and contour plot of equipotential lines")
    plt.show()

### Define a function to draw the desired plot for part b)
def plotPotential():
    ## Define parameters for the charges
    r = [[-1., -1.], [1., 1.]]  # Positions of the charges
    q = 1                       # Strength of the charges
    
    ## Define plot parameters
    xMin = yMin = -2.
    xMax = yMax = 2.
    numPtsXpotential = 50
    numPtsYpotential = 50

    ## Define some variables to store results
    x = linspace(xMin, xMax, numPtsXpotential)
    y = linspace(yMin, yMax, numPtsYpotential)
    X, Y = meshgrid(x, y)

    ## Define a variable to store results
    V = zeros(shape(X)) 

    ## Calculating the potential
    for i in range(len(r)):
        V += calcPotential(r[i], q, X, Y)

    ## Setting up the plot
    plt.figure(2)
    ax = plt.gca(projection = "3d")
    ax.plot_surface(X, Y, V, cmap = cm.coolwarm, linewidth = 1)

    ## Plot title
    plt.title("Surface plot of potential")
    plt.show()
    print(V)

plotEquipotentialAndField()
plotPotential()