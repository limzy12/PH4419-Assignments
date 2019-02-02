from scipy import *
import matplotlib.pyplot as plt

### Define a function to calculate the potential due to a single point charge
def calcPotential(r, Q, X, Y):
    # r = position vector of the point charge
    # Q = charge magnitude of the point charge
    # X, Y = meshgrid coordinates
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

x = y = linspace(-2, 2, 100)
X, Y = meshgrid(x,y)
print(calcPotential([1., 1.], 1., X, Y))
    