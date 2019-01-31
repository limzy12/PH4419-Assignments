from scipy import *
import matplotlib.pyplot as plt

from scipy import integrate # We require the numerical integration library from scipy
import math

### We define a function to convert Cartesian coordinates to Cylindrical coordinates ###
def cart2cyl(r): # The input argument is a 3D vector r in Cartesian coordinates
    x = r[0]
    y = r[1]
    z = r[2]
    rho = sqrt(x**2 + y**2)
    if x == 0 and y == 0:
        phi = 0.
    elif x == 0 and y > 0:
        phi = math.pi / 2
    elif x == 0 and y < 0:
        phi = - math.pi / 2
    elif x > 0 and y == 0:
        phi = 0.
    elif x < 0 and y == 0:
        phi = math.pi
    else:
        phi = arctan(y / x)

    return [rho, phi, z]

### We define a function to convert Cylindrical coordinates to Cartesian coordinates ###
def cyl2cart(r): # The input argument is a 3D vector r in Cylindrical coordinates
    rho = r[0]
    phi = r[1]
    z = r[2]
    x = rho * cos(phi)
    y = rho * sin(phi)

    return [x, y, z]

### We define the integrals to be evaluated
def integrandRho(phi, rho, z, R): # this is the phi component of the integral
    denominator = ((R ** 2) + (rho ** 2) + (z ** 2) - (2 * rho * R * sin(phi))) ** (3/2)
    numerator = z * sin(phi)

    return numerator / denominator

def integrandZ(phi, rho, z, R): # this is the z component of the integral

    denominator = ((R ** 2) + (rho ** 2) + (z ** 2) - (2 * rho * R * sin(phi))) ** (3/2)
    numerator = R - rho * sin(phi)

    return numerator / denominator

def Coil_Bfield(r, current, radius):
    [rho, phi, z] = cart2cyl(r) # convert the given coordinates into cylindrical coordinates

    # Setting the relevant constants
    I = current
    R = radius

    # Defining the integration limits 
    uppLimit = 2 * math.pi
    lowLimit = 0.

    # Calculating coefficient
    coeff = I * R / (4 * math.pi)


    # Computing the integrals component-wise
    integralRho = integrate.quad(integrandRho, lowLimit, uppLimit, (rho, z, R))
    integralZ = integrate.quad(integrandZ, lowLimit, uppLimit, (rho, z ,R))

    BRho = coeff * integralRho[0]
    BZ = coeff * integralZ[0]

    # Converting the result into Cartesian co-ordinates
    resultVec = [BRho, 0., BZ]
    resultCart = cyl2cart(resultVec)

    return resultCart

def plotBZ():
    # Defining the parameters of the plot
    zMin = -5.
    zMax = 5.
    numPts = 100
    I = 1.
    R = 1.
    x = 0.
    y = 0.

    # Setting up the plot values
    z = linspace(zMin, zMax, numPts)
    BZ = zeros(len(z))

    for i in range(len(z)):
        r = [x, y, z[i]]
        B = Coil_Bfield(r, I, R)
        BZ[i] = B[2]
    
    plt.figure(1)
    plt.plot(z, BZ)
    title = "Plot of $B_z$ vs $z$ at x = " + str(x) + " and y = " + str(y) 
    plt.title(title)
    plt.xlabel("$z$")
    plt.ylabel("$B_z$")
    plt.show()

plotBZ()