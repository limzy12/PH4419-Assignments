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
    phi = arctan2(y, x)

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
    resultVec = [BRho, 0., BZ]      ### TODO: angular coordinate?
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

    # Calculation of B_z values
    for i in range(len(z)):
        r = [x, y, z[i]]
        B = Coil_Bfield(r, I, R)
        BZ[i] = B[2]

    # Plotting the graph
    plt.figure(1)
    plt.plot(z, BZ)

    # Plot titles and labels
    title = "Plot of $B_z$ vs $z$ at x = " + str(x) + " and y = " + str(y) 
    plt.title(title)
    plt.xlabel("$z$")
    plt.ylabel("$B_z$")
    plt.show()

def Helmholtz_Bfield(z, current, radius, distance):
    # to vary rho, it suffices to vary x since the result has no explicit phi dependence.
    # Defining the parameters of the plot
    xMin = -2. * radius
    xMax =  2. * radius
    numPts = 100

    # Setting up the plot values
    x = linspace(xMin, xMax, numPts)
    BZ = zeros(len(x))

    # Calculation of B_z values
    # We calculate the B_z contribution by each loop individually, making appropriate coordinate transforms to center the loop at the origin.
    for i in range(len(x)):
        r1 = [x[i], 0., z + distance / 2] 
        r2 = [x[i], 0., z - distance / 2]
        B = Coil_Bfield(r1, current, radius) + Coil_Bfield(r2, current, radius)
        BZ[i] = B[2]

    # Setting up the plot
    plt.figure(2)
    plt.plot(x, BZ)

    # Plot titles and labels
    title = r"Plot of $B_z$ vs $\rho$ at $R$ = " + str(radius) + " and $L$ = " + str(distance)
    plt.title(title)
    plt.xlabel(r"$\rho$")
    plt.ylabel("$B_z$")
    plt.show()

Coil_Bfield([1., 1., 1.], 1., .5)

plotBZ()
Helmholtz_Bfield(0., 1., 1., 1.,)

"""
Discussion on the behaviour of the magnetic field.

Making the relevant plots, we see that for z = 0, when R = L, there is a region where B_z is constant. This region is bounded by the radius of the loop itself. i.e. rho < R. 

When R > L, we observe an 'M'-shaped curve, where there are two peaks. The peaks occur where rho = R and in between the peaks there is a local minimum at rho = 0. On the other hand when R < L, There is only one global maximum of B_z occurring at rho = 0. 

When z != 0, in general, the maximal value of B_z decreases. In the R > L case, the difference in the maxima and the local minima is decreased.
"""