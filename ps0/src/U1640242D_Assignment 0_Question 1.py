from scipy import *
import matplotlib.pyplot as plt

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

### Define a function to calculate the field components due to a single point charge
def calcField(r, Q, X, Y): 
    # r = position vector of the point charge
    # Q = charge magnitude of the point charge
    # X, Y = meshgrid coordinates
    rx = r[0]
    ry = r[1]

    # Sanity checks to prevent errors
    assert X.ndim == Y.ndim == 2
    assert shape(X) == shape(Y)

    # Calculating required coefficients
    distSq = hypot((X - rx), (Y - ry)) ** 2
    coeff = 1 / (4 * pi) # Since we take permittivity = 1

    # Setting up result matrices
    E = zeros(shape(X)) # Magnitude of E-field
    EX = EY = zeros(shape(X)) # x- , y-components of E-field

    # Calculate E-field
    E = coeff * Q / distSq
    
    # Calculate components
    EX = E * cos(arctan2(Y - ry, X - rx))
    EY = E * sin(arctan2(Y - ry, X - rx))

    # Normalise the arrows
    scale = hypot(EX, EY)
    EX = 3 * EX / scale
    EY = 3 * EY / scale

    return [EX, EY]

### Define a function to draw the required plots for part a)
def plotEquipotentialAndField():
    r = [[-1., -1.], [1., 1.]]
    q = 1

    xMin = yMin = -2.
    xMax = yMax = 2.
    numPtsXpotential = 1000
    numPtsYpotential = 1000
    numPtsXfield = 40
    numPtsYfield = 40

    x1 = linspace(xMin, xMax, numPtsXpotential)
    y1 = linspace(yMin, yMax, numPtsYpotential)
    x2 = linspace(xMin, xMax, numPtsXfield)
    y2 = linspace(yMin, yMax, numPtsYfield)

    X1, Y1 = meshgrid(x1, y1)
    X2, Y2 = meshgrid(x2, y2)

    V = zeros(shape(X1))
    EX = zeros(shape(X2))
    EY = zeros(shape(X2))
    Efield = zeros(shape(X2)) 

    for i in range(len(r)):
        V = V + calcPotential(r[i], q, X1, Y1)
        Efield = calcField(r[i], q, X2, Y2)
        EX += Efield[0]
        EY += Efield[1]
    
    plt.figure(1)
    plt.contour(X1, Y1, V)
    plt.quiver(X2, Y2, EX, EY)
    plt.show()

plotEquipotentialAndField()
    