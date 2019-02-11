from scipy import *
import matplotlib.pyplot as plt

import numpy.random as rndm
from mpl_toolkits.mplot3d import Axes3D

### Define a function to estimate the volume of a torus via Monte-Carlo integration
def torus_estimate(r, R, N):
    ## Initialising some variables
    total_pts = 0
    count = 0
    points = zeros([N, 3])
    
    ## Monte Carlo integration
    while count < N:
        ## Generation of random point
        x = -1.5 + 3. * rndm.rand() 
        y = -1.5 + 3. * rndm.rand()
        z = -1 + 2. * rndm.rand()

        ## Check if point is in the torus
        if (sqrt(x ** 2 + y ** 2) - R) ** 2 + z ** 2 <= r ** 2 : 
            points[count] = [x, y, z]
            count += 1

        total_pts += 1
    
    ## Compute the estimated volume
    est_vol = count / total_pts * 18

    return est_vol, points


### Define a function to estimate the error and make a log-log plot against N (Bonus Question)
def error_plot(r, R):
    ## Compute the actual volume
    actualVol = pi * (r ** 2) * 2 * pi * R

    ## Setting up plot variables
    minLimit = 1
    maxLimit = 6
    numPts = 10
    N = logspace(minLimit, maxLimit, numPts)
    error = zeros(len(N))

    ## Estimate the volume and area for each N
    for i in range(len(N)):
        estVol = torus_estimate(r, R, int(N[i]))[0]
        error[i] = abs(estVol - actualVol)

    ## Setting up the plot
    plt.figure(1)
    plt.loglog(N, error)

    ## Plot titles and labels
    plt.title(r"Log-log plot of estimation error against $N$")
    plt.xlabel(r"$N$")
    plt.ylabel("Estimation error")
    plt.show()

## Estimate the volume of a torus and obtain the points within the torus
est_vol, points = torus_estimate(0.5, 1., 1000)

## Plot the points that are within the torus
plt.figure(1)
ax = plt.gca(projection = "3d")
ax.scatter(points[:,0], points[:,1], points[:,2])

## We scale the display to be able to see the torus clearly
ax.set_zlim3d(-1.5, 1.5)

## Plot titles and labels
plt.title("Scatter plot of points generated in the torus")
plt.show()

## Function call to plot log-log plot of error vs N (Bonus question)
error_plot(0.5, 1.0)