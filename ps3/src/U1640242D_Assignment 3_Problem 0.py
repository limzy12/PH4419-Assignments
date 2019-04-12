from scipy import *
import matplotlib.pyplot as plt

import numpy.random as rndm

## Function to return the transformed point given the chosen transformation.
def transform(i, v):
    if i == 1:
        A = [[0., 0.], [0., 0.16]]
        b = [0., 0.]
    elif i == 2:
        A = [[0.85, 0.04], [-0.04, 0.85]]
        b = [0., 1.6]
    elif i == 3:
        A = [[0.2, -0.26], [0.23, 0.22]]
        b = [0., 1.6]
    elif i == 4:
        A = [[-0.15, 0.28], [0.26, 0.24]]
        b = [0., 0.44]
    return matmul(A, v) + b

## Function to choose the transformations and plot the results
def fractal(xmin = -2.2, xmax = 2.7, ymin = 0, ymax = 10, N = 100000):
    points = zeros((N + 1, 2))

    for idx in range(N):
        i = rndm.choice([1, 2, 3, 4], p = [0.01, 0.85, 0.07, 0.07]) ## The choice function selects values from the given array with the corresponding probabilty denoted in p
        points[idx + 1] = transform(i, points[idx])

    # Plotting
    plt.figure(1)
    plt.scatter(points[:, 0], points[:, 1], s = 1, c = 'g', marker = '.')
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.title(r'Generated Barnsley fern using N = ' + str(N))
    plt.legend(['Points generated'], loc = 'best')
    plt.show()

## Plot the fractal for N points
fractal(N = 100000)