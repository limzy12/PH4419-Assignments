from scipy import *
import matplotlib.pyplot as plt

import numpy.random as rndm
from mpl_toolkits.mplot3d import Axes3D

def torus_estimate(r, R, N):
    total_pts = 0
    count = 0
    points = zeros([N, 3])
    while count < N:
        x = -1.5 + 3. * rndm.rand() 
        y = -1.5 + 3. * rndm.rand()
        z = -1 + 2. * rndm.rand()

        if (sqrt(x ** 2 + y ** 2) - R) ** 2 + z ** 2 <= r ** 2 : 
            points[count] = [x, y, z]
            count += 1

        total_pts += 1
    
    est_vol = count / total_pts * 18

    return est_vol, points

est_vol, points = torus_estimate(0.5, 1., 1000)
plt.figure(1)
ax = plt.gca(projection = "3d")
ax.scatter(points[:,0], points[:,1], points[:,2])
ax.set_zlim3d(-1.5, 1.5)
plt.show()

