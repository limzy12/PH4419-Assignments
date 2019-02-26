### Classical Oscillator ###

from scipy import *
import matplotlib.pyplot as plt

import scipy.integrate as integrate

### We define a function that computes the time derivatives at t
def oscillator(u, t, k, D, m):
    x = u[0]
    v = u[1]
    udot = [v, -(D * v + k * x)/m]
    return udot

### We define a function to plot the required plots
def oscillator_plot(t, m, k, D, u0):
    ### Solve for the trajectory
    sol = integrate.odeint(oscillator, u0, t, (k, D, m))
    trajectoryX = sol[:, 0]
    trajectoryV = sol[:, 1]

    ### Setup for phase portrait plot
    xmax = 1.5 * max(trajectoryX)
    vmax = 1.5 * max(trajectoryV)
    numPts = 30
    x = linspace(-xmax, xmax, numPts)
    v = linspace(-vmax, vmax, numPts)
    X, V = meshgrid(x,v)

    dX = zeros(shape(X))
    dV = zeros(shape(X))

    for i in range(numPts):
        for j in range(numPts):
            u = [X[i,j], V[i,j]]
            udot = oscillator(u, 0., k, D, m)
            dX[i,j] = udot[0]
            dV[i,j] = udot[1]
    
    ## Normalise the length of the arrows 
    norm = hypot(dX, dV)
    dX = dX / norm
    dV = dV / norm

    ## Plotting phase portrait and trajectory
    plt.figure(1)
    plt.plot(trajectoryX, trajectoryV)
    plt.quiver(X, V, dX, dV)
    plt.xlabel(r'Position, $x(t)$')
    plt.ylabel(r'Velocity, $v(t)$')
    plt.title(r'Phase portrait and trajectory of $m\ddot{x} + D\dot{x} + kx = 0$;' + '\n m = ' + str(m) + ', D = ' + str(D) + ', k = '+ str(k) + r', $u_0$ = ' + str(u0))
    plt.show()

u0 = [0., pi/2]
t = linspace(0., 10., 1000)
oscillator_plot(t, 1., 1., 0., u0)

