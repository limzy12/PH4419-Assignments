### Classical Oscillator ###

from scipy import *
import matplotlib.pyplot as plt

import scipy.integrate as integrate

##############
## Part (a) ##
##############

### We define a function that computes the time derivatives at t
def oscillator(u, t, k, D, m): # k ,D, m are as defined in the differential equation
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
    plt.legend(['Trajectory'])
    plt.show()

u0 = [0., pi/2]
t = linspace(0., 10., 1000)
oscillator_plot(t, 1., 1., 0., u0)

##############
## Part (b) ##
##############

### Define the differential equation describing the dynamics of the double well
def doubleWell(u, t, m):
    x = u[0]
    v = u[1]
    return [v, (x - x ** 3) / m]

### Function to plot the required plots
def double_potential_plot(m = 1): 

    ### Setting up variables for phase portrait plot
    xmax = 2
    vmax = 2
    numPts = 40
    x = linspace(-xmax, xmax, numPts)
    v = linspace(-vmax, vmax, numPts)
    X, V = meshgrid(x, v)

    dX = zeros(shape(X))
    dV = zeros(shape(X))

    for i in range(numPts):
        for j in range(numPts):
            u = [X[i,j], V[i,j]]
            udot = doubleWell(u, 0., m)
            dX[i,j] = udot[0]
            dV[i,j] = udot[1]

    ## Normalise the length of the arrows 
    norm = hypot(dX, dV)
    dX = dX / norm
    dV = dV / norm

    ### Calculate three closed trajectories
    ## Setting up the integration interval
    tmax = 10.
    numTPts = 1000
    t = linspace(0, tmax, numTPts)

    ## Defining initial points for the three trajectories
    u1 = [0, 0.25 * xmax]
    u2 = [0, 0.5 * xmax]
    u3 = [0, 0.75 * xmax]
    traj1 = integrate.odeint(doubleWell, u1, t, (m,))
    traj2 = integrate.odeint(doubleWell, u2, t, (m,))
    traj3 = integrate.odeint(doubleWell, u3, t, (m,))

    ## Plotting phase portrait and trajectory
    plt.figure(2)
    plt.plot(traj1[:,0], traj1[:,1])
    plt.plot(traj2[:,0], traj2[:,1])
    plt.plot(traj3[:,0], traj3[:,1])
    plt.quiver(X, V, dX, dV)
    plt.xlabel(r'Position, $x(t)$')
    plt.ylabel(r'Velocity, $v(t)$')
    plt.ylim([-1.25 * vmax, 1.25 * vmax])
    plt.title(r'Phase portrait of double potential well with $m = $' + str(m) + '\n and three closed trajectories')
    plt.legend(['Trajectory of initial point ' + str(u1), 'Trajectory of initial point ' + str(u2), 'Trajectory of initial point ' + str(u3)], loc = 1)
    plt.show()

double_potential_plot()

### How do the different initial conditions affect the harmonic motion of the particle?
'''
With different initial conditions, the particle starts out with different amounts of energy. This energy is conserved throughout the harmonic motion since there are no dissipative forces. The higher the initial energy of the particle, the less it is affected by the step in the potential. This is observed in the phase portrait as a more regular, oblong trajectory. Conversely, particles with less initial energy are more sensitive to the step in the potential. This results in the figure 8 trajectory as seen in the phase portrait. 
'''