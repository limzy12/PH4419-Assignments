from scipy import *
import matplotlib.pyplot as plt

import numpy.random as rndm
import time

def metropolisProb(dE, T):
    if dE < 0:
        return True
    elif rndm.rand() < exp(-dE / T):
        return True
    else: 
        return False

def computeTotalEnergy(state, J): ## This only works for systems 3 x 3 or greater
    (Nx, Ny) = state.shape
    totalEnergy = 0
    for i in range(Nx):
        for j in range(Ny):
            totalEnergy += (-J * state[i,j] * (state[(i + 1) % Nx, j] + state[i, (j + 1) % Ny]))
    return totalEnergy

def energyChange(state, J, i, j):
    (Nx, Ny) = state.shape
    ## By flipping the state, the change in energy is simply the initial energy of the particle with its nearest neighbours multiplied by -2
    energy = -J * state[i,j] * (state[(i + 1) % Nx, j] + state[i, (j + 1) % Ny] + state[(i - 1) % Nx, j] + state[i, (j - 1) % Ny])
    return -2 * energy
    

def ising_mc(J = 1., Nx = 16, Ny = 16, nsteps = 50000):
    Tmax = max(1/J, 10/J)
    Tmin = min(1/J, 10/J)
    numTsteps = 50
    threshold = 5000
    T = linspace(Tmin, Tmax, numTsteps)

    avgEnergyPoints = []
    avgMagPoints = []
    heatCapPoints = []
    magSusPoints = []
    algoStart = time.time()

    for temp in T:
        stepStart = time.time()
        ## Start with uniformly distributed state for (hopefully) quicker convergence
        state = rndm.choice([-1, 1], size = (Nx, Ny))
        totalEnergy = computeTotalEnergy(state, J)
        totalMag = mean(state)
        energyPoints = array([totalEnergy])
        magPoints = array([totalMag])

        for i in range(nsteps): 
            ## Choose a random site to flip
            flipX = rndm.randint(0, Nx)
            flipY = rndm.randint(0, Ny)

            ## Check the energy change and whether the flip is accepted
            dE = energyChange(state, J, flipX, flipY)
            if metropolisProb(dE, temp):
                totalMag += ((-2 * state[flipX,flipY]) / (Nx * Ny))
                totalEnergy += dE
                state[flipX, flipY] *= -1 
                
            ## Calculate energy and magnetisation
            energyPoints = append(energyPoints, array([totalEnergy]))
            magPoints = append(magPoints, array([totalMag]))

        ## Calculate plotting values from cutoff timestep (threshold)
        avgE = mean(energyPoints[threshold:])
        avgM = mean(magPoints[threshold:])
        heatCap = var(energyPoints[threshold:]) / temp
        magSus = var(magPoints[threshold:])

        avgEnergyPoints.append(avgE)
        avgMagPoints.append(avgM)
        heatCapPoints.append(heatCap)
        magSusPoints.append(magSus)
        
        stepEnd = time.time()
        print('Step temperature ' + str(temp) + ' runtime = ' + str(stepEnd - stepStart))

    algoEnd = time.time()
    print('Overall runtime = ' + str(algoEnd - algoStart) + '; Average time per step = ' + str((algoEnd - algoStart) / numTsteps))

    fig = plt.figure(1) 
    ax1 = fig.add_subplot(2, 2, 1)
    ax2 = fig.add_subplot(2, 2, 2)
    ax3 = fig.add_subplot(2, 2, 3)
    ax4 = fig.add_subplot(2, 2, 4)
    fig.suptitle('2D Ising model Monte-Carlo Simulation, ' + str(numTsteps) + ' temperature steps with ' + str(nsteps) + ' Monte-Carlo moves')
    ax1.scatter(T, avgEnergyPoints, c = 'r', marker = '.', label = r'$<E>$')
    ax2.scatter(T, avgMagPoints, c = 'b', marker = '.', label = r'$<M>$')
    ax3.scatter(T, heatCapPoints, c = 'r', marker = '.', label = r'$C_B$')
    ax4.scatter(T, magSusPoints, c = 'b', marker = '.', label = r'$\chi_m$')
    ax1.set_xlabel(r'Temperature, $T$ (arb. units)')
    ax2.set_xlabel(r'Temperature, $T$ (arb. units)')
    ax3.set_xlabel(r'Temperature, $T$ (arb. units)')
    ax4.set_xlabel(r'Temperature, $T$ (arb. units)')
    ax1.set_ylabel(r'Average Energy, $<E>$ (arb. units)')
    ax2.set_ylabel(r'Average Magnetisation, $<M>$ (arb. units)')
    ax3.set_ylabel(r'Heat Capacity, $C_B$ (arb. units)')
    ax4.set_ylabel(r'Magnetic Susceptibility, $\chi_m$ (arb. units)')
    ax1.set_title('Average energy against temperature')
    ax2.set_title('Average magnetisation against temperature')
    ax3.set_title('Heat capacity against temperature')
    ax4.set_title('Magnetic susceptibility against temperature')
    ax1.legend(loc = 'best')
    ax2.legend(loc = 'best')
    ax3.legend(loc = 'best')
    ax4.legend(loc = 'best')
    plt.show()

ising_mc()