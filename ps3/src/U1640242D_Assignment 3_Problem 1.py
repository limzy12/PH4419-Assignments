from scipy import *
import matplotlib.pyplot as plt

import numpy.random as rndm

dt = 1

## Function to calculate the decay probability given the half life
def decayProb(halfLife):
    return 1 - (2 ** (-dt / halfLife))

## Define some constants
N = 10000
PbHalfLife = 3.3 * 60
TlHalfLife = 2.2 * 60
BiHalfLife = 46. * 60

stableBiCount = []
activeBiCount = []
TlCount = []
PbCount = []

## Define the starting amounts of each particle
Bi213 = 10000
Tl209 = 0
Pb209 = 0
Bi209 = 0

activeBiCount.append(Bi213)
TlCount.append(Tl209)
PbCount.append(Pb209)
stableBiCount.append(Bi209)

for i in range(N):
    ## Decay from Pb-209 to Bi-213
    PbDecay = sum(rndm.rand(Pb209) < decayProb(PbHalfLife))
    Pb209 -= PbDecay 
    Bi209 += PbDecay

    ## Decay from Tl-209 to Pb-209
    TlDecay = sum(rndm.rand(Tl209) < decayProb(TlHalfLife))
    Tl209 -= TlDecay 
    Pb209 += TlDecay
     
    ## Decay from Bi-213
    BiDecay = sum(rndm.rand(Bi213) < decayProb(BiHalfLife))
    decayToPb = sum(rndm.rand(BiDecay) < 0.9791)
    decayToTl = BiDecay - decayToPb
    Bi213 -= BiDecay
    Pb209 += decayToPb
    Tl209 += decayToTl

    ## Append the current number of particles into an array
    activeBiCount.append(Bi213)
    TlCount.append(Tl209)
    PbCount.append(Pb209)
    stableBiCount.append(Bi209)

## Plotting
time = arange(0, N + 1)
plt.figure(1)
plt.plot(time, activeBiCount, label = 'Bi-213')
plt.plot(time, TlCount, label = 'Tl-209')
plt.plot(time, PbCount, label = 'Pb-209')
plt.plot(time, stableBiCount, label = 'Bi-209')
plt.ylabel('Number of particles')
plt.xlabel('Time')
plt.title('Decay of Bismuth-213')
plt.legend(loc = 'best')
plt.show()

'''
The curves of the two species of Bi are as expected. The curve of Pb is small due to its relatively short half-life, leading to quicker decay of the Pb atoms compared to its accumulation from the decay of Bi-213. It is similar for the Tl, but due to the low probability of decay via the Tl route, the curve looks to be relatively flat at around zero.

We have to implement the decay from bottom-up in order to avoid having cases of a 'double decay' in a particular time step. If we proceeded from the top-down, using the approach above, it is fully possible that a Bi-213 might decay into say, a Pb-209, and that Pb-209 then again decays into Bi-209. As such, it would seem as if the particle decayed directly from Bi-213 to Bi-209.

However, I would argue that it is possible to run a correct simulation from top-down, if one assigns the resultant number of particles in each decay into another variable, which is then to be considered for the next time step. i.e. one has to be careful to not allow scenarios of 'double decay'

'''
