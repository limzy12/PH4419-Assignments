### 1D Schrodinger Equation ###

from scipy import *
import matplotlib.pyplot as plt

import scipy.integrate as integrate

##############
## Part (a) ##
##############

### Define the system of first order DEs
def tise(u, x, nu0, epsilon):
    psi = u[0]
    psidx = u[1]
    if abs(x) < 0.5: ## Piecewise definition due to the nature of the potential 
        return [psidx, -2 * epsilon * psi]
    else:
        return [psidx, 2 * (nu0 - epsilon) * psi]

## Define a function to plot the probability density vs x for energy eigenvals
def FPW_demo(psi0, nu0, x, epsilon):
    for i in range(len(epsilon)):
        u0 = [psi0[0], 0.001] # Set the leftmost boundary condition and some non-trivial slope
        dx = x[1] - x[0]
        psiSolution = integrate.odeint(tise, u0, x, (nu0, epsilon[i])) # "Shoot" for a solution 
        psiSquared = abs(psiSolution[:,0]) ** 2 
        normConst = sum(psiSquared * dx) # Constant to normalise probability density
    
        axes = plt.figure(1).add_subplot(2,2,i+1)
        plt.tight_layout()
        axes.plot(x, psiSquared / normConst)
        axes.set_xlabel(r'$x$')
        axes.set_ylabel(r'$|\psi(x)|^2$')
        axes.set_title(r'Probability density for $\epsilon = $' + str(epsilon[i]))

    plt.show()

x = linspace(-1.2, 1.2, 1000)
epsilonGuess = [3.414, 13.476, 29.452, 48.144] # Values of first four eigenvalues obtained after trial and error
FPW_demo([0.001, 0.001], 50., x, epsilonGuess)

FPW_demo([0.001, 0.001], 50., x, [60.]) ## epsilon > 50
'''
When epsilon > 50, the finite potential well cannot contain the particle as the particle is more energetic than the well depth. This results in a free state that does not decay as x tends to +/- infinity. 
'''
        