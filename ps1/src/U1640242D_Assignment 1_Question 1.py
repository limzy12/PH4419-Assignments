### 1D Schrodinger Equation ###

from scipy import *
import matplotlib.pyplot as plt

import scipy.integrate as integrate
import scipy.linalg as la

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

##############
## Part (b) ##
##############
a = 1.
b = 4.

### Define the potential function
def potentialWell(x, V0 = 50.):
    if abs(x - b / 2) <= a / 2:
        return 0
    elif x < 0 or x > b:
        return inf
    else:
        return V0

### Define the perturbation integrand
def integrand(x, m, n):
    psi_m = sqrt(2 / b) * sin(m * pi * x / b) ## No need for complex conjugate since it is a real value
    psi_n = sqrt(2 / b) * sin(n * pi * x / b) 
    return psi_m * potentialWell(x) * psi_n

### Define a function to calculate the potential matrix element V_mn
def Potential(m, n):
    V_mn = integrate.quad(integrand, 0, b, (m, n))
    return V_mn[0]

### Define a function to return the k smallest eigenvalues
def eigenvalue_solver(A, N, k):
    eigvals = la.eigvalsh(A) ## Use eigvalsh since the matrix will be an observable which is Hermitian
    eigvals.sort()
    return eigvals[0:k]

# Set up the Hamilitonian matrix
N = array(range(1,50))
E_n = (N ** 2) * (pi ** 2) / (2 * (b ** 2))
A = diag(E_n, 0)
for i in range(len(N)): 
    for j in range(len(N)):
        A[i,j] += Potential(i+1,j+1)

# Print the first 4 eigenvalues
print(eigenvalue_solver(A, len(N), 4))

'''
To improve the accuracy of the numerical eigenvalues, we can consider expanding the size of the Hamilitonian matrix. 
'''