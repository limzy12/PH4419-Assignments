from scipy import *
import matplotlib.pyplot as plt

### We define a function that converts Cartesian coordinates to Cylindrical coordinates
def cart2cyl(r): # The input argument is a 3D vector r in Cartesian coordinates
    x = r[0]
    y = r[1]
    z = r[2]
    rho = sqrt(x**2 + y**2)
    theta = arctan(y / x)

    return [rho, theta, z]

vec = [1., 1., 0.]
print(cart2cyl(vec))