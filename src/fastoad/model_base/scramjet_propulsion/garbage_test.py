import numpy as np
from scipy.optimize import root_scalar as root

M1 = 5
p1 = 7172
gamma = 1.4
theta = np.pi/6


def beta_callable(beta):
    zero = 2/np.tan(beta)*((M1**2)*(np.sin(beta)**2)-1)/((M1**2)*(gamma+np.cos(2*beta))+2) - np.tan(theta)
    return zero
    
beta = root(beta_callable, x0 = np.pi/6, x1 = np.pi/3)
print(beta.root)
print(np.pi/6)



