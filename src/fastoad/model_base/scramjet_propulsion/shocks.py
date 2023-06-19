import numpy as np
from scipy.optimize import root_scalar as root

def oblique(M1, p1, r1, T1, theta, gamma=1.4):
    
    def beta_callable(beta):
        zero = 2/np.tan(beta)*((M1**2)*(np.sin(beta)**2)-1)/((M1**2)*(gamma+np.cos(2*beta))+2) - np.tan(theta)
        return zero
    
    beta = root(beta_callable, x0 = np.pi/6, x1 = np.pi/3).root

    p2 = p1*(1 + 2*gamma/(gamma+1)*((M1**2)*(np.sin(beta)**2)-1))
    r2 = r1*(gamma+1)*(M1**2)*(np.sin(beta)**2) / ((gamma-1)*(M1**2)*(np.sin(beta)**2)+2)
    T2 = T1*(p2/p1)*(r1/r2)

    M2 = (1/np.sin(beta-theta)) * ( (1 + (gamma-1)/2 * (M1**2) * (np.sin(beta)**2)) / (gamma*(M1**2)*(np.sin(beta)**2) -(gamma-1)/2 )  )**0.5


    return M2, beta, p2, r2, T2

# M2, beta2, p2, r2, T2 = oblique(5, 7172, 0.115, 217, np.pi/6)

# theta2 = 3*np.pi/12 - np.pi/6

# M3, beta3, p3, r3, T3 = oblique(M2, p2, r2, T2, theta2)

# print(M3, beta3, p3, r3, T3)


def pm_expansion(M1, p1, r1, T1, theta, gamma=1.4):
    def prandtl_meyer(M):
        vM = (((gamma+1)/(gamma-1))**0.5) * np.arctan( ((gamma-1)/(gamma+1)*((M**2)-1))**0.5 ) - np.arctan(((M**2)-1)**0.5)
        return vM

    vM2 = theta + prandtl_meyer(M1)

    def prandtl_meyer_zero(M):
        returnable = (((gamma+1)/(gamma-1))**0.5) * np.arctan( ((gamma-1)/(gamma+1)*((M**2)-1))**0.5 ) - np.arctan(((M**2)-1)**0.5) - vM2
        return returnable
    
    M2 = root(prandtl_meyer_zero, x0 = 2, x1=3).root

    T2 = T1 * ((1+(gamma-1)/2*(M1**2)) / (1+(gamma-1)/2*(M2**2)))
    p2 = p1 * ((1+(gamma-1)/2*(M1**2)) / (1+(gamma-1)/2*(M2**2))) ** (gamma/(gamma-1))
    r2 = r1 * ((1+(gamma-1)/2*(M1**2)) / (1+(gamma-1)/2*(M2**2))) ** (1/(gamma-1))

    return M2, p2, r2, T2

# Mout, pout, rout, Tout = expansion(M3, p3, r3, T3, np.pi/12)
# Mout, pout, rout, Tout = pm_expansion(3, 300000, 0.05, 2500, np.pi/12)
# print(Mout)
# print(pout)
# print(Tout)


###########################
# Functions for use in system solver
