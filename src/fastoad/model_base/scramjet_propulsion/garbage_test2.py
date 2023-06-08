import numpy as np
from scipy.optimize import root_scalar


def get_rayleigh_pressure_ratio(mach_1, mach_2, gamma=1.4):
    '''Return p2/p1 for Rayleigh flow'''
    return ((1 + gamma*mach_1**2) / (1 + gamma*mach_2**2))

def get_rayleigh_temperature_ratio(mach_1, mach_2, gamma=1.4):
    '''Return T2/T1 for Rayleigh flow'''
    return (
        ((1 + gamma*mach_1**2) / (1 + gamma*mach_2**2))**2 *
        (mach_2**2 / mach_1**2)
        )
            
def get_rayleigh_density_ratio(mach_1, mach_2, gamma=1.4):
    '''Return rho2/rho1 for Rayleigh flow'''
    return (
        ((1 + gamma*mach_2**2) / (1 + gamma*mach_1**2)) *
        (mach_1**2 / mach_2**2)
        )
            
def get_rayleigh_stag_temperature_ratio(mach_1, mach_2, gamma=1.4):
    '''Return Tt2/Tt1 for Rayleigh flow'''
    return (
        ((1 + gamma*mach_1**2) / (1 + gamma*mach_2**2))**2 * 
            (mach_2 / mach_1)**2 * (
            (1 + 0.5*(gamma-1)*mach_2**2) /
            (1 + 0.5*(gamma-1)*mach_1**2)
            )
        )

def get_rayleigh_stag_pressure_ratio(mach_1, mach_2, gamma=1.4):
    '''Return pt2/pt1 for Rayleigh flow'''
    return (
        ((1 + gamma*mach_1**2) / (1 + gamma*mach_2**2)) * (
            (1 + 0.5*(gamma-1)*mach_2**2) /
            (1 + 0.5*(gamma-1)*mach_1**2)
            )**(gamma / (gamma - 1))
        )



def get_stagnation_temperature_ratio(mach, gamma=1.4):
    '''Calculates T/Tt stagnation relationship.'''
    return (1.0 / (1 + 0.5*(gamma - 1) * mach**2))



gamma = 1.4
R = 287
cp = 1000
q = 120

T1 = 220
p1 = 70
V1 = 122

a1 = np.sqrt(gamma * R * T1)
M1 = V1 / a1
print(M1)
T1_Tt1 = get_stagnation_temperature_ratio(M1, gamma)
Tt1 = T1 / T1_Tt1

delta_Tt = q / cp
Tt2 = Tt1 + delta_Tt



def solve_rayleigh_stagnation_temperature(M2, M1, Tt2_Tt1, gamma):
    '''Used to find unknown M2 given M1, stagnation temperature ratio, and gamma.'''
    return (
        Tt2_Tt1 - get_rayleigh_stag_temperature_ratio(M1, M2, gamma)
        )


Tt2_Tt1 = Tt2 / Tt1
root = root_scalar(
    solve_rayleigh_stagnation_temperature, 
    x0=0.5, x1=0.7, args=(M1, Tt2_Tt1, gamma)
    )
M2 = root.root
print(M2)

p2 = p1 * get_rayleigh_pressure_ratio(M1, M2, gamma)
T2 = T1 * get_rayleigh_temperature_ratio(M1, M2, gamma)
print(T2)
print(p2)