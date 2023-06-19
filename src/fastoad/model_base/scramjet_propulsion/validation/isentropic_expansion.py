import numpy as np
from scipy.optimize import root_scalar as root

# def M_2(P_1, P_2, M_1, T_1, gamma=1.4, R=287.05):
    
#     def zero_func(M_2):
#         zero = P_1*((1 + (((gamma-1)/2)*(M_1*M_1)))**(gamma/(gamma-1))) - P_2*((1 + (((gamma-1)/2)*(M_2*M_2)))**(gamma/(gamma-1)))
#         return zero
    
#     M_2 = root(zero_func, x0=2, x1=7).root

#     def zero_func(T_2):
#         zero = T_1*((1 + (((gamma-1)/2)*(M_1*M_1)))) - T_2*((1 + (((gamma-1)/2)*(M_2*M_2))))
#         return zero
    
#     T_2 = root(zero_func, x0=1, x1=4000).root

#     def zero_func(r_2):
#         zero = r_2 - (P_2/(R*T_2))
#         return zero
    
#     r_2 = root(zero_func, x0=0, x1=5000).root

#     return M_2, T_2, P_2, r_2


def M_2(P_1, P_2, M_1, T_1, gamma=1.4, R=287.05):
    gamma_1 = 2/(gamma-1)
    P_1_2 = (P_2/P_1)**((1-gamma)/gamma)
    M_2 = (((gamma_1)*(-(1-((P_1_2)*(1 + ((M_1*M_1)/(gamma_1))   )))))**(0.5))
    T_2 = (T_1*(1 + ((M_1*M_1)/(gamma_1)))) / ((1 + ((M_2*M_2)/(gamma_1))))
    r_2 = P_2 / (R*T_2)
    return M_2, T_2, P_2, r_2