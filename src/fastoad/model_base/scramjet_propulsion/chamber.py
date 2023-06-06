import numpy as np
from scipy.optimize import fsolve

def chamber_solver(Min, pin, rin, Tin, rf, vf, Tf, cpin, cpf, cpc, A, mdotf, hf, gamma=1.4, R=287.05):
    
    if Min < 1:
        raise ValueError('Mach inside combustion chamber cannot be lower than 1')
    
    ain = (gamma*R*Tin)
    vin = Min*ain

    mdotin = rin*A*vin

    
    def solver_func(x):
        p2, v2, r2, T2 = x[0], x[1], x[2], x[3]

        zero1 = r2*v2 - rin*vin - rf*vf

        zero2 = p2 + r2*(v2**2) - pin - rin*(vin**2)

        Tt1 = Tin + (vin**2)/(2*cpin)
        Tt2 = T2 + (v2**2)/(2*cpc)
        Ttf = Tf + (vf**2)/(2*cpf)
        zero3 = (mdotin+mdotf)*cpc*Tt2 - mdotin*cpin*Tt1 - mdotf*hf - mdotf*cpf*Ttf
        
        zero4 = r2*R*T2 - p2

        return np.array([zero1, zero2, zero3, zero4])
    
    x = fsolve(solver_func, np.array([1000, 1000, 20, 500]))

    return x[0], x[1], x[2], x[3]



