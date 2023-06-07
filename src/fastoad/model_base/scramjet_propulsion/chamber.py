import numpy as np
from scipy.optimize import root

def chamber_solver(Min, pin, rin, Tin, rf, vf, Tf, cpin, cpf, cpc, A, mdotf, hf, gamma=1.4, R=287.05):
    
    if Min < 1:
        raise ValueError('Mach inside combustion chamber cannot be lower than 1')
    
    ain = (gamma*R*Tin)**0.5
    vin = Min*ain

    mdotin = rin*A*vin

    
    def solver_func(x):
        p2, v2, r2, T2 = x[0], x[1], x[2], x[3]

        zero1 = r2*v2 - rin*vin #- rf*vf

        zero2 = p2 + r2*(v2**2) - pin - rin*(vin**2)

        Tt1 = Tin + (vin**2)/(2*cpin)
        Tt2 = T2 + (v2**2)/(2*cpc)
        Ttf = Tf + (vf**2)/(2*cpf)
        zero3 = (mdotin+mdotf)*cpc*Tt2 - mdotin*cpin*Tt1 - 0.5*mdotf*hf  - mdotf*cpf*Ttf
        
        zero4 = r2*R*T2 - p2

        return np.array([zero1, zero2, zero3, zero4])
    
    y = root(solver_func, np.array([1000, 1000, 20, 500]), method='hybr').x

    return y[0], y[1], y[2], y[3]

#print(chamber_solver(2, 1000, 0.1, 300, 42, 1000, 150, 1004, 1004, 1300, 0.1, 15, 120e6))

