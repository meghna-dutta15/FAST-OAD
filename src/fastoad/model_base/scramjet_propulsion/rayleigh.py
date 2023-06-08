from scipy.optimize import root_scalar

def mach_2_calculation(Min, pin, rin, Tin, rf, vf, Tf, cpin, cpf, cpc, A, mdotf, hf, gamma=1.4, R=287.05): 
    ain = (gamma*R*Tin)**0.5
    vin = Min*ain
    mdotin = rin*A*vin
    Ttf = Tf + (vf**2)/(2*cpf)
    T_t_1 = Tin + (vin**2)/(2*cpin)
    q = mdotf*hf + mdotf*cpf*Ttf
    #T_t_2 = (q/cpf) + T_t_1
    T_t_2 = (mdotin*cpin*T_t_1 + mdotf*hf  + mdotf*cpf*Ttf) / ((mdotin+mdotf)*cpc)

    def get_rayleigh_stag_temperature_ratio(mach_1, mach_2, gamma=1.4):
        '''Return Tt2/Tt1 for Rayleigh flow'''
        return (
            ((1 + gamma*mach_1**2) / (1 + gamma*mach_2**2))**2 * 
                (mach_2 / mach_1)**2 * (
                (1 + 0.5*(gamma-1)*mach_2**2) /
                (1 + 0.5*(gamma-1)*mach_1**2)
                )
            )


    def solve_rayleigh_stagnation_temperature(M2, M1, Tt2_Tt1, gamma):
        '''Used to find unknown M2 given M1, stagnation temperature ratio, and gamma.'''
        return (Tt2_Tt1 - get_rayleigh_stag_temperature_ratio(M1, M2, gamma))
    
    Tt2_Tt1 = T_t_2 / T_t_1
    root = root_scalar(solve_rayleigh_stagnation_temperature, x0=2, x1=3, args=(Min, Tt2_Tt1, gamma))
    M2 = root.root

    def get_rayleigh_pressure_ratio(mach_1, mach_2, gamma=1.4):
        '''Return p2/p1 for Rayleigh flow'''
        p_2_in = ((1 + gamma*mach_1**2) / (1 + gamma*mach_2**2))
        return p_2_in

    def get_rayleigh_temperature_ratio(mach_1, mach_2, gamma=1.4):
        '''Return T2/T1 for Rayleigh flow'''
        t_2_in = ( ((1 + gamma*mach_1**2) / (1 + gamma*mach_2**2))**2 * (mach_2**2 / mach_1**2) )
        return t_2_in

            
    def get_rayleigh_density_ratio(mach_1, mach_2, gamma=1.4):
        '''Return rho2/rho1 for Rayleigh flow'''
        r_2_in = (((1 + gamma*mach_2**2) / (1 + gamma*mach_1**2)) *(mach_1**2 / mach_2**2))
        return r_2_in


    def get_rayleigh_stag_pressure_ratio(mach_1, mach_2, gamma=1.4):
        '''Return pt2/pt1 for Rayleigh flow'''
        return (
            ((1 + gamma*mach_1**2) / (1 + gamma*mach_2**2)) * (
                (1 + 0.5*(gamma-1)*mach_2**2) /
                (1 + 0.5*(gamma-1)*mach_1**2)
                )**(gamma / (gamma - 1))
            )
    
    p_2_in = get_rayleigh_pressure_ratio(Min, M2, gamma)
    t_2_in = get_rayleigh_temperature_ratio(Min, M2, gamma)
    r_2_in = get_rayleigh_density_ratio(Min, M2, gamma)


    def rayleigh_flow(pin, rin, tin, p_2_in, t_2_in, r_2_in):
        p2 = p_2_in * pin
        t2 = t_2_in * tin
        r2 = r_2_in * rin
        return p2, t2, r2

    p2, T2, r2 = rayleigh_flow(pin, rin, Tin, p_2_in, t_2_in, r_2_in)
    
    return p2, M2, r2, T2