from chamber import chamber_solver
from isentropic_expansion import M_2
from shocks import oblique, pm_expansion


#Idea is to have this parent class take an input of a dictionary with all parameters needed. It will call the functions created
# previously and perform the necessary calculations to find thrust and everything else. Perhaps then we can have a child class
# which is the interface between this and the other WPs. 

#List of param keys used: M_freestream, p_freestream, r_freestream, T_freestream, theta_1, gamma_inlet, theta_2, theta_3, R, r_fuel, v_fuel, T_fuel, cp_inlet, cp_fuel, cp_combustor, ER, hf, gamma_combustor, theta_outlet, gamma_outlet, A


class scramjet:
    def __init__(self, input_dict):
        self.params = input_dict

    def inlet(self):
        #First deflection
        M2, beta2, p2, r2, T2 = oblique(self.params['M_freestream'], self.params['p_freestream'], self.params['r_freestream'], self.params['T_freestream'], self.params['theta_1'], self.params['gamma_inlet'])
        #Second deflection
        M3, beta3, p3, r3, T3 = oblique(M2, p2, r2, T2, self.params['theta_2'], self.params['gamma_inlet'])
        #Third deflection, into inlet, assumed to be total input to mixer
        M4, beta4, p4, r4, T4 = oblique(M3, p3, r3, T3, self.params['theta_3'], self.params['gamma_inlet'])

        # a_out = (self.params['gamma_inlet']*self.params['R']*T4)**0.5
        # v_out = a_out*M4
        # mdot_out = r4*v_out*self.params['A']

        self.inlet_mach = M4
        self.inlet_pressure = p4
        self.inlet_density = r4
        self.inlet_temperature = T4
        a = (self.params['gamma_inlet'] * self.params['R'] * T4)**0.5
        vin = M4*a
        self.inlet_velocity = vin
        mdot_inlet = r4*self.params['A']*vin
        self.inlet_mdot = mdot_inlet

    def combustor(self):
        p2, v2, r2, T2 = chamber_solver(self.inlet_mach, self.inlet_pressure, self.inlet_density, self.inlet_temperature, self.params['r_fuel'], self.params['v_fuel'], self.params['T_fuel'], self.params['cp_inlet'], self.params['cp_fuel'], self.params['cp_combustor'], self.params['A'], self.params['ER']*self.inlet_mdot, self.params['hf'], self.params['gamma_combustor'], self.params['R'] )
        self.outlet_pressure = p2
        self.outlet_velocity = v2
        self.outlet_density = r2
        self.outlet_temperature = T2
        a = (self.params['gamma_combustor'] * self.params['R'] * T2)**0.5
        M2 = v2 / a
        self.outlet_mach = M2
        self.mdot_fuel = self.params['ER']*self.inlet_mdot

    def nozzle(self):
        M2, p2, r2, T2 = pm_expansion(self.outlet_mach, self.outlet_pressure, self.outlet_density, self.outlet_temperature, self.params['theta_outlet'], self.params['gamma_combustor'])
        Mout, Tout, pout, rout = M_2(p2, self.params['p_freestream'], M2, T2, self.params['gamma_outlet'], self.params['R'])
        self.exhaust_mach = M2
        self.exhaust_temperature = Tout
        self.exhaust_pressure = pout
        self.exhaust_density = rout
        a = (self.params['gamma_outlet'] * self.params['R'] * Tout)**0.5
        vout = Mout*a
        self.exhaust_velocity = vout

    def thrust(self):
        mdot_inlet = self.inlet_density*self.params['A']*self.inlet_velocity
        a_freestream = (self.params['gamma_inlet'] * self.params['R'] * self.params['T_freestream'])**0.5
        v_freestream = self.params['M_freestream']*a_freestream
        F = (mdot_inlet + self.mdot_fuel)*self.exhaust_velocity - mdot_inlet*v_freestream
        self.engine_thrust = F

