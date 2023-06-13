from chamber import chamber_solver
from isentropic_expansion import M_2
from shocks import oblique, pm_expansion
from rayleigh import mach_2_calculation
import numpy as np


#Idea is to have this parent class take an input of a dictionary with all parameters needed. It will call the functions created
# previously and perform the necessary calculations to find thrust and everything else. Perhaps then we can have a child class
# which is the interface between this and the other WPs. 

#List of param keys used: M_freestream, p_freestream, r_freestream, T_freestream, theta_1, gamma_inlet, theta_2, theta_3, R, r_fuel, v_fuel, T_fuel, cp_inlet, cp_fuel, cp_combustor, ER, hf, gamma_combustor, theta_outlet, gamma_outlet, A


class scramjet_calculations():
    def __init__(self, scale_factors, input_vals):
        self.geo = scale_factors
        self.vals = input_vals

    def define_geometry(self):
        x_default = np.array([0, 2.68165, 4.05500, 3.75976, 5.65000, 8.00000])
        y_default = np.array([0, 0.26383, 0.73307, 1, 0.73307, 0])

        x = x_default*self.geo['x_scale']
        y = y_default*self.geo['y_scale']

        theta1 = np.arctan2(y[1], x[1])
        theta2 = np.arctan2(y[2]-y[1], x[2]-x[1]) - theta1
        theta3 = np.arctan2(y[2]-y[1], x[2]-x[1])
        theta_outlet = np.arctan2(y[4]-y[5], x[5]-x[4])
        A = y[3] - y[2]

        self.params = {'M_freestream': self.vals['M_freestream'], 'p_freestream':self.vals['p_freestream'], 'r_freestream':self.vals['r_freestream'], 'T_freestream':self.vals['T_freestream'], 'theta_1':theta1, 'gamma_inlet':self.vals['gamma_inlet'], 'theta_2':theta2, 'theta_3':theta3, 'R':self.vals['R'], 'r_fuel':self.vals['r_fuel'], 'v_fuel':self.vals['v_fuel'], 'T_fuel':self.vals['T_fuel'], 'cp_inlet':self.vals['cp_inlet'], 'cp_fuel':self.vals['cp_fuel'], 'cp_combustor':self.vals['cp_combustor'], 'ER':self.vals['ER'], 'hf':self.vals['hf'], 'gamma_combustor':self.vals['gamma_combustor'], 'theta_outlet':theta_outlet, 'gamma_outlet':self.vals['gamma_outlet'], 'A':A}

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
        # a_freestream = (self.params['gamma_inlet'] * self.params['R'] * self.params['T_freestream'])**0.5
        # v_freestream = self.params['M_freestream']*a_freestream
        # mdot_inlet = self.params['r_freestream']*(v_freestream)*0.787
        self.inlet_mdot = mdot_inlet

    def combustor(self):
        # p2, v2, r2, T2 = chamber_solver(self.inlet_mach, self.inlet_pressure, self.inlet_density, self.inlet_temperature, self.params['r_fuel'], self.params['v_fuel'], self.params['T_fuel'], self.params['cp_inlet'], self.params['cp_fuel'], self.params['cp_combustor'], self.params['A'], self.params['ER']*self.inlet_mdot, self.params['hf'], self.params['gamma_combustor'], self.params['R'] )
        # self.outlet_pressure = p2
        # self.outlet_velocity = v2
        # self.outlet_density = r2
        # self.outlet_temperature = T2
        # a = (self.params['gamma_combustor'] * self.params['R'] * T2)**0.5
        # M2 = v2 / a
        # self.outlet_mach = M2

        p2, M2, r2, T2 = mach_2_calculation(self.inlet_mach, self.inlet_pressure, self.inlet_density, self.inlet_temperature, self.params['r_fuel'], self.params['v_fuel'], self.params['T_fuel'], self.params['cp_inlet'], self.params['cp_fuel'], self.params['cp_combustor'], self.params['A'], self.params['ER']*self.inlet_mdot, self.params['hf'], self.params['gamma_combustor'], self.params['R'] )
        self.outlet_pressure = p2
        self.outlet_mach = M2
        self.outlet_density = r2
        self.outlet_temperature = T2
        a = (self.params['gamma_combustor'] * self.params['R'] * T2)**0.5
        self.outlet_velocity = M2*a

        self.mdot_fuel = self.params['ER']*self.inlet_mdot

    def nozzle(self):
        M2, p2, r2, T2 = pm_expansion(self.outlet_mach, self.outlet_pressure, self.outlet_density, self.outlet_temperature, self.params['theta_outlet'], self.params['gamma_combustor'])
        Mout, Tout, pout, rout = M_2(p2, self.params['p_freestream'], M2, T2, self.params['gamma_outlet'], self.params['R'])
        self.exhaust_mach = Mout
        self.exhaust_temperature = Tout
        self.exhaust_pressure = pout
        self.exhaust_density = rout
        a = (self.params['gamma_outlet'] * self.params['R'] * Tout)**0.5
        vout = Mout*a
        self.exhaust_velocity = vout

    def thrust(self):
        a_freestream = (self.params['gamma_inlet'] * self.params['R'] * self.params['T_freestream'])**0.5
        v_freestream = self.params['M_freestream']*a_freestream
        #F = self.exhaust_density*0.787*self.exhaust_velocity**2 - self.params['r_freestream']*(v_freestream**2)*0.787
        self.returnable = self.params['r_freestream']*(v_freestream)*0.787
        F = (self.inlet_mdot + self.mdot_fuel)*self.exhaust_velocity - self.inlet_mdot*v_freestream
        self.engine_thrust = F

class scramjet(scramjet_calculations):
    def __init__(self, scale_factors, input_vals):
      super().__init__(scale_factors, input_vals)

    def run(self):
        super().define_geometry()
        super().inlet()
        super().combustor()
        super().nozzle()
        super().thrust()



