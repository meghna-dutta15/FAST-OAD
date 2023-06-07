from scramjet import scramjet

input_dict = {'M_freestream':8 , 'p_freestream':1090.16, 'r_freestream':0.0167207, 'T_freestream':227.130, 'theta_1':0.09806, 'gamma_inlet':1.401, 'theta_2':0.23117, 'theta_3':0.32923, 'R':287.05, 'r_fuel':42, 'v_fuel':209.87, 'T_fuel':150, 'cp_inlet':1006, 'cp_fuel':12820, 'cp_combustor':1287.5, 'ER':0.029, 'hf':120e6, 'gamma_combustor':1.2869, 'theta_outlet':0.30238, 'gamma_outlet':1.31, 'A':0.05393}

data = scramjet(input_dict)
data.inlet()
data.combustor()
data.nozzle()
data.thrust()

print(data.inlet_temperature)
print(data.inlet_density)
print(data.inlet_mach)
print(data.inlet_velocity)
print(data.inlet_pressure)
print(data.inlet_mdot)
print('')
print(data.outlet_density)
print(data.outlet_mach)
print(data.outlet_pressure)
print(data.outlet_temperature)
print(data.outlet_velocity)
#print('')
#print(data.engine_thrust)
