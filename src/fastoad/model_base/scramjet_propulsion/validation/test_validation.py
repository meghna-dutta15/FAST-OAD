from scramjet_validation import scramjet

#input_dict = {'M_freestream': 8, 'p_freestream':1090.16, 'r_freestream':0.0167207, 'T_freestream':227.130, 'theta_1':0.09806, 'gamma_inlet':1.401, 'theta_2':0.23117, 'theta_3':0.32923, 'R':287.05, 'r_fuel':42, 'v_fuel':209.87, 'T_fuel':150, 'cp_inlet':1006, 'cp_fuel':12820, 'cp_combustor':1287.5, 'ER':0.029, 'hf':120e6, 'gamma_combustor':1.2869, 'theta_outlet':0.30238, 'gamma_outlet':1.31, 'A':0.05393}
#input_dict = {'M_freestream': 8, 'p_freestream':1090.16, 'r_freestream':0.0167207, 'T_freestream':227.130, 'theta_1':0.09806, 'gamma_inlet':1.401, 'theta_2':0.23117, 'theta_3':0.32923, 'R':287.05, 'r_fuel':42, 'v_fuel':209.87, 'T_fuel':150, 'cp_inlet':1006, 'cp_fuel':12820, 'cp_combustor':1287.5, 'ER':0.029, 'hf':120e6, 'gamma_combustor':1.2869, 'theta_outlet':0.30238, 'gamma_outlet':1.31, 'A':0.25}
input_vals = {'M_freestream': 8, 'p_freestream':1090.16, 'r_freestream':0.0167207, 'T_freestream':227.130, 'gamma_inlet':1.401, 'R':287.05, 'r_fuel':42, 'v_fuel':209.87, 'T_fuel':150, 'cp_inlet':1006, 'cp_fuel':12820, 'cp_combustor':1287.5, 'ER':0.029, 'hf':120e6, 'gamma_combustor':1.2869, 'gamma_outlet':1.31, 'alpha':0}
scale_factors = {'x_scale':1, 'y_scale':1}

data = scramjet(scale_factors, input_vals)
data.run()

print('Inlet Values:')
print(data.inlet_temperature)
print(data.inlet_density)
print(data.inlet_mach)
print(data.inlet_velocity)
print(data.inlet_pressure)
print(data.inlet_mdot)
print('')

print('Outlet Values:')
print(data.outlet_temperature)
print(data.outlet_density)
print(data.outlet_mach)
print(data.outlet_velocity)
print(data.outlet_pressure)

print('')
print('Exhaust Values:')
print(data.exhaust_temperature)
print(data.exhaust_density)
print(data.exhaust_mach)
print(data.exhaust_velocity)
print(data.exhaust_pressure)

print('')
print('Thrust:')
print(data.engine_thrust)
print(data.returnable)

print('')
print('')

cfd_p1 = 3160
spread_p1 = 3360
fast_p1 = data.ramp1_vals['p']
p1_diff = abs(fast_p1-cfd_p1)/cfd_p1*100

cfd_p2 = 17400
spread_p2 = 17800
fast_p2 = data.ramp2_vals['p']
p2_diff = abs(fast_p2-cfd_p2)/cfd_p2*100

cfd_pinlet = 102000 #Very approximate 
spread_pinlet = 102000
fast_pinlet = data.inlet_vals['p']
pinlet_diff = abs(fast_pinlet-cfd_pinlet)/cfd_pinlet*100

cfd_poutlet = 330000
spread_poutlet = 322000
fast_poutlet = data.outlet_vals['p']
poutlet_diff = abs(fast_poutlet-cfd_poutlet)/cfd_poutlet*100

print(fast_p1)
print(p1_diff)
print('')

print(fast_p2)
print(p2_diff)
print('')

print(fast_pinlet)
print(pinlet_diff)
print('')

print(fast_poutlet)
print(poutlet_diff)



