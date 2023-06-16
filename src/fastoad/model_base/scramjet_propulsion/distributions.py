from scramjet import scramjet
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d as interp
import numpy as np

#input_dict = {'M_freestream': 8, 'p_freestream':1090.16, 'r_freestream':0.0167207, 'T_freestream':227.130, 'theta_1':0.09806, 'gamma_inlet':1.401, 'theta_2':0.23117, 'theta_3':0.32923, 'R':287.05, 'r_fuel':42, 'v_fuel':209.87, 'T_fuel':150, 'cp_inlet':1006, 'cp_fuel':12820, 'cp_combustor':1287.5, 'ER':0.029, 'hf':120e6, 'gamma_combustor':1.2869, 'theta_outlet':0.30238, 'gamma_outlet':1.31, 'A':0.05393}
#input_dict = {'M_freestream': 8, 'p_freestream':1090.16, 'r_freestream':0.0167207, 'T_freestream':227.130, 'theta_1':0.09806, 'gamma_inlet':1.401, 'theta_2':0.23117, 'theta_3':0.32923, 'R':287.05, 'r_fuel':42, 'v_fuel':209.87, 'T_fuel':150, 'cp_inlet':1006, 'cp_fuel':12820, 'cp_combustor':1287.5, 'ER':0.029, 'hf':120e6, 'gamma_combustor':1.2869, 'theta_outlet':0.30238, 'gamma_outlet':1.31, 'A':0.25}
input_vals = {'M_freestream': 8, 'p_freestream':1090.16, 'r_freestream':0.0167207, 'T_freestream':227.130, 'gamma_inlet':1.401, 'R':287.05, 'r_fuel':42, 'v_fuel':209.87, 'T_fuel':150, 'cp_inlet':1006, 'cp_fuel':12820, 'cp_combustor':1287.5, 'ER':0.029, 'hf':120e6, 'gamma_combustor':1.2869, 'gamma_outlet':1.31, 'alpha':0}
scale_factors = {'x_scale':1, 'y_scale':1}

data = scramjet(scale_factors, input_vals)
data.run()

x_coordinates = data.x_coordinates
y_coordinates = data.y_coordinates
ramp1_vals = data.ramp1_vals
ramp2_vals = data.ramp2_vals
inlet_vals = data.inlet_vals
outlet_vals = data.outlet_vals
fan_vals = data.fan_vals
exhaust_vals = data.exhaust_vals

plt.style.use('bmh')
plt.plot([x_coordinates[0], x_coordinates[0]], [data.params['p_freestream']/1000, ramp1_vals['p']/1000], color='black', linestyle='dashed')
plt.plot([x_coordinates[0], x_coordinates[1]], [ramp1_vals['p']/1000, ramp1_vals['p']/1000], label='Ramp 1')
plt.plot([x_coordinates[1], x_coordinates[1]], [ramp1_vals['p']/1000, ramp2_vals['p']/1000], color='black', linestyle='dashed')

plt.plot([x_coordinates[1], x_coordinates[2]], [ramp2_vals['p']/1000, ramp2_vals['p']/1000], label='Ramp 2')
plt.plot([x_coordinates[2], x_coordinates[2]], [ramp2_vals['p']/1000, inlet_vals['p']/1000], color='black', linestyle='dashed')

plt.plot([x_coordinates[2], (x_coordinates[4]-x_coordinates[2])/5+x_coordinates[2]], [inlet_vals['p']/1000, inlet_vals['p']/1000], label='Inlet')
plt.plot([(x_coordinates[4]-x_coordinates[2])/5+x_coordinates[2], (x_coordinates[4]-x_coordinates[2])/5+x_coordinates[2]], [inlet_vals['p']/1000, outlet_vals['p']/1000], color='black', linestyle='dashed')

plt.plot([(x_coordinates[4]-x_coordinates[2])/5+x_coordinates[2], x_coordinates[4]], [outlet_vals['p']/1000, outlet_vals['p']/1000], label='Combustor')
# plt.plot([x_coordinates[4], x_coordinates[4]], [outlet_vals['p'], exhaust_vals['p']], color='black', linestyle='dashed')

# plt.plot([x_coordinates[4], x_coordinates[5]], [exhaust_vals['p'], exhaust_vals['p']], label='Nozzle')

nozzle_x = [x_coordinates[4], (x_coordinates[5]-x_coordinates[4])/10+x_coordinates[4], x_coordinates[5]]
nozzle_y = [outlet_vals['p']/1000, fan_vals['p']/1000, exhaust_vals['p']/1000]
nozzle_func = interp(nozzle_x, nozzle_y, kind = 'linear')

nozzle_x_dense = np.linspace(x_coordinates[4], x_coordinates[5], num=100)
nozzle_y_dense = nozzle_func(nozzle_x_dense)

# nozzle_func = interp(nozzle_x, nozzle_y, kind = 'quadratic')

# nozzle_x_dense = np.linspace(x_coordinates[4], x_coordinates[5], num=100)
# nozzle_y_dense = nozzle_func(nozzle_x_dense)

plt.plot(nozzle_x_dense, nozzle_y_dense, label = 'Nozzle')

plt.xlabel('X Position Along Engine', fontsize = 15)
plt.ylabel('Static Pressure (kPa)', fontsize = 15)
plt.title('Pressure Distribution Throughout Engine', fontsize = 15)

plt.legend()
plt.show()


###########################################################
plt.style.use('bmh')
plt.plot([x_coordinates[0], x_coordinates[0]], [data.params['T_freestream'], ramp1_vals['T']], color='black', linestyle='dashed')
plt.plot([x_coordinates[0], x_coordinates[1]], [ramp1_vals['T'], ramp1_vals['T']], label='Ramp 1')
plt.plot([x_coordinates[1], x_coordinates[1]], [ramp1_vals['T'], ramp2_vals['T']], color='black', linestyle='dashed')

plt.plot([x_coordinates[1], x_coordinates[2]], [ramp2_vals['T'], ramp2_vals['T']], label='Ramp 2')
plt.plot([x_coordinates[2], x_coordinates[2]], [ramp2_vals['T'], inlet_vals['T']], color='black', linestyle='dashed')

plt.plot([x_coordinates[2], (x_coordinates[4]-x_coordinates[2])/5+x_coordinates[2]], [inlet_vals['T'], inlet_vals['T']], label='Inlet')
plt.plot([(x_coordinates[4]-x_coordinates[2])/5+x_coordinates[2], (x_coordinates[4]-x_coordinates[2])/5+x_coordinates[2]], [inlet_vals['T'], outlet_vals['T']], color='black', linestyle='dashed')

plt.plot([(x_coordinates[4]-x_coordinates[2])/5+x_coordinates[2], x_coordinates[4]], [outlet_vals['T'], outlet_vals['T']], label='Combustor')
# plt.plot([x_coordinates[4], x_coordinates[4]], [outlet_vals['p'], exhaust_vals['p']], color='black', linestyle='dashed')

# plt.plot([x_coordinates[4], x_coordinates[5]], [exhaust_vals['p'], exhaust_vals['p']], label='Nozzle')

nozzle_x = [x_coordinates[4], (x_coordinates[5]-x_coordinates[4])/10+x_coordinates[4], x_coordinates[5]]
nozzle_y = [outlet_vals['T'], fan_vals['T'], exhaust_vals['T']]
nozzle_func = interp(nozzle_x, nozzle_y, kind = 'linear')

nozzle_x_dense = np.linspace(x_coordinates[4], x_coordinates[5], num=100)
nozzle_y_dense = nozzle_func(nozzle_x_dense)

# nozzle_func = interp(nozzle_x, nozzle_y, kind = 'quadratic')

# nozzle_x_dense = np.linspace(x_coordinates[4], x_coordinates[5], num=100)
# nozzle_y_dense = nozzle_func(nozzle_x_dense)

plt.plot(nozzle_x_dense, nozzle_y_dense, label = 'Nozzle')

plt.xlabel('X Position Along Engine', fontsize = 15)
plt.ylabel('Static Temperature (K)', fontsize = 15)
plt.title('Temperature Distribution Throughout Engine', fontsize = 15)

plt.legend()
plt.show()


###########################################################
plt.style.use('bmh')
plt.plot([x_coordinates[0], x_coordinates[0]], [data.params['M_freestream'], ramp1_vals['M']], color='black', linestyle='dashed')
plt.plot([x_coordinates[0], x_coordinates[1]], [ramp1_vals['M'], ramp1_vals['M']], label='Ramp 1')
plt.plot([x_coordinates[1], x_coordinates[1]], [ramp1_vals['M'], ramp2_vals['M']], color='black', linestyle='dashed')

plt.plot([x_coordinates[1], x_coordinates[2]], [ramp2_vals['M'], ramp2_vals['M']], label='Ramp 2')
plt.plot([x_coordinates[2], x_coordinates[2]], [ramp2_vals['M'], inlet_vals['M']], color='black', linestyle='dashed')

plt.plot([x_coordinates[2], (x_coordinates[4]-x_coordinates[2])/5+x_coordinates[2]], [inlet_vals['M'], inlet_vals['M']], label='Inlet')
plt.plot([(x_coordinates[4]-x_coordinates[2])/5+x_coordinates[2], (x_coordinates[4]-x_coordinates[2])/5+x_coordinates[2]], [inlet_vals['M'], outlet_vals['M']], color='black', linestyle='dashed')

plt.plot([(x_coordinates[4]-x_coordinates[2])/5+x_coordinates[2], x_coordinates[4]], [outlet_vals['M'], outlet_vals['M']], label='Combustor')
# plt.plot([x_coordinates[4], x_coordinates[4]], [outlet_vals['p'], exhaust_vals['p']], color='black', linestyle='dashed')

# plt.plot([x_coordinates[4], x_coordinates[5]], [exhaust_vals['p'], exhaust_vals['p']], label='Nozzle')

nozzle_x = [x_coordinates[4], (x_coordinates[5]-x_coordinates[4])/10+x_coordinates[4], x_coordinates[5]]
nozzle_y = [outlet_vals['M'], fan_vals['M'], exhaust_vals['M']]
nozzle_func = interp(nozzle_x, nozzle_y, kind = 'linear')

nozzle_x_dense = np.linspace(x_coordinates[4], x_coordinates[5], num=100)
nozzle_y_dense = nozzle_func(nozzle_x_dense)

# nozzle_func = interp(nozzle_x, nozzle_y, kind = 'quadratic')

# nozzle_x_dense = np.linspace(x_coordinates[4], x_coordinates[5], num=100)
# nozzle_y_dense = nozzle_func(nozzle_x_dense)

plt.plot(nozzle_x_dense, nozzle_y_dense, label = 'Nozzle')

plt.xlabel('X Position Along Engine', fontsize = 15)
plt.ylabel('Mach #', fontsize = 15)
plt.title('Flow Mach # Throughout Engine', fontsize = 15)

plt.legend()
plt.show()

