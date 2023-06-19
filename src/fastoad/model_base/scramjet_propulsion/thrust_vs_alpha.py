from parametrize_scramjet import parametrize_scramjet
import numpy as np
import matplotlib.pyplot as plt

#input_dict = {'M_freestream': 8, 'p_freestream':1090.16, 'r_freestream':0.0167207, 'T_freestream':227.130, 'theta_1':0.09806, 'gamma_inlet':1.401, 'theta_2':0.23117, 'theta_3':0.32923, 'R':287.05, 'r_fuel':42, 'v_fuel':209.87, 'T_fuel':150, 'cp_inlet':1006, 'cp_fuel':12820, 'cp_combustor':1287.5, 'ER':0.029, 'hf':120e6, 'gamma_combustor':1.2869, 'theta_outlet':0.30238, 'gamma_outlet':1.31, 'A':0.05393}
#input_dict = {'M_freestream': 8, 'p_freestream':1090.16, 'r_freestream':0.0167207, 'T_freestream':227.130, 'theta_1':0.09806, 'gamma_inlet':1.401, 'theta_2':0.23117, 'theta_3':0.32923, 'R':287.05, 'r_fuel':42, 'v_fuel':209.87, 'T_fuel':150, 'cp_inlet':1006, 'cp_fuel':12820, 'cp_combustor':1287.5, 'ER':0.029, 'hf':120e6, 'gamma_combustor':1.2869, 'theta_outlet':0.30238, 'gamma_outlet':1.31, 'A':0.25}
input_vals = {'M_freestream': 8, 'p_freestream':1090.16, 'r_freestream':0.0167207, 'T_freestream':227.130, 'gamma_inlet':1.401, 'R':287.05, 'r_fuel':42, 'v_fuel':209.87, 'T_fuel':150, 'cp_inlet':1006, 'cp_fuel':12820, 'cp_combustor':1287.5, 'ER':0.029, 'hf':120e6, 'gamma_combustor':1.2869, 'gamma_outlet':1.31, 'alpha':0}
scale_factors = {'x_scale':1, 'y_scale':1}

alpha_array = np.linspace(-0.165,0.08,num=100)

data = parametrize_scramjet(scale_factors, input_vals, alpha_array=alpha_array)
data2 = data.parametrize_alpha()

thrust_array2 = []

for i in data2:
    thrust_array2.append(i.engine_thrust)

#print(thrust_array2)

thrust_array2 = np.array(thrust_array2)


plt.style.use('bmh')
plt.xlabel('Angle of Attack (rad)', fontsize = 15)
plt.ylabel('Net Thrust (kN)', fontsize = 15)
plt.title('Thrust vs \u03B1', fontsize = 15)
plt.plot(alpha_array, thrust_array2/1000)
plt.show()



