from fastoad.model_base.scramjet_propulsion.scramjet import scramjet
from fastoad.model_base.aerothermodynamics.aerodynamics import WingGeometry, Aerodynamics
from fastoad.model_base.atmosphere import Atmosphere

a = Atmosphere(100000)

input_vals = {'M_freestream': 8, 'p_freestream':a.pressure, 'r_freestream':a.density, 'T_freestream':a.temperature, 'gamma_inlet':1.401, 'R':287.05, 'r_fuel':42, 'v_fuel':209.87, 'T_fuel':150, 'cp_inlet':1006, 'cp_fuel':12820, 'cp_combustor':1287.5, 'ER':0.029, 'hf':120e6, 'gamma_combustor':1.2869, 'gamma_outlet':1.31, 'alpha':0}
scale_factors = {'x_scale':1, 'y_scale':.8}

data = scramjet(scale_factors, input_vals)
data.run()

length = data.x_coordinates[-1]
newWing =  WingGeometry(length, 70)
aero = Aerodynamics(data, newWing, 100000, 8, 1.2)

aero.calculateWingLift(newWing, newWing.leadingAngle)
aero.calculateForebodyPressure(newWing, newWing.leadingAngle)
aero.calculateAfterbodyPressure(newWing, newWing.leadingAngle)
aero.calculateEngineLiftDrag(data, 6)
aero.calculateBodyBluntnessDrag(newWing)
aero.calculateInducedDrag(newWing, newWing.leadingAngle)
aero.calculateTotalDrag(newWing, newWing.leadingAngle)

print(aero.wingDrag[newWing.leadingAngle] + aero.engineDrag)