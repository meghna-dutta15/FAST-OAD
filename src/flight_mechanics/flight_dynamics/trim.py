import openmdao as om
import numpy as np
from math import sin
from math import asin
from math import cos
from math import tan
from math import atan
from math import sec
from math import sqrt

# import atmospheric parameters
# import force and moment coefficients
# import weight distribution parameters (inertia tensor, center of gravity, vehicle mass, and more)

# General thoughts:
# - Where is the origin of the aircraft?
# - The functions must be simplified, i.e. make matrices and performs matrix multiplications etc, if possible
# - Nicely import relevant parameters and their value from other WPs and name them, how?

# Parameters TODO: will be changed for when looping
# intertias defined as per mass [kg/m^2/kg]
inertia = ("xx": 2342, "yy": 2344, "zz": 3454, "xz": 50, "zx": 50)

# vehile properties
vehicle_mass = 150*10**3
vehicle_planform_area = 1500
vehicle_chord = 40

# general values
rho = 0.5
g = 9.81

# passed on values to this module
force_coefficients = {"D": 0.3, "N": 0.1, "L": 0.5}
moment_coefficients = {"l": 0.2, "m": 0.1, "n": 0.4}
thrust = {"x": 300*10**3, "y": 0, "z": 0}
position_center_of_gravity = {"x": 0, "y": 0, "z": 0}
position_aerodynamic_center = {"x": -5, "y": 0, "z": 1}
position_center_of_thrust = {"x": -30, "y": 0, "z": -2}
angles = {"phi", "theta", "psi"}
velocity = {"u": 0, "v": 0, "w": 0}


class MatrixSystem(om.Group):
    def setup(self):

        # Input Vector
        # x = [u, v, w, phi, theta, psi, p, q, r]
        self.add_subsystem('inputs', om.IndepVarComp('x', shape=9), promotes=['*'])
        self.add_subsystem('diff_eqs', DiffEqs())
        self.connect('x', 'diff_eqs.x')
        
        self.add_subsystem('diff_eqs_group', om.Group())
        self.diff_eqs_group.add_subsystem('diff_eqs_comp', self.diff_eqs)
        self.add_subsystem('diff_eqs_group', self.diff_eqs_group)

        self.add_subsystem('objective', om.ObjectiveComponent(), promotes=['*'])

class DiffEqs(om.ExplicitComponent):
        def setup(self):
            self.add_input('x', shape=9)
            self.add_output('y', shape=9)
        
        def compute(self, inputs, outputs):
            x = inputs['x']
            y = outputs['y']
            l = moment_coefficients['l']

            ## TODO: add the functions
            u_dot = -g*sin(theta)+1/vehicle_mass*(1/2*rho*sqrt(u**2+w**2)*vehicle_planform_area*(-force_coefficients.D*(w/u)) + thrust.x)
            w_dot = g*cos(phi)*cos(theta)+1/vehicle_mass*1/2*rho*sqrt(u**2+w**2)*vehicle_planform_area*(-force_coefficients.L*atan(w/u))
            q_dot = 1/vehicle_mass*(inertia.yy*1/2*rho*sqrt(u**2+w**2)*vehicle_chord*vehicle_planform_area* \
            (8.1576*10**(-8)*(w/u)**4 - 3.3297*10**(-6)*(w/u)**3 - 9.7364*10**(-6)*(w/u)**2 - 0.00051529*(w/u) + 0.00016104) \
            + 1/2*rho*sqrt(u**2+w**2)*vehicle_planform_area*((position_aerodynamic_center.x-position_center_of_gravity.x)*force_coefficients.L \
            - (position_aerodynamic_center.z-position_center_of_gravity.z)*force_coefficients.D) *(position_center_of_thrust.z-position_center_of_gravity.z) * thrust.x)
            ## angular velocities horizontal to the earth frame

            # theta_dot = q
            # h_dot = 
            

            # v_dot = g*sin(phi)*cos(theta)+1/vehicle_mass*1/2*rho*sqrt(u**2+v**2+w**2)*vehicle_planform_area*(-force_coefficients.N*asin(v/sqrt(u**2+v**2+w**2)))
            # fi_dot = p + q*sin(phi)*tan(theta) + r*cos(phi)*tan(theta)
            # theta_dot = q*cos(phi) - r*sin(phi)
            # psi_dot = q*sin(phi)*sec(theta) + r*cos(phi)*sec(theta)
            # p_dot = inertia.xx*1/vehicle_mass* ( 1/2*rho*sqrt(u**2+v**2+w**2)*vehicle_chord*vehicle_planform_area*moment_coefficients.l*atan(w/u) + 1/2*rho*sqrt(u**2+v**2+w**2)*vehicle_planform_area*((position_aerodynamic_center.z-position_center_of_gravity.z)*force_coefficients.N*asin(v/sqrt(u**2+v**2+w**2)) - (position_aerodynamic_center.y-position_center_of_gravity.y)*force_coefficients.L*atan(w/u))) + inertia.xz*1/vehicle_mass* ( 1/2*rho*sqrt(u**2+v**2+w**2)*vehicle_chord*vehicle_planform_area*moment_coefficients.n*atan(w/u) + 1/2*rho*sqrt(u**2+v**2+w**2)*vehicle_planform_area*((position_aerodynamic_center.y-position_center_of_gravity.y)*force_coefficients.D*atan(w/u) - (position_aerodynamic_center.x-position_center_of_gravity.x)*force_coefficients.N*asin(v/sqrt(u**2+v**2+w**2))) - (position_center_of_thrust.y-position_center_of_gravity.y)*thrust.x)
            # q_dot = inertia.yy*1/vehicle_mass* ( 1/2*rho*sqrt(u**2+v**2+w**2)*vehicle_chord*vehicle_planform_area*moment_coefficients.m*asin(v/sqrt(u**2+v**2+w**2)) + 1/2*rho*sqrt(u**2+v**2+w**2)*vehicle_planform_area*(-(position_aerodynamic_center.z-position_center_of_gravity.z)*force_coefficients.D*atan(w/u) + (position_aerodynamic_center.x-position_center_of_gravity.x)*force_coefficients.L*atan(w/u)) + (position_center_of_thrust.z-position_center_of_gravity.z)*thrust.x)
            # r_dot = inertia.xz*1/vehicle_mass* ( 1/2*rho*sqrt(u**2+v**2+w**2)*vehicle_chord*vehicle_planform_area*moment_coefficients.l*atan(w/u) + 1/2*rho*sqrt(u**2+v**2+w**2)*vehicle_planform_area*((position_aerodynamic_center.z-position_center_of_gravity.z)*force_coefficients.N*asin(v/sqrt(u**2+v**2+w**2)) - (position_aerodynamic_center.y-position_center_of_gravity.y)*force_coefficients.L*atan(w/u))) + inertia.zz*1/vehicle_mass* ( 1/2*rho*sqrt(u**2+v**2+w**2)*vehicle_chord*vehicle_planform_area*moment_coefficients.n*atan(w/u) + 1/2*rho*sqrt(u**2+v**2+w**2)*vehicle_planform_area*((position_aerodynamic_center.y-position_center_of_gravity.y)*force_coefficients.D*atan(w/u) - (position_aerodynamic_center.x-position_center_of_gravity.x)*force_coefficients.N*asin(v/sqrt(u**2+v**2+w**2))) - (position_center_of_thrust.y-position_center_of_gravity.y)*thrust.x)


            y[0] = u_dot
            y[1] = w_dot
            y[2] = q_dot
            y[3] = h_dot
            


class ObjectiveComponent(om.ExplicitComponent):
    def setup(self):
        self.add_input('y', shape=9)
        self.add_output('obj')

    def compute(self, inputs, outputs):
        y = inputs['y']
        outputs['obj'] = y

    def runSolver():
        prob = om.Problem()
        prob.model = om.NonlinearSystem()

        # Set the initial values for the variables
        prob['x'] = np.zeros(4)

        prob.driver = om.pyOptSparseDriver()
        prob.driver.options['optimizer'] = 'SLSQP'

        trim = [0, 0, 0, 0]
        for i, target in enumerate(trim):
            prob.model.add_objective('objective.obj[{}]'.format(i), scaler=1.0, target=target)

        # Run the optimization
        prob.setup()
        prob.run_model()
        prob.run_driver()

        # Print the results
        print('Optimal solution:')
        print('x:', prob['x'])
        print('y:', prob['y'])
        print('Objective:', prob['obj'])
