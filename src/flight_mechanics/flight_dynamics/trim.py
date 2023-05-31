import openmdao as om
import numpy as np
from math import sin
from math import asin
from math import sqrt

# import atmospheric parameters
# import intertia tensor
# import force and moment coefficients

# Parameters TODO: will be changed for when looping
# intertias defined as per mass [kg/m^2/kg]
inertia = ("xx": 2342, "yy": 2344, "zz": 3454, "xz": 50, "zx": 50)

# vehile properties
vehicle_mass = 150*10**3
vehicle_planform_area = 1500
vehicle_chord = 40

# general values
rho = 0.5
gravity = 9.81

# passed on values
force_coefficients = ("D": 0.3, "N": 0.1, "L": 0.5)
moment_coefficients = ("l": 0.2, "m": 0.1, "n": 0.4)
thrust = ("x": 300*10**3, "y": 0, "z": 0)
position_center_of_gravity = ("x": 0, "y": 0, "z": 0)
position_aerodynamic_center = ("x": -5, "y": 0, "z": 1)
position_center_of_thrust = ("x": -30, "y": 0, "z": -2)



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


            ## TODO: add the functions
            u_dot = -g*sin(theta)+1/vehicle_mass*(1/2*rho*sqrt(u**2+v**2+w**2)*vehicle_planform_area *(-force_coefficients.D*asin(w/u)))
            u_dot = 
            v_dot = 
            w_dot = 
            fi_dot = 
            theta_dot = 
            psi_dot = 
            p_dot = 
            q_dot = 
            r_dot = 


            y[0] = u_dot
            y[1] = v_dot
            y[1] = w_dot
            y[1] = fi_dot
            y[1] = theta_dot
            y[1] = psi_dot
            y[1] = p_dot
            y[1] = q_dot
            y[1] = r_dot



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
        prob['x'] = np.zeros(9)

        prob.driver = om.pyOptSparseDriver()
        prob.driver.options['optimizer'] = 'SLSQP'

        trim = [0, 0, 0, 0, 0, 0]
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
