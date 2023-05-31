import openmdao as om
import numpy as np

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
            dy1_dt = x[0] * y[0] + x[1] * y[1] + x[2] * y[2]
            dy2_dt = x[3] * y[3] + x[4] * y[4] + x[5] * y[5]
            dy3_dt = x[6] * y[6] + x[7] * y[7] + x[8] * y[8]


            y[0] = dy1_dt
            y[1] = dy2_dt
            y[2] = dy3_dt


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
