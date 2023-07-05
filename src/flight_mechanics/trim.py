import openmdao.api as om
import numpy as np
from math import sin, cos, atan, sqrt, log 

# import atmospheric parameters
# import force and moment coefficients
# import weight distribution parameters (inertia tensor, center of gravity, vehicle mass, and more)

# General thoughts:
# - Where is the origin of the aircraft?
# - The functions must be simplified, i.e. make matrices and performs matrix multiplications etc, if possible
# - Nicely import relevant parameters and their value from other WPs and name them, how?

# Parameters TODO: will be changed for when looping
# intertias defined as per mass [kg/m^2/kg]
inertia = {"xx": 2342, "yy": 23440, "zz": 3454, "xz": 50, "zx": 50}

# vehile properties
vehicle_mass = 150*10**3
vehicle_planform_area = 1500
vehicle_chord = 40

# general values
rho = 0.5
g = 9.81

L = 30

# passed on values to this module
force_coefficients = {"D": 0.3, "N": 0.1, "L": 0.5}
moment_coefficients = {"l": 0.2, "m": 0.1, "n": 0.4}
thrust = {"x": 300*10**3, "y": 0, "z": 0}
position_center_of_gravity = {"x": 0, "y": 0, "z": 0}
position_aerodynamic_center = {"x": -5, "y": 0, "z": 1}
position_center_of_thrust = {"x": -30, "y": 0, "z": -2}
angles = {"phi":0,"theta":0, "psi":0}
init_altitude = 20000
velocity = {"u": 3000, "v": 0, "w": 0}



class DiffEqs(om.ExplicitComponent):
        def setup(self):
            
            self.add_input('theta')
            self.add_input('u')
            self.add_input('w')
            self.add_input('q')
            self.add_output('theta_dot')
            self.add_output('u_dot')
            self.add_output('w_dot')
            self.add_output('q_dot')
            #self.add_objective('objective')
            
            self.declare_partials(of=['u_dot', 'w_dot','objective'], wrt=['theta'], method='fd')
            self.declare_partials(of=['theta_dot','objective'], wrt=['q'], method='fd')
            self.declare_partials(of=['u_dot', 'w_dot', 'q_dot', 'objective'], wrt=['u', 'w'], method='fd')
            self.add_output('objective')
        
        def compute(self, inputs, outputs):
            theta = inputs['theta']
            u = inputs['u']
            w = inputs['w']
            q = inputs['q']

            ## TODO: add the functions
            #h_dot = sqrt(u**2+w**2)*sin(theta-(w/u))
            outputs['theta_dot'] = q
            outputs['u_dot'] = -g*np.sin(theta)+1/vehicle_mass*(1/2*rho*np.sqrt(u**2+w**2)*vehicle_planform_area*(-force_coefficients['D']*(w/u)) + thrust['x'])
            outputs['w_dot'] = g*np.cos(theta)+1/vehicle_mass*1/2*rho*np.sqrt(u**2+w**2)*vehicle_planform_area*(-force_coefficients['L']*np.arctan(w/u))
            outputs['q_dot'] = 1/vehicle_mass*(inertia['yy']*1/2*rho*np.sqrt(u**2+w**2)*vehicle_chord*vehicle_planform_area \
            * (8.1576*10**(-8)*(w/u)**4 - 3.3297*10**(-6)*(w/u)**3 - 9.7364*10**(-6)*(w/u)**2 - 0.00051529*(w/u) + 0.00016104) \
            + 1/2*rho*np.sqrt(u**2+w**2)*vehicle_planform_area*((position_aerodynamic_center['x']-position_center_of_gravity['x'])*force_coefficients['L'] \
            - (position_aerodynamic_center['z']-position_center_of_gravity['z'])*force_coefficients['D']) * (position_center_of_thrust['z']-position_center_of_gravity['z']) * thrust['x'])

            # v_dot = g*sin(phi)*cos(theta)+1/vehicle_mass*1/2*rho*sqrt(u**2+v**2+w**2)*vehicle_planform_area*(-force_coefficients.N*asin(v/sqrt(u**2+v**2+w**2)))
            # fi_dot = p + q*sin(phi)*tan(theta) + r*cos(phi)*tan(theta)
            # theta_dot = q*cos(phi) - r*sin(phi)
            # psi_dot = q*sin(phi)*sec(theta) + r*cos(phi)*sec(theta)
            # p_dot = inertia.xx*1/vehicle_mass* ( 1/2*rho*sqrt(u**2+v**2+w**2)*vehicle_chord*vehicle_planform_area*moment_coefficients.l*atan(w/u) + 1/2*rho*sqrt(u**2+v**2+w**2)*vehicle_planform_area*((position_aerodynamic_center.z-position_center_of_gravity.z)*force_coefficients.N*asin(v/sqrt(u**2+v**2+w**2)) - (position_aerodynamic_center.y-position_center_of_gravity.y)*force_coefficients.L*atan(w/u))) + inertia.xz*1/vehicle_mass* ( 1/2*rho*sqrt(u**2+v**2+w**2)*vehicle_chord*vehicle_planform_area*moment_coefficients.n*atan(w/u) + 1/2*rho*sqrt(u**2+v**2+w**2)*vehicle_planform_area*((position_aerodynamic_center.y-position_center_of_gravity.y)*force_coefficients.D*atan(w/u) - (position_aerodynamic_center.x-position_center_of_gravity.x)*force_coefficients.N*asin(v/sqrt(u**2+v**2+w**2))) - (position_center_of_thrust.y-position_center_of_gravity.y)*thrust.x)
            # q_dot = inertia.yy*1/vehicle_mass* ( 1/2*rho*sqrt(u**2+v**2+w**2)*vehicle_chord*vehicle_planform_area*moment_coefficients.m*asin(v/sqrt(u**2+v**2+w**2)) + 1/2*rho*sqrt(u**2+v**2+w**2)*vehicle_planform_area*(-(position_aerodynamic_center.z-position_center_of_gravity.z)*force_coefficients.D*atan(w/u) + (position_aerodynamic_center.x-position_center_of_gravity.x)*force_coefficients.L*atan(w/u)) + (position_center_of_thrust.z-position_center_of_gravity.z)*thrust.x)
            # r_dot = inertia.xz*1/vehicle_mass* ( 1/2*rho*sqrt(u**2+v**2+w**2)*vehicle_chord*vehicle_planform_area*moment_coefficients.l*atan(w/u) + 1/2*rho*sqrt(u**2+v**2+w**2)*vehicle_planform_area*((position_aerodynamic_center.z-position_center_of_gravity.z)*force_coefficients.N*asin(v/sqrt(u**2+v**2+w**2)) - (position_aerodynamic_center.y-position_center_of_gravity.y)*force_coefficients.L*atan(w/u))) + inertia.zz*1/vehicle_mass* ( 1/2*rho*sqrt(u**2+v**2+w**2)*vehicle_chord*vehicle_planform_area*moment_coefficients.n*atan(w/u) + 1/2*rho*sqrt(u**2+v**2+w**2)*vehicle_planform_area*((position_aerodynamic_center.y-position_center_of_gravity.y)*force_coefficients.D*atan(w/u) - (position_aerodynamic_center.x-position_center_of_gravity.x)*force_coefficients.N*asin(v/sqrt(u**2+v**2+w**2))) - (position_center_of_thrust.y-position_center_of_gravity.y)*thrust.x)

            # y[0] = theta_dot
            # y[1] = u_dot
            # y[2] = w_dot
            # y[3] = q_dot
            #objective = np.sum(outputs['theta_dot'] ** 2 + outputs['u_dot'] ** 2 + outputs['w_dot'] ** 2 + outputs['q_dot'] ** 2)
            #self.set_objective('objective', scaler=1.0, adder=0.0)
            #outputs['objective'] = objective
            outputs['objective'] = np.sum(outputs['theta_dot']**2 + outputs['u_dot']**2 + outputs['w_dot']**2 + outputs['q_dot']**2)

def runSolver():
    prob = om.Problem()
    prob.model = om.Group()

    #prob.model.add_subsystem('inputs', om.IndepVarComp('x', shape=4, val=[0, velocity['u'], velocity['w'], 0]), promotes=['*'])
    prob.model.add_subsystem('ode', DiffEqs(), promotes_inputs= ['theta', 'u', 'w', 'q'], promotes_outputs=['theta_dot', 'u_dot', 'w_dot', 'q_dot'])
    '''prob.model.add_subsystem('objective', ObjectiveComponent(), promotes=['obj'])

    prob.model.connect('theta_dot', 'objective.theta_dot')
    prob.model.connect('u_dot', 'objective.u_dot')
    prob.model.connect('w_dot', 'objective.w_dot')
    prob.model.connect('q_dot', 'objective.q_dot')
    '''
    prob.driver = om.ScipyOptimizeDriver()
    prob.driver.options['optimizer'] = 'SLSQP'

    prob.model.add_design_var('u', lower = 2500, upper = 3500)
    prob.model.add_design_var('w')
    prob.model.add_design_var('q')


    # Set the initial values for the variables

    
    prob.model.add_objective('ode.objective')
    # Run the optimization
    prob.setup()
    prob.set_val('theta', 2*np.pi/ 180)
    prob.set_val('u', velocity['u'])
    prob.set_val('w', velocity['w'])
    prob.set_val('q', 0)
    #prob.run_model()
    prob.run_driver()

    # Print the results
    prob.model.list_inputs()
    prob.model.list_outputs()
    return prob.model
    # Print the solution
    
class ObjectiveComponent(om.ExplicitComponent):
    def setup(self):
        self.add_input('theta_dot')
        self.add_input('u_dot')
        self.add_input('w_dot')
        self.add_input('q_dot')
        self.add_output('obj')

    def compute(self, inputs, outputs):
        theta_dot = inputs['theta_dot']
        u_dot = inputs['u_dot']
        w_dot = inputs['w_dot']
        q_dot = inputs['q_dot']

        outputs['obj'] = np.sum(theta_dot**2 + u_dot**2 + w_dot**2 + q_dot**2)
# This part of the code is for construction of the dynamics matrix, you (Meghna) said that
# the matrix components are inputed seperatley to a function that calculates the eigenvalues.
# I have therefore listed all the components here:

# ---Necessary parameters (hard brackets denote parameter [name])---
# - [gamma] - ratio of specific heats
# - [mach] - free stream Mach number
# - [pressure_freestream] - free stream static pressure
# - [Cpn] - pressure coefficient
# - [dP_dH] - free stream pressure derived w.r.t altitude
# - [theta_L] - nominal-flow deflection angle (from Propulsion WP): equal to first panel deflection!(?) or alpha
# - [h] - vehicle height
# - [alpha_0] - trim alpha
# - [delta_0] - trim elevator deflection
# - [S_es] - combined elevator area
# - [b] - vehicle width parameter
# - [rho] - free stream density
# - [V_infty] - free stream velocity
# - [z_bar] - z distance to vehicle mass center from z_min of vehicle
# - [x_bar] - x distance to vehicle mass center from x_max of vehicle
# - [l_1] - inlet ramp length
# - [L_1] - inlet ramp length projected on x
# - [l_2] - outlet ramp length
# - [L_2] - outlet ramp length projected on x
# - [dTh_dM] - thrust derived w.r.t Mach number
# - [dTh_dtheta_L] - thrust derived w.r.t nominal-flow deflection angle
# - [dPe_dM] - internal nozzle exit pressure derived w.r.t to Mach number
# - [dPe_dtheta_L] - internal nozzle exit pressure derived w.r.t to nominal-flow deflection angle
# - [P_bar] - pressure ratio, Pe/pressure_freeastram
# - [pressure_freestram] - free stream static pressure
# - [Pe] - pressure internal nozzle/engine exit conditions
# - [I_1] - function of P_bar
# - [I_2] - function of P_bar

# ---A bit further notations mainly for me (Marcus) to keep up---
# Cpn: pressure coefficient used in stability derivateves calculations, always equal to 2 §Chavez
# theta_L_ equal to first panel deflection!(?) or alpha
# h refers to vehicle height and H to flight altitude
# b: is used i §Chavez as the vehicle width parameter, all moments, forces, intertias, and masses are normalized w it
# V_infty: u_0 (u vel projected on x_vel) or sqrt(u^2+w^2)
# x_bar and z_bar, need to define origin of aircraft
# l_1 and L_1, assumed to be one panel for aircraft instead of two which it actually is
# dTh_ thrust
# dP: pressure at internal nozzle exit
# I_1 = [(P_bar-1)-ln(P_bar)]/(P_bar-1)^2
# I_2 = [(P_bar+1)*ln(P_bar)-2*(P_bar-1)]/(P_bar-1)^3

# example values used for writing of the stability derivative equations (defining variables)
u0 = 2400
w0 = -50
gamma = 1.4
M = 8
pressure_freestram = 1090   #[Pa]
Cpn = 2
Dp_dH = -0.1765             #[pa/m]
theta_L = 0.2618            #[rad] = 15deg
h = 15                      #[m]
alpha_0 = -0.0349           #[rad] = -2deg
delta_0 = 0                 #[rad] = 0deg
S_es = 20                   #[m^2] ~0.6% of planform area
x_es = -39                  #[m] elevator x position w.r.t cog
z_es = 0                    #[m] elevator z position w.r.t cog
b = 60                      #[m]
rho = 0.018                 #[kg/m^3] at 30000 meters
V_infty = 2414.2            #[m/s] at 30000 meters, V_infty = M*sqrt(gamma*R*T)
z_bar = h/2                 #[m]---negative?
x_bar = 70                  #[m]---negative?
l_1 = 64                    #[m] should be: l_1 = (aircraft length)*(inlet fraction)/cos(arctan(h/L*(inlet fraction))) where L is length of aircraft
L_1 = 60                    #[m] should be: L_1 = (aircraft length)*(inlet fraction)
l_2 = 50                    #[m] should be: l_2 = (aircraft length)*(outlet fraction)/cos(arctan(h/L*(outlet fraction)))
L_2 = 40                    #[m] should be: L_1 = (aircraft length)*(inlet fraction)
# The four deriivatives below can be found by calling WP Propulsion and
# evaluating an approx derivative---
dTh_dM = 50000              # determine from WP Prop---
dTh_dtheta_L = -100         # determine from WP Prop---
dPe_dM = 0.5                # determine from WP Prop---
dPe_dtheta_L = 0.05         # determine from WP Prop---
P_bar = 290                 # (combustion chamber exit pressure)/(free stream pressure)
Pe = P_bar*pressure_freestram
dPfree_dH = -0.15
I_1 = ((P_bar-1)-log(P_bar))/(P_bar-1)**2
I_2 = ((P_bar+1)*log(P_bar)-2*(P_bar-1))/(P_bar-1)**3

q = 0.5*rho*(sqrt(u0**2+w0**2))**2
V_inf = 0.5*rho*(sqrt(u0**2+w0**2))**2

# External nozzle stability derivatives
#X_E_Minf = h*I_1*dPe_dM
#X_E_alpha = h*I_1*dPe_dtheta_L
#X_E_q = -h/V_inf*((h-z_bar)*sin(alpha_0)-(L_1-x_bar)*cos(alpha_0))*I_1*dPe_dtheta_L

#Z_E_Minf = -L_2*I_1*dPe_dM
#Z_E_alpha = -L_2*I_1*dPe_dtheta_L
#Z_E_q = L_2/V_inf*((h-z_bar)*sin(alpha_0)-(L_1-x_bar)*cos(alpha_0))*I_1*dPe_dtheta_L

#M_E_Minf = (((h-z_bar)*h-(L_1-x_bar)*L_2)*I_1-l_2**2*I_2)*dPe_dM
#M_E_alpha = (((h-z_bar)*h-(L_1-x_bar)*L_2)*I_1-l_2**2*I_2)*dPe_dtheta_L
#M_E_q = -1/V_inf*((h-z_bar)*sin(alpha_0)-(L_1-x_bar)*cos(alpha_0))* \
#(((h-z_bar)*h-(L_1-x_bar)*L_2)*I_1-l_2**2*I_2)*dPe_dtheta_L

# Aerodynamic stability derivatives
X_A_Minf = -gamma*M*pressure_freestram*Cpn*(sin(theta_L)**2*h+sin(alpha_0+delta_0)**2*sin(delta_0)*S_es/b)

X_A_alpha = -q*Cpn*(sin(2*theta_L)*h+sin(2*(alpha_0+delta_0))*sin(delta_0)*S_es/b)

X_A_q = q*Cpn*sin(2*theta_L)*h/V_inf*((h-z_bar)*sin(alpha_0)- \
(L_1-x_bar)*cos(alpha_0) + 0.5*l_1*cos(theta_L))

Z_A_Minf = -gamma*M*pressure_freestram*Cpn*(sin(theta_L)**2*L_1+sin(alpha_0+delta_0)**2*cos(delta_0)*S_es/b)

Z_A_alpha = -q*Cpn*(sin(2*theta_L)*L_1+sin(2*(alpha_0+delta_0))*cos(delta_0)*S_es/b)

Z_A_q = q*Cpn*sin(2*theta_L)*L_1/V_inf*((h-z_bar)*sin(alpha_0)- \
(L_1-x_bar)*cos(alpha_0) + 0.5*l_1*cos(theta_L))

M_A_Minf = pressure_freestram*gamma*M*Cpn*( sin(theta_L)**2*(0.5*l_1**2-(L_1-x_bar)*L_1-(h-z_bar)*h)+ \
sin(alpha_0+delta_0)**2*(x_es*cos(delta_0)-z_es*sin(delta_0))*S_es/b )

M_A_alpha = q*Cpn*( sin(2*theta_L)*(0.5*l_1**2-(L_1-x_bar)*L_1-(h-z_bar)*h)+ \
sin(2*(alpha_0+delta_0))*(x_es*cos(delta_0)-z_es*sin(delta_0))*S_es/b )

M_A_q = -0.5*q*Cpn*sin(2*theta_L)*( l_1**2/V_inf*(2/3*l_1*cos(theta_L)- \
(L_1-x_bar)*cos(alpha_0)+(h-z_bar)*sin(alpha_0))-z_bar/V_inf*(L_1*(L_1-x_bar)+ \
h*(h-z_bar))*(0.5*l_1*cos(theta_L)-(L-x_bar)*cos(alpha_0)+(h-z_bar)*sin(alpha_0)) )

# Engine-thrust stability derivatives
X_T_Minf = dTh_dM
X_T_alpha = dTh_dtheta_L
X_T_q = -1/V_inf*((h*z_bar)*sin(alpha_0)-(L_1-x_bar)*cos(alpha_0))*dTh_dtheta_L

Z_T_Minf = 0
Z_T_alpha = 0
Z_T_q = 0

M_T_Minf = (h-z_bar)*dTh_dM
M_T_alpha = (h-z_bar)*dTh_dtheta_L
M_T_q = -1/V_inf*((h-z_bar)*sin(alpha_0)-(L_1-x_bar)*cos(alpha_0))*(h-z_bar)*dTh_dtheta_L

# correction coefficient for dynamics matrix
c1 = M/pressure_freestram*dPfree_dH

# Components to dynamics matrix A
a11 = c1 * (X_A_Minf + X_T_Minf) / (alpha_trim * vehicle_mass)
a21 = c1 * (Z_A_Minf) / (alpha_trim * vehicle_mass)
a31 = 0
a41 = c1 * (M_A_Minf + M_T_Minf) / (alpha_trim * inertia.yy)

a12 = c1 * (X_A_alpha + X_T_alpha) / (vehicle_mass)
a22 = c1 * (Z_A_alpha) / (vehicle_mass)
a32 = 0
a42 = c1 * (M_A_alpha + M_T_alpha) / (inertia.yy)

a13 = -g
a23 = 0
a33 = 0
a43 = 0

a14 = c1 * (X_A_q + X_T_q) / (vehicle_mass)
a24 = c1 * (Z_A_q) / (vehicle_mass)
a34 = 1
a44 = c1 * (M_A_q + M_T_q) / (inertia.yy)