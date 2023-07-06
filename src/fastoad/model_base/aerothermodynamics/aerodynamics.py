from math import log10, sqrt, exp, sin, pi, tan, cos, atan
from scipy.constants import R, atmosphere
from fastoad.model_base.atmosphere import Atmosphere
from fastoad.model_base.scramjet_propulsion.shocks import oblique, pm_expansion

###################################################################################
## Aerodynamic Calculations
""" 
    @param pressureDict: dictionary with forward location as indices and pressure as value
    @param engineGeometry: acquired from propulsion

"""



class WingGeometry():
    def __init__(self, length, sweep):
        self.length = length
        self.sweep = sweep * pi / 180
        self.AR = 4 / tan(self.sweep)
        self.width = self.length * self.AR / 2
        self.planformArea = self.width**2 / self.AR 
        self.max_cross_a = 1
        self.max_cross_b = 2
        self.max_cross_area = pi * self.max_cross_a * self.max_cross_b
        self.thickness = self.max_cross_a / self.max_cross_b
        self.cL = {}
        self.breakpoint = 0.8
        self.leadingAngle = atan(self.max_cross_b / (self.length * self.breakpoint))
        self.trailingAngle = atan(self.max_cross_b / (self.length * (1 - self.breakpoint)))
        self.noseEmissivity = 0.8
        self.maxTemp = 4350
        self.noseRadius = 0

    def setNoseRadius(self, radius):
        self.noseRadius = radius
    
    def setCrossSection(self, a, b):
        self.max_cross_a = a
        self.max_cross_b = b
        self.thickness = a/b

    def setBreakpoint(self, breakpoint):
        self.breakpoint = breakpoint

          
         

class Aerodynamics():
    def __init__(self, engine, wing, altitude, mach, AR):
        self.lift  = 0
        self.AR = AR
        self.length = engine.x_coordinates[-1]
        self.wing = wing
        #self.target_thrust = engine.target_thrust
        self.mach = mach
        self.criticalMach = sqrt(1 + (4/self.AR)**2)
        self.beta = (self.mach**2-1) ** 0.5

        self.bluntnessDrag = 0

        # All properties dependent on aoa
        self.wingLift = {}
        self.inducedDrag = {}
        self.wingDrag = {}
        self.forebodyPressure = {}
        self.afterbodyPressure = {}

        # Define atmospheric properties
        self.atmosphere = Atmosphere(altitude)
        self.processAtmosphere()

    
    def calculateEngineLiftDrag(self, engine):
        p1 = engine.ramp1_vals['p'] - self.pressure
        p2 = engine.ramp2_vals['p'] - self.pressure
        p3 = engine.inlet_pressure - self.pressure
        p4 = engine.outlet_pressure - self.pressure
        p5 = engine.exhaust_pressure - self.pressure

        pressures = [p1,p2,p3,p4,p5] 

        self.engineLift = 0
        self.engineDrag = 0
        lastLocation = 0

        coords = engine.x_coordinates

        pressureDict = {coords[0]: p1, coords[1]: p2, coords[2]: p3, coords[3]: p4, coords[4]:p5}

        angles = engine.angles
        count = 0
    
        for location, pressure in pressureDict.items():
            if (location != 0):
                self.engineLift += (location - lastLocation) * pressure * sin(angles[count])
                self.engineDrag += (location - lastLocation) * pressure * cos(angles[count])
                count += 1
            lastLocation = location

        self.thrust = engine.engine_thrust
        return [self.engineLift, self.engineDrag]
    
    def processAtmosphere(self):
        atmosphere = self.atmosphere
        self.temperature = atmosphere.temperature
        self.pressure = atmosphere.pressure
        self.density = atmosphere.density

        self.speedOfSound = sqrt(R * 1.4 * atmosphere.temperature)
        self.velocity = self.mach * self.speedOfSound

    
    def calculateWingLift(self, wing, aoa):

        if (self.beta < 4/wing.AR):
            c1 = pi * wing.AR - 0.153*self.beta*(wing.AR)**2
            interp = self.beta / (4/wing.AR)
        else:
            c1 = 4.17/self.beta - 0.13
            interp = 1
        c2 = interp * exp(0.955 - (4.35 / self.mach))

        wing.cl[aoa] = c1 * sin(aoa) + c2 * sin(aoa)**2
        
        
        self.wingLift[aoa] = wing.cl[aoa] * self.density * self.velocity ** 2 * wing.planformArea / 2

    def calculateInducedDrag(self, wing, aoa):
        self.inducedDrag[aoa] = wing.cl[aoa] * tan(aoa)

    def calculateForebodyPressure(self, wing, aoa):
        # Van Dyke Cone
        a = wing.max_cross_a
        b = wing.max_cross_b
        vanDyke = a * b * (2 * log10(4 / (self.beta * (a + b) ))-1)\
            + self.beta**2 * a**2 * b**2 * (3/2 - (a**2 + b**2) / (a*b) + (3/2 * (a**2 + b**2) / (a*b)-2)*log10(4/(self.beta*(a+b))))
        # Oblique Shock
        M2, beta, p2, r2, T2 = oblique(self.mach, self.pressure, self.density, self.temperature, aoa)
        self.forebodyPressure[aoa] = p2
        areaForebody = (self.length * wing.breakpoint) * (wing.breakpoint * wing.width) / 2
        newton = (p2 * cos(aoa) - self.pressure) / (1/2 * self.density * self.velocity **2 * wing.planformArea)

        newtonFraction = self.mach / (12 - self.criticalMach)
        diffCoeff = newton - vanDyke
        self.forebodyPressure[aoa] = newton/30 #vanDyke + newtonFraction ** 2 * diffCoeff
    
    def calculateAfterbodyPressure(self, wing, aoa):
        turningAngle = wing.leadingAngle + wing.trailingAngle
        M2, p2, r2, T2 = pm_expansion(self.mach, self.pressure, self.density, self.temperature, turningAngle)

        self.afterbodyPressure[aoa] = (- p2 *cos(wing.trailingAngle - (aoa - wing.leadingAngle))- self.pressure) / (1/2 * self.density * self.velocity ** 2)


    def calculateSkinFrictionCoefficientLaminar():
        pass
    
    def calculateSkinFrictionCoefficientTurbulent():
        pass

    def calculateSkinFrictionDrag():
        pass

    def sizeNoseRadius(self, wing):
        wing.setNoseRadius(1820*(self.pressure /atmosphere) ** 0.5 * (self.mach*self.velocity*10 **-4) **3.15 / (wing.noseEmmisivity * (wing.maxTemp / 1000)**4))

    def calculateBodyBluntnessDrag(self, wing):
        self.bluntnessDrag = wing.noseRadius **2 *pi / wing.planformArea
        

    def calculateTotalDrag(self, wing, aoa):
        self.cd = self.inducedDrag[aoa] + self.forebodyPressure[aoa] + self.afterbodyPressure[aoa]/100 + self.bluntnessDrag
        self.wingDrag[aoa] = self.cd * self.density * self.velocity ** 2 * wing.planformArea / 2

    def calculateEngineWidth(self, aoa):
        #totalDrag = self.wingDrag[aoa] + self.engineDrag # * engineWidth /self.thrust
        #self.engineWidth = totalDrag / self.thrust

        self.engineWidth = self.thrust * self.wingDrag[aoa] / (self.thrust + 1)
        self.scaledEngineLift = self.engineLift * self.engineWidth
        self.scaledEngineDrag = self.engineDrag * self.engineWidth
 
        return self.engineWidth
    
    def calculateTotalLift(self, aoa):
        self.totalLift = self.wingLift[aoa] + self.scaledEngineLift
        self.projectedMass = self.totalLift / 9.81

