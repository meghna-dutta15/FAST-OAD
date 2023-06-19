from math import log10, sqrt, exp, sin, pi, tan

###################################################################################
## Aerodynamic Calculations
""" 
    @param pressureDict: dictionary with forward location as indices and pressure as value
    @param engineGeometry: acquired from propulsion

"""
def integrateBottomSurfaceEngine(pressureDict, engineGeometry):
    totalLift = 0
    lastLocation = 0
    
    for location, pressure in pressureDict.items():
        if (location != 0):
            totalLift += (location - lastLocation) * engineGeometry.width * pressure

        lastLocation = location

    return totalLift


def calculateLiftNasaHypersonicPaper(mach, AR, dynamicPressure, alpha, area):
    compressibilityBeta = sqrt(mach^2 - 1)
    c1 = 4.17/compressibilityBeta - 0.153 *AR^2* compressibilityBeta
    c2 = exp(0.955- (0.435/mach))

    return c1 * sin(alpha * pi / 180) + c2 * sin(alpha*pi/180)^2 * dynamicPressure *area

class WingGeometry():
    def __init__(self, length):
        self.length = length
        self.sweep = 60
        self.AR = 4 / tan(self.sweep)
        self.width = self.length * self.AR / 2
        self.planform_area = self.width**2 / self.AR 
        self.max_cross_a = 1
        self.max_cross_b = 2
        self.max_cross_area = pi * self.max_cross_a * self.max_cross_b
        self.thickness = self.cross_a / self.cross_b
        self.MAC = 

class Aerodynamics():
    def __init__(self, target_thrust, length, ):
        self.lift  = 0
        self.AR = 1
        self.alpha = 0
        self.

    