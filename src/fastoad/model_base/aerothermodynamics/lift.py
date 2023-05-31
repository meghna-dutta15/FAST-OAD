from math import log10, sqrt, exp, sin, pi

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

