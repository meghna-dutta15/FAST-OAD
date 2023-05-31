from math import log10

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


def calculateCoefficientFriction(reynolds):
    return 0.37*log10(reynolds)**2.584