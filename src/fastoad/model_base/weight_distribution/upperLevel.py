from fastoad.model_base.atmosphere import Atmosphere
from scipy.constants import R

class UpperLevel():
    '''
    Distance (km)
    '''
    def __init__(self, pax, distance, mach, altitude, gamma = 1.4):
        self.pax = pax
        self.distance = distance
        self.mach = mach
        self.altitude = altitude
        self.temperature = Atmosphere(altitude).temperature
        self.gamma = gamma
        self.velocity = 0
        self.time = 0
        self.paxMass = 0

    def calculateVelocity(self):
        speedOfSound = (self.gamma * R * self.temperature) ** 0.5
        self.velocity = speedOfSound * self.mach

        return self.velocity
    
    def flightTime(self):
        vel = self.calculateVelocity()
        self.time = self.distance / vel

        return self.time
    
    def calculatePaxMass(self, massPerPax = 95):
        self.paxMass = self.pax * massPerPax
        return self.paxMass
    
    