from fastoad.model_base.aerothermodynamics.aerodynamics import WingGeometry
from fastoad.model_base.scramjet_propulsion.scramjet import scramjet

class WeightDistribution():
    def __init__(self,wing, engine, upperParams):
        self.CoGWing = 0
        self.CoGEngineX = 0
        self.CoGEngineZ = 0

        self.wing = wing
        self.engine = engine

        self.upperLevel = upperParams

    def calculateCoGWing(self):
        self.CoGWing = self.wing.length * 2 / 3

    def calculateCoGEngine(self):
        x_coords = self.engine.x_coordinates
        z_coords = self.engine.y_coordinates

        inletRamp1CoGx = 2 * (x_coords[1] - x_coords[0]) / 3
        inletRamp2CoGx = x_coords[2] - (x_coords[2] - x_coords[1]) /3 * (z_coords[1]+2*z_coords[2]) / (z_coords[1]+z_coords[2])
        combusterCoGx = x_coords[2] + (x_coords[3] - x_coords[2])/2
        nozzleCoGx = x_coords[3] + (x_coords[4]-x_coords[3]) / 3

        inletRamp1Area = (x_coords[1] - x_coords[0])*(z_coords[1]+z_coords[0]) / 2
        inletRamp2Area = (x_coords[2] - x_coords[1])*(z_coords[2]+z_coords[1]) / 2
        combusterArea = (x_coords[3] - x_coords[2])*(z_coords[3]+z_coords[2]) / 2
        nozzleArea = (x_coords[4] - x_coords[3])*(z_coords[4]+z_coords[3]) / 2

        inletRamp1CoGz = (z_coords[1] - z_coords[0])/3
        inletRamp2CoGz = ((z_coords[1]/2)*z_coords[1] * (x_coords[2]-x_coords[1]) + ((z_coords[2] - z_coords[1])/3 + \
            z_coords[1])) / inletRamp2Area
        combusterCoGz = z_coords[2] / 2
        nozzleCoGz = z_coords[3] / 3

        self.engineArea = inletRamp1Area + inletRamp2Area + combusterArea + nozzleArea

        self.CoGEngineX = (inletRamp1CoGx * inletRamp1Area + inletRamp2CoGx * inletRamp2Area + \
            combusterCoGx * combusterArea + nozzleCoGx * nozzleArea) / self.engineArea
        
        self.CoGEngineZ = (inletRamp1CoGz * inletRamp1Area + inletRamp2CoGz * inletRamp2Area + \
            combusterCoGz * combusterArea + nozzleCoGz * nozzleArea) / self.engineArea
        
    def calculatePreFuelCoG(self):
        self.massWing = self.wing.density * self.wing.max_cross_area * self.wing.length ** 2 /3
        self.massEngine = self.engine.density * self.engine.width * self.engineArea

        self.totalCoGX = self.massWing * self.CoGWing + self.massEngine * self.CoGEngineX / (self.massWing + self.massEngine)
        self.totalCoGZ = self.CoGEngineZ/(self.massWing + self.massEngine)

    def postFuelCoG(self):
        self.calculatePreFuelCoG()
        self.massFuel = self.engine.mdot_fuel * self.upperLevel.time

        desiredCoGLocation = self.wing.ac - self.wing.length * 0.1
        # massFuel * locationFuel + massPlane * planeCoG / (massFuel + massPlane)
        (desiredCoGLocation * ((self.massFuel + self.massWing) + self.massEngine) - ((self.massWing + self.massEngine) * self.totalCoGX))/self.massFuel


    def placePassengers(self, pax):
        pass
        


