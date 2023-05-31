import numpy as np

# import atmospheric parameters
# import position paramteters, such as geodetic latitude, longitude, and height

class dynamicsMatrixGeneration:
    def __init__(
            self,
            eulerAngles: float[3]
    ):
        self.eulerAngles = eulerAngles
        self.phi = eulerAngles[0]
        self.theta = eulerAngles[1]
        self.psi = eulerAngles[2]

    @property
    def 