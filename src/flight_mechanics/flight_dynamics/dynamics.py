import numpy as np
import openmdao as om

# import atmospheric parameters
# import intertia tensor
# import force and moment coefficients
# 

class dynamicsMatrixGeneration:
    def __init__(
            self,
            eulerAngles: float[3]
    ):
        self.eulerAngles = eulerAngles
        self.phi = eulerAngles[0]
        self.theta = eulerAngles[1]
        self.psi = eulerAngles[2]

