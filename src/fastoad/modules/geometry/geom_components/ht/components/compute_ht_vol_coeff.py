"""
    Estimation of horizontal tail volume coefficient
"""

#  This file is part of FAST : A framework for rapid Overall Aircraft Design
#  Copyright (C) 2019  ONERA/ISAE
#  FAST is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.
import numpy as np
from openmdao.core.explicitcomponent import ExplicitComponent

from fastoad.utils.physics.atmosphere import Atmosphere


class ComputeHTVolCoeff(ExplicitComponent):
    # TODO: Document equations. Cite sources
    """ Horizontal tail volume coefficient estimation """

    def setup(self):

        self.add_input('cg_airframe:A51', val=np.nan, units='m')
        self.add_input('cg_airframe:A52', val=np.nan, units='m')
        self.add_input('weight:MTOW', val=np.nan, units='kg')
        self.add_input('geometry:wing_area', val=np.nan, units='m**2')
        self.add_input('geometry:wing_l0', val=np.nan, units='m')
        self.add_input('cg:required_cg_range', val=np.nan)

        self.add_output('delta_lg')
        self.add_output('geometry:ht_vol_coeff')

        self.declare_partials('delta_lg',
                              ['cg_airframe:A51', 'cg_airframe:A52'],
                              method='fd')
        self.declare_partials('geometry:ht_vol_coeff', '*', method='fd')

    def compute(self, inputs, outputs):
        cg_a51 = inputs['cg_airframe:A51']
        cg_a52 = inputs['cg_airframe:A52']
        mtow = inputs['weight:MTOW']
        wing_area = inputs['geometry:wing_area']
        l0_wing = inputs['geometry:wing_l0']
        required_cg_range = inputs['cg:required_cg_range']

        delta_lg = cg_a51 - cg_a52
        atm = Atmosphere(0.)
        rho = atm.density
        sos = atm.speed_of_sound
        vspeed = sos * 0.2  # assume the corresponding Mach of VR is 0.2

        cm_wheel = 0.08 * delta_lg * mtow * 9.81 / \
            (0.5 * rho * vspeed**2 * wing_area * l0_wing)
        delta_cm = mtow * l0_wing * required_cg_range * \
            9.81 / (0.5 * rho * vspeed**2 * wing_area * l0_wing)
        ht_vol_coeff = cm_wheel + delta_cm
        outputs['delta_lg'] = delta_lg
        outputs['geometry:ht_vol_coeff'] = ht_vol_coeff