"""
    FAST - Copyrigvt (c) 2016 ONERA ISAE
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

import math

import numpy as np
from openmdao.core.explicitcomponent import ExplicitComponent


class Cd0VerticalTail(ExplicitComponent):

    def initialize(self):
        self.options.declare('low_speed_aero', default=False, types=bool)

    def setup(self):
        self.low_speed_aero = self.options['low_speed_aero']

        self.add_input('geometry:vertical_tail:length', val=np.nan, units='m')
        self.add_input('geometry:vertical_tail:thickness_ratio', val=np.nan)
        self.add_input('geometry:vertical_tail:sweep_25', val=np.nan, units='deg')
        self.add_input('geometry:vertical_tail:wetted_area', val=np.nan, units='m**2')
        self.add_input('geometry:wing:area', val=np.nan, units='m**2')
        if self.low_speed_aero:
            self.add_input('reynolds_low_speed', val=np.nan)
            self.add_input('Mach_low_speed', val=np.nan)
            self.add_output('cd0_vt_low_speed')
        else:
            self.add_input('reynolds_high_speed', val=np.nan)
            self.add_input('TLAR:cruise_mach', val=np.nan)
            self.add_output('cd0_vt_high_speed')

    def compute(self, inputs, outputs):
        el_vt = inputs['geometry:vertical_tail:thickness_ratio']
        vt_length = inputs['geometry:vertical_tail:length']
        sweep_25_vt = inputs['geometry:vertical_tail:sweep_25']
        wet_area_vt = inputs['geometry:vertical_tail:wetted_area']
        wing_area = inputs['geometry:wing:area']
        if self.low_speed_aero:
            mach = inputs['Mach_low_speed']
            re_hs = inputs['reynolds_low_speed']
        else:
            mach = inputs['TLAR:cruise_mach']
            re_hs = inputs['reynolds_high_speed']

        ki_arrow_cd0 = 0.04

        cf_vt_hs = 0.455 / (
                (1 + 0.144 * mach ** 2) ** 0.65 * (math.log10(re_hs * vt_length)) ** 2.58)
        ke_cd0_vt = 4.688 * el_vt ** 2 + 3.146 * el_vt
        k_phi_cd0_vt = 1 - 0.000178 * sweep_25_vt ** 2 - 0.0065 * sweep_25_vt
        cd0_vt_hs = (ke_cd0_vt * k_phi_cd0_vt + ki_arrow_cd0 / 8 + 1
                     ) * cf_vt_hs * wet_area_vt / wing_area

        if self.low_speed_aero:
            outputs['cd0_vt_low_speed'] = cd0_vt_hs
        else:
            outputs['cd0_vt_high_speed'] = cd0_vt_hs
