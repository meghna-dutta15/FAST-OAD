"""
Test module for OpenMDAO versions of RubberEngine
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
import openmdao.api as om

from fastoad.constants import FlightPhase
from fastoad.modules.propulsion.fuel_engine.rubber_engine.openmdao import OMRubberEngine
from tests.testing_utilities import run_system


def test_OMRubberEngine():
    """ Tests ManualRubberEngine component """
    # Same test as in test_rubber_engine.test_compute_flight_points
    engine = OMRubberEngine(flight_point_count=(2, 5))

    machs = [0, 0.3, 0.3, 0.8, 0.8]
    altitudes = [0, 0, 0, 10000, 13000]
    thrust_rates = [0.8, 0.5, 0.5, 0.4, 0.7]
    thrusts = [0.955 * 0.8, 0.389, 0.357, 0.0967, 0.113]
    phases = [FlightPhase.TAKEOFF, FlightPhase.TAKEOFF,
              FlightPhase.CLIMB, FlightPhase.IDLE,
              FlightPhase.CRUISE.value]  # mix FlightPhase with integers
    expected_sfc = [0.993e-5, 1.35e-5, 1.35e-5, 1.84e-5, 1.60e-5]

    ivc = om.IndepVarComp()
    ivc.add_output('bypass_ratio', 5)
    ivc.add_output('overall_pressure_ratio', 30)
    ivc.add_output('turbine_inlet_temperature', 1500, units='K')
    ivc.add_output('mto_thrust', 1, units='N')
    ivc.add_output('maximum_mach', 0.95)
    ivc.add_output('design_altitude', 10000, units='m')

    ivc.add_output('mach', [machs, machs])
    ivc.add_output('altitude', [altitudes, altitudes], units='m')
    ivc.add_output('phase', [phases, phases])
    ivc.add_output('use_thrust_rate', [[True] * 5, [False] * 5])
    ivc.add_output('required_thrust_rate', [thrust_rates, [0] * 5])
    ivc.add_output('required_thrust', [[0] * 5, thrusts])

    problem = run_system(engine, ivc)

    np.testing.assert_allclose(problem['SFC'], [expected_sfc, expected_sfc], rtol=1e-2)
    np.testing.assert_allclose(problem['thrust_rate'], [thrust_rates, thrust_rates], rtol=1e-2)
    np.testing.assert_allclose(problem['thrust'], [thrusts, thrusts], rtol=1e-2)