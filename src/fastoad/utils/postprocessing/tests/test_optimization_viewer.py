"""
Tests for FAST-OAD optimization viewer
"""
#  This file is part of FAST : A framework for rapid Overall Aircraft Design
#  Copyright (C) 2020  ONERA & ISAE-SUPAERO
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

import os.path as pth
import pytest
from shutil import rmtree

from fastoad.io.configuration import FASTOADProblem
from .. import OptimizationViewer
from ..exceptions import FastMissingFile
from fastoad.cmd import api

DATA_FOLDER_PATH = pth.join(pth.dirname(__file__), "data")
RESULTS_FOLDER_PATH = pth.join(pth.dirname(__file__), "results")


@pytest.fixture(scope="module")
def cleanup():
    rmtree(RESULTS_FOLDER_PATH, ignore_errors=True)


def test_optimization_viewer_load(cleanup):
    """
    Basic tests for testing the OptimizationViewer load method.
    """
    filename = pth.join(DATA_FOLDER_PATH, "valid_sellar.toml")

    # The problem has not yet been ran
    problem = FASTOADProblem()
    problem.configure(filename)

    optim_viewer = OptimizationViewer()

    # No input file exists
    with pytest.raises(FastMissingFile):
        optim_viewer.load(problem)

    api.generate_inputs(filename)

    # Input file exist
    optim_viewer.load(problem)

    # We run the problem
    api.optimize_problem(filename, overwrite=True)

    # Load the results
    optim_viewer.load(problem)
