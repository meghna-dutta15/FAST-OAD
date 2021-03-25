"""
Defines the analysis and plotting functions for postprocessing regarding the mission
"""
#  This file is part of FAST-OAD : A framework for rapid Overall Aircraft Design
#  Copyright (C) 2021  ONERA & ISAE-SUPAERO
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
from typing import Union
import pandas as pd
import ipywidgets as widgets
import plotly
import plotly.graph_objects as go
from IPython.display import display, clear_output

COLS = plotly.colors.DEFAULT_PLOTLY_COLORS

BASE_UNITS = {
    "altitude": "m",
    "true_airspeed": "m/s",
    "equivalent_airspeed": "m/s",
    "range": "m",
    "time": "s",
    "ground_distance": "m",
    "mass": "kg",
    "drag": "N",
    "thrust": "N",
    "sfc": "kg/N",
    "slope_angle": "rad",
    "acceleration": "m/s²",
}


class MissionPostprocessing:
    """
    A class for facilitating the post-processing of mission and trajectories
    """

    def __init__(self):
        # The dataframes containing each mission
        self.missions = {}

        # The figure displayed
        self._fig = None

        # The x selector
        self._x_widget = None

        # The y selector
        self._y_widget = None

    def add_mission(self, mission_data: Union[str, pd.DataFrame], name=None):
        """
        Adds the mission to the mission database (self.missions)
        :param mission_data: path of the mission file or Dataframe containing the mission data
        :param name: name to give to the mission
        """
        if type(mission_data) == str and mission_data.endswith(".csv") and pth.exists(mission_data):
            self.missions[name] = pd.read_csv(mission_data, index_col=0)
        elif type(mission_data) == pd.DataFrame:
            self.missions[name] = mission_data
        else:
            raise TypeError("Unknown type for mission data, please use .csv of DataFrame")

        # Initialize widgets when first mission is added
        if len(self.missions) == 1:
            self._initialize_widgets()

    def _initialize_widgets(self):
        """
        Initializes the widgets for selecting x and y
        """

        key = list(self.missions)[0]
        keys = self.missions[key].keys()

        self._x_widget = widgets.Dropdown(value=keys[0], options=keys)
        self._x_widget.observe(self.display, "value")
        self._y_widget = widgets.Dropdown(value=keys[1], options=keys)
        self._y_widget.observe(self.display, "value")

    def _build_plots(self):
        """
        Add a plot of the mission
        """

        x_name = self._x_widget.value
        y_name = self._y_widget.value

        for name in self.missions:
            if self._fig is None:
                self._fig = go.Figure()

            x = self.missions[name][x_name]
            y = self.missions[name][y_name]

            scatter = go.Scatter(x=x, y=y, mode="lines", name=name)

            self._fig.add_trace(scatter)

            self._fig = go.FigureWidget(self._fig)

        # If the variable has no unit we use "[-]"
        x_unit = BASE_UNITS.get(x_name) if BASE_UNITS.get(x_name) else "-"
        y_unit = BASE_UNITS.get(y_name) if BASE_UNITS.get(y_name) else "-"

        self._fig.update_layout(
            title_text="Mission",
            title_x=0.5,
            xaxis_title=x_name + " [" + x_unit + "]",
            yaxis_title=y_name + " [" + y_unit + "]",
        )

    def display(self, change=None) -> display:
        """
        Display the user interface
        :return the display object
        """
        clear_output(wait=True)
        self._update_plots()
        toolbar = widgets.HBox(
            [widgets.Label(value="x:"), self._x_widget, widgets.Label(value="y:"), self._y_widget]
        )
        ui = widgets.VBox([toolbar, self._fig])
        return display(ui)

    def _update_plots(self):
        """
        Update the plots
        """
        self._fig = None
        self._build_plots()
