"""Python wrapper for CCSM3_SIR_DE"""
import numpy as np

import default_input

class SeaIceRT():

    def __init__(self):
        self.day_of_year = default_input.day_of_year
        self.latitude = default_input.latitude
        self.level = default_input.level
        self.pressure = default_input.pressure
        self.air_temperature = default_input.air_temperature
        self.water_vapor_mixing_ratio = default_input.water_vapor_mixing_ratio
        self.ozone_mixing_ratio = default_input.ozone_mixing_ratio
        self.cloud_cover = default_input.cloud_cover
        self.cloud_liquid_water_path = default_input.cloud_liquid_water_path
        self.surface_pressure = default_input.surface_pressure
        self.co2_volume_mixing_ratio = default_input.co2_volume_mixing_ratio
        self.surface_air_temperature = default_input.surface_air_temperature
        self.ground_temperature = default_input.ground_temperature
        self.snow_depth = default_input.snow_depth
        self.snow_density = default_input.snow_density
        self.snow_grain_radius = default_input.snow_grain_radius
        self.pond_depth = default_input.pond_depth
        self.pond_tuning_parameter = default_input.pond_tuning_parameter
        self.sea_ice_thickness = default_input.sea_ice_thickness
        self.sea_ice_tuning_parameter = default_input.sea_ice_tuning_parameter


    def __repr__(self):
        print(f"{self.day_of_year:7.3f} day of year ",
              "(1=january 1; from 1 to 365) ",
              "(140.477 > mu0 = 0.5) at 80 lat")
        print(f"{self.latitude:5.2f}   latitude    (from +90 to -90)")
        
        
