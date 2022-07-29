"""Python wrapper for CCSM3_SIR_DE"""
from dataclasses import dataclass

import numpy as np

from interface import input_common, output_common, seaice_common, test
import default_input


class SeaIceRT():
    """Main class for running CCSM3_SIR_DE

    TBD: add example for initialization and running

    Need methods to 
       - initialize structure
       - initialize model - part of run
       - repr input parameters
       - get results
       - run model
            - copy initialization to model
            - run model
            - copy results to sructure/variables
       - repr results
       - get inputs from file/ERA5 or some other reanalysis

    """

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

    def get_parameters(self):
        for attr, value in self.__dict__.items():
            print(f"{attr} = {value}")

    def run(self):
        """Run sea ice radiative transfer model"""
        set_model_input(self.__dict__)
        return

    def __repr__(self):
        return (
            "SeaIceRT input parameters:\n" + \
            f"{self.day_of_year:7.3f} day of year " + \
            "(1=january 1; from 1 to 365) " + \
            "(140.477 > mu0 = 0.5) at 80 lat" + "\n" + \
            f"{self.latitude:5.2f}   latitude    (from +90 to -90)" + "\n" + \
            "\n"
        )


def set_model_input(param_dict):
    """Assigns model parameters to input_common structure

    :param_dict: __dict__ object containing instance attributes of SeaIceRT
    """
    input_common.dayyr_in = param_dict["day_of_year"]
    input_common.rlat_in = param_dict["latitude"]
    input_common.ps_in = param_dict["surface_pressure"]
    input_common.co2mix_in = param_dict["co2_volume_mixing_ratio"]
    input_common.ts_in = param_dict["surface_air_temperature"]
    input_common.tg_in = param_dict["ground_temperature"]
    input_common.sndpth_in = param_dict["snow_depth"]
    input_common.rhos_in = param_dict["snow_density"]
    input_common.rs_in = param_dict["snow_grain_radius"]
    input_common.hpnd_in = param_dict["pond_depth"]
    input_common.R_pnd_in = param_dict["pond_tuning_parameter"]
    input_common.hice_in = param_dict["sea_ice_thickness"]
    input_common.R_ice_in = param_dict["sea_ice_tuning_parameter"]
    
    input_common.lev_in = to_array(param_dict["level"], np.int32)
    input_common.pmidm1_in = to_array(param_dict["pressure"], np.float32)
    input_common.tm1_in = to_array(param_dict["air_temperature"], np.float32)
    input_common.qm1_in = to_array(param_dict["water_vapor_mixing_ratio"],
                                   np.float32)
    input_common.o3mmr_in = to_array(param_dict["ozone_mixing_ratio"], np.float32)
    input_common.cld_in = to_array(param_dict["cloud_cover"], np.float32)
    input_common.clwp_in = to_array(param_dict["cloud_liquid_water_path"],
                                    np.float32)
    return

   
@dataclass
class DEInput:
    """
    Input dataclass for Delta Eddington model

    day_of_year: float  decimal day of year 1.0 to 365.99 1=January 1
    latitude: float latitude of single column
    """
    day_of_year: float
    latitude: float
    
