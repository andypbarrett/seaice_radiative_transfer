"""Python wrapper for CCSM3_SIR_DE"""
from dataclasses import dataclass
from itertools import chain

import numpy as np

from seaicert.interface import input_common, output_common, seaice_common, test
import seaicert.default_input as default_input


class SeaIceRT():
    """Main class for running CCSM3_SIR_DE

    Model Parameters
    ----------------

    Spatio-temporal parameters:

    :day_of_year: day of year, 1..365, where day 1 = January 1
    :latitude: latitude (-90 to 90)  (test=80.)

    Surface characteristics:
    :surface_pressure:  Surface pressure in mb (test=1008 mb)
    :co2_volume_mixing_ratio:  CO2 volume mixing ratio (test 3.7e-04)
    :surface_air_temperature:  Surface air temperature (K) (test=273.16 K)
    :ground_temperature:  Surface skin temperature (K) (test=273.17 K)
    :snow_depth:  Physical snow depth in meters (test=0 m)
    :snow_density:  Snow density (kg/m3) (test=330 kg/m3)
    :snow_grain_radius:  Snow grain radius in microns (um) (test=50. um)
    :pond_depth:  Physical pond depth in meters (test=0.5 m)
    :pond_tuning_parameter:  Pond tuning parameter in standard deviations (test=-1.)
    :sea_ice_thickness:  Physical ice thickness in meters (test=1.5 m)
    :sea_ice_tuning_parameter:  Sea ice tuning parameter in standard deviations (test=0.)

    Atmospheric Profile - 18 element array-like objects
    :level: number id of level 1..18
    :pressure: Pressure in mb
    :air_temperature: air temperature in Kelvin
    :water_vapor_mixing_ratio:  Water vapour mixing ration (g/g)
    :ozone_mixing_ratio:  Ozone mixing ration (g/g)
    :cloud_cover:  Cloud cover - non-dimension 0.-1.
    :cloud_liquid_water_path: Cloud liquid water path (g/m2)


    Running the model
    In [1]: from ccsm3_sir_de import SeaIceRT

    In [2]: model = SeaIceRT()

    In [3]: model.run()

    # For a single column run, results can be printed to standard out
    In [17]: model.print_results()
    ----------------------------------------------------------------------
    CCSM3 Sea Ice Delta Eddington calculation
    ----------------------------------------------------------------------
    ----------------------------------------------------------------------
    Visible and near-ir direct and diffuse albedos
    Visible: 0.2 to 0.7 micrometers
    Near-IR: 0.7 to 5.0 micrometers
    ----------------------------------------------------------------------
    Albedo shortwave direct: 0.16
    Albedo shortwave diffuse: 0.19
    Albedo longwave direct: 0.04
    Albedo longwave diffuse: 0.06
    
    ...
    snip
    ...
    
    ----------------------------------------------------------------------
    Snow/Sea ice transmitted flux (Tr fraction) and absorption (Q Wm-2)
    ----------------------------------------------------------------------
    Level      depth Tr_vs  Q_vs   Tr_ni  Q_ni   Q_total
    ----------------------------------------------------------------------
    0 surface                  30.67         77.73 108.40
              0.000 1.0000        1.0000
    1 pond                     14.55         76.65  91.20
              0.250 0.9494        0.0130
    2 pond                     13.78          1.06  14.84
              0.500 0.8625        0.0002
    3 ice                       2.34          0.02   2.36
              0.050 0.7889        0.0000
    4 ice                      12.37          0.00  12.37
              0.375 0.6230        0.0000
    5 ice                      10.74          0.00  10.74
              0.750 0.4731        0.0000
    6 ice                       7.74          0.00   7.74
              1.125 0.3552        0.0000
    7 ice                       5.26          0.00   5.26
              1.500 0.2613        0.0000
    8 ocean                    31.46          0.00  31.46

    # Results can be returned as a dictionary
    results = model.get_results()

    # Individual fields can be accessed using normal dictionary access methods
    In [19]: results["downwelling_shortwave_flux_absorbed_by_ocean"]
    Out[19]: 31.46452522277832

    TBD: add example for initialization and running

    Need methods to 
       - repr input parameters
       - get inputs from file/ERA5 or some other reanalysis

    """

    def __init__(self):
        self.set_default_parameters()

    def print_parameters(self):
        for attr, value in self.__dict__.items():
            print(f"{attr} = {value}")

    def run(self):
        """Run sea ice radiative transfer model"""
        set_model_input(self.__dict__)
        test()
#        print_output()
        return

    def set_default_parameters(self):
        """Sets model parameters to default"""
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
        
    def get_results(self):
        """Returns a dictionary of model results"""
        return get_model_results()

    def print_results(self):
        print_output()

    def __repr__(self):
        return (
            "SeaIceRT input parameters:\n" + \
            f"{self.day_of_year:7.3f} day of year " + \
            "(1=january 1; from 1 to 365) " + \
            "(140.477 > mu0 = 0.5) at 80 lat" + "\n" + \
            f"{self.latitude:5.2f}   latitude    (from +90 to -90)" + "\n" + \
            "\n"
        )


def to_array(obj, dtype):
    """convert list to array

    :obj: list object
    :dtype: ctype equivalent numpy data type

    :returns: ctype array
    """
    return np.ctypeslib.as_ctypes(np.array(obj).astype(dtype))


def set_model_input(param_dict):
    """Assigns model parameters to input_common structure

    :param_dict: __dict__ object containing instance attributes of SeaIceRT
    """
    if param_dict["snow_depth"] < 0.:
        raise ValueError("snow_depth must be greater than or equal to zero!")
    if param_dict["pond_depth"] < 0.:
        raise ValueError("pond_depth must be greater than or equal to zero!")
    if param_dict["sea_ice_thickness"] < 0.:
        raise ValueError("sea_ice_thickness must be greater than or equal to zero!")
    if (param_dict["snow_depth"] > 0.) & (param_dict["pond_depth"] > 0.):
        raise ValueError("snow_depth and pond_depth are both greater than zero!\n"
                         "Only one of these can be greater than zero.\n"
                         f"snow_depth: {param_dict['snow_depth']}, pond_depth: {param_dict['pond_depth']}")
    
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

   
def get_model_results():
    """Returns a dictionary of results"""
    results = {}
    results["surface_direct_shortwave_reflectance"] = output_common.asdir[0]
    results["surface_diffuse_shortwave_reflectance"] = output_common.asdif[0]
    results["surface_direct_longwave_reflectance"] = output_common.aldir[0]
    results["surface_diffuse_longwave_reflectance"] = output_common.aldif[0]
    results["downwelling_shortwave_flux_absorbed_by_ocean"] = output_common.F_SW_ocn_vs
    results["downwelling_longwave_flux_absorbed_by_ocean"] = output_common.F_SW_ocn_ni
    results["downwelling_shortwave_flux_absorbed_by_seaice_layer"] = output_common.Q_SW_vs_out[:]
    results["downwelling_longwave_flux_absorbed_by_seaice_layer"] = output_common.Q_SW_ni_out[:]
    results["downwelling_radiative_flux_absorbed_by_seaice_layer"] = output_common.Q_SW_total_out[:]
    results["surface_downwelling_direct_shortwave_flux"] = output_common.sols[0]
    results["surface_downwelling_diffuse_shortwave_flux"] = output_common.solsd[0]
    results["fraction_of_direct_shortwave_at_surface"] = output_common.vsfdir
    results["surface_downwelling_direct_longwave_flux"] = output_common.soll[0]
    results["surface_downwelling_diffuse_longwave_flux"] = output_common.solld[0]
    results["fraction_of_direct_longwave_at_surface"] = output_common.nifdir
    results["surface_downwelling_radiative_flux"] = output_common.fsds
    results["fraction_of_downwelling_radiative_flux_as_shortwave"] = output_common.frs
    results["surface_albedo"] = output_common.albsrf
    results["downwelling_shortwave_absorbed_by_seaice"] = output_common.F_SW_vs
    results["downwelling_longwave_absorbed_seaice"] = output_common.F_SW_ni
    results["downwelling_radiative_flux_absorbed_by_seaice_surface_layer"] = output_common.F_SW_srf
    results["seaice_layer_type"] = [t[:].decode().strip() for t in output_common.layer_type[:]]
    results["fraction_of_shortwave_flux_transmitted_through_seaice_surface"] = seaice_common.I_vs
    results["fraction_of_longwave_flux_transmitted_through_seaice_surface"] = seaice_common.I_vs
    results["depth_of_seaice_layer_interface"] = seaice_common.zd[:]
    results["fraction_of_surface_shortwave_flux_transmitted_to_layer"] = seaice_common.Tri_vs[:]
    results["fraction_of_surface_longwave_flux_transmitted_to_layer"] = seaice_common.Tri_ni[:]
    results["fraction_of_surface_shortwave_flux_transmitted_to_ocean"] = seaice_common.Tro_vs
    results["fraction_of_surface_longwave_flux_transmitted_to_ocean"] = seaice_common.Tro_ni
    return results

    
@dataclass
class DEInput:
    """
    Input dataclass for Delta Eddington model

    day_of_year: float  decimal day of year 1.0 to 365.99 1=January 1
    latitude: float latitude of single column
    """
    day_of_year: float
    latitude: float


def print_output():
    print("-"*70)
    print("CCSM3 Sea Ice Delta Eddington calculation")
    print("-"*70)
    print("-"*70)    
    print("Visible and near-ir direct and diffuse albedos")
    print("   Visible: 0.2 to 0.7 micrometers")
    print("   Near-IR: 0.7 to 5.0 micrometers")
    print("-"*70)    
    print(f"Albedo shortwave direct: {output_common.asdir[0]:4.2f}")
    print(f"Albedo shortwave diffuse: {output_common.asdif[0]:4.2f}")
    print(f"Albedo longwave direct: {output_common.aldir[0]:4.2f}")
    print(f"Albedo longwave diffuse: {output_common.aldif[0]:4.2f}")
    print(" ")
    print("-"*70)
    print("Surface ansorption and Albedos")
    print("-"*70)
    print(f"Visible solar absorbed by ocean: {output_common.F_SW_ocn_vs}")
    print(f"Near-IR absorbed by ocean: {output_common.F_SW_ocn_ni}")
    print('-'*70)
    print('Surface absorption ad albedos')
    print('-'*70)
    print(f"Solar vs direct surface irradiance: {output_common.sols[0]:6.2f} Wm-2")
    print(" ")
    print("-"*70)
    print("Snow/Sea ice transmitted flux (Tr fraction) and absorption (Q Wm-2)")
    print("-"*70)
    print(f"{' '*2} {'Level':10s} {'depth':5s} {'Tr_vs':6s} {'Q_vs':6s} {'Tr_ni':6s} {'Q_ni':6s} {'Q_total':7s}")      
    print("-"*70)
    # Make flux table strings
    zipped = zip(output_common.layer_type[:],
                 output_common.Q_SW_vs_out[:],
                 output_common.Q_SW_ni_out[:],
                 output_common.Q_SW_total_out[:])
    flux_string = []
    for i, (t, qvs, qni, qtt) in enumerate(zipped):
        flux_string.append(
            f"{i:2d} {t[:].decode():10s} {' '*12} {qvs:6.2f} {' '*6} {qni:6.2f} {qtt:6.2f}"
            )
    # Make trasnmission string
    zipped = zip(
        seaice_common.zd[:],
        seaice_common.Tri_vs[:],
        seaice_common.Tri_ni[:]
    )
    trans_string = []
    for zd, tri_vs, tri_ni in zipped:
        trans_string.append(
            f"{' '*13} {zd:5.3f} {tri_vs:6.4f} {' '*6} {tri_ni:6.4f}"
        )
    trans_string.append(None)
    
    for line in list(chain(*zip(flux_string, trans_string)))[:-1]:
        print(line)

#    print(f"Up vs flux direct: {output_common.Fdirup_vs[:][0]}")
