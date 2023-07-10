"""Main function to run Delta-Eddington radiative transfer model"""
import warnings
import datetime as dt

import numpy as np

from seaicert.ccsm3_sir_de import SeaIceRT


RESULTS_TO_RETURN = [
    'surface_albedo',
    'downwelling_shortwave_flux_absorbed_by_ocean',
    'downwelling_longwave_flux_absorbed_by_ocean',
    'surface_downwelling_radiative_flux',
]


def check_isarray(x):
    """Checks that x is numpy.ndarray.  If not returns array."""
    return np.asarray([x]) if np.ndim(x) == 0 else np.asarray(x)


def deciday(thistime):
    """Returns decimal day of year"""
    seconds_in_day = 86400.
    seconds = thistime.second + (thistime.minute * 60) + (thistime.hour * 3600)
    return seconds / seconds_in_day


def decidayofyear(thisdatetime: dt.datetime):
    """Returns decimal day of year"""
    day_of_year = thisdatetime.timetuple().tm_yday
    return day_of_year + deciday(thisdatetime.time())


def run_model(latitude: float,
              localtime: float,
              ice_thickness: float,
              snow_depth: float,
              skin_temperature: float,
              sea_ice_concentration: float,
              pond_depth: float,
              pond_fraction: float,):
    
    """Runs Delta-Eddington RT model

    Arguments
    ---------
    :latitude: (float) latitude in decimal degrees (scalar or array-like)
    :localtime: (datetime object) local time (scalar or array-like) 
    :ice_thickness: (float) ice thickness in meters (scalar or array-like)
    :snow_depth: (float) snow depth in meters (scalar or array-like)
    :skin_temperature: Skin temperature in degrees C (scalar or array-like)
    :sea_ice_concentration: Sea ice concentration [0-1] (scalar or array-like)
    :pond_depth: (float) pond depth in meters (scalar or array-like)
    :pond_fraction: (float) fraction covered by ponds (scalar or array-like)

    :returns: radiative_flux_absorbed_by_the_ocean,
              par_absorbed_by_the_ocean,
              downwelling_radiative_flux_at_the_surface,
              surface_albedo,
              total_transmittance,
    """

    latitude_a = check_isarray(latitude)
    localtime_a = check_isarray(localtime)
    ice_thickness_a = check_isarray(ice_thickness)
    snow_depth_a = check_isarray(snow_depth)
    skin_temperature_a = check_isarray(skin_temperature)
    sea_ice_concentration_a = check_isarray(sea_ice_concentration)
    pond_depth_a = check_isarray(pond_depth)
    pond_fraction_a = check_isarray(pond_fraction)

    shape = ice_thickness_a.shape
    
    # Check all inputs have the same dimensions: except pond_depth and pond_fraction
    for i, arr in enumerate([latitude_a, localtime_a, ice_thickness_a, snow_depth_a,
                             skin_temperature_a, sea_ice_concentration_a,
                             pond_depth_a, pond_fraction_a]):
        if arr.shape != shape:
            raise ValueError("One or more input arrays have mismatched shaped. "
                             f"Expects {shape}, got {arr.shape} for input {i}")

    zipped = zip(
        latitude_a,
        localtime_a,
        ice_thickness_a,
        snow_depth_a,
        skin_temperature_a,
        sea_ice_concentration_a,
        pond_depth_a,
        pond_fraction_a,
        )
    surface_albedo = []
    surface_flux = []
    ocean_flux = []
    transmittance = []
    for lat, tim, hice, hsnow, tskin, sic, hpond, fpond in zipped:
        alb, sfc_flux, ocn_flux, trns = calculate_flux_and_par(lat, tim, hice, hsnow,
                                                               tskin, sic, hpond, fpond)
        surface_albedo.append(alb)
        surface_flux.append(sfc_flux)
        ocean_flux.append(ocn_flux)
        transmittance.append(trns)

    return (np.array(surface_albedo), np.array(surface_flux),
            np.array(ocean_flux), np.array(transmittance))


def calculate_flux_and_par(
        latitude: float,
        localtime: float,
        ice_thickness: float,
        snow_depth: float,
        skin_temperature: float,
        sea_ice_concentration: float,
        pond_depth: float,
        pond_fraction: float,
):
    """Calculate radiative flux and par at ocean, also returns surface radiative flux,
       surface albedo and transmittance as calculated by Delta-Eddington model"""

    for i, sclr in enumerate([latitude, localtime, ice_thickness,
                              snow_depth, skin_temperature,
                              sea_ice_concentration, pond_depth, pond_fraction]):
        if not np.ndim(sclr) == 0:
            warnings.warn(f"Inputs {i} is not scalar: shape: {sclr.shape} "
                          "This may cause unexpected results", UserWarning)

    model = SeaIceRT()

    model.latitude = latitude
    model.localtime = decidayofyear(localtime)
    model.sea_ice_thickness = ice_thickness
    model.snow_depth = snow_depth
    model.surface_air_temperature = skin_temperature
    model.ground_temperature = skin_temperature
    model.pond_depth = 0.

    model.run()

    result = model.get_results()
    output = process_output(result)

    return output


def extract_results(result):
    """Returns values defined in RESULTS_TO_RETURN"""
    return [result.get(name, np.nan) for name in RESULTS_TO_RETURN]


def process_output(result):
    """Extracts results and derive/calculates total absorbed radiative flux and
       transmittance

    :result: dictionary of output parameters from SeaIceRT

    :returns: surface_albedo, surface_downwelling_radiative_flux, 
              downwelling_radiative_flux_absorbed_by_ocean,
              fraction_of_downwelling_radiative_flux_absorbed_by_ocean
    """
    albedo, sw_ocean, lw_ocean, srfc_flux = extract_results(result)
    total_ocean_flux = lw_ocean + sw_ocean
    transmittance = total_ocean_flux / srfc_flux
    return albedo, srfc_flux, total_ocean_flux, transmittance
