"""Tests for run_model"""

import datetime as dt
import numpy as np

from seaicert.model import run_model

# Define inputs to match default inputs
latitude = 80.0
thisdatetime = dt.datetime(2023, 5, 20, 11, 26, 53, 232000)
ice_thickness = 1.5
snow_depth = 0.0
skin_temperature = 273.16
pond_depth = 0.5
sea_ice_concentration = 1.
pond_fraction = 1.
expected_albedo = 0.14001
expected_surface_flux = 179.19
expected_ocean_flux = 27.57
expected_transmittance = 0.1538

def test_default_run():
    """Tests default inputs"""
    result = run_model(latitude, thisdatetime,
                       ice_thickness, snow_depth,
                       skin_temperature,
                       sea_ice_concentration,
                       pond_depth, pond_fraction)
    expected = np.array([expected_albedo, expected_surface_flux,
                               expected_ocean_flux, expected_transmittance])
    result = np.concatenate(result)
    assert np.allclose(result, expected, atol=1e-2)
