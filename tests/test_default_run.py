"""Tests model run with default parameters"""
import pytest

import numpy as np

from seaicert.ccsm3_sir_de import SeaIceRT

expected = {
    "surface_direct_shortwave_reflectance": 0.17,
    }

def test_default():
    model = SeaIceRT()
    model.run()
    actual = model.get_results()
    assert expected["surface_direct_shortwave_reflectance"] == np.round(actual["surface_direct_shortwave_reflectance"], 2)


# Input values are sea_ice_thickness, snow_depth, pond_depth
@pytest.mark.parametrize(
    "test_input,expected",
    [
        ([1.5, -1., 0.], (ValueError, r"snow_depth must be greater than or equal to zero!")),
        ([1.5, 0., -1.], (ValueError, r"pond_depth must be greater than or equal to zero!")),
        ([-1., 0., 0.], (ValueError, r"sea_ice_thickness must be greater than or equal to zero!")),
        ([1.5, 0.3, 0.3], (ValueError, r"snow_depth and pond_depth are both greater than zero!*"))
    ]
)
def test_valueerror(test_input, expected):
    """Checks that ValueError raised for snow_depth, pond_depth 
    and sea_ice_thickness cases"""
    model = SeaIceRT()
    model.sea_ice_thickness = test_input[0]
    model.snow_depth = test_input[1]
    model.pond_depth = test_input[2]
    with pytest.raises(expected[0], match=expected[1]):
        model.run()
    
