"""Tests model run with default parameters"""

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
