"""Test for checking Delta-Eddington code is running correctly

Output stored in 1D_DE_test is compared against data taken from corresponding 
files in 1D_DE.  The output files in 1D_DE came in the code package zipfile.
"""

import pytest

from parse_ccsm3_sir_output import DeltaEdOutput
from filepath import TEST_OUTPUT_DIRPATH


def test_bare_ice_transmission():
    filepath = TEST_OUTPUT_DIRPATH / 'ccsm3_sir_de_output_fig8bare15m.dat'
    bare_ice_1_5m = DeltaEdOutput(filepath)
    print(bare_ice_1_5m.
