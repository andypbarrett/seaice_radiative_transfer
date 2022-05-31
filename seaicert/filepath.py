"""Contains directory and file paths"""

from pathlib import Path

HOME = Path.home()
PACKAGE_DIRPATH = HOME / 'src' / 'seaice_radiative_transfer' / '1D_dE_CCSM'
TEST_OUTPUT_DIRPATH = PACKAGE_DIRPATH / '1D_DE_test'
