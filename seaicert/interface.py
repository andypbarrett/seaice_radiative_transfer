"""Defines interface with ccsm3_sir_de.for fortran"""
import ctypes

# These must be the same as in init.f
PLON = 1
PLEV = 18
NXPT = 1
plevp = PLEV + 1
plond = PLON + 1 + 2*NXPT


class InputCom(ctypes.Structure):
    """Defines c types for input common block that is
       used to set input parameters"""
    _fields_ = [ ("dayyr_in", ctypes.c_float),
                 ("rlat_in", ctypes.c_float),
                 ("lev_in", ctypes.c_int * PLEV),
                 ("pmidm1_in", ctypes.c_float * PLEV),
                 ("tm1_in", ctypes.c_float * PLEV),
                 ("qm1_in", ctypes.c_float * PLEV),
                 ("o3mmr_in", ctypes.c_float * PLEV),
                 ("cld_in", ctypes.c_float * PLEV),
                 ("clwp_in", ctypes.c_float * PLEV),
                 ("ps_in", ctypes.c_float),
                 ("co2mix_in", ctypes.c_float),
                 ("ts_in", ctypes.c_float),
                 ("tg_in", ctypes.c_float),
                 ("sndpth_in", ctypes.c_float),
                 ("rhos_in", ctypes.c_float),
                 ("rs_in", ctypes.c_float),
                 ("hpnd_in", ctypes.c_float),
                 ("R_pnd_in", ctypes.c_float),
                 ("hice_in", ctypes.c_float),
                 ("R_ice_in", ctypes.c_float),]


class OutputCom(ctypes.Structure):
    """Defines c-types for output common black used to return
       results"""
    pass


# Assign library and common blocks for interface
crmlib = ctypes.CDLL("../1D_dE_CCSM/libtest.so")  # need generic definition

# Common blocks
input_common = InputCom.in_dll(crmlib, "input_")

# Alias for main program function
test = crmlib.test_parameters_

