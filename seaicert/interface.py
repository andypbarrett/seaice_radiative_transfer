"""
Defines interface with ccsm3_sir_de.for fortran

ToDo: 
- add field descriptions for output common structures
- Debug and fix transmission and absorption arrays
- in test - write print statements to mimic output file tables
- Think about creating pandas dataframes or other structures to hold tables
- Think about work around for class attributes - what is the best number
- make sure wave band definitions are consistent visible, nir

visible (vs) - 0.2 to 0.7 micro-meters
near-ir (ni) - 0.7 to 5.0 micrometers
"""

import ctypes

# These must be the same as in init.f
PLON = 1
PLEV = 18
NXPT = 1
plevp = PLEV + 1
plond = PLON + 1 + 2*NXPT
ksnow = 1
kseaice = 4
klev = ksnow + kseaice + 1
klevp = klev + 1 + 1  # required because defs are 0:klevp


class OutputCom(ctypes.Structure):
    """Defines c-types for output common black used to return
       results"""
    _fields_ = [("asdir", ctypes.c_float * plond),
                ("asdif", ctypes.c_float * plond),
                ("aldir", ctypes.c_float * plond),
                ("aldif", ctypes.c_float * plond),
                ("F_SW_ocn_vs", ctypes.c_float),  # Solar vs absorbed in ocean
                ("F_SW_ocn_ni", ctypes.c_float),  # Solar no absorbed in ocean
                ("layer_type", ctypes.c_char * 10 * klevp),  # Name of layer
                ("Q_SW_vs_out", ctypes.c_float * klevp),  # Solar vs absorbed in layer
                ("Q_SW_ni_out", ctypes.c_float * klevp),  # Solar ni absorbed in layer
                ("Q_SW_total_out", ctypes.c_float * klevp),  # Total solar absorbed in layer
                ("sols", ctypes.c_float * plond),  # Downward solar onto surface sw direct
                ("solsd", ctypes.c_float * plond),  # Downward solar onto surface sw diffuse
                ("vsfdir", ctypes.c_float),  # Visible fraction direct
                ("soll", ctypes.c_float * plond),  # Downward solar on surface ni direct
                ("solld", ctypes.c_float * plond),  # Downward solar on surface ni diffuse
                ("nifdir", ctypes.c_float),  # Near-ir fraction direct
                ("fsds", ctypes.c_float),  # Total solar surface irradiance Wm-2
                ("vsfrac", ctypes.c_float),  # Total visible fraction of irradiance
                ("frs", ctypes.c_float),  # Solar absorbed at surface
                ("albsrf", ctypes.c_float),  # Broadband surface albedo
                ("F_SW_vs", ctypes.c_float),  # Visible absorbed in snow/sea ice Wm-2
                ("F_SW_ni", ctypes.c_float),  # Near-IR absorbed in snow/sea ice Wm-2
                ("F_SW_srf", ctypes.c_float),  # Total absorbed in sea ice surface layer Wm-2
                ]


class RadFluxSeaiceCom(ctypes.Structure):
    """Defines c-types for seaice radiation flux variables"""
    _fields_ = [
        ("hi_ssl", ctypes.c_float),
        ("hs_ssl", ctypes.c_float),
        ("ksrf", ctypes.c_int),
        ("Fdirup_vs", ctypes.c_float * klevp * plond),
        ("Fdirdn_vs", ctypes.c_float * klevp * plond),
        ("Fdifup_vs", ctypes.c_float * klevp * plond),
        ("Fdifdn_vs", ctypes.c_float * klevp * plond),
        ("Fdirup_ni", ctypes.c_float * klevp * plond),
        ("Fdirdn_ni", ctypes.c_float * klevp * plond),
        ("Fdifup_ni", ctypes.c_float * klevp * plond),
        ("Fdifdn_ni", ctypes.c_float * klevp * plond),
        ]


class SeaiceCom(ctypes.Structure):
    """Defines c-types for seaice common block L426"""
    _fields_ = [
        ("I_vs", ctypes.c_float),  # Frac. trns vs thru sea ice sfc
        ("I_ni", ctypes.c_float),  # Frac. trns ni thru sea ice sfc
        ("zd", ctypes.c_float * klevp),  # Interface depths for snow/pond, sea ice
        ("Tri_vs", ctypes.c_float * klevp),  # Frac trns vs sfc to sea ice layers
        ("Tri_ni", ctypes.c_float * klevp),  # Frac trns no sfc to sea ice layers
        ("Tro_vs", ctypes.c_float),  # Frac trns vs to ocean
        ("Tro_ni", ctypes.c_float),  # Frac trns ni to ocean
        ]


class InputCom(ctypes.Structure):
    """Defines c types for input common block that is
       used to set input parameters"""
    _fields_ = [("dayyr_in", ctypes.c_float),
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
                ("R_ice_in", ctypes.c_float),
                ]


# Assign library and common blocks for interface
crmlib = ctypes.CDLL("../1D_dE_CCSM/libcrm.so")  # need generic definition

# Common blocks
input_common = InputCom.in_dll(crmlib, "input_")
output_common = OutputCom.in_dll(crmlib, "output_")
seaice_common = SeaiceCom.in_dll(crmlib, "seaice_")

# Alias for main program function
test = crmlib.crm_

