import ctypes

# These must be the same as in init.f
PLON = 1
PLEV = 18
NXPT = 1
plevp = PLEV + 1
plond = PLON + 1 + 2*NXPT


def print_parameters(mycom):
    print(f"Day of year (1..365): {mycom.dayyr_in}")
    print(f"Latitude (-90..+90): {mycom.rlat_in}")
    print(f"Levels: {mycom.lev_in[:]}")
    print(f"Pressure at model mid levels: {mycom.pmidm1_in[:]}")
    print(f"Atmospheric temperature: {mycom.tm1_in[:]}")
    print(f"Atmospheric moisture: {mycom.qm1_in[:]}")
    print(f"O3 mass mixing ratio: {mycom.o3mmr_in[:]}")
    print(f"Cloud fraction: {mycom.cld_in[:]}")
    print(f"Cloud liquid water path: {mycom.clwp_in[:]}")
    print(f"Surface pressure: {mycom.ps_in}")
    print(f"CO2 mixing ratio: {mycom.co2mix_in}")
    print(f"Surface temperature: {mycom.ts_in}")
    print(f"Ground temperature: {mycom.tg_in}")
    print(f"Snowdepth: {mycom.sndpth_in}")
    print(f"Snow density: {mycom.rhos_in}")
    print(f"Snow scaling factor: {mycom.rs_in}")
    print(f"Pond depth: {mycom.hpnd_in}")
    print(f"Pond scaling factor: {mycom.R_pnd_in}")
    print(f"Ice thickness: {mycom.hice_in}")
    print(f"Ice scaling factor: {mycom.R_ice_in}")
    return


class Inputcom(ctypes.Structure):
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


crmlib = ctypes.CDLL("../1D_dE_CCSM/libtest.so")

input_com = Inputcom.in_dll(crmlib, "input_")

# Initialize parameters from default values

test = crmlib.test_parameters_
test()

#print_parameters(input_com)

