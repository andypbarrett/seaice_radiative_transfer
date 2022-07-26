"""Test routines for CRM"""

import numpy as np

from interface import input_common, test
import default_input


def to_array(obj, dtype):
    """convert list to array

    :obj: list object
    :dtype: ctype equivalent numpy data type

    :returns: ctype array
    """
    return np.ctypeslib.as_ctypes(np.array(obj).astype(dtype))


def print_parameters(mycom):
    print(f"Day of year (1..365): {mycom.dayyr_in}")
    print(f"Latitude (-90..+90): {mycom.rlat_in}")

    zipped = zip(
        mycom.lev_in[:],
        mycom.pmidm1_in[:],
        mycom.tm1_in[:],
        mycom.qm1_in[:],
        mycom.o3mmr_in[:],
        mycom.cld_in[:],
        mycom.clwp_in[:],
        )
    """level p(mb)  t(k) h2ommr(g/g) o3mmr(g/g) cld cvr  cld lwp (g/m2)
    18    2.   273.   4.0E-06   7.0E-06    0.0E+00  0.0E+00
"""
    print("level p(mb)  t(k) h2ommr(g/g) o3mmr(g/g) cld cvr  cld lwp (g/m2)")
    for l, p, t, q, o, c, w in zipped:
        print(f"   {l:2d} {p:4.0f} {t:3.0f} {q:3.1E} {o:3.1E} {c:3.1E} {w:3.1E}")
    print("")
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


def print_output(output_common):
    print(f"Albedo shortwave direct: {output_common.asdir[0]}")
    print(f"Albedo shortwave diffuse: {output_common.asdif[0]}")
    print(f"Albedo longwave direct: {output_common.aldir[0]}")
    print(f"Albedo longwave diffuse: {output_common.aldif[0]}")
    print(f"Visible solar absorbed by ocean: {output_common.F_SW_ocn_vs}")
    print(f"Near-IR absorbed by ocean: {output_common.F_SW_ocn_ni}")


# Initialize parameters from default values
input_common.dayyr_in = default_input.day_of_year
input_common.rlat_in = default_input.latitude
input_common.ps_in = default_input.surface_pressure
input_common.co2mix_in = default_input.co2_volume_mixing_ratio
input_common.ts_in = default_input.surface_air_temperature
input_common.tg_in = default_input.ground_temperature
input_common.sndpth_in = default_input.snow_depth
input_common.rhos_in = default_input.snow_density
input_common.rs_in = default_input.snow_grain_radius
input_common.hpnd_in = default_input.pond_depth
input_common.R_pnd_in = default_input.pond_tuning_parameter
input_common.hice_in = default_input.sea_ice_thickness
input_common.R_ice_in = default_input.sea_ice_tuning_parameter

input_common.lev_in = to_array(default_input.level, np.int32)
input_common.pmidm1_in = to_array(default_input.pressure, np.float32)
input_common.tm1_in = to_array(default_input.air_temperature, np.float32)
input_common.qm1_in = to_array(default_input.water_vapor_mixing_ratio,
                               np.float32)
input_common.o3mmr_in = to_array(default_input.ozone_mixing_ratio, np.float32)
input_common.cld_in = to_array(default_input.cloud_cover, np.float32)
input_common.clwp_in = to_array(default_input.cloud_liquid_water_path,
                                np.float32)

test()
print_output(output_common)
