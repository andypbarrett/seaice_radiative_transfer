"""Test routines for CRM"""
from itertools import chain

import numpy as np

from interface import input_common, output_common, seaice_common, test
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


def print_output():
    print("-"*70)
    print("CCSM3 Sea Ice Delta Eddington calculation")
    print("-"*70)
    print("-"*70)    
    print("Visible and near-ir direct and diffuse albedos")
    print("   Visible: 0.2 to 0.7 micrometers")
    print("   Near-IR: 0.7 to 5.0 micrometers")
    print("-"*70)    
    print(f"Albedo shortwave direct: {output_common.asdir[0]:4.2f}")
    print(f"Albedo shortwave diffuse: {output_common.asdif[0]:4.2f}")
    print(f"Albedo longwave direct: {output_common.aldir[0]:4.2f}")
    print(f"Albedo longwave diffuse: {output_common.aldif[0]:4.2f}")
    print(" ")
    print("-"*70)
    print("Surface ansorption and Albedos")
    print("-"*70)
    print(f"Visible solar absorbed by ocean: {output_common.F_SW_ocn_vs}")
    print(f"Near-IR absorbed by ocean: {output_common.F_SW_ocn_ni}")
    print('-'*70)
    print('Surface absorption ad albedos')
    print('-'*70)
    print(f"Solar vs direct surface irradiance: {output_common.sols[0]:6.2f} Wm-2")
    print(" ")
    print("-"*70)
    print("Snow/Sea ice transmitted flux (Tr fraction) and absorption (Q Wm-2)")
    print("-"*70)
    print(f"{' '*2} {'Level':10s} {'depth':5s} {'Tr_vs':6s} {'Q_vs':6s} {'Tr_ni':6s} {'Q_ni':6s} {'Q_total':7s}")      
    print("-"*70)
    # Make flux table strings
    zipped = zip(output_common.layer_type[:],
                 output_common.Q_SW_vs_out[:],
                 output_common.Q_SW_ni_out[:],
                 output_common.Q_SW_total_out[:])
    flux_string = []
    for i, (t, qvs, qni, qtt) in enumerate(zipped):
        flux_string.append(
            f"{i:2d} {t[:].decode():10s} {' '*12} {qvs:6.2f} {' '*6} {qni:6.2f} {qtt:6.2f}"
            )
    # Make trasnmission string
    zipped = zip(
        seaice_common.zd[:],
        seaice_common.Tri_vs[:],
        seaice_common.Tri_ni[:]
    )
    trans_string = []
    for zd, tri_vs, tri_ni in zipped:
        trans_string.append(
            f"{' '*13} {zd:5.3f} {tri_vs:6.4f} {' '*6} {tri_ni:6.4f}"
        )
    trans_string.append(None)
    
    for fs, ts in list(zip(flux_string, trans_string)):
        print(fs)
        print(ts)
#    print(f"Up vs flux direct: {output_common.Fdirup_vs[:][0]}")


def print_seaice():
    
    print("-"*70)
    print("CCSM3 Sea Ice Transmission")
    print("-"*70)
    print(f"I_vs: {seaice_common.I_vs}")
    print(f"I_ni: {seaice_common.I_ni}")
    print("level   depth Tr_vs  Q_vs  Tr_ni  Q_ni  Q_total")
    for i, (zd, tri_vs, tri_ni) in enumerate(zipped):
        print(f"{i:2d} {zd:5.3f} {tri_vs:6.4f}       {tri_ni:6.4f}")
    print(f"Tro_vs: {seaice_common.Tro_vs}")
    print(f"Tro_ni: {seaice_common.Tro_ni}")


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
print_output()
#print_seaice()
