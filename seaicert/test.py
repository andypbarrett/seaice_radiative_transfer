"""Test routines for CRM"""

from interface import input_common, test

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

# Initialize parameters from default values
test()

#print_parameters(input_com)

