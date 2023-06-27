"""Default inputs for RT"""

atmospheric_data_table = """  level p(mb)  t(k) h2ommr(g/g) o3mmr(g/g) cld cvr  cld lwp (g/m2)
    18    2.   273.   4.0E-06   7.0E-06    0.0E+00  0.0E+00
    17    5.   251.   4.0E-06   1.3E-05    0.0E+00  0.0E+00
    16   15.   234.   4.0E-06   1.0E-05    0.0E+00  0.0E+00
    15   35.   226.   4.0E-06   5.5E-06    0.0E+00  0.0E+00
    14   60.   225.   4.0E-06   4.2E-06    0.0E+00  0.0E+00
    13  105.   225.   4.0E-06   2.2E-06    0.0E+00  0.0E+00
    12  160.   225.   6.4E-06   1.0E-06    0.0E+00  0.0E+00
    11  235.   225.   2.6E-05   5.0E-07    0.0E+00  0.0E+00
    10  320.   234.   1.2E-04   2.0E-07    0.0E+00  0.0E+00
     9  420.   247.   5.2E-04   1.4E-07    0.0E+00  0.0E+00
     8  520.   257.   1.1E-03   1.0E-07    0.0E+00  0.0E+00
     7  610.   265.   2.0E-03   8.0E-08    0.0E+00  0.0E+00
     6  710.   272.   3.1E-03   7.0E-08    0.0E+00  0.0E+00
     5  800.   277.   4.2E-03   6.0E-08    0.0E+00  0.0E+00
     4  870.   280.   5.1E-03   5.5E-08    0.0E+00  0.0E+00
     3  930.   281.   5.9E-03   5.0E-08    1.0E+00 60.0E+00
     2  970.   278.   4.0E-03   4.5E-08    0.0E+00  0.0E+00
     1 1000.   276.   3.0E-03   4.0E-08    0.0E+00  0.0E+00"""



day_of_year = 140.477  # 1=january1; from 1 to 365 (140.477 > mu0 = 0.5 at 80 lat)
latitude = 80.00

# Load atmospheric data table
level = []
pressure = []
air_temperature = []
water_vapor_mixing_ratio = []
ozone_mixing_ratio = []
cloud_cover = []
cloud_liquid_water_path = []
for line in atmospheric_data_table.split('\n')[1:]:
    lev, prs, ta, wvmr, o3mr, clcv, cldlwp = line.split()
    level.append(float(lev))
    pressure.append(float(prs))
    air_temperature.append(float(ta))
    water_vapor_mixing_ratio.append(float(wvmr))
    ozone_mixing_ratio.append(float(o3mr))
    cloud_cover.append(float(clcv))
    cloud_liquid_water_path.append(float(cldlwp))


surface_pressure = 1008.  # surface pressure (mb)
co2_volume_mixing_ratio = 3.7e-4  # co2 volume mixing ratio
surface_air_temperature = 273.16  # surface air temperature (K)
ground_temperature = 273.16  # ground temperature (K)
snow_depth = 0.000  # snow physical depth (m)
snow_density = 330.000  # snow density (kg/m3)
snow_grain_radius = 180.000  # snow grain radius (microns)
pond_depth = 0.500  # pond physical depth (m)
pond_tuning_parameter = -1.000  # pond tuning parameter
sea_ice_thickness = 1.500  # sea ice thickness (m)
sea_ice_tuning_parameter = 0.000  # sea ice tuning parameter
