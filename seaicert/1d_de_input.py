"""
CCSM3 Sea Ice Radiation with Delta Eddington Input from ccsm3_sir_de_input.dat

Notes on the input:

1) Day of year includes the fraction that determines the solar elevation
angle; for example, if day 15.0 was used, the local time of calculation
would be midnight; if day 15.5, conditions would be local noon.  Day of 
year also includes the change in earth/sun distance. The solar constant 
is set by a data statement in routine 'radcsw'. mu0 is the coszine solar
zenith angle.

2) Latitude runs from +90 to -90 degrees.

3) Total number of levels must be the same as the 'plev = ' parameter
in the radiation code.
 
4) Pressure data is at the same levels as temperature and other fields.
Pressure data is in milli-bars.

5) Temperatures are in K, h2o and o3 are mass mixing ratios, cloud cover
as fraction, and cloud liquid water path in g/m2. 

6) The last pressure is the surface pressure.

7) CO2 is assumed to be well-mixed in the atmospheric column.

8) The surface temperatures are given as the air temperature in contact
with the surface, and the actual skin, or ground temperature of the surface.

9) Snow physical depth is in meters. Snow and pond are mutually exclusive;
if both are non-zero, an error occurs. If snow physical depth is non-zero,
it is assumed to completely overlie sea ice. The number of snow layers is 
set by an internal parameter at run time. Snow layers are evenly spaced.

10) Snow density in kg/m3.

11) Snow grain radius in micro-meters.

12) Pond physical depth in meters. If pond  physical depth is non-zero,
it is assumed to completely overlie sea ice.

13) Pond tuning parameter (see Briegleb and Light, 2007)

14) Sea ice thickness in meters. The number of sea ice layers is determined
by an internal parameter at run time. Sea ice layers are evenly spaced.

15) Sea ice tuning parameter (see Briegleb and Light, 2007)


Briegleb, B. P., and B. Light (2007): A Delta-Eddington Multiple
   Scattering Parameterization for Solar Radiation in the Sea Ice
   Component of the Community Climate System Model, NCAR Technical
   Note  NCAR/TN-472+STR  To be published in: February 2007

Example:
  140.477  day of year (1=january 1; from 1 to 365) (140.477 > mu0 = 0.5) at 80 lat
   80.00   latitude    (from +90 to -90)
  level p(mb)  t(k) h2ommr(g/g) o3mmr(g/g) cld cvr  cld lwp (g/m2)
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
     1 1000.   276.   3.0E-03   4.0E-08    0.0E+00  0.0E+00
                   1008.    surface pressure (mb)
                   3.7e-4   co2 volume mixing ratio
                   273.16   surface air temperature (K)
                   273.16   ground temperature (K)
                    0.000   snow physical depth (m)
                  330.000   snow density (kg/m3)
                   50.000   snow grain radius (microns)
                    0.500   pond physical depth (m)
                   -1.000   pond tuning parameter
                    1.500   sea ice thickness (m)
                    0.000   sea ice tuning parameter

"""

# Default profile data
level =  [18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
pressure = [2.0, 5.0, 15.0, 35.0, 60.0, 105.0, 160.0, 235.0, 320.0, 420.0,
            520.0, 610.0, 710.0, 800.0, 870.0, 930.0, 970.0, 1000.0]
temperature = [273.0, 251.0, 234.0, 226.0, 225.0, 225.0, 225.0, 225.0,
               234.0, 247.0, 257.0, 265.0, 272.0, 277.0, 280.0, 281.0,
               278.0, 276.0]
h20mmr = [4e-06, 4e-06, 4e-06, 4e-06, 4e-06, 4e-06, 6.4e-06, 2.6e-05, 0.00012,
          0.00052, 0.0011, 0.002, 0.0031, 0.0042, 0.0051, 0.0059, 0.004, 0.003]
03mmr = [7e-06, 1.3e-05, 1e-05, 5.5e-06, 4.2e-06, 2.2e-06, 1e-06, 5e-07, 2e-07,
         1.4e-07, 1e-07, 8e-08, 7e-08, 6e-08, 5.5e-08, 5e-08, 4.5e-08, 4e-08]
cloud_cover = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
cloud_liquid_water_path = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 60.0, 0.0, 0.0]


class 1D_DeltaEddington():
    """Class to hold input for model"""

    def __init__(day_of_year=140.477,
                 latitude=80.00,
                 surface_pressure=1008.,
                 co2_volume_mixing_ratio=3.7e-4,
                 surface_air_temperature=273.16,
                 ground_temperature=273.16,
                 snow_physical_depth=0.000,
                 snow_density=330.000,
                 snow_grain_radius=50.000,
                 pond_physical_depth=0.500,
                 pond_tuning_parameter=-1.000,
                 sea_ice_thickness=1.500,
                 sea_ice_tuning_parameter=0.000):

        self.day_of_year = day_of_year
        self.latitude = latitude
        self.surface_pressure = surface_pressure
        self.co2_volume_mixing_ratio = co2_volume_mixing_ratio
        self.surface_air_temperature = surface_air_temperature
        self.ground_temperature = ground_temperature
        self.snow_physical_depth = snow_physical_depth
        self.snow_density = snow_density
        self.snow_grain_radius = snow_grain_radius
        self.pond_physical_depth = pond_physical_depth
        self.pond_tuning_parameter = pond_tuning_parameter
        self.sea_ice_thickness = sea_ice_thickness
        self.sea_ice_tuning_parameter = sea_ice_tuning_parameter

        self.level = level
        self.pressure = pressure
        self.temperature = temperature
        self.h2ommr = h2ommr
        self.o3mmr = o3mmr
        self.cloud_cover = cloud_cover
        self.cloud_liquid_water_path = cloud_liquid_water_path

    def write_input(input_file):
        """Writes input to file used by 1D_DE"""

        print("-----------------------------------------------------------------")
        print("CCSM3 Sea Ice Radiation with Delta Eddington Input")
        print("-----------------------------------------------------------------")
        print(f"  {self.day_of_year:7.3f}  day of year (1=january 1; from 1 to 365) (140.477 > mu0 = 0.5) at 80 lat")
        print(f"  {self.latitude:5.2f}   latitude    (from +90 to -90)")
        print("  level p(mb)  t(k) h2ommr(g/g) o3mmr(g/g) cld cvr  cld lwp (g/m2)")
        for i in len(self.level):
            print(f"    {self.level[i]:2d} {pressure[i]:4.0f}.   {self.temperature[i]:.0f}.   {self.h2ommr[i]:.1E}   {self.o3mmr[i]:.1E}    {self.cloud_cover[i]:.1E}  {self.cloud_liquid_water_path:.1E}")
        print(f"       {self.surface_pressure:.0f}.                surface pressure (mb)")
        print(f"                   {self.co2_volume_mixing_ratio:.1e}   co2 volume mixing ratio")
        print(f"                   {self.surface_air_temperature:.2f}   surface air temperature (K)")
        print(f"                   {self.ground_temperature:.2f}   ground temperature (K)")
        print(f"                    {self.snow_physical_depth:.3f}   snow physical depth (m)")
        print(f"                  {self.snow_density:.3f}   snow density (kg/m3)")
        print(f"                   {self.snow_grain_radius:.3f}   snow grain radius (microns)")
        print(f"                    {self.pond_physical_depth:.3f}   pond physical depth (m)")
        print(f"                   {self.pond_tuning_parameter: .3f}   pond tuning parameter")
        print(f"                    {self.sea_ice_thickness:.3f}   sea ice thickness (m)")
        print(f"                    {self.sea_ice_tuning_parameter:.3f}   sea ice tuning parameter")


