# SeaIceRT: A python interface for CESM3 sea ice radiative transfer code.

SeaIceRT is a python interface for CESM3 Delta Eddington radiative
transfer for sea ice.  The radiative transfer code is written in
Fortran 77.  The python wrapper allows the sea ice parameters to be
set, the code run and output returned.

The `seaicert` directory contains the class SeaIceRT, which
initializes parameters and call the code.

The fortran code is in the `1D_SIR_DE` directory.  The main fortran
routine has been modified to allow the python wrapper to set and get
parameters.  The main fortran calling routine `crm` has been changed
from a Fortan 77 `PROGRAM` to a `SUBROUTINE`.  This allows the fortran
code to be built as a library, which is called from python using the
`ctypes` package.

A full description of the original model can be found [here](https://opensky.ucar.edu/islandora/object/technotes:484)

The easiest way to install the wrapper and model is using git.  

```
git clone git@github.com:andypbarrett/seaice_radiative_transfer.git

cd seaice_radiative_transfer
```

Running the model can be done from a python IDE, either `python` or `ipython`.

```
(base) nsidc-abarrett-442:seaice_radiative_transfer$ cd seaicert/
(base) nsidc-abarrett-442:seaicert$ ipython
Python 3.7.6 | packaged by conda-forge | (default, Mar  5 2020, 15:27:18) 
Type 'copyright', 'credits' or 'license' for more information
IPython 7.17.0 -- An enhanced Interactive Python. Type '?' for help.

In [1]: from ccsm3_sir_de import SeaIceRT

In [2]: model = SeaIceRT()

In [3]: model.run()

In [4]: model.print_results()
----------------------------------------------------------------------
CCSM3 Sea Ice Delta Eddington calculation
----------------------------------------------------------------------
----------------------------------------------------------------------
Visible and near-ir direct and diffuse albedos
   Visible: 0.2 to 0.7 micrometers
   Near-IR: 0.7 to 5.0 micrometers
----------------------------------------------------------------------
Albedo shortwave direct: 0.17
Albedo shortwave diffuse: 0.19
Albedo longwave direct: 0.06
Albedo longwave diffuse: 0.06
 
----------------------------------------------------------------------
Surface ansorption and Albedos
----------------------------------------------------------------------
Visible solar absorbed by ocean: 27.5656681060791
Near-IR absorbed by ocean: 0.0
----------------------------------------------------------------------
Surface absorption ad albedos
----------------------------------------------------------------------
Solar vs direct surface irradiance:   0.12 Wm-2
 
----------------------------------------------------------------------
Snow/Sea ice transmitted flux (Tr fraction) and absorption (Q Wm-2)
----------------------------------------------------------------------
   Level      depth Tr_vs  Q_vs   Tr_ni  Q_ni   Q_total
----------------------------------------------------------------------
 0 surface                  26.88         68.01  94.89
              0.000 1.0000        1.0000
 1 pond                     12.76         67.06  79.82
              0.250 0.9494        0.0130
 2 pond                     12.08          0.93  13.01
              0.500 0.8625        0.0002
 3 ice                       2.05          0.02   2.07
              0.050 0.7888        0.0000
 4 ice                      10.84          0.00  10.84
              0.375 0.6228        0.0000
 5 ice                       9.41          0.00   9.41
              0.750 0.4730        0.0000
 6 ice                       6.78          0.00   6.78
              1.125 0.3551        0.0000
 7 ice                       4.61          0.00   4.61
              1.500 0.2612        0.0000
 8 ocean                    27.57          0.00  27.57
```
Model parameters can be changed simply by modify `model` parameters, for example:

```
In [11]: model.pond_depth
Out[11]: 0.5

In [12]: model.pond_depth = 1.
```

A list of model parameters is given below.

  #### Spatio-temporal parameters:
  - :day_of_year: day of year, 1..365, where day 1 = January 1
  - :latitude: latitude (-90 to 90)  (test=80.)

  #### Surface characteristics:
  - :surface_pressure:  Surface pressure in mb (test=1008 mb)
  - :co2_volume_mixing_ratio:  CO2 volume mixing ratio (test 3.7e-04)
  - :surface_air_temperature:  Surface air temperature (K) (test=273.16 K)
  - :ground_temperature:  Surface skin temperature (K) (test=273.17 K)
  - :snow_depth:  Physical snow depth in meters (test=0 m)
  - :snow_density:  Snow density (kg/m3) (test=330 kg/m3)
  - :snow_grain_radius:  Snow grain radius in microns (um) (test=50. um)
  -:pond_depth:  Physical pond depth in meters (test=0.5 m)
  - :pond_tuning_parameter:  Pond tuning parameter in standard deviations (test=-1.)
  - :sea_ice_thickness:  Physical ice thickness in meters (test=1.5 m)
  - :sea_ice_tuning_parameter:  Sea ice tuning parameter in standard deviations (test=0.)

  #### Atmospheric Profile - 18 element array-like objects
  - :level: number id of level 1..18
  - :pressure: Pressure in mb
  - :air_temperature: air temperature in Kelvin
  - :water_vapor_mixing_ratio:  Water vapour mixing ration (g/g)
  - :ozone_mixing_ratio:  Ozone mixing ration (g/g)
  - :cloud_cover:  Cloud cover - non-dimension 0.-1.
  - :cloud_liquid_water_path: Cloud liquid water path (g/m2)

## Sensitivity Analysis
A sensitivity analsysis of selected parameters can be found [here](sensitivity_analysis.html)
