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

