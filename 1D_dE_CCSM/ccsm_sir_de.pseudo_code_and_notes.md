# Notes for ccsm_sir_de.for

## Order of routine

Parameter values are set throughout the code.

1. Open output file for writing

2. L524 Call `radini` to convert 4 constants from MKS
(meter-kilogram-second) to CGS (centimeter-gramme-second).  Routine
defined L9049.

3. L530 getdat - reads input from ccsm3_sir_de_input.dat and sets
other variables.

4. L539 Call zenith - calculate solar zenith angle (L11257)
   takes fractional calendar day, diurnal averaging boolean, current
   latitude
   returns: cosine of solar zenith angle

5. L534 Get albedo for land points - can this be ignored?

7. L551-561 Estimate fraction of direct to total nir - used in
albocean

8. L562 Get albedo for ocean - albocean (L2450) - does calculations
and writes to output file.
   Output:
     albs - surface albedo for direct rad (0.2-0.7)
     albl - surface albedo for direct rad (0.7-5.0)
     albsd - surface albedo for diffuse (0.2-0.7)
     albld - surface albedo for diffuse (0.7-5.0)
   Writes to output:
     albs
     albsd
     albl
     albld
     ...by level
     Fdirdn_vs - direct Flux at model interface for vs band
     Fdirup_vs
     Fdifdn_vs
     Fdifup_vs
     klmbda - irradiance extinction coefficient

9. L568 Get cloud particle and size and fraction of ice L4636
   Returns: 
     rel - liquid effective drop size
     rei - ice effective drop size
     fice - fractional ice content

10. L572. Calculate cloud emmisivity cldems L4776
   Returns:
     emis - cloud emissivity for atmospheric levels

11. L574 Calculate effoctive cloud cover

12. No cloud cover at the surface

13. L592 Call main radiation driving routine radctl L1336

14. The rest is a bunch of writing to output file - this could be
returned as variables.

--------------------------------------------------------------------------------------------
## Surface Radiative Transfer Model Notes

From Bruce Briegleb notes in ccsm3_sir_de.for

```
C Code was originally written for longitude bands of points,
C but has only been run for a single point for each submission.
C
C The main driver is program crm. The calling tree to the surface
C radiative transfer model is:
C
C    crm
C      albocean
C        simcsw
C          simded
C
C subroutine albocean    
C Sea ice/ocean shortwave radiation calculation of albedos
C and transmission/absorption in sea ice/ocean.
C
C subroutine simcsw
C Set up optical property profiles, based on snow, water and sea ice
C IOPs from Briegleb and Light, 2007 above.
C
C subroutine simded
C Delta-Eddington solution for snow/air/pond over sea ice
C
c To change number of snow or sea ice layers, change the line:
c     $          ksnow   = 1, kseaice =  4, klev = ksnow + kseaice + 1,
c
c where ksnow = number of thermodynamic snow layers, and
c kseaice = number of thermodynamic sea ice layers.
```

### Dependencies

To run `albocean` it looks like `radinit`, `getdat`, `zenith` need to be used.


### Common blocks in sea ice
radflux_seaice
      common/radflux_seaice/
     &              hi_ssl, hs_ssl
     &,             Fdirup_vs(plond,0:klevp),Fdirdn_vs(plond,0:klevp)
     &,             Fdifup_vs(plond,0:klevp),Fdifdn_vs(plond,0:klevp)
     &,             Fdirup_ni(plond,0:klevp),Fdirdn_ni(plond,0:klevp)
     &,             Fdifup_ni(plond,0:klevp),Fdifdn_ni(plond,0:klevp)
     &,             ksrf

seaice
      common/seaice/I_vs,I_ni,zd(0:klevp)
     &             ,Tri_vs(0:klevp),Tri_ni(0:klevp)
     &             ,Tro_vs,Tro_ni

albocean returns:
albs  -> asdir
albl  -> aldir
albsd -> asdif
albld -> aldif
** These are used by albland and albocean.  Results are written from within albocean
