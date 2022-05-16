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

