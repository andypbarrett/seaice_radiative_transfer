C----------------------------------------------------------------------
C     get_parameters: gets input parameters for ccsm2_sir_de sea ice
C     radiative transfer model.  Used as an interface with ctypes python
C     library
C
C
C     Inputs
C     ------
C     dayyr
C     rlat
C     lev        level input
C     pmidm1     pressure at model mid-levels
C     tm1        atmospheric temperature
C     qm1        moisture field
C     o3mmr      O3 mass mixing ratio
C     cld        cloud fraction - cannot be greater than 0.99999
C     clwp       cloud liquid water path (g/m**2)
C
C     ps         model surface pressure field
C     co2mix     
C     ts
C     tg
C     sndpth
C     rhos
C     rs
C     hpnd
C     R_pnd
C     hice
C     R_ice
C
C----------------------------------------------------------------------
      subroutine get_parameters()

      implicit none

      integer
     $     plon,
     $     plev,
     $     plevp,
     $     plond,   ! slt extended domain longitude
     $     nxpt    ! no.of points outside active domain for interpolant

      parameter(plon = 1,
     $          plev = 18,
     $          nxpt = 1,
     $          plevp = plev + 1,
     $          plond = plon + 1 + 2*nxpt)
      
      real 
     $     dayyr_in(plond),        ! day of year
     $     rlat_in(plond)          ! latitude input

      real co2mix               ! co2 volume mixing ratio read in

      integer
     $     i,                   ! longitude point counter
     $     k,                   ! Level counter
     $     lev_in(plev)            ! level input

      real
     $     cld_in(plond,plev),
     $     clwp_in(plond,plev),
     $     hice_in(plond),
     $     hpnd_in(plond),
     $     o3mmr_in(plond,plev),
     $     pmidm1_in(plond,plev),
     $     ps_in(plond),
     $     qm1_in(plond,plev),
     $     rhos_in(plond),
     $     rs_in(plond),
     $     R_ice_in(plond),
     $     R_pnd_in(plond),
     $     sndpth_in(plond),
     $     tg_in(plond),
     $     tm1_in(plond,plev),
     $     ts_in(plond)

      common /input/ dayyr_in, rlat_in, lev_in, pmidm1_in, tm1_in,
     $     qm1_in, o3mmr_in, cld_in, clwp_in, ps_in, co2mix_in,
     $     ts_in, tg_in, sndpth_in, rhos_in, rs_in, hpnd_in, R_pnd_in,
     $     hice_in, R_ice_in
      
      data lev_in /18.0, 17.0, 16.0, 15.0, 14.0, 13.0, 12.0, 11.0, 10.0,
     $     9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0/
      data pmidm1_in(1,:) /2.0, 5.0, 15.0, 35.0, 60.0, 105.0, 160.0,
     $     235.0, 320.0, 420.0, 520.0, 610.0, 710.0, 800.0, 870.0,
     $     930.0, 970.0, 1000.0/
      data tm1_in(1,:) /273.0, 251.0, 234.0, 226.0, 225.0, 225.0, 225.0,
     $     225.0, 234.0, 247.0, 257.0, 265.0, 272.0, 277.0, 280.0,
     $     281.0, 278.0, 276.0/
      data qm1_in(1,:) /4e-06, 4e-06, 4e-06, 4e-06, 4e-06, 4e-06,
     $     6.4e-06, 2.6e-05, 0.00012, 0.00052, 0.0011, 0.002, 0.0031,
     $     0.0042, 0.0051, 0.0059, 0.004, 0.003/
      data o3mmr_in(1,:) /7e-06, 1.3e-05, 1e-05, 5.5e-06, 4.2e-06,
     $     2.2e-06, 1e-06, 5e-07, 2e-07, 1.4e-07, 1e-07, 8e-08, 7e-08,
     $     6e-08, 5.5e-08, 5e-08, 4.5e-08, 4e-08/
      data cld_in(1,:) /0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     $     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0/
      data clwp_in(1,:) /0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     $     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 60.0, 0.0, 0.0/
c     
C  Make sure this check is done in getdat if( cld(i,k) .gt. 0.99999 ) cld(i,k) = .99999
c     
c
      i = 1
      
      dayyr_in(i) = 140.77
      rlat_in(i) = 80.

      ps_in(i) = 1008.
      co2mix_in = 3.7e-4
      ts_in(i) = 273.16
      tg_in(i) = 273.16
      sndpth_in(i) = 0.000
      rhos_in(i) = 330.000
      rs_in(i) = 50.00
      hpnd_in(i) = 0.500
      R_pnd_in(i) = -1.000
      hice_in(i) = 1.5
      R_ice_in(i) = 0.000

      return
      end
