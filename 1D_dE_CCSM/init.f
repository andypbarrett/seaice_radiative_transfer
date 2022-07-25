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
      subroutine init_parameters()

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
     $     dayyr_in,        ! day of year
     $     rlat_in          ! latitude input

      real co2mix_in               ! co2 volume mixing ratio read in

      integer
     $     i,                   ! longitude point counter
     $     k,                   ! Level counter
     $     lev_in(plev)            ! level input

      real
     $     cld_in(plev),
     $     clwp_in(plev),
     $     hice_in,
     $     hpnd_in,
     $     o3mmr_in(plev),
     $     pmidm1_in(plev),
     $     ps_in,
     $     qm1_in(plev),
     $     rhos_in,
     $     rs_in,
     $     R_ice_in,
     $     R_pnd_in,
     $     sndpth_in,
     $     tg_in,
     $     tm1_in(plev),
     $     ts_in

      common /input/ dayyr_in, rlat_in, lev_in, pmidm1_in, tm1_in,
     $     qm1_in, o3mmr_in, cld_in, clwp_in, ps_in, co2mix_in,
     $     ts_in, tg_in, sndpth_in, rhos_in, rs_in, hpnd_in, R_pnd_in,
     $     hice_in, R_ice_in
      
      real 
     $     dayyr(plond),        ! day of year
     $     rlat(plond)          ! latitude input

      real co2mix               ! co2 volume mixing ratio read in

      integer
     $     lev(plev)            ! level input

      real
     $     cld(plond,plev),
     $     clwp(plond,plev),
     $     hice(plond),
     $     hpnd(plond),
     $     o3mmr(plond,plev),
     $     pmidm1(plond,plev),
     $     ps(plond),
     $     qm1(plond,plev),
     $     rhos(plond),
     $     rs(plond),
     $     R_ice(plond),
     $     R_pnd(plond),
     $     sndpth(plond),
     $     tg(plond),
     $     tm1(plond,plev),
     $     ts(plond)

      common /used/ dayyr, rlat, lev, pmidm1, tm1, qm1, o3mmr, cld,
     $     clwp, ps, co2mix, ts, tg, sndpth, rhos, rs, hpnd, R_pnd,
     $     hice, R_ice

C     Copy input variables to variables used in crm.  This will replace
C     the file read in get_data
      dayyr(1) = dayyr_in
      rlat(1) = rlat_in
      hice(1) = hice_in
      hpnd(1) = hpnd_in
      ps(1) = ps_in
      rhos(1) = rhos_in
      rs(1) = rs_in
      R_ice(1) = R_ice_in
      R_pnd(1) = R_pnd_in
      sndpth(1) = sndpth_in
      tg(1) = tg_in
      ts(1) = ts_in

      do 100 k=1, plev
         cld(1,k) = cld_in(k)
         clwp(1,k) = clwp_in(k)
         o3mmr(1,k) = o3mmr_in(k)
         pmidm1(1,k) = pmidm1_in(k)
         qm1(1,k) = qm1_in(k)
         tm1(1,k) = tm1_in(k)
 100  continue
      
      return
      end
