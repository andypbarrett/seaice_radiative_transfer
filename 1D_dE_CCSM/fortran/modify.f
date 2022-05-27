C----------------------------------------------------------------------
C     Modifies parameters in the input common block
C----------------------------------------------------------------------      
      subroutine modify

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
     $     dayyr(plond),        ! day of year
     $     rlat(plond)          ! latitude input

      real co2mix               ! co2 volume mixing ratio read in
      real dummy
      
      integer
     $     i,                   ! longitude point counter
     $     k,                   ! Level counter
     $     lev(plev)            ! level input

      real
     $     cld(plond,plevp),
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

      common /input/ dayyr, rlat, lev, pmidm1, tm1, qm1, o3mmr, cld,
     $     clwp, ps, co2mix, ts, tg, sndpth, rhos, rs,
     $     hpnd, R_pnd, hice, R_ice
      
      i = 1

      dummy = 9999.99
      
      do 10 k=1, plev
         cld(i,k) = dummy
         clwp(i,k) = dummy
         o3mmr(i,k) = dummy
         pmidm1(i,k) = dummy
         qm1(i,k) = dummy
         tm1(i,k) = dummy
 10   continue

      hice(i) = dummy
      hpnd(i) = dummy
      ps(i) = dummy
      rhos(i) = dummy
      rs(i) = dummy
      R_ice(i) = dummy
      R_pnd(i) = dummy
      sndpth(i) = dummy
      tg(i) = dummy
      ts(i) = dummy

      return
      end
      
