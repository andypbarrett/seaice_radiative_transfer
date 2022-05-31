      subroutine print_parameters

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

      write(*,*) ' day of year (1..365)  = ',dayyr(i)
      write(*,*) ' latitude (-90 to +90) = ',rlat(i)
      do 200 k=1,plev
         write(6,99) k   ,pmidm1(i,k),tm1(i,k),qm1(i,k),o3mmr(i,k)
     +        ,cld(i,k),clwp(i,k)
 99      format(1x,i3,1x,6(1pe10.3,1x))
 200  continue

      write(*,571) ps(i)
 571  format('  surface pressure         = ',f7.2)
      write(*,572) co2mix
 572  format('  atmospheric co2 vmr      = ',1pe10.3)
      write(*,573)    ts(i)
 573  format('  surface air temperature  = ',f7.2)
      write(*,574)    tg(i)
 574  format('  surface skin temperature = ',f7.2)
      write(*,575)    sndpth(i)
 575  format('  snow physical depth (m)        = ',f9.4)
      write(*,576)    rhos(i)
 576  format('  snow density (kg/m3)           = ',f9.4)
      write(*,577)    rs(i)
 577  format('  snow grain radius (microns)    = ',f9.4)
      write(*,578)    hpnd(i)
 578  format('  pond physical depth (m)        = ',f9.4)
      write(*,579)    R_pnd(i)
 579  format('  pond tuning parameter          = ',f9.4)
      write(*,580)    hice(i)
 580  format('  sea ice thickness (m)          = ',f9.4)
      write(*,581)    R_ice(i)
 581  format('  sea ice tuning parameter       = ',f9.4)

      return
      end
      
      
