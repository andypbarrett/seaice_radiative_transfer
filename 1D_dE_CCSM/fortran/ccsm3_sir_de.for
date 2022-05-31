C Column Radiation Code used for the development of the
C CCSM Sea Ice Radiation Delta-Eddington treatment.
C   23 January 2007
C   Bruce P. Briegleb
C   NCAR
C
C This is an atmosphere radiative transfer model over
C the Delta-Eddington snow/bare ice/ponded ice radiative
C transfer model.
C
C The code was taken originally from the CCM (atmospheric) code,
C so it has lots of ancient parameters that are not active. The
C longwave (thermal) code is CCM, and still runs; it gives good
C fluxes, but not the latest from CCSM (CCM = Community Climate
C Model, the ancestor of the Community Atmosphere Model CAM).
C
C The Delta-Eddington snow/bare ice/ponded ice radiative 
C transfer model is described in:
C
C Briegleb, B. P., and B. Light (2007): A Delta-Eddington Multiple
C    Scattering Parameterization for Solar Radiation in the Sea Ice
C    Component of the Community Climate System Model, NCAR Technical
C    Note  NCAR/TN-472+STR  To be published in: February 2007
C
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
c
c 23 January 2007 Bruce P. Briegleb
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c-----------------------------NOTICE------------------------------------
c
c            NCAR COMMUNITY CLIMATE MODEL, VERSION 3.0
c            COPYRIGHT (C) 1996
c            UNIVERSITY CORPORATION FOR ATMOSPHERIC RESEARCH
c            ALL RIGHTS RESERVED
c
c               ------------ ----- --- ---------- ------
c  ********** | distribution terms and conditions notice | ************
c               ------------ ----- --- ---------- ------
c
c (c) copyright 1996 university corporation for atmospheric research/
c national center for atmospheric research/
c climate and global dynamics division
c
c this software, the community climate model (ccm), version ccm3, was
c developed by the climate and global dynamics division (cgd) climate
c modeling section (cms) of the national center for atmospheric research
c (ncar), which is operated by the university corporation for
c atmospheric research (ucar) and sponsored by the national science
c foundation (nsf).
c
c access and use of this software shall impose the following obligations
c and understandings on the user.  the user is granted the right,
c without any fee or cost, to use, copy, modify, alter, enhance and
c distribute this software, and any derivative works thereof, and its
c supporting documentation for any purpose whatsoever, except commercial
c sales, provided that this entire notice appears in all copies of the
c software, derivative works and supporting documentation.  further, the
c user agrees to credit ucar/ncar/cgd in any publications that result
c from the use of this software or in any software package that includes
c this software.  the names ucar/ncar/cgd, however, may not be used in
c any advertising or publicity to endorse or promote any products or
c commercial entity unless specific written permission is obtained from
c ucar/ncar/cgd.
c
c the ccm3 materials are made available with the understanding that
c ucar/ncar/cgd is not obligated to provide (and will not provide) the
c user with any support, consulting, training, or assistance of any kind
c with regard to the use, operation and performance of this software, nor
c to provide the user with any updates, revisions, new versions, or "bug
c fixes."
c
c this software is provided by ucar/ncar/cgd "as is" and any express or
c implied warranties, including but not limited to, the implied
c warranties of merchantability and fitness for a particular purpose are
c disclaimed.  in no event shall ucar/ncar/cgd be liable for any
c special, indirect or consequential damages or any damages whatsoever,
c including but not limited to claims associated with the loss of data
c or profits, which may result from an action in contract, negligence or
c other tortious claim that arises out of or in connection with the
c access, use or performance of this software.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
      subroutine crm()
c     
c     ccm3 radiation column model
c     Jeffrey Kiehl, Bruce Briegleb, and Charlie Zender  
c     May 1996
c
c     All routines except crm(), getdat(), radctl(), and four 
c     dummy routines (readric(), writeric(), outfld(), and radozn()) are 
c     included straight from CCM3 code. The purpose of the (non-dummy) 
c     routines is:
c
c     crm(): the main() routine. This routine takes the place of tphysbc()
c     in the ccm3 calling tree. It places the required calls to the cloud
c     routines cldefr() and cldems() directly, rather than invoking cldint().
c
c     getdat(): reads the stdin file to parse the user-specified column
c     profile. overwrites the co2vmr previsouly set by radini() with the
c     user specified value. computes o3vmr (instead of in radozn()).
c     sets the calday variable in comtim.h used in zenith() and radinp().
c
c     radctl(): main radiation driving routine, same as ccm3 version except:
c     receives additional input variable o3vmr, and passes out
c     additional diagnostic radiation quantities and eccf usually local to 
c     radctl(). does all cgs-->mks conversion for all fluxes.
c
c
c     The files that need to be either -included-, specified in the 
c     compilation statement, or linked from a ccm3 library are...
c
c     aermix.F    fmrgrid.F    radclr.F   radinp.F    trcabn.F    whenfgt.F
c     albland.F   intmax.F     radclw.F   radoz2.F    trcems.F    whenflt.F
c     albocean.F  isrchfgt.F   radcsw.F   radtpl.F    trcmix.F    whenne.F
c     blkdat.F    isrchfle.F   radded.F   resetr.F    trcplk.F    zenith.F
c     cldefr.F    myhandler.F  radems.F   torgrid.F   trcpth.F
c     cldems.F    radabs.F     radini.F   trcab.F     wheneq.F
c
c     
c     Standard input file is a text (ascii) file. At the end of the included
c     input file are notes describing more fully the input. The standard
c     output file is also a text (ascii) file. 
c
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
C------------------------------Parameters-------------------------------
c
c $Id: prgrid.h,v 1.2 1995/02/10 01:09:10 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation resolution and I/O parameters
C
      integer
     $     plon,    ! number of longitudes
     $     plev,    ! number of vertical levels
     $     plat,    ! number of latitudes
     $     pcnst,   ! number of constituents (including water vapor)
     $     plevmx,  ! number of subsurface levels
     $     plevp,   ! plev + 1
     $     nxpt,    ! no.of points outside active domain for interpolant
     $     jintmx,  ! number of extra latitudes in polar region
     $     plond,   ! slt extended domain longitude
     $     platd,   ! slt extended domain lat.
     $     plevd    ! fold plev,pcnst indices into one
      integer
     $     plngbuf  ! length of absorptivity/emissivity record
C
      parameter(plon    = 1,
     $          plev    = 18,
     $          plat    = 1,
     $          pcnst   = 1,
     $          plevmx  = 4,
     $          plevp   = plev + 1,
     $          nxpt    = 1,
     $          jintmx  = 1,
     $          plond   = plon + 1 + 2*nxpt,
     $          platd   = plat + 2*nxpt + 2*jintmx,
     $          plevd   = plev*(3 + pcnst))
      parameter(plngbuf = 512*((plond*plevp*plevp + plond*plev*4 +
     $                          plond*plevp)/512 + 1))
C
C Sea Ice Resolution parameters
C
      integer ksnow   ! Number of snow layers in CCSM
      integer kseaice ! Number of sea ice layers in CCSM
      integer klev    ! Number of layers minus one (0 layer at top) for rad calculation
      integer klevp   ! klev + 1
C
      parameter( ksnow   = 1, kseaice =  4, klev = ksnow + kseaice + 1,
     $           klevp   = klev + 1)
C
      integer kstrt   ! ice starting index for printout
      integer kend    ! ice ending index for printout
C------------------------------Commons----------------------------------
c
c $Id: comtim.h,v 1.1.1.1 1995/02/09 23:26:44 ccm2 Exp $
c $Author: ccm2 $
c
C
C Model time variables
C
      common/comtim/calday  ,dtime   ,twodt   ,nrstrt  ,nstep   ,
     $              nstepr  ,nestep  ,nelapse ,nstop   ,mdbase  ,
     $              msbase  ,mdcur   ,mscur   ,mbdate  ,mbsec   ,
     $              mcdate  ,mcsec   ,nndbas  ,nnsbas  ,nnbdat  ,
     $              nnbsec  ,doabsems,dosw    ,dolw
C
      real calday,   ! Current calendar day = julian day + fraction
     $     dtime,    ! Time step in seconds (delta t)
     $     twodt     ! 2 * delta t 
      integer
     $     nrstrt,   ! Starting time step of restart run (constant) 
     $     nstep,    ! Current time step
     $     nstepr,   ! Current time step of restart run(updated w/nstep)
     $     nestep,   ! Time step on which to stop run
     $     nelapse,  ! Requested elapsed time for model run
     $     nstop,    ! nestep + 1
     $     mdbase,   ! Base day of run
     $     msbase,   ! Base seconds of base day
     $     mdcur,    ! Current day of run
     $     mscur,    ! Current seconds of current day
     $     mbdate,   ! Base date of run (yymmdd format)
     $     mbsec,    ! Base seconds of base date
     $     mcdate,   ! Current date of run (yymmdd format)
     $     mcsec,    ! Current seconds of current date
     $     nndbas,   ! User input base day
     $     nnsbas,   ! User input base seconds of input base day
     $     nnbdat,   ! User input base date (yymmdd format)
     $     nnbsec    ! User input base seconds of input base date
      logical
     $     doabsems, ! True => abs/emiss calculation this timestep
     $     dosw,     ! True => shortwave calculation this timestep
     $     dolw      ! True => longwave calculation this timestep
C
c
c $Id: comctl.h,v 1.1.1.1 1995/02/09 23:26:41 ccm2 Exp $
c $Author: ccm2 $
c
C
C Model control variables
C
      common/comctl/itsst   ,nsrest  ,iradsw  ,iradlw  ,iradae  ,
     $              anncyc  ,nlend   ,nlres   ,nlhst   ,lbrnch  ,
     $              ldebug  ,aeres   ,ozncyc  ,sstcyc  ,dodiavg ,
     $              aeregen ,cpuchek
      integer
     $     itsst,   ! Sea surf. temp. update freq. (iters)
     $     nsrest,  ! Restart flag
     $     iradsw,  ! Iteration frequency for shortwave radiation computation
     $     iradlw,  ! Iteration frequency for longwave radiation computation
     $     iradae   ! Iteration freq. for absorptivity/emissivity comp
      logical
     $     anncyc,  ! Do annual cycle (otherwise perpetual)
     $     nlend,   ! Flag for end of run
     $     nlres,   ! If true, continuation run
     $     nlhst,   ! If true, regeneration run
     $     lbrnch,  ! If true, branch run
     $     ldebug,  ! If in debug mode, link output files to /usr/tmp
C                   !    before mswrite, and remove all but last file
     $     aeres,   ! If true, a/e data will be stored on restart file
     $     ozncyc,  ! If true, cycle ozone dataset
     $     sstcyc,  ! If true, cycle sst dataset
     $     dodiavg, ! true => diurnal averaging
     $     aeregen, ! true => absor/emis part of regeneration data
     $     cpuchek  ! If true, check remaining cpu time at each writeup
C
C------------------------------Arguments--------------------------------
c
c     Fields specified by the user in getdat()
c
      integer 
     $     ioro(plond)          ! land/ocean/sea ice flag
c
      real 
     $     clat,                ! Current latitude (radians)
     $     cld(plond,plevp),    ! fractional cloud cover
     $     clwp(plond,plev),    ! cloud liquid water path
     $     coslat,              ! cosine latitude
c
c     NB: o3mmr and o3vmr should be dimensioned (plond,plevr) if a different 
c     size radiation grid is used. Clashes between prgrid.h and ptrrgrid.h
c     (they both define plngbuf) prevent us from dimensioning anything by
c     plevr in this top level crm() routine.
c
     $     o3mmr(plond,plev),  ! Ozone volume mixing ratio
     $     o3vmr(plond,plev),  ! Ozone volume mixing ratio
     $     pilnm1(plond,plevp), ! natural log of pintm1
     $     pintm1(plond,plevp), ! model interface pressures
     $     pmidm1(plond,plev),  ! model level pressures
     $     pmlnm1(plond,plev),  ! natural log of pmidm1
     $     ps(plond),           ! surface pressure
     $     qm1(plond,plev),     ! model level specific humidity
     $     sndpth(plond),       ! snow physical depth 
     $     rhos(plond),         ! snow density
     $     rs(plond),           ! snow grain radius
     $     R_ice(plond),        ! sea ice standard deviation tuning parameter
     $     R_pnd(plond),        ! ponded ice standard deviation tuning parameter
     $     tg(plond),           ! surface (skin) temperature
     $     hpnd(plond),         ! pond depth (m)
     $     hice(plond),         ! sea ice thickness
     $     tm1(plond,plev),     ! model level temperatures
     $     ts(plond)            ! surface air temperature
c     
c     Fields computed from user input
c
      real 
     $     effcld(plond,plevp), ! effective cloud=cld*emis
     $     emis(plond,plev),    ! cloud emissivity
     $     fice(plond,plev),    ! fractional amount of ice
     $     rei(plond,plev),     ! ice particle size
     $     rel(plond,plev)      ! liquid effective drop size (microns)
      real 
     $     coszrs(plond),       ! cosine solar zenith angle
     $     eccf,                ! earth/sun distance factor 
     $     loctim(plond),       ! local time of solar computation
     $     srfrad(plond)        ! srf radiative heat flux
c
      real asdir(plond),        ! albedo: shortwave, direct
     $     asdif(plond),        ! albedo: shortwave, diffuse
     $     aldir(plond),        ! albedo: longwave, direct
     $     aldif(plond)         ! albedo: longwave, diffuse
c
C     
C     Output longwave arguments from radctl()
C
      real 
     $     flwds(plond),        ! Surface down longwave flux
     $     lwup(plond),         ! Surface up longwave flux from coupler
     $     qrl(plond,plev)      ! Longwave cooling rate
C     
C     Output shortwave arguments from radctl()
C
      real 
     $     fsns(plond),         ! Surface absorbed solar flux
     $     qrs(plond,plev),     ! Solar heating rate
     $     soll(plond),         ! Downward solar rad onto surface (lw direct)
     $     solld(plond),        ! Downward solar rad onto surface (lw diffuse)
     $     sols(plond),         ! Downward solar rad onto surface (sw direct)
     $     solsd(plond)         ! Downward solar rad onto surface (sw diffuse)
c
c     Additional CRM diagnostic output from radctl()
c
      real 
     $     flns(plond),         ! srf longwave cooling (up-dwn) flux
     $     flnsc(plond),        ! clr sky lw flx at srf (up-dwn)
     $     flnt(plond),         ! net outgoing lw flx at model top
     $     flntc(plond),        ! clr sky lw flx at model top
     $     fsnsc(plond),        ! clr sky surface abs solar flux
     $     fsnt(plond),         ! total column absorbed solar flux
     $     fsntc(plond),        ! clr sky total column abs solar flux
     $     solin(plond)         ! solar incident flux
C
C special for sea ice albedos:  27 October 2004
C
      real fnidr(plond)         ! fraction of direct to total nir down srf flux
      real tclrsf(plond,plevp)  ! Product of clr-sky fractions from top
c
c     Local workspace: These variables are not saved.
c
      real
     $     alb,                 ! for radiation output
     $     albc,                ! for radiation output
     $     cfl,                 ! for radiation output
     $     cfn,                 ! for radiation output
     $     clatm,               ! for radiation output
     $     cls,                 ! for radiation output
     $     clt,                 ! for radiation output
     $     fla,                 ! for radiation output
     $     fld,                 ! for radiation output
     $     flw,                 ! for radiation output
     $     frs,                 ! for radiation output
     $     fsnstm,              ! for radiation output
     $     hbuf,                ! history buffer
     $     pie,                 ! for radiation output
     $     pmb,                 ! for radiation output
     $     qrlday,              ! for radiation output
     $     qrsday,              ! for radiation output
     $     rlat,                ! for radiation output
     $     sab,                 ! for radiation output
     $     scf,                 ! for radiation output
     $     sol                  ! for radiation output
      real
     $     fsds,                ! for radiation output
     $     vsfdir,              ! visible fraction direct
     $     nifdir,              ! near-ir fraction direct
     $     vsfrac,              ! visible fraction total srf irradiance
     $     albsrf               ! broadband surface albedo
c     
      integer
     $     i,                   ! longitude index
     $     k,                   ! level index
     $     lat                  ! latitude row index 
c
c     Fundamental constants needed by radini()
c
      real 
     $     cpairx,              ! heat capacity dry air at constant prs (J/kg/K)
     $     epsilox,             ! ratio mean mol weight h2o to dry air
     $     gravx,               ! gravitational acceleration (m/s**2)
     $     stebolx              ! Sefan-Boltzmann constant (W/m**2/K**4)
C-----------------------------------------------------------------------
C Absorption data for sea ice
C
      real
     &   I_vs     ! frac transmission vs through sea ice surface layer
     &,  I_ni     ! frac transmission ni through sea ice surface layer
     &,  Tri_vs   ! frac transmission vs, surface to sea ice layer interface
     &,  Tri_ni   ! frac transmission ni, surface to sea ice layer interface
     &,  Tro_vs   ! frac transmission vs to ocean
     &,  Tro_ni   ! frac transmission ni to ocean
     &,  zd       ! interface depths for snow/pond and sea ice (from its own surface)
      common/seaice/I_vs,I_ni,zd(0:klevp)
     &             ,Tri_vs(0:klevp),Tri_ni(0:klevp)
     &             ,Tro_vs,Tro_ni
C
      real
     &  F_SW_vs     ! solar vs absorbed in sea ice
     &, F_SW_ni     ! solar ni absorbed in sea ice
     &, F_SW_srf_vs ! vs solar absorbed in sea ice surface layer
     &, F_SW_srf_ni ! ni solar absorbed in sea ice surface layer
     &, F_SW_ocn_vs ! vs solar absorbed in underlying ocean
     &, F_SW_ocn_ni ! ni solar absorbed in underlying ocean
     &, F_SW_srf    ! total solar absorbed in sea ice surface layer
     &, Q_SW_vs     ! solar vs absorbed in sea ice layer
     &, Q_SW_ni     ! solar ni absorbed in sea ice layer
     &, z           ! depth in sea ice (m)
     &, z0          ! offset depth in sea ice (m)
     &, dz          ! thickness of thermodynamic ice layer
     &, frac_vs     ! fraction of direct visible flux
     &, frac_ni     ! fraction of direct near-ir flux
      integer 
     &  nmb         ! counter for ice layers
c
C-----------------------------------------------------------------------
C Fluxes for sea ice
C
      real hi_ssl     ! sea ice surface scattering layer thickness (m)
      real hs_ssl     ! snow surface scattering layer thickness (m)
      integer ksrf    ! interface index for surface absorption
C
      real Fdirup_vs  ! Up   flux to dir beam at model interface vs band
      real Fdirdn_vs  ! Down flux to dir beam at model interface vs band
      real Fdifup_vs  ! Up   flux to dif beam at model interface vs band
      real Fdifdn_vs  ! Down flux to dif beam at model interface vs band
C
      real Fdirup_ni  ! Up   flux to dir beam at model interface ni band
      real Fdirdn_ni  ! Down flux to dir beam at model interface ni band
      real Fdifup_ni  ! Up   flux to dif beam at model interface ni band
      real Fdifdn_ni  ! Down flux to dif beam at model interface ni band
C
      common/radflux_seaice/
     &              hi_ssl, hs_ssl
     &,             Fdirup_vs(plond,0:klevp),Fdirdn_vs(plond,0:klevp)
     &,             Fdifup_vs(plond,0:klevp),Fdifdn_vs(plond,0:klevp)
     &,             Fdirup_ni(plond,0:klevp),Fdirdn_ni(plond,0:klevp)
     &,             Fdifup_ni(plond,0:klevp),Fdifdn_ni(plond,0:klevp)
     &,             ksrf
C
c------------------------------Externals--------------------------------
c
      external albland
      external albocean
      external blkdat
      external cldefr
      external cldems
      external getdat
      external radctl
      external radini
      external zenith
c
c-----------------------------------------------------------------------
c
      open(unit=6, file='ccsm3_sir_de_output.dat',status='unknown')
      write(6,3333)
 3333 format(' .... Begin CCSM3 Sea Ice Radiation',
     $       ' with Delta Eddington ....')
c
c     Set latitude index to 1
c
      lat=1
c
c     Set parameters in common block comtim.h: nstep,dosw,dolw,doabsems
c
      nstep  =  0
      dosw  =  .true.
      dolw  =  .true.
      doabsems = .true.
c
c     Set parameters in common block comctl.h: anncyc,dodiavg,iradsw,iradlw,iradae
c
      anncyc = .true.
      dodiavg = .false.
      iradsw = 1
      iradlw = 1
      iradae = 1
c
c     Set parameters required for call to radini():
c
      gravx   =   9.80616
      cpairx  =   1.00464e3
      epsilox =   0.622
      stebolx =   5.67e-8
c     
c     Given these four constants in MKS, radini() will define their CGS
c     equivalents, as well as setting many radiation parameters stored
c     in common blocks. radini() must be called before getdat(), because
c     the co2 mixing ratio set (by the user) in getdat() should overwrite
c     the default CCM3 co2 mixing ratio set by radini().
c
      call radini(gravx,cpairx,epsilox,stebolx)
c
      write(6,*) '.... read in profile data ....'
c
c     getdat() also sets calday (used in zenith() and radinp()).
c
      call getdat(
     $     clat,      cld,   clwp,  coslat,   ioro,
     $     loctim,  o3mmr,  o3vmr,  pilnm1, pintm1,
     $     pmidm1, pmlnm1,     ps,     qm1, sndpth,
     $     rhos,       rs,  R_ice,   R_pnd,   hpnd,
     $     hice,       tg,    tm1,      ts)     
c
c     Get coszrs: needed as input to albland(), albocean(), radctl()
c
      call zenith (calday  ,dodiavg ,clat    ,coszrs  )
c
c     Find the albedo for land points
c
      call albland(lat     ,ioro    ,sndpth  ,coszrs  ,asdir   ,
     $     aldir   ,asdif   ,aldif   )
c
c     Find the albedo for ocean/sea-ice points
c
c special: estimate fraction of direct to total nir surface downwards
c flux as equal to geometrical clear sky:
c
      do i=1,plon
         tclrsf(i,1) = 1.
      end do
      do k=1,plev
         do i=1,plon
            tclrsf(i,k+1) = tclrsf(i,k)*(1.-cld(i,k+1))
         end do
      end do
      do i=1,plon
         fnidr(i) = tclrsf(i,plev)
      end do
      call albocean(lat ,ioro ,fnidr ,sndpth ,rhos, rs, coszrs ,tg ,
     $     hpnd  ,R_pnd ,hice ,R_ice ,
     $     asdir ,aldir ,asdif ,aldif  )
C
C Cloud particle size and fraction of ice
C
      call cldefr(ioro, tm1, rel, rei, fice, ps, pmidm1)
C
C Cloud emissivity
C
      call cldems(clwp, fice, rei, emis)
C
C Effective cloud cover
C
      do k=1,plev
         do i=1,plon
            effcld(i,k) = cld(i,k)*emis(i,k)
         end do
      end do
C
C Cloud cover at surface interface always zero (for safety's sake)
C
      do i=1,plon
         effcld(i,plevp) = 0.
         cld(i,plevp)    = 0.
      end do
C
C     Main radiation driving routine. 
C     NB: All fluxes returned from radctl() have already been converted to MKS.
C
      call radctl(hbuf    ,clat    ,coslat  ,lat     ,ts      ,
     $     pmidm1  ,pintm1  ,pmlnm1  ,pilnm1  ,tm1     ,
     $     qm1     ,cld     ,effcld  ,clwp    ,coszrs  ,
     $     asdir   ,asdif   ,aldir   ,aldif   ,fsns    ,
     $     qrs     ,qrl     ,flwds   ,lwup    ,rel     ,
     $     rei     ,fice    ,sols    ,soll    ,solsd   ,
     $     solld   ,
c++csz
     $     fsnt,fsntc,fsnsc,flnt,flns,flntc,flnsc,solin, ! output 
     $     eccf,
     $     o3vmr) ! input
c--csz
c
      do i=1,plon
         srfrad(i)=fsns(i)+flwds(i)
      end do
c
c write out final results:
c
      write(6,888) 
 888  format(/'  .... Atmosphere/Surface Radiation Calculation ....')

      write(6,*) ' ---- solar and longwave results ---- '
      write(6,*) ' ------------------ '
c
      pie  = 4.*atan(1.)
      rlat = 180.*clat/pie
c
      do 100 i=1,1
c
        write(6,330) calday
 330    format('  calendar day of year      = ',f10.5)
        write(6,331) eccf
 331    format('  earth-sun distance factor = ',f10.5)
        write(6,332) rlat
 332    format('  earth latitude            = ',f10.5)
        write(6,333) coszrs(i)
 333    format('  cosine solar zenith angle = ',f10.5)
        write(6,334) loctim(i)
 334    format('  local solar time          = ',f10.5)
        sol = solin(i)
        write(6,335) sol
 335    format('  solar toa insolation      = ',f7.2,' wm-2')
        sab = fsnt(i)
        if( sol .le. 0. ) then
          alb = -999.
        else
          alb = (sol-sab)/sol   
        endif
        write(6,336) alb
 336    format('  solar toa albedo          = ',f10.5)
        write(6,337) sab
 337    format('  solar toa absorbed        = ',f7.2,' wm-2')
        frs = fsns(i)
        fsnstm = sab - frs
        write(6,338) fsnstm
 338    format('  solar absorbed atmosphere = ',f7.2,' wm-2')
        write(6,339) frs
 339    format('  solar absorbed surface    = ',f7.2,' wm-2')
        clt = fsntc(i)
        write(6,340) clt
 340    format('  solar clear toa absorbed  = ',f7.2,' wm-2')
        cls = fsnsc(i)
        clatm = clt - cls
        write(6,341) clatm
 341    format('  solar clear atm absorbed  = ',f7.2,' wm-2')
        write(6,342) cls
 342    format('  solar clear srf absorbed  = ',f7.2,' wm-2')
        scf = sab - clt
        write(6,343) scf
 343    format('  solar cloud forcing       = ',f7.2,' wm-2')
        if( sol .gt. 0. ) then
          albc = (sol - clt) / sol
          write(6,344) albc
 344      format('  solar clear sky albedo    = ',f10.5)
        endif
        write(6,*) ' ------------------ '
c
        flw = flnt(i)
        write(6,345) flw
 345    format('  longwave net up toa       = ',f7.2,' wm-2')
        fla = flns (i)
        write(6,346) fla
 346    format('  longwave net surface up   = ',f7.2,' wm-2')
        fld = flwds (i)
        write(6,347) fld
 347    format('  longwave down at surface  = ',f7.2,' wm-2')
        clt = flntc(i)
        write(6,348) clt
 348    format('  longwave clear outgoing   = ',f7.2,' wm-2')
        cls = flnsc(i)
        write(6,349) cls
 349    format('  longwave clear net srf    = ',f7.2,' wm-2')
        cfl = clt - flw
        write(6,350) cfl
 350    format('  longwave cloud forcing    = ',f7.2,' wm-2')
c
        write(6,*) ' ------------------ '
        cfn = scf + cfl
        write(6,351) cfn
 351    format('  net cloud forcing         = ',f7.2,' wm-2')
c
        write(6,*) ' ------------------------------ '
        write(6,*) ' ---- cloud particle radii and ice fraction ---- '
        write(6,*) '   level   pr(mb)  r_liq(um)  r_ice(um) ',
     $       ' f(ice)'
c
        do 300 k=1,plev
          pmb    = pmidm1(i,k) / 100.
          write(6,399) k,pmb,rel(i,k),rei(i,k),fice(i,k)
  399     format(2x,i4,2x,f8.3,1x,3(f9.3,2x))
  300   continue
c
        write(6,*) ' ------------------------------ '
        write(6,*) ' ---- heating rates in K/day ---- '
        write(6,*) '   level  pr(mb)    qrs        qrl',
     $       '        qnet'
c
        do 200 k=1,plev
          qrsday = qrs(i,k) * 86400.
          qrlday = qrl(i,k) * 86400.
          pmb    = pmidm1(i,k) / 100.
          write(6,199) k,pmb,qrsday,qrlday,qrsday+qrlday
  199     format(2x,i4,2x,f8.3,1x,3(f9.3,2x))
  200   continue
c
c surface data
c
        write(6,*) ' ------------------------------ '
        write(6,*) ' ---- surface absorption and albedos ---- '
c
        write(6,380) sols(1)
 380    format('  solar vs direct surface irradiance  = ',f7.2,' wm-2')
        write(6,381) solsd(1)
 381    format('  solar vs diffuse surface irradiance = ',f7.2,' wm-2')
        vsfdir = sols(1)/(sols(1)+solsd(1)) 
        write(6,384) vsfdir
 384    format('  vs fraction of direct irradiance    = ',f10.5)
        write(6,386) (1.-vsfdir)
 386    format('  vs fraction of diffuse irradiance   = ',f10.5)
        write(6,382) soll(1)
 382    format('  solar ni direct surface irradiance  = ',f7.2,' wm-2')
        write(6,383) solld(1)
 383    format('  solar ni diffuse surface irradiance = ',f7.2,' wm-2')
        nifdir = soll(1)/(soll(1)+solld(1)) 
        write(6,385) nifdir
 385    format('  ni fraction of direct irradiance    = ',f10.5)
        write(6,387) (1.-nifdir)
 387    format('  ni fraction of diffuse irradiance   = ',f10.5)

        fsds = sols(1)+solsd(1)+soll(1)+solld(1)
        write(6,388) fsds
 388    format('  total solar surface irradiance      = ',f7.2,' wm-2')
        vsfrac = (sols(1)+solsd(1))/fsds
        write(6,361) vsfrac
 361    format('  vs fraction of total irradiance     = ',f10.5)
        write(6,362) (1.-vsfrac)
 362    format('  ni fraction of total irradiance     = ',f10.5)

        write(6,389) frs
 389    format('  solar absorbed at surface           = ',f7.2,' wm-2')
        albsrf = 1.-(frs/fsds)
        write(6,390) albsrf
 390    format('  broad band surface albedo           = ',f10.5)
C
C Total absorption in sea ice
C
        F_SW_vs  = sols(1)*(1.-asdir(1)) + 
     &             solsd(1)*(1.-asdif(1))
        F_SW_ni  = soll(1)*(1.-aldir(1)) + 
     &             solld(1)*(1.-aldif(1))
C
C Compute I,Tr:
C
        F_SW_srf_vs =
     $        ((Fdirdn_vs(1,0)-Fdirup_vs(1,0))*sols(1) +
     $         (Fdifdn_vs(1,0)-Fdifup_vs(1,0))*solsd(1)) -
     $        ((Fdirdn_vs(1,ksrf)-Fdirup_vs(1,ksrf))*sols(1) +
     $         (Fdifdn_vs(1,ksrf)-Fdifup_vs(1,ksrf))*solsd(1))
        I_vs = 1. - (F_SW_srf_vs/F_SW_vs)
        frac_vs = sols(1) / (sols(1) + solsd(1))
        do k=0,klevp
          Tri_vs(k) = Fdirdn_vs(1,k)*frac_vs + 
     $                Fdifdn_vs(1,k)*(1.-frac_vs)
        enddo
        Tro_vs = (Fdirdn_vs(1,klevp)-Fdirup_vs(1,klevp))*frac_vs + 
     $           (Fdifdn_vs(1,klevp)-Fdifup_vs(1,klevp))*(1.-frac_vs)
C
        F_SW_srf_ni =
     $        ((Fdirdn_ni(1,0)-Fdirup_ni(1,0))*soll(1) +
     $         (Fdifdn_ni(1,0)-Fdifup_ni(1,0))*solld(1)) -
     $        ((Fdirdn_ni(1,ksrf)-Fdirup_ni(1,ksrf))*soll(1) +
     $         (Fdifdn_ni(1,ksrf)-Fdifup_ni(1,ksrf))*solld(1))
        I_ni = 1. - (F_SW_srf_ni/F_SW_ni)
        frac_ni = soll(1) / (soll(1) + solld(1))
        do k=0,klevp
          Tri_ni(k) = Fdirdn_ni(1,k)*frac_ni + 
     $                Fdifdn_ni(1,k)*(1.-frac_ni)
        enddo
        Tro_ni = (Fdirdn_ni(1,klevp)-Fdirup_ni(1,klevp))*frac_ni + 
     $           (Fdifdn_ni(1,klevp)-Fdifup_ni(1,klevp))*(1.-frac_ni)
C
        F_SW_srf = F_SW_vs*(1.-I_vs) + F_SW_ni*(1.-I_ni)
        write(6,391) F_SW_vs
 391    format('  solar vs absorbed in snow/sea ice   = ',f7.2,' wm-2')
        write(6,392) F_SW_ni
 392    format('  solar ni absorbed in snow/sea ice   = ',f7.2,' wm-2')
        write(6,393) F_SW_srf
 393    format('  total solar absorbed in srfc layer  = ',f7.2,' wm-2')
        write(6,394) I_vs
 394    format('  I_vs (frac vs abs pentrating srf)   = ',f10.5)
        write(6,395) I_ni  
 395    format('  I_ni (frac ni abs pentrating srf)   = ',f10.5)

        write(6,*) ' ------------------------------ '
        write(6,7788) 
 7788   format('  ---- snow/sea ice transmitted flux (Tr) ',
     $         'and absorption (Q) ---- ')
        write(6,*) 
     & ' level   depth  Tr_vs     Q_vs   Tr_ni    Q_ni     Q_total' 
        write(6,101) F_SW_vs*(1.-I_vs),F_SW_ni*(1.-I_ni),
     &               F_SW_vs*(1.-I_vs)+F_SW_ni*(1.-I_ni)
 101    format('  surface ',15x,f6.2,10x,f6.2,5x,f6.2)
c
c print out column layer absorption and interface transmission
c
        do k=0,klev
            write(6,111) zd(k),Tri_vs(k),Tri_ni(k)
 111        format(2x,8x,f5.3,2x,2(f6.4,11x))
            Q_SW_vs =
     $       ((Fdirdn_vs(1,k)-Fdirup_vs(1,k))*sols(1) + 
     $       (Fdifdn_vs(1,k)-Fdifup_vs(1,k))*solsd(1))
     $      -((Fdirdn_vs(1,k+1)-Fdirup_vs(1,k+1))*sols(1) + 
     $       (Fdifdn_vs(1,k+1)-Fdifup_vs(1,k+1))*solsd(1))
C
            Q_SW_ni =
     $       ((Fdirdn_ni(1,k)-Fdirup_ni(1,k))*soll(1) + 
     $       (Fdifdn_ni(1,k)-Fdifup_ni(1,k))*solld(1))
     $      -((Fdirdn_ni(1,k+1)-Fdirup_ni(1,k+1))*soll(1) + 
     $       (Fdifdn_ni(1,k+1)-Fdifup_ni(1,k+1))*solld(1))
C
            if( sndpth(1) .gt. .00001 ) then
              if( k .le. ksnow ) then
                write(6,322) k,Q_SW_vs,Q_SW_ni,
     $                         Q_SW_vs+Q_SW_ni
 322            format('  snow ',i2,16x,f6.2,10x,f6.2,5x,f6.2)
              else
                write(6,325) k,Q_SW_vs,Q_SW_ni,
     $                         Q_SW_vs+Q_SW_ni
              endif
            else 
              if( hpnd(1) .gt. .00001 ) then
                if( k .le. ksnow ) then
                  write(6,323) k,Q_SW_vs,Q_SW_ni,
     $                         Q_SW_vs+Q_SW_ni
 323              format('  pond ',i2,16x,f6.2,10x,f6.2,5x,f6.2)
                else
                  write(6,325) k,Q_SW_vs,Q_SW_ni,
     $                         Q_SW_vs+Q_SW_ni
                endif
              else
                if( k .le. ksnow ) then
                  write(6,324) k,Q_SW_vs,Q_SW_ni,
     $                         Q_SW_vs+Q_SW_ni
 324              format('  air  ',i2,16x,f6.2,10x,f6.2,5x,f6.2)
                else
                  write(6,325) k,Q_SW_vs,Q_SW_ni,
     $                         Q_SW_vs+Q_SW_ni
 325              format('  ice  ',i2,16x,f6.2,10x,f6.2,5x,f6.2)
                endif
              endif
            endif
        enddo
        write(6,111) zd(klevp),Tri_vs(klevp),Tri_ni(klevp)
        F_SW_ocn_vs =
     $        ((Fdirdn_vs(1,klevp)-Fdirup_vs(1,klevp))*sols(1) +
     $         (Fdifdn_vs(1,klevp)-Fdifup_vs(1,klevp))*solsd(1))
        F_SW_ocn_ni =
     $        ((Fdirdn_ni(1,klevp)-Fdirup_ni(1,klevp))*soll(1) +
     $         (Fdifdn_ni(1,klevp)-Fdifup_ni(1,klevp))*solld(1))
        write(6,223) F_SW_ocn_vs,F_SW_ocn_ni,F_SW_ocn_vs+F_SW_ocn_ni
 223    format(2x,'ocean',18x,f6.2,10x,f6.2,5x,f6.2)
c
 100  continue
c
c     Exiting the program without disabling the traps will cause an
c     IEEE note to appear in stderr even if no exceptions were detected
c     
c.. 061003 BPB      rcode=ieee_handler('clear','common',myhandler)
c     
      write(6,*) ' ------------------------------ '
      write(6,*) ' .... Column Radiation Calculation Completed ....'
c     
      return
      end
      subroutine getdat(
     $     clat,
     $     cld,
     $     clwp,
     $     coslat,
     $     ioro,
     $     loctim,
     $     o3mmr,
     $     o3vmr,
     $     pilnm1,
     $     pintm1,
     $     pmidm1,
     $     pmlnm1,
     $     ps,
     $     qm1, 
     $     sndpth,
     $     rhos,
     $     rs,
     $     R_ice,
     $     R_pnd,
     $     hpnd,
     $     hice,
     $     tg,
     $     tm1,
     $     ts)     

c-----------------------------------------------------------------------
c
c interface routine for column model that both initializes
c certain constants and reads external data:
c 
c o3 mass mixing ratios are read in, but the model also requires the
c path lengths; they are computed here
c
c also, from the cloud input (fraction and liquid water path), the
c cloud longwave emissivity must be computed; this is done here
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
c------------------------------Parameters-------------------------------
c
c $Id: prgrid.h,v 1.2 1995/02/10 01:09:10 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation resolution and I/O parameters
C
      integer
     $     plon,    ! number of longitudes
     $     plev,    ! number of vertical levels
     $     plat,    ! number of latitudes
     $     pcnst,   ! number of constituents (including water vapor)
     $     plevmx,  ! number of subsurface levels
     $     plevp,   ! plev + 1
     $     nxpt,    ! no.of points outside active domain for interpolant
     $     jintmx,  ! number of extra latitudes in polar region
     $     plond,   ! slt extended domain longitude
     $     platd,   ! slt extended domain lat.
     $     plevd    ! fold plev,pcnst indices into one
      integer
     $     plngbuf  ! length of absorptivity/emissivity record
C
      parameter(plon    = 1,
     $          plev    = 18,
     $          plat    = 1,
     $          pcnst   = 1,
     $          plevmx  = 4,
     $          plevp   = plev + 1,
     $          nxpt    = 1,
     $          jintmx  = 1,
     $          plond   = plon + 1 + 2*nxpt,
     $          platd   = plat + 2*nxpt + 2*jintmx,
     $          plevd   = plev*(3 + pcnst))
      parameter(plngbuf = 512*((plond*plevp*plevp + plond*plev*4 +
     $                          plond*plevp)/512 + 1))
C
c------------------------------Commons----------------------------------
c
c $Id: comtim.h,v 1.1.1.1 1995/02/09 23:26:44 ccm2 Exp $
c $Author: ccm2 $
c
C
C Model time variables
C
      common/comtim/calday  ,dtime   ,twodt   ,nrstrt  ,nstep   ,
     $              nstepr  ,nestep  ,nelapse ,nstop   ,mdbase  ,
     $              msbase  ,mdcur   ,mscur   ,mbdate  ,mbsec   ,
     $              mcdate  ,mcsec   ,nndbas  ,nnsbas  ,nnbdat  ,
     $              nnbsec  ,doabsems,dosw    ,dolw
C
      real calday,   ! Current calendar day = julian day + fraction
     $     dtime,    ! Time step in seconds (delta t)
     $     twodt     ! 2 * delta t 
      integer
     $     nrstrt,   ! Starting time step of restart run (constant) 
     $     nstep,    ! Current time step
     $     nstepr,   ! Current time step of restart run(updated w/nstep)
     $     nestep,   ! Time step on which to stop run
     $     nelapse,  ! Requested elapsed time for model run
     $     nstop,    ! nestep + 1
     $     mdbase,   ! Base day of run
     $     msbase,   ! Base seconds of base day
     $     mdcur,    ! Current day of run
     $     mscur,    ! Current seconds of current day
     $     mbdate,   ! Base date of run (yymmdd format)
     $     mbsec,    ! Base seconds of base date
     $     mcdate,   ! Current date of run (yymmdd format)
     $     mcsec,    ! Current seconds of current date
     $     nndbas,   ! User input base day
     $     nnsbas,   ! User input base seconds of input base day
     $     nnbdat,   ! User input base date (yymmdd format)
     $     nnbsec    ! User input base seconds of input base date
      logical
     $     doabsems, ! True => abs/emiss calculation this timestep
     $     dosw,     ! True => shortwave calculation this timestep
     $     dolw      ! True => longwave calculation this timestep
C
c
c $Id: crdalb.h,v 1.1.1.1 1995/02/09 23:26:46 ccm2 Exp $
c $Author: ccm2 $
c
C
C Surface albedo data
C
C The albedos are computed for a model grid box by ascribing values to
C 1x1 degree points of a vegetation dataset, then linearly averaging
C for each grid box; ocean and land values are averaged together along
C coastlines; the fraction of every grid box that has strong zenith
C angle dependence is included also (see Briegleb, Bruce P., 1992:
C Delta-Eddington Approximation for Solar Radiation in the NCAR
C Community Climate Model, Journal of Geophysical Research, Vol 97, D7,
C pp7603-7612).
C
      common/crdalb/albvss(plond,plat),albvsw(plond,plat),
     $              albnis(plond,plat),albniw(plond,plat),
     $              frctst(plond,plat)
C
C vis  = 0.2 - 0.7 micro-meters wavelength range
C nir  = 0.7 - 5.0 micro-meters wavelength range
C
C szad = strong zenith angle dependent
C wzad = weak   zenith angle dependent
C
      real albvss, ! Grid box alb for vis over szad surfaces
     $     albvsw, ! Grid box alb for vis over wzad surfaces
     $     albnis, ! Grid box alb for nir over szad surfaces
     $     albniw, ! Grid box alb for nir over wzad surfaces
     $     frctst  ! Fraction of area in grid box with szad surfaces
C
C Surface boundary data
C
C Vegtyp is used to specify the thermal properites of the surface, as
C well as determine the location of permanent land ice points; it is the
C dominant surface type within the model grid box based on the 1x1
C degree resolution vegetation dataset; it is encoded in the following
C manner:
C
C   1        ocean
C   2        sea ice
C   3        permanent land ice
C   4        tropical evergreen forest
C   5        deciduous forest
C   6        grassland/tundra
C   7        desert
C
C Rghnss is the aerodynamic roughness length for the grid box, computed
C by linear averaging of the values ascribed to the 1x1 degree
C resolution vegetation dataset; ocean and land values are averaged
C together along coastlines.
C
C Evapf is the ratio of actual to potential evaporation, and is computed
C from the 1x1 degree resolution vegetation dataset in a manner similar
C to the aerodynamic roughness.
C
C Vevapf allows for variable snow cover, where the underlying
C evaporability factor is modified.
C
C Snwjan and snwjly are mean climatological snow depths (liquid water
C equivalent) used to compute the prescribed daily values of snow cover.
C
      common/crdsrf/vegtyp(plond,plat),rghnss(plond,plat),
     $              evapf (plond,plat),vevapf(plond,plat),
     $              snwjan(plond,plat),snwjly(plond,plat)
C
      real vegtyp, ! Surface thermal type, based on veg type
     $     rghnss, ! Aerodynamic roughness length
     $     evapf , ! Constant surface evaporability
     $     vevapf, ! Variable surface evaporability
     $     snwjan, ! Snow cover (liq water equiv) for January
     $     snwjly  ! Snow cover (liq water equiv) for July
C
c
c $Id: crdcon.h,v 1.1.1.1 1995/02/09 23:26:46 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation constants
C
      common/crdcon/gravit  ,rga     ,cpair   ,epsilo  ,sslp    ,
     $              stebol  ,rgsslp  ,co2vmr  ,dpfo3   ,dpfco2  ,
     $              dayspy  ,pie
C
      real gravit,    ! Acceleration of gravity
     $     rga,       ! 1./gravit
     $     cpair,     ! Specific heat of dry air
     $     epsilo,    ! Ratio of mol. wght of H2O to dry air
     $     sslp,      ! Standard sea-level pressure
     $     stebol,    ! Stefan-Boltzmann's constant
     $     rgsslp,    ! 0.5/(gravit*sslp)
     $     co2vmr,    ! CO2 volume mixing ratio
     $     dpfo3,     ! Voigt correction factor for O3
     $     dpfco2,    ! Voigt correction factor for CO2
     $     dayspy,    ! Number of days per 1 year
     $     pie        ! 3.14.....
C
c-----------------------------------------------------------------------
c
c output arguments
c
      integer 
     $     ioro(plond)          ! land surface flag
c     
      real 
     $     clat,                ! model latitude in radians
     $     cld(plond,plevp),    ! cloud fraction
     $     clwp(plond,plev),    ! cloud liquid water path (g/m**2)
     $     coslat,              ! cosine latitude
     $     loctim(plond),       ! local time of solar computation
     $     o3mmr(plond,plev),   ! o3 mass mixing ratio
     $     o3vmr(plond,plev),   ! o3 mass mixing ratio
     $     pilnm1(plond,plevp), ! ln(pintm1)
     $     pintm1(plond,plevp), ! pressure at model interfaces 
     $     pmidm1(plond,plev),  ! pressure at model mid-levels 
     $     pmlnm1(plond,plev),  ! ln(pmidm1)
     $     ps(plond),           ! model surface pressure field
     $     qm1(plond,plev),     ! moisture field
     $     sndpth(plond),       ! snow physical depth 
     $     rhos(plond),         ! snow density
     $     rs(plond),           ! snow grain radius
     $     R_ice(plond),        ! sea ice standard deviation tuning parameter
     $     R_pnd(plond),        ! ponded ice standard deviation tuning parameter
     $     hpnd(plond),         ! pond depth (m)
     $     tg(plond),           ! surface (skin) temperature
     $     hice(plond),         ! sea ice thickness
     $     tm1(plond,plev),     ! atmospheric temperature
     $     ts(plond)            ! surface (air)  temperature
c     
c     local workspace
c     
      real 
     $     dayyr(plond),        ! day of year
     $     rlat(plond)          ! latitude input
c     
      real co2mix               ! co2 volume mixing ratio read in
c     
      integer 
     $     i,                   ! longitude index
     $     k,                   ! level  index
     $     lev(plev)            ! level input
c     
      character*80 label
c     
      real 
     $     amd,                 ! effective molecular weight of dry air (g/mol)
     $     amo,                 ! molecular weight of ozone (g/mol)
     $     vmmr                 ! ozone volume mixing ratio
c     
      data amd   /  28.9644   /
      data amo   /  48.0000   /
c
c-----------------------------------------------------------------------
      open(unit=5,file='ccsm3_sir_de_input.dat',status='old')
c     
c     begin read of data:
c     
      do 100 i=1,1
c     
         read(5,101)  label
 101     format(a80)
         write(6,*)   label
c     
         read(5,101)  label
         write(6,*)   label
c     
         read(5,101)  label
         write(6,*)   label
c     
         read(5,*)    dayyr(i)
         write(6,*) ' day of year (1..365)  = ',dayyr(i)
c     
         read(5,*)    rlat(i)
         write(6,*) ' latitude (-90 to +90) = ',rlat(i)
c     
         read(5,101)  label
         write(6,*)   label
c     
         do 200 k=1,plev
            read(5,*) lev(k),pmidm1(i,k),tm1(i,k),qm1(i,k),o3mmr(i,k)
     +           ,cld(i,k),clwp(i,k)
            write(6,99) k   ,pmidm1(i,k),tm1(i,k),qm1(i,k),o3mmr(i,k)
     +           ,cld(i,k),clwp(i,k)
 99         format(1x,i3,1x,6(1pe10.3,1x))
            if( cld(i,k) .gt. 0.99999 ) cld(i,k) = .99999
 200     continue
c
         read(5,*)     ps(i)
         write(6,571) ps(i)
 571     format('  surface pressure         = ',f7.2)
c     
c..........convert pressures from mb to pascals and define 
c..........interface pressures:
c     
         ps(i) = ps(i) * 100.
         do 125 k=1,plev
c     
            pmidm1(i,k) = pmidm1(i,k) * 100.
            pmlnm1(i,k) = alog(pmidm1(i,k))
c     
 125     continue
         do 150 k=1,plevp
c     
            if( k .eq. 1 ) then
               pintm1(i,k) = pmidm1(i,k) / 2.0
            else if ( k .gt. 1 .and. k .le. plev ) then
               pintm1(i,k) = 0.5 * (pmidm1(i,k-1) + pmidm1(i,k))
            else if ( k .eq. plevp ) then
               pintm1(i,k) = ps(i)
            endif
            pilnm1(i,k) = alog(pintm1(i,k))
c     
 150     continue
c     
         read(5,*)       co2mix
         write(6,572) co2mix
 572     format('  atmospheric co2 vmr      = ',1pe10.3)
c
c     NB: co2vmr is also set in radini() so this user value should override 
c     that setting.
c
         co2vmr = co2mix
c
         read(5,*)       ts(i)
         write(6,573)    ts(i)
 573     format('  surface air temperature  = ',f7.2)
c     
         read(5,*)       tg(i)
         write(6,574)    tg(i)
 574     format('  surface skin temperature = ',f7.2)
c     
c set these internally         read(5,*)       ioro(i)
c                              read(5,*)       rghnss(i,1)
c 2=sea ice
         ioro(i) = 2
c .04m
         rghnss(i,1) = .04
c
         read(5,*)       sndpth(i)
         write(6,575)    sndpth(i)
 575     format('  snow physical depth (m)        = ',f9.4)
c
         read(5,*)       rhos(i)
         write(6,576)    rhos(i)
 576     format('  snow density (kg/m3)           = ',f9.4)
c
         read(5,*)       rs(i)
         write(6,577)    rs(i)
 577     format('  snow grain radius (microns)    = ',f9.4)
c
         read(5,*)       hpnd(i)
         write(6,578)    hpnd(i)
 578     format('  pond physical depth (m)        = ',f9.4)
c
         read(5,*)       R_pnd(i)
         write(6,579)    R_pnd(i)
 579     format('  pond tuning parameter          = ',f9.4)
c     
         read(5,*)       hice(i)
         write(6,580)    hice(i)
 580     format('  sea ice thickness (m)          = ',f9.4)
c
         read(5,*)       R_ice(i)
         write(6,581)    R_ice(i)
 581     format('  sea ice tuning parameter       = ',f9.4)
c     
c..         read(5,*)       albvss(i,1)
c..         read(5,*)       albvsw(i,1)
c..         read(5,*)       albnis(i,1)
c..         read(5,*)       albniw(i,1)
c..        read(5,*)       frctst(i,1)
c
c 24 Oct 2006  Just set these to reasonable values:
c
       albvss(i,1) = 0.1
       albvsw(i,1) = 0.1
       albnis(i,1) = 0.1
       albniw(i,1) = 0.1
       frctst(i,1) = 0.0
c     
 100  continue
c     
      calday = dayyr(1)
      loctim(1) = (calday-aint(calday))*24.
      pie    = 4.*atan(1.)
      clat   = rlat(1)*(pie/180.)
      coslat = cos(clat)
C     
C     Convert ozone mass mixing ratio to ozone volume mixing ratio:
C     
      vmmr = amo/amd
      do k=1,plev
         do i=1,plon
c     o3mmr(i,k) = vmmr*o3vmr(i,k)
            o3vmr(i,k) = o3mmr(i,k)/vmmr
         end do
      end do
c     
      write(6,*) ' .... end of profile data input .... '
c     
      return
      end
      subroutine writeric(nabem,absems,lngbuf,nrow)
c
c........ dummy routine for write of ab/em data
c
      return 
      end
      subroutine readric(nabem,absems,lngbuf,nrow)
c
c........ dummy routine for read of ab/em data
c
      return 
      end
      subroutine outfld(name,tmp ,plond,lat,hbuf)
c
c........ dummy routine for history tape write
c
      return 
      end
      subroutine radozn(lat     ,pmid    ,o3vmr   )
c
c     Do nothing routine called from radctl. o3vmr is coming from
c     getdat() instead.
c
      return
      end
c
      subroutine radctl(hbuf    ,clat    ,coslat  ,lat     ,ts      ,
     $                  pmid    ,pint    ,pmln    ,piln    ,t       ,
     $                  h2ommr  ,cld     ,effcld  ,clwp    ,coszrs  ,
     $                  albs    ,albsd   ,albl    ,albld   ,fsns    ,
     $                  qrs     ,qrl     ,flwds   ,lwup    ,rel     ,
     $                  rei     ,fice    ,sols    ,soll    ,solsd   ,
     $                  solld   ,
c++csz
     $     fsnt,fsntc,fsnsc,flnt,flns,flntc,flnsc,solin, ! output
     $     eccf, ! output
     $     o3vmr) ! input
c--csz
C-----------------------------------------------------------------------
C
C Driver for radiation computation.
C
C Radiation uses cgs units, so conversions must be done from
C model fields to radiation fields.
C
C Calling sequence:
C
C     radinp      Converts units of model fields and computes ozone
C                 mixing ratio for solar scheme
C
C     radcsw      Performs solar computation
C       radalb    Computes surface albedos
C       radded    Computes delta-Eddington solution
C       radclr    Computes diagnostic clear sky fluxes
C
C     radclw      Performs longwave computation
C
C       radtpl    Computes path quantities
C       radems    Computes emissivity
C       radabs    Computes absorptivity
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Kiehl, B. Briegleb, August 1992
C
C Modified:          B. Briegleb, March 1995 to add aerosol
C                    to shortwave code
C
C-----------------------------------------------------------------------
c
c $Id: crm.F,v 1.9 1995/11/10 17:39:12 zender Exp $
c $Author: zender $
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
C------------------------------Parameters-------------------------------
c
c $Id: pmgrid.h,v 1.2 1995/02/10 01:09:06 ccm2 Exp $
c $Author: ccm2 $
c
C
C Basic grid point resolution parameters
C
      integer
     $     plon,    ! number of longitudes
     $     plev,    ! number of vertical levels
     $     plat,    ! number of latitudes
     $     pcnst,   ! number of constituents (including water vapor)
     $     plevmx,  ! number of subsurface levels
     $     plevp,   ! plev + 1
     $     nxpt,    ! no.of points outside active domain for interpolant
     $     jintmx,  ! number of extra latitudes in polar region
     $     plond,   ! slt extended domain longitude
     $     platd,   ! slt extended domain lat.
     $     plevd    ! fold plev,pcnst indices into one
C
      parameter(plon   = 1,
     $          plev   = 18,
     $          plat   = 1,
     $          pcnst  = 1,
     $          plevmx = 4,
     $          plevp  = plev + 1,
     $          nxpt   = 1,
     $          jintmx = 1,
     $          plond  = plon + 1 + 2*nxpt,
     $          platd  = plat + 2*nxpt + 2*jintmx,
     $          plevd  = plev*(3 + pcnst))
C
c
c $Id: ptrrgrid.h,v 1.1.1.1 1995/02/09 23:26:59 ccm2 Exp $
c $Author: ccm2 $
c
C
C Define radiation vertical grid and buffer length for abs/ems out-of-core file
C
      integer
     $     plevr,   ! number of vertical levels
     $     plevrp,  ! plevr + 1
     $     plngbuf  ! length of absorptivity/emissivity record
C
      parameter(plevr = 18,
     $          plevrp = plevr + 1,
     $          plngbuf = 512*((plond*plevrp*plevrp + plond*plevr*4 +
     $                          plond*plevrp)/512 + 1))
C
c
c $Id: pagrid.h,v 1.3 1995/03/03 17:47:17 bonan Exp $
c $Author: bonan $
c
C
C Model grid point resolution parameters.
C
      integer
     $     plnlv,    ! Length of multilevel field slice
     $     plndlv,   ! Length of multilevel 3-d field slice
     $     pbflnb,   ! Length of buffer 1
     $     pbflna,   ! Length of buffer 2
     $     pbflnm1,  ! Length of buffer m1
     $     pflenb,   ! Length of buffer 1, padded for unblocked I/O
     $     pflena,   ! Length of buffer 2, padded for unblocked I/O
     $     plenalcl, ! Length of buffer 2, needed in SPEGRD
     $     ptifld,   ! Number of fields on time-invariant boundary dataset
     $     ptvsfld,  ! Number of fields on time-variant boundary dataset
     $     ptvofld,  ! Number of fields on ozone dataset
     $     plenhi,   ! Length of integer header record
     $     plenhc,   ! Length of character header record
     $     plenhr,   ! Length of real header record
     $     plexbuf,  ! Length of communication buffer for flux coupling
     $     ptapes,   ! Maximum number of history tapes allowed
     $     pflds     ! Number of fields in master field list
      integer
     $     ptileni,  ! Length of time-invariant integer header
     $     ptilenc,  ! Length of time-invariant character header
     $     ptvoleni, ! Length of ozone integer header
     $     ptvolenc, ! Length of ozone character header
     $     ptvsleni, ! Length of time-variant integer header
     $     ptvslenc  ! Length of time-variant character header
      integer
     $     plenhis,  ! Length of integer header scalars
     $     plenhcs,  ! Length of character header scalars
     $     ptilenis, ! Length of time-invariant integer scalars
     $     ptilencs, ! Length of time-invariant character scalars
     $     ptolenis, ! Length of ozone integer header scalars
     $     ptolencs, ! Length of ozone character header scalars
     $     ptslenis, ! Length of time-variant integer header scalars
     $     ptslencs  ! Length of time-variant character header scalars
C
      parameter(plnlv=plon*plev,plndlv=plond*plev)
C
C In pbflnb, 9 multi-level fields include the plev levels of plol and
C plos. 2 multi-level fields are pcnst-dependent.
C There are plevmx sub-surface temperature fields. (See User's Guide 
C for complete buffer description)
C There are 4 single-level fields to hold albedos
C
      parameter(pbflnb=(7 + 2*pcnst)*plndlv + (15+plevmx+pcnst)*plond,
C
C In pbflna, there are 3 multi-level and 3 single-level fields.
C
     $          pbflna = (3 + 3*plev)*plond,
     $          pbflnm1 = (1 + 2*plev)*plond,
     $          pflenb = ((pbflnb + pbflnm1)/512 + 1)*512,
     $          pflena = (pbflna/512 + 1)*512,
C
C plenalcl is the buffer size as required in SPEGRD.  
C Only pflena is read/written.
C
     $          plenalcl = ((pbflna + 3*plndlv + plond)/512 + 1)*512,
     $          plexbuf = (((1 + 7*plev)*plond)/512+1)*512,
     $          ptapes = 6,
C
C 8 fields in master list are pcnst-dependent 2 fields occur only
C if pcnst > 1
C
     $          pflds=92+8*pcnst+2*(pcnst-1)+plevmx)
      parameter(ptifld = 11, ptvsfld = 1, ptvofld = 2)
C
C There are 37 scalar words in the integer header and 89 scalar words
C in the character header
C
      parameter(plenhis=37,plenhcs=89,
     $          plenhi=plenhis+3*pflds,plenhc=plenhcs+2*pflds,
     $          plenhr=3*(2*plev + 1) + 2*plat,
     $          ptilenis=plenhis, ptilencs=plenhcs,
     $          ptileni=ptilenis+3*ptifld, ptilenc=ptilencs+2*ptifld,
     $          ptolenis=plenhis, ptolencs=plenhcs,
     $          ptvoleni=ptolenis+3*ptvofld,ptvolenc=ptolencs+2*ptvofld,
     $          ptslenis=plenhis, ptslencs=plenhcs,
     $          ptvsleni=ptslenis+3*ptvsfld,ptvslenc=ptslencs+2*ptvsfld)







C------------------------------Commons----------------------------------
c
c $Id: comctl.h,v 1.1.1.1 1995/02/09 23:26:41 ccm2 Exp $
c $Author: ccm2 $
c
C
C Model control variables
C
      common/comctl/itsst   ,nsrest  ,iradsw  ,iradlw  ,iradae  ,
     $              anncyc  ,nlend   ,nlres   ,nlhst   ,lbrnch  ,
     $              ldebug  ,aeres   ,ozncyc  ,sstcyc  ,dodiavg ,
     $              aeregen ,cpuchek
      integer
     $     itsst,   ! Sea surf. temp. update freq. (iters)
     $     nsrest,  ! Restart flag
     $     iradsw,  ! Iteration frequency for shortwave radiation computation
     $     iradlw,  ! Iteration frequency for longwave radiation computation
     $     iradae   ! Iteration freq. for absorptivity/emissivity comp
      logical
     $     anncyc,  ! Do annual cycle (otherwise perpetual)
     $     nlend,   ! Flag for end of run
     $     nlres,   ! If true, continuation run
     $     nlhst,   ! If true, regeneration run
     $     lbrnch,  ! If true, branch run
     $     ldebug,  ! If in debug mode, link output files to /usr/tmp
C                   !    before mswrite, and remove all but last file
     $     aeres,   ! If true, a/e data will be stored on restart file
     $     ozncyc,  ! If true, cycle ozone dataset
     $     sstcyc,  ! If true, cycle sst dataset
     $     dodiavg, ! true => diurnal averaging
     $     aeregen, ! true => absor/emis part of regeneration data
     $     cpuchek  ! If true, check remaining cpu time at each writeup
C
C-----------------------------------------------------------------------
c
c $Id: comhst.h,v 1.1.1.1 1995/02/09 23:26:42 ccm2 Exp $
c $Author: ccm2 $
c
C
C Integer and logical variables related to history tapes
C
      integer pichsum  ! Max. value of 4*ichar(character)
      parameter (pichsum=508)
C
      common /comhst/
     $   nhtfrq(ptapes)    ,mfilt(ptapes) ,nlfilt             ,
     $   ndens(ptapes)     ,nflds(ptapes) ,nfils(ptapes)      ,
     $   hunit(ptapes)     ,nrlen(ptapes) ,nplen(ptapes)      ,
     $   sunit             ,stfnum        ,mtapes             ,
     $   nexcl             ,nincl         ,hbufpt(ptapes)     ,
     $   nacs(pflds,plat)  ,iflds(3,pflds),nupnt(pflds,ptapes),
     $   npnt(pflds,ptapes),ndcurf(ptapes),ncdatf(ptapes)     ,
     $   nscurf(ptapes)    ,ncsecf(ptapes),nfldsc(0:pichsum,ptapes),
     $   islocc(0:pichsum,ptapes)         ,hstwr(ptapes)      ,
     $   rstwr             ,nacsav(pflds,plat)
C
      integer nhtfrq,  ! Array of write frequencies
     $        mfilt    ! Number of write-ups per volume
      logical nlfilt   ! Flag for extra file on 1st vol (ktape=1)
      logical hstwr    ! Flag for history writes
      logical rstwr    ! Flag for restart writes
      integer ndens,   ! Array of input packing densities
     $        nflds,   ! Array of total fields on tape
     $        nfils,   ! Array of current files on the volume
     $        hunit,   ! History tape disk units
     $        nrlen,   ! Record length
     $        nplen,   ! Packed record length,
     $        sunit,   ! History tape SSD unit
     $        stfnum,  ! Starting number for history tape naming
     $        mtapes,  ! Actual number of tapes requested
     $        nexcl,   ! Actual number of excluded fields
     $        nincl,   ! Actual number of included primary tape fields
     $        hbufpt,  ! Ptrs to start of fields for each tape in hbuf
     $        nacs,    ! Number of accumulations for field
     $        nacsav,  ! Saved accumulations for restart
     $        iflds,   ! Integer portion of master field list
     $        nupnt,   ! Array of unpacked field pointers
     $        npnt,    ! Array of packed field pointers
     $        ndcurf,  ! First "current" day for each tape
     $        ncdatf,  ! First "current" date for each tape
     $        nscurf,  ! First "current" second of day for each tape
     $        ncsecf,  ! First "current" second of date for each tape
     $        nfldsc,  ! Number of fields starting with given ichar(1-4)
     $        islocc   ! Index of starting location for each ichar sum
C
C  Character variables related to history tapes
C
      common /comhtc/
     $   nfpath(ptapes)     ,ppath(ptapes)       ,cpath(ptapes)       ,
     $   nhfil(ptapes)      ,ninavg(ptapes)      ,caseid              ,
     $   ctitle             ,fieldn(2,pflds)     ,exclude(pflds)      ,
     $   primary(pflds)     ,aux(pflds,ptapes-1)
C
      character*80 nfpath,    ! Array of first pathnames, for header
     $             ppath,     ! Array of previous pathnames, for header
     $             cpath      ! Array of current pathnames
      character    nhfil*6,   ! Array of current file names
     $             ninavg*1,  ! Tape fields instantaneous or averaged
     $             caseid*8   ! Case identifier
      character*80 ctitle     ! Case title
      character*8  fieldn,    ! Character portion of master field list
     $             exclude,   ! List of fields to rm from primary tape
     $             primary,   ! List of fields to add to primary tape
     $             aux        ! Lists of fields for auxiliary tapes
C
C-----------------------------------------------------------------------
c
c $Id: comtim.h,v 1.1.1.1 1995/02/09 23:26:44 ccm2 Exp $
c $Author: ccm2 $
c
C
C Model time variables
C
      common/comtim/calday  ,dtime   ,twodt   ,nrstrt  ,nstep   ,
     $              nstepr  ,nestep  ,nelapse ,nstop   ,mdbase  ,
     $              msbase  ,mdcur   ,mscur   ,mbdate  ,mbsec   ,
     $              mcdate  ,mcsec   ,nndbas  ,nnsbas  ,nnbdat  ,
     $              nnbsec  ,doabsems,dosw    ,dolw
C
      real calday,   ! Current calendar day = julian day + fraction
     $     dtime,    ! Time step in seconds (delta t)
     $     twodt     ! 2 * delta t 
      integer
     $     nrstrt,   ! Starting time step of restart run (constant) 
     $     nstep,    ! Current time step
     $     nstepr,   ! Current time step of restart run(updated w/nstep)
     $     nestep,   ! Time step on which to stop run
     $     nelapse,  ! Requested elapsed time for model run
     $     nstop,    ! nestep + 1
     $     mdbase,   ! Base day of run
     $     msbase,   ! Base seconds of base day
     $     mdcur,    ! Current day of run
     $     mscur,    ! Current seconds of current day
     $     mbdate,   ! Base date of run (yymmdd format)
     $     mbsec,    ! Base seconds of base date
     $     mcdate,   ! Current date of run (yymmdd format)
     $     mcsec,    ! Current seconds of current date
     $     nndbas,   ! User input base day
     $     nnsbas,   ! User input base seconds of input base day
     $     nnbdat,   ! User input base date (yymmdd format)
     $     nnbsec    ! User input base seconds of input base date
      logical
     $     doabsems, ! True => abs/emiss calculation this timestep
     $     dosw,     ! True => shortwave calculation this timestep
     $     dolw      ! True => longwave calculation this timestep
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real hbuf(*)              ! history tape buffer,for calls to outfld
      integer lat               ! Latitude row index
      real ts(plond),           ! Surface (skin) temperature
     $     pmid(plond,plev),    ! Model level pressures
     $     pint(plond,plevp),   ! Model interface pressures
     $     pmln(plond,plev),    ! Natural log of pmid
     $     rel(plond,plev),     ! liquid cloud particle effective radius
     $     rei(plond,plev),     ! ice effective drop size (microns)
     $     fice(plond,plev),    ! fractional ice content within cloud
     $     piln(plond,plevp),   ! Natural log of pint
     $     t(plond,plev),       ! Model level temperatures
     $     h2ommr(plond,plev),  ! Model level specific humidity
     $     cld(plond,plevp),    ! Fractional cloud cover
     $     effcld(plond,plevp), ! Effective fractional cloud cover
     $     clwp(plond,plev)     ! Cloud liquid water path
      real coszrs(plond),       ! Cosine solar zenith angle
     $     albs(plond),
     $     albsd(plond),
     $     albl(plond),
     $     albld(plond)
      real clat,                  ! current latitude(radians)
     $     coslat                 ! cosine latitude
C
C Output solar arguments
C
      real fsns(plond),         ! Surface absorbed solar flux
     $     sols(plond),         ! Downward solar rad onto surface (sw direct)
     $     soll(plond),         ! Downward solar rad onto surface (lw direct)
     $     solsd(plond),        ! Downward solar rad onto surface (sw diffuse)
     $     solld(plond),        ! Downward solar rad onto surface (lw diffuse)
     $     qrs(plond,plev)      ! Solar heating rate
C
C Output longwave arguments
C
      real qrl(plond,plev),     ! Longwave cooling rate
     $     flwds(plond),        ! Surface down longwave flux
     $     lwup(plond)          ! Surface up longwave flux from coupler
C
C---------------------------Local variables-----------------------------
C
      integer i                 ! index
      real solin(plond),        ! Solar incident flux
     $     fsnt(plond),         ! Net column abs solar flux at model top
     $     fsntc(plond),        ! Clear sky total column abs solar flux
     $     fsnsc(plond)         ! Clear sky surface abs solar flux
      real flnt(plond),         ! Net outgoing lw flux at model top
     $     flns(plond),         ! Srf longwave cooling (up-down) flux
     $     flntc(plond),        ! Clear sky lw flux at model top
     $     flnsc(plond)         ! Clear sky lw flux at srf (up-down)
      real pbr(plond,plevr),    ! Model mid-level pressures (dynes/cm2)
     $     pnm(plond,plevrp),   ! Model interface pressures (dynes/cm2)
     $     o3vmr(plond,plevr),  ! Ozone volume mixing ratio
     $     o3mmr(plond,plevr),  ! Ozone mass mixing ratio
     $     plco2(plond,plevrp), ! Prs weighted CO2 path
     $     plh2o(plond,plevrp), ! Prs weighted H2O path
     $     tclrsf(plond,plevrp),! Total clear sky fraction, level to space
     $     eccf                 ! Earth/sun distance factor
      real tmp(plond)           ! Temporary for outfld
      real n2o(plond,plev),        ! nitrous oxide mass mixing ratio
     $     ch4(plond,plev),        ! methane mass mixing ratio
     $     cfc11(plond,plev),      ! cfc11 mass mixing ratio
     $     cfc12(plond,plev)       ! cfc12 mass mixing ratio
      real aermmr(plond,plevr), ! level aerosol mass mixing ratio
     $     rh(plond,plevr)      ! level relative humidity (fraction)
C
C Declare local arrays to which model input arrays are interpolated here.
C Current default is none since radiation grid = model grid.
C
C Externals.
C
      external radinp,          ! Computes latitude dependent radiation input
     $         aermix,          ! Specifies aerosol mass mixing ratio
     $         radcsw,          ! Computes solar radiation
     $         radozn,          ! Computes ozone volume mixing ratio
     $         radclw,          ! Computes longwave radiation
     $         torgrid,         ! Interpolate model variables to radiation grid
     $         fmrgrid          ! Interpolate radiation variables to model grid
C--------------------------------------------------------------------------
C
C Interpolate model input arrays to radiation vertical grid.  Currently this 
C is a do-nothing routine because radiation grid = model grid.
C
      call torgrid(pmid    ,pint    ,pmln    ,piln    ,t       ,
     $             h2ommr  ,cld     ,effcld  ,clwp    ,
     $             pmid    ,pint    ,pmln    ,piln    ,t       ,
     $             h2ommr  ,cld     ,effcld  ,clwp    )
C
C Interpolate ozone volume mixing ratio to model levels
C
c++csz
c
c     This is a do nothing routine in the CRM.
c     Instead of interpolating the o3vmr from the time-interpolated
c     values, we pass compute o3vmr in getdat() and pass it directly
c     into radctl(). o3mmr will be computed in radinp().
c
      call radozn(lat     ,pmid    ,o3vmr   )
c--csz
C
C Set latitude dependent radiation input
C
      call radinp(pmid    ,pint    ,h2ommr  ,cld     ,o3vmr   ,
     $            pbr     ,pnm     ,plco2   ,plh2o   ,tclrsf  ,
     $            eccf    ,o3mmr   )
C
C Solar radiation computation
C
      if (dosw) then
C
C Specify aerosol mass mixing ratio
C
         call aermix(pnm     ,pbr     ,h2ommr  ,t       ,aermmr  ,
     $               rh      )
C
         call radcsw(pnm     ,h2ommr  ,cld     ,clwp    ,o3mmr   ,
     $               eccf    ,coszrs  ,albs    ,albsd   ,albl    ,
     $               albld   ,solin   ,qrs     ,fsns    ,fsnt    ,
     $               fsnsc   ,fsntc   ,rel     ,rei     ,fice    ,
     $               sols    ,soll    ,solsd   ,solld   ,aermmr  ,
     $               rh      )
C
C Convert units of shortwave fields needed by rest of model from CGS to MKS
C
         do i=1,plon
            solin(i) = solin(i)*1.e-3
            fsnt(i)  = fsnt(i) *1.e-3
            fsns(i)  = fsns(i) *1.e-3
            fsntc(i) = fsntc(i)*1.e-3
            fsnsc(i) = fsnsc(i)*1.e-3
         end do
C
C Calculate/outfld albedo and clear sky albedo
C
         if (ninavg(1).eq.'Q') then
            do i=1,plon
               if (solin(i).gt.0.) then
                  tmp(i) = (solin(i) - fsnt(i)) / solin(i)
               else
                  tmp(i) = 0.
               end if
            end do
C           call outfld('ALB     ',tmp ,plond,lat,hbuf)
C
            do i=1,plon
               if (solin(i).gt.0.) then
                  tmp(i) = (solin(i) - fsntc(i)) / solin(i)
               else
                  tmp(i) = 0.
               end if
            end do
C           call outfld('ALBCLR  ',tmp ,plond,lat,hbuf)
         end if
C
C Dump shortwave radiation information to history tape buffer (diagnostics)
C
C        call outfld('SOLIN   ',solin ,plond,lat,hbuf)
C        call outfld('FSNT    ',fsnt  ,plond,lat,hbuf)
C        call outfld('FSNS    ',fsns  ,plond,lat,hbuf)
C        call outfld('FSNTC   ',fsntc ,plond,lat,hbuf)
C        call outfld('FSNSC   ',fsnsc ,plond,lat,hbuf)
      end if
C
C Longwave radiation computation
C
      if (dolw) then
c
c Specify trace gas mixing ratios    
c
         call trcmix(pmid, clat, coslat, n2o, ch4, cfc11, cfc12)
c
         call radclw(lat     ,ts      ,t       ,h2ommr  ,o3vmr   ,
     $               pbr     ,pnm     ,pmln    ,piln    ,plco2   ,
     $               plh2o   ,n2o     ,ch4     ,cfc11   ,cfc12   ,
     $               effcld  ,tclrsf  ,qrl     ,flns    ,flnt    ,
     $               flnsc   ,flntc   ,flwds   ,lwup    )
C
C Convert units of longwave fields needed by rest of model from CGS to MKS
C
         do i=1,plon
            flnt(i)  = flnt(i)*1.e-3
            flns(i)  = flns(i)*1.e-3
            flntc(i) = flntc(i)*1.e-3
            flnsc(i) = flnsc(i)*1.e-3
            flwds(i) = flwds(i)*1.e-3
         end do
C
C Dump longwave radiation information to history tape buffer (diagnostics)
C
C        call outfld('FLNT    ',flnt  ,plond,lat,hbuf)
C        call outfld('FLNTC   ',flntc ,plond,lat,hbuf)
C        call outfld('FLNS    ',flns  ,plond,lat,hbuf)
C        call outfld('FLNSC   ',flnsc ,plond,lat,hbuf)
      end if
C
C Interpolate radiation output arrays to model vertical grid.  Currently this 
C is a do-nothing routine because radiation grid = model grid.
C
      call fmrgrid(qrs     ,qrl     ,
     $             qrs     ,qrl     )
C
      return
      end
      subroutine aermix(pint  ,pmid   ,h2ommr  ,t   ,aermmr  ,
     $                  rh    )
C-----------------------------------------------------------------------
C
C Specify aerosol mixing ratio and compute relative humidity
C for later adjustment of aerosol optical properties
C
C Currently (March 1995), aerosol mass mixing ratio is specified in
C such a manner that the column visible aerosol optical depth is a
C specified global number (tauvis). This means that the actual mixing
C ratio depends on pressure thickness of the lowest three atmospheric
C layers near the surface.
C
C---------------------------Code history--------------------------------
C
C Original version:  B. Briegleb  March 1995
C
C-----------------------------------------------------------------------
c
c $Id: aermix.F,v 1.1 1995/03/17 18:54:06 ccm2 Exp $
c $Author: ccm2 $
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
C------------------------------Parameters-------------------------------
c
c $Id: pmgrid.h,v 1.2 1995/02/10 01:09:06 ccm2 Exp $
c $Author: ccm2 $
c
C
C Basic grid point resolution parameters
C
      integer
     $     plon,    ! number of longitudes
     $     plev,    ! number of vertical levels
     $     plat,    ! number of latitudes
     $     pcnst,   ! number of constituents (including water vapor)
     $     plevmx,  ! number of subsurface levels
     $     plevp,   ! plev + 1
     $     nxpt,    ! no.of points outside active domain for interpolant
     $     jintmx,  ! number of extra latitudes in polar region
     $     plond,   ! slt extended domain longitude
     $     platd,   ! slt extended domain lat.
     $     plevd    ! fold plev,pcnst indices into one
C
      parameter(plon   = 1,
     $          plev   = 18,
     $          plat   = 1,
     $          pcnst  = 1,
     $          plevmx = 4,
     $          plevp  = plev + 1,
     $          nxpt   = 1,
     $          jintmx = 1,
     $          plond  = plon + 1 + 2*nxpt,
     $          platd  = plat + 2*nxpt + 2*jintmx,
     $          plevd  = plev*(3 + pcnst))
C
c
c $Id: ptrrgrid.h,v 1.1.1.1 1995/02/09 23:26:59 ccm2 Exp $
c $Author: ccm2 $
c
C
C Define radiation vertical grid and buffer length for abs/ems out-of-core file
C
      integer
     $     plevr,   ! number of vertical levels
     $     plevrp,  ! plevr + 1
     $     plngbuf  ! length of absorptivity/emissivity record
C
      parameter(plevr = 18,
     $          plevrp = plevr + 1,
     $          plngbuf = 512*((plond*plevrp*plevrp + plond*plevr*4 +
     $                          plond*plevrp)/512 + 1))
C
c
c $Id: pagrid.h,v 1.3 1995/03/03 17:47:17 bonan Exp $
c $Author: bonan $
c
C
C Model grid point resolution parameters.
C
      integer
     $     plnlv,    ! Length of multilevel field slice
     $     plndlv,   ! Length of multilevel 3-d field slice
     $     pbflnb,   ! Length of buffer 1
     $     pbflna,   ! Length of buffer 2
     $     pbflnm1,  ! Length of buffer m1
     $     pflenb,   ! Length of buffer 1, padded for unblocked I/O
     $     pflena,   ! Length of buffer 2, padded for unblocked I/O
     $     plenalcl, ! Length of buffer 2, needed in SPEGRD
     $     ptifld,   ! Number of fields on time-invariant boundary dataset
     $     ptvsfld,  ! Number of fields on time-variant boundary dataset
     $     ptvofld,  ! Number of fields on ozone dataset
     $     plenhi,   ! Length of integer header record
     $     plenhc,   ! Length of character header record
     $     plenhr,   ! Length of real header record
     $     plexbuf,  ! Length of communication buffer for flux coupling
     $     ptapes,   ! Maximum number of history tapes allowed
     $     pflds     ! Number of fields in master field list
      integer
     $     ptileni,  ! Length of time-invariant integer header
     $     ptilenc,  ! Length of time-invariant character header
     $     ptvoleni, ! Length of ozone integer header
     $     ptvolenc, ! Length of ozone character header
     $     ptvsleni, ! Length of time-variant integer header
     $     ptvslenc  ! Length of time-variant character header
      integer
     $     plenhis,  ! Length of integer header scalars
     $     plenhcs,  ! Length of character header scalars
     $     ptilenis, ! Length of time-invariant integer scalars
     $     ptilencs, ! Length of time-invariant character scalars
     $     ptolenis, ! Length of ozone integer header scalars
     $     ptolencs, ! Length of ozone character header scalars
     $     ptslenis, ! Length of time-variant integer header scalars
     $     ptslencs  ! Length of time-variant character header scalars
C
      parameter(plnlv=plon*plev,plndlv=plond*plev)
C
C In pbflnb, 9 multi-level fields include the plev levels of plol and
C plos. 2 multi-level fields are pcnst-dependent.
C There are plevmx sub-surface temperature fields. (See User's Guide 
C for complete buffer description)
C There are 4 single-level fields to hold albedos
C
      parameter(pbflnb=(7 + 2*pcnst)*plndlv + (15+plevmx+pcnst)*plond,
C
C In pbflna, there are 3 multi-level and 3 single-level fields.
C
     $          pbflna = (3 + 3*plev)*plond,
     $          pbflnm1 = (1 + 2*plev)*plond,
     $          pflenb = ((pbflnb + pbflnm1)/512 + 1)*512,
     $          pflena = (pbflna/512 + 1)*512,
C
C plenalcl is the buffer size as required in SPEGRD.  
C Only pflena is read/written.
C
     $          plenalcl = ((pbflna + 3*plndlv + plond)/512 + 1)*512,
     $          plexbuf = (((1 + 7*plev)*plond)/512+1)*512,
     $          ptapes = 6,
C
C 8 fields in master list are pcnst-dependent 2 fields occur only
C if pcnst > 1
C
     $          pflds=92+8*pcnst+2*(pcnst-1)+plevmx)
      parameter(ptifld = 11, ptvsfld = 1, ptvofld = 2)
C
C There are 37 scalar words in the integer header and 89 scalar words
C in the character header
C
      parameter(plenhis=37,plenhcs=89,
     $          plenhi=plenhis+3*pflds,plenhc=plenhcs+2*pflds,
     $          plenhr=3*(2*plev + 1) + 2*plat,
     $          ptilenis=plenhis, ptilencs=plenhcs,
     $          ptileni=ptilenis+3*ptifld, ptilenc=ptilencs+2*ptifld,
     $          ptolenis=plenhis, ptolencs=plenhcs,
     $          ptvoleni=ptolenis+3*ptvofld,ptvolenc=ptolencs+2*ptvofld,
     $          ptslenis=plenhis, ptslencs=plenhcs,
     $          ptvsleni=ptslenis+3*ptvsfld,ptvslenc=ptslencs+2*ptvsfld)







C------------------------------Commons----------------------------------
c
c $Id: crdcon.h,v 1.1.1.1 1995/02/09 23:26:46 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation constants
C
      common/crdcon/gravit  ,rga     ,cpair   ,epsilo  ,sslp    ,
     $              stebol  ,rgsslp  ,co2vmr  ,dpfo3   ,dpfco2  ,
     $              dayspy  ,pie
C
      real gravit,    ! Acceleration of gravity
     $     rga,       ! 1./gravit
     $     cpair,     ! Specific heat of dry air
     $     epsilo,    ! Ratio of mol. wght of H2O to dry air
     $     sslp,      ! Standard sea-level pressure
     $     stebol,    ! Stefan-Boltzmann's constant
     $     rgsslp,    ! 0.5/(gravit*sslp)
     $     co2vmr,    ! CO2 volume mixing ratio
     $     dpfo3,     ! Voigt correction factor for O3
     $     dpfco2,    ! Voigt correction factor for CO2
     $     dayspy,    ! Number of days per 1 year
     $     pie        ! 3.14.....
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real pint(plond,plevrp),  ! Radiation level interface pressures (dynes/cm2)
     $     pmid(plond,plevr),   ! Radiation level mid-level pressures (dynes/cm2)
     $     h2ommr(plond,plev),  ! Radiation level specific humidity   (g/g)
     $     t(plond,plev)        ! Radiation level temperatures        (K)
C
C Output arguments
C
      real aermmr(plond,plevr), ! Radiation level aerosol mass mixing ratio
     $     rh(plond,plevr)      ! Radiation level relative humidity (fraction)
C
C---------------------------Local variables-----------------------------
C
      integer i,     ! longitude index
     $        k,     ! level index
     $        mxaerl ! max nmbr aerosol levels counting up from surface
C
      real tauvis,  ! visible optical depth
     $     kaervs,  ! visible extinction coefficiant of aerosol (m2/g)
     $     omgvis,  ! visible omega0
     $     gvis,    ! visible forward scattering asymmetry parameter
     $     rhcnst   ! constant relative humidity factor
C
C Relative humidity factor
C
      real rhfac,              ! multiplication factor for kaer
     $     rhpc,               ! level relative humidity in %
     $     a0,                 ! constant in rh mult factor
     $     a1,                 ! constant in rh mult factor
     $     a2,                 ! constant in rh mult factor
     $     a3                  ! constant in rh mult factor
c
      data a0 / -9.2906106183    /
      data a1 /  0.52570211505   /
      data a2 / -0.0089285760691 /
      data a3 /  5.0877212432e-05/
C
      data mxaerl /  3  /
      data tauvis / .12 /
      data kaervs / 5.3012 /
      data omgvis / 0.999999 /
      data gvis   / 0.694889 /
      data rhcnst / .80 /
C
C--------------------------------------------------------------------------
C
C Set relative humidity and factor; then aerosol amount:
C
      do i=1,plon
        do k=1,plevr
C
          rh(i,k) = rhcnst
C
C Compute relative humidity factor:
C
          if( rh(i,k) .gt. .90 ) then
            rhfac = 2.8
          else if (rh(i,k) .lt. .60 ) then
            rhfac = 1.0
          else
            rhpc  = 100. * rh(i,k)
            rhfac = (a0 + a1*rhpc + a2*rhpc**2 + a3*rhpc**3)
          endif
C
C Find constant aerosol mass mixing ratio for specified levels
C in the column, converting units where appropriate
C
          if( k .ge. plevrp-mxaerl ) then
            aermmr(i,k) = gravit*tauvis /
     $              (1.e4*kaervs*rhfac*(1.-omgvis*gvis*gvis)*
     $              (pint(i,plevrp)-pint(i,plevrp-mxaerl)))
          else
            aermmr(i,k) = 0.0
          endif
C
        enddo
      enddo        
C
C
      return
      end
      subroutine albland(lat     ,ioro    ,sndpth  ,coszrs  ,albs    ,
     $                   albl    ,albsd   ,albld   )
C-----------------------------------------------------------------------
C
C Compute surface albedos over land
C
C Computes surface albedos for direct/diffuse incident radiation for
C two spectral intervals:
C   s = 0.2-0.7 micro-meters
C   l = 0.7-5.0 micro-meters
C
C Uses knowledge of surface type to specify albedo, as follows:
C
C Land without    Albedos specified by two dimensional surface albedo
C     snow        fields, which distinguish surfaces with strong solar
C                 zenith angle dependence from those with weaker solar
C                 zenith angle dependence; alb independent of surface
C                 moisture or other physical factors.
C
C Land with snow  Snow depth (liquid water equivalent) used, along with
C                 aerodynamic roughness to define a horizontal fraction
C                 of land surface covered with snow; snow albedos are
C                 comptd as functions of solar zenith angle; these snow
C                 albedos are then weighted by the horizontal fraction
C                 of coverage with the underlying surface albedos
C                 computed above to produce total grid mean albedo.
C
C land with ice   Surface albedos specified as functions of spectral
C                 interval; combined with overlying snow in a similar
C                 manner to the case of land with snow.
C
C Note, the code collects together surfaces of the same type for various
C computations in order to vectorize longitude loops.
C
C For more details , see Briegleb, Bruce P., 1992: Delta-Eddington
C Approximation for Solar Radiation in the NCAR Community Climate Model,
C Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
C
C The details of the land surface albedo arrays can be found in the
C common block description below.
C
C---------------------------Code history--------------------------------
C
C Original version:        CCM1
C Standardized:            J. Rosinski, June 1992
C Reviewed:                J. Kiehl, B. Briegleb, August 1992
C Rewritten for land only: J. Rosinski, May, 1994
C
C-----------------------------------------------------------------------
c
c $Id: albland.F,v 1.1.1.1 1995/02/09 23:26:33 ccm2 Exp $
c $Author: ccm2 $
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
C------------------------------Parameters-------------------------------
c
c $Id: prgrid.h,v 1.2 1995/02/10 01:09:10 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation resolution and I/O parameters
C
      integer
     $     plon,    ! number of longitudes
     $     plev,    ! number of vertical levels
     $     plat,    ! number of latitudes
     $     pcnst,   ! number of constituents (including water vapor)
     $     plevmx,  ! number of subsurface levels
     $     plevp,   ! plev + 1
     $     nxpt,    ! no.of points outside active domain for interpolant
     $     jintmx,  ! number of extra latitudes in polar region
     $     plond,   ! slt extended domain longitude
     $     platd,   ! slt extended domain lat.
     $     plevd    ! fold plev,pcnst indices into one
      integer
     $     plngbuf  ! length of absorptivity/emissivity record
C
      parameter(plon    = 1,
     $          plev    = 18,
     $          plat    = 1,
     $          pcnst   = 1,
     $          plevmx  = 4,
     $          plevp   = plev + 1,
     $          nxpt    = 1,
     $          jintmx  = 1,
     $          plond   = plon + 1 + 2*nxpt,
     $          platd   = plat + 2*nxpt + 2*jintmx,
     $          plevd   = plev*(3 + pcnst))
      parameter(plngbuf = 512*((plond*plevp*plevp + plond*plev*4 +
     $                          plond*plevp)/512 + 1))
C
C-----------------------------------------------------------------------
c
c $Id: albedo.h,v 1.1.1.1 1995/02/09 23:26:33 ccm2 Exp $
c $Author: ccm2 $
c
      real snws,           ! Snow albedo for 0.2-0.7 micro-meters
     $     snwl,           ! Snow albedo for 0.7-5.0 micro-meters
     $     sices,          ! Sea ice albedo for 0.2-0.7 micro-meters
     $     sicel           ! Sea ice albedo for 0.7-5.0 micro-meters
      parameter (snws  = 0.95,
     $           snwl  = 0.70,
     $           sices = 0.70,
     $           sicel = 0.50)
C------------------------------Commons----------------------------------
c
c $Id: crdalb.h,v 1.1.1.1 1995/02/09 23:26:46 ccm2 Exp $
c $Author: ccm2 $
c
C
C Surface albedo data
C
C The albedos are computed for a model grid box by ascribing values to
C 1x1 degree points of a vegetation dataset, then linearly averaging
C for each grid box; ocean and land values are averaged together along
C coastlines; the fraction of every grid box that has strong zenith
C angle dependence is included also (see Briegleb, Bruce P., 1992:
C Delta-Eddington Approximation for Solar Radiation in the NCAR
C Community Climate Model, Journal of Geophysical Research, Vol 97, D7,
C pp7603-7612).
C
      common/crdalb/albvss(plond,plat),albvsw(plond,plat),
     $              albnis(plond,plat),albniw(plond,plat),
     $              frctst(plond,plat)
C
C vis  = 0.2 - 0.7 micro-meters wavelength range
C nir  = 0.7 - 5.0 micro-meters wavelength range
C
C szad = strong zenith angle dependent
C wzad = weak   zenith angle dependent
C
      real albvss, ! Grid box alb for vis over szad surfaces
     $     albvsw, ! Grid box alb for vis over wzad surfaces
     $     albnis, ! Grid box alb for nir over szad surfaces
     $     albniw, ! Grid box alb for nir over wzad surfaces
     $     frctst  ! Fraction of area in grid box with szad surfaces
C
C Surface boundary data
C
C Vegtyp is used to specify the thermal properites of the surface, as
C well as determine the location of permanent land ice points; it is the
C dominant surface type within the model grid box based on the 1x1
C degree resolution vegetation dataset; it is encoded in the following
C manner:
C
C   1        ocean
C   2        sea ice
C   3        permanent land ice
C   4        tropical evergreen forest
C   5        deciduous forest
C   6        grassland/tundra
C   7        desert
C
C Rghnss is the aerodynamic roughness length for the grid box, computed
C by linear averaging of the values ascribed to the 1x1 degree
C resolution vegetation dataset; ocean and land values are averaged
C together along coastlines.
C
C Evapf is the ratio of actual to potential evaporation, and is computed
C from the 1x1 degree resolution vegetation dataset in a manner similar
C to the aerodynamic roughness.
C
C Vevapf allows for variable snow cover, where the underlying
C evaporability factor is modified.
C
C Snwjan and snwjly are mean climatological snow depths (liquid water
C equivalent) used to compute the prescribed daily values of snow cover.
C
      common/crdsrf/vegtyp(plond,plat),rghnss(plond,plat),
     $              evapf (plond,plat),vevapf(plond,plat),
     $              snwjan(plond,plat),snwjly(plond,plat)
C
      real vegtyp, ! Surface thermal type, based on veg type
     $     rghnss, ! Aerodynamic roughness length
     $     evapf , ! Constant surface evaporability
     $     vevapf, ! Variable surface evaporability
     $     snwjan, ! Snow cover (liq water equiv) for January
     $     snwjly  ! Snow cover (liq water equiv) for July
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer lat          ! Lat index for two dimensional data arrays
      integer ioro(plond)  ! Surface type flag (ocean, land, sea ice)
      real sndpth(plond),  ! Snow physical depth
     $     coszrs(plond)   ! Cosine solar zenith angle
C
C Output arguments
C
      real albs(plond),    ! Srf alb for direct rad   0.2-0.7 micro-ms
     $     albl(plond),    ! Srf alb for direct rad   0.7-5.0 micro-ms
     $     albsd(plond),   ! Srf alb for diffuse rad  0.2-0.7 micro-ms
     $     albld(plond)    ! Srf alb for diffuse rad  0.7-5.0 micro-ms
C
C---------------------------Local variables-----------------------------
C
      integer i,ii         ! Longitude indices
      integer indx(plond)  ! Longitude index array (land)
      integer npts         ! Number of land points
      real rs,             ! Empirical fact strng znth angl dependence
     $     rw,             ! Empirical fact weak  znth angl dependence
     $     frsnow,         ! Horizontal fraction of snow cover
     $     snwhgt,         ! Physical snow height
     $     rghsnw          ! Roughness for horizontal snow cover fractn
      real salbs(plond),   ! Snow alb for direct rad  0.2-0.7 micro-ms
     $     salbl(plond),   ! Snow alb for direct rad  0.7-5.0 micro-ms
     $     salbsd(plond),  ! Snow alb for diffuse rad  0.2-0.7 micro-ms
     $     salbld(plond)   ! Snow alb for diffuse rad  0.7-5.0 micro-ms
C
C Externals
C
      external wheneq
C
C-----------------------------------------------------------------------
C
C Find land surfaces
C
      call wheneq(plon,ioro,1,1,indx,npts)
CDIR$ IVDEP
      do ii=1,npts
         i = indx(ii)
         if (coszrs(i).gt.0.0) then         ! Sun above horizon
C
C Use empirical factors to adjust surface albedos for zenith angle
C effects, distinguishing between strong and weakly dependent surfaces:
C
            rs = 1.4/(1. + .8*coszrs(i))
            rw = 1.1/(1. + .2*coszrs(i))
            albs(i)  = albvss(i,lat)*frctst(i,lat)*rs +
     $                 albvsw(i,lat)*(1. - frctst(i,lat))*rw
            albl(i)  = albnis(i,lat)*frctst(i,lat)*rs +
     $                 albniw(i,lat)*(1. - frctst(i,lat))*rw
            albsd(i) = albvss(i,lat)*frctst(i,lat) +
     $                 albvsw(i,lat)*(1. - frctst(i,lat))
            albld(i) = albnis(i,lat)*frctst(i,lat) +
     $                 albniw(i,lat)*(1. - frctst(i,lat))
            if (sndpth(i).gt.0.) then
               salbsd(i) = snws
               salbld(i) = snwl
            end if
         else
C
C Sun below horizon: set land albedos to zero
C
            albs(i) = 0.
            albl(i) = 0.
            albsd(i) = 0.
            albld(i) = 0.
         end if
      end do
CDIR$ IVDEP
      do ii=1,npts
         i = indx(ii)
         if (sndpth(i).gt.0. .and. coszrs(i).gt.0.) then
            if (coszrs(i).lt.0.5) then
C
C Zenith angle regime 1 ( coszrs < 0.5 ).
C Set direct snow albedos (limit to 0.98 max)
C
               salbs(i) = amin1(0.98,salbsd(i) + (1. - salbsd(i))*0.5*
     $                          ((3./(1. + 4.*coszrs(i))) - 1.))
               salbl(i) = amin1(0.98,salbld(i) + (1. - salbld(i))*0.5*
     $                          ((3./(1. + 4.*coszrs(i))) - 1.))
            else
C
C Zenith angle regime 2 ( coszrs >= 0.5 )
C
               salbs(i) = snws
               salbl(i) = snwl
            end if
C
C Compute both diffuse and direct total albedos
C
            snwhgt = sndpth(i)
            rghsnw = amax1(rghnss(i,lat),0.25)
            frsnow = snwhgt/(rghsnw + snwhgt)
            albs(i)  = albs(i) *(1. - frsnow) + salbs(i) *frsnow
            albl(i)  = albl(i) *(1. - frsnow) + salbl(i) *frsnow
            albsd(i) = albsd(i)*(1. - frsnow) + salbsd(i)*frsnow
            albld(i) = albld(i)*(1. - frsnow) + salbld(i)*frsnow
         end if
      end do
C   
      return
      end
      subroutine albocean(lat  ,ioro  ,fnidr  ,sndpth  ,rhos, rs, mu0  ,
     $                    tg   ,hpnd  ,R_pnd  ,hice    ,R_ice   ,
     $                    albs   ,albl    ,albsd  ,albld  )
C-----------------------------------------------------------------------
C Sea ice/ocean shortwave radiation calculation of albedos
C and transmission/absorption in sea ice/ocean.  
C Bruce P. Briegleb  26 October 2006 for snow or pond or ice
C
C Computes surface albedos for direct/diffuse incident radiation for
C two spectral intervals:
C   s = 0.2-0.7 micro-meters
C   l = 0.7-5.0 micro-meters
C-----------------------------------------------------------------------
      implicit none
C------------------------------Parameters-------------------------------
C
C Resolution parameters
C
      integer
     $     plon,      ! number of longitudes
     $     nxpt,      ! no.of points outside active domain for interpolant
     $     plond      ! slt extended domain longitude
      integer ksnow   ! Number of snow layers in CCSM
      integer kseaice ! Number of sea ice layers in CCSM
      integer klev    ! Number of layers minus one (0 layer at top) for rad calculation
      integer klevp   ! klev + 1
      parameter(plon    = 1,
     $          nxpt    = 1,
     $          plond   = plon + 1 + 2*nxpt,
     $          ksnow   = 1, kseaice =  4, klev = ksnow + kseaice + 1,
     $          klevp   = klev + 1)
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer lat           ! Lat index for two dimensional data arrays
      integer ioro(plond)   ! Surface type flag (ocean, land, sea ice)
      real fnidr(plond)     ! fraction of direct to total nir down srf flux
      real sndpth(plond)    ! Snow physical depth (m)
      real rhos(plond)      ! snow density
      real rs(plond)        ! snow grain radius
      real R_ice(plond)     ! sea ice standard deviation tuning parameter
      real R_pnd(plond)     ! ponded ice standard deviation tuning parameter
      real mu0(plond)       ! Cosine solar zenith angle
      real tg(plond)        ! surface skin temperature (K)
      real hpnd(plond)      ! pond thickness (m)
      real hice(plond)      ! sea ice thickness (m)
C
C Output arguments
C
      real albs(plond)      ! Srf alb for direct rad   0.2-0.7 micro-ms
      real albl(plond)      ! Srf alb for direct rad   0.7-5.0 micro-ms
      real albsd(plond)     ! Srf alb for diffuse rad  0.2-0.7 micro-ms
      real albld(plond)     ! Srf alb for diffuse rad  0.7-5.0 micro-ms
C
C---------------------------Local variables-----------------------------
C
      integer i,ii,         ! Longitude indices
     $        indx(plond),  ! Indices for computation points (ocean or sea ice)
     $        indxo(plond), ! Indices for computation points (ocean)
     $        indxsi(plond),! Indices for computation points (sea ice)
     $        npts,         ! Number of ocean/sea ice points
     $        nptso,        ! Number of ocean points
     $        nptssi        ! Number of sea ice points
      integer k             ! level index
C
C Externals
C
      external  wheneq      ! When equal funct gives indices for condtn
      external  whenne      ! When not equal funct gives indices for condtn
C
C-----------------------------------------------------------------------
C 
      real 
     &   dT_mlt    ! change in temp to give melt pond albedo change
     &,  Timelt    ! sea ice melting temperature
     &,  c0        ! 0.0
     &,  c1        ! 1.0

      parameter (
     &   dT_mlt    = -1.      ! change in temp to give melting albedo change
     &,  Timelt    = 0.00     ! ice melting temperature
     &,  c0        = 0.00     ! 0.0
     &,  c1        = 1.00 )   ! 1.0

      real 
     &   hs     ! snow thickness (m)
     &,  fi     ! bare ice area fraction
     &,  fp     ! melt pond fraction
     &,  hp     ! melt pond depth
     &,  z      ! depth in ice (m)
     &,  Klmbda ! irradiance extinction coefficient
      real dz,dz_ssl  ! Sea ice layer thickness and surface layer thickness
C
C-----------------------------------------------------------------------
C Fluxes for sea ice
C
      real hi_ssl     ! sea ice surface scattering layer thickness (m)
      real hs_ssl     ! snow surface scattering layer thickness (m)
      integer ksrf    ! interface index for surface absorption
C
      real Fdirup_vs  ! Up   flux to dir beam at model interface vs band
      real Fdirdn_vs  ! Down flux to dir beam at model interface vs band
      real Fdifup_vs  ! Up   flux to dif beam at model interface vs band
      real Fdifdn_vs  ! Down flux to dif beam at model interface vs band
C
      real Fdirup_ni  ! Up   flux to dir beam at model interface ni band
      real Fdirdn_ni  ! Down flux to dir beam at model interface ni band
      real Fdifup_ni  ! Up   flux to dif beam at model interface ni band
      real Fdifdn_ni  ! Down flux to dif beam at model interface ni band
C
      common/radflux_seaice/
     &              hi_ssl, hs_ssl
     &,             Fdirup_vs(plond,0:klevp),Fdirdn_vs(plond,0:klevp)
     &,             Fdifup_vs(plond,0:klevp),Fdifdn_vs(plond,0:klevp)
     &,             Fdirup_ni(plond,0:klevp),Fdirdn_ni(plond,0:klevp)
     &,             Fdifup_ni(plond,0:klevp),Fdifdn_ni(plond,0:klevp)
     &,             ksrf
C-----------------------------------------------------------------------
C Absorption data for sea ice
C
      real  
     &   I_vs     ! frac transmission vs through sea ice surface layer
     &,  I_ni     ! frac transmission ni through sea ice surface layer
     &,  Tri_vs   ! frac transmission vs, surface to sea ice layer interface 
     &,  Tri_ni   ! frac transmission ni, surface to sea ice layer interface 
     &,  Tro_vs   ! frac transmission vs to ocean
     &,  Tro_ni   ! frac transmission ni to ocean
     &,  zd       ! interface depths for snow/pond and sea ice (from its own surface)
      common/seaice/I_vs,I_ni,zd(0:klevp)
     &             ,Tri_vs(0:klevp),Tri_ni(0:klevp)
     &             ,Tro_vs,Tro_ni
C
C-----------------------------------------------------------------------
C
C Find ocean and sea ice surfaces.
C
      call whenne(plon,ioro,1,1,indx,npts)      ! Ocean or sea ice
      call wheneq(plon,ioro,1,0,indxo,nptso)    ! Open ocean
      call wheneq(plon,ioro,1,2,indxsi,nptssi)  ! Sea ice
C
C Initialize all ocean/sea ice surface albedos to zero
C
      do ii=1,npts
         i = indx(ii)
C
         albs(i)   = 0.
         albl(i)   = 0.
         albsd(i)  = 0.
         albld(i)  = 0.
C
      end do
C
C Start print out input parameters and computed output albedos:
C
      write(6,777) 
 777  format(
     $ /'  .... Begin CCSM3 Sea Ice Delta Eddington calculation ....')
c     
      if( ioro(i) .eq. 0 ) then
        write(6,*) ' surface type = ocean'
      else if( ioro(i) .eq. 2 ) then
c...        write(6,*) ' surface type = sea ice'
      else 
        write(6,*) ' unknown surface type '
      endif
      write(6,*) ' cosine solar zenith angle    = ',mu0(1)
c..      write(6,*) ' surface skin temp C       = ',
c..     $           (tg(1)-273.16)
C
C Skip if no sea ice points
C
      if( nptssi .gt. 0 ) then
C 
C Sea Ice solar radiation
C
      call simcsw(indxsi  ,nptssi  ,fnidr   ,mu0    ,sndpth ,
     $            rhos    ,rs      ,hice    ,R_ice  ,hpnd   ,
     $            R_pnd   ,albs    ,albl    ,albsd  ,albld )
      do ii=1,nptssi
         i = indxsi(ii)
         if (mu0(i).gt.0.) then
           if( sndpth(i) .gt. .00001 ) then
             ksrf = 1
           else
             ksrf = 1 + ksnow + 1
           endif
         endif
      enddo
C
C End of sea ice albedo calculation.
C
      endif 
C
C Ice-free ocean albedos function of solar zenith angle only, and
C independent of spectral interval:
C
      do ii=1,nptso
         i = indxo(ii)
         if (mu0(i).gt.0.) then
            albl(i)  = (.026/(mu0(i)**1.7 + .065)) +
     $                 (.15*(mu0(i) - 0.10)*
     $                      (mu0(i) - 0.50)*
     $                      (mu0(i) - 1.00)  )
            albs(i)  = albl(i)
            albld(i) = 0.06
            albsd(i) = 0.06
         end if
      end do
C
C Print output albedos:
C
c
        write(6,779) 
 779    format('  visible and near-ir direct and diffuse albedos')


        write(6,591) albs(1)
 591    format('  albvs_dir   = ',f8.5)
        write(6,592) albsd(1)
 592    format('  albvs_dif   = ',f8.5)
        write(6,593) albl(1)
 593    format('  albni_dir   = ',f8.5)
        write(6,594) albld(1)
 594    format('  albni_dif   = ',f8.5)
c
        write(6,*) 
     $     ' up/down normalized fluxes vs, ni into snow/sea ice '
        write(6,*) ' Klamda = irradiance extinction coefficent per m'
        write(6,*) 
     $' level   Fdirdn_vs  Fdirup_vs Fdifdn_vs Fdifup_vs  '
        dz_ssl    = hi_ssl
        dz     = hice(1)/real(klevp-ksnow)
        if( dz_ssl .gt. dz/2.0 ) dz_ssl = dz/2.0
        do k=0,klevp
          if( k.lt.klevp ) then
            Klmbda = -2.*(Fdifdn_vs(1,k+1)-Fdifdn_vs(1,k))
     $               /((Fdifdn_vs(1,k+1)+Fdifdn_vs(1,k)+.0001)
     $                *(zd(k+1)-zd(k)+.0001))
          endif
          if( k.lt.klevp) then
              write(6,111) k,Fdirdn_vs(1,k),Fdirup_vs(1,k),
     $                     Fdifdn_vs(1,k),Fdifup_vs(1,k),
     $                     Klmbda
 111          format(2x,i3,4x,2x,4(f8.5,2x)/50x,'Klamda = ',f8.5)
          else
            write(6,112) k,Fdirdn_vs(1,k),Fdirup_vs(1,k),
     $                   Fdifdn_vs(1,k),Fdifup_vs(1,k)
          endif
        enddo

        write(6,*) 
     $' level   Fdirdn_ni  Fdirup_ni Fdifdn_ni Fdifup_ni  '
        do k=0,klevp
          write(6,112) k,Fdirdn_ni(1,k),Fdirup_ni(1,k),
     $                   Fdifdn_ni(1,k),Fdifup_ni(1,k)
 112      format(2x,i3,4x,2x,4(f8.5,2x))
        enddo
c
      write(6,778)
 778  format('  .... CCSM3 Sea Ice Delta Eddington',
     $       ' calculation completed ....')
c     
C
C done with simalb
C
      return
      end
      subroutine simcsw(indxsi  ,nptssi  ,fnidr   ,mu0     ,hs     ,
     $                  rhos    ,rs      ,hi      ,R_ice   ,hp     ,
     $                  R_pnd   ,albs    ,albl    ,albsd   ,albld)
C-----------------------------------------------------------------------
C /home/bruceb/ccsm/ccsm3/snw_ice_albedo/final_1D_Delt_Edd_solar/src
C
C Set up optical property profiles, based on snow, water and sea ice
C IOPs from Briegleb and Light, 2006 Delta-Eddington solar radiation
C treatment for sea ice in the Community Climate System Model (CCSM).
C     3 November 2006 Bruce P. Briegleb
C
C Computes column Delta-Eddington radiation solution for specific
C surface type: either snow over sea ice, bare sea ice, or ponded sea ice.
C
C Divides solar spectrum into 3 intervals: 0.2-0.7, 0.7-1.19, and
C 1.19-5.0 micro-meters. The latter two are added (using an assumed
C partition of incident shortwave) to give the final output of 0.2-0.7 
C visible and 0.7-5.0 infrared albedos and fluxes.
C
C Specifies vertical layer optical properties based on input snow depth, 
C density and grain radius, along with ice and pond depths, then computes 
C layer by layer delta-Eddington reflectivity, transmissivity and combines
C layers (done by calling routine simded). Finally, surface albedos
C and internal fluxes/flux divergences are evaluated.
C
C  Description of the level and layer index conventions. This is 
C  for the standard case of one snow layer and four sea ice layers.
C
C  Please read the following; otherwise, there is 99.9% chance you
C  will be confused about indices at some point in time........ :)
C
C  CCSM3 has four evenly-spaced sea ice layers, and the new CICE4.0
C  snow treatment has one snow layer above the sea ice. This snow
C  layer has finite heat capacity, so that surface absorption must
C  be distinguished from internal. The Delta-Eddington solar radiation
C  needs to add extra surface scattering layers to both snow and sea
C  ice. Note that in the following, we assume a fixed vertical layer 
C  structure for the radiation calculation. In other words, we always
C  have the structure shown below for one snow and four sea ice layers,
C  but for ponded ice have a pond over the sea ice, and for bare ice
C  just treat these top layers over sea ice as transparent air.
C  
C  SSL = surface scattering layer for either snow or sea ice
C  DL  = drained layer for sea ice
C  INT = interior layers for sea ice
C
C  Notice that the radiation level starts with 0 at the top. Thus,
C  the total number radiation layers is klev+1, where klev is the
C  sum of ksnow, the number of CCSM snow layers, and kseaice, the
C  number of CCSM sea ice layers, plus the sea ice SSL: 
C  klev = 1 + ksnow + kseaice. 
C
C  For the standard case illustrated below, ksnow=1, kseaice=4,
C  and klev=6, with the number of layer interfaces klevp=klev+1.
C  Layer interfaces are the surfaces on which reflectivities,
C  transmissivities and fluxes are evaluated.
C
C  CCSM3 Sea Ice Model            Delta-Eddington Solar Radiation
C                                     Layers and Interfaces
C
C                             Layer Index             Interface Index
C    ---------------------            ---------------------  0
C                                  0  \\\   snow SSL    \\\
C       snow layer 1                  ---------------------  1
C                                  1    rest of snow layer
C    +++++++++++++++++++++            +++++++++++++++++++++  2
C                                  2  \\\ sea ice SSL   \\\
C      sea ice layer 1                ---------------------  3
C                                  3      sea ice  DL
C    ---------------------            ---------------------  4
C
C      sea ice layer 2             4      sea ice INT
C
C    ---------------------            ---------------------  5
C
C      sea ice layer 3             5      sea ice INT
C
C    ---------------------            ---------------------  6
C
C      sea ice layer 4             6      sea ice INT
C
C    ---------------------            ---------------------  7
C
C When snow lies over sea ice, the radiation absorbed in the 
C snow SSL is used for surface heating, and that in the rest
C of the snow layer for its internal heating. For sea ice in
C this case, all of the radiant heat absorbed in both the
C sea ice SSL and the DL are used for sea ice layer 1 heating.
C
C When pond lies over sea ice, and for bare sea ice, all of the
C radiant heat absorbed within and above the sea ice SSL is used
C for surface heating, and that absorbed in the sea ice DL is
C used for sea ice layer 1 heating.
C
C This routine in particular specifies the vertical profiles of
C optical properties required for the radiation solution accomplished
C in simded. 
C        call simded(indxsi  ,nptssi  ,mu0     ,srf     ,tau     ,
C     $              w       ,g       ,albodr  ,albodf  ,trndir  ,
C     $              trntdr  ,trndif  ,rupdir  ,rupdif  ,rdndif)
C
C Basically, vertical profiles of tau, w and g are required over 
C the klev+1 layers, where klev = 1 + ksnow + kseaice.
C
C-----------------------------------------------------------------------
      implicit none
C------------------------------Parameters-------------------------------
C
C Resolution parameters
C
      integer plon      ! number of longitudes
      integer nxpt      ! no.of points outside active domain for interpolant
      integer plond     ! slt extended domain longitude
      integer ksnow     ! Number of snow layers in CCSM
      integer kseaice   ! Number of sea ice layers in CCSM
      integer klev      ! Number of layers minus one (0 layer at top) for rad calculation
      integer klevp     ! klev + 1
      parameter(plon    = 1,
     $          nxpt    = 1,
     $          plond   = plon + 1 + 2*nxpt,
     $          ksnow   = 1, kseaice =  4, klev = ksnow + kseaice + 1,
     $          klevp   = klev + 1)
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer indxsi(plond) ! Indices for computation points (sea ice)
      integer nptssi        ! Number of sea ice points
C
      real fnidr(plond)     ! fraction of direct to total nir down srf flux
      real mu0(plond)       ! cosine solar zenith angle
      real hs(plond)        ! Snow physical depth (m)
      real rhos(plond)      ! snow density (kg/m3)
      real rs(plond)        ! snow grain radius (micro-meters)
      real hi(plond)        ! sea ice thickness (m)
      real hp(plond)        ! Melt pond depth (m)
      real R_ice(plond)     ! sea ice standard deviation tuning parameter
      real R_pnd(plond)     ! ponded ice standard deviation tuning parameter
C
C Output arguments
C
      real albs(plond)         ! 0.2-0.7 micro-meter srfc alb to direct rad
      real albl(plond)         ! 0.7-5.0 micro-meter srfc alb to direct rad
      real albsd(plond)        ! 0.2-0.7 micro-meter srfc alb to diffuse rad
      real albld(plond)        ! 0.7-5.0 micro-meter srfc alb to diffuse rad
C
C------------------------------Externals--------------------------------
C
      integer   isrchfgt       ! Search for first array element > 0
      integer   isrchfle       ! Search for first array element < 0
      external  simded         ! Computes delta-eddington solution
      external  isrchfgt       ! Search for first array element > 0
      external  isrchfle       ! Search for first array element < 0
C
C---------------------------Local variables-----------------------------
C
      integer       ns            ! Spectral loop index
      integer       i,ii          ! Longitude loop indices
      integer       k             ! Level loop index
      integer       ktop          ! Total number of levels overlying sea ice minus 1
      integer       kice          ! Level starting index for sea ice
C
      integer srf(plond)          ! Type of layer over ice: 0=air,1=snow,2=water
      real tau(plond,0:klev)      ! Layer extinction optical depth
      real w(plond,0:klev)        ! Layer single scattering albedo
      real g(plond,0:klev)        ! Layer asymmetery paramter
      real dz,dz_ssl  ! Sea ice layer thickness and surface layer thickness
      real fs      ! DL scaling factor to reduce (kseaice<4) or increase (kseaice>4) tau_dl
      real c0      ! 0
      real c1      ! 1
      parameter (
     &   c0        = 0.00     ! 0
     &,  c1        = 1.00 )   ! 1
C
C Spectral properties
C
      integer   nspint            ! Number of solar spectral intervals
      parameter ( nspint = 3 )
      real wt(nspint)             ! Spectral near-ir weights
C
C Optical properties: k in m^{-1}; Q, w and g are dimensionless
C
C Snow optical properties
C
      real Qs(nspint)           ! Snow extinction efficiency
      real ks(nspint)           ! Snow extinction coefficient
      real ws(nspint)           ! Snow single scattering albedo
      real gs(nspint)           ! Snow asymmetry parameter
C
C Sea ice optical properties, distinguished by spectral interval
C mn = mean values
C
      real ki_ssl_mn(nspint)     ! Surface-layer ice extinction coefficient
      real wi_ssl_mn(nspint)     ! Surface-layer ice single scattering albedo
      real gi_ssl_mn(nspint)     ! Surface-layer ice asymmetry parameter
C
      real ki_dl_mn(nspint)      ! Drained-layer ice extinction coefficient
      real wi_dl_mn(nspint)      ! Drained-layer ice single scattering albedo
      real gi_dl_mn(nspint)      ! Drained-layer ice asymmetry parameter
C
      real ki_int_mn(nspint)     ! Interior-layer ice extinction coefficient
      real wi_int_mn(nspint)     ! Interior-layer ice single scattering albedo
      real gi_int_mn(nspint)     ! Interior-layer ice asymmetry parameter
C
      real ki_p_ssl_mn(nspint)   ! Ice under pond srf scat layer extinction coefficient
      real wi_p_ssl_mn(nspint)   ! Ice under pond srf scat layer single scattering albedo
      real gi_p_ssl_mn(nspint)   ! Ice under pond srf scat layer asymmetry parameter
C
      real ki_p_int_mn(nspint)   ! Ice under pond interior extinction coefficient
      real wi_p_int_mn(nspint)   ! Ice under pond interior single scattering albedo
      real gi_p_int_mn(nspint)   ! Ice under pond interior asymmetry parameter
c
      real fp_ice                 ! ice fraction of scat coeff for + stn dev in alb
      real fm_ice                 ! ice fraction of scat coeff for - stn dev in alb
c
      real fp_pnd                 ! ponded ice fraction of scat coeff for + stn dev in alb
      real fm_pnd                 ! ponded ice fraction of scat coeff for - stn dev in alb
C
C actual used in calculation
C
      real ki_ssl(nspint)         ! Surface-layer ice extinction coefficient
      real wi_ssl(nspint)         ! Surface-layer ice single scattering albedo
      real gi_ssl(nspint)         ! Surface-layer ice asymmetry parameter
C
      real ki_dl(nspint)          ! Drained-layer ice extinction coefficient
      real wi_dl(nspint)          ! Drained-layer ice single scattering albedo
      real gi_dl(nspint)          ! Drained-layer ice asymmetry parameter
C
      real ki_int(nspint)         ! Interior-layer ice extinction coefficient
      real wi_int(nspint)         ! Interior-layer ice single scattering albedo
      real gi_int(nspint)         ! Interior-layer ice asymmetry parameter
C
      real kalg                   ! algae absorption coefficient for 0.5 m thick layer
      data kalg / 0.6 /           ! for path of 75 mg Chl a / m2
C
      real ki_p_ssl(nspint)       ! Ice under pond srf scat layer extinction coefficient
      real wi_p_ssl(nspint)       ! Ice under pond srf scat layer single scattering albedo
      real gi_p_ssl(nspint)       ! Ice under pond srf scat layer asymmetry parameter
C
      real ki_p_int(nspint)       ! Ice under pond extinction coefficient
      real wi_p_int(nspint)       ! Ice under pond single scattering albedo
      real gi_p_int(nspint)       ! Ice under pond asymmetry parameter
C
C Water optical properties
C
      real kw(nspint)             ! Water extinction coefficient
      real ww(nspint)             ! Water single scattering albedo
      real gw(nspint)             ! Water asymmetry parameter
C
C Albedos for underlying ocean
C
      real albodr(plond)          ! Spectral ocean albedo to direct rad
      real albodf(plond)          ! Spectral ocean albedo to diffuse rad
c
c For tuning sea ice and pond iops 
c
      real sig,kabs,sigp
      data fp_ice / 0.15 /
      data fm_ice / 0.15 /
c
      data fp_pnd / 2.00 /
      data fm_pnd / 0.50 /
c
c set Qs, ws and gs table input
c
      integer nmbrad             ! number of snow grain radii for table input
      integer nr                 ! loop index for grain number
      parameter( nmbrad = 32 )
      real rs_tab(nmbrad)        ! snow grain radius for each table entry (micro-meters)
      real delr                  ! snow grain radius interpolation parameter
      real Qs_tab(nspint,nmbrad) ! extinction efficiency for each snow grain radius
      real ws_tab(nspint,nmbrad) ! single scatter albedo for each snow grain radius
      real gs_tab(nspint,nmbrad) ! assymetry parameter   for each snow grain radius
c
      data rs_tab/ 
     $                5.,    7.,   10.,   15.,   20.,
     $               30.,   40.,   50.,   65.,   80.,
     $              100.,  120.,  140.,  170.,  200.,
     $              240.,  290.,  350.,  420.,  500.,
     $              570.,  660.,  760.,  870., 1000.,
     $             1100., 1250., 1400., 1600., 1800.,
     $             2000., 2500./
c
      data Qs_tab/
     $    2.131798,    2.187756,    2.267358,
     $    2.104499,    2.148345,    2.236078,
     $    2.081580,    2.116885,    2.175067,
     $    2.062595,    2.088937,    2.130242,
     $    2.051403,    2.072422,    2.106610,
     $    2.039223,    2.055389,    2.080586,
     $    2.032383,    2.045751,    2.066394,
     $    2.027920,    2.039388,    2.057224,
     $    2.023444,    2.033137,    2.048055,
     $    2.020412,    2.028840,    2.041874,
     $    2.017608,    2.024863,    2.036046,
     $    2.015592,    2.022021,    2.031954,
     $    2.014083,    2.019887,    2.028853,
     $    2.012368,    2.017471,    2.025353,
     $    2.011092,    2.015675,    2.022759,
     $    2.009837,    2.013897,    2.020168,
     $    2.008668,    2.012252,    2.017781,
     $    2.007627,    2.010813,    2.015678,
     $    2.006764,    2.009577,    2.013880,
     $    2.006037,    2.008520,    2.012382,
     $    2.005528,    2.007807,    2.011307,
     $    2.005025,    2.007079,    2.010280,
     $    2.004562,    2.006440,    2.009333,
     $    2.004155,    2.005898,    2.008523,
     $    2.003794,    2.005379,    2.007795,
     $    2.003555,    2.005041,    2.007329,
     $    2.003264,    2.004624,    2.006729,
     $    2.003037,    2.004291,    2.006230,
     $    2.002776,    2.003929,    2.005700,
     $    2.002590,    2.003627,    2.005276,
     $    2.002395,    2.003391,    2.004904,
     $    2.002071,    2.002922,    2.004241/
c
      data ws_tab/
     $   0.9999994,  0.9999673,  0.9954589,
     $   0.9999992,  0.9999547,  0.9938576,
     $   0.9999990,  0.9999382,  0.9917989,
     $   0.9999985,  0.9999123,  0.9889724,
     $   0.9999979,  0.9998844,  0.9866190,
     $   0.9999970,  0.9998317,  0.9823021,
     $   0.9999960,  0.9997800,  0.9785269,
     $   0.9999951,  0.9997288,  0.9751601,
     $   0.9999936,  0.9996531,  0.9706974,
     $   0.9999922,  0.9995783,  0.9667577,
     $   0.9999903,  0.9994798,  0.9621007,
     $   0.9999885,  0.9993825,  0.9579541,
     $   0.9999866,  0.9992862,  0.9541924,
     $   0.9999838,  0.9991434,  0.9490959,
     $   0.9999810,  0.9990025,  0.9444940,
     $   0.9999772,  0.9988171,  0.9389141,
     $   0.9999726,  0.9985890,  0.9325819,
     $   0.9999670,  0.9983199,  0.9256405,
     $   0.9999605,  0.9980117,  0.9181533,
     $   0.9999530,  0.9976663,  0.9101540,
     $   0.9999465,  0.9973693,  0.9035031,
     $   0.9999382,  0.9969939,  0.8953134,
     $   0.9999289,  0.9965848,  0.8865789,
     $   0.9999188,  0.9961434,  0.8773350,
     $   0.9999068,  0.9956323,  0.8668233,
     $   0.9998975,  0.9952464,  0.8589990,
     $   0.9998837,  0.9946782,  0.8476493,
     $   0.9998699,  0.9941218,  0.8367318,
     $   0.9998515,  0.9933966,  0.8227881,
     $   0.9998332,  0.9926888,  0.8095131,
     $   0.9998148,  0.9919968,  0.7968620,
     $   0.9997691,  0.9903277,  0.7677887/
c
      data gs_tab /
     $    0.859913,    0.848003,    0.824415,
     $    0.867130,    0.858150,    0.848445,
     $    0.873381,    0.867221,    0.861714,
     $    0.878368,    0.874879,    0.874036,
     $    0.881462,    0.879661,    0.881299,
     $    0.884361,    0.883903,    0.890184,
     $    0.885937,    0.886256,    0.895393,
     $    0.886931,    0.887769,    0.899072,
     $    0.887894,    0.889255,    0.903285,
     $    0.888515,    0.890236,    0.906588,
     $    0.889073,    0.891127,    0.910152,
     $    0.889452,    0.891750,    0.913100,
     $    0.889730,    0.892213,    0.915621,
     $    0.890026,    0.892723,    0.918831,
     $    0.890238,    0.893099,    0.921540,
     $    0.890441,    0.893474,    0.924581,
     $    0.890618,    0.893816,    0.927701,
     $    0.890762,    0.894123,    0.930737,
     $    0.890881,    0.894397,    0.933568,
     $    0.890975,    0.894645,    0.936148,
     $    0.891035,    0.894822,    0.937989,
     $    0.891097,    0.895020,    0.939949,
     $    0.891147,    0.895212,    0.941727,
     $    0.891189,    0.895399,    0.943339,
     $    0.891225,    0.895601,    0.944915,
     $    0.891248,    0.895745,    0.945950,
     $    0.891277,    0.895951,    0.947288,
     $    0.891299,    0.896142,    0.948438,
     $    0.891323,    0.896388,    0.949762,
     $    0.891340,    0.896623,    0.950916,
     $    0.891356,    0.896851,    0.951945,
     $    0.891386,    0.897399,    0.954156/
c
      real rhoi                   ! pure ice density (kg/m3)
      data rhoi /  917.0 /
C
C ice, ponded ice and water IOPs
C
      data ki_ssl_mn /  1000.1, 1003.7, 7042. /
      data wi_ssl_mn /  .9999,  .9963,  .9088 /
      data gi_ssl_mn /   .94,     .94,    .94 /
C
      data ki_dl_mn  /   100.2, 107.7,  1309. /
      data wi_dl_mn  /   .9980,  .9287, .0305 /
      data gi_dl_mn  /   .94,     .94,    .94 /
C
      data ki_int_mn /    20.2,  27.7,  1445. /
      data wi_int_mn /   .9901, .7223,  .0277 /
      data gi_int_mn /   .94,    .94,     .94 /
C
      data ki_p_ssl_mn /  70.2,  77.7,  1309. /
      data wi_p_ssl_mn / .9972, .9009,  .0305 /
      data gi_p_ssl_mn / .94,   .94,    .94   /
C
      data ki_p_int_mn /  20.2,  27.7,  1445. /
      data wi_p_int_mn / .9901, .7223,  .0277 /
      data gi_p_int_mn / .94,   .94,    .94   /
C
      data kw   /    0.20,   12.0,   729. /
      data ww   /    0.0,     0.0,   0.0  /
      data gw   /    0.0,     0.0,   0.0  /
C
C These arrays are defined at model interfaces; 0 is the top of the
C layer above the sea ice; klevp is the sea ice/ocean interface:
C
      real trndir(plond,0:klevp)  ! Solar beam down transm from top
      real trntdr(plond,0:klevp)  ! Total transmission to direct beam for layers above
      real trndif(plond,0:klevp)  ! Diffuse transmission to diffuse beam for layers above
      real rupdir(plond,0:klevp)  ! Ref to dir rad for layers below
      real rupdif(plond,0:klevp)  ! Ref to dif rad for layers below
      real rdndif(plond,0:klevp)  ! Ref to dif rad for layers above
      real refk                   ! Interface multiple scattering k
C
      real fdirup(plond,0:klevp)  ! Up   flux to direct beam at model interface 
      real fdirdn(plond,0:klevp)  ! Down flux to direct beam at model interface
      real fdifup(plond,0:klevp)  ! Up   flux to diffuse beam at model interface 
      real fdifdn(plond,0:klevp)  ! Down flux to diffuse beam at model interface
C
C-----------------------------------------------------------------------
C Fluxes for sea ice
C
      real hi_ssl     ! sea ice surface scattering layer thickness (m)
      real hs_ssl     ! snow surface scattering layer thickness (m)
      integer ksrf    ! interface index for surface absorption
C
      real Fdirup_vs  ! Up   flux to dir beam at model interface vs band
      real Fdirdn_vs  ! Down flux to dir beam at model interface vs band
      real Fdifup_vs  ! Up   flux to dif beam at model interface vs band
      real Fdifdn_vs  ! Down flux to dif beam at model interface vs band
C
      real Fdirup_ni  ! Up   flux to dir beam at model interface ni band
      real Fdirdn_ni  ! Down flux to dir beam at model interface ni band
      real Fdifup_ni  ! Up   flux to dif beam at model interface ni band
      real Fdifdn_ni  ! Down flux to dif beam at model interface ni band
C
      common/radflux_seaice/
     &              hi_ssl, hs_ssl
     &,             Fdirup_vs(plond,0:klevp),Fdirdn_vs(plond,0:klevp)
     &,             Fdifup_vs(plond,0:klevp),Fdifdn_vs(plond,0:klevp)
     &,             Fdirup_ni(plond,0:klevp),Fdirdn_ni(plond,0:klevp)
     &,             Fdifup_ni(plond,0:klevp),Fdifdn_ni(plond,0:klevp)
     &,             ksrf
C
      data hi_ssl / .05 /     ! sea ice surface scattering layer thickness (m)
      data hs_ssl / .04 /     ! snow surface scattering layer thickness (m)
C
C-----------------------------------------------------------------------
C Absorption data for sea ice
C
      real
     &   I_vs     ! frac transmission vs through sea ice surface layer
     &,  I_ni     ! frac transmission ni through sea ice surface layer
     &,  Tri_vs   ! frac transmission vs, surface to sea ice layer interface
     &,  Tri_ni   ! frac transmission ni, surface to sea ice layer interface
     &,  Tro_vs   ! frac transmission vs to ocean
     &,  Tro_ni   ! frac transmission ni to ocean
     &,  zd       ! interface depths for snow/pond and sea ice (from its own surface)
      common/seaice/I_vs,I_ni,zd(0:klevp)
     &             ,Tri_vs(0:klevp),Tri_ni(0:klevp)
     &             ,Tro_vs,Tro_ni
C-----------------------------------------------------------------------
C Data for melt ponds transition to bare sea ice
C
      real hpmin    ! minimum allowed melt pond depth
      real hp0      ! melt pond depth below which iops wghted bare ice/ponded ice
      real sig_i    ! ice scattering coefficient
      real sig_p    ! pond scattering coefficient
      real kext     ! weighted extinction coefficient
C
      data hpmin / .005 /
      data hp0   / .200 /
C-----------------------------------------------------------------------
C if below min pond depth, set pond depth to 0
      do ii=1,nptssi
         i = indxsi(ii)
         if( mu0(i).gt.0.) then
           if( hp(i) .lt. hpmin ) then
             hp(i) = 0.
           endif
         endif
      enddo
C-----------------------------------------------------------------------
C Set surface type for layer over sea ice, based on snow and melt pond 
C depths; set initial albedos, fluxes to zero and set spectral weights:
C
      do ii=1,nptssi
         i = indxsi(ii)
         if( mu0(i).gt.0.) then
           if( hs(i) .eq. 0.0 .and. hp(i) .eq. 0.0 ) then
C air
             srf(i)  = 0 
           else if( hs(i) .gt. 0.0 .and. hp(i) .eq. 0.0 ) then
C snow
             srf(i)  = 1
           else if( hs(i) .eq. 0.0 .and. hp(i) .gt. 0.0 ) then
C melt pond
             srf(i)  = 2
           else
C unknown
             write(6,5) hs(i),hp(i)
 5           format(/' simcsw  unknown combination hs and hp '/
     $               ' hs = ',f8.5,' hp = ',f8.5)
             stop 
           endif
C initialize
           albs(i)  = 0.
           albl(i)  = 0.
           albsd(i) = 0.
           albld(i) = 0.
C
           do k=0,klevp
             Fdirup_vs(i,k)  = 0.0
             Fdirdn_vs(i,k)  = 0.0
             Fdifup_vs(i,k)  = 0.0
             Fdifdn_vs(i,k)  = 0.0
             Fdirup_ni(i,k)  = 0.0
             Fdirdn_ni(i,k)  = 0.0
             Fdifup_ni(i,k)  = 0.0
             Fdifdn_ni(i,k)  = 0.0
           enddo
C
           wt(1) = 1.0
           wt(2) = .67 + (.78-.67)*(1.-fnidr(i))
           wt(3) = .33 + (.22-.33)*(1.-fnidr(i))
         endif
      enddo
C
C set sea ice iops with tuning
C
      do ii=1,nptssi
        i = indxsi(ii)
        if (mu0(i).gt.0.) then
         if( R_ice(i) .ge. 0.0 ) then
          do ns=1,nspint
            kabs = ki_ssl_mn(ns)*(1.0-wi_ssl_mn(ns))
            sig  = ki_ssl_mn(ns)*wi_ssl_mn(ns)
            sigp = sig*(1.0+fp_ice*R_ice(i))
            ki_ssl(ns) = sigp+kabs
            wi_ssl(ns) = sigp/ki_ssl(ns)

            kabs = ki_dl_mn(ns)*(1.0-wi_dl_mn(ns))
            sig  = ki_dl_mn(ns)*wi_dl_mn(ns)
            sigp = sig*(1.0+fp_ice*R_ice(i))
            ki_dl(ns) = sigp+kabs
            wi_dl(ns) = sigp/ki_dl(ns)

            kabs = ki_int_mn(ns)*(1.0-wi_int_mn(ns))
            sig  = ki_int_mn(ns)*wi_int_mn(ns)
            sigp = sig*(1.0+fp_ice*R_ice(i))
            ki_int(ns) = sigp+kabs
            wi_int(ns) = sigp/ki_int(ns)

            gi_ssl(ns) = gi_ssl_mn(ns)
            gi_dl(ns)  = gi_dl_mn(ns)
            gi_int(ns) = gi_int_mn(ns)
          enddo
         else
          do ns=1,nspint
            kabs = ki_ssl_mn(ns)*(1.0-wi_ssl_mn(ns))
            sig  = ki_ssl_mn(ns)*wi_ssl_mn(ns)
            sigp = sig*(1.0+fm_ice*R_ice(i))
            if( sigp .lt. 0. ) sigp = 0.0
            ki_ssl(ns) = sigp+kabs
            wi_ssl(ns) = sigp/ki_ssl(ns)

            kabs = ki_dl_mn(ns)*(1.0-wi_dl_mn(ns))
            sig  = ki_dl_mn(ns)*wi_dl_mn(ns)
            sigp = sig*(1.0+fm_ice*R_ice(i))
            if( sigp .lt. 0. ) sigp = 0.0
            ki_dl(ns) = sigp+kabs
            wi_dl(ns) = sigp/ki_dl(ns)

            kabs = ki_int_mn(ns)*(1.0-wi_int_mn(ns))
            sig  = ki_int_mn(ns)*wi_int_mn(ns)
            sigp = sig*(1.0+fm_ice*R_ice(i))
            if( sigp .lt. 0. ) sigp = 0.0
            ki_int(ns) = sigp+kabs
            wi_int(ns) = sigp/ki_int(ns)

            gi_ssl(ns) = gi_ssl_mn(ns)
            gi_dl(ns)  = gi_dl_mn(ns)
            gi_int(ns) = gi_int_mn(ns)
          enddo
         endif
        endif
      enddo
C
C set ponded sea ice iops with tuning
C
      do ii=1,nptssi
        i = indxsi(ii)
        if (mu0(i).gt.0.) then
         if( R_pnd(i) .ge. 0.0 ) then
          do ns=1,nspint
            kabs = ki_p_ssl_mn(ns)*(1.0-wi_p_ssl_mn(ns))
            sig  = ki_p_ssl_mn(ns)*wi_p_ssl_mn(ns)
            sigp = sig*(1.0+fp_pnd*R_pnd(i))
            ki_p_ssl(ns) = sigp+kabs
            wi_p_ssl(ns) = sigp/ki_p_ssl(ns)
            gi_p_ssl(ns) = gi_p_ssl_mn(ns)

            kabs = ki_p_int_mn(ns)*(1.0-wi_p_int_mn(ns))
            sig  = ki_p_int_mn(ns)*wi_p_int_mn(ns)
            sigp = sig*(1.0+fp_pnd*R_pnd(i))
            ki_p_int(ns) = sigp+kabs
            wi_p_int(ns) = sigp/ki_p_int(ns)
            gi_p_int(ns) = gi_p_int_mn(ns)
          enddo
         else
          do ns=1,nspint
            kabs = ki_p_ssl_mn(ns)*(1.0-wi_p_ssl_mn(ns))
            sig  = ki_p_ssl_mn(ns)*wi_p_ssl_mn(ns)
            sigp = sig*(1.0+fm_pnd*R_pnd(i))
            if( sigp .lt. 0. ) sigp = 0.0
            ki_p_ssl(ns) = sigp+kabs
            wi_p_ssl(ns) = sigp/ki_p_ssl(ns)
            gi_p_ssl(ns) = gi_p_ssl_mn(ns)

            kabs = ki_p_int_mn(ns)*(1.0-wi_p_int_mn(ns))
            sig  = ki_p_int_mn(ns)*wi_p_int_mn(ns)
            sigp = sig*(1.0+fm_pnd*R_pnd(i))
            if( sigp .lt. 0. ) sigp = 0.0
            ki_p_int(ns) = sigp+kabs
            wi_p_int(ns) = sigp/ki_p_int(ns)
            gi_p_int(ns) = gi_p_int_mn(ns)
          enddo
         endif
        endif
      enddo
C
C Begin spectral loop
C
      do 100 ns=1,nspint
C
C Set optical properties of air/snow/pond/ice
C
C top layer
C
        ktop = ksnow
C       sea ice 
        do ii=1,nptssi
           i = indxsi(ii)
C          sun above horizon
           if (mu0(i).gt.0.) then
C            layers over sea ice
C air
             if( srf(i) .eq. 0 ) then
               zd(0)    = 0.0
               do k=0,ktop
                 tau(i,k) = 0.0
                 w(i,k)   = 0.0
                 g(i,k)   = 0.0
                 zd(k+1)  = 0.0
               enddo
C snow
             else if( srf(i) .eq. 1 ) then
               dz_ssl    = hs_ssl
               dz     = hs(i)/real(ksnow)
               if( dz_ssl .gt. dz/2.0 ) dz_ssl = dz/2.0
               zd(0)  = 0.0
               do k=0,ktop
c
c find snow iops using input snow density and snow grain radius:
c
                 if( rs(i) .lt. rs_tab(1) ) then
                   Qs(ns) = Qs_tab(ns,1)
                   ws(ns) = ws_tab(ns,1)
                   gs(ns) = gs_tab(ns,1)
                 else if( rs(i) .ge. rs_tab(nmbrad) ) then
                   Qs(ns) = Qs_tab(ns,nmbrad)
                   ws(ns) = ws_tab(ns,nmbrad)
                   gs(ns) = gs_tab(ns,nmbrad)
                 else
                   do nr=2,nmbrad
                     if( rs_tab(nr-1) .le. rs(i) .and.
     $                   rs(i) .lt. rs_tab(nr)) then
                           delr = (rs(i) - rs_tab(nr-1)) /
     $                            (rs_tab(nr) - rs_tab(nr-1))
                           Qs(ns) = Qs_tab(ns,nr-1)*(1.-delr) + 
     $                              Qs_tab(ns,nr)*delr
                           ws(ns) = ws_tab(ns,nr-1)*(1.-delr) + 
     $                              ws_tab(ns,nr)*delr
                           gs(ns) = gs_tab(ns,nr-1)*(1.-delr) + 
     $                              gs_tab(ns,nr)*delr
                     endif
                   enddo
                 endif
c
                 ks(ns)   = Qs(ns)*((rhos(i)/rhoi)*3./
     $                              (4.*rs(i)*1.0e-6))
                 if( k.eq.0 ) then
                   tau(i,k) = ks(ns)*dz_ssl
                   w(i,k)   = ws(ns)
                   g(i,k)   = gs(ns)
                   zd(k+1)  = dz_ssl
                 else if( k.eq.1 ) then
                   tau(i,k) = ks(ns)*(dz-dz_ssl)
                   w(i,k)   = ws(ns)
                   g(i,k)   = gs(ns)
                   zd(k+1)  = zd(k) + dz-dz_ssl
                 else if( k.ge.2 ) then
                   tau(i,k) = ks(ns)*dz
                   w(i,k)   = ws(ns)
                   g(i,k)   = gs(ns)
                   zd(k+1)  = zd(k) + dz
                 endif 
               enddo
C pond
             else if( srf(i) .eq. 2 ) then
               dz = hp(i)/real(ktop+1)
               zd(0)    = 0.0
               do k=0,ktop
                 tau(i,k) = kw(ns)*dz
                 w(i,k)   = ww(ns)
                 g(i,k)   = gw(ns)
                 zd(k+1)  = zd(k) + dz
               enddo
             endif
C            end sun above horizon
           endif
C          end sea ice
        enddo
C
C sea ice layers; start at index kice
C
        kice = ktop + 1
        do ii=1,nptssi
           i = indxsi(ii)
           if (mu0(i).gt.0.) then
             dz_ssl = hi_ssl
             dz     = hi(i)/real(kseaice)
c
c make dz_ssl for thin ice hi/30
c
             if( hi(i) .le. 1.5 ) dz_ssl = hi(i)/30.
             if( dz_ssl .gt. dz/2.0 ) dz_ssl = dz/2.0
c
             do k=kice,klev
               if( k.eq.kice ) then
                 tau(i,k) = ki_ssl(ns)*dz_ssl
                 w(i,k)   = wi_ssl(ns)
                 g(i,k)   = gi_ssl(ns)
                 zd(k+1)  = dz_ssl
               else if( k.eq.kice+1 ) then
c scale dz for dl relative to 4 even-layer 1.5m case  
                 fs = real(kseaice)/real(4)
                 tau(i,k) = ki_dl(ns)*(dz-dz_ssl)*fs
                 w(i,k)   = wi_dl(ns)
                 g(i,k)   = gi_dl(ns)
                 zd(k+1)  = zd(k) + dz-dz_ssl
               else if( k.ge.kice+2 .and. k .lt. klev ) then
                 tau(i,k) = ki_int(ns)*dz
                 w(i,k)   = wi_int(ns)
                 g(i,k)   = gi_int(ns)
                 zd(k+1)  = zd(k) + dz
               else if( k.eq.klev ) then
c add algae to lowest sea ice layer, visible only:
                 kabs = ki_int(ns)*(1.0-wi_int(ns))
                 if( ns .eq. 1 ) then
                   kabs = kabs + kalg*(0.50/dz)
                 endif
                 sig  = ki_int(ns)*wi_int(ns)
                 tau(i,k) = (kabs+sig)*dz
                 w(i,k)   = (sig/(sig+kabs))
                 g(i,k)   = gi_int(ns)
                 zd(k+1)  = zd(k) + dz
               endif 
             enddo  
C sea ice layers over melt-ponds: set ice optical properties explicitly
             if( srf(i) .eq. 2 ) then
               do k=kice,klev
                 if( k.eq.kice ) then
                   tau(i,k) = ki_p_ssl(ns)*dz_ssl
                   w(i,k)   = wi_p_ssl(ns)
                   g(i,k)   = gi_p_ssl(ns)
                 else if( k.eq.kice+1 ) then
                   tau(i,k) = ki_p_int(ns)*(dz-dz_ssl)
                   w(i,k)   = wi_p_int(ns)
                   g(i,k)   = gi_p_int(ns)
                 else if( k.gt.kice+1 ) then
                   tau(i,k) = ki_p_int(ns)*dz
                   w(i,k)   = wi_p_int(ns)
                   g(i,k)   = gi_p_int(ns)
                 endif 
               enddo
c
c adjust pond iops if pond depth within specified range 
c
               if( hpmin .le. hp(i) .and. hp(i) .le. hp0 ) then
                 do k=kice,klev
                   if( k.eq.kice ) then
                     sig_i = ki_ssl(ns)*wi_ssl(ns)
                     sig_p = ki_p_ssl(ns)*wi_p_ssl(ns)
                     sig   = sig_i + (sig_p-sig_i)*(hp(i)/hp0)
                     kext  = sig + ki_p_ssl(ns)*(1.0-wi_p_ssl(ns))
                   else if( k.eq.kice+1 ) then
                     sig_i = ki_dl(ns)*wi_dl(ns)
                     sig_p = ki_p_int(ns)*wi_p_int(ns)
                     sig   = sig_i + (sig_p-sig_i)*(hp(i)/hp0)
                     kext  = sig + ki_p_int(ns)*(1.0-wi_p_int(ns))
                   else if( k.ge.kice+2 ) then
                     sig_i = ki_int(ns)*wi_int(ns)
                     sig_p = ki_p_int(ns)*wi_p_int(ns)
                     sig   = sig_i + (sig_p-sig_i)*(hp(i)/hp0)
                     kext  = sig + ki_p_int(ns)*(1.0-wi_p_int(ns))
                   endif
                   if( k.eq.kice ) then
                     tau(i,k) = kext*dz_ssl
                   else if( k.eq.kice+1 ) then
                     tau(i,k) = kext*(dz-dz_ssl)
                   else if( k.gt.kice+1 ) then
                     tau(i,k) = kext*dz
                   endif
                   w(i,k)   = sig/kext
                   g(i,k)   = gi_p_int(ns)
                 enddo
               endif
             endif
C
C          end sun above horizon
           endif
C       end sea ice points
        enddo
C
C Set reflectivities for ocean underlying sea ice
C
        if(ns .eq. 1) then
           do ii=1,nptssi
              i = indxsi(ii)
              if (mu0(i).gt.0.) then
                 albodr(i) = .01
                 albodf(i) = .01
              endif
           enddo
        else if(ns .eq. 2) then
           do ii=1,nptssi
              i = indxsi(ii)
              if (mu0(i).gt.0.) then
                 albodr(i) = .0
                 albodf(i) = .0
              endif
           enddo
        else if(ns .eq. 3) then
           do ii=1,nptssi
              i = indxsi(ii)
              if (mu0(i).gt.0.) then
                 albodr(i) = .0
                 albodf(i) = .0
              endif
           enddo
        endif
C
C Layer input properties now completely specified; compute the
C Delta-Eddington solution reflectivities and transmissivities
C for each layer, starting from the top and working downwards;
C then, add the layers going downwards accounting for multiple
C scattering between layers, then starting from the surface and
C adding successive layers upwards to the top:
C
        call simded(indxsi  ,nptssi  ,mu0     ,srf     ,tau     ,
     $              w       ,g       ,albodr  ,albodf  ,trndir  ,
     $              trntdr  ,trndif  ,rupdir  ,rupdif  ,rdndif)
C
C Compute up and down fluxes for each interface, using the 
C added layer properties at each interface:
C
C              layers       interface
C
C       ---------------------  k
C                 k
C       --------------------- 
C
        do k=0,klevp
           do ii=1,nptssi
             i = indxsi(ii)
             if (mu0(i).gt.0.) then
C
C interface scattering
C
               refk = 1./(1. - rdndif(i,k)*rupdif(i,k))
C
C dir tran ref from below times interface scattering, plus diff
C tran and ref from below times interface scattering
C
               fdirup(i,k) = (trndir(i,k)*rupdir(i,k) +
     $                (trntdr(i,k)-trndir(i,k))*rupdif(i,k))*refk
C
C dir tran plus total diff trans times interface scattering plus
C dir tran with up dir ref and down dif ref times interface scattering 
C
               fdirdn(i,k) = trndir(i,k) + (trntdr(i,k) - trndir(i,k) +
     $                trndir(i,k)*rupdir(i,k)*rdndif(i,k))*refk
C
C diffuse tran ref from below times interface scattering
C
               fdifup(i,k) = trndif(i,k)*rupdif(i,k)*refk
C
C diffuse tran times interface scattering
C
               fdifdn(i,k) = trndif(i,k)*refk
C
             endif
           enddo
        enddo
C
C Set surface albedos and fluxes
C
        if(ns .eq. 1) then
           do ii=1,nptssi
              i = indxsi(ii)
              if (mu0(i).gt.0.) then
                 albs(i)  = rupdir(i,0)
                 albsd(i) = rupdif(i,0)
                 do k=0,klevp
                   Fdirup_vs(i,k)  = fdirup(i,k)
                   Fdirdn_vs(i,k)  = fdirdn(i,k)
                   Fdifup_vs(i,k)  = fdifup(i,k)
                   Fdifdn_vs(i,k)  = fdifdn(i,k)
                 enddo
              endif
           enddo
        else if(ns .gt. 1) then
           do ii=1,nptssi
              i = indxsi(ii)
              if (mu0(i).gt.0.) then
                 albl(i)  = albl(i)  + rupdir(i,0)*wt(ns)
                 albld(i) = albld(i) + rupdif(i,0)*wt(ns)
                 do k=0,klevp
                   Fdirup_ni(i,k)  = Fdirup_ni(i,k) + fdirup(i,k)*wt(ns)
                   Fdirdn_ni(i,k)  = Fdirdn_ni(i,k) + fdirdn(i,k)*wt(ns)
                   Fdifup_ni(i,k)  = Fdifup_ni(i,k) + fdifup(i,k)*wt(ns)
                   Fdifdn_ni(i,k)  = Fdifdn_ni(i,k) + fdifdn(i,k)*wt(ns)
                 enddo
              endif
           enddo
        endif
C
C special print of the spectral albedos:
C
        write(6,4321) ns,rupdir(1,0),rupdif(1,0)
 4321   format('  spectral albedo for interval = ',i3,
     $         ' direct and diffuse = ',2(f7.5,1x))
C
  100 continue              ! End of spectral interval loop
C
C done with simcsw
C
      return
      end
      subroutine simded(indxsi  ,nptssi  ,mu0_in  ,srf_in  ,tau_in  ,
     $                  w_in    ,g_in    ,albodr  ,albodf  ,trndir  ,
     $                  trntdr  ,trndif  ,rupdir  ,rupdif  ,rdndif)
C-----------------------------------------------------------------------
C /home/bruceb/ccsm/ccsm3/snw_ice_albedo/final_1D_Delt_Edd_solar/src
C
C Delta-Eddington solution for snow/air/pond over sea ice
C 26 October 2006   Bruce P. Briegleb
C
C Generic solution for a snow/air/pond input column of klev+1 layers,
C with srf_in determining at what interface fresnel refraction occurs.
C
C Computes layer reflectivities and transmissivities, from the top down
C to the lowest interface using the Delta-Eddington solutions for each 
C layer; adds layers from top down to lowest interface, and from the
C lowest interface (underlying ocean) up to the top of the column.
C
C Note that the diffuse reflectivity and transmissivity are computed
C by integrating the direct over several gaussian angles. This is
C because the diffuse reflectivity expression sometimes is negative.
C
C Assumes monochromatic or spectrally uniform properties across a band
C for the input optical parameters.
C
C If total transmission to the interface above a particular layer is
C less than trmin, then no further Delta-Eddington solutions are
C evaluated for layers below.
C
C The following describes how refraction is handled in the calculation.
C
C First, we assume that radiation is refracted when entering either 
C sea ice at the base of the surface scattering layer, or water (i.e. melt 
C pond) from either air or snow; we assume that radiation does not refract 
C when entering snow, nor upon entering sea ice from a melt pond, nor 
C upon entering the underlying ocean from sea ice. 
C
C To handle refraction, we define a "fresnel" layer, which physically 
C is of neglible vertical extent and is non-absorbing, which can be
C added to any sea ice layer or top of melt pond. 
C
C A fresnel layer is added to the top of a melt pond or to the surface
C scattering layer of sea ice if no melt pond lies over it. The fresnel 
C layer accounts for refraction of direct beam and associated reflection
C and transmission for solar radiation. It is assumed that radiation 
C crossing the air/snow, melt pond/sea ice and sea ice/ocean interfaces 
C is not refracted.
C
C Some caution must be exercised for the fresnel layer, because any layer
C to which it is added is no longer a homogeneous layer, as are all other
C individual layers. For all other layers for example, the direct and diffuse 
C reflectivities/transmissivities (R/T) are the same for radiation above or 
C below the layer. This is the meaning of homogeneous! But for the fresnel 
C layer this is not so. Thus, the R/T for this layer must be distinguished 
C for radiation above from that from radiation below. For generality, we 
C treat all layers to be added as inhomogeneous.
C
C-----------------------------------------------------------------------
      implicit none
C------------------------------Parameters-------------------------------
C
C Resolution parameters
C
      integer plon          ! number of longitudes
      integer nxpt          ! no.of points outside active domain for interpolant
      integer plond         ! slt extended domain longitude
      integer ksnow         ! Number of snow layers in CCSM
      integer kseaice       ! Number of sea ice layers in CCSM
      integer klev          ! Number of layers minus one (0 layer at top) for rad calculation
      integer klevp         ! klev + 1
      parameter(plon    = 1,
     $          nxpt    = 1,
     $          plond   = plon + 1 + 2*nxpt,
     $          ksnow   = 1, kseaice =  4, klev = ksnow + kseaice + 1,
     $          klevp   = klev + 1)
C
C------------------------------Input Arguments--------------------------
C
      integer indxsi(plond)       ! Indices for computation points (sea ice)
      integer nptssi              ! Number of sea ice points
      real    mu0_in(plond)       ! Cosine zenith angle
      integer srf_in(plond)       ! Type of layer over ice: 0=air,1=snow,2=water
      real tau_in(plond,0:klev)   ! Layer extinction optical depth
      real w_in(plond,0:klev)     ! Layer single scattering albedo
      real g_in(plond,0:klev)     ! Layer asymmetry parameter
      real albodr(plond)          ! Ocean albedo to direct rad
      real albodf(plond)          ! Ocean albedo to diffuse rad
C
C------------------------------Output Arguments-------------------------
C
C  Note that the following variables are defined on interfaces, with the
C  index k referring to the top interface of the kth layer. For example, 
C  trntdr(k=4) refers to the total transmission to direct radiation to the 
C  top interface of the 4th layer. klevp refers to the lowest sea ice layer 
C  interface with the ocean.
C
      real trndir(plond,0:klevp)  ! Solar beam down transm from top
      real trntdr(plond,0:klevp)  ! Total transmission to direct beam for layers above
      real trndif(plond,0:klevp)  ! Diffuse transmission to diffuse beam for layers above
      real rupdir(plond,0:klevp)  ! Ref to dir rad for layers below
      real rupdif(plond,0:klevp)  ! Ref to dif rad for layers below
      real rdndif(plond,0:klevp)  ! Ref to dif rad for layers above
C
C------------------------------Externals--------------------------------
C
      external  resetr      ! Resets array elements to zero
      external  whenfgt     ! Collect indices greater than condition
C
C---------------------------Local variables-----------------------------
C
      integer    kfrsnl     ! vertical index of layer whose top is Fresnel layer
C
C Following variables are defined for each layer; 0 refers to layer above 
C sea ice. In general we must distinguish directions above and below in 
C the diffuse reflectivity and transmissivity, as layers are not assumed
C to be homogeneous (apart from the single layer delta-Edd solutions); 
C the direct is always from above.
C
      real rdir(plond,0:klev)     ! Layer reflectivity to direct rad
      real rdif_a(plond,0:klev)   ! Layer reflectivity to diffuse rad from above
      real rdif_b(plond,0:klev)   ! Layer reflectivity to diffuse rad from below
      real tdir(plond,0:klev)     ! Layer transmission to direct rad
      real tdif_a(plond,0:klev)   ! Layer transmission to diffuse rad from above
      real tdif_b(plond,0:klev)   ! Layer transmission to diffuse rad from below
      real trnlay(plond,0:klev)   ! Solar beam transm for layer
C
C Minimum total transmission below which no layer computation are done:
C
      real trmin            ! Minimum total transmission allowed
      parameter (trmin = 1.e-3)
C
      integer i             ! Longitude index
      integer k             ! Level index
      integer ii            ! Longitude index
      integer nval          ! Number of long values satisfying criteria
      integer index(plond)  ! Array of longitude indices
      real tautot           ! Layer optical depth
      real wtot             ! Layer single scatter albedo
      real gtot             ! Layer asymmetry parameter
      real ftot             ! Layer forward scatter fraction
      real ts               ! Layer scaled extinction optical depth
      real ws               ! Layer scaled single scattering albedo
      real gs               ! Layer scaled asymmetry parameter
      real rintfc           ! Reflection (multiple) at an interface
      real refkp1           ! Interface multiple scattering for k+1
      real refkm1           ! Interface multiple scattering for k-1
      real tdrrdir          ! Direct tran times layer direct ref 
      real tdndif           ! Total down diffuse = tot tran - direct tran
      real R1               ! Perpendicular polarization reflection amplitude
      real R2               ! Parallel polarization reflection amplitude
      real T1               ! Perpendicular polarization transmission amplitude
      real T2               ! Parallel polarization transmission amplitude
      real Rf_dir_a         ! Fresnel reflection to direct radiation
      real Tf_dir_a         ! Fresnel transmission to direct radiation
      real Rf_dif_a         ! Fresnel reflection to diff radiation from above
      real Rf_dif_b         ! Fresnel reflection to diff radiation from below
      real Tf_dif_a         ! Fresnel transmission to diff radiation from above
      real Tf_dif_b         ! Fresnel transmission to diff radiation from below
      real refindx          ! Refractive index of ice (used for water also)
      parameter (refindx  = 1.31)
C
C----------------Statement functions and other local variables----------
C
C Statement functions and other local variables
C
      real alpha            ! Term in direct reflect and transmissivity
      real gamma            ! Term in direct reflect and transmissivity
      real el               ! Term in alpha,gamma,n,u
      real taus             ! Scaled extinction optical depth
      real omgs             ! Scaled single particle scattering albedo
      real asys             ! Scaled asymmetry parameter
      real u                ! Term in diffuse reflect and transmissivity
      real n                ! Term in diffuse reflect and transmissivity
      real lm               ! Temporary for el
      real mu0              ! cosine solar zenith angle incident
      real mu0n             ! cosine solar zenith angle in medium
      real mu               ! cosine solar zenith for either snow or water
      real ne               ! Temporary for n
      real w                ! Dummy argument for statement function
      real uu               ! Dummy argument for statement function
      real g                ! Dummy argument for statement function
      real e                ! Dummy argument for statement function
      real f                ! Dummy argument for statement function
      real t                ! Dummy argument for statement function
      real et               ! Dummy argument for statement function
C
C Intermediate terms for delta-eddington solution
C
      real alp              ! Temporary for alpha
      real gam              ! Temporary for gamma
      real ue               ! Temporary for u
      real arg              ! Exponential argument
      real extins           ! Extinction
      real amg              ! Alp - gam
      real apg              ! Alp + gam
C
C Gaussian integration parameters
C
      integer ng            ! gaussian integration index
      integer ngmax         ! number of gaussian angles
      parameter ( ngmax = 8 )
      real gauspt(ngmax)    ! gaussian angles (points)
      real gauswt(ngmax)    ! gaussian weights
      real wt               ! gaussian weight
      real swt              ! sum of weights
      real trn              ! layer transmission
      real rdr              ! rdir for gaussian integration
      real tdr              ! tdir for gaussian integration
      real smr              ! accumulator for rdif gaussian integration
      real smt              ! accumulator for tdif gaussian integration
C
      data gauspt / .9894009,  .9445750,  .8656312,  .7554044,
     $              .6178762,  .4580168,  .2816036,  .0950125 /
      data gauswt / .0271525,  .0622535,  .0951585,  .1246290,
     $              .1495960,  .1691565,  .1826034,  .1894506 /
C
C Delta-Eddington solution expressions
C
      alpha(w,uu,g,e) = .75*w*uu*((1. + g*(1-w))/(1. - e*e*uu*uu))
      gamma(w,uu,g,e) = .50*w*((3.*g*(1.-w)*uu*uu + 1.)/(1.-e*e*uu*uu))
      el(w,g)         = sqrt(3.*(1-w)*(1. - w*g))
      taus(w,f,t)     = (1. - w*f)*t
      omgs(w,f)       = (1. - f)*w/(1. - w*f)
      asys(g,f)       = (g - f)/(1. - f)
      u(w,g,e)        = 1.5*(1. - w*g)/e
      n(uu,et)        = ((uu+1.)*(uu+1.)/et ) - ((uu-1.)*(uu-1.)*et)
C
C-----------------------------------------------------------------------
C
C Initialize all total transmission values to 0, so that nighttime 
C values from previous computations are not used:
C
      call resetr(trntdr,plond*(1+klevp),0.)
C
C Other initialization
C
      do ii=1,nptssi
         i = indxsi(ii)
         if (mu0_in(i).gt.0.) then
C
C Initialize top interface of top layer:
C
          trndir(i,0) =   1.0
          trntdr(i,0) =   1.0
          trndif(i,0) =   1.0
          rdndif(i,0) =   0.0
C
C Compute solar zenith angle in medium
C
          mu0  = mu0_in(i)
          mu0n = sqrt(1.-((1.-mu0*mu0)/(refindx*refindx)))
C
C Compute level of Fresnel refraction
C
          if( srf_in(i) .lt. 2 ) then
            kfrsnl = 1 + ksnow + 1
          else
            kfrsnl = 0
          endif
C
C Compute fresnel reflection and transmission amplitudes
C for two polarizations: 1=perpendicular and 2=parallel to
C the plane containing incident, reflected and refracted rays.
C
          R1 = (mu0 - refindx*mu0n) / 
     $         (mu0 + refindx*mu0n)
          R2 = (refindx*mu0 - mu0n) / 
     $         (refindx*mu0 + mu0n)
          T1 = 2.0*mu0 / 
     $         (mu0 + refindx*mu0n)
          T2 = 2.0*mu0 / 
     $         (refindx*mu0 + mu0n)
C
C Unpolarized light for direct beam
C
          Rf_dir_a = 0.5 * (R1*R1 + R2*R2)
          Tf_dir_a = 0.5 * (T1*T1 + T2*T2)*refindx*mu0n/mu0
C
C Precalculated diffuse reflectivities and transmissivities
C for incident radiation above and below fresnel layer, using
C the direct albedos and accounting for complete internal
C reflection from below; precalculated because high order
C number of gaussian points (~256) is required for convergence:
C
C above
          Rf_dif_a = 0.063
          Tf_dif_a = 1.0 - Rf_dif_a
C below
          Rf_dif_b = 0.455
          Tf_dif_b = 1.0 - Rf_dif_b
C
         endif
C End of initialization
      enddo
C
C Proceed  down one layer at a time; if the total transmission to
C the interface just above a given layer is less than trmin, then no
C Delta-Eddington computation for that layer is done:
C
      do k=0,klev
C
C Initialize current layer properties to zero; only if total
C transmission to the top interface of the current layer exceeds the
C minimum, will these values be computed below:
C
         do ii=1,nptssi
            i = indxsi(ii)
            if (mu0_in(i).gt.0.) then
C
C Calculates the solar beam transmission, total transmission, and
C reflectivity for diffuse radiation from below at the top of the
C current layer:
C
C              layers       interface
C
C       ---------------------  k-1 
C                k-1
C       ---------------------  k
C                 k
C       ---------------------  
C
              if ( k .gt. 0 ) then
                trndir(i,k) = trndir(i,k-1)*trnlay(i,k-1)
                refkm1      = 1./(1. - rdndif(i,k-1)*rdif_a(i,k-1))
                tdrrdir     = trndir(i,k-1)*rdir(i,k-1)
                tdndif      = trntdr(i,k-1) - trndir(i,k-1)
                trntdr(i,k) = trndir(i,k-1)*tdir(i,k-1) + 
     $            (tdndif + tdrrdir*rdndif(i,k-1))*refkm1*tdif_a(i,k-1)
                rdndif(i,k) = rdif_b(i,k-1)  +
     $            (tdif_b(i,k-1)*rdndif(i,k-1)*refkm1*tdif_a(i,k-1))
                trndif(i,k) = trndif(i,k-1)*tdif_a(i,k-1)*refkm1
              endif
C
            endif
         enddo
C
C Compute next layer Delta-Eddington solution only if total transmission
C of radiation to the interface just above the layer exceeds trmin.
C
         call whenfgt(plond,trntdr(1,k),1,trmin,index,nval)
         if(nval.gt.0) then
C
            do ii=1,nval
               i=index(ii)
C
               tautot  = tau_in(i,k)
               wtot    = w_in(i,k)
               gtot    = g_in(i,k)
               ftot    = gtot*gtot
C
               ts   = taus(wtot,ftot,tautot)
               ws   = omgs(wtot,ftot)
               gs   = asys(gtot,ftot)
               lm   = el(ws,gs)
               ue   = u(ws,gs,lm)
C mu0 in medium
               mu0  = mu0_in(i)
               mu0n = sqrt(1.-((1.-mu0*mu0)/(refindx*refindx)))
C do not change if above fresnel layer and in non-pond medium
               if( srf_in(i).ne.2 .and. k.lt.kfrsnl ) mu0n = mu0
C
C Limit arguments of exponentials to 25:
C
               arg  = amin1(lm*ts,25.)
               extins = exp(-arg)
               ne = n(ue,extins)
C
C First calculation of rdif, tdif:
C
               rdif_a(i,k) = (ue+1.)*(ue-1.)*(1./extins - extins)/ne
               tdif_a(i,k) = 4.*ue/ne
C
C Evaluate rdir,tdir:
C
               arg       = amin1(ts/mu0n,25.)
               trnlay(i,k) = exp(-arg)
               alp = alpha(ws,mu0n,gs,lm)
               gam = gamma(ws,mu0n,gs,lm)
               apg = alp + gam
               amg = alp - gam
               rdir(i,k) = amg*(tdif_a(i,k)*trnlay(i,k) - 1.) + 
     $                     apg*rdif_a(i,k)
               tdir(i,k) = apg*tdif_a(i,k) +
     $                     (amg*rdif_a(i,k) - (apg-1.))*trnlay(i,k)
C
C Recalculate rdif,tdif using gaussian integration over rdir,tdir:
C
               swt = 0.0
               smr = 0.0
               smt = 0.0
               do ng=1,ngmax
                 mu  = gauspt(ng)
                 wt  = gauswt(ng)
                 swt = swt + mu*wt
                 arg = amin1(ts/mu,25.)
                 trn = exp(-arg)
                 alp = alpha(ws,mu,gs,lm)
                 gam = gamma(ws,mu,gs,lm)
                 apg = alp + gam
                 amg = alp - gam
                 rdr = amg*(tdif_a(i,k)*trn-1.) + 
     $                 apg*rdif_a(i,k)
                 tdr = apg*tdif_a(i,k) +
     $                 (amg*rdif_a(i,k)-(apg-1.))*trn
                 smr = smr + mu*rdr*wt
                 smt = smt + mu*tdr*wt
               enddo
               rdif_a(i,k) = smr/swt
               tdif_a(i,k) = smt/swt
C
C Homogeneous layer
C
               rdif_b(i,k) = rdif_a(i,k)
               tdif_b(i,k) = tdif_a(i,k)
C
C Add fresnel layer to top of desired layer if either 
C air or snow overlies ice; we ignore refraction in ice 
C if a melt pond overlies it:
C
               if( k .eq. kfrsnl ) then
C
            rintfc      = 1./(1. - Rf_dif_b*rdif_a(i,kfrsnl))
            tdir(i,kfrsnl)   = Tf_dir_a*tdir(i,kfrsnl) + 
     $        Tf_dir_a*rdir(i,kfrsnl)*Rf_dif_b*rintfc*tdif_a(i,kfrsnl)
C
            rdir(i,kfrsnl)   = Rf_dir_a + 
     $        Tf_dir_a*rdir(i,kfrsnl)*rintfc*Tf_dif_b
            rdif_a(i,kfrsnl) = Rf_dif_a + 
     $        Tf_dif_a*rdif_a(i,kfrsnl)*rintfc*Tf_dif_b
            rdif_b(i,kfrsnl) = rdif_b(i,kfrsnl) + 
     $        tdif_b(i,kfrsnl)*Rf_dif_b*rintfc*tdif_a(i,kfrsnl)
            tdif_a(i,kfrsnl) = Tf_dif_a*rintfc*tdif_a(i,kfrsnl)
            tdif_b(i,kfrsnl) = tdif_b(i,kfrsnl)*rintfc*Tf_dif_b
C
C update trnlay to include fresnel transmission
C
            trnlay(i,kfrsnl) = Tf_dir_a*trnlay(i,kfrsnl)
C
               endif
            enddo
         else
            rdir(i,k)   =  0.0
            rdif_a(i,k) =  0.0
            rdif_b(i,k) =  0.0
            tdir(i,k)   =  0.0
            tdif_a(i,k) =  0.0
            tdif_b(i,k) =  0.0
            trnlay(i,k) =  0.0
C end transmission large enough
         endif
C end level k loop
      enddo
C
C Compute total direct beam transmission, total transmission, and
C reflectivity for diffuse radiation (from below) for all layers
C above the surface; note that we ignore refraction between ice
C and underlying ocean:
C
C       For k = klevp
C
C              layers       interface
C
C       ---------------------  k-1 
C                k-1
C       ---------------------  k
C       \\\\\\\ ocean \\\\\\\
C
      k = klevp
      do ii=1,nptssi
         i = indxsi(ii)
         if (mu0_in(i).gt.0.) then
           trndir(i,k) = trndir(i,k-1)*trnlay(i,k-1)
           refkm1      = 1./(1. - rdndif(i,k-1)*rdif_a(i,k-1))
           tdrrdir     = trndir(i,k-1)*rdir(i,k-1)
           tdndif      = trntdr(i,k-1) - trndir(i,k-1)
           trntdr(i,k) = trndir(i,k-1)*tdir(i,k-1) + 
     $       (tdndif + tdrrdir*rdndif(i,k-1))*refkm1*tdif_a(i,k-1)
           rdndif(i,k) = rdif_b(i,k-1)  +
     $       (tdif_b(i,k-1)*rdndif(i,k-1)*refkm1*tdif_a(i,k-1))
           trndif(i,k) = trndif(i,k-1)*tdif_a(i,k-1)*refkm1
         endif
      enddo
C
C Compute reflectivity to direct and diffuse radiation for layers below by
C adding succesive layers starting from the surface and working upwards:
C
C              layers       interface
C
C       ---------------------  k
C                 k
C       ---------------------  k+1
C                k+1
C       ---------------------
C
      do ii=1,nptssi
        i = indxsi(ii)
        if (mu0_in(i).gt.0.) then
          rupdir(i,klevp) = albodr(i)
          rupdif(i,klevp) = albodf(i)
          do k=klev,0,-1
C
C interface scattering
C
            refkp1 = 1./( 1. - rdif_b(i,k)*rupdif(i,k+1))
C
C dir from top layer plus exp tran ref from lower layer, interface
C scattered and tran thru top layer from below, plus diff tran ref
C from lower layer with interface scattering tran thru top from below
C
            rupdir(i,k) = rdir(i,k) + ( trnlay(i,k)*rupdir(i,k+1) +
     $        (tdir(i,k)-trnlay(i,k))*rupdif(i,k+1) ) *
     $        refkp1*tdif_b(i,k)
C
C dif from top layer from above, plus dif tran upwards reflected and
C interface scattered which tran top from below
C
            rupdif(i,k) = rdif_a(i,k) +
     $        tdif_a(i,k)*rupdif(i,k+1)*refkp1*tdif_b(i,k)
          enddo
        endif
      enddo
C
C done with simded
C
      return
      end
      block data blkdat
C-----------------------------------------------------------------------
C
C Contains all data statements for common variables
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Rosinski
C Standardized:      L. Bath, June 1992
C
C-----------------------------------------------------------------------
c
c $Id: blkdat.F,v 1.1.1.1 1995/02/09 23:26:36 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
C------------------------------Parameters-------------------------------
c
c $Id: pmgrid.h,v 1.2 1995/02/10 01:09:06 ccm2 Exp $
c $Author: ccm2 $
c
C
C Basic grid point resolution parameters
C
      integer
     $     plon,    ! number of longitudes
     $     plev,    ! number of vertical levels
     $     plat,    ! number of latitudes
     $     pcnst,   ! number of constituents (including water vapor)
     $     plevmx,  ! number of subsurface levels
     $     plevp,   ! plev + 1
     $     nxpt,    ! no.of points outside active domain for interpolant
     $     jintmx,  ! number of extra latitudes in polar region
     $     plond,   ! slt extended domain longitude
     $     platd,   ! slt extended domain lat.
     $     plevd    ! fold plev,pcnst indices into one
C
      parameter(plon   = 1,
     $          plev   = 18,
     $          plat   = 1,
     $          pcnst  = 1,
     $          plevmx = 4,
     $          plevp  = plev + 1,
     $          nxpt   = 1,
     $          jintmx = 1,
     $          plond  = plon + 1 + 2*nxpt,
     $          platd  = plat + 2*nxpt + 2*jintmx,
     $          plevd  = plev*(3 + pcnst))
C
c
c $Id: pagrid.h,v 1.3 1995/03/03 17:47:17 bonan Exp $
c $Author: bonan $
c
C
C Model grid point resolution parameters.
C
      integer
     $     plnlv,    ! Length of multilevel field slice
     $     plndlv,   ! Length of multilevel 3-d field slice
     $     pbflnb,   ! Length of buffer 1
     $     pbflna,   ! Length of buffer 2
     $     pbflnm1,  ! Length of buffer m1
     $     pflenb,   ! Length of buffer 1, padded for unblocked I/O
     $     pflena,   ! Length of buffer 2, padded for unblocked I/O
     $     plenalcl, ! Length of buffer 2, needed in SPEGRD
     $     ptifld,   ! Number of fields on time-invariant boundary dataset
     $     ptvsfld,  ! Number of fields on time-variant boundary dataset
     $     ptvofld,  ! Number of fields on ozone dataset
     $     plenhi,   ! Length of integer header record
     $     plenhc,   ! Length of character header record
     $     plenhr,   ! Length of real header record
     $     plexbuf,  ! Length of communication buffer for flux coupling
     $     ptapes,   ! Maximum number of history tapes allowed
     $     pflds     ! Number of fields in master field list
      integer
     $     ptileni,  ! Length of time-invariant integer header
     $     ptilenc,  ! Length of time-invariant character header
     $     ptvoleni, ! Length of ozone integer header
     $     ptvolenc, ! Length of ozone character header
     $     ptvsleni, ! Length of time-variant integer header
     $     ptvslenc  ! Length of time-variant character header
      integer
     $     plenhis,  ! Length of integer header scalars
     $     plenhcs,  ! Length of character header scalars
     $     ptilenis, ! Length of time-invariant integer scalars
     $     ptilencs, ! Length of time-invariant character scalars
     $     ptolenis, ! Length of ozone integer header scalars
     $     ptolencs, ! Length of ozone character header scalars
     $     ptslenis, ! Length of time-variant integer header scalars
     $     ptslencs  ! Length of time-variant character header scalars
C
      parameter(plnlv=plon*plev,plndlv=plond*plev)
C
C In pbflnb, 9 multi-level fields include the plev levels of plol and
C plos. 2 multi-level fields are pcnst-dependent.
C There are plevmx sub-surface temperature fields. (See User's Guide 
C for complete buffer description)
C There are 4 single-level fields to hold albedos
C
      parameter(pbflnb=(7 + 2*pcnst)*plndlv + (15+plevmx+pcnst)*plond,
C
C In pbflna, there are 3 multi-level and 3 single-level fields.
C
     $          pbflna = (3 + 3*plev)*plond,
     $          pbflnm1 = (1 + 2*plev)*plond,
     $          pflenb = ((pbflnb + pbflnm1)/512 + 1)*512,
     $          pflena = (pbflna/512 + 1)*512,
C
C plenalcl is the buffer size as required in SPEGRD.  
C Only pflena is read/written.
C
     $          plenalcl = ((pbflna + 3*plndlv + plond)/512 + 1)*512,
     $          plexbuf = (((1 + 7*plev)*plond)/512+1)*512,
     $          ptapes = 6,
C
C 8 fields in master list are pcnst-dependent 2 fields occur only
C if pcnst > 1
C
     $          pflds=92+8*pcnst+2*(pcnst-1)+plevmx)
      parameter(ptifld = 11, ptvsfld = 1, ptvofld = 2)
C
C There are 37 scalar words in the integer header and 89 scalar words
C in the character header
C
      parameter(plenhis=37,plenhcs=89,
     $          plenhi=plenhis+3*pflds,plenhc=plenhcs+2*pflds,
     $          plenhr=3*(2*plev + 1) + 2*plat,
     $          ptilenis=plenhis, ptilencs=plenhcs,
     $          ptileni=ptilenis+3*ptifld, ptilenc=ptilencs+2*ptifld,
     $          ptolenis=plenhis, ptolencs=plenhcs,
     $          ptvoleni=ptolenis+3*ptvofld,ptvolenc=ptolencs+2*ptvofld,
     $          ptslenis=plenhis, ptslencs=plenhcs,
     $          ptvsleni=ptslenis+3*ptvsfld,ptvslenc=ptslencs+2*ptvsfld)







C------------------------------Commons----------------------------------
c
c $Id: comhst.h,v 1.1.1.1 1995/02/09 23:26:42 ccm2 Exp $
c $Author: ccm2 $
c
C
C Integer and logical variables related to history tapes
C
      integer pichsum  ! Max. value of 4*ichar(character)
      parameter (pichsum=508)
C
      common /comhst/
     $   nhtfrq(ptapes)    ,mfilt(ptapes) ,nlfilt             ,
     $   ndens(ptapes)     ,nflds(ptapes) ,nfils(ptapes)      ,
     $   hunit(ptapes)     ,nrlen(ptapes) ,nplen(ptapes)      ,
     $   sunit             ,stfnum        ,mtapes             ,
     $   nexcl             ,nincl         ,hbufpt(ptapes)     ,
     $   nacs(pflds,plat)  ,iflds(3,pflds),nupnt(pflds,ptapes),
     $   npnt(pflds,ptapes),ndcurf(ptapes),ncdatf(ptapes)     ,
     $   nscurf(ptapes)    ,ncsecf(ptapes),nfldsc(0:pichsum,ptapes),
     $   islocc(0:pichsum,ptapes)         ,hstwr(ptapes)      ,
     $   rstwr             ,nacsav(pflds,plat)
C
      integer nhtfrq,  ! Array of write frequencies
     $        mfilt    ! Number of write-ups per volume
      logical nlfilt   ! Flag for extra file on 1st vol (ktape=1)
      logical hstwr    ! Flag for history writes
      logical rstwr    ! Flag for restart writes
      integer ndens,   ! Array of input packing densities
     $        nflds,   ! Array of total fields on tape
     $        nfils,   ! Array of current files on the volume
     $        hunit,   ! History tape disk units
     $        nrlen,   ! Record length
     $        nplen,   ! Packed record length,
     $        sunit,   ! History tape SSD unit
     $        stfnum,  ! Starting number for history tape naming
     $        mtapes,  ! Actual number of tapes requested
     $        nexcl,   ! Actual number of excluded fields
     $        nincl,   ! Actual number of included primary tape fields
     $        hbufpt,  ! Ptrs to start of fields for each tape in hbuf
     $        nacs,    ! Number of accumulations for field
     $        nacsav,  ! Saved accumulations for restart
     $        iflds,   ! Integer portion of master field list
     $        nupnt,   ! Array of unpacked field pointers
     $        npnt,    ! Array of packed field pointers
     $        ndcurf,  ! First "current" day for each tape
     $        ncdatf,  ! First "current" date for each tape
     $        nscurf,  ! First "current" second of day for each tape
     $        ncsecf,  ! First "current" second of date for each tape
     $        nfldsc,  ! Number of fields starting with given ichar(1-4)
     $        islocc   ! Index of starting location for each ichar sum
C
C  Character variables related to history tapes
C
      common /comhtc/
     $   nfpath(ptapes)     ,ppath(ptapes)       ,cpath(ptapes)       ,
     $   nhfil(ptapes)      ,ninavg(ptapes)      ,caseid              ,
     $   ctitle             ,fieldn(2,pflds)     ,exclude(pflds)      ,
     $   primary(pflds)     ,aux(pflds,ptapes-1)
C
      character*80 nfpath,    ! Array of first pathnames, for header
     $             ppath,     ! Array of previous pathnames, for header
     $             cpath      ! Array of current pathnames
      character    nhfil*6,   ! Array of current file names
     $             ninavg*1,  ! Tape fields instantaneous or averaged
     $             caseid*8   ! Case identifier
      character*80 ctitle     ! Case title
      character*8  fieldn,    ! Character portion of master field list
     $             exclude,   ! List of fields to rm from primary tape
     $             primary,   ! List of fields to add to primary tape
     $             aux        ! Lists of fields for auxiliary tapes
C
C-----------------------------------------------------------------------
c
c $Id: comlun.h,v 1.1.1.1 1995/02/09 23:26:42 ccm2 Exp $
c $Author: ccm2 $
c
C
C Logical unit numbers and related variables
C
      integer
     $     pnrg1,             ! maximum number of primary 
C                             !  regeneration files
     $     pnrg2,             ! maximum number of secondary 
C                             !  regeneration files
     $     pnrg3              ! maximum number of secondary 

      parameter (pnrg1 = 5)
      parameter (pnrg2 = 5)
      parameter (pnrg3 = 5)
C
      common/comlun/nsds    ,nrg     ,nrg1(pnrg1)      ,nrg2(pnrg2),
     $              nrg3(pnrg3,ptapes)        ,nra1    ,nrb1    ,
     $              ninit   ,nbndti  ,nozone  ,nsst    ,nabem   ,
     $              nsplit,lutag(99)
      common/comlun/rg1lat(pnrg1+1)  ,rg1siz(pnrg1)    ,rg1buf  ,nnrg1,
     $              rg2lat(pnrg2+1)  ,rg2siz(pnrg2)    ,rg2buf  ,nnrg2,
     $              rg3lat(pnrg3+1,ptapes)    ,rg3siz(pnrg3,ptapes)   ,
     $              rg3buf(ptapes)   ,nnrg3(ptapes)    ,mxszrg  ,
     $              nrefrq  ,rgnht(ptapes)    ,rg3num  ,mresfq
      common/comlunc/rg1ext(pnrg1)   ,rg2ext(pnrg2)    ,
     $               rg3ext(pnrg3,ptapes)
C
      integer nsds,     ! restart dataset unit
     $        nrg,      ! master regeneration dataset unit
     $        nrg1,     ! primary regeneration dataset units
     $        nrg2,     ! secondary regeneration dataset units
     $        nrg3,     ! hbuf regeneration dataset units
     $        nra1,     ! a work file
     $        nrb1,     ! b work file
     $        ninit,    ! initial dataset unit
     $        nbndti,   ! time-invariant boundary dataset
     $        nozone,   ! ozone dataset
     $        nsst,     ! sst dataset
     $        nabem,    ! absorptivity/emissivity work file
     $        nsplit    ! communication between LINEMS1 and LINEMS2
C
      logical
     $     lutag        ! list of flags marking logical units in use
      integer
     $     rg1lat,      ! latitude list for primary regen datasets
     $     rg1siz,      ! file sizes for preallocation
     $     rg1buf,      ! buffer length for assign
     $     nnrg1,       ! number of primary regen files written
     $     rg2lat,      ! lat list for secondary regen datasets
     $     rg2siz,      ! file size for preallocation
     $     rg2buf,      ! buffer length for assign
     $     nnrg2,       ! number of secondary regen files written
     $     rg3lat,      ! latitude list for hbuf regen datasets
     $     rg3siz,      ! file sizes for preallocation
     $     rg3buf,      ! buffer length for assign
     $     nnrg3,       ! number of hbuf regen files written
     $     mxszrg,      ! max size of a regen file (megabytes)
     $     nrefrq,      ! frequency of regeneration file writes
     $     mresfq,      ! frequency of mnthly avg regen file writes
     $     rg3num       ! number of temporary secondary regen files
      logical
     $     rgnht        ! set true if regeneration file for a h-tape exists
      character*2
     $     rg1ext,      ! file extension for primary regen files
     $     rg2ext       ! file extension for secondary regen files
      character*5
     $     rg3ext       ! file extension for secondary regen files
C
C-----------------------------------------------------------------------
c
c $Id: compbl.h,v 1.1.1.1 1995/02/09 23:26:43 ccm2 Exp $
c $Author: ccm2 $
c
C
C Pbl constants
C
      common /compbl/ betam   ,betas   ,betah   ,fak     ,g       ,
     $                onet    ,fakn    ,ricr    ,sffrac  ,vk      ,
     $                ccon    ,binm    ,binh
C
      real betam,  ! constant in wind gradient expression
     $     betas,  ! constant in surface layer gradient expression
     $     betah,  ! constant in temperature gradient expression 
     $     fak,    ! constant in surface temperature excess         
     $     g,      ! gravitational acceleration
     $     onet,   ! 1/3 power in wind gradient expression
     $     fakn,   ! constant in turbulent prandtl number
     $     ricr,   ! critical richardson number
     $     sffrac, ! surface layer fraction of boundary layer
     $     vk,     ! von karman's constant
     $     ccon,   ! fak * sffrac * vk
     $     binm,   ! betam * sffrac
     $     binh    ! betah * sffrac
C
C-----------------------------------------------------------------------
c
c $Id: crdcae.h,v 1.1.1.1 1995/02/09 23:26:46 ccm2 Exp $
c $Author: ccm2 $
c
C
C Water vapor narrow band constants for longwave radiation computations
C
      common/crdcae/realk(2), st(2), a1(2), a2(2), b1(2), b2(2),
     $              coefa(3,4),coefb(4,4),coefc(3,4),coefd(4,4),
     $              coefe(3,4),coeff(6,2),coefg(2,4),coefh(2,4),
     $              coefi(6,2),coefj(3,2),coefk(3,2),
     $              c1(4),c2(4),c3(4),c4(4),c5(4),c6(4),c7(4),c8,c9,
     $              c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,
     $              c22,c23,c24,c25,c26,c27,c28,c29,c30,c31,
     $              fwcoef,fwc1,fwc2,fc1,cfa1
C
      real realk,     ! H2O narrow band parameter
     $     st,        ! H2O narrow band parameter
     $     a1,a2,     ! Temperature correction terms for H2O path
     $     b1,b2      ! Temperature correction terms for H2O path
C
C Constant coefficients for water vapor absorptivity and emissivity
C
      real coefa,coefb,coefc,coefd,coefe,coeff,
     $     coefg,coefh,coefi,coefj,coefk,
     $      c1, c2, c3, c4, c5, c6, c7,c8 ,c9 ,c10,
     $     c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,
     $     c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,c31
C
C Farwing correction constants for narrow-band emissivity model,
C introduced to account for the deficiencies in narrow-band model
C used to derive the emissivity; tuned with Arking's line-by-line
C calculations.
C
      real fwcoef,
     $     fwc1,fwc2,
     $     fc1,
     $     cfa1
C
C-----------------------------------------------------------------------
C
C Logical units from /comlun/
C
      data lutag/99*.FALSE./
C
C /COMHST/
C
      data mtapes/0/
C
C /PBL/
C
      data vk, fak, fakn / 0.40, 8.5 , 7.2/
      data betam, betas, betah, sffrac / 15.0, 5.0, 15.0, 0.1 /
      data  ricr /  0.3 /
C
C /CRDCAE/
C H2O EMISSIVITY AND ABSORTIVITY COEFFICIENTS
C
      data coefa/1.01400e+00, 6.41695e-03, 2.85787e-05,
     $           1.01320e+00, 6.86400e-03, 2.96961e-05,
     $           1.02920e+00, 1.01680e-02, 5.30226e-05,
     $           1.02743e+00, 9.85113e-03, 5.00233e-05/
C
      data coefb/8.85675e+00,-3.51620e-02, 2.38653e-04,-1.71439e-06,
     $           5.73841e+00,-1.91919e-02, 1.65993e-04,-1.54665e-06,
     $           6.64034e+00, 1.56651e-02,-9.73357e-05, 0.0,
     $           7.09281e+00, 1.40056e-02,-1.15774e-04, 0.0/
C
      data coefc/9.90127e-01, 1.22475e-03, 4.90135e-06,
     $           9.89753e-01, 1.97081e-03, 3.42046e-06,
     $           9.75230e-01, 1.03341e-03, 0.0,
     $           9.77366e-01, 8.60014e-04, 0.0/
C
      data coefd/7.03047e-01,-2.63501e-03,-1.57023e-06, 0.0,
     $           5.29269e-01,-3.14754e-03, 4.39595e-06, 0.0,
     $           7.88193e-02, 1.31290e-03, 4.25827e-06,-1.23982e-08,
     $           1.62744e-01, 2.22847e-03, 2.60102e-06,-4.30133e-08/
C
      data coefe/3.93137e-02,-4.34341e-05, 3.74545e-07,
     $           3.67785e-02,-3.10794e-05, 2.94436e-07,
     $           7.42500e-02, 3.97397e-05, 0.0,
     $           7.52859e-02, 4.18073e-05, 0.0/
C
      data coeff/2.2037 e-01,1.39719e-03,-7.32011e-06,
     $          -1.40262e-08,2.13638e-10,-2.35955e-13,
     $           3.07431e-01,8.27225e-04,-1.30067e-05,
     $           3.49847e-08,2.07835e-10,-1.98937e-12/
C
      data coefg/9.04489e+00,-9.56499e-03,
     $           1.80898e+01,-1.91300e-02,
     $           8.72239e+00,-9.53359e-03,
     $           1.74448e+01,-1.90672e-02/
C
      data coefh/5.46557e+01,-7.30387e-02,
     $           1.09311e+02,-1.46077e-01,
     $           5.11479e+01,-6.82615e-02,
     $           1.02296e+02,-1.36523e-01/
C
      data coefi/3.31654e-01,-2.86103e-04,-7.87860e-06,
     $           5.88187e-08,-1.25340e-10,-1.37731e-12,
     $           3.14365e-01,-1.33872e-03,-2.15585e-06,
     $           6.07798e-08,-3.45612e-10,-9.34139e-15/
C
      data coefj/2.82096e-02,2.47836e-04,1.16904e-06,
     $           9.27379e-02,8.04454e-04,6.88844e-06/
C
      data coefk/2.48852e-01,2.09667e-03,2.60377e-06,
     $           1.03594e+00,6.58620e-03,4.04456e-06/
C
C Narrow band data for H2O
C 200CM data for 800-1000 CM-1 and 1000-1200 CM-1.
C
      data realk/  0.18967069430426e-04, 0.70172244841851e-04   /
      data   st /  0.31930234492350e-03, 0.97907319939060e-03   /
      data   a1 /  0.28775403075736e-01, 0.23236701470511e-01   /
      data   a2 / -0.57966222388131e-04,-0.95105504388411e-04   /
      data   b1 /  0.29927771523756e-01, 0.21737073577293e-01   /
      data   b2 / -0.86322071248593e-04,-0.78543550629536e-04   /
      end
      subroutine cldefr(ioro, t, rel, rei, fice, ps, pmid)
C-----------------------------------------------------------------------
C
C Compute cloud drop size
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Kiehl, January 1993
C
C-----------------------------------------------------------------------
c
c $Id: cldefr.F,v 1.2 1995/02/17 21:28:07 jhack Exp $
c $Author: jhack $
c
      implicit none
C------------------------------Parameters-------------------------------
c
c $Id: pmgrid.h,v 1.2 1995/02/10 01:09:06 ccm2 Exp $
c $Author: ccm2 $
c
C
C Basic grid point resolution parameters
C
      integer
     $     plon,    ! number of longitudes
     $     plev,    ! number of vertical levels
     $     plat,    ! number of latitudes
     $     pcnst,   ! number of constituents (including water vapor)
     $     plevmx,  ! number of subsurface levels
     $     plevp,   ! plev + 1
     $     nxpt,    ! no.of points outside active domain for interpolant
     $     jintmx,  ! number of extra latitudes in polar region
     $     plond,   ! slt extended domain longitude
     $     platd,   ! slt extended domain lat.
     $     plevd    ! fold plev,pcnst indices into one
C
      parameter(plon   = 1,
     $          plev   = 18,
     $          plat   = 1,
     $          pcnst  = 1,
     $          plevmx = 4,
     $          plevp  = plev + 1,
     $          nxpt   = 1,
     $          jintmx = 1,
     $          plond  = plon + 1 + 2*nxpt,
     $          platd  = plat + 2*nxpt + 2*jintmx,
     $          plevd  = plev*(3 + pcnst))
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer ioro(plond)       ! nint(oro(i))
      real t(plond,plev)        ! Temperature
      real ps(plond),           ! surface pressure
     $     pmid(plond,plev)     ! midpoint pressures
C
C Output arguments
C
      real rel(plond,plev),     ! liquid effective drop size (microns)
     $     rei(plond,plev),     ! ice effective drop size (microns)
     $     fice(plond,plev)     ! fractional ice content within cloud
      real pirnge,              ! nrmlzd pres range for ice particle changes
     $     picemn,              ! normalized pressure below which rei=reimax
     $     rirnge,              ! range of ice radii (reimax - 10 microns)
     $     reimax,              ! maximum ice effective radius
     $     pnrml,               ! normalized pressure
     $     weight               ! coef. for determining rei as fn of P/PS
C
C---------------------------Local workspace-----------------------------
C
      integer i,k               ! longitude, level indices
      real rliq                 ! temporary liquid drop size
      real pi                   ! pi
C
C-----------------------------------------------------------------------
C
      pi = 4. * atan(1.)
      do k=1,plev
         do i=1,plon
c
c Define liquid drop size
c
            if(ioro(i).ne.1) then
c
c     Effective liquid radius over ocean and sea ice
c
c
c changed 27 April 2004 to agree with 7.0 micron 
c liquid drop size used in Jin et al. 1994
c
              rliq =   7.0
c...              rliq =  10.0
            else
c
c     Effective liquid radius over land
c
              rliq  = 5.0 + 
     $                5.0*amin1(1.0,amax1(0.0,(263.16-t(i,k))*0.05))
            endif
c
            rel(i,k) = rliq
C+            rei(i,k) = 30.0
c
c     Determine rei as function of normalized pressure
c
            reimax   = 30.0
            rirnge   = 20.0 
            pirnge   = 0.4
            picemn   = 0.4
c
            pnrml    = pmid(i,k)/ps(i)
            weight   = amax1(amin1((pnrml-picemn)/pirnge,1.0),0.)
            rei(i,k) = reimax - rirnge*weight
c
c Define fractional amount of cloud that is ice
c
c if warmer than -10 degrees C then water phase
c
             if(t(i,k).gt.263.16) fice(i,k) = 0.0
c
c if colder than -10 degrees C but warmer than -30 C mixed phase
c
             if(t(i,k).le.263.16.and.t(i,k).ge.243.16) then
                   fice(i,k) =(263.16-t(i,k)) / 20.0
             endif
c
c if colder than -30 degrees C then ice phase
c
             if(t(i,k).lt.243.16) fice(i,k) = 1.0
c
c Turn off ice radiative properties by setting fice = 0.0
c
C+             fice(i,k) = 0.0
c
         end do
      end do
C
      return
      end
      subroutine cldems(clwp, fice, rei, emis)
C-----------------------------------------------------------------------
C
C Compute cloud emissivity using cloud liquid water path (g/m**2)
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Kiehl
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Hack, J. Kiehl, August 1992
C
C-----------------------------------------------------------------------
c
c $Id: cldems.F,v 1.2 1995/02/17 21:28:09 jhack Exp $
c $Author: jhack $
c
      implicit none
C------------------------------Parameters-------------------------------
c
c $Id: pmgrid.h,v 1.2 1995/02/10 01:09:06 ccm2 Exp $
c $Author: ccm2 $
c
C
C Basic grid point resolution parameters
C
      integer
     $     plon,    ! number of longitudes
     $     plev,    ! number of vertical levels
     $     plat,    ! number of latitudes
     $     pcnst,   ! number of constituents (including water vapor)
     $     plevmx,  ! number of subsurface levels
     $     plevp,   ! plev + 1
     $     nxpt,    ! no.of points outside active domain for interpolant
     $     jintmx,  ! number of extra latitudes in polar region
     $     plond,   ! slt extended domain longitude
     $     platd,   ! slt extended domain lat.
     $     plevd    ! fold plev,pcnst indices into one
C
      parameter(plon   = 1,
     $          plev   = 18,
     $          plat   = 1,
     $          pcnst  = 1,
     $          plevmx = 4,
     $          plevp  = plev + 1,
     $          nxpt   = 1,
     $          jintmx = 1,
     $          plond  = plon + 1 + 2*nxpt,
     $          platd  = plat + 2*nxpt + 2*jintmx,
     $          plevd  = plev*(3 + pcnst))
C
      real kabs                   ! longwave absorption coeff (m**2/g)
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real clwp(plond,plev),      ! cloud liquid water path (g/m**2)
     $     rei(plond,plev),       ! ice effective drop size (microns)
     $     fice(plond,plev)       ! fractional ice content within cloud
C
C Output arguments
C
      real emis(plond,plev)       ! cloud emissivity (fraction)
C
C---------------------------Local workspace-----------------------------
C
      integer i,k                 ! longitude, level indices
      real kabsl,                 ! longwave absorption coeff (m**2/g)
     $     kabsi                  ! ice absorption coefficient
      parameter (kabsl = 0.090361)
C
C-----------------------------------------------------------------------
C
      do k=1,plev
         do i=1,plon
            kabsi = 0.005 + 1./rei(i,k)
            kabs = kabsl*(1.-fice(i,k))+kabsi*fice(i,k)
            emis(i,k) = 1. - exp(-1.66*kabs*clwp(i,k))
         end do
      end do
C
      return
      end
      subroutine fmrgrid(qrs     ,qrl     ,
     $                   qrsm    ,qrlm    )
C-----------------------------------------------------------------------
C
C Interpolate model arrays to radiation vertical grid.
C
C------------------------------Parameters-------------------------------
c
c $Id: fmrgrid.F,v 1.1.1.1 1995/02/09 23:26:49 ccm2 Exp $
c $Author: ccm2 $
c
c
c $Id: pmgrid.h,v 1.2 1995/02/10 01:09:06 ccm2 Exp $
c $Author: ccm2 $
c
C
C Basic grid point resolution parameters
C
      integer
     $     plon,    ! number of longitudes
     $     plev,    ! number of vertical levels
     $     plat,    ! number of latitudes
     $     pcnst,   ! number of constituents (including water vapor)
     $     plevmx,  ! number of subsurface levels
     $     plevp,   ! plev + 1
     $     nxpt,    ! no.of points outside active domain for interpolant
     $     jintmx,  ! number of extra latitudes in polar region
     $     plond,   ! slt extended domain longitude
     $     platd,   ! slt extended domain lat.
     $     plevd    ! fold plev,pcnst indices into one
C
      parameter(plon   = 1,
     $          plev   = 18,
     $          plat   = 1,
     $          pcnst  = 1,
     $          plevmx = 4,
     $          plevp  = plev + 1,
     $          nxpt   = 1,
     $          jintmx = 1,
     $          plond  = plon + 1 + 2*nxpt,
     $          platd  = plat + 2*nxpt + 2*jintmx,
     $          plevd  = plev*(3 + pcnst))
C
c
c $Id: ptrrgrid.h,v 1.1.1.1 1995/02/09 23:26:59 ccm2 Exp $
c $Author: ccm2 $
c
C
C Define radiation vertical grid and buffer length for abs/ems out-of-core file
C
      integer
     $     plevr,   ! number of vertical levels
     $     plevrp,  ! plevr + 1
     $     plngbuf  ! length of absorptivity/emissivity record
C
      parameter(plevr = 18,
     $          plevrp = plevr + 1,
     $          plngbuf = 512*((plond*plevrp*plevrp + plond*plevr*4 +
     $                          plond*plevrp)/512 + 1))
C
C-----------------------------------------------------------------------
C
C Input arguments (radiation grid)
C
      real qrs(plond,plevr),     ! Shortwave heating rate
     $     qrl(plond,plevr)      ! Longwave heating rate
C
C Output arguments (model grid)
C
      real qrsm(plond,plev),     ! Shortwave heating rate
     $     qrlm(plond,plev)      ! Longwave heating rate
C
C Code to interpolate goes here.  Do nothing could be coded as a memory
C transfer, but left out for efficiency considerations.
C
      return
      end








      subroutine radabs(pbr    ,pnm     ,co2em    ,co2eml  ,tplnka  ,
     $                  s2c    ,s2t     ,w        ,h2otr   ,plco2   ,
     $                  plh2o  ,co2t    ,tint     ,tlayr   ,plol    ,
     $                  plos   ,pmln    ,piln     ,ucfc11  ,ucfc12  , 
     $                  un2o0  ,un2o1   ,uch4     ,uco211  ,uco212  ,
     $                  uco213 ,uco221  ,uco222   ,uco223  ,uptype  ,
     $                  bn2o0  ,bn2o1   ,bch4    ,abplnk1  ,abplnk2 ,
     $                  abstot ,absnxt  )
C-----------------------------------------------------------------------
C
C Compute absorptivities for h2o, co2, and o3
C
C h2o  ....  Uses nonisothermal emissivity for water vapor from
C            Ramanathan, V. and  P.Downey, 1986: A Nonisothermal
C            Emissivity and Absorptivity Formulation for Water Vapor
C            Journal of Geophysical Research, vol. 91., D8, pp 8649-8666
C
C co2  ....  Uses absorptance parameterization of the 15 micro-meter
C            (500 - 800 cm-1) band system of Carbon Dioxide, from
C            Kiehl, J.T. and B.P.Briegleb, 1991: A New Parameterization
C            of the Absorptance Due to the 15 micro-meter Band System
C            of Carbon Dioxide Jouranl of Geophysical Research,
C            vol. 96., D5, pp 9013-9019
C
C o3   ....  Uses absorptance parameterization of the 9.6 micro-meter
C            band system of ozone, from Ramanathan, V. and R.Dickinson,
C            1979: The Role of stratospheric ozone in the zonal and
C            seasonal radiative energy balance of the earth-troposphere
C            system. Journal of the Atmospheric Sciences, Vol. 36,
C            pp 1084-1104
C
C Computes individual absorptivities for non-adjacent layers, accounting
C for band overlap, and sums to obtain the total; then, computes the
C nearest layer contribution.
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Kiehl, B. Briegleb, August 1992
C
C-----------------------------------------------------------------------
c
c $Id: radabs.F,v 1.2 1995/02/17 21:28:25 jhack Exp $
c $Author: jhack $
c
      implicit none
C------------------------------Parameters-------------------------------
c
c $Id: prgrid.h,v 1.2 1995/02/10 01:09:10 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation resolution and I/O parameters
C
      integer
     $     plon,    ! number of longitudes
     $     plev,    ! number of vertical levels
     $     plat,    ! number of latitudes
     $     pcnst,   ! number of constituents (including water vapor)
     $     plevmx,  ! number of subsurface levels
     $     plevp,   ! plev + 1
     $     nxpt,    ! no.of points outside active domain for interpolant
     $     jintmx,  ! number of extra latitudes in polar region
     $     plond,   ! slt extended domain longitude
     $     platd,   ! slt extended domain lat.
     $     plevd    ! fold plev,pcnst indices into one
      integer
     $     plngbuf  ! length of absorptivity/emissivity record
C
      parameter(plon    = 1,
     $          plev    = 18,
     $          plat    = 1,
     $          pcnst   = 1,
     $          plevmx  = 4,
     $          plevp   = plev + 1,
     $          nxpt    = 1,
     $          jintmx  = 1,
     $          plond   = plon + 1 + 2*nxpt,
     $          platd   = plat + 2*nxpt + 2*jintmx,
     $          plevd   = plev*(3 + pcnst))
      parameter(plngbuf = 512*((plond*plevp*plevp + plond*plev*4 +
     $                          plond*plevp)/512 + 1))
C
C------------------------------Commons----------------------------------
c
c $Id: crdcae.h,v 1.1.1.1 1995/02/09 23:26:46 ccm2 Exp $
c $Author: ccm2 $
c
C
C Water vapor narrow band constants for longwave radiation computations
C
      common/crdcae/realk(2), st(2), a1(2), a2(2), b1(2), b2(2),
     $              coefa(3,4),coefb(4,4),coefc(3,4),coefd(4,4),
     $              coefe(3,4),coeff(6,2),coefg(2,4),coefh(2,4),
     $              coefi(6,2),coefj(3,2),coefk(3,2),
     $              c1(4),c2(4),c3(4),c4(4),c5(4),c6(4),c7(4),c8,c9,
     $              c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,
     $              c22,c23,c24,c25,c26,c27,c28,c29,c30,c31,
     $              fwcoef,fwc1,fwc2,fc1,cfa1
C
      real realk,     ! H2O narrow band parameter
     $     st,        ! H2O narrow band parameter
     $     a1,a2,     ! Temperature correction terms for H2O path
     $     b1,b2      ! Temperature correction terms for H2O path
C
C Constant coefficients for water vapor absorptivity and emissivity
C
      real coefa,coefb,coefc,coefd,coefe,coeff,
     $     coefg,coefh,coefi,coefj,coefk,
     $      c1, c2, c3, c4, c5, c6, c7,c8 ,c9 ,c10,
     $     c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,
     $     c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,c31
C
C Farwing correction constants for narrow-band emissivity model,
C introduced to account for the deficiencies in narrow-band model
C used to derive the emissivity; tuned with Arking's line-by-line
C calculations.
C
      real fwcoef,
     $     fwc1,fwc2,
     $     fc1,
     $     cfa1
C
C-----------------------------------------------------------------------
c
c $Id: crdcon.h,v 1.1.1.1 1995/02/09 23:26:46 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation constants
C
      common/crdcon/gravit  ,rga     ,cpair   ,epsilo  ,sslp    ,
     $              stebol  ,rgsslp  ,co2vmr  ,dpfo3   ,dpfco2  ,
     $              dayspy  ,pie
C
      real gravit,    ! Acceleration of gravity
     $     rga,       ! 1./gravit
     $     cpair,     ! Specific heat of dry air
     $     epsilo,    ! Ratio of mol. wght of H2O to dry air
     $     sslp,      ! Standard sea-level pressure
     $     stebol,    ! Stefan-Boltzmann's constant
     $     rgsslp,    ! 0.5/(gravit*sslp)
     $     co2vmr,    ! CO2 volume mixing ratio
     $     dpfo3,     ! Voigt correction factor for O3
     $     dpfco2,    ! Voigt correction factor for CO2
     $     dayspy,    ! Number of days per 1 year
     $     pie        ! 3.14.....
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real pbr(plond,plev),           ! Prssr at mid-levels (dynes/cm2)
     $     pnm(plond,plevp),          ! Prssr at interfaces (dynes/cm2)
     $     co2em(plond,plevp),        ! Co2 emissivity function
     $     co2eml(plond,plev),        ! Co2 emissivity function
     $     tplnka(plond,plevp),       ! Planck fnctn level temperature
     $     s2c(plond,plevp),          ! H2o continuum path length
     $     s2t(plond,plevp),          ! H2o tmp and prs wghted path
     $     w(plond,plevp),            ! H2o prs wghted path
     $     h2otr(plond,plevp),        ! H2o trnsmssn fnct for o3 overlap
     $     plco2(plond,plevp),        ! Co2 prs wghted path length
     $     plh2o(plond,plevp),        ! H2o prs wfhted path length
     $     co2t(plond,plevp),         ! Tmp and prs wghted path length
     $     tint(plond,plevp),         ! Interface temperatures
     $     tlayr(plond,plevp),        ! K-1 level temperatures
     $     plol(plond,plevp),         ! Ozone prs wghted path length
     $     plos(plond,plevp)          ! Ozone path length
      real pmln(plond,plev),          ! Ln(pmidm1)
     $     piln(plond,plevp)          ! Ln(pintm1)
c
c   Trace gas variables
c
      real ucfc11(plond,plevp), ! CFC11 path length
     $     ucfc12(plond,plevp), ! CFC12 path length
     $     un2o0(plond,plevp),  ! N2O path length
     $     un2o1(plond,plevp),  ! N2O path length (hot band)
     $     uch4(plond,plevp),   ! CH4 path length
     $     uco211(plond,plevp), ! CO2 9.4 micron band path length
     $     uco212(plond,plevp), ! CO2 9.4 micron band path length
     $     uco213(plond,plevp), ! CO2 9.4 micron band path length
     $     uco221(plond,plevp), ! CO2 10.4 micron band path length
     $     uco222(plond,plevp), ! CO2 10.4 micron band path length
     $     uco223(plond,plevp), ! CO2 10.4 micron band path length
     $     uptype(plond,plevp), ! continuum path length
     $     bn2o0(plond,plevp),  ! pressure factor for n2o
     $     bn2o1(plond,plevp),  ! pressure factor for n2o
     $     bch4(plond,plevp)    ! pressure factor for ch4
      real abplnk1(14,plond,plevp), ! non-nearest layer Plack factor
     $     abplnk2(14,plond,plevp)  ! nearest layer factor
      real abstrc(plond)            ! total trace gas absorptivity
      real bplnk(14,plond,4) ! Planck functions for sub-divided layers
C
C Output arguments
C
      real abstot(plond,plevp,plevp), ! Total absorptivity
     $     absnxt(plond,plev,4)       ! Total nearest layer absorptivity
C
C---------------------------Local variables-----------------------------
C
      integer i,            ! Longitude index
     $        k,            ! Level index
     $        k1,           ! Level index
     $        k2,           ! Level index
     $        kn,           ! Nearest level index
     $        iband         ! Band  index
C
      real pnew(plond),     ! Effective pressure for H2O vapor linewidth
     $     trline(plond,2), ! Transmission due to H2O lines in window
     $     u(plond),        ! Pressure weighted H2O path length
     $     tbar(plond,4),   ! Mean layer temperature
     $     emm(plond,4),    ! Mean co2 emissivity
     $     o3emm(plond,4),  ! Mean o3 emissivity
     $     o3bndi,          ! Ozone band parameter
     $     temh2o(plond,4), ! Mean layer temperature equivalent to tbar
     $     k21,             ! Exponential coefficient used to calculate
C                           !  rotation band transmissvty in the 650-800
C                           !  cm-1 region (tr1)
     $     k22,             ! Exponential coefficient used to calculate
C                           !  rotation band transmissvty in the 500-650
C                           !  cm-1 region (tr2)
     $     uc1(plond)       ! H2o continuum pathlength in 500-800 cm-1
      real to3h2o(plond),   ! H2o trnsmsn for overlap with o3
     $     pi,              ! For co2 absorptivity computation
     $     sqti(plond),     ! Used to store sqrt of mean temperature
     $     et,              ! Co2 hot band factor
     $     et2,             ! Co2 hot band factor squared
     $     et4,             ! Co2 hot band factor to fourth power
     $     omet,            ! Co2 stimulated emission term
     $     f1co2,           ! Co2 central band factor
     $     f2co2(plond),    ! Co2 weak band factor
     $     f3co2(plond),    ! Co2 weak band factor
     $     t1co2(plond),    ! Overlap factr weak bands on strong band
     $     sqwp,            ! Sqrt of co2 pathlength
     $     f1sqwp(plond)    ! Main co2 band factor
      real oneme,           ! Co2 stimulated emission term
     $     alphat,          ! Part of the co2 stimulated emission term
     $     wco2,            ! Constants used to define co2 pathlength
     $     posqt,           ! Effective pressure for co2 line width
     $     u7,              ! Co2 hot band path length
     $     u8,              ! Co2 hot band path length
     $     u9,              ! Co2 hot band path length
     $     u13,             ! Co2 hot band path length
     $     rbeta7,          ! Inverse of co2 hot band line width par
     $     rbeta8,          ! Inverse of co2 hot band line width par
     $     rbeta9,          ! Inverse of co2 hot band line width par
     $     rbeta13          ! Inverse of co2 hot band line width par
      real tpatha(plond),     ! For absorptivity computation
     $     a,                 ! Eq(2) in table A3a of R&D
     $     abso(plond,6),     ! Absorptivity for various gases/bands
     $     dtp(plond),        ! Path temp minus 300 K used in h2o
C                             !  rotation band absorptivity
     $     dtx(plond),        ! Planck temperature minus 250 K
     $     dty(plond),        ! Path temperature minus 250 K
     $     dtz(plond),        ! Planck temperature minus 300 K
     $     term1(plond,4),    ! Equation(5) in table A3a of R&D(1986)
     $     term2(plond,4)     ! Delta a(Te) in table A3a of R&D(1986)
      real term3(plond,4),    ! DB/dT function for rotation and
C                             !  vibration-rotation band absorptivity
     $     term4(plond,4),    ! Equation(6) in table A3a of R&D(1986)
     $     term5(plond,4),    ! Delta a(Tp) in table A3a of R&D(1986)
     $     term6(plond,plevp),! DB/dT function for window region
     $     term7(plond,2),    ! Kl_inf(i) in eq(8) of table A3a of R&D
     $     term8(plond,2),    ! Delta kl_inf(i) in eq(8)
     $     term9(plond,plevp),! DB/dT function for 500-800 cm-1 region
     $     tr1,               ! Eqn(6) in table A2 of R&D for 650-800
     $     tr10(plond),       ! Eqn (6) times eq(4) in table A2
C                             !  of R&D for 500-650 cm-1 region
     $     tr2                ! Eqn(6) in table A2 of R&D for 500-650
      real tr5,               ! Eqn(4) in table A2 of R&D for 650-800
     $     tr6,               ! Eqn(4) in table A2 of R&D for 500-650
     $     tr9(plond),        ! Equation (6) times eq(4) in table A2
C                             !  of R&D for 650-800 cm-1 region
     $     uc(plond)          ! Y + 0.002U in eq(8) of table A2 of R&D
      real sqrtu(plond),      ! Sqrt of pressure weighted h20 pathlength
     $     fwk(plond),        ! Equation(33) in R&D far wing correction
     $     fwku(plond),       ! GU term in eqs(1) and (6) in table A2
     $     r2st(2),           ! 1/(2*beta) in eq(10) in table A2
     $     dtyp15(plond),     ! DeltaTp in eqs(11) & (12) in table A3a
     $     dtyp15sq(plond),   ! (DeltaTp)^2 in eqs(11) & (12) table A3a
     $     to3co2(plond),     ! P weighted temp in ozone band model
     $     dpnm(plond),       ! Pressure difference between two levels
     $     pnmsq(plond,plevp),! Pressure squared
     $     dw(plond),         ! Amount of h2o between two levels
     $     uinpl(plond,4),    ! Nearest layer subdivision factor
     $     winpl(plond,4),    ! Nearest layer subdivision factor
     $     zinpl(plond,4),    ! Nearest layer subdivision factor
     $     pinpl(plond,4),    ! Nearest layer subdivision factor
     $     dplh2o(plond)      ! Difference in press weighted h2o amount
      real r80257,            ! Conversion factor for h2o pathlength
     $     r293,              ! 1/293
     $     r250,              ! 1/250
     $     r3205,             ! Line width factor for o3 (see R&Di)
     $     r300,              ! 1/300
     $     rsslp,             ! Reciprocal of sea level pressure
     $     r2sslp             ! 1/2 of rsslp
      real  ds2c,     ! Y in eq(7) in table A2 of R&D
     $      a11,      ! A1 in table A3b for rotation band absorptivity
     $      a31,      ! A3 in table A3b for rotation band absorptivity
     $      a21,      ! First part in numerator of A2 in table A3b
     $      a22,      ! Second part in numerator of A2 in table A3b
     $      a23,      ! Denominator of A2 in table A3b (rotation band)
     $      t1t4,     ! Eq(3) in table A3a of R&D
     $      t2t5,     ! Eq(4) in table A3a of R&D
     $      rsum,     ! Eq(1) in table A2 of R&D
     $      a41,      ! Numerator in A2 in Vib-rot abstivity(table A3b)
     $      a51,      ! Denominator in A2 in Vib-rot (table A3b)
     $      a61       ! A3 factor for Vib-rot band in table A3b
      real  phi,      ! Eq(11) in table A3a of R&D
     $      psi,      ! Eq(12) in table A3a of R&D
     $      cf812,    ! Eq(11) in table A2 of R&D
     $      ubar,     ! H2o scaled path see comment for eq(10) table A2
     $      pbar,     ! H2o scaled pres see comment for eq(10) table A2
     $      g4        ! Arguement in exp() in eq(10) table A2
C
      real  dplos,    ! Ozone pathlength eq(A2) in R&Di
     $      dplol,    ! Presure weighted ozone pathlength
     $      tlocal,   ! Local interface temperature
     $      beta,     ! Ozone mean line parameter eq(A3) in R&Di
C                       (includes Voigt line correction factor)
     $      rphat,    ! Effective pressure for ozone beta
     $      tcrfac,   ! Ozone temperature factor table 1 R&Di
     $      tmp1,     ! Ozone band factor see eq(A1) in R&Di
     $      u1,       ! Effective ozone pathlength eq(A2) in R&Di
     $      realnu,   ! 1/beta factor in ozone band model eq(A1)
     $      tmp2,     ! Ozone band factor see eq(A1) in R&Di
     $      u2,       ! Effective ozone pathlength eq(A2) in R&Di
     $      rsqti     ! Reciprocal of sqrt of path temperature
C
      real  tpath,    ! Path temperature used in co2 band model
     $      tmp3,     ! Weak band factor see K&B
     $      rdpnmsq,  ! Reciprocal of difference in press^2
     $      rdpnm,    ! Reciprocal of difference in press
     $      p1,       ! Mean pressure factor
     $      p2,       ! Mean pressure factor
     $      dtym10,   ! T - 260 used in eq(9) and (10) table A3a
     $      dplco2,   ! Co2 pathlength
     $      corfac,   ! Correction factors in table A3b
     $      g2,       ! Part of arguement in eq(10) in table A2
     $      te,       ! A_0 T factor in ozone model table 1 of R&Di
     $      denom     ! Denominator in eq(8) of table A3a of R&D
c
c
      real th2o(plond),
     $     tco2(plond),
     $     to3(plond)
      integer wvl
C
C
C Transmission terms for various spectral intervals:
C
      real trab1(plond),  ! H2o     0 -  800 cm-1
     $     trab2(plond),  ! H2o   500 -  800 cm-1
     $     trab3(plond),  ! Co2   band system
     $     trab4(plond),  ! H2o   800 - 1000 cm-1
     $     trab5(plond),  ! 9.6 micrometer band
     $     trab6(plond),  ! H2o  1000 - 1200 cm-1
     $     trab7(plond)   ! H2o  1200 - 2200 cm-1
C
      real bndfct, ! Band absorptance parameter for co2
     $     absbnd  ! Proportional to co2 band absorptance
C
      real dbvtit(plond,plevp),      ! Intrfc drvtv plnck fnctn for o3
     $     dbvtly(plond,plev)        ! Level drvtv plnck fnctn for o3
C
C--------------------------Statement function---------------------------
C
      real dbvt,t     ! Planck fnctn tmp derivative for o3
C
      dbvt(t)=(-2.8911366682e-4+(2.3771251896e-6+1.1305188929e-10*t)*t)/
     $  (1.0+(-6.1364820707e-3+1.5550319767e-5*t)*t)
C
C-----------------------------------------------------------------------
C
C Initialize
C
      do k=1,plev
         do i=1,plon
            dbvtly(i,k) = dbvt(tlayr(i,k+1))
            dbvtit(i,k) = dbvt(tint(i,k))
         end do
      end do
      do i=1,plon
         dbvtit(i,plevp) = dbvt(tint(i,plevp))
      end do
C
      r80257  = 1./8.0257e-04
      r293    = 1./293.
      r250    = 1./250.
      r3205   = 1./.3205
      r300    = 1./300.
      rsslp   = 1./sslp
      r2sslp  = 1./(2.*sslp)
      r2st(1) = 1./(2.*st(1))
      r2st(2) = 1./(2.*st(2))
      bndfct  = 2.0*22.18/(sqrt(196.)*300.)
C
C Non-adjacent layer absorptivity:
C
C abso(i,1)     0 -  800 cm-1   h2o rotation band
C abso(i,2)  1200 - 2200 cm-1   h2o vibration-rotation band
C abso(i,3)   800 - 1200 cm-1   h2o window
C abso(i,4)   500 -  800 cm-1   h2o rotation band overlap with co2
C abso(i,5)   o3  9.6 micrometer band (nu3 and nu1 bands)
C abso(i,6)   co2 15  micrometer band system
C
      do k=1,plevp
         do i=1,plon
            pnmsq(i,k) = pnm(i,k)**2
            dtx(i) = tplnka(i,k) - 250.
            term6(i,k) = coeff(1,2) + coeff(2,2)*dtx(i)*
     $                   (1. +  c9*dtx(i)*(1. + c11*dtx(i)*
     $                   (1. + c13*dtx(i)*(1. + c15*dtx(i)))))
            term9(i,k) = coefi(1,2) + coefi(2,2)*dtx(i)*
     $                   (1. + c19*dtx(i)*(1. + c21*dtx(i)*
     $                   (1. + c23*dtx(i)*(1. + c25*dtx(i)))))
         end do
      end do
C
C Non-nearest layer level loops
C
      do 200 k1=plevp,1,-1
         do 100 k2=plevp,1,-1
            if(k1.eq.k2) go to 100
            do i=1,plon
               dplh2o(i) = plh2o(i,k1) - plh2o(i,k2)
               u(i)      = abs(dplh2o(i))
               sqrtu(i)  = sqrt(u(i))
               ds2c      = abs(s2c(i,k1) - s2c(i,k2))
               dw(i)     = abs(w(i,k1) - w(i,k2))
               uc1(i)    = (ds2c + 1.7e-3*u(i))*(1. +  2.*ds2c)/
     $                                          (1. + 15.*ds2c)
               uc(i)     = ds2c + 2.e-3*u(i)
            end do
            do i=1,plon
               pnew(i)   = u(i)/dw(i)
               tpatha(i) = (s2t(i,k1) - s2t(i,k2))/dplh2o(i)
               dtx(i)      = tplnka(i,k2) - 250.
               dty(i)      = tpatha(i)    - 250.
               dtyp15(i)   = dty(i) + 15.
               dtyp15sq(i) = dtyp15(i)**2
               dtz(i)      = dtx(i) - 50.
               dtp(i)      = dty(i) - 50.
            end do
            do iband=2,4,2
               do i=1,plon
                  term1(i,iband) = coefe(1,iband) + coefe(2,iband)*
     $                             dtx(i)*(1. + c1(iband)*dtx(i))
                  term2(i,iband) = coefb(1,iband) + coefb(2,iband)*
     $                             dtx(i)*(1. + c2(iband)*dtx(i)*
     $                                     (1. + c3(iband)*dtx(i)))
                  term3(i,iband) = coefd(1,iband) + coefd(2,iband)*
     $                             dtx(i)*(1. + c4(iband)*dtx(i)*
     $                                     (1. + c5(iband)*dtx(i)))
                  term4(i,iband) = coefa(1,iband) + coefa(2,iband)*
     $                             dty(i)*(1. + c6(iband)*dty(i))
                  term5(i,iband) = coefc(1,iband) + coefc(2,iband)*
     $                             dty(i)*(1. + c7(iband)*dty(i))
               end do
            end do
C
C abso(i,1)     0 -  800 cm-1   h2o rotation band
C
            do i=1,plon
               a11 = 0.44 + 3.380e-4*dtz(i) - 1.520e-6*dtz(i)*dtz(i)
               a31 = 1.05 - 6.000e-3*dtp(i) + 3.000e-6*dtp(i)*dtp(i)
               a21 = 1.00 + 1.717e-3*dtz(i) - 1.133e-5*dtz(i)*dtz(i)
               a22 = 1.00 + 4.443e-3*dtp(i) + 2.750e-5*dtp(i)*dtp(i)
               a23 = 1.00 + 3.600*sqrtu(i)
               corfac  = a31*(a11 + ((2.*a21*a22)/a23))
               t1t4    = term1(i,2)*term4(i,2)
               t2t5    = term2(i,2)*term5(i,2)
               a       = t1t4 + t2t5/(1. + t2t5*sqrtu(i)*corfac)
               fwk(i)  = fwcoef + fwc1/(1. + fwc2*u(i))
               fwku(i) = fwk(i)*u(i)
               rsum    = exp(-a*(sqrtu(i) + fwku(i)))
               abso(i,1) = (1. - rsum)*term3(i,2)
C               trab1(i)  = rsum
            end do
C
C abso(i,2)  1200 - 2200 cm-1   h2o vibration-rotation band
C
            do i=1,plon
               a41   = 1.75 - 3.960e-03*dtz(i)
               a51   = 1.00 + 1.3*sqrtu(i)
               a61   = 1.00 + 1.250e-03*dtp(i) + 6.250e-05*dtp(i)*dtp(i)
               corfac = .29*(1. + a41/a51)*a61
               t1t4   = term1(i,4)*term4(i,4)
               t2t5   = term2(i,4)*term5(i,4)
               a      = t1t4 + t2t5/(1. + t2t5*sqrtu(i)*corfac)
               rsum   = exp(-a*(sqrtu(i) + fwku(i)))
               abso(i,2) = (1. - rsum)*term3(i,4)
C               trab7(i)  = rsum
            end do
C
C Line transmission in 800-1000 and 1000-1200 cm-1 intervals
C
            do k=1,2
               do i=1,plon
                  phi   = exp(a1(k)*dtyp15(i) + a2(k)*dtyp15sq(i))
                  psi   = exp(b1(k)*dtyp15(i) + b2(k)*dtyp15sq(i))
                  ubar  = dw(i)*phi*1.66*r80257
                  pbar  = pnew(i)*(psi/phi)
                  cf812 = cfa1 + (1. - cfa1)/(1. + ubar*pbar*10.)
                  g2    = 1. + ubar*4.0*st(k)*cf812/pbar
                  g4    = realk(k)*pbar*r2st(k)*(sqrt(g2) - 1.)
                  trline(i,k) = exp(-g4)
               end do
            end do
            do i=1,plon
               term7(i,1) = coefj(1,1) + coefj(2,1)*dty(i)*
     $                                   (1. + c16*dty(i))
               term8(i,1) = coefk(1,1) + coefk(2,1)*dty(i)*
     $                                   (1. + c17*dty(i))
               term7(i,2) = coefj(1,2) + coefj(2,2)*dty(i)*
     $                                   (1. + c26*dty(i))
               term8(i,2) = coefk(1,2) + coefk(2,2)*dty(i)*
     $                                   (1. + c27*dty(i))
            end do
C
C abso(i,3)   800 - 1200 cm-1   h2o window
C abso(i,4)   500 -  800 cm-1   h2o rotation band overlap with co2
C
            do i=1,plon
               k21    = term7(i,1) + term8(i,1)/
     $             (1. + (c30 + c31*(dty(i)-10.)*(dty(i)-10.))*sqrtu(i))
               k22    = term7(i,2) + term8(i,2)/
     $             (1. + (c28 + c29*(dty(i)-10.))*sqrtu(i))
               tr1    = exp(-(k21*(sqrtu(i) + fc1*fwku(i))))
               tr2    = exp(-(k22*(sqrtu(i) + fc1*fwku(i))))
               tr5    = exp(-((coefh(1,3) + coefh(2,3)*dtx(i))*uc1(i)))
               tr6    = exp(-((coefh(1,4) + coefh(2,4)*dtx(i))*uc1(i)))
               tr9(i)   = tr1*tr5
               tr10(i)  = tr2*tr6
               th2o(i) = tr10(i)
               trab2(i) = 0.65*tr9(i) + 0.35*tr10(i)
               trab4(i) = exp(-(coefg(1,3) + coefg(2,3)*dtx(i))*uc(i))
               trab6(i) = exp(-(coefg(1,4) + coefg(2,4)*dtx(i))*uc(i))
               abso(i,3) = term6(i,k2)*(1. - .5*trab4(i)*trline(i,2) -
     $                                       .5*trab6(i)*trline(i,1))
               abso(i,4) = term9(i,k2)*.5*(tr1 - tr9(i) + tr2 - tr10(i))
            end do
            if(k2.lt.k1) then
               do i=1,plon
                  to3h2o(i) = h2otr(i,k1)/h2otr(i,k2)
               end do
            else
               do i=1,plon
                  to3h2o(i) = h2otr(i,k2)/h2otr(i,k1)
               end do
            end if
C
C abso(i,5)   o3  9.6 micrometer band (nu3 and nu1 bands)
C
            do i=1,plon
               dpnm(i)  = pnm(i,k1) - pnm(i,k2)
               to3co2(i)=(pnm(i,k1)*co2t(i,k1) - pnm(i,k2)*co2t(i,k2))/
     $                   dpnm(i)
               te       = (to3co2(i)*r293)**.7
               dplos    = plos(i,k1) - plos(i,k2)
               dplol    = plol(i,k1) - plol(i,k2)
               u1       = 18.29*abs(dplos)/te
               u2       = .5649*abs(dplos)/te
               rphat    = dplol/dplos
               tlocal   = tint(i,k2)
               tcrfac   = sqrt(tlocal*r250)*te
               beta     = r3205*(rphat + dpfo3*tcrfac)
               realnu   = te/beta
               tmp1     = u1/sqrt(4. + u1*(1. + realnu))
               tmp2     = u2/sqrt(4. + u2*(1. + realnu))
               o3bndi    = 74.*te*alog(1. + tmp1 + tmp2)
               abso(i,5) = o3bndi*to3h2o(i)*dbvtit(i,k2)
               to3(i)   = 1.0/(1. + 0.1*tmp1 + 0.1*tmp2)
C               trab5(i)  = 1.-(o3bndi/(1060-980.))
            end do
C
C abso(i,6)      co2 15  micrometer band system
C
            do i=1,plon
               sqwp      = sqrt(abs(plco2(i,k1) - plco2(i,k2)))
               et        = exp(-480./to3co2(i))
               sqti(i)   = sqrt(to3co2(i))
               rsqti     = 1./sqti(i)
               et2       = et*et
               et4       = et2*et2
               omet      = 1. - 1.5*et2
               f1co2     = 899.70*omet*
     $                   (1. + 1.94774*et + 4.73486*et2)*rsqti
               f1sqwp(i) = f1co2*sqwp
               t1co2(i)  = 1./(1. + (245.18*omet*sqwp*rsqti))
               oneme     = 1. - et2
               alphat    = oneme**3*rsqti
               pi        = abs(dpnm(i))
               wco2      =  2.5221*co2vmr*pi*rga
               u7        =  4.9411e4*alphat*et2*wco2
               u8        =  3.9744e4*alphat*et4*wco2
               u9        =  1.0447e5*alphat*et4*et2*wco2
               u13       = 2.8388e3*alphat*et4*wco2
               tpath     = to3co2(i)
               tlocal    = tint(i,k2)
               tcrfac    = sqrt(tlocal*r250*tpath*r300)
               posqt     = ((pnm(i,k2) + pnm(i,k1))*r2sslp +
     $                     dpfco2*tcrfac)*rsqti
               rbeta7    = 1./(5.3228*posqt)
               rbeta8    = 1./(10.6576*posqt)
               rbeta9    = rbeta7
               rbeta13   = rbeta9
               f2co2(i)  = (u7/sqrt(4. + u7*(1. + rbeta7))) +
     $                     (u8/sqrt(4. + u8*(1. + rbeta8))) +
     $                     (u9/sqrt(4. + u9*(1. + rbeta9)))
               f3co2(i)  = u13/sqrt(4. + u13*(1. + rbeta13))
            end do
            if (k2.ge.k1) then
               do i=1,plon
                  sqti(i) = sqrt(tlayr(i,k2))
               end do
            end if
C
            do i=1,plon
               tmp1      = alog(1. + f1sqwp(i))
               tmp2      = alog(1. + f2co2(i))
               tmp3      = alog(1. + f3co2(i))
               absbnd    = (tmp1 + 2.*t1co2(i)*tmp2 + 2.*tmp3)*sqti(i)
               abso(i,6) = trab2(i)*co2em(i,k2)*absbnd
               tco2(i)=1./(1.0+10.0*(u7/sqrt(4. + u7*(1. + rbeta7))))
C               trab3(i)  = 1. - bndfct*absbnd
            end do
c
c     Calculate absorptivity due to trace gases
c
            call trcab(k1, k2, ucfc11, ucfc12, un2o0,  un2o1,
     $                         uch4,   uco211, uco212, uco213, 
     $                         uco221, uco222, uco223, bn2o0, 
     $                         bn2o1,  bch4,   to3co2, pnm,
     $                         dw,     pnew,   s2c,    uptype,
     $                         u,      abplnk1,tco2,   th2o,   
     $                         to3,    abstrc)
C
C Sum total absorptivity
C
            do i=1,plon
               abstot(i,k1,k2) = abso(i,1) + abso(i,2) + abso(i,3) +
     $                           abso(i,4) + abso(i,5) + abso(i,6)
     $                           + abstrc(i)
            end do
  100    continue
  200 continue        ! End of non-nearest layer level loops
C
C Non-adjacent layer absorptivity:
C
C abso(i,1)     0 -  800 cm-1   h2o rotation band
C abso(i,2)  1200 - 2200 cm-1   h2o vibration-rotation band
C abso(i,3)   800 - 1200 cm-1   h2o window
C abso(i,4)   500 -  800 cm-1   h2o rotation band overlap with co2
C abso(i,5)   o3  9.6 micrometer band (nu3 and nu1 bands)
C abso(i,6)   co2 15  micrometer band system
C
C Nearest layer level loop
C
      do 500 k2=plev,1,-1
         do i=1,plon
            tbar(i,1)   = 0.5*(tint(i,k2+1) + tlayr(i,k2+1))
            emm(i,1)    = 0.5*(co2em(i,k2+1) + co2eml(i,k2))
            tbar(i,2)   = 0.5*(tlayr(i,k2+1) + tint(i,k2))
            emm(i,2)    = 0.5*(co2em(i,k2) + co2eml(i,k2))
            tbar(i,3)   = 0.5*(tbar(i,2) + tbar(i,1))
            emm(i,3)    = emm(i,1)
            tbar(i,4)   = tbar(i,3)
            emm(i,4)    = emm(i,2)
            o3emm(i,1)  = 0.5*(dbvtit(i,k2+1) + dbvtly(i,k2))
            o3emm(i,2)  = 0.5*(dbvtit(i,k2) + dbvtly(i,k2))
            o3emm(i,3)  = o3emm(i,1)
            o3emm(i,4)  = o3emm(i,2)
            temh2o(i,1) = tbar(i,1)
            temh2o(i,2) = tbar(i,2)
            temh2o(i,3) = tbar(i,1)
            temh2o(i,4) = tbar(i,2)
            dpnm(i)     = pnm(i,k2+1) - pnm(i,k2)
         end do
c----------------------------------------------------------
c  Weighted Planck functions for trace gases
c
      do wvl = 1,14
         do i = 1,plon
            bplnk(wvl,i,1) = 0.5*(abplnk1(wvl,i,k2+1) + 
     $                          abplnk2(wvl,i,k2))
            bplnk(wvl,i,2) = 0.5*(abplnk1(wvl,i,k2) + 
     $                          abplnk2(wvl,i,k2))
            bplnk(wvl,i,3) = bplnk(wvl,i,1)
            bplnk(wvl,i,4) = bplnk(wvl,i,2)
         end do
      end do
c---------------------------------------------------------
         do i=1,plon
            rdpnmsq    = 1./(pnmsq(i,k2+1) - pnmsq(i,k2))
            rdpnm      = 1./dpnm(i)
            p1         = .5*(pbr(i,k2) + pnm(i,k2+1))
            p2         = .5*(pbr(i,k2) + pnm(i,k2  ))
            uinpl(i,1) =  (pnmsq(i,k2+1) - p1**2)*rdpnmsq
            uinpl(i,2) = -(pnmsq(i,k2  ) - p2**2)*rdpnmsq
            uinpl(i,3) = -(pnmsq(i,k2  ) - p1**2)*rdpnmsq
            uinpl(i,4) =  (pnmsq(i,k2+1) - p2**2)*rdpnmsq
            winpl(i,1) = (.5*( pnm(i,k2+1) - pbr(i,k2)))*rdpnm
            winpl(i,2) = (.5*(-pnm(i,k2  ) + pbr(i,k2)))*rdpnm
            winpl(i,3) = (.5*( pnm(i,k2+1) + pbr(i,k2)) - pnm(i,k2  ))*
     $                   rdpnm
            winpl(i,4) = (.5*(-pnm(i,k2  ) - pbr(i,k2)) + pnm(i,k2+1))*
     $                   rdpnm
            tmp1       = 1./(piln(i,k2+1) - piln(i,k2))
            tmp2       = piln(i,k2+1) - pmln(i,k2)
            tmp3       = piln(i,k2  ) - pmln(i,k2)
            zinpl(i,1) = (.5*tmp2          )*tmp1
            zinpl(i,2) = (        - .5*tmp3)*tmp1
            zinpl(i,3) = (.5*tmp2 -    tmp3)*tmp1
            zinpl(i,4) = (   tmp2 - .5*tmp3)*tmp1
            pinpl(i,1) = 0.5*(p1 + pnm(i,k2+1))
            pinpl(i,2) = 0.5*(p2 + pnm(i,k2  ))
            pinpl(i,3) = 0.5*(p1 + pnm(i,k2  ))
            pinpl(i,4) = 0.5*(p2 + pnm(i,k2+1))
         end do
         do 400 kn=1,4
            do i=1,plon
               u(i)     = uinpl(i,kn)*abs(plh2o(i,k2) - plh2o(i,k2+1))
               sqrtu(i) = sqrt(u(i))
               dw(i)    = abs(w(i,k2) - w(i,k2+1))
               pnew(i)  = u(i)/(winpl(i,kn)*dw(i))
               ds2c     = abs(s2c(i,k2) - s2c(i,k2+1))
               uc1(i)   = uinpl(i,kn)*ds2c
               uc1(i)   = (uc1(i) + 1.7e-3*u(i))*(1. +  2.*uc1(i))/
     $                                           (1. + 15.*uc1(i))
               uc(i)    = uinpl(i,kn)*ds2c + 2.e-3*u(i)
            end do
            do i=1,plon
               dtx(i)      = temh2o(i,kn) - 250.
               dty(i)      = tbar(i,kn) - 250.
               dtyp15(i)   = dty(i) + 15.
               dtyp15sq(i) = dtyp15(i)**2
               dtz(i)      = dtx(i) - 50.
               dtp(i)      = dty(i) - 50.
            end do
            do iband=2,4,2
               do i=1,plon
                  term1(i,iband) = coefe(1,iband) + coefe(2,iband)*
     $                             dtx(i)*(1. + c1(iband)*dtx(i))
                  term2(i,iband) = coefb(1,iband) + coefb(2,iband)*
     $                             dtx(i)*(1. + c2(iband)*dtx(i)*
     $                                     (1. + c3(iband)*dtx(i)))
                  term3(i,iband) = coefd(1,iband) + coefd(2,iband)*
     $                             dtx(i)*(1. + c4(iband)*dtx(i)*
     $                                     (1. + c5(iband)*dtx(i)))
                  term4(i,iband) = coefa(1,iband) + coefa(2,iband)*
     $                             dty(i)*(1. + c6(iband)*dty(i))
                  term5(i,iband) = coefc(1,iband) + coefc(2,iband)*
     $                             dty(i)*(1. + c7(iband)*dty(i))
               end do
            end do
C
C abso(i,1)     0 -  800 cm-1   h2o rotation band
C
            do i=1,plon
               a11 = 0.44 + 3.380e-4*dtz(i) - 1.520e-6*dtz(i)*dtz(i)
               a31 = 1.05 - 6.000e-3*dtp(i) + 3.000e-6*dtp(i)*dtp(i)
               a21 = 1.00 + 1.717e-3*dtz(i) - 1.133e-5*dtz(i)*dtz(i)
               a22 = 1.00 + 4.443e-3*dtp(i) + 2.750e-5*dtp(i)*dtp(i)
               a23 = 1.00 + 3.600*sqrtu(i)
               corfac    = a31*(a11 + ((2.*a21*a22)/a23))
               t1t4      = term1(i,2)*term4(i,2)
               t2t5      = term2(i,2)*term5(i,2)
               a         = t1t4 + t2t5/(1. + t2t5*sqrtu(i)*corfac)
               fwk(i)    = fwcoef + fwc1/(1. + fwc2*u(i))
               fwku(i)   = fwk(i)*u(i)
               rsum      = exp(-a*(sqrtu(i) + fwku(i)))
               abso(i,1) = (1. - rsum)*term3(i,2)
C               trab1(i) = rsum
            end do
C
C abso(i,2)  1200 - 2200 cm-1   h2o vibration-rotation band
C
            do i=1,plon
               a41   = 1.75 - 3.960e-03*dtz(i)
               a51   = 1.00 + 1.3*sqrtu(i)
               a61   = 1.00 + 1.250e-03*dtp(i) + 6.250e-05*dtp(i)*dtp(i)
               corfac = .29*(1. + a41/a51)*a61
               t1t4   = term1(i,4)*term4(i,4)
               t2t5   = term2(i,4)*term5(i,4)
               a      = t1t4 + t2t5/(1. + t2t5*sqrtu(i)*corfac)
               rsum   = exp(-a*(sqrtu(i) + fwku(i)))
               abso(i,2) = (1. - rsum)*term3(i,4)
C               trab7(i) = rsum
            end do
C
C Line transmission in 800-1000 and 1000-1200 cm-1 intervals
C
            do k=1,2
               do i=1,plon
                  phi   = exp(a1(k)*dtyp15(i) + a2(k)*dtyp15sq(i))
                  psi   = exp(b1(k)*dtyp15(i) + b2(k)*dtyp15sq(i))
                  ubar  = dw(i)*phi*winpl(i,kn)*1.66*r80257
                  pbar  = pnew(i)*(psi/phi)
                  cf812 = cfa1 + (1. - cfa1)/(1. + ubar*pbar*10.)
                  g2    = 1. + ubar*4.0*st(k)*cf812/pbar
                  g4    = realk(k)*pbar*r2st(k)*(sqrt(g2) - 1.)
                  trline(i,k) = exp(-g4)
               end do
            end do
            do i=1,plon
               term7(i,1) = coefj(1,1) + coefj(2,1)*dty(i)*
     $                                   (1. + c16*dty(i))
               term8(i,1) = coefk(1,1) + coefk(2,1)*dty(i)*
     $                                   (1. + c17*dty(i))
               term7(i,2) = coefj(1,2) + coefj(2,2)*dty(i)*
     $                                   (1. + c26*dty(i))
               term8(i,2) = coefk(1,2) + coefk(2,2)*dty(i)*
     $                                   (1. + c27*dty(i))
            end do
C
C abso(i,3)   800 - 1200 cm-1   h2o window
C abso(i,4)   500 -  800 cm-1   h2o rotation band overlap with co2
C
            do i=1,plon
               dtym10     = dty(i) - 10.
               denom      = 1. + (c30 + c31*dtym10*dtym10)*sqrtu(i)
               k21        = term7(i,1) + term8(i,1)/denom
               denom      = 1. + (c28 + c29*dtym10       )*sqrtu(i)
               k22        = term7(i,2) + term8(i,2)/denom
               term9(i,2) = coefi(1,2) + coefi(2,2)*dtx(i)*
     $                     (1. + c19*dtx(i)*(1. + c21*dtx(i)*
     $                      (1. + c23*dtx(i)*(1. + c25*dtx(i)))))
               tr1     = exp(-(k21*(sqrtu(i) + fc1*fwku(i))))
               tr2     = exp(-(k22*(sqrtu(i) + fc1*fwku(i))))
               tr5     = exp(-((coefh(1,3) + coefh(2,3)*dtx(i))*uc1(i)))
               tr6     = exp(-((coefh(1,4) + coefh(2,4)*dtx(i))*uc1(i)))
               tr9(i)  = tr1*tr5
               tr10(i) = tr2*tr6
               trab2(i)= 0.65*tr9(i) + 0.35*tr10(i)
               th2o(i) = tr10(i)
               trab4(i)= exp(-(coefg(1,3) + coefg(2,3)*dtx(i))*uc(i))
               trab6(i)= exp(-(coefg(1,4) + coefg(2,4)*dtx(i))*uc(i))
               term6(i,2) = coeff(1,2) + coeff(2,2)*dtx(i)*
     $                     (1. + c9*dtx(i)*(1. + c11*dtx(i)*
     $                     (1. + c13*dtx(i)*(1. + c15*dtx(i)))))
               abso(i,3)  = term6(i,2)*(1. - .5*trab4(i)*trline(i,2) -
     $                                       .5*trab6(i)*trline(i,1))
               abso(i,4)  = term9(i,2)*.5*(tr1 - tr9(i) + tr2 - tr10(i))
            end do
C
C abso(i,5)  o3  9.6 micrometer (nu3 and nu1 bands)
C
            do i=1,plon
               te        = (tbar(i,kn)*r293)**.7
               dplos     = abs(plos(i,k2+1) - plos(i,k2))
               u1        = zinpl(i,kn)*18.29*dplos/te
               u2        = zinpl(i,kn)*.5649*dplos/te
               tlocal    = tbar(i,kn)
               tcrfac    = sqrt(tlocal*r250)*te
               beta      = r3205*(pinpl(i,kn)*rsslp + dpfo3*tcrfac)
               realnu    = te/beta
               tmp1      = u1/sqrt(4. + u1*(1. + realnu))
               tmp2      = u2/sqrt(4. + u2*(1. + realnu))
               o3bndi    = 74.*te*alog(1. + tmp1 + tmp2)
               abso(i,5) = o3bndi*o3emm(i,kn)*
     $                     (h2otr(i,k2+1)/h2otr(i,k2))
               to3(i)    = 1.0/(1. + 0.1*tmp1 + 0.1*tmp2)
C               trab5(i) = 1.-(o3bndi/(1060-980.))
            end do
C
C abso(i,6)   co2 15  micrometer band system
C
            do 300 i=1,plon
               dplco2   = plco2(i,k2+1) - plco2(i,k2)
               sqwp     = sqrt(uinpl(i,kn)*dplco2)
               et       = exp(-480./tbar(i,kn))
               sqti(i)  = sqrt(tbar(i,kn))
               rsqti    = 1./sqti(i)
               et2      = et*et
               et4      = et2*et2
               omet     = (1. - 1.5*et2)
               f1co2    = 899.70*omet*
     $                    (1. + 1.94774*et + 4.73486*et2)*rsqti
               f1sqwp(i)= f1co2*sqwp
               t1co2(i) = 1./(1. + (245.18*omet*sqwp*rsqti))
               oneme    = 1. - et2
               alphat   = oneme**3*rsqti
               pi       = abs(dpnm(i))*winpl(i,kn)
               wco2     = 2.5221*co2vmr*pi*rga
               u7       = 4.9411e4*alphat*et2*wco2
               u8       = 3.9744e4*alphat*et4*wco2
               u9       = 1.0447e5*alphat*et4*et2*wco2
               u13      = 2.8388e3*alphat*et4*wco2
               tpath    = tbar(i,kn)
               tlocal   = tbar(i,kn)
               tcrfac   = sqrt((tlocal*r250)*(tpath*r300))
               posqt    = (pinpl(i,kn)*rsslp + dpfco2*tcrfac)*rsqti
               rbeta7   = 1./(5.3228*posqt)
               rbeta8   = 1./(10.6576*posqt)
               rbeta9   = rbeta7
               rbeta13  = rbeta9
               f2co2(i) = u7/sqrt(4. + u7*(1. + rbeta7)) +
     $                    u8/sqrt(4. + u8*(1. + rbeta8)) +
     $                    u9/sqrt(4. + u9*(1. + rbeta9))
               f3co2(i) = u13/sqrt(4. + u13*(1. + rbeta13))
               tmp1     = alog(1. + f1sqwp(i))
               tmp2     = alog(1. + f2co2(i))
               tmp3     = alog(1. + f3co2(i))
               absbnd   = (tmp1 + 2.*t1co2(i)*tmp2 + 2.*tmp3)*sqti(i)
               abso(i,6)= trab2(i)*emm(i,kn)*absbnd
               tco2(i)=1.0/(1.0+ 10.0*u7/sqrt(4. + u7*(1. + rbeta7)))
C               trab3(i) = 1. - bndfct*absbnd
  300       continue
c
c   Calculate trace gas absorptivity for nearest layer
c
            call trcabn(k2, kn, ucfc11, ucfc12, un2o0,  un2o1,
     $                          uch4,   uco211, uco212, uco213, 
     $                          uco221, uco222, uco223, bn2o0, 
     $                          bn2o1,  bch4,   tbar,   bplnk, 
     $                          winpl,  pinpl,  tco2,   th2o,
     $                          to3,    uptype, dw,     s2c, 
     $                          u,      pnew,   abstrc)
C
C Total next layer absorptivity:
C
            do i=1,plon
               absnxt(i,k2,kn) = abso(i,1) + abso(i,2) + abso(i,3) +
     $                           abso(i,4) + abso(i,5) + abso(i,6)
     $                           + abstrc(i)
            end do
  400    continue
  500 continue                  !  end of nearest layer level loop
C
      return
      end
      subroutine radclr(coszrs  ,trayoslp,pflx    ,abh2o   ,abo3    ,
     $                  abco2   ,abo2    ,uth2o   ,uto3    ,utco2   ,
     $                  uto2    ,tauaer  ,waer    ,gaer    ,faer    ,
     $                  nloop   ,is      ,ie      ,rdir    ,rdif    ,
     $                  tdir    ,tdif    ,explay  ,exptdn  ,rdndif  ,
     $                  tottrn  )
C-----------------------------------------------------------------------
C
C Delta-Eddington solution for special clear sky computation
C
C Computes total reflectivities and transmissivities for two atmospheric
C layers: an overlying purely ozone absorbing layer, and the rest of the
C column below.
C
C For more details , see Briegleb, Bruce P., 1992: Delta-Eddington
C Approximation for Solar Radiation in the NCAR Community Climate Model,
C Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
C
C---------------------------Code history--------------------------------
C
C Original version:  B. Briegleb
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Kiehl, B. Briegleb, August 1992
C
C-----------------------------------------------------------------------
c
c $Id: radclr.F,v 1.2 1995/03/17 18:54:08 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
C------------------------------Parameters-------------------------------
c
c $Id: prgrid.h,v 1.2 1995/02/10 01:09:10 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation resolution and I/O parameters
C
      integer
     $     plon,    ! number of longitudes
     $     plev,    ! number of vertical levels
     $     plat,    ! number of latitudes
     $     pcnst,   ! number of constituents (including water vapor)
     $     plevmx,  ! number of subsurface levels
     $     plevp,   ! plev + 1
     $     nxpt,    ! no.of points outside active domain for interpolant
     $     jintmx,  ! number of extra latitudes in polar region
     $     plond,   ! slt extended domain longitude
     $     platd,   ! slt extended domain lat.
     $     plevd    ! fold plev,pcnst indices into one
      integer
     $     plngbuf  ! length of absorptivity/emissivity record
C
      parameter(plon    = 1,
     $          plev    = 18,
     $          plat    = 1,
     $          pcnst   = 1,
     $          plevmx  = 4,
     $          plevp   = plev + 1,
     $          nxpt    = 1,
     $          jintmx  = 1,
     $          plond   = plon + 1 + 2*nxpt,
     $          platd   = plat + 2*nxpt + 2*jintmx,
     $          plevd   = plev*(3 + pcnst))
      parameter(plngbuf = 512*((plond*plevp*plevp + plond*plev*4 +
     $                          plond*plevp)/512 + 1))
C
C-----------------------------------------------------------------------
C
C Minimum total transmission below which no layer computation are done:
C
      real  trmin,          ! Minimum total transmission allowed
     $      wray,           ! Rayleigh single scatter albedo
     $      gray,           ! Rayleigh asymetry parameter
     $      fray            ! Rayleigh forward scattered fraction
      parameter (trmin = 1.e-3,
     $           wray = 0.999999,
     $           gray = 0.0,
     $           fray = 0.1)
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real coszrs(plond),         ! Cosine zenith angle
     $     trayoslp,              ! Tray/sslp
     $     pflx(plond,0:plevp),   ! Interface pressure
     $     abh2o,                 ! Absorption coefficiant for h2o
     $     abo3 ,                 ! Absorption coefficiant for o3
     $     abco2,                 ! Absorption coefficiant for co2
     $     abo2 ,                 ! Absorption coefficiant for o2
     $     uth2o(plond),          ! Total column absorber amount of h2o
     $     uto3(plond),           ! Total column absorber amount of  o3
     $     utco2(plond),          ! Total column absorber amount of co2
     $     uto2(plond)            ! Total column absorber amount of  o2
      real tauaer(plond),         ! Total column aerosol extinction
     $     waer(plond),           ! Aerosol single scattering albedo
     $     gaer(plond),           ! Aerosol asymmetry parameter
     $     faer(plond)            ! Aerosol forward scattering fraction
      integer nloop,              ! Number of loops (1 or 2)
     $        is(2),              ! Starting index for 1 or 2 loops
     $        ie(2)               ! Ending index for 1 or 2 loops
C
C Input/Output arguments
C
C Following variables are defined for each layer; note, we use layer 0 
C to refer to the entire atmospheric column:
C
      real rdir(plond,0:plev),    ! Layer reflectivity to direct rad
     $     rdif(plond,0:plev),    ! Layer refflectivity to diffuse rad
     $     tdir(plond,0:plev),    ! Layer transmission to direct rad
     $     tdif(plond,0:plev),    ! Layer transmission to diffuse rad
     $     explay(plond,0:plev)   ! Solar beam exp transmn for layer
C
C Note that the following variables are defined on interfaces, with
C the index k referring to the top interface of the kth layer:
C exptdn,rdndif,tottrn; for example, tottrn(k=5) refers to the total
C transmission to the top interface of the 5th layer.
C
      real exptdn(plond,0:plevp), ! Solar beam exp down transmn from top
     $     rdndif(plond,0:plevp), ! Added dif ref for layers above
     $     tottrn(plond,0:plevp)  ! Total transmission for layers above
C
      external  resetr,     ! Resets array elements to zero
     $          whenfgt     ! Collect indices for greater than condition
C
C---------------------------Local variables-----------------------------
C
      integer i,            ! Longitude index
     $        k,            ! Level index
     $        nn,           ! Index of longitude loops (max=nloop)
     $        ii,           ! Longitude index
     $        nval,         ! Number of long values satisfying criteria
     $        index(plond)  ! Array of longitude indices
C
      real taugab(plond),   ! Total column gas absorption optical depth
     $     tauray(plond),   ! Column rayleigh optical depth
     $     tautot       ,   ! Total column optical depth
     $       wtot       ,   ! Total column single scatter albedo
     $       gtot       ,   ! Total column asymmetry parameter
     $       ftot           ! Total column forward scatter fraction
      real   ts,            ! Column scaled extinction optical depth
     $       ws,            ! Column scaled single scattering albedo
     $       gs             ! Column scaled asymmetry parameter
      real rdenom,          ! Mulitiple scattering term
     $     rdirexp,         ! Layer direct ref times exp transmission
     $     tdnmexp          ! Total transmission minus exp transmission
C
C---------------------------Statement functions-------------------------
C
C Statement functions for delta-Eddington solution; for detailed
C explanation of individual terms, see the routine 'radded'.
C
      real alpha,gamma,el,taus,omgs,asys,u,n,lm,ne
      real w,uu,g,e,f,t,et
C
C Intermediate terms for delta-Eddington solution
C
      real alp,gam,ue,arg,extins,amg,apg
C
      alpha(w,uu,g,e) = .75*w*uu*((1. + g*(1-w))/(1. - e*e*uu*uu))
      gamma(w,uu,g,e) = .50*w*((3.*g*(1.-w)*uu*uu + 1.)/(1.-e*e*uu*uu))
      el(w,g)         = sqrt(3.*(1-w)*(1. - w*g))
      taus(w,f,t)     = (1. - w*f)*t
      omgs(w,f)       = (1. - f)*w/(1. - w*f)
      asys(g,f)       = (g - f)/(1. - f)
      u(w,g,e)        = 1.5*(1. - w*g)/e
      n(uu,et)        = ((uu+1.)*(uu+1.)/et ) - ((uu-1.)*(uu-1.)*et)
C
C-----------------------------------------------------------------------
C
C Initialize all total transmimission values to 0, so that nighttime 
C values from previous computations are not used:
C
      call resetr(tottrn,plond*2,0.)
C
C Compute total direct beam transmission, total transmission, and
C reflectivity for diffuse radiation (from below) for all layers
C above each interface by starting from the top and adding layers
C down:
C
C The top layer is assumed to be a purely absorbing ozone layer, and
C that the mean diffusivity for diffuse transmission is 1.66:
C
      do nn=1,nloop
         do i=is(nn),ie(nn)
C
            taugab(i) = abo3*uto3(i)
C
C Limit argument of exponential to 25, in case coszrs is very small:
C
            arg         = amin1(taugab(i)/coszrs(i),25.)
            explay(i,0) = exp(-arg)
            tdir(i,0)   = explay(i,0)
C
C Same limit for diffuse transmission:
C
            arg         = amin1(1.66*taugab(i),25.)
            tdif(i,0)   = exp(-arg)
C
            rdir(i,0)   = 0.0
            rdif(i,0)   = 0.0
C
C Initialize top interface of extra layer:
C
            exptdn(i,0) =   1.0
            rdndif(i,0) =   0.0
            tottrn(i,0) =   1.0
C
            rdndif(i,1) = rdif(i,0)
            tottrn(i,1) = tdir(i,0)
C
         end do
      end do
C
C Now, complete the rest of the column; if the total transmission
C through the top ozone layer is less than trmin, then no
C delta-Eddington computation for the underlying column is done:
C
      do 200 k=1,1
C
C Initialize current layer properties to zero;only if total transmission
C to the top interface of the current layer exceeds the minimum, will
C these values be computed below:
C
         do nn=1,nloop
            do i=is(nn),ie(nn)
C
               rdir(i,k)   =  0.0
               rdif(i,k)   =  0.0
               tdir(i,k)   =  0.0
               tdif(i,k)   =  0.0
               explay(i,k) =  0.0
C
C Calculates the solar beam transmission, total transmission, and
C reflectivity for diffuse radiation from below at the top of the
C current layer:
C
               exptdn(i,k) = exptdn(i,k-1)*explay(i,k-1)
               rdenom      = 1./(1. - rdif(i,k-1)*rdndif(i,k-1))
               rdirexp     = rdir(i,k-1)*exptdn(i,k-1)
               tdnmexp     = tottrn(i,k-1) - exptdn(i,k-1)
               tottrn(i,k) = exptdn(i,k-1)*tdir(i,k-1) + tdif(i,k-1)*
     $                      (tdnmexp + rdndif(i,k-1)*rdirexp)*rdenom
               rdndif(i,k) = rdif(i,k-1)  +
     $                (rdndif(i,k-1)*tdif(i,k-1))*(tdif(i,k-1)*rdenom)
C
            end do
         end do
C
C Compute next layer delta-Eddington solution only if total transmission
C of radiation to the interface just above the layer exceeds trmin.
C
         call whenfgt(plon,tottrn(1,k),1,trmin,index,nval)
         if(nval.gt.0) then
CDIR$ IVDEP
            do 100 ii=1,nval
               i=index(ii)
C
C Remember, no ozone absorption in this layer:
C
               tauray(i) = trayoslp*pflx(i,plevp)
               taugab(i) = abh2o*uth2o(i) +
     $                     abco2*utco2(i) + abo2*uto2(i)
C
               tautot    = tauray(i) + taugab(i) + tauaer(i)
C
               wtot      = (wray*tauray(i) + waer(i)*tauaer(i))
     $                             /tautot
C
               gtot      = (gray*wray*tauray(i)  +
     $                      gaer(i)*waer(i)*tauaer(i)) 
     $                             / (wtot*tautot)
C
               ftot      = (fray*wray*tauray(i)  +
     $                      faer(i)*waer(i)*tauaer(i)) 
     $                             / (wtot*tautot)
C
               ts        = taus(wtot,ftot,tautot)
               ws        = omgs(wtot,ftot)
               gs        = asys(gtot,ftot)
               lm        = el(ws,gs)
               alp       = alpha(ws,coszrs(i),gs,lm)
               gam       = gamma(ws,coszrs(i),gs,lm)
               ue        = u(ws,gs,lm)
C
C Limit argument of exponential to 25, in case lm very large:
C
               arg       = amin1(lm*ts,25.)
               extins    = exp(-arg)
               ne        = n(ue,extins)
C
               rdif(i,k) = (ue+1.)*(ue-1.)*(1./extins - extins)/ne
               tdif(i,k) =   4.*ue/ne
C
C Limit argument of exponential to 25, in case coszrs is very small:
C
               arg       = amin1(ts/coszrs(i),25.)
               explay(i,k) = exp(-arg)
C
               apg       = alp + gam
               amg       = alp - gam
               rdir(i,k) = amg*(tdif(i,k)*explay(i,k) - 1.) +
     $                     apg*rdif(i,k)
               tdir(i,k) = apg*tdif(i,k) +
     $                     (amg*rdif(i,k) - (apg-1.))*explay(i,k)
C
C Under rare conditions, reflectivies and transmissivities can be
C negative; zero out any negative values
C
               rdir(i,k) = amax1(rdir(i,k),0.0)
               tdir(i,k) = amax1(tdir(i,k),0.0)
               rdif(i,k) = amax1(rdif(i,k),0.0)
               tdif(i,k) = amax1(tdif(i,k),0.0)
  100       continue
         end if
C
  200 continue
C
C Compute total direct beam transmission, total transmission, and
C reflectivity for diffuse radiation (from below) for both layers
C above the surface:
C
      k = 2
      do nn=1,nloop
         do i=is(nn),ie(nn)
            exptdn(i,k) = exptdn(i,k-1)*explay(i,k-1)
            rdenom      = 1./(1. - rdif(i,k-1)*rdndif(i,k-1))
            rdirexp     = rdir(i,k-1)*exptdn(i,k-1)
            tdnmexp     = tottrn(i,k-1) - exptdn(i,k-1)
            tottrn(i,k) = exptdn(i,k-1)*tdir(i,k-1) + tdif(i,k-1)*
     $                   (tdnmexp + rdndif(i,k-1)*rdirexp)*rdenom
            rdndif(i,k) = rdif(i,k-1)  +
     $                  (rdndif(i,k-1)*tdif(i,k-1))*(tdif(i,k-1)*rdenom)
         end do
      end do
C 
      return
      end
      subroutine radclw(lat     ,ts      ,tnm     ,qnm     ,o3vmr   ,
     $                  pmid    ,pint    ,pmln    ,piln    ,plco2   ,
     $                  plh2o   ,n2o     ,ch4     ,cfc11   ,cfc12   ,
     $                  cld     ,tclrsf  ,qrl     ,flns    ,flnt    ,
     $                  flnsc   ,flntc   ,flwds   ,lwup    )
C-----------------------------------------------------------------------
C
C Compute longwave radiation heating rates and boundary fluxes
C
C Uses broad band absorptivity/emissivity method to compute clear sky;
C assumes randomly overlapped clouds with variable cloud emissivity to
C include effects of clouds.
C
C Computes clear sky absorptivity/emissivity at lower frequency (in
C general) than the model radiation frequency; uses previously computed
C and stored values for efficiency
C
C Note: This subroutine contains vertical indexing which proceeds
C       from bottom to top rather than the top to bottom indexing 
C       used in the rest of the model.
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Kiehl, B. Briegleb, August 1992
C
C-----------------------------------------------------------------------
c
c $Id: radclw.F,v 1.3 1995/03/03 17:44:43 bonan Exp $
c $Author: bonan $
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
C------------------------------Parameters-------------------------------
c
c $Id: prgrid.h,v 1.2 1995/02/10 01:09:10 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation resolution and I/O parameters
C
      integer
     $     plon,    ! number of longitudes
     $     plev,    ! number of vertical levels
     $     plat,    ! number of latitudes
     $     pcnst,   ! number of constituents (including water vapor)
     $     plevmx,  ! number of subsurface levels
     $     plevp,   ! plev + 1
     $     nxpt,    ! no.of points outside active domain for interpolant
     $     jintmx,  ! number of extra latitudes in polar region
     $     plond,   ! slt extended domain longitude
     $     platd,   ! slt extended domain lat.
     $     plevd    ! fold plev,pcnst indices into one
      integer
     $     plngbuf  ! length of absorptivity/emissivity record
C
      parameter(plon    = 1,
     $          plev    = 18,
     $          plat    = 1,
     $          pcnst   = 1,
     $          plevmx  = 4,
     $          plevp   = plev + 1,
     $          nxpt    = 1,
     $          jintmx  = 1,
     $          plond   = plon + 1 + 2*nxpt,
     $          platd   = plat + 2*nxpt + 2*jintmx,
     $          plevd   = plev*(3 + pcnst))
      parameter(plngbuf = 512*((plond*plevp*plevp + plond*plev*4 +
     $                          plond*plevp)/512 + 1))
C
c
c $Id: pagrid.h,v 1.3 1995/03/03 17:47:17 bonan Exp $
c $Author: bonan $
c
C
C Model grid point resolution parameters.
C
      integer
     $     plnlv,    ! Length of multilevel field slice
     $     plndlv,   ! Length of multilevel 3-d field slice
     $     pbflnb,   ! Length of buffer 1
     $     pbflna,   ! Length of buffer 2
     $     pbflnm1,  ! Length of buffer m1
     $     pflenb,   ! Length of buffer 1, padded for unblocked I/O
     $     pflena,   ! Length of buffer 2, padded for unblocked I/O
     $     plenalcl, ! Length of buffer 2, needed in SPEGRD
     $     ptifld,   ! Number of fields on time-invariant boundary dataset
     $     ptvsfld,  ! Number of fields on time-variant boundary dataset
     $     ptvofld,  ! Number of fields on ozone dataset
     $     plenhi,   ! Length of integer header record
     $     plenhc,   ! Length of character header record
     $     plenhr,   ! Length of real header record
     $     plexbuf,  ! Length of communication buffer for flux coupling
     $     ptapes,   ! Maximum number of history tapes allowed
     $     pflds     ! Number of fields in master field list
      integer
     $     ptileni,  ! Length of time-invariant integer header
     $     ptilenc,  ! Length of time-invariant character header
     $     ptvoleni, ! Length of ozone integer header
     $     ptvolenc, ! Length of ozone character header
     $     ptvsleni, ! Length of time-variant integer header
     $     ptvslenc  ! Length of time-variant character header
      integer
     $     plenhis,  ! Length of integer header scalars
     $     plenhcs,  ! Length of character header scalars
     $     ptilenis, ! Length of time-invariant integer scalars
     $     ptilencs, ! Length of time-invariant character scalars
     $     ptolenis, ! Length of ozone integer header scalars
     $     ptolencs, ! Length of ozone character header scalars
     $     ptslenis, ! Length of time-variant integer header scalars
     $     ptslencs  ! Length of time-variant character header scalars
C
      parameter(plnlv=plon*plev,plndlv=plond*plev)
C
C In pbflnb, 9 multi-level fields include the plev levels of plol and
C plos. 2 multi-level fields are pcnst-dependent.
C There are plevmx sub-surface temperature fields. (See User's Guide 
C for complete buffer description)
C There are 4 single-level fields to hold albedos
C
      parameter(pbflnb=(7 + 2*pcnst)*plndlv + (15+plevmx+pcnst)*plond,
C
C In pbflna, there are 3 multi-level and 3 single-level fields.
C
     $          pbflna = (3 + 3*plev)*plond,
     $          pbflnm1 = (1 + 2*plev)*plond,
     $          pflenb = ((pbflnb + pbflnm1)/512 + 1)*512,
     $          pflena = (pbflna/512 + 1)*512,
C
C plenalcl is the buffer size as required in SPEGRD.  
C Only pflena is read/written.
C
     $          plenalcl = ((pbflna + 3*plndlv + plond)/512 + 1)*512,
     $          plexbuf = (((1 + 7*plev)*plond)/512+1)*512,
     $          ptapes = 6,
C
C 8 fields in master list are pcnst-dependent 2 fields occur only
C if pcnst > 1
C
     $          pflds=92+8*pcnst+2*(pcnst-1)+plevmx)
      parameter(ptifld = 11, ptvsfld = 1, ptvofld = 2)
C
C There are 37 scalar words in the integer header and 89 scalar words
C in the character header
C
      parameter(plenhis=37,plenhcs=89,
     $          plenhi=plenhis+3*pflds,plenhc=plenhcs+2*pflds,
     $          plenhr=3*(2*plev + 1) + 2*plat,
     $          ptilenis=plenhis, ptilencs=plenhcs,
     $          ptileni=ptilenis+3*ptifld, ptilenc=ptilencs+2*ptifld,
     $          ptolenis=plenhis, ptolencs=plenhcs,
     $          ptvoleni=ptolenis+3*ptvofld,ptvolenc=ptolencs+2*ptvofld,
     $          ptslenis=plenhis, ptslencs=plenhcs,
     $          ptvsleni=ptslenis+3*ptvsfld,ptvslenc=ptslencs+2*ptvsfld)
C-----------------------------------------------------------------------
      integer plevp2,plevp3,plevp4
      parameter (plevp2=plev+2,plevp3=plev+3,plevp4=plev+4)
C------------------------------Commons----------------------------------
c
c $Id: comlun.h,v 1.1.1.1 1995/02/09 23:26:42 ccm2 Exp $
c $Author: ccm2 $
c
C
C Logical unit numbers and related variables
C
      integer
     $     pnrg1,             ! maximum number of primary 
C                             !  regeneration files
     $     pnrg2,             ! maximum number of secondary 
C                             !  regeneration files
     $     pnrg3              ! maximum number of secondary 

      parameter (pnrg1 = 5)
      parameter (pnrg2 = 5)
      parameter (pnrg3 = 5)
C
      common/comlun/nsds    ,nrg     ,nrg1(pnrg1)      ,nrg2(pnrg2),
     $              nrg3(pnrg3,ptapes)        ,nra1    ,nrb1    ,
     $              ninit   ,nbndti  ,nozone  ,nsst    ,nabem   ,
     $              nsplit,lutag(99)
      common/comlun/rg1lat(pnrg1+1)  ,rg1siz(pnrg1)    ,rg1buf  ,nnrg1,
     $              rg2lat(pnrg2+1)  ,rg2siz(pnrg2)    ,rg2buf  ,nnrg2,
     $              rg3lat(pnrg3+1,ptapes)    ,rg3siz(pnrg3,ptapes)   ,
     $              rg3buf(ptapes)   ,nnrg3(ptapes)    ,mxszrg  ,
     $              nrefrq  ,rgnht(ptapes)    ,rg3num  ,mresfq
      common/comlunc/rg1ext(pnrg1)   ,rg2ext(pnrg2)    ,
     $               rg3ext(pnrg3,ptapes)
C
      integer nsds,     ! restart dataset unit
     $        nrg,      ! master regeneration dataset unit
     $        nrg1,     ! primary regeneration dataset units
     $        nrg2,     ! secondary regeneration dataset units
     $        nrg3,     ! hbuf regeneration dataset units
     $        nra1,     ! a work file
     $        nrb1,     ! b work file
     $        ninit,    ! initial dataset unit
     $        nbndti,   ! time-invariant boundary dataset
     $        nozone,   ! ozone dataset
     $        nsst,     ! sst dataset
     $        nabem,    ! absorptivity/emissivity work file
     $        nsplit    ! communication between LINEMS1 and LINEMS2
C
      logical
     $     lutag        ! list of flags marking logical units in use
      integer
     $     rg1lat,      ! latitude list for primary regen datasets
     $     rg1siz,      ! file sizes for preallocation
     $     rg1buf,      ! buffer length for assign
     $     nnrg1,       ! number of primary regen files written
     $     rg2lat,      ! lat list for secondary regen datasets
     $     rg2siz,      ! file size for preallocation
     $     rg2buf,      ! buffer length for assign
     $     nnrg2,       ! number of secondary regen files written
     $     rg3lat,      ! latitude list for hbuf regen datasets
     $     rg3siz,      ! file sizes for preallocation
     $     rg3buf,      ! buffer length for assign
     $     nnrg3,       ! number of hbuf regen files written
     $     mxszrg,      ! max size of a regen file (megabytes)
     $     nrefrq,      ! frequency of regeneration file writes
     $     mresfq,      ! frequency of mnthly avg regen file writes
     $     rg3num       ! number of temporary secondary regen files
      logical
     $     rgnht        ! set true if regeneration file for a h-tape exists
      character*2
     $     rg1ext,      ! file extension for primary regen files
     $     rg2ext       ! file extension for secondary regen files
      character*5
     $     rg3ext       ! file extension for secondary regen files
C
C-----------------------------------------------------------------------
c
c $Id: comtim.h,v 1.1.1.1 1995/02/09 23:26:44 ccm2 Exp $
c $Author: ccm2 $
c
C
C Model time variables
C
      common/comtim/calday  ,dtime   ,twodt   ,nrstrt  ,nstep   ,
     $              nstepr  ,nestep  ,nelapse ,nstop   ,mdbase  ,
     $              msbase  ,mdcur   ,mscur   ,mbdate  ,mbsec   ,
     $              mcdate  ,mcsec   ,nndbas  ,nnsbas  ,nnbdat  ,
     $              nnbsec  ,doabsems,dosw    ,dolw
C
      real calday,   ! Current calendar day = julian day + fraction
     $     dtime,    ! Time step in seconds (delta t)
     $     twodt     ! 2 * delta t 
      integer
     $     nrstrt,   ! Starting time step of restart run (constant) 
     $     nstep,    ! Current time step
     $     nstepr,   ! Current time step of restart run(updated w/nstep)
     $     nestep,   ! Time step on which to stop run
     $     nelapse,  ! Requested elapsed time for model run
     $     nstop,    ! nestep + 1
     $     mdbase,   ! Base day of run
     $     msbase,   ! Base seconds of base day
     $     mdcur,    ! Current day of run
     $     mscur,    ! Current seconds of current day
     $     mbdate,   ! Base date of run (yymmdd format)
     $     mbsec,    ! Base seconds of base date
     $     mcdate,   ! Current date of run (yymmdd format)
     $     mcsec,    ! Current seconds of current date
     $     nndbas,   ! User input base day
     $     nnsbas,   ! User input base seconds of input base day
     $     nnbdat,   ! User input base date (yymmdd format)
     $     nnbsec    ! User input base seconds of input base date
      logical
     $     doabsems, ! True => abs/emiss calculation this timestep
     $     dosw,     ! True => shortwave calculation this timestep
     $     dolw      ! True => longwave calculation this timestep
C
C-----------------------------------------------------------------------
c
c $Id: crdcon.h,v 1.1.1.1 1995/02/09 23:26:46 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation constants
C
      common/crdcon/gravit  ,rga     ,cpair   ,epsilo  ,sslp    ,
     $              stebol  ,rgsslp  ,co2vmr  ,dpfo3   ,dpfco2  ,
     $              dayspy  ,pie
C
      real gravit,    ! Acceleration of gravity
     $     rga,       ! 1./gravit
     $     cpair,     ! Specific heat of dry air
     $     epsilo,    ! Ratio of mol. wght of H2O to dry air
     $     sslp,      ! Standard sea-level pressure
     $     stebol,    ! Stefan-Boltzmann's constant
     $     rgsslp,    ! 0.5/(gravit*sslp)
     $     co2vmr,    ! CO2 volume mixing ratio
     $     dpfo3,     ! Voigt correction factor for O3
     $     dpfco2,    ! Voigt correction factor for CO2
     $     dayspy,    ! Number of days per 1 year
     $     pie        ! 3.14.....
C
C-----------------------------------------------------------------------
c
c $Id: comctl.h,v 1.1.1.1 1995/02/09 23:26:41 ccm2 Exp $
c $Author: ccm2 $
c
C
C Model control variables
C
      common/comctl/itsst   ,nsrest  ,iradsw  ,iradlw  ,iradae  ,
     $              anncyc  ,nlend   ,nlres   ,nlhst   ,lbrnch  ,
     $              ldebug  ,aeres   ,ozncyc  ,sstcyc  ,dodiavg ,
     $              aeregen ,cpuchek
      integer
     $     itsst,   ! Sea surf. temp. update freq. (iters)
     $     nsrest,  ! Restart flag
     $     iradsw,  ! Iteration frequency for shortwave radiation computation
     $     iradlw,  ! Iteration frequency for longwave radiation computation
     $     iradae   ! Iteration freq. for absorptivity/emissivity comp
      logical
     $     anncyc,  ! Do annual cycle (otherwise perpetual)
     $     nlend,   ! Flag for end of run
     $     nlres,   ! If true, continuation run
     $     nlhst,   ! If true, regeneration run
     $     lbrnch,  ! If true, branch run
     $     ldebug,  ! If in debug mode, link output files to /usr/tmp
C                   !    before mswrite, and remove all but last file
     $     aeres,   ! If true, a/e data will be stored on restart file
     $     ozncyc,  ! If true, cycle ozone dataset
     $     sstcyc,  ! If true, cycle sst dataset
     $     dodiavg, ! true => diurnal averaging
     $     aeregen, ! true => absor/emis part of regeneration data
     $     cpuchek  ! If true, check remaining cpu time at each writeup
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer lat                  ! Model latitude index
      real ts(plond)               ! Ground (skin) temperature
C
C Input arguments which are only passed to other routines
C
      real tnm(plond,plev),        ! Level temperature
     $     qnm(plond,plev),        ! Level moisture field
     $     o3vmr(plond,plev),      ! ozone volume mixing ratio
     $     pmid(plond,plev),       ! Level pressure
     $     pint(plond,plevp),      ! Model interface pressure
     $     pmln(plond,plev),       ! Ln(pmid)
     $     piln(plond,plevp),      ! Ln(pint)
     $     plco2(plond,plevp),     ! Path length co2
     $     plh2o(plond,plevp)      ! Path length h2o
      real n2o(plond,plev),        ! nitrous oxide mass mixing ratio
     $     ch4(plond,plev),        ! methane mass mixing ratio
     $     cfc11(plond,plev),      ! cfc11 mass mixing ratio
     $     cfc12(plond,plev)       ! cfc12 mass mixing ratio
C
C Input/Output arguments
C
      real cld(plond,plevp),       ! Cloud cover
     $     tclrsf(plond,plevp)     ! Clear sky fraction
C
C Output arguments
C
      real qrl(plond,plev),        ! Longwave heating rate
     $     flns(plond),            ! Surface cooling flux
     $     flnt(plond),            ! Net outgoing flux
     $     flnsc(plond),           ! Clear sky surface cooing
     $     flntc(plond),           ! Net clear sky outgoing flux
     $     flwds(plond),           ! Down longwave flux at surface
     $     lwup(plond)             ! Surface longwave up flux from coupler
C
C---------------------------Local variables-----------------------------
C
      integer     i,    ! Longitude index
     $            k,    ! Level index
     $           k1,    ! Level index
     $           k2,    ! Level index
     $           k3,    ! Level index
     $           km,    ! Level index
     $          km1,    ! Level index
     $          km2,    ! Level index
     $          km3,    ! Level index
     $          km4     ! Level index
C
      real tmp(plond),
     $     tmp1,                ! Temporary 1
     $     absbt(plond)         ! Downward emission at model top
      real plol(plond,plevp),   ! O3 pressure wghted path length
     $     plos(plond,plevp)    ! O3 path length
C
      real co2em(plond,plevp),  ! Layer co2 normalized planck funct. derivative
     $     co2eml(plond,plev),  ! Interface co2 normalized planck funct. deriv.
     $     delt(plond),         ! Diff t**4 mid layer to top interface
     $     delt1(plond),        ! Diff t**4 lower intrfc to mid layer
     $     bk1(plond),          ! Absrptvty for vertical quadrature
     $     bk2(plond),          ! Absrptvty for vertical quadrature
     $     ful(plond,plevp),    ! Total upwards longwave flux
     $     fsul(plond,plevp),   ! Clear sky upwards longwave flux
     $     fdl(plond,plevp),    ! Total downwards longwave flux
     $     fsdl(plond,plevp),   ! Clear sky downwards longwv flux
     $     fclb4(plond,plev),   ! Sig t**4 for cld bottom interfc
     $     fclt4(plond,plev),   ! Sig t**4 for cloud top interfc
     $     s(plond,plevp,plevp) ! Flx integral sum
      real tplnka(plond,plevp), ! Planck fnctn temperature
     $     s2c(plond,plevp),    ! H2o cont amount
     $     s2t(plond,plevp),    ! H2o cont temperature
     $     w(plond,plevp),      ! H2o path
     $     tplnke(plond)        ! Planck fnctn temperature
      real h2otr(plond,plevp),  ! H2o trnmsn for o3 overlap
     $     co2t(plond,plevp),   ! Prs wghted temperature path
     $     tint(plond,plevp),   ! Interface temperature
     $     tint4(plond,plevp),  ! Interface temperature**4
     $     tlayr(plond,plevp),  ! Level temperature
     $     tlayr4(plond,plevp)  ! Level temperature**4
      real rtclrsf(plond,plevp),! 1./tclrsf(i,k)
     $     gocp                 ! gravit/cpair
      integer klov(plond),      ! Cloud lowest level index
     $        khiv(plond),      ! Cloud highest level index
     $        khivm(plond)      ! khiv(i) - 1
      integer indx(plond),
     $        npts,
     $        ii,
     $        khighest
      logical done(plond),
     $        start(plond)
      real absems(plngbuf)      ! Absorbs's and emiss's in buffer
      pointer (pabsnxt,absnxt(plond,plev,4)),
     $        (pabstot,abstot(plond,plevp,plevp)),
     $        (pemstot,emstot(plond,plevp))
      real absnxt,               ! Nearest layer absorptivities
     $     abstot,               ! Non-adjacent layer absorptivites
     $     emstot                ! Total emissivity
c
c Trace gas variables
c
      real ucfc11(plond,plevp), ! CFC11 path length
     $     ucfc12(plond,plevp), ! CFC12 path length
     $     un2o0(plond,plevp),  ! N2O path length
     $     un2o1(plond,plevp),  ! N2O path length (hot band)
     $     uch4(plond,plevp),   ! CH4 path length
     $     uco211(plond,plevp), ! CO2 9.4 micron band path length
     $     uco212(plond,plevp), ! CO2 9.4 micron band path length
     $     uco213(plond,plevp), ! CO2 9.4 micron band path length
     $     uco221(plond,plevp), ! CO2 10.4 micron band path length
     $     uco222(plond,plevp), ! CO2 10.4 micron band path length
     $     uco223(plond,plevp), ! CO2 10.4 micron band path length
     $     bn2o0(plond,plevp),  ! pressure factor for n2o
     $     bn2o1(plond,plevp),  ! pressure factor for n2o
     $     bch4(plond,plevp),   ! pressure factor for ch4
     $     uptype(plond,plevp)  ! p-type continuum path length
      real abplnk1(14,plond,plevp), ! non-nearest layer Plack factor
     $     abplnk2(14,plond,plevp)  ! nearest layer factor
C
C
C------------------------------Externals--------------------------------
C
      integer intmax
      external intmax
      external radtpl,           ! Compute path lengths
     $         radems,           ! H2o,co2,o3 emissivity
     $         radabs,           ! H2o,co2,o3 absorptivity
     $         whenne,
     $         whenflt
      external writeric,         ! Write for abs/ems
     $         readric           ! Read  for abs/ems
C-----------------------------------------------------------------------
C
C Set pointer variables
C
      pabstot = loc(absems(1                                   ))
      pabsnxt = loc(absems(1 + plond*plevp*plevp               ))
      pemstot = loc(absems(1 + plond*plevp*plevp + plond*plev*4))
C
C Initialize and recompute the tclrsf array
C
      do i=1,plon
         rtclrsf(i,1) = 1.0/tclrsf(i,1)
      end do 
C
      do k=1,plev
         do i=1,plon
            fclb4(i,k) = 0.
            fclt4(i,k) = 0.
            tclrsf(i,k+1) = tclrsf(i,k)*(1. - cld(i,k+1))
            rtclrsf(i,k+1) = 1./tclrsf(i,k+1)
         end do
      end do
C
C Calculate some temperatures needed to derive absorptivity and
C emissivity, as well as some h2o path lengths
C
      call radtpl(tnm     ,ts      ,qnm     ,pint    ,plh2o   ,
     $            tplnka  ,s2c     ,s2t     ,w       ,tplnke  ,
     $            tint    ,tint4   ,tlayr   ,tlayr4  ,pmln    ,
     $            piln    )
      if (doabsems) then
C
C Compute ozone path lengths at frequency of a/e calculation.
C
         call radoz2(o3vmr   ,pint    ,plol    ,plos    )
c
c Compute trace gas path lengths
c
         call trcpth(tnm, pint, cfc11, cfc12, n2o, ch4, qnm,
     $               ucfc11, ucfc12, un2o0,  un2o1,  uch4,
     $               uco211, uco212, uco213, uco221, uco222,
     $               uco223, bn2o0,  bn2o1,  bch4,   uptype)
C
C
C Compute total emissivity:
C
         call radems(s2c     ,s2t     ,w       ,tplnke  ,plh2o   ,
     $               pint    ,plco2   ,tint    ,tint4   ,tlayr   ,
     $               tlayr4  ,plol    ,plos    ,ucfc11  ,ucfc12  , 
     $               un2o0   ,un2o1   ,uch4    ,uco211  ,uco212  ,
     $               uco213  ,uco221  ,uco222  ,uco223  ,uptype  ,
     $               bn2o0   ,bn2o1   ,bch4    ,co2em   ,co2eml  ,
     $               co2t    ,h2otr   ,abplnk1 ,abplnk2 ,emstot  )
C
C Compute total absorptivity:
C
         call radabs(pmid    ,pint    ,co2em   ,co2eml  ,tplnka  ,
     $               s2c     ,s2t     ,w       ,h2otr   ,plco2   ,
     $               plh2o   ,co2t    ,tint    ,tlayr   ,plol    ,
     $               plos    ,pmln    ,piln    ,ucfc11  ,ucfc12  , 
     $               un2o0   ,un2o1   ,uch4    ,uco211  ,uco212  ,
     $               uco213  ,uco221  ,uco222  ,uco223  ,uptype  ,
     $               bn2o0   ,bn2o1   ,bch4    ,abplnk1 ,abplnk2 ,
     $               abstot  ,absnxt  )
C
C Write abs/ems info to ssd.  Note: Need not be done with INCORERAD mods
C
         call writeric(nabem   ,absems(1),plngbuf ,lat     )
      else
C
C Get total abs/ems info from ssd.  Note: Need not be done with INCORERAD mods
C
         call readric(nabem   ,absems(1),plngbuf ,lat     )
      end if
C
C Find the lowest and highest level cloud for each grid point
C Note: Vertical indexing here proceeds from bottom to top
C
      do i=1,plon
         klov(i) = 0
         done(i) = .false.
      end do
      do k=1,plev
         do i=1,plon
            if (.not.done(i) .and. cld(i,plevp2-k).gt.0.0) then
               done(i) = .true.
               klov(i) = k
            end if
         end do
      end do
      call whenne(plon,klov,1,0,indx,npts)
      do i=1,plon
         khiv(i) = klov(i)
         done(i) = .false.
      end do
      do k=plev,1,-1
CDIR$ IVDEP
         do ii=1,npts
            i=indx(ii)
            if (.not.done(i) .and. cld(i,plevp2-k).gt.0.0) then
               done(i) = .true.
               khiv(i) = k
            end if
         end do
      end do
      do i=1,plon
         khivm(i) = khiv(i) - 1
      end do
C
C Note: Vertical indexing here proceeds from bottom to top
C
      do ii=1,npts
         i=indx(ii)
         do k=klov(i),khiv(i)
            fclt4(i,plevp-k) = stebol*tint4(i,plevp2-k)
            fclb4(i,plevp-k) = stebol*tint4(i,plevp3-k)
         end do
      end do
C
C Compute sums used in integrals (all longitude points)
C
C Definition of bk1 & bk2 depends on finite differencing.  for
C trapezoidal rule bk1=bk2. trapezoidal rule applied for nonadjacent
C layers only.
C
C delt=t**4 in layer above current sigma level km.
C delt1=t**4 in layer below current sigma level km.
C
      do i=1,plon
         delt(i) = tint4(i,plev) - tlayr4(i,plevp)
         delt1(i) = tlayr4(i,plevp) - tint4(i,plevp)
         s(i,plevp,plevp) = stebol*(delt1(i)*absnxt(i,plev,1) +
     $                              delt (i)*absnxt(i,plev,4))
         s(i,plev,plevp)  = stebol*(delt (i)*absnxt(i,plev,2) +
     $                              delt1(i)*absnxt(i,plev,3))
      end do
      do k=1,plev-1
         do i=1,plon
            bk2(i) = (abstot(i,k,plev) + abstot(i,k,plevp))*0.5
            bk1(i) = bk2(i)
            s(i,k,plevp) = stebol*(bk2(i)*delt(i) + bk1(i)*delt1(i))
         end do
      end do
C
C All k, km>1
C
      do km=plev,2,-1
         do i=1,plon
            delt(i)  = tint4(i,km-1) - tlayr4(i,km)
            delt1(i) = tlayr4(i,km) - tint4(i,km)
         end do
         do k=plevp,1,-1
            if (k.eq.km) then
               do i=1,plon
                  bk2(i) = absnxt(i,km-1,4)
                  bk1(i) = absnxt(i,km-1,1)
               end do
            else if(k.eq.km-1) then
               do i=1,plon
                  bk2(i) = absnxt(i,km-1,2)
                  bk1(i) = absnxt(i,km-1,3)
               end do
            else
               do i=1,plon
                  bk2(i) = (abstot(i,k,km-1) + abstot(i,k,km))*0.5
                  bk1(i) = bk2(i)
               end do
            end if
            do i=1,plon
               s(i,k,km) = s(i,k,km+1) + stebol*
     $                    (bk2(i)*delt(i) + bk1(i)*delt1(i))
            end do
         end do
      end do
C
C Computation of clear sky fluxes always set first level of fsul
C
      do i=1,plon
         fsul(i,plevp) = stebol*(ts(i)**4)
      end do
C
C Downward clear sky fluxes store intermediate quantities in down flux
C Initialize fluxes to clear sky values.
C
      do i=1,plon
         tmp(i) = fsul(i,plevp) - stebol*tint4(i,plevp)
         fsul(i,1) = fsul(i,plevp) - abstot(i,1,plevp)*tmp(i) + s(i,1,2)
         fsdl(i,1) = stebol*(tplnke(i)**4)*emstot(i,1)
         ful(i,1) = fsul(i,1)
         fdl(i,1) = fsdl(i,1)
      end do
C
C fsdl(i,plevp) assumes isothermal layer
C
      do k=2,plev
         do i=1,plon
            fsul(i,k) = fsul(i,plevp) - abstot(i,k,plevp)*tmp(i) +
     $                  s(i,k,k+1)
            ful(i,k) = fsul(i,k)
            fsdl(i,k) = stebol*(tplnke(i)**4)*emstot(i,k) -
     $                  (s(i,k,2) - s(i,k,k+1))
            fdl(i,k) = fsdl(i,k)
         end do
      end do
C
C Store the downward emission from level 1 = total gas emission * sigma
C t**4.  fsdl does not yet include all terms
C
      do i=1,plon
         ful(i,plevp) = fsul(i,plevp)
         absbt(i) = stebol*(tplnke(i)**4)*emstot(i,plevp)
         fsdl(i,plevp) = absbt(i) - s(i,plevp,2)
         fdl(i,plevp) = fsdl(i,plevp)
      end do
C
C Modifications for clouds
C
C Further qualify longitude subset for computations.  Select only those
C locations where there are clouds (total cloud fraction <= 1.e-3 treated 
C as clear)
C
      call whenflt(plon,tclrsf(1,plevp),1,0.999,indx,npts)
C
C Compute downflux at level 1 for cloudy sky
C
      do ii=1,npts
         i=indx(ii)
C
C First clear sky flux plus flux from cloud at level 1
C
         fdl(i,plevp) = fsdl(i,plevp)*tclrsf(i,plev)*
     $         rtclrsf(i,plevp-khiv(i)) + fclb4(i,plev-1)*cld(i,plev)
      end do
C
C Flux emitted by other layers
C Note: Vertical indexing here proceeds from bottom to top
C
      khighest = khiv(intmax(plon,khiv,1))
      do km=3,khighest
         km1 = plevp - km
         km2 = plevp2 - km
         km4 = plevp4 - km
CDIR$ IVDEP
         do ii=1,npts
            i=indx(ii)
            if (km.le.khiv(i)) then
               tmp1 = cld(i,km2)*tclrsf(i,plev)*rtclrsf(i,km2)
               fdl(i,plevp) = fdl(i,plevp) +
     $                       (fclb4(i,km1) - s(i,plevp,km4))*tmp1
            end if
         end do
      end do
C
C Note: Vertical indexing here proceeds from bottom to top
C
      do k=1,khighest-1
         k1 = plevp - k
         k2 = plevp2 - k
         k3 = plevp3 - k
CDIR$ IVDEP
         do ii=1,npts
            i=indx(ii)
            if (k.ge.klov(i) .and. k.le.khivm(i)) then
               ful(i,k2) = fsul(i,k2)*(tclrsf(i,plevp)*rtclrsf(i,k1))
            end if
         end do
         do km=1,k
            km1 = plevp - km
            km2 = plevp2 - km
            km3 = plevp3 - km
CDIR$ IVDEP
            do ii=1,npts
               i=indx(ii)
               if (k.le.khivm(i) .and. km.ge.klov(i) .and.
     $             km.le.khivm(i)) then
C
                  ful(i,k2) = ful(i,k2) +
     $                 (fclt4(i,km1) + s(i,k2,k3) - s(i,k2,km3))*
     $                 cld(i,km2)*(tclrsf(i,km1)*rtclrsf(i,k1))
               end if
            end do
         end do            ! km=1,k
      end do               ! k=1,khighest-1
C
      do k=1,plevp
         k2 = plevp2 - k
         k3 = plevp3 - k
         do i=1,plon
            start(i) = .false.
         end do
CDIR$ IVDEP
         do ii=1,npts
            i=indx(ii)
            if (k.ge.khiv(i)) then
               start(i) = .true.
               ful(i,k2) = fsul(i,k2)*tclrsf(i,plevp)*
     $                     rtclrsf(i,plevp-khiv(i))
            end if
         end do
         do km=1,khighest
            km1 = plevp - km
            km2 = plevp2 - km
            km3 = plevp3 - km
CDIR$ IVDEP
            do ii=1,npts
               i=indx(ii)
               if(start(i) .and. km.ge.klov(i) .and. km.le.khiv(i)) then
                  ful(i,k2) = ful(i,k2)  +
     $              (cld(i,km2)*tclrsf(i,km1)*rtclrsf(i,plevp-khiv(i)))*
     $              (fclt4(i,km1) + s(i,k2,k3) - s(i,k2,km3))
               end if
            end do
         end do         ! km=1,khighest
      end do            ! k=1,plevp
C
C Computation of the downward fluxes
C
      do k=2,khighest-1
         k1 = plevp - k
         k2 = plevp2 - k
         k3 = plevp3 - k
CDIR$ IVDEP
         do ii=1,npts
            i=indx(ii)
            if (k.le.khivm(i)) fdl(i,k2) = 0.
         end do
         do km=k+1,khighest
            km1 = plevp - km
            km2 = plevp2 - km
            km4 = plevp4 - km
CDIR$ IVDEP
            do ii=1,npts
               i=indx(ii)
               if (k.le.khiv(i) .and. km.ge.max0(k+1,klov(i)) .and.
     $             km.le.khiv(i)) then
C
                  fdl(i,k2) = fdl(i,k2) +
     $                 (cld(i,km2)*tclrsf(i,k1)*rtclrsf(i,km2))*
     $                 (fclb4(i,km1) - s(i,k2,km4) + s(i,k2,k3))
               end if
            end do
         end do            ! km=k+1,khighest
CDIR$ IVDEP
         do ii=1,npts
            i=indx(ii)
            if (k.le.khivm(i)) then
               fdl(i,k2) = fdl(i,k2) + fsdl(i,k2)*
     $              (tclrsf(i,k1)*rtclrsf(i,plevp-khiv(i)))
            end if
         end do
      end do               ! k=1,khighest-1
C
C End cloud modification loops
C
C All longitudes: store history tape quantities
C
      do i=1,plon
C
C Downward longwave flux
C
         flwds(i) = fdl(i,plevp)
C
C Net flux
C
         flns(i) = ful(i,plevp) - fdl(i,plevp)
C
C Clear sky flux at top of atmosphere
C
         flntc(i) = fsul(i,1)
         flnsc(i) = fsul(i,plevp) - fsdl(i,plevp)
C
C Outgoing ir
C
         flnt(i) = ful(i,1) - fdl(i,1)
      end do
C
C Computation of longwave heating (k per sec)
C
      gocp = gravit/cpair
      do k=1,plev
         do i=1,plon
            qrl(i,k) = (ful(i,k) - fdl(i,k) - ful(i,k+1) + fdl(i,k+1))*
     $                 gocp/((pint(i,k) - pint(i,k+1)))
         end do
      end do
C
      return
      end
      subroutine radcsw(pint    ,h2ommr  ,cld     ,clwp    ,o3mmr   ,
     $                  eccf    ,coszrs  ,albs    ,albsd   ,albl    ,
     $                  albld   ,solin   ,qrs     ,fsns    ,fsnt    ,
     $                  fsnsc   ,fsntc   ,rel     ,rei     ,fice    ,
     $                  sols    ,soll    ,solsd   ,solld   ,aermmr  ,
     $                  rh      )
C-----------------------------------------------------------------------
C
C Solar radiation code
C
C Computes incident solar flux, solar heating rate, surface absorbed
C solar flux, and total column absorbed solar flux
C
C Uses the delta-eddington method
C
C Divides solar spectrum into 18 intervals from 0.2-5.0 micro-meters.
C solar flux fractions specified for each interval. allows for
C seasonally and diurnally varying solar input.  Includes molecular,
C cloud, and surface scattering, along with h2o,o3,co2,o2,cloud, and
C surface absorption. Computes delta-eddington reflections and
C transmissions assuming homogeneously mixed layers. Adds the layers 
C assuming scattering between layers to be isotropic, and distinguishes 
C direct solar beam from scattered radiation.
C
C Longitude loops are broken into 1 or 2 sections, so that only daylight
C (i.e. coszrs > 0) computations are done.
C
C Note that an extra layer above the model top layer is added.
C
C cgs units are used.
C
C Special diagnostic calc of the clear sky surface and total column
C absorbed flux is also done; this calculation does not effect the rest
C of the model, but is included for cloud forcing diagnostics.
C
C For more details , see Briegleb, Bruce P., 1992: Delta-Eddington
C Approximation for Solar Radiation in the NCAR Community Climate Model,
C Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
C
C---------------------------Code history--------------------------------
C
C Modified March 1995 to add aerosols
C Original version:  B. Briegleb
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Kiehl, B. Briegleb, August 1992
C
C-----------------------------------------------------------------------
c
c $Id: radcsw.F,v 1.4 1995/03/17 18:54:09 ccm2 Exp $
c $Author: ccm2 $
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
C------------------------------Parameters-------------------------------
c
c $Id: prgrid.h,v 1.2 1995/02/10 01:09:10 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation resolution and I/O parameters
C
      integer
     $     plon,    ! number of longitudes
     $     plev,    ! number of vertical levels
     $     plat,    ! number of latitudes
     $     pcnst,   ! number of constituents (including water vapor)
     $     plevmx,  ! number of subsurface levels
     $     plevp,   ! plev + 1
     $     nxpt,    ! no.of points outside active domain for interpolant
     $     jintmx,  ! number of extra latitudes in polar region
     $     plond,   ! slt extended domain longitude
     $     platd,   ! slt extended domain lat.
     $     plevd    ! fold plev,pcnst indices into one
      integer
     $     plngbuf  ! length of absorptivity/emissivity record
C
      parameter(plon    = 1,
     $          plev    = 18,
     $          plat    = 1,
     $          pcnst   = 1,
     $          plevmx  = 4,
     $          plevp   = plev + 1,
     $          nxpt    = 1,
     $          jintmx  = 1,
     $          plond   = plon + 1 + 2*nxpt,
     $          platd   = plat + 2*nxpt + 2*jintmx,
     $          plevd   = plev*(3 + pcnst))
      parameter(plngbuf = 512*((plond*plevp*plevp + plond*plev*4 +
     $                          plond*plevp)/512 + 1))
C
C-----------------------------------------------------------------------
      real scon                ! Solar constant
      parameter (scon = 1.367e6)
C------------------------------Commons----------------------------------
c
c $Id: crdcon.h,v 1.1.1.1 1995/02/09 23:26:46 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation constants
C
      common/crdcon/gravit  ,rga     ,cpair   ,epsilo  ,sslp    ,
     $              stebol  ,rgsslp  ,co2vmr  ,dpfo3   ,dpfco2  ,
     $              dayspy  ,pie
C
      real gravit,    ! Acceleration of gravity
     $     rga,       ! 1./gravit
     $     cpair,     ! Specific heat of dry air
     $     epsilo,    ! Ratio of mol. wght of H2O to dry air
     $     sslp,      ! Standard sea-level pressure
     $     stebol,    ! Stefan-Boltzmann's constant
     $     rgsslp,    ! 0.5/(gravit*sslp)
     $     co2vmr,    ! CO2 volume mixing ratio
     $     dpfo3,     ! Voigt correction factor for O3
     $     dpfco2,    ! Voigt correction factor for CO2
     $     dayspy,    ! Number of days per 1 year
     $     pie        ! 3.14.....
C
      real rel(plond,plev),    ! liquid effective drop size (microns)
     $     rei(plond,plev),    ! ice effective drop size (microns)
     $     fice(plond,plev)    ! fractional ice content within cloud
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real pint(plond,plevp),  ! Interface pressure
     $     h2ommr(plond,plev), ! Specific humidity (h2o mass mix ratio)
     $     cld(plond,plevp),   ! Fractional cloud cover
     $     clwp(plond,plev),   ! Layer liquid water path
     $     o3mmr(plond,plev),  ! Ozone mass mixing ratio
     $     eccf,               ! Eccentricity factor(1./earth-sun dist. squared)
     $     coszrs(plond)       ! Cosine solar zenith angle
      real albs(plond),        ! 0.2-0.7 micro-meter srfc alb to direct rad
     $     albl(plond),        ! 0.7-5.0 micro-meter srfc alb to direct rad
     $     albsd(plond),       ! 0.2-0.7 micro-meter srfc alb to diffuse rad
     $     albld(plond)        ! 0.7-5.0 micro-meter srfc alb to diffuse rad
      real aermmr(plond,plev), ! level aerosol mass mixing ratio
     $     rh(plond,plev)      ! level relative humidity (fraction)
C
C
C Output arguments
C
      real solin(plond),       ! Incident solar flux
     $     qrs(plond,plev),    ! Solar heating rate
     $     fsns(plond),        ! Surface absorbed solar flux
     $     fsnt(plond),        ! Total column absorbed solar flux
     $     fsnsc(plond),       ! Clear sky surface absorbed solar flux
     $     fsntc(plond),       ! Clear sky total column absorbed solar flx
     $     sols(plond),        ! direct solar rad incident on surface (< 0.7)
     $     soll(plond),        ! direct solar rad incident on surface (>= 0.7)
     $     solsd(plond),       ! diffuse solar rad incident on surface (< 0.7)
     $     solld(plond)        ! diffuse solar rad incident on surface (>= 0.7)

C
C------------------------------Externals--------------------------------
C
      integer   isrchfgt,      ! Search for first array element > 0
     $          isrchfle       ! Search for first array element < 0
      external    radded,      ! Computes delta-eddington solution
     $            radclr,      ! Computes clear sky delta-edd solution
     $          isrchfgt,      ! Search for first array element > 0
     $          isrchfle       ! Search for first array element < 0
C
C---------------------------Local variables-----------------------------
C
      integer       ns,        ! Spectral loop index
     $               i,        ! Longitude loop index
     $               k,        ! Level loop index
     $               n,        ! Loop index for daylight
     $           nloop,        ! Number of daylight loops
     $           is(2),        ! Daytime start indices
     $           ie(2),        ! Daytime end indices
     $          indxsl         ! Index for cloud particle properties
C
C A. Slingo's data for cloud particle radiative properties (from 'A GCM
C Parameterization for the Shortwave Properties of Water Clouds' JAS
C vol. 46 may 1989 pp 1419-1427)
C
      real abarl(4), ! A coefficient for extinction optical depth
     $     bbarl(4), ! B coefficiant for extinction optical depth
     $     cbarl(4), ! C coefficiant for single particle scat albedo
     $     dbarl(4), ! D coefficiant for single particle scat albedo
     $     ebarl(4), ! E coefficiant for asymmetry parameter
     $     fbarl(4)  ! F coefficiant for asymmetry parameter
      save abarl, bbarl, cbarl, dbarl, ebarl, fbarl
C
      data abarl/ 2.817e-02, 2.682e-02,2.264e-02,1.281e-02/
      data bbarl/ 1.305    , 1.346    ,1.454    ,1.641    /
      data cbarl/-5.62e-08 ,-6.94e-06 ,4.64e-04 ,0.201    /
      data dbarl/ 1.63e-07 , 2.35e-05 ,1.24e-03 ,7.56e-03 /
      data ebarl/ 0.829    , 0.794    ,0.754    ,0.826    /
      data fbarl/ 2.482e-03, 4.226e-03,6.560e-03,4.353e-03/
C
      real abarli,   ! A coefficiant for current spectral interval
     $     bbarli,   ! B coefficiant for current spectral interval
     $     cbarli,   ! C coefficiant for current spectral interval
     $     dbarli,   ! D coefficiant for current spectral interval
     $     ebarli,   ! E coefficiant for current spectral interval
     $     fbarli    ! F coefficiant for current spectral interval
C
C Caution... A. Slingo recommends no less than 4.0 micro-meters nor
C greater than 20 micro-meters
C
C
c    ice water coefficients (Ebert and Curry,1992, JGR, 97, 3831-3836)
c
      real abari(4), ! a coefficient for extinction optical depth
     $     bbari(4), ! b coefficiant for extinction optical depth
     $     cbari(4), ! c coefficiant for single particle scat albedo
     $     dbari(4), ! d coefficiant for single particle scat albedo
     $     ebari(4), ! e coefficiant for asymmetry parameter
     $     fbari(4)  ! f coefficiant for asymmetry parameter
      save abari, bbari, cbari, dbari, ebari, fbari
c
      data abari/ 3.448e-03, 3.448e-03,3.448e-03,3.448e-03/
      data bbari/ 2.431    , 2.431    ,2.431    ,2.431    /
      data cbari/ 1.00e-05 , 1.10e-04 ,1.861e-02,.46658   /
      data dbari/ 0.0      , 1.405e-05,8.328e-04,2.05e-05 /
      data ebari/ 0.7661   , 0.7730   ,0.794    ,0.9595   /
      data fbari/ 5.851e-04, 5.665e-04,7.267e-04,1.076e-04/
C
      real abarii,   ! A coefficiant for current spectral interval
     $     bbarii,   ! B coefficiant for current spectral interval
     $     cbarii,   ! C coefficiant for current spectral interval
     $     dbarii,   ! D coefficiant for current spectral interval
     $     ebarii,   ! E coefficiant for current spectral interval
     $     fbarii    ! F coefficiant for current spectral interval
C
C      real cldefr ! Universal cloud effective radius in micro-meters
C      data cldefr / 10.0 /
C
      real delta           ! Pressure (atmospheres) for stratospheric h2o limit
      data delta  /  1.70e-3 /
C
      real o2mmr           ! O2 mass mixing ratio:
      save delta, o2mmr
      data o2mmr / .23143 /
C
C CO2 info:
C
      real mmwair,         ! Mean molecular weight of air
     $     mmwco2,         ! Mean molecular weight of co2
     $     co2mmr          ! Co2 mass mixing ratio
      save mmwair, mmwco2
      data mmwair / 28.9644 /
      data mmwco2 / 44.0000 /
C
      real albdir(plond),  ! Current spc intrvl srf alb to direct rad
     $     albdif(plond)   ! Current spc intrvl srf alb to diffuse rad
C
      integer nspint  ! Num of spectral intervals across solar spectrum
      parameter ( nspint = 18 )
C
C Next series depends on spectral interval
C
      real frcsol(nspint),  ! Fraction of solar flux in each spectral interval
     $     wavmin(nspint),  ! Min wavelength (micro-meters) of interval
     $     wavmax(nspint),  ! Max wavelength (micro-meters) of interval
     $     raytau(nspint),  ! Rayleigh scattering optical depth
     $     abh2o(nspint),   ! Absorption coefficiant for h2o (cm2/g)
     $     abo3 (nspint),   ! Absorption coefficiant for o3  (cm2/g)
     $     abco2(nspint),   ! Absorption coefficiant for co2 (cm2/g)
     $     abo2 (nspint),   ! Absorption coefficiant for o2  (cm2/g)
     $     ph2o(nspint),    ! Weight of h2o in spectral interval
     $     pco2(nspint),    ! Weight of co2 in spectral interval
     $     po2 (nspint)     ! Weight of o2  in spectral interval
      save frcsol ,wavmin ,wavmax ,raytau ,abh2o ,abo3 ,
     $     abco2  ,abo2   ,ph2o   ,pco2  ,po2
C
      data frcsol / .001488, .001389, .001290, .001686, .002877,
     $              .003869, .026336, .426131, .526861, .526861,
     $              .526861, .526861, .526861, .526861, .526861,
     $              .006239, .001834, .001834/
C
      data wavmin / .200,  .245,  .265,  .275,  .285,
     $              .295,  .305,  .350,  .700,  .701,
     $              .701,  .701,  .701,  .702,  .702,
     $             2.630, 4.160, 4.160/
      data wavmax / .245,  .265,  .275,  .285,  .295,
     $              .305,  .350,  .700, 5.000, 5.000,
     $             5.000, 5.000, 5.000, 5.000, 5.000,
     $             2.860, 4.550, 4.550/
C
      data raytau / 4.020, 2.180, 1.700, 1.450, 1.250,
     $              1.085, 0.730, 0.135, 0.020, .0001,
     $              .0001, .0001, .0001, .0001, .0001,
     $              .0001, .0001, .0001/
C
C Absorption coefficiants
C
      data abh2o /    .000,     .000,    .000,    .000,    .000,
     $                .000,     .000,    .000,    .002,    .035,
     $                .377,    1.950,   9.400,  44.600, 190.000,
     $                .000,     .000,    .000/
C
      data abo3  /
     $ 5.370e+04, 13.080e+04,  9.292e+04, 4.530e+04, 1.616e+04,
     $ 4.441e+03,  1.775e+02,  2.101e+01,      .000,      .000,
     $  .000    ,   .000    ,   .000    ,      .000,      .000,
     $  .000    ,   .000    ,   .000    /
C
      data abco2  /    .000,     .000,    .000,    .000,    .000,
     $                 .000,     .000,    .000,    .000,    .000,
     $                 .000,     .000,    .000,    .000,    .000,
     $                 .094,     .196,   1.963/
C
      data abo2  /    .000,     .000,    .000,    .000,    .000,
     $                .000,     .000,1.11e-05,6.69e-05,    .000,
     $                .000,     .000,    .000,    .000,    .000,
     $                .000,    .000,    .000/
C
C Spectral interval weights
C
      data ph2o  /    .000,     .000,    .000,    .000,    .000,
     $                .000,     .000,    .000,    .505,    .210,
     $                .120,     .070,    .048,    .029,    .018,
     $                .000,     .000,    .000/
      data pco2  /    .000,     .000,    .000,    .000,    .000,
     $                .000,     .000,    .000,    .000,    .000,
     $                .000,     .000,    .000,    .000,    .000,
     $               1.000,     .640,    .360/
      data po2   /    .000,     .000,    .000,    .000,    .000,
     $                .000,     .000,   1.000,   1.000,    .000,
     $                .000,     .000,    .000,    .000,    .000,
     $                .000,     .000,    .000/
C
C Diagnostic and accumulation arrays; note that sfltot, fswup, and
C fswdn are not used in the computation,but are retained for future use.
C
      real solflx(plond),         ! Solar flux in current interval
     $     sfltot(plond),         ! Spectrally summed total solar flux
     $     totfld(plond,0:plev),  ! Spectrally summed flux divergence
     $     fswup(plond,0:plevp),  ! Spectrally summed up flux
     $     fswdn(plond,0:plevp)   ! Spectrally summed down flux
C
C Cloud radiative property arrays
C
      real tauxcl(plond,0:plev),  ! water cloud extinction optical depth
     $     tauxci(plond,0:plev),  ! ice cloud extinction optical depth
     $         wcl(plond,0:plev), ! liquid cloud single scattering albedo
     $         gcl(plond,0:plev), ! liquid cloud asymmetry parameter
     $         fcl(plond,0:plev), ! liquid cloud forward scattered fraction
     $         wci(plond,0:plev), ! ice cloud single scattering albedo
     $         gci(plond,0:plev), ! ice cloud asymmetry parameter
     $         fci(plond,0:plev)  ! ice cloud forward scattered fraction
C
C Aerosol radiative property arrays
C
      real tauxar(plond,0:plev),     ! aerosol extinction optical depth
     $         wa(plond,0:plev),     ! aerosol single scattering albedo
     $         ga(plond,0:plev),     ! aerosol assymetry parameter
     $         fa(plond,0:plev)      ! aerosol forward scattered fraction
C
      real tauaer(plond),            ! total column aerosol extinction
     $     waer(plond),              ! aerosol single scattering albedo
     $     gaer(plond),              ! aerosol asymmetry parameter
     $     faer(plond)               ! aerosol forward scattering fraction
C
C Sulphate aerosol properties from August 1992
C
      real ksa(nspint),   ! aerosol spectral mass absorption coeff (m2/g)
     $     wsa(nspint),   ! aerosol spectral single scattering albedo
     $     gsa(nspint)    ! aerosol spectral asymmetry parameter
C
      data ksa /11.1163, 10.5472, 10.2468, 10.0392,  9.8292,
     $           9.6199,  9.0407,  5.3012,  1.9169,  0.3780,
     $           0.3780,  0.3780,  0.3780,  0.5704,  0.5704,
     $           0.5704,  0.5704,  0.5704 /
C
      data wsa / .999999, .999999, .999999, .999999, .999999,
     $           .999999, .999999, .999999, .999991, .989772,
     $           .989772, .989772, .989772, .847061, .847061,
     $           .847061, .847061, .847061 /
C
      data gsa / .719161, .719012, .718453, .717820, .716997,
     $           .715974, .712743, .694889, .618115, .485286,
     $           .485286, .485286, .485286, .295557, .295557,
     $           .295557, .295557, .295557 /
C
C Other variables and arrays needed for aerosol:
C
      real rhfac,              ! multiplication factor for kaer
     $     rhpc,               ! level relative humidity in %
     $     a0,                 ! constant in rh mult factor
     $     a1,                 ! constant in rh mult factor
     $     a2,                 ! constant in rh mult factor
     $     a3                  ! constant in rh mult factor
c
      data a0 / -9.2906106183    /
      data a1 /  0.52570211505   /
      data a2 / -0.0089285760691 /
      data a3 /  5.0877212432e-05/
C
C Various arrays and other constants:
C
      real pflx(plond,0:plevp),   ! Interface press, including extra layer
     $     zenfac(plond),         ! Square root of cos solar zenith angle
     $     sqrco2,                ! Square root of the co2 mass mixg ratio
     $     tmp1,                  ! Temporary constant array
     $     tmp2,                  ! Temporary constant array
     $     pdel,                  ! Pressure difference across layer
     $     path,                  ! Mass path of layer
     $     ptop,                  ! Lower interface pressure of extra layer
     $     ptho2,                 ! Used to compute mass path of o2
     $     ptho3,                 ! Used to compute mass path of o3
     $     pthco2,                ! Used to compute mass path of co2
     $     pthh2o,                ! Used to compute mass path of h2o
     $     h2ostr,                ! Inverse square root h2o mass mixing ratio
     $     wavmid,                ! Spectral interval middle wavelength
     $     trayoslp               ! Rayleigh optical depth/standard pressure
      real tmp1l,                 ! Temporary constant array
     $     tmp2l,                 ! Temporary constant array
     $     tmp3l,                 ! Temporary constant array
     $     tmp1i,                 ! Temporary constant array
     $     tmp2i,                 ! Temporary constant array
     $     tmp3i                  ! Temporary constant array
      real rdenom,                ! Multiple scattering term
     $     psf,                   ! Frac of solar flux in spect interval
     $     gocp                   ! Gravity/cp
C
C Layer absorber amounts; note that 0 refers to the extra layer added
C above the top model layer
C
      real uh2o(plond,0:plev),    ! Layer absorber amount of h2o
     $      uo3(plond,0:plev),    ! Layer absorber amount of  o3
     $     uco2(plond,0:plev),    ! Layer absorber amount of co2
     $      uo2(plond,0:plev),    ! Layer absorber amount of  o2
     $      uaer(plond,0:plev)    ! Layer aerosol amount 
C
C Total column absorber amounts:
C
      real uth2o(plond),          ! Total column  absorber amount of  h2o
     $     uto3(plond),           ! Total column  absorber amount of  o3
     $     utco2(plond),          ! Total column  absorber amount of  co2
     $     uto2(plond),           ! Total column  absorber amount of  o2
     $     utaer(plond)           ! Total column  aerosol
C
C These arrays are defined for plev model layers; 0 refers to the extra
C layer on top:
C
      real rdir(plond,0:plev),    ! Layer reflectivity to direct rad
     $     rdif(plond,0:plev),    ! Layer reflectivity to diffuse rad
     $     tdir(plond,0:plev),    ! Layer transmission to direct rad
     $     tdif(plond,0:plev),    ! Layer transmission to diffuse rad
     $     explay(plond,0:plev),  ! Solar beam exp transmission for layer
     $     flxdiv(plond,0:plev)   ! Flux divergence for layer
C
C These arrays are defined at model interfaces; 0 is the top of the
C extra layer above the model top; plevp is the earth surface:
C
      real rupdir(plond,0:plevp), ! Ref to dir rad for layers below
     $     rupdif(plond,0:plevp), ! Ref to dif rad for layers below
     $     rdndif(plond,0:plevp), ! Ref to dif rad for layers above
     $     exptdn(plond,0:plevp), ! Solar beam exp down transm from top
     $     tottrn(plond,0:plevp), ! Total transmission for layers above
     $     fluxup(plond,0:plevp), ! Up   flux at model interface
     $     fluxdn(plond,0:plevp)  ! Down flux at model interface
C
C-----------------------------------------------------------------------
C
C Initialize output fields:
C
      do i=1, plon
         fsnt(i)  = 0.0
         fsns(i)  = 0.0
         solin(i) = 0.0
         fsnsc(i) = 0.0
         fsntc(i) = 0.0
         sols(i) = 0.0
         soll(i) = 0.0
         solsd(i) = 0.0
         solld(i) = 0.0
      end do
      do k=1, plev
         do i=1, plon
            qrs(i,k) = 0.0
         end do
      end do
C
C Compute starting, ending daytime loop indices:
C
      nloop = 0
      is(1) = isrchfgt(plon,coszrs,1,0.0)
C
C If night everywhere, return:
C
      if(is(1).gt.plon) return
      ie(1) = isrchfle(plon-is(1),coszrs(is(1)+1),1,0.0) + is(1) - 1
      nloop = 1
C
C Possibly 2 daytime loops needed:
C
      if(ie(1).ne.plon) then
         is(2) = isrchfgt(plon-ie(1),coszrs(ie(1)+1),1,0.0) + ie(1)
         if(is(2).le.plon) then
            nloop = 2
            ie(2) = plon
         end if
      end if
C
C Define solar incident radiation and interface pressures:
C
      do n=1,nloop
         do i=is(n),ie(n)
            solin(i) = scon*eccf*coszrs(i)
            pflx(i,0) = 0.
         end do
      end do
      do k=1,plevp
         do n=1,nloop
            do i=is(n),ie(n)
               pflx(i,k) = pint(i,k)
            end do
         end do
      end do
C
C Compute optical paths:
C
      tmp1   = 0.5/(gravit*sslp)
      co2mmr = co2vmr*(mmwco2/mmwair)
      sqrco2 = sqrt(co2mmr)
      do n=1,nloop
         do i=is(n),ie(n)
            ptop      = pflx(i,1)
            ptho2     = o2mmr * ptop / gravit
            ptho3     = o3mmr(i,1) * ptop / gravit
            pthco2    = sqrco2 * (ptop / gravit)
            h2ostr    = sqrt( 1. / h2ommr(i,1) )
            zenfac(i) = sqrt(coszrs(i))
            pthh2o    = ptop**2*tmp1 +
     $                (ptop*rga)*(h2ostr*zenfac(i)*delta)
            uh2o(i,0) = h2ommr(i,1)*pthh2o
            uco2(i,0) = zenfac(i)*pthco2
            uo2 (i,0) = zenfac(i)*ptho2
            uo3 (i,0) = ptho3
            uaer(i,0) = 0.0
         end do
      end do
C
      tmp2 = delta/gravit
      do k=1,plev
         do n=1,nloop
            do i=is(n),ie(n)
               pdel   = pflx(i,k+1) - pflx(i,k)
               path   = pdel / gravit
               ptho2  = o2mmr * path
               ptho3  = o3mmr(i,k) * path
               pthco2 = sqrco2 * path
               h2ostr = sqrt(1.0/h2ommr(i,k))
               pthh2o = (pflx(i,k+1)**2 - pflx(i,k)**2)*tmp1 +
     $                  pdel*h2ostr*zenfac(i)*tmp2
               uh2o(i,k) = h2ommr(i,k)*pthh2o
               uco2(i,k) = zenfac(i)*pthco2
               uo2 (i,k) = zenfac(i)*ptho2
               uo3 (i,k) = ptho3
C
C Adjust aerosol amount by relative humidity factor:
C
               if( rh(i,k) .gt. .90 ) then
                 rhfac = 2.8
               else if (rh(i,k) .lt. .60 ) then
                 rhfac = 1.0
               else
                 rhpc  = 100. * rh(i,k)
                 rhfac = (a0 + a1*rhpc + a2*rhpc**2 + a3*rhpc**3)
               endif
               uaer(i,k) = aermmr(i,k)*rhfac*path
C
            end do
         end do
      end do
C
C Compute column absorber amounts for the clear sky computation:
C
      do n=1,nloop
         do i=is(n),ie(n)
            uth2o(i) = 0.0
            uto3(i)  = 0.0
            utco2(i) = 0.0
            uto2(i)  = 0.0
            utaer(i) = 0.0
         end do
      end do
      do k=1,plev
         do n=1,nloop
            do i=is(n),ie(n)
               uth2o(i) = uth2o(i) + uh2o(i,k)
               uto3(i)  = uto3(i)  + uo3(i,k)
               utco2(i) = utco2(i) + uco2(i,k)
               uto2(i)  = uto2(i)  + uo2(i,k)
               utaer(i) = utaer(i) + uaer(i,k)
            end do
         end do
      end do
C
C Initialize spectrally integrated totals:
C
      do k=0,plev
         do i=1,plon
            totfld(i,k) = 0.0
            fswup (i,k) = 0.0
            fswdn (i,k) = 0.0
         end do
      end do
      do i=1,plon
         sfltot(i)       = 0.0
         fswup (i,plevp) = 0.0
         fswdn (i,plevp) = 0.0
      end do
C
C Set cloud properties for top (0) layer; so long as tauxcl is zero,
C there is no cloud above top of model; the other cloud properties
C are arbitrary:
C
      do n=1,nloop
         do i=is(n),ie(n)
            tauxcl(i,0) = 0.
            wcl(i,0)     = 0.999999
            gcl(i,0)     = 0.85
            fcl(i,0)     = 0.725
C
            tauxci(i,0) = 0.
            wci(i,0)     = 0.999999
            gci(i,0)     = 0.85
            fci(i,0)     = 0.725
C
C Aerosol 
C
            tauxar(i,0) = 0.
            wa(i,0)      = 0.925
            ga(i,0)      = 0.850
            fa(i,0)      = 0.7225
C
         end do
      end do
C
C Begin spectral loop
C
      do 100 ns=1,nspint
C
C Set index for cloud particle properties based on the wavelength,
C according to A. Slingo (1989) equations 1-3:
C Use index 1 (0.25 to 0.69 micrometers) for visible
C Use index 2 (0.69 - 1.19 micrometers) for near-infrared
C Use index 3 (1.19 to 2.38 micrometers) for near-infrared
C Use index 4 (2.38 to 4.00 micrometers) for near-infrared
C
C Note that the minimum wavelength is encoded (with .001, .002, .003)
C in order to specify the index appropriate for the near-infrared
C cloud absorption properties
C
        if(wavmax(ns) .le. 0.7) then
           indxsl = 1
        else if(wavmin(ns) .eq. 0.700) then
           indxsl = 2
        else if(wavmin(ns) .eq. 0.701) then
           indxsl = 3
        else if(wavmin(ns) .eq. 0.702 .or. wavmin(ns) .gt. 2.38) then
           indxsl = 4
        end if
C
C Set cloud extinction optical depth, single scatter albedo,
C asymmetry parameter, and forward scattered fraction:
C
        abarli = abarl(indxsl)
        bbarli = bbarl(indxsl)
        cbarli = cbarl(indxsl)
        dbarli = dbarl(indxsl)
        ebarli = ebarl(indxsl)
        fbarli = fbarl(indxsl)
c
        abarii = abari(indxsl)
        bbarii = bbari(indxsl)
        cbarii = cbari(indxsl)
        dbarii = dbari(indxsl)
        ebarii = ebari(indxsl)
        fbarii = fbari(indxsl)
        do k=1,plev

           do n=1,nloop
              do i=is(n),ie(n)
c   liquid
                 tmp1l = abarli + bbarli/rel(i,k)
                 tmp2l = 1. - cbarli - dbarli*rel(i,k)
                 tmp3l = fbarli*rel(i,k)
c   ice
                 tmp1i = abarii + bbarii/rei(i,k)
                 tmp2i = 1. - cbarii - dbarii*rei(i,k)
                 tmp3i = fbarii*rei(i,k)
C
C Cloud fraction incorporated into cloud extinction optical depth
C
                 tauxcl(i,k) = clwp(i,k)*tmp1l*(1.-fice(i,k))
     $                         *cld(i,k)*sqrt(cld(i,k))
                 tauxci(i,k) = clwp(i,k)*tmp1i*fice(i,k)
     $                         *cld(i,k)*sqrt(cld(i,k))
C
C Do not let single scatter albedo be 1; delta-eddington solution
C for non-conservative case:
C
                 wcl(i,k) = amin1(tmp2l,.999999)
                 gcl(i,k) = ebarli + tmp3l
                 fcl(i,k) = gcl(i,k)*gcl(i,k)
C
                 wci(i,k) = amin1(tmp2i,.999999)
                 gci(i,k) = ebarii + tmp3i
                 fci(i,k) = gci(i,k)*gci(i,k)
C
C Set aerosol properties:
C
                 tauxar(i,k) = 1.e4 * ksa(ns) * uaer(i,k)
C
                 wa(i,k)     = wsa(ns)
                 ga(i,k)     = gsa(ns)
                 fa(i,k)     = gsa(ns)*gsa(ns)
C
                 waer(i)     = wa(i,k)
                 gaer(i)     = ga(i,k)
                 faer(i)     = fa(i,k)
C
              end do
           end do
        end do
C
C Set reflectivities for surface based on mid-point wavelength
C
        wavmid = 0.5*(wavmin(ns) + wavmax(ns))
C
C Wavelength less  than 0.7 micro-meter
C
        if(wavmid .lt. 0.7 ) then
           do n=1,nloop
              do i=is(n),ie(n)
                 albdir(i) = albs(i)
                 albdif(i) = albsd(i)
              end do
           end do
C
C Wavelength greater than 0.7 micro-meter
C
        else
           do n=1,nloop
              do i=is(n),ie(n)
                 albdir(i) = albl(i)
                 albdif(i) = albld(i)
              end do
           end do
        end if
        trayoslp = raytau(ns)/sslp
C
C Layer input properties now completely specified; compute the
C delta-Eddington solution reflectivities and transmissivities
C for each layer, starting from the top and working downwards:
C
        call radded(coszrs   ,trayoslp,pflx   ,abh2o(ns),abo3(ns),
     $              abco2(ns),abo2(ns),uh2o   ,uo3      ,uco2    ,
     $              uo2      ,tauxcl  ,wcl    ,gcl      ,fcl     ,
     $              tauxci   ,wci     ,gci    ,fci      ,tauxar  ,
     $              wa       ,ga      ,fa     ,nloop    ,is      ,
     $              ie       ,rdir    ,rdif   ,tdir     ,tdif    ,
     $              explay   ,exptdn  ,rdndif ,tottrn   )
C
C Compute reflectivity to direct and diffuse radiation for layers below
C by adding succesive layers starting from the surface and working
C upwards:
C
        do n=1,nloop
           do i=is(n),ie(n)
              rupdir(i,plevp) = albdir(i)
              rupdif(i,plevp) = albdif(i)
           end do
        end do
        do k=plev,0,-1
           do n=1,nloop
              do i=is(n),ie(n)
                 rdenom = 1./( 1. - rdif(i,k)*rupdif(i,k+1))
                 rupdir(i,k) = rdir(i,k) + tdif(i,k)*
     $                 (rupdir(i,k+1)*explay(i,k) +
     $                  rupdif(i,k+1)*(tdir(i,k)-explay(i,k)))*rdenom
                 rupdif(i,k) = rdif(i,k) +
     $                         rupdif(i,k+1)*tdif(i,k)**2*rdenom
              end do
           end do
        end do
C
C Compute up and down fluxes for each interface, using the added
C atmospheric layer properties at each interface:
C
        do k=0,plevp
           do n=1,nloop
              do i=is(n),ie(n)
                 rdenom = 1./(1. - rdndif(i,k)*rupdif(i,k))
                 fluxup(i,k) = (exptdn(i,k)*rupdir(i,k) +
     $                  (tottrn(i,k)-exptdn(i,k))*rupdif(i,k))*rdenom
                 fluxdn(i,k)=exptdn(i,k) + (tottrn(i,k) - exptdn(i,k) +
     $                  exptdn(i,k)*rupdir(i,k)*rdndif(i,k))*rdenom
              end do
           end do
        end do
C
C Compute flux divergence in each layer using the interface up and down
C fluxes:
C
        do k=0,plev
           do n=1,nloop
              do i=is(n),ie(n)
                 flxdiv(i,k) = (fluxdn(i,k  ) - fluxdn(i,k+1)) +
     $                         (fluxup(i,k+1) - fluxup(i,k  ))
              end do
           end do
        end do
C
C Monochromatic computation completed; accumulate in totals; adjust
C fraction within spectral interval to allow for the possibility of
C sub-divisions within a particular interval:
C
        psf = 1.0
        if(ph2o(ns).ne.0.) psf = psf*ph2o(ns)
        if(pco2(ns).ne.0.) psf = psf*pco2(ns)
        if(po2 (ns).ne.0.) psf = psf*po2 (ns)
        do n=1,nloop
           do i=is(n),ie(n)
              solflx(i)  = solin(i)*frcsol(ns)*psf
              fsnt(i) = fsnt(i) + solflx(i)*(fluxdn(i,1) - fluxup(i,1))
              fsns(i) = fsns(i) + solflx(i)*
     $                  (fluxdn(i,plevp) - fluxup(i,plevp))
              sfltot(i)  = sfltot(i) + solflx(i)
              fswup(i,0) = fswup(i,0) + solflx(i)*fluxup(i,0)
              fswdn(i,0) = fswdn(i,0) + solflx(i)*fluxdn(i,0)
              if (wavmid .lt. 0.7) then
                 sols(i)=sols(i) + exptdn(i,plevp)*solflx(i)*0.001
                 solsd(i)=solsd(i) + (fluxdn(i,plevp) -
     $                    exptdn(i,plevp)) * solflx(i)*0.001 
              else
                 soll(i)=soll(i) + exptdn(i,plevp)*solflx(i)*0.001
                 solld(i)=solld(i) + (fluxdn(i,plevp) -
     $                    exptdn(i,plevp)) * solflx(i)*0.001 
c turn off print  23 Oct 2006 BPB
c..                 write(6,1233) ns,exptdn(i,plevp)
c..     $                         * solflx(i)*0.001 
c.. 1233            format(2x,'ns=',i3,1x,'flx_dir_dwn=',f8.4)
c..                 write(6,1234) ns,(fluxdn(i,plevp)-exptdn(i,plevp))
c..     $                         * solflx(i)*0.001 
c.. 1234            format(2x,'ns=',i3,1x,'flx_dif_dwn=',f8.4)
              end if
           end do
        end do
        do k=0,plev
           do n=1,nloop
              do i=is(n),ie(n)
                 totfld(i,k)  = totfld(i,k)  + solflx(i)*flxdiv(i,k)
                 fswup(i,k+1) = fswup(i,k+1) + solflx(i)*fluxup(i,k+1)
                 fswdn(i,k+1) = fswdn(i,k+1) + solflx(i)*fluxdn(i,k+1)
              end do
           end do
        end do
C
C
C Following code is the diagnostic clear sky computation:
C
C Compute delta-Eddington solution reflectivities and transmissivities
C for the entire column; note, for convenience, we use the same
C reflectivity and transmissivity arrays as for the full calculation
C above, where 0 for layer quantities refers to the entire atmospheric
C column, and where 0 for interface quantities refers to top of atmos-
C phere, while 1 refers to the surface:
C
C
C Compute total column aerosol optical depth:
C
        do n=1,nloop
           do i=is(n),ie(n)
             tauaer(i) = 1.e4 * ksa(ns) * utaer(i)
           enddo
        enddo
C
        call radclr(coszrs   ,trayoslp,pflx    ,abh2o(ns),abo3(ns) ,
     $              abco2(ns),abo2(ns),uth2o   ,uto3     ,utco2    ,
     $              uto2     ,tauaer  ,waer    ,gaer     ,faer     ,
     $              nloop    ,is      ,ie      ,rdir     ,rdif     ,
     $              tdir     ,tdif    ,explay  ,exptdn   ,rdndif   ,
     $              tottrn   )
C
C Compute reflectivity to direct and diffuse radiation for entire
C column; 0,1 on layer quantities refers to two effective layers
C overlying surface; 0 on interface quantities refers to top of column;
C 2 on interface quantities refers to the surface:
C
        do n=1,nloop
           do i=is(n),ie(n)
              rupdir(i,2) = albdir(i)
              rupdif(i,2) = albdif(i)
           end do
        end do
C
        do k=1,0,-1
           do n=1,nloop
              do i=is(n),ie(n)
                 rdenom = 1./( 1. - rdif(i,k)*rupdif(i,k+1))
                 rupdir(i,k) = rdir(i,k) + tdif(i,k)*
     $                 (rupdir(i,k+1)*explay(i,k) +
     $                  rupdif(i,k+1)*(tdir(i,k)-explay(i,k)))*rdenom
                 rupdif(i,k) = rdif(i,k) +
     $                        rupdif(i,k+1)*tdif(i,k)**2*rdenom
              end do
           end do
        end do
C
C Compute up and down fluxes for each interface, using the added
C atmospheric layer properties at each interface:
C
        do k=0,2
           do n=1,nloop
              do i=is(n),ie(n)
                 rdenom = 1./(1. - rdndif(i,k)*rupdif(i,k))
                 fluxup(i,k) = (exptdn(i,k)*rupdir(i,k) +
     $                  (tottrn(i,k)-exptdn(i,k))*rupdif(i,k))*rdenom
                 fluxdn(i,k)=exptdn(i,k) + (tottrn(i,k) - exptdn(i,k) +
     $                  exptdn(i,k)*rupdir(i,k)*rdndif(i,k))*rdenom
              end do
           end do
        end do
C
        do n=1,nloop
           do i=is(n),ie(n)
              fsntc(i) = fsntc(i) + solflx(i)*(fluxdn(i,0)-fluxup(i,0))
              fsnsc(i) = fsnsc(i) + solflx(i)*(fluxdn(i,2)-fluxup(i,2))
           end do
        end do
C
C End of clear sky calculation
C
  100 continue              ! End of spectral interval loop
C
C Compute solar heating rate (k/s)
C
      gocp = gravit/cpair
      do k=1,plev
         do n=1,nloop
           do i=is(n),ie(n)
              qrs(i,k) = -gocp*totfld(i,k)/(pint(i,k) - pint(i,k+1))
           end do
        end do
      end do
C
      return
      end
      subroutine radded(coszrs  ,trayoslp,pflx    ,abh2o   ,abo3    ,
     $                  abco2   ,abo2    ,uh2o    ,uo3     ,uco2    ,
     $                  uo2     ,tauxcl  ,wcl     ,gcl     ,fcl     ,
     $                  tauxci  ,wci     ,gci     ,fci     ,tauxar  ,
     $                  wa      ,ga      ,fa      ,nloop   ,is      ,
     $                  ie      ,rdir    ,rdif    ,tdir    ,tdif    ,
     $                  explay  ,exptdn  ,rdndif  ,tottrn  )
C-----------------------------------------------------------------------
C
C Computes layer reflectivities and transmissivities, from the top down
C to the surface using the delta-Eddington solutions for each layer;
C adds layers from top down to surface as well.
C
C If total transmission to the interface above a particular layer is
C less than trmin, then no further delta-Eddington solutions are
C evaluated for layers below
C
C For more details , see Briegleb, Bruce P., 1992: Delta-Eddington
C Approximation for Solar Radiation in the NCAR Community Climate Model,
C Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
C
C---------------------------Code history--------------------------------
C
C Original version:  B. Briegleb
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Kiehl, B. Briegleb, August 1992
C
C-----------------------------------------------------------------------
c
c $Id: radded.F,v 1.3 1995/03/17 18:54:11 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
C------------------------------Parameters-------------------------------
c
c $Id: prgrid.h,v 1.2 1995/02/10 01:09:10 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation resolution and I/O parameters
C
      integer
     $     plon,    ! number of longitudes
     $     plev,    ! number of vertical levels
     $     plat,    ! number of latitudes
     $     pcnst,   ! number of constituents (including water vapor)
     $     plevmx,  ! number of subsurface levels
     $     plevp,   ! plev + 1
     $     nxpt,    ! no.of points outside active domain for interpolant
     $     jintmx,  ! number of extra latitudes in polar region
     $     plond,   ! slt extended domain longitude
     $     platd,   ! slt extended domain lat.
     $     plevd    ! fold plev,pcnst indices into one
      integer
     $     plngbuf  ! length of absorptivity/emissivity record
C
      parameter(plon    = 1,
     $          plev    = 18,
     $          plat    = 1,
     $          pcnst   = 1,
     $          plevmx  = 4,
     $          plevp   = plev + 1,
     $          nxpt    = 1,
     $          jintmx  = 1,
     $          plond   = plon + 1 + 2*nxpt,
     $          platd   = plat + 2*nxpt + 2*jintmx,
     $          plevd   = plev*(3 + pcnst))
      parameter(plngbuf = 512*((plond*plevp*plevp + plond*plev*4 +
     $                          plond*plevp)/512 + 1))
C
C-----------------------------------------------------------------------
C
C Minimum total transmission below which no layer computation are done:
C
      real  trmin,          ! Minimum total transmission allowed
     $      wray,           ! Rayleigh single scatter albedo
     $      gray,           ! Rayleigh asymetry parameter
     $      fray            ! Rayleigh forward scattered fraction
      parameter (trmin = 1.e-3,
     $           wray = 0.999999,
     $           gray = 0.0,
     $           fray = 0.1)
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real coszrs(plond),         ! Cosine zenith angle
     $     trayoslp,              ! Tray/sslp
     $     pflx(plond,0:plevp),   ! Interface pressure
     $     abh2o,                 ! Absorption coefficiant for h2o
     $     abo3 ,                 ! Absorption coefficiant for o3
     $     abco2,                 ! Absorption coefficiant for co2
     $     abo2 ,                 ! Absorption coefficiant for o2
     $     uh2o(plond,0:plev),    ! Layer absorber amount of h2o
     $     uo3(plond,0:plev),     ! Layer absorber amount of  o3
     $     uco2(plond,0:plev),    ! Layer absorber amount of co2
     $     uo2(plond,0:plev)      ! Layer absorber amount of  o2
      real tauxcl(plond,0:plev),  ! Cloud extinction optical depth
     $     wcl(plond,0:plev),      ! Cloud single scattering albedo
     $     gcl(plond,0:plev),      ! Cloud assymetry parameter
     $     fcl(plond,0:plev),       ! Cloud forward scattered fraction
     $     tauxci(plond,0:plev),  ! Cloud extinction optical depth
     $     wci(plond,0:plev),      ! Cloud single scattering albedo
     $     gci(plond,0:plev),      ! Cloud assymetry parameter
     $     fci(plond,0:plev)
      real tauxar(plond,0:plev), ! Aerosol extinction optical depth
     $     wa(plond,0:plev),     ! Aerosol single scattering albedo
     $     ga(plond,0:plev),     ! Aerosol assymetry parameter
     $     fa(plond,0:plev)      ! Aerosol forward scattered fraction
      integer nloop,              ! Number of loops (1 or 2)
     $        is(2),              ! Starting index for 1 or 2 loops
     $        ie(2)               ! Ending index for 1 or 2 loops
C
C Input/Output arguments
C
C Following variables are defined for each layer; 0 refers to extra
C layer above top of model:
C
      real rdir(plond,0:plev),    ! Layer reflectivity to direct rad
     $     rdif(plond,0:plev),    ! Layer refflectivity to diffuse rad
     $     tdir(plond,0:plev),    ! Layer transmission to direct rad
     $     tdif(plond,0:plev),    ! Layer transmission to diffuse rad
     $     explay(plond,0:plev)   ! Solar beam exp transm for layer
C
C (Note that the following variables are defined on interfaces, with the
C  index k referring to the top interface of the kth layer:
C  exptdn,rdndif,tottrn; for example, tottrn(k=5) refers to the total
C  transmission to the top interface of the 5th layer; plevp refers to
C  the earth surface)
C
      real rdndif(plond,0:plevp), ! Added dif ref for layers above
     $     exptdn(plond,0:plevp), ! Solar beam exp down transm from top
     $     tottrn(plond,0:plevp)  ! Total transmission for layers above
C
C------------------------------Externals--------------------------------
C
      external  resetr,           ! Resets array elements to zero
     $          whenfgt           ! Collect indices greater than conditn
C
C---------------------------Local variables-----------------------------
C
      integer i,            ! Longitude index
     $        k,            ! Level index
     $        nn,           ! Index of longitude loops (max=nloop)
     $        ii,           ! Longitude index
     $        nval,         ! Number of long values satisfying criteria
     $        index(plond)  ! Array of longitude indices
C
      real taugab(plond),   ! Layer total gas absorption optical depth
     $     tauray(plond),   ! Layer rayleigh optical depth
     $     taucsc,          ! Layer cloud scattering optical depth
     $     tautot,          ! Total layer optical depth
     $     wtot,            ! Total layer single scatter albedo
     $     gtot,            ! Total layer asymmetry parameter
     $     ftot             ! Total layer forward scatter fraction
C
      real wtau,            !  rayleigh layer scattering optical depth
     $     wt,              !  layer total single scattering albedo
     $     ts,              !  layer scaled extinction optical depth
     $     ws,              !  layer scaled single scattering albedo
     $     gs               !  layer scaled asymmetry parameter
C
      real rdenom,          !  mulitiple scattering term
     $     rdirexp,         !  layer direct ref times exp transmission
     $     tdnmexp          !  total transmission minus exp transmission
C
C---------------------------Statement functions-------------------------
C
C Statement functions and other local variables
C
      real alpha,           ! Term in direct reflect and transmissivity
     $     gamma,           ! Term in direct reflect and transmissivity
     $     el,              ! Term in alpha,gamma,n,u
     $     taus,            ! Scaled extinction optical depth
     $     omgs,            ! Scaled single particle scattering albedo
     $     asys,            ! Scaled asymmetry parameter
     $     u,               ! Term in diffuse reflect and transmissivity
     $     n,               ! Term in diffuse reflect and transmissivity
     $     lm,              ! Temporary for el
     $     ne               ! Temporary for n
      real w,               ! Dummy argument for statement function
     $     uu,              ! Dummy argument for statement function
     $     g,               ! Dummy argument for statement function
     $     e,               ! Dummy argument for statement function
     $     f,               ! Dummy argument for statement function
     $     t,               ! Dummy argument for statement function
     $     et               ! Dummy argument for statement function
C
C Intermediate terms for delta-eddington solution
C
      real alp,             ! Temporary for alpha
     $     gam,             ! Temporary for gamma
     $     ue,              ! Temporary for u
     $     arg,             ! Exponential argument
     $     extins,          ! Extinction
     $     amg,             ! Alp - gam
     $     apg              ! Alp + gam
C
      alpha(w,uu,g,e) = .75*w*uu*((1. + g*(1-w))/(1. - e*e*uu*uu))
      gamma(w,uu,g,e) = .50*w*((3.*g*(1.-w)*uu*uu + 1.)/(1.-e*e*uu*uu))
      el(w,g)         = sqrt(3.*(1-w)*(1. - w*g))
      taus(w,f,t)     = (1. - w*f)*t
      omgs(w,f)       = (1. - f)*w/(1. - w*f)
      asys(g,f)       = (g - f)/(1. - f)
      u(w,g,e)        = 1.5*(1. - w*g)/e
      n(uu,et)        = ((uu+1.)*(uu+1.)/et ) - ((uu-1.)*(uu-1.)*et)
C
C-----------------------------------------------------------------------
C
C Initialize all total transmission values to 0, so that nighttime 
C values from previous computations are not used:
C
      call resetr(tottrn,plond*plevp,0.)
C
C Compute total direct beam transmission, total transmission, and
C reflectivity for diffuse radiation (from below) for all layers above
C each interface by starting from the top and adding layers down:
C
C For the extra layer above model top:
C
      do 200 nn=1,nloop
         do 100 i=is(nn),ie(nn)
C
            tauray(i) = trayoslp*(pflx(i,1)-pflx(i,0))
            taugab(i) = abh2o*uh2o(i,0) + abo3*uo3(i,0) +
     $                  abco2*uco2(i,0) + abo2*uo2(i,0)
C
            tautot  = tauxcl(i,0) + tauxci(i,0) + tauray(i) + taugab(i)
     $                            + tauxar(i,0)
            taucsc  = tauxcl(i,0)*wcl(i,0)+tauxci(i,0)*wci(i,0)
     $                            + tauxar(i,0)*wa(i,0)
            wtau    = wray*tauray(i) 
            wt      = wtau + taucsc
            wtot = wt/tautot
            gtot = (wtau*gray + gcl(i,0)*tauxcl(i,0)*wcl(i,0) +
     $                          gci(i,0)*tauxci(i,0)*wci(i,0) +
     $                          ga(i,0) *tauxar(i,0)*wa(i,0))/wt
            ftot = (wtau*fray + fcl(i,0)*tauxcl(i,0)*wcl(i,0) +
     $                          fci(i,0)*tauxci(i,0)*wci(i,0) +
     $                          fa(i,0) *tauxar(i,0)*wa(i,0))/wt
C
            ts   = taus(wtot,ftot,tautot)
            ws   = omgs(wtot,ftot)
            gs   = asys(gtot,ftot)
            lm   = el(ws,gs)
            alp  = alpha(ws,coszrs(i),gs,lm)
            gam  = gamma(ws,coszrs(i),gs,lm)
            ue   = u(ws,gs,lm)
C
C Limit argument of exponential to 25, in case lm*ts very large:
C
            arg  = amin1(lm*ts,25.)
            extins = exp(-arg)
            ne = n(ue,extins)
C
            rdif(i,0) = (ue+1.)*(ue-1.)*(1./extins - extins)/ne
            tdif(i,0) = 4.*ue/ne
C
C Limit argument of exponential to 25, in case coszrs is very small:
C
            arg       = amin1(ts/coszrs(i),25.)
            explay(i,0) = exp(-arg)
C
            apg = alp + gam
            amg = alp - gam
            rdir(i,0) = amg*(tdif(i,0)*explay(i,0) - 1.) + apg*rdif(i,0)
            tdir(i,0) = apg*tdif(i,0) +
     $                  (amg*rdif(i,0) - (apg-1.))*explay(i,0)
C
C Under rare conditions, reflectivies and transmissivities can be
C negative; zero out any negative values
C
            rdir(i,0) = amax1(rdir(i,0),0.0)
            tdir(i,0) = amax1(tdir(i,0),0.0)
            rdif(i,0) = amax1(rdif(i,0),0.0)
            tdif(i,0) = amax1(tdif(i,0),0.0)
C
C Initialize top interface of extra layer:
C
            exptdn(i,0) =   1.0
            rdndif(i,0) =   0.0
            tottrn(i,0) =   1.0
C
            rdndif(i,1) = rdif(i,0)
            tottrn(i,1) = tdir(i,0)
C
  100    continue
  200 continue
C
C Now, continue down one layer at a time; if the total transmission to
C the interface just above a given layer is less than trmin, then no
C delta-eddington computation for that layer is done:
C
      do 400 k=1,plev
C
C Initialize current layer properties to zero; only if total
C transmission to the top interface of the current layer exceeds the
C minimum, will these values be computed below:
C
         do nn=1,nloop
            do i=is(nn),ie(nn)
C
               rdir(i,k)   =  0.0
               rdif(i,k)   =  0.0
               tdir(i,k)   =  0.0
               tdif(i,k)   =  0.0
               explay(i,k) =  0.0
C
C Calculates the solar beam transmission, total transmission, and
C reflectivity for diffuse radiation from below at the top of the
C current layer:
C
               exptdn(i,k) = exptdn(i,k-1)*explay(i,k-1)
               rdenom      = 1./(1. - rdif(i,k-1)*rdndif(i,k-1))
               rdirexp     = rdir(i,k-1)*exptdn(i,k-1)
               tdnmexp     = tottrn(i,k-1) - exptdn(i,k-1)
               tottrn(i,k) = exptdn(i,k-1)*tdir(i,k-1) + tdif(i,k-1)*
     $                      (tdnmexp + rdndif(i,k-1)*rdirexp)*rdenom
               rdndif(i,k) = rdif(i,k-1)  +
     $                (rdndif(i,k-1)*tdif(i,k-1))*(tdif(i,k-1)*rdenom)
C
            end do
         end do
C
C Compute next layer delta-eddington solution only if total transmission
C of radiation to the interface just above the layer exceeds trmin.
C
         call whenfgt(plon,tottrn(1,k),1,trmin,index,nval)
         if(nval.gt.0) then
CDIR$ IVDEP
            do 300 ii=1,nval
               i=index(ii)
C
               tauray(i) = trayoslp*(pflx(i,k+1)-pflx(i,k))
               taugab(i) = abh2o*uh2o(i,k) + abo3*uo3(i,k) +
     $                     abco2*uco2(i,k) + abo2*uo2(i,k)
C
               tautot = tauxcl(i,k) + tauxci(i,k) + 
     $                  tauray(i) + taugab(i) + tauxar(i,k)
               taucsc    = tauxcl(i,k)*wcl(i,k) + tauxci(i,k)*wci(i,k)
     $                            + tauxar(i,k)*wa(i,k)
               wtau      = wray*tauray(i)
               wt        = wtau + taucsc
               wtot   = wt/tautot
               gtot   = (wtau*gray + gcl(i,k)*wcl(i,k)*tauxcl(i,k)
     $                   + gci(i,k)*wci(i,k)*tauxci(i,k)
     $                   + ga(i,k) *wa(i,k) *tauxar(i,k))/wt
               ftot   = (wtau*fray + fcl(i,k)*wcl(i,k)*tauxcl(i,k)
     $                   + fci(i,k)*wci(i,k)*tauxci(i,k)
     $                   + fa(i,k) *wa(i,k) *tauxar(i,k))/wt
C
               ts   = taus(wtot,ftot,tautot)
               ws   = omgs(wtot,ftot)
               gs   = asys(gtot,ftot)
               lm   = el(ws,gs)
               alp  = alpha(ws,coszrs(i),gs,lm)
               gam  = gamma(ws,coszrs(i),gs,lm)
               ue   = u(ws,gs,lm)
C
C Limit argument of exponential to 25, in case lm very large:
C
               arg  = amin1(lm*ts,25.)
               extins = exp(-arg)
               ne = n(ue,extins)
C
               rdif(i,k) = (ue+1.)*(ue-1.)*(1./extins - extins)/ne
               tdif(i,k)   =   4.*ue/ne
C
C Limit argument of exponential to 25, in case coszrs is very small:
C
               arg       = amin1(ts/coszrs(i),25.)
               explay(i,k) = exp(-arg)
C
               apg = alp + gam
               amg = alp - gam
               rdir(i,k) = amg*(tdif(i,k)*explay(i,k) - 1.) +
     $                     apg*rdif(i,k)
               tdir(i,k) = apg*tdif(i,k) +
     $                     (amg*rdif(i,k) - (apg-1.))*explay(i,k)
C
C Under rare conditions, reflectivies and transmissivities can be
C negative; zero out any negative values
C
               rdir(i,k) = amax1(rdir(i,k),0.0)
               tdir(i,k) = amax1(tdir(i,k),0.0)
               rdif(i,k) = amax1(rdif(i,k),0.0)
               tdif(i,k) = amax1(tdif(i,k),0.0)
  300       continue
         end if
C
  400 continue
C
C Compute total direct beam transmission, total transmission, and
C reflectivity for diffuse radiation (from below) for all layers
C above the surface:
C
      k = plevp
      do nn=1,nloop
         do i=is(nn),ie(nn)
            exptdn(i,k) = exptdn(i,k-1)*explay(i,k-1)
            rdenom = 1./(1. - rdif(i,k-1)*rdndif(i,k-1))
            rdirexp = rdir(i,k-1)*exptdn(i,k-1)
            tdnmexp = tottrn(i,k-1) - exptdn(i,k-1)
            tottrn(i,k) = exptdn(i,k-1)*tdir(i,k-1) + tdif(i,k-1)*
     $                   (tdnmexp + rdndif(i,k-1)*rdirexp)*rdenom
            rdndif(i,k) = rdif(i,k-1)  +
     $               (rdndif(i,k-1)*tdif(i,k-1))*(tdif(i,k-1)*rdenom)
         end do
      end do
C
      return
      end

      subroutine radems(s2c     ,s2t     ,w       ,tplnke  ,plh2o   ,
     $                  pnm     ,plco2   ,tint    ,tint4   ,tlayr   ,
     $                  tlayr4  ,plol    ,plos    ,ucfc11  ,ucfc12  , 
     $                  un2o0   ,un2o1   ,uch4    ,uco211  ,uco212  ,
     $                  uco213  ,uco221  ,uco222  ,uco223  ,uptype  ,
     $                  bn2o0   ,bn2o1   ,bch4    ,co2em   ,co2eml  ,
     $                  co2t    ,h2otr   ,abplnk1 ,abplnk2 ,emstot  )
C-----------------------------------------------------------------------
C
C Compute emissivity for H2O, CO2, O3
C
C H2O  ....  Uses nonisothermal emissivity for water vapor from
C            Ramanathan, V. and  P.Downey, 1986: A Nonisothermal
C            Emissivity and Absorptivity Formulation for Water Vapor
C            Jouranl of Geophysical Research, vol. 91., D8, pp 8649-8666
C
C
C CO2  ....  Uses absorptance parameterization of the 15 micro-meter
C            (500 - 800 cm-1) band system of Carbon Dioxide, from
C            Kiehl, J.T. and B.P.Briegleb, 1991: A New Parameterization
C            of the Absorptance Due to the 15 micro-meter Band System
C            of Carbon Dioxide Jouranl of Geophysical Research,
C            vol. 96., D5, pp 9013-9019
C
C O3   ....  Uses absorptance parameterization of the 9.6 micro-meter
C            band system of ozone, from Ramanathan, V. and R. Dickinson,
C            1979: The Role of stratospheric ozone in the zonal and
C            seasonal radiative energy balance of the earth-troposphere
C            system. Journal of the Atmospheric Sciences, Vol. 36,
C            pp 1084-1104
C
C Computes individual emissivities, accounting for band overlap, and
C sums to obtain the total.
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Kiehl, B. Briegleb, August 1992
C
C-----------------------------------------------------------------------
c
c $Id: radems.F,v 1.2 1995/02/17 21:28:41 jhack Exp $
c $Author: jhack $
c
      implicit none
C------------------------------Parameters-------------------------------
c
c $Id: prgrid.h,v 1.2 1995/02/10 01:09:10 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation resolution and I/O parameters
C
      integer
     $     plon,    ! number of longitudes
     $     plev,    ! number of vertical levels
     $     plat,    ! number of latitudes
     $     pcnst,   ! number of constituents (including water vapor)
     $     plevmx,  ! number of subsurface levels
     $     plevp,   ! plev + 1
     $     nxpt,    ! no.of points outside active domain for interpolant
     $     jintmx,  ! number of extra latitudes in polar region
     $     plond,   ! slt extended domain longitude
     $     platd,   ! slt extended domain lat.
     $     plevd    ! fold plev,pcnst indices into one
      integer
     $     plngbuf  ! length of absorptivity/emissivity record
C
      parameter(plon    = 1,
     $          plev    = 18,
     $          plat    = 1,
     $          pcnst   = 1,
     $          plevmx  = 4,
     $          plevp   = plev + 1,
     $          nxpt    = 1,
     $          jintmx  = 1,
     $          plond   = plon + 1 + 2*nxpt,
     $          platd   = plat + 2*nxpt + 2*jintmx,
     $          plevd   = plev*(3 + pcnst))
      parameter(plngbuf = 512*((plond*plevp*plevp + plond*plev*4 +
     $                          plond*plevp)/512 + 1))
C
C------------------------------Commons----------------------------------
c
c $Id: crdcae.h,v 1.1.1.1 1995/02/09 23:26:46 ccm2 Exp $
c $Author: ccm2 $
c
C
C Water vapor narrow band constants for longwave radiation computations
C
      common/crdcae/realk(2), st(2), a1(2), a2(2), b1(2), b2(2),
     $              coefa(3,4),coefb(4,4),coefc(3,4),coefd(4,4),
     $              coefe(3,4),coeff(6,2),coefg(2,4),coefh(2,4),
     $              coefi(6,2),coefj(3,2),coefk(3,2),
     $              c1(4),c2(4),c3(4),c4(4),c5(4),c6(4),c7(4),c8,c9,
     $              c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,
     $              c22,c23,c24,c25,c26,c27,c28,c29,c30,c31,
     $              fwcoef,fwc1,fwc2,fc1,cfa1
C
      real realk,     ! H2O narrow band parameter
     $     st,        ! H2O narrow band parameter
     $     a1,a2,     ! Temperature correction terms for H2O path
     $     b1,b2      ! Temperature correction terms for H2O path
C
C Constant coefficients for water vapor absorptivity and emissivity
C
      real coefa,coefb,coefc,coefd,coefe,coeff,
     $     coefg,coefh,coefi,coefj,coefk,
     $      c1, c2, c3, c4, c5, c6, c7,c8 ,c9 ,c10,
     $     c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,
     $     c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,c31
C
C Farwing correction constants for narrow-band emissivity model,
C introduced to account for the deficiencies in narrow-band model
C used to derive the emissivity; tuned with Arking's line-by-line
C calculations.
C
      real fwcoef,
     $     fwc1,fwc2,
     $     fc1,
     $     cfa1
C
C-----------------------------------------------------------------------
c
c $Id: crdcon.h,v 1.1.1.1 1995/02/09 23:26:46 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation constants
C
      common/crdcon/gravit  ,rga     ,cpair   ,epsilo  ,sslp    ,
     $              stebol  ,rgsslp  ,co2vmr  ,dpfo3   ,dpfco2  ,
     $              dayspy  ,pie
C
      real gravit,    ! Acceleration of gravity
     $     rga,       ! 1./gravit
     $     cpair,     ! Specific heat of dry air
     $     epsilo,    ! Ratio of mol. wght of H2O to dry air
     $     sslp,      ! Standard sea-level pressure
     $     stebol,    ! Stefan-Boltzmann's constant
     $     rgsslp,    ! 0.5/(gravit*sslp)
     $     co2vmr,    ! CO2 volume mixing ratio
     $     dpfo3,     ! Voigt correction factor for O3
     $     dpfco2,    ! Voigt correction factor for CO2
     $     dayspy,    ! Number of days per 1 year
     $     pie        ! 3.14.....
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real s2c(plond,plevp),    ! H2o continuum path length
     $     s2t(plond,plevp),    ! Tmp and prs wghted h2o path length
     $     w(plond,plevp),      ! H2o path length
     $     tplnke(plond),       ! Layer planck temperature
     $     plh2o(plond,plevp),  ! H2o prs wghted path length
     $     pnm(plond,plevp),    ! Model interface pressure
     $     plco2(plond,plevp),  ! Prs wghted path of co2
     $     tint(plond,plevp),   ! Model interface temperatures
     $     tint4(plond,plevp),  ! Tint to the 4th power
     $     tlayr(plond,plevp),  ! K-1 model layer temperature
     $     tlayr4(plond,plevp), ! Tlayr to the 4th power
     $     plol(plond,plevp),   ! Pressure wghtd ozone path
     $     plos(plond,plevp)    ! Ozone path
c
c   Trace gas variables
c
      real ucfc11(plond,plevp), ! CFC11 path length
     $     ucfc12(plond,plevp), ! CFC12 path length
     $     un2o0(plond,plevp),  ! N2O path length
     $     un2o1(plond,plevp),  ! N2O path length (hot band)
     $     uch4(plond,plevp),   ! CH4 path length
     $     uco211(plond,plevp), ! CO2 9.4 micron band path length
     $     uco212(plond,plevp), ! CO2 9.4 micron band path length
     $     uco213(plond,plevp), ! CO2 9.4 micron band path length
     $     uco221(plond,plevp), ! CO2 10.4 micron band path length
     $     uco222(plond,plevp), ! CO2 10.4 micron band path length
     $     uco223(plond,plevp), ! CO2 10.4 micron band path length
     $     bn2o0(plond,plevp),  ! pressure factor for n2o
     $     bn2o1(plond,plevp),  ! pressure factor for n2o
     $     bch4(plond,plevp),   ! pressure factor for ch4
     $     uptype(plond,plevp)  ! p-type continuum path length
C
C Output arguments
C
      real emstot(plond,plevp),  ! Total emissivity
     $     co2em(plond,plevp),   ! Layer co2 normalzd plnck funct drvtv
     $     co2eml(plond,plev),   ! Intrfc co2 normalzd plnck func drvtv
     $     co2t(plond,plevp),    ! Tmp and prs weighted path length
     $     h2otr(plond,plevp)    ! H2o transmission over o3 band
      real emplnk(14,plond),        ! emissivity Planck factor
     $     abplnk1(14,plond,plevp), ! non-nearest layer Plack factor
     $     abplnk2(14,plond,plevp)  ! nearest layer factor
      real emstrc(plond,plevp)  ! total trace gas emissivity

C
C---------------------------Local variables-----------------------------
C
      integer
     $     i,                  ! Longitude index
     $     k,                  ! Level index]
     $     k1,                 ! Level index
     $     iband               ! H2o band index
C
C Local variables for H2O:
C
      real h2oems(plond,plevp),! H2o emissivity
     $     tpathe(plond),      ! Used to compute h2o emissivity
     $     a(plond),           ! Eq(2) in table A3a of R&D
     $     corfac(plond),      ! Correction factors in table A3b
     $     dtp(plond),         ! Path temperature minus 300 K used in 
C                                h2o rotation band absorptivity
     $     dtx(plond),         ! Planck temperature minus 250 K
     $     dty(plond),         ! Path temperature minus 250 K
     $     dtz(plond),         ! Planck temperature minus 300 K
     $     emis(plond,4),      ! Total emissivity (h2o+co2+o3)
     $     rsum(plond),        ! Eq(1) in table A2 of R&D
     $     term1(plond,4),     ! Equation(5) in table A3a of R&D(1986)
     $     term2(plond,4)      ! Delta a(Te) in table A3a of R&D(1986)
      real term3(plond,4),     ! B(T) function for rotation and
C                                vibration-rotation band emissivity
     $     term4(plond,4),     ! Equation(6) in table A3a of R&D(1986)
     $     term5(plond,4),     ! Delta a(Tp) in table A3a of R&D(1986)
     $     term6(plond,2),     ! B(T) function for window region
     $     term7(plond,2),     ! Kl_inf(i) in eq(8) of table A3a of R&D
     $     term8(plond,2),     ! Delta kl_inf(i) in eq(8)
     $     term9(plond,2),     ! B(T) function for 500-800 cm-1 region
     $     tr1(plond),         ! Equation(6) in table A2 for 650-800
     $     tr2(plond),         ! Equation(6) in table A2 for 500-650
     $     tr3(plond)          ! Equation(4) in table A2 for 650-800
      real tr4(plond),         ! Equation(4),table A2 of R&D for 500-650
     $     tr7(plond),         ! Equation (6) times eq(4) in table A2
C                              !   of R&D for 650-800 cm-1 region
     $     tr8(plond),         ! Equation (6) times eq(4) in table A2
C                              !   of R&D for 500-650 cm-1 region
     $     uc(plond),          ! Y + 0.002U in eq(8) of table A2 of R&D
     $     pnew(plond),        ! Effective pressure for h2o linewidth
     $     trline(plond,2),    ! Transmission due to H2O lines in window
     $     k21(plond),         ! Exponential coefficient used to calc
C                              !  rot band transmissivity in the 650-800
C                              !  cm-1 region (tr1)
     $     k22(plond),         ! Exponential coefficient used to calc
C                              !  rot band transmissivity in the 500-650
C                              !  cm-1 region (tr2)
     $     u(plond),           ! Pressure weighted H2O path length
     $     uc1(plond),         ! H2o continuum pathlength 500-800 cm-1
     $     r80257              ! Conversion factor for h2o pathlength
      real a11,                ! A1 in table A3b for rotation band emiss
     $     a31,                ! A3 in table A3b for rotation band emiss
     $     a21,                ! First part in numerator of A2 table A3b
     $     a22,                ! Second part in numertor of A2 table A3b
     $     a23,                ! Denominator of A2 table A3b (rot band)
     $     t1t4,               ! Eq(3) in table A3a of R&D
     $     t2t5,               ! Eq(4) in table A3a of R&D
     $     fwk,                ! Equation(33) in R&D far wing correction
     $     a41,                ! Numerator in A2 in Vib-rot (table A3b)
     $     a51,                ! Denominator in A2 in Vib-rot(table A3b)
     $     a61,                ! A3 factor for Vib-rot band in table A3b
     $     phi,                ! Eq(11) in table A3a of R&D
     $     psi,                ! Eq(12) in table A3a of R&D
     $     ubar,               ! H2o scaled path comment eq(10) table A2
     $     g1,                 ! Part of eq(10) table A2
     $     pbar,               ! H2o scaled pres comment eq(10) table A2
     $     g3,                 ! Part of eq(10) table A2
     $     g2,                 ! Part of arguement in eq(10) in table A2
     $     g4,                 ! Arguement in exp() in eq(10) table A2
     $     cf812               ! Eq(11) in table A2 of R&D
      real troco2(plond,plevp) ! H2o overlap factor for co2 absorption
C
C Local variables for CO2:
C
      real co2ems(plond,plevp), ! Co2 emissivity
     $     co2plk(plond),       ! Used to compute co2 emissivity
     $     sum(plond),          ! Used to calculate path temperature
     $     t1i,                 ! Co2 hot band temperature factor
     $     sqti,                ! Sqrt of temperature
     $     pi,                  ! Pressure used in co2 mean line width
     $     et,                  ! Co2 hot band factor
     $     et2,                 ! Co2 hot band factor
     $     et4,                 ! Co2 hot band factor
     $     omet,                ! Co2 stimulated emission term
     $     ex                   ! Part of co2 planck function
      real f1co2,               ! Co2 weak band factor
     $     f2co2,               ! Co2 weak band factor
     $     f3co2,               ! Co2 weak band factor
     $     t1co2,               ! Overlap factor weak bands strong band
     $     sqwp,                ! Sqrt of co2 pathlength
     $     f1sqwp,              ! Main co2 band factor
     $     oneme,               ! Co2 stimulated emission term
     $     alphat,              ! Part of the co2 stimulated emiss term
     $     wco2,                ! Consts used to define co2 pathlength
     $     posqt,               ! Effective pressure for co2 line width
     $     rbeta7,              ! Inverse of co2 hot band line width par
     $     rbeta8,              ! Inverse of co2 hot band line width par
     $     rbeta9,              ! Inverse of co2 hot band line width par
     $     rbeta13              ! Inverse of co2 hot band line width par
      real tpath,               ! Path temp used in co2 band model
     $     tmp1,                ! Co2 band factor
     $     tmp2,                ! Co2 band factor
     $     tmp3,                ! Co2 band factor
     $     tlayr5,              ! Temperature factor in co2 Planck func
     $     rsqti,               ! Reciprocal of sqrt of temperature
     $     exm1sq               ! Part of co2 Planck function
      real u7,    ! Absorber amount for various co2 band systems
     $     u8,    ! Absorber amount for various co2 band systems
     $     u9,    ! Absorber amount for various co2 band systems
     $     u13    ! Absorber amount for various co2 band systems
      real r250,  ! Inverse 250K
     $     r300,  ! Inverse 300K
     $     rsslp  ! Inverse standard sea-level pressure
C
C Local variables for O3:
C
      real o3ems(plond,plevp),  ! Ozone emissivity
     $     dbvtt(plond),        ! Tmp drvtv of planck fctn for tplnke
     $     te,                  ! Temperature factor
     $     u1,                  ! Path length factor
     $     u2,                  ! Path length factor
     $     phat,                ! Effecitive path length pressure
     $     tlocal,              ! Local planck function temperature
     $     tcrfac,              ! Scaled temperature factor
     $     beta,                ! Absorption funct factor voigt effect
     $     realnu,              ! Absorption function factor
     $     o3bndi               ! Band absorption factor
C
C Transmission terms for various spectral intervals:
C
      real trem1(plond),        ! H2o     0 -  800 cm-1
     $     trem2(plond),        ! H2o   500 -  800 cm-1
     $     trem3(plond),        ! Co2   500 -  800 cm-1
     $     trem4(plond),        ! H2o   800 - 1000 cm-1
     $     trem5(plond),        ! O3     9.6 micro-meter band
     $     trem6(plond),        ! H2o  1000 - 1200 cm-1
     $     trem7(plond)         ! H2o  1200 - 2200 cm-1
      real bndfct,              ! Band absorptance parameter for co2
     $     absbnd               ! Proportional to co2 band absorptance
      real tco2(plond),         ! co2 overlap factor
     $     th2o(plond),         ! h2o overlap factor
     $     to3(plond)           ! o3 overlap factor
C
C---------------------------Statement functions-------------------------
C
C Statement functions
C Derivative of planck function at 9.6 micro-meter wavelength, and
C an absorption function factor:
C
      real dbvt,fo3,t,ux,vx
C
      dbvt(t)=(-2.8911366682e-4+(2.3771251896e-6+1.1305188929e-10*t)*t)/
     $  (1.0+(-6.1364820707e-3+1.5550319767e-5*t)*t)
C
      fo3(ux,vx)=ux/sqrt(4.+ux*(1.+vx))
C
C-----------------------------------------------------------------------
C
C Initialize
C
      r80257  = 1./8.0257e-04
C
      r250  = 1./250.
      r300  = 1./300.
      rsslp = 1./sslp
C
C Planck function for co2
C
      do i=1,plon
         ex        = exp(960./tplnke(i))
         co2plk(i) = 5.e8/((tplnke(i)**4)*(ex - 1.))
         co2t(i,1) = tplnke(i)
         sum(i)    = co2t(i,1)*pnm(i,1)
      end do
      k = 1
      do k1=plevp,2,-1
         k = k + 1
         do i=1,plon
            sum(i)         = sum(i) + tlayr(i,k)*(pnm(i,k)-pnm(i,k-1))
            ex             = exp(960./tlayr(i,k1))
            tlayr5         = tlayr(i,k1)*tlayr4(i,k1)
            co2eml(i,k1-1) = 1.2e11*ex/(tlayr5*(ex - 1.)**2)
            co2t(i,k)      = sum(i)/pnm(i,k)
         end do
      end do
      bndfct = 2.0*22.18/(sqrt(196.)*300.)
C
C Initialize planck function derivative for O3
C
      do i=1,plon
         dbvtt(i) = dbvt(tplnke(i))
      end do
c
c   Calculate trace gas Planck functions
c
      call trcplk(tint, tlayr, tplnke, emplnk, abplnk1, abplnk2)
C
C Interface loop
C
      do 200 k1=1,plevp
C
C H2O emissivity
C
C emis(i,1)     0 -  800 cm-1   rotation band
C emis(i,2)  1200 - 2200 cm-1   vibration-rotation band
C emis(i,3)   800 - 1200 cm-1   window
C emis(i,4)   500 -  800 cm-1   rotation band overlap with co2
C
C For the p type continuum
C
         do i=1,plon
            uc(i)     = s2c(i,k1) + 2.e-3*plh2o(i,k1)
            u(i)      = plh2o(i,k1)
            pnew(i)   = u(i)/w(i,k1)
C
C Apply scaling factor for 500-800 continuum
C
            uc1(i)    = (s2c(i,k1) + 1.7e-3*plh2o(i,k1))*
     $                 (1. + 2.*s2c(i,k1))/(1. + 15.*s2c(i,k1))
            tpathe(i) = s2t(i,k1)/plh2o(i,k1)
         end do
         do i=1,plon
            dtx(i) = tplnke(i) - 250.
            dty(i) = tpathe(i) - 250.
            dtz(i) = dtx(i) - 50.
            dtp(i) = dty(i) - 50.
         end do
         do iband=1,3,2
            do i=1,plon
               term1(i,iband) = coefe(1,iband) + coefe(2,iband)*
     $                          dtx(i)*(1. + c1(iband)*dtx(i))
               term2(i,iband) = coefb(1,iband) + coefb(2,iband)*
     $                          dtx(i)*(1. + c2(iband)*dtx(i)*
     $                                  (1. + c3(iband)*dtx(i)))
               term3(i,iband) = coefd(1,iband) + coefd(2,iband)*
     $                          dtx(i)*(1. +  c4(iband)*dtx(i)*
     $                                  (1. + c5(iband)*dtx(i)))
               term4(i,iband) = coefa(1,iband) + coefa(2,iband)*
     $                          dty(i)*(1. + c6(iband)*dty(i))
               term5(i,iband) = coefc(1,iband) + coefc(2,iband)*
     $                          dty(i)*(1. + c7(iband)*dty(i))
            end do
         end do
C
C emis(i,1)     0 -  800 cm-1   rotation band
C
         do i=1,plon
            a11  = .37 - 3.33e-5*dtz(i) + 3.33e-6*dtz(i)*dtz(i)
            a31  = 1.07 - 1.00e-3*dtp(i) + 1.475e-5*dtp(i)*dtp(i)
            a21  = 1.3870 + 3.80e-3*dtz(i) - 7.8e-6*dtz(i)*dtz(i)
            a22  = 1.0 - 1.21e-3*dtp(i) - 5.33e-6*dtp(i)*dtp(i)
            a23  = 0.9 + 2.62*sqrt(u(i))
            corfac(i) = a31*(a11 + ((a21*a22)/a23))
            t1t4 = term1(i,1)*term4(i,1)
            t2t5 = term2(i,1)*term5(i,1)
            a(i) = t1t4 + t2t5/(1. + t2t5*sqrt(u(i))*corfac(i))
            fwk  = fwcoef + fwc1/(1. + fwc2*u(i))
            rsum(i)   = exp(-a(i)*(sqrt(u(i)) + fwk*u(i)))
            emis(i,1) = (1. - rsum(i))*term3(i,1)
C            trem1(i)  = rsum(i)
C
C emis(i,2)  1200 - 2200 cm-1   vibration-rotation band
C
            a41      = 1.75 - 3.96e-3*dtz(i)
            a51      = 1.00 + 1.3*sqrt(u(i))
            a61      = 1.00 + 1.25e-3*dtp(i) + 6.25e-5*dtp(i)*dtp(i)
            corfac(i)= .3*(1. + (a41)/(a51))*a61
            t1t4     = term1(i,3)*term4(i,3)
            t2t5     = term2(i,3)*term5(i,3)
            a(i)     = t1t4 + t2t5/(1. + t2t5*sqrt(u(i))*corfac(i))
            fwk      = fwcoef + fwc1/(1. + fwc2*u(i))
            rsum(i)  = exp(-a(i)*(sqrt(u(i)) + fwk*u(i)))
            emis(i,2)= (1. - rsum(i))*term3(i,3)
C            trem7(i) = rsum(i)
         end do
C
C Line transmission in 800-1000 and 1000-1200 cm-1 intervals
C
         do k=1,2
            do i=1,plon
               phi  = a1(k)*(dty(i) + 15.) + a2(k)*(dty(i) + 15.)**2
               psi  = b1(k)*(dty(i) + 15.) + b2(k)*(dty(i) + 15.)**2
               phi  = exp(phi)
               psi  = exp(psi)
               ubar = w(i,k1)*phi
               ubar = (ubar*1.66)*r80257
               pbar = pnew(i)*(psi/phi)
               cf812 = cfa1 + ((1.-cfa1)/(1. + ubar*pbar*10.))
               g1   = (realk(k)*pbar)/(2.*st(k))
               g2   = 1. + (ubar*4.0*st(k)*cf812)/pbar
               g3   = sqrt(g2) - 1.
               g4   = g1*g3
               trline(i,k) = exp(-g4)
            end do
         end do
         do i=1,plon
            term7(i,1) = coefj(1,1) + coefj(2,1)*dty(i)*(1.+c16*dty(i))
            term8(i,1) = coefk(1,1) + coefk(2,1)*dty(i)*(1.+c17*dty(i))
            term7(i,2) = coefj(1,2) + coefj(2,2)*dty(i)*(1.+c26*dty(i))
            term8(i,2) = coefk(1,2) + coefk(2,2)*dty(i)*(1.+c27*dty(i))
         end do
C
C emis(i,3)   800 - 1200 cm-1   window
C
         do i=1,plon
            term6(i,1) = coeff(1,1) + coeff(2,1)*dtx(i)*
     $                  (1. +  c8*dtx(i)*(1. + c10*dtx(i)*
     $                  (1. + c12*dtx(i)*(1. + c14*dtx(i)))))
            trem4(i)  = exp(-(coefg(1,1)+coefg(2,1)*dtx(i))*uc(i))
     $                  *trline(i,2)
            trem6(i)  = exp(-(coefg(1,2)+coefg(2,2)*dtx(i))*uc(i))
     $                  *trline(i,1)
            emis(i,3) = term6(i,1)*(1. - .5*trem4(i) -.5*trem6(i))
C
C emis(i,4)   500 -  800 cm-1   rotation band overlap with co2
C
            k21(i) = term7(i,1) + term8(i,1)/
     $           (1. + (c30 + c31*(dty(i)-10.)*(dty(i)-10.))*sqrt(u(i)))
            k22(i) = term7(i,2) + term8(i,2)/
     $           (1. + (c28 + c29*(dty(i)-10.))*sqrt(u(i)))
            term9(i,1) = coefi(1,1) + coefi(2,1)*dtx(i)*
     $                  (1. + c18*dtx(i)*(1. + c20*dtx(i)*
     $                   (1. + c22*dtx(i)*(1. + c24*dtx(i)))))
            fwk    = fwcoef + fwc1/(1.+fwc2*u(i))
            tr1(i) = exp(-(k21(i)*(sqrt(u(i)) + fc1*fwk*u(i))))
            tr2(i) = exp(-(k22(i)*(sqrt(u(i)) + fc1*fwk*u(i))))
            tr3(i) = exp(-((coefh(1,1) + coefh(2,1)*dtx(i))*uc1(i)))
            tr4(i) = exp(-((coefh(1,2) + coefh(2,2)*dtx(i))*uc1(i)))
            tr7(i) = tr1(i)*tr3(i)
            tr8(i) = tr2(i)*tr4(i)
            emis(i,4) = term9(i,1)*.5*(tr1(i)-tr7(i) + tr2(i)-tr8(i))
            h2oems(i,k1) = emis(i,1)+emis(i,2)+emis(i,3)+emis(i,4)
            troco2(i,k1) = 0.65*tr7(i) + 0.35*tr8(i)
            th2o(i) = tr8(i)
C            trem2(i)     = troco2(i,k1)
         end do
C
C CO2 emissivity for 15 micron band system
C
         do 100 i=1,plon
            t1i    = exp(-480./co2t(i,k1))
            sqti   = sqrt(co2t(i,k1))
            rsqti  = 1./sqti
            et     = t1i
            et2    = et*et
            et4    = et2*et2
            omet   = 1. - 1.5*et2
            f1co2  = 899.70*omet*(1. + 1.94774*et + 4.73486*et2)*rsqti
            sqwp   = sqrt(plco2(i,k1))
            f1sqwp = f1co2*sqwp
            t1co2  = 1./(1. + 245.18*omet*sqwp*rsqti)
            oneme  = 1. - et2
            alphat = oneme**3*rsqti
            wco2   = 2.5221*co2vmr*pnm(i,k1)*rga
            u7     = 4.9411e4*alphat*et2*wco2
            u8     = 3.9744e4*alphat*et4*wco2
            u9     = 1.0447e5*alphat*et4*et2*wco2
            u13    = 2.8388e3*alphat*et4*wco2
C
            tpath  = co2t(i,k1)
            tlocal = tplnke(i)
            tcrfac = sqrt((tlocal*r250)*(tpath*r300))
            pi     = pnm(i,k1)*rsslp + 2.*dpfco2*tcrfac
            posqt  = pi/(2.*sqti)
            rbeta7 =  1./( 5.3288*posqt)
            rbeta8 = 1./ (10.6576*posqt)
            rbeta9 = rbeta7
            rbeta13= rbeta9
            f2co2  = (u7/sqrt(4. + u7*(1. + rbeta7))) +
     $               (u8/sqrt(4. + u8*(1. + rbeta8))) +
     $               (u9/sqrt(4. + u9*(1. + rbeta9)))
            f3co2  = u13/sqrt(4. + u13*(1. + rbeta13))
            tmp1   = alog(1. + f1sqwp)
            tmp2   = alog(1. +  f2co2)
            tmp3   = alog(1. +  f3co2)
            absbnd = (tmp1 + 2.*t1co2*tmp2 + 2.*tmp3)*sqti
            tco2(i)=1.0/(1.0+10.0*(u7/sqrt(4. + u7*(1. + rbeta7))))
            co2ems(i,k1)  = troco2(i,k1)*absbnd*co2plk(i)
            ex     = exp(960./tint(i,k1))
            exm1sq = (ex - 1.)**2
            co2em(i,k1) = 1.2e11*ex/(tint(i,k1)*tint4(i,k1)*exm1sq)
C            trem3(i) = 1. - bndfct*absbnd
  100    continue
C
C O3 emissivity
C
         do i=1,plon
            h2otr(i,k1) = exp(-12.*s2c(i,k1))
            te          = (co2t(i,k1)/293.)**.7
            u1          = 18.29*plos(i,k1)/te
            u2          = .5649*plos(i,k1)/te
            phat        = plos(i,k1)/plol(i,k1)
            tlocal      = tplnke(i)
            tcrfac      = sqrt(tlocal*r250)*te
            beta        = (1./.3205)*((1./phat) + (dpfo3*tcrfac))
            realnu      = (1./beta)*te
            o3bndi      = 74.*te*(tplnke(i)/375.)*
     $         alog(1. + fo3(u1,realnu) + fo3(u2,realnu))
            o3ems(i,k1) = dbvtt(i)*h2otr(i,k1)*o3bndi
            to3(i)=1.0/(1. + 0.1*fo3(u1,realnu) + 0.1*fo3(u2,realnu))
C            trem5(i)    = 1.-(o3bndi/(1060-980.))
         end do
c
c   Calculate trace gas emissivities
c
      call trcems(k1,     co2t,   pnm,    ucfc11, ucfc12, un2o0,  
     $            un2o1,  bn2o0,  bn2o1,  uch4,   bch4,   uco211,
     $            uco212, uco213, uco221, uco222, uco223, uptype,
     $            w,      s2c,    u,      emplnk, th2o,   tco2,   
     $            to3,    emstrc)
C
C Total emissivity:
C
         do i=1,plon
            emstot(i,k1) = h2oems(i,k1) + co2ems(i,k1) + o3ems(i,k1)
     $                     + emstrc(i,k1)
         end do
  200 continue          ! End of interface loop
C
      return
      end
      subroutine radini(gravx   ,cpairx  ,epsilox ,stebolx )
C-----------------------------------------------------------------------
C
C Initialize various constants for radiation scheme; note that
C the radiation scheme uses cgs units.
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Kiehl, B. Briegleb, August 1992
C
C-----------------------------------------------------------------------
c
c $Id: radini.F,v 1.2 1995/02/17 21:28:44 jhack Exp $
c $Author: jhack $
c
      implicit none
C------------------------------Parameters-------------------------------
c
c $Id: prgrid.h,v 1.2 1995/02/10 01:09:10 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation resolution and I/O parameters
C
      integer
     $     plon,    ! number of longitudes
     $     plev,    ! number of vertical levels
     $     plat,    ! number of latitudes
     $     pcnst,   ! number of constituents (including water vapor)
     $     plevmx,  ! number of subsurface levels
     $     plevp,   ! plev + 1
     $     nxpt,    ! no.of points outside active domain for interpolant
     $     jintmx,  ! number of extra latitudes in polar region
     $     plond,   ! slt extended domain longitude
     $     platd,   ! slt extended domain lat.
     $     plevd    ! fold plev,pcnst indices into one
      integer
     $     plngbuf  ! length of absorptivity/emissivity record
C
      parameter(plon    = 1,
     $          plev    = 18,
     $          plat    = 1,
     $          pcnst   = 1,
     $          plevmx  = 4,
     $          plevp   = plev + 1,
     $          nxpt    = 1,
     $          jintmx  = 1,
     $          plond   = plon + 1 + 2*nxpt,
     $          platd   = plat + 2*nxpt + 2*jintmx,
     $          plevd   = plev*(3 + pcnst))
      parameter(plngbuf = 512*((plond*plevp*plevp + plond*plev*4 +
     $                          plond*plevp)/512 + 1))
C
C------------------------------Commons----------------------------------
c
c $Id: comozp.h,v 1.1.1.1 1995/02/09 23:27:16 ccm2 Exp $
c $Author: ccm2 $
c
C
C ozone mixing ratios, pressures and times
C
      integer pnoz,     ! Maximum number of levels in ozone input data
     $        pozlon    ! Number of ozone longitudes
      parameter (pnoz=100)
      parameter (pozlon=1)
      common /comozp/nyroz   ,ldoyoz  ,ndoyoz  ,cplos   ,cplol   ,
     $               ozmix(plat,pnoz),ozmixm(pozlon,pnoz,plat,2),
     $               pin(pnoz),koz
C
      integer nyroz  ! Year of ozone data just after current date
      real ldoyoz,   ! Day of yr of ozone data prior to current date
     $     ndoyoz,   ! Day of yr of ozone data after current date
     $     cplos,    ! Constant for ozone path length calculation
     $     cplol,    ! Constant for pressure-weighted o3 path length calc.
     $     ozmix,    ! Time-interpolated mixing ratios
     $     ozmixm,   ! Two consecutive time slices of ozone mixing ratios
     $     pin       ! Pressures at model interfaces
      integer koz    ! Actual number of levels in ozone input data
C
C-----------------------------------------------------------------------
c
c $Id: crdcae.h,v 1.1.1.1 1995/02/09 23:26:46 ccm2 Exp $
c $Author: ccm2 $
c
C
C Water vapor narrow band constants for longwave radiation computations
C
      common/crdcae/realk(2), st(2), a1(2), a2(2), b1(2), b2(2),
     $              coefa(3,4),coefb(4,4),coefc(3,4),coefd(4,4),
     $              coefe(3,4),coeff(6,2),coefg(2,4),coefh(2,4),
     $              coefi(6,2),coefj(3,2),coefk(3,2),
     $              c1(4),c2(4),c3(4),c4(4),c5(4),c6(4),c7(4),c8,c9,
     $              c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,
     $              c22,c23,c24,c25,c26,c27,c28,c29,c30,c31,
     $              fwcoef,fwc1,fwc2,fc1,cfa1
C
      real realk,     ! H2O narrow band parameter
     $     st,        ! H2O narrow band parameter
     $     a1,a2,     ! Temperature correction terms for H2O path
     $     b1,b2      ! Temperature correction terms for H2O path
C
C Constant coefficients for water vapor absorptivity and emissivity
C
      real coefa,coefb,coefc,coefd,coefe,coeff,
     $     coefg,coefh,coefi,coefj,coefk,
     $      c1, c2, c3, c4, c5, c6, c7,c8 ,c9 ,c10,
     $     c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,
     $     c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,c31
C
C Farwing correction constants for narrow-band emissivity model,
C introduced to account for the deficiencies in narrow-band model
C used to derive the emissivity; tuned with Arking's line-by-line
C calculations.
C
      real fwcoef,
     $     fwc1,fwc2,
     $     fc1,
     $     cfa1
C
C-----------------------------------------------------------------------
c
c $Id: crdcon.h,v 1.1.1.1 1995/02/09 23:26:46 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation constants
C
      common/crdcon/gravit  ,rga     ,cpair   ,epsilo  ,sslp    ,
     $              stebol  ,rgsslp  ,co2vmr  ,dpfo3   ,dpfco2  ,
     $              dayspy  ,pie
C
      real gravit,    ! Acceleration of gravity
     $     rga,       ! 1./gravit
     $     cpair,     ! Specific heat of dry air
     $     epsilo,    ! Ratio of mol. wght of H2O to dry air
     $     sslp,      ! Standard sea-level pressure
     $     stebol,    ! Stefan-Boltzmann's constant
     $     rgsslp,    ! 0.5/(gravit*sslp)
     $     co2vmr,    ! CO2 volume mixing ratio
     $     dpfo3,     ! Voigt correction factor for O3
     $     dpfco2,    ! Voigt correction factor for CO2
     $     dayspy,    ! Number of days per 1 year
     $     pie        ! 3.14.....
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real gravx,     ! Acceleration of gravity (MKS)
     $     cpairx,    ! Specific heat of dry air (MKS)
     $     epsilox,   ! Ratio of mol. wght of H2O to dry air
     $     stebolx    ! Stefan-Boltzmann's constant (MKS)
C
C---------------------------Local variables-----------------------------
C
      integer iband   ! H2O band index
      real v0,        ! Volume of a gas at stp (m**3/kmol)
     $     p0,        ! Standard pressure (pascals)
     $     amd,       ! Effective molecular weight of dry air (kg/kmol)
     $     goz        ! Acceleration of gravity (m/s**2)
C
C-----------------------------------------------------------------------
C
C Set general radiation consts; convert to cgs units where appropriate:
C
      gravit  =  100.*gravx
      rga     =  1./gravit
      cpair   =  1.e4*cpairx
      epsilo  =  epsilox
      sslp    =  1.013250e6
      stebol  =  1.e3*stebolx
      rgsslp  =  0.5/(gravit*sslp)
      co2vmr  =  3.55e-4
      dpfo3   =  2.5e-3
      dpfco2  =  5.0e-3
      dayspy  =  365.
      pie     =  4.*atan(1.)
C
C Coefficients for h2o emissivity and absorptivity.
C
      do iband=1,4
         c1(iband) = coefe(3,iband)/coefe(2,iband)
         c2(iband) = coefb(3,iband)/coefb(2,iband)
         c3(iband) = coefb(4,iband)/coefb(3,iband)
         c4(iband) = coefd(3,iband)/coefd(2,iband)
         c5(iband) = coefd(4,iband)/coefd(3,iband)
         c6(iband) = coefa(3,iband)/coefa(2,iband)
         c7(iband) = coefc(3,iband)/coefc(2,iband)
      end do
      c8   = coeff(3,1)/coeff(2,1)
      c9   = coeff(3,2)/coeff(2,2)
      c10  = coeff(4,1)/coeff(3,1)
      c11  = coeff(4,2)/coeff(3,2)
      c12  = coeff(5,1)/coeff(4,1)
      c13  = coeff(5,2)/coeff(4,2)
      c14  = coeff(6,1)/coeff(5,1)
      c15  = coeff(6,2)/coeff(5,2)
      c16  = coefj(3,1)/coefj(2,1)
      c17  = coefk(3,1)/coefk(2,1)
      c18  = coefi(3,1)/coefi(2,1)
      c19  = coefi(3,2)/coefi(2,2)
      c20  = coefi(4,1)/coefi(3,1)
      c21  = coefi(4,2)/coefi(3,2)
      c22  = coefi(5,1)/coefi(4,1)
      c23  = coefi(5,2)/coefi(4,2)
      c24  = coefi(6,1)/coefi(5,1)
      c25  = coefi(6,2)/coefi(5,2)
      c26  = coefj(3,2)/coefj(2,2)
      c27  = coefk(3,2)/coefk(2,2)
      c28  = .5
      c29  = .002053
      c30  = .1
      c31  = 3.0e-5
      cfa1 = .61
C
C Initialize further longwave constants referring to far wing
C correction; R&D refers to:
C
C            Ramanathan, V. and  P.Downey, 1986: A Nonisothermal
C            Emissivity and Absorptivity Formulation for Water Vapor
C            Journal of Geophysical Research, vol. 91., D8, pp 8649-8666
C
      fwcoef = .1           ! See eq(33) R&D
      fwc1   = .30          ! See eq(33) R&D
      fwc2   = 4.5          ! See eq(33) and eq(34) in R&D
      fc1    = 2.6          ! See eq(34) R&D
C
C Initialize ozone data.
C
      v0  = 22.4136         ! Volume of a gas at stp (m**3/kmol)
      p0  = 0.1*sslp        ! Standard pressure (pascals)
      amd = 28.9644         ! Molecular weight of dry air (kg/kmol)
      goz = gravx           ! Acceleration of gravity (m/s**2)
C
C Constants for ozone path integrals (multiplication by 100 for unit
C conversion to cgs from mks):
C
      cplos = v0/(amd*goz)       *100.0
      cplol = v0/(amd*goz*p0)*0.5*100.0
C
      return
      end
      subroutine radinp(pmid    ,pint    ,h2ommr  ,cld     ,o3vmr   ,
     $                  pmidrd  ,pintrd  ,plco2   ,plh2o   ,tclrsf  ,
     $                  eccf    ,o3mmr   )
C-----------------------------------------------------------------------
C
C Set latitude and time dependent arrays for input to solar
C and longwave radiation.
C
C Convert model pressures to cgs, compute path length arrays needed for the 
C longwave radiation, and compute ozone mixing ratio, needed for the solar
C radiation.
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Kiehl, B. Briegleb, August 1992
C
C-----------------------------------------------------------------------
c
c $Id: radinp.F,v 1.1.1.1 1995/02/09 23:27:02 ccm2 Exp $
c $Author: ccm2 $
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
C------------------------------Parameters-------------------------------
c
c $Id: prgrid.h,v 1.2 1995/02/10 01:09:10 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation resolution and I/O parameters
C
      integer
     $     plon,    ! number of longitudes
     $     plev,    ! number of vertical levels
     $     plat,    ! number of latitudes
     $     pcnst,   ! number of constituents (including water vapor)
     $     plevmx,  ! number of subsurface levels
     $     plevp,   ! plev + 1
     $     nxpt,    ! no.of points outside active domain for interpolant
     $     jintmx,  ! number of extra latitudes in polar region
     $     plond,   ! slt extended domain longitude
     $     platd,   ! slt extended domain lat.
     $     plevd    ! fold plev,pcnst indices into one
      integer
     $     plngbuf  ! length of absorptivity/emissivity record
C
      parameter(plon    = 1,
     $          plev    = 18,
     $          plat    = 1,
     $          pcnst   = 1,
     $          plevmx  = 4,
     $          plevp   = plev + 1,
     $          nxpt    = 1,
     $          jintmx  = 1,
     $          plond   = plon + 1 + 2*nxpt,
     $          platd   = plat + 2*nxpt + 2*jintmx,
     $          plevd   = plev*(3 + pcnst))
      parameter(plngbuf = 512*((plond*plevp*plevp + plond*plev*4 +
     $                          plond*plevp)/512 + 1))
C
C------------------------------Commons----------------------------------
c
c $Id: comtim.h,v 1.1.1.1 1995/02/09 23:26:44 ccm2 Exp $
c $Author: ccm2 $
c
C
C Model time variables
C
      common/comtim/calday  ,dtime   ,twodt   ,nrstrt  ,nstep   ,
     $              nstepr  ,nestep  ,nelapse ,nstop   ,mdbase  ,
     $              msbase  ,mdcur   ,mscur   ,mbdate  ,mbsec   ,
     $              mcdate  ,mcsec   ,nndbas  ,nnsbas  ,nnbdat  ,
     $              nnbsec  ,doabsems,dosw    ,dolw
C
      real calday,   ! Current calendar day = julian day + fraction
     $     dtime,    ! Time step in seconds (delta t)
     $     twodt     ! 2 * delta t 
      integer
     $     nrstrt,   ! Starting time step of restart run (constant) 
     $     nstep,    ! Current time step
     $     nstepr,   ! Current time step of restart run(updated w/nstep)
     $     nestep,   ! Time step on which to stop run
     $     nelapse,  ! Requested elapsed time for model run
     $     nstop,    ! nestep + 1
     $     mdbase,   ! Base day of run
     $     msbase,   ! Base seconds of base day
     $     mdcur,    ! Current day of run
     $     mscur,    ! Current seconds of current day
     $     mbdate,   ! Base date of run (yymmdd format)
     $     mbsec,    ! Base seconds of base date
     $     mcdate,   ! Current date of run (yymmdd format)
     $     mcsec,    ! Current seconds of current date
     $     nndbas,   ! User input base day
     $     nnsbas,   ! User input base seconds of input base day
     $     nnbdat,   ! User input base date (yymmdd format)
     $     nnbsec    ! User input base seconds of input base date
      logical
     $     doabsems, ! True => abs/emiss calculation this timestep
     $     dosw,     ! True => shortwave calculation this timestep
     $     dolw      ! True => longwave calculation this timestep
C
C-----------------------------------------------------------------------
c
c $Id: crdcon.h,v 1.1.1.1 1995/02/09 23:26:46 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation constants
C
      common/crdcon/gravit  ,rga     ,cpair   ,epsilo  ,sslp    ,
     $              stebol  ,rgsslp  ,co2vmr  ,dpfo3   ,dpfco2  ,
     $              dayspy  ,pie
C
      real gravit,    ! Acceleration of gravity
     $     rga,       ! 1./gravit
     $     cpair,     ! Specific heat of dry air
     $     epsilo,    ! Ratio of mol. wght of H2O to dry air
     $     sslp,      ! Standard sea-level pressure
     $     stebol,    ! Stefan-Boltzmann's constant
     $     rgsslp,    ! 0.5/(gravit*sslp)
     $     co2vmr,    ! CO2 volume mixing ratio
     $     dpfo3,     ! Voigt correction factor for O3
     $     dpfco2,    ! Voigt correction factor for CO2
     $     dayspy,    ! Number of days per 1 year
     $     pie        ! 3.14.....
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real pmid(plond,plev),   ! Pressure at model mid-levels (pascals)
     $     pint(plond,plevp),  ! Pressure at model interfaces (pascals)
     $     h2ommr(plond,plev), ! H2o mass mixing ratio
     $     cld(plond,plevp),   ! Fractional cloud cover
     $     o3vmr(plond,plev)   ! ozone volume mixing ratio
C
C Output arguments
C
      real pmidrd(plond,plev), ! Pressure at mid-levels (dynes/cm*2)
     $     pintrd(plond,plevp),! Pressure at interfaces (dynes/cm*2)
     $     plco2(plond,plevp), ! Vert. pth lngth of co2 (prs-weighted)
     $     plh2o(plond,plevp), ! Vert. pth lngth h2o vap.(prs-weighted)
     $     tclrsf(plond,plevp) ! Product of clr-sky fractions from top
C                              ! of atmosphere to level.
      real eccf,               ! Earth-sun distance factor
     $     o3mmr(plond,plev)   ! Ozone mass mixing ratio
C
C---------------------------Local variables-----------------------------
C
      integer i,    ! Longitude loop index
     $        k     ! Vertical loop index
      real theta,   ! Earth orbit seasonal angle in radians
     $     p0 ,     ! Standard pressure (dynes/cm**2)
     $     amd,     ! Effective molecular weight of dry air (g/mol)
     $     amo,     ! Molecular weight of ozone (g/mol)
     $     amco2,   ! Molecular weight of co2   (g/mol)
     $     cpwpl,   ! Const in co2 mixing ratio to path length conversn
     $     vmmr     ! Ozone volume mixing ratio
      save     p0   ,amd   ,amo  ,amco2
C
      data p0    /  1.01325e6 /
      data amd   /  28.9644   /
      data amo   /  48.0000   /
      data amco2 /  44.0000   /
C
C-----------------------------------------------------------------------
C
C Compute solar distance factor and cosine solar zenith angle usi
C day value where a round day (such as 213.0) refers to 0z at
C Greenwich longitude.
C
C Use formulas from Paltridge, G.W. and C.M.R. Platt 1976: Radiative
C Processes in Meterology and Climatology, Elsevier Scientific
C Publishing Company, New York  p. 57, p. 62,63.
C
C Compute eccentricity factor (sun-earth distance factor)
C
      theta = 2.*pie*calday/dayspy
      eccf  = 1.000110 + .034221*cos(theta) + .001280*sin(theta) +
     $         .000719*cos(2.*theta) + .000077*sin(2.*theta)
C
C Convert pressure from pascals to dynes/cm2
C
      do k=1,plev
         do i=1,plon
            pmidrd(i,k) = pmid(i,k)*10.0
            pintrd(i,k) = pint(i,k)*10.0
         end do
      end do
      do i=1,plon
         pintrd(i,plevp) = pint(i,plevp)*10.0
      end do
C
C Compute path quantities used in the longwave radiation:
C
      vmmr  = amco2/amd
      cpwpl = vmmr*0.5/(gravit*p0)
      do i=1,plon
         plh2o(i,1)  = rgsslp*h2ommr(i,1)*pintrd(i,1)*pintrd(i,1)
         plco2(i,1)  = co2vmr*cpwpl*pintrd(i,1)*pintrd(i,1)
         tclrsf(i,1) = 1.
      end do
      do k=1,plev
         do i=1,plon
            plh2o(i,k+1)  = plh2o(i,k) + rgsslp*
     $             (pintrd(i,k+1)**2 - pintrd(i,k)**2)*h2ommr(i,k)
            plco2(i,k+1)  = co2vmr*cpwpl*pintrd(i,k+1)**2
            tclrsf(i,k+1) = tclrsf(i,k)*(1.-cld(i,k+1))
         end do
      end do
C
C Convert ozone volume mixing ratio to mass mixing ratio:
C
      vmmr = amo/amd
      do k=1,plev
         do i=1,plon
            o3mmr(i,k) = vmmr*o3vmr(i,k)
         end do
      end do
C
      return
      end





      subroutine radoz2(o3vmr   ,pint    ,plol    ,plos    )
C-----------------------------------------------------------------------
C
C Computes the path length integrals to the model interfaces given the 
C ozone volume mixing ratio
C
C---------------------------Code history--------------------------------
C
C Original version:     CCM1
C Standardized:         J. Rosinski, June 1992
C Reviewed:             J. Kiehl, B. Briegleb, August 1992
C Mixing ratio version: Bruce Biegleb, September 1992
C
C-----------------------------------------------------------------------
c
c $Id: radoz2.F,v 1.1.1.1 1995/02/09 23:27:03 ccm2 Exp $
c $Author: ccm2 $
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
C-----------------------------------------------------------------------
c
c $Id: prgrid.h,v 1.2 1995/02/10 01:09:10 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation resolution and I/O parameters
C
      integer
     $     plon,    ! number of longitudes
     $     plev,    ! number of vertical levels
     $     plat,    ! number of latitudes
     $     pcnst,   ! number of constituents (including water vapor)
     $     plevmx,  ! number of subsurface levels
     $     plevp,   ! plev + 1
     $     nxpt,    ! no.of points outside active domain for interpolant
     $     jintmx,  ! number of extra latitudes in polar region
     $     plond,   ! slt extended domain longitude
     $     platd,   ! slt extended domain lat.
     $     plevd    ! fold plev,pcnst indices into one
      integer
     $     plngbuf  ! length of absorptivity/emissivity record
C
      parameter(plon    = 1,
     $          plev    = 18,
     $          plat    = 1,
     $          pcnst   = 1,
     $          plevmx  = 4,
     $          plevp   = plev + 1,
     $          nxpt    = 1,
     $          jintmx  = 1,
     $          plond   = plon + 1 + 2*nxpt,
     $          platd   = plat + 2*nxpt + 2*jintmx,
     $          plevd   = plev*(3 + pcnst))
      parameter(plngbuf = 512*((plond*plevp*plevp + plond*plev*4 +
     $                          plond*plevp)/512 + 1))
C
C------------------------------Commons----------------------------------
c
c $Id: comozp.h,v 1.1.1.1 1995/02/09 23:27:16 ccm2 Exp $
c $Author: ccm2 $
c
C
C ozone mixing ratios, pressures and times
C
      integer pnoz,     ! Maximum number of levels in ozone input data
     $        pozlon    ! Number of ozone longitudes
      parameter (pnoz=100)
      parameter (pozlon=1)
      common /comozp/nyroz   ,ldoyoz  ,ndoyoz  ,cplos   ,cplol   ,
     $               ozmix(plat,pnoz),ozmixm(pozlon,pnoz,plat,2),
     $               pin(pnoz),koz
C
      integer nyroz  ! Year of ozone data just after current date
      real ldoyoz,   ! Day of yr of ozone data prior to current date
     $     ndoyoz,   ! Day of yr of ozone data after current date
     $     cplos,    ! Constant for ozone path length calculation
     $     cplol,    ! Constant for pressure-weighted o3 path length calc.
     $     ozmix,    ! Time-interpolated mixing ratios
     $     ozmixm,   ! Two consecutive time slices of ozone mixing ratios
     $     pin       ! Pressures at model interfaces
      integer koz    ! Actual number of levels in ozone input data
C
C------------------------------Input arguments--------------------------
C
      real o3vmr(plond,plev)   ! ozone volume mixing ratio
      real pint(plond,plevp)   ! Model interface pressures
C
C----------------------------Output arguments---------------------------
C
      real plol(plond,plevp),  ! Ozone prs weighted path length (cm)
     $     plos(plond,plevp)   ! Ozone path length (cm)
C
C---------------------------Local workspace-----------------------------
C
      integer   i,             ! longitude index
     $          k              ! level index
C
C-----------------------------------------------------------------------
C
C Evaluate the ozone path length integrals to interfaces; 
C factors of .1 and .01 to convert pressures from cgs to mks:
C
C Bug fix, 24 May 1996:  the 0.5 and 0.25 factors removed.
C
      do i=1,plon
         plos(i,1) = 0.1*cplos*o3vmr(i,1)*pint(i,1)
         plol(i,1) = 0.01*cplol*o3vmr(i,1)*pint(i,1)*pint(i,1)
      end do
      do k=2,plevp
         do i=1,plon
            plos(i,k) = plos(i,k-1) + 0.1*cplos*o3vmr(i,k-1)*
     $                  (pint(i,k) - pint(i,k-1))
            plol(i,k) = plol(i,k-1) + 0.01*cplol*o3vmr(i,k-1)*
     $                  (pint(i,k)*pint(i,k) - pint(i,k-1)*pint(i,k-1))
         end do
      end do
C
      return
      end
      subroutine radtpl(tnm    ,ts     ,qnm    ,pnm    ,plh2o  ,
     $                  tplnka ,s2c    ,s2t    ,w      ,tplnke ,
     $                  tint   ,tint4  ,tlayr  ,tlayr4 ,pmln   ,
     $                  piln   )
C-----------------------------------------------------------------------
C
C Compute temperatures and path lengths for longwave radiation
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      L. Buja, June 1992
C Reviewed:          J. Kiehl, B. Briegleb, August 1992
C
C-----------------------------------------------------------------------
c
c $Id: radtpl.F,v 1.1.1.1 1995/02/09 23:27:03 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
C------------------------------Parameters-------------------------------
c
c $Id: prgrid.h,v 1.2 1995/02/10 01:09:10 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation resolution and I/O parameters
C
      integer
     $     plon,    ! number of longitudes
     $     plev,    ! number of vertical levels
     $     plat,    ! number of latitudes
     $     pcnst,   ! number of constituents (including water vapor)
     $     plevmx,  ! number of subsurface levels
     $     plevp,   ! plev + 1
     $     nxpt,    ! no.of points outside active domain for interpolant
     $     jintmx,  ! number of extra latitudes in polar region
     $     plond,   ! slt extended domain longitude
     $     platd,   ! slt extended domain lat.
     $     plevd    ! fold plev,pcnst indices into one
      integer
     $     plngbuf  ! length of absorptivity/emissivity record
C
      parameter(plon    = 1,
     $          plev    = 18,
     $          plat    = 1,
     $          pcnst   = 1,
     $          plevmx  = 4,
     $          plevp   = plev + 1,
     $          nxpt    = 1,
     $          jintmx  = 1,
     $          plond   = plon + 1 + 2*nxpt,
     $          platd   = plat + 2*nxpt + 2*jintmx,
     $          plevd   = plev*(3 + pcnst))
      parameter(plngbuf = 512*((plond*plevp*plevp + plond*plev*4 +
     $                          plond*plevp)/512 + 1))
C
C------------------------------Commons----------------------------------
c
c $Id: crdcon.h,v 1.1.1.1 1995/02/09 23:26:46 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation constants
C
      common/crdcon/gravit  ,rga     ,cpair   ,epsilo  ,sslp    ,
     $              stebol  ,rgsslp  ,co2vmr  ,dpfo3   ,dpfco2  ,
     $              dayspy  ,pie
C
      real gravit,    ! Acceleration of gravity
     $     rga,       ! 1./gravit
     $     cpair,     ! Specific heat of dry air
     $     epsilo,    ! Ratio of mol. wght of H2O to dry air
     $     sslp,      ! Standard sea-level pressure
     $     stebol,    ! Stefan-Boltzmann's constant
     $     rgsslp,    ! 0.5/(gravit*sslp)
     $     co2vmr,    ! CO2 volume mixing ratio
     $     dpfo3,     ! Voigt correction factor for O3
     $     dpfco2,    ! Voigt correction factor for CO2
     $     dayspy,    ! Number of days per 1 year
     $     pie        ! 3.14.....
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real tnm(plond,plev),     ! Model level temperatures
     $     ts(plond),           ! Surface skin temperature
     $     qnm(plond,plev),     ! Model level specific humidity
     $     pnm(plond,plevp),    ! Pressure at model interfaces (dynes/cm2)
     $     plh2o(plond,plevp)   ! Pressure weighted h2o path
C
C Output arguments
C
      real tplnka(plond,plevp), ! Level temperature from interface temperatures
     $     s2c(plond,plevp),    ! H2o continuum path length
     $     s2t(plond,plevp),    ! H2o tmp and prs wghtd path length
     $     w(plond,plevp),      ! H2o path length
     $     tplnke(plond),       ! Equal to tplnka
     $     tint(plond,plevp),   ! Layer interface temperature
     $     tint4(plond,plevp),  ! Tint to the 4th power
     $     tlayr(plond,plevp),  ! K-1 level temperature
     $     tlayr4(plond,plevp), ! Tlayr to the 4th power
     $     pmln(plond,plev),    ! Ln(pmidm1)
     $     piln(plond,plevp)    ! Ln(pintm1)
C
C---------------------------Local variables-----------------------------
C
      integer   i,              ! Longitude index
     $          k               ! Level index
      real r296,                ! Inverse stand temp for h2o continuum
     $     repsil,              ! Inver ratio mol weight h2o to dry air
     $         dy,              ! Thickness of layer for tmp interp
     $       dpnm,              ! Pressure thickness of layer
     $     dpnmsq,              ! Prs squared difference across layer
     $       rtnm               ! Inverse level temperature
C
C-----------------------------------------------------------------------
C
      r296   = 1./296.
      repsil = 1./epsilo
C
C Set the top and bottom intermediate level temperatures,
C top level planck temperature and top layer temp**4.
C
C Tint is lower interface temperature
C (not available for bottom layer, so use ground temperature)
C
      do i=1,plon
         tint(i,plevp)  = ts(i)
         tint4(i,plevp) = tint(i,plevp)**4
         tplnka(i,1)    = tnm(i,1)
         tint(i,1)      = tplnka(i,1)
         tlayr4(i,1)    = tplnka(i,1)**4
         tint4(i,1)     = tlayr4(i,1)
      end do
C
C Intermediate level temperatures are computed using temperature
C at the full level below less dy*delta t,between the full level
C
      do k=2,plev
         do i=1,plon
            dy = (piln(i,k) - pmln(i,k))/(pmln(i,k-1) - pmln(i,k))
            tint(i,k)  = tnm(i,k) - dy*(tnm(i,k)-tnm(i,k-1))
            tint4(i,k) = tint(i,k)**4
         end do
      end do
C
C Now set the layer temp=full level temperatures and establish a
C planck temperature for absorption (tplnka) which is the average
C the intermediate level temperatures.  Note that tplnka is not
C equal to the full level temperatures.
C
      do k=2,plevp
         do i=1,plon
            tlayr(i,k)  = tnm(i,k-1)
            tlayr4(i,k) = tlayr(i,k)**4
            tplnka(i,k) = .5*(tint(i,k) + tint(i,k-1))
         end do
      end do
C
C Calculate tplank for emissivity calculation.
C Assume isothermal tplnke i.e. all levels=ttop.
C
      do i=1,plon
         tplnke(i)  = tplnka(i,1)
         tlayr(i,1) = tint(i,1)
      end do
C
C Now compute h2o path fields:
C
      do i=1,plon
         s2t(i,1) = plh2o(i,1) * tnm(i,1)
         w(i,1)   = (plh2o(i,1)*2.) / pnm(i,1)
         s2c(i,1) = plh2o(i,1) * qnm(i,1) * repsil
      end do
      do k=1,plev
         do i=1,plon
            dpnm       = pnm(i,k+1) - pnm(i,k)
            dpnmsq     = pnm(i,k+1)**2 - pnm(i,k)**2
            rtnm       = 1./tnm(i,k)
            s2t(i,k+1) = s2t(i,k) + rgsslp*dpnmsq*qnm(i,k)*tnm(i,k)
            w(i,k+1)   = w(i,k)   + rga*qnm(i,k)*dpnm
            s2c(i,k+1) = s2c(i,k) + rgsslp*dpnmsq*qnm(i,k)*
     $             exp(1800.*(rtnm - r296))*qnm(i,k)*repsil
         end do
      end do
C 
      return
      end
      subroutine resetr(pa      ,kdim    ,pvalue  )
C-----------------------------------------------------------------------
C
C Reset array pa(kdim) to pvalue
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      L. Bath, June 1992
C
C-----------------------------------------------------------------------
c
c $Id: resetr.F,v 1.1.1.1 1995/02/09 23:27:04 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer kdim      ! Dimension of array pa
      real pvalue       ! Value to store in pa
C
C Output arguments
C
      real pa(kdim)     ! Array to reset
C
C---------------------------Local variable------------------------------
C
      integer j         ! Loop index
C
C-----------------------------------------------------------------------
C
      do j=1,kdim
         pa(j) = pvalue
      end do
C
      return
      end
      subroutine torgrid(pmidm   ,pintm   ,pmlnm   ,pilnm   ,tm      ,
     $                   h2ommrm ,cldm    ,effcldm ,clwpm   ,
     $                   pmid    ,pint    ,pmln    ,piln    ,t       ,
     $                   h2ommr  ,cld     ,effcld  ,clwp    )
C-----------------------------------------------------------------------
C
C Interpolate model arrays to radiation vertical grid.
C
C------------------------------Parameters-------------------------------
c
c $Id: torgrid.F,v 1.1.1.1 1995/02/09 23:27:10 ccm2 Exp $
c $Author: ccm2 $
c
c
c $Id: pmgrid.h,v 1.2 1995/02/10 01:09:06 ccm2 Exp $
c $Author: ccm2 $
c
C
C Basic grid point resolution parameters
C
      integer
     $     plon,    ! number of longitudes
     $     plev,    ! number of vertical levels
     $     plat,    ! number of latitudes
     $     pcnst,   ! number of constituents (including water vapor)
     $     plevmx,  ! number of subsurface levels
     $     plevp,   ! plev + 1
     $     nxpt,    ! no.of points outside active domain for interpolant
     $     jintmx,  ! number of extra latitudes in polar region
     $     plond,   ! slt extended domain longitude
     $     platd,   ! slt extended domain lat.
     $     plevd    ! fold plev,pcnst indices into one
C
      parameter(plon   = 1,
     $          plev   = 18,
     $          plat   = 1,
     $          pcnst  = 1,
     $          plevmx = 4,
     $          plevp  = plev + 1,
     $          nxpt   = 1,
     $          jintmx = 1,
     $          plond  = plon + 1 + 2*nxpt,
     $          platd  = plat + 2*nxpt + 2*jintmx,
     $          plevd  = plev*(3 + pcnst))
C
c
c $Id: ptrrgrid.h,v 1.1.1.1 1995/02/09 23:26:59 ccm2 Exp $
c $Author: ccm2 $
c
C
C Define radiation vertical grid and buffer length for abs/ems out-of-core file
C
      integer
     $     plevr,   ! number of vertical levels
     $     plevrp,  ! plevr + 1
     $     plngbuf  ! length of absorptivity/emissivity record
C
      parameter(plevr = 18,
     $          plevrp = plevr + 1,
     $          plngbuf = 512*((plond*plevrp*plevrp + plond*plevr*4 +
     $                          plond*plevrp)/512 + 1))
C
C-----------------------------------------------------------------------
C
C Input arguments
C
      real pmidm(plond,plev),
     $     pintm(plond,plevp),
     $     pmlnm(plond,plev),
     $     pilnm(plond,plevp),
     $     tm(plond,plev),
     $     h2ommrm(plond,plev),
     $     cldm(plond,plevp),
     $     effcldm(plond,plevp),
     $     clwpm(plond,plev)
C
C Output arguments
C
      real pmid(plond,plevr),
     $     pint(plond,plevrp),
     $     pmln(plond,plevr),
     $     piln(plond,plevrp),
     $     t(plond,plevr),
     $     h2ommr(plond,plevr),
     $     cld(plond,plevrp),
     $     effcld(plond,plevrp),
     $     clwp(plond,plevr)
C
C Code to interpolate goes here.  Do nothing could be coded as a memory
C transfer, but left out for efficiency considerations.
C
      return
      end
      subroutine trcab(k1, k2, ucfc11, ucfc12, un2o0,  un2o1,
     $                         uch4,   uco211, uco212, uco213, 
     $                         uco221, uco222, uco223, bn2o0, 
     $                         bn2o1,  bch4,   to3co2, pnm,
     $                         dw,     pnew,   s2c,    uptype,
     $                         dplh2o, abplnk1,tco2,   th2o,   
     $                         to3,    abstrc)
c----------------------------------------------------------------------
c Calculate absorptivity for non nearest layers for CH4, N2O, CFC11 and
c CFC12.
c
c             Coded by J.T. Kiehl November 21, 1994
c-----------------------------------------------------------------------
c
c $Id: trcab.F,v 1.2 1995/02/17 21:28:52 jhack Exp $
c $Author: jhack $
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
C------------------------------Parameters-------------------------------
c
c $Id: prgrid.h,v 1.2 1995/02/10 01:09:10 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation resolution and I/O parameters
C
      integer
     $     plon,    ! number of longitudes
     $     plev,    ! number of vertical levels
     $     plat,    ! number of latitudes
     $     pcnst,   ! number of constituents (including water vapor)
     $     plevmx,  ! number of subsurface levels
     $     plevp,   ! plev + 1
     $     nxpt,    ! no.of points outside active domain for interpolant
     $     jintmx,  ! number of extra latitudes in polar region
     $     plond,   ! slt extended domain longitude
     $     platd,   ! slt extended domain lat.
     $     plevd    ! fold plev,pcnst indices into one
      integer
     $     plngbuf  ! length of absorptivity/emissivity record
C
      parameter(plon    = 1,
     $          plev    = 18,
     $          plat    = 1,
     $          pcnst   = 1,
     $          plevmx  = 4,
     $          plevp   = plev + 1,
     $          nxpt    = 1,
     $          jintmx  = 1,
     $          plond   = plon + 1 + 2*nxpt,
     $          platd   = plat + 2*nxpt + 2*jintmx,
     $          plevd   = plev*(3 + pcnst))
      parameter(plngbuf = 512*((plond*plevp*plevp + plond*plev*4 +
     $                          plond*plevp)/512 + 1))
C
C------------------------------Commons----------------------------------
c
c $Id: crdcon.h,v 1.1.1.1 1995/02/09 23:26:46 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation constants
C
      common/crdcon/gravit  ,rga     ,cpair   ,epsilo  ,sslp    ,
     $              stebol  ,rgsslp  ,co2vmr  ,dpfo3   ,dpfco2  ,
     $              dayspy  ,pie
C
      real gravit,    ! Acceleration of gravity
     $     rga,       ! 1./gravit
     $     cpair,     ! Specific heat of dry air
     $     epsilo,    ! Ratio of mol. wght of H2O to dry air
     $     sslp,      ! Standard sea-level pressure
     $     stebol,    ! Stefan-Boltzmann's constant
     $     rgsslp,    ! 0.5/(gravit*sslp)
     $     co2vmr,    ! CO2 volume mixing ratio
     $     dpfo3,     ! Voigt correction factor for O3
     $     dpfco2,    ! Voigt correction factor for CO2
     $     dayspy,    ! Number of days per 1 year
     $     pie        ! 3.14.....
C
C------------------------------Arguments--------------------------------
      integer k1,
     $        k2
      real to3co2(plond),       ! pressure weighted temperature
     $     pnm(plond,plevp),    ! interface pressures
     $     ucfc11(plond,plevp), ! CFC11 path length
     $     ucfc12(plond,plevp), ! CFC12 path length
     $     un2o0(plond,plevp),  ! N2O path length
     $     un2o1(plond,plevp),  ! N2O path length (hot band)
     $     uch4(plond,plevp),   ! CH4 path length
     $     uco211(plond,plevp), ! CO2 9.4 micron band path length
     $     uco212(plond,plevp), ! CO2 9.4 micron band path length
     $     uco213(plond,plevp), ! CO2 9.4 micron band path length
     $     uco221(plond,plevp), ! CO2 10.4 micron band path length
     $     uco222(plond,plevp), ! CO2 10.4 micron band path length
     $     uco223(plond,plevp), ! CO2 10.4 micron band path length
     $     bn2o0(plond,plevp),  ! pressure factor for n2o
     $     bn2o1(plond,plevp),  ! pressure factor for n2o
     $     bch4(plond,plevp)    ! pressure factor for ch4
      real dw(plond),           ! h2o path length
     $     pnew(plond),         ! pressure
     $     s2c(plond,plevp),    ! continuum path length
     $     uptype(plond,plevp), ! p-type h2o path length
     $     dplh2o(plond)        ! p squared h2o path length
      real abplnk1(14,plond,plevp) ! Planck factor
      real tco2(plond),         ! co2 transmission factor
     $     th2o(plond),         ! h2o transmission factor
     $     to3(plond)           ! o3 transmission factor
c
c  Output Arguments
c
      real abstrc(plond)  ! total trace gas absorptivity
c
c  Local Variables
c
      real sqti(plond),         ! square root of mean temp
     $     du1,                 ! cfc11 path length
     $     du2,                 ! cfc12 path length
     $     acfc1,               ! cfc11 absorptivity 798 cm-1
     $     acfc2,               ! cfc11 absorptivity 846 cm-1
     $     acfc3,               ! cfc11 absorptivity 933 cm-1
     $     acfc4,               ! cfc11 absorptivity 1085 cm-1
     $     acfc5,               ! cfc12 absorptivity 889 cm-1
     $     acfc6,               ! cfc12 absorptivity 923 cm-1
     $     acfc7,               ! cfc12 absorptivity 1102 cm-1
     $     acfc8,               ! cfc12 absorptivity 1161 cm-1
     $     du01,                ! n2o path length
     $     dbeta01,             ! n2o pressure factor
     $     dbeta11,             !         "
     $     an2o1,               ! absorptivity of 1285 cm-1 n2o band
     $     du02,                ! n2o path length
     $     dbeta02,             ! n2o pressure factor
     $     an2o2,               ! absorptivity of 589 cm-1 n2o band
     $     du03                 ! n2o path length
      real dbeta03,             ! n2o pressure factor
     $     an2o3,               ! absorptivity of 1168 cm-1 n2o band
     $     duch4,               ! ch4 path length
     $     dbetac,              ! ch4 pressure factor
     $     ach4,                ! absorptivity of 1306 cm-1 ch4 band
     $     du11,                ! co2 path length
     $     du12,                !       "
     $     du13,                !       "
     $     dbetc1,              ! co2 pressure factor
     $     dbetc2,              ! co2 pressure factor
     $     aco21,               ! absorptivity of 1064 cm-1 band
     $     du21,                ! co2 path length
     $     du22,                !       "
     $     du23,                !       "
     $     aco22                ! absorptivity of 961 cm-1 band
      real tt(plond),           ! temp. factor for h2o overlap factor
     $     psi1,                !                 "
     $     phi1,                !                 "
     $     p1,                  ! h2o overlap factor
     $     w1,                  !        "
     $     ds2c(plond),         ! continuum path length
     $     duptyp(plond),       ! p-type path length
     $     tw(plond,6),         ! h2o transmission factor
     $     g1(6),               !         "
     $     g2(6),               !         "
     $     g3(6),               !         "
     $     g4(6),               !         "
     $     ab(6),               ! h2o temp. factor
     $     bb(6),               !         "
     $     abp(6),              !         "
     $     bbp(6)               !         "
      real tcfc3,               ! transmission for cfc11 band
     $     tcfc4,               ! transmission for cfc11 band
     $     tcfc6,               ! transmission for cfc12 band
     $     tcfc7,               ! transmission for cfc12 band
     $     tcfc8,               ! transmission for cfc12 band
     $     tlw,                 ! h2o transmission
     $     tch4                 ! ch4 transmission
c
      data g1 /0.0468556,0.0397454,0.0407664,0.0304380,0.0540398,
     $         0.0321962/
      data g2 /14.4832,4.30242,5.23523,3.25342,0.698935,16.5599/
      data g3 /26.1898,18.4476,15.3633,12.1927,9.14992,8.07092/
      data g4 /0.0261782,0.0369516,0.0307266,0.0243854,0.0182932,
     $         0.0161418/
      data ab /3.0857e-2,2.3524e-2,1.7310e-2,2.6661e-2,2.8074e-2,
     $         2.2915e-2/
      data bb /-1.3512e-4,-6.8320e-5,-3.2609e-5,-1.0228e-5,
     $         -9.5743e-5,-1.0304e-4/
      data abp/2.9129e-2,2.4101e-2,1.9821e-2,2.6904e-2,2.9458e-2,
     $         1.9892e-2/
      data bbp/-1.3139e-4,-5.5688e-5,-4.6380e-5,-8.0362e-5,
     $         -1.0115e-4,-8.8061e-5/
      integer i,l
c------------------------------------------------------------------
      real func, u, b
      func(u,b) = u/sqrt(4.0 + u*(1.0 + 1.0 / b))
c
      do i = 1,plon
         sqti(i) = sqrt(to3co2(i))
c h2o transmission 
         tt(i) = abs(to3co2(i) - 250.0)
         ds2c(i) = abs(s2c(i,k1) - s2c(i,k2))
         duptyp(i) = abs(uptype(i,k1) - uptype(i,k2))
      end do
c
      do l = 1,6
            do i = 1,plon
               psi1 = exp(abp(l)*tt(i)+bbp(l)*tt(i)*tt(i))
               phi1 = exp(ab(l)*tt(i)+bb(l)*tt(i)*tt(i))
               p1 = pnew(i) * (psi1/phi1) / sslp
               w1 = dw(i) * phi1
               tw(i,l) = exp(- g1(l)*p1*(sqrt(1.0+g2(l)*(w1/p1))-1.0)
     $                     - g3(l)*ds2c(i)-g4(l)*duptyp(i))
            end do
      end do
c
      do i = 1,plon
            du1 = abs(ucfc11(i,k1) - ucfc11(i,k2))
            du2 = abs(ucfc12(i,k1) - ucfc12(i,k2))
c cfc transmissions
            tcfc3 = exp(-175.005*du1)
            tcfc4 = exp(-1202.18*du1)
            tcfc6 = exp(-5786.73*du2)
            tcfc7 = exp(-2873.51*du2)
            tcfc8 = exp(-2085.59*du2)  
c  Absorptivity for CFC11 bands          
            acfc1 = 50.0*(1.0 - exp(-54.09*du1))*tw(i,1)*abplnk1(7,i,k2)
            acfc2 = 60.0*(1.0 - exp(-5130.03*du1))*tw(i,2)
     $                                            *abplnk1(8,i,k2)
            acfc3 = 60.0*(1.0 - tcfc3) * tw(i,4)*tcfc6*abplnk1(9,i,k2)
            acfc4 = 100.0*(1.0 - tcfc4)* tw(i,5) * abplnk1(10,i,k2)
c  Absorptivity for CFC12 bands
            acfc5 = 45.0*(1.0 - exp(-1272.35*du2))*tw(i,3)*
     $                                             abplnk1(11,i,k2)
            acfc6 = 50.0*(1.0 - tcfc6)* tw(i,4) * abplnk1(12,i,k2)
            acfc7 = 80.0*(1.0 - tcfc7)* tw(i,5) * tcfc4*abplnk1(13,i,k2)
            acfc8 = 70.0*(1.0 - tcfc8)* tw(i,6) * abplnk1(14,i,k2)
c  Emissivity for CH4 band 1306 cm-1
            tlw = exp(-1.0*sqrt(dplh2o(i)))
            duch4 = abs(uch4(i,k1) - uch4(i,k2))
            dbetac = abs(bch4(i,k1) - bch4(i,k2))/duch4
            ach4 = 6.00444*sqti(i)*alog(1.0 + func(duch4,dbetac)) *
     $           tlw * abplnk1(3,i,k2) 
            tch4 = 1.0/(1.0 + 0.02*func(duch4,dbetac))
c  Absorptivity for N2O bands 
            du01 = abs(un2o0(i,k1) - un2o0(i,k2))
            du11 = abs(un2o1(i,k1) - un2o1(i,k2))
            dbeta01 = abs(bn2o0(i,k1) - bn2o0(i,k2))/du01
            dbeta11 = abs(bn2o1(i,k1) - bn2o1(i,k2))/du11
c     1285 cm-1 band
            an2o1 = 2.35558*sqti(i)*alog(1.0 + func(du01,dbeta01)
     $            +  func(du11,dbeta11)) * tlw* tch4*abplnk1(4,i,k2)
            du02 = 0.100090*du01
            du12 = 0.0992746*du11
            dbeta02 = 0.964282*dbeta01
c     589 cm-1 band
            an2o2 = 2.65581*sqti(i)*alog(1.0 + func(du02,dbeta02)
     $            +  func(du12,dbeta02)) * th2o(i) * tco2(i) *
     $               abplnk1(5,i,k2)
            du03 = 0.0333767*du01
            dbeta03 = 0.982143*dbeta01
c     1168 cm-1 band
            an2o3 = 2.54034*sqti(i)*alog(1.0 + func(du03,dbeta03)) *
     $           tw(i,6) * tcfc8 * abplnk1(6,i,k2)
c  Emissivity for 1064 cm-1 band of CO2
            du11 = abs(uco211(i,k1) - uco211(i,k2))
            du12 = abs(uco212(i,k1) - uco212(i,k2))
            du13 = abs(uco213(i,k1) - uco213(i,k2))            
            dbetc1 = 2.97558*abs(pnm(i,k1) + pnm(i,k2))/
     $                                           (2.0*sslp*sqti(i))
            dbetc2 = 2.0 * dbetc1
            aco21 = 3.7571*sqti(i)*alog(1.0 + func(du11,dbetc1)
     $         + func(du12,dbetc2) + func(du13,dbetc2))
     $         * to3(i) * tw(i,5) * tcfc4 * tcfc7 * abplnk1(2,i,k2)
c  Emissivity for 961 cm-1 band
            du21 = abs(uco221(i,k1) - uco221(i,k2))
            du22 = abs(uco222(i,k1) - uco222(i,k2))
            du23 = abs(uco223(i,k1) - uco223(i,k2))
            aco22 = 3.8443*sqti(i)*alog(1.0 + func(du21,dbetc1)
     $         + func(du22,dbetc1) + func(du23,dbetc2))
     $         * tw(i,4) * tcfc3 * tcfc6 * abplnk1(1,i,k2)
c total trace gas absorptivity
            abstrc(i) = acfc1 + acfc2 + acfc3 + acfc4 + acfc5 + acfc6
     $            +  acfc7 + acfc8 + an2o1 + an2o2 + an2o3 + ach4
     $            +  aco21 + aco22
      end do
      return
      end


      subroutine trcabn(k2, kn, ucfc11, ucfc12, un2o0,  un2o1,
     $                          uch4,   uco211, uco212, uco213, 
     $                          uco221, uco222, uco223, bn2o0, 
     $                          bn2o1,  bch4,   tbar,   bplnk, 
     $                          winpl,  pinpl,  tco2,   th2o,
     $                          to3,    uptype, dw,     s2c,    
     $                          up2,    pnew,   abstrc)
c----------------------------------------------------------------------
c Calculate nearest layer absorptivity due to CH4, N2O, CFC11 and CFC12
c
c         Coded by J.T. Kiehl November 21, 1994
c-----------------------------------------------------------------------
c
c $Id: trcabn.F,v 1.2 1995/02/17 21:28:54 jhack Exp $
c $Author: jhack $
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
C------------------------------Parameters-------------------------------
c
c $Id: prgrid.h,v 1.2 1995/02/10 01:09:10 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation resolution and I/O parameters
C
      integer
     $     plon,    ! number of longitudes
     $     plev,    ! number of vertical levels
     $     plat,    ! number of latitudes
     $     pcnst,   ! number of constituents (including water vapor)
     $     plevmx,  ! number of subsurface levels
     $     plevp,   ! plev + 1
     $     nxpt,    ! no.of points outside active domain for interpolant
     $     jintmx,  ! number of extra latitudes in polar region
     $     plond,   ! slt extended domain longitude
     $     platd,   ! slt extended domain lat.
     $     plevd    ! fold plev,pcnst indices into one
      integer
     $     plngbuf  ! length of absorptivity/emissivity record
C
      parameter(plon    = 1,
     $          plev    = 18,
     $          plat    = 1,
     $          pcnst   = 1,
     $          plevmx  = 4,
     $          plevp   = plev + 1,
     $          nxpt    = 1,
     $          jintmx  = 1,
     $          plond   = plon + 1 + 2*nxpt,
     $          platd   = plat + 2*nxpt + 2*jintmx,
     $          plevd   = plev*(3 + pcnst))
      parameter(plngbuf = 512*((plond*plevp*plevp + plond*plev*4 +
     $                          plond*plevp)/512 + 1))
C
C------------------------------Commons----------------------------------
c
c $Id: crdcon.h,v 1.1.1.1 1995/02/09 23:26:46 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation constants
C
      common/crdcon/gravit  ,rga     ,cpair   ,epsilo  ,sslp    ,
     $              stebol  ,rgsslp  ,co2vmr  ,dpfo3   ,dpfco2  ,
     $              dayspy  ,pie
C
      real gravit,    ! Acceleration of gravity
     $     rga,       ! 1./gravit
     $     cpair,     ! Specific heat of dry air
     $     epsilo,    ! Ratio of mol. wght of H2O to dry air
     $     sslp,      ! Standard sea-level pressure
     $     stebol,    ! Stefan-Boltzmann's constant
     $     rgsslp,    ! 0.5/(gravit*sslp)
     $     co2vmr,    ! CO2 volume mixing ratio
     $     dpfo3,     ! Voigt correction factor for O3
     $     dpfco2,    ! Voigt correction factor for CO2
     $     dayspy,    ! Number of days per 1 year
     $     pie        ! 3.14.....
C
C------------------------------Arguments--------------------------------
      integer k2,
     $        kn
      real tbar(plond,4),       ! pressure weighted temperature
     $     ucfc11(plond,plevp), ! CFC11 path length
     $     ucfc12(plond,plevp), ! CFC12 path length
     $     un2o0(plond,plevp),  ! N2O path length
     $     un2o1(plond,plevp),  ! N2O path length (hot band)
     $     uch4(plond,plevp),   ! CH4 path length
     $     uco211(plond,plevp), ! CO2 9.4 micron band path length
     $     uco212(plond,plevp), ! CO2 9.4 micron band path length
     $     uco213(plond,plevp), ! CO2 9.4 micron band path length
     $     uco221(plond,plevp), ! CO2 10.4 micron band path length
     $     uco222(plond,plevp), ! CO2 10.4 micron band path length
     $     uco223(plond,plevp), ! CO2 10.4 micron band path length
     $     bn2o0(plond,plevp),  ! pressure factor for n2o
     $     bn2o1(plond,plevp),  ! pressure factor for n2o
     $     bch4(plond,plevp),   ! pressure factor for ch4
     $     bplnk(14,plond,4),   ! weighted Planck function for absorptivity
     $     winpl(plond,4),      ! fractional path length
     $     pinpl(plond,4)       ! pressure factor for subdivided layer
      real tco2(plond),         ! co2 transmission 
     $     th2o(plond),         ! h2o transmission
     $     to3(plond)           ! o3 transmission
      real dw(plond),           ! h2o path length
     $     pnew(plond),         ! pressure factor
     $     s2c(plond,plevp),    ! h2o continuum factor
     $     uptype(plond,plevp), ! p-type path length
     $     up2(plond)           ! p squared path length
c
c  Output Arguments
c
      real abstrc(plond)        ! total trace gas absorptivity
c
c  Local Variables
c
      real sqti(plond),         ! square root of mean temp
     $     rsqti(plond),        ! reciprocal of sqti
     $     du1,                 ! cfc11 path length
     $     du2,                 ! cfc12 path length
     $     acfc1,               ! absorptivity of cfc11 798 cm-1 band
     $     acfc2,               ! absorptivity of cfc11 846 cm-1 band
     $     acfc3,               ! absorptivity of cfc11 933 cm-1 band
     $     acfc4,               ! absorptivity of cfc11 1085 cm-1 band
     $     acfc5,               ! absorptivity of cfc11 889 cm-1 band
     $     acfc6,               ! absorptivity of cfc11 923 cm-1 band
     $     acfc7,               ! absorptivity of cfc11 1102 cm-1 band
     $     acfc8,               ! absorptivity of cfc11 1161 cm-1 band
     $     du01,                ! n2o path length
     $     dbeta01,             ! n2o pressure factors
     $     dbeta11              !        "
      real  an2o1,              ! absorptivity of the 1285 cm-1 n2o band
     $     du02,                ! n2o path length
     $     dbeta02,             ! n2o pressure factor
     $     an2o2,               ! absorptivity of the 589 cm-1 n2o band
     $     du03,                ! n2o path length
     $     dbeta03,             ! n2o pressure factor
     $     an2o3,               ! absorptivity of the 1168 cm-1 n2o band
     $     duch4,               ! ch4 path length
     $     dbetac,              ! ch4 pressure factor
     $     ach4,                ! absorptivity of the 1306 cm-1 ch4 band
     $     du11,                ! co2 path length
     $     du12,                !       "
     $     du13,                !       "
     $     dbetc1,              ! co2 pressure factor
     $     dbetc2,              ! co2 pressure factor
     $     aco21,               ! absorptivity of the 1064 cm-1 co2 band
     $     du21,                ! co2 path length
     $     du22,                !       "
     $     du23,                !       "
     $     aco22                ! absorptivity of the 961 cm-1 co2 band
      real tt(plond),           ! temp. factor for h2o overlap
     $     psi1,                !          "
     $     phi1,                !          "
     $     p1,                  ! factor for h2o overlap
     $     w1,                  !          "
     $     ds2c(plond),         ! continuum path length
     $     duptyp(plond),       ! p-type path length
     $     tw(plond,6),         ! h2o transmission overlap
     $     g1(6),               ! h2o overlap factor
     $     g2(6),               !         "
     $     g3(6),               !         "
     $     g4(6),               !         "
     $     ab(6),               ! h2o temp. factor
     $     bb(6),               !         "
     $     abp(6),              !         "  
     $     bbp(6)               !         "
      real tcfc3,               ! transmission of cfc11 band
     $     tcfc4,               ! transmission of cfc11 band
     $     tcfc6,               ! transmission of cfc12 band
     $     tcfc7,               !         "
     $     tcfc8,               !         "
     $     tlw,                 ! h2o transmission
     $     tch4                 ! ch4 transmission
      data g1 /0.0468556,0.0397454,0.0407664,0.0304380,0.0540398,
     $         0.0321962/
      data g2 /14.4832,4.30242,5.23523,3.25342,0.698935,16.5599/
      data g3 /26.1898,18.4476,15.3633,12.1927,9.14992,8.07092/
      data g4 /0.0261782,0.0369516,0.0307266,0.0243854,0.0182932,
     $         0.0161418/
      data ab /3.0857e-2,2.3524e-2,1.7310e-2,2.6661e-2,2.8074e-2,
     $         2.2915e-2/
      data bb /-1.3512e-4,-6.8320e-5,-3.2609e-5,-1.0228e-5,
     $         -9.5743e-5,-1.0304e-4/
      data abp/2.9129e-2,2.4101e-2,1.9821e-2,2.6904e-2,2.9458e-2,
     $         1.9892e-2/
      data bbp/-1.3139e-4,-5.5688e-5,-4.6380e-5,-8.0362e-5,
     $         -1.0115e-4,-8.8061e-5/
      integer i,l
c------------------------------------------------------------------
      real func, u, b
      func(u,b) = u/sqrt(4.0 + u*(1.0 + 1.0 / b))
c
      do i = 1,plon
         sqti(i) = sqrt(tbar(i,kn))
         rsqti(i) = 1. / sqti(i)
c h2o transmission
         tt(i) = abs(tbar(i,kn) - 250.0)
         ds2c(i) = abs(s2c(i,k2+1) - s2c(i,k2))
         duptyp(i) = abs(uptype(i,k2+1) - uptype(i,k2))
      end do
c
      do l = 1,6
         do i = 1,plon
            psi1 = exp(abp(l)*tt(i)+bbp(l)*tt(i)*tt(i))
            phi1 = exp(ab(l)*tt(i)+bb(l)*tt(i)*tt(i))
            p1 = pnew(i) * (psi1/phi1) / sslp
            w1 = dw(i) * winpl(i,kn) * phi1
            tw(i,l) = exp(- g1(l)*p1*(sqrt(1.0+g2(l)*(w1/p1))-1.0)
     $                  - g3(l)*ds2c(i)-g4(l)*duptyp(i))
         end do
      end do
c
      do i = 1,plon
c
            du1 = abs(ucfc11(i,k2+1) - ucfc11(i,k2)) * winpl(i,kn)
            du2 = abs(ucfc12(i,k2+1) - ucfc12(i,k2)) * winpl(i,kn)
c cfc transmissions
            tcfc3 = exp(-175.005*du1)
            tcfc4 = exp(-1202.18*du1)
            tcfc6 = exp(-5786.73*du2)
            tcfc7 = exp(-2873.51*du2)
            tcfc8 = exp(-2085.59*du2)
c  Absorptivity for CFC11 bands
            acfc1 = 50.0*(1.0 - exp(-54.09*du1)) * tw(i,1)*bplnk(7,i,kn)
            acfc2 = 60.0*(1.0 - exp(-5130.03*du1))*tw(i,2)*bplnk(8,i,kn)
            acfc3 = 60.0*(1.0 - tcfc3)*tw(i,4)*tcfc6 * bplnk(9,i,kn)
            acfc4 = 100.0*(1.0 - tcfc4)* tw(i,5) * bplnk(10,i,kn)
c  Absorptivity for CFC12 bands
            acfc5 = 45.0*(1.0 - exp(-1272.35*du2))*tw(i,3)
     $                                            *bplnk(11,i,kn)
            acfc6 = 50.0*(1.0 - tcfc6)*tw(i,4)*bplnk(12,i,kn)
            acfc7 = 80.0*(1.0 - tcfc7)* tw(i,5)*tcfc4 *bplnk(13,i,kn)
            acfc8 = 70.0*(1.0 - tcfc8)*tw(i,6)*bplnk(14,i,kn)
c  Emissivity for CH4 band 1306 cm-1
            tlw = exp(-1.0*sqrt(up2(i)))
            duch4 = abs(uch4(i,k2+1) - uch4(i,k2)) * winpl(i,kn)
            dbetac = 2.94449 * pinpl(i,kn) * rsqti(i) / sslp
            ach4 = 6.00444*sqti(i)*alog(1.0 + func(duch4,dbetac)) *
     $           tlw * bplnk(3,i,kn) 
            tch4 = 1.0/(1.0 + 0.02*func(duch4,dbetac)) 
c  Absorptivity for N2O bands 
            du01 = abs(un2o0(i,k2+1) - un2o0(i,k2)) * winpl(i,kn)
            du11 = abs(un2o1(i,k2+1) - un2o1(i,k2)) * winpl(i,kn)
            dbeta01 = 19.399 *  pinpl(i,kn) * rsqti(i) / sslp
            dbeta11 = dbeta01
c     1285 cm-1 band
            an2o1 = 2.35558*sqti(i)*alog(1.0 + func(du01,dbeta01)
     $            + func(du11,dbeta11)) * tlw * tch4 * bplnk(4,i,kn)
            du02 = 0.100090*du01
            du12 = 0.0992746*du11
            dbeta02 = 0.964282*dbeta01
c     589 cm-1 band
            an2o2 = 2.65581*sqti(i)*alog(1.0 + func(du02,dbeta02)
     $             +  func(du12,dbeta02)) * tco2(i) * th2o(i) *
     $                bplnk(5,i,kn)
            du03 = 0.0333767*du01
            dbeta03 = 0.982143*dbeta01
c     1168 cm-1 band
            an2o3 = 2.54034*sqti(i)*alog(1.0 + func(du03,dbeta03)) *
     $           tw(i,6) * tcfc8 * bplnk(6,i,kn)
c  Emissivity for 1064 cm-1 band of CO2
            du11 = abs(uco211(i,k2+1) - uco211(i,k2)) * winpl(i,kn)
            du12 = abs(uco212(i,k2+1) - uco212(i,k2)) * winpl(i,kn)
            du13 = abs(uco213(i,k2+1) - uco213(i,k2)) * winpl(i,kn)
            dbetc1 = 2.97558 * pinpl(i,kn) * rsqti(i) / sslp
            dbetc2 = 2.0 * dbetc1
            aco21 = 3.7571*sqti(i)*alog(1.0 + func(du11,dbetc1)
     $         + func(du12,dbetc2) + func(du13,dbetc2))
     $         * to3(i) * tw(i,5) * tcfc4 * tcfc7 * bplnk(2,i,kn)
c  Emissivity for 961 cm-1 band of co2
            du21 = abs(uco221(i,k2+1) - uco221(i,k2)) * winpl(i,kn)
            du22 = abs(uco222(i,k2+1) - uco222(i,k2)) * winpl(i,kn)
            du23 = abs(uco223(i,k2+1) - uco223(i,k2)) * winpl(i,kn)
            aco22 = 3.8443*sqti(i)*alog(1.0 + func(du21,dbetc1)
     $         + func(du22,dbetc1) + func(du23,dbetc2))
     $         * tw(i,4) * tcfc3 * tcfc6 * bplnk(1,i,kn)
c total trace gas absorptivity
            abstrc(i) = acfc1 + acfc2 + acfc3 + acfc4 + acfc5 + acfc6
     $                + acfc7 + acfc8 + an2o1 + an2o2 + an2o3 + ach4
     $                + aco21 + aco22 
      end do
      return
      end

      subroutine trcems(k,      co2t,   pnm,    ucfc11, ucfc12, un2o0, 
     $                  un2o1,  bn2o0,  bn2o1,  uch4,   bch4,   uco211,
     $                  uco212, uco213, uco221, uco222, uco223, uptype,
     $                  w,      s2c,    up2,    emplnk, th2o,   tco2,   
     $                  to3,    emstrc )
c----------------------------------------------------------------------
c  Calculate emissivity for CH4, N2O, CFC11 and CFC12 bands.
c            Coded by J.T. Kiehl November 21, 1994
c-----------------------------------------------------------------------
c
c $Id: trcems.F,v 1.2 1995/02/17 21:28:56 jhack Exp $
c $Author: jhack $
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
C------------------------------Parameters-------------------------------
c
c $Id: prgrid.h,v 1.2 1995/02/10 01:09:10 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation resolution and I/O parameters
C
      integer
     $     plon,    ! number of longitudes
     $     plev,    ! number of vertical levels
     $     plat,    ! number of latitudes
     $     pcnst,   ! number of constituents (including water vapor)
     $     plevmx,  ! number of subsurface levels
     $     plevp,   ! plev + 1
     $     nxpt,    ! no.of points outside active domain for interpolant
     $     jintmx,  ! number of extra latitudes in polar region
     $     plond,   ! slt extended domain longitude
     $     platd,   ! slt extended domain lat.
     $     plevd    ! fold plev,pcnst indices into one
      integer
     $     plngbuf  ! length of absorptivity/emissivity record
C
      parameter(plon    = 1,
     $          plev    = 18,
     $          plat    = 1,
     $          pcnst   = 1,
     $          plevmx  = 4,
     $          plevp   = plev + 1,
     $          nxpt    = 1,
     $          jintmx  = 1,
     $          plond   = plon + 1 + 2*nxpt,
     $          platd   = plat + 2*nxpt + 2*jintmx,
     $          plevd   = plev*(3 + pcnst))
      parameter(plngbuf = 512*((plond*plevp*plevp + plond*plev*4 +
     $                          plond*plevp)/512 + 1))
C
C------------------------------Commons----------------------------------
c
c $Id: crdcon.h,v 1.1.1.1 1995/02/09 23:26:46 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation constants
C
      common/crdcon/gravit  ,rga     ,cpair   ,epsilo  ,sslp    ,
     $              stebol  ,rgsslp  ,co2vmr  ,dpfo3   ,dpfco2  ,
     $              dayspy  ,pie
C
      real gravit,    ! Acceleration of gravity
     $     rga,       ! 1./gravit
     $     cpair,     ! Specific heat of dry air
     $     epsilo,    ! Ratio of mol. wght of H2O to dry air
     $     sslp,      ! Standard sea-level pressure
     $     stebol,    ! Stefan-Boltzmann's constant
     $     rgsslp,    ! 0.5/(gravit*sslp)
     $     co2vmr,    ! CO2 volume mixing ratio
     $     dpfo3,     ! Voigt correction factor for O3
     $     dpfco2,    ! Voigt correction factor for CO2
     $     dayspy,    ! Number of days per 1 year
     $     pie        ! 3.14.....
C
C------------------------------Arguments--------------------------------
      real co2t(plond,plevp),   ! pressure weighted temperature
     $     pnm(plond,plevp),    ! interface pressure
     $     ucfc11(plond,plevp), ! CFC11 path length
     $     ucfc12(plond,plevp), ! CFC12 path length
     $     un2o0(plond,plevp),  ! N2O path length
     $     un2o1(plond,plevp),  ! N2O path length (hot band)
     $     uch4(plond,plevp),   ! CH4 path length
     $     uco211(plond,plevp), ! CO2 9.4 micron band path length
     $     uco212(plond,plevp), ! CO2 9.4 micron band path length
     $     uco213(plond,plevp), ! CO2 9.4 micron band path length
     $     uco221(plond,plevp), ! CO2 10.4 micron band path length
     $     uco222(plond,plevp), ! CO2 10.4 micron band path length
     $     uco223(plond,plevp), ! CO2 10.4 micron band path length
     $     uptype(plond,plevp), ! continuum path length
     $     bn2o0(plond,plevp),  ! pressure factor for n2o
     $     bn2o1(plond,plevp),  ! pressure factor for n2o
     $     bch4(plond,plevp)    ! pressure factor for ch4
      real emplnk(14,plond)     ! emissivity Planck factor
      real th2o(plond),         ! water vapor overlap factor
     $     tco2(plond),         ! co2 overlap factor
     $     to3(plond)           ! o3 overlap factor
      real s2c(plond,plevp),    ! h2o continuum path length
     $     w(plond,plevp),      ! h2o path length
     $     up2(plond)           ! pressure squared h2o path length
      integer k                 ! level index
c
c  Output Arguments
c
      real emstrc(plond,plevp)  ! total trace gas emissivity
c
c  Local Variables
c
      real sqti(plond),         ! square root of mean temp
     $     ecfc1,               ! emissivity of cfc11 798 cm-1 band
     $     ecfc2,               !     "      "    "   846 cm-1 band
     $     ecfc3,               !     "      "    "   933 cm-1 band
     $     ecfc4,               !     "      "    "   1085 cm-1 band
     $     ecfc5,               !     "      "  cfc12 889 cm-1 band
     $     ecfc6,               !     "      "    "   923 cm-1 band
     $     ecfc7,               !     "      "    "   1102 cm-1 band
     $     ecfc8,               !     "      "    "   1161 cm-1 band
     $     u01,                 ! n2o path length
     $     u11,                 ! n2o path length
     $     beta01,              ! n2o pressure factor
     $     beta11,              ! n2o pressure factor
     $     en2o1,               ! emissivity of the 1285 cm-1 N2O band
     $     u02,                 ! n2o path length
     $     u12,                 ! n2o path length
     $     beta02,              ! n2o pressure factor
     $     en2o2,               ! emissivity of the 589 cm-1 N2O band
     $     u03                  ! n2o path length
      real beta03,              ! n2o pressure factor
     $     en2o3,               ! emissivity of the 1168 cm-1 N2O band
     $     betac,               ! ch4 pressure factor
     $     ech4,                ! emissivity of 1306 cm-1 CH4 band
     $     betac1,              ! co2 pressure factor
     $     betac2,              ! co2 pressure factor
     $     eco21,               ! emissivity of 1064 cm-1 CO2 band
     $     eco22                ! emissivity of 961 cm-1 CO2 band
      real tt(plond),           ! temp. factor for h2o overlap factor
     $     psi1,                ! narrow band h2o temp. factor
     $     phi1,                !             "
     $     p1,                  ! h2o line overlap factor
     $     w1,                  !          "
     $     tw(plond,6),         ! h2o transmission overlap
     $     g1(6),               ! h2o overlap factor
     $     g2(6),               !          "
     $     g3(6),               !          "
     $     g4(6),               !          "
     $     ab(6),               !          "
     $     bb(6),               !          "
     $     abp(6),              !          "
     $     bbp(6)               !          "
      real tcfc3,               ! transmission for cfc11 band
     $     tcfc4,               !          "
     $     tcfc6,               ! transmission for cfc12 band
     $     tcfc7,               !          "
     $     tcfc8,               !          "
     $     tlw,                 ! h2o overlap factor
     $     tch4                 ! ch4 overlap factor
      data g1 /0.0468556,0.0397454,0.0407664,0.0304380,0.0540398,
     $         0.0321962/
      data g2 /14.4832,4.30242,5.23523,3.25342,0.698935,16.5599/
      data g3 /26.1898,18.4476,15.3633,12.1927,9.14992,8.07092/
      data g4 /0.0261782,0.0369516,0.0307266,0.0243854,0.0182932,
     $         0.0161418/
      data ab /3.0857e-2,2.3524e-2,1.7310e-2,2.6661e-2,2.8074e-2,
     $         2.2915e-2/
      data bb /-1.3512e-4,-6.8320e-5,-3.2609e-5,-1.0228e-5,
     $         -9.5743e-5,-1.0304e-4/
      data abp/2.9129e-2,2.4101e-2,1.9821e-2,2.6904e-2,2.9458e-2,
     $         1.9892e-2/
      data bbp/-1.3139e-4,-5.5688e-5,-4.6380e-5,-8.0362e-5,
     $         -1.0115e-4,-8.8061e-5/
      integer i,l
      real func, u, b
      func(u,b) = u/sqrt(4.0 + u*(1.0 + 1.0 / b))
c
      do i = 1,plon
         sqti(i) = sqrt(co2t(i,k))
c Transmission for h2o
         tt(i) = abs(co2t(i,k) - 250.0)
      end do
c
      do l = 1,6
         do i = 1,plon
            psi1 = exp(abp(l)*tt(i)+bbp(l)*tt(i)*tt(i))
            phi1 = exp(ab(l)*tt(i)+bb(l)*tt(i)*tt(i))
            p1 = pnm(i,k) * (psi1/phi1) / sslp
            w1 = w(i,k) * phi1
            tw(i,l) = exp(- g1(l)*p1*(sqrt(1.0+g2(l)*(w1/p1))-1.0)
     $                  - g3(l)*s2c(i,k)-g4(l)*uptype(i,k))
         end do
      end do 
c
      do i = 1,plon
c transmission due to cfc bands
            tcfc3 = exp(-175.005*ucfc11(i,k))
            tcfc4 = exp(-1202.18*ucfc11(i,k)) 
            tcfc6 = exp(-5786.73*ucfc12(i,k))
            tcfc7 = exp(-2873.51*ucfc12(i,k))
            tcfc8 = exp(-2085.59*ucfc12(i,k))
c  Emissivity for CFC11 bands
            ecfc1 = 50.0*(1.0 - exp(-54.09*ucfc11(i,k))) * tw(i,1) * 
     $                                                  emplnk(7,i)
            ecfc2 = 60.0*(1.0 - exp(-5130.03*ucfc11(i,k)))* tw(i,2) *
     $                                                  emplnk(8,i)
            ecfc3 = 60.0*(1.0 - tcfc3)*tw(i,4)*tcfc6*emplnk(9,i)
            ecfc4 = 100.0*(1.0 - tcfc4)*tw(i,5)*emplnk(10,i)
c  Emissivity for CFC12 bands
            ecfc5 = 45.0*(1.0 - exp(-1272.35*ucfc12(i,k)))*tw(i,3)*
     $                                                     emplnk(11,i)
            ecfc6 = 50.0*(1.0 - tcfc6)*tw(i,4)*emplnk(12,i)
            ecfc7 = 80.0*(1.0 - tcfc7)*tw(i,5)* tcfc4 * emplnk(13,i)
            ecfc8 = 70.0*(1.0 - tcfc8)*tw(i,6) * emplnk(14,i)
c  Emissivity for CH4 band 1306 cm-1
            tlw = exp(-1.0*sqrt(up2(i)))
            betac = bch4(i,k)/uch4(i,k)
            ech4 = 6.00444*sqti(i)*alog(1.0 + func(uch4(i,k),betac)) *
     $                 tlw * emplnk(3,i)
            tch4 = 1.0/(1.0 + 0.02*func(uch4(i,k),betac))
c  Emissivity for N2O bands 
            u01 = un2o0(i,k)
            u11 = un2o1(i,k)
            beta01 = bn2o0(i,k)/un2o0(i,k)
            beta11 = bn2o1(i,k)/un2o1(i,k)
c     1285 cm-1 band
            en2o1 = 2.35558*sqti(i)*alog(1.0 + func(u01,beta01) +
     $              func(u11,beta11))*tlw*tch4*emplnk(4,i)
            u02 = 0.100090*u01
            u12 = 0.0992746*u11
            beta02 = 0.964282*beta01
c     589 cm-1 band
            en2o2 = 2.65581*sqti(i)*alog(1.0 + func(u02,beta02) +
     $              func(u12,beta02)) * tco2(i) * th2o(i) * 
     $              emplnk(5,i)
            u03 = 0.0333767*u01
            beta03 = 0.982143*beta01
c     1168 cm-1 band
            en2o3 = 2.54034*sqti(i)*alog(1.0 + func(u03,beta03)) *
     $                 tw(i,6) * tcfc8 * emplnk(6,i)
c  Emissivity for 1064 cm-1 band of CO2
            betac1 = 2.97558*pnm(i,k) / (sslp*sqti(i))
            betac2 = 2.0 * betac1
            eco21 = 3.7571*sqti(i)*alog(1.0 + func(uco211(i,k),betac1)
     $         + func(uco212(i,k),betac2) + func(uco213(i,k),betac2))
     $         * to3(i) * tw(i,5) * tcfc4 * tcfc7 * emplnk(2,i)
c  Emissivity for 961 cm-1 band
            eco22 = 3.8443*sqti(i)*alog(1.0 + func(uco221(i,k),betac1)
     $         + func(uco222(i,k),betac1) + func(uco223(i,k),betac2))
     $         * tw(i,4) * tcfc3 * tcfc6 * emplnk(1,i)
c total trace gas emissivity
            emstrc(i,k) = ecfc1 + ecfc2 + ecfc3 + ecfc4 + ecfc5 +ecfc6
     $                + ecfc7 + ecfc8 + en2o1 + en2o2 + en2o3 + ech4
     $                + eco21 + eco22
      end do
      return
      end
      subroutine trcmix(pmid, clat, coslat, n2o, ch4, cfc11, cfc12)
c-------------------------------------------------------------
c Specify zonal mean mass mixing ratios of CH4, N2O, CFC11 and
c CFC12
c          Code: J.T.Kiehl November 21, 1994
C------------------------------Parameters-------------------------------
c
c $Id: trcmix.F,v 1.2 1995/02/17 21:28:58 jhack Exp $
c $Author: jhack $
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
c
c $Id: pmgrid.h,v 1.2 1995/02/10 01:09:06 ccm2 Exp $
c $Author: ccm2 $
c
C
C Basic grid point resolution parameters
C
      integer
     $     plon,    ! number of longitudes
     $     plev,    ! number of vertical levels
     $     plat,    ! number of latitudes
     $     pcnst,   ! number of constituents (including water vapor)
     $     plevmx,  ! number of subsurface levels
     $     plevp,   ! plev + 1
     $     nxpt,    ! no.of points outside active domain for interpolant
     $     jintmx,  ! number of extra latitudes in polar region
     $     plond,   ! slt extended domain longitude
     $     platd,   ! slt extended domain lat.
     $     plevd    ! fold plev,pcnst indices into one
C
      parameter(plon   = 1,
     $          plev   = 18,
     $          plat   = 1,
     $          pcnst  = 1,
     $          plevmx = 4,
     $          plevp  = plev + 1,
     $          nxpt   = 1,
     $          jintmx = 1,
     $          plond  = plon + 1 + 2*nxpt,
     $          platd  = plat + 2*nxpt + 2*jintmx,
     $          plevd  = plev*(3 + pcnst))
C
c------------------------------input------------------------------------
      real pmid(plond,plev),       ! model pressures
     $     clat,                   ! current latitude in radians
     $     coslat                  ! cosine of latitude
c------------------------------output-----------------------------------
      real n2o(plond,plev),        ! nitrous oxide mass mixing ratio
     $     ch4(plond,plev),        ! methane mass mixing ratio
     $     cfc11(plond,plev),      ! cfc11 mass mixing ratio
     $     cfc12(plond,plev)       ! cfc12 mass mixing ratio
c------------------------------local------------------------------------
      integer i,                   ! longitude loop index
     $        k                    ! level index
      real dlat                    ! latitude in degrees
      real xn2o,                   ! pressure scale height for n2o
     $     xch4,                   ! pressure scale height for ch4
     $     xcfc11,                 ! pressure scale height for cfc11
     $     xcfc12                  ! pressure scale height for cfc12
      real ch40,                   ! tropospheric mass mixing ratio for ch4
     $     n2o0,                   ! tropospheric mass mixing ratio for n2o
     $     cfc110,                 ! tropospheric mass mixing ratio for cfc11
     $     cfc120                  ! tropospheric mass mixing ratio for cfc12
      real ptrop,                  ! pressure level of tropopause
     $     pratio                  ! pressure divided by ptrop
c
c tropospheric mass mixing ratios
      ch40 = 0.55241 * 1.714e-6
      n2o0 = 1.51913 * 0.311e-6
      cfc110 = 4.69548 * 0.280e-9
      cfc120 = 4.14307 * 0.503e-9
c set stratospheric scale height factor for gases
      dlat = abs(57.2958 * clat)
      if(dlat.le.45.0) then
        xn2o = 0.3478 + 0.00116 * dlat
        xch4 = 0.2353
        xcfc11 = 0.7273 + 0.00606 * dlat
        xcfc12 = 0.4000 + 0.00222 * dlat
      else
        xn2o = 0.4000 + 0.013333 * (dlat - 45)
        xch4 = 0.2353 + 0.0225489 * (dlat - 45)
        xcfc11 = 1.00 + 0.013333 * (dlat - 45)
        xcfc12 = 0.50 + 0.024444 * (dlat - 45)
      end if
c pressure of tropopause
      ptrop = 250.0e2 - 150.0e2*coslat**2.0
c
      do k = 1,plev
         do i = 1,plon
            if(pmid(i,k).ge.ptrop) then
              ch4(i,k) = ch40
              n2o(i,k) = n2o0
              cfc11(i,k) = cfc110
              cfc12(i,k) = cfc120
            else
              pratio = pmid(i,k)/ptrop
              ch4(i,k) = ch40 * (pratio)**xch4
              n2o(i,k) = n2o0 * (pratio)**xn2o
              cfc11(i,k) = cfc110 * (pratio)**xcfc11
              cfc12(i,k) = cfc120 * (pratio)**xcfc12
            end if
         end do
      end do
      return
      end
      subroutine trcplk(tint, tlayr, tplnke, emplnk, abplnk1, abplnk2)
c----------------------------------------------------------------------
c   Calculate Planck factors for absorptivity and emissivity of
c   CH4, N2O, CFC11 and CFC12
c
c-----------------------------------------------------------------------
c
c $Id: trcplk.F,v 1.2 1995/02/17 21:29:00 jhack Exp $
c $Author: jhack $
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
C------------------------------Parameters-------------------------------
c
c $Id: prgrid.h,v 1.2 1995/02/10 01:09:10 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation resolution and I/O parameters
C
      integer
     $     plon,    ! number of longitudes
     $     plev,    ! number of vertical levels
     $     plat,    ! number of latitudes
     $     pcnst,   ! number of constituents (including water vapor)
     $     plevmx,  ! number of subsurface levels
     $     plevp,   ! plev + 1
     $     nxpt,    ! no.of points outside active domain for interpolant
     $     jintmx,  ! number of extra latitudes in polar region
     $     plond,   ! slt extended domain longitude
     $     platd,   ! slt extended domain lat.
     $     plevd    ! fold plev,pcnst indices into one
      integer
     $     plngbuf  ! length of absorptivity/emissivity record
C
      parameter(plon    = 1,
     $          plev    = 18,
     $          plat    = 1,
     $          pcnst   = 1,
     $          plevmx  = 4,
     $          plevp   = plev + 1,
     $          nxpt    = 1,
     $          jintmx  = 1,
     $          plond   = plon + 1 + 2*nxpt,
     $          platd   = plat + 2*nxpt + 2*jintmx,
     $          plevd   = plev*(3 + pcnst))
      parameter(plngbuf = 512*((plond*plevp*plevp + plond*plev*4 +
     $                          plond*plevp)/512 + 1))
C
C------------------------------Commons----------------------------------
c
c $Id: crdcon.h,v 1.1.1.1 1995/02/09 23:26:46 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation constants
C
      common/crdcon/gravit  ,rga     ,cpair   ,epsilo  ,sslp    ,
     $              stebol  ,rgsslp  ,co2vmr  ,dpfo3   ,dpfco2  ,
     $              dayspy  ,pie
C
      real gravit,    ! Acceleration of gravity
     $     rga,       ! 1./gravit
     $     cpair,     ! Specific heat of dry air
     $     epsilo,    ! Ratio of mol. wght of H2O to dry air
     $     sslp,      ! Standard sea-level pressure
     $     stebol,    ! Stefan-Boltzmann's constant
     $     rgsslp,    ! 0.5/(gravit*sslp)
     $     co2vmr,    ! CO2 volume mixing ratio
     $     dpfo3,     ! Voigt correction factor for O3
     $     dpfco2,    ! Voigt correction factor for CO2
     $     dayspy,    ! Number of days per 1 year
     $     pie        ! 3.14.....
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real tint(plond,plevp),  ! interface temperatures
     $     tlayr(plond,plevp), ! k-1 level temperatures
     $     tplnke(plond)       ! Top Layer temperature
c
c output arguments
c
      real emplnk(14,plond),        ! emissivity Planck factor
     $     abplnk1(14,plond,plevp), ! non-nearest layer Plack factor
     $     abplnk2(14,plond,plevp)  ! nearest layer factor
c
c local workspace
c
      integer wvl                   ! wavelength index
      integer i,k
      real f1(14),                  ! Planck function factor
     $     f2(14),                  !        "
     $     f3(14)                   !        "
c
      data f1 /5.85713e8,7.94950e8,1.47009e9,1.40031e9,1.34853e8,
     $         1.05158e9,3.35370e8,3.99601e8,5.35994e8,8.42955e8,
     $         4.63682e8,5.18944e8,8.83202e8,1.03279e9/
      data f2 /2.02493e11,3.04286e11,6.90698e11,6.47333e11,
     $         2.85744e10,4.41862e11,9.62780e10,1.21618e11,
     $         1.79905e11,3.29029e11,1.48294e11,1.72315e11,
     $         3.50140e11,4.31364e11/
      data f3 /1383.0,1531.0,1879.0,1849.0,848.0,1681.0,
     $         1148.0,1217.0,1343.0,1561.0,1279.0,1328.0,
     $         1586.0,1671.0/
c
c Calculate emissivity Planck factor
c
      do wvl = 1,14
         do i = 1,plon
            emplnk(wvl,i) = f1(wvl)/
     $                   (tplnke(i)**4.0*(exp(f3(wvl)/tplnke(i))-1.0))
         end do
      end do
c
c Calculate absorptivity Planck factor for tint and tlayr temperatures
c
      do wvl = 1,14
         do k = 1, plevp
            do i = 1, plon
c non-nearlest layer function
               abplnk1(wvl,i,k) = (f2(wvl)*exp(f3(wvl)/tint(i,k)))
     $              /(tint(i,k)**5.0*(exp(f3(wvl)/tint(i,k))-1.0)**2.0)
c nearest layer function
               abplnk2(wvl,i,k) = (f2(wvl)*exp(f3(wvl)/tlayr(i,k)))
     $            /(tlayr(i,k)**5.0*(exp(f3(wvl)/tlayr(i,k))-1.0)**2.0)
            end do
         end do
      end do
      return
      end
      subroutine trcpth(tnm, pnm, cfc11, cfc12, n2o, ch4, qnm,
     $                  ucfc11, ucfc12, un2o0,  un2o1,  uch4, 
     $                  uco211, uco212, uco213, uco221, uco222, 
     $                  uco223, bn2o0,  bn2o1,  bch4,   uptype)
c----------------------------------------------------------------------
c Calculate path lengths and pressure factors for CH4, N2O, CFC11
c and CFC12. 
c           Coded by J.T. Kiehl, November 21, 1994.
c
c-----------------------------------------------------------------------
c
c $Id: trcpth.F,v 1.2 1995/02/17 21:29:04 jhack Exp $
c $Author: jhack $
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
C------------------------------Parameters-------------------------------
c
c $Id: prgrid.h,v 1.2 1995/02/10 01:09:10 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation resolution and I/O parameters
C
      integer
     $     plon,    ! number of longitudes
     $     plev,    ! number of vertical levels
     $     plat,    ! number of latitudes
     $     pcnst,   ! number of constituents (including water vapor)
     $     plevmx,  ! number of subsurface levels
     $     plevp,   ! plev + 1
     $     nxpt,    ! no.of points outside active domain for interpolant
     $     jintmx,  ! number of extra latitudes in polar region
     $     plond,   ! slt extended domain longitude
     $     platd,   ! slt extended domain lat.
     $     plevd    ! fold plev,pcnst indices into one
      integer
     $     plngbuf  ! length of absorptivity/emissivity record
C
      parameter(plon    = 1,
     $          plev    = 18,
     $          plat    = 1,
     $          pcnst   = 1,
     $          plevmx  = 4,
     $          plevp   = plev + 1,
     $          nxpt    = 1,
     $          jintmx  = 1,
     $          plond   = plon + 1 + 2*nxpt,
     $          platd   = plat + 2*nxpt + 2*jintmx,
     $          plevd   = plev*(3 + pcnst))
      parameter(plngbuf = 512*((plond*plevp*plevp + plond*plev*4 +
     $                          plond*plevp)/512 + 1))
C
C------------------------------Commons----------------------------------
c
c $Id: crdcon.h,v 1.1.1.1 1995/02/09 23:26:46 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation constants
C
      common/crdcon/gravit  ,rga     ,cpair   ,epsilo  ,sslp    ,
     $              stebol  ,rgsslp  ,co2vmr  ,dpfo3   ,dpfco2  ,
     $              dayspy  ,pie
C
      real gravit,    ! Acceleration of gravity
     $     rga,       ! 1./gravit
     $     cpair,     ! Specific heat of dry air
     $     epsilo,    ! Ratio of mol. wght of H2O to dry air
     $     sslp,      ! Standard sea-level pressure
     $     stebol,    ! Stefan-Boltzmann's constant
     $     rgsslp,    ! 0.5/(gravit*sslp)
     $     co2vmr,    ! CO2 volume mixing ratio
     $     dpfo3,     ! Voigt correction factor for O3
     $     dpfco2,    ! Voigt correction factor for CO2
     $     dayspy,    ! Number of days per 1 year
     $     pie        ! 3.14.....
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real tnm(plond,plev),     ! Model level temperatures
     $     pnm(plond,plevp),    ! Pressure at model interfaces (dynes/cm2)
     $     qnm(plond,plev),     ! h2o specific humidity
     $     cfc11(plond,plev),   ! CFC11 mass mixing ratio
     $     cfc12(plond,plev),   ! CFC12 mass mixing ratio
     $     n2o(plond,plev),     ! N2O mass mixing ratio
     $     ch4(plond,plev)      ! CH4 mass mixing ratio
C
C Output arguments
C
      real ucfc11(plond,plevp), ! CFC11 path length
     $     ucfc12(plond,plevp), ! CFC12 path length
     $     un2o0(plond,plevp),  ! N2O path length
     $     un2o1(plond,plevp),  ! N2O path length (hot band)
     $     uch4(plond,plevp),   ! CH4 path length
     $     uco211(plond,plevp), ! CO2 9.4 micron band path length
     $     uco212(plond,plevp), ! CO2 9.4 micron band path length
     $     uco213(plond,plevp), ! CO2 9.4 micron band path length
     $     uco221(plond,plevp), ! CO2 10.4 micron band path length
     $     uco222(plond,plevp), ! CO2 10.4 micron band path length
     $     uco223(plond,plevp), ! CO2 10.4 micron band path length
     $     bn2o0(plond,plevp),  ! pressure factor for n2o
     $     bn2o1(plond,plevp),  ! pressure factor for n2o
     $     bch4(plond,plevp),   ! pressure factor for ch4
     $     uptype(plond,plevp)  ! p-type continuum path length
C
C---------------------------Local variables-----------------------------
C
      integer   i,              ! Longitude index
     $          k               ! Level index
      real co2fac(plond,1),     ! co2 factor
     $     alpha1(plond),       ! stimulated emission term
     $     alpha2(plond),       ! stimulated emission term
     $     rt(plond),           ! reciprocal of local temperature
     $     rsqrt(plond),        ! reciprocal of sqrt of temp
     $     pbar(plond),         ! mean pressure
     $     dpnm(plond),         ! difference in pressure
     $     co2mmr
      real diff                 ! diffusivity factor
      data diff /1.66/
c-----------------------------------------------------------------------
c  Calculate path lengths for the trace gases
c-----------------------------------------------------------------------
      co2mmr = 1.51913 * co2vmr
      do i = 1,plon
         ucfc11(i,1) = 1.8 * cfc11(i,1) * pnm(i,1) * rga
         ucfc12(i,1) = 1.8 * cfc12(i,1) * pnm(i,1) * rga
         un2o0(i,1) = diff * 1.02346e5 * n2o(i,1) * pnm(i,1) * rga 
     $                       / sqrt(tnm(i,1))
         un2o1(i,1) = diff * 2.01909 * un2o0(i,1) * 
     $                    exp(-847.36/tnm(i,1))
         uch4(i,1) = diff * 8.60957e4 * ch4(i,1) * pnm(i,1) * rga 
     $                       / sqrt(tnm(i,1))
         co2fac(i,1) = diff * co2mmr * pnm(i,1) * rga
         alpha1(i) = (1.0 - exp(-1540.0/tnm(i,1)))**3.0/sqrt(tnm(i,1))
         alpha2(i) = (1.0 - exp(-1360.0/tnm(i,1)))**3.0/sqrt(tnm(i,1))
         uco211(i,1) = 3.42217e3 * co2fac(i,1) * alpha1(i) * 
     $                             exp(-1849.7/tnm(i,1))
         uco212(i,1) = 6.02454e3 * co2fac(i,1) * alpha1(i) * 
     $                             exp(-2782.1/tnm(i,1))
         uco213(i,1) = 5.53143e3 * co2fac(i,1) * alpha1(i) * 
     $                             exp(-3723.2/tnm(i,1))
         uco221(i,1) = 3.88984e3 * co2fac(i,1) * alpha2(i) * 
     $                             exp(-1997.6/tnm(i,1))
         uco222(i,1) = 3.67108e3 * co2fac(i,1) * alpha2(i) * 
     $                             exp(-3843.8/tnm(i,1))
         uco223(i,1) = 6.50642e3 * co2fac(i,1) * alpha2(i) * 
     $                             exp(-2989.7/tnm(i,1))
         bn2o0(i,1) = diff * 19.399 * pnm(i,1)**2.0 * n2o(i,1) * 
     $                  1.02346e5 * rga / (sslp*tnm(i,1))
         bn2o1(i,1) = bn2o0(i,1) * exp(-847.36/tnm(i,1)) *
     $                  2.06646e5
         bch4(i,1) = diff * 2.94449 * ch4(i,1) * pnm(i,1)**2.0 * rga *
     $                  8.60957e4 / (sslp*tnm(i,1))
         uptype(i,1) = diff * qnm(i,1) * pnm(i,1)**2.0 *
     $                   exp(1800.0*(1.0/tnm(i,1) - 1.0/296.0)) *
     $                   rga / sslp
      end do
      do k = 1,plev
         do i = 1,plon
            rt(i) = 1./tnm(i,k)
            rsqrt(i) = sqrt(rt(i))
            pbar(i) = 0.5 * (pnm(i,k+1) + pnm(i,k)) / sslp
            dpnm(i) = (pnm(i,k+1) - pnm(i,k)) * rga
            alpha1(i) = diff * rsqrt(i) * 
     $                         (1.0 - exp(-1540.0/tnm(i,k)))**3.0
            alpha2(i) = diff * rsqrt(i) * 
     $                         (1.0 - exp(-1360.0/tnm(i,k)))**3.0
            ucfc11(i,k+1) = ucfc11(i,k) +  1.8 * cfc11(i,k) * dpnm(i)
            ucfc12(i,k+1) = ucfc12(i,k) +  1.8 * cfc12(i,k) * dpnm(i)
            un2o0(i,k+1) = un2o0(i,k) + diff * 1.02346e5 * 
     $                                  n2o(i,k) * rsqrt(i) * dpnm(i)
            un2o1(i,k+1) = un2o1(i,k) + diff * 2.06646e5 * n2o(i,k) *
     $           rsqrt(i) * exp(-847.36/tnm(i,k)) * dpnm(i)
            uch4(i,k+1) = uch4(i,k) + diff * 8.60957e4 * ch4(i,k) * 
     $           rsqrt(i) * dpnm(i)
            uco211(i,k+1) = uco211(i,k) + 1.15*3.42217e3 * alpha1(i) *
     $            co2mmr * exp(-1849.7/tnm(i,k)) * dpnm(i)
            uco212(i,k+1) = uco212(i,k) + 1.15*6.02454e3 * alpha1(i) *
     $            co2mmr * exp(-2782.1/tnm(i,k)) * dpnm(i)
            uco213(i,k+1) = uco213(i,k) + 1.15*5.53143e3 * alpha1(i) *
     $            co2mmr * exp(-3723.2/tnm(i,k)) * dpnm(i)
            uco221(i,k+1) = uco221(i,k) + 1.15*3.88984e3 * alpha2(i) *
     $            co2mmr * exp(-1997.6/tnm(i,k)) * dpnm(i)
            uco222(i,k+1) = uco222(i,k) + 1.15*3.67108e3 * alpha2(i) *
     $            co2mmr * exp(-3843.8/tnm(i,k)) * dpnm(i)
            uco223(i,k+1) = uco223(i,k) + 1.15*6.50642e3 * alpha2(i) *
     $            co2mmr * exp(-2989.7/tnm(i,k)) * dpnm(i)
            bn2o0(i,k+1) = bn2o0(i,k) + diff * 19.399 * pbar(i) * rt(i)
     $          * 1.02346e5 * n2o(i,k) * dpnm(i)
            bn2o1(i,k+1) = bn2o1(i,k) + diff * 19.399 * pbar(i) * rt(i) 
     $          * 2.06646e5 * exp(-847.36/tnm(i,k)) * n2o(i,k)*dpnm(i)
            bch4(i,k+1) = bch4(i,k) + diff * 2.94449 * rt(i) * pbar(i)
     $            * 8.60957e4 * ch4(i,k) * dpnm(i)
            uptype(i,k+1) = uptype(i,k) + diff *qnm(i,k) * 
     $                   exp(1800.0*(1.0/tnm(i,k) - 1.0/296.0)) *
     $                   pbar(i) * dpnm(i)
         end do
       end do
       return
       end
      subroutine zenith(calday  ,dodiavg ,clat    ,coszrs  )
C-----------------------------------------------------------------------
C
C Compute cosine of solar zenith angle for albedo and radiation computations
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Rosinski, May, 1994
C
C-----------------------------------------------------------------------
c
c $Id: zenith.F,v 1.1.1.1 1995/02/09 23:27:16 ccm2 Exp $
c $Author: ccm2 $
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
C------------------------------Parameters-------------------------------
c
c $Id: pmgrid.h,v 1.2 1995/02/10 01:09:06 ccm2 Exp $
c $Author: ccm2 $
c
C
C Basic grid point resolution parameters
C
      integer
     $     plon,    ! number of longitudes
     $     plev,    ! number of vertical levels
     $     plat,    ! number of latitudes
     $     pcnst,   ! number of constituents (including water vapor)
     $     plevmx,  ! number of subsurface levels
     $     plevp,   ! plev + 1
     $     nxpt,    ! no.of points outside active domain for interpolant
     $     jintmx,  ! number of extra latitudes in polar region
     $     plond,   ! slt extended domain longitude
     $     platd,   ! slt extended domain lat.
     $     plevd    ! fold plev,pcnst indices into one
C
      parameter(plon   = 1,
     $          plev   = 18,
     $          plat   = 1,
     $          pcnst  = 1,
     $          plevmx = 4,
     $          plevp  = plev + 1,
     $          nxpt   = 1,
     $          jintmx = 1,
     $          plond  = plon + 1 + 2*nxpt,
     $          platd  = plat + 2*nxpt + 2*jintmx,
     $          plevd  = plev*(3 + pcnst))
C
C------------------------------Commons----------------------------------
c
c $Id: crdcon.h,v 1.1.1.1 1995/02/09 23:26:46 ccm2 Exp $
c $Author: ccm2 $
c
C
C Radiation constants
C
      common/crdcon/gravit  ,rga     ,cpair   ,epsilo  ,sslp    ,
     $              stebol  ,rgsslp  ,co2vmr  ,dpfo3   ,dpfco2  ,
     $              dayspy  ,pie
C
      real gravit,    ! Acceleration of gravity
     $     rga,       ! 1./gravit
     $     cpair,     ! Specific heat of dry air
     $     epsilo,    ! Ratio of mol. wght of H2O to dry air
     $     sslp,      ! Standard sea-level pressure
     $     stebol,    ! Stefan-Boltzmann's constant
     $     rgsslp,    ! 0.5/(gravit*sslp)
     $     co2vmr,    ! CO2 volume mixing ratio
     $     dpfo3,     ! Voigt correction factor for O3
     $     dpfco2,    ! Voigt correction factor for CO2
     $     dayspy,    ! Number of days per 1 year
     $     pie        ! 3.14.....
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real calday              ! Calendar day, including fraction
      logical dodiavg          ! true => do diurnal averaging
      real clat                ! Current latitude (radians)
C
C Output arguments
C
      real coszrs(plond)       ! Cosine solar zenith angle
C
C---------------------------Local variables-----------------------------
C
      integer i     ! Longitude loop index
      real phi,     ! Greenwich calendar day + local time + long offset
     $     theta,   ! Earth orbit seasonal angle in radians
     $     delta,   ! Solar declination angle  in radians
     $     sinc,    ! Sine   of latitude
     $     cosc,    ! Cosine of latitude
     $     sind,    ! Sine   of declination
     $     cosd     ! Cosine of declination
      real frac,    ! Daylight fraction
     $     arg,     ! ?
     $     tsun,    ! temporary term in diurnal averaging
     $     coszrsu  ! uniform cosine zenith solar angle 
C
C-----------------------------------------------------------------------
C
C Compute solar distance factor and cosine solar zenith angle usi
C day value where a round day (such as 213.0) refers to 0z at
C Greenwich longitude.
C
C Use formulas from Paltridge, G.W. and C.M.R. Platt 1976: Radiative
C Processes in Meterology and Climatology, Elsevier Scientific
C Publishing Company, New York  p. 57, p. 62,63.
C
      theta = 2.*pie*calday/dayspy
C
C Solar declination in radians:
C
      delta = .006918 - .399912*cos(theta) + .070257*sin(theta) -
     $        .006758*cos(2.*theta) + .000907*sin(2.*theta) -
     $        .002697*cos(3.*theta) + .001480*sin(3.*theta)
C
C Compute local cosine solar zenith angle,
C
      sinc = sin(clat)
      sind = sin(delta)
      cosc = cos(clat)
      cosd = cos(delta)
C
C If using diurnal averaging, then compute the average local cosine solar 
C zenith angle using formulas from paltridge and platt 1976  p. 57, p. 62,63.
C
      if (dodiavg) then
         arg = -(sinc/cosc)*(sind/cosd)
         if (arg .lt. -1.) then
            frac = 1.0
         else if (arg .gt. 1.) then
            frac = 0.0
         else
            frac = (1./pie)*acos(arg)
         endif
         tsun = pie*frac
         if (tsun .gt. 0.) then
            coszrsu =  sinc*sind + (cosc*cosd*sin(tsun))/tsun
         else
            coszrsu = 0.0
         endif
         do i=1,plon
            coszrs(i) = coszrsu
         end do
      else                       ! No diurnal averaging
C
C Calday is the calender day for Greenwich, including fraction
C of day; the fraction of the day represents a local time at
C Greenwich; to adjust this to produce a true instantaneous time
C For other longitudes, we must correct for the local time change:
C local time based on the longitude and day of year
C then compute the local cosine solar zenith angle
C
         do i=1,plon
            phi       = calday + (real(i-1)/real(plon))
            coszrs(i) = sinc*sind - cosc*cosd*cos(2.*pie*phi)
         end do
      end if
C
      return
      end





      integer function intmax(n,ix,inc)
c
c $Id: intmax.F,v 1.1.1.1 1995/02/09 23:26:15 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
      integer n,inc
      integer ix(*)
c
      integer i,mx
c
      mx = ix(1)
      intmax = 1
      do i=1+inc,inc*n,inc
         if (ix(i).gt.mx) then
            mx = ix(i)
            intmax = i
         end if
      end do
      return
      end
      integer function isrchfgt(n, array, inc, target)
c
c $Id: isrchfgt.F,v 1.1.1.1 1995/02/09 23:26:15 ccm2 Exp $
c $Author: ccm2 $
c
      integer n,inc
      real array(*),target
      integer i
      if (n.le.0) then
         isrchfgt = 0
         return
      end if
      ind = 1
      do i=1,n
         if (array(ind).gt.target) then
            isrchfgt = i
            return
         else
            ind = ind + inc
         end if
      end do
      isrchfgt = n + 1
      return
      end

      integer function isrchfle(n, array, inc, target)
c
c $Id: isrchfle.F,v 1.1.1.1 1995/02/09 23:26:15 ccm2 Exp $
c $Author: ccm2 $
c
      integer n,inc
      real array(*),target
      integer i
      if (n.le.0) then
         isrchfle = 0
         return
      end if
      ind = 1
      do i=1,n
         if (array(ind).le.target) then
            isrchfle = i
            return
         else
            ind = ind + inc
         end if
      end do
      isrchfle = n + 1
      return
      end

      integer function myhandler( sig, code, context )
      integer sig, code, context(5)
c 061003 BPB      myhandler = abort()
      return
      end
      subroutine wheneq(n,array,inc,target,index,nval)
c
c $Id: wheneq.F,v 1.2 1995/02/28 22:19:52 ccmproc2 Exp $
c $Author: ccmproc2 $
c
      integer index(*),
     $        array(*),
     $        target
      ina=1
      nval=0
      if(inc .lt. 0) ina=(-inc)*(n-1)+1
      do i=1,n
         if(array(ina) .eq. target) then
           nval=nval+1
           index(nval)=i
         end if
         ina=ina+inc
      enddo
      return
      end
c----------
      subroutine whenfgt(n,array,inc,target,index,nval)
c
c $Id: whenfgt.F,v 1.1.1.1 1995/02/09 23:26:17 ccm2 Exp $
c $Author: ccm2 $
c
	dimension index(*), array(*)
	ina=1
	nval=0
	if(inc .lt. 0) ina=(-inc)*(n-1)+1
	do 100 i=1,n
	    if(array(ina) .gt. target) then
	    nval=nval+1
	    index(nval)=i
	    end if
	    ina=ina+inc
 100    continue
      return
      end           
c---------
      subroutine whenflt(n,array,inc,target,index,nval)
c
c $Id: whenflt.F,v 1.1.1.1 1995/02/09 23:26:17 ccm2 Exp $
c $Author: ccm2 $
c
	dimension index(*), array(*)
	ina=1
	nval=0
	if(inc .lt. 0) ina=(-inc)*(n-1)+1
	do 100 i=1,n
	    if(array(ina) .lt. target) then
	    nval=nval+1
	    index(nval)=i
	    end if
	    ina=ina+inc
 100    continue
      return
      end           
      subroutine whenne(n,array,inc,target,index,nval)
c
c $Id: whenne.F,v 1.1.1.1 1995/02/09 23:26:18 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
      integer n,inc
      integer array(*),target
      integer index(*),nval
c
      integer i,ina
c
      ina=1
      nval=0
      if(inc.lt.0)ina=(-inc)*(n-1)+1
      do i=1,n
         if (array(ina).ne.target) then
            nval=nval+1
            index(nval)=i
         end if
         ina = ina+inc
      end do
      return
      end
