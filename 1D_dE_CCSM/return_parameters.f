C----------------------------------------------------------------------
C     Returns dummy parameters for testing
C
C
C----------------------------------------------------------------------
      subroutine return_parameters()

      implicit none

      integer
     $     plon,
     $     plev,
     $     plevp,
     $     plond,   ! slt extended domain longitude
     $     nxpt,    ! no.of points outside active domain for interpolant
     $     i
      
      parameter(plon = 1,
     $          plev = 18,
     $          nxpt = 1,
     $          plevp = plev + 1,
     $          plond = plon + 1 + 2*nxpt)
      
      real asdir(plond),    ! albedo: shortwave, direct
     $     asdif(plond),    ! albedo: shortwave, diffuse
     $     aldir(plond),    ! albedo: longwave, direct
     $     aldif(plond),    ! albedo: longwave, diffuse
     $     F_SW_vs,         ! solar vs absorbed in sea ice
     $     F_SW_ni,         ! solar ni absorbed in sea ice
     $     F_SW_srf_vs,     ! vs solar absorbed in sea ice surface layer
     $     F_SW_srf_ni,     ! ni solar absorbed in sea ice surface layer
     $     F_SW_ocn_vs,     ! vs solar absorbed in underlying ocean
     $     F_SW_ocn_ni,     ! ni solar absorbed in underlying ocean
     $     F_SW_srf,        ! total solar absorbed in sea ice surface layer
     $     Q_SW_vs,         ! solar vs absorbed in sea ice layer
     $     Q_SW_ni          ! solar ni absorbed in sea ice layer

      common /output/
     $     asdir, asdif, aldir, aldif,
     $     F_SW_ocn_vs, F_SW_ocn_ni

      i = 1
      
      asdir(i) = 1.
      asdif(i) = 2.
      aldir(i) = 3.
      aldif(i) = 4.
      F_SW_ocn_vs = 5.
      F_SW_ocn_ni = 6.
      
      return
      end
