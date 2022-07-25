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
     $     nxpt    ! no.of points outside active domain for interpolant

      parameter(plon = 1,
     $          plev = 18,
     $          nxpt = 1,
     $          plevp = plev + 1,
     $          plond = plon + 1 + 2*nxpt)
      
      real asdir(plond),        ! albedo: shortwave, direct
     $     asdif(plond),        ! albedo: shortwave, diffuse
     $     aldir(plond),        ! albedo: longwave, direct
     $     aldif(plond)         ! albedo: longwave, diffuse

      common /output/
     $     asdir, asdif, aldir, aldif
