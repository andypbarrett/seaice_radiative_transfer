      subroutine calling_sub()

      implicit none

      integer n
      real myreal
      real myarray(5)

      common /mycom/ n, myreal, myarray
      
      call print_mycom()

      return
      end


      subroutine print_mycom()
      implicit none

      integer n
      real myreal
      real myarray(5)

      common /mycom/ n, myreal, myarray
      
      print *, 'integer: ', n
      print *, 'myreal: ', myreal
      print *, 'myarray: ', myarray

      return
      end
      
