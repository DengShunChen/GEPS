      subroutine zilch (x,m)
!
!  subroutine to zero out arrays
!
!  ***input***
!
!  x: input array to zilch
!  m: number of elements to zilch
!
!  ***output***
!
!  x: zeroed array
!
! ****************************************************************
!
      implicit none

      integer m,i

      real x(m)
      do 1 i=1,m
      x(i)= 0.0
    1 continue
      return
      end
