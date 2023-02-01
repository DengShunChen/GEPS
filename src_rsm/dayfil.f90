      subroutine dayfil (msg)
!
      use rank

      implicit none

      character*120 msg
!
!     call remark(msg)
      if(myrank .eq. 0) print 100,msg
  100 format(1x,a120)
!
      return
      end
