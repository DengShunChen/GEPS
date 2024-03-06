      subroutine syslbl_w (clrec,idtg,itau,ciflap)
!
      use const, only: ihdgo
      implicit none
      character*4 ciflap
      character*6 clrec
      integer   itau
      integer*8 idtg
!
      write(ihdgo,1) clrec,itau,ciflap,idtg
#ifdef O38K
 1    format(a6,i6.6,a4,i12.12)
#else
 1    format(a6,i4.4,a4,i12.12)
#endif
      return
      end
