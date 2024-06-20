      subroutine syslbl_r (clrec,idtg,itau,ciflap)
!
      use const, only: ihdgi
      implicit none
      character*4 ciflap
      character*6 clrec
      integer   itau
      integer*8 idtg
!
      write(ihdgi,1) clrec,itau,ciflap,idtg
#ifdef I38K
 1    format(a6,i6.6,a4,i12.12)
#else
 1    format(a6,i4.4,a4,i12.12)
#endif
      return
      end
