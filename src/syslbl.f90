      subroutine syslbl (clrec,idtg,itau,ciflap,cirec)
!
      implicit none
      character*4 ciflap
      character*6 clrec
      character*28 cirec
      integer   itau
      integer*8 idtg
!
      write(cirec,1) clrec,itau,ciflap,idtg
 1    format(a6,i6.6,a4,i12.12)
      return
      end
