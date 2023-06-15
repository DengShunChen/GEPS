      subroutine gathv (mn,nx,lev,ipt,fyy,f,tp1)
!
!  sorting routine to find coefficients for cubic splines for each
!  interpolation grid element.  see bicubv for parameter descriptions
!
      use const, only: RTYPE
!
      implicit  none

      integer   mn,nx,lev,ipt(mn)
      real(kind=RTYPE) fyy(nx*lev),f(nx*lev)
      real      tp1(mn,4)
      integer   i,inx
!
      do 20 i=1,mn
      tp1(i,1)= fyy(ipt(i))
      tp1(i,3)= f  (ipt(i))
      inx= ipt(i)-nx
      tp1(i,2)= fyy(inx)
      tp1(i,4)= f  (inx)
   20 continue
!
      return
      end
