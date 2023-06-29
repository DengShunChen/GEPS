      subroutine vterp (nx,my,lmaxp,lev,yr,f,yin,dout,tensy)
!
!          a bicubic spline interpolator to interpolate from a grid
!          with constant i (first dimension) grid spacing and variable
!          j (second dimension) grid spacing t0 a grid with
!          variable grid spacing. all grids are assumed to have point
!          (i,1) in the lower left corner with i increasing to the right
!          and j increasing upward.
!
! **** input ****
!
!  nx: no. of points in e-w direction of input arrays
!  my: no. of points in n-s direction of input arrays
!  lmaxp: no. of levels in input arrays
!  lev: no. of levels in output arrays
!  yr: independent interpolation variable for input grid
!  f: dependent variable for input grid
!  yin: independent interpolation variable for output grid
!  tensy: cubic spline tension factor.
!
! **** output ****
!
!  dout: dependent variable on output grid
!
      use const, only: RTYPE
!
      implicit  none

      integer   nx,my,lmaxp,lev

      real      f(nx,lmaxp,my),dout(nx*lev,my),yr(nx,lmaxp,my) &
      , yin(nx*lev,my),tensy(lmaxp),tp1(nx*lev*4)
!
!sun  include '../include/paramt.h' .. change im,lm to nx,lev
!
      real(kind=RTYPE) fxx(nx,lmaxp),fyy(nx,lmaxp),pjy(nx*lev*4)


      integer   ipt(nx*lev)

      integer   jym2,nxjym2,mn,j,ll,i,k,i1,i2,i3
!
!          compute ipt and pjy
!
      jym2=lmaxp-2
      nxjym2=nx*jym2
      mn= nx*lev
      do 200 j=1,my
!
      call setupv(yr(1,1,j),yin(1,j),mn,nx,lmaxp,pjy,ipt,tp1)
!
!          compute fyy
!
      ll= nxjym2+nx
      do 110 i=1,ll
      fxx(i,2)= yr(i,2,j)-yr(i,1,j)
  110 continue
!
      do 210 i=1,nxjym2
      fyy(i,2)= (fxx(i,3)*(f(i,1,j)-f(i,2,j))+fxx(i,2)*(f(i,3,j) &
       -f(i,2,j)))/(fxx(i,3)*fxx(i,3)*fxx(i,2))
      fxx(i,2)= fxx(i,2)/fxx(i,3)
  210 continue
!
!      call trdivv(nx,jym2,fxx(1,2),tp1,fyy(1,2))
      call trdivv(nx,jym2,fxx(1,2),fyy(1,2))
!
      do 100 i=1,nx
      fyy(i,1)= 0.0
      fyy(i,lmaxp)=0.0
  100 continue
!
!  apply tension
!
      do 5 k=1,lmaxp
      do 5 i=1,nx
      fyy(i,k)= fyy(i,k)*tensy(k)
    5 continue
!
      call gathv(mn,nx,lev,ipt,fyy,f(1,1,j),tp1)
!
      i1= mn
      i2= mn*2
      i3= mn*3
      do 130 i=1,mn
      dout(i,j)=tp1(i)*pjy(i)+tp1(i+i1)*pjy(i+i1)+tp1(i+i2)*pjy(i+i2) &
       +tp1(i+i3)*pjy(i+i3)
  130 continue
  200 continue
      return
      end
