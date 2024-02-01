      subroutine splin1 ( y,f,n,im,il,cof,tension )
!
!     1-d cubic spline interpolation assuming periodic b.c.
!
!     y(n)    : input array to be interpolated from
!     f(im)   : output array to be interpolated to
!     n       : dimension of input array
!     im      : dimension of output array
!     il(im)  : address index of output points by input grid units
!     cof(n*3): grid-related coef. for gaussain elimination
!     tension : =1.0 for full cubic,  =0.0 for linear
!   * note    : arrays il, and cof are generated in splinc
!
      implicit  none
      integer   n,im,i,ii,iip,il(im)

      real      y(n),f(im),cof(n*3),tension
!
!     local work arrays:
      real      d2ydx2(2000), yp1(2000),dx,dm,tem
      real      tem1,tem2,tem13,tem23
!
      dx = 1.0/float(n)
      dm = 1.0/float(im)
!
      tem = tension / (dx*dx)
      do 95 i = 1, n-2
      f(i+1) = ( y(i) - 2.0*y(i+1) + y(i+2) ) * tem
  95  continue
      f(1) = ( y(n) - 2.0*y(1) + y(2) ) * tem
      f(n) = ( y(n-1) - 2.0*y(n) + y(1) ) * tem
!
      call gauslv ( n,f,d2ydx2,cof )
!
      d2ydx2(n+1) = d2ydx2(1)
      do 97 i =1, n
      yp1(i) = y(i)
  97  continue
      yp1(n+1) = y(1)
!
!cc   if (derv )  go to 40
      do 30 i = 1, im
      ii = il(i)
      iip = ii + 1
      tem1 = dx*ii - dm*(i-1.0)
      tem2 = dm*(i-1.0) - dx*(ii-1.0)
      tem13 = tem1**3/dx
      tem23 = tem2**3/dx
      f(i) = d2ydx2(ii)*tem13 + d2ydx2(iip)*tem23  &
           + tem1*(yp1(ii)/dx - dx*d2ydx2(ii))     &
           + tem2*(yp1(iip)/dx- dx*d2ydx2(iip))
  30  continue
!c
!cc   return
!c40  continue
!
!cc   to output first derivative of the input array
!cc
!cc   do 50 i = 1, im
!cc   ii = il(i)
!cc   iip = ii + 1
!cc   tem1 = dx*ii - dm*(i-1.0)
!cc   tem2 = dm*(i-1.0) - dx*(ii-1.0)
!cc   tem12 = tem1*tem1/dx * 3.0
!cc   tem22 = tem2*tem2/dx * 3.0
!cc   f(i) = d2ydx2(iip)*tem22 - d2ydx2(ii)*tem12
!cc  1     + (yp1(iip)-yp1(ii))/dx - (d2ydx2(iip)-d2ydx2(ii))*dx
!c50  continue
!
      return
      end

      subroutine splin1_sp ( y,fout,n,im,il,cof,tension )
!
!     1-d cubic spline interpolation assuming periodic b.c.
!
!     y(n)    : input array to be interpolated from
!     f(im)   : output array to be interpolated to
!     n       : dimension of input array
!     im      : dimension of output array
!     il(im)  : address index of output points by input grid units
!     cof(n*3): grid-related coef. for gaussain elimination
!     tension : =1.0 for full cubic,  =0.0 for linear
!   * note    : arrays il, and cof are generated in splinc
!
      implicit  none
      integer   n,im,i,ii,iip,il(im)

      real      y(n),f(im),cof(n*3),tension
      real*4 fout(im)
!
!     local work arrays:
      real      d2ydx2(2000), yp1(2000),dx,dm,tem
      real      tem1,tem2,tem13,tem23
!
      dx = 1.0/float(n)
      dm = 1.0/float(im)
!
      tem = tension / (dx*dx)
      do 95 i = 1, n-2
      f(i+1) = ( y(i) - 2.0*y(i+1) + y(i+2) ) * tem
  95  continue
      f(1) = ( y(n) - 2.0*y(1) + y(2) ) * tem
      f(n) = ( y(n-1) - 2.0*y(n) + y(1) ) * tem
!
      call gauslv ( n,f,d2ydx2,cof )
!
      d2ydx2(n+1) = d2ydx2(1)
      do 97 i =1, n
      yp1(i) = y(i)
  97  continue
      yp1(n+1) = y(1)
!
!cc   if (derv )  go to 40
      do 30 i = 1, im
      ii = il(i)
      iip = ii + 1
      tem1 = dx*ii - dm*(i-1.0)
      tem2 = dm*(i-1.0) - dx*(ii-1.0)
      tem13 = tem1**3/dx
      tem23 = tem2**3/dx
      fout(i) = d2ydx2(ii)*tem13 + d2ydx2(iip)*tem23  &
           + tem1*(yp1(ii)/dx - dx*d2ydx2(ii))     &
           + tem2*(yp1(iip)/dx- dx*d2ydx2(iip))
  30  continue
      return
      end
