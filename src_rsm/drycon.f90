      subroutine drycon (j,indx,th,dp,nx,lev )
!
!     perform dry convective adjustment on a number of potential
!     temperature soundings using a variational method.
!
!   parameters:
!     th(nx,lev)    : potential temp (k),
!     dp(lev)       : sigma thickness of layers,
!     nx          : horizontal dimension
!     lev           : vertical dimension
!
!     by c.s. liou  6/92
!     modify to f90 by C-H Lee and sort by River Chen in 2015
!
      implicit  none
      integer   j,nx,lev
      integer   indx(nx)
      real      th(nx,lev),dp(lev)
!
!     local work arrays
!
      real      work(nx*lev*2)
      real      dth(nx*(lev-1)),b(lev-1),c(lev-1)
!
      logical   stable(nx*(lev-1))

      integer   lm1,lmd,nlv,nlv1,l,i,il,k,ilmax,ismax
      real      wmax

      lm1 = lev - 1
      lmd = lev * 2
      nlv = nx*lev
      nlv1= nx*lm1
!
!   check for stability of each layer
!
      do 100 l = 1, lm1
      do 100 i = 1, nx
      il=i+nx*(l-1)
      dth(il) = (th(i,l+1) - th(i,l))*dp(l)
      stable(il) = dth(il) .le. 0.0
 100  continue
!
!   compute off diagonal coeff
!
      do 200 l = 1, lm1
      c(l)= dp(l)/dp(l+1)
      b(l)= - c(l) - 1.0
  200 continue
!
!   iterate (20 times max.) to adjust potential temp
!
      do 800 k = 1, 20
!
!   force no adjustment for stable layers with zero forcing
!
      do 300 i = 1, nlv1
      if ( stable(i) )  dth(i) = 0.0
  300 continue
!
!   solve a tri-diagonal system for temperature adjustment
!
      call triadj( dth,b,c,nx,lm1,lmd,work,stable )
!
!   adjust potential temperature
!
      do 400 l = 2, lm1
      do 400 i = 1, nx
      il=i+nx*(l-1)
      work(il) = (dth(il) - dth(il-nx))/dp(l)
  400 continue
      do 410 i = 1, nx
      work(i) = dth(i) / dp(1)
      work(i+nlv1) = - dth(i+nlv1-nx) / dp(lev)
  410 continue
      wmax = 0.0
      do 420 i = 1, nlv
      th(i,1) = th(i,1) - work(i)
      wmax = max ( wmax, abs(work(i)) )
  420 continue
!
!     assume convergence if max correction is less than .05 degrees
!
      if( wmax .lt. 0.05 ) go to 900   
!
      do 500 l = 1, lm1
      do 500 i = 1, nx
      il= i+nx*(l-1)
      dth(il) = ( th(i,l+1) - th(i,l) )*dp(l)
      stable(il) = dth(il) .le. 0.0
  500 continue
      ilmax= ismax(nx*lm1,dth,1)
      ilmax= ilmax-nx*((ilmax-1)/nx)
      ilmax= indx(ilmax)
!
  800 continue
      print 910, wmax,ilmax,j
  910 format(1x,'*** warning: no drycon converge , wmax = ', &
                e13.6 ,',at i,j=',2i4)
!
  900 continue
!
      return
      end
