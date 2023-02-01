      subroutine splinc ( n,im,il,ib,cof )
!
!     prepare constant coefficients for 1-d cubic spline interpolation
!     ( warning: for spherical grids only )
!
!     n       : dimension of input array to be interpolated from
!     im      : dimension of output array to be interpolated to
!     il(im)  : address index of output points by input grid units
!     ib(n)   : address index of input points by output grid units
!     cof(n*3): cubic interpolation coef. for each grid interval
!
      implicit  none
      integer   n,im,i,j,nm2,nt2,ne,nu
      real      fac,facx,v,bn

      integer   il(im),ib(n)
      real      cof(n*3)
!
      fac = float(n) / float(im)
      do 100 i = 1, im
        il(i) = 1.0 + (float(i-1)*fac + 0.00001)
  100 continue
      facx = float(im) / float(n)
      do 110 i = 1, n
        ib(i) = 1.0 + (float(i-1)*facx + 0.00001)
  110 continue
!
      nm2 = n - 2
      nt2 = n*2
!
!     prepare grid-related constant coef. used in gaussian elimination
!
      cof(n+1) = 0.25
      v = 1.0
      cof(1) = cof(n+1)
      cof(nt2+1) = cof(n+1)
      bn = - v*cof(nt2+1) + 4.0
      do 200 j = 2, nm2
         ne = j + n
         cof(ne) = 1.0/(4.0-cof(j-1))
         cof(j) = cof(ne)
         nu = j + nt2
         cof(nu) = - cof(nu-1)*cof(ne)
         v = - v*cof(j-1)
         bn = bn - v*cof(nu)
  200 continue
      v = 1.0 - v*cof(nm2)
      ne = nt2
      cof(ne-1) = 1.0/(4.0-cof(nm2))
      nu = nt2 + n
      cof(n-1) = (1.0-cof(nu-2)) * cof(ne-1)
      cof(ne) = 1.0/(bn - v*cof(n-1))
!
      return
      end
