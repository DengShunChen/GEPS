      subroutine gauslv ( n,y,z,cof )
!
!     gaussian elimination solver with grid-related constnats provided
!
      use paramt

      implicit  none
      integer   n
      real      y(n),z(n),cof(n*3)

      integer   nm,nm2,nt2,j,k,nu,ne
      real      v,yn
!
      nm = n - 1
      nm2 = n - 2
      nt2 = n*2
!
      v = 1.0
      z(1) = y(1) * cof(n+1)
      yn = y(n) - v*z(1)
      do 130 j = 2, nm2
      v = -v*cof(j-1)
      z(j) = (y(j) - z(j-1)) * cof(j+n)
      yn = yn - v*z(j)
  130 continue
      v = 1.0 - v*cof(nm2)
      ne = nt2
      z(nm) = (y(nm) - z(nm2)) * cof(ne-1)
      yn = yn - v*z(nm)
!
!     backward substitude
!
      z(n) = yn * cof(ne)
      z(nm) = z(nm) - cof(nm)*z(n)
      do 140 j = 2, nm
         k = n - j
         nu = k + nt2
         z(k) = z(k) - cof(k)*z(k+1) - cof(nu)*z(n)
  140 continue
      return
      end
