      subroutine drychk (j,tt,plk,dp,nxj,nx,lev,ncon )
!
!     perform dry convective adjustment check and collect
!     all unsatble points together for the adjustment
!
!     tt(nx,lev)    : real temperature (k),
!     plk(nx,lev)   : pl**kapa at model odd levels,
!     dp(lev)       : sigma thickness of layers,
!     nx           : horizontal dimension
!     lev           : vertical dimension
!
!     by c.s. liou  6/92
!     modify to f90 by C-H Lee and sort by River Chen in 2015
!
      use const, only: RTYPE
!              
      implicit  none

      integer   j,nxj,nx,lev,ncon
      real      dp(lev)
      real(kind=RTYPE) tt(nx,lev),plk(nx,lev)
!
!     local work array
!
      integer   indx(nx)
      real      store(nx,lev)
      logical   unstbl(nx)
!
      integer   nxlev,lm1,i,k,l,il
      real      eps

      eps   = 0.1
      nxlev = nx*lev
      lm1   = lev - 1
!
!     compute potential temperature
!
      do 50 k = 1, lev
      do 50 i = 1, nxj
      store(i,k) = tt(i,k)/plk(i,k)
  50  continue
!
!     detect unstable points
!
      do 100 i = 1, nxj
      unstbl(i) = (store(i,2)-store(i,1)) .ge. eps
  100 continue
      do 200 l = 2, lm1
      do 200 i = 1, nxj
      unstbl(i) = unstbl(i) .or. ( (store(i,l+1)-store(i,l)).ge.eps )
  200 continue
      ncon = 0
      do 300 i = 1, nxj
      if ( unstbl(i) ) then
      ncon = ncon + 1
      indx(ncon)= i
      endif
  300 continue
      if ( ncon .eq. 0 )  return
!
!     pack 'tt/plk' of unstable points into work array 'store'
!
      do 400 l=1,lev
      do 400 i = 1, ncon
         il = (l-1)*ncon + i
         store(il,1) = tt(indx(i),l)/plk(indx(i),l)
  400 continue
!
!     adjust real temperature tt to dry adiabatic profile
!
      call drycon (j,indx,store,dp,ncon,lev )
!
      do 500 l=1,lev
      do 500 i = 1, ncon
         il = (l-1)*ncon + i
         tt(indx(i),l) = store(il,1)*plk(indx(i),l)
  500 continue
!
      return
      end
