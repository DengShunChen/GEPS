      subroutine gausl3 (n,xa,xb,wt,ab)
!
! weights and abscissas for nth order gaussian quadrature on (xa,xb).
! input arguments
!
! n  -the order desired
! xa -the left endpoint of the interval of integration
! xb -the right endpoint of the interval of integration
! output arguments
! ab -the n calculated abscissas
! wt -the n calculated weights
!
! modify to f90 by C-H Lee and sort by River Chen in 2015
!
      use rank
      use const, only: RTYPE

      implicit double precision (a-h,o-z)
!
      real  xa,xb
      real(kind=RTYPE) ab(n),wt(n)
!
! machine dependent constants---
!  tol - convergence criterion for double precision iteration
!  pi  - given to 15 significant digits
!  c1  -  1/8                     these are coefficients in mcmahon"s
!  c2  -  -31/(2*3*8**2)          expansions of the kth zero of the
!  c3  -  3779/(2*3*5*8**3)       bessel function j0(x) (cf. abramowitz,
!  c4  -  -6277237/(3*5*7*8**5)   handbook of mathematical functions).
!  u   -  (1-(2/pi)**2)/4
!
!      data tol/1.d-14/,pi/3.14159265358979/,u/.148678816357662/
      data tol/1.d-14/
      data c1,c2,c3,c4/.125,-.080729166666667,.246028645833333, &
                      -1.82443876720609 /
!
! maximum number of iterations before giving up on convergence
!
      data maxit /5/
!
! arithmetic statement function for converting integer to double
!
      dbli(i) = dble(float(i))
!
      pi = 4.*atan(1.)
      u  = (1-(2/pi)**2)/4.
      ddif = .5d0*(dble(xb)-dble(xa))
      dsum = .5d0*(dble(xb)+dble(xa))
      if (n .gt. 1) go to 101
      ab(1) = 0.
      wt(1) = 2.*ddif
      go to 107
  101 continue
      nnp1 = n*(n+1)
      cond = 1./sqrt((.5+float(n))**2+u)
      lim = n/2
!
      do 105 k=1,lim
	 b = (float(k)-.25)*pi
	 bisq = 1./(b*b)
!
! rootbf approximates the kth zero of the bessel function j0(x)
!
	 rootbf = b*(1.+bisq*(c1+bisq*(c2+bisq*(c3+bisq*c4))))
!
!      initial guess for kth root of legendre poly p-sub-n(x)
!
	 dzero = cos(rootbf*cond)
	 do 103 i=1,maxit
!
	    dpm2 = 1.d0
	    dpm1 = dzero
!
!       recursion relation for legendre polynomials
!
	    do 102 nn=2,n
		dp = (dbli(2*nn-1)*dzero*dpm1-dbli(nn-1)*dpm2)/dbli(nn)
		dpm2 = dpm1
		dpm1 = dp
  102       continue
	    dtmp = 1.d0/(1.d0-dzero*dzero)
	    dppr = dbli(n)*(dpm2-dzero*dp)*dtmp
	    dp2pri = (2.d0*dzero*dppr-dbli(nnp1)*dp)*dtmp
	    drat = dp/dppr
!
!       cubically-convergent iterative improvement of root
!
	    dzeri = dzero-drat*(1.d0+drat*dp2pri/(2.d0*dppr))
	    ddum= dabs(dzeri-dzero)
	 if (ddum .le. tol) go to 104
	    dzero = dzeri
  103    continue
!
	 if(myrank .eq. 0) print 504
!
  504    format(1x,' in gausl3, convergence failed')
  104    continue
	 ddifx = ddif*dzero
	 ab(k) = dsum-ddifx
	 wt(k) = 2.d0*(1.d0-dzero*dzero)/(dbli(n)*dpm2)**2*ddif
	 i = n-k+1
	 ab(i) = dsum+ddifx
	 wt(i) = wt(k)
  105 continue
!
      if (mod(n,2) .eq. 0) go to 107
      ab(lim+1) = dsum
      nm1 = n-1
      dprod = n
      do 106 k=1,nm1,2
	 dprod = dbli(nm1-k)*dprod/dbli(n-k)
  106 continue
      wt(lim+1) = 2.d0/dprod**2*ddif
  107 return
      end
