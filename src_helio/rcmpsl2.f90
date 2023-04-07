      subroutine rcmpsl2( ns,nx,lev,u0,as,ww1,ww2,ww,ttau,ffu,ffd )
!      
! **********************************************************************
! The generalized two stream approximation for nonhomgeneous atmospheres
! in  the  solar  wavelengths.  
!
! ++ input/output variables:
!    (1) input variables :
!        ns    : number of grid points to calculate sw flux
!        nx    : horizontal dimension of input array
!        lev   : number of vertical layers
!        u0    : cos(zenith angle)
!        as    : ground surface albedo
!        ww1~2 : first two coeff. of legendre polynomials for
!                the scattering phase function
!        ww    : single scattering albedo
!        ttau(lev+1) : total optical depth from top to reference level
!
!    (2) output variables :
!        ffu
!        ffd
!
! ++ auther :
!    original : fu and liou (1992)
!    modified : c. t. fong (1996)
!
! **********************************************************************
!
      implicit  none
      integer   ns,nx,lev
      real      u0(nx),as(nx),ww1(nx*lev),ww2(nx*lev),ww(nx*lev), &
                ttau(nx*(lev+1))
      real      ffu(nx*(lev+1)),ffd(nx*(lev+1))
!
!   local working arrays
!
      real      w1(nx*lev),w(nx*lev),tau(nx*(lev+1))
      real      lamdan(nx*lev),gamman(nx*lev),caddn(nx*lev)       &
               ,cminn(nx*lev),caddn0(nx*lev),cminn0(nx*lev)       &
               ,aa(nx*lev),bb(nx*lev),expn(nx*lev),g1g2n(nx*lev)  &
               ,ru0(nx*lev),u0t(nx*lev),sol(nx*(lev+1))
      real      f(nx*lev),xn(nx*lev),yn(nx*lev),zn(nx*lev)
      real      etri(nx,lev*2,4)
      logical edding,quadra

      integer   nsn,nsn1,lev2,i,j,jk,jj
      integer   j2n,j2n1,jje,jke,j2m1,j2m
      real      sqrt3,rsqrt3,pi,f0,fw,x,y,z,gamma1,gamma2,gamma3,gamma4
      real      ru0_mod,dtau,ssfc,wm1,wm2,xx
!
      edding = .true.
      quadra = .false.
!
      nsn  = ns*lev
      nsn1 = ns*(lev+1)
      lev2 = lev*2
      sqrt3= sqrt(3.)
      rsqrt3= 1./sqrt3
      pi   = 4.*atan(1.)
      f0   = 1./pi
!
!fong01
! incorporate a delta-function adjustment to account for the forward
! diffraction  peak in the context of the two stream approximations.
! w1(n),  w(n), and tau(n+1) are the adjusted parameters.
!  
      do 6 i =1, nsn
      f(i) = ww2(i) / 5.0
    6 continue
!
      do 7 i = 1, ns
      tau(i) = 0.
      sol(i) = u0(i)
    7 continue
!
      do 10 i = 1, nsn
      fw = 1.0 - f(i) * ww(i) 
      w1(i) = ( ww1(i) - 3.0 * f(i) ) / ( 1.0 - f(i) )
      w(i)  = ( 1.0 - f(i) ) * ww(i) / fw
      u0t(i) = (ttau(i+ns) - ttau(i)) * fw
   10 continue
!
      do 12 j = 1, lev
      jk = j*ns
      jj = jk-ns
      do 12 i = 1, ns
      tau(i+jk) = tau(i+jj) + u0t(i+jj)
   12 continue
!fong01
!     
!fong02
!
      do 20 i = 1, nsn
      w(i) = min( w(i), 0.999999 )
   20 continue
!
      do 22 j = 1, lev
      jk = j*ns
      jj = jk-ns
      do 22 i = 1, ns
      u0t(i+jj) = u0(i)
      ru0(i+jj) = 1./u0(i)
   22 continue
!
      do 25 i = 1, nsn
      sol(i+ns) = u0t(i) * exp( -tau(i+ns) * ru0(i) )
   25 continue
!
!
! The following is used to calculate the Coefficients For Generalized
! Two-Stream scheme.  We can make choices between Eddington, quadrature
! and  hemispheric  mean  schemes  through  logical variables 'edding',
! 'quadra', and 'hemisp'.  The  Eddington  and  quadrature  schemes are 
! discussed in detail by Liou (1992).  The  hemispheric  mean scheme is 
! derived by assuming that the phase function is equal to 1 + g in  the 
! forward scattering hemisphere and 1 - g  in  the  backward scattering 
! hemisphere where g is the asymmetry factor.   The hemispheric mean is
! only used for infrared wavelengths (Toon et al. 1989).
!
      do 40 i = 1, nsn
!
!fong02_1
!cc   if ( edding ) then
	   x = 0.25 * w1(i)     ! w1 = 3g
	   y = w(i) * x
           gamma1 = 1.75 - w(i) - y
	   gamma2 = - 0.25 + w(i) - y
           gamma3 = 0.5 - x * u0t(i)
           gamma4 = 1.0 - gamma3
!fong      ugts1 = 0.5
!
!cc   elseif ( quadra ) then
!          x = 0.5*sqrt3 * w(i)
!          y = 0.5*rsqrt3 * w1(i)  ! rsqrt3 = 1/sqrt3
!          z = y * w(i)
!          gamma1 = sqrt3 - x - z
!          gamma2 = x - z
!          gamma3 = 0.5 - y * u0t(i)
!          gamma4 = 1.0 - gamma3
!fong      ugts1 = rsqrt3
!cc   endif
!fong02_1
      lamdan(i) = exp( 0.5*log((gamma1 + gamma2)*(gamma1 - gamma2)) )
!cc   lamdan(i) = sqrt((gamma1 + gamma2)*(gamma1 - gamma2))
      gamman(i) = gamma2 / ( gamma1 + lamdan(i) )
      g1g2n(i) = gamma1 + gamma2
      fw = pi * f0 * w(i) * exp ( - ru0(i) * tau(i) )
      x = exp ( - ru0(i) * ( tau(i+ns)-tau(i) ) )
      z = lamdan(i) * lamdan(i) - ru0(i) * ru0(i)
      if(abs(z).lt.0.000001)then
       ru0_mod=1./(u0t(i)+0.01*u0t(i))
       z = lamdan(i) * lamdan(i) - ru0_mod * ru0_mod
      end if      
      caddn0(i) = fw * ( ( gamma1 - ru0(i) ) * gamma3 +  &
                     gamma4 * gamma2 ) / z
      cminn0(i) = fw * ( ( gamma1 + ru0(i) ) * gamma4 +  &
                     gamma3 * gamma2 ) / z
      caddn(i) = caddn0(i) * x
      cminn(i) = cminn0(i) * x
   40 continue
!fong02     
!
!fong03
      do 50 i = 1, nsn
      dtau = tau(i+ns) - tau(i)
      expn(i) = exp ( - lamdan(i) * dtau )
      xn(i) = gamman(i) * expn(i)
      yn(i) = ( expn(i) - gamman(i) ) / ( xn(i) - 1.0 )
      zn(i) = ( expn(i) + gamman(i) ) / ( xn(i) + 1.0 )
   50 continue
!fong03
!
!fong04
!
!  etri : tridiagonal matrix for eq(39) in toon et al. (1989)
!  etri(i,k,m) : i=1,ns (horizontal grids), k=1,2*lev (vertical levels)
!                mm=1,4 (1:a, 2:b, 3:d, 4:e, a,b,d,e refer to eq(39)) 
!
      do 110 i = 1, ns        
      etri(i,1,1) = 0.0
      etri(i,1,2) = xn(i) + 1.0
      etri(i,1,3) = xn(i) - 1.0
      etri(i,1,4) = - cminn0(i)
  110 continue
!
      do 150 j = 1, lev-1   
      jk = j*ns
      jj = jk-ns
      j2n = j+j
      j2n1= j2n+1
      do 150 i = 1, ns
      etri(i,j2n,1) = 1.0+xn(i+jj)-yn(i+jk)*(gamman(i+jj)+expn(i+jj))
      etri(i,j2n,2) = 1.0-xn(i+jj)-yn(i+jk)*(gamman(i+jj)-expn(i+jj))
      etri(i,j2n,3) = yn(i+jk)*(1.0+xn(i+jk))-expn(i+jk)-gamman(i+jk)
      etri(i,j2n,4) = caddn0(i+jk)-caddn(i+jj)-yn(i+jk) *  &
                ( cminn0(i+jk)-cminn(i+jj) )
      etri(i,j2n1,1) = gamman(i+jj)-expn(i+jj)-zn(i+jj)*(1.0-xn(i+jj))
      etri(i,j2n1,2) = -1.0-xn(i+jk)+zn(i+jj)*(expn(i+jk)+gamman(i+jk))
      etri(i,j2n1,3) = zn(i+jj)*(expn(i+jk)-gamman(i+jk))-xn(i+jk)+1.0
      etri(i,j2n1,4) = cminn0(i+jk)-cminn(i+jj)-zn(i+jj) *  &
                 ( caddn0(i+jk)-caddn(i+jj) )
  150 continue
!fong04
!
!fong05
      jke = lev*ns
      jje = jke-ns
      do 170 i = 1, ns
      ssfc = pi * u0(i) * exp(-tau(i+jke)/u0(i)) * as(i) * f0
      wm1 = 1.0 - as(i) * gamman(i+jje)
      wm2 = xn(i+jje) - as(i) * expn(i+jje)
      etri(i,lev2,1) = wm1 + wm2
      etri(i,lev2,2) = wm1 - wm2
      etri(i,lev2,3) = 0.0
      etri(i,lev2,4) = as(i) * cminn(i+jje) - caddn(i+jje) + ssfc
  170 continue
!fong05
!
!fong06
!  solve a tridiagonal matrix using gaussian elimination
!
      call trigau ( etri, nx, ns, lev2, lev2, 0 ) 
!fong06
!
!fong07
      do 180 j = 1, lev 
      jk = j*ns
      jj = jk-ns
      j2m1 = j+j-1
      j2m  = j+j
      do 180 i = 1, ns
      aa(i+jj) = etri(i,j2m1,4) + etri(i,j2m,4)
      bb(i+jj) = etri(i,j2m1,4) - etri(i,j2m,4)
  180 continue
!fong07
!
!fong08
      do 200 i = 1, ns
      xx = aa(i) * expn(i)
      ffu(i) = xx + gamman(i) * bb(i) + caddn0(i)
      ffd(i) = gamman(i) * xx + bb(i) + cminn0(i) + sol(i)
  200 continue
!
      do 210 i = 1, nsn
      xx = bb(i) * expn(i)
      ffu(i+ns) = aa(i) + gamman(i) * xx + caddn(i)
      ffd(i+ns) = gamman(i) * aa(i) + xx + cminn(i) + sol(i+ns)
  210 continue
!fong08
      return
      end
