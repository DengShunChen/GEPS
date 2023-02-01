      subroutine rcmpsl2c ( ns,nx,lev,u0,as,ww1,ww2,ww,ttau,cld,cldb    &
                           ,ffu,ffd )
!
! **********************************************************************
! (1)
! The generalized two stream approximation for nonhomgeneous atmospheres
! in  the  solar  wavelengths.
!
! (2)
! The treatment of cloud overlapping for partial cloudinesses is included
! by Fong(1996), which is based on the method proposed by Geleyn and
! Hollingsworth (1979, contrib. atmos. phys. )
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
!    c. t. fong and fu (1996, summer)
!    modify to f90 by C-H Lee in 2015
!
! **********************************************************************
!
      use fj_pad

      implicit  none
      integer   ns,nx,lev

      real      u0(nx),as(nx),ww1(nx*lev,2),ww2(nx*lev,2),ww(nx*lev,2), &
                ttau(nx*(lev+1),2),cld(nx*lev),cldb(nx*lev,4)
      real      ffu(nx*(lev+1)),ffd(nx*(lev+1))
!
!   local working arrays
!
!fong include 'param.h'
!
!     include '../include/fj_pad.h'
      real      w1(nx*lev,2),w(nx*lev,2),tau(nx*(lev+1),2)
      real      lamdan
      real      expdt(nx*lev,2),sol(nx*(lev+1),2)
      real      rl(nx*lev,2),tl(nx*lev,2),rb(nx*lev,2),tb(nx*lev,2)
!fj      dimension aa(nx,4,4),amtx(nx,lev*4,12),u0t(nx*lev)
      real      amtx(nx+npad,lev*4,12),u0t(nx*lev)

      logical edding,quadra
!ibm
!fj      dimension tmp1(ns*lev),tmp2(ns*lev)
      real      tmp1(ns*lev)
!
      integer   nsn,nsn1,lev4,i,j,jk,jj,m,kini
      real      sqrt3,rsqrt3,pi,f0,fxx,fw,tmp2,x,y,z
      real      gamma1,gamma2,gamma3,gamma4,ugts1,gamman
      real      u0_mod,vh,vk,enn,axx,bxx,cxx,solm_f,solm_c,ssfc_f,ssfc_c

      edding = .true.
      quadra = .false.
!
      nsn  = ns*lev
      nsn1 = ns*(lev+1)
      lev4 = lev*4
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
      do 7 i = 1, ns
      sol(i,1) = u0(i)
      sol(i,2) = 0.
    7 continue
!
      do j = 1, lev
      jk = j*ns
      jj = jk-ns
      do i = 1, ns
       u0t(i+jj) = u0(i)
      enddo
      enddo
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

      do 40 m = 1, 2
!ibm
      do i = 1, ns
       tau(i,m) = 0.
      enddo
!
      do i = 1, nsn
       fxx = ww2(i,m) / 5.0
       fw = 1.0 - fxx * ww(i,m)
       w1(i,m) = ( ww1(i,m) - 3.0 * fxx ) / ( 1.0 - fxx )
       w(i,m)  = ( 1.0 - fxx ) * ww(i,m) / fw
       w(i,m) = min ( w(i,m), 0.999999 )
       tmp1(i) = (ttau(i+ns,m) - ttau(i,m)) * fw
      enddo

      do i = 1, nsn
       tau(i+ns,m) = tau(i,m) + tmp1(i)
      enddo
!ibm
!fj      call vrec(tmp2(1),u0t(1),nsn)
!fj      call vexp(expdt(1,m),tmp2(1),nsn)
!ocl norecurrence,simd
      do 40 i = 1, nsn
         tmp2 = 1.d0/u0t(i)
         tmp2 = -tmp1(i)*tmp2
         expdt(i,m) = exp(tmp2)
!ibm
!fj      do 40 i = 1, nsn
!
!fong02_1
!cc   if ( edding ) then
	   x = 0.25 * w1(i,m)     ! w1 = 3g
           y = w(i,m) * x
           gamma1 = 1.75 - w(i,m) - y
           gamma2 = - 0.25 + w(i,m) - y
           gamma3 = 0.5 - x * u0t(i)
           gamma4 = 1.0 - gamma3
           ugts1 = 0.5
!
!cc   elseif ( quadra ) then
!          x = 0.5*sqrt3 * w(i,m)
!          y = 0.5*rsqrt3 * w1(i,m)  ! rsqrt3 = 1/sqrt3
!          z = y * w(i,m)
!          gamma1 = sqrt3 - x - z
!          gamma2 = x - z
!          gamma3 = 0.5 - y * u0t(i)
!          gamma4 = 1.0 - gamma3
!          ugts1 = rsqrt3
!cc   endif
!fong02_1
!ibm  lamdan = exp( 0.5*log((gamma1 + gamma2)*(gamma1 - gamma2)) )
      lamdan = sqrt((gamma1 + gamma2)*(gamma1 - gamma2))
      gamman = gamma2 / ( gamma1 + lamdan )
      z = lamdan * lamdan * u0t(i) * u0t(i) - 1.0
      if(abs(z) .lt. 0.000001)then
       u0_mod=u0t(i)+0.01*u0t(i)
       z = lamdan * lamdan * u0_mod * u0_mod - 1.0
      end if
      vh = ( -gamma3 + u0t(i)*(gamma1*gamma3+gamma2*gamma4) ) / z
      vk = (  gamma4 + u0t(i)*(gamma1*gamma4+gamma2*gamma3) ) / z
      enn = exp ( - lamdan * ( tau(i+ns,m)-tau(i,m) ) )
      axx = gamman * gamman
      bxx = enn * enn
      cxx = axx*bxx - 1.0
      rl(i,m) = gamman * ( bxx - 1.0 ) / cxx
      tl(i,m) = enn * ( axx - 1.0 ) / cxx
      rb(i,m) = w(i,m) * ( -vk*rl(i,m) - vh*expdt(i,m)* tl(i,m) + vh )
      tb(i,m) = w(i,m) * ( -vk*tl(i,m) - vh*expdt(i,m)* rl(i,m) + expdt(i,m)*vk )
      rl(i,m) = max ( 0., rl(i,m) )
   40 continue
!fong02
!
!ibm
      do 26 j = 1, lev
      jk = j*ns
      jj = jk-ns
      do 26 i = 1, ns
      sol(i+jk,1) = expdt(i+jj,1) * ( cldb(i+jj,1)*sol(i+jj,1) +        &
                    (1.0-cldb(i+jj,3))*sol(i+jj,2) )
      sol(i+jk,2) = expdt(i+jj,2) * ( (1.0-cldb(i+jj,1))*sol(i+jj,1) +  &
                    cldb(i+jj,3)*sol(i+jj,2) )
   26 continue
!ibm
!     do 45 i = 1, nx*lev4*12
!     amtx(i,1,1) = 0.
!  45 continue
!
!fong03

!fj
!fj j = 1
!fj
      j = 1
      kini = (j-1) * 4
      jk = j*ns
      jj = jk-ns

!
      do i = 1, ns
         amtx(i,kini+1,3)=0.
         amtx(i,kini+1,4)=0.
         amtx(i,kini+1,5)=0.
         amtx(i,kini+1,7)=0.
         amtx(i,kini+1,8)=0.
         amtx(i,kini+1,9)=0.
         amtx(i,kini+2,3)=0.
         amtx(i,kini+2,4)=0.
         amtx(i,kini+2,5)=0.
         amtx(i,kini+2,7)=0.
         amtx(i,kini+2,8)=0.
         amtx(i,kini+3,3)=0.
         amtx(i,kini+3,4)=0.
         amtx(i,kini+3,5)=0.
         amtx(i,kini+3,7)=0.
         amtx(i,kini+4,3)=0.
         amtx(i,kini+4,4)=0.
         amtx(i,kini+4,5)=0.
         amtx(i,kini+4,9)=0.
         amtx(i,kini+1,6) = 1.0
         amtx(i,kini+1,10) =  - tl(i+jj,1) * cldb(i+jj,2)
         amtx(i,kini+1,11) =  - tl(i+jj,1) * ( 1.0 - cldb(i+jj,4) )
         amtx(i,kini+2,6) = 1.0
         amtx(i,kini+2,9) = - tl(i+jj,2) * ( 1.0 - cldb(i+jj,2) )
         amtx(i,kini+2,10) = - tl(i+jj,2) * cldb(i+jj,4)
         amtx(i,kini+3,6) = 1.0
         amtx(i,kini+3,8) = - rl(i+jj,1) * cldb(i+jj,2)
         amtx(i,kini+3,9) = - rl(i+jj,1) * ( 1.0 - cldb(i+jj,4) )
         amtx(i,kini+4,6) = 1.0
         amtx(i,kini+4,7) = - rl(i+jj,2) * ( 1.0 - cldb(i+jj,2) )
         amtx(i,kini+4,8) = - rl(i+jj,2) * cldb(i+jj,4)
      enddo

      do i = 1, ns
         solm_f = cldb(i+jj,1)*sol(i+jj,1) + ( 1.0-cldb(i+jj,3) ) * sol(i+jj,2)
         solm_c = ( 1.0-cldb(i+jj,1) )*sol(i+jj,1) + cldb(i+jj,3) * sol(i+jj,2)
         amtx(i,kini+1,12) = rb(i+jj,1) * solm_f
         amtx(i,kini+2,12) = rb(i+jj,2) * solm_c
         amtx(i,kini+3,12) = tb(i+jj,1) * solm_f
         amtx(i,kini+4,12) = tb(i+jj,2) * solm_c
      enddo

!fj
!fj j = 2,lev-1
!fj

      do j = 2,lev-1
         kini = (j-1) * 4
         jk = j*ns
         jj = jk-ns
!
!
         do i = 1, ns
            amtx(i,kini+1,3)=0.
            amtx(i,kini+1,7)=0.
            amtx(i,kini+1,8)=0.
            amtx(i,kini+1,9)=0.
            amtx(i,kini+2,5)=0.
            amtx(i,kini+2,7)=0.
            amtx(i,kini+2,8)=0.
            amtx(i,kini+3,4)=0.
            amtx(i,kini+3,5)=0.
            amtx(i,kini+3,7)=0.
            amtx(i,kini+4,3)=0.
            amtx(i,kini+4,4)=0.
            amtx(i,kini+4,5)=0.
            amtx(i,kini+4,9)=0.
         enddo
         do i = 1, ns
            amtx(i,kini+1,4)  = - rl(i+jj,1) * cldb(i+jj,1)
            amtx(i,kini+1,5)  = - rl(i+jj,1) * ( 1.0 - cldb(i+jj,3) )
            amtx(i,kini+1,10) = - tl(i+jj,1) * cldb(i+jj,2)
            amtx(i,kini+1,11) = - tl(i+jj,1) * ( 1.0 - cldb(i+jj,4) )
            amtx(i,kini+2,3)  = - rl(i+jj,2) * ( 1.0 - cldb(i+jj,1) )
            amtx(i,kini+2,4)  = - rl(i+jj,2) * cldb(i+jj,3)
            amtx(i,kini+2,9)  = - tl(i+jj,2) * ( 1.0 - cldb(i+jj,2) )
            amtx(i,kini+2,10) = - tl(i+jj,2) * cldb(i+jj,4)
            amtx(i,kini+3,2)  = - tl(i+jj,1) * cldb(i+jj,1)
            amtx(i,kini+3,3)  = - tl(i+jj,1) * ( 1.0 - cldb(i+jj,3) )
            amtx(i,kini+3,8)  = - rl(i+jj,1) * cldb(i+jj,2)
            amtx(i,kini+3,9)  = - rl(i+jj,1) * ( 1.0 - cldb(i+jj,4) )
            amtx(i,kini+4,1)  = - tl(i+jj,2) * ( 1.0 - cldb(i+jj,1) )
            amtx(i,kini+4,2)  = - tl(i+jj,2) * cldb(i+jj,3)
            amtx(i,kini+4,7)  = - rl(i+jj,2) * ( 1.0 - cldb(i+jj,2) )
            amtx(i,kini+4,8)  = - rl(i+jj,2) * cldb(i+jj,4)
         enddo

         do i = 1, ns
            amtx(i,kini+1,6) = 1.0
            amtx(i,kini+2,6) = 1.0
            amtx(i,kini+3,6) = 1.0
            amtx(i,kini+4,6) = 1.0
         enddo

         do i = 1,ns
            solm_f = cldb(i+jj,1)*sol(i+jj,1) + ( 1.0-cldb(i+jj,3) ) * sol(i+jj,2)
            solm_c = ( 1.0-cldb(i+jj,1) )*sol(i+jj,1) + cldb(i+jj,3) * sol(i+jj,2)
            amtx(i,kini+1,12) = rb(i+jj,1) * solm_f
            amtx(i,kini+2,12) = rb(i+jj,2) * solm_c
            amtx(i,kini+3,12) = tb(i+jj,1) * solm_f
            amtx(i,kini+4,12) = tb(i+jj,2) * solm_c
         enddo
      enddo

!fj
!fj j = lev
!fj
      j = lev
      kini = (j-1) * 4
      jk = j*ns
      jj = jk-ns
!
!
      do i = 1, ns

         amtx(i,kini+1,3)=0.
         amtx(i,kini+1,4)=0.
         amtx(i,kini+1,5)=0.
         amtx(i,kini+1,7)=0.
         amtx(i,kini+1,8)=0.
         amtx(i,kini+1,9)=0.

         amtx(i,kini+2,3)=0.
         amtx(i,kini+2,4)=0.
         amtx(i,kini+2,5)=0.
         amtx(i,kini+2,7)=0.
         amtx(i,kini+2,8)=0.
         amtx(i,kini+2,9)=0.

         amtx(i,kini+3,3)=0.
         amtx(i,kini+3,4)=0.
         amtx(i,kini+3,5)=0.
         amtx(i,kini+3,7)=0.
         amtx(i,kini+3,8)=0.
         amtx(i,kini+3,9)=0.

         amtx(i,kini+4,3)=0.
         amtx(i,kini+4,4)=0.
         amtx(i,kini+4,5)=0.
         amtx(i,kini+4,7)=0.
         amtx(i,kini+4,8)=0.
         amtx(i,kini+4,9)=0.
         amtx(i,kini+1,6) = 1.0
         amtx(i,kini+1,4) =  - rl(i+jj,1) * cldb(i+jj,1)
         amtx(i,kini+1,5) = - rl(i+jj,1) * ( 1.0 - cldb(i+jj,3) )
         amtx(i,kini+1,8) = amtx(i,kini+1,8) + as(i)*(- tl(i+jj,1) * cldb(i+jj,2))
         amtx(i,kini+1,9) = amtx(i,kini+1,9) + as(i)*(- tl(i+jj,1) * ( 1.0 - cldb(i+jj,4) ))
         amtx(i,kini+2,6) = 1.0

         amtx(i,kini+2,3) = - rl(i+jj,2) * ( 1.0 - cldb(i+jj,1) )
         amtx(i,kini+2,4) = - rl(i+jj,2) * cldb(i+jj,3)
         amtx(i,kini+2,7) = amtx(i,kini+2,7) + as(i)*(- tl(i+jj,2) * ( 1.0 - cldb(i+jj,2) ))
         amtx(i,kini+2,8) = amtx(i,kini+2,8) + as(i)*(- tl(i+jj,2) * cldb(i+jj,4))
         amtx(i,kini+3,6) = 1.0

         amtx(i,kini+3,2) = - tl(i+jj,1) * cldb(i+jj,1)
         amtx(i,kini+3,3) = - tl(i+jj,1) * ( 1.0 - cldb(i+jj,3) )

         amtx(i,kini+3,6) = amtx(i,kini+3,6) + as(i)*(- rl(i+jj,1) * cldb(i+jj,2))
         amtx(i,kini+3,7) = amtx(i,kini+3,7) + as(i)*(- rl(i+jj,1) * ( 1.0 - cldb(i+jj,4) ))
         amtx(i,kini+4,6) = 1.0

         amtx(i,kini+4,1) = - tl(i+jj,2) * ( 1.0 - cldb(i+jj,1) )
         amtx(i,kini+4,2) = - tl(i+jj,2) * cldb(i+jj,3)

         amtx(i,kini+4,5) = amtx(i,kini+4,5) + as(i)*(- rl(i+jj,2) * ( 1.0 - cldb(i+jj,2) ))
         amtx(i,kini+4,6) = amtx(i,kini+4,6) + as(i)*(- rl(i+jj,2) * cldb(i+jj,4))
      enddo

        do i = 1, ns
           solm_f = cldb(i+jj,1)*sol(i+jj,1) + ( 1.0-cldb(i+jj,3) ) * sol(i+jj,2)
           solm_c = ( 1.0-cldb(i+jj,1) )*sol(i+jj,1) + cldb(i+jj,3) * sol(i+jj,2)
           amtx(i,kini+1,12) = rb(i+jj,1) * solm_f
           amtx(i,kini+2,12) = rb(i+jj,2) * solm_c
           amtx(i,kini+3,12) = tb(i+jj,1) * solm_f
           amtx(i,kini+4,12) = tb(i+jj,2) * solm_c
        enddo

        do i = 1, ns
           ssfc_f = as(i) * sol(i+nsn,1) ! assume direct albedo = as
           ssfc_c = as(i) * sol(i+nsn,2) 	
           amtx(i,kini+1,12) = amtx(i,kini+1,12) - ( (- tl(i+jj,1) * cldb(i+jj,2))*ssfc_f +  &
                (- tl(i+jj,1) * ( 1.0 - cldb(i+jj,4) ))*ssfc_c )
           amtx(i,kini+2,12) = amtx(i,kini+2,12) - ( (- tl(i+jj,2) * ( 1.0 - cldb(i+jj,2) ))*ssfc_f +  &
                (- tl(i+jj,2) * cldb(i+jj,4))*ssfc_c )
           amtx(i,kini+3,12) = amtx(i,kini+3,12) - ( (- rl(i+jj,1) * cldb(i+jj,2))*ssfc_f +  &
                (- rl(i+jj,1) * ( 1.0 - cldb(i+jj,4) ))*ssfc_c )
           amtx(i,kini+4,12) = amtx(i,kini+4,12) - ( (- rl(i+jj,2) * ( 1.0 - cldb(i+jj,2) ))*ssfc_f +  &
                (- rl(i+jj,2) * cldb(i+jj,4))*ssfc_c )
        enddo
!fong03
!
!  solve matrix
!
      call mxslvc ( amtx, nx, ns, lev4 )
!
!fong04
!
!  compute net upward and downward fluxes at each level
!
!  level 1
!
      do 210 i = 1, ns
      ffu(i) = amtx(i,1,12) + amtx(i,2,12)
      ffd(i) = sol(i,1) + sol(i,2)
  210 continue
!
!  level 2 - lev
!
      do 220 j = 1, lev-1
      kini = 2 + (j-1)*4
      jk = j*ns
      jj = jk-ns
      do 220 i = 1, ns
      ffd(i+jk) = amtx(i,kini+1,12) + amtx(i,kini+2,12) + sol(i+jk,1) + sol(i+jk,2)
      ffu(i+jk) = amtx(i,kini+3,12) + amtx(i,kini+4,12)
  220 continue
!
!  level lev+1 ( ground )
!
      do 250 i = 1, ns
      ffd(i+nsn) = amtx(i,lev4-1,12) + amtx(i,lev4,12) + sol(i+nsn,1) + sol(i+nsn,2)
      ffu(i+nsn) = as(i) * ( amtx(i,lev4-1,12) + amtx(i,lev4,12) )      &
                   + as(i) * ( sol(i+nsn,1) + sol(i+nsn,2) )
  250 continue
!
      return
      end
