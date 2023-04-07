      subroutine landpack(tsat,dfkt,xktk,dfk)
!
!*****************************************
! 1. purpose:
!
!  set up parameters table for 9 soil types
!                            ,depend on 22 levels of soil moisture
!     moisture diffusivity : dfk(22,9)
!     hydrolic conductivity : xktk(22,9)
!     thermal diffusivity for all soil types: dfkt(22,9)
! 2. parameter specification
!
! tsat   : satuation point  (9)
! dfkt   : soil thermal diffusivity     (22,9)    m2/s
! xktk   : soil hydraulic conductivity  (22,9)    m/s
!  dfk   : soil hydraulic diffusivity   (22,9)    m2/s
!
!***************************
!  the 9 soil types are:
!    1  ... loamy sand (coarse)
!    2  ... silty clay loam (medium)
!    3  ... light clay (fine)
!    4  ... sandy loam (coarse-medium)
!    5  ... sandy clay (coarse-fine)
!    6  ... clay loam  (medium-fine)
!    7  ... sandy clay loam (coarse-med-fine)
!    8  ... loam  (organic)
!    9  ... ice (use loamy sand property)
!****************************************
!
      implicit  none
      integer,  parameter :: ntype=9, ngrid=22 
      real      b(ntype),satpsi(ntype),satkt(ntype)
      real      tsat(ntype)
      real      dfkt(ngrid,ntype)
      real      dfk(ngrid,ntype)
      real      xktk(ngrid,ntype)

      integer   k,i
      real      dynw,f1,f2,theta,pf

!
      data b/4.26,8.72,11.55,4.74,10.73,8.17,6.77,5.25,4.26/
      data satpsi/.04,.62,.47,.14,.10,.26,.14,.36,.04/
      data satkt/1.41e-5,.20e-5,.10e-5,.52e-5,.72e-5,  &
                 .25e-5,.45e-5,.34e-5,1.41e-5/

! get dfkt
!
      do k = 1, ntype
        dynw = tsat(k) * .05
        f1 = log10(satpsi(k)) + b(k) * log10(tsat(k)) + 2.
        do i = 1, ngrid
          theta = float(i-1) * dynw
          theta = min(tsat(k),theta)
          if(theta.gt.0.) then
            pf = f1 - b(k) * log10(theta)
          else
            pf = 5.2
          endif
          if(pf.le.5.1) then
            dfkt(i,k) = exp(-(2.7+pf)) * 420.
          else
            dfkt(i,k) = .1744
          endif
        enddo
      enddo
!
! get dfk,ktk
!
      do k = 1, ntype
        dynw = tsat(k) * .05
        f1 = b(k) * satkt(k) * satpsi(k) / tsat(k) ** (b(k) + 3.)
        f2 = satkt(k) / tsat(k) ** (b(k) * 2. + 3.)
!
!  convert from m/s to kg m-2 s-1 unit
!
        f1 = f1 * 1000.
        f2 = f2 * 1000.
        do i = 1, ngrid
          theta = float(i-1) * dynw
          theta = min(tsat(k),theta)
          dfk(i,k) = f1 * theta ** (b(k) + 2.)
          xktk(i,k) = f2 * theta ** (b(k) * 2. + 3.)
        enddo
      enddo
      return
      end
