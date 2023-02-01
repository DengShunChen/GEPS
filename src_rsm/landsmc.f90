      subroutine landsmc(nxj,mn,km,rhscnpy,rhsmc,ai,bi,ci,smc,land      &
                        ,canopy,precip,runof,snowmt                     &
                        ,zsoil,istyp,sigmaf,dt,dth,tsat)
!
!##################################################################
!                     subroutine description
!
! 1. function description
!
!  nmcsfc_gwt is the second part of the soil model that is executed
!  after precipitation for the time step has been calculated.
!  In this subroutine canopy and soil moisture(smc) will be updated.
!
! 2. parameter specification
!
!     mn : dimension of horizontal-direction
!     km : dimension of soil layer
!     dt : time step for ut, vt, tt, qt updates                  (s)
!    dth : time step for snr updates                             (s)
!   smc : volumetric soil moisture content (mn,km)
! canopy : canopy moisture content (mn)                          (mm)
!rhscnpy : right hand side terms for  canopy prediction (mn)     (mm/s)
! rhsmc  : right hand side terms for soil moisture prediction (mn,km)(/s)
! ai     : first column of tri-diagnal matrix for solving smc (mn,km)
! bi     : second column of tri-diagonal matrix for solving smc (mn,km)
! ci     : third column of tri-diagonal matrix for solving smc (mn,km)
! land   : logical array for land (mn)                    (true = land  )
! drain  : drainage from bottom soil layer  (mn)                 (mm/s)
! snowmt : snow meltin rate (mn)                                 (mm/s)
! precip : total precipitation in each call(dth) (mn)            (mm)
! runof  : run off rate from soil layers (mn)                    (mm/s)
! canopy : canopy moisture content (mn)                          (mm)
! sigmaf : vegetation fraction (mn)
! istyp  : soil type(1-9) (mn)
! zsoil  : soil layer depth(mn,km)                                (m)
! tsat   : satuation point of 9 soil types (9)
!
! 3. local variables
!
!  prcp  : precipitation rate  (mn)                           (mm/s)
!  inf   : infiltration rate   (mn)                           (mm/s)
! infmax : maxmun infiltration rate (mn)                      (mm/s)

! 5. calling modules
!
!   pbltke
!
! 6. usage
!
!     call landsmc(mn,km,rhscnpy,rhsmc,ai,bi,ci,smc,land,
!    &       canopy,precip,runof,snowmt,
!    &       zsoil,istyp,sigmaf,dt,dth)
! 7. modules called
!
!
! 8. limitation
!
! 9. date
!
!     from ncep 2002
!     modified  2004
!
! 10. author
!
! 11. reference
!
!    Mahrt, L., and H. -L. Pan,1984: Atwo-layer model of soil hydrology.
!          Boundary Layer Meteorol., 29, 1-20.
!    Pan, H.-L., and L.Mahrt, 1987: Interaction between soil hydrology
!          and boundary-layer development. Boundary Layer Meteorol.,
!          38,185-202.
!#####################################################################
!
      implicit   none

      real, parameter :: scanop=.5,rhoh2o=1000.
      real, parameter :: ctfil1=.5,ctfil2=1.-ctfil1
      real, parameter :: rffact=1.

! input & output variable
      integer   nxj,mn,km
      real      dt,dth

      real      rhscnpy(mn),rhsmc(mn,km)
      real      ai(mn,km),bi(mn,km),ci(mn,km)
      real      smc(mn,km),canopy(mn),precip(mn)
      real      sigmaf(mn),runof(mn)
      real      zsoil(mn,km),snowmt(mn)
      integer   istyp(mn)
      logical   flag(mn),land(mn)
!
      real      prcp(mn),inf(mn),infmax(mn)
      real      sat(mn),dsat(mn),ksat(mn)
      real      smsoil(mn,km),cnpy(mn),dx(mn),conkdt(mn)
! landpack
      integer,  parameter :: ntype=9
      real      tsat(ntype)
!
      integer   i,k
      real      delt2,delt,drip,cc,tdif

      delt2 = dt
      delt  = dth
!
!  precipitation rate is needed in unit of kg m-2 s-1
!
      do i = 1, nxj
        prcp(i) = rhoh2o*( precip(i)/delt + snowmt(i) )/1000.
        runof(i) = 0.
        cnpy(i) = canopy(i)
      enddo
!
!  update canopy water content
!        ( updated by delt ,final blending is not necessary)
!

      do i = 1, nxj
        if(land(i)) then
          rhscnpy(i) = rhscnpy(i) + sigmaf(i) * prcp(i)
          canopy(i) = canopy(i) + delt * rhscnpy(i)
          canopy(i) = max(canopy(i),0.)
          prcp(i) = prcp(i) * (1. - sigmaf(i))
          if(canopy(i).gt.scanop) then
            drip = canopy(i) - scanop
            canopy(i) = scanop
            prcp(i) = prcp(i) + drip / delt
          endif
!
!  calculate infiltration rate
!
          inf(i) = prcp(i)
          sat(i) = tsat(istyp(i))
!         dsat(i) = funcdf(tsat(i),soiltyp(i))
!         ksat(i) = funckt(tsat(i),soiltyp(i))
!         infmax(i) = -dsat(i) * (tsat(i) - smc(i,1))
!    &                / (.5 * zsoil(i,1))
!    &                + ksat(i)
          infmax(i) = (-zsoil(i,1)) *                                 &
                      ((sat(i) - smc(i,1)) / delt - rhsmc(i,1))       &
                      * rhoh2o
          infmax(i) = max(rffact*infmax(i),0.)
!         if(smc(i,1).ge.tsat(i)) infmax(i) = ksat(i)
!         if(smc(i,1).ge.tsat(i)) infmax(i) = zsoil(i,1) * rhsmc(i,1)
!
! new(SWB model)-----------------------------------------------
!         dx(i)=0.
!         dx(i)=(-zsoil(i,1))*(sat(i)-smc(i,1))
!         do k=2,km
!         dx(i)=dx(i)+(zsoil(i,k-1)-zsoil(i,k))*(sat(i)-smc(i,k))
!         enddo
!
! conkdt=dt(in day)*Kdt_ref*Ks/K_ref
!c                  Kdt_ref=3.   K_ref=2.0e-6
!          conkdt(i)=delt/86400.*3.*ksat(i)/2.*1.0e+6
!          xx2=exp(-conkdt(i))
!          print*,'conkdt,xx=',conkdt(i),xx2
!          infmax(i)=inf(i)*dx(i)*( 1.-xx2 )/
!    +                     (inf(i)+dx(i)*(1.-xx2))
!---------------------------------------------------------------
          if(inf(i).gt.infmax(i)) then
            runof(i) = inf(i) - infmax(i)
            inf(i) = infmax(i)
          endif
          inf(i) = inf(i) / rhoh2o
          rhsmc(i,1) = rhsmc(i,1) - inf(i) / zsoil(i,1)
        endif
      enddo
!
!  we currently ignore the effect of rain on sea ice
!
      do i = 1, nxj
        flag(i) = land(i)
      enddo
!
!  solve the tri-diagonal matrix
!
      do k = 1, km
        do i = 1, nxj
          if(flag(i))  then
            rhsmc(i,k) = rhsmc(i,k) * delt2
            ai(i,k) = ai(i,k) * delt2
            bi(i,k) = 1. + bi(i,k) * delt2
            ci(i,k) = ci(i,k) * delt2
          endif
        enddo
      enddo
!  forward elimination
      do i = 1, nxj
        if(flag(i)) then
          ci(i,1) = -ci(i,1) / bi(i,1)
          rhsmc(i,1) = rhsmc(i,1) / bi(i,1)
        endif
      enddo
      do k = 2, km
        do i = 1, nxj
          if(flag(i)) then
            cc = 1. / (bi(i,k) + ai(i,k) * ci(i,k-1))
            ci(i,k) = -ci(i,k) * cc
            rhsmc(i,k) = (rhsmc(i,k) - ai(i,k) * rhsmc(i,k-1)) * cc
          endif
        enddo
      enddo
!  backward substituttion
      do i = 1, nxj
        if(flag(i)) then
          ci(i,km) = rhsmc(i,km)
        endif
      enddo
      do k = km-1, 1,-1
        do i = 1, nxj
          if(flag(i)) then
            ci(i,k) = ci(i,k) * ci(i,k+1) + rhsmc(i,k)
          endif
        enddo
      enddo
 100  continue
!
!  update soil moisture
!
      do k = 1, km
        do i = 1, nxj
          if(flag(i)) then
            smsoil(i,k) = smc(i,k) + ci(i,k)
            smsoil(i,k) = max(smsoil(i,k),0.)
            tdif = max(smsoil(i,k) - sat(i),0.)
            runof(i) = runof(i) -                           &
                      rhoh2o * tdif * zsoil(i,k) / delt
            smsoil(i,k) = smsoil(i,k) - tdif
          endif
        enddo
      enddo
!
! first update smc by delt2
! then blend smsoil and smc to get (t+delt) smc
!
      do k = 1, km
        do i = 1, nxj
          if(flag(i)) then
            smc(i,k) = ctfil1 * smsoil(i,k) + ctfil2 * smc(i,k)
          endif
        enddo
      enddo
!
!     do i = 1, nxj
!       if(flag(i)) then
!         canopy(i) = ctfil1 * canopy(i) + ctfil2 * cnpy(i)
!       endif
!     enddo
!     i = 1
!     print *, ' smc'
!     print 6000, smc(i,1), smc(i,2)
 6000 format(2(f8.5,','))
      return
      end
