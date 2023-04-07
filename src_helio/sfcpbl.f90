      subroutine sfcpbl ( mn,kk,dt,g,r,cp,hltm,tg,z0,ocean,ps,ts,hgt &
                        , u,v,t,q,ut,vt,tt,qt,sfcw,ustar,tstar,qstar &
                        , e,eps,hflux,qflux,qsfc,gwet,zl,itypbl      &
                        , t2,u10,v10 )
!
!##################################################################
!                     subroutine description
!
! 1. function description
!
!       this subroutine calculated surface fluxes of momentum,
!       heat and moisture.
!       the effects of surface fluxes are added to the bottom
!       levels.
!
! 2. parameter specification
!
!     mn : dimension of horizontal-direction
!     kk : dimension of z-direction
!     dt : time step for ut, vt, tt, qt updates               (s)
!      g : gravity                                            (m/s2)
!      r : air gas constant                                   (j/kg)
!     cp : air specific heat constant                         (j/kg)
!   hltm : water vapor latent heat constant                   (j/kg)
!     tg : ground potential temp. (sea or land)    (mn)      (k)
!     z0 : surface roughness length                (mn)      (m)
!  ocean : logical array for ocean  (mn)      (true = open water)
!     ps : surface pressure     (mn)     at current time lvl (mb)
!     ts : surface air temp
!    hgt : height above ground  (mn,kk)  at current time lvl (m)
!      u : x-direction velocity (mn,kk) at current time lvl (m/s)
!      v : y-direction velocity (mn,kk) at current time lvl (m/s)
!      t : potent temp on odd lvl(mn,kk) at current time lvl (k)
!      q : moisture    on odd lvl(mn,kk) at current time lvl (kg/kg)
!     ut : x-direction velocity (mn,kk) at past/future  lvl (m/s)
!     vt : y-direction velocity (mn,kk) at past/future  lvl (m/s)
!     tt : potent temp on odd lvl(mn,kk) at past/future  lvl (k)
!     qt : moisture    on odd lvl(mn,kk) at past/future  lvl (kg/kg)
!   sfcw : surface wind speed                       (mn)     (m/s)
!   ustar: surface friction velocity                (mn)     (m/s)
!   tstar: surface potential temp scale             (mn)     (k)
!   qstar: surface moisture scale                   (mn)     (kg/kg)
!     e  : turbulence kinetic energy                (mn,kk)  (m2/s2)
!    eps : tke dissipation                          (mn)     (m2/s2)
!   hflux: upward surface heat flux, =-ro*cp*ustar*tstar   (mn)  (w/m2)
!   qflux: upward surface moist flux,=-ro*hltm*ustar*qstar (mn)  (w/m2)
!    qsfc: surface satuated specific humidity       (mn)     (kg/kg)
!    gwet: ground wetness                           (mn)
!     zl : height(kk) over monin-obukhov scale height.
!          a stability indicator.                   (mn)
! itypbl : index for surface layer update, =0 stress b,c., =1 direct
!     t2 : surface air temperature at 2m             (mn)   (k)
!    u10 : surface wind in x-direction at 10m        (mn)   (m/s)
!    v10 : surface wind in y-direction at 10m        (mn)   (m/s)
!
! 4. local work arrays
!
!    phim: non-dimensional vertically integrated wind shear function
!          (businger , 1973)                        (mn)
!    phih: non-dimensional vertically integrated temperature/moisture
!          gradient function (businger , 1973)      (mn)
!    sfcw: wind speed on mass point                 (mn)     (m/s)
!    zsz0: log(zs/z0)                               (mn)     (m/m)
!
! 5. calling modules
!
!    pbltke
!
! 6. usage
!
!     call sfcpbl ( mn,kk,dt,g,r,cp,hltm,tg,z0,ocean,ps,ts,hgt,u,v,t
!    1            , q,ut,vt,tt,qt,sfcw,ustar,tstar,qstar,e,eps,hflux
!    2            , qflux,qsfc,gwet,zl,itypbl )
!
! 7. modules called
!
!
! 8. limitation
!
!   ** 1.5e-5 <= z0 (sea) <= 0.01
!   **   -2.0 <=    zl    <= 2.0
!   **            ustar   >= 1.e-7,  <= 1.e+1
!   **            tatrt              <= 5.e+1
!   **            qatrt              <= 5.e+0
!   **              zs    >= 10.0
!   **             sfcw   >= 0.01
!
! 9. date
!
!     created            10/1/1989
!     modified           7/23/1990
!     rewritten          nov. 1991
!
! 10. author
!      c. liou / f. wang / s. chang
!      modify to f90 by C-H Lee and sort by River Chen in 2015
!
! 11. reference
!       businger j.a.,1973: turbulent transfer in the atmospheric
!            surface layer.  workshop in meteorology, d. a. haugen,
!            ed., american meteorology socity, 67-98.
!#####################################################################
      use paramt

      implicit  none
!
! input & output variable
!
      integer   mn,kk,itypbl

      real      tg(mn),z0(mn),ps(mn),hgt(mn,kk)                    &
              , u(mn,kk),v(mn,kk),t(mn,kk),q(mn,kk),ut(mn,kk)      &
              , vt(mn,kk),tt(mn,kk),qt(mn,kk),sfcw(mn),ustar(mn)   &
              , tstar(mn),qstar(mn),e(mn,kk),eps(mn,kk),hflux(mn)  &
              , qflux(mn),qsfc(mn),gwet(mn),zl(mn),ts(mn)

      logical ocean(mn)
!
!  local work arrays
!
      real      phim(im),phih(im),zsz0(im),rosfc(im)               &
              , wu(im),wv(im),wt(im),wq(im)                        &
!cc for output of t2, u10 and v10 ( 2003/01/20 )
              , u10(mn),v10(mn),t2(mn)
!
!
!          cccccccccccccccccccccccccccccccccccccccccccccc
!          c                                            c
!          c         similarity surface layer           c
!          c         paulson 1970, businger 1973        c
!          c                                            c
!          cccccccccccccccccccccccccccccccccccccccccccccc
!
      integer i,it1
      real    xkapa,hltm,cp,r,g,dt,g04,pihalf,xx,yy
      real    btest,bl,al,zl2,zl10,xx10,phim10,phih2,p2
      real    hpbl,ghpbl,wstar,zlx

      xkapa = r/cp
      do 10 i = 1, mn
      rosfc(i) = 100.0*ps(i)/( r*ts(i))
  10  continue
!
!           iterate 3 times to match u*, t*, q*, and z/l
!
      do 900 it1 = 1, 3
!
!           charnock's relation for surface roughness over oceans
!
      do 100 i = 1, mn
      if ( ocean(i) )  then
         z0(i) = 0.032 *ustar(i)*ustar(i)/g
         z0(i)= max(1.5e-5,min(z0(i),0.01))
      endif
! 100 continue
!     do 140 i=1,mn
!
!           compute surface wind
!
      sfcw(i)= max( 1.0e-2,sqrt( ut(i,kk)*ut(i,kk)+vt(i,kk)*vt(i,kk)))
! 140 continue
!
!           compute z/l
!
!           limit minimum depth of surface layer
!           (saved in zsz0) for surface flux calculation
!
!
      g04 = g*0.4
!     do 200 i=1,mn
      zsz0(i) = max(10.0,hgt(i,kk))
      zl(i) = g04*zsz0(i)*tstar(i)*(1.+0.61*qstar(i)) /     &
              (tt(i,kk)*ustar(i)*ustar(i))
      zl(i)= max(-2.0,min(zl(i),2.0))
! 200 continue
!
!           compute phi functions  (see businger 1973, or paulson 1970)
!
      pihalf = 3.141592654 * 0.5
!     do 240 i = 1, mn
      if(zl(i) .lt. 0.0) then
         xx = exp(0.25*log(1.0 - 15.0*zl(i)))
         phim(i) = 2.0*log(0.5*(1.0+xx)) + log(0.5*(1.0+xx*xx))  &
                          - 2.0*atan(xx) + pihalf
         phih(i) = 2.0*log(0.5*(1.0 +sqrt(1.0 -9.0*zl(i))))
      else
         phim(i) = -4.7  * zl(i)
         phih(i) = -6.35 * zl(i)
      endif
! 240 continue
!
!     do 300 i = 1, mn
      zsz0(i) = log(zsz0(i)/z0(i))
      phim(i) = zsz0(i) - phim(i)
      phih(i) = zsz0(i) - phih(i)
      phim(i)= max(phim(i),1.0e-3)
      phih(i)= max(phih(i),1.0e-3)
! 300 continue
!
!     compute new ut, vt, tt, and qt for ustar, tstar, qstar updates
!     using crank-nicolson (bl=0.5) time scheme for weak forcing or
!        or backward (bl=1.) time scheme for strong forcing
!
!     do 340 i = 1, mn
      yy    = dt*ustar(i) * ustar(i) / (2.0*hgt(i,kk)*sfcw(i))
      btest = 0.5 + (yy-1.0)*0.5
      bl    = max (0.5, min (1.0, btest) )
      al    = 1.0 - bl
      wu(i) = ut(i,kk) * (1.0-al*yy) / (1.0+bl*yy)
      wv(i) = vt(i,kk) * (1.0-al*yy) / (1.0+bl*yy)
      xx    = 0.4*dt*ustar(i) / (0.74*2.0*hgt(i,kk)*phih(i))
      wt(i) = ( tt(i,kk)*(1.0-al*xx) + xx*tg(i) )   / (1.0+bl*xx)
      wq(i) = ( qt(i,kk)*(1.0-al*xx) + xx*qsfc(i) ) / (1.0+bl*xx)
! 340 continue
!
!     update ustar, tstar and qstar
!
!     do 360 i = 1, mn
! phon-ju
! add blending for stability
!
      wu(i)=0.5*(ut(i,kk)+wu(i))
      wv(i)=0.5*(vt(i,kk)+wv(i))
      wt(i)=0.5*(tt(i,kk)+wt(i))
      wq(i)=0.5*(qt(i,kk)+wq(i))
! add over
      sfcw(i) = max (1.0e-2, sqrt( wu(i)*wu(i) + wv(i)*wv(i) ))
      ustar(i) = 0.4*sfcw(i) / phim(i)
      tstar(i) = 0.4*        (wt(i)-tg(i))  /(0.74*phih(i))
      qstar(i) = 0.4*        (wq(i)-qsfc(i))/(0.74*phih(i))
      ustar(i)= max(ustar(i),1.0e-7)
      ustar(i)= min(ustar(i),10.000)
      tstar(i)= min(tstar(i),50.00)
      qstar(i)= min(qstar(i),5.000)
! 360 continue
!
!
!     do 400 i = 1, mn
      hflux(i) = -rosfc(i)*cp*ustar(i)*tstar(i)
      qflux(i) = -rosfc(i)*hltm*ustar(i)*qstar(i)
! 400 continue
  100 continue
!
  900 continue
!
!ccc for output of t2, u10 and v10
!
      do i=1,mn
      zl2=zl(i)*2./hgt(i,kk)
      zl10=zl(i)*10./hgt(i,kk)
      if(zl(i) .lt. 0.0) then
         xx10 = exp(0.25*log(1.0 - 15.0*zl10))
         phim10 = 2.0*log(0.5*(1.0+xx10)) + log(0.5*(1.0+xx10*xx10))  &
                          - 2.0*atan(xx10) + pihalf
         phih2 = 2.0*log(0.5*(1.0 +sqrt(1.0 -9.0*zl2)))
      else
         phim10 = -4.7  * zl10
         phih2 = -6.35 * zl2
      endif
      phim10 =log( 10./z0(i)) - phim10
      phih2 =log( 2./z0(i)) - phih2
      phim10= max(phim10,1.0e-3)
      phih2= max(phih2,1.0e-3)
!
      v10(i)=wv(i)*phim10/phim(i)
      u10(i)=wu(i)*phim10/phim(i)
!
      t2(i)=tg(i)*(1.-phih2/phih(i)) + wt(i)*phih2/phih(i)
!
! convert t2 from potential temp. to temp
      p2=(ps(i)*100.-rosfc(i)*g*2.)/100.
      t2(i)=t2(i)*(p2/1000.)**xkapa
      enddo
!cc
!
!     update ut, vt, tt, and qt, if itypbl=1
!
      if ( itypbl .eq. 1 )  then
         do 420 i = 1, mn
         ut(i,kk) = wu(i)
         vt(i,kk) = wv(i)
         tt(i,kk) = wt(i)
         qt(i,kk) = wq(i)
  420    continue
      endif
!
!           compute surface layer e and eps
!
      hpbl  = 450.0
      ghpbl = g*hpbl
      do 500 i = 1, mn
      wstar = ghpbl*(-ustar(i)*tstar(i)) / tt(i,kk)
      wstar = max(wstar,1.0e-20)
      wstar = exp(0.3333*log(wstar))
      e(i,kk)= 3.75*ustar(i)*ustar(i) + 0.2*wstar*wstar
  500 continue
      do 520 i = 1, mn
      zlx      = min(zl(i),-1.0e-20)
      e(i,kk)  = e(i,kk) + ustar(i)*ustar(i)*exp(0.6666*log(-zlx))
      eps(i,kk)= ustar(i)**3 / (0.4*hgt(i,kk))
  520 continue
!
      do 540 i = 1, mn
      e(i,kk)  = min (250.0, max (1.0e-4, e(i,kk)))
      eps(i,kk)= min (1.0,   max (1.0e-7, eps(i,kk)))
  540 continue
!
      return
      end
