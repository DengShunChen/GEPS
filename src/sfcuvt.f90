      subroutine sfcuvt ( nxj,mn,dt,g,r,cp,hltm,tg,z0,ocean,ps,ts,hgt &
                        , ut,vt,tt,qt,ustar,tstar,qstar,hflux,qflux   &
                        , qsfc,t2,q2,rh2,rh10,u10,v10 )
!
!##################################################################
!                     subroutine description
!
! 1. function description
!
!       this subroutine calculated surface t2 u10 v10
!
! 2. parameter specification
!
!     mn : dimension of horizontal-direction
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
!     ut : x-direction velocity (mn,kk) at past/future  lvl (m/s)
!     vt : y-direction velocity (mn,kk) at past/future  lvl (m/s)
!     tt : potent temp on odd lvl(mn,kk) at past/future  lvl (k)
!     qt : moisture    on odd lvl(mn,kk) at past/future  lvl (kg/kg)
!   sfcw : surface wind speed                       (mn)     (m/s)
!   ustar: surface friction velocity                (mn)     (m/s)
!   tstar: surface potential temp scale             (mn)     (k)
!   qstar: surface moisture scale                   (mn)     (kg/kg)
!   hflux: upward surface heat flux, =-ro*cp*ustar*tstar   (mn)  (w/m2)
!   qflux: upward surface moist flux,=-ro*hltm*ustar*qstar (mn)  (w/m2)
!    qsfc: surface satuated specific humidity       (mn)     (kg/kg)
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
!     zl : height(kk) over monin-obukhov scale height.
!          a stability indicator.                   (mn)
!
! 5. calling modules
!
!    getrdy.f
!
! 6. usage
!
!     call sfcuvt ( mn,kk,dt,g,r,cp,hltm,tg,z0,ocean,ps,ts,hgt
!    1            , ut,vt,tt,qt,sfcw,ustar,tstar,qstar,hflux,qflux
!    2            , qsfc,t2,u10,v10 )
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
!
! 11. reference
!       businger j.a.,1973: turbulent transfer in the atmospheric
!            surface layer.  workshop in meteorology, d. a. haugen,
!            ed., american meteorology socity, 67-98.
!#####################################################################

      implicit   none
!
! input & output variable
!
      integer nxj,mn

      real    dt,g,r,cp,hltm

      real    tg(mn),z0(mn),ps(mn),hgt(mn)                           &
              , ut(mn),vt(mn),tt(mn),qt(mn),sfcw(mn)                 &
              , ustar(mn),tstar(mn),qstar(mn)                        &
              , qsfc(mn),ts(mn),hflux(mn),qflux(mn)                  &
              , u10(mn),v10(mn),t2(mn),t10(mn),rh2(mn),rh10(mn)

      logical ocean(mn)
!
!  local work arrays
!
!sun  include '../include/paramt.h' .. change im to mn
!     include '../include/paramt.h'
!
      real      phim(mn),phih(mn),zsz0(mn),rosfc(mn)                 &
              , wu(mn),wv(mn),wt(mn),wq(mn),zl(mn)                   &
              , q2(mn),p2(mn),q10(mn),p10(mn)

      integer   i,it1
      real      xkapa,g04,pihalf,xx,yy,btest,bl,al,zl2,zl10,xx10,phim10 &
                ,phih2,phih10
!
!
!
!          cccccccccccccccccccccccccccccccccccccccccccccc
!          c                                            c
!          c         similarity surface layer           c
!          c         paulson 1970, businger 1973        c
!          c                                            c
!          cccccccccccccccccccccccccccccccccccccccccccccc
!
      xkapa = r/cp
      do 10 i = 1, nxj
      rosfc(i) = 100.0*ps(i)/( r*ts(i))
  10  continue
!
!           iterate 3 times to match u*, t*, q*, and z/l
!
      do 900 it1 = 1, 3
!
!           charnock's relation for surface roughness over oceans
!
      do 100 i = 1, nxj
      if ( ocean(i) )  then
         z0(i) = 0.032 *ustar(i)*ustar(i)/g
         z0(i)= max(1.5e-5,min(z0(i),0.01))
      endif
! 100 continue
!     do 140 i=1,nxj
!
!           compute surface wind
!
      sfcw(i)= max( 1.0e-2,sqrt( ut(i)*ut(i)+vt(i)*vt(i)))
! 140 continue
!
!           compute z/l
!
!           limit minimum depth of surface layer
!           (saved in zsz0) for surface flux calculation
!
!
      g04 = g*0.4
!     do 200 i=1,nxj
      zsz0(i) = max(10.0,hgt(i))
      zl(i) = g04*zsz0(i)*tstar(i)*(1.+0.61*qstar(i)) /      &
              (tt(i)*ustar(i)*ustar(i))
      zl(i)= max(-2.0,min(zl(i),2.0))
! 200 continue
!
!           compute phi functions  (see businger 1973, or paulson 1970)
!
      pihalf = 3.141592654 * 0.5
!     do 240 i = 1, nxj
      if(zl(i) .lt. 0.0) then
         xx = exp(0.25*log(1.0 - 15.0*zl(i)))
         phim(i) = 2.0*log(0.5*(1.0+xx)) + log(0.5*(1.0+xx*xx))     &
                          - 2.0*atan(xx) + pihalf
         phih(i) = 2.0*log(0.5*(1.0 +sqrt(1.0 -9.0*zl(i))))
      else
         phim(i) = -4.7  * zl(i)
         phih(i) = -6.35 * zl(i)
      endif
! 240 continue
!
!     do 300 i = 1, nxj
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
!     do 340 i = 1, nxj
      yy    = dt*ustar(i) * ustar(i) / (2.0*hgt(i)*sfcw(i))
      btest = 0.5 + (yy-1.0)*0.5
      bl    = max (0.5, min (1.0, btest) )
      al    = 1.0 - bl
      wu(i) = ut(i) * (1.0-al*yy) / (1.0+bl*yy)
      wv(i) = vt(i) * (1.0-al*yy) / (1.0+bl*yy)
      xx    = 0.4*dt*ustar(i) / (0.74*2.0*hgt(i)*phih(i))
      wt(i) = ( tt(i)*(1.0-al*xx) + xx*tg(i) )   / (1.0+bl*xx)
      wq(i) = ( qt(i)*(1.0-al*xx) + xx*qsfc(i) ) / (1.0+bl*xx)
! 340 continue
!
!     update ustar, tstar and qstar
!
!     do 360 i = 1, nxj
! phon-ju
! add blending for stability
!
      wu(i)=0.5*(ut(i)+wu(i))
      wv(i)=0.5*(vt(i)+wv(i))
      wt(i)=0.5*(tt(i)+wt(i))
      wq(i)=0.5*(qt(i)+wq(i))
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
!     do 400 i = 1, nxj
      !hflux(i) = -rosfc(i)*cp*ustar(i)*tstar(i)
      !qflux(i) = -rosfc(i)*hltm*ustar(i)*qstar(i)
      hflux(i) = -rosfc(i)*cp*ustar(i)*max(tstar(i),-0.300) ! -3
      qflux(i) = -rosfc(i)*hltm*ustar(i)*min(max(qstar(i),-0.0003),0.0001)
! 400 continue
!
  100 continue
!
  900 continue
!
!ccc for output of t2, u10 and v10
!
      do i=1,nxj
      zl2=zl(i)*2./hgt(i)
      zl10=zl(i)*10./hgt(i)
      if(zl(i) .lt. 0.0) then
         xx10 = exp(0.25*log(1.0 - 15.0*zl10))
         phim10 = 2.0*log(0.5*(1.0+xx10)) + log(0.5*(1.0+xx10*xx10))  &
                          - 2.0*atan(xx10) + pihalf
         phih2 = 2.0*log(0.5*(1.0 +sqrt(1.0 -9.0*zl2)))
         phih10 = 2.0*log(0.5*(1.0 +sqrt(1.0 -9.0*zl10)))
      else
         phim10 = -4.7  * zl10
         phih2 = -6.35 * zl2
         phih10 = -6.35 * zl10
      endif
      phim10 =log( 10./z0(i)) - phim10
      phih2 =log( 2./z0(i)) - phih2
      phih10 =log( 10./z0(i)) - phih10
      phim10= max(phim10,1.0e-3)
      phih2= max(phih2,1.0e-3)
      phih10= max(phih10,1.0e-3)
!
      v10(i)=wv(i)*phim10/phim(i)
      u10(i)=wu(i)*phim10/phim(i)
!
      t2(i)=tg(i)*(1.-phih2/phih(i))+wt(i)*phih2/phih(i)
      t10(i)=tg(i)*(1.-phih10/phih(i))+wt(i)*phih10/phih(i)
!
      q2(i)=qsfc(i)*(1.-phih2/phih(i))+wq(i)*phih2/phih(i)
      q10(i)=qsfc(i)*(1.-phih10/phih(i))+wq(i)*phih10/phih(i)
!
! convert t2 from potential temp. to temp
      p2(i)=(ps(i)*100.-rosfc(i)*g*2.)/100.
      p10(i)=(ps(i)*100.-rosfc(i)*g*10.)/100.
      t2(i)=t2(i)*(p2(i)/1000.)**xkapa
      t10(i)=t10(i)*(p10(i)/1000.)**xkapa
      enddo
!
      call qsatq( nxj, t2, p2, wq)
!
      do i=1,nxj
       q2(i)=min(wq(i),q2(i))
       rh2(i)=max(min(q2(i)/wq(i),1.),0.0)
      enddo
!
      call qsatq( nxj, t10, p10, wq)
!
      do i=1,nxj
       q10(i)=min(wq(i),q10(i))
       rh10(i)=max(min(q10(i)/wq(i),1.),0.0)
      enddo
!
      return
      end
