      subroutine landsfc ( nxj,mn,kk                                    &
                        , dt,g,r,cp,hltm,tg,z0,ocean,ps,ts,hgt          &
                        , u,v,t,q,ut,vt,tt,qt,sfcw,ustar,tstar,qstar    &
                        , e,eps,hflux,qflux,zl,itypbl                   &
                        , as,ss,rld,stbo,pk,pk2,gfx                     &
                        ,km,snr,smc,stc,canopy,sigmaf,istyp,ice,land    &
                        ,tg3,rhscnpy,rhsmc,aim,bim,cim                  &
                        ,drain,snomt,zsoil,dth,ivegtyp,t2,rh2,u10,v10   &
                        ,smwlt,smref,tsat,dfkt,xktk,dfk)
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
!     zl : height(kk) over monin-obukhov scale height.
!          a stability indicator.                   (mn)
! itypbl : index for surface layer update, =0 stress b,c., =1 direct
!
! soil variables
!     km : dimension of soil layer
!    snr : snow depth                              (mn)     (mm)
!    smc : volumetric soil moisture content        (mn,km)
!    stc : soil temperature                        (mn,km)  (k)
! canopy : canopy moisture content                 (mn)     (mm)
! sigmaf : vegetation fraction                     (mn)
! istyp  : soil type(1-9)                          (mn)
!  ice   : logical array for sea ice  (mn)      (true = sea ice)
! land   : logical array for land     (mn)      (true = land  )
!  tg3   : deep layer soil temperature of annual mean (mn)
!rhscnpy : r.h.s. terms for canopy prediction (mn)   (mm/s)
! rhsmc  : r.h.s. terms for soil moisture prediction (mn,km)(s)
! aim    ; first column of tri-diagnal matrix for solving smc (mn,km)
! bim    ; second column of tri-diagonal matrix for solving smc (mn,km)
! cim    ; third column of tri-diagonal matrix for solving smc (mn,km)
! drain  : drainage from bottom soil layer  (mn)     (mm/s)
! snomt  : snow meltin rate (mn)   (mm/s)
! zsoil  : soil layer depth(mn,km)        (m)
! dth    : time step for snr updates      (s)
! ivegtyp: vegetation type(1-13)
! smwlt  : wiltling point of 9 soil types  (9)
! smref  : field capacity of 9 soil types  (9)
! tsat   : satuation point of 9 soil types (9)
! dfkt   : soil thermal diffusivity    (22,9)    m2/s
! xktk   : soil hydraulic conductivity (22,9)    m/s
!  dfk   : soil hydraulic diffusivity  (22,9)    m2/s
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
!     call landsfc ( mn,kk,dt,g,r,cp,hltm,tg,z0,ocean,ps,ts,hgt
!    1                  , u,v,t,q,ut,vt,tt,qt,sfcw,ustar,tstar,qstar
!    2                  , e,eps,hflux,qflux,zl,itypbl
!    3                  ,as,ss,rld,stbo,pk,pk2,gfx
!    4                  ,km,snr,smc,stc,canopy,sigmaf,istyp,ice,land
!    5                  ,tg3,rhscnpy,rhsmc,aim,bim,cim
!    6                  ,drain,snomt,zsoil,dth,ivegtyp,t2,u10,v10
!    7                  ,smwlt,smref,tsat,dfkt,xktk,dfk)
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
!     modified           jan. 2005
!
! 10. author
!      c. liou / f. wang / s. chang
!      jan. 2005 by f. wang
!      modify to f90 by C-H Lee and sort by River Chen in 2015
!
! 11. reference
! surface layer:
!
!       businger j.a.,1973: turbulent transfer in the atmospheric
!            surface layer.  workshop in meteorology, d. a. haugen,
!            ed., american meteorology socity, 67-98.
!
! land process:
!
!       Mahrt, L., and H. -L. Pan,1984: Atwo-layer model of soil hydrology.
!             Boundary Layer Meteorol., 29, 1-20.
!       Pan, H.-L., and L.Mahrt, 1987: Interaction between soil hydrology
!             and boundary-layer development. Boundary Layer Meteorol.,
!             38,185-202.
!#####################################################################
!
      use const, only : RTYPE
!
      implicit  none
!
! input & output variable
!
      integer   nxj,mn,kk,itypbl,km
      real      dt,g,r,cp,hltm,stbo,dth

      real      tg(mn),z0(mn),ps(mn),hgt(mn,kk),                   &
                u(mn,kk),v(mn,kk),t(mn,kk),                        &
                sfcw(mn),ustar(mn),                                &
                tstar(mn),qstar(mn),e(mn,kk),eps(mn,kk),hflux(mn), &
                qflux(mn),qsfc(mn),zl(mn),ts(mn),                  &
! for pm-evap.
                ss(mn),rld(mn),gfx(mn),                            &
                as(mn)
      real(kind=RTYPE) qt(mn,kk),q(mn,kk),tt(mn,kk),               &
                       ut(mn,kk),vt(mn,kk),pk(mn,kk),pk2(mn,kk)
      logical   ocean(mn)
!
      real      phim(mn),phih(mn),zsz0(mn),rosfc(mn),              &
                wu(mn),wv(mn),wt(mn),wq(mn),                       &
                wt1(mn),wt2(mn)
! for landsfc-layer
      real      z0max(mn),ztmax(mn),dtv(mn),adtv(mn),rb(mn)
      real      fm(mn),fh(mn),hlinf(mn),olinf(mn)
      real      hl1(mn),pm(mn),ph(mn),hl110(mn),hl12(mn)
      real      restar(mn),rat(mn),                                &
! for pm-evap.
                p1(mn),qs1(mn),t12(mn),t14(mn),ta(mn),             &
                rch(mn),rnet(mn),rsmall(mn),delta(mn),             &
                terma(mn),termb(mn),termc(mn),rs(mn),              &
                twork(mn),rsi(mn),qs0(mn)
! for soil input and output
       real     zsoil(mn,km),sigmaf(mn),tg3(mn)
       real     snr(mn),smc(mn,km)
       real     canopy(mn),rhscnpy(mn),rhsmc(mn,km),               &
                aim(mn,km),bim(mn,km),cim(mn,km)
       real     drain(mn),stc(mn,km),snomt(mn)
       logical  land(mn),ice(mn),flagsnw(mn),flag(mn)
! local
       real     stsoil(mn,km),tsurf(mn),smcz(mn)
       real     partlnd(mn),factsnw(mn)
       real     evap(mn),ep(mn),etp(mn),edir(mn),ec(mn),           &
                et(mn,km),snoevp(mn)
       real     tref(mn),twilt(mn),dft0(mn),                       &
                canfac(mn),etpfac(mn),dft1(mn),dft2(mn)
       real     gx(mn),fx(mn),rhstc(mn,km)
       real     dew(mn)
       real     dmdz(mn),ddz(mn),dmdz2(mn),ddz2(mn),               &
                df1(mn),df2(mn)
       real     xx(mn),yy(mn),zz(mn)
       real     dtdz1(mn),dtdz2(mn),hcpct(mn)
       real     ai(mn,km),bi(mn,km),ci(mn,km)
       integer  istyp(mn),ivegtyp(mn)
       real     xkt1(mn),xkt2(mn),smcdry(mn),beta(mn)
! landpack
!--------------------------------------------
! soil function
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
!
      integer,  parameter :: ntype=9,ngrid=22
      real      tsat(ntype),dfkt(ngrid,ntype),           &
                xktk(ngrid,ntype),dfk(ngrid,ntype)
      real      smref(ntype),smwlt(ntype)
!
!     data smref/.283,.387,.412,.312,.338,.382,.315,.329,.283/
!     data smwlt/.029,.119,.139,.047,.010,.103,.069,.066,.029/
!
      real      smmax(ntype),smdry(ntype)
!
      data smmax/.421,.464,.468,.434,.406,.465,.404,.439,.421/
      data smdry/.07,.14,.22,.08,.18,.16,.12,.10,.07/
!------------------------------------------------------------
! vegtype(13)
!
!  the 13 vegetation types are:
!
!  1  ...  broadleave-evergreen trees (tropical forest)
!  2  ...  broadleave-deciduous trees
!  3  ...  broadleave and needle leave trees (mixed forest)
!  4  ...  needleleave-evergreen trees
!  5  ...  needleleave-deciduous trees (larch)
!  6  ...  broadleave trees with groundcover (savanna)
!  7  ...  groundcover only (perenial)
!  8  ...  broadleave shrubs with perenial groundcover
!  9  ...  broadleave shrubs with bare soil
! 10  ...  dwarf trees and shrubs with ground cover (trunda)
! 11  ...  bare soil
! 12  ...  cultivations (use parameters from type 7)
! 13  ...  glacial
!
      integer,  parameter :: nvtype=13
      real      rsmax(nvtype),rgl(nvtype),rsmin(nvtype),hs(nvtype)
      data rsmax/13*5000./
      data rsmin/150.,100.,125.,150.,100.,70.,40.,                  &
                 300.,400.,150.,999.,40.,999./
      data rgl/5*30.,65.,4*100.,999.,100.,999./
      data hs/41.69,54.53,51.93,47.35,47.35,54.53,36.35,            &
              3*42.00,999.,36.35,999./
!
      logical   ncepsfc
      data      ncepsfc/.false./
!-----------------------------------------------------------------
! landsfc parameter
      real, parameter :: alpha=5.,a0=-3.975,a1=12.32,b1=-7.755,b2=6.041
      real, parameter :: a0p=-7.941,a1p=24.75,b1p=-8.705,b2p=7.899,vis=1.4e-5
      real, parameter :: aa1=-1.076,bb1=.7045,cc1=-.05808
      real, parameter :: bb2=-.1954,cc2=.009999
      real, parameter :: ca=.4,topt=298.
!
      real, parameter :: cice=1880.*917.
      real, parameter :: dfsnow=0.31
      real, parameter :: scanop=0.5
      real, parameter :: rhoh2o=1000.
      real, parameter :: hfus=3.3358e+5
      real, parameter :: ch2o=4.2e6,csoil=1.26e6
      real, parameter :: zbot=-3.
      real, parameter :: tgice=271.2
      real, parameter :: ctfil1=0.5,ctfil2=0.5
      real, parameter :: cfactr=0.5
      real, parameter :: t0c=273.16
!
!cc for output of t2, u10 and v10 ( 2003/01/20 )
      real  u10(mn),v10(mn),t2(mn),q2(mn),rh2(mn)

      integer i, k, it1
      real    g04,  xkapa, pihalf, xx1,   yy1,    btest,  bl,     al,   zl2,  &
              zl10, xx10,  phim10, phih2, hl0inf, hltinf, aa,     aa0,  bb,   &
              bb0,  fms,   fhs,    hl0,   hlt,    hl1x,   xtsoil, tflx,       &
              rcs,  rct,   rcq,    ff,    rss,    bfact,  cc,     hpbl, ghpbl,&
              zlx, wstar,  funcdf, funckt,p2

!
!          cccccccccccccccccccccccccccccccccccccccccccccc
!          c                                            c
!          c         similarity surface layer           c
!          c         paulson 1970, businger 1973        c
!          c                                            c
!          cccccccccccccccccccccccccccccccccccccccccccccc
!
      g04 = g*0.4
      xkapa = r/cp
      do 10 i = 1, nxj
      rosfc(i) = 100.0*ps(i)/( r*ts(i))
      p1(i)=ps(i)-rosfc(i)*g*hgt(i,kk)/100.
! define stamotal resistance 60 s/m
!cc   sigmaf(i)=0.7
!cc   rs(i)=60.
      if(ivegtyp(i).gt.0.)rs(i)=rsmin(ivegtyp(i))
!
      factsnw(i)=10.
      if(ice(i))factsnw(i)=3.
!
      tsurf(i)=tg(i)
  10  continue
!
      do i=1,nxj
!
! snr(snow depth) in unit mm
!
      flagsnw(i)=.false.
      if(snr(i).gt.1. .or.ice(i))rs(i)=0.
      if(snr(i).gt.1.)then
      flagsnw(i)=.true.
      as(i)=0.7
      endif
!
! partlnd is percent of land in one grid
!
      partlnd(i)=1.
      if(snr(i).gt.0. .and. snr(i).le.1.)then
      partlnd(i)=1.-snr(i)/1.
      endif
!
      enddo
!
! define soil layer depth
!
      do i=1,nxj
      if(ocean(i))then
       zsoil(i,1)=0.
      elseif(land(i))then
       zsoil(i,1)=-.1
      else
       zsoil(i,1)=-3./km
      endif
      enddo
!
      do k=2,km
       do i=1,nxj
       if(ocean(i))then
        zsoil(i,k)=0.
       elseif(land(i))then
        zsoil(i,k)=zsoil(i,k-1)+(-2.-zsoil(i,1))/(km-1)
       else
        zsoil(i,k)=-3.*float(k)/float(km)
       endif
       enddo
      enddo
!
       do k=1,km
        do i=1,nxj
        aim(i,k)=0.
        bim(i,k)=1.
        cim(i,k)=0.
        rhsmc(i,k)=0.
        et(i,k)=0.
        stsoil(i,k)=stc(i,k)
        enddo
       enddo
!
      do i = 1,nxj
        edir(i) = 0.
        ec(i) = 0.
        evap(i) = 0.
        ep(i) = 0.
        etp(i) = 0.
        snomt(i) = 0.
        gfx(i) = 0.
        rhscnpy(i) = 0.
        fx(i) = 0.
        etpfac(i) = 0.
        canfac(i) = 0.
        drain(i) =0.
      enddo
!
      do i=1,nxj
      twork(i)=tt(i,kk)*pk(i,kk)
      enddo
      call qsatq(nxj,twork,p1,qs1)
      call qsatq(nxj,tsurf,ps,qs0)
!
! define surfce specific humidity
      do i=1,nxj
      if(ocean(i).or. istyp(i).eq.0)then
      qsfc(i)=qs0(i)
      else
      beta(i)=(smc(i,1)-smwlt(istyp(i)))/(smref(istyp(i))-smwlt(istyp(i)))
      qsfc(i)=beta(i)*qs0(i)+(1.-beta(i))*qt(i,kk)
      qsfc(i)=min(qs0(i),qsfc(i) )
      endif
      enddo
!
!
!           charnock's relation for surface roughness over oceans
!
      do 100 i = 1, nxj
      if ( ocean(i) )  then
         z0(i) = 0.032 *ustar(i)*ustar(i)/g
         z0(i)= max(1.5e-5,min(z0(i),0.01))
      endif
!
!
!           compute surface wind
!
      sfcw(i)= max( 1.0e-2,sqrt( ut(i,kk)*ut(i,kk)+vt(i,kk)*vt(i,kk)))
  100 continue
!
!           compute surface wind
!
!------------------------------------------------------------
! for ocean
!
!           iterate 3 times to match u*, t*, q*, and z/l
!
!
      do  it1 = 1, 3
!
!           charnock's relation for surface roughness over oceans
!
      do  i = 1, nxj
      if ( ocean(i) )  then
!
         z0(i) = 0.032 *ustar(i)*ustar(i)/g
         z0(i)= max(1.5e-5,min(z0(i),0.01))
      endif
!
!
!
!           compute z/l
!
!           limit minimum depth of surface layer
!           (saved in zsz0) for surface flux calculation
!
!
      g04 = g*0.4
      zsz0(i) = max(10.0,hgt(i,kk))
      zl(i) = g04*zsz0(i)*tstar(i)*(1.+0.61*qstar(i)) /   &
              (tt(i,kk)*ustar(i)*ustar(i))
      zl(i)= max(-2.0,min(zl(i),2.0))
!
!
!           compute phi functions  (see businger 1973, or paulson 1970)
!
      pihalf = 3.141592654 * 0.5
!
      if(zl(i) .lt. 0.0) then
!ibm        xx1 = exp(0.25*log(1.0 - 15.0*zl(i)))
!ibm        phim(i) = 2.0*log(0.5*(1.0+xx1)) + log(0.5*(1.0+xx1*xx1))
!ibm    1                    - 2.0*atan(xx1) + pihalf
         xx1 = sqrt(sqrt(1.0 - 15.0*zl(i)))
         phim(i) = log(((0.5*(1.0+xx1))**2)*(0.5*(1.0+xx1*xx1)))  &
                          - 2.0*atan(xx1) + pihalf
         phih(i) = 2.0*log(0.5*(1.0 +sqrt(1.0 -9.0*zl(i))))
      else
         phim(i) = -4.7  * zl(i)
         phih(i) = -6.35 * zl(i)
      endif
!
      zsz0(i) = log(zsz0(i)/z0(i))
      phim(i) = zsz0(i) - phim(i)
      phih(i) = zsz0(i) - phih(i)
      phim(i)= max(phim(i),1.0e-3)
      phih(i)= max(phih(i),1.0e-3)
!
!     compute new ut, vt, tt, and qt for ustar, tstar, qstar updates
!     using crank-nicolson (bl=0.5) time scheme for weak forcing or
!        or backward (bl=1.) time scheme for strong forcing
!
      yy1   = dt*ustar(i) * ustar(i) / (2.0*hgt(i,kk)*sfcw(i))
      btest = 0.5 + (yy1-1.0)*0.5
      bl    = max (0.5, min (1.0, btest) )
      al    = 1.0 - bl
      wu(i) = ut(i,kk) * (1.0-al*yy1) / (1.0+bl*yy1)
      wv(i) = vt(i,kk) * (1.0-al*yy1) / (1.0+bl*yy1)
      xx1    = 0.4*dt*ustar(i) / (0.74*2.0*hgt(i,kk)*phih(i))
      wt(i) = ( tt(i,kk)*(1.0-al*xx1)+xx1*tg(i)/pk2(i,kk) )/(1.0+bl*xx1)
      wq(i) = ( qt(i,kk)*(1.0-al*xx1) + xx1*qsfc(i) ) / (1.0+bl*xx1)
!
!     update ustar, tstar and qstar
!
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
      tstar(i) = 0.4*        (wt(i)-tg(i)/pk2(i,kk))  /(0.74*phih(i))
      qstar(i) = 0.4*        (wq(i)-qsfc(i))/(0.74*phih(i))
      ustar(i)= max(ustar(i),1.0e-7)
      ustar(i)= min(ustar(i),10.000)
      tstar(i)= min(tstar(i),50.00)
      qstar(i)= min(qstar(i),5.000)
!
      enddo
      enddo
!---- iteration end
!
!ccc for output of t2, u10 and v10
!
      do i=1,nxj
!
      zl2=zl(i)*2./hgt(i,kk)
      zl10=zl(i)*10./hgt(i,kk)
      if(zl(i) .lt. 0.0) then
!ibm        xx10 = exp(0.25*log(1.0 - 15.0*zl10))
!ibm        phim10 = 2.0*log(0.5*(1.0+xx10)) + log(0.5*(1.0+xx10*xx10))
!ibm    1                    - 2.0*atan(xx10) + pihalf
         xx10 = sqrt(sqrt(1.0 - 15.0*zl10))
         phim10 = log(((0.5*(1.0+xx10))**2)*(0.5*(1.0+xx10*xx10)))   &
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
      q2(i)=qsfc(i)*(1.-phih2/phih(i))+wq(i)*phih2/phih(i)
      t2(i)=tg(i)/pk2(i,kk)*(1.-phih2/phih(i))+wt(i)*phih2/phih(i)
!
! convert t2 from potential temp. to temp
      p2=(ps(i)*100.-rosfc(i)*g*2.)/100.
      wt1(i)=p2
      t2(i)=t2(i)*(p2/1000.)**xkapa
!f    t2(i)=t2(i)*exp(xkapa*log(p2/1000.))

! do transfer from cwb-sfc to match ncep-sfc
      phih(i)=0.4*0.4/0.74/phim(i)/phih(i)
      phim(i)=0.4*0.4/phim(i)/phim(i)
!
      enddo
!
      call qsatq( nxj, t2, wt1, wt2)
!
      do i=1,nxj
       q2(i)=min(wt2(i),q2(i))
       rh2(i)=q2(i)/wt2(i)
      enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! surface layer for land or ice
!
!     refer to ncep
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!---- ncepsfc
      if(ncepsfc)then
!
      do  i=1,nxj
      if(.not. ocean(i))then
      zsz0(i) = max(10.0,hgt(i,kk))
      wt1(i)=tt(i,kk)*pk2(i,kk)*(1.+0.61*qt(i,kk))
      wt2(i)=tg(i)*(1.+0.61*qsfc(i))
      rb(i) = g*zsz0(i)*( wt1(i)-wt2(i))/(0.5*(wt1(i)+wt2(i)))  &
              /(sfcw(i)*sfcw(i))
!
!
!  compute stability indices (rb and hlinf)
!
        z0max(i) = min(z0(i),hgt(i,kk))
        ztmax(i) = z0max(i)
        if(ocean(i)) then
          restar(i) = ustar(i) * z0max(i) / vis
          restar(i) = max(restar(i),.000001)
!         restar(i) = alog(restar(i))
!         restar(i) = min(restar(i),5.)
!         restar(i) = max(restar(i),-5.)
!         rat(i) = aa1 + bb1 * restar(i) + cc1 * restar(i) ** 2
!         rat(i) = rat(i) / (1. + bb2 * restar(i)
!    &                       + cc2 * restar(i) ** 2)
! rat taken from zeng,zhao and dickinson 1997
!ibm      rat(i) = 2.67*restar(i)**.25-2.57
          rat(i) = 2.67*sqrt(sqrt(restar(i)))-2.57
          rat(i) = min(rat(i),7.)
          ztmax(i) = z0max(i) * exp(-rat(i))
        endif
!
        dtv(i) = wt1(i) - wt2(i)
        adtv(i) = abs(dtv(i))
        adtv(i) = max(adtv(i),.001)
        dtv(i) = sign(1.,dtv(i)) * adtv(i)
        rb(i) = g * dtv(i) * hgt(i,kk)/ (.5 * (wt1(i) + wt2(i))  &
                * sfcw(i) * sfcw(i))
        rb(i) = max(rb(i),-5000.)
        fm(i) = log((z0max(i)+hgt(i,kk)) / z0max(i))
        fh(i) = log((ztmax(i)+hgt(i,kk)) / ztmax(i))
        hlinf(i) = rb(i) * fm(i) * fm(i) / fh(i)
!
!
!  stable case
!
!
        if(dtv(i).ge.0.) then
          hl1(i) = hlinf(i)
        endif
        if(dtv(i).ge.0..and.hlinf(i).gt..25) then
          hl0inf = z0max(i) * hlinf(i) / hgt(i,kk)
          hltinf = ztmax(i) * hlinf(i) / hgt(i,kk)
          aa = sqrt(1. + 4. * alpha * hlinf(i))
          aa0 = sqrt(1. + 4. * alpha * hl0inf)
          bb = aa
          bb0 = sqrt(1. + 4. * alpha * hltinf)
          pm(i) = aa0 - aa + log((aa + 1.) / (aa0 + 1.))
          ph(i) = bb0 - bb + log((bb + 1.) / (bb0 + 1.))
          fms = fm(i) - pm(i)
          fhs = fh(i) - ph(i)
          hl1(i) = fms * fms * rb(i) / fhs
        endif
!
!  second iteration
!
!
        if(dtv(i).ge.0.) then
          hl0 = z0max(i) * hl1(i) /hgt(i,kk)
          hlt = ztmax(i) * hl1(i) /hgt(i,kk)
          aa = sqrt(1. + 4. * alpha * hl1(i))
          aa0 = sqrt(1. + 4. * alpha * hl0)
          bb = aa
          bb0 = sqrt(1. + 4. * alpha * hlt)
          pm(i) = aa0 - aa + log((aa + 1.) / (aa0 + 1.))
          ph(i) = bb0 - bb + log((bb + 1.) / (bb0 + 1.))
!
!         hl110(i) = hl1(i) * 10. /hgt(i,kk)
!         aa = sqrt(1. + 4. * alpha * hl110(i))
!         pm10(i) = aa0 - aa + log((aa + 1.) / (aa0 + 1.))
!         hl12(i) = hl1(i) * 2. /hgt(i,kk)
!         aa = sqrt(1. + 4. * alpha * hl12(i))
!         bb = sqrt(1. + 4. * alpha * hl12(i))
!         ph2(i) = bb0 - bb + log((bb + 1.) / (bb0 + 1.))
        endif
!
!  unstable case
!
!
!  check for unphysical obukhov length
!
!
        if(dtv(i).lt.0.) then
          olinf(i) = hgt(i,kk) / hlinf(i)
          if(abs(olinf(i)).le.50. * z0max(i)) then
            hlinf(i) = -hgt(i,kk) / (50. * z0max(i))
          endif
        endif
!
!
!  get pm and ph
!
!
        if(dtv(i).lt.0..and.hlinf(i).ge.-.5) then
          hl1(i) = hlinf(i)
          pm(i) = (a0 + a1 * hl1(i)) * hl1(i)                    &
                  / (1. + b1 * hl1(i) + b2 * hl1(i) * hl1(i))
          ph(i) = (a0p + a1p * hl1(i)) * hl1(i)                  &
                  / (1. + b1p * hl1(i) + b2p * hl1(i) * hl1(i))
        endif
        if(dtv(i).lt.0.and.hlinf(i).lt.-.5) then
          hl1(i) = -hlinf(i)
!ibm         pm(i) = log(hl1(i)) + 2. * hl1(i) ** (-.25) - .8776
!ibm         ph(i) = log(hl1(i)) + .5 * hl1(i) ** (-.5) + 1.386
          hl1x=sqrt(hl1(i))
          pm(i) = log(hl1(i)) + 2./(sqrt(hl1x)) - .8776
          ph(i) = log(hl1(i)) + .5/hl1x + 1.386
        endif
!
!  finish the exchange coefficient computation to provide fm and fh
!
!
        fm(i) = fm(i) - pm(i)
        fh(i) = fh(i) - ph(i)
!       fm10(i) = fm10(i) - pm10(i)
!       fh2(i) = fh2(i) - ph2(i)
!
        phim(i) = ca * ca / (fm(i) * fm(i))
        phih(i) = ca * ca / (fm(i) * fh(i))
        ustar(i) = sqrt(phim(i) * sfcw(i) * sfcw(i))
!
       endif
      enddo
!
      do i=1,nxj
      if(.not. ocean(i))then
      yy1   =dt*ustar(i)*ustar(i) / (2.0*hgt(i,kk)*sfcw(i))
      btest = 0.5 + (yy1-1.0)*0.5
      bl    = max (0.5, min (1.0, btest) )
      al    = 1.0 - bl

      wu(i)=ut(i,kk) * (1.-al*yy1) / (1.+bl*yy1)
      wv(i)=vt(i,kk) * (1.-al*yy1) / (1.+bl*yy1)
      endif
      enddo
!
      do i=1,nxj
      if(.not. ocean(i))then
      wu(i) = 0.5*(wu(i)+ut(i,kk))
      wv(i) = 0.5*(wv(i)+vt(i,kk))
      sfcw(i) =max( 1.0e-2,sqrt( wu(i)*wu(i)+wv(i)*wv(i) ))
      ustar(i)=sqrt(phim(i)*sfcw(i)*sfcw(i))
      endif
      enddo
!
      endif
!---- ncepsfc end
!
!---------------------------------------------------------
!  define roughness and evaporation  over ocean
!
      do i=1,nxj
      if(ocean(i))then
       z0(i) = 0.032 *ustar(i)*ustar(i)/g
       z0(i)= max(1.5e-5,min(z0(i),0.01))
!
       evap(i)=rosfc(i)*hltm*phih(i)*sfcw(i)*(qsfc(i)-wq(i))
      endif

      enddo
!
! ...................................
! prepare for penman-monteith method. to calculate potential evaporation
! ...................................
! change wt from potential tmp. to temperature
!
!
! compute soil/snow/ice heat flux(gfx)
!
      do i=1,nxj
      gfx(i)=0.
      if(land(i))then
       smcz(i)=0.5*(smc(i,1)+0.2)
       dft0(i)=xtsoil(smcz(i),istyp(i),tsat,dfkt)
      elseif(ice(i))then
       dft0(i)=2.2
      endif
      enddo
!
       do i=1,nxj
       if(.not. ocean(i))then
        if(flagsnw(i))then
         tflx=min(twork(i),tsurf(i))
         gfx(i)=-dfsnow*(tflx-stsoil(i,1))/             &
                (factsnw(i)*max(snr(i)/1000.,0.001))
        else
         gfx(i)=-dft0(i)*(tsurf(i)-stsoil(i,1))/(-0.5*zsoil(i,1))
        endif
       gfx(i)=max(gfx(i),-200.)
       gfx(i)=min(gfx(i),200.)
       endif
       enddo
!
! calculate potential evaporation: ep(w/m2)
!
      do i=1,nxj
      if(.not.ocean(i))then
      t12(i)=twork(i)*twork(i)
      t14(i)=t12(i)*t12(i)
!
      ta(i)=tt(i,kk)*pk2(i,kk)
!
      rch(i)=rosfc(i)*cp*phih(i)*sfcw(i)
!
! rnet=ss+ld-stbo*t14-gfx-rho*cp*ch*v*(t1-ta)
!
!     gfx(i)=-rosfc(i)*cp*ustar(i)*tstar(i)*0.1
      rnet(i)=ss(i)+rld(i)-stbo*t14(i)+gfx(i)-rch(i)*(twork(i)-ta(i))
!
! rsmall=4*sigma*t1**3/rch +1
!
      rsmall(i)=4.*stbo*twork(i)*t12(i)/rch(i) +1.
!
! delta= hltm/cp *dqs/dt
!
      delta(i)=hltm/cp*(0.622* 2.5e+6 *qs1(i))/(r*t12(i))
!
      terma(i)=rnet(i)*delta(i)
      termb(i)=rsmall(i)*hltm*(rch(i)/cp)*(qs1(i)-qt(i,kk))
!
      ep(i)=termb(i)+terma(i)
      ep(i)=ep(i)/(rsmall(i)+delta(i))
      endif
      enddo
!
! actual evaporation over land in three parts: edir, et, and ec
!
! edir :direct evaporation from soil
!
       do i=1,nxj
       flag(i)=land(i) .and. (ep(i).gt.0.)
       enddo
       do i=1,nxj
        if(flag(i))then
        df1(i)=funcdf(smc(i,1),istyp(i),tsat,dfk)
        xkt1(i)=funckt(smc(i,1),istyp(i),tsat,xktk)
        endif
        if(flag(i).and.stc(i,1).lt.t0c)then
        df1(i)=0.
        xkt1(i)=0.
        endif
!
        if(flag(i))then
!       tref(i)=0.75*tsat(istyp(i))
!       twilt(i)=twlt(istyp(i))
!       fx(i)=-2.*df1(i)*(smc(i,1)-0.23)/zsoil(i,1)-xkt1(i)
        tref(i)=smref(istyp(i))
        twilt(i)=smwlt(istyp(i))
        smcdry(i)=smdry(istyp(i))
        fx(i)=-2.*df1(i)*(smc(i,1)-smcdry(i))/zsoil(i,1)-xkt1(i)
        fx(i)=min(fx(i),ep(i)/hltm)
        fx(i)=max(fx(i),0.)
!
!       edir(i)=fx(i)*(1-sigmaf(i))*partlnd(i)
!---another dir evap =(1-sigmaf)*ep*beta (mahfouf and noilhan 1991)
!
        edir(i)=(1-sigmaf(i))*partlnd(i)*beta(i)*ep(i)/hltm

        endif
       enddo
!
! et: transpiration from all levels of the soil
!
       do i=1,nxj
        if(flag(i))then
!
        rcs = 1.
        rct = 1.
        rcq = 1.
!
!  resistance due to solar effect
!
        ff  = .55* 2.* ss(i)/rgl(ivegtyp(i))
        rcs = (ff + rs(i)/rsmax(ivegtyp(i)))/(1.+ff)
        rcs = max(rcs, 0.0001)
!
!  resistance due to thermal effect
!
!         rct = 1. - .0016 * (topt - ta(i)) ** 2
!         rct = max(rct,.0001)
!
!  resistance due to humidity
!
!         rcq = 1. / (1. + hs(ivegtyp(i)) * (qs1(i) - qt(i,kk)))
!         rcq = max(rcq,.0001)
!
!  compute resistance without the effect of soil moisture
!
         rs(i) = rs(i) / (rcs * rct * rcq)
        endif
!
        if(flag(i))then
        canfac(i)=(canopy(i)/scanop)**cfactr
        etpfac(i)=sigmaf(i)*(1.-canfac(i))/hltm
        gx(i)=(smc(i,1)-twilt(i))/(tref(i)-twilt(i))
        gx(i)=max(gx(i),0.)
        gx(i)=min(gx(i),1.)
!
! resistance due to soil moisture deficit
!
!       rss = 1.0
        rss =gx(i)*(zsoil(i,1)/zsoil(i,km))
        rss =max(rss,0.0001)
        rsi(i)=rs(i)/rss
!
! transpiration la monteith
!
        termc(i)=delta(i)+rsmall(i)*(1.+rsi(i)*sfcw(i)*phih(i))
        etp(i)=(termb(i)+terma(i))/termc(i)
        et(i,1)=etp(i)*etpfac(i)*partlnd(i)
!
        endif
       enddo
!
       do k=2,km
        do i=1,nxj
        if(flag(i))then
        gx(i)=(smc(i,k)-twilt(i))/(tref(i)-twilt(i))
        gx(i)=max(gx(i),0.)
        gx(i)=min(gx(i),1.)
!
! resistance due to soil moisyure deficit
!
!       rss = 1.
        rss =gx(i)*(zsoil(i,k)-zsoil(i,k-1))/zsoil(i,km)
        rss =max(rss,1.e-6)
        rsi(i)=rs(i)/rss
!
! transpiration la monteith
!
        termc(i)=delta(i)+rsmall(i)*(1.+rsi(i)*sfcw(i)*phih(i))
        etp(i)=(termb(i)+terma(i))/termc(i)
        et(i,k)=etp(i)*etpfac(i)*partlnd(i)
        endif
        enddo
       enddo
!
! ec: canopy re-evaporation
!
       do i=1,nxj
        if(flag(i))then
        ec(i)=sigmaf(i)*canfac(i)*ep(i)/hltm
        ec(i)=ec(i)*partlnd(i)
        ec(i)=min(ec(i),canopy(i)/dt)
        endif
       enddo
!
! sum up total evaporation
!
       do i=1,nxj
       if(flag(i))then
       evap(i)=edir(i)+ec(i)
       endif
       enddo
!
       do k=1,km
        do i=1,nxj
        if(flag(i))then
        evap(i)=evap(i)+et(i,k)
        endif
        enddo
       enddo
!
! return unit from kg m-2 s-1 to watts m-2
!
       do i=1,nxj
        if(flag(i))then
        evap(i)=min(evap(i)*hltm,ep(i))
        endif
       enddo
!
! evaporation over bare sea ice
       do i=1,nxj
        if(ice(i))then
        evap(i)=partlnd(i)*ep(i)
        endif
       enddo
!
! treat downward moisture flux situation
        do i=1,nxj
         flag(i)=.not.ocean(i).and.ep(i).le.0.
         dew(i)=0.
        enddo
        do i=1,nxj
         if(flag(i))then
         dew(i)=-ep(i)*dth/(hltm*rhoh2o)
         evap(i)=ep(i)*partlnd(i)
         dew(i)=dew(i)*partlnd(i)
         endif
        enddo
!
! snow covered land and sea ice
!
        do i=1,nxj
         flag(i)=.not.ocean(i).and.snr(i).gt.0.
        enddo
!
! snow evaporation or sublimation
! affect snow depth and evap,tsurf
!
       do i=1,nxj
        if(flag(i))then
        bfact=snr(i)/(1000.*dth*ep(i)/(hltm*rhoh2o))
        bfact=min(bfact,1.)
!snow evaporation
        if(ep(i).le.0.)bfact=1.
        if(snr(i).le.1.)then
         snoevp(i)=(1.-partlnd(i))*bfact*ep(i)
         evap(i)=evap(i)+snoevp(i)
        else
         snoevp(i)=bfact*ep(i)
         evap(i)=snoevp(i)
        endif
!
        tsurf(i)=twork(i)+(rnet(i)-gfx(i)-dfsnow*(twork(i)         &
               -stsoil(i,1))/(factsnw(i)*max(snr(i)/1000.,0.001))  &
               -snoevp(i))/(rsmall(i)*rch(i)+dfsnow/               &
               (factsnw(i)*max(snr(i)/1000.,0.001) ))
!
        snr(i)=snr(i)-snoevp(i)*dth/(rhoh2o*hltm)*1000.
        snr(i)=max(snr(i),1.0e-6)
        endif
       enddo
!snow melt rate
        do i=1,nxj
         flag(i)=.not.ocean(i).and.snr(i).gt.0.
        enddo
        do i=1,nxj
         if(flag(i).and.tsurf(i).gt.273.15)then
         snomt(i)=rch(i)*rsmall(i)*(tsurf(i)-t0c) &
                  /(rhoh2o*hfus)*1000.
         snomt(i)=min(snomt(i),snr(i)/dth)
         snr(i)=snr(i)-snomt(i)*dth
         snr(i)=max(snr(i),1.0e-6)
         tsurf(i)=max(273.15,tsurf(i)-hfus*rhoh2o      &
                  *snomt(i)/(rch(i)*rsmall(i)*1000.) )
! re-evaluate evap because of snow melt
!        call qsatq(1,tsurf(i),ps(i),qss)
!        evap(i)=rosfc(i)*hltm*phih(i)*sfcw(i)*
!    +           (qss-qt(i,kk))
         endif
       enddo
!....................
!
!     update ustar, tstar and qstar
!
      do 350 i = 1, nxj
      qstar(i)=-evap(i)/(rosfc(i)*hltm*ustar(i))
  350 continue
!
      do 360 i=1,nxj
      ustar(i)= max(ustar(i),1.0e-7)
      ustar(i)= min(ustar(i),10.000)
      tstar(i)= min(tstar(i),50.00)
      qstar(i)= min(qstar(i),5.000)
  360 continue
!
!
! prepare rhsmc except precipitation term
!
      do i=1,nxj
      flag(i)=land(i)
      enddo
!
!
      do i=1,nxj
        if(flag(i))then
        smcz(i)=max(smc(i,1),smc(i,2))
        df1(i)=funcdf(smcz(i),istyp(i),tsat,dfk)
        xkt1(i)=funckt(smcz(i),istyp(i),tsat,xktk)
        endif
        if(flag(i).and.stc(i,1).lt.t0c)then
        df1(i)=0.
        xkt1(i)=0.
        endif
      if(flag(i))then
        rhscnpy(i)=-ec(i)+sigmaf(i)*rhoh2o*dew(i)/dth
!
        dmdz(i)=(smc(i,1)-smc(i,2))/(-.5*zsoil(i,2))
        rhsmc(i,1)=(df1(i)*dmdz(i)+xkt1(i)+(edir(i)   &
              +et(i,1)) )/(zsoil(i,1)*rhoh2o)
!
! add new 91version
        rhsmc(i,1)=rhsmc(i,1)-(1.-sigmaf(i))*dew(i)/(zsoil(i,1)*dth)
!
        ddz(i)=1./(-.5*zsoil(i,2))
        aim(i,1)=0.
        bim(i,1)=df1(i)*ddz(i)/(-zsoil(i,1)*rhoh2o)
        cim(i,1)=-bim(i,1)
      endif
      enddo
!
      do k=2,km
       if(k.lt.km)then
        do i=1,nxj
        if(flag(i))then
        smcz(i)=max(smc(i,k),smc(i,k+1))
        df2(i)=funcdf(smcz(i),istyp(i),tsat,dfk)
        xkt2(i)=funckt(smcz(i),istyp(i),tsat,xktk)
        endif
        if(flag(i).and.stc(i,k).lt.t0c)then
        df2(i)=0.
        xkt2(i)=0.
        endif

        if(flag(i))then
        dmdz2(i)=(smc(i,k)-smc(i,k+1))/(0.5*         &
                (zsoil(i,k-1)-zsoil(i,k+1)))
!       smcz(i)=max(smc(i,k),smc(i,k+1))
!       df2(i)=funcdf(smcz(i),istyp(i),tsat,dfk)
!       xkt2(i)=funckt(smcz(i),istyp(i),tsat,xktk)
        rhsmc(i,k)=(df2(i)*dmdz2(i)+xkt2(i)-         &
               df1(i)*dmdz(i)-xkt1(i)+et(i,k))       &
               /(rhoh2o*(zsoil(i,k)-zsoil(i,k-1)))
        ddz2(i)=2./(zsoil(i,k-1)-zsoil(i,k+1))
        cim(i,k)=-df2(i)*ddz2(i)/(rhoh2o*            &
                 (zsoil(i,k-1)-zsoil(i,k)))
        endif
       enddo
       else        !k=km
        do i=1,nxj
           if(flag(i))then
           xkt2(i)=funckt(smc(i,k),istyp(i),tsat,xktk)
           rhsmc(i,k)=(xkt2(i)-df1(i)*dmdz(i)-xkt1(i)      &
              +et(i,k))/(rhoh2o*(zsoil(i,k)-zsoil(i,k-1)))
           drain(i)=xkt2(i)
           cim(i,k)=0.
           endif
        enddo
       endif
!
       do i=1,nxj
        if(flag(i))then
          aim(i,k)=-df1(i)*ddz(i)/(rhoh2o*(zsoil(i,k-1)-zsoil(i,k)))
          bim(i,k)=-(aim(i,k)+cim(i,k))
!
          if(k .lt. km)then
          df1(i)=df2(i)
          xkt1(i)=xkt2(i)
          dmdz(i)=dmdz2(i)
          ddz(i)=ddz2(i)
          endif
        endif
       enddo
      enddo
!
! update stsoil and tsurf
!
      do i=1,nxj
      flag(i)=.not.ocean(i)
      enddo
!
! skin temperature is part of the update when snow is absent
!
!
      do i=1,nxj
      if(flag(i).and. .not.flagsnw(i))then
        yy(i)=twork(i)+(rnet(i)-gfx(i)-evap(i))            &
              /(rsmall(i)*rch(i))
        zz(i)=1.+dft0(i)/(-.5*zsoil(i,1)*rch(i)*rsmall(i))
        xx(i)=dft0(i)*(stsoil(i,1)-yy(i))/                 &
              (0.5*zsoil(i,1)*zz(i))
      endif
!
      if(flag(i).and. flagsnw(i))then
        yy(i)=stsoil(i,1)
        zz(i)=1.
        xx(i)=dfsnow*(stsoil(i,1)-tsurf(i))           &
              /(-factsnw(i)*max(snr(i)/1000.,0.001))
      endif
      enddo
!
      do i=1,nxj
      if(flag(i))then
       smcz(i)=max(smc(i,1),smc(i,2))
       dtdz1(i)=(stsoil(i,1)-stsoil(i,2))/(-.5*zsoil(i,2))
       if(land(i))then
       dft1(i)=xtsoil(smcz(i),istyp(i),tsat,dfkt)
       hcpct(i)=smc(i,1)*ch2o+(1.-smc(i,1))*csoil
       else
        dft1(i)=dft0(i)
        hcpct(i)=cice
       endif
       dft2(i)=dft1(i)
       ddz(i)=1./(-.5*zsoil(i,2))
!
       ai(i,1)=0.
       bi(i,1)=dft1(i)*ddz(i)/(-zsoil(i,1)*hcpct(i))
       ci(i,1)=-bi(i,1)
       bi(i,1)=bi(i,1)+dft0(i)/(.5*zsoil(i,1)**2*hcpct(i)*zz(i))
       rhstc(i,1)=(dft1(i)*dtdz1(i)-xx(i))/(zsoil(i,1)*hcpct(i))
      endif
      enddo
!
      do k=2,km
        do i=1,nxj
        if(land(i))then
        hcpct(i)=smc(i,k)*ch2o+(1.-smc(i,k))*csoil
        elseif(ice(i))then
        hcpct(i)=cice
        endif
        enddo
!
      if(k.lt.km)then
        do i=1,nxj
        if(flag(i))then
         dtdz2(i)=(stsoil(i,k)-stsoil(i,k+1))/(.5*(zsoil(i,k-1)-zsoil(i,k+1)))
         smcz(i)=max(smc(i,k),smc(i,k+1))
         if(land(i))then
       dft2(i)=xtsoil(smcz(i),istyp(i),tsat,dfkt)
         endif
         ddz2(i)=2./(zsoil(i,k-1)-zsoil(i,k+1))
         ci(i,k)=-dft2(i)*ddz2(i)/((zsoil(i,k-1)-zsoil(i,k))*hcpct(i))
        endif
        enddo
      else                !k=km
        do i=1,nxj
        if(land(i))then
         dtdz2(i)=(stsoil(i,k)-tg3(i))/(.5*(zsoil(i,k-1)+zsoil(i,k))-zbot)
       dft2(i)=xtsoil(smc(i,k),istyp(i),tsat,dfkt)
         ci(i,k)=0.
        endif
        if(ice(i))then
         dtdz2(i)=(stsoil(i,k)-tgice)/(.5*zsoil(i,k-1)-.5*zsoil(i,k))
         dft2(i)=dft1(i)
         ci(i,k)=0.
        endif
        enddo
      endif
!
      do i=1,nxj
        if(flag(i))then
        rhstc(i,k)=(dft2(i)*dtdz2(i)-dft1(i)*dtdz1(i))   &
                  /((zsoil(i,k)-zsoil(i,k-1))*hcpct(i))
        ai(i,k)=-dft1(i)*ddz(i)                          &
                /((zsoil(i,k-1)-zsoil(i,k))*hcpct(i))
        bi(i,k)=-(ai(i,k)+ci(i,k))
!
         if(k .lt. km )then
         dft1(i)=dft2(i)
         dtdz1(i)=dtdz2(i)
         ddz(i)=ddz2(i)
         endif
        endif
      enddo
      enddo
!
! solve tri-diagonal matrix
!
      do k=1,km
         do i=1,nxj
         if(flag(i))then
         rhstc(i,k)=rhstc(i,k)*dt
         ai(i,k)=ai(i,k)*dt
         bi(i,k)=1.+bi(i,k)*dt
         ci(i,k)=ci(i,k)*dt
         endif
         enddo
      enddo
! forward elimination
      do i=1,nxj
         if(flag(i))then
         ci(i,1)=-ci(i,1)/bi(i,1)
         rhstc(i,1)=rhstc(i,1)/bi(i,1)
         endif
      enddo
      do k=2,km
         do i=1,nxj
         if(flag(i))then
         cc=1./(bi(i,k)+ai(i,k)*ci(i,k-1))
         ci(i,k)=-ci(i,k)*cc
         rhstc(i,k)=(rhstc(i,k)-ai(i,k)*rhstc(i,k-1))*cc
         endif
         enddo
      enddo
! backward substitution
       do i=1,nxj
          if(flag(i))then
          ci(i,km)=rhstc(i,km)
          endif
       enddo
       do k=km-1,1,-1
        do i=1,nxj
           if(flag(i))then
           ci(i,k)=ci(i,k)*ci(i,k+1)+rhstc(i,k)
           endif
        enddo
       enddo
! update soil and ice temp
!
        do k=1,km
           do i=1,nxj
           if(flag(i))then
           stsoil(i,k)=stsoil(i,k)+ci(i,k)
           endif
           enddo
        enddo
! update surface temp for snow free
!
!
        do i=1,nxj
         if(.not.ocean(i) .and. .not. flagsnw(i))then
         tsurf(i)=(yy(i)+(zz(i)-1.)*stsoil(i,1))/zz(i)
         endif
!
         if(ice(i).and. .not. flagsnw(i))then
         tsurf(i)=min(tsurf(i),t0c)
         endif
        enddo
!
!
        do k=1,km
         do i=1,nxj
         if(ice(i))then
         stc(i,k)=min(stsoil(i,k),t0c)
         endif
         enddo
        enddo
!
! time filter for soil and skin temperature
!
        do i=1,nxj
        if(.not.ocean(i))then
         tg(i)=ctfil1*tsurf(i)+ctfil2*tg(i)
        endif
        enddo
!
        do k=1,km
         do i=1,nxj
          if(.not.ocean(i))then
          stc(i,k)=ctfil1*stsoil(i,k)+ctfil2*stc(i,k)
          endif
         enddo
        enddo
!
! recalculate gflux
       do i=1,nxj
       if(.not. ocean(i))then
        if(flagsnw(i))then
         gfx(i)=-dfsnow*(tg(i)-stc(i,1))/            &
                (factsnw(i)*max(snr(i)/1000.,0.001))
        else
         gfx(i)=-dft0(i)*(tg(i)-stc(i,1))/(-0.5*zsoil(i,1))
        endif
       gfx(i)=max(gfx(i),-200.)
       gfx(i)=min(gfx(i),200.)
       endif
       enddo

! after tg updated ,recalculate tstar which will be
! transfered to mixpbl as lower boundary
!
! for ncep's land version
!     do i=1,nxj
!     if(.not.ocean(i))then
!     bl   =1.
!     al   =1.0-bl
!     xx1   =dt*sfcw(i)*phih(i) / (2.0*hgt(i,kk))
!     wt(i)=( tt(i,kk) * (1.-al*xx1)+ xx1*tg(i)/pk2(i,kk))/(1.+bl*xx1)
!     endif
!     enddo
      do 400 i = 1, nxj
!cc   tstar(i)=phih(i)*(sfcw(i)*(tt(i,kk)*pk2(i,kk)-tg(i)))/ustar(i)
      tstar(i)=phih(i)*(sfcw(i)*(wt(i)*pk2(i,kk)-tg(i)))/ustar(i)
!
      hflux(i) = -rosfc(i)*cp*ustar(i)*tstar(i)
      qflux(i) = -rosfc(i)*hltm*ustar(i)*qstar(i)
  400 continue
!
!     update ut, vt, tt, and qt, if itypbl=1
!
      if ( itypbl .eq. 1 )  then
         itypbl=0
!soil
!        do 420 i = 1, nxj
!        ut(i,kk) = wu(i)
!        vt(i,kk) = wv(i)
!        tt(i,kk) = wt(i)
!        qt(i,kk) = wq(i)
! 420    continue
!soil
      endif
!
!           compute surface layer e and eps
!
      hpbl  = 450.0
      ghpbl = g*hpbl
      do 500 i = 1, nxj
      e(i,kk)= 3.75*ustar(i)*ustar(i)
      eps(i,kk)= ustar(i)**3 / (0.4*hgt(i,kk))
  500 continue
      do 520 i = 1, nxj
      zl(i) =g04*hgt(i,kk)*tstar(i)/(wt(i)*ustar(i)*ustar(i))
      zlx      = min(max(zl(i),-2.0),2.)
!
         if( zlx .lt.0.)then
      wstar = ghpbl*(-ustar(i)*tstar(i)) / wt(i)
      wstar = max(wstar,1.0e-20)
      wstar = wstar**0.3333
!f    wstar = exp(0.3333*log(wstar))
!
      e(i,kk)  = e(i,kk) + ustar(i)*ustar(i)*(-zlx)**0.6666
!f    e(i,kk)  = e(i,kk) + ustar(i)*ustar(i)*exp(0.6666*log(-zlx))
!                + 0.2*wstar*wstar
         endif
  520 continue
!
      do 540 i = 1, nxj
      e(i,kk)  = min (250.0, max (1.0e-4, e(i,kk)))
      eps(i,kk)= min (1.0,   max (1.0e-7, eps(i,kk)))
  540 continue

      return
      end

!---- function area -----
!
      real function xtsoil(theta,ktype,tsat,dfkt)
!
! xtsoil : heat diffusivity (m2/s)
!
      parameter (ntype=9,ngrid=22)
      dimension dfkt(ngrid,ntype),tsat(ntype)
      w = (theta / tsat(ktype)) * 20. + 1.
      kw = w
      kw = min(kw,21)
      kw = max(kw,1)
      xtsoil = dfkt(kw,ktype)+(w-kw)*(dfkt(kw+1,ktype)-dfkt(kw,ktype))
      return
      end

      function funckt(theta,ktype,tsat,xktk)
!
! funckt : hydraulic conductivity (m/s)
!
      parameter (ntype=9,ngrid=22)
      dimension xktk(ngrid,ntype),tsat(ntype)
      w = (theta / tsat(ktype)) * 20. + 1.
      kw = w
      kw = min(kw,21)
      kw = max(kw,1)
      funckt = xktk(kw,ktype)+(w-kw)*(xktk(kw+1,ktype)-xktk(kw,ktype))
      return
      end

      function funcdf(theta,ktype,tsat,dfk)
!
! funcdf : hydraulic diffusivity (m2/s)
!
      parameter (ntype=9,ngrid=22)
      dimension dfk(ngrid,ntype),tsat(ntype)
      w = (theta / tsat(ktype)) * 20. + 1.
      kw = w
      kw = min(kw,21)
      kw = max(kw,1)
      funcdf = dfk(kw,ktype)+(w-kw)*(dfk(kw+1,ktype)-dfk(kw,ktype))
      return
      end

      function twlt(ktype)
!
! twlt : soil wilting point(should depend on soil type)
!
      twlt = .1
      return
      end
