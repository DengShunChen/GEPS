      subroutine pbltke ( nxj,nx,lev                                    &
                        , ktpbl,dt,g,r,cp,xkapa,hltm,ptop,tice          &
                        , hice,tg,z0,land,topo,phi,pss,u,v,t,q,ut       &
                        , vt,tt,qt,pk,pk2,ustar,tstar,qstar,e,eps,hflux &
                        , qflux,gwclim,tgclim,ocean,ice,snr             &
                        , totalp,ss,rs,alb,imx,xkmx,idg,xkmd,itype      &
                        , t2,rh2,u10,v10                                &
                        , rld,stbo                                      &
                        , km,smc,stc,canopy,runoff,sigmaf,istyp,ivegtyp &
                        , wlt,ref,tsat,dfkt,xktk,dfk)
!
!#######################################################################
!                     subroutine description
!
! 1. function description
!
!       this subroutine is the driver of multilayer pbl calculation
!       to compute surface fluxes and vertical mixing due to eddies.
!       the surface layer calculation follows businger (1973)
!       and the eddy mixing coef calculation follows
!       detering and elting (1984)
!       (use ustar, tstar and qstar as lower b.c. for u,v,t,q update)
!
! 2. parameter specification
!
!     nx : dimension of horizontal-direction
!    lev : dimension of z-direction
!  ktpbl : starting level to consider vertical mixing by pbl
!     dt : time step for ut, vt, tt, qt updates                (s)
!      g : gravity                                             (m/s2)
!      r : air gas constant                                    (j/kg)
!     cp : air specific heat constant                          (j/kg)
!  xkapa : r/cp
!   hltm : water vapor latent heat constant                    (j/kg)
!   ptop : model top pressure                                  (mb)
!   tice : temperature of melting point                        (k)
!   hice : latent heat of snow melting                         (j/kg)
!     tg : ground real temp. (sea or land)  (nx)               (k)
!     z0 : surface roughness length    (nx)                    (m)
!   gwet : ground wetness              (nx)
!   land : logic for bare soil land    (nx)  (=true for bare land)
!   topo : terrain geopotential   (nx)                         (m2/s2)
!    phi : odd level geopotential (nx)                         (m2/s2)
!    pss : (terrain pres -ptop)   (nx)     at current time lvl (mb)
!      u : x-direction velocity   (nx,lev) at current time lvl (m/s)
!      v : y-direction velocity   (nx,lev) at current time lvl (m/s)
!      t : potent temp on odd lvl (nx,lev) at current time lvl (k)
!      q : moisture    on odd lvl (nx,lev) at current time lvl (kg/kg)
!     ut : x-direction velocity   (nx,lev) at past/future  lvl (m/s)
!     vt : y-direction velocity   (nx,lev) at past/future  lvl (m/s)
!     tt : real temp.  on odd lvl (nx,lev) at past/future  lvl (k)
!     qt : moisture    on odd lvl (nx,lev) at past/future  lvl (kg/kg)
!     pk : p**capa on odd levels  (exner func)
!    pk2 : p**capa on even levels (exner func)
!   ustar: friction velocity         (nx)                      (m/s)
!   tstar: potential temp scale      (nx)                      (k)
!   qstar: moisture scale            (nx)                      (kg/kg)
!     e  : turbulence kinetic energy (nx,lev) at even lvls     (m2/s2)
!    eps : tke dissipation           (nx,lev)                  (m2/s2)
!   hflux: upward surface heat flux         (nx)               (w/m2)
!   qflux: upward surface moisture flux     (nx)               (w/m2)
!    gwr : soil water amount           (nx)                    (mm)
!  gwclim: soil water climate          (nx)                    (mm)
!  tgclim: climate ground temperature  (nx)                    (k)
!   ocean: logic for open water      (nx)  (=true for open water)
!     ice: logic for ice cover       (nx)  (=true for ice covered)
!     snr: snow depth in water state (nx)                      (mm)
!  totalp: total precipitation       (nx)                      (mm)
!      ss: net short wave radiat. flux down to ground (nx)     (wat/m2)
!      rs: net long  wave radiat. flux up from ground (nx)     (wat/m2)
!     alb: diffusion albedo           (nx)
!     imx: address of the two points with max xkm at the two lowest lvl
!    xkmx: max mixing coef xkm at two lowest levels            (m2/s)
!     idg: i-index address of selected point to print profile
!    xkmd: xkm profile at the selected idg point               (m2/s)
!   itype: surface layer update type, =0 use surface stress as b.c.
!                                     =1 directly update lowest levels
!      t2: 2m air temp derived by similarity theory (nx)     (K)
!     rh2: 2m air rh (nx)
!     u10: 10m air x-direction wind derived by simil. theory(nx)  (m/s)
!     v10: 10m air y-direction wind derived by simil. theory(nx)  (m/s)
!
!soil io
!
!     smc: volumetric soil moisture content (nx,km_soil)
!     stc: soil temperature                 (nx,km_soil)
!  canopy: canopy moisture content(<0.5mm)          (nx)
!  sigmaf: green vegetation fraction        (nx)
!   istyp: soil type(1-9)                   (nx)
! ivegtyp: vegetation type(1-13)            (nx)
!  runoff: accumulate run off water(runof+drain) (nx)             (mm)
!     rld: long wave radiat. flux down to ground (nx)         (wat/m2)
!   wlt  : wiltling point of 9 soil types  (9)
!   ref  : field capacity of 9 soil types  (9)
! tsat   : satuation point of 9 soil types (9)
! dfkt   : soil thermal diffusivity     (22,9)     m2/s
! xktk   : soil hydraulic conductivity  (22,9)     m/s
!  dfk   : soil hydraulic diffusivity   (22,9)     m2/s
!
! 4. local variable
!
!    phi : geopential on odd lvl  (nx,lev) at current time lvl (m2/s2)
!    xkm : eddy mixing coef. for u,v (nx,lev)                  (m2/s)
!    xkh : eddy mixing coef. for t,q (nx,lev)                  (m2/s)
!      ts: surface air temperature    (nx)                     (k)
!    qsfc: surface moisture (specific humidity)      (nx)      (kg/kg)
!     zl : height(lev) over monin-obukhov scale height.
!          a stability indicator.                    (nx)
!     czh: ground thermal capacity                   (nx)      (j/kg.m2)
!    snow: logic for snow cover                      (nx)
!
!soil
!
!  rhscnpy: righ hand side terms of canopy prediction(nx)
!    rhsmc: righ hand side terms of smc prediction(nx,km)
!      aim: working array smc prediction           (nx,km)
!      bim: working array smc prediction           (nx,km)
!      cim: working array smc prediction           (nx,km)
!   drain: drainage water from bottom soil layer   (nx)        (mm/s)
!     snomt: snow melt rate                        (nx)        (mm/s)
!     zsoil: soil depth of soil layer(-0.1,-2.)    (nx,km)        (m)
!     runof: run off water from surface and soil   (nx)        (mm/s)
!      gfx: ground heat flx                        (nx)        (wat/m2)
!     t850: 850mb air temp                         (nx)         (K)
!     gesh: control number,control precipitation   (nx)
!           conver to snow(0) or not(1)
!
!
! 5. calling modules
!
!     main program
!
! 6. usage
!
!     call pbltke ( nx,lev,ktpbl,dt,g,r,cp,xkapa,hltm,ptop,tice
!    1            , hice,tg,z0,gwet,land,topo,phi,pss,u,v,t,q,ut
!    2            , vt,tt,qt,pk,pk2,ustar,tstar,qstar,e,eps,hflux
!    3            , qflux,gwr,gwclim,tgclim,ocean,ice,snr
!    4            , totalp,ss,rs,alb,imx,xkmx,idg,xkmd,itype
!    5            , rld,stbo
!    6            , km,smc,stc,canopy,runoff,sigmaf,istyp,ivegtyp
!    7            , wlt,ref,tsat,dfkt,xktk,dfk)
!
! 7. modules called
!
!    qsatq, lansfc, landsmc, mixpbl
!
! 8. limitation
!
!    see landsfc, landsmc, mixpbl
!
! 9. date
!
!     created       november 1991
!     modified      july     1992,  april   1993
!     modified      jan      2005  by f. wang
!
! 10. author
!
!      c. liou / f. wang
!
! 11. reference
!
!  surface layer:
!       businger j.a.,1973: turbulent transfer in the atmospheric
!            surface layer.  workshop in meteorology, d. a. haugen,
!            ed., american meteorology socity, 67-98.
!  mixpbl:
!       detering, h. w., and d. etling, 1984: application of the
!            e-eps turbulence model to the atmospheric boundary
!            layer. boundary-layer meteorol., 33, 113-133.
!
!  land process:
!       Mahrt, L., and H. -L. Pan,1984: Atwo-layer model of soil hydrology.
!             Boundary Layer Meteorol., 29, 1-20.
!       Pan, H.-L., and L.Mahrt, 1987: Interaction between soil hydrology
!             and boundary-layer development. Boundary Layer Meteorol.,
!             38,185-202.
!#####################################################################

      use paramt
      use const, only: RTYPE

      implicit none

!
! input & output variable
!

      integer  nxj,nx,lev,ktpbl,idg,itype,km
      real     dt,g,r,cp,xkapa,hltm,ptop,tice,hice,stbo

      integer  imx(2)

      real     tg(nx),z0(nx),                                        &
               u(nx,lev),v(nx,lev),t(nx,lev),                        &
               ustar(nx),                                            &
               tstar(nx),qstar(nx),e(nx,lev),eps(nx,lev),hflux(nx),  &
               qflux(nx),pkd(nx),pk2d(nx),gwclim(nx),                &
               tgclim(nx),snr(nx),totalp(nx),                        &
               ss(nx),rs(nx),alb(nx),xkmx(2),xkmd(lev),              &
               t2(nx),rh2(nx),u10(nx),v10(nx)
      real(kind=RTYPE) phi(nx,lev),qt(nx,lev),q(nx,lev),topo(nx),    &
                       pss(nx),tt(nx,lev),ut(nx,lev),vt(nx,lev),     &
                       pk(nx,lev),pk2(nx,lev)
!soil
      real     smc(nx,km),stc(nx,km),canopy(nx),sigmaf(nx),          &
               rld(nx),runoff(nx)

      integer  istyp(nx),ivegtyp(nx)
!
      logical  land(nx),ocean(nx),ice(nx)
!
!  local work arrays
!
!     include '../include/paramt.h'
      real     hgt(im,lm),xkm(im,lm),xkh(im,lm),ts(im),              &
               qsfc(im),zl(im),czh(im),ps(im),sfcw(im),              &
               dhgt(im,lm),ro2(im,lm),dhgtz(im,lm)
!soil
      real     rhscnpy(nx),rhsmc(nx,km),aim(nx,km),bim(nx,km),       &
               cim(nx,km),drain(nx),snomt(nx),zsoil(nx,km),          &
               runof(nx),gfx(nx),t850(nx),gesh(nx)
!
      integer, parameter :: ntype=9, ngrid=22

      real     wlt(ntype),ref(ntype),tsat(ntype),dfkt(ngrid,ntype),  &
               xktk(ngrid,ntype),dfk(ngrid,ntype)
!soil
!
      real     pk2x(nx,lev),pkx(nx,lev)

      logical snow(im)

      integer   i,k,kc
      real      ppp,ttt,ppu,dth,p850,ddd,cc,ppd
!
! --- ensure ktpbl selection is greater than 2
!
      if ( ktpbl .lt. 2 )  ktpbl = 2
!
! --- compute surface pres and surface air temp at current time level
!
      do 50 i = 1, nxj
      ps(i) = pss(i) + ptop
      ts(i) = t(i,lev)*pk2(i,lev)
  50  continue

      do k = 1, lev
      pkd(:) =pk(:,k)
      pk2d(:)=pk2(:,k)
      call vlog(pk2x(1,k),pk2d,nxj)
      call vlog(pkx(1,k), pkd, nxj)
      do i = 1, nxj
        pk2x(i,k)=pk2x(i,k)*(cp/r)
        pkx(i,k) = pkx(i,k)*(cp/r)
      enddo
      call vexp(pk2x(1,k),pk2x(1,k),nxj)
      call vexp(pkx(1,k), pkx(1,k), nxj)
      enddo
!
      do 100 k = 1, lev
      do 100 i = 1, nxj
      hgt(i,k) = (phi(i,k) - topo(i) ) / g
  100 continue
!
      do 105 i = 1, nxj
      ppd = pk2x(i,1) * 1000.
      dhgt(i,1) = ( ppd - ptop ) * 100. / g
      ppp = pkx(i,1) * 1000.
      ttt = t(i,1)*pk(i,1)*(1.+0.608*q(i,1))
      dhgtz(i,1)= dhgt(i,1) * r * ttt / (100.*ppp)
  105 continue
!
      do 110 k=2, lev
      do 110 i=1, nxj
      ppu=pk2x(i,k-1) * 1000.
      ppd=pk2x(i,k) * 1000.
      dhgt(i,k) = ( ppd - ppu ) * 100. / g
      ppp = pkx(i,k) * 1000.
      ttt = t(i,k)*pk(i,k)*(1.+0.608*q(i,k))
      dhgtz(i,k)= dhgt(i,k) * r * ttt / (100.*ppp)
  110 continue
!
      do 120 i = 1, nxj
      ro2(i,lev) = 100.0*ps(i)/( r*ts(i))
  120 continue
!
      do 130 k = 1, lev-1
      do 130 i = 1, nxj
      ppu = pkx(i,k)*1000.
      ppd = pkx(i,k+1)*1000.
      ro2(i,k) = ( ppd - ppu ) * 100. / ( phi(i,k) - phi(i,k+1) )
  130 continue
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                      c
!     determine ground condition and predict ghround temp and wetness  c
!                                                                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
      dth = dt
!
!                                                                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
! soil ....................................
!        to get environment temp at 850mb
!
      p850=850.
      do i=1,nxj
      t850(i)=tt(i,lev)
      enddo
!
      do k=lev,2,-1
      do i=1,nxj
      ppd = pkx(i,k)*1000.
      ppu = pkx(i,k-1)*1000.
      ddd = (p850-ppd)*(p850-ppu)
        if( ddd .lt. 0. )then
         cc=(log(p850/ppd))/(log(ppu/ppd))
         t850(i)=(tt(i,k-1)-tt(i,k))*cc+tt(i,k)
        endif
      enddo
      enddo
!
!soil......................................
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                      c
!     change tt and tg from temp to potential temp                     c
!                                                                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 240 k = 1, lev
      do 240 i = 1, nxj
      tt(i,k) = tt(i,k)/pk(i,k)
  240 continue
!soil do 250 i = 1, nxj
!     tg(i) = tg(i)/pk2(i,lev)
! 250 continue
!soil
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                      c
!     surface layer parameterization:                                  c
!     compute surface fluxes and save the surafce fluxex as            c
!     lower boundary cond for implicit time integration in mixpbl      c
!                                                                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! soil
!
      call landsfc( nxj,nx,lev,dt,g,r,cp,hltm,tg,z0,ocean,ps,ts,hgt,u,v &
                 , t,q,ut,vt,tt,qt,sfcw,ustar,tstar,qstar,e,eps,hflux   &
                 , qflux,zl,itype                                       &
                 , alb,ss,rld,stbo,pk,pk2,gfx                           &
                 , km,snr,smc,stc,canopy,sigmaf,istyp,ice,land          &
                 , tgclim,rhscnpy,rhsmc,aim,bim,cim,drain,snomt         &
                 , zsoil,dth,ivegtyp,t2,rh2,u10,v10                     &
                 , wlt,ref,tsat,dfkt,xktk,dfk)
!
      do i=1,nxj
       if(t850(i).gt.273.16)then
         gesh(i)=1.
       else
         gesh(i)=0.
       endif
       if(.not.ocean(i))then
        snr(i)=snr(i)+totalp(i)*(1.-gesh(i))
       endif
       totalp(i)=gesh(i)*totalp(i)
      enddo
!
      call landsmc(nxj,nx,km,rhscnpy,rhsmc,aim      &
                  ,bim,cim,smc,land                 &
                  ,canopy,totalp,runof,snomt        &
                  ,zsoil,istyp,sigmaf,dt,dth,tsat)
      do i=1,nxj
      runoff(i)=(drain(i)+runof(i))*dth+runoff(i)
      enddo
!soil
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                      c
!     mixing layer parameterization:                                   c
!     compute vertical eddy mixing coeff. and update all level         c
!     u, v, t and q according to the mixing and surface drag           c
!                                                                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      call mixpbl( nxj,nx,lev,ktpbl,dt,g,hgt,u,v,t,q,ut,vt,tt,qt,e,eps  &
                 , xkm,xkh,zl,sfcw,ustar,tstar,qstar,itype              &
                 , dhgt,ro2,dhgtz)
!
      call maxp ( nxj,xkm(1,lev)  ,xkmx(1),imx(1) )
      call maxp ( nxj,xkm(1,lev-1),xkmx(2),imx(2) )
      if ( idg .ne. 0 )  then
        do 260 k = 1, lev
        xkmd(k) = xkm(idg,k)
  260   continue
      endif
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                      c
!     change tt and tg back to temp                                    c
!                                                                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      do 300 k = 1, lev
      do 300 i = 1, nxj
      tt(i,k) = tt(i,k)*pk(i,k)
  300 continue
!soil do 320 i = 1, nxj
!     tg(i) = tg(i)*pk2(i,lev)
! 320 continue
!soil
      return
      end
