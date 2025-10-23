      subroutine pbl_noah (nxj,nx,lev,ktpbl,dt,g,r,cp,xkapa,hltm,ptop   &
                        , tice,hice,tg,z0,land,topo,phi,phii,pss,u,v,t  &
                        , q,ut,vt,tt,qt,pk,pk2,ustar,tstar,qstar,e,eps  &
                        , hflux,qflux,itstp,gwclim,tgclim,ocean,ice     &
                        , sheleg,totalp,ss,rs,alb,imx,xkmx,idg,xkmd     &
                        , itype,t2,q2,rh2,rh10,u10,v10,fm,fh,fm10,fh2   &
                        , srflag,rld,stbo                               &
                        , km,smc,stc,canopy,runoff,sigmaf,istyp,ivegtyp &
                        , ncld,dsigma,islopetyp,slc,sncover,snwdph      &
                        , shdmax,shdmin,snoalb,albedo2                  &
                        , sld,zice,cice,xtice,hpbl,asl,atl,xmu,gfx      &
                        , kpbl,nmpbl,nmmiph,jj,isot,ivegsrc,sfemis_g    &
#ifdef TIMCOMCPL
                        , dudt,dvdt,dtdt,dqdt,ustress,vstress,ssu, ssv)
#else
                        , dudt,dvdt,dtdt,dqdt)
#endif
!
!#######################################################################
!                     subroutine description
!
! 1. function description
!
!    this version is based on the Troen and Mahrt (1986) to deal with
!    the vertical diffusion with account of non-local flux.
!    the code is origially from ncepgfs
!
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
!   itstp: indicated which time step it is
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
!     created       november 1991
!     modified      july     1992,  april   1993
!     modified      jan      2005  by f. wang
!     modified      may      2009  by f. wang
!
!            ed., american meteorology socity, 67-98.
!  mixpbl:
!       troen, i., and l. mahrt, 1986: a simple model of the atmospheric
!          boundary layer: sensitivity to surface evaporation. bound.
!          layer meteor., 37, 129-148.
!
!  land process:
!       Mahrt, L., and H. -L. Pan,1984: Atwo-layer model of soil hydrology.
!             Boundary Layer Meteorol., 29, 1-20.
!       Pan, H.-L., and L.Mahrt, 1987: Interaction between soil hydrology
!             and boundary-layer development. Boundary Layer Meteorol.,
!             38,185-202.
!  modify to f90 by C-H Lee and sort by River Chen in 2015
!
!#####################################################################
!
      use mpe
      use rank
      use index
      use radn,   only:ntcw,ntiw,ntinc,ntrnc,ntoz,ntrw,ntsw,ntgl
      use const,  only:RTYPE,mom4ice 
!ch   use paramt

!
      implicit none
!
! input & output variable
!
      integer  nxj,nx,lev,ktpbl,jj,itype,idg,ncld,km,nmpbl,nmmiph,nc
      real     dt,g,r,cp,xkapa,hltm,ptop,stbo,hice,tice

      integer  imx(2),itstp

      real     tg(nx),z0(nx),                                              &
               u(nx,lev),v(nx,lev),t(nx,lev),                              &
               ustar(nx),                                                  &
               tstar(nx),qstar(nx),e(nx,lev),eps(nx,lev),hflux(nx),        &
               qflux(nx),pkd(nx),pk2d(nx),gwclim(nx),                      &
               tgclim(nx),snr(nx),totalp(nx),                              &
               ss(nx),rs(nx),alb(nx),xkmx(2),xkmd(lev),                    &
#ifdef TIMCOMCPL
               t2(nx),u10(nx),v10(nx),ustress(nx),vstress(nx),ssu(nx),ssv(nx)
#else
               t2(nx),u10(nx),v10(nx)
#endif
      real(kind=RTYPE) qt(nx,lev*ncld),q(nx,lev*ncld),phi(nx,lev),         &
                       topo(nx),pss(nx),ut(nx,lev),vt(nx,lev),tt(nx,lev),  &
                       pk(nx,lev),pk2(nx,lev)
!soil
      real     smc(nx,km),stc(nx,km),canopy(nx),sigmaf(nx),                &
               rld(nx),runoff(nx)

      integer  istyp(nx),ivegtyp(nx)
!
      logical  land(nx),ocean(nx),ice(nx)
!
!  local work arrays
!
!ch   real     hgt(im,lm),xkm(im,lm),xkh(im,lm),ts(im),                    &
!ch            qsfc(im),zl(im),czh(im),ps(im),sfcw(im),                    &
!ch            dhgt(im,lm),ro2(im,lm),dhgtz(im,lm)
      real     hgt(nx),xkm(nx,lev),xkh(nx,lev),                            &
               qsfc(nx),zl(nx),czh(nx),ps(nx),sfcw(nx),                    &
               dhgt,ro2(nx),dhgtz(nx,lev)
!soil
      real     rhscnpy(nx),rhsmc(nx,km),aim(nx,km),bim(nx,km),             &
               cim(nx,km),drain(nx),snomt(nx),zsoil(nx,km),                &
               runof(nx),gfx(nx),t850(nx),gesh(nx)
!
!ch   logical snow(im)
      logical snow(nx)
!
! ncep_pbl
!ch   integer   kpbl(im)
!ch   real      rb(im),fm(im),fh(im),hpbl(im),                          &
!ch             heat(im),evap(im),stress(im),                           &
!ch             prsl(im,lm),prslk(im,lm),phil(im,lm),del(im,lm),        &
!ch             prsi(im,lm+1),phi2(im,lm+1),phii(im,lm+1),              &
!ch             dsigma(lm,2),rcl(im),                                   &
      integer   kpbl(nx)
      real      rb(nx),fm(nx),fh(nx),hpbl(nx),                          &
                heat(nx),evap(nx),stress(nx),                           &
                prsl(nx,lev),prslk(nx,lev),phil(nx,lev),del(nx,lev),        &
!byl                prsi(nx,lev+1),phi2(nx,lev+1),phii(nx,lev+1),              &
                prsi(nx,lev+1),phii(nx,lev+1),                          &
                rcl(nx),                                                &
                u1(nx,lev),v1(nx,lev),t1(nx,lev)
      real(kind=RTYPE) dsigma(lev,2)
      real, dimension(:,:,:), allocatable :: q1
!
      real      pk2x(nx,lev),pkx(nx,lev)
!
! noah  ------------
! input and output
      integer islopetyp(nx),islmsk(nx),io,jo
      real      sld(nx),                                                &
                slc(nx,km),                                             &
                sheleg(nx),snwdph(nx),sncover(nx),                      &
                zice(nx),cice(nx),xtice(nx),                            &
                rh2(nx),rh10(nx),                                       &
                shdmin(nx),shdmax(nx),snoalb(nx),albedo2(nx),sfemis_g(nx)
! local
      real      z0rl(nx),prslki(nx),cd(nx),cdq(nx),                     &
                tsurf(nx),psi(nx),prsl1(nx),                            &
                fm10(nx),fh2(nx),fh10(nx),                              &
                qsurf(nx),evapc(nx),cmm(nx),chh(nx),ep1d(nx),           &
!                radsl(nx) ,tprcp(nx),                                   &
                tprcp(nx),                                              &
!                phy_f2d(nx),q2(nx) 
                ddvel(nx),q2(nx) 
!
      real      srflag(nx),sfemis(nx)
!orig logical   flag_guess(nx),flag_iter(nx),mom4ice(nx)
      logical   flag_guess(nx),flag_iter(nx)
!
!2018 new
      logical   redrag
      integer   ivegsrc,isot
!----------------------------------------------------
! for new pbl: asl,atl,xmu
      real      asl(nx,lev),atl(nx,lev),swh(nx,lev),hlw(nx,lev),xmu(nx)

!----------------------------------------------------
! for fractional step:
      real      dudt(nx,lev),dvdt(nx,lev),dtdt(nx,lev),dqdt(nx,lev) 

      integer  lsm,i,k,iter,kc,ntrac
      real     ppd,ppp,ttt,ppu,dth,p850,ddd,cc,qqq
!
! noah mode
      lsm   =  1
      redrag =.false.
!------------------------------------------------------------
!     io=67
!     jo=251
!     io=589
!     jo=148
!     io=705
!     jo=236
      io=480
      jo=240
!
      ntrac=ncld
!      if ( nmmiph .eq. 8 ) ntrac=ncld-4
!      if ( nmmiph .eq.18 ) ntrac=ncld-4

      allocate(q1(nx,lev,ntrac))
!
! --- ensure ktpbl selection is greater than 2
!
      if ( ktpbl .lt. 2 )  ktpbl = 2
!
! --- compute surface pres and surface air temp at current time level
!
      do k=1,lev
      pkd(:) =pk(:,k)
      pk2d(:)=pk2(:,k)
      call vlog(pk2x(1,k),pk2d,nxj)
      call vlog(pkx(1,k), pkd, nxj)
      do i=1,nxj
        pk2x(i,k)=pk2x(i,k)*(cp/r)
        pkx(i,k) = pkx(i,k)*(cp/r)
      enddo
      call vexp(pk2x(1,k),pk2x(1,k),nxj)
      call vexp(pkx(1,k), pkx(1,k), nxj)
      enddo
!
      do 50 i = 1, nxj
      ps(i)  = pss(i) + ptop
      hgt(i) = (phi(i,lev) - topo(i) ) / g
  50  continue
!
      do 120 i = 1, nxj
      qqq = max(qt(i,lev), 1.0e-8)
      ttt = tt(i,lev)*(1.+0.608*qqq)
      ppp = pkx(i,lev)*1000.
      ro2(i) = 100.0*ppp/( r*ttt)
  120 continue
!
!      do 130 k = 1, lev-1
!      do 130 i = 1, nxj
!      ppu = pkx(i,k)*1000.
!      ppd = pkx(i,k+1)*1000.
!      ro2(i,k) = ( ppd - ppu ) * 100. / ( phi(i,k) - phi(i,k+1) )
!  130 continue
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
! NOAH land soil model ....................................
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!..........................................
! define some variable needed in noah sub
!..........................................
       do i=1,nxj
! build slimsk table for input
        if(ocean(i))islmsk(i)=0
        if(land(i))islmsk(i) =1
        if(ice(i))islmsk(i)  =2
! build roughness from ustar over ocean for input
        if(itstp.eq.1 .and. ocean(i))z0(i)=ustar(i)*ustar(i)*0.014/g ! get z0 from ustar
!
        prslki(i)            =pk2(i,lev)/pk(i,lev)   !(ps/p1)**r/cp
!     sfemis   - real, sfc lw emissivity (fractional)
        sfemis(i)          =sfemis_g(i)
!        radsl(i)           =-ss(i)-rld(i)  ! snet + rld  upward
!
        tprcp(i)            =totalp(i)/1000.  ! dth precip (m)
!
!        radsl(i)            =-ss(i)-rld(i)  ! snet + rld  upward
!
!       shdmin(i)           =0.01
!       shdmax(i)           =0.99
!       snoalb(i)           =0.6
!        phy_f2d(i)          =0.
        ddvel(i)            =0.
        rcl(i)              =1.
       enddo
!..............................................................
!  some variable pass from input and may be changed in this sub
!  so these variables should be put back to input
!.............................................................
!      do k=1,lev
!       do i=1,nx
!       t(i,k)              =t(i,k)*pk(i,k) ! transfer t to temp
!       enddo
!      enddo

!     io=67
!     if(jj.eq.jo)then
!     print*,'before noah t(in  tmpt)-----------------'
!     do k=1,lev
!     print*,k,tt(io,k),pk(io,k)
!     enddo
!     endif
!
       do i=1,nxj
        z0rl(i)             =z0(i)*100.         ! transport from m to cm
!       sheleg(i)           =snr(i)    !snow equiv water depth (mm)
!       snwdph(i)           =snr(i)*10. !snow real depth (mm)
!        psi(i)              =ps(i)*0.1  ! surface pressure (mb to cb)
!        prsl(i,1)           =pkx(i,lev)*100.     !cb
!noah new
        psi(i)              =ps(i)*100.  ! surface pressure (mb to pa)
        prsl1(i)            =pkx(i,lev)*1000.*100.     !pa
       enddo
!...........................
! initialize  variable
!..........................
        do i=1,nxj
        tsurf(i)         = tg(i)
        flag_guess(i)    =.false.
        flag_iter(i)     =.true.
        drain(i)         =  0.
        ep1d(i)          =  0.
        runof(i)         =  0.
        hflux(i)         =  0.
        qflux(i)         =  0.
        enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
       do i=1,nxj
       srflag(i) = 0.
       if(t850(i).le.273.16)then
         srflag(i)=1.
       endif
       enddo
    
!.............................................
! loop needed to remove unstable in calm situation (sfcw < 2m/s)
!............................................
       do iter =1, 2
!
!surface exchange coefficient
!
      call sfc_diff(nxj,nx,psi,ut(1,lev),vt(1,lev),tt(1,lev),qt(1,lev), &
!                    tg,z0rl,cd,cdq,rb,                                  &
!                    rcl,prsl(1,1),work1,slmsk,                          &
!                    ustar,sfcw,phy_f2d,fm10,fh2,                        &
!                    sigmaf,ivegtyp,shdmax,                              &
!                    tsurf,flag_iter)
                    hgt,snwdph,tg,z0rl,cd,cdq,rb,                       &
                    prsl1,prslki,islmsk,                                &
                    stress,fm,fh,                                       &
                    ustar,sfcw,ddvel,fm10,fh2,fh10,                   &
                    sigmaf,ivegtyp,shdmax,ivegsrc,                      &
#ifdef TIMCOMCPL
                    tsurf,flag_iter,redrag,ustress,vstress,ssu,ssv)
#else
                    tsurf,flag_iter,redrag)
#endif
!
!     print*,'pblnoah,diff'
         do i=1,nxj
           if((iter.eq.1) .and. (sfcw(i).lt.2.))                        &
              flag_guess(i) = .true.
         enddo
!
! surface flux over ocean
!
!      call sfc_ocean(nxj,nx,km,psi,ut(1,lev),vt(1,lev),tt(1,lev),qt(1,lev), &
!                     tg,qsurf,evapc,gfx,cd,cdq,rcl,prsl(1,1),work1,     &
!                     slmsk,qflux,hflux,ep1d,phy_f2d,                    &
!                     flag_iter)
      call sfc_ocean(nxj,nx,psi,ut(1,lev),vt(1,lev),tt(1,lev),qt(1,lev), &
                     tg,cd,cdq,prsl1,prslki,islmsk,ddvel,flag_iter,      &
#ifdef TIMCOMCPL
                     qsurf,gfx,qflux,hflux,ep1d,ssu,ssv)
#else
                     qsurf,gfx,qflux,hflux,ep1d)
#endif
!
!     print*,'pblnoah,ocean'
!
!     if(jj.eq.jo)then
!     print*,'before sfcdrv-----------------'
!     print*,'  k     u     v      t        q'
!     do k=1,lev
!     print*,k,ut(io,k),vt(io,k),tt(io,k),qt(io,k)
!     enddo
!     print*,'sheleg sncover snwdph=',sheleg(io),sncover(io),snwdph(io)
!     print*,'istyp,ivegtyp,sigmaf=',istyp(io),ivegtyp(io),sigmaf(io)
!     print*,'-----------------------------'
!     print*,'ss,rld,sld,radsl=',ss(io),rld(io),sld(io),radsl(io)
!     print*,'prsl,work1,slmsk=',prsl(io,1),work1(io),slmsk(io)
!     print*,'tg,tsurf=',tg(io),tsurf(io)
!     print*,'smc     slc      stc='
!     do k=km,1,-1
!     print*,k,smc(io,k),slc(io,k),stc(io,k)
!     enddo
!     endif
!
!      call sfc_drv(nxj,nx,km,psi,ut(1,lev),vt(1,lev),tt(1,lev),qt(1,lev),&
!                     sheleg,sncover,snwdph,tg,qsurf,tprcp,SRFLAG,       &
!                     smc,stc,slc,evapc,istyp,sigmaf,                    &
!                     ivegtyp,canopy,rld,sld,                            &
!                     radsl,dth,tgclim,gfx,cd,cdq,                       &
!                     rcl,prsl(1,1),work1,slmsk,                         &
!                     drain,qflux,hflux,ep1d,phy_f2d,                    &
!                     runof,SLOPETYP,SHDMIN,SHDMAX,SNOALB,alb,           &
!                     tsurf, flag_iter, flag_guess,albedo2,jj,io,jo)
      call sfc_drv(nxj,nx,km,psi,ut(1,lev),vt(1,lev),tt(1,lev),qt(1,lev),   &
                     istyp,ivegtyp,sigmaf,sfemis,rld,sld,ss,dth,tgclim,     &
                     cd,cdq,prsl1,prslki,hgt,islmsk,ddvel,islopetyp,        &
                     shdmin,shdmax,snoalb,alb,flag_iter,flag_guess,         &
                     isot,ivegsrc,                                          &
                     sheleg,snwdph,tg,tprcp,srflag,                         &
                     smc,stc,slc,canopy,tsurf,z0rl,                         &
                     sncover,qsurf,gfx,                                     &
                     drain,qflux,hflux,ep1d,runof,                          &
                     albedo2,jj,io,jo)
!
!       call sfc_sice(nxj,nx,km,psi,ut(1,lev),vt(1,lev),tt(1,lev),qt(1,lev),   &
!                      zice,cice,xtice,sld,        &    ! FOR SEA-ICE - XW Nov04
!                      sheleg,snwdph,tg,qsurf,tprcp,SRFLAG,stc,evapc,    &
!                      rld,radsl,SNOMT,dth,gfx,cd,cdq,                   &
!                      rcl,prsl(1,1),work1,slmsk,                        &
!                      qflux,hflux,ep1d,phy_f2d,flag_iter,               &
!                      mom4ice,lsm)
       call sfc_sice(nxj,nx,km,psi,ut(1,lev),vt(1,lev),tt(1,lev),qt(1,lev),  &
                     dth,sfemis,rld,ss,sld,srflag,cd,cdq,prsl1,prslki,       &
                     islmsk,ddvel,flag_iter,mom4ice,lsm,                     &
                     zice,cice,xtice,sheleg,tg,tprcp,stc,ep1d,               &
                     snwdph,qsurf,snomt,gfx,qflux,hflux) 
!
        do i=1, nxj
          flag_iter(i)  = .False.
          flag_guess(i) = .False.
          if((islmsk(i) .eq. 1) .and. (iter .eq. 1)) then
            if(sfcw(i).lt.2.) flag_iter(i) = .true.
          endif
        enddo
      enddo                                 !!!!! <---- Clu_q2m_iter
!
!** update near surface fields
!
!      call sfc_diag(nxj,nx,km,psi,ut(1,lev),vt(1,lev),tt(1,lev),qt(1,lev), &
!                    tg,qsurf,u10,v10,t2,q2,rcl,work1,slmsk,             &
!                    qflux,fm,fh,fm10,fh2,fh10,rh2,rh10)
      call sfc_diag(nxj,nx,psi,ut(1,lev),vt(1,lev),tt(1,lev),qt(1,lev),  &
                    tg,qsurf,u10,v10,t2,q2,prslki,                       &
                    qflux,fm,fh,fm10,fh2,fh10,rh2,rh10)    

!

!     if(jj.eq.jo)then
!     print*,'after noah -----------------'
!     print*,'k       tt       t2         q2        rh2-----------------'
!     print*,k,tt(io,lev),t2(io),q2(io),rh2(io)
!     print*,'qflux,hflux=',qflux(io),hflux(io)
!     endif
!      do k=1,lev
!      do i=1,nx
!       t(i,k)          =t(i,k)/pk(i,k)   ! transfer t back to poten temp
!      enddo
!      enddo
!

       do i=1,nxj
        z0(i)             =z0rl(i)/100.         ! transport from cm to m
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
!
!      if ( nmmiph.eq.8 .or. nmmiph.eq.18 ) then ! Thompson
!        do k=1,lev
!          kc=lev-k+1
!          do i=1,nxj
!            q1(i,kc,1) = qt(i,             k)
!            q1(i,kc,2) = qt(i,lev*(ntcw-1)+k)
!            q1(i,kc,3) = qt(i,lev*(ntiw-1)+k)
!            q1(i,kc,4) = qt(i,lev*(ntinc-1)+k)
!            q1(i,kc,5) = qt(i,lev*(ntoz-1)+k)
!          enddo
!        enddo
!      else
        do nc=1,ntrac
          do k=1,lev
            kc=lev-k+1
            do i=1,nxj
              q1(i,kc,nc) = qt(i,lev*(nc-1)+k)
            enddo
          enddo
        enddo
!      endif
!
       do k=1,lev
          kc=lev-k+1
          do i=1,nxj
            t1(i,kc) = tt(i,k)
            u1(i,kc) = ut(i,k)
            v1(i,kc) = vt(i,k)
!
           prslk(i,kc)=pk(i,k)
           prsl(i,kc)=pkx(i,k)*100000.     !Pa
           prsi(i,kc)=pk2x(i,k)*100000.    !Pa
           del(i,kc)=(pss(i)*dsigma(k,1)+dsigma(k,2))*100. !Pa
           phil(i,kc)=phi(i,k)-topo(i)
!byl           phii(i,kc+1)=phi2(i,k)
           enddo
       enddo
!
       do i=1,nxj
       prsi(i,lev+1)=ptop*100.     !Pa
!byl       phii(i,1)=phi2(i,lev+1)
       enddo
!
       do i=1,nxj
       heat(i)=hflux(i)
       evap(i)=qflux(i)
!      stress(i)=stress(i)
       rcl(i)=1.
       enddo
!
!     if(jj.eq.jo)then
!     print*,'before mixpbl uvtq-----------------'
!     do k=lev,1,-1
!     print*,k,u1(io,k),v1(io,k),t1(io,k),q1(io,k,1),q1(io,k,2)
!     enddo
!     print*,'before  phii,p2,phil,p-----------------'
!     do k=lev,1,-1
!     print*,k,phii(io,k),prsi(io,k),phil(io,k),prsl(io,k)
!     enddo
!     print*,'before  del prslk -----------------'
!     do k=lev,1,-1
!     print*,k,del(io,k),prslk(io,k)
!     enddo
!     print*,'sfc tg,sfcw,heat,evap,stress--='
!     print*,tg(io),sfcw(io),heat(io),evap(io),stress(io)
!     print*,'sfc rb,fm,fh--='
!     print*,rb(io),fm(io),fh(io)
!     endif
!
       if(nmpbl.eq.1)then
         if(myrank.eq.0)                                               &
          print *,' warning !!!, nmpbl can not be 1, reassign nmpbl=2'
         nmpbl=2
       endif
       if(nmpbl .eq. 2)then
       call mixpbl_n( nx,nxj,lev,ntrac,u1,v1,t1,q1,pk2(1,lev),rb,fm,fh  &
                 , tg,heat,evap,stress,sfcw,kpbl                        &
                 , prsi,del,prsl,prslk,phii,phil,rcl,dt                 &
                 , hpbl,ice,g,cp,hltm,r,jj)
       endif
!
!---
! new_sas pbl
       if(nmpbl .eq. 3)then
       do k=1,lev
          kc=lev-k+1
       do i=1,nxj
         swh(i,kc)=asl(i,k)/86400.
         hlw(i,kc)=atl(i,k)/86400.
       enddo
       enddo
       call mixpbl_new( nx,nxj,lev,ntrac,ntcw,u1,v1,t1,q1,swh,hlw,xmu   &
                 , pk2(1,lev),rb,fm,fh                                  &
                 , tg,heat,evap,stress,sfcw,kpbl                        &
                 , prsi,del,prsl,prslk,phii,phil,rcl,dt                 &
                 , hpbl,ice,g,cp,hltm,r,jj)
!
       endif

!---
! YSU pbl scheme
       if(nmpbl .eq. 5)then
       do k=1,lev
          kc=lev-k+1
       do i=1,nxj
         swh(i,kc)=asl(i,k)/86400.
         hlw(i,kc)=atl(i,k)/86400.
       enddo
       enddo
        call     ysu2d(u1,v1,t1,q1,prsl,prsi,prslk,pk2(1,lev),        &
                   nx,nxj,lev,ntrac,del,cp,g,r,hltm,phii,phil,psi,    &
                   z0,stress,hpbl,kpbl,fm,fh,islmsk,heat,evap,sfcw,rb, &
                   dt,rcl,u10,v10,swh,hlw,xmu,jj)
!
       endif

!
! NCEP GFS moninedmf
       if(nmpbl .eq. 4)then
       do k=1,lev
          kc=lev-k+1
       do i=1,nxj
         swh(i,kc)=asl(i,k)/86400.
         hlw(i,kc)=atl(i,k)/86400.
       enddo
       enddo
       call moninedmf( nx,nxj,lev,ntrac,ntcw,u1,v1,t1,q1,swh,hlw,xmu   &
                 , pk2(1,lev),rb,z0rl,u10,v10,fm,fh                    &
                 , tg,heat,evap,stress,sfcw,kpbl                   &
                 , prsi,del,prsl,prslk,phii,phil,dt,hpbl)
!
       endif
!
!      if ( nmmiph.eq.8 .or. nmmiph.eq.18) then ! Thompson
!        do k=1,lev
!          kc=lev-k+1
!          do i=1,nxj
!            qt(i,lev*(ntcw-1)+k) = q1(i,kc,2)
!            qt(i,lev*(ntiw-1)+k) = q1(i,kc,3)
!            qt(i,lev*(ntinc-1)+k)= q1(i,kc,4)
!            qt(i,lev*(ntoz-1)+k) = q1(i,kc,5)
!          enddo
!        enddo
!      else
        do nc=2,ntrac
          do k=1,lev
            kc=lev-k+1
            do i=1,nxj
              qt(i,lev*(nc-1)+k) = q1(i,kc,nc)
            enddo
          enddo
        enddo
!      endif
!
       do k=1,lev
!jh       do k=ktpbl,lev
          kc=lev-k+1
          do i=1,nxj
! time split
!            tt(i,k) = t1(i,kc)
!            ut(i,k) = u1(i,kc)
!            vt(i,k) = v1(i,kc)
! no time split
            dtdt(i,kc) = (t1(i,kc)-tt(i,k))/dt
            dudt(i,kc) = (u1(i,kc)-ut(i,k))/dt
            dvdt(i,kc) = (v1(i,kc)-vt(i,k))/dt
            dqdt(i,kc) = (q1(i,kc,1)-qt(i,k))/dt
          enddo
       enddo
!
       deallocate(q1)

!     call maxp ( nx,xkm(1,lev)  ,xkmx(1),imx(1) )
!     call maxp ( nx,xkm(1,lev-1),xkmx(2),imx(2) )
!     if ( idg .ne. 0 )  then
!       do 260 k = 1, lev
!       xkmd(k) = xkm(idg,k)
! 260   continue
!     endif
!
       do i=1,nxj
       hflux(i)=hflux(i)*ro2(i)*cp          ! transfer to W/m2
       qflux(i)=qflux(i)*ro2(i)*hltm        ! transfer to W/m2
       tstar(i)=-heat(i)/ustar(i)
       qstar(i)=-evap(i)/ustar(i)
       runoff(i)=(drain(i)+runof(i))*dth + runoff(i) !wei add at 20231012
       enddo
!     if(jj.eq.jo)then
!     print*,'after mixpbl -----------------'
!     do k=lev,1,-1
!     print*,k,u1(io,k),v1(io,k),t1(io,k),q1(io,k,1),q1(io,k,2)
!     enddo
!     print*,'sfc kpbl,hpbl=',kpbl(io),hpbl(io)
!     print*,'hflux,qflux(w/m2)=',hflux(io),qflux(io)
!     endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                      c
!     change tt and tg back to temp                                    c
!                                                                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     do 300 k = 1, lev
!     do 300 i = 1, nx
!     tt(i,k) = tt(i,k)*pk(i,k)
! 300 continue
!soil do 320 i = 1, nx
!     tg(i) = tg(i)*pk2(i,lev)
! 320 continue
!soil
      return
      end
