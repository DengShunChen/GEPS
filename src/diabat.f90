      subroutine diabat ( fwd,docup,dodry,dolsp,dopbl,dorad,doshl,dograv       &
                    , nx,my,my_max,lev,ncld,nmcup,nmpbl,nmland,nmshl,cgw       &
                    , idg,jdg,ldiag,dt,tau,hours,julian                        &
                    , frad,ozon,njump,itypbl,ktcup,ktpbl,ktshl,grav            &
                    , rgas,cp,stbo,s0,evaprh,hltm,ptop,sigma,dsigma,il,ib,cof  &
                    , xlat,xlon,sgeo,z0,alb,land,ocean,ice,snr                 &
                    , tg,tgclim,curate,plcl,cumtop,totalp,raintot,raincu       &
                    , rainlp,raincu6,rainlp6,hflux,qflux,ustar,tstar,qstar,e   &
                    , eps,o3l,dtrad,ss,rs,plt,pk,pk2,ps,up,vp,ttp,qp           &
                    , pst,ut,vt,tt,qt,gwclim,tice,hice,qgini                   &
                    , thdai,tengi,acld,std,qbrwtot,asol,olr,drag               &
                    , ugws,vgws,sdpbl,t2,rh2,u10,v10,gfx                       &
                    , rld,km_soil,smc,stc,canopy,runoff                        &
                    , sigmaf,istyp,ivegtyp,wlt,ref,tsat,dfkt,xktk,dfk          &
                    , ftp,fqp,fpsp,ftp1,fqp1,fpsp1,sd                          &
                    , shdmax,shdmin,snoalb                                     &
                    , slopetyp,sld,slc,zice,cice,xtice,sncover,sndepth         &
                    , ctot,chig,cmid,clow,hpbl,asl,atl,cosz                    &
!--------------------------------------------------------------------------------
                    , fusl,fdsl,fuir,fdir                                      &
                    , fuslr,fdslr,fuirr,fdirr                                  &
                    , asl_clr,atl_clr,clds                                     &
                    , ss_clr,rs_clr,asol_clr,olr_clr,sld_clr,rld_clr           &
                    , alvsf,alvwf,alnsf,alnwf,facsf,facwf                      &
                    , idtg,doo3l,nfxr)
!--------------------------------------------------------------------------------
!#######################################################################
!
!     driver for all physical parametrization of diabatic processes
!
!
!     parameters
!
!     fwd    : logical variable for forward time step or not
!     docup  : logical variable for including cumulus parameterization
!     dodry  : logical variable for including dry convection
!     dolsp  : logical variable for including large scale precipitation
!     dopbl  : logical variable for including pbl parameterization
!     dorad  : logical variable for including radiation parameterization
!     doshl  : logical variable for including shallow convection calcul
!     dograv : logical variable for including gravity wave drag
!     nx     : x-dimension of model grid
!     my     : y-dimension of model grid
!     lev     : total vertical levels of model
!     idg    : i-index of select point for diagnostic prints
!     jdg    : j-index of selcet point for diagnostic prints
!     ldiag  : index to control the level of runtime diagnostic
!            = 1, basic (normal) diagnostics
!            = 2, advanced disgnostics (cupcwb and radtn diag)
!            = 3, extensive diagnostics
!     dt     : time step (for forward time step, in seconds)
!     tau    : current forecast time (in hours)
!     hours  : current zulu hours
!     julian : current julian day
!     frad   : frequency to call radiation package (in hours)
!     ozon   : logical variable for including ozone in radiative calcul
!     njump  : grid interval for long wave radiative calculation
!     itypbl : index for surface layer update, =0 stress b,c., =1 direct
!     ktcup  : highest model level allowed for cup top
!     ktpbl  : highest model level allowed for vertical mixing calcul
!     ktshl  : highest model level aloowed for shallow convection cal
!     grav   : gravity constant (=9.806)
!     rgas   : dry gas constant (=287.)
!     cp     : specific heat coef. for dry air (=1004.)
!     stbo   : stefan-boltzman constant (=5.669e-8)
!     s0     : solar constant (=1368.3)
!     evaprh : relative humidity where lsp condensation begins
!     hltm   : latent heat constant for water vapor (=2.52e+6)
!     ptop   : model top pressure for radiation calculation (=0.01)
!     sigma  : sigma layer
!     dsigma : sigma layer thickness
!     il     : output point address for 1-d spline interp in lw rad
!     ib     : input  point address for 1-d spline interp in lw rad
!     cof    : coeff. array for 1-d spline interpolation in lw rad
!     xlat   : latitude at each model grid in y-direction
!     xlon   : longtitude at each model grid in x-direction
!     sgeo   : surface terrain geopotential  (m2/s2)
!     z0     : surface rougness (m)
!     alb    : surface ground albedo
!     land   : logical variable for bare soil land grid points
!     ocean  : logical variable for open water grid points
!     ice    : logical variables for ice covered grid points
!     snr    : accumulated surface snow depth (m)
!     tg     : ground surface temperature (k)
!     tgclim : climate value for deep soil temperature (k)
!     curate : cumulus precipitation rate for cloud diagnost (mm/s)
!     plcl   : sigma level values at lcl level (cumulus cloud base)
!     cumtop : sigma level values at cumulus cloud top
!     totalp : total precipitation rate for snow accumulation (mm/s)
!     raincu : accumulated cumulus rain in each output period (mm)
!     rainlp : accumulated large-scale rain in each output period (mm)
!     hflux  : upward surface sensible heat flux (w/m2)
!     qflux  : upward surface latent heat flux (w/m2)
!     ustar  : surface friction velocity  (m/s)
!     tstar  : surface friction temperature (k, <0 for upward flux)
!     qstar  : surface friction mixing ratio ( <0 for upward flux)
!     e      : turbulence k.e. at present/future time step (m2/s2)
!     eps    : turbulence k.e. dissipation rate at present/future time
!     o3l    : ozone concentration for radiative transfer calculation
!     dtrad  : temperature change rate due to radiative transfer (k/day)
!     ss     : net downward solar radiation at ground surface (w/m2)
!     rs     : net upward long wave radiation at ground surface (w/m2)
!     ps     : (terrain pressure-ptop) at present time step (mb)
!     plt    : odd-lvl pressure at present time step (mb)
!     pk     : odd-levl p**kcapa at present time step
!     pk2    : even level p**kapa at present time step
!     up     : u-component at present time step (m/s)
!     vp     : v-component at present time step (m/s)
!     ttp    : virtual potential temperature at present time step (k)
!     qp     : specific humidity at present time step (kg/kg)
!     pst    : (terrain pressure-ptop) at future time step (mb)
!     ut     : u-component at past/future time step (m/s)
!     vt     : v-component at past/future time step (m/s)
!     tt     : virtual potential temp at past/future time step (k)
!     qt     : specific humidity at time/future time step (kg/kg)
!     gwclim : soil water climate value (= gwet-clim * 20.0)
!     tice   : critical temp to separate rain and snow (=273.15)
!     hice   : specific heat constant to melt snow     (=3.336e+5)
!     qgini  : initial global moisture amount
!     acld   :
!     std    : standard deviation of terrain
!     asol   : asorbed solar radiative flux by earth
!     olr    : outgoing longwave radiation
!     drag   : gravity wave drag force
!     ugws   : zonal component of surface gravity wave stress 
!     vgws   : meridional component of surface gravity wave stress 
!     sdpbl  : vertical velocity within pbl    
!     t2     : 2-meter temperature
!     rh2    : 2-meter rh
!     u10/v10: 10-meter wind speed
!     gfx    : ground heat flx                        (nx)        (wat/m2)
!soil
!     rld    : long wave radiat. flux down to ground (nx)    (wat/m2)
!     smc    : volumetric soil moisture content (nx,km_soil)
!     stc    : soil temperature                 (nx,km_soil)
!  canopy    : canopy moisture content(<0.5mm)          (nx)
!  runoff    : accumulate run off water(runof+drain) (nx)    (mm)
!  sigmaf    : green vegetation fraction        (nx)
!   istyp    : soil type(1-9)                   (nx)
! ivegtyp    : vegetation type(1-13)            (nx)
!   wlt      : wiltling point of 9 soil types  (9)
!   ref      : field capacity of 9 soil types  (9)
!   tsat     : satuation point of 9 soil types (9)
!   dfkt     : soil thermal diffusivity    (22,9)     m2/s
!   xktk     : soil hydraulic conductivity (22,9)  m/s
!micro
!     ftp, fqp, fpsp : for temperature, specific humidity and surface
!              pressure for past time step for micro cloud scheme
!
! modify to f90 by C-H Lee and sort by River Chen in 2015
!
!#######################################################################
!
      use mpe
      use rank
      use index
!-----------------------------------------------------------------------
      use radn
      use physpara
!-----------------------------------------------------------------------
      implicit  none
!-----------------------------------------------------------------------
      integer nfxr, ntrac, kk
      real slag,sdec,cdec,solcon,dtlw,dtsw,solhr

! --- for radupdat
      integer idat(8),jdat(8)

! --- for new albedo
      real alvsf(nxp,my_max),alvwf(nxp,my_max),alnsf(nxp,my_max), &
           alnwf(nxp,my_max),facsf(nxp,my_max),facwf(nxp,my_max)

! --- for random number generator (thread safe mode)
      integer ixseed(nx,my,2)
      integer icsdlw(nx),icsdsw(nx)

! --- for rrtmg : input
!     logical lsswr,lslwr,lssav,lprnt
      logical lsswr,lslwr,lssav
      real xlonr(nx,my_max)

! --- new variables setting :
      integer*8 idtg
      logical doo3l
      integer ipt,jpt
!-----------------------------------------------------------------------
      logical   fwd,docup,dodry,dolsp,dopbl,dorad,doshl,dograv,ozon, &
                land(nxp,my_max),ocean(nxp,my_max),ice(nxp,my_max)

      integer   nx,my,my_max,lev,ncld,nmcup,nmpbl,nmland,nmshl,idg,  &
                jdg,ldiag,julian,njump,itypbl,ktcup,ktpbl,ktshl,     &
                km_soil

      real      tice,hice,qgini,thdai,tengi,qbrwtot,ptop,            &
                hltm,evaprh,s0,stbo,cp,rgas,grav,frad,               &
                hours,tau,dt,cgw

      integer   il(nx,4),ib(nx,4)

      real      sigma(lev+1,2),dsigma(lev,2),                                         &
                cof(nx*3,4),xlat(my),                                                 &
                xlon(nx,my_max),sgeo(nxp,my_max),z0(nxp,my_max),              &
                alb(nxp,my_max),snr(nxp,my_max),tg(nxp,my_max),           &
                tgclim(nxp,my_max),curate(nxp,my_max),plcl(nxp,my_max),   &
                cumtop(nxp,my_max),totalp(nxp,my_max),raincu(nxp,my_max), &
                hflux(nxp,my_max),qflux(nxp,my_max),ustar(nxp,my_max),    &
                tstar(nxp,my_max),qstar(nxp,my_max),                          &
                e(nxp,lev,my_max),eps(nxp,lev,my_max),                        &
                o3l(nxp,lev,my_max),dtrad(nxp,lev,my_max),ss(nxp,my_max), &
                rs(nxp,my_max),plt(nxp,lev,my_max),pk(nxp,lev,my_max),    &
                pk2(nxp,lev,my_max),ps(nxp,my_max),up(nxp,lev,my_max),    &
                vp(nxp,lev,my_max),ttp(nxp,lev,my_max),                       &
                qp(nxp,lev*ncld,my_max),rainlp(nxp,my_max),                   &
                pst(nxp,my_max),ut(nxp,lev,my_max),vt(nxp,lev,my_max),    &
                tt(nxp,lev,my_max),qt(nxp,lev*ncld,my_max),                   &
                gwclim(nxp,my_max),acld(lev,my),std(nxp,my_max),              &
                asol(nxp,my_max),olr(nxp,my_max),drag(nxp,lev,my_max),    &
                ugws(nxp,my_max),vgws(nxp,my_max),sdpbl(nxp,my_max),      &
                raintot(nxp,my_max),t2(nxp,my_max),rh2(nxp,my_max),       &
                u10(nxp,my_max),v10(nxp,my_max),hpbl(nxp,my_max),         &
                raincu6(nxp,my_max),rainlp6(nxp,my_max)
!soil (2005/01/12)
      integer,  parameter :: ntype=9, ngrid=22
      integer   istyp(nxp,my_max),ivegtyp(nxp,my_max)

      real      smc(nxp,km_soil,my_max),stc(nxp,km_soil,my_max),  &
                canopy(nxp,my_max),runoff(nxp,my_max),            &
                sigmaf(nxp,my_max),rld(nxp,my_max),               &
                wlt(ntype),ref(ntype),tsat(ntype),dfkt(ngrid,ntype),      &
                xktk(ngrid,ntype),dfk(ngrid,ntype)

! new soil
! for noah
      integer   slopetyp(nxp,my_max)

      real      slc(nxp,km_soil,my_max),zice(nxp,my_max),         &
                cice(nxp,my_max),xtice(nxp,my_max),               &
                sld(nxp,my_max),sncover(nxp,my_max),              &
                sndepth(nxp,my_max),gfx(nxp,my_max),              &
                shdmax(nxp,my_max),shdmin(nxp,my_max),            &
                snoalb(nxp,my_max),albedo2(nxp,my_max),heat(nxp),evap(nxp),  &
! for new pbl
                asl(nxp,lev,my_max),atl(nxp,lev,my_max),xmu(nxp,my_max) 
!---------------------------------------------------------------------------
! for new rad
!---------------------------------------------------------------------------
      real      fusl(nxp,lev+1,my_max),fdsl(nxp,lev+1,my_max),   &
                fuir(nxp,lev+1,my_max),fdir(nxp,lev+1,my_max),   &
                fuslr(nxp,lev+1,my_max),fdslr(nxp,lev+1,my_max), &
                fuirr(nxp,lev+1,my_max),fdirr(nxp,lev+1,my_max), &
                asl_clr(nxp,lev,my_max),atl_clr(nxp,lev,my_max), &
                clds(nxp,lev,my_max)
      real      rld_clr(nxp,my_max),sld_clr(nxp,my_max)
      real      asol_clr(nxp,my_max),olr_clr(nxp,my_max),ss_clr(nxp,my_max), &
                rs_clr(nxp,my_max),asr_clr(lev,my),alr_clr(lev,my),           &
                ctot(nxp,my_max),chig(nxp,my_max),cmid(nxp,my_max),clow(nxp,my_max)
!---------------------------------------------------------------------------
      real      avgdrag_u(my,lev),avgdrag_v(my,lev),drag_u(lev),drag_v(lev)
      real      fnor
      logical   donor, upnor
!     data      donor/.true./, fnor/1./
      data      donor/.true./, fnor/0.5/
!     data      donor/.false./,fnor/1./
!
!#######################################################################
!
!     local logical variables and work arrays
!
      logical   fluxcl,doozon,uprad
      integer   ijdg(my),ipblmx(2,my),itlsp(lev),nnlsp(lev),ilsp(lev,my),&
                nlsp(lev,my),ncup(my),ndry(my),nshl(my),icupmx(my),      &
                nlcl(lev,my),nnegl(lev,my),nosat(lev,my),nwork(lev,my),  &
                ntcup(lev,my),nflx(lev,my),lvlwx(my),ilx(nx,my),         &
                ibx(nx,my)

      real      cosl(my),sinl(my),cosz(nxp,my_max),                  &
                rcup(nxp,my_max),rlsp(nxp,my_max),               &
                asr(lev,my),alr(lev,my),xsr(lev,my),xlr(lev,my),         &
                aflxd(lev+2,my),aflxu(lev+2,my),                         &
                dtcupx(my),dtcupz(lev,my),dqcupz(lev,my),dtcupd(lev),    &
                dqcupd(lev),dtcupl(lev),dqcupl(lev),xkmx(2,my),xkmd(lev),&
                qbrrow(ncld,my),phi(nxp,lev),theda(nxp,lev),albx(nxp,my_max), &
                cofx(nx*3,my),dphi(nxp,lev)

      real      wkj(4,my),dsigpp(lev),qt_diff(ncld)

      integer   nlcl_tmp(lev) ,nnegl_tmp(lev),nosat_tmp(lev),    &
                nwork_tmp(lev),ntcup_tmp(lev),nflx_tmp(lev)

      real      adtrad(nxp,lev),work_pr1(9),work_pr2(lev,9)
!--------
! for ncld=2
      real,     parameter :: dxmax=-8.8818363, dxmin=-5.2574954, &
                             dxinv=1.0/(dxmax-dxmin)
!     parameter (rhzbot=0.85, rhztop=0.85)

      real      work1(nx),work2(nx),rhc(nxp,lev),rhckt,coefrhc
      real      del(nxp,lev),prsl(nxp,lev),psfc(nxp)
      real      qtc(nxp,lev), qtr(nxp,lev), ttc(nxp,lev)
      real      ftp(nxp,lev,my_max), fqp(nxp,lev,my_max), fpsp(nxp,my_max)
      real      ftp1(nxp,lev,my_max), fqp1(nxp,lev,my_max), fpsp1(nxp,my_max)

! vertcal rhc
!      real      ct,cs,px

      logical lprnt
!--------
!
! for sascnv
      integer   kbot(nxp,my_max),ktop(nxp,my_max),kuo(nxp,my_max)
      real      sl(lev),delcup(lev),slimsk(nxp)
      real      dotc(nxp,lev),phil(nxp,lev),utc(nxp,lev),vtc(nxp,lev)
      real      cldwrk(nxp,my_max),sd(nxp,lev,my_max),xkt2(nx)
! for new shlcon
      real      rcup2(nxp)

!CWB 2007-09-27 for random number seed >>>
      real*8    rtc,rsecond
      integer   isize(2)
! CWB <<<

      real      wk1(nx,my),wk2(nx,my)

      integer   i,     j,      k,      jj,     nxj,    nxmy,   nxlev, levmy, &
                icrad, iter,   icnor,  njump1, njump2, njump3, nny,   jcap,  &
                kc,    nncup,  nnshl,  nndry,  ixkmkk, ixkmk1, jxkmkk,jxkmk1,&
                idcupx,jdcupx, ikutx,  jutx,   ikut,   isamax, idummy,jpstx, &
                ipstx, jpstn,  ipstn,  jhflx,  ihflx,  jqflx,  iqflx, jtgx,  &
                itgx,  kdradx, jdradx, idradx, isign,  jmax,   kmax,  jmin,  &
                kmin,  jjdg

      real      prevap,etop,   radus,  radsq,  d2r,    ptrad,  xkapa, xkapa1,&
                okapa, p0k,    op0k,   ptopk,  dta,    rainfc, abxlat,xx,    &
                yy,    dtau,   aps,    arcup,  arlsp,  evapor, qglb,  thda,  &
                tke,   tpe,    cosq,   dsigp,  cosw,   teng,   qbrw,  rainbl,&
                qdiff, engdiff,thdadif,topsd,  toplu,  sfcsd,  sfclu, xoj,   &
                xkmkk, xkmk1,  dtcupg, utx,    speed,  pstx,   pstn,  hflmx, &
                qflmx, tgx,    dtradg, dragmax,dragmin,rcuprr, rlsprr,ud,    &
                vd,    ttd,    qqd,    dqcu,   deg_ju, arg,    tem

!xb110>
      real      mdlon
      real      mflon(nx,my)
!xb110<
!CWB2015 
      kuo=0

!CWB2016 
      icsdsw=0
      icsdlw=0

!     define local constants
! for vertical rhc
!      ct=0.7
!      cs=0.9
!      px=4.
!
!-------------------------------
! pbl package > tke-e  ; npbl=1
! pbl package > ncep   ; npbl=2
!     npbl = 2
!
! shlcon use Frank version , nmshl=1
! shlcon use ncep new version(2010) companion with new_sas, nmshl=2
!     nmshl = 2
!-------------------------------
      prevap= 0.2
      etop= 1.0
      nxmy = nx * my
      nxlev= nx * lev
      levmy= lev* my
      radus = 6370000.
      radsq = radus**2
      d2r   = 3.141592654 / 180.0
      ptrad = max(0.01,ptop)
      xkapa = rgas/cp
      xkapa1= 1.0 + xkapa
      okapa = 1.0/xkapa
      p0k   = 1000.0**xkapa
      op0k  = 1.0/p0k
      ptopk = ptop**xkapa
!
!     define local control variables
!
!------------------------------------------------------------------------------
      ipt=518
      jpt=490
!------------------------------------------------------------------------------
!      if (myrank .eq. 0)print *,'### diabat start ###'
!
!     define local control variables
!
      fluxcl = .true.
      if (fwd) then
         dta    = dt
         rainfc = 1.0
         doozon = .true.
      else
         dta    = 2.0*dt
         rainfc = 0.5
         doozon = .false.
      endif
!------------------------------------------------------------------------------
!     set hours, iter, icrad, julian, uprad, doozon
!------------------------------------------------------------------------------

      hours = hours + dt/3600.0
      if ( hours .gt. 24.0 )  then
         hours = mod ( hours,24.0 )
         julian= julian + 1
         if ( julian .gt. 365 ) julian = julian - 365
         doozon = .true.
      endif
      icrad = frad*3600.0/dt + 0.0001 ! frad =1.0 set in block.f
      iter  = tau*3600.0/dt + 0.0001
      uprad = .false.
      if ( (mod(iter,icrad).eq.0) .or. (iter.eq.1) )  uprad = .true.
      doozon = doozon .and. dorad
      uprad  = uprad  .and. dorad
!
! for nonorographic gravity wave drag
!
      icnor = fnor*3600.0/dt + 0.0001
      upnor = .false.
      if ( (mod(iter,icnor).eq.0) .or. (iter.eq.1) )  upnor = .true.
      upnor  = upnor  .and. donor
!
      if (.not. dopbl)  then
        do jj = 1, jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
        do i = 1, nxj
         qflux(i,jj) = 0.
         hflux(i,jj) = 0.
        enddo
        enddo
      endif
!
      do jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do i = 1, nxj
       rcup(i,jj)  = 0.0
       rlsp(i,jj)  = 0.0
!      cosz(i,jj)  = 0.0
       xmu(i,jj)  = 0.0
      enddo

      do k = 1, lev
       asr(k,j)   = 0.0
       alr(k,j)   = 0.0
       xsr(k,j)   = 0.0
       xlr(k,j)   = 0.0
       dtcupz(k,j)= 0.0
       dqcupz(k,j)= 0.0
       nlcl(k,j)  = 0
       nnegl(k,j) = 0
       nosat(k,j) = 0
       nwork(k,j) = 0
       ntcup(k,j) = 0
       nflx(k,j)  = 0
       ilsp(k,j)  = 0
       nlsp(k,j)  = 0
      enddo

      do k = 1, lev+2
       aflxd(k,j) = 0.0
       aflxu(k,j) = 0.0
      enddo

      enddo
!
      do jj = 1, jlistnum
       j=jlist1(jj)
       ijdg(j) = 0
!ch    qbrrow(1,j)= 0.0
!ch    qbrrow(2,j)= 0.0
       qbrrow(1:ncld,j)= 0.0
       ncup(j) = 0
       ndry(j) = 0
       nshl(j) = 0
       icupmx(j) = 0
       ipblmx(1,j) = 0
       ipblmx(2,j) = 0
       xkmx(1,j) = 0.0
       xkmx(2,j) = 0.0
       dtcupx(j) = 0.0
      enddo

      do k = 1, lev
       itlsp(k) = 0
       nnlsp(k) = 0
       dtcupd(k) = 0.0
       dqcupd(k) = 0.0
       dtcupl(k) = 0.0
       dqcupl(k) = 0.0
       xkmd(k)   = 0.0
      enddo
!
      if ( ldiag .ge. 2 )  then
       do 110 jj = 1, jlistnum
        j=jlist1(jj)
        if ( j .eq. jdg ) ijdg(j) = idg
  110  continue
      endif
!
!     compute cos and sin of latitude
!
      do 140 j = 1, my
      cosl(j) = cos(xlat(j)*d2r)
      sinl(j) = sin(xlat(j)*d2r)
  140 continue
!
!     if radiation is to be called, compute cos of solar zenith angular
!                njump, il, ib, and cof for different meridional zones
!
      if ( uprad )  then
!
!        call ccoszen ( nx,my,my_max,julian,hours,xlat,xlon,cosz )
! cosz was modified to be an average of the  calling period(1 hour fo
        call coszenpm ( nx,my,my_max,julian,hours,xlat,xlon,frad,cosz )
!
!!ocl scalar,nounroll
         do 160 jj = 1, jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)

         abxlat = abs(xlat(j))
         njump1 = njump + 1
         if ( mod(nxj,njump1) .ne. 0 )  njump1 = njump
         njump2 = njump + 2
         if ( mod(nxj,njump2) .ne. 0 )  njump2 = njump1
         njump3 = njump + 3
         if ( mod(nxj,njump3) .ne. 0 )  njump3 = njump2

         if( lreduce.eq.1 ) then
           njump1 = njump
           njump2 = njump
           njump3 = njump
         endif

         if ( abxlat .le. 20.0 )  then
!ch         lvlwx(j) = nxjp(j)/njump
            lvlwx(j) = nxp/njump
            nny = 1
         else if ( abxlat .le. 60.0 )  then
!ch         lvlwx(j) = nxjp(j)/njump1
            lvlwx(j) = nxp/njump1
            nny = 2
         else if ( abxlat .le. 80.0 )  then
!ch         lvlwx(j) = nxjp(j)/njump2
            lvlwx(j) = nxp/njump2
            nny = 3
         else
!ch         lvlwx(j) = nxjp(j)/njump3
            lvlwx(j) = nxp/njump3
            nny = 4
         endif
!
         if( lreduce.eq.1 ) then
           call splinc (lvlwx(j),nxj,ilx(1,j),ibx(1,j),cofx(1,j))
         else
!          do 150 i  = 1, nx
           do 150 i  = 1, nxp
            ilx(i,j)  = il(i,nny)
            ibx(i,j)  = ib(i,nny)
  150      continue
           do 155 i  = 1, nx*3
            cofx(i,j) = cof(i,nny)
  155      continue
         endif
  160    continue
      endif
!
! calculate xmu,
! the instantaneous zenith angular at the 'hours'
!------------------------------------------------------------------------------
        call ccoszen ( nx,my,my_max,julian,hours,xlat,xlon,xmu )
!------------------------------------------------------------------------------

! transfer xmu to a fraction of average value(cosz) during frad
      do 170 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 170 i = 1, nxj
        if(xmu(i,jj).gt.0.0001.and.cosz(i,jj).gt.0.0001) then
          xmu(i,jj) = xmu(i,jj) / cosz(i,jj)
        else
          xmu(i,jj)   = 0.
        endif
 170  continue

!     initial albx by climate values of gwet and alb, while it may
!     be updated in grdcon according ground conditions
!
      do 180 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do 180 i = 1, nxj
       albx(i,jj)  = alb(i,jj)
       albedo2(i,jj)  = alb(i,jj)
  180 continue
!
!-----------------------------------------------------------------------
      if (uprad .and. irad .eq. 2) then
!     if (myrank .eq. 0) then
!         print *,'### prerrtmg start !'
!         print *,'### for prerrtmg : iter  =',iter
!         print *,'### for prerrtmg : tau = ',tau
!         print *,'### for prerrtmg : nx, my, idtg, dt=',nx,my,idtg,dt
!         print *,'### for prerrtmg : frad, uprad =',frad, uprad
!     endif
      call prerrtmg(nx,my,my_max,idtg,tau,dt,hours,frad,uprad,        &
                    isubc_sw,isubc_lw,d2r,xlon,myrank,me,             &
                    idat,jdat,solhr,dtsw,dtlw,lsswr,lslwr,            &
                    slag,sdec,cdec,solcon,                            &
                    xlonr,ixseed)

!     if (myrank .eq. 0) then
!         print *,'### for prerrtmg : idat=',idat
!         print *,'### for prerrtmg : jdat=',jdat
!         print *,'### for prerrtmg : solhr=',solhr
!         print *,'### for prerrtmg : dtsw =',dtsw,' dtlw=',dtlw
!         print *,'### for prerrtmg : lsswr, lslwr=',lsswr,lslwr
!         print *,'### for prerrtmg : slag=',slag
!         print *,'### for prerrtmg : sdec=',sdec
!         print *,'### for prerrtmg : cdec=',cdec
!         print *,'### for prerrtmg : solcon=',solcon
!     endif
!      if (myrank .eq. 0) print *,'after prerrtmg ok'
      endif ! for uprad .and. irad. eq. 2

!-----------------------------------------------------------------------
      if(ncld.ge.3)then
      ntrac=3
      do k = 1, lev
        kk = (ntrac-1)*lev+k
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do i = 1, nxj
            o3l(i,k,jj) = qt(i,kk,jj)
          enddo
        enddo
      enddo
      endif
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                      c
!     begin big j-loop for diabatic calculation in each latitude ring  c
!                                                                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      do 290 jj =1, jlistnum
      j=jlist1(jj)
      nxj=nxdef_2d(j)
!
!    fill up negative moisture fields from level below
!
!!    call postq (nxjp(j),nxp,lev,ncld,dsigma,pst(1,jj),qt(1,1,jj),qbrrow(1,j),cosl(j))
!
!    compute new time level p**kapa quantites
!
!  plt= new odd level pressure
!  pk= (plt/1000)**capa
!
      call prexp_hybrid_cwb ( nxjp(j),nxp,lev,ptop,sigma,pst(1,jj),  &
                          pk(1,1,jj),pk2(1,1,jj),plt(1,1,jj) )
!
!     hydrostatic equation
!
      do 200 i = 1, nxj
      phi(i,lev) = sgeo(i,jj)+ &
                   cp*tt(i,lev,jj)*(pk2(i,lev,jj)-pk(i,lev,jj))
  200 continue
!
      do 210 k = lev-1, 1, -1
      do 210 i = 1, nxj
      phi(i,k) = phi(i,k+1) +cp*(tt(i,k,jj)*(pk2(i,k,jj)-pk(i,k,jj))   &
                            + tt(i,k+1,jj)*(pk(i,k+1,jj)-pk2(i,k,jj)))
  210 continue
!
!     deweight u,v by cosl/radus, and
!     change t from virtual potential temperature to real temperature
!     change ttp from virtual potential temperature to potential temperature
!
      xx = radus/cosl(j)
      do 230 k = 1, lev
      do 230 i = 1, nxj
      ut(i,k,jj) = ut(i,k,jj)*xx
      vt(i,k,jj) = vt(i,k,jj)*xx
      up(i,k,jj) = up(i,k,jj)*xx
      vp(i,k,jj) = vp(i,k,jj)*xx
      tt(i,k,jj) = tt(i,k,jj)*pk(i,k,jj) / (1.0+0.608*qt(i,k,jj))
      ttp(i,k,jj) = ttp(i,k,jj) / (1.0+0.608*qp(i,k,jj))
  230 continue
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     start physical process calculation (from long to short time scale)
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     update tt by radiation heating/cooling rate: dtrad (k/day)
!
      do 240 k = 1, lev
      do 240 i = 1, nxj
      tt(i,k,jj) = tt(i,k,jj) + dta*dtrad(i,k,jj)/86400.0
  240 continue
!
!
      if ( dopbl .and. nmpbl.eq.1 .and. nmland.eq.1)                          &
         call pbltke ( nxjp(j),nxp,lev,ktpbl,dta,grav,rgas,cp,xkapa,hltm,ptop &
                     , tice,hice,tg(1,jj),z0(1,jj),land(1,jj)                 &
                     , sgeo(1,jj),phi,pst(1,jj),up(1,1,jj),vp(1,1,jj)         &
                     , ttp(1,1,jj),qp(1,1,jj),ut(1,1,jj),vt(1,1,jj)           &
                     , tt(1,1,jj),qt(1,1,jj),pk(1,1,jj),pk2(1,1,jj)           &
                     , ustar(1,jj),tstar(1,jj),qstar(1,jj),e(1,1,jj)          &
                     , eps(1,1,jj),hflux(1,jj),qflux(1,jj),fwd                &
                     , gwclim(1,jj),tgclim(1,jj),ocean(1,jj),ice(1,jj)        &
                     , snr(1,jj),totalp(1,jj),ss(1,jj),rs(1,jj),albx(1,jj)    &
                     , ipblmx(1,j),xkmx(1,j),ijdg(j),xkmd,itypbl              &
                     , t2(1,jj),rh2(1,jj),u10(1,jj),v10(1,jj)                 &
                     , rld(1,jj),stbo                                         &
                     , km_soil,smc(1,1,jj),stc(1,1,jj),canopy(1,jj)           &
                     , runoff(1,jj),sigmaf(1,jj),istyp(1,jj),ivegtyp(1,jj)    &
                     , wlt,ref,tsat,dfkt,xktk,dfk )
      if ( dopbl .and. nmpbl.eq.2 .and. nmland.eq.1)                          &
         call pbltke_n ( nxjp(j),nxp,lev,ktpbl,dta,grav,rgas,cp,xkapa,hltm,ptop &
                     , tice,hice,tg(1,jj),z0(1,jj),land(1,jj)                 &
                     , sgeo(1,jj),phi,pst(1,jj),up(1,1,jj),vp(1,1,jj)  &
                     , ttp(1,1,jj),qp(1,1,jj),ut(1,1,jj),vt(1,1,jj)           &
                     , tt(1,1,jj),qt(1,1,jj),pk(1,1,jj),pk2(1,1,jj)           &
                     , ustar(1,jj),tstar(1,jj),qstar(1,jj),e(1,1,jj)          &
                     , eps(1,1,jj),hflux(1,jj),qflux(1,jj),fwd                &
                     , gwclim(1,jj),tgclim(1,jj),ocean(1,jj),ice(1,jj)        &
                     , snr(1,jj),totalp(1,jj),ss(1,jj),rs(1,jj),albx(1,jj)    &
                     , ipblmx(1,j),xkmx(1,j),ijdg(j),xkmd,itypbl              &
                     , t2(1,jj),rh2(1,jj),u10(1,jj),v10(1,jj)                 &
                     , rld(1,jj),stbo                                         &
                     , km_soil,smc(1,1,jj),stc(1,1,jj),canopy(1,jj)           &
                     , runoff(1,jj),sigmaf(1,jj),istyp(1,jj),ivegtyp(1,jj)    &
                     , wlt,ref,tsat,dfkt,xktk,dfk,ncld,dsigma,j )
      if ( dopbl .and. nmland.eq.2)                                           &
       call pbl_noah ( nxjp(j),nxp,lev,ktpbl,dta,grav,rgas,cp,xkapa,hltm,ptop &
                     , tice,hice,tg(1,jj),z0(1,jj),land(1,jj)                 &
                     , sgeo(1,jj),phi,pst(1,jj),up(1,1,jj),vp(1,1,jj)  &
                     , ttp(1,1,jj),qp(1,1,jj),ut(1,1,jj),vt(1,1,jj)           &
                     , tt(1,1,jj),qt(1,1,jj),pk(1,1,jj),pk2(1,1,jj)           &
                     , ustar(1,jj),tstar(1,jj),qstar(1,jj),e(1,1,jj)          &
                     , eps(1,1,jj),hflux(1,jj),qflux(1,jj),fwd                &
                     , gwclim(1,jj),tgclim(1,jj),ocean(1,jj),ice(1,jj)        &
                     , snr(1,jj),totalp(1,jj),ss(1,jj),rs(1,jj),albx(1,jj)    &
                     , ipblmx(1,j),xkmx(1,j),ijdg(j),xkmd,itypbl              & 
                     , t2(1,jj),rh2(1,jj),u10(1,jj),v10(1,jj)                 & 
                     , rld(1,jj),stbo                                         &
                     , km_soil,smc(1,1,jj),stc(1,1,jj),canopy(1,jj)           &
                     , runoff(1,jj),sigmaf(1,jj),istyp(1,jj),ivegtyp(1,jj)    &
                     , ncld,dsigma                                            &
                     , slopetyp(1,jj)                                         &
                     , slc(1,1,jj),sncover(1,jj),sndepth(1,jj)                &
                     , shdmax(1,jj),shdmin(1,jj),snoalb(1,jj),albedo2(1,jj)   &
                     , sld(1,jj),zice(1,jj),cice(1,jj),xtice(1,jj)            &
                     , hpbl(1,jj),asl(1,1,jj),atl(1,1,jj),xmu(1,jj),gfx(1,jj) &
                     , nmpbl,j )
!     recompute phi by tt after pbl to ensure consistence of phi & phi2
!
      do 250 k = 1, lev
      do 250 i = 1, nxj
      theda(i,k) = tt(i,k,jj)*(1.0+0.608*qt(i,k,jj)) / pk(i,k,jj)
  250 continue
!
      do 252 i = 1 ,nxj
      phi(i,lev) = sgeo(i,jj)+  &
                   cp*theda(i,lev)*(pk2(i,lev,jj)-pk(i,lev,jj))
  252 continue
      do 255 k = lev-1, 1, -1
      do 255 i = 1,nxj
      phi(i,k) = phi(i,k+1) +cp*(theda(i,k)*(pk2(i,k,jj)-pk(i,k,jj))  &
                            +theda(i,k+1)*(pk(i,k+1,jj)-pk2(i,k,jj)))
  255 continue
!
      if(dograv)                                                       &
       call gwdp (j,nxjp(j),nxp,lev,                                   &
                  ut(1,1,jj),vt(1,1,jj),tt(1,1,jj),qt(1,1,jj),         &
                  plt(1,1,jj),pk(1,1,jj),pk2(1,1,jj),phi,std(1,jj),dta, &
                  grav,rgas,cp,drag(1,1,jj),ugws(1,jj),vgws(1,jj),cgw)    
!
      if ( docup .and. nmcup.eq. 1 )                                         &
       call cupcwb (j,nxjp(j),nxp,my,lev,ktcup,dta,grav,rgas,cp,hltm,etop,prevap  &
                 , sgeo(1,jj),pst(1,jj),plt(1,1,jj),pk(1,1,jj),pk2(1,1,jj)   &
                 , tt(1,1,jj),qt(1,1,jj),phi,plcl(1,jj),cumtop(1,jj)         &
                 , rcup(1,jj),ncup(j),ptop,dsigma,icupmx(j),dtcupx(j)        &
                 , dtcupz(1,j),dqcupz(1,j),ijdg(j),dtcupd,dqcupd             &
                 , nlcl(1,j),nnegl(1,j),nosat(1,j),nwork(1,j)                &
                 , ntcup(1,j),nflx(1,j) )
!
!cyea---->
!c 20120926 for Tiedtke cumulus
      if ( docup .and. nmcup .eq. 4 .and. ncld .ge. 2 ) then
        call cumastr_driv(nxjp(j),nxp,lev,dt,grav,rgas,cp,hltm,ptop &
                       , land(1,jj),sgeo(1,jj),phi,up(1,1,jj)  &
                       , vp(1,1,jj),ttp(1,1,jj),qp(1,1,jj)          &
                       , ut(1,1,jj),vt(1,1,jj),tt(1,1,jj)           &
                       , qt(1,1,jj),rcup(1,jj),pk(1,1,jj)           &
                       , pk2(1,1,jj),sd(1,1,jj),qflux(1,jj)         &
                       , kbot(1,jj),ktop(1,jj),fwd,ncld,sigma,plt(1,1,jj),pst(1,jj),j )

!c for rad input of convection cloud information
!c bottom(plcl) layer and top(cumtop) layer in pressure(mb)
        do i=1,nxj
         if(kbot(i,jj).eq.lev-1 .and. ktop(i,jj).eq.lev-1)then
          plcl(i,jj)=0.
          cumtop(i,jj)=0.
         else
          plcl(i,jj)=plt(i,kbot(i,jj),jj)
          cumtop(i,jj)=plt(i,ktop(i,jj),jj)
         endif
        enddo
!c
        do i=1,nxj
          rcup(i,jj) = rcup(i,jj) * 1000.         ! mm/call
        enddo
!c
!cyea    doshl=.false.
!cyea----<
      endif    !(end if nmcup=4)
!
!xb110>
      if ( docup .and. nmcup .eq. 5 .and. ncld .ge. 2 ) then

        do i=1,nxj
          if(land(i,jj))slimsk(i)=1
          if(ocean(i,jj))slimsk(i)=0
          if(ice(i,jj))slimsk(i)=2
        enddo

        mflon(1,j) = 0.
        mdlon = 0.
          do i=2,nxj
            mflon(i,j)=mflon(1,j)+float(i-1)*360./nxj
          enddo
        mdlon = mflon(nxj,j) - mflon(nxj-1,j)

        call cumastr_driv_n                                               &
               (nxjp(j)    ,nxp       ,lev        ,dta        ,grav      ,&
                rgas       ,cp        ,hltm       ,ptop      ,land(1,jj) ,&
                sgeo(1,jj) ,phi       ,up(1,1,jj) ,vp(1,1,jj),ttp(1,1,jj),&
                qp(1,1,jj) ,ut(1,1,jj),vt(1,1,jj) ,tt(1,1,jj),qt(1,1,jj) ,&
                rcup(1,jj) ,pk(1,1,jj),pk2(1,1,jj),sd(1,1,jj),qflux(1,jj),&
                kbot(1,jj) ,ktop(1,jj),fwd        ,ncld      ,sigma      ,&
                plt(1,1,jj),pst(1,jj) ,j          ,slimsk    ,hflux(1,jj),&
                xlat(j)    ,mdlon)

        do i=1,nxj
         if(kbot(i,jj).eq.lev-1 .and. ktop(i,jj).eq.lev-1)then
          plcl(i,jj)=0.
          cumtop(i,jj)=0.
         else
          plcl(i,jj)=plt(i,kbot(i,jj),jj)
          cumtop(i,jj)=plt(i,ktop(i,jj),jj)
         endif
        enddo
!c
!        do i=1,nx
        do i=1,nxj
          rcup(i,jj) = rcup(i,jj) * 1000.         ! mm/call
        enddo
      endif    !(end if nmcup=5)
!xb110<

      if ( nmshl.eq.2 .or. (nmcup .eq. 2 .or. nmcup .eq. 3) ) then
!
! setting for nmcup=2,3 and new shallow convection
!
!orig       call random_number(XKT2)
! CWB 2007-09-27 change random number seed dynamically >>>
        rsecond=rtc()
        isize(1)=rsecond
        isize(2)=(rsecond-isize(1))*100000000.
        call random_seed(put=isize(1:2))
        call random_number(XKT2)
! CWB <<<
        lprnt=.false.
        jcap = 240
        do i=1,nxj
          psfc(i)  = pst(i,jj)*0.1        ! change to cb
          if(land(i,jj))slimsk(i)=1
          if(ocean(i,jj))slimsk(i)=0
          if(ice(i,jj))slimsk(i)=2
        enddo
        do k=1,lev
          kc=lev-k+1
           sl(kc)    = (sigma(k,1)+sigma(k+1,1))*0.5
          do i=1,nxj
           prsl(i,kc) = plt(i,k,jj)*0.1 ! change to cb
           del(i,kc)  = psfc(i)*dsigma(k,1)+dsigma(k,2)*0.1  !unit cb
!ch move out
!ch        sl(kc)    = (sigma(k,1)+sigma(k+1,1))*0.5
           dotc(i,kc)= sd(i,k,jj)*0.1
!201803           phil(i,kc)= phi(i,k)
           phil(i,kc)= phi(i,k)-sgeo(i,jj)
           qtc(i,kc) = qt(i,k,jj)
           qtr(i,kc) = qt(i,lev+k,jj)
           ttc(i,kc) = tt(i,k,jj)
           utc(i,kc) = ut(i,k,jj)
           vtc(i,kc) = vt(i,k,jj)
          enddo
        enddo
!
      if ( docup .and. (nmcup .eq. 2 .or. nmcup .eq. 3) ) then
!
! old version SAS
         if(docup .and. nmcup .eq. 2)                                 &
         call sascnv(nxjp(j),nxp,lev,jcap,dta,del,sl,psfc,prsl,phil,qtr &
          ,qtc,ttc,utc,vtc,dotc,cldwrk(1,jj) &
          ,rcup(1,jj),kbot(1,jj),ktop(1,jj)                     &
          ,kuo(1,jj),slimsk,xkt2,ncld,dtcupx(j),icupmx(j)      &
          ,grav,cp,hltm,rgas,tice)
!
! new version SAS
         if(docup .and. nmcup .eq. 3)                                 &
         call sascnv_n(nxjp(j),nxp,lev,jcap,dta,del,psfc,prsl,phil,qtr &
          ,qtc,ttc,utc,vtc,dotc,cldwrk(1,jj) &
          ,rcup(1,jj),kbot(1,jj),ktop(1,jj)                     &
          ,kuo(1,jj),slimsk,xkt2,ncld                          &
          ,grav,cp,hltm,rgas,tice)

! for rad input of convection cloud information
! bottom(plcl) layer and top(cumtop) layer in pressure(mb)
        do i=1,nxj
         if(kbot(i,jj).eq.lev+1 .and. ktop(i,jj).eq.0)then
          plcl(i,jj)=0.
          cumtop(i,jj)=0.
         else
          plcl(i,jj)=prsl(i,kbot(i,jj))*10.         !to mb
          cumtop(i,jj)=prsl(i,ktop(i,jj))*10.       !to mb
         endif
        enddo
!
        do i=1,nxj
          rcup(i,jj) = rcup(i,jj) * 1000.         ! mm/call
        enddo
!
      endif  !(end of nmcup=2,3)
!
       if( doshl .and. nmshl.eq.2  ) then
! new version shalcon

        do i=1,nxj
           heat(i)=-ustar(i,jj)*tstar(i,jj)
           evap(i)=-ustar(i,jj)*qstar(i,jj)
        enddo
!
        call shalcnv_new(nxjp(j),nxp,lev,jcap,dta,del,prsl,psfc,phil,qtr &
          ,qtc,ttc,utc,vtc                         &
          ,rcup2,kbot(1,jj),ktop(1,jj)                                &
          ,kuo(1,jj),slimsk,dotc,ncld,hpbl(1,jj),heat,evap &
          ,grav,cp,hltm,rgas,tice)
!
        do i=1,nxj
          rcup(i,jj) = rcup(i,jj)+rcup2(i) * 1000.         ! mm/call
        enddo
!
       endif  !(end of nmshl=2)
!
        do k=1,lev
          kc=lev-k+1
          do i=1,nxj
            qt(i,k    ,jj) = qtc(i,kc)
            qt(i,k+lev,jj) = qtr(i,kc)
            tt(i,k    ,jj) = ttc(i,kc)
            ut(i,k    ,jj) = utc(i,kc)
            vt(i,k    ,jj) = vtc(i,kc)
          enddo
        enddo
!
      endif  !(end of nmshl.eq.2 .or. (nmcup .eq. 2 .or. nmcup .eq. 3))

      if ( doshl .and. nmshl .eq.1)                                            &
         call shlcon ( nxjp(j),nxp,lev,ktshl,dta,grav,rgas,cp,xkapa,hltm,ptop  &
                     , dsigma,tg(1,jj),pk(1,1,jj),pst(1,jj),sgeo(1,jj),phi     &
                     , plt(1,1,jj),tt(1,1,jj), qt(1,1,jj),nshl(j)              &
                     , rcup(1,jj),ncld )
!
      if ( dolsp .and. ncld.eq.1 )                                             &
!        call lspmst ( tt(1,1,jj),qt(1,1,jj),plt(1,1,jj),ptop,pst(1,jj)
!                 , dsigma,grav,nxj,nx,lev,evaprh,rlsp(1,jj),cp,hltm
!                 , nlsp(1,j),ilsp(1,j) )
!
         call lsp ( tt(1,1,jj),qt(1,1,jj),plt(1,1,jj),pst(1,jj),dsigma         &
                  , grav,nxjp(j),nxp,lev,evaprh,rlsp(1,jj),cp,hltm,nlsp(1,j)   &
                  , ilsp(1,j) )
!
      if ( dolsp .and. ncld.ge.2 ) then
!
        deg_ju=23.45*sin(d2r*(360./365.)*(julian+284.))
        arg=xlat(j)-deg_ju
        if(arg.gt.90.)then
          arg=89.9999
        else if(arg.lt.-90.)then
          arg=-89.9999
        endif
!
        lprnt=.false.
        do i=1,nxj
!!cjh          work1(i) = (log(cosl(j)/real(nxj)) - dxmin) * dxinv
!          work1(i) = (log(cosl(j)/real(nx)) - dxmin) * dxinv
!          work2(i) = 1.0 - work1(i)
          psfc(i)  = pst(i,jj)*0.1        ! change to cb
        enddo
        do k=1,lev
          kc=lev-k+1
!!!!          rhckt=35.+10.*max(cos(4.*d2r*xlat(j))**3,0.)
!!!!          coefrhc=min(float(kc)/rhckt,1.)**2
          do i=1,nxj
!!           tem   = (rhztop-rhzbot) / (pk(i,k,jj)-pk(i,lev,jj))
!!           tmprhc = rhzbot + tem * (pk(i,k,jj)-pk(i,lev,jj))
!!           rhc(i,kc) = 0.999 * work1(i) + tmprhc * work2(i)
!            rhc(i,kc) = 0.999 * work1(i) + 0.85 * work2(i)
!
!!!            rhc(i,kc)=0.999-0.08*cos(d2r*arg)**2    !a3
            rhc(i,kc)=0.7+0.15*cos(d2r*xlat(j))**2    !a3
!!!!             rhc(i,kc)=(1.-coefrhc)*(0.7+0.15*cos(d2r*xlat(j))**2)  &
!!!!                     +coefrhc*(0.6+0.1*max(cos(4.*d2r*xlat(j))**3,0.))   !PYL vertical profile
!
! vertical profile
!!            tem=ct+(cs-ct)*exp(1.-(pst(i,jj)/plt(i,k,jj)**px))
!!            rhc(i,kc)=rhc(i,kc)*tem
!
!            tem   = (plt(i,k,jj) / plt(i,lev,jj))**0.1
!            if(tem.le.0.85)tem=0.85
!            tem   = (plt(i,k,jj) / plt(i,lev,jj))**0.3
!            if(tem.le.0.65)tem=0.65
!            rhc(i,kc)=rhc(i,kc)*tem
!
!!            if(rhc(i,kc).ge.0.98)rhc(i,kc)=0.98
            prsl(i,kc) = plt(i,k,jj)*0.1 ! change to cb
            del(i,kc) = (dsigma(k,1)*pst(i,jj)+dsigma(k,2))*0.1  ! change to cb
            qtc(i,kc) = qt(i,k,jj)
            qtr(i,kc) = qt(i,lev+k,jj)
            ttc(i,kc) = tt(i,k,jj)
          enddo
        enddo
        call gscond(nxjp(j),nxp,lev,dta,prsl,psfc,  &
                    qtc,qtr,ttc,           &
                    ftp (1,1,jj),fqp (1,1,jj),fpsp (1,jj),&
                    ftp1(1,1,jj),fqp1(1,1,jj),fpsp1(1,jj),&
                    rhc,lprnt)
        call precpd(nxjp(j),nxp,lev,dta,del,prsl,psfc, &
                    qtc, qtr, ttc,           &
                    rlsp(1,jj), rhc, lprnt)
        do i=1,nxj
          rlsp(i,jj) = rlsp(i,jj) * 1000.         ! mm/call
        enddo
!jh        do k=1,lev
        do k=ktcup,lev
          kc=lev-k+1
          do i=1,nxj
            qt(i,k    ,jj) = qtc(i,kc)
            qt(i,k+lev,jj) = qtr(i,kc)
            tt(i,k    ,jj) = ttc(i,kc)
          enddo
        enddo
      endif
!
      if ( dodry ) then
         do k=1,lev
           dsigpp(k) = dsigma(k,1)+dsigma(k,2)/1000.
         enddo
         call drychk ( j,tt(1,1,jj),pk(1,1,jj),dsigpp,nxjp(j),nxp,lev,ndry(j) )
      endif
!
      if( upnor )then
        call nor_gwdp (j,nxjp(j),nxp,lev,                                 &
                  ut(1,1,jj),vt(1,1,jj),tt(1,1,jj),qt(1,1,jj),    &
                  plt(1,1,jj),pk(1,1,jj),pk2(1,1,jj),phi,dta,&
                  grav,rgas,sinl(j),cosl(j),drag_u,drag_v,cp)
        do k=1,lev
          avgdrag_u(j,k)=drag_u(k)
          avgdrag_v(j,k)=drag_v(k)
        enddo
      endif
!
!-------------------------------------------------------------------------
!     update radiation heating/cooling rate for next time step
!-------------------------------------------------------------------------
!     calculate ozone concentrtion
!-------------------------------------------------------------------------
      if (doo3l) then
!      if (myrank .eq. 0) print *,'### calculate ozone start !'
!      if (myrank .eq. 0) print *,'### doo3l=',doo3l

      if (ntoz .eq. 0) then
         if (doozon) then
         call rozone(nxjp(j),nxp,lev,plt(1,1,jj),o3l(1,1,jj),sinl(j),julian)
!         if (myrank .eq. 0) print *,'### rozone : set climate o3 values !'
         endif
      else
         call rozphys(nxjp(j),nxp,lev,dta,iter,xlat(j),julian,o3l(1,1,jj),&
                      tt(1,1,jj),plt(1,1,jj),ps(1,jj),myrank)
!         if (myrank .eq. 0) print *,'### rozphys ok!'
      endif ! for ntoz
      endif ! for doo3l


!     if (myrank .eq. 0) print *,'*** for o3l(1,1,jj):',' jj=',jj
!     call qmax2d(o3l(1,1,jj),1,1,nxj,lev)
!     if (myrank .eq. 0) print *,'o3l(1,60,jj)=',o3l(1,60,jj),' jj=',jj

!-----------------------------------------------------------------------
!   Radiation scheme
!-----------------------------------------------------------------------
      if (uprad .and. (irad .eq. 1))  then
!      if (myrank .eq. 0) print *,'### use radtn99 scheme'

         do 260 i = 1, nxj
         curate(i,jj) = rcup(i,jj) * 86400.0/dta
  260    continue

         call radtn99 ( fluxcl,ozon,nxjp(j),nxp,lev,ncld,lvlwx(j),julian        &
                    , stbo,s0,grav                                              &
                    , cp,ptrad,dsigma,sinl(j),cosz(1,jj),albedo2(1,jj),tg(1,jj) &
                    , curate(1,jj),pst(1,jj),plt(1,1,jj),tt(1,1,jj)             &
                    , qt(1,1,jj),o3l(1,1,jj)                                    &
                    , plcl(1,jj),cumtop(1,jj),ss(1,jj),rs(1,jj)                 &
                    , asol_clr(1,jj),olr_clr(1,jj),ss_clr(1,jj),rs_clr(1,jj)    &
                    , dtrad(1,1,jj),asr(1,j),alr(1,j)                           &
                    , asr_clr(1,j),alr_clr(1,j)                                 &
                    , xsr(1,j),xlr(1,j),acld(1,j),aflxd(1,j),aflxu(1,j)         &
                    , ilx(1,j),ibx(1,j),cofx(1,j),sdpbl(1,jj),ctot(1,jj)        &
                    , rld(1,jj),sld(1,jj),chig(1,jj),cmid(1,jj),clow(1,jj)      &
                    , asl(1,1,jj),atl(1,1,jj)                                   &
!--------------------------------------------------------------------------------
                    , fusl(1,1,jj),fdsl(1,1,jj)                                 &
                    , fuir(1,1,jj),fdir(1,1,jj)                                 &
                    , fuslr(1,1,jj),fdslr(1,1,jj)                               &
                    , fuirr(1,1,jj),fdirr(1,1,jj)                               &
                    , asl_clr(1,1,jj),atl_clr(1,1,jj)                           &
                    , clds(1,1,jj),rld_clr(1,jj),sld_clr(1,jj))
!--------------------------------------------------------------------------------
         do 265 i = 1, nxj
         asol(i,jj) = plcl(i,jj)
         olr(i,jj)  = cumtop(i,jj)
!cc      tg2 = tg(i,jj)*tg(i,jj)
!cc      rld(i,jj)  = stbo*(tg2*tg2) - rs(i,jj)
  265    continue
         endif  ! for uprad .and. irad=1
!--------------------------------------------------------------------------------
!   RRTMG scheme
!--------------------------------------------------------------------------------
      if (uprad .and. (irad .eq. 2))  then
      if ((isubc_lw .eq. 2) .or. (isubc_sw .eq. 2)) then
!ch       do i = 1 , nxj
          do i = 1 , nxdef(j)
             icsdsw(i) = ixseed(i,j,1)
             icsdlw(i) = ixseed(i,j,2)
          enddo
      endif  !isubc_lw

!      if (myrank .eq. 0) then
!          print *,'### use RRTMG scheme'
!          print *,'### before rrtmg : iter =',iter
!          print *,'### before rrtmg : tau   =',tau
!          print *,'### before rrtmg : solhr =',solhr
!          print *,'### before rrtmg : solcon=',solcon
!      endif
!--------------------------------------------------------------------------------
      call rrtmg                                                           &
!  ---  inputs:
           ( sigma,pst(1,jj),plt(1,1,jj),std(1,jj),                        &
             tt(1,1,jj),qt(1,1,jj),o3l(1,1,jj),sd(1,1,jj),tg(1,jj),        &
             slimsk   ,cice(1,jj),xtice(1,jj),                             &
             snr(1,jj),sncover(1,jj),snoalb(1,jj),z0(1,jj),                &
             alvsf(1,jj),alnsf(1,jj),alvwf(1,jj),                          &
             alnwf(1,jj),facsf(1,jj),facwf(1,jj),                          &
             curate(1,jj),icsdsw(nxjstart(j)),icsdlw(nxjstart(j)),         &
             sinl(j),cosl(j),xlat(j),xlonr(nxjstart(j),jj),jdat,d2r,xkapa, &
             ptrad,dtlw,dtsw,lsswr,lsswr,lssav,                            &
             nfxr,j,                                                       &
             nxp,nxjp(j),lev,ncld,lprnt,ipt,iter,solhr,solcon,             &
!  ---  outputs:
             asol(1,jj),olr(1,jj),ss(1,jj),rs(1,jj),                       &
             sld(1,jj),rld(1,jj),dtrad(1,1,jj),                            &
             ctot(1,jj),chig(1,jj),cmid(1,jj),clow(1,jj),                  &
             clds(1,1,jj),asl(1,1,jj),atl(1,1,jj),                         &
             fusl(1,1,jj),fdsl(1,1,jj),fuir(1,1,jj),fdir(1,1,jj),          &
             fuslr(1,1,jj),fdslr(1,1,jj),fuirr(1,1,jj),fdirr(1,1,jj),      &
             asl_clr(1,1,jj),atl_clr(1,1,jj),                              &
             asol_clr(1,jj),olr_clr(1,jj),ss_clr(1,jj),rs_clr(1,jj),       &
             sld_clr(1,jj),rld_clr(1,jj))
!      if (myrank .eq. 0) then
!          print *,'### rrtmg ok !!'
!      endif
      endif  ! for uprad .and. irad=2
!--------------------------------------------------------------------------------
!     weight back u and v by cosl/radus and
!     change real temp back to virtual potential temperature
!
      yy = cosl(j) / radus
      do 270 k = 1, lev
      do 270 i = 1, nxj
        ut(i,k,jj) = ut(i,k,jj)*yy
        vt(i,k,jj) = vt(i,k,jj)*yy
        tt(i,k,jj) = tt(i,k,jj)*(1.0+0.608*qt(i,k,jj))/pk(i,k,jj)
  270 continue
!
!
!    fill up negative moisture fields from level below
!
      call postq(nxjp(j),nxp,lev,ncld,dsigma,pst(1,jj),qt(1,1,jj),qbrrow(1,j),cosl(j))
!
  290 continue
!
!--------------------------------------------------------------------------------
!     update o3l to qt
!--------------------------------------------------------------------------------
      if(ncld.ge.3)then
      ntrac=3
      do k = 1, lev
         kk = (ntrac-1)*lev+k
      do jj=1, jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
      do i = 1, nxj
         qt(i,kk,jj) = o3l(i,k,jj)
      enddo
      enddo
      enddo
      endif
!
!      call mpe_global_sum(xkmd  ,lev,mpe_double)
!      call mpe_global_sum(dtcupd,lev,mpe_double)
!      call mpe_global_sum(dqcupd,lev,mpe_double)
!
!     compute accumulated rain (mm)
!
      do 300 jj =1, jlistnum
      j=jlist1(jj)
      nxj=nxdef_2d(j)
      do 300 i = 1,nxj
      totalp(i,jj) = (rcup(i,jj)  + rlsp(i,jj))*rainfc
      raincu(i,jj) = raincu(i,jj) + rcup(i,jj)*rainfc
      rainlp(i,jj) = rainlp(i,jj) + rlsp(i,jj)*rainfc
      raincu6(i,jj)= raincu6(i,jj)+ rcup(i,jj)*rainfc
      rainlp6(i,jj)= rainlp6(i,jj)+ rlsp(i,jj)*rainfc
      raintot(i,jj)= raintot(i,jj)+ totalp(i,jj)
  300 continue
!
!     reduce totalp in first hour forecast to avoid apin-up problem
!
      if ( tau .le. 1.001)  then
       dtau = dt/3600.0
!
       do 305 jj =1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
       do 305 i = 1,nxj
        totalp(i,jj) = (tau-dtau)* totalp(i,jj)
  305  continue
!
      endif
!
!      call mpe_unify(totalp,nx,my,2,mpe_double)
!
      if(ldiag.ge.1)then
!
!     print level-1 diagnostics
!                   (aps,psx,psn,nncup,nnshl,nndry,nnlsp,itlsp,arcup
!                   arlsp,xkmx,dtcupx and layer mean dtcup and dqcup )
!
!
      call mpe_global_sum(xkmd  ,lev,mpe_double)
      call mpe_global_sum(dtcupd,lev,mpe_double)
      call mpe_global_sum(dqcupd,lev,mpe_double)

!     call glbmean ( nx,my,my_max,cosl,pst ,aps   )
      call glbmean_2d ( cosl,pst ,aps   )
!
!     compute mean rain rate and change them from mm/call to mm/day
!
!     call glbmean ( nx,my,my_max,cosl,rcup,arcup   )
!     call glbmean ( nx,my,my_max,cosl,rlsp,arlsp   )
      call glbmean_2d ( cosl,rcup,arcup   )
      call glbmean_2d ( cosl,rlsp,arlsp   )
!
      arcup = arcup * 86400.0 / dta
      arlsp = arlsp * 86400.0 / dta
!
!     compute mean evaporation rate rate in mm/day
!
      do 310 jj =1, jlistnum
      j=jlist1(jj)
      nxj=nxdef_2d(j)
      do 310 i = 1,nxj
      albx(i,jj) = qflux(i,jj)/hltm*86400.0
  310 continue
!
!     call glbmean ( nx,my,my_max,cosl,albx,evapor )
      call glbmean_2d ( cosl,albx,evapor )
!
!     compute global moisture budget
!
      qglb = 0.0
      thda = 0.0
      tke=0.0
      tpe=0.0
!
      do j = 1, my*4
      wkj(j,1) = 0.
      enddo
!
      cosq=0.
      do 315 jj =1, jlistnum
      j=jlist1(jj)
      nxj=nxdef_2d(j)
      do 315 k = 1, lev
      do 315 i = 1, nxj
      dsigp = dsigma(k,1)*pst(i,jj) + dsigma(k,2)
      wkj(1,j) = wkj(1,j) + qt(i,k,jj)*dsigp*cosl(j)
      wkj(2,j) = wkj(2,j) + tt(i,k,jj)*dsigp*cosl(j) &
                 /(1.+0.608*qt(i,k,jj))
      wkj(3,j) = wkj(3,j) + (ut(i,k,jj)**2+vt(i,k,jj)**2)*dsigp &
                 *radsq/(2.*cosl(j))
      wkj(4,j) = wkj(4,j) + cp*tt(i,k,jj)*pk(i,k,jj)*dsigp &
                 *cosl(j)
      cosq = cosq + cosl(j)
  315 continue
!
      call mpe_unify(wkj,4,my,4,mpe_double)
!
      cosw=0.
      do 316 j = 1, my
      cosw = cosw + cosl(j)
      qglb = qglb + wkj(1,j)
      thda = thda + wkj(2,j)
      tke  = tke  + wkj(3,j)
      tpe  = tpe  + wkj(4,j)
 316  continue
      teng = tke+tpe
!     teng = teng/cosq
!     thda = thda/cosq
      teng = teng/cosw
      thda = thda/cosw
!
      qbrw = 0.0
!
!     call mpe_unify(qbrrow,ncld,my,2,mpe_double)
      call mpe_unify(qbrrow,ncld,my,4,mpe_double)
!
      do 320 j = 1, my
      qbrw = qbrw + qbrrow(1,j)
  320 continue
!
      cosq = cosw*nx
!
      qglb  = qglb * 100.0/grav/cosq
      qbrw  = qbrw * 100.0/grav/cosq
      rainbl= ( arcup + arlsp - evapor ) * dta / 86400.0
      qdiff = qgini - (qglb - rainbl - qbrw)
      engdiff=tengi-teng
      thdadif=thdai-thda
      qbrwtot=qbrwtot+qbrw
!
      nncup = 0
      nnshl = 0
      nndry = 0
!
      do 330 jj =1, jlistnum
      j=jlist1(jj)
      nncup = nncup + ncup(j)
      nnshl = nnshl + nshl(j)
      nndry = nndry + ndry(j)
  330 continue
!
      call mpe_global_sum(nncup,1,mpe_integer)
      call mpe_global_sum(nnshl,1,mpe_integer)
      call mpe_global_sum(nndry,1,mpe_integer)
!
      if ( uprad ) then
         call glbmean_2d ( cosl,asol,  topsd )
         call glbmean_2d ( cosl,olr,   toplu )
         call glbmean_2d ( cosl,ss,sfcsd )
         call glbmean_2d ( cosl,rs,sfclu )
      endif
!
      xoj = 1.0/cosw
!
      do 340 jj =1, jlistnum
      j=jlist1(jj)
      do 340 k = 1, lev
      nnlsp(k) = nnlsp(k) + nlsp(k,j)
      itlsp(k) = max (itlsp(k), ilsp(k,j))
      dtcupl(k) = dtcupl(k) + dtcupz(k,j)*cosl(j)
      dqcupl(k) = dqcupl(k) + dqcupz(k,j)*cosl(j)
  340 continue
!
      call mpe_global_sum(nnlsp ,lev,mpe_integer)
      call mpe_global_sum(dtcupl,lev,mpe_double)
      call mpe_global_sum(dqcupl,lev,mpe_double)
      call mpe_global_max(itlsp ,lev,mpe_integer)
!
      do 345 k = 1, lev
      dtcupl(k) = dtcupl(k) * xoj
      dqcupl(k) = dqcupl(k) * xoj * 1000.0
  345 continue
!
        xkmkk = -1.797693134862316d+308
        xkmk1 = -1.797693134862316d+308
        ixkmkk = 0
        ixkmk1 = 0
        jxkmkk = 0
        jxkmk1 = 0
        dtcupg = -1.797693134862316d+308
        idcupx = 0
        jdcupx = 0
        utx    = 0.0
        ikutx  = 0
        jutx   = 0
!
      do 360 jj =1, jlistnum
      j=jlist1(jj)
      nxj=nxdef_2d(j)
      if ( xkmkk .lt. xkmx(1,j) )  then
         xkmkk  = xkmx(1,j)
         ixkmkk = ipblmx(1,j)
         jxkmkk = j
      endif
      if ( xkmk1 .lt. xkmx(2,j) )  then
         xkmk1  = xkmx(2,j)
         ixkmk1 = ipblmx(2,j)
         jxkmk1 = j
      endif
      if ( dtcupg .lt. dtcupx(j) )  then
         dtcupg = dtcupx(j)
         idcupx = icupmx(j)
         jdcupx = j
      endif
         do 370 k = 1, lev
         do 380 i = 1,nxj
         speed     = sqrt(ut(i,k,jj)*ut(i,k,jj)+vt(i,k,jj)*vt(i,k,jj))
         dphi(i,k) = speed*radus/cosl(j)
  380    continue
         ikut = isamax ( nxjp(j), dphi(1,k), 1 )
         ikut=row_rank*nxp+ikut
         if ( dphi(ikut,k) .gt. utx )  then
           utx   = dphi(ikut,k)
           ikutx = ikut + nxp*(k-1)
           jutx  = j
         endif
  370    continue
  360 continue
!
      call mpe_global_maxloc(xkmkk, 2,ixkmkk,jxkmkk,idummy,idummy &
                            ,mpe_double)
      call mpe_global_maxloc(xkmk1, 2,ixkmk1,jxkmk1,idummy,idummy &
                            ,mpe_double)
      call mpe_global_maxloc(dtcupg,2,idcupx,jdcupx,idummy,idummy &
                            ,mpe_double)
      call mpe_global_maxloc(utx,   2,ikutx, jutx  ,idummy,idummy &
                            ,mpe_double)
!
      call maxpp2_2d( pst,pstx,ipstx,jpstx )
      call minpp2_2d( pst,pstn,ipstn,jpstn )
      call maxpp2_2d( hflux,hflmx,ihflx,jhflx )
      call maxpp2_2d( qflux,qflmx,iqflx,jqflx )
      call maxpp2_2d( tg,tgx,itgx,jtgx )

      dtradg = 0.0
      idradx = 0 ; jdradx = 0 ; kdradx = 0
      isign = 1
!
      do jj = 1, jlistnum
      j=jlist1(jj)
      nxj=nxdef_2d(j)

      do k = 1, lev
      do i = 1,nxj
       adtrad(i,k)=abs(dtrad(i,k,jj))
       dtradg= max(dtradg,adtrad(i,k))
      end do
      end do

      do k = 1, lev
      do i = 1,nxj
      if( dtradg.eq.adtrad(i,k) )then
        if ( dtrad(i,k,jj).lt.0 ) then
          isign =-1
        else
          isign = 1
        endif
!ch     idradx = nx*(k-1)+i
        idradx = nx*(k-1)+map2to1(i,j)
        jdradx = j
        kdradx = k
      endif
      end do
      end do

      end do
!
      call mpe_global_maxloc(dtradg,4,isign,idradx,jdradx,kdradx &
                            ,mpe_double)
!
      dtradg = dtradg*isign
!
      if(upnor)then
        call mpe_unify(avgdrag_u,my,lev,1,mpe_double)
        call maxpp2l( lev,my,avgdrag_u,dragmax,kmax,jmax )
        call minpp2l( lev,my,avgdrag_u,dragmin,kmin,jmin )
        if(myrank .eq. 0) then
          print *,'norographic gravity wave drag u_max,kmax,jmax,u_min,kmin,jmin'
          print *, dragmax,kmax,jmax,dragmin,kmin,jmin
        endif
        call mpe_unify(avgdrag_v,my,lev,1,mpe_double)
        call maxpp2l( lev,my,avgdrag_v,dragmax,kmax,jmax )
        call minpp2l( lev,my,avgdrag_v,dragmin,kmin,jmin )
        if(myrank .eq. 0) then
          print *,'norographic gravity wave drag v_max,kmax,jmax,v_min,kmin,jmin'
          print *, dragmax,kmax,jmax,dragmin,kmin,jmin
        endif
      endif
!fj
!
!      diagnoctics
!
      if(myrank .eq. 0) then
      print 8011, iter,aps,pstx,ipstx,jpstx,pstn,ipstn,jpstn
      print 8012, dtcupg,idcupx,jdcupx,dtradg,idradx,jdradx &
                , tgx,itgx,jtgx,utx,ikutx,jutx              &
                , xkmkk,ixkmkk,jxkmkk,xkmk1,ixkmk1,jxkmk1   &
                , hflmx,ihflx,jhflx,qflmx,iqflx,jqflx
      print 8013, nndry,nnshl,nncup,arcup,arlsp,evapor
      print 8014, nnlsp
      print 8015, itlsp
      print 8016, dtcupl
      print 8017, dqcupl
      if ( uprad ) print 8018, topsd,toplu,sfcsd,sfclu
      print *,' global moisture budget: qgini, qglb, qbrw, qdiff ='
      print *, qgini,qglb,qbrw,qdiff,qbrwtot
!
 8888 format('global average thdai,thda,thdadif=',e15.9,2(2x,e15.9))
      print *,'global average thdai,thda,thdadif='
      print *, thdai,thda,thdadif
!
 8889 format('global average tengi,teng,engdiff=',e15.9,2(2x,e15.9))
      print *,'global average tengi,teng,engdiff='
      print *, tengi,teng,engdiff
      endif
!
      if ( mod(tau+0.001,12.) .lt. dt/3600.)then
        if(myrank .eq. 0) then
         print 6001
        endif
 6001   format(4x,' qflux  ',' hflux  ',' dtcup18',' qtcup17', &
       ' dtcup16',' dtcup15',' dtcup14',' dtcup13',' dtcup12', &
       ' dtcup11')
!
!        call mpe_unify(hflux,nx,my,2,mpe_double)
!        call mpe_unify(qflux,nx,my,2,mpe_double)
!        call mpe_unify(dtcupz,lev,my,2,mpe_double)
!
        do 601 jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          xkmx(1,j)=0.
          xkmx(2,j)=0.
        do 600 i = 1, nxj
          xkmx(1,j)=xkmx(1,j)+qflux(i,jj)
          xkmx(2,j)=xkmx(2,j)+hflux(i,jj)
 600    continue
        xkmx(1,j)=xkmx(1,j)/nxjp(j)
        xkmx(2,j)=xkmx(2,j)/nxjp(j)
        if( myrank .eq. 0 ) &
         print 6000,j,xkmx(1,j),xkmx(2,j),(dtcupz(k,j),k=18,12,-1)
 6000   format('j=',i3,10f8.3)
 601    continue
      endif
!
      endif   !(end of ldiag.ge.1) 
!
!     print level-2 diagnostics
!                   (1) cupcwb diag, (2) radiation diag, (3) profile
!
!ldiag2 (start) ---
      if ( ldiag .ge. 2 )  then
!
       if(myrank .eq. 0) print *,'            cumulus diagnostics '
!
       nlcl_tmp(1:lev)  = 0.
       nnegl_tmp(1:lev) = 0.
       nosat_tmp(1:lev) = 0.
       nwork_tmp(1:lev) = 0.
       ntcup_tmp(1:lev) = 0.
       nflx_tmp(1:lev)  = 0.
!
       do 400 jj =1, jlistnum
        j=jlist1(jj)
       do 400 k = 1, lev
       nlcl_tmp(k)  = nlcl_tmp(k)  + nlcl(k,j)
       nnegl_tmp(k) = nnegl_tmp(k) + nnegl(k,j)
       nosat_tmp(k) = nosat_tmp(k) + nosat(k,j)
       nwork_tmp(k) = nwork_tmp(k) + nwork(k,j)
       ntcup_tmp(k) = ntcup_tmp(k) + ntcup(k,j)
       nflx_tmp(k)  = nflx_tmp(k)  + nflx(k,j)
  400  continue
!
       call mpe_global_sum(nlcl_tmp ,lev,mpe_integer)
       call mpe_global_sum(nnegl_tmp,lev,mpe_integer)
       call mpe_global_sum(nosat_tmp,lev,mpe_integer)
       call mpe_global_sum(nwork_tmp,lev,mpe_integer)
       call mpe_global_sum(ntcup_tmp,lev,mpe_integer)
       call mpe_global_sum(nflx_tmp ,lev,mpe_integer)
!
!
       if(myrank .eq. 0) then
        print 8020
        print 8022, ( k,nlcl_tmp(k),nnegl_tmp(k),nosat_tmp(k) &
                  ,   nwork_tmp(k),ntcup_tmp(k),nflx_tmp(k), k=1,lev )
       endif
!
       if ( uprad )                                                    &
        call diagrd ( ozon,nx,lev,my,my_max,njump,julian,stbo,s0,cosl &
                    , tt,qt,asol,olr,ss,rs,asr,alr,xsr,xlr            &
                    , acld,aflxd,aflxu,ncld )
!
!     print vertical profile at selected point (idg,jdg)
!
!      call mpe_unify(rcup,nx,my,2,mpe_double)
!      call mpe_unify(rlsp,nx,my,2,mpe_double)
!      call mpe_unify(pst,nx,my,2,mpe_double)
!      call mpe_unify(tg,nx,my,2,mpe_double)
!
!      rcuprr = rcup(idg,jdg)*86400.0/dta
!      rlsprr = rlsp(idg,jdg)*86400.0/dta
!      if(myrank .eq. 0) print 8024, idg,jdg,pst(idg,jdg),tg(idg,jdg),rcuprr,rlsprr
!
       do 425 jj = 1, jlistnum
        j=jlist1(jj)
        if(j .eq. jdg) then
         jjdg=jj
         do 420 k = 1, lev
          ud = ut(idg,k,jjdg)*radus/cosl(jdg)
          vd = vt(idg,k,jjdg)*radus/cosl(jdg)
          ttd= tt(idg,k,jjdg)
          qqd= qt(idg,k,jjdg)*1000.0
          dqcu= dqcupd(k)*1000.0
          work_pr2(k,1)=plt(idg,k,jjdg)
          work_pr2(k,2)=ud
          work_pr2(k,3)=vd
          work_pr2(k,4)=ttd
          work_pr2(k,5)=qt(idg,k,jjdg)
          work_pr2(k,6)=xkmd(k)
          work_pr2(k,7)=dtcupd(k)
          work_pr2(k,8)=dqcu
          work_pr2(k,9)=dtrad(idg,k,jjdg)
  420    continue
         call mpe_send_print(work_pr2, lev*9, jjdg, mpe_double)
        endif
  425  continue
!
       if(myrank .eq. 0) then
         call mpe_recv_print(work_pr2, lev*9, jjdg, mpe_double)
         print 8026, work_pr2
       endif
!
      endif
!ldiag2 (end) ---

 8011 format( 1x,'+iter =',i5,1x,' surf pres:  mean =', f8.2      &
            , ';  max = ', f8.2, 1x,'at (i,j)=', 2i5              &
            ,  ', min = ', f8.2, 1x,'at (i,j)=', 2i5 )
 8012 format( 1x,' max dtcup = ', f9.2, 1x, 'at (ik,j)= ',2i5     &
            , ';   max dtrad = ', f9.2, 1x, 'at (ik,j)= ',2i5,/   &
            , 1x,' max tg    = ', f9.2, 1x, 'at (i,j) = ',2i5     &
            , ';   max u-spd = ', f9.2, 1x, 'at (ik,j)= ',2i5,/   &
            , 1x,' max xkmkk = ', f9.2, 1x, 'at (i,j) = ',2i5     &
            , ';   max xkmk1 = ', f9.2, 1x, 'at (i,j) = ',2i5,/   &
            , 1x,' max hflux = ', f9.2, 1x, 'at (i,j) = ',2i5     &
            , ';   max qflux = ', f9.2, 1x, 'at (i,j) = ',2i5 )
 8013 format( 1x,' # of n-dry, n-shl, n-cup:', 3i6   &
            ,    ';      cup, lsp rain and evap (mm/day)= ', 3f9.4)
 8014 format( 1x,' # of n-lsp (saturated points), and max # of' &
            ,    ' iterations to converge:',/,1x, 30i6 )
 8015 format( 1x, 30i5 )
 8016 format( 1x, ' layer mean  dtcup and  dqcup  (k/day, g/kg/day):' &
            , /,1x, 30f5.1 )
 8017 format( 1x, 30f5.1)
 8018 format( 1x,' topsd, toplu, sfcsd, sfclu= ', 4f9.2 )
 8019 format( 1x,' global moisture budget: qgini, qglb, qbrw, qdiff =' &
            , 4e14.6,' qbrwtot=',e14.6 )
!
 8020 format('  lev  nlcl nnegl nosat nwork ntcup  nflx')
 8022 format( i4, 6i6)
 8024 format( 1x,'   point profile at (i,j)= (',i3,',',i3,') :  '      &
            ,' ps, tg, rcup, rlsp= ',f8.2,f8.2,2f8.3,/,2x              &
            ,'    pp      ut      vt      tt      qt    xkmd   dtcup'  &
            ,'   dqcup   dtrad')
 8026 format( 1x,f8.2,2f8.2,f8.2,f8.4,f8.2,3f8.3)
!
!     print level-3 dianostics:
!                   qprint rcup,rlsp,hflux,qflux,tg,tt(lev)
!
      if ( ldiag .ge. 3 )  then
      if ( myrank .eq.0 )  then
        print *,'2Dmpi version with array(partial_nx,partial_my) not supported for ldiag=3 !'
      endif
!       call qprnt3 (rcup,'rcup  ','mm/call',1,1,1,nx,my,my_max,1,1.0,0.0)
!       call qprnt3 (rlsp,'rlsp  ','mm/call',1,1,1,nx,my,my_max,1,1.0,0.0)
!       call qprnt3 (qflux,'qflux  ','w/m2',1,1,1,nx,my,my_max,1,1.0,0.0)
!       call qprnt3 (hflux,'hflux  ','w/m2',1,1,1,nx,my,my_max,1,1.0,0.0)
!       call qprnt3 (tg,'tg   ','(c) ',1,1,1,nx,my,my_max,1,1.0,-273.16)
!       call qprnt3 (tt,'tt   ','(c) ',1,1,lev,nx,my,my_max,lev,1.0,-273.16)
      endif
!
      return
      end
