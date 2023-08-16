!#define update_dp
      subroutine diabat ( docup,dodry,dolsp,dopbl,dorad,doshl,dograv,tofd      &
                    , nx,my,my_max,lev,ncld,nmcup,nmpbl,nmland,nmshl,cgw       &
                    , idg,jdg,ldiag,dt,tau,hours,julian,year,yrd               &
                    , frad,ozon,njump,itypbl,ktcup,ktpbl,ktshl,grav            &
                    , rgas,cp,stbo,s0,evaprh,hltm,ptop,sigma,dsigma,il,ib,cof  &
                    , xlat,xlon,sgeo,z0,alb,land,ocean,ice,snr                 &
                    , tg,tgclim,curate,plcl,cumtop,totalp,raintot,raincu       &
                    , rainlp,raincu6,rainlp6,raincu3,rainlp3,raincu1,rainlp1   &
                    , hflux,qflux,ustar,tstar,qstar,e                          &
                    , eps,o3l,dtrad,ss,rs,plt,pk,pk2,ps,up,vp,ttp,qp           &
                    , pst,ut,vt,tt,qt,gwclim,tice,hice,qgini                   &
                    , thdai,tengi,acld,std,asol,olr,drag                       &
                    , ugws,vgws,sdpbl,t2,q2,rh2,rh10,u10,v10,gfx               &
                    , fm,fh,fm10,fh2,srflag                                    &
                    , rld,km_soil,smc,stc,canopy,runoff                        &
                    , sigmaf,istyp,ivegtyp,wlt,ref,tsat,dfkt,xktk,dfk          &
                    , ftp,fqp,fpsp,ftp1,fqp1,fpsp1,deltaq,cnvwr,cnvcr,sd       &
                    , shdmax,shdmin,snoalb                                     &
                    , slopetyp,sld,slc,zice,cice,xtice,sncover,sndepth         &
                    , ctot,chig,cmid,clow,hpbl,asl,atl,cosz                    &
                    , nmgwor,nmgwcv,hprime_b,mtnvar,docgrav,nmmiph             &
!--------------------------------------------------------------------------------
                    , fusl,fdsl,fuir,fdir                                      &
                    , fuslr,fdslr,fuirr,fdirr                                  &
                    , asl_clr,atl_clr,clds                                     &
                    , ss_clr,rs_clr,asol_clr,olr_clr,sld_clr,rld_clr           &
                    , alvsf,alvwf,alnsf,alnwf,facsf,facwf                      &
                    , idtg,doo3l,nfxr,sfalb,sfemis,isot,ivegsrc                &
! sit
                    , itimestep,lrun_sitvdiff,ic_sit                           &
!xb110>
                    , flash,tsflw,vvel)
!xb110<
!--------------------------------------------------------------------------------
!#######################################################################
!
!     driver for all physical parametrization of diabatic processes
!
!
!     parameters
!
!     docup  : logical variable for including cumulus parameterization
!     dodry  : logical variable for including dry convection
!     dolsp  : logical variable for including large scale precipitation
!     dopbl  : logical variable for including pbl parameterization
!     dorad  : logical variable for including radiation parameterization
!     doshl  : logical variable for including shallow convection calcul
!     dograv : logical variable for including gravity wave drag
!     dograv : logical variable for including convective gravity wave drag
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
!   lrun_sitvdiff: logical variable for run sit_vdiff
!   ic_sit   : the 'ic_sit'th times run sit_vdiff
!
! modify to f90 by C-H Lee and sort by River Chen in 2015
!
!#######################################################################
!
      use mpe
      use rank
      use index
      use const,                 ONLY:do_sit,ldailyFCTsst,dailyClm_option,      &
                                      pdfcloud,cmbk,cgwd, fsit, dosppt, doshum, dossst, &
                                      use_zmtnblck,ldailyFCTicesndpt,dSITdt_intv, &
                                      weightSIT,bckfile,ggdef,doclx,doslavepp,    &
                                      RTYPE,qmin
      use mod_sitgrid
      USE mod_sit_vdiff,         ONLY:sit_vdiff,ctfreez
      USE mod_sit_control,       ONLY:ftrigsit,ltrigsit,lsitstart,lsftobswt &
                                     ,timebl_option,timebl_start
      USE mod_eos_ocean,         ONLY:AirVaporPressure,CalcSm
      USE mod_sst,               ONLY:time_weights,now1,now2,wgto1,wgto2 &
                                     ,obswtbnmw1,obswtbnmw2,obswtbwgt1   &
                                     ,obswtbwgt2,obswtbold,obswtbnow     &
                                     ,obswtbnew,dFCTsstdt                &
                                     ,tseadiffFCT,tseadiffFCT24,dtseadt  &
                                     ,tseanow,tseaold
      USE mo_netcdf,             ONLY:lkvl
!-----------------------------------------------------------------------
      use radn
      use physpara
!
      use physcons, only :con_rd,con_fvirt,con_rerth,con_rv,con_pi
! for slavepp
      use phygrid,  only :dtcup,ducup,dvcup,dtshl,dushl,dvshl,dtlsp
! for land_noah_new
      use namelist_soilveg, only :MAX_SLOPETYP,MAX_SOILTYP,MAX_VEGTYP
      use mod_stochastic_physics, only : sppt3d, shum3d, ssst3d
!-----------------------------------------------------------------------
      implicit  none
!-----------------------------------------------------------------------
      integer nfxr, ntrac, kk, nk, n
      real    dtlw,dtsw,solhr,rsolhr
!
! for land_noah_new
       real      sfalb(nxp,my_max),sfemis(nxp,my_max)
       integer   isot,ivegsrc
!-----------------------------------------------------------------------
      integer   nx,my,my_max,lev,ncld,nmcup,nmpbl,nmland,nmshl,idg,  &
                jdg,ldiag,julian,njump,itypbl,ktcup,ktpbl,ktshl,     &
                km_soil

      logical   docup,dodry,dolsp,dopbl,dorad,doshl,dograv,ozon,     &
                land(nxp,my_max),ocean(nxp,my_max),ice(nxp,my_max),  &
                docgrav,tofd

      real      tice,hice,qgini,thdai,tengi,ptop,                    &
                hltm,evaprh,s0,stbo,cp,rgas,grav,frad,               &
                hours,tau,dt,cgw

      integer   il(nxp,4),ib(nxp,4)

      real      cof(nxp*3,4),xlat(my),                                    &
                xlon(nx,my_max),z0(nxp,my_max),                           &
                alb(nxp,my_max),snr(nxp,my_max),tg(nxp,my_max),           &
                tgclim(nxp,my_max),curate(nxp,my_max),plcl(nxp,my_max),   &
                cumtop(nxp,my_max),totalp(nxp,my_max),raincu(nxp,my_max), &
                hflux(nxp,my_max),qflux(nxp,my_max),ustar(nxp,my_max),    &
                tstar(nxp,my_max),qstar(nxp,my_max),                      &
                e(nxp,lev,my_max),eps(nxp,lev,my_max),                    &
                dtrad(nxp,lev,my_max),ss(nxp,my_max),                     &
                rs(nxp,my_max),plt(nxp,lev,my_max),                       &
                rainlp(nxp,my_max),                                       &
                gwclim(nxp,my_max),acld(lev,my),std(nxp,my_max),          &
                asol(nxp,my_max),olr(nxp,my_max),drag(nxp,lev,my_max),    &
                ugws(nxp,my_max),vgws(nxp,my_max),                        &
                raintot(nxp,my_max),t2(nxp,my_max),rh2(nxp,my_max),       &
                rh10(nxp,my_max),                                         &
                q2(nxp,my_max),fm(nxp,my_max),fh(nxp,my_max),             &
                fm10(nxp,my_max),fh2(nxp,my_max),srflag(nxp,my_max),      &
                u10(nxp,my_max),v10(nxp,my_max),hpbl(nxp,my_max),         &
                raincu6(nxp,my_max),rainlp6(nxp,my_max),                  &
                raincu3(nxp,my_max),rainlp3(nxp,my_max),                  &
                raincu1(nxp,my_max),rainlp1(nxp,my_max)
      real(kind=RTYPE) qt(nxp,lev*ncld,my_max),qp(nxp,lev*ncld,my_max),   &
                       up(nxp,lev,my_max),vp(nxp,lev,my_max),             &
                       ttp(nxp,lev,my_max),o3l(nxp,lev,my_max),           &
                       sgeo(nxp,my_max),ps(nxp,my_max),pst(nxp,my_max),   &
                       sdpbl(nxp,my_max),sigma(lev+1,2),dsigma(lev,2),    &
                       ut(nxp,lev,my_max),vt(nxp,lev,my_max),             &
                       tt(nxp,lev,my_max),pk(nxp,lev,my_max),             &
                       pk2(nxp,lev,my_max)
!soil (2005/01/12)
      integer,  parameter :: ntype=9, ngrid=22
      integer   istyp(nxp,my_max),ivegtyp(nxp,my_max)

      real      smc(nxp,km_soil,my_max),stc(nxp,km_soil,my_max),          &
                canopy(nxp,my_max),runoff(nxp,my_max),                    &
                sigmaf(nxp,my_max),rld(nxp,my_max),                       &
                wlt(MAX_SOILTYP),ref(MAX_SOILTYP),tsat(MAX_SOILTYP),      &
                dfkt(ngrid,ntype),xktk(ngrid,ntype),dfk(ngrid,ntype)

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
! --- for radupdat
      integer idat(8),jdat(8)

! --- for new albedo
      real alvsf(nxp,my_max),alvwf(nxp,my_max),alnsf(nxp,my_max), &
           alnwf(nxp,my_max),facsf(nxp,my_max),facwf(nxp,my_max)
!for gravity wave drag ====== #
      integer nmgwor,nmgwcv,mtnvar
      real hprime_b(nxp,mtnvar,my_max)
      real tt_bfcnv(nxp,lev)
      real prsi(nxp,lev+1)
      real utgwc(nxp,lev),vtgwc(nxp,lev),                                  &
           dudtc(nxp,lev),dvdtc(nxp,lev),dtdtc(nxp,lev),dqdtc(nxp,lev),    &
           prslk(nxp,lev)
      real oc(nxp),theta(nxp),gamma(nxp),sigmaog(nxp),elvmax(nxp),hprime(nxp),    &
           dlength(nxp),cldf(nxp),cumabs(nxp),work3(nxp),tauctx(nxp),taucty(nxp), &
           facg(lev)
      real oa4(nxp,4),clx(nxp,4),cgwf(2),cdmbgwd(2)
      real ograv
      integer kpbl(nxp,my_max)
      integer latg
      real eng0,eng1

! --- for random number generator (thread safe mode)
      integer ixseed(nx,my,2)
      integer icsdlw(nxp),icsdsw(nxp)

! --- for rrtmg : input
!     logical lsswr,lslwr,lssav,lprnt
      logical lsswr,lslwr,lssav
      real    xlonr(nxp,my_max),sld_adj(nxp),rld_adj(nxp),ss_adj(nxp), &
              rs_adj(nxp),tsflw(nxp,my_max),rstd(nxp)
! --- for slavepp
      real    dttmp,dutmp,dvtmp

! --- new variables setting :
      integer(kind=8)  :: idtg
      logical :: doo3l
      integer :: ipt,jpt
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

      ! for sppt
      real(kind=RTYPE) :: ut_save_sppt(nxp,lev,my_max)     
      real(kind=RTYPE) :: vt_save_sppt(nxp,lev,my_max)        
      real(kind=RTYPE) :: tt_save_sppt(nxp,lev,my_max)        
      real(kind=RTYPE) :: qt_save_sppt(nxp,lev*ncld,my_max)
     !real(kind=RTYPE) :: qt_save_shum(nxp,lev*ncld,my_max)
      real(kind=RTYPE) :: tg_save_ssst(nxp,1,my_max)
#ifdef VERBOSE
      real(kind=RTYPE) :: ut_update(nxp,lev,my_max)     
      real(kind=RTYPE) :: vt_update(nxp,lev,my_max)        
      real(kind=RTYPE) :: tt_update(nxp,lev,my_max)        
      real(kind=RTYPE) :: qt_update(nxp,lev*ncld,my_max)

      real(kind=RTYPE) :: ut_pbl(nxp,lev,my_max)
      real(kind=RTYPE) :: vt_pbl(nxp,lev,my_max)
      real(kind=RTYPE) :: tt_pbl(nxp,lev,my_max)
      real(kind=RTYPE) :: qt_pbl(nxp,lev*ncld,my_max)

      real(kind=RTYPE) :: ut_cmls(nxp,lev,my_max)     
      real(kind=RTYPE) :: vt_cmls(nxp,lev,my_max)        
      real(kind=RTYPE) :: tt_cmls(nxp,lev,my_max)        
      real(kind=RTYPE) :: qt_cmls(nxp,lev*ncld,my_max)

      real(kind=RTYPE) :: ut_sppt(nxp,lev,my_max)     
      real(kind=RTYPE) :: vt_sppt(nxp,lev,my_max)        
      real(kind=RTYPE) :: tt_sppt(nxp,lev,my_max)        
      real(kind=RTYPE) :: qt_sppt(nxp,lev*ncld,my_max)
#endif
      real :: dtradc(nxp,lev,my_max)
      real :: dtradn(nxp,lev)
      integer :: zmtnblck(nxp)
      real :: ru,vru,upert,vpert,tpert,qpert,qnew,dtdtr

      integer itimestep,ii

      ! for MP WSM6 & Thompson
      logical uni_cloud,lmfshal,lmfdeep2
      real    sr(nxp,my_max)
!---------------------------------------------------------------------------
      real      drag_u(lev),drag_v(lev)
      real      fnor
      logical   donor,upnor
      data      donor/.true./,fnor/0.5/

! for microphysics
      real area

!#######################################################################
!
!     local logical variables and work arrays
!
      logical   fluxcl,doozon,uprad
      integer   ijdg(my),ipblmx(2,my),itlsp(lev),nnlsp(lev),ilsp(lev,my),&
                nlsp(lev,my),ncup(my),ndry(my),nshl(my),icupmx(my),      &
                nlcl(lev,my),nnegl(lev,my),nosat(lev,my),nwork(lev,my),  &
                ntcup(lev,my),nflx(lev,my),lvlwx(my_max),ilx(nxp,my_max),&
                ibx(nxp,my_max)

      real      cosl(my),sinl(my),cosz(nxp,my_max),                  &
                rcup(nxp,my_max),rlsp(nxp,my_max),               &
                asr(lev,my),alr(lev,my),xsr(lev,my),xlr(lev,my),         &
                aflxd(lev+2,my),aflxu(lev+2,my),                         &
                dtcupx(my),dtcupz(lev,my),dqcupz(lev,my),dtcupd(lev),    &
                dqcupd(lev),dtcupl(lev),dqcupl(lev),xkmx(2,my),xkmd(lev),&
                albx(nxp,my_max),cofx(nxp*3,my_max),dphi(nxp,lev)
      real(kind=RTYPE) phi(nxp,lev)

      real      wkj(4,my),dsigpp(lev),qt_diff(ncld)

      integer   nlcl_tmp(lev) ,nnegl_tmp(lev),nosat_tmp(lev),    &
                nwork_tmp(lev),ntcup_tmp(lev),nflx_tmp(lev)

      real      adtrad(nxp,lev),work_pr1(9),work_pr2(lev,9)
!--------
! for ncld=2
      real,     parameter :: dxmax=-16.118095651,dxmin=-9.800790154, &
!      real,     parameter :: dxmax=-17.261145789,dxmin=-12.465355243,&
                             dxinv=1.0/(dxmax-dxmin)
!     parameter (rhzbot=0.85, rhztop=0.85)

      real      work1(nxp),work2(nxp),rhc(nxp,lev),rhckt,psautco(nxp)
      real      rhc_mp(nxp,lev)  !for GFDL MP
      real      del(nxp,lev),prsl(nxp,lev),psfc(nxp)
      real      qtc(nxp,lev), qtr(nxp,lev), ttc(nxp,lev)
      real      ftp(nxp,lev,my_max), fqp(nxp,lev,my_max), fpsp(nxp,my_max)
      real      ftp1(nxp,lev,my_max), fqp1(nxp,lev,my_max), fpsp1(nxp,my_max)
!-------
!for pdf cloud
      integer   kdt
      real      sup
      real      cnvw(nxp,lev),cnvc(nxp,lev),deltaq(nxp,lev,my_max)
      real      cnvwr(nxp,lev,my_max),cnvcr(nxp,lev,my_max)

! vertcal rhc
!      real      ct,cs,px

      logical lprnt
!--------
!
! for sascnv
      integer   kbot(nxp,my_max),ktop(nxp,my_max),kuo(nxp,my_max),  &
                islimsk(nxp)
      real      sl(lev),delcup(lev),slimsk(nxp)
      real      dotc(nxp,lev),phil(nxp,lev),utc(nxp,lev),vtc(nxp,lev)

!ch   real      cldwrk(nxp,my_max),sd(nxp,lev+1,my_max),xkt2(nx)
      real      cldwrk(nxp,my_max),                     xkt2(nx)
      real(kind=RTYPE) sd(nxp,lev+1,my_max),vvel(nxp,lev,my_max)

! for new shlcon
      real      rcup2(nxp)
! for scale-aware convection
      real      garea(nxp),tpr,tem1,tem2,jup,jdn,tpi
! for wsm6 & thompson
      integer   nmmiph
      real      phii(nxp,lev+1)
      real      qti(nxp,lev),qtrw(nxp,lev),qtsw(nxp,lev),qtgl(nxp,lev)
      real      icem,ntnc(nxp,lev,2) !1:ice, 2:liquid
      real      ice00(nxp,lev)
      real      qni
! for updating low boundary condition
      integer   ls(nxp,my_max)
      real      sstc(nxp,my_max),z0ocn(nxp,my_max)
      logical   doclxu,iceold(nxp,my_max)


!CWB 2007-09-27 for random number seed >>>
      real*8    rtc,rsecond
      integer   isize(2)
! CWB <<<


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
                tke,   tpe,    cosq,   dsigp,  cosw,   teng,   rainbl,       &
                engdiff,thdadif,topsd,  toplu,  sfcsd,  sfclu, xoj,          &
                xkmkk, xkmk1,  dtcupg, utx,    speed,  pstx,   pstn,  hflmx, &
                qflmx, tgx,    dtradg, dragmax,dragmin,rcuprr, rlsprr,ud,    &
                vd,    ttd,    qqd,    dqcu,   deg_ju, arg,    tem

!xb110>
!for new precpd & nTDK
      real      u0(nxp,lev),v0(nxp,lev),t0(nxp,lev)
      real(kind=RTYPE) q0(nxp,lev*ncld)
      real      upp(nxp,lev),vpp(nxp,lev),tpp(nxp,lev),ttpp(nxp,lev)
      real      pltp(nxp,lev)
      real(kind=RTYPE) pkp(nxp,lev),pk2p(nxp,lev)
!for lightning
      real      flash(nxp,my_max)        !flash density (unit in flashes km^-2 day^-1)
      real      ztenh(nxp,lev),zqenh(nxp,lev),rho(nxp,lev)              &
               ,snow_flxn(nxp,lev),ptun(nxp,lev),pqun(nxp,lev)          &
               ,cnvwn(nxp,lev)
      real      snow_flx(nxp,lev),ptu(nxp,lev)           &
               ,pqu(nxp,lev)
      integer   kbotc(nxp,my_max), ktopc(nxp,my_max)
!xb110<
!ps
      logical lrun_sitvdiff
      integer ic_sit
      real liw, t_surf
      real sitlat(nxp)
      real sitlon(nxp,my_max)
      integer istep,icurrenttau
      real rnuw,z0w,d0w,Pairvapor,rhoa,ps1w,ustra,vstra
      integer*8 idtg_sitvdiff
      real tauleft, dtsit
      real, dimension(:), allocatable :: sstm, rstm,hfluxtm,qfluxtm, &
                                          u10tm,v10tm, rlsptm, rcuptm, &
                                        ustartm, t2tm,  rh2tm,  psttm, &
                                         cicetm,snrtm, zicetm,xticetm, &
                                       obswtbtm, tgtm
      character*12 cdtg
      integer yr, mo, dy, hr, mn, leap, yrd, year
      real tauhr
      real dtx_tau,dtaup
      INTEGER, PARAMETER :: nerr = 6
!xb110>
      ztenh = 0.
      zqenh = 0.
      rho = 0.
      snow_flx  = 0.
      ptu = 0.
      pqu = 0.
      cnvwn = 0.
!xb110<
!ps
!CWB2015 
      kuo=0
      rcup2=0.
      islimsk=0
      kpbl=1

!CWB2016
      lssav=.false.
      doclxu=.false.
      icsdsw=0
      icsdlw=0
      rld_adj=0.
      sld_adj=0.
      ss_adj =0.
! for MP WSM6 & Thompson
      uni_cloud=.false. !if using SHOC scheme, it should be .true.
      lmfshal=( nmshl .eq. 2 .or. nmshl .eq. 3 ) ! .true. if using mass-flux shallow convection
      lmfdeep2=( nmcup .eq. 6 ) ! .true. if using scale-aware deep con


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
      sup=  1.0
      nxmy = nx * my
      nxlev= nx * lev
      levmy= lev* my
      radus = 6371000.
      radsq = radus**2
      tpi   = 4.0*atan(1.0)
      d2r   = tpi / 180.0
      tpr   = 2.*tpi*radus
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
      dta    = dt
      rainfc = 1.0
      if (itimestep .le. 1) then
         doozon = .true.
         kdt    = 1
      else
         doozon = .false.
         kdt    = 0
      endif

!------------------------------------------------------------------------------
!     set hours, iter, icrad, julian, uprad, doozon
!------------------------------------------------------------------------------
      rsolhr = hours
      hours = hours + dt/3600.0
      if ( hours .gt. 24.0 )  then
         hours = mod ( hours,24.0 )
         julian= julian + 1
         if ( julian .gt. yrd ) julian = julian - yrd
         doozon = .true.
         doclxu = .true.
      endif
      leap = mod ( year , 4 )
      yrd = 365
      if ( leap .eq. 0 ) yrd = 366
!
      icrad = frad*3600.0/dt + 0.0001 ! frad =1.0 set in block.f
      iter  = tau*3600.0/dt - 1. + 0.0001
      uprad = .false.
      if ( (mod(iter,icrad).eq.0) .or. iter.eq.1 )  uprad = .true.
      doozon = doozon .and. dorad
      uprad  = uprad  .and. dorad
!
! update low boundary condition
! 
      if ( doclxu .and. doclx ) then
        if (myrank.eq.0)                                                &
           print *,'update low boundary condition at tau= ',tau
        iceold=ice
        z0ocn=z0
!     read climate data
        call readclx( nx,my,my_max,julian,land,ocean,ice,tgclim,gwclim  &
                   ,z0,alb,sstc,bckfile,sigmaf,istyp,ivegtyp,ls         &
                   ,shdmax,shdmin,slopetyp,snoalb,ggdef,isot,ivegsrc )
!
!     read new albedo
!
        if (irad .eq. 2) then
          call readalb(bckfile,nx,my,my_max,julian,ggdef,         &
                     alvsf,alvwf,alnsf,alnwf,facsf,facwf)

          do jj=1,jlistnum
            j=jlist1(jj)
            nxj=nxdef_2d(j)
            do i=1,nxj
              alvsf(i,jj)=alvsf(i,jj)*0.01
              alvwf(i,jj)=alvwf(i,jj)*0.01
              alnsf(i,jj)=alnsf(i,jj)*0.01
              alnwf(i,jj)=alnwf(i,jj)*0.01
              facsf(i,jj)=facsf(i,jj)*0.01
              facwf(i,jj)=facwf(i,jj)*0.01
            enddo
          enddo
        endif
!
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
        do i=1,nxj
          if(ls(i,jj).eq.0) then
!---------------------------------------------------------------------
! (1)  retain surface roughness over ocean
!---------------------------------------------------------------------
            z0(i,jj)=z0ocn(i,jj)
!---------------------------------------------------------------------
! (2)  tg replaced by climate sea surface temperature
!---------------------------------------------------------------------
            if (ocean(i,jj)) tg(i,jj)=sstc(i,jj)
!---------------------------------------------------------------------
! (3)  set ice thickness => not for couple
!---------------------------------------------------------------------
            if(iceold(i,jj)       .and. .not. ice(i,jj)) then
              zice(i,jj)=0.
              cice(i,jj)=0.
              snr(i,jj) =0.
              sndepth(i,jj)=0.
              sncover(i,jj)=0.
              z0(i,jj)=ustar(i,jj)*ustar(i,jj)*0.014/grav
            endif
            if(.not. iceold(i,jj) .and. ice(i,jj)) then
              xtice(i,jj)=tg(i,jj)
              zice(i,jj)=0.15 ! from himin in sfc_sice 
              cice(i,jj)=0.15 ! from cimin in sfc_sice 
              snr(i,jj) =15.
              sndepth(i,jj)=snr(i,jj)*8.
              sncover(i,jj)=min(1., snr(i,jj)/400.)
              z0(i,jj)=0.0002 ! set new ice point to 0.0002
            endif
          endif ! if(ls(i,jj).eq.0) then
        enddo
        enddo

      endif ! doclxu
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

      do k = 1, lev
        do i = 1, nxp
          utgwc(i,k)  = 0.
          vtgwc(i,k)  = 0.
!for pdfcloud
          cnvw(i,k) = 0.
          cnvc(i,k) = 0.
!for hydrometeor
          qtr(i,k)  = 0.
          qtc(i,k)  = 0.
          qti(i,k)  = -999.9
          dudtc(i,k) = 0.
          dvdtc(i,k) = 0.
          dtdtc(i,k) = 0.
          dqdtc(i,k) = 0.
        enddo
      enddo
!
      do jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
      do i = 1, nxp
       rcup(i,jj)  = 0.0
       rlsp(i,jj)  = 0.0
!       cldwrk(i,k) = 0.0
       cldwrk(i,jj) = 0.0
!      cosz(i,jj)  = 0.0
       xmu(i,jj)  = 0.0
! for wsm6
       sr(i,jj)  = 0.0
      enddo
      enddo

      do j = 1, my
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
!
       ijdg(j) = 0
!ch    qbrrow(1,j)= 0.0
!ch    qbrrow(2,j)= 0.0
!!       qbrrow(1:ncld,j)= 0.0
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
!-------------------------------------------------------------------------
!     if radiation is to be called, compute cos of solar zenith angular
!                njump, il, ib, and cof for different meridional zones
!-------------------------------------------------------------------------
      call prerrtmg(nx,my,my_max,idtg,tau,dt,hours,frad,uprad,        &
                    isubc_sw,isubc_lw,d2r,xlon,myrank,me,             &
                    idat,jdat,solhr,dtsw,dtlw,lsswr,lslwr,            &
                    slag,sdec,cdec,solcon,                            &
                    xlonr,ixseed)
      year  = jdat(1)
!-------------------------------------------------------------------------
! cosz was modified to be an average of the  calling period(1 hour fo
      if ( uprad .and. (irad .eq. 1))  then
        call coszenpm ( nx,my,my_max,julian,solhr,sinl,cosl,xlonr,    &
                        frad,slag,sdec,cdec,cosz )
!
!!ocl scalar,nounroll
         do 160 jj = 1, jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)

         abxlat = abs(xlat(j))
         njump1 = njump + 1
         if ( mod(nxp,njump1) .ne. 0 )  njump1 = njump
         njump2 = njump + 2
         if ( mod(nxp,njump2) .ne. 0 )  njump2 = njump1
         njump3 = njump + 3
         if ( mod(nxp,njump3) .ne. 0 )  njump3 = njump2

         if( lreduce.eq.1 ) then
           njump1 = njump
           njump2 = njump
           njump3 = njump
         endif

         if ( abxlat .le. 20.0 )  then
!ch         lvlwx(j) = nxjp(j)/njump
            lvlwx(jj) = nxp/njump
            nny = 1
         else if ( abxlat .le. 60.0 )  then
!ch         lvlwx(j) = nxjp(j)/njump1
            lvlwx(jj) = nxp/njump1
            nny = 2
         else if ( abxlat .le. 80.0 )  then
!ch         lvlwx(j) = nxjp(j)/njump2
            lvlwx(jj) = nxp/njump2
            nny = 3
         else
!ch         lvlwx(j) = nxjp(j)/njump3
            lvlwx(jj) = nxp/njump3
            nny = 4
         endif
!
         if( lreduce.eq.1 ) then
           call splinc (lvlwx(jj),nxp,ilx(1,jj),ibx(1,jj),cofx(1,jj))
         else
!          do 150 i  = 1, nx
           do 150 i  = 1, nxp
            ilx(i,jj)  = il(i,nny)
            ibx(i,jj)  = ib(i,nny)
  150      continue
           do 155 i  = 1, lvlwx(jj)*3
            cofx(i,jj) = cof(i,nny)
  155      continue
         endif
  160    continue
      endif

!
! calculate xmu,
! the instantaneous zenith angular at the 'hours'
!------------------------------------------------------------------------------
!        call ccoszen ( nx,my,my_max,julian,hours,xlat,xlon,xmu )
!------------------------------------------------------------------------------

! transfer xmu to a fraction of average value(cosz) during frad
!      do 170 jj = 1, jlistnum
!       j=jlist1(jj)
!       nxj=nxdef_2d(j)
!      do 170 i = 1, nxj
!        if(xmu(i,jj).gt.0.0001.and.cosz(i,jj).gt.0.0001) then
!          xmu(i,jj) = xmu(i,jj) / cosz(i,jj)
!        else
!          xmu(i,jj)   = 0.
!        endif
! 170  continue

!     initial albx by climate values of gwet and alb, while it may
!     be updated in grdcon according ground conditions
!
      do jj = 1, jlistnum
         j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i = 1, nxj
          albx(i,jj)  = alb(i,jj)
          albedo2(i,jj)  = alb(i,jj)
        enddo
      enddo
!
!-----------------------------------------------------------------------
      if(ntoz.gt.0)then
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do k = 1, lev
            kk = (ntoz-1)*lev+k
            do i = 1, nxj
              o3l(i,k,jj) = qt(i,kk,jj)
            enddo
          enddo
        enddo
      endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                      c
!     begin big j-loop for diabatic calculation in each latitude ring  c
!                                                                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    do 290 jj =1, jlistnum
      j=jlist1(jj)
      nxj=nxdef_2d(j)

      ! for scale-aware
      tem1      = tpr*cosl(j)/float(nxdef(j))
      jup       = min(j+1,my)
      jdn       = max(j-1, 1)
      tem2      = radus*0.5*abs(xlat(jup)-xlat(jdn))*d2r
      do i = 1,nxj
        work1(i)    = (log(cosl(j) / (nxdef(j)*my)) - dxmin) * dxinv
        work1(i)    = max(0.0, min(1.0,work1(i)))
        work2(i)    = 1.0 - work1(i)
        garea(i)    = tem1*tem2
        if(land(i,jj))slimsk(i)=1
        if(ocean(i,jj))slimsk(i)=0
        if(ice(i,jj))slimsk(i)=2
        if(land(i,jj))islimsk(i)=1
        if(ocean(i,jj))islimsk(i)=0
        if(ice(i,jj))islimsk(i)=2
      enddo

    !    compute new time level p**kapa quantites
    !  plt= new odd level pressure
    !  pk= (plt/1000)**capa
      call prexp_hybrid_cwb ( nxjp(j),nxp,lev,ptop,sigma,pst(1,jj), &
                          pk(1,1,jj),pk2(1,1,jj),plt(1,1,jj) )
      call prexp_hybrid_cwb ( nxjp(j),nxp,lev,ptop,sigma,ps(1,jj),  &
                          pkp,pk2p,pltp )
!
!     hydrostatic equation
      do i = 1, nxj
        phi(i,lev) = sgeo(i,jj)+ cp*tt(i,lev,jj)*(pk2(i,lev,jj)-pk(i,lev,jj))
      enddo

      do k = lev-1, 1, -1
        do i = 1, nxj
          phi(i,k) = phi(i,k+1) + cp*( tt(i,k,jj)  *( pk2(i,k,jj) - pk(i,k,jj)  )   &
                                     + tt(i,k+1,jj)*( pk(i,k+1,jj)- pk2(i,k,jj) ) )
        enddo
      enddo

    !-----------------------------------------------------------------------------
    !  deweight u,v by cosl/radus, and
    !  change t from virtual potential temperature to real temperature
    !  change ttp from virtual potential temperature to potential temperature
    !-----------------------------------------------------------------------------
      xx = radus/cosl(j)
      do k = 1, lev
        do i = 1, nxj
          ut(i,k,jj)  = ut(i,k,jj)*xx
          vt(i,k,jj)  = vt(i,k,jj)*xx
          upp(i,k)    = up(i,k,jj)*xx
          vpp(i,k)    = vp(i,k,jj)*xx
          tt(i,k,jj)  = tt(i,k,jj)*pk(i,k,jj) / (1.0+0.608*qt(i,k,jj))
          tpp(i,k)  = ttp(i,k,jj)*pkp(i,k) / (1.0+0.608*qp(i,k,jj))
          ttpp(i,k) = ttp(i,k,jj) / (1.0+0.608*qp(i,k,jj))
        enddo
      enddo
      !-----------------------------------------------------------------------------
      ! save the old control values for SPPT
      !-----------------------------------------------------------------------------
      if (dosppt) then
        ! Save u, v, t, and q for SPPT
        do k=1,lev
          do i=1,nxj
            ut_save_sppt(i,k,jj)=ut(i,k,jj)
            vt_save_sppt(i,k,jj)=vt(i,k,jj)
            tt_save_sppt(i,k,jj)=tt(i,k,jj)
            qt_save_sppt(i,k,jj)=qt(i,k,jj)
          enddo
        enddo
      endif ! end dosppt if stetement

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     start physical process calculation (from long to short time scale)
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!-------------------------------------------------------------------------
!     calculate ozone concentrtion
!-------------------------------------------------------------------------
      if (doo3l) then
        if (ntoz .eq. 0) then
          if (doozon) then
            call rozone(nxjp(j),nxp,lev,plt(1,1,jj),o3l(1,1,jj),sinl(j),julian)
          endif
        else
          call rozphys(nxjp(j),nxp,lev,dta,iter,xlat(j),julian,o3l(1,1,jj),&
                      tt(1,1,jj),plt(1,1,jj),ps(1,jj),myrank)
        endif ! for ntoz
      endif ! for doo3l
!=======================================================================
!   Radiation scheme
!=======================================================================
      if (uprad .and. (irad .eq. 1))  then
         call radtn99 ( fluxcl,ozon,nxjp(j),nxp,lev,ncld,lvlwx(jj),julian       &
                    , stbo,s0,grav                                              &
                    , cp,ptrad,dsigma,sinl(j),cosz(1,jj),albedo2(1,jj),tg(1,jj) &
                    , curate(1,jj),pst(1,jj),plt(1,1,jj),tt(1,1,jj)             &
                    , qt(1,1,jj),o3l(1,1,jj)                                    &
                    , plcl(1,jj),cumtop(1,jj),ss(1,jj),rs(1,jj)                 &
                    , asol_clr(1,jj),olr_clr(1,jj),ss_clr(1,jj),rs_clr(1,jj)    &
                    , dtradn,asr(1,j),alr(1,j)                                  &
                    , asr_clr(1,j),alr_clr(1,j)                                 &
                    , xsr(1,j),xlr(1,j),acld(1,j),aflxd(1,j),aflxu(1,j)         &
                    , ilx(1,jj),ibx(1,jj),cofx(1,jj),sdpbl(1,jj),ctot(1,jj)     &
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
         do i = 1, nxj
         asol(i,jj) = plcl(i,jj)
         olr(i,jj)  = cumtop(i,jj)
!cc      tg2 = tg(i,jj)*tg(i,jj)
!cc      rld(i,jj)  = stbo*(tg2*tg2) - rs(i,jj)
         enddo
         endif  ! for uprad .and. irad=1
!--------------------------------------------------------------------------------
!   RRTMG scheme
!--------------------------------------------------------------------------------
      if (uprad .and. (irad .eq. 2))  then

        if ((isubc_lw .eq. 2) .or. (isubc_sw .eq. 2)) then
            ii=nxjstart(j)
            do i = 1 , nxj
               icsdsw(i) = ixseed(ii,j,1)
               icsdlw(i) = ixseed(ii,j,2)
               ii=ii+1
            enddo
        endif  !isubc_lw
        if ( nmgwor .eq. 1 ) then
          do i = 1, nxj
            rstd(i) = std(i,jj)
          enddo
        else if ( nmgwor .eq. 2 ) then
          do i = 1, nxj
            rstd(i) = hprime_b(i,1,jj)
          enddo
        endif
!    
!
        do k=1,lev
          kc=lev-k+1
          do i = 1, nxj
            dotc(i,kc)=vvel(i,k,jj)
          enddo
        enddo

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
           ( sigma,ps(1,jj),pltp,rstd,                                     &
             tpp,qp(1,1,jj),o3l(1,1,jj),dotc,tg(1,jj),                     &
             slimsk   ,cice(1,jj),xtice(1,jj),                             &
             snr(1,jj),sncover(1,jj),snoalb(1,jj),z0(1,jj),                &
             alvsf(1,jj),alnsf(1,jj),alvwf(1,jj),                          &
             alnwf(1,jj),facsf(1,jj),facwf(1,jj),                          &
             curate(1,jj),icsdsw,icsdlw,                                   &
             sinl(j),cosl(j),xlat(j),xlonr(1,jj),jdat,d2r,xkapa,           &
             ptrad,dtlw,dtsw,lsswr,lslwr,lssav,                            &
             nfxr,j,                                                       &
             nxp,nxjp(j),lev,ncld,lprnt,ipt,kdt,rsolhr,                    &
             uni_cloud,lmfshal,lmfdeep2,                                   &
             deltaq(1,1,jj),sup,cnvwr(1,1,jj),cnvcr(1,1,jj),               &
             ftp(1,1,jj),ftp1(1,1,jj),fqp(1,1,jj),fqp1(1,1,jj),nmmiph,     &
!  ---  outputs:
             asol(1,jj),olr(1,jj),ss(1,jj),rs(1,jj),                       &
             sld(1,jj),rld(1,jj),tsflw(1,jj),                              &
             ctot(1,jj),chig(1,jj),cmid(1,jj),clow(1,jj),                  &
             clds(1,1,jj),asl(1,1,jj),atl(1,1,jj),                         &
             fusl(1,1,jj),fdsl(1,1,jj),fuir(1,1,jj),fdir(1,1,jj),          &
             fuslr(1,1,jj),fdslr(1,1,jj),fuirr(1,1,jj),fdirr(1,1,jj),      &
             asl_clr(1,1,jj),atl_clr(1,1,jj),cosz(1,jj),                   &
             asol_clr(1,jj),olr_clr(1,jj),ss_clr(1,jj),rs_clr(1,jj),       &
             sld_clr(1,jj),rld_clr(1,jj),sfalb(1,jj),sfemis(1,jj))
          do k = 1, lev
            do i = 1, nxj
              dtrad(i,k,jj) = asl(i,k,jj) + atl(i,k,jj)
            enddo
          enddo
      endif  ! for uprad .and. irad=2

      if ( dorad ) then
        call dcyc2t3                                                  &
          !  ---  inputs:
          ( solhr,slag,sdec,cdec,sinl(j),cosl(j),                     &
            xlonr(1,jj),cosz(1,jj),tg(1,jj),tpp(1,lev),tsflw(1,jj),   &
            sld(1,jj),ss(1,jj),rld(1,jj),asl(1,1,jj),atl(1,1,jj),     &
            asl_clr(1,1,jj),atl_clr(1,1,jj),nxp, nxjp(j), lev,        &
!  ---  outputs:
            dtradn,dtradc(1,1,jj),sld_adj,ss_adj,rld_adj,             & 
            rs_adj,xmu(1,jj) )

        do i = 1, nxj
          rld_adj(i) = rld_adj(i) * sfemis(i,jj)
          rs_adj(i) = rs_adj(i) * sfemis(i,jj)
        enddo
      endif

!xb110> save the variables for TDK before doing PBL parameterization
      do k = 1,lev
       do i = 1,nxj
        u0(i,k) = ut(i,k,jj)
        v0(i,k) = vt(i,k,jj)
        t0(i,k) = tt(i,k,jj)
       end do
      end do

      do k = 1,lev*ncld
       do i = 1,nxj
        q0(i,k) = qt(i,k,jj)
       end do
      end do
!xb110<
!
!=======================================================================
! PBL scheme
!=======================================================================
      if ( dopbl .and. (nmpbl.eq.1 .and. nmland.eq.1))                        &
         call pbltke ( nxjp(j),nxp,lev,ktpbl,dta,grav,rgas,cp,xkapa,hltm,ptop &
                     , tice,hice,tg(1,jj),z0(1,jj),land(1,jj)                 &
                     , sgeo(1,jj),phi,pst(1,jj),upp,vpp                       &
                     , ttpp,qp(1,1,jj),ut(1,1,jj),vt(1,1,jj)                  &
                     , tt(1,1,jj),qt(1,1,jj),pk(1,1,jj),pk2(1,1,jj)           &
                     , ustar(1,jj),tstar(1,jj),qstar(1,jj),e(1,1,jj)          &
                     , eps(1,1,jj),hflux(1,jj),qflux(1,jj)                    &
                     , gwclim(1,jj),tgclim(1,jj),ocean(1,jj),ice(1,jj)        &
!                     , snr(1,jj),totalp(1,jj),ss_adj,rs(1,jj),albx(1,jj)      &
                     , snr(1,jj),totalp(1,jj),ss_adj,rs(1,jj),sfalb(1,jj)      &
                     , ipblmx(1,j),xkmx(1,j),ijdg(j),xkmd,itypbl              &
                     , t2(1,jj),rh2(1,jj),u10(1,jj),v10(1,jj)                 &
                     , rld_adj,stbo                                           &
                     , km_soil,smc(1,1,jj),stc(1,1,jj),canopy(1,jj)           &
                     , runoff(1,jj),sigmaf(1,jj),istyp(1,jj),ivegtyp(1,jj)    &
                     , wlt,ref,tsat,dfkt,xktk,dfk )

      if ( dopbl .and. (nmpbl.eq.2 .and. nmland.eq.1))                        &
         call pbltke_n ( nxjp(j),nxp,lev,ktpbl,dta,grav,rgas,cp,xkapa,hltm,ptop &
                     , tice,hice,tg(1,jj),z0(1,jj),land(1,jj)                 &
                     , sgeo(1,jj),phi,pst(1,jj),upp,vpp                       &
                     , ttpp,qp(1,1,jj),ut(1,1,jj),vt(1,1,jj)                  &
                     , tt(1,1,jj),qt(1,1,jj),pk(1,1,jj),pk2(1,1,jj)           &
                     , ustar(1,jj),tstar(1,jj),qstar(1,jj),e(1,1,jj)          &
                     , eps(1,1,jj),hflux(1,jj),qflux(1,jj)                    &
                     , gwclim(1,jj),tgclim(1,jj),ocean(1,jj),ice(1,jj)        &
!                     , snr(1,jj),totalp(1,jj),ss_adj,rs(1,jj),albx(1,jj)      &
                     , snr(1,jj),totalp(1,jj),ss_adj,rs(1,jj),sfalb(1,jj)      &
                     , ipblmx(1,j),xkmx(1,j),ijdg(j),xkmd,itypbl              &
                     , t2(1,jj),rh2(1,jj),u10(1,jj),v10(1,jj)                 &
                     , rld_adj,stbo                                         &
                     , km_soil,smc(1,1,jj),stc(1,1,jj),canopy(1,jj)           &
                     , runoff(1,jj),sigmaf(1,jj),istyp(1,jj),ivegtyp(1,jj)    &
                     , wlt,ref,tsat,dfkt,xktk,dfk,ncld,dsigma,j )
      if ( dopbl .and. (nmland.eq.2))                                         &
       call pbl_noah ( nxjp(j),nxp,lev,ktpbl,dta,grav,rgas,cp,xkapa,hltm,ptop &
                     , tice,hice,tg(1,jj),z0(1,jj),land(1,jj)                 &
                      , sgeo(1,jj),phi,pst(1,jj),upp,vpp                      &
                     , ttpp,qp(1,1,jj),ut(1,1,jj),vt(1,1,jj)                  &
                     , tt(1,1,jj),qt(1,1,jj),pk(1,1,jj),pk2(1,1,jj)           &
                     , ustar(1,jj),tstar(1,jj),qstar(1,jj),e(1,1,jj)          &
                     , eps(1,1,jj),hflux(1,jj),qflux(1,jj),itimestep          &
                     , gwclim(1,jj),tgclim(1,jj),ocean(1,jj),ice(1,jj)        &
!                     , snr(1,jj),totalp(1,jj),ss_adj,rs(1,jj),albx(1,jj)      &
                     , snr(1,jj),totalp(1,jj),ss_adj,rs(1,jj),sfalb(1,jj)      &
                     , ipblmx(1,j),xkmx(1,j),ijdg(j),xkmd,itypbl              &
                     , t2(1,jj),q2(1,jj),rh2(1,jj),rh10(1,jj),u10(1,jj)       &
                     , v10(1,jj),fm(1,jj),fh(1,jj),fm10(1,jj),fh2(1,jj)       &
                     , srflag(1,jj),rld_adj,stbo                            &
                     , km_soil,smc(1,1,jj),stc(1,1,jj),canopy(1,jj)           &
                     , runoff(1,jj),sigmaf(1,jj),istyp(1,jj),ivegtyp(1,jj)    &
                     , ncld,dsigma                                            &
                     , slopetyp(1,jj)                                         &
                     , slc(1,1,jj),sncover(1,jj),sndepth(1,jj)                &
                     , shdmax(1,jj),shdmin(1,jj),snoalb(1,jj),albedo2(1,jj)   &
                     , sld_adj,zice(1,jj),cice(1,jj),xtice(1,jj)            &
                     , hpbl(1,jj),asl(1,1,jj),atl(1,1,jj),xmu(1,jj),gfx(1,jj) &
                     , kpbl(1,jj),nmpbl,nmmiph,j,isot,ivegsrc,sfemis(1,jj)    &
                     , dudtc,dvdtc,dtdtc,dqdtc)


!
!
!     recompute phi by tt after pbl to ensure consistence of phi & phi2
!
      call get_phi(nxjp(j),nxp,lev,ptop,cp,rgas,grav,sgeo(1,jj),      &
                  pk(1,1,jj),pk2(1,1,jj),tt(1,1,jj),qt(1,1,jj),       &
                  phii,phi)
!

#ifdef VERBOSE
        ut_pbl(1:nxj,1:lev,jj)      = ut(1:nxj,1:lev,jj)
        vt_pbl(1:nxj,1:lev,jj)      = vt(1:nxj,1:lev,jj)
        tt_pbl(1:nxj,1:lev,jj)      = tt(1:nxj,1:lev,jj)
        qt_pbl(1:nxj,1:lev*ncld,jj) = qt(1:nxj,1:lev*ncld,jj)
#endif


!=======================================================================
! topograpic gravity wave drag
!=======================================================================
      if(dograv .and. (nmgwor .eq. 1) ) then 
        do k=1,lev
          kc=lev-k+1
          do i=1,nxj
            tt(i,k,jj) = tt(i,k,jj) + dtdtc(i,kc)*dta
            ut(i,k,jj) = ut(i,k,jj) + dudtc(i,kc)*dta
            vt(i,k,jj) = vt(i,k,jj) + dvdtc(i,kc)*dta
            qt(i,k,jj) = qt(i,k,jj) + dqdtc(i,kc)*dta
          enddo
        enddo

        call gwdp (j,nxjp(j),nxp,lev,                                   &
                  ut(1,1,jj),vt(1,1,jj),tt(1,1,jj),qt(1,1,jj),         &
                  plt(1,1,jj),pk(1,1,jj),pk2(1,1,jj),phi,std(1,jj),dta,&
                  grav,rgas,cp,drag(1,1,jj),ugws(1,jj),vgws(1,jj),cgw)
      endif

      if(dograv .and. (nmgwor .eq. 2) ) then
        lprnt = .false.
!        ipr = 1
!
!        cdmbgwd(1)       = 2.00      ! mtn blking and gwd tuning factors
!        cdmbgwd(2)       = 0.25      ! mtn blking and gwd tuning factors
        cdmbgwd(1)       = cmbk      ! mtn blocking tuning factors
        cdmbgwd(2)       = cgwd      ! gwd tuning factors

        do i = 1, nxj
          hprime(i)=hprime_b(i,1,jj)
          oc(i) = hprime_b(i,2,jj)
          theta(i)  = hprime_b(i,11,jj)
          gamma(i)  = hprime_b(i,12,jj)
          sigmaog(i)  = hprime_b(i,13,jj)
          elvmax(i) = hprime_b(i,14,jj)
        enddo

        do k = 1, 4
          do i = 1, nxj
            oa4(i,k) = hprime_b(i,k+2,jj)
            clx(i,k) = hprime_b(i,k+6,jj)
          enddo
        enddo

        do k=1,lev
          kc=lev-k+1
          do i=1,nxj
            prsl(i,kc) = 100.0*plt(i,k,jj) ! pa
            prslk(i,kc)=(plt(i,k,jj)/1000.)**xkapa
            del(i,kc) = 100.0*( dsigma(k,1)*pst(i,jj)+dsigma(k,2))  !  pa
            phil(i,kc) = phi(i,k)-sgeo(i,jj)
            qtc(i,kc) = qt(i,k,jj)
            ttc(i,kc) = tt(i,k,jj)
            utc(i,kc) = ut(i,k,jj)
            vtc(i,kc) = vt(i,k,jj)
          enddo
        enddo
!
!  for interface pressure
        do k = 1,lev+1
          kc= (lev+1)-k+1
          do i =1,nxj
            prsi(i,kc)   = 100.0*( sigma(k,1)*pst(i,jj)+sigma(k,2)+ptop ) !pa
          enddo
        enddo

        call gwdps(nxjp(j), nxp, nxp,  lev,                        &
               dvdtc, dudtc, dtdtc,utc, vtc, ttc,qtc,              &
               kpbl(1,jj),   prsi, del,   prsl, prslk,             &
               phii,  phil, dta,                                   &
               kdt,    hprime, oc, oa4, clx,                       &
!               theta,sigmaog,gamma,elvmax,dusfcg, dvsfcg,          &
               theta,sigmaog,gamma,elvmax,ugws(1,jj),vgws(1,jj),   &
               grav,cp,con_rd,con_rv, nx, mtnvar, cdmbgwd,         &
!              me,zmtnblck)
               me,zmtnblck,garea,hpbl(1,jj),tofd)

        do k=1,lev
          kc=lev-k+1
          do i=1,nxj
            tt(i,k,jj) = ttc(i,kc) + dtdtc(i,kc)*dta
            ut(i,k,jj) = utc(i,kc) + dudtc(i,kc)*dta
            vt(i,k,jj) = vtc(i,kc) + dvdtc(i,kc)*dta
            qt(i,k,jj) = qtc(i,kc) + dqdtc(i,kc)*dta
          enddo
        enddo
      endif  !(end of topo dograv and nmgwor=2)

!
!     update tt by radiation heating/cooling rate: dtrad (k/day)
!
      if ( doslavepp ) then
        do k = 1, lev
          do i = 1, nxj
            tt(i,k,jj) = tt(i,k,jj) + 0.5*dta*(dtrad(i,k,jj)+dtradn(i,k))/86400.0
            dtrad(i,k,jj) = dtradn(i,k)
          enddo
        enddo
      else
        do k = 1, lev
          do i = 1, nxj
            tt(i,k,jj) = tt(i,k,jj) + dta*dtradn(i,k)/86400.0
          enddo
        enddo
      endif
!
!     recompute phi by tt after pbl to ensure consistence of phi & phi2
!
      call get_phi(nxjp(j),nxp,lev,ptop,cp,rgas,grav,sgeo(1,jj),      &
                  pk(1,1,jj),pk2(1,1,jj),tt(1,1,jj),qt(1,1,jj),       &
                  phii,phi)
!
      do k=1,lev
        do i=1,nxj
          tt_bfcnv(i,k) = tt(i,k,jj)
        enddo
      enddo
!=======================================================================
!  cumulus scheme
!=======================================================================
      ! save old array for Thompson
      if ( nmmiph .eq. 18 ) then
        do k=1,lev
          do i = 1, nxj
            ice00(i,k) = qt(i,(ntiw-1)*lev+k,jj)
          enddo
        enddo
      endif

      if ( docup .and. (nmcup.eq. 1) )                               &
        call cupcwb (j,nxjp(j),nxp,my,lev,ktcup,dta,grav,rgas,cp,hltm,etop,prevap  &
                 , sgeo(1,jj),pst(1,jj),plt(1,1,jj),pk(1,1,jj),pk2(1,1,jj)   &
                 , tt(1,1,jj),qt(1,1,jj),phi,plcl(1,jj),cumtop(1,jj)         &
                 , rcup(1,jj),ncup(j),ptop,dsigma,icupmx(j),dtcupx(j)        &
                 , dtcupz(1,j),dqcupz(1,j),ijdg(j),dtcupd,dqcupd             &
                 , nlcl(1,j),nnegl(1,j),nosat(1,j),nwork(1,j)                &
                 , ntcup(1,j),nflx(1,j) )

    !cyea---->
    !c 20120926 for Tiedtke cumulus
      if ( docup .and. (nmcup .eq. 4 .and. ncld .ge. 2) ) then
        do k=1,lev
          do i = 1, nxj
            dotc(i,k)=vvel(i,k,jj)
          enddo
        enddo
        call cumastr_driv(nxjp(j),nxp,lev,dt,grav,rgas,cp,hltm,ptop &
                       , land(1,jj),sgeo(1,jj),phi,upp              &
                       , vpp,ttpp,qp(1,1,jj)                        &
                       , ut(1,1,jj),vt(1,1,jj),tt(1,1,jj)           &
                       , qt(1,1,jj),rcup(1,jj),pk(1,1,jj)           &
                       , pk2(1,1,jj),dotc,qflux(1,jj)               &
                       , kbot(1,jj),ktop(1,jj),ncld,sigma           &
                       , plt(1,1,jj),pst(1,jj),j,kuo(1,jj) )

    ! for rad input of convection cloud information
    ! bottom(plcl) layer and top(cumtop) layer in pressure(mb)
        do i=1,nxj
          if((kbot(i,jj).eq.lev-1) .and. (ktop(i,jj).eq.lev-1))then
            plcl(i,jj)=0.
            cumtop(i,jj)=0.
          else
            plcl(i,jj)=plt(i,kbot(i,jj),jj)
            cumtop(i,jj)=plt(i,ktop(i,jj),jj)
          endif
        enddo

        do i=1,nxj
          rcup(i,jj) = rcup(i,jj) * 1000.         ! mm/call
          kbot(i,jj) = lev - kbot(i,jj) + 1
          ktop(i,jj) = lev - ktop(i,jj) + 1
        enddo
      endif    !(end if nmcup=4)

      if ( docup .and. (nmcup .eq. 5 .and. ncld .ge. 2) ) then
        do k=1,lev
          do i = 1, nxj
            dotc(i,k)=vvel(i,k,jj)
          enddo
        enddo

        call cumastr_driv_n                                               &
               (nxjp(j)    ,nxp       ,lev        ,dta         ,grav     ,&
                rgas       ,cp        ,hltm       ,ptop      ,land(1,jj) ,&
                sgeo(1,jj) ,phi       ,u0         ,v0        ,t0         ,&
                q0         ,ut(1,1,jj),vt(1,1,jj) ,tt(1,1,jj),qt(1,1,jj) ,&
                rcup(1,jj) ,pk(1,1,jj),pk2(1,1,jj),dotc      ,qflux(1,jj),&
                kbot(1,jj) ,ktop(1,jj),ncld      ,sigma      ,            &
                plt(1,1,jj),pst(1,jj) ,j          ,islimsk   ,hflux(1,jj),&
                garea      ,kuo(1,jj) ,flash(1,jj))

        do i=1,nxj
         if((kbot(i,jj).eq.-1) .and. (ktop(i,jj).eq.-1))then
          plcl(i,jj)=0.
          cumtop(i,jj)=0.
         else
          plcl(i,jj)=plt(i,kbot(i,jj),jj)
          cumtop(i,jj)=plt(i,ktop(i,jj),jj)
         endif
        enddo

        do i=1,nxj
          rcup(i,jj) = rcup(i,jj) * 1000.         ! mm/call
          kbot(i,jj) = lev - kbot(i,jj) + 1
          ktop(i,jj) = lev - ktop(i,jj) + 1
        enddo
      endif    !(end if nmcup=5)

      if ( docup .and. (nmcup.eq.2 .or. nmcup.eq.3 .or. nmcup.eq. 6) ) then
        ! setting for nmcup=2,3 and new shallow convection
        ! setting for nmcup=6 and scale-aware shallow convection
        ! original :        
        !    call random_number(XKT2)
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
        enddo
        do k=1,lev
          kc=lev-k+1
          do i = 1, nxj
            dotc(i,kc)=vvel(i,k,jj)*0.1
          enddo
        enddo
        do k=1,lev
          kc=lev-k+1
            sl(kc)    = (sigma(k,1)+sigma(k+1,1))*0.5
          do i=1,nxj
            prsl(i,kc) = plt(i,k,jj)*0.1 ! change to cb
            del(i,kc)  = psfc(i)*dsigma(k,1)+dsigma(k,2)*0.1  !unit cb
            phil(i,kc)= phi(i,k)-sgeo(i,jj)
            qtc(i,kc) = qt(i,k,jj)
            qtr(i,kc) = qt(i,lev+k,jj)
            if ( nmmiph .gt. 2 ) qti(i,kc) = qt(i,(ntiw-1)*lev+k,jj)
            ttc(i,kc) = tt(i,k,jj)
            utc(i,kc) = ut(i,k,jj)
            vtc(i,kc) = vt(i,k,jj)
          enddo
        enddo

        do k = 1,lev
          do i = 1,nxj
            ztenh(i,k) = tt(i,k,jj)
            zqenh(i,k) = qt(i,k,jj)
            rho(i,k) = plt(i,k,jj)*100./ (con_rd*tt(i,k,jj))
          end do
        end do

        ! old version SAS
        if( nmcup .eq. 2)                                 &
          call sascnv(nxjp(j),nxp,lev,jcap,dta,del,sl,psfc,prsl,phil,qtr &
            ,qtc,ttc,utc,vtc,dotc,cldwrk(1,jj) &
            ,rcup(1,jj),kbot(1,jj),ktop(1,jj)                     &
            ,kuo(1,jj),slimsk,xkt2,ncld,dtcupx(j),icupmx(j)      &
            ,grav,cp,hltm,rgas,tice)

        ! new version SAS
        if( nmcup .eq. 3)                                 &
          call sascnv_n(nxjp(j),nxp,lev,jcap,dta,del,psfc,prsl,phil,qtr &
            ,qtc,ttc,utc,vtc,dotc,cldwrk(1,jj) &
            ,rcup(1,jj),kbot(1,jj),ktop(1,jj)                     &
            ,kuo(1,jj),slimsk,ncld                                &
          , grav,cp,hltm,rgas,tice)

        ! scale-aware SAS
        if( nmcup .eq. 6)                                 &
          call samfdeepcnv(nxjp(j),nxp,lev,dta,del,prsl,psfc,phil       &
            ,qtr,qti,qtc,ttc,utc,vtc,cldwrk(1,jj),rcup(1,jj),kbot(1,jj)  &
            ,ktop(1,jj),kuo(1,jj),islimsk,garea,dotc,ncld,cnvw,cnvc      &
            ,snow_flxn,ptun,pqun)

        ! for rad input of convection cloud information
        ! bottom(plcl) layer and top(cumtop) layer in pressure(mb)
        do i=1,nxj
          if((kbot(i,jj).eq.lev+1) .and. (ktop(i,jj).eq.0))then
            plcl(i,jj)=0.
            cumtop(i,jj)=0.
          else
            plcl(i,jj)=prsl(i,kbot(i,jj))*10.         !to mb
            cumtop(i,jj)=prsl(i,ktop(i,jj))*10.       !to mb
          endif
        enddo

        do i=1,nxj
          rcup(i,jj) = rcup(i,jj) * 1000.         ! mm/call
        enddo

        do k=1,lev
          kc=lev-k+1
          do i=1,nxj
            qt(i,k    ,jj) = max(qtc(i,kc),qmin)
            qt(i,k+lev,jj) = max(qtr(i,kc),qmin)
            if ( nmmiph .gt. 2 ) qt(i,(ntiw-1)*lev+k,jj) = qti(i,kc)
            dttmp = ttc(i,kc)-tt(i,k,jj)
            dutmp = utc(i,kc)-ut(i,k,jj)
            dvtmp = vtc(i,kc)-vt(i,k,jj)
            if ( kdt .eq. 1 .or. .not. doslavepp ) then
              tt(i,k,jj) = ttc(i,kc)
              ut(i,k,jj) = utc(i,kc)
              vt(i,k,jj) = vtc(i,kc)
            else
              tt(i,k,jj) = 0.5*( dttmp + dtcup(i,k,jj) ) + tt(i,k,jj)
              ut(i,k,jj) = 0.5*( dutmp + ducup(i,k,jj) ) + ut(i,k,jj)
              vt(i,k,jj) = 0.5*( dvtmp + dvcup(i,k,jj) ) + vt(i,k,jj)
            endif
            dtcup(i,k,jj)  = dttmp
            ducup(i,k,jj)  = dutmp
            dvcup(i,k,jj)  = dvtmp
            cnvwr(i,kc,jj) = cnvw(i,kc)
            cnvcr(i,kc,jj) = cnvc(i,kc)
            cnvw(i,kc)     = 0.
            cnvc(i,kc)     = 0.

            snow_flx(i,k) = snow_flxn(i,kc)
            ptu(i,k)      = ptun(i,kc)
            pqu(i,k)      = pqun(i,kc)
            cnvwn(i,k)    = cnvwr(i,kc,jj)
          enddo
        enddo

        ! for dianostic lightning calculation    
        do i =1,nxj
          kbotc(i,jj)   =lev-kbot(i,jj)+1
          ktopc(i,jj)   =lev-ktop(i,jj)+1
        enddo
        call lightning_ec (nxjp(j),nxp,lev,ptu,pqu,ztenh,zqenh     &
                        ,cnvwn,rho,kbotc(1,jj),ktopc(1,jj)     &
                        ,phi,plt(1,1,jj),snow_flx               &
                        ,flash(1,jj),islimsk,kuo(1,jj))
      endif  !(end of docup .or. (nmcup .eq. 2 .or. nmcup .eq. 3 .or. nmcup .eq. 6))

#ifdef VERBOSE
        ut_cmls(1:nxj,1:lev,jj)      = ut(1:nxj,1:lev,jj)
        vt_cmls(1:nxj,1:lev,jj)      = vt(1:nxj,1:lev,jj)
        tt_cmls(1:nxj,1:lev,jj)      = tt(1:nxj,1:lev,jj)
        qt_cmls(1:nxj,1:lev*ncld,jj) = qt(1:nxj,1:lev*ncld,jj)
#endif
!=======================================================================
! convective gravity wave drag
!=======================================================================
      if( docgrav .and. upnor .and. (nmgwcv .eq. 1) )then
        call nor_gwdp (j,nxjp(j),nxp,lev,                         &
                  ut(1,1,jj),vt(1,1,jj),tt(1,1,jj),qt(1,1,jj),    &
                  plt(1,1,jj),pk(1,1,jj),pk2(1,1,jj),phi,dta,&
                  grav,rgas,sinl(j),cosl(j),drag_u,drag_v,cp,ptop)
      endif ! (end of docgrav .and. nmgwcv.eq.1)

      if( docgrav .and. (nmgwcv.eq.2) )then
        cgwf(1)  = 0.5      ! cloud top fraction for convective gwd scheme
        cgwf(2)  = 0.05     ! cloud top fraction for convective gwd scheme

        do  k = 1,lev+1
          do  i = 1, nxj
            prsi(i,k) = 100.0*( sigma(k,1)*pst(i,jj)+sigma(k,2) + ptop) !pa
          enddo
        enddo

        do k=1,lev
          do i=1,nxj
            prsl(i,k) = 100.0*plt(i,k,jj) ! pa
            del(i,k) = 100.0*( dsigma(k,1)*pst(i,jj)+dsigma(k,2))  !pa
          enddo
        enddo

        do i = 1, nxj
          cumabs(i) = 0.0
          work3(i)  = 0.0
        enddo

        do k = 1, lev
          do i = 1, nxj
            if ((k <= lev-kbot(i,jj)+1) .and. (k >= lev-ktop(i,jj)+1)) then
              cumabs(i) = cumabs(i) + (tt(i,k,jj) - tt_bfcnv(i,k)) * del(i,k)
              work3(i)  = work3(i)  + del(i,k)
            endif
          enddo
        enddo

        do i=1,nxj
          if (work3(i) > 0.0) cumabs(i) = cumabs(i) / (dta*work3(i))
      !   dlength(i)  = sqrt( tem1*tem1+tem2*tem2 )
      !   cldf(i)     = cgwf(1)*work1(i) + cgwf(2)*work2(i)
        enddo

        dlength(1:nxj)  = sqrt( tem1*tem1+tem2*tem2 )
        cldf(1:nxj)     = cgwf(1)*work1(1:nxj) + cgwf(2)*work2(1:nxj)

        call gwdc (nxjp(j),nxp,nxp,lev,u0,v0,                        &
                     t0,q0,prsl,prsi,del,                            &
                     ktop(1,jj),kbot(1,jj),kuo(1,jj),cldf,cumabs,    &
                     grav,cp,con_rd,con_fvirt,dta,dlength,           &
                     utgwc,vtgwc,tauctx,taucty,j)
        do k=1,lev
          do i=1,nxj
            eng0 = 0.5*(ut(i,k,jj)*ut(i,k,jj)+vt(i,k,jj)*vt(i,k,jj))
            ut(i,k,jj) = ut(i,k,jj) + utgwc(i,k) * dta
            vt(i,k,jj) = vt(i,k,jj) + vtgwc(i,k) * dta
            eng1 = 0.5*(ut(i,k,jj)*ut(i,k,jj)+vt(i,k,jj)*vt(i,k,jj))
            tt(i,k,jj) = tt(i,k,jj) + (eng0-eng1)/(dta*cp)
          enddo
        enddo
      endif  !(end of docgrav and nmgwcv=2)

!=======================================================================
! shallow convection 
!=======================================================================
      if( doshl .and. (nmshl.eq.2 .or. nmshl.eq.3) ) then
        do i=1,nxj
          psfc(i)  = pst(i,jj)*0.1        ! change to cb
        enddo
      ! psfc(1:nxj)  = pst(1:nxj,jj)*0.1 ! change to cb

        do k=1,lev
          kc=lev-k+1
          do i = 1, nxj
            dotc(i,kc)=vvel(i,k,jj)*0.1
          enddo
        enddo

        do k=1,lev
          kc=lev-k+1
          sl(kc)    = (sigma(k,1)+sigma(k+1,1))*0.5
          do i=1,nxj
            prsl(i,kc) = plt(i,k,jj)*0.1 ! change to cb
            del(i,kc)  = psfc(i)*dsigma(k,1)+dsigma(k,2)*0.1  !unit cb
            phil(i,kc) = phi(i,k)-sgeo(i,jj)
            qtc(i,kc)  = qt(i,k,jj)
            qtr(i,kc)  = qt(i,lev+k,jj)
            if ( nmmiph .gt. 2 ) qti(i,kc) = qt(i,(ntiw-1)*lev+k,jj)
            ttc(i,kc)  = tt(i,k,jj)
            utc(i,kc)  = ut(i,k,jj)
            vtc(i,kc)  = vt(i,k,jj)
          enddo
        enddo
!
        do i=1,nxj
           heat(i)=-ustar(i,jj)*tstar(i,jj)
           evap(i)=-ustar(i,jj)*qstar(i,jj)
        enddo
      ! heat(1:nxj)=-ustar(1:nxj,jj)*tstar(1:nxj,jj)
      ! evap(1:nxj)=-ustar(1:nxj,jj)*qstar(1:nxj,jj)
        
!
        ! new version shalcon
        if( nmshl.eq.2 ) then
          call shalcnv_new(nxjp(j),nxp,lev,jcap,dta,del,prsl,psfc,phil,qtr &
            ,qtc,ttc,utc,vtc                         &
            ,rcup2,kbot(1,jj),ktop(1,jj)                                &
            ,kuo(1,jj),slimsk,dotc,ncld,hpbl(1,jj),heat,evap &
            ,grav,cp,hltm,rgas,tice)
        endif

        ! scale-aware shalcon
        if( nmshl.eq.3 ) then
          call samfshalcnv(nxjp(j),nxp,lev,dta,del,prsl,psfc,phil,qtr   &
            ,qti,qtc,ttc,utc,vtc,rcup2,kbot(1,jj),ktop(1,jj),kuo(1,jj)  &
            ,islimsk,garea,dotc,ncld,hpbl(1,jj),cnvw,cnvc)
        endif

        do i=1,nxj
          rcup(i,jj) = rcup(i,jj)+rcup2(i) * 1000.         ! mm/call
        enddo
      !rcup(1:nxj,jj) = rcup(1:nxj,jj)+rcup2(1:nxj) * 1000.         ! mm/call

        do k=1,lev
          kc=lev-k+1
          do i=1,nxj
            qt(i,k    ,jj) = max(qtc(i,kc),qmin)
            qt(i,k+lev,jj) = max(qtr(i,kc),qmin)
            if ( nmmiph .gt. 2 ) qt(i,(ntiw-1)*lev+k,jj) = qti(i,kc)
            dttmp = ttc(i,kc)-tt(i,k,jj)
            dutmp = utc(i,kc)-ut(i,k,jj)
            dvtmp = vtc(i,kc)-vt(i,k,jj)
            if ( kdt .eq. 1 .or. .not. doslavepp ) then
              tt(i,k,jj) = ttc(i,kc)
              ut(i,k,jj) = utc(i,kc)
              vt(i,k,jj) = vtc(i,kc)
            else
              tt(i,k,jj) = 0.5*( dttmp + dtshl(i,k,jj) ) + tt(i,k,jj)
              ut(i,k,jj) = 0.5*( dutmp + dushl(i,k,jj) ) + ut(i,k,jj)
              vt(i,k,jj) = 0.5*( dvtmp + dvshl(i,k,jj) ) + vt(i,k,jj)
            endif
            dtshl(i,k,jj)  = dttmp
            dushl(i,k,jj)  = dutmp
            dvshl(i,k,jj)  = dvtmp
            cnvwr(i,kc,jj) = cnvwr(i,kc,jj) + cnvw(i,kc)
            cnvcr(i,kc,jj) = cnvcr(i,kc,jj) + cnvc(i,kc)
          enddo
        enddo
      endif  !(end of doshl .and. (nmshl.eq.2 .or. nmshl.eq.3)

      if ( doshl .and. (nmshl .eq.1))                                          &
         call shlcon ( nxjp(j),nxp,lev,ktshl,dta,grav,rgas,cp,xkapa,hltm,ptop  &
                     , dsigma,tg(1,jj),pk(1,1,jj),pst(1,jj),sgeo(1,jj),phi     &
                     , plt(1,1,jj),tt(1,1,jj), qt(1,1,jj),nshl(j)              &
                     , rcup(1,jj),ncld )

      ! ice number concentration modification for Thompson
      if ( nmmiph .eq. 18 ) then
        icem = 4./3.*con_pi*3.2768*1.e-14*890.
        do k=1,lev
          do i=1,nxj
            qni = (qt(i,(ntiw-1)*lev+k,jj)-ice00(i,k))/icem
            if ( qni .gt. 0. ) then
              qt(i,(ntinc-1)*lev+k,jj) = qt(i,(ntinc-1)*lev+k,jj) + qni
            endif
          enddo
        enddo
      endif

!=======================================================================
! cloud microphysics
!=======================================================================
      if ( dolsp .and. (ntcw.eq.0)) then
         call lsp ( tt(1,1,jj),qt(1,1,jj),plt(1,1,jj),pst(1,jj),dsigma         &
                  , grav,nxjp(j),nxp,lev,evaprh,rlsp(1,jj),cp,hltm,nlsp(1,j)   &
                  , ilsp(1,j) )
      endif

      if ( dolsp .and. (nmmiph.eq.2)) then
        deg_ju=23.45*sin(d2r*(360./365.)*(julian+284.))
        arg=xlat(j)-deg_ju
        if(arg.gt.90.)then
          arg=89.9999
        else if(arg.lt.-90.)then
          arg=-89.9999
        endif

        lprnt=.false.
        psfc(:)  = pst(:,jj)*0.1        ! change to cb
        psautco(:)  = 4.0e-4
!            psautco(i)  = 8.0e-4 * work1(i) + 5.0e-4 * work2(i)

        do k=1,lev
          kc=lev-k+1
          do i=1,nxj
!!           tem   = (rhztop-rhzbot) / (pk(i,k,jj)-pk(i,lev,jj))
!!           tmprhc = rhzbot + tem * (pk(i,k,jj)-pk(i,lev,jj))
!!           rhc(i,kc) = 0.999 * work1(i) + tmprhc * work2(i)
!            rhc(i,kc) = 0.999 * work1(i) + 0.85 * work2(i)
!
!!!            rhc(i,kc)=0.999-0.08*cos(d2r*arg)**2    !a3
            rhc(i,kc)=0.98-0.07*cos(d2r*xlat(j))**2.0    !v3
!            rhc(i,kc)=0.98-0.07*cos(d2r*arg)**2.0    !wsm6
!            tem   = (max(min(plt(i,k,jj),900.)-700.,0.01) / 200.)
!            rhc(i,kc)=tem*rhc(i,kc)+(1.-tem)*0.7
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
          enddo
        enddo
!
        do k=1,lev
          kc=lev-k+1
          do i=1,nxj
            prsl(i,kc) = plt(i,k,jj)*0.1 ! change to cb
!            phil(i,kc) = phi(i,k)-sgeo(i,jj)
            del(i,kc) = (dsigma(k,1)*pst(i,jj)+dsigma(k,2))*0.1  ! change to cb
            qtc(i,kc) = qt(i,k,jj)
            qtr(i,kc) = qt(i,lev+k,jj)
            ttc(i,kc) = tt(i,k,jj)
          enddo
        enddo

        if ( pdfcloud ) then
          call gscondp(nxjp(j),nxp,lev,dta,prsl,psfc,  &
                      qtc,qtr,ttc,           &
                      ftp (1,1,jj),fqp (1,1,jj),fpsp (1,jj),&
                      ftp1(1,1,jj),fqp1(1,1,jj),fpsp1(1,jj),&
                      rhc,deltaq(1,1,jj),sup,lprnt,kdt)
!
          call precpdp(nxjp(j),nxp,lev,dta,del,prsl,psfc, &
                      qtc, qtr, ttc,           &
                      rlsp(1,jj),rhc,deltaq(1,1,jj),psautco,lprnt)
        else
          call gscond(nxjp(j),nxp,lev,dta,prsl,psfc,  &
                      qtc,qtr,ttc,           &
                      ftp (1,1,jj),fqp (1,1,jj),fpsp (1,jj),&
                      ftp1(1,1,jj),fqp1(1,1,jj),fpsp1(1,jj),&
                      rhc,lprnt)
!
          call precpd(nxjp(j),nxp,lev,dta,del,prsl,psfc, &
                      qtc, qtr, ttc,           &
                      rlsp(1,jj),rhc,psautco,lprnt)
!xb110>
! precipitation over mid-latitude perform not very well, especially
! in climatology.
!!        call precpd_n(nxjp(j),nxp,lev,dta,del,prsl,psfc,              &
!!                    qtc, qtr, rm, sm, phil, ttc,                      &
!!                    rlsp(1,jj), slsp(1,jj), rainp, rhc, lprnt)
!
!xb110<
        endif
       !do i=1,nxj
       !  rlsp(i,jj) = rlsp(i,jj) * 1000.         ! mm/call
       !enddo
        rlsp(1:nxj,jj) = rlsp(1:nxj,jj) * 1000.         ! mm/call

        do k=1,lev
          kc=lev-k+1
          do i=1,nxj
            qt(i,k    ,jj) = max(qtc(i,kc),qmin)
            qt(i,k+lev,jj) = max(qtr(i,kc),qmin)
            dttmp = ttc(i,kc)-tt(i,k,jj)
            if ( kdt .eq. 1 .or. .not. doslavepp ) then
              tt(i,k,jj) = ttc(i,kc)
            else
              tt(i,k,jj) = 0.5*( dttmp + dtlsp(i,k,jj) )+tt(i,k,jj)
            endif
            dtlsp(i,k,jj)  = dttmp
          enddo
        enddo
      endif !( dolsp .and. nmmiph.eq.2 )
!
      if ( dolsp .and. (nmmiph.eq.6 .or.       & ! WSM6
           nmmiph.eq.8 .or. nmmiph.eq.18 .or.  & ! Thompson
           nmmiph.eq.11 .or. nmmiph.eq.12 .or. nmmiph.eq.13 .or. & !GFDL MP
           nmmiph.eq.15 .or. nmmiph.eq.16) ) then !Goddard MP

! for microphysics
      area = tem1*tem2  !area of grid box (m^2)

! define rhc for GFDL MP
      rhc_mp = 1.
#ifdef rhc_GFDL
      if ( arg .gt. 45. ) then
        arg = 45.
      elseif ( arg .lt. -45 ) then
        arg = -45.
      endif
      do k=1,lev
        do i=1,nxj
!          rhc_mp(i,k) = 1.0 - 0.02*cos(d2r*arg)**2
          if ( plt(i,k,jj)/plt(i,lev,jj) .lt. 0.4 ) then
             tem = ( plt(i,k,jj)/plt(i,lev,jj) )**0.015
             if ( tem .le. 0.95 ) tem = 0.95
             rhc_mp(i,k) = rhc_mp(i,k)*tem
          endif
        enddo
      enddo
#endif

      call mp_scheme                                                   &
!  ---  inputs:
           ( nmmiph,nxp,nxjp(j),lev,ncld,plt(1,1,jj),ptop,             &
             dsigma,phii,islimsk,q0,kdt,tpi,me,dta,area,jj,            &
             itimestep,sgeo(1,jj),phi,rhc_mp,pk(1,1,jj),               &
             snr(1,jj),                                                &
!  ---  inputs/outputs:
             tt(1,1,jj),qt(1,1,jj),clds(1,1,jj),                       &
             ut(1,1,jj),vt(1,1,jj),vvel(1,1,jj),                       &
             pst(1,jj),                                                &
!  ---  outputs:
             ftp(1,1,jj),ftp1(1,1,jj),fqp(1,1,jj),fqp1(1,1,jj),        &
             rlsp(1,jj),sr(1,jj) )
!
#ifdef update_dp
      if ( nmmiph.eq.12 .or. nmmiph.eq.13 ) then
        ! compute new time step pk, pk2, and plt
        call prexp_hybrid_cwb ( nxjp(j),nxp,lev,ptop,sigma,pst(1,jj), &
                            pk(1,1,jj),pk2(1,1,jj),plt(1,1,jj) )
      endif
#endif
      endif

!
      if ( dodry ) then
         do k=1,lev
           dsigpp(k) = dsigma(k,1)+dsigma(k,2)/1000.
         enddo
         call drychk ( j,tt(1,1,jj),pk(1,1,jj),dsigpp,nxjp(j),nxp,lev,ndry(j) )
      endif

    ! do  i = 1, nxj
    !   curate(i,jj) = rcup(i,jj) * 86400.0/dta
    ! enddo
      curate(1:nxj,jj) = rcup(1:nxj,jj) * 86400.0/dta

!=======================================================================
! SIT scheme, cal. dTsit/dt
!=======================================================================
      if(do_sit) then
        dtsit=dt

        istep= int(tau/(dt/3600.)+0.01)
        tauleft=float(int((tau-int(tau)+0.001)*3600./dt))*dt
        icurrenttau=int(tau)
        if(tauleft == 3600.0 )then
           icurrenttau=icurrenttau+1
           tauleft=0.0
        endif
        call dtgfix12(idtg,idtg_sitvdiff,icurrenttau)
        call time_weights(idtg_sitvdiff,tauleft)
        write(cdtg,'(i12)')idtg_sitvdiff
        read(cdtg,'(i4,i2,i2,i2,i2)')yr,mo,dy,hr,mn
!        ssttau = mod(tau+0.001, 24.)
        tauhr=float(hr)+tauleft/3600.

        
        allocate(sstm(nxp))
        allocate(rstm(nxp))
        allocate(hfluxtm(nxp))
        allocate(qfluxtm(nxp))
        allocate(u10tm(nxp))
        allocate(v10tm(nxp))
        allocate(rlsptm(nxp))
        allocate(rcuptm(nxp))
        allocate(ustartm(nxp))
        allocate(t2tm(nxp))
        allocate(rh2tm(nxp))
        allocate(psttm(nxp))
        allocate(cicetm(nxp))
        allocate(snrtm(nxp))
        allocate(zicetm(nxp))
        allocate(xticetm(nxp))
        allocate(obswtbtm(nxp))
        allocate(tgtm(nxp))

        if (jj .eq. 1) then
          dtfsit=dtfsit+dtsit
        endif
        do ii = 1, nxj
          i=nxjstart(j)+ii-1
          sitlat(ii)=xlat(j)
          IF(xlon(i,jj).LT. 0.) THEN
            sitlon(ii,jj)=xlon(i,jj)+360.
          ELSE
            sitlon(ii,jj)=xlon(i,jj)
          ENDIF

          if (.NOT. lrun_sitvdiff .or. ic_sit .eq. 1) then
            ssfsit(ii,jj)=ssfsit(ii,jj)+ss_adj(ii)*dtsit
            rsfsit(ii,jj)=rsfsit(ii,jj)+rs_adj(ii)*dtsit
            hfluxfsit(ii,jj)=hfluxfsit(ii,jj)+hflux(ii,jj)*dtsit
            qfluxfsit(ii,jj)=qfluxfsit(ii,jj)+qflux(ii,jj)*dtsit
            u10fsit(ii,jj)=u10fsit(ii,jj)+u10(ii,jj)*dtsit
            v10fsit(ii,jj)=v10fsit(ii,jj)+v10(ii,jj)*dtsit
            rlspfsit(ii,jj)=rlspfsit(ii,jj)+rlsp(ii,jj)*dtsit
            rcupfsit(ii,jj)=rcupfsit(ii,jj)+rcup(ii,jj)*dtsit
            ustarfsit(ii,jj)=ustarfsit(ii,jj)+ustar(ii,jj)*dtsit
            t2fsit(ii,jj)=t2fsit(ii,jj)+t2(ii,jj)*dtsit
            rh2fsit(ii,jj)=rh2fsit(ii,jj)+rh2(ii,jj)*dtsit
            pstfsit(ii,jj)=pstfsit(ii,jj)+pst(ii,jj)*dtsit
            cicefsit(ii,jj)=cicefsit(ii,jj)+cice(ii,jj)*dtsit
            snrfsit(ii,jj)=snrfsit(ii,jj)+snr(ii,jj)*dtsit
            zicefsit(ii,jj)=zicefsit(ii,jj)+zice(ii,jj)*dtsit
            xticefsit(ii,jj)=xticefsit(ii,jj)+xtice(ii,jj)*dtsit
            obswtbfsit(ii,jj)=obswtbfsit(ii,jj)+obswtbnow(ii,jj)*dtsit
            tgfsit(ii,jj)=tgfsit(ii,jj)+tg(ii,jj)*dtsit
          endif

          if(lrun_sitvdiff)then
            if(ic_sit .eq. 1) then
              sstm(ii)=ssfsit(ii,jj)/dtfsit
              rstm(ii)=rsfsit(ii,jj)/dtfsit
              hfluxtm(ii)=hfluxfsit(ii,jj)/dtfsit
              qfluxtm(ii)=qfluxfsit(ii,jj)/dtfsit
              u10tm(ii)=u10fsit(ii,jj)/dtfsit
              v10tm(ii)=v10fsit(ii,jj)/dtfsit
              rlsptm(ii)=rlspfsit(ii,jj)/dtfsit
              rcuptm(ii)=rcupfsit(ii,jj)/dtfsit
              ustartm(ii)=ustarfsit(ii,jj)/dtfsit
              t2tm(ii)=t2fsit(ii,jj)/dtfsit
              rh2tm(ii)=rh2fsit(ii,jj)/dtfsit
              psttm(ii)=pstfsit(ii,jj)/dtfsit
              cicetm(ii)=cicefsit(ii,jj)/dtfsit
              snrtm(ii)=snrfsit(ii,jj)/dtfsit
              zicetm(ii)=zicefsit(ii,jj)/dtfsit
              xticetm(ii)=xticefsit(ii,jj)/dtfsit
              obswtbtm(ii)=obswtbfsit(ii,jj)/dtfsit
              tgtm(ii)=tgfsit(ii,jj)/dtfsit


              ssfsit(ii,jj)=0.
              rsfsit(ii,jj)=0.
              hfluxfsit(ii,jj)=0.
              qfluxfsit(ii,jj)=0.
              u10fsit(ii,jj)=0.
              v10fsit(ii,jj)=0.
              rlspfsit(ii,jj)=0.
              rcupfsit(ii,jj)=0.
              ustarfsit(ii,jj)=0.
              t2fsit(ii,jj)=0.
              rh2fsit(ii,jj)=0.
              pstfsit(ii,jj)=0.
              cicefsit(ii,jj)=0.
              snrfsit(ii,jj)=0.
              zicefsit(ii,jj)=0.
              xticefsit(ii,jj)=0.
              obswtbfsit(ii,jj)=0.
              tgfsit(ii,jj)=0.
            else
              sstm(ii)=ss_adj(ii)
              rstm(ii)=rs_adj(ii)
              hfluxtm(ii)=hflux(ii,jj)
              qfluxtm(ii)=qflux(ii,jj)
              u10tm(ii)=u10(ii,jj)
              v10tm(ii)=v10(ii,jj)
              rlsptm(ii)=rlsp(ii,jj)
              rcuptm(ii)=rcup(ii,jj)
              ustartm(ii)=ustar(ii,jj)
              t2tm(ii)=t2(ii,jj)
              rh2tm(ii)=rh2(ii,jj)
              psttm(ii)=pst(ii,jj)
              cicetm(ii)=cice(ii,jj)
              snrtm(ii)=snr(ii,jj)
              zicetm(ii)=zice(ii,jj)
              xticetm(ii)=xtice(ii,jj)
              obswtbtm(ii)=obswtbnow(ii,jj)
              tgtm(ii)=tg(ii,jj)
              dtfsit=dtsit
            endif

            fluxw(ii,jj)=-(1.-cicetm(ii))*(sstm(ii)-rstm(ii)-hfluxtm(ii)-qfluxtm(ii))
            fluxi(ii,jj)=-cicetm(ii)*(sstm(ii)-rstm(ii)-hfluxtm(ii)-qfluxtm(ii))
            soflw(ii,jj)=(1.-cicetm(ii))*sstm(ii)
            sofli(ii,jj)=cicetm(ii)*sstm(ii)
            wind10w(ii,jj)=sqrt(u10tm(ii)**2+v10tm(ii)**2)
            obsseaice(ii,jj)=cicetm(ii)
            evapw(ii,jj)=-qfluxtm(ii)/hltm
            rsf(ii,jj)=(rlsptm(ii)+rcuptm(ii))/dta

            !ustrw, vstrw
            rnuw = 1.14E-6
            z0w = 0.11*rnuw/ustartm(ii)+0.1*SQRT(rnuw*ustartm(ii)/grav) &
                   +0.015*ustartm(ii)**2/grav  ! Kraus & Businger(1994, page 146)
            d0w = 0.

            Pairvapor= AirVaporPressure(t2tm(ii),rh2tm(ii))
            rhoa = 100.*psttm(ii)/(287.04*t2tm(ii))*(1.0-0.378*Pairvapor/psttm(ii))
                 ! 287.04 [k.g-1.K-1]  Gas constant of dry air
            liw = -(0.4*grav*((hfluxtm(ii)/(t2tm(ii)*cp)+0.61*(qfluxtm(ii)/hltm)))/(ustartm(ii)**3.)*rhoa)
            ps1w = log((10.0-d0w)/z0w)-CalcSm((10.0-d0w)*liw,z0w,liw)+CalcSm(z0w*liw,z0w,liw)
            ustra = u10tm(ii)*0.4/ps1w
            vstra = v10tm(ii)*0.4/ps1w
            ustrw(ii,jj) = rhoa*ustra**2
            vstrw(ii,jj) = rhoa*vstra**2

            !in&out
            seaice(ii,jj)=cicetm(ii)
            thickness(ii,jj)=zicetm(ii)
            sni(ii,jj)=snrtm(ii)/1000.   !mm->m
            t_surf   =obswtbtm(ii)
            t_surf = MAX( t_surf, ctfreez )
            tsi(ii,jj)= MIN( t_surf, ctfreez )
            obswtb(ii,jj)=obswtbtm(ii)
            tsw(ii,jj)=tg(ii,jj)
            dtswdt(ii,jj)=0.

          endif    !end (lrun_sitvdiff)
        enddo    !end i=1,nxj

        if(lrun_sitvdiff) then
         call sit_vdiff ( nxjp(j), nxp, jj, istep, dtfsit,             &
              sitlat, sitlon(:,jj), tau, tauhr,                        &
              sitcor(:,jj), slm(:,jj), sitlclass(:,jj),                &
              sitmask(:,jj), bathy(:,jj), wlvl(:,jj),   &
              ocnmask(:,jj), obox_mask(:,jj),                          &
!            ! - same as lake and ml_ocean
              fluxw(:,jj), dfluxs(:,jj), soflw(:,jj),                  &
              fluxi(:,jj), sofli(:,jj),                                &
!            ! - 1D from mo_memory_g3b (wind stress)
              ustrw(:,jj), vstrw(:,jj),                                &
!            ! -    water mass variables(rain, snow, evap and runoff):
              rsf(:,jj), ssf(:,jj), evapw(:,jj), disch(:,jj),          &
!            ! - 1D from mo_memory_g3b
!              t2(:,j), wind10w(:,jj),                                 &
              t2tm(:), wind10w(:,jj),                               &
!            ! - 1D from mo_memory_g3b (sit variables)
!              obsseaice(:,jj), obswtbtm(:,jj), obswsb(:,jj),           &
              obsseaice(:,jj), obswtb(:,jj), obswsb(:,jj),           &
              sitwtb(:,jj), sitwub(:,jj), sitwvb(:,jj),                &
              sitwsb(:,jj), fluxiw(:,jj), pme2(:,jj),                  &
              subfluxw(:,jj), wsubsal(:,jj),                           &
              sitcc(:,jj), sithc(:,jj), engwac(:,jj),                  &
              sc(:,jj), saltwac(:,jj), wtfns(:,jj), wsfns(:,jj),       &
!            ! 3-d SIT vars: snow/ice
              zsi(:,jj,0:1), silw(:,jj,0:1), tsnic(:,jj,0:3),          &
!            ! 3-d SIT vars: water column
              obswt(:,jj,0:lkvl+1),obsws(:,jj,0:lkvl+1),               &
              obswu(:,jj,0:lkvl+1),obswv(:,jj,0:lkvl+1),               &
              sitwt(:,jj,0:lkvl+1),sitwu(:,jj,0:lkvl+1),               &
              sitwv(:,jj,0:lkvl+1),sitww(:,jj,0:lkvl+1),               &
              sitws(:,jj,0:lkvl+1),sitwtke(:,jj,0:lkvl+1),             &
              wlmx(:,jj,0:lkvl+1),wldisp(:,jj,0:lkvl+1),               &
              wkm(:,jj,0:lkvl+1),wkh(:,jj,0:lkvl+1),                   &
              wrho1000(:,jj,0:lkvl+1),wtfn(:,jj,0:lkvl+1),             &
              wsfn(:,jj,0:lkvl+1), wtfn0(:,jj,0:lkvl+1),               &
              wsfn0(:,jj,0:lkvl+1),awufl(:,jj,0:lkvl+1),               &
              awvfl(:,jj,0:lkvl+1),awtfl(:,jj,0:lkvl+1),               &
              awsfl(:,jj,0:lkvl+1),awtfl0(:,jj,0:lkvl+1),              &
              awsfl0(:,jj,0:lkvl+1),awtkefl(:,jj,0:lkvl+1),            &
!             ! final output only
!              cice(:,j), snr(:,j), zice(:,j), xtice(:,j), tg(:,j),    &
              seaice(:,jj),sni(:,jj),thickness(:,jj),xticetm(:),    &
              tsw(:,jj), tsl(:,jj), tslm(:,jj), tslm1(:,jj),           &
              ocu(:,jj), ocv(:,jj),  ctfreez2(:,jj),                   &
!            ! implicit with vdiff
              grndcapc(:,jj), grndhflx(:,jj), grndflux(:,jj),          &
              oldsitwt(:,jj,0:lkvl+1,0:1),oldsitwu(:,jj,0:lkvl+1,0:1), &
              oldsitwv(:,jj,0:lkvl+1,0:1),oldsitww(:,jj,0:lkvl+1,0:1), &
             oldsitws(:,jj,0:lkvl+1,0:1),oldsitwtke(:,jj,0:lkvl+1,0:1),&
              dtswdt(:,jj), sftobswt(:,jj,0:lkvl+1) )
        
          if(jj .eq. 1) then
            dtsittau=dtsittau+dtfsit
            dtsitmon=dtsitmon+dtfsit
            dtsit24=dtsit24+dtfsit
          endif
          call storesittau(nxjp(j),jj,nxp,my_max,lkvl,sitwt,sitws,sitwu,sitwv,dtfsit)
          call storesit24(nxjp(j),jj,nxp,my_max,lkvl,sitwt,sitws,sitwu,sitwv,dtfsit)

        endif  !end lrun_sitvdiff

        deallocate(sstm)
        deallocate(rstm)
        deallocate(hfluxtm)
        deallocate(qfluxtm)
        deallocate(u10tm)
        deallocate(v10tm)
        deallocate(rlsptm)
        deallocate(rcuptm)
        deallocate(ustartm)
        deallocate(t2tm)
        deallocate(rh2tm)
        deallocate(psttm)
        deallocate(cicetm)
        deallocate(snrtm)
        deallocate(zicetm)
        deallocate(xticetm)
        deallocate(obswtbtm)
        deallocate(tgtm)
        if(jj .eq. jlistnum) then
          if (lrun_sitvdiff) dtfsit=0.
        endif

      endif  !end do_sit


!=======================================================================
! Cal. SST tendency (dtsea/dt)
!=======================================================================
      dtx_tau=dt/3600.
      if(ldailyFCTsst .OR. ldailyFCTicesndpt .OR. dailyClm_option.ge.1) then
        if(itimestep .le. 1) then
          do ii = 1,nxj
            i=nxjstart(j)+ii-1
            dtseadt(ii,jj)=0.
            dFCTsstdt(ii,jj)=0. 
            obswtbnow(ii,jj)=dta*dFCTsstdt(ii,jj)+obswtbold(ii,jj)
            if(ocean(ii,jj))then
              tseadiffFCT(ii,jj)=dta*dFCTsstdt(ii,jj)
              tseanow(ii,jj)=dta*dFCTsstdt(ii,jj)+ tseaold(ii,jj)
            endif
          end do
        endif
        if((tau .ge. 24.)) then
          do ii = 1,nxj
            i=nxjstart(j)+ii-1
            dtseadt(ii,jj)=0.
            tseadiffFCT(ii,jj)=0.
            if(ocean(ii,jj))then
              obswtbnew(ii,jj)=dta*dFCTsstdt(ii,jj)+obswtbold(ii,jj)
              obswtbold(ii,jj)=obswtbnew(ii,jj)
              obswtbnow(ii,jj)=obswtbnew(ii,jj)
              ! sea surface temperature tendency 
              dtseadt(ii,jj)=dFCTsstdt(ii,jj)
              tseadiffFCT(ii,jj)=dta*dFCTsstdt(ii,jj)
              if(do_sit) then
                if(sitmask(ii,jj) .EQ. 1.) then
                  if(lrun_sitvdiff .AND. ltrigsit )then
                    tseadiffSIT(ii,jj)=0.
                    sumdSITdt(ii,jj)=sumdSITdt(ii,jj)+dtswdt(ii,jj)
                    countdSITdt(ii,jj)=countdSITdt(ii,jj)+1.
                  endif
                 
                  dtaup = mod(tau+0.001, dSITdt_intv)
                  if( (dSITdt_intv .lt. 0.) .OR. (dtaup .lt. dtx_tau) )then
                    if(countdSITdt(ii,jj) .ge. 1.) then
                      tseadiffFCT(ii,jj)=dta * (1.-weightSIT*ratioSIT(ii,jj)) * dFCTsstdt(ii,jj)
                      tseadiffSIT(ii,jj)=dta * weightSIT*ratioSIT(ii,jj) * (sumdSITdt(ii,jj)/countdSITdt(ii,jj))
                      tseadiffSIT24(ii,jj)=tseadiffSIT24(ii,jj)+tseadiffSIT(ii,jj)/dta*dt
                      ! sea surface temperature tendency
                      dtseadt(ii,jj)=(1.-weightSIT*ratioSIT(ii,jj)) * dFCTsstdt(ii,jj) &
                                  + weightSIT * ratioSIT(ii,jj) * (sumdSITdt(ii,jj)/countdSITdt(ii,jj))
                      sumdSITdt(ii,jj)=0.
                      countdSITdt(ii,jj)=0.
                    endif
                  endif
                endif ! end if sitmask(ii,jj) .EQ. 1
              endif ! end if do_sit

              if (dossst) then 
                dtseadt(ii,jj) = dtseadt(ii,jj) * (ssst3d(ii,1,jj) + 1.)
              endif 
              tseadiffFCT24(ii,jj)=tseadiffFCT24(ii,jj)+tseadiffFCT(ii,jj)/dta*dt
            endif !end if(ocean)
          end do  !end ii
        endif  !end if(tau .ge. 24.)
      endif !end if(ldailyFCTsst....)


!=======================================================================
! Add SPPT pertubation 
!=======================================================================
#ifdef VERBOSE
      ! Save u, v, t, and q for SPPT
      ut_update(1:nxj,1:lev,jj)      = ut(1:nxj,1:lev,jj)
      vt_update(1:nxj,1:lev,jj)      = vt(1:nxj,1:lev,jj)
      tt_update(1:nxj,1:lev,jj)      = tt(1:nxj,1:lev,jj)
      qt_update(1:nxj,1:lev*ncld,jj) = qt(1:nxj,1:lev*ncld,jj)
#endif
      ! sppt tendencies
      if (dosppt) then
        do k=1,lev
          kc=lev-k+1
          do i=1,nxj
            vru=1.
            if (kc.gt.zmtnblck(i)+2) then
               vru=1.0
            endif
            if (kc.le.zmtnblck(i)) then
               vru=0.0
            endif
            if (kc.eq.zmtnblck(i)+1) then
               vru=0.333333
            endif
            if (kc.eq.zmtnblck(i)+2) then
               vru=0.666667
            endif

            ru = sppt3d(i,k,jj) + 1.
            if (use_zmtnblck) ru = (ru - 1.) * vru + 1.

            dtdtr = dtradc(i,k,jj) * dta / 86400.
            upert = ( ut(i,k,jj) - ut_save_sppt(i,k,jj) ) * ru
            vpert = ( vt(i,k,jj) - vt_save_sppt(i,k,jj) ) * ru
            tpert = ( tt(i,k,jj) - tt_save_sppt(i,k,jj) - dtdtr ) * ru
            qpert = ( qt(i,k,jj) - qt_save_sppt(i,k,jj) ) * ru

            ut(i,k,jj) = ut_save_sppt(i,k,jj) + upert
            vt(i,k,jj) = vt_save_sppt(i,k,jj) + vpert

            !negative humidity check  
            qnew = qt_save_sppt(i,k,jj) + qpert
            if ( qnew .ge. qmin ) then
               qt(i,k,jj) = qnew
               tt(i,k,jj) = tt_save_sppt(i,k,jj) + tpert + dtdtr
            endif
          enddo
        enddo
      endif

      ! SHUM process 
      if (doshum) then
          do k=1,lev
            do i=1,nxj
              ru=shum3d(i,k,jj)
              qt(i,k,jj) = qt(i,k,jj)*(1.+ru)
            enddo
          enddo
      endif 

#ifdef VERBOSE
      ! Save u, v, t, and q for SPPT
      ut_sppt(1:nxj,1:lev     ,jj) = ut(1:nxj,1:lev     ,jj)
      vt_sppt(1:nxj,1:lev     ,jj) = vt(1:nxj,1:lev     ,jj)
      tt_sppt(1:nxj,1:lev     ,jj) = tt(1:nxj,1:lev     ,jj)
      qt_sppt(1:nxj,1:lev*ncld,jj) = qt(1:nxj,1:lev*ncld,jj)
#endif

    !--------------------------------------------------------------------------------
    !     weight back u and v by cosl/radus and
    !     change real temp back to virtual potential temperature
    !--------------------------------------------------------------------------------
      yy = cosl(j) / radus
      do  k = 1, lev
        do  i = 1, nxj
          ut(i,k,jj) = ut(i,k,jj)*yy
          vt(i,k,jj) = vt(i,k,jj)*yy
          tt(i,k,jj) = tt(i,k,jj)*(1.0+0.608*qt(i,k,jj))/pk(i,k,jj)
        enddo
      enddo
     
  290 continue  ! end of big j-loop for diabatic calculation

#ifdef VERBOSE
      if (myrank == 0) write(6,'(a18,a3,3(1x,a15))')'u wind','k','mean','variance','std. dev.'
      do k=1,lev
        call check_global_values(ut_save_sppt(:,k,:),k,'diabat before phys')
        call check_global_values(ut_pbl      (:,k,:),k,'diabat after   pbl')
        call check_global_values(ut_cmls     (:,k,:),k,'diabat after  cmls')
        call check_global_values(ut_update   (:,k,:),k,'diabat after  phys')
        call check_global_values(ut_sppt     (:,k,:),k,'diabat after  sppt')
      enddo
      if (myrank == 0) write(6,'(a18,a3,3(1x,a15))')'v wind','k','mean','variance','std. dev.'
      do k=1,lev
        call check_global_values(vt_save_sppt(:,k,:),k,'diabat before phys')
        call check_global_values(vt_pbl      (:,k,:),k,'diabat after   pbl')
        call check_global_values(vt_cmls     (:,k,:),k,'diabat after  cmls')
        call check_global_values(vt_update   (:,k,:),k,'diabat after  phys')
        call check_global_values(vt_sppt     (:,k,:),k,'diabat after  sppt')
      enddo
      if (myrank == 0) write(6,'(a18,a3,3(1x,a15))')'temperature','k','mean','variance','std. dev.'
      do k=1,lev
        call check_global_values(tt_save_sppt(:,k,:),k,'diabat before phys')
        call check_global_values(tt_pbl      (:,k,:),k,'diabat after   pbl')
        call check_global_values(tt_cmls     (:,k,:),k,'diabat after  cmls')
        call check_global_values(tt_update   (:,k,:),k,'diabat after  phys')
        call check_global_values(tt_sppt     (:,k,:),k,'diabat after  sppt')
      enddo
      if (myrank == 0) write(6,'(a18,a3,3(1x,a15))')'specific humidity','k','mean','variance','std. dev.'
      do k=1,lev
        call check_global_values(qt_save_sppt(:,k,:),k,'diabat before phys')
        call check_global_values(qt_pbl      (:,k,:),k,'diabat after   pbl')
        call check_global_values(qt_cmls     (:,k,:),k,'diabat after  cmls')
        call check_global_values(qt_update   (:,k,:),k,'diabat after  phys')
        call check_global_values(qt_sppt     (:,k,:),k,'diabat after  sppt')
      enddo
#endif 

!--------------------------------------------------------------------------------
!     update o3l to qt
!--------------------------------------------------------------------------------
      if(ntoz.gt.0)then
        do jj=1, jlistnum
           j=jlist1(jj)
           nxj=nxdef_2d(j)
           do k = 1, lev
             kk = (ntoz-1)*lev+k
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

      !  compute accumulated rain (mm)
      do jj =1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i = 1,nxj
          totalp(i,jj) = (rcup(i,jj)  + rlsp(i,jj))*rainfc
          raincu(i,jj) = raincu(i,jj) + rcup(i,jj)*rainfc
          rainlp(i,jj) = rainlp(i,jj) + rlsp(i,jj)*rainfc
          raincu1(i,jj)= raincu1(i,jj)+ rcup(i,jj)*rainfc
          rainlp1(i,jj)= rainlp1(i,jj)+ rlsp(i,jj)*rainfc
          raincu3(i,jj)= raincu3(i,jj)+ rcup(i,jj)*rainfc
          rainlp3(i,jj)= rainlp3(i,jj)+ rlsp(i,jj)*rainfc
          raincu6(i,jj)= raincu6(i,jj)+ rcup(i,jj)*rainfc
          rainlp6(i,jj)= rainlp6(i,jj)+ rlsp(i,jj)*rainfc
          raintot(i,jj)= raintot(i,jj)+ totalp(i,jj)
        enddo
      enddo

      !  reduce totalp in first hour forecast to avoid apin-up problem
      if ( tau .le. 1.001)  then
        dtau = dt/3600.0
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do i = 1,nxj
            totalp(i,jj) = (tau-dtau)* totalp(i,jj)
          enddo
       enddo
      endif
!
!      call mpe_unify(totalp,nx,my,2,mpe_double)
!
      if(ldiag.ge.1) then
!     print level-1 diagnostics
!                   (aps,psx,psn,nncup,nnshl,nndry,nnlsp,itlsp,arcup
!                   arlsp,xkmx,dtcupx and layer mean dtcup and dqcup )
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
!      qbrw = 0.0
!
!     call mpe_unify(qbrrow,ncld,my,2,mpe_double)
!      call mpe_unify(qbrrow,ncld,my,4,mpe_double)
!
!      do 320 j = 1, my
!      qbrw = qbrw + qbrrow(1,j)
!  320 continue
!
      cosq = cosw*nx
!
      qglb  = qglb * 100.0/grav/cosq
!      qbrw  = qbrw * 100.0/grav/cosq
      rainbl= ( arcup + arlsp - evapor ) * dta / 86400.0
!      qdiff = qgini - (qglb - rainbl - qbrw)
      engdiff=tengi-teng
      thdadif=thdai-thda
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
      call mpe_global_maxloc(xkmkk, 2,ixkmkk,jxkmkk,idummy,idummy,mpe_double)
      call mpe_global_maxloc(xkmk1, 2,ixkmk1,jxkmk1,idummy,idummy,mpe_double)
      call mpe_global_maxloc(dtcupg,2,idcupx,jdcupx,idummy,idummy,mpe_double)
      call mpe_global_maxloc(utx,   2,ikutx, jutx  ,idummy,idummy,mpe_double)
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
      call mpe_global_maxloc(dtradg,4,isign,idradx,jdradx,kdradx,mpe_double)
!
      dtradg = dtradg*isign
!
!      if(docgrav)then
!        call mpe_unify(avgdrag_u,my,lev,1,mpe_double)
!        call maxpp2l( lev,my,avgdrag_u,dragmax,kmax,jmax )
!        call minpp2l( lev,my,avgdrag_u,dragmin,kmin,jmin )
!        if(myrank .eq. 0) then
!          print *,'norographic gravity wave drag u_max,kmax,jmax,u_min,kmin,jmin'
!          print *, dragmax,kmax,jmax,dragmin,kmin,jmin
!        endif
!        call mpe_unify(avgdrag_v,my,lev,1,mpe_double)
!        call maxpp2l( lev,my,avgdrag_v,dragmax,kmax,jmax )
!        call minpp2l( lev,my,avgdrag_v,dragmin,kmin,jmin )
!        if(myrank .eq. 0) then
!          print *,'norographic gravity wave drag v_max,kmax,jmax,v_min,kmin,jmin'
!          print *, dragmax,kmax,jmax,dragmin,kmin,jmin
!        endif
!      endif
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
      print *,' global moisture budget: qgini, qglb ='
      print *, qgini,qglb
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
 8019 format( 1x,' global moisture budget: qgini, qglb =' &
            , 2e14.6 )
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
  end subroutine diabat
