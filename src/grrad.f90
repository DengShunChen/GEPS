!!!!!  ==========================================================  !!!!!
!!!!!             'module_radiation_driver' descriptions           !!!!!
!!!!!  ==========================================================  !!!!!
!                                                                      !
!   this is the radiation driver module.  it prepares atmospheric      !
!   profiles and invokes main radiation calculations.                  !
!                                                                      !
!   in module 'module_radiation_driver' there are twe externally       !
!   callable subroutine:                                               !
!                                                                      !
!      'radinit'    -- initialization routine                          !
!         input:                                                       !
!           ( si, nlay, me )                                           !
!         output:                                                      !
!           ( none )                                                   !
!                                                                      !
!      'radupdate'  -- update time sensitive data used by radiations   !
!         input:                                                       !
!           ( idate,jdate,deltsw,deltim,lsswr, me )                    !
!         output:                                                      !
!           ( slag,sdec,cdec,solcon )                                  !
!                                                                      !
!      'grrad'      -- setup and invoke main radiation calls           !
!         input:                                                       !
!          ( prsi,prsl,prslk,tgrs,qgrs,tracer,vvl,slmsk,               !
!            xlon,xlat,tsfc,snowd,sncovr,snoalb,zorl,hprim,            !
!            alvsf,alnsf,alvwf,alnwf,facsf,facwf,                      !
!            sinlat,coslat,solhr,jdate,solcon,                         !
!            cv,cvt,cvb,                                               !
!            icsdsw,icsdlw,ntcw,ncld,ntoz,ntrac,nfxr,                 !
!            dtlw,dtsw,lsswr,lslwr,lssav,                              !
!            ix, im, lm, me, lprnt, ipt, kdt,                          !
!         output:                                                      !
!            htrsw,topfsw,sfcfsw,sfalb,coszen,coszdg,                  !
!            htrlw,topflw,sfcflw,tsflw,semis,cldcov,                   !
!         input/output:                                                !
!            fluxr                                                     !
!         optional output:                                             !
!            htrlw0,htrsw0,htrswb,htrlwb                               !
!                                                                      !
!                                                                      !
!   external modules referenced:                                       !
!       'module physpara'                   in 'physpara.f'            !
!       'module funcphys'                   in 'funcphys.f'            !
!       'module physcons'                   in 'physcons.f'            !
!                                                                      !
!       'module module_radiation_gases'     in 'radiation_gases.f'     !
!       'module module_radiation_aerosols'  in 'radiation_aerosols.f'  !
!       'module module_radiation_surface'   in 'radiation_surface.f'   !
!       'module module_radiation_clouds'    in 'radiation_clouds.f'    !
!                                                                      !
!       'module module_radsw_cntr_para'     in 'radsw_xxxx_param.f'    !
!       'module module_radsw_parameters'    in 'radsw_xxxx_param.f'    !
!       'module module_radsw_main'          in 'radsw_xxxx_main.f'     !
!                                                                      !
!       'module module_radlw_cntr_para'     in 'radlw_xxxx_param.f'    !
!       'module module_radlw_parameters'    in 'radlw_xxxx_param.f'    !
!       'module module_radlw_main'          in 'radlw_xxxx_main.f'     !
!                                                                      !
!    where xxxx may vary according to different scheme selection       !
!                                                                      !
!                                                                      !
!   program history log:                                               !
!     mm-dd-yy    ncep         - created program grrad                 !
!     08-12-03    yu-tai hou   - re-written for modulized radiations   !
!     11-06-03    yu-tai hou   - modified                              !
!     01-18-05    s. moorthi   - noah/ice model changes added          !
!     05-10-05    yu-tai hou   - modified module structure             !
!     12-xx-05    s. moorthi   - sfc lw flux adj by mean temperature   !
!     02-20-06    yu-tai hou   - add time variation for co2 data, and  !
!                                solar const. add sfc emiss change     !
!     03-21-06    s. moorthi   - added surface temp over ice           !
!     07-28-06    yu-tai hou   - add stratospheric vocanic aerosols    !
!     03-14-07    yu-tai hou   - add generalized spectral band interp  !
!                                for aerosol optical prop. (sw and lw) !
!     04-10-07    yu-tai hou   - spectral band sw/lw heating rates     !
!     05-04-07    yu-tai hou   - make options for clim based and modis !
!                                based (h. wei and c. marshall) albedo !
!     09-05-08    yu-tai hou   - add the initial date and time 'idate' !
!                    and control param 'ictm' to the passing param list!
!                    to handel different time/date requirements for    !
!                    external data (co2, aeros, solcon, ...)           !
!     10-10-08    yu-tai hou   - add the ictm=-2 option for combining  !
!                    initial condition data with seasonal cycle from   !
!                    climatology.                                      !
!     03-12-09    yu-tai hou   - use two time stamps to keep tracking  !
!                    dates for init cond and fcst time. remove volcanic!
!                    aerosols data in climate hindcast (ictm=-2).      !
!     03-16-09    yu-tai hou   - included sub-column clouds approx.    !
!                    control flags isubcsw/isubclw in initialization   !
!                    subroutine. passed auxiliary cloud control arrays !
!                    icsdsw/icsdlw (if isubcsw/isubclw =2, it will be  !
!                    the user provided permutation seeds) to the sw/lw !
!                    radiation calculation programs. also moved cloud  !
!                    overlapping control flags iovrsw/iovrlw from main !
!                    radiation routines to the initialization routines.!
!     04-02-09    yu-tai hou   - modified surface control flag iems to !
!                    have additional function of if the surface-air    !
!                    interface have the same or different temperature  !
!                    for radiation calculations.                       !
!     04-03-09    yu-tai hou   - modified to add lw surface emissivity !
!                    as output variable. changed the sign of sfcnsw to !
!                    be positive value denote solar flux goes into the !
!                    ground (this is needed to reduce sign confusion   !
!                    in other part of model)                           !
!     09-09-09    fanglin yang (thru s.moorthi) added qme5 qme6 to e-20!
!     01-09-10    sarah lu     - added gocart option, revised grrad for!
!                    gocart coupling. calling argument modifed: ldiag3 !
!                    removed; cldcov/fluxr sequence changed; cldcov is !
!                    changed from accumulative to instant field and    !
!                    from input/output to output field                 !
!     01-24-10    sarah lu     - added aod to fluxr, added prslk and   !
!                    oz to setaer input argument (for gocart coupling),!
!                    added tau_gocart to setaer output argument (for,  !
!                    aerosol diag by index of nv_aod)                  !
!     07-08-10    s.moorthi - updated the nems version for new physics !
!     07-28-10    yu-tai hou   - changed grrad interface to allow all  !
!                    components of sw/lw toa/sfc instantaneous values  !
!                    being passed to the calling program. moved the    !
!                    computaion of sfc net sw flux (sfcnsw) to the     !
!                    calling program. merged carlos' nmmb modification.!
!     07-30-10    s. moorthi - corrected some errors associated with   !
!                    unit changes                                      !
!     12-02-10    s. moorthi/y. hou - removed the use of aerosol flags !
!                    'iaersw' 'iaerlw' from radiations and replaced    !
!                    them by using the runtime variable laswflg and    !
!                    lalwflg defined in module radiation_aerosols.     !
!                    also replaced param nspc in grrad with the use of !
!                    max_num_gridcomp in module radiation_aerosols.    !
!     jun 2012    yu-tai hou   - added sea/land madk 'slmsk' to the    !
!                    argument list of subrotine setaer call for the    !
!                    newly modified horizontal bi-linear interpolation !
!                    in climatological aerosols schem. also moved the  !
!                    virtual temperature calculations in subroutines   !
!                    'radiation_clouds' and 'radiation_aerosols' to    !
!                    'grrad' to reduce repeat comps. renamed var oz as !
!                    tracer to reflect that it carries various prog    !
!                    tracer quantities.                                !
!                              - modified to add 4 compontents of sw   !
!                    surface downward fluxes to the output. (vis/nir;  !
!                    direct/diffused). re-arranged part of the fluxr   !
!                    variable fields and filled the unused slots for   !
!                    the new components.  added check print of select  !
!                    data (co2 value for now).                         !
!                              - changed the initialization subrution  !
!                    'radinit' into two parts: 'radinit' is called at  !
!                    the start of model run to set up radiation related!
!                    fixed parameters; and 'radupdate' is called in    !
!                    the time-loop to update time-varying data sets    !
!                    and module variables.                             !
!     sep 2012    h-m lin/y-t hou added option of extra top layer for  !
!                    models with low toa ceiling. the extra layer will !
!                    help ozone absorption at higher altitude.         !
!     nov 2012    yu-tai hou   - modified control parameters through   !
!                    module 'physpara'.                                !
!     jan 2013    yu-tai hou   - updated subr radupdate for including  !
!                    options of annual/monthly solar constant table.   !
!     mar 2013    h-m lin/y-t hou corrected a bug in extra top layer   !
!                    when using ferrier microphysics.                  !
!     may 2013    s. mooorthi - removed fpkapx                         !
!     jul 2013    r. sun - added pdf cld and convective cloud water and! 
!	             cover for radiation                               ! 
!                                                                      !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!



!========================================!
      module module_radiation_driver     !
!........................................!
!
      use physpara
      use physcons,                 only : con_eps, con_epsm1, con_fvirt&
     &,                                    rocp => con_rocp
!     use funcphys,                 only : fpvs

      use module_radiation_astronomy,only: sol_init, sol_update, coszmn
      use module_radiation_gases,   only : nf_vgas, getgases, getozn,   &
     &                                     gas_init, gas_update
      use module_radiation_aerosols,only : nf_aesw, nf_aelw, setaer,    &
     &                                     aer_init, aer_update
!    &,                                    nspc1                        ! optn for aod output
      use module_radiation_surface, only : nf_albd, sfc_init, setalb,   &
     &                                     setemis
      use module_radiation_clouds,  only : nf_clds, cld_init,           &
     &                                     progcld1, progcld2, progcld3,&
     &					   progcld4, diagcld1,          &
                                           progcld5, progcld5o,         &
                                           progclduni, progcld6,        &
                                           progcld_thompson, progcld_gce

      use module_radsw_parameters,  only : topfsw_type, sfcfsw_type,    &
     &                                     profsw_type,cmpfsw_type,nbdsw
      use module_radsw_main,        only : rswinit,  swrad

      use module_radlw_parameters,  only : topflw_type, sfcflw_type,    &
     &                                     proflw_type, nbdlw
      use module_radlw_main,        only : rlwinit,  lwrad
!
!    implicit   none
!
!
      private

!  ---  version tag and last revision date
      character(40), parameter ::                                       &
     &   vtagrad='ncep-radiation_driver    v5.2  jan 2013 '
!    &   vtagrad='ncep-radiation_driver    v5.1  nov 2012 '
!    &   vtagrad='ncep-radiation_driver    v5.0  aug 2012 '

!  ---  constant values
      real (kind=kind_phys) :: qmin, qme5, qme6, epsq
!     parameter (qmin=1.0e-10, qme5=1.0e-5,  qme6=1.0e-6,  epsq=1.0e-12)
      parameter (qmin=1.0e-10, qme5=1.0e-7,  qme6=1.0e-7,  epsq=1.0e-12)
!     parameter (qmin=1.0e-10, qme5=1.0e-20, qme6=1.0e-20, epsq=1.0e-12)
      real, parameter :: prsmin = 1.0e-6 ! toa pressure minimum value in mb (hpa)

!  ---  control flags set in subr radinit:
      integer :: itsfc  =0            ! flag for lw sfc air/ground interface temp setting

!  ---  data input control variables set in subr radupdate:
      integer :: month0=0,   iyear0=0,   monthd=0
      logical :: loz1st =.true.       ! first-time clim ozone data read flag

!  ---  optional extra top layer on top of low ceiling models
      integer, parameter :: ltp = 0   ! no extra top layer
!     integer, parameter :: ltp = 1   ! add an extra top layer
      logical, parameter :: lextop = (ltp > 0)


!  ---  publicly accessible module programs:

      public radinit, radupdate, grrad


! =================
      contains
! =================


!-----------------------------------
      subroutine radinit                                                &
!...................................

!  ---  inputs:
     &     ( si, nlay, me, myrank )
!  ---  outputs:
!          ( none )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:   radinit     initialization of radiation calculations    !
!                                                                       !
! usage:        call radinit                                            !
!                                                                       !
! attributes:                                                           !
!   language:  fortran 90                                               !
!   machine:   ibm sp                                                   !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
! input parameters:                                                     !
!   nlay             : number of model vertical layers                  !
!   me               : print control flag                               !
!                                                                       !
!  outputs: (none)                                                      !
!                                                                       !
!  external module variables:  (in module physpara)                     !
!   isolar   : solar constant cntrol flag                               !
!              = 0: use the old fixed solar constant in "physcon"       !
!              =10: use the new fixed solar constant in "physcon"       !
!              = 1: use noaa ann-mean tsi tbl abs-scale with cycle apprx!
!              = 2: use noaa ann-mean tsi tbl tim-scale with cycle apprx!
!              = 3: use cmip5 ann-mean tsi tbl tim-scale with cycl apprx!
!              = 4: use cmip5 mon-mean tsi tbl tim-scale with cycl apprx!
!   iaerflg  : 3-digit aerosol flag (abc for volc, lw, sw)              !
!              a:=0 use background stratospheric aerosol                !
!                =1 include stratospheric vocanic aeros                 !
!              b:=0 no topospheric aerosol in lw radiation              !
!                =1 compute tropspheric aero in 1 broad band for lw     !
!                =2 compute tropspheric aero in multi bands for lw      !
!              c:=0 no topospheric aerosol in sw radiation              !
!                =1 include tropspheric aerosols for sw                 !
!   ico2flg  : co2 data source control flag                             !
!              =0: use prescribed global mean co2 (old  oper)           !
!              =1: use observed co2 annual mean value only              !
!              =2: use obs co2 monthly data with 2-d variation          !
!   ictmflg  : =yyyy#, external data ic time/date control flag          !
!              =   -2: same as 0, but superimpose seasonal cycle        !
!                      from climatology data set.                       !
!              =   -1: use user provided external data for the          !
!                      forecast time, no extrapolation.                 !
!              =    0: use data at initial cond time, if not            !
!                      available, use latest, no extrapolation.         !
!              =    1: use data at the forecast time, if not            !
!                      available, use latest and extrapolation.         !
!              =yyyy0: use yyyy data for the forecast time,             !
!                      no further data extrapolation.                   !
!              =yyyy1: use yyyy data for the fcst. if needed, do        !
!                      extrapolation to match the fcst time.            !
!   ioznflg  : ozone data source control flag                           !
!              =0: use climatological ozone profile                     !
!              =1: use interactive ozone profile                        !
!   ialbflg  : albedo scheme control flag                               !
!              =0: climatology, based on surface veg types              !
!              =1: modis retrieval based surface albedo scheme          !
!   iemsflg  : emissivity scheme cntrl flag (ab 2-digit integer)        !
!              a:=0 set sfc air/ground t same for lw radiation          !
!                =1 set sfc air/ground t diff for lw radiation          !
!              b:=0 use fixed sfc emissivity=1.0 (black-body)           !
!                =1 use varying climtology sfc emiss (veg based)        !
!                =2 future development (not yet)                        !
!   icldflg  : cloud optical property scheme control flag               !
!              =0: use diagnostic cloud scheme                          !
!              =1: use prognostic cloud scheme (default)                !
!   icmphys  : cloud microphysics scheme control flag                   !
!              =1 zhao/carr/sundqvist microphysics scheme               !
!              =2 brad ferrier microphysics scheme                      !
!	       =3 zhao/carr/sundqvist microphysics+pdf cloud & cnvc,cnvw!
!   iovrsw   : control flag for cloud overlap in sw radiation           !
!   iovrlw   : control flag for cloud overlap in lw radiation           !
!              =0: random overlapping clouds                            !
!              =1: max/ran overlapping clouds                           !
!   isubcsw  : sub-column cloud approx control flag in sw radiation     !
!   isubclw  : sub-column cloud approx control flag in lw radiation     !
!              =0: with out sub-column cloud approximation              !
!              =1: mcica sub-col approx. prescribed random seed         !
!              =2: mcica sub-col approx. provided random seed           !
!   lcrick   : control flag for eliminating crick                       !
!              =t: apply layer smoothing to eliminate crick             !
!              =f: do not apply layer smoothing                         !
!   lcnorm   : control flag for in-cld condensate                       !
!              =t: normalize cloud condensate                           !
!              =f: not normalize cloud condensate                       !
!   lnoprec  : precip effect in radiation flag (ferrier microphysics)   !
!              =t: snow/rain has no impact on radiation                 !
!              =f: snow/rain has impact on radiation                    !
!   ivflip   : vertical index direction control flag                    !
!              =0: index from toa to surface                            !
!              =1: index from surface to toa                            !
!                                                                       !
!  subroutines called: sol_init, aer_init, gas_init, cld_init,          !
!                      sfc_init, rlwinit, rswinit                       !
!                                                                       !
!  usage:       call radinit                                            !
!                                                                       !
!  ===================================================================  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: nlay, me, myrank
   
      real (kind=kind_phys), intent(in) :: si(:)

!  ---  outputs: (none, to module variables)

!  ---  locals:

!
!===> ...  begin here
!
!  ---  set up control variables
      itsfc  = iemsflg / 10             ! sfc air/ground temp control
      loz1st = (ioznflg == 0)           ! first-time clim ozone data read flag
      month0 = 0
      iyear0 = 0
      monthd = 0

      if (me == 0 .and. myrank == 0) then
!       print *,' new radiation program structures -- sep 01 2004'
        print *,' new radiation program structures became oper. ',      &
     &          '  may 01 2007'
        print *, vtagrad                !print out version tag
        print *,' - selected control flag settings: ictmflg=',ictmflg,  &
     &    ' isolar =',isolar, ' ico2flg=',ico2flg,' iaerflg=',iaerflg,  &
     &    ' ialbflg=',ialbflg,' iemsflg=',iemsflg,' icldflg=',icldflg,  &
     &    ' icmphys=',icmphys,' ioznflg=',ioznflg
        print *,' ivflip=',ivflip,' iovrsw=',iovrsw,' iovrlw=',iovrlw,  &
     &    ' isubcsw=',isubcsw,' isubclw=',isubclw
        print *,' lcrick=',lcrick,' lcnorm=',lcnorm,' lnoprec=',lnoprec
        print *,' ltp =',ltp,', add extra top layer =',lextop

        if ( ictmflg==0 .or. ictmflg==-2 ) then
          print *,'   data usage is limited by initial condition!'
          print *,'   no volcanic aerosols'
        endif

        if ( isubclw == 0 ) then
          print *,' - isubclw=',isubclw,' no mcica, use grid ',         &
     &            'averaged cloud in lw radiation'
        elseif ( isubclw == 1 ) then
          print *,' - isubclw=',isubclw,' use mcica with fixed ',       &
     &            'permutation seeds for lw random number generator'
        elseif ( isubclw == 2 ) then
          print *,' - isubclw=',isubclw,' use mcica with random ',      &
     &            'permutation seeds for lw random number generator'
        else
          print *,' - error!!! isubclw=',isubclw,' is not a ',          &
     &            'valid option '
          stop
        endif

        if ( isubcsw == 0 ) then
          print *,' - isubcsw=',isubcsw,' no mcica, use grid ',         &
     &            'averaged cloud in sw radiation'
        elseif ( isubcsw == 1 ) then
          print *,' - isubcsw=',isubcsw,' use mcica with fixed ',       &
     &            'permutation seeds for sw random number generator'
        elseif ( isubcsw == 2 ) then
          print *,' - isubcsw=',isubcsw,' use mcica with random ',      &
     &            'permutation seeds for sw random number generator'
        else
          print *,' - error!!! isubcsw=',isubcsw,' is not a ',          &
     &            'valid option '
          stop
        endif

        if ( isubcsw /= isubclw ) then
          print *,' - *** notice *** isubcsw /= isubclw !!!',           &
     &            isubcsw, isubclw
        endif
      endif

!  --- ...  call astronomy initialization routine

      if (me == 0 .and. myrank ==0)print *,'call sol_init' 
      call sol_init ( me , myrank )

!  --- ...  call aerosols initialization routine

      if (me == 0 .and. myrank == 0)print *,'call aer_init' 
      call aer_init ( nlay, me , myrank )

!  --- ...  call co2 and other gases initialization routine

      if (me ==0 .and. myrank == 0)print *,'call gas_init' 
      call gas_init ( me , myrank )

!  --- ...  call surface initialization routine

      if (me == 0 .and. myrank == 0)print *,'call sfc_init' 
      call sfc_init ( me , myrank )

!  --- ...  call cloud initialization routine

      if (me == 0 .and. myrank == 0)print *,'call cld_init' 
      call cld_init ( si, nlay, me, myrank )

!  --- ...  call lw radiation initialization routine

      if (me == 0 .and. myrank == 0)print *,'call rlw_init' 
      call rlwinit ( me , myrank )

!  --- ...  call sw radiation initialization routine

      if (me == 0 .and. myrank == 0)print *,'call rsw_init' 
      call rswinit ( me ,myrank)
!
      return
!...................................
      end subroutine radinit
!-----------------------------------


!-----------------------------------
      subroutine radupdate                                              &
!...................................
!  ---  inputs:
     &     ( idate,jdate,deltsw,deltim,lsswr, me, myrank,               &          
!  ---  outputs:
     &       slag,sdec,cdec,solcon                                      &
     &     )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:   radupdate   calls many update subroutines to check and  !
!   update radiation required but time varying data sets and module     !
!   variables.                                                          !
!                                                                       !
! usage:        call radupdate                                          !
!                                                                       !
! attributes:                                                           !
!   language:  fortran 90                                               !
!   machine:   ibm sp                                                   !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
! input parameters:                                                     !
!   idate(8)       : ncep absolute date and time of initial condition   !
!                    (yr, mon, day, t-zone, hr, min, sec, mil-sec)      !
!   jdate(8)       : ncep absolute date and time at fcst time           !
!                    (yr, mon, day, t-zone, hr, min, sec, mil-sec)      !
!   deltsw         : sw radiation calling frequency in seconds          !
!   deltim         : model timestep in seconds                          !
!   lsswr          : logical flags for sw radiation calculations        !
!   me             : print control flag                                 !
!                                                                       !
!  outputs:                                                             !
!   slag           : equation of time in radians                        !
!   sdec, cdec     : sin and cos of the solar declination angle         !
!   solcon         : sun-earth distance adjusted solar constant (w/m2)  !
!                                                                       !
!  external module variables:                                           !
!   isolar   : solar constant cntrl  (in module physpara)               !
!              = 0: use the old fixed solar constant in "physcon"       !
!              =10: use the new fixed solar constant in "physcon"       !
!              = 1: use noaa ann-mean tsi tbl abs-scale with cycle apprx!
!              = 2: use noaa ann-mean tsi tbl tim-scale with cycle apprx!
!              = 3: use cmip5 ann-mean tsi tbl tim-scale with cycl apprx!
!              = 4: use cmip5 mon-mean tsi tbl tim-scale with cycl apprx!
!   ictmflg  : =yyyy#, external data ic time/date control flag          !
!              =   -2: same as 0, but superimpose seasonal cycle        !
!                      from climatology data set.                       !
!              =   -1: use user provided external data for the          !
!                      forecast time, no extrapolation.                 !
!              =    0: use data at initial cond time, if not            !
!                      available, use latest, no extrapolation.         !
!              =    1: use data at the forecast time, if not            !
!                      available, use latest and extrapolation.         !
!              =yyyy0: use yyyy data for the forecast time,             !
!                      no further data extrapolation.                   !
!              =yyyy1: use yyyy data for the fcst. if needed, do        !
!                      extrapolation to match the fcst time.            !
!                                                                       !
!  module variables:                                                    !
!   loz1st   : first-time clim ozone data read flag                     !
!                                                                       !
!  subroutines called: sol_update, aer_update, gas_update               !
!                                                                       !
!  ===================================================================  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: idate(:), jdate(:), me, myrank
      logical, intent(in) :: lsswr

      real (kind=kind_phys), intent(in) :: deltsw, deltim

!  ---  outputs:
      real (kind=kind_phys), intent(out) :: slag, sdec, cdec, solcon

!  ---  locals:
      integer :: iyear, imon, iday, ihour
      integer :: kyear, kmon, kday, khour

      logical :: lmon_chg       ! month change flag
      logical :: lco2_chg       ! cntrl flag for updating co2 data
      logical :: lsol_chg       ! cntrl flag for updating solar constant
!
!===> ...  begin here
!
!  --- ...  time stamp at fcst time

      iyear = jdate(1)
      imon  = jdate(2)
      iday  = jdate(3)
      ihour = jdate(5)

!  --- ...  set up time stamp used for green house gases (** currently co2 only)

      if ( ictmflg==0 .or. ictmflg==-2 ) then  ! get external data at initial condition time
        kyear = idate(1)
        kmon  = idate(2)
        kday  = idate(3)
        khour = idate(5)
      else                           ! get external data at fcst or specified time
        kyear = iyear
        kmon  = imon
        kday  = iday
        khour = ihour
      endif   ! end if_ictmflg_block

      if ( month0 /= imon ) then
        lmon_chg = .true.
        month0 = imon
      else
        lmon_chg = .false.
      endif

!  --- ...  call astronomy update routine, yearly update, no time interpolation

      if (lsswr) then

        if ( isolar == 0 .or. isolar == 10 ) then
          lsol_chg = .false.
        elseif ( iyear0 /= iyear ) then
          lsol_chg = .true.
        else
          lsol_chg = ( isolar==4 .and. lmon_chg )
        endif
        iyear0 = iyear

        if (me == 0 .and. myrank ==0) print *,'call sol_update'
        call sol_update                                                 &
!  ---  inputs:
     &     ( myrank, jdate,kyear,deltsw,deltim,lsol_chg, me,            &
!  ---  outputs:
     &       slag,sdec,cdec,solcon                                      &
     &     )

      endif  ! end_if_lsswr_block

!  --- ...  call aerosols update routine, monthly update, no time interpolation

      if ( lmon_chg ) then
        if (me == 0 .and. myrank ==0) print *,'call aer_update'
        call aer_update ( myrank, iyear, imon, me )
      endif

!  --- ...  call co2 and other gases update routine

      if ( monthd /= kmon ) then
        monthd = kmon
        lco2_chg = .true.
      else
        lco2_chg = .false.
      endif

      if (me == 0 .and. myrank ==0) print *,'call gas_update'
      call gas_update ( myrank, kyear,kmon,kday,khour,loz1st,lco2_chg, me )

      if ( loz1st ) loz1st = .false.

!  --- ...  call surface update routine (currently not needed)
!     call sfc_update ( iyear, imon, me )

!  --- ...  call clouds update routine (currently not needed)
!     call cld_update ( iyear, imon, me )
!
      return
!...................................
      end subroutine radupdate
!-----------------------------------


!-----------------------------------
      subroutine grrad                                                  &
!...................................
!  ---  inputs:
           ( prsi,prsl,prslk,tgrs,qgrs,tracer,vvl,slmsk,                &
             xlon,xlat,tsfc,snowd,sncovr,snoalb,zorl,hprim,             &
             alvsf,alnsf,alvwf,alnwf,facsf,facwf,fice,tisfc,            &
             sinlat,coslat,solhr,jdate,solcon,                          &
             cv,cvt,cvb,                                                &
             icsdsw,icsdlw,ntcw,ncld,ntoz,ntrac,nfxr,                   &
             dtlw,dtsw,lsswr,lslwr,lssav,                               &
             ix,im,lm,me,lprnt,ipt,kdt,myrank,                          &
             ntiw,ntrw,ntsw,ntgl,uni_cloud,lmfshal,lmfdeep2,            &
             deltaq,sup,cnvw,cnvc,phy_f3d,                              &
!  ---  outputs:
             htrsw,sfalb,coszen,coszdg,                                 &
             htrlw,tsflw,semis,                                         &
!  ---  input/output:
             cldcov,fluxr,                                              &
!! ---  optional outputs:
             htrsw0,htrlw0,                                             &
             fusl,fdsl,fuir,fdir,                                       &
             fuslr,fdslr,fuirr,fdirr                                    &
           )

! =================   subprogram documentation block   ================ !
!                                                                       !
!    this program is the driver of radiation calculation subroutines. * !
!    it sets up profile variables for radiation input, including      * !
!    clouds, surface albedos, atmospheric aerosols, ozone, etc.       * !
!                                                                     * !
!    usage:        call grrad                                         * !
!                                                                     * !
!    subprograms called:                                              * !
!                  setalb, setemis, setaer, getozn, getgases,         * !
!                  progcld1, progcld2, diagcds,                       * !
!                  swrad, lwrad, fpvs                                 * !
!                                                                     * !
!    attributes:                                                      * !
!      language:   fortran 90                                         * !
!      machine:    ibm-sp, sgi                                        * !
!                                                                     * !
!                                                                     * !
!  ====================  defination of variables  ====================  !
!                                                                       !
!    input variables:                                                   !
!      prsi  (ix,lm+1) : model level pressure in cb (kpa)               !
!      prsl  (ix,lm)   : model layer mean pressure in cb (kpa)          !
!      prslk (ix,lm)   : exner function = (p/p0)**rocp                  !
!      tgrs  (ix,lm)   : model layer mean temperature in k              !
!      qgrs  (ix,lm)   : layer specific humidity in gm/gm               !
!      tracer(ix,lm,ntrac):layer prognostic tracer amount/mixing-ratio  !
!                        incl: oz, cwc, aeros, etc.                     !
!      vvl   (ix,lm)   : layer mean vertical velocity in cb/sec         !
!      slmsk (im)      : sea/land mask array (sea:0,land:1,sea-ice:2)   !
!      xlon  (im)      : grid longitude in radians, ok for both 0->2pi  !
!                        or -pi -> +pi ranges                           !
!      xlat  (im)      : grid latitude in radians, default to pi/2 ->   !
!                        -pi/2 range, otherwise adj in subr called      !
!      tsfc  (im)      : surface temperature in k                       !
!      snowd (im)      : snow depth water equivalent in mm              !
!      sncovr(im)      : snow cover in fraction                         !
!      snoalb(im)      : maximum snow albedo in fraction                !
!      zorl  (im)      : surface roughness in cm                        !
!      hprim (im)      : topographic standard deviation in m            !
!      alvsf (im)      : mean vis albedo with strong cosz dependency    !
!      alnsf (im)      : mean nir albedo with strong cosz dependency    !
!      alvwf (im)      : mean vis albedo with weak cosz dependency      !
!      alnwf (im)      : mean nir albedo with weak cosz dependency      !
!      facsf (im)      : fractional coverage with strong cosz dependen  !
!      facwf (im)      : fractional coverage with weak cosz dependency  !
!      fice  (im)      : ice fraction over open water grid              !
!      tisfc (im)      : surface temperature over ice fraction          !
!      sinlat(im)      : sine of the grids' corresponding latitudes     !
!      coslat(im)      : cosine of the grids' corresponding latitudes   !
!      solhr           : hour time after 00z at the t-stepe             !
!      jdate (8)       : current forecast date and time                 !
!                        (yr, mon, day, t-zone, hr, min, sec, mil-sec)  !
!      solcon          : solar constant (sun-earth distant adjusted)    !
!      cv    (im)      : fraction of convective cloud                   !
!      cvt, cvb (im)   : convective cloud top/bottom pressure in cb     !
!      fcice           : fraction of cloud ice  (in ferrier scheme)     !
!      frain           : fraction of rain water (in ferrier scheme)     !
!      rrime           : mass ratio of total to unrimed ice ( >= 1 )    !
!      flgmin          : minimim large ice fraction                     !
!      icsdsw/icsdlw   : auxiliary cloud control arrays passed to main  !
!           (im)         radiations. if isubcsw/isubclw (input to init) !
!                        are set to 2, the arrays contains provided     !
!                        random seeds for sub-column clouds generators  !
!      ntcw            : =0 no cloud condensate calculated              !
!                        >0 array index location for cloud condensate   !
!      ncld            : only used when ntcw .gt. 0                     !
!      ntoz            : =0 climatological ozone profile                !
!                        >0 interactive ozone profile                   !
!      ntrac           : dimension veriable for array oz                !
!      nfxr            : second dimension of input/output array fluxr   !
!      dtlw, dtsw      : time duration for lw/sw radiation call in sec  !
!      lsswr, lslwr    : logical flags for sw/lw radiation calls        !
!      lssav           : logical flag for store 3-d cloud field         !
!      ix,im           : horizontal dimention and num of used points    !
!      lm              : vertical layer dimension                       !
!      me              : control flag for parallel process              !
!      lprnt           : control flag for diagnostic print out          !
!      ipt             : index for diagnostic printout point            !
!      kdt             : time-step number                               !
!      deltaq          : half width of uniform total water distribution !
!      sup             : supersaturation in pdf cloud when t is very low!  
!      cnvw            : layer convective cloud water                   !
!      cnvc            : layer convective cloud cover                   !
!                                                                       !
!    output variables:                                                  !
!      htrsw (ix,lm)   : total sky sw heating rate in k/sec             !
!      topfsw(im)      : sw radiation fluxes at toa, components:        !
!                      (check module_radsw_parameters for definition)   !
!       %upfxc           - total sky upward sw flux at toa (w/m**2)     !
!       %dnflx           - total sky downward sw flux at toa (w/m**2)   !
!       %upfx0           - clear sky upward sw flux at toa (w/m**2)     !
!      sfcfsw(im)      : sw radiation fluxes at sfc, components:        !
!                      (check module_radsw_parameters for definition)   !
!       %upfxc           - total sky upward sw flux at sfc (w/m**2)     !
!       %dnfxc           - total sky downward sw flux at sfc (w/m**2)   !
!       %upfx0           - clear sky upward sw flux at sfc (w/m**2)     !
!       %dnfx0           - clear sky downward sw flux at sfc (w/m**2)   !
!      sfalb (im)      : mean surface diffused sw albedo                !
!      coszen(im)      : mean cos of zenith angle over rad call period  !
!      coszdg(im)      : daytime mean cosz over rad call period         !
!      htrlw (ix,lm)   : total sky lw heating rate in k/sec             !
!      topflw(im)      : lw radiation fluxes at top, component:         !
!                        (check module_radlw_paramters for definition)  !
!       %upfxc           - total sky upward lw flux at toa (w/m**2)     !
!       %upfx0           - clear sky upward lw flux at toa (w/m**2)     !
!      sfcflw(im)      : lw radiation fluxes at sfc, component:         !
!                        (check module_radlw_paramters for definition)  !
!       %upfxc           - total sky upward lw flux at sfc (w/m**2)     !
!       %upfx0           - clear sky upward lw flux at sfc (w/m**2)     !
!       %dnfxc           - total sky downward lw flux at sfc (w/m**2)   !
!       %dnfx0           - clear sky downward lw flux at sfc (w/m**2)   !
!      semis (im)      : surface lw emissivity in fraction              !
!      cldcov(ix,lm)   : 3-d cloud fraction                             !
!      tsflw (im)      : surface air temp during lw calculation in k    !
!                                                                       !
!    input and output variables:                                        !
!      fluxr (ix,nfxr) : to save time accumulated 2-d fields defined as:!
!                 1      - toa total sky upwd lw radiation flux         !
!                 2      - toa total sky upwd sw radiation flux         !
!                 3      - sfc total sky upwd sw radiation flux         !
!                 4      - sfc total sky dnwd sw radiation flux         !
!                 5      - high domain cloud fraction                   !
!                 6      - mid  domain cloud fraction                   !
!                 7      - low  domain cloud fraction                   !
!                 8      - high domain mean cloud top pressure          !
!                 9      - mid  domain mean cloud top pressure          !
!                10      - low  domain mean cloud top pressure          !
!                11      - high domain mean cloud base pressure         !
!                12      - mid  domain mean cloud base pressure         !
!                13      - low  domain mean cloud base pressure         !
!                14      - high domain mean cloud top temperature       !
!                15      - mid  domain mean cloud top temperature       !
!                16      - low  domain mean cloud top temperature       !
!                17      - total cloud fraction                         !
!                18      - boundary layer domain cloud fraction         !
!                19      - sfc total sky dnwd lw radiation flux         !
!                20      - sfc total sky upwd lw radiation flux         !
!                21      - sfc total sky dnwd sw uv-b radiation flux    !
!                22      - sfc clear sky dnwd sw uv-b radiation flux    !
!                23      - toa incoming solar radiation flux            !
!                24      - sfc vis beam dnwd sw radiation flux          !
!                25      - sfc vis diff dnwd sw radiation flux          !
!                26      - sfc nir beam dnwd sw radiation flux          !
!                27      - sfc nir diff dnwd sw radiation flux          !
!                28      - toa clear sky upwd lw radiation flux         !
!                29      - toa clear sky upwd sw radiation flux         !
!                30      - sfc clear sky dnwd lw radiation flux         !
!                31      - sfc clear sky upwd sw radiation flux         !
!                32      - sfc clear sky dnwd sw radiation flux         !
!                33      - sfc clear sky upwd lw radiation flux         !
!optional        34      - aeros opt depth at 550nm (all components)    !
!               ....     - optional for test and future use             !
!                                                                       !
!    optional output variables:                                         !
!      htrswb(ix,lm,nbdsw) : spectral band total sky sw heating rate    !
!      htrlwb(ix,lm,nbdlw) : spectral band total sky lw heating rate    !
!                                                                       !
!                                                                       !
!    definitions of internal variable arrays:                           !
!                                                                       !
!     1. fixed gases:         (defined in 'module_radiation_gases')     !
!          gasvmr(:,:,1)  -  co2 volume mixing ratio                    !
!          gasvmr(:,:,2)  -  n2o volume mixing ratio                    !
!          gasvmr(:,:,3)  -  ch4 volume mixing ratio                    !
!          gasvmr(:,:,4)  -  o2  volume mixing ratio                    !
!          gasvmr(:,:,5)  -  co  volume mixing ratio                    !
!          gasvmr(:,:,6)  -  cf11 volume mixing ratio                   !
!          gasvmr(:,:,7)  -  cf12 volume mixing ratio                   !
!          gasvmr(:,:,8)  -  cf22 volume mixing ratio                   !
!          gasvmr(:,:,9)  -  ccl4 volume mixing ratio                   !
!                                                                       !
!     2. cloud profiles:      (defined in 'module_radiation_clouds')    !
!                ---  for  prognostic cloud  ---                        !
!          clouds(:,:,1)  -  layer total cloud fraction                 !
!          clouds(:,:,2)  -  layer cloud liq water path                 !
!          clouds(:,:,3)  -  mean effective radius for liquid cloud     !
!          clouds(:,:,4)  -  layer cloud ice water path                 !
!          clouds(:,:,5)  -  mean effective radius for ice cloud        !
!          clouds(:,:,6)  -  layer rain drop water path                 !
!          clouds(:,:,7)  -  mean effective radius for rain drop        !
!          clouds(:,:,8)  -  layer snow flake water path                !
!          clouds(:,:,9)  -  mean effective radius for snow flake       !
!                ---  for  diagnostic cloud  ---                        !
!          clouds(:,:,1)  -  layer total cloud fraction                 !
!          clouds(:,:,2)  -  layer cloud optical depth                  !
!          clouds(:,:,3)  -  layer cloud single scattering albedo       !
!          clouds(:,:,4)  -  layer cloud asymmetry factor               !
!                                                                       !
!     3. surface albedo:      (defined in 'module_radiation_surface')   !
!          sfcalb( :,1 )  -  near ir direct beam albedo                 !
!          sfcalb( :,2 )  -  near ir diffused albedo                    !
!          sfcalb( :,3 )  -  uv+vis direct beam albedo                  !
!          sfcalb( :,4 )  -  uv+vis diffused albedo                     !
!                                                                       !
!     4. sw aerosol profiles: (defined in 'module_radiation_aerosols')  !
!          faersw(:,:,:,1)-  sw aerosols optical depth                  !
!          faersw(:,:,:,2)-  sw aerosols single scattering albedo       !
!          faersw(:,:,:,3)-  sw aerosols asymmetry parameter            !
!                                                                       !
!     5. lw aerosol profiles: (defined in 'module_radiation_aerosols')  !
!          faerlw(:,:,:,1)-  lw aerosols optical depth                  !
!          faerlw(:,:,:,2)-  lw aerosols single scattering albedo       !
!          faerlw(:,:,:,3)-  lw aerosols asymmetry parameter            !
!                                                                       !
!     6. sw fluxes at toa:    (defined in 'module_radsw_main')          !
!        (topfsw_type -- derived data type for toa rad fluxes)          !
!          topfsw(:)%upfxc  -  total sky upward flux at toa             !
!          topfsw(:)%dnfxc  -  total sky downward flux at toa           !
!          topfsw(:)%upfx0  -  clear sky upward flux at toa             !
!                                                                       !
!     7. lw fluxes at toa:    (defined in 'module_radlw_main')          !
!        (topflw_type -- derived data type for toa rad fluxes)          !
!          topflw(:)%upfxc  -  total sky upward flux at toa             !
!          topflw(:)%upfx0  -  clear sky upward flux at toa             !
!                                                                       !
!     8. sw fluxes at sfc:    (defined in 'module_radsw_main')          !
!        (sfcfsw_type -- derived data type for sfc rad fluxes)          !
!          sfcfsw(:)%upfxc  -  total sky upward flux at sfc             !
!          sfcfsw(:)%dnfxc  -  total sky downward flux at sfc           !
!          sfcfsw(:)%upfx0  -  clear sky upward flux at sfc             !
!          sfcfsw(:)%dnfx0  -  clear sky downward flux at sfc           !
!                                                                       !
!     9. lw fluxes at sfc:    (defined in 'module_radlw_main')          !
!        (sfcflw_type -- derived data type for sfc rad fluxes)          !
!          sfcflw(:)%upfxc  -  total sky upward flux at sfc             !
!          sfcflw(:)%dnfxc  -  total sky downward flux at sfc           !
!          sfcflw(:)%dnfx0  -  clear sky downward flux at sfc           !
!                                                                       !
!! optional radiation outputs:                                          !
!!   10. sw flux profiles:    (defined in 'module_radsw_main')          !
!!       (profsw_type -- derived data type for rad vertical profiles)   !
!!         fswprf(:,:)%upfxc - total sky upward flux                    !
!!         fswprf(:,:)%dnfxc - total sky downward flux                  !
!!         fswprf(:,:)%upfx0 - clear sky upward flux                    !
!!         fswprf(:,:)%dnfx0 - clear sky downward flux                  !
!!                                                                      !
!!   11. lw flux profiles:    (defined in 'module_radlw_main')          !
!!       (proflw_type -- derived data type for rad vertical profiles)   !
!!         flwprf(:,:)%upfxc - total sky upward flux                    !
!!         flwprf(:,:)%dnfxc - total sky downward flux                  !
!!         flwprf(:,:)%upfx0 - clear sky upward flux                    !
!!         flwprf(:,:)%dnfx0 - clear sky downward flux                  !
!!                                                                      !
!!   12. sw sfc components:   (defined in 'module_radsw_main')          !
!!       (cmpfsw_type -- derived data type for component sfc fluxes)    !
!!         scmpsw(:)%uvbfc  -  total sky downward uv-b flux at sfc      !
!!         scmpsw(:)%uvbf0  -  clear sky downward uv-b flux at sfc      !
!!         scmpsw(:)%nirbm  -  total sky sfc downward nir direct flux   !
!!         scmpsw(:)%nirdf  -  total sky sfc downward nir diffused flux !
!!         scmpsw(:)%visbm  -  total sky sfc downward uv+vis direct flx !
!!         scmpsw(:)%visdf  -  total sky sfc downward uv+vis diff flux  !
!                                                                       !
!    external module variables:                                         !
!     ivflip           : control flag for in/out vertical indexing      !
!                        =0 index from toa to surface                   !
!                        =1 index from surface to toa                   !
!     icmphys          : cloud microphysics scheme control flag         !
!                        =1 zhao/carr/sundqvist microphysics scheme     !
!                        =2 brad ferrier microphysics scheme            !
!                        =3 zhao/carr/sundqvist microphysics +pdf cloud !
!                                                                       !
!    module variables:                                                  !
!     itsfc            : =0 use same sfc skin-air/ground temp           !
!                        =1 use diff sfc skin-air/ground temp (not yet) !
!                                                                       !
!  ======================  end of definations  =======================  !
!
      implicit none


!  ---  inputs: (for rank>1 arrays, horizontal dimensioned by ix)
      integer,  intent(in) :: ix,im, lm, ntrac, nfxr, me,          &
     &                        ntoz, ntcw, ncld, ipt, kdt, myrank,  &
                              ntiw, ntrw, ntsw, ntgl
      integer,  intent(in) :: icsdsw(im), icsdlw(im), jdate(8)

      logical,  intent(in) :: lsswr, lslwr, lssav, lprnt

      real (kind=kind_phys), dimension(ix,lm+1), intent(in) ::  prsi

      real (kind=kind_phys), dimension(ix,lm),   intent(in) ::  prsl,   &
!            prslk, tgrs, qgrs, vvl, fcice, frain, rrime, deltaq, cnvw, & 
!            cnvc
             prslk, tgrs, qgrs, vvl, deltaq, cnvw, cnvc
!     real (kind=kind_phys), dimension(im), intent(in) :: flgmin
      real(kind=kind_phys), intent(in) ::sup

      real (kind=kind_phys), dimension(im),      intent(in) ::  slmsk,  &
             xlon, xlat, tsfc, snowd, zorl, hprim, alvsf, alnsf, alvwf, &
             alnwf, facsf, facwf, cv, cvt, cvb, fice, tisfc,            &
             sncovr, snoalb, sinlat, coslat

      real (kind=kind_phys), intent(in) :: solcon, dtlw, dtsw, solhr,   &
             tracer(ix,lm,ntrac)

!  ---  outputs: (horizontal dimensioned by ix)
      real (kind=kind_phys), dimension(ix,lm),intent(out):: htrsw,htrlw

      real (kind=kind_phys), dimension(im),   intent(out):: tsflw,      &
             sfalb, semis, coszen, coszdg

! --- cmy
      real (kind=kind_phys), dimension(ix,lm+1+ltp)::                   &
             fusl,fdsl,fuir,fdir,fuslr,fdslr,fuirr,fdirr
! --- cmy

!     type (topfsw_type), dimension(im), intent(out) :: topfsw
!     type (sfcfsw_type), dimension(im), intent(out) :: sfcfsw
! --- cmy
      type (topfsw_type), dimension(im) :: topfsw
      type (sfcfsw_type), dimension(im) :: sfcfsw
! --- cmy

!     type (topflw_type), dimension(im), intent(out) :: topflw
!     type (sfcflw_type), dimension(im), intent(out) :: sfcflw
! --- cmy
      type (topflw_type), dimension(im) :: topflw
      type (sfcflw_type), dimension(im) :: sfcflw
! --- cmy

!  ---  variables are for both input and output:
      real (kind=kind_phys), intent(inout) :: cldcov(im,lm+ltp)
      real (kind=kind_phys), intent(out) :: fluxr(ix,nfxr)

!! ---  optional outputs:
!     real (kind=kind_phys), dimension(ix,lm,nbdsw), optional,          &
!    &                       intent(out) :: htrswb
!     real (kind=kind_phys), dimension(ix,lm,nbdlw), optional,          &
!    &                       intent(out) :: htrlwb
!     real (kind=kind_phys), dimension(ix,lm), optional,                &
!    &                       intent(out) :: htrlw0
!     real (kind=kind_phys), dimension(ix,lm), optional,                &
!    &                       intent(out) :: htrsw0

      real (kind=kind_phys), dimension(ix,lm,nbdsw) :: htrswb
      real (kind=kind_phys), dimension(ix,lm,nbdlw) :: htrlwb
      real (kind=kind_phys), dimension(ix,lm) :: htrlw0
      real (kind=kind_phys), dimension(ix,lm) :: htrsw0

!  ---  local variables: (horizontal dimensioned by im)
      real (kind=kind_phys), dimension(im,lm+1+ltp):: plvl, tlvl

      real (kind=kind_phys), dimension(im,lm+ltp)  :: plyr, tlyr, qlyr, &
             olyr, rhly, qstl, vvel, clw, prslk1, tem2da, tem2db, tvly
      real (kind=kind_phys), dimension(im,lm+ltp)  :: qst2, rhly2
      real (kind=kind_phys), dimension(im,lm+ltp)  :: es2, qs2
      real (kind=kind_phys), dimension(im,lm+ltp)  :: qa
      real (kind=kind_phys), dimension(im,lm+ltp)  :: cnvw1, cnvc1

      real (kind=kind_phys), dimension(im) :: tsfa, cvt1, cvb1, tem1d,  &
             sfcemis, tsfg, tskn

      real (kind=kind_phys), dimension(im,lm+ltp,nf_clds) :: clouds
      real (kind=kind_phys), dimension(im,lm+ltp,nf_vgas) :: gasvmr
      real (kind=kind_phys), dimension(im,       nf_albd) :: sfcalb
!     real (kind=kind_phys), dimension(im,       nspc1)   :: aerodp      ! optn for aod output
      real (kind=kind_phys), dimension(im,lm+ltp,ntrac)   :: tracer1

      real (kind=kind_phys), dimension(im,lm+ltp,nbdsw,nf_aesw)::faersw
      real (kind=kind_phys), dimension(im,lm+ltp,nbdlw,nf_aelw)::faerlw

      real (kind=kind_phys), dimension(im,lm+ltp) :: htswc
      real (kind=kind_phys), dimension(im,lm+ltp) :: htlwc

      real (kind=kind_phys), dimension(im,lm+ltp) :: gcice, grain, grime

!! ---  may be used for optional sw/lw outputs:
!!      take out "!!" as needed
      real (kind=kind_phys), dimension(im,lm+ltp)   :: htsw0
      type (profsw_type),    dimension(im,lm+1+ltp) :: fswprf
      type (cmpfsw_type),    dimension(im)          :: scmpsw
      real (kind=kind_phys), dimension(im,lm+ltp,nbdsw) :: htswb

      real (kind=kind_phys), dimension(im,lm+ltp)   :: htlw0
      type (proflw_type),    dimension(im,lm+1+ltp) :: flwprf
      real (kind=kind_phys), dimension(im,lm+ltp,nbdlw) :: htlwb

      real (kind=kind_phys) :: raddt, es, qs, tem0d, cldsa(im,5),qss

      integer :: i, j, k, k1, lv, itop, ibtc, nday, idxday(im),         &
             mbota(im,3), mtopa(im,3), lp1, nb, lmk, lmp, kd, lla, llb, &
             lya, lyb, kt, kb
!effective radius for liquid, ice, snow, rain
      real (kind=kind_phys), dimension(im,lm+ltp,5)   :: phy_f3d
      logical uni_cloud,lmfshal,lmfdeep2

!  ---  for debug test use
!     real (kind=kind_phys) :: temlon, temlat, alon, alat
!     integer :: ipt
!     logical :: lprnt1
!
!! ---  logical flags for optional output fields

!     logical :: lhtrswb  = .false.
!     logical :: lhtrsw0  = .false.
!     logical :: lfswprf  = .false.
!     logical :: lscmpsw  = .false.
!
!     logical :: lhtrlwb  = .false.
!     logical :: lhtrlw0  = .false.
!     logical :: lflwprf  = .false.

!
!===> ...  begin here
!
!     lhtrswb  = present( htrswb )
!     lhtrsw0  = present( htrsw0 )
!     lfswprf  = present( fswprf )
!     lscmpsw  = present( scmpsw )
!
!     lhtrlwb  = present( htrlwb )
!     lhtrlw0  = present( htrlw0 )
!     lflwprf  = present( flwprf )
       
!     if (myrank .eq. 0) print *,' #### present (htrswb)=',lhtrswb
!     if (myrank .eq. 0) print *,' #### present (htrsw0)=',lhtrsw0
!     if (myrank .eq. 0) print *,' #### present (fswprf)=',lfswprf
!     if (myrank .eq. 0) print *,' #### present (scmpsw)=',lscmpsw
!     if (myrank .eq. 0) print *,' #### present (htrlwb)=',lhtrlwb
!     if (myrank .eq. 0) print *,' #### present (htrlw0)=',lhtrlw0
!     if (myrank .eq. 0) print *,' #### present (flwprf)=',lflwprf

!
      lp1 = lm + 1               ! num of in/out levels

!  --- ...  set local /level/layer indexes corresponding to in/out variables

      lmk = lm + ltp             ! num of local layers
      lmp = lmk + 1              ! num of local levels

      if ( lextop ) then
        if ( ivflip == 1 ) then    ! vertical from sfc upward
          kd = 0                   ! index diff between in/out and local
          kt = 1                   ! index diff between lyr and upper bound
          kb = 0                   ! index diff between lyr and lower bound
          lla = lmk                ! local index at the 2nd level from top
          llb = lmp                ! local index at toa level
          lya = lm                 ! local index for the 2nd layer from top
          lyb = lp1                ! local index for the top layer
        else                       ! vertical from toa downward
          kd = 1                   ! index diff between in/out and local
          kt = 0                   ! index diff between lyr and upper bound
          kb = 1                   ! index diff between lyr and lower bound
          lla = 2                  ! local index at the 2nd level from top
          llb = 1                  ! local index at toa level
          lya = 2                  ! local index for the 2nd layer from top
          lyb = 1                  ! local index for the top layer
        endif                    ! end if_ivflip_block
      else
        kd = 0
        if ( ivflip == 1 ) then  ! vertical from sfc upward
          kt = 1                   ! index diff between lyr and upper bound
          kb = 0                   ! index diff between lyr and lower bound
        else                     ! vertical from toa downward
          kt = 0                   ! index diff between lyr and upper bound
          kb = 1                   ! index diff between lyr and lower bound
        endif                    ! end if_ivflip_block
      endif   ! end if_lextop_block

      raddt = min(dtsw, dtlw)

! ---------------------------------------------------------------------
        if ( me == 0 .and. myrank == 0) then
            print *,'###################################################' 
            print *,'### In grrad start !! ###' 
            print *,'###################################################' 
            print *,'### ix=',ix,' im=',im,' lm=',lm,' me=',me
            print *,'### ipt=',ipt
            print *,'### iter=',kdt
            print *,'### solhr=',solhr
            print *,'###################################################' 
            print *,'### ncld=',ncld
            print *,'### ntoz=',ntoz
            print *,'### ntcw=',ntcw
            print *,'### ntrac=',ntrac
            print *,'### nfxr=',nfxr
            print *,'### dtsw=',dtsw
            print *,'### dtlw=',dtlw
            print *,'###################################################' 
            print *,'### lsswr=',lsswr
            print *,'### lslwr=',lslwr
            print *,'### lssav=',lssav
            print *,'### lprnt=',lprnt
            print *,'###################################################' 
            print *,'### solcon=',solcon
            print *,'###################################################' 
            print *,'### xlon(ipt)=',xlon(ipt)
            print *,'### xlat(ipt)=',xlat(ipt)
            print *,'### sinlat(ipt)=',sinlat(ipt)
            print *,'### coslat(ipt)=',coslat(ipt)
            print *,'###################################################' 
            print *,'### jdate(1-4)=',jdate(1),jdate(2),jdate(3),jdate(4) 
            print *,'### jdate(5-8)=',jdate(5),jdate(6),jdate(7),jdate(8)
            print *,'###################################################' 
            print *,'### cv(ipt)=',cv(ipt)
            print *,'### cvb(ipt)=',cvb(ipt)
            print *,'### cvt(ipt)=',cvt(ipt)
            print *,'###################################################' 
            print *,'### prsi(ipt,lm+1)=',prsi(ipt,lm+1)
            print *,'### prsl(ipt,lm)=',prsl(ipt,lm)
            print *,'### prslk(ipt,lm)=',prslk(ipt,lm)
            print *,'### vvl(ipt,lm)=',vvl(ipt,lm)
            print *,'###################################################' 
            print *,'### tsfc(ipt)=',tsfc(ipt)
            print *,'### tgrs(ipt,lm)=',tgrs(ipt,lm)
            print *,'### qgrs(ipt,lm)=',qgrs(ipt,lm)
            print *,'##############################################'
            print *,'### tracer(ipt,lm,3) =',tracer(ipt,lm,3)
            print *,'###################################################' 
            print *,'### slmsk(ipt)=',slmsk(ipt)
            print *,'### fice(ipt)=',fice(ipt)
            print *,'### tisfc(ipt)=',tisfc(ipt)
            print *,'###################################################' 
            print *,'### snowd(ipt)=',snowd(ipt)
            print *,'### sncovr(ipt)=',sncovr(ipt)
            print *,'### snoalb(ipt)=',snoalb(ipt)
            print *,'###################################################' 
            print *,'### zorl(ipt)=',zorl(ipt)
            print *,'### hprim(ipt)=',hprim(ipt)
            print *,'###################################################' 
            print *,'### alvsf(ipt)=',alvsf(ipt)
            print *,'### alnsf(ipt)=',alnsf(ipt)
            print *,'### alvwf(ipt)=',alvwf(ipt)
            print *,'### alnwf(ipt)=',alnwf(ipt)
            print *,'### facsf(ipt)=',facsf(ipt)
            print *,'### facwf(ipt)=',facwf(ipt)
            print *,'###################################################' 
            print *,'### icsdsw(ipt)=',icsdsw(ipt)
            print *,'### icsdlw(ipt)=',icsdlw(ipt)
            print *,'###################################################' 
          endif
! ---------------------------------------------------------------------
!  --- ...  for debug test
!     alon = 120.0
!     alat = 29.5
!     ipt = 0
!     do i = 1, im
!       temlon = xlon(i) * 57.29578
!       if (temlon < 0.0) temlon = temlon + 360.0
!       temlat = xlat(i) * 57.29578
!       lprnt1 = abs(temlon-alon) < 1.1 .and. abs(temlat-alat) < 1.1
!       if ( lprnt1 ) then
!         ipt = i
!         exit
!       endif
!     enddo

      if (me == 0 .and. myrank ==0) print *,'### raddt=',raddt
      if (me == 0 .and. myrank ==0) print *,'### itsfc=',itsfc

!  --- ...  setup surface ground temp and ground/air skin temp if required

      if ( itsfc == 0 ) then            ! use same sfc skin-air/ground temp
           
        do i = 1, im
          tskn(i) = tsfc(i)
          tsfg(i) = tsfc(i)
        enddo
      
      else                              ! use diff sfc skin-air/ground temp
        do i = 1, im
!!        tskn(i) = ta  (i)               ! not yet
!!        tsfg(i) = tg  (i)               ! not yet
          tskn(i) = tsfc(i)
          tsfg(i) = tsfc(i)
        enddo
      endif
       if (me == 0 .and. myrank ==0) print *,'### tskn(ipt)=',tskn(ipt),' ipt=',ipt
       if (me == 0 .and. myrank ==0) print *,'### tsfg(ipt)=',tsfg(ipt),' ipt=',ipt

!  --- ...  prepare atmospheric profiles for radiation input
!           convert pressure unit from cb to mb

      do k = 1, lm
        k1 = k + kd

        do i = 1, im
          plvl(i,k1)   = 10.0 * prsi(i,k)   ! cb (kpa) to mb (hpa)
          plyr(i,k1)   = 10.0 * prsl(i,k)   ! cb (kpa) to mb (hpa)
!         plvl(i,k1)   = 0.01 * prsi(i,k)   ! pa to mb (hpa)
!         plyr(i,k1)   = 0.01 * prsl(i,k)   ! pa to mb (hpa)
          tlyr(i,k1)   = tgrs(i,k)
          prslk1(i,k1) = prslk(i,k)
          cnvw1(i,k1)  = cnvw(i,k)
          cnvc1(i,k1)  = cnvc(i,k)

!  --- ...  compute relative humidity
!         es  = min( prsl(i,k), 0.001 * fpvs( tgrs(i,k) ) )   ! fpvs in pa
!         qs  = max( qmin, con_eps * es / (prsl(i,k) + con_epsm1*es) )
!         rhly(i,k1) = max( 0.0, min( 1.0, max(qmin, qgrs(i,k))/qs ) )
!         qstl(i,k1) = qs
!--------------------------------------------------------------------------
          qlyr(i,k1) = max( qme6, qgrs(i,k) )
          call qsatq(1,tlyr(i,k1),plyr(i,k1),qss) !plyr in mb
          rhly(i,k1)= max( 0.0, min( 1.0, max(qmin, qlyr(i,k1))/qss ) )
          qstl(i,k1) = qss
        enddo
!---------------------------------------------------------------------------
        do j = 1, ntrac
          do i = 1, im
             tracer1(i,k1,j) = tracer(i,k,j)
          enddo
        enddo
      enddo

        if ( me == 0 .and. myrank == 0) then
          print *,'###################################################' 
          print *,'### prsl(ipt,1) =',prsl(ipt,1)*10.,' in mb ###'
          print *,'### tgrs(ipt,1) =',tgrs(ipt,1)
          print *,'### qgrs(ipt,1) =',qgrs(ipt,1)
          print *,'### qlyr(ipt,1) =',qlyr(ipt,1)
          print *,'### rhly(ipt,1) =',rhly(ipt,1)
          print *,'### qstl(ipt,1) =',qstl(ipt,1)
          print *,'###################################################' 
          print *,'### prsl(ipt,lm) =',prsl(ipt,lm)*10.,' in mb ###'
          print *,'### tgrs(ipt,lm) =',tgrs(ipt,lm)
          print *,'### qgrs(ipt,lm) =',qgrs(ipt,lm)
          print *,'### qlyr(ipt,lm) =',qlyr(ipt,lm)
          print *,'### rhly(ipt,lm) =',rhly(ipt,lm)
          print *,'### qstl(ipt,lm) =',qstl(ipt,lm)
          print *,'###################################################' 
        endif


      do i = 1, im
        plvl(i,lp1+kd) = 10.0 * prsi(i,lp1)  ! cb (kpa) to mb (hpa)
!       plvl(i,lp1+kd) = 0.01 * prsi(i,lp1)  ! pa to mb (hpa)
      enddo

      if ( lextop ) then                 ! values for extra top layer
        do i = 1, im
          plvl(i,llb) = prsmin
          if ( plvl(i,lla) <= prsmin ) plvl(i,lla) = 2.0*prsmin
          plyr(i,lyb)   = 0.5 * plvl(i,lla)
          tlyr(i,lyb)   = tlyr(i,lya)
          prslk1(i,lyb) = (plyr(i,lyb)*0.001) ** rocp ! plyr in hpa

          rhly(i,lyb)   = rhly(i,lya)
          qstl(i,lyb)   = qstl(i,lya)
        enddo

        do j = 1, ntrac
          do i = 1, im
!  ---  note: may need to take care the top layer amount
             tracer1(i,lyb,j) = tracer1(i,lya,j)
          enddo
        enddo
      endif
!cmy--------------------------------------------------------------------
!  --- ...  extra variables needed for ferrier's microphysics
!     if (icmphys == 2) then
!       do k = 1, lm
!         k1 = k + kd
!         do i = 1, im
!           gcice(i,k1)= fcice(i,k)
!           grain(i,k1)= frain(i,k)
!           grime(i,k1)= rrime(i,k)
!         enddo
!       enddo
!       if ( lextop ) then
!         do i = 1, im
!           gcice(i,lyb) = fcice(i,lya)
!           grain(i,lyb) = frain(i,lya)
!           grime(i,lyb) = rrime(i,lya)
!         enddo
!       endif
!     endif   ! if_icmphys
!cmy--------------------------------------------------------------------

!  --- ...  get layer ozone mass mixing ratio

      if (ntoz > 0) then            ! interactive ozone generation

        do k = 1, lmk
          do i = 1, im
            olyr(i,k) = max( qmin, tracer1(i,k,ntoz) )
          enddo
        enddo
        if (me == 0 .and. myrank ==0) print *, '### ntoz=',ntoz,' ipt=',ipt
        if (me == 0 .and. myrank ==0) print *, '### olyr(i,k)>= qmin, qmin=',qmin
        if (me == 0 .and. myrank ==0) print *, '### olyr(ipt,lm)=',olyr(ipt,lm)

!     else                          ! climatological ozone

!     print *,' in grrad : calling getozn'
!       call getozn                                                     &
!  ---  inputs:
!    &     ( prslk1,xlat,                                               &
!    &       im, lmk,                                                   &
!  ---  outputs:
!    &       olyr                                                       &
!    &     )

      endif                            ! end_if_ntoz

!  --- ...  compute cosin of zenith angle

      if (me == 0 .and. myrank ==0) print *,'### call coszmn'
      call coszmn                                                       &
!  ---  inputs:
     &     ( xlon,sinlat,coslat,solhr, im, me,                          &
!  ---  outputs:
     &       coszen, coszdg                                             &
     &      )

      if ( myrank == 0 .and. me == 0) print *,'### call coszmn ok ! ###'
      if ( myrank == 0 .and. me == 0) then 
          print *,'solhr=',solhr,' ipt=',ipt                          
          print *,'xlon(ipt)=',xlon(ipt)
          print *,'sinlat(ipt)=',sinlat(ipt)
          print *,'coslat(ipt)=',coslat(ipt)
          print *,'coszen(ipt)=',coszen(ipt)
          print *,'coszdg(ipt)=',coszdg(ipt)
      endif
!
!  --- ...  set up non-prognostic gas volume mixing ratioes

      if (me == 0 .and. myrank ==0) print *,'### call getgases'
      call getgases                                                     &
!  ---  inputs:
     &    ( plvl, xlon, xlat,                                           &
     &      im, lmk,                                                    &
!  ---  outputs:
     &      gasvmr                                                      &
     &     )
      if (me == 0 .and. myrank ==0) then
          print *,'###################################################' 
          print *,'### call getgase ok !!'
          print *,'### plvl(ipt,1)=',plvl(ipt,1)
          print *,'### plvl(ipt,2)=',plvl(ipt,2)
          print *,'### plvl(ipt,lm)=',plvl(ipt,lm)
          print *,'### plvl(ipt,lm+1)=',plvl(ipt,lm+1)
          print *,'### xlon(ipt)=',xlon(ipt)
          print *,'### xlat(ipt)=',xlat(ipt)
          print *,'### lmk=',lmk,' im=',im,' ltp=',ltp
          print *,'### gasvmr(ipt,lm,1)_co2 =',gasvmr(ipt,lm,1)
          print *,'### gasvmr(ipt,lm,2)_n2o =',gasvmr(ipt,lm,2)
          print *,'### gasvmr(ipt,lm,3)_ch4 =',gasvmr(ipt,lm,3)
          print *,'### gasvmr(ipt,lm,4)_o2  =',gasvmr(ipt,lm,4)
          print *,'### gasvmr(ipt,lm,5)_co  =',gasvmr(ipt,lm,5)
          print *,'### gasvmr(ipt,lm,6)_cf11=',gasvmr(ipt,lm,6)
          print *,'### gasvmr(ipt,lm,7)_cf12=',gasvmr(ipt,lm,7)
          print *,'### gasvmr(ipt,lm,8)_cf22=',gasvmr(ipt,lm,8)
          print *,'### gasvmr(ipt,lm,9)_ccl4=',gasvmr(ipt,lm,9)
          print *,'###################################################' 
       endif  
!
!  --- ...  get temperature at layer interface, and layer moisture

      do k = 2, lmk
        do i = 1, im
          tem2da(i,k) = log( plyr(i,k) )
          tem2db(i,k) = log( plvl(i,k) )
        enddo
      enddo


      if (ivflip == 0) then              ! input data from toa to sfc

        do i = 1, im
          tem1d (i)   = qme6
          tem2da(i,1) = log( plyr(i,1) )
          tem2db(i,1) = 1.0
          tsfa  (i)   = tlyr(i,lmk)                  ! sfc layer air temp
          tlvl(i,1)   = tlyr(i,1)
          tlvl(i,lmp) = tskn(i)
        enddo

        do k = 1, lm
          k1 = k + kd

          do i = 1, im
            qlyr(i,k1) = max( tem1d(i), qgrs(i,k) )
            tem1d(i)   = min( qme5, qlyr(i,k1) )
            tvly(i,k1) = tgrs(i,k) * (1.0 + con_fvirt*qlyr(i,k1))! virtual temp in k
          enddo
        enddo

        if ( lextop ) then
          do i = 1, im
            qlyr(i,lyb) = qlyr(i,lya)
            tvly(i,lyb) = tvly(i,lya)
          enddo
        endif

        do k = 2, lmk
          do i = 1, im
            tlvl(i,k) = tlyr(i,k) + (tlyr(i,k-1) - tlyr(i,k))           &
     &                * (tem2db(i,k)   - tem2da(i,k))                   &
     &                / (tem2da(i,k-1) - tem2da(i,k))
          enddo
        enddo

      else                               ! input data from sfc to toa

        do i = 1, im
          tem1d (i)   = qme6
          tem2da(i,1) = log( plyr(i,1) )
          tem2db(i,1) = log( plvl(i,1) )
          tsfa  (i)   = tlyr(i,1)                    ! sfc layer air temp
          tlvl(i,1)   = tskn(i)
          tlvl(i,lmp) = tlyr(i,lmk)
        enddo

        do k = lm, 1, -1
          do i = 1, im
            qlyr(i,k) = max( tem1d(i), qgrs(i,k) )
            tem1d(i)  = min( qme5, qlyr(i,k) )
            tvly(i,k) = tgrs(i,k) * (1.0 + con_fvirt*qlyr(i,k)) ! virtual temp in k
          enddo
        enddo

        if ( lextop ) then
          do i = 1, im
            qlyr(i,lyb) = qlyr(i,lya)
            tvly(i,lyb) = tvly(i,lya)
          enddo
        endif

        do k = 1, lmk-1
          do i = 1, im
            tlvl(i,k+1) = tlyr(i,k) + (tlyr(i,k+1) - tlyr(i,k))         &
     &                  * (tem2db(i,k+1) - tem2da(i,k))                 &
     &                  / (tem2da(i,k+1) - tem2da(i,k))
          enddo
        enddo

      endif                              ! end_if_ivflip
      fusl       = 0.
      fdsl       = 0.
      fuslr      = 0.
      fdslr      = 0.

!  --- ...  check for daytime points

      nday = 0
      do i = 1, im
        if (coszen(i) >= 0.0001) then
          nday = nday + 1
          idxday(nday) = i
        endif
      enddo
!      if (myrank == 0 ) print *,'nday=',nday

!  --- ...  setup aerosols property profile for radiation

      if (me == 0 .and. myrank ==0) print *,'### before setaer ###'
      if (myrank == 0 .and. me == 0)then
          print *,'###################################################' 
          print *,'### plvl(ipt,lm+1)=',plvl(ipt,lm+1)
          print *,'### plyr(ipt,lm)=',plyr(ipt,lm)
          print *,'### prslk1(ipt,lm)=',prslk1(ipt,lm)
          print *,'### tvly(ipt,lm)=',tvly(ipt,lm)
          print *,'### tlyr(ipt,lm)=',tlyr(ipt,lm)
          print *,'### qlyr(ipt,lm)=',qlyr(ipt,lm)
          print *,'### rhly(ipt,lm)=',rhly(ipt,lm)
          print *,'###################################################' 
      endif


      call setaer                                                       &
!  ---  inputs:
     &     ( plvl,plyr,prslk1,tvly,rhly,slmsk,tracer1,xlon,xlat,        &
     &       im,lmk,lmp,lsswr,lslwr,me,myrank,                          &
!  ---  outputs:
     &       faersw,faerlw                                              &
!    &       faersw,faerlw,aerodp                                       &
     &     )

       if (me == 0 .and. myrank ==0) then
           print *,'###################################################' 
           print *,'### after call setaer ###'
           print *,'###################################################' 
           print *,'### im=',im,' lmk=',lmk,' lmp=',lmp,' ipt=',ipt
           print *,'### xlon(ipt)=', xlon(ipt)
           print *,'### xlat(ipt)=', xlat(ipt)
           print *,'### slmsk(ipt)=',slmsk(ipt)
           print *,'### plvl(ipt,lm+1)=',plvl(ipt,lm+1)
           print *,'### plyr(ipt,lm)=',plyr(ipt,lm)
           print *,'### prslk1(ipt,lm)=',prslk1(ipt,lm)
           print *,'### tvly(ipt,lm)=',tvly(ipt,lm)
           print *,'### rhly(ipt,lm)=',rhly(ipt,lm)
           print *,'### tracer1(ipt,lm,1)=',tracer1(ipt,lm,1)
           print *,'### tracer1(ipt,lm,2)=',tracer1(ipt,lm,2)
           print *,'### tracer1(ipt,lm,3)=',tracer1(ipt,lm,3)
           print *,'###################################################' 
           print *,'### faersw(ipt,lm,1,1)=sw#1-opd =',faersw(ipt,lm,1,1)
           print *,'### faersw(ipt,lm,1,2)=sw#1-ssa =',faersw(ipt,lm,1,2)
           print *,'### faersw(ipt,lm,1,3)=sw#1-asy =',faersw(ipt,lm,1,3)
           print *,'### faerlw(ipt,lm,1,1)=lw#1-opd =',faerlw(ipt,lm,1,1)
           print *,'### faerlw(ipt,lm,1,2)=lw#1-ssa =',faerlw(ipt,lm,1,2)
           print *,'### faerlw(ipt,lm,1,3)=lw#1-asy =',faerlw(ipt,lm,1,3)
           print *,'###################################################' 
       endif
   

!  --- ...  obtain cloud information for radiation calculations

      if (ntcw > 0) then                   ! prognostic cloud scheme
!
        if (icmphys == 1) then           ! zhao/moorthi's prognostic cloud scheme
!
          do k = 1, lmk
            do i = 1, im
              clw(i,k) = 0.0
            enddo

            do j = 1, ncld
              lv = ntcw + j - 1
              do i = 1, im
!byl                 clw(i,k) = clw(i,k) + tracer1(i,k,lv)   ! cloud condensate amount
                 clw(i,k) = clw(i,k) + tracer1(i,k,lv) + cnvw(i,k)  ! cloud condensate amount
              enddo
            enddo
          enddo

          do k = 1, lmk
            do i = 1, im
              if ( clw(i,k) < epsq ) clw(i,k) = 0.0
            enddo
          enddo
 
          if ( me == 0 .and. myrank == 0 )                              &
            print *,'### call progcld1 -zhao/moorhi ###' 
          call progcld1                                                 &
!  ---  inputs:
     &     ( plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,                    &
     &       xlat,xlon,slmsk,                                           &
     &       im, lmk, lmp, myrank,                                      &
!  ---  outputs:
     &       clouds,cldsa,mtopa,mbota                                   &
     &      )


!       elseif (icmphys == 2) then       ! ferrier's microphysics

!     print *,' in grrad : calling progcld2'
!         call progcld2                                                 &
!  ---  inputs:
!    &     ( plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,                    &
!    &       xlat,xlon,slmsk, gcice,grain,grime,flgmin,                 &
!    &       im, lmk, lmp,                                              &
!  ---  outputs:
!    &       clouds,cldsa,mtopa,mbota                                   &
!    &      )

       elseif(icmphys == 3) then      ! zhao/moorthi's prognostic cloud+pdfcld
!
          do k = 1, lmk
            do i = 1, im
              clw(i,k) = 0.0
            enddo

            do j = 1, ncld
              lv = ntcw + j - 1
              do i = 1, im
                 clw(i,k) = clw(i,k) + tracer1(i,k,lv)   ! cloud condensate amount
              enddo
            enddo
          enddo

          do k = 1, lmk
            do i = 1, im
              if ( clw(i,k) < epsq ) clw(i,k) = 0.0
            enddo
          enddo
!
         if ( me == 0 .and. myrank == 0 )                               &
           print *,'### call progcld3 -zhao/moorhi with PDF cloud###' 
         call progcld3                                                  &
!  ---  inputs:
     &     ( plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,cnvw,cnvc,          &
     &       xlat,xlon,slmsk,                                           &
     &       im, lmk, lmp,                                              &
     &       deltaq, sup,kdt,me,                                        &
!  ---  outputs:
     &       clouds,cldsa,mtopa,mbota                                   &
     &      )
!
       elseif (icmphys == 6 .or. icmphys == 8) then    ! wsm6 & Thompson
         if ( me == 0 .and. myrank == 0 ) then
            if ( icmphys == 6 ) print *,'### call WSM6 cloud###' 
            if ( icmphys == 8 ) print *,'### call Thompson cloud###' 
         endif
         
         if (kdt == 1) then
           phy_f3d(:,:,1) = 10.
           phy_f3d(:,:,2) = 50.
           phy_f3d(:,:,3) = 250.
         endif
!
         call progcld4                               &
!  --- inputs
          ( plyr,plvl,tlyr,qlyr,qstl,rhly,tracer1,   &
            xlat,xlon,slmsk,                         &
            ntrac,ntcw,ntiw,ntrw,ntsw,ntgl,          &
            im, lmk, lmp,                            &
            uni_cloud,lmfshal,lmfdeep2,              &
            cldcov,phy_f3d(:,:,1),                   &
            phy_f3d(:,:,2),phy_f3d(:,:,3),           &
!   --- outputs:
            clouds,cldsa,mtopa,mbota                 &
           )
       elseif ( icmphys == 18 ) then   ! 2M Thompson
         if ( me == 0 .and. myrank == 0 )                               &
           print *,'### call New Thompson cloud ###'

         if (kdt == 1) then
           phy_f3d(:,:,1) = 10.
           phy_f3d(:,:,2) = 50.
           phy_f3d(:,:,3) = 250.
         endif

!         lwp_ex=0.0  !total liquid water path from explicit microphysics
!         iwp_ex=0.0  !total ice water path from explicit microphysics
!         lwp_fc=0.0  !total liquid water path from cloud fraction scheme
!         iwp_fc=0.0  !total ice water path from cloud fraction scheme
         call progcld_thompson                                          &
!  --- inputs
          ( plyr, plvl, tlyr, qlyr, qstl, rhly, tracer1,                &
            xlat, xlon, slmsk,                                          &
            ntrac, ntcw, ntiw, ntrw, ntsw, ntgl,                        &
            im, lmk, lmp,                                               &
            uni_cloud, lmfshal, lmfdeep2, cldcov,                       &
            phy_f3d(:,:,1), phy_f3d(:,:,2), phy_f3d(:,:,3),             &
!            lwp_ex, iwp_ex, lwp_fc, iwp_fc, dzlay,                      &
!            gridkm,                                                     &
!   --- outputs:
!            cld_frac, cld_lwp, cld_reliq, cld_iwp,                      &
!            cld_reice, cld_rwp, cld_rerain, cld_swp, cld_resnow)
            clouds, cldsa, mtopa, mbota )
       elseif ( icmphys == 11 ) then   ! GFDL MP v1
         if ( me == 0 .and. myrank == 0 )                               &
           print *,'### call GFDL cloud ###'
         clw = 0.0
         if ( .not. lgfdlmprad ) then
         do k = 1, lmk
           do i = 1, im
             do j = 1, ncld - 1
               lv = ntcw + j - 1
               clw(i,k) = clw(i,k) + tracer1(i,k,lv)  ! cloud condensate amount
             enddo
             if ( clw(i,k) < epsq ) clw(i,k) = 0.0
           enddo
         enddo
         endif

         if (kdt == 1) then
           phy_f3d(:,:,1) = 10.
           phy_f3d(:,:,2) = 50.
           phy_f3d(:,:,3) = 250.
           phy_f3d(:,:,4) = 1000.
         endif
         if ( .not. lgfdlmprad ) then  ! no consistency between GFDLMP and radiation
           call progcld5                                                &
!    ---  inputs:
             ( plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,cnvw,cnvc,        &
               xlat,xlon,slmsk,im,lmk,lmp,                              &
               cldcov,                                                  &
!    ---  outputs:
               clouds,cldsa,mtopa,mbota                                 &
              ) 
         else
           call progcld5o                                               &
!    ---  inputs:
             ( plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,tracer1,              &
               xlat,xlon,slmsk,                                         &
               ntrac,ntcw,ntiw,ntrw,ntsw,ntgl,cldcov,                   &
               phy_f3d(:,:,1),phy_f3d(:,:,2),phy_f3d(:,:,3),            &
               phy_f3d(:,:,4),effr_in,                                  &
               im,lmk,lmp,                                              &
!    ---  outputs:
               clouds,cldsa,mtopa,mbota                                 &
              ) 
!           endif
         endif
       elseif ( icmphys == 12 .or. icmphys == 13 ) then   ! GFDL MP v2 / v3
         if ( me == 0 .and. myrank == 0 .and. icmphys == 12 )           &
           print *,'### call GFDL v2 cloud ###'
         if ( me == 0 .and. myrank == 0 .and. icmphys == 13 )           &
           print *,'### call GFDL v3 cloud ###'
         qa = 0.  !aerosol mixing ratio (kg/kg)
         call progcld6                                                  &
!    ---  inputs:
             ( plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,cnvw1,cnvc1,          &
               tracer1(:,:,ntcw),tracer1(:,:,ntrw),tracer1(:,:,ntiw),   &
               tracer1(:,:,ntsw),tracer1(:,:,ntgl),                     &
               cldcov,slmsk,                                            &
               phy_f3d(:,:,1),phy_f3d(:,:,2),phy_f3d(:,:,3),            &
               phy_f3d(:,:,4),effr_in,                                  &
               xlat,xlon,im,lmk,lmp,                                    &
!    ---  outputs:
               clouds,cldsa,mtopa,mbota                                 &
              ) 
       elseif ( icmphys == 15 .or. icmphys == 16 ) then   ! Goddard (GCE)
         if (kdt == 1) then
           phy_f3d(:,:,1) = 10.
           phy_f3d(:,:,2) = 50.
           phy_f3d(:,:,3) = 250.
           phy_f3d(:,:,4) = 1000.
         endif
         if ( me == 0 .and. myrank == 0 )                               &
           print *,'### call Goddard (GCE) cloud ###'
           call progcld_gce                                             &
!    ---  inputs:
             ( plyr, plvl, tlyr, tvly, qlyr, qstl, rhly, tracer1,       &
               xlat, xlon, slmsk, ntrac,                                &
               phy_f3d(:,:,1), phy_f3d(:,:,2), phy_f3d(:,:,3),          &
               phy_f3d(:,:,4), effr_in,                                 &
               im, lmk, lmp, lmfshal, lmfdeep2,                         &
!    ---  outputs:
               clouds, cldsa, mtopa, mbota, cldcov                      &
              )
        endif                            ! end if_icmphys

      else                                 ! diagnostic cloud scheme

        do i = 1, im
          cvt1(i) = 10.0 * cvt(i)
          cvb1(i) = 10.0 * cvb(i)
        enddo

        do k = 1, lm
          k1 = k + kd

          do i = 1, im
            vvel(i,k1) = 10.0 * vvl(i,k)
          enddo
        enddo

        if ( lextop ) then
          do i = 1, im
            vvel(i,lyb) = vvel(i,lya)
          enddo
        endif

!  ---  compute diagnostic cloud related quantities

      if (me == 0 .and. myrank ==0)print *,'### call diagcld1 ###' 
        call diagcld1                                                   &
!  ---  inputs:
     &     ( plyr,plvl,tlyr,rhly,vvel,cv,cvt1,cvb1,                     &
     &       xlat,xlon,slmsk,                                           &
     &       im, lmk, lmp,                                              &
!  ---  outputs:
     &       clouds,cldsa,mtopa,mbota                                   &
     &      )

      endif                                ! end_if_ntcw

       if (me == 0 .and. myrank ==0) then
           print *,'###################################################' 
           print *,'###  after diagcld1'
           print *,'###################################################' 
           print *,'### im=',im,' lmk=',lmk,' lmp=',lmp
           print *,'### ntcw=',ntcw
           print *,'### ncld=',ncld
           print *,'###################################################' 
           print *,'### xlon(ipt)=', xlon(ipt)
           print *,'### xlat(ipt)=', xlat(ipt)
           print *,'### slmsk(ipt)=',slmsk(ipt)
           print *,'###################################################' 
           print *,'### plvl(ipt,lm)=',plvl(ipt,lm)
           print *,'### plyr(ipt,lm)=',plyr(ipt,lm)
           print *,'###################################################' 
           print *,'### tlyr(ipt,lm)=',tlyr(ipt,lm)
           print *,'### tvly(ipt,lm)=',tvly(ipt,lm)
           print *,'###################################################' 
           print *,'### qlyr(ipt,lm)=',qlyr(ipt,lm)
           print *,'### qstl(ipt,lm)=',qstl(ipt,lm)
           print *,'### rhly(ipt,lm)=',rhly(ipt,lm)
           print *,'### clw(ipt,lm) =',clw(ipt,lm)
           print *,'###################################################' 
           print *,'### tracer1(ipt,lm,1)=',tracer1(ipt,lm,1)
           print *,'### tracer1(ipt,lm,2)=',tracer1(ipt,lm,2)
           print *,'### tracer1(ipt,lm,3)=',tracer1(ipt,lm,3)
           print *,'###################################################' 
           print *,'### clouds(ipt,lm,1)-total cloud fraction =',clouds(ipt,lm,1)
           print *,'### clouds(ipt,lm,2)-liq water path       =',clouds(ipt,lm,2)
           print *,'### clouds(ipt,lm,3)-liq effective radius =',clouds(ipt,lm,3)
           print *,'### clouds(ipt,lm,4)-ice water path       =',clouds(ipt,lm,4)
           print *,'### clouds(ipt,lm,5)-ice effective radius =',clouds(ipt,lm,5)
           print *,'### clouds(ipt,lm,6)-rain water path      =',clouds(ipt,lm,6)
           print *,'### clouds(ipt,lm,7)-rain effective radius=',clouds(ipt,lm,7)
           print *,'### clouds(ipt,lm,8)-snow water path      =',clouds(ipt,lm,8)
           print *,'### clouds(ipt,lm,9)-snow effective radius=',clouds(ipt,lm,9)
           print *,'###################################################' 
           print *,'### cldsa(ipt,1)-low clouds fraction=',cldsa(ipt,1)
           print *,'### cldsa(ipt,2)-mid clouds fraction=',cldsa(ipt,2)
           print *,'### cldsa(ipt,3)-hig clouds fraction=',cldsa(ipt,3)
           print *,'### cldsa(ipt,4)-tot clouds fraction=',cldsa(ipt,4)
           print *,'### cldsa(ipt,5)-bl  clouds fraction=',cldsa(ipt,5)
           print *,'###################################################' 
           print *,'### mtopa(ipt,1)-low clouds top =',mtopa(ipt,1)
           print *,'### mtopa(ipt,2)-mid clouds top =',mtopa(ipt,2)
           print *,'### mtopa(ipt,3)-hig clouds top =',mtopa(ipt,3)
           print *,'###################################################' 
           print *,'### mbota(ipt,1)-low clouds bottom =',mbota(ipt,1)
           print *,'### mbota(ipt,2)-mid clouds bottom =',mbota(ipt,2)
           print *,'### mbota(ipt,3)-hig clouds bottom =',mbota(ipt,3)
           print *,'###################################################' 
       endif
           

!  --- ...  start radiation calculations 
!           remember to set heating rate unit to k/sec!

      if (lsswr) then

!  ---  setup surface albedo for sw radiation, incl xw (nov04) sea-ice

      if (me == 0 .and. myrank ==0)print *,'### call setalb ###' 
        call setalb                                                     &
!  ---  inputs:
     &     ( slmsk,snowd,sncovr,snoalb,zorl,coszen,tsfg,tsfa,hprim,     &
     &       alvsf,alnsf,alvwf,alnwf,facsf,facwf,fice,tisfc,            &
     &       im,                                                        &
!  ---  outputs:
     &       sfcalb                                                     &
     &     )

!  --- lu [+4l]: derive sfalb from vis- and nir- diffuse surface albedo
        do i = 1, im
          sfalb(i) = max(0.01, 0.5 * (sfcalb(i,2) + sfcalb(i,4)))
        enddo

        if (me == 0 .and. myrank ==0) then
           print *,'###################################################' 
           print *,'### call setalb ok !!'
           print *,'###################################################' 
           print *,'### slmsk(ipt)=',slmsk(ipt)
           print *,'### fice(ipt)=',fice(ipt)
           print *,'### tisfc(ipt)=',tisfc(ipt)
           print *,'###################################################' 
           print *,'### snowd(ipt)=',snowd(ipt)
           print *,'### sncovr(ipt)=',sncovr(ipt)
           print *,'### snoalb(ipt)=',snoalb(ipt)
           print *,'###################################################' 
           print *,'### zorl(ipt)=',zorl(ipt)
           print *,'### hprim(ipt)=',hprim(ipt)
           print *,'###################################################' 
           print *,'### alvsf(ipt)=',alvsf(ipt)
           print *,'### alnsf(ipt)=',alnsf(ipt)
           print *,'### alvwf(ipt)=',alvwf(ipt)
           print *,'### alnwf(ipt)=',alnwf(ipt)
           print *,'### facsf(ipt)=',facsf(ipt)
           print *,'### facwf(ipt)=',facwf(ipt)
           print *,'###################################################' 
           print *,'### coszen(ipt)=',coszen(ipt)
           print *,'### tsfg(ipt)=',tsfg(ipt)
           print *,'### tsfa(ipt)=',tsfa(ipt)
           print *,'###################################################' 
           print *,'### sfcalb(ipt,1)-near ir direct beam albedo'
           print *,'### sfcalb(ipt,2)-near ir diffused beam albedo'
           print *,'### sfcalb(ipt,3)-uv+vis direct beam albedo'
           print *,'### sfcalb(ipt,4)-uv+vis diffused beam albedo'
           print *,'###################################################' 
           print *,'### sfcalb(ipt,1)=',sfcalb(ipt,1)
           print *,'### sfcalb(ipt,2)=',sfcalb(ipt,2)
           print *,'### sfcalb(ipt,3)=',sfcalb(ipt,3)
           print *,'### sfcalb(ipt,4)=',sfcalb(ipt,4)
           print *,'###################################################' 
           print *,'# sfalb(i)-average diffused beam albedo for sw&lw'
           print *,'# sfalb(i)=max(0.01,0.5*(sfcalb(i,2)+sfcalb(i,4)))'
           print *,'###################################################' 
           print *,'### sfalb(ipt)=',sfalb(ipt)
           print *,'###################################################' 
        endif

      if (nday > 0) then
          if (me == 0 .and. myrank ==0) print *,' #### call swrad ####'
          if (me == 0 .and. myrank ==0) print *,' #### nday=',nday
 
!         if ( present(htrswb) .and. present(htrsw0) ) then
!         if ( present(htrsw0) .and. present(fswprf) ) then

            call swrad                                                  &
!  ---  inputs:
     &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                      &
     &       clouds,icsdsw,faersw,sfcalb,                               &
     &       coszen,solcon, nday,idxday,                                &
     &       im, lmk, lmp, lprnt, myrank,                               &
!  ---  outputs:
     &       htswc,topfsw,sfcfsw,                                       &
!! ---  optional:
     &       hsw0=htsw0,hswb=htswb,                                     &
     &       flxprf=fswprf,fdncmp=scmpsw                                &
     &     )

            do k = 1, lm
              k1 = k + kd
              do j = 1, nbdsw
                do i = 1, im
                  htrswb(i,k,j) = htswb(i,k1,j)
                enddo
              enddo
            enddo

!         else if ( present(htrswb) .and. .not. present(htrsw0) ) then

!           call swrad                                                  &
!  ---  inputs:
!    &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                      &
!    &       clouds,icsdsw,faersw,sfcalb,                               &
!    &       coszen,solcon, nday,idxday,                                &
!    &       im, lmk, lmp, lprnt,                                       &
!  ---  outputs:
!    &       htswc,topfsw,sfcfsw                                        &
!! ---  optional:
!    &,      hsw0=htsw0,flxprf=fswprf                                   &
!    &,      hswb=htswb,fdncmp=scmpsw                                   &
!    &     )

!           do k = 1, lm
!             k1 = k + kd
!             do j = 1, nbdsw
!               do i = 1, im
!                 htrswb(i,k,j) = htswb(i,k1,j)
!               enddo
!             enddo
!           enddo

!         else if ( present(htrsw0) .and. .not. present(htrswb) ) then

!           call swrad                                                  &
!  ---  inputs:
!    &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                      &
!    &       clouds,icsdsw,faersw,sfcalb,                               &
!    &       coszen,solcon, nday,idxday,                                &
!    &       im, lmk, lmp, lprnt,                                       &
!  ---  outputs:
!    &       htswc,topfsw,sfcfsw                                        &
!! ---  optional:
!!   &,      hsw0=htsw0,flxprf=fswprf                                   &
!    &,      hsw0=htsw0,fdncmp=scmpsw                                   &
!    &     )

!         else

!           call swrad                                                  &
!  ---  inputs:
!    &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                      &
!    &       clouds,icsdsw,faersw,sfcalb,                               &
!    &       coszen,solcon, nday,idxday,                                &
!    &       im, lmk, lmp, lprnt,                                       &
!  ---  outputs:
!    &       htswc,topfsw,sfcfsw                                        &
!! ---  optional:
!!   &,      hsw0=htsw0,flxprf=fswprf,hswb=htswb                        &
!!   &,      fdncmp=scmpsw                                              &
!    &     )

!         endif

          do k = 1, lm
            k1 = k + kd
            do i = 1, im
              htrsw(i,k) = htswc(i,k1)
            enddo
          enddo
!         if (present(htrsw0)) then
             do k = 1, lm
               k1 = k + kd
               do i = 1, im
                 htrsw0(i,k) = htsw0(i,k1)
               enddo
             enddo
!         endif
!         if (present(fswprf)) then
             do k = 1, lm+1
                k1 = k + kd
               do i = 1, im
                 fusl(i,k)  = fswprf(i,k1)%upfxc
                 fdsl(i,k)  = fswprf(i,k1)%dnfxc
                 fuslr(i,k) = fswprf(i,k1)%upfx0
                 fdslr(i,k) = fswprf(i,k1)%dnfx0
               enddo
             enddo
!         endif


        else                   ! if_nday_block

          do k = 1, lm
            do i = 1, im
              htrsw(i,k) = 0.0
            enddo
          enddo

          sfcfsw = sfcfsw_type( 0.0, 0.0, 0.0, 0.0 )
          topfsw = topfsw_type( 0.0, 0.0, 0.0 )
          scmpsw = cmpfsw_type( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 )

!! ---  optional:
          fswprf= profsw_type( 0.0, 0.0, 0.0, 0.0 )

!         if ( present(htrswb) ) then
            do j = 1, nbdsw
              do k = 1, lm
                do i = 1, im
                  htrswb(i,k,j) = 0.0
                enddo
              enddo
            enddo
!         endif

!         if ( present(htrsw0) ) then
              do k = 1, lm
                do i = 1, im
                  htrsw0(i,k) = 0.0
                enddo
              enddo
!         endif

        endif                  ! end_if_nday

      endif                                ! end_if_lsswr

       if (me == 0 .and. myrank ==0) then
           print *,'###################################################' 
           print *,'### call swrad ok!!'
           print *,'###################################################' 
           print *,'### im=',im,' lmk=',lmk,' lmp=',lmp
           print *,'### nday=',nday
           print *,'### solcon=',solcon
           print *,'###################################################' 
           print *,'### icsdsw(ipt)=',icsdsw(ipt)
           print *,'###################################################' 
           print *,'### plvl(ipt,lm)=',plvl(ipt,lm)
           print *,'### plyr(ipt,lm)=',plyr(ipt,lm)
           print *,'###################################################' 
           print *,'### tlyr(ipt,lm)=',tlyr(ipt,lm)
           print *,'### tlvl(ipt,lm)=',tlvl(ipt,lm)
           print *,'###################################################' 
           print *,'### qlyr(ipt,lm)=',qlyr(ipt,lm)
           print *,'### olyr(ipt,lm)=',olyr(ipt,lm)
           print *,'###################################################' 
           print *,'### gasvmr(ipt,lm,1)_co2 =',gasvmr(ipt,lm,1)
           print *,'### gasvmr(ipt,lm,2)_n2o =',gasvmr(ipt,lm,2)
           print *,'### gasvmr(ipt,lm,3)_ch4 =',gasvmr(ipt,lm,3)
           print *,'### gasvmr(ipt,lm,4)_o2  =',gasvmr(ipt,lm,4)
           print *,'### gasvmr(ipt,lm,5)_co  =',gasvmr(ipt,lm,5)
           print *,'### gasvmr(ipt,lm,6)_cf11=',gasvmr(ipt,lm,6)
           print *,'### gasvmr(ipt,lm,7)_cf12=',gasvmr(ipt,lm,7)
           print *,'### gasvmr(ipt,lm,8)_cf22=',gasvmr(ipt,lm,8)
           print *,'### gasvmr(ipt,lm,9)_ccl4=',gasvmr(ipt,lm,9)
           print *,'###################################################' 
           print *,'### clouds(ipt,lm,1)-total cloud fraction=',clouds(ipt,lm,1)
           print *,'### clouds(ipt,lm,2)-liq water path=',clouds(ipt,lm,2)
           print *,'### clouds(ipt,lm,3)-liq effective radius=',clouds(ipt,lm,3)
           print *,'### clouds(ipt,lm,4)-ice water path=',clouds(ipt,lm,4)
           print *,'### clouds(ipt,lm,5)-ice effective radius=',clouds(ipt,lm,5)
           print *,'### clouds(ipt,lm,6)-rain water path=',clouds(ipt,lm,6)
           print *,'### clouds(ipt,lm,7)-rain effective radius=',clouds(ipt,lm,7)
           print *,'### clouds(ipt,lm,8)-snow water path=',clouds(ipt,lm,8)
           print *,'### clouds(ipt,lm,9)-snow effective radius=',clouds(ipt,lm,9)
           print *,'###################################################' 
           print *,'### faersw(ipt,lm,1,1)=sw#1-opd',faersw(ipt,lm,1,1)
           print *,'### faersw(ipt,lm,1,2)=sw#1-ssa',faersw(ipt,lm,1,2)
           print *,'### faersw(ipt,lm,1,3)=sw#1-asy',faersw(ipt,lm,1,3)
           print *,'###################################################' 
           print *,'### sfcalb(ipt,1)=',sfcalb(ipt,1)
           print *,'### sfcalb(ipt,2)=',sfcalb(ipt,2)
           print *,'### sfcalb(ipt,3)=',sfcalb(ipt,3)
           print *,'### sfcalb(ipt,4)=',sfcalb(ipt,4)
           print *,'###################################################' 
           print *,'### htswc(ipt,lm) =',htswc(ipt,lm)
           print *,'### htrsw(ipt,lm) =',htrsw(ipt,lm)
           print *,'### htrsw0(ipt,lm)=',htrsw0(ipt,lm)
           print *,'###################################################' 
           print *,'### sfcfsw(ipt)%upfxc=',sfcfsw(ipt)%upfxc
           print *,'### sfcfsw(ipt)%dnfxc=',sfcfsw(ipt)%dnfxc
           print *,'### sfcfsw(ipt)%upfx0=',sfcfsw(ipt)%upfx0
           print *,'### sfcfsw(ipt)%dnfx0=',sfcfsw(ipt)%dnfx0
           print *,'###################################################' 
           print *,'### topfsw(ipt)%upfxc=',topfsw(ipt)%upfxc
           print *,'### topfsw(ipt)%dnfxc=',topfsw(ipt)%dnfxc
           print *,'### topfsw(ipt)%upfx0=',topfsw(ipt)%upfx0
           print *,'###################################################' 
           print *,'### scmpsw(ipt)%uvbfc=',scmpsw(ipt)%uvbfc
           print *,'### scmpsw(ipt)%uvbf0=',scmpsw(ipt)%uvbf0
           print *,'### scmpsw(ipt)%nirbm=',scmpsw(ipt)%nirbm
           print *,'### scmpsw(ipt)%nirdf=',scmpsw(ipt)%nirdf
           print *,'### scmpsw(ipt)%visbm=',scmpsw(ipt)%visbm
           print *,'### scmpsw(ipt)%visdf=',scmpsw(ipt)%visdf
           print *,'###################################################' 
           print *,'### fswprf(ipt,1)%upfxc=',fswprf(ipt,1)%upfxc
           print *,'### fswprf(ipt,1)%dnfxc=',fswprf(ipt,1)%dnfxc
           print *,'### fswprf(ipt,1)%upfx0=',fswprf(ipt,1)%upfx0
           print *,'### fswprf(ipt,1)%dnfx0=',fswprf(ipt,1)%dnfx0
           print *,'###################################################' 
           print *,'### fswprf(ipt,61)%upfxc=',fswprf(ipt,61)%upfxc
           print *,'### fswprf(ipt,61)%dnfxc=',fswprf(ipt,61)%dnfxc
           print *,'### fswprf(ipt,61)%upfx0=',fswprf(ipt,61)%upfx0
           print *,'### fswprf(ipt,61)%dnfx0=',fswprf(ipt,61)%dnfx0
           print *,'###################################################' 
           print *,'### fusl(ipt,1)=',fusl(ipt,1)
           print *,'### fdsl(ipt,1)=',fdsl(ipt,1)
           print *,'### fuslr(ipt,1)=',fuslr(ipt,1)
           print *,'### fdslr(ipt,1)=',fdslr(ipt,1)
           print *,'###################################################' 
           print *,'### fusl(ipt,61)=',fusl(ipt,61)
           print *,'### fdsl(ipt,61)=',fdsl(ipt,61)
           print *,'### fuslr(ipt,61)=',fuslr(ipt,61)
           print *,'### fdslr(ipt,61)=',fdslr(ipt,61)
           print *,'###################################################' 
      endif

!----------------------------------------------------------------------
      if (lslwr) then

!  ---  setup surface emissivity for lw radiation

        call setemis                                                    &
!  ---  inputs:
     &     ( xlon,xlat,slmsk,snowd,sncovr,zorl,tsfg,tsfa,hprim,         &
     &       im,                                                        &
!  ---  outputs:
     &       sfcemis                                                    &
     &     )
!
       if (me == 0 .and. myrank ==0) then
           print *,'###################################################' 
           print *,'### call setemis ok!!'
           print *,'###################################################' 
           print *,'### im=',im,' ipt=',ipt
           print *,'###################################################' 
           print *,'### xlon(ipt)=',xlon(ipt)
           print *,'### xlat(ipt)=',xlat(ipt)
           print *,'###################################################' 
           print *,'### slmsk(ipt)=',slmsk(ipt)
           print *,'### snowd(ipt)=',snowd(ipt)
           print *,'### sncovr(ipt)=',sncovr(ipt)
           print *,'### zorl(ipt)=',zorl(ipt)
           print *,'###################################################' 
           print *,'### tsfg(ipt)=',tsfg(ipt)
           print *,'### tsfa(ipt)=',tsfa(ipt)
           print *,'### hprim(ipt)=',hprim(ipt)
           print *,'###################################################' 
           print *,'### sfcemis(ipt)=',sfcemis(ipt)
           print *,'###################################################' 
      endif


!       if ( present(htrlw0) .and. present(flwprf) ) then
        if (me == 0 .and. myrank ==0) print *,'#### call lwrad ####'
          call lwrad                                                    &
!  ---  inputs:
     &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                      &
     &       clouds,icsdlw,faerlw,sfcemis,tsfg,                         &
     &       im, lmk, lmp, lprnt,myrank,                                &
!  ---  outputs:
     &       htlwc,topflw,sfcflw                                        &
!! ---  optional:
     &,      hlw0=htlw0,hlwb=htlwb,flxprf=flwprf                        &
     &     )

          do k = 1, lm
            k1 = k + kd
            do j = 1, nbdlw
              do i = 1, im
                htrlwb(i,k,j) = htlwb(i,k1,j)
              enddo
            enddo
          enddo

!       else if ( present(htrlwb) .and. .not. present(htrlw0) ) then

!         call lwrad                                                    &
!  ---  inputs:
!    &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                      &
!    &       clouds,icsdlw,faerlw,sfcemis,tsfg,                         &
!    &       im, lmk, lmp, lprnt,                                       &
!  ---  outputs:
!    &       htlwc,topflw,sfcflw                                        &
!! ---  optional:
!!   &,      hlw0=htlw0,flxprf=flwprf                                   &
!    &,      hlwb=htlwb                                                 &
!    &     )

!         do k = 1, lm
!           k1 = k + kd

!           do j = 1, nbdlw
!             do i = 1, im
!               htrlwb(i,k,j) = htlwb(i,k1,j)
!             enddo
!           enddo
!         enddo

!       else if ( present(htrlw0) .and. .not. present(htrlwb) ) then

!         !print *,'call lwrad saving clear sky component'
!         call lwrad                                                    &
!  ---  inputs:
!    &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                      &
!    &       clouds,icsdlw,faerlw,sfcemis,tsfg,                         &
!    &       im, lmk, lmp, lprnt,                                       &
!  ---  outputs:
!    &       htlwc,topflw,sfcflw                                        &
!! ---  optional:
!!   &,      hlw0=htlw0,flxprf=flwprf                                   &
!    &,      hlw0=htlw0                                                 &
!    &     )

!       else

!         call lwrad                                                    &
!  ---  inputs:
!    &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                      &
!    &       clouds,icsdlw,faerlw,sfcemis,tsfg,                         &
!    &       im, lmk, lmp, lprnt,                                       &
!  ---  outputs:
!    &       htlwc,topflw,sfcflw                                        &
!! ---  optional:
!!   &,      hlw0=htlw0,flxprf=flwprf,hlwb=htlwb                        &
!    &     )

!       endif

        do i = 1, im
          semis (i) = sfcemis(i)
!  ---  save surface air temp for diurnal adjustment at model t-steps
          tsflw (i) = tsfa(i)
        enddo

        do k = 1, lm
          k1 = k + kd
          do i = 1, im
            htrlw(i,k) = htlwc(i,k1)
          enddo
        enddo

!       if (present(htrlw0)) then
           do k = 1, lm
             k1 = k + kd
             do i = 1, im
               htrlw0(i,k) = htlw0(i,k1)
             enddo
           enddo
!       endif
!         if (present(flwprf)) then
             do k = 1, lm+1
                k1 = k + kd
               do i = 1, im
                 fuir(i,k)  = flwprf(i,k1)%upfxc
                 fdir(i,k)  = flwprf(i,k1)%dnfxc
                 fuirr(i,k) = flwprf(i,k1)%upfx0
                 fdirr(i,k) = flwprf(i,k1)%dnfx0
               enddo
             enddo
!          endif

      endif                                ! end_if_lslwr

!----------------------------------------------------------------------
       if (me == 0 .and. myrank ==0) then
           print *,'###################################################' 
           print *,'### call lwrad ok!!'
           print *,'###################################################' 
           print *,'### im=',im,' lmk=',lmk,' lmp=',lmp,' kd=',kd
           print *,'###################################################' 
           print *,'### icsdlw(ipt)=',icsdlw(ipt)
           print *,'###################################################' 
           print *,'### plvl(ipt,lm)=',plvl(ipt,lm)
           print *,'### plyr(ipt,lm)=',plyr(ipt,lm)
           print *,'###################################################' 
           print *,'### tlyr(ipt,lm)=',tlyr(ipt,lm)
           print *,'### tlvl(ipt,lm)=',tlvl(ipt,lm)
           print *,'###################################################' 
           print *,'### qlyr(ipt,lm)=',qlyr(ipt,lm)
           print *,'### olyr(ipt,lm)=',olyr(ipt,lm)
           print *,'###################################################' 
           print *,'### gasvmr(ipt,lm,1)_co2 =',gasvmr(ipt,lm,1)
           print *,'### gasvmr(ipt,lm,2)_n2o =',gasvmr(ipt,lm,2)
           print *,'### gasvmr(ipt,lm,3)_ch4 =',gasvmr(ipt,lm,3)
           print *,'### gasvmr(ipt,lm,4)_o2  =',gasvmr(ipt,lm,4)
           print *,'### gasvmr(ipt,lm,5)_co  =',gasvmr(ipt,lm,5)
           print *,'### gasvmr(ipt,lm,6)_cf11=',gasvmr(ipt,lm,6)
           print *,'### gasvmr(ipt,lm,7)_cf12=',gasvmr(ipt,lm,7)
           print *,'### gasvmr(ipt,lm,8)_cf22=',gasvmr(ipt,lm,8)
           print *,'### gasvmr(ipt,lm,9)_ccl4=',gasvmr(ipt,lm,9)
           print *,'###################################################' 
           print *,'### clouds(ipt,lm,1)-total cloud fraction=',clouds(ipt,lm,1)
           print *,'### clouds(ipt,lm,2)-liq water path=',clouds(ipt,lm,2)
           print *,'### clouds(ipt,lm,3)-liq effective radius=',clouds(ipt,lm,3)
           print *,'### clouds(ipt,lm,4)-ice water path=',clouds(ipt,lm,4)
           print *,'### clouds(ipt,lm,5)-ice effective radius=',clouds(ipt,lm,5)
           print *,'### clouds(ipt,lm,6)-rain water path=',clouds(ipt,lm,6)
           print *,'### clouds(ipt,lm,7)-rain effective radius=',clouds(ipt,lm,7)
           print *,'### clouds(ipt,lm,8)-snow water path=',clouds(ipt,lm,8)
           print *,'### clouds(ipt,lm,9)-snow effective radius=',clouds(ipt,lm,9)
           print *,'###################################################' 
           print *,'### faerlw(ipt,lm,1,1)=lw#1-opd',faerlw(ipt,lm,1,1)
           print *,'### faerlw(ipt,lm,1,2)=lw#1-ssa',faerlw(ipt,lm,1,2)
           print *,'### faerlw(ipt,lm,1,3)=lw#1-asy',faerlw(ipt,lm,1,3)
           print *,'###################################################' 
           print *,'### sfcemis(ipt)=',sfcemis(ipt)
           print *,'###################################################' 
           print *,'### tsfg(ipt) =',tsfg(ipt)
           print *,'### tsfa(ipt) =',tsfa(ipt)
           print *,'### tsflw(ipt)=',tsflw(ipt)
           print *,'###################################################' 
           print *,'### htswc(ipt,lm) =',htswc(ipt,lm)
           print *,'### htsw0(ipt,lm) =',htsw0(ipt,lm)
           print *,'### htrsw(ipt,lm) =',htrsw(ipt,lm)
           print *,'### htrsw0(ipt,lm)=',htrsw0(ipt,lm)
           print *,'###################################################' 
           print *,'### htlwc(ipt,lm) =',htlwc(ipt,lm)
           print *,'### htlw0(ipt,lm) =',htlw0(ipt,lm)
           print *,'### htrlw(ipt,lm) =',htrlw(ipt,lm)
           print *,'### htrlw0(ipt,lm)=',htrlw0(ipt,lm)
           print *,'###################################################' 
           print *,'### topflw(ipt)%upfxc=',topflw(ipt)%upfxc
           print *,'### topflw(ipt)%upfx0=',topflw(ipt)%upfx0
           print *,'###################################################' 
           print *,'### topfsw(ipt)%upfxc=',topfsw(ipt)%upfxc
           print *,'### topfsw(ipt)%upfx0=',topfsw(ipt)%upfx0
           print *,'### topfsw(ipt)%dnfxc=',topfsw(ipt)%dnfxc
           print *,'###################################################' 
           print *,'### sfcflw(ipt)%upfxc=',sfcflw(ipt)%upfxc
           print *,'### sfcflw(ipt)%upfx0=',sfcflw(ipt)%upfx0
           print *,'### sfcflw(ipt)%dnfxc=',sfcflw(ipt)%dnfxc
           print *,'### sfcflw(ipt)%dnfx0=',sfcflw(ipt)%dnfx0
           print *,'###################################################' 
           print *,'### sfcfsw(ipt)%upfxc=',sfcfsw(ipt)%upfxc
           print *,'### sfcfsw(ipt)%upfx0=',sfcfsw(ipt)%upfx0
           print *,'### sfcfsw(ipt)%dnfxc=',sfcfsw(ipt)%dnfxc
           print *,'### sfcfsw(ipt)%dnfx0=',sfcfsw(ipt)%dnfx0
           print *,'###################################################' 
           print *,'### flwprf(ipt,1)%upfxc=',flwprf(ipt,1)%upfxc
           print *,'### flwprf(ipt,1)%dnfxc=',flwprf(ipt,1)%dnfxc
           print *,'### flwprf(ipt,1)%upfx0=',flwprf(ipt,1)%upfx0
           print *,'### flwprf(ipt,1)%dnfx0=',flwprf(ipt,1)%dnfx0
           print *,'###################################################' 
           print *,'### flwprf(ipt,61)%upfxc=',flwprf(ipt,61)%upfxc
           print *,'### flwprf(ipt,61)%dnfxc=',flwprf(ipt,61)%dnfxc
           print *,'### flwprf(ipt,61)%upfx0=',flwprf(ipt,61)%upfx0
           print *,'### flwprf(ipt,61)%dnfx0=',flwprf(ipt,61)%dnfx0
           print *,'###################################################' 
           print *,'### fuir(ipt,1)=',fuir(ipt,1)
           print *,'### fdir(ipt,1)=',fdir(ipt,1)
           print *,'### fuirr(ipt,1)=',fuirr(ipt,1)
           print *,'### fdirr(ipt,1)=',fdirr(ipt,1)
           print *,'###################################################' 
           print *,'### fuir(ipt,61)=',fuir(ipt,61)
           print *,'### fdir(ipt,61)=',fdir(ipt,61)
           print *,'### fuirr(ipt,61)=',fuirr(ipt,61)
           print *,'### fdirr(ipt,61)=',fdirr(ipt,61)
           print *,'###################################################' 
       endif
!
!
!  ---  save total-sky TOA and SFC fluxes

          do i = 1, im

!  ---  TOA total-sky SW DOWN fluxes
 
            fluxr(i,1 ) = topfsw(i)%dnfxc   ! total sky top sw dn

!  ---  TOA total-sky SW UP fluxes

            fluxr(i,2 ) = topfsw(i)%upfxc   ! total sky top sw up

!  ---  TOA total-sky LW(OLR) fluxes

            fluxr(i,3 ) = topflw(i)%upfxc   ! total sky top lw up

!  ---  SFC total-sky SW DN/UP fluxes

            fluxr(i,4 ) = sfcfsw(i)%dnfxc   ! total sky sfc sw dn
            fluxr(i,5 ) = sfcfsw(i)%upfxc   ! total sky sfc sw up
   
!  ---  SFC total-sky LW DN/UP fluxes

            fluxr(i,6 ) = sfcflw(i)%dnfxc  ! total sky sfc lw dn
            fluxr(i,7 ) = sfcflw(i)%upfxc  ! total sky sfc lw up

          enddo
          

!  ---  save cld frac,toplyr,botlyr and top temp, note that the order
!       of h,m,l cloud is reversed for the fluxr output.
!  ---  save interface pressure (cb) of top/bot

          do j = 1, 3
            do i = 1, im
              tem0d = cldsa(i,j)
              itop  = mtopa(i,j) - kd
              ibtc  = mbota(i,j) - kd
              fluxr(i,11-j) = tem0d  ! cloud fraction(h,m,l,j=3,2,1)
              fluxr(i,14-j) = prsi(i,itop+kt) ! cloud top pressure
              fluxr(i,17-j) = prsi(i,ibtc+kb) ! cloud bot pressure
              fluxr(i,20-j) = tgrs(i,itop)    ! cloud top temp
            enddo
          enddo

!  ---  save total cloud and bl cloud fraction

          do i = 1, im

            fluxr(i,20) = cldsa(i,4) ! total cloud fraction
            fluxr(i,21) = cldsa(i,5) ! BL domain cloud fraction

          enddo

!  ---  save clear-sky TOA and SFC fluxes

          do i = 1, im

!  ---  TOA clear-sky SW/LW UP fluxes

            fluxr(i,22) = topfsw(i)%upfx0   ! clear sky top sw up
            fluxr(i,23) = topflw(i)%upfx0   ! clear sky top lw up
           
!  ---  SFC clear-sky SW fluxes

            fluxr(i,24) = sfcfsw(i)%dnfx0   ! clear sky sfc sw dn
            fluxr(i,25) = sfcfsw(i)%upfx0   ! clear sky sfc sw up
!
!  ---  SFC clear-sky LW fluxes

            fluxr(i,26) = sfcflw(i)%dnfx0   ! clear sky sfc lw dn
            fluxr(i,27) = sfcflw(i)%upfx0   ! clear sky sfc lw up
!
          enddo
!
!  ---  sw uv-b fluxes

          do i = 1, im

            fluxr(i,28) = scmpsw(i)%uvbfc  ! total sky uv-b sw dn
            fluxr(i,29) = scmpsw(i)%uvbf0  ! clear sky uv-b sw dn

!  ---  sw sfc flux components
!
            fluxr(i,30) = scmpsw(i)%visbm  ! uv/vis beam sw dn
            fluxr(i,31) = scmpsw(i)%visdf  ! uv/vis diff sw dn
            fluxr(i,32) = scmpsw(i)%nirbm  ! nir beam sw dn
            fluxr(i,33) = scmpsw(i)%nirdf  ! nir diff sw dn

          enddo


        do k = 1, lm
          k1 = k + kd

          do i = 1, im
            cldcov(i,k) = clouds(i,k1,1)
          enddo
        enddo
 

!  ---  save optional vertically integrated aerosol optical depth at
!       wavelenth of 550nm aerodp(:,1), and other optional aod for
!       individual species aerodp(:,2:nspc1)

!       if ( laswflg ) then
!         if ( nfxr > 33 ) then
!           do i = 1, im
!             fluxr(i,34) = fluxr(i,34) + dtsw*aerodp(i,1)  ! total aod at 550nm (all species)
!           enddo

!           if ( lspcodp ) then
!             do j = 2, nspc1
!               k = 33 + j

!               do i = 1, im
!                 fluxr(i,k) = fluxr(i,k) + dtsw*aerodp(i,j) ! aod at 550nm for indiv species
!               enddo
!             enddo
!           endif     ! end_if_lspcodp
!         else
!           print *,'  !error! need to increase array fluxr size nfxr ',&
!    &              ' to be able to output aerosol optical depth'
!           stop
!         endif     ! end_if_nfxr
!       endif       ! end_if_laswflg

!
      return
!...................................
      end subroutine grrad
!-----------------------------------


!
!........................................!
      end module module_radiation_driver !
!========================================!
