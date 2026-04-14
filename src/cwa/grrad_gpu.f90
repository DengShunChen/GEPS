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
module module_radiation_driver_gpu     !
!........................................!
!
   use physpara
   use physcons,                 only : con_eps, con_epsm1, con_fvirt&
   &,                                   rocp => con_rocp
!  use funcphys,                 only : fpvs
   use radsw_copyin_gpu, only : copyin_radsw_datatb_gpu
   use radlw_copyin_gpu, only : copyin_radlw_datatb_gpu
   use module_radiation_astronomy_gpu,only: sol_init_gpu, sol_update_gpu, coszmn_gpu, &
                                         copyin_radiation_astronomy_gpu
   use module_radiation_gases_gpu,   only : nf_vgas, getgases_gpu, getozn_gpu,   &
   &                                    gas_init_gpu, gas_update_gpu, &
                                        copyin_radiation_gases_gpu
   use module_radiation_aerosols_gpu,only : nf_aesw, nf_aelw, setaer_sw_gpu, setaer_lw_gpu,    &
   &                                    aer_init_gpu, aer_update_gpu, &
                                        copyin_radiation_aerosols_gpu
!    &,                                    nspc1                        ! optn for aod output
   use module_radiation_surface_gpu, only : nf_albd, sfc_init_gpu, setalb_gpu,   &
   &                                    setemis_gpu, copyin_radiation_surface_gpu
   use module_radiation_clouds_gpu,  only : nf_clds, cld_init_gpu,           &
   &                                    progcld1_gpu, progcld2_gpu, progcld3_gpu,&
   &					                      progcld4_gpu, diagcld1_gpu,          &
                                        progcld5_gpu, progcld5o_gpu,         &
                                        progclduni_gpu, progcld6_gpu,        &
                                        progcld_thompson_gpu, progcld_gce_gpu, &
                                        copyin_radiation_clouds_gpu

   use module_radsw_parameters,  only : topfsw_type, sfcfsw_type,    &
   &                                    profsw_type,cmpfsw_type,nbdsw

   use module_radlw_parameters,  only : topflw_type, sfcflw_type,    &
   &                                    proflw_type, nbdlw
   use module_radlw_main_gpu,        only : rlwinit_gpu,  lwrad_gpu, copyin_radlw_main_gpu
   use module_radsw_main_gpu,        only : rswinit_gpu,  swrad_gpu, copyin_radsw_main_gpu
   use param,                    only : my, my_max
   use index,                    only : jlistnum, nxjp, jlist1, nxptot, nxjp_acc
   !use nvtx
!
!    implicit   none
!
!
   private

!  ---  version tag and last revision date
   character(40), parameter ::                                       &
&     vtagrad='ncep-radiation_driver    v5.2  jan 2013 '
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

   public radinit_gpu, radupdate_gpu, grrad_gpu


! =================
contains
! =================


!-----------------------------------
   subroutine radinit_gpu                                                &
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
      call sol_init_gpu ( me , myrank )

!  --- ...  call aerosols initialization routine

      if (me == 0 .and. myrank == 0)print *,'call aer_init' 
      call aer_init_gpu ( nlay, me , myrank )

!  --- ...  call co2 and other gases initialization routine

      if (me ==0 .and. myrank == 0)print *,'call gas_init' 
      call gas_init_gpu ( me , myrank )

!  --- ...  call surface initialization routine

      if (me == 0 .and. myrank == 0)print *,'call sfc_init' 
      call sfc_init_gpu ( me , myrank )

!  --- ...  call cloud initialization routine

      if (me == 0 .and. myrank == 0)print *,'call cld_init' 
      call cld_init_gpu ( si, nlay, me, myrank )

!  --- ...  call lw radiation initialization routine

      if (me == 0 .and. myrank == 0)print *,'call rlw_init' 
      call rlwinit_gpu ( me , myrank )

!  --- ...  call sw radiation initialization routine

      if (me == 0 .and. myrank == 0)print *,'call rsw_init' 
      call rswinit_gpu ( me ,myrank)
      
      call copyin_radiation_clouds_gpu(1)
      call copyin_radiation_surface_gpu(1)
      call copyin_radiation_gases_gpu(1)
      call copyin_radiation_aerosols_gpu(1)
      call copyin_radiation_astronomy_gpu(1)
      call copyin_radsw_datatb_gpu(1)
      call copyin_radlw_datatb_gpu(1)
      call copyin_radlw_main_gpu(1)
      call copyin_radsw_main_gpu(1)
!
      return
!...................................
   end subroutine radinit_gpu
!-----------------------------------


!-----------------------------------
   subroutine radupdate_gpu                                              &
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
         call sol_update_gpu                                                 &
         !  ---  inputs:
         &     ( myrank, jdate,kyear,deltsw,deltim,lsol_chg, me,            &
         !  ---  outputs:
         &       slag,sdec,cdec,solcon                                      &
         &     )

      endif  ! end_if_lsswr_block

!  --- ...  call aerosols update routine, monthly update, no time interpolation

      if ( lmon_chg ) then
         if (me == 0 .and. myrank ==0) print *,'call aer_update'
         call aer_update_gpu ( myrank, iyear, imon, me )
      endif

!  --- ...  call co2 and other gases update routine

      if ( monthd /= kmon ) then
         monthd = kmon
         lco2_chg = .true.
      else
         lco2_chg = .false.
      endif

      if (me == 0 .and. myrank ==0) print *,'call gas_update'
         call gas_update_gpu ( myrank, kyear,kmon,kday,khour,loz1st,lco2_chg, me )

      if ( loz1st ) loz1st = .false.

!  --- ...  call surface update routine (currently not needed)
!     call sfc_update ( iyear, imon, me )

!  --- ...  call clouds update routine (currently not needed)
!     call cld_update ( iyear, imon, me )
!
      return
!...................................
   end subroutine radupdate_gpu
!-----------------------------------


!-----------------------------------
   subroutine grrad_gpu                                                  &
!...................................
!  ---  inputs:
      ( prsi,prsl,prslk,tgrs,qgrs,tracer,vvl,slmsk,                &
      xlon,xlat,tsfc,snowd,sncovr,snoalb,zorl,hprim,             &
      alvsf,alnsf,alvwf,alnwf,facsf,facwf,fice,tisfc,            &
      sinlat,coslat,solhr,jdate,solcon,                          &
      cv,cvt,cvb,                                                &
      icsdsw,icsdlw,ntcw,ncld,ntoz,ntrac,nfxr,                   &
      dtlw,dtsw,lsswr,lslwr,lssav,                               &
      ix,myim,lm,me,lprnt,ipt,kdt,myrank,                          &
      ntiw,ntrw,ntsw,ntgl,uni_cloud,lmfshal,lmfdeep2,            &
      deltaq,sup,cnvw,cnvc,phy_f3d, async_id,                              &
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
      !$acc routine(fpvs_gpu) seq

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
      integer,  intent(in) :: ix,myim(my_max), lm, ntrac, nfxr, me,          &
         ntoz, ntcw, ncld, ipt, kdt, myrank,  &
         ntiw, ntrw, ntsw, ntgl
      integer,  intent(in) :: icsdsw(ix, my_max), icsdlw(ix, my_max), jdate(8)

      logical,  intent(in) :: lsswr, lslwr, lssav, lprnt

      real (kind=kind_phys), dimension(ix,lm+1, my_max), intent(in) ::  prsi

      real (kind=kind_phys), dimension(ix,lm, my_max),   intent(in) ::  prsl,   &
!            prslk, tgrs, qgrs, vvl, fcice, frain, rrime, deltaq, cnvw, & 
!            cnvc
         prslk, tgrs, qgrs, vvl, deltaq, cnvw, cnvc
!     real (kind=kind_phys), dimension(im), intent(in) :: flgmin
      real(kind=kind_phys), intent(in) ::sup

      real (kind=kind_phys), dimension(ix, my_max),      intent(in) ::  slmsk,  &
         xlon, xlat, tsfc, snowd, zorl, hprim, alvsf, alnsf, alvwf, &
         alnwf, facsf, facwf, cv, cvt, cvb, fice, tisfc,            &
         sncovr, snoalb, sinlat, coslat

      real (kind=kind_phys), intent(in) :: solcon, dtlw, dtsw, solhr,   &
         tracer(ix,lm,ntrac, my_max)

!  ---  outputs: (horizontal dimensioned by ix)
      real (kind=kind_phys), dimension(ix,lm, my_max),intent(out):: htrsw,htrlw

      real (kind=kind_phys), dimension(ix, my_max),   intent(out):: tsflw,      &
         sfalb, semis, coszen, coszdg

! --- cmy
      real (kind=kind_phys), dimension(ix,lm+1+ltp, my_max)::                   &
         fusl,fdsl,fuir,fdir,fuslr,fdslr,fuirr,fdirr
! --- cmy

!     type (topfsw_type), dimension(im), intent(out) :: topfsw
!     type (sfcfsw_type), dimension(im), intent(out) :: sfcfsw
! --- cmy
      !real (kind=kind_phys), dimension(ix, my_max, 3) :: topfsw
      !real (kind=kind_phys), dimension(ix, my_max, 4) :: sfcfsw
! --- cmy

!     type (topflw_type), dimension(im), intent(out) :: topflw
!     type (sfcflw_type), dimension(im), intent(out) :: sfcflw
! --- cmy
      !type (topflw_type), dimension(ix, my_max) :: topflw
      !type (sfcflw_type), dimension(ix, my_max) :: sfcflw
! --- cmy

!  ---  variables are for both input and output:
      real (kind=kind_phys), intent(inout) :: cldcov(ix,lm+ltp, my_max)
      real (kind=kind_phys), intent(out) :: fluxr(ix, my_max, nfxr)

!! ---  optional outputs:
!     real (kind=kind_phys), dimension(ix,lm,nbdsw), optional,          &
!    &                       intent(out) :: htrswb
!     real (kind=kind_phys), dimension(ix,lm,nbdlw), optional,          &
!    &                       intent(out) :: htrlwb
!     real (kind=kind_phys), dimension(ix,lm), optional,                &
!    &                       intent(out) :: htrlw0
!     real (kind=kind_phys), dimension(ix,lm), optional,                &
!    &                       intent(out) :: htrsw0

      real (kind=kind_phys), dimension(ix,lm, my_max) :: htrlw0
      real (kind=kind_phys), dimension(ix,lm, my_max) :: htrsw0
      real (kind=kind_phys), dimension(ix,lm+ltp,5,  my_max)   :: phy_f3d
      logical uni_cloud,lmfshal,lmfdeep2
      integer :: async_id


!  ---  local variables: (horizontal dimensioned by im)
      real (kind=kind_phys), dimension(:,:,:,:), allocatable :: htrswb
      real (kind=kind_phys), dimension(:,:,:,:), allocatable :: htrlwb
      real (kind=kind_phys), dimension(nxptot,lm+1+ltp):: plvl, tlvl

      real (kind=kind_phys), dimension(nxptot,lm+ltp)  :: rhly, qstl, prslk1, tvly
      real (kind=kind_phys), dimension(nxptot,lm+ltp)  :: plyr, tlyr, qlyr, &
         olyr, vvel, clw, tem2da, tem2db
      !real (kind=kind_phys), dimension(ix,lm+ltp, my_max)  :: qst2, rhly2
      !real (kind=kind_phys), dimension(ix,lm+ltp, my_max)  :: es2, qs2
      real (kind=kind_phys), dimension(nxptot,lm+ltp)  :: qa
      real (kind=kind_phys), dimension(nxptot,lm+ltp)  :: cnvw1, cnvc1

      real (kind=kind_phys), dimension(nxptot) :: tsfa, cvt1, cvb1, tem1d,  &
         sfcemis, tsfg, tskn

      real (kind=kind_phys), dimension(nxptot,lm+ltp,nf_clds) :: clouds
      real (kind=kind_phys), dimension(nxptot,lm+ltp) :: cldfrc
      real (kind=kind_phys), dimension(nxptot,lm+ltp) :: cwp, rew, &
         cip, rei, crp, rer, csp, res 
      real (kind=kind_phys), dimension(nxptot,lm+ltp) :: gasvmr_co2
      real (kind=kind_phys), dimension(nf_vgas) :: gasvmr_other
      real (kind=kind_phys), dimension(nxptot,nf_albd) :: sfcalb
!     real (kind=kind_phys), dimension(im,       nspc1)   :: aerodp      ! optn for aod output
      real (kind=kind_phys), dimension(nxptot,lm+ltp,ntrac) :: tracer1

      !real (kind=kind_phys), dimension(ix,lm+ltp,nbdsw,nf_aesw, my_max) :: faersw
      !real (kind=kind_phys), dimension(ix,lm+ltp,nbdlw,nf_aelw, my_max) :: faerlw
      !real (kind=kind_phys), dimension(:,:,:), allocatable :: tauaer
      !real (kind=kind_phys), dimension(ix,lm+ltp, my_max) :: htswc
      !real (kind=kind_phys), dimension(ix,lm+ltp, my_max) :: htlwc

      !real (kind=kind_phys), dimension(ix,lm+ltp, my_max) :: gcice, grain, grime

!! ---  may be used for optional sw/lw outputs:
!!      take out "!!" as needed
      !real (kind=kind_phys), dimension(ix,lm+ltp, my_max)   :: htsw0
      !real (kind=kind_phys),    dimension(ix,lm+1+ltp, my_max, 4) :: fswprf
      !real (kind=kind_phys),    dimension(ix, my_max, 6)          :: scmpsw
      !real (kind=kind_phys), dimension(ix,lm+ltp,nbdsw, my_max) :: htswb

      !real (kind=kind_phys), dimension(ix,lm+ltp, my_max)   :: htlw0
      !type (proflw_type),    dimension(ix,lm+1+ltp, my_max) :: flwprf
      !real (kind=kind_phys), dimension(ix,lm+ltp,nbdlw, my_max) :: htlwb

      real (kind=kind_phys) :: raddt, es, qs, tem0d, cldsa(nxptot,5),qss

      integer :: i, j, k, k1, lv, itop, ibtc, nday(my_max), idxday(nxptot), kc,        &
         mbota(nxptot,3), mtopa(nxptot,3), lp1, nb, lmk, lmp, kd, lla, llb, &
         lya, lyb, kt, kb, jj
!effective radius for liquid, ice, snow, rain
      real (kind=kind_phys), dimension(:,:), allocatable :: plvl_im, tlvl_im, &
         plyr_im, tlyr_im, qlyr_im, olyr_im, rhly_im, qstl_im, vvel_im, clw_im, &
         prslk1_im, tem2da_im, tem2db_im, tvly_im, qst2_im, rhly2_im, es2_im, &
         qs2_im, qa_im, cnvw1_im, cnvc1_im, sfcalb_im, htswc_im, htlwc_im, &
         gcice_im, grain_im, grime_im, htsw0_im, htlw0_im, cldsa_im
      integer, dimension(:,:), allocatable :: mbota_im, mtopa_im
      type (proflw_type), dimension(:,:), allocatable :: flwprf_im
      real (kind=kind_phys), dimension(:,:,:,:), allocatable :: faersw_im, faerlw_im
      real (kind=kind_phys), dimension(:,:,:), allocatable :: clouds_im, gasvmr_im, &
         tracer1_im, htswb_im, htlwb_im
      ! for inlined qsatq
      real :: epsm2, qqq, fpvs_gpu
      character(len=4) :: myrank_str
      integer :: n, m
      integer, parameter :: nxpvs = 7501
      real :: c1xpvs,c2xpvs,tbpvs(nxpvs)
      common/fpvscom/ c1xpvs,c2xpvs,tbpvs(nxpvs)
      integer, parameter :: lwrad_block = 4
      integer, parameter :: swrad_block = 8
      integer :: swrad_smalljj, lwrad_smalljj
      ! for GPU register
      real :: tem2dar
      integer, dimension(my_max) :: offset_nday
      integer, dimension(ix*my_max) :: map_nday_ipt, map_nday_jj
      integer, dimension(nxptot) :: map_jj, map_i
      integer :: local_count
      integer :: max_nxjp_acc_length, jjoffset, jbs, jbe, nxjp_acc_length, jb, &
         max_nday_length, jbs_nday, jbe_nday, nday_length


      ! GPU: htrswb and htrlwb are computed but not used outside RRTMG. GPU version
      ! GPU: switch off to save memory.
      logical :: lhtrswb = .false., lhtrlwb = .false.

      ! GPU: variable name changed: CPU - faersw(:,:,:,1), GPU - tauae
      ! GPU: variable name changed: CPU - faersw(:,:,:,2), GPU - ssaae
      ! GPU: variable name changed: CPU - faersw(:,:,:,3), GPU - asyae
      ! GPU: variable name changed: CPU - faerlw(:,:,:,1), GPU - tauaer*
      ! GPU: variable name changed: CPU - faerlw(:,:,:,2), GPU - tauaer*
      ! GPU: variable name changed: CPU - faerlw(:,:,:,3), GPU - <dismissed>
      ! GPU: variable name changed: CPU - tem2da, GPU - <dismissed>
      ! GPU: variable name changed: CPU - tem2db, GPU - <dismissed>
      ! GPU: *: tauaer is computed from faerlw(:,:,:,1) and faerlw(:,:,:,2).
      ! GPU:    tauaer is computed at lwrad in the CPU version, but is computed
      ! GPU:    at setaer_gpu in this GPU verison and passed into subroutine 
      ! GPU:    lwrad_gpu.
      ! GPU: the following variable name changed when using icmphys = 15 or 16 
      ! GPU: (GCE microphysics)
      ! GPU:                        CPU - clouds(:,:,1), GPU - cldfrc
      ! GPU:                        CPU - clouds(:,:,2), GPU - cwp
      ! GPU:                        CPU - clouds(:,:,3), GPU - rew
      ! GPU:                        CPU - clouds(:,:,4), GPU - cip
      ! GPU:                        CPU - clouds(:,:,5), GPU - rei
      ! GPU:                        CPU - clouds(:,:,6), GPU - crp
      ! GPU:                        CPU - clouds(:,:,7), GPU - rer
      ! GPU:                        CPU - clouds(:,:,8), GPU - csp
      ! GPU:                        CPU - clouds(:,:,9), GPU - res
      ! GPU: if you use other microphysics option (icmphys != 15 or 16), modify
      ! GPU: the corresponding cloud information subroutine (progcld*). As 
      ! GPU: the subroutine progcld_gce_gpu, replace variable 'clouds' as the 
      ! GPU: above GPU variables in the subroutine passing argument list, 
      ! GPU: declaration executing and zone inside the subroutine.
      ! GPU: 
!

      swrad_smalljj = ceiling(float((my_max-1))/float(swrad_block))
      lwrad_smalljj = ceiling(float((my_max-1))/float(lwrad_block))
      !$acc data create(plvl, tlvl, plyr, tlyr, qlyr, olyr, &
      !$acc&     tsfa, tem1d, cldfrc, cwp, rew, cip, rei, crp, rer, csp, res, &
      !$acc&     sfcemis, tsfg, tskn, gasvmr_co2, gasvmr_other, sfcalb, &
      !$acc&     map_jj, map_i, &
      !$acc&     cldsa, nday, idxday, mbota, mtopa, map_nday_ipt, map_nday_jj, offset_nday) async(async_id)
      
      
      ! GPU: compute max_nxjp_acc_length using in lwrad_gpu
      max_nxjp_acc_length = 0
      do jb = 1, lwrad_block
         jjoffset = (my_max-1)*(jb-1)/lwrad_block
         jbs = jjoffset+1
         jbe = (my_max-1)*jb/lwrad_block
         nxjp_acc_length = nxjp_acc(jbe+1) - nxjp_acc(jbs)
         if (nxjp_acc_length .gt. max_nxjp_acc_length) max_nxjp_acc_length = nxjp_acc_length
      end do




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
      !$acc parallel loop collapse(2) private(n) async(async_id)
      do jj = 1, jlistnum
         do i = 1, ix
            if (i .le. myim(jj)) then
               n = i + nxjp_acc(jj) - 1
               map_jj(n) = jj
               map_i(n) = i
            end if
         end do
      end do
!  --- ...  compute cosin of zenith angle

      if (me == 0 .and. myrank ==0) print *,'### call coszmn'
      !call nvtxStartRange("coszmn")
      call coszmn_gpu                                                       &
      !  ---  inputs:
      &     ( xlon,sinlat,coslat, &
              solhr, myim, me, ix, map_jj, map_i, nxptot, async_id,                          &
      !  ---  outputs:
      &       coszen, coszdg                                             &
      &      )
      !call nvtxEndRange

!  --- ...  check for daytime points
      !$acc parallel loop gang private(local_count) async(async_id)
      do jj = 1, jlistnum
         local_count = 0
         !$acc loop vector reduction(+:local_count)
         do i = 1, myim(jj)
            if (coszen(i, jj) >= 0.0001) then
               local_count = local_count + 1
            endif
         enddo
         nday(jj) = local_count
      end do

      !$acc serial async(async_id)
      n = 0
      do jj = 1, jlistnum
         offset_nday(jj) = n
         n = n + nday(jj)
      end do
      offset_nday(my_max) = n
      !$acc end serial
      
      !$acc parallel loop gang private(n, m) async(async_id)
      do jj = 1, jlistnum
         n = offset_nday(jj)
         m = 0
         !$acc loop seq
         do i = 1, myim(jj)
            if (coszen(i, jj) >= 0.0001) then
               n = n + 1
               m = m + 1
               idxday(n) = i
               map_nday_ipt(n) = m
               map_nday_jj(n) = jj
            endif
         enddo
      end do
      
      !$acc update self(nday, offset_nday) async(async_id)
      !$acc wait(async_id)
      !if (me .eq. 0) then
         !$acc update self(map_jj, map_i, map_nday_ipt, map_nday_jj) async(async_id)
         !$acc wait(async_id)
      !end if
      ! GPU: compute max_nday_length using in swrad_gpu
      max_nday_length = 0
      do jb = 1, swrad_block
         jjoffset = (my_max-1)*(jb-1)/swrad_block
         jbs = jjoffset+1
         jbe = (my_max-1)*jb/swrad_block
         jbs_nday = offset_nday(jbs)+1
         jbe_nday = offset_nday(jbe+1)
         nday_length = jbe_nday - jbs_nday + 1
         if (nday_length .gt. max_nday_length) max_nday_length = nday_length
      end do
      

      do n = 1, nxptot
         jj = map_jj(n)
         i = map_i(n)
         if (i .eq. ipt) then
            if ( me == 0 .and. myrank == 0) then
               print *,'###################################################' 
               print *,'### In grrad start !! ###' 
               print *,'###################################################' 
               print *,'### ix=',ix,' im=',myim(jj),' lm=',lm,' me=',me
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
               print *,'### xlon(ipt)=',xlon(ipt, jj)
               print *,'### xlat(ipt)=',xlat(ipt, jj)
               print *,'### sinlat(ipt)=',sinlat(ipt, jj)
               print *,'### coslat(ipt)=',coslat(ipt, jj)
               print *,'###################################################' 
               print *,'### jdate(1-4)=',jdate(1),jdate(2),jdate(3),jdate(4) 
               print *,'### jdate(5-8)=',jdate(5),jdate(6),jdate(7),jdate(8)
               print *,'###################################################' 
               print *,'### cv(ipt)=',cv(ipt, jj)
               print *,'### cvb(ipt)=',cvb(ipt, jj)
               print *,'### cvt(ipt)=',cvt(ipt, jj)
               print *,'###################################################' 
               print *,'### prsi(ipt,lm+1)=',prsi(ipt,lm+1, my_max)
               print *,'### prsl(ipt,lm)=',prsl(ipt,lm, my_max)
               print *,'### prslk(ipt,lm)=',prslk(ipt,lm, my_max)
               print *,'### vvl(ipt,lm)=',vvl(ipt,lm, my_max)
               print *,'###################################################' 
               print *,'### tsfc(ipt)=',tsfc(ipt, jj)
               print *,'### tgrs(ipt,lm)=',tgrs(ipt,lm, my_max)
               print *,'### qgrs(ipt,lm)=',qgrs(ipt,lm, my_max)
               print *,'##############################################'
               print *,'### tracer(ipt,lm,3) =',tracer(ipt,lm,3, jj)
               print *,'###################################################' 
               print *,'### slmsk(ipt)=',slmsk(ipt, jj)
               print *,'### fice(ipt)=',fice(ipt, jj)
               print *,'### tisfc(ipt)=',tisfc(ipt, jj)
               print *,'###################################################' 
               print *,'### snowd(ipt)=',snowd(ipt, jj)
               print *,'### sncovr(ipt)=',sncovr(ipt, jj)
               print *,'### snoalb(ipt)=',snoalb(ipt, jj)
               print *,'###################################################' 
               print *,'### zorl(ipt)=',zorl(ipt, jj)
               print *,'### hprim(ipt)=',hprim(ipt, jj)
               print *,'###################################################' 
               print *,'### alvsf(ipt)=',alvsf(ipt, jj)
               print *,'### alnsf(ipt)=',alnsf(ipt, jj)
               print *,'### alvwf(ipt)=',alvwf(ipt, jj)
               print *,'### alnwf(ipt)=',alnwf(ipt, jj)
               print *,'### facsf(ipt)=',facsf(ipt, jj)
               print *,'### facwf(ipt)=',facwf(ipt, jj)
               print *,'###################################################' 
               print *,'### icsdsw(ipt)=',icsdsw(ipt, jj)
               print *,'### icsdlw(ipt)=',icsdlw(ipt, jj)
               print *,'###################################################' 
            endif
         end if
      end do
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
         !$acc parallel loop private(jj, i) async(async_id)
         do n = 1, nxptot
            jj = map_jj(n)
            i = map_i(n)
            tskn(n) = tsfc(i, jj)
            tsfg(n) = tsfc(i, jj)
         end do

      else                              ! use diff sfc skin-air/ground temp
         !$acc parallel loop private(jj, i) async(async_id)
         do n = 1, nxptot
            jj = map_jj(n)
            i = map_i(n)
            !!        tskn(i) = ta  (i)               ! not yet
            !!        tsfg(i) = tg  (i)               ! not yet
            tskn(n) = tsfc(i, jj)
            tsfg(n) = tsfc(i, jj)
         end do
      endif
      do n = 1, nxptot
         jj = map_jj(n)
         i = map_i(n)
         if (i .eq. ipt) then
            if (me == 0 .and. myrank ==0) print *,'### tskn(ipt)=',tskn(n),' ipt=',ipt
            if (me == 0 .and. myrank ==0) print *,'### tsfg(ipt)=',tsfg(n),' ipt=',ipt
         end if
         end do

!  --- ...  prepare atmospheric profiles for radiation input
!           convert pressure unit from cb to mb
      !$acc enter data create(tracer1, rhly, qstl, prslk1, tvly) async(async_id)
      !$acc parallel loop collapse(2) private(k1, jj, i, epsm2, qqq, qss) async(async_id)
      do k = 1, lm
         do n = 1, nxptot
            jj = map_jj(n)
            i = map_i(n)
            k1 = k + kd
            plvl(n,k1)   = 10.0 * prsi(i,k, jj)   ! cb (kpa) to mb (hpa)
            plyr(n,k1)   = 10.0 * prsl(i,k, jj)   ! cb (kpa) to mb (hpa)
            !         plvl(i,k1)   = 0.01 * prsi(i,k)   ! pa to mb (hpa)
            !         plyr(i,k1)   = 0.01 * prsl(i,k)   ! pa to mb (hpa)
            tlyr(n,k1)   = tgrs(i,k, jj)
            prslk1(n,k1) = prslk(i,k, jj)
            !cnvw1(n,k1)  = cnvw(i,k, jj)
            !cnvc1(n,k1)  = cnvc(i,k, jj)

!  --- ...  compute relative humidity
!         es  = min( prsl(i,k), 0.001 * fpvs( tgrs(i,k) ) )   ! fpvs in pa
!         qs  = max( qmin, con_eps * es / (prsl(i,k) + con_epsm1*es) )
!         rhly(i,k1) = max( 0.0, min( 1.0, max(qmin, qgrs(i,k))/qs ) )
!         qstl(i,k1) = qs
!--------------------------------------------------------------------------
            qlyr(n,k1) = max( qme6, qgrs(i,k, jj) )
            !call qsatq(1,tlyr(i,k1, jj),plyr(i,k1, jj),qss) !plyr in mb
            epsm2=0.622-1.
            qqq = min ( plyr(n,k1) , 0.01*fpvs_gpu(tlyr(n,k1),c1xpvs,c2xpvs,tbpvs) )
            qss = 0.622*qqq/(plyr(n,k1)+epsm2*qqq)
            ! end call qsatq (inlined)
            rhly(n,k1)= max( 0.0, min( 1.0, max(qmin, qlyr(n,k1))/qss ) )
            qstl(n,k1) = qss
!---------------------------------------------------------------------------
            !$acc loop seq
            do j = 1, ntrac
               tracer1(n,k1,j) = tracer(i,k,j, jj)
            enddo
         enddo
      end do
      
      do n = 1, nxptot
         jj = map_jj(n)
         i = map_i(n)
         if (i .eq. ipt) then
            if ( me == 0 .and. myrank == 0) then
               print *,'###################################################' 
               print *,'### prsl(ipt,1) =',prsl(ipt,1, jj)*10.,' in mb ###'
               print *,'### tgrs(ipt,1) =',tgrs(ipt,1, jj)
               print *,'### qgrs(ipt,1) =',qgrs(ipt,1, jj)
               print *,'### qlyr(ipt,1) =',qlyr(n,1)
               print *,'### rhly(ipt,1) =',rhly(n,1)
               print *,'### qstl(ipt,1) =',qstl(n,1)
               print *,'###################################################' 
               print *,'### prsl(ipt,lm) =',prsl(ipt,lm, jj)*10.,' in mb ###'
               print *,'### tgrs(ipt,lm) =',tgrs(ipt,lm, jj)
               print *,'### qgrs(ipt,lm) =',qgrs(ipt,lm, jj)
               print *,'### qlyr(ipt,lm) =',qlyr(n,lm)
               print *,'### rhly(ipt,lm) =',rhly(n,lm)
               print *,'### qstl(ipt,lm) =',qstl(n,lm)
               print *,'###################################################' 
            endif
         end if
      end do

      !$acc parallel loop private(jj, i) async(async_id)
      do n = 1, nxptot
         jj = map_jj(n)
         i = map_i(n)
         plvl(n,lp1+kd) = 10.0 * prsi(i,lp1, jj)  ! cb (kpa) to mb (hpa)
         !       plvl(i,lp1+kd) = 0.01 * prsi(i,lp1)  ! pa to mb (hpa)
      end do

      if ( lextop ) then                 ! values for extra top layer
         !$acc parallel loop private(jj, i) async(async_id)
         do n = 1, nxptot
            jj = map_jj(n)
            i = map_i(n)
            plvl(n,llb) = prsmin
            if ( plvl(n,lla) <= prsmin ) plvl(n,lla) = 2.0*prsmin
            plyr(n,lyb)   = 0.5 * plvl(n,lla)
            tlyr(n,lyb)   = tlyr(n,lya)
            prslk1(n,lyb) = (plyr(n,lyb)*0.001) ** rocp ! plyr in hpa

            rhly(n,lyb)   = rhly(n,lya)
            qstl(n,lyb)   = qstl(n,lya)
         end do

         !$acc parallel loop collapse(2) private(jj, i) async(async_id)
         do j = 1, ntrac
            do n = 1, nxptot
               jj = map_jj(n)
               i = map_i(n)
               if (i .le. myim(jj)) then
                  !  ---  note: may need to take care the top layer amount
                  tracer1(n,lyb,j) = tracer1(n,lya,j)
               end if
            enddo
         end do
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
         !$acc parallel loop gang collapse(2) private(jj, i) async(async_id)
         do k = 1, lmk
            do n = 1, nxptot
               jj = map_jj(n)
               i = map_i(n)
               olyr(n,k) = max( qmin, tracer1(n,k,ntoz) )
            enddo
         end do

         do n = 1, nxptot
            jj = map_jj(n)
            i = map_i(n)
            if (i .eq. ipt) then
               if (me == 0 .and. myrank ==0) print *, '### ntoz=',ntoz,' ipt=',ipt
               if (me == 0 .and. myrank ==0) print *, '### olyr(i,k)>= qmin, qmin=',qmin
               if (me == 0 .and. myrank ==0) print *, '### olyr(ipt,lm)=',olyr(n,lm)
            end if
         end do

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


      if ( myrank == 0 .and. me == 0) print *,'### call coszmn ok ! ###'
      do n = 1, nxptot
         jj = map_jj(n)
         i = map_i(n)
         if (i .eq. ipt) then
            if ( myrank == 0 .and. me == 0) then 
               print *,'solhr=',solhr,' ipt=',ipt                          
               print *,'xlon(ipt)=',xlon(ipt, jj)
               print *,'sinlat(ipt)=',sinlat(ipt, jj)
               print *,'coslat(ipt)=',coslat(ipt, jj)
               print *,'coszen(ipt)=',coszen(ipt, jj)
               print *,'coszdg(ipt)=',coszdg(ipt, jj)
            endif
         end if
      end do
!
!  --- ...  set up non-prognostic gas volume mixing ratioes

      if (me == 0 .and. myrank ==0) print *,'### call getgases'
      !call nvtxStartRange("getgases")
      call getgases_gpu                                                     &
      !  ---  inputs:
      &    ( plvl, xlon, xlat,                                           &
      &      myim, lmk, map_jj, map_i, nxptot, async_id,                                                    &
      !  ---  outputs:
      &      gasvmr_co2, gasvmr_other                                                      &
      &     )
      !call nvtxEndRange
      
      do n = 1, nxptot
         jj = map_jj(n)
         i = map_i(n)
         if (i .eq. ipt) then
            if (me == 0 .and. myrank ==0) then
               print *,'###################################################' 
               print *,'### call getgase ok !!'
               print *,'### plvl(ipt,1)=',plvl(n,1)
               print *,'### plvl(ipt,2)=',plvl(n,2)
               print *,'### plvl(ipt,lm)=',plvl(n,lm)
               print *,'### plvl(ipt,lm+1)=',plvl(n,lm+1)
               print *,'### xlon(ipt)=',xlon(ipt, jj)
               print *,'### xlat(ipt)=',xlat(ipt, jj)
               print *,'### lmk=',lmk,' im=',myim(jj),' ltp=',ltp
               print *,'### gasvmr(ipt,lm,1)_co2 =',gasvmr_co2(n,lm)
               print *,'### gasvmr(ipt,lm,2)_n2o =',gasvmr_other(2)
               print *,'### gasvmr(ipt,lm,3)_ch4 =',gasvmr_other(3)
               print *,'### gasvmr(ipt,lm,4)_o2  =',gasvmr_other(4)
               print *,'### gasvmr(ipt,lm,5)_co  =',gasvmr_other(5)
               print *,'### gasvmr(ipt,lm,6)_cf11=',gasvmr_other(6)
               print *,'### gasvmr(ipt,lm,7)_cf12=',gasvmr_other(7)
               print *,'### gasvmr(ipt,lm,8)_cf22=',gasvmr_other(8)
               print *,'### gasvmr(ipt,lm,9)_ccl4=',gasvmr_other(9)
               print *,'###################################################' 
            end if
         endif  
      end do
!
!  --- ...  get temperature at layer interface, and layer moisture

      if (ivflip == 0) then              ! input data from toa to sfc
         !$acc parallel loop collapse(2) private(jj, i) async(async_id)
         do k = 2, lmk
            do n = 1, nxptot
               jj = map_jj(n)
               i = map_i(n)
               tem2da(n,k) = log( plyr(n,k) )
               tem2db(n,k) = log( plvl(n,k) )
            enddo
         end do

         !$acc parallel loop private(jj, i) async(async_id)
         do n = 1, nxptot
            jj = map_jj(n)
            i = map_i(n)
            tem1d (n)   = qme6
            tem2da(n,1) = log( plyr(n,1) )
            tem2db(n,1) = 1.0
            tsfa  (n)   = tlyr(n,lmk)                  ! sfc layer air temp
            tlvl(n,1)   = tlyr(n,1)
            tlvl(n,lmp) = tskn(n)
         end do

         !$acc parallel loop private(k1, jj, i) async(async_id)
         do n = 1, nxptot
            jj = map_jj(n)
            i = map_i(n)
            !$acc loop seq
            do k = 1, lm
               k1 = k + kd
               qlyr(n,k1) = max( tem1d(n), qgrs(i,k, jj) )
               tem1d(n)   = min( qme5, qlyr(n,k1) )
               tvly(n,k1) = tgrs(i,k, jj) * (1.0 + con_fvirt*qlyr(n,k1))! virtual temp in k
            enddo
         end do

         !$acc parallel loop collapse(2) private(jj, i) async(async_id)
         do k = 2, lmk
            do n = 1, nxptot
               jj = map_jj(n)
               i = map_i(n)
               tlvl(n,k) = tlyr(n,k) + (tlyr(n,k-1) - tlyr(n,k))           &
               &                * (tem2db(n,k)   - tem2da(n,k))                   &
               &                / (tem2da(n,k-1) - tem2da(n,k))
            enddo
         end do
         
         if ( lextop ) then
            !$acc parallel loop private(k1, jj, i) async(async_id)
            do n = 1, nxptot
               jj = map_jj(n)
               i = map_i(n)
               qlyr(n,lyb) = qlyr(n,lya)
               tvly(n,lyb) = tvly(n,lya)
            end do
         endif
      else                               ! input data from sfc to toa

         !$acc parallel loop private(jj, i) async(async_id)
         do n = 1, nxptot
            jj = map_jj(n)
            i = map_i(n)
            tem1d (n)   = qme6
            tsfa  (n)   = tlyr(n,1)                    ! sfc layer air temp
            tlvl(n,1)   = tskn(n)
            tlvl(n,lmp) = tlyr(n,lmk)
         end do

         !$acc parallel loop private(jj, i) async(async_id)
         do n = 1, nxptot
            jj = map_jj(n)
            i = map_i(n)
            !$acc loop seq
            do k = lm, 1, -1
               qlyr(n,k) = max( tem1d(n), qgrs(i,k, jj) )
               tem1d(n)  = min( qme5, qlyr(n,k) )
               tvly(n,k) = tgrs(i,k, jj) * (1.0 + con_fvirt*qlyr(n,k)) ! virtual temp in k
            enddo
         end do

         !$acc parallel loop gang collapse(2) private(tem2dar, jj, i) async(async_id)
         do k = 1, lmk-1
            do n = 1, nxptot
               jj = map_jj(n)
               i = map_i(n)
               tem2dar = log(plyr(n,k))
               tlvl(n,k+1) = tlyr(n,k) + (tlyr(n,k+1) - tlyr(n,k))         &
               &                  * (log(plvl(n,k+1)) - tem2dar)                 &
               &                  / (log(plyr(n,k+1)) - tem2dar)
            enddo
         end do
         if ( lextop ) then
            !$acc parallel loop private(jj, i) async(async_id)
            do n = 1, nxptot
               jj = map_jj(n)
               i = map_i(n)
               qlyr(n,lyb) = qlyr(n,lya)
               tvly(n,lyb) = tvly(n,lya)
            end do
         endif
      endif                 ! end_if_ivflip
      
      !$acc parallel loop collapse(3) async(async_id)
      do jj = 1, jlistnum
         do k = 1, lm+1+ltp
            do i = 1, ix
               fusl(i, k, jj)       = 0.
               fdsl(i, k, jj)       = 0.
               fuslr(i, k, jj)      = 0.
               fdslr(i, k, jj)      = 0.
            end do
         end do
      end do

!      if (myrank == 0 ) print *,'nday=',nday

!  --- ...  setup aerosols property profile for radiation

      if (me == 0 .and. myrank ==0) print *,'### before setaer ###'
      do n = 1, nxptot
         jj = map_jj(n)
         i = map_i(n)
         if (i .eq. ipt) then
            if (myrank == 0 .and. me == 0)then
               print *,'###################################################' 
               print *,'### plvl(ipt,lm+1)=',plvl(n,lm+1)
               print *,'### plyr(ipt,lm)=',plyr(n,lm)
               print *,'### prslk1(ipt,lm)=',prslk1(n,lm)
               print *,'### tvly(ipt,lm)=',tvly(n,lm)
               print *,'### tlyr(ipt,lm)=',tlyr(n,lm)
               print *,'### qlyr(ipt,lm)=',qlyr(n,lm)
               print *,'### rhly(ipt,lm)=',rhly(n,lm)
               print *,'###################################################' 
            endif
         end if
      end do

      if (me == 0 .and. myrank ==0) then
         do n = 1, nxptot
            jj = map_jj(n)
            i = map_i(n)
            if (i .eq. ipt) then
               print *,'###################################################' 
               print *,'### after call setaer ###'
               print *,'###################################################' 
               print *,'### im=',myim(jj),' lmk=',lmk,' lmp=',lmp,' ipt=',ipt
               print *,'### xlon(ipt)=', xlon(ipt, jj)
               print *,'### xlat(ipt)=', xlat(ipt, jj)
               print *,'### slmsk(ipt)=',slmsk(ipt, jj)
               print *,'### plvl(ipt,lm+1)=',plvl(n,lm+1)
               print *,'### plyr(ipt,lm)=',plyr(n,lm)
               print *,'### prslk1(ipt,lm)=',prslk1(n,lm)
               print *,'### tvly(ipt,lm)=',tvly(n,lm)
               print *,'### rhly(ipt,lm)=',rhly(n,lm)
               print *,'### tracer1(ipt,lm,1)=',tracer1(n,lm,1)
               print *,'### tracer1(ipt,lm,2)=',tracer1(n,lm,2)
               print *,'### tracer1(ipt,lm,3)=',tracer1(n,lm,3)
               print *,'###################################################' 
               !print *,'### faersw(ipt,lm,1,1)=sw#1-opd =',tauae(n,lm,1)
               !print *,'### faersw(ipt,lm,1,2)=sw#1-ssa =',ssaae(n,lm,1)
               !print *,'### faersw(ipt,lm,1,3)=sw#1-asy =',asyae(n,lm,1)
               !print *,'### faerlw(ipt,lm,1,1)=lw#1-opd =',tauae(ipt,lm,1,jj)
               !print *,'### faerlw(ipt,lm,1,2)=lw#1-ssa =',ssaae(ipt,lm,1,jj)
               !print *,'### faerlw(ipt,lm,1,3)=lw#1-asy =',asyae(ipt,lm,1,jj)
               print *,'###################################################' 
            end if
         end do
      endif
      


!  --- ...  obtain cloud information for radiation calculations
      if (ntcw > 0) then                   ! prognostic cloud scheme
      !
         if (icmphys == 1) then           ! zhao/moorthi's prognostic cloud scheme
         !
            !$acc parallel loop collapse(2) private(jj, i) async(async_id)
            do k = 1, lmk
               do n = 1, nxptot
                  jj = map_jj(n)
                  i = map_i(n)
                  clw(n,k) = 0.0
                  !$acc loop seq
                  do j = 1, ncld
                  lv = ntcw + j - 1
!byl                 clw(i,k) = clw(i,k) + tracer1(i,k,lv)   ! cloud condensate amount
                     clw(n,k) = clw(n,k) + tracer1(n,k,lv) + cnvw(i,k, jj)  ! cloud condensate amount
                  enddo
               enddo
            end do
            !$acc parallel loop gang collapse(2) private(jj, i) async(async_id)
            do k = 1, lmk
               do n = 1, nxptot
                  jj = map_jj(n)
                  i = map_i(n)
                  if ( clw(n,k) < epsq ) clw(n,k) = 0.0
               enddo
            end do

            do jj = 1, jlistnum
               allocate(plvl_im(myim(jj), lm+1+ltp))
               allocate(plyr_im(myim(jj), lm+ltp))
               allocate(tlyr_im(myim(jj), lm+ltp))
               allocate(tvly_im(myim(jj), lm+ltp))
               allocate(qlyr_im(myim(jj), lm+ltp))
               allocate(qstl_im(myim(jj), lm+ltp))
               allocate(rhly_im(myim(jj), lm+ltp))
               allocate(clw_im(myim(jj), lm+ltp))
               allocate(cldsa_im(myim(jj), 5))
               allocate(mtopa_im(myim(jj), 3))
               allocate(mbota_im(myim(jj), 3))
               allocate(clouds_im(myim(jj),lm+ltp,nf_clds))
               
               do i = 1, myim(jj)
                  n = i + nxjp_acc(jj) - 1
                  do k = 1, 3
                     mtopa_im(i, k) = mtopa(n, k)
                     mbota_im(i, k) = mbota(n, k)
                  end do
                  do k = 1, 5
                     cldsa_im(i, k) = cldsa(n, k)
                  end do
                  do k = 1, lm+ltp
                     clw_im(i, k) = clw(n, k)
                     rhly_im(i, k) = rhly(n, k)
                     qstl_im(i, k) = qstl(n, k)
                     qlyr_im(i, k) = qlyr(n, k)
                     tvly_im(i, k) = tvly(n, k)
                     tlyr_im(i, k) = tlyr(n, k)
                     plyr_im(i, k) = plyr(n, k)
                     do j = 1, nf_clds
                        clouds_im(i, k, j) = clouds(n, k, j)
                     end do
                  end do
                  do k = 1, lm+1+ltp
                     plvl_im(i, k) = plvl(n, k)
                  end do
               end do
               if ( me == 0 .and. myrank == 0 )                              &
                  print *,'### call progcld1 -zhao/moorhi ###' 
               call progcld1_gpu                                                 &
               !  ---  inputs:
               &     ( plyr_im,plvl_im,tlyr_im,tvly_im,qlyr_im,qstl_im,rhly_im,clw_im,        &
               &       xlat(1:myim(jj), jj),xlon(1:myim(jj), jj),slmsk(1:myim(jj), jj),      &
               &       myim(jj), lmk, lmp, myrank,                                      &
               !  ---  outputs:
               &       clouds_im,cldsa_im,mtopa_im,mbota_im                                   &
               &      )
               do i = 1, myim(jj)
                  n = i + nxjp_acc(jj) - 1
                  do k = 1, 3
                     mtopa(n, k) = mtopa_im(i, k)
                     mbota(n, k) = mbota_im(i, k)
                  end do
                  do k = 1, 5
                     cldsa(n, k) = cldsa_im(i, k)
                  end do
                  do k = 1, lm+ltp
                     clw(n, k) = clw_im(i, k)
                     rhly(n, k) = rhly_im(i, k)
                     qstl(n, k) = qstl_im(i, k)
                     qlyr(n, k) = qlyr_im(i, k)
                     tvly(n, k) = tvly_im(i, k)
                     tlyr(n, k) = tlyr_im(i, k)
                     plyr(n, k) = plyr_im(i, k)
                     do j = 1, nf_clds
                        clouds(n, k, j) = clouds_im(i, k, j)
                     end do
                  end do
                  do k = 1, lm+1+ltp
                     plvl_im(i, k) = plvl(n, k)
                  end do
               end do
               deallocate(plyr_im,plvl_im,tlyr_im,tvly_im,qlyr_im,qstl_im,rhly_im,&
                  clw_im,cldsa_im,mtopa_im,mbota_im, clouds_im)
            end do


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
            do jj = 1, jlistnum
               do k = 1, lmk
                  do i = 1, myim(jj)
                     n = i + nxjp_acc(jj) - 1
                     clw(n,k) = 0.0
                  enddo

                  do j = 1, ncld
                     lv = ntcw + j - 1
                     do i = 1, myim(jj)
                        n = i + nxjp_acc(jj) - 1
                        clw(n,k) = clw(n,k) + tracer1(n,k,lv)   ! cloud condensate amount
                     enddo
                  enddo
               enddo
            end do
            
            do jj = 1, jlistnum
               do k = 1, lmk
                  do i = 1, myim(jj)
                     n = i + nxjp_acc(jj) - 1
                     if ( clw(n,k) < epsq ) clw(n,k) = 0.0
                  enddo
               enddo
            end do
!
            do jj = 1, jlistnum
               allocate(plvl_im(myim(jj), lm+1+ltp))
               allocate(plyr_im(myim(jj), lm+ltp))
               allocate(tlyr_im(myim(jj), lm+ltp))
               allocate(tvly_im(myim(jj), lm+ltp))
               allocate(qlyr_im(myim(jj), lm+ltp))
               allocate(qstl_im(myim(jj), lm+ltp))
               allocate(rhly_im(myim(jj), lm+ltp))
               allocate(clw_im(myim(jj), lm+ltp))
               allocate(cldsa_im(myim(jj), 5))
               allocate(mtopa_im(myim(jj), 3))
               allocate(mbota_im(myim(jj), 3))
               allocate(clouds_im(myim(jj),lm+ltp,nf_clds))
               do i = 1, myim(jj)
                  n = i + nxjp_acc(jj) - 1
                  do k = 1, 3
                     mtopa_im(i, k) = mtopa(n, k)
                     mbota_im(i, k) = mbota(n, k)
                  end do
                  do k = 1, 5
                     cldsa_im(i, k) = cldsa(n, k)
                  end do
                  do k = 1, lm+ltp
                     clw_im(i, k) = clw(n, k)
                     rhly_im(i, k) = rhly(n, k)
                     qstl_im(i, k) = qstl(n, k)
                     qlyr_im(i, k) = qlyr(n, k)
                     tvly_im(i, k) = tvly(n, k)
                     tlyr_im(i, k) = tlyr(n, k)
                     plyr_im(i, k) = plyr(n, k)
                     do j = 1, nf_clds
                        clouds_im(i, k, j) = clouds(n, k, j)
                     end do
                  end do
                  do k = 1, lm+1+ltp
                     plvl_im(i, k) = plvl(n, k)
                  end do
               end do
               if ( me == 0 .and. myrank == 0 )                               &
                  print *,'### call progcld3 -zhao/moorhi with PDF cloud###' 
               call progcld3_gpu                                                  &
               !  ---  inputs:
               &     ( plyr_im,plvl_im,tlyr_im,tvly_im,qlyr_im,qstl_im,rhly_im,clw_im,&
               cnvw(:, :, jj),cnvc(:, :, jj),          &
               &       xlat(1:myim(jj), jj),xlon(1:myim(jj), jj),slmsk(1:myim(jj), jj),    &
               &       myim(jj), lmk, lmp,                                              &
               &       deltaq(:, :, jj), sup,kdt,me,                                        &
               !  ---  outputs:
               &       clouds_im,cldsa_im,mtopa_im,mbota_im                                   &
               &      )
               do i = 1, myim(jj)
                  n = i + nxjp_acc(jj) - 1
                  do k = 1, 3
                     mtopa(n, k) = mtopa_im(i, k)
                     mbota(n, k) = mbota_im(i, k)
                  end do
                  do k = 1, 5
                     cldsa(n, k) = cldsa_im(i, k)
                  end do
                  do k = 1, lm+ltp
                     clw(n, k) = clw_im(i, k)
                     rhly(n, k) = rhly_im(i, k)
                     qstl(n, k) = qstl_im(i, k)
                     qlyr(n, k) = qlyr_im(i, k)
                     tvly(n, k) = tvly_im(i, k)
                     tlyr(n, k) = tlyr_im(i, k)
                     plyr(n, k) = plyr_im(i, k)
                     do j = 1, nf_clds
                        clouds(n, k, j) = clouds_im(i, k, j)
                     end do
                  end do
                  do k = 1, lm+1+ltp
                     plvl_im(i, k) = plvl(n, k)
                  end do
               end do
               deallocate(plyr_im,plvl_im,tlyr_im,tvly_im,qlyr_im,qstl_im,rhly_im,&
                  clw_im,cldsa_im,mtopa_im,mbota_im, clouds_im)
            end do
!
         elseif (icmphys == 6 .or. icmphys == 8) then    ! wsm6 & Thompson
            if ( me == 0 .and. myrank == 0 ) then
               if ( icmphys == 6 ) print *,'### call WSM6 cloud###' 
               if ( icmphys == 8 ) print *,'### call Thompson cloud###' 
            endif

            if (kdt == 1) then
               do jj = 1, jlistnum
                  phy_f3d(:,:,1, jj) = 10.
                  phy_f3d(:,:,2, jj) = 50.
                  phy_f3d(:,:,3, jj) = 250.
               end do
            endif
!
            do jj = 1, jlistnum
               allocate(plvl_im(myim(jj), lm+1+ltp))
               allocate(plyr_im(myim(jj), lm+ltp))
               allocate(tlyr_im(myim(jj), lm+ltp))
               allocate(qlyr_im(myim(jj), lm+ltp))
               allocate(qstl_im(myim(jj), lm+ltp))
               allocate(rhly_im(myim(jj), lm+ltp))
               allocate(cldsa_im(myim(jj), 5))
               allocate(mtopa_im(myim(jj), 3))
               allocate(mbota_im(myim(jj), 3))
               allocate(clouds_im(myim(jj),lm+ltp,nf_clds))
               allocate(tracer1_im(myim(jj),lm+ltp,ntrac))
               do i = 1, myim(jj)
                  n = i + nxjp_acc(jj) - 1
                  do k = 1, 3
                     mtopa_im(i, k) = mtopa(n, k)
                     mbota_im(i, k) = mbota(n, k)
                  end do
                  do k = 1, 5
                     cldsa_im(i, k) = cldsa(n, k)
                  end do
                  do k = 1, lm+ltp
                     rhly_im(i, k) = rhly(n, k)
                     qstl_im(i, k) = qstl(n, k)
                     qlyr_im(i, k) = qlyr(n, k)
                     tlyr_im(i, k) = tlyr(n, k)
                     plyr_im(i, k) = plyr(n, k)
                     do j = 1, nf_clds
                        clouds_im(i, k, j) = clouds(n, k, j)
                     end do
                     do j = 1, ntrac
                        tracer1_im(i, k, j) = tracer1(n, k, j)
                     end do
                  end do
                  do k = 1, lm+1+ltp
                     plvl_im(i, k) = plvl(n, k)
                  end do
               end do
               call progcld4_gpu                               &
                  !  --- inputs
                  ( plyr_im,plvl_im,tlyr_im,qlyr_im,qstl_im,rhly_im,tracer1_im,   &
                  xlat(1:myim(jj), jj),xlon(1:myim(jj), jj),slmsk(1:myim(jj), jj),     &
                  ntrac,ntcw,ntiw,ntrw,ntsw,ntgl,          &
                  myim(jj), lmk, lmp,                            &
                  uni_cloud,lmfshal,lmfdeep2,              &
                  cldcov(:, :, jj),phy_f3d(:,:,1, jj),                   &
                  phy_f3d(:,:,2, jj),phy_f3d(:,:,3, jj),           &
                  !   --- outputs:
                  clouds_im,cldsa_im,mtopa_im,mbota_im                 &
                  )
               do i = 1, myim(jj)
                  n = i + nxjp_acc(jj) - 1
                  do k = 1, 3
                     mtopa(n, k) = mtopa_im(i, k)
                     mbota(n, k) = mbota_im(i, k)
                  end do
                  do k = 1, 5
                     cldsa(n, k) = cldsa_im(i, k)
                  end do
                  do k = 1, lm+ltp
                     rhly(n, k) = rhly_im(i, k)
                     qstl(n, k) = qstl_im(i, k)
                     qlyr(n, k) = qlyr_im(i, k)
                     tvly(n, k) = tvly_im(i, k)
                     tlyr(n, k) = tlyr_im(i, k)
                     plyr(n, k) = plyr_im(i, k)
                     do j = 1, nf_clds
                        clouds(n, k, j) = clouds_im(i, k, j)
                     end do
                     do j = 1, ntrac
                        tracer1(n, k, j) = tracer1_im(i, k, j)
                     end do
                  end do
                  do k = 1, lm+1+ltp
                     plvl_im(i, k) = plvl(n, k)
                  end do
               end do
               deallocate(plyr_im,plvl_im,tlyr_im,qlyr_im,qstl_im,rhly_im,cldsa_im, &
                  mtopa_im,mbota_im, clouds_im, tracer1_im)
            end do

         elseif ( icmphys == 18 ) then   ! 2M Thompson
            if ( me == 0 .and. myrank == 0 )                               &
               print *,'### call New Thompson cloud ###'

            if (kdt == 1) then
               do jj = 1, jlistnum
                  phy_f3d(:,:,1, jj) = 10.
                  phy_f3d(:,:,2, jj) = 50.
                  phy_f3d(:,:,3, jj) = 250.
               end do
            endif

!         lwp_ex=0.0  !total liquid water path from explicit microphysics
!         iwp_ex=0.0  !total ice water path from explicit microphysics
!         lwp_fc=0.0  !total liquid water path from cloud fraction scheme
!         iwp_fc=0.0  !total ice water path from cloud fraction scheme
            do jj = 1, jlistnum
               allocate(plvl_im(myim(jj), lm+1+ltp))
               allocate(plyr_im(myim(jj), lm+ltp))
               allocate(tlyr_im(myim(jj), lm+ltp))
               allocate(qlyr_im(myim(jj), lm+ltp))
               allocate(qstl_im(myim(jj), lm+ltp))
               allocate(rhly_im(myim(jj), lm+ltp))
               allocate(cldsa_im(myim(jj), 5))
               allocate(mtopa_im(myim(jj), 3))
               allocate(mbota_im(myim(jj), 3))
               allocate(clouds_im(myim(jj),lm+ltp,nf_clds))
               allocate(tracer1_im(myim(jj),lm+ltp,ntrac))
               do i = 1, myim(jj)
                  n = i + nxjp_acc(jj) - 1
                  do k = 1, 3
                     mtopa_im(i, k) = mtopa(n, k)
                     mbota_im(i, k) = mbota(n, k)
                  end do
                  do k = 1, 5
                     cldsa_im(i, k) = cldsa(n, k)
                  end do
                  do k = 1, lm+ltp
                     rhly_im(i, k) = rhly(n, k)
                     qstl_im(i, k) = qstl(n, k)
                     qlyr_im(i, k) = qlyr(n, k)
                     tlyr_im(i, k) = tlyr(n, k)
                     plyr_im(i, k) = plyr(n, k)
                     do j = 1, nf_clds
                        clouds_im(i, k, j) = clouds(n, k, j)
                     end do
                     do j = 1, ntrac
                        tracer1_im(i, k, j) = tracer1(n, k, j)
                     end do
                  end do
                  do k = 1, lm+1+ltp
                     plvl_im(i, k) = plvl(n, k)
                  end do
               end do
               call progcld_thompson_gpu                                          &
                  !  --- inputs
                  ( plyr_im, plvl_im, tlyr_im, qlyr_im, qstl_im, rhly_im, tracer1_im,                &
                  xlat(1:myim(jj), jj), xlon(1:myim(jj), jj), slmsk(1:myim(jj), jj),     &
                  ntrac, ntcw, ntiw, ntrw, ntsw, ntgl,                        &
                  myim(jj), lmk, lmp,                                               &
                  uni_cloud, lmfshal, lmfdeep2, cldcov(:, :, jj),                       &
                  phy_f3d(:,:,1, jj), phy_f3d(:,:,2, jj), phy_f3d(:,:,3, jj),             &
                  !            lwp_ex, iwp_ex, lwp_fc, iwp_fc, dzlay,                      &
                  !            gridkm,                                                     &
                  !   --- outputs:
                  !            cld_frac, cld_lwp, cld_reliq, cld_iwp,                      &
                  !            cld_reice, cld_rwp, cld_rerain, cld_swp, cld_resnow)
                  clouds_im, cldsa_im, mtopa_im, mbota_im )
               do i = 1, myim(jj)
                  n = i + nxjp_acc(jj) - 1
                  do k = 1, 3
                     mtopa(n, k) = mtopa_im(i, k)
                     mbota(n, k) = mbota_im(i, k)
                  end do
                  do k = 1, 5
                     cldsa(n, k) = cldsa_im(i, k)
                  end do
                  do k = 1, lm+ltp
                     rhly(n, k) = rhly_im(i, k)
                     qstl(n, k) = qstl_im(i, k)
                     qlyr(n, k) = qlyr_im(i, k)
                     tvly(n, k) = tvly_im(i, k)
                     tlyr(n, k) = tlyr_im(i, k)
                     plyr(n, k) = plyr_im(i, k)
                     do j = 1, nf_clds
                        clouds(n, k, j) = clouds_im(i, k, j)
                     end do
                     do j = 1, ntrac
                        tracer1(n, k, j) = tracer1_im(i, k, j)
                     end do
                  end do
                  do k = 1, lm+1+ltp
                     plvl_im(i, k) = plvl(n, k)
                  end do
               end do
               deallocate(plyr_im,plvl_im,tlyr_im,qlyr_im,qstl_im,rhly_im,cldsa_im, &
                  mtopa_im,mbota_im, clouds_im, tracer1_im)
            end do

         elseif ( icmphys == 11 ) then   ! GFDL MP v1

            if ( me == 0 .and. myrank == 0 )                               &
               print *,'### call GFDL cloud ###'
            
            do jj = 1, jlistnum
               do k = 1, lm+ltp
                  do i = 1, myim(jj)
                     n = i + nxjp_acc(jj) - 1
                     clw(n, k) = 0.0
                  end do
               end do
            end do
            if ( .not. lgfdlmprad ) then
               do jj = 1, jlistnum
                  do k = 1, lmk
                     do i = 1, myim(jj)
                        do j = 1, ncld - 1
                           n = i + nxjp_acc(jj) - 1
                           lv = ntcw + j - 1
                           clw(n,k) = clw(n,k) + tracer1(n,k,lv)  ! cloud condensate amount
                        enddo
                        if ( clw(n,k) < epsq ) clw(n,k) = 0.0
                     enddo
                  enddo
               end do
            endif

            if (kdt == 1) then
               do jj = 1, jlistnum
                  phy_f3d(:,:,1, jj) = 10.
                  phy_f3d(:,:,2, jj) = 50.
                  phy_f3d(:,:,3, jj) = 250.
                  phy_f3d(:,:,4, jj) = 1000.
               end do
            endif
               
            if ( .not. lgfdlmprad ) then  ! no consistency between GFDLMP and radiation
               do jj = 1, jlistnum
                  allocate(plvl_im(myim(jj), lm+1+ltp))
                  allocate(plyr_im(myim(jj), lm+ltp))
                  allocate(tlyr_im(myim(jj), lm+ltp))
                  allocate(tvly_im(myim(jj), lm+ltp))
                  allocate(qlyr_im(myim(jj), lm+ltp))
                  allocate(qstl_im(myim(jj), lm+ltp))
                  allocate(rhly_im(myim(jj), lm+ltp))
                  allocate(clw_im(myim(jj), lm+ltp))
                  allocate(cldsa_im(myim(jj), 5))
                  allocate(mtopa_im(myim(jj), 3))
                  allocate(mbota_im(myim(jj), 3))
                  allocate(clouds_im(myim(jj),lm+ltp,nf_clds))
                  do i = 1, myim(jj)
                     n = i + nxjp_acc(jj) - 1
                     do k = 1, 3
                        mtopa_im(i, k) = mtopa(n, k)
                        mbota_im(i, k) = mbota(n, k)
                     end do
                     do k = 1, 5
                        cldsa_im(i, k) = cldsa(n, k)
                     end do
                     do k = 1, lm+ltp
                        clw_im(i, k) = clw(n, k)
                        rhly_im(i, k) = rhly(n, k)
                        qstl_im(i, k) = qstl(n, k)
                        qlyr_im(i, k) = qlyr(n, k)
                        tvly_im(i, k) = tvly(n, k)
                        tlyr_im(i, k) = tlyr(n, k)
                        plyr_im(i, k) = plyr(n, k)
                        do j = 1, nf_clds
                           clouds_im(i, k, j) = clouds(n, k, j)
                        end do
                     end do
                     do k = 1, lm+1+ltp
                        plvl_im(i, k) = plvl(n, k)
                     end do
                  end do
                  call progcld5_gpu                                                &
                     !    ---  inputs:
                     ( plyr_im,plvl_im,tlyr_im,tvly_im,qlyr_im,qstl_im,rhly_im,clw_im, &
                     cnvw(:, :, jj),cnvc(:, :, jj),        &
                     xlat(1:myim(jj), jj),xlon(1:myim(jj), jj),slmsk(1:myim(jj), jj),myim(jj),lmk,lmp,   &
                     cldcov(:, :, jj),                                                  &
                     !    ---  outputs:
                     clouds_im,cldsa_im,mtopa_im,mbota_im                                 &
                     ) 
                  do i = 1, myim(jj)
                     n = i + nxjp_acc(jj) - 1
                     do k = 1, 3
                        mtopa(n, k) = mtopa_im(i, k)
                        mbota(n, k) = mbota_im(i, k)
                     end do
                     do k = 1, 5
                        cldsa(n, k) = cldsa_im(i, k)
                     end do
                     do k = 1, lm+ltp
                        clw(n, k) = clw_im(i, k)
                        rhly(n, k) = rhly_im(i, k)
                        qstl(n, k) = qstl_im(i, k)
                        qlyr(n, k) = qlyr_im(i, k)
                        tvly(n, k) = tvly_im(i, k)
                        tlyr(n, k) = tlyr_im(i, k)
                        plyr(n, k) = plyr_im(i, k)
                        do j = 1, nf_clds
                           clouds(n, k, j) = clouds_im(i, k, j)
                        end do
                     end do
                     do k = 1, lm+1+ltp
                        plvl_im(i, k) = plvl(n, k)
                     end do
                  end do
                  deallocate(plyr_im,plvl_im,tlyr_im,tvly_im,qlyr_im,qstl_im,rhly_im,&
                     clw_im,cldsa_im,mtopa_im,mbota_im, clouds_im)
               end do
            else
               do jj = 1, jlistnum
                  allocate(plvl_im(myim(jj), lm+1+ltp))
                  allocate(plyr_im(myim(jj), lm+ltp))
                  allocate(tlyr_im(myim(jj), lm+ltp))
                  allocate(tvly_im(myim(jj), lm+ltp))
                  allocate(qlyr_im(myim(jj), lm+ltp))
                  allocate(qstl_im(myim(jj), lm+ltp))
                  allocate(rhly_im(myim(jj), lm+ltp))
                  allocate(cldsa_im(myim(jj), 5))
                  allocate(mtopa_im(myim(jj), 3))
                  allocate(mbota_im(myim(jj), 3))
                  allocate(clouds_im(myim(jj),lm+ltp,nf_clds))
                  allocate(tracer1_im(myim(jj),lm+ltp,ntrac))
                  do i = 1, myim(jj)
                     n = i + nxjp_acc(jj) - 1
                     do k = 1, 3
                        mtopa_im(i, k) = mtopa(n, k)
                        mbota_im(i, k) = mbota(n, k)
                     end do
                     do k = 1, 5
                        cldsa_im(i, k) = cldsa(n, k)
                     end do
                     do k = 1, lm+ltp
                        rhly_im(i, k) = rhly(n, k)
                        qstl_im(i, k) = qstl(n, k)
                        qlyr_im(i, k) = qlyr(n, k)
                        tvly_im(i, k) = tvly(n, k)
                        tlyr_im(i, k) = tlyr(n, k)
                        plyr_im(i, k) = plyr(n, k)
                        do j = 1, nf_clds
                           clouds_im(i, k, j) = clouds(n, k, j)
                        end do
                        do j = 1, ntrac
                           tracer1_im(i, k, j) = tracer1(n, k, j)
                        end do
                     end do
                     do k = 1, lm+1+ltp
                        plvl_im(i, k) = plvl(n, k)
                     end do
                  end do
                  call progcld5o_gpu                                               &
                     !    ---  inputs:
                     ( plyr_im,plvl_im,tlyr_im,tvly_im,qlyr_im,qstl_im,rhly_im,tracer1_im,              &
                     xlat(1:myim(jj), jj),xlon(1:myim(jj), jj),slmsk(1:myim(jj), jj),        &
                     ntrac,ntcw,ntiw,ntrw,ntsw,ntgl,cldcov(:, :, jj),                   &
                     phy_f3d(:,:,1, jj),phy_f3d(:,:,2, jj),phy_f3d(:,:,3, jj),            &
                     phy_f3d(:,:,4, jj),effr_in,                                  &
                     myim(jj),lmk,lmp,                                              &
                     !    ---  outputs:
                     clouds_im,cldsa_im,mtopa_im,mbota_im                                 &
                     ) 
                  do i = 1, myim(jj)
                     n = i + nxjp_acc(jj) - 1
                     do k = 1, 3
                        mtopa(n, k) = mtopa_im(i, k)
                        mbota(n, k) = mbota_im(i, k)
                     end do
                     do k = 1, 5
                        cldsa(n, k) = cldsa_im(i, k)
                     end do
                     do k = 1, lm+ltp
                        rhly(n, k) = rhly_im(i, k)
                        qstl(n, k) = qstl_im(i, k)
                        qlyr(n, k) = qlyr_im(i, k)
                        tvly(n, k) = tvly_im(i, k)
                        tlyr(n, k) = tlyr_im(i, k)
                        plyr(n, k) = plyr_im(i, k)
                        do j = 1, nf_clds
                           clouds(n, k, j) = clouds_im(i, k, j)
                        end do
                        do j = 1, ntrac
                           tracer1(n, k, j) = tracer1_im(i, k, j)
                        end do
                     end do
                     do k = 1, lm+1+ltp
                        plvl_im(i, k) = plvl(n, k)
                     end do
                  end do
                  deallocate(plyr_im,plvl_im,tlyr_im,tvly_im,qlyr_im,qstl_im,rhly_im, &
                     cldsa_im,mtopa_im,mbota_im, clouds_im, tracer1_im)
               end do
!           endif
            endif
                  
         elseif ( icmphys == 12 .or. icmphys == 13 ) then   ! GFDL MP v2 / v3
            
            if ( me == 0 .and. myrank == 0 .and. icmphys == 12 )           &
               print *,'### call GFDL v2 cloud ###'
            if ( me == 0 .and. myrank == 0 .and. icmphys == 13 )           &
               print *,'### call GFDL v3 cloud ###'
               
            do jj = 1, jlistnum
               do k = 1, lm+ltp
                  do i = 1, myim(jj)
                     n = i + nxjp_acc(jj) - 1
                     qa(n, k) = 0.  !aerosol mixing ratio (kg/kg)
                  end do
               end do
            end do
               
            do jj = 1, jlistnum
               allocate(plvl_im(myim(jj), lm+1+ltp))
               allocate(plyr_im(myim(jj), lm+ltp))
               allocate(tlyr_im(myim(jj), lm+ltp))
               allocate(tvly_im(myim(jj), lm+ltp))
               allocate(qlyr_im(myim(jj), lm+ltp))
               allocate(qstl_im(myim(jj), lm+ltp))
               allocate(rhly_im(myim(jj), lm+ltp))
               allocate(cnvc1_im(myim(jj), lm+ltp))
               allocate(cnvw1_im(myim(jj), lm+ltp))
               allocate(cldsa_im(myim(jj), 5))
               allocate(mtopa_im(myim(jj), 3))
               allocate(mbota_im(myim(jj), 3))
               allocate(clouds_im(myim(jj),lm+ltp,nf_clds))
               allocate(tracer1_im(myim(jj),lm+ltp,ntrac))
               do i = 1, myim(jj)
                  n = i + nxjp_acc(jj) - 1
                  do k = 1, 3
                     mtopa_im(i, k) = mtopa(n, k)
                     mbota_im(i, k) = mbota(n, k)
                  end do
                  do k = 1, 5
                     cldsa_im(i, k) = cldsa(n, k)
                  end do
                  do k = 1, lm+ltp
                     cnvc1_im(i, k) = cnvc1(n, k)
                     cnvw1_im(i, k) = cnvw1(n, k)
                     rhly_im(i, k) = rhly(n, k)
                     qstl_im(i, k) = qstl(n, k)
                     qlyr_im(i, k) = qlyr(n, k)
                     tvly_im(i, k) = tvly(n, k)
                     tlyr_im(i, k) = tlyr(n, k)
                     plyr_im(i, k) = plyr(n, k)
                     do j = 1, nf_clds
                        clouds_im(i, k, j) = clouds(n, k, j)
                     end do
                     do j = 1, ntrac
                        tracer1_im(i, k, j) = tracer1(n, k, j)
                     end do
                  end do
                  do k = 1, lm+1+ltp
                     plvl_im(i, k) = plvl(n, k)
                  end do
               end do
               call progcld6_gpu                                                  &
                  !    ---  inputs:
                  ( plyr_im,plvl_im,tlyr_im,tvly_im,qlyr_im,qstl_im,rhly_im,cnvw1_im,cnvc1_im,          &
                  tracer1_im(:,:,ntcw),tracer1_im(:,:,ntrw),tracer1_im(:,:,ntiw),   &
                  tracer1_im(:,:,ntsw),tracer1_im(:,:,ntgl),                     &
                  cldcov(:, :, jj),slmsk(1:myim(jj), jj),                                            &
                  phy_f3d(:,:,1, jj),phy_f3d(:,:,2, jj),phy_f3d(:,:,3, jj),            &
                  phy_f3d(:,:,4, jj),effr_in,                                  &
                  xlat(1:myim(jj), jj),xlon(1:myim(jj), jj),myim(jj),lmk,lmp,       &
                  !    ---  outputs:
                  clouds_im,cldsa_im,mtopa_im,mbota_im                                 &
                  ) 
               do i = 1, myim(jj)
                  n = i + nxjp_acc(jj) - 1
                  do k = 1, 3
                     mtopa(n, k) = mtopa_im(i, k)
                     mbota(n, k) = mbota_im(i, k)
                  end do
                  do k = 1, 5
                     cldsa(n, k) = cldsa_im(i, k)
                  end do
                  do k = 1, lm+ltp
                     cnvc1(n, k) = cnvc1_im(i, k)
                     cnvw1(n, k) = cnvw1_im(i, k)
                     rhly(n, k) = rhly_im(i, k)
                     qstl(n, k) = qstl_im(i, k)
                     qlyr(n, k) = qlyr_im(i, k)
                     tvly(n, k) = tvly_im(i, k)
                     tlyr(n, k) = tlyr_im(i, k)
                     plyr(n, k) = plyr_im(i, k)
                     do j = 1, nf_clds
                        clouds(n, k, j) = clouds_im(i, k, j)
                     end do
                     do j = 1, ntrac
                        tracer1(n, k, j) = tracer1_im(i, k, j)
                     end do
                  end do
                  do k = 1, lm+1+ltp
                     plvl_im(i, k) = plvl(n, k)
                  end do
               end do
               deallocate(plyr_im,plvl_im,tlyr_im,tvly_im,qlyr_im,qstl_im,rhly_im, &
                  cnvc1_im,cnvw1_im,cldsa_im,mtopa_im,mbota_im, clouds_im, tracer1_im)
            end do

         elseif ( icmphys == 15 .or. icmphys == 16 ) then   ! Goddard (GCE)
            
            if (kdt == 1) then
               !$acc parallel loop gang collapse(3) async(async_id)
               do jj = 1, jlistnum
                  do k = 1, lm+ltp
                     do i = 1, ix
                        phy_f3d(i,k,1, jj) = 10.
                        phy_f3d(i,k,2, jj) = 50.
                        phy_f3d(i,k,3, jj) = 250.
                        phy_f3d(i,k,4, jj) = 1000.
                     end do
                  end do
               end do
            endif
               
            if ( me == 0 .and. myrank == 0 )                               &
               print *,'### call Goddard (GCE) cloud ###'
            !call nvtxStartRange("progcld_gce")
            ! GPU: the following variable name changed when using icmphys = 15 or 16 
            ! GPU: (GCE microphysics)
            ! GPU:                        CPU - clouds(:,:,1), GPU - cldfrc
            ! GPU:                        CPU - clouds(:,:,2), GPU - cwp
            ! GPU:                        CPU - clouds(:,:,3), GPU - rew
            ! GPU:                        CPU - clouds(:,:,4), GPU - cip
            ! GPU:                        CPU - clouds(:,:,5), GPU - rei
            ! GPU:                        CPU - clouds(:,:,6), GPU - crp
            ! GPU:                        CPU - clouds(:,:,7), GPU - rer
            ! GPU:                        CPU - clouds(:,:,8), GPU - csp
            ! GPU:                        CPU - clouds(:,:,9), GPU - res
            call progcld_gce_gpu                                             &
               !    ---  inputs:
               ( plyr, plvl, tlyr, tvly, qlyr, qstl, rhly, tracer1,       &
               xlat, xlon, slmsk, ntrac,                                &
               phy_f3d, effr_in,                                 &
               myim, ix, lmk, lmp, lmfshal, lmfdeep2, map_jj, map_i, nxptot, async_id,                         &
               !    ---  outputs:
               cldfrc, cwp, rew, cip, rei, crp, rer, csp, res, cldsa, mtopa, mbota, cldcov     &
               )
            !call nvtxEndRange
         endif                            ! end if_icmphys
      else                                 ! diagnostic cloud scheme

         !$acc parallel loop private(jj, i) async(async_id)
         do n = 1, nxptot
            jj = map_jj(n)
            i = map_i(n)
            cvt1(n) = 10.0 * cvt(i, jj)
            cvb1(n) = 10.0 * cvb(i, jj)
         end do

         !$acc parallel loop gang collapse(2) private(k1, jj, i) async(async_id)
         do k = 1, lm
            do n = 1, nxptot
               jj = map_jj(n)
               i = map_i(n)
               k1 = k + kd
               vvel(n,k1) = 10.0 * vvl(i,k, jj)
            enddo
         end do

         if ( lextop ) then
            !$acc parallel loop private(jj, i) async(async_id)
            do n = 1, nxptot
               jj = map_jj(n)
               i = map_i(n)
               vvel(n,lyb) = vvel(n,lya)
            end do
         endif

!  ---  compute diagnostic cloud related quantities

         if (me == 0 .and. myrank ==0)print *,'### call diagcld1 ###' 
            
         do jj = 1, jlistnum
            allocate(plvl_im(myim(jj), lm+1+ltp))
            allocate(plyr_im(myim(jj), lm+ltp))
            allocate(tlyr_im(myim(jj), lm+ltp))
            allocate(rhly_im(myim(jj), lm+ltp))
            allocate(vvel_im(myim(jj), lm+ltp))
            allocate(cldsa_im(myim(jj), 5))
            allocate(mtopa_im(myim(jj), 3))
            allocate(mbota_im(myim(jj), 3))
            allocate(clouds_im(myim(jj),lm+ltp,nf_clds))
            
            do i = 1, myim(jj)
               n = i + nxjp_acc(jj) - 1
               do k = 1, 3
                  mtopa_im(i, k) = mtopa(n, k)
                  mbota_im(i, k) = mbota(n, k)
               end do
               do k = 1, 5
                  cldsa_im(i, k) = cldsa(n, k)
               end do
               do k = 1, lm+ltp
                  vvel_im(i, k) = vvel(n, k)
                  rhly_im(i, k) = rhly(n, k)
                  tlyr_im(i, k) = tlyr(n, k)
                  plyr_im(i, k) = plyr(n, k)
                  do j = 1, nf_clds
                     clouds_im(i, k, j) = clouds(n, k, j)
                  end do
               end do
               do k = 1, lm+1+ltp
                  plvl_im(i, k) = plvl(n, k)
               end do
            end do
            if (myrank .eq. 0) write(*,*) "run diagcld1"
            call diagcld1_gpu                                                   &
               !  ---  inputs:
               &     ( plyr_im,plvl_im,tlyr_im,rhly_im,vvel_im,cv(1:myim(jj), jj), &
               cvt1(nxjp_acc(jj):nxjp_acc(jj+1)-1),cvb1(nxjp_acc(jj):nxjp_acc(jj+1)-1),                     &
               &       xlat(1:myim(jj), jj),xlon(1:myim(jj), jj),slmsk(1:myim(jj), jj),      &
               &       myim(jj), lmk, lmp,                                              &
               !  ---  outputs:
               &       clouds_im,cldsa_im,mtopa_im,mbota_im                                   &
               &      )
            do i = 1, myim(jj)
               n = i + nxjp_acc(jj) - 1
               do k = 1, 3
                  mtopa(n, k) = mtopa_im(i, k)
                  mbota(n, k) = mbota_im(i, k)
               end do
               do k = 1, 5
                  cldsa(n, k) = cldsa_im(i, k)
               end do
               do k = 1, lm+ltp
                  vvel(n, k) = vvel_im(i, k)
                  rhly(n, k) = rhly_im(i, k)
                  tlyr(n, k) = tlyr_im(i, k)
                  plyr(n, k) = plyr_im(i, k)
                  do j = 1, nf_clds
                     clouds(n, k, j) = clouds_im(i, k, j)
                  end do
               end do
               do k = 1, lm+1+ltp
                  plvl_im(i, k) = plvl(n, k)
               end do
            end do
            deallocate(plyr_im,plvl_im,tlyr_im,rhly_im,vvel_im,cldsa_im,mtopa_im, &
               mbota_im, clouds_im)
         end do

      endif                                ! end_if_ntcw

      do n = 1, nxptot
         jj = map_jj(n)
         i = map_i(n)
         if (i .eq. ipt) then
            if (me == 0 .and. myrank ==0) then
               print *,'###################################################' 
               print *,'###  after diagcld1'
               print *,'###################################################' 
               print *,'### im=',myim(jj),' lmk=',lmk,' lmp=',lmp
               print *,'### ntcw=',ntcw
               print *,'### ncld=',ncld
               print *,'###################################################' 
               print *,'### xlon(ipt)=', xlon(ipt, jj)
               print *,'### xlat(ipt)=', xlat(ipt, jj)
               print *,'### slmsk(ipt)=',slmsk(ipt, jj)
               print *,'###################################################' 
               print *,'### plvl(ipt,lm)=',plvl(n,lm)
               print *,'### plyr(ipt,lm)=',plyr(n,lm)
               print *,'###################################################' 
               print *,'### tlyr(ipt,lm)=',tlyr(n,lm)
               print *,'### tvly(ipt,lm)=',tvly(n,lm)
               print *,'###################################################' 
               print *,'### qlyr(ipt,lm)=',qlyr(n,lm)
               print *,'### qstl(ipt,lm)=',qstl(n,lm)
               print *,'### rhly(ipt,lm)=',rhly(n,lm)
               print *,'### clw(ipt,lm) =',clw(n,lm)
               print *,'###################################################' 
               print *,'### tracer1(ipt,lm,1)=',tracer1(n,lm,1)
               print *,'### tracer1(ipt,lm,2)=',tracer1(n,lm,2)
               print *,'### tracer1(ipt,lm,3)=',tracer1(n,lm,3)
               print *,'###################################################' 
               print *,'### clouds(ipt,lm,1)-total cloud fraction =',clouds(n,lm,1)
               print *,'### clouds(ipt,lm,2)-liq water path       =',clouds(n,lm,2)
               print *,'### clouds(ipt,lm,3)-liq effective radius =',clouds(n,lm,3)
               print *,'### clouds(ipt,lm,4)-ice water path       =',clouds(n,lm,4)
               print *,'### clouds(ipt,lm,5)-ice effective radius =',clouds(n,lm,5)
               print *,'### clouds(ipt,lm,6)-rain water path      =',clouds(n,lm,6)
               print *,'### clouds(ipt,lm,7)-rain effective radius=',clouds(n,lm,7)
               print *,'### clouds(ipt,lm,8)-snow water path      =',clouds(n,lm,8)
               print *,'### clouds(ipt,lm,9)-snow effective radius=',clouds(n,lm,9)
               print *,'###################################################' 
               print *,'### cldsa(ipt,1)-low clouds fraction=',cldsa(n,1)
               print *,'### cldsa(ipt,2)-mid clouds fraction=',cldsa(n,2)
               print *,'### cldsa(ipt,3)-hig clouds fraction=',cldsa(n,3)
               print *,'### cldsa(ipt,4)-tot clouds fraction=',cldsa(n,4)
               print *,'### cldsa(ipt,5)-bl  clouds fraction=',cldsa(n,5)
               print *,'###################################################' 
               print *,'### mtopa(ipt,1)-low clouds top =',mtopa(n,1)
               print *,'### mtopa(ipt,2)-mid clouds top =',mtopa(n,2)
               print *,'### mtopa(ipt,3)-hig clouds top =',mtopa(n,3)
               print *,'###################################################' 
               print *,'### mbota(ipt,1)-low clouds bottom =',mbota(n,1)
               print *,'### mbota(ipt,2)-mid clouds bottom =',mbota(n,2)
               print *,'### mbota(ipt,3)-hig clouds bottom =',mbota(n,3)
               print *,'###################################################' 
            end if
         endif
      end do


!  --- ...  start radiation calculations 
!           remember to set heating rate unit to k/sec!
      if (lsswr) then

!  ---  setup surface albedo for sw radiation, incl xw (nov04) sea-ice

         if (me == 0 .and. myrank ==0)print *,'### call setalb ###' 
         call setalb_gpu                                                     &
            !  ---  inputs:
            &     ( slmsk,snowd,sncovr,snoalb,&
            zorl,coszen,tsfg, &
            tsfa,hprim,     &
            &       alvsf,alnsf,alvwf,alnwf,&
            facsf,facwf,fice,tisfc,            &
            &       myim, ix, map_jj, map_i, nxptot, async_id,                                                       &
            !  ---  outputs:
            &       sfcalb                                                     &
            &     )
         !call nvtxEndRange

!  --- lu [+4l]: derive sfalb from vis- and nir- diffuse surface albedo
         !$acc parallel loop private(jj, i) async(async_id)
         do n = 1, nxptot
            jj = map_jj(n)
            i = map_i(n)
            sfalb(i, jj) = max(0.01, 0.5 * (sfcalb(n,2) + sfcalb(n,4)))
         end do
      
         do n = 1, nxptot
            jj = map_jj(n)
            i = map_i(n)
            if (i .eq. ipt) then
               if (me == 0 .and. myrank ==0) then
                  print *,'###################################################' 
                  print *,'### call setalb ok !!'
                  print *,'###################################################' 
                  print *,'### slmsk(ipt)=',slmsk(ipt, jj)
                  print *,'### fice(ipt)=',fice(ipt, jj)
                  print *,'### tisfc(ipt)=',tisfc(ipt, jj)
                  print *,'###################################################' 
                  print *,'### snowd(ipt)=',snowd(ipt, jj)
                  print *,'### sncovr(ipt)=',sncovr(ipt, jj)
                  print *,'### snoalb(ipt)=',snoalb(ipt, jj)
                  print *,'###################################################' 
                  print *,'### zorl(ipt)=',zorl(ipt, jj)
                  print *,'### hprim(ipt)=',hprim(ipt, jj)
                  print *,'###################################################' 
                  print *,'### alvsf(ipt)=',alvsf(ipt, jj)
                  print *,'### alnsf(ipt)=',alnsf(ipt, jj)
                  print *,'### alvwf(ipt)=',alvwf(ipt, jj)
                  print *,'### alnwf(ipt)=',alnwf(ipt, jj)
                  print *,'### facsf(ipt)=',facsf(ipt, jj)
                  print *,'### facwf(ipt)=',facwf(ipt, jj)
                  print *,'###################################################' 
                  print *,'### coszen(ipt)=',coszen(ipt, jj)
                  print *,'### tsfg(ipt)=',tsfg(n)
                  print *,'### tsfa(ipt)=',tsfa(n)
                  print *,'###################################################' 
                  print *,'### sfcalb(ipt,1)-near ir direct beam albedo'
                  print *,'### sfcalb(ipt,2)-near ir diffused beam albedo'
                  print *,'### sfcalb(ipt,3)-uv+vis direct beam albedo'
                  print *,'### sfcalb(ipt,4)-uv+vis diffused beam albedo'
                  print *,'###################################################' 
                  print *,'### sfcalb(ipt,1)=',sfcalb(n,1)
                  print *,'### sfcalb(ipt,2)=',sfcalb(n,2)
                  print *,'### sfcalb(ipt,3)=',sfcalb(n,3)
                  print *,'### sfcalb(ipt,4)=',sfcalb(n,4)
                  print *,'###################################################' 
                  print *,'# sfalb(i)-average diffused beam albedo for sw&lw'
                  print *,'# sfalb(i)=max(0.01,0.5*(sfcalb(i,2)+sfcalb(i,4)))'
                  print *,'###################################################' 
                  print *,'### sfalb(ipt)=',sfalb(ipt, jj)
                  print *,'###################################################' 
               endif
            end if
         end do
         do n = 1, nxptot
            jj = map_jj(n)
            i = map_i(n)
            if (i .eq. 1) then
               if (nday(jj) > 0) then
                  if (me == 0 .and. myrank ==0) print *,' #### call swrad ####'
                  if (me == 0 .and. myrank ==0) print *,' #### nday=',nday(jj)
               end if
            end if
         end do

         if (lhtrswb) then
            allocate(htrswb(ix,lm,nbdsw, my_max))
            !$acc enter data create(htrswb) async(async_id)
         end if

         if (lhtrswb) then
            !$acc parallel loop collapse(3) private(k1, jj, i) async(async_id)
            do k = 1, lm
               do j = 1, nbdsw
                  do n = 1, nxptot
                     jj = map_jj(n)
                     i = map_i(n)
                     k1 = k + kd
                     !htrswb(i,k,j, jj) = htswb(i,k1,j, jj)
                     htrswb(i,k,j, jj) = 0
                  enddo
               end do
            enddo
         end if

!         if ( present(htrswb) .and. present(htrsw0) ) then
!         if ( present(htrsw0) .and. present(fswprf) ) then
         !call nvtxStartRange("swrad")
         ! GPU: variable name changed: CPU - htswb, GPU - htrswb
         ! GPU: variable name changed: CPU - htswc, GPU - htrsw
         ! GPU: variable name changed: CPU - htsw0, GPU - htrsw0
         ! GPU: variable name changed: CPU - topfsw%upfxc, GPU - fluxr(i, jj, 2 )
         ! GPU: variable name changed: CPU - topfsw%dnfxc, GPU - fluxr(i, jj, 1 )
         ! GPU: variable name changed: CPU - topfsw%upfx0, GPU - fluxr(i, jj, 22)
         ! GPU: variable name changed: CPU - sfcfsw%upfxc, GPU - fluxr(i, jj, 5 )
         ! GPU: variable name changed: CPU - sfcfsw%dnfxc, GPU - fluxr(i, jj, 4 )
         ! GPU: variable name changed: CPU - sfcfsw%upfx0, GPU - fluxr(i, jj, 25)
         ! GPU: variable name changed: CPU - sfcfsw%dnfx0, GPU - fluxr(i, jj, 24)
         ! GPU: variable name changed: CPU - fswprf%upfxc, GPU - fusl
         ! GPU: variable name changed: CPU - fswprf%dnfxc, GPU - fdsl
         ! GPU: variable name changed: CPU - fswprf%upfx0, GPU - fuslr
         ! GPU: variable name changed: CPU - fswprf%dnfx0, GPU - fdslr
         ! GPU: variable name changed: CPU - scmpsw%uvbfc, GPU - fluxr(i, jj, 28)
         ! GPU: variable name changed: CPU - scmpsw%uvbf0, GPU - fluxr(i, jj, 29)
         ! GPU: variable name changed: CPU - scmpsw%nirbm, GPU - fluxr(i, jj, 32)
         ! GPU: variable name changed: CPU - scmpsw%nirdf, GPU - fluxr(i, jj, 33)
         ! GPU: variable name changed: CPU - scmpsw%visbm, GPU - fluxr(i, jj, 30)
         ! GPU: variable name changed: CPU - scmpsw%visdf, GPU - fluxr(i, jj, 31)
         ! GPU: the following variable name changed when using icmphys = 15 or 16 
         ! GPU: (GCE microphysics)
         ! GPU:                        CPU - clouds(:,:,1), GPU - cldfrc
         ! GPU:                        CPU - clouds(:,:,2), GPU - cwp
         ! GPU:                        CPU - clouds(:,:,3), GPU - rew
         ! GPU:                        CPU - clouds(:,:,4), GPU - cip
         ! GPU:                        CPU - clouds(:,:,5), GPU - rei
         ! GPU:                        CPU - clouds(:,:,6), GPU - crp
         ! GPU:                        CPU - clouds(:,:,7), GPU - rer
         ! GPU:                        CPU - clouds(:,:,8), GPU - csp
         ! GPU:                        CPU - clouds(:,:,9), GPU - res
         !allocate(tauae(nxptot, lm+ltp,nbdsw))
         !allocate(ssaae(nxptot, lm+ltp,nbdsw))
         !allocate(asyae(nxptot, lm+ltp,nbdsw))
         !!$acc enter data create(tauae, ssaae, asyae) async(async_id)
         call swrad_gpu                                                  &
            !  ---  inputs for setaer_sw_gpu:
                  ( prslk1, tvly, rhly, slmsk, tracer1, xlon, xlat, lsswr, lslwr, &
                    me, map_jj, map_i, my_max, ntrac, &
            !  ---  inputs:
            &       plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr_co2, gasvmr_other,                      &
            &       clouds, cldfrc, cwp, rew, cip, rei, crp, rer, csp, res,icsdsw,&
                    sfcalb,       &
            &       coszen,solcon, nday,idxday, map_nday_ipt, map_nday_jj, offset_nday,   &
            &       myim, lmk, lmp, lprnt, myrank, ix, kd,                              &
                    nf_clds, nf_vgas, nf_albd, nf_aesw, max_nday_length, async_id, &
                    my_max, swrad_block, swrad_smalljj, lhtrswb, &
            !  ---  outputs:
            &       htrsw, fluxr(:, :, 2 ), fluxr(:, :, 1 ), fluxr(:, :, 22), &
                    fluxr(:, :, 5 ), fluxr(:, :, 4 ), fluxr(:, :, 25), fluxr(:, :, 24), &
            !! ---  optional:
            &       hsw0=htrsw0,hswb=htrswb,                                     &
            &       flxprf_upfxc=fusl, flxprf_dnfxc=fdsl, flxprf_upfx0=fuslr, flxprf_dnfx0=fdslr, &
                    fdncmp_uvbfc=fluxr(:, :, 28), fdncmp_uvbf0=fluxr(:, :, 29), &
                    fdncmp_nirbm=fluxr(:, :, 32), fdncmp_nirdf=fluxr(:, :, 33), &
                    fdncmp_visbm=fluxr(:, :, 30), fdncmp_visdf=fluxr(:, :, 31)  &
            &     )
         !call nvtxEndRange
         if (lhtrswb) then
            !$acc exit data delete(htrswb) async(async_id)
            deallocate(htrswb)
         end if
      !!$acc exit data delete(tauae, ssaae, asyae) async(async_id)
      !deallocate(tauae)
      !deallocate(ssaae)
      !deallocate(asyae)

      !$acc exit data delete(qstl) async(async_id)
      !$acc exit data delete(prslk1, tracer1, rhly) async(async_id)
      
               

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
!         if (present(htrsw0)) then
!         endif
!         if (present(fswprf)) then
!         endif

         !$acc parallel loop gang collapse(2) private(jj, i) async(async_id)
         do k = 1, lm
            do n = 1, nxptot
               jj = map_jj(n)
               i = map_i(n)
               if (nday(jj) <= 0) then ! if_nday_block
                  htrsw(i,k, jj) = 0.0
               end if
            enddo
         end do
      
!! ---  optional:

!         if ( present(htrswb) ) then
!         endif

!         if ( present(htrsw0) ) then
         !$acc parallel loop collapse(2) private(jj, i) async(async_id)
         do k = 1, lm
            do n = 1, nxptot
               jj = map_jj(n)
               i = map_i(n)
               if (nday(jj) <= 0) then ! if_nday_block
                  htrsw0(i,k, jj) = 0.0
!         endif
                  
               endif                  ! end_if_nday
            enddo
         enddo
      endif                                ! end_if_lsswr

      do n = 1, nxptot
         jj = map_jj(n)
         i = map_i(n)
         if (i .eq. ipt) then
            if (me == 0 .and. myrank ==0) then
               print *,'###################################################' 
               print *,'### call swrad ok!!'
               print *,'###################################################' 
               print *,'### im=',myim(jj),' lmk=',lmk,' lmp=',lmp
               print *,'### nday=',nday(jj)
               print *,'### solcon=',solcon
               print *,'###################################################' 
               print *,'### icsdsw(ipt)=',icsdsw(ipt, jj)
               print *,'###################################################' 
               print *,'### plvl(ipt,lm)=',plvl(n,lm)
               print *,'### plyr(ipt,lm)=',plyr(n,lm)
               print *,'###################################################' 
               print *,'### tlyr(ipt,lm)=',tlyr(n,lm)
               print *,'### tlvl(ipt,lm)=',tlvl(n,lm)
               print *,'###################################################' 
               print *,'### qlyr(ipt,lm)=',qlyr(n,lm)
               print *,'### olyr(ipt,lm)=',olyr(n,lm)
               print *,'###################################################' 
               print *,'### gasvmr(ipt,lm,1)_co2 =',gasvmr_co2(n,lm)
               print *,'### gasvmr(ipt,lm,2)_n2o =',gasvmr_other(2)
               print *,'### gasvmr(ipt,lm,3)_ch4 =',gasvmr_other(3)
               print *,'### gasvmr(ipt,lm,4)_o2  =',gasvmr_other(4)
               print *,'### gasvmr(ipt,lm,5)_co  =',gasvmr_other(5)
               print *,'### gasvmr(ipt,lm,6)_cf11=',gasvmr_other(6)
               print *,'### gasvmr(ipt,lm,7)_cf12=',gasvmr_other(7)
               print *,'### gasvmr(ipt,lm,8)_cf22=',gasvmr_other(8)
               print *,'### gasvmr(ipt,lm,9)_ccl4=',gasvmr_other(9)
               print *,'###################################################' 
               print *,'### clouds(ipt,lm,1)-total cloud fraction=',clouds(n,lm,1)
               print *,'### clouds(ipt,lm,2)-liq water path=',clouds(n,lm,2)
               print *,'### clouds(ipt,lm,3)-liq effective radius=',clouds(n,lm,3)
               print *,'### clouds(ipt,lm,4)-ice water path=',clouds(n,lm,4)
               print *,'### clouds(ipt,lm,5)-ice effective radius=',clouds(n,lm,5)
               print *,'### clouds(ipt,lm,6)-rain water path=',clouds(n,lm,6)
               print *,'### clouds(ipt,lm,7)-rain effective radius=',clouds(n,lm,j)
               print *,'### clouds(ipt,lm,8)-snow water path=',clouds(n,lm,8)
               print *,'### clouds(ipt,lm,9)-snow effective radius=',clouds(n,lm,j)
               print *,'###################################################' 
               !print *,'### faersw(ipt,lm,1,1)=sw#1-opd',tauae(n,lm,1)
               !print *,'### faersw(ipt,lm,1,2)=sw#1-ssa',ssaae(n,lm,1)
               !print *,'### faersw(ipt,lm,1,3)=sw#1-asy',asyae(n,lm,1)
               print *,'###################################################' 
               print *,'### sfcalb(ipt,1)=',sfcalb(n,1)
               print *,'### sfcalb(ipt,2)=',sfcalb(n,2)
               print *,'### sfcalb(ipt,3)=',sfcalb(n,3)
               print *,'### sfcalb(ipt,4)=',sfcalb(n,4)
               print *,'###################################################' 
               print *,'### htswc(ipt,lm) =',htrsw(ipt,lm, jj)
               print *,'### htrsw(ipt,lm) =',htrsw(ipt,lm, jj)
               print *,'### htrsw0(ipt,lm)=',htrsw0(ipt,lm, jj)
               print *,'###################################################' 
               print *,'### sfcfsw(ipt)%upfxc=',fluxr(i, jj, 5)
               print *,'### sfcfsw(ipt)%dnfxc=',fluxr(i, jj, 4)
               print *,'### sfcfsw(ipt)%upfx0=',fluxr(i, jj, 25)
               print *,'### sfcfsw(ipt)%dnfx0=',fluxr(i, jj, 24)
               print *,'###################################################' 
               print *,'### topfsw(ipt)%upfxc=',fluxr(i, jj, 2 )
               print *,'### topfsw(ipt)%dnfxc=',fluxr(i, jj, 1 )
               print *,'### topfsw(ipt)%upfx0=',fluxr(i, jj, 22)
               print *,'###################################################' 
               print *,'### scmpsw(ipt)%uvbfc=',fluxr(i, jj, 28)
               print *,'### scmpsw(ipt)%uvbf0=',fluxr(i, jj, 29)
               print *,'### scmpsw(ipt)%nirbm=',fluxr(i, jj, 32)
               print *,'### scmpsw(ipt)%nirdf=',fluxr(i, jj, 33)
               print *,'### scmpsw(ipt)%visbm=',fluxr(i, jj, 30)
               print *,'### scmpsw(ipt)%visdf=',fluxr(i, jj, 31)
               print *,'###################################################' 
               print *,'### fswprf(ipt,1)%upfxc=',fusl(ipt,1+kd, jj)
               print *,'### fswprf(ipt,1)%dnfxc=',fdsl(ipt,1+kd, jj)
               print *,'### fswprf(ipt,1)%upfx0=',fuslr(ipt,1+kd, jj)
               print *,'### fswprf(ipt,1)%dnfx0=',fdslr(ipt,1+kd, jj)
               print *,'###################################################' 
               print *,'### fswprf(ipt,61)%upfxc=',fusl(ipt,61+kd, jj)
               print *,'### fswprf(ipt,61)%dnfxc=',fdsl(ipt,61+kd, jj)
               print *,'### fswprf(ipt,61)%upfx0=',fuslr(ipt,61+kd, jj)
               print *,'### fswprf(ipt,61)%dnfx0=',fdslr(ipt,61+kd, jj)
               print *,'###################################################' 
               print *,'### fusl(ipt,1)=',fusl(ipt,1, jj)
               print *,'### fdsl(ipt,1)=',fdsl(ipt,1, jj)
               print *,'### fuslr(ipt,1)=',fuslr(ipt,1, jj)
               print *,'### fdslr(ipt,1)=',fdslr(ipt,1, jj)
               print *,'###################################################' 
               print *,'### fusl(ipt,61)=',fusl(ipt,61, jj)
               print *,'### fdsl(ipt,61)=',fdsl(ipt,61, jj)
               print *,'### fuslr(ipt,61)=',fuslr(ipt,61, jj)
               print *,'### fdslr(ipt,61)=',fdslr(ipt,61, jj)
               print *,'###################################################' 
            endif
         end if
      end do

!----------------------------------------------------------------------
      if (lslwr) then

!  ---  setup surface emissivity for lw radiation
         !call nvtxStartRange("setemis")

         call setemis_gpu                                                    &
            !  ---  inputs:
            &     ( xlon,xlat,slmsk,&
            snowd,sncovr, &
            zorl,tsfg,tsfa,hprim,         &
            &       myim, ix, map_jj, map_i, nxptot, async_id,                                                       &
            !  ---  outputs:
            &       sfcemis                                                    &
            &     )
         !call nvtxEndRange
!
         do n = 1, nxptot
            jj = map_jj(n)
            i = map_i(n)
            if (i .eq. ipt) then
               if (me == 0 .and. myrank ==0) then
                  print *,'###################################################' 
                  print *,'### call setemis ok!!'
                  print *,'###################################################' 
                  print *,'### im=',myim(jj),' ipt=',ipt
                  print *,'###################################################' 
                  print *,'### xlon(ipt)=',xlon(ipt, jj)
                  print *,'### xlat(ipt)=',xlat(ipt, jj)
                  print *,'###################################################' 
                  print *,'### slmsk(ipt)=',slmsk(ipt, jj)
                  print *,'### snowd(ipt)=',snowd(ipt, jj)
                  print *,'### sncovr(ipt)=',sncovr(ipt, jj)
                  print *,'### zorl(ipt)=',zorl(ipt, jj)
                  print *,'###################################################' 
                  print *,'### tsfg(ipt)=',tsfg(n)
                  print *,'### tsfa(ipt)=',tsfa(n)
                  print *,'### hprim(ipt)=',hprim(ipt, jj)
                  print *,'###################################################' 
                  print *,'### sfcemis(ipt)=',sfcemis(n)
                  print *,'###################################################' 
               endif
            end if
         end do

         if (lhtrlwb) then
            allocate(htrlwb(ix,lm,nbdlw, my_max))
            !$acc enter data create(htrlwb) async(async_id)
         end if
         !$acc enter data create(tracer1, rhly, qstl, prslk1) async(async_id)
         !$acc parallel loop collapse(2) private(k1, jj, i, epsm2, qqq, qss) async(async_id)
         do k = 1, lm
            do n = 1, nxptot
               jj = map_jj(n)
               i = map_i(n)
               k1 = k + kd
               plyr(n,k1)   = 10.0 * prsl(i,k, jj)   ! cb (kpa) to mb (hpa)
               tlyr(n,k1)   = tgrs(i,k, jj)
               prslk1(n,k1) = prslk(i,k, jj)
               qlyr(n,k1) = max( qme6, qgrs(i,k, jj) )
               epsm2=0.622-1.
               qqq = min ( plyr(n,k1) , 0.01*fpvs_gpu(tlyr(n,k1),c1xpvs,c2xpvs,tbpvs) )
               qss = 0.622*qqq/(plyr(n,k1)+epsm2*qqq)
               rhly(n,k1)= max( 0.0, min( 1.0, max(qmin, qlyr(n,k1))/qss ) )
               qstl(n,k1) = qss
               !$acc loop seq
               do j = 1, ntrac
                  tracer1(n,k1,j) = tracer(i,k,j, jj)
               enddo
            enddo
         end do
         if ( lextop ) then                 ! values for extra top layer
            !$acc parallel loop private(jj, i) async(async_id)
            do n = 1, nxptot
               prslk1(n,lyb) = (plyr(n,lyb)*0.001) ** rocp ! plyr in hpa

               rhly(n,lyb)   = rhly(n,lya)
               qstl(n,lyb)   = qstl(n,lya)
               !$acc loop seq
               do j = 1, ntrac
                  tracer1(n,lyb,j) = tracer1(n,lya,j)
               end do
            end do
         end if


!       if ( present(htrlw0) .and. present(flwprf) ) then
         if (me == 0 .and. myrank ==0) print *,'#### call lwrad ####'
         !call nvtxStartRange("lwrad")
         ! GPU: variable name changed: CPU - htlwb, GPU - htrlwb
         ! GPU: variable name changed: CPU - htlw0, GPU - htrlw0
         ! GPU: variable name changed: CPU - htlwc, GPU - htrlw
         ! GPU: variable name changed: CPU - topflw%upfxc, GPU - fluxr(i, jj, 3)
         ! GPU: variable name changed: CPU - topflw%upfx0, GPU - fluxr(i, jj, 23)
         ! GPU: variable name changed: CPU - sfcflw%upfxc, GPU - fluxr(i, jj, 7)
         ! GPU: variable name changed: CPU - sfcflw%upfx0, GPU - fluxr(i, jj, 27)
         ! GPU: variable name changed: CPU - sfcflw%dnfxc, GPU - fluxr(i, jj, 6)
         ! GPU: variable name changed: CPU - sfcflw%dnfx0, GPU - fluxr(i, jj, 26)
         ! GPU: variable name changed: CPU - flwprf%upfxc, GPU - fuir
         ! GPU: variable name changed: CPU - flwprf%dnfxc, GPU - fdir
         ! GPU: variable name changed: CPU - flwprf%upfx0, GPU - fuirr
         ! GPU: variable name changed: CPU - flwprf%dnfx0, GPU - fdirr
         call lwrad_gpu                                                    &
            !  ---  inputs for setaer_lw_gpu:
            &     ( prslk1, tvly, rhly, slmsk, tracer1, xlon, xlat, lsswr, lslwr, me, ntrac, my_max,        &
            !  ---  inputs:
            &       plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr_co2, gasvmr_other,                      &
            &       clouds, cldfrc, cwp, rew, cip, rei, crp, rer, csp, res,icsdlw, sfcemis,tsfg,  &
            &       myim, map_jj, map_i, lmk, lmp, lprnt,myrank, ix, kd,                                &
                     nf_vgas, nf_clds, nf_aelw, max_nxjp_acc_length, async_id, &
                     my_max, lwrad_block, lwrad_smalljj, lhtrlwb, &
            !  ---  outputs:
            &       htrlw,fluxr(:, :, 3), fluxr(:, :, 23), fluxr(:, :, 7), &
                    fluxr(:, :, 27), fluxr(:, :, 6), fluxr(:, :, 26)                                        &
            !! ---  optional:
            &,      hlw0=htrlw0,hlwb=htrlwb, flxprf_upfxc=fuir, flxprf_dnfxc=fdir, &
                    flxprf_upfx0=fuirr, flxprf_dnfx0=fdirr                        &
            &     )
         !call nvtxEndRange
         if (lhtrlwb) then
            !$acc exit data delete(htrlwb) async(async_id)
            deallocate(htrlwb)
         end if
         !$acc exit data delete(prslk1, tracer1, rhly, tvly, qstl) async(async_id)


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

         !$acc parallel loop private(jj, i) async(async_id)
         do n = 1, nxptot
            jj = map_jj(n)
            i = map_i(n)
            semis (i, jj) = sfcemis(n)
            !  ---  save surface air temp for diurnal adjustment at model t-steps
            tsflw (i, jj) = tsfa(n)
         end do

!       if (present(htrlw0)) then
!       endif
!         if (present(flwprf)) then
!          endif

      endif                                ! end_if_lslwr

!----------------------------------------------------------------------
      do n = 1, nxptot
         jj = map_jj(n)
         i = map_i(n)
         if (i .eq. ipt) then
            if (me == 0 .and. myrank ==0) then
               print *,'###################################################' 
               print *,'### call lwrad ok!!'
               print *,'###################################################' 
               print *,'### im=',myim(jj),' lmk=',lmk,' lmp=',lmp,' kd=',kd
               print *,'###################################################' 
               print *,'### icsdlw(ipt)=',icsdlw(ipt, jj)
               print *,'###################################################' 
               print *,'### plvl(ipt,lm)=',plvl(n,lm)
               print *,'### plyr(ipt,lm)=',plyr(n,lm)
               print *,'###################################################' 
               print *,'### tlyr(ipt,lm)=',tlyr(n,lm)
               print *,'### tlvl(ipt,lm)=',tlvl(n,lm)
               print *,'###################################################' 
               print *,'### qlyr(ipt,lm)=',qlyr(n,lm)
               print *,'### olyr(ipt,lm)=',olyr(n,lm)
               print *,'###################################################' 
               print *,'### gasvmr(ipt,lm,1)_co2 =',gasvmr_co2(n,lm)
               print *,'### gasvmr(ipt,lm,2)_n2o =',gasvmr_other(2)
               print *,'### gasvmr(ipt,lm,3)_ch4 =',gasvmr_other(3)
               print *,'### gasvmr(ipt,lm,4)_o2  =',gasvmr_other(4)
               print *,'### gasvmr(ipt,lm,5)_co  =',gasvmr_other(5)
               print *,'### gasvmr(ipt,lm,6)_cf11=',gasvmr_other(6)
               print *,'### gasvmr(ipt,lm,7)_cf12=',gasvmr_other(7)
               print *,'### gasvmr(ipt,lm,8)_cf22=',gasvmr_other(8)
               print *,'### gasvmr(ipt,lm,9)_ccl4=',gasvmr_other(9)
               print *,'###################################################' 
               print *,'### clouds(ipt,lm,1)-total cloud fraction=',clouds(n,lm,1)
               print *,'### clouds(ipt,lm,2)-liq water path=',clouds(n,lm,2)
               print *,'### clouds(ipt,lm,3)-liq effective radius=',clouds(n,lm,3)
               print *,'### clouds(ipt,lm,4)-ice water path=',clouds(n,lm,4)
               print *,'### clouds(ipt,lm,5)-ice effective radius=',clouds(n,lm,5)
               print *,'### clouds(ipt,lm,6)-rain water path=',clouds(n,lm,6)
               print *,'### clouds(ipt,lm,7)-rain effective radius=',clouds(n,lm,7)
               print *,'### clouds(ipt,lm,8)-snow water path=',clouds(n,lm,8)
               print *,'### clouds(ipt,lm,9)-snow effective radius=',clouds(n,lm,9)
               print *,'###################################################' 
               !print *,'### faerlw(ipt,lm,1,1)=lw#1-opd',faerlw(ipt,lm,1,1, jj)
               !print *,'### faerlw(ipt,lm,1,2)=lw#1-ssa',faerlw(ipt,lm,1,2, jj)
               !print *,'### faerlw(ipt,lm,1,3)=lw#1-asy',faerlw(ipt,lm,1,3, jj)
               print *,'###################################################' 
               print *,'### sfcemis(ipt)=',sfcemis(n)
               print *,'###################################################' 
               print *,'### tsfg(ipt) =',tsfg(n)
               print *,'### tsfa(ipt) =',tsfa(n)
               print *,'### tsflw(ipt)=',tsflw(ipt, jj)
               print *,'###################################################' 
               print *,'### htswc(ipt,lm) =',htrsw(ipt,lm, jj)
               print *,'### htsw0(ipt,lm) =',htrsw0(ipt,lm, jj)
               print *,'### htrsw(ipt,lm) =',htrsw(ipt,lm, jj)
               print *,'### htrsw0(ipt,lm)=',htrsw0(ipt,lm, jj)
               print *,'###################################################' 
               print *,'### htlwc(ipt,lm) =',htrlw(ipt,lm+kd, jj)
               print *,'### htlw0(ipt,lm) =',htrlw0(ipt,lm+kd, jj)
               print *,'### htrlw(ipt,lm) =',htrlw(ipt,lm, jj)
               print *,'### htrlw0(ipt,lm)=',htrlw0(ipt,lm, jj)
               print *,'###################################################' 
               print *,'### topflw(ipt)%upfxc=',fluxr(i, jj, 3)
               print *,'### topflw(ipt)%upfx0=',fluxr(i, jj, 23)
               print *,'###################################################' 
               print *,'### topfsw(ipt)%upfxc=',fluxr(i, jj, 2 )
               print *,'### topfsw(ipt)%upfx0=',fluxr(i, jj, 1 )
               print *,'### topfsw(ipt)%dnfxc=',fluxr(i, jj, 22)
               print *,'###################################################' 
               print *,'### sfcflw(ipt)%upfxc=',fluxr(i, jj, 7)
               print *,'### sfcflw(ipt)%upfx0=',fluxr(i, jj, 27)
               print *,'### sfcflw(ipt)%dnfxc=',fluxr(i, jj, 6)
               print *,'### sfcflw(ipt)%dnfx0=',fluxr(i, jj, 26)
               print *,'###################################################' 
               print *,'### sfcfsw(ipt)%upfxc=',fluxr(i, jj, 5)
               print *,'### sfcfsw(ipt)%upfx0=',fluxr(i, jj, 4)
               print *,'### sfcfsw(ipt)%dnfxc=',fluxr(i, jj, 25)
               print *,'### sfcfsw(ipt)%dnfx0=',fluxr(i, jj, 24)
               print *,'###################################################' 
               print *,'### flwprf(ipt,1)%upfxc=',fuir(ipt,1+kd, jj)
               print *,'### flwprf(ipt,1)%dnfxc=',fdir(ipt,1+kd, jj)
               print *,'### flwprf(ipt,1)%upfx0=',fuirr(ipt,1+kd, jj)
               print *,'### flwprf(ipt,1)%dnfx0=',fdirr(ipt,1+kd, jj)
               print *,'###################################################' 
               print *,'### flwprf(ipt,61)%upfxc=',fuir(ipt,61+kd, jj)
               print *,'### flwprf(ipt,61)%dnfxc=',fdir(ipt,61+kd, jj)
               print *,'### flwprf(ipt,61)%upfx0=',fuirr(ipt,61+kd, jj)
               print *,'### flwprf(ipt,61)%dnfx0=',fdirr(ipt,61+kd, jj)
               print *,'###################################################' 
               print *,'### fuir(ipt,1)=',fuir(ipt,1, jj)
               print *,'### fdir(ipt,1)=',fdir(ipt,1, jj)
               print *,'### fuirr(ipt,1)=',fuirr(ipt,1, jj)
               print *,'### fdirr(ipt,1)=',fdirr(ipt,1, jj)
               print *,'###################################################' 
               print *,'### fuir(ipt,61)=',fuir(ipt,61, jj)
               print *,'### fdir(ipt,61)=',fdir(ipt,61, jj)
               print *,'### fuirr(ipt,61)=',fuirr(ipt,61, jj)
               print *,'### fdirr(ipt,61)=',fdirr(ipt,61, jj)
               print *,'###################################################' 
            endif
         end if
      end do
!
!


!  ---  save cld frac,toplyr,botlyr and top temp, note that the order
!       of h,m,l cloud is reversed for the fluxr output.
!  ---  save interface pressure (cb) of top/bot

      !$acc parallel loop collapse(2) private(tem0d, itop, ibtc, jj, i) async(async_id)
      do j = 1, 3
         do n = 1, nxptot
            jj = map_jj(n)
            i = map_i(n)
            tem0d = cldsa(n,j)
            itop  = mtopa(n,j) - kd
            ibtc  = mbota(n,j) - kd
            fluxr(i, jj, 11-j) = tem0d  ! cloud fraction(h,m,l,j=3,2,1)
            fluxr(i, jj, 14-j) = prsi(i,itop+kt, jj) ! cloud top pressure
            fluxr(i, jj, 17-j) = prsi(i,ibtc+kb, jj) ! cloud bot pressure
            fluxr(i, jj, 20-j) = tgrs(i,itop, jj)    ! cloud top temp
         enddo
      end do

!  ---  save total cloud and bl cloud fraction
      !$acc parallel loop private(jj, i) async(async_id)
      do n = 1, nxptot
         jj = map_jj(n)
         i = map_i(n)
         fluxr(i, jj, 20) = cldsa(n,4) ! total cloud fraction
         fluxr(i, jj, 21) = cldsa(n,5) ! BL domain cloud fraction
      end do

!
!  ---  sw uv-b fluxes

      !!$acc parallel loop gang collapse(2) private(k1) async(async_id)
      !do jj = 1, jlistnum
      !   do k = 1, lm
      !      k1 = k + kd
      !         !$acc loop vector
      !      do i = 1, myim(jj)
      !         cldcov(i,k, jj) = clouds(i,k1,1, jj)
      !      enddo
      !   enddo
      !end do ! jj-loop


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
      !$acc end data


!
      return
!...................................
   end subroutine grrad_gpu
!-----------------------------------


!
!........................................!
end module module_radiation_driver_gpu !
!========================================!
