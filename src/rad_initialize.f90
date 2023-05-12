!-----------------------------------
      subroutine rad_initialize                                         &
!...................................
!  ---  inputs:
     &     ( si,levr,ictm,isol,ico2,iaer,ialb,iems,ntcw,                &
     &       nmmiph,ntoz,iovr_sw,iovr_lw,isubc_sw,isubc_lw,             &
     &       icliq_sw,icice_sw,icliq_lw,icice_lw,                       &
     &       sashal,crick_proof,ccnorm,norad_precip,idate,iflip,me,myrank)
!  ---  outputs: ( none )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:   rad_initialize - a subprogram to initialize radiation   !
!                                                                       !
! usage:        call rad_initialize                                     !
!                                                                       !
! attributes:                                                           !
!   language:  fortran 90                                               !
!                                                                       !
! program history:                                                      !
!   mar 2012  - yu-tai hou   create the program to initialize fixed     !
!                 control variables for radiaion processes.  this       !
!                 subroutine is called at the start of model run.       !
!   nov 2012  - yu-tai hou   modified control parameter through         !
!                 module 'physparam'.                                   !
!   mar 2014  - sarah lu  iaermdl is determined from iaer               !        
!   jun 2014  - h-m lin/y-t hou   added run-time adjustable cloud       !
!                 optical property control parameters.                  !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
! input parameters:                                                     !
!   si               : model vertical sigma interface or equivalence    !
!   levr             : number of model vertical layers                  !
!   ictm             :=yyyy#, external data time/date control flag      !
!                     =   -2: same as 0, but superimpose seasonal cycle !
!                             from climatology data set.                !
!                     =   -1: use user provided external data for the   !
!                             forecast time, no extrapolation.          !
!                     =    0: use data at initial cond time, if not     !
!                             available, use latest, no extrapolation.  !
!                     =    1: use data at the forecast time, if not     !
!                             available, use latest and extrapolation.  !
!                     =yyyy0: use yyyy data for the forecast time,      !
!                             no further data extrapolation.            !
!                     =yyyy1: use yyyy data for the fcst. if needed, do !
!                             extrapolation to match the fcst time.     !
!   isol             := 0: use the old fixed solar constant in "physcon"!
!                     =10: use the new fixed solar constant in "physcon"!
!                     = 1: use noaa ann-mean tsi tbl abs-scale data tabl!
!                     = 2: use noaa ann-mean tsi tbl tim-scale data tabl!
!                     = 3: use cmip5 ann-mean tsi tbl tim-scale data tbl!
!                     = 4: use cmip5 mon-mean tsi tbl tim-scale data tbl!
!   ico2             :=0: use prescribed global mean co2 (old  oper)    !
!                     =1: use observed co2 annual mean value only       !
!                     =2: use obs co2 monthly data with 2-d variation   !
!   iaer             : 4-digit aerosol flag (dabc for aermdl,volc,lw,sw)!
!                     d: =0 or none, opac-climatology aerosol scheme    !                
!                        =1 use gocart climatology aerosol scheme       !  
!                        =2 use gocart progostic aerosol scheme         !  
!                     a: =0 use background stratospheric aerosol        !
!                        =1 incl stratospheric vocanic aeros            !
!                     b: =0 no topospheric aerosol in lw radiation      !
!                        =1 include tropspheric aerosols for lw         !
!                     c: =0 no topospheric aerosol in sw radiation      !
!                        =1 include tropspheric aerosols for sw         !
!   ialb             : control flag for surface albedo schemes          !
!                     =0: climatology, based on surface veg types       !
!                     =1: modis retrieval based surface albedo scheme   !
!   iems             : ab 2-digit control flag                          !
!                     a: =0 set sfc air/ground t same for lw radiation  !
!                        =1 set sfc air/ground t diff for lw radiation  !
!                     b: =0 use fixed sfc emissivity=1.0 (black-body)   !
!                        =1 use varying climtology sfc emiss (veg based)!
!                        =2 future development (not yet)                !
!   ntcw             :=0 no cloud condensate calculated                 !
!                     >0 array index location for cloud condensate      !
!   nmmiph          :=3: ferrier's microphysics cloud scheme           !
!                     =4: zhao/carr/sundqvist microphysics cloud        !
!                     =5: WSM6 & Thompson microphysics cloud            !
!   ntoz             : ozone data control flag                          !
!                     =0: use climatological ozone profile              !
!                     >0: use interactive ozone profile                 !
!   iovr_sw/iovr_lw  : control flag for cloud overlap (sw/lw rad)       !
!                     =0: random overlapping clouds                     !
!                     =1: max/ran overlapping clouds                    !
!   isubc_sw/isubc_lw: sub-column cloud approx control flag (sw/lw rad) !
!                     =0: with out sub-column cloud approximation       !
!                     =1: mcica sub-col approx. prescribed random seed  !
!                     =2: mcica sub-col approx. provided random seed    !
!   icliq_sw         : sw liq-cloud opt prop contol flag (sw rad)       !
!                     =0: diag cld opt depth input, ignore icice_sw.    !
!                     =1: hu & stamnes method for sw liq-cld opt prop.  !
!   icice_sw         : sw ice-cloud opt prop contol flag (sw rad)       !
!                     =1: ebert & curry method for sw ice-cld opt prop. !
!                     =2: streamer v3.0 method for sw ice-cld opt prop. !
!                     =3: fu method for sw ice-cld opt prop.            !
!   icliq_lw         : lw liq-cloud opt prop contol flag (lw rad)       !
!                     =0: diag cld opt depth input, ignore icice_lw.    !
!                     =1: hu & stamnes method for lw liq-cld opt prop.  !
!   icice_lw         : lw ice-cloud opt prop contol flag (lw rad)       !
!                     =1: ebert & curry method for lw ice-cld opt prop. !
!                     =2: streamer method for lw ice-cld opt prop.      !
!                     =3: fu method for lw ice-cld opt prop.            !
!   sashal           : shallow convection scheme flag                   !
!   crick_proof      : control flag for eliminating CRICK               !
!   ccnorm           : control flag for in-cloud condensate mixing ratio!
!   norad_precip     : control flag for not using precip in radiation   !
!   idate(8)         : ncep absolute date and time of initial condition !
!                      (hour, month, day, year)                         !
!   iflip            : control flag for direction of vertical index     !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   me               : print control flag                               !
!                                                                       !
!  subroutines called: radinit                                          !
!                                                                       !
!  ===================================================================  !
!
      use physpara , only : isolar , ictmflg, ico2flg, ioznflg, iaerflg,&
     &             iaermdl, laswflg, lalwflg, lavoflg, icldflg, icmphys,&
     &             iovrsw , iovrlw , lsashal, lcrick , lcnorm , lnoprec,&
     &             ialbflg, iemsflg, isubcsw, isubclw, ivflip , ipsd0,  &
     &             kind_phys, iswcliq, iswcice, ilwcliq, ilwcice
      use const,     only : pdfcloud

      use module_radiation_driver, only : radinit
! 

      implicit   none
      

!  ---  input:
      integer,  intent(in) :: levr, ictm, isol, ico2, iaer,             &
     &       ntcw, ialb, iems, nmmiph, ntoz, iovr_sw, iovr_lw,          &
     &       icliq_sw,icice_sw,icliq_lw,icice_lw,                       &
     &       isubc_sw, isubc_lw, iflip, me, idate(8),myrank

      real (kind=kind_phys), intent(in) :: si(levr+1)

      logical, intent(in) :: sashal, crick_proof, ccnorm, norad_precip

!  ---  output: ( none )

!  ---  local:
      integer :: icld
!
!===> ...  start here
!
!  ---  set up parameters for radiation initialization

      isolar = isol                     ! solar constant control flag

      ictmflg= ictm                     ! data ic time/date control flag
      ico2flg= ico2                     ! co2 data source control flag
      ioznflg= ntoz                     ! ozone data source control flag

      if ( ictm==0 .or. ictm==-2 ) then
        iaerflg = mod(iaer, 100)        ! no volcanic aerosols for clim hindcast
      else
        iaerflg = mod(iaer, 1000)   
      endif
      laswflg= (mod(iaerflg,10) > 0)    ! control flag for sw tropospheric aerosol
      lalwflg= (mod(iaerflg/10,10) > 0) ! control flag for lw tropospheric aerosol
      lavoflg= (iaerflg >= 100)         ! control flag for stratospheric volcanic aeros
      iaermdl = iaer/1000               ! control flag for aerosol scheme selection                              
      if ( iaermdl < 0 .or.  iaermdl > 2) then
      if (myrank.eq.0) print *, ' Error -- IAER flag is incorrect, Abort'
         stop 7777
      endif

      if ( ntcw > 0 ) then
        icldflg = 1                     ! prognostic cloud optical prop scheme

        if ( nmmiph == 2 ) then
          if ( pdfcloud ) then
            icmphys = 3                 ! zhao/moorthi's prognostic with PDF cloud scheme 
          else
            icmphys = 1                 ! zhao/moorthi's prognostic cloud scheme
          endif
        elseif ( nmmiph == 3 ) then
          icmphys = 2                   ! ferrier's microphysics
        elseif ( nmmiph == 6 ) then
          icmphys = 6                   ! WSM6 microphysics
        elseif ( nmmiph == 8 ) then
          icmphys = 8                   ! Thompson microphysics
        elseif ( nmmiph == 18 ) then
          icmphys = 18                  ! 2M Thompson microphysics
        elseif ( nmmiph == 11 ) then
          icmphys = 11                  ! GFDL microphysics version 1
        elseif ( nmmiph == 12 ) then
          icmphys = 12                  ! GFDL microphysics version 2
        elseif ( nmmiph == 13 ) then
          icmphys = 13                  ! GFDL microphysics version 3
        elseif ( nmmiph == 15 ) then
          icmphys = 15                  ! Goddard (GCE) 3ICE microphysics
        elseif ( nmmiph == 16 ) then
          icmphys = 16                  ! Goddard (GCE) 4ICE microphysics
        endif
      else
        icldflg = 0                     ! diagnostic cloud optical prop scheme
      endif
      iovrsw = iovr_sw                  ! cloud overlapping control flag for sw
      iovrlw = iovr_lw                  ! cloud overlapping control flag for lw

      lsashal = sashal                  ! shallow convection scheme flag
      lcrick  = crick_proof             ! control flag for eliminating CRICK 
      lcnorm  = ccnorm                  ! control flag for in-cld condensate 
      lnoprec = norad_precip            ! precip effect on radiation flag (ferrier microphysics)
      isubcsw = isubc_sw                ! sub-column cloud approx flag in sw radiation
      isubclw = isubc_lw                ! sub-column cloud approx flag in lw radiation
      iswcliq = icliq_sw                ! sw liq-cld opt prop scheme flag in sw radiation
      iswcice = icice_sw                ! sw ice-cld opt prop scheme flag in sw radiation
      ilwcliq = icliq_lw                ! lw liq-cld opt prop scheme flag in lw radiation
      ilwcice = icice_lw                ! lw ice-cld opt prop scheme flag in lw radiation

      ialbflg= ialb                     ! surface albedo control flag
      iemsflg= iems                     ! surface emissivity control flag

      ivflip = iflip                    ! vertical index direction control flag

!  ---  assign initial permutation seed for mcica cloud-radiation
      if ( isubc_sw>0 .or. isubc_lw>0 ) then
!       ipsd0 = 17*idate(1)+43*idate(2)+37*idate(3)+23*idate(4) + ipsd0
        ipsd0 = 17*idate(5)+43*idate(2)+37*idate(3)+23*idate(1)
      endif

      if ( me == 0 .and. myrank == 0 ) then
        print *,'  In rad_initialize, before calling radinit'
        print *,' si =',si
        print *,' levr=',levr,' ictm=',ictm,' isol=',isol,' ico2=',ico2,&
                ' iaer=',iaer,' ialb=',ialb,' iems=',iems,' ntcw=',ntcw  
        print *,' nmmiph=',nmmiph,' ntoz=',ntoz,' iovr_sw=',iovr_sw,    &
                ' iovr_lw=',iovr_lw,' isubc_sw=',isubc_sw,              &
                ' isubc_lw=',isubc_lw,' iflip=',iflip,'  me=',me         
        print *,' sashal=',sashal,' crick_proof=',crick_proof,          &
                ' ccnorm=',ccnorm,' norad_precip=',norad_precip
      endif

      call radinit                                                      &
!  ---  inputs:
     &     (  si, levr, me, myrank )
!  ---  outputs:
!          ( none )

      if ( me == 0  .and. myrank  == 0 ) then
        print *,'  Radiation sub-cloud initial seed =',ipsd0,           &
     &          ' IC-idate =',idate
        print *,' return from rad_initialize - after calling radinit'
      endif
!
      return
!...................................
      end subroutine rad_initialize
!-----------------------------------
