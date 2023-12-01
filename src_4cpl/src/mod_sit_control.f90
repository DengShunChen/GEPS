      MODULE mod_sit_control

       IMPLICIT NONE
       	   


  !--------------------------------------------------------------------------
  !   USE mo_time_control,   ONLY: delta_time, lstart, get_time_step, current_date,   &
  !                             write_date, lobs_ocn_rerun, ltrigsit  
  !-----------------------
!ps     INTEGER :: nn             !   max meridional wave number for m=0.
      LOGICAL, SAVE :: lrere = .FALSE. !   .true. for IC (Initial Condition) run (forcasting mode),
      LOGICAL :: lhd      = .FALSE.  !   .true. for hydrologic discharge model
      LOGICAL :: ltrigsit = .TRUE.     ! need check GFS

!--------------------------------------------------------------------------
! from mo_kinds
!----------------------

  ! Number model from which the SELECTED_*_KIND are requested:
  !
  !                   4 byte REAL      8 byte REAL
  !          CRAY:        -            precision =   13
  !                                    exponent  = 2465
  !          IEEE:    precision =  6   precision =   15  
  !                   exponent  = 37   exponent  =  307 
  !
  ! Most likely this are the only possible models.


  ! Floating point section 

  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6,37)  
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)
  !INTEGER, PARAMETER :: dp = 4                            ! single precision
  INTEGER, PARAMETER :: wp = dp   ! working precision

  ! Integer section

  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)

  ! Missing Value

  real, PARAMETER:: xmissing= -9.e+33     ! default missing value
  REAL,     PARAMETER:: real_missing= -9.e+33    ! default missing value for REAL
  INTEGER,  PARAMETER:: int_missing= -99999     ! default missing value for INTEGER  

!--------------------------------------------------------------------------
! from mo_time_event
!-----------------------
  INTEGER, PARAMETER :: STR_LEN_A = 20
  !+
  ! **************** parameters ------------------------------------------------
  !
  ! TIME_INC_*     predefined counter units
  !
  CHARACTER(len=*), PUBLIC, PARAMETER :: &
       TIME_INC_SECONDS = 'seconds'  ,&!
       TIME_INC_MINUTES = 'minutes'  ,&!
       TIME_INC_HOURS   = 'hours'    ,&!
       TIME_INC_DAYS    = 'days'     ,&!
       TIME_INC_MONTHS  = 'months'   ,&!
       TIME_INC_YEARS   = 'years'      !

  ! TRIG_*  type of trigger adjustment of an event
  !
  CHARACTER(len=*), PUBLIC, PARAMETER :: &
       TRIG_FIRST  = 'first'  ,&! trigger in first step of counter unit
       TRIG_LAST   = 'last'   ,&! trigger in last step of counter unit
       TRIG_EXACT  = 'exact'  ,&! trigger without adjustment in side the unit
       TRIG_NONE   = 'off'     ! event trigger non active

  ! **************** structures ------------------------------------------------
  !
  TYPE, PUBLIC :: io_time_event            ! external given event properties
    INTEGER                  :: counter    = 0                ! No. of steps in given unit
    CHARACTER(len=STR_LEN_A) :: unit       = TIME_INC_SECONDS ! counter unit type
    CHARACTER(len=STR_LEN_A) :: adjustment = TRIG_EXACT       ! adjustment in side the unit
    INTEGER                  :: offset     = 0            ! offset to initial date in seconds
  END TYPE io_time_event


!--------------------------------------------------------------------------
! from mo_time_control
!-----------------------

  ! the adjustment of events for RERUN is dependent on the trigger step
  ! the trigger step can be the present date or the next date
  !
  CHARACTER(len=*), PARAMETER, PRIVATE :: &
       EV_TLEV_PRES    = 'present'        ,&! check event with present date
       EV_TLEV_NEXT    = 'next'           ,&! check event with next date
       TIME_INC_STEPS  = 'steps'          ,&! special event interval unit
       TIME_INC_ALWAYS = 'always'           ! special event always used  


  TYPE(io_time_event),SAVE ::    trigsit       &! time interval for triggerring sit_ocean
       = &                                      ! model
       io_time_event(1,TIME_INC_STEPS,TRIG_EXACT,0) ! every one time step
  TYPE(io_time_event),SAVE ::    trigocn       &! time interval for triggerring embedded ocean
       = &                                      ! model (every 1 time step)
       io_time_event(1,TIME_INC_STEPS,TRIG_EXACT,0)
 

  !---------------------------------------------------------------------
  !   USE mo_control,        ONLY: nn,lsit,lssst,lsit_ice,lsit_salt,lhd,                      &
  !                                sit_ice_option,locaf,lgodas,lrere,locn,                    &
  !                                lsice_nudg,lsit_lw,                                        &
  !                                ssit_restore_time,usit_restore_time,dsit_restore_time,     &
  !                                lwarning_msg,ocn_couple_option,                            &
  !                                socn_restore_time,uocn_restore_time,docn_restore_time,     &
  !                                obox_restore_time,obox_nudg_flag,Prw,csiced
  !   

  ! 3.0 SIT variables
  !   3.1 I/O
      INTEGER :: nocn_sv          = 51              ! *nocn_sv*         logical unit for ocn_sv file  
      INTEGER :: nwoa0            = 97   !  *nwoa0*     logical unit for world ocean atlas profile file (bjt), initial ocean field
      INTEGER :: ngodas           = 98   !  *ngodas*    logical unit for GODAS dataset (bjt), time series ocean field
      INTEGER :: nocaf            = 99   !  *nocaf*     logical unit for read_ocaf (bjt), flux correction terms
      INTEGER :: nrere            = 93  !  *nrere*      logical unit for ReReAnalysis run (xl and xi) (bjt)
      !   3.2 logical variables
      LOGICAL, SAVE :: lasia      = .FALSE. ! .true. for using etopo, new land surface data over rice paddy and Tibet
      !!! LOGICAL, SAVE :: lsit       = .FALSE. ! .true. for calculation of upper ocean temperature profile using sit model
      LOGICAL, SAVE :: lsice_nudg = .FALSE. ! .true. for nudging siced and seaice in SIT  
      LOGICAL, SAVE :: lsit_lw    = .FALSE.  ! .true. for turning LW code for water in SIT  
      LOGICAL, SAVE :: lsit_ice   = .TRUE.  ! .true. (default) for turn on the ice module of sit model
      LOGICAL, SAVE :: lsit_salt  = .TRUE.  ! .true. (default) for turn on salinity module of sit model
      INTEGER, SAVE :: zocn_option = 99   !  =1, 1 m appart in top 10 m, 10 m apart within 100-225 m, but coarse (~200 m) in deep-water formation depths (>500 m).
                                          !  =2 conventioanl godas z cord, fine in thermocline, 10 m apart within 100-225 m, but coarse (~200 m) in deep-water formation depths (>500 m).
                                          !  =others, conventional diecast coord.
      real, SAVE :: ocn_tlz=200.    ! depth of bottom z-level (m) of ocean model (5800. after v9.8992, 2013/10/3) (=5000. prior to v9.8992)
      INTEGER, SAVE :: ocn_k1=20            ! vertical dimension parameter (40 after v9.8992, in number of layer interfaces, equal
                                        ! the number of layers or pressure levels, does not include ghost zone)
      LOGICAL, SAVE :: lssst      = .TRUE.  ! .true. for turnon thermocline skin layer, .false. for turnoff thermocline skin layer
      LOGICAL, SAVE :: lgodas     = .TRUE. ! .true. for reading world ocean atlas (woa) data for sit model
      LOGICAL, SAVE :: lamip      = .TRUE. ! .false. for reading climatology world ocean atlas (woa) data for sit model
      LOGICAL, SAVE :: lwoa0      = .TRUE. ! .true. for reading initial ocean profile (unit: 97)
      INTEGER, SAVE :: lwarning_msg = 0     !  or printing warsning message
                                            !   =0, no message
                                            !   =1, basic message
                                            !   =2, medium message
                                            !   =3, many message
      LOGICAL, SAVE :: locaf = .FALSE.             ! .true. for q flux adjustment
      LOGICAL, SAVE :: locaf0 = .FALSE.             ! .true. for q flux adjustment
      real:: ocaf0_add = 0.             ! .true. for q flux adjustment
      LOGICAL, SAVE :: lcool_skin = .FALSE.        ! .true. for cool-skin parameterization
      LOGICAL, SAVE :: lsteady_TKE=.FALSE.          ! .true. = using calc_steady_TKE for TKE 
      LOGICAL, SAVE :: lwave_breaking=.TRUE.      !
  
  !   3.3 other variables  
      INTEGER, SAVE :: sit_ice_option  = 0  !   ice option in sit (i.e., calc. of snow/ice) (=0, off; >=1, on) 
                                            !   0: for coupling with vdiff semi-implicitly for tsi calculation (default)
                                            !   1: explcity coupling with strong security number
                                            !   2: original SIT output, no security. It can be crashed in few time steps
      INTEGER, SAVE :: maskid  = 1          !   0: DIECAST grids only
                                            !   1: Ocean and lakes (default)
                                            !   2: Ocean
                                            !   3: Ocean within 30N-30S
                                            !   4: all the grid                                       
      real:: ssit_restore_time =-99.    ! surface (0 <= ~ <10 m) sit grids restore time scale (s) (default: no nudging)
      real:: usit_restore_time =604800.      ! upper ocean (10 <= ~ <100 m) sit grids restore time scale (s) (default: no nudging)
      real:: dsit_restore_time =0.      ! deep (>=100 m) restore time scale (s)  (default: no nudging)

      real:: ssits_restore_time =-99.    ! surface (0 <= ~ <10 m) sit grids restore time scale for salinity (s) (default: no nudging)
      real:: usits_restore_time =604800.      ! upper ocean (10 <= ~ <100 m) sit grids restore time scale for salinity (s) (default: no nudging)
      real:: dsits_restore_time =0.      ! deep (>=100 m) restore time scale (s) for salinity (default: no nudging)

      real:: ssituv_restore_time =-99.    ! surface (0 <= ~ <10 m) sit grids restore time scale for u, v (s) (default: no nudging)
      real:: usituv_restore_time =604800.      ! upper ocean (10 <= ~ <100 m) sit grids restore time scale for u, v(s) (default: no nudging)
      real:: dsituv_restore_time =0.      ! deep (>=100 m) restore time scale (s) for u, v (default: no nudging)

      INTEGER :: nsit_nudg= 0                  ! number of nudg squares in ocean grids (default = 0, maximun=6)
      real:: sitbox_nudg_w(6)= -999.              ! west coords (lon) of nudging boxes [-180., 360.](deg). There are 6 boxes.
      real:: sitbox_nudg_e(6)= -999.              ! east coords (lon) of nudging boxes [-180., 360.](deg). There are 6 boxes.
      real:: sitbox_nudg_s(6)= -999.              ! south coords (lat) of nudging boxes [-180., 360.](deg). There are 6 boxes.
      real:: sitbox_nudg_n(6)= -999.              ! north coords (lat) of nudging boxes [-180., 360.](deg). There are 6 boxes.
      real:: sitbox_st_restore_time(6)= -999.       ! surface restore time scale of nudging boxes [-180., 360.](deg). There are 6 boxes.
      real:: sitbox_ut_restore_time(6)= -999.       ! upper ocean restore time scale of nudging boxes [-180., 360.](deg). There are 6 boxes.
      real:: sitbox_dt_restore_time(6)= -999.       ! deep ocean restore time scale of nudging boxes [-180., 360.](deg). There are 6 boxes.
      real:: sitbox_ss_restore_time(6)= -999.       ! surface restore time scale of nudging boxes [-180., 360.](deg). There are 6 boxes.
      real:: sitbox_us_restore_time(6)= -999.       ! upper ocean restore time scale of nudging boxes [-180., 360.](deg). There are 6 boxes.
      real:: sitbox_ds_restore_time(6)= -999.       ! deep ocean restore time scale of nudging boxes [-180., 360.](deg). There are 6 boxes.
      real:: sitbox_suv_restore_time(6)= -999.       ! surface restore time scale of nudging boxes [-180., 360.](deg). There are 6 boxes.
      real:: sitbox_uuv_restore_time(6)= -999.       ! upper ocean restore time scale of nudging boxes [-180., 360.](deg). There are 6 boxes.
      real:: sitbox_duv_restore_time(6)= -999.       ! deep ocean restore time scale of nudging boxes [-180., 360.](deg). There are 6 boxes.

    
      ! 4. 3-D Diecast Ocean Model
      !   4.1 I/O
      INTEGER, SAVE :: etopo_nres=1                 ! 1 for etopo1 (1min) (default), 2 for etopo2 (2min), and 5 for etop5 (5min)
      !   4.2 logical variables
      LOGICAL, SAVE :: locn  = .FALSE.              ! .true. for embedded 3-D ocean
      LOGICAL, SAVE :: locn_msg  = .FALSE.          ! .true. for writing embedded 3-D ocean fields
      LOGICAL, SAVE :: lopen_bound=.FALSE.          ! .TRUE. = open boundary condition (set depth at j=jos0,jon0 according to its physical depth
                                                    ! .FALSE.= closed boundary condition (set depth to be 0 for J=jos0,jon0
      LOGICAL, SAVE :: lall_straits=.TRUE.          ! .FALSE.= only open Gibraltar Strait
      LOGICAL, SAVE :: lstrict_channel=.TRUE.       ! .TRUE.  too strick, less vent in Gibraltar
                                                    ! .FALSE. too loose, open Central Ameican, Phillipine  
      !   4.3 other variables
      real:: ratio_dt_o2a=1.                 ! ratio_delta_time_and_ocn_dt (fractional)
      real:: ocn_domain_w  = 0.              ! west coords (lon) of embedded ocean [-180., 360.](deg).
      real:: ocn_domain_e  = 360.            ! east coords (lon) of embedded ocean [-180., 360.](deg).
      real:: ocn_domain_s  = -80.            ! south coords (lat) of embedded ocean [-180., 360.](deg).
      real:: ocn_domain_n  = 80.             ! north coords (lat) of embedded ocean [-180., 360.](deg).
      INTEGER, SAVE :: ocn_lon_factor=1             ! number of grid per dx in ECHAM resolution. 
      INTEGER, SAVE :: ocn_lat_factor=1             ! number of grid per dy in ECHAM resolution. 
      INTEGER, SAVE :: ocn_couple_option = 0        ! 0: put awust2,awvst2 and full-level T,s coupling with 3-D ocean, and T,u,v,s back to SIT.
!                                               ! 1: put afluxs(G0),awust2(ustr),awvst2(vstr),awfre(P-E) of echam to ocean, and all the T,u,v,s back to ECHAMS,
!                                               !    with sea ice correction (fluxiw rather than fluxs)
!                                               ! 2: Same as option 1, but also with sitwkh to ocean
!                                               ! 3: Same as option 2, but also with sitwkm to ocean
!                                               ! 4: Same as option 3, but without secruity number for diffusivity
!                                               ! 5: Same as option 1, but don't feedback anything back to sit.
!                                               ! 6: Same as option 5, but don't put awust2(ustr),awvst2(vstr) of echam to ocean.
!                                               ! 7: put wind stress, 1st layer sst and salinity of echam to ocean,
!                                               !    but don't feedback anything back to sit.
!                                               ! 8: put awust2,awvst2 of echam to ocean, but nothing back to sit.
!                                               ! 9: put wtb,wsb,awust2,awvst2,subfluxw,wsubsal of echam to ocean, and T,u,v,s back to SIT.
!                                               !10: Same as option 0, but skip TX,TY,TZ,SX,SY,SZ in ocn_stepon
!                                               !11: full-level (T,u,v,s) 2-way coupling with 3-D ocean, and T,u,v,s back to SIT.
!                                               !12: put nothing of echam to ocean, and nothing back to sit.
!                                               !13: Same as option 1, but T,S initialized by DIECAST
!                                               !14: Same as option 1, but without sea ice correction (fluxs rather than fluxiw)
!                                               !15: Same as option 10, + TX,SX 
!                                               !16: Same as option 10, + TY,SY
!                                               !17: Same as option 10, + TZ,SZ
      INTEGER :: high_current_killer = 4        ! 1: based on bottom vorticity: damp_high_vel(1:lnlon,1:lnlat,1:nh)=DAMP_high_current*MERGE( 1.,0.,ABS(vor_bottom(1:lnlon,1:lnlat,1:nh)).GT.VOR_cri )
                                                ! 2: based on bottom vertical velocity: damp_high_vel(i,j,ih)=DAMP_high_current*MERGE( 1.,0.,ABS(W(i,j,KB(i,j,ih),ih)).GT.W_cri ) 
                                                ! 3: based on bottom vorticity: damp_high_vel(1:lnlon,1:lnlat,1:nh)=DAMP_high_current*( ABS(vor_bottom(1:lnlon,1:lnlat,1:nh))/VOR_cri )
                                                ! 4: based on bottom vertical velocity: damp_high_vel(i,j,ih)=DAMP_high_current*( ABS(W(i,j,KB(i,j,ih),ih))/W_cri )
      real:: socn_restore_time =xmissing        ! surface (0 <= ~ <10 m) ocn grids restore time scale (s) (default: no nudging)
      real:: uocn_restore_time =xmissing        ! upper ocean (10 <= ~ <100 m) ocn grids restore time scale (s) (default: no nudging)
      real:: docn_restore_time =xmissing        ! deep (>=100 m) restore time scale (s)  (default: no nudging)
    !  
      INTEGER :: nobox_nudg= 0                  ! number of nudg squares in ocean grids (default = 0, maximun=6)
      real:: obox_restore_time=xmissing         ! restore time scale (s) in all depths for iop_ocnmask>0 grids (default: no nudging)
      INTEGER :: obox_nudg_flag= 0              ! 0: no nudging, 1: nudging inside square boxes, 2: nudging outside the square boxes
      real:: obox_nudg_w(6)= -999.              ! west coords (lon) of nudging boxes [-180., 360.](deg). There are 6 boxes.
      real:: obox_nudg_e(6)= -999.              ! east coords (lon) of nudging boxes [-180., 360.](deg). There are 6 boxes.
      real:: obox_nudg_s(6)= -999.              ! south coords (lat) of nudging boxes [-180., 360.](deg). There are 6 boxes.
      real:: obox_nudg_n(6)= -999.              ! north coords (lat) of nudging boxes [-180., 360.](deg). There are 6 boxes.
      real:: kocn_dm0z=1.                       ! =1 (default), ratio to near shore momentum diffusivity (m2/s) for preventing generating near-shore vortex
      INTEGER::  ncarpet=1                      ! number of coastal grid for carpet filter. =0, no carpet filter; =1 for N2 (default); =2 for N4; =3 for N6
      real:: kcsmag=1.                          ! =1, ratio to Smagorinsky horizontal diff coeff.
                                                ! Smagorinsky, J. General circulation experiments with the primitive equations, 
                                                !   I. The basic experiment. Monthly Weather Rev. 1963, 91, 99-164.
      real:: kalbw=1.                           ! =1, ratio to solar albedo over water by 
                                                ! Li J. J. Scinocca M. Lazare N. McFarlane K. von Salzen L. Solheim 2006
                                                !   Ocean Surface Albedo and Its Impact on Radiation Balance in Climate Models.
                                                !   J. Climate 19 6314??333.
                                                ! Taylor, J. P., J. M. Edwards, M. D. Glew, P. Hignett, and A. Slingo, 1996:
                                                !   Studies with a flexible new radiation code. II: Comparisons with aircraft short-wave
                                                !   observations. Quart. J. Roy. Meteor. Soc., 122, 839??61.
      real:: ck=0.1                             ! (0.1 in GASPAR ET AL., 1990) (v9.887)
      real:: ce=0.7                             ! (0.7 in Bougeault and Lacarrere, 1989)
      real:: Prw=1.                             ! =1, Prandtl number in water (sit_ocean.f90)                    
      real:: d0=0.03                            ! zero-displacement (m) (0.02,0.05) (v0.9871)
      real:: csl=-27.                           ! =-27 m (default), the present Caspian Sea Level (CSL) is about -27m and during the medieval time it was about -30m
      real:: por_min=0.1                        ! = minimum porisity for setting as ocean grid
      real:: csiced=0.                          ! =0., threshold sea ice depth, seaice[csiced]=0.5, seaice[2*csiced]=0.84
      REAL:: salti=0.9_dp                       ! the ratio of (ice salinity/sitws) [0,1]
 
! from mo_time_control

!!!  real           :: delta_time    = 0.0_dp ! distance of adjacent times
      real           :: time_step_len = 0.0 ! forecast time step,
                                               ! at beginning equal delta_time
                                               ! else 2*delta_time

  ! dt_initial  corresponds to the initial condition date
  ! dt_start    defines the start of an experiment
  !
!!!  TYPE(time_days),SAVE :: initial_date       ! should set at initial time from file
  !
      INTEGER, PARAMETER   :: INIT_STEP   = 0    ! initial time step
      INTEGER              :: dt_start(6) = 0    ! (runctl) start date of experiment
                                             ! meaning (yr, mo, dy, hr, mi, se)
!!!  TYPE(time_days),SAVE :: start_date         ! transformed start date
      LOGICAL              :: lsitstart    = .FALSE. ! .TRUE. for the first time step
    
      LOGICAL            :: lfirst_day = .TRUE.  ! .TRUE. during the first day
      LOGICAL            :: l2nd_day   = .TRUE.  ! .TRUE. during the first+second day
    
      INTEGER              :: dt_resume(6) = 0   ! user defined restart time
    
      INTEGER              :: dt_stop(6) = 0     ! (runctl) stop experiment here
                                                 ! meaning (yr, mo, dy, hr, mi, se)
!!!  TYPE(time_days),SAVE :: stop_date          ! transformed stop date
      LOGICAL              :: lbreak   = .FALSE. ! .TRUE. at end of one time segment
      LOGICAL              :: lstop    = .FALSE. ! .TRUE. during the last time step
      LOGICAL              :: labort   = .TRUE.  ! .TRUE. return error signal at end
!!!  LOGICAL              :: lwarmstart   = .FALSE.  ! .TRUE. read rerun file but modified with obs. ocean/atm/land data
      LOGICAL              :: lobs_ocn_rerun   = .FALSE.  ! .TRUE. read rerun file but modified with obs. ocean/atm/land data
      LOGICAL              :: lresume = .FALSE.  ! .TRUE. during rerun step
!!!  TYPE(time_days),SAVE ::  resume_date       ! transformed rerun date

!!!  TYPE(time_days),SAVE ::  previous_date     ! date at (time - delta_time)
!!!  TYPE(time_days),SAVE ::   current_date     ! date at (time)
!!!  TYPE(time_days),SAVE ::      next_date     ! date at (time + delta_time)

!!!  TYPE (time_days),SAVE :: radiation_date    ! date corresponding to rad_calc
      LOGICAL           :: l_orbvsop87 = .TRUE.  ! .TRUE. : orbit routine from vsop87

!ps
      LOGICAL, SAVE :: ldailysst = .TRUE. !   .true. for using daily SST and SIC

      LOGICAL, SAVE :: loutsit24    = .FALSE. !write wt,wu,wv,ws daily mean
      LOGICAL, SAVE :: lpre6hr_sit   = .FALSE. ! use lead 6 hours data of pre 6hr
      REAL :: outsitmean = -99.                 !write wt,wu,wv,ws every tau hours mean
      real:: sit_domain_w  = 0.             ! west coords (lon) of sit domain [0., 360.](deg).
      real:: sit_domain_e  = 360.           ! east coords (lon) of sit domain [0., 360.](deg).
      real:: sit_domain_s  = -30.           ! south coords (lat) of sit domain [-90., 90.](deg).
      real:: sit_domain_n  = 30.            ! north coords (lat) of sit domain [-90., 90.](deg).
      real:: sit_domain_extgrd  = 10.       ! extend degree of sit dimain deg).
      real:: ftrigsit= 0.                   ! start trigsit (return sst to atmospheric model) time (hr)
      LOGICAL,SAVE:: ltimeblending=.TRUE.  ! if true=start time blending for nudging
      INTEGER,SAVE:: timebl_option= 1       ! 0: nudging='0d' before timebl_start, gradually increase to nudging='nn' until timebl_allsit
                                            ! 1: nudging='0d' before timebl_start, after that, use sin to control nudging
      real:: timebl_start=0.                ! start time of time blending for nudging (day)
      real:: timebl_allsit=10.              ! always use sit restore_time for nudging after timebl_alsit (day)
      LOGICAL,SAVE:: lmixedlayer=.FALSE.    ! if lmixedlayer=t, read mixed layer depth data (read mixed_layer)
      real:: bathydepth=-200.               !
      LOGICAL,SAVE:: lsftobswt=.TRUE.      ! logical of shift SWT below 10m (10m=obswtb, delete difference
                                            !             between godas and obswtb data)
      integer:: outsitlev= 20               ! output sit level from 0 to outsitlev+1
!ps


  NAMELIST /sit_nml/          &
    loutsit24,                 &! write wt,wu,wv,ws daily mean
    lpre6hr_sit,               &! use lead 6 hours data of pre 6 hr
    outsitmean,                &! write wt,wu,wv,ws every tau hours mean
    outsitlev,                 &! output sit level from 0 to outsitlev+1
    sit_domain_w,              &! west coords (lon) of sit domain [-180., 360.](deg).
    sit_domain_e,              &! east coords (lon) of sit domain [-180., 360.](deg).
    sit_domain_s,              &! south coords (lat) of sit domain [-90., 90.](deg).
    sit_domain_n,              &! north coords (lat) of sit domain [-90., 90.](deg).
    sit_domain_extgrd,         &! extend degree of the sit domain (deg).
    ftrigsit,                  &! start trigsit (return sst to atmospheric model) forecast hour (default: 0. hr)
    ltimeblending,             &! if true=start time blending for nudging
    timebl_option,             &! time blending option for nudging
    timebl_start,              &! start time blending for nudging
    timebl_allsit,             &! always sit restore_time for nudging
    lmixedlayer,               &! if lmixedlayer=t, read mixed layer depth data (read mixed_layer)
    ldailysst,                 &! .true. for using daily SST and SIC
!ps    lrere,                     &! true for IC (Initial Condition) run, false for BC (Boundary Condition) run.
    lobs_ocn_rerun,            &! .TRUE. read rerun file but modified with obs. ocean/atm/land data
    !!! lsit,                      &! switch sit (i.e., calc. of vertical ocean temp. profile) on/off
    trigsit,                   &! coupling interval - trigger embedded ocean-model (default: 1 time setp)
    ltrigsit,                  &
    lsit_ice,                  &! switch on for turnning ice the ice_module of sit for tsi calculation (default: .true.)
    lsit_salt,                 &! switch salt routine in sit (i.e., calc. of salinity ) on/off
    zocn_option,                   &! godas z cord, fine in thermocline, 10 m apart within 100-225 m, but coarse (~200 m)
                                ! in deep-water formation depths (>500 m). (default: .false.)
    ocn_tlz,                   &! depth of bottom z-level (default =5800 m)
    ocn_k1,                    &! vertical dimension parameter (number of layer interfaces, equal
                                ! the number of layers or pressure levels, does not include ghost zone) (=40, default)
    lssst,                     &! .true. for turnon Skin SST, .false. for turnoff Skin SST.
    sit_ice_option,            &! ice option in sit (=0 (default), for coupling with vdiff semi-implicitlyoff; >=1, explicitly)
                                !   2: original SIT output, no security. It can be crashed in few time steps
    maskid,                    &! maskid (=0 (DIECAST only)) 
    lgodas,                    &! switch for using world ocean monthly atlas data (salinity and temp. profile) (n=12) on/off
    lsftobswt,                 &! shift godas 10 m swt to obs sst
    lamip,                     &! when lgodas=T, .FALSE. for using world ocean climatology monthly atlas data (salinity and temp. profile) (n=12) on/off
    lwoa0,                     &! switch for reading initial ocean salinity and temp. profile (n=1) on/off
    lwarning_msg,              &! switch lwarning_msg (i.e., warning message) on/off
    lsice_nudg,                &! spinup time (s) for nudging SIT/DIECAST (default: 1 yr=365.*86400.)
    lsit_lw,                   &! .true. for turning LW code for water in SIT.
    lcool_skin,                &! .true. for cool-skin parameterization (default=.FALSE.)        
    lsteady_TKE,               &! .true. = using calc_steady_TKE for TKE (default=.TRUE.)     
    lwave_breaking,            &! .true. for turning on wave_breaking TKE from ocean surface. (default=.FALSE.)
    ssit_restore_time,         &!  surface [0m,10m) other SIT grids restore time scale (s) (default: no nudging)
    usit_restore_time,         &!  upper [10m,100m) other SIT grids restore time scale (s) (default: no nudging)
    dsit_restore_time,         &!  deep (>=100 m) other SIT grids restore time scale (s)  (default: no nudging)
                                ! 86400 s : 1 d, -9.e33_dp : no nudging (default)
    ssits_restore_time,         &!  surface [0m,10m) SIT grids restore time scale for u, v (s) (default: no nudging)
    usits_restore_time,         &!  upper [10m,100m) SIT grids restore time scale for u, v (s) (default: no nudging)
    dsits_restore_time,         &!  deep (>=100 m) SIT grids restore time scale for u, v (s)  (default: no nudging)

    ssituv_restore_time,         &!  surface [0m,10m) SIT grids restore time scale for u, v (s) (default: no nudging)
    usituv_restore_time,         &!  upper [10m,100m) SIT grids restore time scale for u, v (s) (default: no nudging)
    dsituv_restore_time,         &!  deep (>=100 m) SIT grids restore time scale for u, v (s)  (default: no nudging)
                                ! 86400 s : 1 d, -9.e33_dp : no nudging (default)

    nsit_nudg,                 &! number of nudg squares in ocean grids (default = 0, maximun=6)
    sitbox_nudg_w,             &! west coords (lon) of nudging boxes [-180., 360.](deg). There are 6 boxes.
    sitbox_nudg_e,             &! east coords (lon) of nudging boxes [-180., 360.](deg). There are 6 boxes.
    sitbox_nudg_s,             &! south coords (lat) of nudging boxes [-180., 360.](deg). There are 6 boxes.
    sitbox_nudg_n,             &! north coords (lat) of nudging boxes [-180., 360.](deg). There are 6 boxes.
    sitbox_st_restore_time,    &! surface restore time scale of nudging boxes [-180., 360.](deg). There are 6 boxes.
    sitbox_ut_restore_time,    &! upper ocean restore time scale of nudging boxes [-180., 360.](deg). There are 6 boxes.
    sitbox_dt_restore_time,    &! deep ocean restore time scale of nudging boxes [-180., 360.](deg). There are 6 boxes.
    sitbox_ss_restore_time,    &! surface restore time scale of nudging boxes [-180., 360.](deg). There are 6 boxes.
    sitbox_us_restore_time,    &! upper ocean restore time scale of nudging boxes [-180., 360.](deg). There are 6 boxes.
    sitbox_ds_restore_time,    &! deep ocean restore time scale of nudging boxes [-180., 360.](deg). There are 6 boxes.
    sitbox_suv_restore_time,   &! surface restore time scale of nudging boxes [-180., 360.](deg). There are 6 boxes.
    sitbox_uuv_restore_time,   &! upper ocean restore time scale of nudging boxes [-180., 360.](deg). There are 6 boxes.
    sitbox_duv_restore_time,   &! deep ocean restore time scale of nudging boxes [-180., 360.](deg). There are 6 boxes.


    locaf,                     &! switch using climatological qflux adjustment on/off
    locaf0,                    &! switch using qflux adjustment on/off
    ocaf0_add,                  &!
!                              
    locn,                      &! switch embedded 3-D ocean on/off
    lopen_bound,               &! .TRUE. = open boundary condition (set depth at j=jos0,jon0 according to its physical depth
                                ! .FALSE.= closed boundary condition (set depth to be 0 for J=jos0,jon0
    lall_straits,              &! .FALSE.= only open Gibraltar Strait
    lstrict_channel,           &! .TRUE.  too strick, less vent in Gibraltar
                                ! .FALSE. too loose, open Central Ameican, Phillipine    
    etopo_nres,                &! 1 for etopo1 (1min), 2 for etopo2 (2min), and 5 for etop5 (5min)
    ocn_domain_w,              &! west coords (lon) of the ocean domain [-180., 360.](deg).
    ocn_domain_e,              &! east coords (lon) of the ocean domain [-180., 360.](deg).
    ocn_domain_s,              &! south coords (lat) of the ocean domain [-90., 90.](deg).
    ocn_domain_n,              &! north coords (lat) of the ocean domain [-90., 90.](deg).  
    ratio_dt_o2a,              &! ratio of the time step between ocn_dt and delta_time (ratio_dt_o2a=ocn_dt/delta_time) (default 1.) 
    ocn_couple_option,         &! switch full-level 2-way coupling with embedded 3-D ocean on/off
    high_current_killer,       &! 1: based on bottom vorticity
                                ! 2: based on bottom vertical velocity
    locn_msg,                  &! .true. for writing embedded 3-D ocean fields
    ocn_lon_factor,            &! number of grid per dx in ECHAM resolution
    ocn_lat_factor,            &! number of grid per dy in ECHAM resolution
    trigocn,                   &! coupling interval - trigger embedded ocean-model (default: 1 time setp)
    socn_restore_time,         &!  surface [0m,10m) ocn grids restore time scale (s) (default: no nudging)
    uocn_restore_time,         &!  upper [10m,100m) ocn grids restore time scale (s) (default: no nudging)
    docn_restore_time,         &!  deep (>=100 m) ocn grids restore time scale (s) (default: no nudging)
    nobox_nudg,                &! number of nudg squares in ocean grids (default = 0, maximun=6)
    obox_restore_time,         &!  restore time scale (s) in all depths for obox_mask>0 grids (default: -9.e33_dp) (no nudging)  
    obox_nudg_flag,            &! 0: no nudging, 1: nudging inside square boxes, 2: nudging outside the square boxes
    obox_nudg_w,               &! west coords (lon) of nudging boxes [-180., 360.](deg). There are 6 boxes.
    obox_nudg_e,               &! east coords (lon) of nudging boxes [-180., 360.](deg). There are 6 boxes.
    obox_nudg_s,               &! south coords (lat) of nudging boxes [-90., 90.](deg). There are 6 boxes.
    obox_nudg_n,               &! north coords (lat) of nudging boxes [-90., 90.](deg). There are 6 boxes.
    kocn_dm0z,                 &! =1, rato to near shore momentum diffusivity (m2/s) for preventing generating near-shore vortex
    ncarpet,                   &! number of coastal grid for carpet filter =0, no carpet filter; =1 for N2 (default); =2 for N4; =3 for N6  
    kcsmag,                    &! =1, ratio to Smagorinsky horizontal diff coeff.
                                ! Smagorinsky, J. General circulation experiments with the primitive equations, I. The basic experiment. Monthly Weather Rev. 1963, 91, 99-164.
    kalbw,                     &! =1, ratio to solar albedo over water by 
                                ! Li J. J. Scinocca M. Lazare N. McFarlane K. von Salzen L. Solheim 2006 Ocean Surface Albedo and Its Impact on Radiation Balance in Climate Models. J. Climate 19 6314??333.
                                ! Taylor, J. P., J. M. Edwards, M. D. Glew, P. Hignett, and A. Slingo, 1996: Studies with a flexible new radiation code. II: Comparisons with aircraft short-wave observations. Quart. J.
                                !   Roy. Meteor. Soc., 122, 839??61.
    ck,                        &! ck for TKE (sit_vdiff.f90)                    
    ce,                        &! ce for TKE (sit_vdiff.f90)                    
    Prw,                       &! Prandtl number in water (sit_vdiff.f90)                    
    d0,                        &! zero-displacement (m) (0.02,0.05) (sit_vdiff.f90) (v0.9871)
    csl,                       &! =-27, the present Caspian Sea Level (CSL) is about -27m and during the medieval time it was about -30m
    por_min,                   &! =0.1, minimum porisity for setting as ocean grid
    csiced,                    &! =0., threshold sea ice depth, seaice[csiced]=0.5, seaice[2*csiced]=0.84
    lasia,                     &! switch for etopo, new land surface data over rice paddy and Tibet on/off
                                ! http://www.ngdc.noaa.gov/mgg/global/global.html                  
    bathydepth                  ! bathy depth

  
      END MODULE mod_sit_control
