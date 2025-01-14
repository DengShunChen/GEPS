!#define SL_sedi
!#define new_saturation
!#define sat_predict
!#define use_declination
!#define use_cpm
!WRF:MODEL_LAYER:PHYSICS
!

MODULE module_mp_gsfcgce_3ice_nuwrf


!   USE     module_wrf_error
!NUWRF BEGIN
!   USE     module_dm
#if ( WRF_CHEM == 1)
!JJS 20110525     vvvvv
! for inline Gocart coupling
   USE     module_gocart_coupling
!JJS 20110525     ^^^^^
#endif
!NUWRF END
!   USE module_utility, ONLY: WRFU_Clock, WRFU_Alarm
!   USE module_domain, ONLY : HISTORY_ALARM, Is_alarm_tstep
   USE module_mp_radar

!   LOGICAL, EXTERNAL :: wrf_dm_on_monitor

!JJS 20140117 vvvvv
   PRIVATE   ! privatize all variables/subroutines in this module excepting public parameter below
   PUBLIC :: gsfcgce_3ice_nuwrf
!JJS 20140117 ^^^^^

!JJS 1/3/2008     vvvvv

!  common /bt/
   REAL,    PRIVATE ::          rd1,  rd2,   al,   cp

!  common /cont/
   REAL,    PRIVATE ::          c38, c358, c610, c149,             &
                               c879, c172, c409,  c76,             &
                               c218, c580, c141
!  common /b3cs/
   REAL,    PRIVATE ::           ag,   bg,   as,   bs,             &
                                 aw,   bw,  bgh,  bgq,             &
                                bsh,  bsq,  bwh,  bwq

!  common /size/
   REAL,    PRIVATE ::          tnw,  tns,  tng,                   &
                               roqs, roqg, roqr

!  common /bterv/
!   REAL,    PRIVATE ::           zrc,  zgc,  zsc,                 &
!                                 vrc,  vgc,  vsc

!  common /rterv/
   REAL,    PRIVATE ::           zrc,  zgc,  zsc, vrc,             &  
                                vrc0, vrc1, vrc2, vrc3,            &
                             	 vgc,  vsc

   REAL,    PRIVATE ::         draimax

!  common /bsnw/
!   REAL,    PRIVATE ::          alv,  alf,  als,   t0,   t00,     &
!                                avc,  afc,  asc,  rn1,  bnd1,     &
!                                rn2, bnd2,  rn3,  rn4,   rn5,     &
!                                rn6,  rn7,  rn8,  rn9,  rn10,     &
!                              rn101,rn10a, rn11,rn11a,  rn12
!
!   REAL,    PRIVATE ::         rn14, rn15,rn15a, rn16,  rn17,     &
!                              rn17a,rn17b,rn17c, rn18, rn18a,     &
!                               rn19,rn19a,rn19b, rn20, rn20a,     &
!                              rn20b, bnd3, rn21, rn22,  rn23,     &
!                              rn23a,rn23b, rn25,rn30a, rn30b,     &
!                              rn30c, rn31, beta, rn32

!  common /bsnw/
   REAL,    PRIVATE ::             bnd1, rn11a

   REAL,    PRIVATE ::          rn17, rn19b,                       &
                               	bnd3, rn23a,                   &
                               rn23b, rn30b,                       &
                               rn30c 

!  common /rsnw/
    REAL,    PRIVATE ::          alv,   alf,   als,    t0,   t00,     &
                                 avc,   afc,   asc,   rn1,   rn2,     &
                                bnd2,   rn3,   rn4,   rn5,  rn50,     &
                                rn51,  rn52,  rn53,   rn6,  rn60,     &
                                rn61,  rn62,  rn63,   rn7,   rn8,     &
                                 rn9,  rn10, rn101, rn102, rn10a,     &
                               rn10b, rn10c,  rn11,  rn12,  rn14,     &
                                rn15, rn15a,  rn16, rn171, rn172,     &
                               rn17a, rn17b, rn17c,  rn18, rn18a,     &
                                rn19, rn191, rn192, rn19a,  rn20,     &
                               rn20a, rn20b,  rn30, rn30a,  rn21,     &
                               bnd21,  rn22,  rn23, rn231, rn232,     &
                                rn25,  rn31,  beta,  rn32,  rn33,     &
                               rn331, rn332,  rn34,  rn35,rnn30a,     &
                               rnn191,rnn192

    REAL,    PRIVATE ::          ami50, ami40

   REAL,    PRIVATE, DIMENSION( 31 ) ::    rn12a, rn12b, rn13, rn25a

   REAL,    PRIVATE, DIMENSION( 31 ) ::    BergCon1,  BergCon2,       &
                                           BergCon3,  BergCon4

   ! critical q of hydrometeor characteristics (eg. fall speed, radius)
   REAL,    PRIVATE :: cwmin, cimin, crmin, csmin, cgmin

!
   REAL,    PRIVATE, DIMENSION( 31 )  ::      aa1,  aa2
   DATA aa1/.7939e-7, .7841e-6, .3369e-5, .4336e-5, .5285e-5,         &
           .3728e-5, .1852e-5, .2991e-6, .4248e-6, .7434e-6,          &
           .1812e-5, .4394e-5, .9145e-5, .1725e-4, .3348e-4,          &
           .1725e-4, .9175e-5, .4412e-5, .2252e-5, .9115e-6,          &
           .4876e-6, .3473e-6, .4758e-6, .6306e-6, .8573e-6,          &
           .7868e-6, .7192e-6, .6513e-6, .5956e-6, .5333e-6,          &
           .4834e-6/
   DATA aa2/.4006, .4831, .5320, .5307, .5319,                        &
           .5249, .4888, .3894, .4047, .4318,                         &
           .4771, .5183, .5463, .5651, .5813,                         &
           .5655, .5478, .5203, .4906, .4447,                         &
           .4126, .3960, .4149, .4320, .4506,                         &
           .4483, .4460, .4433, .4413, .4382,                         &
           .4361/

   real, dimension(0:120) :: itble       ! deposition growth coefficients
   data itble /1.000000,0.979490,0.959401,0.939723,0.920450,          &
               0.899498,0.879023,0.857038,0.833681,0.810961,          &
               0.783430,0.755092,0.703072,0.537032,0.467735,          &
               0.524807,0.630957,0.812831,1.096478,1.479108,          &
               1.905461,2.089296,2.290868,2.398833,2.454709,          &
               2.426610,2.371374,2.290868,2.137962,1.995262,          &
               1.862087,1.737801,1.621810,1.513561,1.396368,          &
               1.288250,1.188502,1.096478,1.000000,0.922571,          &
               0.851138,0.785236,0.724436,0.668344,0.616595,          &
               0.575440,0.537032,0.501187,0.467735,0.436516,          &
               0.407380,0.380189,0.354813,0.331131,0.316228,          &
               0.301995,0.291743,0.285102,0.281838,0.278612,          &
               0.275423,0.278612,0.281838,0.285102,0.291743,          &
               0.298538,0.309030,0.319890,0.331131,0.346737,          &
               0.367282,0.393550,0.426580,0.457088,0.489779,          &
               0.524807,0.562341,0.609537,0.660693,0.716143,          &
               0.785236,0.860994,0.954993,1.047129,1.148154,          &
               1.258925,1.380384,1.496236,1.603245,1.698244,          &
               1.778279,1.840772,1.883649,1.905461,1.905461,          &
               1.883649,1.862087,1.840772,1.798871,1.737801,          &
               1.698244,1.640590,1.584893,1.548817,1.513561,          &
               1.475707,1.452112,1.428894,1.412538,1.393157,          &
               1.377209,1.361445,1.348963,1.336596,1.327394,          &
               1.318257,1.309182,1.303167,1.294196,1.288250,1.279381/
   real, private, parameter :: thrd = 1./3.
!JJS 1/3/2008     ^^^^^

!+---+-----------------------------------------------------------------+
!..The following 6 variables moved here to facilitate reflectivity
!.. calculation similar to other MP schemes, because when they get
!.. declared later in the code (now commented out), it makes things
!.. more difficult to integreate with the radar code.
!.. Values will be defined in subroutine all_flux --- JJS 20140116
      REAL ::     xnor, xnos, xnoh, xnog  
      REAL ::     rhowater, rhosnow
      REAL ::     rhohail, rhograul
!      REAL    , PARAMETER ::     xnor = 8.0e6
!      REAL    , PARAMETER ::     xnos = 1.6e7
!      REAL    , PARAMETER ::     xnoh = 2.0e5
!      REAL    , PARAMETER ::     xnog = 4.0e6
!      REAL    , PARAMETER ::     rhohail = 917.
!      REAL    , PARAMETER ::     rhograul = 400.
!+---+-----------------------------------------------------------------+

!.. for terminal velocity and fall_flux :
      real, private :: consta, constb, constc, constd, o6, cdrag, &
                       abar, bbar, rhoe_s

CONTAINS

!-------------------------------------------------------------------
!  NASA/GSFC GCE
!  Tao et al, 2001, Meteo. & Atmos. Phy., 97-137
!-------------------------------------------------------------------
!  SUBROUTINE gsfcgce(  th, th_old,                                &
  SUBROUTINE gsfcgce_3ice_nuwrf( th                                &
                      ,qv, ql                                      &
                      ,qr, qi                                      &
                      ,qs                                          &
!                      ,qvold, qlold                               &
!                      ,qrold, qiold                               &
!                      ,qsold                                      &
                      ,rho, pii, p, dt_in, z                       &
                      ,ht, dz8w, grav, w                           &
!                      ,rhowater, rhosnow                           &
                      ,itimestep, xlat, sdec, xland                &
!                      ,ids,ide, jds,jde, kds,kde                   & ! domain dims
                      ,ims,ime, jms,jme, kms,kme                   & ! memory dims
                      ,its,ite, jts,jte, kts,kte                   & ! tile   dims
                      ,rainnc, icenc, snownc, graupelnc, sr        &
!                      ,rainncv, snowncv, graupelncv                &
!                      ,f_qg, qg, pgold                            &
                      ,f_qg, qg                                    &
                      ,ihail, ice2                                 &
!NUWRF BEGIN
#ifdef EXT_DIAG
                      ,refl_10cm, diagflag, do_radar_ref           &
!                      ,physc, physe, physd, physs, physm, physf    &
!                      ,acphysc, acphyse, acphysd, acphyss, acphysm, acphysf &
                      ,preci3d, precs3d, precg3d, precr3d          &
#endif
                      ,refc, refr, refi, refs, refg                & ! cloud effective radius
#if ( WRF_CHEM == 1)
!JJS 20110525     vvvvv
                      ,aero, icn_diag, nc_diag, gid                &
!JJS 20110525     ^^^^^
! EMK
                      ,chem_opt                                    &
                      ,gsfcgce_gocart_coupling                     &
#endif 
!NUWRF END
                                                                   )

!-------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------
!
! JJS 2/15/2005
!
!  INTEGER,      INTENT(IN   )    ::   ids,ide, jds,jde, kds,kde
  INTEGER,      INTENT(IN   )    ::   ims,ime, jms,jme, kms,kme , &
                                      its,ite, jts,jte, kts,kte 
  INTEGER,      INTENT(IN   )    ::   itimestep, ihail, ice2

  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),                 &
        INTENT(INOUT) ::                                          &
                                                              th, &
                                                              qv, &
                                                              ql, &
                                                              qr, &
                                                              qi, &
                                                              qs, &
                                                              qg

!NUWRF BEGIN
#if ( WRF_CHEM == 1)
! JJS 20110525 vvvvv
! for inline Gocart coupling
  INTEGER, PARAMETER :: num_go = 14  ! number of the gocart aerosol species
  REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_go), intent(in) :: aero
  REAL, DIMENSION( ims:ime, kms:kme, jms:jme), intent(out) :: icn_diag, nc_diag
  INTEGER,      INTENT(IN   )    ::   gid
! JJS 20110525 ^^^^^
  integer,intent(in) :: chem_opt ! EMK
  integer,intent(in) :: gsfcgce_gocart_coupling ! EMK
#endif
!NUWRF END
!
  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),                 &
        INTENT(IN   ) ::                                          &
!                                                         th_old, &
!                                                          qvold, &
!                                                          qlold, &
!                                                          qrold, &
!                                                          qiold, &
!                                                          qsold, &
!                                                          qgold, &
                                                             rho, &
                                                             pii, &
                                                               p, &
                                                            dz8w, &
                                                               z, &
                                                               w

  REAL, DIMENSION( ims:ime , jms:jme ),                           &
        INTENT(INOUT) ::    rainnc, icenc, snownc, graupelnc, sr
!  REAL, DIMENSION( ims:ime , jms:jme ),                           &
!        INTENT(INOUT) ::           rainncv, snowncv, graupelncv

!JJS 20140225   for calculation of effective radius of cloud species
  REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN)   :: XLAND
  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),                 &
        INTENT(INOUT) ::                               refc,      &
                                                       refr,      &
                                                       refi,      &
                                                       refs,      &
                                                       refg
!JJS 20140225  ^^^^^

!+---+-----------------------------------------------------------------+
#ifdef EXT_DIAG
  REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT)::           &  ! GT
                                                       refl_10cm
  LOGICAL, OPTIONAL, INTENT(IN) :: diagflag
  INTEGER, OPTIONAL, INTENT(IN) :: do_radar_ref
#endif
!+---+-----------------------------------------------------------------+

  REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN) ::       ht

  REAL, INTENT(IN   ) ::                                   dt_in, &
                                                            grav
!                                                        rhowater, &
!                                                         rhosnow
  REAL, INTENT(IN   ) :: xlat
  REAL, INTENT(IN   ) :: sdec  ! sine of solar declination angle
  LOGICAL, INTENT(IN), OPTIONAL :: F_QG

!  LOCAL VAR

!
  INTEGER ::  itaobraun, istatmin, new_ice_sat, id
  INTEGER ::  improve

  INTEGER :: i, j, k
  INTEGER :: iskip, ih, icount, ibud, i24h 
  REAL    :: hour
  REAL    , PARAMETER :: cmin=1.e-20
  REAL    :: dth, dqv, dqrest, dqall, dqall1, rhotot, a1, a2 
 
  LOGICAL :: flag_qg

#ifdef EXT_DIAG
  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),                 &
        INTENT(INOUT) ::                                          &
!             physc, physe, physd, physs, physm, physf,                  &
!             acphysc, acphyse, acphysd, acphyss, acphysm, acphysf,      &
             preci3d, precs3d, precg3d, precr3d
#endif
 
  ! for the conversion of q :
  real, dimension(its:ite,kts:kte,jts:jte) :: qtot

  ! for sub-cycle time steps :
  integer :: n, ntimes
  real :: mp_time, dts

!NUWRF BEGIN
#if ( WRF_CHEM == 1)
! EMK 2011/08/25
  logical,save :: first_makelut = .true. ! Initial value, changed later.
#endif
!NUWRF END

!+---+-----------------------------------------------------------------+

      INTEGER :: NCALL=0

!+---+-----------------------------------------------------------------+


!
!c  ihail = 0    for graupel, for tropical region
!c  ihail = 1    for hail, for mid-lat region

! itaobraun: 0 for Tao's constantis, 1 for Braun's constants
!c        if ( itaobraun.eq.1 ) --> betah=0.5*beta=-.46*0.5=-0.23;   cn0=1.e-6
!c        if ( itaobraun.eq.0 ) --> betah=0.5*beta=-.6*0.5=-0.30;    cn0=1.e-8
   itaobraun = 0 ! Set to zero in NUWRF

! Use Steve's new improvement   9/18/2009

    improve = 3
!   improve = -20         ! to use the original codes

!c  ice2 = 0    for 3 ice --- ice, snow and graupel/hail
!c  ice2 = 1    for 2 ice --- ice and snow only
!c  ice2 = 2    for 2 ice --- ice and graupel only, use ihail = 0 only
!c  ice2 = 3    for 0 ice --- no ice, warm only

!  if (ice2 .eq. 2) ihail = 0

  i24h=nint(86400./dt_in)
  if (mod(itimestep,i24h).eq.1) then
!     write(6,*) 'ihail=',ihail,'  ice2=',ice2
     if (ice2.eq.0) then
!        write(6,*) 'Running 3-ice scheme in GSFCGCE with'
        if (ihail.eq.0) then 
!           write(6,*) '     ice, snow and graupel'
        else if (ihail.eq.1) then
!                write(6,*) '     ice, snow and hail'
        else
!             write(6,*) 'ihail has to be either 1 or 0'
!             call wrf_error_fatal ('gsfcgce_2ice in namelist.input has to be either 1 or 0')
             stop
        endif !ihail
     else if (ice2.eq.1) then
!             write(6,*) 'Running 2-ice scheme in GSFCGCE with'
!             write(6,*) '     ice and snow'
     else if (ice2.eq.2) then
!             write(6,*) 'Running 2-ice scheme in GSFCGCE with'
!             write(6,*) '     ice and graupel'
     else if (ice2.eq.3) then
!             write(6,*) 'Running warm rain only scheme in GSFCGCE without any ice'
     else
!             write(6,*) 'gsfcgce_2ice in namelist.input has to be 0, 1, 2, or 3'
!             call wrf_error_fatal ('gsfcgce_2ice in namelist.input has to be 0, 1, 2, or 3')
             stop
     endif !ice2
  endif !itimestep

!c  new_ice_sat = 0, 1, 2, or 3 
    new_ice_sat = 3 ! Set to 3 in NUWRF

!c istatmin
    istatmin = 180

!c id = 0  without in-line staticstics
!c id = 1  with in-line staticstics
    id = 0

!c ibud = 0 no calculation of dth, dqv, dqrest and dqall
!c ibud = 1 yes
    ibud = 0

!jjs   dt=dt_in
!jjs   rhoe_s=1.29
!
!   IF (P_QI .lt. P_FIRST_SCALAR .or. P_QS .lt. P_FIRST_SCALAR) THEN
!      CALL wrf_error_fatal3 ( "module_mp_lin.b" , 130 ,  'module_mp_lin: Improper use of Lin et al scheme; no ice phase. Please chose another one.')
!   ENDIF

!NUWRF BEGIN
#if ( WRF_CHEM == 1)
! JJS 20110525 vvvvv
! Create CCN and IN candidate look-up table. Only need to run once at the first timestep of each run
!    if (itimestep .eq. 1) call makelut_ccn_icn
! EMK 2011/08/25
      if (mod(itimestep,i24h).eq.1)   &
!         print *, 'in mp_gsfcgce_3ice_nuwrf, chem_opt, gsfcgce_gocart_coupling  = ', chem_opt, gsfcgce_gocart_coupling
      if (first_makelut) then
         ! EMK...Only execute when GOCART and coupling turned on.
         if ( (chem_opt == 300 .or. chem_opt == 301 .or. &
               chem_opt == 302 .or. chem_opt == 303) .and. &
               (gsfcgce_gocart_coupling == 1)) then
!            if ( wrf_dm_on_monitor() ) print *,'Goddard MP coupled with Gocart Aerosol'
            call makelut_ccn_icn
            first_makelut = .false. ! This value is saved for subsequent calls.
         end if
      end if
! JJS 20110525 ^^^^^
#endif

   ! convert specific values of q to mixing ratios :
   qtot = 0.
   do j = jts, jte
   do i = its, ite
   do k = kts, kte
     qtot(i,k,j) = qv(i,k,j) + ql(i,k,j) + qr(i,k,j)             &
                   + qi(i,k,j) + qs(i,k,j) + qg(i,k,j)
     qv(i,k,j) = qv(i,k,j)/(1.-qtot(i,k,j))
     ql(i,k,j) = ql(i,k,j)/(1.-qtot(i,k,j))
     qr(i,k,j) = qr(i,k,j)/(1.-qtot(i,k,j))
     qi(i,k,j) = qi(i,k,j)/(1.-qtot(i,k,j))
     qs(i,k,j) = qs(i,k,j)/(1.-qtot(i,k,j))
     qg(i,k,j) = qg(i,k,j)/(1.-qtot(i,k,j))

     ! qtot must be saved as "mixing ratio" :
     qtot(i,k,j) = qv(i,k,j) + ql(i,k,j) + qr(i,k,j)             &
                   + qi(i,k,j) + qs(i,k,j) + qg(i,k,j)
   enddo
   enddo
   enddo

   ! set up constants used internally in GCE
   call consat_s (ihail, itaobraun, improve)

   ! set sub-cycle time step :
   mp_time = 300.                !standard sub-cycle time step
   ntimes = 1                    !number of sub-cycles
   ntimes = max (ntimes, int (dt_in / min(dt_in, mp_time)))
   dts = dt_in / real (ntimes)   !real sub-cycle time step

   do n = 1, ntimes

   ! calculte fallflux and precipiation in MKS system
   call fall_flux(    dts, qv, qr, qi, qs, qg, p,             &
                      rho, th, pii, z, dz8w, ht,              &
                      grav, itimestep,                        &
                      rainnc, icenc, snownc, graupelnc, sr,   &
!                      rainncv, snowncv, graupelncv,           &
#ifdef EXT_DIAG
                      preci3d, precs3d, precg3d, precr3d,     &
#endif
                      xland,                                  &
                      ihail, ice2, improve,                   &
                      ims,ime, jms,jme, kms,kme,              & ! memory dims
                      its,ite, jts,jte, kts,kte               ) ! tile   dims

#ifdef EXT_DIAG
      ! EMK NUWRF...Moved this WRF radar reflectivity initialization to after
      ! fall_flux, as the rhohail and rhograul variables are set in that
      ! subroutine. 
      IF (NCALL .EQ. 0) THEN
!..Set these variables needed for computing radar reflectivity.  These
!.. get used within radar_init to create other variables used in the
!.. radar module.
         xam_r = 3.14159*rhowater/6.
         xbm_r = 3.
         xmu_r = 0.
         xam_s = 3.14159*rhosnow/6.
         xbm_s = 3.
         xmu_s = 0.
         if (ihail .eq. 1) then
            xam_g = 3.14159*rhohail/6.
         else
            xam_g = 3.14159*rhograul/6.
         endif
         xbm_g = 3.
         xmu_g = 0.

         call radar_init
         NCALL = 1
      ENDIF
!-----------------------------------------------------------------------

      IF (NCALL .EQ. 0) THEN
!..Set these variables needed for computing radar reflectivity.  These
!.. get used within radar_init to create other variables used in the
!.. radar module.
         xam_r = 3.14159*rhowater/6.
         xbm_r = 3.
         xmu_r = 0.
         xam_s = 3.14159*rhosnow/6.
         xbm_s = 3.
         xmu_s = 0.
         if (ihail .eq. 1) then
            xam_g = 3.14159*rhohail/6.
         else
            xam_g = 3.14159*rhograul/6.
         endif
         xbm_g = 3.
         xmu_g = 0.

         call radar_init
         NCALL = 1
      ENDIF
#endif

!c Negative values correction

!   iskip = 1
! 
!   if (iskip.eq.0) then
!      call negcor(qv,rho,dz8w,ims,ime,jms,jme,kms,kme, &
!                           itimestep,1,             &
!                           its,ite,jts,jte,kts,kte)
!      call negcor(ql,rho,dz8w,ims,ime,jms,jme,kms,kme, &
!                           itimestep,2,             &
!                           its,ite,jts,jte,kts,kte)
!      call negcor(qr,rho,dz8w,ims,ime,jms,jme,kms,kme, &
!                           itimestep,3,             &
!                           its,ite,jts,jte,kts,kte)
!      call negcor(qi,rho,dz8w,ims,ime,jms,jme,kms,kme, &
!                           itimestep,4,             &
!                           its,ite,jts,jte,kts,kte)
!      call negcor(qs,rho,dz8w,ims,ime,jms,jme,kms,kme, &
!                           itimestep,5,             &
!                           its,ite,jts,jte,kts,kte)
!      call negcor(qg,rho,dz8w,ims,ime,jms,jme,kms,kme, &
!                           itimestep,6,             &
!                           its,ite,jts,jte,kts,kte)
!!   else if (mod(itimestep,i24h).eq.1) then
!!      print *,'no neg correction in mp at timestep=',itimestep
!   endif ! iskip


   ! microphysics in GCE
   call SATICEL_S( dts, IHAIL, itaobraun, ICE2, istatmin,        &
                   new_ice_sat, id, improve, xlat, sdec,         &
                   th, qv, ql, qr,                               &
                   qi, qs, qg,                                   &
                   rho, pii, p, w,                               &
                   itimestep, xland,                             & 
                   refc, refr, refi, refs, refg,                 & ! cloud effective radius
                   ims,ime, jms,jme, kms,kme,                    & ! memory dims
                   its,ite, jts,jte, kts,kte                     & ! tile   dims
#ifdef EXT_DIAG
                   ,refl_10cm, diagflag, do_radar_ref,           & ! GT added for reflectivity calcs
                   physc, physe, physd, physs, physm, physf,     &
                   acphysc, acphyse, acphysd, acphyss, acphysm, acphysf &
#endif
#if ( WRF_CHEM == 1)
                   ,aero, icn_diag, nc_diag, gid,                &
                   chem_opt,                                     &
                   gsfcgce_gocart_coupling                       &
#endif
                   )
   enddo  !end of do n=1,ntimes

   ! convert mixing values of q back to specific values :
   do j = jts, jte
   do i = its, ite
   do k = kts, kte
     qv(i,k,j) = qv(i,k,j)/(1.+qtot(i,k,j))
     ql(i,k,j) = ql(i,k,j)/(1.+qtot(i,k,j))
     qr(i,k,j) = qr(i,k,j)/(1.+qtot(i,k,j))
     qi(i,k,j) = qi(i,k,j)/(1.+qtot(i,k,j))
     qs(i,k,j) = qs(i,k,j)/(1.+qtot(i,k,j))
     qg(i,k,j) = qg(i,k,j)/(1.+qtot(i,k,j))
   enddo
   enddo
   enddo

  END SUBROUTINE gsfcgce_3ice_nuwrf

  SUBROUTINE fall_flux ( dt, qv, qr, qi, qs, qg, p,           &
                      rho, th, pi_mks, z, dz8w, topo,         &
                      grav, itimestep,                        &
                      rainnc, icenc, snownc, graupelnc, sr,   &
!                      rainncv, snowncv, graupelncv,           &
#ifdef EXT_DIAG
                      preci3d, precs3d, precg3d, precr3d,     &
#endif
                      xland,                                  &
                      ihail, ice2, improve,                   &
                      ims,ime, jms,jme, kms,kme,              & ! memory dims
                      its,ite, jts,jte, kts,kte               ) ! tile   dims
!-----------------------------------------------------------------------
! adopted from Jiun-Dar Chern's codes for Purdue Regional Model
! adopted by Jainn J. Shi, 6/10/2005
! modified by Goddard 7/24/2010
! modified by Tao 11/12/2010
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER, INTENT(IN   )               :: ihail, ice2, improve,       &
                                          ims,ime, jms,jme, kms,kme,  &
                                          its,ite, jts,jte, kts,kte 
  INTEGER, INTENT(IN   )               :: itimestep
  REAL,    DIMENSION( ims:ime , kms:kme , jms:jme ),                  &
           INTENT(INOUT)               :: qv, qr, qi, qs, qg      
  REAL,    DIMENSION( ims:ime , kms:kme , jms:jme ),                  &
           INTENT(IN)                  :: th, pi_mks      

  REAL,    DIMENSION( ims:ime , jms:jme ),                            &
           INTENT(INOUT) :: rainnc, icenc, snownc, graupelnc, sr
!  REAL,    DIMENSION( ims:ime , jms:jme ),                            &
!           INTENT(INOUT)               :: rainncv, snowncv, graupelncv
  REAL,    DIMENSION( ims:ime , kms:kme , jms:jme ),                  &
           INTENT(IN   )               :: rho, z, dz8w, p

  REAL,    INTENT(IN   )               :: dt, grav


  REAL,    DIMENSION( ims:ime , jms:jme ),                            &
           INTENT(IN   )               :: topo, xland

#ifdef EXT_DIAG
  REAL,    DIMENSION( ims:ime , kms:kme , jms:jme ),                  &
           INTENT(OUT)                 :: preci3d, precs3d, precg3d, precr3d
#endif

! temperary vars
 
  REAL,    DIMENSION( kts:kte )           :: fv
  REAL                                    :: tmp1, term0
  REAL                                :: pptrain, pptsnow,        &
                                         pptgraul, pptice, pptall
  REAL,    DIMENSION( kts:kte )       :: qvz, qrz, qiz, qsz, qgz, &
                                         zz, dzw, prez, rhoz,     &
                                         orhoz
  REAL,    DIMENSION( kts:kte )       :: rsed, ised, ssed, gsed
  REAL,    DIMENSION( kts:kte )       :: thz, piz, tz

   INTEGER                    :: k, i, j
!

  REAL, DIMENSION( kts:kte )    :: vtr, vts, vtg, vti

  REAL                          :: dtb, pi, gambp4, gamdp4, gam4pt5, gam4bbar

!NUWRF BEGIN
! New local variable
  INTEGER                       :: ic, igce
  REAL                          :: y1, vr, vs, vg
  REAL                          :: tslopes, tslopeg, vgcr, vscf
  REAL                          :: tair, tairc, fexp
  REAL                          :: ftns, ftnsQ, ftng, ftngQ
  REAL                          :: const_vt, const_d, const_m, bb1, bb2
  REAL,    DIMENSION(7)         :: aice, vice
  REAL                          :: cmin, r00
  REAL                          :: ftns0, ftng0

  DATA tslopes/0./, tslopeg/0./
!  DATA const_vt//, const_d//, const_m//
  DATA aice/1.e-6, 1.e-5, 1.e-4, 1.e-3, 0.01, 0.1, 1./  ! g/m**3
  DATA vice/5.,15.,30.,35.,40.,45.,50./  ! cm/s
  DATA igce/0/ 
  DATA ftns/1./, ftng/1./

!NUWRF END

  REAL    , PARAMETER :: p0 = 1.0e5

! for terminal velocity flux
  INTEGER                       :: min_q, max_q
  REAL                          :: t_del_tv, del_tv, flux, fluxin, fluxout ,tmpqrz
  LOGICAL                       :: notlast
#ifdef SL_sedi
  integer                       :: n,ntimes
  real                          :: dtcfl,precip
  real , dimension(kts:kte)     :: qden
#endif
  real , dimension(its:ite,kts:kte,jts:jte) :: vtr3d,vts3d,vtg3d,vti3d

  if (improve.eq.3) igce=1

!  if (itimestep.eq.1) then
!     write(6, *) 'in fall_flux'
!     write(6, *) 'igce=',igce
!     write(6, *) 'improve=',improve
    
!     write(6, *) 'ims=', ims, '  ime=', ime
!     write(6, *) 'jms=', jms, '  jme=', jme
!     write(6, *) 'kms=', kms, '  kme=', kme
!     write(6, *) 'its=', its, '  ite=', ite
!     write(6, *) 'jts=', jts, '  jte=', jte
!     write(6, *) 'kts=', kts, '  kte=', kte
!     write(6, *) 'dt=', dt
!     write(6, *) 'ihail=', ihail
!     write(6, *) 'ICE2=', ICE2
!     write(6, *) 'dt=', dt
!   endif 

!-----------------------------------------------------------------------
!  This program calculates precipitation fluxes due to terminal velocities.
!-----------------------------------------------------------------------

   dtb=dt
   pi=acos(-1.)

!  Gamma function
   gambp4=gammagce(constb+4.)   !ga4b
   gamdp4=gammagce(constd+4.)   !ga4d
   gam4pt5=gammagce(4.5)
   gam4bbar=gammagce(4.+bbar)
!
      cmin=1.e-10
!
!***********************************************************************
! Calculate precipitation fluxes due to terminal velocities.
!***********************************************************************
!
!- Calculate termianl velocity (vt?)  of precipitation q?z
!- Find maximum vt? to determine the small delta t

 j_loop:  do j = jts, jte
 i_loop:  do i = its, ite

   do k = kts, kte
#ifdef EXT_DIAG
      preci3d(i,k,j)=0.
      precs3d(i,k,j)=0.
      precg3d(i,k,j)=0.
      precr3d(i,k,j)=0.
#endif
      ised(k)=0.
      ssed(k)=0.
      gsed(k)=0.
      rsed(k)=0.
      vtr3d(i,k,j)=0.
      vts3d(i,k,j)=0.
      vtg3d(i,k,j)=0.
      vti3d(i,k,j)=0.
   end do

   pptrain = 0.
   pptsnow = 0.
   pptgraul = 0.
   pptice  = 0.

   ! in MKS system
   do k = kts, kte
      qvz(k)=qv(i,k,j)
      qrz(k)=qr(i,k,j)
      qiz(k)=qi(i,k,j)
      qsz(k)=qs(i,k,j)
      rhoz(k)=rho(i,k,j)
      thz(k)=th(i,k,j)
      piz(k)=pi_mks(i,k,j)
      orhoz(k)=1./rhoz(k)
      prez(k)=p(i,k,j)
      fv(k)=sqrt(rhoe_s/rhoz(k))
!      fv(k)=sqrt(rho(i,1,j)/rhoz(k))
      zz(k)=z(i,k,j)
      dzw(k)=dz8w(i,k,j)
      tz(k)=thz(k)*piz(k)
   enddo !k

   IF (ice2.eq.0 .or. ice2.eq.2) THEN
      DO k = kts, kte
         qgz(k)=qg(i,k,j)
      ENDDO
   ELSE
      DO k = kts, kte
         qgz(k)=0.
      ENDDO
   ENDIF


!
!-- rain
!
#ifdef SL_sedi
    ntimes = 1
    dtcfl = dtb/real(ntimes)

    do n = 1 , ntimes
       vtr(:) = 0.
       precip = 0.
       call semi_sedi('qr',ihail,improve,1,kte,dzw,rhoz,qrz,qvz,prez,tz,xland(i,j),vtr,precip,dtcfl)
       pptrain = pptrain + precip
    enddo

    do k = kts , kte
       qr(i,k,j) = qrz(k)
       vtr3d(i,k,j) = vtr(k)
    enddo
#else
    t_del_tv=0.
    del_tv=dtb
    notlast=.true.
    DO while (notlast)
!
      min_q=kte
      max_q=kts-1
!


      do k=kts,kte-1

         vtr(k)=0.

         if (qrz(k) .gt. crmin) then
            min_q=min0(min_q,k)
            max_q=max0(max_q,k)

           call vtr_mks(rhoz(k),qrz(k),tz(k),vtr(k))
!           if (.not. vtr(k) .gt. 0.0) cycle ! EMK NUWRF Bug fix

            if (k .eq. 1) then
               del_tv=dmin1(del_tv,0.9*(zz(k)-topo(i,j))/vtr(k))
            else
               del_tv=dmin1(del_tv,0.9*(zz(k)-zz(k-1))/vtr(k))
            endif
          endif
      enddo

      if (max_q .ge. min_q) then
!
!- Check if the summation of the small delta t >=  big delta t
!             (t_del_tv)          (del_tv)             (dtb)

         t_del_tv=t_del_tv+del_tv
!
         if ( t_del_tv .ge. dtb ) then
              notlast=.false.
              del_tv=dtb+del_tv-t_del_tv
         endif

! use small delta t to calculate the qrz flux
! termi is the qrz flux pass in the grid box through the upper boundary
! termo is the qrz flux pass out the grid box through the lower boundary
!
         fluxin=0.
         do k=max_q,min_q,-1
            fluxout=rhoz(k)*vtr(k)*qrz(k)
            flux=(fluxin-fluxout)/rhoz(k)/dzw(k)
!            tmpqrz=qrz(k)
            qrz(k)=qrz(k)+del_tv*flux
            qrz(k)=dmax1(0.,qrz(k))
            qr(i,k,j)=qrz(k)
            fluxin=fluxout
            rsed(k)=rsed(k)+fluxin
         enddo
         if (min_q .eq. 1) then
            pptrain=pptrain+fluxin*del_tv
         else
            qrz(min_q-1)=qrz(min_q-1)+del_tv*  &
                          fluxin/rhoz(min_q-1)/dzw(min_q-1)
            qr(i,min_q-1,j)=qrz(min_q-1)
         endif
!
      else
         notlast=.false.
      endif
    ENDDO
#endif

!
!-- snow
!
#ifdef SL_sedi
    ntimes = 1
    dtcfl = dtb/real(ntimes)

    do n = 1 , ntimes
       vts(:) = 0.
       precip = 0.
       call semi_sedi('qs',ihail,improve,1,kte,dzw,rhoz,qsz,qvz,prez,tz,xland(i,j),vts,precip,dtcfl)
       pptsnow = pptsnow + precip
    enddo

    do k = kts , kte
       qs(i,k,j) = qsz(k)
       vts3d(i,k,j) = vts(k)
    enddo
#else
    t_del_tv=0.
    del_tv=dtb
    notlast=.true.

    DO while (notlast)
!
      min_q=kte
      max_q=kts-1

!
      do k=kts,kte-1
         vts(k)=0.

         if (qsz(k) .gt. csmin) then
            min_q=min0(min_q,k)
            max_q=max0(max_q,k)

            call vts_mks(improve,rhoz(k),qsz(k),tz(k),vts(k))

            if (k .eq. 1) then
               del_tv=dmin1(del_tv,0.9*(zz(k)-topo(i,j))/vts(k))
            else
               del_tv=dmin1(del_tv,0.9*(zz(k)-zz(k-1))/vts(k))
            endif
         endif
      enddo

      if (max_q .ge. min_q) then
!
!
!- Check if the summation of the small delta t >=  big delta t
!             (t_del_tv)          (del_tv)             (dtb)

         t_del_tv=t_del_tv+del_tv

         if ( t_del_tv .ge. dtb ) then
              notlast=.false.
              del_tv=dtb+del_tv-t_del_tv
         endif

! use small delta t to calculate the qsz flux
! termi is the qsz flux pass in the grid box through the upper boundary
! termo is the qsz flux pass out the grid box through the lower boundary
!
         fluxin=0.
         do k=max_q,min_q,-1
            fluxout=rhoz(k)*vts(k)*qsz(k)
            flux=(fluxin-fluxout)/rhoz(k)/dzw(k)
            qsz(k)=qsz(k)+del_tv*flux
            qsz(k)=dmax1(0.,qsz(k))
            qs(i,k,j)=qsz(k)
            fluxin=fluxout
            ssed(k)=ssed(k)+fluxin
         enddo
         if (min_q .eq. 1) then
            pptsnow=pptsnow+fluxin*del_tv
         else
            qsz(min_q-1)=qsz(min_q-1)+del_tv*  &
                         fluxin/rhoz(min_q-1)/dzw(min_q-1)
            qs(i,min_q-1,j)=qsz(min_q-1)
         endif
!
      else
         notlast=.false.
      endif

    ENDDO
#endif

!
!   ice2=0 --- with hail/graupel 
!   ice2=1 --- without hail/graupel 
!
  if (ice2.eq.0) then 
!
!-- If IHAIL=1, use hail.
!-- If IHAIL=0, use graupel.
!
!    if (ihail .eq. 1) then
!       xnog = xnoh
!       rhograul = rhohail
!    endif

#ifdef SL_sedi
    ntimes = 1
    dtcfl = dtb/real(ntimes)

    do n = 1 , ntimes
       vtg(:) = 0.
       precip = 0.
       call semi_sedi('qg',ihail,improve,1,kte,dzw,rhoz,qgz,qvz,prez,tz,xland(i,j),vtg,precip,dtcfl)
       pptgraul = pptgraul + precip
    enddo

    do k = kts , kte
       qg(i,k,j) = qgz(k)
       vtg3d(i,k,j) = vtg(k)
    enddo
#else
    t_del_tv=0.
    del_tv=dtb
    notlast=.true.
!
    DO while (notlast)
!
      min_q=kte
      max_q=kts-1
!
      do k=kts,kte-1
          vtg(k)=0.

         if (qgz(k) .gt. cgmin) then
            min_q=min0(min_q,k)
            max_q=max0(max_q,k)

            call vtg_mks(ihail,improve,rhoz(k),qgz(k),tz(k),vtg(k))

            if (k .eq. 1) then
               del_tv=dmin1(del_tv,0.9*(zz(k)-topo(i,j))/vtg(k))
            else
               del_tv=dmin1(del_tv,0.9*(zz(k)-zz(k-1))/vtg(k))
            endif 
!
         endif !qgz
      enddo !k

      if (max_q .ge. min_q) then
!
!
!- Check if the summation of the small delta t >=  big delta t
!             (t_del_tv)          (del_tv)             (dtb)

         t_del_tv=t_del_tv+del_tv

         if ( t_del_tv .ge. dtb ) then
              notlast=.false.
              del_tv=dtb+del_tv-t_del_tv
         endif

! use small delta t to calculate the qgz flux
! termi is the qgz flux pass in the grid box through the upper boundary
! termo is the qgz flux pass out the grid box through the lower boundary
!
         fluxin=0.
         do k=max_q,min_q,-1
            fluxout=rhoz(k)*vtg(k)*qgz(k)
            flux=(fluxin-fluxout)/rhoz(k)/dzw(k)
            qgz(k)=qgz(k)+del_tv*flux
            qgz(k)=dmax1(0.,qgz(k))
            qg(i,k,j)=qgz(k)
            fluxin=fluxout
            gsed(k)=gsed(k)+fluxin
         enddo
         if (min_q .eq. 1) then
            pptgraul=pptgraul+fluxin*del_tv
         else
            qgz(min_q-1)=qgz(min_q-1)+del_tv*  &
                         fluxin/rhoz(min_q-1)/dzw(min_q-1)
            qg(i,min_q-1,j)=qgz(min_q-1)
         endif
!
      else
         notlast=.false.
      endif
!
    ENDDO
#endif
 ENDIF !ice2
!
!-- cloud ice  (03/21/02) follow Vaughan T.J. Phillips at GFDL
!

#ifdef SL_sedi
    ntimes = 1
    dtcfl = dtb/real(ntimes)

    do n = 1 , ntimes
       vti(:) = 0.
       precip = 0.
       call semi_sedi('qi',ihail,improve,0,kte,dzw,rhoz,qiz,qvz,prez,tz,xland(i,j),vti,precip,dtcfl)
       pptice = pptice + precip
    enddo

    do k = kts , kte
       qi(i,k,j) = qiz(k)
       vti3d(i,k,j) = vti(k)
    enddo
#else
    t_del_tv=0.
    del_tv=dtb
    notlast=.true.
!
    DO while (notlast)
!
      min_q=kte
      max_q=kts-1
!
      do k=kts,kte-1
          
         vti(k)=0.

         if (qiz(k) .gt. cimin) then
            min_q=min0(min_q,k)
            max_q=max0(max_q,k)

            call vti_mks(improve,rhoz(k),tz(k),qiz(k),qvz(k),prez(k),xland(i,j),vti(k))

            ! EMK:  Avoid division by zero
            if ((vti(k) .gt. 1.0e-20)) then
               if (k .eq. 1) then
                  del_tv=dmin1(del_tv,0.9*(zz(k)-topo(i,j))/vti(k))
               else
                  del_tv=dmin1(del_tv,0.9*(zz(k)-zz(k-1))/vti(k))
               endif
            endif
         endif
      enddo

      if (max_q .ge. min_q) then
!
!
!- Check if the summation of the small delta t >=  big delta t
!             (t_del_tv)          (del_tv)             (dtb)

         t_del_tv=t_del_tv+del_tv

         if ( t_del_tv .ge. dtb ) then
              notlast=.false.
              del_tv=dtb+del_tv-t_del_tv
         endif

! use small delta t to calculate the qiz flux
! termi is the qiz flux pass in the grid box through the upper boundary
! termo is the qiz flux pass out the grid box through the lower boundary
!

         fluxin=0.
         do k=max_q,min_q,-1
            fluxout=rhoz(k)*vti(k)*qiz(k)
            flux=(fluxin-fluxout)/rhoz(k)/dzw(k)
            qiz(k)=qiz(k)+del_tv*flux
            qiz(k)=dmax1(0.,qiz(k))
            qi(i,k,j)=qiz(k)
            fluxin=fluxout
            ised(k)=ised(k)+fluxin
         enddo
         if (min_q .eq. 1) then
            pptice=pptice+fluxin*del_tv
         else
            qiz(min_q-1)=qiz(min_q-1)+del_tv*  &
                         fluxin/rhoz(min_q-1)/dzw(min_q-1)
            qi(i,min_q-1,j)=qiz(min_q-1)
         endif
!
      else
         notlast=.false.
      endif
!
   ENDDO !notlast
#endif

#ifdef EXT_DIAG
   do k = kts, kte
            preci3d(i,k,j)=ised(k)
            precs3d(i,k,j)=ssed(k)
            precg3d(i,k,j)=gsed(k)
            precr3d(i,k,j)=rsed(k)
   end do
#endif

!   prnc(i,j)=prnc(i,j)+pptrain
!   psnowc(i,j)=psnowc(i,j)+pptsnow
!   pgrauc(i,j)=pgrauc(i,j)+pptgraul
!   picec(i,j)=picec(i,j)+pptice
!                     

!   write(6,*) 'i=',i,' j=',j,'   ', pptrain, pptsnow, pptgraul, pptice
!   call flush(6)

   icenc(i,j) = icenc(i,j) + pptice
!   snowncv(i,j) = pptsnow
   snownc(i,j) = snownc(i,j) + pptsnow
!   graupelncv(i,j) = pptgraul
   graupelnc(i,j) = graupelnc(i,j) + pptgraul 
!   RAINNCV(i,j) = pptrain + pptsnow + pptgraul + pptice
   RAINNC(i,j)  = RAINNC(i,j) + pptrain + pptsnow + pptgraul + pptice
   pptall = pptrain + pptsnow + pptgraul + pptice
   sr(i,j) = 0.
   if ( pptall.gt.0. ) sr(i,j) = (pptsnow + pptgraul + pptice) / pptall

  ENDDO i_loop
  ENDDO j_loop

!  if (itimestep.eq.6480) then
!     write(51,*) 'in the end of fallflux, itimestep=',itimestep
!     do j=jts,jte
!        do i=its,ite 
!           if (rainnc(i,j).gt.400.) then
!              write(50,*) 'i=',i,' j=',j,' rainnc=',rainnc
!           endif
!        enddo
!     enddo
!  endif

  RETURN
  END SUBROUTINE fall_flux

!-----------------------------------------------------------------------
!c Correction of negative values  
   SUBROUTINE negcor ( X, rho, dz8w,                         &
                      ims,ime, jms,jme, kms,kme,              & ! memory dims
                      itimestep, ics,                         &
                      its,ite, jts,jte, kts,kte               ) ! tile   dims
!-----------------------------------------------------------------------
  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),                 &
        INTENT(INOUT) ::                                     X   
  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),                 &
        INTENT(IN   ) ::                              rho, dz8w  
  integer, INTENT(IN   ) ::                           itimestep, ics 

!c Local variables
!  REAL, DIMENSION( kts:kte ) ::  Y1, Y2
  REAL   ::   A0, A1, A2

  A1=0.
  A2=0.
  do k=kts,kte
     do j=jts,jte
        do i=its,ite
        A1=A1+max(X(i,k,j), 0.)*rho(i,k,j)*dz8w(i,k,j)
        A2=A2+max(-X(i,k,j), 0.)*rho(i,k,j)*dz8w(i,k,j)
        enddo
     enddo
  enddo

!  A1=0.0
!  A2=0.0
!  do k=kts,kte
!     A1=A1+Y1(k)
!     A2=A2+Y2(k)
!  enddo

  A0=0.0

  if (A1.NE.0.0.and.A1.GT.A2) then 
     A0=(A1-A2)/A1

!  if (mod(itimestep,540).eq.0) then
!     if (ics.eq.1) then
!        write(61,*) 'kms=',kms,'  kme=',kme,'  kts=',kts,'  kte=',kte
!        write(61,*) 'jms=',jms,'  jme=',jme,'  jts=',jts,'  jte=',jte 
!        write(61,*) 'ims=',ims,'  ime=',ime,'  its=',its,'  ite=',ite 
!     endif 
!     if (ics.eq.1) then
!         write(61,*) 'qv timestep=',itimestep
!         write(61,*) '  A1=',A1,'   A2=',A2,'   A0=',A0
!     else if (ics.eq.2) then
!             write(61,*) 'ql timestep=',itimestep
!             write(61,*) '  A1=',A1,'   A2=',A2,'   A0=',A0
!     else if (ics.eq.3) then
!             write(61,*) 'qr timestep=',itimestep
!             write(61,*) '  A1=',A1,'   A2=',A2,'   A0=',A0
!     else if (ics.eq.4) then
!             write(61,*) 'qi timestep=',itimestep
!             write(61,*) '  A1=',A1,'   A2=',A2,'   A0=',A0
!     else if (ics.eq.5) then
!             write(61,*) 'qs timestep=',itimestep
!             write(61,*) '  A1=',A1,'   A2=',A2,'   A0=',A0
!     else if (ics.eq.6) then
!             write(61,*) 'qg timestep=',itimestep
!             write(61,*) '  A1=',A1,'   A2=',A2,'   A0=',A0
!     else
!             write(61,*) 'wrong cloud specieis number'
!     endif 
!  endif 

     do k=kts,kte
        do j=jts,jte
           do i=its,ite
!           X(i,k,j)=A0*AMAX1(X(i,k,j), 0.0)
           X(i,k,j)=A0*DMAX1(X(i,k,j), 0.0)
           enddo
        enddo
     enddo
  endif

  END SUBROUTINE negcor

  SUBROUTINE consat_s (ihail, itaobraun, improve)  

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                      c
!   Tao, W.-K., and J. Simpson, 1989: Modeling study of a tropical     c
!   squall-type convective line. J. Atmos. Sci., 46, 177-202.          c
!                                                                      c
!   Tao, W.-K., J. Simpson and M. McCumber, 1989: An ice-water         c
!   saturation adjustment. Mon. Wea. Rev., 117, 231-235.               c

!                                                                      c
!   Tao, W.-K., and J. Simpson, 1993: The Goddard Cumulus Ensemble     c
!   Model. Part I: Model description. Terrestrial, Atmospheric and     c
!   Oceanic Sciences, 4, 35-72.                                        c
!                                                                      c
!   Tao, W.-K., J. Simpson, D. Baker, S. Braun, M.-D. Chou, B.         c
!   Ferrier,D. Johnson, A. Khain, S. Lang,  B. Lynn, C.-L. Shie,       c
!   D. Starr, C.-H. Sui, Y. Wang and P. Wetzel, 2003: Microphysics,    c
!   radiation and surface processes in the Goddard Cumulus Ensemble    c
!   (GCE) model, A Special Issue on Non-hydrostatic Mesoscale          c
!   Modeling, Meteorology and Atmospheric Physics, 82, 97-137.         c
!                                                                      c
!   Lang, S., W.-K. Tao, R. Cifelli, W. Olson, J. Halverson, S.        c
!   Rutledge, and J. Simpson, 2007: Improving simulations of           c
!   convective system from TRMM LBA: Easterly and Westerly regimes.    c
!   J. Atmos. Sci., 64, 1141-1164.                                     c
!                                                                      c
!   Coded by Tao (1989-2003), modified by S. Lang (2006/07)            c
!                                                                      c
!   Implemented into WRF  by Roger Shi 2006/2007                       c
!   Additional modifications by Tao, Roger and Steve 2009              c
!   July 25 2010                                                       c
!   Tao November 12 2010                                               c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!        itaobraun=0   ! see Tao and Simpson (1993)
!        itaobraun=1   ! see Tao et al. (2003)

 integer :: ihail, itaobraun, improve
 real    :: cn0

!JJS 1/3/2008  vvvvv
!JJS   the following common blocks have been moved to the top of
!JJS   module_mp_gsfcgce.F
!
! real,   dimension (1:31) ::  a1, a2
! data a1/.7939e-7,.7841e-6,.3369e-5,.4336e-5,.5285e-5,.3728e-5,       &
!      .1852e-5,.2991e-6,.4248e-6,.7434e-6,.1812e-5,.4394e-5,.9145e-5, &
!      .1725e-4,.3348e-4,.1725e-4,.9175e-5,.4412e-5,.2252e-5,.9115e-6, &
!      .4876e-6,.3473e-6,.4758e-6,.6306e-6,.8573e-6,.7868e-6,.7192e-6, &
!         .6513e-6,.5956e-6,.5333e-6,.4834e-6/
! data a2/.4006,.4831,.5320,.5307,.5319,.5249,.4888,.3894,.4047, &
!         .4318,.4771,.5183,.5463,.5651,.5813,.5655,.5478,.5203,.4906, &
!         .4447,.4126,.3960,.4149,.4320,.4506,.4483,.4460,.4433,.4413, &
!         .4382,.4361/
!JJS 1/3/2008  ^^^^^

!23456789012345678901234567890123456789012345678901234567890123456789012
!     ******************************************************************
!JJS
      if (improve .eq. 3) ihail = 0
      if (ihail .eq. 1) improve = -20
      al = 2.5e10
      cp = 1.004e7
      rd1 = 1.e-3
      rd2 = 2.2
!JJS
      cpi=4.*atan(1.)
      cpi2=cpi*cpi
      grvt=980.
!
      c76=7.66
      c141=1.414435e7
      c149=1.496286e-5
      c172=17.26939
      c218=21.87456
      c358=35.86
      c409=4098.026
      c580=5807.695
      c610=6.1078e3
      c879=8.794142
!
      cd1=6.e-1
      cd2=4.*grvt/(3.*cd1)
      tca=2.43e3
      dwv=.226
      dva=1.718e-4
      amw=18.016
      ars=8.314e7
      scv=2.2904487
!
      t0=273.16
      t00=238.16
      alv=2.5e10
      alf=3.336e9
      als=2.8336e10
      avc=alv/cp
      afc=alf/cp
      asc=als/cp
      rw=4.615e6
      cw=4.187e7
      ci=2.093e7
!***   DEFINE THE COEFFICIENTS USED IN TERMINAL VELOCITY
!***   DEFINE THE DENSITY AND SIZE DISTRIBUTION OF PRECIPITATION
!**********   HAIL OR GRAUPEL PARAMETERS   **********

! tnw                  !  rain intercept (1/cm**4)
! roqr                 !  rain density (g/cm**3)
! tns                  !  snow intercept (1/cm**4)
! roqs                 !  snow density (g/cm**3)
! tng                  !  graupel intercept (1/cm**4)
! roqg                 !  graupel density (g/cm**3)


      if(ihail .eq. 1) then
         roqg=.9
         ag=sqrt(cd2*roqg)
         bg=.5
         tng=.002
      else
!         roqg=.4
!         ag=351.2
!         bg=.37
         roqg=.3
         ag=330.22    !for bulk density of 0.3
         bg=.36
         tng=.04
      endif

!**********         SNOW PARAMETERS        **********
!                             6/15/02
!      TNS=1.
!      TNS=.08             ! if ice913=1, tao's
      tns=.16              ! if ice913=0, tao's
      roqs=.1
!      AS=152.93
      as=78.63
!      BS=.25
      bs=.11
      if(improve.eq.3) then
        roqs=.05        !  snow density (g/cm**3)
        tns=0.1         !  snow intercept (1/cm**4)
        as=151.01
        bs=0.24
      endif

!**********         RAIN PARAMETERS        **********
      aw=2115.
      bw=.8
      roqr=1.  !not defined in Steve's
      tnw=.08  !not defined in Steve's
!*****************************************************************
      bgh=.5*bg
      bsh=.5*bs
      bwh=.5*bw
      bgq=.25*bg
      bsq=.25*bs
      bwq=.25*bw
!**********GAMMA FUNCTION CALCULATIONS*************
      ga3=2.
      ga4=6.
      ga5=24.
      ga6=120.
      ga7=720.
      ga8=5040.
      ga9=40320.
!
      ga3b  = gammagce(3.+bw)
      ga4b  = gammagce(4.+bw)
      ga6b  = gammagce(6.+bw)
      ga5bh = gammagce((5.+bw)/2.)
      ga3g  = gammagce(3.+bg)
      ga4g  = gammagce(4.+bg)
      ga5gh = gammagce((5.+bg)/2.)
      ga3d  = gammagce(3.+bs)
      ga4d  = gammagce(4.+bs)
      ga5dh = gammagce((5.+bs)/2.)
!        if(bg.eq.0.37) ga4g=9.730877
!        if(bg.eq.0.37) ga3g=2.8875
!        if(bg.eq.0.37) ga5gh=1.526425
!        if(bs.eq.0.57) ga3d=3.59304
!        if(bs.eq.0.57) ga4d=12.82715
!        if(bs.eq.0.57) ga5dh=1.655588
!        if(bs.eq.0.11) ga3d=2.218906
!        if(bs.eq.0.11) ga4d=6.900796
!        if(bs.eq.0.11) ga5dh=1.382792

      if(improve .eq. 3) then
        ga4g=11.63177
        ga3g=3.3233625          
        ga5gh=1.608355         
        if(bg.eq.0.37) ga4g=9.730877 
        if(bg.eq.0.37) ga3g=2.8875
        if(bg.eq.0.37) ga5gh=1.526425
        if(bg.eq.0.36) ga4g=9.599978
        if(bg.eq.0.36) ga3g=2.857136
        if(bg.eq.0.36) ga5gh=1.520402  
          ga3d=2.54925           
          ga4d=8.285063         
          ga5dh=1.456943
          if(bs.eq.0.57) ga3d=3.59304
          if(bs.eq.0.57) ga4d=12.82715
          if(bs.eq.0.57) ga5dh=1.655588
          if(bs.eq.0.24) ga3d=2.523508
          if(bs.eq.0.24) ga4d=8.176166
          if(bs.eq.0.24) ga5dh=1.451396
          if(bs.eq.0.11) ga3d=2.218906
          if(bs.eq.0.11) ga4d=6.900796
          if(bs.eq.0.11) ga5dh=1.382792 
      endif                       
      ga6d=144.93124
        if(bs.eq.0.24) ga6d=181.654791

!
!CCCCC        LIN ET AL., 1983 OR LORD ET AL., 1984   CCCCCCCCCCCCCCCCC
      ac1=aw
      ac2=ag            ! Steve only defines ac1 and ac2
      ac3=as            ! need to talk about these 3 parameters.
!       ac1=as          ! used in Steve's codes
!       ac2=ag          ! used in Steve's codes

      bc1=bw
      cc1=as
      dc1=bs

      zrc=(cpi*roqr*tnw)**0.25
      zsc=(cpi*roqs*tns)**0.25
      zgc=(cpi*roqg*tng)**0.25
      vrc=aw*ga4b/(6.*zrc**bw)

       vrc0=-26.7
       vrc1=20600./zrc
       vrc2=-204500./(zrc*zrc)
       vrc3=906000./(zrc*zrc*zrc)

      vsc=as*ga4d/(6.*zsc**bs)
      vgc=ag*ga4g/(6.*zgc**bg)

!     ****************************
!     RN1=1.E-3
      rn1=9.4e-15                    ! 6/15/02 tao's
      bnd1=6.e-4
      rn2=1.e-3
!     BND2=1.25E-3
!     BND2=1.5E-3                    ! if ice913=1 6/15/02 tao's
      bnd2=2.0e-3                    ! if ice913=0 6/15/02 tao's
!      rn3=.25*cpi*tns*cc1*ga3d      ! cc1 = as
      rn3=.25*cpi*tns*as*ga3d
!
       esi=.1
      if (improve .eq. 3) esi=0.25
      rn3=.25*cpi*tns*as*esi*ga3d

      esw=1.   ! Steve uses esc and ac1
      esc=1.
        if (improve.eq.3) esc=0.45
      rn4=.25*cpi*esw*tns*as*ga3d   
        if (improve .eq. 3) rn4=.25*cpi*esc*tns*as*ga3d
!     ERI=1.
      eri=.1                        ! 6/17/02 tao's ice913=0 (not 1)
      rn5=.25*cpi*eri*tnw*ac1*ga3b
         if (improve .eq. 3) rn5=.25*cpi*eri*tnw
      rn50=-.267e2*ga3
      rn51=5.15e3*ga4
      rn52=-1.0225e4*ga5
      rn53=7.55e3*ga6

!     AMI=1./(24.*4.19E-10)
      ami=1./(24.*6.e-9)            ! 6/15/02 tao's
      rn6=cpi2*eri*tnw*ac1*roqr*ga6b*ami
        if (improve .eq. 3) rn6=cpi2*eri*tnw*roqr*ami
      rn60=-.267e2*ga6
      rn61=5.15e3*ga7
      rn62=-1.0225e4*ga8
      rn63=7.55e3*ga9

!     ESR=1.                       ! also if ice913=1 for tao's
      esr=.5                       ! 6/15/02 for ice913=0 tao's
      if (improve .eq. 3) esr=1.
      rn7=cpi2*esr*tnw*tns*roqs
      esr=1.
      rn8=cpi2*esr*tnw*tns*roqr
!
      egs=.1

      rn9=cpi2*tns*tng*roqs
        if (improve.eq.-1) egs=.01
        if (improve .eq. 3) rn9=cpi2*egs*tns*tng*roqs
!
      rn10=2.*cpi*tns
        If (improve .eq. 3)  rn10=4.*tns
      rn101=.31*ga5dh*sqrt(cc1)
      rn10a=als*als/rw
      If (improve .eq. 3) then
         rn101=.65
         rn102=.44*sqrt(as/dva)*ga5dh
         rn10a=alv*als*amw/(tca*ars)
      endif
      rn10b=alv/tca
      rn10c=ars/(dwv*amw)
!
      rn11=2.*cpi*tns/alf
        If (improve .eq. 3)  rn11=2.*cpi*tns*tca/alf
      rn11a=cw/alf

!     AMI50=1.51e-7
      ami50=3.84e-6               ! 6/15/02 tao's
!     AMI40=2.41e-8
      ami40=3.08e-8               ! 6/15/02 tao's
      ami50=4.8e-7*(100./50.)**3         ! Roger: actually = 3.84e-6 
        if (improve.eq.3) ami50=4.8e-7
!       ami40=2.46e-7
!       ami40=2.46e-7*.5**3              ! Roger: actually = 3.08e-8
        if (improve.eq.3) ami40=2.46e-7
!
      eiw=1.
!     UI50=20.
      ui50=100. ! 6/15/02 tao's
      ri50=2.*5.e-3
        if (improve.eq.3) ri50=5.e-3
!
      cmn=1.05e-15
      rn12=cpi*eiw*ui50*ri50**2
      do k=1,31
         y1=1.-aa2(k)
         rn13(k)=aa1(k)*y1/(ami50**y1-ami40**y1)
         rn12a(k)=rn13(k)/ami50
         rn12b(k)=aa1(k)*ami50**aa2(k)
         rn25a(k)=aa1(k)*cmn**aa2(k)
           BergCon1(k)=6.*aa1(k)*ami50**(aa2(k)-1.)
         BergCon2(k)=-2.*aa1(k)*ami50**aa2(k)*1.2

         BergCon3(k)=6.*aa2(k)/((aa2(k)+1.)*(aa2(k)+2.))        &
                   *aa1(k)*ami50**(aa2(k)-1.)
         BergCon4(k)=2.*(1.-aa2(k))/((aa2(k)+1.)*(aa2(k)+2.))   &
                   *aa1(k)*ami50**aa2(k)*1.2
      enddo
!
      egw=1.    !Roger: ewc in Steve's code

        if(improve.ge.3) egw=0.65

      rn14=.25*cpi*egw*tng*ga3g*ag
!       egc=1.
!       rn14=.25*cpi*egc*ac2*tng*ga3g      ! ac2 = ag in Steve's code
      egi=.1
        if(improve.eq.-1) egi=.001
      rn15=.25*cpi*egi*tng*ga3g*ag
!     rn15=.25*cpi*egi*tng*ac2*ga3g        ! ac2 = ag in Steve's code

      egi=1.
      rn15a=.25*cpi*egi*tng*ga3g*ag
!     rn15a=.25*cpi*egi*tng*ac2*ga3g       ! ac2 = ag in Steve's code 

      egr=1.
      rn16=cpi2*egr*tng*tnw*roqr
      rn17=2.*cpi*tng
      rn17a=.31*ga5gh*sqrt(ag)
      rn17b=cw-ci
      rn17c=cw
      if (improve .eq. 3) then
         rn171=2.*cpi*tng*alv*dwv
         rn172=2.*cpi*tng*tca
         rn17a=.31*ga5gh*sqrt(ag/dva)
      endif
!
      apri=.66
      bpri=1.e-4
      bpri=0.5*bpri                        ! 6/17/02 tao's
      rn18=20.*cpi2*bpri*tnw*roqr
      rn18a=apri
      rn19=2.*cpi*tng/alf
!      if (improve .eq. 3) rn19=2.*cpi*tng*tca/alf
       rn191=.78                     !Are this same as rnn191 (listed below)?
       rn192=.31*ga5gh*sqrt(ag/dva)  !Are this same as rnn192 (listed below)?
      rn19a=.31*ga5gh*sqrt(ag)
      rn19b=cw/alf 
      
!      rnn191=.78
!      rnn192=.31*ga5gh*sqrt(ag/dva)
       if (improve .eq. 3) then
          rn19=2.*cpi*tng*tca/alf
          rn19a=cw/alf
       endif
!
      rn20=2.*cpi*tng
      rn20a=als*als/rw
      rn20b=.31*ga5gh*sqrt(ag)
      rn30a=alv*alv*amw/(tca*ars) ! EMK per Roger's 20110816 code
      rn30=2.*cpi*tng             ! EMK per Roger's 20110816 code
      if (improve .eq. 3) then
         rn20a=als*als*amw/(tca*ars)
         rn20b=als/tca
!         rn30=2.*cpi*tng
!         rn30a=alv*alv*amw/(tca*ars)
      endif
!
      bnd3=2.e-3
      rn21=1.e3*1.569e-12/0.15
         if (improve .eq. 3)  rn21=1.e-3
      bnd21=1.5e-3
!
      erw=1. ! erc in Steve's code
      rn22=.25*cpi*erw*ac1*tnw*ga3b
      if (improve .eq. 3) then
         erc=1.
         rn22=.25*cpi*erc*tnw
      endif
!
      rn23=2.*cpi*tnw
      rn23a=.31*ga5bh*sqrt(ac1)
      rn23b=alv*alv/rw
      if (improve .eq. 3) then
         rn231=.78
         rn232=.31*ga3*sqrt(3.e3/dva)
      endif
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc
!cc        "c0" in routine      "consat" (2d), "consatrh" (3d)
!cc        if ( itaobraun.eq.1 ) --> betah=0.5*beta=-.46*0.5=-0.23;   cn0=1.e-6
!cc        if ( itaobraun.eq.0 ) --> betah=0.5*beta=-.6*0.5=-0.30;    cn0=1.e-8

       if (improve .eq. 3) itaobraun=1

       if (itaobraun .eq. 0) then
         cn0=1.e-8
         beta=-.6
       else if (itaobraun .eq. 1) then
         cn0=1.e-6
         beta=-.46
       endif
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      rn25=cn0
!      rn30a=alv*als*amw/(tca*ars)
      rn30b=alv/tca
      rn30c=ars/(dwv*amw)
      rn31=1.e-17

      rn32=4.*51.545e-4
!!
      rn33=4.*tns
       rn331=.65
       rn332=.44*sqrt(as/dva)*ga5dh 
!        rn332=.44*sqrt(ac1/dva)*ga5dh
!
      if (improve .eq. 3) then
!         rn30=2.*cpi*tng
         amc=1./(24.*4.e-9)
         rn34=cpi2*esc*amc*as*roqs*tns*ga6d  
!         print *,'rn34, cpi2, esc, amc, as, roqs, tns, ga6d = ', rn34, cpi2, esc, amc, as, roqs, tns, ga6d 
         rn35=alv*alv/(cp*rw)
      endif

!----- for teminal velocity and fall_flux --------
      constb = 0.8
!      constd = 0.25
      constd = 0.11
      consta = 2115.0*0.01**(1-constb)   ! =841.9967
!      constc=152.93*0.01**(1-constd)
      constc=78.63*0.01**(1-constd)      ! =11.72

      o6 = 1./6.
      cdrag = 0.6
      abar = 19.3
      bbar = 0.37
      rhoe_s = 1.29           !surface air density (km/m^3)

      xnor = tnw*1.0e8        !intercept parameter of rain (m^-4)
      xnos = tns*1.0e8        !intercept parameter of snow (m^-4)
      xnog = tng*1.0e8        !intercept parameter of graupel (m^-4)
      rhowater = roqr*1000.   !density of rain (kg/m^3)
      rhosnow = roqs*1000.    !density of snow (kg/m^3)
      rhograul = roqg*1000.   !density of graupel (kg/m^3)
      if (ihail.eq.1) then
         xnoh = xnog          !intercept parameter of hail (m^-4)
         rhohail = rhograul   !density of hail (kg/m^3)
      endif

      draimax=0.0500       !maximum rain diameter (cm)
      draimax=draimax**4.*roqr*cpi

      ! critical q of hydrometeor characteristics (eg. fall speed, radius)
      cwmin = 1.e-12
      cimin = 1.e-12
      crmin = 1.e-10
      csmin = 1.e-10
      cgmin = 1.e-10

    return
  END SUBROUTINE consat_s 

!JJS
!JJS      REAL FUNCTION GAMMA(X)
!JJS        Y=GAMMLN(X)
!JJS        GAMMA=EXP(Y)
!JJS      RETURN
!JJS      END
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!JJS      real function GAMMLN (xx)
  real function gammagce (xx)
!**********************************************************************
  real*8 cof(6),stp,half,one,fpf,x,tmp,ser
  data cof,stp /  76.18009173,-86.50532033,24.01409822, &
     -1.231739516,.120858003e-2,-.536382e-5, 2.50662827465 /
  data half,one,fpf / .5, 1., 5.5 /
!
      x=xx-one
      tmp=x+fpf
      tmp=(x+half)*log(tmp)-tmp
      ser=one
      do  j=1,6
         x=x+one
        ser=ser+cof(j)/x
      enddo !j
      gammln=tmp+log(stp*ser)
!JJS
      gammagce=exp(gammln)
!JJS
      return
 END FUNCTION gammagce

! compute base snow/graupel intercept scaling factor 
  subroutine sgmap(isg,qsg,r00,tairc,ftnsg)
  implicit none

!      common/size/ tnw,tns,tng,roqs,roqg,roqr  !defined in the beginning of the module
  integer, intent(in)  :: isg
  real,    intent(in)  :: qsg, r00, tairc
  real,    intent(out) :: ftnsg
  
! LOCAL variables

!  for snomap 
  real ::  xs,sno11,sno00,dsno11,dsno00,sexp11,sexp00,stt,   &
                stexp,sbase,tslopes,dsnomin,slim
!  for grpmap
  real ::  xg,grp11,grp00,dgrp11,dgrp00,gexp11,gexp00,gtt,   &
                gtexp,gbase,tslopeg,dgrpmin, glim

  real :: taird, qsg1, xx, fexp, cpi, cmin
  real :: ftnsT, sno1, dsno1, sexp1, ftnsQ, tnsmax
  real :: ftngT, grp1, dgrp1, gexp1, ftngQ, tngmax

  CPI=4.*ATAN(1.)
  cmin=1.e-20

  tslopes=0.03579323   ! increase tns by 3.5 from 0 to -35C
  tslopeg=0.03138892   ! increase tng by  3 from  0 to -35C

  dsnomin=0.0110       !minimum snow diameter (cm)
  dgrpmin=0.0145       !minimum graupel diameter (cm)

  dsnomin=dsnomin**4.*roqs*cpi
  dgrpmin=dgrpmin**4.*roqg*cpi

  xs=0.97
  sno11=0.65                !cold aloft
  sno00=0.30                !near melting level
  dsno11=1.25               !Tao 02/23/2012 used to be 1.15
  dsno00=0.60               !Tao 02/23/2012 used to be 0.70
  sexp11=1.2
  sexp00=0.7
  stt=-20.
  stexp=0.25
  slim=0.95
  sbase=0.00110

  xg=0.98
  grp11=0.75
  grp00=0.45
  dgrp11=3.95
  dgrp00=0.25
  gexp11=0.35
  gexp00=0.6
  gtt=-20.
  gtexp=0.30
  glim=0.90
  gbase=0.0058

  ftnsg=1.

  if (qsg .gt. cmin) then
      
     taird=min(0.,max(-35.,tairc)+0.0)
     qsg1=qsg*r00*1.e6

     if (isg.eq.1) then                          !snow
        ftnsT=exp(-1.*tslopes*taird)
        sno1=sno11
        dsno1=dsno11
        sexp1=sexp11

        if (taird.gt.stt) then
           sno1=sno00-(sno00-sno11)*(taird/stt)**stexp
           dsno1=dsno00-(dsno00-dsno11)*(taird/stt)**stexp
           sexp1=sexp00-(sexp00-sexp11)*(taird/stt)**stexp
        endif !taird

        xx=xs-xs*min(slim,max(0.,(qsg1-sno1)/dsno1)**sexp1)
        ftnsT=ftnsT**xx
        fexp=xx
        ftnsQ=1.0
        ftnsQ=(qsg1/sbase)**fexp
        ftnsg=ftnsT*ftnsQ
        tnsmax=r00*qsg/dsnomin

        if (ftnsg*tns.gt.tnsmax) ftnsg=tnsmax/tns

     else if (isg.eq.2) then                      !graupel
     
             ftngT=exp(-1.*tslopeg*taird)
             grp1=grp11
             dgrp1=dgrp11
             gexp1=gexp11

             if (taird.gt.gtt) then
                grp1=grp00-(grp00-grp11)*(taird/gtt)**gtexp
                dgrp1=dgrp00-(dgrp00-dgrp11)*(taird/gtt)**gtexp
                gexp1=gexp00-(gexp00-gexp11)*(taird/gtt)**gtexp
             endif !taird

             xx=xg-xg*min(glim,max(0.0,(qsg1-grp1)/dgrp1)**gexp1)
             ftngT=ftngT**xx
             fexp=xx
             ftngQ=1.0
             ftngQ=(qsg1/gbase)**fexp
             ftnsg=ftngT*ftngQ
             tngmax=r00*qsg/dgrpmin
             if (ftnsg*tng.gt.tngmax) ftnsg=tngmax/tng

     endif !isg

  endif !qsg

  return
  end subroutine sgmap

!     compute fall speed of cloud rain and ice
  subroutine vqrqi(isg,improve,r00,fv,qri,ww1)
  implicit none

!      common/size/ tnw,tns,tng,roqs,roqg,roqr    ! defined in the beginning of the module
!      common/rterv/ zrc,zgc,zsc,vr0,vr1,vr2,vr3,vgc,vsc   ! defined in the beginning of the module

  integer, intent(in) :: isg, improve
  real, intent(in) :: r00, fv
  real, intent(inout) :: qri, ww1 

! LOCAL variables
  integer :: ic 
  real  :: y1,vr,vs,vg
  real  :: cmin, cpi, const_vt, const_d, const_m 
  real  :: bb1, bb2
  real, dimension(7) ::  aice, vice
  data aice/1.e-6, 1.e-5, 1.e-4, 1.e-3, 0.01, 0.1, 1./
  data vice/5,15,30,35,40,45,50/

  cmin=1.e-40
  CPI=4.*ATAN(1.)
  if (improve.gt.2) then 
     const_vt=1.49e4      
     const_d=11.9
     const_m=1./5.38e7
  endif 

  y1=r00*qri
  ww1=0.

  if (y1 .gt. cmin) then

    if (isg.eq.1) then                             !  rain

       vs=sqrt( y1 )
       vg=sqrt( vs )
!       vr=vr0+vr1*vg+vr2*vs+vr3*vg*vs
       vr=vrc0+vrc1*vg+vrc2*vs+vrc3*vg*vs
       ww1=max(fv*vr, 0.e0)

    else if (isg.eq.2) then                         ! cloud ice

            y1=1.e6*r00*qri                            ! to g/m**3

            if (y1 .gt. 1.e-6) then
               if (improve.eq.3) then                      ! from Hong et al. (2004)
                  y1=y1*1.e-3
                  bb1=const_m*y1**0.25
                  bb2=const_d*bb1**0.5
                  ww1=max(const_vt*bb2**1.31, 0.0)
                  ww1=ww1*100. !cm/s
                  if (ww1 .gt. 50.) ww1=50.               ! SLang
               else                                      ! from Starr & Cox (1985)
                  if (y1 .gt. aice(7)) then
                     ww1=vice(7)
                  else
                     do ic=1,6
                        if (y1.gt.aice(ic) .and. y1.le.aice(ic+1)) then
                           ww1=vice(ic)+(vice(ic+1)-vice(ic))*    &
                             (y1-aice(ic))/(aice(ic+1)-aice(ic))
                           if (ww1 .le. 0.0) ww1=0.
                        endif
                     enddo
                  endif  !y1
               endif  !imrpove
            endif  !y1
    endif  !isg
  endif !y1

  return
  end subroutine vqrqi

  SUBROUTINE saticel_s (dt, ihail, itaobraun, ice2, istatmin,          &
                       new_ice_sat, id, improve, xlat, sdec,           &
                       ptwrf, qvwrf, qlwrf, qrwrf,                     &
                       qiwrf, qswrf, qgwrf,                            &
                       rho_mks, pi_mks, p0_mks, w_mks,                 &
                       itimestep, xland,                               &
                       refc, refr, refi, refs, refg,                   & ! cloud effective radius
                       ims,ime, jms,jme, kms,kme,                      &
                       its,ite, jts,jte, kts,kte                       &
#ifdef EXT_DIAG
                       ,refl_10cm, diagflag, do_radar_ref,             & ! GT added for reflectivity calcs
                       physc, physe, physd, physs, physm, physf,       &   
                       acphysc, acphyse, acphysd, acphyss, acphysm, acphysf &   
#endif
#if ( WRF_CHEM == 1)
                       ,aero, icn_diag, nc_diag, gid,                  &
                       chem_opt,                                       &
                       gsfcgce_gocart_coupling                         &
#endif
                       )
!-----------------------------------------------------------------------
!  USE module_dm
  IMPLICIT NONE
!-----------------------------------------------------------------------
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                         c
!   History:                                                              c
!                                                                         c
!   Coded by Tao (1989-2003), modified by S. Lang (2006/07)               c
!                                                                         c
!   Implemented into WRF  by Jainn Shi 2006/2007                          c
!   Improved by S. Lang (2008-2009)                                       c
!   Implemented by Tao, Jainn Shi and tested by Jainn Shi                 c
!   Modified by Tao, Jul. 2010                                            c
!   Modified by Tao, Aug. 2010                                            c
!   Added aerosol coupling by Jainn Shi, Jun. 2012
!   Added cloud droplet eff. radius code by Jainn Shi, Mar. 2014          c
!                                                                         c
!   References:                                                           c
!                                                                         c
!   Tao, W.-K., and J. Simpson, 1989: Modeling study of a tropical        c
!   squall-type convective line. J. Atmos. Sci., 46, 177-202.             c
!                                                                         c
!   Tao, W.-K., J. Simpson and M. McCumber, 1989: An ice-water            c
!   saturation adjustment. Mon. Wea. Rev., 117, 231-235.                  c
!                                                                         c
!                                                                         c
!   Tao, W.-K., and J. Simpson, 1993: The Goddard Cumulus Ensemble        c
!   Model. Part I: Model description. Terrestrial, Atmospheric and        c
!   Oceanic Sciences, 4, 35-72.                                           c
!                                                                         c
!   Tao, W.-K., J. Simpson, D. Baker, S. Braun, M.-D. Chou, B.            c
!   Ferrier,D. Johnson, A. Khain, S. Lang,  B. Lynn, C.-L. Shie,          c
!   D. Starr, C.-H. Sui, Y. Wang and P. Wetzel, 2003: Microphysics,       c
!   radiation and surface processes in the Goddard Cumulus Ensemble       c
!   (GCE) model, A Special Issue on Non-hydrostatic Mesoscale             c
!   Modeling, Meteorology and Atmospheric Physics, 82, 97-137.            c
!                                                                         c
!   Lang, S., W.-K. Tao, R. Cifelli, W. Olson, J. Halverson, S.           c
!   Rutledge, and J. Simpson, 2007: Improving simulations of              c
!   convective system from TRMM LBA: Easterly and Westerly regimes.       c
!   J. Atmos. Sci., 64, 1141-1164.                                        c
!                                                                         c
!   Tao, W.-K., J. J. Shi,  S. Lang, C. Peters-Lidard, A. Hou, S.         c
!   Braun, and J. Simpson, 2007: New, improved bulk-microphysical         c
!   schemes for studying precipitation processes in WRF. Part I:          c
!   Comparisons with other schemes.                                       c
!                                                                         c
!   Lang, S., W.-K. Tao, X. Zeng, and Y. Li, 2011: Reducing the Biases in c
!   Simulated Radar Reflectivities from a Bulk Microphysics Scheme:       c
!   Tropical Convective Systems. J. Atmos. Sci., 68, 2306-2320.           c
!                                                                         c
!   Shi, J. J., T. Matsui, W.-K. Tao, C. Peters-Lidard, M. Chin, Q. Tan,  c
!   K. Pickering, N. Guy, S. Lang, and E. Kemp., 2014: Implementation of  c
!   an Aerosol-Cloud Microphysics-Radiation Coupling into the NASA        c
!   Unified WRF:  Simulation Results for the 6-7 August 2006 AMMA Special c
!   Observing Period. Quart. J. Roy. Meteor. Soc., in press.              c
!                                                                         c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!      COMPUTE ICE PHASE MICROPHYSICS AND SATURATION PROCESSES
!
  integer,    parameter ::  nt=2880, nt2=2*nt 

!cc   using scott braun's way for pint, pidep computations
  integer  ::   itaobraun,ice2,ihail,new_ice_sat,id,istatmin
  integer  ::   itimestep, improve
  real     ::   tairccri, cn0, dt
!cc

!JJS      common/bxyz/ n,isec,nran,kt1,kt2
!JJS      common/option/ lipps,ijkadv,istatmin,iwater,itoga,imlifting,lin,
!JJS     1   irf,iadvh,irfg,ismg,id

!JJS      common/timestat/ ndt_stat,idq
!JJS      common/iice/ new_ice_sat
!JJS      common/bt/ dt,d2t,rijl2,dts,f5,rd1,rd2,bound,al,cp,ra,ck,ce,eps,
!JJS     1    psfc,fcor,sec,aminut,rdt

!JJS   the following common blocks have been moved to the top of 
!JJS   module_mp_gsfcgce_driver_instat.F

!      common/bt/ rd1,rd2,al,cp
!
!
!      common/bterv/ zrc,zgc,zsc,vrc,vgc,vsc
!      common/size/ tnw,tns,tng,roqs,roqg,roqr
!      common/cont/ c38,c358,c610,c149,c879,c172,c409,c76,c218,c580,c141
!        common/b3cs/ ag,bg,as,bs,aw,bw,bgh,bgq,bsh,bsq,bwh,bwq
!      common/bsnw/ alv,alf,als,t0,t00,avc,afc,asc,rn1,bnd1,rn2,bnd2,     &
!         rn3,rn4,rn5,rn6,rn7,rn8,rn9,rn10,rn101,rn10a,rn11,rn11a, &
!         rn12,rn12a(31),rn12b(31),rn13(31),rn14,rn15,rn15a,rn16,rn17,    &
!         rn17a,rn17b,rn17c,rn18,rn18a,rn19,rn19a,rn19b,rn20,rn20a,rn20b, &
!         bnd3,rn21,rn22,rn23,rn23a,rn23b,rn25,rn25a(31),rn30a,rn30b,     &
!         rn30c,rn31,beta,rn32
!      common/rsnw1/ rn10b,rn10c,rnn191,rnn192,rn30,rnn30a,rn33,rn331,    &
!                    rn332
!JJS

!JJS 20090623 vvvvv
!  integer ids,ide,jds,jde,kds,kde                                  
  integer ims,ime,jms,jme,kms,kme                                  
  integer its,ite,jts,jte,kts,kte                                  
  integer i,j,k, kp     

  real   ::   a0 ,a1 ,a2 ,afcp ,alvr ,ami100 ,ami40 ,ami50       &
             ,ascp ,avcp ,betah, bg3 ,bgh5 ,bs3 ,bs6 ,bsh5       &
             ,bw3 ,bw6 ,bwh5 ,cmin ,cmin1 ,cmin2 ,cp409          &
             ,cp580 ,cs580 ,cv409 ,d2t ,del ,dwvp ,ee1 ,ee2      &
             ,f00 ,f2 ,f3 ,ft ,fv0 ,fvs ,pi0 ,pir ,pr0 ,qb0      &
             ,r00 ,r0s ,r101f ,r10ar ,r10t ,r11at ,r11rt         &
             ,r12r ,r14f ,r14r ,r15af ,r15ar ,r15f ,r15r         &
             ,r16r ,r17aq ,r17as ,r17r ,r18r ,r19aq ,r19as       &
             ,r19bt ,r19rt ,r20bq ,r20bs ,r20t ,r22f ,r23af      &
             ,r23br ,r23t ,r25a ,r25rt ,r2ice ,r31r ,r32rt       &
             ,r331r ,r332rf ,r34f                                &
             ,r3f ,r4f ,r5f ,r6f ,r7r ,r8r ,r9r ,r_nci ,rft      &
             ,rijl2 ,rp0 ,rr0 ,rrq ,rrs ,rt0 ,scc ,sccc          &
             ,sddd ,see ,seee ,sfff ,smmm ,ssss ,tb0 ,temp       &
             ,ucog ,ucor ,ucos ,uwet ,vgcf ,vgcr ,vrcf           &
             ,vscf ,zgr ,zrr ,zsr ,rdt

  real :: a_1, a_2, a_3
  real :: a_11, a_22, a_33, a_44
  real :: zdry, zwet
!JJS 20090623 ^^^^^

!JJS 20090708 vvvvv
  real :: fact_fit
!JJS 20090708 ^^^^^

  real, dimension (its:ite,jts:jte,kts:kte) ::  fv
  real, dimension (its:ite,jts:jte,kts:kte) ::  dpt, dqv
  real, dimension (its:ite,jts:jte,kts:kte) ::  qcl, qrn,             &
                                                qci, qcs, qcg
!JJS

  real, dimension (ims:ime, kms:kme, jms:jme) ::  ptwrf, qvwrf 
  real, dimension (ims:ime, kms:kme, jms:jme) ::  qlwrf, qrwrf,        &
                                                  qiwrf, qswrf, qgwrf

!JJS in MKS
  real, dimension (ims:ime, kms:kme, jms:jme) ::  rho_mks
  real, dimension (ims:ime, kms:kme, jms:jme) ::  pi_mks
  real, dimension (ims:ime, kms:kme, jms:jme) ::  p0_mks
  real, dimension (ims:ime, kms:kme, jms:jme) ::  w_mks
!JJS

!JJS 20140225
  real, dimension (ims:ime, kms:kme, jms:jme) , INTENT(OUT   )         &
                                              ::  refc, refr,         &
                                                  refi, refs,         &
                                                  refg 

!  real, dimension (its:ite,jts:jte,kts:kte) ::  ww1
!  real, dimension (its:ite,jts:jte,kts:kte) ::  rsw
!  real, dimension (its:ite,jts:jte,kts:kte) ::  rlw

!JJS      COMMON /BADV/
  real, dimension (its:ite,jts:jte) ::                                  &              
           vg,      zg,       &
           ps,      pg,       &
          prn,     psn,       &
        pwacs,   wgacr,       &
        pidep,    pint,       &
          qsi,     ssi,       &
          esi,     esw,       &
          qsw,      pr,       &
          ssw,   pihom,       &
         pidw,   pimlt,       &
        psaut,   qracs,       &
        psaci,   psacw,       &
        qsacw,   praci,       &
        pmlts,   pmltg,       &
        asss,       y1,    y2
!JJS       Y2(its:ite,jts:jte),    DDE(NB)

!JJS      COMMON/BSAT/
  real, dimension (its:ite,jts:jte) ::        &
        praut,   pracw,       &
         psfw,    psfi,       &
        dgacs,   dgacw,       &
        dgaci,   dgacr,       &
        pgacs,   wgacs,       &
        qgacw,   wgaci,       &
        qgacr,   pgwet,       &
        pgaut,   pracs,       &
        psacr,   qsacr,       &
         pgfr,   psmlt,       &
        pgmlt,   psdep,       &
        pgdep,   piacr,       &
           y5,     scv,       &
          tca,     dwv,       &
          egs,      y3,       &
           y4,     ddb

!JJS      COMMON/BSAT1/
  real, dimension (its:ite,jts:jte) ::        &
           pt,      qv,       &
           qc,      qr,       &
           qi,      qs,       &
           qg,    tair,       &
        tairc,   rtair,       &
          dep,      dd,       &
          dd1,     qvs,       &
           dm,      rq,       &
        rsub1,     col,       &
          cnd,     ern,       &
         dlt1,    dlt2,       &
         dlt3,    dlt4,       &
           zr,      vr,       &
           zs,      vs,       &
        pssub,   pgsub,       &
          dda

!JJS      COMMON/B5/
  real, dimension (its:ite,jts:jte,kts:kte) ::  rho
  real, dimension (kts:kte) ::                 & 
           tb,      qb,    rho1,               &
           ta,      qa,     ta1,     qa1,      &
         coef,      z1,      z2,      z3,      &
           am,     am1,      ub,      vb,      &
           wb,     ub1,     vb1,    rrho,      &
        rrho1,     wbx

!JJS      COMMON/B6/
  real, dimension (its:ite,jts:jte,kts:kte) ::  p0, pi, f0, ww1
  real, dimension (kts:kte) ::    & 
           fd,      fe,           &
           st,      sv,           &
           sq,      sc,           &
           se,     sqa

!JJS      COMMON/BRH1/
  real, dimension (kts:kte) ::             & 
         srro,    qrro,    sqc,    sqr,    &
          sqi,     sqs,    sqg,   stqc,    &
         stqr,    stqi,   stqs,   stqg
  real, dimension (nt) ::                  & 
          tqc,     tqr,    tqi,    tqs,    tqg

!JJS     common/bls/ y0(nx,ny),ts0new(nx,ny),qss0new(nx,ny)

!JJS      COMMON/BLS/
  real, dimension (ims:ime,jms:jme) ::     &
           y0,     ts0,   qss0

!JJS      COMMON/BI/ IT(its:ite,jts:jte), ICS(its:ite,jts:jte,4)
  integer, dimension (its:ite,jts:jte) ::        it  
  integer, dimension (its:ite,jts:jte, 4) ::    ics 

  integer :: i24h
  integer :: iwarm
  real :: r2is, r2ig
  

!JJS      COMMON/MICRO/
!NUWRF BEGIN
#ifdef EXT_DIAG
  real, dimension (ims:ime,kms:kme,jms:jme)  ::     &
          physc,   physe,   physd,                  &
          physs,   physm,   physf,                  &
          acphysc,   acphyse,   acphysd,            &
          acphyss,   acphysm,   acphysf
#endif


! EMK NUWRF...Treat this as local
  real, dimension(its:ite,kts:kte,jts:jte) :: dbz

!JJS  9/30/2009 for Steve's new improvement

!  LOGICAL, EXTERNAL :: wrf_dm_on_monitor
  integer  ::  ihalmos
  real     ::  tslopes, tslopeg
  real     ::  xnsplnt, xmsplnt
  real     ::  hmtemp1, hmtemp2, hmtemp3, hmtemp4
  real     ::  ftnsQ, ftngQ, fexp
  real     ::  xssi, fssi 
  real     ::  dmicrons, dmicrong, dvair, alpha
  real,   dimension (its:ite,jts:jte) :: tairN, tairI,    &
                                          ftns,  ftng,    &
                                         pihms, pihmg,    &
                                          pimm,  pcfr

! for Xiping's new dbz code       
  real     :: a_c, a_i            
  real     :: re_c, re_i, re_s
  real     :: w_c, w_i, w_s
  real     :: ze_cld
  
  real,   dimension (its:ite,jts:jte) :: ftns0, ftng0
  real,   dimension (its:ite,jts:jte) :: vi
  integer  :: improve1
  data improve1/-20/


  real     ::  tairc5, hfact, sfact, yy1
  real     ::  xncld, esat, rv, rlapse_m
  real     ::  delT, bhi, cpi
  real     ::  rc, ra, cna
  real     ::  xccld, xknud, cunnf, diffar

  real     :: r7rf, r8rf, r9rf, r16rf
  real     :: r11t, r19t, r19at
  real     :: r30t, r33t
  real     :: r101r, r102rf, r191r, r192rf
  real     :: rn1s
  real     :: r231r, r232rf
  real     :: ami20

#if ( WRF_CHEM == 1)
! JJS 20110525 vvvvv
! for inline Gocart coupling
  INTEGER,      INTENT(IN   )    ::   gid
  INTEGER, PARAMETER :: num_go = 14  ! number of the gocart aerosol species
  REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_go), intent(in) :: aero
  REAL, DIMENSION( ims:ime, kms:kme, jms:jme), intent(out) :: icn_diag !IN concentration [#/Litre]
  REAL, DIMENSION( ims:ime, kms:kme, jms:jme), intent(out) :: nc_diag !cloud concentration [#/cm3]
  integer,intent(in) :: chem_opt ! EMK
  integer,intent(in) :: gsfcgce_gocart_coupling ! EMK

! Local Variables
!  REAL, DIMENSION( its:ite, jts:jte, kts:kte) :: icn_cgs !IN concentration [#/cm3]
!  REAL, DIMENSION( its:ite, jts:jte, kts:kte) :: nc_cgs !cloud concentration [#/cm3]

!  real, dimension( kms:kme, tgmx ) ::  aero_k    ! interporated aerosol mass conc [g/m3]
  real :: e_sat, e_dry  !saturated and dry air water vapor [hPa, mb]
  real :: rh_rad     ! relative humidity [%]
  real :: super_sat  !super saturation [%]
  real :: ccn_out  ! CCN conc [#/cm3]
  real :: icn_out  ! IN conc [#/Litter]
!  real :: rho_dryair ! dry air density [g/m3]  !
!  real :: L_cloud    ! cloud water [g/cm3] !
  real :: P_liu_daum  ! autoconversion rate [g/cm3 s-1]    !
  real :: re_liu_daum ! effective radius of cloud [micron]   !
  real,parameter :: min_icn = 0.01 !minimum # conc of IN [#/Litre]

! JJS 20110525 ^^^^^
#endif

!JJS 20140226  variables for the calculation of effective radius of cloud species
  REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN)   :: XLAND
  real, parameter :: roqi = 0.9179    ! ice density
!  real, parameter :: ccn_over_land = 1500  ! [#/cm3] climatological value
!  real, parameter :: ccn_over_water = 150  ! [#/cm3] climatological value
  real, parameter :: ccn_over_land = 159  ! [#/cm3] climatological value (from MERRA2)
  real, parameter :: ccn_over_water = 66  ! [#/cm3] climatological value (from MERRA2)
  real :: L_cloud    ! cloud water [g/cm3] !
  real :: I_cloud    ! cloud water [g/cm3] !
  real :: mu, ccn_ref, lambda
  real :: gamfac1, gamfac3
!JJS 20140226  ^^^^^
!NUWRF END

!23456789012345678901234567890123456789012345678901234567890123456789012

!
!JJS 1/3/2008  vvvvv
!JJS   the following common blocks have been moved to the top of
!JJS   module_mp_gsfcgce_driver.F

!  real, dimension (31)   ::      aa1,  aa2
!  data aa1/.7939e-7, .7841e-6, .3369e-5, .4336e-5, .5285e-5,     &
!           .3728e-5, .1852e-5, .2991e-6, .4248e-6, .7434e-6,     &
!           .1812e-5, .4394e-5, .9145e-5, .1725e-4, .3348e-4,     &
!           .1725e-4, .9175e-5, .4412e-5, .2252e-5, .9115e-6,     &
!           .4876e-6, .3473e-6, .4758e-6, .6306e-6, .8573e-6,     &
!           .7868e-6, .7192e-6, .6513e-6, .5956e-6, .5333e-6,     &
!           .4834e-6/
!  data aa2/.4006, .4831, .5320, .5307, .5319,      &
!           .5249, .4888, .3894, .4047, .4318,      &
!           .4771, .5183, .5463, .5651, .5813,      &
!           .5655, .5478, .5203, .4906, .4447,      &
!           .4126, .3960, .4149, .4320, .4506,      &
!           .4483, .4460, .4433, .4413, .4382,      &
!           .4361/

!JJS 1/3/2008  ^^^^^

!+---+-----------------------------------------------------------------+
#ifdef EXT_DIAG
  REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT):: refl_10cm  ! GT
! TYPE (WRFU_Clock):: grid_clock
! TYPE (WRFU_Alarm), POINTER:: grid_alarms(:)


!      REAL, DIMENSION(kts:kte):: qv1d, t1d, p1d, qr1d, qs1d, qg1d, dBZ_GT
      LOGICAL, OPTIONAL, INTENT(IN) :: diagflag
      INTEGER, OPTIONAL, INTENT(IN) :: do_radar_ref
#endif
!+---+-----------------------------------------------------------------+
      integer, parameter :: rewflag = 2
      ! 1 : default
      ! 2 : mapping spectrum, different over land and ocean
      real :: mdc1,mdc2,mdc3,mdc4,mdc5,mdc6,mvdc
      real :: efd1,efd2,efd3,efd4,efd5,efd6,efdc
      real :: ltk,ltk2
      real :: lqc,lqc2

      integer, parameter :: reiflag = 5
      ! 1 : default
      ! 2 : Wyser 1998
      ! 3 : Heymsfild et al. 2014
      ! 4 : Mitchell et al. 2011
      ! 5 : mapping spectrum, different over land and ocean
      real :: reimin, reimax
      real :: iwc_0, bw98
      real :: md22, bd22, sigma, cd22, xd22
      real :: lqi,lqi2,efdi

      ! calculate allowable ice supersautration :
      real, intent(in) :: xlat
      real :: d2r, arg

      ! calculate solar declination angle :
      real, intent(in) :: sdec

#ifdef sat_predict
      ! saturation prediction scheme :
      integer :: hid
      real :: rhoair, xlv, xls, xlf, cpm1, dv1, abw, abi
      real :: taui, tauc, taur, tau, atem
      real :: C1, K1, ncloud, nact, qcmax, mvrc
      real :: qimax, nice, inhgr, rhoi, mvdi, mvri
      real :: lqr, lqr2, mvdr, efdr, kmin, kmax, kdxr, afar,  &
              tnr, lzr, bvr, avr, mur, rhoaj, gr2, gbr25
      real :: cnd1, cnd2, dep1, dep2, fez1, fez2, latr, ern1, ern2
      real, dimension (its:ite,jts:jte) :: pact, fez

      ! transition zone :
      real :: dltd, nbd, sbd, nbd1, nbd2, sbd1, sbd2, latint
      real :: fxlat
#endif

      !cp=1.004e7
      !alv=2.5e10 ; alf=3.336e9 ; als=2.8336e10
      !rw=4.615e6 ; cw=4.187e7 ; ci=2.093e7
      !avcp=alv/cp*pir
      real, parameter :: cvap = 1.846e+7    !specific heat capacity of vapor (in CGS)
      real, parameter :: cliq = 4.218e+7    !specific heat capacity of liquid (in CGS)
      real, parameter :: cice = 2.106e+7    !specific heat capacity of ice (in CGS)
      real :: cpm, hlv, hls, hlf
!
!JJS20090623      save  

!    i24h=nint(86400./dt)
!    if (mod(itimestep,i24h).eq.1) then
#ifdef EXT_DIAG
    if (itimestep.eq.1) then
       do k = kts, kte
          do j = jts, jte
             do i= its, ite
!             physc(i,k,j)=0.
!             physe(i,k,j)=0.
!             physd(i,k,j)=0.
!             physs(i,k,j)=0.
!             physf(i,k,j)=0.
!             physm(i,k,j)=0.
             acphysc(i,k,j)=0.
             acphyse(i,k,j)=0.
             acphysd(i,k,j)=0.
             acphyss(i,k,j)=0.
             acphysf(i,k,j)=0.
             acphysm(i,k,j)=0.
             enddo !i
          enddo !j
       enddo !k
!       write(6, *) 'in satice at itimestep=',itimestep
!       write(6, *) ' accumulate latent heating variables have been intialized to 0. '
!       print *,'improve, new_ice_sat = ', improve, new_ice_sat
!       print *,'c610 = ', c610
    endif
#endif

!JJS  convert from mks to cgs, and move from WRF grid to GCE grid
      do k=kts,kte
         do j=jts,jte
         do i=its,ite
         rho(i,j,k)=rho_mks(i,k,j)*0.001
         p0(i,j,k)=p0_mks(i,k,j)*10.0
         pi(i,j,k)=pi_mks(i,k,j)
         ww1(i,j,k)=w_mks(i,k,j)*100.
         dpt(i,j,k)=ptwrf(i,k,j)
         dqv(i,j,k)=qvwrf(i,k,j)
         qcl(i,j,k)=qlwrf(i,k,j)
         qrn(i,j,k)=qrwrf(i,k,j)
         qci(i,j,k)=qiwrf(i,k,j)
         qcs(i,j,k)=qswrf(i,k,j)
         qcg(i,j,k)=qgwrf(i,k,j)
         enddo !i
         enddo !j
      enddo !k

      do k=kts,kte
         do j=jts,jte
         do i=its,ite
            fv(i,j,k)=sqrt(rho(i,j,1)/rho(i,j,k))
         enddo !i
         enddo !j
      enddo !k
!JJS

!
!     ******   THREE CLASSES OF ICE-PHASE   (LIN ET AL, 1983)  *********

         d2t=dt

!C  TAO 2007 START
!   ICE2=0 ! default, 3ice with loud ice, snow and graupel
!              r2is=1., r2ig=1.
!   ICE2=1 ! 2ice with cloud ice and snow (no graupel) - r2iceg=1, r2ice=0.
!              r2is=1., r2ig=0.
!   ICE2=2 ! 2ice with cloud ice and graupel (no snow) - r2ice=1, r2iceg=0.
!              r2is=0., r2ig=1.
!c
!        r2ice=1.
!        r2iceg=1.
         r2ig=1.
         r2is=1.
          if (ice2 .eq. 1) then
!              r2ice=0.
!              r2iceg=1.
              r2ig=0.
              r2is=1.
          endif
          if (ice2 .eq. 2) then
!              r2ice=1.
!              r2iceg=0.
              r2ig=1.
              r2is=0.
          endif

          If (improve .eq. 3) then
              ihail = 0
          endif

!C  TAO 2007 END
     
!JJS  10/7/2008
!   ICE2=3 ! no ice, warm rain only
    iwarm = 0
    if (ice2 .eq. 3) iwarm = 1

      cmin=1.e-19
        If (improve .eq. 3) cmin=1.e-40
        If (improve .eq. 3) cmin=1.e-20 ! EMK NUWRF Bug fix
      cmin1=1.e-20
      cmin2=1.e-40
      ucor=3071.29/tnw**0.75
      ucos=687.97*roqs**0.25/tns**0.75
      ucog=687.97*roqg**0.25/tng**0.75
      uwet=4.464**0.95

      CPI=4.*ATAN(1.)

   ! maximum allowable ice supersaturation :
!      xssi = 0.10   !default
      xssi = 0.05

   ! xssi can be varied with latitude  :
   !  for abs(xlat)<=30  : cos(arg)=0 , xssi=0.10
   !  for abs(xlat)>=60  : cos(arg)=1 , xssi=0.05
!      d2r = CPI/180.0
!      arg = max( min( abs(xlat),60.0 ),30.0 )
!      arg = (60.0 - arg) * 3.0 * d2r
!      xssi = 0.10 - 0.05 * ( cos(arg) )**2.0

!   ??????????
      tslopes=0.10539656   ! increase tns by 40 from -5 to -40C
      tslopeg=0.08559235   ! increase tng by 20 from -5 to -40C


!  HALLET-MOSSOP RIME SPLINTERING parameters
      ihalmos=1
        if (improve.gt.2) ihalmos=1
      xnsplnt=370.     ! peak # splinters per milligram of rime
      xmsplnt=4.4e-8   ! mass of a splinter (from Ferrier 1994)
      hmtemp1=-2.
      hmtemp2=-4.
      hmtemp3=-6.
      hmtemp4=-8.

!      rijl2 = 1. / (ide-ids) / (jde-jds)

      do j=jts,jte
         do i=its,ite
         it(i,j)=1
         enddo
      enddo

      f2=rd1*d2t
      f3=rd2*d2t

      ft=dt/d2t
!      rft=rijl2*ft
!      a0=.5*istatmin*rijl2
      rt0=1./(t0-t00)

      bw3=bw+3.
      bs3=bs+3.
      bg3=bg+3.
      bsh5=2.5+bsh
      bgh5=2.5+bgh
      bwh5=2.5+bwh
      bw6=bw+6.
      bs6=bs+6.
      betah=.5*beta

      r10t=rn10*d2t
      r11at=rn11a*d2t
      r19bt=rn19b*d2t
      r20t=-rn20*d2t
      r23t=-rn23*d2t
      r25a=rn25

!     ami50 for use in PINT
      ami50=3.76e-8
      ami100=1.51e-7
      ami40=2.41e-8

      ami50=3.84e-6                       ! 6/15/02 tao's
      ami40=3.08e-8                       ! 6/15/02 tao's
      ami50=4.8e-7*(100./50.)**3          ! Roger: actually = 3.84e-6 

      if (improve .eq. 3) then

!         itaobraun=0 ! using original way for pint and pidep
!         itaobraun=1 ! using scott braun's way for pint and pidep
       
         itaobraun=1  
         new_ice_sat=3 
         rdt=1./d2t
         r11t=rn11*d2t
         r19t=rn19*d2t
         r19at=rn19a*d2t
         r20t=rn20*d2t
         r23t=rn23*d2t
         r30t=rn30*d2t
         r33t=rn33*d2t
         ami50=4.8e-7
         ami40=2.46e-7
!
!      if ( wrf_dm_on_monitor() .and. itimestep.eq.1 ) then
!      if ( itimestep.eq.1 ) then
!         print *,'inside satice improve=3'
!         print *,'itaobraun = ', itaobraun
!      endif

         Rc=1.e-3               ! cloud droplet radius 10 microns
         Ra=1.e-5               ! aerosol radius 0.1 microns
         Cna=500.               ! contact nuclei conc per cc
         Bhi=1.01e-2            ! pollen (Deihl et al. 2006)

      endif

!       itaobraun=0 ! original pint and pidep & see Tao and Simpson 1993
        itaobraun=1 ! see Tao et al. (2003)
!
       if ( itaobraun.eq.0 ) then
          cn0=1.e-8
!c        beta=-.6
       else if ( itaobraun.eq.1 ) then
          cn0=1.e-6
!         cn0=1.e-8  ! special
!c        beta=-.46
       endif

!C    ******************************************************************

  do 1000 k=kts,kte
       kp=k+1
       tb0=0.
       qb0=0.

     do 2000 j=jts,jte
     do 2000 i=its,ite

! EMK BUG FIX...Initialize variable
#if ( WRF_CHEM == 1)
        ccn_out = 0.e0 
        icn_out = 0.e0 
        if ( (chem_opt == 300 .or. chem_opt == 301 .or. &
              chem_opt == 302 .or. chem_opt == 303) .and. &
              (gsfcgce_gocart_coupling == 1) ) then
           icn_diag(i,k,j) = icn_out  ! #/Litre
           nc_diag(i,k,j) = ccn_out  ! #/cm3
        else
           icn_diag(i,k,j) = 0.
           nc_diag(i,k,j) = 0.
        end if

#endif
! EMK END

         rp0=3.799052e3/p0(i,j,k)
         pi0=pi(i,j,k)
         pir=1./(pi(i,j,k))
         pr0=1./p0(i,j,k)
         r00=rho(i,j,k)
         r0s=sqrt(rho(i,j,k))
         rr0=1./rho(i,j,k)
         rrs=sqrt(rr0)
         rrq=sqrt(rrs)
         f0(i,j,k)=al/cp/pi(i,j,k)
         f00=f0(i,j,k)
         fv0=fv(i,j,k)
         fvs=sqrt(fv(i,j,k))
         zrr=1.e5*zrc*rrq
         zsr=1.e5*zsc*rrq
         zgr=1.e5*zgc*rrq
         cp409=c409*pi0
         cv409=c409*avc
         cp580=c580*pi0
         cs580=c580*asc
         alvr=r00*alv
         afcp=afc*pir
         avcp=avc*pir
         ascp=asc*pir
         vrcf=vrc*fv0
         vscf=vsc*fv0
         vgcf=vgc*fv0
         vgcr=vgc*rrs
         dwvp=c879*pr0

         r3f=rn3*fv0

         r4f=rn4*fv0
         r5f=rn5*fv0
         r6f=rn6*fv0
         r7r=rn7*rr0
         r8r=rn8*rr0
         r9r=rn9*rr0
         r101f=rn101*fvs
         r10ar=rn10a*r00
         r11rt=rn11*rr0*d2t
         r12r=rn12*r00
         r14r=rn14*rrs
         r14f=rn14*fv0
         r15r=rn15*rrs
         r15ar=rn15a*rrs
         r15f=rn15*fv0
         r15af=rn15a*fv0
         r16r=rn16*rr0
         r17r=rn17*rr0
         r17aq=rn17a*rrq
         r17as=rn17a*fvs
         r18r=rn18*rr0
         r19rt=rn19*rr0*d2t
         r19aq=rn19a*rrq
         r19as=rn19a*fvs
         r20bq=rn20b*rrq
         r20bs=rn20b*fvs
         r22f=rn22*fv0
         r23af=rn23a*fvs
         r23br=rn23b*r00
         r25rt=rn25*rr0*d2t
         r31r=rn31*rr0
         r32rt=rn32*d2t*rrs


!JJS added 10/1/2009
         if (improve .eq. 3) then
    
            r7r=rn7*rr0*fv0
            r8r=rn8*rr0*fv0

            r9rf=rn9*rr0*fv0
            r16rf=rn16*rr0*fv0
            r101r=rn101*rr0
            r102rf=rn102*rrs*fvs
            r191r=rn191*rr0
            r192rf=rn192*rrs*fvs

            r331r=rn331*rr0
            r332rf=rn332*rrs*fvs
            r34f=rn34*fv0

	       r231r=rn231*rr0
	       r232rf=rn232*rrs*fvs
!           xccld=xncld*r00               !cloud number concentration

        endif

        pt(i,j)=dpt(i,j,k)
        qv(i,j)=dqv(i,j,k)
        qc(i,j)=qcl(i,j,k)
        qr(i,j)=qrn(i,j,k)
        qi(i,j)=qci(i,j,k)
        qs(i,j)=qcs(i,j,k)
        qg(i,j)=qcg(i,j,k)
!        IF (QV(I,J)+QB0 .LE. 0.) QV(I,J)=-QB0
         if (qc(i,j) .le.  cmin) qc(i,j)=0.0
         if (qr(i,j) .le.  cmin) qr(i,j)=0.0
         if (qi(i,j) .le.  cmin) qi(i,j)=0.0
         if (qs(i,j) .le.  cmin) qs(i,j)=0.0
         if (qg(i,j) .le.  cmin) qg(i,j)=0.0
        tair(i,j)=(pt(i,j)+tb0)*pi0
        tairc(i,j)=tair(i,j)-t0
         zr(i,j)=zrr
         zs(i,j)=zsr
         zg(i,j)=zgr
         vr(i,j)=0.0 
         vs(i,j)=0.0
         vg(i,j)=0.0
         vi(i,j)=0.0

            ftns(i,j)=1.
            ftng(i,j)=1.
            ftns0(i,j)=1.
            ftng0(i,j)=1.

         cnd(i,j)=0.0
         dep(i,j)=0.
         ern(i,j)=0.0
         pint(i,j)=0.0
         pidep(i,j)=0.0

         psdep(i,j)=0.
         pgdep(i,j)=0.
         dd1(i,j)=0.
         dd(i,j)=0.
          pgsub(i,j)=0.
          psmlt(i,j)=0.
          pgmlt(i,j)=0.
          pimlt(i,j)=0.
          psacw(i,j)=0.
          piacr(i,j)=0.

          pssub(i,j)=0.0
          pgsub(i,j)=0.0

          psfw(i,j)=0.0
	     psfi(i,j)=0.0
          pidep(i,j)=0.0

          pgfr(i,j)=0.
          psacr(i,j)=0.
          wgacr(i,j)=0.
          pihom(i,j)=0.
          pidw(i,j)=0.0
          
          psaut(i,j)=0.0
          psaci(i,j)=0.0
          praci(i,j)=0.0
          pwacs(i,j)=0.0
          qsacw(i,j)=0.0
          
          pracs(i,j)=0.0
	     qracs(i,j)=0.0
	     qsacr(i,j)=0.0
	     pgaut(i,j)=0.0
 
          praut(i,j)=0.0
          pracw(i,j)=0.0
          pgfr(i,j)=0.0

          qracs(i,j)=0.0

          pgacs(i,j)=0.0
          qgacw(i,j)=0.0
          dgaci(i,j)=0.0
          dgacs(i,j)=0.0
          wgacs(i,j)=0.0
	     wgaci(i,j)=0.0
          dgacw(i,j)=0.0
          dgacr(i,j)=0.
          pgwet(i,j)=0.0

          qgacr(i,j)=0.0

	     pihom(i,j)=0.0
	     pimlt(i,j)=0.0
	     pidw(i,j)=0.0
          pimm(i,j)=0.0
          pcfr(i,j)=0.0

          pihms(i,j)=0.0 
          pihmg(i,j)=0.0
          ftns(i,j)=1.
          ftng(i,j)=1.
          pmlts(i,j)=0.0
          pmltg(i,j)=0.0

            dlt4(i,j)=0.0
            dlt3(i,j)=0.0
            dlt2(i,j)=0.0

!     ******************************************************************
!     ***   Y1 : DYNAMIC VISCOSITY OF AIR (U)
!     ***   DWV : DIFFUSIVITY OF WATER VAPOR IN AIR (PI)
!     ***   TCA : THERMAL CONDUCTIVITY OF AIR (KA)
!     ***   Y2 : KINETIC VISCOSITY (V)
 
            y1(i,j)=c149*tair(i,j)**1.5/(tair(i,j)+120.)
            dwv(i,j)=dwvp*tair(i,j)**1.81
            tca(i,j)=c141*y1(i,j)
            scv(i,j)=1./((rr0*y1(i,j))**.1666667*dwv(i,j)**.3333333)

!JJS 10/7/2008     vvvvv
    IF (IWARM .EQ. 1) THEN
!JJS   for calculating processes related to warm rain only
                qi(i,j)=0.0
                qs(i,j)=0.0
                qg(i,j)=0.0

                if (qr(i,j) .gt. cmin1) then
                   dd(i,j)=r00*qr(i,j)
                   y1(i,j)=dd(i,j)**.25
                   zr(i,j)=zrc/y1(i,j)
                endif

!                call vqrqi(1,improve,r00,fv0,qr(i,j),vr(i,j))
                call vtr_mks(rho_mks(i,k,j),qr(i,j),tair(i,j),vr(i,j))  !in MKS
                vr(i,j) = vr(i,j) * 100.  !in CGS

!* 21 * PRAUT   AUTOCONVERSION OF QC TO QR                        **21**
!* 22 * PRACW : ACCRETION OF QC BY QR                             **22**
                pracw(i,j)=0.
                praut(i,j)=0.0
                  if (improve .ge. 3 .or. improve .eq. -1) then  !Di Wu
                     praut(i,j)=max(rn21*(qc(i,j)-bnd21),0.0)

                     y1(i,j)=1./zr(i,j)
                     y2(i,j)=y1(i,j)*y1(i,j)
                     y3(i,j)=y1(i,j)*y2(i,j)
                     y4(i,j)=r22f*qc(i,j)*y3(i,j)*(rn50+rn51*y1(i,j)+  &
                             rn52*y2(i,j)+rn53*y3(i,j))
                     pracw(i,j)=max(y4(i,j), 0.0)
                     if (qr(i,j) .le. cmin) pracw(i,j)=0.
                 else                                 ! below is 2007 saticel_s
                     pracw(i,j)=r22f*qc(i,j)/zr(i,j)**bw3
                     y1(i,j)=qc(i,j)-bnd3
                     if (y1(i,j).gt.0.0) then
                        praut(i,j)=r00*y1(i,j)*y1(i,j)/(1.2e-4+rn21/y1(i,j))
                     endif
                  endif

!C********   HANDLING THE NEGATIVE CLOUD WATER (QC)    ******************
                 Y1(I,J)=QC(I,J)/D2T
                 PRAUT(I,J)=MIN(Y1(I,J), PRAUT(I,J))
                 PRACW(I,J)=MIN(Y1(I,J), PRACW(I,J))
                 Y1(I,J)=(PRAUT(I,J)+PRACW(I,J))*D2T
               
               if (qc(i,j) .lt. y1(i,j) .and. y1(i,j) .ge. cmin2) then
                   y2(i,j)=qc(i,j)/(y1(i,j)+cmin2)
                   praut(i,j)=praut(i,j)*y2(i,j)
                   pracw(i,j)=pracw(i,j)*y2(i,j)
                   qc(i,j)=0.0
               else
                  qc(i,j)=qc(i,j)-y1(i,j)
               endif
               
               PR(I,J)=(PRAUT(I,J)+PRACW(I,J))*D2T
               QR(I,J)=QR(I,J)+PR(I,J)
                        
!*****   TAO ET AL (1989) SATURATION TECHNIQUE  ***********************
           
              cnd(i,j)=0.0
              tair(i,j)=(pt(i,j)+tb0)*pi0
              y1(i,j)=1./(tair(i,j)-c358)
              qsw(i,j)=rp0*exp(c172-c409*y1(i,j))
              dd(i,j)=cp409*y1(i,j)*y1(i,j)
              dm(i,j)=qv(i,j)+qb0-qsw(i,j)
              cnd(i,j)=dm(i,j)/(1.+avcp*dd(i,j)*qsw(i,j))
!c    ******   condensation or evaporation of qc  ******
              cnd(i,j)=max(-qc(i,j), cnd(i,j))
              pt(i,j)=pt(i,j)+avcp*cnd(i,j)
              qv(i,j)=qv(i,j)-cnd(i,j)
              qc(i,j)=qc(i,j)+cnd(i,j)

!* 23 * ERN : EVAPORATION OF QR (SUBSATURATION)                   **23**
            ern(i,j)=0.0

            if (qr(i,j).gt.0.0) then
               tair(i,j)=(pt(i,j)+tb0)*pi0
               rtair(i,j)=1./(tair(i,j)-c358)
               qsw(i,j)=rp0*exp(c172-c409*rtair(i,j))
               ssw(i,j)=(qv(i,j)+qb0)/qsw(i,j)-1.0
               dm(i,j)=qv(i,j)+qb0-qsw(i,j)
               rsub1(i,j)=cv409*qsw(i,j)*rtair(i,j)*rtair(i,j)
               dd1(i,j)=max(-dm(i,j)/(1.+rsub1(i,j)), 0.0)
               y1(i,j)=.78/zr(i,j)**2+r23af*scv(i,j)/zr(i,j)**bwh5
               y2(i,j)=r23br/(tca(i,j)*tair(i,j)**2)+1./(dwv(i,j) &
                       *qsw(i,j))
!cccc
               ern(i,j)=r23t*ssw(i,j)*y1(i,j)/y2(i,j)
               ern(i,j)=min(dd1(i,j),qr(i,j),max(ern(i,j),0.))
 
!JJS 20090708 vvvvv
! xli - reducing evaporation rate according to fitq
               if (qr(i,j) .gt. 1.e-6 ) then
                  fact_fit=0.11*(qr(i,j)*1.e3)**(-1.27)+0.98
                  ern(i,j)=ern(i,j)/fact_fit
               endif
!JJS 20090708 ^^^^^

               pt(i,j)=pt(i,j)-avcp*ern(i,j)
               qv(i,j)=qv(i,j)+ern(i,j)
               qr(i,j)=qr(i,j)-ern(i,j)
            endif

  ELSE       ! part of if (iwarm.eq.1) then
!JJS 10/7/2008     ^^^^^

!JJS   for calculating processes related to both ice and warm rain

!     ***   COMPUTE ZR,ZS,ZG,VR,VS,VG      *****************************

            if (qr(i,j) .gt. cmin) then
	            dd(i,j)=r00*qr(i,j)
	            y1(i,j)=sqrt(dd(i,j))
	            y2(i,j)=sqrt(y1(i,j))
	            zr(i,j)=zrc/y2(i,j)
!               if (improve.gt.2) then
!                  call vqrqi(1,improve,r00,fv0,qr(i,j),vr(i,j))
!               else
!                  vr(i,j)=max(vrcf*dd(i,j)**bwq, 0.)
!               endif
               call vtr_mks(rho_mks(i,k,j),qr(i,j),tair(i,j),vr(i,j))  !in MKS
               vr(i,j) = vr(i,j) * 100.  !in CGS
            endif

            if (qs(i,j) .gt. cmin) then
               dd(i,j)=r00*qs(i,j)
               y1(i,j)=dd(i,j)**.25
               ftns(i,j)=1.
               if (improve.gt.2) then
                  call sgmap(1,qs(i,j),r00,tairc(i,j),ftns0(i,j))
                  ftns(i,j)=ftns0(i,j)**0.25
               endif
               zs(i,j)=zsc/y1(i,j)*ftns(i,j)
!               if (improve.gt.2) ftns(i,j)=ftns0(i,j)**bsq
!               vs(I,J)=MAX(vscf*dd(I,J)**bsq/ftns(i,j), 0.)
               call vts_mks(improve,rho_mks(i,k,j),qs(i,j),tair(i,j),vs(i,j))  !in MKS
               vs(i,j) = vs(i,j) * 100.  !in CGS
            endif

            if (qg(i,j) .gt. cmin) then
               dd(i,j)=r00*qg(i,j)
               y1(i,j)=dd(i,j)**.25
               if (ihail .eq. 1) vgcf=vgcr
               ftng(i,j)=1.
               if (improve.gt.2) then
                  call sgmap(2,qg(i,j),r00,tairc(i,j),ftng0(i,j))
                  ftng(i,j)=ftng0(i,j)**0.25
               endif
               zg(i,j)=zgc/y1(i,j)*ftng(i,j)
!               if (improve.gt.2)ftng(i,j)=ftng0(i,j)**bgq
!               vg(i,j)=max(vgcf*dd(i,j)**bgq/ftng(i,j), 0.0)
               call vtg_mks(ihail,improve,rho_mks(i,k,j),qg(i,j),tair(i,j),vg(i,j))  !in MKS
               vg(i,j) = vg(i,j) * 100.  !in CGS
            endif

!            call vqrqi(2,improve,r00,fv0,qi(i,j),vi(i,j))
            call vti_mks(improve,rho_mks(i,k,j),tair(i,j),qi(i,j),qv(i,j),p0_mks(i,k,j),xland(i,j),vi(i,j))  !in MKS
            vi(i,j) = vi(i,j) * 100.  !in CGS

            if (qr(i,j) .le. crmin) vr(i,j)=0.0
            if (qs(i,j) .le. csmin) vs(i,j)=0.0
            if (qg(i,j) .le. cgmin) vg(i,j)=0.0
            if (qi(i,j) .le. cimin) vi(i,j)=0.0

!     ******************************************************************
!     ***   Y1 : DYNAMIC VISCOSITY OF AIR (U)
!     ***   DWV : DIFFUSIVITY OF WATER VAPOR IN AIR (PI)
!     ***   TCA : THERMAL CONDUCTIVITY OF AIR (KA)
!     ***   Y2 : KINETIC VISCOSITY (V)
!
! JJS 20120315:  move the following 4 lines out of "if (iwarm)" 
!            y1(i,j)=c149*tair(i,j)**1.5/(tair(i,j)+120.)
!            dwv(i,j)=dwvp*tair(i,j)**1.81
!            tca(i,j)=c141*y1(i,j)
!            scv(i,j)=1./((rr0*y1(i,j))**.1666667*dwv(i,j)**.3333333)

!*  1 * PSAUT : AUTOCONVERSION OF QI TO QS                        ***1**
!*  3 * PSACI : ACCRETION OF QI TO QS                             ***3**
!*  4 * PSACW : ACCRETION OF QC BY QS (RIMING) (QSACW FOR PSMLT)  ***4**
!*  5 * PRACI : ACCRETION OF QI BY QR                             ***5**
!*  6 * PIACR : ACCRETION OF QR OR QG BY QI                       ***6**
!* 34 * pwacs : collection of qs by qc                            **34**

          pihms(i,j)=0.0 
          pihmg(i,j)=0.0 
          psaut(i,j)=0.0
          psaci(i,j)=0.0
          praci(i,j)=0.0
          piacr(i,j)=0.0
          psacw(i,j)=0.0
          pwacs(i,j)=0.0
          qsacw(i,j)=0.0
          ftns(i,j)=1.
          ftng(i,j)=1.
          if (improve.gt.2) then
             ftns(i,j)=ftns0(i,j)
             ftng(i,j)=ftng0(i,j)
          endif
      if (improve .lt. 2) then ! below is 2007 saticel_s

         dd(i,j)=1./zs(i,j)**bs3
!         rn1=9.4e-15
         if (tair(i,j).lt.t0) then
            esi(i,j)=exp(.025*tairc(i,j))
            psaut(i,j)=r2is*max(rn1*esi(i,j)*(qi(i,j)-bnd1) ,0.0)
            psaci(i,j)=r2is*r3f*esi(i,j)*qi(i,j)*dd(i,j)
!    to cut water to snow accretion by half
!            PSACW(I,J)=R4F*QC(I,J)*DD(I,J)
!            psacw(i,j)=r2is*0.5*r4f*qc(i,j)*dd(i,j)
            psacw(i,j)=r2is*r4f*qc(i,j)*dd(i,j)
            praci(i,j)=r2is*r5f*qi(i,j)/zr(i,j)**bw3
            piacr(i,j)=r2is*r6f*qi(i,j)*(zr(i,j)**(-bw6))
!JJS           PIACR(I,J)=R6F*QI(I,J)/ZR(I,J)**BW6
         else
            qsacw(i,j)=r2is*r4f*qc(i,j)*dd(i,j)
         endif
!
       else
!
          if (tair(i,j).lt.t0) then

#ifdef sat_predict
             rn1s=1.e-3
             bnd1=1.e-4
             esi(i,j)=exp(0.025*tairc(i,j))
             psaut(i,j)=r2is*max(rn1s*esi(i,j)*(qi(i,j)-bnd1*fv0*fv0) ,0.0)
#else
!             y1(i,j)=rdt*(qi(i,j)-r1r*exp(beta*tairc(i,j)))
!             psaut(i,j)=max(y1(i,j),0.0)
             rn1s=1.e-3
             bnd1=6.e-4
             esi(i,j)=exp(.025*tairc(i,j))
             if (improve.gt.2) esi(i,j)=0.15
             psaut(i,j)=r2is*max(rn1s*esi(i,j)*(qi(i,j)-bnd1*fv0*fv0) ,0.0) 
#endif
	     esi(i,j)=1.0 
             dmicrons=(r00*qs(i,j)/roqs/cpi/(tns*ftns(i,j)))**.25*1.e4
             if (improve.gt.2) esi(i,j)=min(1.,(dmicrons/1500.)**4.) ! f(dmicrons)

             y1(i,j)=1.0
             if (improve.gt.2.and.vs(i,j).gt.0.) y1(i,j)=abs((vs(i,j)-vi(i,j))  &
                                                 /vs(i,j))
	     psaci(i,j)=y1(i,j)*r3f*qi(i,j)/zs(i,j)**bs3*ftns(i,j)*esi(i,j)
             psacw(i,j)=r4f*qc(i,j)/zs(i,j)**bs3*ftns(i,j)
!improvepsacw(i,j)=r4f*qc(i,j)/zs(i,j)**bs3*ftns(i,j)*esi(i,j)
             if (ihalmos.eq.1)then
                y2(i,j)=0.
                if((tairc(i,j).le.hmtemp1).and.(tairc(i,j).ge.hmtemp4))  &
                                                         y2(i,j)=0.5
!                if((tairc(i,j).ge.hmtemp2).and.(tairc(i,j).le.hmtemp3))  &
!vvvvvvvv tao 20110722 vvvvvvvvvvvvvvvvvv
                if((tairc(i,j).le.hmtemp2).and.(tairc(i,j).ge.hmtemp3))  &
!^^^^^^^^ Tao 20110722 ^^^^^^^^^^^^^^^^^^
                                                         y2(i,j)=1.
                pihms(i,j)=psacw(i,j)*y2(i,j)*xnsplnt*1000.*xmsplnt
                psacw(i,j)=psacw(i,j)-pihms(i,j)
             endif
             pwacs(i,j)=r34f*qc(i,j)/zs(i,j)**bs6*ftns(i,j)
! improvepwacs(i,j)=r34f*qc(i,j)/zs(i,j)**bs6*ftns(i,j)*esi(i,j)
             y1(i,j)=1./zr(i,j)
             y2(i,j)=y1(i,j)*y1(i,j)
             y3(i,j)=y1(i,j)*y2(i,j)
             y5(i,j)=1.0
           if (improve.gt.2.and.vr(i,j).gt.0.) y5(i,j)=abs((vr(i,j)-vi(i,j))   &
                                                         /vr(i,j))
	     dd(i,j)=y5(i,j)*r5f*qi(i,j)*y3(i,j)*(rn50+rn51*y1(i,j)             &
                                        +rn52*y2(i,j)+rn53*y3(i,j))
             praci(i,j)=max(dd(i,j),0.0)
	     y4(i,j)=y3(i,j)*y3(i,j)
             dd1(i,j)=y5(i,j)*r6f*qi(i,j)*y4(i,j)*(rn60+rn61*y1(i,j)       &
                                     +rn62*y2(i,j)+rn63*y3(i,j))

             piacr(i,j)=max(dd1(i,j),0.0)
          else
             qsacw(i,j)=r4f*qc(i,j)/zs(i,j)**bs3*ftns(i,j)
          endif   !tairc 

       endif                         ! for Processes 1, 3, 4, 5, 6, & 34

!23456789012345678901234567890123456789012345678901234567890123456789012
!* 21 * PRAUT   AUTOCONVERSION OF QC TO QR                        **21**
!* 22 * PRACW : ACCRETION OF QC BY QR                             **22**

!          if (improve .eq. 3) then
                                     ! Steve's new improvement 9/21/2009
!vvvvvvvvvvvvvv Tao 20110722 vvvvvvvvvvvvvv
          if (improve .eq. 3 .or. improve .eq. -1) then
                                     ! Steve's new improvement 9/21/2009
!^^^^^^^^^^^^^^ Tao 20110722 ^^^^^^^^^^^^^^

#if ( WRF_CHEM == 1)
!JJS 20110602 vvvvv
             ! EMK...Only execute when GOCART and coupling are selected
             if ( (chem_opt == 300 .or. chem_opt == 301 .or. &
                   chem_opt == 302 .or. chem_opt == 303) .and. &
                   (gsfcgce_gocart_coupling == 1) ) then 
                !      sat vapor pressure [hPa,mb]
                e_sat = 6.11 * exp( 5423. *( 1.0/273.15 - 1./tair(i,j) ) )
                !      dry air vapor pressure [hPa, mb]
                e_dry = qv(i,j) / ( qv(i,j) + 0.622 ) * p0(i,j,k) * 1.e-3  ! p0 in [g*cm/s2/cm2]
                rh_rad = max(1.e-6, e_dry/e_sat*100.)  ! relative humidity [%]
                super_sat = max(0.001, rh_rad - 100.e0) !super saturation [%]
                !
                ! convert gocart aerosol mass conc to CN
                !
                call mass2ccn(tair(i,j),super_sat,aero(i,k,j,:),ccn_out )
                !             nc_cgs(i,j,k)  = max(100., ccn_out)  !diagnostic cloud droplet conc  [#/cm3]
                ccn_out = max(100., ccn_out)
!                if (mod(itimestep,i24h).eq.1) then
!                   if (ccn_out > 100.) print *,'in satice, i,j,k,ccn_out = ', i,j,k, ccn_out
!                endif
                ! rho_dryair = p0(i,j,k) / ( tair(i,j) * 2.87 * 1.0e6) ! dry air density [g/cm3]
                L_cloud = qc(i,j) * rho(i,j,k)             ! cloud water [g/cm3]
                 !  g/g        g/cm3
                 !             call auto_conversion( L_cloud, nc_cgs(i,j,k), P_liu_daum, re_liu_daum )
                call auto_conversion( L_cloud, ccn_out, P_liu_daum, re_liu_daum )
                praut(i,j) = P_liu_daum / rho(i,j,k)  !autoconversion rate [g/g s-1]
                !             if (i .eq. int((ite+its)/2) .and. j .and. int((jte+jts)/2) )   &
                !             if (rh_rad .ge. 100.)  then
                !                print 1078,'i','j', 'k', 'e_sat', 'e_dry', 'rh_rad' , 'super_sat', 'nc_cgs', 'rho_dryair', &
                !                       'L_cloud', 'P_liu_daum', 'praut'
                !             1078 format(3a4, 9a11)
                !                print 1079, i, j, k, e_sat, e_dry, rh_rad, super_sat, nc_cgs(i,j,k), rho_dryair, &
                !                       L_cloud, P_liu_daum, praut(i,j)
                !             1079 format(3i4, 9e11.4)
                !             endif
             else
                praut(i,j)=max(rn21*(qc(i,j)-bnd21),0.0)
             end if ! if (gsfcgce_gocart_coupling == 1)
#else
             praut(i,j)=max(rn21*(qc(i,j)-bnd21),0.0)
!JJS 20110602 ^^^^^
#endif

             y1(i,j)=1./zr(i,j)
             y2(i,j)=y1(i,j)*y1(i,j)
             y3(i,j)=y1(i,j)*y2(i,j)
             y4(i,j)=r22f*qc(i,j)*y3(i,j)*(rn50+rn51*y1(i,j)+  &
                     rn52*y2(i,j)+rn53*y3(i,j))
             pracw(i,j)=max(y4(i,j), 0.0) 
          else                                 ! below is 2007 saticel_s
             pracw(i,j)=r22f*qc(i,j)/zr(i,j)**bw3
             y1(i,j)=qc(i,j)-bnd3
             if (y1(i,j).gt.0.0) then
                praut(i,j)=r00*y1(i,j)*y1(i,j)/(1.2e-4+rn21/y1(i,j))
             endif
          endif                                  ! for Processes 21 & 22

!* 12 * PSFW : BERGERON PROCESSES FOR QS (KOENING, 1971)          **12**
!* 13 * PSFI : BERGERON PROCESSES FOR QS                          **13**

          pidep(i,j)=0.0

#ifndef sat_predict
!>>> Note that Bergeron processes are concerned in saturation prediction scheme
         if (improve.eq.3) then
!        if (improve1.eq.3) then
                                     ! Steve's new improvement 9/21/2009
            if (tair(i,j) .lt. t0) then
               y1(i,j)=max( min(tairc(i,j), -1.), -31.)
               it(i,j)=int(abs(y1(i,j)))
               y1(i,j)=rn12a(it(i,j))
               y2(i,j)=rn12b(it(i,j))
               y3(i,j)=rn13(it(i,j))
               psfw(i,j)=r2is*max(d2t*y1(i,j)*(y2(i,j)+r12r*qc(i,j))*  &
                            qi(i,j), 0.0)
               psfi(i,j)=r2is*y3(i,j)*qi(i,j)
!
               y4(i,j)=1./(tair(i,j)-c358)
               y5(i,j)=1./(tair(i,j)-c76)
               qsw(i,j)=rp0*exp(c172-c409*y4(i,j))
               qsi(i,j)=rp0*exp(c218-c580*y5(i,j))
#ifdef new_saturation
               esw(i,j) = min(0.99*p0_mks(i,k,j),esw_mks(tair(i,j)))
               esi(i,j) = min(0.99*p0_mks(i,k,j),esi_mks(tair(i,j)))
               if ( esi(i,j) .gt. esw(i,j) ) esi(i,j) = esw(i,j)
               qsw(i,j) = 0.622*esw(i,j)/(p0_mks(i,k,j)-esw(i,j))
               qsi(i,j) = 0.622*esi(i,j)/(p0_mks(i,k,j)-esi(i,j))
               esw(i,j) = esw(i,j)*10.  !in CGS
               esi(i,j) = esi(i,j)*10.  !in CGS
#endif
               ! EMK...Prevent division by zero
!               hfact=(qv(i,j)+qb0-qsi(i,j))/(qsw(i,j)-qsi(i,j))
               hfact=(qv(i,j)+qb0-qsi(i,j))/(qsw(i,j)-qsi(i,j)+cmin1)
               if(hfact.gt.1.) hfact=1.
               sfact=1
               SSI(i,j)=(qv(i,j)+qb0)/qsi(i,j)-1.

               fssi=min(xssi,max(.0,xssi*(tairc(i,j)+44.)/(44.0-38.0)))  !max ssi f(tair)

               fssi=min(ssi(i,j),fssi)

!  STEVE : PLEASE CHECK
!                 xssi=min(ssi(i,j), fssi)
!                  xssi=min(ssi(i,j), 0.20)
               r_nci=min(1.e-3*exp(-.639+12.96*fssi), 1.)  !meyers et al. 
!  STEVE : PLEASE CHECK
               if (tairc(i,j).le.-5.) then                         !meyers
                  r_nci=max(1.e-3*exp(-.639+12.96*fssi),0.528e-3)  !meyers et al. 
               else
                  r_nci=min(1.e-3*exp(-.639+12.96*fssi),0.528e-3)  !meyers et al. 
               endif  !tairc

               if (r_nci.gt.15.) r_nci=15.                  !cap at 15000/liter

#if ( WRF_CHEM == 1)
               ! EMK...Only execute when GOCART and coupling is turned on.
               if ( (chem_opt == 300 .or. chem_opt == 301 .or. &
                     chem_opt == 302 .or. chem_opt == 303) .and. &
                     (gsfcgce_gocart_coupling == 1) ) then
                  !JJS 20110602 vvvvv
                  ! Conversion rate of cloud water to ice in the Bergeron porcess based on Meyer + DeMott formulae
                  !
                  ! convert gocart aerosol mass conc to IN
                  !      p0 need to be converted from g*cm/s2/cm2 to mb (hPa)
                  !                call mass2icn(p0(i,j,k)*0.001,tair(i,j),aero(i,k,j,:), icn_out)
                  call mass2icn(p0(i,j,k)*0.001,tair(i,j),aero(i,k,j,:), icn_out,i,j,k)
                  icn_out = min(1.e3, max(0.01e0 ,  icn_out) )
!                  if ( wrf_dm_on_monitor() .and. itimestep.eq.1 ) &
!                       print *,'mass2icn is called'
                  !                icn_cgs(i,j,k) = max(min_icn, icn_out) * 1.e-3   !IN conc [#/cm3] <-- [#/Litter]
                  !                icn_cgs(i,j,k) = icn_out * 1.e-3 !IN conc [#/cm3]
                  !                if (rh_rad .ge. 100.) then
                  !                   print 1080, 'i', 'j', 'k', 'p0', 'tair', 'icn_cgs'
                  !                1080 format(3a4, 4a12)         
                  !                   print 1081, i,j,k, p0(i,j,k)*0.001,tair(i,j), icn_cgs(i,j,k)
                  !                1081 format(3i4, 3e12.4)
                  !                endif
                  !                r_nci = icn_cgs(i,j,k)  !DeMotto's formuale
                  r_nci = icn_out * 1.e-3  !DeMotto's formuale
                  !JJS 20110602 ^^^^^
               end if ! if (gsfcgce_gocart_coupling == 1)
#endif

               dd(i,j)=min((r00*qi(i,j)/r_nci), ami40)   !mean cloud ice mass
               y4(i,j)=1.-aa2(it(i,j))
               sfact=(AMI50**Y4(i,j)-AMI40**Y4(i,j))/               &
                    (AMI50**Y4(i,j)-dd(i,j)**Y4(i,j))
               if (hfact.gt.0.) then
                  psfi(i,j)=r2is*psfi(i,j)*hfact*sfact
               else
                  psfi(i,j)=0.
               endif
            endif   !tair(i,j)
            
         else                               ! below is 2007 saticel_s

            if (tair(i,j).lt.t0.and.qi(i,j).gt.cmin) then
               y1(i,j)=max( min(tairc(i,j), -1.), -31.)
               it(i,j)=int(abs(y1(i,j)))
               y1(i,j)=rn12a(it(i,j))
               y2(i,j)=rn12b(it(i,j))  	
               psfw(i,j)=r2is*max(d2t*y1(i,j)*(y2(i,j)+r12r*qc(i,j))*qi(i,j),0.0)
               rtair(i,j)=1./(tair(i,j)-c76)
               y2(i,j)=exp(c218-c580*rtair(i,j))
               qsi(i,j)=rp0*y2(i,j)
               esi(i,j)=c610*y2(i,j)
#ifdef new_saturation
               esi(i,j) = min(0.99*p0_mks(i,k,j),esi_mks(tair(i,j)))
               qsi(i,j) = 0.622*esi(i,j)/(p0_mks(i,k,j)-esi(i,j))
               esi(i,j) = esi(i,j)*10.  !in CGS
#endif
               ssi(i,j)=(qv(i,j)+qb0)/qsi(i,j)-1.
               r_nci=min(1.e-6*exp(-.46*tairc(i,j)),1.)
!               R_NCI=min(1.e-8*EXP(-.6*TAIRC(I,J)),1.) ! use Tao's
               dm(i,j)=max( (qv(i,j)+qb0-qsi(i,j)), 0.)
               rsub1(i,j)=cs580*qsi(i,j)*rtair(i,j)*rtair(i,j)
               y3(i,j)=1./tair(i,j)
               dd(i,j)=y3(i,j)*(rn30a*y3(i,j)-rn30b)+rn30c*tair(i,j)/esi(i,j)
               y1(i,j)=206.18*ssi(i,j)/dd(i,j)
               pidep(i,j)=y1(i,j)*sqrt(r_nci*qi(i,j)/r00)
               dep(i,j)=dm(i,j)/(1.+rsub1(i,j))/d2t

               if (dm(i,j).gt.cmin2) then
                  a2=1.
                  if (pidep(i,j).gt.dep(i,j).and.pidep(i,j).gt.cmin2) then
                     a2=dep(i,j)/pidep(i,j)
                     pidep(i,j)=dep(i,j)
                  endif
!                  psfi(i,j)=r2is*a2*.5*qi(i,j)*y1(i,j)/(sqrt(ami100) &
                  psfi(i,j)=r2is*a2*qi(i,j)*y1(i,j)/(sqrt(ami100) &
                          -sqrt(ami40))
               else if (dm(i,j).lt.-cmin2) then
!
!        SUBLIMATION TERMS USED ONLY WHEN SATURATION ADJUSTMENT FOR ICE
!        IS TURNED OFF
!
                        pidep(i,j)=0.
                        psfi(i,j)=0.
                    else
                        pidep(i,j)=0.
                        psfi(i,j)=0.
               endif  !dm(i,j)
!
            endif !tair(i,j)
!        
          endif                                 ! for Processes 12 & 13
!<<< end of ifndef sat_predict
#endif

!TTT***** QG=QG+MIN(PGDRY,PGWET)
!*  9 * PGACS : ACCRETION OF QS BY QG (DGACS,WGACS: DRY AND WET)  ***9**
!* 14 * DGACW : ACCRETION OF QC BY QG (QGACW FOR PGMLT)           **14**
!* 16 * DGACR : ACCRETION OF QR TO QG (QGACR FOR PGMLT)           **16**
!*******PGDRY : DGACW+DGACI+DGACR+DGACS                           ******
!* 15 * DGACI : ACCRETION OF QI BY QG (WGACI FOR WET GROWTH)      **15**
!* 17 * PGWET : WET GROWTH OF QG                                  **17**
!*  Steve turned off PGWET, set PGWET = 0.
!*  Steve turned off wgaci, set wgaci = 0.                           

      if (improve .eq. 3) then

         y1(i,j)=abs( vg(i,j)-vs(i,j) )
         y2(i,j)=zs(i,j)*zg(i,j)
         y3(i,j)=5./y2(i,j)
         y4(i,j)=.08*y3(i,j)*y3(i,j)
         y5(i,j)=.05*y3(i,j)*y4(i,j)
         y2(i,j)=y1(i,j)*(y3(i,j)/zs(i,j)**5+y4(i,j)/zs(i,j)**3        &
                         +y5(i,j)/zs(i,j))

         pgacs(i,j)=r2ig*r2is*r9rf*y2(i,j)*ftns(i,j)*ftng(i,j)
!         dgacs(i,j)=pgacs(i,j)
         dgacs(i,j)=0.0             !Lang et al. 2007
!crh     wgacs(i,j)=10.*r9rf*y2(i,j)
         wgacs(i,j)=0.0
         y1(i,j)=1./zg(i,j)**bg3
         esi(i,j)=1.0 !egc constant in consatrh via r14f/rn14; use esi( ) to make f(T)/f(q)

         dmicrong=(r00*qg(i,j)/roqg/cpi/(tng*ftng(i,j)))**.25*1.e4
         esi(i,j)=min(1.,(dmicrong/500.)**1.1)       ! f(dmicrons)

         dgacw(i,j)= r2ig*esi(i,j)*r14f*qc(i,j)*y1(i,j)*ftng(i,j)
!         dgacw(i,j)=r2ig*r14f*qc(i,j)*y1(i,j)*ftng(i,j)
         y2(i,j)=0.
         if ((tairc(i,j).le.hmtemp1).and.(tairc(i,j).ge.hmtemp4))  &
                                                    y2(i,j)=0.5
!         if ((tairc(i,j).ge.hmtemp2).and.(tairc(i,j).le.hmtemp3))  &
!vvvvvvvvvv Tao 20110722 vvvvvvvvvvvvv
         if ((tairc(i,j).le.hmtemp2).and.(tairc(i,j).ge.hmtemp3))  &
!^^^^^^^^^^ Tao 20110722 ^^^^^^^^^^^^^
                                                    y2(i,j)=1.
         pihmg(i,j)=r2ig*dgacw(i,j)*y2(i,j)*xnsplnt*1000.*xmsplnt
         dgacw(i,j)=r2ig*dgacw(i,j)-pihmg(i,j)
         qgacw(i,j)=r2ig*dgacw(i,j)
         y5(i,j)=1.0
         if (vg(i,j).gt.0.) y5(i,j)=abs((vg(i,j)-vi(i,j))/vg(i,j))
         dgaci(i,j)= r2ig*y5(i,j)*r15f*qi(i,j)*y1(i,j)*ftng(i,j)
!         dgaci(i,j)=r2ig*r15f*qi(i,j)*y1(i,j)*ftng(i,j)
         dgaci(i,j)=0.0
!crh     wgaci(i,j)=r15af*qi(i,j)*y1(i,j)
         wgaci(i,j)=0.0

         y1(i,j)=abs( vg(i,j)-vr(i,j) )
         y2(i,j)=zr(i,j)*zg(i,j)
         y3(i,j)=5./y2(i,j)
         y4(i,j)=.08*y3(i,j)*y3(i,j)
         y5(i,j)=.05*y3(i,j)*y4(i,j)
         dd(i,j)=r16rf*y1(i,j)*(y3(i,j)/zr(i,j)**5+y4(i,j)/zr(i,j)**3 &
                 +y5(i,j)/zr(i,j))*ftng(i,j)
         dgacr(i,j)=r2ig*max(dd(i,j),0.0)
         qgacr(i,j)=dgacr(i,j)

         if (tair(i,j) .ge. t0) then
            dgacs(i,j)=0.0
!crh            wgacs(i,j)=0.0
            dgacw(i,j)=0.0
            dgaci(i,j)=0.0
!crh            wgaci(i,j)=0.0
            dgacr(i,j)=0.0
         else
            pgacs(i,j)=0.0
            qgacw(i,j)=0.0
            qgacr(i,j)=0.0
         endif

         pgwet(i,j)=0.0

      else                           ! below is 2007 saticel_s

         if (qc(i,j)+qr(i,j).lt.1.e-4) then
            ee1=.01
         else
            ee1=1.
         endif
         ee2=0.09
         egs(i,j)=ee1*exp(ee2*tairc(i,j))
!            EGS(I,J)=0.1 ! 6/15/02 tao's
         if (tair(i,j).ge.t0) egs(i,j)=1.0
         y1(i,j)=abs(vg(i,j)-vs(i,j))
         y2(i,j)=zs(i,j)*zg(i,j)
         y3(i,j)=5./y2(i,j)
         y4(i,j)=.08*y3(i,j)*y3(i,j)
         y5(i,j)=.05*y3(i,j)*y4(i,j)
         dd(i,j)=y1(i,j)*(y3(i,j)/zs(i,j)**5+y4(i,j)/zs(i,j)**3 &
                   +y5(i,j)/zs(i,j))
         pgacs(i,j)=r2ig*r2is*r9r*egs(i,j)*dd(i,j)
!
         if (ihail.eq.1) then
            dgacs(i,j)=pgacs(i,j)
         else
            dgacs(i,j)=0.
         endif
!
         wgacs(i,j)=r2ig*r2is*r9r*dd(i,j)
         y1(i,j)=1./zg(i,j)**bg3

         if (ihail .eq. 1) then
            dgacw(i,j)=r2ig*max(r14r*qc(i,j)*y1(i,j), 0.0)

         else
            dgacw(i,j)=r2ig*max(r14f*qc(i,j)*y1(i,j), 0.0)
         endif

            qgacw(i,j)=dgacw(i,j)
            y1(i,j)=abs(vg(i,j)-vr(i,j))
            y2(i,j)=zr(i,j)*zg(i,j)
            y3(i,j)=5./y2(i,j)
            y4(i,j)=.08*y3(i,j)*y3(i,j)
            y5(i,j)=.05*y3(i,j)*y4(i,j)
            dd(i,j)=r16r*y1(i,j)*(y3(i,j)/zr(i,j)**5+y4(i,j)/zr(i,j)**3 &
                    +y5(i,j)/zr(i,j))
            dgacr(i,j)=r2ig*max(dd(i,j), 0.0)
            qgacr(i,j)=dgacr(i,j)

            if (tair(i,j).ge.t0) then
               dgacs(i,j)=0.0
               wgacs(i,j)=0.0
               dgacw(i,j)=0.0
               dgacr(i,j)=0.0
            else
               pgacs(i,j)=0.0
               qgacw(i,j)=0.0
               qgacr(i,j)=0.0
            endif

            dgaci(i,j)=0.0
            wgaci(i,j)=0.0
            pgwet(i,j)=0.0

            if (tair(i,j).lt.t0) then
               y1(i,j)=qi(i,j)/zg(i,j)**bg3
               if (ihail.eq.1) then
                  dgaci(i,j)=r2ig*r15r*y1(i,j)
                  wgaci(i,j)=r2ig*r15ar*y1(i,j)
               else
!                   DGACI(I,J)=r2ig*R15F*Y1(I,J)
                  dgaci(i,j)=0.
                  wgaci(i,j)=r2ig*r15af*y1(i,j)
!                   WGACI(I,J)=0.  ! 6/15/02 tao's
               endif  !ihail
!
               if (tairc(i,j).ge.-50.) then
                  if (alf+rn17c*tairc(i,j) .eq. 0.)  &
                    write(91,*) itimestep, i,j,k, alf, rn17c, tairc(i,j)
                  y1(i,j)=1./(alf+rn17c*tairc(i,j))
                  if (ihail.eq.1) then
                     y3(i,j)=.78/zg(i,j)**2+r17aq*scv(i,j)/zg(i,j)**bgh5
                  else
                     y3(i,j)=.78/zg(i,j)**2+r17as*scv(i,j)/zg(i,j)**bgh5
                  endif
                  y4(i,j)=alvr*dwv(i,j)*(rp0-(qv(i,j)+qb0)) &
                          -tca(i,j)*tairc(i,j)
                  dd(i,j)=y1(i,j)*(r17r*y4(i,j)*y3(i,j) &
                        +(wgaci(i,j)+wgacs(i,j))*(alf+rn17b*tairc(i,j)))
                  pgwet(i,j)=r2ig*max(dd(i,j), 0.0)
               endif !tairc
            endif !tair

      endif                         ! for Processes 9, 14, 15, 16 & 17

 !********   HANDLING THE NEGATIVE CLOUD WATER (QC)    ******************

        y1(i,j)=qc(i,j)/d2t
        psacw(i,j)=min(y1(i,j), psacw(i,j))
        praut(i,j)=min(y1(i,j), praut(i,j))
        pracw(i,j)=min(y1(i,j), pracw(i,j))
        psfw(i,j)= min(y1(i,j), psfw(i,j))
        dgacw(i,j)=min(y1(i,j), dgacw(i,j))
        qsacw(i,j)=min(y1(i,j), qsacw(i,j))
        qgacw(i,j)=min(y1(i,j), qgacw(i,j))
        pihms(i,j)=min(y1(i,j), pihms(i,j))
        pihmg(i,j)=min(y1(i,j), pihmg(i,j))

        y1(i,j)=d2t*(psacw(i,j)+praut(i,j)+pracw(i,j)+psfw(i,j)       &
           +dgacw(i,j)+qsacw(i,j)+qgacw(i,j)+pihms(i,j)+pihmg(i,j))

        qc(i,j)=qc(i,j)-y1(i,j)
!
        if (qc(i,j) .lt. 0.0) then
           y2(i,j)=1.
           if (y1(i,j) .ne. 0.) y2(i,j)=qc(i,j)/y1(i,j)+1.
           psacw(i,j)=psacw(i,j)*y2(i,j)
           praut(i,j)=praut(i,j)*y2(i,j)
           pracw(i,j)=pracw(i,j)*y2(i,j)
           psfw(i,j)=psfw(i,j)*y2(i,j)
           dgacw(i,j)=dgacw(i,j)*y2(i,j)
           qsacw(i,j)=qsacw(i,j)*y2(i,j)
           qgacw(i,j)=qgacw(i,j)*y2(i,j)
           pihms(i,j)=pihms(i,j)*y2(i,j)
           pihmg(i,j)=pihmg(i,j)*y2(i,j)
           qc(i,j)=0.0
        endif
!        wgacr(i,j)=0.

!        if (improve .lt. 2) then
!vvvvvvvvvv Tao 20110722 vvvvvvvvvvvvvvvvvv
          if (improve .eq.3 .or. improve .eq. -1) then
             wgacr(i,j)= qgacr(i,j)+qgacw(i,j)
          else
!^^^^^^^^^^ Tao 20110722 ^^^^^^^^^^^^^^^^^^

!******** SHED PROCESS (WGACR=PGWET-DGACW-WGACI-WGACS)
!c
            wgacr(i,j)=pgwet(i,j)-dgacw(i,j)-wgaci(i,j)-wgacs(i,j)
            y2(i,j)=dgacw(i,j)+dgaci(i,j)+dgacr(i,j)+dgacs(i,j)

            if (pgwet(i,j).ge.y2(i,j)) then
               wgacr(i,j)=0.0
               wgaci(i,j)=0.0
               wgacs(i,j)=0.0
            else
               dgacr(i,j)=0.0
               dgaci(i,j)=0.0
               dgacs(i,j)=0.0
            endif
        endif

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!********   HANDLING THE NEGATIVE CLOUD ICE (QI)      ******************
        y1(i,j)=qi(i,j)/d2t
        psaut(i,j)=min(y1(i,j), psaut(i,j))
        psaci(i,j)=min(y1(i,j), psaci(i,j))
        praci(i,j)=min(y1(i,j), praci(i,j))
        psfi(i,j)= min(y1(i,j), psfi(i,j))
        dgaci(i,j)=min(y1(i,j), dgaci(i,j))
        wgaci(i,j)=min(y1(i,j), wgaci(i,j))

!       Steve: Please check

        y1(i,j)=d2t*(psaut(i,j)+psaci(i,j)+praci(i,j)+psfi(i,j)       &
                    +dgaci(i,j)+wgaci(i,j)-pihms(i,j)-pihmg(i,j))

        qi(i,j)=qi(i,j)-y1(i,j)
!
        if (qi(i,j) .lt. 0.0) then
           y2(i,j)=1.
           if (y1(i,j) .ne. 0.0) y2(i,j)=qi(i,j)/y1(i,j)+1.
           psaut(i,j)=psaut(i,j)*y2(i,j)
           psaci(i,j)=psaci(i,j)*y2(i,j)
           praci(i,j)=praci(i,j)*y2(i,j)
           psfi(i,j)=psfi(i,j)*y2(i,j)
           dgaci(i,j)=dgaci(i,j)*y2(i,j)
           wgaci(i,j)=wgaci(i,j)*y2(i,j)
           qi(i,j)=0.0
        endif

!
        if (improve .eq. 3) then
!            qi(i,j)=qi(i,j)+d2t*(pihms(i,j)+pihmg(i,j))
            wgacr(i,j)=qgacr(i,j)+qgacw(i,j)
            dlt3(i,j)=0.0
              if (qr(i,j) .lt. 1.e-4) dlt3(i,j)=1.
            dlt4(i,j)=1.
!              if (qc(i,j) .gt. 5.e-4) dlt4(i,j)=0.0
!              if (qs(i,j) .le. 1.e-4) dlt4(i,j)=1.
            if (qc(i,j) .gt. 1.e-3) dlt4(i,j)=0.0
	       if (qs(i,j) .le. 1.e-4) dlt4(i,j)=1.

            if (tair(i,j) .ge. t0) then
               dlt3(i,j)=0.0
               dlt4(i,j)=0.0
            endif
            pr(i,j)=d2t*(qsacw(i,j)+praut(i,j)+pracw(i,j)+wgacr(i,j)   &
                     -qgacr(i,j))
            ps(i,j)=d2t*(psaut(i,j)+psaci(i,j)+dlt4(i,j)*psacw(i,j)    &
                   +psfw(i,j)+psfi(i,j)+dlt3(i,j)*praci(i,j))
            pg(i,j)=d2t*((1.-dlt3(i,j))*praci(i,j)+dgaci(i,j)          &
                   +wgaci(i,j)+dgacw(i,j)+(1.-dlt4(i,j))*psacw(i,j))
!
         else
            dlt3(i,j)=0.0
            dlt2(i,j)=0.0
!
!            DLT4(I,J)=1.0
!            if(qc(i,j) .gt. 5.e-4) dlt4(i,j)=0.0
!            if(qs(i,j) .le. 1.e-4) dlt4(i,j)=1.0
!
!            IF (TAIR(I,J).ge.T0) THEN
!               DLT4(I,J)=0.0
!            ENDIF

             if (tair(i,j).lt.t0) then
                if (qr(i,j).lt.1.e-4) then
                   dlt3(i,j)=1.0
                   dlt2(i,j)=1.0
                endif
                if (qs(i,j).ge.1.e-4) then
                   dlt2(i,j)=0.0
                endif
             endif
            
             if (ice2 .eq. 1) then
                dlt3(i,j)=1.0
                dlt2(i,j)=1.0
             endif

            pr(i,j)=(qsacw(i,j)+praut(i,j)+pracw(i,j)+qgacw(i,j))*d2t
            ps(i,j)=(psaut(i,j)+psaci(i,j)+psacw(i,j)+psfw(i,j) &
                    +psfi(i,j)+dlt3(i,j)*praci(i,j))*d2t
            pg(i,j)=((1.-dlt3(i,j))*praci(i,j)+dgaci(i,j)+wgaci(i,j) &
                    +dgacw(i,j))*d2t
         endif

!*  7 * PRACS : ACCRETION OF QS BY QR                             ***7**
!*  8 * PSACR : ACCRETION OF QR BY QS (QSACR FOR PSMLT)           ***8**
!*  2 * PGAUT : AUTOCONVERSION OF QS TO QG                        ***2**
!* 18 * PGFR : FREEZING OF QR TO QG                               **18**

       qracs(i,j)=0.0
       y1(i,j)=abs(vr(i,j)-vs(i,j))
       y2(i,j)=zr(i,j)*zs(i,j)
       y3(i,j)=5./y2(i,j)
       y4(i,j)=.08*y3(i,j)*y3(i,j)
       y5(i,j)=.05*y3(i,j)*y4(i,j)
       pracs(i,j)=r2ig*r2is*r7r*y1(i,j)*(y3(i,j)/zs(i,j)**5 &
                  +y4(i,j)/zs(i,j)**3+y5(i,j)/zs(i,j))*ftns(i,j)
       psacr(i,j)=r2is*r8r*y1(i,j)*(y3(i,j)/zr(i,j)**5 &
                  +y4(i,j)/zr(i,j)**3+y5(i,j)/zr(i,j))*ftns(i,j)
       qsacr(i,j)=psacr(i,j)
!       if (improve.eq.3) qracs(i,j)=min(d2t*pracs(i,j), qs(i,j))
!vvvvvvvvvvvvv Tao 20110722 vvvvvvvvvvvvvvvvvv
       if (improve.eq.3 .or.improve .eq. -1) qracs(i,j)=min(d2t*pracs(i,j), qs(i,j))
!^^^^^^^^^^^^^ Tao 20110722 ^^^^^^^^^^^^^^^^^^

!
          pgaut(i,j)=0.0
          pgfr(i,j)=0.0
       if (tair(i,j) .lt. t0) then
            y2(i,j)=exp(rn18a*(t0-tair(i,j)))
            temp = 1./zr(i,j)
            temp = temp*temp*temp*temp*temp*temp*temp
         pgfr(i,j)=r2ig*max(r18r*(y2(i,j)-1.)*temp, 0.0)
       endif

       if (tair(i,j).ge.t0) then
           pracs(i,j)=0.0
           psacr(i,j)=0.0
       else
           qsacr(i,j)=0.0
           qracs(i,j)=0.0
       endif
!!
!     endif  ! for Processes 2, 7, 8, & 18

!********   HANDLING THE NEGATIVE RAIN WATER (QR)    *******************
!********   HANDLING THE NEGATIVE SNOW (QS)          *******************

!      if (improve .eq. 3) then
!vvvvvvvvvvvvvvv Tao 20110722 vvvvvvvvvvvvvvv
      if (improve .eq. 3 .or. improve .eq. -1) then
!^^^^^^^^^^^^^^^ Tao 20110722 ^^^^^^^^^^^^^^^ 
#ifdef use_cpm
          cpm = cp + cvap*qv(i,j) + cliq*(qc(i,j)+qr(i,j)) + &
                cice*(qi(i,j)+qs(i,j)+qg(i,j))
          hlv = alv - (cliq-cvap)*(tair(i,j)-t0)
          hlf = alf - (cice-cliq)*(tair(i,j)-t0)
          hls = hlv + hlf
          avcp = hlv/cpm*pir
          ascp = hls/cpm*pir
          afcp = hlf/cpm*pir
#endif
         y1(i,j)=qr(i,j)/d2t
          piacr(i,j)=min(y1(i,j), piacr(i,j))
          dgacr(i,j)=min(y1(i,j), dgacr(i,j))
          psacr(i,j)=min(y1(i,j), psacr(i,j))
          pgfr(i,j)= min(y1(i,j), pgfr(i,j))
          y1(i,j)=(piacr(i,j)+dgacr(i,j)+psacr(i,j)+pgfr(i,j))*d2t
          qr(i,j)=qr(i,j)+pr(i,j)+qracs(i,j)-y1(i,j)
          if (qr(i,j) .lt. 0.0) then
             y2(i,j)=1.
             if (y1(i,j) .ne. 0.0) y2(i,j)=qr(i,j)/y1(i,j)+1.
             piacr(i,j)=piacr(i,j)*y2(i,j)
             dgacr(i,j)=dgacr(i,j)*y2(i,j)
             pgfr(i,j)=pgfr(i,j)*y2(i,j)
             psacr(i,j)=psacr(i,j)*y2(i,j)
             qr(i,j)=0.0
          endif
          dlt2(i,j)=1.
          if (qr(i,j) .gt. 1.e-4) dlt2(i,j)=0.
          if (qs(i,j) .le. 1.e-4) dlt2(i,j)=1.
          if (tair(i,j) .ge. t0) dlt2(i,j)=0.
          y1(i,j)=qs(i,j)/d2t
          pgacs(i,j)=min(y1(i,j), pgacs(i,j))
          dgacs(i,j)=min(y1(i,j), dgacs(i,j))
          wgacs(i,j)=min(y1(i,j), wgacs(i,j))
          pgaut(i,j)=min(y1(i,j), pgaut(i,j))
          pracs(i,j)=min(y1(i,j), pracs(i,j))
          pwacs(i,j)=min(y1(i,j), pwacs(i,j))
          prn(i,j)=d2t*((1.-dlt3(i,j))*piacr(i,j)+dgacr(i,j)+pgfr(i,j) &
                      +(1.-dlt2(i,j))*psacr(i,j))
          ps(i,j)=ps(i,j)+d2t*(dlt3(i,j)*piacr(i,j)+dlt2(i,j)*psacr(i,j))
          pracs(i,j)=(1.-dlt2(i,j))*pracs(i,j)
          pwacs(i,j)=(1.-dlt4(i,j))*pwacs(i,j)
 
          psn(i,j)=d2t*(pgacs(i,j)+dgacs(i,j)+wgacs(i,j)+pgaut(i,j)    &
                   +pracs(i,j)+pwacs(i,j))
 
          qs(i,j)=qs(i,j)+ps(i,j)-qracs(i,j)-psn(i,j)
          if (qs(i,j) .lt. 0.0) then
             y2(i,j)=1.
             if (psn(i,j) .ne. 0.) y2(i,j)=qs(i,j)/psn(i,j)+1.
             pgacs(i,j)=pgacs(i,j)*y2(i,j)
             dgacs(i,j)=dgacs(i,j)*y2(i,j)
             wgacs(i,j)=wgacs(i,j)*y2(i,j)
             pgaut(i,j)=pgaut(i,j)*y2(i,j)
             pracs(i,j)=pracs(i,j)*y2(i,j)
             pwacs(i,j)=pwacs(i,j)*y2(i,j)
             qs(i,j)=0.0
          endif
          psn(i,j)=d2t*(pgacs(i,j)+dgacs(i,j)+wgacs(i,j)+pgaut(i,j)    &
                    +pracs(i,j)+pwacs(i,j))
          qg(i,j)=qg(i,j)+pg(i,j)+prn(i,j)+psn(i,j)
          y1(i,j)=d2t*(psacw(i,j)+psfw(i,j)+dgacw(i,j)+piacr(i,j)      &
               +dgacr(i,j)+psacr(i,j)+pgfr(i,j)+pihms(i,j)+pihmg(i,j)) &
               -qracs(i,j)
          pt(i,j)=pt(i,j)+afcp*y1(i,j)
!
       else                                    ! below is 2007 saticel_s
!
          y1(i,j)=qr(i,j)/d2t
          y2(i,j)=-qg(i,j)/d2t
          piacr(i,j)=min(y1(i,j), piacr(i,j))
          dgacr(i,j)=min(y1(i,j), dgacr(i,j))
          wgacr(i,j)=min(y1(i,j), wgacr(i,j))
          wgacr(i,j)=max(y2(i,j), wgacr(i,j))
          psacr(i,j)=min(y1(i,j), psacr(i,j))
          pgfr(i,j)= min(y1(i,j), pgfr(i,j))
          del=0.
          if (wgacr(i,j) .lt. 0.) del=1.
          y1(i,j)=(piacr(i,j)+dgacr(i,j)+(1.-del)*wgacr(i,j) &
                  +psacr(i,j)+pgfr(i,j))*d2t
          qr(i,j)=qr(i,j)+pr(i,j)-y1(i,j)-del*wgacr(i,j)*d2t
          if (qr(i,j) .lt. 0.0) then
             a1=1.
             if (y1(i,j) .ne. 0.) a1=qr(i,j)/y1(i,j)+1.
             piacr(i,j)=piacr(i,j)*a1
             dgacr(i,j)=dgacr(i,j)*a1
             if (wgacr(i,j).gt.0.) wgacr(i,j)=wgacr(i,j)*a1
             pgfr(i,j)=pgfr(i,j)*a1
             psacr(i,j)=psacr(i,j)*a1
             qr(i,j)=0.0
          endif
!
          prn(i,j)=d2t*((1.-dlt3(i,j))*piacr(i,j)+dgacr(i,j) &
                   +wgacr(i,j)+(1.-dlt2(i,j))*psacr(i,j)+pgfr(i,j))
          ps(i,j)=ps(i,j)+d2t*(dlt3(i,j)*piacr(i,j) &
                  +dlt2(i,j)*psacr(i,j))
          pracs(i,j)=(1.-dlt2(i,j))*pracs(i,j)
          y1(i,j)=qs(i,j)/d2t
          pgacs(i,j)=min(y1(i,j), pgacs(i,j))
          dgacs(i,j)=min(y1(i,j), dgacs(i,j))
          wgacs(i,j)=min(y1(i,j), wgacs(i,j))
          pgaut(i,j)=min(y1(i,j), pgaut(i,j))
          pracs(i,j)=min(y1(i,j), pracs(i,j))
          psn(i,j)=d2t*(pgacs(i,j)+dgacs(i,j)+wgacs(i,j) &
                   +pgaut(i,j)+pracs(i,j))
          qs(i,j)=qs(i,j)+ps(i,j)-psn(i,j)
!
          if (qs(i,j).lt.0.0) then
             a2=1.
             if (psn(i,j) .ne. 0.0) a2=qs(i,j)/psn(i,j)+1.
             pgacs(i,j)=pgacs(i,j)*a2
             dgacs(i,j)=dgacs(i,j)*a2
             wgacs(i,j)=wgacs(i,j)*a2
             pgaut(i,j)=pgaut(i,j)*a2
             pracs(i,j)=pracs(i,j)*a2
             psn(i,j)=psn(i,j)*a2
             qs(i,j)=0.0
          endif
!!
!!C           PSN(I,J)=D2T*(PGACS(I,J)+DGACS(I,J)+WGACS(I,J)
!!c                    +PGAUT(I,J)+PRACS(I,J))
          y2(i,j)=d2t*(psacw(i,j)+psfw(i,j)+dgacw(i,j)+piacr(i,j) &
                  +dgacr(i,j)+wgacr(i,j)+psacr(i,j)+pgfr(i,j))
          pt(i,j)=pt(i,j)+afcp*y2(i,j)
          qg(i,j)=qg(i,j)+pg(i,j)+prn(i,j)+psn(i,j)
!
       endif                       ! for handling the nagtive Qr and Qs


!* 11 * PSMLT : MELTING OF QS                                     **11**
!* 19 * PGMLT : MELTING OF QG TO QR                               **19**

       psmlt(i,j)=0.0
       pgmlt(i,j)=0.0
       tair(i,j)=(pt(i,j)+tb0)*pi0
       tairc(i,j)=tair(i,j)-t0

       ftns(i,j)=1.
       ftng(i,j)=1.
       ftns0(i,j)=1.
       ftng0(i,j)=1.

       if (improve.gt.2) then
          call sgmap(1,qs(i,j),r00,tairc(i,j),ftns0(i,j))
          call sgmap(2,qg(i,j),r00,tairc(i,j),ftng0(i,j))
       endif 

       if (tair(i,j).ge.t0) then
          tairc(i,j)=tair(i,j)-t0

          if (improve .eq. 3) then
             dd(i,j)=r11t*tairc(i,j)*(r101r/zs(i,j)**2+r102rf  &
                               /zs(i,j)**bsh5)*ftns0(i,j)
             psmlt(i,j)=r2is*min(qs(i,j),max(dd(i,j),0.0))
!             y2(i,j)=r191r/zg(i,j)**2+r192rf/zg(i,j)**bgh5
!vvvvvvvvvvvvv Tao vvvvvvvvvvvvvvvvvvv
             y2(i,j)=(r191r/zg(i,j)**2+r192rf/zg(i,j)**bgh5)*ftng0(i,j)

!^^^^^^^^^^^^^ Tao ^^^^^^^^^^^^^^^^^^^

!             dd1(i,j)=tairc(i,j)*(r19t*y2(i,j)+r19at*(qgacw(i,j) &
!                                             +qgacr(i,j)))*ftng0(i,j)
!vvvvvvvvvvvvv Tao vvvvvvvvvvvvvvvvvvv
             dd1(i,j)=tairc(i,j)*(r19t*y2(i,j)+r19at*(qgacw(i,j) &
                                             +qgacr(i,j)))
!^^^^^^^^^^^^^ Tao ^^^^^^^^^^^^^^^^^^^
             pgmlt(i,j)=r2ig*min(qg(i,j),max(dd1(i,j),0.0))
          else                                ! below is 2007 saticel_s
             y1(i,j)=tca(i,j)*tairc(i,j)-alvr*dwv(i,j) &
                               *(rp0-(qv(i,j)+qb0))
             y2(i,j)=.78/zs(i,j)**2+r101f*scv(i,j)/zs(i,j)**bsh5
             dd(i,j)=r11rt*y1(i,j)*y2(i,j)+r11at*tairc(i,j) &
                       *(qsacw(i,j)+qsacr(i,j))
             psmlt(i,j)=r2is*max(0.0, min(dd(i,j), qs(i,j)))

             if (ihail.eq.1) then
                y3(i,j)=.78/zg(i,j)**2+r19aq*scv(i,j)/zg(i,j)**bgh5
             else
                y3(i,j)=.78/zg(i,j)**2+r19as*scv(i,j)/zg(i,j)**bgh5
             endif

             dd1(i,j)=r19rt*y1(i,j)*y3(i,j)+r19bt*tairc(i,j) &
                        *(qgacw(i,j)+qgacr(i,j))
             pgmlt(i,j)=r2ig*max(0.0, min(dd1(i,j), qg(i,j)))
          endif  !improve

#ifdef use_cpm
          cpm = cp + cvap*qv(i,j) + cliq*(qc(i,j)+qr(i,j)) + &
                cice*(qi(i,j)+qs(i,j)+qg(i,j))
          hlv = alv - (cliq-cvap)*(tair(i,j)-t0)
          hlf = alf - (cice-cliq)*(tair(i,j)-t0)
          hls = hlv + hlf
          avcp = hlv/cpm*pir
          ascp = hls/cpm*pir
          afcp = hlf/cpm*pir
#endif
          pt(i,j)=pt(i,j)-afcp*(psmlt(i,j)+pgmlt(i,j))
          qr(i,j)=qr(i,j)+psmlt(i,j)+pgmlt(i,j)
          qs(i,j)=qs(i,j)-psmlt(i,j)
          qg(i,j)=qg(i,j)-pgmlt(i,j)

       endif   ! processes 11 & 19
! 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!* 24 * PIHOM : HOMOGENEOUS FREEZING OF QC TO QI (T < T00)        **24**
!* 25 * PIDW : DEPOSITION GROWTH OF QC TO QI ( T0 < T <= T00)     **25**
!* 26 * PIMLT : MELTING OF QI TO QC (T >= T0)                     **26**
!****** PIMM  : IMMERSION FREEZING OF QC TO QI (T < T0)           ******
!****** PCFR  : CONTACT NUCLEATION OF QC TO QI (T < T0)           ******


       if (qc(i,j).le.cmin1) qc(i,j)=0.0
       if (qi(i,j).le.cmin1) qi(i,j)=0.0
       tair(i,j)=(pt(i,j)+tb0)*pi0
       tairc(i,j) = tair(i,j)-t0

       ftns(i,j)=1.
       ftng(i,j)=1.
       ftns0(i,j)=1.
       ftng0(i,j)=1.

       if (improve.gt.2) then
          call sgmap(1,qs(i,j),r00,tairc(i,j),ftns0(i,j))
          call sgmap(2,qg(i,j),r00,tairc(i,j),ftng0(i,j))
       endif 

       if (tair(i,j).le.t00) then
          pihom(i,j)=qc(i,j)
       else
          pihom(i,j)=0.0
       endif 
       if (tair(i,j).ge.t0) then
          pimlt(i,j)=qi(i,j)
       else
          pimlt(i,j)=0.0
       endif
       pidw(i,j)=0.0
           
       if (improve .eq. 3) then
!      if (improve1.eq.3) then
          if (tair(i,j) .lt. t0 .and. tair(i,j) .gt. t00) then
             tairc(i,j) = tair(i,j)-t0
#ifndef sat_predict
!>>> pidw may be already calculated in saturation prediction scheme.
             y1(i,j) = max( min(tairc(i,j), -1.), -31.)
             it(i,j) = int(abs(y1(i,j)))
             y2(i,j)=aa1(it(i,j))
             y3(i,j) = aa2(it(i,j))
             if (tairc(i,j).le.-5.)then                       !  meyers
                y4(i,j) = 1./(tair(i,j)-c358)
                qsw(i,j)=rp0*exp(c172-c409*y4(i,j))
                rtair(i,j)=1./(tair(i,j)-c76)
                y5(i,j)=exp(c218-c580*rtair(i,j))
                qsi(i,j)=rp0*y5(i,j)
#ifdef new_saturation
                esw(i,j) = min(0.99*p0_mks(i,k,j),esw_mks(tair(i,j)))
                esi(i,j) = min(0.99*p0_mks(i,k,j),esi_mks(tair(i,j)))
                if ( esi(i,j) .gt. esw(i,j) ) esi(i,j) = esw(i,j)
                qsw(i,j) = 0.622*esw(i,j)/(p0_mks(i,k,j)-esw(i,j))
                qsi(i,j) = 0.622*esi(i,j)/(p0_mks(i,k,j)-esi(i,j))
                esw(i,j) = esw(i,j)*10.  !in CGS
                esi(i,j) = esi(i,j)*10.  !in CGS
#endif
                SSI(i,j)=(qv(i,j)+qb0)/qsi(i,j)-1.
                fssi=min(xssi,max(.0,xssi*(tairc(i,j)+44.)/(44.0-38.0))) !max ssi f(tair)
                fssi=min(ssi(i,j),fssi)

!                r_nci=min(1.e-3*exp(-.639+12.96*fssi), 0.528e-3)      ! meyers et al.
!              bug fix by JJS 20120429
                r_nci=max(1.e-3*exp(-.639+12.96*fssi), 0.528e-3)      ! meyers et al.
                if (r_nci.gt.15.) r_nci=15.                          !cap at 15000/liter

#if ( WRF_CHEM == 1)
                ! EMK...Only execute if GOCART and coupling is on
                if ( (chem_opt == 300 .or. chem_opt == 301 .or. &
                      chem_opt == 302 .or. chem_opt == 303) .and. &
                      (gsfcgce_gocart_coupling == 1) ) then
                   !JJS 20110602 vvvvv
                   ! Conversion rate of cloud water to ice in the Bergeron porcess based on Meyer + DeMott formulae
                   !        
                   ! convert gocart aerosol mass conc to IN
                   !      p0 need to be converted from g*cm/s2/cm2 to mb (hPa)
                   !                call mass2icn(p0(i,j,k)*0.001,tair(i,j),aero(i,k,j,:), icn_out)
                   !                icn_cgs(i,j,k) = max(min_icn, icn_out) * 1.e-3   !IN conc [#/cm3] <-- [#/Litter]
                   !                icn_cgs(i,j,k) = icn_out * 1.e-3 !IN conc [#/cm3]      
                   !                if (rh_rad .ge. 100.) then
                   !                   print 1080, 'i', 'j', 'k', 'p0', 'tair', 'icn_cgs'
                   !                1080 format(3a4, 4a12)
                   !                   print 1081, i,j,k, p0(i,j,k)*0.001,tair(i,j), icn_cgs(i,j,k)
                   !                1081 format(3i4, 3e12.4) 
                   !                endif
                   !                r_nci = icn_cgs(i,j,k)  !DeMotto's formuale
                   r_nci = icn_out * 1.e-3  !DeMotto's formuale
                   !JJS 20110602 ^^^^^
                end if
#endif
                   
                dd(i,j)=(r00*qi(i,j)/r_nci)**y3(i,j)                  !meyers
                PIDW(i,j)=min(RR0*D2T*y2(i,j)*r_nci*dd(i,j), qc(i,j)) !meyers
             endif  !tairc
!<<< end of ifndef sat_predict
#endif
             pimm(i,j)=0.0
             pcfr(i,j)=0.0

             if (qc(i,j) .gt. 0.0) then
                y4(i,j) = 1./(tair(i,j)-c358)
                qsw(i,j)=rp0*exp(c172-c409*y4(i,j))
#ifdef new_saturation
                esw(i,j) = min(0.99*p0_mks(i,k,j),esw_mks(tair(i,j)))
                qsw(i,j) = 0.622*esw(i,j)/(p0_mks(i,k,j)-esw(i,j))
                esw(i,j) = esw(i,j)*10.  !in CGS
#endif
                xncld=qc(i,j)/4.e-9                         !cloud number
                esat=0.6112*exp(17.67*tairc(i,j)/(tairc(i,j)+243.5))*10.
                rv=0.622*esat/(p0(i,j,k)/1000.-esat)
                rlapse_m=980.616*(1.+2.5e6*rv/287./tair(i,j))/          &     
                       (1004.67+2.5e6*2.5e6*rv*0.622/(287.*tair(i,j)*tair(i,j)))
                delT=rlapse_m*0.5*(ww1(i,j,k)+ww1(i,j,k+1))
                if (delT.lt.0.) delT=0.
!  STEVE: PLAESE CHECK (2 4.e-9)
!               pimm(i,j)=xncld*Bhi*4.e-9*exp(-tairc(i,j))*delT*d2t*4.e-9
                pimm(i,j)=xncld*Bhi*4.e-9*exp(-tairc(i,j))*delT*d2t*4.e-9
 
                xccld=xncld*r00                           !cloud number concentration
                Xknud=7.37*tair(i,j)/(288.*Ra*p0(i,j,k))  ! Knudsen number
                alpha=1.257+0.400*exp(-1.10/Xknud)     ! Cunningham correction (P&Klett)
 
                cunnF=1.+alpha*Xknud                   ! Cunningham correction (P&Klett)

                if (tairc(i,j).ge.0.) then
                    dvair=(1.718+0.0049*tairc(i,j))*1.e-4   
                       !dynamic visc air (Prupp&Klett)
                else 
                    dvair=(1.718+0.0049*tairc(i,j)-1.2e-5*tairc(i,j)**2)*1.e-4  
                       !dynamic visc air (Prupp&Klett)
                endif
                DIFFar=1.3804e-16*tair(i,j)/6./cpi/dvair/Ra*cunnF     !aerosol diffusion via P&Klett
                if (qv(i,j)+qb0-qsw(i,j).lt.0.) then                  !only when cloud evaporating
                   pcfr(i,j)=4.e-9*4.*cpi*Rc*DIFFar*xccld*Cna*rr0*d2t    !Brownian part only via Cotton
                endif
             endif  !qc
          endif  !tair

!      STEVE: PLEASE CHECK
       y1(i,j)=pihom(i,j)+pidw(i,j)+pimm(i,j)+pcfr(i,j)-pimlt(i,j)

         if (y1(i,j).gt.qc(i,j)) then
            y1(i,j)=qc(i,j)
            y2(i,j)=1.
            y3(i,j)=pihom(i,j)+pidw(i,j)+pimm(i,j)+pcfr(i,j)
            if(y3(i,j).ne.0.) y2(i,j)=(qc(i,j)+pimlt(i,j))/y3(i,j)
            pihom(i,j)=pihom(i,j)*y2(i,j)
            pidw(i,j)=pidw(i,j)*y2(i,j)
            pimm(i,j)=pimm(i,j)*y2(i,j)
            pcfr(i,j)=pcfr(i,j)*y2(i,j)
         endif  !y1

#ifdef use_cpm
         cpm = cp + cvap*qv(i,j) + cliq*(qc(i,j)+qr(i,j)) + &
               cice*(qi(i,j)+qs(i,j)+qg(i,j))
         hlv = alv - (cliq-cvap)*(tair(i,j)-t0)
         hlf = alf - (cice-cliq)*(tair(i,j)-t0)
         hls = hlv + hlf
         avcp = hlv/cpm*pir
         ascp = hls/cpm*pir
         afcp = hlf/cpm*pir
#endif
         pt(i,j)=pt(i,j)+afcp*y1(i,j)
         qc(i,j)=qc(i,j)-y1(i,j)
         qi(i,j)=qi(i,j)+y1(i,j)

      else                                  ! below is 2007 saticel_s

#ifdef use_cpm
         cpm = cp + cvap*qv(i,j) + cliq*(qc(i,j)+qr(i,j)) + &
               cice*(qi(i,j)+qs(i,j)+qg(i,j))
         hlv = alv - (cliq-cvap)*(tair(i,j)-t0)
         hlf = alf - (cice-cliq)*(tair(i,j)-t0)
         hls = hlv + hlf
         avcp = hlv/cpm*pir
         ascp = hls/cpm*pir
         afcp = hlf/cpm*pir
#endif
         if (tair(i,j).lt.t0 .and. tair(i,j).gt.t00) then
            tairc(i,j)=tair(i,j)-t0
            y1(i,j)=max( min(tairc(i,j), -1.), -31.)
            it(i,j)=int(abs(y1(i,j)))
            y2(i,j)=aa1(it(i,j))
            y3(i,j)=aa2(it(i,j))
            y4(i,j)=exp(abs(beta*tairc(i,j)))
            y5(i,j)=(r00*qi(i,j)/(r25a*y4(i,j)))**y3(i,j)
            pidw(i,j)=min(r25rt*y2(i,j)*y4(i,j)*y5(i,j), qc(i,j))
         endif  !tair

         y1(i,j)=pihom(i,j)-pimlt(i,j)+pidw(i,j)
         pt(i,j)=pt(i,j)+afcp*y1(i,j)+ascp*(pidep(i,j))*d2t
         qv(i,j)=qv(i,j)-(pidep(i,j))*d2t
         qc(i,j)=qc(i,j)-y1(i,j)
         qi(i,j)=qi(i,j)+y1(i,j)

      endif  ! pocesses 24, 25 & 26 

#ifdef sat_predict
      ! --------------------------------------------------------------------------------
      ! saturation prediction from TCWA 1-moment scheme (Morrison and Milbrandt 2015) :
      !  (in MKS system)
      ! --------------------------------------------------------------------------------

      pact(i,j) = 0.0
      pint(i,j) = 0.0
      cnd(i,j) = 0.0
      dep(i,j) = 0.0
      fez(i,j) = 0.0
      ern(i,j) = 0.0

      ! -------------
      ! pact : cloud water activation
      ! -------------

      tair(i,j)  = (pt(i,j)+tb0)*pi0
      tairc(i,j) = tair(i,j)-t0
      rhoair = rho_mks(i,k,j)
#ifdef use_cpm
      cpm = cp + cvap*qv(i,j) + cliq*(qc(i,j)+qr(i,j)) + &
            cice*(qi(i,j)+qs(i,j)+qg(i,j))      ! specific heat capacity (in CGS)
      cpm1 = cpm*1.e-4                          ! specific heat capacity (in MKS)
      hlv = alv - (cliq-cvap)*(tair(i,j)-t0)    ! latent heat of vaporization (in CGS)
      xlv = hlv*1.e-4                           ! latent heat of vaporization (in MKS)
#else
      cpm1 = 1005.46*(1.+0.887*qv(i,j))         ! specific heat capacity (in MKS)
      xlv  = 3.1484E6-2370.*tair(i,j)           ! latent heat of vaporization (in MKS)
#endif
#ifdef new_saturation
      esw(i,j) = min(0.99*p0_mks(i,k,j),esw_mks(tair(i,j)))
      qsw(i,j) = 0.622*esw(i,j)/(p0_mks(i,k,j)-esw(i,j))
#else
      y1(i,j) = 1./(tair(i,j)-c358)
      qsw(i,j) = rp0*exp(c172-c409*y1(i,j))
#endif
      abw = 1.+xlv**2.*qsw(i,j)/cpm1/(4.61495E2*tair(i,j)**2.)  ! abw=1+dqsdT*(Lv/Cp)

      if ( qv(i,j).gt.qsw(i,j) .and. ww1(i,j,k).gt.1.e-3 ) then
        ! C1 and K1 over land and ocean (Roger and Yau)
        if ( xland(i,j) .eq. 1. ) then  !land
          C1 = 1000. ; K1 = 0.5
        else
          C1 = 120.  ; K1 = 0.4         !ocean
        endif
        ncloud = min(1.E9,1.E6*0.88*C1**(2./(K1+2.))*          &
                 (7.E-2*ww1(i,j,k)**1.5)**(K1/(K1+2.)))
        nact = max(0.,ncloud-qc(i,j)/5.236E-13)             ! CONVERT FROM CM-3 TO M-3 and 5 microm in radius
        qcmax = max((qv(i,j)-qsw(i,j)),0.)/abw
        pact(i,j) = min(max(nact*1.414E-14,0.),qcmax)       ! 1.5 micron in radius
        tair(i,j) = tair(i,j)+pact(i,j)*xlv/cpm1
        tairc(i,j) = tair(i,j)-t0
        qv(i,j) = max(0.,qv(i,j)-pact(i,j))
        qc(i,j) = max(0.,qc(i,j)+pact(i,j))
      endif

      ! -------------
      ! pint : cloud ice initialization
      ! -------------

#ifdef use_cpm
      cpm = cp + cvap*qv(i,j) + cliq*(qc(i,j)+qr(i,j)) + &
            cice*(qi(i,j)+qs(i,j)+qg(i,j))      ! specific heat capacity (in CGS)
      cpm1 = cpm*1.e-4                          ! specific heat capacity (in MKS)
      hlv = alv - (cliq-cvap)*(tair(i,j)-t0)
      hlf = alf - (cice-cliq)*(tair(i,j)-t0)
      hls = hlv + hlf                           ! latent heat of sublimation (in CGS)
      xls = hls*1.e-4                           ! latent heat of sublimation (in MKS)
#else
      cpm1  = 1005.46*(1.+0.887*qv(i,j))        ! specific heat capacity (in MKS)
      xls   = 3.15E6-2370.*tair(i,j)+0.3337E6   ! latent heat of sublimation  (in MKS)
#endif
#ifdef new_saturation
      esw(i,j) = min(0.99*p0_mks(i,k,j),esw_mks(tair(i,j)))
      esi(i,j) = min(0.99*p0_mks(i,k,j),esi_mks(tair(i,j)))
      if ( esi(i,j) .gt. esw(i,j) ) esi(i,j) = esw(i,j)
      qsw(i,j) = 0.622*esw(i,j)/(p0_mks(i,k,j)-esw(i,j))
      qsi(i,j) = 0.622*esi(i,j)/(p0_mks(i,k,j)-esi(i,j))
#else
      y1(i,j) = 1./(tair(i,j)-c358)
      qsw(i,j) = rp0*exp(c172-c409*y1(i,j))
      y2(i,j) = 1./(tair(i,j)-c76)
      qsi(i,j) = rp0*exp(c218-c580*y2(i,j))
#endif
      abi = 1.+xls**2.*qsi(i,j)/cpm1/(4.61495E2*tair(i,j)**2.)  ! abi=1+dqidT*(Ls/Cp)

      if ( qv(i,j).gt.qsi(i,j) .and. tair(i,j).lt.t0 ) then
        ssi(i,j) = qv(i,j)/qsi(i,j)-1.
        nice = 1.e3*exp(1.296E+1*ssi(i,j)-6.39E-1)
        r_nci = max(0.,nice-qi(i,j)/4.19E-10)               ! RHOI = 800; DI = 1.E-4
        qimax = max((qv(i,j)-qsi(i,j)),0.)/abi
        pint(i,j) = min(max(0.,r_nci*1.02e-13),qimax)
        tair(i,j) = tair(i,j)+pint(i,j)*xls/cpm1
        tairc(i,j) = tair(i,j)-t0
        qv(i,j) = max(0.,qv(i,j)-pint(i,j))
        qi(i,j) = max(0.,qi(i,j)+pint(i,j))
      endif

      ! -------------
      ! cnd : condensation of qv to qc
      ! dep : deposition of qv to qi
      ! fez : freezing of qc to qi
      ! ern : evaporation of qr
      ! -------------

#ifdef use_cpm
      cpm = cp + cvap*qv(i,j) + cliq*(qc(i,j)+qr(i,j)) + &
            cice*(qi(i,j)+qs(i,j)+qg(i,j))      ! specific heat capacity (in CGS)
      cpm1 = cpm*1.e-4                          ! specific heat capacity (in MKS)
      hlv = alv - (cliq-cvap)*(tair(i,j)-t0)    ! latent heat of vaporization (in CGS)
      xlv = hlv*1.e-4                           ! latent heat of vaporization (in MKS)
      hlf = alf - (cice-cliq)*(tair(i,j)-t0)    ! latent heat of fusion (in CGS)
      xlf = hlf*1.e-4                           ! latent heat of fusion (in MKS)
      hls = hlv + hlf                           ! latent heat of sublimation (in CGS)
      xls = hls*1.e-4                           ! latent heat of sublimation (in MKS)
#else
      cpm1  = 1005.46*(1.+0.887*qv(i,j))        ! specific heat capacity (in MKS)
      cpm   = 1.00546e+7*(1.+0.887*qv(i,j))     ! specific heat capacity (in CGS)
      xlv   = 3.1484E6-2370.*tair(i,j)          ! latent heat of vaporization (in MKS)
      xls   = 3.15E6-2370.*tair(i,j)+0.3337E6   ! latent heat of sublimation  (in MKS)
      xlf   = alf*1.e-4                         ! latent heat of fusion (in MKS)
#endif
#ifdef new_saturation
      esw(i,j) = min(0.99*p0_mks(i,k,j),esw_mks(tair(i,j)))
      esi(i,j) = min(0.99*p0_mks(i,k,j),esi_mks(tair(i,j)))
      if ( esi(i,j) .gt. esw(i,j) ) esi(i,j) = esw(i,j)
      qsw(i,j) = 0.622*esw(i,j)/(p0_mks(i,k,j)-esw(i,j))
      qsi(i,j) = 0.622*esi(i,j)/(p0_mks(i,k,j)-esi(i,j))
#else
      y1(i,j) = 1./(tair(i,j)-c358)
      qsw(i,j) = rp0*exp(c172-c409*y1(i,j))
      y2(i,j) = 1./(tair(i,j)-c76)
      qsi(i,j) = rp0*exp(c218-c580*y2(i,j))
#endif
      dv1 = 8.794E-4*pr0*tair(i,j)**1.81

      ! psychrometric correction to condensation/evaporation, abw=1+dqsdT*(Lv/Cp) ; dqsdT=Lv*qsw/Rv/T^2
      abw = 1.+xlv**2.*qsw(i,j)/cpm1/(4.61495E2*tair(i,j)**2.)
      ! psychrometric correction to deposition/sublimation,   abi=1+dqidT*(Ls/Cp) ; dqidT=Ls*qsi/Rv/T^2
      abi = 1.+xls**2.*qsi(i,j)/cpm1/(4.61495E2*tair(i,j)**2.)

      ! tauc : supersaturation relaxation timescale of cloud water
      if ( qc(i,j) .ge. cwmin ) then
        ltk  = log(tair(i,j))
        lqc  = -1.*log(qc(i,j)*rhoair)
        ltk2 = ltk*ltk
        lqc2 = lqc*lqc
        if ( xland(i,j) .eq. 1. ) then
          mvdc = exp(5.8936819 - 7.013013*ltk                  &
                     + 1.3178721*lqc + 1.1741987*ltk2          &
                     + 2.6110916E-3*lqc2 - 0.26646396*ltk*lqc)
        else
          mvdc = exp(173.57305 - 64.370929*ltk                 &
                     + 0.36833626*lqc + 6.1389254*ltk2         &
                     + 5.5915321E-3*lqc2 - 0.12488698*ltk*lqc)
        endif
        mvrc = max(5.e-7,min(5.e-5,mvdc/2.e+6))             ! diameter of cloud water
        ncloud = 3.*qc(i,j)*rhoair/(4.*cpi*1000.*mvrc**3.)  ! concentration of cloud water
        ncloud = min(1.e+9,max(0.1,ncloud))
        tauc = 1./(4.*cpi*dv1*ncloud*mvrc)                  ! tauc=1/(4*pi*Dv*Nc*rc)
      else
        tauc = 1.e+10
      endif

      ! taui : supersaturation relaxation timescale of cloud ice
      if ( qi(i,j) .ge. cimin ) then
        ltk  = log(tair(i,j))
        lqi  = -1.*log(qi(i,j)*rhoair)
        ltk2 = ltk*ltk
        lqi2 = lqi*lqi
        if ( xland(i,j) .eq. 1. ) then
          mvdi = exp(154.88767 - 0.25772005*lqi                &
                     + 3.75543E-4*lqi2 - 54.560357*ltk         &
                     + 5.1248879*ltk2 )
        else
          mvdi = exp(155.33396 - 0.25772005*lqi                &
                     + 3.75543E-4*lqi2 - 54.560357*ltk         &
                     + 5.1248879*ltk2 )
        endif
        mvri = min(5.e-4,max(3.e-6,mvdi/2.e+7))
        if ( tairc(i,j) .ge. -40.0 ) then
          hid = max(min(nint(abs(tairc(i,j))/0.25),120),0)
          inhgr = itble(hid)
          rhoi = 900.*exp(-3.*max((qv(i,j)-qsi(i,j))-5.e-5,0.) &
                 /inhgr)
        else
          inhgr = 5.
          ssi(i,j) = qv(i,j)/qsi(i,j)-1.
          if ( ssi(i,j) .gt. 0.267 ) then
            rhoi = -1027.456*ssi(i,j)+1185.834
          else
            rhoi = -32.332*ssi(i,j)+900.
          endif
        endif
        rhoi = min(max(rhoi,50.),900.)
        nice = 3.*qi(i,j)*rhoair/(4.*cpi*rhoi*mvri**3.)
        nice = min(1.e+7,max(0.1,nice))
        taui = 1./(4.*cpi*dv1*nice*mvri)                    ! taui=1/(4*pi*Dv*Ni*ri)
      else
        taui = 1.e+10
      endif

      ! taur : supersaturation relaxation timescale of rain
      if ( qr(i,j) .ge. crmin ) then
        ltk = log(tair(i,j))
        lqr = -1.*log(qr(i,j)*rhoair)
        ltk2 = ltk*ltk
        lqr2 = lqr*lqr
        kmin = 0.2222
        kmax = 0.8396
        mvdr = exp(33.999045 - 13.813047*ltk                   &
                   + 0.29834969*lqr + 1.6622025*ltk2           &
                   + 1.2480129E-3*lqr2 - 0.104641*ltk*lqr)
        efdr = exp(1.6855156 + 0.84147302*log(mvdr)            &
                   - 0.46783369*log(1.e+3)                     &
                   + 1.005305E-2*(log(mvdr))**2.               &
                   + 3.4833332E-2*(log(1.e+3))**2.             &
                   + 1.4727223E-2*log(mvdr)*log(1.e+3))
        efdr = max(min(efdr,mvdr/kmin**thrd),mvdr/kmax**thrd)
        kdxr = max(kmin,min(kmax,(mvdr/efdr)**3.))
        afar = max(0.,(6.*kdxr-3.+sqrt(8.*kdxr+1.))/(2.-2.*kdxr))
!        zr = (afar+3.)/efdr*1.e+6
        lzr = log((afar+3.)/efdr*1.e+6)
        tnr = log(6.*qr(i,j)*rhoair/cpi/1.e+3)                 & ! slope parameter for rain
              +(4.+afar)*lzr-lgamma(afar+4.)
        avr = exp(7.6004532 - 0.7990953*ltk                    & ! coefficient ... for rain (?)
                  + 1.0281818*lqr - 0.16595505*lqr2            &
                  + 1.110037E-2*lqr*lqr2                       &
                  - 2.0925743E-4*lqr2*lqr2)
        bvr = max(0.5,min(2., 1.2218728 - 0.1281004*ltk        & ! exponent ... for rain (?)
              + 2.6088596E-2*lqr - 7.4467639E-3*lqr2           &
              + 7.7592532E-4*lqr*lqr2                          &
              - 1.7056075E-5*lqr2*lqr2))
        mur = 1.496E-6*tair(i,j)**1.5/(tair(i,j)+120.)           ! shape parameter for rain
        rhoaj = sqrt(1.29/rhoair)
        gr2 = lgamma(afar+2.)
        gbr25 = lgamma(bvr*0.5+afar+2.5)
        taur = 1./(2.*cpi*dv1*(0.78*exp(tnr+gr2-(afar+2.)      &
               *lzr)+0.31*sqrt(avr*rhoaj/mur)*(mur/dv1)        &
               **thrd*exp(tnr+gbr25-(bvr*0.5+afar+2.5)         &
               *lzr)))
      else
        taur = 1.e+10
      endif

      cnd1 = 0. ; cnd2 = 0. ; dep1 = 0. ; dep2 = 0.
      fez1 = 0. ; fez2 = 0. ; ern1 = 0. ; ern2 = 0.
      latr = 0.

!>>> TGFS_modify
      ! increase condensation at lower-level tropics :
!      if ( abs(xlat) .le. 23.5 ) then
!        if ( p0(i,j,k).ge.9.1e+5 .and.                         &
!             p0(i,j,k).lt.9.4e+5 .and.                         &
!             ww1(i,j,k).gt.1.e-3 ) then
!          if ( qv(i,j) .gt. qsw(i,j) ) then
!            tauc = 1.
!          elseif ( qv(i,j) .lt. qsw(i,j)) then
!            tauc = 1.e+10
!          endif
!        elseif ( p0(i,j,k) .ge. 9.5e+5 ) then
!          if ( qv(i,j) .gt. qsw(i,j) ) then
!            tauc = 1.e+10
!          elseif ( qv(i,j) .lt. qsw(i,j) ) then
!            tauc = 1.
!          endif
!        endif
!      endif
!<<< TGFS_modify
 
      ! tau : multiphase supersaturation relaxation time scale defined by
      ! 1/tau = 1/tauc + 1/taur + (1+Ls/Cp*dqsl/dT)/(1+Ls/Cp*dqsi/dT)/taui
      tau = 1./(1./tauc+1./taur+1./taui*(1.+qsw(i,j)*          &
            xls*xlv/cpm1/(4.61495E+2*tair(i,j)**2.))/abi)

      ! atem : change in supersaturation due to Bergeron process
      atem = (qsi(i,j)-qsw(i,j))/taui*(1.+qsw(i,j)*xls*xlv     &
             /cpm1/(4.61495E+2*tair(i,j)**2.))/abi

      ! decide values at lower latitude :
      cnd1 = max(-qc(i,j),(atem*dt*tau/tauc+(qv(i,j)           &
             -qsw(i,j)-atem*tau)*(1.-exp(-dt/tau))*tau         &
             /tauc)/abw)
      dep1 = MAX(-qi(i,j),(atem*dt*tau/taui+(qv(i,j)           &
             -qsw(i,j)-atem*tau)*(1.-exp(-dt/tau))*tau         &
             /taui+(qsw(i,j)-qsi(i,j))*dt/taui)/abi)
      ern1 = MAX(-qr(i,j),(atem*dt*tau/taur+(qv(i,j)           &
             -qsw(i,j)-atem*tau)*(1.-exp(-dt/tau))*tau         &
             /taur)/abw)
      fez1 = 0.

#ifdef use_declination
      ! decide values at polar regions :
      ern2 = ern1
      if ( tair(i,j) .ge. 253.16 ) then
        ! T>=-20     : all ice melt; all excess vapor (over ice) condense; no deposition
        cnd2 = max(-qc(i,j),(qv(i,j)-qsi(i,j))/abi)
        fez2 = -qi(i,j)
        dep2 = 0.
      elseif ( tair(i,j).lt.253.16 .and. tair(i,j).ge.t00 ) then
        ! -40<=T<-20 : all water freeze; all excess vapor (over water) deposite; no condensation
        dep2 = max(-qi(i,j),(qv(i,j)-qsw(i,j))/abw)
        fez2 = qc(i,j)
        cnd2 = 0.
      else
        ! T<-40      : all water freeze; all excess vapor (over ice) deposite; no condensation
        dep2 = max(-qi(i,j),(qv(i,j)-qsi(i,j))/abi)
        fez2 = qc(i,j)
        cnd2 = 0.
      endif

      ! create transition zones :
      dltd = asin(sdec)*180.0/cpi      ! solar declination angle (in degree latitude)
!      latint = 23.5*2.0                ! width of the transition zone (in degree latitude)
      latint = 23.5                    ! width of the transition zone (in degree latitude)
      sbd = dltd-90.0                  ! southern boundary of sun (in degree latitude)
      nbd = dltd+90.0                  ! northern boundary of sun (in degree latitude)
      sbd1 = sbd+latint/2.0            ! southern boundary of southern transition zone
      sbd2 = sbd-latint/2.0            ! northern boundary of southern transition zone
      nbd1 = nbd-latint/2.0            ! southern boundary of northern transition zone
      nbd2 = nbd+latint/2.0            ! northern boundary of northern transition zone

      if ( xlat.gt.sbd1 .and. xlat.lt.nbd1 ) then       ! low-latitude
        fxlat = 0.
      elseif ( xlat.lt.sbd2 .or. xlat.gt.nbd2 ) then    ! polar region
        fxlat = 1.
      elseif ( xlat.ge.sbd2 .and. xlat.le.sbd1 ) then   ! southern transition zone
        fxlat = sin(abs(max((xlat-sbd1)*90.0/latint,-90.0)*cpi/180.0))
      elseif ( xlat.ge.nbd1 .and. xlat.le.nbd2 ) then   ! northern transition zone
        fxlat = sin(min((xlat-nbd1)*90.0/latint,90.0)*cpi/180.0)
      endif

      cnd(i,j) = cnd1*(1.0-fxlat) + cnd2*fxlat
      dep(i,j) = dep1*(1.0-fxlat) + dep2*fxlat
      fez(i,j) = fez1*(1.0-fxlat) + fez2*fxlat
      ern(i,j) = ern1*(1.0-fxlat) + ern2*fxlat

#else
!      ! (TGFS_modify) decide values at polar regions :
!      if ( abs(xlat) .ge. 66.0 ) then
!        ern2 = ern1
!        if ( tair(i,j) .ge. 253.16 ) then
!          ! T>=-20     : all ice melt; all excess vapor (over ice) condense; no deposition
!          cnd2 = max(-qc(i,j),(qv(i,j)-qsi(i,j))/abi)
!          fez2 = -qi(i,j)
!          dep2 = 0.
!        elseif ( tair(i,j).lt.253.16 .and. tair(i,j).ge.t00 ) then
!          ! -40<=T<-20 : all water freeze; all excess vapor (over water) deposite; no condensation
!          dep2 = max(-qi(i,j),(qv(i,j)-qsw(i,j))/abw)
!          fez2 = qc(i,j)
!          cnd2 = 0.
!        else
!          ! T<-40      : all water freeze; all excess vapor (over ice) deposite; no condensation
!          dep2 = max(-qi(i,j),(qv(i,j)-qsi(i,j))/abi)
!          fez2 = qc(i,j)
!          cnd2 = 0.
!        endif
!      endif
!
!      ! create a transition zone at 60~67 degree :
!      if ( abs(xlat) .lt. 66.0 ) then
!        cnd(i,j) = cnd1
!        dep(i,j) = dep1
!        fez(i,j) = 0.
!        ern(i,j) = ern1
!      elseif ( abs(xlat).ge.66.0 .and. abs(xlat).lt.67.0 ) then
!        latr = min(1.,max(0.,abs(xlat)-66.0))
!        cnd(i,j) = cnd1*(1.-latr)+cnd2*latr
!        dep(i,j) = dep1*(1.-latr)+dep2*latr
!        fez(i,j) = fez1*(1.-latr)+fez2*latr
!        ern(i,j) = ern1*(1.-latr)+ern2*latr
!      else
!        cnd(i,j) = cnd2
!        dep(i,j) = dep2
!        fez(i,j) = fez2
!        ern(i,j) = ern2
!      endif
      cnd(i,j) = cnd1
      dep(i,j) = dep1
      fez(i,j) = fez1
      ern(i,j) = ern1
#endif

      ! update tair, qv, qc, qi, qr :
!      tair(i,j) = tair(i,j) + (ern(i,j)+cnd(i,j))*xlv/cpm1     &
!                  + dep(i,j)*xls/cpm1 + fez(i,j)*alf/cpm
      tair(i,j) = tair(i,j) + (ern(i,j)+cnd(i,j))*xlv/cpm1     &
                  + dep(i,j)*xls/cpm1 + fez(i,j)*xlf/cpm1
      qv(i,j) = max(0.,qv(i,j)-cnd(i,j)-dep(i,j)-ern(i,j))
      qc(i,j) = max(0.,qc(i,j)+cnd(i,j)-fez(i,j))
      qi(i,j) = max(0.,qi(i,j)+dep(i,j)+fez(i,j))
      qr(i,j) = max(0.,qr(i,j)+ern(i,j))

!>>> TGFS_modify
      ! refreeze at mid-level polar regions :
!      if ( abs(xlat).ge.66.0 .and.                             &
!           p0(i,j,k).ge.8.e+5 .and.                            &
!           qc(i,j).ge.cwmin ) then
!        latr = min(1.,max(0.,abs(xlat)-66.0))
!        fez(i,j) = qc(i,j)*latr
!        tair(i,j) = tair(i,j) + fez(i,j)*alf/cpm
!        qi(i,j) = max(0.,qi(i,j)+fez(i,j))
!        qc(i,j) = max(0.,qc(i,j)-fez(i,j))
!      endif

      ! small cloud water all refreeze at polar regions, 
      ! in order to make the distribution more continuous :
!      if ( abs(xlat).ge.66.5 .and. qc(i,j).lt.1.e-4 ) then
!        tair(i,j) = tair(i,j)+qc(i,j)*alf/cpm
!        qi(i,j) = max(0.,qi(i,j)+qc(i,j))
!        qc(i,j) = 0.
!      endif
!>>> TGFS_modify

      ! update pt :
      pt(i,j) = tair(i,j)/pi0-tb0   !tair=(pt(i,j)+tb0)*pi0
      tairc(i,j) = tair(i,j)-t0

#else
!>>> if not define sat_predict, use saturation adjustment scheme :

!* 31 * pint  : initiation of qi                                  **31**
!* 32 * pidep : deposition of qi                                  **32**
!
!     CALCULATION OF PINT USES DIFFERENT VALUES OF THE INTERCEPT AND SLOPE FOR
!     THE FLETCHER EQUATION. ALSO, ONLY INITIATE MORE ICE IF THE NEW NUMBER
!     CONCENTRATION EXCEEDS THAT ALREADY PRESENT.

      if ( itaobraun .eq. 0 ) then            !  tao's original
         cn0=1.e-8
         beta=-.6
      else if ( itaobraun.eq.1 ) then         ! scott's
              cn0=1.e-6
              beta=-.46
      endif
        
      if (improve .eq. 3) then
         tair(i,j)=(pt(i,j)+tb0)*pi0
#ifdef use_cpm
         cpm = cp + cvap*qv(i,j) + cliq*(qc(i,j)+qr(i,j)) + &
               cice*(qi(i,j)+qs(i,j)+qg(i,j))
         hlv = alv - (cliq-cvap)*(tair(i,j)-t0)
         hlf = alf - (cice-cliq)*(tair(i,j)-t0)
         hls = hlv + hlf
         avcp = hlv/cpm*pir
         ascp = hls/cpm*pir
         afcp = hlf/cpm*pir
#endif
         if (tair(i,j) .lt. t0) then
            if (qi(i,j) .le. cmin) qi(i,j)=0.
!             if (qi(i,j) .le. cmin2) qi(i,j)=0.
            tairc(i,j)=tair(i,j)-t0
            rtair(i,j)=1./(tair(i,j)-c76)
            y2(i,j)=exp(c218-c580*rtair(i,j))
            qsi(i,j)=rp0*y2(i,j)
            esi(i,j)=c610*y2(i,j)
#ifdef new_saturation
            esi(i,j) = min(0.99*p0_mks(i,k,j),esi_mks(tair(i,j)))
            qsi(i,j) = 0.622*esi(i,j)/(p0_mks(i,k,j)-esi(i,j))
            esi(i,j) = esi(i,j)*10.  !in CGS
#endif
            ssi(i,j)=(qv(i,j)+qb0)/qsi(i,j)-1.
            y1(i,j)=1./tair(i,j)
            y3(i,j)=SQRT(qi(i,j))
            dd(i,j)=y1(i,j)*(RN10A*y1(i,j)-RN10B)+RN10C*tair(i,j)/   &
                                                    esi(i,j)
            dm(i,j)=max(qv(i,j)+qb0-qsi(i,j),0.0)
            rsub1(i,j)=cs580*qsi(i,j)*rtair(i,j)*rtair(i,j)
            dep(i,j)=dm(i,j)/(1.+rsub1(i,j))
            if (tairc(i,j).le.-5.) then     
                y4(i,j)=1./(tair(i,j)-c358)
                qsw(i,j)=rp0*exp(c172-c409*y4(i,j))
                fssi=min(xssi,max(.0,xssi*(tairc(i,j)+44.)/(44.-38.)))
                fssi=min(ssi(i,j),fssi)
!               r_nci=min(1.e-3*exp(-.639+12.96*fssi),1.) 
                r_nci=max(1.e-3*exp(-.639+12.96*fssi),0.528e-3)
       	        if (r_nci.gt.15.) r_nci=15.                         

#if ( WRF_CHEM == 1)
                ! EMK...Only execute if GOCART and coupling is on
                if ( (chem_opt == 300 .or. chem_opt == 301 .or. &
                      chem_opt == 302 .or. chem_opt == 303) .and. &
                      (gsfcgce_gocart_coupling == 1) ) then
                   !JJS 20110602 vvvvv
                   ! Conversion rate of cloud water to ice in the Bergeron porcess based on Meyer + DeMott formulae
                   !        
                   ! convert gocart aerosol mass conc to IN
                   !      p0 need to be converted from g*cm/s2/cm2 to mb (hPa)
                   !                call mass2icn(p0(i,j,k)*0.001,tair(i,j),aero(i,k,j,:), icn_out)
                   !                icn_cgs(i,j,k) = max(min_icn, icn_out) * 1.e-3   !IN conc [#/cm3] <-- [#/Litter]
                   !                call mass2icn(p0(i,j,k)*0.001,tair(i,j),aero(i,k,j,:), icn_out,i,j,k)
                   !                icn_out = min(1.e3, max(0.01e0 ,  icn_out) )
                   !                icn_cgs(i,j,k) = icn_out * 1.e-3 !IN conc [#/cm3]
                   !                r_nci = icn_cgs(i,j,k)  !DeMotto's formuale
                   r_nci = icn_out * 1.e-3  !DeMotto's formuale
                   !JJS 20110602 ^^^^^
                end if
#endif  
                   
                pidep(i,j)=max(R32RT*1.e4*fssi*sqrt(r_nci)*y3(i,j)/dd(i,j), 0.)
                dd(i,j)=max(1.e-9*r_nci/r00-qci(i,j,k)*1.e-9/ami50, 0.) 
                pint(i,j)=max(min(dd(i,j),dm(i,j)),0.)
                pint(i,j)=min(pint(i,j)+pidep(i,j), dep(i,j))
                if (pint(i,j).le.cmin) pint(i,j)=0.
                pt(i,j)=pt(i,j)+ascp*pint(i,j)
                qv(i,j)=qv(i,j)-pint(i,j)
                qi(i,j)=qi(i,j)+pint(i,j)
             endif  !taric
          endif  !tair                  
      endif  !improve    
!
          if ( itaobraun.eq.0 ) then
             tair(i,j)=(pt(i,j)+tb0)*pi0
#ifdef use_cpm
            cpm = cp + cvap*qv(i,j) + cliq*(qc(i,j)+qr(i,j)) + &
                  cice*(qi(i,j)+qs(i,j)+qg(i,j))
            hlv = alv - (cliq-cvap)*(tair(i,j)-t0)
            hlf = alf - (cice-cliq)*(tair(i,j)-t0)
            hls = hlv + hlf
            avcp = hlv/cpm*pir
            ascp = hls/cpm*pir
            afcp = hlf/cpm*pir
#endif
             if (tair(i,j) .lt. t0) then
                if (qi(i,j) .le. cmin1) qi(i,j)=0.
                tairc(i,j)=tair(i,j)-t0
                dd(i,j)=r31r*exp(beta*tairc(i,j))
                rtair(i,j)=1./(tair(i,j)-c76)
                y2(i,j)=exp(c218-c580*rtair(i,j))
                qsi(i,j)=rp0*y2(i,j)
                esi(i,j)=c610*y2(i,j)
#ifdef new_saturation
                esi(i,j) = min(0.99*p0_mks(i,k,j),esi_mks(tair(i,j)))
                qsi(i,j) = 0.622*esi(i,j)/(p0_mks(i,k,j)-esi(i,j))
                esi(i,j) = esi(i,j)*10.  !in CGS
#endif
                ssi(i,j)=(qv(i,j)+qb0)/qsi(i,j)-1.
                dm(i,j)=max( (qv(i,j)+qb0-qsi(i,j)), 0.)
                rsub1(i,j)=cs580*qsi(i,j)*rtair(i,j)*rtair(i,j)
                dep(i,j)=dm(i,j)/(1.+rsub1(i,j))
                pint(i,j)=max(min(dd(i,j), dm(i,j)), 0.)
                y1(i,j)=1./tair(i,j)
                y2(i,j)=exp(betah*tairc(i,j))
                y3(i,j)=sqrt(qi(i,j))
                dd(i,j)=y1(i,j)*(rn10a*y1(i,j)-rn10b) &
                     +rn10c*tair(i,j)/esi(i,j)
                pidep(i,j)=max(r32rt*ssi(i,j)*y2(i,j)*y3(i,j)/dd(i,j), 0.)
                pint(i,j)=pint(i,j)+pidep(i,j)

                pint(i,j)=min(pint(i,j),dep(i,j))
!c               if (pint(i,j) .le. cmin2) pint(i,j)=0.
                pt(i,j)=pt(i,j)+ascp*pint(i,j)
                qv(i,j)=qv(i,j)-pint(i,j)
                qi(i,j)=qi(i,j)+pint(i,j)
             endif  !tair
          endif  !if ( itaobraun.eq.0 )

          if ( itaobraun.eq.1 .and. improve .lt. 2) then
             tair(i,j)=(pt(i,j)+tb0)*pi0
#ifdef use_cpm
             cpm = cp + cvap*qv(i,j) + cliq*(qc(i,j)+qr(i,j)) + &
                   cice*(qi(i,j)+qs(i,j)+qg(i,j))
             hlv = alv - (cliq-cvap)*(tair(i,j)-t0)
             hlf = alf - (cice-cliq)*(tair(i,j)-t0)
             hls = hlv + hlf
             avcp = hlv/cpm*pir
             ascp = hls/cpm*pir
             afcp = hlf/cpm*pir
#endif
             if (tair(i,j) .lt. t0) then
                if (qi(i,j) .le. cmin) qi(i,j)=0.
                tairc(i,j)=tair(i,j)-t0
                rtair(i,j)=1./(tair(i,j)-c76)
                y2(i,j)=exp(c218-c580*rtair(i,j))
                qsi(i,j)=rp0*y2(i,j)
                esi(i,j)=c610*y2(i,j)
#ifdef new_saturation
                esi(i,j) = min(0.99*p0_mks(i,k,j),esi_mks(tair(i,j)))
                qsi(i,j) = 0.622*esi(i,j)/(p0_mks(i,k,j)-esi(i,j))
                esi(i,j) = esi(i,j)*10.  !in CGS
#endif
                ssi(i,j)=(qv(i,j)+qb0)/qsi(i,j)-1.
                ami20=3.76e-8
                y1(i,j)=1./tair(i,j)
                tairccri=tairc(i,j)          ! in degree c
                if (tairccri.le.-30.) tairccri=-30.
                y2(i,j)=exp(betah*tairccri)
                y3(i,j)=sqrt(qi(i,j))
                dd(i,j)=y1(i,j)*(rn10a*y1(i,j)-rn10b)+rn10c*tair(i,j)/esi(i,j)
                pidep(i,j)=max(10.*r32rt*ssi(i,j)*y2(i,j)*y3(i,j)/dd(i,j),0.e0)  !fixed SEL
                r_nci=min(rn25*exp(beta*tairc(i,j)),1.)
                dd(i,j)=max(1.e-9*r_nci/r00-qi(i,j)*1.e-9/ami20, 0.)
                dm(i,j)=max( (qv(i,j)+qb0-qsi(i,j)), 0.0)
                rsub1(i,j)=cs580*qsi(i,j)*rtair(i,j)*rtair(i,j)
                dep(i,j)=dm(i,j)/(1.+rsub1(i,j))
                pint(i,j)=max(min(dd(i,j), dm(i,j)), 0.)
                pint(i,j)=min(pint(i,j)+pidep(i,j), dep(i,j))
                if (pint(i,j) .le. cmin) pint(i,j)=0.
                pt(i,j)=pt(i,j)+ascp*pint(i,j)
                qv(i,j)=qv(i,j)-pint(i,j)
                qi(i,j)=qi(i,j)+pint(i,j)
             endif  !tair
          endif  ! if ( itaobraun.eq.1 ) ! scott's

! End of Process 31 & 32

! WRF satice has new_ice_sat option 0, 1 and 2
! Steve's satice has new_ice_sat option 0, 1, 2, 3 and 9
! option 0, 1 and 2 are identical in both satice
! I added option 3 below and wrapped them with "if (improve.eq.3)"

!*****   TAO ET AL (1989) SATURATION TECHNIQUE  ***********************

         if (new_ice_sat .eq. 0) then

               tair(i,j)=(pt(i,j)+tb0)*pi0
#ifdef use_cpm
               cpm = cp + cvap*qv(i,j) + cliq*(qc(i,j)+qr(i,j)) + &
                     cice*(qi(i,j)+qs(i,j)+qg(i,j))
               hlv = alv - (cliq-cvap)*(tair(i,j)-t0)
               hlf = alf - (cice-cliq)*(tair(i,j)-t0)
               hls = hlv + hlf
               avcp = hlv/cpm*pir
               ascp = hls/cpm*pir
               afcp = hlf/cpm*pir
#endif
               cnd(i,j)=rt0*(tair(i,j)-t00)
               dep(i,j)=rt0*(t0-tair(i,j))
               y1(i,j)=1./(tair(i,j)-c358)
               y2(i,j)=1./(tair(i,j)-c76)
               qsw(i,j)=rp0*exp(c172-c409*y1(i,j))
               qsi(i,j)=rp0*exp(c218-c580*y2(i,j))
#ifdef new_saturation
               esw(i,j) = min(0.99*p0_mks(i,k,j),esw_mks(tair(i,j)))
               esi(i,j) = min(0.99*p0_mks(i,k,j),esi_mks(tair(i,j)))
               if ( esi(i,j) .gt. esw(i,j) ) esi(i,j) = esw(i,j)
               qsw(i,j) = 0.622*esw(i,j)/(p0_mks(i,k,j)-esw(i,j))
               qsi(i,j) = 0.622*esi(i,j)/(p0_mks(i,k,j)-esi(i,j))
               esw(i,j) = esw(i,j)*10.  !in CGS
               esi(i,j) = esi(i,j)*10.  !in CGS
#endif
               dd(i,j)=cp409*y1(i,j)*y1(i,j)
               dd1(i,j)=cp580*y2(i,j)*y2(i,j)
               if (qc(i,j).le.cmin) qc(i,j)=cmin
               if (qi(i,j).le.cmin) qi(i,j)=cmin
               if (tair(i,j).ge.t0) then
                  dep(i,j)=0.0
                  cnd(i,j)=1.
                  qi(i,j)=0.0
               endif

               if (tair(i,j).lt.t00) then
                  cnd(i,j)=0.0
                  dep(i,j)=1.
                  qc(i,j)=0.0
               endif

               y5(i,j)=avcp*cnd(i,j)+ascp*dep(i,j)
!               if (qc(i,j) .ge. cmin .or. qi(i,j) .ge. cmin) then
               y1(i,j)=qc(i,j)*qsw(i,j)/(qc(i,j)+qi(i,j))
               y2(i,j)=qi(i,j)*qsi(i,j)/(qc(i,j)+qi(i,j))
               y4(i,j)=dd(i,j)*y1(i,j)+dd1(i,j)*y2(i,j)
               qvs(i,j)=y1(i,j)+y2(i,j)
               rsub1(i,j)=(qv(i,j)+qb0-qvs(i,j))/(1.+y4(i,j)*y5(i,j))
               cnd(i,j)=cnd(i,j)*rsub1(i,j)
               dep(i,j)=dep(i,j)*rsub1(i,j)
               if (qc(i,j).le.cmin) qc(i,j)=0.
               if (qi(i,j).le.cmin) qi(i,j)=0.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    ******   condensation or evaporation of qc  ******

               cnd(i,j)=max(-qc(i,j),cnd(i,j))

!c    ******   deposition or sublimation of qi    ******

               dep(i,j)=max(-qi(i,j),dep(i,j))

               pt(i,j)=pt(i,j)+avcp*cnd(i,j)+ascp*dep(i,j)
               qv(i,j)=qv(i,j)-cnd(i,j)-dep(i,j)
               qc(i,j)=qc(i,j)+cnd(i,j)
               qi(i,j)=qi(i,j)+dep(i,j)

         endif  !if (new_ice_sat .eq. 0)

         if (new_ice_sat .eq. 1) then

               tair(i,j)=(pt(i,j)+tb0)*pi0
#ifdef use_cpm
               cpm = cp + cvap*qv(i,j) + cliq*(qc(i,j)+qr(i,j)) + &
                     cice*(qi(i,j)+qs(i,j)+qg(i,j))
               hlv = alv - (cliq-cvap)*(tair(i,j)-t0)
               hlf = alf - (cice-cliq)*(tair(i,j)-t0)
               hls = hlv + hlf
               avcp = hlv/cpm*pir
               ascp = hls/cpm*pir
               afcp = hlf/cpm*pir
#endif
               cnd(i,j)=rt0*(tair(i,j)-t00)
               dep(i,j)=rt0*(t0-tair(i,j))
               y1(i,j)=1./(tair(i,j)-c358)
               y2(i,j)=1./(tair(i,j)-c76)
               qsw(i,j)=rp0*exp(c172-c409*y1(i,j))
               qsi(i,j)=rp0*exp(c218-c580*y2(i,j))
#ifdef new_saturation
               esw(i,j) = min(0.99*p0_mks(i,k,j),esw_mks(tair(i,j)))
               esi(i,j) = min(0.99*p0_mks(i,k,j),esi_mks(tair(i,j)))
               if ( esi(i,j) .gt. esw(i,j) ) esi(i,j) = esw(i,j)
               qsw(i,j) = 0.622*esw(i,j)/(p0_mks(i,k,j)-esw(i,j))
               qsi(i,j) = 0.622*esi(i,j)/(p0_mks(i,k,j)-esi(i,j))
               esw(i,j) = esw(i,j)*10.  !in CGS
               esi(i,j) = esi(i,j)*10.  !in CGS
#endif
               dd(i,j)=cp409*y1(i,j)*y1(i,j)
               dd1(i,j)=cp580*y2(i,j)*y2(i,j)
               y5(i,j)=avcp*cnd(i,j)+ascp*dep(i,j)
               y1(i,j)=rt0*(tair(i,j)-t00)*qsw(i,j)
               y2(i,j)=rt0*(t0-tair(i,j))*qsi(i,j)
!               IF (QC(I,J).LE.CMIN) QC(I,J)=CMIN
!               IF (QI(I,J).LE.CMIN) QI(I,J)=CMIN

               if (tair(i,j).ge.t0) then
!                 QI(I,J)=0.0
                  dep(i,j)=0.0
                  cnd(i,j)=1.
                  y2(i,j)=0.

                  y1(i,j)=qsw(i,j)
               endif
               if (tair(i,j).lt.t00) then
                  cnd(i,j)=0.0
                  dep(i,j)=1.
                  y2(i,j)=qsi(i,j)
                  y1(i,j)=0.
!                 QC(I,J)=0.0
               endif

!            Y1(I,J)=QC(I,J)*QSW(I,J)/(QC(I,J)+QI(I,J))
!            Y2(I,J)=QI(I,J)*QSI(I,J)/(QC(I,J)+QI(I,J))

               y4(i,j)=dd(i,j)*y1(i,j)+dd1(i,j)*y2(i,j)
               qvs(i,j)=y1(i,j)+y2(i,j)
               rsub1(i,j)=(qv(i,j)+qb0-qvs(i,j))/(1.+y4(i,j)*y5(i,j))
               cnd(i,j)=cnd(i,j)*rsub1(i,j)
               dep(i,j)=dep(i,j)*rsub1(i,j)
!               IF (QC(I,J).LE.CMIN) QC(I,J)=0.
!               IF (QI(I,J).LE.CMIN) QI(I,J)=0.

!C    ******   CONDENSATION OR EVAPORATION OF QC  ******

               cnd(i,j)=max(-qc(i,j),cnd(i,j))

!C    ******   DEPOSITION OR SUBLIMATION OF QI    ******

               dep(i,j)=max(-qi(i,j),dep(i,j))

               pt(i,j)=pt(i,j)+avcp*cnd(i,j)+ascp*dep(i,j)
               qv(i,j)=qv(i,j)-cnd(i,j)-dep(i,j)
               qc(i,j)=qc(i,j)+cnd(i,j)
               qi(i,j)=qi(i,j)+dep(i,j)

         endif  ! if (new_ice_sat .eq. 1)

!c
!
          if (new_ice_sat .eq. 2) then

          dep(i,j)=0.0
          cnd(i,j)=0.0
          tair(i,j)=(pt(i,j)+tb0)*pi0
#ifdef use_cpm
          cpm = cp + cvap*qv(i,j) + cliq*(qc(i,j)+qr(i,j)) + &
                cice*(qi(i,j)+qs(i,j)+qg(i,j))
          hlv = alv - (cliq-cvap)*(tair(i,j)-t0)
          hlf = alf - (cice-cliq)*(tair(i,j)-t0)
          hls = hlv + hlf
          avcp = hlv/cpm*pir
          ascp = hls/cpm*pir
          afcp = hlf/cpm*pir
#endif
          if (tair(i,j) .ge. 253.16) then
              y1(i,j)=1./(tair(i,j)-c358)
              qsw(i,j)=rp0*exp(c172-c409*y1(i,j))
#ifdef new_saturation
              esw(i,j) = min(0.99*p0_mks(i,k,j),esw_mks(tair(i,j)))
              qsw(i,j) = 0.622*esw(i,j)/(p0_mks(i,k,j)-esw(i,j))
              esw(i,j) = esw(i,j)*10.  !in CGS
#endif
              dd(i,j)=cp409*y1(i,j)*y1(i,j)
              dm(i,j)=qv(i,j)+qb0-qsw(i,j)
              cnd(i,j)=dm(i,j)/(1.+avcp*dd(i,j)*qsw(i,j))
!c    ******   condensation or evaporation of qc  ******
              cnd(i,j)=max(-qc(i,j), cnd(i,j))
             pt(i,j)=pt(i,j)+avcp*cnd(i,j)
             qv(i,j)=qv(i,j)-cnd(i,j)
             qc(i,j)=qc(i,j)+cnd(i,j)
         endif
          if (tair(i,j) .le. 258.16) then
!c             cnd(i,j)=0.0
             y2(i,j)=1./(tair(i,j)-c76)
             qsi(i,j)=rp0*exp(c218-c580*y2(i,j))
#ifdef new_saturation
             esi(i,j) = min(0.99*p0_mks(i,k,j),esi_mks(tair(i,j)))
             qsi(i,j) = 0.622*esi(i,j)/(p0_mks(i,k,j)-esi(i,j))
             esi(i,j) = esi(i,j)*10.  !in CGS
#endif
             dd1(i,j)=cp580*y2(i,j)*y2(i,j)
             dep(i,j)=(qv(i,j)+qb0-qsi(i,j))/(1.+ascp*dd1(i,j)*qsi(i,j))
!c    ******   deposition or sublimation of qi    ******
             dep(i,j)=max(-qi(i,j),dep(i,j))
             pt(i,j)=pt(i,j)+ascp*dep(i,j)
             qv(i,j)=qv(i,j)-dep(i,j)
             qi(i,j)=qi(i,j)+dep(i,j)
          endif

      endif                              ! if (new_ice_sat .eq. 2)

!!!!
!!!   new_ice_sat option 3 from  Steve's satice
!!!!

      if (new_ice_sat .eq. 3) then

         dep(i,j)=0.0
         cnd(i,j)=0.0
         tair(i,j)=(pt(i,j)+tb0)*pi0
#ifdef use_cpm
         cpm = cp + cvap*qv(i,j) + cliq*(qc(i,j)+qr(i,j)) + &
               cice*(qi(i,j)+qs(i,j)+qg(i,j))
         hlv = alv - (cliq-cvap)*(tair(i,j)-t0)
         hlf = alf - (cice-cliq)*(tair(i,j)-t0)
         hls = hlv + hlf
         avcp = hlv/cpm*pir
         ascp = hls/cpm*pir
         afcp = hlf/cpm*pir
#endif
         if (tair(i,j).ge.t00) THEN
            y1(i,j)=1./(tair(i,j)-c358)
            qsw(i,j)=rp0*exp(c172-c409*y1(i,j))
#ifdef new_saturation
            esw(i,j) = min(0.99*p0_mks(i,k,j),esw_mks(tair(i,j)))
            qsw(i,j) = 0.622*esw(i,j)/(p0_mks(i,k,j)-esw(i,j))
            esw(i,j) = esw(i,j)*10.  !in CGS
#endif
            dd(i,j)=cp409*y1(i,j)*y1(i,j)
            dm(i,j)=qv(i,j)+qb0-qsw(i,j)
            cnd(i,j)=dm(i,j)/(1.+avcp*dd(i,j)*qsw(i,j))
!c    ******   condensation or evaporation of qc  ******
            cnd(i,j)=max(-qc(i,j),cnd(i,j))
            pt(i,j)=pt(i,j)+avcp*cnd(i,j)
            qv(i,j)=qv(i,j)-cnd(i,j)
            qc(i,j)=qc(i,j)+cnd(i,j)
         endif
         tair(i,j)=(pt(i,j)+tb0)*pi0
#ifdef use_cpm
         cpm = cp + cvap*qv(i,j) + cliq*(qc(i,j)+qr(i,j)) + &
               cice*(qi(i,j)+qs(i,j)+qg(i,j))
         hlv = alv - (cliq-cvap)*(tair(i,j)-t0)
         hlf = alf - (cice-cliq)*(tair(i,j)-t0)
         hls = hlv + hlf
         avcp = hlv/cpm*pir
         ascp = hls/cpm*pir
         afcp = hlf/cpm*pir
#endif
         if (tair(i,j).le.273.16) THEN
!c    ******   deposition or sublimation of qi    ******
            y1(i,j)=1./(tair(i,j)-c358)
            qsw(i,j)=rp0*exp(c172-c409*y1(i,j))
            y2(i,j)=1./(tair(i,j)-c76)
            qsi(i,j)=rp0*exp(c218-c580*y2(i,j))
#ifdef new_saturation
            esw(i,j) = min(0.99*p0_mks(i,k,j),esw_mks(tair(i,j)))
            esi(i,j) = min(0.99*p0_mks(i,k,j),esi_mks(tair(i,j)))
            if ( esi(i,j) .gt. esw(i,j) ) esi(i,j) = esw(i,j)
            qsw(i,j) = 0.622*esw(i,j)/(p0_mks(i,k,j)-esw(i,j))
            qsi(i,j) = 0.622*esi(i,j)/(p0_mks(i,k,j)-esi(i,j))
            esw(i,j) = esw(i,j)*10.  !in CGS
            esi(i,j) = esi(i,j)*10.  !in CGS
#endif

!            fssi=min(0.20,max(0.,0.20*(tair(i,j)-t0+44.0)/(44.0-38.0)))
!vvvvvvvvvvvvv Tao 20110722 vvvvvvvvvvvvvv
            fssi=min(xssi,max(0.,xssi*(tair(i,j)-t0+44.0)/(44.0-38.0)))
!^^^^^^^^^^^^^ Tao 20110722 ^^^^^^^^^^^^^^

            ! max ssi f(tair)  ^^^^^
            y3(i,j)=1.+min((qsw(i,j)-qsi(i,j))/qsi(i,j), fssi)
            ! Lang 2007  ^^^^^
            y4(i,j)=qsi(i,j)*y3(i,j)
            if (tair(i,j).le.268.16.and.(qv(i,j)+qb0.gt.y4(i,j))) then

               dd1(i,j)=cp580*y2(i,j)*y2(i,j)
               dep(i,j)=(qv(i,j)+qb0-y4(i,j))/(1.+ascp*dd1(i,j)*y4(i,j))
            else if (qv(i,j)+qb0.lt.qsi(i,j).and.qi(i,j).gt.cmin) then
               dd1(i,j)=cp580*y2(i,j)*y2(i,j)
               dep(i,j)=(qv(i,j)+qb0-qsi(i,j))/(1.+ascp*dd1(i,j)*qsi(i,j))
               dep(i,j)=max(-qi(i,j),dep(i,j))
            endif
            pt(i,j)=pt(i,j)+ascp*dep(i,j)
            qv(i,j)=qv(i,j)-dep(i,j)
            qi(i,j)=qi(i,j)+dep(i,j)
         endif
    
    endif                                        ! if (new_ice_sat = 3)

!<<< end of ifdef sat_predict
#endif

!* 10 * PSDEP : DEPOSITION OR SUBLIMATION OF QS                   **10**
!* 20 * PGSUB : SUBLIMATION OF QG                                 **20**

        psdep(i,j)=0.0
        pgdep(i,j)=0.0
        pssub(i,j)=0.0
        pgsub(i,j)=0.0
        tair(i,j)=(pt(i,j)+tb0)*pi0
        tairc(I,j)=tair(I,j)-t0

!        if (improve.gt.2) call sgmap(1,qs(i,j),r00,tairc(i,j),ftns0(i,j))
!        if (improve.gt.2) call sgmap(2,qg(i,j),r00,tairc(i,j),ftng0(i,j))

        if (qs(i,j).lt.cmin1) qs(i,j)=0.0
        if (qg(i,j).lt.cmin1) qg(i,j)=0.0
        if (qc(i,j)+qi(i,j).gt.1.e-5) then
           dlt1(i,j)=1.
        else    
           dlt1(i,j)=0.
        endif

!        IF (IMPROVE .NE. 3) THEN ! below is 2007 WRF saticel_s
!vvvvvvvvvvvv Tao 20110722 vvvvvvvvvvvvvvvvvvvv
        IF (IMPROVE .NE. 3 .and. improve .ne. -1) THEN ! below is 2007 WRF saticel_s
!^^^^^^^^^^^^ Tao 20110722 ^^^^^^^^^^^^^^^^^^^^

           if (tair(i,j).lt.t0) then
              rtair(i,j)=1./(tair(i,j)-c76)
              qsi(i,j)=rp0*exp(c218-c580*rtair(i,j))
              ssi(i,j)=(qv(i,j)+qb0)/qsi(i,j)-1.
!
              y1(i,j)=r10ar/(tca(i,j)*tair(i,j)**2)+1./(dwv(i,j) &
                      *qsi(i,j))
              y2(i,j)=.78/zs(i,j)**2+r101f*scv(i,j)/zs(i,j)**bsh5
              psdep(i,j)=r10t*ssi(i,j)*y2(i,j)/y1(i,j)
              pssub(i,j)=psdep(i,j)
              psdep(i,j)=r2is*max(psdep(i,j), 0.)
              pssub(i,j)=r2is*max(-qs(i,j), min(pssub(i,j), 0.))

              if (ihail.eq.1) then
                 y2(i,j)=.78/zg(i,j)**2+r20bq*scv(i,j)/zg(i,j)**bgh5
              else
                 y2(i,j)=.78/zg(i,j)**2+r20bs*scv(i,j)/zg(i,j)**bgh5
              endif

              pgsub(i,j)=r2ig*r20t*ssi(i,j)*y2(i,j)/y1(i,j)
              dm(i,j)=qv(i,j)+qb0-qsi(i,j)
              rsub1(i,j)=cs580*qsi(i,j)*rtair(i,j)*rtair(i,j)

!     ********   DEPOSITION OR SUBLIMATION OF QS  **********************

              y1(i,j)=dm(i,j)/(1.+rsub1(i,j))
              psdep(i,j)=r2is*min(psdep(i,j),max(y1(i,j),0.))
              y2(i,j)=min(y1(i,j),0.)
              pssub(i,j)=r2is*max(pssub(i,j),y2(i,j))

!     ********   SUBLIMATION OF QG   ***********************************

              dd(i,j)=max((-y2(i,j)-qs(i,j)), 0.)
              pgsub(i,j)=r2ig*min(dd(i,j), qg(i,j), max(pgsub(i,j),0.))

              psdep(i,j)=dlt1(i,j)*psdep(i,j)
              pssub(i,j)=(1.-dlt1(i,j))*pssub(i,j)
              pgsub(i,j)=(1.-dlt1(i,j))*pgsub(i,j)

              pt(i,j)=pt(i,j)+ascp*(psdep(i,j)+pssub(i,j)-pgsub(i,j))
              qv(i,j)=qv(i,j)+pgsub(i,j)-pssub(i,j)-psdep(i,j)
              qs(i,j)=qs(i,j)+psdep(i,j)+pssub(i,j)
              qg(i,j)=qg(i,j)-pgsub(i,j)

           endif ! tair

        else   ! IF (IMPROVE .NE. 3)

           if (tair(i,j) .lt. t0) then 

#ifdef use_cpm
              cpm = cp + cvap*qv(i,j) + cliq*(qc(i,j)+qr(i,j)) + &
                    cice*(qi(i,j)+qs(i,j)+qg(i,j))
              hlv = alv - (cliq-cvap)*(tair(i,j)-t0)
              hlf = alf - (cice-cliq)*(tair(i,j)-t0)
              hls = hlv + hlf
              avcp = hlv/cpm*pir
              ascp = hls/cpm*pir
              afcp = hlf/cpm*pir
#endif
              rtair(i,j)=1./(tair(i,j)-c76)
              y2(i,j)=exp(c218-c580*rtair(i,j))
              qsi(i,j)=rp0*y2(i,j)
              esi(i,j)=c610*y2(i,j)
#ifdef new_saturation
              esi(i,j)=min(0.99*p0_mks(i,k,j),esi_mks(tair(i,j)))
              qsi(i,j)=0.622*esi(i,j)/(p0_mks(i,k,j)-esi(i,j))
              esi(i,j) = esi(i,j)*10.  !in CGS
#endif

              SSI(I,J)=(QV(I,J)+QB0)/QSI(I,J)-1.
              IF (DLT1(I,J).EQ.1.) SSI(I,J)=max(SSI(I,J),0.)
              IF (DLT1(I,J).EQ.0.) SSI(I,J)=min(SSI(I,J),0.)
              DM(I,J)=QV(I,J)+QB0-QSI(I,J)
              RSUB1(I,J)=CS580*QSI(I,J)*RTAIR(I,J)*RTAIR(I,J)
              DD1(I,J)=DM(I,J)/(1.+RSUB1(I,J))
              Y3(I,J)=1./TAIR(I,J)
              DD(I,J)=Y3(I,J)*(RN10A*Y3(I,J)-RN10B)+RN10C*TAIR(I,J)/ESI(I,J)
              TAIRC(I,J)=TAIR(I,J)-T0
!
            ftns(i,j)=1.
            ftng(i,j)=1.
            ftns0(i,j)=1.
            ftng0(i,j)=1.
!            IF (IMPROVE .GE. 3) THEN
            call sgmap(1,qs(i,j),r00,tairc(i,j),ftns0(i,j))
            call sgmap(2,qg(i,j),r00,tairc(i,j),ftng0(i,j))
            ftns(i,j)=ftns0(i,j)
            ftng(i,j)=ftng0(i,j)
!            ENDIF
!
              Y4(I,J)=R10T*SSI(I,J)*(R101R/ZS(I,J)**2+R102RF/ZS(I,J)**BSH5)   &
                                /DD(I,J)*ftns(i,j)
              PSDEP(I,J)=max(-QS(I,J), Y4(I,J))
              DD(I,J)=Y3(I,J)*(RN20A*Y3(I,J)-RN20B)+RN10C*TAIR(I,J)/ESI(I,J)
              Y2(I,J)=R191R/ZG(I,J)**2+R192RF/ZG(I,J)**BGH5
              PGDEP(I,J)=MAX(-qg(i,j), R20T*SSI(I,J)*Y2(I,J)/DD(I,J)          &
                                     *ftng(i,j))
!      ******************************************************************
              Y5(I,J)=min(0.,DD1(I,J))
              DD1(I,J)=max(0.,DD1(I,J))
              IF (DLT1(I,J).EQ.1.) THEN
                 Y1(I,J)=MIN(PSDEP(I,J)+PGDEP(I,J), DD1(I,J) )
                 IF (PSDEP(I,J).ge.DD1(I,J)) THEN
                    PSDEP(I,J)=DD1(I,J)
                    PGDEP(I,J)=0.
                 ENDIF
                 IF (DD1(I,J).gt.PSDEP(I,J).and.(PSDEP(I,J)+PGDEP(I,J)).gt.  &
                    DD1(I,J)) PGDEP(I,J)=Y1(I,J)-PSDEP(I,J)
              ENDIF
              IF (DLT1(I,J).EQ.0.) THEN
                 Y1(I,J)=MAX(PSDEP(I,J)+PGDEP(I,J), Y5(I,J) )
                 IF (Y5(I,J).gt.(PSDEP(I,J)+PGDEP(I,J))) THEN
                    Y3(I,J)=(PSDEP(I,J)+PGDEP(I,J))
                    IF (Y3(I,J).ne.0.0) THEN
                       PSDEP(I,J)=PSDEP(I,J)/Y3(I,J)*Y5(I,J)
                       PGDEP(I,J)=PGDEP(I,J)/Y3(I,J)*Y5(I,J)
                    ENDIF
                 ENDIF
              ENDIF

              PSSUB(i,j)=r2is*min(PSDEP(i,j), 0.)
              PSDEP(i,j)=r2is*max(PSDEP(i,j), 0.)
              PGSUB(i,j)=r2ig*min(PGDEP(i,j), 0.)
              PGDEP(i,j)=r2ig*max(PGDEP(i,j), 0.)

              pt(i,j)=pt(i,j)+ascp*y1(i,j)
              qv(i,j)=qv(i,j)-y1(i,j)
              qs(i,j)=qs(i,j)+psdep(i,j)+pssub(i,j)
              qg(i,j)=qg(i,j)+pgdep(i,j)+pgsub(i,j)

           endif   ! if (tair(i,j) .lt. t0)

        endif !  if (improve.ne.3)

!!!  end of Processes 10 and 20

#ifndef sat_predict
!>>> Note that ern is already calculated in saturation prediction scheme.

!* 23 * ERN : EVAPORATION OF QR (SUBSATURATION)                   **23**
! Steve did not make any improvement on this process

       if (qr(i,j) .gt. 0.0) then
          tair(i,j)=(pt(i,j)+tb0)*pi0
          rtair(i,j)=1./(tair(i,j)-c358)
#ifdef use_cpm
          cpm = cp + cvap*qv(i,j) + cliq*(qc(i,j)+qr(i,j)) + &
                cice*(qi(i,j)+qs(i,j)+qg(i,j))
          hlv = alv - (cliq-cvap)*(tair(i,j)-t0)
          hlf = alf - (cice-cliq)*(tair(i,j)-t0)
          hls = hlv + hlf
          avcp = hlv/cpm*pir
          ascp = hls/cpm*pir
          afcp = hlf/cpm*pir
#endif
          if (improve.eq.3) then
	     y2(i,j)=exp( c172-c409*rtair(i,j) )
	     esw(i,j)=c610*y2(i,j)
             qsw(i,j)=rp0*y2(i,j)
#ifdef new_saturation
             esw(i,j)=min(0.99*p0_mks(i,k,j),esw_mks(tair(i,j)))
             qsw(i,j)=0.622*esw(i,j)/(p0_mks(i,k,j)-esw(i,j))
             esw(i,j) = esw(i,j)*10.  !in CGS
#endif
             ssw(i,j)=(qv(i,j)+qb0)/qsw(i,j)-1.
             dm(i,j)=qv(i,j)+qb0-qsw(i,j)
             rsub1(i,j)=cv409*qsw(i,j)*rtair(i,j)*rtair(i,j)
             dd1(i,j)=max(-dm(i,j)/(1.+rsub1(i,j)),0.0)
             y3(i,j)=1./tair(i,j)
             dd(i,j)=y3(i,j)*(rn30a*y3(i,j)-rn10b)+rn10c*tair(i,j)  &
                             /esw(i,j)
             y1(i,j)=-r23t*ssw(i,j)*(r231r/zr(i,j)**2+r232rf/       &
                                      zr(i,j)**3)/dd(i,j)
             ern(i,j)=min(dd1(i,j),qr(i,j),max(y1(i,j),0.0))
          else
!               tair(i,j)=(pt(i,j)+tb0)*pi0
!               rtair(i,j)=1./(tair(i,j)-c358)
             qsw(i,j)=rp0*exp(c172-c409*rtair(i,j))
#ifdef new_saturation
             esw(i,j)=min(0.99*p0_mks(i,k,j),esw_mks(tair(i,j)))
             qsw(i,j)=0.622*esw(i,j)/(p0_mks(i,k,j)-esw(i,j))
             esw(i,j) = esw(i,j)*10.  !in CGS
#endif
             ssw(i,j)=(qv(i,j)+qb0)/qsw(i,j)-1.0
             dm(i,j)=qv(i,j)+qb0-qsw(i,j)
             rsub1(i,j)=cv409*qsw(i,j)*rtair(i,j)*rtair(i,j)
             dd1(i,j)=max(-dm(i,j)/(1.+rsub1(i,j)), 0.0)
             y1(i,j)=.78/zr(i,j)**2+r23af*scv(i,j)/zr(i,j)**bwh5
             y2(i,j)=r23br/(tca(i,j)*tair(i,j)**2)+1./(dwv(i,j)   &
                       *qsw(i,j))
!cccc
             ern(i,j)=r23t*ssw(i,j)*y1(i,j)/y2(i,j)
             ern(i,j)=min(dd1(i,j),qr(i,j),max(ern(i,j),0.))

!JJS 20090708 vvvvv
! xli - reducing evaporation rate according to fitq
!             if (ihail .eq. 1 .and. qr(i,j) .gt. 1.e-6 ) then
              if (qr(i,j) .gt. 1.e-6 ) then  !Tao 20120224 apply to all options
                fact_fit=0.11*(qr(i,j)*1.e3)**(-1.27)+0.98
                ern(i,j)=ern(i,j)/fact_fit
             endif
!JJS 20090708 ^^^^^
          endif
              pt(i,j)=pt(i,j)-avcp*ern(i,j)
              qv(i,j)=qv(i,j)+ern(i,j)
              qr(i,j)=qr(i,j)-ern(i,j)
        endif
!>>> end of ifndef sat_predict
#endif

!!!!
!!  add processes 30 & 33 fpr pmltg and pmlts
!!
!!!!

!* 30 * pmltg : evaporation of melting qg                         **30**
!* 33 * pmlts : evaporation of melting qs                         **33**

          pmlts(i,j)=0.0
          pmltg(i,j)=0.0
	     tair(i,j)=(pt(i,j)+tb0)*pi0
          tairc(i,j)=tair(i,j)-t0
#ifdef use_cpm
          cpm = cp + cvap*qv(i,j) + cliq*(qc(i,j)+qr(i,j)) + &
                cice*(qi(i,j)+qs(i,j)+qg(i,j))
          hlv = alv - (cliq-cvap)*(tair(i,j)-t0)
          hlf = alf - (cice-cliq)*(tair(i,j)-t0)
          hls = hlv + hlf
          avcp = hlv/cpm*pir
          ascp = hls/cpm*pir
          afcp = hlf/cpm*pir
#endif

            ftns0(i,j)=1.
            ftng0(i,j)=1.

          if (tair(i,j) .ge. t0) then
             ftns(i,j)=1.
             ftng(i,j)=1.
             if (improve.gt.2) then
                call sgmap(1,qs(i,j),r00,tairc(i,j),ftns0(i,j))
                call sgmap(2,qg(i,j),r00,tairc(i,j),ftng0(i,j))
                ftns(i,j)=ftns0(i,j)
                ftng(i,j)=ftng0(i,j)
             endif
!             rtair(i,j)=1./(tair(i,j)-c358)
             rtair(i,j)=1./(t0-c358)
             y2(i,j)=exp( c172-c409*rtair(i,j) )
             esw(i,j)=c610*y2(i,j)
             qsw(i,j)=rp0*y2(i,j)
#ifdef new_saturation
             esw(i,j)=min(0.99*p0_mks(i,k,j),esw_mks(tair(i,j)))
             qsw(i,j)=0.622*esw(i,j)/(p0_mks(i,k,j)-esw(i,j))
             esw(i,j) = esw(i,j)*10.  !in CGS
#endif
             ssw(i,j)=1.-(qv(i,j)+qb0)/qsw(i,j)
             dm(i,j)=qsw(i,j)-qv(i,j)-qb0
             rsub1(i,j)=cv409*qsw(i,j)*rtair(i,j)*rtair(i,j)
             dd1(i,j)=max(dm(i,j)/(1.+rsub1(i,j)),0.0)
             y3(i,j)=1./tair(i,j)
             dd(i,j)=y3(i,j)*(rn30a*y3(i,j)-rn10b)+rn10c*tair(i,j)       &
                                             /esw(i,j)
             y1(i,j)=ftng(i,j)*r30t*ssw(i,j)*(r191r/zg(i,j)**2+r192rf    &
                   /zg(i,j)**bgh5)/dd(i,j)
             pmltg(i,j)=min(qg(i,j),max(y1(i,j),0.0))
             y1(i,j)=ftns(i,j)*r33t*ssw(i,j)*(r331r/zs(i,j)**2+r332rf    &
                                           /zs(i,j)**bsh5)/dd(i,j)
             pmlts(i,j)=min(qs(i,j),max(y1(i,j),0.0))
             y1(i,j)=min(pmltg(i,j)+pmlts(i,j),dd1(i,j))
             pmltg(i,j)=y1(i,j)-pmlts(i,j)
             pt(i,j)=pt(i,j)-ascp*y1(i,j)
             qv(i,j)=qv(i,j)+y1(i,j)
             qs(i,j)=qs(i,j)-pmlts(i,j)
             qg(i,j)=qg(i,j)-pmltg(i,j)
          endif
! end   Processes 30 and 33
  
!JJS 10/7/2008     vvvvv
    ENDIF    ! part of if (iwarm.eq.1) then
!JJS 10/7/2008     ^^^^^

!            IF (QV(I,J)+QB0 .LE. 0.) QV(I,J)=-QB0
            if (qc(i,j) .le. cmin) qc(i,j)=0.
            if (qr(i,j) .le. cmin) qr(i,j)=0.
            if (qi(i,j) .le. cmin) qi(i,j)=0.
            if (qs(i,j) .le. cmin) qs(i,j)=0.
            if (qg(i,j) .le. cmin) qg(i,j)=0.
            dpt(i,j,k)=pt(i,j)
            dqv(i,j,k)=qv(i,j)
            qcl(i,j,k)=qc(i,j)
            qrn(i,j,k)=qr(i,j)
            qci(i,j,k)=qi(i,j)
            qcs(i,j,k)=qs(i,j)
            qcg(i,j,k)=qg(i,j)

#ifdef EXT_DIAG
         scc=0.
         see=0.

!            DD(I,J)=MAX(-CND(I,J), 0.)
!            CND(I,J)=MAX(CND(I,J), 0.)
!            DD1(I,J)=MAX(-DEP(I,J), 0.)

!ccshie 2/21/02 shie follow tao
!CC for reference    QI(I,J)=QI(I,J)-Y2(I,J)+PIDEP(I,J)*D2T
!CC for reference    QV(I,J)=QV(I,J)-(PIDEP(I,J))*D2T

!c            DEP(I,J)=MAX(DEP(I,J), 0.)
!            DEP(I,J)=MAX(DEP(I,J), 0.)+PIDEP(I,J)*D2T
!            SCC=SCC+CND(I,J)
!            SEE=SEE+DD(I,J)+ERN(I,J)

!         SC(K)=SCC+SC(K)
!         SE(K)=SEE+SE(K)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     henry:  please take a look  (start)
!JJS modified by JJS on 5/1/2007  vvvvv

            dd(i,j)=max(-cnd(i,j), 0.)
            cnd(i,j)=max(cnd(i,j), 0.)
!            dd1(i,j)=max(-dep(i,j), 0.)+pidep(i,j)*d2t
	    dd1(i,j)=max(-dep(i,j), 0.)  !bug fix by Di
            dep(i,j)=max(dep(i,j), 0.)

!!!!!!!!!!!DDDDDDDDDDDDDD double check by Lang 02/23/2016
              sccc=cnd(i,j)
              seee=dd(i,j) + ern(i,j)
              sddd=dep(i,j) + dmax1(pint(i,j),0.0) + psdep(i,j) + pgdep(i,j)
              ssss=dd1(i,j) - dmin1(pint(i,j),0.0) + pssub(i,j) + pgsub(i,j) + pmlts(i,j) + pmltg(i,j)
              smmm=psmlt(i,j) + pgmlt(i,j) + pimlt(i,j) + qracs(i,j) 
              sfff=psacw(i,j)*d2t + piacr(i,j)*d2t + psfw(i,j)*d2t + pgfr(i,j)*d2t   &
                  +dgacw(i,j)*d2t + dgacr(i,j)*d2t + psacr(i,j)*d2t + pihom(i,j) &
                  +pidw(i,j)+pimm(i,j) + pcfr(i,j) + pihms(i,j)*d2t + pihmg(i,j)*d2t


! Snapshot values (K/s), consistent with declared units in Registry
           physc(i,k,j) = avcp * sccc / d2t
           physe(i,k,j) = avcp * seee / d2t
           physd(i,k,j) = ascp * sddd / d2t
           physs(i,k,j) = ascp * ssss / d2t
           physf(i,k,j) = afcp * sfff / d2t
           physm(i,k,j) = afcp * smmm / d2t
! Accumulated values (K)
           acphysc(i,k,j) = acphysc(i,k,j) + avcp * sccc 
           acphyse(i,k,j) = acphyse(i,k,j) + avcp * seee 
           acphysd(i,k,j) = acphysd(i,k,j) + ascp * sddd 
           acphyss(i,k,j) = acphyss(i,k,j) + ascp * ssss 
           acphysf(i,k,j) = acphysf(i,k,j) + afcp * sfff 
           acphysm(i,k,j) = acphysm(i,k,j) + afcp * smmm 
!JJS modified by JJS on 5/1/2007  ^^^^^

!JJS   2010/10/19  vvvvv
!      radar reflectivity calculation

!        dbz(i,k,j)=0.0

        a_1=1.e6*r00*qr(i,j)
        a_2=1.e6*r00*qs(i,j)
        a_3=1.e6*r00*qg(i,j)

!   if (a_1+a_2+a_3 .ge. 1.e-8)  then

      tair(i,j)=(pt(i,j)+tb0)*pi0
      tairc(i,j)=tair(i,j) - t0
!vvvvvvvvvvvv 20110816 vvvvvvvvvvvvvvvvv
      uwet=4.464**0.95
!^^^^^^^^^^^^ 20110816 ^^^^^^^^^^^^^^^^^
      if (improve .eq. 3) then      !using new dbz formula from Xiping
            ftns(i,j)=1.
            ftns0(i,j)=1.
            ftng(i,j)=1.
            ftng0(i,j)=1.
            call sgmap(1,qs(i,j),r00,tairc(i,j),ftns0(i,j))
            call sgmap(2,qg(i,j),r00,tairc(i,j),ftng0(i,j))
         A_C = 49.6 * 1.E-6
         A_I = 9.4 * 1.E-6
         RE_C = 10.            ! UM
         RE_I = 25.            ! UM
         RE_S = 75.            ! UM
!vvvvvvvvvvvvv Tao 20110722
!         tairc5=min(0.,max(-35.,tairc(i,j)+5.0))
!         fexp=-0.2-tairc5/100.
!^^^^^^^^^^^^^ Tao 20110722 
         ucor=3071.29/tnw**0.75

! Caculate temperature dependent ucog
!         ftng(i,j)=1.
!vvvvvvvvvvvvv Tao 20110722
         ftng(i,j)=ftng0(i,j)**0.75
!         ftng(i,j)=exp(-1.*tslopeg*tairc5)
!         ftngQ=1.0
!         if (qg(i,j).gt.cmin) ftngQ=(qg(i,j)/0.001)**fexp
!         ftng(i,j)=ftng(i,j)*ftngQ
!^^^^^^^^^^^^^ Tao 20110722

! Caculate temperature dependent ucos
!         ftns(i,j)=1.
!vvvvvvvvvvvvv Tao 20110722
         ftns(i,j)=ftns0(i,j)**0.75
!         ftns(i,j)=exp(-1.*tslopes*tairc5)
!         ftnsQ=1.0
!         if (qs(i,j).gt.cmin) ftnsQ=(qs(i,j)/0.001)**fexp
!         ftns(i,j)=ftns(i,j)*ftnsQ
!^^^^^^^^^^^^^ Tao 20110722
!
!vvvvvvvvvvvvvv 20110816 vvvvvvvvvvvvvvvvvvvvvvvv
!          ucos=687.97*roqs**0.25/ftns(i,j)**0.75
!          ucog=687.97*roqg**0.25/ftng(i,j)**0.75
          ucos=687.97*roqs**0.25/tns**0.75/ftns(i,j)
          ucog=687.97*roqg**0.25/tng**0.75/ftng(i,j)
!^^^^^^^^^^^^^^ 20110816 ^^^^^^^^^^^^^^^^^^^^^^^^^

          a_11=ucor*(max(1.e-12,a_1))**1.75
!         a_22=ucos*(max(1.e-12,a_2))**1.75
!         a_33=ucog*(max(1.e-12,a_3))**1.75

!vvvvvvvvvvvvvv 20110816 vvvvvvvvvvvvvvvvvvvvvvvv 
!         a_22=ucos*(max(1.e-12,a_2))**1.75/ftns(i,j)
!         a_33=ucog*(max(1.e-12,a_3))**1.75/ftng(i,j)
         a_22=ucos*(max(1.e-12,a_2))**1.75
         a_33=ucog*(max(1.e-12,a_3))**1.75
!vvvvvvvvvvvvvv 20110816 vvvvvvvvvvvvvvvvvvvvvvvv

         W_C = r00 * QC(I,J) * 1000.   !LIQUID WATER CONTENT IN G/M^3
         W_I = r00 * QI(I,J) * 1000.   !ICE WATER CONTENT IN G/M^3
         W_S = r00 * QS(I,J) * 1000.   !ICE WATER CONTENT IN G/M^3

! Xiping's formula
!         ZE_CLD = A_C * W_C * RE_C**3 + A_I * W_I * RE_I**3   &  !ZE IN  mm^6/m^3
!                 + A_I * W_S * RE_S**3  ! (as cloud particles NOT precipitation)
         ZE_CLD = A_C * W_C * RE_C**3 + A_I * W_I * RE_I**3      !ZE IN  mm^6/m^3
                  ! (use a_22 and ucos instead)
         IF (TAIR(I,J).LT.273.16) THEN
!            ZDRY = MAX(1.e-4,A_11+A_33+ZE_CLD) ! Xiping's  !rain,snow,cloud ice,cloud water,graupel
            ZDRY = MAX(1.e-4,A_11+A_22+A_33+ZE_CLD) !rain,snow,cloud ice,cloud water,graupel
!            DBZ(I,K,J) = 10.*ALOG10(ZDRY)
            DBZ(I,K,J) = 10.*DLOG10(ZDRY)
         ELSE         
!            A_44 = A_11+UWET*(A_22+A_33)**.95         ! old formula
!            A_44 = A_11+UWET*A_33**.95+ZE_CLD         ! Xiping's
            A_44 = A_11+UWET*(A_22+A_33)**.95+ZE_CLD
            ZWET = MAX(1.e-4,A_44)
!            DBZ(I,K,J) = 10.*ALOG10(ZWET)
            DBZ(I,K,J) = 10.*DLOG10(ZWET)
         ENDIF

      else ! using old dbz formula from Tao
         ucor=3071.29/tnw**0.75
         ucos=687.97*roqs**0.25/tns**0.75
         ucog=687.97*roqg**0.25/tng**0.75
         uwet=4.464**0.95
         a_11=ucor*(max(1.e-5,a_1))**1.75
         a_22=ucos*(max(1.e-5,a_2))**1.75
         a_33=ucog*(max(1.e-5,a_3))**1.75
         zdry=max(1.e-20,a_11+a_22+a_33)
         dbz(i,k,j)=10.*log10(zdry)
         if (tair(i,j).ge.273.16) then
            a_44=a_11+uwet*(a_22+a_33)**.95
            zwet=max(1.e-20,a_44)
            dbz(i,k,j)=10.*log10(zwet)
         endif

      endif  !improve
#endif

!   endif

!JJS   2010/10/19  ^^^^^
#if ( WRF_CHEM == 1)
      ! EMK...Nuclei only calculated when coupling
      if ( (chem_opt == 300 .or. chem_opt == 301 .or. &
            chem_opt == 302 .or. chem_opt == 303) .and. &
            (gsfcgce_gocart_coupling == 1) ) then
         icn_diag(i,k,j) = icn_out  ! #/Litre 
         nc_diag(i,k,j) = ccn_out  ! #/cm3
      else
         icn_diag(i,k,j) = 0.
         nc_diag(i,k,j) = 0.
      end if         
#endif

!JJS 20140305 vvvvv  Calculate effective radius for all cloud species
!   eff_rad is a function of the slope parameter (Lambda)
        
    ! rain
      if (qrn(i,j,k) .lt. crmin) then
         refr(i,k,j) = 0.e0  
      else 
         refr(i,k,j) = eff_rad(zr(i,j))    
      endif  
    ! snow
      if (qcs(i,j,k) .lt. csmin) then
         refs(i,k,j) = 0.e0 
      else 
         refs(i,k,j) = eff_rad(zs(i,j))    
      endif  
    ! graupel/hail
      if (qcg(i,j,k) .lt. cgmin) then
         refg(i,k,j) = 0.e0  
      else 
         refg(i,k,j) = eff_rad(zg(i,j))    
      endif  

! for cloud water

   if ( rewflag .eq. 1 ) then
   ! default
   if (qcl(i,j,k) .lt. cmin) then
      refc(i,k,j) = 0.e0
   else
      L_cloud = qcl(i,j,k) * rho(i,j,k)             ! cloud water [g/cm3]
#if ( WRF_CHEM == 1)
   ! when running with WRF_Chem and using aerosol coupling in Goddard MP
           ! cpi: const_pi = 4.*atan(1.)         ~ 3.1415
           ! roqr: 1.0 g/cm**3, liquid water density
           ! roqi: 0.9179, ice density
            if ( (chem_opt == 300 .or. chem_opt == 301 .or. &
                   chem_opt == 302 .or. chem_opt == 303) .and. &
                   (gsfcgce_gocart_coupling == 1) ) then
              ! for cloud water, estimate lambda (slope of gamma distribution)
                   mu = min(15.e0, (1000.E0/ccn_out + 2.e0))
                   gamfac3 = ( gamma_toshi(mu+4.e0) / gamma_toshi(mu+3.e0) )
                   gamfac1 = ( gamma_toshi(mu+4.e0) / gamma_toshi(mu+1.e0) )
                   lambda = (4.e0/3.e0*cpi*roqr*ccn_out/L_cloud*   &
                            gamfac1)**(1.e0/3.e0)  ! [1/cm]
                   refc(i,k,j) = 1.e0/lambda * gamfac3 * 1.e4  !effective radius [micron]
             else
   ! when running with WRF_Chem but no aerosol coupling in Goddard MP
               if (xland(i,j) .eq. 1.0) then
                  ccn_ref = ccn_over_land
               else if (xland(i,j) .eq. 2.0) then
                  ccn_ref = ccn_over_water
               else
                  print *,' xland is not 1. or 2., run stopped'
                  stop
               endif
                 ! for cloud water, estimate lambda (slope of gamma distribution)
                      mu = min(15.e0, (1000.E0/ccn_ref + 2.e0))
                      gamfac3 = ( gamma_toshi(mu+4.e0) / gamma_toshi(mu+3.e0) )
                      gamfac1 = ( gamma_toshi(mu+4.e0) / gamma_toshi(mu+1.e0) )
                      lambda = (4.e0/3.e0*cpi*roqr*ccn_ref/L_cloud*   &
                               gamfac1)**(1.e0/3.e0)  ! [1/cm]
                      refc(i,k,j) = 1.e0/lambda * gamfac3 * 1.e4  !effective radius [micron]
            endif ! chem_opt and gsfcgce_gocart_coupling
#else
  ! Not running with WRF_Chem
            ! ccn_over_land = 1500  ! [#/cm3] climatological value
            ! ccn_over_water = 150  ! [#/cm3] climatological value
            if (xland(i,j) .eq. 1.0) then
               ccn_ref = ccn_over_land
            else if (xland(i,j) .eq. 2.0) then
               ccn_ref = ccn_over_water
            else
               print *,' xland is not 1. or 2., run stopped'
               stop
            endif
           ! for cloud water, estimate lambda (slope of gamma distribution)
                   mu = min(15.e0, (1000.E0/ccn_ref + 2.e0))
                   gamfac3 = ( gamma_toshi(mu+4.e0) / gamma_toshi(mu+3.e0) )
                   gamfac1 = ( gamma_toshi(mu+4.e0) / gamma_toshi(mu+1.e0) )
                   lambda = (4.e0/3.e0*cpi*roqr*ccn_ref/L_cloud*   &
                            gamfac1)**(1.e0/3.e0)  ! [1/cm]
                   refc(i,k,j) = 1.e0/lambda * gamfac3 * 1.e4  !effective radius [micron]
#endif
   endif ! qcl(i,j,k) < cmin test
   endif

   if ( rewflag .eq. 2 ) then
      ! mapping spectrum, different over land and ocean
      if (qcl(i,j,k) .ge. cwmin) then
         if ( xland(i,j) .eq. 1. ) then
            mdc1 = 5.8936819
            mdc2 = -7.013013
            mdc3 = 1.3178721
            mdc4 = 1.1741987
            mdc5 = 2.6110916E-3
            mdc6 = -0.26646396
         else
            mdc1 = 173.57305
            mdc2 = -64.370929
            mdc3 = 0.36833626
            mdc4 = 6.1389254
            mdc5 = 5.5915321E-3
            mdc6 = -0.12488698
         endif
         efd1 = 1.6855156
         efd2 = 0.84147302
         efd3 = -0.46783369
         efd4 = 1.005305E-2
         efd5 = 3.4833332E-2
         efd6 = 1.4727223E-2

         ltk  = log(tair(i,j))
         lqc  = -1.*log(qcl(i,j,k)*rho_mks(i,k,j))
         ltk2 = ltk*ltk
         lqc2 = lqc*lqc
         mvdc = exp(mdc1 + mdc2*ltk + mdc3*lqc + mdc4*ltk2     &
                + mdc5*lqc2 + mdc6*ltk*lqc)
         efdc = exp(efd1 + efd2*log(mvdc) + efd3*log(1000.)    &
                + efd4*log(mvdc)**2. + efd5*log(1000.)**2.     &
                + efd6*log(mvdc)*log(1000.))
         refc(i,k,j) = min(max(efdc/2.,0.5),50)
      else
         refc(i,k,j) = 0.5  !test
      endif
   endif

! for cloud ice
   if ( reiflag .eq. 1 ) then
!   if (qci(i,j,k) .lt. cmin) then
   if (qci(i,j,k) .lt. cimin) then
      refi(i,k,j) = 0.e0
   else
#if ( WRF_CHEM == 1)
!!! No arosol coupling for calaculating eff radius for ice for the time being 20141027 
!  ! when running with WRF_Chem and using aerosol coupling in Goddard MP
!      I_cloud = qci(i,j,k) * rho(i,j,k)             ! cloud ice [g/cm3]
!      if ( (chem_opt == 300 .or. chem_opt == 301 .or. &
!          chem_opt == 302 .or. chem_opt == 303) .and. &
!          (gsfcgce_gocart_coupling == 1) ) then
!        ! for cloud ice, estimate lambda (slope of gamma distribution)
!         mu = min(15.e0, (1000.E0/(icn_out*1.e-3)  + 2.e0))
!         gamfac3 = ( gamma_toshi(mu+4.e0) / gamma_toshi(mu+3.e0) )
!         gamfac1 = ( gamma_toshi(mu+4.e0) / gamma_toshi(mu+1.e0) )
!         lambda = (4.e0/3.e0*cpi*roqi*icn_out*1.e-3/I_cloud*  &
!                  gamfac1)**(1.e0/3.e0)  ! [1/cm]
!         refi(i,k,j) = 1.e0/lambda * gamfac3 * 1.e4  !effective radius [micron]
!         if (mod(itimestep,i24h).eq.1) then
!            if (refi(i,k,j) .gt. 500.) then
!               print *,'i ,j ,k = ', i, j, k  
!               print *,'I_cloud, qci, rho, icn_out = ', I_cloud, icn_out, qci(i,j,k), rho(i,j,k)
!               print *,'mu, gamfac3, gamfac1 = ', mu, gamfac3, gamfac1
!               print *,'lambda, refi = ', lambda, refi 
!            endif
!         endif 
!      else
  ! when running with WRF_Chem but no aerosol coupling in Goddard MP
        ! for cloud ice effective radius depends on temperature profile, formula from GCE
         refi(i,k,j) = 125.e0 +(tair(i,j)-243.16)*5.e0     ! [micron]
         if (tair(i,j) .gt. 243.16) refi(i,k,j) = 125.e0
         if (tair(i,j) .lt. 223.16) refi(i,k,j) = 25.e0
!      endif ! chem_opt and gsfcgce_gocart_coupling
#else
  ! Not running with WRF_Chem
     ! for cloud ice effective radius depends on temperature profile, formula from GCE
      refi(i,k,j) = 125.e0 +(tair(i,j)-243.16)*5.e0     ! [micron]
      if (tair(i,j) .gt. 243.16) refi(i,k,j) = 125.e0
      if (tair(i,j) .lt. 223.16) refi(i,k,j) = 25.e0
#endif
   endif ! qci(i,j,k) < cmin test
   endif

   if ( reiflag .eq. 2 ) then
      ! Wyser 1998 , equation 15 and 35:
      reimin = 10.0
      if ( qci(i,j,k) .ge. cimin ) then
         iwc_0 = 50.e-6  ! note that rho here is in CGS
                         ! IWC_0 (50 g/m^-3) must be converted to g/cm^3
         bw98 = - 2. + 1.e-3 * log10(rho(i,j,k)*qci(i,j,k)/iwc_0)*max(0.0,-tairc(i,j))**1.5
         refi(i,k,j) = 377.4 + bw98 * (203.3 + bw98 * (37.91 + 2.3696 * bw98))
         refi(i,k,j) = max (reimin, refi(i,k,j))
      else
         refi(i,k,j) = 0.0
      endif
   endif

   if ( reiflag .eq. 3 ) then
      ! Heymsfield et al. 2014 , equation 9e :
      reimin = 10.0
      if ( qci(i,j,k) .ge. cimin ) then
         if ( tairc(i,j) >= -56.0 .and. tairc(i,j) < 0.0 ) then
            refi(i,k,j) = 308.4 * exp ( 0.0152 * tairc(i,j) )      ! 131.657 ~ 308.4 micron
         elseif ( tairc(i,j) >= -71.0 .and. tairc(i,j) < -56.0 ) then
            refi(i,k,j) = 9.1744e+4 * exp ( 0.117 * tairc(i,j) )   ! 22.641 ~ 130.942 micron
         elseif ( tairc(i,j) < -71.0 ) then
            refi(i,k,j) = 83.3 * exp ( 0.0184 * tairc(i,j) )       ! 10 ~ 22.557 micron
         endif
         refi(i,k,j) = max (reimin, refi(i,k,j))
      else
         refi(i,k,j) = 0.0
      endif
   endif

   if ( reiflag .eq. 4 ) then
      reimin = 0.0
      ! Mitchell et al. 2011, equation 9
      ! IWC in mg/m^3, T in oC, rei in micron
      if ( qci(i,j,k) .ge. cimin ) then
         refi(i,k,j) = 132.23 + 1.2433 * tairc(i,j) + &
                       8.6629 * (6. + log10(rho(i,j,k) * qci(i,j,k)))
         refi(i,k,j) = max (reimin, refi(i,k,j))
      else
         refi(i,k,j) = 0.0
      endif
   endif

   if ( reiflag .eq. 5 ) then
      ! mapping spectrum, different over land and ocean
      if ( qci(i,j,k) .ge. cimin ) then
         ltk  = log(tair(i,j))
         lqi  = -1.*log(qci(i,j,k)*rho_mks(i,k,j))
         ltk2 = ltk*ltk
         lqi2 = lqi*lqi
         if ( xland(i,j) .eq. 1.0 ) then
            ! over land
            efdi = exp(161.47584 - 0.26232591*lqi + 4.3393883E-4*lqi2 &
                   - 57.057846*ltk + 5.3668153*ltk2)/10.
         else
            ! over ocean
            efdi = exp(161.92213 - 0.26232591*lqi + 4.3393883E-4*lqi2 &
                   - 57.057846*ltk + 5.3668153*ltk2)/10.
         endif
         refi(i,k,j) = efdi/2.
         refi(i,k,j) = min( max (refi(i,k,j),1.5),500.)
      else
         refi(i,k,j) = 1.5  !test
      endif
   endif
!JJS 20140305 ^^^^^  Calculate effective radius for all cloud species

 2000 continue

 1000 continue

!JJS  ****************************************************************
!JJS  convert from GCE grid back to WRF grid
      do k=kts,kte
         do j=jts,jte
         do i=its,ite
         ptwrf(i,k,j) = dpt(i,j,k)
         qvwrf(i,k,j) = dqv(i,j,k)
         qlwrf(i,k,j) = qcl(i,j,k)
         qrwrf(i,k,j) = qrn(i,j,k)
         qiwrf(i,k,j) = qci(i,j,k)
         qswrf(i,k,j) = qcs(i,j,k)
         qgwrf(i,k,j) = qcg(i,j,k)
!         icn_diag(i,k,j) = icn_cgs(i,j,k) * 1000. ! #/Litre <-- #/cm3
!         nc_diag(i,k,j) = nc_cgs(i,j,k)  ! #/cm3
!         if (gid.eq.1)  icn_diag(i,k,j) = 0.
         enddo !i
         enddo !j
      enddo !k

!       do j=jts,jte
!       do i=its,ite
!          comdbz(i,j) =  dbz(i,kts,j)
!          do k=kts,kte
!             if (dbz(i,k,j) .gt. comdbz(i,j)) comdbz(i,j) = dbz(i,k,j)
!          enddo !k
!       enddo !i
!       enddo !j

!     ****************************************************************

!+---+-----------------------------------------------------------------+
#ifdef EXT_DIAG
! EMK NUWRF...Replace Greg Thompson's dBZ calculations with values calculated
! above.
         IF ( PRESENT (diagflag) ) THEN
         if (diagflag .and. do_radar_ref == 1) then
            do j=jts,jte
               do k = kts,kte
                  do i=its,ite                     
                     refl_10cm(i,k,j) = max(-35.,dBZ(i,k,j))
                  end do
               end do
            end do
!             do j=jts,jte
!             do i=its,ite
!                DO K=kts,kte
!                   t1d(k)=ptwrf(i,k,j)*pi_mks(i,k,j)
!                   p1d(k)=p0_mks(i,k,j)
!                   qv1d(k)=qvwrf(i,k,j)
!                   qr1d(k)=qrwrf(i,k,j)
!                ENDDO
!                if (ice2.eq.0) then
!                   DO K=kts,kte
!                      qs1d(k)=qswrf(i,k,j)
!                      qg1d(k)=qgwrf(i,k,j)
!                   ENDDO
!                elseif (ice2.eq.1) then
!                   DO K=kts,kte
!                      qs1d(k)=qswrf(i,k,j)
!                   ENDDO
!                elseif (ice2.eq.2) then
!                   DO K=kts,kte
!                      qs1d(k)=0.
!                      qg1d(k)=qgwrf(i,k,j)
!                   ENDDO
!                elseif (ice2.eq.3) then
!                   DO K=kts,kte
!                      qs1d(k)=0.
!                      qg1d(k)=0.
!                   ENDDO
!                endif
!                call refl10cm_gsfc (qv1d, qr1d, qs1d, qg1d,             &
!                        t1d, p1d, dBZ_GT, kts, kte, i, j, ihail)
!                do k = kts, kte
!                   refl_10cm(i,k,j) = MAX(-35., dBZ_GT(k))
!                enddo
!             enddo
!             enddo
         endif
         ENDIF
#endif
!+---+-----------------------------------------------------------------+

      return

  END SUBROUTINE saticel_s

  SUBROUTINE auto_conversion( L, N, P , re)
  implicit none
!-----------------------------------------------------------------------------------------------------
! Comments: 
!  This subroutine compute auto conversion rate folloing Li and Daum [2004], which account for 
!  total cloud liquid water, particle number concentrations, PSD lambda, broadening parameters. 
! 
! History:
!  08/2010  Toshi Matsui@NASA GSFC : Initial.
!
!
! References:
! Liu, Y. and P. H. Daum, 2004: Parameterization of the autoconversion process. Part I: Analytical 
!   formulation of the Kessler-type parameterizations. J. Atmos. Sci, 61, 1539-1548.
!-----------------------------------------------------------------------------------------------------
 real,intent(in) :: L    ! cloud liquid water [g cm-3]
 real,intent(in) :: N    ! total number concentration [# cm-3]
 real,intent(out) :: P   ! auto conversion rate [g cm-3 s-1]
 real,intent(out) :: re  ! cloud effective radius [micron]

 real :: mu   ! mu of gamma PSD [-]
 real :: eta  ! eta function [cm3 g-2 s-1]
 real :: beta, beta1, beta2     ! beta function [-]
 real :: gamfac , gfac1 , gfac2 ! gamma function [-]
 real :: R6_6power  ! mean radius of the sixth moment [cm]
 real :: R6_thresh ! threshold of  mean radius of the sixth moment [cm]
 real :: R6        ! mean radius of the sixth moment [cm]
 real :: Heaviside_func  ! Heaviside step function (0 or 1) 
 real :: lambda    ! slope of gamma size ditribution [1/cm]
! real :: No        ! intercept  [cm-4]

 real,parameter :: Rc = 10.e0 * 1.e-4 ! threshold of particle radus (10 micron) [cm]
 real,parameter :: const_pi    = 3.14159e0 ! pai 
 real,parameter :: const_kappa = 1.9e11    ! coefficient for water droplet collection kernel [cm-3 s-1] 
                                           ! from Long [1974, JAS].
! real,parameter :: const_kappa = 1.9e11*10000.e0 !10000 is to adjust the order to keseller


 real,parameter :: const_rho_liq = 1.e0    ! density of liquid water [g cm-3] 
 real,parameter :: eta_func = ((3.e0/(4.e0*const_pi*const_rho_liq))**2) * const_kappa  ! eta function [cm3 g-2 s-1]
                                                                                       ! a part of (eq 27b)

 logical,parameter :: no_thresh = .true.  ! logic to choose no threshold parameterization or not.


!
! When no particel, no autoconversion. 
!
! EMK BUG FIX...Prevent overflow for small but non-zero values of L
! if( N <= 0.e0 .or. L <= 0.e0 ) then
 if( N <= 0.e0 .or. L <= 1.0e-32 ) then
   P = 0.e0
   return
 endif

!
! check bad values of N and L
!
! if( N < 0.e0 ) stop 'MSG auto_conversion: N is negative, it must be positive' 
! if( L < 0.e0 ) stop 'MSG auto_conversion: L is negative, it must be positive' 
! if( L > 1.e0 ) stop 'MSG auto_conversion: L is greater than 1g/cm3.'

!
! Empirical fit of mu as a function of total particle number concentrations. 
! From Martin et al. (1994), assign gamma shape parameter mu for cloud
! drops according to general dispersion characteristics. 
! disp=~0.25 for Maritime and 0.45 for Continental.
! Since disp=SQRT((mu+2)/(mu+1) - 1), mu varies from 15 for Maritime (pristine air)
! to 2 for Continental (really dirty air).  if mu = 0 --> expnential distribution (narrow dist)
!

! orig
 mu = MIN(15.e0, (1000.E0/N + 2.e0))


!
! gamma functions
!
 gfac1 = gamma_toshi(mu+4.e0)
 gfac2 = gamma_toshi(mu+1.e0)
 gamfac = (gfac1/gfac2)

!
! estimate lambda (slope of gamma distribution)
!
 lambda = (4.e0/3.e0*const_pi*const_rho_liq*N/L*gamfac)**(1.e0/3.e0)  ! [1/cm]


 THRESH: if( no_thresh ) then !-------------------------------------------

!
! threshold of particle radius (mean radius of the sixth moment )
!
 gfac1 = gamma_toshi(mu+7.e0)
 gfac2 = gamma_toshi(mu+1.e0)
 gamfac = (gfac1/gfac2)

 R6_6power = (1.e0 / lambda)**6.e0 * gamfac   ![cm] (eq. A3)

!
! auto conversion rate (eq. 26a)
!
 P = const_kappa *   N     * R6_6power *    L      ! [g cm-3 s-1 ]
!    [cm-3 s-1]  * [#/cm3] *   [cm6]   * [g/cm3]


 else  !with threshold ---------------------------------------------------

!
! Estimate eta under gamma PSD
!
 beta1 = (6.e0+mu)*(5.e0+mu)*(4.e0+mu) 
 beta2 = (3.e0+mu)*(2.e0+mu)*(1.e0+mu)
 beta  = beta1 / beta2

 eta = eta_func * beta  ! eta function (eq 27b) [cm3 g-2 s-1]

!
! threshold of particle radius (mean radius of the sixth moment )
!
 R6_thresh = beta * Rc  ![cm] (pg 1545)

!
! mean radius of the sixth moment 
!
 gfac1 = gamma_toshi(6.e0+mu+1.e0)
 gfac2 = gamma_toshi(1.e0+mu)
 gamfac = (gfac1/gfac2)**(1.e0/6.e0)

 R6 = (1.e0 / lambda) * gamfac   ![cm] (eq. A3)

!
! Heaviside step function
!
 if    ( R6 - R6_thresh <= 0.e0 ) then
    Heaviside_func = 0.e0
 elseif( R6 - R6_thresh > 0.e0 ) then
    Heaviside_func = 1.e0
 else
   print*, 'MSG: auto_conversion: Strange value of R6= ', R6 ; stop
 endif 

!
! auto conversion rate [g cm-3 s-1 ] (eq. 27a)
!
 P = eta * (1.e0/N) * (L**3) *  Heaviside_func

!    [cm3 g-2 s-1] * [cm3] * [g3/cm9]

 endif THRESH !------------------------------------------------------------


! optional 

!
! estimate effective radius
!
 gfac1 = gamma_toshi(mu+4.e0)
 gfac2 = gamma_toshi(mu+3.e0)
 gamfac = (gfac1/gfac2)

 re = 1.e0/lambda * gamfac * 1.e4  !effective radius [micron] 

!
! estimate No
!
! call gamma_function(mu+1.e0 ,gfac1)
! No = N * (lambda**(mu+1)) / gfac1  

 RETURN
 END subroutine auto_conversion

 real function gamma_toshi(x)

!---------------------------------------------------------------------------------------------------
! Comments:                             
!   compute the gamma function T(x) for single precision floating point. 
!       input :  x  --- argument of a(x)
!                       ( x is not equal to 0,-1,-2,... )
!                                       
! History:
! 09/2009  Toshi Matsui@NASA GSFC ; Adapted to SDSU
!           
! References: 
!----------------------------------------------------------------------------------------------------
 implicit double precision (a-h,o-z)
 dimension g(26)
 data g/1.0d0,0.5772156649015329d0, &
       -0.6558780715202538d0, -0.420026350340952d-1, &
        0.1665386113822915d0,-.421977345555443d-1, &
        -.96219715278770d-2, .72189432466630d-2, &
        -.11651675918591d-2, -.2152416741149d-3, &
        .1280502823882d-3, -.201348547807d-4, &
        -.12504934821d-5, .11330272320d-5, &
        -.2056338417d-6, .61160950d-8, &
         .50020075d-8, -.11812746d-8, &
        .1043427d-9, .77823d-11, &
        -.36968d-11, .51d-12, &
        -.206d-13, -.54d-14, .14d-14, .1d-15/
 real :: x

 pi=3.141592653589793d0
 if (x.eq.int(x)) then
     if (x.gt.0.0d0) then
         ga=1.0d0
         m1=int(x)-1
        do k=2,m1
           ga=ga*k
        enddo
     else
        ga=1.0d+300
     endif
  else
     if (dabs(dble(x)).gt.1.0d0) then
         z=dabs(dble(x))
         m=int(z)
         r=1.0d0
        do k=1,m
           r=r*(z-k)
        enddo
        z=z-m
     else
        z=dble(x)
     endif
     gr=g(26)
     do k=25,1,-1
        gr=gr*z+g(k)
     enddo
     ga=1.0d0/(gr*z)
     if (dabs(dble(x)).gt.1.0d0) then
         ga=ga*r
         if (x.lt.0.0d0) ga=-pi/(x*ga*dsin(pi*x))
     endif
  endif

  gamma_toshi = real(ga)

  return
  end function gamma_toshi

!+---+-----------------------------------------------------------------+

      subroutine refl10cm_gsfc (qv1d, qr1d, qs1d, qg1d,                 &
                       t1d, p1d, dBZ, kts, kte, ii, jj, ihail)

      IMPLICIT NONE

!..Sub arguments
      INTEGER, INTENT(IN):: kts, kte, ii, jj, ihail
      REAL, DIMENSION(kts:kte), INTENT(IN)::                            &
                      qv1d, qr1d, qs1d, qg1d, t1d, p1d
      REAL, DIMENSION(kts:kte), INTENT(INOUT):: dBZ

!..Local variables
      REAL, DIMENSION(kts:kte):: temp, pres, qv, rho
      REAL, DIMENSION(kts:kte):: rr, rs, rg

      DOUBLE PRECISION, DIMENSION(kts:kte):: ilamr, ilams, ilamg
      DOUBLE PRECISION, DIMENSION(kts:kte):: N0_r, N0_s, N0_g
      DOUBLE PRECISION:: lamr, lams, lamg
      LOGICAL, DIMENSION(kts:kte):: L_qr, L_qs, L_qg

      REAL, DIMENSION(kts:kte):: ze_rain, ze_snow, ze_graupel
      DOUBLE PRECISION:: fmelt_s, fmelt_g

      INTEGER:: i, k, k_0, kbot, n
      LOGICAL:: melti

      DOUBLE PRECISION:: cback, x, eta, f_d
      REAL, PARAMETER:: R=287.
      REAL, PARAMETER:: PIx=3.1415926536

!+---+

      do k = kts, kte
         dBZ(k) = -35.0
      enddo

!+---+-----------------------------------------------------------------+
!..Put column of data into local arrays.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         temp(k) = t1d(k)
         qv(k) = MAX(1.E-10, qv1d(k))
         pres(k) = p1d(k)
         rho(k) = 0.622*pres(k)/(R*temp(k)*(qv(k)+0.622))

         if (qr1d(k) .gt. 1.E-9) then
            rr(k) = qr1d(k)*rho(k)
            N0_r(k) = xnor
            lamr = (xam_r*xcrg(3)*N0_r(k)/rr(k))**(1./xcre(1))
            ilamr(k) = 1./lamr
            L_qr(k) = .true.
         else
            rr(k) = 1.E-12
            L_qr(k) = .false.
         endif

         if (qs1d(k) .gt. 1.E-9) then
            rs(k) = qs1d(k)*rho(k)
            N0_s(k) = xnos
            lams = (xam_s*xcsg(3)*N0_s(k)/rs(k))**(1./xcse(1))
            ilams(k) = 1./lams
            L_qs(k) = .true.
         else
            rs(k) = 1.E-12
            L_qs(k) = .false.
         endif

         if (qg1d(k) .gt. 1.E-9) then
            rg(k) = qg1d(k)*rho(k)
            if (ihail.eq.1) then
               N0_g(k) = xnoh
            else
               N0_g(k) = xnog
            endif
            lamg = (xam_g*xcgg(3)*N0_g(k)/rg(k))**(1./xcge(1))
            ilamg(k) = 1./lamg
            L_qg(k) = .true.
         else
            rg(k) = 1.E-12
            L_qg(k) = .false.
         endif
      enddo

!+---+-----------------------------------------------------------------+
!..Locate K-level of start of melting (k_0 is level above).
!+---+-----------------------------------------------------------------+
      melti = .false.
      k_0 = kts
      do k = kte-1, kts, -1
         if ( (temp(k).gt.273.15) .and. L_qr(k)                         &
                                  .and. (L_qs(k+1).or.L_qg(k+1)) ) then
            k_0 = MAX(k+1, k_0)
            melti=.true.
            goto 195
         endif
      enddo
 195  continue

!+---+-----------------------------------------------------------------+
!..Assume Rayleigh approximation at 10 cm wavelength. Rain (all temps)
!.. and non-water-coated snow and graupel when below freezing are
!.. simple. Integrations of m(D)*m(D)*N(D)*dD.
!+---+-----------------------------------------------------------------+

      do k = kts, kte
         ze_rain(k) = 1.e-22
         ze_snow(k) = 1.e-22
         ze_graupel(k) = 1.e-22
         if (L_qr(k)) ze_rain(k) = N0_r(k)*xcrg(4)*ilamr(k)**xcre(4)
         if (L_qs(k)) ze_snow(k) = (0.176/0.93) * (6.0/PIx)*(6.0/PIx)     &
                                 * (xam_s/900.0)*(xam_s/900.0)          &
                                 * N0_s(k)*xcsg(4)*ilams(k)**xcse(4)
         if (L_qg(k)) ze_graupel(k) = (0.176/0.93) * (6.0/PIx)*(6.0/PIx)  &
                                    * (xam_g/900.0)*(xam_g/900.0)       &
                                    * N0_g(k)*xcgg(4)*ilamg(k)**xcge(4)
      enddo


!+---+-----------------------------------------------------------------+
!..Special case of melting ice (snow/graupel) particles.  Assume the
!.. ice is surrounded by the liquid water.  Fraction of meltwater is
!.. extremely simple based on amount found above the melting level.
!.. Uses code from Uli Blahak (rayleigh_soak_wetgraupel and supporting
!.. routines).
!+---+-----------------------------------------------------------------+

      if (melti .and. k_0.ge.kts+1) then
       do k = k_0-1, kts, -1

!..Reflectivity contributed by melting snow
          if (L_qs(k) .and. L_qs(k_0) ) then
           fmelt_s = MAX(0.005d0, MIN(1.0d0-rs(k)/rs(k_0), 0.99d0))
           eta = 0.d0
           lams = 1./ilams(k)
           do n = 1, nrbins
              x = xam_s * xxDs(n)**xbm_s
              call rayleigh_soak_wetgraupel (x,DBLE(xocms),DBLE(xobms), &
                    fmelt_s, melt_outside_s, m_w_0, m_i_0, lamda_radar, &
                    CBACK, mixingrulestring_s, matrixstring_s,          &
                    inclusionstring_s, hoststring_s,                    &
                    hostmatrixstring_s, hostinclusionstring_s)
              f_d = N0_s(k)*xxDs(n)**xmu_s * DEXP(-lams*xxDs(n))
              eta = eta + f_d * CBACK * simpson(n) * xdts(n)
           enddo
           ze_snow(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
          endif

!..Reflectivity contributed by melting graupel

          if (L_qg(k) .and. L_qg(k_0) ) then
           fmelt_g = MAX(0.005d0, MIN(1.0d0-rg(k)/rg(k_0), 0.99d0))
           eta = 0.d0
           lamg = 1./ilamg(k)
           do n = 1, nrbins
              x = xam_g * xxDg(n)**xbm_g
              call rayleigh_soak_wetgraupel (x,DBLE(xocmg),DBLE(xobmg), &
                    fmelt_g, melt_outside_g, m_w_0, m_i_0, lamda_radar, &
                    CBACK, mixingrulestring_g, matrixstring_g,          &
                    inclusionstring_g, hoststring_g,                    &
                    hostmatrixstring_g, hostinclusionstring_g)
              f_d = N0_g(k)*xxDg(n)**xmu_g * DEXP(-lamg*xxDg(n))
              eta = eta + f_d * CBACK * simpson(n) * xdtg(n)
           enddo
           ze_graupel(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
          endif

       enddo
      endif

      do k = kte, kts, -1
         dBZ(k) = 10.*log10((ze_rain(k)+ze_snow(k)+ze_graupel(k))*1.d18)
      enddo

      end subroutine refl10cm_gsfc

!+---+-----------------------------------------------------------------+

!JJS 20140225
! Calculate cloud droplet effective radius
   real function eff_rad(lambda)

      implicit none

!---------------------------------------------------------------------------------------------------
! Comments:
! Compute drop effective radius from slope parameters (lambda) of expoential size distribution.
!
! History:
! 02/2014  Toshi Matsui@NASA GSFC ; Initial
!
! References:
!----------------------------------------------------------------------------------------------------
      real,intent(in)  :: lambda   ! intercept parameter [1/cm]
!      real,intent(out) :: re  !effective radius [micron]

!
! for no particles.
!
!       if ( lambda <= 0.e0  .or. isnan(lambda) ) then
       if ( lambda <= 0.e0 ) then
          eff_rad = 0.e0
          return
       endif

!
! compute drop effective radius for exponential distribution N(D) = N0*exp(-lam*D)
!
       eff_rad = 1.5e0 / (lambda*100.) * 1.0e+6  ! [micron]
   return
   end function eff_rad

!-------------------------------------------------------------------
      subroutine vti_mks(improve,rhoz,tz,qiz,qvz,p,xland,vti)
      implicit none
      integer, intent(in) :: improve
      real, intent(in) :: rhoz  !air density (kg/m^3)
      real, intent(in) :: qiz   !mixing ratio of cloud ice (kg/kg)
      real, intent(in) :: qvz   !specific humidity (kg/kg)
      real, intent(in) :: p     !pressure (Pa)
      real, intent(in) :: tz    !air temperature (K)
      real, intent(in) :: xland !land-sea mask
      real, intent(out):: vti   !terminal velocity (m/s)

      integer, parameter :: vtiflag = 7
      ! 1 : Starr and Cox (1985)        , igce = 1
      ! 2 : Heymsfield and Donner (1990), igce!= 1
      ! 3 : Hong et al. (2004)          , igce = 1 , improve = 3
      ! 4 : Deng and Mace (2008)
      ! 5 : Mitchell et al. (2011)
      ! 6 : assume rhoi=300 kg/m^3
      ! 7 : TCWA1 semi-theoretical approach, different over land and ocean

      real, parameter :: vimax = 0.5     ! max fall speed for cloud ice (m/s)
      real, parameter :: vimin = 0.      ! min fall speed for cloud ice (m/s)

      !local variables
      integer :: ic
      integer :: igce
      real    :: y1
      real    :: const_vt, const_d, const_m, bb1, bb2
      real, dimension(7) ::  aice, vice
      data aice/1.e-6, 1.e-5, 1.e-4, 1.e-3, 0.01, 0.1, 1./
      data vice/5,15,30,35,40,45,50/
      ! for Deng and Mace (2008)
      real, parameter :: aa = - 4.14122e-5
      real, parameter :: bb = - 0.00538922
      real, parameter :: cc = - 0.0516344
      real, parameter :: dd = 0.00216078
      real, parameter :: ee = 1.9714
      real    :: tc
      real    :: h1, h2
      ! for TCWA1 semi-theoretical approach
      integer :: hid
      real    :: qsi,sqrhoz,rhoi,adagr,inhgr,ltk,lqi,ltk2,lqi2,zeta,vishp,viroi,ssi,lroi
#ifdef new_saturation
      real    :: esi
#endif
      real, parameter :: di0 = 6.e-6

      if ( qiz .ge. cimin ) then
         if ( vtiflag .eq. 1 ) then
            ! Starr and Cox 1985 :
            y1 = rhoz * 1000. * qiz           ! y1 in g/m^3
            if ( y1 .lt. 1.e-6 ) then
               vti = 0.
            else
               do ic = 1 , 6
                  if ( y1.gt.aice(ic) .and. y1.le.aice(ic+1) ) then
                     vti = vice(ic) + ( vice(ic+1) - vice(ic) ) *         &
                             ( y1 - aice(ic) ) / ( aice(ic+1) - aice(ic) )
                     if (vti .lt. 0.0) vti = 0. ! EMK per 20110816 code
                  endif
                  vti = vti * 0.01  ! convert back to MKS
               enddo
            endif

         elseif ( vtiflag .eq. 2 ) then
            ! Heymsfield and Donner 1990 :
            vti= 3.29 * ( rhoz * qiz ) ** 0.16

         elseif ( vtiflag .eq. 3 ) then
            ! Hong et al. 2004, same as HD90 (in a*D**b form)
            y1 = rhoz * 1000. * qiz                  ! y1 in g/m^3
            if ( y1 .lt. 1.e-6 ) then
               vti = 0.
            else
               ! fallspeed     (m/s)  : V=1.49e4*D**1.31
               ! diameter      (m)    : D=11.9*M**0.5
               ! concentration (m^-3) : N=5.38e7*(rho*q)**0.75
               ! mass          (kg)   : M=rho*q/N
               !                         =1/5.38e7*(rho*q)**0.25
               const_vt = 1.49e4
               const_d = 11.9
               const_m = 1./5.38e7
               y1 = y1 * 1.e-3                       ! y1 in kg/m^3
               bb1 = const_m * y1**0.25
               bb2 = const_d * bb1**0.5
               vti = max(const_vt*bb2**1.31, 0.0)    ! vti in m/s
            endif

         elseif ( vtiflag .eq. 4 ) then
            ! Deng and Mace 2008 :
            ! IWC in g/m^3 , tc in oC , and vti in cm/s
            tc = tz - t0
            vti = (3. + log10(qiz * rhoz)) * &
                  (tc * (aa * tc + bb) + cc) + dd * tc + ee
            vti = exp(log(10.) * vti)
            vti = vti * 0.01    ! convert back to MKS

         elseif ( vtiflag .eq. 5 ) then
            ! Mitchell et al. 2011 , equation 10 :
            ! IWC in mg/m^3 , T in oC , and vti in cm/s
            vti = 82.082 + 1.0121 * tc + &
                  6.6303 * (6. + log10(qiz * rhoz))
            vti = vti * 0.01    ! convert back to MKS

         elseif ( vtiflag .eq. 6 ) then
            ! assume constant density (rhoi=300 kg/m^3) :
            tc = tz - t0
            sqrhoz = sqrt(rhoe_s/rhoz)

            ! shape parameter :
            hid   = max(min(nint(abs(tc)/0.25),120),0)
            adagr = itble(hid)
            zeta  = (adagr-1.)/(adagr+2.)
            if ( zeta .gt. 0. ) then
               vishp = di0**(zeta/2.)
            elseif ( zeta .lt. 0. ) then
               vishp = di0**(-zeta)
            else
               vishp = 1.
            endif

            ltk   = log(tz)
            ltk2  = ltk*ltk
            lqi   = -1.*log(rhoz*qiz)
            lqi2  = lqi*lqi
            vti = exp(-1.1100279E-2 + 0.47727519*lqi              &
                  - 8.8757389E-2*lqi2 + 3.6732918E-3*lqi*lqi2     &
                  - 4.5748034E-5*lqi2*lqi2 + 1.3864255*ltk)       &
                  /1.e+3*vishp*sqrhoz

         elseif ( vtiflag .eq. 7 ) then
            ! TCWA1 adopted semi-theoretical approach from NTU 4ICE-3M scheme
            ! with prescribed ice properties(shape and density) :
            tc = tz - t0
            qsi = f_qsi(tz,p)
#ifdef new_saturation
            esi = min(0.99*p,esi_mks(tz))
            qsi = 0.622*esi/(p-esi)
#endif

            ! deposition density :
            if ( tc .ge. -40. ) then
               ! Chen and Lamb 1994a ; Chen and Tsai 2016
               hid = max(min(nint(abs(tc)/0.25),120),0)
               inhgr = itble(hid)  !inherent growth ratio
               rhoi = 900.*exp(-3.*max(qvz-qsi-5.e-5,0.)/inhgr)
            else
               ! Pokrifka et al. 2023, equation 18
               inhgr = 3.          !inherent growth ratio (FIG. 16)
               ssi = qvz/qsi-1.
               if ( ssi .gt. 0.267 ) then
                  rhoi = -1027.456*ssi+1185.834
               else
                  rhoi = -32.332*ssi+900.
               endif
            endif

            ! shape parameter (ice aspect ratio) :
!            adagr = inhgr**thrd
! >>> reduce upper-level vti :
            if ( tc .ge. -40. ) then
              adagr = inhgr**thrd
            else
#ifdef sat_predict
              adagr = inhgr**thrd
#else
              adagr = inhgr**0.9
#endif
            endif
! <<<
            ltk   = log(tz)
            ltk2  = ltk*ltk
            lqi   = -1.*log(rhoz*qiz)
            lqi2  = lqi*lqi
            lroi  = log(min(max(rhoi,50.),900.))
            zeta  = (adagr-1.)/(adagr+2.)
            if ( zeta .gt. 0. ) then
               vishp = di0**(zeta/2.)
            elseif ( zeta .lt. 0. ) then
               vishp = di0**(-zeta)
            else
               vishp = 1.
            endif

            ! density parameter :
            viroi = 2.6795546 - 0.010732829*lqi - 1.176491*lroi         &
                    + 3.2512268E-4*lqi2 + 0.13920371*lroi**2.           &
                    - 5.6243169E-4*lqi*lroi
            viroi = min(1.,viroi)

            ! the divide between land and ocean :
            if ( xland .eq. 1. ) then   ! land
               vti = exp(265.16805 - 0.28802545*lqi - 3.754874E-3*lqi2  &
                     - 93.843132*ltk + 8.8315122*ltk2)                  &
                     /1.e+6*vishp*viroi*(1.0837/rhoz)**0.35
            else                        ! ocean
               vti = exp(252.86312 - 0.23844613*lqi - 4.6185936E-3*lqi2 &
                     - 89.119066*ltk + 8.3851678*ltk2)                  &
                     /1.e+6*vishp*viroi*(1.0837/rhoz)**0.35
            endif
         endif

         vti = min ( vimax , max ( vimin , vti ) )
      else
         vti = vimin
      endif  !end of if qiz

      return
      end subroutine vti_mks
!-------------------------------------------------------------------

!-------------------------------------------------------------------
      subroutine vtr_mks(rhoz,qrz,tz,vtr)
      implicit none
      real, intent(in) :: rhoz  !air density (kg/m^3)
      real, intent(in) :: qrz   !mixing ratio (kg/kg)
      real, intent(in) :: tz    !temperature (K)
      real, intent(out):: vtr   !terminal velocity (m/s)

      real, parameter :: vrmax = 12.0    !(m/s)
      real, parameter :: vrmin = 0.0     !(m/s)

      integer, parameter :: vtrflag = 2
      ! 1 : Lin et al. (1983)   ( igce != 1 )
      ! 2 : Steve's             ( igce = 1 )
      ! 3 : Lang et al. (2014) , Steve's with bin rain evaporation correction

      !local variables
!      integer :: igce
      real    :: pi, gambp4
      real    :: y1, tmp1, vs, vg, vr, fv, zr
      real    :: ftnw, bin_factor, ftnwmin

      pi = acos(-1.)
      gambp4 = gammagce(constb+4.)

      if ( qrz .ge. crmin ) then
         fv = sqrt( rhoe_s/rhoz )
!         igce = 1
         if ( vtrflag .eq. 1 ) then !if igce .ne. 1
            ! Lin et al. (1983) ; in MKS
            ! rhowater=1000., consta=0.8, constb=841.9967
!            tmp1 = sqrt(pi*rhowater*xnor/rhoz/qrz)
!            tmp1 = sqrt(tmp1)
!            vtr = consta*gambp4*fv/tmp1**constb
!            vtr = vtr/6.
            zr = (pi*rhowater*xnor/(rhoz*qrz))**0.25   !slope parameter of rain (in MKS)
            gambp4 = gammagce(constb+4.)
            vtr = consta*gambp4/(6.*zr**constb)*fv

         elseif ( vtrflag .eq. 2 ) then !if igce=1
            ! new codes from Steve's in CGS :
            y1 = rhoz*qrz*0.001  !in CGS
            vs = sqrt( y1 )
            vg = sqrt( vs )
            vr = vrc0+vrc1*vg+vrc2*vs+vrc3*vg*vs
            vtr = max(fv*vr, 0.e0)
            vtr = vtr*0.01  !convert back to MKS

         elseif ( vtrflag .eq. 3 ) then
            ! same as Steve's, but with bin rain evaporation correction (Lang et al. 2014)
            ! (in CGS)
            ftnw = 1.
            if ( tz .gt. t0 ) then
               bin_factor = 0.11*(1000.*qrz)**(-1.27) + 0.98
               bin_factor = min(bin_factor,1.30)
               ftnw = 1./bin_factor**3.35
               ftnwmin = rhoz*0.001*qrz/draimax
               if ( qrz .le. 0.001 ) ftnw = max(ftnw,ftnwmin/tnw)
            endif
            zr = (pi*roqr*tnw/(rhoz*0.001*qrz))**0.25   !slope parameter of rain (in CGS)
            zr = zr/ftnw
            vtr = (-26.7 + 2.06e+4/zr - 2.045e+5/zr**2. + 9.06e+5/zr**3.)*fv
            vtr = vtr*0.01  !convert back to MKS

         endif  !end of if igce
         vtr = min ( vrmax , max ( vrmin , vtr ) )
      else
         vtr = vrmin
      endif  !end of if qrz

      return
      end subroutine vtr_mks
!-------------------------------------------------------------------

!-------------------------------------------------------------------
      subroutine vts_mks(improve,rhoz,qsz,tz,vts)
      implicit none
      integer, intent(in) :: improve
      real, intent(in) :: rhoz  !air density (kg/m^3)
      real, intent(in) :: qsz   !mixing ratio (kg/kg)
      real, intent(in) :: tz    !temperature (K)
      real, intent(out):: vts   !terminal velocity (m/s)

      real, parameter :: vsmax = 5.0     !(m/s)
      real, parameter :: vsmin = 0.0     !(m/s)

      !local variables
      integer :: igce

      real    :: y1, fv, r00, ftns, ftns0, tzc, vscf
      real    :: tmp1, pi, gamdp4

      pi = acos(-1.)
      gamdp4 = gammagce(constd+4.)

      if ( qsz .ge. csmin ) then
         y1 = rhoz*qsz*0.001  !in CGS
         fv = sqrt( rhoe_s/rhoz )
         igce = 1
         if ( igce .ne. 1 ) then
            ! old codes from Chern's in MKS
            tmp1 = sqrt(pi*rhosnow*xnos/rhoz/qsz)
            tmp1 = sqrt(tmp1)
            vts = constc*gamdp4*fv/tmp1**constd
            vts = vts/6.
         else
            ! new codes from Steve's in cgs
!            y1 = rhoz * 0.001 *qsz             ! rhoz need to be in CGS
            y1 = qsz
            r00 = rhoz * 0.001   ! rhoz need to be in CGS
            vscf = vsc*fv

            ftns = 1.
            ftns0 = 1.
            if ( improve .eq. 3 )then
               tzc = tz - t0
               call sgmap(1,y1,r00,tzc,ftns0)
               ftns = ftns0**bsq
            endif
            vts = max(vscf*(r00*y1)**bsq/ftns, 0.0)
            vts = vts * 0.01  ! convert back to MKS
         endif  !end of if igce
         vts = min ( vsmax , max ( vsmin , vts ) )
      else
         vts = vsmin
      endif  !end of if qsz

      return
      end subroutine vts_mks
!-------------------------------------------------------------------

!-------------------------------------------------------------------
      subroutine vtg_mks(ihail,improve,rhoz,qgz,tz,vtg)
      implicit none
      integer, intent(in) :: ihail, improve
      real, intent(in) :: rhoz  !air density (kg/m^3)
      real, intent(in) :: qgz   !mixing ratio (kg/kg)
      real, intent(in) :: tz    !temperature (K)
      real, intent(out):: vtg   !terminal velocity (m/s)

      real, parameter :: vgmax = 8.0     !(m/s)
      real, parameter :: vgmin = 0.0     !(m/s)

      !local variables
      integer :: igce
      real    :: tmp1, pi, gam4pt5, term0, grav, r00, vgcr
      real    :: y1, fv, ftng, ftng0, tzc

      pi = acos(-1.)
      grav = 9.80665e+0
      gam4pt5 = gammagce(4.5)

      if ( qgz .ge. cgmin ) then
         if (ihail .eq. 1) then
            ! for hail, based on Lin et al (1983)
            tmp1 = sqrt(pi*rhohail*xnoh/rhoz/qgz)
            tmp1 = sqrt(tmp1)
            term0 = sqrt(4.*grav*rhohail/3./rhoz/cdrag)
            vtg = gam4pt5*term0*sqrt(1./tmp1)
            vtg = vtg/6.
         else
            ! for graupel, based on RH (1984)
            igce = 1
!            if ( igce .ne. 1 ) then
!               ! old codes from Chern's in MKS; prez(:) needed.
!               tmp1=sqrt(pi*rhograul*xnog/rhoz/qgz)
!               tmp1=sqrt(tmp1)
!               tmp1=tmp1**bbar
!               tmp1=1./tmp1
!               term0=abar*gam4bbar/6.
!               vtg=term0*tmp1*(p0/prez)**0.4
!            else
               ! new codes from Steve's in CGS
               y1 = qgz
               r00 = rhoz*0.001  !in CGS
               fv = sqrt( rhoe_s/rhoz )
               vgcr = vgc*fv
               ftng = 1.
               ftng0 = 1.
               if ( improve .gt. 2 ) then
                  tzc = tz - t0
                  call sgmap(2,y1,r00,tzc,ftng0)
                  ftng = ftng0**bgq
               endif
               vtg = dmax1(vgcr*(r00*y1)**bgq/ftng, 0.0)
               vtg = vtg * 0.01  !convert back to MKS

!               !>>> from GFDL MP :
!!               vcong = 87.2382675
!!               rhof = sqrt(min(10.,1.2/rhoz))
!!               normg = 5026548245.74367
!!               vtg = vcong * rhof * sqrt(sqrt(sqrt(qgz*rhoz/normg)))
!               vtg = 87.2382675 * sqrt(min(10.,rhoe_s/rhoz)) * &
!                     sqrt(sqrt(sqrt(qgz*rhoz/5026548245.74367)))
!               !<<<
!            endif  !end of if igce
            vtg = min ( vgmax , max ( vgmin , vtg ) )
         endif  !end of if ihail
      else
         vtg = vgmin
      endif   !end of if qgz

      return
      end subroutine vtg_mks
!-------------------------------------------------------------------

!-------------------------------------------------------------------
      SUBROUTINE semi_sedi(qvar,ihail,improve,iter,km,dzl,rho,qc,qv,p,tz,xland,ww,precip,dt)
!-------------------------------------------------------------------
!
! This routine is a semi-Lagrangain forward advection for hydrometeors
! with mass conservation and positive definite advection
! 2nd order interpolation with monotonic piecewise parabolic method is used.
! This routine is under assumption of decfl < 1 for semi_Lagrangian
!
! dzl    depth of model layer (m)
! rho    dry air density (kg/m^3)
! qc     dry mixing ratio of condensate (kg/kg)
! qv     specific humidity (kg/kg)
! p      pressure (Pa)
! tz     air temperature (K)
! ww     terminal velocity (m/s)
! precip total precipitation at surface (mm)
! dt     time step (sec)
! iter   how many time to guess mean terminal velocity: 0 pure forward.
!        0 : use departure wind for advection
!        1 : use mean wind for advection
!        > 1 : use mean wind after iter-1 iterations
!
! author: hann-ming henry juang <henry.juang@noaa.gov>
!         implemented by song-you hong
! reference: Juang, H.-M., and S.-Y. Hong, 2010: Forward semi-Lagrangian advection
!         with mass conservation and positive definiteness for falling
!         hydrometeors. *Mon.  Wea. Rev.*, *138*, 1778-1791
!
      implicit none

      character(len=2) :: qvar
      integer, intent(in) :: ihail, improve, km, iter
      real, intent(in) ::  dt
      real, intent(in) :: dzl(km), rho(km), tz(km), qv(km), p(km)
      real, intent(in) :: xland
      real, intent(out) :: ww(km)
      real, intent(out) :: precip
      real, intent(inout) :: qc(km)
      integer  k,m,kk,kb,kt,n
      real  tl,tl2,qql,dql,qqd
      real  th,th2,qqh,dqh
      real  zsum,qsum,dim,dip,con1,fa1,fa2
      real  allold, decfl
      real  dz(km), qq(km)
      real  wi(km+1), zi(km+1), za(km+2)
      real  qn(km)
      real  dza(km+1), qa(km+1), qmi(km+1), qpi(km+1)
      real  wd(km),wa(km),was(km)
!
      precip = 0.0
      qa(:) = 0.0
      qq(:) = 0.0
      ww(:) = 0.0
      wa(:) = 0.0
      was(:) = 0.0
      dz(:) = dzl(:)
      do k = 1,km
        if ( qvar .eq. 'qr' ) call vtr_mks(rho(k),qc(k),tz(k),ww(k))
        if ( qvar .eq. 'qs' ) call vts_mks(improve,rho(k),qc(k),tz(k),ww(k))
        if ( qvar .eq. 'qg' ) call vtg_mks(ihail,improve,rho(k),qc(k),tz(k),ww(k))
        if ( qvar .eq. 'qi' ) call vti_mks(improve,rho(k),tz(k),qc(k),qv(k),p(k),xland,ww(k))
      enddo
      do k = 1,km
        qq(k) = qc(k)*rho(k)
      enddo

! skip for no precipitation for all layers
      allold = 0.0
      do k = 1,km
        allold = allold + qq(k)
      enddo
      if ( allold .le. 0.0 ) then
         return 
      endif
!
! compute interface values
      zi(1) = 0.0
      do k = 1,km
        zi(k+1) = zi(k)+dz(k)
      enddo
!====================================
! save departure wind
      wd(:) = ww(:)
      n = 1
 100  continue
! plm is 2nd order, we can use 2nd order wi or 3rd order wi
! 2nd order interpolation to get wi
      wi(1) = ww(1)
      wi(km+1) = ww(km)
      do k = 2,km
        wi(k) = (ww(k)*dz(k-1)+ww(k-1)*dz(k))/(dz(k-1)+dz(k))
      enddo
! 3rd order interpolation to get wi
      fa1 = 9./16.
      fa2 = 1./16.
      wi(1) = ww(1)
      wi(2) = 0.5*(ww(2)+ww(1))
      do k = 3,km-1
        wi(k) = fa1*(ww(k)+ww(k-1))-fa2*(ww(k+1)+ww(k-2))
      enddo
      wi(km) = 0.5*(ww(km)+ww(km-1))
      wi(km+1) = ww(km)

! terminate of top of raingroup
      do k = 2,km
        if( ww(k).eq.0.0 ) wi(k) = ww(k-1)
      enddo

! diffusivity of wi
      con1 = 0.05
      do k = km,1,-1
        decfl = (wi(k+1)-wi(k))*dt/dz(k)
        if( decfl .gt. con1 ) then
          wi(k) = wi(k+1) - con1*dz(k)/dt
        endif
      enddo
! compute arrival point
      do k = 1,km+1
        za(k) = zi(k) - wi(k)*dt
      enddo
      za(km+2) = zi(km+1)

      do k = 1,km+1
        dza(k) = za(k+1) - za(k)
      enddo

! computer deformation at arrival point
      do k = 1,km
        qa(k) = qq(k)*dz(k)/dza(k)
        qc(k) = qa(k)/rho(k)
      enddo
      qa(km+1) = 0.0

! compute arrival terminal velocity, and estimate mean terminal velocity
! then back to use mean terminal velocity
      if ( n.le.iter ) then
        do k = 1,km
          if ( qvar .eq. 'qr' ) call vtr_mks(rho(k),qc(k),tz(k),wa(k))
          if ( qvar .eq. 'qs' ) call vts_mks(improve,rho(k),qc(k),tz(k),wa(k))
          if ( qvar .eq. 'qg' ) call vtg_mks(ihail,improve,rho(k),qc(k),tz(k),wa(k))
          if ( qvar .eq. 'qi' ) call vti_mks(improve,rho(k),tz(k),qc(k),qv(k),p(k),xland,wa(k))
        enddo

        do k = 1,km
           if ( n.ge.2 ) wa(k) = 0.5*(wa(k)+was(k))
           ! mean wind is average of departure and new arrival winds
           ww(k) = 0.5*(wd(k)+wa(k))
        enddo
        was(:) = wa(:)
        n = n + 1
        go to 100
      endif
!====================================
! estimate values at arrival cell interface with monotone
      do k = 2,km
        dip = (qa(k+1)-qa(k))/(dza(k+1)+dza(k))
        dim = (qa(k)-qa(k-1))/(dza(k-1)+dza(k))
        if ( dip*dim.le.0.0 ) then
          qmi(k) = qa(k)
          qpi(k) = qa(k)
        else
          qpi(k) = qa(k)+0.5*(dip+dim)*dza(k)
          qmi(k) = 2.0*qa(k)-qpi(k)
          if ( qpi(k).lt.0.0 .or. qmi(k).lt.0.0 ) then
            qpi(k) = qa(k)
            qmi(k) = qa(k)
          endif
        endif
      enddo
      qpi(1) = qa(1)
      qmi(1) = qa(1)
      qmi(km+1) = qa(km+1)
      qpi(km+1) = qa(km+1)

! interpolation to regular point
      qn = 0.0
      kb = 1
      kt = 1
      intp : do k = 1,km
             kb = max(kb-1,1)
             kt = max(kt-1,1)
! find kb and kt
             if ( zi(k).ge.za(km+1) ) then
               exit intp
             else
               find_kb : do kk = kb,km
                         if ( zi(k).le.za(kk+1) ) then
                           kb = kk
                           exit find_kb
                         else
                           cycle find_kb
                         endif
               enddo find_kb
               find_kt : do kk = kt,km+2
                         if( zi(k+1).le.za(kk) ) then
                           kt = kk
                           exit find_kt
                         else
                           cycle find_kt
                         endif
               enddo find_kt
               kt = kt - 1
! compute q with piecewise constant method
               if ( kt.eq.kb ) then
                 tl = (zi(k)-za(kb))/dza(kb)
                 th = (zi(k+1)-za(kb))/dza(kb)
                 tl2 = tl*tl
                 th2 = th*th
                 qqd = 0.5*(qpi(kb)-qmi(kb))
                 qqh = qqd*th2+qmi(kb)*th
                 qql = qqd*tl2+qmi(kb)*tl
                 qn(k) = (qqh-qql)/(th-tl)
               else if ( kt.gt.kb ) then
                 tl = (zi(k)-za(kb))/dza(kb)
                 tl2 = tl*tl
                 qqd = 0.5*(qpi(kb)-qmi(kb))
                 qql = qqd*tl2+qmi(kb)*tl
                 dql = qa(kb)-qql
                 zsum  = (1.-tl)*dza(kb)
                 qsum  = dql*dza(kb)
                 if ( kt-kb.gt.1 ) then
                   do m = kb+1,kt-1
                     zsum = zsum + dza(m)
                     qsum = qsum + qa(m) * dza(m)
                   enddo
                 endif
                 th = (zi(k+1)-za(kt))/dza(kt)
                 th2 = th*th
                 qqd = 0.5*(qpi(kt)-qmi(kt))
                 dqh = qqd*th2+qmi(kt)*th
                 zsum  = zsum + th*dza(kt)
                 qsum  = qsum + dqh*dza(kt)
                 qn(k) = qsum/zsum
               endif
               cycle intp
             endif
       enddo intp

! rain out (unit:kg/m^2=mm)
       sum_precip: do k = 1,km
          if ( za(k).lt.0.0 .and. za(k+1).le.0.0 ) then
             precip = precip + qa(k)*dza(k)
             cycle sum_precip
          else if ( za(k).lt.0.0 .and. za(k+1).gt.0.0 ) then
             ! from Thompson MP :
             th = (0.0-za(k))/dza(k)
             th2 = th*th
             qqd = 0.5*(qpi(k)-qmi(k))
             qqh = qqd*th2+qmi(k)*th
             precip = precip + qqh*dza(k)
!             ! from WSM6 MP : 
!             precip = precip + qa(k)*(0.-za(k))
             exit sum_precip
          endif
          exit sum_precip
       enddo sum_precip

! replace the new values
       do k = 1,km
          qc(k) = qn(k)/rho(k)
       enddo

       END SUBROUTINE semi_sedi

       REAL FUNCTION f_qsi(tair,p0)
       IMPLICIT NONE
       REAL :: tair   !real temperature (K)
       REAL :: p0     !pressure (Pa)
       REAL :: rp0

       rp0 = 3.799052e3 / (p0 * 10.)
       f_qsi = rp0 * exp(c218 - c580 / (tair - c76) )

       END FUNCTION f_qsi

#ifdef new_saturation
       REAL FUNCTION esi_mks(tair)
       IMPLICIT NONE
       REAL :: tair   !real temperature (K)
!
!  COMPUTE SATURATION VAPOR PRESSURE POLYSVP RETURNED IN UNITS OF PA. T IS INPUT IN UNITS OF K.
!  REPLACE GOFF-GRATCH WITH FASTER FORMULATION FROM FLATAU ET AL. 1992, TABLE 4 (RIGHT-HAND COLUMN)
       REAL :: DT
       REAL :: a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i
       DATA a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i /6.11147274,0.503160820, &
            0.188439774E-1,0.420895665E-3,0.615021634E-5,0.602588177E-7, &
            0.385852041E-9,0.146898966E-11,0.252751365E-14/

       DT = MAX(-80.,tair - 273.16)
       esi_mks = a0i + DT*(a1i + DT*(a2i + DT*(a3i + DT*(a4i + DT*(a5i + &
                 DT*(a6i + DT*(a7i + a8i*DT)))))))
       esi_mks = esi_mks*100.  !convert to Pa
!  
!  Goff-Gratch equation (Goff and Gratch 1945)
!       esi_mks = c610 * exp( c218 - c580 / (tair - c76) )
!       esi_mks = esi_mks / 10.  !convert to Pa

       END FUNCTION esi_mks

       REAL FUNCTION esw_mks(tair)
       IMPLICIT NONE
       REAL :: tair   !real temperature (K)
!
!  COMPUTE SATURATION VAPOR PRESSURE POLYSVP RETURNED IN UNITS OF PA. T IS INPUT IN UNITS OF K.
!  REPLACE GOFF-GRATCH WITH FASTER FORMULATION FROM FLATAU ET AL. 1992, TABLE 4 (RIGHT-HAND COLUMN)
       REAL :: DT
       REAL :: a0,a1,a2,a3,a4,a5,a6,a7,a8
       DATA a0,a1,a2,a3,a4,a5,a6,a7,a8 /6.11239921,0.443987641,          &
            0.142986287E-1,0.264847430E-3,0.302950461E-5,0.206739458E-7, &
            0.640689451E-10,-0.952447341E-13,-0.976195544E-15/

       DT = MAX(-80.,tair - 273.16)
       esw_mks = a0 + DT*(a1 + DT*(a2 + DT*(a3 + DT*(a4 + DT*(a5 + DT*   &
                 (a6 + DT*(a7 + a8*DT)))))))
       esw_mks = esw_mks*100.  !convert to Pa
!
!  Goff-Gratch equation (Goff and Gratch 1945)
!       esw_mks = c610 * exp( c172 - c409 / (tair - c358) )
!       esw_mks = esw_mks / 10.   !convert to Pa

       END FUNCTION esw_mks
#endif
END MODULE  module_mp_gsfcgce_3ice_nuwrf

