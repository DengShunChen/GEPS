!#define new_saturation
!#define sat_predict
!#define use_declination
!#define use_cpm

MODULE module_mp_gsfcgce_3ice_nuwrf_gpu
   !use nvtx
   use rank, only: myrank
#ifdef Readaeroclx
   USE module_gocart_coupling_gpu, only: mass2ccn_gpu, mass2icn_gpu, &
                    nlut, nsuso, nsoot, ninso, nwaso, nssam, nsscm, &
                    nminm, nmiam, nmicm, lut_in_025, lut_ccn, pts_t, pts_s, &
                    mxpts_t, mxpts_s
   USE const, only: naso4, nadu1, nadu2, nadu3, nadu4, nadu5, &
                    nass1, nass2, nass3, nass4, nass5, nablc, &
                    nabbc, naolc, naobc, namsa, nadms, naso2
#endif
   USE module_mp_radar
   use param, only: my, my_max

   PRIVATE   ! privatize all variables/subroutines in this module excepting public parameter below
   PUBLIC :: gsfcgce_3ice_nuwrf_gpu
   PUBLIC :: consat_s_gpu, fall_flux_gpu, saticel_s_gpu, semi_sedi_gpu, vti_mks_gpu, &
             vtr_mks_gpu, vts_mks_gpu, vtg_mks_gpu, sgmap_gpu, gammagce_gpu, &
             gamma_toshi_gpu, eff_rad_gpu, f_qsi_gpu, esi_mks_gpu, esw_mks_gpu, &
             gce_table_copyin_gpu, gce_table_delete_gpu

   REAL, PRIVATE ::          rd1, rd2, al, cp

   REAL, PRIVATE ::          c38, c358, c610, c149, &
                             c879, c172, c409, c76, &
                             c218, c580, c141
   REAL, PRIVATE ::          ag, bg, as, bs, &
                             aw, bw, bgh, bgq, &
                             bsh, bsq, bwh, bwq

   REAL, PRIVATE ::          tnw, tns, tng, &
                             roqs, roqg, roqr

   REAL, PRIVATE ::          zrc, zgc, zsc, vrc, &
                             vrc0, vrc1, vrc2, vrc3, &
                             vgc, vsc

   REAL, PRIVATE ::          draimax

   REAL, PRIVATE ::          bnd1, rn11a

   REAL, PRIVATE ::          rn17, rn19b, &
                             bnd3, rn23a, &
                             rn23b, rn30b, &
                             rn30c

   REAL, PRIVATE ::          alv, alf, als, t0, t00, &
                             avc, afc, asc, rn1, rn2, &
                             bnd2, rn3, rn4, rn5, rn50, &
                             rn51, rn52, rn53, rn6, rn60, &
                             rn61, rn62, rn63, rn7, rn8, &
                             rn9, rn10, rn101, rn102, rn10a, &
                             rn10b, rn10c, rn11, rn12, rn14, &
                             rn15, rn15a, rn16, rn171, rn172, &
                             rn17a, rn17b, rn17c, rn18, rn18a, &
                             rn19, rn191, rn192, rn19a, rn20, &
                             rn20a, rn20b, rn30, rn30a, rn21, &
                             bnd21, rn22, rn23, rn231, rn232, &
                             rn25, rn31, beta, rn32, rn33, &
                             rn331, rn332, rn34, rn35, rnn30a, &
                             rnn191, rnn192, cn0

   REAL, PRIVATE ::          ami50, ami40, ami100

   REAL, PRIVATE, DIMENSION(31) ::    rn12a, rn12b, rn13, rn25a

   REAL, PRIVATE, DIMENSION(31) ::    BergCon1, BergCon2, &
                                      BergCon3, BergCon4

   ! critical q of hydrometeor characteristics (eg. fall speed, radius)
   REAL, PRIVATE :: cwmin, cimin, crmin, csmin, cgmin

   REAL, PRIVATE, DIMENSION(31)  ::      aa1, aa2
   DATA aa1/.7939e-7, .7841e-6, .3369e-5, .4336e-5, .5285e-5, &
      .3728e-5, .1852e-5, .2991e-6, .4248e-6, .7434e-6, &
      .1812e-5, .4394e-5, .9145e-5, .1725e-4, .3348e-4, &
      .1725e-4, .9175e-5, .4412e-5, .2252e-5, .9115e-6, &
      .4876e-6, .3473e-6, .4758e-6, .6306e-6, .8573e-6, &
      .7868e-6, .7192e-6, .6513e-6, .5956e-6, .5333e-6, &
      .4834e-6/
   DATA aa2/.4006, .4831, .5320, .5307, .5319, &
      .5249, .4888, .3894, .4047, .4318, &
      .4771, .5183, .5463, .5651, .5813, &
      .5655, .5478, .5203, .4906, .4447, &
      .4126, .3960, .4149, .4320, .4506, &
      .4483, .4460, .4433, .4413, .4382, &
      .4361/

   real, dimension(0:120) :: itble       ! deposition growth coefficients
   data itble/1.000000, 0.979490, 0.959401, 0.939723, 0.920450, &
      0.899498, 0.879023, 0.857038, 0.833681, 0.810961, &
      0.783430, 0.755092, 0.703072, 0.537032, 0.467735, &
      0.524807, 0.630957, 0.812831, 1.096478, 1.479108, &
      1.905461, 2.089296, 2.290868, 2.398833, 2.454709, &
      2.426610, 2.371374, 2.290868, 2.137962, 1.995262, &
      1.862087, 1.737801, 1.621810, 1.513561, 1.396368, &
      1.288250, 1.188502, 1.096478, 1.000000, 0.922571, &
      0.851138, 0.785236, 0.724436, 0.668344, 0.616595, &
      0.575440, 0.537032, 0.501187, 0.467735, 0.436516, &
      0.407380, 0.380189, 0.354813, 0.331131, 0.316228, &
      0.301995, 0.291743, 0.285102, 0.281838, 0.278612, &
      0.275423, 0.278612, 0.281838, 0.285102, 0.291743, &
      0.298538, 0.309030, 0.319890, 0.331131, 0.346737, &
      0.367282, 0.393550, 0.426580, 0.457088, 0.489779, &
      0.524807, 0.562341, 0.609537, 0.660693, 0.716143, &
      0.785236, 0.860994, 0.954993, 1.047129, 1.148154, &
      1.258925, 1.380384, 1.496236, 1.603245, 1.698244, &
      1.778279, 1.840772, 1.883649, 1.905461, 1.905461, &
      1.883649, 1.862087, 1.840772, 1.798871, 1.737801, &
      1.698244, 1.640590, 1.584893, 1.548817, 1.513561, &
      1.475707, 1.452112, 1.428894, 1.412538, 1.393157, &
      1.377209, 1.361445, 1.348963, 1.336596, 1.327394, &
      1.318257, 1.309182, 1.303167, 1.294196, 1.288250, 1.279381/
   real, private, parameter :: thrd = 1./3.

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


   subroutine gce_table_copyin_gpu
   implicit none
   integer :: async_id = -1


   !$acc enter data copyin(itble, aa2, rn12a, rn12b, rn13, aa1) async(async_id)
   !$acc wait(async_id)


   end subroutine
   
   subroutine gce_table_delete_gpu
   implicit none
   integer :: async_id = -1


   !$acc exit data delete(itble, aa2, rn12a, rn12b, rn13, aa1) async(async_id)
   !$acc wait(async_id)


   end subroutine

!-------------------------------------------------------------------
!  NASA/GSFC GCE
!  Tao et al, 2001, Meteo. & Atmos. Phy., 97-137
!-------------------------------------------------------------------
   SUBROUTINE gsfcgce_3ice_nuwrf_gpu(myim, th &
                                 , qv, ql &
                                 , qr, qi &
                                 , qs &
                                 , rho, pii, p, dt_in, z &
                                 , ht, dz8w, grav, w &
                                 !                      ,rhowater, rhosnow                           &
                                 , itimestep, xlat, sdec, xland &
                                 !                      ,ids,ide, jds,jde, kds,kde                   & ! domain dims
                                 , ims, ime, jms, jme, kms, kme & ! memory dims
                                 , its, ite, jts, jte, kts, kte & ! tile   dims
                                 , rainnc, icenc, snownc, graupelnc, sr &
                                 !                      ,rainncv, snowncv, graupelncv                &
                                 , f_qg, qg &
#ifdef Readaeroclx
                                 , aeroclx, naero &
#endif
                                 , ihail, ice2 &
#ifdef EXT_DIAG
                                 , refl_10cm, diagflag, do_radar_ref &
                                 !                      ,physc, physe, physd, physs, physm, physf    &
                                 !                      ,acphysc, acphyse, acphysd, acphyss, acphysm, acphysf &
                                 , preci3d, precs3d, precg3d, precr3d &
#endif
                                 , refc, refr, refi, refs, refg & ! cloud effective radius
                                 , SL_sedi, sat_predict, new_saturation, use_cpm, use_declination, &
                                 benchmark, async_id) !flags

!-------------------------------------------------------------------
      IMPLICIT NONE
!-------------------------------------------------------------------
      INTEGER, INTENT(IN)    ::   ims, ime, jms, jme, kms, kme, &
                                  its, ite, jts, jte, kts, kte
      INTEGER, INTENT(IN)    ::   itimestep, ihail, ice2

      integer, intent(in) :: myim(my_max)

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(INOUT) :: th, qv, ql, qr, qi, qs, qg
!
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(IN) :: rho, pii, p, dz8w, z, w

      REAL, DIMENSION(ims:ime, jms:jme), &
         INTENT(INOUT) ::    rainnc, icenc, snownc, graupelnc, sr
!  REAL, DIMENSION( ims:ime , jms:jme ),                           &
!        INTENT(INOUT) ::           rainncv, snowncv, graupelncv

!JJS 20140225   for calculation of effective radius of cloud species
      REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN)   :: XLAND
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(INOUT) :: refc, refr, refi, refs, refg

!+---+-----------------------------------------------------------------+
#ifdef EXT_DIAG
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT):: &  ! GT
         refl_10cm
      LOGICAL, OPTIONAL, INTENT(IN) :: diagflag
      INTEGER, OPTIONAL, INTENT(IN) :: do_radar_ref
#endif
!+---+-----------------------------------------------------------------+

      REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) ::       ht

      REAL, INTENT(IN) ::                                   dt_in, grav
!                                                        rhowater, &
!                                                         rhosnow
      REAL, INTENT(IN) :: xlat(my_max)
      REAL, INTENT(IN) :: sdec  ! sine of solar declination angle
      LOGICAL, INTENT(IN), OPTIONAL :: F_QG
#ifdef Readaeroclx
      ! aerosol climatology
      integer, intent(in) :: naero
      real, dimension(ims:ime, kms:kme, naero, jms:jme), intent(in) :: aeroclx
#endif

!flags
      logical, intent(in) :: SL_sedi, sat_predict, new_saturation, use_cpm, &
                             use_declination

!  LOCAL VAR
      INTEGER ::  itaobraun, istatmin, new_ice_sat, id
      INTEGER ::  improve

      INTEGER :: i, j, k
      INTEGER :: iskip, ih, icount, ibud, i24h
      REAL    :: hour
      REAL, PARAMETER :: cmin = 1.e-20
      REAL    :: dth, dqv, dqrest, dqall, dqall1, rhotot, a1, a2

      LOGICAL :: flag_qg

#ifdef EXT_DIAG
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(INOUT) :: &
         !             physc, physe, physd, physs, physm, physf,                  &
         !             acphysc, acphyse, acphysd, acphyss, acphysm, acphysf,      &
         preci3d, precs3d, precg3d, precr3d
#endif

      ! for the conversion of q :
      real, dimension(its:ite, kts:kte, jts:jte) :: qtot

      ! for sub-cycle time steps :
      integer :: n, ntimes
      real :: mp_time, dts
      
      character(len=4) :: myrank_str
      
      ! for GPU porting
      integer, intent(in) :: async_id
      logical, intent(in) :: benchmark
      real :: qtotr, qvr, qlr, qrr, qir, qsr, qgr



!+---+-----------------------------------------------------------------+

      INTEGER :: NCALL = 0

!+---+-----------------------------------------------------------------+
     !if (benchmark .eq. .true.) call nvtxStartRange("GPU:GCE_3ice")

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

      i24h = nint(86400./dt_in)

!c  new_ice_sat = 0, 1, 2, or 3
      new_ice_sat = 3 ! Set to 3 in NUWRF

      istatmin = 180

!c id = 0  without in-line staticstics
!c id = 1  with in-line staticstics
      id = 0

!c ibud = 0 no calculation of dth, dqv, dqrest and dqall
!c ibud = 1 yes
      ibud = 0

      ! convert specific values of q to mixing ratios :
      !$acc data create(qtot) async(async_id)
      !qtot = 0.
      !$acc parallel loop gang collapse(2) async(async_id) &
      !$acc&         private(qtotr, qvr, qlr, qrr, qir, qsr, qgr)
      do j = jts, jte
         do k = kts, kte
            !$acc loop vector
            do i = its, myim(j)
               qvr = qv(i, k, j)
               qlr = ql(i, k, j)
               qrr = qr(i, k, j)
               qir = qi(i, k, j)
               qsr = qs(i, k, j)
               qgr = qg(i, k, j)
               qtotr = qvr + qlr + qrr + qir + qsr + qgr
               qvr = qvr/(1.-qtotr)
               qlr = qlr/(1.-qtotr)
               qrr = qrr/(1.-qtotr)
               qir = qir/(1.-qtotr)
               qsr = qsr/(1.-qtotr)
               qgr = qgr/(1.-qtotr)
               qv(i, k, j) = qvr
               ql(i, k, j) = qlr
               qr(i, k, j) = qrr
               qi(i, k, j) = qir
               qs(i, k, j) = qsr
               qg(i, k, j) = qgr

               ! qtot must be saved as "mixing ratio" :
               qtot(i, k, j) = qvr + qlr + qrr + qir + qsr + qgr
            end do
         end do
      end do

      ! set up constants used internally in GCE
      if (itimestep .eq. 1) then
         call consat_s_gpu(ihail, itaobraun, improve)
         call gce_table_copyin_gpu
      end if

      ! set sub-cycle time step :
      mp_time = 300.                !standard sub-cycle time step
      ntimes = 1                    !number of sub-cycles
      ntimes = max(ntimes, int(dt_in/min(dt_in, mp_time)))
      dts = dt_in/real(ntimes)   !real sub-cycle time step

      do n = 1, ntimes
            ! calculte fallflux and precipiation in MKS system
            call fall_flux_gpu(myim, dts, qv, qr, qi, &
                           qs, qg, p, &
                           rho, th, pii, &
                           z, dz8w, ht, &
                           grav, itimestep, &
                           rainnc, icenc, snownc, graupelnc, sr, &
                           !                      rainncv, snowncv, graupelncv,           &
   #ifdef EXT_DIAG
                           preci3d, precs3d, &
                           precg3d, precr3d, &
   #endif
                           xland, &
                           ihail, ice2, improve, &
                           ims, ime, jms, jme, kms, kme, & ! memory dims
                           its, ite, jts, jte, kts, kte, & ! tile   dims
                           SL_sedi, sat_predict, new_saturation, benchmark, async_id) !flags

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
               end if
               xbm_g = 3.
               xmu_g = 0.

               call radar_init
               NCALL = 1
            END IF
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
               end if
               xbm_g = 3.
               xmu_g = 0.

               call radar_init
               NCALL = 1
            END IF
   #endif

            ! microphysics in GCE
            call SATICEL_S_gpu(myim, dts, IHAIL, itaobraun, ICE2, istatmin, &
                           new_ice_sat, id, improve, xlat, sdec, &
                           th, qv, ql, qr, &
                           qi, qs, qg, &
                           rho, pii, p, w, &
                           itimestep, xland, &
                           refc, refr, refi, &
                           refs, refg, & ! cloud effective radius
                           ims, ime, jms, jme, kms, kme, & ! memory dims
                           its, ite, jts, jte, kts, kte, & ! tile   dims
   #ifdef Readaeroclx
                           aeroclx, naero, &
   #endif
   #ifdef EXT_DIAG
                           refl_10cm, diagflag, do_radar_ref, & ! GT added for reflectivity calcs
                           physc, physe, physd, physs, physm, physf, &
                           acphysc, acphyse, acphysd, acphyss, acphysm, acphysf, &
   #endif
                           sat_predict, new_saturation, use_cpm, use_declination, &
                           benchmark, async_id) !flags
      end do  !end of do n=1,ntimes

      ! convert mixing values of q back to specific values :
      !$acc parallel loop gang collapse(2) async(async_id) &
      !$acc&         private(qtotr)
      do j = jts, jte
         do k = kts, kte
            !$acc loop vector
            do i = its, myim(j)
               qtotr = qtot(i, k, j)
               qv(i, k, j) = qv(i, k, j)/(1.+qtotr)
               ql(i, k, j) = ql(i, k, j)/(1.+qtotr)
               qr(i, k, j) = qr(i, k, j)/(1.+qtotr)
               qi(i, k, j) = qi(i, k, j)/(1.+qtotr)
               qs(i, k, j) = qs(i, k, j)/(1.+qtotr)
               qg(i, k, j) = qg(i, k, j)/(1.+qtotr)
            end do
         end do
      end do
      !$acc end data
      !if (benchmark .eq. .true.) call nvtxEndRange

   END SUBROUTINE gsfcgce_3ice_nuwrf_gpu

   SUBROUTINE fall_flux_gpu(myim, dt, qv, qr, qi, qs, qg, p, &
                        rho, th, pi_mks, z, dz8w, topo, &
                        grav, itimestep, &
                        rainnc, icenc, snownc, graupelnc, sr, &
                        !                      rainncv, snowncv, graupelncv,           &
#ifdef EXT_DIAG
                        preci3d, precs3d, precg3d, precr3d, &
#endif
                        xland, &
                        ihail, ice2, improve, &
                        ims, ime, jms, jme, kms, kme, & ! memory dims
                        its, ite, jts, jte, kts, kte, & ! tile   dims
                        SL_sedi, sat_predict, new_saturation, benchmark, async_id) !flags
      !$acc routine(vti_mks_gpu) seq
      !$acc routine(vtg_mks_gpu) seq
      !$acc routine(vts_mks_gpu) seq
      !$acc routine(vtr_mks_gpu) seq

!-----------------------------------------------------------------------
! adopted from Jiun-Dar Chern's codes for Purdue Regional Model
! adopted by Jainn J. Shi, 6/10/2005
! modified by Goddard 7/24/2010
! modified by Tao 11/12/2010
!-----------------------------------------------------------------------

      IMPLICIT NONE
      integer, intent(in)               :: myim(my_max)
      INTEGER, INTENT(IN)               :: ihail, ice2, improve, &
                                           ims, ime, jms, jme, kms, kme, &
                                           its, ite, jts, jte, kts, kte
      INTEGER, INTENT(IN)               :: itimestep
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(INOUT)               :: qv, qr, qi, qs, qg
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(IN)                  :: th, pi_mks

      REAL, DIMENSION(ims:ime, jms:jme), &
         INTENT(INOUT) :: rainnc, icenc, snownc, graupelnc, sr
!  REAL,    DIMENSION( ims:ime , jms:jme ),                            &
!           INTENT(INOUT)               :: rainncv, snowncv, graupelncv
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(IN)               :: rho, z, dz8w, p

      REAL, INTENT(IN)               :: dt, grav

      REAL, DIMENSION(ims:ime, jms:jme), &
         INTENT(IN)               :: topo, xland

#ifdef EXT_DIAG
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(OUT)                 :: preci3d, precs3d, precg3d, precr3d
#endif

!flags
      logical, intent(in) :: SL_sedi, sat_predict, new_saturation


! temperary vars
      REAL, DIMENSION(its:ite, kts:kte, jts:jte) :: z1
      REAL, DIMENSION(its:ite, kts:kte, jts:jte) :: fv
      REAL                                       :: tmp1, term0
      REAL, DIMENSION(its:ite, jts:jte)          :: pptrain, pptsnow, &
                                                    pptgraul, pptice, pptall
      REAL, DIMENSION(its:ite, kts:kte, jts:jte) :: orhoz
      !REAL, DIMENSION(its:ite, kts:kte, jts:jte) :: rsed, ised, ssed, gsed
      REAL, DIMENSION(its:ite, kts:kte, jts:jte) :: tz

      INTEGER                                    :: k, i, j
!

      REAL, DIMENSION(its:ite, kts:kte, jts:jte) :: vtr, vts, vtg, vti

      REAL                          :: dtb, pi, gambp4, gamdp4, gam4pt5, gam4bbar

!NUWRF BEGIN
! New local variable
      INTEGER                       :: ic, igce
      REAL                          :: y1, vr, vs, vg
      REAL                          :: tslopes, tslopeg, vgcr, vscf
      REAL                          :: tair, tairc, fexp
      REAL                          :: ftns, ftnsQ, ftng, ftngQ
      REAL                          :: const_vt, const_d, const_m, bb1, bb2
      REAL, DIMENSION(7)            :: aice, vice
      REAL                          :: cmin, r00
      REAL                          :: ftns0, ftng0

      DATA tslopes/0./, tslopeg/0./
      DATA aice/1.e-6, 1.e-5, 1.e-4, 1.e-3, 0.01, 0.1, 1./  ! g/m**3
      DATA vice/5., 15., 30., 35., 40., 45., 50./  ! cm/s
      DATA igce/0/
      DATA ftns/1./, ftng/1./

!NUWRF END

      REAL, PARAMETER :: p0 = 1.0e5

! for terminal velocity flux
      INTEGER                          :: min_q, max_q
      REAL                             :: t_del_tv, del_tv, flux, fluxin, fluxout, tmpqrz
      LOGICAL                          :: notlast
      integer                          :: n, ntimes          !SL_sedi
      real                             :: dtcfl      !SL_sedi
      real, allocatable, dimension(:)  :: qden         !SL_sedi
      !real, dimension(its:ite, kts:kte, jts:jte) :: vtr3d, vts3d, vtg3d, vti3d
      
      ! for GPU porting
      integer, intent(in) :: async_id
      integer :: async_id_qr, async_id_qs, async_id_qg, async_id_qi
      integer :: nhydro
      real    :: rhor, qsr, qir, qrr, qgr, vtrr, vtgr, vtsr, vtir, xlandr
      logical, intent(in) :: benchmark



      !if (SL_sedi) allocate(qden(kts:kte))
      !if (benchmark .eq. .true.) call nvtxStartRange("GPU:fall_flux")

      if (improve .eq. 3) igce = 1

!-----------------------------------------------------------------------
!  This program calculates precipitation fluxes due to terminal velocities.
!-----------------------------------------------------------------------

      dtb = dt
      pi = acos(-1.)
!  Gamma function
      gambp4 = gammagce_gpu(constb + 4.)   !ga4b
      gamdp4 = gammagce_gpu(constd + 4.)   !ga4d
      gam4pt5 = gammagce_gpu(4.5)
      gam4bbar = gammagce_gpu(4.+bbar)
!
      cmin = 1.e-10
!
!***********************************************************************
! Calculate precipitation fluxes due to terminal velocities.
!***********************************************************************
!
!- Calculate termianl velocity (vt?)  of precipitation q?z
!- Find maximum vt? to determine the small delta t

!!!!!!!!!!!! fall_flux_A_loop
      !$acc enter data create(fv, pptrain, pptsnow, pptgraul, pptice, &
      !$acc&      orhoz, tz) async(async_id)
!      !$acc parallel loop gang collapse(2) async(async_id)
!      do j = jts, jte
!         do k = kts, kte
!            !$acc loop vector
!            do i = its, myim(j)
!   #ifdef EXT_DIAG
!               preci3d(i, k, j) = 0.
!               precs3d(i, k, j) = 0.
!               precg3d(i, k, j) = 0.
!               precr3d(i, k, j) = 0.
!   #endif
!               ised(i, k, j) = 0.
!               ssed(i, k, j) = 0.
!               gsed(i, k, j) = 0.
!               rsed(i, k, j) = 0.
!            end do
!         end do
!      end do
      !$acc enter data create(vtr, vts, vtg, vti, z1, pptall) &
      !$acc&      async(async_id)

            ! in MKS system
      !$acc parallel loop gang collapse(2) async(async_id)
      do j = jts, jte
         do k = kts, kte
            !$acc loop vector
            do i = its, myim(j)
               orhoz(i, k, j) = 1./rho(i, k, j)
               fv(i, k, j) = sqrt(rhoe_s/rho(i, k, j))
   !      fv(k)=sqrt(rho(i,1,j)/rhoz(k))
               tz(i, k, j) = th(i, k, j)*pi_mks(i, k, j)
               if (k .ne. 1) z1(i, k, j) = 0.9*(z(i, k, j) - z(i, k-1, j))
            end do !k
         end do
      end do
      
      !$acc parallel loop gang vector collapse(2) async(async_id)
      do j = jts, jte
         do i = its, ite
            if (i .le. myim(j)) then
               pptrain(i, j) = 0.
               pptsnow(i, j) = 0.
               pptgraul(i, j) = 0.
               pptice(i, j) = 0.
               z1(i, 1, j) = 0.9*(z(i, 1, j) - topo(i, j))
            end if
         end do
      end do

      IF ((ice2 .eq. 0 .or. ice2 .eq. 2) .eq. .false.) THEN
         !$acc parallel loop gang collapse(2) async(async_id)
         do j = jts, jte
            do k = kts, kte
               !$acc loop vector
               do i = its, myim(j)
                  qg(i, k, j) = 0.
               END DO
            end do
         end do
      END IF
   !
   !-- rain
   !
      if (SL_sedi) then
         if (async_id .eq. -1) then
            async_id_qr = -1
            async_id_qs = -1
            async_id_qg = -1
            async_id_qi = -1
         else
            async_id_qr = 2
            async_id_qs = 3
            async_id_qg = 4
            async_id_qi = 5
            !$acc wait(async_id)
         end if
         ntimes = 1
         dtcfl = dtb/real(ntimes)

         do n = 1, ntimes
            call semi_sedi_gpu(myim, 'qr', its, ite, jts, jte, &
                ihail, improve, 1, kte, dz8w, rho, qr, qv, p, tz, xland, &
                vtr, pptrain, dtcfl, sat_predict, new_saturation, async_id_qr)
            call semi_sedi_gpu(myim, 'qs', its, ite, jts, jte, &
                 ihail, improve, 1, kte, dz8w, rho, qs, qv, p, tz, xland, &
                 vts, pptsnow, dtcfl, sat_predict, new_saturation, async_id_qs)
   !   ice2=0 --- with hail/graupel
   !   ice2=1 --- without hail/graupel
            if (ice2 .eq. 0) then
!-- If IHAIL=1, use hail.
!-- If IHAIL=0, use graupel.
               call semi_sedi_gpu(myim, 'qg', its, ite, jts, jte, &
                  ihail, improve, 1, kte, dz8w, rho, qg, qv, p, tz, xland, &
                  vtg, pptgraul, dtcfl, sat_predict, new_saturation, async_id_qg)
            end if
            call semi_sedi_gpu(myim, 'qi', its, ite, jts, jte, &
                  ihail, improve, 0, kte, dz8w, rho, qi, qv, p, tz, xland, &
                  vti, pptice, dtcfl, sat_predict, new_saturation, async_id_qi)
         end do
         if (async_id .ne. -1) then
            !$acc wait(async_id_qr)
            !$acc wait(async_id_qs)
            !$acc wait(async_id_qg)
            !$acc wait(async_id_qi)
         end if
      else !SL_sedi
         !$acc parallel loop gang vector collapse(3) async(async_id) &
         !$acc&         private(t_del_tv, del_tv, notlast, min_q, max_q, &
         !$acc&         fluxout, flux, fluxin, rhor, qsr, qir, qrr, qgr, &
         !$acc&         vtrr, vtsr, vtgr, vtir, xlandr)
         do nhydro = 1, 4
            do j = jts, jte
               do i = its, ite
                  if (i .le. myim(j)) then
                     if (nhydro .eq. 1) then ! rain
                        t_del_tv = 0.
                        del_tv = dtb
                        notlast = .true.
                        DO while (notlast)
               !
                           min_q = kte
                           max_q = kts - 1
               !

                           do k = kts, kte - 1

                              vtrr = 0.
                              qrr = qr(i, k, j)
                              if (qrr .gt. crmin) then
                                 min_q = min0(min_q, k)
                                 max_q = max0(max_q, k)

                                 call vtr_mks_gpu(rho(i, k, j), qrr, &
                                    tz(i, k, j), vtrr, constb, consta, &
                                    rhoe_s, crmin, rhowater, xnor, vrc0, vrc1, &
                                    vrc2, vrc3, draimax, tnw, roqr, t0)
               !           if (.not. vtr(k) .gt. 0.0) cycle ! EMK NUWRF Bug fix

                                 del_tv = dmin1(del_tv, z1(i, k, j)/vtrr)
                              end if
                              vtr(i, k, j) = vtrr
                           end do

                           if (max_q .ge. min_q) then
               !
               !- Check if the summation of the small delta t >=  big delta t
               !             (t_del_tv)          (del_tv)             (dtb)

                              t_del_tv = t_del_tv + del_tv
               !
                              if (t_del_tv .ge. dtb) then
                                 notlast = .false.
                                 del_tv = dtb + del_tv - t_del_tv
                              end if

               ! use small delta t to calculate the qrz flux
               ! termi is the qrz flux pass in the grid box through the upper boundary
               ! termo is the qrz flux pass out the grid box through the lower boundary
               !
                              fluxin = 0.
                              !$acc loop seq
                              do k = max_q, min_q, -1
                                 rhor = rho(i, k, j)
                                 qrr = qr(i, k, j)
                                 fluxout = rhor*vtr(i, k, j)*qrr
                                 flux = (fluxin - fluxout)/rhor/dz8w(i, k, j)
               !            tmpqrz=qrz(k)
                                 qrr = qrr + del_tv*flux
                                 qr(i, k, j) = dmax1(0., qrr)
                                 fluxin = fluxout
                                 !rsed(i, k, j) = rsed(i, k, j) + fluxin
                              end do
                              if (min_q .eq. 1) then
                                 pptrain(i, j) = pptrain(i, j) + fluxin*del_tv
                              else
                                 qr(i, min_q - 1, j) = qr(i, min_q - 1, j) &
                                    + del_tv*fluxin/rho(i, min_q - 1, j) &
                                    /dz8w(i, min_q - 1, j)
                              end if
               !
                           else
                              notlast = .false.
                           end if
                        END DO
            
      !
      !-- snow
      !
                     elseif (nhydro .eq. 2) then
                        t_del_tv = 0.
                        del_tv = dtb
                        notlast = .true.

                        DO while (notlast)
               !
                           min_q = kte
                           max_q = kts - 1

               !
                           do k = kts, kte - 1
                              vtsr = 0.
                              qsr = qs(i, k, j)

                              if (qsr .gt. csmin) then
                                 min_q = min0(min_q, k)
                                 max_q = max0(max_q, k)

                                 call vts_mks_gpu(improve, rho(i, k, j), &
                                    qsr, tz(i, k, j), vtsr, &
                                    rhoe_s, constd, csmin, rhosnow, xnos, &
                                    constc, vsc, t0, roqs, roqg, tns, tng, bsq)

                                 del_tv = dmin1(del_tv, z1(i, k, j)/vtsr)
                              end if
                              vts(i, k, j) = vtsr
                           end do

                           if (max_q .ge. min_q) then
               !
               !
               !- Check if the summation of the small delta t >=  big delta t
               !             (t_del_tv)          (del_tv)             (dtb)

                              t_del_tv = t_del_tv + del_tv

                              if (t_del_tv .ge. dtb) then
                                 notlast = .false.
                                 del_tv = dtb + del_tv - t_del_tv
                              end if

               ! use small delta t to calculate the qsz flux
               ! termi is the qsz flux pass in the grid box through the upper boundary
               ! termo is the qsz flux pass out the grid box through the lower boundary
               !
                              fluxin = 0.
                              do k = max_q, min_q, -1
                                 rhor = rho(i, k, j)
                                 qsr = qs(i, k, j)
                                 fluxout = rhor*vts(i, k, j)*qsr
                                 flux = (fluxin - fluxout)/rhor/dz8w(i, k, j)
                                 qsr = qsr + del_tv*flux
                                 qs(i, k, j) = dmax1(0., qsr)
                                 fluxin = fluxout
                                 !ssed(i, k, j) = ssed(i, k, j) + fluxin
                              end do
                              if (min_q .eq. 1) then
                                 pptsnow(i, j) = pptsnow(i, j) + fluxin*del_tv
                              else
                                 qs(i, min_q - 1, j) = qs(i, min_q - 1, j) &
                                    + del_tv*fluxin/rho(i, min_q - 1, j) &
                                    /dz8w(i, min_q - 1, j)
                              end if
               !
                           else
                              notlast = .false.
                           end if

                        END DO
            
      !
      !   ice2=0 --- with hail/graupel
      !   ice2=1 --- without hail/graupel
      !
                     elseif (nhydro .eq. 3) then
                        if (ice2 .eq. 0) then

      !
      !-- If IHAIL=1, use hail.
      !-- If IHAIL=0, use graupel.
      !
                           t_del_tv = 0.
                           del_tv = dtb
                           notlast = .true.
               !
                           DO while (notlast)
               !
                              min_q = kte
                              max_q = kts - 1
               !
                              do k = kts, kte - 1
                                 vtgr = 0.
                                 qgr = qg(i, k, j)
                                 if (qgr .gt. cgmin) then
                                    min_q = min0(min_q, k)
                                    max_q = max0(max_q, k)

                                    call vtg_mks_gpu(ihail, improve, &
                                       rho(i, k, j), qgr, tz(i, k, j), &
                                       vtgr, cgmin, rhohail, xnoh, &
                                       cdrag, rhoe_s, vgc, roqs, roqg, &
                                       tns, tng, bgq, t0)

                                    del_tv = dmin1(del_tv, z1(i, k, j)/vtgr)
               !
                                 end if !qgz
                                 vtg(i, k, j) = vtgr
                              end do !k

                              if (max_q .ge. min_q) then
               !
               !
               !- Check if the summation of the small delta t >=  big delta t
               !             (t_del_tv)          (del_tv)             (dtb)

                                 t_del_tv = t_del_tv + del_tv

                                 if (t_del_tv .ge. dtb) then
                                    notlast = .false.
                                    del_tv = dtb + del_tv - t_del_tv
                                 end if

               ! use small delta t to calculate the qgz flux
               ! termi is the qgz flux pass in the grid box through the upper boundary
               ! termo is the qgz flux pass out the grid box through the lower boundary
               !
                                 fluxin = 0.
                                 do k = max_q, min_q, -1
                                    rhor = rho(i, k, j)
                                    qgr = qg(i, k, j)
                                    fluxout = rhor*vtg(i, k, j)*qgr
                                    flux = (fluxin - fluxout)/rhor/dz8w(i, k, j)
                                    qgr = qgr + del_tv*flux
                                    qg(i, k, j) = dmax1(0., qgr)
                                    fluxin = fluxout
                                    !gsed(i, k, j) = gsed(i, k, j) + fluxin
                                 end do
                                 if (min_q .eq. 1) then
                                    pptgraul(i, j) = pptgraul(i, j) + fluxin*del_tv
                                 else
                                    qg(i, min_q - 1, j) = qg(i, min_q - 1, j) &
                                       + del_tv*fluxin/rho(i, min_q - 1, j) &
                                       /dz8w(i, min_q - 1, j)
                                 end if
               !
                              else
                                 notlast = .false.
                              end if
               !
                           END DO
                        end if !ice2 .eq. 0
      !
      !-- cloud ice  (03/21/02) follow Vaughan T.J. Phillips at GFDL
      !
                     elseif (nhydro .eq. 4) then
                        t_del_tv = 0.
                        del_tv = dtb
                        notlast = .true.
               !
                        DO while (notlast)
               !
                           min_q = kte
                           max_q = kts - 1
               !
                           xlandr = xland(i, j)
                           do k = kts, kte - 1

                              vtir = 0.
                              qir = qi(i, k, j)
                              if (qir .gt. cimin) then
                                 min_q = min0(min_q, k)
                                 max_q = max0(max_q, k)

                                 call vti_mks_gpu(improve, rho(i, k, j), &
                                    tz(i, k, j), qir, qv(i, k, j), &
                                    p(i, k, j), xlandr, vtir, &
                                    sat_predict, new_saturation, cimin, t0, &
                                    rhoe_s, itble, thrd, c218, c580, c76)

                                 ! EMK:  Avoid division by zero
                                 if ((vtir .gt. 1.0e-20)) then
                                    del_tv = dmin1(del_tv, z1(i, k, j)/vtir)
                                 end if
                              end if
                              vti(i, k, j) = vtir
                           end do

                           if (max_q .ge. min_q) then
               !
               !
               !- Check if the summation of the small delta t >=  big delta t
               !             (t_del_tv)          (del_tv)             (dtb)

                              t_del_tv = t_del_tv + del_tv

                              if (t_del_tv .ge. dtb) then
                                 notlast = .false.
                                 del_tv = dtb + del_tv - t_del_tv
                              end if

               ! use small delta t to calculate the qiz flux
               ! termi is the qiz flux pass in the grid box through the upper boundary
               ! termo is the qiz flux pass out the grid box through the lower boundary
               !

                              fluxin = 0.
                              do k = max_q, min_q, -1
                                 rhor = rho(i, k, j)
                                 qir = qi(i, k, j)
                                 fluxout = rhor*vti(i, k, j)*qir
                                 flux = (fluxin - fluxout)/rhor/dz8w(i, k, j)
                                 qir = qir + del_tv*flux
                                 qi(i, k, j) = dmax1(0., qir)
                                 fluxin = fluxout
                                 !ised(i, k, j) = ised(i, k, j) + fluxin
                              end do
                              if (min_q .eq. 1) then
                                 pptice(i, j) = pptice(i, j) + fluxin*del_tv
                              else
                                 qi(i, min_q - 1, j) = qi(i, min_q - 1, j) &
                                    + del_tv*fluxin/rho(i, min_q - 1, j) &
                                    /dz8w(i, min_q - 1, j)
                              end if
               !
                           else
                              notlast = .false.
                           end if
               !
                        END DO !notlast
                     end if
                  end if
               end do
            end do
         end do
      end if !SL_sedi

!   #ifdef EXT_DIAG
!      do j = jts, jte
!         do i = its, myim(j)
!            do k = kts, kte
!               preci3d(i, k, j) = ised(i, k, j)
!               precs3d(i, k, j) = ssed(i, k, j)
!               precg3d(i, k, j) = gsed(i, k, j)
!               precr3d(i, k, j) = rsed(i, k, j)
!            end do
!         end do
!      end do
!   #endif

   !   prnc(i,j)=prnc(i,j)+pptrain
   !   psnowc(i,j)=psnowc(i,j)+pptsnow
   !   pgrauc(i,j)=pgrauc(i,j)+pptgraul
   !   picec(i,j)=picec(i,j)+pptice
   !
      !$acc parallel loop gang vector collapse(2) async(async_id)
      do j = jts, jte
         do i = its, ite
            if (i .le. myim(j)) then

               icenc(i, j) = icenc(i, j) + pptice(i, j)
      !   snowncv(i,j) = pptsnow
               snownc(i, j) = snownc(i, j) + pptsnow(i, j)
      !   graupelncv(i,j) = pptgraul
               graupelnc(i, j) = graupelnc(i, j) + pptgraul(i, j)
      !   RAINNCV(i,j) = pptrain + pptsnow + pptgraul + pptice
               RAINNC(i, j) = RAINNC(i, j) + pptrain(i, j) + pptsnow(i, j) &
                              + pptgraul(i, j) + pptice(i, j)
               pptall(i, j) = pptrain(i, j) + pptsnow(i, j) &
                              + pptgraul(i, j) + pptice(i, j)
               sr(i, j) = 0.
               if (pptall(i, j) .gt. 0.) sr(i, j) = (pptsnow(i, j) + pptgraul(i, j) &
                                              + pptice(i, j))/pptall(i, j)
            end if
         END DO
      END DO
      !$acc exit data delete(fv, pptrain, pptsnow, pptgraul, pptice, &
      !$acc&     orhoz, tz) async(async_id)
      !$acc exit data delete(vtr, vts, vtg, vti, z1, pptall) &
      !$acc&     async(async_id)
!!!!!!!!!!!! end fall_flux_A_loop
      !if (benchmark .eq. .true.) call nvtxEndRange


      !if (SL_sedi) deallocate(qden)


      RETURN
   END SUBROUTINE fall_flux_gpu

!-----------------------------------------------------------------------
!c Correction of negative values
   SUBROUTINE negcor(X, rho, dz8w, &
                     ims, ime, jms, jme, kms, kme, & ! memory dims
                     itimestep, ics, &
                     its, ite, jts, jte, kts, kte) ! tile   dims
!-----------------------------------------------------------------------
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(INOUT) ::                                     X
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(IN) ::                              rho, dz8w
      integer, INTENT(IN) ::                           itimestep, ics

!c Local variables
      REAL   ::   A0, A1, A2

      A1 = 0.
      A2 = 0.
      do k = kts, kte
         do j = jts, jte
            do i = its, ite
               A1 = A1 + max(X(i, k, j), 0.)*rho(i, k, j)*dz8w(i, k, j)
               A2 = A2 + max(-X(i, k, j), 0.)*rho(i, k, j)*dz8w(i, k, j)
            end do
         end do
      end do

      A0 = 0.0

      if (A1 .NE. 0.0 .and. A1 .GT. A2) then
         A0 = (A1 - A2)/A1

         do k = kts, kte
            do j = jts, jte
               do i = its, ite
                  X(i, k, j) = A0*DMAX1(X(i, k, j), 0.0)
               end do
            end do
         end do
      end if

   END SUBROUTINE negcor

   SUBROUTINE consat_s_gpu(ihail, itaobraun, improve)

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

      implicit none

      integer :: ihail, itaobraun, improve

      integer :: k
      real :: cpi, cpi2, grvt, cd1, cd2
      real :: ga3, ga4, ga5, ga6, ga7, ga8, ga9, ga3b, ga4b, ga6b, ga5bh, &
              ga3g, ga4g, ga5gh, ga3d, ga4d, ga5dh, ga6d
      real :: tca, dwv, dva, scv
      real :: ac1, ac2, ac3, bc1, cc1, dc1
      real :: esi, esw, esc, egs, erc, ehs, ehg, egc, eri, erw, esr, eiw, &
              egw, egi, egr, ehw, ehi, ehr
      real :: sc13
      real :: amw, ami, ars, amc
      real :: rw, cw, ci
      real :: ui50, ri50, cmn, y1, apri, bpri

!      if (improve .eq. 3) ihail = 0
!      if (ihail .eq. 1) improve = -20
!      al = 2.5e10
      rhoe_s = 1.29           !surface air density (kg/m^3)
      rd1 = 1.e-3
      rd2 = 2.2

      cpi = 4.*atan(1.)
      cpi2 = cpi*cpi
      grvt = 980.

      ! Fit of Goff-Gratch equation :
      c76 = 7.66
      c141 = 1.414435e7
      c149 = 1.496286e-5
      c172 = 17.26939
      c218 = 21.87456
      c358 = 35.86
      c409 = 4098.026
      c580 = 5807.695
      c610 = 6.1078e3
      c879 = 8.794142
!
      cd1 = 6.e-1
      cd2 = 4.*grvt/(3.*cd1)
      tca = 2.43e3
      dwv = 0.226
      dva = 1.718e-4
      amw = 18.016
      ars = 8.314e7
      scv = 2.2904487
!
      t0 = 273.16
      t00 = 238.16

      ! heat capacity and latent heat (in CGS) :
      rw = 4.615e+6     !gas constant of vapor
      cw = 4.187e+7     !specific heat capacity of liquid water
      ci = 2.093e+7     !specific heat capacity of ice
      cp = 1.004e7      !specific heat capacity of dry air
      alv = 2.5e+10     !latent heat of evaporation/condensation
      alf = 3.336e+9    !latent heat of melting/freezing
      als = 2.8336e+10  !latent heat of sublimation/deposition
      avc = alv/cp
      afc = alf/cp
      asc = als/cp

!***   DEFINE THE DENSITY AND SIZE DISTRIBUTION OF PRECIPITATION
      ! tng  : intercept parameter of graupel/hail (1/cm**4)
      ! roqg : density of graupel/hail (g/cm**3)
      ! ag, bg : coefficients used in terminal velocity
      if (ihail .eq. 1) then
         roqg = 0.9
         ag = sqrt(cd2*roqg)
         bg = 0.5
         tng = 0.002
      else
         roqg = 0.3
         ag = 330.22    !for bulk density of 0.3
         bg = 0.36
         tng = 0.04
      end if

      ! tns  : intercept parameter of snow (1/cm**4)
      ! roqs : density of snow (g/cm**3)
      ! as, bs : coefficients used in terminal velocity
      roqs = 0.05
      tns = 0.1
      as = 151.01
      bs = 0.24

      ! tnw  : intercept parameter of rain (1/cm**4)
      ! roqr : density of rain (g/cm**3)
      ! aw, bw : coefficients used in terminal velocity
      aw = 2115.
      bw = 0.8
      roqr = 1.  !not defined in Steve's
      tnw = 0.08  !not defined in Steve's

      ! Lin et al. 1983 :
      constb = 0.8
      constd = 0.11
      consta = 2115.0*0.01**(1 - constb)   ! =841.9967
      constc = 78.63*0.01**(1 - constd)      ! =11.72
      o6 = 1./6.
      cdrag = 0.6
      abar = 19.3
      bbar = 0.37

      bgh = 0.5*bg
      bsh = 0.5*bs
      bwh = 0.5*bw
      bgq = 0.25*bg
      bsq = 0.25*bs
      bwq = 0.25*bw

      ! unit in MKS :
      xnor = tnw*1.0e8        !intercept parameter of rain (m^-4)
      xnos = tns*1.0e8        !intercept parameter of snow (m^-4)
      xnog = tng*1.0e8        !intercept parameter of graupel (m^-4)
      rhowater = roqr*1000.   !density of rain (kg/m^3)
      rhosnow = roqs*1000.    !density of snow (kg/m^3)
      rhograul = roqg*1000.   !density of graupel (kg/m^3)
      if (ihail .eq. 1) then
         xnoh = xnog          !intercept parameter of hail (m^-4)
         rhohail = rhograul   !density of hail (kg/m^3)
      end if

!***   GAMMA FUNCTION CALCULATIONS   ********
      ga3 = 2.
      ga4 = 6.
      ga5 = 24.
      ga6 = 120.
      ga7 = 720.
      ga8 = 5040.
      ga9 = 40320.

      ga3b = gammagce_gpu(3.+bw)        !not used
      ga4b = gammagce_gpu(4.+bw)        !not used
      ga6b = gammagce_gpu(6.+bw)        !not used
      ga5bh = gammagce_gpu((5.+bw)/2.)  !not used
      ga3g = gammagce_gpu(3.+bg)
      ga4g = gammagce_gpu(4.+bg)
      ga5gh = gammagce_gpu((5.+bg)/2.)
      ga3d = gammagce_gpu(3.+bs)
      ga4d = gammagce_gpu(4.+bs)
      ga5dh = gammagce_gpu((5.+bs)/2.)
      ga6d = gammagce_gpu(6.+bs)

!      if (improve .eq. 3) then
!         ga4g = 11.63177     !bg=0.5
!         ga3g = 3.3233625    !bg=0.5
!         ga5gh = 1.608355    !bg=0.5
!         if (bg .eq. 0.37) ga4g = 9.730877
!         if (bg .eq. 0.37) ga3g = 2.8875
!         if (bg .eq. 0.37) ga5gh = 1.526425
!         if (bg .eq. 0.36) ga4g = 9.599978
!         if (bg .eq. 0.36) ga3g = 2.857136
!         if (bg .eq. 0.36) ga5gh = 1.520402
!         ga3d = 2.54925      !bs=0.25
!         ga4d = 8.285063     !bs=0.25
!         ga5dh = 1.456943    !bs=0.25
!         if (bs .eq. 0.57) ga3d = 3.59304
!         if (bs .eq. 0.57) ga4d = 12.82715
!         if (bs .eq. 0.57) ga5dh = 1.655588
!         if (bs .eq. 0.24) ga3d = 2.523508
!         if (bs .eq. 0.24) ga4d = 8.176166
!         if (bs .eq. 0.24) ga5dh = 1.451396
!         if (bs .eq. 0.11) ga3d = 2.218906
!         if (bs .eq. 0.11) ga4d = 6.900796
!         if (bs .eq. 0.11) ga5dh = 1.382792
!      end if
!      ga6d = 144.93124       !bs=0.11
!      if (bs .eq. 0.24) ga6d = 181.654791

!CCCCC        LIN ET AL., 1983 OR LORD ET AL., 1984   CCCCCCCCCCCCCCCCC
      ac1 = aw
      ac2 = ag            ! Steve only defines ac1 and ac2
      ac3 = as            ! need to talk about these 3 parameters.

      bc1 = bw
      cc1 = as
      dc1 = bs

      ! slope parameter :
      zrc = (cpi*roqr*tnw)**0.25
      zsc = (cpi*roqs*tns)**0.25
      zgc = (cpi*roqg*tng)**0.25

      ! terminal velocity of rain (Steve's) :
      vrc = aw*ga4b/(6.*zrc**bw)
      vrc0 = -26.7
      vrc1 = 20600./zrc
      vrc2 = -204500./(zrc*zrc)
      vrc3 = 906000./(zrc*zrc*zrc)

      ! terminal velocity of snow and graupel (Lin et al. 1983) :
      vsc = as*ga4d/(6.*zsc**bs)
      vgc = ag*ga4g/(6.*zgc**bg)

!     ****************************
      rn1 = 9.4e-15
      bnd1 = 6.e-4   !for psaut
      rn2 = 1.e-3
      bnd2 = 2.0e-3
      rn3 = 0.25*cpi*tns*as*ga3d

!      esi=.1
      esi = 0.25 !improve=3
      rn3 = 0.25*cpi*tns*as*esi*ga3d

      esw = 1.   ! Steve uses esc and ac1
!      esc=1.
      esc = 0.45  !improve=3
!      rn4=.25*cpi*esw*tns*as*ga3d
      rn4 = 0.25*cpi*esc*tns*as*ga3d  !improve=3
      eri = 0.1
!      rn5=.25*cpi*eri*tnw*ac1*ga3b
      rn5 = 0.25*cpi*eri*tnw  !improve=3
      rn50 = -.267e2*ga3
      rn51 = 5.15e3*ga4
      rn52 = -1.0225e4*ga5
      rn53 = 7.55e3*ga6

      ami = 1./(24.*6.e-9)
!      rn6=cpi2*eri*tnw*ac1*roqr*ga6b*ami
      rn6 = cpi2*eri*tnw*roqr*ami  !improve=3
      rn60 = -.267e2*ga6
      rn61 = 5.15e3*ga7
      rn62 = -1.0225e4*ga8
      rn63 = 7.55e3*ga9

!      esr=.5
      esr = 1.  !improve=3
      rn7 = cpi2*esr*tnw*tns*roqs
      rn8 = cpi2*esr*tnw*tns*roqr

      egs = 0.1

!      rn9=cpi2*tns*tng*roqs
      rn9 = cpi2*egs*tns*tng*roqs  !improve=3

!      rn10=2.*cpi*tns
      rn10 = 4.*tns  !improve=3
!      rn101 = 0.31*ga5dh*sqrt(cc1)
!      rn10a = als*als/rw
      rn101 = 0.65  !improve=3
      rn102 = 0.44*sqrt(as/dva)*ga5dh  !improve=3
      rn10a = alv*als*amw/(tca*ars)  !improve=3
      rn10b = alv/tca
      rn10c = ars/(dwv*amw)

!      rn11=2.*cpi*tns/alf
      rn11 = 2.*cpi*tns*tca/alf  !improve=3
      rn11a = cw/alf

!      ami50=3.84e-6
!      ami40=3.08e-8
!      ami50=4.8e-7*(100./50.)**3         ! Roger: actually = 3.84e-6
!      ami40=2.41e-8
!      ami50=3.76e-8
      ami50 = 4.8e-7  !improve=3
      ami40 = 2.46e-7  !improve=3
      ami100 = 1.51e-7

      eiw = 1.
      ui50 = 100. ! 6/15/02 tao's
!      ri50=2.*5.e-3
      ri50 = 5.e-3  !improve=3

      cmn = 1.05e-15
      rn12 = cpi*eiw*ui50*ri50**2
      do k = 1, 31
         y1 = 1.-aa2(k)
         rn13(k) = aa1(k)*y1/(ami50**y1 - ami40**y1)
         rn12a(k) = rn13(k)/ami50
         rn12b(k) = aa1(k)*ami50**aa2(k)
         rn25a(k) = aa1(k)*cmn**aa2(k)
         BergCon1(k) = 6.*aa1(k)*ami50**(aa2(k) - 1.)
         BergCon2(k) = -2.*aa1(k)*ami50**aa2(k)*1.2
         BergCon3(k) = 6.*aa2(k)/((aa2(k) + 1.)*(aa2(k) + 2.)) &
                       *aa1(k)*ami50**(aa2(k) - 1.)
         BergCon4(k) = 2.*(1.-aa2(k))/((aa2(k) + 1.)*(aa2(k) + 2.)) &
                       *aa1(k)*ami50**aa2(k)*1.2
      end do

!      egw=1.
      egw = 0.65  !improve=3
      rn14 = 0.25*cpi*egw*tng*ga3g*ag

      egi = 0.1
!      egi=.001
      rn15 = 0.25*cpi*egi*tng*ga3g*ag
      egi = 1.
      rn15a = 0.25*cpi*egi*tng*ga3g*ag

      egr = 1.
      rn16 = cpi2*egr*tng*tnw*roqr
      rn17 = 2.*cpi*tng
!      rn17a = 0.31*ga5gh*sqrt(ag)
      rn17b = cw - ci
      rn17c = cw
      rn171 = 2.*cpi*tng*alv*dwv  !improve=3
      rn172 = 2.*cpi*tng*tca  !improve=3
      rn17a = 0.31*ga5gh*sqrt(ag/dva)  !improve=3

      apri = 0.66
      bpri = 1.e-4
      bpri = 0.5*bpri                        ! 6/17/02 tao's
      rn18 = 20.*cpi2*bpri*tnw*roqr
      rn18a = apri
!      rn19=2.*cpi*tng/alf
      rn191 = 0.78                     !Are this same as rnn191 (listed below)?
      rn192 = 0.31*ga5gh*sqrt(ag/dva)  !Are this same as rnn192 (listed below)?
!      rn19a=.31*ga5gh*sqrt(ag)
      rn19b = cw/alf

      rn19 = 2.*cpi*tng*tca/alf  !improve=3
      rn19a = cw/alf  !improve=3

      rn20 = 2.*cpi*tng
!      rn20a=als*als/rw
!      rn20b=.31*ga5gh*sqrt(ag)
      rn30a = alv*alv*amw/(tca*ars) ! EMK per Roger's 20110816 code
      rn30 = 2.*cpi*tng             ! EMK per Roger's 20110816 code
      rn20a = als*als*amw/(tca*ars)  !improve=3
      rn20b = als/tca  !improve=3

      bnd3 = 2.e-3
!      rn21=1.e3*1.569e-12/0.15
      rn21 = 1.e-3  !improve=3
      bnd21 = 1.5e-3

!      erw=1. ! erc in Steve's code
!      rn22=.25*cpi*erw*ac1*tnw*ga3b
      erc = 1.  !improve=3
      rn22 = 0.25*cpi*erc*tnw  !improve=3

      rn23 = 2.*cpi*tnw
      rn23a = 0.31*ga5bh*sqrt(ac1)
      rn23b = alv*alv/rw
      rn231 = 0.78  !improve=3
      rn232 = 0.31*ga3*sqrt(3.e3/dva)  !improve=3

!      if (improve .eq. 3) itaobraun=1
      if (itaobraun .eq. 0) then
         cn0 = 1.e-8
         beta = -0.6
      else if (itaobraun .eq. 1) then
         cn0 = 1.e-6
         beta = -0.46
      end if

      rn25 = cn0
      rn30b = alv/tca
      rn30c = ars/(dwv*amw)
      rn31 = 1.e-17

      rn32 = 4.*51.545e-4
      rn33 = 4.*tns
      rn331 = 0.65
      rn332 = 0.44*sqrt(as/dva)*ga5dh

      amc = 1./(24.*4.e-9)  !improve=3
      rn34 = cpi2*esc*amc*as*roqs*tns*ga6d  !improve=3
      rn35 = alv*alv/(cp*rw)  !improve=3

      draimax = 0.0500       !maximum rain diameter (cm)
      draimax = draimax**4.*roqr*cpi

      ! critical q of hydrometeor characteristics (eg. fall speed, radius)
      cwmin = 1.e-12
      cimin = 1.e-12
      crmin = 1.e-10
      csmin = 1.e-10
      cgmin = 1.e-10

      return
   END SUBROUTINE consat_s_gpu

!JJS
!JJS      REAL FUNCTION GAMMA(X)
!JJS        Y=GAMMLN(X)
!JJS        GAMMA=EXP(Y)
!JJS      RETURN
!JJS      END
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!JJS      real function GAMMLN (xx)
   real function gammagce_gpu(xx)
   !$acc routine seq
!**********************************************************************
      real*8 cof(6), stp, half, one, fpf, x, tmp, ser
      data cof, stp/76.18009173, -86.50532033, 24.01409822, &
         -1.231739516, .120858003e-2, -.536382e-5, 2.50662827465/
      data half, one, fpf/.5, 1., 5.5/
!
      x = xx - one
      tmp = x + fpf
      tmp = (x + half)*log(tmp) - tmp
      ser = one
      do j = 1, 6
         x = x + one
         ser = ser + cof(j)/x
      end do !j
      gammln = tmp + log(stp*ser)
!JJS
      gammagce_gpu = exp(gammln)
!JJS
      return
   END FUNCTION gammagce_gpu

! compute base snow/graupel intercept scaling factor
   subroutine sgmap_gpu(isg, qsg, r00, tairc, ftnsg, roqs, roqg, tns, tng)
      !$acc routine seq
      implicit none

!      common/size/ tnw,tns,tng,roqs,roqg,roqr  !defined in the beginning of the module
      integer, intent(in)  :: isg
      real, intent(in)  :: qsg, r00, tairc
      real, intent(out) :: ftnsg
      
      real, intent(in) :: roqs, roqg, tns, tng ! acc routine

! LOCAL variables

!  for snomap
      real ::  xs, sno11, sno00, dsno11, dsno00, sexp11, sexp00, stt, &
              stexp, sbase, tslopes, dsnomin, slim
!  for grpmap
      real ::  xg, grp11, grp00, dgrp11, dgrp00, gexp11, gexp00, gtt, &
              gtexp, gbase, tslopeg, dgrpmin, glim

      real :: taird, qsg1, xx, fexp, cpi, cmin
      real :: ftnsT, sno1, dsno1, sexp1, ftnsQ, tnsmax
      real :: ftngT, grp1, dgrp1, gexp1, ftngQ, tngmax

      CPI = 4.*ATAN(1.)
      cmin = 1.e-20

      tslopes = 0.03579323   ! increase tns by 3.5 from 0 to -35C
      tslopeg = 0.03138892   ! increase tng by  3 from  0 to -35C

      dsnomin = 0.0110       !minimum snow diameter (cm)
      dgrpmin = 0.0145       !minimum graupel diameter (cm)

      dsnomin = dsnomin**4.*roqs*cpi
      dgrpmin = dgrpmin**4.*roqg*cpi

      xs = 0.97
      sno11 = 0.65                !cold aloft
      sno00 = 0.30                !near melting level
      dsno11 = 1.25               !Tao 02/23/2012 used to be 1.15
      dsno00 = 0.60               !Tao 02/23/2012 used to be 0.70
      sexp11 = 1.2
      sexp00 = 0.7
      stt = -20.
      stexp = 0.25
      slim = 0.95
      sbase = 0.00110

      xg = 0.98
      grp11 = 0.75
      grp00 = 0.45
      dgrp11 = 3.95
      dgrp00 = 0.25
      gexp11 = 0.35
      gexp00 = 0.6
      gtt = -20.
      gtexp = 0.30
      glim = 0.90
      gbase = 0.0058

      ftnsg = 1.

      if (qsg .gt. cmin) then

         taird = min(0., max(-35., tairc) + 0.0)
         qsg1 = qsg*r00*1.e6

         if (isg .eq. 1) then                          !snow
            ftnsT = exp(-1.*tslopes*taird)
            sno1 = sno11
            dsno1 = dsno11
            sexp1 = sexp11

            if (taird .gt. stt) then
               sno1 = sno00 - (sno00 - sno11)*(taird/stt)**stexp
               dsno1 = dsno00 - (dsno00 - dsno11)*(taird/stt)**stexp
               sexp1 = sexp00 - (sexp00 - sexp11)*(taird/stt)**stexp
            end if !taird

            xx = xs - xs*min(slim, max(0., (qsg1 - sno1)/dsno1)**sexp1)
            ftnsT = ftnsT**xx
            fexp = xx
            ftnsQ = 1.0
            ftnsQ = (qsg1/sbase)**fexp
            ftnsg = ftnsT*ftnsQ
            tnsmax = r00*qsg/dsnomin

            if (ftnsg*tns .gt. tnsmax) ftnsg = tnsmax/tns

         else if (isg .eq. 2) then                      !graupel

            ftngT = exp(-1.*tslopeg*taird)
            grp1 = grp11
            dgrp1 = dgrp11
            gexp1 = gexp11

            if (taird .gt. gtt) then
               grp1 = grp00 - (grp00 - grp11)*(taird/gtt)**gtexp
               dgrp1 = dgrp00 - (dgrp00 - dgrp11)*(taird/gtt)**gtexp
               gexp1 = gexp00 - (gexp00 - gexp11)*(taird/gtt)**gtexp
            end if !taird

            xx = xg - xg*min(glim, max(0.0, (qsg1 - grp1)/dgrp1)**gexp1)
            ftngT = ftngT**xx
            fexp = xx
            ftngQ = 1.0
            ftngQ = (qsg1/gbase)**fexp
            ftnsg = ftngT*ftngQ
            tngmax = r00*qsg/dgrpmin
            if (ftnsg*tng .gt. tngmax) ftnsg = tngmax/tng

         end if !isg

      end if !qsg

      return
   end subroutine sgmap_gpu

   SUBROUTINE saticel_s_gpu(myim, dt, ihail, itaobraun, ice2, istatmin, &
                        new_ice_sat, id, improve, xlat, sdec, &
                        ptwrf, qvwrf, qlwrf, qrwrf, &
                        qiwrf, qswrf, qgwrf, &
                        rho_mks, pi_mks, p0_mks, w_mks, &
                        itimestep, xland, &
                        refc, refr, refi, refs, refg, & ! cloud effective radius
                        ims, ime, jms, jme, kms, kme, &
                        its, ite, jts, jte, kts, kte, &
#ifdef Readaeroclx
                        aeroclx, naero, &
#endif
#ifdef EXT_DIAG
                        refl_10cm, diagflag, do_radar_ref, & ! GT added for reflectivity calcs
                        physc, physe, physd, physs, physm, physf, &
                        acphysc, acphyse, acphysd, acphyss, acphysm, acphysf, &
#endif
                        sat_predict, new_saturation, use_cpm, use_declination, &
                        benchmark, async_id) !flags
      !$acc routine(esw_mks_gpu) seq
      !$acc routine(esi_mks_gpu) seq
      !$acc routine(vtr_mks_gpu) seq
      !$acc routine(sgmap_gpu) seq
      !$acc routine(vts_mks_gpu) seq
      !$acc routine(vti_mks_gpu) seq
      !$acc routine(vtg_mks_gpu) seq
      !$acc routine(mass2ccn_gpu) seq
      !$acc routine(mass2icn_gpu) seq
      !$acc routine(eff_rad_gpu) seq
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
      integer, parameter ::  nt = 2880, nt2 = 2*nt

!cc   using scott braun's way for pint, pidep computations
      integer  ::   itaobraun, ice2, ihail, new_ice_sat, id, istatmin
      integer  ::   itimestep, improve
      real     ::   tairccri, dt
!cc

!  integer ids,ide,jds,jde,kds,kde
      integer, intent(in) :: myim(my_max)
      integer ims, ime, jms, jme, kms, kme
      integer its, ite, jts, jte, kts, kte
      integer i, j, k, kp, n

      real, dimension(ims:ime, kms:kme, jms:jme) :: afcp, ascp, avcp, &
         pi0, pir, r00, rp0
      real   ::   a0, a1, a2,  alvr, &
                  betah, bg3, bgh5, bs3, bs6, bsh5, &
                  bw3, bw6, bwh5, cmin, cmin1, cmin2,  &
                  d2t, del, ee1, ee2, &
                  f00, f2, f3, ft, fvs, qb0, &
                  r0s, r101f, r10ar, r10t, r11at, r11rt, &
                  r14r, r15af, r15ar, r15r, &
                  r16r, r17aq, r17as, r17r, r19aq, r19as, &
                  r19bt, r19rt, r20bq, r20bs, r20t, &
                  r23t, r25a, r25rt, r2ice, r31r, &
                  r9r, r_nci, rft, &
                  rijl2, rrq, rrs, rt0, scc, sccc, &
                  sddd, see, seee, sfff, smmm, ssss, tb0, temp, &
                  ucog, ucor, ucos, uwet, vrcf, &
                  rdt

      real :: a_1, a_2, a_3
      real :: a_11, a_22, a_33, a_44
      real :: zdry, zwet

      real :: fact_fit

      real, dimension(its:ite, kts:kte, jts:jte) ::  fv
      real, dimension(its:ite, kts:kte, jts:jte) ::  qiwrf_old

      real, dimension(ims:ime, kms:kme, jms:jme) ::  ptwrf, qvwrf
      real, dimension(ims:ime, kms:kme, jms:jme) ::  qlwrf, qrwrf, &
                                                     qiwrf, qswrf, qgwrf

      real, dimension(ims:ime, kms:kme, jms:jme) ::  rho_mks
      real, dimension(ims:ime, kms:kme, jms:jme) ::  pi_mks
      real, dimension(ims:ime, kms:kme, jms:jme) ::  p0_mks
      real, dimension(ims:ime, kms:kme, jms:jme) ::  w_mks

      real, dimension(ims:ime, kms:kme, jms:jme), INTENT(OUT) :: refc, refr, &
         refi, refs, refg

      real :: cv409, cs580

      real :: y1, y2, y3, y4, y5, qsw, ssw, rtair, dm, rsub1, esi, esw, qsi, ssi, &
         prn, psn, xlv, cpm1, xls, xlf, qvs, dd1, dd, ern, pmlts, pmltg, &
         psdep, pgdep, pssub, pgsub, qracs, cnd, dep, pint, psmlt, pgmlt, &
         pimlt, pgfr, psacr, pihom, pidw, f0, zrr, zsr, zgr, vscf, vgcf, &
         dwvp, pidep, wgacr, pracs, qsacr, pgaut, dgaci, wgaci, pgwet, &
         dlt2, pact, fez, abw, abi, dv1, tauc, taui, r3f, r4f, r5f, r6f, &
         cp409, cp580, r7r, r8r, vgcr, r12r, r14f, r15f, r18r, r22f, &
         r23af, r23br, r32rt, r9rf, r16rf, r101r, r102rf, r191r, r192rf, &
         r331r, r332rf, r34f, r231r, r232rf, dwv, tca, scv, rr0, pr0, pr, ps, pg, &
         pwacs, psaut, ftns, ftng,  pihms, piacr, pihmg, praut, pracw, &
         psfw, psfi, dgacs, dgacw, dgacr, pgacs, wgacs, qgacw, &
         psaci, psacw, qsacw, praci, qgacr, pimm, pcfr, tairc 
      
      ! GPU register for argument variables
      real :: xlandr, ptwrfr, qvwrfr, qlwrfr, qrwrfr, qiwrfr, qswrfr, qgwrfr, &
         ftns0r, ftng0r, dlt1, dlt3, dlt4, rhoair, taur, tairr, ascpr, afcpr, &
         avcpr, pi0r, pirr, pr0r, rp0r, r00r, vrr, vsr, vgr, vir, p0r1, refcr, &
         refrr, refir, refsr, refgr

      !real, dimension(its:ite, jts:jte) :: asss


      !real, dimension(its:ite, jts:jte) :: egs, ddb

      real, dimension(its:ite, kts:kte, jts:jte) :: tair, &
         zr, vr, zs, vs, vg, zg, vi

      !real, dimension(its:ite, jts:jte) :: rq, col, dda

      real, dimension(its:ite, kts:kte, jts:jte) ::  rho
      !real, dimension(kts:kte) :: &
      !   tb, qb, rho1, &
      !   ta, qa, ta1, qa1, &
      !   coef, z1, z2, z3, &
      !   am, am1, ub, vb, &
      !   wb, ub1, vb1, rrho, &
      !   rrho1, wbx

      !real, dimension(kts:kte) :: &
      !   fd, fe, &
      !   st, sv, &
      !   sq, sc, &
      !   se, sqa

      !real, dimension(kts:kte) :: &
      !   srro, qrro, sqc, sqr, &
      !   sqi, sqs, sqg, stqc, &
      !   stqr, stqi, stqs, stqg
      !real, dimension(nt) :: &
      !   tqc, tqr, tqi, tqs, tqg

      !real, dimension(ims:ime, jms:jme) :: &
      !   y0, ts0, qss0

      !integer, dimension(its:ite, jts:jte, 4) ::    ics

      integer :: i24h, it
      integer :: iwarm
      real :: r2is, r2ig

#ifdef EXT_DIAG
      real, dimension(ims:ime, kms:kme, jms:jme)  :: &
         physc, physe, physd, &
         physs, physm, physf, &
         acphysc, acphyse, acphysd, &
         acphyss, acphysm, acphysf
      real, dimension(its:ite, kts:kte, jts:jte) :: dbz
#endif


      integer  ::  ihalmos
      real     ::  tslopes, tslopeg
      real     ::  xnsplnt, xmsplnt
      real     ::  hmtemp1, hmtemp2, hmtemp3, hmtemp4
      real     ::  ftnsQ, ftngQ, fexp
      real     ::  xssi, fssi
      real     ::  efsi
      real     ::  dmicrons, dmicrong, fdms, fdmg, dvair, alpha
      !real, dimension(its:ite, jts:jte) :: tairN, tairI

! for Xiping's new dbz code
      real     :: a_c, a_i
      real     :: re_c, re_i, re_s
      real     :: w_c, w_i, w_s
      real     :: ze_cld

      real, dimension(its:ite, kts:kte, jts:jte) :: ftns0, ftng0
      integer  :: improve1
      data improve1/-20/

      real     ::  tairc5, hfact, sfact, yy1
      real     ::  xncld, esat, rv, rlapse_m
      real     ::  delT, bhi, cpi
      real     ::  rc, ra, cna
      real     ::  xccld, xknud, cunnf, diffar

      real     :: r7rf, r8rf
      real     :: r11t, r19t, r19at
      real     :: r30t, r33t
      real     :: rn1s
      real     :: ami20

      REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN)   :: XLAND
      real, parameter :: roqi = 0.9179    ! ice density
!  real, parameter :: ccn_over_land = 1500  ! [#/cm3] climatological value
!  real, parameter :: ccn_over_water = 150  ! [#/cm3] climatological value
      real, parameter :: ccn_over_land = 159  ! [#/cm3] climatological value (from MERRA2)
      real, parameter :: ccn_over_water = 66  ! [#/cm3] climatological value (from MERRA2)
      real :: L_cloud    ! cloud water [g/cm3] !
      real :: I_cloud    ! cloud water [g/cm3] !
      real :: mu, ccn_ref, lambda
      real :: gamfac1, gamfac3

!+---+-----------------------------------------------------------------+
#ifdef EXT_DIAG
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT):: refl_10cm  ! GT
      LOGICAL, OPTIONAL, INTENT(IN) :: diagflag
      INTEGER, OPTIONAL, INTENT(IN) :: do_radar_ref
#endif
!+---+-----------------------------------------------------------------+
      integer, parameter :: rewflag = 2
      ! 1 : default
      ! 2 : mapping spectrum, different over land and ocean
      real :: mdc1, mdc2, mdc3, mdc4, mdc5, mdc6, mvdc
      real :: efd1, efd2, efd3, efd4, efd5, efd6, efdc
      real :: ltk, ltk2
      real :: lqc, lqc2

      integer, parameter :: reiflag = 5
      ! 1 : default
      ! 2 : Wyser 1998
      ! 3 : Heymsfild et al. 2014
      ! 4 : Mitchell et al. 2011
      ! 5 : mapping spectrum, different over land and ocean
      real :: reimin, reimax
      real :: iwc_0, bw98
      real :: md22, bd22, sigma, cd22, xd22
      real :: lqi, lqi2, efdi

      integer, parameter :: ccnflag = 1
      ! 1 : default
      ! 2 : WRF GOCART mass2ccn (LUT approach)

      integer, parameter :: inflag = 1
      ! 1 : default (Meyers et al. 1992)
      ! 2 : WRF GOCART mass2icn (DeMotto et al. 2010)
      ! 3 : Hong et al. 2004
      ! 4 : Cooper curve (Cooper 1986 ; Chern et al. 2016)
      ! 5 : Fletcher 1962

      ! calculate allowable ice supersautration :
      real, intent(in) :: xlat(my_max)
      real :: d2r, arg

      ! calculate solar declination angle :
      real, intent(in) :: sdec

#ifdef Readaeroclx
      ! aerosol climatology :
      integer, intent(in) :: naero
      real, dimension(ims:ime, kms:kme, naero, jms:jme), intent(in) :: aeroclx

      ! WRF GOCART coupling :
!      integer, parameter :: ngo = 14
      real, dimension(nlut, ims:ime, kms:kme, jms:jme) :: aerog
      real :: rhw, ssrw, p_mb, mr2mc

      real :: P_liu_daum  ! autoconversion rate [g/cm3 s-1]
      real :: re_liu_daum ! effective radius of cloud [micron]
      real,dimension(nlut, its:ite, kts:kte, jts:jte)  :: lut_ccn_interp
#endif

!flags
      logical, intent(in) :: sat_predict, new_saturation, use_cpm, use_declination


!sat_predict
      ! saturation prediction scheme :
      integer :: hid
      real :: tau, atem
      real :: C1, K1, ncloud, nact, qcmax, mvrc
      real :: qimax, nice, inhgr, rhoi, mvdi, mvri
      real :: lqr, lqr2, mvdr, efdr, kmin, kmax, kdxr, afar, &
              tnr, lzr, bvr, avr, mur, rhoaj, gr2, gbr25
      real :: cnd1, cnd2, dep1, dep2, fez1, fez2, latr, ern1, ern2

      ! transition zone :
      real :: dltd, nbd, sbd, nbd1, nbd2, sbd1, sbd2, latint
      real :: fxlat
      

      !cp=1.004e7
      !alv=2.5e10 ; alf=3.336e9 ; als=2.8336e10
      !rw=4.615e6 ; cw=4.187e7 ; ci=2.093e7
      !avcp=alv/cp*pir
      real, parameter :: cvap = 1.846e+7    !specific heat capacity of vapor (in CGS)
      real, parameter :: cliq = 4.218e+7    !specific heat capacity of liquid (in CGS)
      real, parameter :: cice = 2.106e+7    !specific heat capacity of ice (in CGS)
      real :: cpm, hlv, hls, hlf
      
      !for GPU porting
      integer, intent(in) :: async_id
      integer, parameter :: i_check = 50600, k_check = 1100, j_check = 1600
      logical, intent(in) :: benchmark
      real :: pir1, p0r, rhor, pirr, rr0r, fv0r

      
      !if (benchmark .eq. .true.) call nvtxStartRange("GPU:saticel_s")


#ifdef EXT_DIAG
      if (itimestep .eq. 1) then
         !$acc parallel loop gang collapse(2) async(async_id)
         do k = kts, kte
            do j = jts, jte
               !$acc loop vector
               do i = its, myim(j)
!             physc(i,k,j)=0.
!             physe(i,k,j)=0.
!             physd(i,k,j)=0.
!             physs(i,k,j)=0.
!             physf(i,k,j)=0.
!             physm(i,k,j)=0.
                  acphysc(i, k, j) = 0.
                  acphyse(i, k, j) = 0.
                  acphysd(i, k, j) = 0.
                  acphyss(i, k, j) = 0.
                  acphysf(i, k, j) = 0.
                  acphysm(i, k, j) = 0.
               end do !i
            end do !j
         end do !k
         !$acc wait(async_id)
      end if
#endif


#ifdef Readaeroclx
      !$acc data create( tair, zr, zs, zg, ascp, pi0, afcp, avcp, pir, rp0, &
      !$acc&     ftns0, ftng0, r00, rho, fv, qiwrf_old, vr, vs, vg, vi, aerog, &
      !$acc&     lut_ccn_interp) async(async_id)
#else
      !$acc data create( tair, zr, zs, zg, ascp, pi0, afcp, avcp, pir, rp0, &
      !$acc&     ftns0, ftng0, r00, rho, fv, qiwrf_old, vr, vs, vg, vi) &
      !$acc&     async(async_id)

#endif

!
!     ******   THREE CLASSES OF ICE-PHASE   (LIN ET AL, 1983)  *********

      d2t = dt

!   ICE2=0 ! default, 3ice with cloud ice, snow and graupel
!              r2is=1., r2ig=1.
!   ICE2=1 ! 2ice with cloud ice and snow (no graupel) - r2iceg=1, r2ice=0.
!              r2is=1., r2ig=0.
!   ICE2=2 ! 2ice with cloud ice and graupel (no snow) - r2ice=1, r2iceg=0.
!              r2is=0., r2ig=1.
!   ICE2=3 ! no ice, warm rain only

      r2ig = 1.
      r2is = 1.
      iwarm = 0
      if (ice2 .eq. 1) then
         r2ig = 0.
         r2is = 1.
      elseif (ice2 .eq. 2) then
         r2ig = 1.
         r2is = 0.
      elseif (ice2 .eq. 3) then
         r2ig = 0.
         r2is = 0.
         iwarm = 1
      end if

!      If (improve .eq. 3) then
!         ihail = 0
!      endif

      cmin = 1.e-19
      If (improve .eq. 3) cmin = 1.e-20 ! EMK NUWRF Bug fix
      cmin1 = 1.e-20
      cmin2 = 1.e-40
      ucor = 3071.29/tnw**0.75
      ucos = 687.97*roqs**0.25/tns**0.75
      ucog = 687.97*roqg**0.25/tng**0.75
      uwet = 4.464**0.95

      CPI = 4.*ATAN(1.)
      
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
      tslopes = 0.10539656   ! increase tns by 40 from -5 to -40C
      tslopeg = 0.08559235   ! increase tng by 20 from -5 to -40C

!  HALLET-MOSSOP RIME SPLINTERING parameters
      ihalmos = 1
!      if (improve.gt.2) ihalmos=1
      xnsplnt = 370.     ! peak # splinters per milligram of rime
      xmsplnt = 4.4e-8   ! mass of a splinter (from Ferrier 1994)
      hmtemp1 = -2.
      hmtemp2 = -4.
      hmtemp3 = -6.
      hmtemp4 = -8.


      f2 = rd1*d2t
      f3 = rd2*d2t

      ft = dt/d2t
      rt0 = 1./(t0 - t00)

      bw3 = bw + 3.
      bs3 = bs + 3.
      bg3 = bg + 3.
      bsh5 = 2.5 + bsh
      bgh5 = 2.5 + bgh
      bwh5 = 2.5 + bwh
      bw6 = bw + 6.
      bs6 = bs + 6.
      betah = .5*beta

      r10t = rn10*d2t
      r11at = rn11a*d2t
      r19bt = rn19b*d2t
      r20t = -rn20*d2t
      r23t = -rn23*d2t
      r25a = rn25

      if (improve .eq. 3) then
!         itaobraun=1
!         new_ice_sat=3
         rdt = 1./d2t
         r11t = rn11*d2t
         r19t = rn19*d2t
         r19at = rn19a*d2t
         r20t = rn20*d2t
         r23t = rn23*d2t
         r30t = rn30*d2t
         r33t = rn33*d2t

         Rc = 1.e-3               ! cloud droplet radius 10 microns
         Ra = 1.e-5               ! aerosol radius 0.1 microns
         Cna = 500.               ! contact nuclei conc per cc
         Bhi = 1.01e-2            ! pollen (Deihl et al. 2006)
      end if

!C    ******************************************************************

!!!!!!!!!!!! saticel_s_B_loop
      tb0 = 0.
      qb0 = 0.
      cv409 = c409*avc
      cs580 = c580*asc
      
!JJS
#ifdef Readaeroclx
      !$acc parallel loop gang collapse(2) async(async_id) private(r0s, rrs, &
      !$acc&         rrq, fvs, pir1, p0r, rhor, pirr, rr0r, mr2mc)
#else
      !$acc parallel loop gang collapse(2) async(async_id) private(r0s, rrs, &
      !$acc&         rrq, fvs, pir1, p0r, rhor, pirr, rr0r)
#endif
      do j = jts, jte
         do k = kts, kte
            !$acc loop vector
            do i = its, myim(j)
!JJS  convert from mks to cgs, and move from WRF grid to GCE grid
               rhor = rho_mks(i, k, j)*0.001
               w_mks(i, k, j) = w_mks(i, k, j)*100.
               qiwrf_old(i, k, j) = qiwrf(i, k, j)
               
               pir1 = pi_mks(i, k, j)
               p0r = p0_mks(i, k, j)*10.0
               rho(i, k, j) = rhor
               rp0(i, k, j) = 3.799052e3/p0r
               pi0(i, k, j) = pir1
               pirr = 1./pir1
               pir(i, k, j) = pirr
               r00(i, k, j) = rhor
               r0s = sqrt(rhor)
               rr0r = 1./rhor
               rrs = sqrt(rr0r)
               rrq = sqrt(rrs)
               zr(i, k, j) = 1.e5*zrc*rrq
               zs(i, k, j) = 1.e5*zsc*rrq
               zg(i, k, j) = 1.e5*zgc*rrq
               afcp(i, k, j) = afc*pirr
               avcp(i, k, j) = avc*pirr
               ascp(i, k, j) = asc*pirr


   !            if (improve .eq. 3) then
   !           xccld=xncld*r00               !cloud number concentration
   !            end if
               if (qlwrf(i, k, j) .le. cmin) qlwrf(i, k, j) = 0.0
               if (qrwrf(i, k, j) .le. cmin) qrwrf(i, k, j) = 0.0
               if (qiwrf(i, k, j) .le. cmin) qiwrf(i, k, j) = 0.0
               if (qswrf(i, k, j) .le. cmin) qswrf(i, k, j) = 0.0
               if (qgwrf(i, k, j) .le. cmin) qgwrf(i, k, j) = 0.0
               tair(i, k, j) = (ptwrf(i, k, j) + tb0)*pir1
               
   #ifdef Readaeroclx
               ! -------------
               ! aerosol-aware :
               ! -------------

               if ((ccnflag .eq. 2) .or. (inflag .eq. 2)) then
                  !!$acc loop seq
                  !do n = 1, nlut
                  !   aerog(n, i, k, j) = 0.
                  !end do
                  mr2mc = rhor*1.e+6  !mass mixing ratio (kg/kg) to mass concentration (g/m^-3)
   #ifdef LUT_aero
                  if (naero .lt. nlut) stop 'naero must not be smaller than nlut!'
                  if (nsuso .le. nlut) aerog(nsuso, i, k, j) = &
                                          max(aeroclx(i, k, nsuso, j), 0.)*mr2mc
                  if (nsoot .le. nlut) aerog(nsoot, i, k, j) = &
                                          max(aeroclx(i, k, nsoot, j), 0.)*mr2mc
                  if (ninso .le. nlut) aerog(ninso, i, k, j) = &
                                          max(aeroclx(i, k, ninso, j), 0.)*mr2mc
                  if (nwaso .le. nlut) aerog(nwaso, i, k, j) = &
                                          max(aeroclx(i, k, nwaso, j), 0.)*mr2mc
                  if (nssam .le. nlut) aerog(nssam, i, k, j) = &
                                          max(aeroclx(i, k, nssam, j), 0.)*mr2mc
                  if (nsscm .le. nlut) aerog(nsscm, i, k, j) = &
                                          max(aeroclx(i, k, nsscm, j), 0.)*mr2mc
                  if (nminm .le. nlut) aerog(nminm, i, k, j) = &
                                          max(aeroclx(i, k, nminm, j), 0.)*mr2mc
                  if (nmiam .le. nlut) aerog(nmiam, i, k, j) = &
                                          max(aeroclx(i, k, nmiam, j), 0.)*mr2mc
                  if (nmicm .le. nlut) aerog(nmicm, i, k, j) = &
                                          max(aeroclx(i, k, nmicm, j), 0.)*mr2mc
   #else
                  if (naero .lt. 15) stop 'not enough aerosol types!'
                  ! convert from MERRA2-aerotype to GOCART-aerotype
                  aerog(1, i, k, j) = max(aeroclx(i, k, naso4, j), 0.)*mr2mc    !sulfur and its precure    (SO4)
                  aerog(2, i, k, j) = max(aeroclx(i, k, nablc, j) + &          !soot                      (BLC
                                 aeroclx(i, k, nabbc, j), 0.)*mr2mc    !                          +BBC)
                  aerog(3, i, k, j) = max(aeroclx(i, k, naobc, j), 0.)*mr2mc    !non-hygroscopic OC        (OBC)
                  aerog(4, i, k, j) = max(aeroclx(i, k, naolc, j), 0.)*mr2mc    !hygroscopic OC            (OLC)
                  aerog(5, i, k, j) = max(aeroclx(i, k, nass1, j), 0.)*mr2mc    !sea salt accumulated mode (SS1)
                  aerog(6, i, k, j) = max(aeroclx(i, k, nass2, j) + &          !sea salt coarse mode      (SS2
                                 aeroclx(i, k, nass3, j) + &          !                          +SS3
                                 aeroclx(i, k, nass4, j), 0.)*mr2mc    !                          +SS4)
                  aerog(7, i, k, j) = max(aeroclx(i, k, nadu1, j), 0.)*mr2mc    !dust mode 1               (DU1)
                  aerog(8, i, k, j) = max(aeroclx(i, k, nadu1, j), 0.)*mr2mc    !dust mode 2               (DU1)
                  aerog(9, i, k, j) = max(aeroclx(i, k, nadu1, j), 0.)*mr2mc    !dust mode 3               (DU1)
                  aerog(10, i, k, j) = max(aeroclx(i, k, nadu1, j), 0.)*mr2mc    !dust mode 4               (DU1)
                  aerog(11, i, k, j) = max(aeroclx(i, k, nadu2, j), 0.)*mr2mc    !dust mode 5               (DU2)
                  aerog(12, i, k, j) = max(aeroclx(i, k, nadu3, j), 0.)*mr2mc    !dust mode 6               (DU3)
                  aerog(13, i, k, j) = max(aeroclx(i, k, nadu4, j), 0.)*mr2mc    !dust mode 7               (DU4)
                  aerog(14, i, k, j) = max(aeroclx(i, k, nadu5, j), 0.)*mr2mc    !dust mode 8               (DU5)
   #endif
               end if
   #endif

            end do
         end do
      end do

      !$acc parallel loop gang collapse(2) async(async_id)
      do j = jts, jte
         do k = kts, kte
            !$acc loop vector
            do i = its, myim(j)
               fv(i, k, j) = sqrt(rho(i, 1, j)/rho(i, k, j))
            end do
         end do
      end do

      
      IF (IWARM .EQ. 1) THEN
         !$acc parallel loop gang collapse(2) async(async_id) private(fact_fit, &
         !$acc&         pr0, rr0r, y1, dwvp, dwv, tca, scv, cp409, r22f, r23af, &
         !$acc&         r23br, dd, pracw, praut, y2, y3, y4, pr, cnd, qsw, dm, &
         !$acc&         cnd, rtair, ssw, tairr, zrr, avcpr, pi0r, rp0r, rhor, fv0r, &
         !$acc&         vrr, ptwrfr, qvwrfr, qiwrfr, qswrfr, qgwrfr, qrwrfr, qlwrfr)
         do j = jts, jte
            do k = kts, kte
               !$acc loop vector
               do i = its, myim(j)
   !     ******************************************************************
   !     ***   Y1 : DYNAMIC VISCOSITY OF AIR (U)
   !     ***   DWV : DIFFUSIVITY OF WATER VAPOR IN AIR (PI)
   !     ***   TCA : THERMAL CONDUCTIVITY OF AIR (KA)
   !     ***   Y2 : KINETIC VISCOSITY (V)
                  ptwrfr = ptwrf(i, k, j)
                  qvwrfr = qvwrf(i, k, j)
                  qrwrfr = qrwrf(i, k, j)
                  qlwrfr = qlwrf(i, k, j)
                  qiwrfr = 0.0
                  qswrfr = 0.0
                  qgwrfr = 0.0
                  rhor = rho(i, k, j)
                  rp0r = rp0(i, k, j)
                  pi0r = pi0(i, k, j)
                  fv0r = fv(i, k, j)
                  avcpr = avcp(i, k, j)
                  zrr = zr(i, k, j)
                  pr0 = 1./(p0_mks(i, k, j)*10.0)
                  rr0r = 1./rhor
                  tairr = tair(i, k, j)
                  vrr = 0.
                  y1 = c149*tairr**1.5/(tairr + 120.)
                  dwvp = c879*pr0
                  dwv = dwvp*tairr**1.81
                  tca = c141*y1
                  scv = 1./((rr0r*y1)**.1666667*dwv**.3333333)
                  
                  ! for calculating processes related to warm rain only
                  cp409 = c409*pi_mks(i, k, j)
                  r22f = rn22*fv0r
                  r23af = rn23a*sqrt(fv0r)
                  r23br = rn23b*rhor

                  if (qrwrfr .gt. cmin1) then
                     dd = r00(i, k, j)*qrwrfr
                     y1 = dd**.25
                     zrr = zrc/y1
                     zr(i, k, j) = zrr
                  end if

                  call vtr_mks_gpu(rho_mks(i, k, j), qrwrfr, tairr, &
                     vrr, constb, consta, rhoe_s, crmin, rhowater, xnor, &
                     vrc0, vrc1, vrc2, vrc3, draimax, tnw, roqr, t0)  !in MKS
                  vrr = vrr*100.  !in CGS
                  vr(i, k, j) = vrr

   !* 21 * PRAUT   AUTOCONVERSION OF QC TO QR                        **21**
   !* 22 * PRACW : ACCRETION OF QC BY QR                             **22**
                  pracw = 0.
                  praut = 0.0
                  praut = max(rn21*(qlwrfr - bnd21), 0.0)

                  y1 = 1./zr(i, k, j)
                  y2 = y1*y1
                  y3 = y1*y2
                  y4 = r22f*qlwrfr*y3*(rn50 + rn51*y1 + &
                                                     rn52*y2 + rn53*y3)
                  pracw = max(y4, 0.0)
                  if (qrwrfr .le. cmin) pracw = 0.

   !C********   HANDLING THE NEGATIVE CLOUD WATER (QC)    ******************
                  Y1 = qlwrf(I, k, J)/D2T
                  PRAUT = MIN(Y1, PRAUT)
                  PRACW = MIN(Y1, PRACW)
                  Y1 = (PRAUT + PRACW)*D2T

                  if (qlwrfr .lt. y1 .and. y1 .ge. cmin2) then
                     y2 = qlwrfr/(y1 + cmin2)
                     praut = praut*y2
                     pracw = pracw*y2
                     qlwrfr = 0.0
                  else
                     qlwrfr = qlwrfr - y1
                  end if

                  PR = (PRAUT + PRACW)*D2T
                  qrwrfr = qrwrfr + PR

   !*****   TAO ET AL (1989) SATURATION TECHNIQUE  ***********************
                  cnd = 0.0
                  tairr = (ptwrfr + tb0)*pi0r
                  tair(i ,k, j) = tairr
                  y1 = 1./(tairr - c358)
                  qsw = rp0r*exp(c172 - c409*y1)
                  dd = cp409*y1*y1
                  dm = qvwrfr + qb0 - qsw
                  cnd = dm/(1.+avcpr*dd*qsw)
   !c    ******   condensation or evaporation of qc  ******
                  cnd = max(-qlwrfr, cnd)
                  ptwrfr = ptwrfr + avcpr*cnd
                  qvwrfr = qvwrfr - cnd
                  qlwrfr = qlwrfr + cnd

   !* 23 * ERN : EVAPORATION OF QR (SUBSATURATION)                   **23**
                  ern = 0.0

                  if (qrwrfr .gt. 0.0) then
                     tairr = (ptwrfr + tb0)*pi0r
                     tair(i ,k, j) = tairr
                     rtair = 1./(tairr - c358)
                     qsw = rp0r*exp(c172 - c409*rtair)
                     ssw = (qvwrfr + qb0)/qsw - 1.0
                     dm = qvwrfr + qb0 - qsw
                     rsub1 = cv409*qsw*rtair*rtair
                     dd1 = max(-dm/(1.+rsub1), 0.0)
                     y1 = .78/zrr**2 + r23af*scv/zrr**bwh5
                     y2 = r23br/(tca*tairr**2) + 1./(dwv*qsw)
                     ern = r23t*ssw*y1/y2
                     ern = min(dd1, qrwrfr, max(ern, 0.))

                     ! reducing evaporation rate according to fitq
                     if (qrwrfr .gt. 1.e-6) then
                        fact_fit = 0.11*(qrwrfr*1.e3)**(-1.27) + 0.98
                        ern = ern/fact_fit
                     end if

                     ptwrfr = ptwrfr - avcpr*ern
                     qvwrfr = qvwrfr + ern
                     qrwrfr = qrwrfr - ern
                  end if
                  ptwrf(i, k, j) = ptwrfr
                  qvwrf(i, k, j) = qvwrfr
                  qiwrf(i, k, j) = qiwrfr
                  qlwrf(i, k, j) = qlwrfr
                  qrwrf(i, k, j) = qrwrfr
                  qswrf(i, k, j) = qswrfr
                  qgwrf(i, k, j) = qgwrfr
               end do
            end do
         end do
      end if
                      ! part of if (iwarm.eq.1) then

   !JJS   for calculating processes related to both ice and warm rain

   !     ***   COMPUTE ZR,ZS,ZG,VR,VS,VG      *****************************
      IF (IWARM .ne. 1) THEN
         !$acc parallel loop gang collapse(2) async(async_id) private(tairc, &
         !$acc&         vgcr, dd, y1, y2, ftns, vgcf, ftng, tairr, r00r, vrr, &
         !$acc&         vsr, vgr, vir, qrwrfr, qswrfr, qgwrfr, qiwrfr)

         do j = jts, jte
            do k = kts, kte
               !$acc loop vector
               do i = its, myim(j)
                  qrwrfr = qrwrf(i, k, j)
                  qswrfr = qswrf(i, k, j)
                  qgwrfr = qgwrf(i, k, j)
                  qiwrfr = qiwrf(i, k, j)
                  tairr = tair(i, k, j)
                  tairc = tairr - t0
                  vgcr = vgc*sqrt(1./rho(i, k, j))
                  r00r = r00(i, k, j)
                  vrr = 0.
                  vsr = 0.
                  vgr = 0.
                  vir = 0.
                  ftns0(i, k, j) = 1.
                  ftng0(i, k, j) = 1.
                  if (qrwrfr .gt. cmin) then
                     dd = r00r*qrwrfr
                     y1 = sqrt(dd)
                     y2 = sqrt(y1)
                     zr(i, k, j) = zrc/y2
                     call vtr_mks_gpu(rho_mks(i, k, j), qrwrfr, tairr, &
                        vrr, constb, consta, rhoe_s, crmin, rhowater, &
                        xnor, vrc0, vrc1, vrc2, vrc3, draimax, tnw, roqr, t0)  !in MKS
                     vrr = vrr*100.  !in CGS
                  end if

                  if (qswrfr .gt. cmin) then
                     dd = r00r*qswrfr
                     y1 = dd**.25
                     ftns = 1.
                     if (improve .gt. 2) then
                        call sgmap_gpu(1, qswrfr, r00r, tairc, &
                                       ftns0(i, k, j), roqs, roqg, tns, tng)
                        ftns = ftns0(i, k, j)**0.25
                     end if
                     zs(i, k, j) = zsc/y1*ftns
                     call vts_mks_gpu(improve, rho_mks(i, k, j), qswrfr, &
                                      tairr, vsr, rhoe_s, constd, &
                                      csmin, rhosnow, xnos, constc, vsc, t0, &
                                      roqs, roqg, tns, tng, bsq)  !in MKS
                     vsr = vsr*100.  !in CGS
                  end if

                  if (qgwrfr .gt. cmin) then
                     dd = r00r*qgwrfr
                     y1 = dd**.25
                     if (ihail .eq. 1) vgcf = vgcr
                     ftng = 1.
                     if (improve .gt. 2) then
                        call sgmap_gpu(2, qgwrfr, r00r, tairc, &
                                       ftng0(i, k, j), roqs, roqg, tns, tng)
                        ftng = ftng0(i, k, j)**0.25
                     end if
                     zg(i, k, j) = zgc/y1*ftng
                     call vtg_mks_gpu(ihail, improve, rho_mks(i, k, j), &
                        qgwrfr, tairr, vgr, cgmin, rhohail, &
                        xnoh, cdrag, rhoe_s, vgc, roqs, roqg, tns, tng, bgq, t0)  !in MKS
                     vgr = vgr*100.  !in CGS
                  end if

                  call vti_mks_gpu(improve, rho_mks(i, k, j), tairr, &
                               qiwrfr, qvwrf(i, k, j), p0_mks(i, k, j), &
                               xland(i, j), vir, sat_predict, &
                               new_saturation, cimin, t0, rhoe_s, itble, thrd, &
                               c218, c580, c76)  !in MKS
                  vir = vir*100.  !in CGS

                  if (qrwrfr .le. crmin) vrr = 0.0
                  if (qswrfr .le. csmin) vsr = 0.0
                  if (qgwrfr .le. cgmin) vgr = 0.0
                  if (qiwrfr .le. cimin) vir = 0.0
                  vr(i, k, j) = vrr
                  vs(i, k, j) = vsr
                  vg(i, k, j) = vgr
                  vi(i, k, j) = vir
               end do
            end do
         end do
      end if
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
      IF (IWARM .ne. 1) THEN
         !$acc parallel loop gang collapse(2) async(async_id) private(tairc, &
         !$acc&         psfw, psfi, pihms, psaut, psaci, praci,piacr, psacw, &
         !$acc&         pwacs, qsacw, ftns, ftng, fv0r, r3f, r4f, r5f, r6f, &
         !$acc&         r12r, r22f, r34f, rn1s, bnd1, esi, dmicrons, y1, y2, &
         !$acc&         y3, y4, y5, dd, dd1, praut, pracw, pidep, it, qsw, &
         !$acc&         qsi, esw, hfact, sfact, ssi, fssi, r_nci, rr0r, r14f, &
         !$acc&         r15f, r9rf, r16rf, pgacs, dgacs, wgacs, dmicrong, &
         !$acc&         dgacw, pihmg, dgaci, wgaci, dgacr, qgacr, qgacw, wgacr, &
         !$acc&         dlt3, dlt4, pr, ps, pg, r7r, r8r, r18r, rrs, fvs, &
         !$acc&         r101r, r102rf, r191r, r192rf, qracs, pracs, psacr, &
         !$acc&         qsacr, pgaut, pgfr, temp, cpm, hlv, hlf, hls, dlt2, &
         !$acc&         prn, psn, psmlt, pgmlt, ftns0r, ftng0r, pgwet, tairr, &
         !$acc&         zrr, zsr, zgr, afcpr, pirr, rp0r, r00r, rhor, vrr, vsr, &
         !$acc&         vgr, vir, ptwrfr, qvwrfr, qlwrfr, qiwrfr, qrwrfr, qswrfr, &
         !$acc&         qgwrfr, p0r1)
         do j = jts, jte
            do k = kts, kte
               !$acc loop vector
               do i = its, myim(j)
                  ptwrfr = ptwrf(i, k, j)
                  qvwrfr = qvwrf(i, k, j)
                  qiwrfr = qiwrf(i, k, j)
                  qlwrfr = qlwrf(i, k, j)
                  qrwrfr = qrwrf(i, k, j)
                  qswrfr = qswrf(i, k, j)
                  qgwrfr = qgwrf(i, k, j)
                  p0r1 = p0_mks(i, k, j)
                  rhor = rho(i, k, j)
                  tairr = tair(i, k, j)
                  tairc = tairr - t0
                  psfw = 0.0
                  psfi = 0.0
                  pihms = 0.0
                  psaut = 0.0
                  psaci = 0.0
                  praci = 0.0
                  piacr = 0.0
                  psacw = 0.0
                  pwacs = 0.0
                  qsacw = 0.0
                  ftns = ftns0(i, k, j)
                  ftng = ftng0(i, k, j)
                  fv0r = fv(i, k, j)
                  r3f = rn3*fv0r
                  r4f = rn4*fv0r
                  r5f = rn5*fv0r
                  r6f = rn6*fv0r
                  r12r = rn12*rhor
                  r22f = rn22*fv0r
                  zrr = zr(i, k, j)
                  zsr = zs(i, k, j)
                  zgr = zg(i, k, j)
                  afcpr = afcp(i, k, j)
                  pirr = pir(i, k, j)
                  rp0r = rp0(i, k, j)
                  r00r = r00(i, k, j)
                  vrr = vr(i, k, j)
                  vsr = vs(i, k, j)
                  vgr = vg(i, k, j)
                  vir = vi(i, k, j)
                  if (improve .eq. 3) then
                     r34f = rn34*fv0r
                  end if

                  if (tairr .lt. t0) then

!                     if (sat_predict) then
                        rn1s = 1.e-3
                        bnd1 = 1.e-4
                        efsi = exp(0.025*tairc)
                        psaut = r2is*max(rn1s*efsi*(qiwrfr &
                                - bnd1*fv0r*fv0r), 0.0)
!                     else !sat_predict
   !             y1(i,j)=rdt*(qi(i,j)-r1r*exp(beta*tairc(i,j)))
   !             psaut(i,j)=max(y1(i,j),0.0)
!                        rn1s = 1.e-3
!                        bnd1 = 6.e-4
!                        esi = exp(.025*tairc)
!                        if (improve .gt. 2) esi = 0.15
!                        psaut = r2is*max(rn1s*esi*(qiwrfr &
!                                - bnd1*fv0r*fv0r), 0.0)
!                     end if !sat_predict
!                     esi = 1.0
                     dmicrons = (r00r*qswrfr &
                                /roqs/cpi/(tns*ftns))**.25*1.e4
                     fdms = min(1., (dmicrons/1500.)**4.) ! f(dmicrons)

                     y1 = 1.0
                     if (vsr .gt. 0.) y1 = abs((vsr - vir) &
                                                   /vsr)
                     psaci = y1*r3f*qiwrfr/zsr**bs3*ftns*fdms
                     psacw = r4f*qlwrfr/zsr**bs3*ftns
                     if (ihalmos .eq. 1) then
                        y2 = 0.
                        if ((tairc .le. hmtemp1) .and. (tairc .ge. hmtemp4)) &
                           y2 = 0.5
                        if ((tairc .le. hmtemp2) .and. (tairc .ge. hmtemp3)) &
                           y2 = 1.
                        pihms = psacw*y2*xnsplnt*1000.*xmsplnt
                        psacw = psacw - pihms
                     end if
                     pwacs = r34f*qlwrfr/zsr**bs6*ftns
                     y1 = 1./zrr
                     y2 = y1*y1
                     y3 = y1*y2
                     y5 = 1.0
                     if (vrr .gt. 0.) y5 = abs((vrr - vir) &
                                                   /vrr)
                     dd = y5*r5f*qiwrfr*y3*(rn50 + rn51*y1 &
                                                                + rn52*y2 + rn53*y3)
                     praci = max(dd, 0.0)
                     y4 = y3*y3
                     dd1 = y5*r6f*qiwrfr*y4*(rn60 + rn61*y1 &
                                                                 + rn62*y2 + rn63*y3)

                     piacr = max(dd1, 0.0)
                  else
                     qsacw = r4f*qlwrfr/zsr**bs3*ftns
                  end if   !tairc

   !23456789012345678901234567890123456789012345678901234567890123456789012
   !* 21 * PRAUT   AUTOCONVERSION OF QC TO QR                        **21**
   !* 22 * PRACW : ACCRETION OF QC BY QR                             **22**

                  praut = max(rn21*(qlwrfr - bnd21), 0.0)
   !            if ( ccnflag .eq. 1 ) then
   !               praut(i,j)=max(rn21*(qc(i,j)-bnd21),0.0)
   !#ifdef Readaeroclx
   !            elseif ( ccnflag .eq. 2 ) then
   !               esw(i,j) = min(0.99*p0_mks(i,k,j),esw_mks(tair(i,j))) !in MKS
   !               qsw(i,j) = 0.622*esw(i,j)/(p0_mks(i,k,j)-esw(i,j))
   !               if ( qv(i,j).gt.qsw(i,j) ) then
   !                  rhw = max(1.e-6, qv(i,j)/qsw(i,j)*100.)      !relative humidity (%)
   !                  ssrw = max(0.001, rhw - 100.e0)              !super saturation rate over water (%)
   !                  call mass2ccn(tair(i,j),ssrw,aerog,ncloud)
   ! >>> Khairoutdinov and Kogan (2000) :
   !                  if ( ncloud .gt. 0. ) then
   !                     ! convert NCCN from m^-3 to cm^-3 :
   !                     praut(i,j) = 1350.0*max(qc(i,j),0.)**2.47*     &
   !                                  (ncloud*1.e+6)**(-1.79)
   !                  else
   !                     praut(i,j) = 0.
   !                  endif
   ! >>> Liu and Daum (2004) :
   !                  ncloud = max(100., ncloud)
   !                  L_cloud = qc(i,j)*rho(i,j,k)                 !cloud water content (g cm^-3)
   !                  call auto_conversion(L_cloud,ncloud,P_liu_daum,re_liu_daum)
   !                  praut(i,j) = P_liu_daum/rho(i,j,k)           !autoconversion rate (s^-1)
   !               endif
   !#endif
   !            else
   !               stop 'ccnflag error!!!'
   !            endif
   
                  y1 = 1./zrr
                  y2 = y1*y1
                  y3 = y1*y2
                  y4 = r22f*qlwrfr*y3*(rn50 + rn51*y1 + &
                                                     rn52*y2 + rn53*y3)
                  pracw = max(y4, 0.0)

   !* 12 * PSFW : BERGERON PROCESSES FOR QS (KOENING, 1971)          **12**
   !* 13 * PSFI : BERGERON PROCESSES FOR QS                          **13**

                  pidep = 0.0

                  if (sat_predict .eq. .false.) then
      !>>> Note that Bergeron processes are concerned in saturation prediction scheme
                     if (tairr .lt. t0) then
                        y1 = max(min(tairc, -1.), -31.)
                        it = int(abs(y1))
                        y1 = rn12a(it)
                        y2 = rn12b(it)
                        y3 = rn13(it)
                        psfw = r2is*max(d2t*y1*(y2 + r12r*qlwrfr)* &
                                              qiwrfr, 0.0)
                        psfi = r2is*y3*qiwrfr
      !
                        y4 = 1./(tairr - c358)
                        y5 = 1./(tairr - c76)
                        qsw = rp0r*exp(c172 - c409*y4)
                        qsi = rp0r*exp(c218 - c580*y5)
                        if (new_saturation) then
                           esw = min(0.99*p0r1, esw_mks_gpu(tairr))
                           esi = min(0.99*p0r1, esi_mks_gpu(tairr))
                           if (esi .gt. esw) esi = esw
                           qsw = 0.622*esw/(p0r1 - esw)
                           qsi = 0.622*esi/(p0r1 - esi)
                           esw = esw*10.  !in CGS
                           esi = esi*10.  !in CGS
                        end if
                        ! EMK...Prevent division by zero
      !               hfact=(qv(i,j)+qb0-qsi(i,j))/(qsw(i,j)-qsi(i,j))
                        hfact = (qvwrfr + qb0 - qsi)/(qsw - qsi + cmin1)
                        if (hfact .gt. 1.) hfact = 1.
                        sfact = 1
                        SSI = (qvwrfr + qb0)/qsi - 1.

                        fssi = min(xssi, max(.0, xssi*(tairc + 44.)/(44.0 - 38.0)))  !max ssi f(tair)

                        fssi = min(ssi, fssi)

      !  STEVE : PLEASE CHECK
      !                 xssi=min(ssi(i,j), fssi)
      !                  xssi=min(ssi(i,j), 0.20)
                        r_nci = min(1.e-3*exp(-.639 + 12.96*fssi), 1.)  !meyers et al. 1992 (cm^-3)
      !  STEVE : PLEASE CHECK
                        if (tairc .le. -5.) then                         !meyers
                           r_nci = max(1.e-3*exp(-.639 + 12.96*fssi), 0.528e-3)  !meyers et al.
                        else
                           r_nci = min(1.e-3*exp(-.639 + 12.96*fssi), 0.528e-3)  !meyers et al.
                        end if  !tairc

                        if (r_nci .gt. 15.) r_nci = 15.                  !cap at 15000/liter

                        dd = min((r00r*qiwrfr/r_nci), ami40)   !mean cloud ice mass
                        y4 = 1.-aa2(it)
                        sfact = (AMI50**Y4 - AMI40**Y4)/ &
                                (AMI50**Y4 - dd**Y4)
                        if (hfact .gt. 0.) then
                           psfi = r2is*psfi*hfact*sfact
                        else
                           psfi = 0.
                        end if
                     end if   !tair(i,j)
                  endif !sat_predict
   !TTT***** QG=QG+MIN(PGDRY,PGWET)
   !*  9 * PGACS : ACCRETION OF QS BY QG (DGACS,WGACS: DRY AND WET)  ***9**
   !* 14 * DGACW : ACCRETION OF QC BY QG (QGACW FOR PGMLT)           **14**
   !* 16 * DGACR : ACCRETION OF QR TO QG (QGACR FOR PGMLT)           **16**
   !*******PGDRY : DGACW+DGACI+DGACR+DGACS                           ******
   !* 15 * DGACI : ACCRETION OF QI BY QG (WGACI FOR WET GROWTH)      **15**
   !* 17 * PGWET : WET GROWTH OF QG                                  **17**
   !*  Steve turned off PGWET, set PGWET = 0.
   !*  Steve turned off wgaci, set wgaci = 0.
                  rr0r = 1./rhor
                  r14f = rn14*fv0r
                  r15f = rn15*fv0r
                  if (improve .eq. 3) then
                     r9rf = rn9*fv0r*rr0r
                     r16rf = rn16*fv0r*rr0r
                  end if
                  y1 = abs(vgr - vsr)
                  y2 = zsr*zgr
                  y3 = 5./y2
                  y4 = .08*y3*y3
                  y5 = .05*y3*y4
                  y2 = y1*(y3/zsr**5 + y4/zsr**3 &
                                       + y5/zsr)

                  pgacs = r2ig*r2is*r9rf*y2*ftns*ftng
   !            dgacs(i,j)=pgacs(i,j)
                  dgacs = 0.0             !Lang et al. 2007
                  wgacs = 0.0
                  y1 = 1./zgr**bg3
!                  esi = 1.0 !egc constant in consatrh via r14f/rn14; use esi( ) to make f(T)/f(q)

                  dmicrong = (r00r*qgwrfr/roqg/cpi/(tng*ftng))**.25*1.e4
                  fdmg = min(1., (dmicrong/500.)**1.1)       ! f(dmicrons)
                  dgacw = r2ig*fdmg*r14f*qlwrfr*y1*ftng
   !            dgacw(i,j)=r2ig*r14f*qc(i,j)*y1(i,j)*ftng(i,j)
                  y2 = 0.
                  if ((tairc .le. hmtemp1) .and. (tairc .ge. hmtemp4)) &
                     y2 = 0.5
                  if ((tairc .le. hmtemp2) .and. (tairc .ge. hmtemp3)) &
                     y2 = 1.
                  pihmg = r2ig*dgacw*y2*xnsplnt*1000.*xmsplnt
                  dgacw = r2ig*dgacw - pihmg
                  qgacw = r2ig*dgacw
                  y5 = 1.0
                  if (vgr .gt. 0.) y5 = abs((vgr - vir) &
                                                /vgr)
                  dgaci = r2ig*y5*r15f*qiwrfr*y1*ftng
   !            dgaci(i,j)=r2ig*r15f*qi(i,j)*y1(i,j)*ftng(i,j)
!                  dgaci = 0.0
                  wgaci = 0.0
                  y1 = abs(vgr - vrr)
                  y2 = zrr*zgr
                  y3 = 5./y2
                  y4 = .08*y3*y3
                  y5 = .05*y3*y4
                  dd = r16rf*y1*(y3/zrr**5 + y4/zrr**3 &
                                             + y5/zrr)*ftng
                  dgacr = r2ig*max(dd, 0.0)
                  qgacr = dgacr

                  if (tair(i, k, j) .ge. t0) then
                     dgacs = 0.0
                     dgacw = 0.0
                     dgaci = 0.0
                     dgacr = 0.0
                  else
                     pgacs = 0.0
                     qgacw = 0.0
                     qgacr = 0.0
                  end if

                  pgwet = 0.0
                  !********   HANDLING THE NEGATIVE CLOUD WATER (QC)    ******************
                  y1 = qlwrfr/d2t
                  psacw = min(y1, psacw)
                  praut = min(y1, praut)
                  pracw = min(y1, pracw)
                  psfw = min(y1, psfw)
                  dgacw = min(y1, dgacw)
                  qsacw = min(y1, qsacw)
                  qgacw = min(y1, qgacw)
                  pihms = min(y1, pihms)
                  pihmg = min(y1, pihmg)

                  y1 = d2t*(psacw + praut + pracw + psfw &
                                  + dgacw + qsacw + qgacw + pihms + pihmg)

                  qlwrfr = qlwrfr - y1
   !
                  if (qlwrfr .lt. 0.0) then
                     y2 = 1.
                     if (y1 .ne. 0.) y2 = qlwrfr/y1 + 1.
                     !if (abs(y1(i, k, j)) .gt. 1e-30) y2(i, k, j) = qc(i, k, j)/y1(i, k, j) + 1.
                     psacw = psacw*y2
                     praut = praut*y2
                     pracw = pracw*y2
                     psfw = psfw*y2
                     dgacw = dgacw*y2
                     qsacw = qsacw*y2
                     qgacw = qgacw*y2
                     pihms = pihms*y2
                     pihmg = pihmg*y2
                     qlwrfr = 0.0
                  end if
   !            wgacr(i,j)=0.
                  wgacr = qgacr + qgacw
   !******** SHED PROCESS (WGACR=PGWET-DGACW-WGACI-WGACS)
   !c
   !            wgacr(i,j)=pgwet(i,j)-dgacw(i,j)-wgaci(i,j)-wgacs(i,j)
   !            y2(i,j)=dgacw(i,j)+dgaci(i,j)+dgacr(i,j)+dgacs(i,j)
   !
   !            if (pgwet(i,j).ge.y2(i,j)) then
   !               wgacr(i,j)=0.0
   !               wgaci(i,j)=0.0
   !               wgacs(i,j)=0.0
   !            else
   !               dgacr(i,j)=0.0
   !               dgaci(i,j)=0.0
   !               dgacs(i,j)=0.0
   !            endif

   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !********   HANDLING THE NEGATIVE CLOUD ICE (QI)      ******************
                  y1 = qiwrfr/d2t
                  psaut = min(y1, psaut)
                  psaci = min(y1, psaci)
                  praci = min(y1, praci)
                  psfi = min(y1, psfi)
                  dgaci = min(y1, dgaci)
                  wgaci = min(y1, wgaci)

   !       Steve: Please check

                  y1 = d2t*(psaut + psaci + praci + psfi &
                                  + dgaci + wgaci - pihms - pihmg)

                  qiwrfr = qiwrfr - y1
                  if (qiwrfr .lt. 0.0) then
                     y2 = 1.
                     if (y1 .ne. 0.0) y2 = qiwrfr/y1 + 1.
                     psaut = psaut*y2
                     psaci = psaci*y2
                     praci = praci*y2
                     psfi = psfi*y2
                     dgaci = dgaci*y2
                     wgaci = wgaci*y2
                     qiwrfr = 0.0
                  end if

   !            qi(i,j)=qi(i,j)+d2t*(pihms(i,j)+pihmg(i,j))
                  wgacr = qgacr + qgacw
                  dlt3 = 0.0
                  if (qrwrfr .lt. 1.e-4) dlt3 = 1.
                  dlt4 = 1.
   !              if (qc(i,j) .gt. 5.e-4) dlt4(i,j)=0.0
   !              if (qs(i,j) .le. 1.e-4) dlt4(i,j)=1.
                  if (qlwrfr .gt. 1.e-3) dlt4 = 0.0
                  if (qswrfr .le. 1.e-4) dlt4 = 1.

                  if (tairr .ge. t0) then
                     dlt3 = 0.0
                     dlt4 = 0.0
                  end if
                  pr = d2t*(qsacw + praut + pracw + wgacr - qgacr)
                  ps = d2t*(psaut + psaci + dlt4*psacw + psfw + psfi + dlt3*praci)
                  pg = d2t*((1.-dlt3)*praci + dgaci + wgaci + dgacw + (1.-dlt4)*psacw)
   !*  7 * PRACS : ACCRETION OF QS BY QR                             ***7**
   !*  8 * PSACR : ACCRETION OF QR BY QS (QSACR FOR PSMLT)           ***8**
   !*  2 * PGAUT : AUTOCONVERSION OF QS TO QG                        ***2**
   !* 18 * PGFR : FREEZING OF QR TO QG                               **18**
                  r7r = rn7*rr0r
                  r8r = rn8*rr0r
                  r18r = rn18*rr0r
                  if (improve .eq. 3) then
                     rrs = sqrt(rr0r)
                     fvs = sqrt(fv0r)
                     r7r = rn7*rr0r*fv0r
                     r8r = rn8*rr0r*fv0r
                     r101r = rn101*rr0r
                     r102rf = rn102*rrs*fvs
                     r191r = rn191*rr0r
                     r192rf = rn192*rrs*fvs
                  end if
                  
                  qracs = 0.0
                  y1 = abs(vrr - vsr)
                  y2 = zrr*zsr
                  y3 = 5./y2
                  y4 = .08*y3*y3
                  y5 = .05*y3*y4
                  pracs = r2ig*r2is*r7r*y1*(y3/zsr**5 &
                          + y4/zsr**3 + y5/zsr)*ftns
                  psacr = r2is*r8r*y1*(y3/zrr**5 &
                          + y4/zrr**3 + y5/zrr)*ftns
                  qsacr = psacr
                  qracs = min(d2t*pracs, qswrfr)

                  pgaut = 0.0
                  if (qswrfr .gt. 2.e-3) then
                     pgaut = r2is*max(1.e-3*exp(0.09*tairc)*qswrfr - 2.e-3, 0.0)
                  endif
                  pgfr = 0.0
                  if (tairr .lt. t0) then
                     y2 = exp(rn18a*(t0 - tairr))
                     temp = 1./zrr
                     temp = temp*temp*temp*temp*temp*temp*temp
                     pgfr = r2ig*max(r18r*(y2 - 1.)*temp, 0.0)
                  end if

                  if (tair(i, k, j) .ge. t0) then
                     pracs = 0.0
                     psacr = 0.0
                  else
                     qsacr = 0.0
                     qracs = 0.0
                  end if

   !********   HANDLING THE NEGATIVE RAIN WATER (QR)    *******************
   !********   HANDLING THE NEGATIVE SNOW (QS)          *******************
                  if (use_cpm) then
                     cpm = cp*(1.0 - qvwrfr - qlwrfr - qiwrfr - qrwrfr - qswrfr - qgwrfr) &
                           + cvap*qvwrfr &
                           + cliq*(qlwrfr + qrwrfr) &
                           + cice*(qiwrfr + qswrfr + qgwrfr)
                     hlv = alv - (cliq - cvap)*(tairr - t0)
                     hlf = alf - (cice - cliq)*(tairr - t0)
                     hls = hlv + hlf
                     avcp(i, k, j) = hlv/cpm*pirr
                     ascp(i, k, j) = hls/cpm*pirr
                     afcpr = hlf/cpm*pirr
                     afcp(i, k, j) = afcpr
                  end if
                  y1 = qrwrfr/d2t
                  piacr = min(y1, piacr)
                  dgacr = min(y1, dgacr)
                  psacr = min(y1, psacr)
                  pgfr = min(y1, pgfr)
                  y1 = (piacr + dgacr + psacr + pgfr)*d2t
                  qrwrfr = qrwrfr + pr + qracs - y1

                  if (qrwrfr .lt. 0.0) then
                     y2 = 1.
                     if (y1 .ne. 0.0) y2 = qrwrfr/y1 + 1.
                     !if (abs(y1(i, k, j)) .gt. 1e-30) y2(i, k, j) = qr(i, k, j)/y1(i, k, j) + 1.
                     piacr = piacr*y2
                     dgacr = dgacr*y2
                     pgfr = pgfr*y2
                     psacr = psacr*y2
                     qrwrfr = 0.0
                  end if
                  dlt2 = 1.
                  if (qrwrfr .gt. 1.e-4) dlt2 = 0.
                  if (qswrfr .le. 1.e-4) dlt2 = 1.
                  if (tairr .ge. t0) dlt2 = 0.
                  y1 = qswrfr/d2t
                  pgacs = min(y1, pgacs)
                  dgacs = min(y1, dgacs)
                  wgacs = min(y1, wgacs)
                  pgaut = min(y1, pgaut)
                  pracs = min(y1, pracs)
                  pwacs = min(y1, pwacs)
                  prn = d2t*((1.-dlt3)*piacr + dgacr + pgfr + (1.-dlt2)*psacr)
                  ps = ps + d2t*(dlt3*piacr + dlt2*psacr)
                  pracs = (1.-dlt2)*pracs
                  pwacs = (1.-dlt4)*pwacs
                  psn = d2t*(pgacs + dgacs + wgacs + pgaut + pracs + pwacs)

                  qswrfr = qswrfr + ps - qracs - psn
                  if (qswrfr .lt. 0.0) then
                     y2 = 1.
                     if (psn .ne. 0.) y2 = qswrfr/psn + 1.
                     !if (abs(psn(i, k, j)) .gt. 1e-30) y2(i, k, j) = qs(i, k, j)/psn(i, k, j) + 1.
                     pgacs = pgacs*y2
                     dgacs = dgacs*y2
                     wgacs = wgacs*y2
                     pgaut = pgaut*y2
                     pracs = pracs*y2
                     pwacs = pwacs*y2
                     qswrfr = 0.0
                  end if
                  psn = d2t*(pgacs + dgacs + wgacs + pgaut + pracs + pwacs)
                  qgwrfr = qgwrfr + pg + prn + psn
                  y1 = d2t*(psacw + psfw + dgacw + piacr &
                          + dgacr + psacr + pgfr + pihms + pihmg) &
                             - qracs
                  ptwrfr = ptwrfr + afcpr*y1

   !* 11 * PSMLT : MELTING OF QS                                     **11**
   !* 19 * PGMLT : MELTING OF QG TO QR                               **19**

                  psmlt = 0.0
                  pgmlt = 0.0
                  tairr = (ptwrfr + tb0)*pi0(i, k, j)
                  tair(i ,k, j) = tairr
                  tairc = tairr - t0

                  ftns = 1.
                  ftng = 1.
                  ftns0r = 1.
                  ftng0r = 1.
                  call sgmap_gpu(1, qswrfr, r00r, tairc, &
                                 ftns0r, roqs, roqg, tns, tng)
                  call sgmap_gpu(2, qgwrfr, r00r, tairc, &
                                 ftng0r, roqs, roqg, tns, tng)

                  if (tairr .ge. t0) then
                     tairc = tairr - t0

                     dd = r11t*tairc*(r101r/zsr**2 + r102rf &
                                                  /zsr**bsh5)*ftns0r
                     psmlt = r2is*min(qswrfr, max(dd, 0.0))
                     y2 = (r191r/zgr**2 + r192rf/zgr**bgh5)*ftng0r
   !               dd1(i,j)=tairc(i,j)*(r19t*y2(i,j)+r19at*(qgacw(i,j) &
   !                                             +qgacr(i,j)))*ftng0(i,j)
                     dd1 = tairc*(r19t*y2 + r19at*(qgacw &
                                                                     + qgacr))
                     pgmlt = r2ig*min(qgwrfr, max(dd1, 0.0))
                     if (use_cpm) then
                        cpm = cp*(1.0 - qvwrfr - qlwrfr - qiwrfr - qrwrfr - qswrfr - qgwrfr) &
                              + cvap*qvwrfr &
                              + cliq*(qlwrfr + qrwrfr) &
                              + cice*(qiwrfr + qswrfr + qgwrfr)
                        hlv = alv - (cliq - cvap)*(tairr - t0)
                        hlf = alf - (cice - cliq)*(tairr - t0)
                        hls = hlv + hlf
                        avcp(i, k, j) = hlv/cpm*pirr
                        ascp(i, k, j) = hls/cpm*pirr
                        afcpr = hlf/cpm*pirr
                        afcp(i, k, j) = afcpr
                     end if
                     ptwrfr = ptwrfr - afcpr*(psmlt + pgmlt)
                     qrwrfr = qrwrfr + psmlt + pgmlt
                     qswrfr = qswrfr - psmlt
                     qgwrfr = qgwrfr - pgmlt
                  end if   ! processes 11 & 19
                  ptwrf(i, k, j) = ptwrfr
                  qlwrf(i, k, j) = qlwrfr
                  qiwrf(i, k, j) = qiwrfr
                  qrwrf(i, k, j) = qrwrfr
                  qswrf(i, k, j) = qswrfr
                  qgwrf(i, k, j) = qgwrfr
               end do
            end do
         end do
      end if
   !
   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !* 24 * PIHOM : HOMOGENEOUS FREEZING OF QC TO QI (T < T00)        **24**
   !* 25 * PIDW : DEPOSITION GROWTH OF QC TO QI ( T0 < T <= T00)     **25**
   !* 26 * PIMLT : MELTING OF QI TO QC (T >= T0)                     **26**
   !****** PIMM  : IMMERSION FREEZING OF QC TO QI (T < T0)           ******
   !****** PCFR  : CONTACT NUCLEATION OF QC TO QI (T < T0)           ******
     IF (IWARM .ne. 1) THEN
         !$acc parallel loop gang collapse(2) async(async_id) private(fssi, r_nci, &
         !$acc&         xncld, esat, rv, rlapse_m, delT, xccld, Xknud, alpha, &
         !$acc&         cunnF, dvair, DIFFar, cpm, hlv, hlf, hls, pimm, pcfr, &
         !$acc&         rr0, tairc, ftns, ftng, ftns0r, ftng0r, pihom, pimlt, &
         !$acc&         pidw, y1, it, y2, y3, y4, qsw, rtair, y5, qsi, esw, &
         !$acc&         esi, ssi, dd, tairr, pirr, rp0r, r00r, ptwrfr, qvwrfr, &
         !$acc&         qlwrfr, qiwrfr, qrwrfr, qswrfr, qgwrfr, p0r1)

         do j = jts, jte
            do k = kts, kte
               !$acc loop vector
               do i = its, myim(j)
                  pimm = 0.0
                  pcfr = 0.0
                  rr0 = 1./rho(i, k, j)
                  ptwrfr = ptwrf(i, k, j)
                  qvwrfr = qvwrf(i, k, j)
                  qlwrfr = qlwrf(i, k, j)
                  qiwrfr = qiwrf(i, k, j)
                  qrwrfr = qrwrf(i, k, j)
                  qswrfr = qswrf(i, k, j)
                  qgwrfr = qgwrf(i, k, j)
                  p0r1 = p0_mks(i, k, j)
                  if (qlwrfr .le. cmin1) qlwrfr = 0.0
                  if (qiwrfr .le. cmin1) qiwrfr = 0.0
                  tairr = (ptwrfr + tb0)*pi0(i, k, j)
                  tair(i ,k, j) = tairr
                  tairc = tairr - t0
                  pirr = pir(i, k, j)
                  rp0r = rp0(i, k, j)
                  r00r = r00(i, k, j)

                  ftns = 1.
                  ftng = 1.
                  ftns0r = 1.
                  ftng0r = 1.
                  call sgmap_gpu(1, qswrfr, r00r, tairc, &
                                 ftns0r, roqs, roqg, tns, tng)
                  call sgmap_gpu(2, qgwrfr, r00r, tairc, &
                                 ftng0r, roqs, roqg, tns, tng)

                  if (tairr .le. t00) then
                     pihom = qlwrfr
                  else
                     pihom = 0.0
                  end if
                  if (tairr .ge. t0) then
                     pimlt = qiwrfr
                  else
                     pimlt = 0.0
                  end if
                  pidw = 0.0

                  if (tairr .lt. t0 .and. tairr .gt. t00) then
                     tairc = tairr - t0
                     if (sat_predict .eq. .false.) then
   !>>> pidw may be already calculated in saturation prediction scheme.
                        y1 = max(min(tairc, -1.), -31.)
                        it = int(abs(y1))
                        y2 = aa1(it)
                        y3 = aa2(it)
                        if (tairc .le. -5.) then                       !  meyers
                           y4 = 1./(tairr - c358)
                           qsw = rp0r*exp(c172 - c409*y4)
                           rtair = 1./(tairr - c76)
                           y5 = exp(c218 - c580*rtair)
                           qsi = rp0r*y5
                           if (new_saturation) then
                              esw = min(0.99*p0r1, esw_mks_gpu(tairr))
                              esi = min(0.99*p0r1, esi_mks_gpu(tairr))
                              if (esi .gt. esw) esi = esw
                              qsw = 0.622*esw/(p0r1 - esw)
                              qsi = 0.622*esi/(p0r1 - esi)
                              esw = esw*10.  !in CGS
                              esi = esi*10.  !in CGS
                           end if
                           SSI = (qvwrfr + qb0)/qsi - 1.
                           fssi = min(xssi, max(.0, xssi*(tairc + 44.)/(44.0 - 38.0))) !max ssi f(tair)
                           fssi = min(ssi, fssi)

                           r_nci = max(1.e-3*exp(-.639 + 12.96*fssi), 0.528e-3)  ! Meyers et al. 1992 (cm^-3)
                           if (r_nci .gt. 15.) r_nci = 15.                       ! cap at 15000/liter

                           dd = (r00r*qiwrfr/r_nci)**y3                  !meyers
                           PIDW = min(RR0*D2T*y2*r_nci*dd, qlwrfr) !meyers
                        end if  !tairc
                     end if !sat_predict
                     
                     if (qlwrfr .gt. 0.0) then
                        y4 = 1./(tairr - c358)
                        qsw = rp0r*exp(c172 - c409*y4)
                        if (new_saturation) then
                           esw = min(0.99*p0r1, esw_mks_gpu(tairr))
                           qsw = 0.622*esw/(p0r1 - esw)
                           esw = esw*10.  !in CGS
                        end if
                        xncld = qlwrfr/4.e-9                         !cloud number
                        esat = 0.6112*exp(17.67*tairc/(tairc + 243.5))*10.
                        rv = 0.622*esat/((p0r1*10.0)/1000.-esat)
                        rlapse_m = 980.616*(1.+2.5e6*rv/287./tairr) &
                                   /(1004.67 + 2.5e6*2.5e6*rv*0.622 &
                                   /(287.*tairr*tairr))
                        delT = rlapse_m*w_mks(i, k, j)
                        if (delT .lt. 0.) delT = 0.
   !  STEVE: PLAESE CHECK (2 4.e-9)
   !                  pimm(i,j)=xncld*Bhi*4.e-9*exp(-tairc(i,j))*delT*d2t*4.e-9
                        pimm = xncld*Bhi*4.e-9*exp(-tairc)*delT*d2t*4.e-9

                        xccld = xncld*r00r                           !cloud number concentration
                        Xknud = 7.37*tairr/(288.*Ra*p0r1*10.0)  ! Knudsen number
                        alpha = 1.257 + 0.400*exp(-1.10/Xknud)     ! Cunningham correction (P&Klett)

                        cunnF = 1.+alpha*Xknud                   ! Cunningham correction (P&Klett)

                        if (tairc .ge. 0.) then
                           dvair = (1.718 + 0.0049*tairc)*1.e-4
                           !dynamic visc air (Prupp&Klett)
                        else
                           dvair = (1.718 + 0.0049*tairc - 1.2e-5*tairc**2)*1.e-4
                           !dynamic visc air (Prupp&Klett)
                        end if
                        DIFFar = 1.3804e-16*tairr/6./cpi/dvair/Ra*cunnF     !aerosol diffusion via P&Klett
                        if (qvwrfr + qb0 - qsw .lt. 0.) then                  !only when cloud evaporating
                           pcfr = 4.e-9*4.*cpi*Rc*DIFFar*xccld*Cna*rr0*d2t    !Brownian part only via Cotton
                        end if
                     end if  !qc
                  end if  !tair

   !      STEVE: PLEASE CHECK
                  y1 = pihom + pidw + pimm + pcfr - pimlt

                  if (y1 .gt. qlwrfr) then
                     y1 = qlwrfr
                     y2 = 1.
                     y3 = pihom + pidw + pimm + pcfr
                     if (y3 .ne. 0.) y2 = (qlwrfr + pimlt)/y3
                     pihom = pihom*y2
                     pidw = pidw*y2
                     pimm = pimm*y2
                     pcfr = pcfr*y2
                  end if  !y1

                  if (use_cpm) then
                     cpm = cp*(1.0 - qvwrfr - qlwrfr - qiwrfr - qrwrfr - qswrfr - qgwrfr) &
                           + cvap*qvwrfr &
                           + cliq*(qlwrfr + qrwrfr) &
                           + cice*(qiwrfr + qswrfr + qgwrfr)
                     hlv = alv - (cliq - cvap)*(tairr - t0)
                     hlf = alf - (cice - cliq)*(tairr - t0)
                     hls = hlv + hlf
                     avcp(i, k, j) = hlv/cpm*pirr
                     ascp(i, k, j) = hls/cpm*pirr
                     afcp(i, k, j) = hlf/cpm*pirr
                  end if
                  ptwrf(i, k, j) = ptwrfr + afcp(i, k, j)*y1
                  qlwrf(i, k, j) = qlwrfr - y1
                  qiwrf(i, k, j) = qiwrfr + y1
               end do
            end do
         end do
      end if

      IF (IWARM .ne. 1) THEN
         if (sat_predict) then
#ifdef Readaeroclx
            !$acc parallel loop gang collapse(2) async(async_id) private(cpm, &
            !$acc&         hlv, C1, K1, rhw, nact, qcmax, hlf, hls, r_nci, &
            !$acc&         qimax, pact, pint, cnd, dep, fez, ern, tairc, &
            !$acc&         rhoair, cpm1, xlv, esw, qsw, y1, abw, xls, esi, &
            !$acc&         qsi, y2, abi, ssi, rp0r, xlandr, qvwrfr)
#else
            !$acc parallel loop gang collapse(2) async(async_id) private(cpm, &
            !$acc&         hlv, C1, K1, nact, qcmax, hlf, hls, r_nci, &
            !$acc&         qimax, pact, pint, cnd, dep, fez, ern, tairc, &
            !$acc&         rhoair, cpm1, xlv, esw, qsw, y1, abw, xls, esi, &
            !$acc&         qsi, y2, abi, ssi, nice, ncloud, tairr, rp0r, xlandr, &
            !$acc&         qvwrfr, qlwrfr, qiwrfr, qrwrfr, qswrfr, qgwrfr, p0r1)
#endif
            do j = jts, jte
               do k = kts, kte
#ifdef Readaeroclx
                  !$acc loop vector private(p_mb, nice, ncloud, ssrw, tairr, qlwrfr, &
                  !$acc&     qiwrfr, qrwrfr, qswrfr, qgwrfr, p0r1)
#endif
                  do i = its, myim(j)
                     ! --------------------------------------------------------------------------------
                     ! saturation prediction from TCWA 1-moment scheme (Morrison and Milbrandt 2015) :
                     !  (in MKS system)
                     ! --------------------------------------------------------------------------------

                     pact = 0.0
                     pint = 0.0
                     cnd = 0.0
                     dep = 0.0
                     fez = 0.0
                     ern = 0.0

                     ! -------------
                     ! pact : cloud water activation
                     ! -------------

                     tairr = (ptwrf(i, k, j) + tb0)*pi0(i, k, j)
                     tair(i ,k, j) = tairr
                     tairc = tairr - t0
                     rhoair = rho_mks(i, k, j)
                     rp0r = rp0(i, k, j)
                     p0r1 = p0_mks(i, k, j)
                     xlandr = xland(i, j)
                     qvwrfr = qvwrf(i, k, j)
                     qlwrfr = qlwrf(i, k, j)
                     qiwrfr = qiwrf(i, k, j)
                     qrwrfr = qrwrf(i, k, j)
                     qswrfr = qswrf(i, k, j)
                     qgwrfr = qgwrf(i, k, j)
                     if (use_cpm) then
                        cpm = cp*(1.0 - qvwrfr - qlwrfr - qiwrfr - qrwrfr - qswrfr - qgwrfr) &
                              + cvap*qvwrfr &
                              + cliq*(qlwrfr + qrwrfr) &
                              + cice*(qiwrfr + qswrfr + qgwrfr)      ! specific heat capacity (in CGS)
                        cpm1 = cpm*1.e-4                          ! specific heat capacity (in MKS)
                        hlv = alv - (cliq - cvap)*(tairr - t0)    ! latent heat of vaporization (in CGS)
                        xlv = hlv*1.e-4                           ! latent heat of vaporization (in MKS)
                     else
                        cpm1 = cp*1.e-4*(1. + 0.887*qvwrfr)       ! specific heat capacity (in MKS)
                        xlv = 3.1484E6 - 2370.*tairr           ! latent heat of vaporization (in MKS)
                     end if
                     if (new_saturation) then
                        esw = min(0.99*p0r1, esw_mks_gpu(tairr))
                        qsw = 0.622*esw/(p0r1 - esw)
                     else
                        y1 = 1./(tairr - c358)
                        qsw = rp0r*exp(c172 - c409*y1)
                     end if
                     abw = 1.+xlv**2.*qsw/cpm1/(4.61495E2*tairr**2.)  ! abw=1+dqsdT*(Lv/Cp)

                     if (ccnflag .eq. 1) then
                        if (qvwrfr .gt. qsw .and. w_mks(i, k, j) .gt. 1.e-3) then
                           ! C1 and K1 over land and ocean (Roger and Yau)
                           if (xlandr .eq. 1.) then  !land
                              C1 = 1000.; K1 = 0.5
                           else
                              C1 = 120.; K1 = 0.4         !ocean
                           end if
                           ncloud = min(1.E9, 1.E6*0.88*C1**(2./(K1 + 2.))* &
                                        (7.E-2*w_mks(i, k, j)**1.5)**(K1/(K1 + 2.)))
                        else
                           ncloud = 0.0
                        end if
      #ifdef Readaeroclx
                     elseif (ccnflag .eq. 2) then
                        if (qvwrfr .gt. qsw) then
                           rhw = max(1.e-6, qvwrfr/qsw*100.)       !relative humidity (%)
                           ssrw = max(0.001, rhw - 100.0)                  !super saturation rate over water (%)
                           call mass2ccn_gpu(tairr, ssrw, aerog(1, i, k, j), &
                              ncloud, lut_ccn, pts_t, pts_s, mxpts_t, mxpts_s, &
                              nlut, lut_ccn_interp(1, i, k, j))
                           ncloud = ncloud*1.e+6                        !CCN (convert from cm^-3 to m^-3)
                        else
                           ncloud = 0.0
                        end if
      #endif
                     else
                        stop 'ccnflag error!!!'
                     end if
                     if (ncloud .gt. 0.) then
                        nact = max(0., ncloud/rhoair - qlwrfr/5.236E-13)      ! CONVERT FROM CM-3 TO M-3 and 5 microm in radius
                        qcmax = max((qvwrfr - qsw), 0.)/abw
                        pact = min(max(nact*1.414E-14, 0.), qcmax)       ! 1.5 micron in radius
                        tairr = tairr + pact*xlv/cpm1
                        tair(i ,k, j) = tairr
                        tairc = tairr - t0
                        qvwrfr = max(0., qvwrfr - pact)
                        qlwrfr = max(0., qlwrfr + pact)
                     end if

                     ! -------------
                     ! pint : cloud ice initialization
                     ! -------------
                     if (use_cpm) then
                        cpm = cp*(1.0 - qvwrfr - qlwrfr - qiwrfr - qrwrfr - qswrfr - qgwrfr) &
                              + cvap*qvwrfr &
                              + cliq*(qlwrfr + qrwrfr) &
                              + cice*(qiwrfr + qswrfr + qgwrfr)      ! specific heat capacity (in CGS)
                        cpm1 = cpm*1.e-4                          ! specific heat capacity (in MKS)
                        hlv = alv - (cliq - cvap)*(tairr - t0)
                        hlf = alf - (cice - cliq)*(tairr - t0)
                        hls = hlv + hlf                           ! latent heat of sublimation (in CGS)
                        xls = hls*1.e-4                           ! latent heat of sublimation (in MKS)
                     else
                        cpm1 = cp*1.e-4*(1. + 0.887*qvwrfr)       ! specific heat capacity (in MKS)
                        xls = 3.15E6 - 2370.*tairr + 0.3337E6    ! latent heat of sublimation  (in MKS)
                     endif
                     if (new_saturation) then
                        esw = min(0.99*p0r1, esw_mks_gpu(tairr))
                        esi = min(0.99*p0r1, esi_mks_gpu(tairr))
                        if (esi .gt. esw) esi = esw
                        qsw = 0.622*esw/(p0r1 - esw)
                        qsi = 0.622*esi/(p0r1 - esi)
                     else
                        y1 = 1./(tairr - c358)
                        qsw = rp0r*exp(c172 - c409*y1)
                        y2 = 1./(tairr - c76)
                        qsi = rp0r*exp(c218 - c580*y2)
                     end if
                     abi = 1.+xls**2.*qsi/cpm1/(4.61495E2*tairr**2.)  ! abi=1+dqidT*(Ls/Cp)

                     if (qvwrfr .gt. qsi .and. tairr .lt. t0) then
                        if (inflag .eq. 1) then       ! Meyers et al. 1992 (m^-3)
                           ssi = qvwrfr/qsi - 1.
                           nice = 1.e3*exp(1.296E+1*ssi - 6.39E-1)  ! IN (convert from L^-1 to m^-3)
      #ifdef Readaeroclx
                        elseif (inflag .eq. 2) then   ! GOCART mass2icn
                           p_mb = p0r1*10.0*1.e-3                       ! pressure (hPa, equal to mbar)
                           call mass2icn_gpu(p_mb, tairr, &
                              aerog(1, i, k, j), nice, nlut, lut_in_025)
                           nice = nice*1.e+3                            ! IN (convert from L^-1 to m^-3)
!                           nice = min(nice, rhoair*qiwrfr/4.71E-10)     ! cap IN to the amount corresponding to rhoi=900, Ri=50
!                           nice = min(nice, rhoair*qiwrfr/3.77E-9)      ! cap IN to the amount corresponding to rhoi=900, Ri=100
      #endif
                        elseif (inflag .eq. 3) then  ! Hong et al. 2004
                           nice = 5.38e+7*exp(0.75*log(qiwrfr*rhoair))   ! IN (m^-3)
                        elseif (inflag .eq. 4) then  ! Cooper curve ; Chern et al. 2016
                           nice = 5.e-3*exp(0.304*abs(max(tairc, -40.0)))*1.e+3   ! IN (m^-3)
                           nice = min(nice, rhoair*qiwrfr/3.77e-9)      ! cap IN to the amount corresponding to rhoi=900, Ri=100
                        elseif (inflag .eq. 5) then  ! Fletcher 1962
                           nice = min(cn0*exp(beta*tairc), 1.0)*1.e+3
                        else
                           stop 'inflag error!!!'
                        end if
      !              r_nci = max(0.,nice/rhoair-qi(i,j)/4.19E-10)          ! RHOI = 800; DI = 1.E-4
                        if (xlandr .eq. 1.) then  !land
                           r_nci = max(0., nice/rhoair - qiwrfr/2.28E-10)      ! RHOI = 850; DI = 8.E-5
                        else
                           r_nci = max(0., nice/rhoair - qiwrfr/4.19E-10)      ! RHOI = 800; DI = 1.E-4
                        end if
                        qimax = max((qvwrfr - qsi), 0.)/abi
                        pint = min(max(0., r_nci*1.02e-13), qimax)
                        tairr = tairr + pint*xls/cpm1
                        tair(i ,k, j) = tairr
                        qvwrfr = max(0., qvwrfr - pint)
                        qiwrfr = max(0., qiwrfr + pint)
                     end if
                     qvwrf(i, k, j) = qvwrfr
                     qlwrf(i, k, j) = qlwrfr
                     qiwrf(i, k, j) = qiwrfr
                  end do
               end do
            end do
         end if
      end if
                     ! -------------
                     ! cnd : condensation of qv to qc
                     ! dep : deposition of qv to qi
                     ! fez : freezing of qc to qi
                     ! ern : evaporation of qr
                     ! -------------
      IF (IWARM .ne. 1) THEN
         if (sat_predict) then
            !$acc parallel loop gang collapse(2) async(async_id) private(cpm, hlv, &
            !$acc&         hlf, hls, ltk, lqc, ltk2, lqc2, mvdc, mvrc, ncloud, &
            !$acc&         lqi, lqi2, mvdi, mvri, hid, inhgr, rhoi, nice, lqr, &
            !$acc&         lqr2, kmin, kmax, mvdr, efdr, kdxr, afar, lzr, tnr, &
            !$acc&         avr, bvr, mur, rhoaj, gr2, gbr25, cnd1, cnd2, dep1, &
            !$acc&         dep2, fez1, fez2, ern1, ern2, latr, tau, atem, dltd, &
            !$acc&         latint, sbd, nbd, sbd1, sbd2, nbd1, nbd2, fxlat, &
            !$acc&         cpm1, xlv, xlf, xls, esw, esi, qsw, qsi, y1, y2, &
            !$acc&         pr0, rhoair, dv1, abw, abi, tauc, tairc, ssi, taui, &
            !$acc&         taur, cnd, dep, fez, ern, tairr, rp0r, xlandr, qvwrfr, &
            !$acc&         qlwrfr, qiwrfr, qrwrfr, qswrfr, qgwrfr, p0r1)
            do j = jts, jte
               do k = kts, kte
                  !$acc loop vector
                  do i = its, myim(j)
                     p0r1 = p0_mks(i, k, j)
                     pr0 = 1./(p0r1*10.0)
                     rhoair = rho_mks(i, k, j)
                     tairr = tair(i, k, j)
                     rp0r = rp0(i, k, j)
                     xlandr = xland(i, j)
                     qvwrfr = qvwrf(i, k, j)
                     qlwrfr = qlwrf(i, k, j)
                     qiwrfr = qiwrf(i, k, j)
                     qrwrfr = qrwrf(i, k, j)
                     qswrfr = qswrf(i, k, j)
                     qgwrfr = qgwrf(i, k, j)
                     if (use_cpm) then
                        cpm = cp*(1.0 - qvwrfr - qlwrfr - qiwrfr - qrwrfr - qswrfr - qgwrfr) &
                              + cvap*qvwrfr &
                              + cliq*(qlwrfr + qrwrfr) &
                              + cice*(qiwrfr + qswrfr + qgwrfr)      ! specific heat capacity (in CGS)
                        cpm1 = cpm*1.e-4                          ! specific heat capacity (in MKS)
                        hlv = alv - (cliq - cvap)*(tairr - t0)    ! latent heat of vaporization (in CGS)
                        xlv = hlv*1.e-4                           ! latent heat of vaporization (in MKS)
                        hlf = alf - (cice - cliq)*(tairr - t0)    ! latent heat of fusion (in CGS)
                        xlf = hlf*1.e-4                           ! latent heat of fusion (in MKS)
                        hls = hlv + hlf                           ! latent heat of sublimation (in CGS)
                        xls = hls*1.e-4                           ! latent heat of sublimation (in MKS)
                     else
                        cpm = cp*(1. + 0.887*qvwrfr)     ! specific heat capacity (in CGS)
                        cpm1 = cpm*1.e-4                 ! specific heat capacity (in MKS)
                        xlv = 3.1484E6 - 2370.*tairr          ! latent heat of vaporization (in MKS)
                        xls = 3.15E6 - 2370.*tairr + 0.3337E6   ! latent heat of sublimation  (in MKS)
                        xlf = alf*1.e-4                         ! latent heat of fusion (in MKS)
                     end if
                     if (new_saturation) then
                        esw = min(0.99*p0r1, esw_mks_gpu(tairr))
                        esi = min(0.99*p0r1, esi_mks_gpu(tairr))
                        if (esi .gt. esw) esi = esw
                        qsw = 0.622*esw/(p0r1 - esw)
                        qsi = 0.622*esi/(p0r1 - esi)
                     else
                        y1 = 1./(tairr - c358)
                        qsw = rp0r*exp(c172 - c409*y1)
                        y2 = 1./(tairr - c76)
                        qsi = rp0r*exp(c218 - c580*y2)
                     end if
                     dv1 = 8.794E-4*pr0*tairr**1.81

                     ! psychrometric correction to condensation/evaporation, abw=1+dqsdT*(Lv/Cp) ; dqsdT=Lv*qsw/Rv/T^2
                     abw = 1.+xlv**2.*qsw/cpm1/(4.61495E2*tairr**2.)
                     ! psychrometric correction to deposition/sublimation,   abi=1+dqidT*(Ls/Cp) ; dqidT=Ls*qsi/Rv/T^2
                     abi = 1.+xls**2.*qsi/cpm1/(4.61495E2*tairr**2.)

                     ! tauc : supersaturation relaxation timescale of cloud water
                     if (qlwrfr .ge. cwmin) then
                        ltk = log(tairr)
                        lqc = -1.*log(qlwrfr*rhoair)
                        ltk2 = ltk*ltk
                        lqc2 = lqc*lqc
                        if (xlandr .eq. 1.) then
                           mvdc = exp(5.8936819 - 7.013013*ltk &
                                      + 1.3178721*lqc + 1.1741987*ltk2 &
                                      + 2.6110916E-3*lqc2 - 0.26646396*ltk*lqc)
                        else
!                           mvdc = exp(173.57305 - 64.370929*ltk &
                           mvdc = exp(173.27305 - 64.370929*ltk &   ! decrease diameter
                                      + 0.36833626*lqc + 6.1389254*ltk2 &
                                      + 5.5915321E-3*lqc2 - 0.12488698*ltk*lqc)
                        end if
                        mvrc = max(5.e-7, min(5.e-5, mvdc/2.e+6))          ! volume-weighted mean radius of cloud water (m)
                        ncloud = 3.*qlwrfr*rhoair/(4.*cpi*1000.*mvrc**3.)  ! number concentration of cloud water (m^-3)
                        ncloud = min(1.e+9, max(0.1, ncloud))
                        tauc = 1./(4.*cpi*dv1*ncloud*mvrc)                  ! tauc=1/(4*pi*Dv*Nc*rc)
                     else
                        tauc = 1.e+10
                     end if
                     ! taui : supersaturation relaxation timescale of cloud ice
                     if (qiwrfr .ge. cimin) then
                        ltk = log(tairr)
                        lqi = -1.*log(qiwrfr*rhoair)
                        ltk2 = ltk*ltk
                        lqi2 = lqi*lqi
                        if (xlandr .eq. 1.) then
!                           mvdi = exp(154.88767 - 0.25772005*lqi &
                           mvdi = exp(155.18767 - 0.25772005*lqi &   ! increase diameter
                                      + 3.75543E-4*lqi2 - 54.560357*ltk &
                                      + 5.1248879*ltk2)
                        else
!                           mvdi = exp(155.33396 - 0.25772005*lqi &
!                           mvdi = exp(155.03396 - 0.25772005*lqi &   ! decrease diameter
                           mvdi = exp(155.23396 - 0.25772005*lqi &   ! slightly decrease diameter
                                      + 3.75543E-4*lqi2 - 54.560357*ltk &
                                      + 5.1248879*ltk2)
                        end if
                        mvri = min(5.e-4, max(3.e-6, mvdi/2.e+7))   ! volume-weighted mean radius of cloud ice (m)
                        tairc = tairr - t0
                        if (tairc .ge. -40.0) then
                           hid = max(min(nint(abs(tairc)/0.25), 120), 0)
                           inhgr = itble(hid)
                           rhoi = 900.*exp(-3.*max((qvwrfr - qsi) - 5.e-5, 0.) &
                                           /inhgr)
                        else
                           inhgr = 5.
                           ssi = qvwrfr/qsi - 1.
                           if (ssi .gt. 0.267) then
                              rhoi = -1027.456*ssi + 1185.834
                           else
                              rhoi = -32.332*ssi + 900.
                           end if
                        end if
                        rhoi = min(max(rhoi, 50.), 900.)
                        nice = 3.*qiwrfr*rhoair/(4.*cpi*rhoi*mvri**3.)  ! number concentration of cloud ice (m^-3)
                        nice = min(1.e+7, max(0.1, nice))
                        taui = 1./(4.*cpi*dv1*nice*mvri)                    ! taui=1/(4*pi*Dv*Ni*ri)
                     else
                        taui = 1.e+10
                     end if

                     ! taur : supersaturation relaxation timescale of rain
                     if (qrwrfr .ge. crmin) then
                        ltk = log(tairr)
                        lqr = -1.*log(qrwrfr*rhoair)
                        ltk2 = ltk*ltk
                        lqr2 = lqr*lqr
                        kmin = 0.2222
                        kmax = 0.8396
                        mvdr = exp(33.999045 - 13.813047*ltk &
                                   + 0.29834969*lqr + 1.6622025*ltk2 &
                                   + 1.2480129E-3*lqr2 - 0.104641*ltk*lqr)
                        efdr = exp(1.6855156 + 0.84147302*log(mvdr) &
                                   - 0.46783369*log(1.e+3) &
                                   + 1.005305E-2*(log(mvdr))**2. &
                                   +3.4833332E-2*(log(1.e+3))**2. &
                                   +1.4727223E-2*log(mvdr)*log(1.e+3))
                        efdr = max(min(efdr, mvdr/kmin**thrd), mvdr/kmax**thrd)
                        kdxr = max(kmin, min(kmax, (mvdr/efdr)**3.))
                        afar = max(0., (6.*kdxr - 3.+sqrt(8.*kdxr + 1.))/(2.-2.*kdxr))
                        lzr = log((afar + 3.)/efdr*1.e+6)
                        tnr = log(6.*qrwrfr*rhoair/cpi/1.e+3) & ! slope parameter for rain
                              + (4.+afar)*lzr - log_gamma(afar + 4.)
                        avr = exp(7.6004532 - 0.7990953*ltk &
                                  + 1.0281818*lqr - 0.16595505*lqr2 &
                                  + 1.110037E-2*lqr*lqr2 &
                                  - 2.0925743E-4*lqr2*lqr2)
                        bvr = max(0.5, min(2., 1.2218728 - 0.1281004*ltk &
                                           + 2.6088596E-2*lqr - 7.4467639E-3*lqr2 &
                                           + 7.7592532E-4*lqr*lqr2 &
                                           - 1.7056075E-5*lqr2*lqr2))
                        mur = 1.496E-6*tairr**1.5/(tairr + 120.)           ! shape parameter for rain
                        rhoaj = sqrt(1.29/rhoair)
                        gr2 = log_gamma(afar + 2.)
                        gbr25 = log_gamma(bvr*0.5 + afar + 2.5)
                        taur = 1./(2.*cpi*dv1*(0.78*exp(tnr + gr2 - (afar + 2.) &
                               *lzr) + 0.31*sqrt(avr*rhoaj/mur)*(mur/dv1) &
                               **thrd*exp(tnr + gbr25 - (bvr*0.5 + afar + 2.5) &
                               *lzr)))
                     else
                        taur = 1.e+10
                     end if

                     cnd1 = 0.; cnd2 = 0.; dep1 = 0.; dep2 = 0.
                     fez1 = 0.; fez2 = 0.; ern1 = 0.; ern2 = 0.
                     latr = 0.

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

                     ! tau : multiphase supersaturation relaxation time scale defined by
                     ! 1/tau = 1/tauc + 1/taur + (1+Ls/Cp*dqsl/dT)/(1+Ls/Cp*dqsi/dT)/taui
                     tau = 1./(1./tauc + 1./taur + 1./taui*(1.+qsw* &
                           xls*xlv/cpm1/(4.61495E+2*tairr**2.))/abi)

                     ! atem : change in supersaturation due to Bergeron process
                     atem = (qsi - qsw)/taui*(1.+qsw*xls*xlv &
                            /cpm1/(4.61495E+2*tairr**2.))/abi

                     ! decide values at lower latitude :
                     cnd1 = max(-qlwrfr, (atem*dt*tau/tauc + (qvwrfr &
                            - qsw - atem*tau)*(1.-exp(-dt/tau))*tau /tauc)/abw)
                     dep1 = MAX(-qiwrfr, (atem*dt*tau/taui + (qvwrfr &
                            - qsw - atem*tau)*(1.-exp(-dt/tau))*tau &
                            /taui + (qsw - qsi)*dt/taui)/abi)
                     ern1 = MAX(-qrwrfr, (atem*dt*tau/taur + (qvwrfr &
                            - qsw - atem*tau)*(1.-exp(-dt/tau))*tau/taur)/abw)
                     fez1 = 0.

                     if (use_declination) then
                        ! decide values at polar regions :
                        ern2 = ern1
                        if (tairr .ge. 253.16) then
                           ! T>=-20     : all ice melt; all excess vapor (over ice) condense; no deposition
                           cnd2 = max(-qlwrfr, (qvwrfr - qsi)/abi)
                           fez2 = -qiwrfr
                           dep2 = 0.
                        elseif (tairr .lt. 253.16 &
                                .and. tairr .ge. t00) then
                           ! -40<=T<-20 : all water freeze; all excess vapor (over water) deposite; no condensation
                           dep2 = max(-qiwrfr, (qvwrfr - qsw)/abw)
                           fez2 = qlwrfr
                           cnd2 = 0.
                        else
                           ! T<-40      : all water freeze; all excess vapor (over ice) deposite; no condensation
                           dep2 = max(-qiwrfr, (qvwrfr - qsi)/abi)
                           fez2 = qlwrfr
                           cnd2 = 0.
                        end if

                        ! create transition zones :
                        dltd = asin(sdec)*180.0/cpi      ! solar declination angle (in degree latitude)
                        latint = 23.5                    ! width of the transition zone (in degree latitude)
                        sbd = dltd - 90.0                  ! southern boundary of sun (in degree latitude)
                        nbd = dltd + 90.0                  ! northern boundary of sun (in degree latitude)
                        sbd1 = sbd + latint/2.0            ! southern boundary of southern transition zone
                        sbd2 = sbd - latint/2.0            ! northern boundary of southern transition zone
                        nbd1 = nbd - latint/2.0            ! southern boundary of northern transition zone
                        nbd2 = nbd + latint/2.0            ! northern boundary of northern transition zone

                        if (xlat(j) .gt. sbd1 .and. xlat(j) .lt. nbd1) then       ! low-latitude
                           fxlat = 0.
                        elseif (xlat(j) .lt. sbd2 .or. xlat(j) .gt. nbd2) then    ! polar region
                           fxlat = 1.
                        elseif (xlat(j) .ge. sbd2 .and. xlat(j) .le. sbd1) then   ! southern transition zone
                           fxlat = sin(abs(max((xlat(j) - sbd1)*90.0/latint, -90.0) &
                                   *cpi/180.0))
                        elseif (xlat(j) .ge. nbd1 .and. xlat(j) .le. nbd2) then   ! northern transition zone
                           fxlat = sin(min((xlat(j) - nbd1)*90.0/latint, 90.0)*cpi/180.0)
                        end if

                        cnd = cnd1*(1.0 - fxlat) + cnd2*fxlat
                        dep = dep1*(1.0 - fxlat) + dep2*fxlat
                        fez = fez1*(1.0 - fxlat) + fez2*fxlat
                        ern = ern1*(1.0 - fxlat) + ern2*fxlat

                     else
                        cnd = cnd1
                        dep = dep1
                        fez = fez1
                        ern = ern1
                     end if

                     ! update tair, qv, qc, qi, qr :
                     tairr = tair(i, k, j) + (ern + cnd)*xlv/cpm1 &
                                  + dep*xls/cpm1 + fez*xlf/cpm1
                     tair(i, k, j) = tairr
                     qvwrf(i, k, j) = max(0., qvwrf(i, k, j) - cnd - dep - ern)
                     qlwrf(i, k, j) = max(0., qlwrf(i, k, j) + cnd - fez)
                     qiwrf(i, k, j) = max(0., qiwrf(i, k, j) + dep + fez)
                     qrwrf(i, k, j) = max(0., qrwrf(i, k, j) + ern)

                     ! refreeze at mid-level polar regions :
      !           if ( abs(xlat).ge.66.0 .and.                             &
      !                p0(i,j,k).ge.8.e+5 .and.                            &
      !                qc(i,j).ge.cwmin ) then
      !              latr = min(1.,max(0.,abs(xlat)-66.0))
      !              fez(i,j) = qc(i,j)*latr
      !              tair(i,j) = tair(i,j) + fez(i,j)*alf/cpm
      !              qi(i,j) = max(0.,qi(i,j)+fez(i,j))
      !              qc(i,j) = max(0.,qc(i,j)-fez(i,j))
      !           endif

                     ! small cloud water all refreeze at polar regions,
                     ! in order to make the distribution more continuous :
      !           if ( abs(xlat).ge.66.5 .and. qc(i,j).lt.1.e-4 ) then
      !              tair(i,j) = tair(i,j)+qc(i,j)*alf/cpm
      !              qi(i,j) = max(0.,qi(i,j)+qc(i,j))
      !              qc(i,j) = 0.
      !           endif

                     ! update pt :
                     ptwrf(i, k, j) = tairr/pi0(i, k, j) - tb0   !tair=(pt(i,j)+tb0)*pi0
                     tairc = tairr - t0
                  end do
               end do
            end do
         end if
      end if

      !else !sat_predict
   !>>> if not define sat_predict, use saturation adjustment scheme :

   !* 31 * pint  : initiation of qi                                  **31**
   !* 32 * pidep : deposition of qi                                  **32**
   !
   !     CALCULATION OF PINT USES DIFFERENT VALUES OF THE INTERCEPT AND SLOPE FOR
   !     THE FLETCHER EQUATION. ALSO, ONLY INITIATE MORE ICE IF THE NEW NUMBER
   !     CONCENTRATION EXCEEDS THAT ALREADY PRESENT.
      
      IF (IWARM .ne. 1) THEN
         if (sat_predict .eq. .false.) then
            !$acc parallel loop gang collapse(2) async(async_id) private(cpm, &
            !$acc&         hlv, hlf, hls, fssi, r_nci, cp409, cp580, r32rt, &
            !$acc&         rtair, y3, dd, dm, rsub1, y4, tairc, y2, qsi, esi, &
            !$acc&         ssi, y1, dep, qsw, pidep, pint, cnd, esw, dd1, y5, &
            !$acc&         qvs, tairr, ascpr, avcpr, pi0r, pirr, rp0r, ptwrfr, &
            !$acc&         qvwrfr, qlwrfr, qiwrfr, qrwrfr, qswrfr, qgwrfr, pir1, &
            !$acc&         p0r1)
            do j = jts, jte
               do k = kts, kte
                  !$acc loop vector
                  do i = its, myim(j)
                     ptwrfr = ptwrf(i, k, j)
                     qvwrfr = qvwrf(i, k, j)
                     qlwrfr = qlwrf(i, k, j)
                     qiwrfr = qiwrf(i, k, j)
                     qrwrfr = qrwrf(i, k, j)
                     qgwrfr = qgwrf(i, k, j)
                     qswrfr = qswrf(i, k, j)
                     pir1 = pi_mks(i, k, j)
                     p0r1 = p0_mks(i, k, j)
                     cp409 = c409*pir1
                     cp580 = c580*pir1
                     r32rt = rn32*d2t*sqrt(1./rho(i, k, j))
                     pi0r = pi0(i, k, j)
                     pirr = pir(i, k, j)
                     tairr = (ptwrfr + tb0)*pi0r
                     tair(i, k, j) = tairr
                     ascpr = ascp(i, k, j)
                     avcpr = avcp(i, k, j)
                     rp0r = rp0(i, k, j)
                     if (use_cpm) then
                        cpm = cp*(1.0 - qvwrfr - qlwrfr - qiwrfr - qrwrfr - qswrfr - qgwrfr) &
                              + cvap*qvwrfr &
                              + cliq*(qlwrfr + qrwrfr) &
                              + cice*(qiwrfr + qswrfr + qgwrfr)
                        hlv = alv - (cliq - cvap)*(tairr - t0)
                        hlf = alf - (cice - cliq)*(tairr - t0)
                        hls = hlv + hlf
                        avcpr = hlv/cpm*pirr
                        ascpr = hls/cpm*pirr
                        afcp(i, k, j) = hlf/cpm*pirr
                        ascp(i, k, j) = ascpr
                        avcp(i, k, j) = avcpr
                     end if

                     if (tairr .lt. t0) then
                        if (qiwrfr .le. cmin) qiwrfr = 0.
                        tairc = tairr - t0
                        rtair = 1./(tairr - c76)
                        y2 = exp(c218 - c580*rtair)
                        qsi = rp0r*y2
                        esi = c610*y2
                        if (new_saturation) then
                           esi = min(0.99*p0r1, esi_mks_gpu(tairr))
                           qsi = 0.622*esi/(p0r1 - esi)
                           esi = esi*10.  !in CGS
                        end if
                        ssi = (qvwrfr + qb0)/qsi - 1.
                        y1 = 1./tairr
                        y3 = SQRT(qiwrfr)
                        dd = y1*(RN10A*y1 - RN10B) + RN10C*tairr/ &
                                   esi
                        dm = max(qvwrfr + qb0 - qsi, 0.0)
                        rsub1 = cs580*qsi*rtair*rtair
                        dep = dm/(1.+rsub1)
                        if (tairc .le. -5.) then
                           y4 = 1./(tairr - c358)
                           qsw = rp0r*exp(c172 - c409*y4)
                           fssi = min(xssi, max(.0, xssi*(tairc + 44.)/(44.-38.)))
                           fssi = min(ssi, fssi)
      !                 r_nci=min(1.e-3*exp(-.639+12.96*fssi),1.)
                           r_nci = max(1.e-3*exp(-.639 + 12.96*fssi), 0.528e-3)  ! Meyers et al. 1992 (cm^-3)
                           if (r_nci .gt. 15.) r_nci = 15.

                           pidep = max(R32RT*1.e4*fssi*sqrt(r_nci)*y3/dd, 0.)
                           dd = max(1.e-9*r_nci/r00(i, k, j) &
                                    - qiwrf_old(i, k, j)*1.e-9/ami50, 0.)
                           pint = max(min(dd, dm), 0.)
                           pint = min(pint + pidep, dep)
                           if (pint .le. cmin) pint = 0.
                           ptwrfr = ptwrfr + ascpr*pint
                           qvwrfr = qvwrfr - pint
                           qiwrfr = qiwrfr + pint
                        end if  !taric
                     end if  !tair

      ! End of Process 31 & 32

      ! WRF satice has new_ice_sat option 0, 1 and 2
      ! Steve's satice has new_ice_sat option 0, 1, 2, 3 and 9
      ! option 0, 1 and 2 are identical in both satice
      ! I added option 3 below and wrapped them with "if (improve.eq.3)"

      !*****   TAO ET AL (1989) SATURATION TECHNIQUE  ***********************
                     if (new_ice_sat .eq. 0) then

                        tairr = (ptwrfr + tb0)*pi0r
                        tair(i, k, j) = tairr
                        if (use_cpm) then
                           cpm = cp*(1.0 - qvwrfr - qlwrfr - qiwrfr - qrwrfr - qswrfr - qgwrfr) &
                                 + cvap*qvwrfr &
                                 + cliq*(qlwrfr + qrwrfr) &
                                 + cice*(qiwrfr + qswrfr + qgwrfr)
                           hlv = alv - (cliq - cvap)*(tairr - t0)
                           hlf = alf - (cice - cliq)*(tairr - t0)
                           hls = hlv + hlf
                           avcpr = hlv/cpm*pirr
                           ascpr = hls/cpm*pirr
                           afcp(i, k, j) = hlf/cpm*pirr
                           ascp(i, k, j) = ascpr
                           avcp(i, k, j) = avcpr
                        end if
                        cnd = rt0*(tairr - t00)
                        dep = rt0*(t0 - tairr)
                        y1 = 1./(tairr - c358)
                        y2 = 1./(tairr - c76)
                        qsw = rp0r*exp(c172 - c409*y1)
                        qsi = rp0r*exp(c218 - c580*y2)

                        if (new_saturation) then
                           esw = min(0.99*p0r1, esw_mks_gpu(tairr))
                           esi = min(0.99*p0r1, esi_mks_gpu(tairr))
                           if (esi .gt. esw) esi = esw
                           qsw = 0.622*esw/(p0r1 - esw)
                           qsi = 0.622*esi/(p0r1 - esi)
                           esw = esw*10.  !in CGS
                           esi = esi*10.  !in CGS
                        end if
                        dd = cp409*y1*y1
                        dd1 = cp580*y2*y2
                        if (qlwrfr .le. cmin) qlwrfr = cmin
                        if (qiwrfr .le. cmin) qiwrfr = cmin
                        if (tairr .ge. t0) then
                           dep = 0.0
                           cnd = 1.
                           qiwrfr = 0.0
                        end if

                        if (tairr .lt. t00) then
                           cnd = 0.0
                           dep = 1.
                           qlwrfr = 0.0
                        end if

                        y5 = avcpr*cnd + ascpr*dep
                        y1 = qlwrfr*qsw/(qlwrfr + qiwrfr)
                        y2 = qiwrfr*qsi/(qlwrfr + qiwrfr)
                        y4 = dd*y1 + dd1*y2
                        qvs = y1 + y2
                        rsub1 = (qvwrfr + qb0 - qvs)/(1.+y4*y5)
                        cnd = cnd*rsub1
                        dep = dep*rsub1
                        if (qlwrfr .le. cmin) qlwrfr = 0.
                        if (qiwrfr .le. cmin) qiwrfr = 0.

                        cnd = max(-qlwrfr, cnd)
                        dep = max(-qiwrfr, dep)

                        ptwrfr = ptwrfr + avcpr*cnd &
                                         + ascpr*dep
                        qvwrfr = qvwrfr - cnd - dep
                        qlwrfr = qlwrfr + cnd
                        qiwrfr = qiwrfr + dep

                     end if  !if (new_ice_sat .eq. 0)

                     if (new_ice_sat .eq. 1) then

                        tairr = (ptwrfr + tb0)*pi0r
                        tair(i, k, j) = tairr
                        if (use_cpm) then
                           cpm = cp*(1.0 - qvwrfr - qlwrfr - qiwrfr - qrwrfr - qswrfr - qgwrfr) &
                                 + cvap*qvwrfr &
                                 + cliq*(qlwrfr + qrwrfr) &
                                 + cice*(qiwrfr + qswrfr + qgwrfr)
                           hlv = alv - (cliq - cvap)*(tairr - t0)
                           hlf = alf - (cice - cliq)*(tairr - t0)
                           hls = hlv + hlf
                           avcpr = hlv/cpm*pirr
                           ascpr = hls/cpm*pirr
                           afcp(i, k, j) = hlf/cpm*pirr
                           ascp(i, k, j) = ascpr
                           avcp(i, k, j) = avcpr
                        end if
                        cnd = rt0*(tairr - t00)
                        dep = rt0*(t0 - tairr)
                        y1 = 1./(tairr - c358)
                        y2 = 1./(tairr - c76)
                        qsw = rp0r*exp(c172 - c409*y1)
                        qsi = rp0r*exp(c218 - c580*y2)

                        if (new_saturation) then
                           esw = min(0.99*p0r1, esw_mks_gpu(tairr))
                           esi = min(0.99*p0r1, esi_mks_gpu(tairr))
                           if (esi .gt. esw) esi = esw
                           qsw = 0.622*esw/(p0r1 - esw)
                           qsi = 0.622*esi/(p0r1 - esi)
                           esw = esw*10.  !in CGS
                           esi = esi*10.  !in CGS
                        end if
                        dd = cp409*y1*y1
                        dd1 = cp580*y2*y2
                        y5 = avcpr*cnd + ascpr*dep
                        y1 = rt0*(tairr - t00)*qsw
                        y2 = rt0*(t0 - tairr)*qsi

                        if (tairr .ge. t0) then
                           dep = 0.0
                           cnd = 1.
                           y2 = 0.
                           y1 = qsw
                        end if
                        if (tairr .lt. t00) then
                           cnd = 0.0
                           dep = 1.
                           y2 = qsi
                           y1 = 0.
                        end if
                        y4 = dd*y1 + dd1*y2
                        qvs = y1 + y2
                        rsub1 = (qvwrfr + qb0 - qvs)/(1.+y4*y5)
                        cnd = cnd*rsub1
                        dep = dep*rsub1

                        cnd = max(-qlwrfr, cnd)
                        dep = max(-qiwrfr, dep)

                        ptwrfr = ptwrfr + avcpr*cnd &
                                         + ascpr*dep
                        qvwrfr = qvwrfr - cnd - dep
                        qlwrfr = qlwrfr + cnd
                        qiwrfr = qiwrfr + dep

                     end if  ! if (new_ice_sat .eq. 1)

                     if (new_ice_sat .eq. 2) then

                        dep = 0.0
                        cnd = 0.0
                        tairr = (ptwrfr + tb0)*pi0r
                        tair(i, k, j) = tairr
                        if (use_cpm) then
                           cpm = cp*(1.0 - qvwrfr - qlwrfr - qiwrfr - qrwrfr - qswrfr - qgwrfr) &
                                 + cvap*qvwrfr &
                                 + cliq*(qlwrfr + qrwrfr) &
                                 + cice*(qiwrfr + qswrfr + qgwrfr)
                           hlv = alv - (cliq - cvap)*(tairr - t0)
                           hlf = alf - (cice - cliq)*(tairr - t0)
                           hls = hlv + hlf
                           avcpr = hlv/cpm*pirr
                           ascpr = hls/cpm*pirr
                           afcp(i, k, j) = hlf/cpm*pirr
                           ascp(i, k, j) = ascpr
                           avcp(i, k, j) = avcpr
                        end if

                        if (tairr .ge. 253.16) then
                           y1 = 1./(tairr - c358)
                           qsw = rp0r*exp(c172 - c409*y1)
                           if (new_saturation) then
                              esw = min(0.99*p0r1, esw_mks_gpu(tairr))
                              qsw = 0.622*esw/(p0r1 - esw)
                              esw = esw*10.  !in CGS
                           end if
                           dd = cp409*y1*y1
                           dm = qvwrfr + qb0 - qsw
                           cnd = dm/(1.+avcpr*dd*qsw)

                           cnd = max(-qlwrfr, cnd)

                           ptwrfr = ptwrfr + avcpr*cnd
                           qvwrfr = qvwrfr - cnd
                           qlwrfr = qlwrfr + cnd
                        end if

                        if (tairr .le. 258.16) then
                           y2 = 1./(tairr - c76)
                           qsi = rp0r*exp(c218 - c580*y2)
                           if (new_saturation) then
                              esi = min(0.99*p0r1, esi_mks_gpu(tairr))
                              qsi = 0.622*esi/(p0r1 - esi)
                              esi = esi*10.  !in CGS
                           end if
                           dd1 = cp580*y2*y2
                           dep = (qvwrfr + qb0 - qsi) &
                                 /(1.+ascpr*dd1*qsi)

                           dep = max(-qiwrfr, dep)

                           ptwrfr = ptwrfr + ascpr*dep
                           qvwrfr = qvwrfr - dep
                           qiwrfr = qiwrfr + dep
                        end if
                     end if                              ! if (new_ice_sat .eq. 2)

                     if (new_ice_sat .eq. 3) then

                        dep = 0.0
                        cnd = 0.0
                        tairr = (ptwrfr + tb0)*pi0r
                        tair(i, k, j) = tairr
                        if (use_cpm) then
                           cpm = cp*(1.0 - qvwrfr - qlwrfr - qiwrfr - qrwrfr - qswrfr - qgwrfr) &
                                 + cvap*qvwrfr &
                                 + cliq*(qlwrfr + qrwrfr) &
                                 + cice*(qiwrfr + qswrfr + qgwrfr)
                           hlv = alv - (cliq - cvap)*(tairr - t0)
                           hlf = alf - (cice - cliq)*(tairr - t0)
                           hls = hlv + hlf
                           avcpr = hlv/cpm*pirr
                           ascpr = hls/cpm*pirr
                           afcp(i, k, j) = hlf/cpm*pirr
                           ascp(i, k, j) = ascpr
                           avcp(i, k, j) = avcpr
                        end if

                        if (tairr .ge. t00) THEN
                           y1 = 1./(tairr - c358)
                           qsw = rp0r*exp(c172 - c409*y1)
                           if (new_saturation) then
                              esw = min(0.99*p0r1, esw_mks_gpu(tairr))
                              qsw = 0.622*esw/(p0r1 - esw)
                              esw = esw*10.  !in CGS
                           end if
                           dd = cp409*y1*y1
                           dm = qvwrfr + qb0 - qsw
                           cnd = dm/(1.+avcpr*dd*qsw)

                           cnd = max(-qlwrfr, cnd)

                           ptwrfr = ptwrfr + avcpr*cnd
                           qvwrfr = qvwrfr - cnd
                           qlwrfr = qlwrfr + cnd
                        end if

                        tairr = (ptwrfr + tb0)*pi0r
                        tair(i, k, j) = tairr
                        if (use_cpm) then
                           cpm = cp*(1.0 - qvwrfr - qlwrfr - qiwrfr - qrwrfr - qswrfr - qgwrfr) &
                                 + cvap*qvwrfr &
                                 + cliq*(qlwrfr + qrwrfr) &
                                 + cice*(qiwrfr + qswrfr + qgwrfr)
                           hlv = alv - (cliq - cvap)*(tairr - t0)
                           hlf = alf - (cice - cliq)*(tairr - t0)
                           hls = hlv + hlf
                           avcpr = hlv/cpm*pirr
                           ascpr = hls/cpm*pirr
                           afcp(i, k, j) = hlf/cpm*pirr
                           ascp(i, k, j) = ascpr
                           avcp(i, k, j) = avcpr
                        end if
                        if (tairr .le. 273.16) THEN
                           y1 = 1./(tairr - c358)
                           qsw = rp0r*exp(c172 - c409*y1)
                           y2 = 1./(tairr - c76)
                           qsi = rp0r*exp(c218 - c580*y2)
                           if (new_saturation) then
                              esw = min(0.99*p0r1, esw_mks_gpu(tairr))
                              esi = min(0.99*p0r1, esi_mks_gpu(tairr))
                              if (esi .gt. esw) esi = esw
                              qsw = 0.622*esw/(p0r1 - esw)
                              qsi = 0.622*esi/(p0r1 - esi)
                              esw = esw*10.  !in CGS
                              esi = esi*10.  !in CGS
                           end if
      !                 fssi=min(0.20,max(0.,0.20*(tair(i,j)-t0+44.0)/(44.0-38.0)))
                           fssi = min(xssi, max(0., xssi*(tairr - t0 + 44.0) &
                                 /(44.0 - 38.0)))

                           y3 = 1.+min((qsw - qsi)/qsi, fssi)
                           y4 = qsi*y3
                           if (tairr .le. 268.16 &
                               .and. (qvwrfr + qb0 .gt. y4)) then
                              dd1 = cp580*y2*y2
                              dep = (qvwrfr + qb0 - y4) &
                                    /(1.+ascpr*dd1*y4)
                           elseif (qvwrfr + qb0 .lt. qsi &
                                   .and. qiwrfr .gt. cmin) then
                              dd1 = cp580*y2*y2
                              dep = (qvwrfr + qb0 - qsi) &
                                    /(1.+ascpr*dd1*qsi)
                              dep = max(-qiwrfr, dep)
                           end if

                           ptwrfr = ptwrfr + ascpr*dep
                           qvwrfr = qvwrfr - dep
                           qiwrfr = qiwrfr + dep
                        end if

                     end if                                        ! if (new_ice_sat = 3)
                     ptwrf(i, k, j) = ptwrfr
                     qvwrf(i, k, j) = qvwrfr
                     qlwrf(i, k, j) = qlwrfr
                     qiwrf(i, k, j) = qiwrfr
                  end do 
               end do
            end do
         end if!sat_predict
      end if

   !* 10 * PSDEP : DEPOSITION OR SUBLIMATION OF QS                   **10**
   !* 20 * PGSUB : SUBLIMATION OF QG                                 **20**
      IF (IWARM .ne. 1) THEN
         !$acc parallel loop gang collapse(2) async(async_id) &
         !$acc&         private(cpm, hlv, hlf, hls, psdep, pgdep, pssub, pgsub, &
         !$acc&         rr0r, rrs, fvs, r101r, r102rf, r191r, r192rf, tairc, &
         !$acc&         dlt1, rtair, y2, qsi, esi, ssi, dm, rsub1, dd1, y3, &
         !$acc&         dd, ftns, ftng, ftns0r, ftng0r, y4, y5, y1, tairr, zsr, &
         !$acc&         zgr, ascpr, pirr, r00r, ptwrfr, qvwrfr, qlwrfr, qiwrfr, &
         !$acc&         qrwrfr, qswrfr, qgwrfr, p0r1)
         do j = jts, jte
            do k = kts, kte
               !$acc loop vector
               do i = its, myim(j)
                  ptwrfr = ptwrf(i, k, j)
                  qvwrfr = qvwrf(i, k, j)
                  qlwrfr = qlwrf(i, k, j)
                  qiwrfr = qiwrf(i, k, j)
                  qrwrfr = qrwrf(i, k, j)
                  qswrfr = qswrf(i, k, j)
                  qgwrfr = qgwrf(i, k, j)
                  p0r1 = p0_mks(i, k, j)
                  psdep = 0.0
                  pgdep = 0.0
                  pssub = 0.0
                  pgsub = 0.0
                  zsr = zs(i, k, j)
                  zgr = zg(i, k, j)
                  ascpr = ascp(i, k, j)
                  pirr = pir(i, k, j)
                  r00r = r00(i, k, j)
                  if (improve .eq. 3) then
                     rr0r = 1./rho(i, k, j)
                     rrs = sqrt(rr0r)
                     fvs = sqrt(fv(i, k, j))
                     r101r = rn101*rr0r
                     r102rf = rn102*rrs*fvs
                     r191r = rn191*rr0r
                     r192rf = rn192*rrs*fvs
                  end if
                  tairr = (ptwrfr + tb0)*pi0(i, k, j)
                  tair(i, k, j) = tairr
                  tairc = tairr - t0
   !           if (improve.gt.2) call sgmap(1,qs(i,j),r00,tairc(i,j),ftns0(i,j))
   !           if (improve.gt.2) call sgmap(2,qg(i,j),r00,tairc(i,j),ftng0(i,j))

                  if (qswrfr .lt. cmin1) qswrfr = 0.0
                  if (qgwrfr .lt. cmin1) qgwrfr = 0.0
                  if (qlwrfr + qiwrfr .gt. 1.e-5) then
                     dlt1 = 1.
                  else
                     dlt1 = 0.
                  end if

                  if (tairr .lt. t0) then
                     if (use_cpm) then
                        cpm = cp*(1.0 - qvwrfr - qlwrfr - qiwrfr - qrwrfr - qswrfr - qgwrfr) &
                              + cvap*qvwrfr &
                              + cliq*(qlwrfr + qrwrfr) &
                              + cice*(qiwrfr + qswrfr + qgwrfr)
                        hlv = alv - (cliq - cvap)*(tairr - t0)
                        hlf = alf - (cice - cliq)*(tairr - t0)
                        hls = hlv + hlf
                        avcp(i, k, j) = hlv/cpm*pirr
                        ascpr = hls/cpm*pirr
                        afcp(i, k, j) = hlf/cpm*pirr
                        ascp(i, k, j) = ascpr
                     end if
                     rtair = 1./(tairr - c76)
                     y2 = exp(c218 - c580*rtair)
                     qsi = rp0(i, k, j)*y2
                     esi = c610*y2
                     if (new_saturation) then
                        esi = min(0.99*p0r1, esi_mks_gpu(tairr))
                        qsi = 0.622*esi/(p0r1 - esi)
                        esi = esi*10.  !in CGS
                     end if

                     SSI = (qvwrfr + QB0)/QSI - 1.
                     IF (DLT1 .EQ. 1.) SSI = max(SSI, 0.)
                     IF (DLT1 .EQ. 0.) SSI = min(SSI, 0.)
                     DM = qvwrfr + QB0 - QSI
                     RSUB1 = CS580*QSI*RTAIR*RTAIR
                     DD1 = DM/(1.+RSUB1)
                     Y3 = 1./TAIRr
                     DD = Y3*(RN10A*Y3 - RN10B) + RN10C*TAIRr/ESI
                     TAIRC = TAIRr - T0

                     ftns = 1.
                     ftng = 1.
                     ftns0r = 1.
                     ftng0r = 1.
                     call sgmap_gpu(1, qswrfr, r00r, tairc, &
                                    ftns0r, roqs, roqg, tns, tng)
                     call sgmap_gpu(2, qgwrfr, r00r, tairc, &
                                    ftng0r, roqs, roqg, tns, tng)
                     ftns = ftns0r
                     ftng = ftng0r

                     Y4 = R10T*SSI*(R101R/ZSr**2 + R102RF/ZSr**BSH5) &
                          /DD*ftns
                     PSDEP = max(-qswrfr, Y4)
                     DD = Y3*(RN20A*Y3 - RN20B) + RN10C*TAIRr/ESI
                     Y2 = R191R/ZGr**2 + R192RF/ZGr**BGH5
                     PGDEP = MAX(-qgwrfr, R20T*SSI*Y2/DD &
                                       *ftng)

                     Y5 = min(0., DD1)
                     DD1 = max(0., DD1)
                     IF (DLT1 .EQ. 1.) THEN
                        Y1 = MIN(PSDEP + PGDEP, DD1)
                        IF (PSDEP .ge. DD1) THEN
                           PSDEP = DD1
                           PGDEP = 0.
                        END IF
                        IF (DD1 .gt. PSDEP .and. (PSDEP + PGDEP) .gt. &
                            DD1) PGDEP = Y1 - PSDEP
                     END IF
                     IF (DLT1 .EQ. 0.) THEN
                        Y1 = MAX(PSDEP + PGDEP, Y5)
                        IF (Y5 .gt. (PSDEP + PGDEP)) THEN
                           Y3 = (PSDEP + PGDEP)
                           IF (Y3 .ne. 0.0) THEN
                              PSDEP = PSDEP/Y3*Y5
                              PGDEP = PGDEP/Y3*Y5
                           END IF
                        END IF
                     END IF

                     PSSUB = r2is*min(PSDEP, 0.)
                     PSDEP = r2is*max(PSDEP, 0.)
                     PGSUB = r2ig*min(PGDEP, 0.)
                     PGDEP = r2ig*max(PGDEP, 0.)

                     ptwrfr = ptwrfr + ascpr*y1
                     qvwrfr = qvwrfr - y1
                     qswrfr = qswrfr + psdep + pssub
                     qgwrfr = qgwrfr + pgdep + pgsub

                  end if   ! if (tair(i,j) .lt. t0)
                  ptwrf(i, k, j) = ptwrfr
                  qvwrf(i, k, j) = qvwrfr
                  qswrf(i, k, j) = qswrfr
                  qgwrf(i, k, j) = qgwrfr
               end do
            end do
         end do
      end if

   !!!  end of Processes 10 and 20
   !>>> Note that ern is already calculated in saturation prediction scheme.

   !* 23 * ERN : EVAPORATION OF QR (SUBSATURATION)                   **23**
      IF (IWARM .ne. 1) THEN
         !$acc parallel loop gang collapse(2) async(async_id) private(cpm, hlv, &
         !$acc&         hlf, hls, rr0r, fv0r, rrs, fvs, r191r, r192rf, r331r, &
         !$acc&         r332rf, r231r, r232rf, rtair, y2, esw, qsw, ssw, dm, &
         !$acc&         rsub1, dd1, y3, dd, y1, ern, pmlts, pmltg, tairc, &
         !$acc&         ftns0r, ftng0r, ftns, ftng, tairr, zrr, zsr, zgr, &
         !$acc&         ascpr, avcpr, pi0r, pirr, rp0r, r00r, ptwrfr, qvwrfr, &
         !$acc&         qlwrfr, qiwrfr, qrwrfr, qswrfr, qgwrfr, p0r1)
         do j = jts, jte
            do k = kts, kte
               !$acc loop vector
               do i = its, myim(j)
                  ptwrfr = ptwrf(i, k, j)
                  qvwrfr = qvwrf(i, k, j)
                  qlwrfr = qlwrf(i, k, j)
                  qiwrfr = qiwrf(i, k, j)
                  qrwrfr = qrwrf(i, k, j)
                  qswrfr = qswrf(i, k, j)
                  qgwrfr = qgwrf(i, k, j)
                  zrr = zr(i, k, j)
                  zsr = zs(i, k, j)
                  zgr = zg(i, k, j)
                  ascpr = ascp(i, k, j)
                  avcpr = avcp(i, k, j)
                  pi0r = pi0(i, k, j)
                  pirr = pir(i, k, j)
                  r00r = r00(i, k, j)
                  p0r1 = p0_mks(i, k, j)
                  if (improve .eq. 3) then
                     rp0r = rp0(i, k, j)
                     rr0r = 1./rho(i, k, j)
                     fv0r = fv(i, k, j)
                     rrs = sqrt(rr0r)
                     fvs = sqrt(fv0r)
                     r191r = rn191*rr0r
                     r192rf = rn192*rrs*fvs
                     r331r = rn331*rr0r
                     r332rf = rn332*rrs*fvs
                     r231r = rn231*rr0r
                     r232rf = rn232*rrs*fvs
                  end if
                  if (sat_predict .eq. .false.) then
                     if (qrwrfr .gt. 0.0) then
                        tairr = (ptwrfr + tb0)*pi0r
                        tair(i, k, j) = tairr
                        rtair = 1./(tairr - c358)
                        if (use_cpm) then
                           cpm = cp*(1.0 - qvwrfr - qlwrfr - qiwrfr - qrwrfr - qswrfr - qgwrfr) &
                                 + cvap*qvwrfr &
                                 + cliq*(qlwrfr + qrwrfr) &
                                 + cice*(qiwrfr + qswrfr + qgwrfr)
                           hlv = alv - (cliq - cvap)*(tairr - t0)
                           hlf = alf - (cice - cliq)*(tairr - t0)
                           hls = hlv + hlf
                           avcpr = hlv/cpm*pirr
                           ascpr = hls/cpm*pirr
                           afcp(i, k, j) = hlf/cpm*pirr
                           ascp(i, k, j) = ascpr
                           avcp(i, k, j) = avcpr
                        end if
                        y2 = exp(c172 - c409*rtair)
                        esw = c610*y2
                        qsw = rp0r*y2
                        if (new_saturation) then
                           esw = min(0.99*p0r1, esw_mks_gpu(tairr))
                           qsw = 0.622*esw/(p0r1 - esw)
                           esw = esw*10.  !in CGS
                        end if
                        ssw = (qvwrfr + qb0)/qsw - 1.
                        dm = qvwrfr + qb0 - qsw
                        rsub1 = cv409*qsw*rtair*rtair
                        dd1 = max(-dm/(1.+rsub1), 0.0)
                        y3 = 1./tair(i, k, j)
                        dd = y3*(rn30a*y3 - rn10b) + rn10c*tairr &
                                   /esw
                        y1 = -r23t*ssw*(r231r/zrr**2 + r232rf/ &
                                                    zrr**3)/dd
                        ern = min(dd1, qrwrfr, max(y1, 0.0))

                        ptwrfr = ptwrfr - avcpr*ern
                        qvwrfr = qvwrfr + ern
                        qrwrfr = qrwrfr - ern
                     end if
                  end if !sat_predict

   !* 30 * pmltg : evaporation of melting qg                         **30**
   !* 33 * pmlts : evaporation of melting qs                         **33**
                  pmlts = 0.0
                  pmltg = 0.0
                  tairr = (ptwrfr + tb0)*pi0r
                  tair(i, k, j) = tairr
                  tairc = tairr - t0
                  if (use_cpm) then
                     cpm = cp*(1.0 - qvwrfr - qlwrfr - qiwrfr - qrwrfr - qswrfr - qgwrfr) &
                           + cvap*qvwrfr &
                           + cliq*(qlwrfr + qrwrfr) &
                           + cice*(qiwrfr + qswrfr + qgwrfr)
                     hlv = alv - (cliq - cvap)*(tairr - t0)
                     hlf = alf - (cice - cliq)*(tairr - t0)
                     hls = hlv + hlf
                     avcpr = hlv/cpm*pirr
                     ascpr = hls/cpm*pirr
                     afcp(i, k, j) = hlf/cpm*pirr
                     ascp(i, k, j) = ascpr
                     avcp(i, k, j) = avcpr
                  end if

                  ftns0r = 1.
                  ftng0r = 1.

                  if (tairr .ge. t0) then
                     ftns = 1.
                     ftng = 1.
                     call sgmap_gpu(1, qswrfr, r00r, tairc, &
                                    ftns0r, roqs, roqg, tns, tng)
                     call sgmap_gpu(2, qgwrfr, r00r, tairc, &
                                    ftng0r, roqs, roqg, tns, tng)
                     ftns = ftns0r
                     ftng = ftng0r

   !              rtair(i,j)=1./(tair(i,j)-c358)
                     rtair = 1./(t0 - c358)
                     y2 = exp(c172 - c409*rtair)
                     esw = c610*y2
                     qsw = rp0r*y2
                     if (new_saturation) then
                     esw = min(0.99*p0r1, esw_mks_gpu(tairr))
                     qsw = 0.622*esw/(p0r1 - esw)
                     esw = esw*10.  !in CGS
                     end if
                     
                     ssw = 1.-(qvwrfr + qb0)/qsw
                     dm = qsw - qvwrfr - qb0
                     rsub1 = cv409*qsw*rtair*rtair
                     dd1 = max(dm/(1.+rsub1), 0.0)
                     y3 = 1./tairr
                     dd = y3*(rn30a*y3 - rn10b) + rn10c*tairr/esw
                     y1 = ftng*r30t*ssw*(r191r/zgr**2 + r192rf &
                          /zgr**bgh5)/dd
                     pmltg = min(qgwrfr, max(y1, 0.0))
                     y1 = ftns*r33t*ssw*(r331r/zsr**2 + r332rf &
                          /zsr**bsh5)/dd
                     pmlts = min(qswrfr, max(y1, 0.0))
                     y1 = min(pmltg + pmlts, dd1)
                     pmltg = y1 - pmlts
                     ptwrfr = ptwrfr - ascpr*y1
                     qvwrfr = qvwrfr + y1
                     qswrfr = qswrfr - pmlts
                     qgwrfr = qgwrfr - pmltg
                  end if
   ! end   Processes 30 and 33
                  ptwrf(i, k, j) = ptwrfr
                  qvwrf(i, k, j) = qvwrfr
                  qrwrf(i, k, j) = qrwrfr
                  qswrf(i, k, j) = qswrfr
                  qgwrf(i, k, j) = qgwrfr
               end do
            end do
         end do
      END IF    ! part of if (iwarm.eq.1) then

   #ifdef EXT_DIAG
      do j = jts, jte
         do k = kts, kte
            do i = its, myim(j)
               scc = 0.
               see = 0.

               dd = max(-cnd, 0.)
               cnd = max(cnd, 0.)
   !        dd1(i,j)=max(-dep(i,j), 0.)+pidep(i,j)*d2t
               dd1 = max(-dep, 0.)  !bug fix by Di
               dep = max(dep, 0.)

               sccc = cnd
               seee = dd + ern
               sddd = dep + dmax1(pint, 0.0) + psdep + &
                      pgdep
               ssss = dd1 - dmin1(pint, 0.0) + pssub + &
                      pgsub + pmlts + pmltg
               smmm = psmlt + pgmlt + pimlt + qracs
               sfff = psacw*d2t + piacr*d2t + psfw*d2t + pgfr*d2t &
                      + dgacw*d2t + dgacr*d2t + psacr*d2t + pihom &
                      + pidw + pimm + pcfr + pihms*d2t + pihmg*d2t

               ! Snapshot values (K/s), consistent with declared units in Registry
               physc(i, k, j) = avcp(i, k, j)*sccc/d2t
               physe(i, k, j) = avcp(i, k, j)*seee/d2t
               physd(i, k, j) = ascp(i, k, j)*sddd/d2t
               physs(i, k, j) = ascp(i, k, j)*ssss/d2t
               physf(i, k, j) = afcp(i, k, j)*sfff/d2t
               physm(i, k, j) = afcp(i, k, j)*smmm/d2t
               ! Accumulated values (K)
               acphysc(i, k, j) = acphysc(i, k, j) + avcp(i, k, j)*sccc
               acphyse(i, k, j) = acphyse(i, k, j) + avcp(i, k, j)*seee
               acphysd(i, k, j) = acphysd(i, k, j) + ascp(i, k, j)*sddd
               acphyss(i, k, j) = acphyss(i, k, j) + ascp(i, k, j)*ssss
               acphysf(i, k, j) = acphysf(i, k, j) + afcp(i, k, j)*sfff
               acphysm(i, k, j) = acphysm(i, k, j) + afcp(i, k, j)*smmm

   !        radar reflectivity calculation :
   !        dbz(i,k,j)=0.0
               a_1 = 1.e6*r00(i, k, j)*qrwrf(i, k, j)
               a_2 = 1.e6*r00(i, k, j)*qswrf(i, k, j)
               a_3 = 1.e6*r00(i, k, j)*qgwrf(i, k, j)
   !        if (a_1+a_2+a_3 .ge. 1.e-8)  then
               tairr = (ptwrf(i, k, j) + tb0)*pi0(i, k, j)
               tair(i, k, j) = tairr
               tairc = tairr - t0
               uwet = 4.464**0.95

               ftns = 1.
               ftns0r = 1.
               ftng = 1.
               ftng0r = 1.
               call sgmap_gpu(1, qswrf(i, k, j), r00(i, k, j), tairc, &
                              ftns0r, roqs, roqg, tns, tng)
               call sgmap_gpu(2, qgwrf(i, k, j), r00(i, k, j), tairc, &
                              ftng0r, roqs, roqg, tns, tng)
               A_C = 49.6*1.E-6
               A_I = 9.4*1.E-6
               RE_C = 10.            ! UM
               RE_I = 25.            ! UM
               RE_S = 75.            ! UM
   !        tairc5=min(0.,max(-35.,tairc(i,j)+5.0))
   !        fexp=-0.2-tairc5/100.
               ucor = 3071.29/tnw**0.75

   ! Caculate temperature dependent ucog
   !        ftng(i,j)=1.
               ftng = ftng0r**0.75
   !        ftng(i,j)=exp(-1.*tslopeg*tairc5)
   !        ftngQ=1.0
   !        if (qg(i,j).gt.cmin) ftngQ=(qg(i,j)/0.001)**fexp
   !        ftng(i,j)=ftng(i,j)*ftngQ
   ! Caculate temperature dependent ucos
   !        ftns(i,j)=1.
               ftns = ftns0r**0.75
   !        ftns(i,j)=exp(-1.*tslopes*tairc5)
   !        ftnsQ=1.0
   !        if (qs(i,j).gt.cmin) ftnsQ=(qs(i,j)/0.001)**fexp
   !        ftns(i,j)=ftns(i,j)*ftnsQ
   !
   !         ucos=687.97*roqs**0.25/ftns(i,j)**0.75
   !         ucog=687.97*roqg**0.25/ftng(i,j)**0.75
               ucos = 687.97*roqs**0.25/tns**0.75/ftns
               ucog = 687.97*roqg**0.25/tng**0.75/ftng

               a_11 = ucor*(max(1.e-12, a_1))**1.75
   !        a_22=ucos*(max(1.e-12,a_2))**1.75
   !        a_33=ucog*(max(1.e-12,a_3))**1.75

   !        a_22=ucos*(max(1.e-12,a_2))**1.75/ftns(i,j)
   !        a_33=ucog*(max(1.e-12,a_3))**1.75/ftng(i,j)
               a_22 = ucos*(max(1.e-12, a_2))**1.75
               a_33 = ucog*(max(1.e-12, a_3))**1.75

               W_C = r00(i, k, j)*qlwrf(I, k, J)*1000.   !LIQUID WATER CONTENT IN G/M^3
               W_I = r00(i, k, j)*qiwrf(I, k, J)*1000.   !ICE WATER CONTENT IN G/M^3
               W_S = r00(i, k, j)*qswrf(I, k, J)*1000.   !ICE WATER CONTENT IN G/M^3

   ! Xiping's formula
   !        ZE_CLD = A_C * W_C * RE_C**3 + A_I * W_I * RE_I**3   &  !ZE IN  mm^6/m^3
   !                + A_I * W_S * RE_S**3  ! (as cloud particles NOT precipitation)
               ZE_CLD = A_C*W_C*RE_C**3 + A_I*W_I*RE_I**3      !ZE IN  mm^6/m^3
               ! (use a_22 and ucos instead)
               IF (TAIRr .LT. 273.16) THEN
   !           ZDRY = MAX(1.e-4,A_11+A_33+ZE_CLD) ! Xiping's  !rain,snow,cloud ice,cloud water,graupel
                  ZDRY = MAX(1.e-4, A_11 + A_22 + A_33 + ZE_CLD) !rain,snow,cloud ice,cloud water,graupel
   !           DBZ(I,K,J) = 10.*ALOG10(ZDRY)
                  DBZ(I, K, J) = 10.*DLOG10(ZDRY)
               ELSE
   !           A_44 = A_11+UWET*(A_22+A_33)**.95         ! old formula
   !           A_44 = A_11+UWET*A_33**.95+ZE_CLD         ! Xiping's
                  A_44 = A_11 + UWET*(A_22 + A_33)**.95 + ZE_CLD
                  ZWET = MAX(1.e-4, A_44)
   !           DBZ(I,K,J) = 10.*ALOG10(ZWET)
                  DBZ(I, K, J) = 10.*DLOG10(ZWET)
               END IF
            end do
         end do
      end do
   !end of EXT_DIAG
   #endif

   !   endif
      !$acc parallel loop gang collapse(2) async(async_id) private(L_cloud, &
      !$acc&         ccn_ref, mu, gamfac3, gamfac1, lambda, mdc1, mdc2, mdc3, &
      !$acc&         mdc4, mdc5, mdc6, efd1, efd2, efd3, efd4, efd5, efd6, &
      !$acc&         ltk, lqc, ltk2, lqc2, mvdc, efdc, reimin, iwc_0, bw98, &
      !$acc&         lqi, lqi2, efdi, tairc, tairr, rhor, xlandr, qlwrfr, &
      !$acc&         qiwrfr, qrwrfr, qswrfr, qgwrfr, rhoair, refcr, refrr, &
      !$acc&         refir, refsr, refgr)
      do j = jts, jte
         do k = kts, kte
            !$acc loop vector
            do i = its, myim(j)
               tairr = tair(i, k, j)
               tairc = tairr - t0
               rhor = rho(i, k, j)
               xlandr = xland(i, j)
               qlwrfr = qlwrf(i, k, j)
               qiwrfr = qiwrf(i, k, j)
               qrwrfr = qrwrf(i, k, j)
               qswrfr = qswrf(i, k, j)
               qgwrfr = qgwrf(i, k, j)
               rhoair = rho_mks(i, k, j)
               refir = refi(i, k, j)
   !        IF (QV(I,J)+QB0 .LE. 0.) QV(I,J)=-QB0
               if (qlwrfr .le. cmin) qlwrfr = 0.
               if (qrwrfr .le. cmin) qrwrfr = 0.
               if (qiwrfr .le. cmin) qiwrfr = 0.
               if (qswrfr .le. cmin) qswrfr = 0.
               if (qgwrfr .le. cmin) qgwrfr = 0.


   !   eff_rad is a function of the slope parameter (Lambda)
               ! effective radii of rain :
               if (qrwrfr .lt. crmin) then
                  refrr = 0.e0
               else
                  refrr = eff_rad_gpu(zr(i, k, j))
               end if

               ! effective radii of snow
               if (qswrfr .lt. csmin) then
                  refsr = 0.e0
               else
                  refsr = eff_rad_gpu(zs(i, k, j))
               end if

               ! effective radii of graupel/hail
               if (qgwrfr .lt. cgmin) then
                  refgr = 0.e0
               else
                  refgr = eff_rad_gpu(zg(i, k, j))
               end if

               ! effective radii of  cloud water
               if (rewflag .eq. 1) then
                  ! default
                  if (qlwrfr .lt. cmin) then
                     refcr = 0.e0
                  else
                     L_cloud = qlwrfr*rhor             ! cloud water [g/cm3]
                     if (xlandr .eq. 1.0) then
                        ccn_ref = ccn_over_land
                     elseif (xlandr .eq. 2.0) then
                        ccn_ref = ccn_over_water
                     else
                        print *, ' xland is not 1. or 2., run stopped'
                        stop
                     end if
                     ! for cloud water, estimate lambda (slope of gamma distribution)
                     mu = min(15.e0, (1000.E0/ccn_ref + 2.e0))
                     gamfac3 = (gamma_toshi_gpu(mu + 4.e0)/gamma_toshi_gpu(mu + 3.e0))
                     gamfac1 = (gamma_toshi_gpu(mu + 4.e0)/gamma_toshi_gpu(mu + 1.e0))
                     lambda = (4.e0/3.e0*cpi*roqr*ccn_ref/L_cloud* &
                               gamfac1)**(1.e0/3.e0)  ! [1/cm]
                     refcr = 1.e0/lambda*gamfac3*1.e4  !effective radius [micron]
                  end if ! qcl(i,j,k) < cmin test
               end if  !end of if rewflag=1

               if (rewflag .eq. 2) then
                  ! mapping spectrum, different over land and ocean
                  if (qlwrfr .ge. cwmin) then
                     if (xlandr .eq. 1.) then
                        mdc1 = 5.8936819
                        mdc2 = -7.013013
                        mdc3 = 1.3178721
                        mdc4 = 1.1741987
                        mdc5 = 2.6110916E-3
                        mdc6 = -0.26646396
                     else
!                        mdc1 = 173.57305
                        mdc1 = 173.27305   ! decrease diameter
                        mdc2 = -64.370929
                        mdc3 = 0.36833626
                        mdc4 = 6.1389254
                        mdc5 = 5.5915321E-3
                        mdc6 = -0.12488698
                     end if
                     efd1 = 1.6855156
                     efd2 = 0.84147302
                     efd3 = -0.46783369
                     efd4 = 1.005305E-2
                     efd5 = 3.4833332E-2
                     efd6 = 1.4727223E-2

                     ltk = log(tairr)
                     lqc = -1.*log(qlwrfr*rhoair)
                     ltk2 = ltk*ltk
                     lqc2 = lqc*lqc
                     mvdc = exp(mdc1 + mdc2*ltk + mdc3*lqc + mdc4*ltk2 &
                                + mdc5*lqc2 + mdc6*ltk*lqc)
                     efdc = exp(efd1 + efd2*log(mvdc) + efd3*log(1000.) &
                                + efd4*log(mvdc)**2.+efd5*log(1000.)**2. &
                                +efd6*log(mvdc)*log(1000.))
                     refcr = min(max(efdc/2., 0.5), 50.)
                  else
                     refcr = 0.5  !test
                  end if  !end of if qcl>=cwmin
               end if  !end of if rewflag=2

               ! radiative radii of cloud ice :
               if (reiflag .eq. 1) then
                  if (qiwrfr .lt. cimin) then
                     refir = 0.e0
                  else
                     ! for cloud ice effective radius depends on temperature profile, formula from GCE
                     refir = 125.e0 + (tairr - 243.16)*5.e0     ! [micron]
                     if (tairr .gt. 243.16) refir = 125.e0
                     if (tairr .lt. 223.16) refir = 25.e0
                  end if !end of if qci(i,j,k) < cimin
               end if  !end of if reiflag=1

               if (reiflag .eq. 2) then
                  ! Wyser 1998 , equation 15 and 35:
                  reimin = 10.0
                  if (qiwrfr .ge. cimin) then
                     iwc_0 = 50.e-6  ! note that rho here is in CGS
                     ! IWC_0 (50 g/m^-3) must be converted to g/cm^3
                     bw98 = -2.+1.e-3*log10(rhor*qiwrfr/iwc_0) &
                            *max(0.0, -tairc)**1.5
                     refir = 377.4 + bw98*(203.3 + bw98*(37.91 + 2.3696*bw98))
                     refir = max(reimin, refir)
                  else
                     refir = 0.0
                  end if
               end if  !end of if reiflag=2

               if (reiflag .eq. 3) then
                  ! Heymsfield et al. 2014 , equation 9e :
                  reimin = 10.0
                  if (qiwrfr .ge. cimin) then
                     if (tairc >= -56.0 .and. tairc < 0.0) then
                        refir = 308.4*exp(0.0152*tairc)      ! 131.657 ~ 308.4 micron
                     elseif (tairc >= -71.0 .and. tairc < -56.0) then
                        refir = 9.1744e+4*exp(0.117*tairc)   ! 22.641 ~ 130.942 micron
                     elseif (tairc < -71.0) then
                        refir = 83.3*exp(0.0184*tairc)       ! 10 ~ 22.557 micron
                     end if
                     refir = max(reimin, refir)
                  else
                     refir = 0.0
                  end if
               end if  !end of if reiflag=3

               if (reiflag .eq. 4) then
                  reimin = 0.0
                  ! Mitchell et al. 2011, equation 9
                  ! IWC in mg/m^3, T in oC, rei in micron
                  if (qiwrfr .ge. cimin) then
                     refir = 132.23 + 1.2433*tairc + &
                                     8.6629*(6.+log10(rhor*qiwrfr))
                     refir = max(reimin, refir)
                  else
                     refir = 0.0
                  end if
               end if  !end of if reiflag=4

               if (reiflag .eq. 5) then
                  ! mapping spectrum, different over land and ocean
                  if (qiwrfr .ge. cimin) then
                     ltk = log(tairr)
                     lqi = -1.*log(qiwrfr*rhoair)
                     ltk2 = ltk*ltk
                     lqi2 = lqi*lqi
                     if (xlandr .eq. 1.0) then
                        ! over land
!                        efdi = exp(161.47584 - 0.26232591*lqi + 4.3393883E-4*lqi2 &
                        efdi = exp(161.77584 - 0.26232591*lqi + 4.3393883E-4*lqi2 &  ! increase diameter
                                   - 57.057846*ltk + 5.3668153*ltk2)/10.
                     else
                        ! over ocean
!                        efdi = exp(161.92213 - 0.26232591*lqi + 4.3393883E-4*lqi2 &
!                        efdi = exp(161.62213 - 0.26232591*lqi + 4.3393883E-4*lqi2 &   ! decrease diameter
                        efdi = exp(161.82213 - 0.26232591*lqi + 4.3393883E-4*lqi2 &   ! slightly decrease diameter
                                   - 57.057846*ltk + 5.3668153*ltk2)/10.
                     end if
                     refir = efdi/2.
                     refir = min(max(refir, 1.5), 500.)
                  else
                     refir = 1.5  !test
                  end if
               end if  !end of if reiflag=5
               qlwrf(i, k, j) = qlwrfr
               qiwrf(i, k, j) = qiwrfr
               qrwrf(i, k, j) = qrwrfr
               qswrf(i, k, j) = qswrfr
               qgwrf(i, k, j) = qgwrfr
               refc(i, k, j) = refcr
               refr(i, k, j) = refrr
               refi(i, k, j) = refir
               refs(i, k, j) = refsr
               refg(i, k, j) = refgr


!!!!!!!!!!!! end saticel_s_B_loop


! ****************************************************************
! convert from GCE grid back to WRF grid
               w_mks(i, k, j) = w_mks(i, k, j)*0.01
   !         icn_diag(i,k,j) = icn_cgs(i,j,k) * 1000. ! #/Litre <-- #/cm3
   !         nc_diag(i,k,j) = nc_cgs(i,j,k)  ! #/cm3
   !         if (gid.eq.1)  icn_diag(i,k,j) = 0.
            end do !i
         end do !j
      end do !k
      !$acc end data
!       do j=jts,jte
!       do i=its,ite
!          comdbz(i,j) =  dbz(i,kts,j)
!          do k=kts,kte
!             if (dbz(i,k,j) .gt. comdbz(i,j)) comdbz(i,j) = dbz(i,k,j)
!          enddo !k
!       enddo !i
!       enddo !j

! ****************************************************************

#ifdef EXT_DIAG
      IF (PRESENT(diagflag)) THEN
      if (diagflag .and. do_radar_ref == 1) then
         do j = jts, jte
            do k = kts, kte
               do i = its, ite
                  refl_10cm(i, k, j) = max(-35., dBZ(i, k, j))
               end do
            end do
         end do
      end if
      END IF
#endif
!+---+-----------------------------------------------------------------+
      !if (benchmark .eq. .true.) call nvtxEndRange

      return

      END SUBROUTINE saticel_s_gpu

            SUBROUTINE auto_conversion(L, N, P, re)
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
               real, intent(in) :: L    ! cloud liquid water [g cm-3]
               real, intent(in) :: N    ! total number concentration [# cm-3]
               real, intent(out) :: P   ! auto conversion rate [g cm-3 s-1]
               real, intent(out) :: re  ! cloud effective radius [micron]

               real :: mu   ! mu of gamma PSD [-]
               real :: eta  ! eta function [cm3 g-2 s-1]
               real :: betaf, beta1, beta2     ! beta function [-]
               real :: gamfac, gfac1, gfac2 ! gamma function [-]
               real :: R6_6power  ! mean radius of the sixth moment [cm]
               real :: R6_thresh ! threshold of  mean radius of the sixth moment [cm]
               real :: R6        ! mean radius of the sixth moment [cm]
               real :: Heaviside_func  ! Heaviside step function (0 or 1)
               real :: lambda    ! slope of gamma size ditribution [1/cm]
! real :: No        ! intercept  [cm-4]

               real, parameter :: Rc = 10.e0*1.e-4 ! threshold of particle radus (10 micron) [cm]
               real, parameter :: const_pi = 3.14159e0 ! pai
               real, parameter :: const_kappa = 1.9e11    ! coefficient for water droplet collection kernel [cm-3 s-1]
               ! from Long [1974, JAS].
! real,parameter :: const_kappa = 1.9e11*10000.e0 !10000 is to adjust the order to keseller

               real, parameter :: const_rho_liq = 1.e0    ! density of liquid water [g cm-3]
               real, parameter :: eta_func = ((3.e0/(4.e0*const_pi*const_rho_liq))**2)*const_kappa  ! eta function [cm3 g-2 s-1]
               ! a part of (eq 27b)

               logical, parameter :: no_thresh = .true.  ! logic to choose no threshold parameterization or not.

!
! When no particel, no autoconversion.
!
! EMK BUG FIX...Prevent overflow for small but non-zero values of L
! if( N <= 0.e0 .or. L <= 0.e0 ) then
               if (N <= 0.e0 .or. L <= 1.0e-32) then
                  P = 0.e0
                  return
               end if

! orig
               mu = MIN(15.e0, (1000.E0/N + 2.e0))

!
! gamma functions
!
               gfac1 = gamma_toshi_gpu(mu + 4.e0)
               gfac2 = gamma_toshi_gpu(mu + 1.e0)
               gamfac = (gfac1/gfac2)

!
! estimate lambda (slope of gamma distribution)
!
               lambda = (4.e0/3.e0*const_pi*const_rho_liq*N/L*gamfac)**(1.e0/3.e0)  ! [1/cm]

               THRESH: if (no_thresh) then !-------------------------------------------

!
! threshold of particle radius (mean radius of the sixth moment )
!
                  gfac1 = gamma_toshi_gpu(mu + 7.e0)
                  gfac2 = gamma_toshi_gpu(mu + 1.e0)
                  gamfac = (gfac1/gfac2)

                  R6_6power = (1.e0/lambda)**6.e0*gamfac   ![cm] (eq. A3)

!
! auto conversion rate (eq. 26a)
!
                  P = const_kappa*N*R6_6power*L      ! [g cm-3 s-1 ]
!    [cm-3 s-1]  * [#/cm3] *   [cm6]   * [g/cm3]

               else  !with threshold ---------------------------------------------------

!
! Estimate eta under gamma PSD
!
                  beta1 = (6.e0 + mu)*(5.e0 + mu)*(4.e0 + mu)
                  beta2 = (3.e0 + mu)*(2.e0 + mu)*(1.e0 + mu)
                  betaf = beta1/beta2

                  eta = eta_func*betaf  ! eta function (eq 27b) [cm3 g-2 s-1]

!
! threshold of particle radius (mean radius of the sixth moment )
!
                  R6_thresh = betaf*Rc  ![cm] (pg 1545)

!
! mean radius of the sixth moment
!
                  gfac1 = gamma_toshi_gpu(6.e0 + mu + 1.e0)
                  gfac2 = gamma_toshi_gpu(1.e0 + mu)
                  gamfac = (gfac1/gfac2)**(1.e0/6.e0)

                  R6 = (1.e0/lambda)*gamfac   ![cm] (eq. A3)

!
! Heaviside step function
!
                  if (R6 - R6_thresh <= 0.e0) then
                     Heaviside_func = 0.e0
                  elseif (R6 - R6_thresh > 0.e0) then
                     Heaviside_func = 1.e0
                  else
                     print *, 'MSG: auto_conversion: Strange value of R6= ', R6; stop
                  end if

!
! auto conversion rate [g cm-3 s-1 ] (eq. 27a)
!
                  P = eta*(1.e0/N)*(L**3)*Heaviside_func

!    [cm3 g-2 s-1] * [cm3] * [g3/cm9]

               end if THRESH !------------------------------------------------------------

! optional

!
! estimate effective radius
!
               gfac1 = gamma_toshi_gpu(mu + 4.e0)
               gfac2 = gamma_toshi_gpu(mu + 3.e0)
               gamfac = (gfac1/gfac2)

               re = 1.e0/lambda*gamfac*1.e4  !effective radius [micron]

!
! estimate No
!
! call gamma_function(mu+1.e0 ,gfac1)
! No = N * (lambda**(mu+1)) / gfac1

               RETURN
            END subroutine auto_conversion

            real function gamma_toshi_gpu(x)

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
               implicit double precision(a - h, o - z)
               dimension g(26)
               data g/1.0d0, 0.5772156649015329d0, &
                  -0.6558780715202538d0, -0.420026350340952d-1, &
                  0.1665386113822915d0, -.421977345555443d-1, &
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

               pi = 3.141592653589793d0
               if (x .eq. int(x)) then
                  if (x .gt. 0.0d0) then
                     ga = 1.0d0
                     m1 = int(x) - 1
                     do k = 2, m1
                        ga = ga*k
                     end do
                  else
                     ga = 1.0d+300
                  end if
               else
                  if (dabs(dble(x)) .gt. 1.0d0) then
                     z = dabs(dble(x))
                     m = int(z)
                     r = 1.0d0
                     do k = 1, m
                        r = r*(z - k)
                     end do
                     z = z - m
                  else
                     z = dble(x)
                  end if
                  gr = g(26)
                  do k = 25, 1, -1
                     gr = gr*z + g(k)
                  end do
                  ga = 1.0d0/(gr*z)
                  if (dabs(dble(x)) .gt. 1.0d0) then
                     ga = ga*r
                     if (x .lt. 0.0d0) ga = -pi/(x*ga*dsin(pi*x))
                  end if
               end if

               gamma_toshi_gpu = real(ga)

               return
            end function gamma_toshi_gpu

!+---+-----------------------------------------------------------------+

            subroutine refl10cm_gsfc(qv1d, qr1d, qs1d, qg1d, &
                                     t1d, p1d, dBZ, kts, kte, ii, jj, ihail)

               IMPLICIT NONE

!..Sub arguments
               INTEGER, INTENT(IN):: kts, kte, ii, jj, ihail
               REAL, DIMENSION(kts:kte), INTENT(IN):: &
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
               REAL, PARAMETER:: R = 287.
               REAL, PARAMETER:: PIx = 3.1415926536

!+---+

               do k = kts, kte
                  dBZ(k) = -35.0
               end do

!+---+-----------------------------------------------------------------+
!..Put column of data into local arrays.
!+---+-----------------------------------------------------------------+
               do k = kts, kte
                  temp(k) = t1d(k)
                  qv(k) = MAX(1.E-10, qv1d(k))
                  pres(k) = p1d(k)
                  rho(k) = 0.622*pres(k)/(R*temp(k)*(qv(k) + 0.622))

                  if (qr1d(k) .gt. 1.E-9) then
                     rr(k) = qr1d(k)*rho(k)
                     N0_r(k) = xnor
                     lamr = (xam_r*xcrg(3)*N0_r(k)/rr(k))**(1./xcre(1))
                     ilamr(k) = 1./lamr
                     L_qr(k) = .true.
                  else
                     rr(k) = 1.E-12
                     L_qr(k) = .false.
                  end if

                  if (qs1d(k) .gt. 1.E-9) then
                     rs(k) = qs1d(k)*rho(k)
                     N0_s(k) = xnos
                     lams = (xam_s*xcsg(3)*N0_s(k)/rs(k))**(1./xcse(1))
                     ilams(k) = 1./lams
                     L_qs(k) = .true.
                  else
                     rs(k) = 1.E-12
                     L_qs(k) = .false.
                  end if

                  if (qg1d(k) .gt. 1.E-9) then
                     rg(k) = qg1d(k)*rho(k)
                     if (ihail .eq. 1) then
                        N0_g(k) = xnoh
                     else
                        N0_g(k) = xnog
                     end if
                     lamg = (xam_g*xcgg(3)*N0_g(k)/rg(k))**(1./xcge(1))
                     ilamg(k) = 1./lamg
                     L_qg(k) = .true.
                  else
                     rg(k) = 1.E-12
                     L_qg(k) = .false.
                  end if
               end do

!+---+-----------------------------------------------------------------+
!..Locate K-level of start of melting (k_0 is level above).
!+---+-----------------------------------------------------------------+
               melti = .false.
               k_0 = kts
               do k = kte - 1, kts, -1
                  if ((temp(k) .gt. 273.15) .and. L_qr(k) &
                      .and. (L_qs(k + 1) .or. L_qg(k + 1))) then
                     k_0 = MAX(k + 1, k_0)
                     melti = .true.
                     goto 195
                  end if
               end do
195            continue

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
                  if (L_qs(k)) ze_snow(k) = (0.176/0.93)*(6.0/PIx)*(6.0/PIx) &
                                            *(xam_s/900.0)*(xam_s/900.0) &
                                            *N0_s(k)*xcsg(4)*ilams(k)**xcse(4)
                  if (L_qg(k)) ze_graupel(k) = (0.176/0.93)*(6.0/PIx)*(6.0/PIx) &
                                               *(xam_g/900.0)*(xam_g/900.0) &
                                               *N0_g(k)*xcgg(4)*ilamg(k)**xcge(4)
               end do

!+---+-----------------------------------------------------------------+
!..Special case of melting ice (snow/graupel) particles.  Assume the
!.. ice is surrounded by the liquid water.  Fraction of meltwater is
!.. extremely simple based on amount found above the melting level.
!.. Uses code from Uli Blahak (rayleigh_soak_wetgraupel and supporting
!.. routines).
!+---+-----------------------------------------------------------------+

               if (melti .and. k_0 .ge. kts + 1) then
                  do k = k_0 - 1, kts, -1

!..Reflectivity contributed by melting snow
                     if (L_qs(k) .and. L_qs(k_0)) then
                        fmelt_s = MAX(0.005d0, MIN(1.0d0 - rs(k)/rs(k_0), 0.99d0))
                        eta = 0.d0
                        lams = 1./ilams(k)
                        do n = 1, nrbins
                           x = xam_s*xxDs(n)**xbm_s
                           call rayleigh_soak_wetgraupel(x, DBLE(xocms), DBLE(xobms), &
                                                         fmelt_s, melt_outside_s, m_w_0, m_i_0, lamda_radar, &
                                                         CBACK, mixingrulestring_s, matrixstring_s, &
                                                         inclusionstring_s, hoststring_s, &
                                                         hostmatrixstring_s, hostinclusionstring_s)
                           f_d = N0_s(k)*xxDs(n)**xmu_s*DEXP(-lams*xxDs(n))
                           eta = eta + f_d*CBACK*simpson(n)*xdts(n)
                        end do
                        ze_snow(k) = SNGL(lamda4/(pi5*K_w)*eta)
                     end if

!..Reflectivity contributed by melting graupel

                     if (L_qg(k) .and. L_qg(k_0)) then
                        fmelt_g = MAX(0.005d0, MIN(1.0d0 - rg(k)/rg(k_0), 0.99d0))
                        eta = 0.d0
                        lamg = 1./ilamg(k)
                        do n = 1, nrbins
                           x = xam_g*xxDg(n)**xbm_g
                           call rayleigh_soak_wetgraupel(x, DBLE(xocmg), DBLE(xobmg), &
                                                         fmelt_g, melt_outside_g, m_w_0, m_i_0, lamda_radar, &
                                                         CBACK, mixingrulestring_g, matrixstring_g, &
                                                         inclusionstring_g, hoststring_g, &
                                                         hostmatrixstring_g, hostinclusionstring_g)
                           f_d = N0_g(k)*xxDg(n)**xmu_g*DEXP(-lamg*xxDg(n))
                           eta = eta + f_d*CBACK*simpson(n)*xdtg(n)
                        end do
                        ze_graupel(k) = SNGL(lamda4/(pi5*K_w)*eta)
                     end if

                  end do
               end if

               do k = kte, kts, -1
                  dBZ(k) = 10.*log10((ze_rain(k) + ze_snow(k) + ze_graupel(k))*1.d18)
               end do

            end subroutine refl10cm_gsfc

!+---+-----------------------------------------------------------------+

!JJS 20140225
! Calculate cloud droplet effective radius
            real function eff_rad_gpu(lambda)
            !$acc routine seq

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
               real, intent(in)  :: lambda   ! intercept parameter [1/cm]
!      real,intent(out) :: re  !effective radius [micron]

!
! for no particles.
!
!       if ( lambda <= 0.e0  .or. isnan(lambda) ) then
               if (lambda <= 0.e0) then
                  eff_rad_gpu = 0.e0
                  return
               end if

!
! compute drop effective radius for exponential distribution N(D) = N0*exp(-lam*D)
!
               eff_rad_gpu = 1.5e0/(lambda*100.)*1.0e+6  ! [micron]
               return
            end function eff_rad_gpu

!-------------------------------------------------------------------
            subroutine vti_mks_gpu(improve, rhoz, tz, qiz, qvz, p, xland, vti, &
                                   sat_predict, new_saturation, cimin, t0, &
                                   rhoe_s, itble, thrd, c218, c580, c76)
               !$acc routine seq
               !$acc routine(f_qsi_gpu) seq
               !$acc routine(esi_mks_gpu) seq
               implicit none
               integer, intent(in) :: improve
               real, intent(in) :: rhoz  !air density (kg/m^3)
               real, intent(in) :: qiz   !mixing ratio of cloud ice (kg/kg)
               real, intent(in) :: qvz   !specific humidity (kg/kg)
               real, intent(in) :: p     !pressure (Pa)
               real, intent(in) :: tz    !air temperature (K)
               real, intent(in) :: xland !land-sea mask
               real, intent(out):: vti   !terminal velocity (m/s)
               logical, intent(in) :: sat_predict, new_saturation
               real, intent(in) :: cimin, t0, rhoe_s, itble(0:120), &
                                   thrd, c218, c580, c76

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
               data vice/5, 15, 30, 35, 40, 45, 50/
               ! for Deng and Mace (2008)
               real, parameter :: aa = -4.14122e-5
               real, parameter :: bb = -0.00538922
               real, parameter :: cc = -0.0516344
               real, parameter :: dd = 0.00216078
               real, parameter :: ee = 1.9714
               real    :: tc
               real    :: h1, h2
               ! for TCWA1 semi-theoretical approach
               integer :: hid
               real    :: qsi, sqrhoz, rhoi, adagr, inhgr, ltk, lqi, &
                          ltk2, lqi2, zeta, vishp, viroi, ssi, lroi
               real    :: esi !new_saturation

               real, parameter :: di0 = 6.e-6

               if (qiz .ge. cimin) then
                  if (vtiflag .eq. 1) then
                     ! Starr and Cox 1985 :
                     y1 = rhoz*1000.*qiz           ! y1 in g/m^3
                     if (y1 .lt. 1.e-6) then
                        vti = 0.
                     else
                        do ic = 1, 6
                           if (y1 .gt. aice(ic) .and. y1 .le. aice(ic + 1)) then
                              vti = vice(ic) + (vice(ic + 1) - vice(ic))* &
                                    (y1 - aice(ic))/(aice(ic + 1) - aice(ic))
                              if (vti .lt. 0.0) vti = 0. ! EMK per 20110816 code
                           end if
                           vti = vti*0.01  ! convert back to MKS
                        end do
                     end if

                  elseif (vtiflag .eq. 2) then
                     ! Heymsfield and Donner 1990 :
                     vti = 3.29*(rhoz*qiz)**0.16

                  elseif (vtiflag .eq. 3) then
                     ! Hong et al. 2004, same as HD90 (in a*D**b form)
                     y1 = rhoz*1000.*qiz                  ! y1 in g/m^3
                     if (y1 .lt. 1.e-6) then
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
                        y1 = y1*1.e-3                       ! y1 in kg/m^3
                        bb1 = const_m*y1**0.25
                        bb2 = const_d*bb1**0.5
                        vti = max(const_vt*bb2**1.31, 0.0)    ! vti in m/s
                     end if

                  elseif (vtiflag .eq. 4) then
                     ! Deng and Mace 2008 :
                     ! IWC in g/m^3 , tc in oC , and vti in cm/s
                     tc = tz - t0
                     vti = (3.+log10(qiz*rhoz))* &
                           (tc*(aa*tc + bb) + cc) + dd*tc + ee
                     vti = exp(log(10.)*vti)
                     vti = vti*0.01    ! convert back to MKS

                  elseif (vtiflag .eq. 5) then
                     ! Mitchell et al. 2011 , equation 10 :
                     ! IWC in mg/m^3 , T in oC , and vti in cm/s
                     vti = 82.082 + 1.0121*tc + &
                           6.6303*(6.+log10(qiz*rhoz))
                     vti = vti*0.01    ! convert back to MKS

                  elseif (vtiflag .eq. 6) then
                     ! assume constant density (rhoi=300 kg/m^3) :
                     tc = tz - t0
                     sqrhoz = sqrt(rhoe_s/rhoz)

                     ! shape parameter :
                     hid = max(min(nint(abs(tc)/0.25), 120), 0)
                     adagr = itble(hid)
                     zeta = (adagr - 1.)/(adagr + 2.)
                     if (zeta .gt. 0.) then
                        vishp = di0**(zeta/2.)
                     elseif (zeta .lt. 0.) then
                        vishp = di0**(-zeta)
                     else
                        vishp = 1.
                     end if

                     ltk = log(tz)
                     ltk2 = ltk*ltk
                     lqi = -1.*log(rhoz*qiz)
                     lqi2 = lqi*lqi
                     vti = exp(-1.1100279E-2 + 0.47727519*lqi &
                               - 8.8757389E-2*lqi2 + 3.6732918E-3*lqi*lqi2 &
                               - 4.5748034E-5*lqi2*lqi2 + 1.3864255*ltk) &
                           /1.e+3*vishp*sqrhoz

                  elseif (vtiflag .eq. 7) then
                     ! TCWA1 adopted semi-theoretical approach from NTU 4ICE-3M scheme
                     ! with prescribed ice properties(shape and density) :
                     tc = tz - t0
                     qsi = f_qsi_gpu(tz, p, c218, c580, c76)
                     if (new_saturation) then
                        esi = min(0.99*p, esi_mks_gpu(tz))
                        qsi = 0.622*esi/(p - esi)
                     end if

                     ! deposition density :
                     if (tc .ge. -40.) then
                        ! Chen and Lamb 1994a ; Chen and Tsai 2016
                        hid = max(min(nint(abs(tc)/0.25), 120), 0)
                        inhgr = itble(hid)  !inherent growth ratio
                        rhoi = 900.*exp(-3.*max(qvz - qsi - 5.e-5, 0.)/inhgr)
                     elseif (tc .lt. -65.) then
                        inhgr = 1.
                        rhoi = 900.
                     else
                        ! Pokrifka et al. 2023, equation 18
!                        inhgr = 3.          !inherent growth ratio (FIG. 16)
                        if ( tc .ge. -50. ) then  !-50<=Tc<-40
                           inhgr = 3.
                        else                      !-65<=Tc<-50
                           inhgr = 2.
                        endif
                        ssi = qvz/qsi - 1.
                        if (ssi .gt. 0.267) then
                           rhoi = -1027.456*ssi + 1185.834
                        else
                           rhoi = -32.332*ssi + 900.
                        end if
                     end if

                     ! shape parameter (ice aspect ratio) :
                     adagr = inhgr**thrd
! >>> reduce upper-level vti :
!                     if (tc .ge. -40.) then
!                        adagr = inhgr**thrd
!                     else
!                        if (sat_predict) then
!                           adagr = inhgr**thrd
!                        else !sat_predict
!                           adagr = inhgr**0.9
!                        end if !sat_predict
!                     end if
! <<<
                     ltk = log(tz)
                     ltk2 = ltk*ltk
                     lqi = -1.*log(rhoz*qiz)
                     lqi2 = lqi*lqi
                     lroi = log(min(max(rhoi, 50.), 900.))
                     zeta = (adagr - 1.)/(adagr + 2.)
                     if (zeta .gt. 0.) then
                        vishp = di0**(zeta/2.)
                     elseif (zeta .lt. 0.) then
                        vishp = di0**(-zeta)
                     else
                        vishp = 1.
                     end if

                     ! density parameter :
                     viroi = 2.6795546 - 0.010732829*lqi - 1.176491*lroi &
                             + 3.2512268E-4*lqi2 + 0.13920371*lroi**2. &
                             -5.6243169E-4*lqi*lroi
                     viroi = min(1., viroi)

                     ! "NO" divide between land and ocean :
                     vti = exp(259.25629 - 0.26367743*lqi - 4.184759E-3*lqi2 &
                           - 91.567622*ltk + 8.6164869*ltk2) &
                           /1.e+6*vishp*viroi*(1.0837/rhoz)**0.35

                     ! the divide between land and ocean :
!                     if (xland .eq. 1.) then   ! land
!                        vti = exp(265.16805 - 0.28802545*lqi - 3.754874E-3*lqi2 &
!                                  - 93.843132*ltk + 8.8315122*ltk2) &
!                              /1.e+6*vishp*viroi*(1.0837/rhoz)**0.35
!                     else                        ! ocean
!                        vti = exp(252.86312 - 0.23844613*lqi - 4.6185936E-3*lqi2 &
!                                  - 89.119066*ltk + 8.3851678*ltk2) &
!                              /1.e+6*vishp*viroi*(1.0837/rhoz)**0.35
!                     end if
                  end if

!                  vti = min(vimax, max(vimin, vti))
                  vti = max(vimin, vti)   ! remove upper limit of vti
               else
                  vti = vimin
               end if  !end of if qiz

               return
            end subroutine vti_mks_gpu
!-------------------------------------------------------------------

!-------------------------------------------------------------------
            subroutine vtr_mks_gpu(rhoz, qrz, tz, vtr, constb, consta, &
                                   rhoe_s, crmin, rhowater, xnor, vrc0, vrc1, &
                                   vrc2, vrc3, draimax, tnw, roqr, t0)
               !$acc routine seq
               !$acc routine(gammagce_gpu) seq
               implicit none
               real, intent(in) :: rhoz  !air density (kg/m^3)
               real, intent(in) :: qrz   !mixing ratio (kg/kg)
               real, intent(in) :: tz    !temperature (K)
               real, intent(out):: vtr   !terminal velocity (m/s)
               real, intent(in) :: constb, consta, rhoe_s, crmin, rhowater, &
                  xnor, vrc0, vrc1, vrc2, vrc3, draimax, tnw, roqr, t0 ! acc routine

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
               gambp4 = gammagce_gpu(constb + 4.)

               if (qrz .ge. crmin) then
                  fv = sqrt(rhoe_s/rhoz)
!         igce = 1
                  if (vtrflag .eq. 1) then !if igce .ne. 1
                     ! Lin et al. (1983) ; in MKS
                     ! rhowater=1000., consta=0.8, constb=841.9967
!            tmp1 = sqrt(pi*rhowater*xnor/rhoz/qrz)
!            tmp1 = sqrt(tmp1)
!            vtr = consta*gambp4*fv/tmp1**constb
!            vtr = vtr/6.
                     zr = (pi*rhowater*xnor/(rhoz*qrz))**0.25   !slope parameter of rain (in MKS)
                     gambp4 = gammagce_gpu(constb + 4.)
                     vtr = consta*gambp4/(6.*zr**constb)*fv

                  elseif (vtrflag .eq. 2) then !if igce=1
                     ! new codes from Steve's in CGS :
                     y1 = rhoz*qrz*0.001  !in CGS
                     vs = sqrt(y1)
                     vg = sqrt(vs)
                     vr = vrc0 + vrc1*vg + vrc2*vs + vrc3*vg*vs
                     vtr = max(fv*vr, 0.e0)
                     vtr = vtr*0.01  !convert back to MKS

                  elseif (vtrflag .eq. 3) then
                     ! same as Steve's, but with bin rain evaporation correction (Lang et al. 2014)
                     ! (in CGS)
                     ftnw = 1.
                     if (tz .gt. t0) then
                        bin_factor = 0.11*(1000.*qrz)**(-1.27) + 0.98
                        bin_factor = min(bin_factor, 1.30)
                        ftnw = 1./bin_factor**3.35
                        ftnwmin = rhoz*0.001*qrz/draimax
                        if (qrz .le. 0.001) ftnw = max(ftnw, ftnwmin/tnw)
                     end if
                     zr = (pi*roqr*tnw/(rhoz*0.001*qrz))**0.25   !slope parameter of rain (in CGS)
                     zr = zr/ftnw
                     vtr = (-26.7 + 2.06e+4/zr - 2.045e+5/zr**2.+9.06e+5/zr**3.)*fv
                     vtr = vtr*0.01  !convert back to MKS

                  end if  !end of if igce
                  vtr = min(vrmax, max(vrmin, vtr))
               else
                  vtr = vrmin
               end if  !end of if qrz

               return
            end subroutine vtr_mks_gpu
!-------------------------------------------------------------------

!-------------------------------------------------------------------
            subroutine vts_mks_gpu(improve, rhoz, qsz, tz, vts, rhoe_s, constd, &
                                   csmin, rhosnow, xnos, constc, vsc, t0, roqs, &
                                   roqg, tns, tng, bsq)
               !$acc routine seq
               !$acc routine(gammagce_gpu) seq
               !$acc routine(sgmap_gpu) seq
               implicit none
               integer, intent(in) :: improve
               real, intent(in) :: rhoz  !air density (kg/m^3)
               real, intent(in) :: qsz   !mixing ratio (kg/kg)
               real, intent(in) :: tz    !temperature (K)
               real, intent(out):: vts   !terminal velocity (m/s)
               real, intent(in) :: rhoe_s, constd, csmin, rhosnow, xnos, constc, &
                  vsc, t0, roqs, roqg, tns, tng, bsq ! acc routine


               real, parameter :: vsmax = 5.0     !(m/s)
               real, parameter :: vsmin = 0.0     !(m/s)

               integer, parameter :: vtsflag = 2
               ! 1 : Lin et al. (1983)   ( igce != 1 )
               ! 2 : Steve's with size mapping  (Lang el al. 2011)

               real    :: y1, fv, r00, ftns, ftns0, tzc, vscf
               real    :: tmp1, pi, gamdp4

               pi = acos(-1.)
               gamdp4 = gammagce_gpu(constd + 4.)

               if (qsz .ge. csmin) then
                  y1 = rhoz*qsz*0.001  !in CGS
                  fv = sqrt(rhoe_s/rhoz)
!         igce = 1
                  if (vtsflag .eq. 1) then
                     ! Lin et al. (1983) ; in MKS
                     tmp1 = sqrt(pi*rhosnow*xnos/rhoz/qsz)
                     tmp1 = sqrt(tmp1)
                     vts = constc*gamdp4*fv/tmp1**constd
                     vts = vts/6.
                  elseif (vtsflag .eq. 2) then
                     ! Steve's with size mapping ; in CGS
                     y1 = qsz
                     r00 = rhoz*0.001   ! rhoz need to be in CGS
                     vscf = vsc*fv
                     tzc = tz - t0

                     ftns = 1.
                     ftns0 = 1.
                     call sgmap_gpu(1, y1, r00, tzc, ftns0, roqs, roqg, tns, tng)
                     ftns = ftns0**bsq

                     vts = max(vscf*(r00*y1)**bsq/ftns, 0.0)
                     vts = vts*0.01  ! convert back to MKS
                  end if  !end of if vtsflag
                  vts = min(vsmax, max(vsmin, vts))
               else
                  vts = vsmin
               end if  !end of if qsz

               return
            end subroutine vts_mks_gpu
!-------------------------------------------------------------------

!-------------------------------------------------------------------
            subroutine vtg_mks_gpu(ihail, improve, rhoz, qgz, tz, vtg, cgmin, &
                                   rhohail, xnoh, cdrag, rhoe_s, vgc, roqs, &
                                   roqg, tns, tng, bgq, t0)
               !$acc routine seq
               !$acc routine(gammagce_gpu) seq
               !$acc routine(sgmap_gpu) seq
               implicit none
               integer, intent(in) :: ihail, improve
               real, intent(in) :: rhoz  !air density (kg/m^3)
               real, intent(in) :: qgz   !mixing ratio (kg/kg)
               real, intent(in) :: tz    !temperature (K)
               real, intent(out):: vtg   !terminal velocity (m/s)
               real, intent(in) :: cgmin, rhohail, xnoh, cdrag, rhoe_s, vgc, &
                  roqs, roqg, tns, tng, bgq, t0 ! acc routine

               real, parameter :: vgmax = 8.0     !(m/s)
               real, parameter :: vgmin = 0.0     !(m/s)

               !local variables
               integer :: igce
               real    :: tmp1, pi, gam4pt5, term0, grav, r00, vgcr
               real    :: y1, fv, ftng, ftng0, tzc

               pi = acos(-1.)
               grav = 9.80665e+0
               gam4pt5 = gammagce_gpu(4.5)

               if (qgz .ge. cgmin) then
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
                     fv = sqrt(rhoe_s/rhoz)
                     vgcr = vgc*fv
                     ftng = 1.
                     ftng0 = 1.
                     if (improve .gt. 2) then
                        tzc = tz - t0
                        call sgmap_gpu(2, y1, r00, tzc, ftng0, roqs, roqg, tns, tng)
                        ftng = ftng0**bgq
                     end if
                     vtg = dmax1(vgcr*(r00*y1)**bgq/ftng, 0.0)
                     vtg = vtg*0.01  !convert back to MKS

!            endif  !end of if igce
                     vtg = min(vgmax, max(vgmin, vtg))
                  end if  !end of if ihail
               else
                  vtg = vgmin
               end if   !end of if qgz

               return
            end subroutine vtg_mks_gpu
!-------------------------------------------------------------------

!-------------------------------------------------------------------
            SUBROUTINE semi_sedi_gpu(myim, qvar, its, ite, jts, jte, &
                       ihail, improve, iter, km, dzl, rho, qc, &
                       qv, p, tz, xland, ww, precip, dt, sat_predict, &
                       new_saturation, async_id)
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
               integer, intent(in) :: myim(my_max)
               integer, intent(in) :: ihail, improve, km, iter, its, ite, jts, jte
               real, intent(in) ::  dt
               real, intent(in) :: dzl(its:ite, km, jts:jte), &
                                   rho(its:ite, km, jts:jte), &
                                   tz(its:ite, km, jts:jte), &
                                   qv(its:ite, km, jts:jte), &
                                   p(its:ite, km, jts:jte)
               real, intent(in) :: xland(its:ite, jts:jte)
               real, intent(out) :: ww(its:ite, km, jts:jte)
               real, intent(inout) :: precip(its:ite, jts:jte)
               real, intent(inout) :: qc(its:ite, km, jts:jte)
               logical, intent(in) :: sat_predict, new_saturation
               integer k, m, kk, kb, kt, n, i, j
               real tl, tl2, qql, dql, qqd
               real th, th2, qqh, dqh
               real zsum, qsum, dim, dip, con1, fa1, fa2
               real allold, decfl
               real, dimension(its:ite, km, jts:jte) :: qq, qn, wd, wa, was, qcold
               real za(its:ite, km + 2, jts:jte)
               real, dimension(its:ite, km+1, jts:jte) :: wi, zi, dza, qa, qmi, qpi
               ! for GPU porting
               integer, intent(in) :: async_id
               real :: dzarb, dzart, qmirb, qmirt, zir, zir1
!
               fa1 = 9./16.
               fa2 = 1./16.
               con1 = 0.05
               !$acc data create(qq, qn, wd, wa, was, za, wi, zi, &
               !$acc&     dza, qa, qmi, qpi, qcold) async(async_id)
               !$acc parallel loop gang collapse(2) async(async_id)
               do j = jts, jte
                  do k = 1, km
                     !$acc loop vector
                     do i = its, myim(j)
                        if (k .eq. km) qa(i, km+1, j) = 0.0
                        qa(i, k, j) = 0.0
                        qq(i, k, j) = 0.0
                        wa(i, k, j) = 0.0
                        was(i, k, j) = 0.0
                        ww(i, k, j) = 0.0
                        qcold(i, k, j) = qc(i, k, j)
                        
                        if (qvar .eq. 'qr') call vtr_mks_gpu(rho(i, k, j), &
                           qc(i, k, j), tz(i, k, j), ww(i, k, j), constb, &
                           consta, rhoe_s, crmin, rhowater, xnor, vrc0, vrc1, &
                           vrc2, vrc3, draimax, tnw, roqr, t0)
                        if (qvar .eq. 'qs') call vts_mks_gpu(improve, rho(i, k, j), &
                           qc(i, k, j), tz(i, k, j), ww(i, k, j), rhoe_s, constd, &
                           csmin, rhosnow, xnos, constc, vsc, t0, roqs, roqg, &
                           tns, tng, bsq)
                        if (qvar .eq. 'qg') call vtg_mks_gpu(ihail, improve, &
                           rho(i, k, j), qc(i, k, j), tz(i, k, j), ww(i, k, j), &
                           cgmin, rhohail, xnoh, cdrag, rhoe_s, vgc, roqs, &
                           roqg, tns, tng, bgq, t0)
                        if (qvar .eq. 'qi') call vti_mks_gpu(improve, rho(i, k, j), &
                           tz(i, k, j), qc(i, k, j), qv(i, k, j), p(i, k, j), &
                           xland(i, j), ww(i, k, j), sat_predict, new_saturation, &
                           cimin, t0, rhoe_s, itble, thrd, c218, c580, c76)
                        qq(i, k, j) = qc(i, k, j)*rho(i, k, j)
                        wd(i, k, j) = ww(i, k, j)
                     end do
                  end do
               end do

      ! plm is 2nd order, we can use 2nd order wi or 3rd order wi
      ! 2nd order interpolation to get wi
                     !wi(i, 1, j) = ww(i, 1, j)
                     !wi(i, km + 1, j) = ww(i, km, j)
                     !do k = 2, km
                     !   wi(i, k, j) = (ww(i, k, j)*dzl(i, k - 1, j) + &
                     !                 ww(i, k - 1, j)*dzl(i, k, j))/(dzl(i, k - 1, j) + dzl(i, k, j))
                     !end do
                     
               !$acc parallel loop gang collapse(2) async(async_id)
               do j = jts, jte
                  do k = 2, km
                     !$acc loop vector
                     do i = its, myim(j)
                        !====================================
                        ! save departure wind
                        ! 3rd order interpolation to get wi
                        if (k .eq. 2) then
                           wi(i, 1, j) = ww(i, 1, j)
                           wi(i, 2, j) = 0.5*(ww(i, 2, j) + ww(i, 1, j))
                        elseif (k .eq. km) then
                           wi(i, km, j) = 0.5*(ww(i, km, j) + ww(i, km - 1, j))
                           wi(i, km + 1, j) = ww(i, km, j)
                        else
                           wi(i, k, j) = fa1*(ww(i, k, j) + ww(i, k - 1, j)) &
                                         - fa2*(ww(i, k + 1, j) + ww(i, k - 2, j))
                        end if
                        ! terminate of top of raingroup
                        if (k .ge. 2) then
                           if (ww(i, k, j) .eq. 0.0) wi(i, k, j) = ww(i, k - 1, j)
                        end if
                     end do
                  end do
               end do

               !$acc parallel loop gang vector collapse(2) async(async_id) &
               !$acc&         private(decfl)
               do j = jts, jte
                  do i = its, ite
                     if (i .le. myim(j)) then
      !
      ! compute interface values
                        zi(i, 1, j) = 0.0
                        !$acc loop seq
                        do k = 1, km
                           zi(i, k + 1, j) = zi(i, k, j) + dzl(i, k, j)
                        end do
                        
         ! diffusivity of wi
                        !$acc loop seq
                        do k = km, 1, -1
                           decfl = (wi(i, k + 1, j) - wi(i, k, j))*dt/dzl(i, k, j)
                           if (decfl .gt. con1) then
                              wi(i, k, j) = wi(i, k + 1, j) - con1*dzl(i, k, j)/dt
                           end if
                        end do
                     end if
                  end do
               end do
      ! compute arrival point
               !$acc parallel loop gang collapse(2) async(async_id)
               do j = jts, jte
                  do k = 1, km+1
                     !$acc loop vector
                     do i = its, myim(j)
                        za(i, k, j) = zi(i, k, j) - wi(i, k, j)*dt
                        if (k .eq. (km + 1)) za(i, km + 2, j) = zi(i, km + 1, j)
                     end do
                  end do
               end do
                     
               !$acc parallel loop gang collapse(2) async(async_id)
               do j = jts, jte
                  do k = 1, km+1
                     !$acc loop vector
                     do i = its, myim(j)
                        dza(i, k, j) = za(i, k + 1, j) - za(i, k, j)
      ! computer deformation at arrival point
                        if (k .le. km) then
                           qa(i, k, j) = qq(i, k, j)*dzl(i, k, j)/dza(i, k, j)
                           qc(i, k, j) = qa(i, k, j)/rho(i, k, j)
                        else
                           qa(i, km + 1, j) = 0.0
                        end if
                     end do
                  end do
               end do


      ! compute arrival terminal velocity, and estimate mean terminal velocity
      ! then back to use mean terminal velocity
               do n = 1, iter
                  !$acc parallel loop gang collapse(2) async(async_id)
                  do j = jts, jte
                     do k = 1, km
                        !$acc loop vector
                        do i = its, myim(j)
                           if (qvar .eq. 'qr') call vtr_mks_gpu(rho(i, k, j), &
                              qc(i, k, j), tz(i, k, j), wa(i, k, j), constb, &
                              consta, rhoe_s, crmin, rhowater, xnor, vrc0, &
                              vrc1, vrc2, vrc3, draimax, tnw, roqr, t0)
                           if (qvar .eq. 'qs') call vts_mks_gpu(improve, &
                              rho(i, k, j), qc(i, k, j), tz(i, k, j), wa(i, k, j), &
                              rhoe_s, constd, csmin, rhosnow, xnos, constc, vsc, &
                              t0, roqs, roqg, tns, tng, bsq)
                           if (qvar .eq. 'qg') call vtg_mks_gpu(ihail, improve, &
                              rho(i, k, j), qc(i, k, j), tz(i, k, j), wa(i, k, j), &
                              cgmin, rhohail, xnoh, cdrag, rhoe_s, vgc, roqs, &
                              roqg, tns, tng, bgq, t0)
                           if (qvar .eq. 'qi') call vti_mks_gpu(improve, rho(i, k, j), &
                              tz(i, k, j), qc(i, k, j), qv(i, k, j), p(i, k, j), &
                              xland(i, j), wa(i, k, j), sat_predict, new_saturation, &
                              cimin, t0, rhoe_s, itble, thrd, c218, c580, c76)
                           if (n .ge. 2) wa(i, k, j) = 0.5*(wa(i, k, j) + was(i, k, j))
                           ! mean wind is average of departure and new arrival winds
                           ww(i, k, j) = 0.5*(wd(i, k, j) + wa(i, k, j))
                           was(i, k, j) = wa(i, k, j)
                        end do
                     end do
                  end do

                  !$acc parallel loop gang collapse(2) async(async_id)
                  do j = jts, jte
                     do k = 1, km
                        !$acc loop vector
                        do i = its, myim(j)
                           ! 3rd order interpolation to get wi
                           if (k .eq. 1) then
                              wi(i, 1, j) = ww(i, 1, j)
                           elseif (k .eq. 2) then
                              wi(i, 2, j) = 0.5*(ww(i, 2, j) + ww(i, 1, j))
                           elseif (k .eq. km) then
                              wi(i, km, j) = 0.5*(ww(i, km, j) + ww(i, km - 1, j))
                              wi(i, km + 1, j) = ww(i, km, j)
                           else
                              wi(i, k, j) = fa1*(ww(i, k, j) + ww(i, k - 1, j)) &
                                            - fa2*(ww(i, k + 1, j) + ww(i, k - 2, j))
                           end if
                           ! terminate of top of raingroup
                           if (k .ge. 2) then
                              if (ww(i, k, j) .eq. 0.0) wi(i, k, j) = ww(i, k - 1, j)
                           end if
                        end do
                     end do
                  end do

                  !$acc parallel loop gang vector collapse(2) async(async_id) &
                  !$acc&         private(decfl)
                  do j = jts, jte
                     do i = its, ite
                        if (i .le. myim(j)) then
         ! diffusivity of wi
                           !$acc loop seq
                           do k = km, 1, -1
                              decfl = (wi(i, k + 1, j) - wi(i, k, j))*dt/dzl(i, k, j)
                              if (decfl .gt. con1) then
                                 wi(i, k, j) = wi(i, k + 1, j) - con1*dzl(i, k, j)/dt
                              end if
                           end do
                        end if
                     end do
                  end do
         ! compute arrival point
                  !$acc parallel loop gang collapse(2) async(async_id)
                  do j = jts, jte
                     do k = 1, km+1
                        !$acc loop vector
                        do i = its, myim(j)
                           za(i, k, j) = zi(i, k, j) - wi(i, k, j)*dt
                           if (k .eq. (km + 1)) za(i, km + 2, j) = zi(i, km + 1, j)
                        end do
                     end do
                  end do
                        

                  !$acc parallel loop gang collapse(2) async(async_id)
                  do j = jts, jte
                     do k = 1, km+1
                        !$acc loop vector
                        do i = its, myim(j)
                           dza(i, k, j) = za(i, k + 1, j) - za(i, k, j)
         ! computer deformation at arrival point
                           if (k .le. km) then
                              qa(i, k, j) = qq(i, k, j)*dzl(i, k, j)/dza(i, k, j)
                              qc(i, k, j) = qa(i, k, j)/rho(i, k, j)
      ! interpolation to regular point
                              qn(i, k, j) = 0.0
                           else
                              qa(i, km + 1, j) = 0.0
                           end if
                        end do
                     end do
                  end do
               end do
      !====================================
               !$acc parallel loop gang collapse(2) async(async_id) &
               !$acc&         private(dip, dim)
               do j = jts, jte
                  do k = 1, km+1
                     !$acc loop vector
                     do i = its, myim(j)
                        if (k .eq. 1) then
                           qpi(i, 1, j) = qa(i, 1, j)
                           qmi(i, 1, j) = qa(i, 1, j)
                        elseif (k .eq. km+1) then
                           qmi(i, km + 1, j) = qa(i, km + 1, j)
                           qpi(i, km + 1, j) = qa(i, km + 1, j)
                        else
         ! estimate values at arrival cell interface with monotone
                           dip = (qa(i, k + 1, j) - qa(i, k, j)) &
                                 /(dza(i, k + 1, j) + dza(i, k, j))
                           dim = (qa(i, k, j) - qa(i, k - 1, j)) &
                                 /(dza(i, k - 1, j) + dza(i, k, j))
                           if (dip*dim .le. 0.0) then
                              qmi(i, k, j) = qa(i, k, j)
                              qpi(i, k, j) = qa(i, k, j)
                           else
                              qpi(i, k, j) = qa(i, k, j) + 0.5*(dip + dim)*dza(i, k, j)
                              qmi(i, k, j) = 2.0*qa(i, k, j) - qpi(i, k, j)
                              if (qpi(i, k, j) .lt. 0.0 .or. qmi(i, k, j) .lt. 0.0) then
                                 qpi(i, k, j) = qa(i, k, j)
                                 qmi(i, k, j) = qa(i, k, j)
                              end if
                           end if
                        end if
                     end do
                  end do
               end do
               
               !$acc parallel loop gang vector collapse(2) async(async_id) &
               !$acc&         private(kb, kt, tl, th, tl2, th2, qqd, qqh, qql, &
               !$acc&         dql, zsum, qsum, dqh, dzarb, dzart, qmirb, qmirt, &
               !$acc&         zir, zir1)
               do j = jts, jte
                  do i = its, ite
                     if (i .le. myim(j)) then
                        kb = 1
                        kt = 1
                        !$acc loop seq
                        intp: do k = 1, km
                           kb = max(kb - 1, 1)
                           kt = max(kt - 1, 1)
         ! find kb and kt
                           zir = zi(i, k, j)
                           zir1 = zi(i, k + 1, j)
                           if (zir .ge. za(i, km + 1, j)) then
                              exit intp
                           else
                              find_kb: do kk = kb, km
                                 if (zir .le. za(i, kk + 1, j)) then
                                    kb = kk
                                    exit find_kb
                                 else
                                    cycle find_kb
                                 end if
                              end do find_kb
                              find_kt: do kk = kt, km + 2
                                 if (zir1 .le. za(i, kk, j)) then
                                    kt = kk
                                    exit find_kt
                                 else
                                    cycle find_kt
                                 end if
                              end do find_kt
                              kt = kt - 1
         ! compute q with piecewise constant method
                              dzarb = dza(i, kb, j)
                              dzart = dza(i, kt, j)
                              qmirb = qmi(i, kb, j)
                              qmirt = qmi(i, kt, j)
                              if (kt .eq. kb) then
                                 tl = (zir - za(i, kb, j))/dzarb
                                 th = (zir1 - za(i, kb, j))/dzarb
                                 tl2 = tl*tl
                                 th2 = th*th
                                 qqd = 0.5*(qpi(i, kb, j) - qmirb)
                                 qqh = qqd*th2 + qmirb*th
                                 qql = qqd*tl2 + qmirb*tl
                                 qn(i, k, j) = (qqh - qql)/(th - tl)
                              else if (kt .gt. kb) then
                                 tl = (zir - za(i, kb, j))/dzarb
                                 tl2 = tl*tl
                                 qqd = 0.5*(qpi(i, kb, j) - qmirb)
                                 qql = qqd*tl2 + qmirb*tl
                                 dql = qa(i, kb, j) - qql
                                 zsum = (1.-tl)*dzarb
                                 qsum = dql*dzarb
                                 if (kt - kb .gt. 1) then
                                    do m = kb + 1, kt - 1
                                       zsum = zsum + dza(i, m, j)
                                       qsum = qsum + qa(i, m, j)*dza(i, m, j)
                                    end do
                                 end if
                                 th = (zir1 - za(i, kt, j))/dzart
                                 th2 = th*th
                                 qqd = 0.5*(qpi(i, kt, j) - qmirt)
                                 dqh = qqd*th2 + qmirt*th
                                 zsum = zsum + th*dzart
                                 qsum = qsum + dqh*dzart
                                 qn(i, k, j) = qsum/zsum
                              end if
                              cycle intp
                           end if
                        end do intp

         ! skip for no precipitation for all layers
                        allold = 0.0
                        !$acc loop seq
                        do k = 1, km
                           allold = allold + qq(i, k, j)
                        end do

      ! rain out (unit:kg/m^2=mm)
                        if (allold .gt. 0.0) then
                           !$acc loop seq
                           sum_precip: do k = 1, km
                              if (za(i, k, j) .lt. 0.0 &
                                  .and. za(i, k + 1, j) .le. 0.0) then
                                 precip(i, j) = precip(i, j) + qa(i, k, j)*dza(i, k, j)
                                 cycle sum_precip
                              else if (za(i, k, j) .lt. 0.0 &
                                       .and. za(i, k + 1, j) .gt. 0.0) then
                                 ! from Thompson MP :
                                 th = (0.0 - za(i, k, j))/dza(i, k, j)
                                 th2 = th*th
                                 qqd = 0.5*(qpi(i, k, j) - qmi(i, k, j))
                                 qqh = qqd*th2 + qmi(i, k, j)*th
                                 precip(i, j) = precip(i, j) + qqh*dza(i, k, j)
            !             ! from WSM6 MP :
            !             precip = precip + qa(k)*(0.-za(k))
                                 exit sum_precip
                              end if
                              exit sum_precip
                           end do sum_precip
                           !$acc loop seq
                           do k = 1, km
                              qc(i, k, j) = qn(i, k, j)/rho(i, k, j) ! replace the new values
                           enddo
                        else
                           !$acc loop seq
                           do k = 1, km
                              ww(i, k, j) = 0.0 ! cycle, do not modify out/inout argument
                              qc(i, k, j) = qcold(i, k, j)
                           enddo
                        end if
                     end if
                  end do
               end do

      
               !$acc end data

            END SUBROUTINE semi_sedi_gpu

            REAL FUNCTION f_qsi_gpu(tair, p0, c218, c580, c76)
               !$acc routine seq
               IMPLICIT NONE
               REAL :: tair   !real temperature (K)
               REAL :: p0     !pressure (Pa)
               REAL :: rp0
               real, intent(in) :: c218, c580, c76

               rp0 = 3.799052e3/(p0*10.)
               f_qsi_gpu = rp0*exp(c218 - c580/(tair - c76))

            END FUNCTION f_qsi_gpu


            REAL FUNCTION esi_mks_gpu(tair)
            !$acc routine seq
               IMPLICIT NONE
               REAL :: tair   !real temperature (K)
!
!  COMPUTE SATURATION VAPOR PRESSURE POLYSVP RETURNED IN UNITS OF PA. T IS INPUT IN UNITS OF K.
!  REPLACE GOFF-GRATCH WITH FASTER FORMULATION FROM FLATAU ET AL. 1992, TABLE 4 (RIGHT-HAND COLUMN)
               REAL :: DT
               REAL :: a0i, a1i, a2i, a3i, a4i, a5i, a6i, a7i, a8i
               DATA a0i, a1i, a2i, a3i, a4i, a5i, a6i, a7i, a8i/6.11147274, 0.503160820, &
                  0.188439774E-1, 0.420895665E-3, 0.615021634E-5, 0.602588177E-7, &
                  0.385852041E-9, 0.146898966E-11, 0.252751365E-14/

               DT = MAX(-80., tair - 273.16)
               esi_mks_gpu = a0i + DT*(a1i + DT*(a2i + DT*(a3i + DT*(a4i + DT*(a5i + &
                                             DT*(a6i + DT*(a7i + a8i*DT)))))))
               esi_mks_gpu = esi_mks_gpu*100.  !convert to Pa
!
!  Goff-Gratch equation (Goff and Gratch 1945)
!       esi_mks = c610 * exp( c218 - c580 / (tair - c76) )
!       esi_mks = esi_mks / 10.  !convert to Pa

            END FUNCTION esi_mks_gpu

            REAL FUNCTION esw_mks_gpu(tair)
            !$acc routine seq
               IMPLICIT NONE
               REAL :: tair   !real temperature (K)
!
!  COMPUTE SATURATION VAPOR PRESSURE POLYSVP RETURNED IN UNITS OF PA. T IS INPUT IN UNITS OF K.
!  REPLACE GOFF-GRATCH WITH FASTER FORMULATION FROM FLATAU ET AL. 1992, TABLE 4 (RIGHT-HAND COLUMN)
               REAL :: DT
               REAL :: a0, a1, a2, a3, a4, a5, a6, a7, a8
               DATA a0, a1, a2, a3, a4, a5, a6, a7, a8/6.11239921, 0.443987641, &
                  0.142986287E-1, 0.264847430E-3, 0.302950461E-5, 0.206739458E-7, &
                  0.640689451E-10, -0.952447341E-13, -0.976195544E-15/

               DT = MAX(-80., tair - 273.16)
               esw_mks_gpu = a0 + DT*(a1 + DT*(a2 + DT*(a3 + DT*(a4 + DT*(a5 + &
                                           DT*(a6 + DT*(a7 + a8*DT)))))))
               esw_mks_gpu = esw_mks_gpu*100.  !convert to Pa

               ! to be closer to function fpvs :
               if ( DT.gt.-5.0 ) esw_mks_gpu = esw_mks_gpu*max(1.0-(DT+5.0)*0.0025/35.0,0.9975)
!
!  Goff-Gratch equation (Goff and Gratch 1945)
!       esw_mks = c610 * exp( c172 - c409 / (tair - c358) )
!       esw_mks = esw_mks / 10.   !convert to Pa

            END FUNCTION esw_mks_gpu

            END MODULE module_mp_gsfcgce_3ice_nuwrf_gpu

