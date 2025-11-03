#define EffectRad_GCE3
!--------------------------
      subroutine mp_init_gpu                                               &
!--------------------------
!  ---  inputs:
           ( nmmiph ,myrank)
!  ---  outputs: ( none )

! for WSM6
      use module_mp_wsm6,     only : wsm6init
! for Thompson
      use module_mp_thompson, only : thompson_init
! for 2M Thompson
      use module_mp_thompson_new,                                       &
                              only : new_thompson_init
! for GFDL MP v1
      use module_mp_gfdl,     only : gfdl_cloud_microphys_init
! for GFDL MP v2
      use module_mp_gfdl_v2,  only : gfdlv2_init => gfdl_cld_mp_init
! for GFDL MP v3
      use module_mp_gfdl_v3,  only : gfdlv3_init => gfdl_cld_mp_init
      use physpara, only : is_aerosol_aware,merra2_aerosol_aware
#ifdef Readaeroclx
! for GOCART coupling
      use module_gocart_coupling, only : makelut_ccn_icn
#endif

      implicit none
!  ---  input:
      integer,  intent(in) :: nmmiph,myrank
!  ---  local:
      integer   ntrac_req
      integer   errflg
      character errmsg
!-----------------------------------------------------------------------
!  for cloud microphysics
!-----------------------------------------------------------------------
! WSM6
        if ( nmmiph .eq. 6 ) then
          call wsm6init()
          if ( myrank .eq. 0 )                                         &
             print *,'WSM6 cloud microphysics initialized'
        endif
! Thompson
        if ( nmmiph .eq. 8 ) then
          call thompson_init()
          if ( myrank .eq. 0 )                                         &
             print *,'Thompson cloud microphysics initialized'
        endif
! 2M Thompson
        if ( nmmiph .eq. 18 ) then
          call new_thompson_init ( is_aerosol_aware,                   &
                   merra2_aerosol_aware, myrank, 0, errmsg, errflg )
          if ( myrank .eq. 0 )                                         &
             print *,'New Thompson cloud microphysics initialized'
        endif
! GFDL MP v1
        if ( nmmiph .eq. 11 ) then
          call gfdl_cloud_microphys_init()
          if ( myrank .eq. 0 )                                         &
             print *, 'GFDL cloud microphysics initialized'
        endif
! GFDL MP v2
        if ( nmmiph .eq. 12 ) then
          call gfdlv2_init()
          if ( myrank .eq. 0 )                                         &
             print *, 'GFDL cloud microphysics version 2 initialized'
        endif
! GFDL MP v3
        if ( nmmiph .eq. 13 ) then
          call gfdlv3_init()
          if ( myrank .eq. 0 )                                         &
             print *, 'GFDL cloud microphysics version 3 initialized'
        endif
! Goddard (GCE) 3ICE MP
        if ( nmmiph .eq. 15 ) then
          if ( myrank .eq. 0 )                                         &
             print *,'Goddard (GCE) 3ICE cloud microphysics initialized'
#ifdef Readaeroclx
          call makelut_ccn_icn
#endif
        endif
! Goddard (GCE) 4ICE MP
        if ( nmmiph .eq. 16 ) then
          if ( myrank .eq. 0 )                                         &
             print *,'Goddard (GCE) 4ICE cloud microphysics initialized'
        endif

      return
!--------------------------
      end subroutine mp_init_gpu
!--------------------------
!--------------------------
      subroutine mp_scheme_gpu                                         &
!--------------------------
!  ---  inputs:
           ( nmmiph,nx,nxjp,lev,ncld,plt,ptop,                         &
             dsigma,phii,islimsk,q0,kdt,tpi,me,dta,area,jjj,           &
             itimestep,sgeo,phi,rhc_mp,pk,                             &
             snr,xlat,sdec,ivegtyp,                                    &
#ifdef Readaeroclx
             aeroclx,naero,                                            &
#endif
!  ---  inputs/outputs:
             tt,qt,qa,ut,vt,vvel,pst,                                  &
!  ---  outputs:
             re_cloud,re_ice,re_snow,re_rain,                          &
             rlsp,rlspi,rlsps,rlspg,sr,                                &
!  ---  flags:
             SL_sedi, sat_predict, new_saturation,                     &
             use_cpm, use_declination)

      use index,               only: jlistnum, jlist1
      use param,               only: my, my_max
      use rank
      use radn,                only: ntcw,ntiw,ntrw,ntsw,ntgl,nthl,     &
                                     ntinc,ntrnc
! for wsm6
      use module_mp_wsm6,      only: wsm6
! for thompson
      use module_mp_thompson,  only: thompson_driver => mp_gt_driver
! for 2M Thompson
      use module_mp_thompson_new,                                       &
                               only: new_thompson_driver => mp_gt_driver&
                                     , cal_cldfra3, cfflag_thom
      use module_mp_thompson_make_number_concentrations,                &
                               only: make_IceNumber, make_RainNumber
! for GFDL MP v1
      use module_mp_gfdl,      only: gfdl_cloud_microphys_driver,       &
                                     cloud_diagnosis
! for GFDL MP v2
      use module_mp_gfdl_v2,   only: gfdlv2_driver =>gfdl_cld_mp_driver,&
                                     cloud_diagnosis_v2,                &
                                     sedi_w_v2 => do_sedi_w
! for GFDL MP v3
      use module_mp_gfdl_v3,   only: gfdlv3_driver =>gfdl_cld_mp_driver,&
                                     cloud_diagnosis_v3,                &
                                     sedi_w_v3 => do_sedi_w
! for Goddard (GCE) 3ICE MP
      use module_mp_gsfcgce_3ice_nuwrf_gpu, only: gsfcgce_3ice_nuwrf_gpu
!      use module_mp_gsfcgce_3ice_cwb  , only: gsfcgce_3ice_cwb
      use module_mp_gsfcgce,   only: gsfcgce
! for Goddard (GCE) 4ICE MP
      use module_mp_gce4ice,   only: gsfcgce_4ice_nuwrf
      use physcons,            only: con_rd,con_fvirt,con_g,con_eps
      use physpara,            only: effr_in
      use const,               only: RTYPE

      implicit none

!  ---  inputs:
      integer,  intent(in)    :: nmmiph,nx,nxjp(my),jjj,lev,ncld,kdt,me
!      integer,  intent(in)    :: ntcw,ntrw,ntiw,ntsw,ntgl,ntinc,ntrnc
      integer,  intent(in)    :: islimsk(nx,my_max)
      integer,  intent(in)    :: itimestep
      integer,  intent(in)    :: ivegtyp(nx,my_max)
      real,     intent(in)    :: tpi,dta,xlat(my),sdec
      real,     intent(in)    :: phii(nx,lev+1,my_max)
      real,     intent(in)    :: area(my_max)
      real,     intent(in)    :: ptop
      real,     intent(in)    :: rhc_mp(nx,lev,my_max)
      real,     intent(in)    :: snr(nx,my_max)
      real,     intent(in)    :: plt(nx,lev,my_max)
      real(kind=RTYPE), intent(in):: phi(nx,lev,my_max),pk(nx,lev,my_max), &
                                     sgeo(nx,my_max)
      real(kind=RTYPE), intent(in):: q0(nx,lev*ncld,my_max),dsigma(lev,2)
#ifdef Readaeroclx
      integer,  intent(in)    :: naero
      real(kind=RTYPE), intent(in):: aeroclx(nx,lev*naero,my_max)
#endif
!  ---  inputs/outputs:
      real, intent(inout) :: tt(nx,lev,my_max)
      real,     intent(inout) :: qa(nx,lev,my_max)
      real(kind=RTYPE), intent(inout) :: vvel(nx,lev,my_max) !mb/s
      real, intent(inout) :: ut(nx,lev,my_max),vt(nx,lev,my_max)
      real(kind=RTYPE), intent(inout) :: qt(nx,lev*ncld,my_max)
      real(kind=RTYPE), intent(in   ) :: pst(nx,my_max)
!  ---  outputs:
      real,     intent(inout)   :: re_cloud(nx,lev,my_max),re_ice(nx,lev,my_max),   &
                                   re_snow(nx,lev,my_max),re_rain(nx,lev,my_max)
      real,     intent(inout)   :: rlsp(nx,my_max),sr(nx,my_max)
      real,     intent(inout)   :: rlspi(nx,my_max),rlsps(nx,my_max), &
                                   rlspg(nx,my_max)
!  ---  local arrays:
      integer   kc,k,i,n,nxj,myim(my_max),jj,j,m
      real      prsl(nx,lev,my_max),del(nx,lev,my_max)
      real      ttc(nx,lev,my_max),qtc(nx,lev,my_max),qtr(nx,lev,my_max), &
                qtrw(nx,lev,my_max), qti(nx,lev,my_max),qtsw(nx,lev,my_max), &
                qtgl(nx,lev,my_max),ntnc(nx,lev,2,my_max),  &
                refl10(nx,lev,my_max)
      real      rainncv(nx,my_max),snowncv(nx,my_max),graupelncv(nx,my_max)
      real      xlat_myim(my_max)
      real      icem,tem
      logical   lradar
      logical   convert_dry_q,q_remove_cond
      real,dimension(:,:),allocatable ::                                  &
              land1d
      real,dimension(:,:,:),allocatable ::                                &
              u2d,v2d,p2d
      real, dimension(:,:,:), allocatable :: qni2d, qnr2d, &
              qv2d,qc2d,qr2d,qi2d,qs2d,qg2d,qnc2d, rew2d, rei2d, res2d,   &
              t2d, dz2d, cld2d, w2d, rer2d, reg2d, dp2d, rhc2d
! 2M Thompson MP
      real,dimension(:,:,:),allocatable ::                                &
              nwfa,nifa,pfils,pflls,vt_dbz_wt
      real,dimension(:,:),allocatable :: nwfasfc,nifasfc,rainnc,snownc,   &
              icenc,graupelnc,icencv,gridkm
      real,dimension(:,:,:),allocatable :: rand_pert
      real,dimension(:,:),allocatable :: spp_prt_list,spp_stddev_cutoff
      character(len=3),dimension(:,:),allocatable :: spp_var_list
      real    dt_inner
      real    rho
      logical sedi_semi,ext_diag,first_time_step,reset_dBZ,aero_ind_fdb,  &
              diagflag,convert_dry_rho
      integer do_radar_ref,rand_perturb_on,has_reqc,has_reqi,has_reqs,    &
              kme_stoch,istep,nsteps,errflg,decfl
      character errmsg
      integer, parameter :: n_var_spp=1
! GFDLMP
      real, parameter ::                                                  &
                rainmin=1.0e-10 !(mm)
      logical   hydrostatic,phys_hydrostatic,sedi_w
     !GFDL MP v2 & v3
      real, dimension(:,:), allocatable ::                                &
                gsize,hs,water1d,rain1d,snow1d,ice1d,graupel1d,snr1d
      real, dimension(:,:,:), allocatable ::                              &
                q_con,cappa,te,dp2d_ef
#ifdef EXT_DIAG
      real, dimension(:,:), allocatable ::                                &
                cond0,dep0,evap0,sub0
      real, dimension(:,:,:), allocatable ::                              &
                prefluxw,prefluxr,prefluxi, prefluxs,prefluxg
#endif
      logical   consv_te,last_step,do_inline_mp
     !GFDL MP v1
      integer, dimension(:,:), allocatable ::                             &
                mask1d
      real, dimension(:,:), allocatable ::                              &
                garea,land2d,rain2d,snow2d,ice2d,graupel2d
      real, dimension(:,:,:), allocatable ::                            &
                qv3d,qc3d,qr3d,qi3d,qs3d,qg3d,cld3d,qnc3d,t3d,w3d,u3d,  &
                v3d,dp3d,dz3d,qvten3d,qcten3d,qrten3d,qiten3d,qsten3d,  &
                qgten3d,cldten3d,uten3d,vten3d,tten3d,rho2d
! Goddard (GCE) MP
      integer nmgce3
      integer ihail,ICE2
      integer ntimes
      real    mp_time,dts
      real    rhowater,rhosnow,dx
      real,dimension(:,:),allocatable ::                                &
              ht,hail2d,sr2d
      real,dimension(:,:,:),allocatable ::                              &
              th3d,qh3d,rho3d,pii3d,p3d,z3d,rew3d,rer3d,rei3d,res3d,    &
              reg3d,reh3d,refl_10cm
      real,dimension(:,:,:,:),allocatable :: aero4d
#ifdef EXT_DIAG
      real,dimension(:,:,:),allocatable ::                              &
              physc, physe, physd, physs, physm, physf,                 &
              acphysc, acphyse, acphysd, acphyss, acphysm, acphysf,     &
              preci3d, precs3d, precg3d, precr3d, prech3d
#endif
      logical, intent(in) :: SL_sedi, sat_predict, new_saturation, &
              use_cpm, use_declination
      logical :: benchmark = .false.
      integer :: async_id = 1
      ! allocate temporary arrays (GCE) on GPU
      !$acc enter data create(prsl, del, ttc, qtc, qtr, qtrw, qti, qtsw, &
      !$acc&      qtgl, refl10, ntnc, rainncv, snowncv, graupelncv) async(async_id)
      !$acc wait(async_id)
      !$acc parallel loop private(j, nxj) async(async_id)
      do jj = 1, jlistnum
         j = jlist1(jj)
         nxj = nxjp(j)
         myim(jj) = nxjp(j)
         xlat_myim(jj) = xlat(j)
      end do
      
      convert_dry_q = .true.
      q_remove_cond = .false.
   !
   ! reset all value to zero
      !$acc parallel loop gang collapse(3) async(async_id)
      do jj = 1, jlistnum
         do k = 1, lev
            do i = 1, nx
               prsl(i,k,jj)  = 0.
               del(i,k,jj)   = 0.
               ttc(i,k,jj)   = 0.
               qtc(i,k,jj)   = 0.
               qtr(i,k,jj)   = 0.
               qtrw(i,k,jj)  = 0.
               qti(i,k,jj)   = 0.
               qtsw(i,k,jj)  = 0.
               qtgl(i,k,jj)  = 0.
               refl10(i,k,jj)= 0.
               ntnc(i,k,1,jj)  = 0.
               ntnc(i,k,2,jj)  = 0.
            end do
         end do
      end do
      !$acc parallel loop collapse(2) async(async_id)
      do jj = 1, jlistnum
         do i = 1, nx
            rainncv(i,jj)=0.
            snowncv(i,jj)=0.
            graupelncv(i,jj)=0.
         end do
      end do
   !
      lradar= .false.
      icem  =  4./3.*tpi*3.2768*1.e-14*890.
   !
      if ( nmmiph.eq.6 .or. nmmiph.eq.8 ) then
         do jj = 1, jlistnum
            do k=1,lev
               kc=lev-k+1
               do i=1,myim(jj)
                  prsl(i,kc,jj)= plt(i,k,jj)*100. ! change to Pa
                  del(i,kc,jj) = (dsigma(k,1)*pst(i,jj)+dsigma(k,2))*100. !change to Pa
                  qtc(i,kc,jj) = qt(i,             k,jj)
                  qtr(i,kc,jj) = qt(i,(ntcw-1)*lev+k,jj)
                  qtrw(i,kc,jj)= qt(i,(ntrw-1)*lev+k,jj)
                  qti(i,kc,jj) = qt(i,(ntiw-1)*lev+k,jj)
                  qtsw(i,kc,jj)= qt(i,(ntsw-1)*lev+k,jj)
                  qtgl(i,kc,jj)= qt(i,(ntgl-1)*lev+k,jj)
                  ttc(i,kc,jj) = tt(i,             k,jj)
               enddo
            enddo
         end do
   !
   !        WSM6
         if ( nmmiph .eq. 6 ) then 
            do jj = 1, jlistnum
               call wsm6(ttc(1,1,jj),phii(1,1,jj),qtc(1,1,jj),qtr(1,1,jj), &
                         qtrw(1,1,jj),qti(1,1,jj),qtsw(1,1,jj),qtgl(1,1,jj), &
                         prsl(1,1,jj),del(1,1,jj),dta, &
                         rainncv(1,jj),sr(1,jj),islimsk(1,jj),re_cloud(1,1,jj), &
                         re_ice(1,1,jj),re_snow(1,1,jj),       &
                         1,nx,1,lev,1,myim(jj),1,lev,snowncv(1,jj),graupelncv(1,jj))
            end do
         end if
   !        Thompson
         if ( nmmiph .eq. 8 )then
            allocate( qni2d(nx,lev,my_max),qnr2d(nx,lev,my_max) )
            do jj = 1, jlistnum
               do k=1,lev
                  kc=lev-k+1
                  do i=1,myim(jj)
                     qni2d(i,kc,jj) = qt(i,(ntinc-1)*lev+k,jj)                   &
                          +max(0.,qti(i,kc,jj)-q0(i,(ntiw-1)*lev+k,jj))/icem
                     qnr2d(i,kc,jj) = qt(i,(ntrnc-1)*lev+k,jj)
                  enddo
               enddo
            end do
            do jj = 1, jlistnum
               call thompson_driver(1,nx,1,lev,1,myim(jj),1,lev,qtc(1,1,jj), &
                       qtr(1,1,jj),qtrw(1,1,jj), qti(1,1,jj),qtsw(1,1,jj), &
                       qtgl(1,1,jj),qni2d(1,1,jj),qnr2d(1,1,jj),ttc(1,1,jj),                    &
                       prsl(1,1,jj),del(1,1,jj),dta,kdt,rainncv(1,jj),sr(1,jj), &
                       islimsk(1,jj),refl10(1,1,jj),lradar,re_cloud(1,1,jj), &
                       re_ice(1,1,jj),re_snow(1,1,jj),me,phii(1,1,jj))
            end do
            do jj = 1, jlistnum
               do k=1,lev
                  kc=lev-k+1
                  do i=1,myim(jj)
                     qt(i,(ntinc-1)*lev+k,jj)=qni2d(i,kc,jj)
                     qt(i,(ntrnc-1)*lev+k,jj)=qnr2d(i,kc,jj)
                  enddo
               enddo
            end do
         endif
   !
         do jj = 1, jlistnum
            do i=1,myim(jj)
               rlsp(i,jj) = rainncv(i,jj) + snowncv(i,jj) + graupelncv(i,jj)
               rlspi(i,jj) = 0.   !ice not precipitating
               rlsps(i,jj) = snowncv(i,jj)
               rlspg(i,jj) = graupelncv(i,jj)
            enddo
         end do
   !
         do jj = 1, jlistnum
            do k=1,lev
               kc=lev-k+1
               do i=1,myim(jj)
                  qt(i,             k,jj) = qtc(i,kc,jj)
                  qt(i,(ntcw-1)*lev+k,jj) = qtr(i,kc,jj)
                  qt(i,(ntrw-1)*lev+k,jj) = qtrw(i,kc,jj)
                  qt(i,(ntiw-1)*lev+k,jj) = qti(i,kc,jj)
                  qt(i,(ntsw-1)*lev+k,jj) = qtsw(i,kc,jj)
                  qt(i,(ntgl-1)*lev+k,jj) = qtgl(i,kc,jj)
                  tt(i,             k,jj) = ttc(i,kc,jj)
                  
               enddo
            enddo
         end do
         if ( nmmiph .eq. 8 ) deallocate( qni2d,qnr2d )

      endif

   !     2M Thompson
      if ( nmmiph .eq. 18 ) then
         allocate                                                             &
          ( t2d(nx,lev,my_max),qv2d(nx,lev,my_max),qc2d(nx,lev,my_max),       &
            qr2d(nx,lev,my_max),qi2d(nx,lev,my_max),qs2d(nx,lev,my_max),      &
            qg2d(nx,lev,my_max),qni2d(nx,lev,my_max),qnr2d(nx,lev,my_max),    &
            qnc2d(nx,lev,my_max),rew2d(nx,lev,my_max),rei2d(nx,lev,my_max),   &
            res2d(nx,lev,my_max),nwfa(nx,lev,my_max),nifa(nx,lev,my_max),     &
            dz2d(nx,lev,my_max),w2d(nx,lev,my_max),pfils(nx,lev,my_max),      &
            pflls(nx,lev,my_max),vt_dbz_wt(nx,lev,my_max),nwfasfc(nx,my_max), &
            nifasfc(nx,my_max),rainnc(nx,my_max),snownc(nx,my_max),           &
            icenc(nx,my_max),graupelnc(nx,my_max),icencv(nx,my_max),          &
            rand_pert(nx,1,my_max),spp_prt_list(n_var_spp,my_max),            &
            spp_stddev_cutoff(n_var_spp,my_max),spp_var_list(n_var_spp,my_max) )
         if ( cfflag_thom .eq. 2 ) allocate                              &
          ( cld2d(nx,lev,my_max),land1d(nx,my_max),gridkm(nx,my_max) )
         
         sedi_semi=.false.        !use Semi-Lagrangian sedimentation for rain and graupel
         ext_diag=.false.         !extended diagnostics, array pointers only associated if ext_diag is .true.
         first_time_step=.false.  !(not sure)
         reset_dBZ=.false.        !if true, set melti=.true.
         aero_ind_fdb=.false.     !(not sure)
         diagflag=.false.         !if diagflag=true and do_radar_ref=1, call calc_refl10cm
         convert_dry_rho=.false.
         do_radar_ref=0
         rand_perturb_on=0        !if!=0, use SPP
         if ( effr_in ) then
            has_reqc=1             !calculate effective radii of cloud water
            has_reqi=1             !calculate effective radii of cloud ice
            has_reqs=1             !calculate effective radii of snow
         endif
         
         dt_inner=dta     !inner time step
         decfl=1          !if .not. sedi_semi
         kme_stoch=1
         istep=1          !current step
         nsteps=1         !maximum number of steps
         
         do jj = 1, jlistnum
            do k = 1, lev
               do i = 1, nx
                  t2d(i,k,jj)=0.
                  qv2d(i,k,jj)=0.
                  qc2d(i,k,jj)=0.
                  qr2d(i,k,jj)=0.
                  qi2d(i,k,jj)=0.
                  qs2d(i,k,jj)=0.
                  qg2d(i,k,jj)=0.
                  rew2d(i,k,jj)=0.
                  rei2d(i,k,jj)=0.
                  res2d(i,k,jj)=0.
                  qni2d(i,k,jj)=0.         !number concentracion of ice
                  qnr2d(i,k,jj)=0.         !number concentracion of rain
                  qnc2d(i,k,jj)=0.         !number concentracion of cloud droplet
                  nwfa(i,k,jj)=0.          !number concentration of water friendly aerosol
                  nifa(i,k,jj)=0.          !number concentration of ice friendly aerosol
                  w2d(i,k,jj)=0.
                  dz2d(i,k,jj)=0.
                  pfils(i,k,jj)=0.
                  pflls(i,k,jj)=0.
                  vt_dbz_wt(i,k,jj)=0.
               end do
            end do
            do i = 1, nx
               nwfasfc(i,jj)=0.       !nwfa at surface
               nifasfc(i,jj)=0.       !nifa at surface
               rainnc(i,jj)=0.        !number concentracion of precipitating rain
               snownc(i,jj)=0.        !number concentracion of precipitating snow
               icenc(i,jj)=0.         !number concentracion of precipitating ice
               graupelnc(i,jj)=0.     !number concentracion of precipitating graupel
               icencv(i,jj)=0.        !amount of precipitating ice
               rand_pert(i,1,jj)=0.
            end do
            do i = 1, n_var_spp
               spp_stddev_cutoff(i, jj)=0.
            end do
         end do
            
         do jj = 1, jlistnum
            do k = 1, lev
               kc = lev - k + 1
               do i = 1, myim(jj)
                  qv2d (i,k,jj) = qt(i,              kc,jj)
                  qc2d (i,k,jj) = qt(i,(ntcw-1) *lev+kc,jj)
                  qr2d (i,k,jj) = qt(i,(ntrw-1) *lev+kc,jj)
                  qi2d (i,k,jj) = qt(i,(ntiw-1) *lev+kc,jj)
                  qs2d (i,k,jj) = qt(i,(ntsw-1) *lev+kc,jj)
                  qg2d (i,k,jj) = qt(i,(ntgl-1) *lev+kc,jj)
                  qni2d(i,k,jj) = qt(i,(ntinc-1)*lev+kc,jj)
                  qnr2d(i,k,jj) = qt(i,(ntrnc-1)*lev+kc,jj)
                  prsl (i,k,jj) = plt(i,kc,jj)*100.                         !layer pressure (Pa)
                  t2d  (i,k,jj) = tt(i,kc,jj)
                  dz2d (i,k,jj) = (phii(i,k+1,jj)-phii(i,k,jj))/con_g          !layer width (m)
                  
                  !> - Convert specific humidity to water vapor mixing ratio.
                  !> - Also, hydrometeor variables are mass or number mixing ratio
                  !> - either kg of species per kg of dry air, or per kg of (dry + vapor).
                  if ( convert_dry_q ) then
                     if ( q_remove_cond ) then
                        ! q convert to : hydrometeor mass / ( dry mass )
                        tem = qv2d(i,k,jj) + qc2d(i,k,jj) + qr2d(i,k,jj) + qi2d(i,k,jj) +   &
                              qs2d(i,k,jj) + qg2d(i,k,jj)
                     else
                        ! q convert to : hydrometeor mass / ( dry mass + condensates )
                        tem = qv2d(i,k,jj)
                     endif
                     qv2d (i,k,jj) = qv2d (i,k,jj)/(1. - tem)
                     qc2d (i,k,jj) = qc2d (i,k,jj)/(1. - tem)
                     qr2d (i,k,jj) = qr2d (i,k,jj)/(1. - tem)
                     qi2d (i,k,jj) = qi2d (i,k,jj)/(1. - tem)
                     qs2d (i,k,jj) = qs2d (i,k,jj)/(1. - tem)
                     qg2d (i,k,jj) = qg2d (i,k,jj)/(1. - tem)
                     qni2d(i,k,jj) = qni2d(i,k,jj)/(1. - tem)
                     qnr2d(i,k,jj) = qnr2d(i,k,jj)/(1. - tem)
                  endif
                  
                  ! Ensure non-negative mass mixing ratios of all water variables
                  if ( qv2d(i,k,jj) .lt. 0. ) qv2d(i,k,jj) = 1.E-10  !should *never* be identically zero
                  if ( qc2d(i,k,jj) .lt. 0. ) qc2d(i,k,jj) = 0.
                  if ( qr2d(i,k,jj) .lt. 0. ) qr2d(i,k,jj) = 0.
                  if ( qi2d(i,k,jj) .lt. 0. ) qi2d(i,k,jj) = 0.
                  if ( qs2d(i,k,jj) .lt. 0. ) qs2d(i,k,jj) = 0.
                  if ( qg2d(i,k,jj) .lt. 0. ) qg2d(i,k,jj) = 0.
                  
                  if ( convert_dry_q ) then
                     ! dry air density : rho = 0.622*p/(Rd*T*(0.622+qv))
                     rho = con_eps*prsl(i,k,jj)/                               &
                           (con_rd*t2d(i,k,jj)*(con_eps+qv2d(i,k,jj)))
                  else
                     ! moist air density : rho = p/(Rd*T*(1+0.608*qv))
                     rho = prsl(i,k,jj)/(con_rd*t2d(i,k,jj)*(1+con_fvirt*qv2d(i,k,jj)))
                  endif
                  w2d(i,k,jj) = - vvel(i,kc,jj)*100./(rho*con_g)          !vertical velocity (m/s)
             
                  ! 1st guess number concentration where mass non-zero
   !               if ( kdt .eq. 1 ) then
   !                 if ( qi2d(i,k) .gt. 0. )                                  &
   !                   qni2d(i,k) = make_IceNumber(qi2d(i,k)*rho,t2d(i,k))/rho
   !                 if ( qr2d(i,k) .gt. 0. )                                  &
   !                   qnr2d(i,k) = make_RainNumber(qr2d(i,k)*rho,t2d(i,k))/rho
   !               endif
               enddo
            enddo
         end do
         
         do jj = 1, jlistnum
            call new_thompson_driver                                            &
                       ( qv2d(1,1,jj),qc2d(1,1,jj),qr2d(1,1,jj),qi2d(1,1,jj),   &
                         qs2d(1,1,jj),qg2d(1,1,jj),qni2d(1,1,jj),qnr2d(1,1,jj), &
                         qnc2d(1,1,jj),nwfa(1,1,jj),nifa(1,1,jj),nwfasfc(1,jj), &
                         nifasfc(1,jj),                                         &!optional
                         t2d(1,1,jj),&!th,pii,                                  &!optional ??
                         prsl(1,1,jj),w2d(1,1,jj),dz2d(1,1,jj),dta,dt_inner,    &
                         sedi_semi,decfl,islimsk(i, jj),                        &
                         rainnc(1,jj),rainncv(1,jj),                            &
                         snownc(1,jj),snowncv(1,jj),icenc(1,jj),icencv(1,jj),   &
                         graupelnc(1,jj),graupelncv(1,jj),                      &!optional
                         sr(i, jj),                                             &
                         refl10(1,1,jj),diagflag,do_radar_ref,                  &!optional
                         vt_dbz_wt(1,1,jj),                                     &!optional
                         first_time_step,                                       &
                         rew2d(1,1,jj),rei2d(1,1,jj),res2d(1,1,jj),             &!optional
                         has_reqc,has_reqi,has_reqs,                            &
                         aero_ind_fdb,                                          &!optional
                         rand_perturb_on,                                       &
                         kme_stoch,                                             &
                         rand_pert(:,:,jj),spp_prt_list(:,jj),spp_var_list(:,jj),&
                         spp_stddev_cutoff(:,jj),n_var_spp,                     &
                         1,nx,1,lev,                                            &
                         1,myim(jj),1,lev,                                      &
                         reset_dBZ,istep,nsteps,                                &
                         errmsg,errflg,                                         &!optional
                         ext_diag,                                              &
   #ifdef EXT_DIAG
                         prw_vcdc, prw_vcde, tpri_inu, tpri_ide_d,              &
                         tpri_ide_s, tprs_ide, tprs_sde_d,                      &
                         tprs_sde_s, tprg_gde_d,                                &
                         tprg_gde_s, tpri_iha, tpri_wfz,                        &
                         tpri_rfz, tprg_rfz, tprs_scw, tprg_scw,                &
                         tprg_rcs, tprs_rcs, tprr_rci, tprg_rcg,                &
                         tprw_vcd_c, tprw_vcd_e, tprr_sml,                      &
                         tprr_gml, tprr_rcg,                                    &
                         tprr_rcs, tprv_rev, tten3, qvten3,                     &
                         qrten3, qsten3, qgten3, qiten3, niten3,                &
                         nrten3, ncten3, qcten3,                                &
   #endif   
                         pfils(1,1,jj), pflls(1,1,jj) )
         end do
      
         if ( cfflag_thom .eq. 2 ) then
            do jj = 1, jlistnum
               do i = 1, nx
                  land1d(i,jj)=0.
                  gridkm(i,jj)=0.
               end do
               do k = 1, lev
                  do i = 1, nx
                     cld2d(i,k,jj)=0.
                  end do
               end do
            end do
            do jj = 1, jlistnum
               do i = 1, myim(jj)
                  if ( islimsk(i,jj) .eq. 1 ) then
                     land1d(i,jj) = 1.      !land
                  else
                     land1d(i,jj) = 2.      !ocean or seaice
                  endif
                  gridkm(i,jj) = sqrt(area(jj))/1000.                !grid length (km)
               enddo
            end do
            
            do jj = 1, jlistnum
               call cal_cldfra3                                                   &
                        ( 1, nx, 1, myim(jj), 1, lev, 1, lev,                     &
                          cld2d(1,1,jj), qv2d(1,1,jj), qc2d(1,1,jj),              &
                          qi2d(1,1,jj), qs2d(1,1,jj), dz2d(1,1,jj),               &
                          prsl(1,1,jj), t2d(1,1,jj), land1d(1,jj), gridkm(1,jj),  &
                          .false., 1.5, .false. )
            end do
            do jj = 1, jlistnum
               do k = 1, lev
                  kc = lev - k + 1
                  do i = 1, myim(jj)
                     qa(i,k,jj) = cld2d(i,kc,jj)
                  enddo
               enddo
            end do
         endif
         
         do jj = 1, jlistnum
            do k = 1, lev
               kc = lev - k + 1
               do i = 1, myim(jj)
                  if ( convert_dry_q ) then
                     if ( q_remove_cond ) then
                        tem = qv2d(i,kc,jj) + qc2d(i,kc,jj) + qr2d(i,kc,jj) +  &
                              qi2d(i,kc,jj) + qs2d(i,kc,jj) + qg2d(i,kc,jj)
                     else
                        tem = qv2d(i,kc,jj)
                     endif
                     qv2d (i,kc,jj) = qv2d (i,kc,jj)/(1. + tem)
                     qc2d (i,kc,jj) = qc2d (i,kc,jj)/(1. + tem)
                     qr2d (i,kc,jj) = qr2d (i,kc,jj)/(1. + tem)
                     qi2d (i,kc,jj) = qi2d (i,kc,jj)/(1. + tem)
                     qs2d (i,kc,jj) = qs2d (i,kc,jj)/(1. + tem)
                     qg2d (i,kc,jj) = qg2d (i,kc,jj)/(1. + tem)
                     qni2d(i,kc,jj) = qni2d(i,kc,jj)/(1. + tem)
                     qnr2d(i,kc,jj) = qnr2d(i,kc,jj)/(1. + tem)
                  endif
                  
                  qt(i,              k,jj) = max( qv2d (i,kc,jj) , 0. )
                  qt(i,(ntcw-1) *lev+k,jj) = max( qc2d (i,kc,jj) , 0. )
                  qt(i,(ntrw-1) *lev+k,jj) = max( qr2d (i,kc,jj) , 0. )
                  qt(i,(ntiw-1) *lev+k,jj) = max( qi2d (i,kc,jj) , 0. )
                  qt(i,(ntsw-1) *lev+k,jj) = max( qs2d (i,kc,jj) , 0. )
                  qt(i,(ntgl-1) *lev+k,jj) = max( qg2d (i,kc,jj) , 0. )
                  qt(i,(ntinc-1)*lev+k,jj) = max( qni2d(i,kc,jj) , 0. )
                  qt(i,(ntrnc-1)*lev+k,jj) = max( qnr2d(i,kc,jj) , 0. )
                  tt(i,k,jj) = t2d(i,kc,jj)
                  
                  re_cloud(i,k,jj) = rew2d(i,k,jj)*1.E+6   ! m to micron
                  re_ice  (i,k,jj) = rei2d(i,k,jj)*1.E+6   ! m to micron
                  re_snow (i,k,jj) = res2d(i,k,jj)*1.E+6   ! m to micron
               enddo
            enddo
         end do
         do jj = 1, jlistnum
            do i = 1, myim(jj)
               rlsp(i,jj) = rainncv(i,jj)  !total large scale precipitation (kg/m^2=mm)
               rlspi(i,jj) = icencv(i,jj)      !not sure
               rlsps(i,jj) = snowncv(i,jj)     !not sure
               rlspg(i,jj) = graupelncv(i,jj)  !not sure
            enddo
         end do
            
            deallocate                                                      &
             ( t2d,qv2d,qc2d,qr2d,qi2d,qs2d,qg2d,qni2d,qnr2d,qnc2d,         &
               nwfa,nifa,dz2d,w2d,pfils,pflls,vt_dbz_wt,nwfasfc,nifasfc,    &
               rew2d,rei2d,res2d,rainnc,snownc,icenc,graupelnc,icencv,      &
               rand_pert,spp_prt_list,spp_stddev_cutoff,spp_var_list )
            if ( cfflag_thom .eq. 2 ) deallocate ( cld2d,land1d,gridkm )
         endif  !end of if nmmiph=18

   !     GFDL MP v1
      do jj = 1, jlistnum
         if ( nmmiph .eq. 11 ) then
            allocate                                                            &
             ( qv3d(myim(jj),1,lev),qc3d(myim(jj),1,lev),qr3d(myim(jj),1,lev),  &
               qi3d(myim(jj),1,lev),qs3d(myim(jj),1,lev),qg3d(myim(jj),1,lev),  &
               cld3d(myim(jj),1,lev),qnc3d(myim(jj),1,lev),t3d(myim(jj),1,lev), &
               w3d(myim(jj),1,lev),u3d(myim(jj),1,lev),v3d(myim(jj),1,lev),     &
               dp3d(myim(jj),1,lev),dz3d(myim(jj),1,lev),                       &
               qvten3d(myim(jj),1,lev),qcten3d(myim(jj),1,lev),                 &
               qrten3d(myim(jj),1,lev), qiten3d(myim(jj),1,lev),                &
               qsten3d(myim(jj),1,lev),qgten3d(myim(jj),1,lev),                 &
               cldten3d(myim(jj),1,lev),uten3d(myim(jj),1,lev),                 &
               vten3d(myim(jj),1,lev), tten3d(myim(jj),1,lev),                  &
               rew2d(nx,lev,my_max),rei2d(nx,lev,my_max),rer2d(nx,lev,my_max),  &
               res2d(nx,lev,my_max), reg2d(nx,lev,my_max),land2d(myim(jj),1),   &
               rain2d(myim(jj),1),snow2d(myim(jj),1), ice2d(myim(jj),1),        &
               graupel2d(myim(jj),1),garea(myim(jj),1),rhc2d(nx,lev,my_max) )
            if ( effr_in ) allocate                                             &
               ( dp2d(nx,lev,my_max),rho2d(nx,lev,my_max),qc2d(nx,lev,my_max),  &
                 qr2d(nx,lev,my_max),qi2d(nx,lev,my_max),qs2d(nx,lev,my_max),   &
                 qg2d(nx,lev,my_max),mask1d(nx,my_max), t2d(nx,lev,my_max) )
            land2d = 0.
            garea = 0.
            qvten3d = 0.
            qcten3d = 0.
            qrten3d = 0.
            qiten3d = 0.
            qsten3d = 0.
            qgten3d = 0.
            cldten3d = 0.
            tten3d = 0.
            uten3d = 0.
            vten3d = 0.
            rain2d = 0.
            snow2d = 0.
            ice2d  = 0.
            graupel2d = 0.
            
            hydrostatic = .false.       !flag for hydrostatic solver
            phys_hydrostatic = .true.   !flag for hydrostatic heating from physics 
            sedi_w = .false.
            
            do i = 1, myim(jj)
               if( islimsk(i,jj) .eq. 1 ) land2d(i,1) = 1.  !land fraction
               if( effr_in ) mask1d(i,jj) = islimsk(i,jj)      !land-sea mask
               garea(i,1) = area(jj)                         !area of grid box (m^-2)
            enddo
            
            do k = 1, lev
               kc = lev - k + 1
               do i = 1, myim(jj)
                  if ( sedi_w ) then
                     prsl(i,k,jj) = 100.0 * plt(i,k,jj)             !layer mean pressure (from mb to Pa)
                     w3d(i,1,k) = -vvel(i,k,jj)*(1.+con_fvirt*qt(i,k,jj))*tt(i,k,jj)    &
                                /prsl(i,k,jj)*con_rd/con_g      !vertical velocity (m/s)
                  else
                     w3d(i,1,k) = 0.
                  endif
                  
                  qv3d (i,1,k) = qt(i,             k,jj)
                  qc3d (i,1,k) = qt(i,(ntcw-1)*lev+k,jj)
                  qr3d (i,1,k) = qt(i,(ntrw-1)*lev+k,jj)
                  qi3d (i,1,k) = qt(i,(ntiw-1)*lev+k,jj)
                  qs3d (i,1,k) = qt(i,(ntsw-1)*lev+k,jj)
                  qg3d (i,1,k) = qt(i,(ntgl-1)*lev+k,jj)
                  qnc3d(i,1,k) = 0.                      ! =0. for prog_ccn=.false. (cm^-3)
                  cld3d(i,1,k) = 0.                      !layer cloud fraction (=0. for do_qa=.false.)
                  t3d  (i,1,k) = tt(i,k,jj)                 !temperature
                  u3d  (i,1,k) = ut(i,k,jj)                 !zonal wind (m/s)
                  v3d  (i,1,k) = vt(i,k,jj)                 !meridional wind (m/s)
                  dp3d (i,1,k) = (dsigma(k,1)*pst(i,jj)+dsigma(k,2))*100. !difference of interface pressure (Pa)
                  dz3d (i,1,k) = (phii(i,kc,jj)-phii(i,kc+1,jj))/con_g       !differences of height (m), dz<0
                  rhc2d(i,  k,jj) = rhc_mp(i,k,jj)
                  if ( effr_in ) dp2d(i,k,jj) = dp3d(i,1,k)
               enddo
            enddo
            
            call gfdl_cloud_microphys_driver                                &
                    ( qv3d, qc3d, qr3d, qi3d, qs3d, qg3d, cld3d, qnc3d,     &
                      qvten3d, qcten3d, qrten3d, qiten3d, qsten3d, qgten3d, &
                      cldten3d, tten3d, t3d, w3d, u3d, v3d, uten3d, vten3d, &
                      dz3d, dp3d, garea, dta, land2d,                       &
                      rain2d, snow2d, ice2d, graupel2d,                     &
                      rhc2d(:,:,jj), hydrostatic, phys_hydrostatic,         &
                      1, myim(jj), 1, 1, 1, lev, 1, lev )
            
            do k = 1, lev
               do i = 1, myim(jj)
                  qc2d(i,k,jj) = qc3d(i,1,k) + qcten3d(i,1,k) * dta
                  qr2d(i,k,jj) = qr3d(i,1,k) + qrten3d(i,1,k) * dta
                  qi2d(i,k,jj) = qi3d(i,1,k) + qiten3d(i,1,k) * dta
                  qs2d(i,k,jj) = qs3d(i,1,k) + qsten3d(i,1,k) * dta
                  qg2d(i,k,jj) = qg3d(i,1,k) + qgten3d(i,1,k) * dta
                  
                  qt(i,             k,jj) = qv3d(i,1,k) + qvten3d(i,1,k) * dta
                  qt(i,(ntcw-1)*lev+k,jj) = qc2d(i,k,jj)
                  qt(i,(ntrw-1)*lev+k,jj) = qr2d(i,k,jj)
                  qt(i,(ntiw-1)*lev+k,jj) = qi2d(i,k,jj)
                  qt(i,(ntsw-1)*lev+k,jj) = qs2d(i,k,jj)
                  qt(i,(ntgl-1)*lev+k,jj) = qg2d(i,k,jj)
                  qa(i,k,jj)  = cld3d(i,1,k) + cldten3d(i,1,k) * dta
                  tt(i,k,jj)  = t3d(i,1,k) + tten3d(i,1,k) * dta
                  ut(i,k,jj)  = u3d(i,1,k) + uten3d(i,1,k) * dta
                  vt(i,k,jj)  = v3d(i,1,k) + vten3d(i,1,k) * dta
                  
                  if ( sedi_w ) then
                     vvel(i,k,jj) = -w3d(i,1,k)*prsl(i,k,jj)*con_g/con_rd       &
                                 /((1+con_fvirt*qt(i,k,jj))*tt(i,k,jj))
                  endif
                  
                  if ( effr_in ) then
                     rho2d(i,k,jj) = 0.622*prsl(i,k,jj)                         &
                                 /( con_rd*tt(i,k,jj)*(qt(i,k,jj)+0.622) ) !air density (kg/m^3)
                  endif
               enddo
            enddo
            
            if ( effr_in ) then
               do k = 1, lev
                  do i = 1, myim(jj)
                     t2d(i,k,jj) = tt(i,k,jj)
                  enddo
               enddo
               call cloud_diagnosis                                           &
                    ( 1, myim(jj), 1, lev, rho2d(1,1,jj), dp2d(1,1,jj),       &
                      mask1d(1,jj), qc2d(1,1,jj), qi2d(1,1,jj), qr2d(1,1,jj), &
                      qs2d(1,1,jj), qg2d(1,1,jj), t2d(1,1,jj), rew2d(1,1,jj), &
                      rei2d(1,1,jj), rer2d(1,1,jj), res2d(1,1,jj), reg2d(1,1,jj) )
               do k = 1, lev
                  kc = lev - k + 1
                  do i = 1, myim(jj)
                     re_cloud(i,k,jj)   = rew2d(i,kc,jj)   !(micron)
                     re_ice(i,k,jj)     = rei2d(i,kc,jj)   !(micron)
                     re_rain(i,k,jj)    = rer2d(i,kc,jj)   !(micron)
                     re_snow(i,k,jj)    = res2d(i,kc,jj)   !(micron)
                  enddo
               enddo
            endif
   !        
            do i = 1, myim(jj)
               if ( rain2d(i,1)    .lt. rainmin ) rain2d(i,1)    = 0.0
               if ( ice2d(i,1)     .lt. rainmin ) ice2d(i,1)     = 0.0
               if ( snow2d(i,1)    .lt. rainmin ) snow2d(i,1)    = 0.0
               if ( graupel2d(i,1) .lt. rainmin ) graupel2d(i,1) = 0.0
               
               rlsp(i,jj) = rain2d(i,1)+snow2d(i,1)+ice2d(i,1)+graupel2d(i,1)  !total large scale precipitation (mm)
               rlspi(i,jj) = ice2d(i,1)
               rlsps(i,jj) = snow2d(i,1)
               rlspg(i,jj) = graupel2d(i,1)
               if ( rlsp(i,jj) .gt. rainmin ) then
                  sr(i,jj) = (snow2d(i,1)+ice2d(i,1)+graupel2d(i,1))             &
                         /(rain2d(i,1)+snow2d(i,1)+ice2d(i,1)+graupel2d(i,1))  !snow ratio
               else
                  sr(i,jj) = 0.0
               endif
            enddo
            
            deallocate                                                      &
             ( qv3d,qc3d,qr3d,qi3d,qs3d,qg3d,cld3d,qnc3d,t3d,w3d,u3d,v3d,   &
               dp3d,dz3d,qvten3d,qcten3d,qrten3d,qiten3d,qsten3d,qgten3d,   &
               cldten3d,uten3d,vten3d,tten3d,rew2d,rei2d,rer2d,res2d,reg2d, &
               land2d,rain2d,snow2d,ice2d,graupel2d,garea,rhc2d )
            if ( effr_in ) deallocate                                       &
               ( dp2d,rho2d,qc2d,qr2d,qi2d,qs2d,qg2d,mask1d,t2d )
         endif  ! end of nmmiph.eq.11
      end do

   !     GFDL MP v2 & v3
      if ( nmmiph .eq. 12 .or. nmmiph .eq. 13 ) then
         do jj = 1, jlistnum
            nxj = myim(jj)
            allocate                                                               &
             ( qv2d(nxj,lev,1),qc2d(nxj,lev,1),qr2d(nxj,lev,1),        &
               qi2d(nxj,lev,1),qs2d(nxj,lev,1),qg2d(nxj,lev,1),        &
               cld2d(nxj,lev,1),qnc2d(nxj,lev,1),qni2d(nxj,lev,1),     &
               w2d(nxj,lev,1),t2d(nxj,lev,1),u2d(nxj,lev,1),           &
               v2d(nxj,lev,1),dp2d(nxj,lev,1),dz2d(nxj,lev,1),         &
               hs(nxj,1),land1d(nxj,1),gsize(nxj,1),rain1d(nxj,1), &
               snow1d(nxj,1),ice1d(nxj,1),graupel1d(nxj,1),            &
               water1d(nxj,1),rhc2d(nxj,lev,1) )
            if ( effr_in ) allocate                                                &
             ( rew2d(nxj,lev,1),rei2d(nxj,lev,1),rer2d(nxj,lev,1),     &
               res2d(nxj,lev,1),reg2d(nxj,lev,1),snr1d(nxj,1),         &
               p2d(nxj,lev,1),dp2d_ef(nxj,lev,1) )
   #ifdef EXT_DIAG
            allocate                                                               &
             ( prefluxr(nxj,lev,1),prefluxi(nxj,lev,1),prefluxs(nxj,lev,1), &
               prefluxg(nxj,lev,1),prefluxw(nxj,lev,1),                    &
               q_con(nxj,lev,1),cappa(nxj,lev,1),te(nxj,lev,1) )
            allocate                                                               &
             ( cond0(nxj,1),dep0(nxj,1),evap0(nxj,1),sub0(nxj,1) )
   #endif   
            
            hydrostatic = .true.       !flag for hydrostatic solver (isobaric assumption)
            phys_hydrostatic = .true.   !flag for hydrostatic heating from physics 
            consv_te = .false.          !flag for energy conservation
            last_step = .true.          !flag for final clean-up (not sure)
            do_inline_mp = .false.      !flag for inline GFDLMP (so far must be false) :
                                        !  = false : run in physical driver, and use real temperature
                                        !  = true  : run in dynamical core, and use virtual temperature
            if(nmmiph.eq.12) sedi_w = sedi_w_v2  !flag for w momentum transportation during sedimentation
            if(nmmiph.eq.13) sedi_w = sedi_w_v3  !flag for w momentum transportation during sedimentation
            
            do i = 1, myim(jj)
               gsize(i,1) = 0.0
               land1d(i,1) = 0.0
               rain1d(i,1) = 0.
               snow1d(i,1) = 0.
               ice1d(i,1)  = 0.
               graupel1d(i,1) = 0.
               water1d(i,1) = 0.
   #ifdef EXT_DIAG
               cond0(i,1) = 0.0
               dep0(i,1)  = 0.0
               evap0(i,1) = 0.0
               sub0(i,1)  = 0.0
   #endif   
            end do
   #ifdef EXT_DIAG
            do k = 1, lev
               do i = 1, myim(jj)
                  te(i,k,1)    = 0.0
                  q_con(i,k,1) = 0.0  !not sure
                  cappa(i,k,1) = 0.0  !not sure
                  prefluxw(i,k,1)  = 0.0
                  prefluxr(i,k,1)  = 0.0
                  prefluxi(i,k,1)  = 0.0
                  prefluxs(i,k,1)  = 0.0
                  prefluxg(i,k,1)  = 0.0
               end do
            end do
   #endif   
            do i = 1, myim(jj)
               gsize(i,1) = sqrt(area(jj))     !square root of grid area (m)
               hs(i,1)    = sgeo(i,jj)        !terrain geopotential (m^2 s^-2)
               if( islimsk(i,jj) .eq. 1 ) land1d(i,1) = 1.  !land fraction
            enddo
            do k = 1, lev
               kc = lev - k + 1
               do i = 1, myim(jj)
                  qv2d (i,k,1) = qt(i,             k,jj)
                  qc2d (i,k,1) = qt(i,(ntcw-1)*lev+k,jj)
                  qr2d (i,k,1) = qt(i,(ntrw-1)*lev+k,jj)
                  qi2d (i,k,1) = qt(i,(ntiw-1)*lev+k,jj)
                  qs2d (i,k,1) = qt(i,(ntsw-1)*lev+k,jj)
                  qg2d (i,k,1) = qt(i,(ntgl-1)*lev+k,jj)
                  qnc2d(i,k,1) = 0.
                  qni2d(i,k,1) = 0.
                  cld2d(i,k,1) = 0.                      !layer cloud fracion, =0 for do_qa=.false.
                  t2d  (i,k,1) = tt(i,k,jj)                 !temperature
                  u2d  (i,k,1) = ut(i,k,jj)                 !zonal wind (m/s)
                  v2d  (i,k,1) = vt(i,k,jj)                 !meridional wind (m/s)
                  dp2d (i,k,1) = (dsigma(k,1)*pst(i,jj)+dsigma(k,2))*100. !difference of interface pressure (Pa)
                  dz2d (i,k,1) = (phii(i,kc,jj)-phii(i,kc+1,jj))/con_g       !differences of height (m), dz<0
                  rhc2d(i,k,1) = rhc_mp(i,k,jj)
                  
                  if ( sedi_w ) then
                     prsl(i,k,jj) = 100.0 * plt(i,k,jj)            !layer mean pressure (from mb to Pa)
                     ! use moist air density : rho = p/(Rd*T*(1+0.608*qv))
                     rho = prsl(i,k,jj)/(con_rd*t2d(i,k,1)*(1+con_fvirt*qv2d(i,k,1)))
                     w2d(i,k,1) = - vvel(i,k,jj)*100./(rho*con_g) !vertical velocity (m/s)
                  else
                     w2d(i,k,1)  = 0.
                  endif
               enddo
            enddo
            
            ! GFDL MP v2
            if ( nmmiph .eq. 12 ) then
               call gfdlv2_driver                                               &
                     ( qv2d, qc2d, qr2d, qi2d,  &
                       qs2d, qg2d,                              &
                       cld2d, qnc2d,qni2d,              &
                       t2d, w2d, u2d, v2d,      &
                       dz2d, dp2d, gsize, dta, hs,  &
                       land1d, rain1d, snow1d, ice1d,   &
                       graupel1d, hydrostatic,                            &
                       1, myim(jj), 1, lev, consv_te,                           &
   #ifdef EXT_DIAG
                       q_con, cappa, te,                &
                       prefluxr, prefluxi, prefluxs,    &
                       prefluxg,                                        &
                       cond0, dep0, evap0, sub0,        &
   #endif      
                       rhc2d, last_step, do_inline_mp )
            end if
            
            ! GFDL MP v3
            if ( nmmiph .eq. 13 ) then
               call gfdlv3_driver                                               &
                     ( qv2d, qc2d, qr2d, qi2d,  &
                       qs2d, qg2d,                              &
                       cld2d, qnc2d, qni2d,             &
                       t2d, w2d, u2d, v2d,      &
                       dz2d, dp2d, gsize, dta, hs,  &
                       land1d, water1d, rain1d, snow1d, &
                       ice1d, graupel1d, hydrostatic,               &
                       1, myim(jj), 1, lev, consv_te,                           &
   #ifdef EXT_DIAG
                       q_con, cappa, te,                &
                       prefluxw, prefluxr, prefluxi,    &
                       prefluxs, prefluxg,                      &
                       cond0, dep0, evap0, sub0,        &
   #endif      
                       last_step, do_inline_mp )
            end if
            
            do k = 1, lev
               do i = 1, myim(jj)
                  qt(i,             k,jj) = qv2d(i,k,1)
                  qt(i,(ntcw-1)*lev+k,jj) = qc2d(i,k,1)
                  qt(i,(ntrw-1)*lev+k,jj) = qr2d(i,k,1)
                  qt(i,(ntiw-1)*lev+k,jj) = qi2d(i,k,1)
                  qt(i,(ntsw-1)*lev+k,jj) = qs2d(i,k,1)
                  qt(i,(ntgl-1)*lev+k,jj) = qg2d(i,k,1)
                  qa(i,k,jj)  = cld2d(i,k,1)
                  tt(i,k,jj)  = t2d  (i,k,1)
                  ut(i,k,jj)  = u2d  (i,k,1)
                  vt(i,k,jj)  = v2d  (i,k,1)
                  
                  if ( sedi_w ) then
                     ! use moist air density : rho = p/(Rd*T*(1+0.608*qv))
                     rho = prsl(i,k,jj)/(con_rd*tt(i,k,jj)*(1+con_fvirt*qt(i,k,jj)))
                     vvel(i,k,jj) = - w2d(i,k,1)*rho*con_g/100.  !convert back to hPa/s
                  endif
               enddo
            enddo
            
            ! calculate cloud effective radii
            if ( effr_in ) then
               do k = 1, lev
                  do i = 1, myim(jj)
                     p2d    (i,k,1) = 100.0 * plt(i,k,jj)                      !layer mean pressure (from mb to Pa)
                     dp2d_ef(i,k,1) = (dsigma(k,1)*pst(i,jj)+dsigma(k,2))*100. !difference of interface pressure (Pa)
                     snr1d  (i,1)   = snr(i,jj)                                !snow depth (mm)
                  enddo
               enddo
               
            if ( nmmiph .eq. 12 ) then
                  call cloud_diagnosis_v2                                          &
                     ( 1, myim(jj), 1, lev, land1d, p2d,             &
                       dp2d_ef, t2d, qc2d, qi2d,   &
                       qr2d, qs2d, qg2d, qnc2d,    &
                       rew2d, rei2d, rer2d, res2d, &
                       reg2d, snr1d )
            end if
            
            if ( nmmiph .eq. 13 ) then
                  call cloud_diagnosis_v3                                          &
                     ( 1, myim(jj), 1, lev, land1d, p2d,             &
                       dp2d_ef, t2d, qc2d, qi2d,   &
                       qr2d, qs2d, qg2d, qnc2d,    &
                       rew2d, rei2d, rer2d, res2d, &
                       reg2d, snr1d )
            end if
               do k = 1, lev
                  kc = lev - k + 1
                  do i = 1, myim(jj)
                     re_cloud(i,k,jj)   = rew2d(i,kc,1)   !(micron)
                     re_ice(i,k,jj)     = rei2d(i,kc,1)   !(micron)
                     re_rain(i,k,jj)    = rer2d(i,kc,1)   !(micron)
                     re_snow(i,k,jj)    = res2d(i,kc,1)   !(micron)
                  enddo
               enddo
            endif
            
            do i = 1, myim(jj)
               rlsp(i,jj) = water1d(i,1)+rain1d(i,1)+snow1d(i,1) &
                            +ice1d(i,1)+graupel1d(i,1)     !total large scale precipitation (mm)
               rlspi(i,jj) = ice1d(i,1)
               rlsps(i,jj) = snow1d(i,1)
               rlspg(i,jj) = graupel1d(i,1)
               if ( rlsp(i,jj) .gt. rainmin ) then
                  sr(i,jj) = (snow1d(i,1)+ice1d(i,1)+graupel1d(i,1))  &
                          /(water1d(i,1)+rain1d(i,1)+snow1d(i,1)      &
                          +ice1d(i,1)+graupel1d(i,1))  !snow ratio
               else
                  sr(i,jj) = 0.0
               endif
            enddo
         
            
            deallocate                                                      &
             ( qv2d,qc2d,qr2d,qi2d,qs2d,qg2d,qnc2d,qni2d,cld2d,w2d,t2d,u2d, &
               v2d,dp2d,dz2d,hs,gsize,rain1d,snow1d,ice1d,graupel1d,water1d,&
               land1d,rhc2d )
            if ( effr_in ) deallocate                                       &
             ( rew2d,rei2d,rer2d,res2d,reg2d,snr1d,p2d,dp2d_ef )
   #ifdef EXT_DIAG
            deallocate                                                      &
             ( prefluxw,prefluxr,prefluxi,prefluxs,prefluxg,cond0,dep0,     &
               evap0,sub0,q_con,cappa,te )
   #endif
         end do
      endif  !end if nmmiph.eq.12 .or nmmiph.eq.13



   !     Goddard (GCE) MP
      if ( (nmmiph .eq. 15) .or. (nmmiph .eq. 16) ) then
         allocate                                                                     &
          ( th3d(nx,lev,my_max),qv3d(nx,lev,my_max),qc3d(nx,lev,my_max),              &
            qr3d(nx,lev,my_max),qi3d(nx,lev,my_max),qs3d(nx,lev,my_max),              &
            qg3d(nx,lev,my_max),rho3d(nx,lev,my_max),pii3d(nx,lev,my_max),            &
            p3d(nx,lev,my_max),z3d(nx,lev,my_max),dz3d(nx,lev,my_max),                &
            rain2d(nx,my_max),snow2d(nx,my_max),graupel2d(nx,my_max),sr2d(nx,my_max), &
            ice2d(nx,my_max),ht(nx,my_max),land2d(nx,my_max),w3d(nx,lev,my_max) )
         allocate                                                                     &
          ( rew3d(nx,lev,my_max),rer3d(nx,lev,my_max),rei3d(nx,lev,my_max),           &
            res3d(nx,lev,my_max),reg3d(nx,lev,my_max) )
#ifdef EXT_DIAG
         allocate                                                                     &
          ( preci3d(nx,lev,my_max),precs3d(nx,lev,my_max),precg3d(nx,lev,my_max),     &
            precr3d(nx,lev,my_max),prech3d(nx,lev,my_max),refl_10cm(nx,lev,my_max) )
#endif   
         
         if ( nmmiph.eq.16 ) allocate                                                 &
          ( qh3d(nx,lev,my_max),hail2d(nx,my_max),reh3d(nx,lev,my_max) )

#ifdef Readaeroclx
      allocate ( aero4d(nx,lev,naero,my_max) )
      !$acc enter data copyin(aeroclx) async(async_id)
      !$acc enter data create(aero4d) async(async_id)
#endif
         !$acc enter data create(th3d, qv3d, qc3d, qr3d, &
         !$acc&      qi3d, qs3d, rain2d, ice2d, snow2d, &
         !$acc&      graupel2d, sr2d, qg3d, rew3d, rer3d, &
         !$acc&      rei3d, res3d, reg3d, myim, xlat_myim) async(async_id)
         !$acc enter data create(rho3d, pii3d, p3d, z3d, ht, dz3d, &
         !$acc&      w3d, land2d) async(async_id)
         
         !$acc parallel loop gang vector collapse(3) async(async_id)
         do jj = 1, jlistnum
            do k = 1, lev
               do i = 1, nx
                  th3d(i,k,jj) = 0.
                  qv3d(i,k,jj) = 0.
                  qc3d(i,k,jj) = 0.
                  qr3d(i,k,jj) = 0.
                  qi3d(i,k,jj) = 0.
                  qs3d(i,k,jj) = 0.
                  qg3d(i,k,jj) = 0.
                  rho3d(i,k,jj) = 0.
                  pii3d(i,k,jj) = 0.
                  p3d(i,k,jj) = 0.
                  z3d(i,k,jj) = 0.
                  dz3d(i,k,jj) = 0.
                  rew3d(i,k,jj) = 0.
                  rer3d(i,k,jj) = 0.
                  rei3d(i,k,jj) = 0.
                  res3d(i,k,jj) = 0.
                  reg3d(i,k,jj) = 0.
                  w3d(i,k,jj) = 0.
#ifdef EXT_DIAG
                  preci3d(i,k,jj) = 0.
                  precs3d(i,k,jj) = 0.
                  precg3d(i,k,jj) = 0.
                  precr3d(i,k,jj) = 0.
                  prech3d(i,k,jj) = 0.
                  refl_10cm(i,k,jj) = 0.
#endif   
#ifdef Readaeroclx
                  !$acc loop seq
                  do n = 1, naero
                     aero4d(i,k,n,jj) = 0.
                  end do
#endif   
               end do
            end do
         end do
         !$acc parallel loop gang vector collapse(2) async(async_id)
         do jj = 1, jlistnum
            do i = 1, nx
               ht(i,jj) = 0.
               rain2d(i,jj) = 0.
               ice2d(i,jj) = 0.
               snow2d(i,jj) = 0.
               graupel2d(i,jj) = 0.
               sr2d(i,jj) = 0.
               land2d(i,jj) = 0.
            end do
         end do
         if ( nmmiph .eq. 16 ) then
            do jj = 1, jlistnum
               do k = 1, lev
                  do i = 1, nx
                     qh3d(i,k,jj) = 0.
                     reh3d(i,k,jj) = 0.
                  end do
               end do
            end do
            do jj = 1, jlistnum
               do i = 1, nx
                  hail2d(i,jj) = 0.
               end do
            end do
         endif
            
         nmgce3 = 2   ! nmgce3=1 : WRF
                      ! nmgce3=2 : NASA Unified WRF
                      ! nmgce3=3 : CWB WRF(not used)
         
         ihail = 0  !run gsfcgce with graupel option
         ICE2  = 0  !run gsfcgce with snow, ice and hail/graupel
         
         diagflag = .false.          !if diagflag=true and do_radar_ref=1, call calc_refl10cm
         do_radar_ref = 0
         
         rhowater = 1000.            !water density (kg/m^3), but not used
         rhosnow = 100.              !snow density (kg/m^3), but not used
         !$acc parallel loop gang collapse(2) async(async_id) private(kc, n)
         do jj = 1, jlistnum
            do k = 1, lev
               kc = lev - k + 1
               !$acc loop vector
               do i = 1, myim(jj)
                  qv3d (i,k,jj) = qt(i,             kc,jj)
                  qc3d (i,k,jj) = qt(i,(ntcw-1)*lev+kc,jj)
                  qr3d (i,k,jj) = qt(i,(ntrw-1)*lev+kc,jj)
                  qi3d (i,k,jj) = qt(i,(ntiw-1)*lev+kc,jj)
                  qs3d (i,k,jj) = qt(i,(ntsw-1)*lev+kc,jj)
                  qg3d (i,k,jj) = qt(i,(ntgl-1)*lev+kc,jj)
                  if(nmmiph.eq.16) qh3d (i,k,jj) = qt(i,(nthl-1)*lev+kc,jj)
                  
                  p3d  (i,k,jj) = 100.0*plt(i,kc,jj)                   !layer mean pressure (from mb to Pa)
                  pii3d(i,k,jj) = pk(i,kc,jj)                          !exner function, =(plt/1000)**(Rd/cp)
                  th3d (i,k,jj) = tt(i,kc,jj)/pk(i,kc,jj)                 !potential temperature (K)
                  z3d  (i,k,jj) = phi(i,kc,jj)/con_g                   !layer geopotential height above sea level (m)
                  dz3d (i,k,jj) = (phii(i,k+1,jj)-phii(i,k,jj))/con_g     !layer thickness (m)
                  rho3d(i,k,jj) = p3d(i,k,jj)/                                  &
                                 (con_rd*tt(i,kc,jj)*(1+con_fvirt*qv3d(i,k,jj)))  !air density (kg/m^3)
                  w3d  (i,k,jj) = -vvel(i,kc,jj)*100./(rho3d(i,k,jj)*con_g)        !vertical velocity (m/s)
               
   #ifdef Readaeroclx
                 !$acc loop seq
                 do n = 1, naero
                    aero4d(i,k,n,jj) = aeroclx(i,(n-1)*lev+kc,jj)
                 enddo
   #endif      
               enddo
            enddo
         end do
         
         !$acc parallel loop gang vector collapse(2) async(async_id)
         do jj = 1, jlistnum
            do i = 1, nx
               if (i .le. myim(jj)) then
                  ht(i,jj) = sgeo(i,jj)/con_g   !terrain geopotential height above sea level (m)
                  if( islimsk(i,jj) .eq. 1 ) then
                     land2d(i,jj) = 1.        !land
                     if ( ivegtyp(i,jj) .eq. 15 ) land2d(i,jj) = 2.   !glacial is seen as ocean
                  else
                     land2d(i,jj) = 2.        !ocean & seaice
                  endif
               end if
            enddo
         end do
         if ( nmmiph .eq. 15 ) then 
            if ( nmgce3 .eq. 2 ) then
                  call gsfcgce_3ice_nuwrf_gpu                                   &
                         ( myim, th3d, qv3d, qc3d, qr3d, qi3d, qs3d,            &
                           rho3d, pii3d, p3d, dta, z3d,                         &
                           ht, dz3d, con_g, w3d,                                &
                           itimestep, xlat_myim, sdec, land2d,                  &
                           1, nx , 1, jlistnum, 1, lev,                         & ! memory dims
                           1, nx, 1, jlistnum, 1, lev,                          & ! tile   dims
                           rain2d, ice2d, snow2d, graupel2d, sr2d,              &
                           .false., qg3d,                                       &
      #ifdef Readaeroclx
                           aero4d, naero,                                       &
      #endif      
                           ihail, ice2,                                         &
      #ifdef EXT_DIAG
                           refl_10cm, diagflag, do_radar_ref,                   &
                           preci3d, precs3d, precg3d, precr3d,                  &
      #endif      
                           rew3d, rer3d, rei3d, res3d, reg3d, &
                           SL_sedi, sat_predict, new_saturation, &
                           use_cpm, use_declination, benchmark, async_id)

            end if
            
            if ( nmgce3 .eq. 1 ) then
               do jj = 1, jlistnum
                  call gsfcgce                                                  &
                         ( th3d(1,1,jj), qv3d(1,1,jj), qc3d(1,1,jj),            &
                           qr3d(1,1,jj), qi3d(1,1,jj), qs3d(1,1,jj),            &
                           rho3d(1,1,jj), pii3d(1,1,jj), p3d(1,1,jj), dta,      &
                           z3d(1,1,jj), ht(1,jj), dz3d(1,1,jj), con_g,          &
                           rhowater, rhosnow, itimestep,                        &
                           1, nx , 1, 1, 1, lev,                                & ! memory dims
                           1, myim(jj), 1, 1, 1, lev,                           & ! tile   dims
                           rain2d(1,jj), snow2d(1,jj),                          &
                           graupel2d(1,jj), sr2d(1,jj),                         &
      #ifdef EXT_DIAG
                           refl_10cm(1,1,jj), diagflag, do_radar_ref,           &
      #endif      
      #ifdef EffectRad_GCE3
                           land2d(1,jj),                                        &
                           rew3d(1,1,jj), rer3d(1,1,jj), rei3d(1,1,jj),         &
                           res3d(1,1,jj), reg3d(1,1,jj),                        &
      #endif      
                           .false., qg3d(1,1,jj), ihail, ice2 )
               end do
            end if
            
!            if ( nmgce3 .eq. 3 )                                          &
!            call gsfcgce_3ice_cwb                                         &
!                   ( th3d, qv3d, qc3d, qr3d, qi3d, qs3d, qg3d,            &
!                     rho3d, pii3d, p3d, dta, z3d,                         &
!                     ht, dz3d, itimestep,                                 &
!                     1, nx , 1, 1, 1, lev,                                & ! memory dims
!                     1, myim(jj), 1, 1, 1, lev,                                & ! tile   dims
!                     sr2d, rainnc2d, rain2d,                              &
!                     snownc2d, snow2d, graupelnc2d, graupel2d,            &
!                     land2d )
         
         endif
            
         if ( nmmiph .eq. 16 ) then
            do jj = 1, jlistnum
               dx = sqrt(area(jj))             !grid length (m)
               call gsfcgce_4ice_nuwrf                                          &
                      ( th3d(1,1,jj), qv3d(1,1,jj), qc3d(1,1,jj), qr3d(1,1,jj), &
                        qi3d(1,1,jj), qs3d(1,1,jj), qh3d(1,1,jj), qg3d(1,1,jj), &
                        rho3d(1,1,jj), pii3d(1,1,jj), p3d(1,1,jj), dta,         &
                        z3d(1,1,jj), ht(1,jj), dz3d(1,1,jj), con_g, w3d(1,1,jj),&
                        rhowater, rhosnow, itimestep, land2d(1,jj), dx,         &
                        1, nx , 1, 1, 1, lev,                                   & ! memory dims
                        1, myim(jj), 1, 1, 1, lev,                              & ! tile   dims
                        rain2d(1,jj), snow2d(1,jj), graupel2d(1,jj),            &
                        hail2d(1,jj), sr2d(1,jj),                               &
                        rew3d(1,1,jj), rer3d(1,1,jj), rei3d(1,1,jj),            &
                        res3d(1,1,jj), reg3d(1,1,jj), reh3d(1,1,jj),            &
   #ifdef EXT_DIAG
                        refl_10cm(1,1,jj), diagflag, do_radar_ref,              &
                        physc, physe, physd, physs, physm, physf,               & !simultaneous diabatic heating rate
                        acphysc, acphyse, acphysd, acphyss, acphysm,acphysf,    & !accumulated diabatic heating rate
                        preci3d(1,1,jj), precs3d(1,1,jj), precg3d(1,1,jj),      &
                        prech3d(1,1,jj), precr3d(1,1,jj),                       & !precitation
   #endif   
                        .false. )
            end do
         end if
         
         !$acc parallel loop gang collapse(2) async(async_id)
         do jj = 1, jlistnum
            do k = 1, lev
               kc = lev - k + 1
               !$acc loop vector
               do i = 1, myim(jj)
                  qt(i,             k,jj) = qv3d(i,kc,jj)
                  qt(i,(ntcw-1)*lev+k,jj) = qc3d(i,kc,jj)
                  qt(i,(ntrw-1)*lev+k,jj) = qr3d(i,kc,jj)
                  qt(i,(ntiw-1)*lev+k,jj) = qi3d(i,kc,jj)
                  qt(i,(ntsw-1)*lev+k,jj) = qs3d(i,kc,jj)
                  qt(i,(ntgl-1)*lev+k,jj) = qg3d(i,kc,jj)
                  tt(i,k,jj) = th3d(i,kc,jj)*pk(i,k,jj)  !convert back to real temperature
                  
                  if ( nmmiph .eq. 15 ) then
                     if ( nmgce3 .eq. 3 ) then
                     ! Goddard 3ICE with constant effective radii
                        re_cloud(i,k,jj) = 10.0          !micron
                        re_rain (i,k,jj) = 1000.0        !micron
                        re_ice  (i,k,jj) = 50.0          !micron
                        re_snow (i,k,jj) = 250.0         !micron
                     else
                        re_cloud(i,k,jj) = rew3d(i,k,jj)  !micron
                        re_rain (i,k,jj) = rer3d(i,k,jj)  !micron
                        re_ice  (i,k,jj) = rei3d(i,k,jj)  !micron
                        re_snow (i,k,jj) = res3d(i,k,jj)  !micron
                     endif
                  endif
                  
                  if ( nmmiph .eq. 16 ) then
                     qt(i,(nthl-1)*lev+k,jj) = qh3d(i,kc,jj)
                     re_cloud(i,k,jj) = rew3d(i,k,jj)  !micron
                     re_rain (i,k,jj) = rer3d(i,k,jj)  !micron
                     re_ice  (i,k,jj) = rei3d(i,k,jj)  !micron
                     re_snow (i,k,jj) = res3d(i,k,jj)  !micron
                  endif
               
               enddo
            enddo
         end do
         
         !$acc parallel loop gang vector collapse(2) async(async_id)
         do jj = 1, jlistnum
            do i = 1, nx
               if (i .le. myim(jj)) then
                  rlsp(i,jj) = rain2d(i,jj)  !total large scale precipitation (kg/m^2=mm)
                  rlspi(i,jj) = ice2d(i,jj)
                  rlsps(i,jj) = snow2d(i,jj)
                  rlspg(i,jj) = graupel2d(i,jj)
                  if ( nmmiph .eq. 16 ) rlspg(i,jj) = rlspg(i,jj) + hail2d(i,jj)
                  sr(i,jj)   = sr2d(i,jj)
               end if
            enddo
         end do
         ! deallocate GPU temporary arrays (GCE)
         !$acc exit data delete(prsl, del, ttc, qtc, qtr, qtrw, qti, qtsw, &
         !$acc&     qtgl, refl10, ntnc, rainncv, snowncv, graupelncv) async(async_id)
         !$acc exit data delete(th3d, qv3d, qc3d, qr3d, &
         !$acc&     qi3d, qs3d, rain2d, ice2d, snow2d, &
         !$acc&     graupel2d, sr2d, qg3d, rew3d, rer3d, &
         !$acc&     rei3d, res3d, reg3d, myim, xlat_myim) async(async_id)
         !$acc exit data delete(rho3d, pii3d, p3d, z3d, ht, dz3d, &
         !$acc&     w3d, land2d) async(async_id)
#ifdef Readaeroclx
         !$acc exit data delete(aero4d) async(async_id)
         !$acc exit data delete(aeroclx) async(async_id)
         deallocate ( aero4d )
#endif
         !$acc wait(async_id)
         deallocate                                                      &
          ( th3d,qv3d,qc3d,qr3d,qs3d,qi3d,qg3d,pii3d,p3d,z3d,dz3d,rho3d, &
            rain2d,ice2d,snow2d,graupel2d,                               &
            ht,sr2d,land2d,w3d,rew3d,rer3d,rei3d,res3d,reg3d )
#ifdef EXT_DIAG
         deallocate                                                      &
          ( preci3d,precs3d,precg3d,precr3d,prech3d,refl_10cm )
#endif   
         if ( nmmiph .eq. 16 ) deallocate                                &
          ( qh3d,hail2d,reh3d )

      endif  !end of if nmmiph=15 .or. nmmiph=16

      return
!--------------------------
      end subroutine mp_scheme_gpu
!--------------------------
