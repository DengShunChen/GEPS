#define EffectRad_GCE3
!--------------------------
      subroutine mp_init                                               &
!--------------------------
!  ---  inputs:
           ( nmmiph ,myrank)
!  ---  outputs: ( none )
      use const,              only : SL_sedi, sat_predict, &
                                     new_saturation, use_cpm, &
                                     use_declination
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
#ifdef USE_CUDA
      use module_gocart_coupling_gpu, only : makelut_ccn_icn_gpu
      use module_mp_gsfcgce_3ice_nuwrf_gpu, only: gce_table_copyin_gpu
#endif
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
          call new_thompson_init ( is_aerosol_aware,                    &
                   merra2_aerosol_aware, myrank, 0, errmsg, errflg )
          if ( myrank .eq. 0 )                                          &
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
#ifdef Readaeroclx
#ifdef USE_CUDA
          call makelut_ccn_icn_gpu
#else
          call makelut_ccn_icn
#endif
#endif
          if ( myrank .eq. 0 ) then
#ifdef USE_CUDA
       print *,'Goddard (GCE) 3ICE cloud microphysics initialized (GPU)'
#else
             print *,'Goddard (GCE) 3ICE cloud microphysics initialized'
#endif
             print *,'Flag settings:'
             print *,'SL_sedi = ', SL_sedi
             print *,'sat_predict = ', sat_predict 
             print *,'new_saturation = ', new_saturation
             print *,'use_cpm = ', use_cpm 
             print *,'use_declination = ', use_declination
          endif
        endif
! Goddard (GCE) 4ICE MP
        if ( nmmiph .eq. 16 ) then
          if ( myrank .eq. 0 )                                         &
             print *,'Goddard (GCE) 4ICE cloud microphysics initialized'
        endif

      return
!--------------------------
      end subroutine mp_init
!--------------------------
!--------------------------
      subroutine mp_scheme                                             &
!--------------------------
!  ---  inputs:
           ( nmmiph,nx,nxj,lev,ncld,plt,ptop,                          &
             dsigma,phii,islimsk,q0,kdt,tpi,me,dta,area,jj,            &
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
             use_cpm, use_declination )

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
      use module_mp_gsfcgce_3ice_nuwrf, only: gsfcgce_3ice_nuwrf
!      use module_mp_gsfcgce_3ice_cwb  , only: gsfcgce_3ice_cwb
      use module_mp_gsfcgce,   only: gsfcgce
! for Goddard (GCE) 4ICE MP
      use module_mp_gce4ice,   only: gsfcgce_4ice_nuwrf
      use physcons,            only: con_rd,con_fvirt,con_g,con_eps
      use physpara,            only: effr_in
      use const,               only: RTYPE

      implicit none

!  ---  inputs:
      integer,  intent(in)    :: nmmiph,nx,nxj,jj,lev,ncld,kdt,me
!      integer,  intent(in)    :: ntcw,ntrw,ntiw,ntsw,ntgl,ntinc,ntrnc
      integer,  intent(in)    :: islimsk(nx)
      integer,  intent(in)    :: itimestep
      integer,  intent(in)    :: ivegtyp(nx)
      real,     intent(in)    :: tpi,dta,xlat,sdec
      real,     intent(in)    :: phii(nx,lev+1)
      real,     intent(in)    :: area
      real,     intent(in)    :: ptop
      real,     intent(in)    :: rhc_mp(nx,lev)
      real,     intent(in)    :: snr(nx)
      real,     intent(in)    :: plt(nx,lev)
      real(kind=RTYPE), intent(in):: phi(nx,lev),pk(nx,lev),sgeo(nx)
      real(kind=RTYPE), intent(in):: q0(nx,lev*ncld),dsigma(lev,2)
#ifdef Readaeroclx
      integer,  intent(in)    :: naero
      real(kind=RTYPE), intent(in):: aeroclx(nx,lev*naero)
#endif
!  ---  inputs/outputs:
      real, intent(inout) :: tt(nx,lev)
      real,     intent(inout) :: qa(nx,lev)
      real(kind=RTYPE), intent(inout) :: vvel(nx,lev) !mb/s
      real, intent(inout) :: ut(nx,lev),vt(nx,lev)
      real(kind=RTYPE), intent(inout) :: qt(nx,lev*ncld)
      real(kind=RTYPE), intent(in   ) :: pst(nx)
!  ---  outputs:
      real,     intent(inout)   :: re_cloud(nx,lev),re_ice(nx,lev),   &
                                   re_snow(nx,lev),re_rain(nx,lev)
      real,     intent(inout)   :: rlsp(nx),sr(nx)
      real,     intent(inout)   :: rlspi(nx),rlsps(nx),rlspg(nx)
!  ---  local arrays:
      integer   kc,k,i,n
      real      prsl(nx,lev),del(nx,lev)
      real      ttc(nx,lev),qtc(nx,lev),qtr(nx,lev),qtrw(nx,lev),      &
                qti(nx,lev),qtsw(nx,lev),qtgl(nx,lev),ntnc(nx,lev,2),  &
                refl10(nx,lev)
      real      rainncv(nx),snowncv(nx),graupelncv(nx)
      real      icem,tem
      logical   lradar
      logical   convert_dry_q,q_remove_cond
      real,dimension(:),allocatable ::                                  &
              land1d
      real,dimension(:,:),allocatable ::                                &
              qv2d,qc2d,qr2d,qi2d,qs2d,qg2d,qnc2d,qni2d,qnr2d,          &
              rew2d,rer2d,rei2d,res2d,reg2d,                            &
              t2d,dp2d,dz2d,cld2d,w2d,u2d,v2d,rhc2d,p2d
! 2M Thompson MP
      real,dimension(:,:),allocatable ::                                &
              nwfa,nifa,pfils,pflls,vt_dbz_wt
      real,dimension(:),allocatable :: nwfasfc,nifasfc,rainnc,snownc,   &
              icenc,graupelnc,icencv,gridkm
      real,dimension(:,:),allocatable :: rand_pert
      real,dimension(:),allocatable :: spp_prt_list,spp_stddev_cutoff
      character(len=3),dimension(:),allocatable :: spp_var_list
      real    dt_inner
      real    rho
      logical sedi_semi,ext_diag,first_time_step,reset_dBZ,aero_ind_fdb,&
              diagflag,convert_dry_rho
      integer do_radar_ref,rand_perturb_on,has_reqc,has_reqi,has_reqs,  &
              kme_stoch,istep,nsteps,errflg,decfl
      character errmsg
      integer, parameter :: n_var_spp=1
! GFDLMP
      real, parameter ::                                                &
                rainmin=1.0e-10 !(mm)
      logical   hydrostatic,phys_hydrostatic,sedi_w
     !GFDL MP v2 & v3
      real, dimension(:), allocatable ::                                &
                gsize,hs,water1d,rain1d,snow1d,ice1d,graupel1d,snr1d
      real, dimension(:,:), allocatable ::                              &
                q_con,cappa,te,dp2d_ef
#ifdef EXT_DIAG
      real, dimension(:), allocatable ::                                &
                cond0,dep0,evap0,sub0
      real, dimension(:,:), allocatable ::                              &
                prefluxw,prefluxr,prefluxi, prefluxs,prefluxg
#endif
      logical   consv_te,last_step,do_inline_mp
     !GFDL MP v1
      integer, dimension(:), allocatable ::                             &
                mask1d
      real, dimension(:,:), allocatable ::                              &
                garea,land2d,rain2d,snow2d,ice2d,graupel2d,rho2d
      real, dimension(:,:,:), allocatable ::                            &
                qv3d,qc3d,qr3d,qi3d,qs3d,qg3d,cld3d,qnc3d,t3d,w3d,u3d,  &
                v3d,dp3d,dz3d,qvten3d,qcten3d,qrten3d,qiten3d,qsten3d,  &
                qgten3d,cldten3d,uten3d,vten3d,tten3d
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
      logical, intent(in) :: SL_sedi, sat_predict, new_saturation,      &
                 use_cpm, use_declination
      logical :: benchmark = .true.
      integer :: jtest = 0
      convert_dry_q = .true.
      q_remove_cond = .false.
!
! reset all value to zero
      prsl  = 0.
      del   = 0.
      ttc   = 0.
      qtc   = 0.
      qtr   = 0.
      qtrw  = 0.
      qti   = 0.
      qtsw  = 0.
      qtgl  = 0.
      ntnc  = 0.
      refl10= 0.
      rainncv=0.
      snowncv=0.
      graupelncv=0.
!
      lradar= .false.
      icem  =  4./3.*tpi*3.2768*1.e-14*890.
!
      if ( nmmiph.eq.6 .or. nmmiph.eq.8 ) then
        do k=1,lev
          kc=lev-k+1
          do i=1,nxj
            prsl(i,kc)= plt(i,k)*100. ! change to Pa
            del(i,kc) = (dsigma(k,1)*pst(i)+dsigma(k,2))*100. !change to Pa
            qtc(i,kc) = qt(i,             k)
            qtr(i,kc) = qt(i,(ntcw-1)*lev+k)
            qtrw(i,kc)= qt(i,(ntrw-1)*lev+k)
            qti(i,kc) = qt(i,(ntiw-1)*lev+k)
            qtsw(i,kc)= qt(i,(ntsw-1)*lev+k)
            qtgl(i,kc)= qt(i,(ntgl-1)*lev+k)
            ttc(i,kc) = tt(i,             k)
          enddo
        enddo
!
!          WSM6
           if ( nmmiph .eq. 6 )                                        &
           call wsm6(ttc,phii,qtc,qtr,qtrw,qti,qtsw,qtgl,prsl,del,dta, &
                     rainncv,sr,islimsk,re_cloud,re_ice,re_snow,       &
                     1,nx,1,lev,1,nxj,1,lev,snowncv,graupelncv)
!          Thompson
           if ( nmmiph .eq. 8 )then
             allocate( qni2d(nx,lev),qnr2d(nx,lev) )
             do k=1,lev
               kc=lev-k+1
               do i=1,nxj
                 qni2d(i,kc) = qt(i,(ntinc-1)*lev+k)                   &
                      +max(0.,qti(i,kc)-q0(i,(ntiw-1)*lev+k))/icem
                 qnr2d(i,kc) = qt(i,(ntrnc-1)*lev+k)
               enddo
             enddo

             call thompson_driver(1,nx,1,lev,1,nxj,1,lev,qtc,qtr,qtrw, &
                     qti,qtsw,qtgl,qni2d,qnr2d,ttc,                    &
                     prsl,del,dta,kdt,rainncv,sr,islimsk,refl10,       &
                     lradar,re_cloud,re_ice,re_snow,me,phii)

             do k=1,lev
               kc=lev-k+1
               do i=1,nxj
                 qt(i,(ntinc-1)*lev+k)=qni2d(i,kc)
                 qt(i,(ntrnc-1)*lev+k)=qnr2d(i,kc)
               enddo
             enddo
           endif
!
        do i=1,nxj
          rlsp(i) = rainncv(i) + snowncv(i) + graupelncv(i)
          rlspi(i) = 0.   !ice not precipitating
          rlsps(i) = snowncv(i)
          rlspg(i) = graupelncv(i)
        enddo
!
        do k=1,lev
          kc=lev-k+1
          do i=1,nxj
            qt(i,             k) = qtc(i,kc)
            qt(i,(ntcw-1)*lev+k) = qtr(i,kc)
            qt(i,(ntrw-1)*lev+k) = qtrw(i,kc)
            qt(i,(ntiw-1)*lev+k) = qti(i,kc)
            qt(i,(ntsw-1)*lev+k) = qtsw(i,kc)
            qt(i,(ntgl-1)*lev+k) = qtgl(i,kc)
            tt(i,             k) = ttc(i,kc)

          enddo
        enddo
        if ( nmmiph .eq. 8 ) deallocate( qni2d,qnr2d )

      endif

!     2M Thompson
      if ( nmmiph .eq. 18 ) then
        allocate                                                        &
         ( t2d(nx,lev),qv2d(nx,lev),qc2d(nx,lev),qr2d(nx,lev),          &
           qi2d(nx,lev),qs2d(nx,lev),qg2d(nx,lev),qni2d(nx,lev),        &
           qnr2d(nx,lev),qnc2d(nx,lev),                                 &
           rew2d(nx,lev),rei2d(nx,lev),res2d(nx,lev),                   &
           nwfa(nx,lev),nifa(nx,lev),dz2d(nx,lev),                      &
           w2d(nx,lev),pfils(nx,lev),pflls(nx,lev),vt_dbz_wt(nx,lev),   &
           nwfasfc(nx),nifasfc(nx),rainnc(nx),snownc(nx),icenc(nx),     &
           graupelnc(nx),icencv(nx),                                    &
           rand_pert(nx,1),spp_prt_list(n_var_spp),                     &
           spp_stddev_cutoff(n_var_spp),spp_var_list(n_var_spp) )
        if ( cfflag_thom .eq. 2 ) allocate                              &
         ( cld2d(nx,lev),land1d(nx),gridkm(nx) )

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

        t2d=0.
        qv2d=0.
        qc2d=0.
        qr2d=0.
        qi2d=0.
        qs2d=0.
        qg2d=0.
        rew2d=0.
        rei2d=0.
        res2d=0.
        qni2d=0.         !number concentracion of ice
        qnr2d=0.         !number concentracion of rain
        qnc2d=0.         !number concentracion of cloud droplet
        nwfa=0.          !number concentration of water friendly aerosol
        nifa=0.          !number concentration of ice friendly aerosol
        nwfasfc=0.       !nwfa at surface
        nifasfc=0.       !nifa at surface
        w2d=0.
        dz2d=0.
        pfils=0.
        pflls=0.
        vt_dbz_wt=0.
        rainnc=0.        !number concentracion of precipitating rain
        snownc=0.        !number concentracion of precipitating snow
        icenc=0.         !number concentracion of precipitating ice
        graupelnc=0.     !number concentracion of precipitating graupel
        icencv=0.        !amount of precipitating ice
        rand_pert=0.
        spp_stddev_cutoff=0.

        do k = 1, lev
          kc = lev - k + 1
          do i = 1, nxj
            qv2d (i,k) = qt(i,              kc)
            qc2d (i,k) = qt(i,(ntcw-1) *lev+kc)
            qr2d (i,k) = qt(i,(ntrw-1) *lev+kc)
            qi2d (i,k) = qt(i,(ntiw-1) *lev+kc)
            qs2d (i,k) = qt(i,(ntsw-1) *lev+kc)
            qg2d (i,k) = qt(i,(ntgl-1) *lev+kc)
            qni2d(i,k) = qt(i,(ntinc-1)*lev+kc)
            qnr2d(i,k) = qt(i,(ntrnc-1)*lev+kc)
            prsl (i,k) = plt(i,kc)*100.                         !layer pressure (Pa)
            t2d  (i,k) = tt(i,kc)
            dz2d (i,k) = (phii(i,k+1)-phii(i,k))/con_g          !layer width (m)

            !> - Convert specific humidity to water vapor mixing ratio.
            !> - Also, hydrometeor variables are mass or number mixing ratio
            !> - either kg of species per kg of dry air, or per kg of (dry + vapor).
            if ( convert_dry_q ) then
              if ( q_remove_cond ) then
                ! q convert to : hydrometeor mass / ( dry mass )
                tem = qv2d(i,k) + qc2d(i,k) + qr2d(i,k) + qi2d(i,k) +   &
                      qs2d(i,k) + qg2d(i,k)
              else
                ! q convert to : hydrometeor mass / ( dry mass + condensates )
                tem = qv2d(i,k)
              endif
              qv2d (i,k) = qv2d (i,k)/(1. - tem)
              qc2d (i,k) = qc2d (i,k)/(1. - tem)
              qr2d (i,k) = qr2d (i,k)/(1. - tem)
              qi2d (i,k) = qi2d (i,k)/(1. - tem)
              qs2d (i,k) = qs2d (i,k)/(1. - tem)
              qg2d (i,k) = qg2d (i,k)/(1. - tem)
              qni2d(i,k) = qni2d(i,k)/(1. - tem)
              qnr2d(i,k) = qnr2d(i,k)/(1. - tem)
            endif

            ! Ensure non-negative mass mixing ratios of all water variables
            if ( qv2d(i,k) .lt. 0. ) qv2d(i,k) = 1.E-10  !should *never* be identically zero
            if ( qc2d(i,k) .lt. 0. ) qc2d(i,k) = 0.
            if ( qr2d(i,k) .lt. 0. ) qr2d(i,k) = 0.
            if ( qi2d(i,k) .lt. 0. ) qi2d(i,k) = 0.
            if ( qs2d(i,k) .lt. 0. ) qs2d(i,k) = 0.
            if ( qg2d(i,k) .lt. 0. ) qg2d(i,k) = 0.

            if ( convert_dry_q ) then
              ! dry air density : rho = 0.622*p/(Rd*T*(0.622+qv))
              rho = con_eps*prsl(i,k)/                                  &
                    (con_rd*t2d(i,k)*(con_eps+qv2d(i,k)))
            else
              ! moist air density : rho = p/(Rd*T*(1+0.608*qv))
              rho = prsl(i,k)/(con_rd*t2d(i,k)*(1+con_fvirt*qv2d(i,k)))
            endif
            w2d(i,k) = - vvel(i,kc)*100./(rho*con_g)          !vertical velocity (m/s)

            ! 1st guess number concentration where mass non-zero
!            if ( kdt .eq. 1 ) then
!              if ( qi2d(i,k) .gt. 0. )                                  &
!                qni2d(i,k) = make_IceNumber(qi2d(i,k)*rho,t2d(i,k))/rho
!              if ( qr2d(i,k) .gt. 0. )                                  &
!                qnr2d(i,k) = make_RainNumber(qr2d(i,k)*rho,t2d(i,k))/rho
!            endif
          enddo
        enddo

        call new_thompson_driver                                        &
                   ( qv2d,qc2d,qr2d,qi2d,qs2d,qg2d,qni2d,qnr2d,         &
                     qnc2d,nwfa,nifa,nwfasfc,nifasfc,                   &!optional
                     t2d,&!th,pii,                                      &!optional ??
                     prsl,w2d,dz2d,dta,dt_inner,                        &
                     sedi_semi,decfl,islimsk,                           &
                     rainnc,rainncv,                                    &
                     snownc,snowncv,icenc,icencv,graupelnc,graupelncv,  &!optional
                     sr,                                                &
                     refl10,diagflag,do_radar_ref,                      &!optional
                     vt_dbz_wt,                                         &!optional
                     first_time_step,                                   &
                     rew2d,rei2d,res2d,                                 &!optional
                     has_reqc,has_reqi,has_reqs,                        &
                     aero_ind_fdb,                                      &!optional
                     rand_perturb_on,                                   &
                     kme_stoch,                                         &
                     rand_pert,spp_prt_list,spp_var_list,               &
                     spp_stddev_cutoff,n_var_spp,                       &
                     1,nx,1,lev,                                        &
                     1,nxj,1,lev,                                       &
                     reset_dBZ,istep,nsteps,                            &
                     errmsg,errflg,                                     &!optional
                     ext_diag,                                          &
#ifdef EXT_DIAG
                     prw_vcdc, prw_vcde, tpri_inu, tpri_ide_d,          &
                     tpri_ide_s, tprs_ide, tprs_sde_d,                  &
                     tprs_sde_s, tprg_gde_d,                            &
                     tprg_gde_s, tpri_iha, tpri_wfz,                    &
                     tpri_rfz, tprg_rfz, tprs_scw, tprg_scw,            &
                     tprg_rcs, tprs_rcs, tprr_rci, tprg_rcg,            &
                     tprw_vcd_c, tprw_vcd_e, tprr_sml,                  &
                     tprr_gml, tprr_rcg,                                &
                     tprr_rcs, tprv_rev, tten3, qvten3,                 &
                     qrten3, qsten3, qgten3, qiten3, niten3,            &
                     nrten3, ncten3, qcten3,                            &
#endif
                     pfils, pflls )

        if ( cfflag_thom .eq. 2 ) then
          land1d=0.
          gridkm=0.
          cld2d=0.
          do i = 1, nxj
            if ( islimsk(i) .eq. 1 ) then
              land1d(i) = 1.      !land
            else
              land1d(i) = 2.      !ocean or seaice
            endif
            gridkm(i) = sqrt(area)/1000.                !grid length (km)
          enddo

          call cal_cldfra3                                              &
                   ( 1, nx, 1, nxj, 1, lev, 1, lev,                     &
                     cld2d, qv2d, qc2d, qi2d, qs2d, dz2d,               &
                     prsl, t2d, land1d, gridkm,                         &
                     .false., 1.5, .false. )

          do k = 1, lev
            kc = lev - k + 1
            do i = 1, nxj
              qa(i,k) = cld2d(i,kc)
            enddo
          enddo
        endif

        do k = 1, lev
          kc = lev - k + 1
          do i = 1, nxj
            if ( convert_dry_q ) then
              if ( q_remove_cond ) then
                tem = qv2d(i,kc) + qc2d(i,kc) + qr2d(i,kc) +            &
                      qi2d(i,kc) + qs2d(i,kc) + qg2d(i,kc)
              else
                tem = qv2d(i,kc)
              endif
              qv2d (i,kc) = qv2d (i,kc)/(1. + tem)
              qc2d (i,kc) = qc2d (i,kc)/(1. + tem)
              qr2d (i,kc) = qr2d (i,kc)/(1. + tem)
              qi2d (i,kc) = qi2d (i,kc)/(1. + tem)
              qs2d (i,kc) = qs2d (i,kc)/(1. + tem)
              qg2d (i,kc) = qg2d (i,kc)/(1. + tem)
              qni2d(i,kc) = qni2d(i,kc)/(1. + tem)
              qnr2d(i,kc) = qnr2d(i,kc)/(1. + tem)
            endif

            qt(i,              k) = max( qv2d (i,kc) , 0. )
            qt(i,(ntcw-1) *lev+k) = max( qc2d (i,kc) , 0. )
            qt(i,(ntrw-1) *lev+k) = max( qr2d (i,kc) , 0. )
            qt(i,(ntiw-1) *lev+k) = max( qi2d (i,kc) , 0. )
            qt(i,(ntsw-1) *lev+k) = max( qs2d (i,kc) , 0. )
            qt(i,(ntgl-1) *lev+k) = max( qg2d (i,kc) , 0. )
            qt(i,(ntinc-1)*lev+k) = max( qni2d(i,kc) , 0. )
            qt(i,(ntrnc-1)*lev+k) = max( qnr2d(i,kc) , 0. )
            tt(i,k) = t2d(i,kc)

            re_cloud(i,k) = rew2d(i,k)*1.E+6   ! m to micron
            re_ice  (i,k) = rei2d(i,k)*1.E+6   ! m to micron
            re_snow (i,k) = res2d(i,k)*1.E+6   ! m to micron
          enddo
        enddo
        do i = 1, nxj
          rlsp(i) = rainncv(i)  !total large scale precipitation (kg/m^2=mm)
          rlspi(i) = icencv(i)      !not sure
          rlsps(i) = snowncv(i)     !not sure
          rlspg(i) = graupelncv(i)  !not sure
        enddo

        deallocate                                                      &
         ( t2d,qv2d,qc2d,qr2d,qi2d,qs2d,qg2d,qni2d,qnr2d,qnc2d,         &
           nwfa,nifa,dz2d,w2d,pfils,pflls,vt_dbz_wt,nwfasfc,nifasfc,    &
           rew2d,rei2d,res2d,rainnc,snownc,icenc,graupelnc,icencv,      &
           rand_pert,spp_prt_list,spp_stddev_cutoff,spp_var_list )
        if ( cfflag_thom .eq. 2 ) deallocate ( cld2d,land1d,gridkm )
      endif  !end of if nmmiph=18

!     GFDL MP v1
      if ( nmmiph .eq. 11 ) then
        allocate                                                        &
         ( qv3d(nxj,1,lev),qc3d(nxj,1,lev),qr3d(nxj,1,lev),             &
           qi3d(nxj,1,lev),qs3d(nxj,1,lev),qg3d(nxj,1,lev),             &
           cld3d(nxj,1,lev),qnc3d(nxj,1,lev),t3d(nxj,1,lev),            &
           w3d(nxj,1,lev),u3d(nxj,1,lev),v3d(nxj,1,lev),                &
           dp3d(nxj,1,lev),dz3d(nxj,1,lev),                             &
           qvten3d(nxj,1,lev),qcten3d(nxj,1,lev),qrten3d(nxj,1,lev),    &
           qiten3d(nxj,1,lev),qsten3d(nxj,1,lev),qgten3d(nxj,1,lev),    &
           cldten3d(nxj,1,lev),uten3d(nxj,1,lev),vten3d(nxj,1,lev),     &
           tten3d(nxj,1,lev),                                           &
           rew2d(nxj,lev),rei2d(nxj,lev),rer2d(nxj,lev),res2d(nxj,lev), &
           reg2d(nxj,lev),land2d(nxj,1),rain2d(nxj,1),snow2d(nxj,1),    &
           ice2d(nxj,1),graupel2d(nxj,1),garea(nxj,1),rhc2d(nxj,lev) )
        if ( effr_in ) allocate                                         &
           ( dp2d(nxj,lev),rho2d(nxj,lev),qc2d(nxj,lev),qr2d(nxj,lev),  &
             qi2d(nxj,lev),qs2d(nxj,lev),qg2d(nxj,lev),mask1d(nxj),     &
             t2d(nxj,lev) )
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

        do i = 1, nxj
          if( islimsk(i) .eq. 1 ) land2d(i,1) = 1.  !land fraction
          if( effr_in ) mask1d(i) = islimsk(i)      !land-sea mask
          garea(i,1) = area                         !area of grid box (m^-2)
        enddo

        do k = 1, lev
          kc = lev - k + 1
          do i = 1, nxj
            if ( sedi_w ) then
              prsl(i,k) = 100.0 * plt(i,k)             !layer mean pressure (from mb to Pa)
              w3d(i,1,k) = -vvel(i,k)*(1.+con_fvirt*qt(i,k))*tt(i,k)    &
                          /prsl(i,k)*con_rd/con_g      !vertical velocity (m/s)
            else
              w3d(i,1,k) = 0.
            endif

            qv3d (i,1,k) = qt(i,             k)
            qc3d (i,1,k) = qt(i,(ntcw-1)*lev+k)
            qr3d (i,1,k) = qt(i,(ntrw-1)*lev+k)
            qi3d (i,1,k) = qt(i,(ntiw-1)*lev+k)
            qs3d (i,1,k) = qt(i,(ntsw-1)*lev+k)
            qg3d (i,1,k) = qt(i,(ntgl-1)*lev+k)
            qnc3d(i,1,k) = 0.                      ! =0. for prog_ccn=.false. (cm^-3)
            cld3d(i,1,k) = 0.                      !layer cloud fraction (=0. for do_qa=.false.)
            t3d  (i,1,k) = tt(i,k)                 !temperature
            u3d  (i,1,k) = ut(i,k)                 !zonal wind (m/s)
            v3d  (i,1,k) = vt(i,k)                 !meridional wind (m/s)
            dp3d (i,1,k) = (dsigma(k,1)*pst(i)+dsigma(k,2))*100. !difference of interface pressure (Pa)
            dz3d (i,1,k) = (phii(i,kc)-phii(i,kc+1))/con_g       !differences of height (m), dz<0
            rhc2d(i,  k) = rhc_mp(i,k)
            if ( effr_in ) dp2d(i,k) = dp3d(i,1,k)
          enddo
        enddo

        call gfdl_cloud_microphys_driver                                &
                ( qv3d, qc3d, qr3d, qi3d, qs3d, qg3d, cld3d, qnc3d,     &
                  qvten3d, qcten3d, qrten3d, qiten3d, qsten3d, qgten3d, &
                  cldten3d, tten3d, t3d, w3d, u3d, v3d, uten3d, vten3d, &
                  dz3d, dp3d, garea, dta, land2d,                       &
                  rain2d, snow2d, ice2d, graupel2d,                     &
                  rhc2d, hydrostatic, phys_hydrostatic,                 &
                  1, nxj, 1, 1, 1, lev, 1, lev )

        do k = 1, lev
          do i = 1, nxj
            qc2d(i,k) = qc3d(i,1,k) + qcten3d(i,1,k) * dta
            qr2d(i,k) = qr3d(i,1,k) + qrten3d(i,1,k) * dta
            qi2d(i,k) = qi3d(i,1,k) + qiten3d(i,1,k) * dta
            qs2d(i,k) = qs3d(i,1,k) + qsten3d(i,1,k) * dta
            qg2d(i,k) = qg3d(i,1,k) + qgten3d(i,1,k) * dta

            qt(i,             k) = qv3d(i,1,k) + qvten3d(i,1,k) * dta
            qt(i,(ntcw-1)*lev+k) = qc2d(i,k)
            qt(i,(ntrw-1)*lev+k) = qr2d(i,k)
            qt(i,(ntiw-1)*lev+k) = qi2d(i,k)
            qt(i,(ntsw-1)*lev+k) = qs2d(i,k)
            qt(i,(ntgl-1)*lev+k) = qg2d(i,k)
            qa(i,k)  = cld3d(i,1,k) + cldten3d(i,1,k) * dta
            tt(i,k)  = t3d(i,1,k) + tten3d(i,1,k) * dta
            ut(i,k)  = u3d(i,1,k) + uten3d(i,1,k) * dta
            vt(i,k)  = v3d(i,1,k) + vten3d(i,1,k) * dta

            if ( sedi_w ) then
              vvel(i,k) = -w3d(i,1,k)*prsl(i,k)*con_g/con_rd            &
                          /((1+con_fvirt*qt(i,k))*tt(i,k))
            endif

            if ( effr_in ) then
              rho2d(i,k) = 0.622*prsl(i,k)                              &
                          /( con_rd*tt(i,k)*(qt(i,k)+0.622) ) !air density (kg/m^3)
            endif
          enddo
        enddo

        if ( effr_in ) then
          do k = 1, lev
            do i = 1, nxj
              t2d(i,k) = tt(i,k)
            enddo
          enddo
          call cloud_diagnosis                                          &
               ( 1, nxj, 1, lev, rho2d, dp2d, mask1d,                   &
                 qc2d, qi2d, qr2d, qs2d, qg2d, t2d,                     &
                 rew2d, rei2d, rer2d, res2d, reg2d )
          do k = 1, lev
            kc = lev - k + 1
            do i = 1, nxj
              re_cloud(i,k)   = rew2d(i,kc)   !(micron)
              re_ice(i,k)     = rei2d(i,kc)   !(micron)
              re_rain(i,k)    = rer2d(i,kc)   !(micron)
              re_snow(i,k)    = res2d(i,kc)   !(micron)
            enddo
          enddo
        endif
!
        do i = 1, nxj
          if ( rain2d(i,1)    .lt. rainmin ) rain2d(i,1)    = 0.0
          if ( ice2d(i,1)     .lt. rainmin ) ice2d(i,1)     = 0.0
          if ( snow2d(i,1)    .lt. rainmin ) snow2d(i,1)    = 0.0
          if ( graupel2d(i,1) .lt. rainmin ) graupel2d(i,1) = 0.0

          rlsp(i) = rain2d(i,1)+snow2d(i,1)+ice2d(i,1)+graupel2d(i,1)  !total large scale precipitation (mm)
          rlspi(i) = ice2d(i,1)
          rlsps(i) = snow2d(i,1)
          rlspg(i) = graupel2d(i,1)
          if ( rlsp(i) .gt. rainmin ) then
            sr(i) = (snow2d(i,1)+ice2d(i,1)+graupel2d(i,1))             &
                   /(rain2d(i,1)+snow2d(i,1)+ice2d(i,1)+graupel2d(i,1))  !snow ratio
          else
            sr(i) = 0.0
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

!     GFDL MP v2 & v3
      if ( nmmiph .eq. 12 .or. nmmiph .eq. 13 ) then
        allocate                                                        &
         ( qv2d(nxj,lev),qc2d(nxj,lev),qr2d(nxj,lev),qi2d(nxj,lev),     &
           qs2d(nxj,lev),qg2d(nxj,lev),cld2d(nxj,lev),qnc2d(nxj,lev),   &
           qni2d(nxj,lev),w2d(nxj,lev),t2d(nxj,lev),u2d(nxj,lev),       &
           v2d(nxj,lev),dp2d(nxj,lev),dz2d(nxj,lev),                    &
           hs(nxj),land1d(nxj),gsize(nxj),rain1d(nxj),snow1d(nxj),      &
           ice1d(nxj),graupel1d(nxj),water1d(nxj),rhc2d(nxj,lev) )
        if ( effr_in ) allocate                                         &
         ( rew2d(nxj,lev),rei2d(nxj,lev),rer2d(nxj,lev),res2d(nxj,lev), &
           reg2d(nxj,lev),snr1d(nxj),p2d(nxj,lev),dp2d_ef(nxj,lev) )
#ifdef EXT_DIAG
        allocate                                                        &
         ( prefluxr(nxj,lev),prefluxi(nxj,lev),prefluxs(nxj,lev),       &
           prefluxg(nxj,lev),prefluxw(nxj,lev),                         &
           q_con(nxj,lev),cappa(nxj,lev),te(nxj,lev) )
        allocate                                                        &
         ( cond0(nxj),dep0(nxj),evap0(nxj),sub0(nxj) )
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

        gsize = 0.0
        land1d = 0.0
        rain1d = 0.
        snow1d = 0.
        ice1d  = 0.
        graupel1d = 0.
        water1d = 0.
#ifdef EXT_DIAG
        te    = 0.0
        q_con = 0.0  !not sure
        cappa = 0.0  !not sure
        cond0 = 0.0
        dep0  = 0.0
        evap0 = 0.0
        sub0  = 0.0
        prefluxw  = 0.0
        prefluxr  = 0.0
        prefluxi  = 0.0
        prefluxs  = 0.0
        prefluxg  = 0.0
#endif

        do i = 1, nxj
          gsize(i) = sqrt(area)     !square root of grid area (m)
          hs(i)    = sgeo(i)        !terrain geopotential (m^2 s^-2)
          if( islimsk(i) .eq. 1 ) land1d(i) = 1.  !land fraction
        enddo
        do k = 1, lev
          kc = lev - k + 1
          do i = 1, nxj
            qv2d (i,k) = qt(i,             k)
            qc2d (i,k) = qt(i,(ntcw-1)*lev+k)
            qr2d (i,k) = qt(i,(ntrw-1)*lev+k)
            qi2d (i,k) = qt(i,(ntiw-1)*lev+k)
            qs2d (i,k) = qt(i,(ntsw-1)*lev+k)
            qg2d (i,k) = qt(i,(ntgl-1)*lev+k)
            qnc2d(i,k) = 0.
            qni2d(i,k) = 0.
            cld2d(i,k) = 0.                      !layer cloud fracion, =0 for do_qa=.false.
            t2d  (i,k) = tt(i,k)                 !temperature
            u2d  (i,k) = ut(i,k)                 !zonal wind (m/s)
            v2d  (i,k) = vt(i,k)                 !meridional wind (m/s)
            dp2d (i,k) = (dsigma(k,1)*pst(i)+dsigma(k,2))*100. !difference of interface pressure (Pa)
            dz2d (i,k) = (phii(i,kc)-phii(i,kc+1))/con_g       !differences of height (m), dz<0
            rhc2d(i,k) = rhc_mp(i,k)

            if ( sedi_w ) then
              prsl(i,k) = 100.0 * plt(i,k)            !layer mean pressure (from mb to Pa)
              ! use moist air density : rho = p/(Rd*T*(1+0.608*qv))
              rho = prsl(i,k)/(con_rd*t2d(i,k)*(1+con_fvirt*qv2d(i,k)))
              w2d(i,k) = - vvel(i,k)*100./(rho*con_g) !vertical velocity (m/s)
            else
              w2d(i,k)  = 0.
            endif
          enddo
        enddo

        ! GFDL MP v2
        if ( nmmiph .eq. 12 )                                           &
          call gfdlv2_driver                                            &
                ( qv2d, qc2d, qr2d, qi2d, qs2d, qg2d,                   &
                  cld2d, qnc2d,qni2d,                                   &
                  t2d, w2d, u2d, v2d, dz2d, dp2d, gsize, dta, hs,       &
                  land1d,                                               &
                  rain1d, snow1d, ice1d, graupel1d, hydrostatic,        &
                  1, nxj, 1, lev, consv_te,                             &
#ifdef EXT_DIAG
                  q_con, cappa, te,                                     &
                  prefluxr, prefluxi, prefluxs, prefluxg,               &
                  cond0, dep0, evap0, sub0,                             &
#endif
                  rhc2d, last_step, do_inline_mp )

        ! GFDL MP v3
        if ( nmmiph .eq. 13 )                                           &
          call gfdlv3_driver                                            &
                ( qv2d, qc2d, qr2d, qi2d, qs2d, qg2d,                   &
                  cld2d, qnc2d, qni2d,                                  &
                  t2d, w2d, u2d, v2d, dz2d, dp2d, gsize, dta, hs,       &
                  land1d, water1d,                                      &
                  rain1d, snow1d, ice1d, graupel1d, hydrostatic,        &
                  1, nxj, 1, lev, consv_te,                             &
#ifdef EXT_DIAG
                  q_con, cappa, te,                                     &
                  prefluxw, prefluxr, prefluxi, prefluxs, prefluxg,     &
                  cond0, dep0, evap0, sub0,                             &
#endif
                  last_step, do_inline_mp )

        do k = 1, lev
          do i = 1, nxj
            qt(i,             k) = qv2d(i,k)
            qt(i,(ntcw-1)*lev+k) = qc2d(i,k)
            qt(i,(ntrw-1)*lev+k) = qr2d(i,k)
            qt(i,(ntiw-1)*lev+k) = qi2d(i,k)
            qt(i,(ntsw-1)*lev+k) = qs2d(i,k)
            qt(i,(ntgl-1)*lev+k) = qg2d(i,k)
            qa(i,k)  = cld2d(i,k)
            tt(i,k)  = t2d  (i,k)
            ut(i,k)  = u2d  (i,k)
            vt(i,k)  = v2d  (i,k)

            if ( sedi_w ) then
              ! use moist air density : rho = p/(Rd*T*(1+0.608*qv))
              rho = prsl(i,k)/(con_rd*tt(i,k)*(1+con_fvirt*qt(i,k)))
              vvel(i,k) = - w2d(i,k)*rho*con_g/100.  !convert back to hPa/s
            endif
          enddo
        enddo

        ! calculate cloud effective radii
        if ( effr_in ) then
          do k = 1, lev
            do i = 1, nxj
              p2d    (i,k) = 100.0 * plt(i,k)                      !layer mean pressure (from mb to Pa)
              dp2d_ef(i,k) = (dsigma(k,1)*pst(i)+dsigma(k,2))*100. !difference of interface pressure (Pa)
              snr1d  (i)   = snr(i)                                !snow depth (mm)
            enddo
          enddo

          if ( nmmiph .eq. 12 )                                         &
            call cloud_diagnosis_v2                                     &
               ( 1, nxj, 1, lev, land1d, p2d, dp2d_ef, t2d,             &
                 qc2d, qi2d, qr2d, qs2d, qg2d, qnc2d,                   &
                 rew2d, rei2d, rer2d, res2d, reg2d, snr1d )

          if ( nmmiph .eq. 13 )                                         &
            call cloud_diagnosis_v3                                     &
               ( 1, nxj, 1, lev, land1d, p2d, dp2d_ef, t2d,             &
                 qc2d, qi2d, qr2d, qs2d, qg2d, qnc2d,                   &
                 rew2d, rei2d, rer2d, res2d, reg2d, snr1d )

          do k = 1, lev
            kc = lev - k + 1
            do i = 1, nxj
              re_cloud(i,k)   = rew2d(i,kc)   !(micron)
              re_ice(i,k)     = rei2d(i,kc)   !(micron)
              re_rain(i,k)    = rer2d(i,kc)   !(micron)
              re_snow(i,k)    = res2d(i,kc)   !(micron)
            enddo
          enddo
        endif

        do i = 1, nxj
          rlsp(i) = water1d(i)+rain1d(i)+snow1d(i)+ice1d(i)+graupel1d(i)     !total large scale precipitation (mm)
          rlspi(i) = ice1d(i)
          rlsps(i) = snow1d(i)
          rlspg(i) = graupel1d(i)
          if ( rlsp(i) .gt. rainmin ) then
            sr(i) = (snow1d(i)+ice1d(i)+graupel1d(i))                   &
                    /(water1d(i)+rain1d(i)+snow1d(i)+ice1d(i)+graupel1d(i))  !snow ratio
          else
            sr(i) = 0.0
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
      endif  !end if nmmiph.eq.12 .or nmmiph.eq.13

!     Goddard (GCE) MP
      if ( (nmmiph .eq. 15) .or. (nmmiph .eq. 16) ) then
        allocate                                                        &
         ( th3d(nx,lev,1),qv3d(nx,lev,1),qc3d(nx,lev,1),qr3d(nx,lev,1), &
           qi3d(nx,lev,1),qs3d(nx,lev,1),qg3d(nx,lev,1),rho3d(nx,lev,1),&
           pii3d(nx,lev,1),p3d(nx,lev,1),z3d(nx,lev,1),dz3d(nx,lev,1),  &
           rain2d(nx,1),snow2d(nx,1),graupel2d(nx,1),sr2d(nx,1),        &
           ice2d(nx,1),ht(nx,1),land2d(nx,1),w3d(nx,lev,1) )
        allocate                                                        &
         ( rew3d(nx,lev,1),rer3d(nx,lev,1),rei3d(nx,lev,1),             &
           res3d(nx,lev,1),reg3d(nx,lev,1) )
#ifdef EXT_DIAG
        allocate                                                        &
         ( preci3d(nx,lev,1),precs3d(nx,lev,1),precg3d(nx,lev,1),       &
           precr3d(nx,lev,1),prech3d(nx,lev,1),refl_10cm(nx,lev,1) )
#endif

        if ( nmmiph.eq.16 ) allocate                                    &
         ( qh3d(nx,lev,1),hail2d(nx,1),reh3d(nx,lev,1) )

        th3d = 0.
        qv3d = 0.
        qc3d = 0.
        qr3d = 0.
        qi3d = 0.
        qs3d = 0.
        qg3d = 0.
        rho3d = 0.
        pii3d = 0.
        p3d = 0.
        z3d = 0.
        ht = 0.
        dz3d = 0.
        rain2d = 0.
        ice2d = 0.
        snow2d = 0.
        graupel2d = 0.
        sr2d = 0.
        land2d = 0.
        rew3d = 0.
        rer3d = 0.
        rei3d = 0.
        res3d = 0.
        reg3d = 0.
        w3d = 0.
#ifdef Readaeroclx
        allocate ( aero4d(nx,lev,1,naero) )
        aero4d = 0.
#endif

#ifdef EXT_DIAG
        preci3d = 0.
        precs3d = 0.
        precg3d = 0.
        precr3d = 0.
        prech3d = 0.
        refl_10cm = 0.
#endif

        if ( nmmiph .eq. 16 ) then
          qh3d = 0.
          hail2d = 0.
          reh3d = 0.
        endif

        nmgce3 = 2   ! nmgce3=1 : WRF
                     ! nmgce3=2 : NASA Unified WRF
                     ! nmgce3=3 : CWB WRF(not used)

        ihail = 0  !run gsfcgce with graupel option
        ICE2  = 0  !run gsfcgce with snow, ice and hail/graupel

        diagflag = .false.          !if diagflag=true and do_radar_ref=1, call calc_refl10cm
        do_radar_ref = 0

        dx = sqrt(area)             !grid length (m)
        rhowater = 1000.            !water density (kg/m^3), but not used
        rhosnow = 100.              !snow density (kg/m^3), but not used

        do k = 1, lev
          kc = lev - k + 1
          do i = 1, nxj
            qv3d (i,k,1) = qt(i,             kc)
            qc3d (i,k,1) = qt(i,(ntcw-1)*lev+kc)
            qr3d (i,k,1) = qt(i,(ntrw-1)*lev+kc)
            qi3d (i,k,1) = qt(i,(ntiw-1)*lev+kc)
            qs3d (i,k,1) = qt(i,(ntsw-1)*lev+kc)
            qg3d (i,k,1) = qt(i,(ntgl-1)*lev+kc)
            if(nmmiph.eq.16) qh3d (i,k,1) = qt(i,(nthl-1)*lev+kc)

            p3d  (i,k,1) = 100.0*plt(i,kc)                   !layer mean pressure (from mb to Pa)
            pii3d(i,k,1) = pk(i,kc)                          !exner function, =(plt/1000)**(Rd/cp)
            th3d (i,k,1) = tt(i,kc)/pk(i,kc)                 !potential temperature (K)
            z3d  (i,k,1) = phi(i,kc)/con_g                   !layer geopotential height above sea level (m)
            dz3d (i,k,1) = (phii(i,k+1)-phii(i,k))/con_g     !layer thickness (m)
            rho3d(i,k,1) = p3d(i,k,1)/                                  &
                           (con_rd*tt(i,kc)*(1+con_fvirt*qv3d(i,k,1)))  !air density (kg/m^3)
            w3d  (i,k,1) = -vvel(i,kc)*100./(rho3d(i,k,1)*con_g)        !vertical velocity (m/s)

#ifdef Readaeroclx
            do n = 1, naero
              aero4d(i,k,1,n) = aeroclx(i,(n-1)*lev+kc)
            enddo
#endif
          enddo
        enddo

        do i = 1, nxj
          ht(i,1) = sgeo(i)/con_g   !terrain geopotential height above sea level (m)
          if( islimsk(i) .eq. 1 ) then
            land2d(i,1) = 1.        !land
            if ( ivegtyp(i) .eq. 15 ) land2d(i,1) = 2.   !glacial is seen as ocean
          else
            land2d(i,1) = 2.        !ocean & seaice
          endif
        enddo

        if ( nmmiph .eq. 15 ) then 
          if ( nmgce3 .eq. 2 )                                          &
          call gsfcgce_3ice_nuwrf                                       &
                 ( th3d, qv3d, qc3d, qr3d, qi3d, qs3d,                  &
                   rho3d, pii3d, p3d, dta, z3d,                         &
                   ht, dz3d, con_g, w3d,                                &
                   itimestep, xlat, sdec, land2d,                       &
                   1, nx , 1, 1, 1, lev,                                & ! memory dims
                   1, nxj, 1, 1, 1, lev,                                & ! tile   dims
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
                   SL_sedi, sat_predict, new_saturation, use_cpm, &
                   use_declination, benchmark )

          if ( nmgce3 .eq. 1 )                                          &
          call gsfcgce                                                  &
                 ( th3d, qv3d, qc3d, qr3d, qi3d, qs3d,                  &
                   rho3d, pii3d, p3d, dta, z3d,                         &
                   ht, dz3d, con_g,                                     &
                   rhowater, rhosnow,                                   &
                   itimestep,                                           &
                   1, nx , 1, 1, 1, lev,                                & ! memory dims
                   1, nxj, 1, 1, 1, lev,                                & ! tile   dims
                   rain2d, snow2d, graupel2d, sr2d,                     &
#ifdef EXT_DIAG
                   refl_10cm, diagflag, do_radar_ref,                   &
#endif
#ifdef EffectRad_GCE3
                   land2d,                                              &
                   rew3d, rer3d, rei3d,                                 &
                   res3d, reg3d,                                        &
#endif
                   .false., qg3d,                                       &
                   ihail, ice2 )

!          if ( nmgce3 .eq. 3 )                                          &
!          call gsfcgce_3ice_cwb                                         &
!                 ( th3d, qv3d, qc3d, qr3d, qi3d, qs3d, qg3d,            &
!                   rho3d, pii3d, p3d, dta, z3d,                         &
!                   ht, dz3d, itimestep,                                 &
!                   1, nx , 1, 1, 1, lev,                                & ! memory dims
!                   1, nxj, 1, 1, 1, lev,                                & ! tile   dims
!                   sr2d, rainnc2d, rain2d,                              &
!                   snownc2d, snow2d, graupelnc2d, graupel2d,            &
!                   land2d )

        endif

        if ( nmmiph .eq. 16 )                                           &
          call gsfcgce_4ice_nuwrf                                       &
                 ( th3d, qv3d, qc3d, qr3d, qi3d, qs3d, qh3d, qg3d,      &
                   rho3d, pii3d, p3d, dta, z3d,                         &
                   ht, dz3d, con_g, w3d,                                &
                   rhowater, rhosnow,                                   &
                   itimestep, land2d, dx,                               &
                   1, nx , 1, 1, 1, lev,                                & ! memory dims
                   1, nxj, 1, 1, 1, lev,                                & ! tile   dims
                   rain2d, snow2d, graupel2d, hail2d, sr2d,             &
                   rew3d, rer3d, rei3d,                                 &
                   res3d, reg3d, reh3d,                                 &
#ifdef EXT_DIAG
                   refl_10cm, diagflag, do_radar_ref,                   &
                   physc, physe, physd, physs, physm, physf,            & !simultaneous diabatic heating rate
                   acphysc, acphyse, acphysd, acphyss, acphysm,acphysf, & !accumulated diabatic heating rate
                   preci3d, precs3d, precg3d, prech3d, precr3d,         & !precitation
#endif
                   .false. )

        do k = 1, lev
          kc = lev - k + 1
          do i = 1, nxj
            qt(i,             k) = qv3d(i,kc,1)
            qt(i,(ntcw-1)*lev+k) = qc3d(i,kc,1)
            qt(i,(ntrw-1)*lev+k) = qr3d(i,kc,1)
            qt(i,(ntiw-1)*lev+k) = qi3d(i,kc,1)
            qt(i,(ntsw-1)*lev+k) = qs3d(i,kc,1)
            qt(i,(ntgl-1)*lev+k) = qg3d(i,kc,1)
            tt(i,k) = th3d(i,kc,1)*pk(i,k)  !convert back to real temperature

            if ( nmmiph .eq. 15 ) then
              if ( nmgce3 .eq. 3 ) then
              ! Goddard 3ICE with constant effective radii
                re_cloud(i,k) = 10.0          !micron
                re_rain (i,k) = 1000.0        !micron
                re_ice  (i,k) = 50.0          !micron
                re_snow (i,k) = 250.0         !micron
              else
                re_cloud(i,k) = rew3d(i,k,1)  !micron
                re_rain (i,k) = rer3d(i,k,1)  !micron
                re_ice  (i,k) = rei3d(i,k,1)  !micron
                re_snow (i,k) = res3d(i,k,1)  !micron
              endif
            endif

            if ( nmmiph .eq. 16 ) then
              qt(i,(nthl-1)*lev+k) = qh3d(i,kc,1)
              re_cloud(i,k) = rew3d(i,k,1)  !micron
              re_rain (i,k) = rer3d(i,k,1)  !micron
              re_ice  (i,k) = rei3d(i,k,1)  !micron
              re_snow (i,k) = res3d(i,k,1)  !micron
            endif

          enddo
        enddo
        do i = 1, nxj
          rlsp(i) = rain2d(i,1)  !total large scale precipitation (kg/m^2=mm)
          rlspi(i) = ice2d(i,1)
          rlsps(i) = snow2d(i,1)
          rlspg(i) = graupel2d(i,1)
          if ( nmmiph .eq. 16 ) rlspg(i) = rlspg(i) + hail2d(i,1)
          sr(i)   = sr2d(i,1)
        enddo

        deallocate                                                      &
         ( th3d,qv3d,qc3d,qr3d,qs3d,qi3d,qg3d,pii3d,p3d,z3d,dz3d,rho3d, &
           rain2d,ice2d,snow2d,graupel2d,                               &
           ht,sr2d,land2d,w3d,rew3d,rer3d,rei3d,res3d,reg3d )
#ifdef Readaeroclx
        deallocate ( aero4d )
#endif
#ifdef EXT_DIAG
        deallocate                                                      &
         ( preci3d,precs3d,precg3d,precr3d,prech3d,refl_10cm )
#endif
        if ( nmmiph .eq. 16 ) deallocate                                &
         ( qh3d,hail2d,reh3d )

      endif  !end of if nmmiph=15 .or. nmmiph=16

      return
!--------------------------
      end subroutine mp_scheme
!--------------------------
