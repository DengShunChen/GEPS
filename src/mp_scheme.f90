!--------------------------
      subroutine mp_init                                               &
!--------------------------
!  ---  inputs:
           ( nmmiph ,myrank)
!  ---  outputs: ( none )

! for WSM6
      use module_mp_wsm6,     only : wsm6init
! for Thompson
      use module_mp_thompson, only : thompson_init => thompson_init
! for New Thompson
      use module_mp_thompson_new,                                       &
                              only : new_thompson_init => thompson_init
! for GFDL MP v1
      use module_mp_gfdl,     only : gfdl_cloud_microphys_init
! for GFDL MP v2
      use module_mp_gfdl_v2,  only : gfdlv2_init => gfdl_cld_mp_init
! for GFDL MP v3
      use module_mp_gfdl_v3,  only : gfdlv3_init => gfdl_cld_mp_init
      use physpara, only : is_aerosol_aware,merra2_aerosol_aware

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
! New Thompson
        if ( nmmiph .eq. 9 ) then
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
           ( nmmiph,nx,nxj,lev,ncld,plt,                               &
             pst,dsigma,phii,islimsk,q0,kdt,tpi,me,dta,area,jj,        &
             itimestep,sgeo,phi,rhc,                                   &
!  ---  inputs/outputs:
             tt,qt,qa,ut,vt,vvel,                                      &
!  ---  outputs:
             re_cloud,re_ice,re_snow,re_rain,                          &
             rlsp,sr )

      use rank
      use radn,                only: ntcw,ntiw,ntrw,ntsw,ntgl,nthl,     &
                                     ntinc,ntrnc
! for wsm6
      use module_mp_wsm6,      only: wsm6
! for thompson
      use module_mp_thompson,  only: thompson_driver => mp_gt_driver
! for New Thompson
      use module_mp_thompson_new,                                       &
                               only: new_thompson_driver => mp_gt_driver&
                                     , cal_cldfra3, cfflag_thom
! for GFDL MP v1
      use module_mp_gfdl,      only: gfdl_cloud_microphys_driver,      &
                                     cloud_diagnosis
! for GFDL MP v2
      use module_mp_gfdl_v2,   only: gfdlv2_driver => gfdl_cld_mp_driver
! for GFDL MP v3
      use module_mp_gfdl_v3,   only: gfdlv3_driver => gfdl_cld_mp_driver
! for Goddard (GCE) 4ICE MP
      use module_mp_gce4ice,   only: gsfcgce_4ice_nuwrf
      use physcons,            only: con_rd,con_fvirt,con_g
      use physpara,            only: effr_in
      use const,               only: RTYPE

      implicit none

!  ---  inputs:
      integer,  intent(in)    :: nmmiph,nx,nxj,lev,ncld,kdt,me
!      integer,  intent(in)    :: ntcw,ntrw,ntiw,ntsw,ntgl,ntinc,ntrnc
      integer,  intent(in)    :: islimsk(nx)
      integer,  intent(in)    :: itimestep
      real,     intent(in)    :: tpi,dta,jj
      real,     intent(in)    :: plt(nx,lev),phii(nx,lev+1),phi(nx,lev)
      real,     intent(in)    :: area
      real,     intent(in)    :: rhc(nxj,lev)
      real,     intent(in)    :: sgeo(nx)
      real,     intent(inout) :: vvel(nx,lev) !mb/s
      real(kind=RTYPE), intent(in):: q0(nx,lev*ncld),pst(nx),          &
                                     dsigma(lev,2)
!  ---  inputs/outputs:
      real(kind=RTYPE), intent(inout) :: tt(nx,lev)
      real,     intent(inout) :: qa(nx,lev)
      real(kind=RTYPE), intent(inout) :: ut(nx,lev),vt(nx,lev)
      real(kind=RTYPE), intent(inout):: qt(nx,lev*ncld)
!  ---  outputs:
      real,     intent(inout)   :: re_cloud(nx,lev),re_ice(nx,lev),   &
                                   re_snow(nx,lev),re_rain(nx,lev)
      real,     intent(inout)   :: rlsp(nx),sr(nx)
!  ---  local arrays:
      integer   kc,k,i
      real      prsl(nx,lev),del(nx,lev)
      real      ttc(nx,lev),qtc(nx,lev),qtr(nx,lev),qtrw(nx,lev),      &
                qti(nx,lev),qtsw(nx,lev),qtgl(nx,lev),ntnc(nx,lev,2),  &
                refl10(nx,lev)
      real      rainncv(nx),snowncv(nx),graupelncv(nx)
      real      icem
      logical   lradar
      real,dimension(:),allocatable ::                                  &
              land1d
      real,dimension(:,:),allocatable ::                                &
              qv2d,qc2d,qr2d,qi2d,qs2d,qg2d,qnc2d,qni2d,qnr2d,          &
              rew2d,rer2d,rei2d,res2d,reg2d,                            &
              t2d,dp2d,dz2d,cld2d,w2d,u2d,v2d
      real    qmin, qnmin
! New Thompson MP
      real,dimension(:,:),allocatable ::                                &
              nwfa,nifa,pfils,pflls,vt_dbz_wt
      real,dimension(:),allocatable :: nwfasfc,nifasfc,rainnc,snownc,   &
              icenc,graupelnc,icencv,gridkm
      real,dimension(:,:),allocatable :: rand_pert
      real,dimension(:),allocatable :: spp_prt_list,spp_stddev_cutoff
      character(len=3),dimension(:),allocatable :: spp_var_list
      real    dt_inner
      logical sedi_semi,ext_diag,first_time_step,reset_dBZ,aero_ind_fdb,&
              diagflag
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
                gsize,hs,water1d,rain1d,snow1d,ice1d,graupel1d
      real, dimension(:,:), allocatable ::                              &
                q_con,cappa,te
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
! Goddard (GCE) 4ICE MP
      real    rhowater,rhosnow,dx
      real,dimension(:,:),allocatable ::                                &
              ht,hail2d,rainnc2d,snownc2d,graupelnc2d,hailnc2d,sr2d
      real,dimension(:,:,:),allocatable ::                              &
              th3d,qh3d,rho3d,pii3d,p3d,z3d,rew3d,rer3d,rei3d,res3d,    &
              reg3d,reh3d,refl_10cm
#ifdef EXT_DIAG
      real,dimension(:,:,:),allocatable ::                              &
              physc, physe, physd, physs, physm, physf,                 &
              acphysc, acphyse, acphysd, acphyss, acphysm, acphysf,     &
              preci3d, precs3d, precg3d, prech3d, precr3d
#endif
!
! define rhc for GFDL MP v1 & v2
!     rhc = 1.0                 ! default
!     rhc = 1.0 - 0.02*cosl**2  ! =0.98 at equator; =1.0 at pole
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
             ntnc=0.
             do k=1,lev
               kc=lev-k+1
               do i=1,nxj
                 ntnc(i,kc,1) = qt(i,(ntinc-1)*lev+k)                  &
                      +max(0.,qti(i,kc)-q0(i,(ntiw-1)*lev+k))/icem
                 ntnc(i,kc,2) = qt(i,(ntrnc-1)*lev+k)
               enddo
             enddo

             call thompson_driver(1,nx,1,lev,1,nxj,1,lev,qtc,qtr,qtrw, &
                     qti,qtsw,qtgl,ntnc(1,1,1),ntnc(1,1,2),ttc,        &
                     prsl,del,dta,kdt,rainncv,sr,islimsk,refl10,       &
                     lradar,re_cloud,re_ice,re_snow,me,phii)

             do k=1,lev
               kc=lev-k+1
               do i=1,nxj
                 qt(i,(ntinc-1)*lev+k)=ntnc(i,kc,1)
                 qt(i,(ntrnc-1)*lev+k)=ntnc(i,kc,2)
               enddo
             enddo
           endif
!
        do i=1,nxj
          rlsp(i) = rainncv(i) + snowncv(i) + graupelncv(i)
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

      endif

!     New Thompson
      if ( nmmiph .eq. 9 ) then
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

        qmin=1.0e-12     !minimum of q (kg/kg)
        qnmin=1.0e-6     !minimum of qn (m^-3)

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
            prsl (i,k) = plt(i,kc)*100.                 !layer pressure (Pa)
            t2d  (i,k) = tt(i,kc)
            w2d  (i,k) = - vvel(i,k)*100.*                              &
                         (1.+con_fvirt*qt(i,kc))*tt(i,kc)/              &
                         prsl(i,k)*con_rd/con_g         !vertical velocity (m/s)
            dz2d (i,k) = (phii(i,k+1)-phii(i,k))/con_g  !layer depth (m)
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
            if ( qc2d(i,kc) .lt. qmin ) qc2d(i,kc) = qmin
            if ( qr2d(i,kc) .lt. qmin ) qr2d(i,kc) = qmin
            if ( qi2d(i,kc) .lt. qmin ) qi2d(i,kc) = qmin
            if ( qs2d(i,kc) .lt. qmin ) qs2d(i,kc) = qmin
            if ( qg2d(i,kc) .lt. qmin ) qg2d(i,kc) = qmin
            if ( qni2d(i,kc) .lt. qnmin ) qni2d(i,kc) = qnmin
            if ( qnr2d(i,kc) .lt. qnmin ) qnr2d(i,kc) = qnmin

            qt(i,              k) = qv2d (i,kc)
            qt(i,(ntcw-1) *lev+k) = qc2d (i,kc)
            qt(i,(ntrw-1) *lev+k) = qr2d (i,kc)
            qt(i,(ntiw-1) *lev+k) = qi2d (i,kc)
            qt(i,(ntsw-1) *lev+k) = qs2d (i,kc)
            qt(i,(ntgl-1) *lev+k) = qg2d (i,kc)
            qt(i,(ntinc-1)*lev+k) = qni2d(i,kc)
            qt(i,(ntrnc-1)*lev+k) = qnr2d(i,kc)
            tt(i,              k) = t2d  (i,kc)

            re_cloud(i,k) = rew2d(i,k)*1.E+6   ! m to micron
            re_ice  (i,k) = rei2d(i,k)*1.E+6   ! m to micron
            re_snow (i,k) = res2d(i,k)*1.E+6   ! m to micron
          enddo
        enddo
        do i = 1, nxj
          rlsp(i) = icencv(i)+rainncv(i)+snowncv(i)+graupelncv(i)  !total large scale precipitation (mm)
        enddo

        deallocate                                                      &
         ( t2d,qv2d,qc2d,qr2d,qi2d,qs2d,qg2d,qni2d,qnr2d,qnc2d,         &
           nwfa,nifa,dz2d,w2d,pfils,pflls,vt_dbz_wt,nwfasfc,nifasfc,    &
           rew2d,rei2d,res2d,rainnc,snownc,icenc,graupelnc,icencv,      &
           rand_pert,spp_prt_list,spp_stddev_cutoff,spp_var_list )
        if ( cfflag_thom .eq. 2 ) deallocate ( cld2d,land1d,gridkm )
      endif  !end of if nmmiph=9

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
           ice2d(nxj,1),graupel2d(nxj,1),garea(nxj,1) )
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
            if ( effr_in ) dp2d(i,k) = dp3d(i,1,k)
          enddo
        enddo

        call gfdl_cloud_microphys_driver                                &
                ( qv3d, qc3d, qr3d, qi3d, qs3d, qg3d, cld3d, qnc3d,     &
                  qvten3d, qcten3d, qrten3d, qiten3d, qsten3d, qgten3d, &
                  cldten3d, tten3d, t3d, w3d, u3d, v3d, uten3d, vten3d, &
                  dz3d, dp3d, garea, dta, land2d,                       &
                  rain2d, snow2d, ice2d, graupel2d,                     &
                  rhc, hydrostatic, phys_hydrostatic,               &
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
           land2d,rain2d,snow2d,ice2d,graupel2d,garea )
        if ( effr_in ) deallocate                                       &
           ( dp2d,rho2d,qc2d,qr2d,qi2d,qs2d,qg2d,mask1d,t2d )
      endif  ! end of nmmiph.eq.11

!     GFDL MP v2 & v3
      if ( nmmiph .eq. 12 .or. nmmiph .eq. 13) then
        allocate                                                        &
         ( qv2d(nxj,lev),qc2d(nxj,lev),qr2d(nxj,lev),qi2d(nxj,lev),     &
           qs2d(nxj,lev),qg2d(nxj,lev),cld2d(nxj,lev),qnc2d(nxj,lev),   &
           qni2d(nxj,lev),w2d(nxj,lev),t2d(nxj,lev),u2d(nxj,lev),       &
           v2d(nxj,lev),dp2d(nxj,lev),dz2d(nxj,lev),                    &
           q_con(nxj,lev),cappa(nxj,lev),te(nxj,lev),                   &
           hs(nxj),land1d(nxj),gsize(nxj),rain1d(nxj),snow1d(nxj),      &
           ice1d(nxj),graupel1d(nxj),water1d(nxj) )
#ifdef EXT_DIAG
        allocate                                                        &
         ( prefluxr(nxj,lev),prefluxi(nxj,lev),prefluxs(nxj,lev),       &
           prefluxg(nxj,lev),prefluxw(nxj,lev) )
        allocate                                                        &
         ( cond0(nxj),dep0(nxj),evap0(nxj),sub0(nxj) )
#endif

        hydrostatic = .false.       !flag for hydrostatic solver
        phys_hydrostatic = .true.   !flag for hydrostatic heating from physics 
        sedi_w = .false.
        consv_te = .false.          !flag for energy conservation
        last_step = .true.          !flag for final clean-up (not sure)
        do_inline_mp = .false.      !flag for inline GFDLMP

        te    = 0.0
        q_con = 0.0  !not sure
        cappa = 0.0  !not sure
        gsize = 0.0
        land1d = 0.0
        rain1d = 0.
        snow1d = 0.
        ice1d  = 0.
        graupel1d = 0.
        water1d = 0.
#ifdef EXT_DIAG
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
            if ( sedi_w ) then
              prsl(i,k) = 100.0 * plt(i,k)            !layer mean pressure (from mb to Pa)
              w2d(i,k)  = -vvel(i,k)*(1.+con_fvirt*qt(i,k))*tt(i,k)     &
                         /prsl(i,k)*con_rd/con_g      !vertical velocity (m/s)
            else
              w2d(i,k)  = 0.
            endif

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
                  1, nxj, 1, lev, q_con, cappa, consv_te, te,           &
#ifdef EXT_DIAG
                  prefluxr, prefluxi, prefluxs, prefluxg,               &
                  cond0, dep0, evap0, sub0,                             &
#endif
                  rhc, last_step, do_inline_mp )

        ! GFDL MP v3
        if ( nmmiph .eq. 13 )                                           &
          call gfdlv3_driver                                            &
                ( qv2d, qc2d, qr2d, qi2d, qs2d, qg2d,                   &
                  cld2d, qnc2d, qni2d,                                  &
                  t2d, w2d, u2d, v2d, dz2d, dp2d, gsize, dta, hs,       &
                  land1d, water1d,                                      &
                  rain1d, snow1d, ice1d, graupel1d, hydrostatic,        &
                  1, nxj, 1, lev, q_con, cappa, consv_te, te,           &
#ifdef EXT_DIAG
                  prefluxw, prefluxr, prefluxi, prefluxs, prefluxg,     &
                  cond0, dep0, evap0, sub0,                             &
#endif
                  rhc, last_step, do_inline_mp )

        qmin = 1.0e-15     !minimum of q (kg/kg)
        do k = 1, lev
          do i = 1, nxj
            qt(i,             k) = qv2d(i,k)
            qt(i,(ntcw-1)*lev+k) = max( qc2d(i,k) , qmin )
            qt(i,(ntrw-1)*lev+k) = max( qr2d(i,k) , qmin )
            qt(i,(ntiw-1)*lev+k) = max( qi2d(i,k) , qmin )
            qt(i,(ntsw-1)*lev+k) = max( qs2d(i,k) , qmin )
            qt(i,(ntgl-1)*lev+k) = max( qg2d(i,k) , qmin )
            qa(i,k)  = cld2d(i,k)
            tt(i,k)  = t2d  (i,k)
            ut(i,k)  = u2d  (i,k)
            vt(i,k)  = v2d  (i,k)

            if ( sedi_w ) then
              vvel(i,k)  = -w2d(i,k)*prsl(i,k)*con_g/con_rd             &
                           /((1+con_fvirt*qt(i,k))*tt(i,k))
            endif
          enddo
        enddo

        do i = 1, nxj
          rlsp(i) = water1d(i)+rain1d(i)+snow1d(i)+ice1d(i)+graupel1d(i)     !total large scale precipitation (mm)
          if ( rlsp(i) .gt. rainmin ) then
            sr(i) = (snow1d(i)+ice1d(i)+graupel1d(i))                   &
                    /(water1d(i)+rain1d(i)+snow1d(i)+ice1d(i)+graupel1d(i))  !snow ratio
          else
            sr(i) = 0.0
          endif
        enddo

        deallocate                                                      &
         ( qv2d,qc2d,qr2d,qi2d,qs2d,qg2d,qnc2d,qni2d,cld2d,w2d,t2d,u2d, &
           v2d,dp2d,dz2d,q_con,cappa,te,hs,gsize,rain1d,snow1d,ice1d,   &
           graupel1d,water1d,land1d )
#ifdef EXT_DIAG
        deallocate                                                      &
         ( prefluxw,prefluxr,prefluxi,prefluxs,prefluxg,cond0,dep0,     &
           evap0,sub0 )
#endif
      endif  !end if nmmiph.eq.12 .or nmmiph.eq.13

!     Goddard (GCE) 4ICE MP
      if ( nmmiph .eq. 16 ) then
        allocate                                                        &
         ( th3d(nx,lev,1),qv3d(nx,lev,1),qc3d(nx,lev,1),qr3d(nx,lev,1), &
           qi3d(nx,lev,1),qs3d(nx,lev,1),qg3d(nx,lev,1),qh3d(nx,lev,1), &
           rho3d(nx,lev,1),pii3d(nx,lev,1),p3d(nx,lev,1),z3d(nx,lev,1), &
           ht(nx,1),dz3d(nx,1,lev),w3d(nx,1,lev),rainnc2d(nx,1),        &
           snownc2d(nx,1),graupelnc2d(nx,1),hailnc2d(nx,1),rain2d(nx,1),&
           snow2d(nx,1),graupel2d(nx,1),hail2d(nx,1),sr2d(nx,1),        &
           rew3d(nx,lev,1),rer3d(nx,lev,1),rei3d(nx,lev,1),             &
           res3d(nx,lev,1),reg3d(nx,lev,1),reh3d(nx,lev,1),             &
           land2d(nx,1),refl_10cm(nx,lev,1) )

        th3d = 0.
        qv3d = 0.
        qc3d = 0.
        qr3d = 0.
        qi3d = 0.
        qs3d = 0.
        qg3d = 0.
        qh3d = 0.
        rho3d = 0.
        pii3d = 0.
        p3d = 0.
        z3d = 0.
        ht = 0.
        dz3d = 0.
        w3d = 0.
        rain2d = 0.
        snow2d = 0.
        graupel2d = 0.
        hail2d = 0.
        rainnc2d = 0.
        snownc2d = 0.
        graupelnc2d = 0.
        hailnc2d = 0.
        sr2d = 0.
        rew3d = 0.
        rer3d = 0.
        rei3d = 0.
        res3d = 0.
        reg3d = 0.
        reh3d = 0.
        land2d = 0.
        refl_10cm = 0.

        diagflag = .false.          !if diagflag=true and do_radar_ref=1, call calc_refl10cm
        do_radar_ref = 0

        dx = sqrt(area)             !grid length (m)
        rhowater = 1000.            !water density (kg/m^3), but not used
        rhosnow = 100.              !snow density (kg/m^3), but not used

        do i = 1, nxj
          ht(i,1) = sgeo(i)/con_g   !terrain geopotential height above sea level (m)
          if( islimsk(i) .eq. 1 ) then
            land2d(i,1) = 1.        !land
          else
            land2d(i,1) = 2.        !ocean & seaice
          endif
        enddo

        do k = 1, lev
          kc = lev - k + 1
          do i = 1, nxj
            qv3d (i,k,1) = qt(i,             kc)
            qc3d (i,k,1) = qt(i,(ntcw-1)*lev+kc)
            qr3d (i,k,1) = qt(i,(ntrw-1)*lev+kc)
            qi3d (i,k,1) = qt(i,(ntiw-1)*lev+kc)
            qs3d (i,k,1) = qt(i,(ntsw-1)*lev+kc)
            qg3d (i,k,1) = qt(i,(ntgl-1)*lev+kc)
            qh3d (i,k,1) = qt(i,(nthl-1)*lev+kc)
            p3d  (i,k,1) = 100.0*plt(i,kc)                   !layer mean pressure (from mb to Pa)
!            pii3d(i,k,1) = pk(i,kc)                          !exner function, =(p/psfc)**(Rd/cp)
            pii3d(i,k,1) = 1.
!            th3d (i,k,1) = tt(i,kc)*pk(i,kc)                 !potential temperature (K)
            th3d (i,k,1) = tt(i,kc)                          !temperature (K)
            z3d  (i,k,1) = phi(i,kc)/con_g                   !layer geopotential height above sea level (m)
            dz3d (i,k,1) = (phii(i,k+1)-phii(i,k))/con_g     !layer thickness (m)
            ! use virtural temperature : Tv = (1+(Rv/Rd-1)*q)*T = (1+con_fvirt*q)*T
            rho3d(i,k,1) = p3d(i,k,1)/con_rd/tt(i,kc)                   &
                          /(1+con_fvirt*qt(i,kc))            !density of air (kg/m^3)
            w3d  (i,k,1) = -vvel(i,k)*100.*(1.+con_fvirt*qt(i,kc))      &
                          *tt(i,kc)/p3d(i,k,1)*con_rd/con_g  !vertical velocity (m/s)
          enddo
        enddo

        call gsfcgce_4ice_nuwrf                                         &
                 ( th3d, qv3d, qc3d, qr3d, qi3d, qs3d, qh3d, qg3d,      &
                   rho3d, pii3d, p3d, dta, z3d,                         &
                   ht, dz3d, con_g, w3d,                                &
                   rhowater, rhosnow,                                   &
                   itimestep, land2d, dx,                               &
!                   ids,ide, jds,jde, kds,kde,                           & ! domain dims
                   1, nx , 1, 1, 1, lev,                                & ! memory dims
                   1, nxj, 1, 1, 1, lev,                                & ! tile   dims
                   rainnc2d, rain2d,                                    &
                   snownc2d, snow2d, sr2d,                              &
                   graupelnc2d, graupel2d,                              &
                   hailnc2d, hail2d,                                    &
                   refl_10cm, diagflag, do_radar_ref,                   &
                   rew3d, rer3d, rei3d,                                 &
                   res3d, reg3d, reh3d,                                 &
#ifdef EXT_DIAG
                   physc, physe, physd, physs, physm, physf,            &
                   acphysc, acphyse, acphysd, acphyss, acphysm, acphysf,&
                   preci3d, precs3d, precg3d, prech3d, precr3d,         &
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
            qt(i,(nthl-1)*lev+k) = qh3d(i,kc,1)
            tt(i,k) = th3d(i,kc,1)

            re_cloud(i,k) = rew3d(i,k,1)  !micron
            re_rain (i,k) = rer3d(i,k,1)  !micron
            re_ice  (i,k) = rei3d(i,k,1)  !micron
            re_snow (i,k) = res3d(i,k,1)  !micron
          enddo
        enddo
        do i = 1, nxj
          rlsp(i) = rain2d(i,1)+snow2d(i,1)+graupel2d(i,1)+hail2d(i,1)  !total large scale precipitation (mm)
          sr(i)   = sr2d(i,1)
        enddo

        deallocate                                                      &
         ( th3d,qv3d,qc3d,qr3d,qs3d,qi3d,qg3d,qh3d,rho3d,pii3d,p3d,z3d, &
           ht,dz3d,w3d,rainnc2d,snownc2d,graupelnc2d,hailnc2d,rain2d,   &
           snow2d,graupel2d,hail2d,sr2d,rew3d,rer3d,rei3d,res3d,reg3d,  &
           reh3d,land2d,refl_10cm )
      endif  !end of if nmmiph=16

      return
!--------------------------
      end subroutine mp_scheme
!--------------------------
