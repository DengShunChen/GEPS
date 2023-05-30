!--------------------------
      subroutine mp_init                                               &
!--------------------------
!  ---  inputs:
           ( nmmiph ,myrank)
!  ---  outputs: ( none )

! for WSM6
      use module_mp_wsm6,     only : wsm6init
! for Thompson
      use module_mp_thompson, only : thompson_init
! for GFDLMP
#if defined (GFDLMP_v2)
      use module_mp_gfdl_v2,  only : gfdl_cld_mp_init
#else
      use module_mp_gfdl,     only : gfdl_cloud_microphys_init
#endif

      implicit none
!  ---  input:
      integer,  intent(in) :: nmmiph,myrank
!  ---  local:
      integer   ntrac_req
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
! GFDLMP
        if ( nmmiph .eq. 11 ) then
#if defined (GFDLMP_v2)
          call gfdl_cld_mp_init()
#else
          call gfdl_cloud_microphys_init()
#endif
          if ( myrank .eq. 0 )                                         &
             print *, 'GFDL cloud microphysics initialized'
        endif

      return
!--------------------------
      end subroutine mp_init
!--------------------------
!--------------------------
      subroutine mp_scheme                                             &
!--------------------------
!  ---  inputs:
           ( nmmiph,nx,nxj,lev,ncld,plt,&!prsi,                          &
             pst,dsigma,phii,islimsk,q0,kdt,ntcw,ntrw,ntiw,ntsw,       &
             ntgl,ntinc,ntrnc,tpi,me,dta,area,jj,                      &
#if defined (GFDLMP_v2)
             sgeo,                                                     &
#endif
!  ---  inputs/outputs:
             tt,qt,qa,ut,vt,sd,                                        &
!  ---  outputs:
             re_cloud,re_ice,re_snow,re_rain,                          &
             rlsp,sr )

! for wsm6
      use module_mp_wsm6,      only: wsm6
! for thompson
      use module_mp_thompson,  only: mp_gt_driver
! for GFDLMP
#if defined (GFDLMP_v2)
      use module_mp_gfdl_v2,   only: gfdl_cld_mp_driver
#else
      use module_mp_gfdl,      only: gfdl_cloud_microphys_driver,      &
                                     cloud_diagnosis
#endif
      use physcons,            only: con_rd,con_fvirt,con_g
      use physpara,            only: effr_in
      use const,               only: RTYPE

      implicit none

!  ---  inputs:
      integer,  intent(in)    :: nmmiph,nx,nxj,lev,ncld,kdt,me
      integer,  intent(in)    :: ntcw,ntrw,ntiw,ntsw,ntgl,ntinc,ntrnc
      integer,  intent(in)    :: islimsk(nx)
      real,     intent(in)    :: tpi,dta,jj
      real,     intent(in)    :: plt(nx,lev),phii(nx,lev+1)!,       &
!                                 prsi(nx,lev+1)
      real,     intent(in)    :: area(nx)  ! area of grid box (m^2)
#if defined (GFDLMP_v2)
      real,     intent(in)    :: sgeo(nx)
#endif
!      real,     intent(in)    :: sd(nx,lev+1)
      real(kind=RTYPE), intent(inout) :: sd(nx,lev+1)
      real(kind=RTYPE), intent(in):: q0(nx,lev*ncld),pst(nx),          &
                                     dsigma(lev,2)
!  ---  inputs/outputs:
      real(kind=RTYPE), intent(inout) :: tt(nx,lev)
      real,     intent(inout) :: qa(nx,lev)  ! only changed in GFDL MP
      real(kind=RTYPE), intent(inout) :: ut(nx,lev),vt(nx,lev)
      real(kind=RTYPE), intent(inout):: qt(nx,lev*ncld)
!  ---  outputs:
      real,     intent(inout)   :: re_cloud(nx,lev),re_ice(nx,lev),    &
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
! GFDLMP
      real, parameter ::                                                &
                rainmin=1.0e-10 !(mm)
      real, dimension(:,:), allocatable ::                              &
                dot,dp
      logical   hydrostatic,phys_hydrostatic,sedi_w 
#if defined (GFDLMP_v2)
      real, dimension(:,:), allocatable ::                              &
                qv,ql,qr,qi,qs,qg,cldcov,qnl,qni,w,delp,dz,q_con,cappa, &
                te,pt,uin,vin,prefluxr,prefluxi, prefluxs,prefluxg
      real, dimension(:), allocatable ::                                &
                gsize,hs,rain0,snow0,ice0,graupel0,                     &
                cond0,dep0,evap0,sub0
      logical   consv_te,last_step,do_inline_mp
#else
      real, dimension(:,:), allocatable ::                              &
                garea,frland,rain0,snow0,ice0,graupel0
      real, dimension(:,:,:), allocatable ::                            &
                qv1,ql1,qr1,qi1,qs1,qg1,qa1,qn1,pt,w,uin,vin,delp,dz,   &
                qv_dt,ql_dt,qr_dt,qi_dt,qs_dt,qg_dt,qa_dt,udt,vdt,pt_dt
      real, dimension(:,:), allocatable ::                              &
                ql2,qr2,qi2,qs2,qg2,rho,                                &
                re_graupel,rew,rei,rer,res,reg
#endif
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

      if ( nmmiph .eq. 11 ) then
#if defined (GFDLMP_v2)
        allocate                                                        &
         ( qv(nxj,lev),ql(nxj,lev),qr(nxj,lev),qi(nxj,lev),qs(nxj,lev), &
           qg(nxj,lev),cldcov(nxj,lev),qnl(nxj,lev),qni(nxj,lev),       &
           w(nxj,lev),pt(nxj,lev),uin(nxj,lev),vin(nxj,lev),            &
           delp(nxj,lev),dz(nxj,lev),prefluxr(nxj,lev),                 &
           prefluxi(nxj,lev),prefluxs(nxj,lev),prefluxg(nxj,lev),       &
           q_con(nxj,lev),cappa(nxj,lev),te(nxj,lev) )
        allocate                                                        &
         ( hs(nxj),gsize(nxj),rain0(nxj),snow0(nxj),ice0(nxj),          &
           graupel0(nxj),cond0(nxj),dep0(nxj),evap0(nxj),sub0(nxj) )
#else
        allocate                                                        &
         ( re_graupel(nxj,lev),rew(nxj,lev),                            &
           rei(nxj,lev),rer(nxj,lev),res(nxj,lev),reg(nxj,lev),         &
           frland(nxj,1),rain0(nxj,1),snow0(nxj,1),ice0(nxj,1),         &
           graupel0(nxj,1),garea(nxj,1),                                &
           qv1(nxj,1,lev),ql1(nxj,1,lev),qr1(nxj,1,lev),qi1(nxj,1,lev), &
           qs1(nxj,1,lev),qg1(nxj,1,lev),qa1(nxj,1,lev),qn1(nxj,1,lev), &
           pt(nxj,1,lev),w(nxj,1,lev),uin(nxj,1,lev),vin(nxj,1,lev),    &
           delp(nxj,1,lev),dz(nxj,1,lev),                               &
           qv_dt(nxj,1,lev),ql_dt(nxj,1,lev),qr_dt(nxj,1,lev),          &
           qi_dt(nxj,1,lev),qs_dt(nxj,1,lev),qg_dt(nxj,1,lev),          &
           qa_dt(nxj,1,lev),udt(nxj,1,lev),vdt(nxj,1,lev),              &
           pt_dt(nxj,1,lev) )
        if ( effr_in ) allocate                                         &
           ( dp(nxj,lev),rho(nxj,lev),ql2(nxj,lev),qr2(nxj,lev),        &
             qi2(nxj,lev),qs2(nxj,lev),qg2(nxj,lev) )
#endif
        if ( sedi_w ) allocate ( dot(nxj,lev) )
#if defined (GFDLMP_v2)
        qnl   = 0.0
        qni   = 0.0
        pt    = 0.0
        uin   = 0.0
        vin   = 0.0
        te    = 0.0
        gsize = 0.0
        cond0 = 0.0
        dep0  = 0.0
        evap0 = 0.0
        sub0  = 0.0
        q_con = 0.0  !not sure
        cappa = 0.0  !not sure
        cldcov    = 0.0  !for do_qa=.false.
        prefluxr  = 0.0
        prefluxi  = 0.0
        prefluxs  = 0.0
        prefluxg  = 0.0
#else
        frland = 0.
        garea = 0.
        qv_dt = 0.
        ql_dt = 0.
        qr_dt = 0.
        qi_dt = 0.
        qs_dt = 0.
        qg_dt = 0.
        qa_dt = 0.
        pt_dt = 0.
        udt   = 0.
        vdt   = 0.
#endif
        rain0 = 0.
        snow0 = 0.
        ice0  = 0.
        graupel0 = 0.
      endif
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
             do k=1,lev
               kc=lev-k+1
               do i=1,nxj
                 ntnc(i,kc,1) = qt(i,(ntinc-1)*lev+k)                  &
                      +max(0.,qti(i,kc)-q0(i,(ntiw-1)*lev+k))/icem
                 ntnc(i,kc,2) = qt(i,(ntrnc-1)*lev+k)
               enddo
             enddo
             call mp_gt_driver(1,nx,1,lev,1,nxj,1,lev,qtc,qtr,qtrw,    &
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

!     GFDLMP
      if ( nmmiph .eq. 11 ) then
        hydrostatic = .false.       !flag for hydrostatic solver
        phys_hydrostatic = .true.   !flag for hydrostatic heating from physics 
        sedi_w = .false.

#if defined (GFDLMP_v2)
        consv_te = .false.          !flag for energy conservation
        last_step = .true.          !flag for final clean-up (not sure)
        do_inline_mp = .false.      !flag for inline GFDLMP

        do i = 1, nxj
          gsize(i) = sqrt(area(i))  !square root of grid area (m)
          hs(i)    = sgeo(i)        !terrain geopotential (gpm) (not sure)
        enddo
        do k = 1, lev
          kc = lev - k + 1
          do i = 1, nxj
            if ( sedi_w ) then
              prsl(i,k) = 100.0 * plt(i,k)            !layer mean pressure (from mb to Pa)
              dot(i,k) = 0.5*(sd(i,k)+sd(i,k+1))*100. !vertical velocity (from mb/s to Pa/s)
              w(i,k)   = -dot(i,k)*(1.+con_fvirt*qt(i,k))*tt(i,k)      &
                         /prsl(i,k)*con_rd/con_g      !vertical velocity (m/s)
            else
              w(i,k)   = 0.
            endif

            qv(i,k)   = qt(i,             k)
            ql(i,k)   = qt(i,(ntcw-1)*lev+k)
            qr(i,k)   = qt(i,(ntrw-1)*lev+k)
            qi(i,k)   = qt(i,(ntiw-1)*lev+k)
            qs(i,k)   = qt(i,(ntsw-1)*lev+k)
            qg(i,k)   = qt(i,(ntgl-1)*lev+k)
            pt(i,k)   = tt(i,k)                 !temperature
            uin(i,k)  = ut(i,k)                 !zonal wind (m/s)
            vin(i,k)  = vt(i,k)                 !meridional wind (m/s)
            delp(i,k) = (dsigma(k,1)*pst(i)+dsigma(k,2))*100. !difference of interface pressure (Pa)
            dz(i,k)   = (phii(i,kc)-phii(i,kc+1))/con_g       !differences of height (m), dz<0
          enddo
        enddo
          
        call gfdl_cld_mp_driver                                         &
                ( qv, ql, qr, qi, qs, qg, cldcov, qnl, qni,             &
                  pt, w, uin, vin, dz, delp, gsize, dta, hs,            &
                  rain0, snow0, ice0, graupel0, hydrostatic,            &
                  1, nxj, 1, lev, q_con, cappa, consv_te, te,           &
                  prefluxr, prefluxi, prefluxs, prefluxg,               &
                  cond0, dep0, evap0, sub0,                             &
                  last_step, do_inline_mp )

        do k = 1, lev
          do i = 1, nxj
            qt(i,             k) = qv(i,k)
            qt(i,(ntcw-1)*lev+k) = ql(i,k)
            qt(i,(ntrw-1)*lev+k) = qr(i,k)
            qt(i,(ntiw-1)*lev+k) = qi(i,k)
            qt(i,(ntsw-1)*lev+k) = qs(i,k)
            qt(i,(ntgl-1)*lev+k) = qg(i,k)
            qa(i,k)  = cldcov(i,k)
            tt(i,k)  = pt (i,k)
            ut(i,k)  = uin(i,k)
            vt(i,k)  = vin(i,k)

            if ( sedi_w ) then
              dot(i,k)  = -w(i,k)*prsl(i,k)*con_g/con_rd                &
                          /((1+con_fvirt*qt(i,k))*tt(i,k))
              sd(i,k+1) = dot(i,k)/100./0.5-sd(i,k)
            endif
          enddo
        enddo

        do i = 1, nxj
          rain0(i)    = rain0(i)
          ice0(i)     = ice0(i)
          snow0(i)    = snow0(i)
          graupel0(i) = graupel0(i)

          rlsp(i) = rain0(i)+snow0(i)+ice0(i)+graupel0(i)     !total large scale precipitation (mm)
          if ( rlsp(i) .gt. rainmin ) then
            sr(i) = (snow0(i)+ice0(i)+graupel0(i))                      &
                    /(rain0(i)+snow0(i)+ice0(i)+graupel0(i))  !snow ratio
          else
            sr(i) = 0.0
          endif
        enddo

        deallocate                                                      &
         ( qv,ql,qr,qi,qs,qg,qnl,qni,cldcov,w,pt,uin,vin,delp,dz,       &
           prefluxr,prefluxi,prefluxs,prefluxg,q_con,cappa,te)
        deallocate                                                      &
         ( hs,gsize,rain0,snow0,ice0,graupel0,cond0,dep0,evap0,sub0 )
#else
        do i = 1, nxj
          if( islimsk(i) .eq. 1 ) frland(i,1) = 1.  !land fraction
          garea(i,1) = area(i)                      !area of grid box (m^-2)
        enddo
         
        do k = 1, lev
          kc = lev - k + 1
          do i = 1, nxj
            if ( sedi_w ) then
              prsl(i,k) = 100.0 * plt(i,k)             !layer mean pressure (from mb to Pa)
              dot(i,k)  = 0.5*(sd(i,k)+sd(i,k+1))*100. !vertical velocity (from mb/s to Pa/s)
              w(i,1,k)  = -dot(i,k)*(1.+con_fvirt*qt(i,k))*tt(i,k)      &
                          /prsl(i,k)*con_rd/con_g      !vertical velocity (m/s)
            else
              w(i,1,k)  = 0.
            endif

            qv1(i,1,k)  = qt(i,             k)
            ql1(i,1,k)  = qt(i,(ntcw-1)*lev+k)
            qr1(i,1,k)  = qt(i,(ntrw-1)*lev+k)
            qi1(i,1,k)  = qt(i,(ntiw-1)*lev+k)
            qs1(i,1,k)  = qt(i,(ntsw-1)*lev+k)
            qg1(i,1,k)  = qt(i,(ntgl-1)*lev+k)
            qn1(i,1,k)  = 0.                      ! =0. for prog_ccn=.false. (cm^-3)
            qa1(i,1,k)  = 0.                      !layer cloud fraction (should set to zero for do_qa=.false.)
            pt(i,1,k)   = tt(i,k)                 !temperature
            uin(i,1,k)  = ut(i,k)                 !zonal wind (m/s)
            vin(i,1,k)  = vt(i,k)                 !meridional wind (m/s)
            delp(i,1,k) = (dsigma(k,1)*pst(i)+dsigma(k,2))*100. !difference of interface pressure (Pa)
            dz(i,1,k)   = (phii(i,kc)-phii(i,kc+1))/con_g !differences of height (m), dz<0
            if ( effr_in ) dp(i,k) = delp(i,1,k)
          enddo
        enddo

        call gfdl_cloud_microphys_driver                                &
                ( qv1, ql1, qr1, qi1, qs1, qg1, qa1, qn1,               &
                  qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt, qa_dt,      &
                  pt_dt, pt, w, uin, vin, udt, vdt,                     &
                  dz, delp, garea, dta, frland,                         &
                  rain0, snow0, ice0, graupel0,                         &
                  hydrostatic, phys_hydrostatic,                        &
                  1, nxj, 1, 1, 1, lev, 1, lev )

        do k = 1, lev
          do i = 1, nxj
            ql2(i,k) = ql1(i,1,k) + ql_dt(i,1,k) * dta
            qr2(i,k) = qr1(i,1,k) + qr_dt(i,1,k) * dta
            qi2(i,k) = qi1(i,1,k) + qi_dt(i,1,k) * dta
            qs2(i,k) = qs1(i,1,k) + qs_dt(i,1,k) * dta
            qg2(i,k) = qg1(i,1,k) + qg_dt(i,1,k) * dta

            qa(i,k)  = qa1(i,1,k) + qa_dt(i,1,k) * dta
            qt(i,             k) = qv1(i,1,k) + qv_dt(i,1,k) * dta
            qt(i,(ntcw-1)*lev+k) = ql2(i,k)
            qt(i,(ntrw-1)*lev+k) = qr2(i,k)
            qt(i,(ntiw-1)*lev+k) = qi2(i,k)
            qt(i,(ntsw-1)*lev+k) = qs2(i,k)
            qt(i,(ntgl-1)*lev+k) = qg2(i,k)
            tt(i,k)  = pt(i,1,k)  + pt_dt(i,1,k) * dta
            ut(i,k)  = uin(i,1,k) + udt(i,1,k)   * dta
            vt(i,k)  = vin(i,1,k) + vdt(i,1,k)   * dta

            if ( sedi_w ) then
              dot(i,k) = -w(i,1,k)*prsl(i,k)*con_g/con_rd               &
                          /((1+con_fvirt*qt(i,k))*tt(i,k))
              sd(i,k+1)= dot(i,k)/100./0.5-sd(i,k)
            endif

            if ( effr_in ) then
              rho(i,k) = 0.622*prsl(i,k)                                &
                         /( con_rd*tt(i,k)*(qt(i,k)+0.622) ) !air density (kg/m^3)
            endif
          enddo
        enddo

        if ( effr_in ) then
          call cloud_diagnosis                                          &
!               ( 1, nx, 1, lev, rho, qtr, qti, qtrw, qtsw, qtgl, tt,    &  ! module_mp_gfdl_fv3.f90
!                 rew, rei, rer, res, reg )
               ( 1, nxj, 1, lev, rho, dp, islimsk,                      &  ! module_mp_gfdl_fv3_v16.f90
                 ql2, qi2, qr2, qs2, qg2, tt,                           &
                 rew, rei, rer, res, reg )
          do k = 1, lev
            kc = lev - k + 1
            do i = 1, nxj
              re_cloud(i,k)   = rew(i,kc)   !(micron)
              re_ice(i,k)     = rei(i,kc)   !(micron)
              re_rain(i,k)    = rer(i,kc)   !(micron)
              re_snow(i,k)    = res(i,kc)   !(micron)
              re_graupel(i,k) = reg(i,kc)   !(micron)
            enddo
          enddo
        endif
!
        do i = 1, nxj
          if ( rain0(i,1)    .lt. rainmin ) rain0(i,1)    = 0.0
          if ( ice0(i,1)     .lt. rainmin ) ice0(i,1)     = 0.0
          if ( snow0(i,1)    .lt. rainmin ) snow0(i,1)    = 0.0
          if ( graupel0(i,1) .lt. rainmin ) graupel0(i,1) = 0.0

          rlsp(i) = rain0(i,1)+snow0(i,1)+ice0(i,1)+graupel0(i,1)  !total large scale precipitation (mm)
          if ( rlsp(i) .gt. rainmin ) then
            sr(i) = (snow0(i,1)+ice0(i,1)+graupel0(i,1))                &
                    /(rain0(i,1)+snow0(i,1)+ice0(i,1)+graupel0(i,1))  !snow ratio
          else
            sr(i) = 0.0
          endif
        enddo

        deallocate                                                      &
          ( re_graupel,rew,rei,rer,res,reg,                             &
            frland,rain0,snow0,ice0,graupel0,garea,                     &
            qv1,ql1,qr1,qi1,qs1,qg1,qa1,qn1,pt,w,uin,vin,delp,dz,       &
            qv_dt,ql_dt,qr_dt,qi_dt,qs_dt,qg_dt,qa_dt,udt,vdt,pt_dt )
        if ( effr_in ) deallocate ( dp,rho,ql2,qi2,qr2,qs2,qg2 )
#endif
        if ( sedi_w ) deallocate ( dot )

      endif  ! end of nmmiph.eq.11

      return
!--------------------------
      end subroutine mp_scheme
!--------------------------
