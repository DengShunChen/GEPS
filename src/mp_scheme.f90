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
      use module_mp_gfdl,     only : gfdl_cloud_microphys_init

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
          call gfdl_cloud_microphys_init()
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
      use module_mp_gfdl,      only: gfdl_cloud_microphys_driver,      &
                                     cloud_diagnosis
      use physcons,            only: con_rd,con_fvirt,con_g
      use physpara,            only: effr_in
      use const,               only: RTYPE

      implicit none

!  ---  inputs:
      integer,  intent(in)    :: nmmiph,nx,nxj,lev,ncld,kdt,me
      integer,  intent(in)    :: ntcw,ntrw,ntiw,ntsw,ntgl,ntinc,ntrnc
      integer,  intent(in)    :: islimsk(nx)
      real,     intent(in)    :: tpi,dta,jj
      real,     intent(in)    :: plt(nx,lev),dsigma(lev,2),            &
                                 phii(nx,lev+1)!,       &
!                                 prsi(nx,lev+1)
      real,     intent(in)    :: area(nx,1)  ! area of grid box (m^2)
!      real,     intent(in)    :: sd(nx,lev+1)
      real,     intent(inout) :: sd(nx,lev+1)
      real(kind=RTYPE), intent(in):: q0(nx,lev*ncld),pst(nx)
!  ---  inputs/outputs:
      real,     intent(inout) :: tt(nx,lev)
      real,     intent(inout) :: qa(nx,lev)  ! only changed in GFDL MP
      real,     intent(inout) :: ut(nx,lev),vt(nx,lev)
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
!                rainmin=1.0e-10  !(mm?)
                rainmin=1.0e-20  !test
      real, dimension(:,:), allocatable ::                              &
                dot,rho,re_graupel,rew,rei,rer,res,reg,dp
      real, dimension(:,:), allocatable ::                              &
                frland,rain0,snow0,ice0,graupel0
      real, dimension(:,:,:), allocatable ::                            &
                qv1,ql1,qr1,qi1,qs1,qg1,qa1,qn1,pt,w,uin,vin,delp,dz,   &
                qv_dt,ql_dt,qr_dt,qi_dt,qs_dt,qg_dt,qa_dt,udt,vdt,pt_dt
      logical   hydrostatic,phys_hydrostatic,sedi_w 
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
        allocate                                                        &
         ( re_graupel(nx,lev),rew(nx,lev),                              &
           rei(nx,lev),rer(nx,lev),res(nx,lev),reg(nx,lev),             &
           frland(nx,1),rain0(nx,1),snow0(nx,1),ice0(nx,1),             &
           graupel0(nx,1),                                              &
           qv1(nx,1,lev),ql1(nx,1,lev),qr1(nx,1,lev),qi1(nx,1,lev),     &
           qs1(nx,1,lev),qg1(nx,1,lev),qa1(nx,1,lev),qn1(nx,1,lev),     &
           pt(nx,1,lev),w(nx,1,lev),uin(nx,1,lev),vin(nx,1,lev),        &
           delp(nx,1,lev),dz(nx,1,lev),                                 &
           qv_dt(nx,1,lev),ql_dt(nx,1,lev),qr_dt(nx,1,lev),             &
           qi_dt(nx,1,lev),qs_dt(nx,1,lev),qg_dt(nx,1,lev),             &
           qa_dt(nx,1,lev),udt(nx,1,lev),vdt(nx,1,lev),pt_dt(nx,1,lev) )
        if ( effr_in ) allocate ( dp(nx,lev),rho(nx,lev) )
        if ( sedi_w ) allocate ( dot(nx,lev) )
        frland = 0.
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

        do i = 1, nxj
          if( islimsk(i) == 1 ) frland(i,1) = 1.  !land fraction
        enddo
         
        do k = 1, lev
          kc = lev - k + 1
          do i = 1, nxj
            prsl(i,k) = 100.0 * plt(i,k)               !layer mean pressure (from mb to Pa)
            if ( sedi_w ) then
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
                  dz, delp, area, dta, frland,                          &
                  rain0, snow0, ice0, graupel0,                         &
                  hydrostatic, phys_hydrostatic,                        &
                  1, nx, 1, 1, 1, lev, 1, lev )

        do k = 1, lev
          do i = 1, nxj
            qtc(i,k)  = qv1(i,1,k) + qv_dt(i,1,k) * dta
            qtr(i,k)  = ql1(i,1,k) + ql_dt(i,1,k) * dta
            qtrw(i,k) = qr1(i,1,k) + qr_dt(i,1,k) * dta
            qti(i,k)  = qi1(i,1,k) + qi_dt(i,1,k) * dta
            qtsw(i,k) = qs1(i,1,k) + qs_dt(i,1,k) * dta
            qtgl(i,k) = qg1(i,1,k) + qg_dt(i,1,k) * dta
            qa(i,k)   = qa1(i,1,k) + qa_dt(i,1,k) * dta
            qt(i,             k) = qtc(i,k)
            qt(i,(ntcw-1)*lev+k) = qtr(i,k)
            qt(i,(ntrw-1)*lev+k) = qtrw(i,k)
            qt(i,(ntiw-1)*lev+k) = qti(i,k)
            qt(i,(ntsw-1)*lev+k) = qtsw(i,k)
            qt(i,(ntgl-1)*lev+k) = qtgl(i,k)
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
               ( 1, nx, 1, lev, rho, dp, islimsk,                       &  ! module_mp_gfdl_fv3_v16.f90
                 qtr, qti, qtrw, qtsw, qtgl, tt,                        &
                 rew, rei, rer, res, reg )
          do k = 1, lev
            kc = lev - k + 1
            do i = 1, nxj
              re_cloud(i,k)   = rew(i,kc)   !(micron)
              re_ice(i,k)     = rei(i,kc)   !(micron)
              re_rain(i,k)    = rer(i,kc)   !(micron)
              re_snow(i,k)    = res(i,kc)   !(micron)
              re_graupel(i,k) = reg(i,kc)
            enddo
          enddo
        endif
!
        do i = 1, nxj
          if ( rain0(i,1)    < rainmin ) rain0(i,1)    = 0.0 
          if ( ice0(i,1)     < rainmin ) ice0(i,1)     = 0.0
          if ( snow0(i,1)    < rainmin ) snow0(i,1)    = 0.0
          if ( graupel0(i,1) < rainmin ) graupel0(i,1) = 0.0

          rlsp(i) = rain0(i,1)+snow0(i,1)+ice0(i,1)+graupel0(i,1)  !total large scale precipitation (mm)
          if ( rlsp(i) > rainmin ) then                         
            sr(i) = (snow0(i,1)+ice0(i,1)+graupel0(i,1))/rlsp(i)   !snow ratio
          else
            sr(i) = 0.0
          endif
        enddo

        deallocate                                                      &
          ( re_graupel,rew,rei,rer,res,reg,                             &
            frland,rain0,snow0,ice0,graupel0,                           &
            qv1,ql1,qr1,qi1,qs1,qg1,qa1,qn1,pt,w,uin,vin,delp,dz,       &
            qv_dt,ql_dt,qr_dt,qi_dt,qs_dt,qg_dt,qa_dt,udt,vdt,pt_dt )
        if ( effr_in ) deallocate ( dp,rho )
        if ( sedi_w ) deallocate ( dot )

      endif  ! end of nmmiph.eq.11

      return
!--------------------------
      end subroutine mp_scheme
!--------------------------
