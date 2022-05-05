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

      return
!--------------------------
      end subroutine mp_init
!--------------------------
!--------------------------
      subroutine mp_scheme                                             &
!--------------------------
!  ---  inputs:
           ( nmmiph,nx,nxj,lev,ncld,plt,pst,dsigma,                    &
             phii,islimsk,q0,kdt,ntcw,ntrw,ntiw,ntsw,ntgl,             &
             ntinc,ntrnc,tpi,me,dta,jj,                                &
!  ---  inputs/outputs:
             tt,qt,                                                    &
!  ---  outputs:
             re_cloud,re_ice,re_snow,rlsp,sr )

! for wsm6
      use module_mp_wsm6,      only: wsm6
! for thompson
      use module_mp_thompson,  only: mp_gt_driver

      implicit none

!  ---  inputs:
      integer,  intent(in)    :: nmmiph,nx,nxj,lev,ncld,kdt,me
      integer,  intent(in)    :: ntcw,ntrw,ntiw,ntsw,ntgl,ntinc,ntrnc
      integer,  intent(in)    :: islimsk(nx)
      real,     intent(in)    :: tpi,dta,jj
      real,     intent(in)    :: plt(nx,lev),pst(nx),dsigma(lev,2),    &
                                 phii(nx,lev+1),q0(nx,lev*ncld)
!  ---  inputs/outputs:
      real,     intent(inout) :: tt(nx,lev),qt(nx,lev*ncld)
!  ---  outputs:
      real,     intent(inout)   :: re_cloud(nx,lev),re_ice(nx,lev),    &
                                   re_snow(nx,lev)
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

      return
!--------------------------
      end subroutine mp_scheme
!--------------------------
