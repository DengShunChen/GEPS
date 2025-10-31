      subroutine gwdps_gpu(im, ix, iy, km, a, b, c, u1, v1, t1, q1, kpbl, &
                           prsi, del, prsl, prslk, phii, phil, deltim, kdt, &
                           hprime, oc, oa4, clx4, theta, sigma, gamma, elvmax, &
                           dusfc, dvsfc, g, cp, rd, rv, imx, &
                           nmtvr, cdmbgwd, me, idxzb, dxmet, hpbl, tofd)
!
!   ********************************************************************
! ----->  i m p l e m e n t a t i o n    v e r s i o n   <----------
!
!          --- not in this code --  history of gwdp at ncep----
!              ----------------     -----------------------
!  version 3  modified for gravity waves, location: .fr30(v3gwd)  *j*
!---       3.1 includes variable saturation flux profile cf isigst
!---       3.g includes ps combined w/ ph (glas and gfdl)
!-----         also included is ri  smooth over a thick lower layer
!-----         also included is decrease in de-acc at top by 1/2
!-----     the nmc gwd incorporating both glas(p&s) and gfdl(migwd)
!-----        mountain induced gravity wave drag
!-----    code from .fr30(v3monnx) for monin3
!-----        this version (06 mar 1987)
!-----        this version (26 apr 1987)    3.g
!-----        this version (01 may 1987)    3.9
!-----    change to fortran 77 (feb 1989)     --- h. juang
!-----    20070601 elvmax bug fix (*j*)
!-----    add dissipation heating (jan 2014)  --- h. juang and f. yang
!
!   version 4
!                ----- this code -----
!
!-----   modified to implement the enhanced low tropospheric gravity
!-----   wave drag developed by kim and arakawa(jas, 1995).
!        orographic std dev (hprime), convexity (oc), asymmetry (oa4)
!        and lx (clx4) are input topographic statistics needed.
!
!-----   programmed and debugged by hong, alpert and kim --- jan 1996.
!-----   debugged again - moorthi and iredell --- may 1998.
!-----
!       further cleanup, optimization and modification
!                                       - s. moorthi may 98, march 99.
!-----   modified for usgs orography data (ncep office note 424)
!        and with several bugs fixed  - moorthi and hong --- july 1999.
!
!-----   modified & implemented into nrl nogaps
!                                       - young-joon kim, july 2000
!-----
!   version lm mb  (6): oz fix 8/2003
!                ----- this code -----
!
!------   changed to include the lott and miller mtn blocking
!         with some modifications by (*j*)  4/02
!        from a principal coordinate calculation using the
!        hi res 8 minute orography, the angle of the
!        mtn with that to the east (x) axis is theta, the slope
!        parameter sigma. the anisotropy is in gamma - all  are input
!        topographic statistics needed.  these are calculated off-line
!        as a function of model resolution in the fortran code ml01rg2.f,
!        with script mlb2.sh.   (*j*)
!-----   gwdps_mb.f version (following lmi) elvmax < hncrit (*j*)
!        mb3a expt to enhance elvmax mtn hgt see sigfac & hncrit
!        gwdps_gwdfix_v6.f fixgwd gf6.0 20070608 sigfac=4.
!-----
!----------------------------------------------------------------------c
!    use
!        routine is called from gbphys  (after call to monnin)
!
!    purpose
!        using the gwd parameterizations of ps-glas and ph-
!        gfdl technique.  the time tendencies of u v
!        are altered to include the effect of mountain induced
!        gravity wave drag from sub-grid scale orography including
!        convective breaking, shear breaking and the presence of
!        critical levels
!
!  input
!        a(iy,km)  non-lin tendency for v wind component
!        b(iy,km)  non-lin tendency for u wind component
!        c(iy,km)  non-lin tendency for temperature
!        u1(ix,km) zonal wind m/sec  at t0-dt
!        v1(ix,km) meridional wind m/sec at t0-dt
!        t1(ix,km) temperature deg k at t0-dt
!        q1(ix,km) specific humidity at t0-dt
!
!        deltim  time step    secs
!        si(n)   p/psfc at base of layer n
!        sl(n)   p/psfc at middle of layer n
!        del(n)  positive increment of p/psfc across layer n
!        kpbl(im) is the index of the top layer of the pbl
!        ipr & lprnt for diagnostics
!
!  output
!        a, b    as augmented by tendency due to gwdps
!                other input variables unmodified.
!   ********************************************************************
!lzl +add === #
         use rank
         use index
         use openacc
         use cudafor
         use param, only: my, my_max
!lzl -end === #
!lzl      use machine , only : kind_phys
         implicit none
         integer im(my), iy, ix, km, imx, kdt, ipr, me

!lzl      integer kpbl(im)                 ! index for the pbl top layer!
         integer kpbl(ix, my_max)                 ! index for the pbl top layer!

!lzl      real(kind=kind_phys) deltim, g, cp, rd, rv,      cdmbgwd(2)
!!lzl     real(kind=kind_phys) a(iy,km),    b(iy,km),      pstar(im)
!lzl      real(kind=kind_phys) a(iy,km),    b(iy,km),       c(iy,km),
!lzl     &                     u1(ix,km),   v1(ix,km),     t1(ix,km),
!lzl     &                     q1(ix,km),   prsi(ix,km+1), del(ix,km),
!lzl     &                     prsl(ix,km), prslk(ix,km),  phil(ix,km),
!lzl     &                     phii(ix,km+1)
!lzl      real(kind=kind_phys) oc(im),     oa4(iy,4), clx4(iy,4)
!lzl     &,                    hprime(im)
!=lzl+add-------------------------------------------------------
         real deltim, g, cp, rd, rv, cdmbgwd(2)
!     real a(iy,km),    b(iy,km),      pstar(im)
         real a(iy, km, my_max), b(iy, km, my_max), c(iy, km, my_max), &
            u1(ix, km, my_max), v1(ix, km, my_max), t1(ix, km, my_max), &
            q1(ix, km, my_max), prsi(ix, km + 1, my_max), del(ix, km, my_max), &
            prsl(ix, km, my_max), prslk(ix, km, my_max), phil(ix, km, my_max), &
            phii(ix, km + 1, my_max)
!lzl      real oc(im),     oa4(iy,4), clx4(iy,4),                            &
!lzl           hprime(im)
         real oc(ix, my_max), oa4(ix, 4, my_max), clx4(ix, 4, my_max), &
            hprime(ix, my_max)

! for lm mtn blocking
!lzl      real(kind=kind_phys) elvmax(im),theta(im),sigma(im),gamma(im)
!lzl      real(kind=kind_phys) wk(im)
!lzl      real(kind=kind_phys) bnv2lm(im,km),pe(im),ek(im),zbk(im),up(im)
!lzl      real(kind=kind_phys) db(im,km),ang(im,km),uds(im,km)
!lzl      real(kind=kind_phys) zlen, dbtmp, r, phiang, cdmb, dbim
!lzl      real(kind=kind_phys) eng0, eng1
!
!=lzl+add---------------------------------------------------------
         real elvmax(ix, my_max), theta(ix, my_max), sigma(ix, my_max), gamma(ix, my_max)
         real wk(ix, my_max)
         real bnv2lm(ix, km, my_max), pe, ek, up
         real db(ix, km, my_max), ang(km), uds(km)
         real zlen, dbtmp, r, phiang, cdmb, dbim
!     real eng0, eng1, eng2
         real eng0, eng1, eng2
!xb118---for TOFD
         logical tofd
         real utendform(ix, km, my_max), vtendform(ix, km, my_max), za, &
            hpbl(ix, my_max), dxmet(ix, my_max)
         real wsp, H_efold, a1, a2, var_temp, dxmeter, ss_taper, ro_tmp, utmp, vtmp
         real varmax_fd, beta_fd, a1_coeff, a2_coeff, TOFD_coeff, Hefold_nom, &
            dxmin_ss, dxmax_ss
         parameter(varmax_fd=160.)
         parameter(beta_fd=0.2)
         parameter(a1_coeff=0.00026615161)  ! Coefficient for TOFD from Beljaars et al. (2004)
         parameter(a2_coeff=0.005363)       !  ""
         parameter(TOFD_coeff=0.0759)       !  ""
         parameter(Hefold_nom=1500.)        ! Nominal TOFD e-folding height (m)
! Small-scale GWD + turbulent form drag
         parameter(dxmin_ss=1000., dxmax_ss=12000.)  ! min,max range of tapering (m)
!xb118---
!     some constants
!
!lzl      real(kind=kind_phys) pi, dw2min, rimin, ric, bnv2min, efmin
!lzl     &,                    efmax,hpmax,hpmin, rad_to_deg, deg_to_rad
!=lzl+add--
         real pi, dw2min, rimin, ric, bnv2min, efmin, &
            efmax, hpmax, hpmin, rad_to_deg, deg_to_rad
!lzl-end
!byl      parameter (pi=3.1415926535897931)
         parameter(pi=4.*atan(1.))
         parameter(rad_to_deg=180.0/pi, deg_to_rad=pi/180.0)
         parameter(dw2min=1., rimin=-100., ric=0.25, bnv2min=1.0e-5)
!     parameter (efmin=0.0, efmax=10.0, hpmax=200.0)
         parameter(efmin=0.0, efmax=10.0, hpmax=2400.0, hpmin=1.0)
!
!lzl      real(kind=kind_phys) frc,    ce,     ceofrc, frmax, cg, gmax
!lzl     &,                    veleps, factop, rlolev, rdi
!lzl      real(kind=kind_phys) critac
!=lzl+add--------------------------------------------------------------
         real frc, ce, ceofrc, frmax, cg, gmax, &
            veleps, factop, rlolev, rdi
         real critac
!lzl-end
         parameter(frc=1.0, ce=0.8, ceofrc=ce/frc, frmax=100., cg=0.5)
         parameter(gmax=1.0, veleps=1.0, factop=0.5)
!     parameter (critac=5.0e-4)
         parameter(rlolev=50000.0)
!     parameter (rlolev=500.0)
!     parameter (rlolev=0.5)
!
!lzl       real(kind=kind_phys) dpmin,hminmt,hncrit,minwnd,sigfac
!
!lzl+add--------------------------------------------------------------
         real dpmin, hminmt, hncrit, minwnd, sigfac

! --- for lm mtn blocking
!     parameter (cdmb = 1.0)     ! non-dim sub grid mtn drag amp (*j*)
         parameter(hncrit=8000.)   ! max value in meters for elvmax (*j*)
!  hncrit set to 8000m and sigfac added to enhance elvmax mtn hgt
         parameter(sigfac=4.0)     ! mb3a expt test for elvmax factor (*j*)
         parameter(hminmt=50.)     ! min mtn height (*j*)
         parameter(minwnd=0.1)     ! min wind component (*j*)

!     parameter (dpmin=00.0)     ! minimum thickness of the reference layer
!!    parameter (dpmin=05.0)     ! minimum thickness of the reference layer
!     parameter (dpmin=20.0)     ! minimum thickness of the reference layer
         ! in centibars
         parameter(dpmin=5000.0)   ! minimum thickness of the reference layer
         ! in pa
!
!lzl      real(kind=kind_phys) fdir
!lzl+add----------------------------------------------------------------
         real fdir
         integer mdir
         parameter(mdir=8, fdir=mdir/(pi + pi))
         integer nwdir(mdir)
         data nwdir/6, 7, 5, 8, 2, 3, 1, 4/
         save nwdir
!
         logical icrilv
!
!----   mountain induced gravity wave drag
!!
!      real(kind=kind_phys) taub(im),  xn(im),     yn(im),    ubar(im)
!     &,                    vbar(im),  ulow(im),   oa(im),    clx(im)
!     &,                    roll(im),  uloi(im),   dusfc(im), dvsfc(im)
!     &,                    dtfac(im), xlinv(im),  delks(im), delks1(im)
!!
!      real(kind=kind_phys) bnv2(im,km),  taup(im,km+1), ri_n(im,km)
!     &,                    taud(im,km),  ro(im,km),     vtk(im,km)
!     &,                    vtj(im,km),   scor(im),      velco(im,km-1)
!     &,                    bnv2bar(im)
!
!lzl+add-----------------------------------------------------------------
!
         real taub, xn(ix, my_max), yn(ix, my_max), ubar(ix, my_max) &
            , vbar(ix, my_max), ulow, oa(ix, my_max), clx(ix, my_max) &
            , roll(ix, my_max), uloi, dusfc(ix, my_max), dvsfc(ix, my_max) &
            , dtfac(ix, my_max), xlinv, delks(ix, my_max), delks1(ix, my_max)
!
         real bnv2(ix, km, my_max), taup(ix, km + 1, my_max), ri_n(ix, km, my_max) &
            , taud(ix, km, my_max), ro(ix, km, my_max) &
            , scor, velco(ix, km - 1, my_max) &
            , bnv2bar(ix, my_max)
!lzl-end-----------------------------------------------------------------
!
!     real(kind=kind_phys) velko(km-1)
         integer kref(ix, my_max), kint(ix, my_max), iwk(ix, my_max), ipt(ix, my_max)
! for lm mtn blocking
         integer kreflm, iwklm(ix, my_max)
         integer idxzb(ix, my_max), ktrial, klevm1, nmtvr
!
!
!      real(kind=kind_phys) gor,    gocp,  fv,    gr2,  bnv,  fr
!     &,                    brvf,   cleff, tem,   tem1,  tem2, temc, temv
!     &,                    wdir,   ti,    rdz,   dw2,   shr2, bvf2
!     &,                    rdelks, efact, coefm, gfobnv
!     &,                    scork,  rscor, hd,    fro,   rim,  sira
!     &,                    dtaux,  dtauy, pkp1log, pklog
!
!lzl+add-----------------------------------------------------------------
         real gor, gocp, fv, gr2, bnv, fr &
            , brvf, cleff, tem, tem1, tem2, temc, temv &
            , wdir, ti, rdz, dw2, shr2, bvf2 &
            , rdelks, efact, coefm, gfobnv &
            , scork, rscor, hd, fro, rim, sira &
            , dtaux, dtauy, pkp1log, pklog &
            , ubartmp, vbartmp, rolltmp, bnv2bartmp, kend
!

         integer kmm1, kmm2, lcap, lcapp1, kbps(my_max), kbpsp1, kbpsm1 &
            , kmps(my_max), idir, nwd, i, j, k, klcap, kp1, kmpbl, npt(my_max), npr &
            , kmll, jj, j1
!    &, kmll,kmds,ihit,jhit
         logical lprnt
         real dutmp, dvtmp, taudtmp, deltmp, vtjp, taup1, &
            vtj0, vtj1, vtk0, vtk1, delp, dtfac_min, t1k0, t1k1, prslk1, &
            atmp, btmp, hprimetmp, velcor, u1r, v1r, kmin

         real time1, time2, time3   !timer
         integer istat, async_id, nptr
         integer(kind=cuda_stream_kind) stream

         async_id = 1
         stream = acc_get_cuda_stream(async_id)

!       parameter (cdmb = 1.0)     ! non-dim sub grid mtn drag amp (*j*)
!   non-dim sub grid mtn drag amp (*j*)
!       cdmb = 1.0/float(imx/192)
!       cdmb = 192.0/float(imx)
         !$acc enter data create(db, ang, uds, taup, npt, ipt, &
         !$acc&      iwklm, ro, bnv2lm, delks, delks1, ubar, vbar, &
         !$acc&      roll, bnv2bar, ri_n, bnv2, iwk, kbps, kmps, kref, &
         !$acc&      oa, clx, xn, yn, dtfac, &
         !$acc&      velco, kint, taud,utendform, vtendform) &
         !$acc& async(async_id)

         nwdir = (/6, 7, 5, 8, 2, 3, 1, 4/)
         !$acc enter data copyin(nwdir) async(async_id)
         cdmb = 4.0*192.0/float(imx)
         if (cdmbgwd(1) >= 0.0) cdmb = cdmb*cdmbgwd(1)
!
         npr = 0
         !jlistnum = 1
         !$acc parallel loop gang private(j1) async(async_id)
         do jj = 1, jlistnum
            j1 = jlist1(jj)
            !$acc loop vector
            do i = 1, im(j1)
               dusfc(i, jj) = 0.
               dvsfc(i, jj) = 0.
            end do
         end do

!
         !$acc host_data use_device(db, ang, uds, taup, utendform, vtendform)
         istat = cudaMemsetAsync(db, 0.0, size(db), stream)
         istat = cudaMemsetAsync(ang, 0.0, size(ang), stream)
         istat = cudaMemsetAsync(uds, 0.0, size(uds), stream)
         istat = cudaMemsetAsync(taup, 0.0, size(taup), stream)
         istat = cudaMemsetAsync(utendform, 0.0, size(utendform), stream)
         istat = cudaMemsetAsync(vtendform, 0.0, size(vtendform), stream)
         !$acc end host_data
         !db=0.
         !ang=0.
         !uds=0.
         !taup=0.
!
         rdi = 1.0/rd
         gor = g/rd
         gr2 = g*gor
         gocp = g/cp
         fv = rv/rd - 1
!
!     ncnt   = 0
         kmm1 = km - 1
         kmm2 = km - 2
         lcap = km
         lcapp1 = lcap + 1

         if ((nmtvr .eq. 14) .or. tofd) then
            !$acc parallel loop private(j1,i,jj,nptr) async(async_id)
            do jj = 1, jlistnum
               j1 = jlist1(jj)
!
!xb118---for TOFD
!---   calculate scale-aware tapering factors

               ipt(:, jj) = 0
               nptr = 0
! ----    for lm and gwd calculation points
               !$acc loop seq
               do i = 1, im(j1)
                  if ((elvmax(i, jj) .gt. hminmt) &
                      .and. (hprime(i, jj) .gt. hpmin)) then
                     nptr = nptr + 1
                     ipt(nptr, jj) = i
!               if (ipr .eq. i) npr = npt
                  end if
               end do
               npt(jj) = nptr
               !if (npt .eq. 0) return     ! no gwd/mb calculation done!
            end do
         end if

         if (tofd) then
            !$acc parallel loop collapse(2) gang private(i,j,k,var_temp,a1,a2,H_efold,za,vtjp,&
            !$acc&                           utmp,vtmp,wsp,ro_tmp,hprimetmp,&
            !$acc&                           atmp,btmp,dxmeter) async(async_id)
            do jj = 1, jlistnum
               do k = 1, km
                  !$acc loop vector
                  do i = 1, npt(jj)
                     j = ipt(i, jj)
                     dxmeter = sqrt(dxmet(j, jj))
                     if (dxmeter .ge. dxmax_ss) then
                        ss_taper = 1.
                     else
                        if (dxmeter .le. dxmin_ss) then
                           ss_taper = 0.
                        else
                           ss_taper = dxmax_ss*(1.-dxmin_ss/dxmeter)/(dxmax_ss - dxmin_ss)
                        end if
                     end if

                     if (ss_taper .gt. 1.e-2) then
                        atmp = a(j, k, jj)
                        btmp = b(j, k, jj)
                        hprimetmp = hprime(j, jj)
                        var_temp = MIN(hprimetmp, varmax_fd) + &
                                   MAX(0., beta_fd*(hprimetmp - varmax_fd))
                        a1 = a1_coeff*var_temp**2
                        a2 = a1*a2_coeff
                        ! Revise e-folding height based on PBL height and topographic std. dev. -- M. Toy 3/12/2018
                        H_efold = max(2*hprimetmp, hpbl(j, jj))
                        H_efold = min(H_efold, Hefold_nom)
                        za = 0.5*(phii(j, k, jj) + phii(j, k + 1, jj))/g
                        !              wsp=SQRT(u1(j,k)**2 + v1(j,k)**2)
                        utmp = u1(j, k, jj) + btmp*deltim
                        vtmp = v1(j, k, jj) + atmp*deltim
                        wsp = SQRT(utmp*utmp + vtmp*vtmp)
                        vtjp = t1(j, k, jj)*(1.+fv*q1(j, k, jj))
                        ro_tmp = rdi*prsl(j, k, jj)/vtjp ! density tons/m**3
                        ! Eqn. (16) of Beljaars et al. (2004)
                        utendform(j, k, jj) = -TOFD_coeff*wsp*utmp*ro_tmp* &
                                              EXP(-(za/H_efold)**1.5)*a2*za**(-1.2)*ss_taper
                        vtendform(j, k, jj) = -TOFD_coeff*wsp*vtmp*ro_tmp* &
                                              EXP(-(za/H_efold)**1.5)*a2*za**(-1.2)*ss_taper
                        a(j, k, jj) = vtendform(j, k, jj) + atmp
                        b(j, k, jj) = utendform(j, k, jj) + btmp
                     end if ! ss_taper
                  end do  ! enddo k
               end do ! enddo i
            end do
         end if ! tofd
!!xb118---
!
!
!         if (lprnt) print *,' npt=',npt,' npr=',npr,' ipr=',ipr,' im=',im
!
         !
!         if (lprnt)
!      &  print *,' in gwdps_lm.f npt,im,ix,iy,km,me=',npt,im,ix,iy,km,me
!
!
!   start lm mtn blocking (mb) section
!
!..............................
!..............................
!
!    (*j*)  11/03:  test upper limit on kmll=km - 1
!        then do not need hncrit -- test with large hncrit first.
!         kmll  = km / 2 ! maximum mtnlm height : # of vertical levels / 2
         kmll = kmm1
! ---   no mtn should be as high as kmll (so we do not have to start at
! ---   the top of the model but could do calc for all levels).
!
!
         if (nmtvr .eq. 14) then
            !$acc parallel loop collapse(2) gang vector private(i,j,k,rdz,pkp1log,pklog,vtjp) async(async_id)
            do jj = 1, jlistnum
               do i = 1, ix
                  if (i .le. npt(jj)) then
                     j = ipt(i, jj)
                     ! ---   iwklm is the level above the height of the of the mountain.
                     ! ---   idxzb is the level of the dividing streamline.
                     !   initialize dividing streamline (ds) control vector
                     iwklm(i, jj) = 2
                     idxzb(i, jj) = 0
                     elvmax(j, jj) = min(elvmax(j, jj) + sigfac*hprime(j, jj), hncrit)
                     !$acc loop seq
                     do k = 1, kmll
                        ! ---   interpolate to max mtn height for index, iwklm(i) wk[gz]
                        ! ---   elvmax is limited to hncrit because to hi res topo30 orog.
                        pkp1log = phil(j, k + 1, jj)/g
                        pklog = phil(j, k, jj)/g
                        !              if(myrank.eq.0)print*,'lzl gwdps pkp1log pklog= ',pkp1log,pklog
  !!!-------       elvmax(j) = min (elvmax(j) + sigfac * hprime(j), hncrit)
                        if ((elvmax(j, jj) .le. pkp1log) .and. &
                            (elvmax(j, jj) .ge. pklog)) then
                           !       print *,' in gwdps_lm.f 1  =',k,elvmax(j),pklog,pkp1log,me
                           ! ---          wk for diags but can be saved and reused.
                           iwklm(i, jj) = max(iwklm(i, jj), k + 1)
                           !       print *,' in gwdps_lm.f 2 npt=',npt,i,j,wk(i),iwklm(i),me
                        end if
                        !
                        ! ---          find at prsl levels large scale environment variables
                        ! ---          these cover all possible mtn max heights
                        vtj0 = t1(j, k, jj)*(1.+fv*q1(j, k, jj))
                        ro(i, k, jj) = rdi*prsl(j, k, jj)/vtj0 ! density kg/m**3
                        !
                        !   testing for highest model level of mountain top
                        !
                        !           ihit = 2
                        !           jhit = 0
                        !          do i = 1, npt
                        !          j=ipt(i)
                        !            if ( iwklm(i) .gt. ihit ) then
                        !              ihit = iwklm(i)
                        !              jhit = j
                        !            endif
                        !          enddo
                        !       print *, ' mb: kdt,max(iwklm),jhit,phil,me=',
                        !      &          kdt,ihit,jhit,phil(jhit,ihit),me
                        if (k .ne. kmll) then
                           vtk0 = vtj0/prslk(j, k, jj)
                           vtj1 = t1(j, k + 1, jj)*(1.+fv*q1(j, k + 1, jj))
                           vtk1 = vtj1/prslk(j, k + 1, jj)

                           rdz = g/(phil(j, k + 1, jj) - phil(j, k, jj))
                           ! ---                                 brunt-vaisala frequency
                           bnv2lm(i, k, jj) = max((g + g)*rdz*(vtk1 - vtk0)/(vtk1 + vtk0), &
                                                  bnv2min)
                        end if
                     end do
                  end if
               end do
            end do
         end if
!      print *,' in gwdps_lm.f 3 npt=',npt,j,rdz,me
!

!       print *,' in gwdps_lm.f 4 npt=',npt,kreflm(npt),me
!
! ---   in the layer kreflm(i) to 1 find pe (which needs n, elvmax)
! ---    make averages, guess dividing stream (ds) line layer.
! ---    this is not used in the first cut except for testing and
! ---   is the vert ave of quantities from the surface to mtn top.
!

         if (nmtvr .eq. 14) then
            !$acc parallel loop collapse(2) private(i,j,k,rdelks,ubartmp,vbartmp,&
            !$acc&                          rolltmp,bnv2bartmp,kreflm,phiang,&
            !$acc&                          zlen,r,dbtmp,ang,uds,pe) async(async_id)
            do jj = 1, jlistnum
               do i = 1, ix
                  if (i .le. npt(jj)) then
                     j = ipt(i, jj)
                     delks(i, jj) = 1.0/(prsi(j, 1, jj) - prsi(j, iwklm(i, jj), jj))
                     delks1(i, jj) = 1.0/(prsl(j, 1, jj) - prsl(j, iwklm(i, jj), jj))
                     ubar(i, jj) = 0.0
                     vbar(i, jj) = 0.0
                     roll(i, jj) = 0.0
                     pe = 0.0
                     ek = 0.0
                     bnv2bar(i, jj) = (prsl(j, 1, jj) - prsl(j, 2, jj))*delks1(i, jj)*bnv2lm(i, 1, jj)

                     ubartmp = 0.
                     vbartmp = 0.
                     rolltmp = 0.
                     bnv2bartmp = 0.
! ---   find the dividing stream line height
! ---   starting from the level above the max mtn downward
! ---   iwklm(i) is the k-index of mtn elvmax elevation
                     kreflm = 0
                     !$acc loop seq
                     do ktrial = kmll, 1, -1
                        if (ktrial .lt. iwklm(i, jj) .and. kreflm .eq. 0) then
                           kreflm = ktrial
                        end if
                     end do
                     !$acc loop seq
                     do k = 1, kreflm
                        rdelks = del(j, k, jj)*delks(i, jj)
                        ubartmp = ubartmp + rdelks*u1(j, k, jj) ! trial mean u below
                        vbartmp = vbartmp + rdelks*v1(j, k, jj) ! trial mean v below
                        rolltmp = rolltmp + rdelks*ro(i, k, jj) ! trial mean ro below
                        rdelks = (prsl(j, k, jj) - prsl(j, k + 1, jj))*delks1(i, jj)
                        bnv2bartmp = bnv2bartmp + bnv2lm(i, k, jj)*rdelks
                        ! ---   these vert ave are for diags, testing and gwd to follow (*j*).
                     end do
                     ubar(i, jj) = ubartmp ! trial mean u below
                     vbar(i, jj) = vbartmp ! trial mean v below
                     roll(i, jj) = rolltmp ! trial mean ro below
                     bnv2bar(i, jj) = bnv2bar(i, jj) + bnv2bartmp
                     !       print *,' in gwdps_lm.f 5  =',i,kreflm(npt),bnv2bar(npt),me

!---------------------------------------
! ---   integrate to get pe in the trial layer.
! ---   need the first layer where pe>ek - as soon as
! ---   idxzb is not 0 we have a hit and zb is found.
!
                     !$acc loop seq
                     do k = iwklm(i, jj), 1, -1
                        phiang = atan2(v1(j, k, jj), u1(j, k, jj))*rad_to_deg
                        ang(k) = (theta(j, jj) - phiang)
                        if (ang(k) .gt. 90.) ang(k) = ang(k) - 180.
                        if (ang(k) .lt. -90.) ang(k) = ang(k) + 180.
                        ang(k) = ang(k)*deg_to_rad
                        !
                        uds(k) = &
                           max(sqrt(u1(j, k, jj)*u1(j, k, jj) + v1(j, k, jj)*v1(j, k, jj)), minwnd)
                        ! ---   test to see if we found zb previously
                        if (idxzb(i, jj) .eq. 0) then
                           pe = pe + bnv2lm(i, k, jj)* &
                                (g*elvmax(j, jj) - phil(j, k, jj))* &
                                (phii(j, k + 1, jj) - phii(j, k, jj))/(g*g)
                           ! ---   ke
                           ! ---   wind projected on the line perpendicular to mtn range, u(zb(k)).
                           ! ---   kenetic energy is at the layer zb
                           ! ---   theta ranges from -+90deg |_ to the mtn "largest topo variations"
                           up = uds(k)*cos(ang(k))
                           ek = 0.5*up*up

                           ! ---   dividing stream line  is found when pe =exceeds ek.
                           if (pe .ge. ek) then
                              idxzb(i, jj) = k
                           end if
                           ! ---   then mtn blocked flow is between zb=k(idxzb(i)) and surface
                           !
                        end if
                     end do
                     !zbk(i,jj) = elvmax(j,jj)                                             &
                     !    - sqrt(ubar(i,jj)*ubar(i,jj) + vbar(i,jj)*vbar(i,jj))/bnv2bar(i,jj)

!-----------------------------------
!
!       print *,' in gwdps_lm.f 6  =',phiang,theta(ipt(npt)),me
!       print *,' in gwdps_lm.f 7  =',idxzb(npt),pe(npt)
!
!       if (lprnt .and. npr .gt. 0) then
!         print *,' bnv2bar,bnv2lm=',bnv2bar(npr),bnv2lm(npr,1:klevm1)
!         print *,' npr,idxzb,uds=',npr,idxzb(npr),uds(npr,:)
!         print *,' pe,up,ek=',pe(npr),up(npr),ek(npr)
!       endif
!
!
!       if (lprnt .and. npr .gt. 0) then
!         print *,' iwklm,zbk=',iwklm(npr),zbk(npr),idxzb(npr)
!         print *,' zb=',phil(ipr),idxzb(npr))/g
!       print *,' in gwdps_lm.f 8 npt =',npt,zbk(npt),up(npt),me
!       endif
!
! ---   the drag for mtn blocked flow
!
                     !        print *,' in gwdps_lm.f 9  =',i,j,idxzb(i),me
                     if (idxzb(i, jj) .gt. 0) then
                        !$acc loop seq
                        do k = idxzb(i, jj), 1, -1
                           if (phil(j, idxzb(i, jj), jj) .gt. phil(j, k, jj)) then
                              zlen = sqrt((phil(j, idxzb(i, jj), jj) - phil(j, k, jj))/ &
                                          (phil(j, k, jj) + g*hprime(j, jj)))
                              ! ---   lm eq 14:
                              r = (cos(ang(k))**2 + gamma(j, jj)*sin(ang(k))**2)/ &
                                  (gamma(j, jj)*cos(ang(k))**2 + sin(ang(k))**2)
                              ! ---   (negitive of db -- see sign at tendency)
                              dbtmp = 0.25*cdmb* &
                                      max(2.-1./r, 0.)*sigma(j, jj)* &
                                      max(cos(ang(k)), gamma(j, jj)*sin(ang(k)))* &
                                      zlen/hprime(j, jj)
                              db(i, k, jj) = dbtmp*uds(k)*ro(i, k, jj)
                              !
                              !                 if(lprnt .and. i .eq. npr) then
                              !                   print *,' in gwdps_lmi.f 10 npt=',npt,i,j,idxzb(i)
                              !      &,           dbtmp,r' ang=',ang(i,k),' gamma=',gamma(j),' k=',k
                              !                   print *,' in gwdps_lmi.f 11   k=',k,zlen,cos(ang(i,k))
                              !                   print *,' in gwdps_lmi.f 12  db=',db(i,k),sin(ang(i,k))
                              !                 endif
                           end if
                        end do
                        !           if(lprnt) print *,' @k=1,zlen,dbtmp=',k,zlen,dbtmp
                     end if
                  end if
               end do
            end do
         end if

!
!.............................
!.............................
! end  mtn blocking section
!
         if (nmtvr .ne. 14) then
            !$acc parallel loop private(j1,i,jj,nptr) async(async_id)
            do jj = 1, jlistnum
               j1 = jlist1(jj)
! ----    for mb not present and  gwd (nmtvr .ne .14)
               ipt(:, jj) = 0
               nptr = 0
               !$acc loop seq
               do i = 1, im(j1)
                  if (hprime(i, jj) .gt. hpmin) then
                     nptr = nptr + 1
                     ipt(npt(jj), jj) = i
!               if (ipr .eq. i) npr = npt
                  end if
               end do
               npt(jj) = nptr
               !if (npt .eq. 0) return     ! no gwd/mb calculation done!
            end do
         end if
!
!          if(myrank.eq.0)print *,' in gwdps npr=',npr,' npt=',npt,' ipr=',ipr
!         if (lprnt) print *,' npr=',npr,' npt=',npt,' ipr=',ipr
!        &,' ipt(npt)=',ipt(npt)
!

!
!.............................
!.............................
!
         kmpbl = km/2 ! maximum pbl height : # of vertical levels / 2
!
!    scale cleff between im=384*2 and 192*2 for t126/t170 and t62
!

         if (imx .gt. 0) then
!       cleff = 1.0e-5 * sqrt(float(imx)/384.0) !  this is inverse of cleff!
!       cleff = 1.0e-5 * sqrt(float(imx)/192.0) !  this is inverse of cleff!
!       cleff = 0.5e-5 * sqrt(float(imx)/192.0) !  this is inverse of cleff!
!       cleff = 1.0e-5 * sqrt(float(imx)/192)/float(imx/192)
!       cleff = 1.0e-5 / sqrt(float(imx)/192.0) !  this is inverse of cleff!
            cleff = 0.5e-5/sqrt(float(imx)/192.0) !  this is inverse of cleff!
!       cleff = 2.0e-5 * sqrt(float(imx)/192.0) !  this is inverse of cleff!
!       cleff = 2.5e-5 * sqrt(float(imx)/192.0) !  this is inverse of cleff!
         end if
         if (cdmbgwd(2) >= 0.0) cleff = cleff*cdmbgwd(2)
!

         !$acc parallel loop collapse(3) private(k,i,j,ti,tem,rdz,tem1,tem2,&
         !$acc&                    dw2,shr2,bvf2,vtj0,vtj1,vtk0,vtk1, &
         !$acc&                    t1k0,t1k1,prslk1) async(async_id)
         do jj = 1, jlistnum
            do k = 1, kmm1
               do i = 1, ix
                  if (i .le. npt(jj)) then
                     j = ipt(i, jj)
                     t1k0 = t1(j, k, jj)
                     t1k1 = t1(j, k + 1, jj)
                     prslk1 = prsl(j, k + 1, jj)

                     vtj0 = t1k0*(1.+fv*q1(j, k, jj))
                     vtk0 = vtj0/prslk(j, k, jj)
                     vtj1 = t1k1*(1.+fv*q1(j, k + 1, jj))
                     vtk1 = vtj1/prslk(j, k + 1, jj)
                     if (k .eq. 1) ro(i, k, jj) = rdi*prsl(j, k, jj)/vtj0 ! density tons/m**3
                     ro(i, k + 1, jj) = rdi*prslk1/vtj1 ! density tons/m**3
                     ti = 2.0/(t1k0 + t1k1)
                     tem = ti/(prsl(j, k, jj) - prslk1)
                     rdz = g/(phil(j, k + 1, jj) - phil(j, k, jj))
                     tem1 = u1(j, k, jj) - u1(j, k + 1, jj)
                     tem2 = v1(j, k, jj) - v1(j, k + 1, jj)
                     dw2 = tem1*tem1 + tem2*tem2
                     shr2 = max(dw2, dw2min)*rdz*rdz
                     bvf2 = g*(gocp + rdz*(vtj1 - vtj0))*ti
                     ri_n(i, k, jj) = max(bvf2/shr2, rimin)   ! richardson number
                     !                                                brunt-vaisala frequency
                     !           tem       = gr2 * (prsl(j,k)+prsl(j,k+1)) * tem
                     !           bnv2(i,k) = tem * (vtk(i,k+1)-vtk(i,k))/(vtk(i,k+1)+vtk(i,k))
                     bnv2(i, k, jj) = max((g + g)*rdz*(vtk1 - vtk0)/(vtk1 + vtk0), &
                                          bnv2min)
                  end if
               end do
            end do
         end do

!        if(myrank.eq.0)print *,' in gwdps npt,kmm1,bnv2=',npt,kmm1,bnv2(npt,kmm1)
!        print *,' in gwdps_lm.f gwd:14  =',npt,kmm1,bnv2(npt,kmm1)
!
!       apply 3 point smoothing on bnv2
!
!       do k=1,km
!         do i=1,im
!           vtk(i,k) = bnv2(i,k)
!         enddo
!       enddo
!       do k=2,kmm1
!         do i=1,im
!           bnv2(i,k) = 0.25*(vtk(i,k-1)+vtk(i,k+1)) + 0.5*vtk(i,k)
!         enddo
!       enddo
!
!       finding the first interface index above 50 hpa level
!

         !$acc parallel loop collapse(2) private(i,k,j,tem) async(async_id)
         do jj = 1, jlistnum
            do i = 1, ix
               if (i .le. npt(jj)) then
                  if (nmtvr .ne. 14) idxzb(i, jj) = 0
                  iwk(i, jj) = 2
                  !$acc loop seq
                  do k = 3, kmpbl
                     j = ipt(i, jj)
                     tem = (prsi(j, 1, jj) - prsi(j, k, jj))
                     !            if(myrank.eq.0)print*,'lzl gwdps tem= ',tem
                     if (tem .lt. dpmin) iwk(i, jj) = k
                  end do
               end if
            end do
         end do

!

         !$acc parallel loop gang private(i,j) async(async_id)
         do jj = 1, jlistnum
            kbps(jj) = 1
            kmps(jj) = km
            !$acc loop vector
            do i = 1, npt(jj)
               j = ipt(i, jj)
               kref(i, jj) = max(iwk(i, jj), kpbl(j, jj) + 1) ! reference level
               delks(i, jj) = 1.0/(prsi(j, 1, jj) - prsi(j, kref(i, jj), jj))
               delks1(i, jj) = 1.0/(prsl(j, 1, jj) - prsl(j, kref(i, jj), jj))
               ubar(i, jj) = 0.0
               vbar(i, jj) = 0.0
               roll(i, jj) = 0.0
               bnv2bar(i, jj) = (prsl(j, 1, jj) - prsl(j, 2, jj))*delks1(i, jj)*bnv2(i, 1, jj)
            end do
            !$acc loop seq
            do i = 1, npt(jj)
               kbps(jj) = max(kbps(jj), kref(i, jj))
               kmps(jj) = min(kmps(jj), kref(i, jj))
            end do
         end do

         !$acc parallel loop collapse(2) gang private(i,j,k,rdelks,ubartmp,vbartmp, &
         !$acc&                          rolltmp,bnv2bartmp,velcor,bnv,fr,efact,&
         !$acc&                          coefm,tem,gfobnv,kmin,kp1,temv,scork,&
         !$acc&                          rscor,brvf,tem1,hd,fro,tem2,rim,temc,&
         !$acc&                          taup1,klcap,sira,taub,ulow,uloi,xlinv,scor) async(async_id)
         do jj = 1, jlistnum
!        if(myrank.eq.0)print *,' in gwdps kbps,kmps=',kbps,kmps
!        print *,' in gwdps_lm.f gwd:15  =',kbps,kmps
            do i = 1, ix
               if (i .le. npt(jj)) then
                  j = ipt(i, jj)
                  ubartmp = 0.
                  vbartmp = 0.
                  rolltmp = 0.
                  bnv2bartmp = 0.
                  kend = min(kbps(jj), kref(i, jj) - 1)
            !!!$acc loop vector reduction(+:ubartmp) &
            !!!$acc&            reduction(+:vbartmp) &
            !!!$acc&            reduction(+:rolltmp) &
            !!!$acc&            reduction(+:bnv2bartmp)
                  !$acc loop seq
                  do k = 1, kend
                     rdelks = del(j, k, jj)*delks(i, jj)
                     ubartmp = ubartmp + rdelks*u1(j, k, jj)   ! mean u below kref
                     vbartmp = vbartmp + rdelks*v1(j, k, jj)   ! mean v below kref
                     rolltmp = rolltmp + rdelks*ro(i, k, jj)   ! mean ro below kref
                     rdelks = (prsl(j, k, jj) - prsl(j, k + 1, jj))*delks1(i, jj)
                     bnv2bartmp = bnv2bartmp + bnv2(i, k, jj)*rdelks
                  end do
                  ubar(i, jj) = ubar(i, jj) + ubartmp
                  vbar(i, jj) = vbar(i, jj) + vbartmp
                  roll(i, jj) = roll(i, jj) + rolltmp
                  bnv2bar(i, jj) = bnv2bar(i, jj) + bnv2bartmp

!       figure out low-level horizontal wind direction and find 'oa'
!
!               nwd  1   2   3   4   5   6   7   8
!                wd  w   s  sw  nw   e   n  ne  se
                  wdir = atan2(ubar(i, jj), vbar(i, jj)) + pi
                  idir = mod(nint(fdir*wdir), mdir) + 1
                  nwd = nwdir(idir)
                  oa(i, jj) = (1 - 2*int((nwd - 1)/4))*oa4(j, mod(nwd - 1, 4) + 1, jj)
                  clx(i, jj) = clx4(j, mod(nwd - 1, 4) + 1, jj)

!
!-----  xn,yn            "low-level" wind projections in zonal
!                                      & meridional directions
!-----  ulow             "low-level" wind magnitude -        (= u)
!-----  bnv2             bnv2 = n**2
!-----  taub             base momentum flux
!-----  = -(ro * u**3/(n*xl)*gf(fr) for n**2 > 0
!-----  = 0.                        for n**2 < 0
!-----  fr               froude    =   n*hprime / u
!-----  g                gmax*fr**2/(fr**2+cg/oc)
!
!-----  initialize some arrays
!
                  xn(i, jj) = 0.0
                  yn(i, jj) = 0.0
                  taub = 0.0
                  ulow = 0.0
                  dtfac(i, jj) = 1.0
!
!----  compute the "low level" wind magnitude (m/s)
!
                  ulow = max(sqrt(ubar(i, jj)*ubar(i, jj) + vbar(i, jj)*vbar(i, jj)), 1.0)
                  uloi = 1.0/ulow

!----------------------------------------------
                  !$acc loop seq
                  do k = 1, kmm1
                     velcor = 0.5*((u1(j, k, jj) + u1(j, k + 1, jj))*ubar(i, jj) &
                                   + (v1(j, k, jj) + v1(j, k + 1, jj))*vbar(i, jj))
                     velco(i, k, jj) = velcor*uloi
!           if ((velco(i,k).lt.veleps) .and. (velco(i,k).gt.0.)) then
!             velco(i,k) = veleps
!           endif
                  end do
!
!       if(lprnt) print *,' ubar=',ubar
!      &,' vbar=',vbar,' ulow=',ulow,' veleps=',veleps
!
!----------------------------------------------
!     find the interface level of the projected wind where
!     low levels & upper levels meet above pbl
!
!     note following not being used (kint reset to kref) so commented out
!       do i=1,npt
!         kint(i) = km
!       enddo
!       do k = 1,kmm1
!         do i = 1,npt
!           if (k .gt. kref(i)) then
!             if(velco(i,k) .lt. veleps .and. kint(i) .eq. km) then
!               kint(i) = k+1
!             endif
!           endif
!         enddo
!       enddo
!    warning  kint = kref !!!!!!!!!
                  kint(i, jj) = kref(i, jj)
                  bnv = sqrt(bnv2bar(i, jj))
                  fr = bnv*uloi*min(hprime(j, jj), hpmax)
                  fr = min(fr, frmax)
                  xn(i, jj) = ubar(i, jj)*uloi
                  yn(i, jj) = vbar(i, jj)*uloi
                  !
                  !       compute the base level stress and store it in taub
                  !       calculate enhancement factor, number of mountains & aspect
                  !       ratio const. use simplified relationship between standard
                  !       deviation & critical hgt
                  !
                  efact = (oa(i, jj) + 2.)**(ceofrc*fr)
                  efact = min(max(efact, efmin), efmax)
                  !
                  coefm = (1.+clx(i, jj))**(oa(i, jj) + 1.)
                  !
                  xlinv = coefm*cleff
                  !
                  tem = fr*fr*oc(j, jj)
                  gfobnv = gmax*tem/((tem + cg)*bnv)  ! g/n0
                  !
                  taub = xlinv*roll(i, jj)*ulow*ulow &
                         *ulow*gfobnv*efact         ! base flux tau0
                  !
                  !           tem      = min(hprime(i),hpmax)
                  !           taub(i)  = xlinv(i) * roll(i) * ulow(i) * bnv * tem * tem
                  !
                  k = max(1, kref(i, jj) - 1)
                  tem = max(velco(i, k, jj)*velco(i, k, jj), 0.1)
                  scor = bnv2(i, k, jj)/tem  ! scorer parameter below ref level

!------------------------------------------------
!       if(lprnt) print *,' taub=',taub
!
!----  set up bottom values of stress
!

                  kmin = min(kbps(jj), kref(i, jj))
                  !$acc loop seq
                  do k = 1, kmin
                     taup(i, k, jj) = taub
                  end do

!-------------------------------------------------
!     now compute vertical structure of the stress.
!

                  !$acc loop seq
                  do k = kref(i, jj), kmm1                   ! vertical level k loop!
                     kp1 = k + 1
                     icrilv = .false. ! initialize critical level control vector

                     !
                     !-----  unstable layer if ri < ric
                     !-----  unstable layer if upper air vel comp along surf vel <=0 (crit lay)
                     !----   at (u-c)=0. crit layer exists and bit vector should be set (.le.)
                     !
                     icrilv = icrilv .or. (ri_n(i, k, jj) .lt. ric) &
                              .or. (velco(i, k, jj) .le. 0.0)

                     if (.not. icrilv .and. taup(i, k, jj) .gt. 0.0) then
                        temv = 1.0/max(velco(i, k, jj), 0.01)
                        !               if (oa(i) .gt. 0. .and.  prsi(ipt(i),kp1).gt.rlolev) then
                        if (oa(i, jj) .gt. 0. .and. kp1 .lt. kint(i, jj)) then
                           scork = bnv2(i, k, jj)*temv*temv
                           rscor = min(1.0, scork/scor)
                           scor = scork
                        else
                           rscor = 1.
                        end if
                        !
                        brvf = sqrt(bnv2(i, k, jj))        ! brunt-vaisala frequency
                        !               tem1 = xlinv(i)*(ro(i,kp1)+ro(i,k))*brvf*velco(i,k)*0.5
                        tem1 = xlinv*(ro(i, kp1, jj) + ro(i, k, jj))*brvf*0.5 &
                               *max(velco(i, k, jj), 0.01)
                        hd = sqrt(taup(i, k, jj)/tem1)
                        fro = brvf*hd*temv
                        !
                        !      rim is the  minimum-richardson number by shutts (1985)
                        !
                        tem2 = sqrt(ri_n(i, k, jj))
                        tem = 1.+tem2*fro
                        rim = ri_n(i, k, jj)*(1.-fro)/(tem*tem)
                        !
                        !      check stability to employ the 'saturation hypothesis'
                        !      of lindzen (1981) except at tropospheric downstream regions
                        !
                        !                                         ----------------------
                        if (rim .le. ric .and. &
                            !      &           (oa(i) .le. 0. .or.  prsi(ipt(i),kp1).le.rlolev )) then
                            (oa(i, jj) .le. 0. .or. kp1 .ge. kint(i, jj))) then
                           temc = 2.0 + 1.0/tem2
                           hd = velco(i, k, jj)*(2.*sqrt(temc) - temc)/brvf
                           taup1 = tem1*hd*hd
                        else
                           taup1 = taup(i, k, jj)*rscor
                        end if
                        taup(i, kp1, jj) = min(taup1, taup(i, k, jj))
                     end if
                  end do

!
!       do i=1,im
!         taup(i,km+1) = taup(i,km)
!       enddo
!
!--------------------------------------
                  if (lcap .le. km) then
                     !$acc loop seq
                     do klcap = lcapp1, km + 1
                        sira = prsi(ipt(i, jj), klcap, jj)/prsi(ipt(i, jj), lcap, jj)
                        taup(i, klcap, jj) = sira*taup(i, lcap, jj)
                     end do
                  end if
               end if
            end do
         end do

!
!       calculate - (g/p*)*d(tau)/d(sigma) and decel terms dtaux, dtauy
!

         !$acc parallel loop collapse(2) gang vector private(i,k,klcap,dtfac_min) async(async_id)
         do jj = 1, jlistnum
            do i = 1, ix
               if (i .le. npt(jj)) then
                  !$acc loop seq
                  do k = 1, km
                     taud(i, k, jj) = g*(taup(i, k + 1, jj) - taup(i, k, jj))/del(ipt(i, jj), k, jj)
                     if (k .ge. lcap) taud(i, k, jj) = taud(i, k, jj)*factop
                  end do
!---------------------------------------
!------  limit de-acceleration (momentum deposition ) at top to 1/2 value
!------  the idea is some stuff must go out the 'top'
!

!
!------  if the gravity wave drag would force a critical line in the
!------  layers below sigma=rlolev during the next deltim timestep,
!------  then only apply drag until that critical line is reached.
!
                  dtfac_min = dtfac(i, jj)
            !!!$acc loop vector reduction(min:dtfac_min)
                  !$acc loop seq
                  do k = kref(i, jj) + 1, kmm1
                     if (prsi(ipt(i, jj), k, jj) .ge. rlolev .and. taud(i, k, jj) .ne. 0.) then
                        tem = deltim*taud(i, k, jj)
                        dtfac_min = min(dtfac_min, abs(velco(i, k, jj)/tem))
                     end if
                  end do
                  dtfac(i, jj) = dtfac_min
               end if
            end do
         end do

!
!       if(lprnt .and. npr .gt. 0) then
!         if(myrank.eq.0)print *,' in gwdps before  a=',a(npr,:)
!         if(myrank.eq.0)print *,' in gwdps before  b=',b(npr,:)
!         print *,' before  a=',a(npr,:)
!         print *,' before  b=',b(npr,:)
!       endif
         !$acc parallel loop collapse(3) gang vector private(i,j,k,dtaux,dtauy,&
         !$acc&                          eng0,dbim,eng1,eng2,taudtmp,u1r,v1r) async(async_id)
         do jj = 1, jlistnum
            do k = 1, km            !org
               do i = 1, ix            !org
                  if (i .le. npt(jj)) then
                     j = ipt(i, jj)   !org
                     taudtmp = taud(i, k, jj)*dtfac(i, jj)
                     dtaux = taudtmp*xn(i, jj)
                     dtauy = taudtmp*yn(i, jj)
                     u1r = u1(j, k, jj)
                     v1r = v1(j, k, jj)
                     eng0 = 0.5*(u1r**2.0 + v1r**2.0)
                     ! ---    lm mb (*j*)  changes overwrite gwd
                     if (k .lt. idxzb(i, jj) .and. idxzb(i, jj) .ne. 0) then
                        dbim = db(i, k, jj)/(1.+db(i, k, jj)*deltim)
                        a(j, k, jj) = -dbim*v1r + a(j, k, jj)
                        b(j, k, jj) = -dbim*u1r + b(j, k, jj)
                        eng1 = eng0*(1.0 - dbim*deltim)**2.0
                        !            if ( abs(dbim * u1(j,k)) .gt. .01 )
                        !      & print *,' in gwdps_lmi.f kdt=',kdt,i,k,db(i,k),
                        !      &                      dbim,idxzb(i),u1(j,k),v1(j,k),me
                     else
                        !
                        a(j, k, jj) = dtauy + a(j, k, jj)
                        b(j, k, jj) = dtaux + b(j, k, jj)
                        eng1 = 0.5*((u1r + dtaux*deltim)**2.0 + &
                                    (v1r + dtauy*deltim)**2.0)
                     end if

                     eng2 = 0.5*((u1r + utendform(j, k, jj)*deltim)**2.0 + &
                                 (v1r + vtendform(j, k, jj)*deltim)**2.0)
                     !         c(j,k) = c(j,k) + max((eng0-eng1-eng2),0.0)/cp/deltim
                     !org          c(j,k) = c(j,k) + max((eng0-eng1),0.0)/cp/deltim   #need to make sure the meaning of the constraint of maximum value
                     c(j, k, jj) = c(j, k, jj) + (eng0 - eng1)/(cp*deltim)
                     c(j, k, jj) = c(j, k, jj) + (eng0 - eng2)/(cp*deltim)
                     !
                     !          u1(j,k) = u1(j,k) + b(j,k) * deltim
                     !          v1(j,k) = v1(j,k) + a(j,k) * deltim
                     !          t1(j,k) = t1(j,k) + c(j,k) * deltim
                     !
                  end if
               end do
            end do
         end do

         tem = -1.0/g
         !$acc parallel loop collapse(2) private(delp,i,j,k,dtaux,dtauy,eng0,dbim,eng1,&
         !$acc&                      dutmp,dvtmp,taudtmp) async(async_id)
         do jj = 1, jlistnum
            do i = 1, ix              !org
               if (i .le. npt(jj)) then
                  j = ipt(i, jj)   !org
                  dutmp = 0.
                  dvtmp = 0.
            !!!$acc loop vector reduction(+:dutmp,dvtmp)
                  !$acc loop seq
                  do k = 1, km            !org
                     delp = del(j, k, jj)
                     taudtmp = taud(i, k, jj)*dtfac(i, jj)
                     dtaux = taudtmp*xn(i, jj)
                     dtauy = taudtmp*yn(i, jj)
                     if (k .lt. idxzb(i, jj) .and. idxzb(i, jj) .ne. 0) then
                        dbim = db(i, k, jj)/(1.+db(i, k, jj)*deltim)
                        dutmp = dutmp - dbim*v1(j, k, jj)*delp
                        dvtmp = dvtmp - dbim*u1(j, k, jj)*delp
                     else
                        dutmp = dutmp + dtaux*delp
                        dvtmp = dvtmp + dtauy*delp
                     end if
                  end do
                  dusfc(j, jj) = tem*dutmp
                  dvsfc(j, jj) = tem*dvtmp
                  !dusfc(j,jj) =  dusfc(j,jj)
                  !dvsfc(j,jj) =  dvsfc(j,jj)
               end if
            end do
         end do

!       if (lprnt) then
!         print *,' in gwdps_lm.f after  a=',a(ipr,:)
!         print *,' in gwdps_lm.f after  b=',b(ipr,:)
!         print *,' db=',db(ipr,:)
!       endif
!
!        do k = 1, km
!          do i = 1, im
!            u1(i,k) = u1(i,k) + b(i,k) * deltim
!            v1(i,k) = v1(i,k) + a(i,k) * deltim
!            t1(i,k) = t1(i,k) + c(i,k) * deltim
!          enddo
!        enddo
!
!      monitor for excessive gravity wave drag tendencies if ncnt>0
!
!       if(ncnt.gt.0) then
!          if(lat.ge.38.and.lat.le.42) then
!  cmic$ guard 37
!             do 92 i = 1,im
!                if(ikount.gt.ncnt) go to 92
!                if(i.lt.319.or.i.gt.320) go to 92
!                do 91 k = 1,km
!                   if(abs(taud(i,k)) .gt. critac) then
!                      if(i.le.im) then
!                         ikount = ikount+1
!                         print 123,i,lat,kdt
!                         print 124,taub(i),bnv(i),ulow(i),
!      1                  gf(i),fr(i),roll(i),hprime(i),xn(i),yn(i)
!                         print 124,(taud(i,kk),kk = 1,km)
!                         print 124,(taup(i,kk),kk = 1,km+1)
!                         print 124,(ri_n(i,kk),kk = 1,km)
!                         do 93 kk = 1,kmm1
!                            velko(kk) =
!      1                  0.5*((u1(i,kk)+u1(i,kk+1))*ubar(i)+
!      2                  (v1(i,kk)+v1(i,kk+1))*vbar(i))*uloi(i)
!  93                     continue
!                         print 124,(velko(kk),kk = 1,kmm1)
!                         print 124,(a    (i,kk),kk = 1,km)
!                         print 124,(dtauy(i,kk),kk = 1,km)
!                         print 124,(b    (i,kk),kk = 1,km)
!                         print 124,(dtaux(i,kk),kk = 1,km)
!                         go to 92
!                      endif
!                   endif
!  91            continue
!  92         continue
!  cmic$ end guard 37
!  123        format('  *** migwd print *** i=',i3,' lat=',i3,' kdt=',i3)
!  124        format(2x,  10e13.6)
!          endif
!       endif
!
!        print *,' in gwdps_lm.f 18  =',a(ipt(1),idxzb(1))
!      &,                          b(ipt(1),idxzb(1)),me
!
         !$acc exit data delete(db, ang, uds, taup, npt, ipt, nwdir, &
         !$acc&     iwklm, ro, bnv2lm, delks, delks1, ubar, vbar, &
         !$acc&     roll, bnv2bar, ri_n, bnv2, iwk, kbps, kmps, kref, &
         !$acc&     oa, clx, xn, yn, dtfac, &
         !$acc&     velco, kint, taud,utendform, vtendform) &
         !$acc&                 async(async_id)
         !$acc wait(async_id)

         return
      end

