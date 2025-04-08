!> \file moninedmf.f
!!  Contains most of the hybrid eddy-diffusivity mass-flux scheme except for the
!!  subroutine that calculates the mass flux and updraft properties.

!> \defgroup PBL Hybrid Eddy-diffusivity Mass-flux Scheme
!! @{
!!  \brief The Hybrid EDMF scheme is a first-order turbulent transport scheme used for subgrid-scale vertical turbulent mixing in the PBL and above. It blends the traditional first-order approach that has been used and improved over the last several years
!!
!!  The PBL scheme's main task is to calculate tendencies of temperature, moisture, and momentum due to vertical diffusion throughout the column (not just the PBL). The scheme is an amalgamation of decades of work, starting from the initial first-order PB
!!
!!  \section diagram Calling Hierarchy Diagram
!!  \image html Hybrid_EDMF_Flowchart.png "Diagram depicting how the Hybrid EDMF PBL scheme is called from the GSM physics time loop" height=2cm
!!  \section intraphysics Intraphysics Communication
!!  This space is reserved for a description of how this scheme uses information from other scheme types and/or how information calculated in this scheme is used in other scheme types.

!>  \brief This subroutine contains all of logic for the Hybrid EDMF PBL scheme except for the calculation of the updraft properties and mass flux.
!!
!!  The scheme works on a basic level by calculating background diffusion coefficients and updating them according to which processes are occurring in the column. The most important difference in diffusion coefficients occurs between those levels in the P
!!
!!  \param[in] ix horizontal dimension
!!  \param[in] im number of used points
!!  \param[in] km vertical layer dimension
!!  \param[in] ntrac number of tracers
!!  \param[in] ntcw cloud condensate index in the tracer array
!!  \param[in,out] dv v-momentum tendency (\f$ m s^{-2} \f$)
!!  \param[in,out] du u-momentum tendency (\f$ m s^{-2} \f$)
!!  \param[in,out] tau temperature tendency (\f$ K s^{-1} \f$)
!!  \param[in,out] rtg moisture tendency (\f$ kg kg^{-1} s^{-1} \f$)
!!  \param[in] u1 u component of layer wind (\f$ m s^{-1} \f$)
!!  \param[in] v1 v component of layer wind (\f$ m s^{-1} \f$)
!!  \param[in] t1 layer mean temperature (\f$ K \f$)
!!  \param[in] q1 layer mean tracer concentration (units?)
!!  \param[in] swh total sky shortwave heating rate (\f$ K s^-1 \f$)
!!  \param[in] hlw total sky longwave heating rate (\f$ K s^-1 \f$)
!!  \param[in] xmu time step zenith angle adjust factor for shortwave
!!  \param[in] psk Exner function at surface interface?
!!  \param[in] rbsoil surface bulk Richardson number
!!  \param[in] zorl surface roughness (units?)
!!  \param[in] u10m 10-m u wind (\f$ m s^{-1} \f$)
!!  \param[in] v10m 10-m v wind (\f$ m s^{-1} \f$)
!!  \param[in] fm fm parameter from PBL scheme
!!  \param[in] fh fh parameter from PBL scheme
!!  \param[in] tsea ground surface temperature (K)
!!  \param[in] qss surface saturation humidity (units?)
!!  \param[in] heat surface sensible heat flux (units?)
!!  \param[in] evap evaporation from latent heat flux (units?)
!!  \param[in] stress surface wind stress? (\f$ cm*v^2\f$ in sfc_diff subroutine) (units?)
!!  \param[in] spd1 surface wind speed? (units?)
!!  \param[out] kpbl PBL top index
!!  \param[in] prsi pressure at layer interfaces (units?)
!!  \param[in] del pressure difference between level k and k+1 (units?)
!!  \param[in] prsl mean layer pressure (units?)
!!  \param[in] prslk Exner function at layer
!!  \param[in] phii interface geopotential height (units?)
!!  \param[in] phil layer geopotential height (units?)
!!  \param[in] delt physics time step (s)
!!  \param[in] dspheat flag for TKE dissipative heating
!!  \param[out] dusfc surface u-momentum tendency (units?)
!!  \param[out] dvsfc surface v-momentum tendency (units?)
!!  \param[out] dtsfc surface temperature tendency (units?)
!!  \param[out] dqsfc surface moisture tendency (units?)
!!  \param[out] hpbl PBL top height (m)
!!  \param[out] hgamt counter gradient mixing term for temperature (units?)
!!  \param[out] hgamq counter gradient mixing term for moisture (units?)
!!  \param[out] dkt diffusion coefficient for temperature (units?)
!!  \param[in] kinver index location of temperature inversion
!!  \param[in] xkzm_m background vertical diffusion coefficient for momentum (units?)
!!  \param[in] xkzm_h background vertical diffusion coefficeint for heat, moisture (units?)
!!  \param[in] xkzm_s sigma threshold for background momentum diffusion (units?)
!!  \param[in] lprnt flag to print some output
!!  \param[in] ipr index of point to print
!!
!!  \section general General Algorithm
!!  -# Compute preliminary variables from input arguments.
!!  -# Calculate the first estimate of the PBL height ("Predictor step").
!!  -# Calculate Monin-Obukhov similarity parameters.
!!  -# Update thermal properties of surface parcel and recompute PBL height ("Corrector step").
!!  -# Determine whether stratocumulus layers exist and compute quantities needed for enhanced diffusion.
!!  -# Calculate the inverse Prandtl number.
!!  -# Compute diffusion coefficients below the PBL top.
!!  -# Compute diffusion coefficients above the PBL top.
!!  -# If the PBL is convective, call the mass flux scheme to replace the countergradient terms.
!!  -# Compute enhanced diffusion coefficients related to stratocumulus-topped PBLs.
!!  -# Solve for the temperature and moisture tendencies due to vertical mixing.
!!  -# Calculate heating due to TKE dissipation and add to the tendency for temperature.
!!  -# Solve for the horizontal momentum tendencies and add them to output tendency terms.
!!  \section detailed Detailed Algorithm
!!  @{
      subroutine moninedmf_gpu(ix, myim, km, ntrac, ntcw, &
                               u1, v1, t1, q1, swh, hlw, xmu, &
                               psk, rbsoil, zorl, u10m, v10m, fm, fh, &
                               tsea, heat, evap, stress, spd1, kpbl, &
                               prsi, del, prsl, prslk, phii, phil, delt, &
                               hpbl, handle, async_id )
!
         use machine, only: kind_phys
         use physcons, grav => con_g, rd => con_rd, cp => con_cp &
                                                          , hvap => con_hvap, fv => con_fvirt
         use const, only: RTYPE
         use param, only: my, my_max
         use index, only: jlist1, jlistnum, nxptot, nxjp_acc
         use rank, only: myrank
         use openacc
         use cudafor
         !use cusparse

         implicit none
!
!     arguments
!
         logical lprnt
         integer ipr
         integer ix,  myim(my_max), km, ntrac, ntcw, kpbl(ix, my_max), kinver(ix)
!
         real(kind=kind_phys) delt, xkzm_m, xkzm_h, xkzm_s
         real(kind=kind_phys) tau(ix, km, my_max), &!rtg(ix, km, ntrac, my_max), &
            u1(ix, km, my_max), v1(ix, km, my_max), &
            t1(ix, km, my_max), q1(ix, km, ntrac, my_max), &
            swh(ix, km, my_max), hlw(ix, km, my_max), &
            xmu(ix, my_max), &
            rbsoil(ix, my_max), zorl(ix, my_max), &
            u10m(ix, my_max), v10m(ix, my_max), &
            fm(ix, my_max), fh(ix, my_max), &
            tsea(ix, my_max), &
            spd1(ix, my_max), &
            prsi(ix, km + 1, my_max), del(ix, km, my_max), &
            prsl(ix, km, my_max), prslk(ix, km, my_max), &
            phii(ix, km + 1, my_max), phil(ix, km, my_max), &
            dusfc(ix, my_max), dvsfc(ix, my_max), &
            dtsfc(ix, my_max), dqsfc(ix, my_max), &
            hpbl(ix, my_max), hpblx(ix, my_max), &
            hgamt(ix, my_max), hgamq(ix, my_max)
         real(kind=RTYPE) psk(ix, km, my_max)
!
         logical dspheat
!          flag for tke dissipative heating
!
!    locals
!
         integer i, iprt, is, iun, k, kk, km1, kmpbl, latd, lond, jj, j1, m
         integer lcld(ix, my_max), icld(ix, my_max), kcld(ix, my_max), krad(ix, my_max)
         integer kx1(ix, my_max), kpblx(ix, my_max)
!
!     real(kind=kind_phys) betaq(im), betat(im),   betaw(im),
         real(kind=kind_phys) evap(ix, my_max), heat(ix, my_max), phih(ix, my_max), &
            phim(ix, my_max), rbdn(ix, my_max), rbup(ix, my_max), &
            stress(ix, my_max), beta(ix, my_max), sflux(ix, my_max), &
            z0(ix, my_max), crb(ix, my_max), wstar(ix, my_max), &
            zol(ix, my_max), ustmin(ix, my_max), ustar(ix, my_max), &
            thermal(ix, my_max), wscale(ix, my_max), wscaleu(ix, my_max)
!
         real(kind=kind_phys) theta(ix, km, my_max), thvx(ix, km, my_max), thlvx(ix, km, my_max), &
            qlx(ix, km, my_max), thetae(ix, km, my_max), &
            qtx(ix, km, my_max), bf(ix, km - 1, my_max), diss(ix, km, my_max), &
            radx(ix, km - 1, my_max), &
            govrth(ix, my_max), hrad(ix, my_max), &
            !    &                     hradm(im),   radmin(im),   vrad(im),         &
            radmin(ix, my_max), vrad(ix, my_max), &
            zd(ix, my_max), zdd(ix, my_max), thlvx1(ix, my_max)
!
         real(kind=kind_phys) rdzt(ix, km - 1, my_max), dktx(ix, km - 1, my_max), &
            zi(ix, km + 1, my_max), zl(ix, km, my_max), xkzo(ix, km - 1, my_max), &
            dku(ix, km - 1, my_max), dkt(ix, km - 1, my_max), xkzmo(ix, km - 1, my_max), &
            cku(ix, km - 1, my_max), ckt(ix, km - 1, my_max), &
            !al(ix,km-1,my_max),  ad(ix,km,my_max),                     &
            !au(ix,km-1,my_max),  a1(ix,km,my_max),                     &
            !a2(ix,km*ntrac,my_max), &
            ti(ix, km - 1, my_max), shr2(ix, km - 1, my_max)
!
         real(kind=kind_phys) tcko(ix, km, my_max), qcko(ix, km, ntrac, my_max), &
            ucko(ix, km, my_max), vcko(ix, km, my_max), xmf(ix, km, my_max)
!
         real(kind=kind_phys) prinv(ix, my_max), rent(ix, my_max)
         real(kind=kind_phys) alg(nxptot, km), adg(nxptot, km), &
            aug(nxptot, km), a3g(nxptot, km, ntrac+1)
            
         real(kind=kind_phys) a2g(nxptot, km, ntrac), aug_backup(nxptot, km), &
            a1g(nxptot, km) ! cusparse

!
         logical pblflg(ix, my_max), sfcflg(ix, my_max), scuflg(ix, my_max)
         logical ublflg(ix, my_max), pcnvflg(ix, my_max)
!
!  pcnvflg: true for convective(strongly unstable) pbl
!  ublflg: true for unstable but not convective(strongly unstable) pbl
!
         real(kind=kind_phys) aphi16, aphi5, bvf2, wfac, &
            cfac, conq, cont, conw, &
            dk, dkmax, dkmin, &
            dq1, dsdz2, dsdzq, dsdzt, &
            dsdzu, dsdzv, &
            dsig, dt2, dthe1, dtodsd, &
            dtodsu, dw2, dw2min, g, &
            gamcrq, gamcrt, gocp, &
            gravi, f0, &
            prnum, prmax, prmin, pfac, crbcon, &
            qmin, tdzmin, qtend, crbmin, crbmax, &
            rbint, rdt, rdz, qlmin, &
            ri, rimin, rl2, rlam, rlamun, &
            rone, rzero, sfcfrac, &
            spdk2, sri, zol1, zolcr, zolcru, &
            robn, ttend, &
            utend, vk, vk2, &
            ust3, wst3, &
            vtend, zfac, vpert, cteit, &
            rentf1, rentf2, radfac, &
            zfmin, zk, tem, tem1, tem2, &
            xkzm, xkzmu, xkzminv, &
            ptem, ptem1, ptem2, tx1(ix, my_max), tx2(ix, my_max)
!
         real(kind=kind_phys) moninq_fac
!
         real(kind=kind_phys) zstblmax, h1, h2, qlcr, actei, &
            cldtime
         ! registers for GPU
         real(kind=kind_phys) dtodsu_1, rdz_1, dsig_1, tem1_1, dsdz2_1, al_1, &
            tem2_1, ptem_1, ptem2_1, dttmp, dqtmp, tem3, &
            time1, time2, aur, alr, k1, q1r, t1r, thetar, a3gr1, a3gr2, &
            crbr, rbupr, rbupr1, kpblr, rbdnr, thvxr, zolr, shr2r, adgr, &
            wscaleur, wscaler, ustarr, dkur, dktr, kradr, qtxr, tx1r, tx2r
         !type(cusparseHandle) :: handle
         integer :: handle
         logical flgr
         integer :: accui, ii
!cc
         parameter(gravi=1.0/grav)
         parameter(g=grav)
         parameter(gocp=g/cp)
         parameter(cont=cp/g, conq=hvap/g, conw=1.0/g)               ! for del in pa
!     parameter(cont=1000.*cp/g,conq=1000.*hvap/g,conw=1000./g) ! for del in kpa
         parameter(rlam=30.0, vk=0.4, vk2=vk*vk)
         parameter(prmin=0.25, prmax=4., zolcr=0.2, zolcru=-0.5)
         parameter(dw2min=0.0001, dkmin=0.0, dkmax=1000., rimin=-100.)
         parameter(crbcon=0.25, crbmin=0.15, crbmax=0.35)
         parameter(wfac=7.0, cfac=6.5, pfac=2.0, sfcfrac=0.1)
!     parameter(qmin=1.e-8,xkzm=1.0,zfmin=1.e-8,aphi5=5.,aphi16=16.)
         parameter(qmin=1.e-8, zfmin=1.e-8, aphi5=5., aphi16=16.)
         parameter(tdzmin=1.e-3, qlmin=1.e-12, f0=1.e-4)
         parameter(h1=0.33333333, h2=0.66666667)
!     parameter(cldtime=500.,xkzminv=0.3)
         parameter(cldtime=500.)
!     parameter(cldtime=500.,xkzmu=3.0,xkzminv=0.3)
!     parameter(gamcrt=3.,gamcrq=2.e-3,rlamun=150.0)
         parameter(xkzm_m=1.0, xkzm_h=1.0, xkzm_s=1.0, xkzminv=0.3, moninq_fac=1.0)
         parameter(gamcrt=3., gamcrq=0., rlamun=150.0)
         parameter(rentf1=0.2, rentf2=1.0, radfac=0.85)
         parameter(iun=84)
!
!     parameter (zstblmax = 2500., qlcr=1.0e-5)
!     parameter (zstblmax = 2500., qlcr=3.0e-5)
!     parameter (zstblmax = 2500., qlcr=3.5e-5)
!     parameter (zstblmax = 2500., qlcr=1.0e-4)
         parameter(zstblmax=2500., qlcr=3.5e-5)
!     parameter (actei = 0.23)
         parameter(actei=0.7)
         integer(kind=cuda_stream_kind) stream
         integer istat, async_id

         stream = acc_get_cuda_stream(async_id)

!
!byl      dspheat=.false. reference from NCEP fv3GFS, it should be .true.
         dspheat = .true.
!
!c-----------------------------------------------------------------------
!
601      format(1x, ' moninp lat lon step hour ', 3i6, f6.1)
602      format(1x, '    k', '        z', '        t', '       th', &
                '      tvh', '        q', '        u', '        v', &
                '       sp')
603      format(1x, i5, 8f9.1)
604      format(1x, '  sfc', 9x, f9.1, 18x, f9.1)
605      format(1x, '    k      zl    spd2   thekv   the1v' &
                , ' thermal    rbup')
606      format(1x, i5, 6f8.2)
607      format(1x, ' kpbl    hpbl      fm      fh   hgamt', &
                '   hgamq      ws   ustar      cd      ch')
608      format(1x, i5, 9f8.2)
609      format(1x, ' k pr dkt dku ', i5, 3f8.2)
610      format(1x, ' k pr dkt dku ', i5, 3f8.2, ' l2 ri t2', &
                ' sr2  ', 2f8.2, 2e10.2)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!>  ## Compute preliminary variables from input arguments
!

         !!$acc enter data create(a1g, a2g, aug_backup) async(async_id) !cusparse

         !!acc enter data create(rtg) async(async_id)
         !!$acc host_data use_device(rtg)
         !istat = cudaMemsetAsync(rtg, 0.0, size(rtg), stream)
         !!$acc end host_data

!
! compute preliminary variables
!
         !if (ix .lt. im) stop
!
!     iprt = 0
!     if(iprt.eq.1) then
!cc   latd = 0
!     lond = 0
!     else
!cc   latd = 0
!     lond = 0
!     endif
!
         dt2 = delt
         rdt = 1./dt2
         km1 = km - 1
         kmpbl = km/2

         !$acc enter data create(zi, zl) async(async_id)
         !>  - Compute physical height of the layer centers and interfaces from the geopotential height (zi and zl)
         !$acc parallel loop gang collapse(2) private(jj,k,i) async(async_id)
         do jj = 1, jlistnum
            do k = 1, km
               !$acc loop vector
               do i = 1, myim(jj)
                  zi(i, k, jj) = phii(i, k, jj)*gravi
                  zl(i, k, jj) = phil(i, k, jj)*gravi
                  if (k .eq. 1) zi(i, km + 1, jj) = phii(i, km + 1, jj)*gravi
               end do
            end do
         end do


         !$acc enter data create(ckt, cku, dku, dkt, dktx, qlx, qtx, radx, &
         !$acc&      theta, thetae, thlvx, thvx) async(async_id)
         !$acc parallel loop gang collapse(2) private(jj,k,i,ptem, ptem1, ptem2, &
         !$acc&     q1r,t1r,thetar,qtxr) async(async_id)
         do jj = 1, jlistnum
            do k = 1, km
                !$acc loop vector
                do i = 1, myim(jj)
         !>  - Compute \f$\theta\f$ (theta), \f$q_l\f$ (qlx), \f$q_t\f$ (qtx), \f$\theta_e\f$ (thetae), \f$\theta_v\f$ (thvx), \f$\theta_{l,v}\f$ (thlvx)
                     q1r = q1(i, k, 1, jj)
                     t1r = t1(i, k, jj)
                     thetar = t1r*psk(i, km, jj)/prslk(i, k, jj)
                     ptem = max(q1(i, k, ntcw, jj), qlmin)
                     qtxr = max(q1r, qmin) + ptem
                     ptem1 = hvap*max(q1r, qmin)/(cp*t1r)
                     thetae(i, k, jj) = thetar*(1.+ptem1)
                     thvx(i, k, jj) = thetar*(1.+fv*max(q1r, qmin) - ptem)
                     ptem2 = thetar - (hvap/cp)*ptem
                     thlvx(i, k, jj) = ptem2*(1.+fv*qtxr)
                     if (k .ne. km) then
!>  - Initialize diffusion coefficients to 0 and calculate the total radiative heating rate (dku, dkt, radx)
                        dku(i, k, jj) = 0.
                        dkt(i, k, jj) = 0.
                        dktx(i, k, jj) = 0.
                        cku(i, k, jj) = 0.
                        ckt(i, k, jj) = 0.
                        tem = zi(i, k + 1, jj) - zi(i, k, jj)
                        radx(i, k, jj) = tem*(swh(i, k, jj)*xmu(i, jj) + hlw(i, k, jj))
                     end if
                     theta(i, k, jj) = thetar
                     qlx(i, k, jj) = ptem
                     qtx(i, k, jj) = qtxr
               end do
            end do
         end do
         !$acc enter data create(beta,bf,diss, &
         !$acc&                  dusfc,dvsfc,dtsfc,dqsfc, &
         !$acc&                  crb,govrth,hgamt,hgamq, &
         !$acc&                  hpblx,hrad,icld,kcld,krad,kpblx, &
         !$acc&                  kx1,lcld,pblflg,pcnvflg,phih,phim, &
         !$acc&                  prinv,qcko,radmin,rent,rbdn, &
         !$acc&                  rbup,rdzt,sfcflg,scuflg,shr2,sflux, &
         !$acc&                  tau,tcko,thermal, &
         !$acc&                  thlvx1,ti,ublflg,ucko, &
         !$acc&                  ustar,ustmin,vcko,vrad,wscale,wscaleu,wstar, &
         !$acc&                  xmf,z0,zd,zdd,zol,xkzo,xkzmo,tx1,tx2, &
         !$acc&                  alg,adg,aug,a3g) async(async_id)
         
         !$acc parallel loop gang vector collapse(2) private(jj,k,i,ptem,tem1,tx1r, tx2r) async(async_id)
         do jj = 1, jlistnum
            do i = 1, ix
               if (i .le. myim(jj)) then
         !>  - Compute reciprocal of pressure (tx1, tx2)
                  kx1(i, jj) = 1
                  tx1r = 1.0/prsi(i, 1, jj)
                  tx2r = tx1r
         !>  - Compute background vertical diffusivities for scalars and momentum (xkzo and xkzmo)
                  !$acc loop seq
                  do k = 1, km1
                     !>  - Compute reciprocal of \f$ \Delta z \f$ (rdzt)
                     rdzt(i, k, jj) = 1.0/(zl(i, k + 1, jj) - zl(i, k, jj))
                     !byl          if (k < kinver(i)) then
                     !                                  vertical background diffusivity
                     ptem = prsi(i, k + 1, jj)*tx1r
                     tem1 = 1.0 - ptem
                     tem1 = tem1*tem1*10.0
                     xkzo(i, k, jj) = xkzm_h*min(1.0, exp(-tem1))

                     !                                  vertical background diffusivity for momentum
                     if (ptem >= xkzm_s) then
                        xkzmo(i, k, jj) = xkzm_m
                        kx1(i, jj) = k + 1
                     else
                        if (k == kx1(i, jj) .and. k > 1) tx2r = 1.0/prsi(i, k, jj)
                        tem1 = 1.0 - prsi(i, k + 1, jj)*tx2r
                        tem1 = tem1*tem1*5.0
                        xkzmo(i, k, jj) = xkzm_m*min(1.0, exp(-tem1))
                     end if
                     !byl          endif
                  end do
                  tx1(i, jj) = tx1r
                  tx2(i, jj) = tx2r
               end if
            end do
         end do
         !     if (lprnt) then
         !       print *,' xkzo=',(xkzo(ipr,k),k=1,km1)
         !       print *,' xkzmo=',(xkzmo(ipr,k),k=1,km1)
         !     endif
         !
         
         ! diffusivity in the inversion layer is set to be xkzminv (m^2/s)
         !>  - The background scalar vertical diffusivity is limited to be less than or equal to xkzminv
         !$acc parallel loop gang collapse(2) private(jj,k,i,tem1,rdz,dw2) async(async_id)
         do jj = 1, jlistnum
            do k = 1, km1
               !$acc loop vector
               do i = 1, myim(jj)
                  rdz = rdzt(i, k, jj)
                  bf(i, k, jj) = (thvx(i, k + 1, jj) - thvx(i, k, jj))*rdz
                  ti(i, k, jj) = 2./(t1(i, k, jj) + t1(i, k + 1, jj))
                  dw2 = (u1(i, k, jj) - u1(i, k + 1, jj))**2 &
                        + (v1(i, k, jj) - v1(i, k + 1, jj))**2
                  shr2(i, k, jj) = max(dw2, dw2min)*rdz*rdz
                  if (k .le. kmpbl) then
                     !         if(zi(i,k+1) > 200..and.zi(i,k+1) < zstblmax) then
                     if (zi(i, k + 1, jj) > 250.) then
                        tem1 = (t1(i, k + 1, jj) - t1(i, k, jj))*rdz
                        if (tem1 > 1.e-5) then
                           xkzo(i, k, jj) = min(xkzo(i, k, jj), xkzminv)
                        end if
                     end if
                  end if
               end do
            end do
         end do

         !$acc parallel loop gang vector collapse(2) &
         !$acc&              private(jj,k,i,spdk2,tem,robn,tem1,crbr,rbupr, &
         !$acc&                      rbupr1,kpblr,rbdnr,flgr,thvxr,rbint,zol1, &
         !$acc&                      wst3,ust3,vpert,zolr,wscaleur,wscaler,&
         !$acc&                      ustarr,kk,kradr) async(async_id)
         do jj = 1, jlistnum
            do i = 1, ix
               if (i .le. myim(jj)) then
         !>  - Some output variables and logical flags are initialized
                  z0(i, jj) = 0.01*zorl(i, jj)
                  dusfc(i, jj) = 0.
                  dvsfc(i, jj) = 0.
                  dtsfc(i, jj) = 0.
                  dqsfc(i, jj) = 0.
                  wscale(i, jj) = 0.
                  wscaleu(i, jj) = 0.
                  kpbl(i, jj) = 1
                  hpbl(i, jj) = zi(i, 1, jj)
                  hpblx(i, jj) = zi(i, 1, jj)
                  pblflg(i, jj) = .true.
                  sfcflg(i, jj) = .true.
                  if (rbsoil(i, jj) > 0.) sfcflg(i, jj) = .false.
                  ublflg(i, jj) = .false.
                  pcnvflg(i, jj) = .false.
                  scuflg(i, jj) = .true.
                  if (scuflg(i, jj)) then
                     radmin(i, jj) = 0.
                     rent(i, jj) = rentf1
                     hrad(i, jj) = zi(i, 1, jj)
                     !          hradm(i) = zi(i,1)
                     krad(i, jj) = 1
                     icld(i, jj) = 0
                     lcld(i, jj) = km1
                     kcld(i, jj) = km1
                     zd(i, jj) = 0.
                  end if
                  
         !  compute virtual potential temp gradient (bf) and winshear square
         !>  - Compute \f$\frac{\partial \theta_v}{\partial z}\f$ (bf) and the wind shear squared (shr2)
                  govrth(i, jj) = g/theta(i, 1, jj)
                  beta(i, jj) = dt2/(zi(i, 2, jj) - zi(i, 1, jj))
                  ustar(i, jj) = sqrt(stress(i, jj))
                  sflux(i, jj) = heat(i, jj) + evap(i, jj)*fv*theta(i, 1, jj)
                  if (.not. sfcflg(i, jj) .or. sflux(i, jj) <= 0.) pblflg(i, jj) = .false.

                  !>  - Set lcld to first index above 2.5km
                  flgr = scuflg(i, jj)
                  !$acc loop seq
                  do k = 1, km1
                     if (flgr .and. zl(i, k, jj) >= zstblmax) then
                        lcld(i, jj) = k
                        flgr = .false.
                     end if
                  end do
!---------------------------------------------
                  !>  - Calculate \f$\frac{g}{\theta}\f$ (govrth), \f$\beta = \frac{\Delta t}{\Delta z}\f$ (beta), \f$u_*\f$ (ustar), total surface flux (sflux), and set pblflag to false if the total surface energy flux is into the surface
                  !>  ## Calculate the first estimate of the PBL height (``Predictor step")
  !!  The calculation of the boundary layer height follows Troen and Mahrt (1986) \cite troen_and_mahrt_1986 section 3. The approach is to find the level in the column where a modified bulk Richardson number exceeds a critical value.
  !!
  !!  The temperature of the thermal is of primary importance. For the initial estimate of the PBL height, the thermal is assumed to have one of two temperatures. If the boundary layer is stable, the thermal is assumed to have a temperature equal to the sur
                  !  compute the pbl height
                  !

                  flgr = .false.
                  rbupr = rbsoil(i, jj)
                  thvxr = thvx(i, 1, jj)
                  !
                  if (pblflg(i, jj)) then
                     thermal(i, jj) = thvxr
                     crbr = crbcon
                  else
                     thermal(i, jj) = tsea(i, jj)*(1.+fv*max(q1(i, 1, 1, jj), qmin))
                     tem = sqrt(u10m(i, jj)**2 + v10m(i, jj)**2)
                     tem = max(tem, 1.)
                     robn = tem/(f0*z0(i, jj))
                     tem1 = 1.e-7*robn
                     crbr = 0.16*(tem1**(-0.18))
                     crbr = max(min(crbr, crbmax), crbmin)
                  end if
                  !>  Given the thermal's properties and the critical Richardson number, a loop is executed to find the first level above the surface where the modified Richardson number is greater than the critical Richardson number, using equation 10a from Troen and Mahr
  !!  \f[
  !!  h = Ri\frac{T_0\left|\vec{v}(h)\right|^2}{g\left(\theta_v(h) - \theta_s\right)}
  !!  \f]
  !!  where \f$h\f$ is the PBL height, \f$Ri\f$ is the Richardson number, \f$T_0\f$ is the virtual potential temperature near the surface, \f$\left|\vec{v}\right|\f$ is the wind speed, and \f$\theta_s\f$ is for the thermal. Rearranging this equation to calc
  !!  \f[
  !!  Ri_k = gz(k)\frac{\left(\theta_v(k) - \theta_s\right)}{\theta_v(1)*\vec{v}(k)}
  !!  \f]
                  !$acc loop seq
                  do k = 1, kmpbl
                     if (.not. flgr) then
                        rbupr1 = rbupr
                        spdk2 = max((u1(i, k, jj)**2 + v1(i, k, jj)**2), 1.)
                        rbupr = (thvx(i, k, jj) - thermal(i, jj))* &
                                (g*zl(i, k, jj)/thvxr)/spdk2
                        kpblr = k
                        flgr = rbupr > crbr
                     end if
                  end do
                  crb(i, jj) = crbr
                  rbup(i, jj) = rbupr
                  rbdn(i, jj) = rbupr1
                  kpbl(i, jj) = kpblr

!------------------------------------
                  !>  Once the level is found, some linear interpolation is performed to find the exact height of the boundary layer top (where \f$Ri = Ri_{cr}\f$) and the PBL height and the PBL top index are saved (hpblx and kpblx, respectively)
                  if (kpbl(i, jj) > 1) then
                     k = kpbl(i, jj)
                     if (rbdn(i, jj) >= crb(i, jj)) then
                        rbint = 0.
                     elseif (rbup(i, jj) <= crb(i, jj)) then
                        rbint = 1.
                     else
                        rbint = (crb(i, jj) - rbdn(i, jj))/(rbup(i, jj) - rbdn(i, jj))
                     end if
                     hpbl(i, jj) = zl(i, k - 1, jj) + rbint*(zl(i, k, jj) - zl(i, k - 1, jj))
                     if (hpbl(i, jj) < zi(i, kpbl(i, jj), jj)) kpbl(i, jj) = kpbl(i, jj) - 1
                  else
                     hpbl(i, jj) = zl(i, 1, jj)
                     kpbl(i, jj) = 1
                  end if
                  kpblx(i, jj) = kpbl(i, jj)
                  hpblx(i, jj) = hpbl(i, jj)
!--------------------------------------------
                  !
                  !  compute similarity parameters
                  !>  ## Calculate Monin-Obukhov similarity parameters
    !!  Using the initial guess for the PBL height, Monin-Obukhov similarity parameters are calculated. They are needed to refine the PBL height calculation and for calculating diffusion coefficients.
    !!
    !!  First, calculate the Monin-Obukhov nondimensional stability parameter, commonly referred to as \f$\zeta\f$ using the following equation from Businger et al. (1971) \cite businger_et_al_1971 (equation 28):
    !!  \f[
    !!  \zeta = Ri_{sfc}\frac{F_m^2}{F_h} = \frac{z}{L}
    !!  \f]
    !!  where \f$F_m\f$ and \f$F_h\f$ are surface Monin-Obukhov stability functions calculated in sfc_diff.f and \f$L\f$ is the Obukhov length. Then, the nondimensional gradients of momentum and temperature (phim and phih) are calculated using equations 5 and
    !!  \f[
    !!  w_* = \left(\frac{g}{\theta_0}h\overline{w'\theta_0'}\right)^{1/3}
    !!  \f]
    !!  and the mixed layer velocity scale is then calculated with equation 6 from Troen and Mahrt (1986) \cite troen_and_mahrt_1986
    !!  \f[
    !!  w_s = (u_*^3 + 7\epsilon k w_*^3)^{1/3}
    !!  \f]
                  zolr = max(rbsoil(i, jj)*fm(i, jj)*fm(i, jj)/fh(i, jj), rimin)
                  ustarr = ustar(i, jj)
                  if (sfcflg(i, jj)) then
                     zolr = min(zolr, -zfmin)
                  else
                     zolr = max(zolr, zfmin)
                  end if
                  zol1 = zolr*sfcfrac*hpbl(i, jj)/zl(i, 1, jj)
                  if (sfcflg(i, jj)) then
                     !          phim(i) = (1.-aphi16*zol1)**(-1./4.)
                     !          phih(i) = (1.-aphi16*zol1)**(-1./2.)
                     tem = 1.0/(1.-aphi16*zol1)
                     phih(i, jj) = sqrt(tem)
                     phim(i, jj) = sqrt(phih(i, jj))
                  else
                     phim(i, jj) = 1.+aphi5*zol1
                     phih(i, jj) = phim(i, jj)
                  end if
                  wscaler = ustarr/phim(i, jj)
                  ustmin(i, jj) = ustarr/aphi5
                  wscale(i, jj) = max(wscaler, ustmin(i, jj))
                  zol(i, jj) = zolr

!---------------------------------------
                  if (pblflg(i, jj)) then
                     if (zolr < zolcru .and. kpbl(i, jj) > 1) then
                        pcnvflg(i, jj) = .true.
                     else
                        ublflg(i, jj) = .true.
                     end if
                     wst3 = govrth(i, jj)*sflux(i, jj)*hpbl(i, jj)
                     wstar(i, jj) = wst3**h1
                     ust3 = ustarr**3.
                     wscaleur = (ust3 + wfac*vk*wst3*sfcfrac)**h1
                     wscaleur = max(wscaleur, ustmin(i, jj))
                     wscaleu(i, jj) = wscaleur
                  end if

!------------------------------------------
                  !
                  ! compute counter-gradient mixing term for heat and moisture
                  !>  ## Update thermal properties of surface parcel and recompute PBL height ("Corrector step").
    !!  Next, the counter-gradient terms for temperature and humidity are calculated using equation 4 of Hong and Pan (1996) \cite hong_and_pan_1996 and are used to calculate the "scaled virtual temperature excess near the surface" (equation 9 in Hong and Pan
                  if (ublflg(i, jj)) then
                     hgamt(i, jj) = min(cfac*heat(i, jj)/wscaleur, gamcrt)
                     hgamq(i, jj) = min(cfac*evap(i, jj)/wscaleur, gamcrq)
                     vpert = hgamt(i, jj) + hgamq(i, jj)*fv*theta(i, 1, jj)
                     vpert = min(vpert, gamcrt)
                     thermal(i, jj) = thermal(i, jj) + max(vpert, 0.)
                     hgamt(i, jj) = max(hgamt(i, jj), 0.0)
                     hgamq(i, jj) = max(hgamq(i, jj), 0.0)
                  end if

!--------------------------------------
                  !  enhance the pbl height by considering the thermal excess
                  !>  The PBL height calculation follows the same procedure as the predictor step, except that it uses an updated virtual potential temperature for the thermal.
                  flgr = .true.
                  if (ublflg(i, jj)) then
                     flgr = .false.
                     rbupr = rbsoil(i, jj)
                  end if
                  crbr = crb(i, jj)
                  !$acc loop seq
                  do k = 2, kmpbl
                     if (.not. flgr) then
                        rbupr1 = rbupr
                        spdk2 = max((u1(i, k, jj)**2 + v1(i, k, jj)**2), 1.)
                        rbupr = (thvx(i, k, jj) - thermal(i, jj))* &
                                (g*zl(i, k, jj)/thvxr)/spdk2
                        kpbl(i, jj) = k
                        flgr = rbupr > crbr
                     end if
                  end do
!------------------------------------------
                  if (ublflg(i, jj)) then
                     rbup(i, jj) = rbupr
                     rbdn(i, jj) = rbupr1
                     k = kpbl(i, jj)
                     if (rbdn(i, jj) >= crb(i, jj)) then
                        rbint = 0.
                     elseif (rbupr <= crb(i, jj)) then
                        rbint = 1.
                     else
                        rbint = (crb(i, jj) - rbdn(i, jj))/(rbupr - rbdn(i, jj))
                     end if
                     hpbl(i, jj) = zl(i, k - 1, jj) + rbint*(zl(i, k, jj) - zl(i, k - 1, jj))
                     if (hpbl(i, jj) < zi(i, kpbl(i, jj), jj)) kpbl(i, jj) = kpbl(i, jj) - 1
                     if (kpbl(i, jj) <= 1) then
                        ublflg(i, jj) = .false.
                        pblflg(i, jj) = .false.
                     end if
                  end if
!------------------------------------------
                  !
                  !  look for stratocumulus
                  !>  ## Determine whether stratocumulus layers exist and compute quantities needed for enhanced diffusion
  !!  - Starting at the PBL top and going downward, if the level is less than 2.5 km and \f$q_l>q_{l,cr}\f$ then set kcld = k (find the cloud top index in the PBL). If no cloud water above the threshold is found, scuflg is set to F.
                  flgr = scuflg(i, jj)
                  !$acc loop seq
                  do k = kmpbl, 1, -1
                     if (flgr .and. k <= lcld(i, jj)) then
                        if (qlx(i, k, jj) .ge. qlcr) then
                           kcld(i, jj) = k
                           flgr = .false.
                        end if
                     end if
                  end do
                  if (scuflg(i, jj) .and. kcld(i, jj) == km1) scuflg(i, jj) = .false.
!------------------------------------------
                  !>  - Starting at the PBL top and going downward, if the level is less than the cloud top, find the level of the minimum radiative heating rate within the cloud. If the level of the minimum is the lowest model level or the minimum radiative heating rate i
                  flgr = scuflg(i, jj)
                  !$acc loop seq
                  do k = kmpbl, 1, -1
                     if (flgr .and. k <= kcld(i, jj)) then
                        if (qlx(i, k, jj) >= qlcr) then
                           if (radx(i, k, jj) < radmin(i, jj)) then
                              radmin(i, jj) = radx(i, k, jj)
                              kradr = k
                           end if
                        else
                           flgr = .false.
                        end if
                     end if
                  end do
                  if (scuflg(i, jj) .and. kradr <= 1) scuflg(i, jj) = .false.
                  if (scuflg(i, jj) .and. radmin(i, jj) >= 0.) scuflg(i, jj) = .false.
!------------------------------------------
                  !>  - Starting at the PBL top and going downward, count the number of levels below the minimum radiative heating rate level that have cloud water above the threshold. If there are none, then set the scuflg to F.
                  flgr = scuflg(i, jj)
                  !$acc loop seq
                  do k = kmpbl, 2, -1
                     if (flgr .and. k <= kradr) then
                        if (qlx(i, k, jj) >= qlcr) then
                           icld(i, jj) = icld(i, jj) + 1
                        else
                           flgr = .false.
                        end if
                     end if
                  end do
                  if (scuflg(i, jj) .and. icld(i, jj) < 1) scuflg(i, jj) = .false.

!------------------------------------------
                  !>  - Find the height of the interface where the minimum in radiative heating rate is located. If this height is less than the second model interface height, then set the scuflg to F.
                  if (scuflg(i, jj)) then
                     hrad(i, jj) = zi(i, kradr + 1, jj)
                     !          hradm(i)= zl(i,krad(i))
                  end if
                  if (scuflg(i, jj) .and. hrad(i, jj) < zi(i, 2, jj)) scuflg(i, jj) = .false.
                  !>  - Calculate the hypothetical \f$\theta_v\f$ at the minimum radiative heating level that a parcel would reach due to radiative cooling after a typical cloud turnover time spent at that level.
                  if (scuflg(i, jj)) then
                     k = kradr
                     tem = zi(i, k + 1, jj) - zi(i, k, jj)
                     tem1 = cldtime*radmin(i, jj)/tem
                     thlvx1(i, jj) = thlvx(i, k, jj) + tem1
                     !         if(thlvx1(i) > thlvx(i,k-1)) scuflg(i)=.false.
                  end if
                  !-----------------------------------------------
                  !>  - Determine the distance that a parcel would sink downwards starting from the level of minimum radiative heating rate by comparing the hypothetical minimum \f$\theta_v\f$ calculated above with the environmental \f$\theta_v\f$.
                  flgr = scuflg(i, jj)
                  !$acc loop seq
                  do k = kmpbl, 1, -1
                     if (flgr .and. k <= kradr) then
                        if (thlvx1(i, jj) <= thlvx(i, k, jj)) then
                           tem = zi(i, k + 1, jj) - zi(i, k, jj)
                           zd(i, jj) = zd(i, jj) + tem
                        else
                           flgr = .false.
                        end if
                     end if
                  end do
                  !>  - Calculate the cloud thickness, where the cloud top is the in-cloud minimum radiative heating level and the bottom is determined previously.
                  if (scuflg(i, jj)) then
                     kk = max(1, kradr + 1 - icld(i, jj))
                     zdd(i, jj) = hrad(i, jj) - zi(i, kk, jj)
                  end if
                  !------------------------------------------------

                  !>  - Find the largest between the cloud thickness and the distance of a sinking parcel, then determine the smallest of that number and the height of the minimum in radiative heating rate. Set this number to \f$zd\f$. Using \f$zd\f$, calculate the charact
                  if (scuflg(i, jj)) then
                     zd(i, jj) = max(zd(i, jj), zdd(i, jj))
                     zd(i, jj) = min(zd(i, jj), hrad(i, jj))
                     tem = govrth(i, jj)*zd(i, jj)*(-radmin(i, jj))
                     vrad(i, jj) = tem**h1
                  end if
                  !     compute inverse prandtl number
                  !>  ## Calculate the inverse Prandtl number
    !!  For an unstable PBL, the Prandtl number is calculated according to Hong and Pan (1996) \cite hong_and_pan_1996, equation 10, whereas for a stable boundary layer, the Prandtl number is simply \f$Pr = \frac{\phi_h}{\phi_m}\f$.
                  if (ublflg(i, jj)) then
                     tem = phih(i, jj)/phim(i, jj) + cfac*vk*sfcfrac
                  else
                     tem = phih(i, jj)/phim(i, jj)
                  end if
                  prinv(i, jj) = 1.0/tem
                  prinv(i, jj) = min(prinv(i, jj), prmax)
                  prinv(i, jj) = max(prinv(i, jj), prmin)
                  if (zol(i, jj) > zolcr) then
                     kpbl(i, jj) = 1
                  end if
                  krad(i, jj) = kradr
               end if !myim
            end do
         end do
         !
         !
         !     compute diffusion coefficients below pbl
         !>  ## Compute diffusion coefficients below the PBL top
  !!  Below the PBL top, the diffusion coefficients (\f$K_m\f$ and \f$K_h\f$) are calculated according to equation 2 in Hong and Pan (1996) \cite hong_and_pan_1996 where a different value for \f$w_s\f$ (PBL vertical velocity scale) is used depending on the
         !$acc parallel loop gang vector collapse(3) private(jj,k,i,zfac,tem,tem1, &
         !$acc&                               dkur,dktr) async(async_id)
         do jj = 1, jlistnum
            do k = 1, kmpbl
               do i = 1, ix
                  if (i .le. myim(jj)) then
                     if (k < kpbl(i, jj)) then
                        !           zfac = max((1.-(zi(i,k+1)-zl(i,1))/
                        !    1             (hpbl(i)-zl(i,1))), zfmin)
                        zfac = max((1.-zi(i, k + 1, jj)/hpbl(i, jj)), zfmin)
                        tem = zi(i, k + 1, jj)*(zfac**pfac)*moninq_fac ! lmh suggested by kg
                        if (pblflg(i, jj)) then
                           tem1 = vk*wscaleu(i, jj)*tem
                           !             dku(i,k) = xkzmo(i,k) + tem1
                           !             dkt(i,k) = xkzo(i,k)  + tem1 * prinv(i)
                           dkur = tem1
                           dktr = tem1*prinv(i, jj)
                        else
                           tem1 = vk*wscale(i, jj)*tem
                           !             dku(i,k) = xkzmo(i,k) + tem1
                           !             dkt(i,k) = xkzo(i,k)  + tem1 * prinv(i)
                           dkur = tem1
                           dktr = tem1*prinv(i, jj)
                        end if
                        dkur = min(dkur, dkmax)
                        dkur = max(dkur, xkzmo(i, k, jj))
                        dktr = min(dktr, dkmax)
                        dktr = max(dktr, xkzo(i, k, jj))
                        dktx(i, k, jj) = dktr
                        dku(i, k, jj) = dkur
                        dkt(i, k, jj) = dktr
                     end if
                  end if
               end do
            end do
         end do

         !
         ! compute diffusion coefficients based on local scheme above pbl
         !>  ## Compute diffusion coefficients above the PBL top
  !!  Diffusion coefficients above the PBL top are computed as a function of local stability (gradient Richardson number), shear, and a length scale from Louis (1979) \cite louis_1979 :
  !!  \f[
  !!  K_{m,h}=l^2f_{m,h}(Ri_g)\left|\frac{\partial U}{\partial z}\right|
  !!  \f]
  !!  The functions used (\f$f_{m,h}\f$) depend on the local stability. First, the gradient Richardson number is calculated as
  !!  \f[
  !!  Ri_g=\frac{\frac{g}{T}\frac{\partial \theta_v}{\partial z}}{\frac{\partial U}{\partial z}^2}
  !!  \f]
  !!  where \f$U\f$ is the horizontal wind. For the unstable case (\f$Ri_g < 0\f$), the Richardson number-dependent functions are given by
  !!  \f[
  !!  f_h(Ri_g) = 1 + \frac{8\left|Ri_g\right|}{1 + 1.286\sqrt{\left|Ri_g\right|}}\\
  !!  \f]
  !!  \f[
  !!  f_m(Ri_g) = 1 + \frac{8\left|Ri_g\right|}{1 + 1.746\sqrt{\left|Ri_g\right|}}\\
  !!  \f]
  !!  For the stable case, the following formulas are used
  !!  \f[
  !!  f_h(Ri_g) = \frac{1}{\left(1 + 5Ri_g\right)^2}\\
  !!  \f]
  !!  \f[
  !!  Pr = \frac{K_h}{K_m} = 1 + 2.1Ri_g
  !!  \f]
  !!  The source for the formulas used for the Richardson number-dependent functions is unclear. They are different than those used in Hong and Pan (1996) \cite hong_and_pan_1996 as the previous documentation suggests. They follow equation 14 of Louis (1979
  !!  \f[
  !!  \frac{1}{l} = \frac{1}{kz} + \frac{1}{l_0}\\
  !!  \f]
  !!  \f[
  !!  or\\
  !!  \f]
  !!  \f[
  !!  l=\frac{l_0kz}{l_0+kz}
  !!  \f]
  !!  where \f$l_0\f$ is currently 30 m for stable conditions and 150 m for unstable. Finally, the diffusion coefficients are kept in a range bounded by the background diffusion and the maximum allowable values.
         !$acc parallel loop gang collapse(2) &
         !$acc&         private(jj,k,i,bvf2,ri,zk,rl2,dk,sri,tem1,prnum, &
         !$acc&                 dkur,dktr,shr2r) async(async_id)
         do jj = 1, jlistnum
            do k = 1, km1
               !$acc loop vector
               do i = 1, myim(jj)
                  if (k >= kpbl(i, jj)) then
                     shr2r = shr2(i, k, jj)
                     bvf2 = g*bf(i, k, jj)*ti(i, k, jj)
                     ri = max(bvf2/shr2r, rimin)
                     zk = vk*zi(i, k + 1, jj)
                     if (ri < 0.) then ! unstable regime
                        rl2 = zk*rlamun/(rlamun + zk)
                        dk = rl2*rl2*sqrt(shr2r)
                        sri = sqrt(-ri)
!               dku(i,k) = xkzmo(i,k) + dk*(1+8.*(-ri)/(1+1.746*sri))
!               dkt(i,k) = xkzo(i,k)  + dk*(1+8.*(-ri)/(1+1.286*sri))
                        dkur = dk*(1 + 8.*(-ri)/(1 + 1.746*sri))
                        dktr = dk*(1 + 8.*(-ri)/(1 + 1.286*sri))
                     else             ! stable regime
                        rl2 = zk*rlam/(rlam + zk)
!!              tem      = rlam * sqrt(0.01*prsi(i,k))
!!              rl2      = zk*tem/(tem+zk)
                        dk = rl2*rl2*sqrt(shr2r)
                        tem1 = dk/(1 + 5.*ri)**2
!
                        if (k >= kpblx(i, jj)) then
                           prnum = 1.0 + 2.1*ri
                           prnum = min(prnum, prmax)
                        else
                           prnum = 1.0
                        end if
!               dku(i,k) = xkzmo(i,k) + tem1 * prnum
!               dkt(i,k) = xkzo(i,k)  + tem1
                        dkur = tem1*prnum
                        dktr = tem1
                     end if
!
                     dkur = min(dkur, dkmax)
                     dkur = max(dkur, xkzmo(i, k, jj))
                     dktr = min(dktr, dkmax)
                     dktr = max(dktr, xkzo(i, k, jj))
                     dku(i, k, jj) = dkur
                     dkt(i, k, jj) = dktr
!
                  end if
!
               end do
            end do
         end do

         !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !  compute components for mass flux mixing by large thermals
         !>  ## If the PBL is convective, call the mass flux scheme to replace the countergradient terms.
  !!  If the PBL is convective, the updraft properties are initialized to be the same as the state variables and the subroutine mfpbl is called.
         !$acc parallel loop gang vector collapse(3) private(jj,k,i) async(async_id)
         do jj = 1, jlistnum
            do k = 1, km
               do i = 1, ix
                  if (i .le. myim(jj)) then
                  if (pcnvflg(i, jj)) then
                     tcko(i, k, jj) = t1(i, k, jj)
                     ucko(i, k, jj) = u1(i, k, jj)
                     vcko(i, k, jj) = v1(i, k, jj)
                     xmf(i, k, jj) = 0.
                     !$acc loop seq
                     do kk = 1, ntrac
                        qcko(i, k, kk, jj) = q1(i, k, kk, jj)
                     end do
                  end if
                  end if
               end do
            end do
         end do
         !>  For details of the mfpbl subroutine, step into its documentation ::mfpbl
         call mfpbl_gpu(myim, ix, km, ntrac, dt2, pcnvflg, &
                        zl, zi, thvx, q1, t1, u1, v1, hpbl, kpbl, &
                        sflux, ustar, wstar, xmf, tcko, qcko, ucko, vcko, async_id)

         !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !  compute diffusion coefficients for cloud-top driven diffusion
         !  if the condition for cloud-top instability is met,
         !    increase entrainment flux at cloud top
         !
         !>  ## Compute enhanced diffusion coefficients related to stratocumulus-topped PBLs
  !!  If a stratocumulus layer has been identified in the PBL, the diffusion coefficients in the PBL are modified in the following way.
  !!
  !!  -# First, the criteria for CTEI is checked, using the threshold from equation 13 of Macvean and Mason (1990) \cite macvean_and_mason_1990. If the criteria is met, the cloud top diffusion is increased:
  !!  \f[
  !!  K_h^{Sc} = -c\frac{\Delta F_R}{\rho c_p}\frac{1}{\frac{\partial \theta_v}{\partial z}}
  !!  \f]
  !!  where the constant \f$c\f$ is set to 0.2 if the CTEI criterion is not met and 1.0 if it is.
  !!
  !!  -# Calculate the diffusion coefficients due to stratocumulus mixing according to equation 5 in Lock et al. (2000) \cite lock_et_al_2000 for every level below the stratocumulus top using the characteristic stratocumulus velocity scale previously calcul
         !$acc parallel loop gang vector collapse(2) &
         !$acc&              private(jj,k,i,tem1,tem2,ptem,dkur,dktr) async(async_id)
         do jj = 1, jlistnum
            do i = 1, ix
               if (i .le. myim(jj) .and. scuflg(i, jj)) then
                  m = krad(i, jj)
                  tem = thetae(i, m, jj) - thetae(i, m + 1, jj)
                  tem1 = qtx(i, m, jj) - qtx(i, m + 1, jj)
                  if (tem > 0. .and. tem1 > 0.) then
                     cteit = cp*tem/(hvap*tem1)
                     if (cteit > actei) rent(i, jj) = rentf2
                  end if
                  tem1 = max(bf(i, m, jj), tdzmin)
                  ckt(i, m, jj) = -rent(i, jj)*radmin(i, jj)/tem1
                  cku(i, m, jj) = ckt(i, m, jj)

                  !$acc loop seq
                  do k = 1, kmpbl
                     if (k < krad(i, jj)) then
                        !if (1) then    !maybe have some error
                        tem1 = hrad(i, jj) - zd(i, jj)
                        tem2 = zi(i, k + 1, jj) - tem1
                        if (tem2 > 0.) then
                           ptem = tem2/zd(i, jj)
                           if (ptem .ge. 1.) ptem = 1.
                           ptem = tem2*ptem*sqrt(1.-ptem)
                           ckt(i, k, jj) = radfac*vk*vrad(i, jj)*ptem
                           cku(i, k, jj) = 0.75*ckt(i, k, jj)
                           ckt(i, k, jj) = max(ckt(i, k, jj), dkmin)
                           ckt(i, k, jj) = min(ckt(i, k, jj), dkmax)
                           cku(i, k, jj) = max(cku(i, k, jj), dkmin)
                           cku(i, k, jj) = min(cku(i, k, jj), dkmax)
                        end if
                     end if
                     !>  After \f$K_h^{Sc}\f$ has been determined from the surface to the top of the stratocumulus layer, it is added to the value for the diffusion coefficient calculated previously using surface-based mixing [see equation 6 of Lock et al. (2000) \cite lock_e
                     dktr = dkt(i, k, jj) + ckt(i, k, jj)
                     dkur = dku(i, k, jj) + cku(i, k, jj)
                     dktr = min(dktr, dkmax)
                     dkur = min(dkur, dkmax)
                     dkt(i, k, jj) = dktr
                     dku(i, k, jj) = dkur
                  end do
               end if
            end do
         end do
         !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         !

         if (.false.) then
            !$acc parallel loop gang vector collapse(2) private(dtodsd,&
            !$acc&                   dtodsu,dsig,rdz,aur,alr, &
            !$acc&                   tem1,dsdz2,tem2,ptem,ptem1,ptem2,&
            !$acc&                   dsdzt,dsdzq,tem,jj,i,k,accui) async(async_id)
            do jj = 1, jlistnum
               do i = 1, ix
                  if (i .le. myim(jj)) then
                     accui = nxjp_acc(jj) - 1
                     alg(accui + i, 1) = 0.
                     aug(accui + i, km) = 0.
                     adg(accui + i, 1) = 1.
                     a1g(accui + i, 1) = t1(i, 1, jj) + beta(i, jj)*heat(i, jj)
                     a2g(accui + i, 1) = q1(i, 1, 1, jj) + beta(i, jj)*evap(i, jj)
                     !$acc loop seq
                     do k = 1, km1
                        dtodsd = dt2/del(i, k, jj)
                        dtodsu = dt2/del(i, k + 1, jj)
                        dsig = prsl(i, k, jj) - prsl(i, k + 1, jj)
                        rdz = rdzt(i, k, jj)
                        tem1 = dsig*dkt(i, k, jj)*rdz
                        dsdz2 = tem1*rdz
                        aur = -dtodsd*dsdz2
                        alr = -dtodsu*dsdz2
                        !
                        if (pcnvflg(i, jj) .and. k < kpbl(i, jj)) then
                           tem2 = dsig*rdz
                           ptem = 0.5*tem2*xmf(i, k, jj)
                           ptem1 = dtodsd*ptem
                           ptem2 = dtodsu*ptem
                           adg(accui + i, k) = adg(accui + i, k) - aur - ptem1
                           adg(accui + i, k + 1) = 1.-alr + ptem2
                           aur = aur - ptem1
                           alr = alr + ptem2
                           ptem = tcko(i, k, jj) + tcko(i, k + 1, jj)
                           dsdzt = tem1*gocp
                           a1g(accui + i, k) = a1g(accui + i, k) + dtodsd*dsdzt - ptem1*ptem
                           a1g(accui + i, k + 1) = t1(i, k + 1, jj) - dtodsu*dsdzt + ptem2*ptem
                           ptem = qcko(i, k, 1, jj) + qcko(i, k + 1, 1, jj)
                           a2g(accui + i, k) = a2g(accui + i, k) - ptem1*ptem
                           a2g(accui + i, k + 1) = q1(i, k + 1, 1, jj) + ptem2*ptem
                        elseif (ublflg(i, jj) .and. k < kpbl(i, jj)) then
                           ptem1 = dsig*dktx(i, k, jj)*rdz
                           tem = 1.0/hpbl(i, jj)
                           dsdzt = tem1*gocp - ptem1*hgamt(i, jj)*tem
                           dsdzq = -ptem1*hgamq(i, jj)*tem
                           adg(accui + i, k) = adg(accui + i, k) - aur
                           adg(accui + i, k + 1) = 1.-alr
                           a1g(accui + i, k) = a1g(accui + i, k) + dtodsd*dsdzt
                           a1g(accui + i, k + 1) = t1(i, k + 1, jj) - dtodsu*dsdzt
                           a2g(accui + i, k) = a2g(accui + i, k) + dtodsd*dsdzq
                           a2g(accui + i, k + 1) = q1(i, k + 1, 1, jj) - dtodsu*dsdzq
                        else
                           adg(accui + i, k) = adg(accui + i, k) - aur
                           adg(accui + i, k + 1) = 1.-alr
                           dsdzt = tem1*gocp
                           a1g(accui + i, k) = a1g(accui + i, k) + dtodsd*dsdzt
                           a1g(accui + i, k + 1) = t1(i, k + 1, jj) - dtodsu*dsdzt
                           a2g(accui + i, k + 1) = q1(i, k + 1, 1, jj)
                        end if
                        aug(accui + i, k) = aur
                        aug_backup(accui + i, k) = aur
                        alg(accui + i, k + 1) = alr
                     end do
                  end if
               end do
            end do
            call cusparseDgtsvInterleavedBatch_async(handle, 0, km, alg, adg, aug, a1g, nxptot, async_id)
            !$acc parallel loop gang collapse(2) private(jj,i,k,accui) async(async_id)
            do jj = 1, jlistnum
               do k = 1, km
                  !$acc loop vector
                  do i = 1, myim(jj)
                     accui = nxjp_acc(jj) - 1
                     if (k .eq. km) then
                        aug(accui + i, k) = 0.
                     else
                        aug(accui + i, k) = aug_backup(accui + i, k)
                     end if
                  end do
               end do
            end do
            call cusparseDgtsvInterleavedBatch_async(handle, 0, km, alg, adg, aug, a2g, nxptot, async_id)
            !$acc parallel loop gang collapse(2) private(jj,k,i,accui,qtend,ttend) async(async_id)
            do jj = 1, jlistnum
               do k = 1, km
                  !$acc loop vector
                  do i = 1, myim(jj)
                     accui = nxjp_acc(jj) - 1
                     !
                     !     recover tendencies of heat and moisture
                     !
                     !>  After returning with the solution, the tendencies for temperature and moisture are recovered.
                     qtend = (a2g(accui + i, k) - q1(i, k, 1, jj))*rdt
                     ttend = (a1g(accui + i, k) - t1(i, k, jj))*rdt
                     tau(i, k, jj) = ttend
  !!            tau(i,k)   = tau(i,k)+ttend
  !!            rtg(i,k,1) = rtg(i,k,1)+qtend

                     !dtsfc(i,jj) = dtsfc(i,jj) + cont*del(i,k,jj)*ttend
                     !dqsfc(i,jj) = dqsfc(i,jj) + conq*del(i,k,jj)*qtend
                     !rtg(i, k, 1, jj) = qtend
                     q1(i, k, 1, jj) = a2g(accui + i, k)
                  end do
               end do
            end do
            !
            if (ntrac >= 2) then
               do kk = 2, ntrac
                  !$acc parallel loop gang vector collapse(2) private(jj,i,k, &
                  !$acc&                   dtodsu,dsig,tem,ptem2,&
                  !$acc&                   tem1,dtodsd,ptem1) async(async_id)
                  do jj = 1, jlistnum
                     do i = 1, ix
                        if (i .le. myim(jj)) then
                           accui = nxjp_acc(jj) - 1
                           a2g(accui + i, 1) = q1(i, 1, kk, jj)
                           !$acc loop seq
                           do k = 1, km1
                              if (pcnvflg(i, jj) .and. k < kpbl(i, jj)) then
                                 dtodsu = dt2/del(i, k + 1, jj)
                                 dtodsd = dt2/del(i, k, jj)
                                 dsig = prsl(i, k, jj) - prsl(i, k + 1, jj)
                                 tem = dsig*rdzt(i, k, jj)
                                 ptem = 0.5*tem*xmf(i, k, jj)
                                 ptem1 = dtodsd*ptem
                                 ptem2 = dtodsu*ptem
                                 tem1 = qcko(i, k, kk, jj) + qcko(i, k + 1, kk, jj)
                                 a2g(accui + i, k) = a2g(accui + i, k) - ptem1*tem1
                                 a2g(accui + i, k + 1) = q1(i, k + 1, kk, jj) + ptem2*tem1
                              else
                                 a2g(accui + i, k + 1) = q1(i, k + 1, kk, jj)
                              end if
                              aug(accui + i, k) = aug_backup(accui + i, k)
                           end do
                           aug(accui + i, km) = 0.
                        end if
                     end do
                  end do

                  !
                  !     solve tridiagonal problem for heat and moisture
                  !
                  !>  The tridiagonal system is solved by calling the internal ::tridin subroutine.
                  call cusparseDgtsvInterleavedBatch_async(handle, 0, km, alg, adg, aug, a2g, nxptot, async_id)
                  !$acc parallel loop gang collapse(2) private(jj,k,i,accui,qtend,ttend) async(async_id)
                  do jj = 1, jlistnum
                     do k = 1, km
                        !$acc loop vector
                        do i = 1, myim(jj)
                           accui = nxjp_acc(jj) - 1
                           !
                           !     recover tendencies of heat and moisture
                           !
                           !>  After returning with the solution, the tendencies for temperature and moisture are recovered.
                           qtend = (a2g(accui + i, k) - q1(i, k, kk, jj))*rdt
                           !rtg(i, k, kk, jj) = qtend
                           q1(i, k, kk, jj) = a2g(accui + i, k)
                        end do
                     end do
                  end do
               end do
            end if
         else
         !     compute tridiagonal matrix elements for heat and moisture
         !
         !>  ## Solve for the temperature and moisture tendencies due to vertical mixing.
  !!  The tendencies of heat, moisture, and momentum due to vertical diffusion are calculated using a two-part process. First, a solution is obtained using an implicit time-stepping scheme, then the time tendency terms are "backed out". The tridiagonal matr
            !$acc parallel loop gang vector collapse(2) private(dtodsd,&
            !$acc&                   dtodsu,dsig,rdz,aur,alr,a3gr2,a3gr1, &
            !$acc&                   tem1,dsdz2,tem2,ptem,ptem1,ptem2,&
            !$acc&                   dsdzt,dsdzq,tem,jj,i,k,accui,adgr) async(async_id)
            do jj = 1, jlistnum
               do i = 1, ix
                  if (i .le. myim(jj)) then
                     ii = nxjp_acc(jj) - 1 + i
                     alg(ii, 1) = 0.
                     aug(ii, km) = 0.
                     adgr = 1.
                     a3gr1 = t1(i, 1, jj) + beta(i, jj)*heat(i, jj)
                     a3gr2 = q1(i, 1, 1, jj) + beta(i, jj)*evap(i, jj)
                     !$acc loop seq
                     do k = 1, km1
                        dtodsd = dt2/del(i, k, jj)
                        dtodsu = dt2/del(i, k + 1, jj)
                        dsig = prsl(i, k, jj) - prsl(i, k + 1, jj)
                        rdz = rdzt(i, k, jj)
                        tem1 = dsig*dkt(i, k, jj)*rdz
                        dsdz2 = tem1*rdz
                        aur = -dtodsd*dsdz2
                        alr = -dtodsu*dsdz2
                        !
                        if (pcnvflg(i, jj) .and. k < kpbl(i, jj)) then
                           tem2 = dsig*rdz
                           ptem = 0.5*tem2*xmf(i, k, jj)
                           ptem1 = dtodsd*ptem
                           ptem2 = dtodsu*ptem
                           adg(ii, k) = adgr - aur - ptem1
                           adgr = 1.-alr + ptem2
                           aur = aur - ptem1
                           alr = alr + ptem2
                           ptem = tcko(i, k, jj) + tcko(i, k + 1, jj)
                           dsdzt = tem1*gocp
                           a3g(ii, k, 1) = a3gr1 + dtodsd*dsdzt - ptem1*ptem
                           a3gr1 = t1(i, k + 1, jj) - dtodsu*dsdzt + ptem2*ptem
                           ptem = qcko(i, k, 1, jj) + qcko(i, k + 1, 1, jj)
                           a3g(ii, k, 2) = a3gr2 - ptem1*ptem
                           a3gr2 = q1(i, k + 1, 1, jj) + ptem2*ptem
                        elseif (ublflg(i, jj) .and. k < kpbl(i, jj)) then
                           ptem1 = dsig*dktx(i, k, jj)*rdz
                           tem = 1.0/hpbl(i, jj)
                           dsdzt = tem1*gocp - ptem1*hgamt(i, jj)*tem
                           dsdzq = -ptem1*hgamq(i, jj)*tem
                           adg(ii, k) = adgr - aur
                           adgr = 1.-alr
                           a3g(ii, k, 1) = a3gr1 + dtodsd*dsdzt
                           a3gr1 = t1(i, k + 1, jj) - dtodsu*dsdzt
                           a3g(ii, k, 2) = a3gr2 + dtodsd*dsdzq
                           a3gr2 = q1(i, k + 1, 1, jj) - dtodsu*dsdzq
                        else
                           adg(ii, k) = adgr - aur
                           adgr = 1.-alr
                           dsdzt = tem1*gocp
                           a3g(ii, k, 1) = a3gr1 + dtodsd*dsdzt
                           a3gr1 = t1(i, k + 1, jj) - dtodsu*dsdzt
                           a3g(ii, k, 2) = a3gr2
                           a3gr2 = q1(i, k + 1, 1, jj)
                        end if
                        aug(ii, k) = aur
                        alg(ii, k + 1) = alr
                     end do
                     adg(ii, km) = adgr
                     a3g(ii, km, 1) = a3gr1
                     a3g(ii, km, 2) = a3gr2
                  end if
               end do
            end do
            !
            !$acc parallel loop gang vector collapse(2) private(jj,i,k, &
            !$acc&                   dtodsu,dsig,tem2,ptem2,ptem,&
            !$acc&                   tem1,dtodsd,ptem1,tem3) async(async_id)
            do jj = 1, jlistnum
               do i = 1, ix
                  if (i .le. myim(jj)) then
                     ii = nxjp_acc(jj) - 1 + i
                        !$acc loop seq
                        do kk = 2, ntrac
                           a3g(ii, 1, kk+1) = q1(i, 1, kk, jj)
                        end do
                        !$acc loop seq
                        do k = 1, km1
                           if (pcnvflg(i, jj) .and. k < kpbl(i, jj)) then
                              rdz = rdzt(i, k, jj)
                              dtodsu = dt2/del(i, k + 1, jj)
                              dtodsd = dt2/del(i, k, jj)
                              dsig = prsl(i, k, jj) - prsl(i, k + 1, jj)
                              tem2 = dsig*rdz
                              ptem = 0.5*tem2*xmf(i, k, jj)
                              ptem1 = dtodsd*ptem
                              ptem2 = dtodsu*ptem
                              !$acc loop seq
                              do kk = 2, ntrac
                                 tem3 = qcko(i, k, kk, jj) + qcko(i, k + 1, kk, jj)
                                 a3g(ii, k, kk+1) = a3g(ii, k, kk+1) - ptem1*tem3
                                 a3g(ii, k + 1, kk+1) = q1(i, k + 1, kk, jj) + ptem2*tem3
                              end do
                           else
                              !$acc loop seq
                              do kk = 2, ntrac
                                 a3g(ii, k + 1, kk+1) = q1(i, k + 1, kk, jj)
                              end do
                           end if
                        end do
                     end if
                  end do
               end do

            !
            !     solve tridiagonal problem for heat and moisture
            !
            !>  The tridiagonal system is solved by calling the internal ::tridin subroutine.
            call tridin_gpu(nxptot, km, ntrac, alg, adg, aug, a3g, aug, a3g, async_id)
            !$acc parallel loop gang vector collapse(3) private(m,jj,k,i,accui,qtend,ttend) async(async_id)
            do jj = 1, jlistnum
               do k = 1, km
                  do i = 1, ix
                     if (i .le. myim(jj)) then
                        accui = nxjp_acc(jj) - 1
                        !
                        !     recover tendencies of heat and moisture
                        !
                        !>  After returning with the solution, the tendencies for temperature and moisture are recovered.
                        ttend = (a3g(accui + i, k, 1) - t1(i, k, jj))*rdt
                        tau(i, k, jj) = ttend
     !!            tau(i,k)   = tau(i,k)+ttend
     !!            rtg(i,k,1) = rtg(i,k,1)+qtend

                        !dtsfc(i,jj) = dtsfc(i,jj) + cont*del(i,k,jj)*ttend
                        !dqsfc(i,jj) = dqsfc(i,jj) + conq*del(i,k,jj)*qtend
                        
                        !$acc loop seq
                        do kk = 1, ntrac
                           m = kk + 1
                           !qtend = (a3g(accui + i, k, m) - q1(i, k, kk, jj))*rdt
                           !rtg(i, k, kk, jj) = qtend
                           q1(i, k, kk, jj) = a3g(accui + i, k, m)
                        end do
                     end if
                  end do
               end do
            end do
         end if

         !
         !   compute tke dissipation rate
         !
         !>  ## Calculate heating due to TKE dissipation and add to the tendency for temperature
  !!  Following Han et al. (2015) \cite han_et_al_2015 , turbulence dissipation contributes to the tendency of temperature in the following way. First, turbulence dissipation is calculated by equation 17 of Han et al. (2015) \cite han_et_al_2015 for the PBL
         if (dspheat) then
            !
            !$acc parallel loop gang collapse(2) private(jj,k,i) async(async_id)
            do jj = 1, jlistnum
               do k = 1, km1
                  !$acc loop vector
                  do i = 1, myim(jj)
                     diss(i, k, jj) = dku(i, k, jj)*shr2(i, k, jj) - g*ti(i, k, jj)*dkt(i, k, jj)*bf(i, k, jj)
                     !           diss(i,k) = dku(i,k)*shr2(i,k)
                  end do
               end do
            end do
         end if
         !
         !
         !
         !$acc parallel loop gang collapse(2) private(jj,k,i,tem,ttend,tem1,tem2) async(async_id)
         do jj = 1, jlistnum
            do k = 1, km
               !$acc loop vector
               do i = 1, myim(jj)
                  if (k .ne. km) then
                     if (dspheat) then
                        if (k .eq. 1) then
         !     add dissipative heating at the first model layer
         !
         !>  Next, the temperature tendency is updated following equation 14.
                           tem = govrth(i, jj)*sflux(i, jj)
                           tem1 = tem + stress(i, jj)*spd1(i, jj)/zl(i, 1, jj)
                           tem2 = 0.5*(tem1 + diss(i, 1, jj))
                        else
                           tem2 = 0.5*(diss(i, k - 1, jj) + diss(i, k, jj))
                        end if
         !     add dissipative heating above the first model layer
                        tem2 = max(tem2, 0.)
                        ttend = tem2/cp
                        tau(i, k, jj) = tau(i, k, jj) + 0.5*ttend
                     end if
                  end if
                  !    update temperature tendency to model layer mean temperature
                  t1(i, k, jj) = t1(i, k, jj) + tau(i, k, jj)*dt2
               end do
            end do
         end do
         !
         !
         !     compute tridiagonal matrix elements for momentum
         !

         !
         !
         if (.false.) then
            !may cause some absolute error
            !$acc parallel loop gang vector collapse(2) private(jj,i,k,&
            !$acc&                   dtodsd,dtodsu,dsig,rdz,aur,alr,adgr,&
            !$acc&                   tem1,dsdz2,tem2,ptem,ptem1,ptem2,accui) async(async_id)
            do jj = 1, jlistnum
               do i = 1, ix
                  if (i .le. myim(jj)) then
                     accui = nxjp_acc(jj) - 1
                     alg(accui + i, 1) = 0.
                     aug(accui + i, km) = 0.
                     adg(accui + i, 1) = 1.0 + beta(i, jj)*stress(i, jj)/spd1(i, jj)
                     a1g(accui + i, 1) = u1(i, 1, jj)
                     a2g(accui + i, 1) = v1(i, 1, jj)
                     !$acc loop seq
                     do k = 1, km1
                        dtodsd = dt2/del(i, k, jj)
                        dtodsu = dt2/del(i, k + 1, jj)
                        dsig = prsl(i, k, jj) - prsl(i, k + 1, jj)
                        rdz = rdzt(i, k, jj)
                        tem1 = dsig*dku(i, k, jj)*rdz
                        dsdz2 = tem1*rdz
                        aur = -dtodsd*dsdz2
                        alr = -dtodsu*dsdz2
                        !
                        if (pcnvflg(i, jj) .and. k < kpbl(i, jj)) then
                           tem2 = dsig*rdz
                           ptem = 0.5*tem2*xmf(i, k, jj)
                           ptem1 = dtodsd*ptem
                           ptem2 = dtodsu*ptem
                           adg(accui + i, k) = adg(accui + i, k) - aur - ptem1
                           adg(accui + i, k + 1) = 1.-alr + ptem2
                           aur = aur - ptem1
                           alr = alr + ptem2
                           ptem = ucko(i, k, jj) + ucko(i, k + 1, jj)
                           a1g(accui + i, k) = a1g(accui + i, k) - ptem1*ptem
                           a1g(accui + i, k + 1) = u1(i, k + 1, jj) + ptem2*ptem
                           ptem = vcko(i, k, jj) + vcko(i, k + 1, jj)
                           a2g(accui + i, k) = a2g(accui + i, k) - ptem1*ptem
                           a2g(accui + i, k + 1) = v1(i, k + 1, jj) + ptem2*ptem
                        else
                           adg(accui + i, k) = adg(accui + i, k) - aur
                           adg(accui + i, k + 1) = 1.-alr
                           a1g(accui + i, k + 1) = u1(i, k + 1, jj)
                           a2g(accui + i, k + 1) = v1(i, k + 1, jj)
                        end if
                        aug(accui + i, k) = aur
                        aug_backup(accui + i, k) = aur
                        alg(accui + i, k + 1) = alr
                     end do
                  end if
               end do
            end do

            call cusparseDgtsvInterleavedBatch_async(handle, 0, km, alg, adg, aug, a1g, nxptot, async_id)
            !$acc parallel loop gang collapse(2) private(jj,i,k,accui) async(async_id)
            do jj = 1, jlistnum
               do k = 1, km
                  !$acc loop vector
                  do i = 1, myim(jj)
                     accui = nxjp_acc(jj) - 1
                     if (k .eq. km) then
                        aug(accui + i, k) = 0.
                     else
                        aug(accui + i, k) = aug_backup(accui + i, k)
                     end if
                  end do
               end do
            end do
            call cusparseDgtsvInterleavedBatch_async(handle, 0, km, alg, adg, aug, a2g, nxptot, async_id)
         else
         !>  ## Solve for the horizontal momentum tendencies and add them to the output tendency terms
  !!  As with the temperature and moisture tendencies, the horizontal momentum tendencies are calculated by solving tridiagonal matrices after the matrices are prepared in this section.
            !$acc parallel loop gang vector collapse(2) private(jj,i,k,&
            !$acc&                   dtodsd,dtodsu,dsig,rdz,aur,alr,adgr,a3gr1,a3gr2,&
            !$acc&                   tem1,dsdz2,tem2,ptem,ptem1,ptem2,accui) async(async_id)
            do jj = 1, jlistnum
               do i = 1, ix
                  if (i .le. myim(jj)) then
                     accui = nxjp_acc(jj) - 1
                     alg(accui + i, 1) = 0.
                     aug(accui + i, km) = 0.
                     adgr = 1.0 + beta(i, jj)*stress(i, jj)/spd1(i, jj)
                     a3gr1 = u1(i, 1, jj)
                     a3gr2 = v1(i, 1, jj)
                     !$acc loop seq
                     do k = 1, km1
                        dtodsd = dt2/del(i, k, jj)
                        dtodsu = dt2/del(i, k + 1, jj)
                        dsig = prsl(i, k, jj) - prsl(i, k + 1, jj)
                        rdz = rdzt(i, k, jj)
                        tem1 = dsig*dku(i, k, jj)*rdz
                        dsdz2 = tem1*rdz
                        aur = -dtodsd*dsdz2
                        alr = -dtodsu*dsdz2
                        !
                        if (pcnvflg(i, jj) .and. k < kpbl(i, jj)) then
                           tem2 = dsig*rdz
                           ptem = 0.5*tem2*xmf(i, k, jj)
                           ptem1 = dtodsd*ptem
                           ptem2 = dtodsu*ptem
                           adg(accui + i, k) = adgr - aur - ptem1
                           adgr = 1.-alr + ptem2
                           aur = aur - ptem1
                           alr = alr + ptem2
                           ptem = ucko(i, k, jj) + ucko(i, k + 1, jj)
                           a3g(accui + i, k, 1) = a3gr1 - ptem1*ptem
                           a3gr1 = u1(i, k + 1, jj) + ptem2*ptem
                           ptem = vcko(i, k, jj) + vcko(i, k + 1, jj)
                           a3g(accui + i, k, 2) = a3gr2 - ptem1*ptem
                           a3gr2 = v1(i, k + 1, jj) + ptem2*ptem
                        else
                           adg(accui + i, k) = adgr - aur
                           adgr = 1.-alr
                           a3g(accui + i, k, 1) = a3gr1
                           a3gr1 = u1(i, k + 1, jj)
                           a3g(accui + i, k, 2) = a3gr2
                           a3gr2 = v1(i, k + 1, jj)
                        end if
                        aug(accui + i, k) = aur
                        alg(accui + i, k + 1) = alr
                     end do
                        adg(accui + i, km) = adgr
                        a3g(accui + i, km, 1) = a3gr1
                        a3g(accui + i, km, 2) = a3gr2
                  end if
               end do
            end do
         !     solve tridiagonal problem for momentum

            call tridi2_gpu(nxptot, km, alg, adg, aug, a3g, aug, a3g, async_id)
         end if
!     recover tendencies of momentum
!
!>  Finally, the tendencies are recovered from the tridiagonal solutions.

         !$acc parallel loop gang collapse(2) private(jj,k,i,accui) async(async_id)
         do jj = 1, jlistnum
            do k = 1, km
               !$acc loop vector
               do i = 1, myim(jj)
                  accui = nxjp_acc(jj) - 1
                  !a1(i,k,jj) = a1g(accui+i,k)
                  !a2(i,k,jj) = a2g(accui+i,k)
                  !     recover tendencies of momentum
                  !
                  !>  Finally, the tendencies are recovered from the tridiagonal solutions.
                  !utend = (a1(i,k,jj)-u1(i,k,jj))*rdt
                  !vtend = (a2(i,k,jj)-v1(i,k,jj))*rdt
  !!            du(i,k)  = du(i,k)  + utend
  !!            dv(i,k)  = dv(i,k)  + vtend
                  !du(i,k,jj)  = utend
                  !dv(i,k,jj)  = vtend
                  !dusfc(i,jj) = dusfc(i,jj) + conw*del(i,k,jj)*utend
                  !dvsfc(i,jj) = dvsfc(i,jj) + conw*del(i,k,jj)*vtend
                  u1(i, k, jj) = a3g(accui + i, k, 1)   ! cusparse a1g -> a1 -> u1
                  v1(i, k, jj) = a3g(accui + i, k, 2)   ! cusparse a2g -> a2 -> v1
                  !
                  !  for dissipative heating for ecmwf model
                  !
                  !           tem1 = 0.5*(a1(i,k)+u1(i,k))
                  !           tem2 = 0.5*(a2(i,k)+v1(i,k))
                  !           diss(i,k) = -(tem1*utend+tem2*vtend)
                  !           diss(i,k) = max(diss(i,k),0.)
                  !           ttend = diss(i,k) / cp
                  !           tau(i,k) = tau(i,k) + ttend
                  !
               end do
            end do
         end do

         !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         !$acc parallel loop gang vector collapse(2) private(jj,i) async(async_id)
         do jj = 1, jlistnum
            do i = 1, ix
               if (i .le. myim(jj)) then
                  hpbl(i, jj) = hpblx(i, jj)
                  kpbl(i, jj) = kpblx(i, jj)
               end if
            end do
         end do
         !$acc exit data delete(beta,bf,ckt,cku,diss, &
         !$acc&                 dusfc,dvsfc,dtsfc,dqsfc,dku,dkt,dktx, &
         !$acc&                 crb,govrth,hgamt,hgamq, &
         !$acc&                 hpblx,hrad,icld,kcld,krad,kpblx, &
         !$acc&                 kx1,lcld,pblflg,pcnvflg,phih,phim, &
         !$acc&                 prinv,qcko,qlx,qtx,radmin,radx,rent,rbdn, &
         !$acc&                 rbup,rdzt,sfcflg,scuflg,shr2,sflux, &
         !$acc&                 tau,tcko,theta,thetae,thermal, &
         !$acc&                 thlvx1,thlvx,thvx,ti,ublflg,ucko, &
         !$acc&                 ustar,ustmin,vcko,vrad,wscale,wscaleu,wstar, &
         !$acc&                 xmf,z0,zd,zdd,zi,zl,zol,xkzo,xkzmo,tx1,tx2, &
         !$acc&                 alg,adg,aug,a3g) async(async_id)
         !!$acc exit data delete(rtg) async(async_id)
         !!$acc exit data delete(a1g, a2g, aug_backup) async(async_id) !cusparse
         !!$acc wait(async_id)

         return
      end subroutine moninedmf_gpu

