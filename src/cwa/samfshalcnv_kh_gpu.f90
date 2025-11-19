!> \defgroup SAMF_shal Scale-Aware Mass-Flux Shallow Convection
!! @{
!!  \brief The scale-aware mass-flux shallow (SAMF_shal) convection scheme is an updated version of the previous mass-flux shallow convection scheme with scale and aerosol awareness and parameterizes the effect of shallow convection on the environment.  T
!!
!!  The previous version of the shallow convection scheme (shalcnv.f) is described in Han and Pan (2011) \cite han_and_pan_2011 and differences between the shallow and deep convection schemes are presented in Han and Pan (2011) \cite han_and_pan_2011 and
!!
!!  \section diagram Calling Hierarchy Diagram
!!  \image html SAMF_shal_Flowchart.png "Diagram depicting how the SAMF shallow convection scheme is called from the FV3GFS physics time loop" height=2cm
!!  \section intraphysics Intraphysics Communication
!!  This space is reserved for a description of how this scheme uses information from other scheme types and/or how information calculated in this scheme is used in other scheme types.

!> \file samfshalcnv.f
!!  Contains the entire SAMF shallow convection scheme.

!>  \brief This subroutine contains the entirety of the SAMF shallow convection scheme.
!!
!!  This routine follows the \ref SAMF deep scheme quite closely, although it can be interpreted as only having the "static" and "feedback" control portions, since the "dynamic" control is not necessary to find the cloud base mass flux. The algorithm is s
!!
!!  \param[in] im number of used points
!!  \param[in] ix horizontal dimension
!!  \param[in] km vertical layer dimension
!!  \param[in] jcap number of spectral wave trancation
!!  \param[in] delt physics time step in seconds
!!  \param[in] delp pressure difference between level k and k+1 (Pa)
!!  \param[in] prslp mean layer presure (Pa)
!!  \param[in] psp surface pressure (Pa)
!!  \param[in] phil layer geopotential (\f$m^s/s^2\f$)
!!  \param[inout] ql cloud water (kg/kg)
!!  \param[inout] qi cloud ice (kg/kg)
!!  \param[inout] q1 updated tracers (kg/kg)
!!  \param[inout] t1 updated temperature (K)
!!  \param[inout] u1 updated zonal wind (\f$m s^{-1}\f$)
!!  \param[inout] v1 updated meridional wind (\f$m s^{-1}\f$)
!!  \param[out] rn convective rain (m)
!!  \param[out] kbot index for cloud base
!!  \param[out] ktop index for cloud top
!!  \param[out] kcnv flag to denote deep convection (0=no, 1=yes)
!!  \param[in] islimsk sea/land/ice mask (=0/1/2)
!!  \param[in] dot layer mean vertical velocity (Pa/s)
!!  \param[in] ncloud number of cloud species
!!  \param[in] hpbl PBL height (m)
!!  \param[in] heat surface sensible heat flux (K m/s)
!!  \param[in] evap surface latent heat flux (kg/kg m/s)
!!  \param[out] ud_mf updraft mass flux multiplied by time step (\f$kg/m^2\f$)
!!  \param[out] dt_mf ud_mf at cloud top (\f$kg/m^2\f$)
!!  \param[out] cnvw convective cloud water (kg/kg)
!!  \param[out] cnvc convective cloud cover (unitless)
!!  \param[in] clam coefficient for entrainment rate
!!  \param[in] c0s convective rain conversion parameter (1/m)
!!  \param[in] c1 conversion parameter of detrainment from liquid water into grid-scale cloud water (1/m)
!!  \param[in] pgcon reduction factor in momentum transport due to convection induced pressure gradient force
!!  \param[in] asolfac aerosol-aware parameter inversely proportional to CCN number concentraion
!!
!!  \section general General Algorithm
!!  -# Compute preliminary quantities needed for the static and feedback control portions of the algorithm.
!!  -# Perform calculations related to the updraft of the entraining/detraining cloud model ("static control").
!!  -# The cloud base mass flux is obtained using the cumulus updraft velocity averaged ove the whole cloud depth.
!!  -# Calculate the tendencies of the state variables (per unit cloud base mass flux) and the cloud base mass flux.
!!  -# For the "feedback control", calculate updated values of the state variables by multiplying the cloud base mass flux and the tendencies calculated per unit cloud base mass flux from the static control.
!!  \section detailed Detailed Algorithm
!!  @{
      subroutine samfshalcnv_kh_gpu(nxjp, ix, km, delt, delp, prslp, psp, phil, ql, &
                                   qi, q1, t1, u1, v1, rn, kbot, ktop, kcnv, islimsk, garea, &
                                   dot, ncloud, hpbl, cnvw, cnvc)
      !$acc routine(fpvs_gpu) seq
!byl           dot,ncloud,hpbl,ud_mf,dt_mf,cnvw,cnvc,
!    &     dot,ncloud,hpbl,ud_mf,dt_mf,cnvw,cnvc,me) &
!byl           clam,c0s,c1,pgcon,asolfac)
!
         use machine, only: kind_phys
!byl      use funcphys , only : fpvs
         use param, only: my, my_max
         use index, only: jlistnum, jlist1
         use rank
         use physcons, grav => con_g, cp => con_cp, hvap => con_hvap &
            , rv => con_rv, fv => con_fvirt, t0c => con_t0c &
            , rd => con_rd, cvap => con_cvap, cliq => con_cliq &
            , eps => con_eps, epsm1 => con_epsm1
         implicit none
!
         real fpvs, fpvs_gpu
!
         integer im, ix, km, ncloud, &
            kbot(ix, my_max), ktop(ix, my_max), kcnv(ix, my_max)
!    &,                  me
         real(kind=kind_phys) delt
         real(kind=kind_phys) psp(ix, my_max), delp(ix, km, my_max), prslp(ix, km, my_max)
         real(kind=kind_phys) ps(ix, my_max), del(ix, km, my_max), prsl(ix, km, my_max), &
            ql(ix, km, my_max), q1(ix, km, my_max), t1(ix, km, my_max), &
            u1(ix, km, my_max), v1(ix, km, my_max), &
            !    &                     u1(ix,km),  v1(ix,km),   rcs(ix), &
            rn(ix, my_max), garea(ix, my_max), &
            dot(ix, km, my_max), phil(ix, km, my_max), hpbl(ix, my_max), &
            cnvw(ix, km, my_max), cnvc(ix, km, my_max) &
            ! hchuang code change mass flux output &
            , ud_mf(ix, km, my_max), dt_mf(ix, km, my_max), qi(ix, km, my_max)
!
         integer i, j, indx, k, kk, km1, n, jj, jjj
         integer kpbl(ix, my_max), nxjp(my), myim(my_max)
         integer, dimension(ix, my_max), intent(in) :: islimsk
!
         real(kind=kind_phys) dellat, delta, &
            c0l, c0s, d0, &
            c1, asolfac, &
            desdt, dp, &
            dq, dqsdp, dqsdt, dt, &
            dt2, dtmax, dtmin, dxcrt, &
            dv1h, dv2h, dv3h, &
            dv1q, dv2q, dv3q, &
            dz, dz1, e1, clam, &
            el2orc, elocp, aafac, cm, &
            es, etah, h1, &
            evef, evfact, evfactl, fact1, &
            fact2, factor, dthk, &
            g, gamma, pprime, betaw, &
            qlk, qrch, qs, &
            rfact, shear, tfac, &
            val, val1, val2, &
            w1, w1l, w1s, w2, &
            w2l, w2s, w3, w3l, &
            w3s, w4, w4l, w4s, &
            rho, tem, tem1, tem2, &
            ptem, ptem1, &
            pgcon
!
         integer kb(ix, my_max), kbcon(ix, my_max), kbcon1(ix, my_max), &
            ktcon(ix, my_max), ktcon1(ix, my_max), ktconn(ix, my_max), &
            kbm(ix, my_max), kmax(ix, my_max)
!
         real(kind=kind_phys) aa1, cina, &
            umean, tauadv, gdx(ix, my_max), &
            delhbar(ix, my_max), delq, delq2, &
            delqbar(ix, my_max), delqev, deltbar(ix, my_max), &
            deltv, dtconv, edt(ix, my_max), &
            pdot(ix, my_max), po, &
            qcond, qevap, hmax, &
            rntot, vshear, &
            xlamud(ix, my_max), xmb(ix, my_max), xmbmax, &
            delubar(ix, my_max), delvbar(ix, my_max)
!
         real(kind=kind_phys) c0(ix, my_max)
!
         real(kind=kind_phys) crtlamd
!
         real(kind=kind_phys) cinpcr, cinpcrmx, cinpcrmn, &
            cinacr, cinacrmx, cinacrmn
!
!  parameters for updraft velocity calculation
         real(kind=kind_phys) bet1, cd1, f1, gam1, &
            bb1, bb2
!    &                     bb1,     bb2,     wucb
!cc
!  physical parameters
         parameter(g=grav, asolfac=0.958)
!byl      parameter(g=grav)
         parameter(elocp=hvap/cp, &
                   el2orc=hvap*hvap/(rv*cp))
         parameter(c0s=0.002, c1=5.e-4, d0=.01)
!byl      parameter(d0=.01)
!     parameter(c0l=c0s*asolfac)
!
! asolfac: aerosol-aware parameter based on Lim & Hong (2012)
!      asolfac= cx / c0s(=.002)
!      cx = min([-0.7 ln(Nccn) + 24]*1.e-4, c0s)
!      Nccn: CCN number concentration in cm^(-3)
!      Until a realistic Nccn is provided, Nccns are assumed
!      as Nccn=100 for sea and Nccn=7000 for land
!
         parameter(cm=1.0, delta=fv)
         parameter(fact1=(cvap - cliq)/rv, fact2=hvap/rv - fact1*t0c)
         parameter(dthk=25.)
         parameter(cinpcrmn=120.)
!      parameter(cinpcrmx=180.,cinpcrmn=120.)
!     parameter(cinacrmx=-120.,cinacrmn=-120.)
         parameter(cinacrmx=-120., cinacrmn=-80.)
         parameter(crtlamd=3.e-4)
         parameter(dtmax=10800., dtmin=600.)
         parameter(bet1=1.875, cd1=.506, f1=2.0, gam1=.5)
         parameter(betaw=.03, dxcrt=15.e3)
         parameter(h1=0.33333333)
!  local variables and arrays
         real(kind=kind_phys) pfld(ix, km, my_max), to(ix, km, my_max), qo(ix, km, my_max), &
            uo(ix, km, my_max), vo(ix, km, my_max), qeso(ix, km, my_max)
!  for updraft velocity calculation
         real(kind=kind_phys) wu2(ix, km, my_max), buo(ix, km, my_max), drag(ix, km, my_max)
         real(kind=kind_phys) wc(ix, my_max), scaldfunc, sigmagfm
!
!  cloud water
!     real(kind=kind_phys) qlko_ktcon(im), dellal(im,km), tvo(im,km),
         real(kind=kind_phys) qlko_ktcon(ix, my_max), dellal(ix, km, my_max), &
            dbyo(ix, km, my_max), zo(ix, km, my_max), xlamue(ix, km, my_max), &
            heo(ix, km, my_max), heso(ix, km, my_max), frh(ix, km, my_max), &
            dellah(ix, km, my_max), dellaq(ix, km, my_max), &
            dellau(ix, km, my_max), dellav(ix, km, my_max), hcko(ix, km,  my_max), &
            ucko(ix, km, my_max), vcko(ix, km, my_max), qcko(ix, km, my_max), &
            qrcko(ix, km, my_max), eta(ix, km, my_max), &
            zi(ix, km, my_max), pwo(ix, km, my_max), c0t(ix, km, my_max), &
            sumx, tx1(ix, my_max), cnvwt(ix, km, my_max)
!
         logical totflg, cnvflg(ix, my_max), flg(ix, my_max)
!
         real(kind=kind_phys) tf, tcr, tcrf
         parameter(tf=233.16, tcr=263.16, tcrf=1.0/(tcr - tf))
!<---for scale-aware parameterization (Kwon and Hong 2017)
         real(kind=kind_phys) cinpcri, frh_sum, wbar(ix, my_max), &
            clear(ix, my_max), sigma(ix, my_max)
         real(kind=kind_phys) po1, po2
         real(kind=kind_phys) sigma_con, pi, dx1km, dx5km, dx250m
         parameter(pi=3.14159)
         parameter(dx1km=1000., dx5km=5000., dx250m=250.)

         !for GPU porting
         integer, parameter :: nxpvs = 7501
         real :: c1xpvs, c2xpvs, tbpvs(nxpvs)
         common/fpvscom/c1xpvs, c2xpvs, tbpvs(nxpvs)
         
         !Register
         integer :: kbr, kmaxr, ktconr, kbmr, kbcon1r, kbconr, ktcon1r, ktconnr
         real :: heor, pfldr1, pfldr, zor, zor1, qesor, tor1, qor, por, uor, &
            uor1, vor, vor1, zir, zir1, dragr, dragr1, buor, buor1, wu2r, &
            qckor, etar, etar1, c0tr, qor, qor1, xlamuer, xlamuer1, hckor1, &
            uckor1, vckor1, heor1, xlamudr, wcr, wbarr, pfldrm1, pfldrp1, &
            wu2r1, dotr, dotr1, tor, rnr, xmbr, edtr, islimskr, delqbarr, deltbarr

         logical :: flgr

         integer :: async_id = 1

!
!c-----------------------------------------------------------------------
!
!************************************************************************
!     convert input Pa terms to Cb terms  -- Moorthi
!>  ## Compute preliminary quantities needed for the static and feedback control portions of the algorithm.
!>  - Convert input pressure terms to centibar units.
!byl      ps   = psp   * 0.001
!byl      prsl = prslp * 0.001
!byl      del  = delp  * 0.001
         
         !$acc data create(myim, ps, prsl, del, cnvflg, ktop, kbcon, ktcon, &
         !$acc&     ktconn, kb, pdot, qlko_ktcon, edt, gdx, &
         !$acc&     sigma, c0, c0t, ud_mf, dt_mf, kbm, kmax, &
         !$acc&     tx1, zo, zi, xlamue, flg, kpbl, pfld, eta, frh, hcko, qcko, &
         !$acc&     qrcko, ucko, vcko, dbyo, pwo, dellal, to, qo, uo, vo, &
         !$acc&     wu2, buo, drag, cnvwt, qeso, heo, heso, &
         !$acc&     xlamud, kbcon1, ktcon1, wc, wbar, clear, &
         !$acc&     dellah, dellaq, dellau, dellav, xmb, &
         !$acc&     delhbar, delqbar, deltbar, delubar, delvbar &
         !$acc&     ) async(async_id)
         
         !$acc parallel loop async(async_id) private(j)
         do jj = 1, jlistnum
            j = jlist1(jj)
            myim(jj) = nxjp(j)
         end do
         
         cinpcrmx = 240.
         !$acc parallel loop gang collapse(2) async(async_id)
         do jj = 1, jlistnum
            do k = 1, km
               !$acc loop vector
               do i = 1, myim(jj)
                  prsl(i, k, jj) = prslp(i, k, jj)
                  del(i, k, jj) = delp(i, k, jj)
               end do
            end do
         end do
   !************************************************************************
   !
            km1 = km - 1
         !$acc parallel loop gang vector collapse(2) async(async_id)
         do jj = 1, jlistnum
            do i = 1, ix
               if (i .le. myim(jj)) then
                  ps(i, jj) = psp(i, jj)
   !
   !  initialize arrays
   !
   !>  - Initialize column-integrated and other single-value-per-column variable arrays.
                  cnvflg(i, jj) = .true.
                  if (kcnv(i, jj) == 1) cnvflg(i, jj) = .false.
                  if (cnvflg(i, jj)) then
                     kbot(i, jj) = km + 1
                     ktop(i, jj) = 0
                  end if
                  rn(i, jj) = 0.
                  kbcon(i, jj) = km
                  ktcon(i, jj) = 1
                  ktconn(i, jj) = 1
                  kb(i, jj) = km
                  pdot(i, jj) = 0.
                  qlko_ktcon(i, jj) = 0.
                  edt(i, jj) = 0.
                  gdx(i, jj) = sqrt(garea(i, jj))

                  sigma(i, jj) = 0.
   !!
   !>  - Return to the calling routine if deep convection is present or the surface buoyancy flux is negative.
            !totflg = .true.
            !do i = 1, myim(jj)
            !   totflg = totflg .and. (.not. cnvflg(i, jj))
            !end do
            !if (totflg) cycle
   !!
   !>  - determine aerosol-aware rain conversion parameter over land
         
                  if (islimsk(i, jj) == 1) then
                     c0(i, jj) = c0s*asolfac
                  else
                     c0(i, jj) = c0s
                  end if
               end if
            end do
         end do
   !
   !>  - determine rain conversion parameter above the freezing level which exponentially decreases with decreasing temperature from Han et al.'s (2017) \cite han_et_al_2017 equation 8.
         
         !$acc parallel loop gang collapse(2) async(async_id) private(tem, tem1)
         do jj = 1, jlistnum
            do k = 1, km
               !$acc loop vector
               do i = 1, myim(jj)
                  if (t1(i, k, jj) > 273.16) then
                     c0t(i, k, jj) = c0(i, jj)
                  else
                     tem = d0*(t1(i, k, jj) - 273.16)
                     tem1 = exp(tem)
                     c0t(i, k, jj) = c0(i, jj)*tem1
                  end if
               end do
            end do
         end do
   !
   !
            dt2 = delt
   !
   !  model tunable parameters are all here
            clam = .3
            aafac = .1
   !     evef    = 0.07
            evfact = 0.3
            evfactl = 0.3
   !
   !     pgcon   = 0.7     ! Gregory et al. (1997, QJRMS)
            pgcon = 0.55    ! Zhang & Wu (2003,JAS)
   !
            w1l = -8.e-3
            w2l = -4.e-2
            w3l = -5.e-3
            w4l = -5.e-4
            w1s = -2.e-4
            w2s = -2.e-3
            w3s = -1.e-3
            w4s = -2.e-5
   !
   !  define top layer for search of the downdraft originating layer
   !  and the maximum thetae for updraft
   !
   !>  - Determine maximum indices for the parcel starting point (kbm) and cloud top (kmax).
         
         !$acc parallel loop gang vector collapse(2) async(async_id)
         do jj = 1, jlistnum
            do i = 1, ix
               if (i .le. myim(jj)) then
                  kbm(i, jj) = km
                  kmax(i, jj) = km
                  tx1(i, jj) = 1.0/ps(i, jj)
   !
   
                  !$acc loop seq
                  do k = 1, km
                     if (prsl(i, k, jj)*tx1(i, jj) > 0.70) kbm(i, jj) = k + 1
                     if (prsl(i, k, jj)*tx1(i, jj) > 0.60) kmax(i, jj) = k + 1
                  end do
                  kbm(i, jj) = min(kbm(i, jj), kmax(i, jj))
               end if
            end do
         end do
   !
   !  hydrostatic height assume zero terr and compute
   !  updraft entrainment rate as an inverse function of height
   !
   !>  - Calculate hydrostatic height at layer centers assuming a flat surface (no terrain) from the geopotential.
         !$acc parallel loop gang collapse(2) async(async_id)
         do jj = 1, jlistnum
            do k = 1, km
               !$acc loop vector
               do i = 1, myim(jj)
                  zo(i, k, jj) = phil(i, k, jj)/g
               end do
            end do
         end do
   !>  - Calculate interface height and the entrainment rate as an inverse function of height.
         !$acc parallel loop gang collapse(2) async(async_id)
         do jj = 1, jlistnum
            do k = 1, km1
               !$acc loop vector
               do i = 1, myim(jj)
                  zi(i, k, jj) = 0.5*(zo(i, k, jj) + zo(i, k + 1, jj))
                  xlamue(i, k, jj) = clam/zi(i, k, jj)
               end do
            end do
         end do
         
         !$acc parallel loop gang vector collapse(2) async(async_id)
         do jj = 1, jlistnum
            do i = 1, ix
               if (i .le. myim(jj)) then
                  xlamue(i, km, jj) = xlamue(i, km1, jj)
   !
   !  pbl height
   !
   !>  - Find the index for the PBL top using the PBL height; enforce that it is lower than the maximum parcel starting level.
                  flg(i, jj) = cnvflg(i, jj)
                  kpbl(i, jj) = 1
                  !$acc loop seq
                  do k = 2, km1
                     if (flg(i, jj) .and. zo(i, k, jj) <= hpbl(i, jj)) then
                        kpbl(i, jj) = k
                     else
                        flg(i, jj) = .false.
                     end if
                  end do

                  kpbl(i, jj) = min(kpbl(i, jj), kbm(i, jj))
               end if
            end do
         end do
   !
   !c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !   convert surface pressure to mb from cb
   !
   !>  - Convert prsl from centibar to millibar, set normalized mass flux to 1, cloud properties to 0, and save model state variables (after advection/turbulence).
         !$acc parallel loop gang collapse(2) async(async_id)
         do jj = 1, jlistnum
            do k = 1, km
               !$acc loop vector
               do i = 1, myim(jj)
                  if (cnvflg(i, jj) .and. k <= kmax(i, jj)) then
                     pfld(i, k, jj) = prsl(i, k, jj)*10.0
                     eta(i, k, jj) = 1.
                     hcko(i, k, jj) = 0.
                     qcko(i, k, jj) = 0.
                     qrcko(i, k, jj) = 0.
                     ucko(i, k, jj) = 0.
                     vcko(i, k, jj) = 0.
                     dbyo(i, k, jj) = 0.
                     pwo(i, k, jj) = 0.
                     dellal(i, k, jj) = 0.
                     to(i, k, jj) = t1(i, k, jj)
                     qo(i, k, jj) = q1(i, k, jj)
                     uo(i, k, jj) = u1(i, k, jj)
                     vo(i, k, jj) = v1(i, k, jj)
   !           uo(i,k)   = u1(i,k) * rcs(i)
   !           vo(i,k)   = v1(i,k) * rcs(i)
                     wu2(i, k, jj) = 0.
                     buo(i, k, jj) = 0.
                     drag(i, k, jj) = 0.
                     cnvwt(i, k, jj) = 0.
                  end if
               end do
            end do
         end do
   !>  - Calculate saturation specific humidity and enforce minimum moisture values.
         !$acc parallel loop gang collapse(2) async(async_id) private(val1, val2, &
         !$acc&         tem, qesor)
         do jj = 1, jlistnum
            do k = 1, km
               !$acc loop vector
               do i = 1, myim(jj)
                  if (cnvflg(i, jj) .and. k <= kmax(i, jj)) then
                     qesor = 0.01*fpvs_gpu(to(i, k, jj),c1xpvs,c2xpvs,tbpvs)      ! fpvs is in pa
   !byl            call qsatq_cwb(1,to(i,k),pfld(i,k),qeso(i,k))
                     qesor = eps*qesor/(pfld(i, k, jj) + epsm1*qesor)
                     val1 = 1.e-8
                     qesor = max(qesor, val1)
                     val2 = 1.e-10
                     qo(i, k, jj) = max(qo(i, k, jj), val2)
   !           qo(i,k)   = min(qo(i,k),qeso(i,k))
   !           tvo(i,k)  = to(i,k) + delta * to(i,k) * qo(i,k)
   !
   !  compute moist static energy
   !
   !>  - Calculate moist static energy (heo) and saturation moist static energy (heso).
   !           tem       = g * zo(i,k) + cp * to(i,k)
                     tem = phil(i, k, jj) + cp*to(i, k, jj)
                     heo(i, k, jj) = tem + hvap*qo(i, k, jj)
                     heso(i, k, jj) = tem + hvap*qesor
                     qeso(i, k, jj) = qesor
   !           heo(i,k)  = min(heo(i,k),heso(i,k))
                  end if
               end do
            end do
         end do
   !
   !  determine level with largest moist static energy within pbl
   !  this is the level where updraft starts
   !
   !> ## Perform calculations related to the updraft of the entraining/detraining cloud model ("static control").
   !> - Search in the PBL for the level of maximum moist static energy to start the ascending parcel.
         !$acc parallel loop gang vector collapse(2) async(async_id) private( &
         !$acc&         dz, dp, es, pprime, qs, dqsdp, desdt, dqsdt, gamma, dt, &
         !$acc&         dq, val1, val2, ptem, hmax, kbr, kmaxr, pfldr1, pfldr, &
         !$acc&         zor, zor1, qesor, tor1, qor, por, uor, uor1, vor, vor1)
         do jj = 1, jlistnum
            do i = 1, ix
               if (i .le. myim(jj)) then
                  if (cnvflg(i, jj)) then
                     hmax = heo(i, 1, jj)
                     kbr = 1
                     kmaxr = kmax(i, jj)

                     !$acc loop seq
                     do k = 2, km
                        heor = heo(i, k, jj)
                        if (heor > hmax .and. k <= kpbl(i, jj)) then
                           kbr = k
                           hmax = heor
                        end if
                     end do
                     kb(i, jj) = kbr
   !
   !> - Calculate the temperature, water vapor mixing ratio, and pressure at interface levels.
                     pfldr = pfld(i, 1, jj)
                     zor = zo(i, 1, jj)
                     uor = uo(i, 1, jj)
                     vor = vo(i, 1, jj)
                     frh(i, kmaxr, jj) = 0.
                     !$acc loop seq
                     do k = 1, km1
                        if (k <= kmaxr - 1) then
                           pfldr1 = pfld(i, k + 1, jj)
                           zor1 = zo(i, k + 1, jj)
                           tor1 = to(i, k + 1, jj)
                           qor = qo(i, k + 1, jj)
                           uor1 = uo(i, k + 1, jj)
                           vor1 = vo(i, k + 1, jj)
                           dz = .5*(zor1 - zor)
                           dp = .5*(pfldr1 - pfldr)
                           es = 0.01*fpvs_gpu(tor1,c1xpvs,c2xpvs,tbpvs)      ! fpvs is in pa
         !byl            call qsatq_cwb(1,to(i,k+1),pfld(i,k+1),es)
                           pprime = pfldr1 + epsm1*es
                           qs = eps*es/pprime
                           dqsdp = -qs/pprime
                           desdt = es*(fact1/tor1 + fact2/(tor1**2))
                           dqsdt = qs*pfldr1*desdt/(es*pprime)
                           gamma = el2orc*qeso(i, k + 1, jj)/(tor1**2)
                           dt = (g*dz + hvap*dqsdp*dp)/(cp*(1.+gamma))
                           dq = dqsdt*dt + dqsdp*dp
                           tor1 = tor1 + dt
                           to(i, k, jj) = tor1
                           qor = qor + dq
                           por = .5*(pfldr + pfldr1)
                           pfldr = pfldr1
   !
   !> - Recalculate saturation specific humidity, moist static energy, saturation moist static energy, and horizontal momentum on interface levels. Enforce minimum specific humidity.
         
                           qesor = 0.01*fpvs_gpu(tor1,c1xpvs,c2xpvs,tbpvs)      ! fpvs is in pa
         !byl            call qsatq_cwb(1,to(i,k),po(i,k),qeso(i,k))
                           qesor = eps*qesor/(por + epsm1*qesor)
                           val1 = 1.e-8
                           qesor = max(qesor, val1)
                           val2 = 1.e-10
                           qor = max(qor, val2)
                           qo(i, k, jj) = qor
         !           qo(i,k)   = min(qo(i,k),qeso(i,k))
                           frh(i, k, jj) = 1.-min(qor/qesor, 1.)
                           heo(i, k, jj) = .5*g*(zor + zor1) + &
                                       cp*tor1 + hvap*qor
                           heso(i, k, jj) = .5*g*(zor + zor1) + &
                                        cp*tor1 + hvap*qesor
                           uo(i, k, jj) = .5*(uor + uor1)
                           vo(i, k, jj) = .5*(vor + vor1)
                           zor = zor1
                           qeso(i, k ,jj) = qesor
                           uor = uor1
                           vor = vor1
                        end if
                     end do
                  end if
               end if
            end do
         end do
   !
   !  look for the level of free convection as cloud base
   !
   !> - Search below the index "kbm" for the level of free convection (LFC) where the condition \f$h_b > h^*\f$ is first met, where \f$h_b, h^*\f$ are the state moist static energy at the parcel's starting level and saturation moist static energy, respective
         !$acc parallel loop gang vector collapse(2) async(async_id) private( &
         !$acc&         sigma_con, cinpcr, tem1, flgr, kbr, kbconr, xlamudr, &
         !$acc&         kbmr, ktconnr, zir1, zir, xlamuer1, xlamuer, etar1, kk)
         do jj = 1, jlistnum
            do i = 1, ix
               if (i .le. myim(jj)) then
                  kbconr = kbcon(i, jj)
                  kbr = kb(i, jj)
                  kbmr = kbm(i, jj)
                  kmaxr = kmax(i, jj)
                  if (cnvflg(i, jj)) then
                     flgr = .true.
                     if (flgr) kbconr = kmaxr

                     !$acc loop seq
                     do k = 2, km1
                        if (flgr .and. k < kbmr) then
                           if (k > kbr .and. heo(i, kbr, jj) > heso(i, k, jj)) then
                              kbconr = k
                              flgr = .false.
                           end if
                        end if
                     end do
   !
                     if (kbconr == kmaxr) cnvflg(i, jj) = .false.
                  end if
   !!
   !> - If no LFC, return to the calling routine without modifying state variables.
            !totflg = .true.
            !do i = 1, myim(jj)
            !   totflg = totflg .and. (.not. cnvflg(i, jj))
            !end do
            !if (totflg) cycle
   !!
   !> - Determine the vertical pressure velocity at the LFC. After Han and Pan (2011) \cite han_and_pan_2011 , determine the maximum pressure thickness between a parcel's starting level and the LFC. If a parcel doesn't reach the LFC within the critical thick
                  if (cnvflg(i, jj)) then
                     pdot(i, jj) = 10.*dot(i, kbconr, jj)
      !          pdot(i)  = 0.01 * dot(i,kbcon(i)) ! Now dot is in Pa/s

                     frh_sum = 0.
                     !$acc loop seq
                     do k = 1, km1
                        if (k .ge. kbr .and. k .le. kbconr) then
                           frh_sum = frh_sum + (1 - frh(i, k, jj))
                        end if
                     end do
   !
   !   turn off convection if pressure depth between parcel source level
   !      and cloud base is larger than a critical value, cinpcr
   !
                     sigma_con = tan(0.4*pi)/(dx5km - dx1km)                     !7.7 e-4 m-1
                     sigma(i, jj) = (1.-1./pi*(atan(sigma_con*(gdx(i, jj) - dx5km)) + pi/2.)) !1(1km),0.1(10km)
                     if (gdx(i, jj) .lt. dx5km) then
                        sigma(i, jj) = min(sigma(i, jj) - 0.01684*gdx(i, jj)/1000.+0.0842, 1.0)
                     end if
                     cinpcr = cinpcrmn + 0.5*(cinpcrmx - cinpcrmn)*(1.-sigma(i, jj))
                     cinpcri = cinpcr*frh_sum/(kbconr - kbr + 1)
                     tem1 = pfld(i, kbr, jj) - pfld(i, kbconr, jj)
                     if (tem1 > cinpcri) then
                        cnvflg(i, jj) = .false.
                     end if
                  end if
   !!
            !totflg = .true.
            !do i = 1, myim(jj)
            !   totflg = totflg .and. (.not. cnvflg(i, jj))
            !end do
            !if (totflg) cycle
   !!
   !
   !  specify the detrainment rate for the updrafts
   !
   !> - The updraft detrainment rate is set constant and equal to the entrainment rate at cloud base.
         
                  if (cnvflg(i, jj)) then
                     xlamudr = xlamue(i, kbconr, jj)
                     xlamud(i, jj) = xlamudr
      !         xlamud(i) = crtlamd
   !
   !  determine updraft mass flux for the subcloud layers
   !
   !> - Calculate the normalized mass flux for subcloud and in-cloud layers according to Pan and Wu (1995) \cite pan_and_wu_1995 equation 1:
   !!  \f[
   !!  \frac{1}{\eta}\frac{\partial \eta}{\partial z} = \lambda_e - \lambda_d
   !!  \f]
   !!  where \f$\eta\f$ is the normalized mass flux, \f$\lambda_e\f$ is the entrainment rate and \f$\lambda_d\f$ is the detrainment rate. The normalized mass flux increases upward below the cloud base and decreases upward above.
                     
                     kk = min(kbconr, km1)
                     kk = max(kk, 1)
                     zir1 = zi(i, kk, jj)
                     xlamuer1 = xlamue(i, kk, jj)
                     etar1 = eta(i, kk, jj)
                     !$acc loop seq
                     do k = km1, 1, -1
                        if (k < kbconr .and. k >= kbr) then
                           zir = zi(i, k, jj)
                           xlamuer = xlamue(i, k, jj)
                           dz = zir1 - zir
                           ptem = 0.5*(xlamuer + xlamuer1) - xlamudr
                           etar1 = etar1/(1.+ptem*dz)
                           eta(i, k, jj) = etar1
                           zir1 = zir
                           xlamuer1 = xlamuer
                        end if
                     end do
   !
   !  compute mass flux above cloud base
   !
                     flgr = .true.
                     kk = max(2, kbconr)
                     kk = min(kk, km1)
                     zir1 = zi(i, kk, jj)
                     xlamuer1 = xlamue(i, kk, jj)
                     etar1 = eta(i, kk, jj)
                     !$acc loop seq
                     do k = 2, km1
                        if (flgr) then
                           if (k > kbconr .and. k < kmaxr) then
                              zir = zi(i, k, jj)
                              xlamuer = xlamue(i, k, jj)
                              dz = zir - zir1
                              ptem = 0.5*(xlamuer + xlamuer1) - xlamudr
                              etar1 = etar1*(1 + ptem*dz)
                              eta(i, k, jj) = etar1
                              if (etar1 <= 0.) then
                                 kmaxr = k
                                 ktconnr = k
                                 kbmr = min(kbmr, kmaxr)
                                 flgr = .false.
                              end if
                              zir1 = zir
                              xlamuer1 = xlamuer
                           end if
                        end if
                     end do
                  end if
                  kbm(i, jj) = kbmr
                  ktconn(i, jj) = ktconnr
                  kmax(i, jj) = kmaxr
                  kbcon(i, jj) = kbconr
               end if
            end do
         end do
   !
   !  compute updraft cloud property
   !
   !> - Set cloud properties equal to the state variables at updraft starting level (kb).
         !$acc parallel loop gang vector collapse(2) async(async_id) private( &
         !$acc&         kbr, flgr, kbcon1r, kbmr, kbconr, hckor1, uckor1, vckor1, &
         !$acc&         zir, zir1, xlamuer, xlamuer1, vor, vor1, uor, uor1, heor, &
         !$acc&         heor1, gamma, rfact, qlk, val, ktconr, &
         !$acc&         w1, w2, w3, w4, tem, val1, val2, tem1, cinacr, k, dp, &
         !$acc&         dz, qrch, factor, dq, etah, ptem, dz1, ptem1, &
         !$acc&         ktcon1r, kk, dragr, dragr1, buor, buor1, &
         !$acc&         wu2r, qckor, etar, etar1, c0tr, qor, qor1, &
         !$acc&         shear, e1, sumx, po1, po2, dotr, dotr1, tor, tor1, &
         !$acc&         wcr, wbarr, pfldrm1, pfldr, pfldrp1, wu2r1, vshear)
         do jj = 1, jlistnum
            do i = 1, ix
               if (i .le. myim(jj)) then
                  kbr = kb(i, jj)
                  if (cnvflg(i, jj)) then
                     heor1 = heo(i, kbr, jj)
                     uor1 = uo(i, kbr, jj)
                     vor1 = vo(i, kbr, jj)
                     hckor1 = heor1
                     uckor1 = uor1
                     vckor1 = vor1
                     hcko(i, kbr, jj) = hckor1
                     ucko(i, kbr, jj) = uckor1
                     vcko(i, kbr, jj) = vckor1
   !
   !  cm is an enhancement factor in entrainment rates for momentum
   !
   !> - Calculate the cloud properties as a parcel ascends, modified by entrainment and detrainment. Discretization follows Appendix B of Grell (1993) \cite grell_1993 . Following Han and Pan (2006) \cite han_and_pan_2006, the convective momentum transport i
                     zir1 = zi(i, kbr, jj)
                     xlamuer1 = xlamue(i, kbr, jj)
                     !$acc loop seq
                     do k = 2, km1
                        if (k > kbr .and. k < kmax(i, jj)) then
                           xlamuer = xlamue(i, k, jj)
                           zir = zi(i, k, jj)
                           heor = heo(i, k, jj)
                           uor = uo(i, k, jj)
                           vor = vo(i, k, jj)
                           
                           dz = zir - zir1
                           tem = 0.5*(xlamuer + xlamuer1)*dz
                           tem1 = 0.5*xlamud(i, jj)*dz
                           factor = 1.+tem - tem1
                           hckor1 = ((1.-tem1)*hckor1 + tem*0.5* &
                                         (heor + heor1))/factor
                           dbyo(i, k, jj) = hckor1 - heso(i, k, jj)
                           hcko(i, k, jj) = hckor1
      !
                           tem = 0.5*cm*tem
                           factor = 1.+tem
                           ptem = tem + pgcon
                           ptem1 = tem - pgcon
                           uckor1 = ((1.-tem)*uckor1 + ptem*uor &
                                         + ptem1*uor1)/factor
                           vckor1 = ((1.-tem)*vckor1 + ptem*vor &
                                         + ptem1*vor1)/factor
                                         
                           ucko(i, k, jj) = uckor1
                           vcko(i, k, jj) = vckor1
                           zir1 = zir
                           xlamuer1 = xlamuer
                           heor1 = heor
                           uor1 = uor
                           vor1 = vor
                        end if
                     end do
   !
   !   taking account into convection inhibition due to existence of
   !    dry layers below cloud base
   !
   !> - With entrainment, recalculate the LFC as the first level where buoyancy is positive. The difference in pressure levels between LFCs calculated with/without entrainment must be less than a threshold (currently 25 hPa). Otherwise, convection is inhibit
                     flgr = .true.
                     kbmr = kbm(i, jj)
                     kbcon1r = kmax(i, jj)
                     kbconr = kbcon(i, jj)
                     
                  
                     !$acc loop seq
                     do k = 2, km1
                        if (flgr .and. k < kbmr) then
                           if (k >= kbconr .and. dbyo(i, k, jj) > 0.) then
                              kbcon1r = k
                              flgr = .false.
                           end if
                        end if
                     end do
                     kbcon1(i, jj) = kbcon1r

                     if (kbcon1r == kmax(i, jj)) cnvflg(i, jj) = .false.
                  end if

                  if (cnvflg(i, jj)) then
                     tem = pfld(i, kbconr, jj) - pfld(i, kbcon1r, jj)
                     if (tem > dthk) then
                        cnvflg(i, jj) = .false.
                     end if
                  end if
   !!
            !totflg = .true.
            !do i = 1, myim(jj)
            !   totflg = totflg .and. (.not. cnvflg(i, jj))
            !end do
            !if (totflg) cycle
   !!
   !
   !  calculate convective inhibition
   !
   !> - Calculate additional trigger condition of the convective inhibition (CIN) according to Han et al.'s (2017) \cite han_et_al_2017 equation 13.
                  ktconr = ktcon(i, jj)
                  if (cnvflg(i, jj)) then
                     kbmr = kbm(i, jj)
                     kbcon1r = kbcon1(i, jj)
                     kbconr = kbcon(i, jj)
                     cina = 0.
                     !$acc loop seq
                     do k = 2, km1
                        if (k > kbr .and. k < kbcon1r) then
                           dz1 = zo(i, k + 1, jj) - zo(i, k, jj)
                           gamma = el2orc*qeso(i, k, jj)/(to(i, k, jj)**2)
                           rfact = 1.+delta*cp*gamma &
                                   *to(i, k, jj)/hvap
                           cina = cina + &
                                     !    &                 dz1 * eta(i,k) * (g / (cp * to(i,k))) &
                                     dz1*(g/(cp*to(i, k, jj))) &
                                     *dbyo(i, k, jj)/(1.+gamma) &
                                     *rfact
                           val = 0.
                           cina = cina + &
                                     !    &                 dz1 * eta(i,k) * g * delta * &
                                     dz1*g*delta* &
                                     max(val, (qeso(i, k, jj) - qo(i, k, jj)))
                        end if
                     end do
   !> - Turn off convection if the CIN is less than a critical value (cinacr) which is inversely proportional to the large-scale vertical velocity.
                     if (islimsk(i, jj) == 1) then
                        w1 = w1l
                        w2 = w2l
                        w3 = w3l
                        w4 = w4l
                     else
                        w1 = w1s
                        w2 = w2s
                        w3 = w3s
                        w4 = w4s
                     end if
                     if (pdot(i, jj) <= w4) then
                        tem = (pdot(i, jj) - w4)/(w3 - w4)
                     elseif (pdot(i, jj) >= -w4) then
                        tem = -(pdot(i, jj) + w4)/(w4 - w3)
                     else
                        tem = 0.
                     end if

                     val1 = -1.
                     tem = max(tem, val1)
                     val2 = 1.
                     tem = min(tem, val2)
                     tem = 1.-tem
                     tem1 = .5*(cinacrmx - cinacrmn)
                     cinacr = cinacrmx - tem*tem1
      !
      !         cinacr = cinacrmx
                     if (cina < cinacr) cnvflg(i, jj) = .false.
                  end if
   !!
            !totflg = .true.
            !do i = 1, myim(jj)
            !   totflg = totflg .and. (.not. cnvflg(i, jj))
            !end do
            !if (totflg) cycle
   !!
   !
   !  determine first guess cloud top as the level of zero buoyancy
   !    limited to the level of P/Ps=0.7
   !
   !> - Calculate the cloud top as the first level where parcel buoyancy becomes negative; the maximum possible value is at \f$p=0.7p_{sfc}\f$.
                  
                  if (cnvflg(i, jj)) then
                     flgr = .true.
                     if (flgr) ktconr = kbmr

                     !$acc loop seq
                     do k = 2, km1
                        if (flgr .and. k < kbmr) then
                           if (k > kbcon1r .and. dbyo(i, k, jj) < 0.) then
                              ktconr = k
                              flgr = .false.
                           end if
                        end if
                     end do
   !
   !  specify upper limit of mass flux at cloud base
   !
   !> - Calculate the maximum value of the cloud base mass flux using the CFL-criterion-based formula of Han and Pan (2011) \cite han_and_pan_2011, equation 7.
      !         xmbmax(i) = .1
      !
                     !k = kbconr
                     !dp = 1000.*del(i, k, jj)
                     !xmbmax(i, jj) = dp/(g*dt2)
      !
      !         tem = dp / (g * dt2)
      !         xmbmax(i) = min(tem, xmbmax(i))
   !
   !  compute cloud moisture property and precipitation
   !
   !> - Set cloud moisture property equal to the enviromental moisture at updraft starting level (kb).
                     qor1 = qo(i, kbr, jj)
                     qcko(i, kbr, jj) = qor1
                     qrcko(i, kbr, jj) = qor1
                     zir1 = zi(i, kbr, jj)
                     etar1 = eta(i, kbr, jj)
                     qckor = qor1
   !> - Calculate the moisture content of the entraining/detraining parcel (qcko) and the value it would have if just saturated (qrch), according to equation A.14 in Grell (1993) \cite grell_1993 . Their difference is the amount of convective cloud water (ql
                     xlamuer1 = xlamue(i, kbr, jj)
                     !$acc loop seq
                     do k = 2, km1
                        if (k > kbr .and. k < ktconr) then
                           buor = 0.
                           xlamuer = xlamue(i, k, jj)
                           qor = qo(i, k, jj)
                           zir = zi(i, k, jj)
                           etar = eta(i, k, jj)
                           dz = zir - zir1
                           gamma = el2orc*qeso(i, k, jj)/(to(i, k, jj)**2)
                           qrch = qeso(i, k, jj) &
                                  + gamma*dbyo(i, k, jj)/(hvap*(1.+gamma))
      !cj
                           tem = 0.5*(xlamuer + xlamuer1)*dz
                           tem1 = 0.5*xlamud(i, jj)*dz
                           factor = 1.+tem - tem1
                           qckor = ((1.-tem1)*qckor + tem*0.5* &
                                         (qor + qor1))/factor
                           qrcko(i, k, jj) = qckor
      !cj
                           dq = etar*(qckor - qrch)
      !
      !             rhbar(i) = rhbar(i) + qo(i,k) / qeso(i,k)
      !
      !  below lfc check if there is excess moisture to release latent heat
      !
                           if (k >= kbconr .and. dq > 0.) then
                              etah = .5*(etar + etar1)
                              dp = 1000.*del(i, k, jj)
                              if (ncloud > 0) then
                                 ptem = c0t(i, k, jj) + c1
                                 qlk = dq/(etar + etah*ptem*dz)
                                 dellal(i, k, jj) = etah*c1*dz*qlk*g/dp
                              else
                                 qlk = dq/(etar + etah*c0t(i, k, jj)*dz)
                              end if
                              buor = buor - g*qlk
                              qckor = qlk + qrch
                              pwo(i, k, jj) = etah*c0t(i, k, jj)*dz*qlk
                              cnvwt(i, k, jj) = etah*qlk*g/dp
                           end if
      !
      !  compute buoyancy and drag for updraft velocity
      !
                           if (k >= kbconr) then
                              rfact = 1.+delta*cp*gamma &
                                      *to(i, k, jj)/hvap
                              buor = buor + (g/(cp*to(i, k, jj))) &
                                          *dbyo(i, k, jj)/(1.+gamma) &
                                          *rfact
                              val = 0.
                              buor = buor + g*delta* &
                                          max(val, (qeso(i, k, jj) - qor))
                              drag(i, k, jj) = max(xlamuer, xlamud(i, jj))
                           end if
                           xlamuer1 = xlamuer
                           qor1 = qor
                           zir1 = zir
                           etar1 = etar
                           qcko(i, k, jj) = qckor
                           buo(i, k, jj) = buor
      !
                        end if
                     end do
   !
   !  calculate cloud work function
   !
   !     do k = 2, km1
   !       do i = 1, im
   !         if (cnvflg(i)) then
   !           if(k >= kbcon(i) .and. k < ktcon(i)) then
   !             dz1 = zo(i,k+1) - zo(i,k)
   !             gamma = el2orc * qeso(i,k) / (to(i,k)**2)
   !             rfact =  1. + delta * cp * gamma
   !    &                 * to(i,k) / hvap
   !             aa1(i) = aa1(i) +
   !!   &                 dz1 * eta(i,k) * (g / (cp * to(i,k)))
   !    &                 dz1 * (g / (cp * to(i,k)))
   !    &                 * dbyo(i,k) / (1. + gamma)
   !    &                 * rfact
   !             val = 0.
   !             aa1(i) = aa1(i) +
   !!   &                 dz1 * eta(i,k) * g * delta *
   !    &                 dz1 * g * delta *
   !    &                 max(val,(qeso(i,k) - qo(i,k)))
   !           endif
   !         endif
   !       enddo
   !     enddo
   !     do i = 1, im
   !       if(cnvflg(i) .and. aa1(i) <= 0.) cnvflg(i) = .false.
   !     enddo
   !
   !  calculate cloud work function
   !
   !> - Calculate the cloud work function according to Pan and Wu (1995) \cite pan_and_wu_1995 equation 4:
   !!  \f[
   !!  A_u=\int_{z_0}^{z_t}\frac{g}{c_pT(z)}\frac{\eta}{1 + \gamma}[h(z)-h^*(z)]dz
   !!  \f]
   !! (discretized according to Grell (1993) \cite grell_1993 equation B.10 using B.2 and B.3 of Arakawa and Schubert (1974) \cite arakawa_and_schubert_1974 and assuming \f$\eta=1\f$) where \f$A_u\f$ is the updraft cloud work function, \f$z_0\f$ and \f$z_t\f
                     aa1 = 0.

                     kk = max(2, kbconr)
                     kk = min(kk, km1)
                     zor = zo(i, kk, jj)
                     !$acc loop seq
                     do k = 2, km1
                        if (k >= kbconr .and. k < ktconr) then
                           zor1 = zo(i, k + 1, jj)
                           dz1 = zor1 - zor
                           aa1 = aa1 + buo(i, k, jj)*dz1
                           zor = zor1
                        end if
                     end do

                     if (cnvflg(i, jj) .and. aa1 <= 0.) cnvflg(i, jj) = .false.
                  end if
   !!
   !> - If the updraft cloud work function is negative, convection does not occur, and the scheme returns to the calling routine.
            !totflg = .true.
            !do i = 1, myim(jj)
            !   totflg = totflg .and. (.not. cnvflg(i, jj))
            !end do
            !if (totflg) cycle
   !!
   !
   !  estimate the onvective overshooting as the level
   !    where the [aafac * cloud work function] becomes zero,
   !    which is the final cloud top
   !    limited to the level of P/Ps=0.7
   !
   !> - Continue calculating the cloud work function past the point of neutral buoyancy to represent overshooting according to Han and Pan (2011) \cite han_and_pan_2011 . Convective overshooting stops when \f$ cA_u < 0\f$ where \f$c\f$ is currently 10%, or w
                  wc(i, jj) = 0.
                  wbar(i, jj) = 0.
                  sumx = 0.
                  if (cnvflg(i, jj)) then
                     aa1 = aafac*aa1

                     flgr = .true.
                     ktcon1r = kbmr

                     kk = max(2, ktconr)
                     kk = min(kk, km1)
                     zor = zo(i, kk, jj)
                     !$acc loop seq
                     do k = 2, km1
                        if (flgr) then
                           if (k >= ktconr .and. k < kbmr) then
                              zor1 = zo(i, k + 1, jj)
                              tor1 = to(i, k, jj)
                              dz1 = zor1 - zor
                              gamma = el2orc*qeso(i, k, jj)/(tor1**2)
                              rfact = 1.+delta*cp*gamma &
                                      *tor1/hvap
                              aa1 = aa1 + &
                                       !    &                 dz1 * eta(i,k) * (g / (cp * to(i,k))) &
                                       dz1*(g/(cp*tor1)) &
                                       *dbyo(i, k, jj)/(1.+gamma) &
                                       *rfact
         !             val = 0.
         !             aa1(i) = aa1(i) +
         !!   &                 dz1 * eta(i,k) * g * delta *
         !    &                 dz1 * g * delta *
         !    &                 max(val,(qeso(i,k) - qo(i,k)))
                              if (aa1 < 0.) then
                                 ktcon1r = k
                                 flgr = .false.
                              end if
                              zor = zor1
                           end if
                        end if
                     end do
                     ktcon1(i, jj) = ktcon1r
   !
   !  compute cloud moisture property, detraining cloud water
   !    and precipitation in overshooting layers
   !
   !> - For the overshooting convection, calculate the moisture content of the entraining/detraining parcel as before. Partition convective cloud water and precipitation and detrain convective cloud water in the overshooting layers.

                     kk = max(2, ktconr)
                     kk = min(kk, km1)
                     zir1 = zi(i, kk - 1, jj)
                     qckor = qcko(i, kk - 1, jj)
                     etar1 = eta(i, kk - 1, jj)
                     qor1 = qo(i, kk - 1, jj)
                     xlamuer1 = xlamue(i, kk - 1, jj)
                     !$acc loop seq
                     do k = 2, km1
                        if (k >= ktconr .and. k < ktcon1r) then
                           xlamuer = xlamue(i, k, jj)
                           etar = eta(i, k, jj)
                           zir = zi(i, k, jj)
                           qor = qo(i, k, jj)
                           qesor = qeso(i, k, jj)
                           dz = zir - zir1
                           gamma = el2orc*qesor/(to(i, k, jj)**2)
                           qrch = qesor &
                                  + gamma*dbyo(i, k, jj)/(hvap*(1.+gamma))
      !cj
                           tem = 0.5*(xlamuer + xlamuer1)*dz
                           tem1 = 0.5*xlamud(i, jj)*dz
                           factor = 1.+tem - tem1
                           qckor = ((1.-tem1)*qckor + tem*0.5* &
                                         (qor + qor1))/factor
                           qrcko(i, k, jj) = qckor
      !cj
                           dq = etar*(qckor - qrch)
      !
      !  check if there is excess moisture to release latent heat
      !
                           if (dq > 0.) then
                              c0tr = c0t(i, k, jj)
                              etah = .5*(etar + etar1)
                              dp = 1000.*del(i, k, jj)
                              if (ncloud > 0) then
                                 ptem = c0tr + c1
                                 qlk = dq/(etar + etah*ptem*dz)
                                 dellal(i, k, jj) = etah*c1*dz*qlk*g/dp
                              else
                                 qlk = dq/(etar + etah*c0tr*dz)
                              end if
                              qckor = qlk + qrch
                              pwo(i, k, jj) = etah*c0tr*dz*qlk
                              cnvwt(i, k, jj) = etah*qlk*g/dp
                           end if
                           qcko(i, k, jj) = qckor
                           zir1 = zir
                           qor1 = qor
                           etar1 = etar
                           xlamuer1 = xlamuer
                        end if
                     end do
   !
   !  compute updraft velocity square(wu2)
   !> - Calculate updraft velocity square(wu2) according to Han et al.'s (2017) \cite han_et_al_2017 equation 7.
   !
   !     bb1 = 2. * (1.+bet1*cd1)
   !     bb2 = 2. / (f1*(1.+gam1))
   !
   !     bb1 = 3.9
   !     bb2 = 0.67
   !
   !     bb1 = 2.0
   !     bb2 = 4.0
   !
            bb1 = 4.0
            bb2 = 0.8
   !
   !     do i = 1, im
   !       if (cnvflg(i)) then
   !         k = kbcon1(i)
   !         tem = po(i,k) / (rd * to(i,k))
   !         wucb = -0.01 * dot(i,k) / (tem * g)
   !         if(wucb > 0.) then
   !           wu2(i,k) = wucb * wucb
   !         else
   !           wu2(i,k) = 0.
   !         endif
   !       endif
   !     enddo

                     kk = max(2, kbcon1r)
                     kk = min(kk, km1)
                     zir1 = zi(i, kk, jj)
                     dragr1 = drag(i, kk, jj)
                     buor1 = buo(i, kk, jj)
                     !$acc loop seq
                     do k = 2, km1
                        if (k > kbcon1r .and. k < ktconr) then
                           wu2r = wu2(i, k - 1, jj)
                           zir = zi(i, k, jj)
                           dragr = drag(i, k, jj)
                           buor = buo(i, k, jj)
                           dz = zir - zir1
                           tem = 0.25*bb1*(dragr + dragr1)*dz
                           tem1 = 0.5*bb2*(buor + buor1)*dz
                           ptem = (1.-tem)*wu2r
                           ptem1 = 1.+tem
                           wu2r = (ptem + tem1)/ptem1
                           wu2r = max(wu2r, 0.)
                           wu2(i, k, jj) = wu2r
                           zir1 = zir
                           dragr1 = dragr
                           buor1 = buor
                        end if
                     end do
                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   !  compute updraft velocity averaged over the whole cumulus
   !
   !> - Calculate the mean updraft velocity within the cloud (wc).
                     wcr = wc(i, jj)
                     wbarr = wbar(i, jj)
                     kk = max(2, kbcon1r)
                     kk = min(kk, km1)
                     pfldrm1 = pfld(i, kk, jj)
                     pfldr = pfld(i, kk + 1, jj)
                     zir1 = zi(i, kk, jj)
                     wu2r1 = wu2(i, kk, jj)
                     dotr1 = dot(i, kk, jj)
                     tor1 = to(i, kk, jj)
                     !$acc loop seq
                     do k = 2, km1
                        if (k > kbcon1r .and. k < ktconr) then
                           pfldrp1 = pfld(i, k + 1, jj)
                           zir = zi(i, k, jj)
                           wu2r = wu2(i, k, jj)
                           dotr = dot(i, k, jj)
                           tor = to(i, k, jj)
                           po1 = .5*(pfldr + pfldrp1)
                           po2 = .5*(pfldrm1 + pfldr)
                           dz = zir - zir1
                           tem = 0.5*(sqrt(wu2r) + sqrt(wu2r1))
                           wcr = wcr + tem*dz
                           tem = 10.*dotr*tor/po1
                           tem1 = 10.*dotr1*tor1/po2
                           wbarr = wbarr + (-0.5*rd/g)*(tem + tem1)*dz !grid-scale vertical velocity
                           sumx = sumx + dz
                           
                           pfldrm1 = pfldr
                           pfldr = pfldrp1
                           zir1 = zir
                           wu2r1 = wu2r
                           dotr1 = dotr
                           tor1 = tor
                        end if
                     end do
                     wbar(i, jj) = wbarr

                     if (sumx == 0.) then
                        cnvflg(i, jj) = .false.
                     else
                        wcr = wcr/sumx
                     end if
                     val = 1.e-4
                     if (wcr < val) cnvflg(i, jj) = .false.
                     wc(i, jj) = wcr
                  end if

                  if (cnvflg(i, jj)) then
                     if (sumx > 0.) then
                        wbarr = wbarr/sumx
                        wbar(i, jj) = wbarr
                     end if
   ! compute mean cloud core fraction
   ! assume mean cloud core fraction to be the ratio of
   ! mean grid-scale vertical velocity (wbar) and mean updraft velocity

                     tem = wbarr/wcr
                     tem = max(tem, 0.)
                     clear(i, jj) = 1.-tem
                     clear(i, jj) = max(min(clear(i, jj), 1.0), 0.)
                     if (wbarr .gt. 0. .and. wbarr .gt. wcr) then
                        cnvflg(i, jj) = .false.
                     end if
                  end if
   !
   ! exchange ktcon with ktcon1
   !
                  if (cnvflg(i, jj)) then
                     kk = ktconr
                     ktconr = ktcon1(i, jj)
                     ktcon1(i, jj) = kk
   !
   !  this section is ready for cloud water
   !
                     if (ncloud > 0) then
   !
   !  compute liquid and vapor separation at cloud top
   !
   !> - => Separate the total updraft cloud water at cloud top into vapor and condensate.
                        k = ktconr - 1
                        gamma = el2orc*qeso(i, k, jj)/(to(i, k, jj)**2)
                        qrch = qeso(i, k, jj) &
                               + gamma*dbyo(i, k, jj)/(hvap*(1.+gamma))
                        dq = qcko(i, k, jj) - qrch
      !
      !  check if there is excess moisture to release latent heat
      !
                        if (dq > 0.) then
                           qlko_ktcon(i, jj) = dq
                           qcko(i, k, jj) = qrch
                        end if
                     end if
   !
   !c--- compute precipitation efficiency in terms of windshear
   !
   !! - Calculate the wind shear and precipitation efficiency according to equation 58 in Fritsch and Chappell (1980) \cite fritsch_and_chappell_1980 :
   !! \f[
   !! E = 1.591 - 0.639\frac{\Delta V}{\Delta z} + 0.0953\left(\frac{\Delta V}{\Delta z}\right)^2 - 0.00496\left(\frac{\Delta V}{\Delta z}\right)^3
   !! \f]
   !! where \f$\Delta V\f$ is the integrated horizontal shear over the cloud depth, \f$\Delta z\f$, (the ratio is converted to units of \f$10^{-3} s^{-1}\f$). The variable "edt" is \f$1-E\f$ and is constrained to the range \f$[0,0.9]\f$.
                     vshear = 0.
         
                     !$acc loop seq
                     do k = 2, km
                        if (k > kbr .and. k <= ktconr) then
                           shear = sqrt((uo(i, k, jj) - uo(i, k - 1, jj))**2 &
                                        + (vo(i, k, jj) - vo(i, k - 1, jj))**2)
                           vshear = vshear + shear
                        end if
                     end do
                  
                     vshear = 1.e3*vshear/(zi(i, ktconr, jj) - zi(i, kbr, jj))
                     e1 = 1.591 - .639*vshear &
                          + .0953*(vshear**2) - .00496*(vshear**3)
                     edt(i, jj) = 1.-e1
                     val = .9
                     edt(i, jj) = min(edt(i, jj), val)
                     val = .0
                     edt(i, jj) = max(edt(i, jj), val)
                  end if
                  ktcon(i, jj) = ktconr
               end if
            end do
         end do
   !
   !c--- what would the change be, that a cloud with unit mass
   !c--- will do to the environment?
   !
   !> ## Calculate the tendencies of the state variables (per unit cloud base mass flux) and the cloud base mass flux.
   !> - Calculate the change in moist static energy, moisture mixing ratio, and horizontal winds per unit cloud base mass flux for all layers below cloud top from equations B.14 and B.15 from Grell (1993) \cite grell_1993, and for the cloud top from B.16 and
         !$acc parallel loop gang collapse(2) async(async_id)
         do jj = 1, jlistnum
            do k = 1, km
               !$acc loop vector
               do i = 1, myim(jj)
                  if (cnvflg(i, jj) .and. k <= kmax(i, jj)) then
                     dellah(i, k, jj) = 0.
                     dellaq(i, k, jj) = 0.
                     dellau(i, k, jj) = 0.
                     dellav(i, k, jj) = 0.
                  end if
               end do
            end do
         end do
   !
   !c--- changed due to subsidence and entrainment
   !
         !$acc parallel loop gang collapse(2) async(async_id) private(dp, &
         !$acc&         dz, dv1h, dv2h, dv3h, dv1q, dv2q, dv3q, tem, tem1)
         do jj = 1, jlistnum
            do k = 2, km1
               !$acc loop vector
               do i = 1, myim(jj)
                  if (cnvflg(i, jj)) then
                     if (k > kb(i, jj) .and. k < ktcon(i, jj)) then
                        dp = 1000.*del(i, k, jj)
                        dz = zi(i, k, jj) - zi(i, k - 1, jj)
   !
                        dv1h = heo(i, k, jj)
                        dv2h = .5*(heo(i, k, jj) + heo(i, k - 1, jj))
                        dv3h = heo(i, k - 1, jj)
                        dv1q = qo(i, k, jj)
                        dv2q = .5*(qo(i, k, jj) + qo(i, k - 1, jj))
                        dv3q = qo(i, k - 1, jj)
   !
                        tem = 0.5*(xlamue(i, k, jj) + xlamue(i, k - 1, jj))
                        tem1 = xlamud(i, jj)
   !cj
                        dellah(i, k, jj) = &
                                       (eta(i, k, jj)*dv1h - eta(i, k - 1, jj)*dv3h &
                                        - tem*eta(i, k - 1, jj)*dv2h*dz &
                                        + tem1*eta(i, k - 1, jj)*.5*(hcko(i, k, jj) &
                                        + hcko(i, k - 1, jj))*dz)*g/dp
   !cj
                        dellaq(i, k, jj) = &
                                       (eta(i, k, jj)*dv1q - eta(i, k - 1, jj)*dv3q &
                                        - tem*eta(i, k - 1, jj)*dv2q*dz &
                                        + tem1*eta(i, k - 1, jj)*.5*(qrcko(i, k, jj) &
                                        + qcko(i, k - 1, jj))*dz)*g/dp
   !cj
                        tem1 = eta(i, k, jj)*(uo(i, k, jj) - ucko(i, k, jj))
                        tem2 = eta(i, k - 1, jj)*(uo(i, k - 1, jj) - ucko(i, k - 1, jj))
                        dellau(i, k, jj) = (tem1 - tem2)*g/dp
   !cj
                        tem1 = eta(i, k, jj)*(vo(i, k, jj) - vcko(i, k, jj))
                        tem2 = eta(i, k - 1, jj)*(vo(i, k - 1, jj) - vcko(i, k - 1, jj))
                        dellav(i, k, jj) = (tem1 - tem2)*g/dp
   !cj
                     end if
                  end if
               end do
            end do
         end do
   !
   !c------- cloud top
   !
         !$acc parallel loop gang vector collapse(2) async(async_id) private( &
         !$acc&         indx, dp, dv1h, dv1q, tem, tfac, dz, val, sumx, umean, &
         !$acc&         po, dp, sigmagfm, scaldfunc, xmbmax, qesor)
         do jj = 1, jlistnum
            do i = 1, ix
               if (i .le. myim(jj)) then
                  if (cnvflg(i, jj)) then
                     indx = ktcon(i, jj)
                     dp = 1000.*del(i, indx, jj)
                     dv1h = heo(i, indx - 1, jj)
                     dellah(i, indx, jj) = eta(i, indx - 1, jj)* &
                                       (hcko(i, indx - 1, jj) - dv1h)*g/dp
                     dv1q = qo(i, indx - 1, jj)
                     dellaq(i, indx, jj) = eta(i, indx - 1, jj)* &
                                       (qcko(i, indx - 1, jj) - dv1q)*g/dp
                     dellau(i, indx, jj) = eta(i, indx - 1, jj)* &
                                       (ucko(i, indx - 1, jj) - uo(i, indx - 1, jj))*g/dp
                     dellav(i, indx, jj) = eta(i, indx - 1, jj)* &
                                       (vcko(i, indx - 1, jj) - vo(i, indx - 1, jj))*g/dp
      !
      !  cloud water
      !
                     dellal(i, indx, jj) = eta(i, indx - 1, jj)* &
                                       qlko_ktcon(i, jj)*g/dp
   !
   !
   !  compute convective turn-over time
   !
   !> - Following Bechtold et al. (2008) \cite bechtold_et_al_2008, calculate the convective turnover time using the mean updraft velocity (wc) and the cloud depth. It is also proportional to the grid size (gdx).
                     tem = zi(i, ktcon1(i, jj), jj) - zi(i, kbcon1(i, jj), jj)
                     dtconv = tem/wc(i, jj)
                     tfac = 1.+gdx(i, jj)/75000.
      ! reference from eq.3 in Zheng et al. 2016
      !byl          tfac = 1. + log(25000./gdx(i))
                     dtconv = tfac*dtconv
                     dtconv = max(dtconv, dtmin)
                     dtconv = max(dtconv, dt2)
                     dtconv = min(dtconv, dtmax)
   !
   !> - Calculate advective time scale (tauadv) using a mean cloud layer wind speed.
                     sumx = 0.
                     umean = 0.
                  
                     !$acc loop seq
                     do k = 2, km1
                        if (k >= kbcon1(i, jj) .and. k < ktcon1(i, jj)) then
                           dz = zi(i, k, jj) - zi(i, k - 1, jj)
                           tem = sqrt(u1(i, k, jj)*u1(i, k, jj) + v1(i, k, jj)*v1(i, k, jj))
                           umean = umean + tem*dz
                           sumx = sumx + dz
                        end if
                     end do

                     umean = umean/sumx
                     umean = max(umean, 1.)
                     tauadv = gdx(i, jj)/umean
   !
   !  compute cloud base mass flux as a function of the mean
   !     updraft velcoity
   !
   !> - From Han et al.'s (2017) \cite han_et_al_2017 equation 6, calculate cloud base mass flux as a function of the mean updraft velcoity.
   !!  As discussed in Han et al. (2017) \cite han_et_al_2017 , when dtconv is larger than tauadv, the convective mixing is not fully conducted before the cumulus cloud is advected out of the grid cell. In this case, therefore, the cloud base mass flux is fu
                     k = kbcon(i, jj)
                     po = .5*(pfld(i, k, jj) + pfld(i, k + 1, jj))
                     rho = po*100./(rd*to(i, k, jj))
                     tfac = tauadv/dtconv
                     tfac = min(tfac, 1.)
                     xmb(i, jj) = tfac*betaw*rho*wc(i, jj)
   !> - Calculate the maximum value of the cloud base mass flux using the CFL-criterion-based formula of Han and Pan (2011) \cite han_and_pan_2011, equation 7.
                     dp = 1000.*del(i, k, jj)
                     xmbmax = dp/(g*dt2)
   !
   !> - For scale-aware parameterization, the updraft fraction (sigmagfm) is first computed as a function of the lateral entrainment rate at cloud base (see Han et al.'s (2017) \cite han_et_al_2017 equation 4 and 5), following the study by Grell and Freitas
                     sigmagfm = wbar(i, jj)/wc(i, jj)
                     sigmagfm = tem1/garea(i, jj)
                     sigmagfm = max(sigmagfm, 0.001)
                     sigmagfm = min(sigmagfm, 0.999)
   !
   !> - Then, calculate the reduction factor (scaldfunc) of the vertical convective eddy transport of mass flux as a function of updraft fraction from the studies by Arakawa and Wu (2013) \cite arakawa_and_wu_2013 (also see Han et al.'s (2017) \cite han_et_a
                     scaldfunc = clear(i, jj)*(1.-sigma(i, jj))
                     scaldfunc = max(min(scaldfunc, 1.0), 0.)
                     xmb(i, jj) = xmb(i, jj)*scaldfunc
                     xmb(i, jj) = min(xmb(i, jj), xmbmax)
   !> ## For the "feedback control", calculate updated values of the state variables by multiplying the cloud base mass flux and the tendencies calculated per unit cloud base mass flux from the static control.
   !! - Recalculate saturation specific humidity.
   !
   !c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
                     !$acc loop seq
                     do k = 1, km
                        if (k <= kmax(i, jj)) then
                           qesor = 0.01*fpvs_gpu(t1(i, k, jj),c1xpvs,c2xpvs,tbpvs)      ! fpvs is in pa
         !byl            call qsatq_cwb(1,t1(i,k),pfld(i,k),qeso(i,k))
                           qesor = eps*qesor/(pfld(i, k, jj) + epsm1*qesor)
                           val = 1.e-8
                           qeso(i, k, jj) = max(qesor, val)
                        end if
                     end do
                  end if
               end if
            end do
         end do
   !c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   !> - Calculate the temperature tendency from the moist static energy and specific humidity tendencies.
   !> - Update the temperature, specific humidity, and horiztonal wind state variables by multiplying the cloud base mass flux-normalized tendencies by the cloud base mass flux.
   !> - Accumulate column-integrated tendencies.
         !$acc parallel loop gang vector collapse(2) async(async_id)
         do jj = 1, jlistnum
            do i = 1, ix
               if (i .le. myim(jj)) then
                  delhbar(i, jj) = 0.
                  delqbar(i, jj) = 0.
                  deltbar(i, jj) = 0.
                  delubar(i, jj) = 0.
                  delvbar(i, jj) = 0.
               end if
            end do
         end do
         
         !$acc parallel loop gang collapse(2) async(async_id) private(dellat, val)
         do jj = 1, jlistnum
            do k = 1, km
               !$acc loop vector
               do i = 1, myim(jj)
                  if (cnvflg(i, jj)) then
                     if (k > kb(i, jj) .and. k <= ktcon(i, jj)) then
                        dellat = (dellah(i, k, jj) - hvap*dellaq(i, k, jj))/cp
                        t1(i, k, jj) = t1(i, k, jj) + dellat*xmb(i, jj)*dt2
                        q1(i, k, jj) = q1(i, k, jj) + dellaq(i, k, jj)*xmb(i, jj)*dt2
   !             tem = 1./rcs(i)
   !             u1(i,k) = u1(i,k) + dellau(i,k) * xmb(i) * dt2 * tem
   !             v1(i,k) = v1(i,k) + dellav(i,k) * xmb(i) * dt2 * tem
                        u1(i, k, jj) = u1(i, k, jj) + dellau(i, k, jj)*xmb(i, jj)*dt2
                        v1(i, k, jj) = v1(i, k, jj) + dellav(i, k, jj)*xmb(i, jj)*dt2
                        dp = 1000.*del(i, k, jj)
                        !$acc atomic
                        delhbar(i, jj) = delhbar(i, jj) + dellah(i, k, jj)*xmb(i, jj)*dp/g
                        !$acc atomic
                        delqbar(i, jj) = delqbar(i, jj) + dellaq(i, k, jj)*xmb(i, jj)*dp/g
                        !$acc atomic
                        deltbar(i, jj) = deltbar(i, jj) + dellat*xmb(i, jj)*dp/g
                        !$acc atomic
                        delubar(i, jj) = delubar(i, jj) + dellau(i, k, jj)*xmb(i, jj)*dp/g
                        !$acc atomic
                        delvbar(i, jj) = delvbar(i, jj) + dellav(i, k, jj)*xmb(i, jj)*dp/g
   !
   !> - Recalculate saturation specific humidity using the updated temperature.
                        qeso(i, k, jj) = 0.01*fpvs_gpu(t1(i, k, jj),c1xpvs,c2xpvs,tbpvs)      ! fpvs is in pa
   !byl              call qsatq_cwb(1,t1(i,k),pfld(i,k),qeso(i,k))
                        qeso(i, k, jj) = eps*qeso(i, k, jj)/(pfld(i, k, jj) + epsm1*qeso(i, k, jj))
                        val = 1.e-8
                        qeso(i, k, jj) = max(qeso(i, k, jj), val)
                     end if
                  end if
               end do
            end do
         end do
   !
   !> - Add up column-integrated convective precipitation by multiplying the normalized value by the cloud base mass flux.
         !$acc parallel loop gang vector collapse(2) async(async_id) private( &
         !$acc&         evef, dp, tem, tem1, flgr, rntot, delqev, delq2, deltv, &
         !$acc&         delq, qevap, ktconr, kbr, rnr, xmbr, edtr, islimskr, &
         !$acc&         qcond, delqbarr, deltbarr)
         do jj = 1, jlistnum
            do i = 1, ix
               if (i .le. myim(jj)) then
                  qcond = 0.
                  if (cnvflg(i, jj)) then
                     rntot = 0.
                     delqev = 0.
                     delq2 = 0.
                     flgr = .true.
                     ktconr = ktcon(i, jj)
                     kbr = kb(i, jj)
                     rnr = rn(i, jj)
                     xmbr = xmb(i, jj)
                     edtr = edt(i, jj)
                     islimskr = islimsk(i, jj)
                     delqbarr = delqbar(i, jj)
                     deltbarr = deltbar(i, jj)

                     !$acc loop seq
                     do k = km1, 1, -1
                        if (k < ktconr .and. k > kbr) then
                           rntot = rntot + pwo(i, k, jj)*xmbr*.001*dt2
                        end if
                     end do
   !
   ! evaporating rain
   !
   !> - Determine the evaporation of the convective precipitation and update the integrated convective precipitation.
   !> - Update state temperature and moisture to account for evaporation of convective precipitation.
   !> - Update column-integrated tendencies to account for evaporation of convective precipitation.
                     delqbarr = 0.
                     delqbarr = 0.
                     !$acc loop seq
                     do k = km1, 1, -1
                        if (k <= kmax(i, jj)) then
                           deltv = 0.
                           delq = 0.
                           qevap = 0.
                           if (k < ktconr .and. k > kbr) then
                              rnr = rnr + pwo(i, k, jj)*xmbr*.001*dt2
                           end if
                           if (flgr .and. k < ktconr) then
                              evef = edtr*evfact
                              if (islimskr == 1) evef = edtr*evfactl
         !             if(islimsk(i) == 1) evef=.07
         !             if(islimsk(i) == 1) evef = 0.
                              qcond = evef*(q1(i, k, jj) - qeso(i, k, jj)) &
                                         /(1.+el2orc*qeso(i, k, jj)/t1(i, k, jj)**2)
                              dp = 1000.*del(i, k, jj)
                              if (rnr > 0. .and. qcond < 0.) then
                                 qevap = -qcond*(1.-exp(-.32*sqrt(dt2*rnr)))
                                 qevap = min(qevap, rnr*1000.*g/dp)
                                 delq2 = delqev + .001*qevap*dp/g
                              end if
                              if (rnr > 0. .and. qcond < 0. .and. &
                                  delq2 > rntot) then
                                 qevap = 1000.*g*(rntot - delqev)/dp
                                 flgr = .false.
                              end if
                              if (rnr > 0. .and. qevap > 0.) then
                                 tem = .001*dp/g
                                 tem1 = qevap*tem
                                 if (tem1 > rnr) then
                                    qevap = rnr/tem
                                    rnr = 0.
                                 else
                                    rnr = rnr - tem1
                                 end if
                                 q1(i, k, jj) = q1(i, k, jj) + qevap
                                 t1(i, k, jj) = t1(i, k, jj) - elocp*qevap
                                 deltv = -elocp*qevap/dt2
                                 delq = +qevap/dt2
                                 delqev = delqev + .001*dp*qevap/g
                              end if
                              delqbarr = delqbarr + delq*dp/g
                              deltbarr = deltbarr + deltv*dp/g
                           end if
                        end if
                     end do
                     delqbar(i, jj) = delqbarr
                     deltbar(i, jj) = deltbarr
                     
   !cj
   !     do i = 1, im
   !     if(me == 31 .and. cnvflg(i)) then
   !     if(cnvflg(i)) then
   !       print *, ' shallow delhbar, delqbar, deltbar = ',
   !    &             delhbar(i),hvap*delqbar(i),cp*deltbar(i)
   !       print *, ' shallow delubar, delvbar = ',delubar(i),delvbar(i)
   !       print *, ' precip =', hvap*rn(i)*1000./dt2
   !       print*,'pdif= ',pfld(i,kbcon(i))-pfld(i,ktcon(i))
   !     endif
   !     enddo
   !cj
                     if (rnr < 0. .or. .not. flgr) rnr = 0.
                     if (rnr < 1.e-13) rnr = 0.
                     rn(i, jj) = rnr
                     ktop(i, jj) = ktconr
                     kbot(i, jj) = kbcon(i, jj)
                     kcnv(i, jj) = 2
                  end if
               end if
            end do
         end do
   !
   !  convective cloud water
   !
   !> - Calculate shallow convective cloud water.
         !$acc parallel loop gang collapse(2) async(async_id) private(tem, tem1)
         do jj = 1, jlistnum
            do k = 1, km
               !$acc loop vector
               do i = 1, myim(jj)
   !>  - Initialize convective cloud water and cloud cover to zero.
   ! hchuang code change
                  cnvw(i, k, jj) = 0.
                  cnvc(i, k, jj) = 0.
   !>  - Initialize updraft mass fluxes to zero.
                  ud_mf(i, k, jj) = 0.
                  dt_mf(i, k, jj) = 0.
                  if (cnvflg(i, jj)) then
                     if (k >= kbcon(i, jj) .and. k < ktcon(i, jj)) then
                        cnvw(i, k, jj) = cnvwt(i, k, jj)*xmb(i, jj)*dt2

   !
   !  convective cloud cover
   !
   !> - Calculate convective cloud cover, which is used when pdf-based cloud fraction is used (i.e., pdfcld=.true.).
                        tem = 0.04*log(1.+675.*eta(i, k, jj)*xmb(i, jj))
                        tem = min(tem, 0.2)
                        cnvc(i, k, jj) = max(tem, 0.0)
                     end if
   !
   ! hchuang code change
   !
   !> - Calculate and retain the updraft mass flux for dust transport by cumulus convection.
   !
   !> - Calculate the updraft convective mass flux.
                     if (k >= kb(i, jj) .and. k < ktop(i, jj)) then
                        ud_mf(i, k, jj) = eta(i, k, jj)*xmb(i, jj)*dt2
                     end if
   !
   !  cloud water
   !
   !> - Separate detrained cloud water into liquid and ice species as a function of temperature only.
   !
                     if (ncloud > 0) then
   !           if (k > kb(i) .and. k <= ktcon(i)) then
                        if (k >= kbcon(i, jj) .and. k <= ktcon(i, jj)) then
                           tem = dellal(i, k, jj)*xmb(i, jj)*dt2
                           tem1 = max(0.0, min(1.0, (tcr - t1(i, k, jj))*tcrf))
                           if (qi(i, k, jj) > -999.0) then
                              qi(i, k, jj) = qi(i, k, jj) + tem*tem1            ! ice
                              ql(i, k, jj) = ql(i, k, jj) + tem*(1.0 - tem1)       ! water
                           else
                              ql(i, k, jj) = ql(i, k, jj) + tem
                           end if
                        end if
                     end if
                  end if
               end do
            end do
         end do
   !
   !> - save the updraft convective mass flux at cloud top.
         !$acc parallel loop gang vector collapse(2) async(async_id) private(k)
         do jj = 1, jlistnum
            do i = 1, ix
               if (i .le. myim(jj)) then
                  if (cnvflg(i, jj)) then
                     k = ktop(i, jj) - 1
                     dt_mf(i, k, jj) = ud_mf(i, k, jj)
                  end if
               end if
            end do
         end do
         !$acc end data
!!
         return
      end
