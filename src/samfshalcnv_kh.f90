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
      subroutine samfshalcnv_kh(im,ix,km,delt,delp,prslp,psp,phil,ql, &
           qi,q1,t1,u1,v1,rn,kbot,ktop,kcnv,islimsk,garea, &
           dot,ncloud,hpbl,cnvw,cnvc)
!byl           dot,ncloud,hpbl,ud_mf,dt_mf,cnvw,cnvc,
!    &     dot,ncloud,hpbl,ud_mf,dt_mf,cnvw,cnvc,me) &
!byl           clam,c0s,c1,pgcon,asolfac)
!
      use machine , only : kind_phys
!byl      use funcphys , only : fpvs
      use physcons, grav => con_g, cp => con_cp, hvap => con_hvap &
      ,             rv => con_rv, fv => con_fvirt, t0c => con_t0c &
      ,             rd => con_rd, cvap => con_cvap, cliq => con_cliq &
      ,             eps => con_eps, epsm1 => con_epsm1
      implicit none
!
      real fpvs
!
      integer            im, ix,  km, ncloud, &
                         kbot(im), ktop(im), kcnv(im)
!    &,                  me
      real(kind=kind_phys) delt
      real(kind=kind_phys) psp(im),    delp(ix,km), prslp(ix,km)
      real(kind=kind_phys) ps(im),     del(ix,km),  prsl(ix,km), &
                           ql(ix,km),  q1(ix,km),   t1(ix,km),   &
                           u1(ix,km),  v1(ix,km),                &
!    &                     u1(ix,km),  v1(ix,km),   rcs(im), &
                           rn(im),     garea(im),                &
                           dot(ix,km), phil(ix,km), hpbl(im),    &
                           cnvw(ix,km),cnvc(ix,km)               &
! hchuang code change mass flux output &
      ,                    ud_mf(im,km),dt_mf(im,km),qi(ix,km)
!
      integer              i,j,indx, k, kk, km1, n
      integer              kpbl(im)
      integer, dimension(im), intent(in) :: islimsk
!
      real(kind=kind_phys) dellat,  delta, &
                           c0l,     c0s,     d0, &
                           c1,      asolfac, &
                           desdt,   dp, &
                           dq,      dqsdp,   dqsdt,   dt, &
                           dt2,     dtmax,   dtmin,   dxcrt, &
                           dv1h,    dv2h,    dv3h, &
                           dv1q,    dv2q,    dv3q, &
                           dz,      dz1,     e1,      clam, &
                           el2orc,  elocp,   aafac,   cm, &
                           es,      etah,    h1, &
                           evef,    evfact,  evfactl, fact1, &
                           fact2,   factor,  dthk, &
                           g,       gamma,   pprime,  betaw, &
                           qlk,     qrch,    qs, &
                           rfact,   shear,   tfac, &
                           val,     val1,    val2, &
                           w1,      w1l,     w1s,     w2, &
                           w2l,     w2s,     w3,      w3l, &
                           w3s,     w4,      w4l,     w4s, &
                           rho,     tem,     tem1,    tem2, &
                           ptem,    ptem1, &
                           pgcon
!
      integer              kb(im), kbcon(im), kbcon1(im), &
                           ktcon(im), ktcon1(im), ktconn(im), &
                           kbm(im), kmax(im)
!
      real(kind=kind_phys) aa1(im),     cina(im), &
                           umean(im),  tauadv(im),  gdx(im), &
                           delhbar(im), delq(im),   delq2(im), &
                           delqbar(im), delqev(im), deltbar(im), &
                           deltv(im),   dtconv(im), edt(im), &
                           pdot(im),    po(im,km), &
                           qcond(im),   qevap(im),  hmax(im), &
                           rntot(im),   vshear(im), &
                           xlamud(im),  xmb(im),    xmbmax(im), &
                           delubar(im), delvbar(im)
!
      real(kind=kind_phys) c0(im)
!
      real(kind=kind_phys) crtlamd
!
      real(kind=kind_phys) cinpcr,  cinpcrmx,  cinpcrmn, &
                           cinacr,  cinacrmx,  cinacrmn
!
!  parameters for updraft velocity calculation
      real(kind=kind_phys) bet1,    cd1,     f1,      gam1, &
                           bb1,     bb2
!    &                     bb1,     bb2,     wucb
!cc
!  physical parameters
      parameter(g=grav,asolfac=0.958)
!byl      parameter(g=grav)
      parameter(elocp=hvap/cp, &
                el2orc=hvap*hvap/(rv*cp))
      parameter(c0s=0.002,c1=5.e-4,d0=.01)
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
      parameter(cm=1.0,delta=fv)
      parameter(fact1=(cvap-cliq)/rv,fact2=hvap/rv-fact1*t0c)
      parameter(dthk=25.)
      parameter(cinpcrmn=120.)
!      parameter(cinpcrmx=180.,cinpcrmn=120.)
!     parameter(cinacrmx=-120.,cinacrmn=-120.)
      parameter(cinacrmx=-120.,cinacrmn=-80.)
      parameter(crtlamd=3.e-4)
      parameter(dtmax=10800.,dtmin=600.)
      parameter(bet1=1.875,cd1=.506,f1=2.0,gam1=.5)
      parameter(betaw=.03,dxcrt=15.e3)
      parameter(h1=0.33333333)
!  local variables and arrays
      real(kind=kind_phys) pfld(im,km),    to(im,km),     qo(im,km), &
                           uo(im,km),      vo(im,km),     qeso(im,km)
!  for updraft velocity calculation
      real(kind=kind_phys) wu2(im,km),     buo(im,km),    drag(im,km)
      real(kind=kind_phys) wc(im),         scaldfunc(im), sigmagfm(im)
!
!  cloud water
!     real(kind=kind_phys) qlko_ktcon(im), dellal(im,km), tvo(im,km),
      real(kind=kind_phys) qlko_ktcon(im), dellal(im,km), &
                           dbyo(im,km),    zo(im,km),     xlamue(im,km), &
                           heo(im,km),     heso(im,km),   frh(im,km),&
                           dellah(im,km),  dellaq(im,km), &
                           dellau(im,km),  dellav(im,km), hcko(im,km), &
                           ucko(im,km),    vcko(im,km),   qcko(im,km), &
                           qrcko(im,km),   eta(im,km), &
                           zi(im,km),      pwo(im,km),    c0t(im,km), &
                           sumx(im),       tx1(im),       cnvwt(im,km)
!
      logical totflg, cnvflg(im), flg(im)
!
      real(kind=kind_phys) tf, tcr, tcrf
      parameter (tf=233.16, tcr=263.16, tcrf=1.0/(tcr-tf))
!<---for scale-aware parameterization (Kwon and Hong 2017)
      real(kind=kind_phys) cinpcri(im), frh_sum(im), wbar(im),      &
                           clear(im),sigma(im)
      real(kind=kind_phys) po1(im,km)
      real(kind=kind_phys) sigma_con, pi, dx1km, dx5km, dx250m
      parameter (pi=3.14159)
      parameter (dx1km = 1000., dx5km = 5000., dx250m = 250.)
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
      do i=1,im
        ps(i) = psp(i)
      enddo
      do k=1,km
        do i=1,im
          prsl(i,k) = prslp(i,k)
          del(i,k)  = delp(i,k)
        enddo
      enddo
!************************************************************************
!
      km1 = km - 1
!
!  initialize arrays
!
!>  - Initialize column-integrated and other single-value-per-column variable arrays.
      do i=1,im
        cnvflg(i) = .true.
        if(kcnv(i) == 1) cnvflg(i) = .false.
        if(cnvflg(i)) then
          kbot(i)=km+1
          ktop(i)=0
        endif
        rn(i)=0.
        kbcon(i)=km
        ktcon(i)=1
        ktconn(i)=1
        kb(i)=km
        pdot(i) = 0.
        qlko_ktcon(i) = 0.
        edt(i)  = 0.
        aa1(i)  = 0.
        cina(i) = 0.
        vshear(i) = 0.
        gdx(i) = sqrt(garea(i))
      enddo

      cinpcrmx=240.
      do i = 1,im
       sigma(i) = 0.
       frh_sum(i) = 0.
       sigmagfm(i) = 0.
      end do
!!
!>  - Return to the calling routine if deep convection is present or the surface buoyancy flux is negative.
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!>  - determine aerosol-aware rain conversion parameter over land
      do i=1,im
        if(islimsk(i) == 1) then
           c0(i) = c0s*asolfac
        else
           c0(i) = c0s
        endif
      enddo
!
!>  - determine rain conversion parameter above the freezing level which exponentially decreases with decreasing temperature from Han et al.'s (2017) \cite han_et_al_2017 equation 8.
      do k = 1, km
        do i = 1, im
          if(t1(i,k) > 273.16) then
            c0t(i,k) = c0(i)
          else
            tem = d0 * (t1(i,k) - 273.16)
            tem1 = exp(tem)
            c0t(i,k) = c0(i) * tem1
          endif
        enddo
      enddo
!
!>  - Initialize convective cloud water and cloud cover to zero.
      do k = 1, km
        do i = 1, im
          cnvw(i,k) = 0.
          cnvc(i,k) = 0.
        enddo
      enddo
! hchuang code change
!>  - Initialize updraft mass fluxes to zero.
      do k = 1, km
        do i = 1, im
          ud_mf(i,k) = 0.
          dt_mf(i,k) = 0.
        enddo
      enddo
!
      dt2   = delt
!
!  model tunable parameters are all here
      clam    = .3
      aafac   = .1
!     evef    = 0.07
      evfact  = 0.3
      evfactl = 0.3
!
!     pgcon   = 0.7     ! Gregory et al. (1997, QJRMS)
      pgcon   = 0.55    ! Zhang & Wu (2003,JAS)
!
      w1l     = -8.e-3
      w2l     = -4.e-2
      w3l     = -5.e-3
      w4l     = -5.e-4
      w1s     = -2.e-4
      w2s     = -2.e-3
      w3s     = -1.e-3
      w4s     = -2.e-5
!
!  define top layer for search of the downdraft originating layer
!  and the maximum thetae for updraft
!
!>  - Determine maximum indices for the parcel starting point (kbm) and cloud top (kmax).
      do i=1,im
        kbm(i)   = km
        kmax(i)  = km
        tx1(i)   = 1.0 / ps(i)
      enddo
!
      do k = 1, km
        do i=1,im
          if (prsl(i,k)*tx1(i) > 0.70) kbm(i)   = k + 1
          if (prsl(i,k)*tx1(i) > 0.60) kmax(i)  = k + 1
        enddo
      enddo
      do i=1,im
        kbm(i)   = min(kbm(i),kmax(i))
      enddo
!
!  hydrostatic height assume zero terr and compute
!  updraft entrainment rate as an inverse function of height
!
!>  - Calculate hydrostatic height at layer centers assuming a flat surface (no terrain) from the geopotential.
      do k = 1, km
        do i=1,im
          zo(i,k) = phil(i,k) / g
        enddo
      enddo
!>  - Calculate interface height and the entrainment rate as an inverse function of height.
      do k = 1, km1
        do i=1,im
          zi(i,k) = 0.5*(zo(i,k)+zo(i,k+1))
          xlamue(i,k) = clam / zi(i,k)
        enddo
      enddo
      do i=1,im
        xlamue(i,km) = xlamue(i,km1)
      enddo
!
!  pbl height
!
!>  - Find the index for the PBL top using the PBL height; enforce that it is lower than the maximum parcel starting level.
      do i=1,im
        flg(i) = cnvflg(i)
        kpbl(i)= 1
      enddo
      do k = 2, km1
        do i=1,im
          if (flg(i) .and. zo(i,k) <= hpbl(i)) then
            kpbl(i) = k
          else
            flg(i) = .false.
          endif
        enddo
      enddo
      do i=1,im
        kpbl(i)= min(kpbl(i),kbm(i))
      enddo
!
!c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   convert surface pressure to mb from cb
!
!>  - Convert prsl from centibar to millibar, set normalized mass flux to 1, cloud properties to 0, and save model state variables (after advection/turbulence).
      do k = 1, km
        do i = 1, im
          if (cnvflg(i) .and. k <= kmax(i)) then
            pfld(i,k) = prsl(i,k) * 10.0
            eta(i,k)  = 1.
            frh(i,k)  = 0.
            hcko(i,k) = 0.
            qcko(i,k) = 0.
            qrcko(i,k)= 0.
            ucko(i,k) = 0.
            vcko(i,k) = 0.
            dbyo(i,k) = 0.
            pwo(i,k)  = 0.
            dellal(i,k) = 0.
            to(i,k)   = t1(i,k)
            qo(i,k)   = q1(i,k)
            uo(i,k)   = u1(i,k)
            vo(i,k)   = v1(i,k)
!           uo(i,k)   = u1(i,k) * rcs(i)
!           vo(i,k)   = v1(i,k) * rcs(i)
            wu2(i,k)  = 0.
            buo(i,k)  = 0.
            drag(i,k) = 0.
            cnvwt(i,k) = 0.
          endif
        enddo
      enddo
!>  - Calculate saturation specific humidity and enforce minimum moisture values.
      do k = 1, km
        do i=1,im
          if (cnvflg(i) .and. k <= kmax(i)) then
            qeso(i,k) = 0.01 * fpvs(to(i,k))      ! fpvs is in pa
!byl            call qsatq_cwb(1,to(i,k),pfld(i,k),qeso(i,k))
            qeso(i,k) = eps * qeso(i,k) / (pfld(i,k) + epsm1*qeso(i,k))
            val1      =             1.e-8
            qeso(i,k) = max(qeso(i,k), val1)
            val2      =           1.e-10
            qo(i,k)   = max(qo(i,k), val2 )
!           qo(i,k)   = min(qo(i,k),qeso(i,k))
!           tvo(i,k)  = to(i,k) + delta * to(i,k) * qo(i,k)
          endif
        enddo
      enddo
!
!  compute moist static energy
!
!>  - Calculate moist static energy (heo) and saturation moist static energy (heso).
      do k = 1, km
        do i=1,im
          if (cnvflg(i) .and. k <= kmax(i)) then
!           tem       = g * zo(i,k) + cp * to(i,k)
            tem       = phil(i,k) + cp * to(i,k)
            heo(i,k)  = tem  + hvap * qo(i,k)
            heso(i,k) = tem  + hvap * qeso(i,k)
!           heo(i,k)  = min(heo(i,k),heso(i,k))
          endif
        enddo
      enddo
!
!  determine level with largest moist static energy within pbl
!  this is the level where updraft starts
!
!> ## Perform calculations related to the updraft of the entraining/detraining cloud model ("static control").
!> - Search in the PBL for the level of maximum moist static energy to start the ascending parcel.
      do i=1,im
         if (cnvflg(i)) then
            hmax(i) = heo(i,1)
            kb(i) = 1
         endif
      enddo
      do k = 2, km
        do i=1,im
          if (cnvflg(i) .and. k <= kpbl(i)) then
            if(heo(i,k) > hmax(i)) then
              kb(i)   = k
              hmax(i) = heo(i,k)
            endif
          endif
        enddo
      enddo
!
!> - Calculate the temperature, water vapor mixing ratio, and pressure at interface levels.
      do k = 1, km1
        do i=1,im
          if (cnvflg(i) .and. k <= kmax(i)-1) then
            dz      = .5 * (zo(i,k+1) - zo(i,k))
            dp      = .5 * (pfld(i,k+1) - pfld(i,k))
            es      = 0.01 * fpvs(to(i,k+1))      ! fpvs is in pa
!byl            call qsatq_cwb(1,to(i,k+1),pfld(i,k+1),es)
            pprime  = pfld(i,k+1) + epsm1 * es
            qs      = eps * es / pprime
            dqsdp   = - qs / pprime
            desdt   = es * (fact1 / to(i,k+1) + fact2 / (to(i,k+1)**2))
            dqsdt   = qs * pfld(i,k+1) * desdt / (es * pprime)
            gamma   = el2orc * qeso(i,k+1) / (to(i,k+1)**2)
            dt      = (g * dz + hvap * dqsdp * dp) / (cp * (1. + gamma))
            dq      = dqsdt * dt + dqsdp * dp
            to(i,k) = to(i,k+1) + dt
            qo(i,k) = qo(i,k+1) + dq
            po(i,k) = .5 * (pfld(i,k) + pfld(i,k+1))
          endif
        enddo
      enddo
!
!> - Recalculate saturation specific humidity, moist static energy, saturation moist static energy, and horizontal momentum on interface levels. Enforce minimum specific humidity.
      do k = 1, km1
        do i=1,im
          if (cnvflg(i) .and. k <= kmax(i)-1) then
            qeso(i,k) = 0.01 * fpvs(to(i,k))      ! fpvs is in pa
!byl            call qsatq_cwb(1,to(i,k),po(i,k),qeso(i,k))
            qeso(i,k) = eps * qeso(i,k) / (po(i,k) + epsm1*qeso(i,k))
            val1      =             1.e-8
            qeso(i,k) = max(qeso(i,k), val1)
            val2      =           1.e-10
            qo(i,k)   = max(qo(i,k), val2 )
!           qo(i,k)   = min(qo(i,k),qeso(i,k))
            frh(i,k)  = 1. - min(qo(i,k)/qeso(i,k), 1.)
            heo(i,k)  = .5 * g * (zo(i,k) + zo(i,k+1)) + &
                        cp * to(i,k) + hvap * qo(i,k)
            heso(i,k) = .5 * g * (zo(i,k) + zo(i,k+1)) + &
                        cp * to(i,k) + hvap * qeso(i,k)
            uo(i,k)   = .5 * (uo(i,k) + uo(i,k+1))
            vo(i,k)   = .5 * (vo(i,k) + vo(i,k+1))
          endif
        enddo
      enddo
!
!  look for the level of free convection as cloud base
!
!> - Search below the index "kbm" for the level of free convection (LFC) where the condition \f$h_b > h^*\f$ is first met, where \f$h_b, h^*\f$ are the state moist static energy at the parcel's starting level and saturation moist static energy, respective
      do i=1,im
        flg(i)   = cnvflg(i)
        if(flg(i)) kbcon(i) = kmax(i)
      enddo
      do k = 2, km1
        do i=1,im
          if (flg(i) .and. k < kbm(i)) then
            if(k > kb(i) .and. heo(i,kb(i)) > heso(i,k)) then
              kbcon(i) = k
              flg(i)   = .false.
            endif
          endif
        enddo
      enddo
!
      do i=1,im
        if(cnvflg(i)) then
          if(kbcon(i) == kmax(i)) cnvflg(i) = .false.
        endif
      enddo
!!
!> - If no LFC, return to the calling routine without modifying state variables.
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!> - Determine the vertical pressure velocity at the LFC. After Han and Pan (2011) \cite han_and_pan_2011 , determine the maximum pressure thickness between a parcel's starting level and the LFC. If a parcel doesn't reach the LFC within the critical thick
      do i=1,im
        if(cnvflg(i)) then
          pdot(i)  = 10.* dot(i,kbcon(i))
!          pdot(i)  = 0.01 * dot(i,kbcon(i)) ! Now dot is in Pa/s
        endif
      enddo

      do k = 1, km1
       do i = 1, im
         if (cnvflg(i)) then
          if (k.ge.kb(i).and.k.le.kbcon(i)) then
           frh_sum(i) = frh_sum(i) + (1-frh(i,k))
          endif
         endif
       enddo
      enddo
!
!   turn off convection if pressure depth between parcel source level
!      and cloud base is larger than a critical value, cinpcr
!
      do i=1,im
        if(cnvflg(i)) then
          sigma_con = tan(0.4*pi)/(dx5km-dx1km)                     !7.7 e-4 m-1
          sigma(i)  = (1.-1./pi*(atan(sigma_con*(gdx(i)-dx5km))+pi/2.)) !1(1km),0.1(10km)
            if (gdx(i).lt.dx5km) then
             sigma(i) = min(sigma(i)-0.01684*gdx(i)/1000.+0.0842,1.0)
            endif
            cinpcr = cinpcrmn + 0.5*(cinpcrmx-cinpcrmn) * (1.-sigma(i))
            cinpcri(i) = cinpcr * frh_sum(i)/(kbcon(i)-kb(i)+1)
            tem1 = pfld(i,kb(i)) - pfld(i,kbcon(i))
            if(tem1 > cinpcri(i)) then
               cnvflg(i) = .false.
            endif
        endif
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!
!  specify the detrainment rate for the updrafts
!
!> - The updraft detrainment rate is set constant and equal to the entrainment rate at cloud base.
      do i = 1, im
        if(cnvflg(i)) then
          xlamud(i) = xlamue(i,kbcon(i))
!         xlamud(i) = crtlamd
        endif
      enddo
!
!  determine updraft mass flux for the subcloud layers
!
!> - Calculate the normalized mass flux for subcloud and in-cloud layers according to Pan and Wu (1995) \cite pan_and_wu_1995 equation 1:
!!  \f[
!!  \frac{1}{\eta}\frac{\partial \eta}{\partial z} = \lambda_e - \lambda_d
!!  \f]
!!  where \f$\eta\f$ is the normalized mass flux, \f$\lambda_e\f$ is the entrainment rate and \f$\lambda_d\f$ is the detrainment rate. The normalized mass flux increases upward below the cloud base and decreases upward above.
      do k = km1, 1, -1
        do i = 1, im
          if (cnvflg(i)) then
            if(k < kbcon(i) .and. k >= kb(i)) then
              dz       = zi(i,k+1) - zi(i,k)
              ptem     = 0.5*(xlamue(i,k)+xlamue(i,k+1))-xlamud(i)
              eta(i,k) = eta(i,k+1) / (1. + ptem * dz)
            endif
          endif
        enddo
      enddo
!
!  compute mass flux above cloud base
!
      do i = 1, im
        flg(i) = cnvflg(i)
      enddo
      do k = 2, km1
        do i = 1, im
         if(flg(i))then
           if(k > kbcon(i) .and. k < kmax(i)) then
              dz       = zi(i,k) - zi(i,k-1)
              ptem     = 0.5*(xlamue(i,k)+xlamue(i,k-1))-xlamud(i)
              eta(i,k) = eta(i,k-1) * (1 + ptem * dz)
              if(eta(i,k) <= 0.) then
                kmax(i) = k
                ktconn(i) = k
                kbm(i) = min(kbm(i),kmax(i))
                flg(i) = .false.
              endif
           endif
         endif
        enddo
      enddo
!
!  compute updraft cloud property
!
!> - Set cloud properties equal to the state variables at updraft starting level (kb).
      do i = 1, im
        if(cnvflg(i)) then
          indx         = kb(i)
          hcko(i,indx) = heo(i,indx)
          ucko(i,indx) = uo(i,indx)
          vcko(i,indx) = vo(i,indx)
        endif
      enddo
!
!  cm is an enhancement factor in entrainment rates for momentum
!
!> - Calculate the cloud properties as a parcel ascends, modified by entrainment and detrainment. Discretization follows Appendix B of Grell (1993) \cite grell_1993 . Following Han and Pan (2006) \cite han_and_pan_2006, the convective momentum transport i
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k > kb(i) .and. k < kmax(i)) then
              dz   = zi(i,k) - zi(i,k-1)
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              hcko(i,k) = ((1.-tem1)*hcko(i,k-1)+tem*0.5* &
                           (heo(i,k)+heo(i,k-1)))/factor
              dbyo(i,k) = hcko(i,k) - heso(i,k)
!
              tem  = 0.5 * cm * tem
              factor = 1. + tem
              ptem = tem + pgcon
              ptem1= tem - pgcon
              ucko(i,k) = ((1.-tem)*ucko(i,k-1)+ptem*uo(i,k) &
                           +ptem1*uo(i,k-1))/factor
              vcko(i,k) = ((1.-tem)*vcko(i,k-1)+ptem*vo(i,k) &
                           +ptem1*vo(i,k-1))/factor
            endif
          endif
        enddo
      enddo
!
!   taking account into convection inhibition due to existence of
!    dry layers below cloud base
!
!> - With entrainment, recalculate the LFC as the first level where buoyancy is positive. The difference in pressure levels between LFCs calculated with/without entrainment must be less than a threshold (currently 25 hPa). Otherwise, convection is inhibit
      do i=1,im
        flg(i) = cnvflg(i)
        kbcon1(i) = kmax(i)
      enddo
      do k = 2, km1
      do i=1,im
        if (flg(i) .and. k < kbm(i)) then
          if(k >= kbcon(i) .and. dbyo(i,k) > 0.) then
            kbcon1(i) = k
            flg(i)    = .false.
          endif
        endif
      enddo
      enddo
      do i=1,im
        if(cnvflg(i)) then
          if(kbcon1(i) == kmax(i)) cnvflg(i) = .false.
        endif
      enddo
      do i=1,im
        if(cnvflg(i)) then
          tem = pfld(i,kbcon(i)) - pfld(i,kbcon1(i))
          if(tem > dthk) then
             cnvflg(i) = .false.
          endif
        endif
      enddo
!!
      totflg = .true.
      do i = 1, im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!
!  calculate convective inhibition
!
!> - Calculate additional trigger condition of the convective inhibition (CIN) according to Han et al.'s (2017) \cite han_et_al_2017 equation 13.
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k > kb(i) .and. k < kbcon1(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              rfact =  1. + delta * cp * gamma &
                       * to(i,k) / hvap
              cina(i) = cina(i) +                 &
!    &                 dz1 * eta(i,k) * (g / (cp * to(i,k))) &
                       dz1 * (g / (cp * to(i,k))) &
                       * dbyo(i,k) / (1. + gamma) &
                       * rfact
              val = 0.
              cina(i) = cina(i) +                      &
!    &                 dz1 * eta(i,k) * g * delta * &
                       dz1 * g * delta *               &
                       max(val,(qeso(i,k) - qo(i,k)))
            endif
          endif
        enddo
      enddo
!> - Turn off convection if the CIN is less than a critical value (cinacr) which is inversely proportional to the large-scale vertical velocity.
      do i = 1, im
        if(cnvflg(i)) then
!
          if(islimsk(i) == 1) then
            w1 = w1l
            w2 = w2l
            w3 = w3l
            w4 = w4l
          else
            w1 = w1s
            w2 = w2s
            w3 = w3s
            w4 = w4s
          endif
          if(pdot(i) <= w4) then
            tem = (pdot(i) - w4) / (w3 - w4)
          elseif(pdot(i) >= -w4) then
            tem = - (pdot(i) + w4) / (w4 - w3)
          else
            tem = 0.
          endif

          val1    =            -1.
          tem = max(tem,val1)
          val2    =             1.
          tem = min(tem,val2)
          tem = 1. - tem
          tem1= .5*(cinacrmx-cinacrmn)
          cinacr = cinacrmx - tem * tem1
!
!         cinacr = cinacrmx
          if(cina(i) < cinacr) cnvflg(i) = .false.
        endif
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!
!  determine first guess cloud top as the level of zero buoyancy
!    limited to the level of P/Ps=0.7
!
!> - Calculate the cloud top as the first level where parcel buoyancy becomes negative; the maximum possible value is at \f$p=0.7p_{sfc}\f$.
      do i = 1, im
        flg(i) = cnvflg(i)
        if(flg(i)) ktcon(i) = kbm(i)
      enddo
      do k = 2, km1
      do i=1,im
        if (flg(i) .and. k < kbm(i)) then
          if(k > kbcon1(i) .and. dbyo(i,k) < 0.) then
             ktcon(i) = k
             flg(i)   = .false.
          endif
        endif
      enddo
      enddo
!
!  specify upper limit of mass flux at cloud base
!
!> - Calculate the maximum value of the cloud base mass flux using the CFL-criterion-based formula of Han and Pan (2011) \cite han_and_pan_2011, equation 7.
      do i = 1, im
        if(cnvflg(i)) then
!         xmbmax(i) = .1
!
          k = kbcon(i)
          dp = 1000. * del(i,k)
          xmbmax(i) = dp / (g * dt2)
!
!         tem = dp / (g * dt2)
!         xmbmax(i) = min(tem, xmbmax(i))
        endif
      enddo
!
!  compute cloud moisture property and precipitation
!
!> - Set cloud moisture property equal to the enviromental moisture at updraft starting level (kb).
      do i = 1, im
        if (cnvflg(i)) then
          aa1(i) = 0.
          qcko(i,kb(i)) = qo(i,kb(i))
          qrcko(i,kb(i)) = qo(i,kb(i))
        endif
      enddo
!> - Calculate the moisture content of the entraining/detraining parcel (qcko) and the value it would have if just saturated (qrch), according to equation A.14 in Grell (1993) \cite grell_1993 . Their difference is the amount of convective cloud water (ql
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k > kb(i) .and. k < ktcon(i)) then
              dz    = zi(i,k) - zi(i,k-1)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              qrch = qeso(i,k) &
                   + gamma * dbyo(i,k) / (hvap * (1. + gamma))
!cj
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              qcko(i,k) = ((1.-tem1)*qcko(i,k-1)+tem*0.5* &
                           (qo(i,k)+qo(i,k-1)))/factor
              qrcko(i,k) = qcko(i,k)
!cj
              dq = eta(i,k) * (qcko(i,k) - qrch)
!
!             rhbar(i) = rhbar(i) + qo(i,k) / qeso(i,k)
!
!  below lfc check if there is excess moisture to release latent heat
!
              if(k >= kbcon(i) .and. dq > 0.) then
                etah = .5 * (eta(i,k) + eta(i,k-1))
                dp = 1000. * del(i,k)
                if(ncloud > 0) then
                  ptem = c0t(i,k) + c1
                  qlk = dq / (eta(i,k) + etah * ptem * dz)
                  dellal(i,k) = etah * c1 * dz * qlk * g / dp
                else
                  qlk = dq / (eta(i,k) + etah * c0t(i,k) * dz)
                endif
                buo(i,k) = buo(i,k) - g * qlk
                qcko(i,k)= qlk + qrch
                pwo(i,k) = etah * c0t(i,k) * dz * qlk
                cnvwt(i,k) = etah * qlk * g / dp
              endif
!
!  compute buoyancy and drag for updraft velocity
!
              if(k >= kbcon(i)) then
                rfact =  1. + delta * cp * gamma &
                         * to(i,k) / hvap
                buo(i,k) = buo(i,k) + (g / (cp * to(i,k))) &
                         * dbyo(i,k) / (1. + gamma) &
                         * rfact
                val = 0.
                buo(i,k) = buo(i,k) + g * delta * &
                           max(val,(qeso(i,k) - qo(i,k)))
                drag(i,k) = max(xlamue(i,k),xlamud(i))
              endif
!
            endif
          endif
        enddo
      enddo
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
      do i = 1, im
        if (cnvflg(i)) then
          aa1(i) = 0.
        endif
      enddo
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k >= kbcon(i) .and. k < ktcon(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              aa1(i) = aa1(i) + buo(i,k) * dz1
            endif
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i) .and. aa1(i) <= 0.) cnvflg(i) = .false.
      enddo
!!
!> - If the updraft cloud work function is negative, convection does not occur, and the scheme returns to the calling routine.
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!
!  estimate the onvective overshooting as the level
!    where the [aafac * cloud work function] becomes zero,
!    which is the final cloud top
!    limited to the level of P/Ps=0.7
!
!> - Continue calculating the cloud work function past the point of neutral buoyancy to represent overshooting according to Han and Pan (2011) \cite han_and_pan_2011 . Convective overshooting stops when \f$ cA_u < 0\f$ where \f$c\f$ is currently 10%, or w
      do i = 1, im
        if (cnvflg(i)) then
          aa1(i) = aafac * aa1(i)
        endif
      enddo
!
      do i = 1, im
        flg(i) = cnvflg(i)
        ktcon1(i) = kbm(i)
      enddo
      do k = 2, km1
        do i = 1, im
          if (flg(i)) then
            if(k >= ktcon(i) .and. k < kbm(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              rfact =  1. + delta * cp * gamma &
                       * to(i,k) / hvap
              aa1(i) = aa1(i) +                   &
!    &                 dz1 * eta(i,k) * (g / (cp * to(i,k))) &
                       dz1 * (g / (cp * to(i,k))) &
                       * dbyo(i,k) / (1. + gamma) &
                       * rfact
!             val = 0.
!             aa1(i) = aa1(i) +
!!   &                 dz1 * eta(i,k) * g * delta *
!    &                 dz1 * g * delta *
!    &                 max(val,(qeso(i,k) - qo(i,k)))
              if(aa1(i) < 0.) then
                ktcon1(i) = k
                flg(i) = .false.
              endif
            endif
          endif
        enddo
      enddo
!
!  compute cloud moisture property, detraining cloud water
!    and precipitation in overshooting layers
!
!> - For the overshooting convection, calculate the moisture content of the entraining/detraining parcel as before. Partition convective cloud water and precipitation and detrain convective cloud water in the overshooting layers.
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k >= ktcon(i) .and. k < ktcon1(i)) then
              dz    = zi(i,k) - zi(i,k-1)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              qrch = qeso(i,k) &
                   + gamma * dbyo(i,k) / (hvap * (1. + gamma))
!cj
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              qcko(i,k) = ((1.-tem1)*qcko(i,k-1)+tem*0.5* &
                           (qo(i,k)+qo(i,k-1)))/factor
              qrcko(i,k) = qcko(i,k)
!cj
              dq = eta(i,k) * (qcko(i,k) - qrch)
!
!  check if there is excess moisture to release latent heat
!
              if(dq > 0.) then
                etah = .5 * (eta(i,k) + eta(i,k-1))
                dp = 1000. * del(i,k)
                if(ncloud > 0) then
                  ptem = c0t(i,k) + c1
                  qlk = dq / (eta(i,k) + etah * ptem * dz)
                  dellal(i,k) = etah * c1 * dz * qlk * g / dp
                else
                  qlk = dq / (eta(i,k) + etah * c0t(i,k) * dz)
                endif
                qcko(i,k) = qlk + qrch
                pwo(i,k) = etah * c0t(i,k) * dz * qlk
                cnvwt(i,k) = etah * qlk * g / dp
              endif
            endif
          endif
        enddo
      enddo
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
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k > kbcon1(i) .and. k < ktcon(i)) then
              dz    = zi(i,k) - zi(i,k-1)
              tem  = 0.25 * bb1 * (drag(i,k)+drag(i,k-1)) * dz
              tem1 = 0.5 * bb2 * (buo(i,k)+buo(i,k-1)) * dz
              ptem = (1. - tem) * wu2(i,k-1)
              ptem1 = 1. + tem
              wu2(i,k) = (ptem + tem1) / ptem1
              wu2(i,k) = max(wu2(i,k), 0.)
            endif
          endif
        enddo
      enddo
!
!  compute updraft velocity averaged over the whole cumulus
!
!> - Calculate the mean updraft velocity within the cloud (wc).
      do i = 1, im
        wc(i) = 0.
        sumx(i) = 0.
        wbar(i) = 0.
      enddo

      do k = 1,km
        do i = 1,im
          po1(i,k) = .5 * (pfld(i,k) + pfld(i,k+1))
        enddo
      enddo

      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k > kbcon1(i) .and. k < ktcon(i)) then
              dz = zi(i,k) - zi(i,k-1)
              tem = 0.5 * (sqrt(wu2(i,k)) + sqrt(wu2(i,k-1)))
              wc(i) = wc(i) + tem * dz
              tem  = 10. * dot(i,k)   * to(i,k)   / po1(i,k)
              tem1 = 10. * dot(i,k-1) * to(i,k-1) / po1(i,k-1)
              wbar(i) = wbar(i) + (-0.5 * rd/g) * (tem + tem1) * dz !grid-scale vertical velocity
              sumx(i) = sumx(i) + dz
            endif
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i)) then
          if(sumx(i) == 0.) then
             cnvflg(i)=.false.
          else
             wc(i) = wc(i) / sumx(i)
          endif
          val = 1.e-4
          if (wc(i) < val) cnvflg(i)=.false.
        endif
      enddo

      do i = 1, im
        if(cnvflg(i) .and. sumx(i) > 0.) then
          wbar(i) = wbar(i) / sumx(i)
        endif
      enddo
! compute mean cloud core fraction
! assume mean cloud core fraction to be the ratio of
! mean grid-scale vertical velocity (wbar) and mean updraft velocity
      do i = 1,im
        if (cnvflg(i)) then
          tem = wbar(i) / wc(i)
          tem = max(tem, 0.)
          clear(i) = 1. - tem
          clear(i) = max(min(clear(i), 1.0), 0.)
          if (wbar(i).gt.0. .and. wbar(i).gt.wc(i)) then
            cnvflg(i) = .false.
          endif
        endif
      enddo
!
! exchange ktcon with ktcon1
!
      do i = 1, im
        if(cnvflg(i)) then
          kk = ktcon(i)
          ktcon(i) = ktcon1(i)
          ktcon1(i) = kk
        endif
      enddo
!
!  this section is ready for cloud water
!
      if(ncloud > 0) then
!
!  compute liquid and vapor separation at cloud top
!
!> - => Separate the total updraft cloud water at cloud top into vapor and condensate.
      do i = 1, im
        if(cnvflg(i)) then
          k = ktcon(i) - 1
          gamma = el2orc * qeso(i,k) / (to(i,k)**2)
          qrch = qeso(i,k) &
               + gamma * dbyo(i,k) / (hvap * (1. + gamma))
          dq = qcko(i,k) - qrch
!
!  check if there is excess moisture to release latent heat
!
          if(dq > 0.) then
            qlko_ktcon(i) = dq
            qcko(i,k) = qrch
          endif
        endif
      enddo
      endif
!
!c--- compute precipitation efficiency in terms of windshear
!
!! - Calculate the wind shear and precipitation efficiency according to equation 58 in Fritsch and Chappell (1980) \cite fritsch_and_chappell_1980 :
!! \f[
!! E = 1.591 - 0.639\frac{\Delta V}{\Delta z} + 0.0953\left(\frac{\Delta V}{\Delta z}\right)^2 - 0.00496\left(\frac{\Delta V}{\Delta z}\right)^3
!! \f]
!! where \f$\Delta V\f$ is the integrated horizontal shear over the cloud depth, \f$\Delta z\f$, (the ratio is converted to units of \f$10^{-3} s^{-1}\f$). The variable "edt" is \f$1-E\f$ and is constrained to the range \f$[0,0.9]\f$.
      do i = 1, im
        if(cnvflg(i)) then
          vshear(i) = 0.
        endif
      enddo
      do k = 2, km
        do i = 1, im
          if (cnvflg(i)) then
            if(k > kb(i) .and. k <= ktcon(i)) then
              shear= sqrt((uo(i,k)-uo(i,k-1)) ** 2 &
                        + (vo(i,k)-vo(i,k-1)) ** 2)
              vshear(i) = vshear(i) + shear
            endif
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i)) then
          vshear(i) = 1.e3 * vshear(i) / (zi(i,ktcon(i))-zi(i,kb(i)))
          e1=1.591-.639*vshear(i) &
             +.0953*(vshear(i)**2)-.00496*(vshear(i)**3)
          edt(i)=1.-e1
          val =         .9
          edt(i) = min(edt(i),val)
          val =         .0
          edt(i) = max(edt(i),val)
        endif
      enddo
!
!c--- what would the change be, that a cloud with unit mass
!c--- will do to the environment?
!
!> ## Calculate the tendencies of the state variables (per unit cloud base mass flux) and the cloud base mass flux.
!> - Calculate the change in moist static energy, moisture mixing ratio, and horizontal winds per unit cloud base mass flux for all layers below cloud top from equations B.14 and B.15 from Grell (1993) \cite grell_1993, and for the cloud top from B.16 and
      do k = 1, km
        do i = 1, im
          if(cnvflg(i) .and. k <= kmax(i)) then
            dellah(i,k) = 0.
            dellaq(i,k) = 0.
            dellau(i,k) = 0.
            dellav(i,k) = 0.
          endif
        enddo
      enddo
!
!c--- changed due to subsidence and entrainment
!
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k > kb(i) .and. k < ktcon(i)) then
              dp = 1000. * del(i,k)
              dz = zi(i,k) - zi(i,k-1)
!
              dv1h = heo(i,k)
              dv2h = .5 * (heo(i,k) + heo(i,k-1))
              dv3h = heo(i,k-1)
              dv1q = qo(i,k)
              dv2q = .5 * (qo(i,k) + qo(i,k-1))
              dv3q = qo(i,k-1)
!
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1))
              tem1 = xlamud(i)
!cj
              dellah(i,k) = dellah(i,k) + &
           ( eta(i,k)*dv1h - eta(i,k-1)*dv3h &
          -  tem*eta(i,k-1)*dv2h*dz &
          +  tem1*eta(i,k-1)*.5*(hcko(i,k)+hcko(i,k-1))*dz &
               ) *g/dp
!cj
              dellaq(i,k) = dellaq(i,k) + &
           ( eta(i,k)*dv1q - eta(i,k-1)*dv3q &
          -  tem*eta(i,k-1)*dv2q*dz &
          +  tem1*eta(i,k-1)*.5*(qrcko(i,k)+qcko(i,k-1))*dz &
               ) *g/dp
!cj
              tem1=eta(i,k)*(uo(i,k)-ucko(i,k))
              tem2=eta(i,k-1)*(uo(i,k-1)-ucko(i,k-1))
              dellau(i,k) = dellau(i,k) + (tem1-tem2) * g/dp
!cj
              tem1=eta(i,k)*(vo(i,k)-vcko(i,k))
              tem2=eta(i,k-1)*(vo(i,k-1)-vcko(i,k-1))
              dellav(i,k) = dellav(i,k) + (tem1-tem2) * g/dp
!cj
            endif
          endif
        enddo
      enddo
!
!c------- cloud top
!
      do i = 1, im
        if(cnvflg(i)) then
          indx = ktcon(i)
          dp = 1000. * del(i,indx)
          dv1h = heo(i,indx-1)
          dellah(i,indx) = eta(i,indx-1) * &
                           (hcko(i,indx-1) - dv1h) * g / dp
          dv1q = qo(i,indx-1)
          dellaq(i,indx) = eta(i,indx-1) * &
                           (qcko(i,indx-1) - dv1q) * g / dp
          dellau(i,indx) = eta(i,indx-1) * &
                   (ucko(i,indx-1) - uo(i,indx-1)) * g / dp
          dellav(i,indx) = eta(i,indx-1) * &
                   (vcko(i,indx-1) - vo(i,indx-1)) * g / dp
!
!  cloud water
!
          dellal(i,indx) = eta(i,indx-1) * &
                           qlko_ktcon(i) * g / dp
        endif
      enddo
!
!
!  compute convective turn-over time
!
!> - Following Bechtold et al. (2008) \cite bechtold_et_al_2008, calculate the convective turnover time using the mean updraft velocity (wc) and the cloud depth. It is also proportional to the grid size (gdx).
      do i= 1, im
        if(cnvflg(i)) then
          tem = zi(i,ktcon1(i)) - zi(i,kbcon1(i))
          dtconv(i) = tem / wc(i)
          tfac = 1. + gdx(i) / 75000.
! reference from eq.3 in Zheng et al. 2016
!byl          tfac = 1. + log(25000./gdx(i))
          dtconv(i) = tfac * dtconv(i)
          dtconv(i) = max(dtconv(i),dtmin)
          dtconv(i) = max(dtconv(i),dt2)
          dtconv(i) = min(dtconv(i),dtmax)
        endif
      enddo
!
!> - Calculate advective time scale (tauadv) using a mean cloud layer wind speed.
      do i= 1, im
        if(cnvflg(i)) then
          sumx(i) = 0.
          umean(i) = 0.
        endif
      enddo
      do k = 2, km1
        do i = 1, im
          if(cnvflg(i)) then
            if(k >= kbcon1(i) .and. k < ktcon1(i)) then
              dz = zi(i,k) - zi(i,k-1)
              tem = sqrt(u1(i,k)*u1(i,k)+v1(i,k)*v1(i,k))
              umean(i) = umean(i) + tem * dz
              sumx(i) = sumx(i) + dz
            endif
          endif
        enddo
      enddo
      do i= 1, im
        if(cnvflg(i)) then
           umean(i) = umean(i) / sumx(i)
           umean(i) = max(umean(i), 1.)
           tauadv(i) = gdx(i) / umean(i)
        endif
      enddo
!
!  compute cloud base mass flux as a function of the mean
!     updraft velcoity
!
!> - From Han et al.'s (2017) \cite han_et_al_2017 equation 6, calculate cloud base mass flux as a function of the mean updraft velcoity.
!!  As discussed in Han et al. (2017) \cite han_et_al_2017 , when dtconv is larger than tauadv, the convective mixing is not fully conducted before the cumulus cloud is advected out of the grid cell. In this case, therefore, the cloud base mass flux is fu
      do i= 1, im
        if(cnvflg(i)) then
          k = kbcon(i)
          rho = po(i,k)*100. / (rd*to(i,k))
          tfac = tauadv(i) / dtconv(i)
          tfac = min(tfac, 1.)
          xmb(i) = tfac*betaw*rho*wc(i)
        endif
      enddo
!
!> - For scale-aware parameterization, the updraft fraction (sigmagfm) is first computed as a function of the lateral entrainment rate at cloud base (see Han et al.'s (2017) \cite han_et_al_2017 equation 4 and 5), following the study by Grell and Freitas
      do i = 1, im
        if(cnvflg(i)) then
          sigmagfm(i) = wbar(i) / wc(i)
          sigmagfm(i) = tem1 / garea(i)
          sigmagfm(i) = max(sigmagfm(i), 0.001)
          sigmagfm(i) = min(sigmagfm(i), 0.999)
        endif
      enddo
!
!> - Then, calculate the reduction factor (scaldfunc) of the vertical convective eddy transport of mass flux as a function of updraft fraction from the studies by Arakawa and Wu (2013) \cite arakawa_and_wu_2013 (also see Han et al.'s (2017) \cite han_et_a
      do i = 1, im
        if(cnvflg(i)) then
          scaldfunc(i) = clear(i) * (1.-sigma(i))
          scaldfunc(i) = max(min(scaldfunc(i), 1.0), 0.)
          xmb(i) = xmb(i) * scaldfunc(i)
          xmb(i) = min(xmb(i),xmbmax(i))
        endif
      enddo
!> ## For the "feedback control", calculate updated values of the state variables by multiplying the cloud base mass flux and the tendencies calculated per unit cloud base mass flux from the static control.
!! - Recalculate saturation specific humidity.
!
!c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      do k = 1, km
        do i = 1, im
          if (cnvflg(i) .and. k <= kmax(i)) then
            qeso(i,k) = 0.01 * fpvs(t1(i,k))      ! fpvs is in pa
!byl            call qsatq_cwb(1,t1(i,k),pfld(i,k),qeso(i,k))
            qeso(i,k) = eps * qeso(i,k) / (pfld(i,k) + epsm1*qeso(i,k))
            val     =             1.e-8
            qeso(i,k) = max(qeso(i,k), val )
          endif
        enddo
      enddo
!c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!> - Calculate the temperature tendency from the moist static energy and specific humidity tendencies.
!> - Update the temperature, specific humidity, and horiztonal wind state variables by multiplying the cloud base mass flux-normalized tendencies by the cloud base mass flux.
!> - Accumulate column-integrated tendencies.
      do i = 1, im
        delhbar(i) = 0.
        delqbar(i) = 0.
        deltbar(i) = 0.
        delubar(i) = 0.
        delvbar(i) = 0.
        qcond(i) = 0.
      enddo
      do k = 1, km
        do i = 1, im
          if (cnvflg(i)) then
            if(k > kb(i) .and. k <= ktcon(i)) then
              dellat = (dellah(i,k) - hvap * dellaq(i,k)) / cp
              t1(i,k) = t1(i,k) + dellat * xmb(i) * dt2
              q1(i,k) = q1(i,k) + dellaq(i,k) * xmb(i) * dt2
!             tem = 1./rcs(i)
!             u1(i,k) = u1(i,k) + dellau(i,k) * xmb(i) * dt2 * tem
!             v1(i,k) = v1(i,k) + dellav(i,k) * xmb(i) * dt2 * tem
              u1(i,k) = u1(i,k) + dellau(i,k) * xmb(i) * dt2
              v1(i,k) = v1(i,k) + dellav(i,k) * xmb(i) * dt2
              dp = 1000. * del(i,k)
              delhbar(i) = delhbar(i) + dellah(i,k)*xmb(i)*dp/g
              delqbar(i) = delqbar(i) + dellaq(i,k)*xmb(i)*dp/g
              deltbar(i) = deltbar(i) + dellat*xmb(i)*dp/g
              delubar(i) = delubar(i) + dellau(i,k)*xmb(i)*dp/g
              delvbar(i) = delvbar(i) + dellav(i,k)*xmb(i)*dp/g
            endif
          endif
        enddo
      enddo
!
!> - Recalculate saturation specific humidity using the updated temperature.
      do k = 1, km
        do i = 1, im
          if (cnvflg(i)) then
            if(k > kb(i) .and. k <= ktcon(i)) then
              qeso(i,k) = 0.01 * fpvs(t1(i,k))      ! fpvs is in pa
!byl              call qsatq_cwb(1,t1(i,k),pfld(i,k),qeso(i,k))
              qeso(i,k) = eps * qeso(i,k)/(pfld(i,k) + epsm1*qeso(i,k))
              val     =             1.e-8
              qeso(i,k) = max(qeso(i,k), val )
            endif
          endif
        enddo
      enddo
!
!> - Add up column-integrated convective precipitation by multiplying the normalized value by the cloud base mass flux.
      do i = 1, im
        rntot(i) = 0.
        delqev(i) = 0.
        delq2(i) = 0.
        flg(i) = cnvflg(i)
      enddo
      do k = km, 1, -1
        do i = 1, im
          if (cnvflg(i)) then
            if(k < ktcon(i) .and. k > kb(i)) then
              rntot(i) = rntot(i) + pwo(i,k) * xmb(i) * .001 * dt2
            endif
          endif
        enddo
      enddo
!
! evaporating rain
!
!> - Determine the evaporation of the convective precipitation and update the integrated convective precipitation.
!> - Update state temperature and moisture to account for evaporation of convective precipitation.
!> - Update column-integrated tendencies to account for evaporation of convective precipitation.
      do k = km, 1, -1
        do i = 1, im
          if (k <= kmax(i)) then
            deltv(i) = 0.
            delq(i) = 0.
            qevap(i) = 0.
            if(cnvflg(i)) then
              if(k < ktcon(i) .and. k > kb(i)) then
                rn(i) = rn(i) + pwo(i,k) * xmb(i) * .001 * dt2
              endif
            endif
            if(flg(i) .and. k < ktcon(i)) then
              evef = edt(i) * evfact
              if(islimsk(i) == 1) evef=edt(i) * evfactl
!             if(islimsk(i) == 1) evef=.07
!             if(islimsk(i) == 1) evef = 0.
              qcond(i) = evef * (q1(i,k) - qeso(i,k)) &
                       / (1. + el2orc * qeso(i,k) / t1(i,k)**2)
              dp = 1000. * del(i,k)
              if(rn(i) > 0. .and. qcond(i) < 0.) then
                qevap(i) = -qcond(i) * (1.-exp(-.32*sqrt(dt2*rn(i))))
                qevap(i) = min(qevap(i), rn(i)*1000.*g/dp)
                delq2(i) = delqev(i) + .001 * qevap(i) * dp / g
              endif
              if(rn(i) > 0. .and. qcond(i) < 0. .and. &
                 delq2(i) > rntot(i)) then
                qevap(i) = 1000.* g * (rntot(i) - delqev(i)) / dp
                flg(i) = .false.
              endif
              if(rn(i) > 0. .and. qevap(i) > 0.) then
                tem  = .001 * dp / g
                tem1 = qevap(i) * tem
                if(tem1 > rn(i)) then
                  qevap(i) = rn(i) / tem
                  rn(i) = 0.
                else
                  rn(i) = rn(i) - tem1
                endif
                q1(i,k) = q1(i,k) + qevap(i)
                t1(i,k) = t1(i,k) - elocp * qevap(i)
                deltv(i) = - elocp*qevap(i)/dt2
                delq(i) =  + qevap(i)/dt2
                delqev(i) = delqev(i) + .001*dp*qevap(i)/g
              endif
              delqbar(i) = delqbar(i) + delq(i)*dp/g
              deltbar(i) = deltbar(i) + deltv(i)*dp/g
            endif
          endif
        enddo
      enddo
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
      do i = 1, im
        if(cnvflg(i)) then
          if (rn(i) < 1.e-13) rn(i) = 0.
          if(rn(i) < 0. .or. .not.flg(i)) rn(i) = 0.
          ktop(i) = ktcon(i)
          kbot(i) = kbcon(i)
          kcnv(i) = 2
        endif
      enddo
!
!  convective cloud water
!
!> - Calculate shallow convective cloud water.
      do k = 1, km
        do i = 1, im
          if (cnvflg(i)) then
            if (k >= kbcon(i) .and. k < ktcon(i)) then
              cnvw(i,k) = cnvwt(i,k) * xmb(i) * dt2
            endif
          endif
        enddo
      enddo

!
!  convective cloud cover
!
!> - Calculate convective cloud cover, which is used when pdf-based cloud fraction is used (i.e., pdfcld=.true.).
      do k = 1, km
        do i = 1, im
          if (cnvflg(i)) then
            if (k >= kbcon(i) .and. k < ktcon(i)) then
              cnvc(i,k) = 0.04 * log(1. + 675. * eta(i,k) * xmb(i))
              cnvc(i,k) = min(cnvc(i,k), 0.2)
              cnvc(i,k) = max(cnvc(i,k), 0.0)
            endif
          endif
        enddo
      enddo
!
!  cloud water
!
!> - Separate detrained cloud water into liquid and ice species as a function of temperature only.
      if (ncloud > 0) then
!
      do k = 1, km1
        do i = 1, im
          if (cnvflg(i)) then
!           if (k > kb(i) .and. k <= ktcon(i)) then
            if (k >= kbcon(i) .and. k <= ktcon(i)) then
              tem  = dellal(i,k) * xmb(i) * dt2
              tem1 = max(0.0, min(1.0, (tcr-t1(i,k))*tcrf))
              if (qi(i,k) > -999.0) then
                qi(i,k) = qi(i,k) + tem * tem1            ! ice
                ql(i,k) = ql(i,k) + tem *(1.0-tem1)       ! water
              else
                ql(i,k) = ql(i,k) + tem
              endif
            endif
          endif
        enddo
      enddo
!
      endif
!
! hchuang code change
!
!> - Calculate and retain the updraft mass flux for dust transport by cumulus convection.
!
!> - Calculate the updraft convective mass flux.
      do k = 1, km
        do i = 1, im
          if(cnvflg(i)) then
            if(k >= kb(i) .and. k < ktop(i)) then
              ud_mf(i,k) = eta(i,k) * xmb(i) * dt2
            endif
          endif
        enddo
      enddo
!> - save the updraft convective mass flux at cloud top.
      do i = 1, im
        if(cnvflg(i)) then
           k = ktop(i)-1
           dt_mf(i,k) = ud_mf(i,k)
        endif
      enddo
!!
      return
      end
