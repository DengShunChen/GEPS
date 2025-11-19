!> \defgroup SAMF Scale-Aware Mass-Flux Deep Convection
!! @{
!!  \brief The scale-aware mass-flux (SAMF) deep convection scheme is an updated version of the previous Simplified Arakawa-Schubert (SAS) scheme with scale and aerosol awareness and parameterizes the effect of deep convection on the environment (represen
!!
!!  The SAS scheme uses the working concepts put forth in Arakawa and Schubert (1974) \cite arakawa_and_schubert_1974 but includes modifications and simplifications from Grell (1993) \cite grell_1993 such as saturated downdrafts and only one cloud type (t
!!
!! The SAMF scheme updates the SAS scheme with scale- and aerosol-aware parameterizations from Han et al. (2017) \cite han_et_al_2017 based on the studies by Arakawa and Wu (2013) \cite arakawa_and_wu_2013 and Grell and Freitas (2014) \cite grell_and_frei
!!
!!  \section diagram Calling Hierarchy Diagram
!!  \image html SAMF_Flowchart.png "Diagram depicting how the SAMF deep convection scheme is called from the FV3GFS physics time loop" height=2cm
!!  \section intraphysics Intraphysics Communication
!!  This space is reserved for a description of how this scheme uses information from other scheme types and/or how information calculated in this scheme is used in other scheme types.

!> \file samfdeepcnv.f
!!  Contains the entire SAMF deep convection scheme.

!>  \brief This subroutine contains the entirety of the SAMF deep convection scheme.
!!
!!  For grid sizes larger than threshold value, as in Grell (1993) \cite grell_1993 , the SAMF deep convection scheme can be described in terms of three types of "controls": static, dynamic, and feedback. The static control component consists of the simpl
!!
!! For grid sizes smaller than threshold value, the cloud base mass flux in the SAMF scheme is determined by the cumulus updraft velocity averaged ove the whole cloud depth (Han et al., 2017 \cite han_et_al_2017 ), which in turn, determines changes of the
!!
!!  \param[in] im number of used points
!!  \param[in] ix horizontal dimension
!!  \param[in] km vertical layer dimension
!!  \param[in] delt physics time step in seconds
!!  \param[in] delp pressure difference between level k and k+1 (Pa)
!!  \param[in] prslp mean layer presure (Pa)
!!  \param[in] psp surface pressure (Pa)
!!  \param[in] phil layer geopotential (\f$m^2/s^2\f$)
!!  \param[inout] ql cloud water (kg/kg)
!!  \param[inout] qi cloud ice (kg/kg)
!!  \param[inout] q1 updated tracers (kg/kg)
!!  \param[inout] t1 updated temperature (K)
!!  \param[inout] u1 updated zonal wind (\f$m s^{-1}\f$)
!!  \param[inout] v1 updated meridional wind (\f$m s^{-1}\f$)
!!  \param[out] cldwrk cloud workfunction (\f$m^2/s^2\f$)
!!  \param[out] rn convective rain (m)
!!  \param[out] kbot index for cloud base
!!  \param[out] ktop index for cloud top
!!  \param[out] kcnv flag to denote deep convection (0=no, 1=yes)
!!  \param[in] islimsk sea/land/ice mask (=0/1/2)
!!  \param[in] garea area of grid box (\f$m^2\f$)
!!  \param[in] dot layer mean vertical velocity (Pa/s)
!!  \param[in] ncloud number of cloud species
!!  \param[out] ud_mf updraft mass flux multiplied by time step (\f$kg/m^2\f$)
!!  \param[out] dd_mf downdraft mass flux multiplied by time step (\f$kg/m^2\f$)
!!  \param[out] dt_mf ud_mf at cloud top (\f$kg/m^2\f$)
!!  \param[out] cnvw convective cloud water (kg/kg)
!!  \param[out] cnvc convective cloud cover (unitless)
!!  \param[in] clam coefficient for entrainment rate
!!  \param[in] c0s convective rain conversion parameter (1/m)
!!  \param[in] c1 conversion parameter of detrainment from liquid water into grid-scale cloud water (1/m)
!!  \param[in] betal fraction factor of downdraft air mass reaching ground surface over land
!!  \param[in] betas fraction factor of downdraft air mass reaching ground surface over sea
!!  \param[in] evfact evaporation factor from convective rain
!!  \param[in] evfactl evaporation factor from convective rain over land
!!  \param[in] pgcon reduction factor in momentum transport due to convection induced pressure gradient force
!!  \param[in] asolfac aerosol-aware parameter inversely proportional to CCN number concentraion
!!
!!  \section general General Algorithm
!!  -# Compute preliminary quantities needed for static, dynamic, and feedback control portions of the algorithm.
!!  -# Perform calculations related to the updraft of the entraining/detraining cloud model ("static control").
!!  -# Perform calculations related to the downdraft of the entraining/detraining cloud model ("static control").
!!
!!  -# For grid sizes larger than the threshold value (currently 8 km):
!!  -# 1) Using the updated temperature and moisture profiles that were modified by the convection on a short time-scale, recalculate the total cloud work function to determine the change in the cloud work function due to convection, or the stabilizing ef
!!  -# 2) For the "dynamic control", using a reference cloud work function, estimate the change in cloud work function due to the large-scale dynamics. Following the quasi-equilibrium assumption, calculate the cloud base mass flux required to keep the lar
!!
!!  -# For grid sizes smaller than the threshold value (currently 8 km):
!!  -# 1) compute the cloud base mass flux using the cumulus updraft velocity averaged ove the whole cloud depth.
!!
!!  -# For scale awareness, the updraft fraction (sigma) is obtained as a function of cloud base entrainment. Then, the final cloud base mass flux is obtained by the original mass flux multiplied by the \f$(1-sigma)^2\f$.
!!  -# For the "feedback control", calculate updated values of the state variables by multiplying the cloud base mass flux and the tendencies calculated per unit cloud base mass flux from the static control.
!!  \section detailed Detailed Algorithm
!!  @{
      subroutine samfdeepcnv_kh_gpu(nxjp,ix,km,delt,delp,prslp,psp,phil,ql, &
           qi,q1,t1,u1,v1,cldwrk,rn,kbot,ktop,kcnv,islimsk,garea, &
           dot,ncloud,cnvw,cnvc, &
!xb110>> for flash parameterization
           snow_flx,ptu,pqu, &
!xb99>> test dissipation
           dissdc1, myrank, icheck, kcheck, jjcheck, rcheck)
!xb110<<
!,ud_mf,dd_mf,dt_mf,cnvw,cnvc, &
!           clam,c0s,c1,betal,betas,evfact,evfactl,pgcon,asolfac)
!
      !$acc routine(fpvs_gpu) seq
      use machine , only : kind_phys
      use param, only: my, my_max
      use index, only: jlistnum, jlist1
!      use mpe
      use rank
!      use index
!byl      use funcphys , only : fpvs
      use physcons, grav => con_g, cp => con_cp, hvap => con_hvap &
      ,             rv => con_rv, fv => con_fvirt, t0c => con_t0c &
      ,             rd => con_rd, cvap => con_cvap, cliq => con_cliq &
      ,             eps => con_eps, epsm1 => con_epsm1
      implicit none
!
      real fpvs, fpvs_gpu
!
      integer, intent(in) :: myrank, icheck, kcheck, jjcheck, rcheck
      integer, intent(in)  :: ix,  km, ncloud, nxjp(my)
      integer, intent(in)  :: islimsk(ix, my_max)
      real(kind=kind_phys), intent(in) ::  delt
      real(kind=kind_phys), intent(in) :: psp(ix, my_max), delp(ix,km, my_max), &
         prslp(ix,km, my_max),  garea(ix, my_max), dot(ix,km, my_max), phil(ix,km, my_max)

      integer, intent(inout)  :: kcnv(ix, my_max)
      real(kind=kind_phys), intent(inout) ::   ql(ix,km, my_max), &
         q1(ix,km, my_max), t1(ix,km, my_max),   u1(ix,km, my_max), v1(ix,km, my_max),    &
         qi(ix,km, my_max)

      integer, intent(out) :: kbot(ix, my_max), ktop(ix, my_max)
      real(kind=kind_phys), intent(out) :: cldwrk(ix, my_max),rn(ix, my_max)
      real(kind=kind_phys) cnvw(ix,km, my_max),  cnvc(ix,km, my_max),  &
         ud_mf(ix,km, my_max),dd_mf(ix,km, my_max), dt_mf(ix,km, my_max)
!-----test dissipation from convective
      real(kind=kind_phys), intent(inout):: dissdc1(ix,km, my_max)

!
!------local variables
      integer              i, indx, jmn, k, kk, km1, n, j, jj
      integer              myim(my_max)
!     integer              latd,lond
!
      real(kind=kind_phys) clam,    cxlamu,  cxlamd, &
                           xlamde,  xlamdd, &
                           crtlamu, crtlamd
!
!     real(kind=kind_phys) detad
      real(kind=kind_phys) adw,     aup,     aafac, &
                           beta,    betal,   betas, &
                           c0l,     c0s,     d0, &
                           c1,      asolfac, &
                           dellat,  delta,   desdt,   dg, &
                           dh,      dhh,     dp, &
                           dq,      dqsdp,   dqsdt,   dt, &
                           dt2,     dtmax,   dtmin, &
                           dxcrtas, dxcrtuf, &
                           dv1h,    dv2h,    dv3h, &
                           dv1q,    dv2q,    dv3q, &
                           dz,      dz1,     e1,      edtmax, &
                           edtmaxl, edtmaxs, el2orc,  elocp, &
                           es,      etah, &
                           cthk,    dthk, &
                           evef,    evfact,  evfactl, fact1, &
                           fact2,   factor, &
                           g,       gamma,   pprime,  cm, &
                           qlk,     qrch,    qs, &
                           rain,    rfact,   shear,   tfac, &
                           val,     val1,    val2, &
                           w1,      w1l,     w1s,     w2, &
                           w2l,     w2s,     w3,      w3l, &
                           w3s,     w4,      w4l,     w4s, &
                           rho,     betaw, &
                           xdby,    xpw,     xpwd,          &
!    &                     xqrch,   mbdt,    tem, &
                           xqrch,   tem,     tem1,    tem2, &
                           ptem,    ptem1,   ptem2, &
                           pgcon
!
      integer              kb(ix, my_max), kbcon(ix, my_max), kbcon1(ix, my_max), &
                           ktcon(ix, my_max), ktcon1(ix, my_max), ktconn(ix, my_max), &
                           jmin(ix, my_max), lmin, kbmax(ix, my_max), &
                           kbm(ix, my_max), kmax(ix, my_max)
!
!     real(kind=kind_phys) aa1(im),     acrt(im),   acrtfct(im),
      real(kind=kind_phys) aa1(ix, my_max), &
                           ps(ix, my_max), del(ix,km, my_max),  prsl(ix,km, my_max), &
                           umean,   tauadv(ix, my_max), gdx(ix, my_max), &
                           delhbar(ix, my_max), delq,   delq2, &
                           delqbar(ix, my_max), delqev(ix, my_max), deltbar(ix, my_max), &
                           deltv,   dtconv(ix, my_max), edt(ix, my_max), &
                           edto(ix, my_max),    edtx(ix, my_max),   fld(ix, my_max), &
                           hcdo(ix,km, my_max), hmax(ix, my_max),   hmin, &
                           ucdo(ix,km, my_max), vcdo(ix,km, my_max),aa2, &
                           pdot(ix, my_max),    po, &
                           pwavo(ix, my_max),   pwevo,  mbdt(ix, my_max), &
                           qcdo(ix,km, my_max), qcond,  qevap, &
                           rntot(ix, my_max),   vshear, xaa0(ix, my_max), &
                           xk(ix, my_max),      xlamd(ix, my_max),  cina(ix, my_max), &
                           xmb(ix, my_max),     xmbmax(ix, my_max), xpwav(ix, my_max), &
                           xpwev(ix, my_max),   xlamx, &
                           delubar(ix, my_max),delvbar(ix, my_max)
!
      real(kind=kind_phys) c0(ix, my_max)
!cj
      real(kind=kind_phys) cinpcr,  cinpcrmx,  cinpcrmn, &
                           cinacr,  cinacrmx,  cinacrmn
!cj
!
!  parameters for updraft velocity calculation
      real(kind=kind_phys) bet1,    cd1,     f1,      gam1, &
                           bb1,     bb2
!    &                     bb1,     bb2,     wucb
!
!  physical parameters
      parameter(g=grav,asolfac=0.958)
!byl      parameter(g=grav)
      parameter(elocp=hvap/cp,el2orc=hvap*hvap/(rv*cp))
!      parameter(c0s=.002,c1=.0015,d0=.01)     !K4
!      parameter(c0s=.001,c1=.0015,d0=.01)     !K5
      parameter(c0s=.0025,c1=.0015,d0=.01)     !K6
!byl      parameter(d0=.01)
!     parameter(c0l=c0s*asolfac)
!
! asolfac: aerosol-aware parameter based on Lim (2011)
!      asolfac= cx / c0s(=.002)
!      cx = min([-0.7 ln(Nccn) + 24]*1.e-4, c0s)
!      Nccn: CCN number concentration in cm^(-3)
!      Until a realistic Nccn is provided, Nccns are assumed
!      as Nccn=100 for sea and Nccn=1000 for land
!
      parameter(cm=1.0,delta=fv)
      parameter(fact1=(cvap-cliq)/rv,fact2=hvap/rv-fact1*t0c)
      parameter(cthk=200.,dthk=25.)
      parameter(cinpcrmn=120.)
!      parameter(cinpcrmx=180.,cinpcrmn=120.)
!     parameter(cinacrmx=-120.,cinacrmn=-120.)
      parameter(cinacrmx=-120.,cinacrmn=-80.)
      parameter(bet1=1.875,cd1=.506,f1=2.0,gam1=.5)
      parameter(betaw=.03,dxcrtas=8.e3,dxcrtuf=15.e3)
!
!  local variables and arrays
      real(kind=kind_phys) pfld(ix,km, my_max),    to(ix,km, my_max),     qo(ix,km, my_max), &
                           uo(ix,km, my_max),      vo(ix,km, my_max),     qeso(ix,km, my_max)
!  for updraft velocity calculation
      real(kind=kind_phys) wu2(ix,km, my_max),     buo(ix,km, my_max),    drag, dragm1
      real(kind=kind_phys) wc(ix, my_max),         scaldfunc(ix, my_max), sigmagfm(ix, my_max)
!
!  cloud water
!     real(kind=kind_phys) tvo(im,km)
      real(kind=kind_phys) qlko_ktcon(ix, my_max), dellal(ix,km, my_max), tvo(ix,km), &
                           dbyo(ix,km, my_max),    zo(ix,km, my_max), &
                           xlamue(ix,km, my_max),  xlamud(ix,km, my_max), &
                           fent1,   fent2,  frh(ix,km, my_max), &
                           heo(ix,km, my_max),     heso(ix,km, my_max), &
                           qrcd(ix,km, my_max),    dellah(ix,km, my_max), dellaq(ix,km, my_max), &
                           dellau(ix,km, my_max),  dellav(ix,km, my_max), hcko(ix,km, my_max), &
                           ucko(ix,km, my_max),    vcko(ix,km, my_max),   qcko(ix,km, my_max), &
                           eta(ix,km, my_max),     etad(ix,km, my_max),   zi(ix,km, my_max), &
                           qrcko(ix,km, my_max),   qrcdo(ix,km, my_max), &
                           pwo(ix,km, my_max),     pwdo(ix,km, my_max),   c0t(ix,km, my_max), &
                           tx1(ix, my_max),        sumx,      cnvwt(ix,km, my_max)
!    &,                    rhbar(im)
!
      logical totflg, cnvflg(ix, my_max), asqecflg(ix, my_max), flg(ix, my_max)
!
!    asqecflg: flag for the quasi-equilibrium assumption of Arakawa-Schubert
!
!     real(kind=kind_phys) pcrit(15), acritt(15), acrit(15)
!!    save pcrit, acritt
!     data pcrit/850.,800.,750.,700.,650.,600.,550.,500.,450.,400.,
!    &           350.,300.,250.,200.,150./
!     data acritt/.0633,.0445,.0553,.0664,.075,.1082,.1521,.2216,
!    &           .3151,.3677,.41,.5255,.7663,1.1686,1.6851/
!  gdas derived acrit
!     data acritt/.203,.515,.521,.566,.625,.665,.659,.688,
!    &            .743,.813,.886,.947,1.138,1.377,1.896/
      real(kind=kind_phys) tf, tcr, tcrf
      parameter (tf=233.16, tcr=263.16, tcrf=1.0/(tcr-tf))
!
!xb110>> for lightning parameterization
      real(kind=kind_phys) snow_flx(ix,km, my_max),rainf,&
                           ptu(ix,km, my_max),pqu(ix,km, my_max)
!<---for scale-aware parameterization (Kwon and Hong 2017)
      real(kind=kind_phys) cinpcri, frh_sum(ix, my_max), wbar(ix, my_max),      &
                           clear(ix, my_max),sigma(ix, my_max)
      real(kind=kind_phys) po1, po1m1
      real(kind=kind_phys) sigma_con, pi, dx1km, dx5km, dx250m
      parameter (pi=3.14159)
      parameter (dx1km = 1000., dx5km = 5000., dx250m = 250.)
      
      
      ! for GPU porting
      ! register
      real :: qesor, qor, tor, tem3, tem4, tem5, tem6, tem7, tem8, &
         etar, etar1, etadr, etadr1, uor, uor1, vor, vor1, zir1, zir, &
         xlamudr1, xlamudr, xlamuer1, xlamuer, hckor1, heor1, uckor1, uor1, &
         vckor1, vor1, uor, vor, heor, heor1, frh_sumr, cinar, hmaxr, zor1, &
         zor, pfldr1, pfldr, tor1, qor1, uor, uor1, vor, vor1, pdotr, aa1r, &
         wu2r1, wu2r, zir1, zir, buor1, buor, pwavor, qckor1, qckor, qrckor1, &
         qrckor, c0tr, wcr, wbarr

      integer :: kmaxr, kbconr, kbr, kbmaxr, kbcon1r, jminr, ktconr, ktconnr, ktcon1r
      logical :: flgr
      
      integer, parameter :: nxpvs = 7501
      real :: c1xpvs,c2xpvs,tbpvs(nxpvs)
      common/fpvscom/ c1xpvs,c2xpvs,tbpvs(nxpvs)
      
      integer :: async_id = 1
      !character(len=4):: device_str = 'gpu'
!xb110<<
!c-----------------------------------------------------------------------
!>  ## Compute preliminary quantities needed for static, dynamic, and feedback control portions of the algorithm.
!>  - Convert input pressure terms to centibar units.
!************************************************************************
!     convert input Pa terms to Cb terms  -- Moorthi
!byl      ps   = psp   * 0.001
!byl      prsl = prslp * 0.001
!byl      del  = delp  * 0.001
      !$acc data create(myim, ps, prsl, del, cnvflg, &
      !$acc&     mbdt, kbcon, ktcon, ktconn, dtconv, pdot, jmin, &
      !$acc&     qlko_ktcon, edt, edto, edtx, aa1, xaa0, cina, pwavo, &
      !$acc&     xpwav, xpwev, gdx, sigma, frh_sum, &
      !$acc&     sigmagfm, c0, c0t, ud_mf, dd_mf, dt_mf, kbmax, kbm, kmax, &
      !$acc&     tx1, zo, zi, xlamue, pfld, eta, frh, hcko, &
      !$acc&     qcko, qrcko, ucko, vcko, etad, hcdo, qcdo, ucdo, vcdo, &
      !$acc&     qrcd, qrcdo, dbyo, pwo, pwdo, dellal, to, qo, uo, vo, wu2, &
      !$acc&     buo, cnvwt, qeso, heo, heso, hmax, kb, flg, &
      !$acc&     xlamud, kbcon1, xmbmax, ktcon1, wc, wbar, &
      !$acc&     xlamd, dellah, dellaq, dellau, dellav, asqecflg, &
      !$acc&     tauadv, xmb, fld, xk, scaldfunc, delhbar, delqbar, deltbar, &
      !$acc&     delubar, delvbar, rntot, delqev, clear) async(async_id)
      !$acc parallel loop async(async_id) private(j)
      do jj = 1, jlistnum
         j = jlist1(jj)
         myim(jj) = nxjp(j)
      end do

      
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, jlistnum
         do k = 1, km
            !$acc loop vector
            do i = 1, myim(jj)
              prsl(i, k, jj) = prslp(i, k, jj)
              del(i, k, jj)  = delp(i, k, jj)
            enddo
         enddo
      end do

   !xb110>
      !$acc parallel loop gang vector collapse(3) async(async_id)
      do jj = 1, jlistnum
         do k = 1, km
            do i = 1, ix
               ptu(i, k, jj) = 0.
               pqu(i, k, jj) = 0.
               snow_flx(i, k, jj) = 0. !xb110, snow flux (kgm-2s-1)
            end do
         end do
      end do
   !xb110<
   !************************************************************************
   !
   !
         km1 = km - 1
         cinpcrmx=240.
      !$acc parallel loop gang vector collapse(2) async(async_id)
      do jj = 1, jlistnum
         do i = 1, ix
            if (i .le. myim(jj)) then
               ps(i, jj) = psp(i, jj)
   !>  - Initialize column-integrated and other single-value-per-column variable arrays.
   !
   !  initialize arrays
   !
               cnvflg(i, jj) = .true.
               rn(i, jj)=0.
               mbdt(i, jj)=10.
               kbot(i, jj)=km+1
               ktop(i, jj)=0
               kbcon(i, jj)=km
               ktcon(i, jj)=1
               ktconn(i, jj)=1
               dtconv(i, jj) = 3600.
               cldwrk(i, jj) = 0.
               pdot(i, jj) = 0.
               jmin(i, jj) = 1
               qlko_ktcon(i, jj) = 0.
               edt(i, jj)  = 0.
               edto(i, jj) = 0.
               edtx(i, jj) = 0.
      !        acrt(i) = 0.
      !        acrtfct(i) = 1.
               aa1(i, jj)  = 0.
               xaa0(i, jj) = 0.
               cina(i, jj) = 0.
               xpwav(i, jj)= 0.
               xpwev(i, jj)= 0.
               gdx(i, jj) = sqrt(garea(i, jj))

               sigma(i, jj) = 0.
               frh_sum(i, jj) = 0.
               sigmagfm(i, jj) = 0.
   !
   !>  - determine aerosol-aware rain conversion parameter over land
               if(islimsk(i, jj) == 1) then
                  c0(i, jj) = c0s*asolfac
               else
                  c0(i, jj) = c0s
               endif
              
   !
   !  define top layer for search of the downdraft originating layer
   !  and the maximum thetae for updraft
   !
   !>  - Determine maximum indices for the parcel starting point (kbm), LFC (kbmax), and cloud top (kmax).
               kbmax(i, jj) = km
               kbm(i, jj)   = km
               kmax(i, jj)  = km
               tx1(i, jj)   = 1.0 / ps(i, jj)
              
               !$acc loop seq
               do k = 1, km
                  if (prsl(i, k, jj)*tx1(i, jj) > 0.04) kmax(i, jj)  = k + 1
                  if (prsl(i, k, jj)*tx1(i, jj) > 0.45) kbmax(i, jj) = k + 1
                  if (prsl(i, k, jj)*tx1(i, jj) > 0.70) kbm(i, jj)   = k + 1
               enddo
              
               kmax(i, jj)  = min(km, kmax(i, jj))
               kbmax(i, jj) = min(kbmax(i, jj), kmax(i, jj))
               kbm(i, jj)   = min(kbm(i, jj), kmax(i, jj))
            end if
         enddo
      end do
   !>  - determine rain conversion parameter above the freezing level which exponentially decreases with decreasing temperature from Han et al.'s (2017) \cite han_et_al_2017 equation 8.
      !$acc parallel loop gang collapse(2) async(async_id) private(tem, tem1)
      do jj = 1, jlistnum
         do k = 1, km
            !$acc loop vector
            do i = 1, myim(jj)
               if(t1(i, k, jj) > 273.16) then
                  c0t(i, k, jj) = c0(i, jj)
               else
                  tem = d0 * (t1(i, k, jj) - 273.16)
                  tem1 = exp(tem)
                  c0t(i, k, jj) = c0(i, jj) * tem1
               endif
            enddo
         enddo
      end do
   !>  - Initialize convective cloud water and cloud cover to zero.
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, jlistnum
         do k = 1, km
            !$acc loop vector
            do i = 1, myim(jj)
               cnvw(i, k, jj) = 0.
               cnvc(i, k, jj) = 0.
   ! hchuang code change
   !>  - Initialize updraft and downdraft mass fluxes to zero.
               ud_mf(i, k, jj) = 0.
               dd_mf(i, k, jj) = 0.
               dt_mf(i, k, jj) = 0.
               dissdc1(i, k, jj) = 0.
            enddo
         enddo
      end do
!
!     do k = 1, 15
!       acrit(k) = acritt(k) * (975. - pcrit(k))
!     enddo
!
      dt2 = delt
!     val   =         1200.
      val   =         600.
      dtmin = max(dt2, val )
!     val   =         5400.
      val   =         10800.
      dtmax = max(dt2, val )
!  model tunable parameters are all here
      edtmaxl = .3
      edtmaxs = .3
      clam    = .1
      aafac   = .1
!     betal   = .15
!     betas   = .15
      betal   = .05
      betas   = .05
!     evef    = 0.07
      evfact  = 0.3
      evfactl = 0.3
!
      crtlamu = 1.0e-4
      crtlamd = 1.0e-4
!
      cxlamu  = 1.0e-3
      cxlamd  = 1.0e-4
      xlamde  = 1.0e-4
      xlamdd  = 1.0e-4
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
   !  hydrostatic height assume zero terr and initially assume
   !    updraft entrainment rate as an inverse function of height
   !
   !>  - Calculate hydrostatic height at layer centers assuming a flat surface (no terrain) from the geopotential.
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, jlistnum
         do k = 1, km
            !$acc loop vector
            do i = 1, myim(jj)
               zo(i, k, jj) = phil(i, k, jj) / g
            enddo
         enddo
      end do
   !>  - Calculate interface height and the initial entrainment rate as an inverse function of height.
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, jlistnum
         do k = 1, km1
            !$acc loop vector
            do i = 1, myim(jj)
               zi(i, k, jj) = 0.5*(zo(i, k, jj)+zo(i, k+1, jj))
               xlamue(i, k, jj) = clam / zi(i, k, jj)
   !           xlamue(i,k) = max(xlamue(i,k), crtlamu)
            enddo
         enddo
      end do
   !
   !c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !   convert surface pressure to mb from cb
   !
   !>  - Convert prsl from centibar to millibar, set normalized mass fluxes to 1, cloud properties to 0, and save model state variables (after advection/turbulence).
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, jlistnum
         do k = 1, km
            !$acc loop vector
            do i = 1, myim(jj)
               if (k <= kmax(i, jj)) then
                  pfld(i, k, jj) = prsl(i, k, jj) * 10.0
                  eta(i, k, jj)  = 1.
                  frh(i, k, jj)  = 0.
                  hcko(i, k, jj) = 0.
                  qcko(i, k, jj) = 0.
                  qrcko(i ,k, jj)= 0.
                  ucko(i, k, jj) = 0.
                  vcko(i, k, jj) = 0.
                  etad(i, k, jj) = 1.
                  hcdo(i, k, jj) = 0.
                  qcdo(i, k, jj) = 0.
                  ucdo(i, k, jj) = 0.
                  vcdo(i, k, jj) = 0.
                  qrcd(i, k, jj) = 0.
                  qrcdo(i ,k, jj)= 0.
                  dbyo(i, k, jj) = 0.
                  pwo(i, k, jj)  = 0.
                  pwdo(i, k, jj) = 0.
                  dellal(i, k, jj) = 0.
                  to(i, k, jj)   = t1(i, k, jj)
                  qo(i, k, jj)   = q1(i, k, jj)
                  uo(i, k, jj)   = u1(i, k, jj)
                  vo(i, k, jj)   = v1(i, k, jj)
   !              uo(i,k)   = u1(i,k) * rcs(i)
   !              vo(i,k)   = v1(i,k) * rcs(i)
                  wu2(i, k, jj)  = 0.
                  buo(i, k, jj)  = 0.
                  cnvwt(i, k, jj)= 0.
               endif
            enddo
         enddo
      end do
   !>  - Calculate saturation specific humidity and enforce minimum moisture values.
      !$acc parallel loop gang collapse(2) async(async_id) private(val1, val2, tem)
      do jj = 1, jlistnum
         do k = 1, km
            !$acc loop vector
            do i = 1, myim(jj)
               if (k <= kmax(i, jj)) then
                 
                  qesor = 0.01 * fpvs_gpu(to(i, k, jj),c1xpvs,c2xpvs,tbpvs)      ! fpvs is in pa
   !byl               call qsatq_cwb(1,to(i,k),pfld(i,k),qeso(i,k))
                  qesor = eps * qesor / (pfld(i, k, jj) + epsm1*qesor)
                  val1      =             1.e-8
                  qesor = max(qesor, val1)
                  val2      =           1.e-10
                  qo(i, k, jj)   = max(qo(i, k, jj), val2 )
   !              qo(i,k)   = min(qo(i,k),qeso(i,k))
   !              tvo(i,k)  = to(i,k) + delta * to(i,k) * qo(i,k)
   !
   !  compute moist static energy
   !
   !>  - Calculate moist static energy (heo) and saturation moist static energy (heso).
   !              tem       = g * zo(i,k) + cp * to(i,k)
                  tem       = phil(i, k, jj) + cp * to(i, k, jj)
                  heo(i, k, jj)  = tem  + hvap * qo(i, k, jj)
                  heso(i, k, jj) = tem  + hvap * qesor
                  qeso(i, k, jj) = qesor
   !              heo(i,k)  = min(heo(i,k),heso(i,k))
               endif
            enddo
         enddo
      end do
      val1      =             1.e-8
      val2      =           1.e-10
      !$acc parallel loop gang vector collapse(2) async(async_id) private(dz, dp, &
      !$acc&         es, pprime, qs, dqsdp, desdt, dqsdt, gamma, dt, dq, val1, &
      !$acc&         val2, qesor, qor1, tem, tor, kbr, hmaxr, kmaxr, zor1, zor, &
      !$acc&         pfldr1, pfldr, tor1, uor, uor1, vor, vor1, po)
      do jj = 1, jlistnum
         do i = 1, ix
            if (i .le. myim(jj)) then
   !
   !  determine level with largest moist static energy
   !  this is the level where updraft starts
   !
   !> ## Perform calculations related to the updraft of the entraining/detraining cloud model ("static control").
   !> - Search below index "kbm" for the level of maximum moist static energy.
   !
               hmaxr = heo(i, 1, jj)
               kbr   = 1
               
               !$acc loop seq
               do k = 2, km
                  if (k <= kbm(i, jj)) then
                     zor = heo(i, k, jj)
                     if(zor > hmaxr) then
                        kbr   = k
                        hmaxr = zor
                     endif
                  endif
               enddo
               kb(i, jj) = kbr
               hmax(i, jj) = hmaxr
              
   !> - Calculate the temperature, specific humidity, and pressure at interface levels.
               kmaxr = kmax(i, jj)
               zor = zo(i, 1, jj)
               uor = uo(i, 1, jj)
               vor = vo(i, 1, jj)
               pfldr = pfld(i, 1, jj)
               !$acc loop seq
               do k = 1, km1
                  if (k <= kmaxr-1) then
                     zor1 = zo(i, k+1, jj)
                     pfldr1 = pfld(i, k+1, jj)
                     tor1 = to(i, k+1, jj)
                     qor1 = qo(i, k+1, jj)
                     uor1 = uo(i, k+1, jj)
                     vor1 = vo(i, k+1, jj)
                     dz      = .5 * (zor1 - zor)
                     dp      = .5 * (pfldr1 - pfldr)
                     es      = 0.01 * fpvs_gpu(tor1,c1xpvs,c2xpvs,tbpvs)      ! fpvs is in pa
      !byl               call qsatq_cwb(1,to(i,k+1),pfld(i,k+1),es)
                     pprime  = pfldr1 + epsm1 * es
                     qs      = eps * es / pprime
                     dqsdp   = - qs / pprime
                     desdt   = es * (fact1 / tor1 + fact2 / (tor1**2))
                     dqsdt   = qs * pfldr1 * desdt / (es * pprime)
                     gamma   = el2orc * qeso(i, k+1, jj) / (tor1**2)
                     dt      = (g * dz + hvap * dqsdp * dp) / (cp * (1. + gamma))
                     dq      = dqsdt * dt + dqsdp * dp
                     tor1 = tor1 + dt
                     to(i, k, jj) = tor1
                     qor1 = qor1 + dq
                     po = .5 * (pfldr + pfldr1)
   !
   !> - Recalculate saturation specific humidity, moist static energy, saturation moist static energy, and horizontal momentum on interface levels. Enforce minimum specific humidity and calculate \f$(1 - RH)\f$.
                     qesor = 0.01 * fpvs_gpu(tor1,c1xpvs,c2xpvs,tbpvs)      ! fpvs is in pa
      !byl               call qsatq_cwb(1,to(i,k),po(i,k),qeso(i,k))
                     qesor = eps * qesor / (po + epsm1*qesor)
                     qesor = max(qesor, val1)
                     qor1   = max(qor1, val2 )
                     qo(i, k, jj) = qor1
                     tem = .5 * g * (zor + zor1) + &
                                 cp * tor1 
      !              qo(i,k)   = min(qo(i,k),qeso(i,k))
                     frh(i, k, jj)  = 1. - min(qor1/qesor, 1.)
                     heo(i, k, jj)  = tem + hvap * qor1
                     heso(i, k, jj) = tem + hvap * qesor
                     qeso(i, k, jj) = qesor
                     uo(i, k, jj)   = .5 * (uor + uor1)
                     vo(i, k, jj)   = .5 * (vor + vor1)
                     zor = zor1
                     pfldr = pfldr1
                     uor = uor1
                     vor = vor1
                  endif
               enddo
            end if
         enddo
      end do
      
      !$acc parallel loop gang vector collapse(2) async(async_id) private(flgr, &
      !$acc&         sigma_con, cinpcr, tem1, tem, ptem, dz, factor, ptem1, &
      !$acc&         kmaxr, kbconr, kbr, kbmaxr, kbcon1r, fent1, fent2, xlamx, &
      !$acc&         cinpcri, etar1, zir1, zir, xlamudr1, xlamudr, xlamuer1, xlamuer, &
      !$acc&         hckor1, heor1, uckor1, uor1, vckor1, vor1, uor, vor, heor, heor1, &
      !$acc&         frh_sumr)
      do jj = 1, jlistnum
         do i = 1, ix
            if (i .le. myim(jj)) then
   !
   !  look for the level of free convection as cloud base
   !
   !> - Search below the index "kbmax" for the level of free convection (LFC) where the condition \f$h_b > h^*\f$ is first met, where \f$h_b, h^*\f$ are the state moist static energy at the parcel's starting level and saturation moist static energy, respecti
               flgr   = .true.
               kmaxr = kmax(i, jj)
               kbconr = kmaxr
               kbr = kb(i, jj)
               kbmaxr = kbmax(i, jj)
               !$acc loop seq
               do k = 1, km1
                  if (flgr .and. k <= kbmaxr) then
                     if(k > kbr .and. heo(i, kbr, jj) > heso(i, k, jj)) then
                        kbconr = k
                        flgr   = .false.
                     endif
                  endif
               enddo
   !
   !> - If no LFC, return to the calling routine without modifying state variables.
               if(kbconr == kmaxr) cnvflg(i, jj) = .false.
   !!
   !!
   !> - Determine the vertical pressure velocity at the LFC. After Han and Pan (2011) \cite han_and_pan_2011 , determine the maximum pressure thickness between a parcel's starting level and the LFC. If a parcel doesn't reach the LFC within the critical thick
               if(cnvflg(i, jj)) then
                  pdot(i, jj)  = 10.* dot(i, kbconr, jj)
      !           pdot(i)  = 0.01 * dot(i,kbcon(i)) ! Now dot is in Pa/s
                  frh_sumr = frh_sum(i, jj)
                  !$acc loop seq
                  do k = 1, km1
                     if (k.ge.kbr.and.k.le.kbconr) then
                        ! frh_sum can use reduction
                        frh_sumr = frh_sumr + (1-frh(i, k, jj))
                     endif
                  enddo

   !
   !   turn off convection if pressure depth between parcel source level
   !      and cloud base is larger than a critical value, cinpcr
   !
                  sigma_con = tan(0.4*pi)/(dx5km-dx1km)           !7.7e-4 m-1
                  sigma(i, jj)  = (1.-1./pi*(atan(sigma_con*(gdx(i, jj)-dx5km))+pi/2.)) !1(1km),0.1(10km)
                  if (gdx(i, jj).lt.dx5km) then
                     sigma(i, jj) = min(sigma(i, jj)-0.01684*gdx(i, jj)/1000.+0.0842,1.0)
                  endif
                  cinpcr = cinpcrmn + 0.5*(cinpcrmx-cinpcrmn) * (1.-sigma(i, jj))
                  cinpcri = cinpcr * frh_sumr/(kbconr-kbr+1)
                  tem1 = pfld(i, kbr, jj) - pfld(i, kbconr, jj)
                  if(tem1 > cinpcri) then
                     cnvflg(i, jj) = .false.
                  endif
                  frh_sum(i, jj) = frh_sumr
               endif
   !!
   !!
   !
   !  assume that updraft entrainment rate above cloud base is
   !    same as that at cloud base
   !
   !> - Calculate the entrainment rate according to Han and Pan (2011) \cite han_and_pan_2011 , equation 8, after Bechtold et al. (2008) \cite bechtold_et_al_2008, equation 2 given by:
   !!  \f[
   !!  \epsilon = \epsilon_0F_0 + d_1\left(1-RH\right)F_1
   !!  \f]
   !!  where \f$\epsilon_0\f$ is the cloud base entrainment rate, \f$d_1\f$ is a tunable constant, and \f$F_0=\left(\frac{q_s}{q_{s,b}}\right)^2\f$ and \f$F_1=\left(\frac{q_s}{q_{s,b}}\right)^3\f$ where \f$q_s\f$ and \f$q_{s,b}\f$ are the saturation specific
               if(cnvflg(i, jj)) then
                  xlamx = xlamue(i, kbconr, jj)
                  if(1 < kmaxr) then
                     xlamud(i, 1, jj) = xlamx
                  end if
               
                  !$acc loop seq
                  do k = 2, km1
      !
      !  specify detrainment rate for the updrafts
      !
      !> - The updraft detrainment rate is set constant and equal to the entrainment rate at cloud base.
                     if(k < kmaxr) then
                        xlamud(i, k, jj) = xlamx
         !              xlamud(i,k) = crtlamd
                     endif
                     
                     if(k > kbconr .and. k < kmaxr) then
                          xlamue(i, k, jj) = xlamx
      !
      !  functions rapidly decreasing with height, mimicking a cloud ensemble
      !    (Bechtold et al., 2008)
      !
                        tem = qeso(i, k, jj)/qeso(i, kbconr, jj)
                        fent1 = tem**2
                        fent2 = tem**3
      !
      !  final entrainment and detrainment rates as the sum of turbulent part and
      !    organized entrainment depending on the environmental relative humidity
      !    (Bechtold et al., 2008)
      !
                        tem = cxlamu * frh(i, k, jj) * fent2
                        xlamue(i, k, jj) = xlamue(i, k, jj)*fent1 + tem
         !              tem1 = cxlamd * frh(i,k)
         !              xlamud(i,k) = xlamud(i,k) + tem1
                     endif
                  enddo
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   !  determine updraft mass flux for the subcloud layers
   !
   !> - Calculate the normalized mass flux for subcloud and in-cloud layers according to Pan and Wu (1995) \cite pan_and_wu_1995 equation 1:
   !!  \f[
   !!  \frac{1}{\eta}\frac{\partial \eta}{\partial z} = \lambda_e - \lambda_d
   !!  \f]
   !!  where \f$\eta\f$ is the normalized mass flux, \f$\lambda_e\f$ is the entrainment rate and \f$\lambda_d\f$ is the detrainment rate.
                  
                  zir1 = zi(i, kbconr, jj)
                  xlamudr1 = xlamud(i, kbconr, jj)
                  xlamuer1 = xlamue(i, kbconr, jj)
                  etar1 = eta(i, kbconr, jj)
                  !$acc loop seq
                  do k = km1, 1, -1
                     if(k < kbconr .and. k >= kbr) then
                        zir = zi(i, k, jj)
                        xlamudr = xlamud(i, k, jj)
                        xlamuer = xlamue(i, k, jj)
                        dz       = zir1 - zir
                        tem      = 0.5*(xlamudr+xlamudr1)
                        ptem     = 0.5*(xlamuer+xlamuer1)-tem
                        etar1 = etar1 / (1. + ptem * dz)
                        eta(i, k, jj) = etar1
                        zir1 = zir
                        xlamudr1 = xlamudr
                        xlamuer1 = xlamuer
                     endif
                  enddo
      !
      !  compute mass flux above cloud base
      !
                  flgr = .true.
                  zir1 = zi(i, kbconr, jj)
                  xlamudr1 = xlamud(i, kbconr, jj)
                  xlamuer1 = xlamue(i, kbconr, jj)
                  etar1 = eta(i, kbconr, jj)
                  !$acc loop seq
                  do k = 2, km1
                     if(flgr)then
                        if(k > kbconr .and. k < kmaxr) then
                           zir = zi(i, k, jj)
                           xlamudr = xlamud(i, k, jj)
                           xlamuer = xlamue(i, k, jj)
                           dz       = zir - zir1
                           tem      = 0.5*(xlamudr+xlamudr1)
                           ptem     = 0.5*(xlamuer+xlamuer1)-tem
                           etar1 = etar1 * (1 + ptem * dz)
                           eta(i, k, jj) = etar1
                           zir1 = zir
                           xlamudr1 = xlamudr
                           xlamuer1 = xlamuer
                           if(etar1 <= 0.) then
                              kmaxr = k
                              ktconn(i, jj) = k
                              flgr   = .false.
                           endif
                        endif
                     endif
                  enddo
   !
   !  compute updraft cloud properties
   !
   !> - Set cloud properties equal to the state variables at updraft starting level (kb).
                  heor1 = heo(i, kbr, jj)
                  uor1 = uo(i, kbr, jj)
                  vor1 = vo(i, kbr, jj)
                  hcko(i, kbr, jj) = heor1
                  ucko(i, kbr, jj) = uor1
                  vcko(i, kbr, jj) = vor1
                  hckor1 = heor1
                  uckor1 = uor1
                  vckor1 = vor1
   !
   !  cloud property is modified by the entrainment process
   !
   !  cm is an enhancement factor in entrainment rates for momentum
   !
   !> - Calculate the cloud properties as a parcel ascends, modified by entrainment and detrainment. Discretization follows Appendix B of Grell (1993) \cite grell_1993 . Following Han and Pan (2006) \cite han_and_pan_2006, the convective momentum transport i
                  
                  zir1 = zi(i, kbr, jj)
                  xlamudr1 = xlamud(i, kbr, jj)
                  xlamuer1 = xlamue(i, kbr, jj)
                  !$acc loop seq
                  do k = 2, km1
                     if(k > kbr .and. k < kmaxr) then
                        heor = heo(i, k, jj)
                        zir = zi(i, k, jj)
                        xlamudr = xlamud(i, k, jj)
                        xlamuer = xlamue(i, k, jj)
                        uor = uo(i, k, jj)
                        vor = vo(i, k, jj)
                        dz   = zir - zir1
                        tem  = 0.5 * (xlamuer+xlamuer1) * dz
                        tem1 = 0.25 * (xlamudr+xlamudr1) * dz
                        factor = 1. + tem - tem1
                        hckor1 = ((1.-tem1)*hckor1+tem*0.5* &
                                     (heor+heor1))/factor
                        dbyo(i, k, jj) = hckor1 - heso(i, k, jj)
                        hcko(i, k, jj) = hckor1
      !         
                        tem  = 0.5 * cm * tem
                        factor = 1. + tem
                        ptem = tem + pgcon
                        ptem1= tem - pgcon
                        uckor1 = ((1.-tem)*uckor1+ptem*uor &
                                     +ptem1*uor1)/factor
                        vckor1 = ((1.-tem)*vckor1+ptem*vor &
                                     +ptem1*vor1)/factor
                        ucko(i, k, jj) = uckor1
                        vcko(i, k, jj) = vckor1
                        zir1 = zir
                        xlamudr1 = xlamudr
                        xlamuer1 = xlamuer
                        uor1 = uor
                        vor1 = vor
                        heor1 = heor
                     endif
                  enddo
   !
   !   taking account into convection inhibition due to existence of
   !    dry layers below cloud base
   !
   !> - With entrainment, recalculate the LFC as the first level where buoyancy is positive. The difference in pressure levels between LFCs calculated with/without entrainment must be less than a threshold (currently 25 hPa). Otherwise, convection is inhibit
                  flgr = .true.
                  kbcon1r = kmaxr
                  
                  !$acc loop seq
                  do k = 2, km1
                     if (flgr .and. k < kmaxr) then
                        if(k >= kbconr .and. dbyo(i, k, jj) > 0.) then
                           kbcon1r = k
                           flgr    = .false.
                        endif
                     endif
                  enddo
               
                  if(kbcon1r == kmaxr) cnvflg(i, jj) = .false.
               endif
               
               if(cnvflg(i, jj)) then
                  tem = pfld(i, kbconr, jj) - pfld(i, kbcon1r, jj)
                  if(tem > dthk) then
                     cnvflg(i, jj) = .false.
                  endif
               endif
               kmax(i, jj) = kmaxr
               kbcon(i, jj) = kbconr
               kbcon1(i, jj) = kbcon1r
            end if
         enddo
      end do
   !!
   !!
   !
   !  calculate convective inhibition
   !
   !> - Calculate additional trigger condition of the convective inhibition (CIN) according to Han et al.'s (2017) \cite han_et_al_2017 equation 13.
      ! cina can use reduction
      !$acc parallel loop gang collapse(2) async(async_id) private(dz1, gamma, &
      !$acc&         rfact, val, cinar, tor)
      do jj = 1, jlistnum
         do k = 2, km1
            !$acc loop vector
            do i = 1, myim(jj)
               if (cnvflg(i, jj)) then
                  if(k > kb(i, jj) .and. k < kbcon1(i, jj)) then
                     tor = to(i, k, jj)
                     dz1 = zo(i, k+1, jj) - zo(i, k, jj)
                     gamma = el2orc * qeso(i, k, jj) / (tor**2)
                     rfact =  1. + delta * cp * gamma &
                              * tor / hvap
                     cinar =                  &
   !    &                     dz1 * eta(i,k) * (g / (cp * to(i,k))) &
                              dz1 * (g / (cp * tor)) &
                              * dbyo(i, k, jj) / (1. + gamma) &
                              * rfact
                     val = 0.
                     cinar = cinar +                      &
   !    &                     dz1 * eta(i,k) * g * delta * &
                              dz1 * g * delta *              &
                              max(val,(qeso(i, k, jj) - qo(i, k, jj)))
                     !$acc atomic
                     cina(i, jj) = cina(i, jj) + cinar
                  endif
               endif
            enddo
         enddo
      end do
   !> - Turn off convection if the CIN is less than a critical value (cinacr) which is inversely proportional to the large-scale vertical velocity.
      !$acc parallel loop gang vector collapse(2) async(async_id) private(w1, &
      !$acc&         w2, w3, w4, tem, val1, val2, tem1, cinacr, k ,dp, dz, &
      !$acc&         gamma, qrch, factor, dq, etah, ptem, qlk, rfact, val, dz1, &
      !$acc&         ptem1, flgr, pdotr, ktconr, kmaxr, kbcon1r, ktconnr, &
      !$acc&         kbconr, kbr, kbmaxr, hmin, lmin, jminr, ktcon1r, heor, &
      !$acc&         aa1r, drag, dragm1, wu2r1, wu2r, zir1, zir, buor1, buor, &
      !$acc&         qckor1, qor, qor1, qesor, tor, kk, &
      !$acc&         xlamuer1, xlamuer, xlamudr1, xlamudr, etar, etar1, c0tr)
      do jj = 1, jlistnum
         do i = 1, ix
            if (i .le. myim(jj)) then
               if(cnvflg(i, jj)) then
      !
                  if(islimsk(i, jj) == 1) then
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
                  pdotr = pdot(i, jj)
                  if(pdotr <= w4) then
                     tem = (pdotr - w4) / (w3 - w4)
                  elseif(pdotr >= -w4) then
                     tem = - (pdotr + w4) / (w4 - w3)
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
      !           cinacr = cinacrmx
                  if(cina(i, jj) < cinacr) cnvflg(i, jj) = .false.
               endif
   !!
   !!
   !
   !  determine first guess cloud top as the level of zero buoyancy
   !
   !> - Calculate the cloud top as the first level where parcel buoyancy becomes negative. If the thickness of the calculated convection is less than a threshold (currently 200 hPa), then convection is inhibited, and the scheme returns to the calling routine
               ktconr = 1
               kmaxr = kmax(i, jj)
               kbcon1r = kbcon1(i, jj)
               ktconnr = ktconn(i, jj)
               kbconr = kbcon(i, jj)
               kbr = kb(i, jj)
               kbmaxr = kbmax(i, jj)
               pwavor = 0.
               if(cnvflg(i, jj)) then
                  flgr = .true.

                  !$acc loop seq
                  do k = 2, km1
                     if (flgr .and. k < kmaxr) then
                        if(k > kbcon1r .and. dbyo(i, k, jj) < 0.) then
                           ktconr = k
                           flgr   = .false.
                        endif
                     endif
                  enddo
   !
                  if(ktconr == 1 .and. ktconnr > 1) then
                     ktconr = ktconnr
                  endif
                  tem = pfld(i, kbconr, jj)-pfld(i, ktconr, jj)
                  if(tem < cthk) cnvflg(i, jj) = .false.
               endif
   !!
   !
   !  search for downdraft originating level above theta-e minimum
   !
   !> - To originate the downdraft, search for the level above the minimum in moist static energy. Return to the calling routine without modification if this level is determined to be outside of the convective cloud layers.
               if(cnvflg(i, jj)) then
                  hmin = heo(i, kbcon1r, jj)
                  lmin = 1
                  lmin = kbmaxr
                  jminr = kbmaxr
                  
                  !$acc loop seq
                  do k = 2, km1
                     if (k <= kbmaxr) then
                        heor = heo(i, k, jj)
                        if(k > kbcon1r .and. heor < hmin) then
                           lmin = k + 1
                           hmin = heor
                        endif
                     endif
                  enddo
   !
   !  make sure that jmin is within the cloud
   !
                  jminr = min(lmin, ktconr-1)
                  jminr = max(jminr, kbcon1r+1)
                  if(jminr >= ktconr) cnvflg(i, jj) = .false.
                  jmin(i, jj) = jminr
               endif
   !
   !  specify upper limit of mass flux at cloud base
   !
   !> - Calculate the maximum value of the cloud base mass flux using the CFL-criterion-based formula of Han and Pan (2011) \cite han_and_pan_2011, equation 7.
               if(cnvflg(i, jj)) then
      !           xmbmax(i) = .1
      !         
                  k = kbconr
                  dp = 1000. * del(i, k, jj)
                  xmbmax(i, jj) = dp / (g * dt2)
      !         
      !           mbdt(i) = 0.1 * dp / g
      !         
      !           tem = dp / (g * dt2)
      !           xmbmax(i) = min(tem, xmbmax(i))
   !
   !  compute cloud moisture property and precipitation
   !
   !> - Set cloud moisture property equal to the enviromental moisture at updraft starting level (kb).
      !           aa1(i) = 0.
                  qckor1 = qo(i, kbr, jj)
                  qrcko(i, kbr, jj) = qo(i, kbr, jj)
                  qcko(i, kbr, jj) = qckor1
      !           rhbar(i) = 0.
   !> - Calculate the moisture content of the entraining/detraining parcel (qcko) and the value it would have if just saturated (qrch), according to equation A.14 in Grell (1993) \cite grell_1993 . Their difference is the amount of convective cloud water (ql
                  jminr = jmin(i, jj)
                  zir1 = zi(i, kbr, jj)
                  xlamuer1 = xlamue(i, kbr, jj)
                  xlamudr1 = xlamud(i, kbr, jj)
                  qor1 = qo(i, kbr, jj)
                  etar1 = eta(i, kbr, jj)
                  !$acc loop seq
                  do k = 2, km1
                     if(k > kbr .and. k < ktconr) then
                        qesor = qeso(i, k, jj)
                        tor = to(i, k, jj)
                        zir = zi(i, k, jj)
                        xlamuer = xlamue(i, k, jj)
                        xlamudr = xlamud(i, k, jj)
                        qor = qo(i, k, jj)
                        etar = eta(i, k, jj)
                        c0tr = c0t(i, k, jj)
                        dz    = zir - zir1
                        gamma = el2orc * qesor / (tor**2)
                        qrch = qesor &
                             + gamma * dbyo(i, k, jj) / (hvap * (1. + gamma))
      !cj         
                        tem  = 0.5 * (xlamuer+xlamuer1) * dz
                        tem1 = 0.25 * (xlamudr+xlamudr1) * dz
                        factor = 1. + tem - tem1
                        qckor1 = ((1.-tem1)*qckor1+tem*0.5* &
                                     (qor+qor1))/factor
                        qrcko(i, k, jj) = qckor1
      !cj         
                        dq = etar * (qckor1 - qrch)
      !           
      !                 rhbar(i) = rhbar(i) + qo(i,k) / qeso(i,k)
      !
      !  check if there is excess moisture to release latent heat
      !
                        buor = buo(i, k, jj)
                        if(k >= kbconr .and. dq > 0.) then
                           etah = .5 * (etar + etar1)
                           dp = 1000. * del(i, k, jj)
                           if(ncloud > 0 .and. k > jminr) then
                              ptem = c0tr + c1
                              qlk = dq / (etar + etah * ptem * dz)
                              dellal(i, k, jj) = etah * c1 * dz * qlk * g / dp
                           else
                              qlk = dq / (etar + etah * c0tr * dz)
                           endif
      !                    aa1(i) = aa1(i) - dz * g * qlk * etah
      !                    aa1(i) = aa1(i) - dz * g * qlk
                           buor = buor - g * qlk
                           qckor1 = qlk + qrch
                           pwo(i, k, jj) = etah * c0tr * dz * qlk
                           pwavor = pwavor + pwo(i, k, jj)
      !                    cnvwt(i,k) = (etah*qlk + pwo(i,k)) * g / dp
                           cnvwt(i, k, jj) = etah * qlk * g / dp
                        endif
      !
      !  compute buoyancy and drag for updraft velocity
      !
                        if(k >= kbconr) then
                           rfact =  1. + delta * cp * gamma &
                                    * tor / hvap
                           buor = buor + (g / (cp * tor)) &
                                    * dbyo(i, k, jj) / (1. + gamma) &
                                    * rfact
                           val = 0.
                           buor = buor + g * delta * &
                                      max(val,(qesor - qor))
                        endif
                        zir1 = zir
                        xlamuer1 = xlamuer
                        xlamudr1 = xlamudr
                        qor1 = qor
                        etar1 = etar
                        buo(i, k, jj) = buor
                        qcko(i, k, jj) = qckor1
                     endif
                  enddo
   !
   !     do i = 1, im
   !       if(cnvflg(i)) then
   !         indx = ktcon(i) - kb(i) - 1
   !         rhbar(i) = rhbar(i) / float(indx)
   !       endif
   !     enddo
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
   !
   !  calculate cloud work function
   !
   !> - Calculate the cloud work function according to Pan and Wu (1995) \cite pan_and_wu_1995 equation 4:
   !!  \f[
   !!  A_u=\int_{z_0}^{z_t}\frac{g}{c_pT(z)}\frac{\eta}{1 + \gamma}[h(z)-h^*(z)]dz
   !!  \f]
   !! (discretized according to Grell (1993) \cite grell_1993 equation B.10 using B.2 and B.3 of Arakawa and Schubert (1974) \cite arakawa_and_schubert_1974 and assuming \f$\eta=1\f$) where \f$A_u\f$ is the updraft cloud work function, \f$z_0\f$ and \f$z_t\f
                  aa1r = 0.
                  ! aa1 can use reduction
                  kk = max(2, kbconr)
                  kk = min(kk, km1)
                  zir = zo(i, kk, jj)
                  !$acc loop seq
                  do k = 2, km1
                     if(k >= kbconr .and. k < ktconr) then
                        zir1 = zo(i, k+1, jj)
                        dz1 = zir1 - zir
      !                 aa1(i) = aa1(i) + buo(i,k) * dz1 * eta(i,k)
                        aa1r = aa1r + buo(i, k, jj) * dz1
                        zir = zir1
                     endif
                  enddo
                  aa1(i, jj) = aa1r
   !
   !> - If the updraft cloud work function is negative, convection does not occur, and the scheme returns to the calling routine.
               if(aa1r <= 0.) cnvflg(i, jj) = .false.
            endif

   !!
   !
   !  estimate the onvective overshooting as the level
   !    where the [aafac * cloud work function] becomes zero,
   !    which is the final cloud top
   !
   !> - Continue calculating the cloud work function past the point of neutral buoyancy to represent overshooting according to Han and Pan (2011) \cite han_and_pan_2011 . Convective overshooting stops when \f$ cA_u < 0\f$ where \f$c\f$ is currently 10%, or w
               ktcon1r = kmaxr
               if (cnvflg(i, jj)) then
                  aa2 = aafac * aa1(i, jj)
   !
                  flgr = .true.
                  kk = max(2, ktconr)
                  kk = min(kk, km1)
                  zir = zo(i, kk, jj)
                  !$acc loop seq
                  do k = 2, km1
                     if (flgr) then
                        if(k >= ktconr .and. k < kmaxr) then
                           tor = to(i, k, jj)
                           zir1 = zo(i, k+1, jj)
                           dz1 = zir1 - zir
                           gamma = el2orc * qeso(i, k, jj) / (tor**2)
                           rfact =  1. + delta * cp * gamma &
                                    * tor / hvap
                           aa2 = aa2 +                   &
         !    &                     dz1 * eta(i,k) * (g / (cp * to(i,k))) &
                                    dz1 * (g / (cp * tor)) &
                                    * dbyo(i, k, jj) / (1. + gamma) &
                                    * rfact
         !                 val = 0.
         !                 aa2(i) = aa2(i) +
         !!   &                     dz1 * eta(i,k) * g * delta *
         !    &                     dz1 * g * delta *
         !    &                     max(val,(qeso(i,k) - qo(i,k)))
                           if(aa2 < 0.) then
                              ktcon1r = k
                              flgr = .false.
                           endif
                           zir = zir1
                        endif
                     endif
                  enddo
   !
   !  compute cloud moisture property, detraining cloud water
   !    and precipitation in overshooting layers
   !
   !> - For the overshooting convection, calculate the moisture content of the entraining/detraining parcel as before. Partition convective cloud water and precipitation and detrain convective cloud water above the mimimum in moist static energy.
                  kk = max(2, ktconr)
                  kk = min(kk, km1)
                  zir1 = zi(i, kk-1, jj)
                  xlamuer1 = xlamue(i, kk-1, jj)
                  xlamudr1 = xlamud(i, kk-1, jj)
                  etar1 = eta(i, kk-1, jj)
                  qor1 = qo(i, kk-1, jj)
                  qckor1 = qcko(i, kk-1, jj)
                  !$acc loop seq
                  do k = 2, km1
                     if(k >= ktconr .and. k < ktcon1r) then
                        zir = zi(i, k, jj)
                        xlamuer = xlamue(i, k, jj)
                        xlamudr = xlamud(i, k, jj)
                        etar = eta(i, k, jj)
                        qor = qo(i, k, jj)
                        qesor = qeso(i, k, jj)
                        dz    = zir - zir1
                        gamma = el2orc * qesor / (to(i, k, jj)**2)
                        qrch = qesor &
                             + gamma * dbyo(i, k, jj) / (hvap * (1. + gamma))
      !cj         
                        tem  = 0.5 * (xlamuer+xlamuer1) * dz
                        tem1 = 0.25 * (xlamudr+xlamudr1) * dz
                        factor = 1. + tem - tem1
                        qckor1 = ((1.-tem1)*qckor1+tem*0.5* &
                                     (qor+qor1))/factor
                        qrcko(i, k, jj) = qckor1
      !cj         
                        dq = etar * (qckor1 - qrch)
      !
      !  check if there is excess moisture to release latent heat
      !
                        if(dq > 0.) then
                           c0tr = c0t(i, k, jj)
                           etah = .5 * (etar + etar1)
                           dp = 1000. * del(i, k, jj)
                           if(ncloud > 0) then
                             ptem = c0tr + c1
                             qlk = dq / (etar + etah * ptem * dz)
                             dellal(i, k, jj) = etah * c1 * dz * qlk * g / dp
                           else
                             qlk = dq / (etar + etah * c0tr * dz)
                           endif
                           qckor1 = qlk + qrch
                           pwo(i, k, jj) = etah * c0tr * dz * qlk
                           pwavor = pwavor + pwo(i, k, jj)
      !                    cnvwt(i,k) = (etah*qlk + pwo(i,k)) * g / dp
                           cnvwt(i, k, jj) = etah * qlk * g / dp
                        endif
                        zir1 = zir
                        xlamuer1 = xlamuer
                        xlamudr1 = xlamudr
                        qor1 = qor
                        etar1 = etar
                        qcko(i, k, jj) = qckor1
                     endif
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
   !      bb1 = 2.0
   !      bb2 = 4.0
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
                  zir1 = zi(i, kbcon1r, jj)
                  buor1 = buo(i, kbcon1r, jj)
                  wu2r1 = wu2(i, kbcon1r, jj)
                  dragm1 = 0.
                  !$acc loop seq
                  do k = 2, km1
                     drag = 0.
                     if (k > kbr .and. k < ktconr) then
                        if (k >= kbconr) then
                           drag = max(xlamue(i, k, jj),xlamud(i, k, jj))
                        end if
                     end if
                     if(k > kbcon1r .and. k < ktconr) then
                        zir = zi(i, k, jj)
                        buor = buo(i, k, jj)
                        dz    = zir - zir1
                        tem  = 0.25 * bb1 * (drag+dragm1) * dz
                        tem1 = 0.5 * bb2 * (buor+buor1) * dz
                        ptem = (1. - tem) * wu2r1
                        ptem1 = 1. + tem
                        wu2r = (ptem + tem1) / ptem1
                        wu2r = max(wu2r, 0.)
                        wu2(i, k, jj) = wu2r
                        wu2r1 = wu2r
                        zir1 = zir
                        buor1 = buor
                     endif
                     dragm1 = drag
                  enddo
               end if

               ktcon(i, jj) = ktconr
               ktcon1(i, jj) = ktcon1r
               pwavo(i, jj) = pwavor
            end if
         enddo
      end do

      !$acc parallel loop gang collapse(2) async(async_id) private(dz, tem, &
      !$acc&         tem1, val, kk, k, gamma, qrch, dq, shear, e1, beta, &
      !$acc&         factor, ptem, ptem1, edtmax, po1, po1m1, sumx, wcr, wbarr, &
      !$acc&         vshear)
      do jj = 1, jlistnum
         do i = 1, ix
            if (i .le. myim(jj)) then
   !
   !  compute updraft velocity average over the whole cumulus
   !
   !> - Calculate the mean updraft velocity within the cloud (wc).
               wcr = 0.
               sumx = 0.
               wbarr = 0.
   !        ptem = -0.5 * rd/g
               if (cnvflg(i, jj)) then
                  !$acc loop seq
                  do k = 2, km1
                     if(k > kbcon1(i, jj) .and. k < ktcon(i, jj)) then
                        po1 = .5 * (pfld(i, k, jj) + pfld(i, k+1, jj))
                        po1m1 = .5 * (pfld(i, k-1, jj) + pfld(i, k, jj))
                        dz = zi(i, k, jj) - zi(i, k-1, jj)
                        tem = 0.5 * (sqrt(wu2(i, k, jj)) + sqrt(wu2(i, k-1, jj)))
                        wcr = wcr + tem * dz
                        tem=10.0*dot(i, k, jj)*to(i, k, jj)/po1
                        tem1=10.0*dot(i, k-1, jj)*to(i, k-1, jj)/po1m1
                        wbarr=wbarr+(-0.5*rd/g)*(tem+tem1)*dz 
                        sumx = sumx + dz
                     endif
                  enddo
               
                  if(sumx == 0.) then
                     cnvflg(i, jj)=.false.
                  else
                     wcr = wcr / sumx
                  endif
                  val = 1.e-4
                  if (wcr < val) cnvflg(i, jj)=.false.
               endif
               
               if(cnvflg(i, jj)) then
                  if(sumx > 0.) then
                     wbarr = wbarr / sumx   
                  endif
   ! compute mean cloud core fraction
   ! assume mean cloud core fraction to be the ratio of
   ! mean grid-scale vertical velocity (wbar) and mean updraft velocity
      
                  tem = wbarr / wcr
                  tem = max(tem, 0.)
                  clear(i, jj) = 1. - tem
                  clear(i, jj) = max(min(clear(i, jj), 1.0), 0.)
                  if (wbarr.gt.0. .and. wbarr.gt.wcr) then
                     cnvflg(i, jj) = .false.
                  endif
               endif
   !
   ! exchange ktcon with ktcon1
   !
   !> - Swap the indices of the convective cloud top (ktcon) and the overshooting convection top (ktcon1) to use the same cloud top level in the calculations of \f$A^+\f$ and \f$A^*\f$.
               if(cnvflg(i, jj)) then
                  kk = ktcon(i, jj)
                  ktcon(i, jj) = ktcon1(i, jj)
                  ktcon1(i, jj) = kk
   !
   !  this section is ready for cloud water
   !
   !> - Separate the total updraft cloud water at cloud top into vapor and condensate.
                  if(ncloud > 0) then
   !
   !  compute liquid and vapor separation at cloud top
   !
                     k = ktcon(i, jj) - 1
                     gamma = el2orc * qeso(i, k, jj) / (to(i, k, jj)**2)
                     qrch = qeso(i, k, jj) &
                          + gamma * dbyo(i, k, jj) / (hvap * (1. + gamma))
                     dq = qcko(i, k, jj) - qrch
         !
         !  check if there is excess moisture to release latent heat
         !
                     if(dq > 0.) then
                        qlko_ktcon(i, jj) = dq
                        qcko(i, k, jj) = qrch
                     endif
                  endif
   !
   !ccccc if(lat.==.latd.and.lon.==.lond.and.cnvflg(i)) then
   !ccccc   print *, ' aa1(i) before dwndrft =', aa1(i)
   !ccccc endif
   !
   !c------- downdraft calculations
   !
   !c--- compute precipitation efficiency in terms of windshear
   !
   !> ## Perform calculations related to the downdraft of the entraining/detraining cloud model ("static control").
   !! - First, in order to calculate the downdraft mass flux (as a fraction of the updraft mass flux), calculate the wind shear and precipitation efficiency according to equation 58 in Fritsch and Chappell (1980) \cite fritsch_and_chappell_1980 :
   !! \f[
   !! E = 1.591 - 0.639\frac{\Delta V}{\Delta z} + 0.0953\left(\frac{\Delta V}{\Delta z}\right)^2 - 0.00496\left(\frac{\Delta V}{\Delta z}\right)^3
   !! \f]
   !! where \f$\Delta V\f$ is the integrated horizontal shear over the cloud depth, \f$\Delta z\f$, (the ratio is converted to units of \f$10^{-3} s^{-1}\f$). The variable "edto" is \f$1-E\f$ and is constrained to the range \f$[0,0.9]\f$.
                  vshear = 0.
      ! vshear can use reduction
                  !$acc loop seq
                  do k = 2, km
                     if(k > kb(i, jj) .and. k <= ktcon(i, jj)) then
                        shear= sqrt((uo(i, k, jj)-uo(i, k-1, jj)) ** 2 &
                                  + (vo(i, k, jj)-vo(i, k-1, jj)) ** 2)
                        vshear = vshear + shear
                     endif
                  enddo
                  vshear = 1.e3 * vshear / (zi(i, ktcon(i, jj), jj)-zi(i, kb(i, jj), jj))
                  e1=1.591-.639*vshear &
                     +.0953*(vshear**2)-.00496*(vshear**3)
                  edt(i, jj)=1.-e1
                  val =         .9
                  edt(i, jj) = min(edt(i, jj),val)
                  val =         .0
                  edt(i, jj) = max(edt(i, jj),val)
                  edto(i, jj)=edt(i, jj)
                  edtx(i, jj)=edt(i, jj)
   !
   !  determine detrainment rate between 1 and kbcon
   !
   !> - Next, calculate the variable detrainment rate between the surface and the LFC according to:
   !! \f[
   !! \lambda_d = \frac{1-\beta^{\frac{1}{k_{LFC}}}}{\overline{\Delta z}}
   !! \f]
   !! \f$\lambda_d\f$ is the detrainment rate, \f$\beta\f$ is a constant currently set to 0.05, implying that only 5% of downdraft mass flux at LFC reaches the ground surface due to detrainment, \f$k_{LFC}\f$ is the vertical index of the LFC level, and \f$\o
                  sumx = 0.
      
      ! sumx can use reduction
                  !$acc loop seq
                  do k = 1, km1
                     if(k >= 1 .and. k < kbcon(i, jj)) then
                        dz = zi(i, k+1, jj) - zi(i, k, jj)
                        sumx = sumx + dz
                     endif
                  enddo
              
                  beta = betas
                  if(islimsk(i, jj) == 1) beta = betal
                  dz  = (sumx+zi(i, 1, jj))/float(kbcon(i, jj))
                  tem = 1./float(kbcon(i, jj))
                  xlamd(i, jj) = (1.-beta**tem)/dz
   !
   !  determine downdraft mass flux
   !
   !> - Calculate the normalized downdraft mass flux from equation 1 of Pan and Wu (1995) \cite pan_and_wu_1995 . Downdraft entrainment and detrainment rates are constants from the downdraft origination to the LFC.

                  !$acc loop seq
                  do k = km1, 1, -1
                     if (k <= kmax(i, jj)-1) then
                        if(k < jmin(i, jj) .and. k >= kbcon(i, jj)) then
                           dz        = zi(i, k+1, jj) - zi(i, k, jj)
                           ptem      = xlamdd - xlamde
                           etad(i, k, jj) = etad(i, k+1, jj) * (1. - ptem * dz)
                        else if(k < kbcon(i, jj)) then
                           dz        = zi(i, k+1, jj) - zi(i, k, jj)
                           ptem      = xlamd(i, jj) + xlamdd - xlamde
                           etad(i, k, jj) = etad(i, k+1, jj) * (1. - ptem * dz)
                        endif
                     endif
                  enddo
   !
   !c--- downdraft moisture properties
   !
   !> - Set initial cloud downdraft properties equal to the state variables at the downdraft origination level.
                  jmn = jmin(i, jj)
                  hcdo(i,jmn, jj) = heo(i,jmn, jj)
                  qcdo(i,jmn, jj) = qo(i,jmn, jj)
                  qrcdo(i,jmn, jj)= qo(i,jmn, jj)
                  ucdo(i,jmn, jj) = uo(i,jmn, jj)
                  vcdo(i,jmn, jj) = vo(i,jmn, jj)
                  pwevo = 0.
   !cj
   !> - Calculate the cloud properties as a parcel descends, modified by entrainment and detrainment. Discretization follows Appendix B of Grell (1993) \cite grell_1993 .
                  !$acc loop seq
                  do k = km1, 1, -1
                     if (k < jmin(i, jj)) then
                        dz = zi(i, k+1, jj) - zi(i, k, jj)
                        if(k >= kbcon(i, jj)) then
                           tem  = xlamde * dz
                           tem1 = 0.5 * xlamdd * dz
                        else
                           tem  = xlamde * dz
                           tem1 = 0.5 * (xlamd(i, jj)+xlamdd) * dz
                        endif
                        factor = 1. + tem - tem1
                        hcdo(i, k, jj) = ((1.-tem1)*hcdo(i, k+1, jj)+tem*0.5* &
                                     (heo(i, k, jj)+heo(i, k+1, jj)))/factor
                        dbyo(i, k, jj) = hcdo(i, k, jj) - heso(i, k, jj)
         !         
                        tem  = 0.5 * cm * tem
                        factor = 1. + tem
                        ptem = tem - pgcon
                        ptem1= tem + pgcon
                        ucdo(i, k, jj) = ((1.-tem)*ucdo(i, k+1, jj)+ptem*uo(i, k+1, jj) &
                                     +ptem1*uo(i, k, jj))/factor
                        vcdo(i, k, jj) = ((1.-tem)*vcdo(i, k+1, jj)+ptem*vo(i, k+1, jj) &
                                     +ptem1*vo(i, k, jj))/factor
                     endif
                  enddo
   !
   !> - Compute the amount of moisture that is necessary to keep the downdraft saturated.
                  !$acc loop seq
                  do k = km1, 1, -1
                     if (k < jmin(i, jj)) then
                        gamma      = el2orc * qeso(i, k, jj) / (to(i, k, jj)**2)
                        qrcdo(i, k, jj) = qeso(i, k, jj)+ &
                                (1./hvap)*(gamma/(1.+gamma))*dbyo(i, k, jj)
         !              detad      = etad(i,k+1) - etad(i,k)
         !cj      
                        dz = zi(i, k+1, jj) - zi(i, k, jj)
                        if(k >= kbcon(i, jj)) then
                           tem  = xlamde * dz
                           tem1 = 0.5 * xlamdd * dz
                        else
                           tem  = xlamde * dz
                           tem1 = 0.5 * (xlamd(i, jj)+xlamdd) * dz
                        endif
                        factor = 1. + tem - tem1
                        qcdo(i, k, jj) = ((1.-tem1)*qrcdo(i, k+1, jj)+tem*0.5* &
                                     (qo(i, k, jj)+qo(i, k+1, jj)))/factor
         !cj      
         !              pwdo(i,k)  = etad(i,k+1) * qcdo(i,k+1) -
         !    &                      etad(i,k) * qrcdo(i,k)
         !              pwdo(i,k)  = pwdo(i,k) - detad *
         !    &                     .5 * (qrcdo(i,k) + qrcdo(i,k+1))
         !cj      
                        pwdo(i, k, jj)  = etad(i, k, jj) * (qcdo(i, k, jj) - qrcdo(i, k, jj))
                        pwevo   = pwevo + pwdo(i, k, jj)
                     endif
                  enddo
   !
   !c--- final downdraft strength dependent on precip
   !c--- efficiency (edt), normalized condensate (pwav), and
   !c--- evaporate (pwev)
   !
   !> - Update the precipitation efficiency (edto) based on the ratio of normalized cloud condensate (pwavo) to normalized cloud evaporate (pwevo).
                  edtmax = edtmaxl
                  if(islimsk(i, jj) == 0) edtmax = edtmaxs
                  if(pwevo < 0.) then
                     edto(i, jj) = -edto(i, jj) * pwavo(i, jj) / pwevo
                     edto(i, jj) = min(edto(i, jj),edtmax)
                  else
                     edto(i, jj) = 0.
                  endif
   !
   !c--- downdraft cloudwork functions
   !
   !> - Calculate downdraft cloud work function (\f$A_d\f$) according to equation A.42 (discretized by B.11) in Grell (1993) \cite grell_1993 . Add it to the updraft cloud work function, \f$A_u\f$.
                  !$acc loop seq
                  do k = km1, 1, -1
                     if (k < jmin(i, jj)) then
                        gamma = el2orc * qeso(i, k, jj) / to(i, k, jj)**2
                        dhh=hcdo(i, k, jj)
                        dt=to(i, k, jj)
                        dg=gamma
                        dh=heso(i, k, jj)
                        dz=-1.*(zo(i, k+1, jj)-zo(i, k, jj))
         !              aa1(i)=aa1(i)+edto(i)*dz*etad(i,k)
                        aa1(i, jj)=aa1(i, jj)+edto(i, jj)*dz &
                               *(g/(cp*dt))*((dhh-dh)/(1.+dg)) &
                               *(1.+delta*cp*dg*dt/hvap)
                        val=0.
         !              aa1(i)=aa1(i)+edto(i)*dz*etad(i,k)
                        aa1(i, jj)=aa1(i, jj)+edto(i, jj)*dz &
                               *g*delta*max(val,(qeso(i, k, jj)-qo(i, k, jj)))
                     endif
                  end do
   !> - Check for negative total cloud work function; if found, return to calling routine without modifying state variables.
                  if(aa1(i, jj) <= 0.) then
                     cnvflg(i, jj) = .false.
                  endif
               end if
               wc(i, jj) = wcr
               wbar(i, jj) = wbarr
            end if
         enddo
      end do
   !!
   !!
   !
   !c--- what would the change be, that a cloud with unit mass
   !c--- will do to the environment?
   !
   !> - Calculate the change in moist static energy, moisture mixing ratio, and horizontal winds per unit cloud base mass flux near the surface using equations B.18 and B.19 from Grell (1993) \cite grell_1993, for all layers below cloud top from equations B.
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, jlistnum
         do k = 1, km
            !$acc loop vector
            do i = 1, myim(jj)
               if(cnvflg(i, jj) .and. k <= kmax(i, jj)) then
                  dellah(i, k, jj) = 0.
                  dellaq(i, k, jj) = 0.
                  dellau(i, k, jj) = 0.
                  dellav(i, k, jj) = 0.
               endif
            enddo
         enddo
      end do
   !
   !c--- changed due to subsidence and entrainment
   !
      !$acc parallel loop gang collapse(2) async(async_id) private(aup, adw, &
      !$acc&         dp, dz, dv1h, dv2h, dv3h, dv1q, dv2q, dv3q, tem, tem1, &
      !$acc&         ptem, ptem1, tem2, ptem2, tem3, tem4, tem5, tem6, tem7, tem8, &
      !$acc&         etar, etar1, etadr, etadr1, uor, uor1, vor, vor1)
      do jj = 1, jlistnum
         do k = 2, km1
            !$acc loop vector
            do i = 1, myim(jj)
               if (cnvflg(i, jj) .and. k < ktcon(i, jj)) then
                  etar = eta(i, k, jj)
                  etar1 = eta(i, k-1, jj)
                  etadr = etad(i, k, jj)
                  etadr1 = etad(i, k-1, jj)
                  uor = uo(i, k, jj)
                  uor1 = uo(i, k-1, jj)
                  vor = vo(i, k, jj)
                  vor1 = vo(i, k-1, jj)
                  aup = 1.
                  if(k <= kb(i, jj)) aup = 0.
                  adw = 1.
                  if(k > jmin(i, jj)) adw = 0.
                  dp = 1000. * del(i, k, jj)
                  dz = zi(i, k, jj) - zi(i, k-1, jj)
   !         
                  dv1h = heo(i, k, jj)
                  dv3h = heo(i, k-1, jj)
                  dv2h = .5 * (dv1h + dv3h)
                  dv1q = qo(i, k, jj)
                  dv3q = qo(i, k-1, jj)
                  dv2q = .5 * (dv1q + dv3q)
   !
                  tem  = 0.5 * (xlamue(i, k, jj)+xlamue(i, k-1, jj))
                  tem1 = 0.5 * (xlamud(i, k, jj)+xlamud(i, k-1, jj))
   !          
                  if(k <= kbcon(i, jj)) then
                     ptem  = xlamde
                     ptem1 = xlamd(i, jj)+xlamdd
                  else
                     ptem  = xlamde
                     ptem1 = xlamdd
                  endif
                  tem3 = adw*edto(i, jj)
                  tem4 = aup*etar-tem3*etadr
                  tem5 = aup*etar1-tem3*etadr1
                  tem6 = aup*tem*etar1+tem3*ptem*etadr
                  tem7 = aup*tem1*etar1*.5
                  tem8 = tem3*ptem1*etadr*.5
   !cj
                  dellah(i, k, jj) = (tem4*dv1h - tem5*dv3h - tem6*dv2h*dz &
                   +  tem7*(hcko(i, k, jj)+hcko(i, k-1, jj))*dz &
                   +  tem8*(hcdo(i, k, jj)+hcdo(i, k-1, jj))*dz &
                   ) *g/dp
   !cj     
                  dellaq(i, k, jj) = (tem4*dv1q - tem5*dv3q - tem6*dv2q*dz &
                   +  tem7*(qrcko(i, k, jj)+qcko(i, k-1, jj))*dz &
                   +  tem8*(qrcdo(i, k, jj)+qcdo(i, k-1, jj))*dz &
                   ) *g/dp
   !cj
                 tem1=etar*(uor-ucko(i, k, jj))
                 tem2=etar1*(uor1-ucko(i, k-1, jj))
                 ptem1=etadr*(uor-ucdo(i, k, jj))
                 ptem2=etadr1*(uor1-ucdo(i, k-1, jj))
                 dellau(i, k, jj) = (aup*(tem1-tem2)-tem3*(ptem1-ptem2))*g/dp
   !cj
                 tem1=etar*(vor-vcko(i, k, jj))
                 tem2=etar1*(vor1-vcko(i, k-1, jj))
                 ptem1=etadr*(vor-vcdo(i, k, jj))
                 ptem2=etadr1*(vor1-vcdo(i, k-1, jj))
                 dellav(i, k, jj) = (aup*(tem1-tem2)-tem3*(ptem1-ptem2))*g/dp
   !cj
               endif
            enddo
         enddo
      end do
   !
   !c------- cloud top
   !
      !$acc parallel loop gang vector collapse(2) async(async_id) private(indx, &
      !$acc&         dp, dv1h, dv1q)
      do jj = 1, jlistnum
         do i = 1, ix
            if (i .le. myim(jj)) then
               if(cnvflg(i, jj)) then
                  dp = 1000. * del(i, 1, jj)
                  dellah(i, 1, jj) = edto(i, jj) * etad(i, 1, jj) * (hcdo(i, 1, jj) &
                                 - heo(i, 1, jj)) * g / dp
                  dellaq(i, 1, jj) = edto(i, jj) * etad(i, 1, jj) * (qrcdo(i, 1, jj) &
                                 - qo(i, 1, jj)) * g / dp
                  dellau(i, 1, jj) = edto(i, jj) * etad(i, 1, jj) * (ucdo(i, 1, jj) &
                                 - uo(i, 1, jj)) * g / dp
                  dellav(i, 1, jj) = edto(i, jj) * etad(i, 1, jj) * (vcdo(i, 1, jj) &
                                 - vo(i, 1, jj)) * g / dp
               endif
               if(cnvflg(i, jj)) then
                  indx = ktcon(i, jj)
                  dp = 1000. * del(i, indx, jj)
                  dv1h = heo(i, indx-1, jj)
                  dellah(i, indx, jj) = eta(i, indx-1, jj) * &
                                   (hcko(i, indx-1, jj) - dv1h) * g / dp
                  dv1q = qo(i, indx-1, jj)
                  dellaq(i, indx, jj) = eta(i, indx-1, jj) * &
                                   (qcko(i, indx-1, jj) - dv1q) * g / dp
                  dellau(i, indx, jj) = eta(i, indx-1, jj) * &
                           (ucko(i, indx-1, jj) - uo(i, indx-1, jj)) * g / dp
                  dellav(i, indx, jj) = eta(i, indx-1, jj) * &
                           (vcko(i, indx-1, jj) - vo(i, indx-1, jj)) * g / dp
      !
      !  cloud water
      !
                  dellal(i,indx, jj) = eta(i,indx-1, jj) * &
                                   qlko_ktcon(i, jj) * g / dp
               endif
   
   !
   !c------- final changed variable per unit mass flux
   !
   !> - If grid size is less than a threshold value (dxcrtas: currently 8km), the quasi-equilibrium assumption of Arakawa-Schubert is not used any longer.
   !
               asqecflg(i, jj) = cnvflg(i, jj)
               if(asqecflg(i, jj) .and. gdx(i, jj) < dxcrtas) then
                  asqecflg(i, jj) = .false.
               endif
            end if
         enddo
      end do

   !
   !> - If grid size is larger than the threshold value (i.e., asqecflg=.true.), the quasi-equilibrium assumption is used to obtain the cloud base mass flux. To begin with, calculate the change in the temperature and moisture profiles per unit cloud base mas
      !$acc parallel loop gang collapse(2) async(async_id) private(val, dellat)
      do jj = 1, jlistnum
         do k = 1, km
            !$acc loop vector
            do i = 1, myim(jj)
               if (asqecflg(i, jj) .and. k <= kmax(i, jj)) then
                  if(k > ktcon(i, jj)) then
                     qo(i, k, jj) = q1(i, k, jj)
                     to(i, k, jj) = t1(i, k, jj)
                  endif
                  if(k <= ktcon(i, jj)) then
                     qo(i, k, jj) = dellaq(i, k, jj) * mbdt(i, jj) + q1(i, k, jj)
                     dellat = (dellah(i, k, jj) - hvap * dellaq(i, k, jj)) / cp
                     to(i, k, jj) = dellat * mbdt(i, jj) + t1(i, k, jj)
                     val   =           1.e-10
                     qo(i, k, jj) = max(qo(i, k, jj), val  )
                  endif
               endif
   !c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   !c--- the above changed environment is now used to calulate the
   !c--- effect the arbitrary cloud (with unit mass flux)
   !c--- would have on the stability,
   !c--- which then is used to calculate the real mass flux,
   !c--- necessary to keep this change in balance with the large-scale
   !c--- destabilization.
   !
   !c--- environmental conditions again, first heights
   !
   !> ## Using the updated temperature and moisture profiles that were modified by the convection on a short time-scale, recalculate the total cloud work function to determine the change in the cloud work function due to convection, or the stabilizing effect
   !! - Using notation from Pan and Wu (1995) \cite pan_and_wu_1995, the previously calculated cloud work function is denoted by \f$A^+\f$. Now, it is necessary to use the entraining/detraining cloud model ("static control") to determine the cloud work funct
   !! - Recalculate saturation specific humidity.
               if(asqecflg(i, jj) .and. k <= kmax(i, jj)) then
                  qeso(i, k, jj) = 0.01 * fpvs_gpu(to(i, k, jj),c1xpvs,c2xpvs,tbpvs)      ! fpvs is in pa
   !byl               call qsatq_cwb(1,to(i,k),pfld(i,k),qeso(i,k))
                  qeso(i, k, jj) = eps * qeso(i, k, jj) / (pfld(i, k, jj)+epsm1*qeso(i, k, jj))
                  val       =             1.e-8
                  qeso(i, k, jj) = max(qeso(i, k, jj), val )
   !              tvo(i,k)  = to(i,k) + delta * to(i,k) * qo(i,k)
               endif
            enddo
         enddo
      end do
   !
   !c--- moist static energy
   !
   !! - Recalculate moist static energy and saturation moist static energy.
      !$acc parallel loop gang vector collapse(2) async(async_id) private(dz, k, &
      !$acc&         dp, es, pprime, qs, dqsdp, desdt,dqsdt,gamma, dt, dq, &
      !$acc&         val1, val2, indx, tem, tem1, factor, po)
      do jj = 1, jlistnum
         do i = 1, ix
            if (i .le. myim(jj)) then
               !$acc loop seq
               do k = 1, km1
                  if(asqecflg(i, jj) .and. k <= kmax(i, jj)-1) then
                     dz = .5 * (zo(i, k+1, jj) - zo(i, k, jj))
                     dp = .5 * (pfld(i, k+1, jj) - pfld(i, k, jj))
                     es = 0.01 * fpvs_gpu(to(i, k+1, jj),c1xpvs,c2xpvs,tbpvs)      ! fpvs is in pa
      !byl               call qsatq_cwb(1,to(i,k+1),pfld(i,k+1),es)
                     pprime = pfld(i, k+1, jj) + epsm1 * es
                     qs = eps * es / pprime
                     dqsdp = - qs / pprime
                     desdt = es * (fact1 / to(i, k+1, jj) + fact2 / (to(i, k+1, jj)**2))
                     dqsdt = qs * pfld(i, k+1, jj) * desdt / (es * pprime)
                     gamma = el2orc * qeso(i, k+1, jj) / (to(i, k+1, jj)**2)
                     dt = (g * dz + hvap * dqsdp * dp) / (cp * (1. + gamma))
                     dq = dqsdt * dt + dqsdp * dp
                     to(i, k, jj) = to(i, k+1, jj) + dt
                     qo(i, k, jj) = qo(i, k+1, jj) + dq
                  endif
                  if(asqecflg(i, jj) .and. k <= kmax(i, jj)-1) then
                     po = .5 * (pfld(i, k, jj) + pfld(i, k+1, jj))
                     qeso(i, k, jj) = 0.01 * fpvs_gpu(to(i, k, jj),c1xpvs,c2xpvs,tbpvs)      ! fpvs is in pa
      !byl               call qsatq_cwb(1,to(i,k),po(i,k),qeso(i,k))
                     qeso(i, k, jj) = eps * qeso(i, k, jj) / (po + epsm1 * qeso(i, k, jj))
                     val1      =             1.e-8
                     qeso(i, k, jj) = max(qeso(i, k, jj), val1)
                     val2      =           1.e-10
                     qo(i, k, jj)   = max(qo(i, k, jj), val2 )
      !              qo(i,k)   = min(qo(i,k),qeso(i,k))
                     heo(i, k, jj)   = .5 * g * (zo(i, k, jj) + zo(i, k+1, jj)) + &
                                   cp * to(i, k, jj) + hvap * qo(i, k, jj)
                     heso(i, k, jj) = .5 * g * (zo(i, k, jj) + zo(i, k+1, jj)) + &
                                 cp * to(i, k, jj) + hvap * qeso(i, k, jj)
                  endif
               enddo
               if(asqecflg(i, jj)) then
                  k = kmax(i, jj)
                  heo(i, k, jj) = g * zo(i, k, jj) + cp * to(i, k, jj) + hvap * qo(i, k, jj)
                  heso(i, k, jj) = g * zo(i, k, jj) + cp * to(i, k, jj) + hvap * qeso(i, k, jj)
      !           heo(i,k) = min(heo(i,k),heso(i,k))
               endif
   !
   !c**************************** static control
   !
   !c------- moisture and cloud work functions
   !
   !> - As before, recalculate the updraft cloud work function.
               if(asqecflg(i, jj)) then
                  xaa0(i, jj) = 0.
                  xpwav(i, jj) = 0.
               endif
   !
               if(asqecflg(i, jj)) then
                  indx = kb(i, jj)
                  hcko(i, indx, jj) = heo(i, indx, jj)
                  qcko(i, indx, jj) = qo(i, indx, jj)
               endif

               !$acc loop seq
               do k = 2, km1
                  if (asqecflg(i, jj)) then
                     if(k > kb(i, jj) .and. k <= ktcon(i, jj)) then
                        dz = zi(i, k, jj) - zi(i, k-1, jj)
                        tem  = 0.5 * (xlamue(i, k, jj)+xlamue(i, k-1, jj)) * dz
                        tem1 = 0.25 * (xlamud(i, k, jj)+xlamud(i, k-1, jj)) * dz
                        factor = 1. + tem - tem1
                        hcko(i, k, jj) = ((1.-tem1)*hcko(i, k-1, jj)+tem*0.5* &
                                     (heo(i, k, jj)+heo(i, k-1, jj)))/factor
                     endif
                  endif
               enddo
            end if
         enddo
      end do
      
      !$acc parallel loop gang vector collapse(2) async(async_id) private( &
      !$acc&         dz, gamma, xdby, xqrch, tem, tem1, factor, dq, etah, ptem, &
      !$acc&         qlk, xpw, dz1, rfact, val, edtmax, jmn, dt, dh, xpwd)
      do jj = 1, jlistnum
         do i = 1, ix
            if (i .le. myim(jj)) then
               !$acc loop seq
               do k = 2, km1
                  if (asqecflg(i, jj)) then
                     if(k > kb(i, jj) .and. k < ktcon(i, jj)) then
                        dz = zi(i, k, jj) - zi(i, k-1, jj)
                        gamma = el2orc * qeso(i, k, jj) / (to(i, k, jj)**2)
                        xdby = hcko(i, k, jj) - heso(i, k, jj)
                        xqrch = qeso(i, k, jj) &
                              + gamma * xdby / (hvap * (1. + gamma))
      !cj         
                        tem  = 0.5 * (xlamue(i, k, jj)+xlamue(i, k-1, jj)) * dz
                        tem1 = 0.25 * (xlamud(i, k, jj)+xlamud(i, k-1, jj)) * dz
                        factor = 1. + tem - tem1
                        qcko(i, k, jj) = ((1.-tem1)*qcko(i, k-1, jj)+tem*0.5* &
                                     (qo(i, k, jj)+qo(i, k-1, jj)))/factor
      !cj         
                        dq = eta(i, k, jj) * (qcko(i, k, jj) - xqrch)
      !           
                        if(k >= kbcon(i, jj) .and. dq > 0.) then
                           etah = .5 * (eta(i, k, jj) + eta(i, k-1, jj))
                           if(ncloud > 0 .and. k > jmin(i, jj)) then
                              ptem = c0t(i, k, jj) + c1
                              qlk = dq / (eta(i, k, jj) + etah * ptem * dz)
                           else
                              qlk = dq / (eta(i, k, jj) + etah * c0t(i, k, jj) * dz)
                           endif
                           if(k < ktcon1(i, jj)) then
         !                    xaa0(i) = xaa0(i) - dz * g * qlk * etah
                              xaa0(i, jj) = xaa0(i, jj) - dz * g * qlk
                           endif
                           qcko(i, k, jj) = qlk + xqrch
                           xpw = etah * c0t(i, k, jj) * dz * qlk
                           xpwav(i, jj) = xpwav(i, jj) + xpw
                        endif
                     endif
                     if(k >= kbcon(i, jj) .and. k < ktcon1(i, jj)) then
                        dz1 = zo(i, k+1, jj) - zo(i, k, jj)
                        gamma = el2orc * qeso(i, k, jj) / (to(i, k, jj)**2)
                        rfact =  1. + delta * cp * gamma &
                                 * to(i, k, jj) / hvap
                        xaa0(i, jj) = xaa0(i, jj)                       &
      !    &                    + dz1 * eta(i,k) * (g / (cp * to(i,k))) &
                                + dz1 * (g / (cp * to(i, k, jj))) &
                                * xdby / (1. + gamma) &
                                * rfact
                        val=0.
                        xaa0(i, jj) = xaa0(i, jj) +                     &
      !    &                     dz1 * eta(i,k) * g * delta * &
                                 dz1 * g * delta * &
                                 max(val,(qeso(i, k, jj) - qo(i, k, jj)))
                     endif
                  endif
               enddo
   !
   !c------- downdraft calculations
   !
   !c--- downdraft moisture properties
   !
   !> - As before, recalculate the downdraft cloud work function.
               if(asqecflg(i, jj)) then
                  jmn = jmin(i, jj)
                  hcdo(i, jmn, jj) = heo(i, jmn, jj)
                  qcdo(i, jmn, jj) = qo(i, jmn, jj)
                  qrcd(i, jmn, jj) = qo(i, jmn, jj)
                  xpwev(i, jj) = 0.
               endif
   !cj
               !$acc loop seq
               do k = km1, 1, -1
                  if (asqecflg(i, jj) .and. k < jmin(i, jj)) then
                       dz = zi(i, k+1, jj) - zi(i, k, jj)
                       if(k >= kbcon(i, jj)) then
                           tem  = xlamde * dz
                           tem1 = 0.5 * xlamdd * dz
                       else
                           tem  = xlamde * dz
                           tem1 = 0.5 * (xlamd(i, jj)+xlamdd) * dz
                       endif
                       factor = 1. + tem - tem1
                       hcdo(i, k, jj) = ((1.-tem1)*hcdo(i, k+1, jj)+tem*0.5* &
                                    (heo(i, k, jj)+heo(i, k+1, jj)))/factor
                  endif
               enddo
   !cj
               !$acc loop seq
               do k = km1, 1, -1
                  if (asqecflg(i, jj) .and. k < jmin(i, jj)) then
                     dq = qeso(i, k, jj)
                     dt = to(i, k, jj)
                     gamma    = el2orc * dq / dt**2
                     dh       = hcdo(i, k, jj) - heso(i, k, jj)
                     qrcd(i, k, jj)=dq+(1./hvap)*(gamma/(1.+gamma))*dh
      !              detad    = etad(i,k+1) - etad(i,k)
      !cj      
                     dz = zi(i, k+1, jj) - zi(i, k, jj)
                     if(k >= kbcon(i, jj)) then
                        tem  = xlamde * dz
                        tem1 = 0.5 * xlamdd * dz
                     else
                        tem  = xlamde * dz
                        tem1 = 0.5 * (xlamd(i, jj)+xlamdd) * dz
                     endif
                     factor = 1. + tem - tem1
                     qcdo(i, k, jj) = ((1.-tem1)*qrcd(i, k+1, jj)+tem*0.5* &
                                  (qo(i, k, jj)+qo(i, k+1, jj)))/factor
      !cj      
      !              xpwd     = etad(i,k+1) * qcdo(i,k+1) -
      !    &                    etad(i,k) * qrcd(i,k)
      !              xpwd     = xpwd - detad *
      !    &                  .5 * (qrcd(i,k) + qrcd(i,k+1))
      !cj      
                     xpwd     = etad(i, k, jj) * (qcdo(i, k, jj) - qrcd(i, k, jj))
                     xpwev(i, jj) = xpwev(i, jj) + xpwd
                  endif
               enddo
   !
               edtmax = edtmaxl
               if(islimsk(i, jj) == 0) edtmax = edtmaxs
               if(asqecflg(i, jj)) then
                  if(xpwev(i, jj) >= 0.) then
                     edtx(i, jj) = 0.
                  else
                     edtx(i, jj) = -edtx(i, jj) * xpwav(i, jj) / xpwev(i, jj)
                     edtx(i, jj) = min(edtx(i, jj),edtmax)
                  endif
               endif
            end if
         enddo
      end do
   !
   !
   !c--- downdraft cloudwork functions
   !
   !
      ! xaa0 can use reduction
      !$acc parallel loop gang vector collapse(2) async(async_id) private(gamma, dhh, &
      !$acc&         dt, dg, dh, dz, val, tem, tfac, k, rho, po, sumx, umean)
      do jj = 1, jlistnum
         do i = 1, ix
            if (i .le. myim(jj)) then
               !$acc loop seq
               do k = km1, 1, -1
                  if (asqecflg(i, jj) .and. k < jmin(i, jj)) then
                     gamma = el2orc * qeso(i, k, jj) / to(i, k, jj)**2
                     dhh=hcdo(i, k, jj)
                     dt= to(i, k, jj)
                     dg= gamma
                     dh= heso(i, k, jj)
                     dz=-1.*(zo(i, k+1, jj)-zo(i, k, jj))
      !              xaa0(i)=xaa0(i)+edtx(i)*dz*etad(i,k)
                     xaa0(i, jj)=xaa0(i, jj)+edtx(i, jj)*dz &
                             *(g/(cp*dt))*((dhh-dh)/(1.+dg)) &
                             *(1.+delta*cp*dg*dt/hvap)
                     val=0.
      !              xaa0(i)=xaa0(i)+edtx(i)*dz*etad(i,k)
                     xaa0(i, jj)=xaa0(i, jj)+edtx(i, jj)*dz &
                             *g*delta*max(val,(qeso(i, k, jj)-qo(i, k, jj)))
                  endif
               enddo
   !
   !  calculate critical cloud work function
   !
   !     do i = 1, im
   !       if(cnvflg(i)) then
   !         if(pfld(i,ktcon(i)) < pcrit(15))then
   !           acrt(i)=acrit(15)*(975.-pfld(i,ktcon(i)))
   !    &              /(975.-pcrit(15))
   !         else if(pfld(i,ktcon(i)) > pcrit(1))then
   !           acrt(i)=acrit(1)
   !         else
   !           k =  int((850. - pfld(i,ktcon(i)))/50.) + 2
   !           k = min(k,15)
   !           k = max(k,2)
   !           acrt(i)=acrit(k)+(acrit(k-1)-acrit(k))*
   !    &           (pfld(i,ktcon(i))-pcrit(k))/(pcrit(k-1)-pcrit(k))
   !         endif
   !       endif
   !     enddo
   !     do i = 1, im
   !       if(cnvflg(i)) then
   !         if(islimsk(i) == 1) then
   !           w1 = w1l
   !           w2 = w2l
   !           w3 = w3l
   !           w4 = w4l
   !         else
   !           w1 = w1s
   !           w2 = w2s
   !           w3 = w3s
   !           w4 = w4s
   !         endif
   !
   !  modify critical cloud workfunction by cloud base vertical velocity
   !
   !         if(pdot(i) <= w4) then
   !           acrtfct(i) = (pdot(i) - w4) / (w3 - w4)
   !         elseif(pdot(i) >= -w4) then
   !           acrtfct(i) = - (pdot(i) + w4) / (w4 - w3)
   !         else
   !           acrtfct(i) = 0.
   !         endif
   !         val1    =            -1.
   !         acrtfct(i) = max(acrtfct(i),val1)
   !         val2    =             1.
   !         acrtfct(i) = min(acrtfct(i),val2)
   !         acrtfct(i) = 1. - acrtfct(i)
   !
   !  modify acrtfct(i) by colume mean rh if rhbar(i) is greater than 80 percent
   !
   !         if(rhbar(i) >= .8) then
   !           acrtfct(i) = acrtfct(i) * (.9 - min(rhbar(i),.9)) * 10.
   !         endif
   !
   !  modify adjustment time scale by cloud base vertical velocity
   !
   !         dtconv(i) = dt2 + max((1800. - dt2),0.) *
   !    &                (pdot(i) - w2) / (w1 - w2)
   !         dtconv(i) = max(dtconv(i), dt2)
   !         dtconv(i) = 1800. * (pdot(i) - w2) / (w1 - w2)
   !
   !         dtconv(i) = max(dtconv(i),dtmin)
   !         dtconv(i) = min(dtconv(i),dtmax)
   !
   !       endif
   !     enddo
   !
   !  compute convective turn-over time
   !
   !> - Following Bechtold et al. (2008) \cite bechtold_et_al_2008, the convective adjustment time (dtconv) is set to be proportional to the convective turnover time, which is computed using the mean updraft velocity (wc) and the cloud depth. It is also prop
               if(cnvflg(i, jj)) then
                  tem = zi(i, ktcon1(i, jj), jj) - zi(i, kbcon1(i, jj), jj)
                  dtconv(i, jj) = tem / wc(i, jj)
                  tfac = 1. + gdx(i, jj) / 75000.
      ! reference from eq.3 in Zheng et al. 2016
      !byl          tfac = 1. + log(25000./gdx(i))
                  dtconv(i, jj) = tfac * dtconv(i, jj)
                  dtconv(i, jj) = max(dtconv(i, jj), dtmin)
                  dtconv(i, jj) = min(dtconv(i, jj), dtmax)
   !
   !> - Calculate advective time scale (tauadv) using a mean cloud layer wind speed.
                  sumx = 0.
                  umean = 0.

                  !$acc loop seq
                  do k = 2, km1
                     if(k >= kbcon1(i, jj) .and. k < ktcon1(i, jj)) then
                        dz = zi(i, k, jj) - zi(i, k-1, jj)
                        tem = sqrt(u1(i, k, jj)*u1(i, k, jj)+v1(i, k, jj)*v1(i, k, jj))
                        umean = umean + tem * dz
                        sumx = sumx + dz
                     endif
                  enddo

                  umean = umean / sumx
                  umean = max(umean, 1.)
                  tauadv(i, jj) = gdx(i, jj) / umean
   !> - From Han et al.'s (2017) \cite han_et_al_2017 equation 6, calculate cloud base mass flux as a function of the mean updraft velcoity for the grid sizes where the quasi-equilibrium assumption of Arakawa-Schubert is not valid any longer.
   !!  As discussed in Han et al. (2017) \cite han_et_al_2017 , when dtconv is larger than tauadv, the convective mixing is not fully conducted before the cumulus cloud is advected out of the grid cell. In this case, therefore, the cloud base mass flux is fu
      
                  if(.not.asqecflg(i, jj)) then
                     k = kbcon(i, jj)
                     po = .5 * (pfld(i, k, jj) + pfld(i, k+1, jj))
                     rho = po*100. / (rd*to(i, k, jj))
                     tfac = tauadv(i, jj) / dtconv(i, jj)
                     tfac = min(tfac, 1.)
                     xmb(i, jj) = tfac*betaw*rho*wc(i, jj)
                  end if
               endif
   !> - For the cases where the quasi-equilibrium assumption of Arakawa-Schubert is valid, first calculate the large scale destabilization as in equation 5 of Pan and Wu (1995) \cite pan_and_wu_1995 :
   !! \f[
   !!  \frac{\partial A}{\partial t}_{LS}=\frac{A^+-cA^0}{\Delta t_{LS}}
   !! \f]
   !! Here \f$A^0\f$ is set to zero following  Han et al.'s (2017) \cite han_et_al_2017 , implying that the instability is completely eliminated after the convective adjustment time, \f$\Delta t_{LS}\f$.
      
               if(asqecflg(i, jj)) then
      !           fld(i)=(aa1(i)-acrt(i)*acrtfct(i))/dtconv(i)
                  fld(i, jj)=aa1(i, jj)/dtconv(i, jj)
                  if(fld(i, jj) <= 0.) then
                     asqecflg(i, jj) = .false.
                     cnvflg(i, jj) = .false.
                  endif
               endif
      !> - Calculate the stabilization effect of the convection (per unit cloud base mass flux) as in equation 6 of Pan and Wu (1995) \cite pan_and_wu_1995 :
      !! \f[
      !! \frac{\partial A}{\partial t}_{cu}=\frac{A^*-A^+}{\Delta t_{cu}}
      !! \f]
      !! \f$\Delta t_{cu}\f$ is the short timescale of the convection.
               if(asqecflg(i, jj)) then
      !           xaa0(i) = max(xaa0(i),0.)
                  xk(i, jj) = (xaa0(i, jj) - aa1(i, jj)) / mbdt(i, jj)
                  if(xk(i, jj) >= 0.) then
                     asqecflg(i, jj) = .false.
                     cnvflg(i, jj) = .false.
                  endif
               endif
      !
      !c--- kernel, cloud base mass flux
      !
      !> - The cloud base mass flux (xmb) is then calculated from equation 7 of Pan and Wu (1995) \cite pan_and_wu_1995
      !! \f[
      !! M_c=\frac{-\frac{\partial A}{\partial t}_{LS}}{\frac{\partial A}{\partial t}_{cu}}
      !! \f]
      !!
      !!  Again when dtconv is larger than tauadv, the cloud base mass flux is further reduced in proportion to the ratio of tauadv to dtconv.
               if(asqecflg(i, jj)) then
                  tfac = tauadv(i, jj) / dtconv(i, jj)
                  tfac = min(tfac, 1.)
                  xmb(i, jj) = -tfac * fld(i, jj) / xk(i, jj)
      !           xmb(i) = min(xmb(i),xmbmax(i))
               endif
   !!
   !> - If the large scale destabilization is less than zero, or the stabilization by the convection is greater than zero, then the scheme returns to the calling routine without modifying the state variables.
      !   totflg = .true.
      !do jj = 1, jlistnum
      !   do i=1,myim(jj)
      !     totflg = totflg .and. (.not. cnvflg(i, jj))
      !   enddo
      !end do
         !if(totflg) cycle
   !!
   !
   !> - For scale-aware parameterization, the updraft fraction (sigmagfm) is first computed as a function of the lateral entrainment rate at cloud base (see Han et al.'s (2017) \cite han_et_al_2017 equation 4 and 5), following the study by Grell and Freitas
               if(cnvflg(i, jj)) then
                  sigmagfm(i, jj) = wbar(i, jj) / wc(i, jj) !scale-aware parameterization (Kwon and Hung 2017)
                  sigmagfm(i, jj) = max(sigmagfm(i, jj), 0.001)
                  sigmagfm(i, jj) = min(sigmagfm(i, jj), 0.999)
   !
   !> - Then, calculate the reduction factor (scaldfunc) of the vertical convective eddy transport of mass flux as a function of updraft fraction from the studies by Arakawa and Wu (2013) \cite arakawa_and_wu_2013 (also see Han et al.'s (2017) \cite han_et_a
      
                  scaldfunc(i, jj) = clear(i, jj) * (1.-sigma(i, jj))
                  scaldfunc(i, jj) = max(min(scaldfunc(i, jj), 1.0), 0.)
                  xmb(i, jj) = xmb(i, jj) * scaldfunc(i, jj)
                  xmb(i, jj) = min(xmb(i, jj),xmbmax(i, jj))
   !
   !  restore to,qo,uo,vo to t1,q1,u1,v1 in case convection stops
   !
                  !$acc loop seq
                  do k = 1, km
                     if (k <= kmax(i, jj)) then
                        to(i, k, jj) = t1(i, k, jj)
                        qo(i, k, jj) = q1(i, k, jj)
                        uo(i, k, jj) = u1(i, k, jj)
                        vo(i, k, jj) = v1(i, k, jj)
                        qeso(i, k, jj) = 0.01 * fpvs_gpu(t1(i, k, jj),c1xpvs,c2xpvs,tbpvs)      ! fpvs is in pa
         !byl               call qsatq_cwb(1,t1(i,k),pfld(i,k),qeso(i,k))
                        qeso(i, k, jj) = eps * qeso(i, k, jj) / (pfld(i, k, jj) + epsm1*qeso(i, k, jj))
                        val     =             1.e-8
                        qeso(i, k, jj) = max(qeso(i, k, jj), val )
                     endif
                  enddo
               end if
            end if
         enddo
      end do
   !c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   !c--- feedback: simply the changes from the cloud with unit mass flux
   !c---           multiplied by  the mass flux necessary to keep the
   !c---           equilibrium with the larger-scale.
   !
   !> ## For the "feedback" control, calculate updated values of the state variables by multiplying the cloud base mass flux and the tendencies calculated per unit cloud base mass flux from the static control.
   !> - Calculate the temperature tendency from the moist static energy and specific humidity tendencies.
   !> - Update the temperature, specific humidity, and horiztonal wind state variables by multiplying the cloud base mass flux-normalized tendencies by the cloud base mass flux.
   !> - Accumulate column-integrated tendencies.
      !$acc parallel loop gang vector collapse(2) async(async_id) private(dp, dellat)
      do jj = 1, jlistnum
         do i = 1, ix
            if (i .le. myim(jj)) then
               delhbar(i, jj) = 0.
               delqbar(i, jj) = 0.
               deltbar(i, jj) = 0.
               delubar(i, jj) = 0.
               delvbar(i, jj) = 0.
   !
   !> - Add up column-integrated convective precipitation by multiplying the normalized value by the cloud base mass flux.
               rntot(i, jj) = 0.
               delqev(i, jj) = 0.
               flg(i, jj) = cnvflg(i, jj)
            end if
         enddo
      end do
      
      !$acc parallel loop gang collapse(2) async(async_id) private(dp, dellat, val)
      do jj = 1, jlistnum
         do k = 1, km
            !$acc loop vector
            do i = 1, myim(jj)
               if (cnvflg(i, jj) .and. k <= kmax(i, jj)) then
                  if(k <= ktcon(i, jj)) then
                     dellat = (dellah(i, k, jj) - hvap * dellaq(i, k, jj)) / cp
                     t1(i, k, jj) = t1(i, k, jj) + dellat * xmb(i, jj) * dt2
                     q1(i, k, jj) = q1(i, k, jj) + dellaq(i, k, jj) * xmb(i, jj) * dt2
   !                 tem = 1./rcs(i)
   !                 u1(i,k) = u1(i,k) + dellau(i,k) * xmb(i) * dt2 * tem
   !                 v1(i,k) = v1(i,k) + dellav(i,k) * xmb(i) * dt2 * tem
                     u1(i, k, jj) = u1(i, k, jj) + dellau(i, k, jj) * xmb(i, jj) * dt2
                     v1(i, k, jj) = v1(i, k, jj) + dellav(i, k, jj) * xmb(i, jj) * dt2
                     dp = 1000. * del(i, k, jj)
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
   !> - Recalculate saturation specific humidity using the updated temperature.
                    qeso(i, k, jj) = 0.01 * fpvs_gpu(t1(i, k, jj),c1xpvs,c2xpvs,tbpvs)      ! fpvs is in pa
   !byl                 call qsatq_cwb(1,t1(i,k),pfld(i,k),qeso(i,k))
                    qeso(i, k, jj) = eps * qeso(i, k, jj)/(pfld(i, k, jj) + epsm1*qeso(i, k, jj))
                    val     =             1.e-8
                    qeso(i, k, jj) = max(qeso(i, k, jj), val )
                 endif
              endif
            enddo
         enddo
      end do


   !> - Determine the evaporation of the convective precipitation and update the integrated convective precipitation.
   !> - Update state temperature and moisture to account for evaporation of convective precipitation.
   !> - Update column-integrated tendencies to account for evaporation of convective precipitation.
      !$acc parallel loop gang vector collapse(2) async(async_id) private( &
      !$acc&         deltv, delq, qevap, delq2, aup, adw, rain, evef, qcond, &
      !$acc&         dp, delq2, rainf, tem1)
      do jj = 1, jlistnum
         do i = 1, ix
            if (i .le. myim(jj)) then
               !$acc loop seq
               do k = km1, 1, -1
                  if (cnvflg(i, jj) .and. k <= kmax(i, jj)) then
                     if(k < ktcon(i, jj)) then
                        aup = 1.
                        if(k <= kb(i, jj)) aup = 0.
                        adw = 1.
                        if(k >= jmin(i, jj)) adw = 0.
                        rain =  aup * pwo(i, k, jj) + adw * edto(i, jj) * pwdo(i, k, jj)
                        rainf = rain * xmb(i, jj)    !xb110
                        rntot(i, jj) = rntot(i, jj) + rain * xmb(i, jj) * .001 * dt2
      !xb110>  for lightning parameterization
                        tem1 = max(0.0, min(1.0, (tcr-t1(i, k, jj))*tcrf))
                        snow_flx(i, k, jj) = rainf*tem1 !snow flux
      !xb110<
                     endif
                  endif
               end do
               !$acc loop seq
               do k = km1, 1, -1
                  if (k <= kmax(i, jj)) then
                     deltv = 0.
                     delq = 0.
                     qevap = 0.
                     delq2 = 0.
                     if(cnvflg(i, jj) .and. k < ktcon(i, jj)) then
                        aup = 1.
                        if(k <= kb(i, jj)) aup = 0.
                        adw = 1.
                        if(k >= jmin(i, jj)) adw = 0.
                        rain =  aup * pwo(i, k, jj) + adw * edto(i, jj) * pwdo(i, k, jj)
                        rn(i, jj) = rn(i, jj) + rain * xmb(i, jj) * .001 * dt2
                     endif
                     if(flg(i, jj) .and. k < ktcon(i, jj)) then
                        evef = edt(i, jj) * evfact
                        if(islimsk(i, jj) == 1) evef=edt(i, jj) * evfactl
      !                 if(islimsk(i) == 1) evef=.07
      !                 if(islimsk(i) == 1) evef = 0.
                        qcond = evef * (q1(i, k, jj) - qeso(i, k, jj)) &
                                 / (1. + el2orc * qeso(i, k, jj) / t1(i, k, jj)**2)
                        dp = 1000. * del(i, k, jj)
                        if(rn(i, jj) > 0. .and. qcond < 0.) then
                           qevap = -qcond * (1.-exp(-.32*sqrt(dt2*rn(i, jj))))
                           qevap = min(qevap, rn(i, jj)*1000.*g/dp)
                           delq2 = delqev(i, jj) + .001 * qevap * dp / g
                        endif
                        if(rn(i, jj) > 0. .and. qcond < 0. .and. &
                           delq2 > rntot(i, jj)) then
                           qevap = 1000.* g * (rntot(i, jj) - delqev(i, jj)) / dp
                           flg(i, jj) = .false.
                        endif
                        if(rn(i, jj) > 0. .and. qevap > 0.) then
                           q1(i, k, jj) = q1(i, k, jj) + qevap
                           t1(i, k, jj) = t1(i, k, jj) - elocp * qevap
                           rn(i, jj) = rn(i, jj) - .001 * qevap * dp / g
                           deltv = - elocp*qevap/dt2
                           delq =  + qevap/dt2
                           delqev(i, jj) = delqev(i, jj) + .001*dp*qevap/g
                        endif
                        delqbar(i, jj) = delqbar(i, jj) + delq*dp/g
                        deltbar(i, jj) = deltbar(i, jj) + deltv*dp/g
                     endif
                  endif
              enddo
   !cj
   !     do i = 1, im
   !     if(me == 31 .and. cnvflg(i)) then
   !     if(cnvflg(i)) then
   !       print *, ' deep delhbar, delqbar, deltbar = ',
   !    &             delhbar(i),hvap*delqbar(i),cp*deltbar(i)
   !       print *, ' deep delubar, delvbar = ',delubar(i),delvbar(i)
   !       print *, ' precip =', hvap*rn(i)*1000./dt2
   !       print*,'pdif= ',pfld(i,kbcon(i))-pfld(i,ktcon(i))
   !     endif
   !     enddo
   !
   !  precipitation rate converted to actual precip
   !  in unit of m instead of kg
   !
               if(cnvflg(i, jj)) then
      !
      !  in the event of upper level rain evaporation and lower level downdraft
      !    moistening, rn can become negative, in this case, we back out of the
      !    heating and the moistening
      !
                  if(rn(i, jj) < 0. .and. .not.flg(i, jj)) rn(i, jj) = 0.
                  if(rn(i, jj) < 1.e-13) rn(i, jj) = 0.
                  if(rn(i, jj) <= 0.) then
                     rn(i, jj) = 0.
                  else
                     ktop(i, jj) = ktcon(i, jj)
                     kbot(i, jj) = kbcon(i, jj)
                     kcnv(i, jj) = 1
                     cldwrk(i, jj) = aa1(i, jj)
                  endif
               endif
            end if
         enddo
      end do
      !$acc parallel loop gang collapse(2) async(async_id) private(tem, tem1, &
      !$acc&         dz, rho, ptem, po1)
      do jj = 1, jlistnum
         do k = 1, km
            !$acc loop vector
            do i = 1, myim(jj)
   !
   !  convective cloud water
   !
   !> - Calculate convective cloud water.
               if (cnvflg(i, jj) .and. rn(i, jj) > 0.) then
                  if (k >= kbcon(i, jj) .and. k < ktcon(i, jj)) then
                     cnvw(i, k, jj) = cnvwt(i, k, jj) * xmb(i, jj) * dt2
                  endif
               endif
   !
   !  convective cloud cover
   !
   !> - Calculate convective cloud cover, which is used when pdf-based cloud fraction is used (i.e., pdfcld=.true.).
               if (cnvflg(i, jj) .and. rn(i, jj) > 0.) then
                  if (k >= kbcon(i, jj) .and. k < ktcon(i, jj)) then
                     cnvc(i, k, jj) = 0.04 * log(1. + 675. * eta(i, k, jj) * xmb(i, jj))
                     cnvc(i, k, jj) = min(cnvc(i, k, jj), 0.6)
                     cnvc(i, k, jj) = max(cnvc(i, k, jj), 0.0)
                  endif
               endif
   !
   !  cloud water
   !
   !> - Separate detrained cloud water into liquid and ice species as a function of temperature only.
               if (ncloud > 0) then
                  if (cnvflg(i, jj) .and. rn(i, jj) > 0.) then
      !              if (k > kb(i) .and. k <= ktcon(i)) then
                     if (k >= kbcon(i, jj) .and. k <= ktcon(i, jj)) then
                        tem  = dellal(i, k, jj) * xmb(i, jj) * dt2
                        tem1 = max(0.0, min(1.0, (tcr-t1(i, k, jj))*tcrf))
                        if (qi(i, k, jj) > -999.0) then
                           qi(i, k, jj) = qi(i, k, jj) + tem * tem1            ! ice
                           ql(i, k, jj) = ql(i, k, jj) + tem *(1.0-tem1)       ! water
                        else
                           ql(i, k, jj) = ql(i, k, jj) + tem
                        endif
                     endif
                  endif
               endif
   !
   !> - If convective precipitation is zero or negative, reset the updated state variables back to their original values (negating convective changes).
               if(cnvflg(i, jj) .and. rn(i, jj) <= 0.) then
                  if (k <= kmax(i, jj)) then
                     t1(i, k, jj) = to(i, k, jj)
                     q1(i, k, jj) = qo(i, k, jj)
                     u1(i, k, jj) = uo(i, k, jj)
                     v1(i, k, jj) = vo(i, k, jj)
                  endif
               endif
   !
   ! hchuang code change
   !
   !> - Calculate and retain the updraft and downdraft mass fluxes for dust transport by cumulus convection.
   !
   !> - Calculate the updraft convective mass flux.
               if(cnvflg(i, jj) .and. rn(i, jj) > 0.) then
                  if(k >= kb(i, jj) .and. k < ktop(i, jj)) then
                     ud_mf(i, k, jj) = eta(i, k, jj) * xmb(i, jj) * dt2
                  endif
               endif
   !
   !
   !> - Calculate the downdraft convective mass flux.
               if(cnvflg(i, jj) .and. rn(i, jj) > 0.) then
                  if(k >= 1 .and. k <= jmin(i, jj)) then
                     dd_mf(i, k, jj) = edto(i, jj) * etad(i, k, jj) * xmb(i, jj) * dt2
                  endif
               endif
             
               if(cnvflg(i, jj) .and. rn(i, jj) > 0.) then
                  if (k .ge. 2 .and. k .le. km1) then
                     if(k < ktop(i, jj)) then
      !(test) dissipation 
                        dz = zi(i, k, jj) - zi(i, k-1, jj)
                        po1 = .5 * (pfld(i, k, jj) + pfld(i, k+1, jj))
                        rho = po1*100. / (rd*to(i, k, jj))
      !                  tem=(ud_mf(i,k)+dd_mf(i,k))/dt2    !(kg/m^2)
                        tem=ud_mf(i, k, jj)/dt2    !(kg/m^2)
                        ptem= 0.5*(xlamud(i, k, jj)+xlamud(i, k+1, jj))   !(1/m)
      !                  ptem2= 0.5*(xlamue(i,k)+xlamue(i,k+1))   !(1/m)
                        dissdc1(i, k, jj)=(tem**2.)*(ptem)/(rho**3.)*1.e4/(2.6*2.6)*15.  !from Berner_2009*15.
                     endif
                  end if
               endif
   !!
   !xb110>> storage the updarft T and q for flash parameterization
               if (cnvflg(i, jj))then
                  if(k >= kbcon(i, jj) .and. k < ktcon(i, jj)) then
                     ptu(i, k, jj) = (hcko(i, k, jj)-hvap*qcko(i, k, jj)-phil(i, k, jj))/cp
                     pqu(i, k, jj) = qcko(i, k, jj)
                  end if
               end if
            enddo
         enddo
      end do

   !> - save the updraft convective mass flux at cloud top.
      !$acc parallel loop gang vector collapse(2) async(async_id) private(k)
      do jj = 1, jlistnum
         do i = 1, ix
            if (i .le. myim(jj)) then
               if(cnvflg(i, jj) .and. rn(i, jj) > 0.) then
                  k = ktop(i, jj)-1
                  dt_mf(i, k, jj) = ud_mf(i, k, jj)
               endif
            end if
         enddo
      end do
      !$acc end data
      !$acc wait async(async_id)
!      if (myrank .eq.0) print*,"max ptu:",maxval(ptu)
!      if (myrank .eq.0) print*,"max pqu:",maxval(pqu)
!b110<<
      return
      end
