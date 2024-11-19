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
      subroutine moninedmf_gpu(ix,im,km,ntrac,ntcw,                         &
         u1,v1,t1,q1,swh,hlw,xmu,                                       &
         psk,rbsoil,zorl,u10m,v10m,fm,fh,                               &
         tsea,heat,evap,stress,spd1,kpbl,                               &
         prsi,del,prsl,prslk,phii,phil,delt,hpbl)
!
      !$acc routine(tridi2_gpu) vector
      !$acc routine(tridin_gpu) vector
      use machine  , only : kind_phys
      use physcons, grav => con_g, rd => con_rd, cp => con_cp &
      ,             hvap => con_hvap, fv => con_fvirt
      use const    , only : RTYPE
      use param, only : my, my_max
      use index, only: jlist1, jlistnum
      use rank, only: myrank
      use openacc
      use cudafor

      implicit none
!
!     arguments
!
      logical lprnt
      integer ipr
      integer ix, im(my), myim(my_max), km, ntrac, ntcw, kpbl(ix,my_max), kinver(ix)
!
      real(kind=kind_phys) delt, xkzm_m, xkzm_h, xkzm_s
      real(kind=kind_phys) tau(ix,km,my_max),    rtg(ix,km,ntrac,my_max),             &
                           u1(ix,km,my_max),     v1(ix,km,my_max),                    &
                           t1(ix,km,my_max),     q1(ix,km,ntrac,my_max),              &
                           swh(ix,km,my_max),    hlw(ix,km,my_max),                   &
                           xmu(ix,my_max),                                     &
                           rbsoil(ix,my_max),    zorl(ix,my_max),                     &
                           u10m(ix,my_max),      v10m(ix,my_max),                     &
                           fm(ix,my_max),        fh(ix,my_max),                       &
                           tsea(ix,my_max),                                    &
                                          spd1(ix,my_max),                     &
                           prsi(ix,km+1,my_max), del(ix,km,my_max),                   &
                           prsl(ix,km,my_max),   prslk(ix,km,my_max),                 &
                           phii(ix,km+1), phil(ix,km,my_max),                  &
                           dusfc(ix,my_max),     dvsfc(ix,my_max),                    &
                           dtsfc(ix,my_max),     dqsfc(ix,my_max),                    &
                           hpbl(ix,my_max),      hpblx(ix,my_max),                    &
                           hgamt(ix,my_max),     hgamq(ix,my_max)
      real(kind=RTYPE)     psk(ix,km,my_max)
!
      logical dspheat
!          flag for tke dissipative heating
!
!    locals
!
      integer i,iprt,is,iun,k,kk,km1,kmpbl,latd,lond,jj,j1
      integer lcld(ix,my_max),icld(ix,my_max),kcld(ix,my_max),krad(ix,my_max)
      integer kx1(ix,my_max), kpblx(ix,my_max)
!
!     real(kind=kind_phys) betaq(im), betat(im),   betaw(im),
      real(kind=kind_phys) evap(ix,my_max),  heat(ix,my_max),    phih(ix,my_max),            &
                           phim(ix,my_max),  rbdn(ix,my_max),    rbup(ix,my_max),            &
                           stress(ix,my_max),beta(ix,my_max),    sflux(ix,my_max),           &
                           z0(ix,my_max),    crb(ix,my_max),     wstar(ix,my_max),           &
                           zol(ix,my_max),   ustmin(ix,my_max),  ustar(ix,my_max),           &
                           thermal(ix,my_max),wscale(ix,my_max), wscaleu(ix,my_max)
!
      real(kind=kind_phys) theta(ix,km,my_max),thvx(ix,km,my_max),  thlvx(ix,km,my_max),     &
                           qlx(ix,km,my_max),  thetae(ix,km,my_max),                  &
                           qtx(ix,km,my_max),  bf(ix,km-1,my_max),  diss(ix,km,my_max),      &
                           radx(ix,km-1,my_max),                               &
                           govrth(ix,my_max),  hrad(ix,my_max),                       &
!    &                     hradm(im),   radmin(im),   vrad(im),         &
                           radmin(ix,my_max),  vrad(ix,my_max),                       &
                           zd(ix,my_max),      zdd(ix,my_max),      thlvx1(ix,my_max)
!
      real(kind=kind_phys) rdzt(ix,km-1,my_max),dktx(ix,km-1,my_max),                 &
                           zi(ix,km+1,my_max),  zl(ix,km,my_max),    xkzo(ix,km-1,my_max),   &
                           dku(ix,km-1,my_max), dkt(ix,km-1,my_max), xkzmo(ix,km-1,my_max),  &
                           cku(ix,km-1,my_max), ckt(ix,km-1,my_max),                  &
                           ti(ix,km-1,my_max),  shr2(ix,km-1,my_max),                 &
                           al(ix,km-1,my_max),  ad(ix,km,my_max),                     &
                           au(ix,km-1,my_max),  a1(ix,km,my_max),                     &
                           a2(ix,km*ntrac,my_max)
!
      real(kind=kind_phys) tcko(ix,km,my_max),  qcko(ix,km,ntrac,my_max),             &
                           ucko(ix,km,my_max),  vcko(ix,km,my_max),  xmf(ix,km,my_max)
!
      real(kind=kind_phys) prinv(ix,my_max), rent(ix,my_max)
!
      logical  pblflg(ix,my_max), sfcflg(ix,my_max), scuflg(ix,my_max), flg(ix,my_max)
      logical  ublflg(ix,my_max), pcnvflg(ix,my_max)
!
!  pcnvflg: true for convective(strongly unstable) pbl
!  ublflg: true for unstable but not convective(strongly unstable) pbl
!
      real(kind=kind_phys) aphi16,  aphi5,  bvf2,   wfac, &
                           cfac,    conq,   cont,   conw, &
                           dk,      dkmax,  dkmin, &
                           dq1,     dsdz2,  dsdzq,  dsdzt, &
                           dsdzu,   dsdzv, &
                           dsig,    dt2,    dthe1,  dtodsd, &
                           dtodsu,  dw2,    dw2min, g, &
                           gamcrq,  gamcrt, gocp, &
                           gravi,   f0, &
                           prnum,   prmax,  prmin,  pfac,  crbcon, &
                           qmin,    tdzmin, qtend,  crbmin,crbmax, &
                           rbint,   rdt,    rdz,    qlmin, &
                           ri,      rimin,  rl2,    rlam,  rlamun, &
                           rone,    rzero,  sfcfrac, &
                           spdk2,   sri,    zol1,   zolcr, zolcru, &
                           robn,    ttend, &
                           utend,   vk,     vk2, &
                           ust3,    wst3, &
                           vtend,   zfac,   vpert,  cteit, &
                           rentf1,  rentf2, radfac, &
                           zfmin,   zk,     tem,    tem1,  tem2, &
                           xkzm,    xkzmu,  xkzminv, &
                           ptem,    ptem1,  ptem2, tx1(ix,my_max), tx2(ix,my_max),fk
!
      real(kind=kind_phys) moninq_fac
!
      real(kind=kind_phys) zstblmax,h1,     h2,     qlcr,  actei, &
                           cldtime
      real(kind=kind_phys) dtodsu_1,rdz_1,dsig_1,tem1_1,dsdz2_1,al_1, &
                           tem2_1,ptem_1,ptem2_1,dttmp,dqtmp,fk(ix,my_max),fkk(ix,km-2,my_max)
!cc
      parameter(gravi=1.0/grav)
      parameter(g=grav)
      parameter(gocp=g/cp)
      parameter(cont=cp/g,conq=hvap/g,conw=1.0/g)               ! for del in pa
!     parameter(cont=1000.*cp/g,conq=1000.*hvap/g,conw=1000./g) ! for del in kpa
      parameter(rlam=30.0,vk=0.4,vk2=vk*vk)
      parameter(prmin=0.25,prmax=4.,zolcr=0.2,zolcru=-0.5)
      parameter(dw2min=0.0001,dkmin=0.0,dkmax=1000.,rimin=-100.)
      parameter(crbcon=0.25,crbmin=0.15,crbmax=0.35)
      parameter(wfac=7.0,cfac=6.5,pfac=2.0,sfcfrac=0.1)
!     parameter(qmin=1.e-8,xkzm=1.0,zfmin=1.e-8,aphi5=5.,aphi16=16.)
      parameter(qmin=1.e-8,         zfmin=1.e-8,aphi5=5.,aphi16=16.)
      parameter(tdzmin=1.e-3,qlmin=1.e-12,f0=1.e-4)
      parameter(h1=0.33333333,h2=0.66666667)
!     parameter(cldtime=500.,xkzminv=0.3)
      parameter(cldtime=500.)
!     parameter(cldtime=500.,xkzmu=3.0,xkzminv=0.3)
!     parameter(gamcrt=3.,gamcrq=2.e-3,rlamun=150.0)
      parameter(xkzm_m=1.0,xkzm_h=1.0,xkzm_s=1.0,xkzminv=0.3,moninq_fac=1.0)
      parameter(gamcrt=3.,gamcrq=0.,rlamun=150.0)
      parameter(rentf1=0.2,rentf2=1.0,radfac=0.85)
      parameter(iun=84)
!
!     parameter (zstblmax = 2500., qlcr=1.0e-5)
!     parameter (zstblmax = 2500., qlcr=3.0e-5)
!     parameter (zstblmax = 2500., qlcr=3.5e-5)
!     parameter (zstblmax = 2500., qlcr=1.0e-4)
      parameter (zstblmax = 2500., qlcr=3.5e-5)
!     parameter (actei = 0.23)
      parameter (actei = 0.7)
      integer(kind=cuda_stream_kind) stream
      integer istat, async_id

      async_id = -1
      stream = acc_get_cuda_stream(async_id)

!
!byl      dspheat=.false. reference from NCEP fv3GFS, it should be .true.
      dspheat=.true.
!
!c-----------------------------------------------------------------------
!
 601  format(1x,' moninp lat lon step hour ',3i6,f6.1)
 602      format(1x,'    k','        z','        t','       th', &
           '      tvh','        q','        u','        v', &
           '       sp')
 603      format(1x,i5,8f9.1)
 604      format(1x,'  sfc',9x,f9.1,18x,f9.1)
 605      format(1x,'    k      zl    spd2   thekv   the1v' &
               ,' thermal    rbup')
 606      format(1x,i5,6f8.2)
 607      format(1x,' kpbl    hpbl      fm      fh   hgamt', &
               '   hgamq      ws   ustar      cd      ch')
 608      format(1x,i5,9f8.2)
 609      format(1x,' k pr dkt dku ',i5,3f8.2)
 610      format(1x,' k pr dkt dku ',i5,3f8.2,' l2 ri t2', &
               ' sr2  ',2f8.2,2e10.2)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!>  ## Compute preliminary variables from input arguments
!
      !jlistnum = 2

      if (.false.) then
      !!!com133 for testing
      kinver = 0
      tau = 0.
      rtg = 0.
      dusfc = 0.
      dvsfc = 0.
      dtsfc = 0.
      dqsfc = 0.
      hpblx = 0.
      hgamt = 0.
      hgamq = 0.
      lcld = 0
      icld = 0
      kcld = 0
      krad = 0
      kx1 = 0
      kpblx = 0
      phih = 0.
      phim = 0.
      rbdn = 0.
      rbup = 0.
      beta = 0.
      sflux = 0.
      z0 = 0.
      crb = 0.
      wstar = 0.
      zol = 0.
      ustmin = 0.
      ustar = 0.
      thermal = 0.
      wscale = 0.
      wscaleu = 0.
      theta = 0.
      thvx = 0.
      thlvx = 0.
      qlx = 0.
      thetae = 0.
      qtx = 0.
      bf = 0.
      diss = 0.
      radx = 0.
      govrth = 0.
      hrad = 0.
      radmin = 0.
      vrad = 0.
      zd = 0.
      zdd = 0.
      thlvx1 = 0.
      rdzt = 0.
      dktx = 0.
      zi = 0.
      zl = 0.
      xkzo = 0.
      dku = 0.
      dkt = 0.
      xkzmo = 0.
      cku = 0.
      ckt = 0.
      ti = 0.
      shr2 = 0.
      al = 0.
      ad = 0.
      au = 0.
      a1 = 0.
      a2 = 0.
      tcko = 0.
      qcko = 0.
      ucko = 0.
      vcko = 0.
      xmf = 0.
      prinv = 0.
      rent = 0.
      pblflg = .false.
      sfcflg = .false.
      scuflg = .false.
      flg = .false.
      ublflg = .false.
      pcnvflg = .false.
      tx1 = 0.
      tx2 = 0.
      !!! end com133 for testing    
      endif	  


      !$acc enter data create(a1,a2,ad,al,au,beta,bf,ckt,cku,diss, &
      !$acc&                  dusfc,dvsfc,dtsfc,dqsfc,dku,dkt,dktx, &
      !$acc&                  crb,flg,govrth,hgamt,hgamq,fk,fkk, &
      !$acc&                  hpblx,hrad,icld,kcld,krad,kpblx, &
      !$acc&                  kx1,lcld,pblflg,pcnvflg,phih,phim, &
      !$acc&                  prinv,qcko,qlx,qtx,radmin,radx,rent,rbdn, &
      !$acc&                  rbup,rdzt,rtg,sfcflg,scuflg,shr2,sflux, &
      !$acc&                  tau,tcko,theta,thetae,thermal,myim, &
      !$acc&                  thlvx1,thlvx,thvx,ti,ublflg,ucko, &
      !$acc&                  ustar,ustmin,vcko,vrad,wscale,wscaleu,wstar, &
      !$acc&                  xmf,z0,zd,zdd,zi,zl,zol,xkzo,xkzmo,tx1,tx2)
      !$acc enter data copyin(jlist1)
      
      !$acc host_data use_device(rtg)
      istat = cudaMemsetAsync(rtg, 0.0, size(rtg), stream)
      !$acc end host_data

      !rtg   = 0.

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
      dt2   = delt
      rdt   = 1. / dt2
      km1   = km - 1
      kmpbl = km / 2
      
      !$acc parallel loop private(j1,jj)
      do jj = 1, jlistnum
        j1 = jlist1(jj)
        myim(jj) = im(j1)
      enddo
  !>  - Compute physical height of the layer centers and interfaces from the geopotential height (zi and zl)
      !$acc parallel loop gang collapse(2) private(jj,k,i)
      do jj = 1, jlistnum
        do k=1,km
          !$acc loop vector
          do i=1,myim(jj)
            zi(i,k,jj) = phii(i,k) * gravi
            zl(i,k,jj) = phil(i,k,jj) * gravi
          enddo
        enddo
      enddo
      
      !$acc parallel loop gang private(jj,i)
      do jj = 1, jlistnum
        !$acc loop vector
        do i=1,myim(jj)
           zi(i,km+1,jj) = phii(i,km+1) * gravi
        enddo
      enddo
  !>  - Compute reciprocal of \f$ \Delta z \f$ (rdzt)
      !$acc parallel loop gang collapse(2) private(jj,k,i)
      do jj = 1, jlistnum
        do k = 1,km1
          !$acc loop vector
          do i=1,myim(jj)
            rdzt(i,k,jj) = 1.0 / (zl(i,k+1,jj) - zl(i,k,jj))
          enddo
        enddo
      enddo
  !>  - Compute reciprocal of pressure (tx1, tx2)
      !$acc parallel loop gang private(jj,i)
      do jj = 1, jlistnum
        !$acc loop vector
        do i=1,myim(jj)
          kx1(i,jj) = 1
          tx1(i,jj) = 1.0 / prsi(i,1,jj)
          tx2(i,jj) = tx1(i,jj)
        enddo
      enddo
  !>  - Compute background vertical diffusivities for scalars and momentum (xkzo and xkzmo)
      !$acc parallel loop gang private(jj,k,i,ptem,tem1)
      do jj = 1, jlistnum
        !$acc loop vector
        do i=1,myim(jj)
          !$acc loop seq
          do k = 1,km1
            xkzo(i,k,jj)  = 0.0
            xkzmo(i,k,jj) = 0.0
  !byl          if (k < kinver(i)) then
  !                                  vertical background diffusivity
              ptem      = prsi(i,k+1,jj) * tx1(i,jj)
              tem1      = 1.0 - ptem
              tem1      = tem1 * tem1 * 10.0
              xkzo(i,k,jj) = xkzm_h * min(1.0, exp(-tem1))

  !                                  vertical background diffusivity for momentum
              if (ptem >= xkzm_s) then
                xkzmo(i,k,jj) = xkzm_m
                kx1(i,jj)     = k + 1
              else
                if (k == kx1(i,jj) .and. k > 1) tx2(i,jj) = 1.0 / prsi(i,k,jj)
                tem1 = 1.0 - prsi(i,k+1,jj) * tx2(i,jj)
                tem1 = tem1 * tem1 * 5.0
                xkzmo(i,k,jj) = xkzm_m * min(1.0, exp(-tem1))
              endif
  !byl          endif
          enddo
        enddo
      enddo
  !     if (lprnt) then
  !       print *,' xkzo=',(xkzo(ipr,k),k=1,km1)
  !       print *,' xkzmo=',(xkzmo(ipr,k),k=1,km1)
  !     endif
  !
  ! diffusivity in the inversion layer is set to be xkzminv (m^2/s)
  !>  - The background scalar vertical diffusivity is limited to be less than or equal to xkzminv
      !$acc parallel loop gang collapse(2) private(jj,k,i,tem1)
      do jj = 1, jlistnum
        do k = 1,kmpbl
          !$acc loop vector
          do i=1,myim(jj)
  !         if(zi(i,k+1) > 200..and.zi(i,k+1) < zstblmax) then
            if(zi(i,k+1,jj) > 250.) then
              tem1 = (t1(i,k+1,jj)-t1(i,k,jj)) * rdzt(i,k,jj)
              if(tem1 > 1.e-5) then
                 xkzo(i,k,jj)  = min(xkzo(i,k,jj),xkzminv)
              endif
            endif
          enddo
        enddo
      enddo
  !>  - Some output variables and logical flags are initialized
      !$acc parallel loop gang private(jj,i)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1,myim(jj)
           z0(i,jj)    = 0.01 * zorl(i,jj)
           dusfc(i,jj) = 0.
           dvsfc(i,jj) = 0.
           dtsfc(i,jj) = 0.
           dqsfc(i,jj) = 0.
           wscale(i,jj)= 0.
           wscaleu(i,jj)= 0.
           kpbl(i,jj)  = 1
           hpbl(i,jj)  = zi(i,1,jj)
           hpblx(i,jj) = zi(i,1,jj)
           pblflg(i,jj)= .true.
           sfcflg(i,jj)= .true.
           if(rbsoil(i,jj) > 0.) sfcflg(i,jj) = .false.
           ublflg(i,jj)= .false.
           pcnvflg(i,jj)= .false.
           scuflg(i,jj)= .true.
           if(scuflg(i,jj)) then
             radmin(i,jj)= 0.
             rent(i,jj)  = rentf1
             hrad(i,jj)  = zi(i,1,jj)
  !          hradm(i) = zi(i,1)
             krad(i,jj)  = 1
             icld(i,jj)  = 0
             lcld(i,jj)  = km1
             kcld(i,jj)  = km1
             zd(i,jj)    = 0.
          endif
        enddo
      enddo
  !>  - Compute \f$\theta\f$ (theta), \f$q_l\f$ (qlx), \f$q_t\f$ (qtx), \f$\theta_e\f$ (thetae), \f$\theta_v\f$ (thvx), \f$\theta_{l,v}\f$ (thlvx)
      !$acc parallel loop gang collapse(2) private(jj,k,i,ptem, ptem1, ptem2)
      do jj = 1, jlistnum
        do k = 1,km
          !$acc loop vector
          do i = 1,myim(jj)
            theta(i,k,jj) = t1(i,k,jj) * psk(i,km,jj) / prslk(i,k,jj)
            qlx(i,k,jj)   = max(q1(i,k,ntcw,jj),qlmin)
            qtx(i,k,jj)   = max(q1(i,k,1,jj),qmin)+qlx(i,k,jj)
            ptem       = qlx(i,k,jj)
            ptem1      = hvap*max(q1(i,k,1,jj),qmin)/(cp*t1(i,k,jj))
            thetae(i,k,jj)= theta(i,k,jj)*(1.+ptem1)
            thvx(i,k,jj)  = theta(i,k,jj)*(1.+fv*max(q1(i,k,1,jj),qmin)-ptem)
            ptem2      = theta(i,k,jj)-(hvap/cp)*ptem
            thlvx(i,k,jj) = ptem2*(1.+fv*qtx(i,k,jj))
          enddo
        enddo
      enddo
  !>  - Initialize diffusion coefficients to 0 and calculate the total radiative heating rate (dku, dkt, radx)
      !$acc parallel loop gang collapse(2) private(jj,k,i,tem)
      do jj = 1, jlistnum
        do k = 1,km1
          !$acc loop vector
          do i = 1,myim(jj)
            dku(i,k,jj)  = 0.
            dkt(i,k,jj)  = 0.
            dktx(i,k,jj) = 0.
            cku(i,k,jj)  = 0.
            ckt(i,k,jj)  = 0.
            tem       = zi(i,k+1,jj)-zi(i,k,jj)
            radx(i,k,jj) = tem*(swh(i,k,jj)*xmu(i,jj)+hlw(i,k,jj))
          enddo
        enddo
      enddo
  !>  - Set lcld to first index above 2.5km
      !$acc parallel loop gang private(jj,i)
      do jj = 1, jlistnum
        !$acc loop vector
        do i=1,myim(jj)
           flg(i,jj)  = scuflg(i,jj)
        enddo
      enddo
      !$acc parallel loop gang private(jj,k,i)
      do jj = 1, jlistnum
      !$acc loop vector
        do i=1,myim(jj)
          !$acc loop seq
          do k = 1, km1
            if(flg(i,jj).and.zl(i,k,jj) >= zstblmax) then
               lcld(i,jj)=k
               flg(i,jj)=.false.
            endif
          enddo
        enddo
      enddo
  !
  !  compute virtual potential temp gradient (bf) and winshear square
  !>  - Compute \f$\frac{\partial \theta_v}{\partial z}\f$ (bf) and the wind shear squared (shr2)
      !$acc parallel loop gang collapse(2) private(jj,k,i,rdz,dw2)
      do jj = 1, jlistnum
        do k = 1, km1
          !$acc loop vector
          do i = 1, myim(jj)
             rdz  = rdzt(i,k,jj)
             bf(i,k,jj) = (thvx(i,k+1,jj)-thvx(i,k,jj))*rdz
             ti(i,k,jj) = 2./(t1(i,k,jj)+t1(i,k+1,jj))
             dw2  = (u1(i,k,jj)-u1(i,k+1,jj))**2 &
                  + (v1(i,k,jj)-v1(i,k+1,jj))**2
             shr2(i,k,jj) = max(dw2,dw2min)*rdz*rdz
          enddo
        enddo
      enddo
    !>  - Calculate \f$\frac{g}{\theta}\f$ (govrth), \f$\beta = \frac{\Delta t}{\Delta z}\f$ (beta), \f$u_*\f$ (ustar), total surface flux (sflux), and set pblflag to false if the total surface energy flux is into the surface
      !$acc parallel loop gang private(jj,i)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1,myim(jj)
          govrth(i,jj) = g/theta(i,1,jj)
        enddo
      enddo
  !
      !$acc parallel loop gang private(jj,i)
      do jj = 1, jlistnum
        !$acc loop vector
        do i=1,myim(jj)
           beta(i,jj)  = dt2 / (zi(i,2,jj)-zi(i,1,jj))
        enddo
      enddo
  !
      !$acc parallel loop gang private(jj,i)
      do jj = 1, jlistnum
        !$acc loop vector
        do i=1,myim(jj)
           ustar(i,jj) = sqrt(stress(i,jj))
        enddo
      enddo
  !
      !$acc parallel loop gang private(jj,i)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1,myim(jj)
           sflux(i,jj)  = heat(i,jj) + evap(i,jj)*fv*theta(i,1,jj)
           if(.not.sfcflg(i,jj) .or. sflux(i,jj) <= 0.) pblflg(i,jj)=.false.
        enddo
      enddo
  !>  ## Calculate the first estimate of the PBL height (``Predictor step")
  !!  The calculation of the boundary layer height follows Troen and Mahrt (1986) \cite troen_and_mahrt_1986 section 3. The approach is to find the level in the column where a modified bulk Richardson number exceeds a critical value.
  !!
  !!  The temperature of the thermal is of primary importance. For the initial estimate of the PBL height, the thermal is assumed to have one of two temperatures. If the boundary layer is stable, the thermal is assumed to have a temperature equal to the sur
  !  compute the pbl height
  !
      !$acc parallel loop gang private(jj,i,tem,robn,tem1)
      do jj = 1, jlistnum
        !$acc loop vector
        do i=1,myim(jj)
           flg(i,jj) = .false.
           rbup(i,jj) = rbsoil(i,jj)
  !
           if(pblflg(i,jj)) then
             thermal(i,jj) = thvx(i,1,jj)
             crb(i,jj) = crbcon
           else
             thermal(i,jj) = tsea(i,jj)*(1.+fv*max(q1(i,1,1,jj),qmin))
             tem = sqrt(u10m(i,jj)**2+v10m(i,jj)**2)
             tem = max(tem, 1.)
             robn = tem / (f0 * z0(i,jj))
             tem1 = 1.e-7 * robn
             crb(i,jj) = 0.16 * (tem1 ** (-0.18))
             crb(i,jj) = max(min(crb(i,jj), crbmax), crbmin)
           endif
        enddo
      enddo
  !>  Given the thermal's properties and the critical Richardson number, a loop is executed to find the first level above the surface where the modified Richardson number is greater than the critical Richardson number, using equation 10a from Troen and Mahr
  !!  \f[
  !!  h = Ri\frac{T_0\left|\vec{v}(h)\right|^2}{g\left(\theta_v(h) - \theta_s\right)}
  !!  \f]
  !!  where \f$h\f$ is the PBL height, \f$Ri\f$ is the Richardson number, \f$T_0\f$ is the virtual potential temperature near the surface, \f$\left|\vec{v}\right|\f$ is the wind speed, and \f$\theta_s\f$ is for the thermal. Rearranging this equation to calc
  !!  \f[
  !!  Ri_k = gz(k)\frac{\left(\theta_v(k) - \theta_s\right)}{\theta_v(1)*\vec{v}(k)}
  !!  \f]
      !$acc parallel loop gang private(jj,k,i,spdk2)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1, myim(jj)
          !$acc loop seq
          do k = 1, kmpbl
            if(.not.flg(i,jj)) then
              rbdn(i,jj) = rbup(i,jj)
              spdk2   = max((u1(i,k,jj)**2+v1(i,k,jj)**2),1.)
              rbup(i,jj) = (thvx(i,k,jj)-thermal(i,jj))* &
                        (g*zl(i,k,jj)/thvx(i,1,jj))/spdk2
              kpbl(i,jj) = k
              flg(i,jj)  = rbup(i,jj) > crb(i,jj)
            endif
          enddo
        enddo
      enddo
  !>  Once the level is found, some linear interpolation is performed to find the exact height of the boundary layer top (where \f$Ri = Ri_{cr}\f$) and the PBL height and the PBL top index are saved (hpblx and kpblx, respectively)
      !$acc parallel loop gang private(jj,i,k,rbint)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1,myim(jj)
          if(kpbl(i,jj) > 1) then
            k = kpbl(i,jj)
            if(rbdn(i,jj) >= crb(i,jj)) then
              rbint = 0.
            elseif(rbup(i,jj) <= crb(i,jj)) then
              rbint = 1.
            else
              rbint = (crb(i,jj)-rbdn(i,jj))/(rbup(i,jj)-rbdn(i,jj))
            endif
            hpbl(i,jj) = zl(i,k-1,jj) + rbint*(zl(i,k,jj)-zl(i,k-1,jj))
            if(hpbl(i,jj) < zi(i,kpbl(i,jj),jj)) kpbl(i,jj) = kpbl(i,jj) - 1
          else
            hpbl(i,jj) = zl(i,1,jj)
            kpbl(i,jj) = 1
          endif
          kpblx(i,jj) = kpbl(i,jj)
          hpblx(i,jj) = hpbl(i,jj)
        enddo
      enddo
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
      !$acc parallel loop gang private(i,jj,zol1,tem)
      do jj = 1, jlistnum
        !$acc loop vector
        do i=1,myim(jj)
           zol(i,jj) = max(rbsoil(i,jj)*fm(i,jj)*fm(i,jj)/fh(i,jj),rimin)
           if(sfcflg(i,jj)) then
             zol(i,jj) = min(zol(i,jj),-zfmin)
           else
             zol(i,jj) = max(zol(i,jj),zfmin)
           endif
           zol1 = zol(i,jj)*sfcfrac*hpbl(i,jj)/zl(i,1,jj)
           if(sfcflg(i,jj)) then
  !          phim(i) = (1.-aphi16*zol1)**(-1./4.)
  !          phih(i) = (1.-aphi16*zol1)**(-1./2.)
             tem     = 1.0 / (1. - aphi16*zol1)
             phih(i,jj) = sqrt(tem)
             phim(i,jj) = sqrt(phih(i,jj))
           else
             phim(i,jj) = 1. + aphi5*zol1
             phih(i,jj) = phim(i,jj)
           endif
           wscale(i,jj) = ustar(i,jj)/phim(i,jj)
           ustmin(i,jj) = ustar(i,jj)/aphi5
           wscale(i,jj) = max(wscale(i,jj),ustmin(i,jj))
        enddo
      enddo
      !$acc parallel loop gang private(jj,i,wst3,ust3)
      do jj = 1, jlistnum
        !$acc loop vector
        do i=1,myim(jj)
          if(pblflg(i,jj)) then
            if(zol(i,jj) < zolcru .and. kpbl(i,jj) > 1) then
              pcnvflg(i,jj) = .true.
            else
              ublflg(i,jj) = .true.
            endif
            wst3 = govrth(i,jj)*sflux(i,jj)*hpbl(i,jj)
            wstar(i,jj)= wst3**h1
            ust3 = ustar(i,jj)**3.
            wscaleu(i,jj) = (ust3+wfac*vk*wst3*sfcfrac)**h1
            wscaleu(i,jj) = max(wscaleu(i,jj),ustmin(i,jj))
          endif
        enddo
      enddo
  !
  ! compute counter-gradient mixing term for heat and moisture
  !>  ## Update thermal properties of surface parcel and recompute PBL height ("Corrector step").
  !!  Next, the counter-gradient terms for temperature and humidity are calculated using equation 4 of Hong and Pan (1996) \cite hong_and_pan_1996 and are used to calculate the "scaled virtual temperature excess near the surface" (equation 9 in Hong and Pan
      !$acc parallel loop gang private(jj,i,vpert)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1,myim(jj)
           if(ublflg(i,jj)) then
             hgamt(i,jj)  = min(cfac*heat(i,jj)/wscaleu(i,jj),gamcrt)
             hgamq(i,jj)  = min(cfac*evap(i,jj)/wscaleu(i,jj),gamcrq)
             vpert     = hgamt(i,jj) + hgamq(i,jj)*fv*theta(i,1,jj)
             vpert     = min(vpert,gamcrt)
             thermal(i,jj)= thermal(i,jj)+max(vpert,0.)
             hgamt(i,jj)  = max(hgamt(i,jj),0.0)
             hgamq(i,jj)  = max(hgamq(i,jj),0.0)
           endif
        enddo
      enddo
  !
  !  enhance the pbl height by considering the thermal excess
  !>  The PBL height calculation follows the same procedure as the predictor step, except that it uses an updated virtual potential temperature for the thermal.
      !$acc parallel loop gang private(jj,i)
      do jj = 1, jlistnum
        !$acc loop vector
        do i=1,myim(jj)
           flg(i,jj)  = .true.
           if(ublflg(i,jj)) then
             flg(i,jj)  = .false.
             rbup(i,jj) = rbsoil(i,jj)
           endif
        enddo
      enddo
      !$acc parallel loop gang private(jj,k,i,spdk2)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1, myim(jj)
          !$acc loop seq
          do k = 2, kmpbl
            if(.not.flg(i,jj)) then
              rbdn(i,jj) = rbup(i,jj)
              spdk2   = max((u1(i,k,jj)**2+v1(i,k,jj)**2),1.)
              rbup(i,jj) = (thvx(i,k,jj)-thermal(i,jj))* &
                        (g*zl(i,k,jj)/thvx(i,1,jj))/spdk2
              kpbl(i,jj) = k
              flg(i,jj)  = rbup(i,jj) > crb(i,jj)
            endif
          enddo
        enddo
      enddo
      !$acc parallel loop gang private(jj,i,k,rbint)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1,myim(jj)
          if(ublflg(i,jj)) then
             k = kpbl(i,jj)
             if(rbdn(i,jj) >= crb(i,jj)) then
                rbint = 0.
             elseif(rbup(i,jj) <= crb(i,jj)) then
                rbint = 1.
             else
                rbint = (crb(i,jj)-rbdn(i,jj))/(rbup(i,jj)-rbdn(i,jj))
             endif
             hpbl(i,jj) = zl(i,k-1,jj) + rbint*(zl(i,k,jj)-zl(i,k-1,jj))
             if(hpbl(i,jj) < zi(i,kpbl(i,jj),jj)) kpbl(i,jj) = kpbl(i,jj) - 1
             if(kpbl(i,jj) <= 1) then
                ublflg(i,jj) = .false.
                pblflg(i,jj) = .false.
             endif
          endif
        enddo
      enddo
  !
  !  look for stratocumulus
  !>  ## Determine whether stratocumulus layers exist and compute quantities needed for enhanced diffusion
  !!  - Starting at the PBL top and going downward, if the level is less than 2.5 km and \f$q_l>q_{l,cr}\f$ then set kcld = k (find the cloud top index in the PBL). If no cloud water above the threshold is found, scuflg is set to F.
      !$acc parallel loop gang private(jj,i)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1, myim(jj)
          flg(i,jj)=scuflg(i,jj)
        enddo
      enddo
      
      !$acc parallel loop gang private(jj,k,i)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1, myim(jj)
          !$acc loop seq
          do k = kmpbl,1,-1
            if(flg(i,jj) .and. k <= lcld(i,jj)) then
              if(qlx(i,k,jj).ge.qlcr) then
                 kcld(i,jj)=k
                 flg(i,jj)=.false.
              endif
            endif
          enddo
        enddo
      enddo
      !$acc parallel loop gang private(jj,i)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1, myim(jj)
          if(scuflg(i,jj) .and. kcld(i,jj)==km1) scuflg(i,jj)=.false.
        enddo
      enddo
  !>  - Starting at the PBL top and going downward, if the level is less than the cloud top, find the level of the minimum radiative heating rate within the cloud. If the level of the minimum is the lowest model level or the minimum radiative heating rate i
      !$acc parallel loop gang private(jj,i)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1, myim(jj)
          flg(i,jj)=scuflg(i,jj)
        enddo
      enddo
      !$acc parallel loop gang private(jj,k,i)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1, myim(jj)
          !$acc loop seq
          do k = kmpbl,1,-1
            if(flg(i,jj) .and. k <= kcld(i,jj)) then
              if(qlx(i,k,jj) >= qlcr) then
                if(radx(i,k,jj) < radmin(i,jj)) then
                  radmin(i,jj)=radx(i,k,jj)
                  krad(i,jj)=k
                endif
              else
                flg(i,jj)=.false.
              endif
            endif
          enddo
        enddo
      enddo
      !$acc parallel loop gang private(jj,i)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1, myim(jj)
          if(scuflg(i,jj) .and. krad(i,jj) <= 1) scuflg(i,jj)=.false.
          if(scuflg(i,jj) .and. radmin(i,jj)>=0.) scuflg(i,jj)=.false.
        enddo
      enddo
  !>  - Starting at the PBL top and going downward, count the number of levels below the minimum radiative heating rate level that have cloud water above the threshold. If there are none, then set the scuflg to F.
      !$acc parallel loop gang private(jj,i)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1, myim(jj)
          flg(i,jj)=scuflg(i,jj)
        enddo
      enddo
      !$acc parallel loop gang collapse(2) private(jj,k,i)
      do jj = 1, jlistnum
        do k = kmpbl,2,-1
          !$acc loop vector
          do i = 1, myim(jj)
            if(flg(i,jj) .and. k <= krad(i,jj)) then
              if(qlx(i,k,jj) >= qlcr) then
                !$acc atomic
                icld(i,jj)=icld(i,jj)+1
              else
                flg(i,jj)=.false.
              endif
            endif
          enddo
        enddo
      enddo
      !$acc parallel loop gang private(jj,i)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1, myim(jj)
          if(scuflg(i,jj) .and. icld(i,jj) < 1) scuflg(i,jj)=.false.
        enddo
      enddo
  !>  - Find the height of the interface where the minimum in radiative heating rate is located. If this height is less than the second model interface height, then set the scuflg to F.
      !$acc parallel loop gang private(jj,i)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1, myim(jj)
          if(scuflg(i,jj)) then
             hrad(i,jj) = zi(i,krad(i,jj)+1,jj)
  !          hradm(i)= zl(i,krad(i))
          endif
        enddo
      enddo
  !
      !$acc parallel loop gang private(jj,i)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1, myim(jj)
          if(scuflg(i,jj) .and. hrad(i,jj)<zi(i,2,jj)) scuflg(i,jj)=.false.
        enddo
      enddo
  !>  - Calculate the hypothetical \f$\theta_v\f$ at the minimum radiative heating level that a parcel would reach due to radiative cooling after a typical cloud turnover time spent at that level.
      !$acc parallel loop gang private(jj,i,k,tem,tem1)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1, myim(jj)
          if(scuflg(i,jj)) then
          !if (1) then
            k    = krad(i,jj)
            tem  = zi(i,k+1,jj)-zi(i,k,jj)
            tem1 = cldtime*radmin(i,jj)/tem
            thlvx1(i,jj) = thlvx(i,k,jj)+tem1
  !         if(thlvx1(i) > thlvx(i,k-1)) scuflg(i)=.false.
          endif
        enddo
      enddo
  !>  - Determine the distance that a parcel would sink downwards starting from the level of minimum radiative heating rate by comparing the hypothetical minimum \f$\theta_v\f$ calculated above with the environmental \f$\theta_v\f$.
      !$acc parallel loop gang private(jj,i,k,tem,tem1)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1, myim(jj)
           flg(i,jj)=scuflg(i,jj)
        enddo
      enddo
      !$acc parallel loop gang private(jj,k,i,tem)
      do jj = 1, jlistnum  
        !$acc loop vector      
        do i = 1, myim(jj)
          !$acc loop seq
          do k = kmpbl,1,-1
            if(flg(i,jj) .and. k <= krad(i,jj))then
            !if (1) then
              if(thlvx1(i,jj) <= thlvx(i,k,jj))then
                 tem=zi(i,k+1,jj)-zi(i,k,jj)
                 zd(i,jj)=zd(i,jj)+tem
              else
                 flg(i,jj)=.false.
              endif
            endif
          enddo
        enddo
      enddo
  !>  - Calculate the cloud thickness, where the cloud top is the in-cloud minimum radiative heating level and the bottom is determined previously.
      !$acc parallel loop gang private(jj,i,kk)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1, myim(jj)
          if(scuflg(i,jj))then
          !if (1) then
            kk = max(1, krad(i,jj)+1-icld(i,jj))
            zdd(i,jj) = hrad(i,jj)-zi(i,kk,jj)
          endif
        enddo
      enddo
  !>  - Find the largest between the cloud thickness and the distance of a sinking parcel, then determine the smallest of that number and the height of the minimum in radiative heating rate. Set this number to \f$zd\f$. Using \f$zd\f$, calculate the charact
      !$acc parallel loop gang private(jj,i,tem)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1, myim(jj)
          if(scuflg(i,jj))then
          !if (1) then
            zd(i,jj) = max(zd(i,jj),zdd(i,jj))
            zd(i,jj) = min(zd(i,jj),hrad(i,jj))
            tem   = govrth(i,jj)*zd(i,jj)*(-radmin(i,jj))
            vrad(i,jj)= tem**h1
          endif
        enddo
      enddo
  !
  !     compute inverse prandtl number
  !>  ## Calculate the inverse Prandtl number
  !!  For an unstable PBL, the Prandtl number is calculated according to Hong and Pan (1996) \cite hong_and_pan_1996, equation 10, whereas for a stable boundary layer, the Prandtl number is simply \f$Pr = \frac{\phi_h}{\phi_m}\f$.
      !$acc parallel loop gang private(jj,i,tem)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1, myim(jj)
          if(ublflg(i,jj)) then
            tem = phih(i,jj)/phim(i,jj)+cfac*vk*sfcfrac
          else
            tem = phih(i,jj)/phim(i,jj)
          endif
          prinv(i,jj) =  1.0 / tem
          prinv(i,jj) = min(prinv(i,jj),prmax)
          prinv(i,jj) = max(prinv(i,jj),prmin)
        enddo
      enddo
      !$acc parallel loop gang private(jj,i)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1, myim(jj)
          if(zol(i,jj) > zolcr) then
            kpbl(i,jj) = 1
          endif
        enddo
      enddo
  !
  !     compute diffusion coefficients below pbl
  !>  ## Compute diffusion coefficients below the PBL top
  !!  Below the PBL top, the diffusion coefficients (\f$K_m\f$ and \f$K_h\f$) are calculated according to equation 2 in Hong and Pan (1996) \cite hong_and_pan_1996 where a different value for \f$w_s\f$ (PBL vertical velocity scale) is used depending on the
      !$acc parallel loop gang collapse(2) private(jj,k,i,zfac,tem,tem1)
      do jj = 1, jlistnum
        do k = 1, kmpbl
          !$acc loop vector
          do i=1,myim(jj)
             if(k < kpbl(i,jj)) then
    !           zfac = max((1.-(zi(i,k+1)-zl(i,1))/
    !    1             (hpbl(i)-zl(i,1))), zfmin)
                zfac = max((1.-zi(i,k+1,jj)/hpbl(i,jj)), zfmin)
                tem = zi(i,k+1,jj) * (zfac**pfac) * moninq_fac ! lmh suggested by kg
                if(pblflg(i,jj)) then
                  tem1 = vk * wscaleu(i,jj) * tem
    !             dku(i,k) = xkzmo(i,k) + tem1
    !             dkt(i,k) = xkzo(i,k)  + tem1 * prinv(i)
                  dku(i,k,jj) = tem1
                  dkt(i,k,jj) = tem1 * prinv(i,jj)
                else
                  tem1 = vk * wscale(i,jj) * tem
    !             dku(i,k) = xkzmo(i,k) + tem1
    !             dkt(i,k) = xkzo(i,k)  + tem1 * prinv(i)
                  dku(i,k,jj) = tem1
                  dkt(i,k,jj) = tem1 * prinv(i,jj)
                endif
                dku(i,k,jj) = min(dku(i,k,jj),dkmax)
                dku(i,k,jj) = max(dku(i,k,jj),xkzmo(i,k,jj))
                dkt(i,k,jj) = min(dkt(i,k,jj),dkmax)
                dkt(i,k,jj) = max(dkt(i,k,jj),xkzo(i,k,jj))
                dktx(i,k,jj)= dkt(i,k,jj)
             endif
          enddo
        enddo
      enddo
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
      !$acc&         private(jj,k,i,bvf2,ri,zk,rl2,dk,sri,tem1,prnum)
      do jj = 1, jlistnum
        do k = 1, km1
          !$acc loop vector
          do i=1,myim(jj)
            if(k >= kpbl(i,jj)) then
              bvf2 = g*bf(i,k,jj)*ti(i,k,jj)
              ri   = max(bvf2/shr2(i,k,jj),rimin)
              zk   = vk*zi(i,k+1,jj)
              if(ri < 0.) then ! unstable regime
                rl2      = zk*rlamun/(rlamun+zk)
                dk       = rl2*rl2*sqrt(shr2(i,k,jj))
                sri      = sqrt(-ri)
!               dku(i,k) = xkzmo(i,k) + dk*(1+8.*(-ri)/(1+1.746*sri))
!               dkt(i,k) = xkzo(i,k)  + dk*(1+8.*(-ri)/(1+1.286*sri))
                dku(i,k,jj) = dk*(1+8.*(-ri)/(1+1.746*sri))
                dkt(i,k,jj) = dk*(1+8.*(-ri)/(1+1.286*sri))
              else             ! stable regime
                rl2      = zk*rlam/(rlam+zk)
!!              tem      = rlam * sqrt(0.01*prsi(i,k))
!!              rl2      = zk*tem/(tem+zk)
                dk       = rl2*rl2*sqrt(shr2(i,k,jj))
                tem1     = dk/(1+5.*ri)**2
!
                if(k >= kpblx(i,jj)) then
                  prnum = 1.0 + 2.1*ri
                  prnum = min(prnum,prmax)
                else
                  prnum = 1.0
                endif
!               dku(i,k) = xkzmo(i,k) + tem1 * prnum
!               dkt(i,k) = xkzo(i,k)  + tem1
                dku(i,k,jj) = tem1 * prnum
                dkt(i,k,jj) = tem1
              endif
!
              dku(i,k,jj) = min(dku(i,k,jj),dkmax)
              dku(i,k,jj) = max(dku(i,k,jj),xkzmo(i,k,jj))
              dkt(i,k,jj) = min(dkt(i,k,jj),dkmax)
              dkt(i,k,jj) = max(dkt(i,k,jj),xkzo(i,k,jj))
!
            endif
!
          enddo
        enddo
      enddo
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  compute components for mass flux mixing by large thermals
  !>  ## If the PBL is convective, call the mass flux scheme to replace the countergradient terms.
  !!  If the PBL is convective, the updraft properties are initialized to be the same as the state variables and the subroutine mfpbl is called.
      !$acc parallel loop gang collapse(2) private(jj,k,i)
      do jj = 1, jlistnum
        do k = 1, km
          !$acc loop vector
          do i = 1, myim(jj)
            if(pcnvflg(i,jj)) then
              tcko(i,k,jj) = t1(i,k,jj)
              ucko(i,k,jj) = u1(i,k,jj)
              vcko(i,k,jj) = v1(i,k,jj)
              xmf(i,k,jj) = 0.
            endif
          enddo
        enddo
      enddo
      !$acc parallel loop gang collapse(3) private(jj,kk,k,i)
      do jj = 1, jlistnum
        do kk = 1, ntrac
          do k = 1, km
            !$acc loop vector
            do i = 1, myim(jj)
              if(pcnvflg(i,jj)) then
                qcko(i,k,kk,jj) = q1(i,k,kk,jj)
              endif
            enddo
          enddo
        enddo
      enddo
  !>  For details of the mfpbl subroutine, step into its documentation ::mfpbl
      call mfpbl_gpu(myim,ix,km,ntrac,dt2,pcnvflg, &
               zl,zi,thvx,q1,t1,u1,v1,hpbl,kpbl, &
               sflux,ustar,wstar,xmf,tcko,qcko,ucko,vcko)

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
      !$acc parallel loop gang private(jj,i,k,tem,tem1,cteit)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1, myim(jj)
          if(scuflg(i,jj)) then
          !if (1) then
             k = krad(i,jj)
             tem = thetae(i,k,jj) - thetae(i,k+1,jj)
             tem1 = qtx(i,k,jj) - qtx(i,k+1,jj)
             if (tem > 0. .and. tem1 > 0.) then
               cteit= cp*tem/(hvap*tem1)
               if(cteit > actei) rent(i,jj) = rentf2
             endif
          endif
        enddo
      enddo
      !$acc parallel loop gang private(jj,i,k,tem1)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1, myim(jj)
          if(scuflg(i,jj)) then
          !if (1) then
             k = krad(i,jj)
             tem1  = max(bf(i,k,jj),tdzmin)
             ckt(i,k,jj) = -rent(i,jj)*radmin(i,jj)/tem1
             cku(i,k,jj) = ckt(i,k,jj)
          endif
        enddo
      enddo
  !
      !$acc parallel loop gang &
      !$acc&              private(jj,k,i,tem1,tem2,ptem)
      do jj = 1, jlistnum
	    !$acc loop vector
        do i=1,myim(jj)
          !$acc loop seq
           do k = 1, kmpbl
              if(scuflg(i,jj) .and. k < krad(i,jj)) then
              !if (1) then    !maybe have some error
                 tem1=hrad(i,jj)-zd(i,jj)
                 tem2=zi(i,k+1,jj)-tem1
                 if(tem2 > 0.) then
                    ptem= tem2/zd(i,jj)
                    if(ptem.ge.1.) ptem= 1.
                    ptem= tem2*ptem*sqrt(1.-ptem)
                    ckt(i,k,jj) = radfac*vk*vrad(i,jj)*ptem
                    cku(i,k,jj) = 0.75*ckt(i,k,jj)
                    ckt(i,k,jj) = max(ckt(i,k,jj),dkmin)
                    ckt(i,k,jj) = min(ckt(i,k,jj),dkmax)
                    cku(i,k,jj) = max(cku(i,k,jj),dkmin)
                    cku(i,k,jj) = min(cku(i,k,jj),dkmax)
                 endif
              endif
           enddo
        enddo
      enddo
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !>  After \f$K_h^{Sc}\f$ has been determined from the surface to the top of the stratocumulus layer, it is added to the value for the diffusion coefficient calculated previously using surface-based mixing [see equation 6 of Lock et al. (2000) \cite lock_e
      !$acc parallel loop gang collapse(2) private(jj,k,i)
      do jj = 1, jlistnum
        do k = 1, kmpbl
          !$acc loop vector
          do i=1,myim(jj)
            if(scuflg(i,jj)) then
            !if (1) then
               dkt(i,k,jj) = dkt(i,k,jj)+ckt(i,k,jj)
               dku(i,k,jj) = dku(i,k,jj)+cku(i,k,jj)
               dkt(i,k,jj) = min(dkt(i,k,jj),dkmax)
               dku(i,k,jj) = min(dku(i,k,jj),dkmax)
            endif
          enddo
        enddo
      enddo
  !
  !     compute tridiagonal matrix elements for heat and moisture
  !
  !>  ## Solve for the temperature and moisture tendencies due to vertical mixing.
  !!  The tendencies of heat, moisture, and momentum due to vertical diffusion are calculated using a two-part process. First, a solution is obtained using an implicit time-stepping scheme, then the time tendency terms are "backed out". The tridiagonal matr

        if(ntrac >= 2) then
          !$acc parallel loop gang collapse(2) private(jj,k,i,is)
          do jj = 1, jlistnum
            do k = 2, ntrac
              is = (k-1) * km
              !$acc loop vector
              do i = 1, myim(jj)
                a2(i,1+is,jj) = q1(i,1,k,jj)
              enddo
            enddo
          enddo
        endif
        
      !$acc parallel loop gang private(jj,i)
      do jj = 1, jlistnum
        !$acc loop vector
        do i=1,myim(jj)
           ad(i,1,jj) = 1.
           a1(i,1,jj) = t1(i,1,jj)   + beta(i,jj) * heat(i,jj)
           a2(i,1,jj) = q1(i,1,1,jj) + beta(i,jj) * evap(i,jj)
           !if (jj .eq. 5 .and. i .eq. 3) write(*,*) 'a3',beta(i,jj), heat(i,jj),a1(i,1,jj)

        enddo
      enddo

        
      !$acc parallel loop gang private(dtodsd,dtodsu,dsig,rdz, &
      !$acc&                   tem1,dsdz2,tem2,ptem,ptem1,ptem2,&
      !$acc&                   dsdzt,dsdzq,tem,jj,i,k)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1,myim(jj)
          !$acc loop seq
          do k = 1,km1
            dtodsd = dt2/del(i,k,jj)
            dtodsu = dt2/del(i,k+1,jj)
            dsig   = prsl(i,k,jj)-prsl(i,k+1,jj)
            rdz    = rdzt(i,k,jj)
            tem1   = dsig * dkt(i,k,jj) * rdz
            dsdz2     = tem1 * rdz
            au(i,k,jj)   = -dtodsd*dsdz2
            al(i,k,jj)   = -dtodsu*dsdz2
  !
            if(pcnvflg(i,jj) .and. k < kpbl(i,jj)) then
               tem2      = dsig * rdz
               ptem      = 0.5 * tem2 * xmf(i,k,jj)
               ptem1     = dtodsd * ptem
               ptem2     = dtodsu * ptem
               !if (jj .eq. 5 .and. i .eq. 3 .and. k .eq. 1) write(*,*) 'a10',ptem,tem2,xmf(i,k,jj)
               !if (k .gt. 1) then
               !  dtodsu_1 = dt2/del(i,k,jj)
               !  rdz_1    = rdzt(i,k-1,jj)
               !  dsig_1   = prsl(i,k-1,jj)-prsl(i,k,jj)
               !  tem1_1   = dsig_1 * dkt(i,k-1,jj) * rdz_1
               !  dsdz2_1 = tem1_1 * rdz_1
               !  al_1 = -dtodsu_1*dsdz2_1
               !  tem2_1 = dsig_1 * rdz_1
               !  ptem_1 = 0.5 * tem2_1 * xmf(i,k-1,jj)
               !  ptem2_1 = dtodsu_1 * ptem_1
               !  write(*,*) 'a',jj,i,k,al_1,ptem2_1
               !endif
               ad(i,k,jj)   = ad(i,k,jj)-au(i,k,jj)-ptem1
               ad(i,k+1,jj) = 1.-al(i,k,jj)+ptem2
               !write(*,*) 'b',jj,i,k,al(i,k,jj),ptem2
               au(i,k,jj)   = au(i,k,jj)-ptem1
               al(i,k,jj)   = al(i,k,jj)+ptem2
               ptem      = tcko(i,k,jj) + tcko(i,k+1,jj)
               dsdzt     = tem1 * gocp
               !if (jj .eq. 5 .and. i .eq. 3 .and. k .eq. 1) write(*,*) 'a1',a1(i,k+1,jj), a1(i,k,jj)
               !if (jj .eq. 5 .and. i .eq. 3 .and. k .eq. 1) write(*,*) 'a00',dtodsd,dsdzt !ok
               !if (jj .eq. 5 .and. i .eq. 3 .and. k .eq. 1) write(*,*) 'a01',ptem1,ptem
               !if (jj .eq. 5 .and. i .eq. 3 .and. k .eq. 1) write(*,*) 'a02',dtodsu,ptem2
               a1(i,k,jj)   = a1(i,k,jj)+dtodsd*dsdzt-ptem1*ptem
               a1(i,k+1,jj) = t1(i,k+1,jj)-dtodsu*dsdzt+ptem2*ptem
               ptem      = qcko(i,k,1,jj) + qcko(i,k+1,1,jj)
               a2(i,k,jj)   = a2(i,k,jj) - ptem1 * ptem
               a2(i,k+1,jj) = q1(i,k+1,1,jj) + ptem2 * ptem
               !if (jj .eq. 5 .and. i .eq. 3 .and. k .eq. 1) write(*,*) 'a2',a1(i,k+1,jj), a1(i,k,jj)
            elseif(ublflg(i,jj) .and. k < kpbl(i,jj)) then
               ptem1 = dsig * dktx(i,k,jj) * rdz
               tem   = 1.0 / hpbl(i,jj)
               dsdzt = tem1 * gocp - ptem1 * hgamt(i,jj) * tem
               dsdzq = - ptem1 * hgamq(i,jj) * tem
               ad(i,k,jj)   = ad(i,k,jj)-au(i,k,jj)
               ad(i,k+1,jj) = 1.-al(i,k,jj)
               a1(i,k,jj)   = a1(i,k,jj)+dtodsd*dsdzt
               a1(i,k+1,jj) = t1(i,k+1,jj)-dtodsu*dsdzt
               a2(i,k,jj)   = a2(i,k,jj)+dtodsd*dsdzq
               a2(i,k+1,jj) = q1(i,k+1,1,jj)-dtodsu*dsdzq
            else
               ad(i,k,jj)   = ad(i,k,jj)-au(i,k,jj)
               ad(i,k+1,jj) = 1.-al(i,k,jj)
               dsdzt     = tem1 * gocp
               a1(i,k,jj)   = a1(i,k,jj)+dtodsd*dsdzt
               a1(i,k+1,jj) = t1(i,k+1,jj)-dtodsu*dsdzt
               a2(i,k+1,jj) = q1(i,k+1,1,jj)
            endif
  !
               !if (jj .eq. 5 .and. i .eq. 3 .and. k .eq. 1) write(*,*) 'a6',t1(i,k+1,jj),q1(i,k+1,1,jj)

          enddo
        enddo
      enddo

  !
        if(ntrac >= 2) then
          !$acc parallel loop gang collapse(2) private(jj,kk,i,k,is, &
          !$acc&                   dtodsd,dtodsu,dsig,tem,ptem1,ptem2,&
          !$acc&                   tem1)
          do jj = 1, jlistnum
            do kk = 2, ntrac
              is = (kk-1) * km
              !$acc loop vector
              do i = 1, myim(jj)
                !$acc loop seq
                do k = 1, km1
                  if(pcnvflg(i,jj) .and. k < kpbl(i,jj)) then
                    dtodsd = dt2/del(i,k,jj)
                    dtodsu = dt2/del(i,k+1,jj)
                    dsig  = prsl(i,k,jj)-prsl(i,k+1,jj)
                    tem   = dsig * rdzt(i,k,jj)
                    ptem  = 0.5 * tem * xmf(i,k,jj)
                    ptem1 = dtodsd * ptem
                    ptem2 = dtodsu * ptem
                    tem1  = qcko(i,k,kk,jj) + qcko(i,k+1,kk,jj)
                    a2(i,k+is,jj) = a2(i,k+is,jj) - ptem1*tem1
                    a2(i,k+1+is,jj)= q1(i,k+1,kk,jj) + ptem2*tem1
                  else
                    a2(i,k+1+is,jj) = q1(i,k+1,kk,jj)
                  endif
                enddo
              enddo
            enddo
          enddo
        endif
  !
  !     solve tridiagonal problem for heat and moisture
  !
  !>  The tridiagonal system is solved by calling the internal ::tridin subroutine.
      !write(*,*) 'a4',a1(3,2,5)
      !$acc parallel loop gang
      do jj = 1, jlistnum
        call tridin_gpu(ix,myim(jj),km,ntrac,al(1,1,jj),ad(1,1,jj),&
                        au(1,1,jj),a1(1,1,jj),a2(1,1,jj),&
                        au(1,1,jj),a1(1,1,jj),a2(1,1,jj),fk(1,jj),fkk(1,1,jj))
      enddo

  !
  !     recover tendencies of heat and moisture
  !
  !>  After returning with the solution, the tendencies for temperature and moisture are recovered.
      
      !$acc parallel loop gang private(jj,i,k,ttend,qtend,dttmp,dqtmp)&
      !$acc&         reduction(+:dttmp) reduction(+:dqtmp)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1,myim(jj)
          dttmp = 0.
          dqtmp = 0.
          !$acc loop seq
          do k = 1,km
            ttend      = (a1(i,k,jj)-t1(i,k,jj)) * rdt
            qtend      = (a2(i,k,jj)-q1(i,k,1,jj))*rdt
!!            tau(i,k)   = tau(i,k)+ttend
!!            rtg(i,k,1) = rtg(i,k,1)+qtend
            tau(i,k,jj)   = ttend
            rtg(i,k,1,jj) = qtend
            dttmp = dttmp + cont*del(i,k,jj)*ttend
            dqtmp = dqtmp + conq*del(i,k,jj)*qtend
            q1(i,k,1,jj)  = a2(i,k,jj)
            !if (jj .eq. 5 .and. i .eq. 3 .and. k .eq. 2) write(*,*) 'a5',a1(i,k,jj),t1(i,k,jj)
          enddo
          dtsfc(i,jj)   = dtsfc(i,jj) + dttmp
          dqsfc(i,jj)   = dqsfc(i,jj) + dqtmp
        enddo
      enddo
        if(ntrac >= 2) then
          !$acc parallel loop gang collapse(3) private(jj,kk,k,i,is,qtend)
          do jj = 1, jlistnum
            do kk = 2, ntrac
              do k = 1, km
                is = (kk-1) * km
                !$acc loop vector
                do i = 1, myim(jj)
                  qtend = (a2(i,k+is,jj)-q1(i,k,kk,jj))*rdt
                  rtg(i,k,kk,jj) = rtg(i,k,kk,jj)+qtend
                  q1(i,k,kk,jj)  = a2(i,k+is,jj)
                enddo
              enddo
            enddo
          enddo
        endif
  !
  !   compute tke dissipation rate
  !
  !>  ## Calculate heating due to TKE dissipation and add to the tendency for temperature
  !!  Following Han et al. (2015) \cite han_et_al_2015 , turbulence dissipation contributes to the tendency of temperature in the following way. First, turbulence dissipation is calculated by equation 17 of Han et al. (2015) \cite han_et_al_2015 for the PBL
      
        if(dspheat) then
  !
          !$acc parallel loop gang collapse(2) private(jj,k,i)
          do jj = 1, jlistnum
            do k = 1,km1
              !$acc loop vector
              do i = 1,myim(jj)
                diss(i,k,jj) = dku(i,k,jj)*shr2(i,k,jj)-g*ti(i,k,jj)*dkt(i,k,jj)*bf(i,k,jj)
    !           diss(i,k) = dku(i,k)*shr2(i,k)
              enddo
            enddo
          enddo
  !
  !     add dissipative heating at the first model layer
  !
  !>  Next, the temperature tendency is updated following equation 14.
          !$acc parallel loop gang private(jj,i,tem,tem1,tem2,ttend)
          do jj = 1, jlistnum
            !$acc loop vector
            do i = 1,myim(jj)
               tem   = govrth(i,jj)*sflux(i,jj)
               tem1  = tem + stress(i,jj)*spd1(i,jj)/zl(i,1,jj)
               tem2  = 0.5 * (tem1+diss(i,1,jj))
               tem2  = max(tem2, 0.)
               ttend = tem2 / cp
               tau(i,1,jj) = tau(i,1,jj)+0.5*ttend
            enddo
          enddo
  !
  !     add dissipative heating above the first model layer
  ! 
          !$acc parallel loop gang collapse(2) private(jj,k,i,tem,ttend)
          do jj = 1, jlistnum
            do k = 2,km1
              !$acc loop vector
              do i = 1,myim(jj)
                tem = 0.5 * (diss(i,k-1,jj)+diss(i,k,jj))
                tem  = max(tem, 0.)
                ttend = tem / cp
                !if (jj .eq. 5 .and. i .eq. 3 .and. k .eq. 2) write(*,*) 'a7',diss(i,k-1,jj),diss(i,k,jj)
                tau(i,k,jj) = tau(i,k,jj) + 0.5*ttend
              enddo
            enddo
          enddo
  !
        endif
  !    update temperature tendency to model layer mean temperature
      !$acc parallel loop gang collapse(2) private(jj,k,i)
      do jj = 1, jlistnum
        do  k = 1,km
           !$acc loop vector
           do i = 1,myim(jj)
              !if (jj .eq. 5 .and. i .eq. 3 .and. k .eq. 2) write(*,*) 'a8',t1(i,k,jj)
              t1(i,k,jj) = t1(i,k,jj) + tau(i,k,jj) * dt2
              !if (jj .eq. 5 .and. i .eq. 3 .and. k .eq. 2) write(*,*) 'a9',t1(i,k,jj),tau(i,k,jj)
           enddo
        enddo
      enddo
  !
  !     compute tridiagonal matrix elements for momentum
  !
  !>  ## Solve for the horizontal momentum tendencies and add them to the output tendency terms
  !!  As with the temperature and moisture tendencies, the horizontal momentum tendencies are calculated by solving tridiagonal matrices after the matrices are prepared in this section.
      !$acc parallel loop gang private(jj,i)
      do jj = 1, jlistnum
        !$acc loop vector
        do i=1,myim(jj)
           ad(i,1,jj) = 1.0 + beta(i,jj) * stress(i,jj) / spd1(i,jj)
           a1(i,1,jj) = u1(i,1,jj)
           a2(i,1,jj) = v1(i,1,jj)
        enddo
      enddo
  !
      !may cause some absolute error
      !$acc parallel loop gang private(jj,i,k,dtodsd,dtodsu,dsig,rdz,&
      !$acc&                   tem1,dsdz2,tem2,ptem,ptem1,ptem2)
      do jj = 1, jlistnum
        !$acc loop vector
        do i=1,myim(jj)
          !$acc loop seq
          do k = 1,km1
            dtodsd  = dt2/del(i,k,jj)
            dtodsu  = dt2/del(i,k+1,jj)
            dsig    = prsl(i,k,jj)-prsl(i,k+1,jj)
            rdz     = rdzt(i,k,jj)
            tem1    = dsig*dku(i,k,jj)*rdz
            dsdz2   = tem1 * rdz
            au(i,k,jj) = -dtodsd*dsdz2
            al(i,k,jj) = -dtodsu*dsdz2
  !
            if(pcnvflg(i,jj) .and. k < kpbl(i,jj)) then
               tem2      = dsig * rdz
               ptem      = 0.5 * tem2 * xmf(i,k,jj)
               ptem1     = dtodsd * ptem
               ptem2     = dtodsu * ptem
               ad(i,k,jj)   = ad(i,k,jj)-au(i,k,jj)-ptem1
               ad(i,k+1,jj) = 1.-al(i,k,jj)+ptem2
               au(i,k,jj)   = au(i,k,jj)-ptem1
               al(i,k,jj)   = al(i,k,jj)+ptem2
               ptem      = ucko(i,k,jj) + ucko(i,k+1,jj)
               a1(i,k,jj)   = a1(i,k,jj) - ptem1 * ptem
               a1(i,k+1,jj) = u1(i,k+1,jj) + ptem2 * ptem
               ptem      = vcko(i,k,jj) + vcko(i,k+1,jj)
               a2(i,k,jj)   = a2(i,k,jj) - ptem1 * ptem
               a2(i,k+1,jj) = v1(i,k+1,jj) + ptem2 * ptem
            else
               ad(i,k,jj)   = ad(i,k,jj)-au(i,k,jj)
               ad(i,k+1,jj) = 1.-al(i,k,jj)
               a1(i,k+1,jj) = u1(i,k+1,jj)
               a2(i,k+1,jj) = v1(i,k+1,jj)
            endif
  !
          enddo
        enddo
      enddo
  !
  !     solve tridiagonal problem for momentum
  !
      !$acc parallel loop gang
      do jj = 1, jlistnum
        call tridi2_gpu(ix,myim(jj),km,al(1,1,jj),ad(1,1,jj),&
                        au(1,1,jj),a1(1,1,jj),a2(1,1,jj),&
                        au(1,1,jj),a1(1,1,jj),a2(1,1,jj))
      enddo
        !do i=1,im(j1)
        !  fk      = 1./ad(i,1)
        !  au(i,1) = fk*au(i,1)
        !  a1(i,1) = fk*a1(i,1)
        !  a2(i,1) = fk*a2(i,1)
        !enddo
        !do k=2,km-1
        !  do i=1,im(j1)
        !    fk      = 1./(ad(i,k)-al(i,k-1)*au(i,k-1))
        !    au(i,k) = fk*au(i,k)
        !    a1(i,k) = fk*(a1(i,k)-al(i,k-1)*a1(i,k-1))
        !    a2(i,k) = fk*(a2(i,k)-al(i,k-1)*a2(i,k-1))
        !  enddo
        !enddo
        !do i=1,im(j1)
        !  fk      = 1./(ad(i,km)-al(i,km-1)*au(i,km-1))
        !  a1(i,km) = fk*(a1(i,km)-al(i,km-1)*a1(i,km-1))
        !  a2(i,km) = fk*(a2(i,km)-al(i,km-1)*a2(i,km-1))
        !enddo
        !do k=km-1,1,-1
        !  do i=1,im(j1)
        !    a1(i,k) = a1(i,k)-au(i,k)*a1(i,k+1)
        !    a2(i,k) = a2(i,k)-au(i,k)*a2(i,k+1)
        !  enddo
        !enddo


  !
  !     recover tendencies of momentum
  !
  !>  Finally, the tendencies are recovered from the tridiagonal solutions.
      !$acc parallel loop gang private(jj,i,k,utend,vtend)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1,myim(jj)
           !$acc loop seq
           do k = 1,km
              utend = (a1(i,k,jj)-u1(i,k,jj))*rdt
              vtend = (a2(i,k,jj)-v1(i,k,jj))*rdt
  !!            du(i,k)  = du(i,k)  + utend
  !!            dv(i,k)  = dv(i,k)  + vtend
              !du(i,k,jj)  = utend
              !dv(i,k,jj)  = vtend
              dusfc(i,jj) = dusfc(i,jj) + conw*del(i,k,jj)*utend
              dvsfc(i,jj) = dvsfc(i,jj) + conw*del(i,k,jj)*vtend
              u1(i,k,jj) = a1(i,k,jj)
              v1(i,k,jj) = a2(i,k,jj)
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
           enddo
        enddo
      enddo

  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
      !$acc parallel loop gang private(jj,i)
      do jj = 1, jlistnum
        !$acc loop vector
        do i = 1, myim(jj)
           hpbl(i,jj) = hpblx(i,jj)
           kpbl(i,jj) = kpblx(i,jj)
        enddo
      enddo
      !$acc exit data delete(a1,a2,ad,al,au,beta,bf,ckt,cku,diss, &
      !$acc&                 dusfc,dvsfc,dtsfc,dqsfc,dku,dkt,dktx, &
      !$acc&                 crb,flg,govrth,hgamt,hgamq,fk,fkk, &
      !$acc&                 hpblx,hrad,icld,jlist1,kcld,krad,kpblx, &
      !$acc&                 kx1,lcld,pblflg,pcnvflg,phih,phim, &
      !$acc&                 prinv,qcko,qlx,qtx,radmin,radx,rent,rbdn, &
      !$acc&                 rbup,rdzt,rtg,sfcflg,scuflg,shr2,sflux, &
      !$acc&                 tau,tcko,theta,thetae,thermal, &
      !$acc&                 thlvx1,thlvx,thvx,ti,ublflg,ucko,myim, &
      !$acc&                 ustar,ustmin,vcko,vrad,wscale,wscaleu,wstar, &
      !$acc&                 xmf,z0,zd,zdd,zi,zl,zol,xkzo,xkzmo,tx1,tx2)

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      return
      end
      
