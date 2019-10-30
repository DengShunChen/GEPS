      subroutine cuasc_n &
     &    (nxj,      klon,     klev,     klevp1,   klevm1,   ptenh,&
     &     pqenh,    pxenh,    puen,     pven,     pten,     pqen,&
     &     pqsen,    pgeo,     pgeoh,    pap,      paph,&
     &     pqte,     pverv,    klwmin,   ldcum,    phcbase,&
     &     ktype,    klab,     ptu,      pqu,      plu,&
     &     puu,      pvu,      pmfu,     pmfub,    pmfus,  &
     &     pmfuq,    pmful,    plude,    pdmfup,   kcbot,  &
     &     kctop,    kctop0,   kcum,     ztmst,    pqsenh, &
     &     plglac,   lndj,     wup,      wbase,    kdpl,   &
     &     pmfude_rate, zoentr1)
!     this routine does the calculations for cloud ascents
!     for cumulus parameterization
!     m.tiedtke         e.c.m.w.f.     7/86 modif.  12/89
!     y.wang            iprc           11/01 modif.
!     c.zhang           iprc           05/12 modif.
!***purpose.
!   --------
!          to produce cloud ascents for cu-parametrization
!          (vertical profiles of t,q,l,u and v and corresponding
!           fluxes as well as precipitation rates)
!***interface
!   ---------
!          this routine is called from *cumastr*.
!***method.
!  --------
!          lift surface air dry-adiabatically to cloud base
!          and then calculate moist ascent for
!          entraining/detraining plume.
!          entrainment and detrainment rates differ for
!          shallow and deep cumulus convection.
!          in case there is no penetrative or shallow convection
!          check for possibility of mid level convection
!          (cloud base values calculated in *cubasmc*)
!***externals
!   ---------
!          *cuadjtqn* adjust t and q due to condensation in ascent
!          *cuentrn*  calculate entrainment/detrainment rates
!          *cubasmcn* calculate cloud base values for midlevel convection
!***reference
!   ---------
!          (tiedtke,1989)
!***input variables:
!       ptenh [ztenh] - environ temperature on half levels. (cuini)
!       pqenh [zqenh] - env. specific humidity on half levels. (cuini)
!       puen - environment wind u-component. (mssflx)
!       pven - environment wind v-component. (mssflx)
!       pten - environment temperature. (mssflx)
!       pqen - environment specific humidity. (mssflx)
!       pqsen - environment saturation specific humidity. (mssflx)
!       pgeo - geopotential. (mssflx)
!       pgeoh [zgeoh] - geopotential on half levels, (mssflx)
!       pap - pressure in pa.  (mssflx)
!       paph - pressure of half levels. (mssflx)
!       pqte - moisture convergence (delta q/delta t). (mssflx)
!       pverv - large scale vertical velocity (omega). (mssflx)
!       klwmin [ilwmin] - level of minimum omega. (cuini)
!       klab [ilab] - level label - 1: sub-cloud layer.
!                                   2: condensation level (cloud base)
!       pmfub [zmfub] - updraft mass flux at cloud base. (cumastr)
!***variables modified by cuasc:
!       ldcum - logical denoting profiles. (cubase)
!       ktype - convection type - 1: penetrative  (cumastr)
!                                 2: stratocumulus (cumastr)
!                                 3: mid-level  (cuasc)
!       ptu - cloud temperature.
!       pqu - cloud specific humidity.
!       plu - cloud liquid water (moisture condensed out)
!       puu [zuu] - cloud momentum u-component.
!       pvu [zvu] - cloud momentum v-component.
!       pmfu - updraft mass flux.
!       pmfus [zmfus] - updraft flux of dry static energy. (cubasmc)
!       pmfuq [zmfuq] - updraft flux of specific humidity.
!       pmful [zmful] - updraft flux of cloud liquid water.
!       plude - liquid water returned to environment by detrainment.
!       pdmfup [zmfup] -
!       kcbot - cloud base level. (cubase)
!       kctop - cloud top level
!       kctop0 [ictop0] - estimate of cloud top. (cumastr)
!       kcum [icum] - flag to control the call
!-------------------------------------------------------------------
  USE mo_constants,     ONLY: g,       & ! gravity acceleration
                              cpd,     & ! specific heat at constant pressure
                              rd,      & ! gas constant for dry air
                              rv,      & !   idem      for water vapour
                              alv,     & ! latent heat for vaporisation
                              vtmpc1,  & ! vtmpc1=rv/rd-1
                              rcpd,    & ! rcpd=1./cpd
                              zrg,     & ! 1.0/g
                              ralfdcp, &
                              rtber,   &
                              rtice 
  USE mo_cumulus_flux,  ONLY: lmfdudv, & !
                              lmfmid,  & !
                              cmfcmin, & !
                              cprcon_n,  & !xb110
                              cmfctop

!-------------------------------------------------------------------
      implicit none

      integer  klev,klon,klevp1,klevm1,nxj
      real     ptenh(klon,klev),       pqenh(klon,klev), &
     &         pxenh(klon,klev),                       &
     &         puen(klon,klev),        pven(klon,klev),&
     &         pten(klon,klev),        pqen(klon,klev),&
     &         pgeo(klon,klev),        pgeoh(klon,klevp1),&
     &         pap(klon,klev),         paph(klon,klevp1),&
     &         pqsen(klon,klev),       pqte(klon,klev),&
     &         pverv(klon,klev),       pqsenh(klon,klev)
      real     ptu(klon,klev),         pqu(klon,klev),&
     &         puu(klon,klev),         pvu(klon,klev),&
     &         pmfu(klon,klev),        zph(klon),&
     &         pmfub(klon),            &
     &         pmfus(klon,klev),       pmfuq(klon,klev),&
     &         plu(klon,klev),         plude(klon,klev),&
     &         pmful(klon,klev),       pdmfup(klon,klev)
      real     zdmfen(klon),           zdmfde(klon),&
     &         zmfuu(klon),            zmfuv(klon),&
     &         zpbase(klon),           zqold(klon)              
      real     phcbase(klon),          zluold(klon)
      real     zprecip(klon),          zlrain(klon,klev)
      real     zbuo(klon,klev),        kup(klon,klev)
      real     wup(klon)
      real     wbase(klon),            zodetr(klon,klev)
      real     plglac(klon,klev)

      real     eta(klon),dz(klon)

      integer  klwmin(klon),           ktype(klon),&
     &         klab(klon,klev),        kcbot(klon),&
     &         kctop(klon),            kctop0(klon)
      integer  lndj(klon)
      logical  ldcum(klon),            loflag(klon)
      logical  llo2,llo3,              llo1(klon)

      integer  kdpl(klon)
      real     zoentr(klon),           zdpmean(klon)
      real     pdmfen(klon,klev),      pmfude_rate(klon,klev)
      real     ushear(klon,klev),      zoentr1(klon,klev)      
      integer  jl,jk
      integer  ikb,icum,itopm2,ik,icall,is,kcum,jlm,jll
      integer  jlx(klon)
      real     ztmst,zcons2,zfacbuo,zprcdgw,z_cwdrag,z_cldmax,z_cwifrac,z_cprc2
      real     zmftest,zmfmax,zqeen,zseen,zscde,zqude
      real     zmfusk,zmfuqk,zmfulk
      real     zbc,zbe,zkedke,zmfun,zwu,zprcon,zdt,zcbf,zzco
      real     zlcrit,zdfi,zc,zd,zint,zlnew,zvw,zvi,zalfaw,zrold
      real     zrnew,zz,zdmfeu,zdmfdu,dp
      real     zfac,zbuoc,zdkbuo,zdken,zvv,zarg,zchange,zxe,zxs,zdshrd
      real     atop1,atop2,abot
      real     foealfa,ushmax                 
!--------------------------------
!*    1.       specify parameters
!--------------------------------
      do jk = 1, klev
        do jl = 1,nxj
          zoentr1(jl,jk) = 0.
          ushear(jl,jk) = 0.
        enddo
      enddo
      zcons2=3./(g*ztmst)
      zfacbuo = 0.5/(1.+0.5)
      zprcdgw = cprcon_n*zrg
      z_cldmax = 5.e-3
      z_cwifrac = 0.5
      z_cprc2 = 0.5
      z_cwdrag = (3.0/8.0)*0.506/0.2
!---------------------------------
!     2.        set default values
!---------------------------------
!xb110>
      kctop = 0.
      jlx   = 0
      zph   = 0.
      zbuoc = 0.
!xb110<
      llo3 = .false.
       do jl = 1, nxj         
        zluold(jl)=0.
        wup(jl)=0.
        zdpmean(jl)=0.
        zoentr(jl)=0.
        if(.not.ldcum(jl)) then
          ktype(jl)=0
          kcbot(jl) = -1
          pmfub(jl) = 0.
          pqu(jl,klev) = 0.
        end if
      end do

! initialize variout quantities     
      do jk=1,klev
       do jl = 1, nxj         
          if(jk.ne.kcbot(jl)) plu(jl,jk)=0.
          pmfu(jl,jk)=0.
          pmfus(jl,jk)=0.
          pmfuq(jl,jk)=0.
          pmful(jl,jk)=0.
          plude(jl,jk)=0.
          plglac(jl,jk)=0.
          pdmfup(jl,jk)=0.
          zlrain(jl,jk)=0.
          zbuo(jl,jk)=0.
          kup(jl,jk)=0.
          pdmfen(jl,jk) = 0.
          pmfude_rate(jl,jk) = 0.
          if(.not.ldcum(jl).or.ktype(jl).eq.3) klab(jl,jk)=0
          if(.not.ldcum(jl).and.paph(jl,jk).lt.4.e4) kctop0(jl)=jk
      end do
      end do

       do jl = 1, nxj         
        if ( ktype(jl) == 3 ) ldcum(jl) = .false.
      end do
!------------------------------------------------
!     3.0      initialize values at cloud base level
!------------------------------------------------
       do jl = 1, nxj         
        kctop(jl)=kcbot(jl)
        if(ldcum(jl)) then
          ikb = kcbot(jl)
          kup(jl,ikb) = 0.5*wbase(jl)**2
          pmfu(jl,ikb) = pmfub(jl)
          pmfus(jl,ikb) = pmfub(jl)*(cpd*ptu(jl,ikb)+pgeoh(jl,ikb))
          pmfuq(jl,ikb) = pmfub(jl)*pqu(jl,ikb)
          pmful(jl,ikb) = pmfub(jl)*plu(jl,ikb)
        end if
      end do
!
!-----------------------------------------------------------------
!     4.       do ascent: subcloud layer (klab=1) ,clouds (klab=2)
!              by doing first dry-adiabatic ascent and then
!              by adjusting t,q and l accordingly in *cuadjtqn*,
!              then check for buoyancy and set flags accordingly
!-----------------------------------------------------------------
!
      do jk=klevm1,3,-1
!             specify cloud base values for midlevel convection
!             in *cubasmc* in case there is not already convection
! ---------------------------------------------------------------------
      ik=jk
      call cubasmc_n&
     &    (nxj,      klon,     klev,     klevm1,   ik,      pten,&
     &     pqen,     pqsen,    puen,     pven,    pverv,&
     &     pgeo,     pgeoh,    ldcum,   ktype,   klab,  zlrain,&
     &     pmfu,     pmfub,    kcbot,   ptu,&
     &     pqu,      plu,      puu,     pvu,      pmfus,&
     &     pmfuq,    pmful,    pdmfup)
      is = 0
      jlm = 0
       do jl = 1, nxj         
        loflag(jl) = .false.
        zprecip(jl) = 0.
        llo1(jl) = .false.
        is = is + klab(jl,jk+1)
        if ( klab(jl,jk+1) == 0 ) klab(jl,jk) = 0
        if ( (ldcum(jl) .and. klab(jl,jk+1) == 2) .or.   &
             (ktype(jl) == 3 .and. klab(jl,jk+1) == 1) ) then
          loflag(jl) = .true.
          jlm = jlm + 1
          jlx(jlm) = jl
        end if
        zph(jl) = paph(jl,jk)
        if ( ktype(jl) == 3 .and. jk == kcbot(jl) ) then
          zmfmax = (paph(jl,jk)-paph(jl,jk-1))*zcons2
          if ( pmfub(jl) > zmfmax ) then
            zfac = zmfmax/pmfub(jl)
            pmfu(jl,jk+1) = pmfu(jl,jk+1)*zfac
            pmfus(jl,jk+1) = pmfus(jl,jk+1)*zfac
            pmfuq(jl,jk+1) = pmfuq(jl,jk+1)*zfac
            pmfub(jl) = zmfmax
          end if
          pmfub(jl)=min(pmfub(jl),zmfmax)
        end if
      end do

      if(is.gt.0) llo3 = .true.
!
!*     specify entrainment rates in *cuentr*
! -------------------------------------
      ik=jk
      call cuentr_n(nxj,klon,klev,ik,kcbot,ldcum,llo3, &
                  pgeoh,pmfu,zdmfen,zdmfde)
!
!      do adiabatic ascent for entraining/detraining plume
      if(llo3) then
! -------------------------------------------------------
!
       do jl = 1, nxj         
          zqold(jl) = 0.
        end do
        do jll = 1 , jlm
          jl = jlx(jll)
          zdmfde(jl) = min(zdmfde(jl),0.75*pmfu(jl,jk+1))
          if ( jk == kcbot(jl) ) then
            zoentr(jl) = -1.75e-3*(min(1.,pqen(jl,jk)/pqsen(jl,jk)) - &
                         1.)*(pgeoh(jl,jk)-pgeoh(jl,jk+1))*zrg
            zoentr(jl) = min(0.4,zoentr(jl))*pmfu(jl,jk+1)
          end if
          if ( jk < kcbot(jl) ) then
            zmfmax = (paph(jl,jk)-paph(jl,jk-1))*zcons2
            zxs = max(pmfu(jl,jk+1)-zmfmax,0.)
            wup(jl) = wup(jl) + kup(jl,jk+1)*(pap(jl,jk+1)-pap(jl,jk))
            zdpmean(jl) = zdpmean(jl) + pap(jl,jk+1) - pap(jl,jk)
            zdmfen(jl) = zoentr(jl)
            if ( ktype(jl) >= 2 ) then
              zdmfen(jl) = 2.0*zdmfen(jl)
              zdmfde(jl) = zdmfen(jl)
            end if
            zdmfde(jl) = zdmfde(jl) * &
                         (1.6-min(1.,pqen(jl,jk)/pqsen(jl,jk)))
            zmftest = pmfu(jl,jk+1) + zdmfen(jl) - zdmfde(jl)
            zchange = max(zmftest-zmfmax,0.)
            zxe = max(zchange-zxs,0.)
            zdmfen(jl) = zdmfen(jl) - zxe
            zchange = zchange - zxe
            zdmfde(jl) = zdmfde(jl) + zchange
          end if
          pdmfen(jl,jk) = zdmfen(jl) - zdmfde(jl)
          pmfu(jl,jk) = pmfu(jl,jk+1) + zdmfen(jl) - zdmfde(jl)
          zqeen = pqenh(jl,jk+1)*zdmfen(jl)
          zseen = (cpd*ptenh(jl,jk+1)+pgeoh(jl,jk+1))*zdmfen(jl)
          zscde = (cpd*ptu(jl,jk+1)+pgeoh(jl,jk+1))*zdmfde(jl)
          zqude = pqu(jl,jk+1)*zdmfde(jl)
          plude(jl,jk) = plu(jl,jk+1)*zdmfde(jl)
          zmfusk = pmfus(jl,jk+1) + zseen - zscde
          zmfuqk = pmfuq(jl,jk+1) + zqeen - zqude
          zmfulk = pmful(jl,jk+1) - plude(jl,jk)
          plu(jl,jk) = zmfulk*(1./max(cmfcmin,pmfu(jl,jk)))
          pqu(jl,jk) = zmfuqk*(1./max(cmfcmin,pmfu(jl,jk)))
          ptu(jl,jk) = (zmfusk * &
            (1./max(cmfcmin,pmfu(jl,jk)))-pgeoh(jl,jk))*rcpd
          ptu(jl,jk) = max(100.,ptu(jl,jk))
          ptu(jl,jk) = min(400.,ptu(jl,jk))
          zqold(jl) = pqu(jl,jk)
          zlrain(jl,jk) = zlrain(jl,jk+1)*(pmfu(jl,jk+1)-zdmfde(jl)) * &
                          (1./max(cmfcmin,pmfu(jl,jk)))
          zluold(jl) = plu(jl,jk)
        end do
! reset to environmental values if below departure level
       do jl = 1, nxj         
          if ( jk > kdpl(jl) ) then
            ptu(jl,jk) = ptenh(jl,jk)
            pqu(jl,jk) = pqenh(jl,jk)
            plu(jl,jk) = 0.
            zluold(jl) = plu(jl,jk)
          end if
        end do
!*             do corrections for moist ascent
!*             by adjusting t,q and l in *cuadjtq*
!------------------------------------------------
      ik=jk
      icall=1
!
      if ( jlm > 0 ) then
        call cuadjtq_n(nxj,klon,klev,ik,zph,ptu,pqu,loflag,icall)
      end if
! compute the upfraft speed in cloud layer
        do jll = 1 , jlm
          jl = jlx(jll)
          if ( pqu(jl,jk) /= zqold(jl) ) then
            plglac(jl,jk) = plu(jl,jk) * &
                           ((1.-foealfa(ptu(jl,jk)))- &
                            (1.-foealfa(ptu(jl,jk+1))))
            ptu(jl,jk) = ptu(jl,jk) + ralfdcp*plglac(jl,jk)
          end if
        end do
        do jll = 1 , jlm
          jl = jlx(jll)
          if ( pqu(jl,jk) /= zqold(jl) ) then
            klab(jl,jk) = 2
            plu(jl,jk) = plu(jl,jk) + zqold(jl) - pqu(jl,jk)
            zbc = ptu(jl,jk)*(1.+vtmpc1*pqu(jl,jk)-plu(jl,jk+1) - &
              zlrain(jl,jk+1))
            zbe = ptenh(jl,jk)*(1.+vtmpc1*pqenh(jl,jk))
            zbuo(jl,jk) = zbc - zbe
! set flags for the case of midlevel convection
            if ( ktype(jl) == 3 .and. klab(jl,jk+1) == 1 ) then
              if ( zbuo(jl,jk) > -0.5 ) then
                ldcum(jl) = .true.
                kctop(jl) = jk
                kup(jl,jk) = 0.5
              else
                klab(jl,jk) = 0
                pmfu(jl,jk) = 0.
                plude(jl,jk) = 0.
                plu(jl,jk) = 0.
              end if
            end if
            if ( klab(jl,jk+1) == 2 ) then
              if ( zbuo(jl,jk) < 0. ) then
                ptenh(jl,jk) = 0.5*(pten(jl,jk)+pten(jl,jk-1))
                pqenh(jl,jk) = 0.5*(pqen(jl,jk)+pqen(jl,jk-1))
                zbuo(jl,jk) = zbc - ptenh(jl,jk)*(1.+vtmpc1*pqenh(jl,jk))
              end if
              zbuoc = (zbuo(jl,jk) / &
                (ptenh(jl,jk)*(1.+vtmpc1*pqenh(jl,jk)))+zbuo(jl,jk+1) / &
                (ptenh(jl,jk+1)*(1.+vtmpc1*pqenh(jl,jk+1))))*0.5
              zdkbuo = (pgeoh(jl,jk)-pgeoh(jl,jk+1))*zfacbuo*zbuoc
! mixing and "pressure" gradient term in upper troposphere
              if ( zdmfen(jl) > 0. ) then
                zdken = min(1.,(1.+z_cwdrag)*zdmfen(jl) / &
                        max(cmfcmin,pmfu(jl,jk+1)))
              else
                zdken = min(1.,(1.+z_cwdrag)*zdmfde(jl) / &
                        max(cmfcmin,pmfu(jl,jk+1)))
              end if
              kup(jl,jk) = (kup(jl,jk+1)*(1.-zdken)+zdkbuo) / &
                           (1.+zdken)
              if ( zbuo(jl,jk) < 0. ) then
                zkedke = kup(jl,jk)/max(1.e-10,kup(jl,jk+1))
                zkedke = max(0.,min(1.,zkedke))
                zmfun = sqrt(zkedke)*pmfu(jl,jk+1)
                zdmfde(jl) = max(zdmfde(jl),pmfu(jl,jk+1)-zmfun)
                plude(jl,jk) = plu(jl,jk+1)*zdmfde(jl)
                pmfu(jl,jk) = pmfu(jl,jk+1) + zdmfen(jl) - zdmfde(jl)
              end if
              if ( zbuo(jl,jk) > 0.  ) then
                ikb = kcbot(jl)
                zoentr(jl) = 1.75e-3*(0.3-(min(1.,pqen(jl,jk-1) /    &
                  pqsen(jl,jk-1))-1.))*(pgeoh(jl,jk-1)-pgeoh(jl,jk)) * &
                  zrg*min(1.,pqsen(jl,jk)/pqsen(jl,ikb))**3
                zoentr1(jl,jk)=zoentr(jl)
                zoentr(jl) = min(0.4,zoentr(jl))*pmfu(jl,jk)
              else
                zoentr(jl) = 0.
              end if
! erase values if below departure level
              if ( jk > kdpl(jl) ) then
                pmfu(jl,jk) = pmfu(jl,jk+1)
                kup(jl,jk) = 0.5
              end if
              if ( kup(jl,jk) > 0. .and. pmfu(jl,jk) > 0. ) then
                kctop(jl) = jk
                llo1(jl) = .true.
              else
                klab(jl,jk) = 0
                pmfu(jl,jk) = 0.
                kup(jl,jk) = 0.
                zdmfde(jl) = pmfu(jl,jk+1)
                plude(jl,jk) = plu(jl,jk+1)*zdmfde(jl)
              end if
! save detrainment rates for updraught
              if ( pmfu(jl,jk+1) > 0. ) pmfude_rate(jl,jk) = zdmfde(jl)
            end if
          else if ( ktype(jl) == 2 .and. pqu(jl,jk) == zqold(jl) ) then
            klab(jl,jk) = 0
            pmfu(jl,jk) = 0.
            kup(jl,jk) = 0.
            zdmfde(jl) = pmfu(jl,jk+1)
            plude(jl,jk) = plu(jl,jk+1)*zdmfde(jl)
            pmfude_rate(jl,jk) = zdmfde(jl)
          end if
        end do

       do jl = 1, nxj         
          if ( llo1(jl) ) then
! conversions only proceeds if plu is greater than a threshold liquid water
! content of 0.3 g/kg over water and 0.5 g/kg over land to prevent precipitation
! generation from small water contents.
            if ( lndj(jl).eq.1 ) then
              zdshrd = 5.e-4
            else
              zdshrd = 3.e-4
            end if
            ikb=kcbot(jl)
            if ( plu(jl,jk) > zdshrd )then
              zwu = min(15.0,sqrt(2.*max(0.1,kup(jl,jk+1))))
              zprcon = zprcdgw/(0.75*zwu)
! PARAMETERS FOR BERGERON-FINDEISEN PROCESS (T < -5C)
              zdt = min(rtber-rtice,max(rtber-ptu(jl,jk),0.))
              zcbf = 1. + z_cprc2*sqrt(zdt)
              zzco = zprcon*zcbf
              zlcrit = zdshrd/zcbf
              zdfi = pgeoh(jl,jk) - pgeoh(jl,jk+1)
              zc = (plu(jl,jk)-zluold(jl))
              zarg = (plu(jl,jk)/zlcrit)**2
              if ( zarg < 25.0 ) then
                zd = zzco*(1.-exp(-zarg))*zdfi
              else
                zd = zzco*zdfi
              end if
              zint = exp(-zd)
              zlnew = zluold(jl)*zint + zc/zd*(1.-zint)
              zlnew = max(0.,min(plu(jl,jk),zlnew))
              zlnew = min(z_cldmax,zlnew)
              zprecip(jl) = max(0.,zluold(jl)+zc-zlnew)
              pdmfup(jl,jk) = zprecip(jl)*pmfu(jl,jk)
              zlrain(jl,jk) = zlrain(jl,jk) + zprecip(jl)
              plu(jl,jk) = zlnew
            end if
          end if
        end do
       do jl = 1, nxj         
          if ( llo1(jl) ) then
            if ( zlrain(jl,jk) > 0. ) then
              zvw = 21.18*zlrain(jl,jk)**0.2
              zvi = z_cwifrac*zvw
              zalfaw = foealfa(ptu(jl,jk))
              zvv = zalfaw*zvw + (1.-zalfaw)*zvi
              zrold = zlrain(jl,jk) - zprecip(jl)
              zc = zprecip(jl)
              zwu = min(15.0,sqrt(2.*max(0.1,kup(jl,jk))))
              zd = zvv/zwu
              zint = exp(-zd)
              zrnew = zrold*zint + zc/zd*(1.-zint)
              zrnew = max(0.,min(zlrain(jl,jk),zrnew))
              zlrain(jl,jk) = zrnew
            end if
          end if
        end do
        do jll = 1 , jlm
          jl = jlx(jll)
          pmful(jl,jk) = plu(jl,jk)*pmfu(jl,jk)
          pmfus(jl,jk) = (cpd*ptu(jl,jk)+pgeoh(jl,jk))*pmfu(jl,jk)
          pmfuq(jl,jk) = pqu(jl,jk)*pmfu(jl,jk)
        end do
      end if
      end do
!----------------------------------------------------------------------
! 5.       final calculations
! ------------------
       do jl = 1, nxj         

       if ( kctop(jl) == -1 ) ldcum(jl) = .false.
        kcbot(jl) = max(kcbot(jl),kctop(jl))
!xb110> 20190711
       if (.not.ldcum(jl)) then
          kcbot(jl) = -1
          kctop(jl) = -1
       end if
!xb110< 20190711

        if ( ldcum(jl) ) then
          wup(jl) = max(1.e-2,wup(jl)/max(1.,zdpmean(jl)))
          wup(jl) = sqrt(2.*wup(jl))
        end if
      end do

      return
      end subroutine cuasc_n
