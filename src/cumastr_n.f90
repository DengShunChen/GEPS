!***********************************************************
!           subroutine cumastrn
!***********************************************************
SUBROUTINE cumastr_n  &
     &      (nxj   ,klon ,klev ,klevp1,klevm1 ,&
     &       pten  ,pqen ,pxen ,puen  ,pven   ,&
     &       ldland,pverv,pqsen,pqhfl ,pap    ,&
     &       paph  ,pgeo ,ptte ,pqte  ,pxtec  ,&
     &       pvom  ,pvol ,prsfc,pssfc ,kcbot  ,&
     &       kctop ,ztmst,jin  ,ptu   ,pqu    ,&
     &       pmfu  ,pmfd ,prain,pcte  ,phhfl  ,&
     &       lndj  ,ldcum,dx   ,kcnv   ,&
     &       flash)
!
!***cumastrn*  master routine for cumulus massflux-scheme
!     m.tiedtke      e.c.m.w.f.      1986/1987/1989
!     modifications
!     y.wang         i.p.r.c         2001
!     c.zhang                        2012
!***purpose
!   -------
!          this routine computes the physical tendencies of the
!     prognostic variables t,q,u and v due to convective processes.
!     processes considered are: convective fluxes, formation of
!     precipitation, evaporation of falling rain below cloud base,
!     saturated cumulus downdrafts.
!***method
!   ------
!     parameterization is done using a massflux-scheme.
!        (1) define constants and parameters
!        (2) specify values (t,q,qs...) at half levels and
!            initialize updraft- and downdraft-values in 'cuinin'
!        (3) calculate cloud base in 'cutypen', calculate cloud types in cutypen,
!            and specify cloud base massflux
!        (4) do cloud ascent in 'cuascn' in absence of downdrafts
!        (5) do downdraft calculations:
!              (a) determine values at lfs in 'cudlfsn'
!              (b) determine moist descent in 'cuddrafn'
!              (c) recalculate cloud base massflux considering the
!                  effect of cu-downdrafts
!        (6) do final adjusments to convective fluxes in 'cuflxn',
!            do evaporation in subcloud layer
!        (7) calculate increments of t and q in 'cudtdqn'
!        (8) calculate increments of u and v in 'cududvn'
!***externals.
!   ----------
!       cuinin:  initializes values at vertical grid used in cu-parametr.
!       cutypen: cloud bypes, 1: deep cumulus 2: shallow cumulus
!       cuascn:  cloud ascent for entraining plume
!       cudlfsn: determines values at lfs for downdrafts
!       cuddrafn:does moist descent for cumulus downdrafts
!       cuflxn:  final adjustments to convective fluxes (also in pbl)
!       cudqdtn: updates tendencies for t and q
!       cududvn: updates tendencies for u and v
!***switches.
!   --------
!          lmfmid=.t.   midlevel convection is switched on
!          lmfdd=.t.    cumulus downdrafts switched on
!          lmfdudv=.t.  cumulus friction switched on
!***
!     model parameters (defined in subroutine cuparam)
!     ------------------------------------------------
!     entrdd     entrainment rate for cumulus downdrafts
!     cmfcmax    maximum massflux value allowed for
!     cmfcmin    minimum massflux value (for safety)
!     cmfdeps    fractional massflux for downdrafts at lfs
!     cprcon     coefficient for conversion from cloud water to rain
!***reference.
!   ----------
!         1. paper on massflux scheme (Tiedtke,1989)
!         2. Bechtold,P.et.al.(2014):Representing equilibrium and
!            nonequilibrium convection in large-scale models.
!         3. ECMWF,23-25 May 2016 NWP Training Course: Parametrization
!            of subgrid physical processes.
!            link:
!     https://software.ecmwf.int/wiki/display/OPTR/NWP+Training+Material+2016
!-----------------------------------------------------------------
USE mo_constants,     only: rd,      &! gas constant for dry air
                              c4les,   &!
                              c5les,   &!
                              cpd,     &! specific heat at constant pressure
                              g,       &! gravity acceleration
                              vtmpc1,  &! vtmpc1=rv/rd-1
                              alv,     &! latent heat for vaporisation
                              zrg,     &! 1.0/g
                              rcpd
USE mo_cumulus_flux,  only: lmfdudv, &! true if cum. friction is switched on
                              cmfdeps, &! fractional massflux for downdr. at lfs
                              lmfdd,   &! true if cum. downdraft is switched on
                              entrpen, &! entrainment rate for penetrative conv.
                              entrscv, &! entrainment rate for shallow conv.
                              zdnoprc, &
                              momtrans,&
                              lmfpen,  &
                              cmfcmin, &
                              lmfpen,  &
                              lmfscv

      use mpe
      use rank
      use index                              

      implicit none
      integer klev,klon,klevp1,klevm1,nxj,jin
      logical ldland(klon)
      real     pten(klon,klev),        pqen(klon,klev),& 
     &         puen(klon,klev),        pven(klon,klev),&
     &         ptte(klon,klev),        pqte(klon,klev),&
     &         pvom(klon,klev),        pvol(klon,klev),&
     &         pqsen(klon,klev),       pgeo(klon,klev),&
     &         pap(klon,klev),         paph(klon,klevp1),&
     &         pverv(klon,klev),       pqhfl(klon),&
     &         phhfl(klon)
      real     pxtec(klon,klev),       pxen(klon,klev)            
      real     ptu(klon,klev),         pqu(klon,klev),&
     &         plu(klon,klev),         plude(klon,klev),&
     &         pmfu(klon,klev),        pmfd(klon,klev),&
     &         prain(klon),&
     &         prsfc(klon),            pssfc(klon)
      real     ztenh(klon,klev),       zqenh(klon,klev),&
     &         zgeoh(klon,klevp1),     zqsenh(klon,klev),&
     &         ztd(klon,klev),         zqd(klon,klev),&
     &         zmfus(klon,klev),       zmfds(klon,klev),&
     &         zmfuq(klon,klev),       zmfdq(klon,klev),&
     &         zdmfup(klon,klev),      zdmfdp(klon,klev),&
     &         zmful(klon,klev),       zrfl(klon),&
     &         zuu(klon,klev),         zvu(klon,klev),&
     &         zud(klon,klev),         zvd(klon,klev),&
     &         zlglac(klon,klev),                     &
               zxenh(klon,klev)
      real     pmflxr(klon,klevp1),    pmflxs(klon,klevp1)
      real     zhcbase(klon),& 
     &         zmfub(klon),            zmfub1(klon),&
     &         zdhpbl(klon)          
      real     zsfl(klon),             zdpmel(klon,klev),&
     &         pcte(klon,klev),        zcape(klon),&
     &         zcape1(klon),           zcape2(klon),&
     &         ztauc(klon),            ztaubl(klon),&
     &         zheat(klon)            
      real     wup(klon),              zdqcv(klon)            
      real     wbase(klon),            zmfuub(klon)
      real     upbl(klon)
      real     xmin,xmax
      real     pmfude_rate(klon,klev), pmfdde_rate(klon,klev)
      real     zmfuus(klon,klev),      zmfdus(klon,klev)
      real     zuv2(klon,klev),ztenu(klon,klev),ztenv(klon,klev)
      real     zmfuvb(klon),zsum12(klon),zsum22(klon)
      integer  ilab(klon,klev),        idtop(klon),&
     &         ictop0(klon),           ilwmin(klon)
      integer  kdpl(klon)
      integer  kcbot(klon),            kctop(klon),&
     &         ktype(klon),            lndj(klon)
      logical  ldcum(klon)
      logical  loddraf(klon),          llo1,   llo2(klon)
      integer  p950,p650
      

!  local varaiables
      real     zcons,zcons2,zqumqe,zdqmin,zdh,zmfmax,xlat
      real     zalfaw,zalv,zqalv,zc5ldcp,zc4les,zhsat,zgam,zzz,zhhat
      real     zpbmpt,zro,zdz,zdp,zeps,zfac,wspeed
      integer  jl,jk,ik
      integer  ikb,ikt,icum,itopm2,jim
      real     ztmst,ztau,zerate,zderate,zmfa
      real     zmfs(klon),pmean(klev),zlon,pi,rad,torad,dlon,dlat
      real     apha,capa,klvl
      real     zduten,zdvten,ztdis,pgf_u,pgf_v
      real     zoentr1(klon,klev),tmp2(klon,klev),tmp3(klon,klev) &
     &        ,tmp4(klon,klev),tmp5(klon,klev)
!xb110>
      real     mdlon,gdx,re,rr
      integer  kcnv(klon)
      real     sumpap(klon)
      logical  adj
!for lightning parameterization
      real     pluu(klon,klev),pf(klon,klev),rho(klon,klev),papn(klon,klev)
      real     flash(klon),dx(klon)
!xb110<
!-------------------------------------------
!     1.    specify constants and parameters
!-------------------------------------------
!xb110>
     pmean = 0.
     zlon  = 0.
     zrfl  = 0.
     sumpap = 0.
!xb110<
      zcons=1./(g*ztmst)
      zcons2=3./(g*ztmst)

      zlon = real(klon)
      do jk = klev , 1 , -1
!xb110> bug fixed, pap has an undefined value
!        pmean(jk) = sum(pap(:,jk))/zlon
      do jl = 1,nxj                      
        sumpap(jk) = sumpap(jk)+pap(jl,jk)
      end do
!xb110<
        pmean(jk) = sumpap(jk)/zlon
      end do
      p950 = klev-2
      p650 = klev-2  
      do jk = klev , 3 , -1
        if ( pmean(jk)/pmean(klev)*1.013250e5 >  950.e2 ) p950 = jk
        if ( pmean(jk)/pmean(klev)*1.013250e5 >  650.e2 ) p650 = jk
      end do
      p950 = min(klev-2,p950)
!--------------------------------------------------------------
!*    2.    initialize values at vertical grid points in 'cuini'
!--------------------------------------------------------------
      call cuini_n &
     &    (nxj,      klon,     klev,     klevp1,   klevm1,&
           pten,     pqen,     pqsen,    pxen,     puen,     pven,  &
           pverv,    pgeo,     paph,     zgeoh,    ztenh, &
           zqenh,    zqsenh,   zxenh,    ilwmin,   ptu,      pqu,   &
           ztd,      zqd,      zuu,      zvu,      zud,   &
           zvd,      pmfu,     pmfd,     zmfus,    zmfds, &
           zmfuq,    zmfdq,    zdmfup,   zdmfdp,   zdpmel,&
           plu,      plude,    ilab,     jin)
!----------------------------------
!*    3.0   cloud base calculations
!----------------------------------
!*         (a) determine cloud base values in 'cutypen',
!              and the cumulus type 1 or 2
!          -------------------------------------------
       call cutype_n &
     &     (  nxj,     klon,     klev,     klevp1,   klevm1,     pqen,&
     &      ztenh,    zqenh,     zqsenh,    zgeoh,     paph,&
     &      phhfl,    pqhfl,       pgeo,    pqsen,      pap,&
     &       pten,     lndj,        ptu,      pqu,     ilab,&
     &      ldcum,    kcbot,     ictop0,    ktype,    wbase,&
              plu,     kdpl,     p650,      jin)

!*         (b) assign the first guess mass flux at cloud base
!              ------------------------------------------
       do jl = 1, nxj                  
         zdhpbl(jl)=0.0
         upbl(jl) = 0.0
         idtop(jl)=0
         zmfub(jl) = 0.
       end do

       do jk=2,klev
       do jl = 1, nxj                   
         if(jk.ge.kcbot(jl) .and. ldcum(jl)) then
            zdhpbl(jl)=zdhpbl(jl)+(alv*pqte(jl,jk)+cpd*ptte(jl,jk))&
     &                 *(paph(jl,jk+1)-paph(jl,jk))
            wspeed = sqrt(puen(jl,jk)**2 + pven(jl,jk)**2)
            upbl(jl) = upbl(jl) + wspeed*(paph(jl,jk+1)-paph(jl,jk))
         end if
       end do
       end do

        do jl = 1, nxj                
        if(ldcum(jl)) then
           ikb=kcbot(jl)
           zmfmax = (paph(jl,ikb)-paph(jl,ikb-1))*zcons2
           if(ktype(jl) == 1) then
             zmfub(jl)= 0.1*zmfmax
           else if ( ktype(jl) == 2 ) then
             zqumqe = pqu(jl,ikb) + plu(jl,ikb) - zqenh(jl,ikb)
             zdqmin = max(0.01*zqenh(jl,ikb),1.e-10)
             zdh = cpd*(ptu(jl,ikb)-ztenh(jl,ikb)) + alv*zqumqe
             zdh = g*max(zdh,1.e5*zdqmin)
             if ( zdhpbl(jl) > 0. ) then
               zmfub(jl) = zdhpbl(jl)/zdh
               zmfub(jl) = min(zmfub(jl),zmfmax)
             else
               zmfub(jl) = 0.1*zmfmax
               ldcum(jl) = .false.
             end if
            end if  
        else
           zmfub(jl) = 0.
        end if
      end do
!------------------------------------------------------
!*    4.0   determine cloud ascent for entraining plume
!------------------------------------------------------
!*    (a) do ascent in 'cuasc'in absence of downdrafts
!----------------------------------------------------------
      call cuasc_n &
     &    (nxj,      klon,     klev,     klevp1,   klevm1,   ztenh,  &
     &     zqenh,    zxenh,    puen,     pven,     pten,     pqen,   &
     &     pqsen,    pgeo,     zgeoh,    pap,      paph,   &
     &     pqte,     pverv,    ilwmin,   ldcum,    zhcbase,&
     &     ktype,    ilab,     ptu,      pqu,      plu,    &
     &     zuu,      zvu,      pmfu,     zmfub,    zmfus,  &
     &     zmfuq,    zmful,    plude,    zdmfup,   kcbot,  &
     &     kctop,    ictop0,   icum,     ztmst,    zqsenh, &
     &     zlglac,   lndj,     wup,      wbase,    kdpl,   &
     &     pmfude_rate,zoentr1 )

!xb110>
!      !storing the cloud water (kg/kg) value within the convective updraft
!       do jk = 1,klev
!        do jl = 1,nxj
!          pluu(jl,jk) = plu(jl,jk)          ! cloud liquid water within updraft
!        enddo
!       enddo
!xb110<

!*     (b) check cloud depth and change entrainment rate accordingly
!          calculate precipitation rate (for downdraft calculation)
!------------------------------------------------------------------
        do jl = 1, nxj               
        if ( ldcum(jl) ) then
          ikb = kcbot(jl)
          itopm2 = kctop(jl)
          zpbmpt = paph(jl,ikb) - paph(jl,itopm2)
          if ( ktype(jl) == 1 .and. zpbmpt <  zdnoprc ) ktype(jl) = 2
          if ( ktype(jl) == 2 .and. zpbmpt >= zdnoprc ) ktype(jl) = 1
          ictop0(jl) = kctop(jl)
        end if
        zrfl(jl)=zdmfup(jl,1)
      end do

      do jk=2,klev
        do jl = 1, nxj              
          zrfl(jl)=zrfl(jl)+zdmfup(jl,jk)
        end do
      end do

      do jk = 1,klev
       do jl = 1, nxj        
        pmfd(jl,jk) = 0.
        zmfds(jl,jk) = 0.
        zmfdq(jl,jk) = 0.
        zdmfdp(jl,jk) = 0.
        zdpmel(jl,jk) = 0.
      end do
      end do
!-----------------------------------------
!*    6.0   cumulus downdraft calculations
!-----------------------------------------
      if(lmfdd) then
!*      (a) determine lfs in 'cudlfsn'
!--------------------------------------
        call cudlfs_n    &
     &    (nxj,      klon,     klev,                        &
     &     kcbot,    kctop,    lndj,   ldcum,               &
     &     ztenh,    zqenh,    puen,     pven,              &
     &     pten,     pqsen,    pgeo,                        &
     &     zgeoh,    paph,     ptu,      pqu,      plu,     &
     &     zuu,      zvu,      zmfub,    zrfl,              &
     &     ztd,      zqd,      zud,      zvd,               &
     &     pmfd,     zmfds,    zmfdq,    zdmfdp,            &
     &     idtop,    loddraf)
!*     (b)  determine downdraft t,q and fluxes in 'cuddrafn'
!------------------------------------------------------------
        call cuddraf_n                                         &
     &   ( nxj,       klon,     klev,     loddraf,             &
     &     ztenh,    zqenh,    puen,     pven,                 &
     &     pgeo,     zgeoh,    paph,     zrfl,                  &
     &     ztd,      zqd,      zud,      zvd,      pmfu,         &
     &     pmfd,     zmfds,    zmfdq,    zdmfdp,   pmfdde_rate )
!-----------------------------------------------------------
      end if
!
!-----------------------------------------------------------------------
!* 6.0          closure and clean work
! ------
!-- 6.1 recalculate cloud base massflux from a cape closure
!       for deep convection (ktype=1) 
!
!xb110>
        zheat=0.0
        zcape=0.0
        zcape1=0.0
        zcape2=0.0
        zmfub1=0.0
        ztauc=0.0
        ztaubl=0.0
        zmfs=0.0
!xb110<
      do jl = 1, nxj
      if(ldcum(jl) .and. ktype(jl) .eq. 1) then
        ikb = kcbot(jl)
        ikt = kctop(jl)
!xb110>
!        zheat(jl)=0.0
!        zcape(jl)=0.0
!        zcape1(jl)=0.0
!        zcape2(jl)=0.0
!xb110<
        zmfub1(jl)=zmfub(jl)
        upbl(jl) = max(2.,upbl(jl)/(paph(jl,klev+1)-paph(jl,ikb)))
        ztauc(jl)  = (zgeoh(jl,ikt)-zgeoh(jl,ikb)) / &
                   ((2.+ min(15.0,wup(jl)))*g)
!> xb110
!        if(lndj(jl) .eq. 1) then 
!          ztaubl(jl) = ztauc(jl)
!        else
!          ztaubl(jl) = (zgeoh(jl,ikb)-zgeoh(jl,klev+1))*zrg/upbl(jl)
!        end if
          ztaubl(jl) = ztauc(jl)
!< xb110
      end if    
      end do
!

      do jk = 1 , klev
       do jl = 1, nxj            
        llo1 = ldcum(jl) .and. ktype(jl) .eq. 1
        if ( llo1 .and. jk <= kcbot(jl) .and. jk > kctop(jl) ) then
          ikb = kcbot(jl)
          zdz = pgeo(jl,jk-1)-pgeo(jl,jk)
          zdp = pap(jl,jk)-pap(jl,jk-1)
          zheat(jl) = zheat(jl) + ((pten(jl,jk-1)-pten(jl,jk)+zdz*rcpd) / &
                      ztenh(jl,jk)+vtmpc1*(pqen(jl,jk-1)-pqen(jl,jk))) * &
                      (g*(pmfu(jl,jk)+pmfd(jl,jk)))
          zcape1(jl) = zcape1(jl) + ((ptu(jl,jk)-ztenh(jl,jk))/ztenh(jl,jk) + &
                      vtmpc1*(pqu(jl,jk)-zqenh(jl,jk))-plu(jl,jk))*zdp
        end if

        if ( llo1 .and. jk >= kcbot(jl) ) then
          ikb = kcbot(jl)
          if(paph(jl,klev+1)-paph(jl,ikb) <= 50.e2) then
            zdp = paph(jl,jk+1)-paph(jl,jk)
            zcape2(jl) = zcape2(jl) + ztaubl(jl)* &
                      (ptte(jl,jk)+vtmpc1*pten(jl,jk)*pqte(jl,jk))*zdp 
          end if
        end if
      end do
      end do

       do jl = 1, nxj            
       if(ldcum(jl).and.ktype(jl).eq.1) then
           ikb = kcbot(jl)
           ikt = kctop(jl)
           ztau = ztauc(jl) * (1.+1.33e-5*dx(jl))       
           ztau = max(ztmst,ztau)
           ztau = max(720.,ztau)
           ztau = min(10800.,ztau)
           zcape2(jl)= max(0.,zcape2(jl))
           zcape(jl) = max(0.,min(zcape1(jl)-zcape2(jl),5000.))
           zheat(jl) = max(1.e-4,zheat(jl))
           zmfub1(jl) = (zcape(jl)*zmfub(jl))/(zheat(jl)*ztau)
           zmfub1(jl) = max(zmfub1(jl),0.001)
           zmfmax=(paph(jl,ikb)-paph(jl,ikb-1))*zcons2
           zmfub1(jl)=min(zmfub1(jl),zmfmax)
       end if
      end do
!
!*  6.2   recalculate convective fluxes due to effect of
!         downdrafts on boundary layer moist static energy budget (ktype=2)
!--------------------------------------------------------
        do jl=1,nxj            
         if(ldcum(jl) .and. ktype(jl) .eq. 2) then
           ikb=kcbot(jl)
           if(pmfd(jl,ikb).lt.0.0 .and. loddraf(jl)) then
              zeps=-pmfd(jl,ikb)/max(zmfub(jl),cmfcmin)
           else
              zeps=0.
           endif
           zqumqe=pqu(jl,ikb)+plu(jl,ikb)-  &
     &            zeps*zqd(jl,ikb)-(1.-zeps)*zqenh(jl,ikb)
           zdqmin=max(0.01*zqenh(jl,ikb),cmfcmin)
           zmfmax=(paph(jl,ikb)-paph(jl,ikb-1))*zcons2
!  using moist static engergy closure instead of moisture closure
           zdh=cpd*(ptu(jl,ikb)-zeps*ztd(jl,ikb)- &
     &       (1.-zeps)*ztenh(jl,ikb))+alv*zqumqe
           zdh=g*max(zdh,1.e5*zdqmin)
           if(zdhpbl(jl).gt.0.)then
             zmfub1(jl)=zdhpbl(jl)/zdh
           else
             zmfub1(jl) = zmfub(jl)
           end if
           zmfub1(jl) = min(zmfub1(jl),zmfmax)
         end if

!*  6.3   mid-level convection - nothing special
!---------------------------------------------------------
         if(ldcum(jl) .and. ktype(jl) .eq. 3 ) then
            zmfub1(jl) = zmfub(jl)
         end if

       end do

!*  6.4   scaling the downdraft mass flux
!---------------------------------------------------------
       do jk=1,klev
        do jl = 1, nxj           
        if( ldcum(jl) ) then
           zfac=zmfub1(jl)/max(zmfub(jl),cmfcmin)
           pmfd(jl,jk)=pmfd(jl,jk)*zfac
           zmfds(jl,jk)=zmfds(jl,jk)*zfac
           zmfdq(jl,jk)=zmfdq(jl,jk)*zfac
           zdmfdp(jl,jk)=zdmfdp(jl,jk)*zfac
           pmfdde_rate(jl,jk) = pmfdde_rate(jl,jk)*zfac
        end if
       end do
       end do

!*  6.5   scaling the updraft mass flux
! --------------------------------------------------------
      do jl = 1, nxj        
        if ( ldcum(jl) ) zmfs(jl) = zmfub1(jl)/max(cmfcmin,zmfub(jl))
      end do
      do jk = 2 , klev
       do jl = 1, nxj      
        if ( ldcum(jl) .and. jk >= kctop(jl)-1 ) then
          ikb = kcbot(jl)
          if ( jk>ikb ) then
            zdz = ((paph(jl,klev+1)-paph(jl,jk))/(paph(jl,klev+1)-paph(jl,ikb)))
            pmfu(jl,jk) = pmfu(jl,ikb)*zdz
          end if
          zmfmax = (paph(jl,jk)-paph(jl,jk-1))*zcons2
          if ( pmfu(jl,jk)*zmfs(jl) > zmfmax ) then
            zmfs(jl) = min(zmfs(jl),zmfmax/pmfu(jl,jk))
          end if
        end if
      end do
      end do
      do jk = 2 , klev
        do jl = 1, nxj        
        if ( ldcum(jl) .and. jk <= kcbot(jl) .and. jk >= kctop(jl)-1 ) then
          pmfu(jl,jk) = pmfu(jl,jk)*zmfs(jl)
          zmfus(jl,jk) = zmfus(jl,jk)*zmfs(jl)
          zmfuq(jl,jk) = zmfuq(jl,jk)*zmfs(jl)
          zmful(jl,jk) = zmful(jl,jk)*zmfs(jl)
          zdmfup(jl,jk) = zdmfup(jl,jk)*zmfs(jl)
          plude(jl,jk) = plude(jl,jk)*zmfs(jl)
          pmfude_rate(jl,jk) = pmfude_rate(jl,jk)*zmfs(jl)
        end if
      end do
      end do

!*    6.6  if ktype = 2, kcbot=kctop is not allowed
! ---------------------------------------------------
      do jl = 1, nxj            
      if ( ktype(jl) == 2 .and. &
           kcbot(jl) == kctop(jl) .and. kcbot(jl) >= klev-1 ) then
        ldcum(jl) = .false.
        ktype(jl) = 0
      end if
      end do

      if ( .not. lmfscv .or. .not. lmfpen ) then
       do jl = 1, nxj 
        llo2(jl) = .false.
        if ( (.not. lmfscv .and. ktype(jl) == 2) .or. &
             (.not. lmfpen .and. ktype(jl) == 1) ) then
          llo2(jl) = .true.
          ldcum(jl) = .false.
        end if
      end do
      end if

!*   6.7  set downdraft mass fluxes to zero above cloud top
!----------------------------------------------------
      do jl = 1, nxj         
      if ( loddraf(jl) .and. idtop(jl) <= kctop(jl) ) then
        idtop(jl) = kctop(jl) + 1
      end if
      end do
      do jk = 2 , klev
        do jl = 1, nxj        
        if ( loddraf(jl) ) then
          if ( jk < idtop(jl) ) then
            pmfd(jl,jk) = 0.
            zmfds(jl,jk) = 0.
            zmfdq(jl,jk) = 0.
            pmfdde_rate(jl,jk) = 0.
            zdmfdp(jl,jk) = 0.
          else if ( jk == idtop(jl) ) then
            pmfdde_rate(jl,jk) = 0.
          end if
        end if
      end do
      end do

!----------------------------------------------------------
!*    7.0      determine final convective fluxes in 'cuflx'
!----------------------------------------------------------
       call cuflx_n                                        &
     &  (  nxj,      klon,     klev,     ztmst             &   
     &  ,  pten,     pqen,     pqsen,    ztenh,   zqenh    &        
     &  ,  paph,     pap,      zgeoh,    lndj,    ldcum    &      
     &  ,  kcbot,    kctop,    idtop,    itopm2            &      
     &  ,  ktype,    loddraf                               &       
     &  ,  pmfu,     pmfd,     zmfus,    zmfds             &      
     &  ,  zmfuq,    zmfdq,    zmful,    plude             &     
     &  ,  zdmfup,   zdmfdp,   zdpmel,   zlglac            &   
     &  ,  prain,    pmfdde_rate, pmflxr, pmflxs )    

!xb110, WRF 4.0 version>
! some adjustments needed
adj=.true.
if (adj)then
    do jl=1,nxj
      zmfs(jl) = 1.
      zmfuub(jl)=0.
    end do
    do jk = 2 , klev
      do jl = 1,nxj
        if ( loddraf(jl) .and. jk >= idtop(jl)-1 ) then   ! below cloud top level
          zmfmax = pmfu(jl,jk)*0.98
          if ( pmfd(jl,jk)+zmfmax+1.e-15 < 0. ) then
            zmfs(jl) = min(zmfs(jl),-zmfmax/pmfd(jl,jk))   ! factor for downdraft adjustment
          end if
        end if
      end do
    end do

    do jk = 2 , klev
      do jl = 1 , nxj
        if ( zmfs(jl) < 1. .and. jk >= idtop(jl)-1 ) then
          pmfd(jl,jk) = pmfd(jl,jk)*zmfs(jl)
          zmfds(jl,jk) = zmfds(jl,jk)*zmfs(jl)
          zmfdq(jl,jk) = zmfdq(jl,jk)*zmfs(jl)
          pmfdde_rate(jl,jk) = pmfdde_rate(jl,jk)*zmfs(jl)
          zmfuub(jl) = zmfuub(jl) - (1.-zmfs(jl))*zdmfdp(jl,jk)
          pmflxr(jl,jk+1) = pmflxr(jl,jk+1) + zmfuub(jl)
          zdmfdp(jl,jk) = zdmfdp(jl,jk)*zmfs(jl)
        end if
      end do
    end do

    do jk = 2 , klev - 1
      do jl = 1, nxj
        if ( loddraf(jl) .and. jk >= idtop(jl)-1 ) then
          zerate = -pmfd(jl,jk) + pmfd(jl,jk-1) + pmfdde_rate(jl,jk)
          if ( zerate < 0. ) then
            pmfdde_rate(jl,jk) = pmfdde_rate(jl,jk) - zerate
          end if
        end if
        if ( ldcum(jl) .and. jk >= kctop(jl)-1 ) then
          zerate = pmfu(jl,jk) - pmfu(jl,jk+1) + pmfude_rate(jl,jk)
          if ( zerate < 0. ) then
            pmfude_rate(jl,jk) = pmfude_rate(jl,jk) - zerate
          end if
          zdmfup(jl,jk) = pmflxr(jl,jk+1) + pmflxs(jl,jk+1) - &
                          pmflxr(jl,jk) - pmflxs(jl,jk)
          zdmfdp(jl,jk) = 0.
        end if
      end do
    end do
end if       !end for adjustment

! avoid negative humidities at ddraught top
    do jl = 1,nxj
      if ( loddraf(jl) ) then
        jk = idtop(jl)
        ik = min(jk+1,klev)
        if ( zmfdq(jl,jk) < 0.3*zmfdq(jl,ik) ) then
            zmfdq(jl,jk) = 0.3*zmfdq(jl,ik)
        end if
      end if
    end do

! avoid negative humidities near cloud top because gradient of precip flux
! and detrainment / liquid water flux are too large
    do jk = 2 , klev
      do jl = 1, nxj
        if ( ldcum(jl) .and. jk >= kctop(jl)-1 .and. jk < kcbot(jl) ) then
          zdz = ztmst*g/(paph(jl,jk+1)-paph(jl,jk))
          zmfa = zmfuq(jl,jk+1) + zmfdq(jl,jk+1) - &
                 zmfuq(jl,jk) - zmfdq(jl,jk) + &
                 zmful(jl,jk+1) - zmful(jl,jk) + zdmfup(jl,jk)
          zmfa = (zmfa-plude(jl,jk))*zdz
          if ( pqen(jl,jk)+zmfa < 0. ) then
            plude(jl,jk) = plude(jl,jk) + 2.*(pqen(jl,jk)+zmfa)/zdz
          end if
          if ( plude(jl,jk) < 0. ) plude(jl,jk) = 0.
        end if
        if ( .not. ldcum(jl) ) pmfude_rate(jl,jk) = 0.
        if ( abs(pmfd(jl,jk-1)) < 1.0e-20 ) pmfdde_rate(jl,jk) = 0.
      end do
    end do
!xb110, WRF<

      do jl = 1, nxj           
        prsfc(jl) = pmflxr(jl,klev+1)
        pssfc(jl) = pmflxs(jl,klev+1)
      end do

!xb110>
!for lightning parameterization from ECMWF
      pf = 0.
      rho = 0.

      do jk = 1,klev
       do jl = 1,nxj
        pf(jl,jk) = pmflxs(jl,jk)
        rho(jl,jk) = pap(jl,jk)/(rd*ztenh(jl,jk))
!environment air density
        papn(jl,jk) = pap(jl,jk)* 0.01 !Pa to mb
       end do
      end do
    
!      if (myrank .eq. 0)then
!        print*, "max sflx:",maxval(pf),"min sflx:",minval(pf)
!      end if

!xb110<

!----------------------------------------------------------------
!*    8.0      update tendencies for t and q in subroutine cudtdq
!----------------------------------------------------------------
      call cudtdq_n(nxj,klon,klev,itopm2,kctop,idtop,ldcum,loddraf, &
                 ztmst,paph,zgeoh,pgeo,pten,ztenh,pqen,zqenh,pqsen,     &
                 zlglac,plude,pmfu,pmfd,zmfus,zmfds,zmfuq,zmfdq,zmful,   &
                 zdmfup,zdmfdp,zdpmel,ptte,pqte,pcte,pxtec)
!----------------------------------------------------------------
!*    9.0      update tendencies for u and u in subroutine cududv
!----------------------------------------------------------------
      if(lmfdudv) then
      do jk = klev-1 , 2 , -1
        ik = jk + 1
        do jl = 1, nxj          
          if ( ldcum(jl) ) then
            if ( jk == kcbot(jl) .and. ktype(jl) < 3 ) then
              ikb = kdpl(jl)
              zuu(jl,jk) = puen(jl,ikb-1)
              zvu(jl,jk) = pven(jl,ikb-1)
            else if ( jk == kcbot(jl) .and. ktype(jl) == 3 ) then
              zuu(jl,jk) = puen(jl,jk-1)
              zvu(jl,jk) = pven(jl,jk-1)
            end if
            if ( jk < kcbot(jl) .and. jk >= kctop(jl) ) then
            if(momtrans .eq. 1)then
              zfac = 0.
              if ( ktype(jl) == 1 .or. ktype(jl) == 3 ) zfac = 2.
              if ( ktype(jl) == 1 .and. jk <= kctop(jl)+2 ) zfac = 3.
              zerate = pmfu(jl,jk) - pmfu(jl,ik) + &
                (1.+zfac)*pmfude_rate(jl,jk)
              zderate = (1.+zfac)*pmfude_rate(jl,jk)
              zmfa = 1./max(cmfcmin,pmfu(jl,jk))
              zuu(jl,jk) = (zuu(jl,ik)*pmfu(jl,ik) + &
                zerate*puen(jl,jk)-zderate*zuu(jl,ik))*zmfa
              zvu(jl,jk) = (zvu(jl,ik)*pmfu(jl,ik) + &
                zerate*pven(jl,jk)-zderate*zvu(jl,ik))*zmfa
            else
              if(ktype(jl) == 1 .or. ktype(jl) == 3) then
                pgf_u   =  -0.7*0.5*(pmfu(jl,ik)*(puen(jl,ik)-puen(jl,jk))+&
                                   pmfu(jl,jk)*(puen(jl,jk)-puen(jl,jk-1)))
                pgf_v   =  -0.7*0.5*(pmfu(jl,ik)*(pven(jl,ik)-pven(jl,jk))+&
                                   pmfu(jl,jk)*(pven(jl,jk)-pven(jl,jk-1)))
              else
                pgf_u   = 0.
                pgf_v   = 0.
              end if
              zerate = pmfu(jl,jk) - pmfu(jl,ik) + pmfude_rate(jl,jk)
              zderate = pmfude_rate(jl,jk)
              zmfa = 1./max(cmfcmin,pmfu(jl,jk))
              zuu(jl,jk) = (zuu(jl,ik)*pmfu(jl,ik) + &
                zerate*puen(jl,jk)-zderate*zuu(jl,ik)+pgf_u)*zmfa
              zvu(jl,jk) = (zvu(jl,ik)*pmfu(jl,ik) + &
                zerate*pven(jl,jk)-zderate*zvu(jl,ik)+pgf_v)*zmfa
            end if
            end if
          end if
        end do
      end do

      if(lmfdd) then
      do jk = 3 , klev
        ik = jk - 1
        do jl = 1, nxj
          if ( ldcum(jl) ) then
            if ( jk == idtop(jl) ) then
              zud(jl,jk) = 0.5*(zuu(jl,jk)+puen(jl,ik))
              zvd(jl,jk) = 0.5*(zvu(jl,jk)+pven(jl,ik))
            else if ( jk > idtop(jl) ) then
              zerate = -pmfd(jl,jk) + pmfd(jl,ik) + pmfdde_rate(jl,jk)
              zmfa = 1./min(-cmfcmin,pmfd(jl,jk))
              zud(jl,jk) = (zud(jl,ik)*pmfd(jl,ik) - &
                zerate*puen(jl,ik)+pmfdde_rate(jl,jk)*zud(jl,ik))*zmfa
              zvd(jl,jk) = (zvd(jl,ik)*pmfd(jl,ik) - &
                zerate*pven(jl,ik)+pmfdde_rate(jl,jk)*zvd(jl,ik))*zmfa
            end if
          end if
        end do
      end do
      end if
!   --------------------------------------------------
!   rescale massfluxes for stability in Momentum
!------------------------------------------------------------------------
      zmfs(:) = 1.
      do jk = 2 , klev
        do jl = 1, nxj          
          if ( ldcum(jl) .and. jk >= kctop(jl)-1 ) then
            zmfmax = (paph(jl,jk)-paph(jl,jk-1))*zcons
            if ( pmfu(jl,jk) > zmfmax .and. jk >= kctop(jl) ) then
              zmfs(jl) = min(zmfs(jl),zmfmax/pmfu(jl,jk))
            end if
          end if
        end do
      end do
      do jk = 1 , klev
        do jl = 1, nxj           
          zmfuus(jl,jk) = pmfu(jl,jk)
          zmfdus(jl,jk) = pmfd(jl,jk)
          if ( ldcum(jl) .and. jk >= kctop(jl)-1 ) then
            zmfuus(jl,jk) = pmfu(jl,jk)*zmfs(jl)
            zmfdus(jl,jk) = pmfd(jl,jk)*zmfs(jl)
          end if
        end do
      end do
!*  9.1          update u and v in subroutine cududvn
!-------------------------------------------------------------------
     do jk = 1 , klev
       do jl = 1, nxj        
          ztenu(jl,jk) = pvom(jl,jk)
          ztenv(jl,jk) = pvol(jl,jk)
        end do
      end do

      call cududv_n(nxj, klon,klev,itopm2,ktype,kcbot,kctop, &
                  ldcum,ztmst,paph,puen,pven,zmfuus,zmfdus,zuu,  &
                  zud,zvu,zvd,pvom,pvol,jin)
!  calculate KE dissipation
      do jl = 1, nxj          
        zsum12(jl) = 0.
        zsum22(jl) = 0.
      end do
        do jk = 1 , klev
          do jl = 1, nxj             
            zuv2(jl,jk) = 0.
            if ( ldcum(jl) .and. jk >= kctop(jl)-1 ) then
              zdz = (paph(jl,jk+1)-paph(jl,jk))
              zduten = pvom(jl,jk) - ztenu(jl,jk)
              zdvten = pvol(jl,jk) - ztenv(jl,jk)
              zuv2(jl,jk) = sqrt(zduten**2+zdvten**2)
              zsum22(jl) = zsum22(jl) + zuv2(jl,jk)*zdz
              zsum12(jl) = zsum12(jl) - &
                (puen(jl,jk)*zduten+pven(jl,jk)*zdvten)*zdz
            end if
          end do
        end do
        do jk = 1 , klev
          do jl = 1, nxj           
            if ( ldcum(jl) .and. jk>=kctop(jl)-1 ) then
              ztdis = rcpd*zsum12(jl)*zuv2(jl,jk)/max(1.e-15,zsum22(jl))
              ptte(jl,jk) = ptte(jl,jk) + ztdis
            end if
          end do
        end do

      end if
!----------------------------------------------------------------------
!*   10.           IN CASE THAT EITHER DEEP OR SHALLOW IS SWITCHED OFF
! NEED TO SET SOME VARIABLES A POSTERIORI TO ZERO
! ---------------------------------------------------
      if ( .not. lmfscv .or. .not. lmfpen ) then
      do jk = 2 , klev
        do jl = 1, nxj          
          if ( llo2(jl) .and. jk >= kctop(jl)-1 ) then
            ptu(jl,jk) = pten(jl,jk)
            pqu(jl,jk) = pqen(jl,jk)
            plu(jl,jk) = 0.
            pmfude_rate(jl,jk) = 0.
            pmfdde_rate(jl,jk) = 0.
          end if
        end do
      end do
       do jl = 1, nxj         
        if ( llo2(jl) ) then
          kctop(jl) = -1
          kcbot(jl) = -1
        end if
      end do
      end if

!xb110>
      kcnv = 0
      do jl = 1,nxj
       if (ktype(jl) .eq. 1 .or. ktype(jl) .eq. 3 ) kcnv(jl) = 1
       if (ktype(jl) .eq. 2 ) kcnv(jl) = 0
      end do
      
      call lightning_ec (nxj,klon,klev,ptu,pqu,ztenh,zqenh     &
                        ,plu,rho,kcbot,kctop,pgeo,papn,pf  &
                        ,flash,lndj,kcnv)
!xb110<

      return
      end subroutine cumastr_n
