      subroutine cutype_n & 
     &   ( nxj,     klon,    klev,     klevp1,   klevm1,     pqen,&
     &     ptenh,   pqenh,     pqsenh,    pgeoh,     paph,& 
     &       hfx,     qfx,       pgeo,    pqsen,      pap,&
     &      pten,    lndj,       cutu,     cuqu,    culab,&
     &     ldcum,   cubot,      cutop,    ktype,    wbase,    culu,   kdpl,&
     &     p650,     jin    )
!      zhang & wang      iprc           2011-2013
!***purpose.
!   --------
!          to produce first guess updraught for cu-parameterizations
!          calculates condensation level, and sets updraught base variables and
!          first guess cloud type
!***interface
!   ---------
!          this routine is called from *cumastr*.
!          input are environm. values of t,q,p,phi at half levels.
!          it returns cloud types as follows;
!                 ktype=1 for deep cumulus
!                 ktype=2 for shallow cumulus
!***method.
!  --------
!          based on a simplified updraught equation
!            partial(hup)/partial(z)=eta(h - hup)
!            eta is the entrainment rate for test parcel
!            h stands for dry static energy or the total water specific humidity
!            references: christian jakob, 2003: a new subcloud model for
!            mass-flux convection schemes
!                        influence on triggering, updraft properties, and model
!                        climate, mon.wea.rev.
!                        131, 2765-2778
!            and
!                        ifs documentation - cy36r1,cy38r1 
!***input variables:
!       ptenh [ztenh] - environment temperature on half levels. (cuini)
!       pqenh [zqenh] - env. specific humidity on half levels. (cuini)
!       pgeoh [zgeoh] - geopotential on half levels, (mssflx)
!       paph - pressure of half levels. (mssflx)
!       rho  - density of the lowest model level
!       qfx  - net upward moisture flux at the surface (kg/m^2/s)
!       hfx  - net upward heat flux at the surface (w/m^2)
!***variables output by cutype:
!       ktype - convection type - 1: penetrative  (cumastr)
!                                 2: stratocumulus (cumastr)
!                                 3: mid-level  (cuasc)
!       information for updraft parcel (ptu,pqu,plu,kcbot,klab,kdpl...)
! ----------------------------------------------------------------
  USE mo_constants,    ONLY: cpd,    &! specific heat at constant pressure
                             rcpd,   &! rcpd=1./cpd
                             vtmpc1, &! vtmpc1=rv/rd-1
                             zrg,    &! 1.0/g
                             ralfdcp,&
                             rd,     &
                             g,      &
                             c4les,  &
                             c5les,  &
                             c4ies,  &
                             c5ies
  USE mo_cumulus_flux, ONLY: lmfdudv, & ! true if cumulus friction is switched on
                             zdnoprc

!-------------------------------------------------------------------
      implicit none
!-------------------------------------------------------------------
      integer  klon, klev, klevp1, klevm1,nxj,p650,jin
      real     foeewm,foealfa
      real     ptenh(klon,klev),       pqenh(klon,klev),& 
     &         pqsen(klon,klev),       pqsenh(klon,klev),&
     &         pgeoh(klon,klevp1),     paph(klon,klevp1),&
     &         pap(klon,klev),         pqen(klon,klev)
      real     pten(klon,klev)
      real     ptu(klon,klev),pqu(klon,klev),plu(klon,klev)
      real     pgeo(klon,klev)
      integer  klab(klon,klev)
      integer  kctop(klon),kcbot(klon)

      real     qfx(klon),hfx(klon)
      real     zph(klon)
      integer  lndj(klon)
      logical  loflag(klon), deepflag(klon), resetflag(klon)

! output variables
      real     cutu(klon,klev), cuqu(klon,klev), culu(klon,klev)
      integer  culab(klon,klev)
      real     wbase(klon)
      integer  ktype(klon),cubot(klon),cutop(klon),kdpl(klon)
      logical  ldcum(klon)

! local variables
      real     zqold(klon)
      real     rho, part1, part2, root, conw, deltt, deltq
      real     eta(klon),dz(klon),coef(klon)
      real     dhen(klon,klev), dh(klon,klev)
      real     plude(klon,klev)
      real     kup(klon,klev)
      real     vptu(klon,klev),vten(klon,klev)
      real     zbuo(klon,klev),abuoy(klon,klev)
  
      real     zz,zdken,zdq
      real     fscale,crirh1,pp
      real     atop1,atop2,abot
      real     tmix,zmix,qmix,pmix
      real     zlglac,dp,t13
      integer  nk,is,ikb,ikt

      real     zqsu,zcor,zdp,zesdp,zalfaw,zfacw,zfaci,zfac,zdsdp,zdqsdt,zdtdp
      real     zpdifftop, zpdiffbot
      integer  zcbase(klon)
      integer  jl,jk,ik,icall,levels
      logical  needreset, lldcum(klon)
!--------------------------------------------------------------
      t13 = 1.0/3.0
!xb110>
      eta = 0.
      dz  = 0.
      coef = 0.
      zph  = 0.
      zcbase = 0.

      cubot = -1
      cutop = -1
!xb110<
!
      do jl = 1, nxj            
        kcbot(jl)=klev
        kctop(jl)=klev
        kdpl(jl) =klev
        ktype(jl)=0
        wbase(jl)=0.
        ldcum(jl)=.false.
      end do

!-----------------------------------------------------------
! let's do test,and check the shallow convection first
! the first level is klev
! define deltat and deltaq
!-----------------------------------------------------------
      do jk=1,klev
       do jl = 1, nxj            
          plu(jl,jk)=culu(jl,jk)  ! parcel liquid water
          ptu(jl,jk)=cutu(jl,jk)  ! parcel temperature
          pqu(jl,jk)=cuqu(jl,jk)  ! parcel specific humidity
          klab(jl,jk)=culab(jl,jk)
           dh(jl,jk)=0.0  ! parcel dry static energy
         dhen(jl,jk)=0.0  ! environment dry static energy
          kup(jl,jk)=0.0  ! updraught kinetic energy for parcel
         vptu(jl,jk)=0.0  ! parcel virtual temperature considering water-loading
         vten(jl,jk)=0.0  ! environment virtual temperature
         zbuo(jl,jk)=0.0  ! parcel buoyancy
         abuoy(jl,jk)=0.0
      end do
      end do

       do jl = 1, nxj             
         zqold(jl) = 0.
         lldcum(jl) = .false.
         loflag(jl) = .true.
      end do

! check the levels from lowest level to second top level
      do jk=klevm1,2,-1

! define the variables at the first level      
      if(jk .eq. klevm1) then
       do jl = 1, nxj             
        rho=paph(jl,klev+1)/ &
     &         (rd*(pten(jl,klev)*(1.+vtmpc1*pqen(jl,klev))))
        part1 = 1.5*0.4*(pgeoh(jl,klev)-pgeoh(jl,klev+1))/ &
     &              (rho*pten(jl,klev))
        part2 = -hfx(jl)*rcpd-vtmpc1*pten(jl,klev)*qfx(jl)
        root  = 0.001-part1*part2
        if(part2 .lt. 0.) then
           conw  = 1.2*(root)**t13
           deltt = max(1.5*hfx(jl)/(rho*cpd*conw),0.)
           deltq = max(1.5*qfx(jl)/(rho*conw),0.)
           kup(jl,klev) = 0.5*(conw**2)
           pqu(jl,klev)= pqenh(jl,klev) + deltq
          dhen(jl,klev)= pgeoh(jl,klev) + ptenh(jl,klev)*cpd
           dh(jl,klev) = dhen(jl,klev)  + deltt*cpd
          ptu(jl,klev) = (dh(jl,klev)-pgeoh(jl,klev))*rcpd 
          vptu(jl,klev)=ptu(jl,klev)*(1.+vtmpc1*pqu(jl,klev))
          vten(jl,klev)=ptenh(jl,klev)*(1.+vtmpc1*pqenh(jl,klev))
          zbuo(jl,klev)=(vptu(jl,klev)-vten(jl,klev))/vten(jl,klev)
          klab(jl,klev) = 1
        else
          loflag(jl) = .false.
        end if
      end do
      end if
 
      is=0
       do jl = 1, nxj          
         if(loflag(jl))then
            is=is+1
         endif
      enddo
      if(is.eq.0) exit

! the next levels, we use the variables at the first level as initial values
       do jl = 1, nxj           
      if(loflag(jl)) then
!         eta(jl) = 0.55/((pgeoh(jl,jk)-pgeoh(jl,klev+1))*zrg)+1.0e-4
        eta(jl) = 0.8/(pgeo(jl,jk)*zrg)+2.e-4             !xb110, for more stable
        dz(jl)  = (pgeoh(jl,jk)-pgeoh(jl,jk+1))*zrg
        coef(jl)= 0.5*eta(jl)*dz(jl)
        dhen(jl,jk) = pgeoh(jl,jk) + cpd*ptenh(jl,jk)
        dh(jl,jk) = (coef(jl)*(dhen(jl,jk+1)+dhen(jl,jk))&
     &              +(1.-coef(jl))*dh(jl,jk+1))/(1.+coef(jl))
        pqu(jl,jk) =(coef(jl)*(pqenh(jl,jk+1)+pqenh(jl,jk))&
     &              +(1.-coef(jl))*pqu(jl,jk+1))/(1.+coef(jl))
        ptu(jl,jk) = (dh(jl,jk)-pgeoh(jl,jk))*rcpd
        zqold(jl) = pqu(jl,jk)
        zph(jl)=paph(jl,jk)
      end if
      end do
! check if the parcel is saturated
      ik=jk
      icall=1
      call cuadjtq_n(nxj,klon,klev,ik,zph,ptu,pqu,loflag,icall)
        do jl = 1, nxj          
        if( loflag(jl) ) then
          pqu(jl,jk) =max(pqu(jl,jk),0.) !xb110, to make sure that pqu is positive
          zdq = max((zqold(jl) - pqu(jl,jk)),0.)
          plu(jl,jk) = plu(jl,jk+1) + zdq
          zlglac=zdq*((1.-foealfa(ptu(jl,jk))) - &
                      (1.-foealfa(ptu(jl,jk+1))))
          plu(jl,jk) = min(plu(jl,jk),5.e-3)
          dh(jl,jk) =  pgeoh(jl,jk) + cpd*(ptu(jl,jk)+ralfdcp*zlglac)
! compute the updraft speed
          vptu(jl,jk) = ptu(jl,jk)*(1.+vtmpc1*pqu(jl,jk)-plu(jl,jk))+&
                        ralfdcp*zlglac
          vten(jl,jk) = ptenh(jl,jk)*(1.+vtmpc1*pqenh(jl,jk))
          zbuo(jl,jk) = (vptu(jl,jk) - vten(jl,jk))/vten(jl,jk)
          abuoy(jl,jk)=(zbuo(jl,jk)+zbuo(jl,jk+1))*0.5*g
          atop1 = 1.0 - 2.*coef(jl)
          atop2 = 2.0*dz(jl)*abuoy(jl,jk)
          abot =  1.0 + 2.*coef(jl)
          kup(jl,jk)  = (atop1*kup(jl,jk+1) + atop2) / abot

! let's find the exact cloud base
         if ( plu(jl,jk) > 0. .and. klab(jl,jk+1) == 1 ) then
              ik = jk + 1
              zqsu = foeewm(ptu(jl,ik))/paph(jl,ik)
              zqsu = min(0.5,zqsu)
              zcor = 1./(1.-vtmpc1*zqsu)
              zqsu = zqsu*zcor
              zdq = min(0.,pqu(jl,ik)-zqsu)
              zalfaw = foealfa(ptu(jl,ik))
              zfacw = c5les/((ptu(jl,ik)-c4les)**2)
              zfaci = c5ies/((ptu(jl,ik)-c4ies)**2)
              zfac = zalfaw*zfacw + (1.-zalfaw)*zfaci
              zesdp = foeewm(ptu(jl,ik))/paph(jl,ik)
              zcor = 1./(1.-vtmpc1*zesdp)
              zdqsdt = zfac*zcor*zqsu
              zdtdp = rd*ptu(jl,ik)/(cpd*paph(jl,ik))
              zdp = zdq/(zdqsdt*zdtdp)
              zcbase(jl) = paph(jl,ik) + zdp
! chose nearest half level as cloud base (jk or jk+1)
              zpdifftop = zcbase(jl) - paph(jl,jk)
              zpdiffbot = paph(jl,jk+1) - zcbase(jl)
              if ( zpdifftop > zpdiffbot .and. kup(jl,jk+1) > 0. ) then
                ikb = min(klev-1,jk+1)
                klab(jl,ikb) = 2
                klab(jl,jk) = 2
                kcbot(jl) = ikb
                plu(jl,jk+1) = 1.0e-8
              else if ( zpdifftop <= zpdiffbot .and.kup(jl,jk) > 0. ) then
                klab(jl,jk) = 2
                kcbot(jl) = jk
              end if
          end if

          if(kup(jl,jk) .lt. 0.)then
            loflag(jl) = .false.
            if(plu(jl,jk+1) .gt. 0.) then
              kctop(jl) = jk
              lldcum(jl) = .true.
            else
              lldcum(jl) = .false.
            end if
          else if(plu(jl,jk) .gt. 0.)then
              klab(jl,jk)=2
          else
              klab(jl,jk)=1
          end if
        end if
      end do

      end do ! end all the levels

      do jl = 1, nxj           
        ikb = kcbot(jl)
        ikt = kctop(jl)
        if(paph(jl,ikb) - paph(jl,ikt) > zdnoprc) lldcum(jl) = .false.
        if(lldcum(jl)) then
            ktype(jl) = 2
            ldcum(jl) = .true.
            wbase(jl) = sqrt(max(2.*kup(jl,ikb),0.))
            cubot(jl) = ikb
            cutop(jl) = ikt
            kdpl(jl)  = klev
        else
            cutop(jl) = -1
            cubot(jl) = -1
            kdpl(jl)  = klev - 1
            ldcum(jl) = .false.
            wbase(jl) = 0.
        end if
      end do

      do jk=klev,1,-1
        do jl = 1, nxj          
           ikt = kctop(jl)
           if(jk .ge. ikt)then
             culab(jl,jk) = klab(jl,jk)
             cutu(jl,jk)  = ptu(jl,jk)
             cuqu(jl,jk)  = pqu(jl,jk)
             culu(jl,jk)  = plu(jl,jk)
           end if
       end do
      end do
      
!-----------------------------------------------------------
! next, let's check the deep convection
! the first level is klevm1-1
! define deltat and deltaq
!----------------------------------------------------------
! we check the parcel starting level by level
! assume the mix-layer is 60hPa
      deltt = 0.2
      deltq = 1.0e-4
      do jl = 1, nxj           
        deepflag(jl) = .false.
      end do

      do levels=klevm1-1,p650,-1 ! loop starts
        do jk=1,klev
           do jl = 1, nxj              
             plu(jl,jk)=0.0  ! parcel liquid water
             ptu(jl,jk)=0.0  ! parcel temperature
             pqu(jl,jk)=0.0  ! parcel specific humidity
             dh(jl,jk)=0.0   ! parcel dry static energy
             dhen(jl,jk)=0.0  ! environment dry static energy
             kup(jl,jk)=0.0   ! updraught kinetic energy for parcel
             vptu(jl,jk)=0.0  ! parcel virtual temperature consideringwater-loading
             vten(jl,jk)=0.0  ! environment virtual temperature
             abuoy(jl,jk)=0.0
             zbuo(jl,jk)=0.0
             klab(jl,jk)=0
          end do
        end do

         do jl = 1, nxj              
           kcbot(jl)    =  levels
           kctop(jl)    =  levels
           zqold(jl)    = 0.
           lldcum(jl)   = .false.
           resetflag(jl)= .false.
           loflag(jl)   = .not. deepflag(jl)
        end do

! start the inner loop to search the deep convection points
      do jk=levels,2,-1
        is=0
         do jl = 1, nxj              
         if(loflag(jl))then
            is=is+1
         endif
        enddo
        if(is.eq.0) exit

! define the variables at the departure level 
        if(jk .eq. levels) then
         do jl = 1,nxj           
          if(loflag(jl)) then
            if((paph(jl,klev+1)-paph(jl,jk)) < 60.e2) then
              tmix=0.
              qmix=0.
              zmix=0.
              pmix=0.
              do nk=jk+2,jk,-1
              if(pmix < 50.e2) then
                dp = paph(jl,nk) - paph(jl,nk-1)
                tmix=tmix+dp*ptenh(jl,nk)
                qmix=qmix+dp*pqenh(jl,nk)
                zmix=zmix+dp*pgeoh(jl,nk)
                pmix=pmix+dp
              end if
              end do
              tmix=tmix/pmix
              qmix=qmix/pmix
              zmix=zmix/pmix
            else
              tmix=ptenh(jl,jk+1)
              qmix=pqenh(jl,jk+1)
              zmix=pgeoh(jl,jk+1)
            end if

            pqu(jl,jk+1) = qmix + deltq
            dhen(jl,jk+1)= zmix + tmix*cpd
            dh(jl,jk+1)  = dhen(jl,jk+1) + deltt*cpd
            ptu(jl,jk+1) = (dh(jl,jk+1)-pgeoh(jl,jk+1))*rcpd
            kup(jl,jk+1) = 0.5
            klab(jl,jk+1)= 1
            vptu(jl,jk+1)=ptu(jl,jk+1)*(1.+vtmpc1*pqu(jl,jk+1))
            vten(jl,jk+1)=ptenh(jl,jk+1)*(1.+vtmpc1*pqenh(jl,jk+1))
            zbuo(jl,jk+1)=(vptu(jl,jk+1)-vten(jl,jk+1))/vten(jl,jk+1)
          end if
          end do
        end if

! the next levels, we use the variables at the first level as initial values
         do jl = 1, nxj                   
           if(loflag(jl)) then
! define the fscale
             fscale = min(1.,(pqsen(jl,jk)/pqsen(jl,levels))**3)
             eta(jl) = 1.75e-3*fscale
             dz(jl)  = (pgeoh(jl,jk)-pgeoh(jl,jk+1))*zrg
             coef(jl)= 0.5*eta(jl)*dz(jl)
             dhen(jl,jk) = pgeoh(jl,jk) + cpd*ptenh(jl,jk)
             dh(jl,jk) = 0.5*coef(jl)*(dhen(jl,jk+1)+dhen(jl,jk))&
     &                  +(1.-coef(jl))*dh(jl,jk+1)
             pqu(jl,jk) =0.5*coef(jl)*(pqenh(jl,jk+1)+pqenh(jl,jk))&
     &                  +(1.-coef(jl))*pqu(jl,jk+1)
             ptu(jl,jk) = (dh(jl,jk)-pgeoh(jl,jk))*rcpd
             zqold(jl) = pqu(jl,jk)
             zph(jl)=paph(jl,jk)
           end if
        end do
! check if the parcel is saturated
        ik=jk
        icall=1
        call cuadjtq_n(nxj,klon,klev,ik,zph,ptu,pqu,loflag,icall)
 
        do jl = 1, nxj              
        if( loflag(jl) ) then
          pqu(jl,jk) =max(pqu(jl,jk),0.) !xb110, to make sure that pqu is positive
          zdq = max((zqold(jl) - pqu(jl,jk)),0.)
          plu(jl,jk) = plu(jl,jk+1) + zdq
          zlglac=zdq*((1.-foealfa(ptu(jl,jk))) - &
                      (1.-foealfa(ptu(jl,jk+1))))
          plu(jl,jk) = 0.5*plu(jl,jk)
          dh(jl,jk) =  pgeoh(jl,jk) + cpd*(ptu(jl,jk)+ralfdcp*zlglac)
! compute the updraft speed
          vptu(jl,jk) = ptu(jl,jk)*(1.+vtmpc1*pqu(jl,jk)-plu(jl,jk))+&
                        ralfdcp*zlglac
          vten(jl,jk) = ptenh(jl,jk)*(1.+vtmpc1*pqenh(jl,jk))
          zbuo(jl,jk) = (vptu(jl,jk) - vten(jl,jk))/vten(jl,jk)
          abuoy(jl,jk)=(zbuo(jl,jk)+zbuo(jl,jk+1))*0.5*g
          atop1 = 1.0 - 2.*coef(jl)
          atop2 = 2.0*dz(jl)*abuoy(jl,jk)
          abot =  1.0 + 2.*coef(jl)
          kup(jl,jk)  = (atop1*kup(jl,jk+1) + atop2) / abot
! let's find the exact cloud base
          if ( plu(jl,jk) > 0. .and. klab(jl,jk+1) == 1 ) then
              ik = jk + 1
              zqsu = foeewm(ptu(jl,ik))/paph(jl,ik)
              zqsu = min(0.5,zqsu)
              zcor = 1./(1.-vtmpc1*zqsu)
              zqsu = zqsu*zcor
              zdq = min(0.,pqu(jl,ik)-zqsu)
              zalfaw = foealfa(ptu(jl,ik))
              zfacw = c5les/((ptu(jl,ik)-c4les)**2)
              zfaci = c5ies/((ptu(jl,ik)-c4ies)**2)
              zfac = zalfaw*zfacw + (1.-zalfaw)*zfaci
              zesdp = foeewm(ptu(jl,ik))/paph(jl,ik)
              zcor = 1./(1.-vtmpc1*zesdp)
              zdqsdt = zfac*zcor*zqsu
              zdtdp = rd*ptu(jl,ik)/(cpd*paph(jl,ik))
              zdp = zdq/(zdqsdt*zdtdp)
              zcbase(jl) = paph(jl,ik) + zdp
! chose nearest half level as cloud base (jk or jk+1)
              zpdifftop = zcbase(jl) - paph(jl,jk)
              zpdiffbot = paph(jl,jk+1) - zcbase(jl)
              if ( zpdifftop > zpdiffbot .and. kup(jl,jk+1) > 0. ) then
                ikb = min(klev-1,jk+1)
                klab(jl,ikb) = 2
                klab(jl,jk) = 2
                kcbot(jl) = ikb
                plu(jl,jk+1) = 1.0e-8
              else if ( zpdifftop <= zpdiffbot .and.kup(jl,jk) > 0. ) then
                klab(jl,jk) = 2
                kcbot(jl) = jk
              end if
          end if

          if(kup(jl,jk) .lt. 0.)then
            loflag(jl) = .false.
            if(plu(jl,jk+1) .gt. 0.) then
              kctop(jl) = jk
              lldcum(jl) = .true.
            else
              lldcum(jl) = .false.
            end if
          else if(plu(jl,jk) .gt. 0.)then
              klab(jl,jk)=2
          else
              klab(jl,jk)=1
          end if
        end if
      end do

      end do ! end all the levels

      needreset = .false.
      do jl = 1, nxj              
        ikb = kcbot(jl)
        ikt = kctop(jl)
        if(paph(jl,ikb) - paph(jl,ikt) < zdnoprc) lldcum(jl) = .false.
        if(lldcum(jl)) then      
         ktype(jl)    = 1
         ldcum(jl)    = .true.
         deepflag(jl) = .true.
         wbase(jl)    = sqrt(max(2.*kup(jl,ikb),0.))
         cubot(jl)    = ikb
         cutop(jl)    = ikt
         kdpl(jl)     = levels+1
         needreset    = .true.
         resetflag(jl)= .true.
        end if
      end do

      if(needreset) then
      do jk=klev,1,-1
         do jl = 1, nxj                 
          if(resetflag(jl)) then
            ikt = kctop(jl)
            ikb = kdpl(jl)
            if(jk .le. ikb .and. jk .ge. ikt )then
             culab(jl,jk) = klab(jl,jk)
             cutu(jl,jk)  = ptu(jl,jk)
             cuqu(jl,jk)  = pqu(jl,jk)
             culu(jl,jk)  = plu(jl,jk)
            else
             culab(jl,jk) = 1
             cutu(jl,jk)  = ptenh(jl,jk)
             cuqu(jl,jk)  = pqenh(jl,jk)
             culu(jl,jk)  = 0.
            end if
            if ( jk .lt. ikt ) culab(jl,jk) = 0
          end if
        end do
      end do
      end if

      end do ! end all cycles

      return
      end subroutine cutype_n
