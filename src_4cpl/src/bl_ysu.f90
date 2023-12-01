!    call     ysu2d(u1,v1,t1,q1,prsl,prsi,prslk,pk2(1,lve)         &
!                  nx,nxj,lev,ntrac=2,del                          &
!                  cp,g,r,hltm,                                    & !constant
!                  phil,psi,                                       &
!                  z0,stress,hpbl,fm,fh,                           &
!                  slmsk,heat,evap,sfcw,rb,                        &
!                  dt,rcl,                                    &
!                  u10,v10,                                        &
!                  swh,hlw,xmu                                     &
!                   )
   subroutine ysu2d(ux,vx,tx,qxt,p2d,p2di,prslk,psk,              &
                  ix,im,km,ndiff,del,                             &
                  cp,g,rd,xlv,                                    &!const
!                  dz8w2d,psfcpa,                                  &
                  phii,phil,psfcpa,                                  &
                  znt,stress,hpbl,kpbl,psim,psih,                 &
                  ls,heat,evap,wspd,br,                           &
                  dt,rcl,                                         &
                  u10,v10,                                        &
                  swh,hlw,xmu,jj)
!-------------------------------------------------------------------------------
    use rank
    use const, only: RTYPE

   implicit none
!-------------------------------------------------------------------------------
!
!     this code is a revised vertical diffusion package ("ysupbl")
!     with a nonlocal turbulent mixing in the pbl after "mrfpbl".
!     the ysupbl (hong et al. 2006) is based on the study of noh
!     et al.(2003) and accumulated realism of the behavior of the
!     troen and mahrt (1986) concept implemented by hong and pan(1996).
!     the major ingredient of the ysupbl is the inclusion of an explicit
!     treatment of the entrainment processes at the entrainment layer.
!     this routine uses an implicit approach for vertical flux
!     divergence and does not require "miter" timesteps.
!     it includes vertical diffusion in the stable atmosphere
!     and moist vertical diffusion in clouds.
!
!     mrfpbl:
!     coded by song-you hong (ncep), implemented by jimy dudhia (ncar)
!              fall 1996
!
!     ysupbl:
!     coded by song-you hong (yonsei university) and implemented by
!              song-you hong (yonsei university) and jimy dudhia (ncar)
!              summer 2002
!
!     further modifications :
!              an enhanced stable layer mixing, april 2008
!              ==> increase pbl height when sfc is stable (hong 2010)
!              pressure-level diffusion, april 2009
!              ==> negligible differences
!              implicit forcing for momentum with clean up, july 2009
!              ==> prevents model blowup when sfc layer is too low
!              incresea of lamda, maximum (30, 0.1 x del z) feb 2010
!              ==> prevents model blowup when delz is extremely large
!              revised prandtl number at surface, peggy lemone, feb 2010
!              ==> increase kh, decrease mixing due to counter-gradient term
!              revised thermal, shin et al. mon. wea. rev. , songyou hong, aug 2011
!              ==> reduce the thermal strength when z1 < 0.1 h
!              revised prandtl number for free convection, dudhia, mar 2012
!              ==> pr0 = 1 + bke (=0.272) when neutral, kh is reduced
!              minimum kzo = 0.01, lo = min (30m,delz), hong, mar 2012
!              ==> weaker mixing when stable, and les resolution in vertical
!              gz1oz0 is removed, and psim psih are ln(z1/z0)-psim,h, hong, mar 2012
!              ==> consider thermal z0 when differs from mechanical z0
!              a bug fix in wscale computation in stable bl, sukanta basu, jun 2012
!              ==> wscale becomes small with height, and less mixing in stable bl
!              revision in background diffusion (kzo), jan 2016
!              ==> kzo = 0.1 for momentum and = 0.01 for mass to account for
!                  internal wave mixing of large et al. (1994), songyou hong, feb 2016
!              ==> alleviate superious excessive mixing when delz is large
!
!     references:
!
!        hong (2010) quart. j. roy. met. soc
!        hong, noh, and dudhia (2006), mon. wea. rev.
!        hong and pan (1996), mon. wea. rev.
!        noh, chun, hong, and raasch (2003), boundary layer met.
!        troen and mahrt (1986), boundary layer met.
!
!-------------------------------------------------------------------------------
!
   real,parameter    ::  xkzminm = 0.1,xkzminh = 0.01
   real,parameter    ::  xkzmin = 0.01,xkzmax = 1000.,rimin = -100.
   real,parameter    ::  rlam = 30.,prmin = 0.25,prmax = 4.
!   real,parameter    ::  brcr_ub = 0.0,brcr_sb = 0.25,cori = 1.e-4
   real,parameter    ::  brcr_ub = 0.25,brcr_sb = 0.25,cori = 1.e-4    !xb118
   real,parameter    ::  afac = 6.8,bfac = 6.8,pfac = 2.0,pfac_q = 2.0
!   real,parameter    ::  afac = 6.8,bfac = 6.8,pfac = 3.5,pfac_q = 3.4  !xb118
   real,parameter    ::  phifac = 8.,sfcfrac = 0.1
!   real,parameter    ::  phifac = 7.,sfcfrac = 0.1     !xb118
   real,parameter    ::  d1 = 0.02, d2 = 0.05, d3 = 0.001
   real,parameter    ::  h1 = 0.33333335, h2 = 0.6666667
   real,parameter    ::  zfmin = 1.e-8,aphi5 = 5.,aphi16 = 16.
   real,parameter    ::  tmin=1.e-2
   real,parameter    ::  gamcrt = 3.,gamcrq = 2.e-3
   real,parameter    ::  xka = 2.4e-5
   integer,parameter ::  imvdif = 1
!
!   integer,  intent(in   )   ::  ids,ide, jds,jde, kds,kde,            &
!                                 ims,ime, jms,jme, kms,kme,            &
!                                 its,ite, jts,jte, kts,kte,            &
!                                 j,ndiff
!===
!in
   integer           ::     im,ix,km,ndiff,jj

   integer,parameter ::     ysu_topdown_pblmix = 1
!   integer,parameter ::     ysu_topdown_pblmix = 0
!
   real              ::     dt,cp,g,xlv,rd
!
   real              ::     rovcp,rovg,ep1,ep2,karman,rv
!
   real,     dimension( ix, km )       ::                dz8w2d,prslk
   real,     dimension( ix, km+1 )     ::                        phii
   real,     dimension( ix, km )       ::                        phil
!
   real,     dimension( ix, km,ndiff ) ::                         qxt
   real,     dimension( im, km*3 )     ::                          qx
!
   real,     dimension( ix, km+1 )     ::                        p2di
!
   real,     dimension( ix, km )       ::                         p2d
!
   real,     dimension( ix, km )       ::        ux,vx,tx,swh,hlw,del
   real,     dimension( im )           ::                     u10,v10
   real,     dimension( im, km )       ::                        pi2d
   real,     dimension( im )           ::                      rcl,br
   real,     dimension( im )           ::                   psim,psih
   real(kind=RTYPE), dimension( im )   ::                         psk
   real,     dimension( im )           ::     xland,hfx,qfx,heat,evap
   integer,  dimension( im )           ::  ls
   real,     dimension( im )           ::                      psfcpa
   real,     dimension( im )           ::                hpbl,znt,xmu
   real,     dimension( im )           ::                 stress,wspd
! inout
   real,     dimension( im )           ::                         ust
   real,     dimension( ix, km )       ::              utnp,vtnp,ttnp
   real,     dimension( ix, km*ndiff ) ::                        qtnp
!
! out
!   integer,  dimension( im )           ::                      kpbl1d
!===
   real,     dimension( im )            ::                ctopo,ctopo2
!   real,     dimension( im )                                         , &
!!             optional                                               , &
!             intent(in   )   ::                         ctopo,ctopo2
!   real,     dimension( im )                                         , &
!             optional                                                , &
!             intent(inout)   ::                               regime
!
! local vars
!
   real,     dimension( im )       ::                              hol
   real,     dimension( im, km+1 ) ::                               zq
!
   real,     dimension( im, km )   ::                   thx,thvx,thlix, &
                                                dza,dzq,xkzom,xkzoh,za
!
   real,    dimension( im )        ::   rhox,govrth,zl1,thermal,wscale, &
                                       hgamt,hgamq,brdn,brup,phim,phih, &
                                         dusfc,dvsfc,dtsfc,dqsfc,prpbl, &
                                                       wspd1,thermalli
!
   real,    dimension( im, km )    ::      xkzm,xkzh,f1,f2,r1,r2,ad,au, &
                                          cu,al,xkzq,zfac,rhox2,hgamt2
!
!===
! inout
! in
   real,    dimension( im )           ::                       uox,vox
!===
   real,    dimension( im )           ::     brcr,sflux,zol1,brcr_sbro
!
   real,    dimension( im, km, ndiff) ::                         r3,f3
   integer, dimension( im )           ::                  kpbl,kpblold
!
   logical, dimension( im )           :: pblflg,sfcflg,stable,cloudflg

   logical ::  definebrup
!
   integer ::  n,i,k,l,ic,is,kk
   integer ::  klpbl, ktrace1, ktrace2, ktrace3
!
!
   real    ::  dt2,rdt,spdk2,fm,fh,hol1,gamfac,vpert,prnum,prnum0
   real    ::  ss,ri,qmean,tmean,alph,chi,zk,rl2,dk,sri
   real    ::  brint,dtodsd,dtodsu,rdz,dsdzt,dsdzq,dsdz2,rlamdz
   real    ::  utend,vtend,ttend,qtend
   real    ::  dtstep,govrthv
   real    ::  cont, conq, conw, conwrc
!

   real, dimension( im, km )     ::                   wscalek,wscalek2
   real, dimension( im )         ::                        wstar,delta
   real, dimension( im, km )     ::         xkzml,xkzhl,zfacent,entfac
   real, dimension( im )         ::               ust3,wstar3,wstar3_2, &
                                                           hgamu,hgamv, &
                                                               wm2, we, &
                                                  bfxpbl,hfxpbl,qfxpbl, &
                                                   ufxpbl,vfxpbl,dthvx
   real  :: prnumfac,bfx0,hfx0,qfx0,delb,dux,dvx,                       &
            dsdzu,dsdzv,wm3,dthx,dqx,wspd10,ross,tem1,dsig,tvcon,conpr, &
            prfac,prfac2,phim8z,radsum,tmp1,templ,rvls,temps,ent_eff,   &
            rcldb,bruptmp,radflux,vconvlim,vconvnew,fluxc,vconvc,vconv
!topo-corr
   real,    dimension( im, km )    ::                             fric, &
                                     tke_ysu,el_ysu,shear_ysu,buoy_ysu
   real,    dimension( im )        ::                 pblh_ysu,vconvfx
!===xb118
   real     :: tem,ho,sfx
!
!------------------------------------------------------------------------
      if (ix .lt. im) stop
!
!for topo_wind=0
  ctopo=1.
  ctopo2=1.

   klpbl = km
!=== xb118
   rovcp=rd/cp
   rovg=rd/g
!   karman=0.4
   karman=0.38
   rv=4.615e+2
   ep1=(rv/rd)-1.
   ep2=rd/rv
   uox=0.0
   vox=0.0
   do i=1,im
   if (ls(i).eq.0) then      ! ls=0:sea ; ls=1:land ; ls=2:sea ice
       xland(i)=2.            ! xland=2:sea & sea ice ; xland=1:land
   else
       xland(i)=real(ls(i))
   endif
   enddo

!     p2di(i,km+1)=prsi(i,km+1)*1000.
!  p2di=prsi*1000.
!  p2d=prsl*1000.
!  del=deli*1000.
!===
   cont=cp/g
   conq=xlv/g
   conw=1./g
   conwrc = conw*sqrt(rcl(1))
   conpr = bfac*karman*sfcfrac
!
!  k-start index for tracer diffusion
!
   ktrace1 = 0
   ktrace2 = 0 + km
   ktrace3 = 0 + km*2

!=== xb118

    if (ix .lt. im) stop

    qx=0.0
   if (ndiff .eq. 3)then
     do k=1,km
      qx(:,k+ktrace1)=qxt(:,k,1)
      qx(:,k+ktrace2)=qxt(:,k,2)
      qx(:,k+ktrace3)=qxt(:,k,3)
     enddo
   elseif (ndiff .eq. 2)then 
    do k=1,km
     do i=1,im
     qx(i,k)=qxt(i,k,1) 
     qx(i,k+ktrace2)=qxt(i,k,ndiff) 
     enddo
    enddo
   else
    print*, "warning !!! Q type must >= 2."
    stop
   endif
!===
!
   do k = 1,km
     do i = 1,im
       pi2d(i,k)=prslk(i,k)/psk(i)
       thx(i,k) = tx(i,k)/pi2d(i,k)
       thlix(i,k) = (tx(i,k)-xlv*qx(i,ktrace2+k)/cp-2.834E6*qx(i,ktrace3+k)/cp)/pi2d(i,k)
     enddo
   enddo
!
   do k = 1,km
     do i = 1,im
       tvcon = (1.+ep1*qx(i,k))
       thvx(i,k) = thx(i,k)*tvcon
     enddo
   enddo
!
   do i = 1,im
     tvcon = (1.+ep1*qx(i,1))
     rhox(i) = (psfcpa(i))/(rd*tx(i,1)*tvcon)
     govrth(i) = g/thx(i,1)
   enddo
!
!-----compute the height of full- and half-sigma levels above ground
!     level, and the layer thicknesses.
!
!   do i = 1,im
!     zq(i,1) = 0.
!   enddo
!
   do k = 1,km
     do i = 1,im
!=== xb118
!       zq(i,k+1) = dz8w2d(i,k)+zq(i,k)
!       zq(i,k+1) = dz8w2d(i,k)*conw
       zq(i,k) = phii(i,k)*conw
       za(i,k)= phil(i,k)*conw
!===
       tvcon = (1.+ep1*qx(i,k))
       rhox2(i,k) = p2d(i,k)/(rd*tx(i,k)*tvcon)
     enddo
   enddo
   do i = 1,im
       zq(i,km+1) = phii(i,km+1)*conw
   enddo
!
   do k = 1,km
     do i = 1,im
!       za(i,k) = 0.5*(zq(i,k)+zq(i,k+1))
       dzq(i,k) = zq(i,k+1)-zq(i,k)
!       del(i,k) = p2di(i,k)-p2di(i,k+1)
     enddo
   enddo
!
   do i = 1,im
     dza(i,1) = za(i,1)
   enddo
!
   do k = 1+1,km
     do i = 1,im
       dza(i,k) = za(i,k)-za(i,k-1)
     enddo
   enddo
!
!
!-----initialize vertical tendencies and
!
   utnp(:,:) = 0.
   vtnp(:,:) = 0.
   ttnp(:,:) = 0.
   qtnp(:,:) = 0.
!
   do i = 1,im
     wspd1(i) = sqrt( (ux(i,1)-uox(i))*(ux(i,1)-uox(i)) + (vx(i,1)-vox(i))*(vx(i,1)-vox(i)) )+1.e-9
   enddo
!
!---- compute vertical diffusion
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     compute preliminary variables
!

   dtstep = dt
!   dt2 = 2.*dtstep
   dt2 = 1.*dtstep
   rdt = 1./dt2
!
   do i = 1,im
     bfxpbl(i) = 0.0
     hfxpbl(i) = 0.0
     qfxpbl(i) = 0.0
     ufxpbl(i) = 0.0
     vfxpbl(i) = 0.0
     hgamu(i)  = 0.0
     hgamv(i)  = 0.0
     delta(i)  = 0.0
     wstar3_2(i) =  0.0
     wstar3(i) =  0.0
     wstar(i) =  0.0
   enddo
!
   do k = 1,klpbl
     do i = 1,im
       xkzh(i,k) = 0.0
       xkzm(i,k) = 0.0
       xkzq(i,k) = 0.0
       wscalek(i,k) = 0.0
       wscalek2(i,k) = 0.0
     enddo
   enddo
!
   do k = 1,klpbl
     do i = 1,im
       zfac(i,k) = 0.0
     enddo
   enddo
   do k = 1,klpbl-1
     do i = 1,im
       xkzom(i,k) = xkzminm
       xkzoh(i,k) = xkzminh
     enddo
   enddo
!
   do i = 1,im
     dusfc(i) = 0.
     dvsfc(i) = 0.
     dtsfc(i) = 0.
     dqsfc(i) = 0.
   enddo
!
   do i = 1,im
     hgamt(i)  = 0.
     hgamq(i)  = 0.
     wscale(i) = 0.
     kpbl(i)   = 1
     hpbl(i)   = zq(i,1)
     zl1(i)    = za(i,1)
     thermal(i)= thvx(i,1)
     thermalli(i) = thlix(i,1)
     pblflg(i) = .true.
     sfcflg(i) = .true.
     hfx(i)=heat(i)*rhox(i)*cp     ! heat[m*K/s] * rho[kg/m^3] * cp[J/kg/K] = [W/m^2]
     qfx(i)=evap(i)*rhox(i)        ! evap[m/s] * rho[kg/m^3] = [kg/m^2/s]
     sflux(i) = hfx(i)/rhox(i)/cp + qfx(i)/rhox(i)*ep1*thx(i,1)   ![m*K/s]
!     if(br(i).gt.0.0) sfcflg(i) = .false.
     if(br(i).gt.0.25) sfcflg(i) = .false.
   enddo
!
!     compute the first guess of pbl height
!
   do i = 1,im
     stable(i) = .false.
     brup(i) = br(i)
     brcr(i) = brcr_ub
   enddo
!
   do k = 2,klpbl
     do i = 1,im
       if(.not.stable(i))then
         brdn(i) = brup(i)
         spdk2   = max(ux(i,k)**2+vx(i,k)**2,1.)
         brup(i) = (thvx(i,k)-thermal(i))*(g*za(i,k)/thvx(i,1))/spdk2
         kpbl(i) = k
         stable(i) = brup(i).gt.brcr(i)
       endif
     enddo
   enddo
!
   do i = 1,im
     k = kpbl(i)
     if(brdn(i).ge.brcr(i))then
       brint = 0.
     elseif(brup(i).le.brcr(i))then
       brint = 1.
     else
       brint = (brcr(i)-brdn(i))/(brup(i)-brdn(i))
     endif
     hpbl(i) = za(i,k-1)+brint*(za(i,k)-za(i,k-1))
     if(hpbl(i).lt.zq(i,2)) kpbl(i) = 1
     if(kpbl(i).le.1) pblflg(i) = .false.
   enddo
!
   do i = 1,im
     fm = psim(i)
     fh = psih(i)
     zol1(i) = max(br(i)*fm*fm/fh,rimin)
     if(sfcflg(i))then
       zol1(i) = min(zol1(i),-zfmin)
     else
       zol1(i) = max(zol1(i),zfmin)
     endif
     hol1 = zol1(i)*hpbl(i)/zl1(i)*sfcfrac
     if(sfcflg(i))then
       phim(i) = (1.-aphi16*hol1)**(-1./4.)
       phih(i) = (1.-aphi16*hol1)**(-1./2.)
       bfx0 = max(sflux(i),0.)
       hfx0 = max(hfx(i)/rhox(i)/cp,0.)
       qfx0 = max(ep1*thx(i,1)*qfx(i)/rhox(i),0.)
       wstar3(i) = (govrth(i)*bfx0*hpbl(i))
       wstar(i) = (wstar3(i))**h1
     else
       phim(i) = (1.+aphi5*hol1)
       phih(i) = phim(i)
       wstar(i)  = 0.
       wstar3(i) = 0.
     endif
     ust(i) = sqrt(stress(i))
     ust3(i)   = ust(i)**3.
     wscale(i) = (ust3(i)+phifac*karman*wstar3(i)*0.5)**h1
     wscale(i) = min(wscale(i),ust(i)*aphi16)
     wscale(i) = max(wscale(i),ust(i)/aphi5)
   enddo
!
!     compute the surface variables for pbl height estimation
!     under unstable conditions
!
   do i = 1,im
     if(sfcflg(i).and.sflux(i).gt.0.0)then
       gamfac   = bfac/rhox(i)/wscale(i)
       hgamt(i) = min(gamfac*hfx(i)/cp,gamcrt)
       hgamq(i) = min(gamfac*qfx(i),gamcrq)
       vpert = (hgamt(i)+ep1*thx(i,1)*hgamq(i))/bfac*afac
       thermal(i) = thermal(i)+max(vpert,0.)*min(za(i,1)/(sfcfrac*hpbl(i)),1.0)
       thermalli(i)= thermalli(i)+max(vpert,0.)*min(za(i,1)/(sfcfrac*hpbl(i)),1.0)
       hgamt(i) = max(hgamt(i),0.0)
       hgamq(i) = max(hgamq(i),0.0)
       brint    = -15.9*ust(i)*ust(i)/wspd(i)*wstar3(i)/(wscale(i)**4.)
!       brint    = -8.*ust(i)*ust(i)/wspd(i)*wstar3(i)/(wscale(i)**4.)  !xb118
       hgamu(i) = brint*ux(i,1)
       hgamv(i) = brint*vx(i,1)
     else
       pblflg(i) = .false.
     endif
   enddo
!     enhance the pbl height by considering the thermal
!
   do i = 1,im
     if(pblflg(i))then
       kpbl(i) = 1
       hpbl(i) = zq(i,1)
     endif
   enddo
!
   do i = 1,im
     if(pblflg(i))then
       stable(i) = .false.
       brup(i) = br(i)
       brcr(i) = brcr_ub
     endif
   enddo
!
   do k = 2,klpbl
     do i = 1,im
       if(.not.stable(i).and.pblflg(i))then
         brdn(i) = brup(i)
         spdk2   = max(ux(i,k)**2+vx(i,k)**2,1.)
         brup(i) = (thvx(i,k)-thermal(i))*(g*za(i,k)/thvx(i,1))/spdk2
         kpbl(i) = k
         stable(i) = brup(i).gt.brcr(i)
       endif
     enddo
   enddo

!
!     enhance pbl by theta-li
!
   if (ysu_topdown_pblmix.eq.1)then
!   if (ysu_topdown_pblmix.eq.2)then
     do i = 1,im
        kpblold(i) = kpbl(i)
        definebrup=.false.
        do k = kpblold(i), km-1
           spdk2   = max(ux(i,k)**2+vx(i,k)**2,1.)
           bruptmp = (thlix(i,k)-thermalli(i))*(g*za(i,k)/thlix(i,1))/spdk2
           stable(i) = bruptmp.ge.brcr(i)
           if (definebrup) then
           kpbl(i) = k
           brup(i) = bruptmp
           definebrup=.false.
           endif
           if (.not.stable(i)) then !overwrim brup brdn values
           brdn(i)=bruptmp
           definebrup=.true.
           pblflg(i)=.true.
           endif
        enddo
     enddo
   endif

   do i = 1,im
     if(pblflg(i)) then
       k = kpbl(i)
       if(brdn(i).ge.brcr(i))then
         brint = 0.
       elseif(brup(i).le.brcr(i))then
         brint = 1.
       else
         brint = (brcr(i)-brdn(i))/(brup(i)-brdn(i))
       endif
       hpbl(i) = za(i,k-1)+brint*(za(i,k)-za(i,k-1))
       if(hpbl(i).lt.zq(i,2)) kpbl(i) = 1
       if(kpbl(i).le.1) pblflg(i) = .false.
     endif
   enddo
!
!     stable boundary layer
!
   do i = 1,im
     if((.not.sfcflg(i)).and.hpbl(i).lt.zq(i,2)) then
       brup(i) = br(i)
       stable(i) = .false.
     else
       stable(i) = .true.
     endif
   enddo
!
   do i = 1,im
     if((.not.stable(i)).and.((xland(i)-1.5).ge.0))then
       wspd10 = u10(i)*u10(i) + v10(i)*v10(i)
       wspd10 = sqrt(wspd10)
       ross = wspd10 / (cori*znt(i))
       brcr_sbro(i) = min(0.16*(1.e-7*ross)**(-0.18),.3)
     endif
   enddo
!
   do i = 1,im
     if(.not.stable(i))then
       if((xland(i)-1.5).ge.0)then
         brcr(i) = brcr_sbro(i)
       else
         brcr(i) = brcr_sb
       endif
     endif
   enddo
!
   do k = 2,klpbl
     do i = 1,im
       if(.not.stable(i))then
         brdn(i) = brup(i)
         spdk2   = max(ux(i,k)**2+vx(i,k)**2,1.)
         brup(i) = (thvx(i,k)-thermal(i))*(g*za(i,k)/thvx(i,1))/spdk2
         kpbl(i) = k
         stable(i) = brup(i).gt.brcr(i)
       endif
     enddo
   enddo
!
   do i = 1,im
     if((.not.sfcflg(i)).and.hpbl(i).lt.zq(i,2)) then
       k = kpbl(i)
       if(brdn(i).ge.brcr(i))then
         brint = 0.
       elseif(brup(i).le.brcr(i))then
         brint = 1.
       else
         brint = (brcr(i)-brdn(i))/(brup(i)-brdn(i))
       endif
       hpbl(i) = za(i,k-1)+brint*(za(i,k)-za(i,k-1))
       if(hpbl(i).lt.zq(i,2)) kpbl(i) = 1
       if(kpbl(i).le.1) pblflg(i) = .false.
     endif
   enddo

!
!     estimate the entrainment parameters
!
   do i = 1,im
   radsum=0.
     cloudflg(i)=.false. 
     if(pblflg(i)) then
       k = kpbl(i) - 1
       wm3       = wstar3(i) + 5. * ust3(i)
       wm2(i)    = wm3**h2
       bfxpbl(i) = -0.15*thvx(i,1)/g*wm3/hpbl(i)
       dthvx(i)  = max(thvx(i,k+1)-thvx(i,k),tmin)
       we(i) = max(bfxpbl(i)/dthvx(i),-sqrt(wm2(i)))
       if((qx(i,ktrace2+k)+qx(i,ktrace3+k)).gt.0.01e-3.and.ysu_topdown_pblmix.eq.1)then
           if ( kpbl(i) .ge. 2) then
                cloudflg(i)=.true. 
                templ=thlix(i,k)*p2di(i,k+1)**rovcp
                !rvls is ws at full level
                rvls=100.*6.112*EXP(17.67*(templ-273.16)/(templ-29.65))*(ep2/p2di(i,k+1))
                temps=templ + ((qx(i,k)+qx(i,ktrace2+k))-rvls)/(cp/xlv  + &
                ep2*xlv*rvls/(rd*templ**2))
                rvls=100.*6.112*EXP(17.67*(temps-273.15)/(temps-29.65))*(ep2/p2di(i,k+1))
                rcldb=max((qx(i,k)+qx(i,ktrace2+k))-rvls,0.)
                !entrainment efficiency
                dthvx(i)  = (thlix(i,k+2)+thx(i,k+2)*ep1*(qx(i,k+2)+qx(i,ktrace2+k+2))) &
                          - (thlix(i,k) + thx(i,k)  *ep1*(qx(i,k)  +qx(i,ktrace2+k)))
                dthvx(i)  = max(dthvx(i),0.1)
                tmp1      = xlv/cp * rcldb/(pi2d(i,k)*dthvx(i))
                ent_eff   = 0.2 * 8. * tmp1 +0.2

                radsum=0.
                do kk = 1,kpbl(i)-1
!=== xb118
!                   radflux=rthraten(i,kk)*pi2d(i,kk) !converts theta/s to temp/s
!                   radflux=dzq(i,kk)*(swh(i,kk)*xmu(i)+hlw(i,kk)) ! SW+LW [temp/s]
                   radflux=swh(i,kk)*xmu(i)+hlw(i,kk)
!                   radflux=radflux*cp/g*(p2diORG(i,kk)-p2diORG(i,kk+1)) ! converts temp/s to W/m^2
                   radflux=radflux*cp/g*(p2di(i,kk)-p2di(i,kk+1)) ! converts temp/s to W/m^2
!===
                   if (radflux < 0.0 ) radsum=abs(radflux)+radsum
                enddo
                radsum=max(radsum,0.0)

                !recompute entrainment from sfc thermals
!                bfx0 = max(max(sflux(i),0.0)-radsum/rhox2(i,k)/cp,0.)
                bfx0 = max(sflux(i),0.0)
                wm3 = (govrth(i)*bfx0*hpbl(i))+5. * ust3(i)
                wm2(i)    = wm3**h2
                bfxpbl(i) = -0.15*thvx(i,1)/g*wm3/hpbl(i)
                dthvx(i)  = max(thvx(i,k+1)-thvx(i,k),0.1)
!                dthvx(i)  = max(thvx(i,k+1)-thvx(i,k),tmin)
                we(i) = max(bfxpbl(i)/dthvx(i),-sqrt(wm2(i)))
                 if (.not.sfcflg(i))  we(i) =0.0

                !entrainment from PBL top thermals
!                bfx0 = max(radsum/rhox2(i,k)/cp-max(sflux(i),0.0),0.)
                bfx0 = max(radsum/rhox2(i,k)/cp,0.)
                wm3       = (g/thvx(i,k)*bfx0*hpbl(i)) ! this is wstar3(i)
                wm2(i)    = wm2(i)+wm3**h2
                bfxpbl(i) = - ent_eff * bfx0
                dthvx(i)  = max(thvx(i,k+1)-thvx(i,k),0.1)
                we(i) = we(i) + max(bfxpbl(i)/dthvx(i),-sqrt(wm3**h2))

!                !wstar3_2
                bfx0 = max(radsum/rhox2(i,k)/cp,0.)
                wstar3_2(i) =  (g/thvx(i,k)*bfx0*hpbl(i))
!                !recompute hgamt 
                wscale(i) = (ust3(i)+phifac*karman*(wstar3(i)+wstar3_2(i))*0.5)**h1
                wscale(i) = min(wscale(i),ust(i)*aphi16)
                wscale(i) = max(wscale(i),ust(i)/aphi5)
                gamfac   = bfac/rhox(i)/wscale(i)
                hgamt(i) = min(gamfac*hfx(i)/cp,gamcrt)
                hgamq(i) = min(gamfac*qfx(i),gamcrq)
                gamfac   = bfac/rhox2(i,k)/wscale(i)
                hgamt2(i,k) = min(gamfac*radsum/cp,gamcrt)
                hgamt(i) = max(hgamt(i),0.0) + max(hgamt2(i,k),0.0)
                brint    = -15.9*ust(i)*ust(i)/wspd(i)*(wstar3(i)+wstar3_2(i))/(wscale(i)**4.)
!                brint    = -8.*ust(i)*ust(i)/wspd(i)*(wstar3(i)+wstar3_2(i))/(wscale(i)**4.) !xb118
                hgamu(i) = brint*ux(i,1)
                hgamv(i) = brint*vx(i,1)
           endif
       endif
       prpbl(i) = 1.0
       dthx  = max(thx(i,k+1)-thx(i,k),tmin)
       dqx   = min(qx(i,k+1)-qx(i,k),0.0)
       hfxpbl(i) = we(i)*dthx
       qfxpbl(i) = we(i)*dqx
!
       dux = ux(i,k+1)-ux(i,k)
       dvx = vx(i,k+1)-vx(i,k)
       if(dux.gt.tmin) then
         ufxpbl(i) = max(prpbl(i)*we(i)*dux,-ust(i)*ust(i))
       elseif(dux.lt.-tmin) then
         ufxpbl(i) = min(prpbl(i)*we(i)*dux,ust(i)*ust(i))
       else
         ufxpbl(i) = 0.0
       endif
       if(dvx.gt.tmin) then
         vfxpbl(i) = max(prpbl(i)*we(i)*dvx,-ust(i)*ust(i))
       elseif(dvx.lt.-tmin) then
         vfxpbl(i) = min(prpbl(i)*we(i)*dvx,ust(i)*ust(i))
       else
         vfxpbl(i) = 0.0
       endif
       delb  = govrth(i)*d3*hpbl(i)
       delta(i) = min(d1*hpbl(i) + d2*wm2(i)/delb,100.)
     endif
   enddo
!
   do k = 1,klpbl
     do i = 1,im
       if(pblflg(i).and.k.ge.kpbl(i))then
         entfac(i,k) = ((zq(i,k+1)-hpbl(i))/delta(i))**2.
       else
         entfac(i,k) = 1.e30
       endif
     enddo
   enddo
!
!     compute diffusion coefficients below pbl
!
   do k = 1,klpbl
     do i = 1,im
       if(k.lt.kpbl(i)) then
         zfac(i,k) = min(max((1.-(zq(i,k+1)-zl1(i))/(hpbl(i)-zl1(i))),zfmin),1.)
         zfacent(i,k) = (1.-zfac(i,k))**3.
!         zfacent(i,k) = (1.-zfac(i,k))**4.
         wscalek(i,k) = (ust3(i)+phifac*karman*wstar3(i)*(1.-zfac(i,k)))**h1
         wscalek2(i,k) = (phifac*karman*wstar3_2(i)*(zfac(i,k)))**h1
         if(sfcflg(i)) then
           prfac = conpr
           prfac2 = 15.9*(wstar3(i)+wstar3_2(i))/ust3(i)/(1.+4.*karman*(wstar3(i)+wstar3_2(i))/ust3(i))
!           prfac2 = 8.*(wstar3(i)+wstar3_2(i))/ust3(i)/(1.+4.*karman*(wstar3(i)+wstar3_2(i))/ust3(i)) !xb118
           prnumfac = -3.*(max(zq(i,k+1)-sfcfrac*hpbl(i),0.))**2./hpbl(i)**2.
         else
           prfac = 0.
           prfac2 = 0.
           prnumfac = 0.
           phim8z = 1.+aphi5*zol1(i)*zq(i,k+1)/zl1(i)
           wscalek(i,k) = ust(i)/phim8z
           wscalek(i,k) = max(wscalek(i,k),0.001)
         endif
         prnum0 = (phih(i)/phim(i)+prfac)
         prnum0 = max(min(prnum0,prmax),prmin)
           xkzm(i,k) = wscalek(i,k) *karman*    zq(i,k+1)      *    zfac(i,k)**pfac+ &
                       wscalek2(i,k)*karman*(hpbl(i)-zq(i,k+1))*(1-zfac(i,k))**pfac
         !Do not include xkzm at kpbl-1 since it changes entrainment
         if (k.eq.kpbl(i)-1.and.cloudflg(i).and.we(i).lt.0.0) then
           xkzm(i,k) = 0.0
         endif
         prnum =  1. + (prnum0-1.)*exp(prnumfac)
         xkzq(i,k) = xkzm(i,k)/prnum*zfac(i,k)**(pfac_q-pfac)
         prnum0 = prnum0/(1.+prfac2*karman*sfcfrac)
         prnum =  1. + (prnum0-1.)*exp(prnumfac)
         xkzh(i,k) = xkzm(i,k)/prnum
         xkzm(i,k) = xkzm(i,k)+xkzom(i,k)
         xkzh(i,k) = xkzh(i,k)+xkzoh(i,k)
         xkzq(i,k) = xkzq(i,k)+xkzoh(i,k)
         xkzm(i,k) = min(xkzm(i,k),xkzmax)
         xkzh(i,k) = min(xkzh(i,k),xkzmax)
         xkzq(i,k) = min(xkzq(i,k),xkzmax)
       endif
     enddo
   enddo

!
!     compute diffusion coefficients over pbl (free atmosphere)
!
   do k = 1,km-1
     do i = 1,im
       if(k.ge.kpbl(i)) then
         ss = ((ux(i,k+1)-ux(i,k))*(ux(i,k+1)-ux(i,k))                         &
              +(vx(i,k+1)-vx(i,k))*(vx(i,k+1)-vx(i,k)))                        &
              /(dza(i,k+1)*dza(i,k+1))+1.e-9
         govrthv = g/(0.5*(thvx(i,k+1)+thvx(i,k)))
         ri = govrthv*(thvx(i,k+1)-thvx(i,k))/(ss*dza(i,k+1))
!         if(imvdif.eq.1.and.ndiff.ge.3)then
         if(imvdif.eq.1.and.ndiff.ge.2)then
           if((qx(i,ktrace2+k)+qx(i,ktrace3+k)).gt.0.01e-3.and.(qx(i           &
             ,ktrace2+k+1)+qx(i,ktrace3+k+1)).gt.0.01e-3)then
!      in cloud
             qmean = 0.5*(qx(i,k)+qx(i,k+1))
             tmean = 0.5*(tx(i,k)+tx(i,k+1))
             alph  = xlv*qmean/rd/tmean
             chi   = xlv*xlv*qmean/cp/rv/tmean/tmean
             ri    = (1.+alph)*(ri-g*g/ss/tmean/cp*((chi-alph)/(1.+chi)))
           endif
         endif
         zk = karman*zq(i,k+1)
         rlamdz = min(max(0.1*dza(i,k+1),rlam),300.)
         rlamdz = min(dza(i,k+1),rlamdz)
         rl2 = (zk*rlamdz/(rlamdz+zk))**2
         dk = rl2*sqrt(ss)
         if(ri.lt.0.)then
! unstable regime
           ri = max(ri, rimin)
           sri = sqrt(-ri)
           xkzm(i,k) = dk*(1+8.*(-ri)/(1+1.746*sri))
           xkzh(i,k) = dk*(1+8.*(-ri)/(1+1.286*sri))
         else
! stable regime
           xkzh(i,k) = dk/(1+5.*ri)**2
           prnum = 1.0+2.1*ri
           prnum = min(prnum,prmax)
           xkzm(i,k) = xkzh(i,k)*prnum
         endif
!
         xkzm(i,k) = xkzm(i,k)+xkzom(i,k)
         xkzh(i,k) = xkzh(i,k)+xkzoh(i,k)
         xkzm(i,k) = min(xkzm(i,k),xkzmax)
         xkzh(i,k) = min(xkzh(i,k),xkzmax)
         xkzml(i,k) = xkzm(i,k)
         xkzhl(i,k) = xkzh(i,k)
       endif
     enddo
   enddo
!
!     compute tridiagonal matrix elements for heat
!
   do k = 1,km
     do i = 1,im
       au(i,k) = 0.
       al(i,k) = 0.
       ad(i,k) = 0.
       f1(i,k) = 0.
     enddo
   enddo
!
   do i = 1,im
     ad(i,1) = 1.
     f1(i,1) = thx(i,1)-300.+hfx(i)/cont/del(i,1)*dt2
   enddo
!
   do k = 1,km-1
     do i = 1,im
       dtodsd = dt2/del(i,k)
       dtodsu = dt2/del(i,k+1)
       dsig   = p2d(i,k)-p2d(i,k+1)
       rdz    = 1./dza(i,k+1)
       tem1   = dsig*xkzh(i,k)*rdz
       if(pblflg(i).and.k.lt.kpbl(i)) then
         dsdzt = tem1*(-hgamt(i)/hpbl(i)-hfxpbl(i)*zfacent(i,k)/xkzh(i,k))
         f1(i,k)   = f1(i,k)+dtodsd*dsdzt
         f1(i,k+1) = thx(i,k+1)-300.-dtodsu*dsdzt
       elseif(pblflg(i).and.k.ge.kpbl(i).and.entfac(i,k).lt.4.6) then
         xkzh(i,k) = -we(i)*dza(i,kpbl(i))*exp(-entfac(i,k))
         xkzh(i,k) = sqrt(xkzh(i,k)*xkzhl(i,k))
         xkzh(i,k) = max(xkzh(i,k),xkzoh(i,k))
         xkzh(i,k) = min(xkzh(i,k),xkzmax)
         f1(i,k+1) = thx(i,k+1)-300.
       else
         f1(i,k+1) = thx(i,k+1)-300.
       endif
       tem1   = dsig*xkzh(i,k)*rdz
       dsdz2     = tem1*rdz
       au(i,k)   = -dtodsd*dsdz2
       al(i,k)   = -dtodsu*dsdz2
       ad(i,k)   = ad(i,k)-au(i,k)
       ad(i,k+1) = 1.-al(i,k)
     enddo
   enddo
!
! copies here to avoid duplicate input args for tridin
!
   do k = 1,km
     do i = 1,im
       cu(i,k) = au(i,k)
       r1(i,k) = f1(i,k)
     enddo
   enddo
!
   call tridin_ysu(al,ad,cu,r1,au,f1,1,im,1,km,1)
!
!     recover tendencies of heat
!
   do k = km,1,-1
     do i = 1,im
       ttend = (f1(i,k)-thx(i,k)+300.)*rdt*pi2d(i,k)
       ttnp(i,k) = ttnp(i,k)+ttend
       dtsfc(i) = dtsfc(i)+ttend*cont*del(i,k)/pi2d(i,k)
     enddo
   enddo

     do k = 1,km
       do i = 1,im
         tx(i,k) = (f1(i,k)+300.)*pi2d(i,k)
       enddo
     enddo
!
!     compute tridiagonal matrix elements for moisture, clouds, and gases
!
   do k = 1,km
     do i = 1,im
       au(i,k) = 0.
       al(i,k) = 0.
       ad(i,k) = 0.
     enddo
   enddo
!
   do ic = 1,ndiff
     do i = 1,im
       do k = 1,km
         f3(i,k,ic) = 0.
       enddo
     enddo
   enddo
!
   do i = 1,im
     ad(i,1) = 1.
     f3(i,1,1) = qx(i,1)+qfx(i)*g/del(i,1)*dt2
   enddo
!
   if(ndiff.ge.2) then
     do ic = 2,ndiff
       is = (ic-1) * km
       do i = 1,im
         f3(i,1,ic) = qx(i,1+is)
       enddo
     enddo
   endif
!
   do k = 1,km-1
     do i = 1,im
       if(k.ge.kpbl(i)) then
         xkzq(i,k) = xkzh(i,k)
       endif
     enddo
   enddo
!
   do k = 1,km-1
     do i = 1,im
       dtodsd = dt2/del(i,k)
       dtodsu = dt2/del(i,k+1)
       dsig   = p2d(i,k)-p2d(i,k+1)
       rdz    = 1./dza(i,k+1)
       tem1   = dsig*xkzq(i,k)*rdz
       if(pblflg(i).and.k.lt.kpbl(i)) then
         dsdzq = tem1*(-qfxpbl(i)*zfacent(i,k)/xkzq(i,k))
         f3(i,k,1) = f3(i,k,1)+dtodsd*dsdzq
         f3(i,k+1,1) = qx(i,k+1)-dtodsu*dsdzq
       elseif(pblflg(i).and.k.ge.kpbl(i).and.entfac(i,k).lt.4.6) then
         xkzq(i,k) = -we(i)*dza(i,kpbl(i))*exp(-entfac(i,k))
         xkzq(i,k) = sqrt(xkzq(i,k)*xkzhl(i,k))
         xkzq(i,k) = max(xkzq(i,k),xkzoh(i,k))
         xkzq(i,k) = min(xkzq(i,k),xkzmax)
         f3(i,k+1,1) = qx(i,k+1)
       else
         f3(i,k+1,1) = qx(i,k+1)
       endif
       tem1   = dsig*xkzq(i,k)*rdz
       dsdz2     = tem1*rdz
       au(i,k)   = -dtodsd*dsdz2
       al(i,k)   = -dtodsu*dsdz2
       ad(i,k)   = ad(i,k)-au(i,k)
       ad(i,k+1) = 1.-al(i,k)
!        if (f3(i,k,1).lt.-1e-3)print*,i,k,'f3(i,k,1)=',f3(i,k,1),'xkzm=',xkzm(i,k),'xkzh=',xkzh(i,k),'xkzq=',xkzq(i,k),'sfcflg=',sfcflg(i),'wscale=',wscale(i),'wscalek=',wscalek(i,k),'wscalek2=',wscalek2(i,k),'ust3=',ust3(i),'wstar3=',wstar3(i),'zfac=',zfac(i,k),'qfxpbl=',qfxpbl(i),'we=',we(i)
     enddo
   enddo
!
   if(ndiff.ge.2) then
     do ic = 2,ndiff
       is = (ic-1) * km
       do k = 1,km-1
         do i = 1,im
           f3(i,k+1,ic) = qx(i,k+1+is)
         enddo
       enddo
     enddo
   endif
!
! copies here to avoid duplicate input args for tridin
!
   do k = 1,km
     do i = 1,im
       cu(i,k) = au(i,k)
     enddo
   enddo
!
   do ic = 1,ndiff
     do k = 1,km
       do i = 1,im
         r3(i,k,ic) = f3(i,k,ic)
       enddo
     enddo
   enddo
!
!     solve tridiagonal problem for moisture, clouds, and gases
!
   call tridin_ysu(al,ad,cu,r3,au,f3,1,im,1,km,ndiff)
!
!     recover tendencies of heat and moisture
!
   do k = km,1,-1
     do i = 1,im
       qtend = (f3(i,k,1)-qx(i,k))*rdt
       qtnp(i,k) = qtnp(i,k)+qtend
       dqsfc(i) = dqsfc(i)+qtend*conq*del(i,k)
     enddo
   enddo
!
   if(ndiff.ge.2) then
     do ic = 2,ndiff
       is = (ic-1) * km
       do k = km,1,-1
         do i = 1,im
           qtend = (f3(i,k,ic)-qx(i,k+is))*rdt
           qtnp(i,k+is) = qtnp(i,k+is)+qtend
         enddo
       enddo
     enddo
   endif

    do k = 1,km
     do i =1,im
      qxt(i,k,1)=f3(i,k,1)
     enddo
    enddo
   if (ndiff .eq. 2 ) then
    do k = 1,km
     do i =1,im
      qxt(i,k,ndiff)=f3(i,k,ndiff)
     enddo
    enddo
   endif
!
!     compute tridiagonal matrix elements for momentum
!
   do i = 1,im
     do k = 1,km
       au(i,k) = 0.
       al(i,k) = 0.
       ad(i,k) = 0.
       f1(i,k) = 0.
       f2(i,k) = 0.
     enddo
   enddo
!
!************************************************************************* xb118
! paj: ctopo=1 if topo_wind=0 (default)
!raquel---paj tke code (could be replaced with shin-hong tke in future
   do i = 1,im
      do k= 1, km-1
        shear_ysu(i,k)=xkzm(i,k)*((-hgamu(i)/hpbl(i)+(ux(i,k+1)-ux(i,k))/dza(i,k+1))*(ux(i,k+1)-ux(i,k))/dza(i,k+1) &
        + (-hgamv(i)/hpbl(i)+(vx(i,k+1)-vx(i,k))/dza(i,k+1))*(vx(i,k+1)-vx(i,k))/dza(i,k+1))
         buoy_ysu(i,k)=xkzh(i,k)*g*(1.0/thx(i,k))*(-hgamt(i)/hpbl(i)+(thx(i,k+1)-thx(i,k))/dza(i,k+1))

       zk = karman*zq(i,k+1)
 !over pbl
       if (k.ge.kpbl(i)) then
        rlamdz = min(max(0.1*dza(i,k+1),rlam),300.)
        rlamdz = min(dza(i,k+1),rlamdz)
       else
 !in pbl
        rlamdz = 150.0
       endif
       el_ysu(i,k) = zk*rlamdz/(rlamdz+zk)
       tke_ysu(i,k)=16.6*el_ysu(i,k)*(shear_ysu(i,k)-buoy_ysu(i,k))
 !q2 when q3 positive
       if(tke_ysu(i,k).le.0) then
        tke_ysu(i,k)=0.0
       else
        tke_ysu(i,k)=(tke_ysu(i,k))**0.66
       endif
      enddo
 !Hybrid pblh of MYNN
 !tke is q2
      CALL GET_PBLH(1,km,pblh_ysu(i),thvx(i,1:km),&
      &    tke_ysu(i,1:km),zq(i,1:km+1),dzq(i,1:km),xland(i))

!--- end of paj tke
! compute vconv
!      Use Beljaars over land
        if (xland(i).lt.1.5) then
        fluxc = max(sflux(i),0.0)
        vconvc=1.
        VCONV = vconvc*(g/thvx(i,1)*pblh_ysu(i)*fluxc)**.33
        else
! for water there is no topo effect so vconv not needed
        VCONV = 0.
        endif
        vconvfx(i) = vconv
!raquel
!ctopo stability correction
      fric(i,1)=ust(i)**2/wspd1(i)*rhox(i)*g/del(i,1)*dt2         &
        *(wspd1(i)/wspd(i))**2
      if(.true.) then
!      if(present(ctopo)) then
        vconvnew=0.9*vconvfx(i)+1.5*(max((pblh_ysu(i)-500)/1000.0,0.0))
        vconvlim = min(vconvnew,1.0)
        ad(i,1) = 1.+fric(i,1)*vconvlim+ctopo(i)*fric(i,1)*(1-vconvlim)
      else
       ad(i,1) = 1.+fric(i,1)
      endif
     f1(i,1) = ux(i,1)+uox(i)*ust(i)**2*rhox(i)*g/del(i,1)*dt2/wspd1(i)*(wspd1(i)/wspd(i))**2
     f2(i,1) = vx(i,1)+vox(i)*ust(i)**2*rhox(i)*g/del(i,1)*dt2/wspd1(i)*(wspd1(i)/wspd(i))**2
   enddo
!********************************************************************************
!
   do k = 1,km-1
     do i = 1,im
       dtodsd = dt2/del(i,k)
       dtodsu = dt2/del(i,k+1)
       dsig   = p2d(i,k)-p2d(i,k+1)
       rdz    = 1./dza(i,k+1)
       tem1   = dsig*xkzm(i,k)*rdz
       if(pblflg(i).and.k.lt.kpbl(i))then
         dsdzu     = tem1*(-hgamu(i)/hpbl(i)-ufxpbl(i)*zfacent(i,k)/xkzm(i,k))
         dsdzv     = tem1*(-hgamv(i)/hpbl(i)-vfxpbl(i)*zfacent(i,k)/xkzm(i,k))
         f1(i,k)   = f1(i,k)+dtodsd*dsdzu
         f1(i,k+1) = ux(i,k+1)-dtodsu*dsdzu
         f2(i,k)   = f2(i,k)+dtodsd*dsdzv
         f2(i,k+1) = vx(i,k+1)-dtodsu*dsdzv
       elseif(pblflg(i).and.k.ge.kpbl(i).and.entfac(i,k).lt.4.6) then
         xkzm(i,k) = prpbl(i)*xkzh(i,k)
         xkzm(i,k) = sqrt(xkzm(i,k)*xkzml(i,k))
         xkzm(i,k) = max(xkzm(i,k),xkzom(i,k))
         xkzm(i,k) = min(xkzm(i,k),xkzmax)
        f1(i,k+1) = ux(i,k+1)
         f2(i,k+1) = vx(i,k+1)
       else
         f1(i,k+1) = ux(i,k+1)
         f2(i,k+1) = vx(i,k+1)
       endif
       tem1   = dsig*xkzm(i,k)*rdz
       dsdz2     = tem1*rdz
       au(i,k)   = -dtodsd*dsdz2
       al(i,k)   = -dtodsu*dsdz2
       ad(i,k)   = ad(i,k)-au(i,k)
       ad(i,k+1) = 1.-al(i,k)
     enddo
   enddo
!
! copies here to avoid duplicate input args for tridin
!
   do k = 1,km
     do i = 1,im
       cu(i,k) = au(i,k)
       r1(i,k) = f1(i,k)
       r2(i,k) = f2(i,k)
     enddo
   enddo
!
!     solve tridiagonal problem for momentum
!
   call tridi1n(al,ad,cu,r1,r2,au,f1,f2,1,im,1,km,1)
!
!     recover tendencies of momentum
!
   do k = km,1,-1
     do i = 1,im
       utend = (f1(i,k)-ux(i,k))*rdt
       vtend = (f2(i,k)-vx(i,k))*rdt
       utnp(i,k) = utnp(i,k)+utend
       vtnp(i,k) = vtnp(i,k)+vtend
       dusfc(i) = dusfc(i) + utend*conwrc*del(i,k)
       dvsfc(i) = dvsfc(i) + vtend*conwrc*del(i,k)
     enddo
   enddo
!
!******************************************************************************* xb118
! paj: ctopo2=1 if topo_wind=0 (default)
!
   do i = 1,im
     if(.true.) then ! mchen for NMM
!     if(present(ctopo).and.present(ctopo2)) then ! mchen for NMM
       u10(i) = ctopo2(i)*u10(i)+(1-ctopo2(i))*ux(i,1)
       v10(i) = ctopo2(i)*v10(i)+(1-ctopo2(i))*vx(i,1)
     endif !mchen
   enddo
!********************************************************************************
!
!---- end of vertical diffusion
!

    do k = 1,km
     do i =1,im
      ux(i,k)=f1(i,k)/rcl(i)
      vx(i,k)=f2(i,k)/rcl(i)
     enddo
    enddo


!   do i = 1,im
!     kpbl1d(i) = kpbl(i)
!   enddo
!
!
   end subroutine ysu2d
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine tridi1n(cl,cm,cu,r1,r2,au,f1,f2,its,ite,kts,kte,nt)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
   integer, intent(in )      ::     its,ite, kts,kte, nt
!
   real, dimension( its:ite, kts+1:kte+1 )                                   , &
         intent(in   )  ::                                                 cl
!
   real, dimension( its:ite, kts:kte )                                       , &
         intent(in   )  ::                                                 cm, &
                                                                           r1
   real, dimension( its:ite, kts:kte,nt )                                    , &
         intent(in   )  ::                                                 r2
!
   real, dimension( its:ite, kts:kte )                                       , &
         intent(inout)  ::                                                 au, &
                                                                           cu, &
                                                                           f1
   real, dimension( its:ite, kts:kte,nt )                                    , &
         intent(inout)  ::                                                 f2
!
   real    :: fk
   integer :: i,k,l,n,it
!
!-------------------------------------------------------------------------------
!
   l = ite
   n = kte
!
   do i = its,l
     fk = 1./cm(i,1)
     au(i,1) = fk*cu(i,1)
     f1(i,1) = fk*r1(i,1)
   enddo
!
   do it = 1,nt
     do i = its,l
       fk = 1./cm(i,1)
       f2(i,1,it) = fk*r2(i,1,it)
     enddo
   enddo
!
   do k = kts+1,n-1
     do i = its,l
       fk = 1./(cm(i,k)-cl(i,k)*au(i,k-1))
       au(i,k) = fk*cu(i,k)
       f1(i,k) = fk*(r1(i,k)-cl(i,k)*f1(i,k-1))
     enddo
   enddo
!
   do it = 1,nt
     do k = kts+1,n-1
       do i = its,l
         fk = 1./(cm(i,k)-cl(i,k)*au(i,k-1))
         f2(i,k,it) = fk*(r2(i,k,it)-cl(i,k)*f2(i,k-1,it))
       enddo
     enddo
   enddo
!
   do i = its,l
     fk = 1./(cm(i,n)-cl(i,n)*au(i,n-1))
     f1(i,n) = fk*(r1(i,n)-cl(i,n)*f1(i,n-1))
   enddo
!
   do it = 1,nt
     do i = its,l
       fk = 1./(cm(i,n)-cl(i,n)*au(i,n-1))
       f2(i,n,it) = fk*(r2(i,n,it)-cl(i,n)*f2(i,n-1,it))
     enddo
   enddo
!
   do k = n-1,kts,-1
     do i = its,l
       f1(i,k) = f1(i,k)-au(i,k)*f1(i,k+1)
     enddo
   enddo
!
   do it = 1,nt
     do k = n-1,kts,-1
       do i = its,l
         f2(i,k,it) = f2(i,k,it)-au(i,k)*f2(i,k+1,it)
       enddo
     enddo
   enddo
!
   end subroutine tridi1n
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine tridin_ysu(cl,cm,cu,r2,au,f2,its,ite,kts,kte,nt)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
   integer, intent(in )      ::     its,ite, kts,kte, nt
!
   real, dimension( its:ite, kts+1:kte+1 )                                   , &
         intent(in   )  ::                                                 cl
!
   real, dimension( its:ite, kts:kte )                                       , &
         intent(in   )  ::                                                 cm
   real, dimension( its:ite, kts:kte,nt )                                    , &
         intent(in   )  ::                                                 r2
!
   real, dimension( its:ite, kts:kte )                                       , &
         intent(inout)  ::                                                 au, &
                                                                           cu
   real, dimension( its:ite, kts:kte,nt )                                    , &
         intent(inout)  ::                                                 f2
!
   real    :: fk
   integer :: i,k,l,n,it
!
!-------------------------------------------------------------------------------
!
   l = ite
   n = kte
!
   do it = 1,nt
     do i = its,l
       fk = 1./cm(i,1)
       au(i,1) = fk*cu(i,1)
       f2(i,1,it) = fk*r2(i,1,it)
     enddo
   enddo
!
   do it = 1,nt
     do k = kts+1,n-1
       do i = its,l
         fk = 1./(cm(i,k)-cl(i,k)*au(i,k-1))
         au(i,k) = fk*cu(i,k)
         f2(i,k,it) = fk*(r2(i,k,it)-cl(i,k)*f2(i,k-1,it))
       enddo
     enddo
   enddo
!
   do it = 1,nt
     do i = its,l
       fk = 1./(cm(i,n)-cl(i,n)*au(i,n-1))
       f2(i,n,it) = fk*(r2(i,n,it)-cl(i,n)*f2(i,n-1,it))
     enddo
   enddo
!
   do it = 1,nt
     do k = n-1,kts,-1
       do i = its,l
         f2(i,k,it) = f2(i,k,it)-au(i,k)*f2(i,k+1,it)
       enddo
     enddo
   enddo
!
   end subroutine tridin_ysu
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine ysuinit(rublten,rvblten,rthblten,rqvblten,                       &
                      rqcblten,rqiblten,p_qi,p_first_scalar,                   &
                      restart, allowed_to_read,                                &
                      ids, ide, jds, jde, kds, kde,                            &
                      ims, ime, jms, jme, ktes, ktee,                            &
                      its, ite, jts, jte, kts, kte                 )
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
   logical , intent(in)          :: restart, allowed_to_read
   integer , intent(in)          ::  ids, ide, jds, jde, kds, kde,             &
                                     ims, ime, jms, jme, ktes, ktee,             &
                                     its, ite, jts, jte, kts, kte
   integer , intent(in)          ::  p_qi,p_first_scalar
   real , dimension( ims:ime , ktes:ktee , jms:jme ), intent(out) ::             &
                                                                      rublten, &
                                                                      rvblten, &
                                                                     rthblten, &
                                                                     rqvblten, &
                                                                     rqcblten, &
                                                                     rqiblten
   integer :: i, j, k, itf, jtf, ktf
!
   jtf = min0(jte,jde-1)
   ktf = min0(kte,kde-1)
   itf = min0(ite,ide-1)
!
   if(.not.restart)then
     do j = jts,jtf
       do k = kts,ktf
         do i = its,itf
            rublten(i,k,j) = 0.
            rvblten(i,k,j) = 0.
            rthblten(i,k,j) = 0.
            rqvblten(i,k,j) = 0.
            rqcblten(i,k,j) = 0.
         enddo
       enddo
     enddo
   endif
!
   if (p_qi .ge. p_first_scalar .and. .not.restart) then
     do j = jts,jtf
       do k = kts,ktf
         do i = its,itf
           rqiblten(i,k,j) = 0.
         enddo
       enddo
     enddo
   endif
!
   end subroutine ysuinit
!-------------------------------------------------------------------------------
! ==================================================================

      SUBROUTINE GET_PBLH(KTS,KTE,zi,thetav1D,qke1D,zw1D,dz1D,landsea)
! Copied from MYNN PBL

      !---------------------------------------------------------------
      !             NOTES ON THE PBLH FORMULATION
      !
      !The 1.5-theta-increase method defines PBL heights as the level at
      !which the potential temperature first exceeds the minimum potential
      !temperature within the boundary layer by 1.5 K. When applied to
      !observed temperatures, this method has been shown to produce PBL-
      !height estimates that are unbiased relative to profiler-based
      !estimates (Nielsen-Gammon et al. 2008). However, their study did not
      !include LLJs. Banta and Pichugina (2008) show that a TKE-based
      !threshold is a good estimate of the PBL height in LLJs. Therefore,
      !a hybrid definition is implemented that uses both methods, weighting
      !the TKE-method more during stable conditions (PBLH < 400 m).
      !A variable tke threshold (TKEeps) is used since no hard-wired
      !value could be found to work best in all conditions.
      !---------------------------------------------------------------

      INTEGER,INTENT(IN) :: KTS,KTE
      REAL, INTENT(OUT) :: zi
      REAL, INTENT(IN) :: landsea
      REAL, DIMENSION(KTS:KTE), INTENT(IN) :: thetav1D, qke1D, dz1D
      REAL, DIMENSION(KTS:KTE+1), INTENT(IN) :: zw1D
      !LOCAL VARS
      REAL ::  PBLH_TKE,qtke,qtkem1,wt,maxqke,TKEeps,minthv
      REAL :: delt_thv   !delta theta-v; dependent on land/sea point
      REAL, PARAMETER :: sbl_lim  = 200. !Theta-v PBL lower limit of trust (m).
      REAL, PARAMETER :: sbl_damp = 400. !Damping range for averaging with TKE-based PBLH (m).
      INTEGER :: I,J,K,kthv,ktke

      !FIND MAX TKE AND MIN THETAV IN THE LOWEST 500 M
      k = kts+1
      kthv = 1
      ktke = 1
      maxqke = 0.
      minthv = 9.E9

      DO WHILE (zw1D(k) .LE. 500.)
        qtke  =MAX(Qke1D(k),0.)   ! maximum QKE
         IF (maxqke < qtke) then
            maxqke = qtke
            ktke = k
         ENDIF
         IF (minthv > thetav1D(k)) then
             minthv = thetav1D(k)
             kthv = k
         ENDIF
         k = k+1
      ENDDO
      !TKEeps = maxtke/20. = maxqke/40.
      TKEeps = maxqke/40.
      TKEeps = MAX(TKEeps,0.025)
      TKEeps = MIN(TKEeps,0.25)

      !FIND THETAV-BASED PBLH (BEST FOR DAYTIME).
      zi=0.
      k = kthv+1
      IF((landsea-1.5).GE.0)THEN
      ! WATER
          delt_thv = 0.75
      ELSE
      ! LAND
          delt_thv = 1.5
      ENDIF

      zi=0.
      k = kthv+1
      DO WHILE (zi .EQ. 0.)
         IF (thetav1D(k) .GE. (minthv + delt_thv))THEN
            zi = zw1D(k) - dz1D(k-1)* &
              & MIN((thetav1D(k)-(minthv + delt_thv))/MAX(thetav1D(k)-thetav1D(k-1),1E-6),1.0)
        ENDIF
        k = k+1
         IF (k .EQ. kte-1) zi = zw1D(kts+1) !EXIT SAFEGUARD
      ENDDO

      !print*,"IN GET_PBLH:",thsfc,zi
      !FOR STABLE BOUNDARY LAYERS, USE TKE METHOD TO COMPLEMENT THE
      !THETAV-BASED DEFINITION (WHEN THE THETA-V BASED PBLH IS BELOW ~0.5 KM).
      !THE TANH WEIGHTING FUNCTION WILL MAKE THE TKE-BASED DEFINITION NEGLIGIBLE
      !WHEN THE THETA-V-BASED DEFINITION IS ABOVE ~1 KM.
      !FIND TKE-BASED PBLH (BEST FOR NOCTURNAL/STABLE CONDITIONS).

      PBLH_TKE=0.
      k = ktke+1
     DO WHILE (PBLH_TKE .EQ. 0.)
        !QKE CAN BE NEGATIVE (IF CKmod == 0)... MAKE TKE NON-NEGATIVE.
         qtke  =MAX(Qke1D(k)/2.,0.)      ! maximum TKE
         qtkem1=MAX(Qke1D(k-1)/2.,0.)
         IF (qtke .LE. TKEeps) THEN
               PBLH_TKE = zw1D(k) - dz1D(k-1)* &
               & MIN((TKEeps-qtke)/MAX(qtkem1-qtke, 1E-6), 1.0)
             !IN CASE OF NEAR ZERO TKE, SET PBLH = LOWEST LEVEL.
             PBLH_TKE = MAX(PBLH_TKE,zw1D(kts+1))
             !print *,"PBLH_TKE:",i,j,PBLH_TKE, Qke1D(k)/2., zw1D(kts+1)
         ENDIF
         k = k+1
         IF (k .EQ. kte-1) PBLH_TKE = zw1D(kts+1) !EXIT SAFEGUARD
      ENDDO

    !BLEND THE TWO PBLH TYPES HERE:

      wt=.5*TANH((zi - sbl_lim)/sbl_damp) + .5
      zi=PBLH_TKE*(1.-wt) + zi*wt

   END SUBROUTINE GET_PBLH
! ==================================================================

!-------------------------------------------------------------------------------
