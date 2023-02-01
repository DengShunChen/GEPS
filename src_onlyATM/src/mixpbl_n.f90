      subroutine mixpbl_n(ix,im,km,ntrac,                               &
           u1,v1,t1,q1,                                                 &
           psk,rbsoil,fm,fh,tsea,heat,evap,stress,spd1,kpbl,            &
           prsi,del,prsl,prslk,phii,phil,rcl,deltim,                    &
           hpbl,ice,grav,cp,hvap,rd,j)
!
!     use machine     , only : kind_phys
!     use physcons, grav => con_g, rd => con_rd, cp => con_cp
!     ,             hvap => con_hvap, rog => con_rog, fv => con_fvirt
      use const, only: RTYPE
!
      implicit none
!
      real grav,cp,hvap,rd
!
!     arguments
!
      integer ix, im, km, ntrac, kpbl(im),j
!
      real deltim
      real dv(im,km),     du(im,km),                                    &
                           tau(im,km),    rtg(im,km,ntrac),             &
                           u1(ix,km),     v1(ix,km),                    &
                           t1(ix,km),     q1(ix,km,ntrac),              &
                           rbsoil(im),                                  &
!                          cd(im),        ch(im),
                           fm(im),        fh(im),                       &
                           tsea(im),      qss(im),                      &
                                          spd1(im),                     &
!                          dphi(im),      spd1(im),
                           prsi(ix,km+1), del(ix,km),                   &
                           prsl(ix,km),   prslk(ix,km),                 &
                           phii(ix,km+1), phil(ix,km),                  &
                           rcl(im),       dusfc(im),                    &
                           dvsfc(im),     dtsfc(im),                    &
                           dqsfc(im),     hpbl(im),                     &
                           hgamt(im),     hgamq(im)
      real(kind=RTYPE)     psk(im)
!
!    locals
!
      integer i,iprt,is,iun,k,kk,kmpbl,lond
      real evap(im),  heat(im),    phih(im),                            &
                           phim(im),  rbdn(im),    rbup(im),            &
                           the1(im),  stress(im),  beta(im),            &
                           the1v(im), thekv(im),   thermal(im),         &
                           thesv(im), ustar(im),   wscale(im)
!
      real rdzt(im,km-1),                                               &
                           zi(im,km+1),     zl(im,km),                  &
                           dku(im,km-1),    dkt(im,km-1),               &
                           al(im,km-1),     ad(im,km),                  &
                           au(im,km-1),     a1(im,km),                  &
                           a2(im,km*ntrac), theta(im,km)
      logical              pblflg(im),   sfcflg(im), stable(im)         &
                          ,ice(im)
!
      real aphi16,  aphi5,  bet1,   bvf2,                               &
                           cfac,    conq,   cont,   conw,               &
                           conwrc,  dk,     dkmax,  dkmin,              &
                           dq1,     dsdz2,  dsdzq,  dsdzt,              &
                           dsig,    dt,     dthe1,  dtodsd,             &
                           dtodsu,  dw2,    dw2min, g,                  &
                           gamcrq,  gamcrt, gocp,   gor, gravi,         &
                           hol,     pfac,   prmax,  prmin, prinv,       &
                           prnum,   qmin,   qtend,  rbcr,               &
                           rbint,   rdt,    rdz,                        &
                           ri,      rimin,  rl2,    rlam,  rlamun,      &
                           rone,   rzero,   sfcfrac,                    &
                           sflux,   shr2,   spdk2,  sri,                &
                           tem,     ti,     ttend,  tvd,                &
                           tvu,     utend,  vk,     vk2,                &
                           vpert,   vtend,  xkzo(im,km),   zfac,        &
                           zfmin,   zk,     tem1, xkzm
!c
      real rv,fv,rog
!     parameter (rv=4.615e+2,fv=rv/rd-1.,rog=rd/grav)
!
!     parameter (gravi=1.0/grav)
!     parameter(g=grav)
!     parameter(gor=g/rd,gocp=g/cp)
!     parameter(cont=1000.*cp/g,conq=1000.*hvap/g,conw=1000./g)
      parameter(rlam=30.0,vk=0.4,vk2=vk*vk,prmin=1.0,prmax=4.)
      parameter(dw2min=0.0001,dkmin=0.0,dkmax=1000.,rimin=-100.)
      parameter(rbcr=0.25,cfac=7.8,pfac=2.0,sfcfrac=0.1)
!test
!     parameter(rbcr=0.5,cfac=7.8,pfac=2.0,sfcfrac=0.1)
!     parameter(qmin=1.e-8,xkzm=3.0,zfmin=1.e-8,aphi5=5.,aphi16=16.)
!     parameter(qmin=1.e-8,xkzm=2.0,zfmin=1.e-8,aphi5=5.,aphi16=16.)
!     parameter(qmin=1.e-8,xkzm=0.5,zfmin=1.e-8,aphi5=5.,aphi16=16.)
!     parameter(qmin=1.e-8,xkzm=0.25,zfmin=1.e-8,aphi5=5.,aphi16=16.)
!     parameter(qmin=1.e-8,xkzm=0.0,zfmin=1.e-8,aphi5=5.,aphi16=16.)
      parameter(qmin=1.e-8,xkzm=1.0,zfmin=1.e-8,aphi5=5.,aphi16=16.)
!test090708
!     parameter(qmin=1.e-8,xkzm=0.2,zfmin=1.e-8,aphi5=5.,aphi16=16.)
!     parameter(gamcrt=3.,gamcrq=2.e-3)
!     parameter(gamcrt=3.,gamcrq=0., rlamun=30.0)
      parameter(gamcrt=3.,gamcrq=0., rlamun=150.0)
      parameter(iun=84)
!
!CWBinit
      rtg=0.
      tau=0.
      dtsfc=0.
      dqsfc=0.
      dusfc=0.
      dvsfc=0.
      du=0.
      dv=0.

      gravi=1.0/grav
      g=grav
      rv=4.615e+2
      fv=rv/rd-1.
      rog=rd/grav
      gor=g/rd
      gocp=g/cp
      cont=1000.*cp/g
      conq=1000.*hvap/g
      conw=1000./g
!
!-----------------------------------------------------------------------
!
 601  format(1x,' moninp lat lon step hour ',3i6,f6.1)
 602      format(1x,'    k','        z','        t','       th',        &
           '      tvh','        q','        u','        v',             &
           '       sp')
 603      format(1x,i5,8f9.1)
 604      format(1x,'  sfc',9x,f9.1,18x,f9.1)
 605      format(1x,'    k      zl    spd2   thekv   the1v'             &
               ,' thermal    rbup')
 606      format(1x,i5,6f8.2)
 607      format(1x,' kpbl    hpbl      fm      fh   hgamt',            &
               '   hgamq      ws   ustar      cd      ch')
 608      format(1x,i5,9f8.2)
 609      format(1x,' k pr dkt dku ',i5,3f8.2)
 610      format(1x,' k pr dkt dku ',i5,3f8.2,' l2 ri t2',              &
               ' sr2  ',2f8.2,2e10.2)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     compute preliminary variables
!
      if (ix .lt. im) stop
!
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
!     dt    = 2. * deltim
!phon
      dt    =  deltim
      rdt   = 1. / dt
      kmpbl = km / 2
!
      do k=1,km
        do i=1,im
          zi(i,k) = phii(i,k) * gravi
          zl(i,k) = phil(i,k) * gravi
        enddo
      enddo
!
      do k=1,kmpbl
        do i=1,im
          theta(i,k) = t1(i,k) * psk(i) / prslk(i,k)
        enddo
      enddo
!
      do k = 1,km-1
        do i=1,im
          rdzt(i,k) = 1.0 / (zl(i,k+1) - zl(i,k))
!         rdzt(i,k) = gor * prsi(i,k+1) / (prsl(i,k) - prsl(i,k+1))
!         if (prsi(i,1) .gt. 60.0) then
!           tem1    = max((prsi(i,k+1)-30.0)/(prsi(i,1)-30.0), 0.0)
!         else
!           tem1    = 0.0
!         endif
!         xkzo(i,k) = xkzm * tem1 * tem1
!!        tem1      = (zi(i,k+1) - zi(i,1)) * 0.002
!!        xkzo(i,k) = xkzm * exp(-tem1)
        enddo
      enddo
!
      do i = 1,im
!p        dusfc(i) = 0.
!p        dvsfc(i) = 0.
!p        dusfc(i) = 0.
!p        dtsfc(i) = 0.
!p        dqsfc(i) = 0.
         hgamt(i) = 0.
         hgamq(i) = 0.
         wscale(i) = 0.
         kpbl(i) = 1
         hpbl(i) = zi(i,2)
         pblflg(i) = .true.
         sfcflg(i) = .true.
         if(rbsoil(i).gt.0.0) sfcflg(i) = .false.
      enddo
!!
      do i=1,im
!        rdzt1    = gor * prsl(i,1) / del(i,1)
!        bet1     = dt*rdzt1*spd1(i)/t1(i,1)
!        beta(i)  = dt*rdzt1/t1(i,1)
         beta(i)  = dt / (zi(i,2)-zi(i,1))
!        betaw(i) = bet1*cd(i)
!        betat(i) = bet1*ch(i)
!        betaq(i) = dphi(i)*betat(i)
      enddo
!
      do i=1,im
!        zl1(i) = 0.-(t1(i,1)+tsea(i))/2.*log(prsl(i,1)/prsi(i,1))*rog
!        ustar(i) = sqrt(cd(i)*spd1(i)**2)
         ustar(i) = sqrt(stress(i))
      enddo
!
      do i=1,im
!        thesv(i)   = tsea(i)*(1.+fv*max(qss(i),qmin))
         the1(i)    = theta(i,1)
         the1v(i)   = the1(i)*(1.+fv*max(q1(i,1,1),qmin))
         thermal(i) = the1v(i)
!        dthe1      = (the1(i)-tsea(i))
!        dq1        = (max(q1(i,1,1),qmin) - max(qss(i),qmin))
!        heat(i)    = -ch(i)*spd1(i)*dthe1
!        evap(i)    = -ch(i)*spd1(i)*dq1
      enddo
!
!
!     compute the first guess of pbl height
!
      do i=1,im
         stable(i) = .false.
!        zl(i,1) = zl1(i)
         rbup(i) = rbsoil(i)
! phon add
      if(ice(i) .and. .not.sfcflg(i))stable(i)=.true.
! phon
      enddo
      do k = 2, kmpbl
        do i = 1, im
          if(.not.stable(i)) then
             rbdn(i)   = rbup(i)
!            zl(i,k)   = zl(i,k-1) - (t1(i,k)+t1(i,k-1))/2 *
!    &                   log(prsl(i,k)/prsl(i,k-1)) * rog
             thekv(i)  = theta(i,k)*(1.+fv*max(q1(i,k,1),qmin))
             spdk2     = max(rcl(i)*(u1(i,k)**2+v1(i,k)**2),1.)
             rbup(i)   = (thekv(i)-the1v(i))*(g*zl(i,k)/the1v(i))/spdk2
             kpbl(i)   = k
             stable(i) = rbup(i).gt.rbcr
          endif
        enddo
      enddo
!
      do i = 1,im
! phon add
      if(kpbl(i).gt.1)then
         k = kpbl(i)
         if(rbdn(i).ge.rbcr) then
            rbint = 0.
         elseif(rbup(i).le.rbcr) then
            rbint = 1.
         else
            rbint = (rbcr-rbdn(i))/(rbup(i)-rbdn(i))
         endif
         hpbl(i) = zl(i,k-1) + rbint*(zl(i,k)-zl(i,k-1))
         if(hpbl(i).lt.zi(i,kpbl(i))) kpbl(i) = kpbl(i) - 1
       endif
! phon add
      enddo
!!
      do i=1,im
           hol = max(rbsoil(i)*fm(i)*fm(i)/fh(i),rimin)
           if(sfcflg(i)) then
              hol = min(hol,-zfmin)
           else
              hol = max(hol,zfmin)
           endif
!
!          hol = hol*hpbl(i)/zl1(i)*sfcfrac
           hol = hol*hpbl(i)/zl(i,1)*sfcfrac
           if(sfcflg(i)) then
!             phim = (1.-aphi16*hol)**(-1./4.)
!             phih = (1.-aphi16*hol)**(-1./2.)
              tem  = 1.0 / (1. - aphi16*hol)
              phih(i) = sqrt(tem)
              phim(i) = sqrt(phih(i))
           else
              phim(i) = (1.+aphi5*hol)
              phih(i) = phim(i)
           endif
           wscale(i) = ustar(i)/phim(i)
           wscale(i) = min(wscale(i),ustar(i)*aphi16)
           wscale(i) = max(wscale(i),ustar(i)/aphi5)
      enddo
!
!     compute the surface variables for pbl height estimation
!     under unstable conditions
!
      do i = 1,im
         sflux  = heat(i) + evap(i)*fv*the1(i)
         if(sfcflg(i).and.sflux.gt.0.0) then
           hgamt(i)   = min(cfac*heat(i)/wscale(i),gamcrt)
           hgamq(i)   = min(cfac*evap(i)/wscale(i),gamcrq)
           vpert      = hgamt(i) + fv*the1(i)*hgamq(i)
           vpert      = min(vpert,gamcrt)
           thermal(i) = thermal(i) + max(vpert,0.)
           hgamt(i)   = max(hgamt(i),0.0)
           hgamq(i)   = max(hgamq(i),0.0)
         else
           pblflg(i) = .false.
         endif
      enddo
!
      do i = 1,im
         if(pblflg(i)) then
            kpbl(i) = 1
            hpbl(i) = zi(i,2)
         endif
      enddo
!
!     enhance the pbl height by considering the thermal
!
      do i = 1, im
         if(pblflg(i)) then
            stable(i) = .false.
            rbup(i) = rbsoil(i)
         endif
      enddo
      do k = 2, kmpbl
        do i = 1, im
          if(.not.stable(i).and.pblflg(i)) then
            rbdn(i)   = rbup(i)
!           zl(i,k)   = zl(i,k-1) - (t1(i,k)+t1(i,k-1))/2 *
!    &                  log(prsl(i,k)/prsl(i,k-1)) * rog
            thekv(i)  = theta(i,k)*(1.+fv*max(q1(i,k,1),qmin))
            spdk2     = max(rcl(i)*(u1(i,k)**2+v1(i,k)**2),1.)
            rbup(i)   = (thekv(i)-thermal(i))*(g*zl(i,k)/the1v(i))/spdk2
            kpbl(i)   = k
            stable(i) = rbup(i).gt.rbcr
          endif
        enddo
      enddo
!
      do i = 1,im
         if(pblflg(i)) then
            k = kpbl(i)
            if(rbdn(i).ge.rbcr) then
               rbint = 0.
            elseif(rbup(i).le.rbcr) then
               rbint = 1.
            else
               rbint = (rbcr-rbdn(i))/(rbup(i)-rbdn(i))
            endif
            hpbl(i) = zl(i,k-1) + rbint*(zl(i,k)-zl(i,k-1))
            if(hpbl(i).lt.zi(i,kpbl(i))) kpbl(i) = kpbl(i) - 1
            if(kpbl(i).le.1) pblflg(i) = .false.
         endif
      enddo
!!
!     do k = 1,km-1
!       do i=1,im
!         if (rbsoil(i) .gt. 0.0 .or. (.not. pblflg(i))) then
!           tem1    = max((prsi(i,k+1)-30.0)/(prsi(i,1)-30.0), 0.0)
!           xkzo(i,k) = xkzm * tem1 * tem1
!!        tem1      = (zi(i,k+1) - zi(i,1)) * 0.002
!!        xkzo(i,k) = xkzm * exp(-tem1)
!         else
!           xkzo(i,k) = 0.0
!
!         if (pblflg(i)) then
!         if (sfcflg(i)) then
!           xkzo(i,k) = 0.0
!         else
!           tem1      = (zi(i,k+1) - zi(i,1)) * 0.0005
!           tem1      = (100.0 - prsi(i,k+1)) * 0.075
!           tem1      = 100.0 - prsi(i,k+1)
!           tem1      = max(0.0, 100.0 - prsi(i,k+1))
!           tem1      = tem1 * tem1 * 0.00075
!           tem1      = tem1 * tem1 * 0.001
!           tem1      = tem1 * tem1 * 0.0011
!           tem1      = tem1 * tem1 * 0.0012
!
!           tem1      = 1.0 - prsi(i,k+1) / prsi(i,1)
!           tem1      = tem1 * tem1 * 5.0
!           tem1      = tem1 * tem1 * 7.5
!           tem1      = tem1 * tem1 * 10.0
!           tem1      = tem1 * tem1 * 12.0
!
!           xkzo(i,k) = xkzm * min(1.0, exp(-tem1))
!           if (xkzo(i,k) .lt. 0.01) xkzo(i,k) = 0.0
!         endif
!       enddo
!     enddo
!
      do k = 1,km-1
        do i=1,im
          tem1      = 1.0 - prsi(i,k+1) / prsi(i,1)
          tem1      = tem1 * tem1 * 10.0
          xkzo(i,k) = xkzm * min(1.0, exp(-tem1))
        enddo
!     write(77,*)'k=',k,'j=',j,'xkzo=',xkzo(360,k)
      enddo
!!
!
!     compute diffusion coefficients below pbl
!
      do k = 1, kmpbl
         do i=1,im
            if(kpbl(i).gt.k) then
               prinv = 1.0 / (phih(i)/phim(i)+cfac*vk*.1)
               prinv = min(prinv,prmax)
               prinv = max(prinv,prmin)
!              zfac = max((1.-(zi(i,k+1)-zl1(i))/
!    1                (hpbl(i)-zl1(i))), zfmin)
               zfac = max((1.-(zi(i,k+1)-zl(i,1))/  &
                      (hpbl(i)-zl(i,1))), zfmin)
               dku(i,k) = xkzo(i,k) + wscale(i)*vk*zi(i,k+1)  &
                               * zfac**pfac
               dkt(i,k) = dku(i,k)*prinv
               dku(i,k) = min(dku(i,k),dkmax)
               dku(i,k) = max(dku(i,k),dkmin)
               dkt(i,k) = min(dkt(i,k),dkmax)
               dkt(i,k) = max(dkt(i,k),dkmin)
            endif
         enddo
      enddo
!
!     compute diffusion coefficients over pbl (free atmosphere)
!
      do k = 1, km-1
         do i=1,im
            if(k.ge.kpbl(i)) then
!              ti   = 0.5*(t1(i,k)+t1(i,k+1))
               ti   = 2.0 / (t1(i,k)+t1(i,k+1))
!              rdz  = rdzt(i,k)/ti
!              rdz  = rdzt(i,k) * ti
               rdz  = rdzt(i,k)

               dw2  = rcl(i)*((u1(i,k)-u1(i,k+1))**2  &
                            + (v1(i,k)-v1(i,k+1))**2)
               shr2 = max(dw2,dw2min)*rdz*rdz
               tvd  = t1(i,k)*(1.+fv*max(q1(i,k,1),qmin))
               tvu  = t1(i,k+1)*(1.+fv*max(q1(i,k+1,1),qmin))
!              bvf2 = g*(gocp+rdz*(tvu-tvd))/ti
               bvf2 = g*(gocp+rdz*(tvu-tvd)) * ti
               ri   = max(bvf2/shr2,rimin)
               zk   = vk*zi(i,k+1)
!              rl2  = (zk*rlam/(rlam+zk))**2
!              dk   = rl2*sqrt(shr2)
!              rl2  = zk*rlam/(rlam+zk)
!              dk   = rl2*rl2*sqrt(shr2)
               if(ri.lt.0.) then ! unstable regime
                  rl2      = zk*rlamun/(rlamun+zk)
                  dk       = rl2*rl2*sqrt(shr2)
                  sri      = sqrt(-ri)
                  dku(i,k) = xkzo(i,k) + dk*(1+8.*(-ri)/(1+1.746*sri))
                  dkt(i,k) = xkzo(i,k) + dk*(1+8.*(-ri)/(1+1.286*sri))
               else             ! stable regime
                  rl2       = zk*rlam/(rlam+zk)
!                 tem       = rlam * sqrt(0.01*prsi(i,k))
!                 rl2       = zk*tem/(tem+zk)
                  dk        = rl2*rl2*sqrt(shr2)
                  dkt(i,k)  = xkzo(i,k) + dk/(1+5.*ri)**2
                  prnum     = 1.0 + 2.1*ri
                  prnum     = min(prnum,prmax)
                  dku(i,k)  = (dkt(i,k)-xkzo(i,k))*prnum + xkzo(i,k)
               endif
!
!
               dku(i,k) = min(dku(i,k),dkmax)
               dku(i,k) = max(dku(i,k),dkmin)
               dkt(i,k) = min(dkt(i,k),dkmax)
               dkt(i,k) = max(dkt(i,k),dkmin)
         if(ice(i) .and.  .not.pblflg(i))dkt(i,k)=min(dkt(i,k), 0.005)
         if(ice(i) .and.  .not.pblflg(i))dku(i,k)=min(dku(i,k), 0.005)
!
!cc   if(i.eq.lond.and.lat.eq.latd) then
!cc   prnum = dku(k)/dkt(k)
!cc   write(iun,610) k,prnum,dkt(k),dku(k),rl2,ri,
!cc   1              bvf2,shr2
!cc   endif
!
            endif
         enddo
      enddo
!
!     compute tridiagonal matrix elements for heat and moisture
!
      do i=1,im
         ad(i,1) = 1.
         a1(i,1) = t1(i,1)   + beta(i) * heat(i)
         a2(i,1) = q1(i,1,1) + beta(i) * evap(i)
!        a1(i,1) = t1(i,1)-betat(i)*(theta(i,1)-tsea(i))
!        a2(i,1) = q1(i,1,1)-betaq(i)*
!    &           (max(q1(i,1,1),qmin)-max(qss(i),qmin))
      enddo
!byl      if(ntrac.eq.2) then
      if(ntrac.ge.2) then
        do k = 2, ntrac
          is = (k-1) * km
          do i = 1, im
            a2(i,1+is) = q1(i,1,k)
          enddo
        enddo
      endif
!
      do k = 1,km-1
        do i = 1,im
          dtodsd = dt/del(i,k)
          dtodsu = dt/del(i,k+1)
          dsig   = prsl(i,k)-prsl(i,k+1)
!         rdz    = rdzt(i,k)*2./(t1(i,k)+t1(i,k+1))
          rdz    = rdzt(i,k)
          tem1   = dsig * dkt(i,k) * rdz
          if(pblflg(i).and.k.lt.kpbl(i)) then
!            dsdzt = dsig*dkt(i,k)*rdz*(gocp-hgamt(i)/hpbl(i))
!            dsdzq = dsig*dkt(i,k)*rdz*(-hgamq(i)/hpbl(i))
             tem   = 1.0 / hpbl(i)
             dsdzt = tem1 * (gocp-hgamt(i)*tem)
             dsdzq = tem1 * (-hgamq(i)*tem)
             a2(i,k)   = a2(i,k)+dtodsd*dsdzq
             a2(i,k+1) = q1(i,k+1,1)-dtodsu*dsdzq
          else
!            dsdzt = dsig*dkt(i,k)*rdz*(gocp)
             dsdzt = tem1 * gocp
             a2(i,k+1) = q1(i,k+1,1)
          endif
!         dsdz2 = dsig*dkt(i,k)*rdz*rdz
          dsdz2     = tem1 * rdz
          au(i,k)   = -dtodsd*dsdz2
          al(i,k)   = -dtodsu*dsdz2
          ad(i,k)   = ad(i,k)-au(i,k)
          ad(i,k+1) = 1.-al(i,k)
          a1(i,k)   = a1(i,k)+dtodsd*dsdzt
          a1(i,k+1) = t1(i,k+1)-dtodsu*dsdzt
        enddo
      enddo
!byl      if(ntrac.eq.2) then
      if(ntrac.ge.2) then
        do kk = 2, ntrac
          is = (kk-1) * km
          do k = 1, km - 1
            do i = 1, im
              a2(i,k+1+is) = q1(i,k+1,kk)
            enddo
          enddo
        enddo
      endif
!
!     solve tridiagonal problem for heat and moisture
!
      call tridin(im,km,ntrac,al,ad,au,a1,a2,au,a1,a2)
!
!     recover tendencies of heat and moisture
!
      do  k = 1,km
         do i = 1,im
            ttend      = (a1(i,k)-t1(i,k))*rdt
            qtend      = (a2(i,k)-q1(i,k,1))*rdt
            tau(i,k)   = tau(i,k)+ttend
            rtg(i,k,1) = rtg(i,k,1)+qtend
            dtsfc(i)   = dtsfc(i)+cont*del(i,k)*ttend
            dqsfc(i)   = dqsfc(i)+conq*del(i,k)*qtend
         enddo
      enddo
!byl      if(ntrac.eq.2) then
      if(ntrac.ge.2) then
        do kk = 2, ntrac
          is = (kk-1) * km
          do k = 1, km
            do i = 1, im
              qtend = (a2(i,k+is)-q1(i,k,kk))*rdt
              rtg(i,k,kk) = rtg(i,k,kk)+qtend
            enddo
          enddo
        enddo
      endif
!
      do  k = 1,km
         do i = 1,im
            t1(i,k)      = a1(i,k)
            q1(i,k,1)    = a2(i,k)
         enddo
      enddo
!
!byl       if(ntrac.eq.2)then
       if(ntrac.ge.2)then
        do kk = 2, ntrac
          is = (kk-1) * km
          do k = 1, km
            do i = 1, im
            q1(i,k,kk) = a2(i,k+is)
            enddo
          enddo
        enddo
       endif
!
!     compute tridiagonal matrix elements for momentum
!
      do i=1,im
!        ad(i,1) = 1.+betaw(i)
         ad(i,1) = 1.0 + beta(i) * stress(i) / spd1(i)
         a1(i,1) = u1(i,1)
         a2(i,1) = v1(i,1)
!        ad(i,1) = 1.0
!        tem     = 1.0 + beta(i) * stress(i) / spd1(i)
!        a1(i,1) = u1(i,1) * tem
!        a2(i,1) = v1(i,1) * tem
      enddo
!
      do k = 1,km-1
        do i=1,im
          dtodsd    = dt/del(i,k)
          dtodsu    = dt/del(i,k+1)
          dsig      = prsl(i,k)-prsl(i,k+1)
!         rdz       = rdzt(i,k)*2./(t1(i,k)+t1(i,k+1))
          rdz       = rdzt(i,k)
          dsdz2     = dsig*dku(i,k)*rdz*rdz
          au(i,k)   = -dtodsd*dsdz2
          al(i,k)   = -dtodsu*dsdz2
          ad(i,k)   = ad(i,k)-au(i,k)
          ad(i,k+1) = 1.-al(i,k)
          a1(i,k+1) = u1(i,k+1)
          a2(i,k+1) = v1(i,k+1)
        enddo
      enddo
!
!     solve tridiagonal problem for momentum
!
      call tridi2(im,km,al,ad,au,a1,a2,au,a1,a2)
!
!     recover tendencies of momentum
!
      do k = 1,km
         do i = 1,im
            conwrc = conw*sqrt(rcl(i))
            utend = (a1(i,k)-u1(i,k))*rdt
            vtend = (a2(i,k)-v1(i,k))*rdt
            du(i,k)  = du(i,k)+utend
            dv(i,k)  = dv(i,k)+vtend
            dusfc(i) = dusfc(i)+conwrc*del(i,k)*utend
            dvsfc(i) = dvsfc(i)+conwrc*del(i,k)*vtend
         enddo
      enddo
!!
      do  k = 1,km
         do i = 1,im
            u1(i,k)      = a1(i,k)
            v1(i,k)      = a2(i,k)
         enddo
      enddo
      return
      end
!fpp$ noconcur r
!-----------------------------------------------------------------------
      subroutine tridi2(l,n,cl,cm,cu,r1,r2,au,a1,a2)
!sela %include dbtridi2;
!c
!     use machine     , only : kind_phys
      implicit none
      integer             k,n,l,i
      real fk
!c
      real cl(l,2:n),cm(l,n),cu(l,n-1),r1(l,n),r2(l,n), &
                au(l,n-1),a1(l,n),a2(l,n)
!-----------------------------------------------------------------------
      do i=1,l
        fk      = 1./cm(i,1)
        au(i,1) = fk*cu(i,1)
        a1(i,1) = fk*r1(i,1)
        a2(i,1) = fk*r2(i,1)
      enddo
      do k=2,n-1
        do i=1,l
          fk      = 1./(cm(i,k)-cl(i,k)*au(i,k-1))
          au(i,k) = fk*cu(i,k)
          a1(i,k) = fk*(r1(i,k)-cl(i,k)*a1(i,k-1))
          a2(i,k) = fk*(r2(i,k)-cl(i,k)*a2(i,k-1))
        enddo
      enddo
      do i=1,l
        fk      = 1./(cm(i,n)-cl(i,n)*au(i,n-1))
        a1(i,n) = fk*(r1(i,n)-cl(i,n)*a1(i,n-1))
        a2(i,n) = fk*(r2(i,n)-cl(i,n)*a2(i,n-1))
      enddo
      do k=n-1,1,-1
        do i=1,l
          a1(i,k) = a1(i,k)-au(i,k)*a1(i,k+1)
          a2(i,k) = a2(i,k)-au(i,k)*a2(i,k+1)
        enddo
      enddo
!-----------------------------------------------------------------------
      return
      end
!fpp$ noconcur r
!-----------------------------------------------------------------------
      subroutine tridin(l,n,nt,cl,cm,cu,r1,r2,au,a1,a2)
!sela %include dbtridi2;
!c
!     use machine     , only : kind_phys
      implicit none
      integer             is,k,kk,n,nt,l,i
      real fk(l)
!c
      real cl(l,2:n), cm(l,n), cu(l,n-1),                               &
                           r1(l,n),   r2(l,n*nt),                       &
                           au(l,n-1), a1(l,n), a2(l,n*nt),              &
                           fkk(l,2:n-1)
!-----------------------------------------------------------------------
      do i=1,l
        fk(i)   = 1./cm(i,1)
        au(i,1) = fk(i)*cu(i,1)
        a1(i,1) = fk(i)*r1(i,1)
      enddo
      do k = 1, nt
        is = (k-1) * n
        do i = 1, l
          a2(i,1+is) = fk(i) * r2(i,1+is)
        enddo
      enddo
      do k=2,n-1
        do i=1,l
          fkk(i,k) = 1./(cm(i,k)-cl(i,k)*au(i,k-1))
          au(i,k)  = fkk(i,k)*cu(i,k)
          a1(i,k)  = fkk(i,k)*(r1(i,k)-cl(i,k)*a1(i,k-1))
        enddo
      enddo
      do kk = 1, nt
        is = (kk-1) * n
        do k=2,n-1
          do i=1,l
            a2(i,k+is) = fkk(i,k)*(r2(i,k+is)-cl(i,k)*a2(i,k+is-1))
          enddo
        enddo
      enddo
      do i=1,l
        fk(i)   = 1./(cm(i,n)-cl(i,n)*au(i,n-1))
        a1(i,n) = fk(i)*(r1(i,n)-cl(i,n)*a1(i,n-1))
      enddo
      do k = 1, nt
        is = (k-1) * n
        do i = 1, l
          a2(i,n+is) = fk(i)*(r2(i,n+is)-cl(i,n)*a2(i,n+is-1))
        enddo
      enddo
      do k=n-1,1,-1
        do i=1,l
          a1(i,k) = a1(i,k) - au(i,k)*a1(i,k+1)
        enddo
      enddo
      do kk = 1, nt
        is = (kk-1) * n
        do k=n-1,1,-1
          do i=1,l
            a2(i,k+is) = a2(i,k+is) - au(i,k)*a2(i,k+is+1)
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
      return
      end
