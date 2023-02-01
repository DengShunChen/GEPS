      subroutine shalcnv_sa(im,ix,km,jcap,delt,del,prsl,ps,phil,ql_all,&
           q1,t1,u1,v1,rn,kbot,ktop,kcnv,slimsk,garea,                 &
           dot,ncloud,hpbl,heat,evap,grav,cp,hvap,rd,t0c)
!          dot,ncloud,hpbl,heat,evap,ud_mf,dt_mf,me)
!
!     use machine , only : kind_phys
!     use funcphys , only : fpvs
!     use physcons, grav => con_g, cp => con_cp, hvap => con_hvap
!    &,             rv => con_rv, fv => con_fvirt, t0c => con_t0c
!    &,             rd => con_rd, cvap => con_cvap, cliq => con_cliq
!    &,             eps => con_eps, epsm1 => con_epsm1
      implicit none
!
      integer            im, ix,  km, jcap, ncloud,                     &
                         kbot(ix), ktop(ix), kcnv(ix)
!                        me
      real delt
      real ps(im),     del(ix,km),  prsl(ix,km),                        &
                           ql(ix,km,2),q1(ix,km),   t1(ix,km),          &
                           u1(ix,km),  v1(ix,km),   rcs(im),            &
                           rn(ix),     slimsk(ix),                      &
                           dot(ix,km), phil(ix,km), hpbl(ix),           &
                           heat(ix),   evap(ix)                         &
                          ,ql_all(ix,km)                                &
! hchuang code change mass flux output
      ,                    ud_mf(im,km),dt_mf(im,km)
!
      integer              i,j,indx, jmn, k, kk, latd, lond, km1, n
      integer              kpbl(im)
!
      real                 c0s,     cpoel,   dellat,  delta,            &
                           desdt,   deta,    detad,   dg,               &
                           dh,      dhh,     dlnsig,  dp,               &
                           dq,      dqsdp,   dqsdt,   dt,               &
                           dt2,     dtmax,   dtmin,   dv1h,             &
                           dv1q,    dv2h,    dv2q,    dv1u,             &
                           dv1v,    dv2u,    dv2v,    dv3q,             &
                           dv3h,    dv3u,    dv3v,    clam,             &
                           dz,      dz1,     e1,                        &
                           el2orc,  elocp,   aafac,                     &
                           es,      etah,    h1,      dthk,             &
                           evef,    evfact,  evfactl, fact1,            &
                           fact2,   factor,  tfac,    rho,              &
                           g,       gamma,   pprime,  betaw,            &
                           qlk,     qrch,    qs,      c1,               &
                           rain,    rfact,   shear,   tem1,             &
                           tem2,    terr,    val,     val1,             &
                           val2,    w1,      w1l,     w1s,              &
                           w2,      w2l,     w2s,     w3,               &
                           w3l,     w3s,     w4,      w4l,              &
                           w4s,     tem,     ptem,    ptem1,            &
                           pgcon,   garea,   dxcrt
!
      integer              kb(im), kbcon(im), kbcon1(im),               & 
                           ktcon(im), ktcon1(im),ktconn(im),            &
                           kbm(im), kmax(im)
!
      real                 aa1(im),     cina(im),                       &
                           umean(im),   tauadv(im), gdx(im),            &
                           delhbar(im), delq(im),   delq2(im),          &
                           delqbar(im), delqev(im), deltbar(im),        &
                           deltv(im),   dtconv(im), edt(im),            &
                           pdot(im),    po(im,km),                      &
                           qcond(im),   qevap(im),  hmax(im),           &
                           rntot(im),   vshear(im),                     &
                           xlamud(im),  xmb(im),    xmbmax(im),         &
                           delubar(im), delvbar(im)
!
      real                 cinpcr, cinpcrmx, cinpcrmn
      real                 cinacr, cinacrmx, cinacrmn
      real                 crtlamd,d0
!
!  parameters for updraft velocity calculation
      real                 bet1,    cd1,     f1,      gam1,             &
                           bb1,     bb2 
!c

!c
!phon
      real grav,cp,hvap,rd,t0c
      real rv,fv,cvap,cliq,eps,epsm1,cm
!phon
      parameter(dthk=25.,cm=1.0,d0=0.01)
      parameter(cinpcrmx=180.,cinpcrmn=120.)
      parameter(cinacrmx=-120.,cinacrmn=-80.)
      parameter(crtlamd=3.e-4)
      parameter(dtmax=10800.,dtmin=600.)
      parameter(bet1=1.875,cd1=.506,f1=2.0,gam1=.5)
      parameter(dxcrt=15.e3)
      parameter(terr=0.,c0s=.002,c1=5.e-4)
      parameter(h1=0.33333333)
!  physical parameters
!     parameter(g=grav)
!     parameter(cpoel=cp/hvap,elocp=hvap/cp,
!    &          el2orc=hvap*hvap/(rv*cp))
!     parameter(terr=0.,c0=.002,c1=5.e-4,delta=fv)
!     parameter(fact1=(cvap-cliq)/rv,fact2=hvap/rv-fact1*t0c)
!  local variables and arrays
      real                 pfld(im,km),    to(im,km),     qo(im,km),    &
                           uo(im,km),      vo(im,km),     qeso(im,km)
!  for updraft velocity calculation
      real                 wu2(im,km),     buo(im,km),    drag(im,km)
      real                 wc(im),         scaldfunc(im), sigmagfm(im)
!  cloud water
!     real                 qlko_ktcon(im), dellal(im,km), tvo(im,km),
      real                 qlko_ktcon(im), dellal(im,km),               &
                           dbyo(im,km),    zo(im,km),     xlamue(im,km),&
                           heo(im,km),     heso(im,km),   qrcko(im,km), &
                           dellah(im,km),  dellaq(im,km),               &
                           dellau(im,km),  dellav(im,km), hcko(im,km),  &
                           ucko(im,km),    vcko(im,km),   qcko(im,km),  &
                           eta(im,km),     zi(im,km),     pwo(im,km),   &
                           tx1(im),        c0t(im,km),    sumx(im),     &
                           c0(im)
                           
!
      logical totflg, cnvflg(im), flg(im)
!
      real                 tf, tcr, tcrf
      parameter (tf=233.16, tcr=263.16, tcrf=1.0/(tcr-tf))
!
!fong--- merge qsatq in here
      real vpsat(191),pqs,qqq,temx,t11
      integer ic
      data vpsat/                                                       &
          9.67165e-05,  1.15983e-04,  1.38819e-04,  1.65835e-04,        &
          1.97736e-04,  2.35339e-04,  2.79584e-04,  3.31553e-04,        &
          3.92489e-04,  4.63820e-04,  5.47177e-04,  6.44430e-04,        &
          7.57710e-04,  8.89450e-04,  1.04242e-03,  1.21975e-03,        &
          1.42503e-03,  1.66230e-03,  1.93614e-03,  2.25172e-03,        &
          2.61488e-03,  3.03222e-03,  3.51113e-03,  4.05995e-03,        &
          4.68804e-03,  5.40589e-03,  6.22523e-03,  7.15922e-03,        &
          8.22253e-03,  9.43153e-03,  1.08045e-02,  1.23617e-02,        &
          1.41258e-02,  1.61219e-02,  1.83779e-02,  2.09244e-02,        &
          2.37959e-02,  2.70300e-02,  3.06684e-02,  3.47573e-02,        &
          3.93475e-02,  4.44947e-02,  5.02607e-02,  5.67130e-02,        &
          6.39258e-02,  7.19807e-02,  8.09670e-02,  9.09823e-02,        &
          1.02134e-01,  1.14538e-01,  1.28323e-01,                      &
!
                        1.45280e-01,  1.64189e-01,  1.85241e-01,        &
          2.08643e-01,  2.34615e-01,  2.63398e-01,  2.95248e-01,        &
          3.30441e-01,  3.69270e-01,  4.12053e-01,  4.59124e-01,        &
          5.10843e-01,  5.67591e-01,  6.29773e-01,  6.97819e-01,        &
          7.72185e-01,  8.53352e-01,  9.41827e-01,  1.03814e+00,        &
          1.14287e+00,  1.25659e+00,  1.37992e+00,  1.51352e+00,        &
          1.65806e+00,  1.81424e+00,  1.98279e+00,  2.16447e+00,        &
          2.36006e+00,  2.57039e+00,  2.79628e+00,  3.03858e+00,        &
          3.29819e+00,  3.57599e+00,  3.87289e+00,  4.18982e+00,        &
          4.52773e+00,  4.88753e+00,  5.27019e+00,  5.67664e+00,        &
!
          6.10780e+00,  6.56617e+00,  7.05475e+00,  7.57526e+00,        &
          8.12946e+00,  8.71922e+00,  9.34647e+00,  1.00132e+01,        &
          1.07216e+01,  1.14739e+01,  1.22723e+01,  1.31192e+01,        &
          1.40172e+01,  1.49688e+01,  1.59767e+01,  1.70438e+01,        &
          1.81729e+01,  1.93672e+01,  2.06298e+01,  2.19639e+01,        &
          2.33729e+01,  2.48605e+01,  2.64302e+01,  2.80858e+01,        &
          2.98314e+01,  3.16708e+01,  3.36085e+01,  3.56487e+01,        &
          3.77959e+01,  4.00548e+01,  4.24303e+01,  4.49274e+01,        &
          4.75511e+01,  5.03069e+01,  5.32001e+01,  5.62365e+01,        &
          5.94220e+01,  6.27625e+01,  6.62643e+01,  6.99337e+01,        &
!
          7.37774e+01,  7.78022e+01,  8.20150e+01,  8.64231e+01,        &
          9.10338e+01,  9.58548e+01,  1.00894e+02,  1.06159e+02,        &
          1.11659e+02,  1.17401e+02,  1.23395e+02,  1.29650e+02,        &
          1.36174e+02,  1.42978e+02,  1.50070e+02,  1.57461e+02,        &
          1.65161e+02,  1.73180e+02,  1.81529e+02,  1.90218e+02,        &
          1.99260e+02,  2.08665e+02,  2.18446e+02,  2.28613e+02,        &
          2.39180e+02,  2.50159e+02,  2.61562e+02,  2.73404e+02,        &
          2.85696e+02,  2.98453e+02,  3.11689e+02,  3.25418e+02,        &
          3.39655e+02,  3.54414e+02,  3.69711e+02,  3.85560e+02,        &
          4.01979e+02,  4.18982e+02,  4.36586e+02,  4.54808e+02,        &
!
          4.73665e+02,  4.93175e+02,  5.13354e+02,  5.34221e+02,        &
          5.55795e+02,  5.78093e+02,  6.01135e+02,  6.24940e+02,        &
          6.49527e+02,  6.74918e+02,  7.01131e+02,  7.28188e+02,        &
          7.56110e+02,  7.84918e+02,  8.14633e+02,  8.45278e+02,        &
          8.76876e+02,  9.09448e+02,  9.43018e+02,  9.77609e+02,        &
          1.01325e+03/

!-----------------------------------------------------------------------
!c--------------------------------------------------------------------
      g=grav
      rv=4.615e+2
      fv=rv/rd-1.
      cvap=1.846e+3
      cliq=4.1855e+3
      eps=rd/rv
      epsm1=eps-1
      cpoel=cp/hvap
      elocp=hvap/cp
      el2orc=hvap*hvap/(rv*cp)
      delta=fv
      fact1=(cvap-cliq)/rv
      fact2=hvap/rv-fact1*t0c
!oc--------------------------------
      temx = 1./1.622

!
      km1 = km - 1
!
!  initialize arrays
!
      do i=1,im
        cnvflg(i) = .true.
        if(kcnv(i).eq.1) cnvflg(i) = .false.
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
        gdx(i) = sqrt(garea)
!
        rcs(i)=1.
      enddo
! hchuang code change
      do k = 1, km
        do i = 1, im
          ud_mf(i,k) = 0.
          dt_mf(i,k) = 0.
        enddo
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!
!>  - determine aerosol-aware rain conversion parameter over land
      do i=1,im
        if(slimsk(i) == 1.) then
           c0(i) = c0s*0.958
        else
           c0(i) = c0s
        endif
      enddo
!
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
!!
!
      dt2   = delt
!!      val   =         1200.
!!      dtmin = max(dt2, val )
!!      val   =         3600.
!!      dtmax = max(dt2, val )
!  model tunable parameters are all here
!byl      clam    = .3
      clam    = .3
      aafac   = .1
      betaw   = .03
!     evef    = 0.07
      evfact  = 0.3
      evfactl = 0.3
!
!     pgcon   = 0.7     ! Gregory et al. (1997, QJRMS)
      pgcon   = 0.55    ! Zhang & Wu (2003,JAS)
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
      do i=1,im
        kbm(i)   = km
        kmax(i)  = km
        tx1(i)   = 1.0 / ps(i)
      enddo
!
      do k = 1, km
        do i=1,im
          if (prsl(i,k)*tx1(i) .gt. 0.70) kbm(i)   = k + 1
          if (prsl(i,k)*tx1(i) .gt. 0.60) kmax(i)  = k + 1
        enddo
      enddo
      do i=1,im
        kbm(i)   = min(kbm(i),kmax(i))
      enddo
!
!  hydrostatic height assume zero terr and compute
!  updraft entrainment rate as an inverse function of height
!
      do k = 1, km
        do i=1,im
          zo(i,k) = phil(i,k) / g
        enddo
      enddo
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
      do i=1,im
        flg(i) = cnvflg(i)
        kpbl(i)= 1
      enddo
      do k = 2, km1
        do i=1,im
          if (flg(i).and.zo(i,k).le.hpbl(i)) then
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   convert surface pressure to mb from cb
!
      do k = 1, km
        do i = 1, im
          if (cnvflg(i) .and. k .le. kmax(i)) then
            pfld(i,k) = prsl(i,k) * 10.0
            eta(i,k)  = 1.
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
            uo(i,k)   = u1(i,k) * rcs(i)
            vo(i,k)   = v1(i,k) * rcs(i)
            wu2(i,k)  = 0.
            buo(i,k)  = 0.
            drag(i,k) = 0.
          endif
        enddo
      enddo
!
!  column variables
!  p is pressure of the layer (mb)
!  t is temperature at t-dt (k)..tn
!  q is mixing ratio at t-dt (kg/kg)..qn
!  to is temperature at t+dt (k)... this is after advection and turbulan
!  qo is mixing ratio at t+dt (kg/kg)..q1
!
      do k = 1, km
        do i=1,im
          if (cnvflg(i) .and. k .le. kmax(i)) then
!         qeso(i,k) = 0.01 * fpvs(to(i,k))      ! fpvs is in pa
!         qeso(i,k) = eps * qeso(i,k) / (pfld(i,k) + epsm1*qeso(i,k))
!cwb qsat
          t11 = max(1.00001, min(190.999, to(i,k)-182.16))
          ic = int(t11)
          pqs = pfld(i,k)
          qqq = min(temx*pqs, vpsat(ic)+(vpsat(1+ic)-vpsat(ic)) &
                                     *(t11-float(ic)))
          qeso(i,k) = 0.622*qqq/(pqs+epsm1*qqq)
!cwb
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
      do k = 1, km
        do i=1,im
          if (cnvflg(i) .and. k .le. kmax(i)) then
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
      do i=1,im
         if (cnvflg(i)) then
            hmax(i) = heo(i,1)
            kb(i) = 1
         endif
      enddo
      do k = 2, km
        do i=1,im
          if (cnvflg(i).and.k.le.kpbl(i)) then
            if(heo(i,k).gt.hmax(i)) then
              kb(i)   = k
              hmax(i) = heo(i,k)
            endif
          endif
        enddo
      enddo
!
      do k = 1, km1
        do i=1,im
          if (cnvflg(i) .and. k .le. kmax(i)-1) then
            dz      = .5 * (zo(i,k+1) - zo(i,k))
            dp      = .5 * (pfld(i,k+1) - pfld(i,k))
!           es      = 0.01 * fpvs(to(i,k+1))      ! fpvs is in pa
!cwb qsat
          t11 = max(1.00001, min(190.999, to(i,k+1)-182.16))
          ic = int(t11)
          pqs =pfld(i,k+1)
          es = min(temx*pqs, vpsat(ic)+(vpsat(1+ic)-vpsat(ic)) &
                                     *(t11-float(ic)))
!cwb
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
      do k = 1, km1
        do i=1,im
          if (cnvflg(i) .and. k .le. kmax(i)-1) then
!           qeso(i,k) = 0.01 * fpvs(to(i,k))      ! fpvs is in pa
!           qeso(i,k) = eps * qeso(i,k) / (po(i,k) + epsm1*qeso(i,k))
!cwb qsat
          t11 = max(1.00001, min(190.999, to(i,k)-182.16))
          ic = int(t11)
          pqs =po(i,k)
          qqq = min(temx*pqs, vpsat(ic)+(vpsat(1+ic)-vpsat(ic)) &
                                     *(t11-float(ic)))
          qeso(i,k) = 0.622*qqq/(pqs+epsm1*qqq)
!cwb
            val1      =             1.e-8
            qeso(i,k) = max(qeso(i,k), val1)
            val2      =           1.e-10
            qo(i,k)   = max(qo(i,k), val2 )
!           qo(i,k)   = min(qo(i,k),qeso(i,k))
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
      do i=1,im
        flg(i)   = cnvflg(i)
        if(flg(i)) kbcon(i) = kmax(i)
      enddo
      do k = 2, km1
        do i=1,im
          if (flg(i).and.k.lt.kbm(i)) then
            if(k.gt.kb(i).and.heo(i,kb(i)).gt.heso(i,k)) then
              kbcon(i) = k
              flg(i)   = .false.
            endif
          endif
        enddo
      enddo
!
      do i=1,im
        if(cnvflg(i)) then
          if(kbcon(i).eq.kmax(i)) cnvflg(i) = .false.
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
!  determine critical convective inhibition
!  as a function of vertical velocity at cloud base.
!
      do i=1,im
        if(cnvflg(i)) then
          pdot(i)  = 10.* dot(i,kbcon(i))
        endif
      enddo
      do i=1,im
        if(cnvflg(i)) then
          if(slimsk(i).eq.1.) then
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
          if(pdot(i).le.w4) then
            ptem = (pdot(i) - w4) / (w3 - w4)
          elseif(pdot(i).ge.-w4) then
            ptem = - (pdot(i) + w4) / (w4 - w3)
          else
            ptem = 0.
          endif
          val1    =             -1.
          ptem = max(ptem,val1)
          val2    =             1.
          ptem = min(ptem,val2)
          ptem = 1. - ptem
          ptem1= .5*(cinpcrmx-cinpcrmn)
          cinpcr = cinpcrmx - ptem * ptem1
          tem1 = pfld(i,kb(i)) - pfld(i,kbcon(i))
          if(tem1.gt.cinpcr) then
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
!  assume the detrainment rate for the updrafts to be same as
!  the entrainment rate at cloud base
!
      do i = 1, im
        if(cnvflg(i)) then
          xlamud(i) = xlamue(i,kbcon(i))
        endif
      enddo
!
!  determine updraft mass flux for the subcloud layers
!
      do k = km1, 1, -1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.lt.kbcon(i).and.k.ge.kb(i)) then
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
           if(k.gt.kbcon(i).and.k.lt.kmax(i)) then
              dz       = zi(i,k) - zi(i,k-1)
              ptem     = 0.5*(xlamue(i,k)+xlamue(i,k-1))-xlamud(i)
              eta(i,k) = eta(i,k-1) * (1 + ptem * dz)
              if(eta(i,k) .le. 0.) then
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
      do i = 1, im
        if(cnvflg(i)) then
          indx         = kb(i)
          hcko(i,indx) = heo(i,indx)
          ucko(i,indx) = uo(i,indx)
          vcko(i,indx) = vo(i,indx)
        endif
      enddo
!
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.lt.kmax(i)) then
              dz   = zi(i,k) - zi(i,k-1)
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              hcko(i,k) = ((1.-tem1)*hcko(i,k-1)+tem*0.5* &
                           (heo(i,k)+heo(i,k-1)))/factor
              dbyo(i,k) = hcko(i,k) - heso(i,k)

              tem  = 0.5 * cm * tem
              factor = 1. + tem
              ptem = tem + pgcon
              ptem1= tem - pgcon
              ucko(i,k) = ((1.-tem1)*ucko(i,k-1)+ptem*uo(i,k) &
                           +ptem1*uo(i,k-1))/factor
              vcko(i,k) = ((1.-tem1)*vcko(i,k-1)+ptem*vo(i,k) &
                           +ptem1*vo(i,k-1))/factor
            endif
          endif
        enddo
      enddo
!
!   taking account into convection inhibition due to existence of
!    dry layers below cloud base
!
      do i=1,im
        flg(i) = cnvflg(i)
        kbcon1(i) = kmax(i)
      enddo
      do k = 2, km1
      do i=1,im
        if (flg(i).and.k.lt.kbm(i)) then
          if(k.ge.kbcon(i).and.dbyo(i,k).gt.0.) then
            kbcon1(i) = k
            flg(i)    = .false.
          endif
        endif
      enddo
      enddo
      do i=1,im
        if(cnvflg(i)) then
          if(kbcon1(i).eq.kmax(i)) cnvflg(i) = .false.
        endif
      enddo
      do i=1,im
        if(cnvflg(i)) then
          tem = pfld(i,kbcon(i)) - pfld(i,kbcon1(i))
          if(tem.gt.dthk) then
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
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k > kb(i) .and. k < kbcon1(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              rfact =  1. + delta * cp * gamma                   &
                      * to(i,k) / hvap
              cina(i) = cina(i) +                                &
!                      dz1 * eta(i,k) * (g / (cp * to(i,k)))     &
                       dz1 * (g / (cp * to(i,k)))                &
                       * dbyo(i,k) / (1. + gamma)                &
                       * rfact
              val = 0.
              cina(i) = cina(i) +                                &
!                      dz1 * eta(i,k) * g * delta *              &
                       dz1 * g * delta *                         &
                       max(val,(qeso(i,k) - qo(i,k)))
            endif
          endif
        enddo
      enddo
!
      do i = 1, im
        if(cnvflg(i)) then
!
          if(slimsk(i) == 1) then
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
!    limited to the level of sigma=0.7
!
      do i = 1, im
        flg(i) = cnvflg(i)
        if(flg(i)) ktcon(i) = kbm(i)
      enddo
      do k = 2, km1
      do i=1,im
        if (flg(i).and.k .lt. kbm(i)) then
          if(k.gt.kbcon1(i).and.dbyo(i,k).lt.0.) then
             ktcon(i) = k
             flg(i)   = .false.
          endif
        endif
      enddo
      enddo
!
!  turn off shallow convection if cloud top is less than pbl top
!
!     do i=1,im
!       if(cnvflg(i)) then
!         kk = kpbl(i)+1
!         if(ktcon(i).le.kk) cnvflg(i) = .false.
!       endif
!     enddo
!!
!     totflg = .true.
!     do i = 1, im
!       totflg = totflg .and. (.not. cnvflg(i))
!     enddo
!     if(totflg) return
!!
!
!  specify upper limit of mass flux at cloud base
!
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
      do i = 1, im
        if (cnvflg(i)) then
          aa1(i) = 0.
          qcko(i,kb(i)) = qo(i,kb(i))
          qrcko(i,kb(i)) = qo(i,kb(i))
        endif
      enddo
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.lt.ktcon(i)) then
              dz    = zi(i,k) - zi(i,k-1)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              qrch = qeso(i,k) &
                   + gamma * dbyo(i,k) / (hvap * (1. + gamma))
!j
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              qcko(i,k) = ((1.-tem1)*qcko(i,k-1)+tem*0.5* &
                           (qo(i,k)+qo(i,k-1)))/factor
              qrcko(i,k) = qcko(i,k)
!j
              dq = eta(i,k) * (qcko(i,k) - qrch)
!
!             rhbar(i) = rhbar(i) + qo(i,k) / qeso(i,k)
!
!  below lfc check if there is excess moisture to release latent heat
!
              if(k.ge.kbcon(i).and.dq.gt.0.) then
                etah = .5 * (eta(i,k) + eta(i,k-1))
                dp = 1000. * del(i,k)
                if(ncloud.gt.0.) then
                  ptem = c0t(i,k) + c1
                  qlk = dq / (eta(i,k) + etah * ptem * dz)
                  dellal(i,k) = etah * c1 * dz * qlk * g / dp
                else
                  qlk = dq / (eta(i,k) + etah * c0t(i,k) * dz)
                endif
                buo(i,k) = buo(i,k) - g * qlk
!                aa1(i) = aa1(i) - dz * g * qlk
                qcko(i,k)= qlk + qrch
                pwo(i,k) = etah * c0t(i,k) * dz * qlk
              endif
!
!  compute buoyancy and drag for updraft velocity
!
              if(k .ge. kbcon(i)) then
                rfact =  1. + delta * cp * gamma                 &
                         * to(i,k) / hvap
                buo(i,k) = buo(i,k) + (g / (cp * to(i,k)))       &
                         * dbyo(i,k) / (1. + gamma)              &
                         * rfact
                val = 0.
                buo(i,k) = buo(i,k) + g * delta *                &
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
!!      do k = 2, km1
!!        do i = 1, im
!!          if (cnvflg(i)) then
!!            if(k.ge.kbcon(i).and.k.lt.ktcon(i)) then
!!              dz1 = zo(i,k+1) - zo(i,k)
!!              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
!!              rfact =  1. + delta * cp * gamma &
!!                       * to(i,k) / hvap
!!              aa1(i) = aa1(i) +                        &
!!                       dz1 * (g / (cp * to(i,k)))      &
!!                       * dbyo(i,k) / (1. + gamma)      &
!!                       * rfact
!!              val = 0.
!!              aa1(i)=aa1(i)+                           &
!!                       dz1 * g * delta *               &
!!                       max(val,(qeso(i,k) - qo(i,k)))
!!            endif
!!          endif
!!        enddo
!!      enddo
!!
      do i = 1, im
        if (cnvflg(i)) then
          aa1(i) = 0.
        endif
      enddo
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k .ge. kbcon(i) .and. k .lt. ktcon(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              aa1(i) = aa1(i) + buo(i,k) * dz1
            endif
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i) .and. aa1(i) .le. 0.) cnvflg(i) = .false.
      enddo
!
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
!    limited to the level of sigma=0.7
!
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
            if(k.ge.ktcon(i).and.k.lt.kbm(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              rfact =  1. + delta * cp * gamma    &
                       * to(i,k) / hvap
              aa1(i) = aa1(i) +                   &
                       dz1 * (g / (cp * to(i,k))) &
                       * dbyo(i,k) / (1. + gamma) &
                       * rfact
              if(aa1(i).lt.0.) then
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
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.ge.ktcon(i).and.k.lt.ktcon1(i)) then
              dz    = zi(i,k) - zi(i,k-1)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              qrch = qeso(i,k) + gamma * dbyo(i,k) / (hvap * (1. + gamma))
!j
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              qcko(i,k) = ((1.-tem1)*qcko(i,k-1)+tem*0.5*  &
                           (qo(i,k)+qo(i,k-1)))/factor
              qrcko(i,k) = qcko(i,k)
!j
              dq = eta(i,k) * (qcko(i,k) - qrch)
!
!  check if there is excess moisture to release latent heat
!
              if(dq.gt.0.) then
                etah = .5 * (eta(i,k) + eta(i,k-1))
                dp = 1000. * del(i,k)
                if(ncloud.gt.0.) then
                  ptem = c0t(i,k) + c1
                  qlk = dq / (eta(i,k) + etah * ptem * dz)
                  dellal(i,k) = etah * c1 * dz * qlk * g / dp
                else
                  qlk = dq / (eta(i,k) + etah * c0t(i,k) * dz)
                endif
                qcko(i,k) = qlk + qrch
                pwo(i,k) = etah * c0t(i,k) * dz * qlk
              endif
            endif
          endif
        enddo
      enddo
!
!  compute updraft velocity square(wu2)
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
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k .gt. kbcon1(i) .and. k .lt. ktcon(i)) then
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
      do i = 1, im
        wc(i) = 0.
        sumx(i) = 0.
      enddo
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k .gt. kbcon1(i) .and. k .lt. ktcon(i)) then
              dz = zi(i,k) - zi(i,k-1)
              tem = 0.5 * (sqrt(wu2(i,k)) + sqrt(wu2(i,k-1)))
              wc(i) = wc(i) + tem * dz
              sumx(i) = sumx(i) + dz
            endif
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i)) then
          if(sumx(i) .eq. 0.) then
             cnvflg(i)=.false.
          else
             wc(i) = wc(i) / sumx(i)
          endif
          val = 1.e-4
          if (wc(i) .lt. val) cnvflg(i)=.false.
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
      if(ncloud.gt.0) then
!
!  compute liquid and vapor separation at cloud top
!
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
          if(dq.gt.0.) then
            qlko_ktcon(i) = dq
            qcko(i,k) = qrch
          endif
        endif
      enddo
      endif
!
!--- compute precipitation efficiency in terms of windshear
!
      do i = 1, im
        if(cnvflg(i)) then
          vshear(i) = 0.
        endif
      enddo
      do k = 2, km
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.le.ktcon(i)) then
              shear= sqrt((uo(i,k)-uo(i,k-1)) ** 2  &
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
!--- what would the change be, that a cloud with unit mass
!--- will do to the environment?
!
      do k = 1, km
        do i = 1, im
          if(cnvflg(i) .and. k .le. kmax(i)) then
            dellah(i,k) = 0.
            dellaq(i,k) = 0.
            dellau(i,k) = 0.
            dellav(i,k) = 0.
          endif
        enddo
      enddo
!
!--- changed due to subsidence and entrainment
!
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.lt.ktcon(i)) then
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
! 
              dellah(i,k) = dellah(i,k) +                               &
           ( eta(i,k)*dv1h - eta(i,k-1)*dv3h                            &
          -  tem*eta(i,k-1)*dv2h*dz                                     &
          +  tem1*eta(i,k-1)*.5*(hcko(i,k)+hcko(i,k-1))*dz              &
               ) *g/dp
! 
              dellaq(i,k) = dellaq(i,k) +                               &
           ( eta(i,k)*dv1q - eta(i,k-1)*dv3q                            &
          -  tem*eta(i,k-1)*dv2q*dz                                     &
          +  tem1*eta(i,k-1)*.5*(qcko(i,k)+qcko(i,k-1))*dz              &
               ) *g/dp
! 
              tem1=eta(i,k)*(uo(i,k)-ucko(i,k))
              tem2=eta(i,k-1)*(uo(i,k-1)-ucko(i,k-1))
              dellau(i,k) = dellau(i,k) + (tem1-tem2) * g/dp
!
              tem1=eta(i,k)*(vo(i,k)-vcko(i,k))
              tem2=eta(i,k-1)*(vo(i,k-1)-vcko(i,k-1))
              dellav(i,k) = dellav(i,k) + (tem1-tem2) * g/dp
!
            endif
          endif
        enddo
      enddo
!
!------- cloud top
!
      do i = 1, im
        if(cnvflg(i)) then
          indx = ktcon(i)
          dp = 1000. * del(i,indx)
          dv1h = heo(i,indx-1)
          dellah(i,indx) = eta(i,indx-1) * (hcko(i,indx-1) - dv1h) * g / dp
          dv1q = qo(i,indx-1)
          dellaq(i,indx) = eta(i,indx-1) * (qcko(i,indx-1) - dv1q) * g / dp
          dv1u = uo(i,indx-1)
          dellau(i,indx) = eta(i,indx-1) * (ucko(i,indx-1) - dv1u) * g / dp
          dv1v = vo(i,indx-1)
          dellav(i,indx) = eta(i,indx-1) * (vcko(i,indx-1) - dv1v) * g / dp
!
!  cloud water
!
          dellal(i,indx) = eta(i,indx-1) * qlko_ktcon(i) * g / dp
        endif
      enddo
!
!
!  compute convective turn-over time
!
      do i= 1, im
        if(cnvflg(i)) then
          tem = zi(i,ktcon1(i)) - zi(i,kbcon1(i))
          dtconv(i) = tem / wc(i)
          tfac = 1. + gdx(i) / 75000.
          dtconv(i) = tfac * dtconv(i)
          dtconv(i) = max(dtconv(i),dtmin)
          dtconv(i) = max(dtconv(i),dt2)
          dtconv(i) = min(dtconv(i),dtmax)
        endif
      enddo
!
!> - Calculate advective time scale (tauadv) using a mean cloud layer
!wind speed.
      do i= 1, im
        if(cnvflg(i)) then
          sumx(i) = 0.
          umean(i) = 0.
        endif
      enddo
      do k = 2, km1
        do i = 1, im
          if(cnvflg(i)) then
            if(k .ge. kbcon1(i) .and. k .lt. ktcon1(i)) then
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
!  mass flux at cloud base for shallow convection
!  (Grant, 2001)
!
!!      do i= 1, im
!!        if(cnvflg(i)) then
!!          k = kbcon(i)
!!!         ptem = g*sflx(i)*zi(i,k)/t1(i,1)
!!          ptem = g*sflx(i)*hpbl(i)/t1(i,1)
!!          wstar(i) = ptem**h1
!!          tem = po(i,k)*100. / (rd*t1(i,k))
!!          xmb(i) = betaw*tem*wstar(i)
!!          xmb(i) = min(xmb(i),xmbmax(i))
!!        endif
!!      enddo
!
!
!> - For scale-aware parameterization, the updraft fraction (sigmagfm)
!is first computed as a function of the lateral entrainment rate at
!cloud base (see Han et al.'s (2017) \cite han_et_al_2017 equation 4 and
!5), following the study by Grell and Freitas (2014) \cite
!grell_and_freitus_2014.
      do i = 1, im
        if(cnvflg(i)) then
          tem = min(max(xlamue(i,kbcon(i)), 2.e-4), 6.e-4)
          tem = 0.2 / tem
          tem1 = 3.14 * tem * tem
          sigmagfm(i) = tem1 / garea
          sigmagfm(i) = max(sigmagfm(i), 0.001)
          sigmagfm(i) = min(sigmagfm(i), 0.999)
        endif
      enddo
!
!> - Then, calculate the reduction factor (scaldfunc) of the vertical
!convective eddy transport of mass flux as a function of updraft
!fraction from the studies by Arakawa and Wu (2013) \cite
!arakawa_and_wu_2013 (also see Han et al.'s (2017) \cite han_et_al_2017
!equation 1 and 2). The final cloud base mass flux with scale-aware
!parameterization is obtained from the mass flux when sigmagfm << 1,
!multiplied by the reduction factor (Han et al.'s (2017) \cite
!han_et_al_2017 equation 2).
      do i = 1, im
        if(cnvflg(i)) then
          if (gdx(i) < dxcrt) then
            scaldfunc(i) = (1.-sigmagfm(i)) * (1.-sigmagfm(i))
            scaldfunc(i) = max(min(scaldfunc(i), 1.0), 0.)
          else
            scaldfunc(i) = 1.0
          endif
          xmb(i) = xmb(i) * scaldfunc(i)
          xmb(i) = min(xmb(i),xmbmax(i))
        endif
      enddo

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      do k = 1, km
        do i = 1, im
          if (cnvflg(i) .and. k .le. kmax(i)) then
!           qeso(i,k) = 0.01 * fpvs(t1(i,k))      ! fpvs is in pa
!           qeso(i,k) = eps * qeso(i,k) / (pfld(i,k) + epsm1*qeso(i,k))
!cwb qsat
          t11 = max(1.00001, min(190.999, t1(i,k)-182.16))
          ic = int(t11)
          pqs =pfld(i,k)
          qqq = min(temx*pqs, vpsat(ic)+(vpsat(1+ic)-vpsat(ic)) &
                                     *(t11-float(ic)))
          qeso(i,k) = 0.622*qqq/(pqs+epsm1*qqq)
!cwb
            val     =             1.e-8
            qeso(i,k) = max(qeso(i,k), val )
          endif
        enddo
      enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
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
            if(k.gt.kb(i).and.k.le.ktcon(i)) then
              dellat = (dellah(i,k) - hvap * dellaq(i,k)) / cp
              t1(i,k) = t1(i,k) + dellat * xmb(i) * dt2
              q1(i,k) = q1(i,k) + dellaq(i,k) * xmb(i) * dt2
!              tem = 1./rcs(i)
!              u1(i,k) = u1(i,k) + dellau(i,k) * xmb(i) * dt2 * tem
!              v1(i,k) = v1(i,k) + dellav(i,k) * xmb(i) * dt2 * tem
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
      do k = 1, km
        do i = 1, im
          if (cnvflg(i)) then
            if(k.gt.kb(i).and.k.le.ktcon(i)) then
!             qeso(i,k) = 0.01 * fpvs(t1(i,k))      ! fpvs is in pa
!             qeso(i,k) = eps * qeso(i,k)/(pfld(i,k) + epsm1*qeso(i,k))
!cwb qsat
          t11 = max(1.00001, min(190.999, t1(i,k)-182.16))
          ic = int(t11)
          pqs =pfld(i,k)
          qqq = min(temx*pqs, vpsat(ic)+(vpsat(1+ic)-vpsat(ic)) &
                                     *(t11-float(ic)))
          qeso(i,k) = 0.622*qqq/(pqs+epsm1*qqq)
!cwb
              val     =             1.e-8
              qeso(i,k) = max(qeso(i,k), val )
            endif
          endif
        enddo
      enddo
!
      do i = 1, im
        rntot(i) = 0.
        delqev(i) = 0.
        delq2(i) = 0.
        flg(i) = cnvflg(i)
      enddo
      do k = km, 1, -1
        do i = 1, im
          if (cnvflg(i)) then
            if(k.lt.ktcon(i).and.k.gt.kb(i)) then
              rntot(i) = rntot(i) + pwo(i,k) * xmb(i) * .001 * dt2
            endif
          endif
        enddo
      enddo
!
! evaporating rain
!
      do k = km, 1, -1
        do i = 1, im
          if (k .le. kmax(i)) then
            deltv(i) = 0.
            delq(i) = 0.
            qevap(i) = 0.
            if(cnvflg(i)) then
              if(k.lt.ktcon(i).and.k.gt.kb(i)) then
                rn(i) = rn(i) + pwo(i,k) * xmb(i) * .001 * dt2
              endif
            endif
            if(flg(i).and.k.lt.ktcon(i)) then
              evef = edt(i) * evfact
              if(slimsk(i).eq.1.) evef=edt(i) * evfactl
!             if(slimsk(i).eq.1.) evef=.07
!             if(slimsk(i).ne.1.) evef = 0.
              qcond(i) = evef * (q1(i,k) - qeso(i,k))  &
                       / (1. + el2orc * qeso(i,k) / t1(i,k)**2)
              dp = 1000. * del(i,k)
              if(rn(i).gt.0..and.qcond(i).lt.0.) then
                qevap(i) = -qcond(i) * (1.-exp(-.32*sqrt(dt2*rn(i))))
                qevap(i) = min(qevap(i), rn(i)*1000.*g/dp)
                delq2(i) = delqev(i) + .001 * qevap(i) * dp / g
              endif
              if(rn(i).gt.0..and.qcond(i).lt.0..and.  &
                 delq2(i).gt.rntot(i)) then
                qevap(i) = 1000.* g * (rntot(i) - delqev(i)) / dp
                flg(i) = .false.
              endif
              if(rn(i).gt.0..and.qevap(i).gt.0.) then
                tem  = .001 * dp / g
                tem1 = qevap(i) * tem
                if(tem1.gt.rn(i)) then
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
!              dellaq(i,k) = dellaq(i,k) + delq(i) / xmb(i)
              delqbar(i) = delqbar(i) + delq(i)*dp/g
              deltbar(i) = deltbar(i) + deltv(i)*dp/g
            endif
          endif
        enddo
      enddo
!j
!     do i = 1, im
!     if(me.eq.31.and.cnvflg(i)) then
!     if(cnvflg(i)) then
!       print *, ' shallow delhbar, delqbar, deltbar = ',
!    &             delhbar(i),hvap*delqbar(i),cp*deltbar(i)
!       print *, ' shallow delubar, delvbar = ',delubar(i),delvbar(i)
!       print *, ' precip =', hvap*rn(i)*1000./dt2
!       print*,'pdif= ',pfld(i,kbcon(i))-pfld(i,ktcon(i))
!     endif
!     enddo
!j
      do i = 1, im
        if(cnvflg(i)) then
          if(rn(i).lt.0..or..not.flg(i)) rn(i) = 0.
          ktop(i) = ktcon(i)
          kbot(i) = kbcon(i)
          kcnv(i) = 2
        endif
      enddo
!
!  cloud water
!
      if (ncloud.gt.0) then
!
      do k = 1, km1
        do i = 1, im
          if (cnvflg(i)) then
!            if (k.ge.kb(i).and.k.le.ktcon(i)) then
            if (k.ge.kbcon(i).and.k.le.ktcon(i)) then
              tem  = dellal(i,k) * xmb(i) * dt2
!cwb
              ql_all(i,k)=ql_all(i,k)+tem
!cwb
!             tem1 = max(0.0, min(1.0, (tcr-t1(i,k))*tcrf))
!             if (ql(i,k,2) .gt. -999.0) then
!               ql(i,k,1) = ql(i,k,1) + tem * tem1            ! ice
!               ql(i,k,2) = ql(i,k,2) + tem *(1.0-tem1)       ! water
!             else
!               ql(i,k,1) = ql(i,k,1) + tem
!             endif
            endif
          endif
        enddo
      enddo
!
      endif
!
! hchuang code change
!
      do k = 1, km
        do i = 1, im
          if(cnvflg(i)) then
            if(k.ge.kb(i) .and. k.lt.ktop(i)) then
              ud_mf(i,k) = eta(i,k) * xmb(i) * dt2
            endif
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i)) then
           k = ktop(i)-1
           dt_mf(i,k) = ud_mf(i,k)
        endif
      enddo
!!
      return
      end
