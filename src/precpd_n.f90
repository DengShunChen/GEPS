       subroutine precpd_n (im,ix,km,dt,del,prsl,ps, &
                          q,cwm,qr,qs,phi,t,rn,sr, &
                          rainp,u00k,lprnt)
!xb110>
!                         dzl,psautco,prautco,evpco,wminco,
!     &                    lprnt,jpr,kdt,me)
!xb110<
!e
!
!     ******************************************************************
!     *                                                                *
!     *           subroutine for precipitation processes               *
!     *           from suspended cloud water/ice                       *
!     *                                                                *
!     ******************************************************************
!     *                                                                *
!     *  originally created by  q. zhao                jan. 1995       *
!     *                         -------                                *    
!     *  modified and rewritten by shrinivas moorthi   oct. 1998       *
!     *                            -----------------                   *
!     *  and                       hua-lu pan                          *
!     *                            ----------                          *
!     *                                                                *
!     *  references:                                                   *
!     *                                                                *
!     *  zhao and carr (1997), monthly weather review (august)         *
!     *  sundqvist et al., (1989) monthly weather review. (august)     *
!     *  chuang 2013, modify sr to define frozen precipitation fraction*
!     ******************************************************************
!
!     in this code vertical indexing runs from surface to top of the
!     model
!
!     argument list:
!     --------------
!       im         : inner dimension over which calculation is made
!       ix         : maximum inner dimension
!       km         : number of vertical levels
!       dt         : time step in seconds
!       del(km)    : pressure layer thickness (bottom to top)
!       prsl(km)   : pressure values for model layers (bottom to top)
!       ps(im)     : surface pressure (centibars)
!       q(ix,km)   : specific humidity (updated in the code)
!       cwm(ix,km) : condensate mixing ratio (updated in the code)
!       qr(ix,km) : condensate rain mixing ratio (updated in the code)
!       qs(ix,km) : condensate snow mixing ratio (updated in the code)
!       t(ix,km)   : temperature       (updated in the code)
!       rn(im)     : precipitation over one time-step dt (m/dt)
!old      sr(im)     : index (=-1 snow, =0 rain/snow, =1 rain)
!new    sr(im)     : "snow ratio", ratio of snow to total precipitation
!       cll(ix,km) : cloud cover
!hchuang rn(im) unit in m per time step
!        precipitation rate conversion 1 mm/s = 1 kg/m2/s
!
!xb110>
!      use machine , only : kind_phys
!      use funcphys , only : fpvs
!      use physcons, grav => con_g, hvap => con_hvap, hfus => con_hfus
!     &,             ttp => con_ttp, cp => con_cp
!     &,             eps => con_eps, epsm1 => con_epsm1
!     &,             R   => con_rd 
!xb110<
      implicit none
!xb110>
      real hvap,grav,hfus,ttp,rd,rv,cp,eps,epsm1

      parameter(  hvap = 2.5000e+6 )
      parameter(  grav = 9.80665e+0 )
      parameter(  hfus = 3.3358e+5 )
      parameter(   ttp = 2.7316e+2 )
      parameter(    rd = 2.8705e+2 )
      parameter(    rv = 4.6150e+2 )
      parameter(    cp = 1.0046e+3 )
      parameter(   eps = rd/rv )
      parameter( epsm1 = eps-1 )
!xb110<
!     include 'constant.h'
!
!      real (kind=kind_phys) g,      h1,    h1000    !xb110
      real  g,      h1,    h1000                                  &
      ,                     d00                                   &
      ,                     elwv,   eliv,  row                    &
      ,                     epsq,   eliw                          &
      ,                     rcp,    rrow                          &
!->rsun 
      ,                     vts, vtr, onsr(ix), onss(ix)          
       
      integer               nsr(ix),    nss(ix), nt 
!<-rsun 
       parameter (g=grav,         h1=1.e0,     h1000=1000.0       &
      ,           d00=0.e0                                        &
      ,           elwv=hvap,      eliv=hvap+hfus,   row=1.e3      &
      ,           epsq=2.e-12                                     &
      ,           eliw=eliv-elwv, rcp=h1/cp,   rrow=h1/row)       
!->rsun 
!       parameter( vts = 1.0, vtr = 5.0) 
!       parameter( vts = 1.0, vtr = 4.0) !rsun 
       parameter( vts = 1.0, vtr = 2.0) !xb110
!      real(kind=kind_phys), parameter :: minp=1.e-12
!      real(kind=kind_phys), parameter :: minp=1.e-15 !xb110
      real, parameter :: minp=1.e-15
!      real(kind=kind_phys) watertot1(im), watertot2(im)  !xb110
      real watertot1(im), watertot2(im) 
!<-rsun 
     
!
!      real(kind=kind_phys), parameter :: cons_0=0.0,     cons_p01=0.01    !xb110
      real, parameter :: cons_0=0.0,     cons_p01=0.01                         &
      ,                                  cons_20=20.0                          &
      ,                                  cons_m30=-30.0, cons_50=50.0
!
      integer im, ix, km, jpr,kdt
!      real (kind=kind_phys) q(ix,km),   t(ix,km),    cwm(ix,km)             !xb110
      real                  q(ix,km),   t(ix,km),    cwm(ix,km)          &
      ,                                 del(ix,km),  prsl(ix,km)         &
!    &,                     cll(im,km), del(ix,km),  prsl(ix,km)         &
      ,                     ps(ix),     rn(ix),      sr(ix)              &
      ,                     dt                                           &
!->rsun 
      ,                     qr(ix,km), qs(ix,km), dzl(ix,km),sn(im)      &
!<-rsun 
!hchuang code change [+1l] : add record to record information in vertical in
!                       addition to total column precrl
      ,                     rainp(ix,km), rnp(im),                       &
                            psautco, prautco, evpco, wminco(2)
!
!
!      real (kind=kind_phys) err(im),      ers(im) ! rsun,     precrl(im)        !xb110
      real  err(im),      ers(im)                                        & ! rsun,     precrl(im)                 
!rsun      &,                     precsl(im),   precrl1(im),!precsl1(im) &
      ,                     rq(im),       condt(im)                       &
      ,                     conde(im),    rconde(im),  tmt0(im)           &
      ,                     wmin(im,km),  wmink(im),   pres(im)           &
      ,                     wmini(im,km), ccr(im)                         &
      ,                     tt(im),       qq(im),      ww(im)             &
      ,                     u00k(ix,km)                                   &
      ,                     zaodt                                         &
!->rsun 
      ,                     rho(im,km),seds(im,km),sedr(im,km)            &
      ,                     qqs(im),qqr(im)
!<-rsun 
!       real (kind=kind_phys) cclim(km) !xb110
       real  cclim(km)
!
      integer iw(im,km), ipr(im), iwl(im),     iwl1(im)
!
       logical comput(im)
       logical lprnt
       logical noskip 
!
!      real (kind=kind_phys) ke,   rdt,  us, cclimit, climit, cws, csm1      !xb110
      real  ke,   rdt,  us, cclimit, climit, cws, csm1                      &
      ,                     crs1, crs2, cr, aa2,     dtcp,   c00, cmr       &
      ,                     tem,  c1,   c2, wwn                             &
!    &,                     tem,  c1,   c2, u00b,    u00t,   wwn            &
      ,                     precrk, precsk, pres1,   qk,     qw,  qi        &
      ,                     qint, fiw, wws, cwmk, expf                      &
      ,                     psaut, psaci, amaxcm, tem1, tem2                &
      ,                     tmt0k, psm1, psm2, ppr                          &
      ,                     rprs,  erk,   pps, sid, rid, amaxps             &
      ,                     praut, fi, qc, amaxrq, rqkll                    &
!->rsan 
      ,                     erkdt 
      integer me
!<-rsun 
      integer i, k, ihpr, n

!xb110>
      real      phi(ix,km)

!--- merge qsatq in here
      real vpsat(191),pqs,qqq,temx,t1
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
!
!xb110<
!
!xb110>
      temx = 1.0/1.622
      psautco   = 4.0e-4
      prautco   = 1.0e-4
      evpco     = 2.0e-5
      wminco(1) = 0.5e-5
      wminco(2) = 0.5e-5

      do k = 1,km
       do i = 1,im
         dzl(i,k) = phi(i,k)/g
       end do
      end do
!xb110<
!-----------------------preliminaries ---------------------------------
!
!     do k=1,km
!       do i=1,im
!         cll(i,k) = 0.0
!       enddo
!     enddo
!
      rdt     = h1 / dt
!rsun evpco   = 2.0e-5 
!     ke      = 2.0e-5  ! commented on 09/10/99  -- opr value
!     ke      = 2.0e-6
!     ke      = 1.0e-5
!!!   ke      = 5.0e-5
!!    ke      = 7.0e-5
      ke      = evpco
!     ke      = 7.0e-5
      us      = h1
      cclimit = 1.0e-3
      climit  = 1.0e-20
      cws     = 0.025
!
      zaodt   = 800.0 * rdt
!xb110      zaodt   = 800.0 * rdt*2.0 ! rsun 
!xb110      zaodt   = 800.0 * rdt*4.0 ! rsun 
!
      csm1    = 5.0000e-8   * zaodt
      crs1    = 5.00000e-6  * zaodt
      crs2    = 6.66600e-10 * zaodt
      cr      = 5.0e-4      * zaodt
      aa2     = 1.25e-3     * zaodt
!rsun aa2     = 2.5e-3     * zaodt
!
      ke      = ke * sqrt(rdt)
!     ke      = ke * sqrt(zaodt)
!
      dtcp    = dt * rcp
!
!     c00 = 1.5e-1 * dt
!     c00 = 10.0e-1 * dt
!     c00 = 3.0e-1 * dt          !05/09/2000
!      c00 = 1.0e-4 * dt          !05/09/2000   !xb110
      c00 = prautco * dt         !05/09/2000
      cmr = 1.0 / 3.0e-4
!     cmr = 1.0 / 5.0e-4
!     c1  = 100.0
      c1  = 300.0
      c2  = 0.5
!
!
!--------calculate c0 and cmr using lc at previous step-----------------
!
!xb110>
      prsl = prsl *1000.    ! convert the unit from cb to Pa !xb110
      del  = del  *1000.
!xb110<

      do k=1,km
        do i=1,im
          tem   = (prsl(i,k)*0.00001)  
!         tem   = sqrt(tem)
          iw(i,k)    = 0.0
!         wmin(i,k)  = 1.0e-5 * tem
!         wmini(i,k) = 1.0e-5 * tem       ! testing for ras
!

          wmin(i,k)  = wminco(1) * tem
          wmini(i,k) = wminco(2) * tem

          rainp(i,k) = 0.0
 
          rho(i,k) = 0.622*prsl(i,k)/(rd*t(i,k)*(q(i,k) + 0.622))      

        enddo
      enddo

!       conservation 

      do i = 1, im
        watertot1(i) = 0.0 
        do k = 1, km 
          watertot1(i) = watertot1(i) +                           &
              (q(i,k)+cwm(i,k)+qr(i,k)+qs(i,k))*del(i,k)/grav
        enddo 
      enddo
!->rsun 
      do i=1,im
!       c0(i)  = 1.5e-1
!       cmr(i) = 3.0e-4
!
        iwl1(i)    = 0
! - - - - - - - - - - - 
!rsun        precrl1(i) = d00
!rsun        precsl1(i) = d00
! - - - - - - - - - - - 
        comput(i)  = .false.
        rn(i)      = d00
!->rsun
        sn(i)      = d00 
!<-rsun 
        sr(i)      = d00
        ccr(i)     = d00
!
        rnp(i)     = d00
      enddo
!------------select columns where precip falls --------------
!      do i = 1, im 
!        tem = dzl(i,1)/vtr
!        nsr(i)  = int(dt/tem+1.)
!        tem = dzl(i,1)/vts
!        nss(i)  = int(dt/tem+1.)
!        onsr(i) = 1./nsr(i)
!        onss(i) = 1./nss(i)
!      enddo 
!      
!      if(me == 0 ) then 
!!        print*,'dzl(1,1)',kdt,dzl(1,1)
!!        print*,'nsr,nss,onsr, nnss:',nsr,nss,onsr, onss
!!        print*,'precpd before seds : qs:',kdt,qs
!      endif 
!
!      do i = 1, im 
!        do nt = 1, nss (i)
!          do k = km, 1, -1 
!            if(qs(i,k) .gt. minp) then 
!              seds(i,k)  = vts * qs(i,k) * rho(i,k) 
!            else 
!              seds(i,k)  = 0.0 
!            endif 
!          enddo 
!!rsun          qs(i,km) = qs(i,km)-seds(i,km)/dzl(i,km)*dt/nss(i)/rho(i,km) 
!          qs(i,km) = qs(i,km)-seds(i,km)/
!     &               del(i,km)*dt/nss(i)*grav
!          qs(i,km) = max(qs(i,km),0.0)
!          do k = km-1, 1, -1 
!            tem1=qs(i,k+1)
!            tem =qs(i,k)
!!            qs(i,k) = qs(i,k) + 
!!     &              (seds(i,k+1) - seds(i,k))/rho(i,k)*
!!     &              dt*onss(i)/dzl(i,k)  
!            qs(i,k) = qs(i,k) + 
!     &             (seds(i,k+1) - seds(i,k))/del(i,k)*grav*
!     &              dt*onss(i)
!!rsun            qs(i,k) = max(qs(i,k),0.0)
!!            if(qs(i,k) < 0.0) then 
!!               print*,'neg value qs:', 
!!     &                  k,tem1, tem, rho(i,k+1), rho(i,k), 
!!     &                  seds(i,k+1), seds(i,k),
!!     &                  onss(i),dzl(i,k),qs(i,k) 
!! 
!!            endif 
!          enddo 
!          rn(i) = rn(i)+seds(i,1)*dt*onss(i)*rrow
!          sn(i) = sn(i)+seds(i,1)*dt*onss(i)*rrow
!        enddo 
!      enddo 
!  
!!      if(me==0) then 
!!        print*,'precpd after seds: seds:',kdt,seds
!!        print*,'precpd after seds : rn:',kdt,rn
!!        print*,'precpd after seds : sn:',kdt,sn
!!        print*,'precpd after seds : qs:',kdt,qs
!!        print*, 'precpd before sedr : qr:',kdt, qr 
!!      endif 
!
!      do i = 1, im 
!        do nt = 1, nsr(i)
!          do k = km, 1, -1 
!            if(qr(i,k) .gt. minp) then 
!              sedr(i,k)  = vtr * qr(i,k) * rho(i,k) 
!            else 
!              sedr(i,k)  = 0.0 
!            endif 
!          enddo 
!          qr(i,km) = qr(i,km)-sedr(i,km)/
!     &               del(i,km)*dt/nsr(i)*grav
!          do k = km-1, 1, -1 
!            tem = qr(i,k) 
!            tem1 = qr(i,k+1)
!            
!            qr(i,k) = qr(i,k) + 
!     &             (sedr(i,k+1) - sedr(i,k))/del(i,k)*grav*
!     &              dt*onsr(i)
!!rsun            qr(i,k) = max(qr(i,k),0.0)
!!            if(qr(i,k) < 0.0) then 
!!               print*,'neg value qr:', 
!!     &                  k,tem1, tem, rho(i,k+1), rho(i,k), 
!!     &                  sedr(i,k+1), sedr(i,k),
!!     &                  onsr(i),dzl(i,k),qr(i,k) 
!! 
!!            endif 
!          enddo 
!          rn(i) = rn(i)+sedr(i,1)*dt*onsr(i)*rrow
!
!        enddo 
!      enddo 
!
!!      if(me ==0) then 
!!        print*,'precpd after sedr and dt: sedr:',kdt,sedr
!!        print*,'precpd after sedr  dt: rn:',kdt,rn
!!        print*,'precpd after sedr dt : sn:',kdt,sn
!!        print*,'precpd after sedr dt : qr:',kdt,qr
!!      endif 
!
!      do i = 1, im 
!        rid = rn(i)
!        if (rid < 1.e-13) then
!           sr(i) = 0.
!        else
!           sr(i) = sn(i)/rid
!        endif
!      enddo 
!!->rsun conservation 
!
!      do i = 1, im
!        watertot2(i) = 0.0
!        do k = 1, km
!          watertot2(i) = watertot2(i) +
!     &        (q(i,k)+cwm(i,k)+qr(i,k)+qs(i,k))*del(i,k)/grav
!        enddo
!        watertot2(i) = watertot2(i) + rn(i)*row 
!
!        if(abs(watertot2(i) - watertot1(i)) > 1.e-6) then 
!          print*,'different1:',kdt, watertot1(i), watertot2(i), rn(i),
!     &        abs(watertot2(i) - watertot1(i))
!        endif 
!      enddo

!------------select columns where rain can be produced--------------
      do k=1, km-1
        do i=1,im
          tem = min(wmin(i,k), wmini(i,k))
          if (cwm(i,k) > tem) comput(i) = .true.
        enddo
      enddo
      ihpr = 0
      do i=1,im
        if (comput(i)) then
           ihpr      = ihpr + 1
           ipr(ihpr) = i
        endif
      enddo
!***********************************************************************
!-----------------begining of precipitation calculation-----------------
!***********************************************************************
!     do k=km-1,2,-1
      do k=km,1,-1
        do n=1,ihpr
!rsun          precrl(n) = precrl1(n)
!rsun          precsl(n) = precsl1(n)
          err  (n)  = d00
          ers  (n)  = d00
          iwl  (n)  = 0
!
          i         = ipr(n)
          tt(n)     = t(i,k)
          qq(n)     = q(i,k)
          ww(n)     = cwm(i,k)
          wmink(n)  = wmin(i,k)
          pres(n)   = prsl(i,k)
!->rsun 
          qqs(n)    = qs(i,k) 
          qqr(n)    = qr(i,k) 
!<-rsun 
!
!rsun     precrk = max(cons_0,    precrl1(n))
          precrk = max(cons_0,    qqr(n) * del(i,k)/grav)
!rsun     precsk = max(cons_0,    precsl1(n))
          precsk = max(cons_0,    qqs(n) * del(i,k)/grav) 
          wwn    = max(ww(n), climit)
!         if (wwn .gt. wmink(n) .or. (precrk+precsk) .gt. d00) then
          if (wwn > climit .or. (precrk+precsk) > d00) then
!rsun     if (wwn > climit .or. (seds(i,k)+sedr(i,k)) > d00) then ! rsun ???
            comput(n) = .true.
          else
            comput(n) = .false.
          endif
        enddo
!
!       es(1:ihpr) = fpvs(tt(1:ihpr))
        do n=1,ihpr
          if (comput(n)) then
            i = ipr(n)
            conde(n)  = (dt/g) * del(i,k) 
            condt(n)  = conde(n) * rdt
            rconde(n) = h1 / conde(n)
            qk        = max(epsq,  qq(n))
            tmt0(n)   = tt(n) - 273.16
            wwn       = max(ww(n), climit)
!
!           pl = pres(n) * 0.01
!           call qsatd(tt(n), pl, qc)
!           rq(n) = max(qq(n), epsq) / max(qc, 1.0e-10)
!           rq(n) = max(1.0e-10, rq(n))           ! -- relative humidity---
!
!  the global qsat computation is done in pa
            pres1   = pres(n) 
!           qw      = es(n)
!xb110>
!            qw      = min(pres1, fpvs(tt(n)))
!            qw      = eps * qw / (pres1 + epsm1 * qw)
!            qw      = max(qw,epsq)
          t1 = max(1.00001, min(190.999, tt(n)-182.16))
          ic = int(t1)
          pqs = pres(n)/100.
          qqq = min(temx*pqs, vpsat(ic)+(vpsat(1+ic)-vpsat(ic))  &
                                     *(t1-float(ic)))
!         qw = 0.622*qqq/(pqs-qqq)
          qw = 0.622*qqq/(pqs+epsm1*qqq)

          qw = max(qw,epsq)
!xb110<
!
!           tmt15 = min(tmt0(n), cons_m15)
!           ai    = 0.008855
!           bi    = 1.0
!           if (tmt0(n) .lt. -20.0) then
!             ai = 0.007225
!             bi = 0.9674
!           endif
!           qi   = qw * (bi + ai*min(tmt0(n),cons_0))
!           qint = qw * (1.-0.00032*tmt15*(tmt15+15.))
!
            qi   = qw
            qint = qw
!           if (tmt0(n).le.-40.) qint = qi
!
!-------------------ice-water id number iw------------------------------
            if(tmt0(n) < -15.) then
!rsun       if(tmt0(n) < -12.) then
!rsun       if(tmt0(n) < -10.) then
               fi = qk - u00k(i,k)*qi
               if(fi > d00 .or. wwn > climit) then
                  iwl(n) = 1
               else
                  iwl(n) = 0
               endif
!           endif
            elseif (tmt0(n) >= 0.) then
               iwl(n) = 0
!
!           if(tmt0(n).lt.0.0.and.tmt0(n).ge.-15.0) then
            else
              iwl(n) = 0
              if(iwl1(n) == 1 .and. wwn > climit) iwl(n) = 1
            endif
!
!           if(tmt0(n).ge.0.) then
!              iwl(n) = 0
!           endif
!----------------the satuation specific humidity------------------------
            fiw   = float(iwl(n))
            qc    = (h1-fiw)*qint + fiw*qi
!----------------the relative humidity----------------------------------
            if(qc <= 1.0e-10) then
               rq(n) = d00
            else
               rq(n) = qk / qc
            endif
!----------------cloud cover ratio ccr----------------------------------
            if(rq(n) < u00k(i,k)) then
                   ccr(n) = d00
            elseif(rq(n) >= us) then
                   ccr(n) = us
            else
                 rqkll  = min(us,rq(n))
                 ccr(n) = h1-sqrt((us-rqkll)/(us-u00k(i,k))) !eq.11
            endif
!
          endif
        enddo
!-------------------ice-water id number iwl------------------------------
!       do n=1,ihpr
!         if (comput(n) .and.  (ww(n) .gt. climit)) then
!           if (tmt0(n) .lt. -15.0
!    *         .or. (tmt0(n) .lt. 0.0 .and. iwl1(n) .eq. 1))
!    *                                      iwl(n) = 1
!             cll(ipr(n),k) = 1.0                           ! cloud cover!
!             cll(ipr(n),k) = min(1.0, ww(n)*cclim(k))      ! cloud cover!
!         endif
!       enddo
!
!---   precipitation production --  auto conversion and accretion
!
      noskip=.true.
      if(noskip) then 
        do n=1,ihpr
          if (comput(n) .and. ccr(n) > 0.0) then
!->rsun 
!rsun       i = ipr(n)
!<-rsun 
            wws    = ww(n)
            cwmk   = max(cons_0, wws)
!           amaxcm = max(cons_0, cwmk - wmink(n))
            if (iwl(n) == 1) then                 !  ice phase
               amaxcm = max(cons_0, cwmk - wmini(ipr(n),k))
               expf      = dt * exp(0.025*tmt0(n))          !eq.26
               psaut     = min(cwmk, psautco*expf*amaxcm)   !eq.25
               ww(n)     = ww(n) - psaut
               cwmk      = max(cons_0, ww(n))
!              cwmk      = max(cons_0, ww(n)-wmini(ipr(n),k))
!rsun          psaci     = min(cwmk, aa2*expf*precsl1(n)*cwmk)
!rsun          psaci     = min(cwmk,
!rsun&              aa2*expf*qs(ipr(n),k)*del(ipr(n),k)
!rsun&              /grav*cwmk)

!
!       use a criterion for the snow collecion 
!

               precsk = max(cons_0,    qqs(n) * condt(n)) 
               psaci     = min(cwmk, aa2*expf*precsk*cwmk)  !eq.28, 29

               ww(n)     = ww(n) - psaci
!->rsun           
               tem = wws - ww(n) 
               if(ww(n) > 0.0 .and. ww(n) > wws) then 
                 ww(n) = wws 
                 tem=0.0 
               endif 
               tem1 = qs(ipr(n),k)
!rsun          precsl(n) = precsl(n) + (wws - ww(n)) * condt(n)
!rsun          qs(ipr(n),k) = qs(ipr(n),k) + max(wws - ww(n),0.)
!rsun          qs(ipr(n),k) = qs(ipr(n),k) + tem 
               qqs(n) = qqs(n) + tem 
               if(wws-ww(n) < 0.0) then 
                 print*,'1precpd neg values:',k,wws,ww(n),  &
                      wws - ww(n),qs(ipr(n),k),tem1
               endif 

!->> rsun    test test 
!               ww(n) = wws 
!               qs(ipr(n),k) = tem  
!<<-rsun 

!<-rsun 
            else                                    !  liquid water
!
!          for using sundqvist precip formulation of rain
!
               amaxcm    = max(cons_0, cwmk - wmink(n))
!!             amaxcm    = cwmk
!rsun          tem1      = precsl1(n) + precrl1(n)
!               tem1      = (qs(ipr(n),k)*rho(ipr(n),k)*dzl(ipr(n),k) + 
!     &                      qr(ipr(n),k)*rho(ipr(n),k)*dzl(ipr(n),k))
!rsun          tem1      = max(cons_0,(qs(ipr(n),k)+qr(ipr(n),k))
!rsun&                    *del(ipr(n),k)/grav)
               tem1      = max(cons_0,(qqs(n)+qqr(n))*condt(n))
               tem2      = min(max(cons_0, 268.0-tt(n)), cons_20)
               tem       = (1.0+c1*sqrt(tem1*rdt)) * (1+c2*sqrt(tem2))
!
               tem2      = amaxcm * cmr * tem / max(ccr(n),cons_p01)
               tem2      = min(cons_50, tem2*tem2)
               praut     = c00  * tem * amaxcm * (1.0-exp(-tem2)) !eq.24
               praut     = min(praut, cwmk)
               ww(n)     = ww(n) - praut
!
!          below is for zhao's precip formulation (water)
!
!              amaxcm    = max(cons_0, cwmk - wmink(n))
!              praut     = min(cwmk, c00*amaxcm*amaxcm)
!              ww(n)     = ww(n) - praut
!
!              cwmk      = max(cons_0, ww(n))
!              tem1      = precsl1(n) + precrl1(n)
!              pracw     = min(cwmk, cr*dt*tem1*cwmk)
!              ww(n)     = ww(n) - pracw
!
!rsun          precrl(n) = precrl(n) + (wws - ww(n)) * condt(n)
               tem = wws - ww(n) 
               if(tem < 0.0) then 
                 tem = 0.0 
                 ww(n) = wws 
               endif 
               tem1 = qr(ipr(n),k)
!rsun          qqr(n)    = qqr(n) + max(wws - ww(n),0.)
!rsun          qr(ipr(n),k)    = qr(ipr(n),k) + max(wws - ww(n),0.)
!rsun          qr(ipr(n),k)    = qr(ipr(n),k) + (wws - ww(n)) ! rsun does this do anything to conservation 
!rsun          qr(ipr(n),k)    = qr(ipr(n),k) + tem ! rsun does this do anything to conservation 
               qqr(n)    = qqr(n) + tem ! rsun does this do anything to conservation 
               if(wws-ww(n) < 0.0) then 
                 print*,'2precpd neg values:',wws,ww(n),          &
                      wws - ww(n),qr(ipr(n),k),tem1,qqr(n) 
               endif 
!->> rsun    test test
!               ww(n) = wws
!               qr(ipr(n),k) = tem
!<<-rsun
!
!hchuang code change [+1l] : add record to record information in vertical
! turn rnp in unit of ww (cwm and q, kg/kg ???)
               rnp(n) = rnp(n) + (wws - ww(n))
            endif
          endif
        enddo
      endif ! if(noskip) 
!
!-----evaporation of precipitation-------------------------
!**** err & ers positive--->evaporation-- negtive--->condensation
!
      noskip = .false.
      noskip = .true.
      if(noskip) then 
        do n=1,ihpr
          if (comput(n)) then
!rsun       i      = ipr(n)
            qk     = max(epsq,  qq(n))
            tmt0k  = max(cons_m30, tmt0(n))
!rsun       precrk = max(cons_0,    precrl(n))
!rsun       precrk = max(cons_0,    qr(i,k)*del(i,k)/grav)
            precrk = max(cons_0,    qqr(n)*condt(n))
!rsun       precsk = max(cons_0,    precsl(n))
!rsun       precsk = max(cons_0,    qs(i,k)*del(i,k)/grav)
            precsk = max(cons_0,    qqs(n)*condt(n))
            amaxrq = max(cons_0,    u00k(i,k)-rq(n)) * conde(n)
!----------------------------------------------------------------------
! increase the evaporation for strong/light prec
!----------------------------------------------------------------------
            ppr    = ke * amaxrq * sqrt(precrk)    !evaporation of rain; eq.35
!           ppr    = ke * amaxrq * sqrt(precrk*rdt)
!->rsun 
!rsun       ppr    = min (ppr,qr(i,k)*condt(n)) 
!<-rsun 
            if (tmt0(n) .ge. 0.) then
              pps = 0.
            else
              pps = (crs1+crs2*tmt0k) * amaxrq * precsk / u00k(i,k) !evaporation of snow; eq.36
!->rsun 
!rsun         pps = min (pps,qs(i,k)*condt(n)) 
!<-rsun 
            end if
!---------------correct if over-evapo./cond. occurs--------------------
            erk=precrk+precsk
            erkdt = erk 
            if(rq(n).ge.1.0e-10)  erk = amaxrq * qk * rdt / rq(n)
            if (ppr+pps .gt. abs(erk)) then
!rsun          rprs   = erk / (precrk+precsk)
               rprs   = erk /erkdt 
               ppr    = precrk * rprs
               pps    = precsk * rprs
            endif
            ppr       = min(ppr, precrk)
            pps       = min(pps, precsk)
            err(n)    = ppr * rconde(n)
            ers(n)    = pps * rconde(n)
!rsun       precrl(n) = precrl(n) - ppr
!rsun       qr(i,k)   = max(0.0,qr(i,k) - ppr /del(i,k)*grav) ! rsun potenially negative values 
!rsun       qr(i,k)   = qr(i,k) - ppr /del(i,k)*grav ! rsun potenially negative values 
!rsun       qr(i,k)   = qr(i,k) - err(n) * dt 
            qqr(n)   = qqr(n) - err(n) * dt   ! rsun: make sure this is positive 
!hchuang code change [+1l] : add record to record information in vertical
! use err for kg/kg/dt not the ppr (mm/dt=kg/m2/dt)
!
            rnp(n) = rnp(n) - err(n)
!
!rsun       precsl(n) = precsl(n) - pps
!rsun       qs(i,k)   = max(0.,qs(i,k) - pps /del(i,k)*grav)
!rsun       qs(i,k)   = qs(i,k) - pps /del(i,k)*grav
!rsun       qs(i,k)   = qs(i,k) - ers(n) * dt 
            qqs(n)   = qqs(n) - ers(n) * dt 
          endif
        enddo
      endif ! if(noskip) 
!--------------------melting of the snow--------------------------------
      noskip=.true.
      if(noskip) then  
        do n=1,ihpr
          if (comput(n)) then
            i      = ipr(n)
!rsun        precrk = max(cons_0,    qr(i,k)*del(i,k)/grav) 
             precrk = max(cons_0,    qqr(n) * condt(n)) 
!rsun        precsk = max(cons_0,    qs(i,k)*del(i,k)/grav) 
             precsk = max(cons_0,    qqs(n)*condt(n))
            if (tmt0(n) .gt. 0.) then
!rsun          amaxps = max(cons_0,    precsl(n))
               amaxps = max(cons_0,    precsk)
               psm1   = csm1 * tmt0(n) * tmt0(n) * amaxps !eq.30; continous melting of melting snow
               psm2   = cws * cr * max(cons_0, ww(n)) * amaxps   !eq.31, 32; melting rate of snow
               ppr    = (psm1 + psm2) * conde(n)
               if (ppr .gt. amaxps) then
                 ppr  = amaxps
                 psm1 = amaxps * rconde(n)
               endif
!rsun          precrl(n) = precrl(n) + ppr
!rsun          qr(i,k) = qr(i,k) + ppr/del(i,k)*grav
!rsun          qqr(n) = qqr(n) + ppr/del(i,k)*grav
               qqr(n) = qqr(n) + psm1*dt 
!
!hchuang code change [+1l] : add record to record information in vertical
! turn ppr (mm/dt=kg/m2/dt) to kg/kg/dt -> ppr/air density (kg/m3)
               rnp(n) = rnp(n) + ppr * rconde(n)
!
!rsun          precsl(n) = precsl(n) - ppr
!rsun          qs(i,k) = qs(i,k) - ppr/del(i,k)*grav
!rsun          qqs(n) = qqs(n) - ppr/del(i,k)*grav
               qqs(n) = qqs(n) - psm1*dt 
            else
               psm1 = d00
            endif
!
!---------------update t and q------------------------------------------
            tt(n) = tt(n) - dtcp * (elwv*err(n)+eliv*ers(n)+eliw*psm1)
            qq(n) = qq(n) + dt * (err(n)+ers(n))
          endif
        enddo
      endif ! if(noskip) 
!
        do n=1,ihpr
          iwl1(n)    = iwl(n)
!rsun          precrl1(n) = max(cons_0, precrl(n))
!rsun          precsl1(n) = max(cons_0, precsl(n))
          i          = ipr(n)
          t(i,k)     = tt(n)
          q(i,k)     = qq(n)
          cwm(i,k)   = ww(n)
          iw(i,k)    = iwl(n)
          qr(i,k)    = qqr(n) 
          qs(i,k)    = qqs(n) 
!hchuang code change [+1l] : add record to record information in vertical
! rnp = precrl1*rconde(n) unit in kg/kg/dt
!
          rainp(i,k) = rnp(n)
        enddo
!
!  move water from vapor to liquid should the liquid amount be negative
!
        do i = 1, im
          if (cwm(i,k) < 0.) then
            tem      = q(i,k) + cwm(i,k)
            if (tem >= 0.0) then
              q(i,k)   = tem
              t(i,k)   = t(i,k) - elwv * rcp * cwm(i,k)
              cwm(i,k) = 0.
            elseif (q(i,k) > 0.0) then
              cwm(i,k) = tem
              t(i,k)   = t(i,k) + elwv * rcp * q(i,k)
              q(i,k)   = 0.0
            endif
          endif
        enddo
!
      enddo                               ! k loop ends here!
!->rsun conservation 
     
      do i = 1, im
        watertot2(i) = 0.0
        do k = 1, km
          watertot2(i) = watertot2(i) +                          &
              (q(i,k)+cwm(i,k)+qr(i,k)+qs(i,k))*del(i,k)/grav
        enddo
        watertot2(i) = watertot2(i) + rn(i)*row 

        if(abs(watertot2(i) - watertot1(i)) > 1.e-6) then 
!          print*,'different2:',kdt, watertot1(i), watertot2(i), rn(i), & !xb110
!     &        abs(watertot2(i) - watertot1(i))                    !     !xb110
        endif 
      enddo
!      print*,'after different2'
!<-rsun 
!->rsun sedimentation 
      do i = 1, im 
        tem = dzl(i,1)/vtr
        nsr(i)  = int(dt/tem+1.)
        tem = dzl(i,1)/vts
        nss(i)  = int(dt/tem+1.)
        onsr(i) = 1./nsr(i)
        onss(i) = 1./nss(i)
      enddo 
      
!      if(me == 0 ) then    !xb110
!        print*,'dzl(1,1)',kdt,dzl(1,1)
!        print*,'nsr,nss,onsr, nnss:',nsr,nss,onsr, onss
!        print*,'precpd before seds : qs:',kdt,qs
!      endif                 !xb110

      do i = 1, im 
        do nt = 1, nss (i)
          do k = km, 1, -1 
            if(qs(i,k) .gt. minp) then 
              seds(i,k)  = vts * qs(i,k) * rho(i,k) 
            else 
              seds(i,k)  = 0.0 
            endif 
          enddo 
!rsun          qs(i,km) = qs(i,km)-seds(i,km)/dzl(i,km)*dt/nss(i)/rho(i,km) 
          qs(i,km) = qs(i,km)-seds(i,km)/                         &
                     del(i,km)*dt/nss(i)*grav
          qs(i,km) = max(qs(i,km),0.0)
          do k = km-1, 1, -1 
            tem1=qs(i,k+1)
            tem =qs(i,k)
!            qs(i,k) = qs(i,k) + 
!     &              (seds(i,k+1) - seds(i,k))/rho(i,k)*
!     &              dt*onss(i)/dzl(i,k)  
            qs(i,k) = qs(i,k) +                                 &
                   (seds(i,k+1) - seds(i,k))/del(i,k)*grav*     &
                    dt*onss(i)
!rsun            qs(i,k) = max(qs(i,k),0.0)
!            if(qs(i,k) < 0.0) then 
!               print*,'neg value qs:', 
!     &                  k,tem1, tem, rho(i,k+1), rho(i,k), 
!     &                  seds(i,k+1), seds(i,k),
!     &                  onss(i),dzl(i,k),qs(i,k) 
! 
!            endif 
          enddo 
          rn(i) = rn(i)+seds(i,1)*dt*onss(i)*rrow
          sn(i) = sn(i)+seds(i,1)*dt*onss(i)*rrow
        enddo 
      enddo 
  
!      if(me==0) then 
!        print*,'precpd after seds: seds:',kdt,seds
!        print*,'precpd after seds : rn:',kdt,rn
!        print*,'precpd after seds : sn:',kdt,sn
!        print*,'precpd after seds : qs:',kdt,qs
!        print*, 'precpd before sedr : qr:',kdt, qr 
!      endif 

      do i = 1, im 
        do nt = 1, nsr(i)
          do k = km, 1, -1 
            if(qr(i,k) .gt. minp) then 
              sedr(i,k)  = vtr * qr(i,k) * rho(i,k) 
            else 
              sedr(i,k)  = 0.0 
            endif 
          enddo 
          qr(i,km) = qr(i,km)-sedr(i,km)/                         &
                     del(i,km)*dt/nsr(i)*grav
          do k = km-1, 1, -1 
            tem = qr(i,k) 
            tem1 = qr(i,k+1)
            
            qr(i,k) = qr(i,k) +                                   &
                   (sedr(i,k+1) - sedr(i,k))/del(i,k)*grav*       &
                    dt*onsr(i)
!rsun            qr(i,k) = max(qr(i,k),0.0)
!            if(qr(i,k) < 0.0) then 
!               print*,'neg value qr:', 
!     &                  k,tem1, tem, rho(i,k+1), rho(i,k), 
!     &                  sedr(i,k+1), sedr(i,k),
!     &                  onsr(i),dzl(i,k),qr(i,k) 
! 
!            endif 
          enddo 
          rn(i) = rn(i)+sedr(i,1)*dt*onsr(i)*rrow

        enddo 
      enddo 

!      if(me ==0) then 
!        print*,'precpd after sedr and dt: sedr:',kdt,sedr
!        print*,'precpd after sedr  dt: rn:',kdt,rn
!        print*,'precpd after sedr dt : sn:',kdt,sn
!        print*,'precpd after sedr dt : qr:',kdt,qr
!      endif 

      do i = 1, im 
        rid = rn(i)
        if (rid < 1.e-13) then
           sr(i) = 0.
        else
           sr(i) = sn(i)/rid
        endif
      enddo 
!->rsun conservation 

      do i = 1, im
        watertot2(i) = 0.0
        do k = 1, km
          watertot2(i) = watertot2(i) +                                 &
              (q(i,k)+cwm(i,k)+qr(i,k)+qs(i,k))*del(i,k)/grav
        enddo
        watertot2(i) = watertot2(i) + rn(i)*row 

        if(abs(watertot2(i) - watertot1(i)) > 1.e-6) then 
!          print*,'different1:',kdt, watertot1(i), watertot2(i), rn(i),  &!xb110
!     &        abs(watertot2(i) - watertot1(i))                           !xb110
        endif 
      enddo
!<-rsun end of sedimentation 

!**********************************************************************
!-----------------------end of precipitation processes-----------------
!**********************************************************************
!
!->rsun 
!      do n=1,ihpr
!        i = ipr(n)
!        rn(i) = (precrl1(n)  + precsl1(n)) * rrow  ! prate: flux to rate in DT (not /s) 
!!
!!----sr=1 if sfc prec is rain ; ----sr=-1 if sfc prec is snow
!!----sr=0 for both of them or no sfc prec
!!
!!        rid = 0.
!!        sid = 0.
!!        if (precrl1(n) .ge. 1.e-13) rid = 1.
!!        if (precsl1(n) .ge. 1.e-13) sid = -1.
!!        sr(i) = rid + sid  ! sr=1 --> rain, sr=-1 -->snow, sr=0 -->both
!! chuang, june 2013: change sr to define fraction of frozen precipitation instead
!! because wpc uses it in their winter experiment
!
!        rid = precrl1(n) + precsl1(n)
!        if (rid < 1.e-13) then
!           sr(i) = 0.
!        else
!           sr(i) = precsl1(n)/rid
!        endif 
!      enddo
!<-rsun 
!
      return
      end
