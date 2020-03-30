       subroutine precpd (im,ix,km,dt,del,prsl,ps,q,cwm,t,rn,u00k,  &
                          psautco,lprnt)
!
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
!     *  modified by h juang                           jun, 2004       *
!     *                                                                *
!     *  references:                                                   *
!     *                                                                *
!     *  zhao and carr (1997), monthly weather review (august)         *
!     *  sundqvist et al., (1989) monthly weather review. (august)     *
!     *                                                                *
!     ******************************************************************
!
!     ************ important note *************
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
!       t(ix,km)   : temperature       (updated in the code)
!       rn(im)     : precipitation over one time-step dt (m/dt)
!       sr(im)     : index (=-1 snow, =0 rain/snow, =1 rain)
!       tcw(im)    : vertically integrated liquid water (kg/m**2)
!       cll(ix,km) : cloud cover
!

      implicit none

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

      real fpvs
!
      real                  g,      h1,    h2,   h1000                  &
      ,                     h1000g, d00,   d125, d5                     &
      ,                     elwv,   eliv,  row                          &
      ,                     epsq,   dldt,  tm10, eliw                   &
      ,                     rcp,    rrow
       parameter (g=grav,         h1=1.e0,     h2=2.e0,     h1000=1000.0 &
      ,           h1000g=h1000/g, d00=0.e0,    d125=.125e0, d5=0.5e0    &
      ,           elwv=hvap,      eliv=hvap+hfus,   row=1.e3            &
      ,           epsq=2.e-12,    dldt=2274.e0,tm10=ttp-10.0            &
      ,           eliw=eliv-elwv, rcp=h1/cp,   rrow=h1/row)
!
      integer im, ix, km, lat
      real                  q(ix,km),   t(ix,km),    cwm(ix,km)         &
      ,                                 del(ix,km),  prsl(ix,km)        &
!     ,                     cll(im,km), del(ix,km),  prsl(ix,km)
      ,                     ps(im),     rn(im),      sr(im)             &
      ,                     tcw(im),    dt                              &
      ,                     u00k(ix,km),psautco(im)
!
!
      real                  err(im),      ers(im),     precrl(im)       &
      ,                     precsl(im),   precrl1(im), precsl1(im)      &
      ,                     rq(im),       condt(im)                     &
      ,                     conde(im),    rconde(im),  tmt0(im)         &
      ,                     wmin(im,km),  wmink(im),   pres(im)         &
      ,                     wmini(im,km), ccr(im),     cclim(km)        &
      ,                     tt(im),       qq(im),      ww(im)           &
      ,                     wfix(km),                  es(im)           &
      ,                     zaodt
!
      integer iw(im,km), ipr(im), iwl(im),     iwl1(im)
!
       logical comput(im)
       logical lprnt
!
      real                  ke,   rdt,  us, cclimit, climit, cws, csm1  &
      ,                     crs1, crs2, cr, aa2, aa1, dtcp,   c00, cmr   &
      ,                     tem,  c1,   c2, wwn                         &
!     ,                     tem,  c1,   c2, u00b,    u00t,   wwn
      ,                     precrk, precsk, pres1,   qk,     qw,  qi    &
      ,                     ai,     bi, qint, fiw, wws, cwmk, expf      &
      ,                     psaut, psaci, amaxcm, tem1, tem2            &
      ,                     tmt0k, tmt15, psm1, psm2, ppr               &
      ,                     rprs,  erk,   pps, sid, rid, amaxps         &
      ,                     praut, pracw, fi, qc, amaxrq, rqkll
      integer i, k, ihpr, n
!
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
      temx = 1.0/1.622

!CWB2015 fixed undefine problem in line 389
      ccr=0.0
!-----------------------preliminaries ---------------------------------
!
!     do k=1,km
!       do i=1,im
!         cll(i,k) = 0.0
!       enddo
!     enddo
!
      rdt     = h1 / dt
      ke      = 2.0e-5  ! commented on 09/10/99
!     ke      = 2.0e-6
!     ke      = 1.0e-5
!     ke      = 5.0e-5
      us      = h1
      cclimit = 1.0e-3
      climit  = 1.0e-20
      cws     = 0.025
!
!       zaodt   = 80.0 * rdt
      zaodt   = 800.0 * rdt    !original
!
!      csm1    = 5.0000e-8   * zaodt
!      crs1    = 5.00000e-6  * zaodt
!      crs2    = 6.66600e-10 * zaodt
!      cr      = 5.0e-4      * zaodt
!      aa2     = 1.25e-3     * zaodt
!
      csm1    = 5.0000e-8   * zaodt        !test(for psm1)
      cr      = 5.0e-4      * zaodt        !test(for psm2)
      crs1    = 5.00000e-6  * zaodt        !test(for Ers)
      crs2    = 6.66600e-10 * zaodt        !test(for Ers)
      ke      = ke * sqrt(rdt)             !test(for Err)
      aa2     = 1.25e-3     * zaodt        !test(for Psaci)
!      aa1     = 4.0e-4                     !test(for Psaut) original
!      aa1     = 4.0e-4                     !test(for Psaut) original
!byl      aa1     = 8.0e-4                     !TCo639L72
!      aa1     = 8.0e-4                     !test(for Psaut) test (rrt10)
!      aa1     = 6.0e-4                     !test(for Psaut) test (rrt11)
!      aa1     = 5.0e-4                     !test(for Psaut) test (rrt13)
!
!      ke      = ke * sqrt(rdt)
!     ke      = ke * sqrt(zaodt)
!
      dtcp    = dt * rcp
!
!     c00 = 1.5e-1 * dt
!     c00 = 10.0e-1 * dt
!     c00 = 3.0e-1 * dt          !05/09/2000
!ops      c00 = 1.0e-4 * dt          !05/09/2000
!test rsas8p2
!      c00 = 3.0e-4 * dt          !2008/12/5
!byl      c00 = 1.0e-4 * dt          ! test(for Praut)
      c00 = 1.0e-4 * dt          ! TCo639L72
!
      cmr = 1.0 / 3.0e-4
!     cmr = 1.0 / 5.0e-4
!     c1  = 100.0
      c1  = 300.0
      c2  = 0.5
!
!
!--------calculate c0 and cmr using lc at previous step-----------------
!
      do k=1,km
        do i=1,im
          tem   = (prsl(i,k)*0.01)
!         tem   = sqrt(tem)
          iw(i,k)    = 0.0
          wmin(i,k)  = 1.0e-5 * tem
          wmini(i,k) = 1.0e-5 * tem       ! testing for ras
!         wmini(i,k) = 1.0e-6 * tem       ! for sas
!          wmin(i,k)  = 0.5e-5 * tem
!          wmini(i,k) = 0.5e-5 * tem       ! testing for ras
!          wmin(i,k)  = 1.5e-5 * tem      ! for tem   = sqrt(tem)
!          wmini(i,k) = 1.5e-5 * tem      ! for tem   = sqrt(tem)
        enddo
      enddo
      do i=1,im
!       c0(i)  = 1.5e-1
!       cmr(i) = 3.0e-4
!
        iwl1(i)    = 0
        precrl1(i) = d00
        precsl1(i) = d00
        comput(i)  = .false.
        rn(i)      = d00
        sr(i)      = d00
      enddo
!------------select columns where rain can be produced--------------
      do k=1, km-1
        do i=1,im
          tem = min(wmin(i,k), wmini(i,k))
          if (cwm(i,k) .gt. tem) comput(i) = .true.
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
          precrl(n) = precrl1(n)
          precsl(n) = precsl1(n)
          err  (n)  = d00
          ers  (n)  = d00
          iwl  (n)  = 0
!
          i         = ipr(n)
          tt(n)     = t(i,k)
          qq(n)     = q(i,k)
          ww(n)     = cwm(i,k)
          wmink(n)  = wmin(i,k)
          pres(n)   = h1000 * prsl(i,k)
!
          precrk = max(0.,    precrl1(n))
          precsk = max(0.,    precsl1(n))
          wwn    = max(ww(n), climit)
!         if (wwn .gt. wmink(n) .or. (precrk+precsk) .gt. d00) then
          if (wwn .gt. climit .or. (precrk+precsk) .gt. d00) then
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
            conde(n)  = (h1000*dt/g) * del(i,k)
            condt(n)  = conde(n) * rdt
            rconde(n) = h1 / conde(n)
            qk        = max(epsq,  qq(n))
            tmt0(n)   = tt(n) - 273.16
            wwn       = max(ww(n), climit)
!
!  the global qsat computation is done in pa
            pres1   = pres(n)
!           qw      = es(n)
            qw      = min(pres1, fpvs(tt(n)))
            qw      = eps * qw / (pres1 + epsm1 * qw)
!
!!           call qsatq(1,tt(n),pres(n)/100.,qw)
!          t1 = max(1.00001, min(190.999, tt(n)-182.16))
!         ic = int(t1)
!         pqs = pres(n)/100.
!         qqq = min(temx*pqs, vpsat(ic)+(vpsat(1+ic)-vpsat(ic))  &
!                                    *(t1-float(ic)))
!!         qw = 0.622*qqq/(pqs-qqq)
!         qw = 0.622*qqq/(pqs+epsm1*qqq)
!
          qw = max(qw,epsq)
!
!
!           tmt15 = min(tmt0(n), -15.)
!           ai    = 0.008855
!           bi    = 1.0
!           if (tmt0(n) .lt. -20.0) then
!             ai = 0.007225
!             bi = 0.9674
!           endif
!           qi   = qw * (bi + ai*min(tmt0(n),0.))
!           qint = qw * (1.-0.00032*tmt15*(tmt15+15.))
!
            qi   = qw
            qint = qw
!           if (tmt0(n).le.-40.) qint = qi
!
!-------------------ice-water id number iw------------------------------
            if(tmt0(n).lt.-15.) then
               fi = qk - u00k(i,k)*qi
               if(fi.gt.d00.or.wwn.gt.climit) then
                  iwl(n) = 1
               else
                  iwl(n) = 0
               endif
!           endif
            elseif (tmt0(n).ge.0.) then
               iwl(n) = 0
!
!           if(tmt0(n).lt.0.0.and.tmt0(n).ge.-15.0) then
            else
              iwl(n) = 0
              if(iwl1(n).eq.1.and.wwn.gt.climit) iwl(n)=1
            endif
!
!           if(tmt0(n).ge.0.) then
!              iwl(n) = 0
!           endif
!----------------the satuation specific humidity------------------------
            fiw   = float(iwl(n))
            qc    = (h1-fiw)*qint + fiw*qi
!----------------the relative humidity----------------------------------
            if(qc .le. 1.0e-10) then
               rq(n) = d00
            else
               rq(n) = qk / qc
            endif
!----------------cloud cover ratio ccr----------------------------------
            if(rq(n).lt.u00k(i,k)) then
                   ccr(n)=d00
            elseif(rq(n).ge.us) then
                   ccr(n)=us
            else
                 rqkll=min(us,rq(n))
                 ccr(n)= h1-sqrt((us-rqkll)/(us-u00k(i,k)))      !  eq(11) (b)
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
        do n=1,ihpr
          if (comput(n) .and. ccr(n) .gt. 0.0) then
            wws    = ww(n)
            cwmk   = max(0.0, wws)
!           amaxcm = max(0.0, cwmk - wmink(n))
            if (iwl(n) .eq. 1) then                 !  ice phase
               amaxcm = max(0.0, cwmk - wmini(ipr(n),k))
               expf      = dt * exp(0.025*tmt0(n))               !eq(26)
!              psaut     = min(cwmk, 2.0e-3*expf*amaxcm)
!              psaut     = min(cwmk, 1.0e-3*expf*amaxcm)
!              psaut     = min(cwmk, 5.0e-4*expf*amaxcm)
!cjh               psaut     = min(cwmk, 4.0e-4*expf*amaxcm)
!byl               psaut     = min(cwmk, aa1*expf*amaxcm)            !eq(25) aa1=4.0e-4
               psaut     = min(cwmk, psautco(ipr(n))*expf*amaxcm)            !eq(25) aa1=4.0e-4

               ww(n)     = ww(n) - psaut
               cwmk      = max(0.0, ww(n))
!              cwmk      = max(0.0, ww(n)-wmini(ipr(n),k))
               psaci     = min(cwmk, aa2*expf*precsl1(n)*cwmk)   !eq(28)

               ww(n)     = ww(n) - psaci

               precsl(n) = precsl(n) + (wws - ww(n)) * condt(n)
            else                                                 !  liquid water
!
!          for using sundqvist precip formulation of rain
!
!              amaxcm    = max(0.0, cwmk - wmink(n))
               amaxcm    = cwmk
               tem1      = precsl1(n) + precrl1(n)
               tem2      = min(max(0.0, 268.0-tt(n)), 20.0)
               tem       = (1.0+c1*sqrt(tem1*rdt)) * (1+c2*sqrt(tem2))
!
               tem2      = amaxcm * cmr * tem / max(ccr(n),0.01)
               tem2      = min(50.0, tem2*tem2)
               praut     = c00  * tem * amaxcm * (1.0-exp(-tem2))      !eq(24)
               praut     = min(praut, cwmk)
               ww(n)     = ww(n) - praut
!
!          below is for zhao's precip formulation (water)
!
!              amaxcm    = max(0.0, cwmk - wmink(n))
!              praut     = min(cwmk, c00*amaxcm*amaxcm)
!              ww(n)     = ww(n) - praut
!
!              cwmk      = max(0.0, ww(n))
!              tem1      = precsl1(n) + precrl1(n)
!              pracw     = min(cwmk, cr*dt*tem1*cwmk)                  !eq(27)
!              ww(n)     = ww(n) - pracw
!
               precrl(n) = precrl(n) + (wws - ww(n)) * condt(n)
            endif
          endif
        enddo
!
!-----evaporation of precipitation-------------------------
!**** err & ers positive--->evaporation-- negtive--->condensation
!
        do n=1,ihpr
          if (comput(n)) then
            i      = ipr(n)
            qk     = max(epsq,  qq(n))
            tmt0k  = max(-30.0, tmt0(n))
            precrk = max(0.,    precrl(n))
            precsk = max(0.,    precsl(n))
            amaxrq = max(0.,    u00k(i,k)-rq(n)) * conde(n)
!----------------------------------------------------------------------
! increase the evaporation for strong/light prec
!----------------------------------------------------------------------
            ppr    = ke * amaxrq * sqrt(precrk)                        ! eq(35)(Err)
!           ppr    = ke * amaxrq * sqrt(precrk*rdt)
            if (tmt0(n) .ge. 0.) then
              pps = 0.
            else
              pps = (crs1+crs2*tmt0k) * amaxrq * precsk / u00k(i,k)    !eq(36)(Ers)
            end if
!---------------correct if over-evapo./cond. occurs--------------------
            erk=precrk+precsk
            if(rq(n).ge.1.0e-10)  erk = amaxrq * qk * rdt / rq(n)
            if (ppr+pps .gt. abs(erk)) then
               rprs   = erk / (precrk+precsk)
               ppr    = precrk * rprs
               pps    = precsk * rprs
            endif
            ppr       = min(ppr, precrk)
            pps       = min(pps, precsk)
            err(n)    = ppr * rconde(n)
            ers(n)    = pps * rconde(n)
            precrl(n) = precrl(n) - ppr
            precsl(n) = precsl(n) - pps
          endif
        enddo
!--------------------melting of the snow--------------------------------
        do n=1,ihpr
          if (comput(n)) then
            if (tmt0(n) .gt. 0.) then
               amaxps = max(0.,    precsl(n))
               psm1   = csm1 * tmt0(n) * tmt0(n) * amaxps              !eq(30)
               psm2   = cws * cr * max(0.0, ww(n)) * amaxps            !eq(31,32)
               ppr    = (psm1 + psm2) * conde(n)
               if (ppr .gt. amaxps) then
                 ppr  = amaxps
                 psm1 = amaxps * rconde(n)
               endif
               precrl(n) = precrl(n) + ppr
               precsl(n) = precsl(n) - ppr
            else
               psm1 = d00
            endif
!
!---------------update t and q------------------------------------------
!
            tt(n) = tt(n) - dtcp * (elwv*err(n)+eliv*ers(n)+eliw*psm1)
            qq(n) = qq(n) + dt * (err(n)+ers(n))
          endif
        enddo
!
        do n=1,ihpr
          iwl1(n)    = iwl(n)
          precrl1(n) = max(0.0, precrl(n))
          precsl1(n) = max(0.0, precsl(n))
          i          = ipr(n)
          t(i,k)     = tt(n)
          q(i,k)     = qq(n)
          cwm(i,k)   = ww(n)
          iw(i,k)    = iwl(n)
        enddo
!
!  move water from vapor to liquid should the liquid amount be negative
!
        do i = 1, im
          if (cwm(i,k) .lt. 0.) then
            q(i,k)   = q(i,k) + cwm(i,k)
            t(i,k)   = t(i,k) - elwv * rcp * cwm(i,k)
            cwm(i,k) = 0.
          endif
        enddo
!
      enddo                               ! k loop ends here!
!**********************************************************************
!-----------------------end of precipitation processes-----------------
!**********************************************************************
!
      do n=1,ihpr
        i = ipr(n)
        rn(i) = (precrl1(n)  + precsl1(n)) * rrow  ! precip at surface
!
!----sr=1 if sfc prec is rain ; ----sr=-1 if sfc prec is snow
!----sr=0 for both of them or no sfc prec
!
        rid = 0.
        sid = 0.
        if (precrl1(n) .ge. 1.e-13) rid = 1.
        if (precsl1(n) .ge. 1.e-13) sid = -1.
        sr(i) = rid + sid  ! sr=1 --> rain, sr=-1 -->snow, sr=0 -->both
      enddo
!
      return
      end
