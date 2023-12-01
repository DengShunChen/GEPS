      subroutine gscond (im,ix,km,dt,prsl,ps,q,cwm,t          &
      ,                  tp, qp, psp,tp1,qp1,psp1, u, lprnt, fwd)
!
!     ******************************************************************
!     *                                                                *
!     *  subroutine for grid-scale condensation & evaporation          *
!     *  for the mrf model at ncep.                                    *
!     *                                                                *
!     ******************************************************************
!     *                                                                *
!     *  created by:   q.  zhao         jan. 1995                      *
!     *  modified by:  h.-l. pan        sep. 1998                      *
!     *  modified by:  s. moorthi       aug. 1999, 2000                *
!     *  modified by:  h. juang         jun. 2004                      *
!     *                                                                *
!     *  references:                                                   *
!     *                                                                *
!     ******************************************************************
!
!     ************ important note *************
!     in this code vertical indexing runs from surface to top of the
!     model
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
      implicit none

      real  psat, hvap, grav, hfus, ttp, rd, rv, cp, eps, epsm1
      parameter(  psat = 6.1078e+2 )
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
      real                   g,    h1,   h2, h1000                      &
      ,                     d00,  d125,  d5,   elwv, eliv               &
      ,                     epsq, tm10,  eliw, arcp                     &
      ,                     a1,   r,     cpr,  rcpr, rcp
      parameter (h1=1.e0,       h2=2.e0,               h1000=1000.0     &
      ,          d00=0.e0,      d125=.125e0,           d5=0.5e0         &
      ,          a1=psat                                                &
      ,          elwv=hvap,     eliv=hvap+hfus, g=grav                  &
      ,          epsq=2.e-12,   tm10=ttp-10.,   r=rd                    &
      ,          cpr=cp*r,      rcpr=h1/cpr,    rcp=h1/cp)
!
      integer im, ix, km
      real                  q(ix,km),    t(ix,km),    cwm(ix,km)        &
      ,                     prsl(ix,km), ps(im), dt                     &
      ,                     tp(ix,km),   qp(ix,km),   psp(im)           &
      ,                     tp1(ix,km),  qp1(ix,km),   psp1(im)         &
      ,                     u(ix,km)                                    &
      ,                     t0(ix,km),   q0(ix,km),   ps0(im)
!
      real                   qi(im), qint(im), ccrik, e0                &
      ,                      cond,   rdt, us, cclimit, climit           &
      ,                      u00b,   u00t, tmt0, tmt15, qik, cwmik      &
      ,                      ai, bi, qw, u00ik, tik, pres, pp0, fi      &
      ,                      at, aq, ap, fiw, elv, qc, rqik             &
      ,                      rqikk, tx1, tx2, tx3, es, qs               &
      ,                      tsq, delq, condi, cone0, us00, ccrik1      &
      ,                      aa, ab, ac, ad, ae, af, ag                 &
      ,                      el2orc, albycp, vprs(im), rdt2
      integer iw(im,km), i, k, iwik
      logical lprnt,fwd
!
!--- merge qsatq in here
      real vpsat(191),pqs,qqq,tem,t1
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
      tem = 1./1.622
!

!
!-----------------prepare constants for later uses-----------------
!

      el2orc = hvap*hvap / (rv*cp)
      albycp = hvap / cp
!
! for 2 delta t
      rdt     = h1/dt
!
! for tendency caculation ,2-time level was used, just need 1 dt
      rdt2    = h1/dt*2.
!
      us      = h1
      cclimit = 1.0e-3
      climit  = 1.0e-20
!
      do  i = 1, im
        iw(i,km) = d00
      enddo
!
!  check for first time step
!
      if (tp(1,1) .lt. 1.) then
        do k = 1, km
          do i = 1, im
            tp(i,k) = t(i,k)
            qp(i,k) = max(q(i,k),epsq)
            tp1(i,k) = t(i,k)
            qp1(i,k) = max(q(i,k),epsq)
          enddo
        enddo
        do i = 1, im
          psp(i)  = ps(i)
          psp1(i)  = ps(i)
        enddo
      endif

!-----------------------------------------
!
!*************************************************************
!*******begining of  grid-scale condensation/evap. loop*******
!*************************************************************
!
!     do k = km-1,2,-1
      do k = km,1,-1
!       vprs(:) = 0.001 * fpvs(t(:,k))       ! fpvs in pa
!-----------------------------------------------------------------------
!------------------qw, qi and qint--------------------------------------
        do i = 1, im
          tmt0  = t(i,k)-273.16
          tmt15 = min(tmt0,-15.)
          qik   = max(q(i,k),epsq)
          cwmik = max(cwm(i,k),climit)
!
!  the global qsat computation is done in cb
          pres    = prsl(i,k)
!fpvs
           qw      = min(pres, 0.001 * fpvs(t(i,k)))
           qw      = eps * qw / (pres + epsm1 * qw)
!qsatq
!          call qsatq(1,t(i,k),pres*10.,qw)  !prsl in unit of cb
!         t1 = max(1.00001, min(190.999, t(i,k)-182.16))
!         ic = int(t1)
!         pqs = pres*10.
!         qqq = min(tem*pqs, vpsat(ic)+(vpsat(1+ic)-vpsat(ic))  &
!                                    *(t1-float(ic)))
!!         qw = 0.622*qqq/(pqs-qqq)
!          qw = 0.622*qqq/(pqs+epsm1*qqq)
!
          qw      = max(qw,epsq)
          qi(i)   = qw
          qint(i) = qw
!-------------------ice-water id number iw------------------------------
!
          iw(i,k) = 0
!
         if(tmt0.lt.-15.0) then
           u00ik = u(i,k)
           fi    = qik - u00ik*qi(i)
           if(fi.gt.d00.or.cwmik.gt.climit) then
              iw(i,k) = 1
           else
             iw(i,k) = 0
           end if
         end if

         if(tmt0.ge.0.0) then
           iw(i,k) = 0
         end if

         if (tmt0 .lt. 0.0 .and. tmt0 .ge. -15.0) then
           iw(i,k) = 0
           if (k .lt. km) then
           if (iw(i,k+1) .eq. 1 .and. cwmik .gt. climit) iw(i,k) = 1
           endif
         end if
!
        enddo

!--------------condensation and evaporation of cloud--------------------
        do i = 1, im
!------------------------at, aq and dp/dt-------------------------------
          qik   = max(q(i,k),epsq)
          cwmik = max(cwm(i,k),climit)
          iwik  = iw(i,k)
          u00ik = u(i,k)
          tik   = t(i,k)
          pres  = prsl(i,k)   * h1000
          if ( fwd ) then
            pp0   = (pres / ps(i)) * psp1(i)
            at    = (tik-tp1(i,k)) * rdt
            aq    = (qik-qp1(i,k)) * rdt
          else
           pp0   = (pres / ps(i)) * psp(i)
           at    = (tik-tp(i,k)) * rdt
           aq    = (qik-qp(i,k)) * rdt
          endif
          ap    = (pres-pp0)    * rdt
! tik and tp just across 2-time level
!         at    = (tik-tp(i,k)) * rdt2
!         aq    = (qik-qp(i,k)) * rdt2
!         ap    = (pres-pp0)    * rdt2
!----------------the satuation specific humidity------------------------
          fiw   = float(iwik)
          elv   = (h1-fiw)*elwv    + fiw*eliv
          qc    = (h1-fiw)*qint(i) + fiw*qi(i)
!     if (lprnt) print *,' qc=',qc,' qint=',qint(i),' qi=',qi(i)
!----------------the relative humidity----------------------------------
          if(qc.le.1.0e-10) then
            rqik=d00
          else
            rqik = qik/qc
          endif
!----------------cloud cover ratio ccrik--------------------------------
          if (rqik .lt. u00ik) then
             ccrik = d00
          elseif(rqik.ge.us) then
             ccrik = us
          else
             rqikk  = min(us,rqik)
             ccrik = h1-sqrt((us-rqikk)/(us-u00ik))      ! (b) eq(11)
          endif
!-----------correct ccr if it is too small in large cwm regions--------
!         if(ccrik.ge.0.01.and.ccrik.le.0.2.and
!    &          .cwmik.ge.0.2e-3) then
!          ccrik=min(1.0,cwmik*1.0e3)
!         end if
!----------------------------------------------------------------------
!   if no cloud exists then evaporate any existing cloud condensate
!----------------evaporation of cloud water-----------------------------
          e0 = d00
!
          if (ccrik.le.cclimit.and.cwmik.gt.climit)  then
!  
!!            if(iwik.eq.1 .and. t(i,k).le.248.16)then                         !begin test
!!              tmt0  = t(i,k)-248.16
!!              u00ik=u00ik*exp(0.0015*tmt0)
!!            endif                                       !end test
!
!   first iteration - increment halved
!
            tx1 = tik
            tx3 = qik
!fpvs
!pres in unit pa
             es   = min(pres, fpvs(tx1))
             qs   = u00ik * eps * es / (pres + epsm1*es)
!qsatq
!!           call qsatq(1,tx1,pres/100.,qs)     !pres in unit of pa
!           t1 = max(1.00001, min(190.999, tx1-182.16))
!           ic = int(t1)
!           pqs = pres*0.01
!           qqq = min(tem*pqs, vpsat(ic)+(vpsat(1+ic)-vpsat(ic))  &
!                                    *(t1-float(ic)))
!!            qs = 0.622*qqq/(pqs-qqq)
!!            qs = 0.622*qqq/(pqs+epsm1*qqq)
!            qs   = u00ik*qs
!
            tsq  = tx1 * tx1
            delq = 0.5 * (qs - tx3) * tsq / (tsq + el2orc * qs)
!
            tx2   = delq
            tx1   = tx1 - delq * albycp
            tx3   = tx3 + delq
!
!   second iteration
!
             es   = min(pres, fpvs(tx1))
             qs   = u00ik * eps * es / (pres + epsm1*es)
!
!!           call qsatq(1,tx1,pres/100.,qs)     !pres in unit of pa
!           t1 = max(1.00001, min(190.999, tx1-182.16))
!           ic = int(t1)
!           pqs = pres*0.01
!           qqq = min(tem*pqs, vpsat(ic)+(vpsat(1+ic)-vpsat(ic)) &
!                                    *(t1-float(ic)))
!           qs = 0.622*qqq/(pqs-qqq)
!           qs = 0.622*qqq/(pqs+epsm1*qqq)
!            qs   = u00ik*qs
!
            tsq  = tx1 * tx1
            delq = (qs - tx3) * tsq / (tsq + el2orc * qs)
!
            tx2  = tx2 + delq
            tx1  = tx1 - delq * albycp
            tx3  = tx3 + delq
!
!   third iteration
!fpvs
             es   = min(pres, fpvs(tx1))
             qs   = u00ik * eps * es / (pres + epsm1*es)
!qsatq
!!            call qsatq(1,tx1,pres/100.,qs)     !pres in unit of pa
!           t1 = max(1.00001, min(190.999, tx1-182.16))
!           ic = int(t1)
!           pqs = pres*0.01
!           qqq = min(tem*pqs, vpsat(ic)+(vpsat(1+ic)-vpsat(ic))  &
!                                    *(t1-float(ic)))
!!           qs = 0.622*qqq/(pqs-qqq)
!!           qs = 0.622*qqq/(pqs+epsm1*qqq)
!            qs   = u00ik*qs
!
            tsq  = tx1 * tx1
            delq = (qs - tx3) * tsq / (tsq + el2orc * qs)
            tx2  = tx2 + delq
!
            e0   = max(tx2*rdt, 0.0)

!           e0 = max(qc*(u00ik-rqik)*rdt, 0.0)
            e0 = min(cwmik*rdt,   e0)
            e0 = max(0.0,e0)
          end if
!
!   if cloud cover > 0.2 condense water vapor in to cloud condensate
!-----------the eqs. for cond. has been reorganized to reduce cpu------
!
          cond = d00
!         if (ccrik .gt. 0.20 .and. qc .gt. epsq) then
          if (ccrik .gt. cclimit .and. qc .gt. epsq) then
             us00   = us  - u00ik
             ccrik1 = 1.0 - ccrik
             aa     = eps*elv*pres*qik
             ab     = ccrik*ccrik1*qc*us00
             ac     = ab + 0.5*cwmik
             ad     = ab * ccrik1
             ae     = cpr*tik*tik
             af     = ae * pres
             ag     = aa * elv
             ai     = cp * aa
             cond   = (ac-ad)*(af*aq-ai*at+ae*qik*ap)/(ac*(af+ag))
!-----------check & correct if over condensation occurs-----------------
             condi  = (qik   -u00ik   *qc*1.0)*rdt
             cond   = min(cond, condi)
!----------check & correct if supersatuation is too high----------------
!             qtemp=qik-max(0.,(cond-e0))*dt
!             if(qc.le.1.0e-10) then
!               rqtmp=0.0
!             else
!               rqtmp=qtemp/qc
!             end if
!             if(rqtmp.ge.1.10) then
!               cond=(qik-1.10*qc)*rdt
!             end if
!-----------------------------------------------------------------------
             cond = max(cond, d00)
!-------------------update of t, q and cwm------------------------------
!
          end if
          cone0    = (cond-e0) * dt
          cwm(i,k) = cwm(i,k) + cone0
          t(i,k)   = t(i,k)   + elv*rcp*cone0
          q(i,k)   = q(i,k)   - cone0
        enddo                                  ! end of i-loop!

      enddo                                    ! end of k-loop!

!
!*********************************************************************
!****************end of the condensation/evaporation loop*************
!*********************************************************************
!----------------store t, q, ps for next time step
      do k = 1, km
        do i = 1, im
          tp(i,k) = tp1(i,k)
          qp(i,k) = qp1(i,k)

          tp1(i,k) = t(i,k)
          qp1(i,k) = max(q(i,k),epsq)
        enddo
      enddo
      do i = 1, im
        psp(i) = psp1(i)
        psp1(i) = ps(i)
      enddo

!-----------------------------------------------------------------------
      return
      end
