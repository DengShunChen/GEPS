      subroutine mixpbl( nxj,mn,kk,ktpbl                               &
                       , dt,g,hgt,u,v,t,q,ut,vt,tt,qt,e,eps            &
                       , xkm,xkh,zl,sfcw,ustar,tstar,qstar,itypbl      &
                       ,dhgt,ro2,dhgtz)
!
!##################################################################
!                     subroutine description
!
! 1. function description
!
!       this subroutine calculates eddy transport coefficient by
!       computing e and eps, and update ut, vt, tt, qt due to
!       turbulence mixing effects using backward time scheme.
!
! 2. common block specification
!
! 3. parameter specification
!
!     mn : dimension of horizontal-direction
!     kk : dimension of z-direction
!  ktpbl : starting level to onsider vertical mixing effects
!     dt : time step for ut, vt, tt, qt updates                (s)
!      g : gravity                                              (m/s2)
!    hgt : height above ground    (mn,kk)  at  current time lvl (m)
!      u : x-direction velocity   (mn,kk)  at  current time lvl (m/s)
!      v : y-direction velocity   (mn,kk)  at  current time lvl (m/s)
!      t : potent temp on odd lvl  (mn,kk) at  current time lvl (k)
!      q : moisture    on odd lvl  (mn,kk) at  current time lvl (kg/kg)
!     ut : x-direction velocity   (mn,kk)  at past/future t lvl (m/s)
!     vt : y-direction velocity   (mn,kk)  at past/future t lvl (m/s)
!     tt : potent temp on odd lvl  (mn,kk) at past/future t lvl (k)
!     qt : moisture    on odd lvl  (mn,kk) at past/future t lvl (kg/kg)
!     e  : turbulence kinetic energy (mn,kk) at current/future  (m2/s2)
!    eps : tke dissipation           (mn,kk) at current/future  (m2/s2)
!    xkm : eddy mixing coef for u, v (mn,kk) at current t lvl   (m2/s2)
!    xkh : eddy mixing coef for t, q (mn,kk) at current t lvl   (m2/s2)
!     zl : height(kk) over monin-obukhov scale height.
!          a stability indicator.                   (mn)
!   sfcw : surface wind speed        (mn)                      (m/s)
!   ustar: friction velocity         (mn)                      (m/s)
!   tstar: potential temp scale      (mn)                      (k)
!   qstar: moisture scale            (mn)                      (kg/kg)
! itypbl : index for surface layer update, =0 stress b,c., =1 direct
!
!
! 4. local work arrays
!
!   wkk1 : working array               (mn,kk)
!   wkk2 : working array               (mn,kk)
!   wktri: tri-diagonal matrix solver working array   (mn,kk,4)
!
! 5. calling modules
!
!     pbltke
!
! 6. usage
!
!     call mixpbl( mn,kk,ktpbl,dt,g,hgt,u,v,t,q,ut,vt,tt,qt,e,eps
!    1           , xkm,xkh,zl,sfcw,ustar,tstar,qstar,itypbl )
!
! 7. modules called
!
!     couvtq
!     trigau
!     eeps
!
! 8. limitation
!
!   ** ktpbl  >= 2
!   **     e  >= 1.0e-4,    <= 250.0
!   **    eps >= 1.0e-7,    <= 1.0
!   **    xkm >=0.01 (stable) or 1.0 (unstable);  <= 1200.0
!   **    xkh = xhk*phih/phit
!
! 9. date
!
!     rewritten          november 1991
!
! 10. author
!
!      c. liou / f. wang / s. chang
!      modify to f90 by C-H Lee abd sort by River Chen in 2015
!
! 11. reference
!       detering, h. w., and d. etling, 1984: application of the
!            e-eps turbulence model to the atmospheric boundary
!            layer. boundary-layer meteorol., 33, 113-133.
!
!#####################################################################
      use paramt
      use const, only : RTYPE

      implicit  none
      integer   nxj,mn,kk,ktpbl,itypbl
      real      dt,g
!
! input & output variables
!
      real      hgt(mn,kk),u(mn,kk),v(mn,kk),t(mn,kk),           &
                e(mn,kk),                                        &
                eps(mn,kk),xkm(mn,kk),xkh(mn,kk),zl(mn),         &
                sfcw(mn),ustar(mn),tstar(mn),qstar(mn),          &
                dhgt(mn,kk),ro2(mn,kk),dhgtz(mn,kk)
      real(kind=RTYPE) q(mn,kk),qt(mn,kk),tt(mn,kk),             &
                       ut(mn,kk),vt(mn,kk)
!
!  local work arrays
!
      real      wkk1(im,lm),wkk2(im,lm)
      real      wktri(mn,kk,4)
!
      integer  k,i,kk1,kt,ktpbl1,nlimt,nlimq,nee,nc
      real     ec1,ec2,ec3,ec4,ec5,wgtx,wgt,umin,smin
      real     duu,dvv,beta,tmin,qmin,dth,dtee

      data ec1,ec2,ec3,ec4,ec5 / 2.0, 0.026, 1.38, 1.9, 0.77 /
!
!     compute eddy mixing coef xkm, xkh using current time lvl e, eps.
!
      do 100 k = 1, kk
      do 100 i = 1, nxj
      xkm(i,k) = 0.0
      xkh(i,k) = 0.0
  100 continue
!
      kk1 = kk - 1
      ktpbl1 = ktpbl - 1
      wgtx = 1.0/((kk-ktpbl1)*(kk-ktpbl1))
!
!     put upper and lower limits on xkm (depneding on lvls, stability)
!
      do 120 k = ktpbl, kk
      wgt = (k-ktpbl1)*(k-ktpbl1)*wgtx
      umin = 1.00 * wgt
      smin = 0.01 * wgt
      do 120 i = 1, nxj
      xkm(i,k) = ec2*e(i,k)*e(i,k)/eps(i,k)
!zzzz if (xkm(i,k) .gt. 1000.0)  xkm(i,k) = 1000.0
      if (xkm(i,k) .gt. 1200.0)  xkm(i,k) = 1200.0
      if ((t(i,k-1).le.t(i,k)).and.(xkm(i,k-1).lt.umin)) &
                                    xkm(i,k-1) =  umin
      if ((t(i,k-1).gt.t(i,k)).and.(xkm(i,k-1).lt.smin)) &
                                    xkm(i,k-1) =  smin
  120 continue
!
      do 200 k = ktpbl1, kk1
!
!     (1) compute richardson number : wkk1(,1)
          do 220 i = 1, nxj
          duu = u(i,k) - u(i,k+1)
          dvv = v(i,k) - v(i,k+1)
          wktri(i,1,1) = duu*duu + dvv*dvv
          if ( wktri(i,1,1) .lt. 1.0e-6 )  wktri(i,1,1) = 1.0e-6
! 220     continue
!         do 240 i = 1, nxj
          wkk1(i,1) = g*(t(i,k)-t(i,k+1))*(hgt(i,k)-hgt(i,k+1)) / &
                  ( 0.5*(t(i,k)+t(i,k+1))*wktri(i,1,1) )
! 240     continue
!
!     (2) compute phim/phih : wkk2(,1)
!         do 260 i = 1, nxj
          wktri(i,1,1) = zl(i)*0.5*(hgt(i,k)+hgt(i,k+1)) / hgt(i,kk)
          if( wktri(i,1,1) .gt. 10.0 )  wktri(i,1,1) = 10.0
          if( wktri(i,1,1) .lt. -10.0 ) wktri(i,1,1) = -10.0
! 260     continue
!         do 280 i = 1, nxj
!         if ( zl(i) .gt. 0.0 )  then
!            wkk2(i,1) =(1.0+4.7*wktri(i,1,1))/(0.74+4.7*wktri(i,1,1))
!         else
!            wkk2(i,1) =1.35*sqrt( (1.0-9.0*wktri(i,1,1))/ 
!    1                         sqrt(1.0-15.0*wktri(i,1,1)) )
!         endif
! ---     limit phim/phih = 1.35, if ri > 0.25
!         if ( wkk1(i,1) .gt. 0.25 ) wkk2(i,1) = 1.35
! 280     continue
!
! yamada(1983) formula of xkh/xkm
!       do 280 i=1,nxj
         if(wkk1(i,1) .lt. 0.16)then
           wkk2(i,1)=1.318*(0.2231-wkk1(i,1))/(0.2341-wkk1(i,1))
         else
           wkk2(i,1)=1.12
         endif
!280    continue
!       do 290 i = 1, nxj
        xkh(i,k) = xkm(i,k)*wkk2(i,1)
! 290   continue
!
  220   continue
  200 continue
      do 295 i = 1, nxj
      xkh(i,kk) = 1.35*xkm(i,kk)
  295 continue
!
!     update ut, vt, t, qt due to vertical mixing using backward
!     time scheme by solving a tri-diagonal matrix:
!     if itypbl= 0 (.ne. 1) use surface stress as lower b.c. in forcing
!
! (a) for ut  (assume boussineq approximation)
!     a*u(k+1,t+1) + b*u(k,t+1) + c*u(k-1,t+1) = u(k,t-1)
!     a= -dt*xkm(k+1/2) / (dz(k+1/2)*dz(k)),
!     b= 1 + dt*(xkm(k+1/2)/dz(k+1/2) + xkm(k-1/2)/dz(k-1/2)) / (dz(k)
!     c= -dt*xkm(k-1)/(dz(k-1/2)*dz(k)),    dt = (t+1) - (t-1),
!     dz(k) = 1/2*(dz(k+1/2)+dz(k-1/2)) = 1/2(z(k-1)-z(k+1))
!     k: odd level index (u lvl, k+1/2 --> k for xkm)
!
      kt = kk - ktpbl + 1
      beta = 1.0
!
      call couvtq ( nxj,mn,kk,ktpbl,dt,wktri,xkm,hgt,ut,beta,dhgt,ro2 )
      if ( itypbl .ne. 1 )  then
         do 300 i = 1, nxj
!        sfcwx = max( 1.0e-4, sqrt(ut(i,kk)**2 + vt(i,kk)**2) )
         wkk1(i,1) = ustar(i)/sfcw(i)
!org     wkk1(i,2) = dt*ustar(i)/( 0.5*(hgt(i,kk-1)+hgt(i,kk)) )
         wkk1(i,2) = dt*ustar(i)/( dhgt(i,kk) )
         wktri(i,kt,2) = wktri(i,kt,2) + wkk1(i,1)*wkk1(i,2)*ro2(i,kk)
  300    continue
      endif
      call trigau ( wktri, mn, nxj, kk, kt, 0 )
      do 320 k = ktpbl, kk
      do 320 i = 1, nxj
      ut(i,k) = wktri(i,k-ktpbl+1,4)
!cc   do 320 ik = (ktpbl-1)*mn+1, kk*mn
!cc   ut(ik,1) = wktri(ik-(ktpbl-1)*mn,1,4)
  320 continue
!
! (b) for vt
!
      call couvtq ( nxj,mn,kk,ktpbl,dt,wktri,xkm,hgt,vt,beta,dhgt,ro2 )
      if ( itypbl .ne. 1 )  then
         do 330 i = 1, nxj
         wktri(i,kt,2) = wktri(i,kt,2) + wkk1(i,1)*wkk1(i,2)*ro2(i,kk)
  330    continue
      endif
      call trigau ( wktri, mn, nxj, kk, kt, 0 )
      do 340 k = ktpbl, kk
      do 340 i = 1, nxj
      vt(i,k) = wktri(i,k-ktpbl+1,4)
!cc   do 340 ik = (ktpbl-1)*mn+1, kk*mn
!cc   vt(ik,1) = wktri(ik-(ktpbl-1)*mn,1,4)
  340 continue
!
! (c) for tt
!     limit surface forcing effects to less than 10% of current
!     values to avoid too large forcing
!
      call couvtq ( nxj,mn,kk,ktpbl,dt,wktri,xkh,hgt,tt,beta,dhgt,ro2 )
      if ( itypbl .ne. 1 )  then
       nlimt=0
         do 350 i = 1, nxj
!
         tmin = wktri(i,kt,4)*0.9
         wktri(i,kt,4) = wktri(i,kt,4) - wkk1(i,2)*tstar(i)*ro2(i,kk)
         wktri(i,kt,4) = max ( tmin, wktri(i,kt,4) )
         if(wktri(i,kt,4) .eq. tmin)nlimt=nlimt+1
  350    continue
      endif
      call trigau ( wktri, mn, nxj, kk, kt, 0 )
      do 360 k = ktpbl, kk
      do 360 i = 1, nxj
      tt(i,k) = wktri(i,k-ktpbl+1,4)
!cc   do 360 ik = (ktpbl-1)*mn+1, kk*mn
!cc   tt(ik,1) = wktri(ik-(ktpbl-1)*mn,1,4)
  360 continue
!
! (d) for qt
!     limit surface forcing effects to less than 20% of current
!     values to avoid too large forcing
!
      call couvtq ( nxj,mn,kk,ktpbl,dt,wktri,xkh,hgt,qt,beta,dhgt,ro2 )
      if ( itypbl .ne. 1 )  then
       nlimq=0
         do 370 i = 1, nxj
         qmin = wktri(i,kt,4)*0.8
         wktri(i,kt,4) = wktri(i,kt,4) - wkk1(i,2)*qstar(i)*ro2(i,kk)
         wktri(i,kt,4) = max ( qmin, wktri(i,kt,4) )
         if(wktri(i,kt,4) .eq. qmin)nlimq=nlimq+1
  370    continue
      endif
      call trigau ( wktri, mn, nxj, kk, kt, 0 )
      do 380 k = ktpbl, kk
      do 380 i = 1, nxj
      qt(i,k) = wktri(i,k-ktpbl+1,4)
  380 continue
!
!     if(j.eq.99)print*,'after matrix qq(kk)=',qt(87,kk)
!    +           ,qt(87,kk-1),qt(87,kk-2)
!     update e, eps from current time level to future time level
!
      dth = dt
!  -- limit time step for e eps prediction to be less than 150 sec
!
      nee = 1.0 + (dth-1.0)/150.0
      dtee = dth / nee
      do 400 nc = 1, nee
      call eeps ( nxj,mn,kk,ktpbl,g,dtee,ec1,ec2,ec3,ec4,ec5,u,v,t &
                , hgt,zl,e,eps,xkm,xkh,wkk1,wkk2,wktri,dhgtz )
  400 continue
!
      return
      end
