      subroutine eeps ( nxj,mn,kk,ktpbl,g,dt,ec1,ec2,ec3,ec4,ec5,u,v,t &
                      , hgt,zl,e,eps,xkm,xkh,shr,wkk2,wktri,dhgt )
!##################################################################
!                     subroutine description
!
! 1. function description
!
!       this subroutine predicts e and eps in the future time level
!       for vertical levels from kk-1 to ktpbl.
!       ( e and eps at kk levels are diagnosed in sfcpbl )
!
! 2. common block specification
!
!
! 3. parameter specification
!
!     mn : dimension of horizontal direction
!     kk : dimension of vertical direction
!  ktpbl : starting level to consider vertical mixing effects
!     dt : time step for e and eps updates
!      g : gravity                                               (m/s2)
!     ec1: coef for   e-prediction equation, eq. (15) of gerber (1989)
!     ec2: constant to compute xkm from e**2/eps
! ec3,4,5: coef for eps-prediction equation, eq. (16) of gerber (1989)
!      u : x-direction velocity (m,n,kk) at the current time lvl (m/s)
!      v : y-direction velocity (m,n,kk) at the current time lvl (m/s)
!      t : temp.     on odd lvl (m,n,kk) at the current time lvl (k)
!    hgt : heigh above ground   (m,n,kk) at the current time lvl (m)
!     zl : lowest sigma level height over mixing length  (mn)
!     e  : turbulence kinetic energy (mn,kk) at current/future   (m2/s2)
!    eps : tke dissipation           (mn,kk) at current/future   (m2/s2)
!    xkm : momentum transport coef   (mn,kk) at current time lvl (m2/s2)
!    xkh : heat/moist transpt coef   (mn,kk) at current time lvl (m2/s2)
!
! 4. work arrays
!    shr : working array for coef. calculation            (mn,kk)
!   wkk2 : working array for coef. calculation            (mn,kk)
!  wktri : working array for tri-diagonal system solver   (mn,kk,4)
!
! 4.1 local work array
!     eo : working array to save e at previous time lvl   (mn,kk)
!
! 5. calling modules
!
!      mixpbl
!
! 6. usage
!
!     call eeps ( mn,kk,ktpbl,g,dt,ec1,ec2,ec3,ec4,ec5,u,v,t
!    1          , hgt,zl,e,eps,xkm,xkh,wkk1,wkk2,wktri )
!
! 7. modules called
!
!     trigau
!
! 8. limitation
!
!         e   >= 1.0e-04, e   <= 250.0
!         eps >= 1.0e-07, eps <= 1.0
!
! 9. date
!
!     created              10/01/1989
!     rewritten            nov.  1991
!
!
! 10. author
!
!      c. liou / f. wang / c. chang
!      modify to f90 by C-H Lee and sort by River Chen in 2015
!
! 11. reference
!       detering h.,w., and d. etling, 1984: application of the
!            e-eps turbulence model to the atmospheric boundary
!            layer. boundary-layer meteorol., 33, 113-133.
!
!#####################################################################
!
      use paramt

      implicit  none
      integer   nxj,mn,kk,ktpbl
      real      g,dt,ec1,ec2,ec3,ec4,ec5

!  input & output variable (including work arrays)
      real      u(mn,kk),v(mn,kk),t(mn,kk),hgt(mn,kk),zl(mn), &
                e(mn,kk),eps(mn,kk),xkm(mn,kk),xkh(mn,kk),    &
                shr(mn,kk),wkk2(mn,kk),wktri(mn,kk,4)
      real      dhgt(mn,kk)
!
!  local work array
!
      real      eo(im,lm)
      real      vdiff,tflt
      data      vdiff /0.1/, tflt /0.2/

      integer   kk1,kt1,kt2,mnkt1,kpp1,kpm1,ik,iksp,iks,iksm, &
                k,i,ktpbx,kppx,kpmx
      real      duu,dvv,dzz,dzzu,dzzup,buoy,eee,xkmnew,factor

!
! *note: e, eps prediction are only from ktpbl to kk-1 levels
!
      kk1   = kk - 1
      kt1   = kk1 - ktpbl + 1
      kt2   = kt1 - 1
      mnkt1 = mn*kt1
      kpp1  = ktpbl + 1
      kpm1  = ktpbl - 1
!
!-----------------------------------------------------------------
! (1) prepare tri-diagonal system coef for e-prediction equation -
!-----------------------------------------------------------------
!
      do ik = 1, nxj
      iksp=(ktpbl  )*mn+ik
      iks =(ktpbl-1)*mn+ik
      iksm=(ktpbl-2)*mn+ik
      duu = u(iks,1) - u(iksp,1)
      dvv = v(iks,1) - v(iksp,1)
      dzz = hgt(iks,1) - hgt(iksp,1)
!cc   dzzu = 0.5*( hgt(iksm,1) - hgt(iksp,1) )
      dzzu = dhgt(iks,1)
      dzzup= dhgt(iksp,1)
      shr(ik,1) = (duu*duu+dvv*dvv)*xkm(iks,1)/dzz
      buoy= -g*xkh(iks,1)*(t(iks,1)-t(iksp,1))/(0.5*(t(iks,1) &
            +t(iksp,1)))
      wktri(ik,1,1) = - dt*ec1*0.5*(xkm(iks,1)+xkm(iksm,1))/dzzu
      wktri(ik,1,3) = - dt*ec1*0.5*(xkm(iks,1)+xkm(iksp,1)) &
                      /dzzup
      wktri(ik,1,2) = dzz - wktri(ik,1,1) - wktri(ik,1,3)
      wktri(ik,1,4) = dzz*(e(iks,1)-dt*eps(iks,1)) &
            + dt*(shr(ik,1)+buoy)
      enddo

      do 120 k = 2, kt1
      do 120 i = 1, nxj
      ik = (k-1)*mn + i
      iksp=(ktpbl  )*mn+ik
      iks =(ktpbl-1)*mn+ik
      iksm=(ktpbl-2)*mn+ik
      duu = u(iks,1) - u(iksp,1)
      dvv = v(iks,1) - v(iksp,1)
      dzz = hgt(iks,1) - hgt(iksp,1)
!cc   dzzu = 0.5*( hgt(iksm,1) - hgt(iksp,1) )
      dzzu = dhgt(iks,1)
      dzzup= dhgt(iksp,1)
      shr(ik,1) = (duu*duu+dvv*dvv)*xkm(iks,1)/dzz
      buoy= -g*xkh(iks,1)*(t(iks,1)-t(iksp,1))/(0.5*(t(iks,1) &
            +t(iksp,1)))
      wktri(ik,1,1) = wktri(ik-mn,1,3)
      wktri(ik,1,3) = - dt*ec1*0.5*(xkm(iks,1)+xkm(iksp,1)) &
                      /dzzup
      wktri(ik,1,2) = dzz - wktri(ik,1,1) - wktri(ik,1,3)
      wktri(ik,1,4) = dzz*(e(iks,1)-dt*eps(iks,1)) &
            + dt*(shr(ik,1)+buoy)
  120 continue
!
!     -- do not allow   dt*(total damping) < e(t)   --
!
      do 140 k = 1, kt1
      do 140 i = 1, nxj
      if ( wktri(i,k,4) .lt. 0.0 )  wktri(i,k,4) = 0.0
  140 continue
!
!
!     fix upper b.c (move -e(ktpbl-1)*wktri(1,1) to forcing term)
!
      do 150 i = 1, nxj
      wktri(i,1,4) = wktri(i,1,4) - e(i,ktpbl-1)*wktri(i,1,1)
  150 continue
!
!     put surface forcing to kk level
!
      do 160 i = 1, nxj
      wktri(i,kt1,4) = wktri(i,kt1,4) - e(i,kk)*wktri(i,kt1,3)
  160 continue
!
!     update e to future time level by backward time scheme
!
      call trigau ( wktri, mn, nxj, kk, kt1, 0 )
!
      do 180 k = 1, kt1
      ktpbx=ktpbl+k-1
      do 180 i = 1, nxj
      eo(i,ktpbx) = e(i,ktpbx)
      e(i,ktpbx)  = wktri(i,k,4)
  180 continue
!
!     apply vertical diffusion to e
!
      do 185 k = 1, kt1
      kppx=kpp1+k-1
      kpmx=kpm1+k-1
      ktpbx=ktpbl+k-1
      do 185 i = 1, nxj
      wktri(i,k,4) = 0.5*vdiff*(e(i,kppx)-2.0*e(i,ktpbx)+e(i,kpmx))
  185 continue
      do 187 k = 1, kt1
      ktpbx=ktpbl+k-1
      do 187 i = 1, nxj
      e(i,ktpbx) = e(i,ktpbx) + wktri(i,k,4)
  187 continue
!
!     limit e to 1.0e-04 and 250.0
!
      do 190 k = 1, kt1
      ktpbx=ktpbl+k-1
      do 190 i = 1, nxj
      if ( e(i,ktpbx) .lt. 1.0e-04 )  e(i,ktpbx)  = 1.0e-04
      if ( e(i,ktpbx) .gt. 250.0   )  e(i,ktpbx)  = 250.0
  190 continue
!
!---------------------------------------------------------------------
! (2) prepare tri-diagonal system coef for eps-prediction equation   -
!---------------------------------------------------------------------
!
!     enhance shear effects, if buoyance > 0
!
      do 220 k = 1, kt1
      do 220 i = 1, nxj
      ik = (k-1)*mn + i
      iks =(ktpbl-1)*mn+ik
      iksp=(ktpbl  )*mn+ik
      if( t(iks,1).lt.t(iksp,1) )  shr(ik,1) = shr(ik,1) +  &
        (-g)*xkh(iks,1)*(t(iks,1)-t(iksp,1))                &
        /(0.5*(t(iks,1)+t(iksp,1)))
  220 continue
!
      do ik = 1, nxj
      iks =(ktpbl-1)*mn+ik
      iksp=(ktpbl  )*mn+ik
      iksm=(ktpbl-2)*mn+ik
      dzz = hgt(iks,1) - hgt(iksp,1)
!cc   dzzu = 0.5*( hgt(iksm,1) - hgt(iksp,1) )
      dzzu = dhgt(iks,1)
      dzzup= dhgt(iksp,1)
      eee  = 0.5*( e(iks,1) + eo(iks,1) )
      wktri(ik,1,1) = - dt*ec5*0.5*(xkm(iks,1)+xkm(iksm,1)) &
                      /dzzu
      wktri(ik,1,3) = - dt*ec5*0.5*(xkm(iks,1)+xkm(iksp,1)) &
                      /dzzup
      wktri(ik,1,2) = dzz + dzz*dt*ec4*eps(iks,1)/eee &
                     - wktri(ik,1,1) - wktri(ik,1,3)
      wktri(ik,1,4) = eps(iks,1)*dzz +                &
                           dt*ec3*ec2*eee*shr(ik,1)/xkm(iks,1)
      enddo

      do 240 k = 2, kt1
      do 240 i = 1, nxj
      ik = (k-1)*mn + i
      iks =(ktpbl-1)*mn+ik
      iksp=(ktpbl  )*mn+ik
      iksm=(ktpbl-2)*mn+ik
      dzz = hgt(iks,1) - hgt(iksp,1)
!cc   dzzu = 0.5*( hgt(iksm,1) - hgt(iksp,1) )
      dzzu = dhgt(iks,1)
      dzzup= dhgt(iksp,1)
      eee  = 0.5*( e(iks,1) + eo(iks,1) )
      wktri(ik,1,1) = wktri(ik-mn,1,3)
      wktri(ik,1,3) = - dt*ec5*0.5*(xkm(iks,1)+xkm(iksp,1)) &
                      /dzzup
      wktri(ik,1,2) = dzz + dzz*dt*ec4*eps(iks,1)/eee &
                     - wktri(ik,1,1) - wktri(ik,1,3)
      wktri(ik,1,4) = eps(iks,1)*dzz +                &
                           dt*ec3*ec2*eee*shr(ik,1)/xkm(iks,1)
  240 continue
!
!     fix upper b.c (move -eps(ktpbl-1)*wktri(1,1) to forcing term)
!
      do 250 i = 1, nxj
      wktri(i,1,4) = wktri(i,1,4) - eps(i,ktpbl-1)*wktri(i,1,1)
  250 continue
!
!     put surface forcing to kk level
!
      do 260 i = 1, nxj
      wktri(i,kt1,4) = wktri(i,kt1,4) - wktri(i,kt1,3)*eps(i,kk)
  260 continue
!
!     update eps to future time level by backward time scheme
!
      call trigau ( wktri, mn, nxj, kk, kt1, 0 )
!
      do 280 k = 1, kt1
      ktpbx=ktpbl+k-1
      do 280 i = 1, nxj
      wkk2(i,ktpbx) = eps(i,ktpbx)
      eps(i,ktpbx)  = wktri(i,k,4)
  280 continue
!
!     apply vertical diffusion to eps
!
      do 285 k = 1, kt1
      kppx=kpp1+k-1
      kpmx=kpm1+k-1
      ktpbx=ktpbl+k-1
      do 285 i = 1, nxj
      wktri(i,k,4) = 0.5*vdiff*(eps(i,kppx)-2.0*eps(i,ktpbx) &
                                               +eps(i,kpmx))
  285 continue
      do 287 k = 1, kt1
      ktpbx=ktpbl+k-1
      do 287 i = 1, nxj
      eps(i,ktpbx) = eps(i,ktpbx) + wktri(i,k,4)
  287 continue
!
!     limit eps to 1.0e-07 and 1.0
!
      do 290 k = 1, kt1
      ktpbx=ktpbl+k-1
      do 290 i = 1, nxj
      if ( eps(i,ktpbx) .lt. 1.0e-07 )  eps(i,ktpbx) = 1.0e-07
      if ( eps(i,ktpbx) .gt. 1.0     )  eps(i,ktpbx) = 1.0
  290 continue
!
!     apply two level time filter on e and eps to avoid oscillation
!
      do 300 k = 1, kt1
      ktpbx=ktpbl+k-1
      do 300 i = 1, nxj
      e(i,ktpbx)   = (1.0-tflt)*e(i,ktpbx)   + tflt*eo(i,ktpbx)
      eps(i,ktpbx) = (1.0-tflt)*eps(i,ktpbx) + tflt*wkk2(i,ktpbx)
  300 continue
!
!     re-define xkm, xkh
!
      do 400 k = 1, kt1
      ktpbx=ktpbl+k-1
      do 400 i = 1, nxj
      xkmnew = ec2*e(i,ktpbx)*e(i,ktpbx)/eps(i,ktpbx)
      if ( xkmnew .gt. 1200.0 )  xkmnew = 1200.0
      factor = xkmnew / xkm(i,ktpbx)
      xkm(i,ktpbx) = xkmnew
      xkh(i,ktpbx) = xkh(i,ktpbx) * factor
  400 continue
      return
      end
