      subroutine ozphys_2015 (ix, im, lev,  dt, theta, julian, ozi, &
                              ozo, tin, prsl, delp, me)
!                              prsl, ozplout, pl_coeff, delp, ldiag3d,&
!                              ozp,me)
!
!     this code assumes that both prsl and po3 are from bottom to top
!     as are all other variables
!     This code is specifically for NRL parameterization and
!     climatological T and O3 are in location 5 and 6 of ozplout array
! June 2015 - Shrinivas Moorthi
!
      use machine , only : kind_phys
      use physcons, only : grav => con_g
      use ozne_def
      use const,    only : RTYPE
      implicit none
!
      real, parameter :: gravi=1.0/grav
      integer im, ix, lev, me, julian
      real    theta,dt
      real(kind=kind_phys) po3(levozp),  prsl(ix,lev),  delp(ix,lev)
      real(kind=RTYPE)     tin(ix,lev), ozi(ix,lev), ozo(ix,lev)
!                          ozp(ix,lev,4),  dt
!
      integer k,kk,kmax,kmin,l,ll,i,j,n,n1,lozc,daynum
      real     con1, con2
      logical              ldiag3d, flg(im)
      real(kind=kind_phys) ozwk1(latsozp,levozp,pl_coeff)
      real(kind=kind_phys) ozwk2(levozp,pl_coeff)
      real(kind=kind_phys) ozplout(im,lev,pl_coeff)
      real(kind=kind_phys) pmax, pmin, tem, temp
      real(kind=kind_phys) wk1(im), wk2(im), wk3(im),                  &
                           prod(im,pl_coeff),                          &
                           ozib(im), colo3(im,lev+1), coloz(im,lev+1)
!----------------------------------------------------------------------
      lozc = latsozp*levozp*pl_coeff
!----------------------------------------------------------------------
!  check julian locating on which month is 
!  and linearly interpolate for time
!----------------------------------------------------------------------
      daynum = julian
      if ( daynum .gt. 365 )     daynum = 365.
      if ( daynum .le. pl_time(1) ) daynum = julian+365.

      do n = 1,timeoz
         con1 = (daynum-pl_time(n)) / (pl_time(n+1)-pl_time(n))
         n1=n+1
         if (n1 .gt. 12) n1=n1-12
         if ( (con1.ge.0.0) .and. (con1.le.1.0) ) then
            con2=1.0-con1
            do i=1,lozc
               ozwk1(i,1,1) = con2*ozplin(i,1,1,n) + con1*ozplin(i,1,1,n1)
            enddo
         endif
       enddo
!
!  linearly interpolate for latitude
!
!ch      if ((theta.le.pl_lat(1)) .or. (theta.ge.pl_lat(latsozp))) then
!ch      do 200 i = 1,pl_coeff
!ch      do 210 k = 1,levozp
!ch         ozwk2(k,i)=0.
!ch         if(theta .le. pl_lat(1)) ozwk2(k,i) = ozwk1(1,k,i)
!ch         if(theta .ge. pl_lat(latsozp)) ozwk2(k,i) = ozwk1(latsozp,k,i)
!ch 210  continue
!ch 200  continue
!ch >>
      if (theta.le.pl_lat(1))then
      do i = 1,pl_coeff
      do k = 1,levozp
         ozwk2(k,i) = ozwk1(1,k,i)
      enddo
      enddo
      elseif(theta.ge.pl_lat(latsozp)) then
      do i = 1,pl_coeff
      do k = 1,levozp
         ozwk2(k,i) = ozwk1(latsozp,k,i)
      enddo
      enddo
!ch <<

      else

      ll  = (theta-pl_lat(1))/5.0 + 1 
      con1 = (theta-pl_lat(ll)) / 5.0
      con2 = 1.0-con1
      do 220 i = 1,pl_coeff
      do 230 k = 1,levozp
        ozwk2(k,i) = con2*ozwk1(ll,k,i) + con1*ozwk1(ll+1,k,i)
 230  continue
 220  continue
      endif

!
!  linear interpolate for p
!
!ch      do  300 j = 1, pl_coeff
!ch      do  310 k = 1, lev
!ch      do  320 i = 1, im
!ch          ozplout(i,k,j)=0.
!ch          if ( prsl(i,k) .le. pl_pres(1) )  ozplout(i,k,j) = ozwk2(1,j)
!ch          if ( prsl(i,k) .ge. pl_pres(levozp) ) ozplout(i,k,j) = ozwk2(levozp,j)
!ch 320  continue
!ch 310  continue
!ch 300  continue
!ch >>
      do  310 k = 1, lev
      do  320 i = 1, im 
          if ( prsl(i,k) .le. pl_pres(1) )then
             do j = 1, pl_coeff
               ozplout(i,k,j) = ozwk2(1,j)
             enddo
          elseif ( prsl(i,k) .ge. pl_pres(levozp) )then
             do j = 1, pl_coeff
               ozplout(i,k,j) = ozwk2(levozp,j)
             enddo
          else
             do j = 1, pl_coeff
               ozplout(i,k,j) = 0.
             enddo
          endif
 320  continue
 310  continue
!ch <<

!
      do 350 j  = 1, pl_coeff
      do 360 kk = 1, lev
      do 370 k  = 1, levozp-1
!!ocl simd
      do 380 i  = 1, im
         con1 = (prsl(i,kk)-pl_pres(k)) / (pl_pres(k+1)-pl_pres(k))
         if ((con1.gt.0.0) .and. (con1.le.1.0)) then !prsl(kk) in pl(k,k+1)
            con2=1.0-con1
            ozplout(i,kk,j) = con2*ozwk2(k,j) + con1*ozwk2(k+1,j)
         endif
 380  continue
 370  continue
 360  continue
 350  continue
!
!----------------------------------------------------------------------

      ldiag3d =.true.

!
        colo3(:,lev+1) = 0.0
        coloz(:,lev+1) = 0.0
!
      do l=lev,1,-1
        pmin =  1.0e10
        pmax = -1.0e10
!
        do i=1,im
!          wk1(i) = log(prsl(i,l))
          wk1(i) = prsl(i,l)
          pmin   = min(wk1(i), pmin)
          pmax   = max(wk1(i), pmax)
          prod(i,:) = 0.0
        enddo
        kmax = 1
        kmin = 1
        do k=levozp-1,1,-1
          if (pmin < pl_pres(k)) kmax = k
          if (pmax < pl_pres(k)) kmin = k
        enddo
!
        do k=kmin,kmax
!          temp = 1.0 / (po3(k) - po3(k+1))
          temp = 1.0 / (pl_pres(k+1) - pl_pres(k))
          do i=1,im
            flg(i) = .false.
!            if (wk1(i) < po3(k) .and. wk1(i) >= po3(k+1)) then
            if (wk1(i) < pl_pres(k+1) .and. wk1(i) >= pl_pres(k)) then
              flg(i) = .true.
!              wk2(i) = (wk1(i) - po3(k+1)) * temp
              wk2(i) = (wk1(i) - pl_pres(k)) * temp
              wk3(i) = 1.0 - wk2(i)
            endif
          enddo
          do j=1,pl_coeff
            do i=1,im
              if (flg(i)) then
                prod(i,j)  = wk2(i) * ozplout(i,k+1,j)                 &
                           + wk3(i) * ozplout(i,k,j)
              endif
            enddo
          enddo
        enddo
!
        do j=1,pl_coeff
          do i=1,im
!            if (wk1(i) < po3(levozp)) then
            if (wk1(i) >= pl_pres(levozp)) then
              prod(i,j) = ozplout(i,levozp,j)
            endif
!            if (wk1(i) >= po3(1)) then
            if (wk1(i) < pl_pres(1)) then
              prod(i,j) = ozplout(i,1,j)
            endif
          enddo
        enddo
        do i=1,im
          colo3(i,l) = colo3(i,l+1) + ozi(i,l)  * delp(i,l)*gravi
          coloz(i,l) = coloz(i,l+1) + prod(i,6) * delp(i,l)*gravi
          prod(i,2)  = min(prod(i,2), 0.0)
        enddo
!       write(1000+me,*) ' colo3=',colo3(1,l),' coloz=',coloz(1,l)
!    &,' l=',l
        do i=1,im
          ozib(i)  = ozi(i,l)            ! no filling
          tem      = prod(i,1) - prod(i,2) * prod(i,6)                 &
                   + prod(i,3) * (tin(i,l) - prod(i,5))                &
                   + prod(i,4) * (colo3(i,l)-coloz(i,l))

!     if (me .eq. 0) print *,'ozphys_2015 tem=',tem,' prod=',prod(i,:)
!    &,' ozib=',ozib(i),' l=',l,' tin=',tin(i,l),'colo3=',colo3(i,l+1)

            ozo(i,l) = (ozib(i)  + tem*dt) / (1.0 - prod(i,2)*dt)
        enddo
!        if (ldiag3d) then     !     ozone change diagnostics
!          do i=1,im
!            ozp(i,l,1) = ozp(i,l,1) + (prod(i,1)-prod(i,2)*prod(i,6))*dt
!            ozp(i,l,2) = ozp(i,l,2) + (ozo(i,l) - ozib(i))
!            ozp(i,l,3) = ozp(i,l,3) + prod(i,3)*(tin(i,l)-prod(i,5))*dt
!            ozp(i,l,4) = ozp(i,l,4) + prod(i,4)
!     &                              * (colo3(i,l)-coloz(i,l))*dt
!          enddo
!        endif
      enddo                                ! vertical loop
!
      return
      end
