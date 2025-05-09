      subroutine rozphys(nxj, nx, lev, dt, iter, theta, julian,        &
                         o3l, tt, pp, ps, myrank)
      use machine,  only : kind_phys
      use physcons, only : grav => con_g
      use ozne_def
      use const, only : RTYPE
      integer nxj,nx,lev,julian 
      real (kind=kind_phys) ozwk1(latsozp,levozp,pl_coeff)
      real (kind=kind_phys) ozwk2(levozp,pl_coeff)
      real (kind=kind_phys) ozplout(nx,lev,pl_coeff)
      real (kind=kind_phys) amin(pl_coeff)
      real (kind=kind_phys) amax(pl_coeff)
        
!
      real, parameter :: gravi=1.0/grav
      integer pl_coeff2, kmax(pl_coeff),kmin(pl_coeff)
             

      real(kind=kind_phys) pp(nx,lev),ozp(nx,lev,pl_coeff)
      real(kind=RTYPE)     o3l(nx,lev),tt(nx,lev),ps(nx)
      real(kind=kind_phys) dt
!
      integer k,i,j
      logical ldiag3d
      real(kind=kind_phys) coef(nx,pl_coeff),                    &
                           ozb(nx),  colo3(nx,lev),              &
                           ozo(nx,lev), delp(nx,lev),temp
      ozo=0.
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
             ozwk1(:,:,:) = con2*ozplin(:,:,:,n) + con1*ozplin(:,:,:,n1)
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

!
!  linear interpolate for p
!
!ch      do  300 j = 1, pl_coeff
!ch      do  310 k = 1, lev
!ch      do  320 i = 1, nxj
!ch          ozplout(i,k,j)=0.
!ch          if ( pp(i,k) .le. pl_pres(1) )  ozplout(i,k,j) = ozwk2(1,j)
!ch          if ( pp(i,k) .ge. pl_pres(levozp) ) ozplout(i,k,j) = ozwk2(levozp,j)
!ch 320  continue
!ch 310  continue
!ch 300  continue
!ch >>
      do  310 k = 1, lev
      do  320 i = 1, nxj
          if ( pp(i,k) .le. pl_pres(1) )then
             do j = 1, pl_coeff
               ozplout(i,k,j) = ozwk2(1,j)
             enddo
          elseif ( pp(i,k) .ge. pl_pres(levozp) )then
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
      do 380 i  = 1, nxj
         con1 = (pp(i,kk)-pl_pres(k)) / (pl_pres(k+1)-pl_pres(k)) 
         if ((con1.gt.0.0) .and. (con1.le.1.0)) then !pp(kk) in pl(k,k+1)
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
!     do i=1,nx*lev*pl_coeff
!        ozp(i,1,1)=0.0
!     enddo
         ozp=0.0

      pl_coeff2 = 2
!cmy--------------------------------------------------------------------
!for pl_coeff > 2
!cmy--------------------------------------------------------------------
      if (pl_coeff .gt. 2) then
          colo3(:,lev) = 0.0
!     
          do k=2,lev
          do i=1,nxj
             delp(i,k)=pp(i,k)-pp(i,k-1)
          enddo
          enddo
!
          do i=1,nxj
             delp(i,1)=pp(i,1)
             colo3(i,1)=o3l(i,1) * delp(i,1) * gravi
          enddo
!
          do k=2,lev
          do i=1,nxj
             colo3(i,k) = colo3(i,k-1) + o3l(i,k) * delp(i,k) * gravi
          enddo
          enddo
      endif ! for pl_coeff > 2
!cmy--------------------------------------------------------------------
!  K - loop start
!cmy--------------------------------------------------------------------
      do k=1,lev
        do j=1,pl_coeff
          do i=1,nxj
              coef(i,j) = ozplout(i,k,j)
          enddo
        enddo
!cmy--------------------------------------------------------------------
        if (pl_coeff2 .eq. 2) then 
          do i=1,nxj
            ozb(i)    = o3l(i,k)           ! NO FilliNG
            ozo(i,k)  = (ozb(i) + coef(i,1)*dt) / (1.0 + coef(i,2)*dt)
          enddo
          if (ldiag3d) then     !     Ozone change diagnostics
              do i=1,nxj
                 ozp(i,k,1) = ozp(i,k,1) + coef(i,1)*DT
                 ozp(i,k,2) = ozp(i,k,2) + (ozo(i,k) - ozb(i))
              enddo
          endif ! for ldiag3d
        endif ! for pl_coeff2=2
!cmy--------------------------------------------------------------------
        if (pl_coeff2 .eq. 4) then 
          do i=1,nxj
            ozb(i)   = o3l(i,k)            ! NO FilliNG
            temp     = coef(i,1) + coef(i,3)*tt(i,k) + coef(i,4)*colo3(i,k)
            ozo(i,k) = (OZB(i)  + temp*dt) / (1.0 + coef(i,2)*dt)
          enddo
          if (ldiag3d) then     !     Ozone change diagnostics
              do i=1,nxj
                 OZP(i,k,1) = OZP(i,k,1) + coef(i,1)*DT
                 OZP(i,k,2) = OZP(i,k,2) + (OZO(i,k)-OZB(i))
                 OZP(i,k,3) = OZP(i,k,3) + coef(i,3)*tt(i,k)*DT
                 OZP(i,k,4) = OZP(i,k,4) + coef(i,4)*colo3(i,k+1)*DT
              enddo
          endif !for ldiag3d
      endif !for pl_coeff2=4
!cmy------------------------------------------------------------------
      enddo ! Vertical loop
!cmy------------------------------------------------------------------
!     do i=1,nx*lev
!      o3l(i,1)=ozo(i,1)
!     enddo
       o3l=ozo
!cmy------------------------------------------------------------------
!     if (myrank .eq. 0) print *,'*** ozphys.f end !!' 
!     if (myrank .eq. 0) print *,'*** o3l ***'
!     call qmax2d(o3l,1,1,nx,lev)

      RETURN
      END
