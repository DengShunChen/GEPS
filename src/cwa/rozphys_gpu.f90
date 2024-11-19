      subroutine rozphys_gpu(nxj, nx, lev, dt, iter, theta, julian,        &
                         o3l, tt, pp, ps, myrank)
      use physcons, only : grav => con_g
      use ozne_def
      use const, only : RTYPE
      use param, only : my, my_max
      use index, only: jlist1, jlistnum
      
      implicit none
      
      integer nxj(my),nx,lev,julian 
      real (kind=kind_phys) ozwk1(latsozp,levozp,pl_coeff)
      real (kind=kind_phys) ozwk2(levozp,pl_coeff,my_max)
      real (kind=kind_phys) ozplout(nx,lev,pl_coeff,my_max)
      real (kind=kind_phys) amin(pl_coeff)
      real (kind=kind_phys) amax(pl_coeff)
      integer :: iter, myrank, lozc, n, n1, ll,kk
      real :: theta(my), con1, con2, daynum
        
!
      real, parameter :: gravi=1.0/grav
      integer pl_coeff2, kmax(pl_coeff),kmin(pl_coeff)
             

      real(kind=kind_phys) ps(nx,my_max),                               &
                           pp(nx,lev,my_max),ozp(nx,lev,pl_coeff,my_max)
      real(kind=RTYPE)     o3l(nx,lev,my_max),tt(nx,lev,my_max)
      real(kind=kind_phys) dt
!
      integer k,i,j, jj, j1
      logical ldiag3d
      real(kind=kind_phys) coef(nx,lev,pl_coeff,my_max),                    &
                           ozb(nx,lev,my_max),  colo3(nx,lev,my_max),              &
                           ozo(nx,lev,my_max), delp(nx,lev,my_max),temp
      integer async_id
      real :: time2, time1
      async_id = 1
!----------------------------------------------------------------------
      lozc = latsozp*levozp*pl_coeff
!----------------------------------------------------------------------
!  check julian locating on which month is 
!  and linearly interpolate for time
!----------------------------------------------------------------------
      daynum = julian
      if ( daynum .gt. 365 )     daynum = 365.
      if ( daynum .le. pl_time(1) ) daynum = julian+365.
    
      !!$acc parallel loop seq private(con1,n,n1,i,con2)
      !do n = 1,timeoz
      !   con1 = (daynum-pl_time(n)) / (pl_time(n+1)-pl_time(n))
      !   n1=n+1
      !   if (n1 .gt. 12) n1=n1-12
      !   if ( (con1.ge.0.0) .and. (con1.le.1.0) ) then
      !      con2=1.0-con1
      !      do i=1,lozc
      !         ozwk1(i,1,1) = con2*ozplin(i,1,1,n) + con1*ozplin(i,1,1,n1)
      !      enddo
      !   endif
      !enddo

      !!$acc parallel loop private(i,n,n1,con1,con2)
      do i=1,lozc
        !!$acc loop seq
        do n = 1,timeoz
         con1 = (daynum-pl_time(n)) / (pl_time(n+1)-pl_time(n))
         n1=n+1
         if (n1 .gt. 12) n1=n1-12
         if ( (con1.ge.0.0) .and. (con1.le.1.0) ) then
            con2=1.0-con1
               ozwk1(i,1,1) = con2*ozplin(i,1,1,n) + con1*ozplin(i,1,1,n1)
         endif
        enddo
      enddo

      !$acc enter data copyin(jlist1, pl_lat, ozwk1, pl_pres) async(async_id)
      !$acc enter data create(ozwk2, ozplout, colo3, delp, ozb, ozo, ozp) async(async_id)
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
      !$acc parallel loop gang private(jj,i,k,j1,ll,con1,con2) async(async_id)
      do jj = 1, jlistnum
        j1 = jlist1(jj)

        if (theta(j1).le.pl_lat(1))then
          !$acc loop vector collapse(2)
          do i = 1,pl_coeff
            do k = 1,levozp
               ozwk2(k,i,jj) = ozwk1(1,k,i)
            enddo
          enddo
        elseif(theta(j1).ge.pl_lat(latsozp)) then
          !$acc loop vector collapse(2)
          do i = 1,pl_coeff
            do k = 1,levozp
               ozwk2(k,i,jj) = ozwk1(latsozp,k,i)
            enddo
          enddo
!ch <<

        else
          ll  = (theta(j1)-pl_lat(1))/5.0 + 1 
          con1 = (theta(j1)-pl_lat(ll)) / 5.0
          con2 = 1.0-con1
          !$acc loop vector collapse(2)
          do i = 1,pl_coeff
            do k = 1,levozp
              ozwk2(k,i,jj) = con2*ozwk1(ll,k,i) + con1*ozwk1(ll+1,k,i)
          enddo
          enddo
        endif
      
      enddo ! jj-loop

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
      !$acc parallel loop gang collapse(2) private(jj,j1,i,j,k) async(async_id)
      do jj = 1, jlistnum
        do k = 1, lev
          j1 = jlist1(jj)
          !$acc loop vector
          do i = 1, nxj(j1)
            if ( pp(i,k,jj) .le. pl_pres(1) )then
               do j = 1, pl_coeff
                 ozplout(i,k,j,jj) = ozwk2(1,j,jj)
               enddo
            elseif ( pp(i,k,jj) .ge. pl_pres(levozp) )then
               do j = 1, pl_coeff
                 ozplout(i,k,j,jj) = ozwk2(levozp,j,jj)
               enddo
            else
               do j = 1, pl_coeff
                 ozplout(i,k,j,jj) = 0.
               enddo
            endif
          enddo
        enddo
      enddo

!ch <<

      !$acc parallel loop gang collapse(3) private(jj,j1,kk,k,i,con1,con2) async(async_id)
      do jj = 1, jlistnum
        do j  = 1, pl_coeff
          do kk = 1, lev
            j1 = jlist1(jj)
            !$acc loop vector
            do k  = 1, levozp-1
!!ocl simd
              do i  = 1, nxj(j1)
                 con1 = (pp(i,kk,jj)-pl_pres(k)) / (pl_pres(k+1)-pl_pres(k)) 
                 if ((con1.gt.0.0) .and. (con1.le.1.0)) then !pp(kk) in pl(k,k+1)
                    con2=1.0-con1
                    ozplout(i,kk,j,jj) = con2*ozwk2(k,j,jj) + con1*ozwk2(k+1,j,jj)
                 endif
              enddo
            enddo
          enddo
        enddo
      enddo

!
!----------------------------------------------------------------------
   
        ldiag3d =.true.
!
!     do i=1,nx*lev*pl_coeff
!        ozp(i,1,1)=0.0
!     enddo
           !ozp=0.0

        pl_coeff2 = 2
!cmy--------------------------------------------------------------------
!for pl_coeff > 2
!cmy--------------------------------------------------------------------

      if (pl_coeff .gt. 2) then
        
        !$acc parallel loop collapse(2) private(jj,i) async(async_id)
        do jj = 1,jlistnum
          do i = 1, nx
            colo3(i,lev,jj) = 0.0
          enddo
        enddo
        
        !$acc parallel loop gang collapse(2) private(jj,j1,k) async(async_id)
        do jj = 1, jlistnum
          do k=2,lev
            j1 = jlist1(jj)
            !$acc loop vector 
            do i=1,nxj(j1)
               delp(i,k,jj)=pp(i,k,jj)-pp(i,k-1,jj)
            enddo
          enddo
        enddo
        
        !$acc parallel loop gang private(jj,j1,i,k) async(async_id)
        do jj = 1, jlistnum
          j1 = jlist1(jj)
          !$acc loop vector
          do i=1,nxj(j1)
             delp(i,1,jj)=pp(i,1,jj)
             colo3(i,1,jj)=o3l(i,1,jj) * delp(i,1,jj) * gravi
            !$acc loop seq
            do k=2,lev
               colo3(i,k,jj) = colo3(i,k-1,jj) + o3l(i,k,jj) * delp(i,k,jj) * gravi
            enddo
          enddo
        enddo
      endif ! for pl_coeff > 2
      
!cmy--------------------------------------------------------------------
!  K - loop start
!cmy--------------------------------------------------------------------

      if (pl_coeff2 .eq. 2) then
        !$acc parallel loop gang collapse(2) private(jj,j1,k,i) async(async_id)  
        do jj = 1, jlistnum
          do k=1,lev
            j1 = jlist1(jj)
            !$acc loop vector
            do i=1,nxj(j1)
              ozb(i,k,jj)    = o3l(i,k,jj)           ! NO FilliNG
              ozo(i,k,jj)  = (ozb(i,k,jj) + ozplout(i,k,1,jj)*dt) / (1.0 + ozplout(i,k,2,jj)*dt)
            enddo
            if (ldiag3d) then     !     Ozone change diagnostics
              !$acc loop vector
              do i=1,nxj(j1)
                ozp(i,k,1,jj) = ozp(i,k,1,jj) + ozplout(i,k,1,jj)*DT
                ozp(i,k,2,jj) = ozp(i,k,2,jj) + (ozo(i,k,jj) - ozb(i,k,jj))
              enddo
            endif ! for ldiag3d
          enddo
        enddo
      endif ! for pl_coeff2=2
!cmy--------------------------------------------------------------------
      if (pl_coeff2 .eq. 4) then 
        !$acc parallel loop gang collapse(2) private(jj,j1,k,i,temp) async(async_id)
        do jj = 1, jlistnum
          do k=1,lev
            j1 = jlist1(jj)
            !$acc loop vector
            do i=1,nxj(j1)
              ozb(i,k,jj)   = o3l(i,k,jj)            ! NO FilliNG
              temp     = ozplout(i,k,1,jj) + ozplout(i,k,3,jj)*tt(i,k,jj) + ozplout(i,k,4,jj)*colo3(i,k,jj)
              ozo(i,k,jj) = (OZB(i,k,jj)  + temp*dt) / (1.0 + ozplout(i,k,2,jj)*dt)
            enddo
            if (ldiag3d) then     !     Ozone change diagnostics
              !$acc loop vector
              do i=1,nxj(j1)
                OZP(i,k,1,jj) = OZP(i,k,1,jj) + ozplout(i,k,1,jj)*DT
                OZP(i,k,2,jj) = OZP(i,k,2,jj) + (OZO(i,k,jj)-OZB(i,k,jj))
                OZP(i,k,3,jj) = OZP(i,k,3,jj) + ozplout(i,k,3,jj)*tt(i,k,jj)*DT
                OZP(i,k,4,jj) = OZP(i,k,4,jj) + ozplout(i,k,4,jj)*colo3(i,k+1,jj)*DT
              enddo
            endif !for ldiag3d
          enddo
        enddo
      endif !for pl_coeff2=4

      
!cmy------------------------------------------------------------------
!     do i=1,nx*lev
!      o3l(i,1)=ozo(i,1)
!     enddo
      !$acc parallel loop gang private(jj,j1,k,i) async(async_id)
      do jj = 1, jlistnum
        j1 = jlist1(jj)
        !$acc loop vector
        do k=1,lev
          do i=1,nxj(j1)
            o3l(i,k,jj)=ozo(i,k,jj)
          enddo
        enddo
      enddo
!cmy------------------------------------------------------------------
!     if (myrank .eq. 0) print *,'*** ozphys.f end !!' 
!     if (myrank .eq. 0) print *,'*** o3l ***'
!     call qmax2d(o3l,1,1,nx,lev)
      !$acc exit data delete(jlist1, pl_lat, ozwk1, pl_pres, ozwk2, ozplout,&
      !$acc&                 colo3, delp, ozb, ozo, ozp)
      !$acc wait(async_id)

      RETURN
      END


