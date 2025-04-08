      subroutine rozphys_gpu(nxj, nx, lev, dt, iter, theta, julian, &
                             o3l, tt, pp, ps, myrank)
         use physcons, only: grav => con_g
         use ozne_def
         use const, only: RTYPE
         use param, only: my, my_max
         use index, only: jlist1, jlistnum
         use openacc
         use cudafor

         implicit none

         integer nxj(my), nx, lev, julian
         real(kind=kind_phys) ozwk1(latsozp, levozp, pl_coeff)
         real(kind=kind_phys) ozwk2(levozp, pl_coeff, my_max)
         real(kind=kind_phys) ozplout(nx, lev, pl_coeff, my_max)
         real(kind=kind_phys) amin(pl_coeff)
         real(kind=kind_phys) amax(pl_coeff)
         integer :: iter, myrank, lozc, n, n1, ll, kk, right, left, mid
         real :: theta(my), con1, con2, daynum

!
         real, parameter :: gravi = 1.0/grav
         integer pl_coeff2, kmax(pl_coeff), kmin(pl_coeff), myim(my_max)

         real(kind=kind_phys) ps(nx, my_max), &
            pp(nx, lev, my_max), ozp(nx, lev, pl_coeff, my_max)
         real(kind=RTYPE) o3l(nx, lev, my_max), tt(nx, lev, my_max)
         real(kind=kind_phys) dt
!
         integer k, i, j, jj, j1
         logical ldiag3d
         real(kind=kind_phys) ozb, colo3(nx, lev, my_max), &
            ozo, delp, temp
         integer async_id, istat
         real :: time2, time1, ppr, ozbr, ozor
         integer(kind=cuda_stream_kind) stream

         async_id = 1
         stream = acc_get_cuda_stream(async_id)
!----------------------------------------------------------------------
         lozc = latsozp*levozp*pl_coeff
!----------------------------------------------------------------------
!  check julian locating on which month is
!  and linearly interpolate for time
!----------------------------------------------------------------------
         daynum = julian
         if (daynum .gt. 365) daynum = 365.
         if (daynum .le. pl_time(1)) daynum = julian + 365.

         do n = 1, timeoz
            con1 = (daynum - pl_time(n))/(pl_time(n + 1) - pl_time(n))
            if ((con1 .ge. 0.0) .and. (con1 .le. 1.0)) then
               kk = n
               n1 = n + 1
               con2 = 1.0 - con1
               if (n1 .gt. 12) n1 = n1 - 12
               exit
            end if
         end do

         !$acc enter data create(ozwk1, ozwk2, ozplout, colo3, ozp, myim) async(async_id)

         !$acc parallel loop gang vector collapse(3) private(i,j,k) async(async_id)
         do k = 1, pl_coeff
            do j = 1, levozp
               do i = 1, latsozp
                  ozwk1(i, j, k) = con2*ozplin(i, j, k, n) + con1*ozplin(i, j, k, n1)
               end do
            end do
         end do

         !$acc parallel loop private(j1,jj) async(async_id)
         do jj = 1, jlistnum
            j1 = jlist1(jj)
            myim(jj) = nxj(j1)
         end do

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
         !$acc parallel loop gang vector collapse(3) private(jj,i,k,j1,ll,con1,con2) async(async_id)
         do jj = 1, jlistnum
            do i = 1, pl_coeff
               do k = 1, levozp
                  j1 = jlist1(jj)
                  if (theta(j1) .le. pl_lat(1)) then
                     ozwk2(k, i, jj) = ozwk1(1, k, i)
                  elseif (theta(j1) .ge. pl_lat(latsozp)) then
                     ozwk2(k, i, jj) = ozwk1(latsozp, k, i)
                  else
                     ll = (theta(j1) - pl_lat(1))/5.0 + 1
                     con1 = (theta(j1) - pl_lat(ll))/5.0
                     con2 = 1.0 - con1
                     ozwk2(k, i, jj) = con2*ozwk1(ll, k, i) + con1*ozwk1(ll + 1, k, i)
                  end if
               end do
            end do
         end do ! jj-loop
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
         !$acc parallel loop gang collapse(2) private(jj,i,j,ppr,con1,con2,&
         !$acc&                   right,left,mid) async(async_id)
         do jj = 1, jlistnum
            do kk = 1, lev
               !$acc loop vector
               do i = 1, myim(jj)
                  ppr = pp(i, kk, jj)
                  if (ppr .le. pl_pres(1)) then
                     !$acc loop seq
                     do j = 1, pl_coeff
                        ozplout(i, kk, j, jj) = ozwk2(1, j, jj)
                     end do
                  elseif (ppr .ge. pl_pres(levozp)) then
                     !$acc loop seq
                     do j = 1, pl_coeff
                        ozplout(i, kk, j, jj) = ozwk2(levozp, j, jj)
                     end do
                  else
                     left = 1
                     right = levozp
                     do while (right - left > 1)
                        mid = (left + right)/2
                        con1 = ppr - pl_pres(mid)
                        if (con1 .lt. 0.) then
                           right = mid
                        else
                           left = mid
                        end if
                     end do
                     con1 = (ppr - pl_pres(left))/(pl_pres(right) - pl_pres(left))
                     con2 = 1.0 - con1
                     !$acc loop seq
                     do j = 1, pl_coeff
                        ozplout(i, kk, j, jj) = con2*ozwk2(left, j, jj) + con1*ozwk2(right, j, jj)
                     end do
                  end if
               end do
            end do
         end do

!ch <<

!
!----------------------------------------------------------------------

         ldiag3d = .true.
!
!     do i=1,nx*lev*pl_coeff
!        ozp(i,1,1)=0.0
!     enddo
      !!$acc host_data use_device(ozp)
         !istat = cudaMemsetAsync(ozp, 0.0, size(ozp), stream)
      !!$acc end host_data
         pl_coeff2 = 2
!cmy--------------------------------------------------------------------
!for pl_coeff > 2
!cmy--------------------------------------------------------------------

         if (pl_coeff .gt. 2) then
            !$acc parallel loop gang vector collapse(2) private(jj,i,k,delp) async(async_id)
            do jj = 1, jlistnum
               do i = 1, nx
                  if (i .le. myim(jj)) then
                     delp = pp(i, 1, jj)
                     colo3(i, 1, jj) = o3l(i, 1, jj)*delp*gravi
                     colo3(i, lev, jj) = 0.0
                     !$acc loop seq
                     do k = 2, lev
                        delp = pp(i, k, jj) - pp(i, k - 1, jj)
                        colo3(i, k, jj) = colo3(i, k - 1, jj) + o3l(i, k, jj)*delp*gravi
                     end do
                  end if
               end do
            end do
         end if ! for pl_coeff > 2
!cmy--------------------------------------------------------------------
!  K - loop start
!cmy--------------------------------------------------------------------

         !$acc parallel loop gang collapse(2) private(jj,k,i,temp,ozb,ozo) async(async_id)
         do jj = 1, jlistnum
            do k = 1, lev
               !$acc loop vector
               do i = 1, myim(jj)
                  if (pl_coeff2 .eq. 2) then
                     ozb = o3l(i, k, jj)           ! NO FilliNG
                     ozo = (ozb + ozplout(i, k, 1, jj)*dt)/(1.0 + ozplout(i, k, 2, jj)*dt)
                     o3l(i, k, jj) = ozo
                     !     Ozone change diagnostics
                     if (ldiag3d) then
                        ozp(i, k, 1, jj) = ozplout(i, k, 1, jj)*DT
                        ozp(i, k, 2, jj) = (ozo - ozb)
                     end if ! for ldiag3d
                  end if ! for pl_coeff2=2
!cmy--------------------------------------------------------------------

                  if (pl_coeff2 .eq. 4) then
                     ozb = o3l(i, k, jj)            ! NO FilliNG
                     temp = ozplout(i, k, 1, jj) + ozplout(i, k, 3, jj)*tt(i, k, jj) + ozplout(i, k, 4, jj)*colo3(i, k, jj)
                     ozo = (OZB + temp*dt)/(1.0 + ozplout(i, k, 2, jj)*dt)
                     o3l(i, k, jj) = ozo
                     if (ldiag3d) then     !     Ozone change diagnostics
                        OZP(i, k, 1, jj) = ozplout(i, k, 1, jj)*DT
                        OZP(i, k, 2, jj) = (OZO - OZB)
                        OZP(i, k, 3, jj) = ozplout(i, k, 3, jj)*tt(i, k, jj)*DT
                        OZP(i, k, 4, jj) = ozplout(i, k, 4, jj)*colo3(i, k + 1, jj)*DT
                     end if !for ldiag3d
                  end if !for pl_coeff2=4
               end do
            end do
         end do

!cmy------------------------------------------------------------------
!     do i=1,nx*lev
!      o3l(i,1)=ozo(i,1)
!     enddo

!cmy------------------------------------------------------------------
!     if (myrank .eq. 0) print *,'*** ozphys.f end !!'
!     if (myrank .eq. 0) print *,'*** o3l ***'
!     call qmax2d(o3l,1,1,nx,lev)
         !$acc exit data delete(ozwk1,ozwk2, ozplout,&
         !$acc&                 colo3, ozp, myim) async(async_id)
         !$acc wait(async_id)

         RETURN
      END

