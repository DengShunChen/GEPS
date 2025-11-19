      subroutine ozphys_2015_gpu(ix, nxjp, lev, dt, xlat, julian, ozi, &
                                ozo, tin, prsl, dsigma, pst, me)
!                              prsl, ozplout, pl_coeff, delp, ldiag3d,&
!                              ozp,me)
!
!     this code assumes that both prsl and po3 are from bottom to top
!     as are all other variables
!     This code is specifically for NRL parameterization and
!     climatological T and O3 are in location 5 and 6 of ozplout array
! June 2015 - Shrinivas Moorthi
!
         use machine, only: kind_phys
         use physcons, only: grav => con_g
         use ozne_def
         use const, only: RTYPE
         use param, only: my, my_max
         use index, only: jlist1, jlistnum, nxdef_2d
         implicit none
!
         real, parameter :: gravi = 1.0/grav
         integer ix, lev, me, julian, myim(my_max), nxjp(my), nxj(my_max)
         real theta(my_max), dt, xlat(my)
         real(kind=kind_phys) po3(levozp), prsl(ix, lev, my_max), delp
         real(kind=RTYPE) tin(ix, lev, my_max), ozi(ix, lev, my_max), ozo(ix, lev, my_max), &
            dsigma(lev, 2), pst(ix, my_max)
!                          ozp(ix,lev,4),  dt
!
         integer k, kk, kmax, kmin, &
            l, ll, i, j, n, n1, lozc, daynum, jj
         real con1, con2
         logical ldiag3d
         real(kind=kind_phys) ozwk11, ozwk12
         real(kind=kind_phys) ozwk2(levozp, pl_coeff, my_max)
         real(kind=kind_phys) ozplout(ix, lev, pl_coeff, my_max)
         real(kind=kind_phys) pmax, pmin, tem, temp
         real(kind=kind_phys) wk1, wk2, wk3, &
            prod(ix, pl_coeff, lev, my_max), &
            ozib, colo3, coloz
            
         ! for GPU porting
         integer :: async_id = 1
         ! GPU register
         integer :: mid, left, right
         real :: ppr, colo3r, colozr, tmp1, tmp2
!----------------------------------------------------------------------
         lozc = latsozp*levozp*pl_coeff
!----------------------------------------------------------------------
!  check julian locating on which month is
!  and linearly interpolate for time
!----------------------------------------------------------------------
         daynum = julian
         if (daynum .gt. 365) daynum = 365.
         if (daynum .le. pl_time(1)) daynum = julian + 365.

         do n = 1, timeoz ! 12
            con1 = (daynum - pl_time(n))/(pl_time(n + 1) - pl_time(n))
            n1 = n + 1
            if (n1 .gt. 12) n1 = n1 - 12
            if ((con1 .ge. 0.0) .and. (con1 .le. 1.0)) then
               con2 = 1.0 - con1
               exit
            end if
         end do
         !$acc data create(myim, theta, nxj, ozwk2, ozplout, &
         !$acc&     prod) async(async_id)
         
         
         !$acc parallel loop async(async_id) private(j)
         do jj = 1, jlistnum
            j = jlist1(jj)
            myim(jj) = nxjp(j)
            theta(jj) = xlat(j)
            nxj(jj) = nxdef_2d(j)
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
         !$acc parallel loop async(async_id) collapse(3) private(ll, tmp1, tmp2) &
         !$acc&         firstprivate(con1, con2, n, n1)
         do jj = 1, jlistnum
            do i = 1, pl_coeff
               do k = 1, levozp
                  if (theta(jj) .le. pl_lat(1)) then
                     ozwk2(k, i, jj) = con2*ozplin(1, k, i, n) + con1*ozplin(1, k, i, n1)
                  elseif (theta(jj) .ge. pl_lat(latsozp)) then
                     ozwk2(k, i, jj) = con2*ozplin(latsozp, k, i, n) + &
                                       con1*ozplin(latsozp, k, i, n1)
   !ch <<
                  else
                     ll = (theta(jj) - pl_lat(1))/5.0 + 1
                     ozwk11 = con2*ozplin(ll, k, i, n) + con1*ozplin(ll, k, i, n1)
                     ozwk12 = con2*ozplin(ll + 1, k, i, n) + con1*ozplin(ll + 1, k, i, n1)
                     tmp1 = (theta(jj) - pl_lat(ll))/5.0
                     tmp2 = 1.0 - tmp1
                     ozwk2(k, i, jj) = tmp2*ozwk11 + tmp1*ozwk12
                  end if
               end do
            end do
         end do

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
         !$acc parallel loop async(async_id) collapse(2) gang private(mid, &
         !$acc&         left, right, ppr)
         do jj = 1, jlistnum
            do k = 1, lev
               !$acc loop vector
               do i = 1, myim(jj)
                  ppr = prsl(i, k, jj)
                  if (ppr .le. pl_pres(1)) then
                     !$acc loop seq
                     do j = 1, pl_coeff
                        ozplout(i, k, j, jj) = ozwk2(1, j, jj)
                     end do
                  elseif (ppr .ge. pl_pres(levozp)) then
                     !$acc loop seq
                     do j = 1, pl_coeff
                        ozplout(i, k, j, jj) = ozwk2(levozp, j, jj)
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
                        ozplout(i, k, j, jj) = con2*ozwk2(left, j, jj) + con1*ozwk2(right, j, jj)
                     end do
                  end if
               end do
            end do
         end do
   !ch <<

   !
   !
   !----------------------------------------------------------------------

         ldiag3d = .true.

   !
   !
         !$acc parallel loop async(async_id) collapse(2) gang private(wk2, wk3, &
         !$acc&         temp, kmax, kmin, wk1)
         do jj = 1, jlistnum
            do l = lev, 1, -1
               pmin = 1.0e10
               pmax = -1.0e10
   !
               !$acc loop vector reduction(min:pmin) reduction(max:pmax)
               do i = 1, myim(jj)
   !          wk1(i) = log(prsl(i,l))
                  wk1 = prsl(i, l, jj)
                  pmin = min(wk1, pmin)
                  pmax = max(wk1, pmax)
                  !$acc loop seq
                  do j = 1, pl_coeff
                     prod(i, j, l, jj) = 0.0
                  end do
               end do
               kmax = 1
               kmin = 1
               !$acc loop seq
               do k = levozp - 1, 1, -1
                  if (pmin < pl_pres(k)) kmax = k
                  if (pmax < pl_pres(k)) kmin = k
               end do
   !
               !$acc loop vector
               do i = 1, myim(jj)
                  wk1 = prsl(i, l, jj)
                  !$acc loop seq
                  do k = kmin, kmax
   !          temp = 1.0 / (po3(k) - po3(k+1))
                  temp = 1.0/(pl_pres(k + 1) - pl_pres(k))
   !            if (wk1(i) < po3(k) .and. wk1(i) >= po3(k+1)) then
                     if (wk1 < pl_pres(k + 1) .and. wk1 >= pl_pres(k)) then
   !              wk2(i) = (wk1(i) - po3(k+1)) * temp
                        wk2 = (wk1 - pl_pres(k))*temp
                        wk3 = 1.0 - wk2
                        !$acc loop seq
                        do j = 1, pl_coeff
                           prod(i, j, l, jj) = wk2*ozplout(i, k + 1, j, jj) &
                                        + wk3*ozplout(i, k, j, jj)
                        end do
                     end if
                  end do
!            if (wk1(i) < po3(levozp)) then
                  if (wk1 >= pl_pres(levozp)) then
                     !$acc loop seq
                     do j = 1, pl_coeff
                        prod(i, j, l, jj) = ozplout(i, levozp, j, jj)
                     end do
                  end if
!            if (wk1(i) >= po3(1)) then
                  if (wk1 < pl_pres(1)) then
                     !$acc loop seq
                     do j = 1, pl_coeff
                        prod(i, j, l, jj) = ozplout(i, 1, j, jj)
                     end do
                  end if
               end do
            end do
         end do
            
   !
         !$acc parallel loop async(async_id) collapse(2) private(delp, colo3r, colozr, tem)
         do jj = 1, jlistnum
            do i = 1, ix
               if (i .le. myim(jj)) then
                  colo3r = 0.
                  colozr = 0.
                  !$acc loop seq
                  do l = lev, 1, -1
                     delp = 100.0*(dsigma(l, 1)*pst(i, jj) + dsigma(l, 2))  !  pa
                     colo3r = colo3r + ozi(i, l, jj)*delp*gravi
                     colozr = colozr + prod(i, 6, l, jj)*delp*gravi
                     colo3 = colo3r
                     coloz = colozr
   !       write(1000+me,*) ' colo3=',colo3(1,l),' coloz=',coloz(1,l)
   !    &,' l=',l
                     prod(i, 2, l, jj) = min(prod(i, 2, l, jj), 0.0)
                     tem = prod(i, 1, l, jj) - prod(i, 2, l, jj)*prod(i, 6, l, jj) &
                           + prod(i, 3, l, jj)*(tin(i, l, jj) - prod(i, 5, l, jj)) &
                           + prod(i, 4, l, jj)*(colo3 - coloz)

      !     if (me .eq. 0) print *,'ozphys_2015 tem=',tem,' prod=',prod(i,:)
      !    &,' ozib=',ozib(i),' l=',l,' tin=',tin(i,l),'colo3=',colo3(i,l+1)

                     ozo(i, l, jj) = (ozi(i, l, jj) + tem*dt)/(1.0 - prod(i, 2, l, jj)*dt)
                  end do
   !        if (ldiag3d) then     !     ozone change diagnostics
   !          do i=1,im
   !            ozp(i,l,1) = ozp(i,l,1) + (prod(i,1)-prod(i,2)*prod(i,6))*dt
   !            ozp(i,l,2) = ozp(i,l,2) + (ozo(i,l) - ozib(i))
   !            ozp(i,l,3) = ozp(i,l,3) + prod(i,3)*(tin(i,l)-prod(i,5))*dt
   !            ozp(i,l,4) = ozp(i,l,4) + prod(i,4)
   !     &                              * (colo3(i,l)-coloz(i,l))*dt
   !          enddo
   !        endif
               end if
            end do                                ! vertical loop
         end do
         !$acc end data
!
         return
      end
