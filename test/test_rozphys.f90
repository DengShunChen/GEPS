!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_rozphys
   implicit none
   call mpe_init
   call cons
   call getrdy
   call rozphys_unit
   call mpe_finalize

end program

subroutine assert_real(actual, n_actual, desired, n_desired, rtol, atol, err_msg)
   use const, only: RTYPE

   implicit none

   integer, intent(in) :: n_actual, n_desired
   real(kind=RTYPE), dimension(n_actual), intent(in) :: actual
   real(kind=RTYPE), dimension(n_desired), intent(in) :: desired
   real(kind=RTYPE), intent(in) :: rtol, atol
   character(len=*), intent(in), optional :: err_msg
   real(kind=RTYPE), parameter :: eps = 1e-15
   real(kind=RTYPE) :: rel_diff, abs_diff
   logical :: equal
   integer :: i

   if (n_actual .ne. n_desired) then
      print *, "Arrays have different lengths!"
      if (present(err_msg)) print *, err_msg
      call exit(1)
   end if

   equal = .true.
   rel_diff = 0.0
   abs_diff = 0.0

   rel_diff = maxval(abs((actual - desired)/(desired + eps)))
   abs_diff = maxval(abs(actual - desired))
   if (abs_diff > atol) equal = .false.
   if (rel_diff > rtol) equal = .false.

   if (.not. equal) then
      print *, "Arrays are not close within tolerance rtol =", rtol, "atol =", atol
      print *, "Max relative difference = ", rel_diff
      print *, "Max absolute difference = ", abs_diff
      if (present(err_msg)) print *, err_msg
      call exit(1)
   end if

end subroutine assert_real

subroutine assert_integer(actual, n_actual, desired, n_desired, err_msg)
   use const, only: RTYPE

   implicit none

   integer, intent(in) :: n_actual, n_desired
   integer, dimension(n_actual), intent(in) :: actual
   integer, dimension(n_desired), intent(in) :: desired
   character(len=*), intent(in), optional :: err_msg
   integer :: abs_diff
   logical :: equal
   integer :: i

   if (n_actual .ne. n_desired) then
      print *, "Arrays have different lengths!"
      if (present(err_msg)) print *, err_msg
      call exit(1)
   end if

   equal = .true.
   abs_diff = maxval(abs(actual - desired))
   if (abs_diff .ne. 0) equal = .false.

   if (.not. equal) then
      print *, "Arrays are not close to 0."
      print *, "Max absolute difference = ", abs_diff
      if (present(err_msg)) print *, err_msg
      call exit(1)
   end if

end subroutine assert_integer

subroutine rozphys_unit
   use rank, only: myrank
   use index, only: jlist1, jlistnum, nxp, nxjp
   use param, only: my_max, lev, my
   use const, only: julian, RTYPE
   use phygrid, only: xlat
   use machine, only: kind_phys
   use ozne_def
   !use nvtx

   implicit none

   real(kind=RTYPE) :: o3l(nxp, lev, my_max), tt(nxp, lev, my_max)
   real(kind=RTYPE) :: o3l_gpu(nxp, lev, my_max)
   real(kind=kind_phys) :: plt(nxp, lev, my_max)
   real(kind=kind_phys) :: ps(nxp, my_max)
   real(kind=kind_phys) :: dta, mmax
   real :: iter_r, pi, diff, rel_diff, rtol
   integer :: iter

   integer :: i, j, jj, async_id, n, k
   integer, dimension(34) :: seed

   real :: time1, time2, ct, gt

   async_id = 1
   seed = (/10003, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,&
           &25, 26, 27, 28, 29, 30, 31, 32, 33/)
   !call random_seed()
   call random_seed(put=seed)
   call random_number(o3l)
   call random_number(tt)
   call random_number(plt)
   call random_number(ps)
   call random_number(dta)
   call random_number(iter_r)
   iter = floor(iter_r*6)
   o3l_gpu = o3l

   ct = 0.
   do i = 1, 16
      call cpu_time(time1)
      !if (i .ne. 1) call nvtxStartRange("CPU compute")
      do jj = 1, jlistnum
         j = jlist1(jj)
         call rozphys(nxjp(j), nxp, lev, dta, iter, xlat(j), julian, o3l(1, 1, jj), &
                      tt(1, 1, jj), plt(1, 1, jj), ps(1, jj), myrank)
      end do
      !if (i .ne. 1) call nvtxEndRange
      call cpu_time(time2)
      if (i .ne. 1) ct = ct + time2 - time1
   end do

   if (.true.) then

      !$acc enter data copyin(ozplin,jlist1, pl_lat, pl_pres)
      !$acc enter data copyin(nxjp, nxp, lev, dta, iter, xlat, julian,      &
      !$acc&            tt, plt, ps, myrank)
      !$acc enter data copyin(o3l_gpu)
      gt = 0.
      do i = 1, 16
         call cpu_time(time1)
         !if (i .ne. 1) call nvtxStartRange("GPU compute")
         call rozphys_gpu(nxjp, nxp, lev, dta, iter, xlat, julian, o3l_gpu, &
                          tt, plt, ps, myrank)
         !if (i .ne. 1) call nvtxEndRange
         call cpu_time(time2)
         if (i .ne. 1) gt = gt + time2 - time1
      end do
      !$acc exit data copyout(o3l_gpu)
      !$acc exit data delete(nxjp, nxp, lev, dta, iter, xlat, julian,      &
      !$acc&            tt, plt, ps, myrank)
      !$acc exit data delete(ozplin,jlist1, pl_lat, pl_pres)

      diff = 0.
      do jj = 1, jlistnum
         j = jlist1(jj)
         do k = 1, lev
            do i = 1, nxjp(j)
               rel_diff = (abs(o3l(i, k, jj) - o3l_gpu(i, k, jj))/o3l(i, k, jj) + 1e-15)
               if (diff .lt. rel_diff) then
                  diff = rel_diff
               end if
            end do
         end do
      end do
      rtol = 1e-8
      if (diff .gt. rtol) then
         print *, "Arrays are not close within tolerance rtol =", rtol
         print *, "Max relative difference = ", diff
         call exit(1)
      end if
   end if
   write (*, *) 'Timing for CPU & GPU: ', ct, gt
   write (*, *) 'speedup ratio: ', ct/gt

end subroutine rozphys_unit
!mmax = 0.
!do j = 1, my_max
!do k = 1, lev
!do i = 1, nxp
!  if (mmax .lt. abs(o3l(i,k,j) - o3l_gpu(i,k,j))) then
!    mmax = abs(o3l(i,k,j) - o3l_gpu(i,k,j))
!    write(*,*) 'i,k,j,o3l,o3l_gpu:', i,k,j,o3l(i,k,j),o3l_gpu(i,k,j)
!  endif
!enddo
!enddo
!enddo

!do j = 1, my_max
!do k = 1, lev
!do i = 1, nxp
!  if (abs(o3l(i,k,j) - o3l_gpu(i,k,j)) .gt. 1.e-8) then
!    write(*,*) 'i,k,j,o3l,o3l_gpu:', i,k,j,o3l(i,k,j),o3l_gpu(i,k,j)
!  endif
!enddo
!enddo
!enddo
