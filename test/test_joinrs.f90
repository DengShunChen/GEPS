!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_joinrs
   use param
   use const

   implicit none

   call mpe_init
   call cons
   call joinrs_unit
   call mpe_finalize

end program

subroutine joinrs_unit
   use param
   use const, only: RTYPE
   use index

   implicit none

   integer, parameter :: steps = 10
   real(kind=RTYPE), dimension(nx + 2, levp, 1 + ncld, my_max) :: cc, cc_gpu
   real(kind=RTYPE), dimension(nxp, lev, my_max) :: r1
   real(kind=RTYPE), dimension(nxp, lev*ncld, my_max) :: r2
   real(kind=RTYPE) dummy
   integer :: i
   integer :: async_id

   integer j, k, n, jj, kk, nk, nxj

   async_id = 1

   call random_seed()
   call random_number(r1)
   call random_number(r2)

   cc = 0.
   cc_gpu = 0.

   do i = 1, steps
      call joinrs(cc, r1, r2, dummy, dummy, nx, my_max, lev, jlistnum, 2, ncld)
   end do
   !$acc enter data create(cc_gpu) copyin(r1, r2, jlist1, nxjlen, nxjlen_all) async(async_id)
   do i = 1, steps
      call joinrs_gpu(cc_gpu, r1, r2, dummy, dummy, nx, my_max, lev, jlistnum, 2, ncld)
   end do
   !$acc exit data delete(r1, r2, jlist1, nxjlen, nxjlen_all) copyout(cc_gpu) async(async_id)
   !$acc wait(async_id)

   ! call assert_realspace_close(cc_gpu, cc, nx+2, levp, ncld+1, &
   !      1.e-10, 1.e-10, "joinsr1_sr")
   if (all(abs(cc - cc_gpu) <= 1e-10)) then
      print *, "(all close) test_joinrs passed."
   else
      print *, "(all close) test_joinrs failed."
      call exit(1)
   end if

end subroutine
