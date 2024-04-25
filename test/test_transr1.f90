!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_transr1
   use param
   use const

   implicit none

   call mpe_init
   call cons
   call transr1_unit
   call mpe_finalize

end program

subroutine transr1_unit
   use param
   use const, only: RTYPE
   use index

   implicit none

   integer, parameter :: steps = 16
   integer :: i
   real(kind=RTYPE) :: poly(jtrun, my/2, jtmax)
   real(kind=RTYPE) :: s(jtrun, jtmax, 2)
   real(kind=RTYPE) :: r(nxp, my_max)
   real(kind=RTYPE) :: r_gpu(nxp, my_max)
   integer :: async_id

   async_id = 1

   call random_seed()
   call random_number(poly)
   call random_number(s)
   r = 0.
   r_gpu = 0.

   do i = 1, steps
      call transr1(jtrun, jtmax, nx, my, my_max, poly, s, r, nsizey)
   end do
   !$acc enter data copyin(poly, s, r_gpu) async(async_id)
   !$acc enter data copyin(jlist2, jlist1, mtrundef, nlist, nxjlen, nxjstart, nxjend) async(async_id)
   do i = 1, steps
      call transr1_gpu(jtrun, jtmax, nx, my, my_max, poly, s, r_gpu, nsizey)
   end do
   !$acc exit data delete(jlist2, jlist1, mtrundef, nlist, nxjlen, nxjstart, nxjend) async(async_id)
   !$acc exit data copyout(poly, s, r_gpu) async(async_id)
   !$acc wait(async_id)

   if (all(abs(r - r_gpu) <= 1e-10)) then
      print *, "test_transr1 passed."
   else
      print *, "test_transr1 failed."
      call exit(1)
   end if

end subroutine
