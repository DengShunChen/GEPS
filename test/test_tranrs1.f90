!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_tranrs1
   use param
   use const

   implicit none

   call mpe_init
   call cons
   call tranrs1_unit
   call mpe_finalize

end program

subroutine tranrs1_unit
   use param
   use const, only: RTYPE
   use index

   implicit none

   integer, parameter :: steps = 16
   integer, parameter :: num = 1
   real(kind=RTYPE) poly(jtrun, my/2, jtmax)
   real(kind=RTYPE) w(my)
   real(kind=RTYPE) r(nx, my_max)
   real(kind=RTYPE) s(jtrun, jtmax, 2)
   real(kind=RTYPE) s_gpu(jtrun, jtmax, 2)
   integer :: i, j, k, l, m
   integer :: async_id

   async_id = 1

   call random_seed()
   call random_number(poly)
   call random_number(w)
   call random_number(r)
   s = 0.
   s_gpu = 0.

   do i = 1, steps
      call tranrs1(jtrun, jtmax, nx, my, my_max, poly, w, r, s, nsizey)
   end do

   !$acc enter data copyin(poly, w, r, s_gpu) async(async_id)
   !$acc enter data copyin(jlist1, nxdef, mtrundef, nlist) async(async_id)
   do i = 1, steps
      call tranrs1_gpu(jtrun, jtmax, nx, my, my_max, poly, w, r, s_gpu, nsizey)
   end do
   !$acc exit data delete(jlist1, nxdef, mtrundef, nlist) async(async_id)
   !$acc exit data copyout(poly, w, r, s_gpu) async(async_id)
   !$acc wait(async_id)

   if (all(abs(s - s_gpu) <= 1e-10)) then
      print *, "test_tranrs1 passed."
   else
      print *, "test_tranrs1 failed."
      call exit(1)
   end if

end subroutine
