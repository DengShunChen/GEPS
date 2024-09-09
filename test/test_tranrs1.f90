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
   real(kind=RTYPE), dimension(jtrun, my/2, jtmax) :: poly
   real(kind=RTYPE), dimension(my) :: w
   real(kind=RTYPE), dimension(nx, my_max) :: r
   real(kind=RTYPE), dimension(jtrun, jtmax, 2) :: s, s_gpu
   real(kind=RTYPE), dimension(nx + 2, my_max) :: cc, gwk1
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

   !$acc enter data copyin(poly, w, r, s_gpu, cc, gwk1) async(async_id)
   !$acc enter data copyin(jlist1, nxdef, mtrundef, nlist) async(async_id)
   do i = 1, steps
      call tranrs1_gpu(jtrun, jtmax, nx, my, my_max, poly, w, r, s_gpu, nsizey, cc, gwk1)
   end do
   !$acc exit data delete(jlist1, nxdef, mtrundef, nlist, cc, gwk1) async(async_id)
   !$acc exit data copyout(poly, w, r, s_gpu) async(async_id)
   !$acc wait(async_id)

   call assert_allclose(s_gpu, size(s_gpu), s, size(s), 1e-10, 1e-10, "Array s")

end subroutine tranrs1_unit
