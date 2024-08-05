!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_tranrs
   use param
   use const

   implicit none

   call mpe_init
   call cons
   call tranrs_unit
   call mpe_finalize

end program

subroutine tranrs_unit
   use param
   use const, only: RTYPE
   use index

   implicit none

   integer, parameter :: steps = 10
   integer, parameter :: num = 1
   real(kind=RTYPE), dimension(jtrun, my/2, jtmax) :: poly
   real(kind=RTYPE), dimension(my) :: w
   real(kind=RTYPE), dimension(nx + 2, lev, num, my_max) :: cc, cc_buffer, gwk1
   real(kind=RTYPE), dimension(lev, 2, num, jtrun, jtmax) :: wss, wss_gpu
   integer :: i, j, k, l, m
   integer :: async_id

   async_id = 1

   call random_seed()
   call random_number(poly)
   call random_number(w)
   call random_number(cc)
   wss = 0.
   wss_gpu = 0.

   ! tranrs will modify cc, so a copy cc_buffer is made
   do i = 1, steps
      cc_buffer = cc
      call tranrs(jtrun, jtmax, nx, my, my_max, lev, poly, w, cc_buffer, wss, num, nsizey)
   end do

   !$acc enter data copyin(poly, w, cc_buffer, nlist, jlist2, jlist1, nxdef, cc) create(wss_gpu, gwk1) async(async_id)
   do i = 1, steps
      !$acc kernels async(async_id)
      cc_buffer = cc
      !$acc end kernels
      call tranrs_gpu(jtrun, jtmax, nx, my, my_max, lev, poly, w, cc_buffer, wss_gpu, num, nsizey, gwk1)
   end do
   !$acc exit data delete(poly, w, cc_buffer, nlist, jlist2, jlist1, nxdef, gwk1, cc) copyout(wss_gpu) async(async_id)
   !$acc wait(async_id)

   call assert_allclose(wss_gpu, size(wss_gpu), wss, size(wss), 1e-10, 1e-10, "Array wss")

end subroutine
