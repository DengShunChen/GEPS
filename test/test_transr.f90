!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_tranuv
   use param
   use const

   implicit none

   call mpe_init
   call cons
   call transr_unit
   call mpe_finalize

end program

subroutine transr_unit
   use param
   use const, only: RTYPE
   use index

   implicit none

   integer, parameter :: steps = 10
   integer, parameter :: num = 1
   real(kind=RTYPE), dimension(jtrun, my/2, jtmax) :: poly
   real(kind=RTYPE), dimension(lev, 2, num, jtrun, jtmax) :: wss
   real(kind=RTYPE), dimension(nx + 2, lev, num, my_max) :: cc, cc_gpu, gwk1
   integer :: i, async_id

   async_id = 1

   call random_seed()
   call random_number(poly)
   call random_number(wss)
   cc = 0.
   cc_gpu = 0.

   do i = 1, steps
      call transr(jtrun, jtmax, nx, my, my_max, levp, poly, wss, cc, num, nsizey)
   end do
   !$acc enter data copyin(poly, wss, cc_gpu, jlist2, jlist1, mtrundef, nlist, nxdef, gwk1) async(async_id)
   do i = 1, steps
      call transr_gpu(jtrun, jtmax, nx, my, my_max, levp, poly, wss, cc_gpu, num, nsizey, gwk1)
   end do
   !$acc exit data copyout(poly, wss, cc_gpu, jlist2, jlist1, mtrundef, nlist, nxdef, gwk1) async(async_id)
   !$acc wait(async_id)

   call assert_allclose(cc_gpu, size(cc_gpu), cc, size(cc), 1e-10, 1e-10, "Array cc")

end subroutine
