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
   use const, only: RTYPE, poly, polyf
   use index

   implicit none

   integer, parameter :: steps = 10
   integer, parameter :: num = 1
   ! real(kind=RTYPE), dimension(jtrun, my/2, jtmax) :: poly
   real(kind=RTYPE), dimension(nx + 2, levp, num, my_max) :: cc, cc_gpu, gwk1
   real, dimension(levp, 2, num, jtrun, jtmax) :: s, wss
   real, dimension(levp, 2, num, my, jtmax) :: tcc
   real(kind=RTYPE), dimension(levp, 2, num, jtmax, my_max*nsizey) :: wcc_fk
   integer :: i, async_id

   async_id = 1

   call random_seed()
   ! call random_number(poly)
   call random_number(wss)
   cc = 0.
   cc_gpu = 0.

   do i = 1, steps
      call transr(jtrun, jtmax, nx, my, my_max, levp, poly, s, cc, num, nsizey)
   end do
   !$acc enter data create(wss, tcc, wcc_fk) async(async_id)
   !$acc enter data copyin(polyf, s, cc_gpu, jlist2, jlist1, mtrundef, mlist, nlist, nxdef, tcolt_jlist, poly_mlist, gwk1) async(async_id)
   do i = 1, steps
      ! call transr_gpu(jtrun, jtmax, nx, my, my_max, levp, poly, s, cc_gpu, num, nsizey, gwk1)
      call transr_gpu_cuda_graph(jtrun, jtmax, nx, my, my_max, levp, polyf, s, &
                                 cc_gpu, num, nsizey, &
                                 gwk1, wss, tcc, wcc_fk)

   end do
   !$acc exit data copyout(polyf, wss, cc_gpu, jlist2, jlist1, mtrundef, mlist, nlist, nxdef,  tcolt_jlist, poly_mlist, gwk1) async(async_id)
   !$acc exit data delete(wss, tcc, wcc_fk) async(async_id)
   !$acc wait(async_id)

#ifdef SP
   call assert_rmse(cc_gpu, size(cc_gpu), cc, size(cc), 1e-4_4, "Array cc")
#else
   call assert_allclose(cc_gpu, size(cc_gpu), cc, size(cc), 1e-10, 1e-10, "Array cc")
#endif

end subroutine
