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
   use const, only: RTYPE, poly, polyf, w => weight
   use index

   implicit none

   integer, parameter :: steps = 10
   integer, parameter :: num = 1
   ! real(kind=RTYPE), dimension(jtrun, my/2, jtmax) :: poly
   ! real(kind=RTYPE), dimension(my) :: w
   real(kind=RTYPE), dimension(nx + 2, lev, num, my_max) :: cc, cc_buffer, gwk1
   real(kind=RTYPE), dimension(lev, 2, num, jtrun, jtmax) :: s, s_gpu
   real(kind=RTYPE), dimension(lev, 2, num, jtmax, my_max*nsizey) :: wcc_fk
   real, dimension(lev, 2, num, jtrun, jtmax):: wss
   real, dimension(lev, 2, num, my, jtmax) :: wcc
   real, dimension((jtrun + nsizey)*my/2*jtmax) :: fj_wp
   integer :: i, j, k, l, m
   integer :: async_id

   async_id = 1

   call random_seed()
   ! call random_number(poly)
   ! call random_number(w)
   call random_number(cc)
   s = 0.
   s_gpu = 0.

   ! tranrs will modify cc, so a copy cc_buffer is made
   do i = 1, steps
      cc_buffer = cc
      call tranrs(jtrun, jtmax, nx, my, my_max, lev, poly, w, cc_buffer, s, num, nsizey)
   end do
   !$acc enter data create(wss, wcc, wcc_fk, fj_wp) async(async_id)
   !$acc enter data copyin(polyf, w, cc_buffer, nlist, jlist2, jlist1, nxdef, cc, &
   !$acc& mlist, tcolt_jlist, poly_mlist) create(s_gpu, gwk1) async(async_id)
   do i = 1, steps
      !$acc kernels async(async_id)
      cc_buffer = cc
      !$acc end kernels
      ! call tranrs_gpu(jtrun, jtmax, nx, my, my_max, lev, poly, w, cc_buffer, s_gpu, num, nsizey, gwk1)
      call tranrs_gpu_cuda_graph(jtrun, jtmax, nx, my, my_max, lev, polyf, w, &
                                 cc_buffer, s_gpu, num, nsizey, &
                                 gwk1, wss, wcc, wcc_fk, fj_wp)
   end do
   !$acc exit data delete(polyf, w, cc_buffer, nlist, jlist2, jlist1, nxdef, cc, gwk1, &
   !$acc& mlist, tcolt_jlist, poly_mlist) copyout(s_gpu) async(async_id)
   !$acc exit data delete(wss, wcc, wcc_fk, fj_wp) async(async_id)
   !$acc wait(async_id)

#ifdef SP
   call assert_rmse(s_gpu, size(s_gpu), s, size(s), 1e-4_4, "Array s")
#else
   call assert_allclose(s_gpu, size(s_gpu), s, size(s), 1e-10, 1e-10, "Array s")
#endif

end subroutine
