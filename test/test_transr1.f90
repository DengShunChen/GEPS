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
   use const, only: RTYPE, poly, polyf
   use index

   implicit none

   integer, parameter :: steps = 16
   integer :: i, async_id
   ! real(kind=RTYPE), dimension(jtrun, my/2, jtmax) :: poly
   real(kind=RTYPE), dimension(jtrun, jtmax, 2) :: s
   real(kind=RTYPE), dimension(nxp, my_max) :: r, r_gpu
   real(kind=RTYPE), dimension(nx + 2, my_max) :: cc, gwk1
   real, dimension(jtrun, jtmax, 2):: wss
   real, dimension(my, jtmax, 2):: tcc
   real(kind=RTYPE), dimension(my_max*nsizey, jtmax, 2):: wcc_fk

   async_id = 1

   call random_seed()
   ! call random_number(poly)
   call random_number(s)
   r = 0.
   r_gpu = 0.

   do i = 1, steps
      call transr1(jtrun, jtmax, nx, my, my_max, poly, s, r, nsizey)
   end do
   !$acc enter data create(wss, tcc, wcc_fk) async(async_id)
   !$acc enter data copyin(poly, polyf, s, r_gpu, cc, gwk1) async(async_id)
   !$acc enter data copyin(jlist2, jlist1, mtrundef, mlist, nlist, nxjlen, &
   !$acc& nxjstart,nxjend, nxdef, mlist, tcolt_jlist, poly_mlist) async(async_id)
   do i = 1, steps
      call transr1_gpu(jtrun, jtmax, nx, my, my_max, polyf, s, r_gpu, nsizey, &
                       cc, gwk1)
      ! call transr1_gpu_cuda_graph(jtrun, jtmax, nx, my, my_max, polyf, &
      !      s, r_gpu, nsizey, &
      !      cc, gwk1, wss, tcc, wcc_fk)
   end do
   !$acc exit data delete(jlist2, jlist1, mtrundef, nlist, nxjlen, nxjstart, nxjend, nxdef, mlist, tcolt_jlist, poly_mlist) async(async_id)
   !$acc exit data copyout(polyf, s, r_gpu, cc, gwk1) async(async_id)
   !$acc exit data delete(wss, tcc, wcc_fk) async(async_id)
   !$acc wait(async_id)

#ifdef SP
   call assert_rmse(r_gpu, size(r_gpu), r, size(r), 1e-4_4, "Array r")
#else
   call assert_allclose(r_gpu, size(r_gpu), r, size(r), 1e-10, 1e-10, "Array r")
#endif

end subroutine
