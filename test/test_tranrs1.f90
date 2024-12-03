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
   use const, only: RTYPE, poly, polyf, w => weight
   use index

   implicit none

   integer, parameter :: steps = 16
   integer, parameter :: num = 1
   ! real(kind=RTYPE), dimension(jtrun, my/2, jtmax) :: poly
   ! real(kind=RTYPE), dimension(my) :: w
   real(kind=RTYPE), dimension(nx, my_max) :: r
   real(kind=RTYPE), dimension(jtrun, jtmax, 2) :: s, s_gpu
   real(kind=RTYPE), dimension(nx + 2, my_max) :: cc, gwk1
   real, dimension(jtrun, jtmax, 2):: wss
   real, dimension(my, jtmax, 2) :: wcc
   real(kind=RTYPE), dimension(jtmax, my_max*nsizey, 2) :: wcc_fk
   real, dimension((jtrun + nsizey)*my/2*jtmax) :: fj_wp
   integer :: i, j, k, l, m
   integer :: async_id

   async_id = 1

   call random_seed()
   ! call random_number(poly)
   ! call random_number(w)
   call random_number(r)
   s = 0.
   s_gpu = 0.

   do i = 1, steps
      call tranrs1(jtrun, jtmax, nx, my, my_max, poly, w, r, s, nsizey)
   end do
   !$acc enter data create(wss, wcc, wcc_fk, fj_wp) async(async_id)
   !$acc enter data copyin(polyf, w, r, s_gpu, cc, gwk1) async(async_id)
   !$acc enter data copyin(jlist1, jlist2, nxdef, mtrundef, nlist, &
   !$acc& mlist, tcolt_jlist, poly_mlist) async(async_id)
   do i = 1, steps
      call tranrs1_gpu(jtrun, jtmax, nx, my, my_max, polyf, w, r, s_gpu, nsizey, cc, gwk1)
   end do
   !$acc exit data delete(jlist1, jlist2, nxdef, mtrundef, nlist, mlist, &
   !$acc& tcolt_jlist, poly_mlist, cc, gwk1) async(async_id)
   !$acc exit data copyout(polyf, w, r, s_gpu) async(async_id)
   !$acc exit data delete(wss, wcc, wcc_fk, fj_wp) async(async_id)
   !$acc wait(async_id)

#ifdef SP
   call assert_allclose(s_gpu, size(s_gpu), s, size(s), 1e-4_4, 1e-4_4, "Array s")
#else
   call assert_allclose(s_gpu, size(s_gpu), s, size(s), 1e-10, 1e-10, "Array s")
#endif

end subroutine tranrs1_unit
