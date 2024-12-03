!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_prexp_hybrid_cwb
   use param
   use const

   implicit none

   call mpe_init
   call cons
   call prexp_hybrid_cwb_unit
   call mpe_finalize

end program

subroutine prexp_hybrid_cwb_unit
   !$acc routine(prexp_hybrid_cwb_gpu) vector
   use param
   use const, only: RTYPE, ptop, sigma
   use index

   implicit none

   integer, parameter :: steps = 16
   integer :: i, jj, j, nxj, async_id, seed_size
   integer, allocatable :: seed(:)
   ! real(kind=RTYPE) sigma(lev + 1, 2)
   real(kind=RTYPE) ptm(nxp, my_max)
   real(kind=RTYPE) pk(nxp, lev, my_max), pk_gpu(nxp, lev, my_max)
   real(kind=RTYPE) pk2(nxp, lev, my_max), pk2_gpu(nxp, lev, my_max)
   real(kind=RTYPE) plt(nxp, lev, my_max), plt_gpu(nxp, lev, my_max)

   async_id = 1

   call random_seed(size=seed_size)
   allocate (seed(seed_size))
   seed = 123
   call random_seed(put=seed)
   !! call random_number(sigma)
   call random_number(ptm)
   ptm = ptm*ptop
   pk = 0.
   pk_gpu = 0.
   pk2 = 0.
   pk2_gpu = 0.
   plt = 0.
   plt_gpu = 0.

   do i = 1, steps
      do jj = 1, jlistnum
         j = jlist1(jj)
         call prexp_hybrid_cwb(nxjp(j), nxp, lev, ptop, sigma, ptm(1, jj), pk(1, 1, jj), pk2(1, 1, jj), plt(1, 1, jj))
      end do
   end do

   !$acc enter data copyin(sigma, ptm, pk_gpu, pk2_gpu, plt_gpu,jlistnum, lev, jlist1, nxjp, nxdef_2d) async(async_id)
   !$acc wait(async_id)
   do i = 1, steps
      ! !$acc parallel loop gang async(async_id)
      ! do jj = 1, jlistnum
      !    j = jlist1(jj)
      !  call prexp_hybrid_cwb_gpu(nxjp(j), nxp, lev, ptop, sigma, ptm(1, jj), pk_gpu(1, 1, jj), pk2_gpu(1, 1, jj), plt_gpu(1, 1, jj))
      ! end do
      ! ------------------------------------------------------------
      call prexp_hybrid_cwb_gpu_refactor(nxjp, nxp, lev, ptop, sigma, ptm, pk_gpu, pk2_gpu, plt_gpu)
   end do
   !$acc wait(async_id)
   !$acc exit data copyout(pk_gpu, pk2_gpu, plt_gpu) delete(sigma, ptm, jlist1, nxjp, nxdef_2d) async(async_id)
   !$acc wait(async_id)

#ifdef SP
   call assert_rmse(pk_gpu, size(pk_gpu), pk, size(pk), 1e-4_4, "Array pk")
   call assert_rmse(pk2_gpu, size(pk2_gpu), pk2, size(pk2), 1e-4_4, "Array pk2")
   call assert_rmse(plt_gpu, size(plt_gpu), plt, size(plt), 1e-4_4, "Array plt")
#else
   call assert_allclose(pk_gpu, size(pk_gpu), pk, size(pk), 1e-8, 1e-8, "Array pk")
   call assert_allclose(pk2_gpu, size(pk2_gpu), pk2, size(pk2), 1e-8, 1e-8, "Array pk2")
   call assert_allclose(plt_gpu, size(plt_gpu), plt, size(plt), 1e-8, 1e-8, "Array plt")
#endif

end subroutine
