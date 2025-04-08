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
   use rank, only: myrank
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
   real plt(nxp, lev, my_max), plt_gpu(nxp, lev, my_max)

   real(kind=RTYPE) err(3)
   async_id = 1

   call random_seed(size=seed_size)
   allocate (seed(seed_size))
   seed = 123
   call random_seed(put=seed)
   !! call random_number(sigma)
   call random_number(ptm)
   ptm = ptm*1000.+300.
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

#ifdef VERBOSE
   if (myrank .eq. 0) then
      write (*, '(1X, A8, 1X, A15, 4A15)') &
         "name", "max|Err|/max|A|", "L2-Err", "RMSE", "max|A|", "max|B|"
   end if
   call varErr(err, pk_gpu, pk, lev, 1, 'pk')
   call varErr(err, pk2_gpu, pk2, lev, 1, 'pk2')
#endif

#ifdef SP
   call assert_realspace_close(pk_gpu, pk, nxp, lev, 1, 1e-4_4, 1e-4_4, "Array pk")
   call assert_realspace_close(pk2_gpu, pk2, nxp, lev, 1, 1e-4_4, 1e-4_4, "Array pk2")
   ! use r4 to compute, but store in r8 array.
   ! call assert_allclose_r8(plt_gpu, size(plt_gpu), plt, size(plt), 1.3e-4, 1.3e-4, "Array plt")
#else
   call assert_allclose(pk_gpu, size(pk_gpu), pk, size(pk), 1e-8, 1e-8, "Array pk")
   call assert_allclose(pk2_gpu, size(pk2_gpu), pk2, size(pk2), 1e-8, 1e-8, "Array pk2")
   call assert_allclose_r8(plt_gpu, size(plt_gpu), plt, size(plt), 1e-8, 1e-8, "Array plt")
#endif

end subroutine
subroutine VarErr(Err, a, b, lev, nvar, lab)
   use const, only: RTYPE, numreduce, weight
   use param, only: nx, my_max, my, octahedral
   use index, only: nxp, nxdef, nxdef_2d, jlist1, jlistnum
   use rank, only: myrank
   use mpe
   implicit none
   ! Measure the difference of between two arrays.
   ! There are the three measurement methods.
   ! 1. Maximum norm
   ! 2. L2-norm (The vertical integral weights need to be corrected)
   ! 3. RMSE
   character(len=*), intent(in):: lab
   real(kind=RTYPE), intent(out):: Err(3, nvar)
   integer, intent(in)::lev, nvar
   real(kind=RTYPE), intent(in):: A(nxp, lev, nvar, my_max), &
                                  B(nxp, lev, nvar, my_max)

   integer i, j, k, n, nxj, jj, pts, lons
   real(kind=8) vamax
   real(kind=RTYPE) Err_loc(nvar, 3), avar_amax(nvar, 2), sum_local(nvar), pi, tmp

   if (octahedral) then
      pts = (20 + nx)*my*lev
   elseif (numreduce == -99) then
      pts = nx*my*lev
   end if

   pi = 4.*atan(1.)
   Err = 0.
   Err_loc = 0.
   avar_amax = 0.

   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef_2d(j)
      lons = nxdef(j)
      sum_local = 0.
      do n = 1, nvar
         do k = 1, lev
            do i = 1, nxj
               avar_amax(n, 1) = max(avar_amax(n, 1), abs(A(i, k, n, jj)))
               avar_amax(n, 2) = max(avar_amax(n, 2), abs(B(i, k, n, jj)))

               tmp = (A(i, k, n, jj) - B(i, k, n, jj))/B(i, k, n, jj)

               Err_loc(n, 1) = max(Err_loc(n, 1), abs(tmp))
               sum_local(n) = sum_local(n) + tmp**2
            end do
         end do
         Err_loc(n, 2) = Err_loc(n, 2) + sum_local(n)*(2.*pi/lons)*weight(j)/lev
         Err_loc(n, 3) = Err_loc(n, 3) + sum_local(n)
      end do
   end do

   call mpe_global_max(Err_loc(1, 1), nvar, RTYPE)
   call mpe_global_sum(Err_loc(1, 2), nvar*2, RTYPE)
   call mpe_global_max(avar_amax(1, 1), nvar*2, RTYPE)

   do n = 1, nvar
      Err_loc(n, 2) = sqrt(Err_loc(n, 2))
      Err_loc(n, 3) = sqrt(Err_loc(n, 3)/pts/nvar)
   end do

   if (myrank .eq. 0) then
      do n = 1, nvar
         Err_loc(n, 1) = Err_loc(n, 1)/avar_amax(n, 1)
         write (*, '(1X, A6, i0.2, 1X, 1pe15.7, 4(1pe15.7))') &
            trim(lab), n, Err_loc(n, 1:3), avar_amax(n, 1:2)
      end do
   end if

   return
end subroutine VarErr
