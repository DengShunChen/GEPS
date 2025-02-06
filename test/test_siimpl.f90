!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_siimpl
   use param
   use const

   implicit none

   call mpe_init
   call cons
   call siimpl_unit
   call mpe_finalize

end program

subroutine siimpl_unit
   use param
   use const, only: RTYPE, dt, itter, ptmeans, alpha, eps4L
   use spec, only: jtwvp
   use index

   implicit none

   integer, parameter :: steps = 16
   real(kind=RTYPE) dtahi
   real(kind=RTYPE) dsigma(lev, 2), spalm(lev), eps4(jtrun, jtmax), eigval(lev), evecin(lev, lev), evectr(lev, lev), arrhyd(lev, lev), arsddt(lev, lev)
   real(kind=RTYPE) temmid(levp, 2, jtrun, jtmax), divmid(levp, 2, jtrun, jtmax), plmid(jtrun, jtmax, 2)
   real(kind=RTYPE) temten(levp, 2, jtrun, jtmax), divten(levp, 2, jtrun, jtmax), plten(jtrun, jtmax, 2)
   real(kind=RTYPE) temten_gpu(levp, 2, jtrun, jtmax), divten_gpu(levp, 2, jtrun, jtmax), plten_gpu(jtrun, jtmax, 2)
   integer i, seed_size
   integer, allocatable :: seed(:)
   integer async_id
   real(kind=RTYPE) abs_err, rel_err

   async_id = 1

   dtahi = dt*0.5/float(itter)
   call random_seed(size=seed_size)
   allocate (seed(seed_size))
   seed = 128
   call random_seed(put=seed)
   call random_number(dsigma)
   call random_number(spalm)
   call random_number(eps4)
   call random_number(eigval)
   call random_number(evecin)
   call random_number(evectr)
   call random_number(arrhyd)
   call random_number(arsddt)
   call random_number(temmid)
   call random_number(divmid)
   call random_number(plmid)
   call random_number(temten)
   call random_number(divten)
   call random_number(plten)
   temten_gpu = temten
   divten_gpu = divten
   plten_gpu = plten

   do i = 1, steps
      call siimpl(jtrun, jtmax, lev, dtahi, ptmeans, dsigma, spalm, eps4, eigval &
                  , evecin, evectr, arrhyd, arsddt, temmid, divmid, plmid &
                  , temmid, divmid, plmid, temten, divten, plten, alpha)
   end do

   !$acc enter data copyin(plten_gpu, plmid, temmid, temten_gpu, divmid, divten_gpu, &
   !$acc& jtwvp, spalm, arrhyd, eps4L, evecin, eigval, evectr, arsddt, dsigma, mlist) async(async_id)
   do i = 1, steps
      call siimpl_gpu(jtrun, jtmax, lev, dtahi, ptmeans, dsigma, spalm, eps4, eigval &
                     , evecin, evectr, arrhyd, arsddt, temmid, divmid, plmid &
                     , temmid, divmid, plmid, temten_gpu, divten_gpu, plten_gpu, alpha)
   end do
   !$acc exit data copyout(plten_gpu, plmid, temmid, temten_gpu, divmid, divten_gpu, &
   !$acc& jtwvp, spalm, arrhyd, eps4L, evecin, eigval, evectr, arsddt, dsigma, mlist) async(async_id)
   !$acc wait(async_id)

#ifdef SP
   ! call assert_rmse(temten_gpu, size(temten_gpu), temten, size(temten), 1e-4_4, "Array temten")
   ! call assert_rmse(divten_gpu, size(divten_gpu), divten, size(divten), 1e-4_4, "Array divten")
   ! call assert_rmse(plten_gpu, size(plten_gpu), plten, size(plten), 1e-4_4, "Array plten")
   ! call assert_allclose(temten_gpu, size(temten_gpu), temten, size(temten), 1e-3_4, 1e-3_4, "Array temten")
   ! call assert_allclose(divten_gpu, size(divten_gpu), divten, size(divten), 1e-3_4, 1e-3_4, "Array divten")
   ! call assert_allclose(plten_gpu, size(plten_gpu), plten, size(plten), 1e-3_4, 1e-3_4, "Array plten")
#else
   call assert_allclose(temten_gpu, size(temten_gpu), temten, size(temten), 1e-10, 1e-10, "Array temten")
   call assert_allclose(divten_gpu, size(divten_gpu), divten, size(divten), 1e-10, 1e-10, "Array divten")
   call assert_allclose(plten_gpu, size(plten_gpu), plten, size(plten), 1e-10, 1e-10, "Array plten")
#endif

end subroutine siimpl_unit
