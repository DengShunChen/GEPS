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
   use const, only: RTYPE, dt, itter, ptmeans, alpha
   use index

   implicit none

   integer, parameter :: steps = 1
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

   call siimpl(jtrun, jtmax, lev, dtahi, ptmeans, dsigma, spalm, eps4, eigval &
               , evecin, evectr, arrhyd, arsddt, temmid, divmid, plmid &
               , temmid, divmid, plmid, temten, divten, plten, alpha)

   call siimpl_gpu(jtrun, jtmax, lev, dtahi, ptmeans, dsigma, spalm, eps4, eigval &
                   , evecin, evectr, arrhyd, arsddt, temmid, divmid, plmid &
                   , temmid, divmid, plmid, temten_gpu, divten_gpu, plten_gpu, alpha)

   abs_err = maxval(abs(temten - temten_gpu))
   rel_err = maxval(abs((temten - temten_gpu)/(temten + 1e-15)))
   print *, abs_err, rel_err
   if (all(abs(temten - temten_gpu) <= 1e-10)) then
      print *, "test_siimpl temten passed."
   else
      print *, "test_siimpl temten failed."
      ! call exit(1)
   end if
   if (all(abs(divten - divten_gpu) <= 1e-10)) then
      print *, "test_siimpl divten passed."
   else
      print *, "test_siimpl divten failed."
      ! call exit(1)
   end if
   if (all(abs(plten - plten_gpu) <= 1e-10)) then
      print *, "test_siimpl plten passed."
   else
      print *, "test_siimpl plten failed."
      ! call exit(1)
   end if

end subroutine siimpl_unit
