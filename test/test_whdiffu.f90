!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_whdiffu
   use param
   use const

   implicit none

   call mpe_init
   call cons
   call whdiffu_unit
   call mpe_finalize

end program

subroutine whdiffu_unit
   use param
   use const
   use index

   implicit none

   integer, parameter :: steps = 16
   real(kind=RTYPE) dtah
   real(kind=RTYPE) hfiltm
   real(kind=RTYPE) um(nxp, lev, my_max)
   real(kind=RTYPE) vm(nxp, lev, my_max)
   real(kind=RTYPE) vormid(levp, 2, jtrun, jtmax), vormid_gpu(levp, 2, jtrun, jtmax), vormid_cpu(levp, 2, jtrun, jtmax)
   real(kind=RTYPE) divmid(levp, 2, jtrun, jtmax), divmid_gpu(levp, 2, jtrun, jtmax), divmid_cpu(levp, 2, jtrun, jtmax)
   real(kind=RTYPE) temmid(levp, 2, jtrun, jtmax), temmid_gpu(levp, 2, jtrun, jtmax), temmid_cpu(levp, 2, jtrun, jtmax)
   real(kind=RTYPE) trefs(levp, 2, jtrun, jtmax)
   integer i, seed_size
   integer, allocatable :: seed(:)
   integer async_id

   async_id = 1

   dtah = 0.5*dt
   hfiltx = 0.8*hfilt
   hfiltm = mwhd*hfiltx
   call random_seed(size=seed_size)
   allocate (seed(seed_size))
   seed = 128
   call random_seed(put=seed)
   call random_number(um)
   call random_number(vm)
   call random_number(vormid)
   call random_number(divmid)
   call random_number(temmid)
   call random_number(trefs)

   do i = 1, steps
      vormid_cpu = vormid
      divmid_cpu = divmid
      temmid_cpu = temmid
      call whdiffu(dtah, my, my_max, nx, jtrun, jtmax, lev, ncld &
                   , hfiltm, rad, cosl, um, vm, vormid_cpu, divmid_cpu, temmid_cpu &
                   , eps4, trefs)
   end do

   !$acc enter data copyin(cosl, um, vm, vormid_gpu, divmid_gpu, temmid_gpu, eps4, trefs, jlist1, nxdef_2d, Llist, hdk2, &
   !$acc& vormid, divmid, temmid, mlist) async(async_id)
   do i = 1, steps
      !$acc kernels async(async_id)
      vormid_gpu = vormid
      divmid_gpu = divmid
      temmid_gpu = temmid
      !$acc end kernels
      call whdiffu_gpu(dtah, my, my_max, nx, jtrun, jtmax, lev, ncld &
                       , hfiltm, rad, cosl, um, vm, vormid_gpu, divmid_gpu, temmid_gpu &
                       , eps4, trefs)
   end do
   !$acc exit data copyout(cosl, um, vm, vormid_gpu, divmid_gpu, temmid_gpu, eps4, trefs, jlist1, nxdef_2d, Llist, hdk2) &
   !$acc& delete(vormid, divmid, temmid, mlist) async(async_id)
   !$acc wait(async_id)

#ifdef SP
   call assert_allclose(vormid_gpu, size(vormid_gpu), vormid_cpu, size(vormid_cpu), 1e-4_4, 1e-4_4, "Array vormid")
   call assert_allclose(divmid_gpu, size(divmid_gpu), divmid_cpu, size(divmid_cpu), 1e-4_4, 1e-4_4, "Array divmid")
   call assert_allclose(temmid_gpu, size(temmid_gpu), temmid_cpu, size(temmid_cpu), 1e-4_4, 1e-4_4, "Array temmid")
#else
   call assert_allclose(vormid_gpu, size(vormid_gpu), vormid_cpu, size(vormid_cpu), 1e-10, 1e-10, "Array vormid")
   call assert_allclose(divmid_gpu, size(divmid_gpu), divmid_cpu, size(divmid_cpu), 1e-10, 1e-10, "Array divmid")
   call assert_allclose(temmid_gpu, size(temmid_gpu), temmid_cpu, size(temmid_cpu), 1e-10, 1e-10, "Array temmid")
#endif

end subroutine whdiffu_unit
