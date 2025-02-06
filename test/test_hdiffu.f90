!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_hdiffu
   use param
   use const

   implicit none

   call mpe_init
   call cons
   call hdiffu_unit
   call mpe_finalize

end program

subroutine hdiffu_unit
   use param
   use const
   use index

   implicit none

   integer, parameter :: steps = 16
   real(kind=RTYPE) dta
   real(kind=RTYPE) um(nxp, lev, my_max)
   real(kind=RTYPE) vm(nxp, lev, my_max)
   real(kind=RTYPE) vornow(levp, 2, jtrun, jtmax), vornow_gpu(levp, 2, jtrun, jtmax), vornow_cpu(levp, 2, jtrun, jtmax)
   real(kind=RTYPE) divnow(levp, 2, jtrun, jtmax), divnow_gpu(levp, 2, jtrun, jtmax), divnow_cpu(levp, 2, jtrun, jtmax)
   real(kind=RTYPE) temnow(levp, 2, jtrun, jtmax), temnow_gpu(levp, 2, jtrun, jtmax), temnow_cpu(levp, 2, jtrun, jtmax)
   real(kind=RTYPE) trefs(levp, 2, jtrun, jtmax)
   integer i, seed_size
   integer, allocatable :: seed(:)
   integer async_id

   async_id = 1

   dta = dt
   hfiltx = 0.8*hfilt
   call random_seed(size=seed_size)
   allocate (seed(seed_size))
   seed = 128
   call random_seed(put=seed)
   call random_number(um)
   call random_number(vm)
   call random_number(vornow)
   call random_number(divnow)
   call random_number(temnow)
   call random_number(trefs)

   do i = 1, steps
      vornow_cpu = vornow
      divnow_cpu = divnow
      temnow_cpu = temnow
      call hdiffu(dta, my, my_max, nx, jtrun, jtmax, lev, ncld &
                  , hfiltx, rad, cosl, um, vm, vornow_cpu, divnow_cpu, temnow_cpu &
                  , eps4, trefs)
   end do

   !$acc enter data copyin(cosl, um, vm, vornow_gpu, divnow_gpu, temnow_gpu, &
   !$acc& eps4, trefs, jlist1, nxdef_2d, Llist, hdk2, vornow, divnow, temnow, mlist) async(async_id)
   do i = 1, steps
      !$acc kernels async(async_id)
      vornow_gpu = vornow
      divnow_gpu = divnow
      temnow_gpu = temnow
      !$acc end kernels
      ! Present on device: cosl, ut, vt, vornow, divnow, temnow, eps4, trefs
      ! Present on device: jlist1, nxdef_2d, Llist, hdk2
      call hdiffu_gpu(dta, my, my_max, nx, jtrun, jtmax, lev, ncld &
                      , hfiltx, rad, cosl, um, vm, vornow_gpu, divnow_gpu, temnow_gpu &
                      , eps4, trefs)
   end do
   !$acc exit data copyout(cosl, um, vm, vornow_gpu, divnow_gpu, temnow_gpu, &
   !$acc& eps4, trefs, jlist1, nxdef_2d, Llist, hdk2, vornow, divnow, temnow, mlist) async(async_id)
   !$acc wait(async_id)

#ifdef SP
   call assert_allclose(vornow_gpu, size(vornow_gpu), vornow_cpu, size(vornow_cpu), 1e-4_4, 1e-4_4, "Array vornow")
   call assert_allclose(divnow_gpu, size(divnow_gpu), divnow_cpu, size(divnow_cpu), 1e-4_4, 1e-4_4, "Array divnow")
   call assert_allclose(temnow_gpu, size(temnow_gpu), temnow_cpu, size(temnow_cpu), 1e-4_4, 1e-4_4, "Array temnow")
#else
   call assert_allclose(vornow_gpu, size(vornow_gpu), vornow_cpu, size(vornow_cpu), 1e-10, 1e-10, "Array vornow")
   call assert_allclose(divnow_gpu, size(divnow_gpu), divnow_cpu, size(divnow_cpu), 1e-10, 1e-10, "Array divnow")
   call assert_allclose(temnow_gpu, size(temnow_gpu), temnow_cpu, size(temnow_cpu), 1e-10, 1e-10, "Array temnow")
#endif

end subroutine hdiffu_unit
