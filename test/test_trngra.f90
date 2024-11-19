!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_trngra
   use param
   use const

   implicit none

   call mpe_init
   call cons
   call trngra_unit
   call mpe_finalize

end program

subroutine trngra_unit
   use param
   use const, only: RTYPE
   use index

   implicit none

   integer, parameter :: steps = 16
   real(kind=RTYPE), dimension(jtmax) :: cim
   real(kind=RTYPE), dimension(jtrun, my/2, jtmax) :: poly, dpoly
   real(kind=RTYPE), dimension(jtrun, jtmax, 2) :: s
   real(kind=RTYPE), dimension(nxp, my_max) :: dlpl, dtpl, dlpl_gpu, dtpl_gpu
   real(kind=RTYPE), dimension(nx + 2, 2, my_max) :: cc, gwk1
   integer :: i, async_id

   async_id = 1

   call random_seed()
   call random_number(cim)
   call random_number(poly)
   call random_number(dpoly)
   call random_number(s)

   dlpl = 0.
   dtpl = 0.
   dlpl_gpu = 0.
   dtpl_gpu = 0.

   do i = 1, steps
      call trngra(jtrun, jtmax, nx, my, my_max, cim, poly, dpoly, s, dlpl, dtpl, nsizey)
   end do

   !$acc enter data copyin(cim, poly, dpoly, s, dlpl_gpu, dtpl_gpu, cc, gwk1) async(async_id)
   !$acc enter data copyin(mlist, mtrundef, jlist2, jlist1, nlist, nxjlen, nxjstart, nxjend, nxdef) async(async_id)
   do i = 1, steps
      call trngra_gpu(jtrun, jtmax, nx, my, my_max, cim, poly, dpoly, s, dlpl_gpu, dtpl_gpu, nsizey, cc, gwk1)
   end do
   !$acc exit data delete(mlist, mtrundef, jlist2, jlist1, nlist, nxjlen, nxjstart, nxjend, nxdef) async(async_id)
   !$acc exit data copyout(dlpl_gpu, dtpl_gpu) delete(cim, poly, dpoly, s, cc, gwk1) async(async_id)
   !$acc wait(async_id)

#ifdef SP
   call assert_rmse(dlpl_gpu, size(dlpl_gpu), dlpl, size(dlpl), 1e-3_4, "Array dlpl")
   call assert_rmse(dtpl_gpu, size(dtpl_gpu), dtpl, size(dtpl), 1e-3_4, "Array dtpl")
#else
   call assert_allclose(dlpl_gpu, size(dlpl_gpu), dlpl, size(dlpl), 1e-10, 1e-10, "Array dlpl")
   call assert_allclose(dtpl_gpu, size(dtpl_gpu), dtpl, size(dtpl), 1e-10, 1e-10, "Array dtpl")
#endif

end subroutine trngra_unit
