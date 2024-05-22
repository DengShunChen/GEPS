!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_trngra3
   use param
   use const

   implicit none

   call mpe_init
   call cons
   call trngra3_unit
   call mpe_finalize

end program

subroutine trngra3_unit
   use param
   use const, only: RTYPE
   use index

   implicit none

   integer, parameter :: steps = 16
   real(kind=RTYPE), dimension(jtmax) :: cim
   real(kind=RTYPE), dimension(jtrun, my/2, jtmax) :: poly, dpoly
   real(kind=RTYPE), dimension(lev, 2, jtrun, jtmax) :: s
   real(kind=RTYPE), dimension(nxp, levF, my_max) :: dlpl, dtpl, dlpl_gpu, dtpl_gpu
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
      call trngra3(jtrun, jtmax, nx, levp, my, my_max, cim, poly, dpoly, s, dlpl, dtpl, nsizey)
   end do
   
   !$acc enter data copyin(mtrundef, jlist1, jlist2, nlist, mlist) async(async_id)
   !$acc enter data copyin(cim, poly, dpoly, s, dlpl_gpu, dtpl_gpu) async(async_id)
   do i = 1, steps
      call trngra3_gpu(jtrun, jtmax, nx, levp, my, my_max, cim, poly, dpoly, s, dlpl_gpu, dtpl_gpu, nsizey)
   end do
   !$acc exit data copyout(dlpl_gpu, dtpl_gpu) delete(cim, poly, dpoly, s) async(async_id)
   !$acc exit data delete(mtrundef, jlist1, jlist2, nlist, mlist) async(async_id)
   !$acc wait(async_id)

   call assert_allclose(dlpl_gpu, size(dlpl_gpu), dlpl, size(dlpl), 1e-10, 1e-10, "Array dlpl")
   call assert_allclose(dtpl_gpu, size(dtpl_gpu), dtpl, size(dtpl), 1e-10, 1e-10, "Array dtpl")
   
end subroutine trngra3_unit
