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
   real(kind=RTYPE) cim(jtmax)
   real(kind=RTYPE) poly(jtrun, my/2, jtmax), dpoly(jtrun, my/2, jtmax)
   real(kind=RTYPE) s(lev, 2, jtrun, jtmax)
   real(kind=RTYPE) dlpl(nxp, levF, my_max), dtpl(nxp, levF, my_max)
   real(kind=RTYPE) dlpl_gpu(nxp, levF, my_max), dtpl_gpu(nxp, levF, my_max)
   integer i
   integer async_id

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

   do i = 1, steps
      !$acc enter data copyin(cim, poly, dpoly, s, dlpl_gpu, dtpl_gpu) async(async_id)
      !$acc wait(async_id)
      call trngra3_gpu(jtrun, jtmax, nx, levp, my, my_max, cim, poly, dpoly, s, dlpl_gpu, dtpl_gpu, nsizey)
      !$acc wait(async_id)
      !$acc exit data copyout(dlpl_gpu, dtpl_gpu) delete(cim, poly, dpoly, s) async(async_id)
      !$acc wait(async_id)
   end do

   if (all(abs(dlpl - dlpl_gpu) <= 1e-10) .AND. all(abs(dtpl - dtpl_gpu) <= 1e-10)) then
      PRINT *, "test_trngra3 passed."
   else
      PRINT *, "test_trngra3 failed."
      call exit(1)
   end if

end subroutine trngra3_unit
