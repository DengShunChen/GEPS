!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_ujoinsr
   use param
   use const

   implicit none

   call mpe_init
   call cons
   call ujoinsr_unit
   call mpe_finalize

end program

subroutine ujoinsr_unit
   use param
   use const, only: RTYPE
   use index

   implicit none

   integer, parameter :: steps = 10
   real(kind=RTYPE) cc(nx + 2, levp, 1, my_max)
   real(kind=RTYPE) r1(nxp, lev, my_max)
   real(kind=RTYPE) r1_gpu(nxp, lev, my_max)
   real(kind=RTYPE) dummy
   integer async_id
   integer i

   async_id = 1

   call random_seed()
   call random_number(cc)
   r1 = 0.
   r1_gpu = 0.

   do i = 1, steps
      call ujoinsr(cc, r1, dummy, dummy, dummy, nx, my_max, lev, jlistnum, 1, 1)
   end do

   do i = 1, steps
      !$acc enter data copyin(cc, r1_gpu) async(async_id)
      call ujoinsr_gpu(cc, r1_gpu, dummy, dummy, dummy, nx, my_max, lev, jlistnum, 1, 1)
      !$acc wait(async_id)
      !$acc exit data copyout(r1_gpu) delete(cc) async(async_id)
      !$acc wait(async_id)
   end do
   if (all(abs(r1 - r1_gpu) <= 1e-10)) then
      PRINT *, "(all close) test_ujoinsr passed."
   else
      PRINT *, "(all close) test_ujoinsr failed."
      call exit(1)
   end if

end subroutine ujoinsr_unit
