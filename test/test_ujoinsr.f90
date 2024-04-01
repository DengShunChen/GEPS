!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_tranuv
   use param
   use const
   
   implicit none
   
   call mpe_init
   call cons
   call tranuv_unit
   call mpe_finalize

end program

subroutine tranuv_unit
   use param
   use const, only: RTYPE
   use index
   
   implicit none
   
   real(kind=RTYPE) cc(nx + 2, levp, 1, my_max)
   real(kind=RTYPE) r1(nxp, lev, my_max)
   real(kind=RTYPE) r1_gpu(nxp, lev, my_max)
   real(kind=RTYPE) dummy
   
   call random_seed()
   call random_number(cc)
   r1 = 0.
   r1_gpu = 0.
   
   call ujoinsr(cc, r1, dummy, dummy, dummy, nx, my_max, lev, jlistnum, 1, 1)
   call ujoinsr_gpu(cc, r1_gpu, dummy, dummy, dummy, nx, my_max, lev, jlistnum, 1, 1)
   print *, r1_gpu(123, 4, 5), r1(123, 4, 5)
   if (all(abs(r1 - r1_gpu) <= 1e-10)) then
      PRINT *, "test_tranuv passed."
   else
      PRINT *, "test_tranuv failed."
      call exit(1)
   end if
   
end subroutine tranuv_unit
