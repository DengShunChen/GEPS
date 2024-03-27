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
   call transr_unit
   call mpe_finalize

end program

subroutine transr_unit
   use param
   use const, only: RTYPE
   use index

   implicit none

   integer, parameter :: steps = 10
   integer, parameter :: num = 1
   integer :: i
   real(kind=RTYPE) :: poly(jtrun, my/2, jtmax)
   real(kind=RTYPE) :: wss(lev, 2, num, jtrun, jtmax)
   real(kind=RTYPE) :: cc(nx + 2, lev, num, my_max)
   real(kind=RTYPE) :: cc_gpu(nx + 2, lev, num, my_max)

   call random_seed()
   call random_number(poly)
   call random_number(wss)
   cc = 0.
   cc_gpu = 0.

   do i = 1, steps
      call transr(jtrun, jtmax, nx, my, my_max, levp, poly, wss, cc, num, nsizey)
   end do
   do i = 1, steps
      call transr_gpu(jtrun, jtmax, nx, my, my_max, levp, poly, wss, cc_gpu, num, nsizey)
   end do

   if (all(abs(cc - cc_gpu) <= 1e-10)) then
      print *, "test_transr passed."
   else
      print *, "test_transr failed."
      call exit(1)
   end if

end subroutine
