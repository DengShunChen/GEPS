!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_trandv
   use param
   use const

   implicit none

   call mpe_init
   call cons
   call trandv_unit
   call mpe_finalize

end program

subroutine trandv_unit
   use param
   use const, only: RTYPE
   use index

   implicit none

   integer, parameter :: steps = 10
   real(kind=RTYPE), dimension(nxp, lev, my_max)      :: ut, vt
   real(kind=RTYPE), dimension(my)                    :: w, onocos
   real(kind=RTYPE), dimension(jtmax)                 :: cim
   real(kind=RTYPE), dimension(jtrun, my/2, jtmax)    :: poly, dpoly
   real(kind=RTYPE), dimension(levp, 2, jtrun, jtmax) :: vor, div
   real(kind=RTYPE), dimension(levp, 2, jtrun, jtmax) :: vor_gpu, div_gpu
   integer :: i

   call random_seed()
   call random_number(ut)
   call random_number(vt)
   call random_number(w)
   call random_number(onocos)
   call random_number(cim)
   call random_number(poly)
   call random_number(dpoly)
   vor = 0.
   div = 0.
   vor_gpu = 0.
   div_gpu = 0.

   do i = 1, steps
      call trandv(jtrun, jtmax, nx, my, my_max, lev, ut, vt, w, cim, onocos, poly, dpoly, vor, div, nsizey)
   end do
   do i = 1, steps
      call trandv_gpu(jtrun, jtmax, nx, my, my_max, lev, ut, vt, w, cim, onocos, poly, dpoly, vor_gpu, div_gpu, nsizey)
   end do

   if (all(abs(vor - vor_gpu) <= 1e-10) .AND. all(abs(div - div_gpu) <= 1e-10)) then
      print *, "test_trandv passed."
   else
      print *, "test_trandv failed."
      call exit(1)
   end if

end subroutine
