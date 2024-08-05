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

   integer, parameter :: steps = 16
   real(kind=RTYPE), dimension(my) :: onocos
   real(kind=RTYPE), dimension(jtrun, jtmax) :: wcfac, wdfac
   real(kind=RTYPE), dimension(jtrun, my/2, jtmax) :: poly, dpoly
   real(kind=RTYPE), dimension(lev, 2, jtrun, jtmax) :: vor, div
   real(kind=RTYPE), dimension(nxp, levF, my_max) :: ut, vt, ut_gpu, vt_gpu
   real(kind=RTYPE), dimension(nx + 2, levp, 2, my_max) :: cc, gwk1
   integer :: i, async_id

   async_id = 1

   call random_seed()
   call random_number(onocos)
   call random_number(wcfac)
   call random_number(wdfac)
   call random_number(poly)
   call random_number(dpoly)
   call random_number(vor)
   call random_number(div)

   ut = 0.
   vt = 0.
   ut_gpu = 0.
   vt_gpu = 0.

   do i = 1, steps
      call tranuv(jtrun, jtmax, nx, my, my_max, levp, onocos, wcfac, wdfac, poly, dpoly, vor, div, ut, vt, nsizey)
   end do
   !$acc enter data copyin(onocos, wcfac, wdfac, poly, dpoly, vor, div, ut_gpu, vt_gpu, &
   !$acc& jlist1, nlist, mtrundef, mlist, jlist2, nxdef, cc, gwk1) async(async_id)
   do i = 1, steps
      call tranuv_gpu(jtrun, jtmax, nx, my, my_max, levp, onocos, wcfac, wdfac, &
                      poly, dpoly, vor, div, ut_gpu, vt_gpu, nsizey, cc, gwk1)
   end do
   !$acc exit data copyout(onocos, wcfac, wdfac, poly, dpoly, vor, div, ut_gpu, vt_gpu, &
   !$acc& jlist1, nlist, mtrundef, mlist, jlist2, nxdef, cc, gwk1) async(async_id)
   !$acc wait(async_id)

   call assert_allclose(ut_gpu, size(ut_gpu), ut, size(ut), 1e-8, 1e-8, "Array ut")
   call assert_allclose(vt_gpu, size(vt_gpu), vt, size(vt), 1e-8, 1e-8, "Array vt")

end subroutine tranuv_unit
