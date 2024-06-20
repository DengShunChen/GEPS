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

   real(kind=RTYPE) onocos(my)
   real(kind=RTYPE) wcfac(jtrun, jtmax), wdfac(jtrun, jtmax)
   real(kind=RTYPE) poly(jtrun, my/2, jtmax), dpoly(jtrun, my/2, jtmax)
   real(kind=RTYPE) vor(lev, 2, jtrun, jtmax), div(lev, 2, jtrun, jtmax)
   real(kind=RTYPE) ut(nxp, levF, my_max), vt(nxp, levF, my_max)
   real(kind=RTYPE) ut_gpu(nxp, levF, my_max), vt_gpu(nxp, levF, my_max)
   integer i, async_id

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

   do i = 1, 16
      call tranuv(jtrun, jtmax, nx, my, my_max, levp, onocos, wcfac, wdfac, poly, dpoly, vor, div, ut, vt, nsizey)
   end do
   do i = 1, 16
      !$acc enter data copyin(onocos, wcfac, wdfac, poly, dpoly, vor, div, ut_gpu, vt_gpu, jlist1, nlist, mtrundef, mlist, jlist2) async(async_id)
      call tranuv_gpu(jtrun, jtmax, nx, my, my_max, levp, onocos, wcfac, wdfac, poly, dpoly, vor, div, ut_gpu, vt_gpu, nsizey)
      !$acc exit data copyout(onocos, wcfac, wdfac, poly, dpoly, vor, div, ut_gpu, vt_gpu, jlist1, nlist, mtrundef, mlist, jlist2) async(async_id)
      !$acc wait(async_id)
   end do

   call assert_allclose(ut_gpu, size(ut_gpu), ut, size(ut), 1e-8, 1e-8, "Array ut")
   call assert_allclose(vt_gpu, size(vt_gpu), vt, size(vt), 1e-8, 1e-8, "Array vt")

end subroutine tranuv_unit
