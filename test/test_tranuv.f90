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
   use const, only: RTYPE, onocos, wcfac, wdfac, poly, dpoly, polyf, dpolyf, coslr
   use index

   implicit none

   integer, parameter :: steps = 5
   ! real(kind=RTYPE), dimension(my) :: onocos
   ! real(kind=RTYPE), dimension(jtrun, jtmax) :: wcfac, wdfac
   ! real(kind=RTYPE), dimension(jtrun, my/2, jtmax) :: poly, dpoly
   real(kind=RTYPE), dimension(lev, 2, jtrun, jtmax) :: vor, div
   real(kind=RTYPE), dimension(nxp, levF, my_max) :: ut, vt, ut_gpu, vt_gpu
   real(kind=RTYPE), dimension(nx + 2, levp, 2, my_max) :: cc, gwk1
   real, dimension(lev, 2, 2, jtrun, jtmax):: ws1, ws2
   real, dimension(lev, 2, 2, my, jtmax) :: tcc
   real(kind=RTYPE), dimension(lev, 2, 2, jtmax, my_max*nsizey):: wcc_fk
   real, dimension(jtrun + nsizey, my/2, jtmax) :: fj_wc, fj_wd
   integer :: i, async_id

   async_id = 1

   call random_seed()
   ! call random_number(onocos)
   ! call random_number(wcfac)
   ! call random_number(wdfac)
   ! call random_number(poly)
   ! call random_number(dpoly)
   call random_number(vor)
   call random_number(div)

   ut = 0.
   vt = 0.
   ut_gpu = 0.
   vt_gpu = 0.

   do i = 1, steps
      call tranuv(jtrun, jtmax, nx, my, my_max, levp, onocos, wcfac, wdfac, poly, dpoly, vor, div, ut, vt, nsizey)
   end do
   !$acc enter data copyin(mlistnum,jtrun,lev)
   !$acc enter data copyin(coslr, wcfac, wdfac, polyf, dpolyf, vor, div, ut_gpu, vt_gpu, &
   !$acc& jlist1, nlist, mtrundef, mlist, jlist2, nxdef, tcolt_jlist, poly_mlist, cc, gwk1) async(async_id)
   !$acc enter data create(ws1, ws2, tcc, wcc_fk , fj_wc, fj_wd) async(async_id)
   do i = 1, steps
      call tranuv_gpu_cuda_graph(jtrun, jtmax, nx, my, my_max, levp, &
                                 coslr, wcfac, wdfac, &
                                 polyf, dpolyf, vor, div, ut_gpu, vt_gpu, nsizey, &
                                 cc, gwk1, ws1, ws2, tcc, wcc_fk, fj_wc, fj_wd)
   end do
   !$acc exit data copyout(ut_gpu, vt_gpu) async(async_id)
   !$acc exit data delete(coslr, wcfac, wdfac, polyf, dpolyf, vor, div, &
   !$acc& jlist1, nlist, mtrundef, mlist, jlist2, nxdef, tcolt_jlist, poly_mlist, cc, gwk1) async(async_id)
   !$acc exit data delete(ws1, ws2, tcc, wcc_fk , fj_wc, fj_wd) async(async_id)
   !$acc wait(async_id)

#ifdef SP
   call assert_rmse(ut_gpu, size(ut_gpu), ut, size(ut), 1e-5_4, "Array ut")
   call assert_rmse(vt_gpu, size(vt_gpu), vt, size(vt), 1e-5_4, "Array vt")
#else
   call assert_allclose(ut_gpu, size(ut_gpu), ut, size(ut), 1e-8, 1e-8, "Array ut")
   call assert_allclose(vt_gpu, size(vt_gpu), vt, size(vt), 1e-8, 1e-8, "Array vt")
#endif

end subroutine tranuv_unit
