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
   use const, only: RTYPE, w => weight, cim, onocos, poly, dpoly, polyf, dpolyf
   use index

   implicit none

   integer, parameter :: steps = 10
   real(kind=RTYPE), dimension(nxp, lev, my_max)      :: ut, vt
   real(kind=RTYPE), dimension(levp, 2, jtrun, jtmax) :: vor, div
   real(kind=RTYPE), dimension(levp, 2, jtrun, jtmax) :: vor_gpu, div_gpu
   real(kind=RTYPE) cc(nx + 2, lev, 2, my_max), gwk1(nx + 2, lev, 2, my_max)
   real(kind=RTYPE) :: wcc_fk(levp*2*jtmax*my_max*nsizey*2)
   real :: wc(levp*2*my*jtmax*2, 2)
   real :: ws(levp*2*jtrun*jtmax*2)
   real :: fj_weight((jtrun + nsizey)*(my/2)*jtmax, 2)
   integer i, async_id
   async_id = 1

   call random_seed()
   call random_number(ut)
   call random_number(vt)
   vor = 0.
   div = 0.
   vor_gpu = 0.
   div_gpu = 0.

   do i = 1, steps
      call trandv(jtrun, jtmax, nx, my, my_max, lev, ut, vt, w, cim, onocos, poly, dpoly, vor, div, nsizey)
   end do
   !$acc enter data create(cc, gwk1, ws, wc, wcc_fk, fj_weight) async(async_id)
   !$acc enter data copyin(ut, vt, w, cim, onocos, polyf, dpolyf,&
   !$acc& vor_gpu, div_gpu, nlist, jlist2, mlist, &
   !$acc& mtrundef, jlist1, nxjlen, nxjlen_all, nxdef, tcolt_jlist, poly_mlist) async(async_id)
   !$acc wait(async_id)
   do i = 1, steps
      call trandv_gpu_cuda_graph(jtrun, jtmax, nx, my, my_max, levp, &
                                 ut, vt, w, cim, onocos, polyf, dpolyf, &
                                 vor_gpu, div_gpu, nsizey, &
                                 cc, gwk1, ws, wc(1, 1), wc(1, 2), &
                                 wcc_fk, fj_weight(1, 1), fj_weight(1, 2))
      !$acc wait(async_id)
   end do
   !$acc exit data copyout(vor_gpu, div_gpu) &
   !$acc& delete(vt, ut, w, cim, onocos, polyf, dpolyf, nlist, jlist2, mlist,&
   !$acc& mtrundef, jlist1, nxjlen, nxjlen_all, nxdef, tcolt_jlist, poly_mlist) async(async_id)
   !$acc exit data delete(cc, gwk1, ws, wc, wcc_fk, fj_weight) async(async_id)
   !$acc wait(async_id)

#ifdef SP
   call assert_rmse(div_gpu, size(div_gpu), div, size(div), 1e-4_4, "Array div")
   call assert_rmse(vor_gpu, size(vor_gpu), vor, size(vor), 1e-4_4, "Array vor")
#else
   call assert_allclose(div_gpu, size(div_gpu), div, size(div), 1e-10, 1e-10, "Array div")
   call assert_allclose(vor_gpu, size(vor_gpu), vor, size(vor), 1e-10, 1e-10, "Array vor")
#endif

end subroutine
