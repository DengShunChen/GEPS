!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_gridnl_hybrid_ndsl
   use param
   use const

   implicit none

   call mpe_init
   call cons
   call gridnl_hybrid_ndsl_unit
   call mpe_finalize

end program

subroutine gridnl_hybrid_ndsl_unit
   !$acc routine(gridnl_hybrid_ndsl_gpu) vector
   use param
   use const
   use index
   use grid, only: latpart

   implicit none

   integer, parameter :: steps = 5
   integer :: i, jj, j, nxj, async_id, seed_size
   integer, allocatable :: seed(:)
   real(kind=RTYPE) um(nxp, lev, my_max)
   real(kind=RTYPE) vm(nxp, lev, my_max)
   real(kind=RTYPE) rdivm(nxp, lev, my_max)
   real(kind=RTYPE) tm(nxp, lev, my_max)
   real(kind=RTYPE) qt(nxp, lev*ncld, my_max)
   real(kind=RTYPE) phi(nxp, lev, my_max), phi_gpu(nxp, lev, my_max)
   real(kind=RTYPE) ptm(nxp, my_max)
   real(kind=RTYPE) dtpl(nxp, my_max)
   real(kind=RTYPE) dlpl(nxp, my_max)
   real(kind=RTYPE) pk(nxp, lev, my_max)
   real(kind=RTYPE) pk2(nxp, lev, my_max)
   real(kind=RTYPE), dimension(nxp, lev, my_max) :: diveng, diveng_gpu
   real(kind=RTYPE) vdmerdg(nxp, lev, my_max), vdmerdg_gpu(nxp, lev, my_max)
   real(kind=RTYPE) vdzonlg(nxp, lev, my_max), vdzonlg_gpu(nxp, lev, my_max)
   real(kind=RTYPE) pten(nxp, lev, my_max)
   real(kind=RTYPE) deldm(nxp, my_max), deldm_gpu(nxp, my_max)
   real(kind=RTYPE) sdpbl(nxp, my_max), sdpbl_gpu(nxp, my_max)
   real(kind=RTYPE) sd(nxp, lev, my_max), sd_gpu(nxp, lev, my_max)
   !! << >>
   real(kind=RTYPE) pdot(nxp, lev + 1, latpart), pdot_gpu(nxp, lev + 1, latpart)
   !! << >>
   real(kind=RTYPE) vvel(nxp, lev, my_max), vvel_gpu(nxp, lev, my_max)
   real(kind=RTYPE) sgeo(nxp, my_max)
   real(kind=RTYPE) abs_diff, rel_diff

   async_id = 1

   call random_seed(size=seed_size)
   allocate (seed(seed_size))
   seed = 123
   call random_seed(put=seed)
   call random_number(um)
   call random_number(vm)
   call random_number(rdivm)
   call random_number(tm)
   call random_number(qt)
   call random_number(ptm)
   call random_number(dtpl)
   call random_number(dlpl)
   call random_number(pk)
   call random_number(pk2)
   call random_number(pten)
   call random_number(vvel)
   call random_number(sgeo)
   deldm = 0.
   sd = 0.
   pdot = 0.
   vvel = 0.
   diveng = 0.
   vdmerdg = 0.
   vdzonlg = 0.
   sdpbl = 0.
   phi = 0.
   deldm_gpu = 0.
   sd_gpu = 0.
   pdot_gpu = 0.
   vvel_gpu = 0.
   diveng_gpu = 0.
   vdmerdg_gpu = 0.
   vdzonlg_gpu = 0.
   sdpbl_gpu = 0.
   phi_gpu = 0.

   do i = 1, steps
      do jj = 1, jlistnum
         j = jlist1(jj)
         nxj = nxdef_2d(j)
         call gridnl_hybrid_ndsl(nxjp(j), nxp, lev, ncld &
                                 , cp, radsq, um(1, 1, jj), vm(1, 1, jj), rdivm(1, 1, jj), tm(1, 1, jj) &
                                 , qt(1, 1, jj), phi(1, 1, jj), ptm(1, jj), dtpl(1, jj), dlpl(1, jj), sinl(j) &
                                 , pk(1, 1, jj), pk2(1, 1, jj), dsigma, sigma, onocos(j), cor(j) &
                                 , diveng(1, 1, jj), vdmerdg(1, 1, jj), vdzonlg(1, 1, jj), pten(1, 1, jj) &
                                 , deldm(1, jj), sdpbl(1, jj), sd(1, 1, jj), pdot(1, 1, jj), vvel(1, 1, jj) &
                                 , sgeo(1, jj))
      end do
   end do

   !$acc enter data copyin(um, vm, rdivm, tm, qt, phi, phi_gpu, ptm, dtpl, dlpl, sinl, pk, pk2, onocos, cor, diveng_gpu, vdmerdg_gpu, vdzonlg_gpu, pten, deldm_gpu, sdpbl_gpu, sd_gpu, pdot_gpu, vvel_gpu, sgeo, dsigma, sigma, nxjp, jlist1) async(async_id)
   do i = 1, steps
      ! !$acc parallel loop gang async(async_id) private(j, nxj)
      ! do jj = 1, jlistnum
      !    j = jlist1(jj)
      !    nxj = nxdef_2d(j)
      !    call gridnl_hybrid_ndsl_gpu(nxjp(j), nxp, lev, ncld &
      !                                , cp, radsq, um(1, 1, jj), vm(1, 1, jj), rdivm(1, 1, jj), tm(1, 1, jj) &
      !                                , qt(1, 1, jj), phi_gpu(1, 1, jj), ptm(1, jj), dtpl(1, jj), dlpl(1, jj), sinl(j) &
      !                                , pk(1, 1, jj), pk2(1, 1, jj), dsigma, sigma, onocos(j), cor(j) &
      !                                , diveng_gpu(1, 1, jj), vdmerdg_gpu(1, 1, jj), vdzonlg_gpu(1, 1, jj), pten(1, 1, jj) &
      !                               , deldm_gpu(1, jj), sdpbl_gpu(1, jj), sd_gpu(1, 1, jj), pdot_gpu(1, 1, jj), vvel_gpu(1, 1, jj) &
      !                                , sgeo(1, jj))
      ! end do
      ! ------------------------------------------------------------
      call gridnl_hybrid_ndsl_gpu_refactor(nxjp, nxp, lev, ncld &
           , cp, radsq, um, vm, rdivm, tm &
           , qt, phi_gpu, ptm, dtpl, dlpl, sinl &
           , pk, pk2, dsigma, sigma, onocos, cor &
           , diveng_gpu, vdmerdg_gpu, vdzonlg_gpu, pten &
           , deldm_gpu, sdpbl_gpu, sd_gpu, pdot_gpu, vvel_gpu &
           , sgeo)
   end do
   !$acc exit data delete(um, vm, rdivm, tm, qt, phi, ptm, dtpl, dlpl, sinl, pk, pk2, onocos, cor, pten, sgeo, nxjp, jlist1) &
   !$acc& copyout(phi_gpu, deldm_gpu, sd_gpu, pdot_gpu, vvel_gpu, diveng_gpu, vdmerdg_gpu, vdzonlg_gpu, dsigma, sigma, sdpbl_gpu) async(async_id)
   !$acc wait(async_id)

#ifdef SP
   call assert_rmse(sd_gpu, size(sd_gpu), sd, size(sd), 1e-4_4, "Array sd")
   call assert_rmse(deldm_gpu, size(deldm_gpu), deldm, size(deldm), 1e-4_4, "Array deldm")
   call assert_rmse(pdot_gpu, size(pdot_gpu), pdot, size(pdot), 1e-4_4, "Array pdot")
   call assert_rmse(vvel_gpu, size(vvel_gpu), vvel, size(vvel), 1e-4_4, "Array vvel")
   call assert_rmse(sdpbl_gpu, size(sdpbl_gpu), sdpbl, size(sdpbl), 1e-4_4, "Array sdpbl")
   call assert_rmse(vdmerdg_gpu, size(vdmerdg_gpu), vdmerdg, size(vdmerdg), 1e-4_4, "Array vdmerdg")
   call assert_rmse(vdzonlg_gpu, size(vdzonlg_gpu), vdzonlg, size(vdzonlg), 1e-4_4, "Array vdzonlg")
   call assert_rmse(diveng_gpu, size(diveng_gpu), diveng, size(diveng), 1e-4_4, "Array diveng")
#else
   call assert_allclose(sd_gpu, size(sd_gpu), sd, size(sd), 1e-10, 1e-10, "Array sd")
   call assert_allclose(deldm_gpu, size(deldm_gpu), deldm, size(deldm), 1e-10, 1e-10, "Array deldm")
   call assert_allclose(pdot_gpu, size(pdot_gpu), pdot, size(pdot), 1e-10, 1e-10, "Array pdot")
   call assert_allclose(vvel_gpu, size(vvel_gpu), vvel, size(vvel), 1e-10, 1e-10, "Array vvel")
   call assert_allclose(sdpbl_gpu, size(sdpbl_gpu), sdpbl, size(sdpbl), 1e-10, 1e-10, "Array sdpbl")
   call assert_allclose(vdmerdg_gpu, size(vdmerdg_gpu), vdmerdg, size(vdmerdg), 1e-10, 1e-10, "Array vdmerdg")
   call assert_allclose(vdzonlg_gpu, size(vdzonlg_gpu), vdzonlg, size(vdzonlg), 1e-10, 1e-10, "Array vdzonlg")
   call assert_allclose(diveng_gpu, size(diveng_gpu), diveng, size(diveng), 1e-10, 1e-10, "Array diveng")
#endif

end subroutine
