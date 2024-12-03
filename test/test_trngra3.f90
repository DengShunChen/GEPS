!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_trngra3

   call mpe_init
   call cons
   call trngra3_unit
   call mpe_finalize

end program

subroutine trngra3_unit
   use const, only: RTYPE, polyf, dpolyf, poly, dpoly, cim
   use index
   use param

   use rank
   implicit none

   integer, parameter :: steps = 16
   ! real(kind=RTYPE), dimension(jtmax) :: cim
   ! real(kind=RTYPE), dimension(jtrun, my/2, jtmax) :: poly, dpoly
   real(kind=RTYPE), dimension(levp, 2, jtrun, jtmax) :: s
   real(kind=RTYPE), dimension(nx + 2, levp, 2, my_max) :: cc, gwk1
   real, dimension(levp, 2, jtrun, jtmax, 2) :: ws
   real, dimension(levp, 2, my, jtmax, 2) :: tcc
   real, dimension(levp, 2, 2, jtmax, my_max*nsize) :: wcc_fk
   real(kind=RTYPE), dimension(nxp, levF, my_max) :: dlpl, dtpl, dlpl_gpu, dtpl_gpu
   integer :: i, async_id, m, mf, l, k
   integer :: m_str, j_str, jlistnum_fj, llistnum_fj, ind, j, jj

   async_id = 1

   call init_seed()
   call random_number(cim)
   call random_number(poly)
   call random_number(dpoly)
   call random_number(s)

   do m = 1, mlistnum
      m_str = poly_mlist(m)
      jlistnum_fj = tcolt_jlist(1, m)
      j_str = tcolt_jlist(2, m)

      mf = mlist(m)
      llistnum_fj = jtrun - mf + 1
      do j = 1, jlistnum_fj
         do l = mf, jtrun
            ind = (l - mf + 1) + (j - 1)*llistnum_fj + m_str - 1
            jj = j_str + j - 1
            polyf(ind) = poly(l, jj, m)
            dpolyf(ind) = dpoly(l, jj, m)

            ind = ind + jlistnum_fj*llistnum_fj
            jj = my/2 - j + 1
            polyf(ind) = (-1)**(l - mf)*poly(l, jj, m)
            dpolyf(ind) = (-1)**(l - mf + 1)*dpoly(l, jj, m)
         end do

      end do
   end do

   dlpl = 0.
   dtpl = 0.
   dlpl_gpu = 0.
   dtpl_gpu = 0.

   do i = 1, steps
      call trngra3(jtrun, jtmax, nx, levp, my, my_max, cim, poly, dpoly, s, dlpl, dtpl, nsizey)
   end do

   !$acc enter data copyin(mtrundef, jlist1, jlist2, nlist, mlist, nxdef, &
   !$acc& tcolt_jlist, poly_mlist) async(async_id)
   !$acc enter data copyin(cim, poly, dpoly, polyf, dpolyf, s, dlpl_gpu, dtpl_gpu) async(async_id)
   !$acc enter data create(cc, gwk1, ws, tcc, wcc_fk) async(async_id)
   !$acc wait(async_id)
   do i = 1, steps
      ! call trngra3_gpu(jtrun, jtmax, nx, levp, my, my_max, &
      !                             cim, poly, dpoly, &
      !                             s, dlpl_gpu, dtpl_gpu, nsizey, &
      !                             cc, gwk1)
      call trngra3_gpu_cuda_graph(jtrun, jtmax, nx, levp, my, my_max, &
                                  cim, polyf, dpolyf, &
                                  s, dlpl_gpu, dtpl_gpu, nsizey, &
                                  cc, gwk1, ws, tcc, wcc_fk)
   end do
   !$acc exit data copyout(dlpl_gpu, dtpl_gpu) delete(cim, poly, dpoly, polyf, dpolyf, s, cc, gwk1) async(async_id)
   !$acc exit data delete(mtrundef, jlist1, jlist2, nlist, mlist, nxdef, tcolt_jlist, poly_mlist, &
   !$acc& ws, tcc, wcc_fk) async(async_id)
   !$acc wait(async_id)

#ifdef SP
   call assert_rmse(dlpl_gpu, size(dlpl_gpu), dlpl, size(dlpl), 1e-3_4, "Array dlpl")
   call assert_rmse(dtpl_gpu, size(dtpl_gpu), dtpl, size(dtpl), 1e-3_4, "Array dtpl")
#else
   call assert_allclose(dlpl_gpu, size(dlpl_gpu), dlpl, size(dlpl), 1e-10, 1e-10, &
                        "(allclose)Array dlpl")
   call assert_allclose(dtpl_gpu, size(dtpl_gpu), dtpl, size(dtpl), 1e-10, 1e-10, &
                        "(allclose)Array dtpl")
#endif

end subroutine trngra3_unit
