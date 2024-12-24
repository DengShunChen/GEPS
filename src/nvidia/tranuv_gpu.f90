!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define CUDACHECK(ierr) call cuda_check_helper(ierr, __FILE__, __LINE__)

subroutine tranuv_gpu(jtrun, jtmax, nx, my, my_max, lev, coslr, wcfac &
                      , wdfac, poly, dpoly, vor, div, ut, vt, nsize, cc, gwk1)
! Present on device: coslr, wcfac, wdfac, poly, dpoly, vor, div, ut, vt, cc, gwk1
! jlist1, nlist, mtrundef, mlist, jlist2, nxdef
!  subroutine to transform vorticity and divergence to velocity
!  components
!
! *** input ***
!
!  jtrun: zonal wavenumber truncation limit
!  jtmax: maximum amount of zonal waves located in each pe
!  nx: e-w dimension no.
!  my: n-s dimension no.
!  my_max: maximum amount of n-s grids located in each pe
!  lev: number of veritical levels to transform
!  coslr: cos(lat)**2
!  wcfac: constants defined in cons
!  wdfac: constants defined in cons
!  poly: legendre polynomials
!  dpoly: d(poly)/d(sin(lat))
!  vor: spectral vorticity
!  div: spectral divergence
!
! *** output ***
!
!  ut: e-w velocity component
!  vt: n-s velocity component
!
!  ****************************************
!
   use const, only: RTYPE

   use index
   use fftcom
   use spec_cuda_graph, only: fft_cg => tranuv_fft_cg
   use openacc
   use cudafor

   implicit none

   integer jtrun, jtmax, nx, my, my_max, lev, nsize
   integer myhalf, lev2, mlx, j, m, mf, l, k, nb, jchk
   integer jje, jlistnum_fj, l_fj, j_fj, llistnum_fj
   integer kk, ll, jj, jx, j2, i, jtrunj, mchk, mm, mp, mlst
   integer mm1, mp1, mlst1, mm2, mp2, mlst2, mm3, mp3, mlst3, nxj, ierr

   real(kind=RTYPE) sa00, sa10, sa20, sa30, dummy

   real(kind=RTYPE) poly(jtrun, my/2, jtmax), dpoly(jtrun, my/2, jtmax), &
      coslr(my)
   real(kind=RTYPE) vor(lev, 2, jtrun, jtmax), div(lev, 2, jtrun, jtmax)
   real(kind=RTYPE) ut(nxp, levF, my_max), vt(nxp, levF, my_max)
   real(kind=RTYPE) gwk1(nx + 2, lev, 2, my_max)
   real(kind=RTYPE) wcc_fk(lev, 2, 2, jtmax, my_max*nsize)
   real(kind=RTYPE) twcc_fk(lev, 2, 2, jtmax*nsize, my_max)
   real(kind=RTYPE) cc(nx + 2, lev, 2, my_max)
   real(kind=RTYPE) tcc(lev, 2, 2, my)
   real(kind=RTYPE) ws3(lev, 2, 2, jtrun)
   real(kind=RTYPE) ws4(lev, 2, 2, jtrun)
   real(kind=RTYPE) wcfac(jtrun, jtmax), wdfac(jtrun, jtmax)

   real(kind=RTYPE) tc2(lev, 2, 2, my)
   real(kind=RTYPE) wc(jtrun, my/2), wd(jtrun, my/2)

   integer jlist_fj(my/2, mlistnum)
   real fj_ws3(lev*2*2, jtrun)
   real fj_ws4(lev*2*2, jtrun)
   real fj_tcc(lev*2*2, my)
   real fj_wc(jtrun, my/2), fj_wd(jtrun, my/2)
   real fj_tc2(lev*2*2, my)
   integer async_id, istat
   integer(kind=cuda_stream_kind) stream
   integer jlistnum_fj_array(mlistnum)

   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   myhalf = my/2
   lev2 = lev*2
   mlx = (jtrun/2)*((jtrun + 1)/2)
   nb = 32
   jchk = iand(myhalf, 1)
   jje = myhalf - jchk
   !$acc enter data create(jlist_fj, jlistnum_fj_array) async(async_id)
   !$acc parallel loop private(mf, jlistnum_fj) async(async_id)
   do m = 1, mlistnum
      mf = mlist(m)
      jlistnum_fj = 0
      !$acc loop seq
      do j = 1, jje
         if (mf .le. mtrundef(j)) then
            jlistnum_fj = jlistnum_fj + 1
            jlist_fj(jlistnum_fj, m) = j
         end if
      end do
      jlistnum_fj_array(m) = jlistnum_fj
   end do
   !$acc exit data copyout(jlistnum_fj_array) wait(async_id)

   !$acc enter data &
   !$acc& create(twcc_fk, ws3, ws4, tcc, tc2, wcc_fk, fj_ws3, fj_ws4, fj_tcc, fj_wc, fj_wd, fj_tc2, wc, wd) async(async_id)
   !$acc host_data use_device(wcc_fk)
   CUDACHECK(cudaMemsetAsync(wcc_fk, real(0.0, RTYPE), size(wcc_fk), stream))
   !$acc end host_data
   do m = 1, mlistnum
      mf = mlist(m)
      !$acc parallel loop collapse(2) async(async_id)
      do j = 1, myhalf
         do l = 1, jtrun
            mf = mlist(m)
            if ((mf .le. mtrundef(j)) .AND. (l .ge. mf)) then
               wc(l, j) = wcfac(l, m)*poly(l, j, m)
               wd(l, j) = wdfac(l, m)*dpoly(l, j, m)*coslr(j)
            end if
         end do
      end do
      !$acc parallel loop collapse(2) async(async_id)
      do l = 1, jtrun
         do k = 1, lev
            if (l .ge. mf) then
               ws3(k, 1, 1, l) = +div(k, 2, l, m)
               ws3(k, 2, 1, l) = -div(k, 1, l, m)
               ws3(k, 1, 2, l) = +vor(k, 2, l, m)
               ws3(k, 2, 2, l) = -vor(k, 1, l, m)

               ws4(k, 1, 1, l) = +vor(k, 1, l, m)
               ws4(k, 2, 1, l) = +vor(k, 2, l, m)
               ws4(k, 1, 2, l) = -div(k, 1, l, m)
               ws4(k, 2, 2, l) = -div(k, 2, l, m)
            end if
         end do
      end do
      !$acc host_data use_device(tcc, fj_ws3, fj_ws4, fj_tcc, fj_wc, fj_wd)
      CUDACHECK(cudaMemsetAsync(tcc, real(0.0, RTYPE), size(tcc), stream))
      CUDACHECK(cudaMemsetAsync(fj_ws3, 0.0, size(fj_ws3), stream))
      CUDACHECK(cudaMemsetAsync(fj_ws4, 0.0, size(fj_ws4), stream))
      CUDACHECK(cudaMemsetAsync(fj_tcc, 0.0, size(fj_tcc), stream))
      CUDACHECK(cudaMemsetAsync(fj_wc, 0.0, size(fj_wc), stream))
      CUDACHECK(cudaMemsetAsync(fj_wd, 0.0, size(fj_wd), stream))
      !$acc end host_data
      !$acc parallel loop collapse(2) private(l_fj) async(async_id)
      do l = 1, jtrun
         do k = 1, lev2*2
            if (l .ge. mf) then
               l_fj = l - mf + 1
               fj_ws3(k, l_fj) = ws3(k, 1, 1, l)
            end if
         end do
      end do
      jlistnum_fj = jlistnum_fj_array(m)
      !$acc parallel loop collapse(2) private(j, l_fj) async(async_id)
      do j_fj = 1, jlistnum_fj
         do l = 1, jtrun
            if (l .ge. mf) then
               j = jlist_fj(j_fj, m)
               l_fj = l - mf + 1
               fj_wc(l_fj, j_fj) = wc(l, j)
            end if
         end do
      end do

      llistnum_fj = jtrun - mf + 1
      call dgemm_async('n', 'n', lev*2*2, jlistnum_fj, llistnum_fj, &
                       1.0d+0, fj_ws3, lev*2*2, fj_wc, &
                       jtrun, 1.0d+0, fj_tcc, lev*2*2, async_id)

      !$acc parallel loop collapse(2) private(j) async(async_id)
      do j_fj = 1, jlistnum_fj
         do k = 1, lev2*2
            j = jlist_fj(j_fj, m)
            tcc(k, 1, 1, j) = tcc(k, 1, 1, j) + fj_tcc(k, j_fj)
         end do
      end do

      !$acc parallel loop collapse(2) private(l_fj) async(async_id)
      do l = 1, jtrun
         do k = 1, lev2*2
            if (l .ge. mf) then
               l_fj = l - mf + 1
               fj_ws4(k, l_fj) = ws4(k, 1, 1, l)
            end if
         end do
      end do

      !$acc parallel loop collapse(2) async(async_id)
      do j_fj = 1, jlistnum_fj
         do l = 1, jtrun
            if (l .ge. mf) then
               j = jlist_fj(j_fj, m)
               l_fj = l - mf + 1
               fj_wd(l_fj, j_fj) = wd(l, j)
            end if
         end do
      end do
      !$acc host_data use_device(fj_tcc)
      CUDACHECK(cudaMemsetAsync(fj_tcc, 0.0, size(fj_tcc), stream))
      !$acc end host_data

      llistnum_fj = jtrun - mf + 1
      call dgemm_async('n', 'n', lev*2*2, jlistnum_fj, llistnum_fj, 1.0d+0, fj_ws4, lev*2*2, fj_wd, &
                       jtrun, 1.0d+0, fj_tcc, lev*2*2, async_id)

      !$acc parallel loop collapse(2) private(j) async(1)
      do j_fj = 1, jlistnum_fj
         do k = 1, lev2*2
            j = jlist_fj(j_fj, m)
            tcc(k, 1, 1, j) = tcc(k, 1, 1, j) + fj_tcc(k, j_fj)
         end do
      end do
!
! odd number
!
      if (jchk .eq. 1) then
         j = myhalf
         if (mf .le. mtrundef(j)) then
            !$acc parallel loop collapse(2) private(sa00, sa10, sa20, sa30) async(async_id)
            do kk = 1, lev2*2, nb
            do ll = mf, jtrun, nb
               do k = kk, min(kk + nb - 1, lev2*2), 4
                  sa00 = tcc(k, 1, 1, j)
                  sa10 = tcc(k + 1, 1, 1, j)
                  sa20 = tcc(k + 2, 1, 1, j)
                  sa30 = tcc(k + 3, 1, 1, j)
                  do l = ll, min(ll + nb - 1, jtrun)
                     sa00 = sa00 + ws3(k, 1, 1, l)*wc(l, j) + ws4(k, 1, 1, l)*wd(l, j)
                     sa10 = sa10 + ws3(k + 1, 1, 1, l)*wc(l, j) + ws4(k + 1, 1, 1, l)*wd(l, j)
                     sa20 = sa20 + ws3(k + 2, 1, 1, l)*wc(l, j) + ws4(k + 2, 1, 1, l)*wd(l, j)
                     sa30 = sa30 + ws3(k + 3, 1, 1, l)*wc(l, j) + ws4(k + 3, 1, 1, l)*wd(l, j)
                  end do
                  tcc(k, 1, 1, j) = sa00
                  tcc(k + 1, 1, 1, j) = sa10
                  tcc(k + 2, 1, 1, j) = sa20
                  tcc(k + 3, 1, 1, j) = sa30
               end do
            end do
            end do
         end if
      end if

      !$acc parallel loop collapse(2) async(async_id)
      do l = mf, jtrun, 2
      do k = 1, lev
         ws3(k, 1, 1, l) = +div(k, 2, l, m)
         ws3(k, 2, 1, l) = -div(k, 1, l, m)
         ws3(k, 1, 2, l) = +vor(k, 2, l, m)
         ws3(k, 2, 2, l) = -vor(k, 1, l, m)
         ws4(k, 1, 1, l) = -vor(k, 1, l, m)
         ws4(k, 2, 1, l) = -vor(k, 2, l, m)
         ws4(k, 1, 2, l) = +div(k, 1, l, m)
         ws4(k, 2, 2, l) = +div(k, 2, l, m)
      end do
      end do

      !$acc parallel loop collapse(2) async(async_id)
      do l = mf + 1, jtrun, 2
         do k = 1, lev
            ws3(k, 1, 1, l) = -div(k, 2, l, m)
            ws3(k, 2, 1, l) = +div(k, 1, l, m)
            ws3(k, 1, 2, l) = -vor(k, 2, l, m)
            ws3(k, 2, 2, l) = +vor(k, 1, l, m)
            ws4(k, 1, 1, l) = +vor(k, 1, l, m)
            ws4(k, 2, 1, l) = +vor(k, 2, l, m)
            ws4(k, 1, 2, l) = -div(k, 1, l, m)
            ws4(k, 2, 2, l) = -div(k, 2, l, m)
         end do
      end do

      !$acc host_data use_device(tc2, fj_ws3, fj_ws4, fj_tc2, fj_wc, fj_wd)
      CUDACHECK(cudaMemsetAsync(tc2, real(0.0, RTYPE), size(tc2), stream))
      CUDACHECK(cudaMemsetAsync(fj_ws3, 0.0, size(fj_ws3), stream))
      CUDACHECK(cudaMemsetAsync(fj_ws4, 0.0, size(fj_ws4), stream))
      CUDACHECK(cudaMemsetAsync(fj_tc2, 0.0, size(fj_tc2), stream))
      CUDACHECK(cudaMemsetAsync(fj_wc, 0.0, size(fj_wc), stream))
      CUDACHECK(cudaMemsetAsync(fj_wd, 0.0, size(fj_wd), stream))
      !$acc end host_data

      !$acc parallel loop collapse(2) private(l_fj) async(async_id)
      do l = 1, jtrun
         do k = 1, lev2*2
            if (l .ge. mf) then
               l_fj = l - mf + 1
               fj_ws3(k, l_fj) = ws3(k, 1, 1, l)
            end if
         end do
      end do

      !$acc parallel loop collapse(2) private(j, l_fj) async(async_id)
      do j_fj = 1, jlistnum_fj
         do l = 1, jtrun
            if (l .ge. mf) then
               j = jlist_fj(j_fj, m)
               l_fj = l - mf + 1
               fj_wc(l_fj, j_fj) = wc(l, j)
            end if
         end do
      end do

      llistnum_fj = jtrun - mf + 1

      call dgemm_async('n', 'n', lev*2*2, jlistnum_fj, llistnum_fj, 1.0d+0, fj_ws3, lev*2*2, fj_wc, &
                       jtrun, 1.0d+0, fj_tc2, lev*2*2, async_id)

      !$acc parallel loop collapse(2) private(j) async(async_id)
      do j_fj = 1, jlistnum_fj
         do k = 1, lev2*2
            j = jlist_fj(j_fj, m)
            tc2(k, 1, 1, j) = tc2(k, 1, 1, j) + fj_tc2(k, j_fj)
         end do
      end do

      !$acc parallel loop collapse(2) private(l_fj) async(async_id)
      do l = 1, jtrun
         do k = 1, lev2*2
            if (l .ge. mf) then
               l_fj = l - mf + 1
               fj_ws4(k, l_fj) = ws4(k, 1, 1, l)
            end if
         end do
      end do

      !$acc parallel loop collapse(2) private(j, l_fj) async(async_id)
      do j_fj = 1, jlistnum_fj
         do l = 1, jtrun
            if (l .ge. mf) then
               j = jlist_fj(j_fj, m)
               l_fj = l - mf + 1
               fj_wd(l_fj, j_fj) = wd(l, j)
            end if
         end do
      end do

      !$acc host_data use_device(fj_tc2)
      CUDACHECK(cudaMemsetAsync(fj_tc2, 0.0, size(fj_tc2), stream))
      !$acc end host_data

      llistnum_fj = jtrun - mf + 1
      call dgemm_async('n', 'n', lev*2*2, jlistnum_fj, llistnum_fj, 1.0d+0, fj_ws4, lev*2*2, fj_wd, &
                       jtrun, 1.0d+0, fj_tc2, lev*2*2, async_id)

      !$acc parallel loop collapse(2) private(j) async(async_id)
      do j_fj = 1, jlistnum_fj
         do k = 1, lev2*2
            j = jlist_fj(j_fj, m)
            tc2(k, 1, 1, j) = tc2(k, 1, 1, j) + fj_tc2(k, j_fj)
         end do
      end do

!
! odd number
!
      if (jchk .eq. 1) then
         j = myhalf
         if (mf .le. mtrundef(j)) then
            !$acc parallel loop collapse(2) private(sa00, sa10, sa20, sa30) async(async_id)
            do kk = 1, lev2*2, nb
            do ll = mf, jtrun, nb
               do k = kk, min(kk + nb - 1, lev2*2), 4
                  sa00 = tc2(k, 1, 1, j)
                  sa10 = tc2(k + 1, 1, 1, j)
                  sa20 = tc2(k + 2, 1, 1, j)
                  sa30 = tc2(k + 3, 1, 1, j)
                  do l = ll, min(ll + nb - 1, jtrun)
                     sa00 = sa00 + ws3(k, 1, 1, l)*wc(l, j) + ws4(k, 1, 1, l)*wd(l, j)
                     sa10 = sa10 + ws3(k + 1, 1, 1, l)*wc(l, j) + ws4(k + 1, 1, 1, l)*wd(l, j)
                     sa20 = sa20 + ws3(k + 2, 1, 1, l)*wc(l, j) + ws4(k + 2, 1, 1, l)*wd(l, j)
                     sa30 = sa30 + ws3(k + 3, 1, 1, l)*wc(l, j) + ws4(k + 3, 1, 1, l)*wd(l, j)
                  end do
                  tc2(k, 1, 1, j) = sa00
                  tc2(k + 1, 1, 1, j) = sa10
                  tc2(k + 2, 1, 1, j) = sa20
                  tc2(k + 3, 1, 1, j) = sa30
               end do
            end do
            end do
         end if
      end if

      !$acc parallel loop collapse(2) private(jj, jx, j2) async(async_id)
      do j = 1, myhalf
         do k = 1, lev*2*2
            jj = jlist2(j)
            jx = my - j + 1
            j2 = jlist2(jx)
            wcc_fk(k, 1, 1, m, jj) = tcc(k, 1, 1, j)
            wcc_fk(k, 1, 1, m, j2) = tc2(k, 1, 1, j)
         end do
      end do
   end do   ! end of big m loop

   ! Present on device: wcc_fk, twcc_fk
   call mpe_transpose_sr_sp_gpu(wcc_fk, twcc_fk, lev*2*2, jtmax, my_max, &
                                nsize, nccl_col_comm)

   !$acc host_data use_device(cc)
   CUDACHECK(cudaMemSetAsync(cc, real(0.0, RTYPE), size(cc), stream))
   !$acc end host_data

   !$acc parallel loop collapse(3) private(j, jtrunj, mm, mp, mlst) async(async_id)
   do jj = 1, jlistnum
      do m = 1, jtrun
         do k = 1, lev
            j = jlist1(jj)
            jtrunj = mtrundef(j)
            if (m .le. jtrunj) then
               mm = 2*m - 1
               mp = mm + 1
               mlst = nlist(m)
               cc(mm, k, 1, jj) = twcc_fk(k, 1, 1, mlst, jj)
               cc(mp, k, 1, jj) = twcc_fk(k, 2, 1, mlst, jj)
               cc(mm, k, 2, jj) = twcc_fk(k, 1, 2, mlst, jj)
               cc(mp, k, 2, jj) = twcc_fk(k, 2, 2, mlst, jj)
            end if
         end do
      end do
   end do

   if (length_fft .eq. 0 .and. lreduce .eq. 0) then
      call rfftmlt(cc, gwk1, trigs, ifax, 1, nx + 2, nx, lev*jlistnum*2, 1) ! CWB2015
   else
      if (fft_cg%created) then
         CUDACHECK(cudaGraphLaunch(fft_cg%graph_exec, stream))
      else
         call rfftmlt_loop_identical_cuda_graph(cc, gwk1, trigsj, ifaxj, jlist1, nxdef, jlistnum, nx + 2, lev*2, 1, fft_cg%graph)
         CUDACHECK(cudaGraphInstantiate(fft_cg%graph_exec, fft_cg%graph, 0))
         fft_cg%created = .true.
         CUDACHECK(cudaGraphLaunch(fft_cg%graph_exec, stream))
      end if
   end if

   call ujoinsr_gpu(cc, ut, vt, dummy, dummy, nx, my_max, levF, jlistnum, 2, 1)
   !$acc exit data delete(jlist_fj, &
   !$acc& twcc_fk, ws3, ws4, tcc, tc2, wcc_fk, fj_ws3, fj_ws4, fj_tcc, fj_wc, fj_wd, fj_tc2, wc, wd) async(async_id)

   return
end subroutine tranuv_gpu
! ============================================================

subroutine tranuv_gpu_cuda_graph(jtrun, jtmax, nx, my, my_max, lev, coslr, wcfac, &
                                 wdfac, poly, dpoly, vor, div, ut, vt, nsize, &
                                 cc, gwk1, ws1, ws2, wcc, wcc_fk, fj_wp, fj_wd)
! Present on device: coslr, wcfac, wdfac, poly, dpoly, vor, div, ut, vt, cc, gwk1
! jlist1, nlist, mtrundef, mlist, jlist2, nxdef
!  subroutine to transform vorticity and divergence to velocity
!  components
!
! *** input ***
!
!  jtrun: zonal wavenumber truncation limit
!  jtmax: maximum amount of zonal waves located in each pe
!  nx: e-w dimension no.
!  my: n-s dimension no.
!  my_max: maximum amount of n-s grids located in each pe
!  lev: number of veritical levels to transform
!  coslr: cos(lat)**2
!  wcfac: constants defined in cons
!  wdfac: constants defined in cons
!  poly: legendre polynomials
!  dpoly: d(poly)/d(sin(lat))
!  vor: spectral vorticity
!  div: spectral divergence
!
! *** output ***
!
!  ut: e-w velocity component
!  vt: n-s velocity component
!
! *** buffer ***
!
! cc, gwk1, ws1, ws2, wcc, wcc_fk, fj_wp, fj_wd
!
!  ****************************************
!
   use const, only: RTYPE
   use index
   use fftcom
   use spec_cuda_graph, only: fft_cg => tranuv_fft_cg, lt_cg => tranuv_lt_cg
   use openacc
   use cudafor
   use cublas

   implicit none

   ! << output >>
   real(kind=RTYPE), dimension(nxp, levF, my_max):: ut, vt
   ! << input >>
   real(kind=RTYPE), dimension(lev, 2, jtrun, jtmax):: vor, div
   ! << const >>
   real, dimension((jtrun + nsize)*my/2*jtmax) :: poly, dpoly
   real(kind=RTYPE), dimension(jtrun, jtmax):: wcfac, wdfac
   real(kind=RTYPE), dimension(my):: coslr
   ! << buffer >>
   real(kind=RTYPE), dimension(nx + 2, lev, 2, my_max) :: cc, gwk1
   real, dimension((jtrun + nsize)*my/2*jtmax) :: fj_wp, fj_wd
   real, dimension(lev, 2, 2, jtrun, jtmax):: ws1, ws2
   real, dimension(lev, 2, 2, my, jtmax) :: wcc
   real(kind=RTYPE), dimension(lev, 2, 2, jtmax, my_max*nsize):: wcc_fk
   ! << local >>
   real(kind=RTYPE), dimension(lev, 2, 2, jtmax*nsize, my_max):: twcc_fk
   real(kind=RTYPE) dummy

   integer jtrun, jtmax, nx, my, my_max, lev, nsize
   integer myhalf, j, m, mf, l, k, lev4
   integer jlistnum_fj, llistnum_fj, j_str, m_str, ind
   integer ll, jj, jtrunj, mm, mp, mlst

   integer async_id, istat
   integer(kind=cuda_stream_kind) stream, lt_cg_stream(jtmax)
   type(cudaEvent) :: spread_event, pack_event
   type(cublashandle) :: handle

   async_id = 1

   myhalf = my/2
   lev4 = lev*4
   stream = acc_get_cuda_stream(async_id)

   !$acc enter data create(twcc_fk) async(async_id)

   !$acc parallel loop collapse(3) private(mf) async(async_id)
   do m = 1, mlistnum
      do l = 1, jtrun
         do k = 1, lev
            mf = mlist(m)
            if (l .ge. mf) then
               ! -\sqrt{-1} div
               ws1(k, 1, 1, l, m) = +div(k, 2, l, m)
               ws1(k, 1, 2, l, m) = -div(k, 1, l, m)
               ! -\sqrt{-1} vor
               ws1(k, 2, 1, l, m) = +vor(k, 2, l, m)
               ws1(k, 2, 2, l, m) = -vor(k, 1, l, m)

               ! +vor
               ws2(k, 1, 1, l, m) = +vor(k, 1, l, m)
               ws2(k, 1, 2, l, m) = +vor(k, 2, l, m)
               ! -div
               ws2(k, 2, 1, l, m) = -div(k, 1, l, m)
               ws2(k, 2, 2, l, m) = -div(k, 2, l, m)
            end if
         end do
      end do
   end do
   !$acc parallel loop collapse(3) &
   !$acc& private(mf, jlistnum_fj, llistnum_fj, j_str, m_str, ind, ll, jj) async(async_id)
   do m = 1, mlistnum
      do j = 1, myhalf
         do l = 1, jtrun
            mf = mlist(m)
            jlistnum_fj = tcolt_jlist(1, m)
            llistnum_fj = jtrun - mf + 1
            if ((j .le. jlistnum_fj) .AND. (l .le. llistnum_fj)) then
               j_str = tcolt_jlist(2, m)
               m_str = poly_mlist(m)
               ll = mf + l - 1

               ind = l + (j - 1)*llistnum_fj + m_str - 1
               fj_wp(ind) = wcfac(ll, m)*poly(ind)
               fj_wd(ind) = wdfac(ll, m)*dpoly(ind)*coslr(j_str + j - 1)

               ind = ind + jlistnum_fj*llistnum_fj
               jj = myhalf - j + 1
               fj_wp(ind) = wcfac(ll, m)*poly(ind)
               fj_wd(ind) = wdfac(ll, m)*dpoly(ind)*coslr(jj)
            end if
         end do
      end do
   end do

   if (.not. (lt_cg%created)) then
      istat = cudaEventCreate(spread_event)
      istat = cudaEventCreate(pack_event)
      handle = cublasGetHandle()
      do m = 1, mlistnum
         lt_cg_stream(m) = acc_get_cuda_stream(m + 1)
      end do

      CUDACHECK(cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal))

      istat = cudaEventRecord(spread_event, stream)
      !$acc host_data use_device(wcc_fk)
      CUDACHECK(cudaMemsetAsync(wcc_fk, real(0.0, RTYPE), size(wcc_fk), stream))
      !$acc end host_data

      do m = 1, mlistnum
         istat = cudaStreamWaitEvent(lt_cg_stream(m), spread_event, 0)
         mf = mlist(m)

         llistnum_fj = jtrun - mf + 1
         jlistnum_fj = tcolt_jlist(1, m)*2
         m_str = poly_mlist(m)
         !$acc host_data use_device(ws1, ws2, fj_wp, fj_wd, wcc)
         istat = cublasSetStream(handle, lt_cg_stream(m))
         call dgemm('n', 'n', lev4, jlistnum_fj, llistnum_fj, &
                    1.0, ws1(1, 1, 1, mf, m), lev4, &
                    fj_wp(m_str), llistnum_fj, &
                    0.0, wcc(1, 1, 1, 1, m), lev4)

         call dgemm('n', 'n', lev4, jlistnum_fj, llistnum_fj, &
                    1.0, ws2(1, 1, 1, mf, m), lev4, &
                    fj_wd(m_str), llistnum_fj, &
                    1.0, wcc(1, 1, 1, 1, m), lev4)
         istat = cudaEventRecord(pack_event, lt_cg_stream(m))
         istat = cudaStreamWaitEvent(stream, pack_event, 0)
         !$acc end host_data
      end do

      !$acc parallel loop collapse(3) &
      !$acc& private(jj, j_str, jlistnum_fj) async(async_id)
      do m = 1, mlistnum
         do j = 1, my
            do k = 1, lev
               jlistnum_fj = tcolt_jlist(1, m)*2
               if (j .le. jlistnum_fj) then
                  j_str = tcolt_jlist(2, m)
                  jj = jlist2(j_str + j - 1)
                  ! U
                  wcc_fk(k, 1, 1, m, jj) = wcc(k, 1, 1, j, m)
                  wcc_fk(k, 2, 1, m, jj) = wcc(k, 2, 1, j, m)
                  ! V
                  wcc_fk(k, 1, 2, m, jj) = wcc(k, 1, 2, j, m)
                  wcc_fk(k, 2, 2, m, jj) = wcc(k, 2, 2, j, m)
               end if
            end do
         end do
      end do

      CUDACHECK(cudaStreamEndCapture(stream, lt_cg%graph))
      CUDACHECK(cudaGraphInstantiate(lt_cg%graph_exec, lt_cg%graph, 0))
      lt_cg%created = .true.
      istat = cudaEventDestroy(spread_event)
      istat = cudaEventDestroy(pack_event)
   end if
   CUDACHECK(cudaGraphLaunch(lt_cg%graph_exec, stream))

   ! Present on device: wcc_fk, twcc_fk
   call mpe_transpose_sr_sp_gpu(wcc_fk, twcc_fk, lev4, jtmax, my_max, &
                                nsize, nccl_col_comm)

   !$acc host_data use_device(gwk1)
   CUDACHECK(cudaMemSetAsync(gwk1, real(0.0, RTYPE), size(gwk1), stream))
   !$acc end host_data

   !$acc parallel loop collapse(2) private(j, jtrunj, mm, mp, mlst) async(async_id)
   do jj = 1, jlistnum
      do m = 1, jtrun
         j = jlist1(jj)
         jtrunj = mtrundef(j)
         if (m .le. jtrunj) then
            mm = 2*m - 1
            mp = mm + 1
            mlst = nlist(m)
            !$acc loop vector
            do k = 1, lev*2
               gwk1(mm, k, 1, jj) = twcc_fk(k, 1, 1, mlst, jj)
               gwk1(mp, k, 1, jj) = twcc_fk(k, 1, 2, mlst, jj)
            end do
         end if
      end do
   end do

   if (length_fft .eq. 0 .and. lreduce .eq. 0) then
      call rfftmlt(cc, gwk1, trigs, ifax, 1, nx + 2, nx, lev*jlistnum*2, 1) ! CWB2015
   else
      if (.not. (fft_cg%created)) then
         call rfftmlt_loop_identical_cuda_graph(cc, gwk1, trigsj, ifaxj, jlist1, nxdef, jlistnum, nx + 2, lev*2, 1, fft_cg%graph)
         CUDACHECK(cudaGraphInstantiate(fft_cg%graph_exec, fft_cg%graph, 0))
         fft_cg%created = .true.
      end if
      CUDACHECK(cudaGraphLaunch(fft_cg%graph_exec, stream))
   end if

   call ujoinsr_gpu(cc, ut, vt, dummy, dummy, nx, my_max, levF, jlistnum, 2, 1)
   !$acc exit data delete(twcc_fk) async(async_id)

   return
end subroutine tranuv_gpu_cuda_graph
