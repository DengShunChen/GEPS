!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define CUDACHECK(ierr) call cuda_check_helper(ierr, __FILE__, __LINE__)

subroutine transr_gpu(jtrun, jtmax, nx, my, my_max, lev, poly, wss &
                      , cc, num, nsize, gwk1)
! Present on device: poly, wss, cc, jlist2, jlist1, mtrundef, nlist, mlist, nxdef
!
!  subroutine to transform a spectral coefficient field to
!  grid point form
!
! *** input ***
!
!  jtrun: zonal wavenumber resolution limit
!  jtmax: maximum amount of zonal waves located in each pe
!  nx: e-w dimension no.
!  my: n-s dimension no.
!  my_max: maximum amount of n-s grids located in each pe
!  lev: number of levels to transform
!  poly: legendre polynomials
!  wss: spectral coefficient array to transform
!  num: number of variables grouped together
!
! *** output ***
!
!  cc_r8: 3-d output grid point fields
!
!  **************************************
!
   use const, only: RTYPE
   use index
   use fftcom
   use spec_cuda_graph, only: fft_cg => transr_fft_cg
   use openacc
   use cudafor

   implicit none

   integer jtrun, jtmax, nx, my, my_max, lev, num, nsize
   integer mlx, myhalf, lev2, m, mf, lmax, lchk, j, k, l
   integer nb, jchk, jje, jlistnum_fj, j_fj, l_fj, llistnum_fj
   integer kk, ll, jj, jx, j2, ii, i, jtrunj, mchk, mm, mp, mlst
   integer mm1, mp1, mlst1, mm2, mp2, mlst2, mm3, mp3, mlst3, nxj

   real(kind=RTYPE) sa00, sa10, sb00, sb10

   real(kind=RTYPE) poly(jtrun, my/2, jtmax)
   real(kind=RTYPE) cc(nx + 2, lev, num, my_max), wss(lev, 2, num, jtrun, jtmax)

   real(kind=RTYPE) gwk1(nx + 2, lev, num, my_max)

   real(kind=RTYPE) wcc_fk(lev, 2, num, jtmax, my_max*nsize)
   real(kind=RTYPE) twcc_fk(lev, 2, num, jtmax*nsize, my_max)
   real(kind=RTYPE) tcc(lev, 2, num, my/2), tc2(lev, 2, num, my/2)
   real ws2(lev, 2, num, jtrun)

   real fj_poly(jtrun, my/2)
   real fj_tcc(lev*2*num, my/2)
   real fj_tc2(lev*2*num, my/2)
   real fj_ws2(lev*2*num, jtrun)
   real fj_wss(lev*2*num, jtrun)
   integer jlist_fj(my/2, mlistnum)
   integer async_id
   integer(kind=cuda_stream_kind) stream
   integer jlistnum_fj_array(mlistnum)

   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   !$acc enter data create(ws2, tcc, tc2, fj_tcc, fj_tc2, fj_wss, fj_ws2, fj_poly, wcc_fk, twcc_fk, jlist_fj, jlistnum_fj_array) async(async_id)
   !$acc host_data use_device(cc, wcc_fk)
   CUDACHECK(cudaMemsetAsync(cc, real(0.0, RTYPE), size(cc), stream))
   CUDACHECK(cudaMemsetAsync(wcc_fk, real(0.0, RTYPE), size(wcc_fk), stream))
   !$acc end host_data

   mlx = (jtrun/2)*((jtrun + 1)/2)
   myhalf = my/2
   lev2 = lev*2
   nb = 32
   jchk = iand(myhalf, 1)
   jje = myhalf - jchk

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

   !$acc exit data copyout(jlistnum_fj_array) async(async_id)
   !$acc wait(async_id)

   do m = 1, mlistnum
      mf = mlist(m)
      lmax = jtrun - mf + 1
      lchk = iand(lmax, 1)
      jlistnum_fj = jlistnum_fj_array(m)

      !$acc parallel loop collapse(2) async(async_id)
      do l = mf, jtrun - 1, 2
      do k = 1, lev2*num
         ws2(k, 1, 1, l) = wss(k, 1, 1, l, m)
         ws2(k, 1, 1, l + 1) = -wss(k, 1, 1, l + 1, m)
      end do
      end do
      if (lchk .eq. 1) then
         l = jtrun
         !$acc parallel loop async(async_id)
         do k = 1, lev2*num
            ws2(k, 1, 1, l) = wss(k, 1, 1, l, m)
         end do
      end if
      !$acc kernels async(async_id)
      !$acc loop
      do k = 1, lev2*num*myhalf
         tcc(k, 1, 1, 1) = 0.
         tc2(k, 1, 1, 1) = 0.
      end do

      fj_tcc = 0.0
      fj_tc2 = 0.0
      fj_wss = 0.0
      fj_ws2 = 0.0
      fj_poly = 0.0

      !$acc loop gang
      do j_fj = 1, jlistnum_fj
         j = jlist_fj(j_fj, m)
         !$acc loop vector
         do l = mf, jtrun
            l_fj = l - mf + 1
            fj_poly(l_fj, j_fj) = poly(l, j, m)
         end do
      end do
      !$acc loop gang
      do l = mf, jtrun
         l_fj = l - mf + 1
         !$acc loop vector
         do k = 1, lev2*num
            fj_wss(k, l_fj) = wss(k, 1, 1, l, m)
         end do
      end do
      !$acc loop gang
      do l = mf, jtrun
         l_fj = l - mf + 1
         !$acc loop vector
         do k = 1, lev2*num
            fj_ws2(k, l_fj) = ws2(k, 1, 1, l)
         end do
      end do
      !$acc end kernels

      llistnum_fj = jtrun - mf + 1

      call dgemm_async('n', 'n', lev2*num, jlistnum_fj, llistnum_fj &
                       , 1.0d+0, fj_wss, lev2*num, fj_poly, jtrun, 1.0d+0, fj_tcc, lev2*num, async_id)

      !$acc parallel loop gang async(async_id)
      do j_fj = 1, jlistnum_fj
         j = jlist_fj(j_fj, m)
         !$acc loop vector
         do k = 1, lev2*num
            tcc(k, 1, 1, j) = tcc(k, 1, 1, j) + fj_tcc(k, j_fj)
         end do
      end do

      call dgemm_async('n', 'n', lev2*num, jlistnum_fj, llistnum_fj &
                       , 1.0d+0, fj_ws2, lev2*num, fj_poly, jtrun, 1.0d+0, fj_tc2, lev2*num, async_id)

      !$acc parallel loop gang async(async_id)
      do j_fj = 1, jlistnum_fj
         j = jlist_fj(j_fj, m)
         !$acc loop vector
         do k = 1, lev2*num
            tc2(k, 1, 1, j) = tc2(k, 1, 1, j) + fj_tc2(k, j_fj)
         end do
      end do

!
! odd number
!
      if (jchk .eq. 1) then
         j = myhalf
         if (mf .le. mtrundef(j)) then
            !$acc parallel loop gang collapse(2) async(async_id)
            do kk = 1, lev2*num, nb
            do ll = mf, jtrun, nb
               !$acc loop vector
               do k = kk, min(kk + nb - 1, lev2*num), 2
                  sa00 = tcc(k, 1, 1, j)
                  sa10 = tcc(k + 1, 1, 1, j)
                  sb00 = tc2(k, 1, 1, j)
                  sb10 = tc2(k + 1, 1, 1, j)
                  do l = ll, min(ll + nb - 1, jtrun)
                     sa00 = sa00 + poly(l, j, m)*wss(k, 1, 1, l, m)
                     sa10 = sa10 + poly(l, j, m)*wss(k + 1, 1, 1, l, m)
                     sb00 = sb00 + poly(l, j, m)*ws2(k, 1, 1, l)
                     sb10 = sb10 + poly(l, j, m)*ws2(k + 1, 1, 1, l)
                  end do
                  tcc(k, 1, 1, j) = sa00
                  tcc(k + 1, 1, 1, j) = sa10
                  tc2(k, 1, 1, j) = sb00
                  tc2(k + 1, 1, 1, j) = sb10
               end do
            end do
            end do
         end if
      end if

      !$acc parallel loop gang async(async_id)
      do j = 1, myhalf
         jj = jlist2(j)
         jx = my - j + 1
         j2 = jlist2(jx)
         !$acc loop vector
         do k = 1, lev2*num
            wcc_fk(k, 1, 1, m, jj) = tcc(k, 1, 1, j)
            wcc_fk(k, 1, 1, m, j2) = tc2(k, 1, 1, j)
         end do
      end do
   end do

   call mpe_transpose_sr_sp_gpu(wcc_fk, twcc_fk, lev*2*num, jtmax, my_max, nsize, nccl_col_comm)

   !$acc parallel loop gang async(async_id)
   do jj = 1, jlistnum

      !$acc loop vector collapse(3)
      do ii = 1, num
      do k = 1, lev
      do i = 1, nx + 2
         cc(i, k, ii, jj) = 0.
      end do
      end do
      end do

      j = jlist1(jj)
      jtrunj = mtrundef(j)
      mchk = iand(jtrunj, 3)

      !$acc loop worker
      do m = 1, mchk
         mm = 2*m - 1
         mp = mm + 1
         mlst = nlist(m)
         !$acc loop vector collapse(2)
         do ii = 1, num
         do k = 1, lev
            cc(mm, k, ii, jj) = twcc_fk(k, 1, ii, mlst, jj)
            cc(mp, k, ii, jj) = twcc_fk(k, 2, ii, mlst, jj)
         end do
         end do
      end do
      !$acc loop worker
      do m = mchk + 1, jtrunj, 4
         mm = 2*m - 1
         mp = mm + 1
         mlst = nlist(m)
         mm1 = 2*(m + 1) - 1
         mp1 = mm1 + 1
         mlst1 = nlist(m + 1)
         mm2 = 2*(m + 2) - 1
         mp2 = mm2 + 1
         mlst2 = nlist(m + 2)
         mm3 = 2*(m + 3) - 1
         mp3 = mm3 + 1
         mlst3 = nlist(m + 3)
         !$acc loop vector collapse(2)
         do ii = 1, num
         do k = 1, lev
            cc(mm, k, ii, jj) = twcc_fk(k, 1, ii, mlst, jj)
            cc(mp, k, ii, jj) = twcc_fk(k, 2, ii, mlst, jj)
            cc(mm1, k, ii, jj) = twcc_fk(k, 1, ii, mlst1, jj)
            cc(mp1, k, ii, jj) = twcc_fk(k, 2, ii, mlst1, jj)
            cc(mm2, k, ii, jj) = twcc_fk(k, 1, ii, mlst2, jj)
            cc(mp2, k, ii, jj) = twcc_fk(k, 2, ii, mlst2, jj)
            cc(mm3, k, ii, jj) = twcc_fk(k, 1, ii, mlst3, jj)
            cc(mp3, k, ii, jj) = twcc_fk(k, 2, ii, mlst3, jj)
         end do
         end do
      end do

   end do

   if (length_fft .eq. 0 .and. lreduce .eq. 0) then
      call rfftmlt(cc, gwk1, trigs, ifax, 1, nx + 2, nx, lev*jlistnum*num, 1)
   else
      if (fft_cg%created) then
         CUDACHECK(cudaGraphLaunch(fft_cg%graph_exec, stream))
      else
         call rfftmlt_loop_identical_cuda_graph(cc, gwk1, trigsj, ifaxj, jlist1, nxdef, jlistnum, nx + 2, lev*num, 1, fft_cg%graph)
         CUDACHECK(cudaGraphInstantiate(fft_cg%graph_exec, fft_cg%graph, 0))
         fft_cg%created = .true.
         CUDACHECK(cudaGraphLaunch(fft_cg%graph_exec, stream))
      end if
   end if

   !$acc exit data delete(jlist_fj, ws2, tcc, tc2, fj_tcc, fj_tc2, fj_wss, fj_ws2, fj_poly, wcc_fk, twcc_fk) async(async_id)

   return
end subroutine transr_gpu
! ------------------------------------------------------------

subroutine transr_gpu_cuda_graph(jtrun, jtmax, nx, my, my_max, lev, poly, s &
                                 , cc, num, nsize, gwk1, wss, tcc, wcc_fk)
! Present on device: poly, wss, cc, jlist2, jlist1, mtrundef, nlist, mlist, nxdef
!
!  subroutine to transform a spectral coefficient field to
!  grid point form
!
! *** input ***
!
!  jtrun: zonal wavenumber resolution limit
!  jtmax: maximum amount of zonal waves located in each pe
!  nx: e-w dimension no.
!  my: n-s dimension no.
!  my_max: maximum amount of n-s grids located in each pe
!  lev: number of levels to transform
!  poly: legendre polynomials
!  wss: spectral coefficient array to transform
!  num: number of variables grouped together
!
! *** output ***
!
!  cc_r8: 3-d output grid point fields
!
!  **************************************
!
   use const, only: RTYPE
   use index
   use fftcom
   use spec_cuda_graph, only: fft_cg => transr_fft_cg, lt_cg => transr_lt_cg
   use openacc
   use cudafor
   use cublas

   implicit none

   ! << output >>
   real(kind=RTYPE), dimension(nx + 2, lev, num, my_max):: cc
   ! << input >>
   real(kind=RTYPE), dimension(lev, 2, num, jtrun, jtmax) :: s
   ! << const >>
   real, dimension((jtrun + nsize)*my/2*jtmax) :: poly
   ! << buffer >>
   real(kind=RTYPE), dimension(nx + 2, lev, num, my_max):: gwk1
   real, dimension(lev, 2, num, jtrun, jtmax) :: wss
   real, dimension(lev, 2, num, my, jtmax) :: tcc
   real(kind=RTYPE), dimension(lev, 2, num, jtmax, my_max*nsize):: wcc_fk
   ! << local >>
   real(kind=RTYPE), dimension(lev, 2, num, jtmax*nsize, my_max):: twcc_fk

   integer jtrun, jtmax, nx, my, my_max, lev, num, nsize
   integer myhalf, nlev2, m, mf, j, k, l
   integer jlistnum_fj, llistnum_fj, j_str, m_str
   integer j, jj, ii, i, jtrunj, mm, mp, mlst

   integer :: async_id, istat, id_str
   integer(kind=cuda_stream_kind) :: stream, lt_cg_stream(jtmax)
   type(cudaEvent) :: spread_event, pack_event
   type(cublashandle) :: handle

   async_id = 1
   stream = acc_get_cuda_stream(async_id)
   nlev2 = lev*2*num
   !$acc enter data create(twcc_fk) async(async_id)
   !$acc host_data use_device(gwk1)
   CUDACHECK(cudaMemsetAsync(gwk1, real(0.0, RTYPE), size(gwk1), stream))
   !$acc end host_data

   !$acc parallel loop collapse(3) private(mf) async(async_id)
   do m = 1, mlistnum
      do L = 1, jtrun
         do k = 1, nlev2
            mf = mlist(m)
            if (L .ge. mf) then
               wss(k, 1, 1, L, m) = s(k, 1, 1, L, m)
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
      CUDACHECK(cudaMemSetAsync(wcc_fk, real(0.0, RTYPE), size(wcc_fk), stream))
      !$acc end host_data

      do m = 1, mlistnum
         istat = cudaStreamWaitEvent(lt_cg_stream(m), spread_event, 0)

         mf = mlist(m)
         jlistnum_fj = tcolt_jlist(1, m)*2
         llistnum_fj = jtrun - mf + 1
         m_str = poly_mlist(m)
         !$acc host_data use_device(wss, tcc, poly)
         ! --------------------
         istat = cublasSetStream(handle, lt_cg_stream(m))
         call dgemm('n', 'n', nlev2, jlistnum_fj, llistnum_fj, &
                    1.0, wss(1, 1, 1, mf, m), nlev2, &
                    poly(m_str), llistnum_fj, &
                    0.0, tcc(1, 1, 1, 1, m), nlev2)
         istat = cudaEventRecord(pack_event, lt_cg_stream(m))
         istat = cudaStreamWaitEvent(stream, pack_event, 0)
         !$acc end host_data
      end do

      !$acc parallel loop collapse(3) private(jj, jlistnum_fj, j_str) async(async_id)
      do m = 1, mlistnum
         do j = 1, my
            do k = 1, nlev2
               jlistnum_fj = tcolt_jlist(1, m)*2
               if (j .le. jlistnum_fj) then
                  j_str = tcolt_jlist(2, m)
                  jj = jlist2(j_str + j - 1)

                  wcc_fk(k, 1, 1, m, jj) = tcc(k, 1, 1, j, m)
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

   call mpe_transpose_sr_sp_gpu(wcc_fk, twcc_fk, nlev2, jtmax, my_max, &
                                nsize, nccl_col_comm)

   !$acc parallel loop collapse(3) private(j, jtrunj, mm, mp, mlst) async(async_id)
   do jj = 1, jlistnum
      do m = 1, jtrun
         do ii = 1, num
            j = jlist1(jj)
            jtrunj = mtrundef(j)
            if (m .le. jtrunj) then
               mm = 2*m - 1
               mp = mm + 1
               mlst = nlist(m)
               !$acc loop vector
               do k = 1, lev
                  gwk1(mm, k, ii, jj) = twcc_fk(k, 1, ii, mlst, jj)
                  gwk1(mp, k, ii, jj) = twcc_fk(k, 2, ii, mlst, jj)
               end do
            end if
         end do
      end do
   end do

   if (length_fft .eq. 0 .and. lreduce .eq. 0) then
      call rfftmlt(cc, gwk1, trigs, ifax, 1, nx + 2, nx, lev*jlistnum*num, 1)
   else
      if (.not. (fft_cg%created)) then
         call rfftmlt_loop_identical_cuda_graph(cc, gwk1, trigsj, ifaxj, jlist1, nxdef, jlistnum, nx + 2, lev*num, 1, fft_cg%graph)
         CUDACHECK(cudaGraphInstantiate(fft_cg%graph_exec, fft_cg%graph, 0))
         fft_cg%created = .true.
      end if
      CUDACHECK(cudaGraphLaunch(fft_cg%graph_exec, stream))
   end if

   !$acc exit data delete(twcc_fk) async(async_id)

   return
end subroutine transr_gpu_cuda_graph
