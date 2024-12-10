#define CUDACHECK(ierr) call cuda_check_helper(ierr, __FILE__, __LINE__)

subroutine rstrandz_gpu(jtrun, jtmax, nx, my, my_max, lev &
                        , vdmer, vdzon, w, cim, onocos, poly, dpoly &
                        , hldten, vorten, nsize, cc, gwk1)
   ! Present on device: vdmer, vdzon, w, cim, onocos, poly, dpoly, hldten, vorten
   ! Present on device: nlist, jlist2, mtrundef, jlist1, nxjlen, nxjlen_all, nxdef
   use const, only: RTYPE
   use index
   use paramt
   use fftcom
   use spec_cuda_graph, only: fft_cg => rstrandz_fft_cg
   use openacc
   use cudafor

   implicit none

   integer jtrun, jtmax, nx, my, my_max, lev, nsize
   integer myhalf, lev2, jj, j, nxj, k, i, m, mm, mp, mlst, mf, j2, j1, l
   integer lchk, lle, jlistnum_fj, j_fj, l_fj, llistnum_fj

   real(kind=RTYPE) poly(jtrun, my/2, jtmax), &
      dpoly(jtrun, my/2, jtmax), cim(jtmax), &
      onocos(my), w(my)

   real(kind=RTYPE) hldten(lev, 2, jtrun, jtmax), vorten(lev, 2, jtrun, jtmax)
   real(kind=RTYPE) vdmer(nxp, levf, my_max), vdzon(nxp, levf, my_max), dummy

   real(kind=RTYPE) gwk1(nx + 2, lev, 2, my_max)
   real(kind=RTYPE) wss(lev, 2, 2, jtrun)
   real(kind=RTYPE) wcc_fk(lev, 2, 2, jtmax, my_max*nsize)
   real(kind=RTYPE) twcc_fk(lev, 2, 2, jtmax*nsize, my_max)
   real(kind=RTYPE) cc(nx + 2, lev, 2, my_max)

   real(kind=RTYPE) wcc2(lev, 2, 2, my/2)
   real(kind=RTYPE) wcc3(lev, 2, 2, my/2)
   real(kind=RTYPE) wcc4(lev, 2, 2, my/2)
   real(kind=RTYPE) wcc5(lev, 2, 2, my/2)

   real(kind=RTYPE) wp(my/2, jtrun), wd(my/2, jtrun)

   real fj_wcc2(lev*2*2, my/2)
   real fj_wcc3(lev*2*2, my/2)
   real fj_wcc4(lev*2*2, my/2)
   real fj_wcc5(lev*2*2, my/2)
   real fj_wd2(my/2, jtrun), fj_wp3(my/2, jtrun)
   real fj_wd4(my/2, jtrun), fj_wp5(my/2, jtrun)
   real fj_wss23(lev*2*2, jtrun)
   real fj_wss45(lev*2*2, jtrun)
   integer jlist_fj(my/2, mlistnum)
   integer jlistnum_fj_array(mlistnum)
   integer async_id, istat
   integer(kind=cuda_stream_kind) :: stream
   real(kind=RTYPE) :: sum_local

   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   !$acc enter data create(twcc_fk, wcc_fk) async(async_id)

   !$acc host_data use_device(twcc_fk, gwk1)
   istat = cudaMemSetAsync(twcc_fk, real(0.0, RTYPE), size(twcc_fk), stream)
   istat = cudaMemSetAsync(gwk1, real(0.0, RTYPE), size(gwk1), stream)
   !$acc end host_data

   myhalf = my/2
   lev2 = lev*2

   ! Present on device: cc, vdmer, vdzon, jlist1, nxjlen, nxjlen_all
   call joinrs_gpu(cc, vdmer, vdzon, dummy, dummy, nx, my_max, levf, jlistnum, 2, 1)

   if (length_fft .eq. 0 .and. lreduce .eq. 0) then
      call rfftmlt(cc, gwk1, trigs, ifax, 1, nx + 2, nx, lev*jlistnum*2, -1)
   else
      if (fft_cg%created) then
         CUDACHECK(cudaGraphLaunch(fft_cg%graph_exec, stream))
      else
         call rfftmlt_loop_identical_cuda_graph(cc, gwk1, trigsj, ifaxj, jlist1, nxdef, jlistnum, nx + 2, lev*2, -1, fft_cg%graph)
         CUDACHECK(cudaGraphInstantiate(fft_cg%graph_exec, fft_cg%graph, 0))
         fft_cg%created = .true.
         CUDACHECK(cudaGraphLaunch(fft_cg%graph_exec, stream))
      end if
   end if

   !$acc parallel loop collapse(3) async(async_id)
   do j = 1, jlistnum
      do m = 1, jtrun
         do k = 1, lev
            mm = 2*m - 1
            mp = mm + 1
            mlst = nlist(m)
            twcc_fk(k, 1, 1, mlst, j) = cc(mm, k, 1, j)
            twcc_fk(k, 2, 1, mlst, j) = cc(mp, k, 1, j)
            twcc_fk(k, 1, 2, mlst, j) = cc(mm, k, 2, j)
            twcc_fk(k, 2, 2, mlst, j) = cc(mp, k, 2, j)
         end do
      end do
   end do

#ifdef SP
   call mpe_transpose_rs_sp_gpu(twcc_fk, wcc_fk, lev*2*2, jtmax, my_max, nsize, nccl_col_comm)
#else
   ! Present on device: twcc_fk, wcc_fk
   call mpe_transpose_rs_gpu(twcc_fk, wcc_fk, lev*2*2, jtmax, my_max, nsize, nccl_col_comm)
#endif
   !$acc exit data delete(twcc_fk) async(async_id)

   do m = 1, mlistnum
      mf = mlist(m)
      if (jtrun .gt. mf) then
         jlistnum_fj = 0
         do j = 1, myhalf
            if (mf .le. mtrundef(j)) then
               jlistnum_fj = jlistnum_fj + 1
               jlist_fj(jlistnum_fj, m) = j
            end if
         end do
         jlistnum_fj_array(m) = jlistnum_fj
      end if
   end do

   !$acc enter data copyin(jlist_fj) create(wcc2, wcc3, wcc4, wcc5, wd, wp, wss, fj_wss23, fj_wss45, fj_wcc2, fj_wcc3, fj_wcc4, fj_wcc5, fj_wd2, fj_wp3, fj_wd4, fj_wp5) async(async_id)
   do m = 1, mlistnum
      mf = mlist(m)

      !$acc parallel loop gang async(async_id)
      do j = 1, myhalf
         j1 = jlist2(j); j2 = jlist2(my - j + 1)
         if (mf .le. mtrundef(j)) then
            !$acc loop vector
            do k = 1, lev
               ! +V
               wcc2(k, 1, 1, j) = +(wcc_fk(k, 1, 1, m, j1) - wcc_fk(k, 1, 1, m, j2))
               wcc2(k, 2, 1, j) = +(wcc_fk(k, 2, 1, m, j1) - wcc_fk(k, 2, 1, m, j2))
               ! -U
               wcc2(k, 1, 2, j) = -(wcc_fk(k, 1, 2, m, j1) - wcc_fk(k, 1, 2, m, j2))
               wcc2(k, 2, 2, j) = -(wcc_fk(k, 2, 2, m, j1) - wcc_fk(k, 2, 2, m, j2))
               ! \sqrt{-i} x U
               wcc3(k, 1, 1, j) = -(wcc_fk(k, 2, 2, m, j1) + wcc_fk(k, 2, 2, m, j2))
               wcc3(k, 2, 1, j) = +(wcc_fk(k, 1, 2, m, j1) + wcc_fk(k, 1, 2, m, j2))
               ! \sqrt{-i} x V
               wcc3(k, 1, 2, j) = -(wcc_fk(k, 2, 1, m, j1) + wcc_fk(k, 2, 1, m, j2))
               wcc3(k, 2, 2, j) = +(wcc_fk(k, 1, 1, m, j1) + wcc_fk(k, 1, 1, m, j2))

               wcc4(k, 1, 1, j) = +(wcc_fk(k, 1, 1, m, j1) + wcc_fk(k, 1, 1, m, j2))
               wcc4(k, 2, 1, j) = +(wcc_fk(k, 2, 1, m, j1) + wcc_fk(k, 2, 1, m, j2))
               wcc4(k, 1, 2, j) = -(wcc_fk(k, 1, 2, m, j1) + wcc_fk(k, 1, 2, m, j2))
               wcc4(k, 2, 2, j) = -(wcc_fk(k, 2, 2, m, j1) + wcc_fk(k, 2, 2, m, j2))

               wcc5(k, 1, 1, j) = -(wcc_fk(k, 2, 2, m, j1) - wcc_fk(k, 2, 2, m, j2))
               wcc5(k, 2, 1, j) = +(wcc_fk(k, 1, 2, m, j1) - wcc_fk(k, 1, 2, m, j2))
               wcc5(k, 1, 2, j) = -(wcc_fk(k, 2, 1, m, j1) - wcc_fk(k, 2, 1, m, j2))
               wcc5(k, 2, 2, j) = +(wcc_fk(k, 1, 1, m, j1) - wcc_fk(k, 1, 1, m, j2))

            end do
         end if
      end do

      !$acc parallel loop collapse(2) async(async_id)
      do l = mf, jtrun
      do j = 1, myhalf
         wd(j, l) = w(j)*dpoly(l, j, m)
         wp(j, l) = w(j)*onocos(j)*cim(m)*poly(l, j, m)
      end do
      end do

      !$acc host_data use_device(wss)
      istat = cudaMemSetAsync(wss, real(0.0, RTYPE), size(wss), stream)
      !$acc end host_data

      lchk = iand(jtrun - mf + 1, 1)
      lle = jtrun - lchk

      if (jtrun .gt. mf) then
         !$acc host_data use_device(fj_wss23, fj_wss45, fj_wcc2, fj_wcc3, fj_wcc4, fj_wcc5, fj_wd2, fj_wp3, fj_wd4, fj_wp5)
         istat = cudaMemSetAsync(fj_wss23, 0.0, size(fj_wss23), stream)
         istat = cudaMemSetAsync(fj_wss45, 0.0, size(fj_wss45), stream)
         istat = cudaMemSetAsync(fj_wcc2, 0.0, size(fj_wcc2), stream)
         istat = cudaMemSetAsync(fj_wcc3, 0.0, size(fj_wcc3), stream)
         istat = cudaMemSetAsync(fj_wcc4, 0.0, size(fj_wcc4), stream)
         istat = cudaMemSetAsync(fj_wcc5, 0.0, size(fj_wcc5), stream)
         istat = cudaMemSetAsync(fj_wd2, 0.0, size(fj_wd2), stream)
         istat = cudaMemSetAsync(fj_wp3, 0.0, size(fj_wp3), stream)
         istat = cudaMemSetAsync(fj_wd4, 0.0, size(fj_wd4), stream)
         istat = cudaMemSetAsync(fj_wp5, 0.0, size(fj_wp5), stream)
         !$acc end host_data

         jlistnum_fj = jlistnum_fj_array(m)

         !$acc parallel loop gang async(async_id)
         do j_fj = 1, jlistnum_fj
            j = jlist_fj(j_fj, m)
            !$acc loop vector
            do k = 1, lev*2*2
               fj_wcc2(k, j_fj) = wcc2(k, 1, 1, j)
               fj_wcc3(k, j_fj) = wcc3(k, 1, 1, j)
               fj_wcc4(k, j_fj) = wcc4(k, 1, 1, j)
               fj_wcc5(k, j_fj) = wcc5(k, 1, 1, j)
            end do
         end do

         !$acc parallel loop gang async(async_id)
         do l = mf, lle, 2
            l_fj = (l - mf + 2)/2
            !$acc loop vector
            do j_fj = 1, jlistnum_fj
               j = jlist_fj(j_fj, m)
               fj_wd2(j_fj, l_fj) = wd(j, l)
               fj_wp3(j_fj, l_fj) = wp(j, l)
               fj_wd4(j_fj, l_fj) = wd(j, l + 1)
               fj_wp5(j_fj, l_fj) = wp(j, l + 1)
            end do
         end do

         do l = mf, lle, 2
            l_fj = (l - mf + 2)/2
         end do

         llistnum_fj = l_fj

         call dgemm_async('n', 'n', lev*2*2, llistnum_fj, jlistnum_fj, 1.0d+0, fj_wcc2, lev*2*2, fj_wd2, myhalf, &
                          1.0d+0, fj_wss23, lev*2*2, async_id)

         call dgemm_async('n', 'n', lev*2*2, llistnum_fj, jlistnum_fj, 1.0d+0, fj_wcc3, lev*2*2, fj_wp3, myhalf, &
                          1.0d+0, fj_wss23, lev*2*2, async_id)

         !$acc parallel loop gang async(async_id)
         do L = mf, lle, 2
            l_fj = (l - mf + 2)/2
            sum_local = 0.0
            !$acc loop vector
            do k = 1, lev*2*2
               wss(k, 1, 1, L) = wss(k, 1, 1, L) + fj_wss23(k, l_fj)
            end do
         end do

         call dgemm_async('n', 'n', lev*2*2, llistnum_fj, jlistnum_fj, 1.0d+0, fj_wcc4, lev*2*2, fj_wd4, myhalf, &
                          1.0d+0, fj_wss45, lev*2*2, async_id)

         call dgemm_async('n', 'n', lev*2*2, llistnum_fj, jlistnum_fj, 1.0d+0, fj_wcc5, lev*2*2, fj_wp5, myhalf, &
                          1.0d+0, fj_wss45, lev*2*2, async_id)

         !$acc parallel loop gang async(async_id)
         do L = mf, lle, 2
            l_fj = (l - mf + 2)/2
            !$acc loop vector
            do k = 1, lev*2*2
               wss(k, 1, 1, L + 1) = wss(k, 1, 1, L + 1) + fj_wss45(k, l_fj)
            end do
         end do

      end if

      if (lchk .eq. 1) then
         L = jtrun
         !$acc parallel loop async(async_id)
         do k = 1, lev*2*2
            !$acc loop seq
            do j = 1, myhalf
               if (mf .le. mtrundef(j)) then
                  wss(k, 1, 1, L) = wss(k, 1, 1, L) + wcc2(k, 1, 1, j)*wd(j, L) + wcc3(k, 1, 1, j)*wp(j, L)
               end if
            end do
         end do
      end if

      !$acc parallel loop collapse(2) async(async_id)
      do L = mf, jtrun
      do k = 1, lev*2
         hldten(k, 1, L, m) = wss(k, 1, 1, L)
         vorten(k, 1, L, m) = -wss(k, 1, 2, L)
      end do
      end do

   end do
   !$acc exit data delete(wcc_fk, wcc2, wcc3, wcc4, wcc5, wd, wp, wss, jlist_fj, fj_wss23, fj_wss45, fj_wcc2, fj_wcc3, fj_wcc4, fj_wcc5, fj_wd2, fj_wp3, fj_wd4, fj_wp5) async(async_id)
   return
end subroutine rstrandz_gpu
! ============================================================
subroutine rstrandz_gpu_cuda_graph(jtrun, jtmax, nx, my, my_max, lev, &
                                   vdmer, vdzon, w, cim, onocos, poly, dpoly, &
                                   hldten, vorten, nsize, &
                                   cc, gwk1, wss, wc1, wc2, wcc_fk, fj_wp, fj_wd)
   ! Present on device: vdmer, vdzon, w, cim, onocos, poly, dpoly, hldten, vorten
   ! Present on device: nlist, jlist2, mtrundef, jlist1, nxjlen, nxjlen_all, nxdef
   ! Presnet on device: cc, gwk1, wss, wc1, wc2, wcc_fk, fj_wp, fj_wd
   ! Presnet on device: tcolt_jlist, poly_mlist
   use const, only: RTYPE
   use index
   use paramt
   use fftcom
   use spec_cuda_graph, only: fft_cg => rstrandz_fft_cg, lt_cg => rstrandz_lt_cg
   use openacc
   use cudafor
   use cublas

   implicit none

   integer jtrun, jtmax, nx, my, my_max, lev, nsize
   ! << output >>
   real(kind=RTYPE), dimension(lev, 2, jtrun, jtmax):: hldten, vorten
   ! << input >>
   real(kind=RTYPE), dimension(nxp, levf, my_max)::vdmer, vdzon
   ! << const >>
   real(kind=RTYPE), dimension(jtmax):: cim
   real(kind=RTYPE), dimension(my):: onocos, w
   real, dimension((jtrun + nsize)*my/2*jtmax):: poly, dpoly
   ! << buffer >>
   real(kind=RTYPE), dimension(nx + 2, lev, 2, my_max):: cc, gwk1
   real, dimension(lev, 2, 2, jtrun, jtmax):: wss
   real, dimension(lev, 2, 2, my, jtmax):: wc1, wc2
   real(kind=RTYPE), dimension(lev, 2, 2, jtmax, my_max*nsize):: wcc_fk
   real, dimension((jtrun + nsize)*my/2*jtmax):: fj_wd, fj_wp
   ! << local >>
   real(kind=RTYPE), dimension(lev, 2, 2, jtmax*nsize, my_max):: twcc_fk
   real(kind=RTYPE) dummy

   integer myhalf, lev2, jj, j, k, i, m, mm, mp, mlst, mf, j2, j1, l
   integer jlistnum_fj, llistnum_fj, jtrunj, j_str, m_str, ind

   integer async_id, istat
   integer(kind=cuda_stream_kind):: stream, lt_cg_stream(jtmax)
   type(cudaEvent):: spread_event, pack_event
   type(cublashandle):: handle

   async_id = 1
   stream = acc_get_cuda_stream(async_id)
   myhalf = my/2
   lev2 = lev*2

   !$acc enter data create(twcc_fk) async(async_id)

   !$acc host_data use_device(twcc_fk, gwk1)
   istat = cudaMemSetAsync(twcc_fk, real(0.0, RTYPE), size(twcc_fk), stream)
   istat = cudaMemSetAsync(gwk1, real(0.0, RTYPE), size(gwk1), stream)
   !$acc end host_data

   ! Present on device: cc, vdmer, vdzon, jlist1, nxjlen, nxjlen_all
   call joinrs_gpu(cc, vdmer, vdzon, dummy, dummy, nx, my_max, levf, jlistnum, 2, 1)

   if (length_fft .eq. 0 .and. lreduce .eq. 0) then
      call rfftmlt(cc, gwk1, trigs, ifax, 1, nx + 2, nx, lev*jlistnum*2, -1)
   else
      if (fft_cg%created) then
         CUDACHECK(cudaGraphLaunch(fft_cg%graph_exec, stream))
      else
         call rfftmlt_loop_identical_cuda_graph(cc, gwk1, trigsj, ifaxj, jlist1, nxdef, jlistnum, nx + 2, lev*2, -1, fft_cg%graph)
         CUDACHECK(cudaGraphInstantiate(fft_cg%graph_exec, fft_cg%graph, 0))
         fft_cg%created = .true.
         CUDACHECK(cudaGraphLaunch(fft_cg%graph_exec, stream))
      end if
   end if

   !$acc parallel loop collapse(3) private(mm, mp, mlst) async(async_id)
   do j = 1, jlistnum
      do m = 1, jtrun
         ! jj = jlist1(j)
         ! jtrunj = mtrundef(jj)
         ! if (m .le. jtrunj) then
         do k = 1, lev
            mm = 2*m - 1
            mp = mm + 1
            mlst = nlist(m)
            ! V
            twcc_fk(k, 1, 1, mlst, j) = cc(mm, k, 1, j)
            twcc_fk(k, 1, 2, mlst, j) = cc(mp, k, 1, j)
            ! U
            twcc_fk(k, 2, 1, mlst, j) = cc(mm, k, 2, j)
            twcc_fk(k, 2, 2, mlst, j) = cc(mp, k, 2, j)
         end do
         ! end if
      end do
   end do

#ifdef SP
   call mpe_transpose_rs_sp_gpu(twcc_fk, wcc_fk, lev*2*2, jtmax, my_max, nsize, nccl_col_comm)
#else
   ! Present on device: twcc_fk, wcc_fk
   call mpe_transpose_rs_gpu(twcc_fk, wcc_fk, lev*2*2, jtmax, my_max, nsize, nccl_col_comm)
#endif
   !$acc exit data delete(twcc_fk) async(async_id)

   !$acc parallel loop collapse(3) private(jlistnum_fj, j1, j2) async(async_id)
   do j = 1, my
      do m = 1, mlistnum
         do k = 1, lev
            jlistnum_fj = tcolt_jlist(1, m)*2
            if (j .le. jlistnum_fj) then
               j1 = tcolt_jlist(2, m) + j - 1
               j2 = jlist2(j1)
               ! +V
               wc1(k, 1, 1, j, m) = wcc_fk(k, 1, 1, m, j2)
               wc1(k, 1, 2, j, m) = wcc_fk(k, 1, 2, m, j2)
               ! +U
               wc1(k, 2, 1, j, m) = wcc_fk(k, 2, 1, m, j2)
               wc1(k, 2, 2, j, m) = wcc_fk(k, 2, 2, m, j2)
               ! \sqrt{-1} x U
               wc2(k, 1, 1, j, m) = -wcc_fk(k, 2, 2, m, j2)
               wc2(k, 1, 2, j, m) = +wcc_fk(k, 2, 1, m, j2)
               ! -\sqrt{-1} x V
               wc2(k, 2, 1, j, m) = +wcc_fk(k, 1, 2, m, j2)
               wc2(k, 2, 2, j, m) = -wcc_fk(k, 1, 1, m, j2)
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

      !$acc parallel loop collapse(3) &
      !$acc& private(mf, jlistnum_fj, llistnum_fj, j_str, m_str, ind, jj) async(async_id)
      do m = 1, mlistnum
         do j = 1, myhalf
            do l = 1, jtrun
               mf = mlist(m)
               jlistnum_fj = tcolt_jlist(1, m)
               llistnum_fj = jtrun - mf + 1
               if ((j .le. jlistnum_fj) .AND. (l .le. llistnum_fj)) then
                  j_str = tcolt_jlist(2, m)
                  m_str = poly_mlist(m)

                  ind = l + (j - 1)*llistnum_fj + m_str - 1
                  jj = j_str + j - 1
                  fj_wd(ind) = w(jj)*dpoly(ind)
                  fj_wp(ind) = w(jj)*onocos(jj)*cim(m)*poly(ind)

                  ind = ind + jlistnum_fj*llistnum_fj
                  jj = myhalf - j + 1
                  fj_wd(ind) = w(jj)*dpoly(ind)
                  fj_wp(ind) = w(jj)*onocos(jj)*cim(m)*poly(ind)
               end if
            end do
         end do
      end do

      istat = cudaEventRecord(spread_event, stream)
      do m = 1, mlistnum
         istat = cudaStreamWaitEvent(lt_cg_stream(m), spread_event, 0)

         mf = mlist(m)
         llistnum_fj = jtrun - mf + 1
         jlistnum_fj = tcolt_jlist(1, m)*2
         m_str = poly_mlist(m)

         !$acc host_data use_device(fj_wp, fj_wd, wc1, wc2, wss)
         istat = cublasSetStream(handle, lt_cg_stream(m))
         call dgemm('n', 't', lev*4, llistnum_fj, jlistnum_fj, &
                    1.0, wc1(1, 1, 1, 1, m), lev*4, &
                    fj_wd(m_str), llistnum_fj, &
                    0.0, wss(1, 1, 1, mf, m), lev*4)

         call dgemm('n', 't', lev*4, llistnum_fj, jlistnum_fj, &
                    1.0, wc2(1, 1, 1, 1, m), lev*4, &
                    fj_wp(m_str), llistnum_fj, &
                    1.0, wss(1, 1, 1, mf, m), lev*4)
         !$acc end host_data
         istat = cudaEventRecord(pack_event, lt_cg_stream(m))
         istat = cudaStreamWaitEvent(stream, pack_event, 0)
      end do

      CUDACHECK(cudaStreamEndCapture(stream, lt_cg%graph))
      CUDACHECK(cudaGraphInstantiate(lt_cg%graph_exec, lt_cg%graph, 0))
      lt_cg%created = .true.
      istat = cudaEventDestroy(spread_event)
      istat = cudaEventDestroy(pack_event)
   end if
   CUDACHECK(cudaGraphLaunch(lt_cg%graph_exec, stream))

   !$acc parallel loop collapse(3) private(mf) async(async_id)
   do m = 1, mlistnum
      do L = 1, jtrun
         do k = 1, lev
            mf = mlist(m)
            if (L .ge. mf) then
               hldten(k, 1, L, m) = wss(k, 1, 1, L, m)
               hldten(k, 2, L, m) = wss(k, 1, 2, L, m)

               vorten(k, 1, L, m) = wss(k, 2, 1, L, m)
               vorten(k, 2, L, m) = wss(k, 2, 2, L, m)
            end if
         end do
      end do
   end do

   return
end subroutine rstrandz_gpu_cuda_graph
