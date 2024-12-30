#define CUDACHECK(ierr) call cuda_check_helper(ierr, __FILE__, __LINE__)

subroutine trngra3_gpu(jtrun, jtmax, nx, lev, my, my_max, cim, poly, dpoly, s, dlpl, dtpl, nsize, cc, gwk1)
! Present on device: cim, poly, dpoly, s, dlpl, dtpl, mtrundef, jlist1, jlist2, nlist, mlist, nxdef, cc, gwk1
!  subroutine to transform spectral terrain pressure to grid point
!  fields of zonal and meridional derivatives of terrain pressure
!
! *** input ***
!
!  jtrun: zonal wavenumber resolution limit
!  jtmax: maximum amount of zonal waves located in each pe
!  nx: e-w dimension no.
!  my: n-s dimension no.
!  my_max: maximum amount of n-s grids located in each pe
!  cim: zonal wavenumber array
!  poly: legendre polynomials
!  dpoly: d(poly)/d(sin(lat))
!  s: spectral coefficient array
!  nsize: pe number
!
! **** output ****
!
!  dlpl: d(pt)/d(longitude)
!  dtpl: d(pt)/d(sin(lat))
!
! ****************************************************
!
   use const, only: RTYPE
   use index
   use fftcom
   use spec_cuda_graph, only: fft_cg => trngra3_fft_cg
   use openacc
   use cudafor

   implicit none

   integer :: jtrun, jtmax, nx, lev, my, my_max, nsize
   real(kind=RTYPE), dimension(jtrun, my/2, jtmax) :: poly, dpoly
   real(kind=RTYPE) :: p_t, dp_t, s1_t, s2_t
   real(kind=RTYPE), dimension(lev, 2, jtrun, jtmax) :: s
   real(kind=RTYPE), dimension(jtmax) :: cim
   real(kind=RTYPE), dimension(nxp, levF, my_max) :: dlpl, dtpl
   integer :: myhalf, k, m, mf, l, j, jj, i, jtrunj, mm, mp, mlst, nxj, ii
   real(kind=RTYPE), dimension(nx + 2, lev, 2, my_max) :: cc, gwk1
   ! real(kind=RTYPE), dimension(my_max, jtmax*nsize, lev, 2) :: twcc_fk, twdd_fk
   ! real(kind=RTYPE), dimension(my_max*nsize, jtmax, lev, 2) :: wcu_fk, wcv_fk
   real, dimension(lev, 2, 2, jtmax, my_max*nsize) :: wcc_fk
   real, dimension(lev, 2, 2, jtmax*nsize, my_max) :: twcc_fk
   real(kind=RTYPE), dimension(2, 2) :: wcu_fk_t, wcv_fk_t
   real(kind=RTYPE) :: dummy
   integer :: async_id, istat
   integer(kind=cuda_stream_kind) :: stream

   async_id = 1
   stream = acc_get_cuda_stream(async_id)
   myhalf = my/2

   !$acc enter data create(wcc_fk, twcc_fk) async(async_id)

   !$acc host_data use_device(wcc_fk)
   CUDACHECK(cudaMemSetAsync(wcc_fk, 0.0, size(wcc_fk), stream))
   !$acc end host_data

   !$acc parallel loop collapse(3) private(mf, jj, wcu_fk_t, wcv_fk_t, s1_t, s2_t, p_t, dp_t) async(async_id)
   do m = 1, mlistnum
      do j = 1, myhalf
         do k = 1, lev
            mf = mlist(m)
            if (mf .le. mtrundef(j)) then
               wcu_fk_t = 0.0
               wcv_fk_t = 0.0
               !$acc loop seq
               do L = mf, jtrun, 2
                  s1_t = s(k, 1, L, m)
                  s2_t = s(k, 2, L, m)
                  p_t = poly(L, j, m)
                  dp_t = dpoly(L, j, m)
                  wcu_fk_t(1, 1) = wcu_fk_t(1, 1) + s2_t*p_t
                  wcu_fk_t(2, 1) = wcu_fk_t(2, 1) - s1_t*p_t
                  wcv_fk_t(1, 1) = wcv_fk_t(1, 1) - s1_t*dp_t
                  wcv_fk_t(2, 1) = wcv_fk_t(2, 1) - s2_t*dp_t
                  wcu_fk_t(1, 2) = wcu_fk_t(1, 2) + s2_t*p_t
                  wcu_fk_t(2, 2) = wcu_fk_t(2, 2) - s1_t*p_t
                  wcv_fk_t(1, 2) = wcv_fk_t(1, 2) + s1_t*dp_t
                  wcv_fk_t(2, 2) = wcv_fk_t(2, 2) + s2_t*dp_t
               end do
               !$acc loop seq
               do L = mf + 1, jtrun, 2
                  s1_t = s(k, 1, L, m)
                  s2_t = s(k, 2, L, m)
                  p_t = poly(L, j, m)
                  dp_t = dpoly(L, j, m)
                  wcu_fk_t(1, 1) = wcu_fk_t(1, 1) + s2_t*p_t
                  wcu_fk_t(2, 1) = wcu_fk_t(2, 1) - s1_t*p_t
                  wcv_fk_t(1, 1) = wcv_fk_t(1, 1) - s1_t*dp_t
                  wcv_fk_t(2, 1) = wcv_fk_t(2, 1) - s2_t*dp_t
                  wcu_fk_t(1, 2) = wcu_fk_t(1, 2) - s2_t*p_t
                  wcu_fk_t(2, 2) = wcu_fk_t(2, 2) + s1_t*p_t
                  wcv_fk_t(1, 2) = wcv_fk_t(1, 2) - s1_t*dp_t
                  wcv_fk_t(2, 2) = wcv_fk_t(2, 2) - s2_t*dp_t
               end do
               jj = jlist2(j)
               wcc_fk(k, 1, 1, m, jj) = wcu_fk_t(1, 1)*cim(m)
               wcc_fk(k, 1, 2, m, jj) = wcu_fk_t(2, 1)*cim(m)
               wcc_fk(k, 2, 1, m, jj) = wcv_fk_t(1, 1)
               wcc_fk(k, 2, 2, m, jj) = wcv_fk_t(2, 1)

               jj = jlist2(my - j + 1)
               wcc_fk(k, 1, 1, m, jj) = wcu_fk_t(1, 2)*cim(m)
               wcc_fk(k, 1, 2, m, jj) = wcu_fk_t(2, 2)*cim(m)
               wcc_fk(k, 2, 1, m, jj) = wcv_fk_t(1, 2)
               wcc_fk(k, 2, 2, m, jj) = wcv_fk_t(2, 2)
            end if
         end do
      end do
   end do

   ! Present on device: wcu_fk, twcc_fk, wcv_fk, twdd_fk
   ! call mpe_transpose_rs1_sp_gpu(wcu_fk, twcc_fk, my_max, jtmax, lev*2, nsize, nccl_col_comm)
   ! call mpe_transpose_rs1_sp_gpu(wcv_fk, twdd_fk, my_max, jtmax, lev*2, nsize, nccl_col_comm)
   call mpe_transpose_sr_sp_gpu(wcc_fk, twcc_fk, lev*4, jtmax, my_max, &
                                nsize, nccl_col_comm, async_id)
   !$acc host_data use_device(gwk1)
   CUDACHECK(cudaMemSetAsync(gwk1, real(0.0, RTYPE), size(gwk1), stream))
   !$acc end host_data

   !$acc parallel loop collapse(2) private(j, jtrunj) async(async_id)
   do jj = 1, jlistnum
   do k = 1, lev
      j = jlist1(jj)
      jtrunj = mtrundef(j)
      !$acc loop private(mm, mp, mlst)
      do m = 1, jtrunj
         mm = 2*m - 1
         mp = mm + 1
         mlst = nlist(m)
         gwk1(mm, k, 1, jj) = twcc_fk(k, 1, 1, mlst, jj)
         gwk1(mp, k, 1, jj) = twcc_fk(k, 1, 2, mlst, jj)
         gwk1(mm, k, 2, jj) = twcc_fk(k, 2, 1, mlst, jj)
         gwk1(mp, k, 2, jj) = twcc_fk(k, 2, 2, mlst, jj)
      end do
   end do
   end do

   if (lreduce .eq. 0) then
      call rfftmlt_gpu(cc, gwk1, trigs, ifax, 1, nx + 2, nx, jlistnum*lev*2, 1)
   else
      if (.not. (fft_cg%created)) then
         call rfftmlt_loop_identical_cuda_graph(cc, gwk1, trigsj, ifaxj, jlist1, nxdef, jlistnum, nx + 2, lev*2, 1, fft_cg%graph)
         CUDACHECK(cudaGraphInstantiate(fft_cg%graph_exec, fft_cg%graph, 0))
         fft_cg%created = .true.
      end if
      CUDACHECK(cudaGraphLaunch(fft_cg%graph_exec, stream))
   end if

   call ujoinsr_gpu(cc, dlpl, dtpl, dummy, dummy, nx, my_max, levF, jlistnum, 2, 1)
   !$acc kernels async(async_id)
   dlpl = -dlpl
   dtpl = -dtpl
   !$acc end kernels
   !$acc exit data delete(twcc_fk, wcc_fk) async(async_id)

end subroutine trngra3_gpu
! ------------------------------------------------------------
subroutine trngra3_gpu_cuda_graph(jtrun, jtmax, nx, lev, my, my_max, &
                                  cim, poly, dpoly, &
                                  s, dlpl, dtpl, nsize, &
                                  cc, gwk1, wss, tcc, wcc_fk)
! Present on device: cim, poly, dpoly, s, dlpl, dtpl, mtrundef, jlist1, jlist2, nlist, mlist, nxdef, cc, gwk1, wcu_t, wcv_t, wcc_fk
!  subroutine to transform spectral terrain pressure to grid point
!  fields of zonal and meridional derivatives of terrain pressure
!
! *** input ***
!
!  jtrun: zonal wavenumber resolution limit
!  jtmax: maximum amount of zonal waves located in each pe
!  nx: e-w dimension no.
!  my: n-s dimension no.
!  my_max: maximum amount of n-s grids located in each pe
!  cim: zonal wavenumber array
!  poly: legendre polynomials
!  dpoly: d(poly)/d(sin(lat))
!  s: spectral coefficient array
!  nsize: pe number
!
! **** output ****
!
!  dlpl: d(pt)/d(longitude)
!  dtpl: d(pt)/d(sin(lat))
!
! ****************************************************
!
   use const, only: RTYPE
   use index
   use fftcom
   use spec_cuda_graph, only: fft_cg => trngra3_fft_cg, lt_cg => trngra3_lt_cg
   use openacc
   use cudafor
   use cublas

   implicit none

   integer :: jtrun, jtmax, nx, lev, my, my_max, nsize
   ! << output >>
   real(kind=RTYPE), dimension(nxp, levF, my_max) :: dlpl, dtpl
   ! << input >>
   real(kind=RTYPE), dimension(lev, 2, jtrun, jtmax) :: s
   ! << const >>
   real, dimension((jtrun + nsize)*my/2*jtmax) :: poly, dpoly
   real(kind=RTYPE), dimension(jtmax) :: cim !! type
   real :: cim_t
   ! << buffer >>
   real(kind=RTYPE), dimension(nx + 2, lev, 2, my_max) :: cc, gwk1
   real, dimension(lev, 2, jtrun, jtmax, 2) :: wss
   real, dimension(lev, 2, my, jtmax, 2) :: tcc
   real(kind=RTYPE), dimension(lev, 2, 2, jtmax, my_max*nsize) :: wcc_fk
   ! << local >>
   real(kind=RTYPE), dimension(lev, 2, 2, jtmax*nsize, my_max) :: twcc_fk
   real(kind=RTYPE) :: dummy

   integer :: myhalf, lev2, k, m, mf, l, j, jj, i, jtrunj, &
              mm, mp, mlst, jlistnum_fj, llistnum_fj, j_str, m_str

   integer :: async_id, istat, id_str
   integer(kind=cuda_stream_kind) :: stream, lt_cg_stream(jtmax*2)
   type(cudaEvent) :: spread_event, pack_event
   type(cublashandle) :: handle

   async_id = 1
   stream = acc_get_cuda_stream(async_id)
   myhalf = my/2
   lev2 = lev*2

   !$acc enter data create(twcc_fk) async(async_id)
   !$acc parallel loop collapse(3) private(mf) async(async_id)
   do m = 1, mlistnum
      do l = 1, jtrun
         do k = 1, lev
            mf = mlist(m)
            if (l .ge. mf) then
               wss(k, 1, l, m, 1) = -s(k, 2, l, m)
               wss(k, 2, l, m, 1) = s(k, 1, l, m)

               wss(k, 1, l, m, 2) = s(k, 1, l, m)
               wss(k, 2, l, m, 2) = s(k, 2, l, m)
            end if
         end do
      end do
   end do

   !$acc host_data use_device(wcc_fk)
   CUDACHECK(cudaMemSetAsync(wcc_fk, real(0.0, RTYPE), size(wcc_fk), stream))
   !$acc end host_data

   if (.not. (lt_cg%created)) then
      istat = cudaEventCreate(spread_event)
      istat = cudaEventCreate(pack_event)
      handle = cublasGetHandle()
      do m = 1, mlistnum*2
         lt_cg_stream(m) = acc_get_cuda_stream(m + 1)
      end do

      CUDACHECK(cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal))
      istat = cudaEventRecord(spread_event, stream)

      do m = 1, mlistnum
         id_str = (m - 1)*2 + 1
         istat = cudaStreamWaitEvent(lt_cg_stream(id_str), spread_event, 0)
         istat = cudaStreamWaitEvent(lt_cg_stream(id_str + 1), spread_event, 0)

         mf = mlist(m)
         jlistnum_fj = tcolt_jlist(1, m)*2
         llistnum_fj = jtrun - mf + 1
         m_str = poly_mlist(m)
         !$acc host_data use_device(wss, tcc, poly, dpoly)
         ! --------------------
         istat = cublasSetStream(handle, lt_cg_stream(id_str))
         cim_t = cim(m)
         call dgemm('n', 'n', lev2, jlistnum_fj, llistnum_fj, &
                    cim_t, wss(1, 1, mf, m, 1), lev2, &
                    poly(m_str), llistnum_fj, &
                    0.0d+0, tcc(1, 1, 1, m, 1), lev2)
         istat = cudaEventRecord(pack_event, lt_cg_stream(id_str))
         istat = cudaStreamWaitEvent(stream, pack_event, 0)
         ! --------------------
         istat = cublasSetStream(handle, lt_cg_stream(id_str + 1))
         call dgemm('n', 'n', lev2, jlistnum_fj, llistnum_fj, &
                    1.0d+0, wss(1, 1, mf, m, 2), lev2, &
                    dpoly(m_str), llistnum_fj, &
                    0.0d+0, tcc(1, 1, 1, m, 2), lev2)
         istat = cudaEventRecord(pack_event, lt_cg_stream(id_str + 1))
         istat = cudaStreamWaitEvent(stream, pack_event, 0)
         !$acc end host_data
      end do

      !$acc parallel loop collapse(3) private(jj, jlistnum_fj, j_str) async(async_id)
      do m = 1, mlistnum
         do j = 1, my
            do k = 1, lev
               jlistnum_fj = tcolt_jlist(1, m)*2
               if (j .le. jlistnum_fj) then
                  j_str = tcolt_jlist(2, m)
                  jj = jlist2(j_str + j - 1)

                  wcc_fk(k, 1, 1, m, jj) = tcc(k, 1, j, m, 1)
                  wcc_fk(k, 1, 2, m, jj) = tcc(k, 2, j, m, 1)

                  wcc_fk(k, 2, 1, m, jj) = tcc(k, 1, j, m, 2)
                  wcc_fk(k, 2, 2, m, jj) = tcc(k, 2, j, m, 2)
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
   call mpe_transpose_sr_sp_gpu(wcc_fk, twcc_fk, lev*4, jtmax, my_max, &
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
            do k = 1, lev2
               gwk1(mm, k, 1, jj) = twcc_fk(k, 1, 1, mlst, jj)
               gwk1(mp, k, 1, jj) = twcc_fk(k, 1, 2, mlst, jj)
            end do
         end if
      end do
   end do

   if (lreduce .eq. 0) then
      call rfftmlt_gpu(cc, gwk1, trigs, ifax, 1, nx + 2, nx, jlistnum*lev*2, 1)
   else
      if (.not. (fft_cg%created)) then
         call rfftmlt_loop_identical_cuda_graph(cc, gwk1, trigsj, ifaxj, jlist1, nxdef, jlistnum, nx + 2, lev*2, 1, fft_cg%graph)
         CUDACHECK(cudaGraphInstantiate(fft_cg%graph_exec, fft_cg%graph, 0))
         fft_cg%created = .true.
      end if
      CUDACHECK(cudaGraphLaunch(fft_cg%graph_exec, stream))
   end if

   call ujoinsr_gpu(cc, dlpl, dtpl, dummy, dummy, nx, my_max, levF, jlistnum, 2, 1)
   !$acc exit data delete(twcc_fk) async(async_id)

end subroutine trngra3_gpu_cuda_graph
