#define CUDACHECK(ierr) call cuda_check_helper(ierr, __FILE__, __LINE__)

subroutine transr1_gpu(jtrun, jtmax, nx, my, my_max, poly, s, r, nsize, cc, gwk1)
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
!  poly: legendre polynomials
!  s: spectral coefficient array to transform
!
! *** output ***
!
!  r: 3-d output grid point fields
!
!  **************************************
!
   use const, only: RTYPE
   use index
   use fftcom
   use fft_cuda_graph
   use openacc
   use cudafor

   implicit none

   integer jtrun, jtmax, nx, my, my_max, nsize

   integer mlx, myhalf, m, mf, l, jlistnum_fj, j, l_fj, j_fj
   integer llistnum_fj, jj, i, jtrunj, mm, mp, mlst, nxj

   real(kind=RTYPE) r(nxp, my_max) ! Present on device
!
   real(kind=RTYPE) s(jtrun, jtmax, 2) ! Present on device
   real(kind=RTYPE) gwk1(nx + 2, my_max)
!
   real(kind=RTYPE) wcc_fk(my_max*nsize, jtmax, 2), twcc_fk(my_max, jtmax*nsize, 2)
   real wss(jtrun, 2)

   real(kind=RTYPE) poly(jtrun, my/2, jtmax) ! Present on device
   real(kind=RTYPE) cc(nx + 2, my_max)
   real tcc(my, 2)
!
   real ws2(jtrun, 2)
!
   integer jlist_fj(my/2, mlistnum)
   integer jlistnum_fj_array(mlistnum)
   real fj_poly(jtrun, my/2)
   real fj_tcc(2, my)
   real fj_wss(2, jtrun)
   real fj_ws2(2, jtrun)
   integer async_id
   integer(kind=cuda_stream_kind) :: stream

   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   mlx = (jtrun/2)*((jtrun + 1)/2)
   myhalf = my/2

   !$acc enter data create(wcc_fk) async(async_id)
   !$acc host_data use_device(wcc_fk)
   CUDACHECK(cudaMemSetAsync(wcc_fk, 0.0, size(wcc_fk), stream))
   !$acc end host_data

   do m = 1, mlistnum
      mf = mlist(m)
      jlistnum_fj = 0
      do j = 1, myhalf
         if (mf .le. mtrundef(j)) then
            jlistnum_fj = jlistnum_fj + 1
            jlist_fj(jlistnum_fj, m) = j
         end if
      end do
      jlistnum_fj_array(m) = jlistnum_fj
   end do

   !$acc enter data create(wss, ws2, tcc, fj_tcc, fj_poly, fj_wss, fj_ws2) copyin(jlist_fj) async(async_id)
   do m = 1, mlistnum
      mf = mlist(m)

      !$acc parallel loop async(async_id)
      do l = 1, jtrun
         wss(l, 1) = s(l, m, 1)
         wss(l, 2) = s(l, m, 2)
      end do

      !$acc parallel loop async(async_id)
      do l = mf, jtrun, 2
         ws2(l, 1) = wss(l, 1)
         ws2(l, 2) = wss(l, 2)
      end do
      !$acc parallel loop async(async_id)
      do l = mf + 1, jtrun, 2
         ws2(l, 1) = -wss(l, 1)
         ws2(l, 2) = -wss(l, 2)
      end do
      !$acc host_data use_device(tcc, fj_tcc, fj_poly, fj_wss, fj_ws2)
      CUDACHECK(cudaMemSetAsync(tcc, 0.0, size(tcc), stream))
      CUDACHECK(cudaMemSetAsync(fj_tcc, 0.0, size(fj_tcc), stream))
      CUDACHECK(cudaMemSetAsync(fj_poly, 0.0, size(fj_poly), stream))
      CUDACHECK(cudaMemSetAsync(fj_wss, 0.0, size(fj_wss), stream))
      CUDACHECK(cudaMemSetAsync(fj_ws2, 0.0, size(fj_ws2), stream))
      !$acc end host_data

      jlistnum_fj = jlistnum_fj_array(m)

      !$acc parallel loop gang async(async_id)
      do l = mf, jtrun
         l_fj = l - mf + 1
         !$acc loop vector
         do j_fj = 1, jlistnum_fj
            j = jlist_fj(j_fj, m)
            fj_poly(l_fj, j_fj) = poly(l, j, m)
         end do
      end do

      !$acc parallel loop async(async_id)
      do l = mf, jtrun
         l_fj = l - mf + 1
         fj_wss(1, l_fj) = wss(l, 1)
         fj_wss(2, l_fj) = wss(l, 2)
      end do

      !$acc parallel loop async(async_id)
      do l = mf, jtrun
         l_fj = l - mf + 1
         fj_ws2(1, l_fj) = ws2(l, 1)
         fj_ws2(2, l_fj) = ws2(l, 2)
      end do

      llistnum_fj = jtrun - mf + 1
      call dgemm_async('n', 'n', 2, jlistnum_fj, llistnum_fj, 1.0d+0, fj_wss, 2 &
                       , fj_poly, jtrun, 1.0d+0, fj_tcc, 2, async_id)

      !$acc parallel loop async(async_id)
      do j_fj = 1, jlistnum_fj
         j = jlist_fj(j_fj, m)
         ! Elements are unique in jlist_fj
         tcc(j, 1) = tcc(j, 1) + fj_tcc(1, j_fj)
         tcc(j, 2) = tcc(j, 2) + fj_tcc(2, j_fj)
      end do

      !$acc host_data use_device(fj_tcc)
      CUDACHECK(cudaMemSetAsync(fj_tcc, 0.0, size(fj_tcc), stream))
      !$acc end host_data

      llistnum_fj = jtrun - mf + 1
      call dgemm_async('n', 'n', 2, jlistnum_fj, llistnum_fj, 1.0d+0, fj_ws2, 2 &
                       , fj_poly, jtrun, 1.0d+0, fj_tcc, 2, async_id)

      !$acc parallel loop async(async_id)
      do j_fj = 1, jlistnum_fj
         j = jlist_fj(j_fj, m)
         jj = my - j + 1
         tcc(jj, 1) = tcc(jj, 1) + fj_tcc(1, j_fj)
         tcc(jj, 2) = tcc(jj, 2) + fj_tcc(2, j_fj)
      end do
      !$acc parallel loop async(async_id)
      do j = 1, my
         wcc_fk(jlist2(j), m, 1) = tcc(j, 1)
         wcc_fk(jlist2(j), m, 2) = tcc(j, 2)
      end do

   end do

!*** r1 start ***

   !$acc enter data create(twcc_fk) async(async_id)
   call mpe_transpose_rs1_sp_gpu(wcc_fk, twcc_fk, my_max, jtmax, 2, nsize, nccl_col_comm, async_id)

   !$acc host_data use_device(cc)
   CUDACHECK(cudaMemSetAsync(cc, 0.0, size(cc), stream))
   !$acc end host_data
   !$acc parallel loop async(async_id)
   do jj = 1, jlistnum
      j = jlist1(jj)
      jtrunj = mtrundef(j)
      !$acc loop vector
      do m = 1, jtrunj
         mm = 2*m - 1
         mp = mm + 1
         mlst = nlist(m)
         cc(mm, jj) = twcc_fk(jj, mlst, 1)
         cc(mp, jj) = twcc_fk(jj, mlst, 2)
      end do
   end do

#ifdef SP
   print *, "Symbol SP is not supported."
#endif

   if (lreduce .eq. 0) then
      call rfftmlt(cc, gwk1, trigs, ifax, 1, nx + 2, nx, jlistnum, 1)
   else
      if (transr1_graph_created) then
         CUDACHECK(cudaGraphLaunch(transr1_graph_exec, stream))
      else
         call rfftmlt_loop_identical_cuda_graph(cc, gwk1, trigsj, ifaxj, jlist1, nxdef, jlistnum, nx + 2, 1, 1, transr1_graph)
         CUDACHECK(cudaGraphInstantiate(transr1_graph_exec, transr1_graph, transr1_error_node, transr1_buffer, transr1_buffer_len))
         transr1_graph_created = .true.
         CUDACHECK(cudaGraphLaunch(transr1_graph_exec, stream))
      end if
   end if

   !$acc parallel loop async(async_id)
   do jj = 1, jlistnum
      j = jlist1(jj)
      r(1:nxjlen(j), jj) = cc(nxjstart(j):nxjend(j), jj)
   end do
   !$acc exit data delete(twcc_fk, jlist_fj, wss, ws2, tcc, fj_tcc, fj_poly, fj_wss, fj_ws2, wcc_fk) async(async_id)

!*** r1  end  ***

   return
end
