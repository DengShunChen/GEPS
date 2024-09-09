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
   use fft_cuda_graph
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
   real(kind=RTYPE), dimension(my_max, jtmax*nsize, lev, 2) :: twcc_fk, twdd_fk
   real(kind=RTYPE), dimension(my_max*nsize, jtmax, lev, 2) :: wcu_fk, wcv_fk
   real(kind=RTYPE), dimension(2, 2) :: wcu_fk_t, wcv_fk_t
   real(kind=RTYPE) :: dummy
   integer :: async_id, istat
   integer(kind=cuda_stream_kind) :: stream

   async_id = 1
   stream = acc_get_cuda_stream(async_id)
   myhalf = my/2

   !$acc enter data create(wcu_fk, wcv_fk, twcc_fk, twdd_fk) async(async_id)

   !$acc host_data use_device(wcu_fk, wcv_fk)
   CUDACHECK(cudaMemSetAsync(wcu_fk, 0.0, size(wcu_fk), stream))
   CUDACHECK(cudaMemSetAsync(wcv_fk, 0.0, size(wcv_fk), stream))
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
               wcu_fk(jj, m, k, 1) = wcu_fk_t(1, 1)*cim(m)
               wcu_fk(jj, m, k, 2) = wcu_fk_t(2, 1)*cim(m)
               wcv_fk(jj, m, k, 1) = wcv_fk_t(1, 1)
               wcv_fk(jj, m, k, 2) = wcv_fk_t(2, 1)
               jj = jlist2(my - j + 1)
               wcu_fk(jj, m, k, 1) = wcu_fk_t(1, 2)*cim(m)
               wcu_fk(jj, m, k, 2) = wcu_fk_t(2, 2)*cim(m)
               wcv_fk(jj, m, k, 1) = wcv_fk_t(1, 2)
               wcv_fk(jj, m, k, 2) = wcv_fk_t(2, 2)
            end if
         end do
      end do
   end do

   ! Present on device: wcu_fk, twcc_fk, wcv_fk, twdd_fk
   call mpe_transpose_rs1_sp_gpu(wcu_fk, twcc_fk, my_max, jtmax, lev*2, nsize, nccl_col_comm)
   call mpe_transpose_rs1_sp_gpu(wcv_fk, twdd_fk, my_max, jtmax, lev*2, nsize, nccl_col_comm)

   !$acc host_data use_device(cc)
   CUDACHECK(cudaMemSetAsync(cc, 0.0, size(cc), stream))
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
         cc(mm, k, 1, jj) = twcc_fk(jj, mlst, k, 1)
         cc(mp, k, 1, jj) = twcc_fk(jj, mlst, k, 2)
         cc(mm, k, 2, jj) = twdd_fk(jj, mlst, k, 1)
         cc(mp, k, 2, jj) = twdd_fk(jj, mlst, k, 2)
      end do
   end do
   end do

#ifdef SP
   print *, "Symbol SP is not supported."
   call exit(1)
#endif

   if (lreduce .eq. 0) then
      call rfftmlt_gpu(cc, gwk1, trigs, ifax, 1, nx + 2, nx, jlistnum*lev*2, 1)
   else
      if (trngra3_graph_created) then
         CUDACHECK(cudaGraphLaunch(trngra3_graph_exec, stream))
      else
         call rfftmlt_loop_identical_cuda_graph(cc, gwk1, trigsj, ifaxj, jlist1, nxdef, jlistnum, nx + 2, lev*2, 1, trngra3_graph)
         CUDACHECK(cudaGraphInstantiate(trngra3_graph_exec, trngra3_graph, trngra3_error_node, trngra3_buffer, trngra3_buffer_len))
         trngra3_graph_created = .true.
         CUDACHECK(cudaGraphLaunch(trngra3_graph_exec, stream))
      end if
   end if

   call ujoinsr_gpu(cc, dlpl, dtpl, dummy, dummy, nx, my_max, levF, jlistnum, 2, 1)
   !$acc kernels async(async_id)
   dlpl = -dlpl
   dtpl = -dtpl
   !$acc end kernels
   !$acc exit data delete(twcc_fk, twdd_fk, wcu_fk, wcv_fk) async(async_id)

end subroutine trngra3_gpu
