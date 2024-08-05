#define CUDACHECK(ierr) call cuda_check_helper(ierr, __FILE__, __LINE__)

subroutine tranrs_gpu(jtrun, jtmax, nx, my, my_max, lev, poly, w, cc &
                      , wss, num, nsize, gwk1)
!  Present on device: poly, w, cc, wss, nlist, jlist2, jlist1, nxdef
!  subroutine to transform a scalar grid point field to spectral
!  coefficients
!
! *** input ***
!
!  jtrun: zonal wavenumber truncation limit
!  jtmax: maximum amount of zonal waves located in each pe
!  nx: e-w dimension no.
!  my: n-s dimension no.
!  my_max: maximum amount of n-s grids located in each pe
!  lev: number of vertical levels to transform
!  poly: legendre polynomials
!  w: gaussian quadrature weights
!  cc_r8: 3-dim input grid pt. field to be transformed
!  num: number of variables grouped together
!
! *** output ***
!
!  wss: spectral coefficient fields
!
!  **********************************
!
   use const, only: RTYPE
   use index
   use fftcom
   use fft_cuda_graph
   use fj_pad
   use openacc
   use cudafor

   implicit none

   integer jtrun, jtmax, nx, my, my_max, lev, num, nsize
   integer mlx, myhalf, lev2, nxj, j, jj, mchk, jtrunj, m, mm, mp
   integer mlst, ii, k, mm1, mp1, mlst1, mm2, mp2, mlst2, mm3, mp3, mlst3, mf
   integer lchk, lle, jlistnum_fj, j_fj, j1, j2, l, l_fj, llistnum_fj

   real(kind=RTYPE) poly(jtrun, my/2, jtmax), w(my)
   real(kind=RTYPE) wss(lev, 2, num, jtrun, jtmax)
   real(kind=RTYPE) gwk1(nx + 2, lev, num, my_max)
   real(kind=RTYPE) wcc_fk(lev, 2, num, jtmax, my_max*nsize)
   real(kind=RTYPE) twcc_fk(lev, 2, num, jtmax*nsize, my_max)
   real(kind=RTYPE) cc(nx + 2, lev, num, my_max)
   real fj_wss_sum(lev*2*num, jtrun)
   real fj_wss_dif(lev*2*num, jtrun)
   real fj_polyw_sum(my/2 + npad, jtrun)
   real fj_polyw_dif(my/2 + npad, jtrun)
   real fj_wccSUM(lev*2*num, my/2)
   real fj_wccDIF(lev*2*num, my/2)
   integer jlist_fj(my/2, mlistnum)
   integer jlistnum_fj_array(mlistnum)
   integer async_id, istat
   integer(kind=cuda_stream_kind) stream

   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   !$acc enter data create(twcc_fk, wcc_fk) async(async_id)

   !$acc host_data use_device(twcc_fk, wss, gwk1)
   istat = cudaMemsetAsync(twcc_fk, 0.0, size(twcc_fk), stream)
   istat = cudaMemsetAsync(wss, 0.0, size(wss), stream)
   istat = cudaMemsetAsync(gwk1, 0.0, size(gwk1), stream)
   !$acc end host_data

   mlx = (jtrun/2)*((jtrun + 1)/2)
   myhalf = my/2
   lev2 = lev*2
!
!  fft for each guassian latitude of 2-d field
!
#ifdef SP
   print *, "Symbol SP is not supported."
   call exit(1)
#endif

   if (length_fft .eq. 0 .and. lreduce .eq. 0) then
      call rfftmlt(cc, gwk1, trigs, ifax, 1, nx + 2, nx, lev*jlistnum*num, -1)
   else
      if (tranrs_graph_created) then
         CUDACHECK(cudaGraphLaunch(tranrs_graph_exec, stream))
      else
         call rfftmlt_loop_identical_cuda_graph(cc, gwk1, trigsj, ifaxj, jlist1, nxdef, jlistnum, nx + 2, lev*num, -1, tranrs_graph)
         CUDACHECK(cudaGraphInstantiate(tranrs_graph_exec, tranrs_graph, tranrs_error_node, tranrs_buffer, tranrs_buffer_len))
         tranrs_graph_created = .true.
         CUDACHECK(cudaGraphLaunch(tranrs_graph_exec, stream))
      end if
   end if

   mchk = iand(jtrun, 3)

   do j = 1, jlistnum
      jj = jlist1(j)
      jtrunj = mtrundef(jj)
      mchk = iand(jtrunj, 3)
      !$acc parallel loop gang async(async_id)
      do m = 1, mchk
         mm = 2*m - 1
         mp = mm + 1
         mlst = nlist(m)
         !$acc loop vector collapse(2)
         do ii = 1, num
            do k = 1, lev
               twcc_fk(k, 1, ii, mlst, j) = cc(mm, k, ii, j)
               twcc_fk(k, 2, ii, mlst, j) = cc(mp, k, ii, j)
            end do
         end do
      end do

      !$acc parallel loop gang async(async_id)
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
               twcc_fk(k, 1, ii, mlst, j) = cc(mm, k, ii, j)
               twcc_fk(k, 2, ii, mlst, j) = cc(mp, k, ii, j)
               twcc_fk(k, 1, ii, mlst1, j) = cc(mm1, k, ii, j)
               twcc_fk(k, 2, ii, mlst1, j) = cc(mp1, k, ii, j)
               twcc_fk(k, 1, ii, mlst2, j) = cc(mm2, k, ii, j)
               twcc_fk(k, 2, ii, mlst2, j) = cc(mp2, k, ii, j)
               twcc_fk(k, 1, ii, mlst3, j) = cc(mm3, k, ii, j)
               twcc_fk(k, 2, ii, mlst3, j) = cc(mp3, k, ii, j)
            end do
         end do
      end do

   end do

   ! Present on device: twcc_fk, wcc_fk
   call mpe_transpose_rs_sp_gpu(twcc_fk, wcc_fk, lev*2*num, jtmax, my_max, nsize, nccl_col_comm)
   !$acc exit data delete(twcc_fk) async(async_id)

   do m = 1, mlistnum
      jlistnum_fj = 0
      do j = 1, myhalf
         if (mf .le. mtrundef(j)) then
            jlistnum_fj = jlistnum_fj + 1
            jlist_fj(jlistnum_fj, m) = j
         end if
      end do
      jlistnum_fj_array(m) = jlistnum_fj
   end do
   !$acc enter data create(fj_wss_sum, fj_wss_dif, fj_polyw_sum, fj_polyw_dif, fj_wccSUM, fj_wccDIF) copyin(jlist_fj) async(async_id)

   do m = 1, mlistnum
      mf = mlist(m)

      !$acc kernels async(async_id)
      fj_wss_sum = 0.0
      fj_wss_dif = 0.0
      fj_polyw_sum = 0.0
      fj_polyw_dif = 0.0
      fj_wccSUM = 0.0
      fj_wccDIF = 0.0
      !$acc end kernels

      lchk = iand(jtrun - mf + 1, 1)
      lle = jtrun - lchk
      jlistnum_fj = jlistnum_fj_array(m)

      !$acc parallel loop gang async(async_id)
      do j_fj = 1, jlistnum_fj
         j = jlist_fj(j_fj, m)
         j1 = jlist2(j)
         j2 = jlist2(my - j + 1)
         !$acc loop vector
         do k = 1, lev2*num
            fj_wccSUM(k, j_fj) = (wcc_fk(k, 1, 1, m, j1) + wcc_fk(k, 1, 1, m, j2))
            fj_wccDIF(k, j_fj) = (wcc_fk(k, 1, 1, m, j1) - wcc_fk(k, 1, 1, m, j2))
         end do
      end do

      if (jtrun .gt. mf) then
         !$acc parallel loop gang async(async_id)
         do l = mf, lle, 2
            l_fj = (l - mf + 2)/2
            !$acc loop vector
            do j_fj = 1, jlistnum_fj
               j = jlist_fj(j_fj, m)
               fj_polyw_sum(j_fj, l_fj) = poly(l, j, m)*w(j)
               fj_polyw_dif(j_fj, l_fj) = poly(l + 1, j, m)*w(j)
            end do
         end do

         do l = mf, lle, 2
            l_fj = (l - mf + 2)/2
         end do
         llistnum_fj = l_fj

         call dgemm_async('n', 'n', lev2*num, llistnum_fj, jlistnum_fj, 1.0d+0, fj_wccSUM, lev*2*num, &
                          fj_polyw_sum, my/2 + npad, 1.0d+0, fj_wss_sum, lev*2*num, async_id)

         !$acc parallel loop gang async(async_id)
         do l = mf, lle, 2
            l_fj = (l - mf + 2)/2
            !$acc loop vector
            do k = 1, lev2*num
               !$acc atomic
               wss(k, 1, 1, l, m) = wss(k, 1, 1, l, m) + fj_wss_sum(k, l_fj)
            end do
         end do

         call dgemm_async('n', 'n', lev2*num, llistnum_fj, jlistnum_fj, 1.0d+0, fj_wccDIF, lev*2*num, &
                          fj_polyw_dif, my/2 + npad, 1.0d+0, fj_wss_dif, lev*2*num, async_id)
         !$acc parallel loop gang async(async_id)
         do l = mf, lle, 2
            l_fj = (l - mf + 2)/2
            !$acc loop vector
            do k = 1, lev2*num
               !$acc atomic
               wss(k, 1, 1, l + 1, m) = wss(k, 1, 1, l + 1, m) + fj_wss_dif(k, l_fj)
            end do
         end do

      end if

      if (lchk .eq. 1) then
         l = jtrun
         !$acc parallel loop async(async_id)
         do j_fj = 1, jlistnum_fj
            j = jlist_fj(j_fj, m)
            fj_polyw_sum(j_fj, l) = poly(l, j, m)*w(j)
         end do
         !$acc parallel loop gang async(async_id)
         do j_fj = 1, jlistnum_fj
            !$acc loop vector
            do k = 1, lev2*num
               !$acc atomic
               wss(k, 1, 1, l, m) = wss(k, 1, 1, l, m) + fj_polyw_sum(j_fj, l)*fj_wccSUM(k, j_fj)
            end do
         end do

      end if
   end do
   !$acc exit data delete(jlist_fj, wcc_fk, fj_wss_sum, fj_wss_dif, fj_polyw_sum, fj_polyw_dif, fj_wccSUM, fj_wccDIF) async(async_id)

   return
end
