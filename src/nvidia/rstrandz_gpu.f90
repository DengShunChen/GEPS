subroutine rstrandz_gpu(jtrun, jtmax, nx, my, my_max, lev &
                        , vdmer, vdzon, w, cim, onocos, poly, dpoly &
                        , hldten, vorten, nsize)
   ! Present on device: vdmer, vdzon, w, cim, onocos, poly, dpoly, hldten, vorten
   ! Present on device: nlist, jlist2, mtrundef, jlist1, nxjlen, nxjlen_all
   use const, only: RTYPE
   use index
   use paramt
   use fftcom
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

   !$acc enter data create(twcc_fk, gwk1, cc, wcc_fk) async(async_id)

   !$acc host_data use_device(twcc_fk, gwk1)
   istat = cudaMemSetAsync(twcc_fk, 0.0, size(twcc_fk), stream)
   istat = cudaMemSetAsync(gwk1, 0.0, size(gwk1), stream)
   !$acc end host_data

   myhalf = my/2
   lev2 = lev*2

   ! Present on device: cc, vdmer, vdzon, jlist1, nxjlen, nxjlen_all
   call joinrs_gpu(cc, vdmer, vdzon, dummy, dummy, nx, my_max, levf, jlistnum, 2, 1)

#ifdef SP
   print *, "Symbol SP is not supported."
   call exit(1)
#endif

   if (length_fft .eq. 0 .and. lreduce .eq. 0) then
      call rfftmlt(cc, gwk1, trigs, ifax, 1, nx + 2, nx, lev*jlistnum*2, -1)
   else
      call rfftmlt_loop_identical(cc, gwk1, trigsj, ifaxj, jlist1, nxdef, jlistnum, nx + 2, lev*2, -1)
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
   call mpe_transpose_rs_sp_gpu(twcc_fk, wcc_fk, lev*2*2, jtmax, my_max, nsize, col_comm)
#else
   ! Present on device: twcc_fk, wcc_fk
   call mpe_transpose_rs_gpu(twcc_fk, wcc_fk, lev*2*2, jtmax, my_max, nsize, nccl_col_comm)
#endif
   !$acc exit data delete(gwk1, cc, twcc_fk) async(async_id)

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

               wcc2(k, 1, 1, j) = +(wcc_fk(k, 1, 1, m, j1) - wcc_fk(k, 1, 1, m, j2))
               wcc2(k, 2, 1, j) = +(wcc_fk(k, 2, 1, m, j1) - wcc_fk(k, 2, 1, m, j2))
               wcc2(k, 1, 2, j) = -(wcc_fk(k, 1, 2, m, j1) - wcc_fk(k, 1, 2, m, j2))
               wcc2(k, 2, 2, j) = -(wcc_fk(k, 2, 2, m, j1) - wcc_fk(k, 2, 2, m, j2))

               wcc3(k, 1, 1, j) = -(wcc_fk(k, 2, 2, m, j1) + wcc_fk(k, 2, 2, m, j2))
               wcc3(k, 2, 1, j) = +(wcc_fk(k, 1, 2, m, j1) + wcc_fk(k, 1, 2, m, j2))
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
      istat = cudaMemSetAsync(wss, 0.0, size(wss), stream)
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
end
