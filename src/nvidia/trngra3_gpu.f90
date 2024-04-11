subroutine trngra3_gpu(jtrun, jtmax, nx, lev, my, my_max, cim, poly, dpoly, s, dlpl, dtpl, nsize)
!
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
   use openacc
   use cudafor

   implicit none

   integer jtrun, jtmax, nx, lev, my, my_max, nsize
   real(kind=RTYPE) poly(jtrun, my/2, jtmax), dpoly(jtrun, my/2, jtmax) ! Present on device
   real(kind=RTYPE) s(lev, 2, jtrun, jtmax) ! Present on device
   real(kind=RTYPE) cim(jtmax) ! Present on device
   real(kind=RTYPE) dlpl(nxp, levF, my_max), dtpl(nxp, levF, my_max) ! Present on device
   integer myhalf, k, m, mf, l, j, jj, i, jtrunj, mm, mp, mlst, nxj
   real(kind=RTYPE) cc(nx + 2, lev, 2, my_max)

   real(kind=RTYPE) gwk1(nx + 2, lev, 2, my_max)

   real(kind=RTYPE) twcc_fk(my_max, jtmax*nsize, lev, 2)
   real(kind=RTYPE) twdd_fk(my_max, jtmax*nsize, lev, 2)

   real(kind=RTYPE) wcu_fk(my_max*nsize, jtmax, lev, 2), wcv_fk(my_max*nsize, jtmax, lev, 2)
   real(kind=RTYPE) wcu_t(lev, 2, my), wcv_t(lev, 2, my)

   real(kind=RTYPE) dummy
   integer async_id, istat
   integer(kind=cuda_stream_kind) :: stream

   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   !$acc enter data copyin(mtrundef, jlist1, jlist2, nlist) async(async_id)
   !$acc enter data create(wcu_t, wcv_t, wcu_fk, wcv_fk, twcc_fk, twdd_fk) async(async_id)

   myhalf = my/2

   !$acc host_data use_device(wcu_fk, wcv_fk)
   istat = cudaMemSetAsync(wcu_fk, 0.0, size(wcu_fk), stream)
   istat = cudaMemSetAsync(wcv_fk, 0.0, size(wcv_fk), stream)
   !$acc end host_data

   do m = 1, mlistnum
      mf = mlist(m)

      !$acc host_data use_device(wcu_t, wcv_t)
      istat = cudaMemSetAsync(wcu_t, 0.0, size(wcu_t), stream)
      istat = cudaMemSetAsync(wcv_t, 0.0, size(wcv_t), stream)
      !$acc end host_data

      !$acc parallel loop collapse(2) async(async_id)
      do j = 1, myhalf
         do k = 1, lev
            !$acc loop seq
            do L = mf, jtrun
               if (mf .le. mtrundef(j)) then
                  wcu_t(k, 1, j) = wcu_t(k, 1, j) + s(k, 2, L, m)*poly(L, j, m)*cim(m)
                  wcu_t(k, 2, j) = wcu_t(k, 2, j) - s(k, 1, L, m)*poly(L, j, m)*cim(m)
                  wcv_t(k, 1, j) = wcv_t(k, 1, j) - s(k, 1, L, m)*dpoly(L, j, m)
                  wcv_t(k, 2, j) = wcv_t(k, 2, j) - s(k, 2, L, m)*dpoly(L, j, m)
               end if
            end do
         end do
      end do

      !$acc parallel loop gang async(async_id)
      do jj = myhalf + 1, my
         j = my - jj + 1
         !$acc loop vector
         do k = 1, lev
            !$acc loop seq
            do L = mf, jtrun, 2
               if (mf .le. mtrundef(j)) then
                  wcu_t(k, 1, jj) = wcu_t(k, 1, jj) + s(k, 2, L, m)*poly(L, j, m)*cim(m)
                  wcu_t(k, 2, jj) = wcu_t(k, 2, jj) - s(k, 1, L, m)*poly(L, j, m)*cim(m)
                  wcv_t(k, 1, jj) = wcv_t(k, 1, jj) + s(k, 1, L, m)*dpoly(L, j, m)
                  wcv_t(k, 2, jj) = wcv_t(k, 2, jj) + s(k, 2, L, m)*dpoly(L, j, m)
               end if
            end do
            !$acc loop seq
            do L = mf + 1, jtrun, 2
               if (mf .le. mtrundef(j)) then
                  wcu_t(k, 1, jj) = wcu_t(k, 1, jj) - s(k, 2, L, m)*poly(L, j, m)*cim(m)
                  wcu_t(k, 2, jj) = wcu_t(k, 2, jj) + s(k, 1, L, m)*poly(L, j, m)*cim(m)
                  wcv_t(k, 1, jj) = wcv_t(k, 1, jj) - s(k, 1, L, m)*dpoly(L, j, m)
                  wcv_t(k, 2, jj) = wcv_t(k, 2, jj) - s(k, 2, L, m)*dpoly(L, j, m)
               end if
            end do
         end do
      end do

      !$acc parallel loop gang async(async_id)
      do j = 1, my
         jj = jlist2(j)
         !$acc loop vector
         do k = 1, lev
            wcu_fk(jj, m, k, 1) = wcu_t(k, 1, j)
            wcu_fk(jj, m, k, 2) = wcu_t(k, 2, j)
            wcv_fk(jj, m, k, 1) = wcv_t(k, 1, j)
            wcv_fk(jj, m, k, 2) = wcv_t(k, 2, j)
         end do
      end do
   end do

   call mpe_transpose_rs1_sp_gpu(wcu_fk, twcc_fk, my_max, jtmax, lev*2, nsize, col_comm)
   call mpe_transpose_rs1_sp_gpu(wcv_fk, twdd_fk, my_max, jtmax, lev*2, nsize, col_comm)
   !$acc exit data delete(wcu_t, wcv_t, jlist2, wcu_fk, wcv_fk) async(async_id)

   !$acc enter data create(cc, gwk1) async(async_id)
   !$acc host_data use_device(cc)
   istat = cudaMemSetAsync(cc, 0.0, size(cc), stream)
   !$acc end host_data

   !$acc parallel loop gang async(async_id)
   do jj = 1, jlistnum
      j = jlist1(jj)
      jtrunj = mtrundef(j)
      !$acc loop vector collapse(2)
      do k = 1, lev
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

   if (lreduce .eq. 0) then
#ifdef SP
      call rfftmlt_sp(cc, gwk1, trigs, ifax, 1, nx + 2, nx, jlistnum*lev*2, 1)
#else
      call rfftmlt_gpu(cc, gwk1, trigs, ifax, 1, nx + 2, nx, jlistnum*lev*2, 1)
#endif
   else
#ifdef SP
!$omp  parallel do default(none)                                &
!$omp  private(jj,j,nxj,gwk1)                                   &
!$omp  shared(jlistnum,jlist1,nxdef,cc,trigsj,ifaxj,nx,lev)     &
!$omp  schedule(dynamic)
      do jj = 1, jlistnum
         j = jlist1(jj)
         nxj = nxdef(j)
         call rfftmlt_sp(cc(1, 1, 1, jj), gwk1(1, 1, 1, jj), trigsj(1, j), ifaxj(1, j), 1, nx + 2, nxj, lev*2, 1)
      end do
!$omp end parallel do
#else
      call rfftmlt_loop_identical(cc, gwk1, trigsj, ifaxj, jlist1, nxdef, jlistnum, nx + 2, lev*2, 1)

#endif
   end if

   call ujoinsr_gpu(cc, dlpl, dtpl, dummy, dummy, nx, my_max, levF, jlistnum, 2, 1)
   !$acc kernels async(async_id)
   dlpl = -dlpl
   dtpl = -dtpl
   !$acc end kernels
   !$acc exit data delete(jlist1, mtrundef, nlist, twcc_fk, twdd_fk, gwk1, cc) async(async_id)
   !$acc wait(async_id)

   return
end
