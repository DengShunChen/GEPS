subroutine trngra3_gpu(jtrun, jtmax, nx, lev, my, my_max, cim, poly, dpoly, s, dlpl, dtpl, nsize)
! Present on device: cim, poly, dpoly, s, dlpl, dtpl, mtrundef, jlist1, jlist2, nlist, mlist
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
   use mpi

   implicit none

   integer jtrun, jtmax, nx, lev, my, my_max, nsize
   real(kind=RTYPE) poly(jtrun, my/2, jtmax), dpoly(jtrun, my/2, jtmax)
   real(kind=RTYPE) s(lev, 2, jtrun, jtmax)
   real(kind=RTYPE) cim(jtmax)
   real(kind=RTYPE) dlpl(nxp, levF, my_max), dtpl(nxp, levF, my_max)
   integer myhalf, k, m, mf, l, j, jj, i, jtrunj, mm, mp, mlst, nxj, ii
   real(kind=RTYPE) cc(nx + 2, lev, 2, my_max)
   real(kind=RTYPE) gwk1(nx + 2, lev, 2, my_max)
   real(kind=RTYPE) twcc_fk(my_max, jtmax*nsize, lev, 2)
   real(kind=RTYPE) twdd_fk(my_max, jtmax*nsize, lev, 2)
   real(kind=RTYPE) wcu_fk(my_max*nsize, jtmax, lev, 2), wcv_fk(my_max*nsize, jtmax, lev, 2)
   real(kind=RTYPE) wcu_t(lev, 2, my), wcv_t(lev, 2, my), wcu_t1, wcu_t2, wcv_t1, wcv_t2
   real(kind=RTYPE) dummy
   integer async_id, istat
   integer(kind=cuda_stream_kind) :: stream

   async_id = 1
   stream = acc_get_cuda_stream(async_id)
   myhalf = my/2

   !$acc enter data create(wcu_t, wcv_t, wcu_fk, wcv_fk, twcc_fk, twdd_fk, cc, gwk1) async(async_id)

   !$acc host_data use_device(wcu_fk, wcv_fk)
   istat = cudaMemSetAsync(wcu_fk, 0.0, size(wcu_fk), stream)
   istat = cudaMemSetAsync(wcv_fk, 0.0, size(wcv_fk), stream)
   !$acc end host_data

   !$acc parallel loop collapse(3) private(mf, jj, wcu_t1, wcu_t2, wcv_t1, wcv_t2, i, ii) async(async_id)
   do m = 1, mlistnum
      do j = 1, my
         do k = 1, lev
            mf = mlist(m)
            if (j .le. myhalf) then
               jj = j
            else
               jj = my - j + 1
            end if
            if (mf .le. mtrundef(jj)) then
               wcu_t1 = 0.0
               wcu_t2 = 0.0
               wcv_t1 = 0.0
               wcv_t2 = 0.0
               if (j .le. myhalf) then
                  !$acc loop seq
                  do L = mf, jtrun
                     wcu_t1 = wcu_t1 + s(k, 2, L, m)*poly(L, jj, m)
                     wcu_t2 = wcu_t2 - s(k, 1, L, m)*poly(L, jj, m)
                     wcv_t1 = wcv_t1 - s(k, 1, L, m)*dpoly(L, jj, m)
                     wcv_t2 = wcv_t2 - s(k, 2, L, m)*dpoly(L, jj, m)
                  end do
               else
                  !$acc loop seq
                  do L = mf, jtrun, 2
                     wcu_t1 = wcu_t1 + s(k, 2, L, m)*poly(L, jj, m)
                     wcu_t2 = wcu_t2 - s(k, 1, L, m)*poly(L, jj, m)
                     wcv_t1 = wcv_t1 + s(k, 1, L, m)*dpoly(L, jj, m)
                     wcv_t2 = wcv_t2 + s(k, 2, L, m)*dpoly(L, jj, m)
                  end do
                  !$acc loop seq
                  do L = mf + 1, jtrun, 2
                     wcu_t1 = wcu_t1 - s(k, 2, L, m)*poly(L, jj, m)
                     wcu_t2 = wcu_t2 + s(k, 1, L, m)*poly(L, jj, m)
                     wcv_t1 = wcv_t1 - s(k, 1, L, m)*dpoly(L, jj, m)
                     wcv_t2 = wcv_t2 - s(k, 2, L, m)*dpoly(L, jj, m)
                  end do
               end if
               jj = jlist2(j)
               wcu_fk(jj, m, k, 1) = wcu_t1*cim(m)
               wcu_fk(jj, m, k, 2) = wcu_t2*cim(m)
               wcv_fk(jj, m, k, 1) = wcv_t1
               wcv_fk(jj, m, k, 2) = wcv_t2
            end if
         end do
      end do
   end do

   ! Present on device: wcu_fk, twcc_fk, wcv_fk, twdd_fk
   call mpe_transpose_rs1_sp_gpu(wcu_fk, twcc_fk, my_max, jtmax, lev*2, nsize, col_comm)
   call mpe_transpose_rs1_sp_gpu(wcv_fk, twdd_fk, my_max, jtmax, lev*2, nsize, col_comm)

   !$acc host_data use_device(cc)
   istat = cudaMemSetAsync(cc, 0.0, size(cc), stream)
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
      call rfftmlt_loop_identical(cc, gwk1, trigsj, ifaxj, jlist1, nxdef, jlistnum, nx + 2, lev*2, 1)
   end if

   call ujoinsr_gpu(cc, dlpl, dtpl, dummy, dummy, nx, my_max, levF, jlistnum, 2, 1)
   !$acc kernels async(async_id)
   dlpl = -dlpl
   dtpl = -dtpl
   !$acc end kernels
   !$acc exit data delete(twcc_fk, twdd_fk, gwk1, cc, wcu_t, wcv_t, wcu_fk, wcv_fk) async(async_id)

   return
end
