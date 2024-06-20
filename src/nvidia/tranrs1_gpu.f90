subroutine tranrs1_gpu(jtrun, jtmax, nx, my, my_max, poly, w, r, s, nsize)
   ! Present on device: poly, w, r, s, jlist1, nxdef, mtrundef, nlist, mlist, jlist2
!
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
!  poly: legendre polynomials
!  w: gaussian quadrature weights
!  r: 3-dim input grid pt. field to be transformed
!
! *** output ***
!
!  s: spectral coefficient fields
!
!  **********************************
!
   use const, only: RTYPE
   use index
   use paramt
   use fftcom
   use openacc
   use cudafor

   implicit none

   integer jtrun, jtmax, nx, my, my_max, nsize
   integer mlx, myhalf, jj, j, nxj, i, jtrunj, m, mm, mp, mlst, mf
   integer l, i1, i2, i3, j1, j2

   real(kind=RTYPE) poly(jtrun, my/2, jtmax), w(my)
   real(kind=RTYPE) r(nx, my_max)
   real(kind=RTYPE) s(jtrun, jtmax, 2)
   real(kind=RTYPE) gwk1(nx + 2, my_max)

   real(kind=RTYPE) wcc_fk(jtmax, my_max*nsize, 2)
   real(kind=RTYPE) twcc_fk(jtmax*nsize, my_max, 2)
   real(kind=RTYPE) wss(jtrun, 2)
   real(kind=RTYPE) cc(nx + 2, my_max)

   real(kind=RTYPE) wccSUM(my, 2), wccSUM1, wccSUM2
   real(kind=RTYPE) wccDIF(my, 2), wccDIF1, wccDIF2
   real(kind=RTYPE) fj_polyw(my/2, jtrun)
   logical wfirst
   data wfirst/.true./
   save wfirst
   integer async_id, istat
   integer(kind=cuda_stream_kind) :: stream

   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   !$acc enter data create(twcc_fk, gwk1, cc, wcc_fk) async(async_id)
   !$acc host_data use_device(twcc_fk, gwk1)
   istat = cudaMemSetAsync(twcc_fk, 0.0, size(twcc_fk), stream)
   istat = cudaMemSetAsync(gwk1, 0.0, size(gwk1), stream)
   !$acc end host_data

   myhalf = my/2

   !$acc parallel loop gang async(async_id) private(j, nxj)
   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef(j)
      !$acc loop vector
      do i = 1, nxj
         cc(i, jj) = r(i, jj)
      end do
   end do

#ifdef SP
   print *, "Symbol SP is not supported."
   call exit(1)
#endif

   if (lreduce .eq. 0) then
      call rfftmlt(cc, gwk1, trigs, ifax, 1, nx + 2, nx, jlistnum, -1)
   else
      call rfftmlt_loop_identical(cc, gwk1, trigsj, ifaxj, jlist1, nxdef, jlistnum, nx + 2, 1, -1)
   end if

   !$acc parallel loop gang async(async_id) private(jj, jtrunj)
   do j = 1, jlistnum
      jj = jlist1(j)
      jtrunj = mtrundef(jj)
      !$acc loop vector private(mm, mp, mlst)
      do m = 1, jtrunj
         mm = 2*m - 1
         mp = mm + 1
         mlst = nlist(m)
         twcc_fk(mlst, j, 1) = cc(mm, j)
         twcc_fk(mlst, j, 2) = cc(mp, j)
      end do
   end do

   ! Present on device: twcc_fk, wcc_fk
   call mpe_transpose_rs1_sp_gpu(twcc_fk, wcc_fk, jtmax, my_max, 2, nsize, nccl_col_comm, async_id)
   !$acc exit data delete(twcc_fk, gwk1, cc) async(async_id)
   !$acc enter data create(wss, fj_polyw, wccSUM, wccDIF) async(async_id)
   !$acc parallel loop gang async(async_id) private(mf, wss, fj_polyw, i1, i2, i3, wccSUM, wccDIF)
   do m = 1, mlistnum
      mf = mlist(m)

      !$acc loop vector
      do L = mf, jtrun
         wss(L, 1) = 0.
         wss(L, 2) = 0.
      end do

      !$acc loop vector collapse(2)
      do j = 1, myhalf
         do L = mf, jtrun
            fj_polyw(j, L) = poly(L, j, m)*w(j)
         end do
      end do

      i1 = (jtrun - mf + 1)/4
      i2 = (jtrun - mf + 1 - i1*4)/2
      i3 = jtrun - mf + 1 - i1*4 - i2*2

      !$acc loop vector private(j1, j2)
      do j = 1, myhalf
         j1 = jlist2(j)
         j2 = jlist2(my - j + 1)
         if (mf .le. mtrundef(j)) then
            wccSUM(j, 1) = wcc_fk(m, j1, 1) + wcc_fk(m, j2, 1)
            wccDIF(j, 1) = wcc_fk(m, j1, 1) - wcc_fk(m, j2, 1)
            wccSUM(j, 2) = wcc_fk(m, j1, 2) + wcc_fk(m, j2, 2)
            wccDIF(j, 2) = wcc_fk(m, j1, 2) - wcc_fk(m, j2, 2)
         end if
      end do

      !$acc loop vector private(L)
      do i = 1, i1
         L = mf + (i - 1)*4
         !$acc loop seq
         do j = 1, myhalf
            jj = my - j + 1
            if (mf .le. mtrundef(j)) then

               wss(L, 1) = wss(L, 1) + fj_polyw(j, L)*wccSUM(j, 1)
               wss(L + 1, 1) = wss(L + 1, 1) + fj_polyw(j, L + 1)*wccDIF(j, 1)
               wss(L + 2, 1) = wss(L + 2, 1) + fj_polyw(j, L + 2)*wccSUM(j, 1)
               wss(L + 3, 1) = wss(L + 3, 1) + fj_polyw(j, L + 3)*wccDIF(j, 1)
               wss(L, 2) = wss(L, 2) + fj_polyw(j, L)*wccSUM(j, 2)
               wss(L + 1, 2) = wss(L + 1, 2) + fj_polyw(j, L + 1)*wccDIF(j, 2)
               wss(L + 2, 2) = wss(L + 2, 2) + fj_polyw(j, L + 2)*wccSUM(j, 2)
               wss(L + 3, 2) = wss(L + 3, 2) + fj_polyw(j, L + 3)*wccDIF(j, 2)
            end if
         end do
      end do
      !$acc loop vector private(L)
      do i = 1, i2
         L = mf + i1*4 + (i - 1)*2
         !$acc loop seq
         do j = 1, myhalf
            jj = my - j + 1
            if (mf .le. mtrundef(j)) then
               wss(L, 1) = wss(L, 1) + fj_polyw(j, L)*wccSUM(j, 1)
               wss(L + 1, 1) = wss(L + 1, 1) + fj_polyw(j, L + 1)*wccDIF(j, 1)
               wss(L, 2) = wss(L, 2) + fj_polyw(j, L)*wccSUM(j, 2)
               wss(L + 1, 2) = wss(L + 1, 2) + fj_polyw(j, L + 1)*wccDIF(j, 2)
            end if
         end do
      end do
      !$acc loop vector private(L)
      do i = 1, i3
         L = jtrun
         !$acc loop seq
         do j = 1, myhalf
            jj = my - j + 1
            if (mf .le. mtrundef(j)) then
               wss(L, 1) = wss(L, 1) + fj_polyw(j, L)*wccSUM(j, 1)
               wss(L, 2) = wss(L, 2) + fj_polyw(j, L)*wccSUM(j, 2)
            end if
         end do
      end do
      !$acc loop vector
      do L = mf, jtrun
         s(L, m, 1) = wss(L, 1)
         s(L, m, 2) = wss(L, 2)
      end do

   end do
   !$acc exit data delete(wcc_fk, wss, fj_polyw, wccSUM, wccDIF) async(async_id)

   return
end
