subroutine tranuv_gpu(jtrun, jtmax, nx, my, my_max, lev, onocos, wcfac &
                      , wdfac, poly, dpoly, vor, div, ut, vt, nsize)
! Present on device: onocos, wcfac, wdfac, poly, dpoly, vor, div, ut, vt, jlist1, nlist, mtrundef, mlist, jlist2
!  subroutine to transform vorticity and divergence to velocity
!  components
!
! *** input ***
!
!  jtrun: zonal wavenumber truncation limit
!  jtmax: maximum amount of zonal waves located in each pe
!  nx: e-w dimension no.
!  my: n-s dimension no.
!  my_max: maximum amount of n-s grids located in each pe
!  lev: number of veritical levels to transform
!  onocos: 1.0/(cos(lat)**2)
!  wcfac: constants defined in cons
!  wdfac: constants defined in cons
!  poly: legendre polynomials
!  dpoly: d(poly)/d(sin(lat))
!  vor: spectral vorticity
!  div: spectral divergence
!
! *** output ***
!
!  ut: e-w velocity component
!  vt: n-s velocity component
!
!  ****************************************
!
   use const, only: RTYPE
   use index
!     use paramt
   use fftcom
   use openacc
   use cudafor

   implicit none

   integer jtrun, jtmax, nx, my, my_max, lev, nsize
   integer myhalf, lev2, mlx, j, m, mf, l, k, nb, jchk
   integer jje, jlistnum_fj, l_fj, j_fj, llistnum_fj
   integer kk, ll, jj, jx, j2, i, jtrunj, mchk, mm, mp, mlst
   integer mm1, mp1, mlst1, mm2, mp2, mlst2, mm3, mp3, mlst3, nxj, ierr

   real(kind=RTYPE) sa00, sa10, sa20, sa30, dummy

   real(kind=RTYPE) poly(jtrun, my/2, jtmax), dpoly(jtrun, my/2, jtmax), &
      onocos(my)
   real(kind=RTYPE) vor(lev, 2, jtrun, jtmax), div(lev, 2, jtrun, jtmax)
   real(kind=RTYPE) ut(nxp, levF, my_max), vt(nxp, levF, my_max)
   real(kind=RTYPE) gwk1(nx + 2, lev, 2, my_max)
   real(kind=RTYPE) wcc_fk(lev, 2, 2, jtmax, my_max*nsize)
   real(kind=RTYPE) twcc_fk(lev, 2, 2, jtmax*nsize, my_max)
   real(kind=RTYPE) cc(nx + 2, lev, 2, my_max)
   real(kind=RTYPE) tcc(lev, 2, 2, my)
   real(kind=RTYPE) ws3(lev, 2, 2, jtrun)
   real(kind=RTYPE) ws4(lev, 2, 2, jtrun)
   real(kind=RTYPE) wcfac(jtrun, jtmax), wdfac(jtrun, jtmax)

   real(kind=RTYPE) tc2(lev, 2, 2, my)
   real(kind=RTYPE) wc(jtrun, my/2), wd(jtrun, my/2)

!CWB2015
!     real      coslr(jm)
!     save coslr
   real(kind=RTYPE), dimension(:), allocatable, save ::  coslr

   logical lfirst
   data lfirst/.true./
   save lfirst

   integer jlist_fj(my/2, mlistnum)
   real fj_ws3(lev*2*2, jtrun)
   real fj_ws4(lev*2*2, jtrun)
   real fj_tcc(lev*2*2, my)
   real fj_wc(jtrun, my/2), fj_wd(jtrun, my/2)
   real fj_tc2(lev*2*2, my)
   integer async_id, istat
   integer(kind=cuda_stream_kind) stream
   integer jlistnum_fj_array(mlistnum)

   async_id = 1
   stream = acc_get_cuda_stream(async_id)

!CWBinit
   wcc_fk = 0.

   myhalf = my/2
   lev2 = lev*2
   mlx = (jtrun/2)*((jtrun + 1)/2)
   nb = 32
   jchk = iand(myhalf, 1)
   jje = myhalf - jchk
   do m = 1, mlistnum
      mf = mlist(m)
      jlistnum_fj = 0
      do j = 1, jje
         if (mf .le. mtrundef(j)) then
            jlistnum_fj = jlistnum_fj + 1
            jlist_fj(jlistnum_fj, m) = j
         end if
      end do
      jlistnum_fj_array(m) = jlistnum_fj
   end do

   if (lfirst) then
!CWB2015 >>>
!       allocate(coslr(jm),stat=ierr)
      allocate (coslr(my), stat=ierr)
      if (ierr /= 0) then
         write (6, *) 'tranuv : allocate fail '
         stop
      end if
!CWB2015 <<<
      do j = 1, my
         coslr(j) = 1./onocos(j)
      end do
      lfirst = .false.
   end if

   !$acc enter data copyin(coslr, jlist_fj) &
   !$acc& create(gwk1, cc, twcc_fk, ws3, ws4, tcc, tc2, wcc_fk, fj_ws3, fj_ws4, fj_tcc, fj_wc, fj_wd, fj_tc2, wc, wd) async(async_id)
   !$acc wait(async_id)
   do m = 1, mlistnum
      mf = mlist(m)
      !$acc parallel loop collapse(2) async(async_id)
      do j = 1, myhalf
         do l = 1, jtrun
            mf = mlist(m)
            if ((mf .le. mtrundef(j)) .AND. (l .ge. mf)) then
               wc(l, j) = wcfac(l, m)*poly(l, j, m)
               wd(l, j) = wdfac(l, m)*dpoly(l, j, m)*coslr(j)
            end if
         end do
      end do
      !$acc parallel loop collapse(2) async(async_id)
      do l = 1, jtrun
         do k = 1, lev
            if (l .ge. mf) then
               ws3(k, 1, 1, l) = +div(k, 2, l, m)
               ws3(k, 2, 1, l) = -div(k, 1, l, m)
               ws3(k, 1, 2, l) = +vor(k, 2, l, m)
               ws3(k, 2, 2, l) = -vor(k, 1, l, m)
               ws4(k, 1, 1, l) = +vor(k, 1, l, m)
               ws4(k, 2, 1, l) = +vor(k, 2, l, m)
               ws4(k, 1, 2, l) = -div(k, 1, l, m)
               ws4(k, 2, 2, l) = -div(k, 2, l, m)
            end if
         end do
      end do
      !$acc host_data use_device(tcc, fj_ws3, fj_ws4, fj_tcc, fj_wc, fj_wd)
      istat = cudaMemsetAsync(tcc, 0.0, size(tcc), stream)
      istat = cudaMemsetAsync(fj_ws3, 0.0, size(fj_ws3), stream)
      istat = cudaMemsetAsync(fj_ws4, 0.0, size(fj_ws4), stream)
      istat = cudaMemsetAsync(fj_tcc, 0.0, size(fj_tcc), stream)
      istat = cudaMemsetAsync(fj_wc, 0.0, size(fj_wc), stream)
      istat = cudaMemsetAsync(fj_wd, 0.0, size(fj_wd), stream)
      !$acc end host_data
      !$acc parallel loop collapse(2) async(async_id)
      do l = 1, jtrun
         do k = 1, lev2*2
            if (l .ge. mf) then
               l_fj = l - mf + 1
               fj_ws3(k, l_fj) = ws3(k, 1, 1, l)
            end if
         end do
      end do
      jlistnum_fj = jlistnum_fj_array(m)
      !$acc parallel loop collapse(2) private(j, l_fj) async(async_id)
      do j_fj = 1, jlistnum_fj
         do l = 1, jtrun
            if (l .ge. mf) then
               j = jlist_fj(j_fj, m)
               l_fj = l - mf + 1
               fj_wc(l_fj, j_fj) = wc(l, j)
            end if
         end do
      end do

      llistnum_fj = jtrun - mf + 1
      call dgemm_async('n', 'n', lev*2*2, jlistnum_fj, llistnum_fj, 1.0d+0, fj_ws3, lev*2*2, fj_wc, &
                       jtrun, 1.0d+0, fj_tcc, lev*2*2, async_id)

      !$acc parallel loop collapse(2) private(j) async(async_id)
      do j_fj = 1, jlistnum_fj
         do k = 1, lev2*2
            j = jlist_fj(j_fj, m)
            tcc(k, 1, 1, j) = tcc(k, 1, 1, j) + fj_tcc(k, j_fj)
         end do
      end do

      !$acc parallel loop collapse(2) private(l_fj) async(async_id)
      do l = 1, jtrun
         do k = 1, lev2*2
            if (l .ge. mf) then
               l_fj = l - mf + 1
               fj_ws4(k, l_fj) = ws4(k, 1, 1, l)
            end if
         end do
      end do

      !$acc parallel loop collapse(2) async(async_id)
      do j_fj = 1, jlistnum_fj
         do l = 1, jtrun
            if (l .ge. mf) then
               j = jlist_fj(j_fj, m)
               l_fj = l - mf + 1
               fj_wd(l_fj, j_fj) = wd(l, j)
            end if
         end do
      end do
      !$acc host_data use_device(fj_tcc)
      istat = cudaMemsetAsync(fj_tcc, 0.0, size(fj_tcc), stream)
      !$acc end host_data

      llistnum_fj = jtrun - mf + 1
      call dgemm_async('n', 'n', lev*2*2, jlistnum_fj, llistnum_fj, 1.0d+0, fj_ws4, lev*2*2, fj_wd, &
                       jtrun, 1.0d+0, fj_tcc, lev*2*2, async_id)

      !$acc parallel loop collapse(2) private(j) async(1)
      do j_fj = 1, jlistnum_fj
         do k = 1, lev2*2
            j = jlist_fj(j_fj, m)
            tcc(k, 1, 1, j) = tcc(k, 1, 1, j) + fj_tcc(k, j_fj)
         end do
      end do
!
! odd number
!
      if (jchk .eq. 1) then
         j = myhalf
         if (mf .le. mtrundef(j)) then
            !$acc parallel loop collapse(2) private(sa00, sa10, sa20, sa30) async(async_id)
            do kk = 1, lev2*2, nb
            do ll = mf, jtrun, nb
               do k = kk, min(kk + nb - 1, lev2*2), 4
                  sa00 = tcc(k, 1, 1, j)
                  sa10 = tcc(k + 1, 1, 1, j)
                  sa20 = tcc(k + 2, 1, 1, j)
                  sa30 = tcc(k + 3, 1, 1, j)
                  do l = ll, min(ll + nb - 1, jtrun)
                     sa00 = sa00 + ws3(k, 1, 1, l)*wc(l, j) + ws4(k, 1, 1, l)*wd(l, j)
                     sa10 = sa10 + ws3(k + 1, 1, 1, l)*wc(l, j) + ws4(k + 1, 1, 1, l)*wd(l, j)
                     sa20 = sa20 + ws3(k + 2, 1, 1, l)*wc(l, j) + ws4(k + 2, 1, 1, l)*wd(l, j)
                     sa30 = sa30 + ws3(k + 3, 1, 1, l)*wc(l, j) + ws4(k + 3, 1, 1, l)*wd(l, j)
                  end do
                  tcc(k, 1, 1, j) = sa00
                  tcc(k + 1, 1, 1, j) = sa10
                  tcc(k + 2, 1, 1, j) = sa20
                  tcc(k + 3, 1, 1, j) = sa30
               end do
            end do
            end do
         end if
      end if

      !$acc parallel loop collapse(2) async(async_id)
      do l = mf, jtrun, 2
      do k = 1, lev
         ws3(k, 1, 1, l) = +div(k, 2, l, m)
         ws3(k, 2, 1, l) = -div(k, 1, l, m)
         ws3(k, 1, 2, l) = +vor(k, 2, l, m)
         ws3(k, 2, 2, l) = -vor(k, 1, l, m)
         ws4(k, 1, 1, l) = -vor(k, 1, l, m)
         ws4(k, 2, 1, l) = -vor(k, 2, l, m)
         ws4(k, 1, 2, l) = +div(k, 1, l, m)
         ws4(k, 2, 2, l) = +div(k, 2, l, m)
      end do
      end do

      !$acc parallel loop collapse(2) async(async_id)
      do l = mf + 1, jtrun, 2
         do k = 1, lev
            ws3(k, 1, 1, l) = -div(k, 2, l, m)
            ws3(k, 2, 1, l) = +div(k, 1, l, m)
            ws3(k, 1, 2, l) = -vor(k, 2, l, m)
            ws3(k, 2, 2, l) = +vor(k, 1, l, m)
            ws4(k, 1, 1, l) = +vor(k, 1, l, m)
            ws4(k, 2, 1, l) = +vor(k, 2, l, m)
            ws4(k, 1, 2, l) = -div(k, 1, l, m)
            ws4(k, 2, 2, l) = -div(k, 2, l, m)
         end do
      end do

      !$acc host_data use_device(tc2, fj_ws3, fj_ws4, fj_tc2, fj_wc, fj_wd)
      istat = cudaMemsetAsync(tc2, 0.0, size(tc2), stream)
      istat = cudaMemsetAsync(fj_ws3, 0.0, size(fj_ws3), stream)
      istat = cudaMemsetAsync(fj_ws4, 0.0, size(fj_ws4), stream)
      istat = cudaMemsetAsync(fj_tc2, 0.0, size(fj_tc2), stream)
      istat = cudaMemsetAsync(fj_wc, 0.0, size(fj_wc), stream)
      istat = cudaMemsetAsync(fj_wd, 0.0, size(fj_wd), stream)
      !$acc end host_data

      !$acc parallel loop collapse(2) private(l_fj) async(async_id)
      do l = 1, jtrun
         do k = 1, lev2*2
            if (l .ge. mf) then
               l_fj = l - mf + 1
               fj_ws3(k, l_fj) = ws3(k, 1, 1, l)
            end if
         end do
      end do

      !$acc parallel loop collapse(2) private(j, l_fj) async(async_id)
      do j_fj = 1, jlistnum_fj
         do l = 1, jtrun
            if (l .ge. mf) then
               j = jlist_fj(j_fj, m)
               l_fj = l - mf + 1
               fj_wc(l_fj, j_fj) = wc(l, j)
            end if
         end do
      end do

      llistnum_fj = jtrun - mf + 1

      call dgemm_async('n', 'n', lev*2*2, jlistnum_fj, llistnum_fj, 1.0d+0, fj_ws3, lev*2*2, fj_wc, &
                       jtrun, 1.0d+0, fj_tc2, lev*2*2, async_id)

      !$acc parallel loop collapse(2) private(j) async(async_id)
      do j_fj = 1, jlistnum_fj
         do k = 1, lev2*2
            j = jlist_fj(j_fj, m)
            tc2(k, 1, 1, j) = tc2(k, 1, 1, j) + fj_tc2(k, j_fj)
         end do
      end do

      !$acc parallel loop collapse(2) private(l_fj) async(async_id)
      do l = 1, jtrun
         do k = 1, lev2*2
            if (l .ge. mf) then
               l_fj = l - mf + 1
               fj_ws4(k, l_fj) = ws4(k, 1, 1, l)
            end if
         end do
      end do

      !$acc parallel loop collapse(2) private(j, l_fj) async(async_id)
      do j_fj = 1, jlistnum_fj
         do l = 1, jtrun
            if (l .ge. mf) then
               j = jlist_fj(j_fj, m)
               l_fj = l - mf + 1
               fj_wd(l_fj, j_fj) = wd(l, j)
            end if
         end do
      end do

      !$acc host_data use_device(fj_tc2)
      istat = cudaMemsetAsync(fj_tc2, 0.0, size(fj_tc2), stream)
      !$acc end host_data

      llistnum_fj = jtrun - mf + 1
      call dgemm_async('n', 'n', lev*2*2, jlistnum_fj, llistnum_fj, 1.0d+0, fj_ws4, lev*2*2, fj_wd, &
                       jtrun, 1.0d+0, fj_tc2, lev*2*2, async_id)

      !$acc parallel loop collapse(2) private(j) async(async_id)
      do j_fj = 1, jlistnum_fj
         do k = 1, lev2*2
            j = jlist_fj(j_fj, m)
            tc2(k, 1, 1, j) = tc2(k, 1, 1, j) + fj_tc2(k, j_fj)
         end do
      end do

!
! odd number
!
      if (jchk .eq. 1) then
         j = myhalf
         if (mf .le. mtrundef(j)) then
            !$acc parallel loop collapse(2) private(sa00, sa10, sa20, sa30) async(async_id)
            do kk = 1, lev2*2, nb
            do ll = mf, jtrun, nb
               do k = kk, min(kk + nb - 1, lev2*2), 4
                  sa00 = tc2(k, 1, 1, j)
                  sa10 = tc2(k + 1, 1, 1, j)
                  sa20 = tc2(k + 2, 1, 1, j)
                  sa30 = tc2(k + 3, 1, 1, j)
                  do l = ll, min(ll + nb - 1, jtrun)
                     sa00 = sa00 + ws3(k, 1, 1, l)*wc(l, j) + ws4(k, 1, 1, l)*wd(l, j)
                     sa10 = sa10 + ws3(k + 1, 1, 1, l)*wc(l, j) + ws4(k + 1, 1, 1, l)*wd(l, j)
                     sa20 = sa20 + ws3(k + 2, 1, 1, l)*wc(l, j) + ws4(k + 2, 1, 1, l)*wd(l, j)
                     sa30 = sa30 + ws3(k + 3, 1, 1, l)*wc(l, j) + ws4(k + 3, 1, 1, l)*wd(l, j)
                  end do
                  tc2(k, 1, 1, j) = sa00
                  tc2(k + 1, 1, 1, j) = sa10
                  tc2(k + 2, 1, 1, j) = sa20
                  tc2(k + 3, 1, 1, j) = sa30
               end do
            end do
            end do
         end if
      end if

      !$acc parallel loop collapse(2) private(jj, jx, j2) async(async_id)
      do j = 1, myhalf
         do k = 1, lev*2*2
            jj = jlist2(j)
            jx = my - j + 1
            j2 = jlist2(jx)
            wcc_fk(k, 1, 1, m, jj) = tcc(k, 1, 1, j)
            wcc_fk(k, 1, 1, m, j2) = tc2(k, 1, 1, j)
         end do
      end do
   end do   ! end of big m loop

   ! Present on device: wcc_fk, twcc_fk
   call mpe_transpose_sr_sp_gpu(wcc_fk, twcc_fk, lev*2*2, jtmax, my_max, nsize, col_comm)

   do jj = 1, jlistnum
      !$acc kernels present(cc) async(async_id)
      do i = 1, (nx + 2)*lev*2
         cc(i, 1, 1, jj) = 0.
      end do
      !$acc end kernels
!
      j = jlist1(jj)
      jtrunj = mtrundef(j)
      mchk = iand(jtrunj, 3)

      !$acc parallel loop gang async(async_id)
      do m = 1, mchk
         mm = 2*m - 1
         mp = mm + 1
         mlst = nlist(m)
         !$acc loop vector
         do k = 1, lev
            cc(mm, k, 1, jj) = twcc_fk(k, 1, 1, mlst, jj)
            cc(mp, k, 1, jj) = twcc_fk(k, 2, 1, mlst, jj)
            cc(mm, k, 2, jj) = twcc_fk(k, 1, 2, mlst, jj)
            cc(mp, k, 2, jj) = twcc_fk(k, 2, 2, mlst, jj)
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
         !$acc loop vector
         do k = 1, lev
            cc(mm, k, 1, jj) = twcc_fk(k, 1, 1, mlst, jj)
            cc(mp, k, 1, jj) = twcc_fk(k, 2, 1, mlst, jj)
            cc(mm, k, 2, jj) = twcc_fk(k, 1, 2, mlst, jj)
            cc(mp, k, 2, jj) = twcc_fk(k, 2, 2, mlst, jj)
            cc(mm1, k, 1, jj) = twcc_fk(k, 1, 1, mlst1, jj)
            cc(mp1, k, 1, jj) = twcc_fk(k, 2, 1, mlst1, jj)
            cc(mm1, k, 2, jj) = twcc_fk(k, 1, 2, mlst1, jj)
            cc(mp1, k, 2, jj) = twcc_fk(k, 2, 2, mlst1, jj)
            cc(mm2, k, 1, jj) = twcc_fk(k, 1, 1, mlst2, jj)
            cc(mp2, k, 1, jj) = twcc_fk(k, 2, 1, mlst2, jj)
            cc(mm2, k, 2, jj) = twcc_fk(k, 1, 2, mlst2, jj)
            cc(mp2, k, 2, jj) = twcc_fk(k, 2, 2, mlst2, jj)
            cc(mm3, k, 1, jj) = twcc_fk(k, 1, 1, mlst3, jj)
            cc(mp3, k, 1, jj) = twcc_fk(k, 2, 1, mlst3, jj)
            cc(mm3, k, 2, jj) = twcc_fk(k, 1, 2, mlst3, jj)
            cc(mp3, k, 2, jj) = twcc_fk(k, 2, 2, mlst3, jj)
         end do
      end do

   end do

#ifdef SP
   print *, "Symbol SP is not supported."
   call exit(1)
#endif

   if (length_fft .eq. 0 .and. lreduce .eq. 0) then
      call rfftmlt(cc, gwk1, trigs, ifax, 1, nx + 2, nx, lev*jlistnum*2, 1) ! CWB2015
   else
      call rfftmlt_loop_identical(cc, gwk1, trigsj, ifaxj, jlist1, nxdef, jlistnum, nx + 2, lev*2, 1)
   end if
!
!     do 22 jj=1,jlistnum
!       j= jlist1(jj)
!       nxj=nxdef(j)
!     do 22 k=1,lev
!     do 22 i=1,nxj
!     ut(i,k,jj)= cc(i,k,1,jj)
!     vt(i,k,jj)= cc(i,k,2,jj)
!  22 continue

!2dMPI
   call ujoinsr_gpu(cc, ut, vt, dummy, dummy, nx, my_max, levF, jlistnum, 2, 1)
   !$acc exit data delete(coslr, jlist_fj, &
   !$acc& gwk1, cc, twcc_fk, ws3, ws4, tcc, tc2, wcc_fk, fj_ws3, fj_ws4, fj_tcc, fj_wc, fj_wd, fj_tc2, wc, wd) async(async_id)

   return
end
