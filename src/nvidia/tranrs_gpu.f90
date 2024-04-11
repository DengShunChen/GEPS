subroutine tranrs_gpu(jtrun, jtmax, nx, my, my_max, lev, poly, w, cc &
                      , wss, num, nsize)
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
   use fj_pad

   implicit none

   integer jtrun, jtmax, nx, my, my_max, lev, num, nsize
   integer mlx, myhalf, lev2, nxj, j, jj, mchk, jtrunj, m, mm, mp
   integer mlst, ii, k, mm1, mp1, mlst1, mm2, mp2, mlst2, mm3, mp3, mlst3, mf
   integer lchk, lle, jlistnum_fj, j_fj, j1, j2, l, l_fj, llistnum_fj

   real(kind=RTYPE) poly(jtrun, my/2, jtmax), w(my) ! present on device
   real(kind=RTYPE) wss(lev, 2, num, jtrun, jtmax) ! present on device
!
   real(kind=RTYPE) gwk1(nx + 2, lev, num, my_max)
!
   real(kind=RTYPE) wcc_fk(lev, 2, num, jtmax, my_max*nsize)
   real(kind=RTYPE) twcc_fk(lev, 2, num, jtmax*nsize, my_max)
   real(kind=RTYPE) cc(nx + 2, lev, num, my_max) ! present on device
!
   real fj_wss_sum(lev*2*num, jtrun)
   real fj_wss_dif(lev*2*num, jtrun)
   real fj_polyw_sum(my/2 + npad, jtrun)
   real fj_polyw_dif(my/2 + npad, jtrun)
   real fj_wccSUM(lev*2*num, my/2)
   real fj_wccDIF(lev*2*num, my/2)
   integer jlist_fj(my/2, mlistnum)
   integer jlistnum_fj_array(mlistnum)
   integer async_id

   async_id = 1

   !$acc enter data copyin(nlist, jlist2) async(async_id)
   !$acc enter data create(twcc_fk, gwk1, wcc_fk) async(async_id)
   !$acc kernels async(async_id)
   twcc_fk = 0.
   wss = 0.
   gwk1 = 0.
   !$acc end kernels

   mlx = (jtrun/2)*((jtrun + 1)/2)
   myhalf = my/2
   lev2 = lev*2
!
!  fft for each guassian latitude of 2-d field
!
   if (length_fft .eq. 0 .and. lreduce .eq. 0) then
#ifdef SP
      call rfftmlt_sp(cc, gwk1, trigs, ifax, 1, nx + 2, nx, lev*jlistnum*num, -1)
#else
      call rfftmlt(cc, gwk1, trigs, ifax, 1, nx + 2, nx, lev*jlistnum*num, -1)
#endif
   else
#ifdef SP
!$omp  parallel do default(none)                                &
!$omp  private(jj,j,nxj,gwk1)                                   &
!$omp  shared(jlistnum,jlist1,nxdef,cc,trigsj,ifaxj,nx,lev,num) &
!$omp  schedule(dynamic)
      do jj = 1, jlistnum
         j = jlist1(jj)
         nxj = nxdef(j)
         call rfftmlt_sp(cc(1, 1, 1, jj), gwk1(1, 1, 1, jj), trigsj(1, j), ifaxj(1, j), &
                         end do
!$omp end parallel do
#else
         call rfftmlt_loop_identical(cc, gwk1, trigsj, ifaxj, jlist1, nxdef, jlistnum, nx + 2, lev*num, -1)
#endif
         end if
         !$acc exit data delete(gwk1) async(async_id)

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
         !$acc exit data delete(cc) async(async_id)

         call mpe_transpose_rs_sp_gpu(twcc_fk, wcc_fk, lev*2*num, jtmax, my_max, nsize, col_comm)
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
!
         end do
         !$acc exit data delete(jlist_fj, wcc_fk, poly, w, fj_wss_sum, fj_wss_dif, fj_polyw_sum, fj_polyw_dif, fj_wccSUM, fj_wccDIF) async(async_id)
         !$acc exit data delete(nlist, jlist2) async(async_id)

         return
      end
!
      subroutine tranrs_dp_gpu(jtrun, jtmax, nx, my, my_max, lev, poly, w, cc &
                               , wss, num, nsize)
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
!     use paramt
         use fftcom
         use fj_pad
!
         implicit none

         integer jtrun, jtmax, nx, my, my_max, lev, num, nsize
         integer mlx, myhalf, lev2, nxj, j, jj, mchk, jtrunj, m, mm, mp
         integer mlst, ii, k, mm1, mp1, mlst1, mm2, mp2, mlst2, mm3, mp3, mlst3, mf
         integer lchk, lle, jlistnum_fj, j_fj, j1, j2, l, l_fj, llistnum_fj

         real(kind=RTYPE) poly(jtrun, my/2, jtmax), w(my)
         real wss(lev, 2, num, jtrun, jtmax)
!
         real gwk1(nx + 2, lev, num, my_max)
!
         real wcc_fk(lev, 2, num, jtmax, my_max*nsize)
         real twcc_fk(lev, 2, num, jtmax*nsize, my_max)
         real(kind=RTYPE) cc(nx + 2, lev, num, my_max)
         real cc8(nx + 2, lev, num, my_max)
!
         real fj_wss_sum(lev*2*num, jtrun)
         real fj_wss_dif(lev*2*num, jtrun)
         real fj_polyw_sum(my/2 + npad, jtrun)
         real fj_polyw_dif(my/2 + npad, jtrun)
         real fj_wccSUM(lev*2*num, my/2)
         real fj_wccDIF(lev*2*num, my/2)
         integer jlist_fj(my/2)
!
!CWBinit
         twcc_fk = 0.
         wss = 0.

!CWB2014
         gwk1 = 0.

!CWB2021 for single precision test

         mlx = (jtrun/2)*((jtrun + 1)/2)
         myhalf = my/2
         lev2 = lev*2
         cc8 = cc
!
!  fft for each guassian latitude of 2-d field
!
         if (length_fft .eq. 0 .and. lreduce .eq. 0) then
!#ifdef SP
!      call rfftmlt_sp(cc,gwk1,trigs,ifax,1,nx+2,nx,lev*jlistnum*num,-1)
!#else
            call rfftmlt(cc8, gwk1, trigs, ifax, 1, nx + 2, nx, lev*jlistnum*num, -1)
!#endif
         else
!$omp  parallel do default(none)                                &
!$omp  private(jj,j,nxj,gwk1)                                   &
!$omp  shared(jlistnum,jlist1,nxdef,cc8,trigsj,ifaxj,nx,lev,num) &
!$omp  schedule(dynamic)
            do jj = 1, jlistnum
               j = jlist1(jj)
               nxj = nxdef(j)
!#ifdef SP
!        call rfftmlt_sp(cc(1,1,1,jj),gwk1(1,1,1,jj),trigsj(1,j),ifaxj(1,j), &
!#else
               call rfftmlt(cc8(1, 1, 1, jj), gwk1(1, 1, 1, jj), trigsj(1, j), ifaxj(1, j), &
                            !#endif
                            1, nx + 2, nxj, lev*num, -1)
            end do
!$omp end parallel do
         end if
!
         mchk = iand(jtrun, 3)

         do j = 1, jlistnum
            jj = jlist1(j)
            jtrunj = mtrundef(jj)
            mchk = iand(jtrunj, 3)

            do m = 1, mchk
               mm = 2*m - 1
               mp = mm + 1
               mlst = nlist(m)
               do ii = 1, num
                  do k = 1, lev
                     twcc_fk(k, 1, ii, mlst, j) = cc8(mm, k, ii, j)
                     twcc_fk(k, 2, ii, mlst, j) = cc8(mp, k, ii, j)
                  end do
               end do
            end do

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
               do ii = 1, num
                  do k = 1, lev
                     twcc_fk(k, 1, ii, mlst, j) = cc8(mm, k, ii, j)
                     twcc_fk(k, 2, ii, mlst, j) = cc8(mp, k, ii, j)
                     twcc_fk(k, 1, ii, mlst1, j) = cc8(mm1, k, ii, j)
                     twcc_fk(k, 2, ii, mlst1, j) = cc8(mp1, k, ii, j)
                     twcc_fk(k, 1, ii, mlst2, j) = cc8(mm2, k, ii, j)
                     twcc_fk(k, 2, ii, mlst2, j) = cc8(mp2, k, ii, j)
                     twcc_fk(k, 1, ii, mlst3, j) = cc8(mm3, k, ii, j)
                     twcc_fk(k, 2, ii, mlst3, j) = cc8(mp3, k, ii, j)
                  end do
               end do
            end do

         end do

!      call mpe_transpose_rs_sp(twcc_fk,wcc_fk,lev*2*num,jtmax,my_max,nsize,col_comm)
         call mpe_transpose_rs(twcc_fk, wcc_fk, lev*2*num, jtmax, my_max, nsize, col_comm)

         do m = 1, mlistnum
            mf = mlist(m)

            fj_wss_sum = 0.0
            fj_wss_dif = 0.0
            fj_polyw_sum = 0.0
            fj_polyw_dif = 0.0
            fj_wccSUM = 0.0
            fj_wccDIF = 0.0

            lchk = iand(jtrun - mf + 1, 1)
            lle = jtrun - lchk

            jlistnum_fj = 0
            do j = 1, myhalf
               if (mf .le. mtrundef(j)) then
                  jlistnum_fj = jlistnum_fj + 1
                  jlist_fj(jlistnum_fj) = j
               end if
            end do

            do j_fj = 1, jlistnum_fj
               j = jlist_fj(j_fj)
               j1 = jlist2(j)
               j2 = jlist2(my - j + 1)
               do k = 1, lev2*num
                  fj_wccSUM(k, j_fj) = (wcc_fk(k, 1, 1, m, j1) + wcc_fk(k, 1, 1, m, j2))
                  fj_wccDIF(k, j_fj) = (wcc_fk(k, 1, 1, m, j1) - wcc_fk(k, 1, 1, m, j2))
               end do
            end do

            if (jtrun .gt. mf) then

               do l = mf, lle, 2
                  l_fj = (l - mf + 2)/2
                  do j_fj = 1, jlistnum_fj
                     j = jlist_fj(j_fj)
                     fj_polyw_sum(j_fj, l_fj) = poly(l, j, m)*w(j)
                     fj_polyw_dif(j_fj, l_fj) = poly(l + 1, j, m)*w(j)
                  end do
               end do

               llistnum_fj = l_fj

               call dgemm('n', 'n', lev2*num, llistnum_fj, jlistnum_fj, 1.0d+0, fj_wccSUM, lev*2*num, &
                          fj_polyw_sum, my/2 + npad, 1.0d+0, fj_wss_sum, lev*2*num)

               do l = mf, lle, 2
                  l_fj = (l - mf + 2)/2
                  do k = 1, lev2*num
                     wss(k, 1, 1, l, m) = wss(k, 1, 1, l, m) + fj_wss_sum(k, l_fj)
                  end do
               end do

               call dgemm('n', 'n', lev2*num, llistnum_fj, jlistnum_fj, 1.0d+0, fj_wccDIF, lev*2*num, &
                          fj_polyw_dif, my/2 + npad, 1.0d+0, fj_wss_dif, lev*2*num)

               do l = mf, lle, 2
                  l_fj = (l - mf + 2)/2
                  do k = 1, lev2*num
                     wss(k, 1, 1, l + 1, m) = wss(k, 1, 1, l + 1, m) + fj_wss_dif(k, l_fj)
                  end do
               end do

            end if

            if (lchk .eq. 1) then
               l = jtrun
               do j_fj = 1, jlistnum_fj
                  j = jlist_fj(j_fj)
                  fj_polyw_sum(j_fj, l) = poly(l, j, m)*w(j)
               end do

               do j_fj = 1, jlistnum_fj
                  do k = 1, lev2*num
                     wss(k, 1, 1, l, m) = wss(k, 1, 1, l, m) + fj_polyw_sum(j_fj, l)*fj_wccSUM(k, j_fj)
                  end do
               end do

            end if
!
         end do

         return
      end
