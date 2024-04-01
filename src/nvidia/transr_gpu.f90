!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine transr_gpu(jtrun, jtmax, nx, my, my_max, lev, poly, wss &
                      , cc, num, nsize)
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
!  lev: number of levels to transform
!  poly: legendre polynomials
!  wss: spectral coefficient array to transform
!  num: number of variables grouped together
!
! *** output ***
!
!  cc_r8: 3-d output grid point fields
!
!  **************************************
!
   use const, only: RTYPE
   use index
   use fftcom

   implicit none

   integer jtrun, jtmax, nx, my, my_max, lev, num, nsize
   integer mlx, myhalf, lev2, m, mf, lmax, lchk, j, k, l
   integer nb, jchk, jje, jlistnum_fj, j_fj, l_fj, llistnum_fj
   integer kk, ll, jj, jx, j2, ii, i, jtrunj, mchk, mm, mp, mlst
   integer mm1, mp1, mlst1, mm2, mp2, mlst2, mm3, mp3, mlst3, nxj

   real(kind=RTYPE) sa00, sa10, sb00, sb10

   real(kind=RTYPE) poly(jtrun, my/2, jtmax)
   real(kind=RTYPE) cc(nx + 2, lev, num, my_max), wss(lev, 2, num, jtrun, jtmax)

   real(kind=RTYPE) gwk1(nx + 2, lev, num, my_max)

   real(kind=RTYPE) wcc_fk(lev, 2, num, jtmax, my_max*nsize)
   real(kind=RTYPE) twcc_fk(lev, 2, num, jtmax*nsize, my_max)
   real(kind=RTYPE) tcc(lev, 2, num, my/2), tc2(lev, 2, num, my/2)
   real ws2(lev, 2, num, jtrun)

   real fj_poly(jtrun, my/2)
   real fj_tcc(lev*2*num, my/2)
   real fj_tc2(lev*2*num, my/2)
   real fj_ws2(lev*2*num, jtrun)
   real fj_wss(lev*2*num, jtrun)
   integer jlist_fj(my/2)
   integer, parameter :: async_id = 1

   !$acc data create(ws2, tcc, tc2, fj_tcc, fj_tc2, fj_wss, fj_ws2, fj_poly, wcc_fk, twcc_fk) copyin(wss, poly, jlist2) copyout(cc) async(async_id)
   !$acc kernels async(async_id)
   cc = 0.
   wcc_fk = 0.
   !$acc end kernels

   mlx = (jtrun/2)*((jtrun + 1)/2)
   myhalf = my/2
   lev2 = lev*2

   do m = 1, mlistnum
      mf = mlist(m)
      lmax = jtrun - mf + 1
      lchk = iand(lmax, 1)

      nb = 32
      jchk = iand(myhalf, 1)
      jje = myhalf - jchk
      !$acc data copy(jlistnum_fj, mtrundef, jlist_fj) async(async_id)
      !$acc kernels async(async_id)
      jlistnum_fj = 0
      !$acc loop seq
      do j = 1, jje
         if (mf .le. mtrundef(j)) then
            jlistnum_fj = jlistnum_fj + 1
            jlist_fj(jlistnum_fj) = j
         end if
      end do
      !$acc end kernels
      !$acc end data

      !$acc data copy(jlist_fj) async(async_id)
      !$acc parallel loop collapse(2) async(async_id)
      do l = mf, jtrun - 1, 2
      do k = 1, lev2*num
         ws2(k, 1, 1, l) = wss(k, 1, 1, l, m)
         ws2(k, 1, 1, l + 1) = -wss(k, 1, 1, l + 1, m)
      end do
      end do
      if (lchk .eq. 1) then
         l = jtrun
         !$acc parallel loop async(async_id)
         do k = 1, lev2*num
            ws2(k, 1, 1, l) = -wss(k, 1, 1, l, m)
         end do
      end if
      !$acc kernels async(async_id)
      !$acc loop
      do k = 1, lev2*num*myhalf
         tcc(k, 1, 1, 1) = 0.
         tc2(k, 1, 1, 1) = 0.
      end do

      fj_tcc = 0.0
      fj_tc2 = 0.0
      fj_wss = 0.0
      fj_ws2 = 0.0
      fj_poly = 0.0

      !$acc loop gang
      do j_fj = 1, jlistnum_fj
         j = jlist_fj(j_fj)
         !$acc loop vector
         do l = mf, jtrun
            l_fj = l - mf + 1
            fj_poly(l_fj, j_fj) = poly(l, j, m)
         end do
      end do
      !$acc loop gang
      do l = mf, jtrun
         l_fj = l - mf + 1
         !$acc loop vector
         do k = 1, lev2*num
            fj_wss(k, l_fj) = wss(k, 1, 1, l, m)
         end do
      end do
      !$acc loop gang
      do l = mf, jtrun
         l_fj = l - mf + 1
         !$acc loop vector
         do k = 1, lev2*num
            fj_ws2(k, l_fj) = ws2(k, 1, 1, l)
         end do
      end do
      !$acc end kernels

      llistnum_fj = jtrun - mf + 1

      call dgemm_async('n', 'n', lev2*num, jlistnum_fj, llistnum_fj &
                       , 1.0d+0, fj_wss, lev2*num, fj_poly, jtrun, 1.0d+0, fj_tcc, lev2*num, async_id)

      !$acc parallel loop gang async(async_id)
      do j_fj = 1, jlistnum_fj
         j = jlist_fj(j_fj)
         !$acc loop vector
         do k = 1, lev2*num
            tcc(k, 1, 1, j) = tcc(k, 1, 1, j) + fj_tcc(k, j_fj)
         end do
      end do

      call dgemm_async('n', 'n', lev2*num, jlistnum_fj, llistnum_fj &
                       , 1.0d+0, fj_ws2, lev2*num, fj_poly, jtrun, 1.0d+0, fj_tc2, lev2*num, async_id)

      !$acc parallel loop gang async(async_id)
      do j_fj = 1, jlistnum_fj
         j = jlist_fj(j_fj)
         !$acc loop vector
         do k = 1, lev2*num
            tc2(k, 1, 1, j) = tc2(k, 1, 1, j) + fj_tc2(k, j_fj)
         end do
      end do

!
! odd number
!
      if (jchk .eq. 1) then
         j = myhalf
         if (mf .le. mtrundef(j)) then
            !$acc parallel loop gang collapse(2) async(async_id)
            do kk = 1, lev2*num, nb
            do ll = mf, jtrun, nb
               !$acc loop vector
               do k = kk, min(kk + nb - 1, lev2*num), 2
                  sa00 = tcc(k, 1, 1, j)
                  sa10 = tcc(k + 1, 1, 1, j)
                  sb00 = tc2(k, 1, 1, j)
                  sb10 = tc2(k + 1, 1, 1, j)
                  do l = ll, min(ll + nb - 1, jtrun)
                     sa00 = sa00 + poly(l, j, m)*wss(k, 1, 1, l, m)
                     sa10 = sa10 + poly(l, j, m)*wss(k + 1, 1, 1, l, m)
                     sb00 = sb00 + poly(l, j, m)*ws2(k, 1, 1, l)
                     sb10 = sb10 + poly(l, j, m)*ws2(k + 1, 1, 1, l)
                  end do
                  tcc(k, 1, 1, j) = sa00
                  tcc(k + 1, 1, 1, j) = sa10
                  tc2(k, 1, 1, j) = sb00
                  tc2(k + 1, 1, 1, j) = sb10
               end do
            end do
            end do
         end if
      end if

      !$acc parallel loop gang async(async_id)
      do j = 1, myhalf
         jj = jlist2(j)
         jx = my - j + 1
         j2 = jlist2(jx)
         !$acc loop vector
         do k = 1, lev2*num
            wcc_fk(k, 1, 1, m, jj) = tcc(k, 1, 1, j)
            wcc_fk(k, 1, 1, m, j2) = tc2(k, 1, 1, j)
         end do
      end do
      !$acc end data
      !$acc wait(async_id)
   end do

   call mpe_transpose_sr_sp_async(wcc_fk, twcc_fk, lev*2*num, jtmax, my_max, nsize, col_comm, async_id)
!      call mpe_transpose_sr(wcc_fk,twcc_fk,lev*2*num,jtmax,my_max,nsize,col_comm)

   !$acc data copy(jlist1, mtrundef, nlist) async(async_id)
   !$acc parallel loop gang async(async_id)
   do jj = 1, jlistnum

      !$acc loop vector collapse(3)
      do ii = 1, num
      do k = 1, lev
      do i = 1, nx + 2
         cc(i, k, ii, jj) = 0.
      end do
      end do
      end do

      j = jlist1(jj)
      jtrunj = mtrundef(j)
      mchk = iand(jtrunj, 3)

      !$acc loop worker
      do m = 1, mchk
         mm = 2*m - 1
         mp = mm + 1
         mlst = nlist(m)
         !$acc loop vector collapse(2)
         do ii = 1, num
         do k = 1, lev
            cc(mm, k, ii, jj) = twcc_fk(k, 1, ii, mlst, jj)
            cc(mp, k, ii, jj) = twcc_fk(k, 2, ii, mlst, jj)
         end do
         end do
      end do
      !$acc loop worker
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
            cc(mm, k, ii, jj) = twcc_fk(k, 1, ii, mlst, jj)
            cc(mp, k, ii, jj) = twcc_fk(k, 2, ii, mlst, jj)
            cc(mm1, k, ii, jj) = twcc_fk(k, 1, ii, mlst1, jj)
            cc(mp1, k, ii, jj) = twcc_fk(k, 2, ii, mlst1, jj)
            cc(mm2, k, ii, jj) = twcc_fk(k, 1, ii, mlst2, jj)
            cc(mp2, k, ii, jj) = twcc_fk(k, 2, ii, mlst2, jj)
            cc(mm3, k, ii, jj) = twcc_fk(k, 1, ii, mlst3, jj)
            cc(mp3, k, ii, jj) = twcc_fk(k, 2, ii, mlst3, jj)
         end do
         end do
      end do

   end do
   !$acc end data
   !$acc wait(async_id)

!
   if (length_fft .eq. 0 .and. lreduce .eq. 0) then
#ifdef SP
      call rfftmlt_sp(cc, gwk1, trigs, ifax, 1, nx + 2, nx, lev*jlistnum*num, 1)
#else
      call rfftmlt(cc, gwk1, trigs, ifax, 1, nx + 2, nx, lev*jlistnum*num, 1)
#endif
   else
#ifdef SP
!$omp  parallel do default(none)                                      &
!$omp  private(jj,j,nxj,gwk1)                                         &
!$omp  shared(jlistnum,jlist1,nxdef,cc,trigsj,ifaxj,nx,lev,num)       &
!$omp  schedule(dynamic)
      do jj = 1, jlistnum
         j = jlist1(jj)
         nxj = nxdef(j)
         call rfftmlt_sp(cc(1, 1, 1, jj), gwk1(1, 1, 1, jj), trigsj(1, j), ifaxj(1, j), &
      end do
!$omp end parallel do
#else
      !$acc data create(gwk1)
      call rfftmlt_loop(cc, gwk1, trigsj, ifaxj, jlist1, nxdef, jlistnum, nx + 2, lev*num, 1)
      !$acc end data
#endif
   end if
   !$acc end data
!
20       continue

!CWB2021

         return
      end
