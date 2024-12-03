program test_cyclic_cell_massadvx
   implicit none

   call mpe_init
   call cons
   call cyclic_cell_massadvx_unit(.true.)
   call cyclic_cell_massadvx_unit(.false.)
   call mpe_finalize
end program test_cyclic_cell_massadvx
! ============================================================
subroutine cyclic_cell_massadvx_unit(forward)
   use const, only: RTYPE, dt, cosl
   use rank, only: myrank
   use param, only: nx, my, my_max, lev, ncld
   use index, only: nxp, levp, levf, myf, jlistnum, jlen, nsizex, row_comm, &
                    nxdef, jlist1, jlist2_2d, nxjlen_all, nxdef
   use grid, only: gglati, fa1, fa2, fa3, fa4
   use mod_ndslfv_monoadv_gpu
   use mpe
   implicit none
   logical forward
   integer, parameter:: steps = 5
   integer, parameter:: nvar = 4

   real(kind=RTYPE):: um(nxp, lev, my_max)
   real(kind=RTYPE):: ut(nxp, lev, my_max), &
                      vt(nxp, lev, my_max), &
                      tt(nxp, lev, my_max), &
                      qt(nxp, lev*ncld, my_max)

   real(kind=RTYPE):: ut_CPU(nxp, lev, my_max), &
                      vt_CPU(nxp, lev, my_max), &
                      tt_CPU(nxp, lev, my_max), &
                      qt_CPU(nxp, lev*ncld, my_max)

   real(kind=RTYPE):: ut_GPU(nxp, lev, my_max), &
                      vt_GPU(nxp, lev, my_max), &
                      tt_GPU(nxp, lev, my_max), &
                      qt_GPU(nxp, lev*ncld, my_max)

   real(kind=RTYPE):: ut_wrk(nxp, lev, my_max), &
                      vt_wrk(nxp, lev, my_max), &
                      tt_wrk(nxp, lev, my_max), &
                      qt_wrk(nxp, lev*ncld, my_max)

   real(kind=RTYPE) dtah, dt_tmp

   integer i, k, j, jj, nxj, n, nn
   real(kind=RTYPE):: err_arr(4, nvar), vamax, ummax, vmmax
   character(len=6):: name(nvar)

   name = (/'t', 'u', 'v', 'q'/)
   dtah = 0.5*dt

   if (myrank .eq. 0) then
      print *, "========================================"
      print *, "    start test cyclic_cell_massadvx_gpu"
      print *, "    with forward = ", forward
      print *, "========================================"
      write (*, '(1X, A5, 1X, f15.2)') "dt=", dt
      write (*, '(1X, A5, 1X, f15.2)') "dtah=", dtah
   end if

   call random_seed()
   call random_number(um)
   call random_number(ut)
   call random_number(vt)
   call random_number(qt)

   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef(j)
      do k = 1, lev
         do i = 1, nxj
            um(i, k, jj) = um(i, k, jj)/dtah*5e-3
            !! um(i, k, jj) = um(i, k, jj)/dtah*0.02
         end do
      end do
   end do

   call gen_cosinebell(tt, 0., 1.)

   ummax = vamax(um, nxp, lev)
   if (myrank .eq. 0) then
      write (*, '(A, 1pe15.7)'), "umamax=", ummax
   end if
   !$acc enter data copyin(nx, my, lev, ncld, nxp, nsizex, &
   !$acc& gglati, fa1, fa2, fa3, fa4, jlist2_2d, nxjlen_all, cosl, nxdef)

   !$acc data create(ut_GPU, vt_GPU, tt_GPU, qt_GPU, &
   !$acc& ut, vt, tt, qt, um)
   do n = 1, steps
      if (myrank .eq. 0) write (*, '("<< ", i3, " >>")') n
      !! << CPU >>
      call massadvx_cpu_wrapper(ut_CPU, vt_CPU, tt_CPU, qt_CPU, &
                                ut, vt, tt, qt, um, dtah, 0, forward, 0)
      !! << GPU >>
      !$acc update device(ut, vt, tt, qt, um)
      call massadvx_gpu_wrapper(ut_GPU, vt_GPU, tt_GPU, qt_GPU, &
                                ut, vt, tt, qt, um, dtah, 0, forward, 1)
      !$acc update self(ut_GPU, vt_GPU, tt_GPU, qt_GPU)
   end do
   !$acc end data

   !$acc exit data delete(nx, my, lev, ncld, nxp, nsizex, &
   !$acc& gglati, fa1, fa2, fa3, fa4, jlist2_2d, nxjlen_all, cosl, nxdef)

   call Varerr(err_arr(1, 1), tt_CPU, nxp, tt_GPU, nxp, lev, 1)
   call Varerr(err_arr(1, 2), ut_CPU, nxp, ut_GPU, nxp, lev, 1)
   call Varerr(err_arr(1, 3), vt_CPU, nxp, vt_GPU, nxp, lev, 1)
   call Varerr(err_arr(1, 4), qt_CPU, nxp, qt_GPU, nxp, lev, ncld)
   if (myrank .eq. 0) then
      do k = 1, nvar
         write (*, '(1X, A6, 1pe23.15, 3(1pe15.7))') &
            name(k), err_arr(1:4, k)
      end do
   end if

   if (all(err_arr(1, 1:nvar) < 1e-10)) then
      if (myrank .eq. 0) write (*, '(A)') &
         "test_cyclic_cell_massadvx passed."
   else
      if (myrank .eq. 0) write (*, '(A)') &
         "test_cyclic_cell_massadvx failed."
      call exit(1)
   end if
end subroutine cyclic_cell_massadvx_unit
! ============================================================
subroutine gen_cosinebell(dat, tim, u0)
   use const, only: RTYPE, sinl, cosl
   use param, only: nx, my_max, my, lev
   use index, only: nxp, nxdef, jlist1, jlistnum
   implicit none

   real(kind=RTYPE), intent(out):: dat(nxp, lev, my_max)
   real(kind=RTYPE), intent(in):: tim, u0

   real(kind=RTYPE):: xlon(nxp, my_max), &
                      xlat(my_max)
   real(kind=RTYPE):: pi, r2d, dx, radius, h0, r0, dist
   integer i, j, k, jj, nxj

   pi = 4.*atan(1.)
   r2d = 180./pi
   radius = 6.37122*1e+6
   r0 = 1./3.
   ! h0 = 1000.
   h0 = 1.
   do jj = 1, jlistnum
      xlat(jj) = asin(sinl(jj))
   end do

   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef(j)
      do i = 1, nxj
         dx = 2.*pi/nxj
         xlon(i, jj) = (i - 1)*dx
      end do
   end do

   dat = 0.
   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef(j)
      do k = 1, lev
         do i = 1, nxj
            dist = acos(cosl(j)*cos(xlon(i, jj) - u0*tim &
                                    - (1.5*pi + (k - 1)/lev*2*pi)))
            if (dist .le. r0) then
               dat(i, k, jj) = 0.5*h0*(1.+cos(pi*dist/r0))
            end if
         end do
      end do
   end do

end subroutine gen_cosinebell
! ============================================================
subroutine massadvx_cpu_wrapper(ut, vt, tt, qt, &
                                ut_sl, vt_sl, tt_sl, qt_sl, &
                                um_sl, deltim, mass, forward, typ)
   use const, only: RTYPE
   use param, only: nx, my, lev, my_max, ncld
   use index
   use grid, only: lonfull, latpart, ndslhvar
   use rank, only: myrank
   implicit none
   integer mass, typ
   real(kind=RTYPE):: deltim
   logical forward
   real(kind=RTYPE), intent(out):: ut(nxp, lev, my_max), &
                                   vt(nxp, lev, my_max), &
                                   tt(nxp, lev, my_max), &
                                   qt(nxp, lev*ncld, my_max)

   real(kind=RTYPE), intent(in):: ut_sl(nxp, lev, my_max), &
                                  vt_sl(nxp, lev, my_max), &
                                  tt_sl(nxp, lev, my_max), &
                                  qt_sl(nxp, lev*ncld, my_max), &
                                  um_sl(nxp, lev, my_max)

   real(kind=RTYPE):: um(nx, levp, my_max), &
                      vdzonl(nx, levp, my_max), &
                      vdmerd(nx, levp, my_max), &
                      ddtemp(nx, levp, my_max), &
                      qvadv(nx, levp, ncld, my_max)

   real(kind=RTYPE):: qqlon(lonfull, levp*ndslhvar, latpart), &
                      uulon(lonfull, levp, latpart)
   integer i, j, k, m, n, ii, jf, nn, kk
   integer lan, lat, lons_lat, ku, kv, kt, kq, kuu, kvv, ktt, kqq

   kuu = 1
   kvv = kuu + levp
   ktt = kvv + levp
   kqq = ktt + levp
   ! ----------------------------------------
   call mpe2d_transpose_ndsl_p2f(um_sl, um, &
                                 nxp, nx, levf, levp, 1, &
                                 myf, my_max, jlistnum, jlen, nsizex, row_comm)

   call mpe2d_transpose_ndsl_p2f(ut_sl, vdzonl, &
                                 nxp, nx, levf, levp, 1, &
                                 myf, my_max, jlistnum, jlen, nsizex, row_comm)

   call mpe2d_transpose_ndsl_p2f(vt_sl, vdmerd, &
                                 nxp, nx, levf, levp, 1, &
                                 myf, my_max, jlistnum, jlen, nsizex, row_comm)

   call mpe2d_transpose_ndsl_p2f(tt_sl, ddtemp, &
                                 nxp, nx, levf, levp, 1, &
                                 myf, my_max, jlistnum, jlen, nsizex, row_comm)

   call mpe2d_transpose_ndsl_p2f(qt_sl, qvadv, &
                                 nxp, nx, levf, levp, ncld, &
                                 myf, my_max, jlistnum, jlen, nsizex, row_comm)

   ! =================================================================
   !   prepare wind and variable in flux form with gaussina weight
   ! =================================================================
   !
   !$omp parallel do &
   !$omp private(lan,lat,lons_lat,i,n,kk,k,kt,kv,ku,kq) &
   !$omp schedule(dynamic)
   do lan = 1, jlistnum
      lat = jlist1(lan)
      lons_lat = nxdef(lat)
      !
      ! u v t at n-1
      !
      do k = 1, levp
         ku = kuu + k - 1
         kv = kvv + k - 1
         kt = ktt + k - 1
         do i = 1, lons_lat
            qqlon(i, ku, lan) = vdzonl(i, k, lan)
            qqlon(i, kv, lan) = vdmerd(i, k, lan)
            qqlon(i, kt, lan) = ddtemp(i, k, lan)
         end do
      end do
      ! rq at n-1
      do n = 1, ncld
         do k = 1, levp
            kq = kqq + k - 1 + (n - 1)*levp
            do i = 1, lons_lat
               qqlon(i, kq, lan) = qvadv(i, k, n, lan)
            end do
         end do
      end do

      ! first set positive advection in east-west direction
      if (forward) then
         call cyclic_cell_massadvxl(lons_lat, lonfull, levp, ndslhvar, deltim, &
                                    um(1, 1, lan), qqlon(1, 1, lan), mass, forward)
      else
         call cyclic_cell_massadvx(lons_lat, lonfull, levp, ndslhvar, deltim, &
                                   um(1, 1, lan), qqlon(1, 1, lan), mass)
      end if

      do k = 1, levp
         ku = kuu + k - 1
         kv = kvv + k - 1
         kt = ktt + k - 1
         do i = 1, lons_lat
            vdzonl(i, k, lan) = qqlon(i, ku, lan)
            vdmerd(i, k, lan) = qqlon(i, kv, lan)
            ddtemp(i, k, lan) = qqlon(i, kt, lan)
         end do
      end do
      ! rq at n-1
      do n = 1, ncld
         do k = 1, levp
            kq = kqq + k - 1 + (n - 1)*levp
            do i = 1, lons_lat
               qvadv(i, k, n, lan) = qqlon(i, kq, lan)
            end do
         end do
      end do
   end do
   !$omp end parallel do

   call mpe2d_transpose_ndsl_f2p(vdzonl, ut, &
                                 nxp, nx, levf, levp, 1, &
                                 myf, my_max, jlistnum, jlen, nsizex, row_comm)

   call mpe2d_transpose_ndsl_f2p(vdmerd, vt, &
                                 nxp, nx, levf, levp, 1, &
                                 myf, my_max, jlistnum, jlen, nsizex, row_comm)

   call mpe2d_transpose_ndsl_f2p(ddtemp, tt, &
                                 nxp, nx, levf, levp, 1, &
                                 myf, my_max, jlistnum, jlen, nsizex, row_comm)

   call mpe2d_transpose_ndsl_f2p(qvadv, qt, &
                                 nxp, nx, levf, levp, ncld, &
                                 myf, my_max, jlistnum, jlen, nsizex, row_comm)

   return
end subroutine massadvx_cpu_wrapper
! ============================================================
subroutine massadvx_gpu_wrapper(ut, vt, tt, qt, &
                                ut_sl, vt_sl, tt_sl, qt_sl, &
                                um_sl, deltim, mass, forward)
   use const, only: RTYPE
   use param, only: nx, my, lev, my_max, ncld
   use index, only: nxp, nsizex, jlist2_2d, nxjlen_all, nxdef
   use grid, only: ndslhvar
   use rank, only: myrank
   use mod_ndslfv_monoadv_gpu
   implicit none
   integer mass
   real(kind=RTYPE):: deltim
   logical forward
   real(kind=RTYPE), intent(out):: ut(nxp, lev, my_max), &
                                   vt(nxp, lev, my_max), &
                                   tt(nxp, lev, my_max), &
                                   qt(nxp, lev*ncld, my_max)

   real(kind=RTYPE), intent(in):: ut_sl(nxp, lev, my_max), &
                                  vt_sl(nxp, lev, my_max), &
                                  tt_sl(nxp, lev, my_max), &
                                  qt_sl(nxp, lev*ncld, my_max), &
                                  um_sl(nxp, lev, my_max)

   real(kind=RTYPE):: um_3df(nx, my, lev)
   integer, parameter:: async_id = -1
   integer i, j, k, m, n, ii, jf, nn, kk
   ! ----------------------------------------
   !$acc host_data use_device(um_sl, ut_sl, vt_sl, tt_sl, qt_sl)
   call advh_gather4GPU_dev(umwrk_d, um_sl, 1, 0)
   ! u v t at n-1
   call advh_gather4GPU_dev(ainp_d(1, 1, 1), ut_sl, 1, 0)
   call advh_gather4GPU_dev(ainp_d(1, 1, 2), vt_sl, 1, 0)
   call advh_gather4GPU_dev(ainp_d(1, 1, 3), tt_sl, 1, 0)
   ! rq at n-1
   call advh_gather4GPU_dev(qpwrk_d, qt_sl, ncld, 0)
   !$acc end host_data
   if (myrank .eq. 0) then
      !$acc data &
      !$acc& create( &
      !$acc& qq_3df(:nx,:my,:lev,:ndslhvar), &
      !$acc& um_3df(:nx,:my,:lev) &
      !$acc& ) &
      !$acc& present(jlist2_2d(nsizex,my), nxjlen_all(nsizex,my) &
      !$acc& ) &
      !$acc& copyin(deltim)
      ! --------------------

      !$acc parallel loop collapse(2) &
      !$acc& deviceptr(ainp_d) &
      !$acc& private(ii,jf,nn,kk,i,m)
      do k = 1, lev
         do j = 1, my
            ii = 0
            !$acc loop seq
            do i = 1, nsizex
               jf = jlist2_2d(i, j)
               nn = nxjlen_all(i, j)
               kk = (k - 1)*nxp
               !$acc loop independent
               do m = 1, nn
                  ! wind at time step n
                  qq_3df(ii + m, j, k, 1) = ainp_d(kk + m, jf, 1)
                  qq_3df(ii + m, j, k, 2) = ainp_d(kk + m, jf, 2)
                  qq_3df(ii + m, j, k, 3) = ainp_d(kk + m, jf, 3)
               end do
               ii = ii + nn
            end do
         end do
      end do

      !$acc parallel loop collapse(3) deviceptr(qpwrk_d)&
      !$acc& private(ii,jf,nn,kk,i,m)
      do n = 1, ncld
         do k = 1, lev
            do j = 1, my
               ii = 0
               !$acc loop seq
               do i = 1, nsizex
                  jf = jlist2_2d(i, j)
                  nn = nxjlen_all(i, j)
                  kk = (k - 1)*nxp + (n - 1)*nxp*lev
                  !$acc loop independent
                  do m = 1, nn
                     qq_3df(ii + m, j, k, 3 + n) = qpwrk_d(kk + m, jf)
                  end do
                  ii = ii + nn
               end do
            end do
         end do
      end do
      !$acc parallel loop collapse(2) async(async_id)&
      !$acc& private(ii,jf,nn,kk,i,m)
      do k = 1, lev
         do j = 1, my
            ii = 0
            !$acc loop seq
            do i = 1, nsizex
               jf = jlist2_2d(i, j)
               nn = nxjlen_all(i, j)
               kk = (k - 1)*nxp
               !$acc loop independent
               do m = 1, nn
                  ! wind at time step n
                  um_3df(ii + m, j, k) = umwrk(kk + m, jf)
               end do
               ii = ii + nn
            end do
         end do
      end do

      call cyclic_cell_massadvx_gpu(qq_3df, um_3df, deltim, mass, forward, &
                                    nxdef, my, nx, my, lev, ndslhvar, -1)
      ! ========================================
      !$acc parallel loop collapse(2) &
      !$acc& deviceptr(aout_d) &
      !$acc& private(ii,jf,nn,kk,i,m)
      do k = 1, lev
         do j = 1, my
            ii = 0
            !$acc loop seq
            do i = 1, nsizex
               jf = jlist2_2d(i, j)
               nn = nxjlen_all(i, j)
               kk = (k - 1)*nxp
               !$acc loop independent
               do m = 1, nn
                  aout_d(kk + m, jf, 1) = qq_3df(ii + m, j, k, 1)
                  aout_d(kk + m, jf, 2) = qq_3df(ii + m, j, k, 2)
                  aout_d(kk + m, jf, 3) = qq_3df(ii + m, j, k, 3)
               end do
               ii = ii + nn
            end do
         end do
      end do

      !$acc parallel loop collapse(3) deviceptr(qpwrk_d)&
      !$acc& private(ii,jf,nn,kk,i,m)
      do n = 1, ncld
         do k = 1, lev
            do j = 1, my
               ii = 0
               !$acc loop seq
               do i = 1, nsizex
                  jf = jlist2_2d(i, j)
                  nn = nxjlen_all(i, j)
                  kk = (k - 1)*nxp + (n - 1)*nxp*lev
                  !$acc loop independent
                  do m = 1, nn
                     qpwrk_d(kk + m, jf) = qq_3df(ii + m, j, k, 3 + n)
                  end do
                  ii = ii + nn
               end do
            end do
         end do
      end do
      !$acc end data
   end if
   ! call mpe_barrier

   ! ----------------------------------------
   ! u v t at n
   !$acc host_data use_device(ut,vt,tt,qt)
   call advh_scatter4GPU_dev(ut, aout_d(1, 1, 1), 1, 0)
   call advh_scatter4GPU_dev(vt, aout_d(1, 1, 2), 1, 0)
   call advh_scatter4GPU_dev(tt, aout_d(1, 1, 3), 1, 0)
   ! rq at n
   call advh_scatter4GPU_dev(qt, qpwrk_d, ncld, 0)
   !$acc end host_data
   ! ----------------------------------------
   return

end subroutine massadvx_gpu_wrapper
! ============================================================
subroutine vcopy(out, inp, n)
   use const, only: RTYPE
   use param, only: nx, my_max, lev
   use index, only: nxp, nxdef, jlist1, jlistnum
   implicit none
   integer, intent(in):: n
   real(kind=RTYPE), intent(out):: out(nxp, lev*n, my_max)
   real(kind=RTYPE), intent(in):: inp(nxp, lev*n, my_max)

   integer i, j, k, nxj, jj

   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef(j)
      do k = 1, lev*n
         do i = 1, nxj
            out(i, k, jj) = inp(i, k, jj)
         end do
      end do
   end do
   return
end subroutine vcopy
! ============================================================
subroutine VarErr(Err, a, lda, b, ldb, lev, nvar)
   use const, only: RTYPE, numreduce, weight
   use param, only: nx, my_max, my, octahedral
   use index, only: nxdef, jlist1, jlistnum
   use rank, only: myrank
   use mpe
   implicit none
   real(kind=RTYPE), intent(out):: Err(4)
   integer, intent(in)::lda, ldb, lev, nvar
   real(kind=RTYPE), intent(in):: A(lda, lev, nvar, my_max), &
                                  B(ldb, lev, nvar, my_max)

   integer i, j, k, n, nxj, jj, pts
   real(kind=RTYPE) tmp, vamax, sum_local, pi
   if (octahedral) then
      pts = (20 + nx)*my*lev
   elseif (numreduce == -99) then
      pts = nx*my*lev
   end if
   pi = 4.*atan(1.)
   Err = 0.
   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef(j)
      sum_local = 0.
      do n = 1, nvar
         do k = 1, lev
            do i = 1, nxj
               tmp = A(i, k, n, jj) - B(i, k, n, jj)
               Err(1) = max(Err(1), abs(tmp))
               sum_local = sum_local + tmp**2
            end do
         end do
      end do
      Err(2) = Err(2) + sum_local*(2.*pi/nxj)/(lev*nvar)*weight(j)
      Err(3) = Err(3) + sum_local
   end do

   call mpe_global_max(Err(1), 1, RTYPE)
   call mpe_global_sum(Err(2), 2, RTYPE)
   Err(2) = sqrt(Err(2))
   Err(3) = sqrt(Err(3)/pts/nvar)
   Err(4) = vamax(a, lda, lev)

   return
end subroutine VarErr
! ============================================================
real(kind=8) function Vamax(a, lda, lev)
   use const, only: RTYPE
   use param, only: nx, my_max, my
   use index, only: nxdef, jlist1, jlistnum, levp
   use rank, only: myrank
   use mpe
   implicit none
   integer, intent(in)::lda, lev
   real(kind=RTYPE), intent(in):: A(lda, lev, my_max)

   integer i, j, k, nxj, jj

   vamax = 0.
   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef(j)
      do k = 1, lev
         do i = 1, nxj
            vamax = max(abs(A(i, k, jj)), vamax)
         end do
      end do
   end do

   call mpe_global_max(vamax, 1, RTYPE)

   return
end function Vamax
