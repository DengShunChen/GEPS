subroutine siimpl_gpu(jtrun, jtmax, lev, dta, ptmean, dsigma, spalm, eps4 &
                      , eigval, evecin, evectr, arrhyd, arsddt, temold, divold, plold &
                      , temnow, divnow, plnow, temten, divten, plten, alpha)
! Present on device: plten, plnow, plold, temold, temnow, temten, divold, divnow, divten,
! Present on device: jtwvp, spalm, arrhyd, eps4L, evecin, eigval, evectr, arsddt, dsigma, mlist
!
!
!  computes corrections to explicit tendencies to convert model to a
!  semi-implicit model
!
!  *** input ***
!
!  dta: time step in seconds
!  ptmean: mean terrain pressure chosen for best stability properties
!  dsigma: thickness of sigma layers
!  spalm: mean energy conversion terms for each level
!  eps4: spherical harmonic laplacian operator
!  eigval: gravity mode phase speed eigenvalues
!  evecin: inverse of gravity mode eigenvector matrix
!  evectr: gravity mode eigenvector matrix
!  arrhyd: linearized hydrostatic matrix
!  arsddt: linerized vertical temperature advection matrix
!  temold: (t-dt) spectral temperature
!  divold: (t-dt) spectral divergence
!  plold: (t-dt) spectral terrain pressure
!  temnow: current time temperature
!  divnow: current time divergence
!  plnow: current time terrain pressure
!  temten: explicit temperature tendency
!  divten: explicit divergence tendency
!  plten: explicit terrain pressure tendency
!
! *** output ***
!
!  temten: semi-implicit temperature tendency
!  divten: semi-implicit divergence tendency
!  plten: semi-implicit terrain pressure tendency
!
! **************************************************
!

!CWB2017 2dMPI version

   use index
   use paramt
   use const, only: eps4L, RTYPE
   use spec, only: plnowL, ploldL, pltenL, jtwvp
   use openacc
   use cudafor

   implicit none
   integer jtrun, jtmax, lev

   real(kind=RTYPE) dsigma(lev, 2), eps4(jtrun, jtmax), eigval(lev), evecin(lev, lev) &
      , evectr(lev, lev), arrhyd(lev, lev), arsddt(lev, lev), spalm(lev)

   real(kind=RTYPE) temold(levp, 2, jtrun, jtmax), divold(levp, 2, jtrun, jtmax) &
      , temnow(levp, 2, jtrun, jtmax), divnow(levp, 2, jtrun, jtmax) &
      , temten(levp, 2, jtrun, jtmax), divten(levp, 2, jtrun, jtmax)
   real(kind=RTYPE) plold(jtrun, jtmax, 2), plnow(jtrun, jtmax, 2), plten(jtrun, jtmax, 2)

   real(kind=RTYPE) divavg(lev, 2, jtlen)

   integer m, mf, k, n, l, j, i
   real alpha
   real(kind=RTYPE) dta, dd, odd, dd2, ptmean, tem, s1, s2, d1, d2, dp, s, d

   real(kind=RTYPE) wrk1(lev, 2, jtp), wrk2(lev, 2, jtp), wrk3(lev, 2, jtp), &
      wrk4(lev, 2, jtp), wrk5(lev, 2, jtp), wrk6(lev, 2, jtp)
   integer async_id, istat
   integer(kind=cuda_stream_kind) :: stream
   real(kind=RTYPE) eps4e, phiave1, phiave2

   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   dd = alpha*dta
   odd = 1.0/dd
   dd2 = dd*dd

#ifdef MULTIPLE
   print *, "Symbol MULTIPLE is not supported."
   call exit(1)
#endif

   !$acc enter data async(async_id) &
   !$acc& create(pltenL, plnowL, ploldL, wrk1, wrk2, wrk3, wrk4, wrk5, wrk6, divavg)

   call mpe2d_reshape_pl_gpu(plten, pltenL)
   call mpe2d_reshape_pl_gpu(plnow, plnowL)
   call mpe2d_reshape_pl_gpu(plold, ploldL)

   call mpe2d_transpose_siimpl_gpu(temold, &
                                   wrk1, &
                                   levp, jtrun, jtmax, lev, jtp, jtf, mlistnum, mlist, nsizex, nccl_row_comm)
   call mpe2d_transpose_siimpl_gpu(temnow, &
                                   wrk2, &
                                   levp, jtrun, jtmax, lev, jtp, jtf, mlistnum, mlist, nsizex, nccl_row_comm)
   call mpe2d_transpose_siimpl_gpu(temten, &
                                   wrk3, &
                                   levp, jtrun, jtmax, lev, jtp, jtf, mlistnum, mlist, nsizex, nccl_row_comm)
   call mpe2d_transpose_siimpl_gpu(divold, &
                                   wrk4, &
                                   levp, jtrun, jtmax, lev, jtp, jtf, mlistnum, mlist, nsizex, nccl_row_comm)
   call mpe2d_transpose_siimpl_gpu(divnow, &
                                   wrk5, &
                                   levp, jtrun, jtmax, lev, jtp, jtf, mlistnum, mlist, nsizex, nccl_row_comm)
   call mpe2d_transpose_siimpl_gpu(divten, &
                                   wrk6, &
                                   levp, jtrun, jtmax, lev, jtp, jtf, mlistnum, mlist, nsizex, nccl_row_comm)

   !$acc parallel loop collapse(3) private(n) async(async_id)
   do m = 1, jtlen
      do i = 1, 2
         do k = 1, lev
            n = jtwvp(m)
            if (n .ne. 1) then
               divavg(k, i, m) = wrk1(k, i, m) + dd*wrk3(k, i, m) - wrk2(k, i, m)
            end if
         end do
      end do
   end do

   !$acc parallel loop collapse(3) private(n, s) async(async_id)
   do m = 1, jtlen
      do i = 1, 2
         do k = 1, lev
            n = jtwvp(m)
            if (n .ne. 1) then
               s = spalm(k)*(ploldL(m, i) + dd*pltenL(m, i) - plnowL(m, i))
               !$acc loop seq
               do L = 1, lev
                  s = s + arrhyd(k, L)*divavg(L, i, m)
               end do
               wrk6(k, i, m) = dd*(eps4L(m)*s + wrk6(k, i, m)) &
                               - wrk5(k, i, m) + wrk4(k, i, m)
            end if
         end do
      end do
   end do

!
! transform time averaged divergence to eigenspace and compute
! semi-implicit values
!
   !$acc parallel loop collapse(3) private(n, d, tem, eps4e) async(async_id)
   do m = 1, jtlen
      do i = 1, 2
         do k = 1, lev
            n = jtwvp(m)
            if (n .ne. 1) then
               d = 0.0
               !$acc loop seq
               do L = 1, lev
                  d = d + evecin(k, L)*wrk6(L, i, m)
               end do
               tem = dd2*eigval(k)
               eps4e = 1.0/(1.0 + tem*eps4L(m))
               divavg(k, i, m) = d*eps4e
            end if
         end do
      end do
   end do

   !$acc parallel loop collapse(3) private(n, s) async(async_id)
   do m = 1, jtlen
      do i = 1, 2
         do k = 1, lev
            n = jtwvp(m)
            s = 0.0
            if (n .ne. 1) then
               s = 0.0
               !$acc loop seq
               do L = 1, lev
                  s = s + evectr(k, L)*divavg(L, i, m)
               end do
            end if
            wrk6(k, i, m) = s
         end do
      end do
   end do
!
! add contributions of time averaged divergence to temperature
! tendency
!
   !$acc parallel loop collapse(3) private(n, s) async(async_id)
   do m = 1, jtlen
      do i = 1, 2
         do k = 1, lev
            n = jtwvp(m)
            if (n .ne. 1) then
               s = wrk3(k, i, m)
               !$acc loop seq
               do L = 1, lev
                  s = s - arsddt(k, L)*wrk6(L, i, m)
               end do
               wrk3(k, i, m) = s
            end if
         end do
      end do
   end do
!
! add contribution of vertically integrated time-averaged divergence
! to surface pressure tendency.  convert time-averaged divergence to
! semi-implicit divergence tendency.
!
   !$acc parallel loop collapse(2) private(n, dp) async(async_id)
   do m = 1, jtlen
      do i = 1, 2
         !$acc loop seq
         do k = 1, lev
            n = jtwvp(m)
            if (n .ne. 1) then
               dp = dsigma(k, 1)*ptmean + dsigma(k, 2)
               pltenL(m, i) = pltenL(m, i) - dp*wrk6(k, i, m)
            end if
         end do
      end do
   end do

   !$acc parallel loop collapse(3) private(n, dp) async(async_id)
   do m = 1, jtlen
      do i = 1, 2
         do k = 1, lev
            n = jtwvp(m)
            if (n .ne. 1) then
               wrk6(k, i, m) = (wrk5(k, i, m) + wrk6(k, i, m) - wrk4(k, i, m))*odd
            end if
         end do
      end do
   end do

   call mpe2d_reshape_pl_back_gpu(pltenL, plten)

   call mpe2d_transpose_siimpl_back_gpu(wrk3, temten, levp, jtrun, jtmax, lev, jtp, jtf, mlistnum, mlist, nsizex, nccl_row_comm)
   call mpe2d_transpose_siimpl_back_gpu(wrk6, divten, levp, jtrun, jtmax, lev, jtp, jtf, mlistnum, mlist, nsizex, nccl_row_comm)

   !$acc exit data async(async_id) &
   !$acc& delete(pltenL, plnowL, ploldL, wrk1, wrk2, wrk3, wrk4, wrk5, wrk6, divavg)
   return
end
