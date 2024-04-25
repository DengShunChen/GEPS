subroutine siimpl_gpu(jtrun, jtmax, lev, dta, ptmean, dsigma, spalm, eps4 &
                      , eigval, evecin, evectr, arrhyd, arsddt, temold, divold, plold &
                      , temnow, divnow, plnow, temten, divten, plten, alpha)
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

   real(kind=RTYPE) divavg(lev, 2)

   integer m, mf, k, n, l, j
   real alpha
   real(kind=RTYPE) dta, dd, odd, dd2, ptmean, tem, s1, s2, d1, d2, dp

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
   call mpe2d_reshape_pl_multi(plten, plnow, plold, pltenL, plnowL, ploldL)

   call mpe2d_transpose_siimpl_multi(temold, temnow, temten, divold, divnow, divten, &
                                     wrk1, wrk2, wrk3, wrk4, wrk5, wrk6, &
                                     levp, jtrun, jtmax, lev, jtp, jtf, mlistnum, mlist, nsizex, row_comm)
#else
   !$acc enter data copyin(plten, plnow, plold) create(pltenL, plnowL, ploldL) async(async_id)
   call mpe2d_reshape_pl_gpu(plten, pltenL)
   call mpe2d_reshape_pl_gpu(plnow, plnowL)
   call mpe2d_reshape_pl_gpu(plold, ploldL)
   !$acc exit data copyout(pltenL, plnowL, ploldL) delete(plten, plnow, plold) async(async_id)
   !$acc wait(async_id)

   !$acc enter data copyin(temold, temnow, temten, divold, divnow, divten) create(wrk1, wrk2, wrk3, wrk4, wrk5, wrk6) async(async_id)
   call mpe2d_transpose_siimpl_gpu(temold, &
                                   wrk1, &
                                   levp, jtrun, jtmax, lev, jtp, jtf, mlistnum, mlist, nsizex, row_comm)
   call mpe2d_transpose_siimpl_gpu(temnow, &
                                   wrk2, &
                                   levp, jtrun, jtmax, lev, jtp, jtf, mlistnum, mlist, nsizex, row_comm)
   call mpe2d_transpose_siimpl_gpu(temten, &
                                   wrk3, &
                                   levp, jtrun, jtmax, lev, jtp, jtf, mlistnum, mlist, nsizex, row_comm)
   call mpe2d_transpose_siimpl_gpu(divold, &
                                   wrk4, &
                                   levp, jtrun, jtmax, lev, jtp, jtf, mlistnum, mlist, nsizex, row_comm)
   call mpe2d_transpose_siimpl_gpu(divnow, &
                                   wrk5, &
                                   levp, jtrun, jtmax, lev, jtp, jtf, mlistnum, mlist, nsizex, row_comm)
   call mpe2d_transpose_siimpl_gpu(divten, &
                                   wrk6, &
                                   levp, jtrun, jtmax, lev, jtp, jtf, mlistnum, mlist, nsizex, row_comm)
   !$acc exit data copyout(wrk1, wrk2, wrk3, wrk4, wrk5, wrk6) delete(temold, temnow, temten, divold, divnow, divten) async(async_id)
   !$acc wait(async_id)
#endif

   do m = 1, jtlen
      n = jtwvp(m)

      if (n .ne. 1) then
         do k = 1, lev
            divavg(k, 1) = wrk1(k, 1, m) + dd*wrk3(k, 1, m) - wrk2(k, 1, m)
            divavg(k, 2) = wrk1(k, 2, m) + dd*wrk3(k, 2, m) - wrk2(k, 2, m)
         end do
         do k = 1, lev
            s1 = spalm(k)*(ploldL(m, 1) + dd*pltenL(m, 1) - plnowL(m, 1))
            s2 = spalm(k)*(ploldL(m, 2) + dd*pltenL(m, 2) - plnowL(m, 2))
            do L = 1, lev
               s1 = s1 + arrhyd(k, L)*divavg(L, 1)
               s2 = s2 + arrhyd(k, L)*divavg(L, 2)
            end do
            wrk6(k, 1, m) = dd*(eps4L(m)*s1 + wrk6(k, 1, m)) &
                            - wrk5(k, 1, m) + wrk4(k, 1, m)
            wrk6(k, 2, m) = dd*(eps4L(m)*s2 + wrk6(k, 2, m)) &
                            - wrk5(k, 2, m) + wrk4(k, 2, m)
         end do
      end if
!
! transform time averaged divergence to eigenspace and compute
! semi-implicit values
!
      if (n .eq. 1) then

         do l = 1, lev
            wrk6(l, 1, m) = 0.0
            wrk6(l, 2, m) = 0.0
         end do
      else
         do L = 1, lev
            d1 = 0.
            d2 = 0.
            do k = 1, lev
               d1 = d1 + evecin(L, k)*wrk6(k, 1, m)
               d2 = d2 + evecin(L, k)*wrk6(k, 2, m)
            end do
            tem = dd2*eigval(L)
            eps4e = 1.0/(1.0 + tem*eps4L(m))
            divavg(L, 1) = d1*eps4e
            divavg(L, 2) = d2*eps4e
         end do

         do k = 1, lev
            s1 = 0.
            s2 = 0.
            do L = 1, lev
               s1 = s1 + evectr(k, L)*divavg(L, 1)
               s2 = s2 + evectr(k, L)*divavg(L, 2)
            end do
            wrk6(k, 1, m) = s1
            wrk6(k, 2, m) = s2
         end do
      end if
!
! add contributions of time averaged divergence to temperature
! tendency
!
      if (n .ne. 1) then
         do k = 1, lev
            s1 = wrk3(k, 1, m)
            s2 = wrk3(k, 2, m)
            do L = 1, lev
               s1 = s1 - arsddt(k, L)*wrk6(L, 1, m)
               s2 = s2 - arsddt(k, L)*wrk6(L, 2, m)
            end do
            wrk3(k, 1, m) = s1
            wrk3(k, 2, m) = s2
         end do
      end if
!
! add contribution of vertically integrated time-averaged divergence
! to surface pressure tendency.  convert time-averaged divergence to
! semi-implicit divergence tendency.
!
      if (n .ne. 1) then
         do k = 1, lev
         do j = 1, 2
            dp = dsigma(k, 1)*ptmean + dsigma(k, 2)
            pltenL(m, j) = pltenL(m, j) - dp*wrk6(k, j, m)
            wrk6(k, j, m) = (wrk5(k, j, m) + wrk6(k, j, m) - wrk4(k, j, m))*odd
         end do
         end do
      end if
   end do ! end of large me loop

   call mpe2d_reshape_pl_back(pltenL, plten)

#ifdef MULTIPLE
   call mpe2d_transpose_siimpl_back_multi(wrk3,wrk6,temten,divten,levp,jtrun,jtmax,lev,jtp,jtf,mlistnum,mlist,nsizex,row_comm)
#else
   call mpe2d_transpose_siimpl_back(wrk3, temten, levp, jtrun, jtmax, lev, jtp, jtf, mlistnum, mlist, nsizex, row_comm)
   call mpe2d_transpose_siimpl_back(wrk6, divten, levp, jtrun, jtmax, lev, jtp, jtf, mlistnum, mlist, nsizex, row_comm)
#endif

   return
end
