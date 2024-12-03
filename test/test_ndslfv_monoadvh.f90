program test_ndslfv_monoadvh
   implicit none
   integer xy

   call mpe_init
   call cons

   ! << advh2_xy >>
   xy = 1
   call ndslfv_monoadvh_unit(xy)
   ! << advh2_yx >>
   xy = -1
   call ndslfv_monoadvh_unit(xy)

   call mpe_finalize
end program test_ndslfv_monoadvh
! ============================================================
subroutine ndslfv_monoadvh_unit(xy)
   use const, only: RTYPE, dt, cosl
   use rank, only: myrank
   use param, only: nx, my, my_max, lev, ncld
   use index, only: nxp, levp, levf, myf, jlistnum, jlen, nsizex, row_comm, &
                    nxdef, jlist1, jlist2_2d, nxjlen_all, nxdef
   use grid, only: ut, vt, tt, qm, qt, ut_sl, vt_sl, &
                   gglati, fa1, fa2, fa3, fa4
   use mod_ndslfv_monoadv_gpu, only: ndslfv_monoadvh_gpu
   use mpe
   implicit none
   integer, intent(in):: xy
   integer, parameter:: steps = 5
   integer, parameter:: nvar = 4

   real(kind=RTYPE):: um(nxp, lev, my_max), &
                      vm(nxp, lev, my_max), &
                      pten(nxp, lev, my_max)

   real(kind=RTYPE):: uum_sl(nx, levp, my_max), &
                      vvm_sl(nx, levp, my_max), &
                      ttm_sl(nx, levp, my_max), &
                      qm_sl(nx, levp*ncld, my_max), &
                      pten_sl(nx, levp, my_max)

   real(kind=RTYPE):: ut_cpu(nxp, lev, my_max), &
                      vt_cpu(nxp, lev, my_max), &
                      tt_cpu(nxp, lev, my_max), &
                      qt_cpu(nxp, lev*ncld, my_max)

   real(kind=RTYPE):: ut_gpu(nxp, lev, my_max), &
                      vt_gpu(nxp, lev, my_max), &
                      tt_gpu(nxp, lev, my_max), &
                      qt_gpu(nxp, lev*ncld, my_max)

   real(kind=RTYPE) dtah, irma, irm2a
   logical forward
   integer i, k, j, jj, nxj

   real(kind=RTYPE):: err_arr(4, nvar), vamax, ummax, vmmax
   character(len=6):: name(nvar)

   name = (/'t', 'u', 'v', 'q'/)
   dtah = 0.5*dt
   forward = .false.

   if (myrank .eq. 0) then
      print *, "========================================"
      if (xy .gt. 0.5) then
         print *, "    start test ndslfv_monoadvh_xy"
      elseif (xy .le. -0.5) then
         print *, "    start test ndslfv_monoadvh_yx"
      end if
      print *, "========================================"
      write (*, '(1X, A5, 1X, f15.2)') "dt=", dt
      write (*, '(1X, A5, 1X, f15.2)') "dtah=", dtah
   end if

   call random_seed()
   call random_number(um)
   call random_number(vm)
   call random_number(ut)
   call random_number(vt)
   call random_number(qm)

   call gen_cosinebell(tt, 0., 1.)

   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef(j)
      irma = cosl(j)
      irm2a = cosl(j)/irma
      do k = 1, lev
         do i = 1, nxj
            um(i, k, jj) = um(i, k, jj)*irma/dt*1e-3
            vm(i, k, jj) = vm(i, k, jj)*irm2a/dt*1e-5
         end do
      end do
   end do
   ummax = vamax(um, nxp, lev)
   vmmax = vamax(vm, nxp, lev)
   if (myrank .eq. 0) then
      write (*, '(A, 1pe15.7)'), "umamax=", ummax
      write (*, '(A, 1pe15.7)'), "vmamax=", vmmax
   end if

   !$acc enter data copyin(nx, my, lev, ncld, nxp, nsizex, &
   !$acc& gglati, fa1, fa2, fa3, fa4, jlist2_2d, nxjlen_all, cosl, nxdef)

   !$acc data create(tt_gpu, ut_gpu, vt_gpu, qt_gpu, &
   !$acc& um, vm, qm)
   do i = 1, steps
      if (myrank .eq. 0) write (*, '("<< ", i3, " >>")') i

      ! << CPU >>
      call mpe2d_transpose_ndsl_p2f(um, ut_sl, &
                                    nxp, nx, levf, levp, 1, myf, &
                                    my_max, jlistnum, jlen, nsizex, row_comm)
      call mpe2d_transpose_ndsl_p2f(vm, vt_sl, &
                                    nxp, nx, levf, levp, 1, myf, &
                                    my_max, jlistnum, jlen, nsizex, row_comm)
      call mpe2d_transpose_ndsl_p2f(ut, uum_sl, &
                                    nxp, nx, levf, levp, 1, myf, &
                                    my_max, jlistnum, jlen, nsizex, row_comm)
      call mpe2d_transpose_ndsl_p2f(vt, vvm_sl, &
                                    nxp, nx, levf, levp, 1, myf, &
                                    my_max, jlistnum, jlen, nsizex, row_comm)
      call mpe2d_transpose_ndsl_p2f(tt, ttm_sl, &
                                    nxp, nx, levf, levp, 1, myf, &
                                    my_max, jlistnum, jlen, nsizex, row_comm)
      call mpe2d_transpose_ndsl_p2f(qm, qm_sl, &
                                    nxp, nx, levf, levp, ncld, myf, &
                                    my_max, jlistnum, jlen, nsizex, row_comm)

      call ndslfv_monoadvh(ttm_sl, qm_sl, pten_sl, uum_sl, vvm_sl, &
                           nxdef, dtah, xy, levp)

      call mpe2d_transpose_ndsl_f2p(ttm_sl, tt_cpu, &
                                    nxp, nx, levf, levp, 1, myf, &
                                    my_max, jlistnum, jlen, nsizex, row_comm)
      call mpe2d_transpose_ndsl_f2p(uum_sl, ut_cpu, &
                                    nxp, nx, levf, levp, 1, myf, &
                                    my_max, jlistnum, jlen, nsizex, row_comm)
      call mpe2d_transpose_ndsl_f2p(vvm_sl, vt_cpu, &
                                    nxp, nx, levf, levp, 1, myf, &
                                    my_max, jlistnum, jlen, nsizex, row_comm)
      call mpe2d_transpose_ndsl_f2p(qm_sl, qt_cpu, &
                                    nxp, nx, levf, levp, ncld, myf, &
                                    my_max, jlistnum, jlen, nsizex, row_comm)

      ! << GPU >>
      call vcopy(tt_gpu, tt, 1)
      call vcopy(ut_gpu, ut, 1)
      call vcopy(vt_gpu, vt, 1)
      !$acc update device(tt_gpu, ut_gpu, vt_gpu)
      !$acc update device(um, vm)
      !$acc update device(qm)
      call ndslfv_monoadvh_gpu(tt_gpu, pten, ut_gpu, vt_gpu, qt_gpu, &
                               qm, um, vm, dtah, xy, forward)
      !$acc update self(tt_gpu, ut_gpu, vt_gpu, qt_gpu)
   end do
   !$acc end data

   !$acc exit data delete(nx, my, lev, ncld, nxp, nsizex, &
   !$acc& gglati, fa1, fa2, fa3, fa4, jlist2_2d, nxjlen_all, cosl, nxdef)

   call Varerr(err_arr(1, 1), tt_gpu, nxp, tt_cpu, nxp, lev, 1)
   call Varerr(err_arr(1, 2), ut_gpu, nxp, ut_cpu, nxp, lev, 1)
   call Varerr(err_arr(1, 3), vt_gpu, nxp, vt_cpu, nxp, lev, 1)
   call Varerr(err_arr(1, 4), qt_gpu, nxp, qt_cpu, nxp, lev, ncld)

   if (myrank .eq. 0) then
      do k = 1, nvar
         write (*, '(1X, A6, 1pe23.15, 3(1pe15.7))') &
            name(k), err_arr(1:4, k)
      end do
   end if

   if (all(err_arr(1, 1:nvar) < 1e-10)) then
      if (myrank .eq. 0) write (*, '(A,i3,A)') &
         "test_ndslfv_monoadvh (xy=", xy, ") passed."
   else
      if (myrank .eq. 0) write (*, '(A,i3,A)') &
         "test_ndslfv_monoadvh (xy=", xy, ") failed."
      call exit(1)
   end if
end subroutine ndslfv_monoadvh_unit
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
subroutine VarErr(Err, a, lda, b, ldb, lev, nvar)
   use const, only: RTYPE, numreduce, weight
   use param, only: nx, my_max, my, octahedral
   use index, only: nxdef, jlist1, jlistnum
   use rank, only: myrank
   use mpe
   implicit none
   real(kind=RTYPE), intent(out):: Err(4)
   integer, intent(in)::lda, ldb, lev, nvar
   real(kind=RTYPE), intent(in):: A(lda, lev*nvar, my_max), &
                                  B(ldb, lev*nvar, my_max)

   integer i, j, k, nxj, jj, pts
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
      do k = 1, lev*nvar
         do i = 1, nxj
            tmp = A(i, k, jj) - B(i, k, jj)
            Err(1) = max(Err(1), abs(tmp))
            sum_local = sum_local + tmp**2
         end do
      end do
      Err(2) = Err(2) + sum_local*(2.*pi/nxj)/(lev*nvar)*weight(j)
      Err(3) = Err(3) + sum_local
   end do

   call mpe_global_max(Err(1), 1, RTYPE)
   call mpe_global_sum(Err(2), 2, RTYPE)
   Err(2) = sqrt(Err(2))
   Err(3) = sqrt(Err(3)/pts/nvar)
   err(4) = vamax(a, lda, lev)

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
