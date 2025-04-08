subroutine nnmi_unit(no)
   use rank, only: myrank
   use const, only: RTYPE, nnmivm, rad, omega, eigval, cutfreq
   use param, only: lev, jtrun, jtmax
   use index, only: levp, mlist, mlistnum, Llist, Lstart, Lend, row_rank
   use cudafor
   use openacc
   use mod_vcopy
   use mpe

   implicit none

   integer, parameter :: steps = 3
   integer, no
   real, dimension(no, 2, 2*jtmax*nnmivm) :: x, x_c
   real, dimension(no*no, 2*jtmax*nnmivm) :: mx, evec, epos
   real, dimension(no, 2*jtmax*nnmivm) :: eval
   real, dimension(jtrun, nnmivm) :: bal_tmp_c, bal_tmp
   real bal

   real wk(no*no), wc(no*2), wd(no*2), ew(no*no)
   real, dimension(jtrun, jtrun, nnmivm):: a, b, c
   real, dimension(jtrun, jtmax, nnmivm):: h
   integer nw(jtrun, jtmax)
   integer L, k, mf, m, ns, na, nbig, j, ii, jj
   integer brank   !root rank of row_broadcast
   logical nnmi_cg_created
   integer nnlist(2*jtmax*nnmivm)
   type(acc_graph_t) nnmi_cg
   real nnmi_buf(no, 2, 2*jtmax*nnmivm)

   real(kind=RTYPE) :: err(3)

   real(kind=8) tm_1, tm_2, tm_use, mpi_wtime
   integer async_id, s

   async_id = 1
   nnmi_cg_created = .false.

   call random_seed()
   call random_number(x)
   x = x*1e-5
   x_c = x
   do L = 1, nnmivm
      do m = 1, mlistnum
         mf = mlist(m)
         nbig = jtrun - mf + 1
         !  calculate the size of symmetric and antisymmetric matrix which
         !  include the gravity and rossby wave
         ns = 2*int((nbig + 1)/2) + int(nbig/2)
         na = 2*int(nbig/2) + int((nbig + 1)/2)
         j = 1 + (m - 1)*2 + (L - 1)*2*jtmax
         nnlist(j) = ns
         nnlist(j + 1) = na
      end do
   end do
   call inicons(rad, omega, eigval, lev, jtrun, jtmax, &
                nw, a, b, c, h, nnmivm)
   do L = 1, nnmivm
      if ((L .ge. Lstart) .and. (L .le. Lend)) then
         k = L - Lstart + 1
         do m = 1, mlistnum
            mf = mlist(m)
            nbig = jtrun - mf + 1
            !  calculate the size of symmetric and antisymmetric matrix which
            !  include the gravity and rossby wave
            ns = 2*int((nbig + 1)/2) + int(nbig/2)
            na = 2*int(nbig/2) + int((nbig + 1)/2)
            j = 1 + (m - 1)*2 + (L - 1)*2*jtmax
            ! symmetric matrix
            call coftrix(mf, L, mx(1:ns*ns, j), ns, &
                         a(1, 1, L), b(1, 1, L), c(1, 1, L), jtrun, lev, +1)
            ! antisymmetric matrix
            call coftrix(mf, L, mx(1:na*na, j + 1), na, &
                         a(1, 1, L), b(1, 1, L), c(1, 1, L), jtrun, lev, -1)
         end do
      end if
   end do
   do L = 1, nnmivm
      brank = 0
      if ((L .ge. Lstart) .and. (L .le. Lend)) then
         brank = row_rank
      end if
      j = 1 + (L - 1)*jtmax*2
      call mpe2d_row_broadcast(mx(1, j), no*no*2*jtmax, brank)
   end do
   do L = 1, nnmivm
      do m = 1, mlistnum
         mf = mlist(m)

         nbig = jtrun - mf + 1
         ns = 2*int((nbig + 1)/2) + int(nbig/2)
         na = 2*int(nbig/2) + int((nbig + 1)/2)

         j = 1 + (m - 1)*2 + (L - 1)*jtmax*2
         call eigen(mx(1:ns*ns, j), ns, &
                    eval(1:ns, j), &
                    evec(1:ns*ns, j), &
                    epos(1:ns*ns, j), wk)

         call eigen(mx(1:na*na, j + 1), na, &
                    eval(1:na, j + 1), &
                    evec(1:na*na, j + 1), &
                    epos(1:na*na, j + 1), wk)
      end do
   end do
   bal_tmp_c = 0.

   ! << CPU >>
   tm_1 = mpi_wtime()
   do L = 1, nnmivm
      do m = 1, mlistnum
         mf = mlist(m)
         nbig = jtrun - mf + 1
         ns = 2*int((nbig + 1)/2) + int(nbig/2)
         na = 2*int(nbig/2) + int((nbig + 1)/2)

         j = 1 + (m - 1)*2 + (L - 1)*jtmax*2

         call nnmi(no, x_c(1:ns, 1:2, j), &
                   epos(1:ns*ns, j), &
                   eval(1:ns, j), &
                   evec(1:ns*ns, j), &
                   ns, wc, wd, bal_tmp_c(mf, L), cutfreq)

         call nnmi(no, x_c(1:na, 1:2, j + 1), &
                   epos(1:na*na, j + 1), &
                   eval(1:na, j + 1), &
                   evec(1:na*na, j + 1), &
                   na, wc, wd, bal_tmp_c(mf, L), cutfreq)

      end do
   end do
   tm_use = mpi_wtime() - tm_1
   if (myrank .eq. 0) write (*, '(A, f15.3, A)') &
      "Elapsed time of CPU: ", tm_use*1e+3, " ms"
   do L = 1, nnmivm
      if (myrank .eq. 0) print *, 'vertical mode l=', L
      !
      ! set convergence variable
      !
      bal = 0.
      call mpe_unify(bal_tmp_c(1, L), 1, jtrun, 3, mpe_double)
      do mf = 1, jtrun
         bal = bal + bal_tmp_c(mf, L)
      end do
      !
      if (myrank .eq. 0) print *, 'bal=', bal
   end do

   ! << GPU >>
   bal_tmp = 0.
   !$acc enter data copyin(lev, jtrun, jtmax, levp, mlistnum, mlist, Llist, &
   !$acc&                  Lstart, Lend, cutfreq, eval, evec, nnlist) async(async_id)
   !$acc enter data create(nnmi_buf) async(async_id)
   !$acc enter data create(x, bal_tmp) async(async_id)
   do s = 1, steps
      !$acc update device(x, bal_tmp) async(async_id)
      !$acc wait(async_id)
      tm_1 = mpi_wtime()
      call nnmi_gpu(x, eval, evec, no, bal_tmp, cutfreq, nnlist, &
                    nnmi_cg_created, nnmi_cg, nnmi_buf)

      !$acc wait(async_id)
      tm_use = mpi_wtime() - tm_1
      if (myrank .eq. 0) write (*, '(A, f15.3, A)') &
         "Elapsed time of GPU: ", tm_use*1e+3, " ms"
   end do
   !$acc exit data delete(lev, jtrun, jtmax, levp, mlistnum, mlist, Llist, &
   !$acc&                 Lstart, Lend, cutfreq, eval, evec, nnlist) async(async_id)
   !$acc exit data delete(nnmi_buf)
   !$acc exit data copyout(x, bal_tmp) async(async_id)
   !$acc wait(async_id)

   do L = 1, nnmivm
      if (myrank .eq. 0) print *, 'vertical mode l=', L
      !
      ! set convergence variable
      !
      bal = 0.
      call mpe_unify(bal_tmp(1, L), 1, jtrun, 3, mpe_double)
      do mf = 1, jtrun
         bal = bal + bal_tmp(mf, L)
      end do
      !
      if (myrank .eq. 0) print *, 'bal=', bal
   end do

   call assert_allclose_r8(x, size(x), x_c, size(x_c), 1e-10, 1e-10, "x")
   call assert_allclose_r8(bal_tmp, size(bal_tmp), bal_tmp_c, size(bal_tmp_c), 1e-10, 1e-10, "bal_tmp")

end subroutine nnmi_unit

program test_nnmi
   use rank, only: myrank
   use param, only: jtrun
   implicit none
   integer no

   call mpe_init
   call cons
   no = 2*((jtrun + 1)/2) + (jtrun/2) + 10
   call nnmi_unit(no)
   call mpe_finalize

end program test_nnmi
