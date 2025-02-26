#define CUSOLVERCHECK(ierr) call cuSOLVER_check_helper(ierr, __FILE__, __LINE__)
subroutine eigen_unit(no)
   use rank, only: myrank
   use const, only: RTYPE, nnmivm, rad, omega, eigval, cutfreq
   use param, only: lev, jtrun, jtmax
   use index, only: levp, mlist, mlistnum, Llist, Lstart, Lend, row_rank
   use cusolverDn
   use cudafor
   use cublas
   use openacc
   use mod_vcopy
   use mpe

   implicit none

   integer, parameter :: steps = 2
   integer, no
   real, dimension(no*no, 2*jtmax*nnmivm) :: mx, evec, epos, evec_c, mx_c
   real, dimension(no, 2*jtmax*nnmivm) :: eval, eval_c
   integer nnlist(2*jtmax*nnmivm)

   real wk(no*no)
   real, dimension(jtrun, jtrun, nnmivm):: a, b, c
   real, dimension(jtrun, jtmax, nnmivm):: h

   integer nw(jtrun, jtmax)
   integer L, k, mf, m, ns, na, nbig, j, ii, jj, i, nn, ind, isym
   integer brank   !root rank of row_broadcast

   real(kind=8) tm_1, tm_2, tm_use, mpi_wtime
   integer async_id, s, istat
   type(cublashandle) :: cublas_h
   type(cusolverDnHandle) cusolverDn_h

   istat = cublasCreate(cublas_h)
   istat = cusolverDNCreate(cusolverDn_h)
   async_id = 1
   evec_c = 0.
   mx = 0.
   eval = 0.

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

   ! << CPU >>
   tm_1 = mpi_wtime()
   do L = 1, nnmivm
      do m = 1, mlistnum
         do i = 1, 2
            j = i + (m - 1)*2 + (L - 1)*jtmax*2
            nn = nnlist(j)
            call eigen(mx(1, j), nn, &
                       eval_c(1, j), evec_c(1, j), epos(1, j), wk)
         end do
      end do
   end do
   tm_use = mpi_wtime() - tm_1
   if (myrank .eq. 0) write (*, '(A, f15.3, A)') &
      "Elapsed time of CPU: ", tm_use*1e+3, " ms"

   ! << GPU >>
   !$acc enter data copyin(mx) async(async_id)
   !$acc enter data create(eval) async(async_id)
   do s = 1, steps
      if (s .gt. 1) then
         !$acc update device(mx) async(async_id)
         !$acc wait(async_id)
      end if
      tm_1 = mpi_wtime()
      call eigen_mx_gpu(mx, eval, nnlist, no)
      !$acc wait(async_id)
      tm_use = mpi_wtime() - tm_1
      if (myrank .eq. 0) write (*, '(A, f15.3, A)') &
         "Elapsed time of GPU: ", tm_use*1e+3, " ms"
   end do
   !$acc exit data copyout(mx, eval) async(async_id)
   !$acc wait(async_id)

   !! align sign
   do L = 1, nnmivm
      do m = 1, mlistnum
         do i = 1, 2
            j = i + (m - 1)*2 + (L - 1)*jtmax*2
            nn = nnlist(j)
            if (nn .gt. 2) then
               do ii = 1, nn
                  do jj = 1, nn
                     ind = jj + (ii - 1)*nn
                     mx(ind, j) = eval(ii, j)*sign(mx(ind, j), evec_c(ind, j))
                     evec_c(ind, j) = eval_c(ii, j)*evec_c(ind, j)
                  end do
               end do
            else
               do ii = 1, nn
                  eval(ii, j) = 0.
                  do jj = 1, nn
                     ind = jj + (ii - 1)*nn
                     mx(ind, j) = 0.
                  end do
               end do
            end if
         end do
      end do
   end do
   ! ========================================

   call assert_allclose_r8(eval, size(eval), eval_c, size(eval_c), 1e-10, 1e-10, "eval")
   call assert_allclose_r8(mx, size(mx), evec_c, size(evec_c), 1e-10, 1e-10, "evec")

end subroutine eigen_unit

program test_eigen
   use rank, only: myrank
   use param, only: jtrun
   implicit none
   integer no

   call mpe_init
   call cons
   no = 2*((jtrun + 1)/2) + (jtrun/2) + 10
   call eigen_unit(no)
   call mpe_finalize

end program test_eigen
