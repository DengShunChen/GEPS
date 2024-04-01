subroutine mpe_transpose_rs1_sp_gpu(sbuf, rbuf, n, m, lev, nsize, comm)

   use const, only: RTYPE, MPI_RTYPE
   use mpi

   implicit none
   integer n, m, lev, nsize, i, j, k, ii, len_tr, ierr, comm

   real(kind=RTYPE) sbuf(n, nsize, m, lev), rbuf(n, m * nsize, lev)
   real(kind=RTYPE) swork(lev, n, m, nsize), rwork(lev, n, m * nsize)

   len_tr = n * m

   !$acc data copyin(sbuf) copyout(rbuf) create(swork, rwork)
   !$acc parallel loop collapse(2) present(swork, sbuf)
   do j = 1, m
   do ii = 1, nsize
   do i = 1, n
   do k = 1, lev
      swork(k, i, j, ii) = sbuf(i, ii, j, k)
   end do
   end do
   end do
   end do
   !$acc host_data use_device(swork, rwork)
   call MPI_ALLTOALL(SWORK, LEN_TR * LEV, MPI_RTYPE, &
                     RWORK, LEN_TR * LEV, MPI_RTYPE, &
                     comm, IERR)
   !$acc end host_data
   !$acc parallel loop collapse(2) present(rbuf, rwork)
   do k = 1, lev
   do j = 1, m * nsize
   do i = 1, n
      rbuf(i, j, k) = rwork(k, i, j)
   end do
   end do
   end do
   !$acc end data

   return
end subroutine mpe_transpose_rs1_sp_gpu
