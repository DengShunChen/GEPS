subroutine mpe_transpose_rs_sp_gpu(sbuf, rbuf, lev, n, m, nsize, comm)
   ! Present on device: sbuf, rbuf
   use const, only: RTYPE, MPI_RTYPE
   use mpi

   implicit none
   integer n, m, lev, nsize, i, j, k, ii, len_tr, ierr, comm

   real(kind=RTYPE) sbuf(lev, n, nsize, m), rbuf(lev, n, m*nsize) ! present on device
   real(kind=RTYPE) swork(lev, n, m, nsize)

   integer async_id

   async_id = 1
   len_tr = m*n

   !$acc enter data create(swork) async(async_id)

   !$acc parallel loop collapse(4) async(async_id)
   do j = 1, m
      do ii = 1, nsize
         do i = 1, n
            do k = 1, lev
               swork(k, i, j, ii) = sbuf(k, i, ii, j)
            end do
         end do
      end do
   end do

   !$acc wait(async_id)

   !$acc host_data use_device(swork, rbuf)
   call MPI_ALLTOALL(SWORK, LEN_TR*LEV, MPI_RTYPE, RBUF, LEN_TR*LEV, MPI_RTYPE, comm, IERR)
   !$acc end host_data

   !$acc exit data delete(swork) async(async_id)

   return
end

subroutine mpe_transpose_rs_gpu(sbuf, rbuf, lev, n, m, nsize, comm)
!
   use mpi

   implicit none
   integer n, m, lev, nsize, i, j, k, ii, len_tr, ierr, comm

   real sbuf(lev, n, nsize, m) ! Present on device
   real rbuf(lev, n, m*nsize) ! Present on device
#ifdef MPISP
#   ifdef UNIFY_ONLY
   real*8 swork(lev, n, m, nsize)
#   else
   real*4 swork(lev, n, m, nsize), rwork(lev, n, m*nsize)
#   endif
#else
   real*8 swork(lev, n, m, nsize)
#endif
   integer async_id

   async_id = 1
!
   len_tr = m*n
!
   !$acc enter data create(swork) async(async_id)
#ifdef MPISP
#ifndef UNIFY_ONLY
   !$acc enter data create(rwork) async(async_id)
#endif
#endif
   !$acc parallel loop collapse(4) async(async_id)
   do j = 1, m
   do ii = 1, nsize
   do i = 1, n
   do k = 1, lev
      swork(k, i, j, ii) = sbuf(k, i, ii, j)
   end do
   end do
   end do
   end do
   !$acc wait(async_id)
!
#ifdef MPISP
#   ifdef UNIFY_ONLY
   !$acc host_data use_device(swork, rbuf)
   call MPI_ALLTOALL(SWORK, LEN_TR*LEV, MPI_REAL8, &
                     RBUF, LEN_TR*LEV, MPI_REAL8, &
                     comm, IERR)
   !$acc end host_data
#   else
   !$acc host_data use_device(swork, rwork)
   call MPI_ALLTOALL(SWORK, LEN_TR*LEV, MPI_REAL4, &
                     rwork, LEN_TR*LEV, MPI_REAL4, &
                     comm, IERR)
   !$acc end host_data
   !$acc kernels async(async_id)
   rbuf = rwork
   !$acc end kernels
#   endif
#else
   !$acc host_data use_device(swork, rbuf)
   call MPI_ALLTOALL(SWORK, LEN_TR*LEV, MPI_REAL8, &
                     RBUF, LEN_TR*LEV, MPI_REAL8, &
                     comm, IERR)
   !$acc end host_data
#endif
   !$acc exit data delete(swork) async(async_id)
#ifdef MPISP
#ifndef UNIFY_ONLY
   !$acc exit data delete(rwork) async(async_id)
#endif
#endif
!
   return
end
