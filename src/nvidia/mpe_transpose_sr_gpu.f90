      subroutine mpe_transpose_sr_gpu(sbuf, rbuf, lev, n, m, nsize, comm)
!
         use mpi

         implicit none
         integer n, m, lev, nsize, i, j, k, ii, len_tr, ierr, comm

         real sbuf(lev, n, m*nsize)
         real rbuf(lev, n, nsize, m)
#ifdef MPISP
#   ifdef UNIFY_ONLY
         real*8 rwork(lev, n, m, nsize)
#   else
         real*4 rwork(lev, n, m, nsize), sbuf1(lev, n, m*nsize)
#   endif
#else
         real*8 rwork(lev, n, m, nsize)
#endif
!
!
         len_tr = m*n

#ifdef MPISP
#   ifdef UNIFY_ONLY
         call MPI_ALLTOALL(SBUF, LEN_TR*LEV, MPI_REAL8, &
                           RWORK, LEN_TR*LEV, MPI_REAL8, &
                           comm, IERR)
#   else
         sbuf1 = sbuf
         call MPI_ALLTOALL(sbuf1, LEN_TR*LEV, MPI_REAL4, &
                           RWORK, LEN_TR*LEV, MPI_REAL4, &
                           comm, IERR)
#   endif
#else
         call MPI_ALLTOALL(SBUF, LEN_TR*LEV, MPI_REAL8, &
                           RWORK, LEN_TR*LEV, MPI_REAL8, &
                           comm, IERR)
#endif
!
         do j = 1, m
         do ii = 1, nsize
         do i = 1, n
         do k = 1, lev
            rbuf(k, i, ii, j) = rwork(k, i, j, ii)
         end do
         end do
         end do
         end do
!
         return
      end
!------------------------------------------------------------
!CWB2021 for single precision test
      subroutine mpe_transpose_sr_sp_gpu(sbuf, rbuf, lev, n, m, nsize, comm)
         ! Present on device: sbuf, rbuf
         use const, only: RTYPE, MPI_RTYPE
         use mpi

         implicit none
         integer n, m, lev, nsize, i, j, k, ii, len_tr, ierr, comm

         real(kind=RTYPE) sbuf(lev, n, m*nsize), rbuf(lev, n, nsize, m)
         real(kind=RTYPE) rwork(lev, n, m, nsize)
         integer async_id

         async_id = 1
         len_tr = m*n

         !$acc enter data create(rwork) async(async_id)
         !$acc wait(async_id)

         !$acc host_data use_device(sbuf, rwork)
         call MPI_ALLTOALL(SBUF, LEN_TR*LEV, MPI_RTYPE, &
                           RWORK, LEN_TR*LEV, MPI_RTYPE, &
                           comm, IERR)

         !$acc end host_data

         !$acc parallel loop collapse(4)
         do j = 1, m
         do ii = 1, nsize
         do i = 1, n
         do k = 1, lev
            rbuf(k, i, ii, j) = rwork(k, i, j, ii)
         end do
         end do
         end do
         end do

         !$acc exit data delete(rwork) async(async_id)

         return
      end

      subroutine mpe_transpose_sr_sp_async(sbuf, rbuf, lev, n, m, nsize, comm, async_id)
!
         use const, only: RTYPE, MPI_RTYPE
         use mpi

         implicit none
         integer n, m, lev, nsize, i, j, k, ii, len_tr, ierr, comm
         integer async_id

         real(kind=RTYPE) sbuf(lev, n, m*nsize), rbuf(lev, n, nsize, m)
         real(kind=RTYPE) rwork(lev, n, m, nsize)

         len_tr = m*n

         !$acc enter data create(rwork) async(async_id)
         !$acc wait(async_id)
         !$acc host_data use_device(sbuf, rwork)
         call MPI_ALLTOALL(SBUF, LEN_TR*LEV, MPI_RTYPE, &
                           RWORK, LEN_TR*LEV, MPI_RTYPE, &
                           comm, IERR)

         !$acc end host_data
         !$acc parallel loop collapse(2) async(async_id)
         do j = 1, m
         do ii = 1, nsize
         do i = 1, n
         do k = 1, lev
            rbuf(k, i, ii, j) = rwork(k, i, j, ii)
         end do
         end do
         end do
         end do
         !$acc exit data delete(rwork) async(async_id)
!
         return
      end
