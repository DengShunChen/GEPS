subroutine mpe_transpose_rs_sp_gpu(sbuf, rbuf, lev, n, m, nsize, comm)
   ! Present on device: sbuf, rbuf
   use const, only: RTYPE
   use nccl

   implicit none
   integer n, m, lev, nsize, i, j, k, ii, len_tr
   type(ncclComm) :: comm

   real(kind=RTYPE) sbuf(lev, n, nsize, m), rbuf(lev, n, m*nsize)
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

   call nccl_alltoall(swork, len_tr*lev, rbuf, len_tr*lev, comm, nsize, async_id)

   !$acc exit data delete(swork) async(async_id)

   return
end

subroutine mpe_transpose_rs_gpu(sbuf, rbuf, lev, n, m, nsize, comm)
   ! Present on device: sbuf, rbuf
   use nccl

   implicit none
   integer n, m, lev, nsize, i, j, k, ii, len_tr, ierr
   type(ncclComm) :: comm

   real sbuf(lev, n, nsize, m)
   real rbuf(lev, n, m*nsize)
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

#ifdef MPISP
#   ifdef UNIFY_ONLY
   call nccl_alltoall(swork, len_tr*lev, rbuf, len_tr*lev, comm, nsize, async_id)
#   else
   call nccl_alltoall(swork, len_tr*lev, rwork, len_tr*lev, comm, nsize, async_id)
   !$acc kernels async(async_id)
   rbuf = rwork
   !$acc end kernels
#   endif
#else
   call nccl_alltoall(swork, len_tr*lev, rbuf, len_tr*lev, comm, nsize, async_id)
#endif
   !$acc exit data delete(swork) async(async_id)
#ifdef MPISP
#ifndef UNIFY_ONLY
   !$acc exit data delete(rwork) async(async_id)
#endif
#endif

   return
end
