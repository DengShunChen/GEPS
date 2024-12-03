subroutine mpe_transpose_sr_sp_gpu(sbuf, rbuf, lev, n, m, nsize, comm)
   ! Present on device: sbuf, rbuf
   use const, only: RTYPE
   use nccl

   implicit none
   integer n, m, lev, nsize, i, j, k, ii, len_tr
   type(ncclComm) :: comm

   real(kind=RTYPE) sbuf(lev, n, m*nsize), rbuf(lev, n, nsize, m)
   real(kind=RTYPE) rwork(lev, n, m, nsize)
   integer async_id

   async_id = 1
   len_tr = m*n

   !$acc enter data create(rwork) async(async_id)

   call nccl_alltoall(sbuf, len_tr*lev, rwork, len_tr*lev, comm, nsize, async_id)

   !$acc parallel loop collapse(4) async(async_id)
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

end
