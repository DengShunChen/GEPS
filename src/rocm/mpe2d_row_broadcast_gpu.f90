! First-class ROCm row broadcast: OpenMP use_device_addr + RCCL.
! 1-rank comm is a no-op (this lab has one full MI300X).

subroutine mpe2d_row_broadcast_gpu(buf, n, brank)
   use iso_c_binding
   use rank, only: nsize
   use index, only: nccl_row_comm
   use openacc
   use cudafor
   use nccl
   implicit none

   integer, intent(in) :: n, brank
   real(kind=8), intent(inout), target :: buf(n)

   integer :: async_id
   integer(kind=cuda_stream_kind) :: stream
   type(ncclResult) :: nerr

   if (nsize <= 1 .or. n <= 0) return

   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   !$omp target data use_device_addr(buf)
   nerr = ncclBroadcast(buf, buf, n, ncclFloat64, brank, nccl_row_comm, stream)
   !$omp end target data

   if (nerr /= ncclSuccess) then
      print *, 'geps_rocm_row_broadcast ncclBroadcast failed', nerr%val, 'n=', n
      error stop
   end if
end subroutine mpe2d_row_broadcast_gpu
