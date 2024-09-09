#define NCCLCHECK(ierr) call nccl_check_helper(ierr, __FILE__, __LINE__)

subroutine cuda_check_helper(ierr, filename, line)
   use cudafor

   implicit none

   integer, intent(in) :: ierr
   character(len=*), intent(in) :: filename
   integer, intent(in) :: line

   if (ierr .ne. cudaSuccess) then
      print '(5g0)', "Failed, CUDA error ", filename, ":", line
      print '(5g0)', cudaGetErrorString(ierr)
   end if

end subroutine cuda_check_helper

subroutine cufft_check_helper(ierr, filename, line)
   use cufft

   implicit none

   integer, intent(in) :: ierr
   character(len=*), intent(in) :: filename
   integer, intent(in) :: line

   if (ierr .ne. CUFFT_SUCCESS) then
      print '(5g0)', "Failed, cuFFT error ", filename, ":", line
      print '(5g0)', "Error : ", ierr
   end if

end subroutine cufft_check_helper

subroutine nccl_check_helper(ierr, filename, line)
   use nccl

   implicit none

   type(ncclResult), intent(in) :: ierr
   character(len=*), intent(in) :: filename
   integer, intent(in) :: line

   if (ierr .ne. ncclSuccess) then
      print '(5g0)', "Failed, NCCL error ", filename, ":", line
      print '(5g0)', ncclGetErrorString(ierr)
   end if

end subroutine nccl_check_helper

subroutine nccl_alltoall_fp64(sendbuf, sendcount, recvbuf, recvcount, comm, nsize, async_id)
   ! Present on device: sendbuf, recvbuf
   use openacc
   use cudafor
   use nccl

   implicit none

   real(8), intent(in) :: sendbuf(sendcount, nsize)
   real(8), intent(out) :: recvbuf(recvcount, nsize)
   integer, intent(in) :: sendcount, recvcount, nsize, async_id
   type(ncclComm), intent(in) :: comm

   integer :: i
   integer(kind=cuda_stream_kind) :: stream

   stream = acc_get_cuda_stream(async_id)

   NCCLCHECK(ncclGroupStart())
   !$acc host_data use_device(sendbuf, recvbuf)
   do i = 1, nsize
      NCCLCHECK(ncclSend(sendbuf(1, i), sendcount, ncclFloat64, i - 1, comm, stream))
      NCCLCHECK(ncclRecv(recvbuf(1, i), recvcount, ncclFloat64, i - 1, comm, stream))
   end do
   !$acc end host_data
   NCCLCHECK(ncclGroupEnd())

end subroutine nccl_alltoall_fp64

subroutine nccl_alltoallv_stride_fp64(sendbuf, sendcount, senddisplace, recvbuf, recvcount, recvdisplace, comm, nsize, async_id)
   ! Currently the displacement is an integer representing the stride of data layout accorss processes
   use cudafor
   use nccl
   use openacc

   implicit none

   real(8), intent(in) :: sendbuf(senddisplace, nsize), recvbuf(recvdisplace, nsize)
   integer, intent(in), dimension(nsize) :: sendcount, recvcount
   integer, intent(in) :: senddisplace, recvdisplace
   type(ncclComm), intent(in) :: comm
   integer, intent(in) :: nsize, async_id
   integer(kind=cuda_stream_kind) :: stream
   integer :: i

   stream = acc_get_cuda_stream(async_id)
   NCCLCHECK(ncclGroupStart())
   !$acc host_data use_device(sendbuf, recvbuf)
   do i = 1, nsize
      NCCLCHECK(ncclSend(sendbuf(1, i), sendcount(i), ncclFloat64, i - 1, comm, stream))
      NCCLCHECK(ncclRecv(recvbuf(1, i), recvcount(i), ncclFloat64, i - 1, comm, stream))
   end do
   !$acc end host_data
   NCCLCHECK(ncclGroupEnd())

end subroutine nccl_alltoallv_stride_fp64
