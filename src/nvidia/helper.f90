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

subroutine cusolver_check_helper(ierr, filename, line)
   use cusolverDN
   implicit none

   integer, intent(in) :: ierr
   character(len=*), intent(in) :: filename
   integer, intent(in) :: line

   if (ierr .ne. CUSOLVER_STATUS_SUCCESS) then
      print '(5g0)', "Failed, cuSOLVER error ", filename, ":", line
      print '(5g0)', "Error : ", ierr
   end if

end subroutine cusolver_check_helper

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

subroutine nccl_alltoall(sendbuf, sendcount, recvbuf, recvcount, comm, nsize, async_id)
   ! Present on device: sendbuf, recvbuf
   use const, only: RTYPE
   use openacc
   use cudafor
   use nccl

   implicit none

   real(kind=RTYPE), intent(in) :: sendbuf(sendcount, nsize)
   real(kind=RTYPE), intent(out) :: recvbuf(recvcount, nsize)
   integer, intent(in) :: sendcount, recvcount, nsize, async_id
   type(ncclComm), intent(in) :: comm
   type(ncclDataType) :: dtype

   integer :: i
   integer(kind=cuda_stream_kind) :: stream

   stream = acc_get_cuda_stream(async_id)

#ifdef SP
   dtype = ncclFloat32
#else
   dtype = ncclFloat64
#endif

   NCCLCHECK(ncclGroupStart())
   !$acc host_data use_device(sendbuf, recvbuf)
   do i = 1, nsize
      NCCLCHECK(ncclSend(sendbuf(1, i), sendcount, dtype, i - 1, comm, stream))
      NCCLCHECK(ncclRecv(recvbuf(1, i), recvcount, dtype, i - 1, comm, stream))
   end do
   !$acc end host_data
   NCCLCHECK(ncclGroupEnd())

end subroutine nccl_alltoall

subroutine nccl_alltoallv_stride(sendbuf, sendcount, senddisplace, recvbuf, recvcount, recvdisplace, comm, nsize, async_id)
   ! Currently the displacement is an integer representing the stride of data layout accorss processes
   use const, only: RTYPE
   use cudafor
   use nccl
   use openacc

   implicit none

   real(kind=RTYPE), intent(in) :: sendbuf(senddisplace, nsize), recvbuf(recvdisplace, nsize)
   integer, intent(in), dimension(nsize) :: sendcount, recvcount
   integer, intent(in) :: senddisplace, recvdisplace
   type(ncclComm), intent(in) :: comm
   type(ncclDataType) :: dtype
   integer, intent(in) :: nsize, async_id
   integer(kind=cuda_stream_kind) :: stream
   integer :: i

   stream = acc_get_cuda_stream(async_id)
#ifdef SP
   dtype = ncclFloat32
#else
   dtype = ncclFloat64
#endif
   NCCLCHECK(ncclGroupStart())
   !$acc host_data use_device(sendbuf, recvbuf)
   do i = 1, nsize
      NCCLCHECK(ncclSend(sendbuf(1, i), sendcount(i), dtype, i - 1, comm, stream))
      NCCLCHECK(ncclRecv(recvbuf(1, i), recvcount(i), dtype, i - 1, comm, stream))
   end do
   !$acc end host_data
   NCCLCHECK(ncclGroupEnd())

end subroutine nccl_alltoallv_stride

subroutine histogram(data, n, num_bins)
   use const, only: RTYPE

   implicit none

   real(kind=RTYPE), dimension(n), intent(in) :: data
   integer, intent(in) :: n
   integer, intent(in) :: num_bins
   real(kind=RTYPE) :: min_val, max_val, bin_width
   real(kind=RTYPE), dimension(num_bins + 1) :: bin_edges
   integer, dimension(num_bins) :: bins
   integer :: bin_index, i

   min_val = minval(data)
   max_val = maxval(data)
   bin_width = (max_val - min_val)/real(num_bins, RTYPE)

   bin_edges(1) = min_val
   do i = 2, num_bins + 1
      bin_edges(i) = min_val + (i - 1)*bin_width
   end do

   bins = 0
   do i = 1, n
      bin_index = min(num_bins, max(1, floor((data(i) - min_val)/bin_width) + 1))
      bins(bin_index) = bins(bin_index) + 1
   end do

   print *, bin_edges
   print *, bins

end subroutine histogram
