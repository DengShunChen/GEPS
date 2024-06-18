#define NCCLCHECK(ierr) call nccl_check_helper(ierr, __FILE__, __LINE__)

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
