!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define NCCLCHECK(ierr) call nccl_check_helper(ierr, __FILE__, __LINE__)

program test_nccl
   use mpi
   use nccl
   use cudafor
   use openacc

   implicit none

   integer :: ierr, rank, nranks, num_device
   type(ncclUniqueId) :: nccl_id
   type(ncclResult) :: nccl_result
   type(ncclComm) :: nccl_comm
   real(4), allocatable :: sendbuf(:), recvbufmpi(:), recvbufnccl(:)
   integer, parameter :: buffer_size = 1000
   integer, parameter :: async_id = 1
   integer(kind=cuda_stream_kind) :: stream

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)

   num_device = acc_get_num_devices(acc_device_nvidia)
   if (num_device .ne. nranks) then
      print *, "Please set the number of MPI process to be equal to number of GPUs."
      print *, "# of MPI processes =", nranks
      print *, "# of GPU devices =", num_device
      call exit(1)
   end if
   call acc_set_device_num(rank, acc_device_nvidia)

   if (rank == 0) then
      NCCLCHECK(ncclGetUniqueId(nccl_id))
   end if

   call MPI_Bcast(nccl_id%internal, size(nccl_id%internal), MPI_BYTE, 0, MPI_COMM_WORLD, ierr)

   NCCLCHECK(ncclCommInitRank(nccl_comm, nranks, nccl_id, rank))

   allocate (sendbuf(buffer_size))
   allocate (recvbufmpi(buffer_size))
   allocate (recvbufnccl(buffer_size))

   call random_seed()
   call random_number(sendbuf)

   call MPI_Allreduce(sendbuf, recvbufmpi, buffer_size, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)

   !$acc enter data copyin(sendbuf, recvbufnccl) async(async_id)

   stream = acc_get_cuda_stream(async_id)
   !$acc host_data use_device(sendbuf, recvbufnccl)
   NCCLCHECK(ncclAllReduce(sendbuf, recvbufnccl, buffer_size, ncclFloat32, ncclSum, nccl_comm, stream))
   !$acc end host_data

   !$acc exit data copyout(sendbuf, recvbufnccl) async(async_id)
   !$acc wait(async_id)

   if (MAXVAL(ABS(recvbufnccl - recvbufmpi)) .ge. 1e-5) then
      print *, "test_nccl failed."
      call exit(1)
   end if

   NCCLCHECK(ncclCommDestroy(nccl_comm))

   deallocate (sendbuf)
   deallocate (recvbufmpi)
   deallocate (recvbufnccl)

   call MPI_Finalize(ierr)

end program
