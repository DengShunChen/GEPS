!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define NCCLCHECK(ierr) call nccl_check_helper(ierr, __FILE__, __LINE__)

subroutine device_init(rank, size)
   use openacc
   use cublas
   use cusparse
   use cusolverDn
   implicit none
   integer :: rank, size, istat
   integer :: num_device, device_id
   type(cublashandle) :: handle
   type(cusparseHandle) :: sparsehandle
   type(cusolverDnHandle) :: solverDnhandle
   num_device = acc_get_num_devices(acc_device_nvidia)
   device_id = mod(rank, num_device)
   call acc_set_device_num(device_id, acc_device_nvidia)
   istat = cublasCreate(handle)
   istat = cusparseCreate(sparsehandle)
   istat = cusolverDnCreate(solverDnhandle)
end subroutine device_init

subroutine nccl_init()
   use mpi
   use nccl
   use rank, only: nsize, myrank, MPI_COMM_gfs, nccl_id, nccl_comm_gfs
   use index, only: nccl_row_comm, nccl_col_comm, mrow, ncol, row_rank, col_rank

   implicit none

   integer :: ierr

   if (myrank == 0) then
      NCCLCHECK(ncclGetUniqueId(nccl_id))
   end if
   call MPI_Bcast(nccl_id%internal, size(nccl_id%internal), MPI_BYTE, 0, MPI_COMM_gfs, ierr)

   NCCLCHECK(ncclCommInitRank(nccl_comm_gfs, nsize, nccl_id, myrank))
   NCCLCHECK(ncclCommSplit(nccl_comm_gfs, mrow, row_rank, nccl_row_comm, c_null_ptr))
   NCCLCHECK(ncclCommSplit(nccl_comm_gfs, ncol, col_rank, nccl_col_comm, c_null_ptr))

end subroutine nccl_init

subroutine nccl_destroy()
   use nccl
   use rank, only: nccl_comm_gfs
   use index, only: nccl_row_comm, nccl_col_comm

   implicit none

   NCCLCHECK(ncclCommDestroy(nccl_row_comm))
   NCCLCHECK(ncclCommDestroy(nccl_col_comm))
   NCCLCHECK(ncclCommDestroy(nccl_comm_gfs))

end subroutine nccl_destroy
