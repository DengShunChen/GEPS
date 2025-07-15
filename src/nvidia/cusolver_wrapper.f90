#define CUSOLVERCHECK(ierr) call cuSOLVER_check_helper(ierr, __FILE__, __LINE__)
#define CUDACHECK(ierr) call cuda_check_helper(ierr, __FILE__, __LINE__)

subroutine DnDsyevd_async(jobz, uplo, n, A, lda, W, &
                          devinfo, async_id)
   use openacc
   use cudafor
   use cusolverDN
   implicit none
   integer :: jobz, uplo
   integer :: n, lda
   real(kind=8), dimension(lda, *) :: A
   real(kind=8), dimension(n) :: W
   integer :: devinfo

   real(kind=8), allocatable, device :: work(:)
   integer:: lwork
   integer async_id

   integer(kind=cuda_stream_kind) :: stream
   type(cusolverDnHandle) handle
   integer istat

   CUSOLVERCHECK(cusolverDnCreate(handle))
   stream = acc_get_cuda_stream(async_id)
   CUSOLVERCHECK(cusolverDnSetStream(handle, stream))

   !$acc host_data use_device(A, W)
   CUSOLVERCHECK(cusolverDnDsyevd_bufferSize(handle, jobz, uplo, n, A, lda, W, lwork))
   !$acc end host_data

   CUDACHECK(cudaMallocAsync(work, lwork, stream))

   !$acc host_data use_device(A, W, devinfo)
   CUSOLVERCHECK(cusolverDnDsyevd(handle, jobz, uplo, n, A, lda, W, work, lwork, devinfo))
   !$acc end host_data
   CUDACHECK(cudaFreeAsync(work, stream))

end subroutine DnDsyevd_Async
! ============================================================
subroutine DnDsyevj_async(handle, jobz, uplo, n, A, lda, W, &
                          devinfo, params, async_id)
   use openacc
   use cudafor
   use cusolverDN
   implicit none
   type(cusolverDnHandle) :: handle
   integer :: jobz, uplo
   integer :: n, lda
   real(kind=8), dimension(lda, *) :: A
   real(kind=8), dimension(n) :: W
   integer :: devinfo
   type(cusolverDnSyevjInfo) :: params

   real(kind=8), allocatable, device :: work(:)
   integer:: lwork
   integer async_id

   integer(kind=cuda_stream_kind) :: stream
   integer istat

   stream = acc_get_cuda_stream(async_id)
   CUSOLVERCHECK(cusolverDnSetStream(handle, stream))

   !$acc host_data use_device(A, W)
   CUSOLVERCHECK(cusolverDnDsyevj_bufferSize(handle, jobz, uplo, n, A, lda, W, lwork, params))
   !$acc end host_data

   CUDACHECK(cudaMallocAsync(work, lwork, stream))

   !$acc host_data use_device(A, W, devinfo)
   CUSOLVERCHECK(cusolverDnDsyevj(handle, jobz, uplo, n, A, lda, W, work, lwork, devinfo, params))
   !$acc end host_data
   CUDACHECK(cudaFreeAsync(work, stream))

end subroutine DnDsyevj_Async
