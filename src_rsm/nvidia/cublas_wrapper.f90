!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dgemm_async(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, async_id)
    use cudafor
    use cublas
    use openacc
    implicit none
    character*1 :: transa, transb
    integer :: m, n, k, lda, ldb, ldc
    real(8), dimension(lda, *) :: a
    real(8), dimension(ldb, *) :: b
    real(8), dimension(ldc, *) :: c
    real(8) :: alpha, beta
    integer :: async_id
    integer(kind=cuda_stream_kind) :: stream
    integer :: istat
    type(cublashandle) :: handle

    handle = cublasGetHandle()
    stream = acc_get_cuda_stream(async_id)
    istat = cublasSetStream(handle, stream)

    !$acc host_data use_device(a, b, c)
    call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    !$acc end host_data
end subroutine dgemm_async

subroutine dgemm_gpu(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    implicit none
    character*1 :: transa, transb
    integer :: m, n, k, lda, ldb, ldc
    real(8), dimension(lda, *) :: a
    real(8), dimension(ldb, *) :: b
    real(8), dimension(ldc, *) :: c
    real(8) :: alpha, beta
    integer :: a_sd, b_sd

    if (transa == 'n') then
        a_sd = k
    else
        a_sd = m
    end if
    if (transb == 'n') then
        b_sd = n
    else
        b_sd = k
    end if

    !$acc data copyin(a(:lda, :a_sd), b(:ldb, :b_sd)) copy(c(:ldc, :n))
    call dgemm_async(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, 1)
    !$acc wait(1)
    !$acc end data
end subroutine dgemm_gpu


