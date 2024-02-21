!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_dgemm_wrapper
        use cudafor
        integer :: driver_ver, runtime_ver, err
        err = cudaDriverGetVersion(driver_ver)
        err = cudaRuntimeGetVersion(runtime_ver)
        PRINT *, driver_ver
        PRINT *, runtime_ver
!   call cublas_dgemm_unit('n', 'n', 32, 64, 128, 32, 128, 128)
!   call cublas_dgemm_unit('n', 'n', 24, 32, 80, 30, 120, 34)
end program test_dgemm_wrapper

subroutine cublas_dgemm_unit(transa, transb, m, n, k, lda, ldb, ldc)
   character*1 :: transa, transb
   integer :: m, n, k, lda, ldb, ldc
   real(8) :: alpha, beta
   real(8), dimension(:, :), allocatable :: a, b
   real(8), dimension(ldc, n) :: c, c_cublas

   if (transa == 'n') then
      allocate (a(lda, k))
   else
      allocate (a(lda, m))
   end if
   if (transb == 'n') then
      allocate (b(ldb, n))
   else
      allocate (b(ldb, k))
   end if

   call random_seed()
   call random_number(a)
   call random_number(b)
   call random_number(c)
   c_cublas = c

   call random_number(alpha)
   call random_number(beta)

   call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
   call dgemm_gpu(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c_cublas, ldc)

   if (all(abs(c - c_cublas) <= 1e-10)) then
      PRINT *, "test_dgemm_wrapper passed."
   else
      PRINT *, "test_dgemm_wrapper failed."
      call exit(1)
   end if

end subroutine cublas_dgemm_unit
