#define CUSOLVERCHECK(ierr) call cuSOLVER_check_helper(ierr, __FILE__, __LINE__)
subroutine DnDsyevd_unit
   use cudafor
   use cusolverDN
   implicit none
   integer, parameter:: lda = 3
   integer:: n
   real(kind=8):: a(lda, lda)
   real(kind=8):: w(lda), lambda(lda)
   integer :: devinfo
   integer :: jobz, uplo

   integer async_id
   async_id = 1
   n = lda

   jobz = CUSOLVER_EIG_MODE_VECTOR; ! compute eigenvectors
   uplo = CUBLAS_FILL_MODE_LOWER; 
   a(1:n, 1) = (/3.5, 0.5, 0/)
   a(1:n, 2) = (/0.5, 3.5, 0/)
   a(1:n, 3) = (/0.0, 0.0, 2/)
   lambda(1:n) = (/2, 3, 4/)

   !$acc enter data copyin(a) async(async_id)
   !$acc enter data create(w, devinfo) async(async_id)
   call DnDsyevd_Async(jobz, uplo, n, A, lda, W, &
                       devinfo, async_id)
   !$acc exit data copyout(a, w, devinfo) async(async_id)
   !$acc wait(async_id)
   print *, 'devinfo=', devinfo
   write (*, '(A, 3(1pe15.7,1X))') "W: ", w(1:n)

   call assert_allclose_r8(W, size(W), lambda, size(lambda), 1e-15, 1e-15, "W")

end subroutine DnDsyevd_unit

program test_DnDsyevj
   call DnDsyevd_unit
end program test_DnDsyevj
