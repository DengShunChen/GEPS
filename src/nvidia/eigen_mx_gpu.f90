#define CUSOLVERCHECK(ierr) call cuSOLVER_check_helper(ierr, __FILE__, __LINE__)

subroutine eigen_mx_gpu(mx, eval, nnlist, no)
   use const, only: nnmivm
   use param, only: jtrun, jtmax
   use index, only: mlistnum
   use cusolverDn
   implicit none
   real(kind=8), intent(inout):: mx(no*no, 2*jtmax*nnmivm)
   real(kind=8), intent(out):: eval(no, 2*jtmax*nnmivm)
   integer, intent(in):: nnlist(2*jtmax*nnmivm)
   integer, intent(in):: no

   integer devinfo(2*jtmax*nnmivm)
   integer jobz, uplo
   type(cusolverDnSyevjInfo) :: params
   real(kind=8):: tol
   integer :: MaxSweeps, SortEig
   integer i, j, m, L, nn, k
   integer async_id
   async_id = 1

   jobz = CUSOLVER_EIG_MODE_VECTOR ! compute eigenvectors
   uplo = CUBLAS_FILL_MODE_UPPER

   ! ------------------------------------------------------------
   ! << for DnDsyevj >>
   tol = 1e-15
   ! MaxSweeps = 15
   SortEig = 0                      ! don't sort eigenvalues
   CUSOLVERCHECK(cusolverDnCreateSyevjInfo(params))
   CUSOLVERCHECK(cusolverDnXsyevjSetTolerance(params, tol))
   ! CUSOLVERCHECK(cusolverDnXsyevjSetMaxSweeps(params, MaxSweeps))
   CUSOLVERCHECK(cusolverDnXsyevjSetSortEig(params, SortEig))
   ! ------------------------------------------------------------

   !$acc enter data create(devinfo) async(async_id)
   do L = 1, nnmivm
      do m = 1, mlistnum
         do i = 1, 2
            j = i + (m - 1)*2 + (L - 1)*jtmax*2
            nn = nnlist(j)

            !$acc wait(async_id) async(j+1)
            ! <<<< -------------------------------------------------------
            ! call DnDsyevd_Async(jobz, uplo, nn, mx(1, j), nn, eval(1, j), &
            !      devinfo(j), j+1)
            ! ------------------------------------------------------------
            call DnDsyevj_Async(jobz, uplo, nn, mx(1, j), nn, eval(1, j), &
                                devinfo(j), params, j + 1)
            ! ------------------------------------------------------- >>>>

         end do
      end do
   end do

   do L = 1, nnmivm
      do m = 1, mlistnum
         do i = 1, 2
            j = i + (m - 1)*2 + (L - 1)*jtmax*2
            !$acc wait(j+1) async(async_id)
         end do
      end do
   end do
   !$acc exit data delete(devinfo) async(async_id)

end subroutine eigen_mx_gpu
