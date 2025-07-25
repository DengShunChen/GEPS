!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define CUSOLVERCHECK(ierr) call cuSOLVER_check_helper(ierr, __FILE__, __LINE__)

subroutine eigen_mx_gpu(mx, eval, nnlist, no)
   use const, only: nnmivm
   use param, only: jtrun, jtmax
   use index, only: mlistnum
   use cusolverDn

   implicit none

   integer, parameter :: max_handles = 8

   real(kind=8), intent(inout) :: mx(no*no, 2*jtmax*nnmivm)
   real(kind=8), intent(out)   :: eval(no, 2*jtmax*nnmivm)
   integer, intent(in)         :: nnlist(2*jtmax*nnmivm)
   integer, intent(in)         :: no

   integer                   :: devinfo(2*jtmax*nnmivm)
   integer                   :: jobz, uplo
   type(cusolverDnSyevjInfo) :: params
   real(kind=8)              :: tol
   integer                   :: SortEig
   integer                   :: i, j, m, L, nn, handle_id
   integer                   :: async_id
   type(cusolverDnHandle)    :: handle(max_handles)

   async_id = 1

   ! Create cuSOLVER handles
   do i = 1, max_handles
      CUSOLVERCHECK(cusolverDnCreate(handle(i)))
   end do

   ! Configure Syevj Solver
   jobz = CUSOLVER_EIG_MODE_VECTOR ! Compute eigenvectors
   uplo = CUBLAS_FILL_MODE_UPPER
   tol = 1e-15
   SortEig = 0                     ! Don't sort eigenvalues
   CUSOLVERCHECK(cusolverDnCreateSyevjInfo(params))
   CUSOLVERCHECK(cusolverDnXsyevjSetTolerance(params, tol))
   CUSOLVERCHECK(cusolverDnXsyevjSetSortEig(params, SortEig))

   !$acc enter data create(devinfo) async(async_id)

   ! Main Eigenvalue Decomposition Loop
   do L = 1, nnmivm
      do m = 1, mlistnum
         do i = 1, 2
            j = i + (m - 1)*2 + (L - 1)*jtmax*2
            nn = nnlist(j)
            handle_id = mod(j, max_handles) + 1
            call DnDsyevj_Async(handle(handle_id), jobz, uplo, nn, mx(1, j), &
                                 nn, eval(1, j), devinfo(j), params, handle_id)
         end do
      end do
   end do

   ! Wait for All Streams and Clean Up Handles
   do i = 1, max_handles
      !$acc wait(i)
      CUSOLVERCHECK(cusolverDnDestroy(handle(i)))
   end do

   !$acc exit data delete(devinfo) async(async_id)

end subroutine eigen_mx_gpu
