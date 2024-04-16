!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine mpe2d_transpose_nx_levp_gpu(ain, aout, nxp, nx, lev, levp, num, my, my_max, jlistnum, jlen, nsizex, comm)

! transpose (nx full,lev partial) to (nx partial,lev full), num variables packed

   use index, only: jlist1, nxjlen_all
   use const, only: RTYPE, MPI_RTYPE
   use openacc
   use cudafor

   implicit none

   include 'mpif.h'
   integer nxp, nx, lev, levp, my, my_max, jlen, nsizex, comm
   real(kind=RTYPE) ain(nx, levp, num, my_max), aout(nxp, lev, num, my_max) ! Present on device
   real(kind=RTYPE) b1(levp, num, jlen, nxp, nsizex), b2(levp, num, jlen, nxp, nsizex)
   integer nlen, j, jj, i, k, KL, ierr, jlistnum, num, n, i1, i2, j1
   integer async_id, istat
   integer(kind=cuda_stream_kind) :: stream

   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   !$acc enter data create(b1, b2) copyin(jlist1) async(async_id)

   !$acc host_data use_device(b1, b2, aout)
   istat = cudaMemSetAsync(b1, 0.0, size(b1), stream)
   istat = cudaMemSetAsync(b2, 0.0, size(b2), stream)
   istat = cudaMemSetAsync(aout, 0.0, size(aout), stream)
   !$acc end host_data

   !$acc parallel loop gang async(async_id)
   do j = 1, jlistnum
      j1 = jlist1(j)
      !$acc loop vector collapse(2)
      do n = 1, num
      do k = 1, levp
         i1 = 1
         do i = 1, nsizex
            i2 = nxjlen_all(i, j1)
            b1(k, n, j, 1:i2, i) = ain(i1:i1 + i2 - 1, k, n, j)
            i1 = i1 + i2
         end do
      end do
      end do
   end do

   nlen = nxp*levp*jlen*num
   !$acc wait(async_id)
   !$acc host_data use_device(b1, b2)
   call MPI_ALLTOALL(b1, nlen, MPI_RTYPE, &
                     b2, nlen, MPI_RTYPE, &
                     comm, IERR)
   !$acc end host_data

   !$acc parallel loop collapse(2) async(async_id)
   do jj = 1, jlistnum
   do n = 1, num
      !$acc loop seq
      do i = 1, nxp
         k = 1
         do j = 1, nsizex
            aout(i, k:k + levp - 1, n, jj) = b2(1:levp, n, jj, i, j)
            k = k + levp
         end do
      end do
   end do
   end do
   !$acc exit data delete(b1, b2, jlist1) async(async_id)

   return
end

subroutine mpe2d_transpose_nxp_lev_gpu(ain, aout, nxp, nx, lev, levp, num, my, my_max, jlistnum, jlen, nsizex, comm)

! transpose (nx partial,lev full) to (nx full,lev partial), num variable packed

   use index, only: jlist1, nxjlen_all, nxjlen
   use const, only: RTYPE, MPI_RTYPE
   use openacc
   use cudafor

   implicit none

   include 'mpif.h'
   integer nx, nxp, lev, levp, my, my_max, jlen, nsizex, comm
   real(kind=RTYPE) ain(nxp, lev, num, my_max), aout(nx, levp, num, my_max) ! present on device
   real(kind=RTYPE) b1(nxp, jlen, num, lev), b2(nxp, jlen, num, levp, nsizex)
   integer nlen, j, i, k, ierr, jlistnum, num, n, i1, i2, j1

   integer async_id, istat
   integer(kind=cuda_stream_kind) :: stream

   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   !$acc enter data copyin(jlist1, nxjlen, nxjlen_all) create(b1, b2) async(async_id)

   !$acc host_data use_device(b1, b2, aout)
   istat = cudaMemSetAsync(b1, 0.0, size(b1), stream)
   istat = cudaMemSetAsync(b2, 0.0, size(b2), stream)
   istat = cudaMemSetAsync(aout, 0.0, size(aout), stream)
   !$acc end host_data

   !$acc parallel loop gang collapse(2) async(async_id)
   do k = 1, lev
   do n = 1, num
      !$acc loop worker
      do j = 1, jlistnum
         j1 = jlist1(j)
         !$acc loop vector
         do i = 1, nxjlen(j1)
            b1(i, j, n, k) = ain(i, k, n, j)
         end do
      end do
   end do
   end do
   !$acc wait(async_id)

   nlen = nxp*levp*jlen*num
   !$acc host_data use_device(b1, b2)
   call MPI_ALLTOALL(b1, nlen, MPI_RTYPE, &
                     b2, nlen, MPI_RTYPE, &
                     comm, IERR)
   !$acc end host_data

   !$acc kernels async(async_id)
   !$acc loop collapse(3) private(j1, i1, i2)
   do n = 1, num
   do k = 1, levp
   do j = 1, jlistnum
      j1 = jlist1(j)
      i1 = 1
      do i = 1, nsizex
         i2 = nxjlen_all(i, j1)
         aout(i1:i1 + i2 - 1, k, n, j) = b2(1:i2, j, n, k, i)
         i1 = i1 + i2
      end do
   end do
   end do
   end do
   !$acc end kernels
   !$acc exit data delete(jlist1, nxjlen, nxjlen_all, b1, b2) async(async_id)

   return
end
