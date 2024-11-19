!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine joinrs_gpu(cc, r1, r2, r3, r4, nx, my_max, lev, jlistnum, num, ncld)
   ! Present on device: cc, r1, r2, r3, r4, jlist1, nxjlen, nxjlen_all

   use const, only: RTYPE

   implicit none
   real(kind=RTYPE) cc(*), r1(*), r2(*), r3(*), r4(*)
   integer nx, my_max, lev, jlistnum, num, ncld

   if (num .eq. 1) call join1rs_gpu(cc, r1, nx, my_max, lev, jlistnum, ncld)
   if (num .eq. 2) call join2rs_gpu(cc, r1, r2, nx, my_max, lev, jlistnum, ncld)
   if (num .eq. 3) call join3rs_gpu(cc, r1, r2, r3, nx, my_max, lev, jlistnum, ncld)

   return
end

subroutine join1rs_gpu(cc, r1, nx, my_max, lev, jnum, ncld)
   use index
   use const, only: RTYPE
   use param, only: npex
   use openacc
   use cudafor

   implicit none
   real(kind=RTYPE) cc(nx + 2, levp, ncld, my_max)
   real(kind=RTYPE) r1(nxp, lev*ncld, my_max)
   real(kind=RTYPE) bufA(nxp, lev, ncld, my_max)
   real(kind=RTYPE) bufB(nx, levp, ncld, my_max)
   integer nx, my_max, lev, jnum, ncld
   integer jj, j, nxj, k, i, n, nk, kk
   integer async_id, istat
   integer(kind=cuda_stream_kind) stream

   async_id = 1
   stream = acc_get_cuda_stream(async_id)
   if (npex .eq. 1) then
      !$acc parallel loop collapse(4) private(j, nxj, kk) async(async_id)
      do jj = 1, jlistnum
         do n = 1, ncld
            do k = 1, lev
               do i = 1, nxp
                  j = jlist1(jj)
                  nxj = nxdef(j)
                  if (i .le. nxj) then
                     kk = (n - 1)*lev + k
                     cc(i, k, n, jj) = r1(i, kk, jj)
                  end if
               end do
            end do
         end do
      end do

   else
      !$acc enter data create(bufA, bufB) async(async_id)

      !$acc host_data use_device(bufA, bufB)
      istat = cudaMemsetAsync(bufA, real(0.0, RTYPE), size(bufA), stream)
      istat = cudaMemsetAsync(bufB, real(0.0, RTYPE), size(bufB), stream)
      !$acc end host_data

      !$acc parallel loop collapse(4) private(kk) async(async_id)
      do jj = 1, jlistnum
         do n = 1, ncld
            do k = 1, lev
               do i = 1, nxp
                  kk = (n - 1)*lev + k
                  bufA(i, k, n, jj) = r1(i, kk, jj)
               end do
            end do
         end do
      end do

      ! Present on device: bufA, bufB, jlist1, nxjlen, nxjlen_all
      call mpe2d_transpose_nxp_lev_gpu(bufA, bufB, nxp, nx, lev, levp, ncld, myf, my_max, jlistnum, jlen, nsizex, nccl_row_comm)

      !$acc parallel loop collapse(4) async(async_id)
      do jj = 1, jlistnum
         do n = 1, ncld
            do k = 1, levp
               do i = 1, nx
                  cc(i, k, n, jj) = bufB(i, k, n, jj)
               end do
            end do
         end do
      end do

      !$acc exit data delete(bufA, bufB) async(async_id)

   end if
   return
end

subroutine join2rs_gpu(cc, r1, r2, nx, my_max, lev, jnum, ncld)

   use index
   use const, only: RTYPE
   use param, only: npex
   use openacc
   use cudafor

   implicit none
   real(kind=RTYPE) cc(nx + 2, levp, 1 + ncld, my_max) ! Present on device
   real(kind=RTYPE) r1(nxp, lev, my_max) ! Present on device
   real(kind=RTYPE) r2(nxp, lev*ncld, my_max) ! Present on device
   real(kind=RTYPE) bufA(nxp, lev, 1 + ncld, my_max)
   real(kind=RTYPE) bufB(nx, levp, 1 + ncld, my_max)
   integer nx, my_max, lev, jnum, ncld
   integer jj, j, nxj, k, i, n, nk, kk, nxj
   integer async_id, istat
   integer(kind=cuda_stream_kind) stream

   async_id = 1
   stream = acc_get_cuda_stream(async_id)
   if (npex .eq. 1) then
      !$acc parallel loop collapse(3) private(j, nxj) async(async_id)
      do jj = 1, jlistnum
         do k = 1, lev
            do i = 1, nxp
               j = jlist1(jj)
               nxj = nxdef(j)
               if (i .le. nxj) then
                  cc(i, k, 1, jj) = r1(i, k, jj)
               end if
            end do
         end do
      end do

      !$acc parallel loop gang collapse(4) private(j, nxj, nk, kk) async(async_id)
      do jj = 1, jlistnum
         do n = 1, ncld
            do k = 1, lev
               do i = 1, nxp
                  j = jlist1(jj)
                  nxj = nxdef(j)
                  if (i .le. nxj) then
                     nk = (n - 1)*lev
                     kk = nk + k
                     cc(i, k, 1 + n, jj) = r2(i, kk, jj)
                  end if
               end do
            end do
         end do
      end do
   else
      !$acc enter data create(bufA, bufB) async(async_id)

      !$acc host_data use_device(bufA, bufB)
      istat = cudaMemsetAsync(bufA, real(0.0, RTYPE), size(bufA), stream)
      istat = cudaMemsetAsync(bufB, real(0.0, RTYPE), size(bufB), stream)
      !$acc end host_data

      !$acc parallel loop collapse(3) async(async_id)
      do jj = 1, jlistnum
         do k = 1, lev
            do i = 1, nxp
               bufA(i, k, 1, jj) = r1(i, k, jj)
            end do
         end do
      end do

      !$acc parallel loop gang collapse(4) private(nk, kk) async(async_id)
      do jj = 1, jlistnum
         do n = 1, ncld
            do k = 1, lev
               do i = 1, nxp
                  nk = (n - 1)*lev
                  kk = nk + k
                  bufA(i, k, 1 + n, jj) = r2(i, kk, jj)
               end do
            end do
         end do
      end do

      ! Present on device: bufA, bufB, jlist1, nxjlen, nxjlen_all
      call mpe2d_transpose_nxp_lev_gpu(bufA, bufB, nxp, nx, lev, levp, 1 + ncld, myf, my_max, jlistnum, jlen, nsizex, nccl_row_comm)

      !$acc parallel loop collapse(4) async(async_id)
      do jj = 1, jlistnum
         do n = 1, 1 + ncld
            do k = 1, levp
               do i = 1, nx
                  cc(i, k, n, jj) = bufB(i, k, n, jj)
               end do
            end do
         end do
      end do

      !$acc exit data delete(bufA, bufB) async(async_id)
   end if
   return
end

subroutine join3rs_gpu(cc, r1, r2, r3, nx, my_max, lev, jnum, ncld)

   use index
   use const, only: RTYPE

   implicit none
   real(kind=RTYPE) cc(nx + 2, levp, 2 + ncld, my_max)
   real(kind=RTYPE) r1(nxp, lev, my_max)
   real(kind=RTYPE) r2(nxp, lev, my_max)
   real(kind=RTYPE) r3(nxp, lev*ncld, my_max)
   real(kind=RTYPE) bufA(nxp, lev, 2 + ncld, my_max)
   real(kind=RTYPE) bufB(nx, levp, 2 + ncld, my_max)
   integer nx, my_max, lev, jnum, ncld
   integer jj, j, nxj, k, i, n, nk, kk

   bufA = 0.
   bufB = 0.

   do jj = 1, jlistnum
   do k = 1, lev
   do i = 1, nxp
      bufA(i, k, 1, jj) = r1(i, k, jj)
      bufA(i, k, 2, jj) = r2(i, k, jj)
   end do
   end do
   end do

   do jj = 1, jlistnum
   do n = 1, ncld
      nk = (n - 1)*lev
      do k = 1, lev
         kk = nk + k
         do i = 1, nxp
            bufA(i, k, 2 + n, jj) = r3(i, kk, jj)
         end do
      end do
   end do
   end do

   call mpe2d_transpose_nxp_lev_gpu(bufA, bufB, nxp, nx, lev, levp, 2 + ncld, myf, my_max, jlistnum, jlen, nsizex, nccl_row_comm)

   do jj = 1, jlistnum
   do n = 1, 2 + ncld
   do k = 1, levp
   do i = 1, nx
      cc(i, k, n, jj) = bufB(i, k, n, jj)
   end do
   end do
   end do
   end do

   return
end
