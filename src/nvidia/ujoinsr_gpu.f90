subroutine ujoinsr_gpu(cc, r1, r2, r3, r4, nx, my_max, lev, jlistnum, num, ncld)
   ! Present on device: cc, r1, r2, r3, r4, jlist1
   use const, only: RTYPE

   implicit none
   real(kind=RTYPE) cc(*), r1(*), r2(*), r3(*), r4(*)
   integer nx, my_max, lev, jlistnum, num, ncld

   if (num .eq. 1) call ujoin1sr_gpu(cc, r1, nx, my_max, lev &
                                     , jlistnum, ncld)
   if (num .eq. 2) call ujoin2sr_gpu(cc, r1, r2, nx, my_max, lev &
                                     , jlistnum, ncld)
   if (num .eq. 3) call ujoin3sr_gpu(cc, r1, r2, r3, nx, my_max, lev &
                                     , jlistnum, ncld)
   if (num .eq. 4) call ujoin4sr_gpu(cc, r1, r2, r3, r4, nx, my_max, lev &
                                     , jlistnum, ncld)

   return
end

subroutine ujoin1sr_gpu(cc, r1, nx, my_max, lev, jnum, ncld)

   use index
   use const, only: RTYPE
   use param, only: npex
   use openacc
   use cudafor

   implicit none
   real(kind=RTYPE) cc(nx + 2, levp, ncld, my_max) ! Present on device
   real(kind=RTYPE) r1(nxp, lev*ncld, my_max) ! Present on device
   real(kind=RTYPE) bufA(nx, levp, ncld, my_max)
   real(kind=RTYPE) bufB(nxp, lev, ncld, my_max)
   integer nx, my_max, lev, jnum, ncld
   integer jj, j, nxj, k, i, nk, kk, n
   integer async_id, istat
   integer(kind=cuda_stream_kind) :: stream

   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   !$acc host_data use_device(r1)
   istat = cudaMemSetAsync(r1, real(0.0, RTYPE), size(r1), stream)
   !$acc end host_data
   if (npex .eq. 1) then
      !$acc parallel loop gang collapse(4) private(j, nxj, kk) async(async_id)
      do jj = 1, jlistnum
         do n = 1, ncld
            do k = 1, lev
               do i = 1, nxp
                  j = jlist1(jj)
                  nxj = nxdef(j)
                  if (i .le. nxj) then
                     kk = (n - 1)*lev + k
                     r1(i, kk, jj) = cc(i, k, n, jj)
                  end if
               end do
            end do
         end do
      end do
   else
      !$acc enter data create(bufA, bufB) async(async_id)

      !$acc host_data use_device(bufA, bufB, r1)
      istat = cudaMemSetAsync(bufA, real(0.0, RTYPE), size(bufA), stream)
      istat = cudaMemSetAsync(bufB, real(0.0, RTYPE), size(bufB), stream)
      !$acc end host_data

      !$acc parallel loop collapse(4) async(async_id)
      do jj = 1, jlistnum
      do n = 1, ncld
      do k = 1, levp
      do i = 1, nx
         bufA(i, k, n, jj) = cc(i, k, n, jj)
      end do
      end do
      end do
      end do

      ! Present on device: jlist1, bufA, bufB
      call mpe2d_transpose_nx_levp_gpu(bufA, bufB, nxp, nx, lev, levp, ncld, myf, my_max, jlistnum, jlen, nsizex, nccl_row_comm)

      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, jlistnum
      do n = 1, ncld
         nk = (n - 1)*lev
         !$acc loop worker
         do k = 1, lev
            kk = nk + k
            !$acc loop vector
            do i = 1, nxp
               r1(i, kk, jj) = bufB(i, k, n, jj)
            end do
         end do
      end do
      end do
      !$acc exit data delete(bufA, bufB) async(async_id)
   end if

   return
end

subroutine ujoin2sr_gpu(cc, r1, r2, nx, my_max, lev, jnum, ncld)

   use index
   use const, only: RTYPE
   use param, only: npex
   use openacc
   use cudafor

   implicit none
   real(kind=RTYPE) cc(nx + 2, levp, 1 + ncld, my_max)
   real(kind=RTYPE) r1(nxp, lev, my_max)
   real(kind=RTYPE) r2(nxp, lev*ncld, my_max)
   real(kind=RTYPE) bufA(nx, levp, 1 + ncld, my_max)
   real(kind=RTYPE) bufB(nxp, lev, 1 + ncld, my_max)
   integer nx, my_max, lev, jnum, ncld
   integer jj, j, nxj, k, i, nk, kk, n
   integer async_id, istat
   integer(kind=cuda_stream_kind) :: stream

   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   !$acc host_data use_device(r1, r2)
   istat = cudaMemSetAsync(r1, real(0.0, RTYPE), size(r1), stream)
   istat = cudaMemSetAsync(r2, real(0.0, RTYPE), size(r2), stream)
   !$acc end host_data

   if (npex .eq. 1) then
      !$acc parallel loop collapse(3) private(j, nxj) async(async_id)
      do jj = 1, jlistnum
         do k = 1, lev
            do i = 1, nxp
               j = jlist1(jj)
               nxj = nxdef(j)
               if (i .le. nxj) then
                  r1(i, k, jj) = cc(i, k, 1, jj)
               end if
            end do
         end do
      end do

      !$acc parallel loop collapse(4) private(j, nxj, kk) async(async_id)
      do jj = 1, jlistnum
         do n = 1, ncld
            do k = 1, lev
               do i = 1, nxp
                  j = jlist1(jj)
                  nxj = nxdef(j)
                  if (i .le. nxj) then
                     kk = (n - 1)*lev + k
                     r2(i, kk, jj) = cc(i, k, 1 + n, jj)
                  end if
               end do
            end do
         end do
      end do
   else
      !$acc enter data create(bufA, bufB) async(async_id)

      !$acc host_data use_device(bufA, bufB)
      istat = cudaMemSetAsync(bufA, real(0.0, RTYPE), size(bufA), stream)
      istat = cudaMemSetAsync(bufB, real(0.0, RTYPE), size(bufB), stream)
      !$acc end host_data

      !$acc parallel loop collapse(4) async(async_id)
      do jj = 1, jlistnum
         do n = 1, 1 + ncld
            do k = 1, levp
               do i = 1, nx
                  bufA(i, k, n, jj) = cc(i, k, n, jj)
               end do
            end do
         end do
      end do

      ! Present on device: bufA, bufB, jlist1
      call mpe2d_transpose_nx_levp_gpu(bufA, bufB, nxp, nx, lev, levp, 1 + ncld, myf, my_max, jlistnum, jlen, nsizex, nccl_row_comm)

      !$acc parallel loop collapse(3) async(async_id)
      do jj = 1, jlistnum
         do k = 1, lev
            do i = 1, nxp
               r1(i, k, jj) = bufB(i, k, 1, jj)
            end do
         end do
      end do

      !$acc parallel loop collapse(2) gang async(async_id)
      do jj = 1, jlistnum
         do n = 1, ncld
            nk = (n - 1)*lev
            !$acc loop worker
            do k = 1, lev
               kk = nk + k
               !$acc loop vector
               do i = 1, nxp
                  r2(i, kk, jj) = bufB(i, k, 1 + n, jj)
               end do
            end do
         end do
      end do
      !$acc exit data delete(bufA, bufB) async(async_id)

   end if
   return
end

subroutine ujoin3sr_gpu(cc, r1, r2, r3, nx, my_max, lev, jnum, ncld)

   use index
   use const, only: RTYPE
   use openacc
   use cudafor

   implicit none
   real(kind=RTYPE) cc(nx + 2, levp, 2 + ncld, my_max) ! Present on device
   real(kind=RTYPE) r1(nxp, lev, my_max) ! Present on device
   real(kind=RTYPE) r2(nxp, lev, my_max) ! Present on device
   real(kind=RTYPE) r3(nxp, lev*ncld, my_max) ! Present on device
   real(kind=RTYPE) bufA(nx, levp, 2 + ncld, my_max)
   real(kind=RTYPE) bufB(nxp, lev, 2 + ncld, my_max)
   integer nx, my_max, lev, jnum, ncld
   integer jj, j, nxj, k, i, nk, kk, n
   integer async_id, istat
   integer(kind=cuda_stream_kind) :: stream

   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   !$acc enter data create(bufA, bufB) async(async_id)

   !$acc host_data use_device(bufA, bufB, r1)
   istat = cudaMemSetAsync(bufA, real(0.0, RTYPE), size(bufA), stream)
   istat = cudaMemSetAsync(bufB, real(0.0, RTYPE), size(bufB), stream)
   istat = cudaMemSetAsync(r1, real(0.0, RTYPE), size(r1), stream)
   !$acc end host_data

   !$acc parallel loop collapse(4) async(async_id)
   do jj = 1, jlistnum
   do n = 1, 2 + ncld
   do k = 1, levp
   do i = 1, nx
      bufA(i, k, n, jj) = cc(i, k, n, jj)
   end do
   end do
   end do
   end do

   call mpe2d_transpose_nx_levp_gpu(bufA, bufB, nxp, nx, lev, levp, 2 + ncld, myf, my_max, jlistnum, jlen, nsizex, nccl_row_comm)

   !$acc parallel loop collapse(3) async(async_id)
   do jj = 1, jlistnum
   do k = 1, lev
   do i = 1, nxp
      r1(i, k, jj) = bufB(i, k, 1, jj)
      r2(i, k, jj) = bufB(i, k, 2, jj)
   end do
   end do
   end do

   !$acc parallel loop gang collapse(2) async(async_id)
   do jj = 1, jlistnum
   do n = 1, ncld
      nk = (n - 1)*lev
      !$acc loop worker
      do k = 1, lev
         kk = nk + k
         !$acc loop vector
         do i = 1, nxp
            r3(i, kk, jj) = bufB(i, k, 2 + n, jj)
         end do
      end do
   end do
   end do
   !$acc exit data delete(bufA, bufB) async(async_id)

   return
end

subroutine ujoin4sr_gpu(cc, r1, r2, r3, r4, nx, my_max, lev, jnum, ncld)

   use index
   use const, only: RTYPE

   implicit none
   real(kind=RTYPE) cc(nx + 2, levp, 3 + ncld, my_max) ! Present on device
   real(kind=RTYPE) r1(nxp, lev, my_max) ! Present on device
   real(kind=RTYPE) r2(nxp, lev, my_max) ! Present on device
   real(kind=RTYPE) r3(nxp, lev, my_max) ! Present on device
   real(kind=RTYPE) r4(nxp, lev*ncld, my_max) ! Present on device
   real(kind=RTYPE) bufA(nx, levp, 3 + ncld, my_max)
   real(kind=RTYPE) bufB(nxp, lev, 3 + ncld, my_max)
   integer nx, my_max, lev, jnum, ncld
   integer jj, j, nxj, k, i, nk, kk, n
   integer async_id

   async_id = 1

   !$acc enter data create(bufA, bufB) async(async_id)
   !$acc parallel loop collapse(4) async(async_id)
   do jj = 1, jlistnum
   do n = 1, 3 + ncld
   do k = 1, levp
   do i = 1, nx
      bufA(i, k, n, jj) = cc(i, k, n, jj)
   end do
   end do
   end do
   end do

   call mpe2d_transpose_nx_levp_gpu(bufA, bufB, nxp, nx, lev, levp, 3 + ncld, myf, my_max, jlistnum, jlen, nsizex, nccl_row_comm)

   !$acc parallel loop collapse(3) async(async_id)
   do jj = 1, jlistnum
   do k = 1, lev
   do i = 1, nxp
      r1(i, k, jj) = bufB(i, k, 1, jj)
      r2(i, k, jj) = bufB(i, k, 2, jj)
      r3(i, k, jj) = bufB(i, k, 3, jj)
   end do
   end do
   end do

   !$acc parallel loop gang collapse(2) async(async_id)
   do jj = 1, jlistnum
   do n = 1, ncld
      nk = (n - 1)*lev
      !$acc loop worker
      do k = 1, lev
         kk = nk + k
         !$acc loop vector
         do i = 1, nxp
            r4(i, kk, jj) = bufB(i, k, 3 + n, jj)
         end do
      end do
   end do
   end do
   !$acc exit data delete(bufA, bufB) async(async_id)

   return
end
