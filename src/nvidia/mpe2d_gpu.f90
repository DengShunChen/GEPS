!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine mpe2d_transpose_nx_levp_gpu(ain, aout, nxp, nx, lev, levp, num, my, my_max, jlistnum, jlen, nsizex, comm)
! Present on device: ain, aout, jlist1
! transpose (nx full,lev partial) to (nx partial,lev full), num variables packed

   use index, only: jlist1, nxjlen_all
   use const, only: RTYPE, MPI_RTYPE
   use openacc
   use cudafor
   use mpi

   implicit none

   integer nxp, nx, lev, levp, my, my_max, jlen, nsizex, comm
   real(kind=RTYPE) ain(nx, levp, num, my_max), aout(nxp, lev, num, my_max)
   real(kind=RTYPE) b1(levp, num, jlen, nxp, nsizex), b2(levp, num, jlen, nxp, nsizex)
   integer nlen, j, jj, i, k, KL, ierr, jlistnum, num, n, i1, i2, j1
   integer async_id, istat
   integer(kind=cuda_stream_kind) :: stream

   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   !$acc enter data create(b1, b2) async(async_id)

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
   !$acc exit data delete(b1, b2) async(async_id)

   return
end

subroutine mpe2d_transpose_nxp_lev_gpu(ain, aout, nxp, nx, lev, levp, num, my, my_max, jlistnum, jlen, nsizex, comm)
   ! Present on device: ain, aout, jlist1, nxjlen, nxjlen_all
   ! transpose (nx partial,lev full) to (nx full,lev partial), num variable packed

   use index, only: jlist1, nxjlen_all, nxjlen
   use const, only: RTYPE, MPI_RTYPE
   use openacc
   use cudafor
   use mpi

   implicit none

   integer nx, nxp, lev, levp, my, my_max, jlen, nsizex, comm
   real(kind=RTYPE) ain(nxp, lev, num, my_max), aout(nx, levp, num, my_max)
   real(kind=RTYPE) b1(nxp, jlen, num, lev), b2(nxp, jlen, num, levp, nsizex)
   integer nlen, j, i, k, ierr, jlistnum, num, n, i1, i2, j1, ii

   integer async_id, istat
   integer(kind=cuda_stream_kind) stream
   integer i1_array(nsizex, jlistnum)

   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   !$acc enter data create(b1, b2) async(async_id)

   !$acc host_data use_device(b1, b2, aout)
   istat = cudaMemSetAsync(b1, 0.0, size(b1), stream)
   istat = cudaMemSetAsync(b2, 0.0, size(b2), stream)
   istat = cudaMemSetAsync(aout, 0.0, size(aout), stream)
   !$acc end host_data

   !$acc parallel loop collapse(4) private(j1) async(async_id)
   do k = 1, lev
      do n = 1, num
         do j = 1, jlistnum
            do i = 1, nxp
               j1 = jlist1(j)
               if (i .le. nxjlen(j1)) then
                  b1(i, j, n, k) = ain(i, k, n, j)
               end if
            end do
         end do
      end do
   end do

   do j = 1, jlistnum
      j1 = jlist1(j)
      i1 = 0
      do i = 1, nsizex
         i1_array(i, j) = i1
         i1 = i1 + nxjlen_all(i, j1)
      end do
   end do
   !$acc enter data copyin(i1_array) async(async_id)

   !$acc wait(async_id)

   nlen = nxp*levp*jlen*num
   !$acc host_data use_device(b1, b2)
   call MPI_ALLTOALL(b1, nlen, MPI_RTYPE, b2, nlen, MPI_RTYPE, comm, IERR)
   !$acc end host_data

   !$acc parallel loop collapse(5) private(i1, i2) async(async_id)
   do n = 1, num
      do k = 1, levp
         do j = 1, jlistnum
            do i = 1, nsizex
               do ii = 1, nxp
                  j1 = jlist1(j)
                  i2 = nxjlen_all(i, j1)
                  if (ii .le. i2) then
                     i1 = i1_array(i, j)
                     aout(i1 + ii, k, n, j) = b2(ii, j, n, k, i)
                  end if
               end do
            end do
         end do
      end do
   end do

   !$acc exit data delete(b1, b2, i1_array) async(async_id)

   return
end

subroutine mpe2d_reshape_pl_gpu(plin, plout)

! reshape spec coef

   use mpi
   use param
   use index
   use const, only: RTYPE
   use openacc
   use cudafor

   implicit none

   integer i, j, m, mf, nl
   integer i_array(mlistnum)
   integer jtsize

   real(kind=RTYPE) plin(jtrun, jtmax, 2) ! Present on device
   real(kind=RTYPE) plout(jtp, 2) ! Present on device
   real(kind=RTYPE) b1(jtf, 2)
   integer async_id, istat
   integer(kind=cuda_stream_kind) :: stream

   async_id = 1
   stream = acc_get_cuda_stream(async_id)
   !$acc enter data copyin(mlist) async(async_id)
   !$acc enter data create(b1) async(async_id)
   !$acc host_data use_device(plout, b1)
   istat = cudaMemSetAsync(plout, 0.0, size(plout), stream)
   istat = cudaMemSetAsync(b1, 0.0, size(b1), stream)
   !$acc end host_data
   i = 1
   do m = 1, mlistnum
      i_array(m) = i
      mf = mlist(m)
      i = i + jtrun - mf + 1
   end do
   !$acc enter data copyin(i_array) async(async_id)
   !$acc parallel loop async(async_id)
   do m = 1, mlistnum
      mf = mlist(m)
      NL = jtrun - mf + 1
      i = i_array(m)
      b1(i:i + nl - 1, 1) = plin(mf:jtrun, m, 1)
      b1(i:i + nl - 1, 2) = plin(mf:jtrun, m, 2)
   end do

   jtsize = jtend - jtstart + 1
   !$acc parallel loop async(async_id)
   do i = 1, jtsize
      plout(i, 1) = b1(jtstart + i - 1, 1)
      plout(i, 2) = b1(jtstart + i - 1, 2)
   end do
   !$acc exit data delete(mlist, b1, i_array) async(async_id)
   return
end

subroutine mpe2d_transpose_siimpl_gpu(ain, aout, &
                                      levp, jtrun, jtmax, lev, jtp, jtf, mlistnum, mlist, nsizex, row_comm)
   ! Present on device: ain, aout, mlist
! transpose siimpl spec

   use mpi
   use const, only: RTYPE, MPI_RTYPE
   use openacc
   use cudafor

   implicit none

   integer levp, jtrun, jtmax, lev, jtp, jtf, mlistnum, mlist(jtrun), nsizex, row_comm
   integer m, mf, nl, i, levp2, nlen, ierr, j, k

   real(kind=RTYPE) ain(levp, 2, jtrun, jtmax) ! Present on device
   real(kind=RTYPE) aout(lev, 2, jtp) ! Present on device
   real(kind=RTYPE) c1(levp, 2, jtf)
   real(kind=RTYPE) c2(levp, 2, jtp, nsizex)
   integer async_id, istat
   integer(kind=cuda_stream_kind) :: stream
   integer i_array(mlistnum)

   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   !$acc enter data create(c1, c2) async(async_id)
   !$acc host_data use_device(aout, c1, c2)
   istat = cudaMemSetAsync(aout, 0.0, size(aout), stream)
   istat = cudaMemSetAsync(c1, 0.0, size(c1), stream)
   istat = cudaMemSetAsync(c2, 0.0, size(c2), stream)
   !$acc end host_data

   levp2 = levp*2
   i = 1
   do m = 1, mlistnum
      i_array(m) = i
      mf = mlist(m)
      i = i + jtrun - mf + 1
   end do
   !$acc enter data copyin(i_array) async(async_id)

   !$acc parallel loop async(async_id)
   do m = 1, mlistnum
      mf = mlist(m)
      NL = jtrun - mf + 1
      i = i_array(m)
      c1(1:levp, 1, i:i + nl - 1) = ain(1:levp, 1, mf:jtrun, m)
      c1(1:levp, 2, i:i + nl - 1) = ain(1:levp, 2, mf:jtrun, m)
   end do

   nlen = levp2*jtp
   !$acc wait(async_id)
   !$acc host_data use_device(c1, c2)
   call MPI_ALLTOALL(c1, nlen, MPI_RTYPE, &
                     c2, nlen, MPI_RTYPE, &
                     row_comm, IERR)
   !$acc end host_data

   !$acc parallel loop collapse(2) async(async_id)
   do j = 1, jtp
      do i = 1, nsizex
         k = levp*(i - 1) + 1
         aout(k:k + levp - 1, 1, j) = c2(1:levp, 1, j, i)
         aout(k:k + levp - 1, 2, j) = c2(1:levp, 2, j, i)
      end do
   end do
   !$acc exit data delete(c1, c2, i_array) async(async_id)

   return
end
subroutine mpe2d_unify_nx_gpu(work, a)
   ! Present ont device: work, a, nxjlen_all
   ! unify a(nx_partial,my_partial) to work(nx_full,my_partial)

   use param
   use index
   use mpi
   use const, only: RTYPE, MPI_RTYPE
   use openacc
   use cudafor

   real(kind=RTYPE) work(nx, my_max)
   real(kind=RTYPE) a(nxp, my_max)
   real(kind=RTYPE) b(nxp, my_max, nsizex)
   integer async_id, istat
   integer(kind=cuda_stream_kind) stream
   integer ii_array(nsizex, jlistnum)

   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   !$acc enter data create(b) async(async_id)

   !$acc host_data use_device(work, b)
   istat = cudaMemsetAsync(work, 0.0, size(work), stream)
   istat = cudaMemsetAsync(b, 0.0, size(b), stream)
   !$acc end host_data

   do jj = 1, jlistnum
      j = jlist1(jj)
      ii = 0
      do i = 1, nsizex
         ii_array(i, jj) = ii
         ii = ii + nxjlen_all(i, j)
      end do
   end do
   !$acc enter data copyin(ii_array) async(async_id)

   !$acc wait(async_id)

   !$acc host_data use_device(a, b)
   call MPI_ALLGATHER(a, nxp*my_max, MPI_RTYPE, &
                      b, nxp*my_max, MPI_RTYPE, &
                      row_comm, IERR)
   !$acc end host_data

   !$acc parallel loop collapse(2) private(j, ii, nn) async(async_id)
   do jj = 1, jlistnum
      do i = 1, nsizex
         j = jlist1(jj)
         ii = ii_array(i, jj) + 1
         nn = nxjlen_all(i, j)
         work(ii:ii + nn - 1, jj) = b(1:nn, jj, 1)
      end do
   end do

   !$acc exit data delete(b, ii_array) async(async_id)

   return
end

subroutine mpe2d_reshape_pl_back_gpu(plin, plout)
   ! Present on device: plin, plout
   ! reshape spec coef back

   use mpi
   use param
   use index
   use const, only: RTYPE, MPI_RTYPE
   use openacc
   use cudafor

   implicit none

   integer n, i, j, m, mf, ierr, q

   real(kind=RTYPE) plin(jtp, 2)
   real(kind=RTYPE) plout(jtrun, jtmax, 2)
   real(kind=RTYPE) b2(jtp, 2, nsizex)
   integer i_array(jtrun, mlistnum)
   integer j_array(jtrun, mlistnum)
   integer async_id, istat
   integer(kind=cuda_stream_kind) stream

   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   !$acc enter data create(b2) async(async_id)

   !$acc host_data use_device(plout, b2)
   istat = cudaMemsetAsync(plout, 0.0, size(plout), stream)
   istat = cudaMemsetAsync(b2, 0.0, size(b2), stream)
   !$acc end host_data

   i = 1
   j = 1
   do m = 1, mlistnum
      mf = mlist(m)
      do n = mf, jtrun
         i_array(n, m) = i
         j_array(n, m) = j
         i = i + 1
         if (i .gt. jtlen_all(j)) then
            i = 1
            j = j + 1
         end if
      end do
   end do
   !$acc enter data copyin(i_array, j_array) async(async_id)

   !$acc wait(async_id)

   !$acc host_data use_device(plin, b2)
   call MPI_ALLGATHER(plin, jtp*2, MPI_RTYPE, &
                      b2, jtp*2, MPI_RTYPE, &
                      row_comm, IERR)
   !$acc end host_data

   !$acc parallel loop collapse(3) private(mf, i, j) async(async_id)
   do q = 1, 2
      do m = 1, mlistnum
         do n = 1, jtrun
            mf = mlist(m)
            if (n .ge. mf) then
               i = i_array(n, m)
               j = j_array(n, m)
               plout(n, m, q) = b2(i, q, j)
            end if
         end do
      end do
   end do

   !$acc exit data delete(b2, i_array, j_array) async(async_id)

   return
end

subroutine mpe2d_transpose_siimpl_back_gpu(ain, aout, levp, jtrun, jtmax, lev, jtp, jtf, mlistnum, mlist, nsizex, row_comm)
   ! Present on device: ain, aout
   ! transpose siimpl spec back

   use mpi
   use index, only: jtlen_all
   use const, only: RTYPE, MPI_RTYPE
   use openacc
   use cudafor

   implicit none

   integer levp, jtrun, jtmax, lev, jtp, jtf, mlistnum, mlist(jtrun), nsizex, row_comm
   integer m, n, mf, nl, i, nlen, ierr, j, k, q

   real(kind=RTYPE) ain(lev, 2, jtp)
   real(kind=RTYPE) aout(levp, 2, jtrun, jtmax)

   real(kind=RTYPE) b1(jtp, 2, lev)
   real(kind=RTYPE) b2(jtp, 2, levp, nsizex)
   integer i_array(jtrun, mlistnum), j_array(jtrun, mlistnum)
   integer async_id, istat
   integer(kind=cuda_stream_kind) stream

   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   !$acc enter data create(b1, b2) async(async_id)

   !$acc host_data use_device(b1, b2, aout)
   istat = cudaMemsetAsync(b1, 0.0, size(b1), stream)
   istat = cudaMemsetAsync(b2, 0.0, size(b2), stream)
   istat = cudaMemsetAsync(aout, 0.0, size(aout), stream)
   !$acc end host_data

   !$acc parallel loop collapse(3) async(async_id)
   do j = 1, jtp
      do q = 1, 2
         do k = 1, lev
            b1(j, q, k) = ain(k, q, j)
         end do
      end do
   end do

   i = 1
   j = 1
   do m = 1, mlistnum
      mf = mlist(m)
      do n = mf, jtrun
         i_array(n, m) = i
         j_array(n, m) = j
         i = i + 1
         if (i .gt. jtlen_all(j)) then
            i = 1
            j = j + 1
         end if
      end do
   end do
   !$acc enter data copyin(i_array, j_array) async(async_id)

   !$acc wait(async_id)

   nlen = jtp*2*levp
   !$acc host_data use_device(b1, b2)
   call MPI_ALLTOALL(b1, nlen, MPI_RTYPE, &
                     b2, nlen, MPI_RTYPE, &
                     row_comm, IERR)
   !$acc end host_data

   !$acc parallel loop collapse(4) private(mf, i, j) async(async_id)
   do m = 1, mlistnum
      do n = 1, jtrun
         do q = 1, 2
            do k = 1, levp
               mf = mlist(m)
               if (n .ge. mf) then
                  i = i_array(n, m)
                  j = j_array(n, m)
                  aout(k, q, n, m) = b2(i, q, k, j)
               end if
            end do
         end do
      end do
   end do

   !$acc exit data delete(b1, b2, i_array, j_array) async(async_id)

   return
end
