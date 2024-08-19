subroutine para_we2ns_gpu(a, b, levs, latg)
   ! Present on device: a, b, lonlen, lonstr, latlen, jlist1_sl
!
! mpi transport from full dimension of west-east to full dimension of
! north-south with latitude shuffl. created by hann-ming henry juang
!
! program log
! 2011 02 20 : henry juang, created for ndsl advection
!
!
   use const, only: RTYPE, MPI_RTYPE
   use grid, only: lonfull, lonhalf, lonpart, lonlenmax, mylonlen, &
                   latfull, lathalf, latpart, latlenmax, mylatlen, &
                   latstr, latlen, lonstr, lonlen
   use rank
   use index
   use mpi
   use cudafor
   use openacc

   implicit none

   integer levs, latg
   real(kind=RTYPE) a(lonfull, levs, latpart)
   real(kind=RTYPE) b(latfull, levs, lonpart)
   real(kind=RTYPE) works(2, levs, lonlenmax*latlenmax, nsizey)
   real(kind=RTYPE) workr(2, levs, lonlenmax*latlenmax, nsizey)
   integer lensend(nsizey), lenrecv(nsizey)
   integer locsend(nsizey), locrecv(nsizey)
   integer i, j, k, n, mn, lat1, lat2, ierr
   integer async_id
   integer(kind=cuda_stream_kind) stream

   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   !$acc enter data create(works, workr) async(async_id)
   !$acc parallel loop collapse(4) private(mn) async(async_id)
   do n = 1, nsizey
      do j = 1, mylatlen
         do i = 1, lonlenmax
            do k = 1, levs
               if (i .le. lonlen(n)) then
                  mn = (j - 1)*lonlen(n) + i
                  works(1, k, mn, n) = a(i + lonstr(n) - 1, k, j)
                  works(2, k, mn, n) = a(i + lonstr(n) - 1 + lonhalf, k, j)
               end if
            end do
         end do
      end do
   end do

   do n = 1, nsizey
      mn = mylatlen*lonlen(n)
      lensend(n) = mn*2*levs
      locsend(n) = (n - 1)*lonlenmax*latlenmax*2*levs
      lenrecv(n) = latlen(n)*mylonlen*2*levs
      locrecv(n) = locsend(n)
   end do

   call nccl_alltoallv_stride_fp64(works, lensend, lonlenmax*latlenmax*2*levs, workr, lenrecv, lonlenmax*latlenmax*2*levs, nccl_col_comm, nsizey, async_id)

   !$acc parallel loop collapse(4) private(lat1, lat2, mn) async(async_id)
   do n = 1, nsizey
      do j = 1, latlenmax
         do i = 1, mylonlen
            do k = 1, levs
               if (j .le. latlen(n)) then
                  lat1 = jlist1_sl(j, n)
                  lat2 = latfull + 1 - lat1
                  mn = (j - 1)*mylonlen + i
                  b(lat1, k, i) = workr(1, k, mn, n)
                  b(lat2, k, i) = workr(2, k, mn, n)
               end if
            end do
         end do
      end do
   end do

   !$acc exit data delete(works, workr) async(async_id)

end

! ======================================================================
subroutine para_ns2we_gpu(a, b, levs, latg)
   ! Present on device: latlen, jlist1_sl, a, lonlen, lonstr, b
!
! mpi transport from full dimension of west-east to full dimension of
! north-south with latitude shuffl.
!
   use const, only: RTYPE, MPI_RTYPE
   use grid, only: lonfull, lonhalf, lonpart, lonlenmax, mylonlen, &
                   latfull, lathalf, latpart, latlenmax, mylatlen, &
                   latstr, latlen, lonstr, lonlen
   use rank
   use index
   use mpi

   implicit none

   integer levs, latg
   real(kind=RTYPE) a(latfull, levs, lonpart)
   real(kind=RTYPE) b(lonfull, levs, latpart)
   real(kind=RTYPE) works(2, levs, lonlenmax*latlenmax, nsizey)
   real(kind=RTYPE) workr(2, levs, lonlenmax*latlenmax, nsizey)
   integer lensend(nsizey), lenrecv(nsizey)
   integer locsend(nsizey), locrecv(nsizey)
   integer i, j, k, n, mn, lat1, lat2, ierr
   integer :: async_id

   async_id = 1

   !$acc enter data create(works, workr) async(async_id)
   !$acc parallel loop collapse(4) private(lat1, lat2, mn) async(async_id)
   do n = 1, nsizey
      do j = 1, latlenmax
         do i = 1, mylonlen
            do k = 1, levs
               if (j .le. latlen(n)) then
                  lat1 = jlist1_sl(j, n)
                  lat2 = latfull + 1 - lat1
                  mn = (j - 1)*mylonlen + i
                  works(1, k, mn, n) = a(lat1, k, i)
                  works(2, k, mn, n) = a(lat2, k, i)
               end if
            end do
         end do
      end do
   end do

   do n = 1, nsizey
      mn = latlen(n)*mylonlen
      lensend(n) = mn*2*levs
      locsend(n) = (n - 1)*lonlenmax*latlenmax*2*levs
      lenrecv(n) = mylatlen*lonlen(n)*2*levs
      locrecv(n) = locsend(n)
   end do

   call nccl_alltoallv_stride_fp64(works, lensend, lonlenmax*latlenmax*2*levs, workr, lenrecv, lonlenmax*latlenmax*2*levs, nccl_col_comm, nsizey, async_id)

   !$acc parallel loop collapse(4) private(mn) async(async_id)
   do n = 1, nsizey
      do j = 1, mylatlen
         do i = 1, lonlenmax
            do k = 1, levs
               if (i .le. lonlen(n)) then
                  mn = (j - 1)*lonlen(n) + i
                  b(i + lonstr(n) - 1, k, j) = workr(1, k, mn, n)
                  b(i + lonstr(n) - 1 + lonhalf, k, j) = workr(2, k, mn, n)
               end if
            end do
         end do
      end do
   end do
   !$acc exit data delete(works, workr) async(async_id)
end
