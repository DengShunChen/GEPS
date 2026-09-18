! NDSL WE<->NS. nsizey==1: alltoall is local, skip works/workr (~14GB at TCo383).
! nsizey>1 (2026-09-16): the *_multi routines below, ported from src/nvidia.

subroutine para_we2ns_gpu(a, b, levs, latg)
   use const, only: RTYPE
   use grid, only: lonfull, lonhalf, lonpart, latfull, latpart, &
                   mylatlen, lonlen, lonstr
   use index, only: nsizey, jlist1_sl
   implicit none
   integer levs, latg
   real(kind=RTYPE) a(lonfull, levs, latpart)
   real(kind=RTYPE) b(latfull, levs, lonpart)
   integer i, j, k, lat1, lat2
   ! #region agent log
   integer(kind=8) :: dbg0, dbg1, dbg2
   interface
      subroutine geps_dbg_vram(hyp, locid, p0, p1, p2)
         integer hyp, locid
         integer(kind=8) p0, p1, p2
      end subroutine
   end interface
   ! #endregion

   ! #region agent log
   dbg0 = int(levs, 8); dbg1 = int(nsizey, 8); dbg2 = int(mylatlen, 8)
   call geps_dbg_vram(2, 37, dbg0, dbg1, dbg2)
   ! #endregion

   if (nsizey > 1) then
      ! 2026-09-16: multi-rank path, a direct OpenMP port of the NVIDIA original
      ! (src/nvidia/ndslfv_para_gpu.f90): pack -> RCCL alltoallv -> unpack.
      ! works/workr are transient (~4 GB each at NPEY=2), so only the 1-rank
      ! path skips them. The wait after the collective is the layer-6 rule.
      call para_we2ns_gpu_multi(a, b, levs, latg)
      return
   end if

   !$omp target teams distribute parallel do collapse(3) private(lat1, lat2)
   do j = 1, mylatlen
      do i = 1, lonlen(1)
         do k = 1, levs
            lat1 = jlist1_sl(j, 1)
            lat2 = latfull + 1 - lat1
            b(lat1, k, i) = a(i + lonstr(1) - 1, k, j)
            b(lat2, k, i) = a(i + lonstr(1) - 1 + lonhalf, k, j)
         end do
      end do
   end do
end subroutine para_we2ns_gpu

subroutine para_ns2we_gpu(a, b, levs, latg)
   use const, only: RTYPE
   use grid, only: lonfull, lonhalf, lonpart, latfull, latpart, &
                   mylatlen, lonlen, lonstr
   use index, only: nsizey, jlist1_sl
   implicit none
   integer levs, latg
   real(kind=RTYPE) a(latfull, levs, lonpart)
   real(kind=RTYPE) b(lonfull, levs, latpart)
   integer i, j, k, lat1, lat2

   if (nsizey > 1) then
      call para_ns2we_gpu_multi(a, b, levs, latg)
      return
   end if

   !$omp target teams distribute parallel do collapse(3) private(lat1, lat2)
   do j = 1, mylatlen
      do i = 1, lonlen(1)
         do k = 1, levs
            lat1 = jlist1_sl(j, 1)
            lat2 = latfull + 1 - lat1
            b(i + lonstr(1) - 1, k, j) = a(lat1, k, i)
            b(i + lonstr(1) - 1 + lonhalf, k, j) = a(lat2, k, i)
         end do
      end do
   end do
end subroutine para_ns2we_gpu

! ======================================================================
! Multi-rank versions (ported from src/nvidia/ndslfv_para_gpu.f90, 2026-09-16).
subroutine para_we2ns_gpu_multi(a, b, levs, latg)
   use const, only: RTYPE
   use grid, only: lonfull, lonhalf, lonpart, lonlenmax, mylonlen, &
                   latfull, lathalf, latpart, latlenmax, mylatlen, &
                   latstr, latlen, lonstr, lonlen
   use index, only: nsizey, jlist1_sl, nccl_col_comm
   implicit none
   integer levs, latg
   real(kind=RTYPE) a(lonfull, levs, latpart)
   real(kind=RTYPE) b(latfull, levs, lonpart)
   real(kind=RTYPE) works(2, levs, lonlenmax*latlenmax, nsizey)
   real(kind=RTYPE) workr(2, levs, lonlenmax*latlenmax, nsizey)
   integer lensend(nsizey), lenrecv(nsizey)
   integer i, j, k, n, mn, lat1, lat2
   integer, parameter :: async_id = 1

   !$omp target enter data map(alloc:works, workr)
   !$omp target teams distribute parallel do collapse(4) private(mn)
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
      lensend(n) = mylatlen*lonlen(n)*2*levs
      lenrecv(n) = latlen(n)*mylonlen*2*levs
   end do

   call nccl_alltoallv_stride(works, lensend, lonlenmax*latlenmax*2*levs, &
                              workr, lenrecv, lonlenmax*latlenmax*2*levs, nccl_col_comm, nsizey, async_id)
   call geps_acc_wait(async_id)

   !$omp target teams distribute parallel do collapse(4) private(lat1, lat2, mn)
   do n = 1, nsizey
      do i = 1, mylonlen
         do k = 1, levs
            do j = 1, latlenmax
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
   !$omp target exit data map(release:works, workr)
end subroutine para_we2ns_gpu_multi

subroutine para_ns2we_gpu_multi(a, b, levs, latg)
   use const, only: RTYPE
   use grid, only: lonfull, lonhalf, lonpart, lonlenmax, mylonlen, &
                   latfull, lathalf, latpart, latlenmax, mylatlen, &
                   latstr, latlen, lonstr, lonlen
   use index, only: nsizey, jlist1_sl, nccl_col_comm
   implicit none
   integer levs, latg
   real(kind=RTYPE) a(latfull, levs, lonpart)
   real(kind=RTYPE) b(lonfull, levs, latpart)
   real(kind=RTYPE) works(2, levs, lonlenmax*latlenmax, nsizey)
   real(kind=RTYPE) workr(2, levs, lonlenmax*latlenmax, nsizey)
   integer lensend(nsizey), lenrecv(nsizey)
   integer i, j, k, n, mn, lat1, lat2
   integer, parameter :: async_id = 1

   !$omp target enter data map(alloc:works, workr)
   !$omp target teams distribute parallel do collapse(4) private(lat1, lat2, mn)
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
      lensend(n) = latlen(n)*mylonlen*2*levs
      lenrecv(n) = mylatlen*lonlen(n)*2*levs
   end do

   call nccl_alltoallv_stride(works, lensend, lonlenmax*latlenmax*2*levs, &
                              workr, lenrecv, lonlenmax*latlenmax*2*levs, nccl_col_comm, nsizey, async_id)
   call geps_acc_wait(async_id)

   !$omp target teams distribute parallel do collapse(4) private(mn)
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
   !$omp target exit data map(release:works, workr)
end subroutine para_ns2we_gpu_multi
