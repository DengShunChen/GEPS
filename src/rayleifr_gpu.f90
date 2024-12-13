#define NCCLCHECK(ierr) call nccl_check_helper(ierr, __FILE__, __LINE__)

subroutine rayleifr_gpu(nx, my, my_max, lev, rad, cosl, dt, ut, vt)
   use mpe
   use rank
   use index
   use const, only: RTYPE
   use openacc
   use cudafor

   implicit none
   integer nx, my, my_max, lev
   real rad

   real(kind=RTYPE) cosl(my)
   real(kind=RTYPE) ut(nxp, lev, my_max), vt(nxp, lev, my_max)
   real wmax(lev), wmaxtmp, wmax_buf(lev), wt
   !
   integer, parameter :: levtop = 6, lev2 = 3

   real ckdx, day, windmax
   !     data ckdx/20./, day/1./
   data ckdx/40./, day/1./

   data windmax/140./
   !     data windmax/130./

   logical dofric

   integer jj, j, nxj, k, i
   real xx, dt, ckdy, frictime, fac, fricd, cdx, cdy
   integer async_id, istat, ierr
   integer(kind=cuda_stream_kind) :: stream
   async_id = 1
   stream = acc_get_cuda_stream(async_id)
   !
   dofric = .false.
   ckdy = ckdx
   frictime = 1./(day*86400.)
   wmaxtmp = windmax
   !$acc enter data create(wmax, wmax_buf) async(async_id)

   !$acc parallel loop gang private(wt) async(async_id)
   do k = 1, lev
      wt = 0.
      !$acc loop vector collapse(2) reduction(max:wt) private(j, nxj, xx)
      do jj = 1, jlistnum
         do i = 1, nxp
            j = jlist1(jj)
            nxj = nxdef_2d(j)
            if (i .le. nxj) then
               xx = rad/cosl(j)
               wt = max(wt, xx*sqrt(ut(i, k, jj)**2 + vt(i, k, jj)**2))
            end if
         end do
      end do
      wmax_buf(k) = wt
   end do
   !$acc host_data use_device(wmax, wmax_buf)
   NCCLCHECK(ncclAllReduce(wmax_buf, wmax, levtop, ncclFloat64, ncclMax, nccl_comm_gfs, stream))
   !$acc end host_data
   !$acc update self(wmax) async(async_id)
   !$acc wait(async_id)
   do k = 1, levtop
      wmaxtmp = max(wmax(k), wmaxtmp)
      if (wmax(k) .gt. windmax) then
         if (myrank .eq. 0) print *, ' (GPU)rayleifr k=', k, ' wmax=', wmax(k)
      end if
   end do

   if (wmaxtmp .gt. windmax) dofric = .true.
   if (myrank .eq. 0) print *, ' (GPU)rayleifr  dofric=', dofric
   if (dofric) then
      !$acc parallel loop collapse(3) private(j, nxj, fac, fricd, cdx, cdy) &
      !$acc& async(async_id)
      do k = 1, levtop
         do jj = 1, jlistnum
            do i = 1, nxp
               j = jlist1(jj)
               nxj = nxdef_2d(j)
               if (i .le. nxj) then
                  fac = 1.+exp(-1.0*(k - lev2)/lev2)
                  if (fac .gt. 3.) fac = 3.
                  fricd = frictime*fac
                  cdx = ckdx*fricd
                  cdy = ckdy*fricd
                  ut(i, k, jj) = ut(i, k, jj) - cdx*ut(i, k, jj)
                  vt(i, k, jj) = vt(i, k, jj) - cdy*vt(i, k, jj)
               end if
            end do
         end do
      end do

   end if
   !$acc exit data delete(wmax, wmax_buf) async(async_id)
   !
   return
end subroutine rayleifr_gpu
