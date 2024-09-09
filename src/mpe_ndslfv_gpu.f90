#define NCCLCHECK(ierr) call nccl_check_helper(ierr, __FILE__, __LINE__)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mpe_nxp2full(aout, ainp, nv)
   ! from (nx partial, lev full, var. packed, my partial)
   ! to   (nx full, my full, lev full, var. packed)
   use const, only: RTYPE, MPI_RTYPE
   use param, only: nx, my, lev, my_max
   use index, only: nxp, nsizex, jlist2_2d, nxjlen_all
   use rank, only: nsize, mpi_comm_gfs
   use cudafor
   use mpi
   implicit none
   integer, intent(in):: nv
   real(kind=RTYPE) aout(nx, my, lev, nv), &
      ainp(nxp, lev, nv, my_max), &
      work(nxp*lev*nv, my_max*nsize)
   integer n, k, ierr, j, ii, i, jf, nn, kk

   call MPI_ALLGATHER(ainp, nxp*lev*nv*my_max, MPI_RTYPE, &
                      work, nxp*lev*nv*my_max, MPI_RTYPE, &
                      MPI_COMM_gfs, IERR)
   do n = 1, nv
      do k = 1, lev
         do j = 1, my
            ii = 1
            do i = 1, nsizex
               jf = jlist2_2d(i, j)
               nn = nxjlen_all(i, j)
               kk = 1 + (k - 1)*nxp + (n - 1)*nxp*lev
               aout(ii:ii + nn - 1, j, k, n) = work(kk:kk + nn - 1, jf)
               ii = ii + nn
            end do
         end do
      end do
   end do

   return
end subroutine mpe_nxp2full
! ============================================================
subroutine mpe_full2nxp(aout, ainp, nv)
   ! from (nx full, my full, lev full, var. packed)
   ! to (nx partial, lev full, var. packed, my partial)
   use const, only: RTYPE, MPI_RTYPE
   use param, only: nx, my, lev, my_max
   use index, only: nxp, nsizex, jlist2_2d, nxjlen_all
   use rank, only: nsize, mpi_comm_gfs, myrank
   use cudafor
   use mpi
   implicit none
   integer, intent(in):: nv
   real(kind=RTYPE) ainp(nx, my, lev, nv), &
      aout(nxp*lev*nv, my_max), &
      work(nxp*lev*nv, my_max*nsize)
   integer n, k, ierr, j, ii, i, jf, nn, kk

   if (myrank .eq. 0) then
      do n = 1, nv
         do k = 1, lev
            do j = 1, my
               ii = 1
               do i = 1, nsizex
                  jf = jlist2_2d(i, j)
                  nn = nxjlen_all(i, j)
                  kk = 1 + (k - 1)*nxp + (n - 1)*nxp*lev
                  work(kk:kk + nn - 1, jf) = ainp(ii:ii + nn - 1, j, k, n)
                  ii = ii + nn
               end do
            end do
         end do
      end do
   end if
   call MPI_SCATTER(work, nxp*lev*nv*my_max, MPI_RTYPE, &
                    aout, nxp*lev*nv*my_max, MPI_RTYPE, &
                    0, MPI_COMM_gfs, IERR)

   return
end subroutine mpe_full2nxp
! ============================================================
subroutine advh_reshape_c2g(aout, ainp, ainp_d)
   use param
   use const, only: RTYPE, cosl
   use index, only: nsizex, jlist2_2d, nxjlen_all
   use grid, only: nxp
   use rank, only: myrank, nsize
   use cudafor
   implicit none
   real(kind=RTYPE), intent(out) :: aout(nx, my, lev)
   real(kind=RTYPE), intent(in)  :: ainp(nxp*lev, my_max*nsize)
   real(kind=RTYPE), device :: ainp_d(nxp*lev, my_max*nsize)

   integer i, j, k, m, ii, jf, kk, nn, ierr

   ierr = cudamemcpy(ainp_d, ainp, &
                     nxp*lev*my_max*nsize, &
                     cudaMemcpyHostToDevice)
   !$acc parallel loop deviceptr(ainp_d) &
   !$acc& copyin(cosl(:my)) &
   !$acc& present( &
   !$acc& my,lev,nxp,nsizex, &
   !$acc& jlist2_2d(nsizex,my), nxjlen_all(nsizex,my), &
   !$acc& aout(:nx,:my,:lev)) &
   !$acc& collapse(2) &
   !$acc& private(ii,jf,nn,kk,i,m)
   do k = 1, lev
      do j = 1, my
         ii = 0
         !$acc loop seq
         do i = 1, nsizex
            jf = jlist2_2d(i, j)
            kk = (k - 1)*nxp
            nn = nxjlen_all(i, j)
            !$acc loop independent
            do m = 1, nn
               aout(ii + m, j, k) = ainp_d(kk + m, jf)
            end do
            ii = ii + nn
         end do
      end do
   end do

end subroutine advh_reshape_c2g
! ============================================================
subroutine advh_reshape_c2g_q(aout, ainp_d, nv, async)
   use param
   use const, only: RTYPE, cosl
   use index, only: nsizex, jlist2_2d, nxjlen_all
   use grid, only: nxp
   use rank, only: myrank, nsize
   use cudafor
   implicit none
   real(kind=RTYPE), intent(out) :: aout(nx, my, lev, nv)
   real(kind=RTYPE), device      :: ainp_d(nxp*lev, my_max*nsize, nv)
   integer, intent(in) :: nv, async

   integer n, i, j, k, m, ii, jf, kk, nn, ierr
   integer jf_i, jf_j, jf_t

   !$acc parallel loop deviceptr(ainp_d) async(async) &
   !$acc& copyin(cosl(:my)) &
   !$acc& present( &
   !$acc& my,lev,nxp,nsizex, &
   !$acc& jlist2_2d(nsizex,my), nxjlen_all(nsizex,my), &
   !$acc& aout(:nx,:my,:lev,:nv)) &
   !$acc& collapse(3) &
   !$acc& private(ii,jf,nn,kk,i,m,jf_i,jf_j,jf_t)
   do n = 1, nv
      do k = 1, lev
         do j = 1, my
            ii = 0
            !$acc loop seq
            do i = 1, nsizex
               jf = jlist2_2d(i, j)
               jf_j = int((jf - 1)/my_max) + 1
               jf_i = jf - (jf_j - 1)*my_max
               jf_t = jf_j + (jf_i - 1)*nsize
               kk = (k - 1)*nxp
               nn = nxjlen_all(i, j)
               !$acc loop independent
               do m = 1, nn
                  aout(ii + m, j, k, n) = ainp_d(kk + m, jf_t, n)
               end do
               ii = ii + nn
            end do
         end do
      end do
   end do
end subroutine advh_reshape_c2g_q
! ============================================================
subroutine advh_reshape_g2c(aout, aout_d, ainp)
   use param
   use const, only: RTYPE, cosl
   use index, only: nsizex, jlist2_2d, nxjlen_all
   use grid, only: nxp
   use rank, only: myrank, nsize
   use cudafor
   implicit none
   real(kind=RTYPE), intent(out) :: aout(nxp*lev, my_max*nsize)
   real(kind=RTYPE), device      :: aout_d(nxp*lev, my_max*nsize)
   real(kind=RTYPE), intent(in)  :: ainp(nx, my, lev)

   integer i, j, k, m, ii, jf, kk, nn, ierr
   !$acc parallel loop deviceptr(aout_d) &
   !$acc& copyin(cosl(:my)) &
   !$acc& present( &
   !$acc& my,lev,nxp,nsizex, &
   !$acc& jlist2_2d(nsizex,my), nxjlen_all(nsizex,my), &
   !$acc& ainp(:nx,:my,:lev)) &
   !$acc& collapse(2) &
   !$acc& private(ii,jf,nn,kk,i,m)
   do k = 1, lev
      do j = 1, my
         ii = 0
         !$acc loop seq
         do i = 1, nsizex
            jf = jlist2_2d(i, j)
            nn = nxjlen_all(i, j)
            kk = (k - 1)*nxp
            !$acc loop independent
            do m = 1, nn
               aout_d(kk + m, jf) = ainp(ii + m, j, k)
            end do
            ii = ii + nn
         end do
      end do
   end do
   ierr = cudamemcpy(aout, aout_d, &
                     nxp*lev*my_max*nsize, &
                     cudaMemcpyDeviceToHost)

end subroutine advh_reshape_g2c
! ============================================================
subroutine advh_reshape_g2c_q(aout_d, ainp, nv, async)
   use param
   use const, only: RTYPE, cosl
   use index, only: nsizex, jlist2_2d, nxjlen_all
   use grid, only: nxp
   use rank, only: myrank, nsize
   use cudafor
   implicit none
   real(kind=RTYPE), device      :: aout_d(nxp*lev, my_max*nsize, nv)
   real(kind=RTYPE), intent(in)  :: ainp(nx, my, lev, nv)
   integer, intent(in) :: nv, async

   integer n, i, j, k, m, ii, jf, kk, nn, ierr
   integer jf_i, jf_j, jf_t

   !$acc parallel loop deviceptr(aout_d) async(async) &
   !$acc& copyin(cosl(:my)) &
   !$acc& present( &
   !$acc& my,lev,nxp,nsizex, &
   !$acc& jlist2_2d(nsizex,my), nxjlen_all(nsizex,my), &
   !$acc& ainp(:nx,:my,:lev,:nv)) &
   !$acc& collapse(3) &
   !$acc& private(ii,jf,nn,kk,i,m,jf_i,jf_j,jf_t)
   do n = 1, nv
      do k = 1, lev
         do j = 1, my
            ii = 0
            !$acc loop seq
            do i = 1, nsizex
               jf = jlist2_2d(i, j)
               jf_j = int((jf - 1)/my_max) + 1
               jf_i = jf - (jf_j - 1)*my_max
               jf_t = jf_j + (jf_i - 1)*nsize
               nn = nxjlen_all(i, j)
               kk = (k - 1)*nxp
               !$acc loop independent
               do m = 1, nn
                  aout_d(kk + m, jf_t, n) = ainp(ii + m, j, k, n)
               end do
               ii = ii + nn
            end do
         end do
      end do
   end do

end subroutine advh_reshape_g2c_q
! ============================================================
subroutine advh_gather4GPU(aout, ainp, nv, id)
   use const, only: RTYPE, MPI_RTYPE
   use param, only: nx, my, lev, my_max
   use index, only: nxp, nsizex, jlist2_2d, nxjlen_all
   use rank, only: nsize, mpi_comm_gfs
   use cudafor
   use mpi
   implicit none
   integer, intent(in):: nv, id
   real(kind=RTYPE) aout(nxp*lev*nv, my_max*nsize)
   real(kind=RTYPE) ainp(nxp, lev, nv, my_max)
   integer ierr

   call MPI_GATHER(ainp, nxp*lev*nv*my_max, MPI_RTYPE, &
                   aout, nxp*lev*nv*my_max, MPI_RTYPE, &
                   id, MPI_COMM_gfs, IERR)
   return
end subroutine advh_gather4GPU
! ============================================================
subroutine advh_allgather4GPU(aout, ainp, nv)
   use const, only: RTYPE, MPI_RTYPE
   use param, only: nx, my, lev, my_max
   use index, only: nxp, nsizex, jlist2_2d, nxjlen_all
   use rank, only: nsize, mpi_comm_gfs
   use cudafor
   use mpi
   implicit none
   integer, intent(in):: nv
   real(kind=RTYPE) aout(nxp*lev*nv, my_max*nsize)
   real(kind=RTYPE) ainp(nxp, lev, nv, my_max)
   integer ierr

   call MPI_ALLGATHER(ainp, nxp*lev*nv*my_max, MPI_RTYPE, &
                      aout, nxp*lev*nv*my_max, MPI_RTYPE, &
                      MPI_COMM_gfs, IERR)
   return
end subroutine advh_allgather4GPU
! ============================================================
subroutine advh_gather4GPU_dev(aout, ainp, nv, id)
   use const, only: RTYPE, MPI_RTYPE
   use param, only: nx, my, lev, my_max
   use index, only: nxp
   use rank, only: nsize, mpi_comm_gfs
   use cudafor
   use mpi
   implicit none
   integer, intent(in):: nv, id
   real(kind=RTYPE), device:: aout(nxp*lev*nv, my_max*nsize)
   real(kind=RTYPE), device:: ainp(nxp*lev*nv, my_max)
   integer pts, ierr

   pts = nxp*lev*nv*my_max
   call MPI_GATHER(ainp, pts, MPI_RTYPE, &
                   aout, pts, MPI_RTYPE, &
                   id, MPI_COMM_gfs, IERR)

   return
end subroutine advh_gather4GPU_dev
! ============================================================
subroutine advh_gather4GPU_nccl(aout, ainp, nv, id)
  use const, only: RTYPE, MPI_RTYPE
  use param, only: nx, my, lev, my_max
  use index, only: nxp
  use rank, only: nsize, nccl_comm_gfs, myrank
  use cudafor
  use openacc
  use nccl
  implicit none
  integer, intent(in):: nv, id
  real(kind=RTYPE), device:: aout(nxp*lev*nv*my_max,nsize)
  real(kind=RTYPE), device:: ainp(nxp*lev*nv*my_max)
  integer pts, i
  integer(kind=cuda_stream_kind) :: stream
  integer async_id
  pts = nxp*lev*nv*my_max

  async_id = 1
  stream = acc_get_cuda_stream(async_id)
  NCCLCHECK(ncclGroupStart())

  NCCLCHECK(ncclSend(ainp, pts, ncclFloat64, id, nccl_comm_gfs, stream))

  if (myrank.eq.0) then
     do i = 1, nsize
        NCCLCHECK(ncclRecv(aout(1, i), pts, ncclFloat64, i-1, nccl_comm_gfs, stream))
     end do
  end if

  NCCLCHECK(ncclGroupEnd())
  !$acc wait(async_id)

  return
end subroutine advh_gather4GPU_nccl
! ============================================================
subroutine advh_gather4GPU_q(aout, ainp, nv, ncld, myrank)
   use const, only: RTYPE, MPI_RTYPE
   use param, only: nx, my, lev, my_max
   use index, only: nxp, nsizex, jlist2_2d, nxjlen_all
   use rank, only: nsize, mpi_comm_gfs
   use cudafor
   use mpi
   implicit none
   integer, intent(in):: nv, ncld, myrank
   real(kind=RTYPE) aout(nxp*lev, nsize, my_max)
   real(kind=RTYPE) ainp(nxp*lev, ncld, my_max)
   integer ierr
   integer j

   do j = 1, my_max
      call MPI_GATHER(ainp(1, nv, j), nxp*lev, MPI_RTYPE, &
                      aout(1, 1, j), nxp*lev, MPI_RTYPE, &
                      myrank, MPI_COMM_gfs, IERR)
   end do
   return
end subroutine advh_gather4GPU_q
! ============================================================
subroutine advh_scatter4GPU(aout, ainp, nv, id)
   use const, only: RTYPE, MPI_RTYPE
   use param, only: nx, my, lev, my_max
   use index, only: nxp
   use rank, only: nsize, mpi_comm_gfs
   use cudafor
   use mpi
   implicit none
   integer, intent(in):: nv, id
   real(kind=RTYPE) aout(nxp*lev*nv, my_max)
   real(kind=RTYPE) ainp(nxp*lev*nv, my_max*nsize)
   integer ierr

   call MPI_SCATTER(ainp, nxp*lev*nv*my_max, MPI_RTYPE, &
                    aout, nxp*lev*nv*my_max, MPI_RTYPE, &
                    id, MPI_COMM_gfs, IERR)

   return
end subroutine advh_scatter4GPU
! ============================================================
subroutine advh_scatter4GPU_dev(aout, ainp, nv, id)
   use const, only: RTYPE, MPI_RTYPE
   use param, only: nx, my, lev, my_max
   use index, only: nxp
   use rank, only: nsize, mpi_comm_gfs
   use cudafor
   use mpi
   implicit none
   integer, intent(in):: nv, id
   real(kind=RTYPE), device:: aout(nxp*lev*nv, my_max)
   real(kind=RTYPE), device:: ainp(nxp*lev*nv, my_max*nsize)
   integer ierr, pts

   pts = nxp*lev*nv*my_max
   call MPI_SCATTER(ainp, pts, MPI_RTYPE, &
                    aout, pts, MPI_RTYPE, &
                    id, MPI_COMM_gfs, IERR)

   return
 end subroutine advh_scatter4GPU_dev
 ! ============================================================
 subroutine advh_scatter4GPU_nccl(aout, ainp, nv, id)
   use const, only: RTYPE, MPI_RTYPE
   use param, only: nx, my, lev, my_max
   use index, only: nxp
   use rank, only: nsize, nccl_comm_gfs, myrank
   use cudafor
   use openacc
   use nccl
   implicit none
   integer, intent(in):: nv, id
   real(kind=RTYPE), device:: aout(nxp*lev*nv*my_max)
   real(kind=RTYPE), device:: ainp(nxp*lev*nv*my_max,nsize)
   integer i, pts
   integer(kind=cuda_stream_kind) :: stream
   integer async_id

   pts = nxp*lev*nv*my_max
   async_id = 1
   stream = acc_get_cuda_stream(async_id)
   NCCLCHECK(ncclGroupStart())

   if (myrank.eq.0) then
      do i = 1, nsize
         NCCLCHECK(ncclSend(ainp(1,i), pts, ncclFloat64, i-1, nccl_comm_gfs, stream))

      end do
   end if
   NCCLCHECK(ncclRecv(aout, pts, ncclFloat64, id, nccl_comm_gfs, stream))

   NCCLCHECK(ncclGroupEnd())
   !$acc wait(async_id)
   return
 end subroutine advh_scatter4GPU_nccl
 ! ============================================================
 subroutine advh_scatter4GPU_q(aout, ainp, nv, ncld, myrank)
   use const, only: RTYPE, MPI_RTYPE
   use param, only: nx, my, lev, my_max
   use index, only: nxp
   use rank, only: nsize, mpi_comm_gfs
   use cudafor
   use mpi
   implicit none
   integer, intent(in):: nv, ncld, myrank
   real(kind=RTYPE) aout(nxp*lev, ncld, my_max)
   real(kind=RTYPE) ainp(nxp*lev, nsize, my_max)
   integer ierr, j

   do j = 1, my_max
      call MPI_SCATTER(ainp(1, 1, j), nxp*lev, MPI_RTYPE, &
                       aout(1, nv, j), nxp*lev, MPI_RTYPE, &
                       myrank, MPI_COMM_gfs, IERR)
   end do

   return
end subroutine advh_scatter4GPU_q
! ============================================================
