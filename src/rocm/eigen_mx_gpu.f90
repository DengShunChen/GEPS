! First-class ROCm eigen_mx: OpenMP target data + hipSOLVER Dsyevj.
! Loop structure matches nvidia/eigen_mx_gpu.f90; no CUDA Fortran host_data.

subroutine eigen_mx_gpu(mx, eval, nnlist, no)
   use iso_c_binding
   use const, only: nnmivm
   use param, only: jtrun, jtmax
   use index, only: mlistnum
   implicit none

   real(kind=8), intent(inout), target :: mx(no*no, 2*jtmax*nnmivm)
   real(kind=8), intent(out), target   :: eval(no, 2*jtmax*nnmivm)
   integer, intent(in)                 :: nnlist(2*jtmax*nnmivm)
   integer, intent(in)                 :: no

   integer, target :: devinfo(2*jtmax*nnmivm)
   ! #region agent log
   integer :: dbg_c, dbg_r, dbg_n, dbg_worstnn
   real(kind=8) :: dbg_s, dbg_v, dbg_worstnorm, dbg_worstmax
   integer(kind=8) :: dbg0, dbg1, dbg2
   interface
      subroutine geps_dbg_vram(hyp, locid, p0, p1, p2)
         integer hyp, locid
         integer(kind=8) p0, p1, p2
      end subroutine
   end interface
   ! #endregion
   integer :: i, j, m, L, nn, ierr

   interface
      function geps_rocm_dsyevj_dev(A, n, W, info) bind(C, name="geps_rocm_dsyevj_dev")
         import c_ptr, c_int
         type(c_ptr), value :: A, W, info
         integer(c_int), value :: n
         integer(c_int) :: geps_rocm_dsyevj_dev
      end function
   end interface

   devinfo = 0
   !$omp target enter data map(alloc: devinfo)

   !$omp target data use_device_addr(mx, eval, devinfo)
   do L = 1, nnmivm
      do m = 1, mlistnum
         do i = 1, 2
            j = i + (m - 1)*2 + (L - 1)*jtmax*2
            nn = nnlist(j)
            if (nn <= 0) cycle
            ierr = geps_rocm_dsyevj_dev(c_loc(mx(1, j)), nn, c_loc(eval(1, j)), &
                                        c_loc(devinfo(j)))
            if (ierr /= 0) then
               print *, 'geps_rocm_dsyevj_dev failed', ierr, 'n=', nn, 'j=', j
               error stop
            end if
         end do
      end do
   end do
   !$omp end target data

   !$omp target exit data map(delete: devinfo)

   ! #region agent log: verify what the solve actually produced.
   ! Only the first nn*nn entries of each column are used - the rest of the
   ! column is padding that is never written, so maxval() over the whole array
   ! is meaningless (that mistake cost us a run; see handoff layer 8).
   ! Eigenvectors of a symmetric matrix must be orthonormal: every column of
   ! the nn x nn result must have sum of squares == 1 and |entry| <= 1.
   !$omp target update from(mx, eval)
   dbg_n = 0
   dbg_worstnorm = 0.0d0
   dbg_worstmax = 0.0d0
   dbg_worstnn = 0
   do L = 1, nnmivm
      do m = 1, mlistnum
         do i = 1, 2
            j = i + (m - 1)*2 + (L - 1)*jtmax*2
            nn = nnlist(j)
            if (nn <= 0) cycle
            ! NO sampling cap. The previous version stopped after 64
            ! matrices, and the loop order (L, m, k) visits the SMALL ones
            ! first (ns = 48, 50, 144, ...), so the large ones - ns up to 530,
            ! which is where the mass of `bal` lives - were never checked.
            ! Third time this sampling mistake bit us today (handoff layer 9).
            dbg_n = dbg_n + 1
            do dbg_c = 1, nn
               dbg_s = 0.0d0
               do dbg_r = 1, nn
                  dbg_v = mx((dbg_c - 1)*nn + dbg_r, j)
                  dbg_s = dbg_s + dbg_v*dbg_v
                  if (abs(dbg_v) > dbg_worstmax) dbg_worstmax = abs(dbg_v)
               end do
               if (abs(dbg_s - 1.0d0) > dbg_worstnorm) then
                  dbg_worstnorm = abs(dbg_s - 1.0d0)
                  dbg_worstnn = nn
               end if
            end do
         end do
      end do
   end do
   if (dbg_worstmax /= dbg_worstmax .or. dbg_worstmax > 1.0d12) then
      dbg0 = -1_8
   else
      dbg0 = nint(dbg_worstmax*1.0d6, kind=8)
   end if
   if (dbg_worstnorm /= dbg_worstnorm .or. dbg_worstnorm > 1.0d12) then
      dbg1 = -1_8
   else
      dbg1 = nint(dbg_worstnorm*1.0d9, kind=8)
   end if
   dbg2 = int(dbg_worstnn, 8)
   call geps_dbg_vram(6, 500, dbg0, dbg1, dbg2)
   ! #endregion
end subroutine eigen_mx_gpu
