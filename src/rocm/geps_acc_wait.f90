! External wait used by acc2omp-translated `!$acc wait` (no module required).
subroutine geps_acc_wait(async_id)
  use iso_c_binding, only: c_int
  implicit none
  integer, intent(in) :: async_id
  interface
      subroutine geps_acc_wait_c(id) bind(C, name="geps_hip_wait")
      import c_int
      integer(c_int), value :: id
    end subroutine
  end interface
  call geps_acc_wait_c(int(async_id, c_int))
end subroutine geps_acc_wait

! Device-wide wait, inserted by acc2omp after cudaMemsetAsync/cudaMemcpyAsync
! runs (see geps_hip_wait_all in hip_compat.cc for why).
subroutine geps_acc_wait_all()
  implicit none
  interface
      subroutine geps_acc_wait_all_c() bind(C, name="geps_hip_wait_all")
    end subroutine
  end interface
  call geps_acc_wait_all_c()
end subroutine geps_acc_wait_all


! wall clock for the acc2omp phase timers (GEPS_ACC2OMP_TIMERS=1)
function geps_wtime() result(t)
   use, intrinsic :: iso_fortran_env, only: int64
   implicit none
   real(kind=8) :: t
   integer(int64) :: c, r
   call system_clock(c, r)
   t = real(c, 8)/real(r, 8)
end function geps_wtime

! MPI rank for the acc2omp phase timers (avoids a `use rank` inside translated units)
function geps_hide_myrank_t() result(r)
   use rank, only: myrank
   implicit none
   integer :: r
   r = myrank
end function geps_hide_myrank_t

! 2026-09-19: turn the shim's per-dgemm stream sync off/on around the Legendre
! dgemm loops (acc2omp _defer_graph_loop_syncs); the loop is followed by
! geps_acc_wait_all().
subroutine geps_blas_defer_sync(on)
  use iso_c_binding, only: c_int
  implicit none
  integer, intent(in) :: on
  interface
      subroutine geps_blas_defer_sync_c(v) bind(C, name="geps_blas_defer_sync_set")
      import c_int
      integer(c_int), value :: v
    end subroutine
  end interface
  call geps_blas_defer_sync_c(int(on, c_int))
end subroutine geps_blas_defer_sync
