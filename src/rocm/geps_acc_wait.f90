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
