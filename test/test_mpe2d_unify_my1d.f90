!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024,
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_mpe2d_unify_my1d
  use param
  use const

  implicit none

  call mpe_init
  call cons
  call mpe2d_unify_my1d_unit
  call mpe_finalize

end program test_mpe2d_unify_my1d

subroutine mpe2d_unify_my1d_unit
  use param
  use const, only: RTYPE
  use index
  use rank, only: myrank

  implicit none

  integer, parameter :: steps = 10
  real(kind=RTYPE), dimension(my_max) :: inp
  real(kind=RTYPE), dimension(my) :: out, out_gpu
  integer :: async_id, i

  async_id = 1
  if (myrank.eq.0) print*, "Start mpe2d_unify_my1d unit test."

  call random_seed()
  call random_number(inp)

  out= 0.
  out_gpu = 0.

  do i = 1, steps
     call mpe2d_unify_my1d(out, inp)
  end do
  !$acc enter data create(out_gpu) copyin(jlist2, inp) async(async_id)
  do i = 1, steps
     call mpe2d_unify_my1d_gpu(out_gpu, inp)
  end do
  !$acc exit data delete(jlist2, inp) copyout(out_gpu) async(async_id)
  !$acc wait(async_id)

  if (all(abs(out - out_gpu) <= 1e-10)) then
     print *, "(all close) test_mpe2d_unify_my1d passed."
  else
     print *, "(all close) test_mpe2d_unify_my1d failed."
     call exit(1)
  end if

end subroutine mpe2d_unify_my1d_unit
