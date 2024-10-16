!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024,
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_mpe2d_unify_nx
  use param
  use const

  implicit none

  call mpe_init
  call cons
  call mpe2d_unify_nx_unit
  call mpe_finalize

end program test_mpe2d_unify_nx

subroutine mpe2d_unify_nx_unit
  use param
  use const, only: RTYPE
  use index
  use rank, only: myrank

  implicit none

  integer, parameter :: steps = 10
  real(kind=RTYPE), dimension(nxp, my_max) :: inp
  real(kind=RTYPE), dimension(nx, my_max) :: out, out_gpu
  integer :: async_id, i

  async_id = 1
  if (myrank.eq.0) print*, "Start mpe2d_unify_nx unit test."
  print*, "nxp: ",nxp

  call random_seed()
  call random_number(inp)
  inp = inp*(myrank+1)

  do i = 1, steps
     call mpe2d_unify_nx(out, inp)
  end do
  !$acc enter data create(out_gpu) copyin(nxjlen_all, inp) async(async_id)
  do i = 1, steps
     call mpe2d_unify_nx_gpu(out_gpu, inp)
  end do
  !$acc exit data delete(nxjlen_all, inp) copyout(out_gpu) async(async_id)
  !$acc wait(async_id)


  write(*,'(2(A, 1pe15.7, 1X))') &
       "CPU: ", maxval(out), &
       "GPU: ", maxval(out_gpu)
  if (all(abs(out - out_gpu) <= 1e-10)) then
     print *, "(all close) test_mpe2d_unify_nx passed."
  else
     print *, "(all close) test_mpe2d_unify_nx failed."
     call exit(1)
  end if

end subroutine mpe2d_unify_nx_unit
