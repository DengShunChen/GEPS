!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024,
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_ptotc
  use param
  use const

  implicit none

  call mpe_init
  call cons
  call ptotc_unit
  call mpe_finalize

end program test_ptotc

subroutine ptotc_unit
  use param
  use const, only: RTYPE, dsigma, cosl
  use index, only: jlist1, nxdef_2d, nxjlen_all, jlist2, nxdef
  use rank, only: myrank
  use grid, only: pt, qt

  implicit none

  integer, parameter :: steps = 10
  real :: pdrym, pdrym_gpu
  real(kind=RTYPE), dimension
  integer :: async_id, i
  logical, parameter :: lprint = .false.

  async_id = 1
  if (myrank.eq.0) print*, "Start ptotc unit test."

  call random_seed()
  call random_number(pt)
  call random_number(qt)

  pdrym = 0.
  pdrym_gpu = 0.

  do i = 1, steps
     call ptotc(pdrym, lprint)
  end do
  !$acc enter data async(async_id) &
  !$acc& copyin(jlist1, nxdef, nxdef_2d, dsigma, pt, qt, nxjlen_all, cosl, jlist2)
  do i = 1, steps
     call ptotc_gpu(pdrym_gpu, lprint)
  end do
  !$acc exit data async(async_id) &
  !$acc& delete(jlist1, nxdef, nxdef_2d, dsigma, pt, qt, nxjlen_all, cosl, jlist2)
  !$acc wait(async_id)

#ifdef SP
  if (abs(pdrym - pdrym_gpu) <= 1e-4) then
#else
  if (abs(pdrym - pdrym_gpu) <= 1e-10) then
#endif
     print *, "(all close) test_ptotc passed."
  else
     print *, "(all close) test_ptotc failed."
     write(*,'(2(A, 1pe15.7, 1X))') &
          "CPU: ", pdrym, &
          "GPU: ", pdrym_gpu

     call exit(1)
  end if

end subroutine ptotc_unit
