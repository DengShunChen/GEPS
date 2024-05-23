!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine assert_allclose(actual, n_actual, desired, n_desired, rtol, atol, err_msg)
   use const, only: RTYPE

   implicit none

   integer, intent(in) :: n_actual, n_desired
   real(kind=RTYPE), dimension(n_actual), intent(in) :: actual
   real(kind=RTYPE), dimension(n_desired), intent(in) :: desired
   real(kind=RTYPE), intent(in) :: rtol, atol
   character(len=*), intent(in), optional :: err_msg
   real(kind=RTYPE), parameter :: eps = 1e-15
   real(kind=RTYPE) :: rel_diff, abs_diff
   logical :: equal
   integer :: i

   if (n_actual .ne. n_desired) then
      print *, "Arrays have different lengths!"
      if (present(err_msg)) print *, err_msg
      call exit(1)
   end if

   equal = .true.
   rel_diff = 0.0
   abs_diff = 0.0
   do i = 1, n_actual
      rel_diff = max(rel_diff, abs((actual(i) - desired(i))/(desired(i) + eps)))
      abs_diff = max(abs_diff, abs(actual(i) - desired(i)))
      if (abs(actual(i) - desired(i)) > atol + rtol*abs(desired(i))) equal = .false.
   end do

   if (.not. equal) then
      print *, "Arrays are not close within tolerance rtol =", rtol, "atol =", atol
      print *, "Max relative difference = ", rel_diff
      print *, "Max absolute difference = ", abs_diff
      if (present(err_msg)) print *, err_msg
      call exit(1)
   end if

end subroutine assert_allclose
