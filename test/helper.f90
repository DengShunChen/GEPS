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
      write (*, '(2(1A, 1pe12.4))') &
         "Arrays are not close within tolerance rtol =", rtol, &
         ", atol =", atol
      write (*, '(1A, 1pe15.7)') "Max relative difference = ", rel_diff
      write (*, '(1A, 1pe15.7)') "Max absolute difference = ", abs_diff
      if (present(err_msg)) print *, err_msg
      call exit(1)
   end if

end subroutine assert_allclose

! ============================================================
subroutine assert_realspace_close(actual, desired, ld, lev, ncld, &
                                  rtol, atol, err_msg)
   use const, only: RTYPE
   use param, only: my_max
   use index, only: nxp, nxdef, jlist1, jlistnum

   implicit none

   integer, intent(in) :: ld, lev, ncld
   real(kind=RTYPE), dimension(ld, lev*ncld, my_max), intent(in) :: actual, desired
   real(kind=RTYPE), intent(in) :: rtol, atol
   character(len=*), intent(in), optional :: err_msg
   real(kind=RTYPE), parameter :: eps = 1e-15
   real(kind=RTYPE) :: rel_diff, abs_diff, loc_diff
   logical :: equal
   integer :: i, j, k, n, nxj, jj, kk, nk

   equal = .true.
   rel_diff = 0.0
   abs_diff = 0.0
   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef(j)
      do n = 1, ncld
         nk = (n - 1)*lev
         do k = 1, lev
            kk = k + nk
            do i = 1, nxj
               loc_diff = abs(actual(i, kk, jj) - desired(i, kk, jj))
               rel_diff = max(rel_diff, loc_diff/abs(desired(i, kk, jj) + eps))
               abs_diff = max(abs_diff, loc_diff)
               if (loc_diff > atol + rtol*abs(desired(i, kk, jj))) then
                  write (*, '(I4,1X,I3,1X,I2,1X,I2,1X,4(1pe15.7, 1X))') &
                     i, j, k, n, &
                     actual(i, kk, jj), desired(i, kk, jj), &
                     loc_diff, loc_diff/abs(desired(i, kk, jj) + eps)
                  equal = .false.
               end if
            end do
         end do
      end do
   end do

   if (.not. equal) then
      write (*, '(2(1A, 1pe12.4))') &
         "Arrays are not close within tolerance rtol =", rtol, &
         ", atol =", atol
      write (*, '(1A, 1pe15.7)') "Max relative difference = ", rel_diff
      write (*, '(1A, 1pe15.7)') "Max absolute difference = ", abs_diff
      if (present(err_msg)) print *, err_msg
      ! call exit(1)
   end if

end subroutine assert_realspace_close
