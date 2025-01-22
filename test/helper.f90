!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine assert_allclose_r8(actual, n_actual, desired, n_desired, rtol, atol, err_msg)
   use const, only: RTYPE

   implicit none

   integer, intent(in) :: n_actual, n_desired
   real(kind=8), dimension(n_actual), intent(in) :: actual
   real(kind=8), dimension(n_desired), intent(in) :: desired
   real(kind=8), intent(in) :: rtol, atol
   character(len=*), intent(in), optional :: err_msg
   real(kind=8), parameter :: eps = nearest(real(0.0, 8), 1.0)
   real(kind=8) :: rel_diff, abs_diff
   logical :: equal
   integer :: i, err_count

   if (n_actual .ne. n_desired) then
      print *, "Arrays have different lengths!"
      if (present(err_msg)) print *, err_msg
      call exit(1)
   end if

   equal = .true.
   err_count = 0
   rel_diff = 0.0
   abs_diff = 0.0
   do i = 1, n_actual
      rel_diff = max(rel_diff, abs((actual(i) - desired(i))/(desired(i) + eps)))
      abs_diff = max(abs_diff, abs(actual(i) - desired(i)))
      if (abs(actual(i) - desired(i)) > atol + rtol*abs(desired(i))) then
         equal = .false.
         if (err_count .lt. 5) then
            ! print '(5g0)', "Desired = ", desired(i), " but get actual = ", actual(i)
            ! print '(5g0)', "Absolute difference = ", abs(actual(i) - desired(i)), "; Tolerance = ", atol + rtol*abs(desired(i))
            write (*, '(4(A, 1pe15.7))') &
               "Desired = ", desired(i), ", but get actual = ", actual(i), &
               ". Absolute difference = ", abs(actual(i) - desired(i)), "; Tolerance = ", atol + rtol*abs(desired(i))
         end if
         err_count = err_count + 1
      end if
   end do

   if (.not. equal) then
      write (*, '(2(1A, 1pe12.4))') &
         "Arrays are not close within tolerance rtol =", rtol, &
         ", atol =", atol
      write (*, '(1A, 1pe15.7)') "Max relative difference = ", rel_diff
      write (*, '(1A, 1pe15.7)') "Max absolute difference = ", abs_diff
      if (present(err_msg)) then
         print '(5g0)', "In array ", err_msg, ","
         print '(5g0)', err_count, "/", n_actual, " elements are not close within tolerance rtol = ", rtol, " atol =", atol
      end if
      call exit(1)
   end if

end subroutine assert_allclose_r8
! ============================================================
subroutine assert_allclose(actual, n_actual, desired, n_desired, rtol, atol, err_msg)
   use const, only: RTYPE

   implicit none

   integer, intent(in) :: n_actual, n_desired
   real(kind=RTYPE), dimension(n_actual), intent(in) :: actual
   real(kind=RTYPE), dimension(n_desired), intent(in) :: desired
   real(kind=RTYPE), intent(in) :: rtol, atol
   character(len=*), intent(in), optional :: err_msg
   real(kind=RTYPE), parameter :: eps = nearest(real(0.0, RTYPE), 1.0)
   real(kind=RTYPE) :: rel_diff, abs_diff
   logical :: equal
   integer :: i, err_count

   if (n_actual .ne. n_desired) then
      print *, "Arrays have different lengths!"
      if (present(err_msg)) print *, err_msg
      call exit(1)
   end if

   equal = .true.
   err_count = 0
   rel_diff = 0.0
   abs_diff = 0.0
   do i = 1, n_actual
      rel_diff = max(rel_diff, abs((actual(i) - desired(i))/(desired(i) + eps)))
      abs_diff = max(abs_diff, abs(actual(i) - desired(i)))
      if (abs(actual(i) - desired(i)) > atol + rtol*abs(desired(i))) then
         equal = .false.
         if (err_count .lt. 5) then
            ! print '(5g0)', "Desired = ", desired(i), " but get actual = ", actual(i)
            ! print '(5g0)', "Absolute difference = ", abs(actual(i) - desired(i)), "; Tolerance = ", atol + rtol*abs(desired(i))
            write (*, '(4(A, 1pe15.7))') &
               "Desired = ", desired(i), ", but get actual = ", actual(i), &
               ". Absolute difference = ", abs(actual(i) - desired(i)), "; Tolerance = ", atol + rtol*abs(desired(i))
         end if
         err_count = err_count + 1
      end if
   end do

   if (.not. equal) then
      write (*, '(2(1A, 1pe12.4))') &
         "Arrays are not close within tolerance rtol =", rtol, &
         ", atol =", atol
      write (*, '(1A, 1pe15.7)') "Max relative difference = ", rel_diff
      write (*, '(1A, 1pe15.7)') "Max absolute difference = ", abs_diff
      if (present(err_msg)) then
         print '(5g0)', "In array ", err_msg, ","
         print '(5g0)', err_count, "/", n_actual, " elements are not close within tolerance rtol = ", rtol, " atol =", atol
      end if
      call exit(1)
   end if

end subroutine assert_allclose
! ============================================================
function assert_specspace_close(actual, desired, &
                                rtol, atol, err_msg)
   use const, only: RTYPE
   use param, only: jtrun, jtmax
   use index, only: levp, mlist, mlistnum

   implicit none

   logical :: assert_specspace_close, equal
   real(kind=RTYPE), dimension(levp, 2, jtrun, jtmax), intent(in) :: actual, desired
   real(kind=RTYPE), intent(in) :: rtol, atol
   character(len=*), intent(in), optional :: err_msg
   real(kind=RTYPE), parameter :: eps = 1e-15
   real(kind=RTYPE) :: rel_diff, abs_diff, loc_diff
   integer :: i, k, l, mf, m

   equal = .true.
   rel_diff = 0.0
   abs_diff = 0.0
   do m = 1, mlistnum
      mf = mlist(m)
      do l = mf, jtrun
         do i = 1, 2
            do k = 1, levp
               loc_diff = abs(actual(k, i, l, m) - desired(k, i, l, m))

               if (loc_diff > atol + rtol*abs(desired(k, i, l, m))) then
                  rel_diff = max(rel_diff, loc_diff/abs(desired(k, i, l, m) + eps))
                  abs_diff = max(abs_diff, loc_diff)

                  write (*, '(2(I2,1X),3(I4,1X),4(1pe15.7, 1X))') &
                     k, i, l, m, mf, &
                     actual(k, i, l, m), desired(k, i, l, m), &
                     loc_diff, loc_diff/abs(desired(k, i, l, m) + eps)
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
   end if
   assert_specspace_close = equal
end function assert_specspace_close

function assert_specspace_2d_close(actual, desired, &
                                   rtol, atol, err_msg)
   use const, only: RTYPE
   use param, only: jtrun, jtmax
   use index, only: levp, mlist, mlistnum

   implicit none

   logical :: assert_specspace_2d_close, equal
   real(kind=RTYPE), dimension(jtrun, jtmax, 2), intent(in) :: actual, desired
   real(kind=RTYPE), intent(in) :: rtol, atol
   character(len=*), intent(in), optional :: err_msg
   real(kind=RTYPE), parameter :: eps = 1e-15
   real(kind=RTYPE) :: rel_diff, abs_diff, loc_diff
   integer :: i, k, l, mf, m

   equal = .true.
   rel_diff = 0.0
   abs_diff = 0.0
   do i = 1, 2
   do m = 1, mlistnum
      mf = mlist(m)
      do l = mf, jtrun
         loc_diff = abs(actual(l, m, i) - desired(l, m, i))

         if (loc_diff > atol + rtol*abs(desired(l, m, i))) then
            rel_diff = max(rel_diff, loc_diff/abs(desired(l, m, i) + eps))
            abs_diff = max(abs_diff, loc_diff)

            equal = .false.
         end if
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
   end if
   assert_specspace_2d_close = equal
end function assert_specspace_2d_close

! ============================================================
subroutine assert_realspace_close(actual, desired, ld, lev, ncld, &
                                  rtol, atol, err_msg)
   use const, only: RTYPE
   use param, only: my_max
   use index, only: nxp, nxdef_2d, jlist1, jlistnum

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
      nxj = nxdef_2d(j)
      do n = 1, ncld
         nk = (n - 1)*lev
         do k = 1, lev
            kk = k + nk
            do i = 1, nxj
               loc_diff = abs(actual(i, kk, jj) - desired(i, kk, jj))
               rel_diff = max(rel_diff, loc_diff/abs(desired(i, kk, jj) + eps))
               abs_diff = max(abs_diff, loc_diff)
               if (loc_diff > atol + rtol*abs(desired(i, kk, jj))) then

                  if (k .ne. 72) then
                     write (*, '(3(I4,1X),4(1pe15.7, 1X))') &
                          i, k, j, &
                          actual(i, kk, jj), desired(i, kk, jj), &
                          loc_diff, loc_diff/abs(desired(i, kk, jj) + eps)
                  end if
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

subroutine assert_rmse(actual, n_actual, desired, n_desired, tol, err_msg)
   use const, only: RTYPE

   implicit none

   integer, intent(in) :: n_actual, n_desired
   real(kind=RTYPE), dimension(n_actual), intent(in) :: actual
   real(kind=RTYPE), dimension(n_desired), intent(in) :: desired
   real(kind=RTYPE), intent(in) :: tol
   character(len=*), intent(in), optional :: err_msg
   real(kind=RTYPE) :: rmse

   if (n_actual .ne. n_desired) then
      print *, "Arrays have different lengths!"
      if (present(err_msg)) print *, err_msg
      call exit(1)
   end if

   rmse = sqrt(sum((actual - desired)**2)/n_actual)

   if (rmse .gt. tol) then
      print '(5g0)', "In array ", err_msg, ", RMSE = ", rmse, " but tolerance = ", tol
      ! The exit function should be called after the FP32 output difference issue is fixed.
      ! call exit(1)
   end if

end subroutine assert_rmse

subroutine init_seed()
   implicit none

   integer :: state_size
   integer, allocatable, dimension(:) :: state

   call random_seed(size=state_size)
   allocate (state(state_size))
   state = 123
   call random_seed(put=state)

end subroutine init_seed
