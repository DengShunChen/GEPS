!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_cufft_plan_cache
   call cufft_plan_cache_unit(1, 720, 70, 120, 1)
end program test_cufft_plan_cache

subroutine cufft_plan_cache_unit(inc, jump, n, m, isign)
   integer :: inc, jump, n, m, isign
   real(8), dimension(jump, m) :: pa
   real(8), dimension(jump, m) :: a
   integer(4) :: plan, cached_plan
   integer(kind=int_ptr_kind()) :: work_size

   interface c_interface
      subroutine find_fft_plan(inc, jump, n, m, isign, plan) bind(C, name="find_fft_plan")
         use iso_c_binding
         implicit none
         integer(c_int), value, intent(in) :: inc, jump, n, m, isign
         integer(c_int), intent(out) :: plan
      end subroutine find_fft_plan
      subroutine cache_fft_plan(inc, jump, n, m, isign, plan) bind(C, name="cache_fft_plan")
         use iso_c_binding
         implicit none
         integer(c_int), value, intent(in) :: inc, jump, n, m, isign, plan
      end subroutine cache_fft_plan
   end interface c_interface

   call find_fft_plan(inc, jump, n, m, isign, plan)
   if (plan .ne. -1) then
      PRINT *, "plan_table should be empty."
      call exit(1)
   end if

   call fft_create_plan(plan, 1)
   call fft_make_plan(inc, jump, n, m, isign, plan, work_size)
   call cache_fft_plan(inc, jump, n, m, isign, plan)

   call find_fft_plan(inc, jump, n, m, isign, cached_plan)
   if (plan .ne. cached_plan) then
      PRINT *, "cached plan id is not correct."
      call exit(1)
   end if
   PRINT *, "test_cufft_plan_cache passed."

end subroutine cufft_plan_cache_unit
