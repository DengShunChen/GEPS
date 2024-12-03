!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_cufft_loop
   call cufft_loop_unit(45, 2560, 720, 120, 1)
end program test_cufft_loop

subroutine cufft_loop_unit(jlistnum, my_max, jump, m, isign)
   use const, only: RTYPE

   integer :: jlistnum, my_max, jump, m, isign
   real(kind=RTYPE), dimension(jump, m, my_max) :: cc
   real(kind=RTYPE), dimension(jump, m, my_max) :: cc_cufft
   real(kind=RTYPE), dimension(jump, m, my_max) :: work
   real, dimension(4096, 2560) :: trigsj
   real, dimension(19, 2560) :: ifaxj
   integer, dimension(2560) :: jlist1
   integer, dimension(2560) :: nxdef
   real :: temp

   do i = 1, jlistnum
      jlist1(i) = jlistnum + 1 - i
   end do

   do i = 1, jlistnum
      call random_number(temp)
      nxdef(i) = FLOOR((jump - 2)*temp)
   end do

   call random_number(cc)
   cc_cufft = cc

   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef(j)
#ifdef SP
      call rfftmlt_sp(cc(1, 1, jj), work(1, 1, jj), trigsj(1, j), ifaxj(1, j), 1, jump, nxj, m, isign)
#else
      call rfftmlt(cc(1, 1, jj), work(1, 1, jj), trigsj(1, j), ifaxj(1, j), 1, jump, nxj, m, isign)
#endif
   end do

   !$acc enter data copyin(cc_cufft) create(work) async(1)
   call rfftmlt_loop(cc_cufft, work, trigsj, ifaxj, jlist1, nxdef, jlistnum, jump, m, isign)
   !$acc exit data copyout(cc_cufft) delete(work) async(1)
   !$acc wait(1)

#ifdef SP
   call assert_allclose(cc_cufft, size(cc_cufft), cc, size(cc), 1e-4_4, 1e-4_4, "Array cc")
#else
   call assert_allclose(cc_cufft, size(cc_cufft), cc, size(cc), 1e-10_8, 1e-10_8, "Array cc")
#endif

end subroutine cufft_loop_unit

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
