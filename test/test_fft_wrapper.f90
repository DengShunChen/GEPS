!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024, NVIDIA CORPORATION & AFFILIATES. All rights reserved.
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_fft_wrapper
   call cufft_rfftmlt_unit(1, 720, 70, 120, -1)
   call cufft_rfftmlt_unit(1, 720, 30, 150, 1)
end program test_fft_wrapper

subroutine cufft_rfftmlt_unit(inc, jump, n, m, isign)
   integer :: inc, jump, n, m, isign, i, j
   real(8), dimension(jump, m) :: a, a_cufft, work, work_gpu
   real(8), dimension(4096) :: trigs
   integer, dimension(19) :: ifax

   call random_seed()
   call random_number(a)
   if (isign .eq. 1) then
      do j = 1, m
         do i = 1, n + 2
            work_gpu(i, j) = a(i, j)
         end do
      end do
   end if
   a_cufft = a

   call rfftmlt(a, work, trigs, ifax, inc, jump, n, m, isign)
   call rfftmlt_gpu(a_cufft, work_gpu, trigs, ifax, inc, jump, n, m, isign)

   if (all(abs(a - a_cufft) <= 1e-10)) then
      PRINT *, "test_fft_wrapper passed."
   else
      PRINT *, "test_fft_wrapper failed."
      call exit(1)
   end if

end subroutine cufft_rfftmlt_unit
