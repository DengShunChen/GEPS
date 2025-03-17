!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024,
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_mpe2d_unify_lev
   use param
   use const

   implicit none

   call mpe_init
   call cons
   call mpe2d_unify_lev_unit
   call mpe2d_unify_spec_lev_zx_unit
   call mpe_finalize

end program test_mpe2d_unify_lev

subroutine mpe2d_unify_lev_unit
   use param
   use const, only: RTYPE
   use index
   use rank, only: myrank

   implicit none

   integer, parameter :: steps = 10
   real(kind=RTYPE), dimension(levp, 2, jtrun, jtmax) :: inp
   real(kind=RTYPE), dimension(lev, 2, jtrun, jtmax) :: out, out_gpu
   real(kind=RTYPE), atol
   integer :: async_id, i

   async_id = 1
   if (myrank .eq. 0) print *, "Start mpe2d_unify_lev unit test."

   call random_seed()
   call random_number(inp)
   inp = inp*(myrank + 1)

   do i = 1, steps
      call mpe2d_unify_lev(inp, out, lev, levp, jtrun, jtmax, mlistnum, &
                           nsizex, row_comm)
   end do
   !$acc enter data create(out_gpu) copyin(inp) async(async_id)
   do i = 1, steps
      call mpe2d_unify_lev_gpu(inp, out_gpu, lev, levp, jtrun, jtmax, mlistnum, &
                               nsizex, nccl_row_comm)
   end do
   !$acc exit data delete(inp) copyout(out_gpu) async(async_id)
   !$acc wait(async_id)

   write (*, '(2(A, i0.3, A, 1pe15.7, "; "))') &
      "CPU-", myrank, ": ", maxval(out), &
      "GPU-", myrank, ": ", maxval(out_gpu)
#ifdef SP
   atol = 1e-4
#else
   atol = 1e-10
#endif
   if (all(abs(out - out_gpu) <= atol)) then
      print *, "(all close) test_mpe2d_unify_lev passed."
   else
      print *, "(all close) test_mpe2d_unify_lev failed."
      call exit(1)
   end if

end subroutine mpe2d_unify_lev_unit

subroutine mpe2d_unify_spec_lev_zx_unit
   use param
   use const, only: RTYPE
   use index
   use rank, only: myrank

   implicit none

   integer, parameter :: steps = 10
   real(kind=RTYPE), dimension(levp, 2, jtrun, jtmax, 3) :: inp
   real(kind=RTYPE), dimension(lev, 2, 3, jtrun, jtmax) :: out, out_gpu
   real(kind=RTYPE), atol
   integer :: async_id, i

   async_id = 1
   if (myrank .eq. 0) print *, "Start mpe2d_unify_spec_lev_zx unit test."

   call random_seed()
   call random_number(inp)
   inp = inp*(myrank + 1)

   do i = 1, steps
      call mpe2d_unify_spec_lev_zx(out, inp(1, 1, 1, 1, 1), &
                                   inp(1, 1, 1, 1, 2), inp(1, 1, 1, 1, 3), &
                                   lev, levp, jtrun, jtmax, mlistnum, &
                                   nsizex, row_comm)
   end do
   !$acc enter data create(out_gpu) copyin(inp) async(async_id)
   do i = 1, steps
      call mpe2d_unify_spec_lev_zx_gpu(out_gpu, inp(1, 1, 1, 1, 1), &
                                       inp(1, 1, 1, 1, 2), inp(1, 1, 1, 1, 3), &
                                       lev, levp, jtrun, jtmax, mlistnum, &
                                       nsizex, nccl_row_comm)
   end do
   !$acc exit data delete(inp) copyout(out_gpu) async(async_id)
   !$acc wait(async_id)

   write (*, '(2(A, i0.3, A, 1pe15.7, "; "))') &
      "CPU-", myrank, ": ", maxval(out), &
      "GPU-", myrank, ": ", maxval(out_gpu)
#ifdef SP
   atol = 1e-4
#else
   atol = 1e-10
#endif
   if (all(abs(out - out_gpu) <= atol)) then
      print *, "(all close) test_mpe2d_unify_spec_lev_zx passed."
   else
      print *, "(all close) test_mpe2d_unify_spec_lev_zx failed."
      call exit(1)
   end if

end subroutine mpe2d_unify_spec_lev_zx_unit
