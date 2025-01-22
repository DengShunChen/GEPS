subroutine mpe2d_row_bcast_unit(no)
   use rank, only: myrank
   use const, only: RTYPE, nnmivm
   use param, only: jtmax
   use index, only: Lstart, Lend, row_rank
   use cudafor
   use openacc
   use mod_vcopy

   implicit none

   integer, parameter :: steps = 3
   integer, no, brank
   real, dimension(no*no, 2*jtmax*nnmivm) :: mx, mx_c
   integer L, j

   real(kind=8) tm_1, tm_2, tm_use, mpi_wtime
   integer async_id, s

   async_id = 1

   call random_seed()
   call random_number(mx)
   mx_c = mx

   ! << CPU >>
   tm_1 = mpi_wtime()
   do L = 1, nnmivm
      brank = 0
      if ((L .ge. Lstart) .and. (L .le. Lend)) then
         brank = row_rank
      end if
      j = 1 + (L - 1)*jtmax*2
      call mpe2d_row_broadcast(mx_c(1, j), no*no*2*jtmax, brank)
   end do

   tm_use = mpi_wtime() - tm_1
   if (myrank .eq. 0) write (*, '(A, f15.3, A)') &
      "Elapsed time of CPU: ", tm_use*1e+6, " us"

   ! << GPU >>
   !$acc enter data copyin(mx) async(async_id)
   !$acc wait(async_id)
   do s = 1, steps
      tm_1 = mpi_wtime()
      do L = 1, nnmivm
         brank = 0
         if ((L .ge. Lstart) .and. (L .le. Lend)) then
            brank = row_rank
         end if
         j = 1 + (L - 1)*jtmax*2
         call mpe2d_row_broadcast_gpu(mx(1, j), no*no*2*jtmax, brank)
      end do

      !$acc wait(async_id)
      tm_use = mpi_wtime() - tm_1
      if (myrank .eq. 0) write (*, '(A, f15.3, A)') &
         "Elapsed time of GPU: ", tm_use*1e+6, " us"
   end do
   !$acc exit data copyout(mx) async(async_id)
   !$acc wait(async_id)

   call assert_allclose_r8(mx, size(mx), mx_c, size(mx_c), 1e-15, 1e-15, "mx")

end subroutine mpe2d_row_bcast_unit

program test_mpe2d_row_bcast
   use rank, only: myrank
   use param, only: jtrun
   implicit none
   integer no

   call mpe_init
   call cons
   no = 2*((jtrun + 1)/2) + (jtrun/2) + 10
   call mpe2d_row_bcast_unit(no)
   call mpe_finalize

end program test_mpe2d_row_bcast
