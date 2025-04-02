program test_fpvs
   implicit none
   call mpe_init
   call cons
   call fpvs_unit
   call mpe_finalize

end program

subroutine assert_real(actual, n_actual, desired, n_desired, rtol, err_msg)
   use const, only: RTYPE

   implicit none

   integer, intent(in) :: n_actual, n_desired
   real(kind=RTYPE), dimension(n_actual), intent(in) :: actual
   real(kind=RTYPE), dimension(n_desired), intent(in) :: desired
   real(kind=RTYPE), intent(in) :: rtol
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

   rel_diff = maxval(abs((actual - desired)/(desired + eps)))
   abs_diff = maxval(abs(actual - desired))
   !if (abs_diff > atol) equal = .false.
   if (rel_diff > rtol) equal = .false.

   if (.not. equal) then
      print *, "Arrays are not close within tolerance rtol =", rtol
      print *, "Max relative difference = ", rel_diff
      !print *, "Max absolute difference = ", abs_diff
      if (present(err_msg)) print *, err_msg
      call exit(1)
   end if

end subroutine assert_real

subroutine assert_integer(actual, n_actual, desired, n_desired, err_msg)
   use const, only: RTYPE

   implicit none

   integer, intent(in) :: n_actual, n_desired
   integer, dimension(n_actual), intent(in) :: actual
   integer, dimension(n_desired), intent(in) :: desired
   character(len=*), intent(in), optional :: err_msg
   integer :: abs_diff
   logical :: equal
   integer :: i

   if (n_actual .ne. n_desired) then
      print *, "Arrays have different lengths!"
      if (present(err_msg)) print *, err_msg
      call exit(1)
   end if

   equal = .true.
   abs_diff = maxval(abs(actual - desired))
   if (abs_diff .ne. 0) equal = .false.

   if (.not. equal) then
      print *, "Arrays are not close to 0."
      print *, "Max absolute difference = ", abs_diff
      if (present(err_msg)) print *, err_msg
      call exit(1)
   end if

end subroutine assert_integer

subroutine fpvs_unit
   !$acc routine(fpvs_gpu) seq
   use radn, only: ntcw
   use index, only: nxjp, nxp, jlistnum, jlist1, nxjp_acc
   use param, only: lev, my_max, my, ncld
   use const, only: grav, cp, nx, mtnvar, tofd, cmbk, cgwd
   use physcons, only: con_rd, con_rv
   use machine, only: kind_phys
   use const, only: RTYPE
   use rank, only: myrank
   !use nvtx

   implicit none

   integer :: ntrac
   integer :: i, ii, j, jj, async_id, n, k
   integer, parameter :: nxpvs=7501
   integer, dimension(34) :: seed

   real :: time1, time2, ct, gt, diff
   real :: c1xpvs,c2xpvs,tbpvs(nxpvs)
   real(kind=kind_phys) fpvs, fpvs_gpu

!-------------- inputs
   integer, dimension(my_max) :: myim
   real, dimension(nxp, my_max) :: t
!-------------- inputs/outputs
!--------------outputs
   real(kind=kind_phys), dimension(nxp, my_max) :: tout, tout_gpu
   
   async_id = 1
   seed = (/10006, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,&
           &25, 26, 27, 28, 29, 30, 31, 32, 33/)
   !call random_seed()
   call random_seed(put=seed)
   ct = 0.
   gt = 0.
   do ii = 1, 16
      call random_number(t)

      t = t*20. + 273.15
      tout = 0.
      tout_gpu = 0.

      
      do jj = 1, jlistnum
         j = jlist1(jj)
         myim(jj) = nxjp(j)
      end do


      call cpu_time(time1)
      !if (ii .ne. 1) call nvtxStartRange("CPU compute")
      do jj = 1, jlistnum
         do i = 1, myim(jj)
            tout(i, jj) = fpvs(t(i, jj))
         end do
      end do
      call cpu_time(time2)
      !if (ii .ne. 1) call nvtxEndRange
      if (ii .ne. 1) ct = ct + time2 - time1


      if (.true.) then
         !$acc enter data copyin(tbpvs, t, tout_gpu) async(async_id)
         call cpu_time(time1)
         call gpvs_gpu(c1xpvs,c2xpvs,tbpvs,async_id)
         !if (ii .ne. 1) call nvtxStartRange("GPU compute")
         !$acc parallel loop gang private(jj,i) async(async_id)
         do jj = 1, jlistnum
            !$acc loop vector
            do i = 1, myim(jj)
               tout_gpu(i, jj) = fpvs_gpu(t(i, jj),c1xpvs,c2xpvs,tbpvs)
            end do
         end do
         !if (ii .ne. 1) call nvtxEndRange
         call cpu_time(time2)
         !$acc exit data copyout(tout_gpu) delete(t, tbpvs) async(async_id)
         !$acc wait(async_id)
         if (ii .ne. 1) gt = gt + time2 - time1


         call assert_real(tout_gpu, size(tout_gpu), tout, size(tout), &
                          1e-12, "Array tout")
      end if
   end do
   write (*, *) 'Timing for CPU & GPU: ', ct, gt
   write (*, *) 'speedup ratio: ', ct/gt, myrank

end subroutine fpvs_unit

