
subroutine assert_real(actual, n_actual, desired, n_desired, rtol, err_msg)
   use const, only: RTYPE
   use rank

   implicit none

   integer, intent(in) :: n_actual, n_desired
   real(kind=RTYPE), dimension(n_actual), intent(in) :: actual
   real(kind=RTYPE), dimension(n_desired), intent(in) :: desired
   real(kind=RTYPE), intent(in) :: rtol
   character(len=*), intent(in), optional :: err_msg
   real(kind=RTYPE), parameter :: eps = 1e-300
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
      if (present(err_msg)) print *, err_msg, myrank
      !call exit(1)
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
      !call exit(1)
   end if

end subroutine assert_integer

subroutine adjptqintp_unit(itimestep, benchmark, ct, gt)
   use index, only: jlistnum, jlist1, nxdef_2d, nxjp_acc, nxjp, nxp
   use param, only: my_max, my, ncld, lev
   use machine, only: kind_phys
   use rank, only: myrank
   use const, only: RTYPE, dsigma
   !use nvtx

   implicit none

   integer :: ntrac, check, myim(my_max)
   integer :: i, ii, j, jj, async_id, n, k, iii, cnt, tcnt, nxj
   real :: dta, ptop
   integer, dimension(34) :: seed
   real :: time1, time2
   character(len=4) :: flag_str
   !flags
   integer, intent(in) :: itimestep
   logical, intent(in) :: benchmark
   real, intent(inout) :: ct, gt
   real(kind=RTYPE), dimension(nxp, my_max) :: pst, ps
   real(kind=RTYPE), dimension(nxp, lev, my_max) :: ut, vt, tt
   real(kind=RTYPE), dimension(nxp, lev*ncld, my_max) :: qt, qtp
   real(kind=RTYPE), dimension(nxp, my_max) :: pst_cpu, ps_cpu
   real(kind=RTYPE), dimension(nxp, lev, my_max) :: ut_cpu, vt_cpu, tt_cpu
   real(kind=RTYPE), dimension(nxp, lev*ncld, my_max) :: qt_cpu
   real(kind=RTYPE), dimension(nxp, my_max) :: pst_gpu, ps_gpu
   real(kind=RTYPE), dimension(nxp, lev, my_max) :: ut_gpu, vt_gpu, tt_gpu
   real(kind=RTYPE), dimension(nxp, lev*ncld, my_max) :: qt_gpu
   
   
   write(flag_str,'(I3)') itimestep
   async_id = 1
   seed = (/10006, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,&
           &25, 26, 27, 28, 29, 30, 31, 32, 33/)
   !call random_seed()
   call random_seed(put=seed)
   call random_number(ut)
   call random_number(vt)
   call random_number(tt)
   call random_number(qt)
   call random_number(pst)
   call random_number(ps)
   call random_number(qtp)
   
   cnt = 0
   tcnt = 0
   
   
   do jj = 1, jlistnum
      j = jlist1(jj)
      myim(jj) = nxjp(j)
   end do
   
   ut = ut*60. - 30.
   vt = vt*40. - 20.
   tt = tt*100. + 200.
   qt = qt*2.e-2 + 2.e-7
   qtp = qtp*2.e-2 + 2.e-7
   ps = ps*1000
   pst = pst*1000

   ut_gpu = ut
   vt_gpu = vt
   tt_gpu = tt
   qt_gpu = qt
   pst_gpu = pst
   ps_gpu = ps
   ut_cpu = ut
   vt_cpu = vt
   tt_cpu = tt
   qt_cpu = qt
   pst_cpu = pst
   ps_cpu = ps
   





   !if (benchmark .eq. .true.) call nvtxStartRange("CPU compute")
      call cpu_time(time1)
      do jj = 1, jlistnum
         j = jlist1(jj)
        call adjptqintp(ut_cpu(1,1,jj),vt_cpu(1,1,jj),tt_cpu(1,1,jj),qt_cpu(1,1,jj),  &
                        qtp(1,1,jj),pst_cpu(1,jj),ps_cpu(1,jj),nxjp(j),nxp,my_max,    &
                        lev,ncld,dta)
      end do
      call cpu_time(time2)
      if (benchmark .eq. .true.) ct = ct + time2 - time1
   !if (benchmark .eq. .true.) call nvtxEndRange


   if (.true.) then
         !$acc wait(async_id)
         !$acc enter data copyin(jlist1, nxjp, dsigma) async(async_id)
         !$acc enter data copyin(ut_gpu, vt_gpu, tt_gpu, qt_gpu, qtp, pst_gpu, ps_gpu) async(async_id)
         call cpu_time(time1)
         !if (benchmark .eq. .true.) call nvtxStartRange("GPU compute")
         call adjptqintp_gpu(ut_gpu, vt_gpu, tt_gpu, &
                         qt_gpu, qtp, pst_gpu, ps_gpu, &
                         nxjp, nxp, my_max, lev, ncld, dta)
         !$acc wait(async_id)
         !if (benchmark .eq. .true.) call nvtxEndRange
         call cpu_time(time2)
         if (benchmark .eq. .true.) gt = gt + time2 - time1
         !$acc exit data copyout(ut_gpu, vt_gpu, tt_gpu, qt_gpu, pst_gpu, ps_gpu) delete(qtp) async(async_id)
         !$acc exit data delete(jlist1, nxjp, dsigma) async(async_id)
         !$acc wait(async_id)
         
   call assert_real(ut_gpu, size(ut_gpu), ut_cpu, size(ut_cpu), &
                    1e-8, "Array ut"//trim(flag_str))
   call assert_real(vt_gpu, size(vt_gpu), vt_cpu, size(vt_cpu), &
                    1e-8, "Array vt"//trim(flag_str))
   call assert_real(tt_gpu, size(tt_gpu), tt_cpu, size(tt_cpu), &
                    1e-8, "Array tt"//trim(flag_str))
   call assert_real(qt_gpu, size(qt_gpu), qt_cpu, size(qt_cpu), &
                    1e-8, "Array qt"//trim(flag_str))
   call assert_real(pst_gpu, size(pst_gpu), pst_cpu, size(pst_cpu), &
                    1e-8, "Array pst"//trim(flag_str))
   call assert_real(ps_gpu, size(ps_gpu), ps_cpu, size(ps_cpu), &
                    1e-8, "Array ps"//trim(flag_str))
   end if
   


end subroutine adjptqintp_unit

program test_adjptqintp
   use index, only: nxp
   use rank, only: myrank
   use param
   use const

   implicit none
   integer :: i, no
   real :: ct, gt
   
   
   ct = 0.
   gt = 0.
   
   call mpe_init
   call cons
   
   call adjptqintp_unit(3, .false., ct, gt)
   do i = 1, 2
      call adjptqintp_unit(3, .true., ct, gt)
   end do
   
   write (*, *) 'Timing for CPU & GPU: ', ct, gt
   write (*, *) 'speedup ratio: ', ct/gt, myrank

   
   call mpe_finalize

end program

