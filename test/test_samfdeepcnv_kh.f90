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

subroutine samfdeepcnv_kh_unit(itimestep, benchmark, ct, gt)
   use index, only: jlistnum, jlist1, nxdef_2d, nxjp_acc, nxjp, nxp
   use param, only: my_max, my, ncld, lev
   use machine, only: kind_phys
   use rank, only: myrank
   use const, only: dsigma, RTYPE, naero
   !use nvtx

   implicit none

   integer, parameter :: icheck = 499, kcheck = 16, jjcheck = 71, rcheck = 3
   integer :: ntrac, check, myim(my_max)
   integer :: i, ii, j, jj, async_id, n, k, iii, cnt, tcnt, nxj
   integer, dimension(34) :: seed
   real :: time1, time2
   character(len=4) :: flag_str
   !flags
   integer, intent(in) :: itimestep
   logical, intent(in) :: benchmark
   real, intent(inout) :: ct, gt
   real :: dta
   
   integer, parameter :: nxpvs = 7501
   real :: c1xpvs,c2xpvs,tbpvs(nxpvs)
   common/fpvscom/ c1xpvs,c2xpvs,tbpvs(nxpvs)
   
   real, dimension(nxp, lev, my_max) :: del, prsl, phil, qtr, qti, qtc, ttc, utc, &
      vtc, dotc, diss_dcc
   real, dimension(nxp, my_max) :: psfc, garea, islimsk_r
   integer, dimension(nxp, my_max) :: islimsk
   real,dimension(nxp, lev, my_max) :: qtr_gpu, qti_gpu, qtc_gpu, ttc_gpu, &
      utc_gpu, vtc_gpu, cnvw_gpu, cnvc_gpu, snow_flxn_gpu, ptun_gpu, pqun_gpu, &
      diss_dcc_gpu
   real,dimension(nxp, my_max) :: cldwrk_gpu, rcup_gpu 
   integer,dimension(nxp, my_max) :: kbot_gpu, ktop_gpu, kuo_gpu
   real,dimension(nxp, lev, my_max) :: qtr_cpu, qti_cpu, qtc_cpu, ttc_cpu, &
      utc_cpu, vtc_cpu, cnvw_cpu, cnvc_cpu, snow_flxn_cpu, ptun_cpu, pqun_cpu, &
      diss_dcc_cpu
   real,dimension(nxp, my_max) :: cldwrk_cpu, rcup_cpu 
   integer,dimension(nxp, my_max) :: kbot_cpu, ktop_cpu, kuo_cpu
   
   real :: sum_cpu, sum_gpu
   
   write(flag_str,'(I3)') itimestep
   async_id = 1
   seed = (/10006, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,&
           &25, 26, 27, 28, 29, 30, 31, 32, 33/)
   !call random_seed()
   call random_seed(put=seed)
   call random_number(del)
   call random_number(prsl)
   call random_number(psfc)
   call random_number(phil)
   call random_number(qtr)
   call random_number(qti)
   call random_number(qtc)
   call random_number(ttc)
   call random_number(utc)
   call random_number(vtc)
   call random_number(islimsk_r)
   call random_number(garea)
   call random_number(dotc)
   call random_number(diss_dcc)
   cnt = 0
   tcnt = 0
   !call read_data(myrank, itimestep)
   
   do jj = 1, jlistnum
      j = jlist1(jj)
      myim(jj) = nxjp(j)
   end do
   !write(*,*) "test itimestep = ", itimestep
   do jj = 1, jlistnum
      do i = 1, nxp
         if (islimsk_r(i, jj) .gt. 0.5) then
            islimsk(i, jj) = 1
         else
            islimsk(i, jj) = 0
         end if
      end do
   end do
   prsl = prsl*100000. + 100.
   del = del*3000. + 60.
   psfc = psfc*5. + 96.
   phil = phil*500000. + 220.
   qtr = qtr*4.e-4 + 1.e-20
   qti = qti*5.e-4 + 1.e-20
   qtc = qtc*2.e-2 + 2.e-7
   utc = utc*60. - 30.
   vtc = vtc*40. - 20.
   ttc = ttc*100. + 200.
   garea = garea*50000000. + 860000000.
   dotc = dotc*2.3e-3 + -1.8e-3 
   diss_dcc = diss_dcc*1.
   
   dta = 600.
   
   
   
   qtr_gpu = qtr
   qti_gpu = qti
   qtc_gpu = qtc
   ttc_gpu = ttc
   utc_gpu = utc
   vtc_gpu = vtc
   cldwrk_gpu = 0.
   rcup_gpu = 0.
   kbot_gpu = 0
   ktop_gpu = 0
   kuo_gpu = 0
   cnvw_gpu = 0.
   cnvc_gpu = 0.
   snow_flxn_gpu = 0.
   ptun_gpu = 0.
   pqun_gpu = 0.
   diss_dcc_gpu = diss_dcc

   qtr_cpu = qtr
   qti_cpu = qti
   qtc_cpu = qtc
   ttc_cpu = ttc
   utc_cpu = utc
   vtc_cpu = vtc
   cldwrk_cpu = 0.
   rcup_cpu = 0.
   kbot_cpu = 0
   ktop_cpu = 0
   kuo_cpu = 0
   cnvw_cpu = 0.
   cnvc_cpu = 0.
   snow_flxn_cpu = 0.
   ptun_cpu = 0.
   pqun_cpu = 0.
   diss_dcc_cpu = diss_dcc

   !if (benchmark .eq. .true.) call nvtxStartRange("CPU compute")
      do jj = 1, jlistnum
         call cpu_time(time1)
         j = jlist1(jj)
         nxj = nxdef_2d(j)
         !write(*,*) nxp, nxjp(j)
         call samfdeepcnv_kh(nxjp(j), nxp, lev, dta, del(1, 1, jj), &
                  prsl(1, 1, jj), psfc(1, jj), phil(1, 1, jj), &
                  qtr_cpu(1, 1, jj), qti_cpu(1, 1, jj), qtc_cpu(1, 1, jj), &
                  ttc_cpu(1, 1, jj), utc_cpu(1, 1, jj), vtc_cpu(1, 1, jj), &
                  cldwrk_cpu(1, jj), rcup_cpu(1, jj), kbot_cpu(1, jj), &
                  ktop_cpu(1, jj), kuo_cpu(1, jj), islimsk(1, jj), &
                  garea(1, jj), dotc(1, 1, jj), ncld, &
                  cnvw_cpu(1, 1, jj), cnvc_cpu(1, 1, jj), &
                  snow_flxn_cpu(1, 1, jj), ptun_cpu(1, 1, jj), & 
                  pqun_cpu(1, 1, jj), diss_dcc_cpu(1, 1, jj), &
                  jj, myrank, icheck, kcheck, jjcheck, rcheck)
         call cpu_time(time2)
         if (benchmark .eq. .true.) ct = ct + time2 - time1
      end do
   !if (benchmark .eq. .true.) call nvtxEndRange


   if (.false.) then
         !$acc wait(async_id)
         !$acc enter data copyin(jlist1, tbpvs) async(async_id)
         !$acc enter data copyin(nxjp, del, prsl, psfc, phil, islimsk, garea, &
         !$acc&      dotc) async(async_id)
         !$acc enter data copyin(qtr_gpu, qti_gpu, qtc_gpu, ttc_gpu, utc_gpu, &
         !$acc&      vtc_gpu, kuo_gpu, diss_dcc_gpu) async(async_id)
         !$acc enter data copyin(cldwrk_gpu, rcup_gpu, kbot_gpu, ktop_gpu, &
         !$acc&      cnvw_gpu, cnvc_gpu, snow_flxn_gpu, ptun_gpu, pqun_gpu) &
         !$acc&      async(async_id)
         call cpu_time(time1)
         !if (benchmark .eq. .true.) call nvtxStartRange("GPU compute")
         call samfdeepcnv_kh_gpu(nxjp, nxp, lev, dta, del, &
                  prsl, psfc, phil, &
                  qtr_gpu, qti_gpu, qtc_gpu, &
                  ttc_gpu, utc_gpu, vtc_gpu, &
                  cldwrk_gpu, rcup_gpu, kbot_gpu, &
                  ktop_gpu, kuo_gpu, islimsk, &
                  garea, dotc, ncld, &
                  cnvw_gpu, cnvc_gpu, &
                  snow_flxn_gpu, ptun_gpu, & 
                  pqun_gpu, diss_dcc_gpu, &
                  myrank, icheck, kcheck, jjcheck, rcheck)
         !$acc wait(async_id)
         !if (benchmark .eq. .true.) call nvtxEndRange
         call cpu_time(time2)
         if (benchmark .eq. .true.) gt = gt + time2 - time1
         !$acc exit data delete(jlist1, tbpvs) async(async_id)
         !$acc exit data delete(nxjp, del, prsl, psfc, phil, islimsk, garea, &
         !$acc&     dotc) async(async_id)
         !$acc exit data copyout(qtr_gpu, qti_gpu, qtc_gpu, ttc_gpu, utc_gpu, &
         !$acc&     vtc_gpu, kuo_gpu, diss_dcc_gpu) async(async_id)
         !$acc exit data copyout(cldwrk_gpu, rcup_gpu, kbot_gpu, ktop_gpu, &
         !$acc&     cnvw_gpu, cnvc_gpu, snow_flxn_gpu, ptun_gpu, pqun_gpu) &
         !$acc&     async(async_id)
         !$acc wait(async_id)
      
      cnt = 0
      tcnt = 0
   !   do jj = 1, jlistnum
   !      do i = 1, nxp
   !         do k = 1, lev
   !         if (abs((cnvw_gpu(i, k, jj) - cnvw_cpu(i, k, jj))/cnvw_cpu(i, k, jj)) .gt. 1e-6) then
   !            tcnt = tcnt + 1
   !            if (i .le. myim(jj)) then
   !               write(*,*) 'aaa', myrank, jj, k, i, cnvw_gpu(i, k, jj), cnvw_cpu(i, k, jj), &
   !                          abs((cnvw_gpu(i, k, jj) - cnvw_cpu(i, k, jj))/cnvw_cpu(i, k, jj))
   !               cnt = cnt + 1
   !            end if
   !         end if
   !         end do
   !      end do
   !   end do
   !   cnt = 0
   !   tcnt = 0
   !   do jj = 1, jlistnum
   !      do i = 1, nxp
   !      if ((kuo_gpu(i, jj) - kuo_cpu(i, jj)) .ne. 0) then
   !         tcnt = tcnt + 1
   !         if (i .le. myim(jj)) then
   !            write(*,*) 'aaa', myrank, jj, i, kuo_gpu(i, jj), kuo_cpu(i, jj)
   !            cnt = cnt + 1
   !         end if
   !      end if
   !      end do
   !   end do
   !   write(*,*) 'check', cnt, tcnt
   call check_real2D(rcup_gpu, rcup_cpu, 1e-8, 'rcup', cnt)
   call check_real3D(cnvw_gpu, cnvw_cpu, 1e-8, 'cnvw', cnt)
   call check_real3D(cnvc_gpu, cnvc_cpu, 1e-8, 'cnvc', cnt)
   call check_real2D(ptun_gpu, ptun_cpu, 1e-8, 'ptun', cnt)
   call check_real2D(pqun_gpu, pqun_cpu, 1e-8, 'pqun', cnt)
   call check_real3D(qtr_gpu, qtr_cpu, 1e-8, 'qtr', cnt)
   call check_real3D(qti_gpu, qti_cpu, 1e-8, 'qti', cnt)
   call check_real3D(qtc_gpu, qtc_cpu, 1e-8, 'qtc', cnt)
   call check_real3D(ttc_gpu, ttc_cpu, 1e-8, 'ttc', cnt)
   call check_real3D(utc_gpu, utc_cpu, 1e-8, 'utc', cnt)
   call check_real3D(vtc_gpu, vtc_cpu, 1e-8, 'vtc', cnt)
   call assert_real(qtr_gpu, size(qtr_gpu), qtr_cpu, size(qtr_cpu), &
                    1e-8, "Array qtr"//trim(flag_str))
   call assert_real(qti_gpu, size(qti_gpu), qti_cpu, size(qti_cpu), &
                    1e-8, "Array qti"//trim(flag_str))
   call assert_real(qtc_gpu, size(qtc_gpu), qtc_cpu, size(qtc_cpu), &
                    1e-8, "Array qtc"//trim(flag_str))
   call assert_real(ttc_gpu, size(ttc_gpu), ttc_cpu, size(ttc_cpu), &
                    1e-8, "Array ttc"//trim(flag_str))
   call assert_real(utc_gpu, size(utc_gpu), utc_cpu, size(utc_cpu), &
                    1e-8, "Array utc"//trim(flag_str))
   call assert_real(vtc_gpu, size(vtc_gpu), vtc_cpu, size(vtc_cpu), &
                    1e-8, "Array vtc"//trim(flag_str))
   call assert_real(cldwrk_gpu, size(cldwrk_gpu), cldwrk_cpu, size(cldwrk_cpu), &
                    1e-10, "Array cldwrk"//trim(flag_str))
   call assert_real(rcup_gpu, size(rcup_gpu), rcup_cpu, size(rcup_cpu), &
                    1e-10, "Array rcup"//trim(flag_str))
   call assert_integer(kbot_gpu, size(kbot_gpu), kbot_cpu, size(kbot_cpu), &
                    "Array kbot"//trim(flag_str))
   call assert_integer(ktop_gpu, size(ktop_gpu), ktop_cpu, size(ktop_cpu), &
                    "Array ktop"//trim(flag_str))
   call assert_integer(kuo_gpu, size(kuo_gpu), kuo_cpu, size(kuo_cpu), &
                    "Array kuo"//trim(flag_str))
   call assert_real(cnvw_gpu, size(cnvw_gpu), cnvw_cpu, size(cnvw_cpu), &
                    1e-10, "Array cnvw"//trim(flag_str))
   call assert_real(cnvc_gpu, size(cnvc_gpu), cnvc_cpu, size(cnvc_cpu), &
                    1e-10, "Array cnvc"//trim(flag_str))
   call assert_real(snow_flxn_gpu, size(snow_flxn_gpu), snow_flxn_cpu, size(snow_flxn_cpu), &
                    1e-10, "Array snow_flxn"//trim(flag_str))
   call assert_real(ptun_gpu, size(ptun_gpu), ptun_cpu, size(ptun_cpu), &
                    1e-10, "Array ptun"//trim(flag_str))
   call assert_real(pqun_gpu, size(pqun_gpu), pqun_cpu, size(pqun_cpu), &
                    1e-10, "Array pqun"//trim(flag_str))
   call assert_real(diss_dcc_gpu, size(diss_dcc_gpu), diss_dcc_cpu, size(diss_dcc_cpu), &
                    1e-10, "Array diss_dcc"//trim(flag_str))
   end if
   !if (myrank .eq. 3) write(*,*) 'kkc', th_cpu(235,:,57)
   !if (myrank .eq. 3) write(*,*) 'kkg', th_gpu(235,:,57)
   
   
   contains
   
      subroutine check_real2D(actual, desired, rtol, err_msg, cnt)
         implicit none

         real(kind=RTYPE), dimension(nxp, my_max), intent(in) :: actual
         real(kind=RTYPE), dimension(nxp, my_max), intent(in) :: desired
         real(kind=RTYPE), intent(in) :: rtol
         character(len=*), intent(in), optional :: err_msg
         real(kind=RTYPE), parameter :: eps = 1e-15
         real(kind=RTYPE) :: reldiff, absdiff, sum_actual, sum_desired, sum_rel, sum_abs
         logical :: equal
         integer :: i, cnt
         
         cnt = 0
         sum_actual = 0.
         sum_desired = 0.
         do jj = 1, jlistnum
            do i = 1, myim(jj)
               absdiff = abs(actual(i, jj) - desired(i, jj))
               reldiff = absdiff/desired(i, jj)
               sum_actual = sum_actual + actual(i, jj)
               sum_desired = sum_desired + desired(i, jj)
               if (reldiff .gt. rtol) then
                  !write(*,*) err_msg, myrank, jj, i, reldiff, absdiff, actual(i, jj), desired(i, jj)
                  cnt = cnt + 1
               end if
            end do
         end do
         sum_abs = abs(sum_actual - sum_desired)
         sum_rel = sum_abs/sum_desired
         write(*,*) err_msg, myrank, cnt, sum_abs, sum_rel, sum_actual, sum_desired


      end subroutine check_real2D
      
      subroutine check_real3D(actual, desired, rtol, err_msg, cnt)
         implicit none

         real(kind=RTYPE), dimension(nxp, lev, my_max), intent(in) :: actual
         real(kind=RTYPE), dimension(nxp, lev, my_max), intent(in) :: desired
         real(kind=RTYPE), intent(in) :: rtol
         character(len=*), intent(in), optional :: err_msg
         real(kind=RTYPE), parameter :: eps = 1e-15
         real(kind=RTYPE) :: reldiff, absdiff, sum_actual, sum_desired, sum_rel, sum_abs
         logical :: equal
         integer :: i, cnt
         
         cnt = 0
         sum_actual = 0.
         sum_desired = 0.
         do jj = 1, jlistnum
            do i = 1, myim(jj)
               do k = 1, lev
                  absdiff = abs(actual(i, k, jj) - desired(i, k, jj))
                  reldiff = (absdiff)/desired(i, k, jj)
                  sum_actual = sum_actual + actual(i, k, jj)
                  sum_desired = sum_desired + desired(i, k, jj)
                  if (reldiff .gt. rtol) then
                     !write(*,*) err_msg, myrank, jj, i, k, reldiff, absdiff
                     cnt = cnt + 1
                  end if
               end do
            end do
         end do
         sum_abs = abs(sum_actual - sum_desired)
         sum_rel = sum_abs/sum_desired
         write(*,*) err_msg, myrank, cnt, sum_abs, sum_rel, sum_actual, sum_desired
         


      end subroutine check_real3D


end subroutine samfdeepcnv_kh_unit

program test_samfdeepcnv_kh
   use param, only: lev, my_max, my
   use const, only: naero
   use index, only: nxp, jlistnum
   use rank, only: myrank

   implicit none
   integer :: itimestep, i
   real :: ct, gt
   
   
   ct = 0.
   gt = 0.
   
   call mpe_init
   call cons
   call samfdeepcnv_kh_unit(4, .false., ct, gt)
   do i = 1, 15
      call samfdeepcnv_kh_unit(4, .true., ct, gt)
   end do
   
   write (*, *) 'Timing for CPU & GPU: ', ct, gt
   write (*, *) 'speedup ratio: ', ct/gt, myrank

   
   call mpe_finalize

end program

