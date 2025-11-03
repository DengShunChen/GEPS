
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

subroutine ozphys_2015_unit(itimestep, benchmark, ct, gt)
   use index, only: jlistnum, jlist1, nxdef_2d, nxjp_acc, nxjp, nxp
   use param, only: my_max, my, ncld, lev
   use machine, only: kind_phys
   use rank, only: myrank
   use const, only: RTYPE, naero, sigma, dsigma, julian
   !use nvtx
   use ozne_def
   use phygrid, only: xlat, o3l
   use grid, only: tt, pk, pk2, pt

   implicit none

   integer, parameter :: icheck = 68, kcheck = 71, jjcheck = 95, rcheck = 3
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
   real(kind=RTYPE), dimension(nxp, my_max) :: pst
   real,dimension(nxp, lev, my_max) :: plt
   real, dimension(nxp, lev, my_max) :: del
   real, dimension(nxp, lev, my_max) :: o3l_gpu, o3l_cpu
   
   
   write(flag_str,'(I3)') itimestep
   async_id = 1
   seed = (/10006, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,&
           &25, 26, 27, 28, 29, 30, 31, 32, 33/)
   !call random_seed()
   call random_seed(put=seed)
   cnt = 0
   tcnt = 0
   
   dta = 600.
   ptop = 0.1000000000000000
   pst = pt
   
   do jj = 1, jlistnum
      j = jlist1(jj)
      myim(jj) = nxjp(j)
   end do
   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef_2d(j)
      !
      !   new p**capa quantities were computed in previous diabat call
      !
      call prexp_hybrid_cwb(nxjp(j), nxp, lev, ptop, sigma, pt(1, jj), &
                          pk(1, 1, jj), pk2(1, 1, jj), plt(1, 1, jj))
   end do !jj = 1,jlistnum
   !write(*,*) "test itimestep = ", itimestep

   o3l_gpu = o3l
   o3l_cpu = o3l



   !if (benchmark .eq. .true.) call nvtxStartRange("CPU compute")
      call cpu_time(time1)
      do jj = 1, jlistnum
         j = jlist1(jj)
         nxj = nxdef_2d(j)
         do k=1,lev
            do i=1,nxj
               del(i,k,jj) = 100.0*( dsigma(k,1)*pst(i,jj)+dsigma(k,2))  !  pa
            enddo
         enddo
         call ozphys_2015 (nxp, nxjp(j), lev , dta, xlat(j), julian, &
                      o3l_cpu(1,1,jj), o3l_cpu(1,1,jj), tt(1,1,jj),          &
                      plt(1,1,jj), del(1,1,jj), myrank, jj, rcheck, jjcheck, kcheck, icheck)
      end do
      call cpu_time(time2)
      if (benchmark .eq. .true.) ct = ct + time2 - time1
   !if (benchmark .eq. .true.) call nvtxEndRange


   if (.false.) then
         !$acc wait(async_id)
         !$acc enter data copyin(jlist1, nxdef_2d, nxjp, xlat, tt, plt, dsigma, pst) &
         !$acc&      async(async_id)
         !$acc enter data copyin(o3l_gpu) async(async_id)
         !$acc enter data copyin(ozplin, pl_lat, pl_pres) async(async_id)
         call cpu_time(time1)
         !if (benchmark .eq. .true.) call nvtxStartRange("GPU compute")
         call ozphys_2015_gpu (nxp, nxjp, lev , dta, xlat, julian, &
                      o3l_gpu, o3l_gpu, tt,          &
                      plt, dsigma, pst, myrank)
         !$acc wait(async_id)
         !if (benchmark .eq. .true.) call nvtxEndRange
         call cpu_time(time2)
         if (benchmark .eq. .true.) gt = gt + time2 - time1
         !$acc wait(async_id)
         !$acc exit data delete(jlist1, nxdef_2d, nxjp, xlat, tt, plt, dsigma, pst) &
         !$acc&      async(async_id)
         !$acc exit data copyout(o3l_gpu) async(async_id)
         !$acc exit data delete(ozplin, pl_lat, pl_pres) async(async_id)
         
   call assert_real(o3l_gpu, size(o3l_gpu), o3l_cpu, size(o3l_cpu), &
                    1e-10, "Array o3l"//trim(flag_str))
   end if
   
   
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
                  write(*,*) err_msg, myrank, jj, i, reldiff, absdiff, actual(i, jj), desired(i, jj)
                  cnt = cnt + 1
               end if
            end do
         end do
         sum_abs = abs(sum_actual - sum_desired)
         sum_rel = sum_abs/sum_desired
         !write(*,*) err_msg, myrank, sum_abs, sum_rel, sum_actual, sum_desired
         write(*,*) err_msg, myrank, cnt


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
         !write(*,*) err_msg, myrank, sum_abs, sum_rel, sum_actual, sum_desired
         


      end subroutine check_real3D


end subroutine ozphys_2015_unit

program test_ozphys_2015
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
   call getrdy
   !call getrdy
   if (donnmi .and. taui .lt. 1.0) then
      no = 2*((jtrun + 1)/2) + (jtrun/2) + 10
      call initial(no, jtrun, jtmax, lev, nx, my, my_max, mlmax)
   end if
   
   call ozphys_2015_unit(3, .false., ct, gt)
   do i = 1, 15
      call ozphys_2015_unit(3, .true., ct, gt)
   end do
   
   write (*, *) 'Timing for CPU & GPU: ', ct, gt
   write (*, *) 'speedup ratio: ', ct/gt, myrank

   
   call mpe_finalize

end program

