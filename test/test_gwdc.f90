
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

subroutine gwdc_unit(itimestep, benchmark, ct, gt)
   use index, only: jlistnum, jlist1, nxdef_2d, nxjp_acc, nxjp, nxp
   use param, only: my_max, my, ncld, lev
   use machine, only: kind_phys
   use rank, only: myrank
   use const, only: dsigma, RTYPE, naero
   use physcons, only: con_rd, con_fvirt
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
   real :: dta, grav, cp
   real, dimension(nxp, lev, my_max) :: u0, v0, t0, prsl, del
   real, dimension(nxp, lev+1, my_max) :: prsi
   real(kind=RTYPE), dimension(nxp, lev*ncld, my_max) :: q0
   integer, dimension(nxp, my_max) :: ktop, kbot, kuo
   real, dimension(nxp, my_max) :: ktop_r, kbot_r, kuo_r
   real, dimension(nxp, my_max) :: cldf, cumabs, dlength
   
   
   real, dimension(nxp, lev, my_max) :: utgwc_gpu, vtgwc_gpu
   real, dimension(nxp, my_max) :: tauctx_gpu, taucty_gpu
   real, dimension(nxp, lev, my_max) :: utgwc_cpu, vtgwc_cpu
   real, dimension(nxp, my_max) :: tauctx_cpu, taucty_cpu
   
   
   write(flag_str,'(I3)') itimestep
   async_id = 1
   seed = (/10006, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,&
           &25, 26, 27, 28, 29, 30, 31, 32, 33/)
   !call random_seed()
   call random_seed(put=seed)
   call random_number(u0)
   call random_number(v0)
   call random_number(t0)
   call random_number(q0)
   call random_number(prsl)
   call random_number(prsi)
   call random_number(del)
   call random_number(ktop_r)
   call random_number(kbot_r)
   call random_number(kuo_r)
   call random_number(cldf)
   call random_number(cumabs)
   call random_number(dlength)
   ktop_r = floor(ktop_r*30) + 1
   kbot_r = 35 + floor(kbot_r*30) + 1
   do j = 1, jlistnum
      do i = 1, nxp
         if (kuo_r(i, jj) .gt. 0.5) then
            kuo(i, jj) = 1
         else
            kuo(i, jj) = 0
         end if
      end do
   end do
   
   u0 = u0*60. - 30.
   v0 = v0*40. - 20.
   t0 = t0*100. + 200.
   q0 = q0*1.8e-2 + 1.e-20
   prsl = prsl*100000. + 100.
   prsi = prsi*100000. + 100.
   del = del*3000. + 60.
   cldf = cldf*0.002 + 0.335
   cumabs = cumabs*5e-4
   dlength = dlength*3000. + 39700
   dta = 600.
   grav = 9.806649999999999
   cp = 1004.600000000000

   
   
   cnt = 0
   tcnt = 0
   
   do jj = 1, jlistnum
      j = jlist1(jj)
      myim(jj) = nxjp(j)
   end do




   utgwc_gpu = 0.
   vtgwc_gpu = 0.
   tauctx_gpu = 0.
   taucty_gpu = 0.

   utgwc_cpu = 0.
   vtgwc_cpu = 0.
   tauctx_cpu = 0.
   taucty_cpu = 0.

   !if (benchmark .eq. .true.) call nvtxStartRange("CPU compute")
      call cpu_time(time1)
         do jj = 1, jlistnum
            j = jlist1(jj)
            nxj = nxdef_2d(j)
            call gwdc(nxjp(j), nxp, nxp, lev, u0(1, 1, jj), v0(1, 1, jj), &
                      t0(1, 1, jj), q0(1, 1, jj), prsl(1, 1, jj), prsi(1, 1, jj), del(1, 1, jj), &
                      ktop(1, jj), kbot(1, jj), kuo(1, jj), cldf(1, jj), cumabs(1, jj), &
                      grav, cp, con_rd, con_fvirt, dta, dlength(1, jj), &
                      utgwc_cpu(1, 1, jj), vtgwc_cpu(1, 1, jj), tauctx_cpu(1, jj), taucty_cpu(1, jj), j)
         end do
      call cpu_time(time2)
      if (benchmark .eq. .true.) ct = ct + time2 - time1
   !if (benchmark .eq. .true.) call nvtxEndRange


   if (.true.) then
         !$acc wait(async_id)
         !$acc enter data copyin(jlist1, nxjp) async(async_id)
         !$acc enter data copyin(u0, v0, t0, q0, prsl, prsi, del, ktop, &
         !$acc&      kbot, kuo, cldf, cumabs, dlength) async(async_id)
         !$acc enter data copyin(utgwc_gpu, vtgwc_gpu, tauctx_gpu, taucty_gpu) &
         !$acc&      async(async_id)
         call cpu_time(time1)
         !if (benchmark .eq. .true.) call nvtxStartRange("GPU compute")
         call gwdc_gpu(nxp, nxp, lev, u0, v0, &
                   t0, q0, prsl, prsi, del, &
                   ktop, kbot, kuo, cldf, cumabs, &
                   grav, cp, con_rd, con_fvirt, dta, dlength, &
                   utgwc_gpu, vtgwc_gpu, tauctx_gpu, taucty_gpu)
         !$acc wait(async_id)
         !if (benchmark .eq. .true.) call nvtxEndRange
         call cpu_time(time2)
         if (benchmark .eq. .true.) gt = gt + time2 - time1
         !$acc exit data delete(jlist1, nxjp) async(async_id)
         !$acc exit data delete(u0, v0, t0, q0, prsl, prsi, del, ktop, &
         !$acc&     kbot, kuo, cldf, cumabs, dlength) async(async_id)
         !$acc exit data copyout(utgwc_gpu, vtgwc_gpu, tauctx_gpu, taucty_gpu) &
         !$acc&     async(async_id)
         !$acc wait(async_id)
      
   call assert_real(utgwc_gpu, size(utgwc_gpu), utgwc_cpu, size(utgwc_cpu), &
                    1e-10, "Array utgwc"//trim(flag_str))
   call assert_real(vtgwc_gpu, size(vtgwc_gpu), vtgwc_cpu, size(vtgwc_cpu), &
                    1e-10, "Array vtgwc"//trim(flag_str))
   call assert_real(tauctx_gpu, size(tauctx_gpu), tauctx_cpu, size(tauctx_cpu), &
                    1e-10, "Array tauctx"//trim(flag_str))
   call assert_real(taucty_gpu, size(taucty_gpu), taucty_cpu, size(taucty_cpu), &
                    1e-10, "Array taucty"//trim(flag_str))
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
         write(*,*) err_msg, myrank, sum_abs, sum_rel, sum_actual, sum_desired
         


      end subroutine check_real3D


end subroutine gwdc_unit

program test_gwdc
   use param, only: lev, my_max, my, ncld
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
   call gwdc_unit(3, .false., ct, gt)
   do i = 1, 15
      call gwdc_unit(3, .true., ct, gt)
   end do
   
   write (*, *) 'Timing for CPU & GPU: ', ct, gt
   write (*, *) 'speedup ratio: ', ct/gt, myrank
   call mpe_finalize
   

end program

