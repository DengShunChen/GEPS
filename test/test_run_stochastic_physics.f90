program test_run_stochastic_physics
   use rank, only: myrank
   use const, only: use_zmtnblck
   use mod_stochastic_physics, only : init_stochastic_physics, &
      destroy_stochastic_physics, sppt3d, shum3d, skeb3du, skeb3dv, ssst3d, &
      ncep_seeds, sppt, sppt_seed, sppt_decort, sppt_lscale, &
      sppt_sigtop1, sppt_sigtop2, sppt_sigbot1, sppt_sigbot2, &
      sppt_sfclimit, sppt_logit, &
      shum, shum_seed, shum_decort, shum_lscale, &
      shum_sigefold, &
      skeb, skeb_seed, skeb_decort, skeb_lscale, &
      skeb_sigtop1, skeb_sigtop2, skeb_sigbot1, skeb_sigbot2, &
      skeb_vdof,skebnorm, skebfilt, &
      ssst, ssst_seed, ssst_decort, ssst_lscale, do_unit_test
   use mod_stochastic_physics_gpu, only : init_stochastic_physics_gpu, &
      destroy_stochastic_physics_gpu, sppt3d_gpu => sppt3d, shum3d_gpu => shum3d, &
      skeb3du_gpu => skeb3du, skeb3dv_gpu => skeb3dv, ssst3d_gpu => ssst3d, &
      ncep_seeds_gpu => ncep_seeds, &
      sppt_gpu => sppt, &
      sppt_seed_gpu => sppt_seed, &
      sppt_decort_gpu => sppt_decort, &
      sppt_lscale_gpu => sppt_lscale, &
      sppt_sigtop1_gpu => sppt_sigtop1, &
      sppt_sigtop2_gpu => sppt_sigtop2, &
      sppt_sigbot1_gpu => sppt_sigbot1, &
      sppt_sigbot2_gpu => sppt_sigbot2, &
      sppt_sfclimit_gpu => sppt_sfclimit, &
      sppt_logit_gpu => sppt_logit, &
      shum_gpu => shum, &
      shum_seed_gpu => shum_seed, &
      shum_decort_gpu => shum_decort, &
      shum_lscale_gpu => shum_lscale, &
      shum_sigefold_gpu => shum_sigefold, &
      skeb_gpu => skeb, &
      skeb_seed_gpu => skeb_seed, &
      skeb_decort_gpu => skeb_decort, &
      skeb_lscale_gpu => skeb_lscale, &
      skeb_sigtop1_gpu => skeb_sigtop1, &
      skeb_sigtop2_gpu => skeb_sigtop2, &
      skeb_sigbot1_gpu => skeb_sigbot1, &
      skeb_sigbot2_gpu => skeb_sigbot2, &
      skeb_vdof_gpu => skeb_vdof, &
      skebnorm_gpu => skebnorm, &
      skebfilt_gpu => skebfilt, &
      ssst_gpu => ssst, &
      ssst_seed_gpu => ssst_seed, &
      ssst_decort_gpu => ssst_decort, &
      ssst_lscale_gpu => ssst_lscale, &
      do_unit_test_gpu => do_unit_test
      
   implicit none
   do_unit_test = .true.
   do_unit_test_gpu = .true.
   call mpe_init
   call cons
      ncep_seeds = ncep_seeds_gpu
      sppt = sppt_gpu
      sppt_seed = sppt_seed_gpu
      sppt_decort = sppt_decort_gpu
      sppt_lscale = sppt_lscale_gpu
      sppt_sigtop1 = sppt_sigtop1_gpu
      sppt_sigtop2 = sppt_sigtop2_gpu
      sppt_sigbot1 = sppt_sigbot1_gpu
      sppt_sigbot2 = sppt_sigbot2_gpu
      sppt_sfclimit = sppt_sfclimit_gpu
      sppt_logit = sppt_logit_gpu
      shum = shum_gpu
      shum_seed = shum_seed_gpu
      shum_decort = shum_decort_gpu
      shum_lscale = shum_lscale_gpu
      shum_sigefold = shum_sigefold_gpu
      skeb = skeb_gpu
      skeb_seed = skeb_seed_gpu
      skeb_decort = skeb_decort_gpu
      skeb_lscale = skeb_lscale_gpu
      skeb_sigtop1 = skeb_sigtop1_gpu
      skeb_sigtop2 = skeb_sigtop2_gpu
      skeb_sigbot1 = skeb_sigbot1_gpu
      skeb_sigbot2 = skeb_sigbot2_gpu
      skeb_vdof = skeb_vdof_gpu
      skebnorm = skebnorm_gpu
      skebfilt = skebfilt_gpu
      ssst = ssst_gpu
      ssst_seed = ssst_seed_gpu
      ssst_decort = ssst_decort_gpu
      ssst_lscale = ssst_lscale_gpu
   call init_stochastic_physics(720.)
   call run_stochastic_physics_unit
   call destroy_stochastic_physics()
   call destroy_stochastic_physics_gpu()
   call mpe_finalize
   

end program

subroutine assert_real(actual, n_actual, desired, n_desired, rtol, err_msg)
   use const, only: RTYPE
   use rank, only: myrank

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
      print *, "Max relative difference = ", rel_diff, myrank
      !print *, "Max absolute difference = ", abs_diff
      if (present(err_msg)) print *, err_msg
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

subroutine run_stochastic_physics_unit
   use mod_stochastic_physics, only : run_stochastic_physics, &
      destroy_stochastic_physics, sppt3d, shum3d, skeb3du, skeb3dv, ssst3d, &
      rpattern_sppt, rpattern_shum, rpattern_skeb, rpattern_ssst, &
      nsppt, nshum, nskeb, nssst, skeblevs, skebest
   use mod_stochastic_physics_gpu, only : run_stochastic_physics_gpu, skebest_gpu, &
      destroy_stochastic_physics_gpu, sppt3d_gpu => sppt3d, shum3d_gpu => shum3d, &
      !skeb3du_gpu => skeb3du, skeb3dv_gpu => skeb3dv, ssst3d_gpu => ssst3d
      skeb3du_gpu => skeb3du, skeb3dv_gpu => skeb3dv, ssst3d_gpu => ssst3d, &
      rpattern_skeb_n2du_gpu => rpattern_skeb_n2du, rpattern_skeb_n2dv_gpu => rpattern_skeb_n2dv
   use radn, only: me
   use index, only: nxjp, nxp, jlistnum, jlist1, nxdef, nlist, tcolt_jlist, &
      poly_mlist, jlist2, mtrundef, nxjlen, nxjstart, mlist, nxdef_2d
   use param, only: lev, my_max, my
   use const, only: grav, cp, nx, mtnvar, tofd, cmbk, cgwd, polyf, RTYPE, onocos, &
      wcfac, wdfac, poly, dpoly, weight, cosl, first_call
   use physcons, only: con_rd, con_rv
   use rank, only: myrank
   use fftcom
   !use nvtx
   use grid, only: ut, vt

   implicit none


   integer :: i, j, jj, async_id, n, k, cnt, myim(my_max)
   integer, dimension(34) :: seed

   real :: time1, time2, ct, gt
   real(kind=RTYPE), allocatable, dimension(:,:,:,:) :: rpattern_skeb_n2du_cpu, rpattern_skeb_n2dv_cpu
   real(kind=RTYPE), dimension(nxp, lev, my_max) :: ut_ori, vt_ori, ut_cpu, vt_cpu, &
      ut_gpu, vt_gpu, um, vm

   do jj = 1, jlistnum
      j = jlist1(jj)
      myim(jj) = nxjp(j)
   end do

   async_id = 1
   seed = (/10004, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,&
           &25, 26, 27, 28, 29, 30, 31, 32, 33/)
   !call random_seed()
   call random_seed(put=seed)
   call random_number(um)
   call random_number(vm)
   call random_number(ut)
   call random_number(vt)
   um = um*2e-7
   vm = vm*2e-7
   ut = ut*2e-5
   vt = vt*2e-5
   ut_ori = ut
   vt_ori = vt

   ct = 0.
   do i = 1, 10
      call cpu_time(time1)
      !if (i .gt. 2) call nvtxStartRange("CPU compute")
      call run_stochastic_physics()
      call skebest(um, vm)
      !if (i .gt. 2) call nvtxEndRange
      call cpu_time(time2)
      if (i .gt. 2) ct = ct + time2 - time1
   end do
   ut_cpu = ut
   vt_cpu = vt
   ut = ut_ori
   vt = vt_ori

   if (.true.) then
      gt = 0.
      !$acc enter data copyin(polyf, trigsj, ifaxj, jlist1, nxdef, onocos, &
      !$acc&      wcfac, wdfac, poly, dpoly, nlist, tcolt_jlist, poly_mlist, &
      !$acc&      jlist2, mtrundef, nxjlen, nxjstart, weight, mlist, nxdef_2d, cosl) async(async_id)
      !$acc enter data copyin(ut, vt, um, vm) async(async_id)
      !$acc wait(async_id)
      do i = 1, 10
         call cpu_time(time1)
         !if (i .gt. 2) call nvtxStartRange("GPU compute")
         call run_stochastic_physics_gpu()
         !$acc wait(async_id)
         !if (i .gt. 2) call nvtxEndRange
         !if (i .gt. 2) call nvtxStartRange("GPU compute")
         call skebest_gpu(um, vm)
         !$acc wait(async_id)
         !if (i .gt. 2) call nvtxEndRange
         call cpu_time(time2)
         if (i .gt. 2) gt = gt + time2 - time1
      end do
      !$acc exit data delete(polyf, trigsj, ifaxj, jlist1, nxdef, onocos, &
      !$acc&      wcfac, wdfac, poly, dpoly, nlist, tcolt_jlist, poly_mlist, cosl, &
      !$acc&      jlist2, mtrundef, nxjlen, nxjstart, weight, um, vm, mlist, nxdef_2d) async(async_id)
      !$acc update self(sppt3d_gpu, shum3d_gpu, ssst3d_gpu, skeb3du_gpu, skeb3dv_gpu, &
      !$acc&       rpattern_skeb_n2du_gpu, rpattern_skeb_n2dv_gpu) async(async_id)
      !$acc exit data copyout(ut, vt) async(async_id)
      !$acc wait(async_id)
      allocate(rpattern_skeb_n2du_cpu(nxp, my_max, skeblevs, nskeb))
      allocate(rpattern_skeb_n2dv_cpu(nxp, my_max, skeblevs, nskeb))
      do n = 1, nskeb
         rpattern_skeb_n2du_cpu(:,:,:,n) = rpattern_skeb(n)%n2du(:,:,:)
         rpattern_skeb_n2dv_cpu(:,:,:,n) = rpattern_skeb(n)%n2dv(:,:,:)
      end do
      ut_gpu = ut
      vt_gpu = vt

      !call check_real3D(rpattern_skeb_n2dv_gpu, rpattern_skeb_n2dv_cpu, 1e-2, 'n2du', cnt)
      !write(*,*) cnt
      call assert_real(sppt3d_gpu, size(sppt3d_gpu), sppt3d, size(sppt3d), &
                       1e-8, "Array sppt3d")
      call assert_real(shum3d_gpu, size(shum3d_gpu), shum3d, size(shum3d), &
                       1e-8, "Array shum3d")
      call assert_real(skeb3du_gpu, size(skeb3du_gpu), skeb3du, size(skeb3du), &
                       1e-12, "Array skeb3du")
      call assert_real(skeb3dv_gpu, size(skeb3dv_gpu), skeb3dv, size(skeb3dv), &
                       1e-12, "Array skeb3dv")
      call assert_real(rpattern_skeb_n2du_gpu, size(rpattern_skeb_n2du_gpu), &
                       rpattern_skeb_n2du_cpu, size(rpattern_skeb_n2du_cpu), &
                       1e-12, "Array rpattern_skeb_n2du")
      call assert_real(rpattern_skeb_n2dv_gpu, size(rpattern_skeb_n2dv_gpu), &
                       rpattern_skeb_n2dv_cpu, size(rpattern_skeb_n2dv_cpu), &
                       1e-12, "Array rpattern_skeb_n2dv")
      call assert_real(ssst3d_gpu, size(ssst3d_gpu), ssst3d, size(ssst3d), &
                       1e-12, "Array ssst3d")
      call assert_real(ut_gpu, size(ut_gpu), ut_cpu, size(ut_cpu), &
                       1e-12, "Array ut")
      call assert_real(vt_gpu, size(vt_gpu), vt_cpu, size(vt_cpu), &
                       1e-12, "Array vt")
      deallocate(rpattern_skeb_n2du_cpu, rpattern_skeb_n2dv_cpu)
   end if
   write (*, *) 'Timing for CPU & GPU: ', ct, gt
   write (*, *) 'speedup ratio: ', ct/gt

   contains
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
                  reldiff = abs((absdiff)/desired(i, k, jj))
                  sum_actual = sum_actual + actual(i, k, jj)
                  sum_desired = sum_desired + desired(i, k, jj)
                  if (reldiff .gt. rtol) then
                     write(*,*) err_msg, myrank, jj, i, k, reldiff, absdiff, actual(i, k, jj), desired(i, k, jj)
                     cnt = cnt + 1
                  end if
               end do
            end do
         end do
         !sum_abs = abs(sum_actual - sum_desired)
         !sum_rel = sum_abs/sum_desired
         !write(*,*) err_msg, myrank, cnt, sum_abs, sum_rel, sum_actual, sum_desired
         


      end subroutine check_real3D

end subroutine run_stochastic_physics_unit
