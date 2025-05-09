program test_sfc_ocean
   implicit none
   call mpe_init
   call cons
   call sfc_ocean_unit
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

subroutine sfc_ocean_unit
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
   integer, dimension(34) :: seed

   real :: time1, time2, ct, gt, diff
   real, dimension(nxp, my_max)  :: islmsk_r, flag_iter_r
   integer, parameter :: nxpvs = 7501
   real :: c1xpvs,c2xpvs,tbpvs(nxpvs)
   common/fpvscom/ c1xpvs,c2xpvs,tbpvs(nxpvs)

!-------------- inputs
   integer, dimension(my_max) :: myim
   integer, dimension(nxp, my_max) :: islmsk
   logical, dimension(nxp, my_max) :: flag_iter
   real(kind=kind_phys), dimension(nxp, my_max) :: psi, tg, cd, cdq, &
                                                   prsl1, prslki, ddvel
   real(kind=kind_phys), dimension(nxp, lev, my_max) :: ut, vt, tt
   real(kind=kind_phys), dimension(nxp, lev*ncld, my_max) :: qt
!-------------- inputs/outputs
!--------------outputs
   real(kind=kind_phys), dimension(nxp, my_max) :: qsurf, gfx, qflux, &
                                                   hflux, ep1d
   real(kind=kind_phys), dimension(nxp, my_max) :: qsurf_gpu, gfx_gpu, &
                                                   qflux_gpu, hflux_gpu, ep1d_gpu
   
   async_id = 1
   seed = (/10006, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,&
           &25, 26, 27, 28, 29, 30, 31, 32, 33/)
   !call random_seed()
   call random_seed(put=seed)
   ct = 0.
   gt = 0.
   call random_number(psi)
   call random_number(ut)
   call random_number(vt)
   call random_number(tt)
   call random_number(qt)
   call random_number(tg)
   call random_number(cd)
   call random_number(cdq)
   call random_number(prsl1)
   call random_number(prslki)
   call random_number(islmsk_r)
   call random_number(ddvel)
   call random_number(flag_iter_r)

   do jj = 1, jlistnum
      do i = 1, nxp
         if (islmsk_r(i, jj) >= 0. .and. islmsk_r(i, jj) < 0.680) then
            islmsk(i, jj) = 0
         elseif (islmsk_r(i, jj) >= 0.680 .and. islmsk_r(i, jj) < 0.963) then
            islmsk(i, jj) = 1
         else
            islmsk(i, jj) = 2
         end if
      end do
   end do

   do jj = 1, jlistnum
      do i = 1, nxp
         if (flag_iter_r(i, jj) .gt. 0.) then
            flag_iter(i, jj) = .true.
         else
            flag_iter(i, jj) = .false.
         end if
      end do
   end do
   psi = psi*30000. + 70000.
   ut = ut*40. - 20.
   vt = vt*40. - 20.
   tt = tt*100. + 200.
   qt = qt*0.001
   tg = tg*30. + 260.
   cd = cd*0.05
   cdq = cdq*0.05
   prsl1 = prsl1*30000. + 70000.
   prslki = prslki*0.001 + 1.
   ddvel = ddvel*1.
   
   ntrac = ncld

   qsurf = 0.
   gfx = 0.
   qflux = 0.
   hflux = 0.
   ep1d = 0.
   qsurf_gpu = 0.
   gfx_gpu = 0.
   qflux_gpu = 0.
   hflux_gpu = 0.
   ep1d_gpu = 0.
   
   do jj = 1, jlistnum
      j = jlist1(jj)
      myim(jj) = nxjp(j)
   end do

   do ii = 1, 16
      call cpu_time(time1)
      !if (ii .ne. 1) call nvtxStartRange("CPU compute")
      do jj = 1, jlistnum
         !write(*,*) 'CPU',myrank, j, nxjp(j)
         call sfc_ocean(myim(jj), nx, psi(1, jj), ut(1, lev, jj), &
                        vt(1, lev, jj), tt(1, lev, jj), &
                        qt(1, lev, jj), tg(1, jj), cd(1, jj), cdq(1, jj), &
                        prsl1(1, jj), prslki(1, jj), islmsk(1, jj), &
                        ddvel(1, jj), flag_iter(1, jj), &
                        qsurf(1, jj), gfx(1, jj), qflux(1, jj), &
                        hflux(1, jj), ep1d(1, jj))
      end do
      call cpu_time(time2)
      !if (ii .ne. 1) call nvtxEndRange
      if (ii .ne. 1) ct = ct + time2 - time1
   end do

   !write(*,*) 12345

   if (.true.) then
      !$acc enter data copyin(tbpvs) async(async_id)
      !$acc enter data copyin(psi, ut, vt, tt, qt, tg, cd, cdq, prsl1, prslki, &
      !$acc&      islmsk, ddvel, flag_iter, myim) async(async_id)
      !$acc enter data copyin(qsurf_gpu, gfx_gpu, qflux_gpu, &
      !$acc&      hflux_gpu, ep1d_gpu) async(async_id)
      do ii = 1, 16
      !$acc wait(async_id)
         call cpu_time(time1)
         !if (ii .ne. 1) call nvtxStartRange("GPU compute")
         call sfc_ocean_gpu(myim, nx, lev, ncld, psi, ut, &
                        vt, tt, qt, tg, cd, cdq, &
                        prsl1, prslki, islmsk, &
                        ddvel, flag_iter, &
                        qsurf_gpu, gfx_gpu, qflux_gpu, &
                        hflux_gpu, ep1d_gpu, async_id)
         !if (ii .ne. 1) call nvtxEndRange
      !$acc wait(async_id)
         call cpu_time(time2)
         if (ii .ne. 1) gt = gt + time2 - time1
      end do
      !$acc exit data delete(tbpvs) async(async_id)
      !$acc exit data delete(psi, ut, vt, tt, qt, tg, cd, cdq, prsl1, prslki, &
      !$acc&     islmsk, ddvel, flag_iter) async(async_id)
      !$acc exit data copyout(qsurf_gpu, gfx_gpu, qflux_gpu, &
      !$acc&     hflux_gpu, ep1d_gpu) async(async_id)
      !$acc wait(async_id)



      call assert_real(qsurf_gpu, size(qsurf_gpu), qsurf, size(qsurf), &
                       1e-12, "Array qsurf")
      call assert_real(gfx_gpu, size(gfx_gpu), gfx, size(gfx), &
                       1e-12, "Array gfx")
      call assert_real(qflux_gpu, size(qflux_gpu), qflux, size(qflux), &
                       1e-12, "Array qflux")
      call assert_real(hflux_gpu, size(hflux_gpu), hflux, size(hflux), &
                       1e-12, "Array hflux")
      call assert_real(ep1d_gpu, size(ep1d_gpu), ep1d, size(ep1d), &
                       1e-12, "Array ep1d")
   end if

   write (*, *) 'Timing for CPU & GPU: ', ct, gt
   write (*, *) 'speedup ratio: ', ct/gt, myrank

end subroutine sfc_ocean_unit

