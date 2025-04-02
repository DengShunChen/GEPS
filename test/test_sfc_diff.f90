program test_sfc_diff
   implicit none
   call mpe_init
   call cons
   call sfc_diff_unit
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

subroutine sfc_diff_unit
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
   real, dimension(nxp, my_max)  :: islmsk_r, ivegtyp_r, flag_iter_r
   integer, parameter :: nxpvs = 7501
   real :: c1xpvs,c2xpvs,tbpvs(nxpvs)
   common/fpvscom/ c1xpvs,c2xpvs,tbpvs(nxpvs)

!-------------- inputs
   integer ivegsrc
   integer, dimension(my_max) :: myim
   integer, dimension(nxp, my_max) :: islmsk, ivegtyp
   logical redrag
   logical, dimension(nxp, my_max) :: flag_iter
   real(kind=kind_phys), dimension(nxp, my_max) :: psi, hgt, snwdph, tg, &
      prsl1, prslki, ddvel, sigmaf, shdmax, tsurf
   real(kind=kind_phys), dimension(nxp, lev, my_max) :: ut, vt, tt
   real(kind=kind_phys), dimension(nxp, lev*ncld, my_max) :: qt
!-------------- inputs/outputs
   real(kind=kind_phys), dimension(nxp, my_max) :: z0rl, ustar
   real(kind=kind_phys), dimension(nxp, my_max) :: z0rl_gpu, ustar_gpu
!--------------outputs
   real(kind=kind_phys), dimension(nxp, my_max) :: cd, cdq, rb, stress, &
      fm, fh, sfcw, fm10, fh2, fh10
   real(kind=kind_phys), dimension(nxp, my_max) :: cd_gpu, cdq_gpu, rb_gpu, &
      stress_gpu, fm_gpu, fh_gpu, sfcw_gpu, fm10_gpu, fh2_gpu, &
      fh10_gpu
   
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
   call random_number(hgt)
   call random_number(snwdph)
   call random_number(tg)
   call random_number(z0rl)
   call random_number(prsl1)
   call random_number(prslki)
   call random_number(islmsk_r)
   call random_number(ustar)
   call random_number(ddvel)
   call random_number(sigmaf)
   call random_number(ivegtyp_r)
   call random_number(shdmax)
   call random_number(tsurf)
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
   ivegtyp = floor(ivegtyp_r*20)
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
   hgt = hgt*3. + 20.
   snwdph = snwdph*1000.
   tg = tg*30. + 260.
   z0rl = z0rl*100.
   prsl1 = prsl1*30000. + 70000.
   prslki = prslki*0.001 + 1.
   ustar = ustar*1.
   ddvel = ddvel*1.
   sigmaf = sigmaf*1.
   shdmax = shdmax*1.
   tsurf = tsurf*30. + 260.
   
   redrag = .false.
   ivegsrc = 1
   ntrac = ncld

   z0rl_gpu = z0rl
   ustar_gpu = ustar
   cd = 0.
   cdq = 0.
   rb = 0.
   stress = 0.
   fm = 0.
   fh = 0.
   sfcw = 0.
   fm10 = 0.
   fh2 = 0.
   fh10 = 0.     
   cd_gpu = 0.
   cdq_gpu = 0.
   rb_gpu = 0.
   stress_gpu = 0.
   fm_gpu = 0.
   fh_gpu = 0.
   sfcw_gpu = 0.
   fm10_gpu = 0.
   fh2_gpu = 0.
   fh10_gpu = 0.   
   
   do jj = 1, jlistnum
      j = jlist1(jj)
      myim(jj) = nxjp(j)
   end do

   do ii = 1, 16
      call cpu_time(time1)
      !if (ii .ne. 1) call nvtxStartRange("CPU compute")
      do jj = 1, jlistnum
         !write(*,*) 'CPU',myrank, j, nxjp(j)
         call sfc_diff(myim(jj), nx, psi(1, jj), ut(1, lev, jj), &
                       vt(1, lev, jj), tt(1, lev, jj), &
                       qt(1, lev, jj), &
                       hgt(1, jj), snwdph(1, jj), tg(1, jj), z0rl(1, jj), cd(1, jj), &
                       cdq(1, jj), rb(1, jj), prsl1(1, jj), prslki(1, jj), islmsk(1, jj), &
                       stress(1, jj), fm(1, jj), fh(1, jj), &
                       ustar(1, jj), sfcw(1, jj), ddvel(1, jj), fm10(1, jj), &
                       fh2(1, jj), fh10(1, jj), sigmaf(1, jj), &
                       ivegtyp(1, jj), shdmax(1, jj), ivegsrc, &
                       tsurf(1, jj), flag_iter(1, jj), redrag)
      end do
      call cpu_time(time2)
      !if (ii .ne. 1) call nvtxEndRange
      if (ii .ne. 1) ct = ct + time2 - time1
   end do

      !write(*,*) 12345

   if (.true.) then
      !$acc enter data copyin(tbpvs) async(async_id)
      !$acc enter data copyin(myim, psi, ut, vt, tt, qt, hgt, snwdph, &
      !$acc&      tg, prsl1, prslki, islmsk, ddvel, sigmaf, ivegtyp, &
      !$acc&      shdmax, tsurf, flag_iter) async(async_id)
      !$acc enter data copyin(z0rl_gpu, cd_gpu, cdq_gpu, rb_gpu, stress_gpu, &
      !$acc&      fm_gpu, fh_gpu, ustar_gpu, &
      !$acc&      sfcw_gpu, fm10_gpu, fh2_gpu, fh10_gpu) async(async_id)
      !$acc wait(async_id)
      do ii = 1, 16
         call cpu_time(time1)
         !if (ii .ne. 1) call nvtxStartRange("GPU compute")
         call sfc_diff_gpu(myim, nx, lev, ncld, psi, ut, &
                       vt, tt, qt, &
                       hgt, snwdph, tg, z0rl_gpu, cd_gpu, &
                       cdq_gpu, rb_gpu, prsl1, prslki, islmsk, &
                       stress_gpu, fm_gpu, fh_gpu, &
                       ustar_gpu, sfcw_gpu, ddvel, fm10_gpu, &
                       fh2_gpu, fh10_gpu, sigmaf, &
                       ivegtyp, shdmax, ivegsrc, &
                       tsurf, flag_iter, redrag, async_id)
         !$acc wait(async_id)
         !if (ii .ne. 1) call nvtxEndRange
         call cpu_time(time2)
         if (ii .ne. 1) gt = gt + time2 - time1
      end do
      !$acc exit data delete(myim, psi, ut, vt, tt, qt, hgt, snwdph, &
      !$acc&     tg, prsl1, prslki, islmsk, ddvel, sigmaf, ivegtyp, &
      !$acc&     shdmax, tsurf, flag_iter) async(async_id)
      !$acc exit data copyout(z0rl_gpu, cd_gpu, cdq_gpu, rb_gpu, stress_gpu, &
      !$acc&     fm_gpu, fh_gpu, ustar_gpu, &
      !$acc&     sfcw_gpu, fm10_gpu, fh2_gpu, fh10_gpu) async(async_id)
      !$acc exit data delete(tbpvs) async(async_id)
      !$acc wait(async_id)
   end if

   !write(*,*) 'NaN Check:', sum(u1), sum(v1), sum(t1), sum(q1),sum(kpbl),sum(hpbl)
   !write(*,*) 'u1 2', u1(1,1,2), u1_gpu(1,1,2)

   call assert_real(cd_gpu, size(cd_gpu), cd, size(cd), &
                    1e-12, "Array cd")
   call assert_real(cdq_gpu, size(cdq_gpu), cdq, size(cdq), &
                    1e-12, "Array cdq")
   call assert_real(rb_gpu, size(rb_gpu), rb, size(rb), &
                    1e-12, "Array rb")
   call assert_real(stress_gpu, size(stress_gpu), stress, size(stress), &
                    1e-12, "Array stress")
   call assert_real(fm_gpu, size(fm_gpu), fm, size(fm), &
                    1e-12, "Array fm")
   call assert_real(fh_gpu, size(fh_gpu), fh, size(fh), &
                    1e-12, "Array fh")
   call assert_real(ustar_gpu, size(ustar_gpu), ustar, size(ustar), &
                    1e-12, "Array ustar")
   call assert_real(sfcw_gpu, size(sfcw_gpu), sfcw, size(sfcw), &
                    1e-12, "Array sfcw")
   call assert_real(fm10_gpu, size(fm10_gpu), fm10, size(fm10), &
                    1e-12, "Array fm10")
   call assert_real(fh2_gpu, size(fh2_gpu), fh2, size(fh2), &
                    1e-12, "Array fh2")
   call assert_real(fh10_gpu, size(fh10_gpu), fh10, size(fh10), &
                    1e-12, "Array fh10")
   call assert_real(z0rl_gpu, size(z0rl_gpu), z0rl, size(z0rl), &
                    1e-12, "Array z0rl")

   write (*, *) 'Timing for CPU & GPU: ', ct, gt
   write (*, *) 'speedup ratio: ', ct/gt, myrank

end subroutine sfc_diff_unit

