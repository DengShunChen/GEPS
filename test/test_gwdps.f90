program test_gwdps
   implicit none
   call mpe_init
   call cons
   call gwdps_unit
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

subroutine gwdps_unit
   use radn, only: me
   use index, only: nxjp, nxp, jlistnum, jlist1
   use param, only: lev, my_max, my
   use const, only: grav, cp, nx, mtnvar, tofd, cmbk, cgwd
   use physcons, only: con_rd, con_rv
   !use nvtx

   implicit none

   real :: dvdtc(nxp, lev, my_max), dudtc(nxp, lev, my_max), dtdtc(nxp, lev, my_max), &
           utc(nxp, lev, my_max), vtc(nxp, lev, my_max), ttc(nxp, lev, my_max), &
           qtc(nxp, lev, my_max), prsi(nxp, lev + 1, my_max), del(nxp, lev, my_max), &
           prsl(nxp, lev, my_max), prslk(nxp, lev, my_max), phii(nxp, lev + 1, my_max), &
           phil(nxp, lev, my_max)
   real :: dvdtc_gpu(nxp, lev, my_max), dudtc_gpu(nxp, lev, my_max), dtdtc_gpu(nxp, lev, my_max)
   integer :: kpbl(nxp, my_max)
   real :: dta
   integer :: kdt
   real :: hprime(nxp, my_max), oc(nxp, my_max), oa4(nxp, 4, my_max), clx(nxp, 4, my_max)
   real :: theta(nxp, my_max), sigmaog(nxp, my_max), gamma(nxp, my_max), &
           elvmax(nxp, my_max), ugws(nxp, my_max), vgws(nxp, my_max)
   real :: elvmax_gpu(nxp, my_max), ugws_gpu(nxp, my_max), vgws_gpu(nxp, my_max)
   real :: cdmbgwd(2)
   integer :: zmtnblck(nxp, my_max)
   integer :: zmtnblck_gpu(nxp, my_max)
   real :: garea(nxp, my_max), hpbl(nxp, my_max)

   integer :: i, j, jj, async_id, n, k
   integer, dimension(34) :: seed

   real :: kpbl_r(nxp, my_max), kdt_r, zmtnblck_r(nxp, my_max)

   real :: time1, time2, ct, gt

   real :: bnv2lm(nxp, lev, my_max), rdelks, ro(nxp, lev, my_max)

   async_id = 1
   seed = (/10004, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,&
           &25, 26, 27, 28, 29, 30, 31, 32, 33/)
   !call random_seed()
   call random_seed(put=seed)
   call random_number(dvdtc)
   call random_number(dudtc)
   call random_number(dtdtc)
   call random_number(utc)
   call random_number(vtc)
   call random_number(ttc)
   call random_number(qtc)
   call random_number(kpbl_r)!
   call random_number(prsi)
   call random_number(del)
   call random_number(prsl)
   call random_number(prslk)
   call random_number(phii)
   call random_number(phil)
   call random_number(dta)
   call random_number(kdt_r)!
   call random_number(hprime)
   call random_number(oc)
   call random_number(oa4)
   call random_number(clx)
   call random_number(theta)
   call random_number(sigmaog)
   call random_number(gamma)
   call random_number(elvmax)
   call random_number(ugws)
   call random_number(vgws)
   call random_number(zmtnblck_r)!
   call random_number(garea)
   call random_number(hpbl)
   kpbl = floor(kpbl_r*6)
   kdt = floor(kdt_r*6)
   zmtnblck = floor(zmtnblck_r*6)

   hprime = hprime*1000
   elvmax = elvmax*1000
   garea = garea*10000000

   dvdtc_gpu = dvdtc
   dudtc_gpu = dudtc
   dtdtc_gpu = dtdtc
   elvmax_gpu = elvmax
   ugws_gpu = ugws
   vgws_gpu = vgws
   zmtnblck_gpu = zmtnblck

   cdmbgwd(1) = cmbk
   cdmbgwd(2) = cgwd
   ct = 0.
   do i = 1, 16
      call cpu_time(time1)
      !if (i .ne. 1) call nvtxStartRange("CPU compute")
      do jj = 1, jlistnum
         j = jlist1(jj)
         call gwdps(nxjp(j), nxp, nxp, lev, dvdtc(1, 1, jj), dudtc(1, 1, jj), dtdtc(1, 1, jj), &
                    utc(1, 1, jj), vtc(1, 1, jj), ttc(1, 1, jj), qtc(1, 1, jj), kpbl(1, jj), prsi(1, 1, jj), &
                    del(1, 1, jj), prsl(1, 1, jj), prslk(1, 1, jj), phii(1, 1, jj), phil(1, 1, jj), dta, kdt, &
                    hprime(1, jj), oc(1, jj), oa4(1, 1, jj), clx(1, 1, jj), &
                    theta(1, jj), sigmaog(1, jj), gamma(1, jj), elvmax(1, jj), ugws(1, jj), vgws(1, jj), &
                    grav, cp, con_rd, con_rv, nx, mtnvar, cdmbgwd, me, &
                    zmtnblck(1, jj), garea(1, jj), hpbl(1, jj), .true.)
      end do
      !if (i .ne. 1) call nvtxEndRange
      call cpu_time(time2)
      if (i .ne. 1) ct = ct + time2 - time1
   end do

   if (.true.) then

      !$acc enter data copyin(nxjp, nxp, lev, utc, vtc, ttc, qtc,         &
      !$acc&            kpbl, prsi, del, prsl, prslk, phii, phil,         &
      !$acc&            dta, kdt, hprime, oc, oa4, clx, theta, sigmaog,   &
      !$acc&            gamma, grav, cp,                                  &
      !$acc&            nx, mtnvar, cdmbgwd, me,                          &
      !$acc&            garea, hpbl, jlist1)
      !$acc enter data copyin(dvdtc_gpu, dudtc_gpu, dtdtc_gpu, elvmax_gpu,      &
      !$acc&            ugws_gpu, vgws_gpu, zmtnblck_gpu)
      gt = 0.
      do i = 1, 16
         call cpu_time(time1)
         !if (i .ne. 1) call nvtxStartRange("GPU compute")
         call gwdps_gpu(nxjp, nxp, nxp, lev, dvdtc_gpu, dudtc_gpu, dtdtc_gpu, &
                        utc, vtc, ttc, qtc, kpbl, prsi, del, &
                        prsl, prslk, phii, phil, dta, kdt, hprime, oc, oa4, clx, &
                        theta, sigmaog, gamma, elvmax_gpu, ugws_gpu, vgws_gpu, &
                        grav, cp, con_rd, con_rv, nx, mtnvar, cdmbgwd, me, &
                        zmtnblck_gpu, garea, hpbl, .true.)
         !if (i .ne. 1) call nvtxEndRange
         call cpu_time(time2)
         if (i .ne. 1) gt = gt + time2 - time1
      end do
      !$acc exit data copyout(dvdtc_gpu, dudtc_gpu,dtdtc_gpu, elvmax_gpu, &
      !$acc&                  ugws_gpu, vgws_gpu, zmtnblck_gpu)
      !$acc exit data delete(nxjp, nxp, lev, utc, vtc, ttc, qtc,         &
      !$acc&            kpbl, prsi, del, prsl, prslk, phii, phil,         &
      !$acc&            dta, kdt, hprime, oc, oa4, clx, theta, sigmaog,   &
      !$acc&            gamma, grav, cp,                                  &
      !$acc&            nx, mtnvar, cdmbgwd, me,                          &
      !$acc&            garea, hpbl, jlist1)

      call assert_real(dvdtc_gpu, size(dvdtc_gpu), dvdtc, size(dvdtc), &
                       1e-12, "Array dvdtc")
      call assert_real(dudtc_gpu, size(dudtc_gpu), dudtc, size(dudtc), &
                       1e-12, "Array dudtc")
      call assert_real(dtdtc_gpu, size(dtdtc_gpu), dtdtc, size(dtdtc), &
                       1e-12, "Array dtdtc")
      call assert_real(elvmax_gpu, size(elvmax_gpu), elvmax, size(elvmax), &
                       1e-12, "Array elvmax")
      call assert_real(ugws_gpu, size(ugws_gpu), ugws, size(ugws), &
                       1e-12, "Array ugws")
      call assert_real(vgws_gpu, size(vgws_gpu), vgws, size(vgws), &
                       1e-12, "Array vgws")
      call assert_integer(zmtnblck_gpu, size(zmtnblck_gpu), zmtnblck, size(zmtnblck), &
                          "Array zmtnblck")
   end if
   write (*, *) 'Timing for CPU & GPU: ', ct, gt
   write (*, *) 'speedup ratio: ', ct/gt
end subroutine gwdps_unit
