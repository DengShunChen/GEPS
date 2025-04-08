program test_moninedmf
   implicit none
   call mpe_init
   call cons
   call mfpbl_unit
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

subroutine mfpbl_unit
   use radn, only: ntcw
   use index, only: nxjp, nxp, jlistnum, jlist1
   use param, only: lev, my_max, my, ncld
   use const, only: grav, cp, nx, mtnvar, tofd, cmbk, cgwd
   use physcons, only: con_rd, con_rv
   use machine, only: kind_phys
   use const, only: RTYPE
   use rank, only: myrank
   !USE, INTRINSIC :: IEEE_ARITHMETIC
   implicit none

   real(kind=kind_phys) xmf(nxp, lev, my_max), xmf_gpu(nxp, lev, my_max)
   real(kind=kind_phys) tcko(nxp, lev, my_max), tcko_gpu(nxp, lev, my_max)
   real(kind=kind_phys) ucko(nxp, lev, my_max), ucko_gpu(nxp, lev, my_max)
   real(kind=kind_phys) vcko(nxp, lev, my_max), vcko_gpu(nxp, lev, my_max)
   real(kind=kind_phys) qcko(nxp, lev, ncld, my_max), qcko_gpu(nxp, lev, ncld, my_max)
   real(kind=kind_phys) zl(nxp, lev, my_max), zm(nxp, lev + 1, my_max), &
      thvx(nxp, lev, my_max), sflx(nxp, my_max), &
      u1(nxp, lev, my_max), v1(nxp, lev, my_max), &
      t1(nxp, lev, my_max), q1(nxp, lev, ncld, my_max), &
      ustar(nxp, my_max), wstar(nxp, my_max)
   logical cnvflg(nxp, my_max)
   real(kind=kind_phys) hpbl(nxp, my_max), hpbl_gpu(nxp, my_max)
   real(kind=kind_phys) kpbl_r(nxp, my_max)
   integer kpbl(nxp, my_max), kpbl_gpu(nxp, my_max)
   integer :: ntrac, myim(my_max)
   integer :: i, ii, j, jj, async_id, n, k, kk, im, j1, cnt, nancnt
   integer, dimension(34) :: seed
   real(kind=kind_phys) dt2, num
   real(kind=kind_phys), allocatable, dimension(:, :) :: zl0, zm0, xmf0, tcko0, ucko0, vcko0, thvx0
   real(kind=kind_phys), allocatable, dimension(:, :, :) :: qcko0

   real :: time1, time2, ct, gt, diff

   async_id = 1
   seed = (/10006, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,&
           &25, 26, 27, 28, 29, 30, 31, 32, 33/)
   !call random_seed()
   call random_seed(put=seed)
   do ii = 1, 16
      call random_number(zl)
      call random_number(zm)
      call random_number(thvx)
      call random_number(q1)
      call random_number(t1)
      call random_number(u1)
      call random_number(v1)
      call random_number(sflx)
      call random_number(ustar)
      call random_number(wstar)
      call random_number(hpbl)
      call random_number(kpbl_r)
      call random_number(tcko)
      call random_number(qcko)
      call random_number(ucko)
      call random_number(vcko)

      if (.true.) then
         zl = zl*40000.
         zm = zm*40000.
         thvx = thvx*100.+200.
         q1 = q1*5.
         t1 = t1*10.
         u1 = u1*60.-30.
         v1 = v1*60.-30.
         hpbl = hpbl*160.+20.
         kpbl_r = kpbl_r*20.+1.
         sflx = -sflx*1e-2 - 4e-2
         ustar = ustar*0.6 + 0.2
         wstar = wstar*1e-5 + 5.5e-5
         tcko = tcko*-2e-5
         qcko = qcko*-3e-5
         ucko = ucko*-2e-5
         vcko = vcko*-2e-5
      end if

      kpbl = floor(kpbl_r)
      kpbl_gpu = kpbl
      cnvflg = .true.
      tcko_gpu = tcko
      qcko_gpu = qcko
      ucko_gpu = ucko
      vcko_gpu = vcko
      hpbl_gpu = hpbl
      xmf = 0.
      xmf_gpu = 0.
      dt2 = 10.

      ntrac = ncld

      ct = 0.
      do jj = 1, jlistnum
         j = jlist1(jj)
         im = nxjp(j)
         allocate (zl0(im, lev))
         allocate (zm0(im, lev + 1))
         allocate (thvx0(im, lev))
         allocate (xmf0(im, lev))
         allocate (tcko0(im, lev))
         allocate (ucko0(im, lev))
         allocate (vcko0(im, lev))
         allocate (qcko0(im, lev, ntrac))

         do k = 1, lev
            do i = 1, im
               zl0(i, k) = zl(i, k, jj)
               zm0(i, k) = zm(i, k, jj)
               thvx0(i, k) = thvx(i, k, jj)
               xmf0(i, k) = xmf(i, k, jj)
               tcko0(i, k) = tcko(i, k, jj)
               ucko0(i, k) = ucko(i, k, jj)
               vcko0(i, k) = vcko(i, k, jj)
               do kk = 1, ntrac
                  qcko0(i, k, kk) = qcko(i, k, kk, jj)
               end do
            end do
         end do
         do i = 1, im
            zm0(i, lev + 1) = zm(i, lev + 1, jj)
         end do

         call cpu_time(time1)
         call mfpbl(im, nxp, lev, ntrac, dt2, cnvflg(1, jj), zl0, zm0, &
                    thvx0, q1(1, 1, 1, jj), t1(1, 1, jj), u1(1, 1, jj), v1(1, 1, jj), &
                    hpbl(1, jj), kpbl(1, jj), sflx(1, jj), ustar(1, jj), wstar(1, jj), &
                    xmf0, tcko0, qcko0, ucko0, vcko0)
         call cpu_time(time2)

         do k = 1, lev
            do i = 1, im
               xmf(i, k, jj) = xmf0(i, k)
               tcko(i, k, jj) = tcko0(i, k)
               ucko(i, k, jj) = ucko0(i, k)
               vcko(i, k, jj) = vcko0(i, k)
               do kk = 1, ntrac
                  qcko(i, k, kk, jj) = qcko0(i, k, kk)
               end do
            end do
         end do

         deallocate (zl0)
         deallocate (zm0)
         deallocate (thvx0)
         deallocate (xmf0)
         deallocate (tcko0)
         deallocate (ucko0)
         deallocate (vcko0)
         deallocate (qcko0)

         if (ii .ne. 1) ct = ct + time2 - time1
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do jj = 1, jlistnum
         j1 = jlist1(jj)
         myim(jj) = nxjp(j1)
      end do

      !$acc enter data copyin(myim,nxp,lev,dt2,cnvflg,zl,zm,thvx, &
      !$acc&                  q1,t1,u1,v1,sflx,ustar,wstar)
      !$acc enter data copyin(hpbl_gpu,kpbl_gpu,xmf_gpu, &
      !$acc&                  tcko_gpu,qcko_gpu,ucko_gpu,vcko_gpu)
      gt = 0.
      call cpu_time(time1)

      call mfpbl_gpu(myim, nxp, lev, ntrac, dt2, cnvflg, zl, zm, thvx, q1, t1, u1, v1, &
                     hpbl_gpu, kpbl_gpu, sflx, ustar, wstar, xmf_gpu, &
                     tcko_gpu, qcko_gpu, ucko_gpu, vcko_gpu, -1)
      call cpu_time(time2)
      if (ii .ne. 1) gt = gt + time2 - time1
      !$acc exit data copyout(hpbl_gpu,kpbl_gpu,xmf_gpu, &
      !$acc&                  tcko_gpu,qcko_gpu,ucko_gpu,vcko_gpu)
      !$acc exit data delete(myim,nxp,lev,dt2,cnvflg,zl,zm,thvx, &
      !$acc&                 q1,t1,u1,v1,sflx,ustar,wstar)


      call assert_integer(kpbl_gpu, size(kpbl_gpu), kpbl, size(kpbl), &
                          "Array kpbl")
      call assert_real(xmf_gpu, size(xmf_gpu), xmf, size(xmf), &
                       1e-8, "Array xmf")
      call assert_real(tcko_gpu, size(tcko_gpu), tcko, size(tcko), &
                       1e-8, "Array tcko")
      call assert_real(qcko_gpu, size(qcko_gpu), qcko, size(qcko), &
                       1e-6, "Array qcko")
      call assert_real(ucko_gpu, size(ucko_gpu), ucko, size(ucko), &
                       1e-8, "Array ucko")
      call assert_real(vcko_gpu, size(vcko_gpu), vcko, size(vcko), &
                       1e-8, "Array vcko")
      call assert_real(hpbl_gpu, size(hpbl_gpu), hpbl, size(hpbl), &
                       1e-8, "Array hpbl")
   end do
   write (*, *) 'Timing for CPU & GPU: ', ct, gt
   write (*, *) 'speedup ratio: ', ct/gt, myrank

end subroutine mfpbl_unit

