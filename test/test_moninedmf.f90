program test_moninedmf
   implicit none
   call mpe_init
   call cons
   call moninedmf_unit
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

subroutine moninedmf_unit
   use radn, only: ntcw
   use index, only: nxjp, nxp, jlistnum, jlist1, nxjp_acc
   use param, only: lev, my_max, my, ncld
   use const, only: grav, cp, nx, mtnvar, tofd, cmbk, cgwd
   use physcons, only: con_rd, con_rv
   use machine, only: kind_phys
   use const, only: RTYPE
   use rank, only: myrank
   !use cusparse
   !use nvtx

   implicit none

   integer :: ntrac, myim(my_max)
   real :: dta
   real(kind=kind_phys) u1_gpu(nxp, lev, my_max), v1_gpu(nxp, lev, my_max)
   real(kind=kind_phys) u1(nxp, lev, my_max), v1(nxp, lev, my_max), &
      swh(nxp, lev, my_max), hlw(nxp, lev, my_max)
   real(kind=kind_phys) t1(nxp, lev, my_max), q1(nxp, lev, ncld, my_max)
   real(kind=kind_phys) t1_gpu(nxp, lev, my_max), q1_gpu(nxp, lev, ncld, my_max)
   real(kind=RTYPE) pk2(nxp, lev, my_max)
   real(kind=kind_phys) xmu(nxp, my_max), u10(nxp, my_max), &
      v10(nxp, my_max), fm(nxp, my_max), &
      fh(nxp, my_max), tg(nxp, my_max)
   real(kind=kind_phys) rb(nxp, my_max), z0rl(nxp, my_max), &
      heat(nxp, my_max), evap(nxp, my_max), &
      stress(nxp, my_max), sfcw(nxp, my_max)
   real(kind=kind_phys) prsi(nxp, lev + 1, my_max), del(nxp, lev, my_max), &
      prsl(nxp, lev, my_max), prslk(nxp, lev, my_max), &
      phil(nxp, lev, my_max), phii(nxp, lev + 1, my_max)
   real(kind=kind_phys) hpbl(nxp, my_max), hpbl_gpu(nxp, my_max)
   integer ::           kpbl(nxp, my_max), kpbl_gpu(nxp, my_max)
   !type(cusparseHandle) :: sparsehandle
   integer :: sparsehandle
   integer :: i, ii, j, jj, async_id, n, k, istat
   integer, dimension(34) :: seed

   real :: time1, time2, ct, gt, diff

   async_id = 1
   seed = (/10006, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,&
           &25, 26, 27, 28, 29, 30, 31, 32, 33/)
   !call random_seed()
   call random_seed(put=seed)
   ct = 0.
   gt = 0.
   !istat = cusparseCreate(sparsehandle)
   do ii = 1, 16
      call random_number(u1)
      call random_number(v1)
      call random_number(swh)
      call random_number(hlw)
      call random_number(t1)
      call random_number(q1)
      call random_number(pk2)
      call random_number(xmu)!
      call random_number(u10)
      call random_number(v10)
      call random_number(fm)
      call random_number(fh)
      call random_number(tg)
      call random_number(rb)
      call random_number(z0rl)
      call random_number(heat)!
      call random_number(evap)
      call random_number(stress)
      call random_number(sfcw)
      call random_number(prsi)
      call random_number(del)
      call random_number(prsl)
      call random_number(prslk)
      call random_number(phii)
      call random_number(phil)
      call random_number(dta)

      u1 = u1*40.-20.
      v1 = v1*40.-20.
      t1 = t1*300.
      q1 = q1*5.
      swh = swh*5e-4 - 1e-4
      hlw = hlw*200.
      xmu = xmu*2.
      pk2 = pk2
      rb = rb*20.-10.
      z0rl = z0rl*200.
      u10 = u10*20.-10.
      v10 = v10*20.-10.
      fm = fm*30.
      fh = fh*30.
      tg = tg*50.+280.
      heat = heat*0.6 - 0.1
      evap = evap*2e-4
      stress = stress*2.
      sfcw = sfcw*20.
      prsi = prsi*100000.
      del = del*3000.
      prsl = prsl*100000.
      prslk = prslk*2.
      phii = phii*600000.
      phil = phil*540000.
      dta = 10.

      u1_gpu = u1
      v1_gpu = v1
      t1_gpu = t1
      q1_gpu = q1
      hpbl = 0.
      hpbl_gpu = 0.
      kpbl = 0
      kpbl_gpu = 0

      ntrac = ncld
      do jj = 1, jlistnum
         j = jlist1(jj)
         myim(jj) = nxjp(j)
      end do
      call cpu_time(time1)
      !if (ii .ne. 1) call nvtxStartRange("CPU compute")
      do jj = 1, jlistnum
         j = jlist1(jj)
         !write(*,*) 'CPU',myrank, j, nxjp(j)
         call moninedmf(nxp, nxjp(j), lev, ntrac, ntcw, u1(1, 1, jj), v1(1, 1, jj), &
                        t1(1, 1, jj), q1(1, 1, 1, jj), swh(1, 1, jj), hlw(1, 1, jj), &
                        xmu(1, jj), pk2(1, lev, jj), rb(1, jj), z0rl(1, jj), &
                        u10(1, jj), v10(1, jj), fm(1, jj), fh(1, jj), tg(1, jj), &
                        heat(1, jj), evap(1, jj), stress(1, jj), sfcw(1, jj), &
                        kpbl(1, jj), prsi(1, 1, jj), del(1, 1, jj), prsl(1, 1, jj), &
                        prslk(1, 1, jj), phii(1, 1, jj), phil(1, 1, jj), dta, hpbl(1, jj))
      end do
      call cpu_time(time2)
      !if (ii .ne. 1) call nvtxEndRange
      if (ii .ne. 1) ct = ct + time2 - time1

      !write(*,*) 12345

      if (.true.) then
         !$acc enter data copyin(nxp, myim, lev, ntrac, ntcw,                &
         !$acc&            swh, hlw, xmu, pk2, rb, z0rl, u10, v10,           &
         !$acc&            fm, fh, tg, heat, evap, stress, sfcw, prsi, del,  &
         !$acc&            prsl, prslk, phii, phil, dta, jlist1, nxjp_acc)   &
         !$acc&            async(async_id)
         !$acc enter data copyin(u1_gpu, v1_gpu, t1_gpu, q1_gpu, kpbl_gpu,   &
         !$acc&            hpbl_gpu) async(async_id)
         call cpu_time(time1)
         !if (ii .ne. 1) call nvtxStartRange("GPU compute")
         call moninedmf_gpu(nxp, myim, lev, ntrac, ntcw, u1_gpu, v1_gpu, &
                            t1_gpu, q1_gpu, swh, hlw, &
                            xmu, pk2, rb, z0rl, &
                            u10, v10, fm, fh, tg, &
                            heat, evap, stress, sfcw, &
                            kpbl_gpu, prsi, del, prsl, &
                            prslk, phii, phil, dta, hpbl_gpu, sparsehandle, async_id)
         !if (ii .ne. 1) call nvtxEndRange
         !$acc wait(async_id)
         call cpu_time(time2)
         if (ii .ne. 1) gt = gt + time2 - time1
         !$acc exit data copyout(u1_gpu, v1_gpu, t1_gpu, q1_gpu, kpbl_gpu,   &
         !$acc&           hpbl_gpu) async(async_id)
         !$acc exit data delete(nxp, myim, lev, ntrac, ntcw,                 &
         !$acc&            swh, hlw, xmu, pk2, rb, z0rl, u10, v10,           &
         !$acc&            fm, fh, tg, heat, evap, stress, sfcw, prsi, del,  &
         !$acc&            prsl, prslk, phii, phil, dta, jlist1, nxjp_acc) &
         !$acc&           async(async_id)
         !$acc wait(async_id)
         do jj = 1, jlistnum
            j = jlist1(jj)
            do i = 1, nxjp(j)
               if (kpbl(i, jj) .ne. kpbl_gpu(i, jj)) then
                  if (myrank .eq. 0) write (*, *) jj, i, kpbl(i, jj), kpbl_gpu(i, jj)
                  !write(*,*) sum(t1(i,:,jj)), sum(t1_gpu(i,:,jj))
               end if
            end do
         end do

         call assert_real(u1_gpu, size(u1_gpu), u1, size(u1), &
                          1e-10, "Array u1")
         call assert_real(v1_gpu, size(v1_gpu), v1, size(v1), &
                          1e-10, "Array v1")
         call assert_real(t1_gpu, size(t1_gpu), t1, size(t1), &
                          1e-10, "Array t1")
         call assert_real(q1_gpu, size(q1_gpu), q1, size(q1), &
                          1e-10, "Array q1")
         call assert_real(hpbl_gpu, size(hpbl_gpu), hpbl, size(hpbl), &
                          1e-10, "Array hpbl")
         call assert_integer(kpbl_gpu, size(kpbl_gpu), kpbl, size(kpbl), &
                             "Array kpbl")
      end if
   end do
   write (*, *) 'Timing for CPU & GPU: ', ct, gt
   write (*, *) 'speedup ratio: ', ct/gt, myrank
   !istat = cusparseDestroy(sparsehandle)

end subroutine moninedmf_unit

