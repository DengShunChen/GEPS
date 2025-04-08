program test_sfc_diag
   implicit none
   call mpe_init
   call cons
   call sfc_diag_unit
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

subroutine sfc_diag_unit
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
   integer, parameter :: km=4
   integer :: i, ii, j, jj, async_id, n, k
   integer, dimension(34) :: seed

   real :: time1, time2, ct, gt, diff
   integer, parameter :: nxpvs = 7501
   real :: c1xpvs,c2xpvs,tbpvs(nxpvs)
   common/fpvscom/ c1xpvs,c2xpvs,tbpvs(nxpvs)

!-------------- inputs
   integer, dimension(my_max) :: myim
   real(kind=kind_phys), dimension(nxp, my_max) :: psi, tg, qsurf, prslki, qflux, &
                                                   fm, fh, fm10, fh2, fh10
   real(kind=kind_phys), dimension(nxp, lev, my_max) :: ut, vt, tt
   real(kind=kind_phys), dimension(nxp, lev*ncld, my_max) :: qt
!-------------- inputs/outputs
!--------------outputs
   real(kind=kind_phys), dimension(nxp, my_max) :: u10, v10, t2, q2, rh2, rh10
   real(kind=kind_phys), dimension(nxp, my_max) :: u10_gpu, v10_gpu, t2_gpu, &
                                                   q2_gpu, rh2_gpu, rh10_gpu
   
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
   call random_number(qsurf)
   call random_number(prslki)
   call random_number(qflux)
   call random_number(fm)
   call random_number(fh)
   call random_number(fm10)
   call random_number(fh2)
   call random_number(fh10)

   psi = psi*30000. + 70000.
   ut = ut*40. - 20.
   vt = vt*40. - 20.
   tt = tt*100. + 200.
   qt = qt*0.001
   tg = tg*30. + 260.
   qsurf = qsurf*5e-2
   fm = fm*50.
   fh = fh*50.
   fm10 = fm10*30.
   fh2 = fh*25.
   fh10 = fh10*40.
   
   
   ntrac = ncld

   u10 = 0.
   v10 = 0.
   t2 = 0.
   q2 = 0.
   rh2 = 0.
   rh10 = 0.
   u10_gpu = 0.
   v10_gpu = 0.
   t2_gpu = 0.
   q2_gpu = 0.
   rh2_gpu = 0.
   rh10_gpu = 0.
   
   do jj = 1, jlistnum
      j = jlist1(jj)
      myim(jj) = nxjp(j)
   end do

   do ii = 1, 16
      call cpu_time(time1)
      !if (ii .ne. 1) call nvtxStartRange("CPU compute")
      do jj = 1, jlistnum
         !write(*,*) 'CPU',myrank, j, nxjp(j)
      call sfc_diag(myim(jj), nx, psi(1, jj), ut(1, lev, jj), vt(1, lev, jj), &
                 tt(1, lev, jj), qt(1, lev, jj), tg(1, jj), &
                 qsurf(1, jj), u10(1, jj), v10(1, jj), t2(1, jj), &
                 q2(1, jj), prslki(1, jj), qflux(1, jj), fm(1, jj), &
                 fh(1, jj), fm10(1, jj), fh2(1, jj), fh10(1, jj), &
                 rh2(1, jj), rh10(1, jj))
      end do
      call cpu_time(time2)
      !if (ii .ne. 1) call nvtxEndRange
      if (ii .ne. 1) ct = ct + time2 - time1
   end do

      !write(*,*) 12345

   if (.true.) then
      !$acc enter data copyin(tbpvs) async(async_id)
      !$acc enter data copyin(myim, psi, ut, vt, tt, qt, tg, qsurf, prslki, qflux, &
      !$acc&      fm, fh, fm10, fh2, fh10) async(async_id)
      !$acc enter data copyin(u10_gpu, v10_gpu, t2_gpu, q2_gpu, rh2_gpu, &
      !$acc&      rh10_gpu) async(async_id)
      !$acc wait(async_id)
      
      do ii = 1, 16
         call cpu_time(time1)
         !if (ii .ne. 1) call nvtxStartRange("GPU compute")
         call sfc_diag_gpu(myim, nx, lev, ncld, psi, ut, vt, &
                    tt, qt, tg, &
                    qsurf, u10_gpu, v10_gpu, t2_gpu, &
                    q2_gpu, prslki, qflux, fm, &
                    fh, fm10, fh2, fh10, &
                    rh2_gpu, rh10_gpu, async_id)
            !if (ii .ne. 1) call nvtxEndRange
         !$acc wait(async_id)
         call cpu_time(time2)
         if (ii .ne. 1) gt = gt + time2 - time1
      end do
      !$acc exit data delete(tbpvs) async(async_id)
      !$acc exit data delete(psi, ut, vt, tt, qt, tg, qsurf, prslki, qflux, &
      !$acc&      fm, fh, fm10, fh2, fh10) async(async_id)
      !$acc exit data copyout(u10_gpu, v10_gpu, t2_gpu, q2_gpu, rh2_gpu, &
      !$acc&     rh10_gpu) async(async_id)
      !$acc wait(async_id)
   end if

   !write(*,*) 'NaN Check:', sum(u1), sum(v1), sum(t1), sum(q1),sum(kpbl),sum(hpbl)
   !write(*,*) 'u1 2', u1(1,1,2), u1_gpu(1,1,2)

   call assert_real(u10_gpu, size(u10_gpu), u10, size(u10), &
                    1e-12, "Array u10")
   call assert_real(v10_gpu, size(v10_gpu), v10, size(v10), &
                    1e-12, "Array v10")
   call assert_real(t2_gpu, size(t2_gpu), t2, size(t2), &
                    1e-12, "Array t2")
   call assert_real(q2_gpu, size(q2_gpu), q2, size(q2), &
                    1e-12, "Array q2")
   call assert_real(rh2_gpu, size(rh2_gpu), rh2, size(rh2), &
                    1e-12, "Array rh2")
   call assert_real(rh10_gpu, size(rh10_gpu), rh10, size(rh10), &
                    1e-12, "Array rh10")

   write (*, *) 'Timing for CPU & GPU: ', ct, gt
   write (*, *) 'speedup ratio: ', ct/gt, myrank

end subroutine sfc_diag_unit

