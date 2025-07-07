program test_sfc_drv
   implicit none
   call mpe_init
   call cons
   call set_soilveg(1,1)
   !call getrdy
   call sfc_drv_unit
   call mpe_finalize

end program

subroutine assert_real(actual, n_actual, desired, n_desired, rtol, err_msg)
   use const, only: RTYPE
   use rank

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
      if (present(err_msg)) print *, err_msg, myrank
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

subroutine sfc_drv_unit
   use radn, only: ntcw
   use index, only: nxjp, nxp, jlistnum, jlist1, nxjp_acc
   use param, only: lev, my_max, my, ncld
   use const, only: grav, cp, nx, mtnvar, tofd, cmbk, cgwd
   use physcons, only: con_rd, con_rv
   use machine, only: kind_phys
   use const, only: RTYPE
   use rank, only: myrank
   use namelist_soilveg
   !use nvtx

   implicit none

   integer :: ntrac, total_count, check
   integer, parameter :: km=4
   integer :: i, ii, j, jj, async_id, n, k, iii, error_count, error, idx(1), idx_check
   integer, dimension(34) :: seed

   real :: time1, time2, ct, gt, diff, tol
   real, dimension(nxp, my_max)  :: islmsk_r, istyp_r, ivegtyp_r, islopetyp_r, &
      flag_iter_r, flag_guess_r
   integer, parameter :: nxpvs = 7501
   real :: c1xpvs,c2xpvs,tbpvs(nxpvs)
   common/fpvscom/ c1xpvs,c2xpvs,tbpvs(nxpvs)
   integer, dimension(nxp*my_max) :: i_idx, jj_idx

!-------------- inputs
   integer isot, ivegsrc, j, io, jo
   integer, dimension(my_max) :: myim
   integer, dimension(nxp, my_max) :: islmsk, istyp, ivegtyp, islopetyp
   logical, dimension(nxp, my_max) :: flag_iter, flag_guess, flag
   real(kind=kind_phys) :: dth
   real(kind=kind_phys), dimension(nxp, my_max) :: psi, sigmaf, sfemis, rld, sld, &
      ss, tgclim, cd, cdq, prsl1, prslki, hgt, ddvel, shdmin, shdmax, snoalb, alb
   real(kind=kind_phys), dimension(nxp, lev, my_max) :: ut, vt, tt
   real(kind=kind_phys), dimension(nxp, lev*ncld, my_max) :: qt
!-------------- inputs/outputs
   real(kind=kind_phys), dimension(nxp, my_max) :: sheleg, snwdph, tg, tprcp, &
      srflag, canopy, tsurf, z0rl
   real(kind=kind_phys), dimension(nxp, my_max) :: sheleg_gpu, snwdph_gpu, &
      tg_gpu, tprcp_gpu, srflag_gpu, canopy_gpu, tsurf_gpu, z0rl_gpu
   real(kind=kind_phys), dimension(nxp, km, my_max) :: smc, stc, slc
   real(kind=kind_phys), dimension(nxp, km, my_max) :: smc_gpu, stc_gpu, slc_gpu
!--------------outputs
   real(kind=kind_phys), dimension(nxp, my_max) :: sncover, qsurf, gfx, drain, &
      qflux, hflux, ep1d, runof, albedo2
   real(kind=kind_phys), dimension(nxp, my_max) :: sncover_gpu, qsurf_gpu, gfx_gpu, &
   drain_gpu, qflux_gpu, hflux_gpu, ep1d_gpu, runof_gpu, albedo2_gpu
   
   async_id = 1
   seed = (/10006, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,&
           &25, 26, 27, 28, 29, 30, 31, 32, 33/)
   !call random_seed()
   call random_seed(put=seed)
   ct = 0.
   gt = 0.
   error_count = 0
   do ii = 1, 16
   call random_number(psi)
   call random_number(ut)
   call random_number(vt)
   call random_number(tt)
   call random_number(qt)
   call random_number(istyp_r)
   call random_number(ivegtyp_r)
   call random_number(sigmaf)
   call random_number(sfemis)
   call random_number(rld)
   call random_number(sld)
   call random_number(ss)
   call random_number(tgclim)
   call random_number(cd)
   call random_number(cdq)
   call random_number(prsl1)
   call random_number(prslki)
   call random_number(hgt)
   call random_number(islmsk_r)
   call random_number(ddvel)
   call random_number(islopetyp_r)
   call random_number(shdmin)
   call random_number(shdmax)
   call random_number(snoalb)
   call random_number(alb)
   call random_number(flag_iter_r)
   call random_number(flag_guess_r)
   call random_number(sheleg)
   call random_number(snwdph)
   call random_number(tg)
   call random_number(tprcp)
   call random_number(srflag)
   call random_number(smc)
   call random_number(stc)
   call random_number(slc)
   call random_number(canopy)
   call random_number(tsurf)
   call random_number(z0rl)

   istyp = floor(istyp_r*15)+1
   ivegtyp = floor(ivegtyp_r*16)+1
   !islmsk = floor(islmsk_r*3)
   islopetyp = floor(islopetyp_r*8)+1
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
         if (flag_guess_r(i, jj) .gt. 0.762) then
            flag_guess(i, jj) = .true.
         else
            flag_guess(i, jj) = .false.
         end if
         if (srflag(i, jj) .gt. 0.5) then
            srflag(i, jj) = 1.
         else
            srflag(i, jj) = 0.
         end if
      end do
   end do
   psi = psi*30000. + 70000.
   ut = ut*40. - 20.
   vt = vt*40. - 20.
   tt = tt*100. + 200.
   qt = qt*0.001
   sigmaf = sigmaf*1.
   sfemis = sfemis*1.
   rld = rld*1.
   sld = sld*1.
   ss = ss*1.
   tgclim = tgclim*30. + 241.
   cd = cd*0.05
   cdq = cdq*0.05
   prsl1 = prsl1*30000. + 70000.
   prslki = prslki*0.001 + 1.
   hgt = hgt*3. + 20.
   ddvel = ddvel*1.
   shdmin = shdmin*1.
   shdmax = shdmax*1.
   snoalb = snoalb*1.
   alb = alb*1.
   sheleg = sheleg*150.
   snwdph = snwdph*1000.
   tg = tg*30. + 260.
   tprcp = tprcp*1.
   !srflag = 1.
   smc = smc*1.
   stc = stc*10. + 260.
   slc = slc*1.
   canopy = canopy*1.
   tsurf = tsurf*30. + 260.
   z0rl= z0rl*100.
   
   isot = 1
   ivegsrc = 1
   dth = 720.

   sheleg_gpu = sheleg
   snwdph_gpu = snwdph
   tg_gpu = tg
   tprcp_gpu = tprcp
   srflag_gpu = srflag
   canopy_gpu = canopy
   tsurf_gpu = tsurf
   z0rl_gpu = z0rl
   smc_gpu = smc
   stc_gpu = stc
   slc_gpu = slc
   
   sncover = 0.
   qsurf = 0.
   gfx = 0.
   drain = 0.
   qflux = 0.
   hflux = 0.
   ep1d = 0.
   runof = 0.
   albedo2 = 0.
   sncover_gpu = 0.
   qsurf_gpu = 0.
   gfx_gpu = 0.
   drain_gpu = 0.
   qflux_gpu = 0.
   hflux_gpu = 0.
   ep1d_gpu = 0.
   runof_gpu = 0.
   albedo2_gpu = 0.
   

   
   do jj = 1, jlistnum
      j = jlist1(jj)
      myim(jj) = nxjp(j)
   end do
   do jj = 1, jlistnum
      do i = 1, myim(jj)
        flag(i, jj) = (islmsk(i, jj) == 1)
      end do
   end do
      call cpu_time(time1)
      !if (ii .ne. 1) call nvtxStartRange("CPU compute")
      do jj = 1, jlistnum
         !write(*,*) 'CPU',myrank, j, nxjp(j)
         call sfc_drv(myim(jj), nx, km, psi(1, jj), ut(1, lev, jj), &
                      vt(1, lev, jj), tt(1, lev, jj), &
                      qt(1, lev, jj), istyp(1, jj), &
                      ivegtyp(1, jj), sigmaf(1, jj), sfemis(1, jj), &
                      rld(1, jj), sld(1, jj), ss(1, jj), dth, &
                      tgclim(1, jj), cd(1, jj), cdq(1, jj), prsl1(1, jj), prslki(1, jj), &
                      hgt(1, jj), islmsk(1, jj), ddvel(1, jj), islopetyp(1, jj), &
                      shdmin(1, jj), shdmax(1, jj), &
                      snoalb(1, jj), alb(1, jj), flag_iter(1, jj), &
                      flag_guess(1, jj), isot, ivegsrc, sheleg(1, jj), &
                      snwdph(1, jj), tg(1, jj), tprcp(1, jj), &
                      srflag(1, jj), smc(1, 1, jj), &
                      stc(1, 1, jj), slc(1, 1, jj), &
                      canopy(1, jj), tsurf(1, jj), z0rl(1, jj), sncover(1, jj), &
                      qsurf(1, jj), gfx(1, jj), drain(1, jj), qflux(1, jj), &
                      hflux(1, jj), ep1d(1, jj), runof(1, jj), &
                      albedo2(1, jj), jj, io, jo)
      end do
      call cpu_time(time2)
      !if (ii .ne. 1) call nvtxEndRange
      if (ii .ne. 1) ct = ct + time2 - time1
   if (.true.) then
      !$acc enter data copyin(tbpvs) async(async_id)
      !$acc enter data copyin(slope_data, bb, drysmc, f11, maxsmc, refsmc, &
      !$acc&      satpsi, satdk, satdw, wltsmc, qtz, rsmtbl, rgltbl, hstbl, &
      !$acc&      snupx, lai_data, nroot_data) async(async_id)
      !$acc enter data copyin(myim, psi, ut, vt, tt, qt, istyp, ivegtyp, &
      !$acc&      sigmaf, sfemis, rld, sld, ss, dth, tgclim, cd, cdq, prsl1, &
      !$acc&      prslki, hgt, islmsk, ddvel, islopetyp, shdmin, shdmax, &
      !$acc&      snoalb, alb, flag_iter, flag_guess) async(async_id)
      !$acc enter data copyin(sheleg_gpu, snwdph_gpu, tg_gpu, tprcp_gpu, &
      !$acc&      srflag_gpu, smc_gpu, stc_gpu, slc_gpu, canopy_gpu, &
      !$acc&      tsurf_gpu, z0rl_gpu, sncover_gpu, qsurf_gpu, gfx_gpu, &
      !$acc&      drain_gpu, qflux_gpu, hflux_gpu, ep1d_gpu, runof_gpu, &
      !$acc&      albedo2_gpu) async(async_id)
      !!$acc wait(async_id)
      call cpu_time(time1)
      !if (ii .ne. 1) call nvtxStartRange("GPU compute")
      call sfc_drv_gpu(myim, nx, km, lev, ncld, psi, ut, &
                   vt, tt, qt, istyp, &
                   ivegtyp, sigmaf, sfemis, &
                   rld, sld, ss, dth, &
                   tgclim, cd, cdq, prsl1, prslki, &
                   hgt, islmsk, ddvel, islopetyp, &
                   shdmin, shdmax, &
                   snoalb, alb, flag_iter, &
                   flag_guess, isot, ivegsrc, sheleg_gpu, &
                   snwdph_gpu, tg_gpu, tprcp_gpu, &
                   srflag_gpu, smc_gpu, &
                   stc_gpu, slc_gpu, &
                   canopy_gpu, tsurf_gpu, z0rl_gpu, sncover_gpu, &
                   qsurf_gpu, gfx_gpu, drain_gpu, qflux_gpu, &
                   hflux_gpu, ep1d_gpu, runof_gpu, &
                   albedo2_gpu, jj, io, jo, async_id)
      !$acc wait(async_id)
      !if (ii .ne. 1) call nvtxEndRange
      call cpu_time(time2)
      if (ii .ne. 1) gt = gt + time2 - time1
      !$acc exit data delete(tbpvs) async(async_id)
      !$acc exit data delete(slope_data, bb, drysmc, f11, maxsmc, refsmc, &
      !$acc&     satpsi, satdk, satdw, wltsmc, qtz, rsmtbl, rgltbl, hstbl, &
      !$acc&     snupx, lai_data, nroot_data) async(async_id)
      !$acc exit data delete(myim, psi, ut, vt, tt, qt, istyp, ivegtyp, &
      !$acc&     sigmaf, sfemis, rld, sld, ss, dth, tgclim, cd, cdq, prsl1, &
      !$acc&     prslki, hgt, islmsk, ddvel, islopetyp, shdmin, shdmax, &
      !$acc&     snoalb, alb, flag_iter, flag_guess) async(async_id)
      !$acc exit data copyout(sheleg_gpu, snwdph_gpu, tg_gpu, tprcp_gpu, &
      !$acc&     srflag_gpu, smc_gpu, stc_gpu, slc_gpu, canopy_gpu, &
      !$acc&     tsurf_gpu, z0rl_gpu, sncover_gpu, qsurf_gpu, gfx_gpu, &
      !$acc&     drain_gpu, qflux_gpu, hflux_gpu, ep1d_gpu, runof_gpu, &
      !$acc&     albedo2_gpu) async(async_id)
      !$acc wait(async_id)
      !write(*,*) check
   end if
      if (.true.) then
      call assert_real(sheleg_gpu, size(sheleg_gpu), sheleg, size(sheleg), &
                       1e-10, "Array sheleg")
      call assert_real(snwdph_gpu, size(snwdph_gpu), snwdph, size(snwdph), &
                       1e-10, "Array snwdph")
      call assert_real(tg_gpu, size(tg_gpu), tg, size(tg), &
                       1e-10, "Array tg")
      call assert_real(tprcp_gpu, size(tprcp_gpu), tprcp, size(tprcp), &
                       1e-10, "Array tprcp")
      call assert_real(srflag_gpu, size(srflag_gpu), srflag, size(srflag), &
                       1e-10, "Array srflag")
      call assert_real(smc_gpu, size(smc_gpu), smc, size(smc), &
                       1e-8, "Array smc")
      call assert_real(stc_gpu, size(stc_gpu), stc, size(stc), &
                       1e-10, "Array stc")
      call assert_real(slc_gpu, size(slc_gpu), slc, size(slc), &
                       1e-8, "Array slc")
      call assert_real(canopy_gpu, size(canopy_gpu), canopy, size(canopy), &
                       1e-10, "Array canopy")
      call assert_real(tsurf_gpu, size(tsurf_gpu), tsurf, size(tsurf), &
                       1e-10, "Array tsurf")
      call assert_real(z0rl_gpu, size(z0rl_gpu), z0rl, size(z0rl), &
                       1e-10, "Array z0rl")
      call assert_real(sncover_gpu, size(sncover_gpu), sncover, size(sncover), &
                       1e-10, "Array sncover")
      call assert_real(qsurf_gpu, size(qsurf_gpu), qsurf, size(qsurf), &
                       1e-10, "Array qsurf")
      call assert_real(gfx_gpu, size(gfx_gpu), gfx, size(gfx), &
                       1e-10, "Array gfx")
      call assert_real(drain_gpu, size(drain_gpu), drain, size(drain), &
                       1e-10, "Array drain")
      call assert_real(qflux_gpu, size(qflux_gpu), qflux, size(qflux), &
                       1e-8, "Array qflux")
      call assert_real(hflux_gpu, size(hflux_gpu), hflux, size(hflux), &
                       1e-8, "Array hflux")
      call assert_real(ep1d_gpu, size(ep1d_gpu), ep1d, size(ep1d), &
                       1e-8, "Array ep1d")
      call assert_real(runof_gpu, size(runof_gpu), runof, size(runof), &
                       1e-8, "Array runof")
      call assert_real(albedo2_gpu, size(albedo2_gpu), albedo2, size(albedo2), &
                       1e-10, "Array albedo2")
      end if
               
      total_count = 0
      error_count = 0
      do jj = 1, jlistnum
         do i = 1, myim(jj)
            !do k = 1, 4
               if (flag_iter(i, jj) .and. flag(i, jj)) then
               total_count = total_count + 1
               if (abs((hflux(i, jj) - hflux_gpu(i, jj))/hflux(i, jj)) .gt. 1e-10) then
                  !if (myrank .eq. 0) write(*,*) jj, i, hflux(i, jj), hflux_gpu(i, jj)
                  error_count = error_count + 1
               end if
               end if
            !end do
         end do
      end do
      !if (myrank .eq. 0) write(*,*) 'error', total_count, error_count
   end do
   write (*, *) 'Timing for CPU & GPU: ', ct, gt
   write (*, *) 'speedup ratio: ', ct/gt, myrank
end subroutine sfc_drv_unit

