program test_sfc_drv
   implicit none
   call mpe_init
   call cons
   call set_soilveg(1,1)
   !call getrdy
   call sflx_unit
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
      !call exit(1)
   end if

end subroutine assert_integer

subroutine sflx_unit
   use index, only: nxjp, nxp, jlistnum, jlist1
   use param, only: lev, my_max, my, ncld
   use machine, only: kind_phys
   use rank, only: myrank
   use namelist_soilveg
   !use nvtx

   implicit none

   integer :: ntrac, check, myim(my_max)
   integer, parameter :: km=4
   integer :: i, ii, j, jj, async_id, n, k, iii, cnt, tcnt
   integer, dimension(34) :: seed
   integer :: islimsk(nxp, my_max)
   real(kind=kind_phys) :: time1, time2, ct, gt, diff, tol
   real(kind=kind_phys), save         :: zsoil_noah(4)
   data zsoil_noah / -0.1, -0.4, -1.0, -2.0 /
   real(kind=kind_phys), dimension(nxp, my_max)  :: islmsk_r, vtype_r, stype_r, slope_r, &
      flag_iter_r, ice_r, couple_r


!-------------- inputs
   integer ivegsrc, nsoil, idx2
   integer, dimension(nxp*my_max) :: jj_idx, i_idx
   integer, dimension(nxp, my_max) :: couple, ice, vtype, stype, slope
   logical, dimension(nxp, my_max) :: flag_iter, flag
   real(kind=kind_phys) :: delt
   real(kind=kind_phys), dimension(nxp, my_max) :: ffrozp, zlvl, swdn, solnet, lwdn, &
      sfcems, sfcprs, sfctmp, sfcspd, prcp, q2, q2sat, dqsdt2, th2, shdmin1d, alb, &
      snoalb1d
   real(kind=kind_phys), dimension(km) :: sldpth
!-------------- inputs/outputs
   real(kind=kind_phys), dimension(nxp, my_max) :: tbot, cmc, tsea, sneqv, chx, cmx, z0, &
      shdfac, snowh
   real(kind=kind_phys), dimension(nxp, my_max) :: tbot_gpu, cmc_gpu, tsea_gpu, sneqv_gpu, &
      chx_gpu, cmx_gpu, z0_gpu, shdfac_gpu, snowh_gpu
   real(kind=kind_phys), dimension(nxp, km, my_max) :: stsoil, smsoil, slsoil
   real(kind=kind_phys), dimension(nxp, km, my_max) :: stsoil_gpu, smsoil_gpu, slsoil_gpu
!--------------outputs
   integer, dimension(nxp, my_max) :: nroot
   integer, dimension(nxp, my_max) :: nroot_gpu
   real(kind=kind_phys), dimension(nxp, my_max) :: albedo, eta, sheat, ec, edir, ett, &
      esnow, drip, dew, beta, etp, ssoil, flx1, flx2, flx3, runoff1, runoff2, runoff3, &
      snomlt, sncovr, rc, pc, rsmin, xlai, rcs, rct, rcq, rcsoil, soilw, soilm, smcwlt, &
      smcdry, smcref, smcmax
   real(kind=kind_phys), dimension(nxp, my_max) :: albedo_gpu, eta_gpu, sheat_gpu, ec_gpu, &
      edir_gpu, ett_gpu, esnow_gpu, drip_gpu, dew_gpu, beta_gpu, etp_gpu, ssoil_gpu, &
      flx1_gpu, flx2_gpu, flx3_gpu, runoff1_gpu, runoff2_gpu, runoff3_gpu, &
      snomlt_gpu, sncovr_gpu, rc_gpu, pc_gpu, rsmin_gpu, xlai_gpu, rcs_gpu, rct_gpu, &
      rcq_gpu, rcsoil_gpu, soilw_gpu, soilm_gpu, smcwlt_gpu, &
      smcdry_gpu, smcref_gpu, smcmax_gpu
   real(kind=kind_phys), dimension(nxp, km, my_max) :: et
   real(kind=kind_phys), dimension(nxp, km, my_max) :: et_gpu
   
   real(kind=kind_phys), dimension(km) :: sldpth2, stsoil2, smsoil2, slsoil2, et2
   
   
   async_id = 1
   seed = (/10006, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,&
           &25, 26, 27, 28, 29, 30, 31, 32, 33/)
   !call random_seed()
   call random_seed(put=seed)
   ct = 0.
   gt = 0.
   cnt = 0
   tcnt = 0
   do ii = 1, 16
      call random_number(islmsk_r)
      call random_number(flag_iter_r)
      call random_number(ffrozp)
      call random_number(zlvl)
      call random_number(swdn)
      call random_number(solnet)
      call random_number(lwdn)
      call random_number(sfcems)
      call random_number(sfcprs)
      call random_number(sfctmp)
      call random_number(sfcspd)
      call random_number(prcp)
      call random_number(q2)
      call random_number(q2sat)
      call random_number(dqsdt2)
      call random_number(th2)
      call random_number(vtype_r)
      call random_number(stype_r)
      call random_number(slope_r)
      call random_number(shdmin1d)
      call random_number(alb)
      call random_number(snoalb1d)
      call random_number(tbot)
      call random_number(cmc)
      call random_number(tsea)
      call random_number(stsoil)
      call random_number(smsoil)
      call random_number(slsoil)
      call random_number(sneqv)
      call random_number(chx)
      call random_number(cmx)
      call random_number(z0)
      call random_number(shdfac)
      call random_number(snowh)
      call random_number(ice_r)
      call random_number(couple_r)
      

      stype = floor(stype_r*15)+1
      vtype = floor(vtype_r*16)+1
      !islimsk = floor(islmsk_r*3)
      slope = floor(slope_r*8)+1
      
      do jj = 1, jlistnum
         do i = 1, nxp
            if (islmsk_r(i, jj) >= 0. .and. islmsk_r(i, jj) < 0.680) then
               islimsk(i, jj) = 0
            elseif (islmsk_r(i, jj) >= 0.680 .and. islmsk_r(i, jj) < 0.963) then
               islimsk(i, jj) = 1
            else
               islimsk(i, jj) = 2
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
            if (ffrozp(i, jj) .gt. 0.5) then
               ffrozp(i, jj) = 1.
            else
               ffrozp(i, jj) = 0.
            end if
            if (ice_r(i, jj) .gt. 0.5) then
               ice(i, jj) = 1
            else
               ice(i, jj) = 0
            end if
            if (couple_r(i, jj) .gt. 0.5) then
               couple(i, jj) = 1
            else
               couple(i, jj) = 0
            end if
            
            sldpth(1) = - zsoil_noah(1)
            do k = 2, km
            sldpth(k) = zsoil_noah(k-1) - zsoil_noah(k)
            enddo
            
         end do
      end do


      nsoil = km
      couple = 1
      ice = 0
      zlvl = zlvl*3. + 20.
      swdn = swdn*1.
      solnet = solnet*1.
      lwdn = lwdn*1.
      sfcems = sfcems*1.
      sfcprs = sfcprs*30000. + 70000.
      sfctmp = sfctmp*100. + 200.
      sfcspd = sfcspd*26. + 1.
      prcp = prcp*1.4
      q2 = q2*1e-3
      q2sat = q2sat*2.7e-2
      dqsdt2 = dqsdt2*1.6e-3
      th2 = th2*100. + 200.
      shdmin1d = shdmin1d*1.
      alb = alb*1.
      snoalb1d = snoalb1d*1.
      tbot = tbot*30. + 241.
      cmc = cmc*1e-3
      tsea = tsea*30. + 260.
      smsoil = smsoil*1.
      stsoil = stsoil*10. + 260.
      slsoil = slsoil*1.
      sneqv = sneqv*0.15
      chx = chx*1.4
      cmx = cmx*1.4
      z0 = z0*1.
      shdfac = shdfac*1.
      snowh = snowh*1.
      ivegsrc = 1
      delt = 720.

      tbot_gpu = tbot
      cmc_gpu = cmc
      tsea_gpu = tsea
      sneqv_gpu = sneqv
      chx_gpu = chx
      cmx_gpu = cmx
      z0_gpu = z0
      shdfac_gpu = shdfac
      snowh_gpu = snowh
      stsoil_gpu = stsoil
      smsoil_gpu = smsoil
      slsoil_gpu = slsoil

      nroot = 0
      albedo = 0.
      eta = 0.
      sheat = 0.
      ec = 0.
      edir = 0.
      ett = 0.
      esnow = 0.
      drip = 0.
      dew = 0.
      beta = 0.
      etp = 0.
      ssoil = 0.
      flx1 = 0.
      flx2 = 0.
      flx3 = 0.
      runoff1 = 0.
      runoff2 = 0.
      runoff3 = 0.
      snomlt = 0.
      sncovr = 0.
      rc = 0.
      pc = 0.
      rsmin = 0.
      xlai = 0.
      rcs = 0.
      rct = 0.
      rcq = 0.
      rcsoil = 0.
      soilw = 0.
      soilm = 0.
      smcwlt = 0.
      smcdry = 0.
      smcref = 0.
      smcmax = 0.
      et = 0.

      nroot_gpu = 0
      albedo_gpu = 0.
      eta_gpu = 0.
      sheat_gpu = 0.
      ec_gpu = 0.
      edir_gpu = 0.
      ett_gpu = 0.
      esnow_gpu = 0.
      drip_gpu = 0.
      dew_gpu = 0.
      beta_gpu = 0.
      etp_gpu = 0.
      ssoil_gpu = 0.
      flx1_gpu = 0.
      flx2_gpu = 0.
      flx3_gpu = 0.
      runoff1_gpu = 0.
      runoff2_gpu = 0.
      runoff3_gpu = 0.
      snomlt_gpu = 0.
      sncovr_gpu = 0.
      rc_gpu = 0.
      pc_gpu = 0.
      rsmin_gpu = 0.
      xlai_gpu = 0.
      rcs_gpu = 0.
      rct_gpu = 0.
      rcq_gpu = 0.
      rcsoil_gpu = 0.
      soilw_gpu = 0.
      soilm_gpu = 0.
      smcwlt_gpu = 0.
      smcdry_gpu = 0.
      smcref_gpu = 0.
      smcmax_gpu = 0.
      et_gpu = 0.


      do jj = 1, jlistnum
         j = jlist1(jj)
         myim(jj) = nxjp(j)
      end do

      do jj = 1, jlistnum
         do i = 1, myim(jj)
            flag(i, jj) = (islimsk(i, jj) == 1)
         end do
      end do
      !if (ii .ne. 1) call nvtxStartRange("CPU compute")
      do jj = 1, jlistnum
         do i = 1, myim(jj)
            if (flag(i, jj) .and. flag_iter(i, jj)) then
               do k = 1, km
                  sldpth2(k) = sldpth(k)
                  stsoil2(k) = stsoil(i, k, jj)
                  smsoil2(k) = smsoil(i, k, jj)
                  slsoil2(k) = slsoil(i, k, jj)
                  et2(k) = 0
               end do
               call cpu_time(time1)
               call sflx                                                     &
                 ( nsoil, couple(i, jj), ice(i, jj), ffrozp(i, jj), delt, zlvl(i, jj), &
                   sldpth2, swdn(i, jj), solnet(i, jj), lwdn(i, jj), &
                   sfcems(i, jj), sfcprs(i, jj), sfctmp(i, jj), sfcspd(i, jj), &
                   prcp(i, jj), q2(i, jj), q2sat(i, jj), dqsdt2(i, jj), th2(i, jj), &
                   ivegsrc, vtype(i, jj), stype(i, jj), slope(i, jj), shdmin1d(i, jj), &
                   alb(i, jj), snoalb1d(i, jj), &
      !  ---  input/outputs: &
                   tbot(i, jj), cmc(i, jj), tsea(i, jj), stsoil2, &
                   smsoil2, slsoil2, sneqv(i, jj), chx(i, jj), &
                   cmx(i, jj), z0(i, jj), &
      !  ---  outputs: &
                   nroot(i, jj), shdfac(i, jj), snowh(i, jj), albedo(i, jj), eta(i, jj), &
                   sheat(i, jj), ec(i, jj), edir(i, jj), et2, ett(i, jj), &
                   esnow(i, jj), drip(i, jj), dew(i, jj), beta(i, jj), etp(i, jj), &
                   ssoil(i, jj), flx1(i, jj), flx2(i, jj), flx3(i, jj), runoff1(i, jj), &
                   runoff2(i, jj), runoff3(i, jj), snomlt(i, jj), sncovr(i, jj), &
                   rc(i, jj), pc(i, jj), rsmin(i, jj), xlai(i, jj), rcs(i, jj), &
                   rct(i, jj), rcq(i, jj), rcsoil(i, jj), soilw(i, jj), soilm(i, jj), &
                   smcwlt(i, jj), smcdry(i, jj), smcref(i, jj), smcmax(i, jj) &
                   )
               call cpu_time(time2)
               if (ii .ne. 1) ct = ct + time2 - time1
               do k = 1, km
                  stsoil(i, k, jj) = stsoil2(k)
                  smsoil(i, k, jj) = smsoil2(k)
                  slsoil(i, k, jj) = slsoil2(k)
                  et(i, k, jj) = et2(k)
               end do
            end if
         end do
      end do
      !if (ii .ne. 1) call nvtxEndRange


      if (.true.) then
         !$acc enter data copyin(slope_data, bb, drysmc, f11, maxsmc, refsmc, &
         !$acc&      satpsi, satdk, satdw, wltsmc, qtz, rsmtbl, rgltbl, hstbl, &
         !$acc&      snupx, lai_data, nroot_data) async(async_id)
         !$acc enter data copyin(myim, flag, flag_iter, couple, ice, ffrozp, &
         !$acc&      zlvl, sldpth, swdn, solnet, lwdn, sfcems, sfcprs, sfctmp, &
         !$acc&      sfcspd, prcp, q2, q2sat, dqsdt2, th2, vtype, stype, &
         !$acc&      slope, shdmin1d, alb, snoalb1d) async(async_id)
         !$acc enter data copyin(tbot_gpu, cmc_gpu, tsea_gpu, stsoil_gpu, &
         !$acc&      smsoil_gpu, slsoil_gpu, sneqv_gpu, chx_gpu, cmx_gpu, z0_gpu, &
         !$acc&      shdfac_gpu, snowh_gpu) async(async_id)
         !$acc enter data copyin(nroot_gpu, albedo_gpu, eta_gpu, sheat_gpu, &
         !$acc&      ec_gpu, edir_gpu, et_gpu, ett_gpu, esnow_gpu, drip_gpu, &
         !$acc&      dew_gpu, beta_gpu, etp_gpu, ssoil_gpu, flx1_gpu, flx2_gpu, &
         !$acc&      flx3_gpu, runoff1_gpu, runoff2_gpu, runoff3_gpu, snomlt_gpu, &
         !$acc&      sncovr_gpu, rc_gpu, pc_gpu, rsmin_gpu, xlai_gpu, rcs_gpu, &
         !$acc&      rct_gpu, rcq_gpu, rcsoil_gpu, soilw_gpu, soilm_gpu, &
         !$acc&      smcwlt_gpu, smcdry_gpu, smcref_gpu, smcmax_gpu) async(async_id)
         !$acc wait(async_id)
         call cpu_time(time1)
         !if (ii .ne. 1) call nvtxStartRange("GPU compute")
         tcnt = tcnt + 1
         idx2 = 0
         !$acc enter data copyin(idx2) async(async_id)
         !$acc parallel loop gang vector collapse(2) private(iii, jj, i) async(async_id)
         do jj = 1, jlistnum
            do i = 1, nxp
               if (i .le. myim(jj)) then
                  if (flag_iter(i, jj) .and. flag(i, jj)) then
                     !$acc atomic capture
                     idx2 = idx2 + 1
                     iii = idx2
                     !$acc end atomic
                     jj_idx(iii) = jj
                     i_idx(iii) = i
                 end if
              end if   ! end if_flag_iter_and_flag_block
            end do   ! end do_i_loop
         end do

         !$acc update self(idx2) async(async_id)
         !$acc wait(async_id)

         call sflx_gpu &
       !  ---  inputs: &
            (nsoil, nxp, myim, flag, flag_iter, &
            couple, ice, ffrozp, delt, zlvl, &
            sldpth, swdn, solnet, lwdn, &
            sfcems, sfcprs, sfctmp, sfcspd, &
            prcp, q2, q2sat, dqsdt2, th2, &
            ivegsrc, vtype, stype, slope, shdmin1d, &
            alb, snoalb1d, &
   !  ---  input/outputs: &
             tbot_gpu, cmc_gpu, tsea_gpu, &
             stsoil_gpu, smsoil_gpu, slsoil_gpu, &
             sneqv_gpu, chx_gpu, cmx_gpu, z0_gpu, &
   !  ---  outputs: &
             nroot_gpu, shdfac_gpu, snowh_gpu, &
             albedo_gpu, eta_gpu, sheat_gpu, &
             ec_gpu, edir_gpu, et_gpu, ett_gpu, &
             esnow_gpu, drip_gpu, dew_gpu, beta_gpu, &
             etp_gpu, ssoil_gpu, flx1_gpu, flx2_gpu, &
             flx3_gpu, runoff1_gpu, runoff2_gpu, &
             runoff3_gpu, snomlt_gpu, sncovr_gpu, &
             rc_gpu, pc_gpu, rsmin_gpu, xlai_gpu, &
             rcs_gpu, rct_gpu, rcq_gpu, rcsoil_gpu, &
             soilw_gpu, soilm_gpu, smcwlt_gpu, &
             smcdry_gpu, smcref_gpu, smcmax_gpu, &
             async_id)
         !if (ii .ne. 1) call nvtxEndRange
         !$acc exit data delete(idx2) async(async_id)
         !$acc wait(async_id)
         call cpu_time(time2)
         !$acc exit data delete(slope_data, bb, drysmc, f11, maxsmc, refsmc, &
         !$acc&     satpsi, satdk, satdw, wltsmc, qtz, rsmtbl, rgltbl, hstbl, &
         !$acc&     snupx, lai_data, nroot_data) async(async_id)
         !$acc exit data delete(myim, flag, flag_iter, couple, ice, ffrozp, &
         !$acc&     zlvl, sldpth, swdn, solnet, lwdn, sfcems, sfcprs, sfctmp, &
         !$acc&     sfcspd, prcp, q2, q2sat, dqsdt2, th2, vtype, stype, &
         !$acc&     slope, shdmin1d, alb, snoalb1d) async(async_id)
         !$acc exit data copyout(tbot_gpu, cmc_gpu, tsea_gpu, stsoil_gpu, &
         !$acc&     smsoil_gpu, slsoil_gpu, sneqv_gpu, chx_gpu, cmx_gpu, z0_gpu, &
         !$acc&     shdfac_gpu, snowh_gpu) async(async_id)
         !$acc exit data copyout(nroot_gpu, albedo_gpu, eta_gpu, sheat_gpu, &
         !$acc&     ec_gpu, edir_gpu, et_gpu, ett_gpu, esnow_gpu, drip_gpu, &
         !$acc&     dew_gpu, beta_gpu, etp_gpu, ssoil_gpu, flx1_gpu, flx2_gpu, &
         !$acc&     flx3_gpu, runoff1_gpu, runoff2_gpu, runoff3_gpu, snomlt_gpu, &
         !$acc&     sncovr_gpu, rc_gpu, pc_gpu, rsmin_gpu, xlai_gpu, rcs_gpu, &
         !$acc&     rct_gpu, rcq_gpu, rcsoil_gpu, soilw_gpu, soilm_gpu, &
         !$acc&     smcwlt_gpu, smcdry_gpu, smcref_gpu, smcmax_gpu) async(async_id)
         !$acc wait(async_id)
         if (ii .ne. 1) gt = gt + time2 - time1
      end if
      
      !write(*,*) check
      if (.true.) then
         do jj = 1, jlistnum
            do i = 1, myim(jj)
               if (flag(i, jj) .and. flag_iter(i, jj)) then
                  do k = 1, 4
                     if ((stsoil(i, k, jj) - stsoil_gpu(i, k, jj))/stsoil(i, k, jj) .gt. 1e-10) then
                        cnt = cnt + 1
                        if (myrank .eq. 0) then
                           write(*,*) k, jj, i, stsoil(i, k, jj), stsoil_gpu(i, k, jj)
                           !write(*,*) 'start', k, jj, i, &
                           !      nsoil, couple(i, jj), ice(i, jj), ffrozp(i, jj), delt, zlvl(i, jj), &
                           !      sldpth(1, i, jj), sldpth(2, i, jj), sldpth(3, i, jj), sldpth(4, i, jj), &
                           !      swdn(i, jj), solnet(i, jj), lwdn(i, jj), &
                           !      sfcems(i, jj), sfcprs(i, jj), sfctmp(i, jj), sfcspd(i, jj), &
                           !      prcp(i, jj), q2(i, jj), q2sat(i, jj), dqsdt2(i, jj), th2(i, jj), &
                           !      ivegsrc, vtype(i, jj), stype(i, jj), slope(i, jj), shdmin1d(i, jj), &
                           !      alb(i, jj), snoalb1d(i, jj), &
                     !  ---  input/outputs: &
                           !       tbot_gpu(i, jj), cmc_gpu(i, jj), tsea_gpu(i, jj), &
                           !       stsoil(1, i, jj), stsoil_gpu(1, i, jj), stsoil(2, i, jj), stsoil_gpu(2, i, jj), &
                           !       stsoil(3, i, jj), stsoil_gpu(3, i, jj), stsoil(4, i, jj), stsoil_gpu(4, i, jj), &
                           !       smsoil_gpu(1, i, jj), smsoil_gpu(2, i, jj), smsoil_gpu(3, i, jj), smsoil_gpu(4, i, jj), &
                           !       slsoil(1, i, jj), slsoil_gpu(1, i, jj), slsoil(2, i, jj), slsoil_gpu(2, i, jj), &
                           !       slsoil(3, i, jj), slsoil_gpu(3, i, jj), slsoil(4, i, jj), slsoil_gpu(4, i, jj), &
                           !       sneqv_gpu(i, jj), chx_gpu(i, jj), cmx_gpu(i, jj), z0_gpu(i, jj), &
                     !  ---  outputs: &
                           !       nroot_gpu(i, jj), shdfac_gpu(i, jj), snowh_gpu(i, jj), 'end'
                           !write(*,*) ' '
                        end if
                     end if
                  end do
               end if
            end do
         end do
      end if
      !if (myrank .eq. 0) write(*,*) 'tcnt, cnt', tcnt, cnt

      !write(*,*) 'NaN Check:', sum(u1), sum(v1), sum(t1), sum(q1),sum(kpbl),sum(hpbl)
      !write(*,*) 'u1 2', u1(1,1,2), u1_gpu(1,1,2)
      call assert_real(tbot_gpu, size(tbot_gpu), tbot, size(tbot), &
                       1e-10, "Array tbot")
      call assert_real(cmc_gpu, size(cmc_gpu), cmc, size(cmc), &
                       1e-10, "Array cmc")
      call assert_real(tsea_gpu, size(tsea_gpu), tsea, size(tsea), &
                       1e-10, "Array tsea")
      call assert_real(stsoil_gpu, size(stsoil_gpu), stsoil, size(stsoil), &
                       1e-10, "Array stsoil")
      call assert_real(smsoil_gpu, size(smsoil_gpu), smsoil, size(smsoil), &
                       1e-10, "Array smsoil")
      call assert_real(slsoil_gpu, size(slsoil_gpu), slsoil, size(slsoil), &
                       1e-10, "Array slsoil")
      call assert_real(sneqv_gpu, size(sneqv_gpu), sneqv, size(sneqv), &
                       1e-10, "Array sneqv")
      call assert_real(chx_gpu, size(chx_gpu), chx, size(chx), &
                       1e-10, "Array chx")
      call assert_real(cmx_gpu, size(cmx_gpu), cmx, size(cmx), &
                       1e-10, "Array cmx")
      call assert_real(z0_gpu, size(z0_gpu), z0, size(z0), &
                       1e-10, "Array z0")
      call assert_integer(nroot_gpu, size(nroot_gpu), nroot, size(nroot), &
                       "Array nroot")
      call assert_real(shdfac_gpu, size(shdfac_gpu), shdfac, size(shdfac), &
                       1e-10, "Array shdfac")
      call assert_real(snowh_gpu, size(snowh_gpu), snowh, size(snowh), &
                       1e-10, "Array snowh")
      call assert_real(albedo_gpu, size(albedo_gpu), albedo, size(albedo), &
                       1e-10, "Array albedo")
      call assert_real(eta_gpu, size(eta_gpu), eta, size(eta), &
                       1e-8, "Array eta")
      call assert_real(sheat_gpu, size(sheat_gpu), sheat, size(sheat), &
                       1e-8, "Array sheat")
      call assert_real(ec_gpu, size(ec_gpu), ec, size(ec), &
                       1e-10, "Array ec")
      call assert_real(edir_gpu, size(edir_gpu), edir, size(edir), &
                       1e-10, "Array edir")
      call assert_real(et_gpu, size(et_gpu), et, size(et), &
                       1e-10, "Array et")
      call assert_real(ett_gpu, size(ett_gpu), ett, size(ett), &
                       1e-10, "Array ett")
      call assert_real(esnow_gpu, size(esnow_gpu), esnow, size(esnow), &
                       1e-8, "Array esnow")
      call assert_real(drip_gpu, size(drip_gpu), drip, size(drip), &
                       1e-10, "Array drip")
      call assert_real(dew_gpu, size(dew_gpu), dew, size(dew), &
                       1e-10, "Array dew")
      call assert_real(beta_gpu, size(beta_gpu), beta, size(beta), &
                       1e-10, "Array beta")
      call assert_real(etp_gpu, size(etp_gpu), etp, size(etp), &
                       1e-8, "Array etp")
      call assert_real(ssoil_gpu, size(ssoil_gpu), ssoil, size(ssoil), &
                       1e-10, "Array ssoil")
      call assert_real(flx1_gpu, size(flx1_gpu), flx1, size(flx1), &
                       1e-10, "Array flx1")
      call assert_real(flx2_gpu, size(flx2_gpu), flx2, size(flx2), &
                       1e-10, "Array flx2")
      call assert_real(flx3_gpu, size(flx3_gpu), flx3, size(flx3), &
                       1e-10, "Array flx3")
      call assert_real(runoff1_gpu, size(runoff1_gpu), runoff1, size(runoff1), &
                       1e-10, "Array runoff1")
      call assert_real(runoff2_gpu, size(runoff2_gpu), runoff2, size(runoff2), &
                       1e-10, "Array runoff2")
      call assert_real(runoff3_gpu, size(runoff3_gpu), runoff3, size(runoff3), &
                       1e-10, "Array runoff3")
      call assert_real(snomlt_gpu, size(snomlt_gpu), snomlt, size(snomlt), &
                       1e-10, "Array snomlt")
      call assert_real(sncovr_gpu, size(sncovr_gpu), sncovr, size(sncovr), &
                       1e-10, "Array sncovr")
      call assert_real(rc_gpu, size(rc_gpu), rc, size(rc), &
                       1e-10, "Array rc")
      call assert_real(pc_gpu, size(pc_gpu), pc, size(pc), &
                       1e-10, "Array pc")
      call assert_real(rsmin_gpu, size(rsmin_gpu), rsmin, size(rsmin), &
                       1e-10, "Array rsmin")
      call assert_real(xlai_gpu, size(xlai_gpu), xlai, size(xlai), &
                       1e-10, "Array xlai")
      call assert_real(rcs_gpu, size(rcs_gpu), rcs, size(rcs), &
                       1e-10, "Array rcs")
      call assert_real(rct_gpu, size(rct_gpu), rct, size(rct), &
                       1e-10, "Array rct")
      call assert_real(rcq_gpu, size(rcq_gpu), rcq, size(rcq), &
                       1e-10, "Array rcq")
      call assert_real(rcsoil_gpu, size(rcsoil_gpu), rcsoil, size(rcsoil), &
                       1e-10, "Array rcsoil")
      call assert_real(soilw_gpu, size(soilw_gpu), soilw, size(soilw), &
                       1e-10, "Array soilw")
      call assert_real(soilm_gpu, size(soilm_gpu), soilm, size(soilm), &
                       1e-10, "Array soilm")
      call assert_real(smcwlt_gpu, size(smcwlt_gpu), smcwlt, size(smcwlt), &
                       1e-10, "Array smcwlt")
      call assert_real(smcdry_gpu, size(smcdry_gpu), smcdry, size(smcdry), &
                       1e-10, "Array smcdry")
      call assert_real(smcref_gpu, size(smcref_gpu), smcref, size(smcref), &
                       1e-10, "Array smcref")
      call assert_real(smcmax_gpu, size(smcmax_gpu), smcmax, size(smcmax), &
                       1e-10, "Array smcmax")
   end do      

   write (*, *) 'Timing for CPU & GPU: ', ct, gt
   write (*, *) 'speedup ratio: ', ct/gt, myrank

end subroutine sflx_unit

