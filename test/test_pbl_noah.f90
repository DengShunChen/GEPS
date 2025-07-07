program test_sfc_drv
   use param
   use const

   implicit none
   integer :: no
   call mpe_init
   call cons
   call getrdy
   if (donnmi .and. taui .lt. 1.0) then
      no = 2*((jtrun + 1)/2) + (jtrun/2) + 10
      call initial(no, jtrun, jtmax, lev, nx, my, my_max, mlmax)
   end if

   call pbl_noah_unit
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

subroutine pbl_noah_unit
   use index, only: nxjp, nxp, jlistnum, jlist1, nxdef_2d, nxjp_acc
   use param, only: lev, my_max, my, ncld
   use machine, only: kind_phys
   use rank, only: myrank
   use noah, only: istyp, ivegtyp, slopetyp, isot, ivegsrc, sfalb, sigmaf, &
      shdmax, shdmin, snoalb, sfemis, canopy, zice, cice, xtice, smc, stc, slc, &
      runoff, sndepth
   use phygrid, only: land, ice, ocean, gwclim, tgclim, totalp, rs, e, eps, &
      asl, atl, ustar, tg, z0, snr, xlat
   use const, only: dsigma, RTYPE
   use grid, only: sgeo, up, vp, ttp, qm, ut, vt, tt, qt, pk, pk2, phi, pt
   use namelist_soilveg
   use nvtx
   !use nvtx

   implicit none

   integer :: ntrac, check, myim(my_max)
   integer, parameter :: km_soil=4
   integer :: i, ii, j, jj, async_id, n, k, iii, cnt, tcnt, nxj
   integer, dimension(34) :: seed
   integer :: islimsk(nxp, my_max)
   real(kind=kind_phys) :: time1, time2, ct, gt, diff, tol, tmp
   real(kind=kind_phys) :: radus, tpi, d2r, cosl(my), xx, xkapa
   integer, parameter :: nxpvs = 7501
   real :: c1xpvs,c2xpvs,tbpvs(nxpvs)
   common/fpvscom/ c1xpvs,c2xpvs,tbpvs(nxpvs)


!-------------- inputs
   integer :: itimestep, itypbl, nmmiph
   integer, dimension(my) :: ijdg
   integer, dimension(2, my) :: ipblmx
   real(kind=kind_phys) :: dta, grav, rgas, cp, xlapa, hltm, ptop, tice, hice, &
      stbo
   real(kind=kind_phys), dimension(lev) :: xkmd
   real(kind=kind_phys), dimension(2, my) :: xkmx
   real(kind=kind_phys), dimension(nxp, my_max) :: ss_adj, rld_adj, sld_adj, xmu
   real(kind=kind_phys), dimension(nxp, lev, my_max) :: upp, vpp, ttpp
   real(kind=kind_phys), dimension(nxp, lev+1, my_max) :: phii
   real(kind=RTYPE), dimension(nxp, my_max) :: pst
   real(kind=RTYPE), dimension(nxp, lev*ncld, my_max) :: qp
!-------------- inputs/outputs
   integer :: ktpbl, nmpbl
   integer :: ktpbl_gpu, nmpbl_gpu
   real(kind=kind_phys), dimension(nxp, my_max) :: tg_gpu, z0_gpu, ustar_gpu, &
      snr_gpu, canopy_gpu, zice_gpu, cice_gpu, xtice_gpu
   real(kind=kind_phys), dimension(nxp, km_soil, my_max) :: smc_gpu, stc_gpu, &
      slc_gpu
   real(kind=RTYPE), dimension(nxp, lev*ncld, my_max) :: qt_gpu
!--------------outputs
   integer, dimension(nxp, my_max) :: kpbl
   integer, dimension(nxp, my_max) :: kpbl_gpu
   real(kind=kind_phys), dimension(nxp, my_max) :: tstar, qstar, hflux, qflux, &
      t2, q2, rh2, rh10, u10, v10, fm, fh, fm10, fh2, srflag, sncover, &
      albedo2, hpbl, gfx
   real(kind=kind_phys), dimension(nxp, my_max) :: tstar_gpu, qstar_gpu, &
      hflux_gpu, qflux_gpu, t2_gpu, q2_gpu, rh2_gpu, rh10_gpu, u10_gpu, v10_gpu, &
      fm_gpu, fh_gpu, fm10_gpu, fh2_gpu, srflag_gpu, runoff_gpu, sncover_gpu, &
      sndepth_gpu, albedo2_gpu, hpbl_gpu, gfx_gpu
   real(kind=kind_phys), dimension(nxp, lev, my_max) :: dudtc, dvdtc, dtdtc, &
      dqdtc
   real(kind=kind_phys), dimension(nxp, lev, my_max) :: dudtc_gpu, dvdtc_gpu, &
      dtdtc_gpu, dqdtc_gpu
   
   
   async_id = 1
   seed = (/10006, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,&
           &25, 26, 27, 28, 29, 30, 31, 32, 33/)
   !call random_seed()
   call random_seed(put=seed)
   ct = 0.
   gt = 0.
   cnt = 0
   tcnt = 0
   do ii = 1, 1
      call random_number(ss_adj)
      call random_number(rld_adj)
      call random_number(sld_adj)
      call random_number(xmu)
      
      do jj = 1, jlistnum
         j = jlist1(jj)
         myim(jj) = nxjp(j)
      end do

      
      ss_adj = ss_adj*900.
      rld_adj = rld_adj*200* + 200.
      sld_adj = sld_adj*1000.
      xmu = xmu*2.
      
      radus = 6371000.
      tpi   = 4.0*atan(1.0)
      d2r   = tpi / 180.0
      do j = 1, my
         cosl(j) = cos(xlat(j)*d2r)
      end do
      
      do jj = 1, jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         xx = radus/cosl(j)
         call get_phi(nxjp(j),nxp,lev,ptop,cp,rgas,grav,sgeo(1,jj),      &
                     pk(1,1,jj),pk2(1,1,jj),tt(1,1,jj),qt(1,1,jj),       &
                     phii(1,1,jj),phi(1,1,jj))
         do k = 1, lev
           do i = 1, nxj
             ut(i,k,jj)  = ut(i,k,jj)*xx
             vt(i,k,jj)  = vt(i,k,jj)*xx
             upp(i,k,jj)    = up(i,k,jj)*xx
             vpp(i,k,jj)    = vp(i,k,jj)*xx
             tt(i,k,jj)  = tt(i,k,jj)*pk(i,k,jj) / (1.0+0.608*qt(i,k,jj))
             ttpp(i,k,jj) = ttp(i,k,jj) / (1.0+0.608*qp(i,k,jj))
           enddo
         enddo
      end do
      pst = pt
      qp = qm

      itypbl = 0
      itimestep = 1
      nmmiph = 12
      dta = 720.
      grav = 9.806649999999999
      rgas = 287.0285714285714
      cp = 1004.600000000000
      xkapa = 0.2857142857142857
      hltm = 2500000.000000000
      ptop = 0.1000000000000000
      tice = 271.2000000000000
      hice = 333580.0000000000
      stbo = 5.6704000000000003E-008
      ktpbl = 2
      nmpbl = 4
      
      !$acc parallel loop gang private(tmp)
      do j = 1, jlistnum
         !$acc loop vector
         do i = 1, myim(jj)
            tmp = gfx(i, jj)
         end do
      end do
      
      ktpbl_gpu = ktpbl
      nmpbl_gpu = nmpbl
      ustar_gpu = ustar
      tg_gpu = tg
      z0_gpu = z0
      snr_gpu = snr
      canopy_gpu = canopy
      zice_gpu = zice
      cice_gpu = cice
      xtice_gpu = xtice
      smc_gpu = smc
      stc_gpu = stc
      slc_gpu = slc
      dudtc_gpu = dudtc
      dvdtc_gpu = dvdtc
      dtdtc_gpu = dtdtc
      dqdtc_gpu = dqdtc
      qt_gpu = qt
      runoff_gpu = runoff
      sndepth_gpu = sndepth
     
      ntrac = ncld 
      
      kpbl = 0
      tstar = 0.
      qstar = 0.
      hflux = 0.
      qflux = 0.
      t2 = 0.
      q2 = 0.
      rh2 = 0.
      rh10 = 0.
      u10 = 0.
      v10 = 0.
      fm = 0.
      fh = 0.
      fm10 = 0.
      fh2 = 0.
      srflag = 0.
      sncover = 0.
      albedo2 = 0.
      hpbl = 0.
      gfx = 0.
      kpbl_gpu = 0
      tstar_gpu = 0.
      qstar_gpu = 0.
      hflux_gpu = 0.
      qflux_gpu = 0.
      t2_gpu = 0.
      q2_gpu = 0.
      rh2_gpu = 0.
      rh10_gpu = 0.
      u10_gpu = 0.
      v10_gpu = 0.
      fm_gpu = 0.
      fh_gpu = 0.
      fm10_gpu = 0.
      fh2_gpu = 0.
      srflag_gpu = 0.
      sncover_gpu = 0.
      albedo2_gpu = 0.
      hpbl_gpu = 0.
      gfx_gpu = 0.

      call cpu_time(time1)
      if (ii .ne. 1) call nvtxStartRange("CPU compute")
      do jj = 1, jlistnum
         j=jlist1(jj)
         call pbl_noah ( nxjp(j),nxp,lev,ktpbl,dta,grav,rgas,cp,xkapa,hltm,ptop &
                       , tice,hice,tg(1,jj),z0(1,jj),land(1,jj)                 &
                       , sgeo(1,jj),phi(1,1,jj),phii(1,1,jj),pst(1,jj),upp(1,1,jj),vpp(1,1,jj)                  &
                       , ttpp(1,1,jj),qp(1,1,jj),ut(1,1,jj),vt(1,1,jj)                  &
                       , tt(1,1,jj),qt(1,1,jj),pk(1,1,jj),pk2(1,1,jj)           &
                       , ustar(1,jj),tstar(1,jj),qstar(1,jj),e(1,1,jj)          &
                       , eps(1,1,jj),hflux(1,jj),qflux(1,jj),itimestep          &
                       , gwclim(1,jj),tgclim(1,jj),ocean(1,jj),ice(1,jj)        &
   !                       , snr(1,jj),totalp(1,jj),ss_adj,rs(1,jj),albx(1,jj)      &
                       , snr(1,jj),totalp(1,jj),ss_adj(1,jj),rs(1,jj),sfalb(1,jj)     &
                       , ipblmx(1,j),xkmx(1,j),ijdg(j),xkmd,itypbl              &
                       , t2(1,jj),q2(1,jj),rh2(1,jj),rh10(1,jj),u10(1,jj)       &
                       , v10(1,jj),fm(1,jj),fh(1,jj),fm10(1,jj),fh2(1,jj)       &
                       , srflag(1,jj),rld_adj(1,jj),stbo                              &
                       , km_soil,smc(1,1,jj),stc(1,1,jj),canopy(1,jj)           &
                       , runoff(1,jj),sigmaf(1,jj),istyp(1,jj),ivegtyp(1,jj)    &
                       , ncld,dsigma                                            &
                       , slopetyp(1,jj)                                         &
                       , slc(1,1,jj),sncover(1,jj),sndepth(1,jj)                &
                       , shdmax(1,jj),shdmin(1,jj),snoalb(1,jj),albedo2(1,jj)   &
                       , sld_adj(1,jj),zice(1,jj),cice(1,jj),xtice(1,jj)              &
                       , hpbl(1,jj),asl(1,1,jj),atl(1,1,jj),xmu(1,jj),gfx(1,jj) &
                       , kpbl(1,jj),nmpbl,nmmiph,j,isot,ivegsrc,sfemis(1,jj)    &
                       , dudtc(1,1,jj),dvdtc(1,1,jj),dtdtc(1,1,jj),dqdtc(1,1,jj))
      end do
      call cpu_time(time2)
      if (ii .ne. 1) call nvtxEndRange
      if (ii .ne. 1) ct = ct + time2 - time1


      if (.true.) then
         !$acc enter data copyin(slope_data, bb, drysmc, f11, maxsmc, refsmc, &
         !$acc&      satpsi, satdk, satdw, wltsmc, qtz, rsmtbl, rgltbl, hstbl, &
         !$acc&      snupx, lai_data, nroot_data) async(async_id)
         !$acc enter data copyin(jlist1, nxjp, nxjp_acc) async(async_id)
         !$acc enter data copyin(tbpvs) async(async_id)
         !$acc enter data copyin(land, sgeo, phi, phii, pst, upp, vpp, ttpp, &
         !$acc&      qp, ut, vt, tt, pk, pk2, e, eps, hflux, qflux, gwclim, tgclim, &
         !$acc&      ocean, ice, totalp, ss_adj, rs, sfalb, ipblmx, xkmx, ijdg, &
         !$acc&      xkmd, rld_adj, sigmaf, istyp, ivegtyp, dsigma, slopetyp, &
         !$acc&      shdmax, shdmin, snoalb, sld_adj, asl, atl, xmu, sfemis) &
         !$acc&      async(async_id)
         !$acc enter data copyin(tg_gpu, z0_gpu, qt_gpu, ustar_gpu, snr_gpu, &
         !$acc&      smc_gpu, stc_gpu, canopy_gpu, runoff_gpu, slc_gpu, &
         !$acc&      sndepth_gpu, zice_gpu, cice_gpu, xtice_gpu) async(async_id)
         !$acc enter data copyin(tstar_gpu, qstar_gpu, hflux_gpu, qflux_gpu, &
         !$acc&      t2_gpu, q2_gpu, rh2_gpu, rh10_gpu, u10_gpu, v10_gpu, &
         !$acc&      fm_gpu, fh_gpu, fm10_gpu, fh2_gpu, srflag_gpu, sncover_gpu, &
         !$acc&      albedo2_gpu, hpbl_gpu, gfx_gpu, kpbl_gpu, dudtc_gpu, &
         !$acc&      dvdtc_gpu, dtdtc_gpu, dqdtc_gpu) async(async_id)
         call cpu_time(time1)
         if (ii .ne. 1) call nvtxStartRange("GPU compute")
         call pbl_noah_gpu(nxjp, nxp, lev, ktpbl_gpu, dta, grav, rgas, cp, xkapa, hltm, ptop, &
                       tice, hice, tg_gpu, z0_gpu, land, &
                       sgeo, phi, phii, pst, upp, vpp, &
                       ttpp, qp, ut, vt, &
                       tt, qt_gpu, pk, pk2, &
                       ustar_gpu, tstar_gpu, qstar_gpu, e, &
                       eps, hflux_gpu, qflux_gpu, itimestep, &
                       gwclim, tgclim, ocean, ice, &
                       !                     , snr(1,jj),totalp(1,jj),ss_adj(1,jj),rs(1,jj),albx(1,jj)      , &
                       snr_gpu, totalp, ss_adj, rs, sfalb, &
                       ipblmx, xkmx, ijdg, xkmd, itypbl, &
                       t2_gpu, q2_gpu, rh2_gpu, rh10_gpu, u10_gpu, &
                       v10_gpu, fm_gpu, fh_gpu, fm10_gpu, fh2_gpu, &
                       srflag_gpu, rld_adj, stbo, &
                       km_soil, smc_gpu, stc_gpu, canopy_gpu, &
                       runoff_gpu, sigmaf, istyp, ivegtyp, &
                       ncld, dsigma, &
                       slopetyp, &
                       slc_gpu, sncover_gpu, sndepth_gpu, &
                       shdmax, shdmin, snoalb, albedo2_gpu, &
                       sld_adj, zice_gpu, cice_gpu, xtice_gpu, &
                       hpbl_gpu, asl, atl, xmu, gfx_gpu, &
                       kpbl_gpu, nmpbl_gpu, nmmiph, isot, ivegsrc, sfemis, &
                       dudtc_gpu, dvdtc_gpu, dtdtc_gpu, dqdtc_gpu, ntrac)
         !$acc wait(async_id)
         if (ii .ne. 1) call nvtxEndRange
         call cpu_time(time2)
         !$acc exit data delete(slope_data, bb, drysmc, f11, maxsmc, refsmc, &
         !$acc&     satpsi, satdk, satdw, wltsmc, qtz, rsmtbl, rgltbl, hstbl, &
         !$acc&     snupx, lai_data, nroot_data) async(async_id)
         !$acc exit data delete(jlist1, nxjp, nxjp_acc) async(async_id)
         !$acc exit data delete(tbpvs) async(async_id)
         !$acc exit data delete(pk, pk2, e, eps, hflux, qflux, gwclim, tgclim, &
         !$acc&     ocean, ice, totalp, ss_adj, rs, sfalb, ipblmx, xkmx, ijdg, &
         !$acc&     xkmd, rld_adj, sigmaf, istyp, ivegtyp, dsigma, slopetyp, &
         !$acc&     shdmax, shdmin, snoalb, sld_adj, asl, atl, xmu, sfemis) &
         !$acc&     async(async_id)
         !$acc exit data copyout(tg_gpu, z0_gpu, qt_gpu, ustar_gpu, snr_gpu, &
         !$acc&     smc_gpu, stc_gpu, canopy_gpu, runoff_gpu, slc_gpu, &
         !$acc&     sndepth_gpu, zice_gpu, cice_gpu, xtice_gpu) async(async_id)
         !$acc exit data copyout(tstar_gpu, qstar_gpu, hflux_gpu, qflux_gpu, &
         !$acc&     t2_gpu, q2_gpu, rh2_gpu, rh10_gpu, u10_gpu, v10_gpu, &
         !$acc&     fm_gpu, fh_gpu, fm10_gpu, fh2_gpu, srflag_gpu, sncover_gpu, &
         !$acc&     albedo2_gpu, hpbl_gpu, gfx_gpu, kpbl_gpu, dudtc_gpu, &
         !$acc&     dvdtc_gpu, dtdtc_gpu, dqdtc_gpu) async(async_id)
         !$acc wait(async_id)
         if (ii .ne. 1) gt = gt + time2 - time1
      end if

      call assert_real(ustar_gpu, size(ustar_gpu), ustar, size(ustar), &
                       1e-10, "Array ustar")
      call assert_real(tg_gpu, size(tg_gpu), tg, size(tg), &
                       1e-10, "Array tg")
      call assert_real(z0_gpu, size(z0_gpu), z0, size(z0), &
                       1e-10, "Array z0")
      call assert_real(snr_gpu, size(snr_gpu), snr, size(snr), &
                       1e-10, "Array snr")
      call assert_real(canopy_gpu, size(canopy_gpu), canopy, size(canopy), &
                       1e-1, "Array canopy")
      call assert_real(zice_gpu, size(zice_gpu), zice, size(zice), &
                       1e-10, "Array zice")
      call assert_real(cice_gpu, size(cice_gpu), cice, size(cice), &
                       1e-10, "Array cice")
      call assert_real(xtice_gpu, size(xtice_gpu), xtice, size(xtice), &
                       1e-10, "Array xtice")
      call assert_real(smc_gpu, size(smc_gpu), smc, size(smc), &
                       1e-10, "Array smc")
      call assert_real(stc_gpu, size(stc_gpu), stc, size(stc), &
                       1e-4, "Array stc")
      call assert_real(slc_gpu, size(slc_gpu), slc, size(slc), &
                       1e-10, "Array slc")
      call assert_real(dudtc_gpu, size(dudtc_gpu), dudtc, size(dudtc), &
                       1e-5, "Array dudtc")
      call assert_real(dvdtc_gpu, size(dvdtc_gpu), dvdtc, size(dvdtc), &
                       1e-5, "Array dvdtc")
      call assert_real(dtdtc_gpu, size(dtdtc_gpu), dtdtc, size(dtdtc), &
                       1e-3, "Array dtdtc")
      call assert_real(dqdtc_gpu, size(dqdtc_gpu), dqdtc, size(dqdtc), &
                       1e-5, "Array dqdtc")
      call assert_real(qt_gpu, size(qt_gpu), qt, size(qt), &
                       1e-9, "Array qt")
      call assert_real(tstar_gpu, size(tstar_gpu), tstar, size(tstar), &
                       1e-8, "Array tstar")
      call assert_real(qstar_gpu, size(qstar_gpu), qstar, size(qstar), &
                       1e-10, "Array qstar")
      call assert_real(hflux_gpu, size(hflux_gpu), hflux, size(hflux), &
                       1e-8, "Array hflux")
      call assert_real(qflux_gpu, size(qflux_gpu), qflux, size(qflux), &
                       1e-10, "Array qflux")
      call assert_real(t2_gpu, size(t2_gpu), t2, size(t2), &
                       1e-10, "Array t2")
      call assert_real(q2_gpu, size(q2_gpu), q2, size(q2), &
                       1e-10, "Array q2")
      call assert_real(rh2_gpu, size(rh2_gpu), rh2, size(rh2), &
                       1e-10, "Array rh2")
      call assert_real(rh10_gpu, size(rh10_gpu), rh10, size(rh10), &
                       1e-10, "Array rh10")
      call assert_real(u10_gpu, size(u10_gpu), u10, size(u10), &
                       1e-10, "Array u10")
      call assert_real(v10_gpu, size(v10_gpu), v10, size(v10), &
                       1e-10, "Array v10")
      call assert_real(fm_gpu, size(fm_gpu), fm, size(fm), &
                       1e-10, "Array fm")
      call assert_real(fh_gpu, size(fh_gpu), fh, size(fh), &
                       1e-10, "Array fh")
      call assert_real(fm10_gpu, size(fm10_gpu), fm10, size(fm10), &
                       1e-10, "Array fm10")
      call assert_real(fh2_gpu, size(fh2_gpu), fh2, size(fh2), &
                       1e-10, "Array fh2")
      call assert_real(srflag_gpu, size(srflag_gpu), srflag, size(srflag), &
                       1e-10, "Array srflag")
      call assert_real(runoff_gpu, size(runoff_gpu), runoff, size(runoff), &
                       1e-9, "Array runoff")
      call assert_real(sncover_gpu, size(sncover_gpu), sncover, size(sncover), &
                       1e-10, "Array sncover")
      call assert_real(sndepth_gpu, size(sndepth_gpu), sndepth, size(sndepth), &
                       1e-10, "Array sndepth")
      call assert_real(albedo2_gpu, size(albedo2_gpu), albedo2, size(albedo2), &
                       1e-10, "Array albedo2")
      call assert_real(hpbl_gpu, size(hpbl_gpu), hpbl, size(hpbl), &
                       1e-10, "Array hpbl")
      call assert_real(gfx_gpu, size(gfx_gpu), gfx, size(gfx), &
                       1e-8, "Array gfx")
      call assert_integer(kpbl_gpu, size(kpbl_gpu), kpbl, size(kpbl), &
                       "Array kpbl")
      if (ktpbl_gpu .ne.ktpbl) then
         write(*,*) ktpbl_gpu - ktpbl, 'ktpbl'
         !exit(1)
      end if
      if (nmpbl_gpu .ne.nmpbl) then
         write(*,*) nmpbl_gpu - nmpbl, 'nmpbl'
         !exit(1)
      end if
   end do      

   write (*, *) 'Timing for CPU & GPU: ', ct, gt
   write (*, *) 'speedup ratio: ', ct/gt, myrank

end subroutine pbl_noah_unit

