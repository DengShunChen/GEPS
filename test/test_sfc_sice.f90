program test_sfc_sice
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
   call sfc_sice_unit
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

subroutine sfc_sice_unit
   use index, only: nxjp, nxp, jlistnum, jlist1, nxdef_2d
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

   implicit none

   integer :: ntrac, itstp, nx
   integer, parameter :: km=4, km_soil = 4
   integer :: i, ii, j, jj, async_id, n, k, nxj
   integer, dimension(34) :: seed

   real :: time1, time2, ct, gt, diff
   integer, parameter :: nxpvs = 7501
   real :: c1xpvs,c2xpvs,tbpvs(nxpvs)
   common/fpvscom/ c1xpvs,c2xpvs,tbpvs(nxpvs)

   real(kind=kind_phys), dimension(nxp, my_max) :: ss_adj, rld_adj, sld_adj
   real(kind=kind_phys) :: radus, tpi, d2r, cosl(my), xx, cp, rgas, grav, ptop, &
      pkxr, r, g, t850, ppd, ppu, ddd, p850, cc
   real(kind=kind_phys), dimension(nxp, lev+1, my_max) :: phii
   real(kind=kind_phys), dimension(nxp, my_max) :: topo, pss, ps, hgt, z0rl, tsurf
   real(kind=kind_phys) :: pkx(nxp, lev, my_max)
   real(kind=kind_phys), dimension(nxp, my_max) :: rb, stress, fm, fh, sfcw, fm10, fh10, fh2


!-------------- inputs
   integer :: mom4ice, lsm
   integer, dimension(my_max) :: myim
   integer, dimension(nxp, my_max) :: islmsk
   logical, dimension(nxp, my_max) :: flag_iter
   real(kind=kind_phys) :: dth
   real(kind=kind_phys), dimension(nxp, my_max) :: psi, rld, ss, &
      sld, srflag, cd, cdq, prsl1, prslki, ddvel
!-------------- inputs/outputs
   real(kind=kind_phys), dimension(nxp, my_max) :: sheleg, tprcp, ep1d
   real(kind=kind_phys), dimension(nxp, my_max) :: zice_gpu, cice_gpu, & 
      xtice_gpu, sheleg_gpu, tg_gpu, tprcp_gpu, ep1d_gpu
   real(kind=kind_phys), dimension(nxp, km, my_max) :: stc_gpu
!--------------outputs
   real(kind=kind_phys), dimension(nxp, my_max) :: snwdph, qsurf, snomt, &
                                                   gfx, qflux, hflux
   real(kind=kind_phys), dimension(nxp, my_max) :: snwdph_gpu, qsurf_gpu, &
      snomt_gpu, gfx_gpu, qflux_gpu, hflux_gpu
   
   async_id = 1
   nx = nxp
   seed = (/10006, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,&
           &25, 26, 27, 28, 29, 30, 31, 32, 33/)
   call random_seed(put=seed)
   ct = 0.
   gt = 0.

   call random_number(ss_adj)
   call random_number(rld_adj)
   call random_number(sld_adj)
   
   do jj = 1, jlistnum
      j = jlist1(jj)
      myim(jj) = nxjp(j)
   end do

   grav = 9.806649999999999
   rgas = 287.0285714285714
   cp = 1004.600000000000
   ptop = 0.1000000000000000
   r = rgas
   g = grav
   itstp = 1
   
   ss_adj = ss_adj*900.
   rld_adj = rld_adj*200* + 200.
   sld_adj = sld_adj*1000.
   
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
            tt(i,k,jj)  = tt(i,k,jj)*pk(i,k,jj) / (1.0+0.608*qt(i,k,jj))
         end do
      end do
   end do
      
   topo = sgeo
   pss = pt
   sheleg = snr
   ss = ss_adj
   rld = rld_adj
   snwdph = sndepth
   sld = sld_adj
      
   dth = 720.
   mom4ice = .true.
   lsm = 1
   ntrac = ncld

   do jj = 1, jlistnum
      do k = 1, lev
         do i = 1, nxp
            if (i .le. myim(jj)) then
               pkxr = pk(i, k, jj)
               pkxr=log(pk(i, k, jj))
               pkxr = pkxr*(cp/r)
               pkxr=exp(pkxr)*1000.
               pkx(i, k, jj) = pkxr
            end if
         end do
      end do
   end do
   p850 = 850.
   do jj = 1, jlistnum
      do i = 1, nxp
         if (i .le. myim(jj)) then
            ps(i, jj) = pss(i, jj) + ptop
            hgt(i, jj) = (phi(i, lev, jj) - topo(i, jj))/g
      
      
            if (ocean(i, jj)) islmsk(i, jj) = 0
            if (land(i, jj)) islmsk(i, jj) = 1
            if (ice(i, jj)) islmsk(i, jj) = 2
            if (itstp .eq. 1 .and. ocean(i, jj)) &
               z0(i, jj) = ustar(i, jj)*ustar(i, jj)*0.014/g ! get z0 from ustar
            prslki(i, jj) = pk2(i, lev, jj)/pk(i, lev, jj)   !(ps/p1)**r/cp
            tprcp(i, jj) = totalp(i, jj)/1000.  ! dth precip (m)
            ddvel(i, jj) = 0.
      
            z0rl(i, jj) = z0(i, jj)*100.         ! transport from m to cm
            psi(i, jj) = ps(i, jj)*100.  ! surface pressure (mb to pa)
            prsl1(i, jj) = pkx(i, lev, jj)*100.     !pa
            
            tsurf(i, jj) = tg(i, jj)
            flag_iter(i, jj) = .true.
            ep1d(i, jj) = 0.
            hflux(i, jj) = 0.
            qflux(i, jj) = 0.
            
            t850 = tt(i, lev, jj)
         do k = lev, 2, -1
            if (i .le. myim(jj)) then
               ppd = pkx(i, k, jj)
               ppu = pkx(i, k - 1, jj)
               ddd = (p850 - ppd)*(p850 - ppu)
               if (ddd .lt. 0.) then
                  cc = (log(p850/ppd))/(log(ppu/ppd))
                  t850 = (tt(i, k - 1, jj) - tt(i, k, jj))*cc &
                     + tt(i, k, jj)
               end if
            end if
         end do
            srflag(i, jj) = 0.
            if (t850 .le. 273.16) then
               srflag(i, jj) = 1.
            end if
         end if
      end do
   end do


   do jj = 1, jlistnum
      call sfc_diff(myim(jj), nx, psi(1, jj), ut(1, lev, jj), &
                    vt(1, lev, jj), tt(1, lev, jj), &
                    qt(1, lev, jj), &
                    hgt(1, jj), snwdph(1, jj), tg(1, jj), z0rl(1, jj), cd(1, jj), &
                    cdq(1, jj), rb(1, jj), prsl1(1, jj), prslki(1, jj), islmsk(1, jj), &
                    stress(1, jj), fm(1, jj), fh(1, jj), &
                    ustar(1, jj), sfcw(1, jj), ddvel(1, jj), fm10(1, jj), &
                    fh2(1, jj), fh10(1, jj), sigmaf(1, jj), &
                    ivegtyp(1, jj), shdmax(1, jj), ivegsrc, &
                    tsurf(1, jj), flag_iter(1, jj), .false.)
   end do


   zice_gpu = zice
   cice_gpu = cice
   xtice_gpu = xtice
   sheleg_gpu = sheleg
   tg_gpu = tg
   tprcp_gpu = tprcp
   stc_gpu = stc
   ep1d_gpu = ep1d

   snwdph = 0.
   qsurf = 0.
   snomt = 0.
   gfx = 0.
   qflux = 0.
   hflux = 0.
   snwdph_gpu = 0.
   qsurf_gpu = 0.
   snomt_gpu = 0.
   gfx_gpu = 0.
   qflux_gpu = 0.
   hflux_gpu = 0.
      
   do jj = 1, jlistnum
      j = jlist1(jj)
      myim(jj) = nxjp(j)
   end do

   do ii = 1, 16
      call cpu_time(time1)
      !if (ii .ne. 1) call nvtxStartRange("CPU compute")
      do jj = 1, jlistnum
         !write(*,*) 'CPU',myrank, j, nxjp(j)
         call sfc_sice(myim(jj), nx, km, psi(1, jj), ut(1, lev, jj), &
                       vt(1, lev, jj), tt(1, lev, jj), &
                       qt(1, lev, jj), dth, sfemis(1, jj), rld(1, jj), &
                       ss(1, jj), sld(1, jj), srflag(1, jj), &
                       cd(1, jj), cdq(1, jj), prsl1(1, jj), prslki(1, jj), &
                       islmsk(1, jj), ddvel(1, jj), &
                       flag_iter(1, jj), mom4ice, lsm, zice(1, jj), &
                       cice(1, jj), xtice(1, jj), sheleg(1, jj), &
                       tg(1, jj), tprcp(1, jj), stc(1, 1, jj), ep1d(1, jj), &
                       snwdph(1, jj), qsurf(1, jj), snomt(1, jj), gfx(1, jj), &
                       qflux(1, jj), hflux(1, jj))
      end do
      call cpu_time(time2)
      !if (ii .ne. 1) call nvtxEndRange
      if (ii .ne. 1) ct = ct + time2 - time1
   end do

      !write(*,*) 12345

   if (.true.) then
      !$acc enter data copyin(tbpvs) async(async_id)
      !$acc enter data copyin(myim, psi, ut, vt, tt, qt, sfemis, rld, ss, &
      !$acc&      sld, srflag, cd, cdq, prsl1, prslki, islmsk, ddvel, &
      !$acc&      flag_iter, nx, dth) async(async_id)
      !$acc enter data copyin(zice_gpu, cice_gpu, xtice_gpu, sheleg_gpu, &
      !$acc&      tg_gpu, tprcp_gpu, stc_gpu, ep1d_gpu, snwdph_gpu, qsurf_gpu, &
      !$acc&      snomt_gpu, gfx_gpu, qflux_gpu, hflux_gpu) async(async_id)
      do ii = 1, 16
         !$acc wait(async_id)
         call cpu_time(time1)
         !if (ii .ne. 1) call nvtxStartRange("GPU compute")
         !write(*,*) 'CPU',myrank, j, nxjp(j)
         call sfc_sice_gpu(myim, nx, km, lev, ncld, psi, ut, &
                       vt, tt, &
                       qt, dth, sfemis, rld, &
                       ss, sld, srflag, &
                       cd, cdq, prsl1, prslki, &
                       islmsk, ddvel, &
                       flag_iter, mom4ice, lsm, zice_gpu, &
                       cice_gpu, xtice_gpu, sheleg_gpu, &
                       tg_gpu, tprcp_gpu, stc_gpu, ep1d_gpu, &
                       snwdph_gpu, qsurf_gpu, snomt_gpu, gfx_gpu, &
                       qflux_gpu, hflux_gpu, async_id) 
         !$acc wait(async_id)
         !if (ii .ne. 1) call nvtxEndRange
         call cpu_time(time2)
         if (ii .ne. 1) gt = gt + time2 - time1
      end do
      !$acc exit data delete(tbpvs) async(async_id)
      !$acc exit data delete(myim, psi, ut, vt, tt, qt, sfemis, rld, ss, &
      !$acc&     sld, srflag, cd, cdq, prsl1, prslki, islmsk, ddvel, &
      !$acc&     flag_iter, nx, dth) async(async_id)
      !$acc exit data copyout(zice_gpu, cice_gpu, xtice_gpu, sheleg_gpu, &
      !$acc&     tg_gpu, tprcp_gpu, stc_gpu, ep1d_gpu, snwdph_gpu, qsurf_gpu, &
      !$acc&     snomt_gpu, gfx_gpu, qflux_gpu, hflux_gpu) async(async_id)
      !$acc wait(async_id)
   end if




   call assert_real(zice_gpu, size(zice_gpu), zice, size(zice), &
                    1e-12, "Array zice")
   call assert_real(cice_gpu, size(cice_gpu), cice, size(cice), &
                    1e-12, "Array cice")
   call assert_real(xtice_gpu, size(xtice_gpu), xtice, size(xtice), &
                    1e-12, "Array xtice")
   call assert_real(sheleg_gpu, size(sheleg_gpu), sheleg, size(sheleg), &
                    1e-12, "Array sheleg")
   call assert_real(tg_gpu, size(tg_gpu), tg, size(tg), &
                    1e-12, "Array tg")
   call assert_real(tprcp_gpu, size(tprcp_gpu), tprcp, size(tprcp), &
                    1e-12, "Array tprcp")
   call assert_real(stc_gpu, size(stc_gpu), stc, size(stc), &
                    1e-12, "Array stc")
   call assert_real(ep1d_gpu, size(ep1d_gpu), ep1d, size(ep1d), &
                    1e-12, "Array ep1d")
   call assert_real(snwdph_gpu, size(snwdph_gpu), snwdph, size(snwdph), &
                    1e-12, "Array snwdph")
   call assert_real(qsurf_gpu, size(qsurf_gpu), qsurf, size(qsurf), &
                    1e-12, "Array qsurf")
   call assert_real(snomt_gpu, size(snomt_gpu), snomt, size(snomt), &
                    1e-12, "Array snomt")
   call assert_real(gfx_gpu, size(gfx_gpu), gfx, size(gfx), &
                    1e-10, "Array gfx")
   call assert_real(qflux_gpu, size(qflux_gpu), qflux, size(qflux), &
                    1e-10, "Array qflux")
   call assert_real(hflux_gpu, size(hflux_gpu), hflux, size(hflux), &
                    1e-12, "Array hflux")

   write (*, *) 'Timing for CPU & GPU: ', ct, gt
   write (*, *) 'speedup ratio: ', ct/gt, myrank

end subroutine sfc_sice_unit

