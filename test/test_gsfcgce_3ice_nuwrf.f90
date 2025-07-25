
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
      !call exit(1)
   end if

end subroutine assert_real

subroutine check_real(actual, n_actual, desired, n_desired, rtol, ntol, err_msg)
   use const, only: RTYPE
   use rank

   implicit none

   integer, intent(in) :: n_actual, n_desired
   real(kind=RTYPE), dimension(n_actual), intent(in) :: actual
   real(kind=RTYPE), dimension(n_desired), intent(in) :: desired
   real(kind=RTYPE), intent(in) :: rtol
   integer, intent(in) :: ntol
   character(len=*), intent(in), optional :: err_msg
   real(kind=RTYPE), parameter :: eps = 1e-15
   real(kind=RTYPE) :: rel_diff, abs_diff
   logical :: equal
   integer :: i, cnt

   if (n_actual .ne. n_desired) then
      print *, "Arrays have different lengths!"
      if (present(err_msg)) print *, err_msg
      call exit(1)
   end if

   equal = .true.
   rel_diff = 0.0
   cnt = 0
   do i = 1, n_actual
      if (abs((actual(i) - desired(i))/(desired(i))) .gt. rtol) then
         cnt = cnt + 1
      end if
   end do
   
   1000 format("number of elements > ", EN7.1, " ", A, ": ", I8)
   2000 format("Too many elements (", I4, ") > ", EN7.1, " ", A, ": ", I8)
   if (present(err_msg)) print 1000, rtol, err_msg, cnt
   if (cnt .gt. ntol) then
      if (present(err_msg)) print 2000,  ntol, rtol, err_msg, cnt
      !call exit(1)
   end if

end subroutine check_real


subroutine fall_flux_unit(SL_sedi, sat_predict, new_saturation, &
                                   use_cpm, use_declination, benchmark, ct, gt)
   use module_mp_gsfcgce_3ice_nuwrf_gpu, only: gsfcgce_3ice_nuwrf_gpu, &
       consat_s_gpu, gce_table_copyin_gpu, fall_flux_gpu
   use module_mp_gsfcgce_3ice_nuwrf, only: gsfcgce_3ice_nuwrf, fall_flux, consat_s
   use phygrid, only: land, ice, ocean, xlat, aeroclxm
   use noah, only: ivegtyp
   use grid, only: sgeo, qt, pk, tt, phi, pk2, pt
   use index, only: nxjp, nxp, jlistnum, jlist1, nxdef_2d, nxjp_acc
   use param, only: lev, my_max, my, ncld
   use machine, only: kind_phys
   use rank, only: myrank
   use const, only: dsigma, RTYPE, naero, sigma
   use radn,                only: ntcw,ntiw,ntrw,ntsw,ntgl,nthl,     &
                                  ntinc,ntrnc
   use physcons, only: con_rd, con_fvirt, con_g

   !use nvtx

   implicit none

   integer :: ntrac, check, myim(my_max)
   integer :: i, ii, j, jj, async_id, n, k, iii, cnt, tcnt, nxj, kc
   integer, dimension(34) :: seed
   real :: time1, time2
   character(len=:),allocatable :: flag_str
   !flags
   logical, intent(in) :: SL_sedi, sat_predict, new_saturation, use_cpm, use_declination
   logical, intent(in) :: benchmark
   real, intent(inout) :: ct, gt
   
   integer :: itimestep, ihail, ice2, islimsk(nxp, my_max)
   real :: dta, sdec, grav, rgas, cp, ptop
   real, dimension(my_max) :: xlat_myim
   real(kind=kind_phys), dimension(nxp, lev+1, my_max) :: phii
   real,dimension(nxp, my_max) :: rain2d, ice2d, snow2d, &
      graupel2d, sr2d, land2d, ht
   real,dimension(nxp, my_max) :: rain2d_cpu, ice2d_cpu, snow2d_cpu, &
      graupel2d_cpu, sr2d_cpu
   real,dimension(nxp, my_max) :: rain2d_gpu, ice2d_gpu, snow2d_gpu, &
      graupel2d_gpu, sr2d_gpu
   real,dimension(nxp, lev, my_max) :: th3d, qv3d, &
      qc3d, qr3d, qi3d, qs3d, qg3d, rew3d, rer3d, &
      rei3d, res3d, reg3d, rho3d, pii3d, p3d, z3d, dz3d, w3d, plt, vvel, tr
   real,dimension(nxp, lev, my_max) :: th3d_cpu, qv3d_cpu, &
      qc3d_cpu, qr3d_cpu, qi3d_cpu, qs3d_cpu, qg3d_cpu, rew3d_cpu, rer3d_cpu, &
      rei3d_cpu, res3d_cpu, reg3d_cpu
   real,dimension(nxp, lev, my_max) :: th3d_gpu, qv3d_gpu, &
      qc3d_gpu, qr3d_gpu, qi3d_gpu, qs3d_gpu, qg3d_gpu, rew3d_gpu, rer3d_gpu, &
      rei3d_gpu, res3d_gpu, reg3d_gpu
   real, dimension(nxp, lev, 1, naero) :: aero3d
   real, dimension(nxp, lev, naero, my_max) :: aero4d
   real, dimension(nxp, lev, my_max, naero) :: aero
   
   flag_str = ''
   if (SL_sedi .eq. .true.) flag_str = flag_str//', SL_sedi'
   if (sat_predict .eq. .true.) flag_str = flag_str//', sat_predict'
   if (new_saturation .eq. .true.) flag_str = flag_str//', new_saturation'
   if (use_cpm .eq. .true.) flag_str = flag_str//', use_cpm'
   if (use_declination .eq. .true.) flag_str = flag_str//', use_declination'
   if (myrank .eq. 0) write(*,*) flag_str
   
   async_id = 1
   seed = (/10006, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,&
           &25, 26, 27, 28, 29, 30, 31, 32, 33/)
   !call random_seed()
   call random_seed(put=seed)
   call random_number(vvel)
   vvel = vvel * 4e-5
   
   cnt = 0
   tcnt = 0

   grav = 9.806649999999999
   rgas = 287.0285714285714
   cp = 1004.600000000000
   ptop = 0.1000000000000000
   sdec = 0.
   ihail = 0
   ice2 = 0
   dta = 600.
   itimestep = 25
   
   do jj = 1, jlistnum
      j = jlist1(jj)
      myim(jj) = nxjp(j)
      xlat_myim(jj) = xlat(j)
   end do
   
   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef_2d(j)
      do i = 1, nxj
         if (land(i, jj)) islimsk(i, jj) = 1
         if (ocean(i, jj)) islimsk(i, jj) = 0
         if (ice(i, jj)) islimsk(i, jj) = 2
      end do
   end do
   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef_2d(j)
      !
      !   new p**capa quantities were computed in previous diabat call
      !
      call prexp_hybrid_cwb(nxjp(j), nxp, lev, ptop, sigma, pt(1, jj), &
                          pk(1, 1, jj), pk2(1, 1, jj), plt(1, 1, jj))
   end do !jj = 1,jlistnum


   do jj = 1, jlistnum
      j=jlist1(jj)
      nxj=nxdef_2d(j)
      call get_phi(nxjp(j),nxp,lev,ptop,cp,rgas,grav,sgeo(1,jj),      &
                  pk(1,1,jj),pk2(1,1,jj),tt(1,1,jj),qt(1,1,jj),       &
                  phii(1,1,jj),phi(1,1,jj))
      do k = 1, lev
        do i = 1, nxj
          tr(i,k,jj)  = tt(i,k,jj)*pk(i,k,jj) / (1.0+0.608*qt(i,k,jj))
        enddo
      enddo
   end do

   do jj = 1, jlistnum
      do k = 1, lev
         kc = lev - k + 1
         do i = 1, myim(jj)
            qv3d (i,k,jj) = qt(i,             kc,jj)
            qc3d (i,k,jj) = qt(i,(ntcw-1)*lev+kc,jj)
            qr3d (i,k,jj) = qt(i,(ntrw-1)*lev+kc,jj)
            qi3d (i,k,jj) = qt(i,(ntiw-1)*lev+kc,jj)
            qs3d (i,k,jj) = qt(i,(ntsw-1)*lev+kc,jj)
            qg3d (i,k,jj) = qt(i,(ntgl-1)*lev+kc,jj)
            
            p3d  (i,k,jj) = 100.0*plt(i,kc,jj)                   !layer mean pressure (from mb to Pa)
            pii3d(i,k,jj) = pk(i,kc,jj)                          !exner function, =(plt/1000)**(Rd/cp)
            th3d (i,k,jj) = tr(i,kc,jj)/pk(i,kc,jj)                 !potential temperature (K)
            z3d  (i,k,jj) = phi(i,kc,jj)/con_g                   !layer geopotential height above sea level (m)
            dz3d (i,k,jj) = (phii(i,k+1,jj)-phii(i,k,jj))/con_g     !layer thickness (m)
            rho3d(i,k,jj) = p3d(i,k,jj)/                                  &
                           (con_rd*tr(i,kc,jj)*(1+con_fvirt*qv3d(i,k,jj)))  !air density (kg/m^3)
            w3d  (i,k,jj) = -vvel(i,kc,jj)*100./(rho3d(i,k,jj)*con_g)        !vertical velocity (m/s)
         
#ifdef Readaeroclx
           do n = 1, naero
              aero4d(i,k,n,jj) = aeroclxm(i,(n-1)*lev+kc,jj)
              aero(i, k, jj, n) = aero4d(i,k,n,jj)
           enddo
#endif      
         enddo
      enddo
   end do
      
   do jj = 1, jlistnum
      do i = 1, myim(jj)
         ht(i,jj) = sgeo(i,jj)/con_g   !terrain geopotential height above sea level (m)
         if( islimsk(i,jj) .eq. 1 ) then
            land2d(i,jj) = 1.        !land
            if ( ivegtyp(i,jj) .eq. 15 ) land2d(i,jj) = 2.   !glacial is seen as ocean
         else
            land2d(i,jj) = 2.        !ocean & seaice
         endif
      enddo
   end do
   
   
   

   rain2d_gpu = 0.
   ice2d_gpu = 0.
   snow2d_gpu = 0.
   graupel2d_gpu = 0.
   sr2d_gpu = 0.
   th3d_gpu = th3d
   qv3d_gpu = qv3d
   qc3d_gpu = qc3d
   qr3d_gpu = qr3d
   qi3d_gpu = qi3d
   qs3d_gpu = qs3d
   qg3d_gpu = qg3d
   rew3d_gpu = 0.
   rer3d_gpu = 0.
   rei3d_gpu = 0.
   res3d_gpu = 0.
   reg3d_gpu = 0.
   rain2d_cpu = 0.
   ice2d_cpu = 0.
   snow2d_cpu = 0.
   graupel2d_cpu = 0.
   sr2d_cpu = 0.
   th3d_cpu = th3d
   qv3d_cpu = qv3d
   qc3d_cpu = qc3d
   qr3d_cpu = qr3d
   qi3d_cpu = qi3d
   qs3d_cpu = qs3d
   qg3d_cpu = qg3d
   rew3d_cpu = 0.
   rer3d_cpu = 0.
   rei3d_cpu = 0.
   res3d_cpu = 0.
   reg3d_cpu = 0.
   
   call consat_s_gpu(ihail, 0, 3)
   call consat_s(ihail, 0, 3)
   call gce_table_copyin_gpu

   !if (benchmark .eq. .true.) call nvtxStartRange("CPU compute")
      do jj = 1, jlistnum
#ifdef Readaeroclx
         do k = 1, lev
            do i = 1, nxp
               do n = 1, naero
                   aero3d(i,k,1,n) = aero4d(i,k,n,jj)
               end do
            end do
         end do
#endif
      
   call cpu_time(time1)
         call fall_flux(dta, qv3d_cpu(1, 1, jj), qr3d_cpu(1, 1, jj), qi3d_cpu(1, 1, jj), &
                        qs3d_cpu(1, 1, jj), qg3d_cpu(1, 1, jj), p3d(1, 1, jj), &
                        rho3d(1, 1, jj), th3d(1, 1, jj), pii3d(1, 1, jj), &
                        z3d(1, 1, jj), dz3d(1, 1, jj), ht(1, jj), &
                        grav, itimestep, &
                        rain2d_cpu(1, jj), ice2d_cpu(1, jj), snow2d_cpu(1, jj), graupel2d_cpu(1, jj), sr2d_cpu(1, jj), &
                        !                      rainncv, snowncv, graupelncv,           &
#ifdef EXT_DIAG
                        preci3d(1, 1, jj), precs3d(1, 1, jj), &
                        precg3d(1, 1, jj), precr3d(1, 1, jj), &
#endif
                        land2d(1, jj), &
                        ihail, ice2, 3, &
                        1, nxp, 1, 1, 1, lev, & ! memory dims
                        1, myim(jj), 1, 1, 1, lev, & ! tile   dims
                        SL_sedi, sat_predict, new_saturation, benchmark) !flags
   call cpu_time(time2)
   !if (benchmark .eq. .true.) call nvtxEndRange
   if (benchmark .eq. .true.) ct = ct + time2 - time1
      end do


   if (.true.) then
      !if (benchmark .eq. .true.) call nvtxStartRange("GPU compute")
      !$acc enter data copyin(th3d_gpu, qv3d_gpu, qc3d_gpu, qr3d_gpu, &
      !$acc&      qi3d_gpu, qs3d_gpu, rain2d_gpu, ice2d_gpu, snow2d_gpu, &
      !$acc&      graupel2d_gpu, sr2d_gpu, qg3d_gpu, rew3d_gpu, rer3d_gpu, &
      !$acc&      rei3d_gpu, res3d_gpu, reg3d_gpu) async(async_id)
      !$acc enter data copyin(myim, rho3d, pii3d, p3d, z3d, ht, dz3d, &
      !$acc&      w3d, xlat_myim, land2d) async(async_id)
#ifdef Readaeroclx
      !$acc enter data copyin(aero4d) async(async_id)
#endif
      !$acc wait(async_id)
      call cpu_time(time1)
      call fall_flux_gpu(myim, dta, qv3d_gpu, qr3d_gpu, qi3d_gpu, &
                     qs3d_gpu, qg3d_gpu, p3d, &
                     rho3d, th3d, pii3d, &
                     z3d, dz3d, ht, &
                     grav, itimestep, &
                     rain2d_gpu, ice2d_gpu, snow2d_gpu, graupel2d_gpu, sr2d_gpu, &
                     !                      rainncv, snowncv, graupelncv,           &
#ifdef EXT_DIAG
                     preci3d, precs3d, &
                     precg3d, precr3d, &
#endif
                     land2d, &
                     ihail, ice2, 3, &
                     1, nxp, 1, jlistnum, 1, lev, & ! memory dims
                     1, nxp, 1, jlistnum, 1, lev, & ! tile   dims
                     SL_sedi, sat_predict, new_saturation, benchmark, async_id) !flags
      !sif (benchmark .eq. .true.) call nvtxEndRange
      !$acc wait(async_id)
      call cpu_time(time2)
      !$acc exit data copyout(th3d_gpu, qv3d_gpu, qc3d_gpu, qr3d_gpu, &
      !$acc&     qi3d_gpu, qs3d_gpu, rain2d_gpu, ice2d_gpu, snow2d_gpu, &
      !$acc&     graupel2d_gpu, sr2d_gpu, qg3d_gpu, rew3d_gpu, rer3d_gpu, &
      !$acc&     rei3d_gpu, res3d_gpu, reg3d_gpu) async(async_id)
      !$acc exit data delete(myim, rho3d, pii3d, p3d, z3d, ht, dz3d, &
      !$acc&     w3d, xlat_myim, land2d) async(async_id)
#ifdef Readaeroclx
      !$acc exit data delete(aero4d) async(async_id)
#endif
      !$acc wait(async_id)
      if (benchmark .eq. .true.) gt = gt + time2 - time1



   call assert_real(th3d_gpu, size(th3d_gpu), th3d_cpu, size(th3d_cpu), &
                    1e-8, "Array th3d"//trim(flag_str))
   call assert_real(qv3d_gpu, size(qv3d_gpu), qv3d_cpu, size(qv3d_cpu), &
                    1e-8, "Array qv3d"//trim(flag_str))
   call assert_real(qc3d_gpu, size(qc3d_gpu), qc3d_cpu, size(qc3d_cpu), &
                    1e-8, "Array qc3d"//trim(flag_str))
   call assert_real(qr3d_gpu, size(qr3d_gpu), qr3d_cpu, size(qr3d_cpu), &
                    1e-8, "Array qr3d"//trim(flag_str))
   call assert_real(qi3d_gpu, size(qi3d_gpu), qi3d_cpu, size(qi3d_cpu), &
                    1e-8, "Array qi3d"//trim(flag_str))
   call assert_real(qs3d_gpu, size(qs3d_gpu), qs3d_cpu, size(qs3d_cpu), &
                    1e-8, "Array qs3d"//trim(flag_str))
   call assert_real(qg3d_gpu, size(qg3d_gpu), qg3d_cpu, size(qg3d_cpu), &
                    1e-8, "Array qg3d"//trim(flag_str))
   call assert_real(rain2d_gpu, size(rain2d_gpu), rain2d_cpu, size(rain2d_cpu), &
                    1e-8, "Array rain2d"//trim(flag_str))
   call assert_real(ice2d_gpu, size(ice2d_gpu), ice2d_cpu, size(ice2d_cpu), &
                    1e-8, "Array ice2d"//trim(flag_str))
   call assert_real(snow2d_gpu, size(snow2d_gpu), snow2d_cpu, size(snow2d_cpu), &
                    1e-8, "Array snow2d"//trim(flag_str))
   call assert_real(graupel2d_gpu, size(graupel2d_gpu), graupel2d_cpu, size(graupel2d_cpu), &
                    1e-8, "Array graupel2d"//trim(flag_str))
   call assert_real(sr2d_gpu, size(sr2d_gpu), sr2d_cpu, size(sr2d_cpu), &
                    1e-8, "Array sr2d"//trim(flag_str))
   end if
   deallocate(flag_str)

  

end subroutine fall_flux_unit

subroutine saticel_s_unit(SL_sedi, sat_predict, new_saturation, &
                                   use_cpm, use_declination, benchmark, ct, gt)
   use module_mp_gsfcgce_3ice_nuwrf_gpu, only: gsfcgce_3ice_nuwrf_gpu, &
       consat_s_gpu, gce_table_copyin_gpu, saticel_s_gpu
   use module_mp_gsfcgce_3ice_nuwrf, only: gsfcgce_3ice_nuwrf, saticel_s, consat_s
   use phygrid, only: land, ice, ocean, xlat, aeroclxm
   use noah, only: ivegtyp
   use grid, only: sgeo, qt, pk, tt, phi, pk2, pt
   use index, only: nxjp, nxp, jlistnum, jlist1, nxdef_2d, nxjp_acc
   use param, only: lev, my_max, my, ncld
   use machine, only: kind_phys
   use rank, only: myrank
   use const, only: dsigma, RTYPE, naero, sigma
   use radn,                only: ntcw,ntiw,ntrw,ntsw,ntgl,nthl,     &
                                  ntinc,ntrnc
   use physcons, only: con_rd, con_fvirt, con_g

   !use nvtx

   implicit none

   integer :: ntrac, check, myim(my_max)
   integer :: i, ii, j, jj, async_id, n, k, iii, cnt, tcnt, nxj, kc
   integer, dimension(34) :: seed
   real :: time1, time2
   character(len=:),allocatable :: flag_str
   !flags
   logical, intent(in) :: SL_sedi, sat_predict, new_saturation, use_cpm, use_declination
   logical, intent(in) :: benchmark
   real, intent(inout) :: ct, gt
   
   integer :: itimestep, ihail, ice2, islimsk(nxp, my_max)
   real :: dta, sdec, grav, rgas, cp, ptop
   real, dimension(my_max) :: xlat_myim
   real(kind=kind_phys), dimension(nxp, lev+1, my_max) :: phii
   real,dimension(nxp, my_max) :: rain2d, ice2d, snow2d, &
      graupel2d, sr2d, land2d, ht
   real,dimension(nxp, my_max) :: rain2d_cpu, ice2d_cpu, snow2d_cpu, &
      graupel2d_cpu, sr2d_cpu
   real,dimension(nxp, my_max) :: rain2d_gpu, ice2d_gpu, snow2d_gpu, &
      graupel2d_gpu, sr2d_gpu
   real,dimension(nxp, lev, my_max) :: th3d, qv3d, &
      qc3d, qr3d, qi3d, qs3d, qg3d, rew3d, rer3d, &
      rei3d, res3d, reg3d, rho3d, pii3d, p3d, z3d, dz3d, w3d, plt, vvel, tr
   real,dimension(nxp, lev, my_max) :: th3d_cpu, qv3d_cpu, &
      qc3d_cpu, qr3d_cpu, qi3d_cpu, qs3d_cpu, qg3d_cpu, rew3d_cpu, rer3d_cpu, &
      rei3d_cpu, res3d_cpu, reg3d_cpu
   real,dimension(nxp, lev, my_max) :: th3d_gpu, qv3d_gpu, &
      qc3d_gpu, qr3d_gpu, qi3d_gpu, qs3d_gpu, qg3d_gpu, rew3d_gpu, rer3d_gpu, &
      rei3d_gpu, res3d_gpu, reg3d_gpu
   real, dimension(nxp, lev, 1, naero) :: aero3d
   real, dimension(nxp, lev, naero, my_max) :: aero4d
   real, dimension(nxp, lev, my_max, naero) :: aero
   
   flag_str = ''
   if (SL_sedi .eq. .true.) flag_str = flag_str//', SL_sedi'
   if (sat_predict .eq. .true.) flag_str = flag_str//', sat_predict'
   if (new_saturation .eq. .true.) flag_str = flag_str//', new_saturation'
   if (use_cpm .eq. .true.) flag_str = flag_str//', use_cpm'
   if (use_declination .eq. .true.) flag_str = flag_str//', use_declination'
   if (myrank .eq. 0) write(*,*) flag_str
   
   async_id = 1
   seed = (/10006, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,&
           &25, 26, 27, 28, 29, 30, 31, 32, 33/)
   !call random_seed()
   call random_seed(put=seed)
   call random_number(vvel)
   vvel = vvel * 4e-5
   
   cnt = 0
   tcnt = 0

   grav = 9.806649999999999
   rgas = 287.0285714285714
   cp = 1004.600000000000
   ptop = 0.1000000000000000
   sdec = 0.
   ihail = 0
   ice2 = 0
   dta = 600.
   itimestep = 25
   
   do jj = 1, jlistnum
      j = jlist1(jj)
      myim(jj) = nxjp(j)
      xlat_myim(jj) = xlat(j)
   end do
   
   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef_2d(j)
      do i = 1, nxj
         if (land(i, jj)) islimsk(i, jj) = 1
         if (ocean(i, jj)) islimsk(i, jj) = 0
         if (ice(i, jj)) islimsk(i, jj) = 2
      end do
   end do
   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef_2d(j)
      !
      !   new p**capa quantities were computed in previous diabat call
      !
      call prexp_hybrid_cwb(nxjp(j), nxp, lev, ptop, sigma, pt(1, jj), &
                          pk(1, 1, jj), pk2(1, 1, jj), plt(1, 1, jj))
   end do !jj = 1,jlistnum


   do jj = 1, jlistnum
      j=jlist1(jj)
      nxj=nxdef_2d(j)
      call get_phi(nxjp(j),nxp,lev,ptop,cp,rgas,grav,sgeo(1,jj),      &
                  pk(1,1,jj),pk2(1,1,jj),tt(1,1,jj),qt(1,1,jj),       &
                  phii(1,1,jj),phi(1,1,jj))
      do k = 1, lev
        do i = 1, nxj
          tr(i,k,jj)  = tt(i,k,jj)*pk(i,k,jj) / (1.0+0.608*qt(i,k,jj))
        enddo
      enddo
   end do

   do jj = 1, jlistnum
      do k = 1, lev
         kc = lev - k + 1
         do i = 1, myim(jj)
            qv3d (i,k,jj) = qt(i,             kc,jj)
            qc3d (i,k,jj) = qt(i,(ntcw-1)*lev+kc,jj)
            qr3d (i,k,jj) = qt(i,(ntrw-1)*lev+kc,jj)
            qi3d (i,k,jj) = qt(i,(ntiw-1)*lev+kc,jj)
            qs3d (i,k,jj) = qt(i,(ntsw-1)*lev+kc,jj)
            qg3d (i,k,jj) = qt(i,(ntgl-1)*lev+kc,jj)
            
            p3d  (i,k,jj) = 100.0*plt(i,kc,jj)                   !layer mean pressure (from mb to Pa)
            pii3d(i,k,jj) = pk(i,kc,jj)                          !exner function, =(plt/1000)**(Rd/cp)
            th3d (i,k,jj) = tr(i,kc,jj)/pk(i,kc,jj)                 !potential temperature (K)
            z3d  (i,k,jj) = phi(i,kc,jj)/con_g                   !layer geopotential height above sea level (m)
            dz3d (i,k,jj) = (phii(i,k+1,jj)-phii(i,k,jj))/con_g     !layer thickness (m)
            rho3d(i,k,jj) = p3d(i,k,jj)/                                  &
                           (con_rd*tr(i,kc,jj)*(1+con_fvirt*qv3d(i,k,jj)))  !air density (kg/m^3)
            w3d  (i,k,jj) = -vvel(i,kc,jj)*100./(rho3d(i,k,jj)*con_g)        !vertical velocity (m/s)
         
#ifdef Readaeroclx
           do n = 1, naero
              aero4d(i,k,n,jj) = aeroclxm(i,(n-1)*lev+kc,jj)
              aero(i, k, jj, n) = aero4d(i,k,n,jj)
           enddo
#endif      
         enddo
      enddo
   end do
      
   do jj = 1, jlistnum
      do i = 1, myim(jj)
         ht(i,jj) = sgeo(i,jj)/con_g   !terrain geopotential height above sea level (m)
         if( islimsk(i,jj) .eq. 1 ) then
            land2d(i,jj) = 1.        !land
            if ( ivegtyp(i,jj) .eq. 15 ) land2d(i,jj) = 2.   !glacial is seen as ocean
         else
            land2d(i,jj) = 2.        !ocean & seaice
         endif
      enddo
   end do
   
   
   

   rain2d_gpu = 0.
   ice2d_gpu = 0.
   snow2d_gpu = 0.
   graupel2d_gpu = 0.
   sr2d_gpu = 0.
   th3d_gpu = th3d
   qv3d_gpu = qv3d
   qc3d_gpu = qc3d
   qr3d_gpu = qr3d
   qi3d_gpu = qi3d
   qs3d_gpu = qs3d
   qg3d_gpu = qg3d
   rew3d_gpu = 0.
   rer3d_gpu = 0.
   rei3d_gpu = 0.
   res3d_gpu = 0.
   reg3d_gpu = 0.
   rain2d_cpu = 0.
   ice2d_cpu = 0.
   snow2d_cpu = 0.
   graupel2d_cpu = 0.
   sr2d_cpu = 0.
   th3d_cpu = th3d
   qv3d_cpu = qv3d
   qc3d_cpu = qc3d
   qr3d_cpu = qr3d
   qi3d_cpu = qi3d
   qs3d_cpu = qs3d
   qg3d_cpu = qg3d
   rew3d_cpu = 0.
   rer3d_cpu = 0.
   rei3d_cpu = 0.
   res3d_cpu = 0.
   reg3d_cpu = 0.
   
   call consat_s_gpu(ihail, 0, 3)
   call consat_s(ihail, 0, 3)
   call gce_table_copyin_gpu

   !if (benchmark .eq. .true.) call nvtxStartRange("CPU compute")
      do jj = 1, jlistnum
#ifdef Readaeroclx
         do k = 1, lev
            do i = 1, nxp
               do n = 1, naero
                   aero3d(i,k,1,n) = aero4d(i,k,n,jj)
               end do
            end do
         end do
#endif
      
   call cpu_time(time1)
         call SATICEL_S(dta, IHAIL, 0, ICE2, 180, &
                        3, 0, 3, xlat(jj), sdec, &
                        th3d_cpu(1, 1, jj), qv3d_cpu(1, 1, jj), qc3d_cpu(1, 1, jj), qr3d_cpu(1, 1, jj), &
                        qi3d_cpu(1, 1, jj), qs3d_cpu(1, 1, jj), qg3d_cpu(1, 1, jj), &
                        rho3d(1, 1, jj), pii3d(1, 1, jj), p3d(1, 1, jj), w3d(1, 1, jj), &
                        itimestep, land2d(1, jj), &
                        rew3d_cpu(1, 1, jj), rer3d_cpu(1, 1, jj), rei3d_cpu(1, 1, jj), &
                        res3d_cpu(1, 1, jj), reg3d_cpu(1, 1, jj), & ! cloud effective radius
                         1, nxp , 1, 1, 1, lev,         & ! memory dims
                         1, myim(jj), 1, 1, 1, lev,     & ! tile   dims
#ifdef Readaeroclx
                        aero3d, naero, &
#endif
#ifdef EXT_DIAG
                        refl_10cm, diagflag, do_radar_ref, & ! GT added for reflectivity calcs
                        physc, physe, physd, physs, physm, physf, &
                        acphysc, acphyse, acphysd, acphyss, acphysm, acphysf, &
#endif
                        sat_predict, new_saturation, use_cpm, use_declination, &
                        benchmark) !flagis
   call cpu_time(time2)
   !if (benchmark .eq. .true.) call nvtxEndRange
   if (benchmark .eq. .true.) ct = ct + time2 - time1
      end do


   if (.true.) then
      !if (benchmark .eq. .true.) call nvtxStartRange("GPU compute")
      !$acc enter data copyin(th3d_gpu, qv3d_gpu, qc3d_gpu, qr3d_gpu, &
      !$acc&      qi3d_gpu, qs3d_gpu, rain2d_gpu, ice2d_gpu, snow2d_gpu, &
      !$acc&      graupel2d_gpu, sr2d_gpu, qg3d_gpu, rew3d_gpu, rer3d_gpu, &
      !$acc&      rei3d_gpu, res3d_gpu, reg3d_gpu) async(async_id)
      !$acc enter data copyin(myim, rho3d, pii3d, p3d, z3d, ht, dz3d, &
      !$acc&      w3d, xlat_myim, land2d) async(async_id)
#ifdef Readaeroclx
      !$acc enter data copyin(aero4d) async(async_id)
#endif
      !$acc wait(async_id)
      call cpu_time(time1)
      call SATICEL_S_gpu(myim, dta, IHAIL, 0, ICE2, 180, &
                     3, 0, 3, xlat, sdec, &
                     th3d_gpu, qv3d_gpu, qc3d_gpu, qr3d_gpu, &
                     qi3d_gpu, qs3d_gpu, qg3d_gpu, &
                     rho3d, pii3d, p3d, w3d, &
                     itimestep, land2d, &
                     rew3d_gpu, rer3d_gpu, rei3d_gpu, &
                     res3d_gpu, reg3d_gpu, & ! cloud effective radius
                     1, nxp, 1, jlistnum, 1, lev, & ! memory dims
                     1, nxp, 1, jlistnum, 1, lev, & ! tile   dims
#ifdef Readaeroclx
                     aero4d, naero, &
#endif
#ifdef EXT_DIAG
                     refl_10cm, diagflag, do_radar_ref, & ! GT added for reflectivity calcs
                     physc, physe, physd, physs, physm, physf, &
                     acphysc, acphyse, acphysd, acphyss, acphysm, acphysf, &
#endif
                     sat_predict, new_saturation, use_cpm, use_declination, &
                     benchmark, async_id) !flags
      !sif (benchmark .eq. .true.) call nvtxEndRange
      !$acc wait(async_id)
      call cpu_time(time2)
      !$acc exit data copyout(th3d_gpu, qv3d_gpu, qc3d_gpu, qr3d_gpu, &
      !$acc&     qi3d_gpu, qs3d_gpu, rain2d_gpu, ice2d_gpu, snow2d_gpu, &
      !$acc&     graupel2d_gpu, sr2d_gpu, qg3d_gpu, rew3d_gpu, rer3d_gpu, &
      !$acc&     rei3d_gpu, res3d_gpu, reg3d_gpu) async(async_id)
      !$acc exit data delete(myim, rho3d, pii3d, p3d, z3d, ht, dz3d, &
      !$acc&     w3d, xlat_myim, land2d) async(async_id)
#ifdef Readaeroclx
      !$acc exit data delete(aero4d) async(async_id)
#endif
      !$acc wait(async_id)
      if (benchmark .eq. .true.) gt = gt + time2 - time1


   call check_real(qv3d_gpu, size(qv3d_gpu), qv3d_cpu, size(qv3d_cpu), &
                   1e-6, 200, "Check qv3d"//trim(flag_str))
   call check_real(qc3d_gpu, size(qc3d_gpu), qc3d_cpu, size(qc3d_cpu), &
                   1e-6, 200, "Check qc3d"//trim(flag_str))
   call check_real(qs3d_gpu, size(qs3d_gpu), qs3d_cpu, size(qs3d_cpu), &
                   1e-6, 200, "Check qs3d"//trim(flag_str))
   call check_real(qr3d_gpu, size(qr3d_gpu), qr3d_cpu, size(qr3d_cpu), &
                   1e-6, 200, "Check qr3d"//trim(flag_str))
   call check_real(qi3d_gpu, size(qi3d_gpu), qi3d_cpu, size(qi3d_cpu), &
                   1e-6, 200, "Check qi3d"//trim(flag_str))
   call check_real(qg3d_gpu, size(qg3d_gpu), qg3d_cpu, size(qg3d_cpu), &
                   1e-6, 200, "Check qg3d"//trim(flag_str))

   call assert_real(th3d_gpu, size(th3d_gpu), th3d_cpu, size(th3d_cpu), &
                    1e-8, "Array th3d"//trim(flag_str))
   call assert_real(qv3d_gpu, size(qv3d_gpu), qv3d_cpu, size(qv3d_cpu), &
                    1e-8, "Array qv3d"//trim(flag_str))
   !call assert_real(qc3d_gpu, size(qc3d_gpu), qc3d_cpu, size(qc3d_cpu), &
   !                 1e-8, "Array qc3d"//trim(flag_str))
   call assert_real(qr3d_gpu, size(qr3d_gpu), qr3d_cpu, size(qr3d_cpu), &
                    1e-6, "Array qr3d"//trim(flag_str))
   !call assert_real(qi3d_gpu, size(qi3d_gpu), qi3d_cpu, size(qi3d_cpu), &
   !                 1e-8, "Array qi3d"//trim(flag_str))
   call assert_real(qs3d_gpu, size(qs3d_gpu), qs3d_cpu, size(qs3d_cpu), &
                    1e-8, "Array qs3d"//trim(flag_str))
   call assert_real(qg3d_gpu, size(qg3d_gpu), qg3d_cpu, size(qg3d_cpu), &
                    1e-8, "Array qg3d"//trim(flag_str))
   call assert_real(rew3d_gpu, size(rew3d_gpu), rew3d_cpu, size(rew3d_cpu), &
                    1e-6, "Array rew3d"//trim(flag_str))
   call assert_real(rer3d_gpu, size(rer3d_gpu), rer3d_cpu, size(rer3d_cpu), &
                    1e-8, "Array rer3d"//trim(flag_str))
   call assert_real(rei3d_gpu, size(rei3d_gpu), rei3d_cpu, size(rei3d_cpu), &
                    1e-8, "Array rei3d"//trim(flag_str))
   call assert_real(res3d_gpu, size(res3d_gpu), res3d_cpu, size(res3d_cpu), &
                    1e-8, "Array res3d"//trim(flag_str))
   call assert_real(reg3d_gpu, size(reg3d_gpu), reg3d_cpu, size(reg3d_cpu), &
                    1e-8, "Array reg3d"//trim(flag_str))
   end if
   deallocate(flag_str)

  

end subroutine saticel_s_unit


subroutine gsfcgce_3ice_nuwrf_unit(SL_sedi, sat_predict, new_saturation, &
                                   use_cpm, use_declination, benchmark, ct, gt)
   use module_mp_gsfcgce_3ice_nuwrf_gpu, only: gsfcgce_3ice_nuwrf_gpu, &
       consat_s_gpu, gce_table_copyin_gpu
   use module_mp_gsfcgce_3ice_nuwrf, only: gsfcgce_3ice_nuwrf
   use phygrid, only: land, ice, ocean, xlat, aeroclxm
   use noah, only: ivegtyp
   use grid, only: sgeo, qt, pk, tt, phi, pk2, pt
   use index, only: nxjp, nxp, jlistnum, jlist1, nxdef_2d, nxjp_acc
   use param, only: lev, my_max, my, ncld
   use machine, only: kind_phys
   use rank, only: myrank
   use const, only: dsigma, RTYPE, naero, sigma
   use radn,                only: ntcw,ntiw,ntrw,ntsw,ntgl,nthl,     &
                                  ntinc,ntrnc
   use physcons, only: con_rd, con_fvirt, con_g

   !use nvtx

   implicit none

   integer :: ntrac, check, myim(my_max)
   integer :: i, ii, j, jj, async_id, n, k, iii, cnt, tcnt, nxj, kc
   integer, dimension(34) :: seed
   real :: time1, time2
   character(len=:),allocatable :: flag_str
   !flags
   logical, intent(in) :: SL_sedi, sat_predict, new_saturation, use_cpm, use_declination
   logical, intent(in) :: benchmark
   real, intent(inout) :: ct, gt
   
   integer :: itimestep, ihail, ice2, islimsk(nxp, my_max)
   real :: dta, sdec, grav, rgas, cp, ptop
   real, dimension(my_max) :: xlat_myim
   real(kind=kind_phys), dimension(nxp, lev+1, my_max) :: phii
   real,dimension(nxp, my_max) :: rain2d, ice2d, snow2d, &
      graupel2d, sr2d, land2d, ht
   real,dimension(nxp, my_max) :: rain2d_cpu, ice2d_cpu, snow2d_cpu, &
      graupel2d_cpu, sr2d_cpu
   real,dimension(nxp, my_max) :: rain2d_gpu, ice2d_gpu, snow2d_gpu, &
      graupel2d_gpu, sr2d_gpu
   real,dimension(nxp, lev, my_max) :: th3d, qv3d, &
      qc3d, qr3d, qi3d, qs3d, qg3d, rew3d, rer3d, &
      rei3d, res3d, reg3d, rho3d, pii3d, p3d, z3d, dz3d, w3d, plt, vvel, tr
   real,dimension(nxp, lev, my_max) :: th3d_cpu, qv3d_cpu, &
      qc3d_cpu, qr3d_cpu, qi3d_cpu, qs3d_cpu, qg3d_cpu, rew3d_cpu, rer3d_cpu, &
      rei3d_cpu, res3d_cpu, reg3d_cpu
   real,dimension(nxp, lev, my_max) :: th3d_gpu, qv3d_gpu, &
      qc3d_gpu, qr3d_gpu, qi3d_gpu, qs3d_gpu, qg3d_gpu, rew3d_gpu, rer3d_gpu, &
      rei3d_gpu, res3d_gpu, reg3d_gpu
   real, dimension(nxp, lev, 1, naero) :: aero3d
   real, dimension(nxp, lev, naero, my_max) :: aero4d
   real, dimension(nxp, lev, my_max, naero) :: aero
   
   flag_str = ''
   if (SL_sedi .eq. .true.) flag_str = flag_str//', SL_sedi'
   if (sat_predict .eq. .true.) flag_str = flag_str//', sat_predict'
   if (new_saturation .eq. .true.) flag_str = flag_str//', new_saturation'
   if (use_cpm .eq. .true.) flag_str = flag_str//', use_cpm'
   if (use_declination .eq. .true.) flag_str = flag_str//', use_declination'
   if (myrank .eq. 0) write(*,*) flag_str
   
   async_id = 1
   seed = (/10006, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,&
           &25, 26, 27, 28, 29, 30, 31, 32, 33/)
   !call random_seed()
   call random_seed(put=seed)
   call random_number(vvel)
   vvel = vvel * 4e-5
   
   cnt = 0
   tcnt = 0

   grav = 9.806649999999999
   rgas = 287.0285714285714
   cp = 1004.600000000000
   ptop = 0.1000000000000000
   sdec = 0.
   ihail = 0
   ice2 = 0
   dta = 600.
   itimestep = 25
   
   do jj = 1, jlistnum
      j = jlist1(jj)
      myim(jj) = nxjp(j)
      xlat_myim(jj) = xlat(j)
   end do
   
   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef_2d(j)
      do i = 1, nxj
         if (land(i, jj)) islimsk(i, jj) = 1
         if (ocean(i, jj)) islimsk(i, jj) = 0
         if (ice(i, jj)) islimsk(i, jj) = 2
      end do
   end do
   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef_2d(j)
      !
      !   new p**capa quantities were computed in previous diabat call
      !
      call prexp_hybrid_cwb(nxjp(j), nxp, lev, ptop, sigma, pt(1, jj), &
                          pk(1, 1, jj), pk2(1, 1, jj), plt(1, 1, jj))
   end do !jj = 1,jlistnum


   do jj = 1, jlistnum
      j=jlist1(jj)
      nxj=nxdef_2d(j)
      call get_phi(nxjp(j),nxp,lev,ptop,cp,rgas,grav,sgeo(1,jj),      &
                  pk(1,1,jj),pk2(1,1,jj),tt(1,1,jj),qt(1,1,jj),       &
                  phii(1,1,jj),phi(1,1,jj))
      do k = 1, lev
        do i = 1, nxj
          tr(i,k,jj)  = tt(i,k,jj)*pk(i,k,jj) / (1.0+0.608*qt(i,k,jj))
        enddo
      enddo
   end do

   do jj = 1, jlistnum
      do k = 1, lev
         kc = lev - k + 1
         do i = 1, myim(jj)
            qv3d (i,k,jj) = qt(i,             kc,jj)
            qc3d (i,k,jj) = qt(i,(ntcw-1)*lev+kc,jj)
            qr3d (i,k,jj) = qt(i,(ntrw-1)*lev+kc,jj)
            qi3d (i,k,jj) = qt(i,(ntiw-1)*lev+kc,jj)
            qs3d (i,k,jj) = qt(i,(ntsw-1)*lev+kc,jj)
            qg3d (i,k,jj) = qt(i,(ntgl-1)*lev+kc,jj)
            
            p3d  (i,k,jj) = 100.0*plt(i,kc,jj)                   !layer mean pressure (from mb to Pa)
            pii3d(i,k,jj) = pk(i,kc,jj)                          !exner function, =(plt/1000)**(Rd/cp)
            th3d (i,k,jj) = tr(i,kc,jj)/pk(i,kc,jj)                 !potential temperature (K)
            z3d  (i,k,jj) = phi(i,kc,jj)/con_g                   !layer geopotential height above sea level (m)
            dz3d (i,k,jj) = (phii(i,k+1,jj)-phii(i,k,jj))/con_g     !layer thickness (m)
            rho3d(i,k,jj) = p3d(i,k,jj)/                                  &
                           (con_rd*tr(i,kc,jj)*(1+con_fvirt*qv3d(i,k,jj)))  !air density (kg/m^3)
            w3d  (i,k,jj) = -vvel(i,kc,jj)*100./(rho3d(i,k,jj)*con_g)        !vertical velocity (m/s)
         
#ifdef Readaeroclx
           do n = 1, naero
              aero4d(i,k,n,jj) = aeroclxm(i,(n-1)*lev+kc,jj)
              aero(i, k, jj, n) = aero4d(i,k,n,jj)
           enddo
#endif      
         enddo
      enddo
   end do
      
   do jj = 1, jlistnum
      do i = 1, myim(jj)
         ht(i,jj) = sgeo(i,jj)/con_g   !terrain geopotential height above sea level (m)
         if( islimsk(i,jj) .eq. 1 ) then
            land2d(i,jj) = 1.        !land
            if ( ivegtyp(i,jj) .eq. 15 ) land2d(i,jj) = 2.   !glacial is seen as ocean
         else
            land2d(i,jj) = 2.        !ocean & seaice
         endif
      enddo
   end do
   
   
   

   rain2d_gpu = 0.
   ice2d_gpu = 0.
   snow2d_gpu = 0.
   graupel2d_gpu = 0.
   sr2d_gpu = 0.
   th3d_gpu = th3d
   qv3d_gpu = qv3d
   qc3d_gpu = qc3d
   qr3d_gpu = qr3d
   qi3d_gpu = qi3d
   qs3d_gpu = qs3d
   qg3d_gpu = qg3d
   rew3d_gpu = 0.
   rer3d_gpu = 0.
   rei3d_gpu = 0.
   res3d_gpu = 0.
   reg3d_gpu = 0.
   rain2d_cpu = 0.
   ice2d_cpu = 0.
   snow2d_cpu = 0.
   graupel2d_cpu = 0.
   sr2d_cpu = 0.
   th3d_cpu = th3d
   qv3d_cpu = qv3d
   qc3d_cpu = qc3d
   qr3d_cpu = qr3d
   qi3d_cpu = qi3d
   qs3d_cpu = qs3d
   qg3d_cpu = qg3d
   rew3d_cpu = 0.
   rer3d_cpu = 0.
   rei3d_cpu = 0.
   res3d_cpu = 0.
   reg3d_cpu = 0.
   
   call consat_s_gpu(ihail, 0, 3)
   call gce_table_copyin_gpu

   !if (benchmark .eq. .true.) call nvtxStartRange("CPU compute")
      do jj = 1, jlistnum
#ifdef Readaeroclx
         do k = 1, lev
            do i = 1, nxp
               do n = 1, naero
                   aero3d(i,k,1,n) = aero4d(i,k,n,jj)
               end do
            end do
         end do
#endif
      
   call cpu_time(time1)
         call gsfcgce_3ice_nuwrf                                       &
                ( th3d_cpu(1,1,jj), qv3d_cpu(1,1,jj), qc3d_cpu(1,1,jj), qr3d_cpu(1,1,jj), &
                  qi3d_cpu(1,1,jj), qs3d_cpu(1,1,jj),                  &
                  rho3d(1,1,jj), pii3d(1,1,jj), p3d(1,1,jj), dta, z3d(1,1,jj),                         &
                  ht(1,jj), dz3d(1,1,jj), con_g, w3d(1,1,jj),                                &
                  itimestep, xlat_myim(jj), sdec, land2d(1,jj),                       &
                  1, nxp , 1, 1, 1, lev,                                & ! memory dims
                  1, myim(jj), 1, 1, 1, lev,                                & ! tile   dims
                  rain2d_cpu(1,jj), ice2d_cpu(1,jj), snow2d_cpu(1,jj), graupel2d_cpu(1,jj), &
                  sr2d_cpu(1,jj),              &
                  .false., qg3d_cpu(1,1,jj),                                       &
#ifdef Readaeroclx
                  aero3d, naero,                                       &
#endif      
                  ihail, ice2,                                         &
#ifdef EXT_DIAG
                  refl_10cm(1,1,jj), diagflag, do_radar_ref,                   &
                  preci3d(1,1,jj), precs3d(1,1,jj), precg3d(1,1,jj), precr3d(1,1,jj),                  &
#endif      
                  rew3d_cpu(1,1,jj), rer3d_cpu(1,1,jj), rei3d_cpu(1,1,jj), res3d_cpu(1,1,jj), &
                  reg3d_cpu(1,1,jj), &
                  SL_sedi, sat_predict, new_saturation, use_cpm, use_declination, &
                  benchmark)
   call cpu_time(time2)
   !if (benchmark .eq. .true.) call nvtxEndRange
   if (benchmark .eq. .true.) ct = ct + time2 - time1
      end do


   if (.true.) then
      !if (benchmark .eq. .true.) call nvtxStartRange("GPU compute")
      !$acc enter data copyin(th3d_gpu, qv3d_gpu, qc3d_gpu, qr3d_gpu, &
      !$acc&      qi3d_gpu, qs3d_gpu, rain2d_gpu, ice2d_gpu, snow2d_gpu, &
      !$acc&      graupel2d_gpu, sr2d_gpu, qg3d_gpu, rew3d_gpu, rer3d_gpu, &
      !$acc&      rei3d_gpu, res3d_gpu, reg3d_gpu) async(async_id)
      !$acc enter data copyin(myim, rho3d, pii3d, p3d, z3d, ht, dz3d, &
      !$acc&      w3d, xlat_myim, land2d) async(async_id)
#ifdef Readaeroclx
      !$acc enter data copyin(aero4d) async(async_id)
#endif
      !$acc wait(async_id)
      call cpu_time(time1)
      call gsfcgce_3ice_nuwrf_gpu                                       &
             ( myim, th3d_gpu, qv3d_gpu, qc3d_gpu, qr3d_gpu, qi3d_gpu, qs3d_gpu,                  &
               rho3d, pii3d, p3d, dta, z3d,                         &
               ht, dz3d, con_g, w3d,                                &
               itimestep, xlat_myim, sdec, land2d,                       &
               1, nxp , 1, jlistnum, 1, lev,                                & ! memory dims
               1, nxp, 1, jlistnum, 1, lev,                                & ! tile   dims
               rain2d_gpu, ice2d_gpu, snow2d_gpu, graupel2d_gpu, sr2d_gpu,              &
               .false., qg3d_gpu,                                       &
#ifdef Readaeroclx
               aero4d, naero,                                       &
#endif      
               ihail, ice2,                                         &
#ifdef EXT_DIAG
               refl_10cm_gpu, diagflag, do_radar_ref,                   &
               preci3d_gpu, precs3d_gpu, precg3d_gpu, precr3d_gpu,                  &
#endif      
               rew3d_gpu, rer3d_gpu, rei3d_gpu, res3d_gpu, reg3d_gpu, &
               SL_sedi, sat_predict, new_saturation, use_cpm, use_declination, &
               benchmark, async_id)
      !sif (benchmark .eq. .true.) call nvtxEndRange
      !$acc wait(async_id)
      call cpu_time(time2)
      !$acc exit data copyout(th3d_gpu, qv3d_gpu, qc3d_gpu, qr3d_gpu, &
      !$acc&     qi3d_gpu, qs3d_gpu, rain2d_gpu, ice2d_gpu, snow2d_gpu, &
      !$acc&     graupel2d_gpu, sr2d_gpu, qg3d_gpu, rew3d_gpu, rer3d_gpu, &
      !$acc&     rei3d_gpu, res3d_gpu, reg3d_gpu) async(async_id)
      !$acc exit data delete(myim, rho3d, pii3d, p3d, z3d, ht, dz3d, &
      !$acc&     w3d, xlat_myim, land2d) async(async_id)
#ifdef Readaeroclx
      !$acc exit data delete(aero4d) async(async_id)
#endif
      !$acc wait(async_id)
      if (benchmark .eq. .true.) gt = gt + time2 - time1


   call check_real(qv3d_gpu, size(qv3d_gpu), qv3d_cpu, size(qv3d_cpu), &
                   1e-6, 200, "Check qv3d"//trim(flag_str))
   call check_real(qc3d_gpu, size(qc3d_gpu), qc3d_cpu, size(qc3d_cpu), &
                   1e-6, 200, "Check qc3d"//trim(flag_str))
   call check_real(qs3d_gpu, size(qs3d_gpu), qs3d_cpu, size(qs3d_cpu), &
                   1e-6, 200, "Check qs3d"//trim(flag_str))
   call check_real(qr3d_gpu, size(qr3d_gpu), qr3d_cpu, size(qr3d_cpu), &
                   1e-6, 200, "Check qr3d"//trim(flag_str))
   call check_real(qi3d_gpu, size(qi3d_gpu), qi3d_cpu, size(qi3d_cpu), &
                   1e-6, 200, "Check qi3d"//trim(flag_str))
   call check_real(qg3d_gpu, size(qg3d_gpu), qg3d_cpu, size(qg3d_cpu), &
                   1e-6, 200, "Check qg3d"//trim(flag_str))

   call assert_real(th3d_gpu, size(th3d_gpu), th3d_cpu, size(th3d_cpu), &
                    1e-8, "Array th3d"//trim(flag_str))
   call assert_real(qv3d_gpu, size(qv3d_gpu), qv3d_cpu, size(qv3d_cpu), &
                    1e-6, "Array qv3d"//trim(flag_str))
   call assert_real(qc3d_gpu, size(qc3d_gpu), qc3d_cpu, size(qc3d_cpu), &
                    1e-4, "Array qc3d"//trim(flag_str))
   call assert_real(qr3d_gpu, size(qr3d_gpu), qr3d_cpu, size(qr3d_cpu), &
                    1e-6, "Array qr3d"//trim(flag_str))
   !call assert_real(qi3d_gpu, size(qi3d_gpu), qi3d_cpu, size(qi3d_cpu), &
   !                 1e-8, "Array qi3d"//trim(flag_str))
   !call assert_real(qs3d_gpu, size(qs3d_gpu), qs3d_cpu, size(qs3d_cpu), &
   !                 1e-8, "Array qs3d"//trim(flag_str))
   call assert_real(rain2d_gpu, size(rain2d_gpu), rain2d_cpu, size(rain2d_cpu), &
                    1e-8, "Array rain2d"//trim(flag_str))
   call assert_real(ice2d_gpu, size(ice2d_gpu), ice2d_cpu, size(ice2d_cpu), &
                    1e-8, "Array ice2d"//trim(flag_str))
   call assert_real(snow2d_gpu, size(snow2d_gpu), snow2d_cpu, size(snow2d_cpu), &
                    1e-8, "Array snow2d"//trim(flag_str))
   call assert_real(graupel2d_gpu, size(graupel2d_gpu), graupel2d_cpu, size(graupel2d_cpu), &
                    1e-8, "Array graupel2d"//trim(flag_str))
   call assert_real(sr2d_gpu, size(sr2d_gpu), sr2d_cpu, size(sr2d_cpu), &
                    1e-8, "Array sr2d"//trim(flag_str))
   call assert_real(qg3d_gpu, size(qg3d_gpu), qg3d_cpu, size(qg3d_cpu), &
                    1e-8, "Array qg3d"//trim(flag_str))
   call assert_real(rew3d_gpu, size(rew3d_gpu), rew3d_cpu, size(rew3d_cpu), &
                    1e-4, "Array rew3d"//trim(flag_str))
   call assert_real(rer3d_gpu, size(rer3d_gpu), rer3d_cpu, size(rer3d_cpu), &
                    1e-8, "Array rer3d"//trim(flag_str))
   call assert_real(rei3d_gpu, size(rei3d_gpu), rei3d_cpu, size(rei3d_cpu), &
                    1e-1, "Array rei3d"//trim(flag_str))
   call assert_real(res3d_gpu, size(res3d_gpu), res3d_cpu, size(res3d_cpu), &
                    1e-8, "Array res3d"//trim(flag_str))
   call assert_real(reg3d_gpu, size(reg3d_gpu), reg3d_cpu, size(reg3d_cpu), &
                    1e-8, "Array reg3d"//trim(flag_str))
   end if
   deallocate(flag_str)

   

end subroutine gsfcgce_3ice_nuwrf_unit

program test_gsfcgce_3ice_nuwrf
   use index, only: nxp
   use rank, only: myrank
   use param
   use const
   use module_gocart_coupling, only : makelut_ccn_icn

   implicit none
   integer :: no
   real :: ct, gt
   
   
   
   call mpe_init
   call cons
   call makelut_ccn_icn
   call getrdy
   if (donnmi .and. taui .lt. 1.0) then
      no = 2*((jtrun + 1)/2) + (jtrun/2) + 10
      call initial(no, jtrun, jtmax, lev, nx, my, my_max, mlmax)
   end if
   ct = 0.
   gt = 0.
   !                            SL_sedi, sat_pre, new_sat, use_cpm, use_dec
   call gsfcgce_3ice_nuwrf_unit(.false., .false., .false., .false., .false., .false., ct, gt)
   call gsfcgce_3ice_nuwrf_unit( .true., .false., .false., .false., .false., .false., ct, gt)
   call gsfcgce_3ice_nuwrf_unit( .true., .false., .false., .false., .false.,  .true., ct, gt)
   call gsfcgce_3ice_nuwrf_unit(.false.,  .true., .false., .false., .false.,  .true., ct, gt)
   call gsfcgce_3ice_nuwrf_unit(.false., .false.,  .true., .false., .false.,  .true., ct, gt)
   call gsfcgce_3ice_nuwrf_unit(.false., .false., .false.,  .true., .false.,  .true., ct, gt)
   call gsfcgce_3ice_nuwrf_unit(.false.,  .true.,  .true., .false., .false.,  .true., ct, gt)
   call gsfcgce_3ice_nuwrf_unit(.false.,  .true., .false.,  .true., .false.,  .true., ct, gt)
   call gsfcgce_3ice_nuwrf_unit(.false.,  .true., .false., .false.,  .true.,  .true., ct, gt)
   
   write (*, *) 'Timing for CPU & GPU: ', ct, gt
   write (*, *) 'speedup ratio: ', ct/gt, myrank
   ct = 0.
   gt = 0.
   !                   SL_sedi, sat_pre, new_sat, use_cpm, use_dec
   call saticel_s_unit(.false., .false., .false., .false., .false., .false., ct, gt)
   call saticel_s_unit(.false.,  .true., .false., .false., .false.,  .true., ct, gt)
   call saticel_s_unit(.false., .false.,  .true., .false., .false.,  .true., ct, gt)
   call saticel_s_unit(.false., .false., .false.,  .true., .false.,  .true., ct, gt)
   call saticel_s_unit(.false.,  .true.,  .true., .false., .false.,  .true., ct, gt)
   call saticel_s_unit(.false.,  .true., .false.,  .true., .false.,  .true., ct, gt)
   call saticel_s_unit(.false.,  .true., .false., .false.,  .true.,  .true., ct, gt)
   write (*, *) 'Timing for CPU & GPU: ', ct, gt
   write (*, *) 'speedup ratio: ', ct/gt, myrank
   
   ct = 0.
   gt = 0.
   call fall_flux_unit(.false., .false., .false., .false., .false., .false., ct, gt)
   call fall_flux_unit( .true., .false., .false., .false., .false., .false., ct, gt)
   call fall_flux_unit( .true., .false., .false., .false., .false.,  .true., ct, gt)
   call fall_flux_unit(.false.,  .true., .false., .false., .false.,  .true., ct, gt)
   call fall_flux_unit(.false., .false.,  .true., .false., .false.,  .true., ct, gt)
   write (*, *) 'Timing for CPU & GPU: ', ct, gt
   write (*, *) 'speedup ratio: ', ct/gt, myrank
   

   
   call mpe_finalize

end program

