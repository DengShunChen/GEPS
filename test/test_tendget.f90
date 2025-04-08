!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024,
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_tendget
   use param
   use const

   implicit none
   integer no
   no = 2*((jtrun + 1)/2) + (jtrun/2) + 10

   call mpe_init
   call cons
   call get_initial_fields
   call tendget_unit(no)
   call mpe_finalize

end program test_tendget

! use param, only: nco, nx, my, lev, jtrun, mlmax, octahedral, ncld, &
!                  jtmax, my_max, npex, mpey
! ------------------------------------------------------------
subroutine read_sgeo_dms
   ! ** output **
   ! sgeo, spgeo
   !
   ! ** input **
   ! ihdgi, ggdef, ksgeo, bckfile
   ! grav, poly, dpoly, weight
   ! jlistnum, jlist1, nxjstart, nxdef, nxdef_2d, nsizey, lreduce
   ! nx, my, my_max, jtrun, jtmax, RTYPE
   !
   use const, only: RTYPE, ihdgi, ggdef, ksgeo, bckfile, &
                    grav, poly, dpoly, weight
   use param, only: nx, my, my_max, jtrun, jtmax
   use index, only: jlistnum, jlist1, nxjstart, nxdef, nxdef_2d, nsizey, &
                    lreduce
   use grid, only: sgeo
   use spec, only: spgeo
   implicit none
   real ww1(nx, my)
   real(kind=RTYPE) ww3(nx, my_max)
   character topohgt*4
   integer nxmy
   integer i, j, ii, jj, nxj
   integer istat

   nxmy = nx*my
   ! read terrain geopotential(sgeo) from data base
   if (ksgeo .eq. 99) then
      topohgt = 'gbkf'
   else
      write (topohgt, '(a3,i1.1)') 'gbk', ksgeo
   end if
#ifdef I38K
   write (ihdgi, '("s00060",2x,a4,a4,12x)') topohgt, ggdef
#else
   write (ihdgi, '("s00060",a4,a4,12x)') topohgt, ggdef
#endif
   call dmsread(nx, my, nxmy, 'H', bckfile, ww1, istat)
   if (istat .ne. 0) then
      call mpe_finalize
      call dmsexit(-1)
   end if
   do jj = 1, jlistnum
      j = jlist1(jj)
      ii = nxjstart(j)
      nxj = nxdef_2d(j)
      if (lreduce .eq. 1) call reducepick(ww1(1, j), nxdef(j), nx, 1)
      do i = 1, nxj
         sgeo(i, jj) = ww1(ii, j)*grav
         ii = ii + 1
      end do
   end do
   ! << spectrally truncate >>
   call mpe2d_unify_nx(ww3, sgeo)
   call tranrs1(jtrun, jtmax, nx, my, my_max, poly, weight, ww3, spgeo, nsizey)
   call transr1(jtrun, jtmax, nx, my, my_max, poly, spgeo, sgeo, nsizey)
end subroutine read_sgeo_dms

! ------------------------------------------------------------
subroutine get_initial_fields
   ! ** output **
   ! tt, ut, vt, qt, pt, plt, pk, pk2, dlpl, dtpl, phi,
   ! pdiff, tsave, t1000,
   ! o3l,
   ! temnow, plnow, vornow, divnow
   !
   ! ** input **
   ! sgeo
   ! idtg, ggdef, gmdef,
   ! cstar, ktrop, ptop, capa, grav, rgas, rad, cp, taup,
   ! poly, dpoly, weight, cim, onocos, wcfac, wdfac, cosl,
   ! sigma, dsigma
   !
   use rank, only: myrank
   use const, only: RTYPE, idtg, ggdef, gmdef, &
                    cstar, ktrop, ptop, capa, grav, rgas, rad, cp, taup, &
                    poly, dpoly, weight, cim, onocos, wcfac, wdfac, cosl, &
                    sigma, dsigma
   use param, only: nx, my, my_max, lev, ncld, jtrun, jtmax
   use index, only: levp, jlistnum, jlist1, nxp, nxjp, nxdef_2d, &
                    nsizey
   use grid, only: tt, ut, vt, qt, pt, rdiv, rvor, &
                   plt, pk, pk2, dlpl, dtpl, &
                   phi, pdiff, tsave, t1000, sgeo
   use spec, only: temnow, plnow, vornow, divnow
   use phygrid, only: o3l

   implicit none
   integer lmax
   real taux
   real(kind=RTYPE) cc(nx + 2, levp, 1, my_max), dummy, ww3(nx, my_max)
   integer j, jj, nxj
   lmax = 26

   if (myrank .eq. 0) print *, "Read terrain geopotential(sgeo) from DMS"
   call read_sgeo_dms
   if (myrank .eq. 0) print *, 'idtg=', idtg
   taux = 0.
   !! input: sgeo
   !! output: phi, tt, ut, vt, qt, pt, pdiff, plt, pk, pk2
   call sigful(nx, my, my_max, lev, ncld, lmax, jtrun, jtmax, &
               cstar, ktrop, idtg, ptop, taux, capa, grav, rgas, rad, &
               cp, weight, poly, sigma, cosl, phi, tt, ut, vt, qt, o3l, pt, sgeo, &
               pdiff, tsave, t1000, plt, pk, pk2, taup, &
               ggdef, gmdef)
   ! << temnow >>
   call joinrs(cc, tt, dummy, dummy, dummy, nx, my_max, lev, jlistnum, 1, 1)
   call tranrs(jtrun, jtmax, nx, my, my_max, levp, poly, weight, cc, &
               temnow, 1, nsizey)
   ! << plnow >>
   call mpe2d_unify_nx(ww3, pt)
   call tranrs1(jtrun, jtmax, nx, my, my_max, poly, weight, ww3, &
                plnow, nsizey)
   ! << vornow and divnow >>
   ! compute vorticity and divergence from u and v
   call trandv(jtrun, jtmax, nx, my, my_max, lev, ut, vt, weight, cim, &
               onocos, poly, dpoly, vornow, divnow, nsizey)
   ! << dlpl, dtpl >>
   ! zonal and meridional gradients of terrain pressure
   call trngra(jtrun, jtmax, nx, my, my_max, cim, poly, dpoly, plnow, &
               dlpl, dtpl, nsizey)
   ! << rvor >>
   call transr(jtrun, jtmax, nx, my, my_max, levp, poly, vornow, cc, 1, nsizey)
   call ujoinsr(cc, rvor, dummy, dummy, dummy, nx, my_max, lev, jlistnum, 1, 1)
   ! << rdiv >>
   call transr(jtrun, jtmax, nx, my, my_max, levp, poly, divnow, cc, 1, nsizey)
   call ujoinsr(cc, rdiv, dummy, dummy, dummy, nx, my_max, lev, jlistnum, 1, 1)
   ! << (spectrally truncated) ut, vt >>
   ! compute the spectrally truncated velocity coefficients
   call tranuv(jtrun, jtmax, nx, my, my_max, levp, onocos, wcfac, wdfac, &
               poly, dpoly, vornow, divnow, ut, vt, nsizey)
   ! << (spectrally truncated) tt >>
   call transr(jtrun, jtmax, nx, my, my_max, levp, poly, temnow, cc, 1, nsizey)
   call ujoinsr(cc, tt, dummy, dummy, dummy, nx, my_max, lev, jlistnum, 1, 1)
   ! << (spectrally trucated) pt >>
   call transr1(jtrun, jtmax, nx, my, my_max, poly, plnow, pt, nsizey)
   ! << pk, pk2, plt >>
   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef_2d(j)
      call prexp_hybrid_cwb(nxjp(j), nxp, lev, ptop, sigma, pt(1, jj), &
                            pk(1, 1, jj), pk2(1, 1, jj), plt(1, 1, jj))
   end do

   if (myrank .eq. 0) print *, "Completed the generation of initail fields."
end subroutine get_initial_fields

subroutine prt_varinfo
   use param, only: lev, ncld
   use grid, only: tt, ut, vt, qt, pt, rdiv, rvor, &
                   plt, pk, pk2, dlpl, dtpl, &
                   phi, t1000, sgeo
   implicit none

   character*6:: var_name(ncld)
   integer j

   var_name(1) = 'sgeo'
   call stats_var(sgeo, var_name, 1, 1)
   var_name(1) = 'pt'
   call stats_var(pt, var_name, 1, 1)
   var_name(1) = 'dlpl'
   call stats_var(dlpl, var_name, 1, 1)
   var_name(1) = 'dtpl'
   call stats_var(dtpl, var_name, 1, 1)

   var_name(1) = 'rvor'
   call stats_var(rvor, var_name, lev, 1)
   var_name(1) = 'rdiv'
   call stats_var(rdiv, var_name, lev, 1)
   var_name(1) = 'tt'
   call stats_var(tt, var_name, lev, 1)
   var_name(1) = 'ut'
   call stats_var(ut, var_name, lev, 1)
   var_name(1) = 'vt'
   call stats_var(vt, var_name, lev, 1)
   do j = 1, ncld
      write (var_name(j), '(A1, i1)') 'q', j
   end do
   call stats_var(qt, var_name, lev, ncld)
   var_name(1) = 'plt'
   call stats_var(plt, var_name, lev, 1)
   var_name(1) = 'pk'
   call stats_var(pk, var_name, lev, 1)
   var_name(1) = 'pk2'
   call stats_var(pk2, var_name, lev, 1)
   var_name(1) = 'phi'
   call stats_var(phi, var_name, lev, 1)
   var_name(1) = 't1000'
   call stats_var(t1000, var_name, lev, 1)
end subroutine prt_varinfo

subroutine stats_var(var, name, lev, num)
   use rank, only: myrank
   use const, only: RTYPE, numreduce
   use param, only: nx, my, my_max, octahedral
   use index, only: nxp, jlistnum, jlist1, nxdef_2d
   use mpe
   implicit none
   real(kind=RTYPE), intent(in):: var(nxp, lev, num, my_max)
   character*6, intent(in)::  name(num)
   integer, intent(in):: num, lev

   real(kind=RTYPE) wrk(nxp, lev, num, my_max)
   real(kind=RTYPE), dimension(num):: varamax, varsum_t1, varsum_t2, &
                                      sum_horizon(lev*num)
   integer npts
   integer i, j, jj, k, n, nxj, kk

   if (octahedral) then
      npts = (20 + nx)*my*lev
   elseif (numreduce == -99) then
      npts = nx*my*lev
   else
      npts = 0
   end if
   varamax = 0.
   varsum_t1 = 0.
   varsum_t2 = 0.
   sum_horizon = 0.
   ! varmax  = 0.
   ! varmin  =
   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef_2d(j)
      do n = 1, num
         do k = 1, lev
            kk = k + (n - 1)*lev
            do i = 1, nxj
               varamax(n) = max(varamax(n), abs(var(i, k, n, jj)))
               varsum_t1(n) = varsum_t1(n) + var(i, k, n, jj)
               sum_horizon(kk) = sum_horizon(kk) + var(i, k, n, jj)
            end do
         end do
      end do
   end do
   call mpe_global_max(varamax, num)
   call mpe_global_sum(varsum_t1, num)
   call mpe_global_sum(sum_horizon, lev*num)
   do n = 1, num
      do k = 1, lev
         kk = k + (n - 1)*lev
         varsum_t2(n) = varsum_t2(n) + sum_horizon(kk)
      end do
   end do

   if (myrank .eq. 0) then
      do n = 1, num
         write (*, '(1X, A6, 1pe15.7, 2(", ", 1pe15.7))') &
            name(n), varamax(n), varsum_t1(n)/npts, varsum_t2(n)/npts
      end do
   end if

end subroutine stats_var

subroutine tendget_unit(no)
   use rank, only: myrank
   use param, only: nx, my, my_max, lev, jtrun, jtmax, &
                    ! ----------------------------------------
                    ncld
   use const, only: RTYPE, mlmax, &
                    ! ----------------------------------------
                    poly, dpoly, polyf, dpolyf, cim, wcfac, wdfac, coslr, weight, &
                    radsq, onocos, cosl, spalm, arrhyd, &
                    ! ========================================
                    cp, sinl, dsigma, sigma, cor, ptop

   use index, only: levp, &
                    ! ----------------------------------------
                    jlistnum, jlist1, jlist2, mlistnum, mtrundef, mlist, nlist, Llist, &
                    nsizey, nsizex, &
                    nxjp, nxp, nxdef, nxdef_2d, nxjp_acc, nxptot, &
                    jlist2_2d, nxjstart, nxjend, jlist1_sl, &
                    tcolt_jlist, poly_mlist, nxjlen_all, nxjlen
   use grid, only: rdiv, pt, tt, ut, vt, qt, plt, pk, pk2, dtpl, dlpl, &
                   sd, vvel, up, vp, ttp, qm, &
                   ut_sl, vt_sl, &
                   phi, dlphi, dtphi, &
                   ! ----------------------------------------
                   sgeo, &
                   lonlen, lonstr, latlen, fa1, fa2, fa3, fa4, gglati

   use spec, only: temten, vorten, divten, hldten, plten

   use spec_cuda_graph, only: cc_cg, gwk1_cg, wcc_fk_cg, &
                              wc_cg, ws_cg, fj_weight_cg, &
                              allocate_spec_cg_buffer, &
                              deallocate_spec_cg_buffer
   use mod_vcopy
   use openacc
   implicit none
   integer, parameter:: nvar = 2
   integer no
   integer steps

   real, dimension(nxp, lev, my_max):: plt_c
   real(kind=RTYPE), dimension(nxp, lev, my_max):: pk_c, pk2_c, phi_c
   real(kind=RTYPE), dimension(nxp, my_max) :: dlpl_c, dtpl_c
   real(kind=RTYPE), dimension(levp, 2, jtrun, jtmax):: phiten_c, phiten_g, &
                                                        vorten_c, divten_c
   real(kind=RTYPE), dimension(lev*6*jtrun*jtmax):: buf
   real(kind=RTYPE) cc(nx + 2, levp, my_max)
   real(kind=RTYPE):: dummy
   real(kind=8) tm_1, tm_2, tm_use, mpi_wtime
   real(kind=RTYPE):: err_arr(3)
   character(len=6):: name

   type(acc_graph_t) matmul_cg
   logical matmul_cg_created, equal
   logical assert_specspace_close
   character*10:: lab
   integer mlmax2

   integer i, k
   integer async_id

   equal = .true.
   matmul_cg_created = .false.
   async_id = 1
   steps = 2
   mlmax2 = mlmax*2
   if (myrank .eq. 0) print *, 'jrtun=', jtrun, ' mlmax=', mlmax
   lab = 'pt'
   call check(pt, nx, my, my_max, lab)
   ! ========================================
   !$acc enter data copyin(nx, my, my_max, jlistnum, mlistnum, jtrun, jtmax, levp, nsizey, nsizex) async(async_id)
   !$acc enter data copyin(nxjp, nxp, lev, ncld) async(async_id)
   !$acc enter data copyin(jlist1, jlist2, nxdef_2d, nxjp_acc, nxptot, &
   !$acc& jlist2_2d, nxdef, mtrundef, mlist, nlist, Llist, &
   !$acc& nxjstart, nxjend) async(async_id)
   !! << const >>
   !$acc enter data copyin(&
   !$acc& poly, dpoly, polyf, dpolyf, cim, wcfac, wdfac, coslr, weight, &
   !$acc& radsq, onocos, cosl, spalm, arrhyd &
   !$acc& ) async(async_id)
   !$acc enter data copyin(cp, sinl, dsigma, sigma, cor, ptop) async(async_id)
   !! << index >>
   !$acc enter data copyin(tcolt_jlist, poly_mlist) async(async_id)
   !$acc enter data copyin(jlist1_sl, nxjlen_all, nxjlen) async(async_id)
   !! << grid >>
   !$acc enter data copyin(sgeo) async(async_id)
   !$acc enter data copyin(lonlen, lonstr, latlen, fa1, fa2, fa3, fa4, gglati) async(async_id)

   ! ========================================
   ! << CPU >>
   tm_1 = mpi_wtime()
   call tendget(phiten_c)
   tm_use = mpi_wtime() - tm_1
   if (myrank .eq. 0) write (*, '(A, f15.3, A)') &
      "Elapsed time of CPU: ", tm_use, " seconds"
   call vcopy_grid(plt_c, plt, 1)
   call vcopy_grid(pk_c, pk, 1)
   call vcopy_grid(pk2_c, pk2, 1)
   call vcopy_spec(vorten_c, vorten)
   call vcopy_spec(divten_c, divten)

   ! << GPU >>

   !$acc enter data create(rdiv, pt, tt, ut, vt, qt, &
   !$acc& plt, pk, pk2, dtpl, dlpl, &
   !$acc& sd, vvel, up, vp, ttp, qm, &
   !$acc& ut_sl, vt_sl, &
   !$acc& phi, dlphi, dtphi) async(async_id)
   !$acc enter data create(temten, vorten, divten, hldten, plten) async(async_id)

   !$acc enter data create(buf) async(async_id)
   !$acc enter data create(phiten_g) async(async_id)

   !$acc update device(pt, ut, vt, tt, qt) async(async_id)
   !$acc update device(rdiv, phi, dtpl, dlpl) async(async_id)
   !$acc wait(async_id)
   do i = 1, steps
      tm_1 = mpi_wtime()
      call tendget_gpu(phiten_g, buf(1), buf(1 + lev*2*jtrun*jtmax), &
                       matmul_cg_created, matmul_cg)
      !$acc wait(async_id)
      tm_use = mpi_wtime() - tm_1
      if (myrank .eq. 0) write (*, '(A, f15.3, A)') &
         "Elapsed time of GPU: ", tm_use, " seconds"
   end do
   !$acc update self(hldten, dlphi, dtphi, &
   !$acc& plten, vorten, divten, &
   !$acc& pk, pk2, plt, sd, vvel, &
   !$acc& phiten_g) async(async_id)

   !$acc exit data delete(rdiv, pt, tt, ut, vt, qt, &
   !$acc& plt, pk, pk2, dtpl, dlpl, &
   !$acc& sd, vvel, up, vp, ttp, qm, &
   !$acc& ut_sl, vt_sl, &
   !$acc& phi, dlphi, dtphi) async(async_id)
   !$acc exit data delete(temten, vorten, divten, hldten, plten) async(async_id)

   !$acc exit data delete(buf) async(async_id)
   !$acc exit data delete(phiten_g) async(async_id)
   !$acc wait(async_id)

#ifdef VERBOSE
   if (myrank .eq. 0) then
      write (*, '(1X, A6, A23, 3A15)') &
         "name", "max|Err|", "L2-Err", "RMSE", "max|Var|"
   end if
   call varerr(err_arr, pk, pk_c, lev, 1, 'pk')
   call varerr(err_arr, pk2, pk2_c, lev, 1, 'pk2')
   call varerr_spec3d(err_arr, vorten, vorten_c, 'vorten')
   call varerr_spec3d(err_arr, divten, divten_c, 'divten')
   call varerr_spec3d(err_arr, phiten_g, phiten_c, 'phiten')
#endif

#ifdef SP
   call assert_realspace_close(pk, pk_c, nxp, lev, 1, 1e-4_4, 1e-4_4, "Array pk")
   call assert_realspace_close(pk2, pk2_c, nxp, lev, 1, 1e-4_4, 1e-4_4, "Array pk2")

   equal = (equal .and. assert_specspace_close(vorten, vorten_c, 1e-4_4, 1e-4_4, &
                                               "vorten"))
   equal = (equal .and. assert_specspace_close(divten, divten_c, 1e-4_4, 1e-4_4, &
                                               "divten"))
   equal = (equal .and. assert_specspace_close(phiten_g, phiten_c, 1e-4_4, 1e-4_4, &
                                               "phiten"))
#else
   call varerr(err_arr, plt, plt_c, lev, 1, 'plt')

   call assert_realspace_close(plt, plt_c, nxp, lev, 1, 1e-10, 1e-10, "Array plt")
   call assert_realspace_close(pk, pk_c, nxp, lev, 1, 1e-10, 1e-10, "Array pk")
   call assert_realspace_close(pk2, pk2_c, nxp, lev, 1, 1e-10, 1e-10, "Array pk2")

   equal = (equal .and. assert_specspace_close(vorten, vorten_c, 1e-10, 1e-10, &
                                               "vorten"))
   equal = (equal .and. assert_specspace_close(divten, divten_c, 1e-10, 1e-10, &
                                               "divten"))
   equal = (equal .and. assert_specspace_close(phiten_g, phiten_c, 1e-10, 1e-10, &
                                               "phiten"))
#endif
   if (.not. equal) call exit(1)

end subroutine tendget_unit
! ============================================================
subroutine VarErr_spec3d(Err, a_spec, b_spec, lab)
   use const, only: RTYPE, poly
   use param, only: nx, my_max, my, lev, jtrun, jtmax
   use index, only: nxp, levp, nxdef, jlist1, jlistnum, &
                    nsizey
   use rank, only: myrank
   use mpe
   implicit none
   character(len=6), intent(in) :: lab
   real(kind=RTYPE), intent(out):: Err(3)
   real(kind=RTYPE), intent(in):: A_spec(levp, 2, jtrun, jtmax), &
                                  B_spec(levp, 2, jtrun, jtmax)

   real(kind=RTYPE), dimension(nxp, lev, my_max):: A_grid, B_grid
   real(kind=RTYPE):: cc(nx + 2, levp, my_max), dummy

   call transr(jtrun, jtmax, nx, my, my_max, levp, poly, A_spec, cc, 1, nsizey)
   call ujoinsr(cc, A_grid, dummy, dummy, dummy, nx, my_max, lev, jlistnum, 1, 1)

   call transr(jtrun, jtmax, nx, my, my_max, levp, poly, B_spec, cc, 1, nsizey)
   call ujoinsr(cc, B_grid, dummy, dummy, dummy, nx, my_max, lev, jlistnum, 1, 1)

   call VarErr(Err, a_grid, b_grid, lev, 1, lab)

   return
end subroutine VarErr_Spec3d

subroutine VarErr(Err, a, b, lev, nvar, lab)
   use const, only: RTYPE, numreduce, weight
   use param, only: nx, my_max, my, octahedral
   use index, only: nxp, nxdef, nxdef_2d, jlist1, jlistnum
   use rank, only: myrank
   use mpe
   implicit none
   ! Measure the difference of between two arrays.
   ! There are the three measurement methods.
   ! 1. Maximum norm
   ! 2. L2-norm (The vertical integral weights need to be corrected)
   ! 3. RMSE
   character(len=*), intent(in):: lab
   real(kind=RTYPE), intent(out):: Err(3, nvar)
   integer, intent(in)::lev, nvar
   real(kind=RTYPE), intent(in):: A(nxp, lev, nvar, my_max), &
                                  B(nxp, lev, nvar, my_max)

   integer i, j, k, n, nxj, jj, pts, lons
   real(kind=RTYPE) Err_loc(nvar, 3), avar_amax(nvar), sum_local(nvar), pi, tmp

   if (octahedral) then
      pts = (20 + nx)*my*lev
   elseif (numreduce == -99) then
      pts = nx*my*lev
   end if

   pi = 4.*atan(1.)
   Err = 0.
   Err_loc = 0.
   avar_amax = 0.

   do jj = 1, jlistnum
      j = jlist1(jj)
      nxj = nxdef_2d(j)
      lons = nxdef(j)
      sum_local = 0.
      do n = 1, nvar
         do k = 1, lev
            do i = 1, nxj
               avar_amax(n) = max(avar_amax(n), abs(A(i, k, n, jj)))

               tmp = A(i, k, n, jj) - B(i, k, n, jj)

               Err_loc(n, 1) = max(Err_loc(n, 1), abs(tmp))
               sum_local(n) = sum_local(n) + tmp**2
            end do
         end do
         Err_loc(n, 2) = Err_loc(n, 2) + sum_local(n)*(2.*pi/lons)*weight(j)/lev
         Err_loc(n, 3) = Err_loc(n, 3) + sum_local(n)
      end do
   end do

   call mpe_global_max(Err_loc(1, 1), nvar, RTYPE)
   call mpe_global_sum(Err_loc(1, 2), nvar*2, RTYPE)
   call mpe_global_max(avar_amax, nvar, RTYPE)

   do n = 1, nvar
      Err_loc(n, 2) = sqrt(Err_loc(n, 2))
      Err_loc(n, 3) = sqrt(Err_loc(n, 3)/pts/nvar)
   end do

   if (myrank .eq. 0) then
      do n = 1, nvar
         write (*, '(1X, A6, 1pe23.15, 3(1pe15.7))') &
            trim(lab), Err_loc(n, 1:3), avar_amax(n)
      end do
   end if

   return
end subroutine VarErr
