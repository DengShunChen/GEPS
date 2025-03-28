!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (c) 2024,
!
! See LICENSE for license information.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mod_ut_initial
   use const, only: RTYPE
   use mod_vcopy
   implicit none

   type diag_var
      real(kind=RTYPE), dimension(:, :, :, :), allocatable:: vornow, divnow, temnow
      real(kind=RTYPE), dimension(:, :, :), allocatable:: plnow
      real(kind=RTYPE), dimension(:, :, :), allocatable:: rvor, rdiv, tt, ut, vt
      real(kind=RTYPE), dimension(:, :), allocatable :: pt, dlpl, dtpl
   end type diag_var

contains
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
      ! var_name(1) = 'plt'
      ! call stats_var(plt, var_name, lev, 1)
      var_name(1) = 'pk'
      call stats_var(pk, var_name, lev, 1)
      var_name(1) = 'pk2'
      call stats_var(pk2, var_name, lev, 1)
      var_name(1) = 'phi'
      call stats_var(phi, var_name, lev, 1)
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
                                         sum_horizon(lev*num), varmin
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
      varmin = varamax
      do jj = 1, jlistnum
         j = jlist1(jj)
         nxj = nxdef_2d(j)
         do n = 1, num
            do k = 1, lev
               kk = k + (n - 1)*lev
               do i = 1, nxj
                  varmin(n) = min(varmin(n), var(i, k, n, jj))
               end do
            end do
         end do
      end do
      call mpe_global_min(varmin, num)


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
            write (*, '(1X, A6, 1pe15.7, 1pe15.7, 2(", ", 1pe15.7))') &
                 name(n), varamax(n), varmin(n), &
                 varsum_t1(n)/npts, varsum_t2(n)/npts
         end do
      end if

   end subroutine stats_var
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
   subroutine get_initial_fields(DOREAD)
      ! ** output **
      ! tt, ut, vt, qt, pt, plt, pk, pk2, dlpl, dtpl, phi,
      ! pdiff, tsave, t1000,
      ! o3l,
      ! temnow, plnow, vornow, divnow
      !
      ! ** input **
      ! sgeo
      ! doincr, idtg, ggdef, gmdef,
      ! cstar, ktrop, ptop, capa, grav, rgas, rad, cp, taup,
      ! poly, dpoly, weight, cim, onocos, wcfac, wdfac, cosl,
      ! sigma, dsigma
      !
      use rank, only: myrank
      use const, only: RTYPE, doincr, idtg, ggdef, gmdef, &
                       cstar, ktrop, ptop, capa, grav, rgas, rad, cp, taup, &
                       poly, dpoly, weight, cim, onocos, wcfac, wdfac, cosl, &
                       sigma, dsigma
      use param, only: nx, my, my_max, lev, ncld, jtrun, jtmax
      use index, only: levp, jlistnum, jlist1, nxp, nxjp, nxdef_2d, &
                       nsizey, mlist, mlistnum, nsizex, nsizex
      use grid, only: tt, ut, vt, qt, pt, rdiv, rvor, &
                      plt, pk, pk2, dlpl, dtpl, &
                      phi, pdiff, tsave, t1000, sgeo
      use spec, only: temnow, vornow, divnow, plnow, &
                      temold, vorold, divold, plold
      use phygrid, only: o3l

      implicit none
      logical DOREAD
      integer lmax
      real taux
      real(kind=RTYPE) cc(nx + 2, levp, my_max), dummy, ww3(nx, my_max)
      integer j, jj, nxj, m, mf, k
      character(len=256) file
      integer fid
      character(len=9) myrank_s
      lmax = 26

      if (DOREAD) then
         if (myrank .eq. 0) print *, "Obatin the initail fields by reading file"
         fid = 200 + myrank
         write (myrank_s, '("x", i1, "y", i0.2, "Id", i0.2)') &
            nsizex, nsizey, myrank
         file = "UT_inp_"
         open (fid, file=trim(file)//trim(myrank_s)//".bin", &
               status='old', ACCESS='stream', form='unformatted')
         read (fid) sgeo, &
            tt, ut, vt, qt, pt, rdiv, &
            plt, pk, pk2, dlpl, dtpl, phi, &
            pdiff, tsave, t1000, o3l, &
            temnow, vornow, divnow, plnow
         ! temold, vorold, divold, plold
         close (fid)
         ! call prt_varinfo
      else
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
         ! << vornow and divnow >>
         ! compute vorticity and divergence from u and v
         call trandv(jtrun, jtmax, nx, my, my_max, lev, ut, vt, weight, cim, &
                     onocos, poly, dpoly, vornow, divnow, nsizey)
         ! << plnow >>
         call mpe2d_unify_nx(ww3, pt)
         call tranrs1(jtrun, jtmax, nx, my, my_max, poly, weight, ww3, &
                      plnow, nsizey)

         do m = 1, mlistnum
            mf = mlist(m)
            if (mf .eq. 1) then
               do k = 1, levp
                  divnow(k, 1, 1, m) = 0.0
                  divnow(k, 2, 1, m) = 0.0
                  vornow(k, 1, 1, m) = 0.0
                  vornow(k, 2, 1, m) = 0.0
               end do
            end if
         end do

         if (.not. doincr) then
            call vcopy_spec(divold, divnow)
            call vcopy_spec(vorold, vornow)
            call vcopy_spec(temold, temnow)
            call vcopy_spec(plold, plnow)
         end if  !end of (not doincr)

         ! << dlpl, dtpl >>
         ! zonal and meridional gradients of terrain pressure
         call trngra(jtrun, jtmax, nx, my, my_max, cim, poly, dpoly, plnow, &
                     dlpl, dtpl, nsizey)
         ! ! << rvor >>
         ! call transr(jtrun, jtmax, nx, my, my_max, levp, poly, vornow, cc, 1, nsizey)
         ! call ujoinsr(cc, rvor, dummy, dummy, dummy, nx, my_max, lev, jlistnum, 1, 1)
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
         fid = 200 + myrank
         write (myrank_s, '("x", i1, "y", i0.2, "Id", i0.2)') &
            nsizex, nsizey, myrank
         file = "UT_inp_"
         open (fid, file=trim(file)//trim(myrank_s)//".bin", &
               status='replace', ACCESS='stream', form='unformatted')
         write (fid) sgeo, &
            tt, ut, vt, qt, pt, rdiv, &
            plt, pk, pk2, dlpl, dtpl, phi, &
            pdiff, tsave, t1000, o3l, &
            temnow, vornow, divnow, plnow
         ! temold, vorold, divold, plold
         close (fid)
      end if
   end subroutine get_initial_fields

   subroutine initial_unit(no)
      use rank, only: myrank
      use param, only: nx, my, lev, jtrun, ncld, mlmax, &
                       my_max, jtmax
      use const, only: RTYPE, nnmiit, &
                       ! ----------
                       poly, dpoly, polyf, dpolyf, weight, cim, &
                       wcfac, wdfac, onocos  ! tranuv, trandv
      use index, only: mlistnum, mlist, nlist, &
                       jlistnum, jlist1, jlist2, jlist2_2d, jlist1_sl, &
                       nxdef, mtrundef, &
                       nsizex, nsizey, tcolt_jlist, poly_mlist, &
                       nxp, levp, nxptot, &
                       Llist, nxjp, nxjstart, nxjend, nxjlen, &
                       nxdef_2d, nxjp_acc, nxjlen_all

      use spec, only: vornow, divnow, temnow, plnow, &
                      ! ----------
                      vorold, divold, temold, plold, &
                      vorten, divten, temten, &
                      hldten, plten
      use grid, only: rvor, rdiv, tt, ut, vt, pt, dlpl, dtpl, &
                      ! ----------
                      qt, &
                      pk, pk2, plt, &
                      sgeo, phi, dlphi, dtphi, &
                      sd, vvel
      implicit none
      type(diag_var):: cpu

      integer no
      integer steps

      real(kind=8) tm_1, tm_2, tm_use, mpi_wtime

      logical DOCPU

      integer i, k, s
      integer async_id

      DOCPU = .true.
      async_id = 1
      steps = 1
      nnmiit = 3
      allocate (cpu%vornow(levp, 2, jtrun, jtmax), &
                cpu%divnow(levp, 2, jtrun, jtmax), &
                cpu%temnow(levp, 2, jtrun, jtmax), &
                cpu%plnow(jtrun, jtmax, 2), &
                cpu%rvor(nxp, lev, my_max), &
                cpu%rdiv(nxp, lev, my_max), &
                cpu%tt(nxp, lev, my_max), &
                cpu%ut(nxp, lev, my_max), &
                cpu%vt(nxp, lev, my_max), &
                cpu%pt(nxp, my_max), &
                cpu%dlpl(nxp, my_max), &
                cpu%dtpl(nxp, my_max) &
                )
      ! << CPU >>
      if (DOCPU) then
         if (myrank .eq. 0) print *, "Run the CPU version."
         call get_initial_fields(doread=.false.)
         ! call get_nnmi_const(doread=.false.)
         tm_1 = mpi_wtime()
         call initial(no, jtrun, jtmax, lev, nx, my, my_max, mlmax)
         tm_use = mpi_wtime() - tm_1
         if (myrank .eq. 0) write (*, '(A, f15.3, A)') &
            "Elapsed time of CPU: ", tm_use, " s"
         call vcopy_spec(cpu%vornow, vornow)
         call vcopy_spec(cpu%divnow, divnow)
         call vcopy_spec(cpu%temnow, temnow)
         call vcopy_spec(cpu%plnow, plnow)
         call vcopy_grid(cpu%rvor, rvor, 1)
         call vcopy_grid(cpu%rdiv, rdiv, 1)
         call vcopy_grid(cpu%tt, tt, 1)
         call vcopy_grid(cpu%ut, ut, 1)
         call vcopy_grid(cpu%vt, vt, 1)
         call vcopy_grid(cpu%pt, pt)
         call vcopy_grid(cpu%dlpl, dlpl)
         call vcopy_grid(cpu%dtpl, dtpl)
         call write_file(cpu)
      else
         if (myrank .eq. 0) print *, "Not running the CUP version. "// &
            "Obtain the cpu version results by reading file."
         call read_file(cpu)
      end if

      ! << GPU >>
      if (myrank .eq. 0) print *, "Run the GPU version. "
      !! << param >>
      !$acc enter data copyin(nx, my, lev, jtrun, ncld, &
      !$acc&                  my_max, jtmax) async(async_id)
      !! << index >>
      !$acc enter data copyin(mlistnum, mlist, nlist, &
      !$acc& jlistnum, jlist1, jlist2, jlist2_2d, jlist1_sl, &
      !$acc& nxdef, mtrundef, &
      !$acc& nsizex, nsizey, tcolt_jlist, poly_mlist, &
      !$acc& nxp, levp, nxptot, &
      !$acc& Llist, nxjp, nxjstart, nxjend, nxjlen, nxdef_2d, nxjp_acc, nxjlen_all &
      !$acc& ) async(async_id)
      !! << const >>
      !$acc enter data copyin(poly, dpoly, polyf, dpolyf, weight, cim, &
      !$acc&                  wcfac, wdfac, onocos) async(async_id)

      !! << spec >>
      !$acc enter data create(temten, vorten, divten, hldten, plten, &
      !$acc&                  temold, vorold, divold, plold, &
      !$acc&                  temnow, vornow, divnow, plnow &
      !$acc& ) async(async_id)
      !! << grid >>
      !$acc enter data create(rdiv, ut, vt, tt, qt, pt, &
      !$acc&                  dtpl, dlpl, &
      !$acc&                  pk, pk2, plt, sgeo, phi, dtphi, dlphi, &
      !$acc&                  sd, vvel &
      !$acc& )async(async_id)

      !$acc update device(sgeo) async(async_id)
      do s = 1, steps
         call get_initial_fields(doread=.true.)
         !$acc update device(pt, ut, vt, tt, qt) async(async_id)
         !$acc update device(rdiv, phi, dtpl, dlpl) async(async_id)
         !$acc update device(vornow, divnow, temnow, plnow) async(async_id)
         !$acc wait(async_id)
         tm_1 = mpi_wtime()
         call initial_gpu(no, jtrun, jtmax, lev, nx, my, my_max, mlmax)
         tm_use = mpi_wtime() - tm_1
         if (myrank .eq. 0) write (*, '(A, f15.3, A)') &
            "Elapsed time of GPU: ", tm_use, " s"
      end do
      !$acc update self(pt, ut, vt, tt, qt) async(async_id)
      !$acc update self(rdiv, dtpl, dlpl) async(async_id)
      !$acc update self(vornow, divnow, temnow, plnow) async(async_id)
      !$acc update self(vorold, divold, temold, plold) async(async_id)
      !$acc update self(vorten, divten, plten) async(async_id)
      !$acc update self(pk, pk2, plt, sd, vvel) async(async_id)
      !$acc update self(hldten, dlphi, dtphi) async(async_id)

      !! << param >>
      !$acc exit data delete(nx, my, lev, jtrun, ncld, &
      !$acc&                  my_max, jtmax) async(async_id)
      !! << index >>
      !$acc exit data delete(mlistnum, mlist, nlist, &
      !$acc& jlistnum, jlist1, jlist2, jlist2_2d, jlist1_sl, &
      !$acc& nxdef, mtrundef, &
      !$acc& nsizex, nsizey, tcolt_jlist, poly_mlist, &
      !$acc& nxp, levp, nxptot, &
      !$acc& Llist, nxjp, nxjstart, nxjend, nxjlen, nxdef_2d, nxjp_acc, nxjlen_all &
      !$acc& ) async(async_id)
      !! << const >>
      !$acc exit data delete(poly, dpoly, polyf, dpolyf, weight, cim, &
      !$acc&                  wcfac, wdfac, onocos) async(async_id)

      !! << spec >>
      !$acc exit data delete(temten, vorten, divten, hldten, plten, &
      !$acc&                  temold, vorold, divold, plold, &
      !$acc&                  temnow, vornow, divnow, plnow &
      !$acc& ) async(async_id)
      !! << grid >>
      !$acc exit data delete(rdiv, ut, vt, tt, qt, pt, &
      !$acc&                  dtpl, dlpl, &
      !$acc&                  pk, pk2, plt, sgeo, phi, dtphi, dlphi, &
      !$acc&                  sd, vvel &
      !$acc& )async(async_id)
      !$acc wait(async_id)
      call assert_output(cpu)

   end subroutine initial_unit
! ============================================================
   subroutine assert_output(cpu)
      use rank, only: myrank
      use const, only: RTYPE
      use index, only: nxp
      use param, only: lev
      use spec, only: vornow, divnow, temnow, plnow
      use grid, only: rvor, rdiv, tt, ut, vt, pt, dlpl, dtpl
      implicit none
      type(diag_var):: cpu
      real(kind=RTYPE) err(3)
      logical equal
      logical assert_specspace_close, assert_specspace_2d_close
      character*6:: lab

      equal = .true.

      if (myrank .eq. 0) then
         write (*, '(1X, A8, 1X, A15, 4A15)') &
            "name", "max|Err|/max|A|", "L2-Err", "RMSE", "max|A|", "max|B|"
      end if
      call varerr_spec3D(err, vornow, cpu%vornow, 'vornow')
      call varerr_spec3D(err, divnow, cpu%divnow, 'divnow')
      call varerr_spec3D(err, temnow, cpu%temnow, 'temnow')
      ! call varerr(err, rvor, cpu%rvor, lev, 1, 'rvor')
      call varerr(err, rdiv, cpu%rdiv, lev, 1, 'rdiv')
      call varerr(err, tt, cpu%tt, lev, 1, 'tt')
      call varerr(err, ut, cpu%ut, lev, 1, 'ut')
      call varerr(err, vt, cpu%vt, lev, 1, 'vt')
      call varerr(err, pt, cpu%pt, 1, 1, 'pt')
      call varerr(err, dlpl, cpu%dlpl, 1, 1, 'dlpl')
      call varerr(err, dtpl, cpu%dtpl, 1, 1, 'dtpl')

#ifdef SP
      equal = assert_specspace_close(vornow, cpu%vornow, 1e-4_4, 1e-4_4, "vornow") &
              .and. equal
      equal = assert_specspace_close(divnow, cpu%divnow, 1e-4_4, 1e-4_4, "divnow") &
              .and. equal
      ! equal = assert_specspace_close(temnow, cpu%temnow, 1e-4_4, 1e-4_4, "temnow") &
      !         .and. equal
      equal = assert_specspace_2d_close(plnow, cpu%plnow, 1e-4_4, 1e-4_4, "plnow") &
              .and. equal

      call assert_realspace_close(rdiv, cpu%rdiv, nxp, lev, 1, &
                                  1.e-4_4, 1.e-4_4, "rdiv")
      call assert_realspace_close(tt, cpu%tt, nxp, lev, 1, &
                                  1.e-4_4, 1.e-4_4, "tt")
      call assert_realspace_close(ut, cpu%ut, nxp, lev, 1, &
                                  1.e-4_4, 1.e-4_4, "ut")
      call assert_realspace_close(vt, cpu%vt, nxp, lev, 1, &
                                  1.e-4_4, 1.e-4_4, "vt")
      call assert_realspace_close(pt, cpu%pt, nxp, 1, 1, &
                                  1.e-4_4, 1.e-4_4, "pt")
#else
      equal = assert_specspace_close(vornow, cpu%vornow, 1e-10, 1e-10, "vornow") &
              .and. equal
      equal = assert_specspace_close(divnow, cpu%divnow, 1e-10, 1e-10, "divnow") &
              .and. equal
      equal = assert_specspace_close(temnow, cpu%temnow, 1e-10, 1e-10, "temnow") &
              .and. equal
      equal = assert_specspace_2d_close(plnow, cpu%plnow, 1e-10, 1e-10, "plnow") &
              .and. equal

      ! call assert_realspace_close(rvor, cpu%rvor, nxp, lev, 1, &
      !                             1.e-10, 1.e-10, "rvor")
      call assert_realspace_close(rdiv, cpu%rdiv, nxp, lev, 1, &
                                  1.e-10, 1.e-10, "rdiv")
      call assert_realspace_close(tt, cpu%tt, nxp, lev, 1, &
                                  1.e-10, 1.e-10, "tt")
      call assert_realspace_close(ut, cpu%ut, nxp, lev, 1, &
                                  1.e-10, 1.e-10, "ut")
      call assert_realspace_close(vt, cpu%vt, nxp, lev, 1, &
                                  1.e-10, 1.e-10, "vt")
      call assert_realspace_close(pt, cpu%pt, nxp, 1, 1, &
                                  1.e-10, 1.e-10, "pt")

      ! call assert_realspace_close(dlpl, cpu%dlpl, nxp, 1, 1, &
      !                             1.e-10, 1.e-10, "dlpl")
      ! call assert_realspace_close(dtpl, cpu%dtpl, nxp, 1, 1, &
      !                             1.e-10, 1.e-10, "dtpl")
#endif
      if (.not. equal) call exit(1)
   end subroutine assert_output
! ============================================================
   subroutine VarErr_spec3D(Err, a_spec, b_spec, lab)
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
    end subroutine VarErr_Spec3D
! ============================================================
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
      real(kind=8) vamax
      real(kind=RTYPE) Err_loc(nvar, 3), avar_amax(nvar, 2), sum_local(nvar), pi, tmp

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
                  avar_amax(n, 1) = max(avar_amax(n, 1), abs(A(i, k, n, jj)))
                  avar_amax(n, 2) = max(avar_amax(n, 2), abs(B(i, k, n, jj)))

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
      call mpe_global_max(avar_amax(1, 1), nvar*2, RTYPE)

      do n = 1, nvar
         Err_loc(n, 2) = sqrt(Err_loc(n, 2))
         Err_loc(n, 3) = sqrt(Err_loc(n, 3)/pts/nvar)
      end do

      if (myrank .eq. 0) then
         do n = 1, nvar
            Err_loc(n, 1) = Err_loc(n, 1)/avar_amax(n, 1)
            write (*, '(1X, A6, i0.2, 1X, 1pe15.7, 4(1pe15.7))') &
               trim(lab), n, Err_loc(n, 1:3), avar_amax(n, 1:2)
         end do
      end if

      return
   end subroutine VarErr

   subroutine write_file(dat)
      use rank, only: myrank
      use index, only: nsizex, nsizey
      implicit none
      type(diag_var):: dat
      character(len=256) file
      integer fid
      character(len=9) myrank_s

      fid = 200 + myrank
      write (myrank_s, '("x", i1, "y", i0.2, "Id", i0.2)') &
         nsizex, nsizey, myrank
      file = "UT_out_"
      open (fid, file=trim(file)//trim(myrank_s)//".bin", &
            status='replace', ACCESS='stream', form='unformatted')
      write (fid) dat%vornow, dat%divnow, dat%temnow, dat%plnow
      write (fid) dat%rvor, dat%rdiv, dat%tt, dat%ut, dat%vt
      write (fid) dat%pt, dat%dlpl, dat%dtpl
      close (fid)

   end subroutine write_file

   subroutine read_file(dat)
      use rank, only: myrank
      use index, only: nsizex, nsizey
      implicit none
      type(diag_var):: dat
      character(len=256) file
      integer fid
      character(len=9) myrank_s

      fid = 200 + myrank
      write (myrank_s, '("x", i1, "y", i0.2, "Id", i0.2)') &
         nsizex, nsizey, myrank
      file = "UT_out_"
      open (fid, file=trim(file)//trim(myrank_s)//".bin", &
            status='old', ACCESS='stream', form='unformatted')
      read (fid) dat%vornow, dat%divnow, dat%temnow, dat%plnow
      read (fid) dat%rvor, dat%rdiv, dat%tt, dat%ut, dat%vt
      read (fid) dat%pt, dat%dlpl, dat%dtpl
      close (fid)

   end subroutine read_file
end module mod_ut_initial

program test_initail
   use param
   use const
   use mod_ut_initial

   implicit none
   integer no
   
   call mpe_init
   call cons
   no = 2*((jtrun + 1)/2) + (jtrun/2) + 10
   call initial_unit(no)
   call mpe_finalize

end program test_initail
