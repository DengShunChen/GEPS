! Tiled vertical advection. Full-londim work arrays (dd/ds/xgrid/xpast/xnext/da)
! are ~5GB plus PPM ~15GB; HIP 1-rank TCo383 has ~0.5GB free after NDSL.
! acc2omp splices this into module mod_ndslfv_monoadv_gpu.

subroutine vertical_cell_advect_gpu(lons, londim, levs, nvars, &
                                     deltim, ssi, wwi, qql, mass, forward, &
                                     async_id, ibase)
   use const, only: RTYPE
   use grid, only: latpart
   use index, only: jlistnum, jlist1, nxjp, nxp, nxjp_acc
   implicit none

   integer, intent(in) :: lons, londim, levs, nvars, mass, async_id
   integer, intent(in), optional :: ibase
   real(kind=RTYPE), intent(in) :: deltim
   real(kind=RTYPE), intent(in) :: ssi(levs + 1, londim)
   real(kind=RTYPE), intent(in) :: wwi(nxp, levs + 1, latpart)
   real(kind=RTYPE), intent(inout) :: qql(levs, nvars, londim)
   logical :: forward
   integer, parameter :: otile = 8192   ! 2026-09-18: was 32 (OOM-era); 32 meant ~9.4k tiles x 8 launches per call + an O(nxptot) pack scan per tile (rocprofv3: 634 s of monoadvv pack kernels per run)
   real(kind=RTYPE) :: dd(levs + 1, otile), ds(levs, otile)
   real(kind=RTYPE) :: xgrid(levs + 1, otile), xpast(levs + 1, otile), xnext(levs + 1, otile)
   real(kind=RTYPE) :: da(levs, nvars, otile), step(10), dd_step, dsfact
   integer :: i, j, k, n, nst, nstep, lat, nxj, istr, ob, ot, nloc, i0, i1, ig
   integer :: nnan, ninf, nzero, ib
   real(kind=RTYPE) :: dsmin, ddmax, hostchk, chk
   ! #region agent log
   integer(kind=8) :: dbg0, dbg1, dbg2
   interface
      subroutine geps_dbg_vram(hyp, locid, p0, p1, p2)
         integer hyp, locid
         integer(kind=8) p0, p1, p2
      end subroutine
   end interface
   ! #endregion

   nnan = 0
   ninf = 0
   ib = 1
   if (present(ibase)) ib = ibase
   ! #region agent log
   if (ib .eq. 1) then
      dbg0 = int(londim, 8); dbg1 = int(levs, 8); dbg2 = int(otile, 8)
      call geps_dbg_vram(4, 249, dbg0, dbg1, dbg2)
   end if
   ! #endregion
   !$omp target enter data map(alloc:dd, ds, xgrid, xpast, xnext, da)
   ! #region agent log
   if (ib .eq. 1) then
      dbg0 = int(londim, 8); dbg1 = int(levs, 8); dbg2 = int(nvars, 8)
      call geps_dbg_vram(4, 250, dbg0, dbg1, dbg2)
   end if
   ! #endregion

   do ob = 1, londim, otile
      nloc = min(otile, londim - ob + 1)

      !$omp target teams distribute parallel do collapse(2) firstprivate(nloc)
      do ot = 1, nloc
         do k = 1, levs + 1
            dd(k, ot) = 0.0_RTYPE
         end do
      end do

      !$omp target teams distribute parallel do collapse(2) firstprivate(ob, nloc)
      do ot = 1, nloc
         do k = 1, levs
            ig = ob + ot - 1
            ds(k, ot) = ssi(k, ig) - ssi(k + 1, ig)
         end do
      end do

      !$omp target teams distribute parallel do private(lat, nxj, istr, i0, i1, i, k) firstprivate(ob, nloc, ib)
      do j = 1, jlistnum
         lat = jlist1(j)
         nxj = nxjp(lat)
         istr = nxjp_acc(j) - 1
         i0 = max(1, ib + ob - 1 - istr)
         i1 = min(nxj, ib + ob + nloc - 2 - istr)
         if (i0 .le. i1) then
            do k = 2, levs
               do i = i0, i1
                  dd(k, istr + i - (ib + ob - 1) + 1) = wwi(i, k, j)*deltim
               end do
            end do
         end if
      end do

      ! #region agent log
      if (ob .eq. 1 .and. ib .eq. 1) then
         nnan = 0
         ninf = 0
         !$omp target update from(dd, ds)
         do ot = 1, nloc
            do k = 1, levs + 1
               if (dd(k, ot) .ne. dd(k, ot)) nnan = nnan + 1
               if (abs(dd(k, ot)) .gt. 1.0e20_RTYPE) ninf = ninf + 1
            end do
         end do
         dbg0 = int(nnan, 8)
         dbg1 = int(ninf, 8)
         dbg2 = nint(ds(1, 1)*1.0d6, kind=8)
         call geps_dbg_vram(4, 251, dbg0, dbg1, dbg2)
         nzero = 0
         dsmin = 1.0e30_RTYPE
         do ot = 1, nloc
            do k = 1, levs
               if (ds(k, ot) .eq. 0.0_RTYPE) nzero = nzero + 1
               if (abs(ds(k, ot)) .lt. dsmin) dsmin = abs(ds(k, ot))
            end do
         end do
         dbg0 = int(nzero, 8)
         dbg1 = nint(dsmin*1.0d9, kind=8)
         dbg2 = int(ib, 8)
         call geps_dbg_vram(4, 252, dbg0, dbg1, dbg2)
         ! Recompute on the host the exact quantity def_cfl_step_gpu_type2
         ! computes on the device, from the same dd/ds we just copied back.
         ! host value ~= device value  -> the input really is bad (physics)
         ! host value small, device big -> device/mapping problem
         ddmax = 0.0_RTYPE
         hostchk = 0.0_RTYPE
         do ot = 1, nloc
            do k = 1, levs + 1
               if (abs(dd(k, ot)) .gt. ddmax) ddmax = abs(dd(k, ot))
            end do
            do k = 1, levs
               if (ds(k, ot) .ne. 0.0_RTYPE) then
                  chk = abs((dd(k + 1, ot) - dd(k, ot))/ds(k, ot))
                  if (chk .gt. hostchk) hostchk = chk
               end if
            end do
         end do
         if (ddmax .ne. ddmax .or. ddmax .gt. 1.0d12) then
            dbg0 = -1_8
         else
            dbg0 = nint(ddmax*1.0d3, kind=8)
         end if
         if (hostchk .ne. hostchk .or. hostchk .gt. 1.0d12) then
            dbg1 = -1_8
         else
            dbg1 = nint(hostchk*1.0d3, kind=8)
         end if
         dbg2 = int(nloc, 8)
         call geps_dbg_vram(4, 257, dbg0, dbg1, dbg2)
      end if
      ! #endregion

      call def_cfl_step_gpu_type2(step, nstep, dd, ds, 'advv', &
                                  levs + 1, nloc, async_id)
      ! #region agent log
      if (ob .eq. 1 .and. ib .eq. 1) then
         dbg0 = int(nstep, 8)
         dbg1 = int(nloc, 8)
         dbg2 = int(ib, 8)
         call geps_dbg_vram(4, 253, dbg0, dbg1, dbg2)
      end if
      ! #endregion

      do nst = 1, nstep
         if (forward) then
            !$omp target teams distribute parallel do collapse(2) firstprivate(ob, nloc, nst)
            do ot = 1, nloc
               do k = 1, levs + 1
                  ig = ob + ot - 1
                  dd_step = dd(k, ot)*step(nst)
                  xgrid(k, ot) = ssi(k, ig)
                  xpast(k, ot) = ssi(k, ig)
                  xnext(k, ot) = ssi(k, ig) + dd_step
               end do
            end do
            !$omp target teams distribute parallel do collapse(3) firstprivate(ob, nloc)
            do ot = 1, nloc
               do n = 1, nvars
                  do k = 1, levs
                     ig = ob + ot - 1
                     da(k, n, ot) = qql(k, n, ig)
                  end do
               end do
            end do
         else
            !$omp target teams distribute parallel do collapse(2) firstprivate(ob, nloc, nst)
            do ot = 1, nloc
               do k = 1, levs + 1
                  ig = ob + ot - 1
                  dd_step = dd(k, ot)*step(nst)
                  xgrid(k, ot) = ssi(k, ig)
                  xpast(k, ot) = ssi(k, ig) - dd_step
                  xnext(k, ot) = ssi(k, ig) + dd_step
               end do
            end do
            ! #region agent log
            if (ob .eq. 1 .and. ib .eq. 1 .and. nst .eq. 1) then
               dbg0 = int(nloc, 8); dbg1 = int(nvars, 8); dbg2 = int(nstep, 8)
               call geps_dbg_vram(4, 254, dbg0, dbg1, dbg2)
            end if
            ! #endregion
            call vertical_cell_ppm_intp_gpu(xgrid, qql(1, 1, ob), xpast, da, &
                                             levs, nvars, nloc, async_id)
            ! #region agent log
            if (ob .eq. 1 .and. ib .eq. 1 .and. nst .eq. 1) then
               dbg0 = int(nloc, 8); dbg1 = int(nvars, 8); dbg2 = 1_8
               call geps_dbg_vram(4, 255, dbg0, dbg1, dbg2)
            end if
            ! #endregion
         end if

         if (mass .eq. 1) then
            !$omp target teams distribute parallel do collapse(3) firstprivate(nloc)
            do ot = 1, nloc
               do n = 1, nvars
                  do k = 1, levs
                     dsfact = (xpast(k, ot) - xpast(k + 1, ot)) &
                              /(xnext(k, ot) - xnext(k + 1, ot))
                     da(k, n, ot) = da(k, n, ot)*dsfact
                  end do
               end do
            end do
         end if

         call vertical_cell_ppm_intp_gpu(xnext, da, xgrid, qql(1, 1, ob), &
                                          levs, nvars, nloc, async_id)
         ! #region agent log
         if (ob .eq. 1 .and. ib .eq. 1 .and. nst .eq. 1) then
            dbg0 = int(nloc, 8); dbg1 = int(nvars, 8); dbg2 = 2_8
            call geps_dbg_vram(4, 256, dbg0, dbg1, dbg2)
         end if
         ! #endregion
      end do
   end do

   !$omp target exit data map(delete:dd, ds, xgrid, xpast, xnext, da)
end subroutine vertical_cell_advect_gpu
