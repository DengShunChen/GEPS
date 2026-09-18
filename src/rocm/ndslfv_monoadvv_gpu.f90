! Tile pack/vertical/unpack so we never GPU-create plev(nxptot) or qqlon(nxptot).
! Those maps are ~3.8GB; 1-rank HIP TCo383 only has ~0.5GB after NDSL.

subroutine ndslfv_monoadvv_gpu(ddtemp, qvadv, vdzonl, vdmerd, pdot, &
                               pt, lonsperlat, deltim, forward)
   use param
   use grid, only: latpart, ndslvvar
   use index, only: jlistnum, jlist1, nxp, nxjp_acc, nxptot
   use rank
   use const
   implicit none

   real(kind=RTYPE), intent(inout) :: ddtemp(nxp, lev, my_max), &
                                      qvadv(nxp, lev, ncld, my_max), &
                                      vdmerd(nxp, lev, my_max), &
                                      vdzonl(nxp, lev, my_max)
   real(kind=RTYPE), intent(in) :: pdot(nxp, lev + 1, latpart), pt(nxp, latpart)
   integer, intent(in) :: lonsperlat(my)
   real(kind=RTYPE), intent(in) :: deltim
   integer, parameter :: otile = 8192   ! 2026-09-18: was 32 (OOM-era); 32 meant ~9.4k tiles x 8 launches per call + an O(nxptot) pack scan per tile (rocprofv3: 634 s of monoadvv pack kernels per run)
   real(kind=RTYPE) :: plev_t(lev + 1, otile), qql_t(lev, ndslvvar, otile)
   integer :: mono, mass, i, n, k, kk, lat, lons_lat, istr, j, ob, nloc, ot, ig
   logical :: forward
   integer, parameter :: async_id = 1
   integer :: nnan
   integer(kind=8) :: dbg0, dbg1, dbg2
   interface
      subroutine geps_dbg_vram(hyp, locid, p0, p1, p2)
         integer hyp, locid
         integer(kind=8) p0, p1, p2
      end subroutine
   end interface

   mono = 1
   mass = 0
   dbg0 = int(nxptot, 8); dbg1 = int(ndslvvar, 8); dbg2 = int(otile, 8)
   call geps_dbg_vram(4, 248, dbg0, dbg1, dbg2)
   !$omp target enter data map(alloc:plev_t, qql_t)
   dbg0 = int(nxptot, 8); dbg1 = int(lev, 8); dbg2 = int(otile, 8)
   call geps_dbg_vram(4, 247, dbg0, dbg1, dbg2)

   do ob = 1, nxptot, otile
      nloc = min(otile, nxptot - ob + 1)

      !$omp target teams distribute parallel do collapse(2) firstprivate(nloc)
      do ot = 1, nloc
         do k = 1, lev + 1
            plev_t(k, ot) = 0.0_RTYPE
         end do
      end do

      !$omp target teams distribute parallel do private(lat, lons_lat, istr, ig, k, kk, n) firstprivate(ob, nloc)
      do j = 1, jlistnum
         lat = jlist1(j)
         lons_lat = lonsperlat(lat)
         istr = nxjp_acc(j) - 1
         do i = 1, lons_lat
            ig = istr + i
            if (ig .ge. ob .and. ig .le. ob + nloc - 1) then
               ot = ig - ob + 1
               do k = lev, 1, -1
                  kk = lev - k + 1
                  plev_t(k, ot) = plev_t(k + 1, ot) + dsigma(kk, 1)*pt(i, j) + dsigma(kk, 2)
               end do
               do k = 1, lev
                  kk = lev - k + 1
                  qql_t(k, 1, ot) = vdzonl(i, kk, j)
                  qql_t(k, 2, ot) = vdmerd(i, kk, j)
                  qql_t(k, 3, ot) = ddtemp(i, kk, j)
               end do
               do n = 1, ncld
                  do k = 1, lev
                     kk = lev - k + 1
                     qql_t(k, n + 3, ot) = qvadv(i, kk, n, j)
                  end do
               end do
            end if
         end do
      end do

      call vertical_cell_advect_gpu(nloc, nloc, lev, ndslvvar, &
                                     deltim, plev_t, pdot, qql_t, mass, forward, &
                                     async_id, ob)

      !$omp target teams distribute parallel do private(lat, lons_lat, istr, ig, k, kk, n) firstprivate(ob, nloc)
      do j = 1, jlistnum
         lat = jlist1(j)
         lons_lat = lonsperlat(lat)
         istr = nxjp_acc(j) - 1
         do i = 1, lons_lat
            ig = istr + i
            if (ig .ge. ob .and. ig .le. ob + nloc - 1) then
               ot = ig - ob + 1
               do k = 1, lev
                  kk = lev - k + 1
                  vdzonl(i, kk, j) = qql_t(k, 1, ot)
                  vdmerd(i, kk, j) = qql_t(k, 2, ot)
                  ddtemp(i, kk, j) = qql_t(k, 3, ot)
               end do
               do n = 1, ncld
                  do k = 1, lev
                     kk = lev - k + 1
                     qvadv(i, kk, n, j) = qql_t(k, n + 3, ot)
                  end do
               end do
            end if
         end do
      end do
   end do

   !$omp target exit data map(delete:plev_t, qql_t)
end subroutine ndslfv_monoadvv_gpu

subroutine ndslfv_monoadvv_fgnl_gpu(vdzonl, vdmerd, ddtemp, pdot, &
                                    pt, lonsperlat, deltim, nvars, forward)
   use param
   use grid, only: latpart
   use index, only: jlistnum, jlist1, nxp, nxjp_acc, nxptot
   use rank
   use const
   implicit none

   real(kind=RTYPE), intent(inout) :: ddtemp(nxp, lev, my_max), &
                                      vdmerd(nxp, lev, my_max), &
                                      vdzonl(nxp, lev, my_max)
   real(kind=RTYPE), intent(in) :: pdot(nxp, lev + 1, latpart), pt(nxp, latpart)
   integer, intent(in) :: lonsperlat(my), nvars
   real(kind=RTYPE), intent(in) :: deltim
   integer, parameter :: otile = 8192   ! 2026-09-18: was 32 (OOM-era); 32 meant ~9.4k tiles x 8 launches per call + an O(nxptot) pack scan per tile (rocprofv3: 634 s of monoadvv pack kernels per run)
   real(kind=RTYPE) :: plev_t(lev + 1, otile), qql_t(lev, nvars, otile)
   integer :: mono, mass, i, k, kk, lat, lons_lat, istr, j, ob, nloc, ot, ig
   logical :: forward
   integer, parameter :: async_id = 1
   integer(kind=8) :: dbg0, dbg1, dbg2
   interface
      subroutine geps_dbg_vram(hyp, locid, p0, p1, p2)
         integer hyp, locid
         integer(kind=8) p0, p1, p2
      end subroutine
   end interface

   mono = 1
   mass = 0
   dbg0 = int(nxptot, 8); dbg1 = int(nvars, 8); dbg2 = int(otile, 8)
   call geps_dbg_vram(4, 246, dbg0, dbg1, dbg2)
   !$omp target enter data map(alloc:plev_t, qql_t)
   call geps_dbg_vram(4, 245, dbg0, dbg1, dbg2)

   do ob = 1, nxptot, otile
      nloc = min(otile, nxptot - ob + 1)

      !$omp target teams distribute parallel do collapse(2) firstprivate(nloc)
      do ot = 1, nloc
         do k = 1, lev + 1
            plev_t(k, ot) = 0.0_RTYPE
         end do
      end do

      !$omp target teams distribute parallel do private(lat, lons_lat, istr, ig, k, kk) firstprivate(ob, nloc)
      do j = 1, jlistnum
         lat = jlist1(j)
         lons_lat = lonsperlat(lat)
         istr = nxjp_acc(j) - 1
         do i = 1, lons_lat
            ig = istr + i
            if (ig .ge. ob .and. ig .le. ob + nloc - 1) then
               ot = ig - ob + 1
               do k = lev, 1, -1
                  kk = lev - k + 1
                  plev_t(k, ot) = plev_t(k + 1, ot) + dsigma(kk, 1)*pt(i, j) + dsigma(kk, 2)
               end do
               do k = 1, lev
                  kk = lev - k + 1
                  qql_t(k, 1, ot) = vdzonl(i, kk, j)
                  qql_t(k, 2, ot) = vdmerd(i, kk, j)
                  qql_t(k, 3, ot) = ddtemp(i, kk, j)
               end do
            end if
         end do
      end do

      call vertical_cell_advect_gpu(nloc, nloc, lev, nvars, &
                                     deltim, plev_t, pdot, qql_t, mass, forward, &
                                     async_id, ob)

      !$omp target teams distribute parallel do private(lat, lons_lat, istr, ig, k, kk) firstprivate(ob, nloc)
      do j = 1, jlistnum
         lat = jlist1(j)
         lons_lat = lonsperlat(lat)
         istr = nxjp_acc(j) - 1
         do i = 1, lons_lat
            ig = istr + i
            if (ig .ge. ob .and. ig .le. ob + nloc - 1) then
               ot = ig - ob + 1
               do k = 1, lev
                  kk = lev - k + 1
                  vdzonl(i, kk, j) = qql_t(k, 1, ot)
                  vdmerd(i, kk, j) = qql_t(k, 2, ot)
                  ddtemp(i, kk, j) = qql_t(k, 3, ot)
               end do
            end if
         end do
      end do
   end do

   !$omp target exit data map(delete:plev_t, qql_t)
end subroutine ndslfv_monoadvv_fgnl_gpu
