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
   integer, parameter :: otile = 32768   ! 2026-09-19: 8192 -> 32768 (131072 was slower) (4x fewer tiles, prof11: NDSL launch-gap bound); 2026-09-18: was 32 (OOM-era); 32 meant ~9.4k tiles x 8 launches per call + an O(nxptot) pack scan per tile (rocprofv3: 634 s of monoadvv pack kernels per run)
   real(kind=RTYPE) :: plev_t(lev + 1, otile), qql_t(lev, ndslvvar, otile)
   integer :: col_i(nxptot), col_j(nxptot)   ! column -> (i, j) map for tile-local pack/unpack
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
   !$omp target enter data map(alloc:col_i, col_j)
   !$omp target teams distribute parallel do private(lat, lons_lat, istr, i)
   do j = 1, jlistnum
      lat = jlist1(j)
      lons_lat = lonsperlat(lat)
      istr = nxjp_acc(j) - 1
      do i = 1, lons_lat
         col_i(istr + i) = i
         col_j(istr + i) = j
      end do
   end do
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

      ! 2026-09-18: tile-local pack via the column map (was a scan over ALL
      ! columns per tile: O(nxptot) per launch, 159 ms x 37 tiles per call).
      !$omp target teams distribute parallel do private(ig, i, j, k, kk) firstprivate(ob, nloc)
      do ot = 1, nloc
         ig = ob + ot - 1
         i = col_i(ig)
         j = col_j(ig)
         do k = lev, 1, -1
            kk = lev - k + 1
            plev_t(k, ot) = plev_t(k + 1, ot) + dsigma(kk, 1)*pt(i, j) + dsigma(kk, 2)
         end do
      end do
      !$omp target teams distribute parallel do collapse(2) private(ig, i, j, kk, n) firstprivate(ob, nloc)
      do ot = 1, nloc
         do k = 1, lev
            ig = ob + ot - 1
            i = col_i(ig)
            j = col_j(ig)
            kk = lev - k + 1
            qql_t(k, 1, ot) = vdzonl(i, kk, j)
            qql_t(k, 2, ot) = vdmerd(i, kk, j)
            qql_t(k, 3, ot) = ddtemp(i, kk, j)
            do n = 1, ncld
               qql_t(k, n + 3, ot) = qvadv(i, kk, n, j)
            end do
         end do
      end do

      call vertical_cell_advect_gpu(nloc, nloc, lev, ndslvvar, &
                                     deltim, plev_t, pdot, qql_t, mass, forward, &
                                     async_id, ob)

      !$omp target teams distribute parallel do collapse(2) private(ig, i, j, kk, n) firstprivate(ob, nloc)
      do ot = 1, nloc
         do k = 1, lev
            ig = ob + ot - 1
            i = col_i(ig)
            j = col_j(ig)
            kk = lev - k + 1
            vdzonl(i, kk, j) = qql_t(k, 1, ot)
            vdmerd(i, kk, j) = qql_t(k, 2, ot)
            ddtemp(i, kk, j) = qql_t(k, 3, ot)
            do n = 1, ncld
               qvadv(i, kk, n, j) = qql_t(k, n + 3, ot)
            end do
         end do
      end do
   end do

   !$omp target exit data map(delete:plev_t, qql_t)
   !$omp target exit data map(delete:col_i, col_j)
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
   integer, parameter :: otile = 32768   ! 2026-09-19: 8192 -> 32768 (131072 was slower) (4x fewer tiles, prof11: NDSL launch-gap bound); 2026-09-18: was 32 (OOM-era); 32 meant ~9.4k tiles x 8 launches per call + an O(nxptot) pack scan per tile (rocprofv3: 634 s of monoadvv pack kernels per run)
   real(kind=RTYPE) :: plev_t(lev + 1, otile), qql_t(lev, nvars, otile)
   integer :: col_i(nxptot), col_j(nxptot)   ! column -> (i, j) map for tile-local pack/unpack
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
   !$omp target enter data map(alloc:col_i, col_j)
   !$omp target teams distribute parallel do private(lat, lons_lat, istr, i)
   do j = 1, jlistnum
      lat = jlist1(j)
      lons_lat = lonsperlat(lat)
      istr = nxjp_acc(j) - 1
      do i = 1, lons_lat
         col_i(istr + i) = i
         col_j(istr + i) = j
      end do
   end do
   call geps_dbg_vram(4, 245, dbg0, dbg1, dbg2)

   do ob = 1, nxptot, otile
      nloc = min(otile, nxptot - ob + 1)

      !$omp target teams distribute parallel do collapse(2) firstprivate(nloc)
      do ot = 1, nloc
         do k = 1, lev + 1
            plev_t(k, ot) = 0.0_RTYPE
         end do
      end do

      !$omp target teams distribute parallel do private(ig, i, j, k, kk) firstprivate(ob, nloc)
      do ot = 1, nloc
         ig = ob + ot - 1
         i = col_i(ig)
         j = col_j(ig)
         do k = lev, 1, -1
            kk = lev - k + 1
            plev_t(k, ot) = plev_t(k + 1, ot) + dsigma(kk, 1)*pt(i, j) + dsigma(kk, 2)
         end do
      end do
      !$omp target teams distribute parallel do collapse(2) private(ig, i, j, kk) firstprivate(ob, nloc)
      do ot = 1, nloc
         do k = 1, lev
            ig = ob + ot - 1
            i = col_i(ig)
            j = col_j(ig)
            kk = lev - k + 1
            qql_t(k, 1, ot) = vdzonl(i, kk, j)
            qql_t(k, 2, ot) = vdmerd(i, kk, j)
            qql_t(k, 3, ot) = ddtemp(i, kk, j)
         end do
      end do

      call vertical_cell_advect_gpu(nloc, nloc, lev, nvars, &
                                     deltim, plev_t, pdot, qql_t, mass, forward, &
                                     async_id, ob)

      !$omp target teams distribute parallel do collapse(2) private(ig, i, j, kk) firstprivate(ob, nloc)
      do ot = 1, nloc
         do k = 1, lev
            ig = ob + ot - 1
            i = col_i(ig)
            j = col_j(ig)
            kk = lev - k + 1
            vdzonl(i, kk, j) = qql_t(k, 1, ot)
            vdmerd(i, kk, j) = qql_t(k, 2, ot)
            ddtemp(i, kk, j) = qql_t(k, 3, ot)
         end do
      end do
   end do

   !$omp target exit data map(delete:plev_t, qql_t)
   !$omp target exit data map(delete:col_i, col_j)
end subroutine ndslfv_monoadvv_fgnl_gpu
