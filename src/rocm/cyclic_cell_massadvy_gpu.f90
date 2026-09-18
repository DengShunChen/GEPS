! Tiled massadvy: full mylonlen maps are ~4GB (gglati_dup/var/dist/ypast/ynext/dyfact).
! HIP 1-rank TCo383 only has ~1.2GB free after NDSL lon arrays.

subroutine cyclic_cell_massadvy_mylonlen_gpu(latfull, levs, nvars, deltim, vv, qq, mass, forward)
   use grid, only: mylonlen, lonpart, gglati, fa1, fa2, fa3, fa4
   use const, only: RTYPE
   implicit none

   integer, intent(in) :: latfull, levs, nvars, mass
   real(kind=RTYPE), intent(in) :: deltim
   real(kind=RTYPE), dimension(latfull, levs, lonpart), intent(in) :: vv
   real(kind=RTYPE), dimension(latfull, nvars, levs, lonpart) :: qq
   integer, parameter :: otile = 8
   integer :: lon, j, k, nst, n, max_nstep, ot, og, ob, nloc
   real(kind=RTYPE) ds(latfull), sc
   real(kind=RTYPE) var(latfull, levs, otile), dist(latfull + 1, levs, otile)
   real(kind=RTYPE) step(10, levs, otile)
   integer :: nstep(levs, otile)
   real(kind=RTYPE) ypast(latfull + 1, levs, otile), ynext(latfull + 1, levs, otile)
   real(kind=RTYPE) dyfact(latfull, levs, otile), dist_step
   integer :: idx11, idx12, idx21, idx22, idxfa
   integer :: outer_index(3, otile), outer_index_def(otile)
   logical :: nstep_less(levs, otile)
   real(kind=RTYPE) gglati_dup(latfull + 1, levs, otile)
   real(kind=RTYPE) ds_dup(latfull + 1, otile)
   logical :: forward
   real(kind=RTYPE) :: vv_t
   integer :: nnan, ninf, jchk
   ! #region agent log
   integer(kind=8) :: dbg0, dbg1, dbg2
   interface
      subroutine geps_dbg_vram(hyp, locid, p0, p1, p2)
         integer hyp, locid
         integer(kind=8) p0, p1, p2
      end subroutine
   end interface
   ! #endregion

   sc = gglati(latfull + 1) - gglati(1)

   do j = 1, latfull
      ds(j) = gglati(j + 1) - gglati(j)
   end do

   ! #region agent log
   dbg0 = int(latfull, 8); dbg1 = int(nvars, 8); dbg2 = int(otile, 8)
   call geps_dbg_vram(2, 240, dbg0, dbg1, dbg2)
   ! #endregion
   !$omp target enter data map(to:ds)
   !$omp target enter data map(alloc:gglati_dup, var, dist, outer_index_def, ds_dup, step, nstep, ypast, ynext, dyfact, outer_index, nstep_less)
   ! #region agent log
   dbg0 = int(latfull, 8); dbg1 = int(nvars, 8); dbg2 = int(mylonlen, 8)
   call geps_dbg_vram(2, 241, dbg0, dbg1, dbg2)
   ! #endregion

   do ob = 1, mylonlen, otile
      nloc = min(otile, mylonlen - ob + 1)

      !$omp target teams distribute parallel do collapse(2) private(og) firstprivate(ob, nloc)
      do ot = 1, nloc
         do k = 1, levs
            og = ob + ot - 1
            gglati_dup(:, k, ot) = gglati
         end do
      end do

      !$omp target teams distribute parallel do collapse(3) private(og, vv_t) firstprivate(ob, nloc)
      do ot = 1, nloc
         do k = 1, levs
            do j = 1, latfull
               og = ob + ot - 1
               vv_t = vv(j, k, og)*deltim
               if (j .le. latfull/2) then
                  var(j, k, ot) = vv_t
               else
                  var(j, k, ot) = -vv_t
               end if
            end do
         end do
      end do

      !$omp target teams distribute parallel do collapse(3) private(idx11, idx12, idx21, idx22, idxfa) firstprivate(nloc)
      do ot = 1, nloc
         do k = 1, levs
            do j = 1, latfull + 1
               idx11 = mod(j - 3 + latfull, latfull) + 1
               idx12 = mod(j - 2 + latfull, latfull) + 1
               idx21 = mod(j - 1, latfull) + 1
               idx22 = mod(j, latfull) + 1
               idxfa = mod(j - 1, latfull) + 1
               dist(j, k, ot) = fa1(idxfa)*var(idx11, k, ot) &
                                + fa2(idxfa)*var(idx12, k, ot) &
                                + fa3(idxfa)*var(idx21, k, ot) &
                                + fa4(idxfa)*var(idx22, k, ot)
            end do
         end do
      end do

      !$omp target teams distribute parallel do firstprivate(nloc)
      do ot = 1, nloc
         outer_index_def(ot) = latfull
         outer_index(1, ot) = latfull
         outer_index(2, ot) = latfull
         outer_index(3, ot) = latfull
         ds_dup(1:latfull, ot) = ds
      end do

      call def_cfl_step_two_loops_gpu(outer_index_def, nloc, levs, dist, ds_dup, step, nstep, max_nstep, 'advy', latfull)
      ! #region agent log
      if (ob .eq. 1) then
         dbg0 = int(max_nstep, 8)
         dbg1 = nint(ds(1)*1.0d6, kind=8)
         dbg2 = nint(ds(max(1, latfull/2))*1.0d6, kind=8)
         call geps_dbg_vram(3, 243, dbg0, dbg1, dbg2)
      end if
      ! #endregion

      do nst = 1, max_nstep
         !$omp target teams distribute parallel do collapse(2) firstprivate(nloc, nst)
         do ot = 1, nloc
            do k = 1, levs
               nstep_less(k, ot) = (nst .le. nstep(k, ot))
            end do
         end do
         !$omp target teams distribute parallel do collapse(3) private(dist_step) firstprivate(nloc, nst, forward)
         do ot = 1, nloc
            do k = 1, levs
               do j = 1, latfull + 1
                  if (nstep_less(k, ot)) then
                     dist_step = dist(j, k, ot)*step(nst, k, ot)
                     if (forward) then
                        ypast(j, k, ot) = gglati(j)
                     else
                        ypast(j, k, ot) = gglati(j) - dist_step
                     end if
                     ynext(j, k, ot) = gglati(j) + dist_step
                  end if
               end do
            end do
         end do
         if (mass .eq. 1) then
            !$omp target teams distribute parallel do collapse(3) firstprivate(nloc)
            do ot = 1, nloc
               do k = 1, levs
                  do j = 1, latfull
                     if (nstep_less(k, ot)) then
                        dyfact(j, k, ot) = (ypast(j + 1, k, ot) - ypast(j, k, ot))/(ynext(j + 1, k, ot) - ynext(j, k, ot))
                     end if
                  end do
               end do
            end do
         end if
         if (.not. forward) then
            call cyclic_cell_ppm_intp_two_loops_gpu(outer_index, nloc, levs, &
                                                    gglati_dup, ypast, qq(1, 1, 1, ob), &
                                                    latfull, nvars, sc, nstep_less)
         end if
         if (mass .eq. 1) then
            !$omp target teams distribute parallel do collapse(4) private(og) firstprivate(ob, nloc)
            do ot = 1, nloc
               do k = 1, levs
                  do n = 1, nvars
                     do j = 1, latfull
                        og = ob + ot - 1
                        if (nstep_less(k, ot)) then
                           qq(j, n, k, og) = qq(j, n, k, og)*dyfact(j, k, ot)
                        end if
                     end do
                  end do
               end do
            end do
         end if
         call cyclic_cell_ppm_intp_two_loops_gpu(outer_index, nloc, levs, &
                                                 ynext, gglati_dup, qq(1, 1, 1, ob), &
                                                 latfull, nvars, sc, nstep_less)
      end do
   end do

   ! #region agent log
   nnan = 0
   ninf = 0
   !$omp target update from(qq(1:latfull, 1, 1, 1))
   do jchk = 1, latfull
      if (qq(jchk, 1, 1, 1) .ne. qq(jchk, 1, 1, 1)) nnan = nnan + 1
      if (abs(qq(jchk, 1, 1, 1)) .gt. 1.0e20_RTYPE) ninf = ninf + 1
   end do
   dbg0 = int(nnan, 8)
   dbg1 = int(ninf, 8)
   dbg2 = nint(sc*1.0d6, kind=8)
   call geps_dbg_vram(3, 242, dbg0, dbg1, dbg2)
   ! #endregion

   !$omp target exit data map(delete:ds, gglati_dup, var, dist, outer_index_def, ds_dup, step, nstep, ypast, ynext, dyfact, outer_index, nstep_less)
end subroutine cyclic_cell_massadvy_mylonlen_gpu
