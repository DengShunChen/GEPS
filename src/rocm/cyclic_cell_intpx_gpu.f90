! Tiled intpx: nvidia/ndslfv_pack_gpu maps full old+new (~2*qqlon, ~14GB) plus
! xpast/xnext. HIP 1-rank TCo383 only has ~2GB free after massadvx. `new` is unused.

subroutine cyclic_cell_intpx_jlist_gpu(levs, nvars, lonfull, qq, to_lonfull)
   use index, only: jlistnum, jlist1, nxdef
   use const, only: RTYPE
   use grid, only: latpart
   implicit none

   integer, intent(in) :: levs, nvars, lonfull
   real(kind=RTYPE), dimension(lonfull, nvars, levs, latpart) :: qq
   logical, intent(in) :: to_lonfull
   integer :: lan, lat, lons_lat, i, imp, imf, k, n, ot, og, ob, nloc
   real(kind=RTYPE) :: pi, two_pi
   integer, parameter :: otile = 8
   real(kind=RTYPE), dimension(lonfull + 1, levs, otile) :: xpast, xnext
   real(kind=RTYPE), dimension(lonfull, nvars, levs, otile) :: old
   integer :: async_id
   integer :: outer_index(3, otile)
   logical :: nstep_less(levs, otile)
   integer :: lons_size(jlistnum)

   ! #region agent log
   integer(kind=8) :: dbg0, dbg1, dbg2
   interface
      subroutine geps_dbg_vram(hyp, locid, p0, p1, p2)
         integer hyp, locid
         integer(kind=8) p0, p1, p2
      end subroutine
   end interface
   ! #endregion

   async_id = 1
   pi = 4.0*atan(1.0)
   two_pi = 2.0*pi

   do i = 1, jlistnum
      lat = jlist1(i)
      lons_size(i) = nxdef(lat)
   end do

   ! #region agent log
   dbg0 = int(lonfull, 8); dbg1 = int(nvars, 8); dbg2 = int(otile, 8)
   call geps_dbg_vram(2, 32, dbg0, dbg1, dbg2)
   ! #endregion
   !$omp target enter data map(alloc:xpast, xnext, old, outer_index, nstep_less)
   !$omp target enter data map(to:lons_size)
   ! #region agent log
   dbg0 = int(lonfull, 8); dbg1 = int(nvars, 8); dbg2 = int(jlistnum, 8)
   call geps_dbg_vram(2, 33, dbg0, dbg1, dbg2)
   ! #endregion

   do ob = 1, jlistnum, otile
      nloc = min(otile, jlistnum - ob + 1)

      !$omp target teams distribute parallel do collapse(2) private(lons_lat, imp, imf, og) firstprivate(ob, nloc, to_lonfull)
      do ot = 1, nloc
         do i = 1, lonfull + 1
            og = ob + ot - 1
            lons_lat = lons_size(og)
            if (to_lonfull) then
               imp = lons_lat
               imf = lonfull
            else
               imp = lonfull
               imf = lons_lat
            end if
            if (i .le. imp + 1) then
               xpast(i, :, ot) = real(i - 1.5)*two_pi/imp
            end if
            if (i .le. imf + 1) then
               xnext(i, :, ot) = real(i - 1.5)*two_pi/imf
            end if
         end do
      end do

      !$omp target teams distribute parallel do collapse(3) private(lons_lat, og) firstprivate(ob, nloc)
      do ot = 1, nloc
         do k = 1, levs
            do n = 1, nvars
               og = ob + ot - 1
               lons_lat = lons_size(og)
               if (lons_lat .eq. lonfull) then
                  do i = 1, lons_lat
                     old(i, n, k, ot) = qq(i, n, k, og)
                  end do
               end if
            end do
         end do
      end do

      !$omp target teams distribute parallel do private(lons_lat, imp, imf, og) firstprivate(ob, nloc, to_lonfull)
      do ot = 1, nloc
         og = ob + ot - 1
         lons_lat = lons_size(og)
         if (to_lonfull) then
            imp = lons_lat
            imf = lonfull
         else
            imp = lonfull
            imf = lons_lat
         end if
         outer_index(1, ot) = lonfull
         outer_index(2, ot) = imp
         outer_index(3, ot) = imf
      end do

      !$omp target
      nstep_less = .true.
      !$omp end target

      call cyclic_cell_ppm_intp_two_loops_gpu(outer_index, nloc, levs, &
                                              xpast, xnext, qq(1, 1, 1, ob), &
                                              lonfull, nvars, two_pi, nstep_less)

      !$omp target teams distribute parallel do collapse(3) private(lons_lat, og) firstprivate(ob, nloc)
      do ot = 1, nloc
         do k = 1, levs
            do n = 1, nvars
               og = ob + ot - 1
               lons_lat = lons_size(og)
               if (lons_lat .eq. lonfull) then
                  do i = 1, lons_lat
                     qq(i, n, k, og) = old(i, n, k, ot)
                  end do
               end if
            end do
         end do
      end do
   end do

   !$omp target exit data map(delete:xpast, xnext, old, outer_index, nstep_less, lons_size)
end subroutine cyclic_cell_intpx_jlist_gpu
