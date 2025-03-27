subroutine cyclic_cell_intpx_jlist_gpu(levs, nvars, lonfull, qq, to_lonfull)
   ! Present on device: qq, jlist1, nxdef
   use index, only: jlistnum, jlist1, nxdef
   use const, only: RTYPE
   use grid, only: latpart

   implicit none

   integer, intent(in) :: levs, nvars, lonfull
   real(kind=RTYPE), dimension(lonfull, nvars, levs, latpart) :: qq
   logical, intent(in) :: to_lonfull
   integer :: lan, lat, lons_lat, i, j, imp, imf, nlevs, k, n
   real(kind=RTYPE) :: pi, two_pi
   real(kind=RTYPE), dimension(lonfull + 1, levs, jlistnum) :: xpast, xnext
   real(kind=RTYPE), dimension(lonfull, nvars, levs, jlistnum) :: old, new
   integer :: async_id

   integer :: outer_index(3, jlistnum)
   logical :: nstep_less(levs, jlistnum)
   integer :: lons_size(jlistnum)

   async_id = 1
   nlevs = levs*nvars
   pi = 4.0*atan(1.0)
   two_pi = 2.0*pi

   do i = 1, jlistnum
      lat = jlist1(i)
      lons_size(i) = nxdef(lat)
   end do

   !$acc enter data create(xpast, xnext, old, outer_index, new, nstep_less) async(async_id)
   !$acc enter data copyin(lons_size) async(async_id)

   !$acc parallel loop collapse(2) private(lat, lons_lat, imp, imf) async(async_id)
   do lan = 1, jlistnum
      do i = 1, lonfull + 1
         lons_lat = lons_size(lan)
         if (to_lonfull) then
            imp = lons_lat
            imf = lonfull
         else
            imp = lonfull
            imf = lons_lat
         end if
         if (i .le. imp + 1) then
            xpast(i, :, lan) = dble(i - 1.5)*two_pi/imp
         end if
         if (i .le. imf + 1) then
            xnext(i, :, lan) = dble(i - 1.5)*two_pi/imf
         end if
      end do
   end do

   !$acc parallel loop collapse(3) private(lat, lons_lat) async(async_id)
   do lan = 1, jlistnum
      do k = 1, levs
         do n = 1, nvars
            lons_lat = lons_size(lan)
            if (lons_lat .eq. lonfull) then
               do i = 1, lons_lat
                  old(i, n, k, lan) = qq(i, n, k, lan)
               end do
            end if
         end do
      end do
   end do

   !$acc parallel loop private(lat, lons_lat, imp, imf) async(async_id)
   do lan = 1, jlistnum
      lons_lat = lons_size(lan)
      if (to_lonfull) then
         imp = lons_lat
         imf = lonfull
      else
         imp = lonfull
         imf = lons_lat
      end if
      outer_index(1, lan) = lonfull
      outer_index(2, lan) = imp
      outer_index(3, lan) = imf
   end do

   !$acc kernels async(async_id)
   nstep_less = .true.
   !$acc end kernels
   call cyclic_cell_ppm_intp_two_loops_gpu(outer_index, jlistnum, levs, &
                                           xpast, xnext, qq, &
                                           lonfull, nvars, two_pi, nstep_less)

   !$acc parallel loop collapse(3) private(lat, lons_lat) async(async_id)
   do lan = 1, jlistnum
      do k = 1, levs
         do n = 1, nvars
            lons_lat = lons_size(lan)
            if (lons_lat .eq. lonfull) then
               do i = 1, lons_lat
                  qq(i, n, k, lan) = old(i, n, k, lan)
               end do
            end if
         end do
      end do
   end do
   !$acc exit data delete(xpast, xnext, old, outer_index, new, nstep_less, lons_size) async(async_id)

end subroutine cyclic_cell_intpx_jlist_gpu

subroutine cyclic_cell_massadvx_jlist_gpu(levs, nvars, lonfull, deltim, uu, qq, mass, forward)
   ! Present on device: uu, qq, jlist1, nxdef
   use index, only: jlistnum, jlist1, nxdef
   use const, only: RTYPE
   use grid, only: latpart

   implicit none

   integer, intent(in) :: levs, nvars, lonfull, mass
   real(kind=RTYPE), intent(in) :: deltim
   real(kind=RTYPE), intent(in), dimension(lonfull, levs, latpart) :: uu
   real(kind=RTYPE), intent(out), dimension(lonfull, nvars, levs, latpart) :: qq
   integer :: lan, lat, lons_lat, im, i, k
   real(kind=RTYPE), dimension(lonfull + 1, jlistnum) :: ds, xreg
   real(kind=RTYPE) :: ds_t
   real(kind=RTYPE) :: pi, sc, dist(lonfull + 1, levs, jlistnum)
   real(kind=RTYPE), parameter :: fa1 = 9./16.
   real(kind=RTYPE), parameter :: fa2 = 1./16.
   real(kind=RTYPE), dimension(lonfull + 1, levs, jlistnum) :: xpast, xnext
   real(kind=RTYPE) :: step(10, levs, jlistnum), dist_step, dxfact(lonfull, levs, jlistnum)
   integer :: nstep(levs, jlistnum), nst, n, idx11, idx12, idx21, idx22, max_nstep
   integer :: async_id
   integer :: outer_index(3, jlistnum)
   logical :: nstep_less(levs, jlistnum)
   real(kind=RTYPE) :: xreg_dup(lonfull + 1, levs, jlistnum)
   logical :: forward
   integer :: nstep_max
   integer :: lons_size(jlistnum)

   async_id = 1

   pi = 4.0*atan(1.0)
   sc = 2.0*pi

   do lan = 1, jlistnum
      lat = jlist1(lan)
      lons_size(lan) = nxdef(lat)
   end do

   !$acc enter data create(ds, xreg, dist, step, nstep, &
   !$acc& xpast, xnext, dxfact, outer_index, xreg_dup, nstep_less) async(async_id)
   !$acc enter data copyin(lons_size) async(async_id)
   !$acc parallel loop gang collapse(2) private(lons_lat) async(async_id)
   do lan = 1, jlistnum
      do k = 1, levs
         lons_lat = lons_size(lan)
         !$acc loop vector private(idx11, idx12, idx21, idx22, ds_t)
         do i = 1, lons_lat + 1
            if (k .eq. 1) then
               ds_t = sc/float(lons_lat)
               ds(i, lan) = ds_t
               xreg(i, lan) = (i - 1.5)*ds_t
            end if
            idx11 = mod(i - 1 + lons_lat, lons_lat) + 1
            idx12 = mod(i - 2 + lons_lat, lons_lat) + 1
            idx21 = mod(i + lons_lat, lons_lat) + 1
            idx22 = mod(i - 3 + lons_lat, lons_lat) + 1
            dist(i, k, lan) = (fa1*(uu(idx11, k, lan) + uu(idx12, k, lan)) - fa2*(uu(idx21, k, lan) + uu(idx22, k, lan)))*deltim
         end do
      end do
   end do
   call def_cfl_step_two_loops_gpu(lons_size, jlistnum, levs, dist, ds, step, nstep, max_nstep, 'advx', lonfull)
   !$acc parallel loop async(async_id)
   do i = 1, 3
      outer_index(i, :) = lons_size
   end do
   !$acc parallel loop async(async_id)
   do k = 1, levs
      xreg_dup(:, k, :) = xreg
   end do
   do nst = 1, max_nstep
      !$acc kernels async(async_id)
      nstep_less = (nst .le. nstep)
      !$acc end kernels
      !$acc parallel loop gang collapse(2) private(lons_lat) async(async_id)
      do lan = 1, jlistnum
         do k = 1, levs
            if (nstep_less(k, lan)) then
               lons_lat = lons_size(lan)
               !$acc loop vector private(dist_step)
               do i = 1, lons_lat + 1
                  dist_step = dist(i, k, lan)*step(nst, k, lan)
                  if (forward) then
                     xpast(i, k, lan) = xreg(i, lan)
                  else
                     xpast(i, k, lan) = xreg(i, lan) - dist_step
                  end if
                  xnext(i, k, lan) = xreg(i, lan) + dist_step
               end do
            end if
         end do
      end do
      if (mass .eq. 1) then
         !$acc parallel loop gang collapse(2) private(lons_lat) async(async_id)
         do lan = 1, jlistnum
            do k = 1, levs
               if (nstep_less(k, lan)) then
                  lons_lat = lons_size(lan)
                  !$acc loop vector
                  do i = 1, lons_lat
                     dxfact(i, k, lan) = (xpast(i + 1, k, lan) - xpast(i, k, lan)) &
                                         /(xnext(i + 1, k, lan) - xnext(i, k, lan))
                  end do
               end if
            end do
         end do
      end if
      if (.not. forward) then
         call cyclic_cell_ppm_intp_two_loops_gpu(outer_index, jlistnum, levs, &
                                                 xreg_dup, xpast, qq, &
                                                 lonfull, nvars, sc, nstep_less)
      end if
      if (mass .eq. 1) then
         !$acc parallel loop gang collapse(3) private(lat, lons_lat) async(async_id)
         do lan = 1, jlistnum
            do k = 1, levs
               do n = 1, nvars
                  if (nstep_less(k, lan)) then
                     lons_lat = lons_size(lan)
                     !$acc loop vector
                     do i = 1, lons_lat
                        qq(i, n, k, lan) = qq(i, n, k, lan)*dxfact(i, k, lan)
                     end do
                  end if
               end do
            end do
         end do
      end if
      call cyclic_cell_ppm_intp_two_loops_gpu(outer_index, jlistnum, levs, &
                                              xnext, xreg_dup, qq, &
                                              lonfull, nvars, sc, nstep_less)
   end do
   !$acc exit data delete(ds, xreg, dist, step, nstep, &
   !$acc& xpast, xnext, dxfact, outer_index, xreg_dup, nstep_less, lons_size) async(async_id)

end subroutine cyclic_cell_massadvx_jlist_gpu

subroutine cyclic_cell_massadvy_mylonlen_gpu(latfull, levs, nvars, deltim, vv, qq, mass, forward)
   ! Present on device: vv, qq, gglati, fa1, fa2, fa3, fa4
   use grid, only: mylonlen, lonpart, gglati, fa1, fa2, fa3, fa4
   use const, only: RTYPE

   implicit none

   integer, intent(in) :: latfull, levs, nvars, mass
   real(kind=RTYPE), intent(in) :: deltim
   real(kind=RTYPE), dimension(latfull, levs, lonpart), intent(in) :: vv
   real(kind=RTYPE), dimension(latfull, nvars, levs, lonpart), intent(out) :: qq
   integer :: lon, j, k, jm, nstep(levs, mylonlen), nst, n, max_nstep
   real(kind=RTYPE) ds(latfull), var(latfull, levs, mylonlen), dist(latfull + 1, levs, mylonlen)
   real(kind=RTYPE) step(10, levs, mylonlen), sc
   real(kind=RTYPE) ypast(latfull + 1, levs, mylonlen), ynext(latfull + 1, levs, mylonlen)
   real(kind=RTYPE) dyfact(latfull, levs, mylonlen), dist_step
   integer :: idx11, idx12, idx21, idx22, idxfa
   integer :: async_id
   integer :: outer_index(3, mylonlen), nstep_max
   logical :: nstep_less(levs, mylonlen)
   real(kind=RTYPE) gglati_dup(latfull + 1, levs, mylonlen)
   integer :: outer_index_def(mylonlen)
   real(kind=RTYPE) ds_dup(latfull + 1, mylonlen)
   logical :: forward
   real(kind=RTYPE) :: vv_t

   async_id = 1

   sc = gglati(latfull + 1) - gglati(1)

   !$acc enter data create(gglati_dup, ds, var, dist, &
   !$acc& outer_index_def, ds_dup, step, nstep, ypast, ynext, dyfact, &
   !$acc& outer_index, nstep_less) async(async_id)

   !$acc parallel loop collapse(2) async(async_id)
   do lon = 1, mylonlen
      do k = 1, levs
         gglati_dup(:, k, lon) = gglati
      end do
   end do

   !$acc parallel loop async(async_id)
   do j = 1, latfull
      ds(j) = gglati(j + 1) - gglati(j)
   end do

   !$acc parallel loop collapse(3) async(async_id)
   do lon = 1, mylonlen
      do k = 1, levs
         do j = 1, latfull
            vv_t = vv(j, k, lon)*deltim
            if (j .le. latfull/2) then
               var(j, k, lon) = vv_t
            else
               var(j, k, lon) = -vv_t
            end if
         end do
      end do
   end do

   !$acc parallel loop collapse(3) private(idx11, idx12, idx21, idx22, idxfa) async(async_id)
   do lon = 1, mylonlen
      do k = 1, levs
         do j = 1, latfull + 1
            idx11 = mod(j - 3 + latfull, latfull) + 1
            idx12 = mod(j - 2 + latfull, latfull) + 1
            idx21 = mod(j - 1, latfull) + 1
            idx22 = mod(j, latfull) + 1
            idxfa = mod(j - 1, latfull) + 1
            dist(j, k, lon) = fa1(idxfa)*var(idx11, k, lon) &
                              + fa2(idxfa)*var(idx12, k, lon) &
                              + fa3(idxfa)*var(idx21, k, lon) &
                              + fa4(idxfa)*var(idx22, k, lon)
         end do
      end do
   end do

   !$acc kernels async(async_id)
   outer_index_def = latfull
   !$acc end kernels
   !$acc kernels async(async_id)
   outer_index = latfull
   !$acc end kernels
   !$acc parallel loop async(async_id)
   do lon = 1, mylonlen
      ds_dup(1:latfull, lon) = ds
   end do
   call def_cfl_step_two_loops_gpu(outer_index_def, mylonlen, levs, dist, ds_dup, step, nstep, max_nstep, 'advy', latfull)
   do nst = 1, max_nstep
      !$acc kernels async(async_id)
      nstep_less = (nst .le. nstep)
      !$acc end kernels
      !$acc parallel loop collapse(3) private(dist_step) async(async_id)
      do lon = 1, mylonlen
         do k = 1, levs
            do j = 1, latfull + 1
               if (nstep_less(k, lon)) then
                  dist_step = dist(j, k, lon)*step(nst, k, lon)
                  if (forward) then
                     ypast(j, k, lon) = gglati(j)
                  else
                     ypast(j, k, lon) = gglati(j) - dist_step
                  end if
                  ynext(j, k, lon) = gglati(j) + dist_step
               end if
            end do
         end do
      end do
      if (mass .eq. 1) then
         !$acc parallel loop collapse(3) async(async_id)
         do lon = 1, mylonlen
            do k = 1, levs
               do j = 1, latfull
                  if (nstep_less(k, lon)) then
                     dyfact(j, k, lon) = (ypast(j + 1, k, lon) - ypast(j, k, lon))/(ynext(j + 1, k, lon) - ynext(j, k, lon))
                  end if
               end do
            end do
         end do
      end if
      if (.not. forward) then
         call cyclic_cell_ppm_intp_two_loops_gpu(outer_index, mylonlen, levs, &
                                                 gglati_dup, ypast, qq, &
                                                 latfull, nvars, sc, nstep_less)
      end if
      if (mass .eq. 1) then
         !$acc parallel loop collapse(4) async(async_id)
         do lon = 1, mylonlen
            do k = 1, levs
               do n = 1, nvars
                  do j = 1, latfull
                     if (nstep_less(k, lon)) then
                        qq(j, n, k, lon) = qq(j, n, k, lon)*dyfact(j, k, lon)
                     end if
                  end do
               end do
            end do
         end do
      end if
      call cyclic_cell_ppm_intp_two_loops_gpu(outer_index, mylonlen, levs, &
                                              ynext, gglati_dup, qq, &
                                              latfull, nvars, sc, nstep_less)
   end do
   !$acc exit data delete(gglati_dup, ds, var, dist, &
   !$acc& outer_index_def, ds_dup, step, nstep, ypast, ynext, dyfact, &
   !$acc& outer_index, nstep_less) async(async_id)

end subroutine cyclic_cell_massadvy_mylonlen_gpu

subroutine def_cfl_step_two_loops_gpu(outer_index, outer_size, inner_size, dist, ds, step, nstep, max_nstep, job, lonn)
   use const, only: RTYPE

   implicit none

   integer :: outer_index(outer_size), outer_size, inner_size, lonn
   real(kind=RTYPE) :: dist(lonn + 1, inner_size, outer_size)
   real(kind=RTYPE) :: ds(lonn + 1, outer_size)
   real(kind=RTYPE) :: step(10, inner_size, outer_size)
   integer :: nstep(inner_size, outer_size), max_nstep, max_nstep_loc
   character*4 :: job
   integer :: outer, inner, im, n, k, nchk(inner_size, outer_size), nchk_local, nstep_loc
   real(kind=RTYPE) :: rstep, check, check_max(inner_size, outer_size), &
                       check_point, safe_step, last_step, check_local, check_max_local, nchk_max_local
   integer :: async_id

   async_id = 1

   check_point = 1.00
   safe_step = 0.99
   max_nstep = 1
   !$acc enter data create(nchk, check_max) async(async_id)

   !$acc parallel loop collapse(2) async(async_id)
   do outer = 1, outer_size
      do inner = 1, inner_size
         nstep(inner, outer) = 1
         step(1, inner, outer) = 1.0
      end do
   end do
   !$acc parallel loop collapse(2) async(async_id)
   do outer = 1, outer_size
      do inner = 1, inner_size
         nchk_local = 0
         check_local = 0.
         !$acc loop private(im, check) reduction(max:nchk_local, check_local)
         do n = 1, lonn
            im = outer_index(outer)
            if (n .le. im) then
               check = abs((dist(n + 1, inner, outer) - dist(n, inner, outer))/ds(n, outer))
               if (check .ge. check_point) then
                  nchk_local = max(nchk_local, n)
                  check_local = max(check_local, check)
               end if
            end if
         end do
         nchk(inner, outer) = nchk_local
         check_max(inner, outer) = check_local
      end do
   end do
   !$acc parallel loop gang private(max_nstep_loc) copy(max_nstep) async(async_id)
   do outer = 1, outer_size
      max_nstep_loc = 0
      !$acc loop vector private(nstep_loc, rstep, last_step) reduction(max:max_nstep_loc)
      do inner = 1, inner_size
         nstep_loc = 1
         if (check_max(inner, outer) .ge. check_point) then
            nstep_loc = int(check_max(inner, outer)/safe_step) + 1
            nstep(inner, outer) = nstep_loc
            ! GPU does not support print (need D2H copy)
            ! if (job .eq. 'advv') then
            !    print *, ' max def_cfl ', check_max(inner, outer), ' needs ', nstep(inner, outer), &
            !       'steps at level', im + 1 - nchk(inner, outer), 'of', im, 'in ', job, ' processing'
            ! else if (job .eq. 'advx' .and. im .gt. 25) then
            !    print *, ' max def_cfl ', check_max(inner, outer), ' needs ', nstep(inner, outer), &
            !       'steps in', nchk(inner, outer), 'of', im, 'at level', k, 'in ', job, ' processing'
            ! else if (job .eq. 'advy') then
            !    print *, ' max def_cfl ', check_max(inner, outer), ' needs ', nstep(inner, outer), &
            !       'steps in', nchk(inner, outer), 'of', im, 'at level', k, 'in ', job, ' processing'
            ! end if
            rstep = safe_step/check_max(inner, outer)
            !$acc loop seq
            do n = 1, nstep_loc - 1
               step(n, inner, outer) = rstep
            end do
            last_step = 1.-(nstep_loc - 1)*rstep
            step(nstep_loc, inner, outer) = last_step
         end if
         max_nstep_loc = max(max_nstep_loc, nstep_loc)
      end do
      !$acc atomic
      max_nstep = max(max_nstep, max_nstep_loc)
   end do
   !$acc wait (async_id)
   if (max_nstep .ne. 1) then
      print *, "max_nstep: ", max_nstep
      !$acc update self(check_max, nstep)
      do outer = 1, outer_size
         do inner = 1, inner_size
            if (nstep(inner, outer) .ne. 1) then
               print *, ' max def_cfl ', check_max(inner, outer), ' needs ', nstep(inner, outer), &
                  'steps of ', outer, ' at level ', inner, ' in ', job, ' processing'
            end if
         end do
      end do
   end if
   !$acc exit data delete(nchk, check_max) async(async_id)

end subroutine def_cfl_step_two_loops_gpu

subroutine locs_pp_map(pp, index, lonn, imp, sc, loc)
   !$acc routine seq
   use const, only: RTYPE

   implicit none

   real(kind=RTYPE), dimension(lonn + 1), intent(in) :: pp
   integer, intent(in) :: index, lonn, imp
   real(kind=RTYPE), intent(in) :: sc
   real(kind=RTYPE), intent(out) :: loc
   integer :: idx1, idx2

   idx1 = (index - 1)/imp
   idx2 = mod(index + imp - 1, imp) + 1
   loc = pp(idx2) + (idx1 - 1)*sc

end subroutine locs_pp_map

subroutine cyclic_cell_ppm_intp_two_loops_gpu(outer_index, outer_size, inner_size, pp, pn, qq, lonn, nvars, sc, nstep_less)
   !$acc routine(locs_pp_map) seq
   use const, only: RTYPE

   implicit none

   integer :: outer_index(3, outer_size), outer_size, inner_size, lonn, nvars
   real(kind=RTYPE), dimension(lonn + 1, inner_size, outer_size) :: pp, pn
   real(kind=RTYPE) :: qq(lonn, nvars, inner_size, outer_size)
   real(kind=RTYPE) :: sc, pn_t
   logical :: nstep_less(inner_size, outer_size)
   real(kind=RTYPE) hh(3*lonn, inner_size, outer_size)
   real(kind=RTYPE) qmi_t, qpi_t, qn_t
   real(kind=RTYPE) dqmono(3*lonn, nvars, inner_size, outer_size)
   real(kind=RTYPE) qi(3*lonn, nvars, inner_size, outer_size)
   integer kkh_array(lonn + 1, inner_size, outer_size)
   real(kind=RTYPE) tl_array(lonn + 1, inner_size, outer_size)
   real(kind=RTYPE) dql_array(lonn + 1, nvars, inner_size, outer_size)
   real(kind=RTYPE) cyclic_length
   integer ik, le, kstr(inner_size, outer_size), kend(inner_size, outer_size)
   real(kind=RTYPE) pnmin, pnmax, locbndmin, locbndmax, pnlocs, locs1, locs2, locs3
   real(kind=RTYPE) dqi, dqimax, dqimin, dqi_p, dqi_t, dqimax_p, dqimax_t, dqimin_p, dqimin_t
   real(kind=RTYPE) tl, tl2, tl3, qql, tlp, tlm, tlc
   real(kind=RTYPE) th, th2, th3, qqh, thp, thm, thc
   real(kind=RTYPE) :: dql_t, dqh_t
   real(kind=RTYPE) dpp, dqq, c1, c2, cc, r3, r6
   real(kind=RTYPE) dqq_array(lonn, nvars, inner_size, outer_size)
   real(kind=RTYPE) rt, dqmono_pre, dqmono_cur, mass_ppre, mass_pre, mass_t, mass_nxt, hh1, hh2
   real(kind=RTYPE) rdthtl
   integer i, j, k, kl, kh, kk, kkl, kkh, n, left, right, mid, kkn
   integer, parameter :: mono = 1

   integer :: outer, inner, im, imp, imf
   integer :: async_id
   logical :: has_error
   logical :: inside, accumulative

   async_id = 1
   has_error = .false.

   !$acc enter data copyin(has_error) async(async_id)
   !$acc parallel loop collapse(3) private(imp, imf, locbndmin, locbndmax, rt, pn_t) async(async_id)
   do outer = 1, outer_size
      do inner = 1, inner_size
         do i = 1, lonn + 1
            if (nstep_less(inner, outer)) then
               imp = outer_index(2, outer)
               imf = outer_index(3, outer)
               if (i .le. imf + 1) then
                  pn_t = pn(1, inner, outer)
                  call locs_pp_map(pp(1, inner, outer), imp + 4, lonn, imp, sc, locbndmin)
                  call locs_pp_map(pp(1, inner, outer), 2*imp - 4, lonn, imp, sc, locbndmax)
                  if (pn_t .lt. locbndmin - sc) then
                     rt = int((locbndmin - pn_t)/sc)*sc
                  else if (pn_t .gt. locbndmax) then
                     rt = (int((locbndmax - pn_t)/sc) - 1)*sc
                  else
                     rt = 0.0
                  end if
                  pn(i, inner, outer) = pn(i, inner, outer) + rt
               end if
            end if
         end do
      end do
   end do
   !$acc enter data create(kstr, kend) async(async_id)

   !$acc parallel loop collapse(2) private(imp, imf, pn_t, pnlocs, locs1, locs2) async(async_id)
   do outer = 1, outer_size
      do inner = 1, inner_size
         if (nstep_less(inner, outer)) then
            imp = outer_index(2, outer)
            imf = outer_index(3, outer)
            pn_t = pn(1, inner, outer) ! pn_min
            ! << kstr >>
            kstr(inner, outer) = 0
            ! call locs_pp_map(pp(1, inner, outer), imp + 1, lonn, imp, sc, pnlocs)
            pnlocs = pp(1, inner, outer)
            if (pn_t .lt. pnlocs) then
               !$acc loop vector private(locs1, locs2)
               do i = imp, 1, -1
                  call locs_pp_map(pp(1, inner, outer), i, lonn, imp, sc, locs1)
                  call locs_pp_map(pp(1, inner, outer), i + 1, lonn, imp, sc, locs2)
                  if ((pn_t .ge. locs1) .and. &
                      (pn_t .lt. locs2)) then
                     kstr(inner, outer) = i
                  end if
               end do
            else
               !$acc loop vector private(locs1, locs2)
               do i = imp + 1, 2*imp
                  call locs_pp_map(pp(1, inner, outer), i, lonn, imp, sc, locs1)
                  call locs_pp_map(pp(1, inner, outer), i + 1, lonn, imp, sc, locs2)
                  if ((pn_t .ge. locs1) .and. &
                      (pn_t .lt. locs2)) then
                     kstr(inner, outer) = i
                  end if
               end do
            end if

            if (kstr(inner, outer) .eq. 0) then
               has_error = .true.
            end if
            kstr(inner, outer) = max(3, kstr(inner, outer))

            ! << kend >>
            pn_t = pn(imf + 1, inner, outer) ! pnmax
            kend(inner, outer) = 0
            ! call locs_pp_map(pp(1, inner, outer), 2*imp + 1, lonn, imp, sc, pnlocs)
            pnlocs = pp(1, inner, outer) + sc
            if (pn_t .lt. pnlocs) then
               !$acc loop vector private(locs1, locs2)
               do i = 2*imp, imp, -1
                  call locs_pp_map(pp(1, inner, outer), i, lonn, imp, sc, locs1)
                  call locs_pp_map(pp(1, inner, outer), i + 1, lonn, imp, sc, locs2)
                  if (pn_t .ge. locs1 .and. pn_t .lt. locs2) then
                     kend(inner, outer) = i + 1
                  end if
               end do
            else
               !$acc loop vector private(locs1, locs2)
               do i = 2*imp + 1, 3*imp - 1
                  call locs_pp_map(pp(1, inner, outer), i, lonn, imp, sc, locs1)
                  call locs_pp_map(pp(1, inner, outer), i + 1, lonn, imp, sc, locs2)
                  if (pn_t .ge. locs1 .and. pn_t .lt. locs2) then
                     kend(inner, outer) = i + 1
                  end if
               end do
            end if

            if (kend(inner, outer) .eq. 0) then
               has_error = .true.
               ! call locs_pp_map(pp(1, inner, outer), imp, lonn, imp, sc, locs1)
               ! call locs_pp_map(pp(1, inner, outer), 3*imp, lonn, imp, sc, locs2)
               ! print *, ' Error: cannot get kend: pnmax locs(lonp) locs(3*lonp)', &
               !    kend(inner, outer), pnmax, locs1, locs2
               ! print *, ' Error: pn(lonn-1) pn(lonn) pn(lonn+1) ', &
               !    pn(imf - 1, inner, outer), pn(imf, inner, outer), pn(imf + 1, inner, outer)
            end if
            kend(inner, outer) = min(3*imp - 2, kend(inner, outer))
         end if
      end do
   end do

   !$acc enter data create(kkh_array) async(async_id)
   !$acc parallel loop collapse(2) private(imp, imf) async(async_id)
   do k = 1, outer_size
      do j = 1, inner_size
         if (nstep_less(j, k)) then
            imp = outer_index(2, k)
            imf = outer_index(3, k)
            kkh_array(1, j, k) = kstr(j, k)
            !$acc loop private(mid, left, right, locs1)
            do i = 1, imf
               left = kstr(j, k)
               right = kend(j, k) + 1
               do while (right - left > 1)
                  mid = (left + right)/2
                  call locs_pp_map(pp(1, j, k), mid, lonn, imp, sc, locs1)
                  if (pn(i + 1, j, k) .lt. locs1) then
                     right = mid
                  else
                     left = mid
                  end if
               end do
               kkh_array(i + 1, j, k) = left
            end do
         end if
      end do
   end do

   !$acc parallel loop collapse(2) private(imp, imf) async(async_id)
   do outer = 1, outer_size
      do inner = 1, inner_size
         if (nstep_less(inner, outer)) then
            imp = outer_index(2, outer)
            imf = outer_index(3, outer)
            !$acc loop vector private(kkl, kkh, pnmin, pnmax)
            do i = 1, imf
               kkl = kkh_array(i, inner, outer)
               kkh = kkh_array(i + 1, inner, outer)
               if (kkh .eq. kend(inner, outer) + 2) then
                  call locs_pp_map(pp(1, inner, outer), kkl + 1, lonn, imp, sc, locs1)
                  call locs_pp_map(pp(1, inner, outer), kend(inner, outer) + 1, lonn, imp, sc, locs2)
                  pnmin = pn(1, inner, outer)
                  pnmax = pn(imf + 1, inner, outer)
                  print *, ' Error in cyclic_cell_ppm_intp location not found '
                  print *, ' lons=', im, ' lonp=', imp, ' lonn=', imf
                  print *, ' pnmin=', pnmin, ' pnmax=', pnmax !
                  print *, ' pn(1)=', pn(1, inner, outer), ' pn(lonn+1)=', pn(imf + 1, inner, outer)
                  print *, ' kstr =', kstr(inner, outer), ' kend =', kend(inner, outer)
                  print *, ' kh=', i + 1, ' pn(kh)=', pn(i + 1, inner, outer)
                  print *, ' kkl +1=', kkl + 1, ' locs(kkl +1)=', locs1
                  print *, ' kend+1=', kend(inner, outer) + 1, ' locs(kend+1)=', locs2
                  has_error = .true.
               end if
               if (kkh .lt. kkl) then
                  call locs_pp_map(pp(1, inner, outer), kkl - 1, lonn, imp, sc, locs1)
                  print *, ' Error in cyclic_cell_ppm_intp location messed up '
                  print *, ' kkl=', kkl, ' kkh=', kkh
                  print *, ' kh=', i + 1, ' pn(kh)=', pn(i + 1, inner, outer)
                  print *, ' kkl-1=', kkl - 1, ' locs(kkl-1)=', locs1
                  has_error = .true.
               end if
            end do
         end if
      end do
   end do

   !$acc enter data create(hh, dqmono, qi) async(async_id)
   !$acc parallel loop gang collapse(3) async(async_id)
   do outer = 1, outer_size
      do inner = 1, inner_size
         do n = 1, nvars
            if (nstep_less(inner, outer)) then
               !$acc loop vector private(imp, mass_pre, mass_t, mass_nxt, dqi, dqimax, dqimin, locs1, locs2)
               do i = kstr(inner, outer) - 2, kend(inner, outer) + 2
                  imp = outer_index(2, outer)
                  if (n .eq. 1) then
                     call locs_pp_map(pp(1, inner, outer), i, lonn, imp, sc, locs1)
                     call locs_pp_map(pp(1, inner, outer), i + 1, lonn, imp, sc, locs2)
                     hh(i, inner, outer) = locs2 - locs1
                  end if
                  mass_pre = qq(mod(i - 2, imp) + 1, n, inner, outer)
                  mass_t = qq(mod(i - 1, imp) + 1, n, inner, outer)
                  mass_nxt = qq(mod(i, imp) + 1, n, inner, outer)
                  dqi = 0.25*(mass_nxt - mass_pre)
                  dqimax = max(mass_pre, mass_t, mass_nxt) - mass_t
                  dqimin = mass_t - min(mass_pre, mass_t, mass_nxt)
                  dqmono(i, n, inner, outer) = sign(min(abs(dqi), dqimin, dqimax), dqi)
               end do
            end if
         end do
      end do
   end do
   !$acc parallel loop gang collapse(2) private(imp) async(async_id)
   do outer = 1, outer_size
      do inner = 1, inner_size
         if (nstep_less(inner, outer)) then
            imp = outer_index(2, outer)
            !$acc loop vector private(dqmono_cur, dqmono_pre, mass_pre, mass_t, hh1, hh2, cc, rt)
            do i = kstr(inner, outer) - 1, kend(inner, outer) + 2
               call locs_pp_map(pp(1, inner, outer), i - 1, lonn, imp, sc, locs1)
               call locs_pp_map(pp(1, inner, outer), i, lonn, imp, sc, locs2)
               call locs_pp_map(pp(1, inner, outer), i + 1, lonn, imp, sc, locs3)
               hh1 = locs3 - locs2
               hh2 = locs2 - locs1
               cc = 1./(hh1 + hh2)
               hh1 = hh1*cc
               hh2 = hh2*cc
               !$acc loop
               do n = 1, nvars
                  mass_pre = qq(mod(i - 2, imp) + 1, n, inner, outer)
                  mass_t = qq(mod(i - 1, imp) + 1, n, inner, outer)
                  dqmono_pre = dqmono(i - 1, n, inner, outer)
                  dqmono_cur = dqmono(i, n, inner, outer)
                  qi(i, n, inner, outer) = mass_pre*hh1 + mass_t*hh2 + (dqmono_pre - dqmono_cur)/3.
               end do
            end do
         end if
      end do
   end do
   !$acc exit data delete(dqmono) async(async_id)
   !$acc enter data create(tl_array) async(async_id)
   !$acc enter data create(dql_array, dqq_array) async(async_id)
   !$acc parallel loop gang collapse(2) private(imp, imf) async(async_id)
   do outer = 1, outer_size
      do inner = 1, inner_size
         if (nstep_less(inner, outer)) then
            imp = outer_index(2, outer)
            imf = outer_index(3, outer)
            !$acc loop vector private(locs1, kkh, th, th2, th3, thp, thm, thc, qmi_t, qpi_t, mass_t, c1, c2, inside, dql_t, dqq, accumulative)
            do i = 1, imf + 1
               kkh = kkh_array(i, inner, outer)
               if (i .ne. imf + 1) then
                  kkn = kkh_array(i + 1, inner, outer)
               end if
               accumulative = (i .ne. imf + 1) .and. (kkh .ne. kkn)
               call locs_pp_map(pp(1, inner, outer), kkh, lonn, imp, sc, locs1)
               th = (pn(i, inner, outer) - locs1)/hh(kkh, inner, outer)
               tl_array(i, inner, outer) = th
               th2 = th*th
               th3 = th2*th
               thp = th3 - th2
               thm = th3 - 2.*th2 + th
               thc = -2.*th3 + 3.*th2
               inside = (mono .eq. 1) .and. (kkh .ge. kstr(inner, outer) - 1) .and. (kkh .le. kend(inner, outer) + 1)
               !$acc loop seq
               do n = 1, nvars
                  qmi_t = 0.0
                  qpi_t = 0.0
                  mass_t = qq(mod(kkh - 1, imp) + 1, n, inner, outer)
                  if (inside) then
                     qmi_t = qi(kkh, n, inner, outer)
                     qpi_t = qi(kkh + 1, n, inner, outer)
                     c1 = qpi_t - mass_t
                     c2 = mass_t - qmi_t
                     if (c1*c2 .le. 0.0) then
                        qmi_t = mass_t
                        qpi_t = mass_t
                     else
                        cc = qpi_t - qmi_t
                        c1 = cc*(mass_t - 0.5*(qpi_t + qmi_t))
                        c2 = cc*cc/6.
                        if (c1 .gt. c2) then
                           qmi_t = 3.*mass_t - 2.*qpi_t
                        else if (c1 .lt. -c2) then
                           qpi_t = 3.*mass_t - 2.*qmi_t
                        end if
                     end if
                  end if
                  dql_t = thp*qpi_t + thm*qmi_t + thc*mass_t
                  dql_array(i, n, inner, outer) = dql_t
                  if (accumulative) then
                     dqq = (mass_t - dql_t)*hh(kkh, inner, outer)
                     !$acc loop seq
                     do kk = kkh + 1, kkn - 1
                        dqq = dqq + qq(mod(kk - 1, imp) + 1, n, inner, outer)*hh(kk, inner, outer)
                     end do
                     dqq_array(i, n, inner, outer) = dqq
                  end if
               end do
            end do
         end if
      end do
   end do
   !$acc parallel loop gang collapse(2) private(imp, imf) async(async_id)
   do outer = 1, outer_size
      do inner = 1, inner_size
         if (nstep_less(inner, outer)) then
            imp = outer_index(2, outer)
            imf = outer_index(3, outer)
            !$acc loop vector private(kkl, kkh, dql_t, dqh_t, dqq, tl, th, dpp)
            do i = 1, imf
               dpp = 1.0
               kkl = kkh_array(i, inner, outer)
               kkh = kkh_array(i + 1, inner, outer)
               tl = tl_array(i, inner, outer)
               th = tl_array(i + 1, inner, outer)
               if (kkh .gt. kkl) then
                  dpp = (1.0 - tl)*hh(kkl, inner, outer) + th*hh(kkh, inner, outer)
                  !$acc loop seq
                  do kk = kkl + 1, kkh - 1
                     dpp = dpp + hh(kk, inner, outer)
                  end do
               end if
               !$acc loop
               do n = 1, nvars
                  dql_t = dql_array(i, n, inner, outer)
                  dqh_t = dql_array(i + 1, n, inner, outer)
                  if (kkh .eq. kkl) then
                     rdthtl = th - tl
                     if (rdthtl .ne. 0.) rdthtl = 1./rdthtl
                     qn_t = (dqh_t - dql_t)*rdthtl
                  else
                     dqq = dqq_array(i, n, inner, outer) + dqh_t*hh(kkh, inner, outer)
                     qn_t = dqq/dpp
                  end if
                  qq(i, n, inner, outer) = qn_t
               end do
            end do
         end if
      end do
   end do
   !$acc exit data delete(qi, kstr, kend) async(async_id)
   !$acc exit data delete(hh, kkh_array, tl_array, dql_array, dqq_array) copyout(has_error) async(async_id)
   !$acc wait(async_id)
   if (has_error) then
      print *, "[ERROR] There is an error in cyclic_cell_ppm_intp_two_loops_gpu."
      print *, "[ERROR] Please check pp and pn."
      call abort
   end if

end subroutine cyclic_cell_ppm_intp_two_loops_gpu
