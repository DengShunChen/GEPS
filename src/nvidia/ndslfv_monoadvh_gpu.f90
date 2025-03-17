subroutine ndslfv_monoadvh_gpu_refactor(ddtemp, qvadv, pten, vdzonl, vdmerd &
                                        , lonsperlat, deltim, xy, levs)
   ! Present on device: ddtemp, qvadv, pten, vdzonl, vdmerd, lonsperlat, ut_sl, vt_sl
   ! Present on device: jlist1, cosl, nxdef, lonlen, lonstr, latlen, jlist1_sl, gglati, fa1, fa2, fa3, fa4
   use param
   use const, only: RTYPE
   implicit none
   real(kind=RTYPE) pten(nx, levs, my_max)
   real(kind=RTYPE) ddtemp(nx, levs, my_max), qvadv(nx, levs*ncld, my_max)
   real(kind=RTYPE) vdmerd(nx, levs, my_max), vdzonl(nx, levs, my_max)
   integer, intent(in):: lonsperlat(my)
   real(kind=RTYPE), intent(in):: deltim
   integer xy, levs

   if (xy .eq. 0) call ndslfv_monoadvh2(ddtemp, qvadv, pten, vdzonl &
                                        , vdmerd, lonsperlat, deltim, levs)
   if (xy .gt. 0.5) call ndslfv_monoadvh2_xy_gpu_refactor(ddtemp, qvadv, pten, vdzonl &
                                                          , vdmerd, lonsperlat, deltim, levs)
   if (xy .lt. -0.5) call ndslfv_monoadvh2_yx_gpu_refactor(ddtemp, qvadv, pten, vdzonl &
                                                           , vdmerd, lonsperlat, deltim, levs)

end subroutine ndslfv_monoadvh_gpu_refactor

subroutine ndslfv_monoadvh2_gpu_refactor(ddtemp, qvadv, pten, vdzonl, vdmerd, &
                                         lonsperlat, deltim, levs)
!
! a routine to do non-iteration semi-Lagrangain advection
! considering advection  with monotonicity in interpolation
! contact: hann-ming henry juang
! program log
! 2011 02 20 : henry juang, created for ndsl advection
! 2013 06 20 : Henry Juang correct wind direction for north-south advection
!
   use param
   use grid
   use index
   use const
   use rank
   use cudafor
   use openacc

   implicit none

   real(kind=RTYPE) pten(nx, levs, my_max)
   real(kind=RTYPE) ddtemp(nx, levs, my_max), qvadv(nx, levs, ncld, my_max)
   real(kind=RTYPE) vdmerd(nx, levs, my_max), vdzonl(nx, levs, my_max)
   integer, intent(in):: lonsperlat(my)
   real(kind=RTYPE), intent(in):: deltim

   real(kind=RTYPE) uulon(lonfull, levs, latpart)
   real(kind=RTYPE) vvlon(lonfull, levs, latpart)
   real(kind=RTYPE) qqlon(lonfull, levs*ndslhvar, latpart)
   real(kind=RTYPE) rrlon(lonfull, levs*ndslhvar, latpart)

   real(kind=RTYPE) vvlat(latfull, levs, lonpart)
   real(kind=RTYPE) qqlat(latfull, levs*ndslhvar, lonpart)
   real(kind=RTYPE) rrlat(latfull, levs*ndslhvar, lonpart)
   real(kind=RTYPE) xr(lonfull, levs)
   real(kind=RTYPE) xcp(lonfull, levs)
   real(kind=RTYPE) sumrq(lonfull, levs)
   real(kind=RTYPE) xkappa(lonfull, levs)
   real(kind=RTYPE) rdt2, rkt, pi, cons0, cons1, rma, rm2a

   integer mono, mass, levs
   integer nlevs, nvars
   integer i, j, n, k, lon, lan, lat, lons_lat, irc, kk, KL
   integer kuu, kvv, kqq, ktt, kup, nqq
   integer ku, kv, kq, kt, kp

   integer :: async_id, istat
   integer(kind=cuda_stream_kind) :: stream

   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   mono = 1
   mass = 0
   cons0 = 0.0
   cons1 = 1.0

   nvars = ndslhvar
   nlevs = nvars*levs

   rdt2 = 0.5/deltim
   !$acc enter data create(uulon, vvlon, qqlon, vvlat, qqlat, rrlon, rrlat) async(async_id)
   !$acc host_data use_device(uulon, vvlon, qqlon, rrlon)
   istat = cudaMemsetAsync(qqlon, real(0.0, RTYPE), size(qqlon), stream)
   istat = cudaMemsetAsync(rrlon, real(0.0, RTYPE), size(qqlon), stream)
   istat = cudaMemsetAsync(uulon, real(0.0, RTYPE), size(uulon), stream)
   istat = cudaMemsetAsync(vvlon, real(0.0, RTYPE), size(vvlon), stream)
   !$acc end host_data

! =================================================================
!   prepare wind and variable in flux form with gaussina weight
! =================================================================

   !$acc parallel loop gang collapse(3) private(lat, lons_lat) async(async_id)
   do lan = 1, jlistnum
      do n = 1, ncld
         do k = 1, levs
            lat = jlist1(lan)
            lons_lat = lonsperlat(lat)
            !$acc loop vector private(rma, kk)
            do i = 1, lons_lat
               if (n .eq. 1) then
                  rma = 1./cosl(lat)
                  kk = (k - 1)*nvars
                  uulon(i, k, lan) = ut_sl(i, k, lan)*rma*rma
                  vvlon(i, k, lan) = vt_sl(i, k, lan)*rma
                  qqlon(i, kk + 1, lan) = vdzonl(i, k, lan)
                  qqlon(i, kk + 2, lan) = vdmerd(i, k, lan)
                  qqlon(i, kk + 3, lan) = ddtemp(i, k, lan)

                  rrlon(i, kk + 1, lan) = vdzonl(i, k, lan)
                  rrlon(i, kk + 2, lan) = vdmerd(i, k, lan)
                  rrlon(i, kk + 3, lan) = ddtemp(i, k, lan)
               end if
               kk = (k - 1)*nvars + 3
               qqlon(i, kk + n, lan) = qvadv(i, k, n, lan)
               rrlon(i, kk + n, lan) = qvadv(i, k, n, lan)
            end do
         end do
      end do
   end do

   call cyclic_cell_massadvx_jlist_gpu(levs, nvars, lonfull, deltim, uulon, qqlon, mass, .false.)
   call cyclic_cell_intpx_jlist_gpu(levs, nvars, lonfull, qqlon, .true.)
   call cyclic_cell_intpx_jlist_gpu(levs, nvars, lonfull, rrlon, .true.)
   call cyclic_cell_intpx_jlist_gpu(levs, 1, lonfull, vvlon, .true.)

! ---------------------------------------------------------------------
! mpi para from east-west full grid to north-south full grid
! ---------------------------------------------------------------------

   call para_we2ns_gpu(vvlon, vvlat, levs, my)
   call para_we2ns_gpu(qqlon, qqlat, nlevs, my)
   call para_we2ns_gpu(rrlon, rrlat, nlevs, my)

! ---------------------------------------------------------------------
! -------------- in north-soutn great circle -------------------
! ---------------------------------------------------------------------

   !$acc parallel loop collapse(3) private(kk) async(async_id)
   do lon = 1, mylonlen
      do k = 1, levs
         do j = 1, lathalf
            kk = (k - 1)*nvars
            qqlat(j, kk + 1, lon) = -qqlat(j, kk + 1, lon)
            qqlat(j, kk + 2, lon) = -qqlat(j, kk + 2, lon)
            rrlat(j, kk + 1, lon) = -rrlat(j, kk + 1, lon)
            rrlat(j, kk + 2, lon) = -rrlat(j, kk + 2, lon)
         end do
      end do
   end do

   call cyclic_cell_massadvy_mylonlen_gpu(latfull, levs, nvars, deltim, vvlat, qqlat, mass, .false.)
   call cyclic_cell_massadvy_mylonlen_gpu(latfull, levs, nvars, deltim, vvlat, rrlat, mass, .false.)

   !$acc parallel loop collapse(3) private(kk) async(async_id)
   do lon = 1, mylonlen
      do k = 1, levs
         do j = 1, lathalf
            kk = (k - 1)*nvars
            qqlat(j, kk + 1, lon) = -qqlat(j, kk + 1, lon)
            qqlat(j, kk + 2, lon) = -qqlat(j, kk + 2, lon)
            rrlat(j, kk + 1, lon) = -rrlat(j, kk + 1, lon)
            rrlat(j, kk + 2, lon) = -rrlat(j, kk + 2, lon)
         end do
      end do
   end do

! ----------------------------------------------------------------------
! mpi para from north-south direction to east-west direeectory
! ----------------------------------------------------------------------
   call para_ns2we_gpu(qqlat, qqlon, nlevs, my)
   call para_ns2we_gpu(rrlat, rrlon, nlevs, my)

! ---------------------------------------------------------------
! ---------------- back to east-west direction ------------------
! ---------------------------------------------------------------

   call cyclic_cell_intpx_jlist_gpu(levs, nvars, lonfull, qqlon, .false.)
   call cyclic_cell_intpx_jlist_gpu(levs, nvars, lonfull, rrlon, .false.)
   call cyclic_cell_massadvx_jlist_gpu(levs, nvars, lonfull, deltim, uulon, rrlon, mass, .false.)

   !$acc parallel loop gang collapse(3) private(lat, lons_lat) async(async_id)
   do lan = 1, jlistnum
      do n = 1, ncld
         do k = 1, levs
            lat = jlist1(lan)
            lons_lat = lonsperlat(lat)
            !$acc loop vector private(kk)
            do i = 1, lons_lat
               if (n .eq. 1) then
                  kk = (k - 1)*nvars
                  vdzonl(i, k, lan) = .5*(qqlon(i, kk + 1, lan) + rrlon(i, kk + 1, lan))
                  vdmerd(i, k, lan) = .5*(qqlon(i, kk + 2, lan) + rrlon(i, kk + 2, lan))
                  ddtemp(i, k, lan) = .5*(qqlon(i, kk + 3, lan) + rrlon(i, kk + 3, lan))
               end if
               kk = (k - 1)*nvars + 3
               qvadv(i, k, n, lan) = .5*(qqlon(i, kk + n, lan) + rrlon(i, kk + n, lan))
            end do
         end do
      end do
   end do

   !$acc exit data delete(uulon, vvlon, qqlon, vvlat, qqlat, rrlon, rrlat) async(async_id)

end subroutine ndslfv_monoadvh2_gpu_refactor

subroutine ndslfv_monoadvh2_xy_gpu_refactor(ddtemp, qvadv, pten, vdzonl, vdmerd &
                                            , lonsperlat, deltim, levs)
!
! a routine to do non-iteration semi-Lagrangain advection
! considering advection  with monotonicity in interpolation
! contact: hann-ming henry juang
! program log
! 2011 02 20 : henry juang, created for ndsl advection
! 2013 06 20 : Henry Juang correct wind direction for north-south advection
!
   use param
   use grid
   use index
   use const
   use rank
   use cudafor
   use openacc

   implicit none

   real(kind=RTYPE) pten(nx, levs, my_max)
   real(kind=RTYPE) ddtemp(nx, levs, my_max), qvadv(nx, levs, ncld, my_max)
   real(kind=RTYPE) vdmerd(nx, levs, my_max), vdzonl(nx, levs, my_max)
   integer, intent(in):: lonsperlat(my)
   real(kind=RTYPE), intent(in):: deltim

   real(kind=RTYPE) uulon(lonfull, levs, latpart)
   real(kind=RTYPE) vvlon(lonfull, levs, latpart)
   real(kind=RTYPE) qqlon(lonfull, levs*ndslhvar, latpart)

   real(kind=RTYPE) vvlat(latfull, levs, lonpart)
   real(kind=RTYPE) qqlat(latfull, levs*ndslhvar, lonpart)
   real(kind=RTYPE) xr(lonfull, levs)
   real(kind=RTYPE) xcp(lonfull, levs)
   real(kind=RTYPE) sumrq(lonfull, levs)
   real(kind=RTYPE) xkappa(lonfull, levs)
   real(kind=RTYPE) rdt2, rkt, pi, cons0, cons1, rma, rm2a

   integer mono, mass, levs
   integer nlevs, nvars
   integer i, j, n, k, lon, lan, lat, lons_lat, irc, kk, KL
   integer kuu, kvv, kqq, ktt, kup, nqq
   integer ku, kv, kq, kt, kp

   integer :: async_id, istat
   integer(kind=cuda_stream_kind) :: stream

   async_id = 1
   stream = acc_get_cuda_stream(async_id)

   mono = 1
   mass = 0
   cons0 = 0.0
   cons1 = 1.0

   nvars = ndslhvar
   nlevs = nvars*levs

   rdt2 = 0.5/deltim
   !$acc enter data create(uulon, vvlon, qqlon, vvlat, qqlat) async(async_id)
   !$acc host_data use_device(uulon, vvlon, qqlon)
   istat = cudaMemsetAsync(qqlon, real(0.0, RTYPE), size(qqlon), stream)
   istat = cudaMemsetAsync(uulon, real(0.0, RTYPE), size(uulon), stream)
   istat = cudaMemsetAsync(vvlon, real(0.0, RTYPE), size(vvlon), stream)
   !$acc end host_data

! =================================================================
!   prepare wind and variable in flux form with gaussina weight
! =================================================================

   !$acc parallel loop gang collapse(3) private(lat, lons_lat) async(async_id)
   do lan = 1, jlistnum
      do n = 1, ncld
         do k = 1, levs
            lat = jlist1(lan)
            lons_lat = lonsperlat(lat)
            !$acc loop vector private(rma, kk)
            do i = 1, lons_lat
               if (n .eq. 1) then
                  rma = 1./cosl(lat)
                  kk = (k - 1)*nvars
                  uulon(i, k, lan) = ut_sl(i, k, lan)*rma*rma
                  vvlon(i, k, lan) = vt_sl(i, k, lan)*rma
                  qqlon(i, kk + 1, lan) = vdzonl(i, k, lan)
                  qqlon(i, kk + 2, lan) = vdmerd(i, k, lan)
                  qqlon(i, kk + 3, lan) = ddtemp(i, k, lan)
               end if
               kk = (k - 1)*nvars + 3
               qqlon(i, kk + n, lan) = qvadv(i, k, n, lan)
            end do
         end do
      end do
   end do

   call cyclic_cell_massadvx_jlist_gpu(levs, nvars, lonfull, deltim, uulon, qqlon, mass, .false.)
   call cyclic_cell_intpx_jlist_gpu(levs, nvars, lonfull, qqlon, .true.)
   call cyclic_cell_intpx_jlist_gpu(levs, 1, lonfull, vvlon, .true.)

! ---------------------------------------------------------------------
! mpi para from east-west full grid to north-south full grid
! ---------------------------------------------------------------------

   call para_we2ns_gpu(vvlon, vvlat, levs, my)
   call para_we2ns_gpu(qqlon, qqlat, nlevs, my)

! ---------------------------------------------------------------------
! -------------- in north-soutn great circle -------------------
! ---------------------------------------------------------------------

   !$acc parallel loop collapse(3) private(kk) async(async_id)
   do lon = 1, mylonlen
      do k = 1, levs
         do j = 1, lathalf
            kk = (k - 1)*nvars
            qqlat(j, kk + 1, lon) = -qqlat(j, kk + 1, lon)
            qqlat(j, kk + 2, lon) = -qqlat(j, kk + 2, lon)
         end do
      end do
   end do

   call cyclic_cell_massadvy_mylonlen_gpu(latfull, levs, nvars, deltim, vvlat, qqlat, mass, .false.)

   !$acc parallel loop collapse(3) private(kk) async(async_id)
   do lon = 1, mylonlen
      do k = 1, levs
         do j = 1, lathalf
            kk = (k - 1)*nvars
            qqlat(j, kk + 1, lon) = -qqlat(j, kk + 1, lon)
            qqlat(j, kk + 2, lon) = -qqlat(j, kk + 2, lon)
         end do
      end do
   end do

! ----------------------------------------------------------------------
! mpi para from north-south direction to east-west direeectory
! ----------------------------------------------------------------------
   call para_ns2we_gpu(qqlat, qqlon, nlevs, my)

! ---------------------------------------------------------------
! ---------------- back to east-west direction ------------------
! ---------------------------------------------------------------

   call cyclic_cell_intpx_jlist_gpu(levs, nvars, lonfull, qqlon, .false.)

   !$acc parallel loop gang collapse(3) private(lat, lons_lat) async(async_id)
   do lan = 1, jlistnum
      do n = 1, ncld
         do k = 1, levs
            lat = jlist1(lan)
            lons_lat = lonsperlat(lat)
            !$acc loop vector private(kk)
            do i = 1, lons_lat
               if (n .eq. 1) then
                  kk = (k - 1)*nvars
                  vdzonl(i, k, lan) = qqlon(i, kk + 1, lan)
                  vdmerd(i, k, lan) = qqlon(i, kk + 2, lan)
                  ddtemp(i, k, lan) = qqlon(i, kk + 3, lan)
               end if
               kk = (k - 1)*nvars + 3
               qvadv(i, k, n, lan) = qqlon(i, kk + n, lan)
            end do
         end do
      end do
   end do

   !$acc exit data delete(uulon, vvlon, qqlon, vvlat, qqlat) async(async_id)

end subroutine ndslfv_monoadvh2_xy_gpu_refactor

subroutine ndslfv_monoadvh2_yx_gpu_refactor(ddtemp, qvadv, pten, vdzonl, vdmerd &
                                            , lonsperlat, deltim, levs)
!
! a routine to do non-iteration semi-Lagrangain advection
! considering advection  with monotonicity in interpolation
! contact: hann-ming henry juang
! program log
! 2011 02 20 : henry juang, created for ndsl advection
! 2013 06 20 : Henry Juang correct wind direction for north-south advection
!
   use param
   use grid
   use index
   use const
   use rank
   use openacc
   use cudafor

   implicit none

   real(kind=RTYPE) pten(nx, levs, my_max)
   real(kind=RTYPE) ddtemp(nx, levs, my_max), qvadv(nx, levs, ncld, my_max)
   real(kind=RTYPE) vdmerd(nx, levs, my_max), vdzonl(nx, levs, my_max)
   integer, intent(in):: lonsperlat(my)
   real(kind=RTYPE), intent(in):: deltim
   real(kind=RTYPE) uulon(lonfull, levs, latpart)
   real(kind=RTYPE) vvlon(lonfull, levs, latpart)
   real(kind=RTYPE) qqlon(lonfull, levs*ndslhvar, latpart)
   real(kind=RTYPE) vvlat(latfull, levs, lonpart)
   real(kind=RTYPE) qqlat(latfull, levs*ndslhvar, lonpart)
   real(kind=RTYPE) xr(lonfull, levs)
   real(kind=RTYPE) xcp(lonfull, levs)
   real(kind=RTYPE) sumrq(lonfull, levs)
   real(kind=RTYPE) xkappa(lonfull, levs)
   real(kind=RTYPE) rdt2, rkt, pi, cons0, cons1, rma, rm2a
   integer mono, mass, levs
   integer nlevs, nvars
   integer i, j, n, k, lon, lan, lat, lons_lat, irc, kk, KL
   integer kuu, kvv, kqq, ktt, kup, nqq
   integer ku, kv, kq, kt, kp
   integer async_id, istat
   integer(kind=cuda_stream_kind) stream

   async_id = 1
   stream = acc_get_cuda_stream(async_id)
   mono = 1
   mass = 0
   cons0 = 0.0
   cons1 = 1.0

   nvars = ndslhvar
   nlevs = nvars*levs
   rdt2 = 0.5/deltim
   !$acc enter data create(uulon, vvlon, qqlon, vvlat, qqlat) async(async_id)
   !$acc host_data use_device(uulon, vvlon, qqlon)
   istat = cudaMemsetAsync(qqlon, real(0.0, RTYPE), size(qqlon), stream)
   istat = cudaMemsetAsync(uulon, real(0.0, RTYPE), size(uulon), stream)
   istat = cudaMemsetAsync(vvlon, real(0.0, RTYPE), size(vvlon), stream)
   !$acc end host_data

! =================================================================
!   prepare wind and variable in flux form with gaussina weight
! =================================================================

   !$acc parallel loop collapse(3) private(lat, lons_lat, rma, rm2a, kk) async(async_id)
   do lan = 1, jlistnum
      do k = 1, levs
         do i = 1, lonfull
            lat = jlist1(lan)
            lons_lat = lonsperlat(lat)
            if (i .le. lons_lat) then
               rma = 1./cosl(lat)
               rm2a = rma*rma
               kk = (k - 1)*nvars

               uulon(i, k, lan) = ut_sl(i, k, lan)*rm2a
               vvlon(i, k, lan) = vt_sl(i, k, lan)*rma
               qqlon(i, kk + 1, lan) = vdzonl(i, k, lan)
               qqlon(i, kk + 2, lan) = vdmerd(i, k, lan)
               qqlon(i, kk + 3, lan) = ddtemp(i, k, lan)
            end if
         end do
      end do
   end do
   !$acc parallel loop collapse(4) private(lat, lons_lat, kk) async(async_id)
   do lan = 1, jlistnum
      do n = 1, ncld
         do k = 1, levs
            do i = 1, lonfull
               lat = jlist1(lan)
               lons_lat = lonsperlat(lat)
               if (i .le. lons_lat) then
                  kk = (k - 1)*nvars + 3
                  qqlon(i, kk + n, lan) = qvadv(i, k, n, lan)
               end if
            end do
         end do
      end do
   end do

   call cyclic_cell_intpx_jlist_gpu(levs, nvars, lonfull, qqlon, .true.)
   call cyclic_cell_intpx_jlist_gpu(levs, 1, lonfull, vvlon, .true.)

! ---------------------------------------------------------------------
! mpi para from east-west full grid to north-south full grid
! ---------------------------------------------------------------------

   call para_we2ns_gpu(vvlon, vvlat, levs, my)
   call para_we2ns_gpu(qqlon, qqlat, nlevs, my)

! ---------------------------------------------------------------------
! -------------- in north-soutn great circle -------------------
! ---------------------------------------------------------------------

   !$acc parallel loop collapse(3) private(kk) async(async_id)
   do lon = 1, mylonlen
      do k = 1, levs
         do j = 1, lathalf
            kk = (k - 1)*nvars
            qqlat(j, kk + 1, lon) = -qqlat(j, kk + 1, lon)
            qqlat(j, kk + 2, lon) = -qqlat(j, kk + 2, lon)
         end do
      end do
   end do
   call cyclic_cell_massadvy_mylonlen_gpu(latfull, levs, nvars, deltim, vvlat, qqlat, mass, .false.)
   !$acc parallel loop collapse(3) private(kk) async(async_id)
   do lon = 1, mylonlen
      do k = 1, levs
         do j = 1, lathalf
            kk = (k - 1)*nvars
            qqlat(j, kk + 1, lon) = -qqlat(j, kk + 1, lon)
            qqlat(j, kk + 2, lon) = -qqlat(j, kk + 2, lon)
         end do
      end do
   end do

! ----------------------------------------------------------------------
! mpi para from north-south direction to east-west direeectory
! ----------------------------------------------------------------------

   call para_ns2we_gpu(qqlat, qqlon, nlevs, my)

! ---------------------------------------------------------------
! ---------------- back to east-west direction ------------------
! ---------------------------------------------------------------

   call cyclic_cell_intpx_jlist_gpu(levs, nvars, lonfull, qqlon, .false.)
   call cyclic_cell_massadvx_jlist_gpu(levs, nvars, lonfull, deltim, uulon, qqlon, mass, .false.)
   !$acc parallel loop collapse(3) private(lat, lons_lat, kk) async(async_id)
   do lan = 1, jlistnum
      do k = 1, levs
         do i = 1, lonfull
            lat = jlist1(lan)
            lons_lat = lonsperlat(lat)
            if (i .le. lons_lat) then
               kk = (k - 1)*nvars
               vdzonl(i, k, lan) = qqlon(i, kk + 1, lan)
               vdmerd(i, k, lan) = qqlon(i, kk + 2, lan)
               ddtemp(i, k, lan) = qqlon(i, kk + 3, lan)
            end if
         end do
      end do
   end do
   !$acc parallel loop collapse(4) private(lat, lons_lat, kk) async(async_id)
   do lan = 1, jlistnum
      do n = 1, ncld
         do k = 1, levs
            do i = 1, lonfull
               lat = jlist1(lan)
               lons_lat = lonsperlat(lat)
               if (i .le. lons_lat) then
                  kk = (k - 1)*nvars + 3
                  qvadv(i, k, n, lan) = qqlon(i, kk + n, lan)
               end if
            end do
         end do
      end do
   end do
   !$acc exit data delete(uulon, vvlon, qqlon, vvlat, qqlat) async(async_id)
end subroutine ndslfv_monoadvh2_yx_gpu_refactor

subroutine ndslfv_monoadvh_fgnl_gpu_refactor(vdzonl, vdmerd, ddtemp &
                                             , lonsperlat, deltim, xy, levs, nvars, forward)
   ! Present on device: vdzonl, vdmerd, ddtemp, lonsperlat, jlist1, cosl, ut_sl, vt_sl, nxdef
   ! Present on device: lonlen, lonstr, latlen, jlist1_sl, gglati, fa1, fa2, fa3, fa4
   use param
   use const, only: RTYPE
   implicit none
   real(kind=RTYPE) vdmerd(nx, levs, my_max), vdzonl(nx, levs, my_max)
   real(kind=RTYPE) ddtemp(nx, levs, my_max)
   integer, intent(in):: lonsperlat(my)
   real(kind=RTYPE), intent(in):: deltim
   integer xy, levs, nvars
   logical forward

   if (xy .eq. 0) call ndslfv_monoadvh2_fgnl(vdzonl, vdmerd, ddtemp &
                                             , lonsperlat, deltim, levs, nvars, forward)
   if (xy .gt. 0.5) call ndslfv_monoadvh2_fgnl_xy_gpu_refactor(vdzonl, vdmerd, ddtemp &
                                                               , lonsperlat, deltim, levs, nvars, forward)
   if (xy .lt. -0.5) call ndslfv_monoadvh2_fgnl_yx_gpu_refactor(vdzonl, vdmerd, ddtemp &
                                                                , lonsperlat, deltim, levs, nvars, forward)
   return
end

subroutine ndslfv_monoadvh2_fgnl_xy_gpu_refactor(vdzonl, vdmerd, ddtemp &
                                                 , lonsperlat, deltim, levs, nvars, forward)
!
! a routine to do non-iteration semi-Lagrangain advection
! considering advection  with monotonicity in interpolation
! contact: hann-ming henry juang
! program log
! 2011 02 20 : henry juang, created for ndsl advection
! 2013 06 20 : Henry Juang correct wind direction for north-south advection
!
   use param
   use grid
   use index
   use const
   use rank

   implicit none

   integer nvars
   real(kind=RTYPE) ddtemp(nx, levs, my_max)
   real(kind=RTYPE) vdmerd(nx, levs, my_max), vdzonl(nx, levs, my_max)
   integer, intent(in):: lonsperlat(my)
   real(kind=RTYPE), intent(in):: deltim
   real(kind=RTYPE) uulon(lonfull, levs, latpart)
   real(kind=RTYPE) vvlon(lonfull, levs, latpart)
   real(kind=RTYPE) qqlon(lonfull, levs*nvars, latpart)
   real(kind=RTYPE) vvlat(latfull, levs, lonpart)
   real(kind=RTYPE) qqlat(latfull, levs*nvars, lonpart)
   real(kind=RTYPE) xr(lonfull, levs)
   real(kind=RTYPE) xcp(lonfull, levs)
   real(kind=RTYPE) sumrq(lonfull, levs)
   real(kind=RTYPE) xkappa(lonfull, levs)
   real(kind=RTYPE) rdt2, rkt, pi, cons0, cons1, rma, rm2a
   integer mono, mass, levs
   integer nlevs
   integer i, j, n, k, lon, lan, lat, lons_lat, irc, kk, KL
   integer kuu, kvv, ktt, kup, nqq
   integer ku, kv, kt, kp
   logical forward
   integer async_id

   async_id = 1
   mono = 1
   mass = 0
   cons0 = 0.0
   cons1 = 1.0
   nlevs = nvars*levs
   rdt2 = 0.5/deltim

! =================================================================
!   prepare wind and variable in flux form with gaussina weight
! =================================================================

   !$acc enter data create(uulon, vvlon, qqlon, vvlat, qqlat) async(async_id)
   !$acc parallel loop collapse(3) private(lat, lons_lat, rma, rm2a, kk, ku, kv, kt) async(async_id)
   do lan = 1, jlistnum
      do k = 1, levs
         do i = 1, lonfull
            lat = jlist1(lan)
            lons_lat = lonsperlat(lat)
            if (i .le. lons_lat) then
               rma = 1./cosl(lat)
               rm2a = rma*rma
               kk = (k - 1)*nvars
               uulon(i, k, lan) = ut_sl(i, k, lan)*rm2a
               vvlon(i, k, lan) = vt_sl(i, k, lan)*rma
               qqlon(i, kk + 1, lan) = vdzonl(i, k, lan)
               qqlon(i, kk + 2, lan) = vdmerd(i, k, lan)
               if (nvars .ge. 3) qqlon(i, kk + 3, lan) = ddtemp(i, k, lan)
            end if
         end do
      end do
   end do
   call cyclic_cell_massadvx_jlist_gpu(levs, nvars, lonfull, deltim, uulon, qqlon, mass, forward)
   call cyclic_cell_intpx_jlist_gpu(levs, nvars, lonfull, qqlon, .true.)
   call cyclic_cell_intpx_jlist_gpu(levs, 1, lonfull, vvlon, .true.)

! ---------------------------------------------------------------------
! mpi para from east-west full grid to north-south full grid
! ---------------------------------------------------------------------

   call para_we2ns_gpu(vvlon, vvlat, levs, my)
   call para_we2ns_gpu(qqlon, qqlat, nlevs, my)

! ---------------------------------------------------------------------
! -------------- in north-soutn great circle -------------------
! ---------------------------------------------------------------------

   !$acc parallel loop collapse(3) private(kk) async(async_id)
   do lon = 1, mylonlen
      do k = 1, levs
         do j = 1, lathalf
            kk = (k - 1)*nvars
            qqlat(j, kk + 1, lon) = -qqlat(j, kk + 1, lon)
            qqlat(j, kk + 2, lon) = -qqlat(j, kk + 2, lon)
         end do
      end do
   end do
   call cyclic_cell_massadvy_mylonlen_gpu(latfull, levs, nvars, deltim, vvlat, qqlat, mass, forward)
   !$acc parallel loop collapse(3) private(kk) async(async_id)
   do lon = 1, mylonlen
      do k = 1, levs
         do j = 1, lathalf
            kk = (k - 1)*nvars
            qqlat(j, kk + 1, lon) = -qqlat(j, kk + 1, lon)
            qqlat(j, kk + 2, lon) = -qqlat(j, kk + 2, lon)
         end do
      end do
   end do

! ----------------------------------------------------------------------
! mpi para from north-south direction to east-west direeectory
! ----------------------------------------------------------------------

   call para_ns2we_gpu(qqlat, qqlon, nlevs, my)
   call cyclic_cell_intpx_jlist_gpu(levs, nvars, lonfull, qqlon, .false.)
   !$acc parallel loop collapse(3) private(lat, lons_lat, kk) async(async_id)
   do lan = 1, jlistnum
      do k = 1, levs
         do i = 1, lonfull
            lat = jlist1(lan)
            lons_lat = lonsperlat(lat)
            if (i .le. lons_lat) then
               kk = (k - 1)*nvars
               vdzonl(i, k, lan) = qqlon(i, kk + 1, lan)
               vdmerd(i, k, lan) = qqlon(i, kk + 2, lan)
               if (nvars .ge. 3) ddtemp(i, k, lan) = qqlon(i, kk + 3, lan)
            end if
         end do
      end do
   end do
   !$acc exit data delete(uulon, vvlon, qqlon, vvlat, qqlat) async(async_id)
end subroutine ndslfv_monoadvh2_fgnl_xy_gpu_refactor

subroutine ndslfv_monoadvh2_fgnl_yx_gpu_refactor(vdzonl, vdmerd, ddtemp &
                                                 , lonsperlat, deltim, levs, nvars, forward)
!
! a routine to do non-iteration semi-Lagrangain advection
! considering advection  with monotonicity in interpolation
! contact: hann-ming henry juang
! program log
! 2011 02 20 : henry juang, created for ndsl advection
! 2013 06 20 : Henry Juang correct wind direction for north-south advection
!
   use param
   use grid
   use index
   use const
   use rank

   implicit none

   integer nvars
   real(kind=RTYPE) ddtemp(nx, levs, my_max)
   real(kind=RTYPE) vdmerd(nx, levs, my_max), vdzonl(nx, levs, my_max)
   integer, intent(in):: lonsperlat(my)
   real(kind=RTYPE), intent(in):: deltim
   real(kind=RTYPE) uulon(lonfull, levs, latpart)
   real(kind=RTYPE) vvlon(lonfull, levs, latpart)
   real(kind=RTYPE) qqlon(lonfull, levs*nvars, latpart)
   real(kind=RTYPE) vvlat(latfull, levs, lonpart)
   real(kind=RTYPE) qqlat(latfull, levs*nvars, lonpart)
   real(kind=RTYPE) xr(lonfull, levs)
   real(kind=RTYPE) xcp(lonfull, levs)
   real(kind=RTYPE) sumrq(lonfull, levs)
   real(kind=RTYPE) xkappa(lonfull, levs)
   real(kind=RTYPE) rdt2, rkt, pi, cons0, cons1, rma, rm2a
   integer mono, mass, levs
   integer nlevs
   integer i, j, n, k, lon, lan, lat, lons_lat, irc, kk, KL
   integer kuu, kvv, ktt, kup, nqq
   integer ku, kv, kt, kp
   logical forward
   integer async_id

   async_id = 1
   mono = 1
   mass = 0
   cons0 = 0.0
   cons1 = 1.0
   nlevs = nvars*levs
   rdt2 = 0.5/deltim

! =================================================================
!   prepare wind and variable in flux form with gaussina weight
! =================================================================
   !$acc enter data create(uulon, vvlon, qqlon, vvlat, qqlat) async(async_id)
   !$acc parallel loop collapse(3) private(lat, lons_lat, rma, rm2a, kk) async(async_id)
   do lan = 1, jlistnum
      do k = 1, levs
         do i = 1, lonfull
            lat = jlist1(lan)
            lons_lat = lonsperlat(lat)
            if (i .le. lons_lat) then
               rma = 1./cosl(lat)
               rm2a = rma/cosl(lat)
               kk = (k - 1)*nvars
               uulon(i, k, lan) = ut_sl(i, k, lan)*rm2a
               vvlon(i, k, lan) = vt_sl(i, k, lan)*rma
               qqlon(i, kk + 1, lan) = vdzonl(i, k, lan)
               qqlon(i, kk + 2, lan) = vdmerd(i, k, lan)
               if (nvars .ge. 3) qqlon(i, kk + 3, lan) = ddtemp(i, k, lan)
            end if
         end do
      end do
   end do
   call cyclic_cell_intpx_jlist_gpu(levs, nvars, lonfull, qqlon, .true.)
   call cyclic_cell_intpx_jlist_gpu(levs, 1, lonfull, vvlon, .true.)

! ---------------------------------------------------------------------
! mpi para from east-west full grid to north-south full grid
! ---------------------------------------------------------------------

   call para_we2ns_gpu(vvlon, vvlat, levs, my)
   call para_we2ns_gpu(qqlon, qqlat, nlevs, my)

! ---------------------------------------------------------------------
! -------------- in north-soutn great circle -------------------
! ---------------------------------------------------------------------

   !$acc parallel loop collapse(3) private(kk) async(async_id)
   do lon = 1, mylonlen
      do k = 1, levs
         do j = 1, lathalf
            kk = (k - 1)*nvars
            qqlat(j, kk + 1, lon) = -qqlat(j, kk + 1, lon)
            qqlat(j, kk + 2, lon) = -qqlat(j, kk + 2, lon)
         end do
      end do
   end do
   call cyclic_cell_massadvy_mylonlen_gpu(latfull, levs, nvars, deltim, vvlat, qqlat, mass, forward)
   !$acc parallel loop collapse(3) private(kk) async(async_id)
   do lon = 1, mylonlen
      do k = 1, levs
         do j = 1, lathalf
            kk = (k - 1)*nvars
            qqlat(j, kk + 1, lon) = -qqlat(j, kk + 1, lon)
            qqlat(j, kk + 2, lon) = -qqlat(j, kk + 2, lon)
         end do
      end do
   end do

! ----------------------------------------------------------------------
! mpi para from north-south direction to east-west direeectory
! ----------------------------------------------------------------------

   call para_ns2we_gpu(qqlat, qqlon, nlevs, my)

! ---------------------------------------------------------------
! ---------------- back to east-west direction ------------------
! ---------------------------------------------------------------

   call cyclic_cell_intpx_jlist_gpu(levs, nvars, lonfull, qqlon, .false.)
   call cyclic_cell_massadvx_jlist_gpu(levs, nvars, lonfull, deltim, uulon, qqlon, mass, forward)
   !$acc parallel loop collapse(3) private(lat, lons_lat, kk) async(async_id)
   do lan = 1, jlistnum
      do k = 1, levs
         do i = 1, lonfull
            lat = jlist1(lan)
            lons_lat = lonsperlat(lat)
            if (i .le. lons_lat) then
               kk = (k - 1)*nvars
               vdzonl(i, k, lan) = qqlon(i, kk + 1, lan)
               vdmerd(i, k, lan) = qqlon(i, kk + 2, lan)
               if (nvars .ge. 3) ddtemp(i, k, lan) = qqlon(i, kk + 3, lan)
            end if
         end do
      end do
   end do
   !$acc exit data delete(uulon, vvlon, qqlon, vvlat, qqlat) async(async_id)

end subroutine ndslfv_monoadvh2_fgnl_yx_gpu_refactor
