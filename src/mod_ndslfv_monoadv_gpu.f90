#define advh_gather4GPU advh_gather4GPU_nccl
#define advh_scatter4GPU advh_scatter4GPU_nccl
module mod_ndslfv_monoadv_gpu
   use const, only: RTYPE
   use cudafor
   use openacc
   implicit none
   real(kind=RTYPE), allocatable:: qq_3df(:, :, :, :), &
                                   vvlat(:, :, :)
   real(kind=RTYPE), allocatable, pinned:: umwrk(:, :), &
                                           vmwrk(:, :), &
                                           advhwrk(:, :, :)

   real(kind=RTYPE), allocatable, device:: umwrk_d(:, :), &
                                           vmwrk_d(:, :), &
                                           qpwrk_d(:, :), &
                                           ainp_d(:, :, :), &
                                           aout_d(:, :, :)

   integer, allocatable:: nydef_loc(:), imflst(:), &
                          wndmdf_all(:), wndmdf(:)
   logical:: firstcomp
   ! ----------------------------------------
   real(kind=RTYPE):: pi, two_pi, &
                      dxf, hfdxf, &
                      sc_advy
   ! ----------------------------------------
   integer:: use_gpu_num, num_ndslhvar_gpu, &
             num_batch, batch_size
   !----------------------------------------
   ! << cuda stream >>
   integer :: num_cudaST
   integer(kind=cuda_stream_kind), allocatable :: cudaST(:)
   type(cudaEvent), allocatable :: cudaEV(:)
   integer :: cudast_D2H, cudast_H2D, cudast_advh_u, cudast_advh_v
contains
   subroutine allocate_ndslfv_array_gpu
      use param
      use index, only: nsizex, jlist2_2d, nxjlen_all, nxdef
      use grid, only: nxp, ndslhvar, &
                      gglati, fa1, fa2, fa3, fa4
      use rank, only: myrank, nsize, mpi_comm_gfs
      use const, only: cosl
      use mpi
      implicit none
      real(kind=RTYPE) one

      integer npts1_size, npts2_size, npts3_size, istat
      integer i
      ! ----------------------------------------
      use_gpu_num = 1
      num_ndslhvar_gpu = ndslhvar/use_gpu_num
      if (num_ndslhvar_gpu < 1) then
         print *, "use too many gpus"
         use_gpu_num = 6
         batch_size = 1
      end if

      num_batch = 1
      batch_size = num_ndslhvar_gpu/num_batch

      allocate (cudaST(-1:num_batch + 4), cudaEV(4))

      do i = 1, num_batch + 4
         istat = cudaStreamCreatewithFlags(cudaST(i), cudastreamnonblocking)
         call acc_set_cuda_stream(i, cudaST(i))
      end do
      ! <<<< -----------------------------------
      ! << for CUDA-AWARE MPI >>
      cudaST(-1) = 0
      call acc_set_cuda_stream(-1, cudaST(-1))
      cudast_d2h = -1
      cudast_advh_u = -1
      cudast_advh_v = -1
      ! ++++++++++++++++++++++++++++++++++++++++
      ! << for batch >>
      ! cudast_d2h    = num_batch+1
      ! cudast_advh_u = num_batch+2
      ! cudast_advh_v = num_batch+3
      ! ----------------------------------- >>>>
      do i = 1, 4
         istat = cudaEventCreate(cudaEV(i))
      end do
      ! ----------------------------------------
      allocate (nydef_loc(nx), imflst(my), &
                wndmdf_all(ndslhvar), &
                wndmdf(num_ndslhvar_gpu))
      nydef_loc = my*2
      imflst = nx
      wndmdf_all = 1
      wndmdf_all(1:2) = (/-1, -1/)

      do i = 1, ndslhvar
         if (myrank .eq. mod(i - 1, use_gpu_num)) then
            wndmdf(mod(i - 1, num_ndslhvar_gpu) + 1) = wndmdf_all(i)
         end if
      end do
      ! ----------------------------------------
      one = 1.0
      npts1_size = nx*my*lev*sizeof(one)
      npts2_size = nxp*my_max*nsize*lev*sizeof(one)
      ! if (myrank .eq. 0) then
      !    print *, "num_gpu:", use_gpu_num
      !    print *, "num_batch: ", num_batch
      !    print *, "batch_size: ", batch_size
      !    print *, "ncld:", ncld
      !    print *, "ndslhvar:", ndslhvar
      ! end if

      allocate ( &
         umwrk(nxp*lev, my_max*nsize), &
         vmwrk(nxp*lev, my_max*nsize), &
         advhwrk(nxp*lev, my_max*nsize, num_ndslhvar_gpu), &
         stat=istat)
      if (istat /= 0) then
         write (6, *) 'mod_ndsl_monoadv_gpu : allocate fail 1'
         stop
      end if

      if (myrank < use_gpu_num) then
         allocate (vvlat(my*2, nx/2, lev), &
                   qq_3df(nx, my, lev, batch_size), &
                   stat=istat)

         allocate ( &
            umwrk_d(nxp*lev, my_max*nsize), &
            vmwrk_d(nxp*lev, my_max*nsize), &
            qpwrk_d(nxp*lev*ncld, my_max*nsize), &
            ainp_d(nxp*lev, my_max*nsize, batch_size), &
            aout_d(nxp*lev, my_max*nsize, batch_size), &
            stat=istat)

         if (istat /= 0) then
            write (6, *) 'mod_ndsl_monoadv_gpu : allocate fail 3'
            stop
         end if
         call acc_map_data(umwrk, umwrk_d, npts2_size)
         call acc_map_data(vmwrk, vmwrk_d, npts2_size)

         pi = 4.0*atan(1.0)
         two_pi = 2.0*pi
         dxf = two_pi/nx
         hfdxf = 0.5*dxf
         sc_advy = gglati(my*2 + 1) - gglati(1)

         !! !$acc enter data &
         !! !$acc& copyin(&
         !! !$acc& nx,my,lev,ncld,nxp,nsizex,batch_size, &
         !! !$acc& gglati(:my*2+1), fa1(:my*2), fa2(:my*2), &
         !! !$acc& fa3(:my*2), fa4(:my*2), &
         !! !$acc& cosl(:my), &
         !! !$acc& jlist2_2d(:nsizex,:nx), nxjlen_all(:nsizex,:nx), &
         !! !$acc& nxdef(:my), imflst(:my), &
         !! !$acc& nydef_loc(:nx), &
         !! !$acc& pi, two_pi, dxf, hfdxf, sc_advy) &
         !! !$acc& create(vvlat(:my*2,:nx/2,:lev))

         !$acc enter data &
         !$acc& copyin(&
         !$acc& batch_size, &
         !$acc& imflst(:my), &
         !$acc& nydef_loc(:nx), &
         !$acc& pi, two_pi, dxf, hfdxf, sc_advy) &
         !$acc& create(vvlat(:my*2,:nx/2,:lev))

      end if

   end subroutine allocate_ndslfv_array_gpu
   ! ============================================================
   subroutine deallocate_ndslfv_array_gpu
      use rank, only: myrank
      implicit none

      if (myrank < use_gpu_num) then
         call acc_unmap_data(umwrk)
         call acc_unmap_data(vmwrk)

         !$acc exit data &
         !$acc& delete( &
         !$acc& batch_size, &
         !$acc& imflst, &
         !$acc& nydef_loc, &
         !$acc& pi, two_pi, dxf, hfdxf, sc_advy, &
         !$acc& vvlat &
         !$acc& )

         deallocate (umwrk_d, vmwrk_d, qpwrk_d, ainp_d, aout_d)
         deallocate (vvlat, qq_3df)
      end if

      deallocate (cudaST, cudaEV)
      deallocate (nydef_loc, imflst, wndmdf_all, wndmdf)
      deallocate (umwrk, vmwrk, advhwrk)

   end subroutine deallocate_ndslfv_array_gpu
   ! ============================================================
   subroutine ndslfv_monoadvh_gpu(tt, pt, ut, vt, qt, qm, &
                                  um, vm, deltim, xy, forward)
      use param
      use rank
      use index
      use const, only: RTYPE, cosl, MPI_RTYPE
      use grid, only: ndslhvar
      use mpi
      use cudafor
      implicit none
      ! input
      real(kind=RTYPE):: um(nxp, lev, my_max), &
                         vm(nxp, lev, my_max), &
                         qm(nxp, lev*ncld, my_max)

      real(kind=RTYPE), intent(in):: deltim
      logical, intent(in):: forward
      ! output
      real(kind=RTYPE):: qt(nxp, lev*ncld, my_max), &
                         pt(nxp, lev, my_max), &
                         ut(nxp, lev, my_max), &
                         vt(nxp, lev, my_max), &
                         tt(nxp, lev, my_max)
      ! local
      real(kind=RTYPE):: um3df(nx, my, lev), &
                         vm3df(nx, my, lev), &
                         rma, rm2a
      ! integer wndmdf(ndslhvar)
      !
      integer xy, istat, i, j, k, n
      integer ierr, ii, jf, nn, kk, m, pts
      ! wndmdf = 1
      ! wndmdf(1:2) = (/-1, -1/)
      firstcomp = .true.
      ! ----------------------------------------
      !$acc host_data use_device(um, vm, ut, vt, tt, qm)
      call advh_gather4GPU(umwrk_d, um, 1, 0)
      call advh_gather4GPU(vmwrk_d, vm, 1, 0)
      ! u v t at n-1
      call advh_gather4GPU(ainp_d(1, 1, 1), ut, 1, 0)
      call advh_gather4GPU(ainp_d(1, 1, 2), vt, 1, 0)
      call advh_gather4GPU(ainp_d(1, 1, 3), tt, 1, 0)
      ! rq at n-1
      call advh_gather4GPU(qpwrk_d, qm, ncld, 0)
      !$acc end host_data
      ! ----------------------------------------
      if (myrank .eq. 0) then
         !$acc data &
         !$acc& create( &
         !$acc& qq_3df(:nx,:my,:lev,:ndslhvar), &
         !$acc& um3df(:nx,:my,:lev), &
         !$acc& vm3df(:nx,:my,:lev) &
         !$acc& ) &
         !$acc& present(jlist2_2d(nsizex,my), nxjlen_all(nsizex,my), &
         !$acc& cosl(:my) &
         !$acc& ) &
         !$acc& copyin(deltim)
         ! --------------------

         !$acc parallel loop collapse(2) &
         !$acc& deviceptr(ainp_d) &
         !$acc& private(ii,jf,nn,kk,i,m,rm2a,rma)
         do k = 1, lev
            do j = 1, my
               ii = 0
               !$acc loop seq
               do i = 1, nsizex
                  jf = jlist2_2d(i, j)
                  nn = nxjlen_all(i, j)
                  kk = (k - 1)*nxp
                  !$acc loop independent
                  do m = 1, nn
                     ! wind at time step n
                     qq_3df(ii + m, j, k, 1) = ainp_d(kk + m, jf, 1)
                     qq_3df(ii + m, j, k, 2) = ainp_d(kk + m, jf, 2)
                     qq_3df(ii + m, j, k, 3) = ainp_d(kk + m, jf, 3)
                  end do
                  ii = ii + nn
               end do
            end do
         end do

         !$acc parallel loop collapse(3) deviceptr(qpwrk_d)&
         !$acc& private(ii,jf,nn,kk,i,m)
         do n = 1, ncld
            do k = 1, lev
               do j = 1, my
                  ii = 0
                  !$acc loop seq
                  do i = 1, nsizex
                     jf = jlist2_2d(i, j)
                     nn = nxjlen_all(i, j)
                     kk = (k - 1)*nxp + (n - 1)*nxp*lev
                     !$acc loop independent
                     do m = 1, nn
                        qq_3df(ii + m, j, k, 3 + n) = qpwrk_d(kk + m, jf)
                     end do
                     ii = ii + nn
                  end do
               end do
            end do
         end do

         if (xy .gt. 0.5) then
            call ndslfv_monoadvh2_xy_gpu(qq_3df, deltim, um3df, vm3df, &
                                         wndmdf, forward, lev, ndslhvar, -1)
         elseif (xy .lt. -0.5) then
            call ndslfv_monoadvh2_yx_gpu(qq_3df, deltim, um3df, vm3df, &
                                         wndmdf, forward, lev, ndslhvar, -1)
         end if
         ! ========================================
         !$acc parallel loop collapse(2) &
         !$acc& deviceptr(aout_d) &
         !$acc& private(ii,jf,nn,kk,i,m,rm2a,rma)
         do k = 1, lev
            do j = 1, my
               ii = 0
               !$acc loop seq
               do i = 1, nsizex
                  jf = jlist2_2d(i, j)
                  nn = nxjlen_all(i, j)
                  kk = (k - 1)*nxp
                  !$acc loop independent
                  do m = 1, nn
                     aout_d(kk + m, jf, 1) = qq_3df(ii + m, j, k, 1)
                     aout_d(kk + m, jf, 2) = qq_3df(ii + m, j, k, 2)
                     aout_d(kk + m, jf, 3) = qq_3df(ii + m, j, k, 3)
                  end do
                  ii = ii + nn
               end do
            end do
         end do

         !$acc parallel loop collapse(3) deviceptr(qpwrk_d)&
         !$acc& private(ii,jf,nn,kk,i,m)
         do n = 1, ncld
            do k = 1, lev
               do j = 1, my
                  ii = 0
                  !$acc loop seq
                  do i = 1, nsizex
                     jf = jlist2_2d(i, j)
                     nn = nxjlen_all(i, j)
                     kk = (k - 1)*nxp + (n - 1)*nxp*lev
                     !$acc loop independent
                     do m = 1, nn
                        qpwrk_d(kk + m, jf) = qq_3df(ii + m, j, k, 3 + n)
                     end do
                     ii = ii + nn
                  end do
               end do
            end do
         end do
         !$acc end data
      end if
      ! call mpe_barrier
      ! ----------------------------------------
      ! u v t at n
      !$acc host_data use_device(ut, vt, tt, qt)
      call advh_scatter4GPU(ut, aout_d(1, 1, 1), 1, 0)
      call advh_scatter4GPU(vt, aout_d(1, 1, 2), 1, 0)
      call advh_scatter4GPU(tt, aout_d(1, 1, 3), 1, 0)
      ! rq at n
      call advh_scatter4GPU(qt, qpwrk_d, ncld, 0)
      !$acc end host_data
      ! ----------------------------------------

      return
   end subroutine ndslfv_monoadvh_gpu
   ! ============================================================
   subroutine ndslfv_monoadvh_fgnl_gpu(vdzonl, vdmerd, ddtemp, &
                                       ut, vt, tt, um, vm, deltim, &
                                       xy, nvars, forward)
      use param
      use rank
      use index
      use const, only: RTYPE, cosl, MPI_RTYPE
      use grid, only: ndslhvar
      use mpi
      use cudafor
      implicit none
      ! input
      real(kind=RTYPE):: ut(nxp, lev, my_max), &
                         vt(nxp, lev, my_max), &
                         tt(nxp, lev, my_max), &
                         um(nxp, lev, my_max), &
                         vm(nxp, lev, my_max)
      real(kind=RTYPE), intent(in):: deltim
      logical, intent(in):: forward
      integer, intent(in):: nvars
      ! output
      real(kind=RTYPE), intent(out):: vdzonl(nxp, lev, my_max), &
                                      vdmerd(nxp, lev, my_max), &
                                      ddtemp(nxp, lev, my_max)
      ! local
      real(kind=RTYPE):: um3df(nx, my, lev), &
                         vm3df(nx, my, lev), &
                         rma, rm2a
      ! integer wndmdf(ndslhvar)
      !
      integer xy, istat, i, j, k, n
      integer ierr, ii, jf, nn, kk, m, pts
      ! wndmdf = 1
      ! wndmdf(1:2) = (/-1, -1/)
      firstcomp = .true.
      ! ----------------------------------------
      !$acc host_data use_device(um, vm, ut, vt, tt)
      call advh_gather4GPU(umwrk_d, um, 1, 0)
      call advh_gather4GPU(vmwrk_d, vm, 1, 0)
      ! u v t at n-1
      call advh_gather4GPU(ainp_d(1, 1, 1), ut, 1, 0)
      call advh_gather4GPU(ainp_d(1, 1, 2), vt, 1, 0)
      call advh_gather4GPU(ainp_d(1, 1, 3), tt, 1, 0)
      !$acc end host_data
      ! ----------------------------------------
      if (myrank .eq. 0) then
         !$acc data &
         !$acc& create( &
         !$acc& qq_3df(:nx,:my,:lev, :ndslhvar), &
         !$acc& um3df(:nx,:my,:lev), &
         !$acc& vm3df(:nx,:my,:lev) &
         !$acc& ) &
         !$acc& present(jlist2_2d(nsizex,my), nxjlen_all(nsizex,my), &
         !$acc& cosl(:my) &
         !$acc& ) &
         !$acc& copyin(deltim)
         ! --------------------

         !$acc parallel loop collapse(2) &
         !$acc& deviceptr(ainp_d) &
         !$acc& private(ii,jf,nn,kk,i,m,rm2a,rma)
         do k = 1, lev
            do j = 1, my
               ii = 0
               !$acc loop seq
               do i = 1, nsizex
                  jf = jlist2_2d(i, j)
                  nn = nxjlen_all(i, j)
                  kk = (k - 1)*nxp
                  !$acc loop independent
                  do m = 1, nn
                     ! wind at time step n
                     qq_3df(ii + m, j, k, 1) = ainp_d(kk + m, jf, 1)
                     qq_3df(ii + m, j, k, 2) = ainp_d(kk + m, jf, 2)
                     qq_3df(ii + m, j, k, 3) = ainp_d(kk + m, jf, 3)
                  end do
                  ii = ii + nn
               end do
            end do
         end do

         if (xy .gt. 0.5) then
            call ndslfv_monoadvh2_xy_gpu(qq_3df, deltim, um3df, vm3df, &
                                         wndmdf, forward, lev, 3, -1)
         elseif (xy .lt. -0.5) then
            call ndslfv_monoadvh2_yx_gpu(qq_3df, deltim, um3df, vm3df, &
                                         wndmdf, forward, lev, 3, -1)
         end if
         ! ========================================
         !$acc parallel loop collapse(2) &
         !$acc& deviceptr(aout_d) &
         !$acc& private(ii,jf,nn,kk,i,m,rm2a,rma)
         do k = 1, lev
            do j = 1, my
               ii = 0
               !$acc loop seq
               do i = 1, nsizex
                  jf = jlist2_2d(i, j)
                  nn = nxjlen_all(i, j)
                  kk = (k - 1)*nxp
                  !$acc loop independent
                  do m = 1, nn
                     aout_d(kk + m, jf, 1) = qq_3df(ii + m, j, k, 1)
                     aout_d(kk + m, jf, 2) = qq_3df(ii + m, j, k, 2)
                     aout_d(kk + m, jf, 3) = qq_3df(ii + m, j, k, 3)
                  end do
                  ii = ii + nn
               end do
            end do
         end do
         !$acc end data
      end if
      ! call mpe_barrier
      ! ----------------------------------------
      ! u v t at n
      !$acc host_data use_device(vdzonl, vdmerd, ddtemp)
      call advh_scatter4GPU(vdzonl, aout_d(1, 1, 1), 1, 0)
      call advh_scatter4GPU(vdmerd, aout_d(1, 1, 2), 1, 0)
      call advh_scatter4GPU(ddtemp, aout_d(1, 1, 3), 1, 0)
      !$acc end host_data
      ! ----------------------------------------

      return
   end subroutine ndslfv_monoadvh_fgnl_gpu
   ! ============================================================
   subroutine ndslfv_monoadvv_gpu(ddtemp, qvadv, vdzonl, vdmerd, pdot, &
                                  pt, lonsperlat, deltim, forward)
      !
      ! a routine to do non-iteration semi-Lagrangain advection
      ! considering advection  with monotonicity in interpolation
      ! contact: hann-ming henry juang
      ! program log:
      ! 2011 02 20 : henry juang, initial implemented into nems as NDSL with mass_dp
      ! 2013 09 30 : henry juang, add option of theta advection, (used later)
      !
      !
      use openacc
      use cudafor
      use param
      use grid, only: latpart, ndslvvar, lonfull
      use index, only: jlistnum, jlist1, nxp, nxjp_acc, nxjp, nxptot
      use rank
      use const
      use mpe
      implicit none

      real(kind=RTYPE), intent(inout):: ddtemp(nxp, lev, my_max), &
                                        qvadv(nxp, lev, ncld, my_max), &
                                        vdmerd(nxp, lev, my_max), &
                                        vdzonl(nxp, lev, my_max)
      real(kind=RTYPE), intent(in):: pdot(nxp, lev + 1, latpart), &
                                     pt(nxp, latpart)
      integer, intent(in):: lonsperlat(my)
      real(kind=RTYPE), intent(in):: deltim

      real(kind=RTYPE):: plev(lev + 1, nxptot), &
                         qqlon(lev, ndslvvar, nxptot)

      integer mono, mass
      integer ii, i, n, k, kk, lon, lat, lons_lat, istr, j, nxj
      logical forward
      integer istat
      integer, parameter :: async_id = 1
      integer(kind=cuda_stream_kind) :: stream
      !$acc data create(plev, qqlon) async(async_id)

      mono = 1
      mass = 0
      !
      stream = acc_get_cuda_stream(async_id)
      !$acc host_data use_device(plev)
      istat = cudaMemsetAsync(plev, real(0.0, RTYPE), size(plev), stream)
      !$acc end host_data

      !$acc parallel loop async(async_id) &
      !$acc& private(lat,lons_lat,istr)
      do j = 1, jlistnum
         lat = jlist1(j)
         lons_lat = lonsperlat(lat)
         istr = nxjp_acc(j) - 1
         !$acc loop seq private(kk)
         do k = lev, 1, -1
            kk = lev - k + 1
            !$acc loop vector
            do i = 1, lons_lat
               plev(k, istr + i) = plev(k + 1, istr + i) &
                                   + dsigma(kk, 1)*pt(i, j) + dsigma(kk, 2)
            end do
         end do
      end do

      ! d z t at n+1*
      !$acc parallel loop collapse(3) async(async_id) &
      !$acc& private(lat,lons_lat,istr,kk)
      do j = 1, jlistnum
         do k = 1, lev
            do i = 1, nxp
               lat = jlist1(j)
               lons_lat = lonsperlat(lat)
               if (i .le. lons_lat) then
                  istr = nxjp_acc(j) - 1
                  kk = lev - k + 1
                  qqlon(kk, 1, istr + i) = vdzonl(i, k, j)
                  qqlon(kk, 2, istr + i) = vdmerd(i, k, j)
                  qqlon(kk, 3, istr + i) = ddtemp(i, k, j)
               end if
            end do
         end do
      end do

      ! rq at n+1*
      !$acc parallel loop collapse(4) async(async_id) &
      !$acc& private(lat,lons_lat,istr,kk)
      do j = 1, jlistnum
         do n = 1, ncld
            do k = 1, lev
               do i = 1, nxp
                  lat = jlist1(j)
                  lons_lat = lonsperlat(lat)
                  if (i .le. lons_lat) then
                     istr = nxjp_acc(j) - 1
                     kk = lev - k + 1
                     qqlon(kk, n + 3, istr + i) = qvadv(i, k, n, j)
                  end if
               end do
            end do
         end do
      end do
      !$acc end parallel loop

      call vertical_cell_advect_gpu(lons_lat, nxptot, lev, ndslvvar, &
                                    deltim, plev, pdot, &
                                    qqlon, mass, forward, async_id)

      ! u v t tendency at n
      !$acc parallel loop collapse(3) independent async(async_id)&
      !$acc& private(lat,lons_lat,istr,kk)
      do j = 1, jlistnum
         do i = 1, nxp
            do k = 1, lev
               lat = jlist1(j)
               lons_lat = lonsperlat(lat)
               if (i .le. lons_lat) then
                  istr = nxjp_acc(j) - 1
                  kk = lev - k + 1
                  vdzonl(i, kk, j) = qqlon(k, 1, istr + i)
                  vdmerd(i, kk, j) = qqlon(k, 2, istr + i)
                  ddtemp(i, kk, j) = qqlon(k, 3, istr + i)
               end if
            end do
         end do
      end do

      ! rq tendency at n
      !$acc parallel loop collapse(4) independent async(async_id)&
      !$acc& private(lat,lons_lat,istr,kk)
      do j = 1, jlistnum
         do i = 1, nxp
            do n = 1, ncld
               do k = 1, lev
                  lat = jlist1(j)
                  lons_lat = lonsperlat(lat)
                  if (i .le. lons_lat) then
                     istr = nxjp_acc(j) - 1
                     kk = lev - k + 1
                     qvadv(i, kk, n, j) = qqlon(k, n + 3, istr + i)
                  end if
               end do
            end do
         end do
      end do
      !$acc end data

      return
   end subroutine ndslfv_monoadvv_gpu
   ! ============================================================
   subroutine ndslfv_monoadvv_fgnl_gpu(vdzonl, vdmerd, ddtemp, pdot, &
                                       pt, lonsperlat, deltim, nvars, forward)
      !
      ! a routine to do non-iteration semi-Lagrangain advection
      ! considering advection  with monotonicity in interpolation
      ! contact: hann-ming henry juang
      ! program log:
      ! 2011 02 20 : henry juang, initial implemented into nems as NDSL with mass_dp
      ! 2013 09 30 : henry juang, add option of theta advection, (used later)
      !
      !
      use openacc
      use cudafor
      use param
      use grid, only: latpart, ndslvvar
      use index, only: jlistnum, jlist1, nxp, nxjp_acc, nxjp, nxptot
      use rank
      use const
      use nvtx
      implicit none

      real(kind=RTYPE), intent(inout):: ddtemp(nxp, lev, my_max), &
                                        vdmerd(nxp, lev, my_max), &
                                        vdzonl(nxp, lev, my_max)
      real(kind=RTYPE), intent(in):: pdot(nxp, lev + 1, latpart), &
                                     pt(nxp, latpart)
      integer, intent(in):: lonsperlat(my), nvars
      real(kind=RTYPE), intent(in):: deltim

      real(kind=RTYPE):: plev(lev + 1, nxptot), &
                         qqlon(lev, nvars, nxptot)

      integer mono, mass
      integer ii, i, n, k, kk, lon, lat, lons_lat, istr, j, nxj
      logical forward
      integer istat
      integer, parameter :: async_id = 1
      integer(kind=cuda_stream_kind) :: stream

      !$acc data async(async_id) &
      !$acc& create(plev, qqlon)
      mono = 1
      mass = 0
      !
      stream = acc_get_cuda_stream(async_id)
      !$acc host_data use_device(plev)
      istat = cudaMemsetAsync(plev, real(0.0, RTYPE), size(plev), stream)
      !$acc end host_data

      !$acc parallel loop async(async_id) &
      !$acc& private(lat,lons_lat,istr)
      do j = 1, jlistnum
         lat = jlist1(j)
         lons_lat = lonsperlat(lat)
         istr = nxjp_acc(j) - 1
         !$acc loop seq private(kk)
         do k = lev, 1, -1
            kk = lev - k + 1
            !$acc loop vector
            do i = 1, lons_lat
               plev(k, istr + i) = plev(k + 1, istr + i) &
                                   + dsigma(kk, 1)*pt(i, j) + dsigma(kk, 2)
            end do
         end do
         ! d z t at n+1*
         !$acc loop vector collapse(2) independent &
         !$acc& private(kk)
         do k = 1, lev
            do i = 1, lons_lat
               kk = lev - k + 1
               qqlon(k, 1, istr + i) = vdzonl(i, kk, j)
               qqlon(k, 2, istr + i) = vdmerd(i, kk, j)
               qqlon(k, 3, istr + i) = ddtemp(i, kk, j)
            end do
         end do
      end do
      !$acc end parallel loop

      call vertical_cell_advect_gpu(lons_lat, nxptot, lev, nvars, &
                                    deltim, plev, pdot, &
                                    qqlon, mass, forward, async_id)

      !$acc parallel loop independent async(async_id)&
      !$acc& private(lat,lons_lat,istr)
      do j = 1, jlistnum
         lat = jlist1(j)
         lons_lat = lonsperlat(lat)
         istr = nxjp_acc(j) - 1

         ! u v t tendency at n
         !$acc loop vector collapse(2) independent &
         !$acc& private(kk)
         do i = 1, lons_lat
            do k = 1, lev
               kk = lev - k + 1
               vdzonl(i, kk, j) = qqlon(k, 1, istr + i)
               vdmerd(i, kk, j) = qqlon(k, 2, istr + i)
               ddtemp(i, kk, j) = qqlon(k, 3, istr + i)
            end do
         end do
      end do
      !$acc end parallel loop
      !$acc end data

      return
   end subroutine ndslfv_monoadvv_fgnl_gpu
   ! ============================================================
   subroutine ndslfv_update_gpu(lonsperlat, vdzonl, vdmerd, &
                                vdzonlr, vdmerdr, deltim, forward)
!  Present on device: lonsperlat, vdzonl, vdmerd, vdzonlr, vdmerdr, jlist1
!  update all horizontal components into momentum eqs
!  for Semi-Lagrangian vertical advection

      use index
      use param, only: nx, my, lev, my_max
      use const, only: RTYPE
      implicit none
      integer, intent(in):: lonsperlat(my)
      real(kind=RTYPE), intent(in):: deltim
      real(kind=RTYPE) vdmerd(nxp, lev, my_max), vdzonl(nxp, lev, my_max)
      real(kind=RTYPE) vdmerdr(nxp, lev, my_max), vdzonlr(nxp, lev, my_max)
      integer i, k, j, lat, lons_lat
      integer dt2
      logical forward
      integer, parameter :: async_id = 1
      if (forward) then
         dt2 = deltim
      else
         dt2 = 2.*deltim
      end if
      !$acc parallel loop collapse(3) async(async_id) &
      !$acc& present(jlist1, lonsperlat, vdzonl, vdzonlr, vdmerd, vdmerdr) &
      !$acc& private(lat, lons_lat)
      do j = 1, jlistnum
      do k = 1, lev
      do i = 1, nxp
         lat = jlist1(j)
         lons_lat = lonsperlat(lat)
         if (i .le. lons_lat) then
            vdzonl(i, k, j) = vdzonlr(i, k, j)*dt2 + vdzonl(i, k, j)
            vdmerd(i, k, j) = vdmerdr(i, k, j)*dt2 + vdmerd(i, k, j)
         end if
      end do
      end do
      end do
      !$acc end parallel loop
      return
   end subroutine ndslfv_update_gpu
   ! ============================================================
   subroutine ndslfv_monoadvh_gpu_batch(ttm_sl, pten_sl, uum_sl, vvm_sl, qm_sl, &
                                        um, vm, deltim, xy, forward)
      use param
      use rank
      use index
      use const, only: RTYPE, cosl, MPI_RTYPE
      use grid, only: ut, vt, tt, qm, ut_sl, vt_sl, ndslhvar
      use mpi
      use cudafor
      implicit none

      ! input
      real(kind=RTYPE) um(nxp, lev, my_max), vm(nxp, lev, my_max)
      real(kind=RTYPE), intent(in):: deltim
      logical, intent(in):: forward
      ! output
      real(kind=RTYPE), intent(out):: qm_sl(nxp, lev*ncld, my_max), &
                                      pten_sl(nxp, lev, my_max), &
                                      uum_sl(nxp, lev, my_max), &
                                      vvm_sl(nxp, lev, my_max), &
                                      ttm_sl(nxp, lev, my_max)

      ! local
      real(kind=RTYPE):: um3df(nx, my, lev), &
                         vm3df(nx, my, lev), &
                         rma, rm2a
      !
      real(kind=RTYPE) uttmp(nxp*lev, my_max*nsize)
      real(kind=RTYPE) qm_sl_wrk(nxp, lev*ncld, my_max)
      real maxerr, locerr
      integer lonp, jj, nxj, l
      !
      integer xy, istat, i, j, k, n
      integer ierr, ii, jf, nn, kk, m, pts
      integer stlev, lb

      firstcomp = .true.
      ! ----------------------------------------
      do i = 0, use_gpu_num - 1
         call advh_gather4GPU(umwrk, um, 1, i)
         call advh_gather4GPU(vmwrk, vm, 1, i)
      end do

      ! u v t at n-1
      call advh_gather4GPU_q(advhwrk(1, 1, mod(0, num_ndslhvar_gpu) + 1), ut, &
                             1, 1, mod(0, use_gpu_num))
      call advh_gather4GPU_q(advhwrk(1, 1, mod(1, num_ndslhvar_gpu) + 1), vt, &
                             1, 1, mod(1, use_gpu_num))
      call advh_gather4GPU_q(advhwrk(1, 1, mod(2, num_ndslhvar_gpu) + 1), tt, &
                             1, 1, mod(2, use_gpu_num))
      ! rq at n-1
      call advh_gather4GPU_q(advhwrk(1, 1, mod(3, num_ndslhvar_gpu) + 1), qm, &
                             1, 3, mod(3, use_gpu_num))
      call advh_gather4GPU_q(advhwrk(1, 1, mod(4, num_ndslhvar_gpu) + 1), qm, &
                             2, 3, mod(4, use_gpu_num))
      call advh_gather4GPU_q(advhwrk(1, 1, mod(5, num_ndslhvar_gpu) + 1), qm, &
                             3, 3, mod(5, use_gpu_num))
      ! ----------------------------------------
      if (myrank < use_gpu_num) then
         pts = nxp*lev*my_max*nsize
         istat = cudamemcpyasync(vmwrk_d, vmwrk, pts, &
                                 cudamemcpyhosttodevice, &
                                 cudast(cudast_advh_v))
         istat = cudamemcpyasync(umwrk_d, umwrk, pts, &
                                 cudamemcpyhosttodevice, &
                                 cudast(cudast_advh_u))
         !$acc data &
         !$acc& create( &
         !$acc& qq_3df(:nx,:my,:lev,:batch_size), &
         !$acc& um3df(:nx,:my,:lev), &
         !$acc& vm3df(:nx,:my,:lev) &
         !$acc& ) &
         !$acc& copyin(deltim)
         ! --------------------
         if (xy .gt. 0.5) then
            do n = 1, num_ndslhvar_gpu, batch_size
               istat = cudamemcpyasync(ainp_d, advhwrk(1, 1, n), &
                                       nxp*lev*my_max*nsize*batch_size, &
                                       cudaMemcpyHostToDevice, cudast(1))
               istat = cudaEventRecord(cudaEV(4), cudast(1))
               istat = cudaStreamWaitEvent(cudast(2), cudaEV(4), 0)
               istat = cudaStreamWaitEvent(cudast(2), cudaEV(1), 0)

               call advh_reshape_c2g_q(qq_3df, ainp_d, &
                                       batch_size, 2)

               call ndslfv_monoadvh2_xy_gpu(qq_3df, deltim, um3df, vm3df, &
                                            wndmdf(n:n + batch_size - 1), forward, &
                                            lev, batch_size, 2)

               call advh_reshape_g2c_q(aout_d, qq_3df, &
                                       batch_size, 2)

               istat = cudaEventRecord(cudaEV(1), cudast(2))
               istat = cudaStreamWaitEvent(cudast(cudast_D2H), cudaEV(1), 0)
               istat = cudamemcpyasync(advhwrk(1, 1, n), aout_d, &
                                       nxp*lev*my_max*nsize*batch_size, &
                                       cudaMemcpyDeviceToHost, cudast(cudast_D2H))
            end do
         elseif (xy .lt. -0.5) then
            do n = 1, num_ndslhvar_gpu, batch_size
               istat = cudamemcpyasync(ainp_d, advhwrk(1, 1, n), &
                                       nxp*lev*my_max*nsize*batch_size, &
                                       cudaMemcpyHostToDevice, cudast(1))
               istat = cudaEventRecord(cudaEV(4), cudast(1))
               istat = cudaStreamWaitEvent(cudast(2), cudaEV(4), 0)
               istat = cudaStreamWaitEvent(cudast(n), cudaEV(1), 0)
               call advh_reshape_c2g_q(qq_3df, ainp_d, &
                                       batch_size, 2)

               !
               call ndslfv_monoadvh2_yx_gpu(qq_3df, deltim, um3df, vm3df, &
                                            wndmdf(n:n + batch_size - 1), forward, &
                                            lev, batch_size, 2)

               call advh_reshape_g2c_q(aout_d, qq_3df, &
                                       batch_size, 2)

               istat = cudaEventRecord(cudaEV(1), cudast(2))
               istat = cudaStreamWaitEvent(cudast(cudast_D2H), cudaEV(1), 0)
               istat = cudamemcpyasync(advhwrk(1, 1, n), aout_d, &
                                       nxp*lev*my_max*nsize*batch_size, &
                                       cudaMemcpyDeviceToHost, cudast(cudast_D2H))

            end do
         end if
         !$acc wait
         !$acc end data
      end if
      ! call mpe_barrier
      ! ----------------------------------------
      ! u v t at n
      call advh_scatter4GPU_q(uum_sl, advhwrk(1, 1, mod(0, num_ndslhvar_gpu) + 1), &
                              1, 1, mod(0, use_gpu_num))
      call advh_scatter4GPU_q(vvm_sl, advhwrk(1, 1, mod(1, num_ndslhvar_gpu) + 1), &
                              1, 1, mod(1, use_gpu_num))
      call advh_scatter4GPU_q(ttm_sl, advhwrk(1, 1, mod(2, num_ndslhvar_gpu) + 1), &
                              1, 1, mod(2, use_gpu_num))
      ! rq at n
      call advh_scatter4GPU_q(qm_sl, advhwrk(1, 1, mod(3, num_ndslhvar_gpu) + 1), &
                              1, 3, mod(3, use_gpu_num))
      call advh_scatter4GPU_q(qm_sl, advhwrk(1, 1, mod(4, num_ndslhvar_gpu) + 1), &
                              2, 3, mod(4, use_gpu_num))
      call advh_scatter4GPU_q(qm_sl, advhwrk(1, 1, mod(5, num_ndslhvar_gpu) + 1), &
                              3, 3, mod(5, use_gpu_num))
      ! ----------------------------------------

      return
   end subroutine ndslfv_monoadvh_gpu_batch
   ! ============================================================
   subroutine ndslfv_monoadvh2_gpu(qqlon, deltim, um, vm, &
                                   wndmdf, forward, levs, nvars, async)
      !
      ! a routine to do non-iteration semi-Lagrangain advection
      ! considering advection  with monotonicity in interpolation
      ! contact: hann-ming henry juang
      ! program log
      ! 2011 02 20 : henry juang, created for ndsl advection
      ! 2013 06 20 : Henry Juang correct wind direction for north-south advection
      !
      use param
      use index
      use const
      use rank, only: nsize
      implicit none

      real(kind=RTYPE), intent(inout) :: qqlon(nx, my, levs, nvars)
      real(kind=RTYPE), intent(in):: deltim
      real(kind=RTYPE):: um(nx, my, levs), vm(nx, my, levs)
      integer, intent(in):: wndmdf(nvars)
      integer, intent(in):: levs, nvars
      logical, intent(in):: forward
      integer, intent(in):: async
      ! local
      real(kind=RTYPE):: rrlon(nx, my, levs, nvars), &
                         rrlat(my*2, nx/2, levs, nvars), &
                         qqlat(my*2, nx/2, levs, nvars)
      real(kind=RTYPE) rma, rm2a
      integer one(1)
      integer mass
      integer i, j, n, k, istat
      integer jf, ii, m, nn, kk, pts
      mass = 0
      !
      ! =================================================================
      !   prepare wind and variable in flux form with gaussina weight
      ! =================================================================
      !
      ! wind at time step n

      !$acc data async(async) &
      !$acc& copyin(levs,nvars,deltim) &
      !$acc& create(rrlon(:nx, :my, :levs, :nvars), &
      !$acc& rrlat(:my*2,:nx/2,:levs,:nvars), &
      !$acc& qqlat(:my*2,:nx/2,:levs,:nvars), &
      !$acc& wndmdf(:nvars),one(:1) ) &
      !$acc& present( &
      !$acc& my,lev,ncld,nxp,nsizex, &
      !$acc& jlist2_2d(nsizex,my), &
      !$acc& nxjlen_all(nsizex,my) )
      !$acc update device(wndmdf(:nvars)) async(async)

      if (firstcomp) then
         one = 1
         !$acc update device(one(1)) async(cudast_advh_v)
         !$acc parallel loop collapse(2) async(cudast_advh_v)&
         !$acc& private(ii,jf,nn,kk,i,m,rm2a,rma)
         do k = 1, lev
            do j = 1, my
               ii = 0
               rma = 1./cosl(j)
               !$acc loop seq
               do i = 1, nsizex
                  jf = jlist2_2d(i, j)
                  nn = nxjlen_all(i, j)
                  kk = (k - 1)*nxp
                  !$acc loop independent
                  do m = 1, nn
                     ! wind at time step n
                     vm(ii + m, j, k) = vmwrk(kk + m, jf)*rma
                  end do
                  ii = ii + nn
               end do
            end do
         end do
         call cyclic_cell_intpx_r2f_gpu(vvlat, vm, one, &
                                        nxdef, my, nx, my, levs, one(1), &
                                        cudast_advh_v)
         istat = cudaEventRecord(cudaEV(3), cudast(cudast_advh_v))

         !$acc parallel loop collapse(2) async(cudast_advh_u)&
         !$acc& private(ii,jf,nn,kk,i,m,rm2a,rma)
         do k = 1, lev
            do j = 1, my
               ii = 0
               rma = 1./cosl(j)
               rm2a = rma/cosl(j)
               !$acc loop seq
               do i = 1, nsizex
                  jf = jlist2_2d(i, j)
                  nn = nxjlen_all(i, j)
                  kk = (k - 1)*nxp
                  !$acc loop independent
                  do m = 1, nn
                     ! wind at time step n
                     um(ii + m, j, k) = umwrk(kk + m, jf)*rm2a
                  end do
                  ii = ii + nn
               end do
            end do
         end do
         istat = cudaEventRecord(cudaEV(2), cudast(cudast_advh_u))

         firstcomp = .false.
      end if

      !$acc kernels present(qqlon(nx,my,levs,nvars), &
      !$acc& rrlon(nx,my,levs,nvars))
      rrlon = qqlon
      !$acc end kernels

      istat = cudaStreamWaitEvent(cudast(async), cudaEV(2), 0)
      call cyclic_cell_massadvx_gpu(rrlon, um, deltim, mass, forward, &
                                    nxdef, my, nx, my, levs, nvars, &
                                    async)

      ! ---------------------------------------------------------------------
      ! from east-west reduced grid to north-south full grid
      ! ---------------------------------------------------------------------
      call cyclic_cell_intpx_r2f_gpu(rrlat, rrlon, wndmdf, &
                                     nxdef, my, nx, my, levs, nvars, &
                                     async)
      call cyclic_cell_intpx_r2f_gpu(qqlat, qqlon, wndmdf, &
                                     nxdef, my, nx, my, levs, nvars, &
                                     async)

      istat = cudaStreamWaitEvent(cudast(async), cudaEV(3), 0)
      call cyclic_cell_massadvy_gpu(rrlat, vvlat, deltim, mass, forward, &
                                    my*2, nx/2, my*2, nx/2, levs, nvars, &
                                    async)
      call cyclic_cell_massadvy_gpu(qqlat, vvlat, deltim, mass, forward, &
                                    my*2, nx/2, my*2, nx/2, levs, nvars, &
                                    async)
      ! ----------------------------------------------------------------------
      ! from north-south full grid to east-west reduced grid
      ! ----------------------------------------------------------------------
      call cyclic_cell_intpx_f2r_gpu(rrlon, qqlat, wndmdf, &
                                     nxdef, my, nx, my, levs, nvars, &
                                     async)
      call cyclic_cell_intpx_f2r_gpu(qqlon, qqlat, wndmdf, &
                                     nxdef, my, nx, my, levs, nvars, &
                                     async)

      call cyclic_cell_massadvx_gpu(qqlon, um, deltim, mass, forward, &
                                    nxdef, my, nx, my, levs, nvars, &
                                    async)

      !$acc kernels present(qqlon(nx,my,levs,nvars), rrlon(nx,my,levs,nvars))
      qqlon = .5*(qqlon + rrlon)
      !$acc end kernels
      !$acc end data

      return
   end subroutine ndslfv_monoadvh2_gpu
   ! ============================================================
   subroutine ndslfv_monoadvh2_xy_gpu(qqlon, deltim, um, vm, &
                                      wndmdf, forward, levs, nvars, async)
      !
      ! a routine to do non-iteration semi-Lagrangain advection
      ! considering advection  with monotonicity in interpolation
      ! contact: hann-ming henry juang
      ! program log
      ! 2011 02 20 : henry juang, created for ndsl advection
      ! 2013 06 20 : Henry Juang correct wind direction for north-south advection
      !
      use param
      use index
      use const
      use rank, only: nsize, myrank
      implicit none

      real(kind=RTYPE), intent(inout) :: qqlon(nx, my, levs, nvars)
      real(kind=RTYPE), intent(in):: deltim
      real(kind=RTYPE):: um(nx, my, levs), vm(nx, my, levs)
      integer, intent(in):: wndmdf(nvars)
      integer, intent(in):: levs, nvars
      logical, intent(in):: forward
      integer, intent(in):: async
      ! local
      real(kind=RTYPE) qqlat(my*2, nx/2, levs, nvars)
      real(kind=RTYPE) rma, rm2a
      integer one(1)
      integer mass
      integer i, j, n, k, istat
      integer jf, ii, m, nn, kk, pts
      mass = 0
      !
      ! =================================================================
      !   prepare wind and variable in flux form with gaussina weight
      ! =================================================================
      !
      ! wind at time step n

      !$acc data async(async) &
      !$acc& copyin(levs,nvars,deltim) &
      !$acc& create(qqlat(:my*2,:nx/2,:levs,:nvars), &
      !$acc& wndmdf(:nvars),one(:1) ) &
      !$acc& present( &
      !$acc& my,lev,ncld,nxp,nsizex, &
      !$acc& jlist2_2d(nsizex,my), &
      !$acc& nxjlen_all(nsizex,my) )
      !$acc update device(wndmdf(:nvars)) async(async)

      if (firstcomp) then
         one = 1
         !$acc update device(one(1)) async(cudast_advh_v)
         !$acc parallel loop collapse(2) async(cudast_advh_v)&
         !$acc& private(ii,jf,nn,kk,i,m,rm2a,rma)
         do k = 1, lev
            do j = 1, my
               ii = 0
               rma = 1./cosl(j)
               !$acc loop seq
               do i = 1, nsizex
                  jf = jlist2_2d(i, j)
                  nn = nxjlen_all(i, j)
                  kk = (k - 1)*nxp
                  !$acc loop independent
                  do m = 1, nn
                     ! wind at time step n
                     vm(ii + m, j, k) = vmwrk(kk + m, jf)*rma
                  end do
                  ii = ii + nn
               end do
            end do
         end do
         call cyclic_cell_intpx_r2f_gpu(vvlat, vm, one, &
                                        nxdef, my, nx, my, levs, one(1), &
                                        cudast_advh_v)
         istat = cudaEventRecord(cudaEV(3), cudast(cudast_advh_v))

         !$acc parallel loop collapse(2) async(cudast_advh_u)&
         !$acc& private(ii,jf,nn,kk,i,m,rm2a,rma)
         do k = 1, lev
            do j = 1, my
               ii = 0
               rma = 1./cosl(j)
               rm2a = rma/cosl(j)
               !$acc loop seq
               do i = 1, nsizex
                  jf = jlist2_2d(i, j)
                  nn = nxjlen_all(i, j)
                  kk = (k - 1)*nxp
                  !$acc loop independent
                  do m = 1, nn
                     ! wind at time step n
                     um(ii + m, j, k) = umwrk(kk + m, jf)*rm2a
                  end do
                  ii = ii + nn
               end do
            end do
         end do
         istat = cudaEventRecord(cudaEV(2), cudast(cudast_advh_u))

         firstcomp = .false.
      end if

      istat = cudaStreamWaitEvent(cudast(async), cudaEV(2), 0)
      call cyclic_cell_massadvx_gpu(qqlon, um, deltim, mass, forward, &
                                    nxdef, my, nx, my, levs, nvars, &
                                    async)

      ! ---------------------------------------------------------------------
      ! from east-west reduced grid to north-south full grid
      ! ---------------------------------------------------------------------
      call cyclic_cell_intpx_r2f_gpu(qqlat, qqlon, wndmdf, &
                                     nxdef, my, nx, my, levs, nvars, &
                                     async)

      istat = cudaStreamWaitEvent(cudast(async), cudaEV(3), 0)
      call cyclic_cell_massadvy_gpu(qqlat, vvlat, deltim, mass, forward, &
                                    my*2, nx/2, my*2, nx/2, levs, nvars, &
                                    async)
      ! ----------------------------------------------------------------------
      ! from north-south full grid to east-west reduced grid
      ! ----------------------------------------------------------------------
      call cyclic_cell_intpx_f2r_gpu(qqlon, qqlat, wndmdf, &
                                     nxdef, my, nx, my, levs, nvars, &
                                     async)
      !$acc end data

      return
   end subroutine ndslfv_monoadvh2_xy_gpu
   ! ============================================================
   subroutine ndslfv_monoadvh2_yx_gpu(qqlon, deltim, um, vm, &
                                      wndmdf, forward, levs, nvars, async)
      !
      ! a routine to do non-iteration semi-Lagrangain advection
      ! considering advection  with monotonicity in interpolation
      ! contact: hann-ming henry juang
      ! program log
      ! 2011 02 20 : henry juang, created for ndsl advection
      ! 2013 06 20 : Henry Juang correct wind direction for north-south advection
      !
      use param
      use index
      use const
      use rank, only: nsize
      implicit none

      real(kind=RTYPE), intent(inout) :: qqlon(nx, my, levs, nvars)
      real(kind=RTYPE), intent(in):: deltim
      real(kind=RTYPE) um(nx, my, levs), vm(nx, my, levs)
      integer, intent(in):: wndmdf(nvars)
      integer, intent(in):: levs, nvars
      logical, intent(in):: forward
      integer, intent(in):: async
      ! local
      real(kind=RTYPE) qqlat(my*2, nx/2, levs, nvars) ! , vvlat(my*2,nx/2,levs)
      real(kind=RTYPE) rma, rm2a
      integer one(1)
      integer mass
      integer i, j, n, k, istat
      integer jf, ii, m, nn, kk, pts
      mass = 0

      !$acc data async(async) &
      !$acc& present_or_copyin(nx,my,levs,nvars) &
      !$acc& create(qqlat(:my*2,:nx/2,:levs,:nvars), &
      !$acc& wndmdf(:nvars),one(:1) )
      !$acc update device(wndmdf(:nvars)) async(async)

      ! ---------------------------------------------------------------------
      ! from east-west reduced grid to north-south full grid
      ! ---------------------------------------------------------------------
      call cyclic_cell_intpx_r2f_gpu(qqlat, qqlon, wndmdf, &
                                     nxdef, my, nx, my, levs, nvars, &
                                     async)

      if (firstcomp) then
         one = 1.
         !$acc update device(one(1)) async(cudast_advh_v)
         !$acc parallel loop collapse(2) async(cudast_advh_v)&
         !$acc& private(ii,jf,nn,kk,i,m,rm2a,rma)
         do k = 1, lev
            do j = 1, my
               ii = 0
               rma = 1./cosl(j)
               rm2a = rma/cosl(j)
               !$acc loop seq
               do i = 1, nsizex
                  jf = jlist2_2d(i, j)
                  nn = nxjlen_all(i, j)
                  kk = (k - 1)*nxp
                  !$acc loop independent
                  do m = 1, nn
                     ! wind at time step n
                     vm(ii + m, j, k) = vmwrk(kk + m, jf)*rma
                  end do
                  ii = ii + nn
               end do
            end do
         end do
         call cyclic_cell_intpx_r2f_gpu(vvlat, vm, one, &
                                        nxdef, my, nx, my, levs, one(1), &
                                        cudast_advh_v)
         istat = cudaEventRecord(cudaEV(3), cudast(cudast_advh_v))

         !$acc parallel loop collapse(2) async(cudast_advh_u)&
         !$acc& private(ii,jf,nn,kk,i,m,rm2a,rma)
         do k = 1, lev
            do j = 1, my
               ii = 0
               rma = 1./cosl(j)
               rm2a = rma/cosl(j)
               !$acc loop seq
               do i = 1, nsizex
                  jf = jlist2_2d(i, j)
                  nn = nxjlen_all(i, j)
                  kk = (k - 1)*nxp
                  !$acc loop independent
                  do m = 1, nn
                     ! wind at time step n
                     um(ii + m, j, k) = umwrk(kk + m, jf)*rm2a
                  end do
                  ii = ii + nn
               end do
            end do
         end do
         istat = cudaEventRecord(cudaEV(2), cudast(cudast_advh_u))

         firstcomp = .false.
      end if

      istat = cudaStreamWaitEvent(cudast(async), cudaEV(3), 0)
      call cyclic_cell_massadvy_gpu(qqlat, vvlat, deltim, mass, forward, &
                                    my*2, nx/2, my*2, nx/2, levs, nvars, &
                                    async)
      ! ----------------------------------------------------------------------
      ! from north-south full grid to east-west reduced grid
      ! ----------------------------------------------------------------------
      call cyclic_cell_intpx_f2r_gpu(qqlon, qqlat, wndmdf, &
                                     nxdef, my, nx, my, levs, nvars, &
                                     async)
    !!     !$acc wait
      istat = cudaStreamWaitEvent(cudast(async), cudaEV(2), 0)
      call cyclic_cell_massadvx_gpu(qqlon, um, deltim, mass, forward, &
                                    nxdef, my, nx, my, levs, nvars, &
                                    async)
      !$acc end data

      return
   end subroutine ndslfv_monoadvh2_yx_gpu
   ! ============================================================
   subroutine cyclic_cell_massadvx_gpu(qq, uc, delt, mass, forward, &
                                       im, jm, lons, lats, levs, nv, async)
      !
      ! compute local positive advection with mass conservation
      ! qq is advected by uc from past to next position
      !
      ! author: hann-ming henry juang 2008
      !
      use const, only: RTYPE
      use rank, only: myrank
      implicit none
      !
      real(kind=RTYPE), intent(inout):: qq(lons, lats, levs, nv)
      real(kind=RTYPE), intent(in)::    uc(lons, lats, levs)
      integer, intent(in):: im(lats), jm, lons, lats, levs, nv, mass
      real(kind=RTYPE), intent(in):: delt
      logical forward
      integer, intent(in):: async
      ! local
      real(kind=RTYPE) da(lons, lats, levs, nv)
      real(kind=RTYPE) dxfact(lons)
      real(kind=RTYPE) xreg(lons + 1, lats, levs)
      real(kind=RTYPE) xpast(lons + 1, lats, levs), xnext(lons + 1, lats, levs)
      real(kind=RTYPE) uint(lons + 1)
      real(kind=RTYPE) dist(lons + 1, lats, levs)
      real(kind=RTYPE) ds(lons + 1, lats), step(10), dist_step
      real(kind=RTYPE) wrk_step(10)
      real(kind=RTYPE), parameter :: fa1 = 9./16.
      real(kind=RTYPE), parameter :: fa2 = 1./16.
      real(kind=RTYPE):: sc

      integer imp
      integer i, j, k, n, nst, nstep, wrk_ns

      ! sc = 2.0 * pi
      !$acc data async(async) &
      !$acc& present_or_copyin(im(:lats),jm,levs,sc,nv,two_pi) &
      !$acc& create( &
      !$acc& da(:lons,:lats,:levs,:nv), sc, &
      !$acc& xpast(:lons+1,:lats,:levs), &
      !$acc& xnext(:lons+1,:lats,:levs), &
      !$acc& xreg(:lons+1,:lats,:levs), &
      !$acc& ds(:lons+1,:lats), dist(:lons+1,:lats,:levs), &
      !$acc& nstep,step(:10)) &
      !$acc& present( &
      !$acc& uc(lons,lats,levs), &
      !$acc& delt &
      !$acc& )
      !----------------------------------------
      !$acc parallel loop collapse(2) async(async) &
      !$acc& private(imp,i,uint(:lons+1))
      do k = 1, levs
         do j = 1, jm
            imp = im(j)
            !
            ! 4th order interpolation from mid point to cell interfaces
            !
            do i = 3, imp - 1
               uint(i) = fa1*(uc(i, j, k) + uc(i - 1, j, k)) &
                         - fa2*(uc(i + 1, j, k) + uc(i - 2, j, k))
            end do
            uint(2) = fa1*(uc(2, j, k) + uc(1, j, k)) &
                      - fa2*(uc(3, j, k) + uc(imp, j, k))
            uint(1) = fa1*(uc(1, j, k) + uc(imp, j, k)) &
                      - fa2*(uc(2, j, k) + uc(imp - 1, j, k))

            uint(imp + 1) = uint(1)
            uint(imp) = fa1*(uc(imp, j, k) + uc(imp - 1, j, k)) &
                        - fa2*(uc(1, j, k) + uc(imp - 2, j, k))

            do i = 1, imp + 1
               dist(i, j, k) = uint(i)*delt
            end do
         end do
      end do
      ! ----------------------------------------
      !$acc parallel loop async(async) &
      !$acc& private(imp,i)
      do j = 1, jm
         imp = im(j)
         do i = 1, imp + 1
            ds(i, j) = two_pi/float(imp)
         end do
      end do

      nstep = 0
      ! <<<< -------------------------------------------------------
      call def_cfl_step_gpu(step, nstep, &
                            dist, ds, 'advx', &
                            im, lons, lats, levs, async)

      ! ============================================================
      ! !$acc update self(dist(:lons+1,:lats,:levs), ds(:lons+1,:lats))
      ! do k = 1, levs
      !    do j = 1, jm
      !       ! print*, "dist=", maxval(dist(1:im(j)+1,j,k)), &
      !       !      maxval(ds(1:im(j),j))
      !       call def_cfl_step(im(j) + 1, dist(1, j, k), ds(1, j), &
      !            wrk_step, wrk_ns, levs + 1 - k, 'advx')
      !       if (nstep .lt. wrk_ns) then
      !          nstep = wrk_ns
      !          step(1:10) = wrk_step(1:10)
      !       end if
      !    end do
      ! end do
      ! !$acc update device(step(:10))
      ! ------------------------------------------------------- >>>>

      do nst = 1, nstep
         if (forward) then
            !$acc parallel loop collapse(2) copyin(nst) async(async)&
            !$acc& private(imp,i,dist_step)
            do k = 1, levs
               do j = 1, jm
                  imp = im(j)
                  do i = 1, imp + 1
                     xreg(i, j, k) = (i - 1.5)*ds(i, j)
                  end do
                  !
                  do i = 1, imp + 1
                     dist_step = dist(i, j, k)*step(nst)
                     xpast(i, j, k) = xreg(i, j, k)
                     xnext(i, j, k) = xreg(i, j, k) + dist_step
                  end do
               end do
            end do
            !$acc kernels async(async)
            da = qq
            !$acc end kernels
         else
            !$acc parallel loop collapse(2) copyin(nst) async(async)&
            !$acc& private(imp,i,dist_step)
            do k = 1, levs
               do j = 1, jm
                  imp = im(j)
                  do i = 1, imp + 1
                     xreg(i, j, k) = (i - 1.5)*ds(i, j)
                  end do
                  !
                  do i = 1, imp + 1
                     dist_step = dist(i, j, k)*step(nst)
                     xpast(i, j, k) = xreg(i, j, k) - dist_step
                     xnext(i, j, k) = xreg(i, j, k) + dist_step
                  end do
               end do
            end do
            call cyclic_cell_ppm_intp_gpu_v2(xreg, qq, xpast, da, &
                                             im, im, jm, lons, lats, levs, nv, two_pi, &
                                             async)
         end if

         if (mass .eq. 1) then
            !$acc parallel loop collapse(3) async(2)&
            !$acc& private(imp,i,dxfact(:lons))
            do n = 1, nv
               do k = 1, levs
                  do j = 1, jm
                     imp = im(j)
                     do i = 1, imp + 1
                        dxfact(i) = (xpast(i + 1, j, k) - xpast(i, j, k)) &
                                    /(xnext(i + 1, j, k) - xnext(i, j, k))

                        da(i, j, k, n) = da(i, j, k, n)*dxfact(i)
                     end do
                  end do
               end do
            end do
         end if
         call cyclic_cell_ppm_intp_gpu(xnext, da, xreg, qq, &
                                       im, im, jm, lons, lats, levs, nv, two_pi, &
                                       async)

      end do
      !$acc end data
      return
   end subroutine cyclic_cell_massadvx_gpu
   ! ============================================================
   subroutine cyclic_cell_massadvx_LBL_gpu(qq, uc, delt, mass, forward, &
                                           im, jm, lons, lats, levs, nv, async, typ)
      !
      ! compute local positive advection with mass conservation
      ! qq is advected by uc from past to next position
      !
      ! author: hann-ming henry juang 2008
      !
      ! --------------------
      ! For debugging
      ! --------------------
      use const, only: RTYPE
      use rank, only: myrank
      implicit none
      !
      real(kind=RTYPE), intent(inout):: qq(lons, lats, levs, nv)
      real(kind=RTYPE), intent(in)::    uc(lons, lats, levs)
      integer, intent(in):: im(lats), jm, lons, lats, levs, nv, mass
      real(kind=RTYPE), intent(in):: delt
      logical forward
      integer, intent(in):: async, typ
      ! local
      real(kind=RTYPE) da(lons, lats, levs, nv)
      real(kind=RTYPE) dxfact(lons)
      real(kind=RTYPE) xreg(lons + 1, lats, levs)
      real(kind=RTYPE) xpast(lons + 1, lats, levs), xnext(lons + 1, lats, levs)
      real(kind=RTYPE) uint(lons + 1)
      real(kind=RTYPE) dist(lons + 1, lats, levs)
      real(kind=RTYPE) ds(lons + 1, lats), step(10), dist_step
      real(kind=RTYPE) wrk_step(10)
      real(kind=RTYPE), parameter :: fa1 = 9./16.
      real(kind=RTYPE), parameter :: fa2 = 1./16.
      real(kind=RTYPE):: sc
      real(kind=RTYPE) qq_wrk(lons, nv), xreg_wrk(lons + 1), xpast_wrk(lons + 1)
      real(kind=RTYPE) da_wrk(lons, nv), xnext_wrk(lons + 1)

      integer imp
      integer i, j, k, n, nst, nstep, wrk_ns

      ! sc = 2.0 * pi
      !$acc data async(async) &
      !$acc& present_or_copyin(im(:lats),jm,levs,sc,nv,two_pi) &
      !$acc& create( &
      !$acc& da(:lons,:lats,:levs,:nv), sc, &
      !$acc& xpast(:lons+1,:lats,:levs), &
      !$acc& xnext(:lons+1,:lats,:levs), &
      !$acc& xreg(:lons+1,:lats,:levs), &
      !$acc& ds(:lons+1,:lats), dist(:lons+1,:lats,:levs), &
      !$acc& qq_wrk, da_wrk, xreg_wrk, xpast_wrk, xnext_wrk, &
      !$acc& nstep,step(:10)) &
      !$acc& present( &
      !$acc& uc(lons,lats,levs), &
      !$acc& delt &
      !$acc& )
      !----------------------------------------
      !$acc parallel loop collapse(2) async(async) &
      !$acc& private(imp,i,uint(:lons+1))
      do k = 1, levs
         do j = 1, jm
            imp = im(j)
            !
            ! 4th order interpolation from mid point to cell interfaces
            !
            do i = 3, imp - 1
               uint(i) = fa1*(uc(i, j, k) + uc(i - 1, j, k)) &
                         - fa2*(uc(i + 1, j, k) + uc(i - 2, j, k))
            end do
            uint(2) = fa1*(uc(2, j, k) + uc(1, j, k)) &
                      - fa2*(uc(3, j, k) + uc(imp, j, k))
            uint(1) = fa1*(uc(1, j, k) + uc(imp, j, k)) &
                      - fa2*(uc(2, j, k) + uc(imp - 1, j, k))

            uint(imp + 1) = uint(1)
            uint(imp) = fa1*(uc(imp, j, k) + uc(imp - 1, j, k)) &
                        - fa2*(uc(1, j, k) + uc(imp - 2, j, k))

            do i = 1, imp + 1
               dist(i, j, k) = uint(i)*delt
            end do
         end do
      end do
      ! ----------------------------------------
      !$acc parallel loop async(async) &
      !$acc& private(imp,i)
      do j = 1, jm
         imp = im(j)
         do i = 1, imp + 1
            ds(i, j) = two_pi/float(imp)
         end do
      end do

      nstep = 0
      if (typ .ne. 0) then
         call def_cfl_step_gpu(step, nstep, &
                               dist, ds, 'advx', &
                               im, lons, lats, levs, async)
      end if

      do k = 1, levs
         do j = 1, jm
            if (typ .eq. 0) then
               call def_cfl_step_gpu(step, nstep, &
                                     dist(1, j, k), ds(1, j), 'advx', &
                                     (/im(j)/), lons, 1, 1, async)
               ! else
               !    nstep = typ
               !    do nst = 1,nstep
               !       step(nst) = 1./nstep
               !    end do
               !    !$acc  update device(step)
            end if

            do nst = 1, nstep
               imp = im(j)
               !$acc parallel loop copyin(imp,nst) &
               !$acc& private(dist_step)
               do i = 1, imp + 1
                  xreg(i, j, k) = (i - 1.5)*ds(i, j)

                  dist_step = dist(i, j, k)*step(nst)
                  xpast(i, j, k) = xreg(i, j, k) - dist_step
                  xnext(i, j, k) = xreg(i, j, k) + dist_step
               end do

               !$acc kernels
               do n = 1, nv
                  do i = 1, imp
                     qq_wrk(i, n) = qq(i, j, k, n)
                  end do

                  do i = 1, imp + 1
                     xreg_wrk(i) = xreg(i, j, k)
                     xpast_wrk(i) = xpast(i, j, k)
                     xnext_wrk(i) = xnext(i, j, k)
                  end do

               end do
               !$acc end kernels
               call cyclic_cell_ppm_intp_gpu_v2(xreg_wrk, qq_wrk, &
                                                xpast_wrk, da_wrk, &
                                                (/imp/), (/imp/), 1, lons, 1, 1, nv, two_pi, &
                                                async)

               call cyclic_cell_ppm_intp_gpu(xnext_wrk, da_wrk, &
                                             xreg_wrk, qq_wrk, &
                                             (/imp/), (/imp/), 1, lons, 1, 1, nv, two_pi, &
                                             async)

               !$acc kernels
               do n = 1, nv
                  do i = 1, imp
                     qq(i, j, k, n) = qq_wrk(i, n)
                  end do
               end do
               !$acc end kernels

            end do ! do nst
         end do ! do jm
      end do ! do levs
      !$acc end data
      return
   end subroutine cyclic_cell_massadvx_LBL_gpu
   ! ============================================================
   subroutine cyclic_cell_massadvy_gpu(qq, vc, delt, mass, forward, &
                                       jm, im, latf, lons, levs, nv, async)
      !
      ! compute local positive advection with mass conserving
      ! qq will be advect by vc from past to next location with 2*delt
      !
      ! author: hann-ming henry juang 2007
      !
      !
      use grid, only: gglati, fa1, fa2, fa3, fa4
      use const, only: RTYPE

      implicit none
      !
      real(kind=RTYPE), intent(inout):: qq(latf, lons, levs, nv)
      real(kind=RTYPE), intent(in):: vc(latf, lons, levs)
      real(kind=RTYPE), intent(in):: delt
      integer, intent(in)::jm, im, latf, lons, levs, nv, mass
      logical, intent(in):: forward
      integer, intent(in):: async
      ! local
      integer jmh
      !
      real(kind=RTYPE):: var(latf), &
                         past(latf, lons, levs, nv), &
                         next(latf, lons, levs, nv), &
                         da(latf, lons, levs, nv), &
                         dyfact, &
                         ypast(latf + 1, lons, levs), &
                         ynext(latf + 1, lons, levs), &
                         dist(latf + 1, lons, levs), &
                         ds(latf + 1, lons), &
                         step(10), dist_step, &
                         wrk_step(10), &
                         sc, &
                         ygrd(latf + 1, lons, levs)

      integer i, n, k, j, latfh, nst, nstep, wrk_ns
      !
      ! preparations ---------------------------
      !
      sc = gglati(jm + 1) - gglati(1)
      !$acc data async(async) &
      !$acc& present_or_copyin(im,jm,latf,lons,levs,nv) &
      !$acc& create( &
      !$acc& da(:latf,:lons,:levs,:nv), &
      !$acc& dist(:latf+1,:lons,:levs), &
      !$acc& ypast(:latf+1,:lons,:levs),ynext(:latf+1,:lons,:levs), &
      !$acc& ygrd(:latf+1,:lons,:levs),ds(:latf+1,:lons), &
      !$acc& step(:10), nstep) &
      !$acc& present( delt, &
      !$acc& vc(latf,lons,levs), &
      !$acc& gglati(latf+1), &
      !$acc& fa1(latf),fa2(latf),fa3(latf),fa4(latf) &
      !$acc& )
      ! ----------------------------------------
      !$acc parallel loop collapse(2) async(async) &
      !$acc& private(j,var(:latf),jmh)
      do k = 1, levs
         do i = 1, im
            jmh = jm/2
            do j = 1, jmh
               var(j) = vc(j, i, k)*delt
               var(j + jmh) = -vc(j + jmh, i, k)*delt
            end do
            ! for Gaussian latitude
            do j = 3, jm - 1
               dist(j, i, k) = fa1(j)*var(j - 2) + fa2(j)*var(j - 1) &
                               + fa3(j)*var(j) + fa4(j)*var(j + 1)
            end do
            ! over pole
            dist(2, i, k) = fa1(2)*var(jm) + fa2(2)*var(1) &
                            + fa3(2)*var(2) + fa4(2)*var(3)
            dist(jm, i, k) = fa1(jm)*var(jm - 2) + fa2(jm)*var(jm - 1) &
                             + fa3(jm)*var(jm) + fa4(jm)*var(1)

            dist(1, i, k) = fa1(1)*var(jm - 1) + fa2(1)*var(jm) &
                            + fa3(1)*var(1) + fa4(1)*var(2)
            dist(jm + 1, i, k) = dist(1, i, k)
         end do
      end do

      !$acc parallel loop collapse(2) async(async) &
      !$acc& present(jm,im,ds(latf+1,lons))
      do i = 1, im
         do j = 1, jm
            ds(j, i) = gglati(j + 1) - gglati(j)
         end do
      end do

      nstep = 0
      ! <<<< -------------------------------------------------------
      call def_cfl_step_gpu(step, nstep, &
                            dist, ds, 'advy', &
                            nydef_loc, latf, lons, levs, async)
      ! ============================================================
      ! !$acc update self(dist(:latf+1,:lons,:levs), ds(:latf+1,:lons))
      ! do k =1,levs
      !    do i =1,im
      !       ! print*, "dist=", maxval(dist(1:jm+1,i,k)), &
      !       !      maxval(ds(1:jm,i))
      !       call def_cfl_step(jm+1,dist(1,i,k),ds, &
      !            wrk_step,wrk_ns,levs+1-k,'advy')
      !       if(nstep.lt.wrk_ns) then
      !          nstep = wrk_ns
      !          step(1:10) = wrk_step(1:10)
      !       end if
      !    end do
      ! enddo
      ! !$acc update device(step(:10))
      ! ------------------------------------------------------- >>>>

      !
      ! advection all in y
      !
      do nst = 1, nstep
         if (forward) then
            !$acc parallel loop collapse(2) copyin(nst) async(async) &
            !$acc& private(j,dist_step)
            do k = 1, levs
               do i = 1, im
                  do j = 1, jm + 1
                     ygrd(j, i, k) = gglati(j)
                  end do

                  do j = 1, jm + 1
                     dist_step = dist(j, i, k)*step(nst)
                     ypast(j, i, k) = ygrd(j, i, k)
                     ynext(j, i, k) = ygrd(j, i, k) + dist_step
                  end do
               end do
            end do
            !$acc kernels
            da = qq
            !$acc end kernels
         else
            !$acc parallel loop collapse(2) copyin(nst) async(async) &
            !$acc& private(j,dist_step)
            do k = 1, levs
               do i = 1, im
                  do j = 1, jm + 1
                     ygrd(j, i, k) = gglati(j)
                  end do

                  do j = 1, jm + 1
                     dist_step = dist(j, i, k)*step(nst)
                     ypast(j, i, k) = ygrd(j, i, k) - dist_step
                     ynext(j, i, k) = ygrd(j, i, k) + dist_step
                  end do
               end do
            end do
            call cyclic_cell_ppm_intp_gpu(ygrd, qq, ypast, da, &
                                          nydef_loc, nydef_loc, &
                                          im, latf, lons, levs, nv, sc_advy, &
                                          async)
         end if
         if (mass .eq. 1) then
            !$acc parallel loop collapse(4) async(async)&
            !$acc& present(im,jm,levs,nv, da(latf,lons,levs,nv)) &
            !$acc& private(dyfact)
            do n = 1, nv
               do k = 1, levs
                  do i = 1, im
                     do j = 1, jm
                        dyfact = (ypast(j + 1, i, k) - ypast(j, i, k)) &
                                 /(ynext(j + 1, i, k) - ynext(j, i, k))
                        da(j, i, k, n) = da(j, i, k, n)*dyfact
                     end do
                  end do
               end do
            end do
         end if
         !
         call cyclic_cell_ppm_intp_gpu(ynext, da, ygrd, qq, &
                                       nydef_loc, nydef_loc, &
                                       im, latf, lons, levs, nv, sc_advy, &
                                       async)
      end do
      !$acc end data
      return
   end subroutine cyclic_cell_massadvy_gpu
   ! ============================================================
   subroutine cyclic_cell_intpx_r2f_gpu(out, inp, wndmdf, &
                                        imrlst, jm, lons, lats, levs, nv, &
                                        async)
      !
      ! do  mass conserving interpolation from different grid at given latitude
      !
      ! author: hann-ming henry juang 2008
      !
      use const, only: RTYPE

      implicit none
      !
      real(kind=RTYPE), intent(out):: out(lats*2, lons/2, levs, nv)
      real(kind=RTYPE):: inp(lons, lats, levs, nv), wrk(lons, lats, levs, nv)
      integer, intent(in):: wndmdf(nv)
      integer, intent(in):: imrlst(lats), jm, lons, lats, levs, nv
      integer, intent(in):: async
      ! local
      real(kind=RTYPE):: xpast(lons + 1, lats, levs), &
                         xnext(lons + 1, lats, levs), &
                         dxp, hfdxp

      integer imf, imp
      ! for ppm
      real(kind=RTYPE):: locs(3*lons, lats, levs)
      integer:: kklist(lons + 1, lats, levs), &
                kstr(lats, levs), kend(lats, levs)
      !----
      real(kind=RTYPE):: hh(3*lons, lats, levs), &
                         mass(3*lons, lats, levs, nv), &
                         qmi(lons + 1, lats, levs, nv), &
                         qpi(lons + 1, lats, levs, nv)
      !
      integer i, j, k, n, ind
      integer nlevs
      ! --------------------
      nlevs = levs*nv
      imf = lons
      ! --------------------
      !$acc data async(async) &
      !$acc& present_or_copyin(imflst(:lats), &
      !$acc& imf,nlevs) &
      !$acc& present( &
      !$acc& jm,levs,nv,lons,lats, &
      !$acc& imrlst(:lats), &
      !$acc& inp(lons,lats,levs,nv), &
      !$acc& out(lats*2,lons/2,levs,nv),wndmdf(nv), &
      !$acc& two_pi, dxf, hfdxf &
      !$acc& ) &
      !$acc& create(xnext(:lons+1,:lats,:levs), &
      !$acc& xpast(:lons+1,:lats,:levs), &
      !$acc& wrk(:lons,:lats,:levs,:nv))
      ! ----------------------------------------
      !$acc parallel loop collapse(2) async(async) &
      !$acc& private(imp,dxp,hfdxp,i)
      do k = 1, levs
         do j = 1, jm
            imp = imrlst(j)
            dxp = two_pi/imp
            hfdxp = 0.5*dxp
            do i = 1, imp + 1
               xpast(i, j, k) = (i - 1)*dxp - hfdxp
            end do

            do i = 1, imf + 1
               xnext(i, j, k) = (i - 1)*dxf - hfdxf
            end do
         end do
      end do
      !
      ! mass conserving interpolation from full grid to reduced grid
      call cyclic_cell_ppm_intp_gpu_v2(xpast, inp, xnext, wrk, &
                                       imrlst, imflst, &
                                       jm, lons, lats, 1, nlevs, two_pi, &
                                       async)

      !$acc parallel loop collapse(3) async(async) &
      !$acc& private(i)
      do n = 1, nv
         do k = 1, levs
            do i = 1, imf/2
               do j = 1, jm
                  out(j, i, k, n) = wrk(i, j, k, n)*wndmdf(n)
               end do
               do j = 1, jm
                  out(jm*2 - j + 1, i, k, n) = wrk(i + imf/2, j, k, n)
               end do
            end do
         end do
      end do
      ! ------------------------------------------------------- >>>>
      ! the data of closed to the equator remains unchanged
      !$acc parallel loop collapse(3) async(async) &
      !$acc& private(i)
      do n = 1, nv
         do k = 1, levs
            do j = 0, 1
               do i = 1, imf/2
                  out(jm/2 + j, i, k, n) = inp(i, jm/2 + j, k, n)*wndmdf(n)
               end do
               do i = 1, imf/2
                  out(3*jm/2 - j + 1, i, k, n) = inp(i + imf/2, jm/2 + j, k, n)
               end do
            end do
         end do
      end do
      !$acc end data
      return
   end subroutine cyclic_cell_intpx_r2f_gpu
   ! ============================================================
   subroutine cyclic_cell_intpx_f2r_gpu(out, inp, wndmdf, &
                                        imrlst, jm, lons, lats, levs, nv, &
                                        async)
      !
      ! do  mass conserving interpolation from different grid at given latitude

      !
      ! author: hann-ming henry juang 2008
      !

      use const, only: RTYPE
      implicit none
      !
      real(kind=RTYPE), intent(out):: out(lons, lats, levs, nv)
      real(kind=RTYPE), intent(in)::  inp(lats*2, lons/2, levs, nv)
      integer, intent(in):: wndmdf(nv)
      integer, intent(in):: imrlst(lats), jm, lons, lats, levs, nv
      integer, intent(in):: async
      ! local
      integer imf, imp
      real(kind=RTYPE):: xpast(lons + 1, lats, levs), &
                         xnext(lons + 1, lats, levs), &
                         dxp, hfdxp
      ! for ppm
      real(kind=RTYPE):: locs(3*lons, lats, levs)
      integer:: kklist(lons + 1, lats, levs), &
                kstr(lats, levs), kend(lats, levs)
      !----
      real(kind=RTYPE):: hh(3*lons, lats, levs), &
                         mass(3*lons, lats, levs, nv), &
                         qmi(lons + 1, lats, levs, nv), &
                         qpi(lons + 1, lats, levs, nv)
      !
      real(kind=RTYPE):: wrk(lons, lats, levs, nv)
      integer i, j, k, n, ind
      ! --------------------
      imf = lons
      ! --------------------
      !$acc data async(async) &
      !$acc& present_or_copyin(imf, &
      !$acc& jm,lons,lats,levs,nv) &
      !$acc& present( &
      !$aacc out(:lons,:lats,:levs,:nv), &
      !$acc& inp(:lats*2,:lons/2,:levs,:nv), &
      !$acc& wndmdf(nv), two_pi, dxf, hfdxf,  &
      !$acc& imflst(:lats),imrlst(:lats) &
      !$acc& ) &
      !$acc& create( &
      !$acc& wrk(:lons,:lats,:levs,:nv), &
      !$acc& xnext(:lons+1,:lats,:levs), &
      !$acc& xpast(:lons+1,:lats,:levs))
      ! -----

      ! ----------------------------------------
      !$acc parallel loop collapse(2) async(async)&
      !$acc& private(imp,dxp,hfdxp,i)
      do k = 1, levs
         do j = 1, jm
            imp = imrlst(j)
            dxp = two_pi/imp
            hfdxp = 0.5*dxp
            do i = 1, imp + 1
               xnext(i, j, k) = (i - 1)*dxp - hfdxp
            end do

            do i = 1, imf + 1
               xpast(i, j, k) = (i - 1)*dxf - hfdxf
            end do
         end do
      end do
      ! change the direction from north-south to east-west and convert wind after advy

      !$acc parallel loop collapse(3) async(async)&
      !$acc& private(i)
      do n = 1, nv
         do k = 1, levs
            do i = 1, jm
               do j = 1, imf/2
                  wrk(j, i, k, n) = inp(i, j, k, n)*wndmdf(n)
               end do
               do j = 1, imf/2
                  wrk(j + imf/2, jm - i + 1, k, n) = inp(i + jm, j, k, n)
               end do
            end do
         end do
      end do
      ! mass conserving interpolation from full grid to reduced grid
      call cyclic_cell_ppm_intp_gpu_v2(xpast, wrk, xnext, out, &
                                       imflst, imrlst, &
                                       jm, lons, lats, 1, levs*nv, &
                                       two_pi, async)

      ! the data of closed to the equator remains unchanged
      !$acc parallel loop collapse(2) async(async)
      do n = 1, nv
         do k = 1, levs
            out(1:imf, jm/2:jm/2 + 1, k, n) = wrk(1:imf, jm/2:jm/2 + 1, k, n)
         end do
      end do

      !$acc end data
      return
   end subroutine cyclic_cell_intpx_f2r_gpu
   ! ============================================================
   subroutine cyclic_cell_ppm_intp_gpu(pp, qq, pn, qn, &
                                       imlst_inp, imlst_out, &
                                       jm, lons, lats, levs, nv, sc, &
                                       run_st)
      !
      ! mass conservation in cyclic bc interpolation: interpolate a group
      ! of grid point  coordiante call pp at interface with quantity qq at
      ! cell averaged to a group of new grid point coordinate call pn at
      ! interface with quantity qn at cell average with ppm spline.
      ! in horizontal with mass conservation is under the condition that
      ! variable value at pp(1)= pp(lons+1)=pn(lons+1)
      !
      ! pp    location at interfac point as input
      ! qq    quantity at averaged-cell as input
      ! pn    location at interface of new grid structure as input
      ! qn    quantity at averaged-cell as output
      ! lons  numer of cells for dimension
      ! imlst_inp  numer of cells for input
      ! imlst_out  numer of cells for output
      ! levs  number of vertical layers
      ! mono  monotonicity o:no, 1:yes
      !
      ! author : henry.juang@noaa.gov
      !
      !
      use const, only: RTYPE
      implicit none
      !
      real(kind=RTYPE), intent(in)::  pp(lons + 1, lats*levs)
      real(kind=RTYPE), intent(in)::  qq(lons, lats, levs, nv)
      real(kind=RTYPE) pn(lons + 1, lats*levs)
      real(kind=RTYPE), intent(out):: qn(lons, lats*levs, nv)
      integer, intent(in):: imlst_inp(lats), imlst_out(lats), jm
      integer, intent(in):: lons, nv, lats, levs
      real(kind=RTYPE), intent(in):: sc
      integer, intent(in) :: run_st
      ! locall
      integer kklist(lons + 1, lats*levs)
      integer kstr(lats*levs), kend(lats*levs)

      real(kind=RTYPE) mass(3*lons, lats, levs, nv)
      real(kind=RTYPE) locs(3*lons, lats*levs)
      real(kind=RTYPE) hh(3*lons, lats*levs)
      real(kind=RTYPE) qmi(lons + 1, lats*levs, nv)
      real(kind=RTYPE) qpi(lons + 1, lats*levs, nv)

      integer n, k, j, i, im
      !$acc data present_or_copyin(imlst_inp(:lats),imlst_out(:lats), &
      !$acc& jm,lons,lats,levs,nv) &
      !$acc& create( mass(:3*lons,:lats,:levs,:nv), &
      !$acc& hh(:3*lons,:lats*levs), qmi(:lons+1,:lats*levs,:nv), &
      !$acc& qpi(:lons+1,:lats*levs,:nv), &
      !$acc& locs(:3*lons,:lats*levs), &
      !$acc& kklist(:lons+1,:lats*levs), kstr(:lats*levs), kend(:lats*levs)) &
      !$acc& async(run_st)
      ! ----------------------------------------
      call ppm_indx(locs, kklist, kstr, kend, pn, pp, sc, &
                    imlst_inp, imlst_out, jm, lons, lats, levs, run_st)
      call ppm_poly(qq, hh, qpi, qmi, mass, locs, kklist, kstr, kend, &
                    imlst_inp, imlst_out, jm, lons, lats, levs, nv, run_st)
      call ppm_intp(qn, pn, locs, kklist, mass, hh, qpi, qmi, &
                    imlst_out, jm, lons, lats, levs, nv, run_st)
      !$acc end data
      return
   end subroutine cyclic_cell_ppm_intp_gpu
   ! ============================================================
   subroutine cyclic_cell_ppm_intp_gpu_v2(pp, qq, pn, qn, &
                                          imlst_inp, imlst_out, &
                                          jm, lons, lats, levs, nv, sc, &
                                          async)
      !
      ! mass conservation in cyclic bc interpolation: interpolate a group
      ! of grid point  coordiante call pp at interface with quantity qq at
      ! cell averaged to a group of new grid point coordinate call pn at
      ! interface with quantity qn at cell average with ppm spline.
      ! in horizontal with mass conservation is under the condition that
      ! variable value at pp(1)= pp(lons+1)=pn(lons+1)
      !
      ! pp    location at interfac point as input
      ! qq    quantity at averaged-cell as input
      ! pn    location at interface of new grid structure as input
      ! qn    quantity at averaged-cell as output
      ! lons  numer of cells for dimension
      ! imlst_inp  numer of cells for input
      ! imlst_out  numer of cells for output
      ! levs  number of vertical layers
      ! mono  monotonicity o:no, 1:yes
      !
      ! author : henry.juang@noaa.gov
      !
      !
      use const, only: RTYPE
      implicit none
      !
      real(kind=RTYPE), intent(in):: pp(lons + 1, lats*levs), &
                                     qq(lons, lats, levs, nv)
      real(kind=RTYPE) pn(lons + 1, lats*levs)
      real(kind=RTYPE), intent(out):: qn(lons, lats*levs, nv)
      integer, intent(in):: imlst_inp(lats), imlst_out(lats), jm, &
                            lons, nv, lats, levs, &
                            async
      real(kind=RTYPE), intent(in):: sc
      ! locall
      integer:: kklist(lons + 1, lats*levs), &
                kstr(lats*levs), kend(lats*levs)

      real(kind=RTYPE):: mass(3*lons, lats, levs, nv), &
                         locs(3*lons, lats*levs), &
                         hh(3*lons, lats*levs), &
                         qmi(lons + 1, lats*levs, nv), &
                         qpi(lons + 1, lats*levs, nv)

      integer n, k, j, i, im
      !$acc data async(async) &
      !$acc& present_or_copyin(imlst_inp(:lats),imlst_out(:lats), &
      !$acc& jm,lons,lats,levs,nv) &
      !$acc& create( mass(:3*lons,:lats,:levs,:nv), &
      !$acc& hh(:3*lons,:lats*levs), qmi(:lons+1,:lats*levs,:nv), &
      !$acc& qpi(:lons+1,:lats*levs,:nv), &
      !$acc& locs(:3*lons,:lats*levs), &
      !$acc& kklist(:lons+1,:lats*levs), kstr(:lats*levs), kend(:lats*levs))
      ! ----------------------------------------
      call ppm_indx_intpsr(locs, kklist, kstr, kend, pn, pp, sc, &
                           imlst_inp, imlst_out, jm, lons, lats, levs, &
                           async)
      call ppm_poly(qq, hh, qpi, qmi, mass, locs, kklist, kstr, kend, &
                    imlst_inp, imlst_out, jm, lons, lats, levs, nv, &
                    async)
      call ppm_intp(qn, pn, locs, kklist, mass, hh, qpi, qmi, &
                    imlst_out, jm, lons, lats, levs, nv, &
                    async)
      !$acc end data
      return
   end subroutine cyclic_cell_ppm_intp_gpu_v2
   ! ------------------------------------------------------------------------
   subroutine vertical_cell_advect_gpu(lons, londim, levs, nvars, &
                                       deltim, ssi, wwi, qql, mass, forward, &
                                       async_id)
      !
      use const, only: RTYPE
      use grid, only: latpart
      use index, only: jlistnum, jlist1, nxjp, nxp, nxjp_acc
      use rank
      implicit none

      real(kind=RTYPE), intent(inout):: qql(levs, nvars, londim)
      integer, intent(in):: londim, levs, nvars, lons, mass, &
                            async_id

      real(kind=RTYPE), intent(in):: deltim, &
                                     ssi(levs + 1, londim), &
                                     wwi(nxp, levs + 1, latpart)

      real(kind=RTYPE):: xgrid(levs + 1, londim), &
                         xpast(levs + 1, londim), &
                         xnext(levs + 1, londim), &
                         dd(levs + 1, londim), &
                         ds(levs, londim), &
                         step(10), dd_step, &
                         dsfact, sstmp, dpdt, check, &
                         da(levs, nvars, londim)

      integer km, i, j, k, n, nst, nstep, nxj, lat, istr
      logical forward
      integer istat
      integer(kind=cuda_stream_kind) :: stream

      stream = acc_get_cuda_stream(async_id)
      !$acc data async(async_id) &
      !$acc& create(dd, ds, xgrid, xpast, xnext, da, &
      !$acc& step, nstep)

      !$acc host_data use_device(dd)
      istat = cudaMemsetAsync(dd, real(0.0, RTYPE), size(dd), stream)
      !$acc end host_data

      !$acc kernels async(async_id)
      do i = 1, londim
         do k = 1, levs
            ds(k, i) = ssi(k, i) - ssi(k + 1, i)
         end do
      end do
      !$acc end kernels

      !$acc parallel loop collapse(3) async(async_id)&
      !$acc& private(lat,nxj,istr)
      do j = 1, jlistnum
         do k = 2, levs
            do i = 1, nxp
               lat = jlist1(j)
               nxj = nxjp(lat)
               if (i .le. nxj) then
                  istr = nxjp_acc(j) - 1
                  dd(k, istr + i) = wwi(i, k, j)*deltim
               end if
            end do
         end do
      end do
      !$acc end parallel loop

      ! <<<< ========================================
      ! !$acc update self(dd, ds) async(async_id)
      ! !$acc wait(async_id)
      ! do i = 1, londim
      !    call def_cfl_step(levs + 1, dd(1, i), ds(1, i), &
      !                      step, nstep, levs + 1, 'advv')
      ! end do
      ! !$acc update device(step) async(async_id)
      ! !$acc wait(async_id)
      ! =============================================
      call def_cfl_step_gpu_type2(step, nstep, dd, ds, 'advv', &
                                  levs + 1, londim, async_id)
      ! ======================================== >>>>
      do nst = 1, nstep
         !hmhj give direction for value larger with k larger
         if (forward) then

            !$acc kernels async(async_id)
            do i = 1, londim
               do k = 1, levs + 1
                  dd_step = dd(k, i)*step(nst)
                  ! for ppm interpolation
                  xgrid(k, i) = ssi(k, i)
                  xpast(k, i) = ssi(k, i)
                  xnext(k, i) = ssi(k, i) + dd_step
               end do
            end do
            do i = 1, londim
               do n = 1, nvars
                  do k = 1, levs
                     da(k, n, i) = qql(k, n, i)
                  end do
               end do
            end do
            !$acc end kernels
         else
            !$acc kernels async(async_id)
            do i = 1, londim
               do k = 1, levs + 1
                  dd_step = dd(k, i)*step(nst)
                  ! for ppm interpolation
                  xgrid(k, i) = ssi(k, i)
                  xpast(k, i) = ssi(k, i) - dd_step
                  xnext(k, i) = ssi(k, i) + dd_step
               end do
            end do
            !$acc end kernels

            call vertical_cell_ppm_intp_gpu(xgrid, qql, xpast, da, &
                                            levs, nvars, londim, async_id)
         end if

         !

         if (mass .eq. 1) then
            !$acc kernels async(async_id)
            do i = 1, londim
               do n = 1, nvars
                  do k = 1, levs
                     dsfact = (xpast(k, i) - xpast(k + 1, i)) &
                              /(xnext(k, i) - xnext(k + 1, i))
                     da(k, n, i) = da(k, n, i)*dsfact
                  end do
               end do
            end do
            !$acc end kernels
         end if

         call vertical_cell_ppm_intp_gpu(xnext, da, xgrid, qql, &
                                         levs, nvars, londim, async_id)
      end do
      !$acc end data

      return
   end subroutine vertical_cell_advect_gpu
   ! ============================================================
   subroutine vertical_cell_ppm_intp_gpu(pp, qq, pn, qn, &
                                         levs, nvars, nxy, async_id)
      !
      ! mass conservation in vertical interpolation: interpolate a group
      ! of grid point  coordiante call pp at interface with quantity qq at
      ! cell averaged to a group of new grid point coordinate call pn at
      ! interface with quantity qn at cell average with ppm spline.
      ! in vertical with mass conservation is under the condition that
      ! pp(1)=pn(1), pp(levs+1)=pn(levs+1)
      !
      ! pp    pressure at interfac level as input
      ! qq    quantity at layer as input
      ! pn    pressure at interface of new grid structure as input
      ! qn    quantity at layer as output
      ! levs  numer of verical layers
      !
      ! author : henry.juang@noaa.gov
      !
      use const, only: RTYPE
      use index, only: nxp
      implicit none
      !
      real(kind=RTYPE) pp(levs + 1, nxy)
      real(kind=RTYPE) qq(levs, nvars, nxy)
      real(kind=RTYPE) pn(levs + 1, nxy)
      real(kind=RTYPE) qn(levs, nvars, nxy)
      integer levs, nvars, nxy
      integer, intent(in):: async_id
      !
      real(kind=RTYPE):: massm, massc, massp, massbot, masstop, &
                         qmi(levs, nvars, nxy), qpi(levs, nvars, nxy), &
                         dql, dqh, dqlist(levs + 1, nvars, nxy), &
                         hh(levs, nxy), &
                         dqi, dqimax, dqimin, dqmono(levs, nvars, nxy), &
                         tl, tl2, tl3, tlp, tlm, tlc, &
                         th(levs + 1, nxy), th2, th3, thp, &
                         thm, thc, &
                         dpp, dqq, c1, c2, qq_t, qmi_t, qpi_t, &
                         hhm, hhc, hhp, dqmono_pre, dqmono_cur

      integer i, k, kl, kh, kk, kkl, kkh, n
      integer kklist(levs + 1, nxy), left, right, mid
      integer, parameter :: mono = 1
      integer istat
      integer(kind=cuda_stream_kind) :: stream

      stream = acc_get_cuda_stream(async_id)

      !$acc data async(async_id) &
      !$acc create(hh,kklist,qmi,qpi,th,dqlist)

      istat = 0
      !$acc parallel loop async(async_id) &
      !$acc& copy(istat)
      do i = 1, nxy
         if ((pp(1, i) .ne. pn(1, i)) .or. &
             (pp(levs + 1, i) .ne. pn(levs + 1, i))) then
            ! print *, ' Error in vertical_cell_ppm_intp for domain values '
            ! print *, "i pp1 pn1", i, pp(1, i), pn(1, i)
            ! print *, "i ppt pnt", i, pp(levs + 1, i), pn(levs + 1, i)
            !$acc atomic
            istat = istat + 1
         end if
      end do
      !$acc end parallel
      !$acc wait(async_id)
      if (istat .ne. 0) then
         print *, "istat=", istat
         call exit(2)
      end if
      ! if (pp(1) .ne. pn(1) .or. pp(levs + 1) .ne. pn(levs + 1)) then
      !    print *, ' Error in vertical_cell_ppm_intp for domain values '
      !    !! print *, ' i pp1 pn1 ppt pnt ', i, &
      !    ! pp(1), pn(1), pp(levs + 1), pn(levs + 1)
      !    call abort
      ! end if

      !$acc host_data use_device(th,dqlist,kklist)
      istat = cudaMemsetAsync(th, real(0.0, RTYPE), size(th), stream)
      istat = cudaMemsetAsync(dqlist, real(0.0, RTYPE), size(dqlist), stream)
      istat = cudaMemsetAsync(kklist, 1, size(kklist), stream)
      !$acc end host_data
      ! ****************************************
      ! prepare thickness for grid
      ! ****************************************
      !$acc kernels async(async_id)
      do i = 1, nxy
         do k = 1, levs
            hh(k, i) = pp(k + 1, i) - pp(k, i)
         end do
      end do
      !$acc end kernels

      ! ****************************************
      ! find kkh
      ! ****************************************
      !$acc parallel loop collapse(2) async(async_id) &
      !$acc& private(left,right,mid)
      do i = 1, nxy
         do k = 1, levs
            left = 1
            right = levs + 1
            do while (right - left > 1)
               mid = (left + right)/2
               if (pn(k + 1, i) .ge. pp(mid, i)) then
                  right = mid
               else
                  left = mid
               end if
            end do
            kklist(k + 1, i) = left
         end do
      end do

      ! ************************************************************
      ! prepare location with monotonic concerns
      ! ************************************************************
      !$acc parallel loop collapse(2) async(async_id)&
      !$acc& create(dqmono)
      do i = 1, nxy
         do n = 1, nvars
            !$acc loop vector private(massbot, masstop, massm, massc, massp, dqi, dqimax, dqimin, hhm, hhc, hhp, dqmono_cur, dqmono_pre)
            do k = 1, levs
               hhc = hh(k, i)
               massc = qq(k, n, i)
               if (k .eq. 1) then
                  hhp = hh(k + 1, i)
                  massp = qq(k + 1, n, i)

                  massbot = (3.*hhc + hhp)*massc &
                            - 2.*hhc*massp
                  massm = massbot/(hhc + hhp)
               else if (k .eq. levs) then
                  hhm = hh(k - 1, i)
                  massm = qq(k - 1, n, i)

                  masstop = (3.*hhc + hhm)*massc &
                            - 2.*hhc*massm
                  massp = masstop/(hhc + hhm)
               else
                  massm = qq(k - 1, n, i)
                  massp = qq(k + 1, n, i)
               end if
               dqi = 0.25*(massp - massm)
               dqimax = max(massm, massc, massp) - massc
               dqimin = massc - min(massm, massc, massp)
               dqmono(k, n, i) = sign(min(abs(dqi), dqimin, dqimax), dqi)
            end do

            ! ************************************************************
            ! compute value at interface with momotone
            ! ************************************************************
            !$acc loop vector private(massc, massm, hhc, hhm, dqmono_cur, dqmono_pre)
            do k = 1, levs
               massc = qq(k, n, i)
               if (k .eq. 1) then
                  qmi(k, n, i) = massc
               else
                  massm = qq(k - 1, n, i)
                  hhc = hh(k, i)
                  hhm = hh(k - 1, i)
                  dqmono_cur = dqmono(k, n, i)
                  dqmono_pre = dqmono(k - 1, n, i)
                  qmi(k, n, i) = (massm*hhc + massc*hhm)/(hhc + hhm) &
                                 + (dqmono_pre - dqmono_cur)/3.0
                  qpi(k - 1, n, i) = qmi(k, n, i)
               end if

               if (k .eq. levs) then
                  qmi(k, n, i) = massc
                  qpi(k, n, i) = massc
               end if
               if (k .eq. 2) then
                  qpi(k - 1, n, i) = massm
               end if
            end do
         end do
      end do
      !$acc end parallel loop

      ! ****************************************
      ! do monotonicity
      ! ****************************************
      if (mono .eq. 1) then
         !$acc parallel loop collapse(3) async(async_id)&
         !$acc private(c1, c2, qq_t, qmi_t, qpi_t)
         do i = 1, nxy
            do n = 1, nvars
               do k = 1, levs
                  qq_t = qq(k, n, i)
                  qmi_t = qmi(k, n, i)
                  qpi_t = qpi(k, n, i)

                  c1 = qpi_t - qq_t
                  c2 = qq_t - qmi_t
                  if (c1*c2 .le. 0.0) then
                     qmi_t = qq_t
                     qpi_t = qq_t
                  end if

                  c1 = (qpi_t - qmi_t) &
                       *(qq_t - 0.5*(qpi_t + qmi_t))
                  c2 = (qpi_t - qmi_t)*(qpi_t - qmi_t)/6.
                  if (c1 .gt. c2) then
                     qmi_t = 3.*qq_t - 2.*qpi_t
                  else if (c1 .lt. -c2) then
                     qpi_t = 3.*qq_t - 2.*qmi_t
                  end if
                  qmi(k, n, i) = qmi_t
                  qpi(k, n, i) = qpi_t
               end do
            end do
         end do
         !$acc end parallel loop
      end if

      ! ************************************************************
      ! start interpolation by integral of ppm spline
      ! ************************************************************
      !$acc parallel loop async(async_id) &
      !$acc& private(th2,th3,kkh,thp,thm,thc)
      do i = 1, nxy
         !$acc loop vector
         do kh = 2, levs + 1
            kkh = kklist(kh, i)
            th(kh, i) = (pn(kh, i) - pp(kkh, i))/hh(kkh, i)
         end do

         !$acc loop vector collapse(2)
         do n = 1, nvars
            do kh = 2, levs + 1
               kkh = kklist(kh, i)

               th2 = th(kh, i)*th(kh, i)
               th3 = th2*th(kh, i)
               thp = th3 - th2
               thm = th3 - 2.*th2 + th(kh, i)
               thc = -2.*th3 + 3.*th2

               dqlist(kh, n, i) = thp*qpi(kkh, n, i) &
                                  + thm*qmi(kkh, n, i) &
                                  + thc*qq(kkh, n, i)
            end do
         end do
      end do
      !$acc end parallel loop

      !$acc parallel loop collapse(3) gang async(async_id) &
      !$acc& private(k,kl,kh,kkh,kkl,tl,dqh,dql,dpp,dqq,kk)
      do i = 1, nxy
         do n = 1, nvars
            do k = 1, levs
               kl = k
               kh = k + 1
               kkh = kklist(kh, i)
               kkl = kklist(k, i)

               tl = th(k, i)
               dqh = dqlist(kh, n, i)
               dql = dqlist(kl, n, i)
               ! ****************************************
               ! mass interpolate
               ! ****************************************
               if (kkh .eq. kkl) then
                  qn(k, n, i) = (dqh - dql)/(th(kh, i) - tl)
               else if (kkh .gt. kkl) then
                  dpp = (1.-tl)*hh(kkl, i) + th(kh, i)*hh(kkh, i)
                  !$acc loop seq
                  do kk = kkl + 1, kkh - 1
                     dpp = dpp + hh(kk, i)
                  end do

                  dql = qq(kkl, n, i) - dql
                  dqq = dql*hh(kkl, i) + dqh*hh(kkh, i)
                  !$acc loop seq
                  do kk = kkl + 1, kkh - 1
                     dqq = dqq + qq(kk, n, i)*hh(kk, i)
                  end do
                  qn(k, n, i) = dqq/dpp
                  ! else
                  !    print *, ' Error in vertical_cell_ppm_intp for lev messed up '
                  !    print *, ' pn ', (pn(kk, i), kk=1, levs + 1)
                  !    print *, ' pp ', (pp(kk, i), kk=1, levs + 1)
               end if
            end do     ! end of k loop
         end do
      end do
      !$acc end parallel loop
      !$acc end data
      return
   end subroutine vertical_cell_ppm_intp_gpu
   ! ============================================================
   subroutine ppm_indx(locs, kklist, kstr, kend, pn, pp, sc, &
                       im_inp, im_out, jm, lons, lats, levs, async)
      ! pp    location at interfac point as input
      ! pn    location at interface of new grid structure as input
      ! lons  numer of cells for dimension
      ! lonp  numer of cells for input
      ! lonn  numer of cells for output
      ! levs  number of vertical layers

      use const, only: RTYPE
      implicit none
      real(kind=RTYPE), intent(out):: locs(3*lons, lats, levs)
      integer, intent(out):: kklist(lons + 1, lats, levs), &
                             kstr(lats, levs), kend(lats, levs)
      real(kind=RTYPE), intent(inout):: pn(lons + 1, lats, levs)
      real(kind=RTYPE), intent(in)   :: pp(lons + 1, lats, levs)
      real(kind=RTYPE), intent(in)   :: sc
      integer, intent(in):: im_inp(lats), im_out(lats), jm
      integer, intent(in):: lons, lats, levs
      integer, intent(in):: async
      !
      integerlonp, lonn
      real(kind=RTYPE) pnmin, pnmax, locbndmin, locbndmax
      integer kk, kh, kkl, left, right, mid
      !
      integer i, j, k
      ! --------------------
      !$acc data copyin(im_inp(:lats),im_out(:lats), jm,levs,sc) &
      !$acc& present( &
      !$acc& locs(3*lons,lats,levs), &
      !$acc& pn(lons+1,lats,levs), &
      !$acc& pp(lons+1,lats,levs), &
      !$acc& kstr(lats,levs), &
      !$acc& kend(lats,levs), &
      !$acc& kklist(lons+1,lats,levs) &
      !$acc& ) async(async)
      ! ----------------------------------------
      !$acc parallel loop collapse(2) async(async) &
      !$acc& private(lonp,lonn,pnmin,pnmax,locbndmin,locbndmax)
      do k = 1, levs
         do j = 1, jm
            lonp = im_inp(j)
            lonn = im_out(j)

            locs(lonp + 1:2*lonp, j, k) = pp(1:lonp, j, k)
            do i = 1, lonp
               locs(i, j, k) = locs(i + lonp, j, k) - sc
               locs(i + 2*lonp, j, k) = locs(i + lonp, j, k) + sc
            end do

            pnmin = pn(1, j, k)
            pnmax = pn(lonn + 1, j, k)

            locbndmin = locs(lonp + 4, j, k)
            locbndmax = locs(2*lonp - 4, j, k)
            if (pnmin .lt. locbndmin - sc) then
               do i = 1, lonn + 1
                  pn(i, j, k) = pn(i, j, k) + int((locbndmin - pnmin)/sc)*sc
               end do
            else if (pnmin .gt. locbndmax) then
               do i = 1, lonn + 1
                  pn(i, j, k) = pn(i, j, k) + (int((locbndmax - pnmin)/sc) - 1)*sc
               end do
            end if
         end do
      end do
      ! ----------------------------------------
      !$acc parallel loop collapse(2) async(async)&
      !$acc& present(im_inp(lats),im_out(lats), jm,levs, &
      !$acc& kstr(lats,levs), locs(3*lons,lats,levs), pn(lons+1,lats,levs)) &
      !$acc& private(i,lonp,lonn,pnmin,pnmax)
      do k = 1, levs
         do j = 1, jm
            lonp = im_inp(j)
            lonn = im_out(j)

            pnmin = pn(1, j, k)
            pnmax = pn(lonn + 1, j, k)
            kstr(j, k) = 0
            if (pnmin .lt. locs(lonp + 1, j, k)) then
               !$acc loop independent
               do i = lonp, 1, -1
                  if ((pnmin .ge. locs(i, j, k)) .and. &
                      (pnmin .lt. locs(i + 1, j, k))) then
                     kstr(j, k) = i
                  end if
               end do
            else
               !$acc loop independent
               do i = lonp + 1, 2*lonp
                  if ((pnmin .ge. locs(i, j, k)) .and. &
                      (pnmin .lt. locs(i + 1, j, k))) then
                     kstr(j, k) = i
                  end if
               end do
            end if
            ! if (kstr(j, k) .eq. 0) then
            !    print *, ' Error: can not find kstr: pnmin locs(1,j,k) locs(2*lonp,j,k) ', &
            !       pnmin, locs(1, j, k), locs(2*lonp, j, k)
            !    print *, ' Error: pn(1) pn(2) pn(3) ', pn(1, j, k), pn(2, j, k), pn(3, j, k)
            ! end if
            kstr(j, k) = max(3, kstr(j, k))
         end do
      end do
      ! ----------------------------------------
      !$acc parallel loop collapse(2) async(async) &
      !$acc& private(i,lonp,lonn,pnmin,pnmax)
      do k = 1, levs
         do j = 1, jm
            lonp = im_inp(j)
            lonn = im_out(j)

            pnmin = pn(1, j, k)
            pnmax = pn(lonn + 1, j, k)
            kend(j, k) = 4*lonp
            if (pnmax .lt. locs(2*lonp + 1, j, k)) then
               !$acc loop independent
               do i = 2*lonp, lonp, -1
                  if (pnmax .ge. locs(i, j, k) .and. pnmax .lt. locs(i + 1, j, k)) then
                     kend(j, k) = i + 1
                  end if
               end do
            else
               !$acc loop independent
               do i = 2*lonp + 1, 3*lonp - 1
                  if (pnmax .ge. locs(i, j, k) .and. pnmax .lt. locs(i + 1, j, k)) then
                     kend(j, k) = i + 1
                  end if
               end do
            end if
            ! if(kend(j,k).eq.4*lonp) then
            !    print *,' Error: cannot get kend: pnmax locs(lonp,j,k) locs(3*lonp,j,k)',&
            !         kend(j,k), pnmax,locs(lonp,j,k),locs(3*lonp,j,k)
            !    print *,' Error: pn(lonn-1) pn(lonn) pn(lonn+1) ',               &
            !         pn(lonn-1,j,k),pn(lonn,j,k),pn(lonn+1,j,k)
            ! end if
            kend(j, k) = min(3*lonp - 2, kend(j, k))
         end do
      end do
    !! ****************************************
    !! << kklist >>
    !! ****************************************
      !$acc parallel loop collapse(2) async(async) &
      !$acc& private(i,lonn,mid,right,left)
      do k = 1, levs
         do j = 1, jm
            lonn = im_out(j)
            ! start interpolation by integral of ppm spline
            kklist(1, j, k) = kstr(j, k)
            !$acc loop private(mid,left,right)
            do i = 1, lonn
               left = kstr(j, k)
               right = kend(j, k) + 1
               do while (right - left > 1)
                  mid = (left + right)/2
                  if (pn(i + 1, j, k) .lt. locs(mid, j, k)) then
                     right = mid
                  else
                     left = mid
                  end if
               end do
               kklist(i + 1, j, k) = left
            end do
         end do
      end do
      !$acc end data
      return
   end subroutine ppm_indx
   ! ============================================================
   subroutine ppm_indx_intpsr(locs, kklist, kstr, kend, pn, pp, sc, &
                              im_inp, im_out, jm, lons, lats, levs, async)
      ! pp    location at interfac point as input
      ! pn    location at interface of new grid structure as input
      ! lons  numer of cells for dimension
      ! lonp  numer of cells for input
      ! lonn  numer of cells for output
      ! levs  number of vertical layers
      use const, only: RTYPE
      implicit none
      real(kind=RTYPE), intent(out):: locs(3*lons, lats, levs)
      integer, intent(out):: kklist(lons + 1, lats, levs), &
                             kstr(lats, levs), kend(lats, levs)
      integer :: kklist_test(lons + 1, lats, levs), k1, k2
      real(kind=RTYPE) tmp(lons + 1, lats, levs)
      real(kind=RTYPE), intent(inout):: pn(lons + 1, lats, levs)
      real(kind=RTYPE), intent(in)   :: pp(lons + 1, lats, levs)
      real(kind=RTYPE), intent(in)   :: sc
      integer, intent(in):: im_inp(lats), im_out(lats), jm
      integer, intent(in):: lons, lats, levs
      integer, intent(in):: async
      !
      integerlonp, lonn
      real(kind=RTYPE) pnmin, pnmax, locbndmin, locbndmax
      integer kk, kh, kkl, left, right, mid
      !
      integer i, j, k
      ! --------------------
      !$acc data async(async) &
      !$acc& present_or_copyin(im_inp(:lats),im_out(:lats), jm,levs,sc) &
      !$acc& present( &
      !$acc& locs(3*lons,lats,levs), &
      !$acc& pn(lons+1,lats,levs), &
      !$acc& pp(lons+1,lats,levs), &
      !$acc& kstr(lats,levs), &
      !$acc& kend(lats,levs), &
      !$acc& kklist(lons+1,lats,levs) &
      !$acc& )
      ! ----------------------------------------
      !$acc parallel loop collapse(2) async(async)&
      !$acc& private(lonp,lonn,pnmin,pnmax,locbndmin,locbndmax)
      do k = 1, levs
         do j = 1, jm
            lonp = im_inp(j)
            lonn = im_out(j)

            locs(lonp + 1:2*lonp, j, k) = pp(1:lonp, j, k)
            do i = 1, lonp
               locs(i, j, k) = locs(i + lonp, j, k) - sc
               locs(i + 2*lonp, j, k) = locs(i + lonp, j, k) + sc
            end do

            pnmin = pn(1, j, k)
            pnmax = pn(lonn + 1, j, k)

            locbndmin = locs(lonp + 4, j, k)
            locbndmax = locs(2*lonp - 4, j, k)
            if (pnmin .lt. locbndmin - sc) then
               do i = 1, lonn + 1
                  pn(i, j, k) = pn(i, j, k) + int((locbndmin - pnmin)/sc)*sc
               end do
            else if (pnmin .gt. locbndmax) then
               do i = 1, lonn + 1
                  pn(i, j, k) = pn(i, j, k) + (int((locbndmax - pnmin)/sc) - 1)*sc
               end do
            end if
         end do
      end do
      !$acc parallel loop collapse(2) async(async)&
      !$acc& private(i,lonp,lonn,pnmax)
      do k = 1, levs
         do j = 1, jm
            lonp = im_inp(j)
            lonn = im_out(j)

            kstr(j, k) = floor((pn(1, j, k) - locs(1, j, k))/(sc/lonp)) + 1
            kstr(j, k) = max(3, kstr(j, k))

            pnmax = pn(lonn + 1, j, k)
            kend(j, k) = floor((pnmax - locs(1, j, k))/(sc/lonp)) + 2
            kend(j, k) = min(3*lonp - 2, kend(j, k))
         end do
      end do
    !! ****************************************
    !! << kklist >>
    !! ****************************************
      !$acc parallel loop collapse(2) async(async)&
      !$acc& private(i,kk,lonp,lonn)
      do k = 1, levs
         do j = 1, jm
            lonp = im_inp(j)
            lonn = im_out(j)

            kklist(1, j, k) = kstr(j, k)
            !$acc loop
            do i = 2, lonn + 1
               kk = floor((pn(i, j, k) - locs(1, j, k))/(sc/lonp)) + 1
               if (abs(pn(i, j, k) - locs(kk + 1, j, k)) .le. 1e-10) then
                  kk = kk + 1
               end if
               if (pn(i, j, k) .lt. locs(kk, j, k)) then
                  kk = kk - 1
               end if
               kklist(i, j, k) = min(kk, kend(j, k) + 1)
            end do
         end do
      end do
      !$acc end data

      return
   end subroutine ppm_indx_intpsr
   ! ============================================================
   subroutine ppm_poly(qq, hh, qpi, qmi, &
                       mass, locs, kklist, kstr, kend, &
                       im_inp, im_out, jm, lons, lats, levs, nv, async)
      use const, only: RTYPE
      implicit none
      !
      real(kind=RTYPE) qq(lons, lats, levs, nv)
      real(kind=RTYPE) mass(3*lons, lats, levs, nv)
      real(kind=RTYPE) locs(3*lons, lats, levs)
      real(kind=RTYPE) hh(3*lons, lats, levs)

      real(kind=RTYPE) fm(3*lons, lats, levs)
      real(kind=RTYPE) fn(3*lons, lats, levs)

      real(kind=RTYPE) qmi(lons + 1, lats, levs, nv)
      real(kind=RTYPE) qpi(lons + 1, lats, levs, nv)
      integer kklist(lons + 1, lats, levs)
      integer kstr(lats, levs), kend(lats, levs)

      !
      integer im_inp(lats), im_out(lats), jm
      integer lons, lonn, lonp, lats, levs, nv
      integer async
      integer, parameter :: mono = 1
      !
      real(kind=RTYPE) dqmono(3*lons)
      real(kind=RTYPE) dqi, dqimax, dqimin
      real(kind=RTYPE) c1, c2, cc
      real(kind=RTYPE), parameter:: r3 = 1./3.
      real(kind=RTYPE), parameter:: r6 = 1./6.
      integer i, j, k, n, kkl
      !
      !$acc data async(async)&
      !$acc& present_or_copyin(im_out(:lats),jm,levs,nv) &
      !$acc& present( &
      !$acc& qq(lons,lats,levs,nv), &
      !$acc& kstr(lats,levs), kend(lats,levs), &
      !$acc& hh(3*lons,lats,levs), locs(3*lons,lats,levs), &
      !$acc& fm(3*lons,lats,levs), fn(3*lons,lats,levs), &
      !$acc& kklist(lons+1,lats,levs), mass(3*lons,lats,levs,nv), &
      !$acc& qmi(lons+1,lats,levs,nv), qpi(lons+1,lats,levs,nv) &
      !$acc& ) &
      !$acc& create(fm(:3*lons,:lats,:levs),fn(:3*lons,:lats,:levs))
      ! ----------------------------------------
      !$acc parallel loop collapse(2) async(async)
      do k = 1, levs
         do j = 1, jm
            !$acc loop independent
            do i = kstr(j, k) - 2, kend(j, k) + 2
               hh(i, j, k) = locs(i + 1, j, k) - locs(i, j, k)
            end do
            !$acc loop independent
            do i = kstr(j, k) - 1, kend(j, k) + 2
               cc = 1./(hh(i, j, k) + hh(i - 1, j, k))
               fm(i, j, k) = hh(i, j, k)*cc
               fn(i, j, k) = hh(i - 1, j, k)*cc
            end do
         end do
      end do
      ! ----------------------------------------
      !$acc parallel loop collapse(3) async(async) &
      !$acc& private(i,lonn,kkl,dqi,dqimax,dqimin, dqmono(:3*lons))
      do n = 1, nv
         do k = 1, levs
            do j = 1, jm
               lonp = im_inp(j)
               lonn = im_out(j)

               mass(1:lonp, j, k, n) = qq(1:lonp, j, k, n)
               mass(1 + lonp:lonp*2, j, k, n) = qq(1:lonp, j, k, n)
               mass(1 + lonp*2:lonp*3, j, k, n) = qq(1:lonp, j, k, n)
               ! prepare location with monotonic concerns
               !
               !$acc loop independent
               do i = kstr(j, k) - 2, kend(j, k) + 2
                  dqi = 0.25*(mass(i + 1, j, k, n) - mass(i - 1, j, k, n))
                  dqimax = max(mass(i - 1, j, k, n), mass(i, j, k, n), mass(i + 1, j, k, n)) - mass(i, j, k, n)
                  dqimin = mass(i, j, k, n) - min(mass(i - 1, j, k, n), mass(i, j, k, n), mass(i + 1, j, k, n))
                  dqmono(i) = sign(min(abs(dqi), dqimin, dqimax), dqi)
               end do
               !
               ! compute value at interface with monotone
               !
               !$acc loop
               do i = 1, lonn + 1
                  kkl = kklist(i, j, k)
                  qmi(i, j, k, n) = mass(kkl - 1, j, k, n)*fm(kkl, j, k) &
                                    + mass(kkl, j, k, n)*fn(kkl, j, k) &
                                    + (dqmono(kkl - 1) - dqmono(kkl))*r3

                  qpi(i, j, k, n) = mass(kkl, j, k, n)*fm(kkl + 1, j, k) &
                                    + mass(kkl + 1, j, k, n)*fn(kkl + 1, j, k) &
                                    + (dqmono(kkl) - dqmono(kkl + 1))*r3
               end do
            end do
         end do
      end do
      !
      ! do monotonicity within cell
      !
      if (mono .eq. 1) then
         !$acc parallel loop collapse(3) async(async)&
         !$acc& private(i,lonn,kkl,c1,c2,cc)
         do n = 1, nv
            do k = 1, levs
               do j = 1, jm
                  lonn = im_out(j)
                  do i = 1, lonn + 1
                     kkl = kklist(i, j, k)

                     c1 = qpi(i, j, k, n) - mass(kkl, j, k, n)
                     c2 = mass(kkl, j, k, n) - qmi(i, j, k, n)
                     if (c1*c2 .le. 0.0) then
                        qmi(i, j, k, n) = mass(kkl, j, k, n)
                        qpi(i, j, k, n) = mass(kkl, j, k, n)
                     else
                        cc = qpi(i, j, k, n) - qmi(i, j, k, n)
                        c1 = cc*(mass(kkl, j, k, n) - 0.5*(qpi(i, j, k, n) + qmi(i, j, k, n)))
                        c2 = cc*cc*r6
                        if (c1 .gt. c2) then
                           qmi(i, j, k, n) = 3.*mass(kkl, j, k, n) - 2.*qpi(i, j, k, n)
                        else if (c1 .lt. -c2) then
                           qpi(i, j, k, n) = 3.*mass(kkl, j, k, n) - 2.*qmi(i, j, k, n)
                        end if
                     end if
                  end do
               end do
            end do
         end do
      end if
      !$acc end data
      return
   end subroutine ppm_poly
   ! ============================================================
   subroutine ppm_intp(out, pn, locs, kklist, &
                       mass, hh, qpi, qmi, &
                       im_out, jm, lons, lats, levs, nv, async)
      use const, only: RTYPE
      implicit none
      real(kind=RTYPE), intent(out):: out(lons, lats, levs, nv)
      real(kind=RTYPE), intent(in) :: pn(lons + 1, lats, levs), &
                                      locs(3*lons, lats, levs), &
                                      hh(3*lons, lats, levs), &
                                      mass(3*lons, lats, levs, nv), &
                                      qmi(lons + 1, lats, levs, nv), &
                                      qpi(lons + 1, lats, levs, nv)
      integer, intent(in):: kklist(lons + 1, lats, levs)
      integer, intent(in):: im_out(lats), jm, lons, lats, levs, nv, async
      ! local
      integer lonn
      real(kind=RTYPE) tl(lons + 1, lats, levs), tl2, tl3
      real(kind=RTYPE) tlp(lons + 1, lats, levs)
      real(kind=RTYPE) tlm(lons + 1, lats, levs)
      real(kind=RTYPE) tlc(lons + 1, lats, levs)
      real(kind=RTYPE) dqlist(lons + 1, lats, levs, nv)

      real(kind=RTYPE) dql, dqh, dqq, dpp
      integer i, j, k, n, kk, kkl, kkh
      !$acc data async(async) &
      !$acc& present_or_copyin(im_out(:lats),jm,lons,lats,levs,nv) &
      !$acc& present( &
      !$acc& kklist(lons+1,lats,levs), pn(lons+1,lats,levs), &
      !$acc& locs(3*lons,lats,levs), hh(3*lons,lats,levs), &
      !$acc& mass(3*lons,lats,levs,nv), &
      !$acc& qpi(lons+1,lats,levs,nv), qmi(lons+1,lats,levs,nv), &
      !$acc& out(lons,lats,levs,nv) &
      !$acc& ) &
      !$acc& create(tl(:lons+1,:lats,:levs),dqlist(:lons+1,:lats,:levs,:nv),&
      !$acc& tlp(:lons+1,:lats,:levs), tlm(:lons+1,:lats,:levs), &
      !$acc& tlc(:lons+1,:lats,:levs))
      ! ----------------------------------------
      !$acc parallel loop collapse(2) async(async) &
      !$acc& private(lonn,i,kkl,tl2,tl3)
      do k = 1, levs
         do j = 1, jm
            lonn = im_out(j)
            do i = 1, lonn + 1
               kkl = kklist(i, j, k)
               tl(i, j, k) = (pn(i, j, k) - locs(kkl, j, k))/hh(kkl, j, k)
               tl2 = tl(i, j, k)*tl(i, j, k)
               tl3 = tl2*tl(i, j, k)

               tlp(i, j, k) = tl3 - tl2
               tlm(i, j, k) = tl3 - 2.*tl2 + tl(i, j, k)
               tlc(i, j, k) = -2.*tl3 + 3.*tl2
            end do
         end do
      end do
      !$acc parallel loop collapse(3) async(async) &
      !$acc& private(lonn,i,kkl)
      do n = 1, nv
         do k = 1, levs
            do j = 1, jm
               lonn = im_out(j)
               do i = 1, lonn + 1
                  kkl = kklist(i, j, k)
                  dqlist(i, j, k, n) = tlp(i, j, k)*qpi(i, j, k, n) &
                                       + tlm(i, j, k)*qmi(i, j, k, n) &
                                       + tlc(i, j, k)*mass(kkl, j, k, n)
               end do
            end do
         end do
      end do

      !$acc parallel loop collapse(3) async(async) &
      !$acc& private(lonn,i,j,k,kkl,kkh,dql,dqh,dpp,dqq,kk)
      do n = 1, nv
         do k = 1, levs
            do j = 1, jm
               lonn = im_out(j)
               !$acc loop
               do i = 1, lonn
                  kkl = kklist(i, j, k)
                  kkh = kklist(i + 1, j, k)

                  dql = dqlist(i, j, k, n)
                  dqh = dqlist(i + 1, j, k, n)

                  if (kkh .eq. kkl) then
                     out(i, j, k, n) = (dqh - dql)/(tl(i + 1, j, k) - tl(i, j, k))
                  else if (kkh .gt. kkl) then
                     dpp = (1.-tl(i, j, k))*hh(kkl, j, k) + tl(i + 1, j, k)*hh(kkh, j, k)
                     !$acc loop seq
                     do kk = kkl + 1, kkh - 1
                        dpp = dpp + hh(kk, j, k)
                     end do
                     dql = mass(kkl, j, k, n) - dql
                     dqq = dql*hh(kkl, j, k) + dqh*hh(kkh, j, k)
                     !$acc loop seq
                     do kk = kkl + 1, kkh - 1
                        dqq = dqq + mass(kk, j, k, n)*hh(kk, j, k)
                     end do
                     out(i, j, k, n) = dqq/dpp
                  else
                     print *, ' Error in cyclic_cell_ppm_intp location messed up '
                     print *, ' levs=', levs, ' j=', j, ' i=', i
                     print *, ' kkl=', kkl, ' kkh=', kkh
                     ! call abort
                  end if
               end do
            end do
         end do
      end do
      !$acc end data
   end subroutine ppm_intp
   ! ============================================================
   subroutine def_cfl_step_gpu(step, nstep, dist, del, job, &
                               im, lons, lats, levs, async)
      !
      ! compute the deformation cfl condition
      ! select the maxima value of the deformation CFL and provide step to
      ! avoid it.
      !
      !
      use const, only: RTYPE
      implicit none
      !

      real(kind=RTYPE), intent(out):: step(10)
      integer, intent(out):: nstep
      real(kind=RTYPE), intent(in) :: dist(lons + 1, lats, levs), del(lons + 1, lats)
      integer, intent(in) :: im(lats), lons, lats, levs
      integer, intent(in) :: async
      ! local
      integer i, j, n, k, nchk, imp
      real(kind=RTYPE) rstep, check, check_max, check_point, loc_max
      real(kind=RTYPE) safe_step, last_step
      character*4 job
      !
      check_point = 1.00
      safe_step = 0.99
      nstep = 1
      step(1) = 1.0

      check_max = 0.0
      loc_max = 0.0
      !$acc data present_or_copyin(im(:lats),lons,lats,levs)
      !----------------------------------------
      !$acc parallel async(async) &
      !$acc& copy(check_max) &
      !$acc& copyin(loc_max, check_point) &
      !$acc& present(dist(lons+1,lats,levs), del(lons+1,lats), &
      !$acc& im(lats),lons,lats,levs )
      !--------------------
      !$acc loop collapse(2) &
      !$acc& private(imp,i,check,loc_max) &
      !$acc& reduction(max:check_max)
      do k = 1, levs
         do j = 1, lats
            imp = im(j)
            loc_max = 0.
            !$acc loop seq
            do i = 1, imp
               check = abs((dist(i + 1, j, k) - dist(i, j, k))/del(i, j))
               if ((check .ge. check_point) .and. (check .gt. loc_max)) then
                  loc_max = check
               end if
            end do
            check_max = max(check_max, loc_max)
         end do
      end do
      !$acc end parallel
      !$acc end data

      if (check_max .ge. check_point) then
         nstep = int(check_max/safe_step) + 1
         ! ----------------------------------------
         write (*, "(1X,A, 1pe15.7, A, i3, 3(1X,A))") &
            " max def_cfl ", check_max, " needs ", nstep, &
            " steps in ", job, "(gpu) processing"
         ! ----------------------------------------
         rstep = safe_step/check_max
         ! rstep = 1./nstep
         do n = 1, nstep - 1
            step(n) = rstep
         end do
         last_step = 1.-(nstep - 1)*rstep
         step(nstep) = last_step

         if (nstep > 10) then
            call exit(1)
            print *, "Error in def_cfl_step_gpu"
         end if
      end if
      !$acc update device(step)

      return
   end subroutine def_cfl_step_gpu
   ! ============================================================
   subroutine def_cfl_step_gpu_type2(step, nstep, dist, del, job, &
                                     lda, len, async_id)
      !
      ! compute the deformation cfl condition
      ! select the maxima value of the deformation CFL and provide step to
      ! avoid it.
      !
      !
      use const, only: RTYPE
      use rank
      implicit none
      !

      real(kind=RTYPE), intent(out):: step(10)
      integer, intent(out):: nstep
      real(kind=RTYPE), intent(in) :: dist(lda, len), del(lda - 1, len)
      integer, intent(in) :: lda, len
      integer, intent(in) :: async_id
      ! local
      integer i, j, n, k, nchk
      real(kind=RTYPE) rstep, check, check_max, check_point, loc_max
      real(kind=RTYPE) safe_step, last_step
      character*4 job
      !
      check_point = 1.00
      safe_step = 0.99
      nstep = 1
      step(1) = 1.0

      check_max = 0.0
      loc_max = 0.0
      !$acc update device(nstep,step(:10)) async(async_id)
      !----------------------------------------
      !$acc kernels async(async_id) &
      !$acc& copy(check_max) &
      !$acc& copyin(loc_max, check_point,safe_step)
      !--------------------
      !$acc loop &
      !$acc& private(i,check,loc_max)
      do j = 1, len
         loc_max = 0.
         !$acc loop seq
         do i = 1, lda - 1
            check = abs((dist(i + 1, j) - dist(i, j))/del(i, j))
            if ((check .ge. check_point) .and. (check .gt. loc_max)) then
               loc_max = check
            end if
         end do
         !$acc atomic
         check_max = max(check_max, loc_max)
      end do
      !$acc end kernels

      if (check_max .ge. check_point) then
         nstep = int(check_max/safe_step) + 1
         rstep = safe_step/check_max
         ! ----------------------------------------
         write (*, "(1X,A, 1pe15.7, 2(1X,A,i3), 3(1X,A))") &
            " max def_cfl ", check_max, "needs ", nstep, &
            "steps in ", myrank, "rank in", job, "(gpu) processing"
         ! ----------------------------------------
         if (nstep > 10) then
            call exit(1)
            print *, "Error in def_cfl_step_gpu_type2"
         end if

         do n = 1, nstep - 1
            step(n) = rstep
         end do
         last_step = 1.-(nstep - 1)*rstep
         step(nstep) = last_step
      end if
      !$acc update device(step) async(async_id)

      return
   end subroutine def_cfl_step_gpu_type2
   ! ============================================================
end module mod_ndslfv_monoadv_gpu
