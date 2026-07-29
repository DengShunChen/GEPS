#define NCCLCHECK(ierr) call nccl_check_helper(ierr, __FILE__, __LINE__)
module mod_stochastic_physics_gpu
   use mpe, only: mpe_bcast, mpe_double, mpe_integer
   use rank, only: myrank
   use index
   use param
   use const, only: aki, bki, first_call, dosppt, doshum, doskeb, dossst, &
                    poly, dpoly, wdfac, wcfac, onocos, radsq, weight, &
                    RTYPE, rad, cosl, doskeb_dc, polyf
   use mersenne_twister_gpu, only: random_setseed, random_gauss, random_stat
   implicit none
   private

   real(kind=RTYPE), public, allocatable, save :: rpattern_sppt_n2du(:, :, :, :), &
                                          rpattern_shum_n2du(:, :, :, :), &
                                          rpattern_skeb_n2du(:, :, :, :), &
                                          rpattern_ssst_n2du(:, :, :, :)
   real(kind=RTYPE), public, allocatable, save :: rpattern_sppt_n2dv(:, :, :, :), &
                                          rpattern_shum_n2dv(:, :, :, :), &
                                          rpattern_skeb_n2dv(:, :, :, :), &
                                          rpattern_ssst_n2dv(:, :, :, :)
   real(kind=RTYPE), public, allocatable, save :: rpattern_sppt_n2d(:, :, :), &
                                          rpattern_shum_n2d(:, :, :), &
                                          rpattern_skeb_n2d(:, :, :), &
                                          rpattern_ssst_n2d(:, :, :)
   real, allocatable, save :: rpattern_sppt_kenorm(:, :, :), &
                              rpattern_shum_kenorm(:, :, :), &
                              rpattern_skeb_kenorm(:, :, :), &
                              rpattern_ssst_kenorm(:, :, :)
   real(kind=RTYPE), allocatable, save :: rpattern_sppt_spec(:, :, :), &
                                          rpattern_shum_spec(:, :, :), &
                                          rpattern_skeb_spec(:, :, :), &
                                          rpattern_ssst_spec(:, :, :)
   real(kind=RTYPE), allocatable, save :: rpattern_sppt_varspec(:, :), &
                                          rpattern_shum_varspec(:, :), &
                                          rpattern_skeb_varspec(:, :), &
                                          rpattern_ssst_varspec(:, :)
   real, allocatable, save :: rpattern_sppt_stdev(:), &
                              rpattern_shum_stdev(:), &
                              rpattern_skeb_stdev(:), &
                              rpattern_ssst_stdev(:) ! stochastic physics tendency amplitude
   real, allocatable, save :: rpattern_sppt_decortau(:), &
                              rpattern_shum_decortau(:), &
                              rpattern_skeb_decortau(:), &
                              rpattern_ssst_decortau(:) ! time scales
   real, allocatable, save :: rpattern_sppt_lenscale(:), &
                              rpattern_shum_lenscale(:), &
                              rpattern_skeb_lenscale(:), &
                              rpattern_ssst_lenscale(:) ! length scales
   real, allocatable, save :: rpattern_sppt_phi(:), &
                              rpattern_shum_phi(:), &
                              rpattern_skeb_phi(:), &
                              rpattern_ssst_phi(:)
   integer, allocatable, save :: rpattern_sppt_mlmax(:), &
                                 rpattern_shum_mlmax(:), &
                                 rpattern_skeb_mlmax(:), &
                                 rpattern_ssst_mlmax(:)
   integer, allocatable, save :: rpattern_sppt_jtrun(:), &
                                 rpattern_shum_jtrun(:), &
                                 rpattern_skeb_jtrun(:), &
                                 rpattern_ssst_jtrun(:)
   integer, save :: rpattern_sppt_mlmax_max, &
                    rpattern_shum_mlmax_max, &
                    rpattern_skeb_mlmax_max, &
                    rpattern_ssst_mlmax_max
   integer, save :: rpattern_sppt_jtrun_max, &
                    rpattern_shum_jtrun_max, &
                    rpattern_skeb_jtrun_max, &
                    rpattern_ssst_jtrun_max
   integer, allocatable, save :: rpattern_sppt_mlsort(:, :, :), &
                                 rpattern_shum_mlsort(:, :, :), &
                                 rpattern_skeb_mlsort(:, :, :), &
                                 rpattern_ssst_mlsort(:, :, :)
   integer, allocatable, save :: rpattern_sppt_msort(:, :), &
                                 rpattern_shum_msort(:, :), &
                                 rpattern_skeb_msort(:, :), &
                                 rpattern_ssst_msort(:, :)
   integer, allocatable, save :: rpattern_sppt_lsort(:, :), &
                                 rpattern_shum_lsort(:, :), &
                                 rpattern_skeb_lsort(:, :), &
                                 rpattern_ssst_lsort(:, :)
   type(random_stat), public, allocatable, save :: rpattern_sppt_rstate(:), &
                                                   rpattern_shum_rstate(:), &
                                                   rpattern_skeb_rstate(:), &
                                                   rpattern_ssst_rstate(:)
   integer, public, allocatable, save :: rpattern_sppt_seed(:), &
                                         rpattern_shum_seed(:), &
                                         rpattern_skeb_seed(:), &
                                         rpattern_ssst_seed(:)

   integer, save :: recnsppt = 1, recnskeb = 1, recnshum = 1
   real ::  dt
   logical, public :: ncep_seeds = .false.
   real, allocatable :: sl(:)


   ! SPPT
   integer :: nsppt = 0
   real(kind=RTYPE), allocatable, save :: sppt3d(:, :, :)
   real :: sppt(5) = -999.             ! amplitude(0.~1.)
   real :: sppt_seed(5) = -999.        ! random seeds
   real :: sppt_decort(5) = -999.      ! time scales(seconds)
   real :: sppt_lscale(5) = -999.      ! length scales(meters)
   real, allocatable, dimension(:) :: vfact_sppt
   real, public :: sppt_sigtop1 = 0.1
   real, public :: sppt_sigtop2 = 0.025
   real, public :: sppt_sigbot1 = 0.975
   real, public :: sppt_sigbot2 = 0.9
   logical, public :: sppt_sfclimit = .false.
   logical, public :: sppt_logit = .true.

   ! SHUM
   integer :: nshum = 0
   real(kind=RTYPE), allocatable, save :: shum3d(:, :, :)
   real :: shum(5) = -999.             ! amplitude(0.~1.)
   real :: shum_seed(5) = -999.         ! random seeds
   real :: shum_decort(5) = -999.      ! time scales(seconds)
   real :: shum_lscale(5) = -999.      ! length scales(meters)
   real, allocatable, dimension(:) :: vfact_shum
   real, public :: shum_sigefold = 0.2
   real(kind=RTYPE), allocatable, save :: shum3d_dq(:, :, :)

   ! SKEB
   integer :: nskeb = 0, skeblevs
   real(kind=RTYPE), allocatable, save :: skeb3du(:, :, :), skeb3dv(:, :, :), diss_est(:, :, :)
   real(kind=RTYPE), allocatable, save :: diss_dc(:, :, :)    !dissipation from deep convection
   real(kind=RTYPE), allocatable, save :: kea(:, :, :), keb(:, :, :)
   real :: skeb(5) = -999.             ! amplitude(0.~1.)
   real :: skeb_seed(5) = -999.        ! random seeds
   real :: skeb_decort(5) = -999.      ! time scales(seconds)
   real :: skeb_lscale(5) = -999.      ! length scales(meters)
   real, allocatable, dimension(:) :: vfact_skeb
   real, public :: skeb_sigtop1 = 0.1
   real, public :: skeb_sigtop2 = 0.025
   real, public :: skeb_sigbot1 = 0.975
   real, public :: skeb_sigbot2 = 0.9
   real, public :: skebnorm = 2
   real, public :: skeb_vdof = 5 ! proxy for vertical correlation, 5 is close to 40 passes of the 1-2-1 filter in the GFS
   real, public :: skebfilt = 12
   real, public, allocatable, dimension(:, :) :: skeb_vwts
   integer, public, allocatable, dimension(:, :) :: skeb_vpts

   ! SSST
   integer :: nssst = 0
   real(kind=RTYPE), allocatable, save :: ssst3d(:, :, :)
   real :: ssst(5) = -999.             ! amplitude(0.~1.)
   real :: ssst_seed(5) = -999.         ! random seeds
   real :: ssst_decort(5) = -999.      ! time scales(seconds)
   real :: ssst_lscale(5) = -999.      ! length scales(meters)
   real, allocatable, dimension(:) :: vfact_ssst

   logical :: do_unit_test = .false.


   public nsppt, sppt, sppt_seed, sppt_decort, sppt_lscale, sppt3d
   public nshum, shum, shum_seed, shum_decort, shum_lscale, shum3d, shum3d_dq
   public nskeb, skeb, skeb_seed, skeb_decort, skeb_lscale, &
      skeb3du, skeb3dv, diss_est, diss_dc, &
      skeblevs, keb, kea
   public nssst, ssst, ssst_seed, ssst_decort, ssst_lscale, ssst3d

   public init_stochastic_physics_gpu, &
      run_stochastic_physics_gpu, &
      destroy_stochastic_physics_gpu

   public spptout, shumout, skebout, skebest_gpu
   public avevar_sppt2d, do_unit_test

contains

   subroutine init_stochastic_physics_gpu(dtau)
      implicit none
      integer :: n, k
      integer :: async_id = 1
      real :: dtau
      allocate (sl(lev))
      ! calculation sigma values
      do k = 1, lev
         sl(k) = 0.5*(aki(k)/1013.0 + bki(k) + aki(k + 1)/1013.0 + bki(k + 1))
      end do
      if (dosppt) then
         call init_sppt(dtau)
         !$acc enter data copyin(sppt3d) copyin(vfact_sppt, rpattern_sppt_n2d, &
         !$acc&      rpattern_sppt_spec, rpattern_sppt_varspec, rpattern_sppt_mlsort) async(async_id)
         !$acc wait(async_id)
      end if
      if (doshum) then
         call init_shum(dtau)
         !$acc enter data create(shum3d) copyin(shum3d_dq, vfact_shum, rpattern_shum_n2d, &
         !$acc&      rpattern_shum_spec, rpattern_shum_varspec, rpattern_shum_mlsort) async(async_id)
         !$acc wait(async_id)
      end if
      if (doskeb) then
         call init_skeb(dtau)
         !$acc enter data create(skeb3du, skeb3dv) copyin(diss_est, diss_dc, &
         !$acc&      keb, kea, vfact_skeb, skeb_vwts, skeb_vpts, rpattern_skeb_n2du, &
         !$acc&      rpattern_skeb_n2dv, rpattern_skeb_spec, rpattern_skeb_varspec, rpattern_skeb_mlsort) async(async_id)
         !$acc wait(async_id)

      end if
      if (dossst) then
         call init_ssst(dtau)
         !$acc enter data create(ssst3d) copyin(vfact_ssst, rpattern_ssst_spec, &
         !$acc&      rpattern_ssst_varspec, rpattern_ssst_mlsort) async(async_id)
         !$acc wait(async_id)

      end if
   end subroutine init_stochastic_physics_gpu

   subroutine run_stochastic_physics_gpu()
      implicit none
      integer :: k, n, i, jj
      integer :: async_id = 1
!    real  :: glob(nx,my),temp(nxp,my_max)
!    real :: aves,vars,stds

      if (dosppt) then
         call get_random_pattern_run(rpattern_sppt_n2du, rpattern_sppt_n2dv, rpattern_sppt_n2d,  &
          rpattern_sppt_kenorm, rpattern_sppt_spec, rpattern_sppt_varspec, &
          rpattern_sppt_stdev, rpattern_sppt_decortau, rpattern_sppt_lenscale, &
          rpattern_sppt_phi, rpattern_sppt_mlmax, rpattern_sppt_jtrun, &
          rpattern_sppt_mlmax_max, rpattern_sppt_jtrun_max, &
          rpattern_sppt_mlsort, rpattern_sppt_msort, rpattern_sppt_lsort, &
          rpattern_sppt_rstate, rpattern_sppt_seed, nsppt)
         call get_stochy_physics(rpattern_sppt_n2du, rpattern_sppt_n2dv, rpattern_sppt_n2d,  &
          rpattern_sppt_kenorm, rpattern_sppt_spec, rpattern_sppt_varspec, &
          rpattern_sppt_stdev, rpattern_sppt_decortau, rpattern_sppt_lenscale, &
          rpattern_sppt_phi, rpattern_sppt_mlmax, rpattern_sppt_jtrun, &
          rpattern_sppt_mlmax_max, rpattern_sppt_jtrun_max, &
          rpattern_sppt_mlsort, rpattern_sppt_msort, rpattern_sppt_lsort, &
          rpattern_sppt_rstate, rpattern_sppt_seed, nsppt, lev, vfact_sppt, sppt3d)
         if (sppt_logit) then
            !$acc parallel loop collapse(3) async(async_id)
            do jj = 1, my_max
               do k = 1, lev
                  do i = 1, nxp
                     sppt3d(i, k, jj) = (2./(1.+exp(sppt3d(i, k, jj)))) - 1.
                  end do
               end do
            end do
         end if
      end if
      if (doshum) then
         !$acc wait(async_id)
         call get_random_pattern_run(rpattern_shum_n2du, rpattern_shum_n2dv, rpattern_shum_n2d,  &
          rpattern_shum_kenorm, rpattern_shum_spec, rpattern_shum_varspec, &
          rpattern_shum_stdev, rpattern_shum_decortau, rpattern_shum_lenscale, &
          rpattern_shum_phi, rpattern_shum_mlmax, rpattern_shum_jtrun, &
          rpattern_shum_mlmax_max, rpattern_shum_jtrun_max, &
          rpattern_shum_mlsort, rpattern_shum_msort, rpattern_shum_lsort, &
          rpattern_shum_rstate, rpattern_shum_seed, nshum)
         call get_stochy_physics(rpattern_shum_n2du, rpattern_shum_n2dv, rpattern_shum_n2d,  &
          rpattern_shum_kenorm, rpattern_shum_spec, rpattern_shum_varspec, &
          rpattern_shum_stdev, rpattern_shum_decortau, rpattern_shum_lenscale, &
          rpattern_shum_phi, rpattern_shum_mlmax, rpattern_shum_jtrun, &
          rpattern_shum_mlmax_max, rpattern_shum_jtrun_max, &
          rpattern_shum_mlsort, rpattern_shum_msort, rpattern_shum_lsort, &
          rpattern_shum_rstate, rpattern_shum_seed, nshum, lev, vfact_shum, shum3d)
         !$acc wait(async_id)
         if (sppt_logit) then
            !$acc parallel loop collapse(3) async(async_id)
            do jj = 1, my_max
               do k = 1, lev
                  do i = 1, nxp
                     shum3d(i, k, jj) = (2./(1.+exp(shum3d(i, k, jj)))) - 1.
                  end do
               end do
            end do
            !$acc wait(async_id)
         end if
      end if
      if (doskeb) then
         if (first_call) then
            do k = skeblevs - 1, 1, -1
               call get_random_pattern_run_vect(rpattern_skeb_n2du, rpattern_skeb_n2dv, rpattern_skeb_n2d,  &
               rpattern_skeb_kenorm, rpattern_skeb_spec, rpattern_skeb_varspec, &
               rpattern_skeb_stdev, rpattern_skeb_decortau, rpattern_skeb_lenscale, &
               rpattern_skeb_phi, rpattern_skeb_mlmax, rpattern_skeb_jtrun, &
               rpattern_skeb_mlmax_max, rpattern_skeb_jtrun_max, &
               rpattern_skeb_mlsort, rpattern_skeb_msort, rpattern_skeb_lsort, &
               rpattern_skeb_rstate, rpattern_skeb_seed, nskeb, k)
            end do
         end if
         !$acc parallel loop collapse(3) async(async_id)
         do n = 1, nskeb
            do jj = 1, my_max
               do i = 1, nxp
                  !$acc loop seq
                  do k = skeblevs, 2, -1
                     rpattern_skeb_n2du(i, jj, k, n) = rpattern_skeb_n2du(i, jj, k - 1, n)
                     rpattern_skeb_n2dv(i, jj, k, n) = rpattern_skeb_n2dv(i, jj, k - 1, n)
                  end do
               end do
            end do
         end do
         call get_random_pattern_run_vect(rpattern_skeb_n2du, rpattern_skeb_n2dv, rpattern_skeb_n2d,  &
               rpattern_skeb_kenorm, rpattern_skeb_spec, rpattern_skeb_varspec, &
               rpattern_skeb_stdev, rpattern_skeb_decortau, rpattern_skeb_lenscale, &
               rpattern_skeb_phi, rpattern_skeb_mlmax, rpattern_skeb_jtrun, &
               rpattern_skeb_mlmax_max, rpattern_skeb_jtrun_max, &
               rpattern_skeb_mlsort, rpattern_skeb_msort, rpattern_skeb_lsort, &
               rpattern_skeb_rstate, rpattern_skeb_seed, nskeb, 1)
         call get_stochy_physics_vect(rpattern_skeb_n2du, rpattern_skeb_n2dv, rpattern_skeb_n2d,  &
               rpattern_skeb_kenorm, rpattern_skeb_spec, rpattern_skeb_varspec, &
               rpattern_skeb_stdev, rpattern_skeb_decortau, rpattern_skeb_lenscale, &
               rpattern_skeb_phi, rpattern_skeb_mlmax, rpattern_skeb_jtrun, &
               rpattern_skeb_mlmax_max, rpattern_skeb_jtrun_max, &
               rpattern_skeb_mlsort, rpattern_skeb_msort, rpattern_skeb_lsort, &
               rpattern_skeb_rstate, rpattern_skeb_seed, nskeb, lev, vfact_skeb, skeb3du, skeb3dv)
         first_call = .false.
      end if
      if (dossst) then
         call get_random_pattern_run(rpattern_ssst_n2du, rpattern_ssst_n2dv, rpattern_ssst_n2d,  &
          rpattern_ssst_kenorm, rpattern_ssst_spec, rpattern_ssst_varspec, &
          rpattern_ssst_stdev, rpattern_ssst_decortau, rpattern_ssst_lenscale, &
          rpattern_ssst_phi, rpattern_ssst_mlmax, rpattern_ssst_jtrun, &
          rpattern_ssst_mlmax_max, rpattern_ssst_jtrun_max, &
          rpattern_ssst_mlsort, rpattern_ssst_msort, rpattern_ssst_lsort, &
          rpattern_ssst_rstate, rpattern_ssst_seed, nssst)
         call get_stochy_physics(rpattern_ssst_n2du, rpattern_ssst_n2dv, rpattern_ssst_n2d,  &
          rpattern_ssst_kenorm, rpattern_ssst_spec, rpattern_ssst_varspec, &
          rpattern_ssst_stdev, rpattern_ssst_decortau, rpattern_ssst_lenscale, &
          rpattern_ssst_phi, rpattern_ssst_mlmax, rpattern_ssst_jtrun, &
          rpattern_ssst_mlmax_max, rpattern_ssst_jtrun_max, &
          rpattern_ssst_mlsort, rpattern_ssst_msort, rpattern_ssst_lsort, &
          rpattern_ssst_rstate, rpattern_ssst_seed, nssst, 1, vfact_ssst, ssst3d)
         if (sppt_logit) ssst3d(:, 1, :) = (2./(1.+exp(ssst3d(:, 1, :)))) - 1.
      end if

   end subroutine run_stochastic_physics_gpu

   subroutine destroy_stochastic_physics_gpu()
      implicit none
      integer :: async_id = 1
      deallocate (sl)

      if (dosppt) then
         !$acc exit data delete(sppt3d, vfact_sppt, rpattern_sppt_n2d, &
         !$acc&     rpattern_sppt_spec, rpattern_sppt_varspec, rpattern_sppt_mlsort) async(async_id)
         !$acc wait(async_id)
         call get_random_pattern_destroy(rpattern_sppt_n2d, rpattern_sppt_n2du, &
                                         rpattern_sppt_n2dv, rpattern_sppt_kenorm, &
                                         rpattern_sppt_spec, rpattern_sppt_varspec, &
                                         rpattern_sppt_msort, rpattern_sppt_lsort, &
                                         rpattern_sppt_mlsort, nsppt, &
                                         rpattern_sppt_mlmax_max, rpattern_sppt_jtrun_max)
         deallocate (sppt3d)
         deallocate (rpattern_sppt_stdev)
         deallocate (rpattern_sppt_decortau)
         deallocate (rpattern_sppt_lenscale)
         deallocate (rpattern_sppt_seed)
         deallocate (vfact_sppt)
         deallocate (rpattern_sppt_rstate)
         deallocate (rpattern_sppt_phi)
      end if
      if (doshum) then
         !$acc exit data delete(shum3d, shum3d_dq, vfact_shum, rpattern_shum_n2d, &
         !$acc&     rpattern_shum_spec, rpattern_shum_varspec, rpattern_shum_mlsort) async(async_id)
         !$acc wait(async_id)
         call get_random_pattern_destroy(rpattern_shum_n2d, rpattern_shum_n2du, &
                                         rpattern_shum_n2dv, rpattern_shum_kenorm, &
                                         rpattern_shum_spec, rpattern_shum_varspec, &
                                         rpattern_shum_msort, rpattern_shum_lsort, &
                                         rpattern_shum_mlsort, nshum, &
                                         rpattern_shum_mlmax_max, rpattern_shum_jtrun_max)
         deallocate (shum3d)
         deallocate (shum3d_dq)
         deallocate (rpattern_shum_stdev)
         deallocate (rpattern_shum_decortau)
         deallocate (rpattern_shum_lenscale)
         deallocate (rpattern_shum_seed)
         deallocate (vfact_shum)
         deallocate (rpattern_shum_rstate)
         deallocate (rpattern_shum_phi)
      end if
      if (doskeb) then
         !$acc exit data delete(skeb3du, skeb3dv, diss_est, diss_dc, keb, kea, &
         !$acc&     vfact_skeb, skeb_vwts, skeb_vpts, rpattern_skeb_n2du, rpattern_skeb_n2dv, &
         !$acc&     rpattern_skeb_spec, rpattern_skeb_varspec, rpattern_skeb_mlsort) async(async_id)
         !$acc wait(async_id)
         call get_random_pattern_destroy(rpattern_skeb_n2d, rpattern_skeb_n2du, &
                                         rpattern_skeb_n2dv, rpattern_skeb_kenorm, &
                                         rpattern_skeb_spec, rpattern_skeb_varspec, &
                                         rpattern_skeb_msort, rpattern_skeb_lsort, &
                                         rpattern_skeb_mlsort, nskeb, &
                                         rpattern_skeb_mlmax_max, rpattern_skeb_jtrun_max)
         deallocate (skeb3du)
         deallocate (skeb3dv)
         deallocate (diss_est)
         deallocate (diss_dc)
         deallocate (rpattern_skeb_stdev)
         deallocate (rpattern_skeb_decortau)
         deallocate (rpattern_skeb_lenscale)
         deallocate (rpattern_skeb_seed)
         deallocate (vfact_skeb)
         deallocate (keb)
         deallocate (kea)
         deallocate (rpattern_skeb_rstate)
         deallocate (rpattern_skeb_phi)
      end if
      if (dossst) then
         !$acc exit data delete(ssst3d, vfact_ssst, rpattern_ssst_spec, rpattern_ssst_varspec, &
         !$acc&     rpattern_ssst_mlsort) async(async_id)
         !$acc wait(async_id)
         call get_random_pattern_destroy(rpattern_ssst_n2d, rpattern_ssst_n2du, &
                                         rpattern_ssst_n2dv, rpattern_ssst_kenorm, &
                                         rpattern_ssst_spec, rpattern_ssst_varspec, &
                                         rpattern_ssst_msort, rpattern_ssst_lsort, &
                                         rpattern_ssst_mlsort, nssst, &
                                         rpattern_ssst_mlmax_max, rpattern_ssst_jtrun_max)
         deallocate (ssst3d)
         deallocate (rpattern_ssst_stdev)
         deallocate (rpattern_ssst_decortau)
         deallocate (rpattern_ssst_lenscale)
         deallocate (rpattern_ssst_seed)
         deallocate (vfact_ssst)
         deallocate (rpattern_ssst_rstate)
         deallocate (rpattern_ssst_phi)
      end if

   end subroutine destroy_stochastic_physics_gpu

   subroutine init_sppt(dtau)
      implicit none
      real :: dtau
      integer :: n, k

      do n = 1, size(sppt)
         if (sppt(n) > 0) then
            nsppt = nsppt + 1
         else
            exit
         end if
      end do

      allocate (rpattern_sppt_stdev(nsppt))
      allocate (rpattern_sppt_decortau(nsppt))
      allocate (rpattern_sppt_lenscale(nsppt))
      allocate (rpattern_sppt_seed(nsppt))
      allocate (rpattern_sppt_rstate(nsppt))
      allocate (rpattern_sppt_phi(nsppt))
      do n = 1, nsppt
         rpattern_sppt_stdev(n) = sppt(n)
         rpattern_sppt_decortau(n) = sppt_decort(n)
         rpattern_sppt_lenscale(n) = sppt_lscale(n)
         rpattern_sppt_seed(n) = int(sppt_seed(n))
         if (myrank .eq. 0) then
            write (6, *) 'mod_stochastic_physics : sppt : stdev  ', sppt(n)
            write (6, *) 'mod_stochastic_physics : sppt : decort ', sppt_decort(n)
            write (6, *) 'mod_stochastic_physics : sppt : lscale ', sppt_lscale(n)
            write (6, *) 'mod_stochastic_physics : sppt : seed   ', sppt_seed(n)
         end if
      end do

      allocate (sppt3d(nxp, lev, my_max))
      call get_random_pattern_init(rpattern_sppt_n2du, rpattern_sppt_n2dv, rpattern_sppt_n2d,  &
          rpattern_sppt_kenorm, rpattern_sppt_spec, rpattern_sppt_varspec, &
          rpattern_sppt_stdev, rpattern_sppt_decortau, rpattern_sppt_lenscale, &
          rpattern_sppt_phi, rpattern_sppt_mlmax, rpattern_sppt_jtrun, &
          rpattern_sppt_mlmax_max, rpattern_sppt_jtrun_max, &
          rpattern_sppt_mlsort, rpattern_sppt_msort, rpattern_sppt_lsort, &
          rpattern_sppt_rstate, rpattern_sppt_seed, nsppt, dtau, .false.)

      ! set up vfact_sppt
      allocate (vfact_sppt(lev))
      do k = 1, lev
         if (sl(k) .lt. sppt_sigtop1 .and. sl(k) .gt. sppt_sigtop2) then
            vfact_sppt(k) = (sl(k) - sppt_sigtop2)/(sppt_sigtop1 - sppt_sigtop2)
         else if (sl(k) .lt. sppt_sigtop2) then
            vfact_sppt(k) = 0.0
         else
            vfact_sppt(k) = 1.0
         end if
      end do

      if (sppt_sfclimit) then
         ! vfact_sppt(lev-1)=vfact_sppt(lev-2)*0.5
         ! vfact_sppt(lev)=0.0
         do k = 1, lev
            if (sl(k) .lt. sppt_sigbot1 .and. sl(k) .gt. sppt_sigbot2) then
               vfact_sppt(k) = 1.-(sl(k) - sppt_sigbot2)/(sppt_sigbot1 - sppt_sigbot2)
            else if (sl(k) .gt. sppt_sigbot1) then
               vfact_sppt(k) = 0.0
            end if
         end do
      end if

      do k = 1, lev
         if (myrank == 0) print 301, 'mod_stochastic_physics : k,sl,vfact_sppt', k, sl(k), vfact_sppt(k)
      end do
301   format(A, I4, 2F20.12)

   end subroutine init_sppt

   subroutine init_shum(dtau)
      implicit none
      real :: dtau
      integer :: n, k
      do n = 1, size(shum)
         if (shum(n) > 0) then
            nshum = nshum + 1
         else
            exit
         end if
      end do

      allocate (rpattern_shum_stdev(nshum))
      allocate (rpattern_shum_decortau(nshum))
      allocate (rpattern_shum_lenscale(nshum))
      allocate (rpattern_shum_seed(nshum))
      allocate (rpattern_shum_rstate(nshum))
      allocate (rpattern_shum_phi(nshum))
      do n = 1, nshum
         rpattern_shum_stdev(n) = shum(n)
         rpattern_shum_decortau(n) = shum_decort(n)
         rpattern_shum_lenscale(n) = shum_lscale(n)
         rpattern_shum_seed(n) = int(shum_seed(n))
         if (myrank .eq. 0) then
            write (6, *) 'mod_stochastic_physics : shum : stdev  ', shum(n)
            write (6, *) 'mod_stochastic_physics : shum : decort ', shum_decort(n)
            write (6, *) 'mod_stochastic_physics : shum : lscale ', shum_lscale(n)
            write (6, *) 'mod_stochastic_physics : shum : seed   ', shum_seed(n)
         end if
      end do

      allocate (shum3d(nxp, lev, my_max))
      allocate (shum3d_dq(nxp, lev, my_max))
      shum3d_dq = 0.
      call get_random_pattern_init(rpattern_shum_n2du, rpattern_shum_n2dv, rpattern_shum_n2d,  &
          rpattern_shum_kenorm, rpattern_shum_spec, rpattern_shum_varspec, &
          rpattern_shum_stdev, rpattern_shum_decortau, rpattern_shum_lenscale, &
          rpattern_shum_phi, rpattern_shum_mlmax, rpattern_shum_jtrun, &
          rpattern_shum_mlmax_max, rpattern_shum_jtrun_max, &
          rpattern_shum_mlsort, rpattern_shum_msort, rpattern_shum_lsort, &
          rpattern_shum_rstate, rpattern_shum_seed, nshum, dtau, .false.)

      allocate (vfact_shum(lev))
      do k = 1, lev
         vfact_shum(k) = exp((sl(k) - 1.)/shum_sigefold)
         if (sl(k) .LT. 2*shum_sigefold) then
            vfact_shum(k) = 0.0
         end if
         if (myrank == 0) print 301, 'mod_stochastic_physics : k,sl,vfact_shum', k, sl(k), vfact_shum(k)
      end do
301   format(A, I4, 2F20.12)
   end subroutine init_shum

   subroutine init_skeb(dtau)
      implicit none
      real :: dtau
      integer :: n, k, k2
      real, allocatable, dimension(:) :: skeb_vloc

      do n = 1, size(skeb)
         if (skeb(n) > 0) then
            nskeb = nskeb + 1
         else
            exit
         end if
      end do

      allocate (rpattern_skeb_stdev(nskeb))
      allocate (rpattern_skeb_decortau(nskeb))
      allocate (rpattern_skeb_lenscale(nskeb))
      allocate (rpattern_skeb_seed(nskeb))
      allocate (rpattern_skeb_rstate(nskeb))
      allocate (rpattern_skeb_phi(nskeb))
      do n = 1, nskeb
         if (skebnorm .eq. 0) then ! stream function norm
            rpattern_skeb_stdev(n) = skeb(n)*1.111e3*sqrt(dtau)
         end if
         if (skebnorm .eq. 1) then ! kinectic energy function norm
            rpattern_skeb_stdev(n) = skeb(n)*0.00222e3*sqrt(dtau)
         end if
         if (skebnorm .eq. 2) then ! vorticity function norm
            rpattern_skeb_stdev(n) = skeb(n)*1.111e-9*sqrt(dtau)
         end if
         rpattern_skeb_decortau(n) = skeb_decort(n)
         rpattern_skeb_lenscale(n) = skeb_lscale(n)
         rpattern_skeb_seed(n) = int(skeb_seed(n))
         if (myrank .eq. 0) then
            write (6, *) 'mod_stochastic_physics : skeb : stdev  ', skeb(n)
            write (6, *) 'mod_stochastic_physics : skeb : decort ', skeb_decort(n)
            write (6, *) 'mod_stochastic_physics : skeb : lscale ', skeb_lscale(n)
            write (6, *) 'mod_stochastic_physics : skeb : seed   ', skeb_seed(n)
         end if
      end do

      allocate (skeb3du(nxp, lev, my_max))
      allocate (skeb3dv(nxp, lev, my_max))
      allocate (diss_est(nxp, lev, my_max))
      allocate (diss_dc(nxp, lev, my_max))
      allocate (keb(nxp, lev, my_max))
      allocate (kea(nxp, lev, my_max))

      diss_est = 0.
      diss_dc = 0.
      kea = 0.
      keb = 0.

      ! for 3 time level scheme to keep skeblevs same and make sure can be
      ! reproduced when restart ( restart function not ready yet ).
      if (first_call) then
         skeblevs = nint(rpattern_skeb_decortau(1)/(2.*dtau)*skeb_vdof)
      else
         skeblevs = nint(rpattern_skeb_decortau(1)/dtau*skeb_vdof)
      end if

      call get_random_pattern_init(rpattern_skeb_n2du, rpattern_skeb_n2dv, rpattern_skeb_n2d,  &
          rpattern_skeb_kenorm, rpattern_skeb_spec, rpattern_skeb_varspec, &
          rpattern_skeb_stdev, rpattern_skeb_decortau, rpattern_skeb_lenscale, &
          rpattern_skeb_phi, rpattern_skeb_mlmax, rpattern_skeb_jtrun, &
          rpattern_skeb_mlmax_max, rpattern_skeb_jtrun_max, &
          rpattern_skeb_mlsort, rpattern_skeb_msort, rpattern_skeb_lsort, &
          rpattern_skeb_rstate, rpattern_skeb_seed, nskeb, dtau, .true.)

      ! set up vfact_skeb
      allocate (vfact_skeb(lev))
      allocate (skeb_vloc(skeblevs))
      allocate (skeb_vwts(lev, 2))
      allocate (skeb_vpts(lev, 2))
      do k = 1, lev
         if (sl(k) .lt. skeb_sigtop1 .and. sl(k) .gt. skeb_sigtop2) then
            vfact_skeb(k) = (sl(k) - skeb_sigtop2)/(skeb_sigtop1 - skeb_sigtop2)
         else if (sl(k) .lt. skeb_sigtop2) then
            vfact_skeb(k) = 0.0
         else
            vfact_skeb(k) = 1.0
         end if
      end do

      do k = 1, lev
         if (myrank == 0) print 301, 'mod_stochastic_physics : k,sl,vfact_skeb', k, sl(k), vfact_skeb(k)
      end do
301   format(A, I4, 2F20.12)
      ! calculate vertical interpolation weights
      do k = 1, skeblevs
         skeb_vloc(k) = sl(lev) - real(skeblevs - k)/real(skeblevs - 1.0)*(sl(lev) - sl(1))
      end do
      ! surface
      skeb_vwts(lev, 2) = 0.
      skeb_vpts(lev, 1) = skeblevs - 2
      ! top
      skeb_vwts(1, 2) = 1.
      skeb_vpts(1, 1) = 1
      ! internal
      do k = 2, lev - 1
         do k2 = 1, skeblevs - 1
            if (sl(k) .LE. skeb_vloc(k2 + 1) .AND. sl(k) .GT. skeb_vloc(k2)) then
               skeb_vpts(k, 1) = k2
               skeb_vwts(k, 2) = (sl(k) - skeb_vloc(k2))/(skeb_vloc(k2 + 1) - skeb_vloc(k2))
            end if
         end do
      end do
      deallocate (skeb_vloc)
      skeb_vwts(:, 1) = 1.0 - skeb_vwts(:, 2)
      skeb_vpts(:, 2) = skeb_vpts(:, 1) + 1
      if (myrank .eq. 0) then
         do k = 1, lev
            print *, 'skeb vpts ', skeb_vpts(k, 1), skeb_vwts(k, 2)
         end do
      end if

   end subroutine init_skeb

   subroutine init_ssst(dtau)
      implicit none
      real :: dtau
      integer :: n, k

      do n = 1, size(ssst)
         if (ssst(n) > 0) then
            nssst = nssst + 1
         else
            exit
         end if
      end do

      allocate (rpattern_ssst_stdev(nssst))
      allocate (rpattern_ssst_decortau(nssst))
      allocate (rpattern_ssst_lenscale(nssst))
      allocate (rpattern_ssst_seed(nssst))
      allocate (rpattern_ssst_rstate(nssst))
      allocate (rpattern_ssst_phi(nssst))
      do n = 1, nssst
         rpattern_ssst_stdev(n) = ssst(n)
         rpattern_ssst_decortau(n) = ssst_decort(n)
         rpattern_ssst_lenscale(n) = ssst_lscale(n)
         rpattern_ssst_seed(n) = int(ssst_seed(n))
         if (myrank .eq. 0) then
            write (6, *) 'mod_stochastic_physics : ssst : stdev  ', ssst(n)
            write (6, *) 'mod_stochastic_physics : ssst : decort ', ssst_decort(n)
            write (6, *) 'mod_stochastic_physics : ssst : lscale ', ssst_lscale(n)
            write (6, *) 'mod_stochastic_physics : ssst : seed   ', ssst_seed(n)
         end if
      end do

      allocate (ssst3d(nxp, 1, my_max))
      call get_random_pattern_init(rpattern_ssst_n2du, rpattern_ssst_n2dv, rpattern_ssst_n2d,  &
          rpattern_ssst_kenorm, rpattern_ssst_spec, rpattern_ssst_varspec, &
          rpattern_ssst_stdev, rpattern_ssst_decortau, rpattern_ssst_lenscale, &
          rpattern_ssst_phi, rpattern_ssst_mlmax, rpattern_ssst_jtrun, &
          rpattern_ssst_mlmax_max, rpattern_ssst_jtrun_max, &
          rpattern_ssst_mlsort, rpattern_ssst_msort, rpattern_ssst_lsort, &
          rpattern_ssst_rstate, rpattern_ssst_seed, nssst, dtau, .false.)

      allocate (vfact_ssst(1))
      vfact_ssst = 1.

   end subroutine init_ssst

   subroutine get_random_pattern_init(rpattern_n2du, rpattern_n2dv, rpattern_n2d,  &
          rpattern_kenorm, rpattern_spec, rpattern_varspec, &
          rpattern_stdev, rpattern_decortau, rpattern_lenscale, &
          rpattern_phi, rpattern_mlmax, rpattern_jtrun, &
          rpattern_mlmax_max, rpattern_jtrun_max, &
          rpattern_mlsort, rpattern_msort, rpattern_lsort, &
          rpattern_rstate, rpattern_seed, nscale, dt, skebrun)
!---- documentation block
!  purpose: To generate 2D SPPT strucutre
!
!  SPPT 2D 500km example
!       parameter (ncx=128,mcy=ncx/2)
!       parameter (jtrun= 2*((1+(ncx-1)/3)/2), mlmax= jtrun*(jtrun+1)/2)
!  SPPT 2D 1000km example
!       parameter (ncx=64,mcy=ncx/2)
!       parameter (jtrun= 2*((1+(ncx-1)/3)/2), mlmax= jtrun*(jtrun+1)/2)
!  SPPT 2D 2000km example
!       parameter (ncx=32,mcy=ncx/2)
!       parameter (jtrun= 2*((1+(ncx-1)/3)/2), mlmax= jtrun*(jtrun+1)/2)
!---- end of documentation block

      implicit none
      integer :: timearray(3), iseed
      integer :: n, k, nscale, ncx, ml, ms, ns, i, j
      real :: rerth, pi, var, correLsq, rkT, rnn1
      real(kind=RTYPE), allocatable :: rpattern_n2du(:, :, :, :)
      real(kind=RTYPE), allocatable :: rpattern_n2dv(:, :, :, :)
      real(kind=RTYPE), allocatable :: rpattern_n2d(:, :, :)
      real, allocatable :: rpattern_kenorm(:, :, :)
      real(kind=RTYPE), allocatable :: rpattern_spec(:, :, :)
      real(kind=RTYPE), allocatable :: rpattern_varspec(:, :)
      real:: rpattern_stdev(nscale)
      real:: rpattern_decortau(nscale)
      real :: rpattern_lenscale(nscale)
      real :: rpattern_phi(nscale)
      integer, allocatable :: rpattern_mlmax(:)
      integer, allocatable :: rpattern_jtrun(:)
      integer :: rpattern_mlmax_max
      integer :: rpattern_jtrun_max
      integer, allocatable :: rpattern_mlsort(:, :, :)
      integer, allocatable :: rpattern_msort(:, :)
      integer, allocatable :: rpattern_lsort(:, :)
      type(random_stat) :: rpattern_rstate(nscale)
      integer :: rpattern_seed(nscale)
      integer :: irand
      real :: dt
      real(kind=RTYPE), allocatable :: noise(:, :)
      integer(8) count, count_rate, count_max, count_trunc
      integer(8) :: iscale = 10000000000_8
      integer :: count4
      logical :: skebrun

      rerth = 6.3712e+6      ! radius of earth (m)
      pi = 4.*atan(1.)
      allocate(rpattern_jtrun(nscale))
      allocate(rpattern_mlmax(nscale))
      rpattern_jtrun_max = 0
      rpattern_mlmax_max = 0
      do n = 1, nscale
         ncx = 2.*pi*rerth/rpattern_lenscale(n)
         !  rpattern(n)_jtrun = 2*((1+(ncx-1)/3)/2)
         rpattern_jtrun(n) = 2*((1 + (4*ncx - 1)/4)/2)
         rpattern_mlmax(n) = rpattern_jtrun(n)*(rpattern_jtrun(n) + 1)/2
         if (rpattern_jtrun(n) .gt. rpattern_jtrun_max) rpattern_jtrun_max = rpattern_jtrun(n)
         if (rpattern_mlmax(n) .gt. rpattern_mlmax_max) rpattern_mlmax_max = rpattern_mlmax(n)
      end do
      allocate (rpattern_n2d(nxp, my_max, nscale))
      if (doskeb) then
         allocate (rpattern_n2du(nxp, my_max, skeblevs, nscale))
         allocate (rpattern_n2dv(nxp, my_max, skeblevs, nscale))
         allocate (rpattern_kenorm(rpattern_mlmax_max, 2, nscale))
      end if
      allocate (rpattern_spec(rpattern_mlmax_max, 2, nscale))
      allocate (rpattern_varspec(rpattern_mlmax_max, nscale))
      allocate (rpattern_msort(rpattern_mlmax_max, nscale))
      allocate (rpattern_lsort(rpattern_mlmax_max, nscale))
      allocate (rpattern_mlsort(rpattern_jtrun_max, rpattern_jtrun_max, nscale))
      !    radsq = rerth*rerth
      do n = 1, nscale
         allocate (noise(rpattern_mlmax(n), 2))
         ncx = 2.*pi*rerth/rpattern_lenscale(n)
         !  rpattern(n)_jtrun = 2*((1+(ncx-1)/3)/2)
         rpattern_jtrun(n) = 2*((1 + (4*ncx - 1)/4)/2)
         rpattern_mlmax(n) = rpattern_jtrun(n)*(rpattern_jtrun(n) + 1)/2

         if (doskeb) then
            rpattern_n2du(:, :, :, n) = 0.0
            rpattern_n2dv(:, :, :, n) = 0.0
         end if

         ! Real random seeds
         if (myrank .eq. 0) then
            if (.not. ncep_seeds) then
               call itime(timearray)
               iseed = irand(0)
               count4 = irand(timearray(1) + (i*11467 - iseed)*timearray(2) &
                              + (iseed - i*23)*timearray(3))
               write(*,*) timearray
            else
               call system_clock(count, count_rate, count_max)
               count_trunc = iscale*(count/iscale)
               count4 = count - count_trunc
            end if
         end if
!     call mpe_bcast(count4,1,0,mpe_double)
         call mpe_bcast(count4, 1, 0, mpe_integer)
         if (rpattern_seed(n) == -999) then
            rpattern_seed(n) = count4
         end if
         if (do_unit_test) rpattern_seed(n) = 1000
         if (myrank .eq. 0) write (6, *) 'scale', n, ' : using seed :', rpattern_seed(n)
         call random_setseed(rpattern_seed(n), rpattern_rstate(n))

         ! horizontal decorrelation
         correLsq = rpattern_lenscale(n)*rpattern_lenscale(n)
         rkT = 0.25*correLsq/radsq
         if (myrank .eq. 0) write (6, *) 'mod_stochastic_physics : n,rkT = ', n, rkT

         ! time decorrelation
         rpattern_phi(n) = exp(-dt/rpattern_decortau(n))

         call sortml_gpu(rpattern_jtrun(n), rpattern_mlmax(n), &
                     rpattern_msort(:, n), rpattern_lsort(:, n), rpattern_mlsort(:,:,n), &
                     rpattern_jtrun_max, rpattern_mlmax_max)

         noise = 0.
         do ml = 1, rpattern_mlmax(n)
            noise(ml, 1) = 1./sqrt(float(2*(rpattern_lsort(ml, n) - 1)) + 1)
            noise(ml, 2) = 1./sqrt(float(2*(rpattern_lsort(ml, n) - 1)) + 1)
            if (rpattern_msort(ml, n) .eq. 1) then
               noise(ml, 1) = sqrt(2.)/sqrt(float(2*(rpattern_lsort(ml, n) - 1)) + 1)
               noise(ml, 2) = 0.
            end if
         end do

         noise(1, 1) = 0.
         noise(1, 2) = 0.
         noise = noise*sqrt(1./float(rpattern_jtrun(n)))

         ! set up the amplitude of noise
         do ml = 1, rpattern_mlmax(n)
            rpattern_varspec(ml, n) = sqrt(float(rpattern_jtrun(n)) &
                                           *exp(-rkT*float(rpattern_lsort(ml, n))*(float(rpattern_lsort(ml, n) - 1))))
         end do
         do ml = 1, rpattern_mlmax(n)
            noise(ml, 1) = noise(ml, 1)*rpattern_varspec(ml, n)
            noise(ml, 2) = noise(ml, 2)*rpattern_varspec(ml, n)
         end do


         ! get specral variance
         var = 0.
         do ml = 1, rpattern_mlmax(n)
            if (rpattern_msort(ml, n) .ne. 1) then
               var = var + (noise(ml, 1)**2 + noise(ml, 2)**2)
            else
               var = var + 0.5*(noise(ml, 1)**2 + noise(ml, 2)**2)
            end if
         end do
         do ml = 1, rpattern_mlmax(n)
            rpattern_varspec(ml, n) = rpattern_varspec(ml, n)/sqrt(var)
         end do

         !if (myrank .eq. 0) write (6, *) 'GACtotal variance =', var
         !if (myrank .eq. 0) write (6, *) 'GADrpattern_varspec(n) =', rpattern_varspec(n)
#ifdef VERBOSE
#endif
         ! initialize spectrum coefficient
         noise = 0.
         call get_noise(rpattern_mlmax(n), rpattern_jtrun(n), &
          rpattern_msort(:,n), rpattern_lsort(:,n), &
          rpattern_rstate(n), rpattern_seed(n), noise, rpattern_jtrun_max, rpattern_mlmax_max)
         do ml = 1, rpattern_mlmax(n)
            rpattern_spec(ml, 1, n) = rpattern_stdev(n)*rpattern_varspec(ml, n)*noise(ml, 1)
            rpattern_spec(ml, 2, n) = rpattern_stdev(n)*rpattern_varspec(ml, n)*noise(ml, 2)
         end do
         rpattern_spec(1, 1, n) = 0.
         rpattern_spec(1, 2, n) = 0.

         deallocate (noise)
!
         if (skebrun) then
            rpattern_kenorm(:, :, n) = 1.
            if (skebnorm .eq. 0) then
               do j = 1, rpattern_jtrun(n)
                  do i = j, rpattern_jtrun(n)
                     ml = rpattern_mlsort(j, i, n)
                     rnn1 = float(i*(i + 1))
                     rpattern_kenorm(ml, 1, n) = rnn1/radsq
                     rpattern_kenorm(ml, 2, n) = rnn1/radsq
                  end do
               end do
               if (myrank .eq. 0) print *, 'using streamfunction ', &
                  maxval(rpattern_kenorm(:, 1, n)), minval(rpattern_kenorm(:, 1, n))
            end if
            if (skebnorm .eq. 1) then
               do j = 1, rpattern_jtrun(n)
                  do i = j, rpattern_jtrun(n)
                     ml = rpattern_mlsort(j, i, n)
                     rnn1 = float(i*(i + 1))
                     rpattern_kenorm(ml, 1, n) = sqrt(rnn1)/rerth
                     rpattern_kenorm(ml, 2, n) = sqrt(rnn1)/rerth
                  end do
               end do
               if (myrank .eq. 0) print *, 'using kenorm ', &
                  maxval(rpattern_kenorm(:, 1, n)), minval(rpattern_kenorm(:, 1, n))
            end if
         end if
      end do

      if (myrank .eq. 0) then
         write (6, *) 'mod_stochastic_physics : dt    = ', dt
         write (6, *) 'mod_stochastic_physics : stdev = ', rpattern_stdev
         write (6, *) 'mod_stochastic_physics : seed  = ', rpattern_seed
         write (6, *) 'mod_stochastic_physics : jtrun = ', rpattern_jtrun
         write (6, *) 'mod_stochastic_physics : mlmax = ', rpattern_mlmax
         write (6, *) 'mod_stochastic_physics : tau   = ', rpattern_decortau
         write (6, *) 'mod_stochastic_physics : phi   = ', rpattern_phi
      end if

   end subroutine get_random_pattern_init

   subroutine get_random_pattern_destroy(rpattern_n2d, rpattern_n2du, rpattern_n2dv, &
      rpattern_kenorm, rpattern_spec, rpattern_varspec, rpattern_msort, rpattern_lsort, &
      rpattern_mlsort, nscale, rpattern_mlmax_max, rpattern_jtrun_max)
      implicit none
      integer :: n, nscale, rpattern_mlmax_max, rpattern_jtrun_max
      real(kind=RTYPE), allocatable :: rpattern_n2du(:,:,:,:)
      real(kind=RTYPE), allocatable :: rpattern_n2dv(:,:,:,:)
      real(kind=RTYPE), allocatable :: rpattern_n2d(:,:,:)
      real, allocatable :: rpattern_kenorm(:,:,:)
      real(kind=RTYPE), allocatable :: rpattern_spec(:,:,:)
      real(kind=RTYPE), allocatable :: rpattern_varspec(:,:)
      integer, allocatable :: rpattern_mlsort(:,:,:)
      integer, allocatable :: rpattern_msort(:,:)
      integer, allocatable :: rpattern_lsort(:,:)

      deallocate (rpattern_n2d)
      if (doskeb) then
         deallocate (rpattern_n2du)
         deallocate (rpattern_n2dv)
         deallocate (rpattern_kenorm)
      end if
      deallocate (rpattern_spec)
      deallocate (rpattern_varspec)
      deallocate (rpattern_msort)
      deallocate (rpattern_lsort)
      deallocate (rpattern_mlsort)

   end subroutine get_random_pattern_destroy

   subroutine get_random_pattern_run(rpattern_n2du, rpattern_n2dv, rpattern_n2d,  &
          rpattern_kenorm, rpattern_spec, rpattern_varspec, &
          rpattern_stdev, rpattern_decortau, rpattern_lenscale, &
          rpattern_phi, rpattern_mlmax, rpattern_jtrun, &
          rpattern_mlmax_max, rpattern_jtrun_max, &
          rpattern_mlsort, rpattern_msort, rpattern_lsort, &
          rpattern_rstate, rpattern_seed, nscale)
      implicit none
      integer :: nscale
      integer :: n, ii, i, jj, j, k, nxj
      real(kind=RTYPE) :: rpattern_n2du(nxp, my_max, skeblevs, nscale)
      real(kind=RTYPE) :: rpattern_n2dv(nxp, my_max, skeblevs, nscale)
      real(kind=RTYPE) :: rpattern_n2d(nxp, my_max, nscale)
      real :: rpattern_kenorm(rpattern_mlmax_max, 2, nscale)
      real(kind=RTYPE) :: rpattern_spec(rpattern_mlmax_max, 2, nscale)
      real(kind=RTYPE) :: rpattern_varspec(rpattern_mlmax_max, nscale)
      real :: rpattern_stdev(nscale)
      real :: rpattern_decortau(nscale)
      real :: rpattern_lenscale(nscale)
      real :: rpattern_phi(nscale)
      integer :: rpattern_mlmax(nscale)
      integer :: rpattern_jtrun(nscale)
      integer :: rpattern_mlmax_max
      integer :: rpattern_jtrun_max
      integer :: rpattern_mlsort(rpattern_jtrun_max, rpattern_jtrun_max, nscale)
      integer :: rpattern_msort(rpattern_mlmax_max, nscale)
      integer :: rpattern_lsort(rpattern_mlmax_max, nscale)
      type(random_stat) :: rpattern_rstate(nscale)
      integer :: rpattern_seed(nscale)

      do n = 1, nscale
         call gen_random_pattern_2d(rpattern_n2du(:,:,:,n), rpattern_n2dv(:,:,:,n), rpattern_n2d(:,:,n),  &
          rpattern_kenorm(:,:,n), rpattern_spec(:,:,n), rpattern_varspec(:,n), &
          rpattern_stdev(n), rpattern_decortau(n), rpattern_lenscale(n), &
          rpattern_phi(n), rpattern_mlmax(n), rpattern_jtrun(n), &
          rpattern_mlsort(:,:,n), rpattern_msort(:,n), rpattern_lsort(:,n), &
          rpattern_rstate(n), rpattern_seed(n), rpattern_jtrun_max, rpattern_mlmax_max)
      end do

   end subroutine get_random_pattern_run

   subroutine get_random_pattern_run_vect(rpattern_n2du, rpattern_n2dv, rpattern_n2d,  &
          rpattern_kenorm, rpattern_spec, rpattern_varspec, &
          rpattern_stdev, rpattern_decortau, rpattern_lenscale, &
          rpattern_phi, rpattern_mlmax, rpattern_jtrun, &
          rpattern_mlmax_max, rpattern_jtrun_max, &
          rpattern_mlsort, rpattern_msort, rpattern_lsort, &
          rpattern_rstate, rpattern_seed, nscale, k)
      implicit none
      integer :: nscale
      integer :: n, ii, i, jj, j, k, nxj, k2
      real(kind=RTYPE) :: rpattern_n2du(nxp, my_max, skeblevs, nscale)
      real(kind=RTYPE) :: rpattern_n2dv(nxp, my_max, skeblevs, nscale)
      real(kind=RTYPE) :: rpattern_n2d(nxp, my_max, nscale)
      real :: rpattern_kenorm(rpattern_mlmax_max, 2, nscale)
      real(kind=RTYPE) :: rpattern_spec(rpattern_mlmax_max, 2, nscale)
      real(kind=RTYPE) :: rpattern_varspec(rpattern_mlmax_max, nscale)
      real :: rpattern_stdev(nscale)
      real :: rpattern_decortau(nscale)
      real :: rpattern_lenscale(nscale)
      real :: rpattern_phi(nscale)
      integer :: rpattern_mlmax(nscale)
      integer :: rpattern_jtrun(nscale)
      integer :: rpattern_mlmax_max
      integer :: rpattern_jtrun_max
      integer :: rpattern_mlsort(rpattern_jtrun_max, rpattern_jtrun_max, nscale)
      integer :: rpattern_msort(rpattern_mlmax_max, nscale)
      integer :: rpattern_lsort(rpattern_mlmax_max, nscale)
      type(random_stat) :: rpattern_rstate(nscale)
      integer :: rpattern_seed(nscale)

      do n = 1, nscale
         call gen_random_pattern_2d_vect(rpattern_n2du(:,:,:,n), rpattern_n2dv(:,:,:,n), rpattern_n2d(:,:,n),  &
          rpattern_kenorm(:,:,n), rpattern_spec(:,:,n), rpattern_varspec(:,n), &
          rpattern_stdev(n), rpattern_decortau(n), rpattern_lenscale(n), &
          rpattern_phi(n), rpattern_mlmax(n), rpattern_jtrun(n), &
          rpattern_mlsort(:,:,n), rpattern_msort(:,n), rpattern_lsort(:,n), &
          rpattern_rstate(n), rpattern_seed(n), rpattern_jtrun_max, rpattern_mlmax_max, k)
      end do

   end subroutine get_random_pattern_run_vect

   subroutine get_stochy_physics(rpattern_n2du, rpattern_n2dv, rpattern_n2d,  &
          rpattern_kenorm, rpattern_spec, rpattern_varspec, &
          rpattern_stdev, rpattern_decortau, rpattern_lenscale, &
          rpattern_phi, rpattern_mlmax, rpattern_jtrun, &
          rpattern_mlmax_max, rpattern_jtrun_max, &
          rpattern_mlsort, rpattern_msort, rpattern_lsort, &
          rpattern_rstate, rpattern_seed, nscale, nlev, vfact, n3d)
!------------------------------------------------------------------------!
!  purpose: To generate 3D SPPT strucutre
!  output: n3d
!------------------------------------------------------------------------!
      implicit none
      integer :: n, ii, i, jj, j, k, nxj, nlev
      integer, intent(in) :: nscale
      real, intent(in) :: vfact(nlev)
      real(kind=RTYPE) :: rpattern_n2du(nxp, my_max, skeblevs, nscale)
      real(kind=RTYPE) :: rpattern_n2dv(nxp, my_max, skeblevs, nscale)
      real(kind=RTYPE) :: rpattern_n2d(nxp, my_max, nscale)
      real :: rpattern_kenorm(rpattern_mlmax_max, 2, nscale)
      real(kind=RTYPE) :: rpattern_spec(rpattern_mlmax_max, 2, nscale)
      real(kind=RTYPE) :: rpattern_varspec(rpattern_mlmax_max, nscale)
      real :: rpattern_stdev(nscale)
      real :: rpattern_decortau(nscale)
      real :: rpattern_lenscale(nscale)
      real :: rpattern_phi(nscale)
      integer :: rpattern_mlmax(nscale)
      integer :: rpattern_jtrun(nscale)
      integer :: rpattern_mlmax_max
      integer :: rpattern_jtrun_max
      integer :: rpattern_mlsort(rpattern_jtrun_max, rpattern_jtrun_max, nscale)
      integer :: rpattern_msort(rpattern_mlmax_max, nscale)
      integer :: rpattern_lsort(rpattern_mlmax_max, nscale)
      type(random_stat) :: rpattern_rstate(nscale)
      integer :: rpattern_seed(nscale)
      real(kind=RTYPE), intent(out) :: n3d(nxp, nlev, my_max)
      real(kind=RTYPE) :: n3dtmp
      integer :: async_id = 1
      !$acc parallel loop collapse(2) private(nxj, j) async(async_id)
      do jj = 1, jlistnum
         do k = 1, nlev
            j = jlist1(jj)
            nxj = nxdef_2d(j)
            !$acc loop vector private(n3dtmp)
            do i = 1, nxj
               n3dtmp = 0.
               !$acc loop seq
               do n = 1, nscale
                  n3dtmp = n3dtmp + rpattern_n2d(i, jj, n)*vfact(k)
               end do
               n3d(i, k, jj) = n3dtmp
            end do
         end do
      end do

   end subroutine get_stochy_physics

   subroutine get_stochy_physics_vect(rpattern_n2du, rpattern_n2dv, rpattern_n2d,  &
          rpattern_kenorm, rpattern_spec, rpattern_varspec, &
          rpattern_stdev, rpattern_decortau, rpattern_lenscale, &
          rpattern_phi, rpattern_mlmax, rpattern_jtrun, &
          rpattern_mlmax_max, rpattern_jtrun_max, &
          rpattern_mlsort, rpattern_msort, rpattern_lsort, &
          rpattern_rstate, rpattern_seed, nscale, nlev, vfact, n3du, n3dv)
!------------------------------------------------------------------------!
!  purpose: To generate 3D SKEB strucutre
!  output: n3du n3dv
!------------------------------------------------------------------------!
      implicit none
      integer :: n, ii, i, jj, j, k, nxj, nlev
      integer, intent(in) :: nscale
      real, intent(in) :: vfact(nlev)
      real(kind=RTYPE) :: rpattern_n2du(nxp, my_max, skeblevs, nscale)
      real(kind=RTYPE) :: rpattern_n2dv(nxp, my_max, skeblevs, nscale)
      real(kind=RTYPE) :: rpattern_n2d(nxp, my_max, nscale)
      real :: rpattern_kenorm(rpattern_mlmax_max, 2, nscale)
      real(kind=RTYPE) :: rpattern_spec(rpattern_mlmax_max, 2, nscale)
      real(kind=RTYPE) :: rpattern_varspec(rpattern_mlmax_max, nscale)
      real :: rpattern_stdev(nscale)
      real :: rpattern_decortau(nscale)
      real :: rpattern_lenscale(nscale)
      real :: rpattern_phi(nscale)
      integer :: rpattern_mlmax(nscale)
      integer :: rpattern_jtrun(nscale)
      integer :: rpattern_mlmax_max
      integer :: rpattern_jtrun_max
      integer :: rpattern_mlsort(rpattern_jtrun_max, rpattern_jtrun_max, nscale)
      integer :: rpattern_msort(rpattern_mlmax_max, nscale)
      integer :: rpattern_lsort(rpattern_mlmax_max, nscale)
      type(random_stat) :: rpattern_rstate(nscale)
      integer :: rpattern_seed(nscale)
      real(kind=RTYPE), intent(out) :: n3du(nxp, nlev, my_max), n3dv(nxp, nlev, my_max)
      integer :: async_id = 1

      !$acc parallel loop collapse(3) async(async_id)
      do jj = 1, my_max
         do k = 1, nlev
            do i = 1, nxp
               n3du(i, k, jj) = 0.
               n3dv(i, k, jj) = 0.
            end do
         end do
      end do
      !$acc parallel loop collapse(2) gang private(j, nxj) async(async_id)
      do jj = 1, jlistnum
         do k = 1, nlev
            j = jlist1(jj)
            nxj = nxdef_2d(j)
            !$acc loop vector
            do i = 1, nxj
               !$acc loop seq
               do n = 1, nscale
                  n3du(i, k, jj) = n3du(i, k, jj) + vfact(k) &
                                   *(skeb_vwts(k, 1)*rpattern_n2du(i, jj, skeb_vpts(k, 1), n) &
                                     + skeb_vwts(k, 2)*rpattern_n2du(i, jj, skeb_vpts(k, 2), n))
                  n3dv(i, k, jj) = n3dv(i, k, jj) + vfact(k) &
                                   *(skeb_vwts(k, 1)*rpattern_n2dv(i, jj, skeb_vpts(k, 1), n) &
                                     + skeb_vwts(k, 2)*rpattern_n2dv(i, jj, skeb_vpts(k, 2), n))
               end do
            end do
         end do
      end do

   end subroutine get_stochy_physics_vect

   subroutine get_noise(rpattern_mlmax, rpattern_jtrun, &
          rpattern_msort, rpattern_lsort, &
          rpattern_rstate, rpattern_seed, noise, rpattern_jtrun_max, rpattern_mlmax_max)
      implicit none
      integer :: ml, ns, ms, rpattern_jtrun_max, rpattern_mlmax_max
      integer :: rpattern_mlmax, rpattern_jtrun, rpattern_msort(rpattern_mlmax_max), &
         rpattern_lsort(rpattern_mlmax_max), rpattern_seed
      type(random_stat) :: rpattern_rstate
      real(kind=RTYPE), intent(out) :: noise(rpattern_mlmax, 2)
      real :: noise_gauss(2*rpattern_mlmax)
      real :: ave, var, std

      ! get white noise  with a gaussian (normal) distribution
      call random_gauss(noise_gauss, rpattern_rstate)
      noise_gauss(1) = 0.; noise_gauss(rpattern_mlmax + 1) = 0.
      noise_gauss = noise_gauss*sqrt(1./float(rpattern_jtrun - 1))

      noise = 0.
      do ml = 1, rpattern_mlmax
         ! set up to red noise
         noise(ml, 1) = noise_gauss(ml)/sqrt(float(2*(rpattern_lsort(ml) - 1) + 1))
         noise(ml, 2) = noise_gauss(rpattern_mlmax + ml)/sqrt(float(2*(rpattern_lsort(ml) - 1) + 1))
         ! zero out the imagenary part when m(zonal wavenumber) equal to 1
         if (rpattern_msort(ml) .eq. 1) then
            noise(ml, 1) = sqrt(2.)*noise(ml, 1)
            noise(ml, 2) = 0.
         end if
      end do
   end subroutine get_noise

   subroutine gen_random_pattern_2d(rpattern_n2du, rpattern_n2dv, rpattern_n2d,  &
          rpattern_kenorm, rpattern_spec, rpattern_varspec, &
          rpattern_stdev, rpattern_decortau, rpattern_lenscale, &
          rpattern_phi, rpattern_mlmax, rpattern_jtrun, &
          rpattern_mlsort, rpattern_msort, rpattern_lsort, &
          rpattern_rstate, rpattern_seed, rpattern_jtrun_max, rpattern_mlmax_max)
      use spec_cuda_graph, only: cc_cg, gwk1_cg
      implicit none
      real(kind=RTYPE) :: rpattern_n2du(nxp, my_max, skeblevs)
      real(kind=RTYPE) :: rpattern_n2dv(nxp, my_max, skeblevs)
      real(kind=RTYPE) :: rpattern_n2d(nxp, my_max)
      real :: rpattern_kenorm(rpattern_mlmax_max, 2)
      real(kind=RTYPE) :: rpattern_spec(rpattern_mlmax_max, 2)
      real(kind=RTYPE) :: rpattern_varspec(rpattern_mlmax_max)
      real :: rpattern_stdev
      real :: rpattern_decortau
      real :: rpattern_lenscale
      real :: rpattern_phi
      integer :: rpattern_mlmax
      integer :: rpattern_jtrun
      integer :: rpattern_mlmax_max
      integer :: rpattern_jtrun_max
      integer :: rpattern_mlsort(rpattern_jtrun_max, rpattern_jtrun_max)
      integer :: rpattern_msort(rpattern_mlmax_max)
      integer :: rpattern_lsort(rpattern_mlmax_max)
      type(random_stat) :: rpattern_rstate
      integer :: rpattern_seed
      !real(kind=RTYPE), intent(out) :: sppt2d(nxp,my_max)
      integer :: ml, ns, ms
      real(kind=RTYPE) :: noise(rpattern_mlmax, 2)
      real(kind=RTYPE) :: specp(jtrun, jtmax, 2), bufr2d(jtrun, jtmax*nsizey, 2)
      integer :: async_id = 1
      real(kind=RTYPE) :: n2d(nxp, my_max)
      integer :: i, j, k, iy, global_j, flat_idx
      real(kind=RTYPE) :: speci(rpattern_mlmax, 2)
      



      !$acc data create(bufr2d, specp, speci) async(async_id)
      if (col_rank .eq. 0) then
         ! get noise
         call get_noise(rpattern_mlmax, rpattern_jtrun, &
          rpattern_msort, rpattern_lsort, &
          rpattern_rstate, rpattern_seed, noise, rpattern_jtrun_max, rpattern_mlmax_max)  ! remain on CPU
         !$acc wait(async_id)
         !  radom pattern advance with first order AR
         !$acc parallel loop collapse(2) copyin(noise) async(async_id)
         do j = 1, 2
            do i = 1, rpattern_mlmax
               rpattern_spec(i, j) = rpattern_phi*rpattern_spec(i, j) + &
                                    sqrt(1.-rpattern_phi**2.)*rpattern_stdev*rpattern_varspec(i)*noise(i, j)
               speci(i, j) = rpattern_spec(i, j)
            end do
         end do

         ! ready for mpi_scatter random pattern
         call spectrun_inp2d(rpattern_jtrun, jtrun, jtmax &
                             , rpattern_mlsort, nsizey &
                             , speci, bufr2d, rpattern_jtrun_max, async_id)

      end if

      !  mpi_scatter random pattern from root
      call nccl_scatter_gpu(bufr2d, specp, jtrun*jtmax*2, nsizey, &
         nccl_col_comm, col_rank, async_id) ! on CPU

      ! transform spectral to physical space
      call transr1_gpu(jtrun, jtmax, nx, my, my_max, polyf, specp, rpattern_n2d, &
                       nsizey, cc_cg, gwk1_cg) ! on GPU
      !$acc end data

   end subroutine gen_random_pattern_2d

   subroutine gen_random_pattern_2d_vect(rpattern_n2du, rpattern_n2dv, rpattern_n2d,  &
          rpattern_kenorm, rpattern_spec, rpattern_varspec, &
          rpattern_stdev, rpattern_decortau, rpattern_lenscale, &
          rpattern_phi, rpattern_mlmax, rpattern_jtrun, &
          rpattern_mlsort, rpattern_msort, rpattern_lsort, &
          rpattern_rstate, rpattern_seed, rpattern_jtrun_max, rpattern_mlmax_max, k)

      implicit none
      real(kind=RTYPE) :: rpattern_n2du(nxp, my_max, skeblevs)
      real(kind=RTYPE) :: rpattern_n2dv(nxp, my_max, skeblevs)
      real(kind=RTYPE) :: rpattern_n2d(nxp, my_max)
      real :: rpattern_kenorm(rpattern_mlmax_max, 2)
      real(kind=RTYPE) :: rpattern_spec(rpattern_mlmax_max, 2)
      real(kind=RTYPE) :: rpattern_varspec(rpattern_mlmax_max)
      real :: rpattern_stdev
      real :: rpattern_decortau
      real :: rpattern_lenscale
      real :: rpattern_phi
      integer :: rpattern_mlmax
      integer :: rpattern_jtrun
      integer :: rpattern_mlmax_max
      integer :: rpattern_jtrun_max
      integer :: rpattern_mlsort(rpattern_jtrun_max, rpattern_jtrun_max)
      integer :: rpattern_msort(rpattern_mlmax_max)
      integer :: rpattern_lsort(rpattern_mlmax_max)
      type(random_stat) :: rpattern_rstate
      integer :: rpattern_seed
      integer :: ml, ns, ms, k, i, j, jj, nxj, kk
      real    :: xx
      real(kind=RTYPE) :: specpv(jtrun, jtmax, 2), specpd(jtrun, jtmax, 2)
      real(kind=RTYPE) :: bufr2d(jtrun, jtmax*nsizey, 2), noise(rpattern_mlmax, 2), specf(rpattern_mlmax, 2)
      integer :: async_id = 1

     
      !$acc data create(specpd, specpv, bufr2d, specf) async(async_id)
      !$acc parallel loop collapse(3) async(async_id)
      do j = 1, 2
          do kk = 1, jtmax
             do i = 1, jtrun
                specpd(i, kk, j) = 0.0
                specpv(i, kk, j) = 0.0
             end do
          end do
      end do
      if (col_rank .eq. 0) then
         ! get noise
         call get_noise(rpattern_mlmax, rpattern_jtrun, &
          rpattern_msort, rpattern_lsort, &
          rpattern_rstate, rpattern_seed, noise, rpattern_jtrun_max, rpattern_mlmax_max)

         !  radom pattern advance with first order AR
         !$acc wait(async_id)
         !$acc parallel loop async(async_id)
         do i = 1, rpattern_mlmax
            rpattern_spec(i, 1) = rpattern_phi*rpattern_spec(i, 1) + &
                                 sqrt(1.-rpattern_phi**2.)*rpattern_stdev*rpattern_varspec(i)*noise(i, 1)
            rpattern_spec(i, 2) = rpattern_phi*rpattern_spec(i, 2) + &
                                 sqrt(1.-rpattern_phi**2.)*rpattern_stdev*rpattern_varspec(i)*noise(i, 2)

            specf(i, 1) = rpattern_kenorm(i, 1)*rpattern_spec(i, 1)
            specf(i, 2) = rpattern_kenorm(i, 2)*rpattern_spec(i, 2)
         end do
         ! ready for mpi_scatter random pattern
         call spectrun_inp2d(rpattern_jtrun, jtrun, jtmax &
                             , rpattern_mlsort, nsizey &
                             , specf, bufr2d, rpattern_jtrun_max, async_id)

      end if

      !  mpi_scatter random pattern from root
      !call mpe_scatter_sppt(bufr2d, specpv, 2*jtrun*jtmax, nsizey)
      call nccl_scatter_gpu(bufr2d, specpv, jtrun*jtmax*2, nsizey, &
         nccl_col_comm, col_rank, async_id) ! on CPU

      !  spectral transform for velocity components
      call tranuv1_gpu(jtrun, jtmax, nx, my, my_max, 1, onocos, wcfac, wdfac &
                   , poly, dpoly, specpv, specpd, rpattern_n2du(:, :, k), rpattern_n2dv(:, :, k), nsizey)
      !$acc parallel loop gang private(j, nxj, xx) async(async_id)
      do jj = 1, jlistnum
         j = jlist1(jj)
         nxj = nxdef_2d(j)
         xx = rad/cosl(j)
         !$acc loop vector
         do i = 1, nxj
            rpattern_n2du(i, jj, k) = rpattern_n2du(i, jj, k)*xx
            rpattern_n2dv(i, jj, k) = rpattern_n2dv(i, jj, k)*xx
         end do
      end do
      !$acc end data
      !$acc wait(async_id)


   end subroutine gen_random_pattern_2d_vect

   SUBROUTINE avevar_sppt(data, n, ave, var, std)
      INTEGER n
      REAL ave, var, data(n)
      INTEGER j
      REAL s, ep, std
      ave = 0.0
      do j = 1, n
         ave = ave + data(j)
      end do
      ave = ave/n
      var = 0.0
      ep = 0.0
      do j = 1, n
         s = data(j) - ave
         ep = ep + s
         var = var + s*s
      end do
      var = (var - ep**2/n)/n
      std = sqrt(var)

   end SUBROUTINE avevar_sppt
!
!
   SUBROUTINE avevar_sppt2d(data2d, n, m, ave, var, std)
      INTEGER :: n, m, nmdim
      REAL :: ave, var, data(n*m)
      REAL(kind=RTYPE) :: data2d(n, m)
      INTEGER :: i, j
      REAL :: s, ep, std

      nmdim = 0
      do j = 1, m
         do i = 1, n
            nmdim = nmdim + 1
            data(nmdim) = data2d(i, j)
         end do
      end do

      ave = 0.0
      do j = 1, nmdim
         ave = ave + data(j)
      end do
      ave = ave/nmdim
      var = 0.0
      ep = 0.0
      do j = 1, nmdim
         s = data(j) - ave
         ep = ep + s
         var = var + s*s
      end do
!   var=(var-ep**2/n)/(n-1)    ! sample numners < 30
      var = (var - ep**2/nmdim)/nmdim
      std = sqrt(var)
   end SUBROUTINE avevar_sppt2d
!
   subroutine spptout(tau)
      implicit none
      integer      :: i, j, k, jj, nxj, ihead, n, async_id = 1
      integer      :: nxmy4
      real         :: tau
      real*4       :: glob4(nx, my)
      real(kind=RTYPE) :: glob(nx, my), temp(nxp, my_max)

      !$acc wait(async_id)
      !$acc update self(rpattern_sppt_n2d, sppt3d) async(async_id)
      !$acc wait(async_id)

      ihead = 15
      nxmy4 = nx*my*4
      if (myrank .eq. 0) then
         open (ihead, file='sppt.dat', access='direct', form='unformatted' &
               , recl=nxmy4, status='unknown', convert='big_endian')
      end if

      do n = 1, nsppt
         call unify_reduceintp(nx, my, my_max, rpattern_sppt_n2d(:,:, n), glob)
         if (myrank .eq. 0) then
            glob4 = glob
            write (ihead, rec=recnsppt) glob4
            recnsppt = recnsppt + 1
         end if
      end do
      do k = lev, 1, -1
         do jj = 1, jlistnum
            j = jlist1(jj)
            nxj = nxdef_2d(j)
            do i = 1, nxj
               temp(i, jj) = sppt3d(i, k, jj)
            end do
         end do
         call unify_reduceintp(nx, my, my_max, temp, glob)
         if (myrank .eq. 0) then
            glob4 = glob
            write (ihead, rec=recnsppt) glob4
            recnsppt = recnsppt + 1
         end if
      end do

      if (myrank .eq. 0) close (ihead)

      !creating ctl file
      call spptctl(tau)

   end subroutine spptout
!----
   subroutine shumout(tau)
      implicit none
      integer      :: i, j, k, jj, nxj, ihead, n, async_id = 1
      integer      :: nxmy4
      real         :: tau
      real*4       :: glob4(nx, my)
      real(kind=RTYPE) :: glob(nx, my), temp(nxp, my_max)

      !$acc wait(async_id)
      !$acc update self(rpattern_shum_n2d, shum3d_dq) async(async_id)
      !$acc wait(async_id)

      ihead = 15
      nxmy4 = nx*my*4
      if (myrank .eq. 0) then
         open (ihead, file='shum.dat', access='direct', form='unformatted' &
               , recl=nxmy4, status='unknown', convert='big_endian')
      end if

      do n = 1, nshum
         call unify_reduceintp(nx, my, my_max, rpattern_shum_n2d(:,:,n), glob)
         if (myrank .eq. 0) then
            glob4 = glob
            write (ihead, rec=recnshum) glob4
            recnshum = recnshum + 1
         end if
      end do
      do k = lev, 1, -1
         do jj = 1, jlistnum
            j = jlist1(jj)
            nxj = nxdef_2d(j)
            do i = 1, nxj
               temp(i, jj) = shum3d_dq(i, k, jj)
            end do
         end do
         call unify_reduceintp(nx, my, my_max, temp, glob)
         if (myrank .eq. 0) then
            glob4 = glob
            write (ihead, rec=recnshum) glob4
            recnshum = recnshum + 1
         end if
      end do

      if (myrank .eq. 0) close (ihead)

      !creating ctl file
      call shumctl(tau)

   end subroutine shumout
!----
   subroutine skebest_gpu(um, vm)
      use grid, only: ut, vt
      use spec_cuda_graph, only: cc_cg, gwk1_cg, wcc_fk_cg, &
                           wc_cg, ws_cg, fj_weight_cg
      implicit none

      integer :: n, ii, i, jj, j, k, nxj
      real    :: xx, axx
      real(kind=RTYPE) :: spectmp(levp, 2, jtrun, jtmax), &
                          cc(nx + 2, levp, 1, my_max), &
                          um(nxp, lev, my_max), vm(nxp, lev, my_max), &
                          dummy, temp
      integer :: async_id = 1
      !$acc data create(diss_est, spectmp) async(async_id)
      ! estimate the dissipation of kinectic energy for SKEB
      !$acc parallel loop collapse(2) gang private(j, nxj, xx) async(async_id)
      do jj = 1, jlistnum
         do k = 1, lev
            j = jlist1(jj)
            nxj = nxdef_2d(j)
            xx = radsq*onocos(j)
            !$acc loop vector
            do i = 1, nxj
               diss_est(i, k, jj) = ((um(i, k, jj)*ut(i, k, jj) &
                                      + vm(i, k, jj)*vt(i, k, jj)) &
                                     + 0.5*(um(i, k, jj)**2.+vm(i, k, jj)**2.)) &
                                    *xx
            end do
         end do
      end do
      call joinrs_gpu(cc_cg, diss_est, dummy, dummy, dummy, nx, my_max, lev &
                  , jlistnum, 1, 1)
      call tranrs_gpu_cuda_graph(jtrun, jtmax, nx, my, my_max, levp, polyf, weight, &
                                 cc_cg, spectmp, 1, nsizey, &
                                 gwk1_cg, ws_cg, wc_cg, wcc_fk_cg, fj_weight_cg)
      ! apply spectral filter of the dissipation of kinetic energy
      call filter_skeb_gpu(jtrun, jtmax, levp, spectmp, skebfilt, async_id)
      call transr_gpu_cuda_graph(jtrun, jtmax, nx, my, my_max, levp, polyf, &
                                 spectmp, cc_cg, 1, nsizey, &
                                 gwk1_cg, ws_cg, wc_cg, wcc_fk_cg)
      call ujoinsr_gpu(cc_cg, diss_est, dummy, dummy, dummy, nx, my_max, lev, jlistnum, 1, 1)
      !
!        if ( myrank .eq. 0 ) print *,'intgrt: diss_est(1,72,1)=',diss_est(1,72,1)
      !$acc parallel loop gang collapse(2) private(j, nxj, xx, axx) async(async_id)
      do jj = 1, jlistnum
         do k = 1, lev
            j = jlist1(jj)
            nxj = nxdef_2d(j)
            xx = rad/cosl(j)
            axx = cosl(j)/rad
            !$acc loop vector private(temp)
            do i = 1, nxj
               temp = diss_est(i, k, jj)
               if (doskeb_dc) then
                  temp = diss_est(i, k, jj) + diss_dc(i, k, jj)
               end if
               !change virtual wind to real wind
               ut(i, k, jj) = ut(i, k, jj)*xx
               vt(i, k, jj) = vt(i, k, jj)*xx
               keb(i, k, jj) = 0.5*(ut(i, k, jj)**2.+vt(i, k, jj)**2.)
!              ut(i,k,jj)=ut(i,k,jj)+skeb3du(i,k,jj)*diss_est(i,k,jj)
!              vt(i,k,jj)=vt(i,k,jj)+skeb3dv(i,k,jj)*diss_est(i,k,jj)
               ut(i, k, jj) = ut(i, k, jj) + skeb3du(i, k, jj)*temp
               vt(i, k, jj) = vt(i, k, jj) + skeb3dv(i, k, jj)*temp
               kea(i, k, jj) = 0.5*(ut(i, k, jj)**2.+vt(i, k, jj)**2.)
               !change back to virtual wind
               ut(i, k, jj) = ut(i, k, jj)*axx
               vt(i, k, jj) = vt(i, k, jj)*axx
            end do
         end do
      end do
      !$acc end data
      !$acc wait(async_id)
!        if ( myrank .eq. 0 ) print *,'intgrt: keb(1,72,1)=',keb(1,72,1)
!        if ( myrank .eq. 0 ) print *,'intgrt: kea(1,72,1)=',kea(1,72,1)

   end subroutine skebest_gpu
!----
   subroutine skebout(tau)
      implicit none
      integer      :: i, j, k, jj, nxj, ihead, n, async_id = 1
      integer      :: nxmy4
      real         :: tau, xx
      real(kind=RTYPE):: glob(nx, my), temp(nxp, my_max)
      real(kind=4) :: glob4(nx, my)

      !$acc wait(async_id)
      !$acc update self(skeb3du, skeb3dv, diss_est, keb, kea, diss_dc) async(async_id)
      !$acc wait(async_id)

      ihead = 15
      nxmy4 = nx*my*4
      if (myrank .eq. 0) then
         open (ihead, file='skeb.dat', access='direct', form='unformatted' &
               , recl=nxmy4, status='unknown', convert='big_endian')
      end if

      do k = lev, 1, -1
         do jj = 1, jlistnum
            j = jlist1(jj)
            nxj = nxdef_2d(j)
            do i = 1, nxj
               temp(i, jj) = skeb3du(i, k, jj)
            end do
         end do
         call unify_reduceintp(nx, my, my_max, temp, glob)
         if (myrank .eq. 0) then
            glob4 = glob
            write (ihead, rec=recnskeb) glob4
            recnskeb = recnskeb + 1
         end if
      end do

      do k = lev, 1, -1
         do jj = 1, jlistnum
            j = jlist1(jj)
            nxj = nxdef_2d(j)
            do i = 1, nxj
               temp(i, jj) = skeb3dv(i, k, jj)
            end do
         end do
         call unify_reduceintp(nx, my, my_max, temp, glob)
         if (myrank .eq. 0) then
            glob4 = glob
            write (ihead, rec=recnskeb) glob4
            recnskeb = recnskeb + 1
         end if
      end do

      do k = lev, 1, -1
         do jj = 1, jlistnum
            j = jlist1(jj)
            nxj = nxdef_2d(j)
            do i = 1, nxj
               temp(i, jj) = diss_est(i, k, jj)
            end do
         end do
         call unify_reduceintp(nx, my, my_max, temp, glob)
         if (myrank .eq. 0) then
            glob4 = glob
            write (ihead, rec=recnskeb) glob4
            recnskeb = recnskeb + 1
         end if
      end do

      do k = lev, 1, -1
         do jj = 1, jlistnum
            j = jlist1(jj)
            nxj = nxdef_2d(j)
            do i = 1, nxj
               temp(i, jj) = keb(i, k, jj)
            end do
         end do
         call unify_reduceintp(nx, my, my_max, temp, glob)
         if (myrank .eq. 0) then
            glob4 = glob
            write (ihead, rec=recnskeb) glob4
            recnskeb = recnskeb + 1
         end if
      end do

      do k = lev, 1, -1
         do jj = 1, jlistnum
            j = jlist1(jj)
            nxj = nxdef_2d(j)
            do i = 1, nxj
               temp(i, jj) = kea(i, k, jj)
            end do
         end do
         call unify_reduceintp(nx, my, my_max, temp, glob)
         if (myrank .eq. 0) then
            glob4 = glob
            write (ihead, rec=recnskeb) glob4
            recnskeb = recnskeb + 1
         end if
      end do

      do k = lev, 1, -1
         do jj = 1, jlistnum
            j = jlist1(jj)
            nxj = nxdef_2d(j)
            do i = 1, nxj
               temp(i, jj) = diss_dc(i, k, jj)
            end do
         end do
         call unify_reduceintp(nx, my, my_max, temp, glob)
         if (myrank .eq. 0) then
            glob4 = glob
            write (ihead, rec=recnskeb) glob4
            recnskeb = recnskeb + 1
         end if
      end do
      if (myrank .eq. 0) close (ihead)

      !creating ctl file
      call skebctl(tau)

   end subroutine skebout
!
!----
   subroutine spptctl(tau)
!
      use const, only: idtg, sinl, sigma
      use rank, only: myrank

      integer :: nxj, j, k, itau, ch, iter, remd, js, je, kk, lev1, mn, n
      real :: pi, r2d, dlon, tau
      real, allocatable :: mlat(:), prsl(:)
      character(len=30) ::  forydef, forzdef

      character yy*4, dd*2, hh*2, mm*2, dtg*12
      character*3 mon(12)
      logical jrem

      data mon/'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep' &
         , 'Oct', 'Nov', 'Dec'/

      allocate (mlat(my), prsl(lev))

      itau = tau
      ch = 10
      lev1 = 1
      pi = 4.0*atan(1.0)
      r2d = 180./pi
      dlon = 360./float(nx)

      do j = 1, my
         mlat(j) = asin(sinl(j))*r2d
      end do

      do k = 1, lev
         kk = lev - k + 1
         prsl(kk) = sigma(k, 2) + sigma(k + 1, 2)
         prsl(kk) = prsl(kk) + (sigma(k, 1) + sigma(k + 1, 1))*1000.
         prsl(kk) = 0.5*prsl(kk)
      end do

      write (forydef, 8) my
      write (forzdef, 9) lev

      write (dtg, '(I12)') idtg
      read (dtg, '(A4,I2,A2,A2,A2)') yy, mn, dd, hh, mm

      if (myrank .eq. 0) then
         OPEN (UNIT=ch, FILE='sppt.ctl', STATUS='UNKNOWN' &
               , ACCESS='SEQUENTIAL')

         write (ch, '(A14)') 'dset ^sppt.dat'
         write (ch, '(A18)') 'options big_endian'
         write (ch, '(A12)') 'undef -999.0'
         write (ch, 12) 'ydef', my, 'levels'
         iter = my/8
         remd = mod(my, 8)
         jrem = (remd .eq. 0)
         write (forydef, 8) remd
         do j = 1, iter
            js = 1 + 8*(j - 1)
            je = js + 7
            write (ch, 10) mlat(js:je)
         end do
         if (.not. jrem) write (10, forydef) mlat(je + 1:my)

         write (ch, 13) 'xdef', nx, 'linear 0.0', dlon
         write (ch, 14) 'tdef', itau, 'linear', hh, 'Z', dd, mon(mn), yy, '1hr'
         write (ch, 12) 'zdef', lev, 'levels '
         iter = lev/8
         remd = mod(lev, 8)
         jrem = (remd .eq. 0)
         write (forzdef, 9) remd
         do j = 1, iter
            js = 1 + 8*(j - 1)
            je = js + 7
            write (ch, 11) prsl(js:je)
         end do
         if (.not. jrem) write (ch, forzdef) prsl(je + 1:my)
         write (ch, '(A4,1X,I2)') 'vars', nsppt + 1
         do n = 1, nsppt
            write (ch, 15) 'scale', n, lev1, '99', rpattern_sppt_lenscale(n)/1000., 'km 2D Random Pattern'
         end do
         write (ch, 16) 'sppt  ', lev, '99', '3D Random Pattern'
         write (ch, '(A7)') 'endvars'

         close (ch)
      end if

8     format("(", I4, "(2x,F11.7))")
9     format("(", I4, "(2x,F10.5))")
10    format(8(2x, F11.7))
11    format(8(2x, F10.5))
12    format(A4, 1X, I4, 1X, A6)
13    format(A4, 1X, I4, 1X, A10, 1X, F10.7)
14    format(A4, 1X, I4, 1X, A6, 1X, A2, A1, A2, A3, A4, 1X, A3)
15    format(A5, I1, 3X, I3, 1X, A2, 1X, F7.1, A20)
16    format(A6, 3X, I3, 1X, A2, 1X, A20)

      deallocate (mlat, prsl)
   end subroutine spptctl
!
!----
!----
   subroutine shumctl(tau)
!
      use const, only: idtg, sinl, sigma
      use rank, only: myrank

      integer :: nxj, j, k, itau, ch, iter, remd, js, je, kk, lev1, mn, n
      real :: pi, r2d, dlon, tau
      real, allocatable :: mlat(:), prsl(:)
      character(len=30) ::  forydef, forzdef

      character yy*4, dd*2, hh*2, mm*2, dtg*12
      character*3 mon(12)
      logical jrem

      data mon/'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep' &
         , 'Oct', 'Nov', 'Dec'/

      allocate (mlat(my), prsl(lev))

      itau = tau
      ch = 10
      lev1 = 1
      pi = 4.0*atan(1.0)
      r2d = 180./pi
      dlon = 360./float(nx)

      do j = 1, my
         mlat(j) = asin(sinl(j))*r2d
      end do

      do k = 1, lev
         kk = lev - k + 1
         prsl(kk) = sigma(k, 2) + sigma(k + 1, 2)
         prsl(kk) = prsl(kk) + (sigma(k, 1) + sigma(k + 1, 1))*1000.
         prsl(kk) = 0.5*prsl(kk)
      end do

      write (forydef, 8) my
      write (forzdef, 9) lev

      write (dtg, '(I12)') idtg
      read (dtg, '(A4,I2,A2,A2,A2)') yy, mn, dd, hh, mm

      if (myrank .eq. 0) then
         OPEN (UNIT=ch, FILE='shum.ctl', STATUS='UNKNOWN' &
               , ACCESS='SEQUENTIAL')

         write (ch, '(A14)') 'dset ^shum.dat'
         write (ch, '(A18)') 'options big_endian'
         write (ch, '(A12)') 'undef -999.0'
         write (ch, 12) 'ydef', my, 'levels'
         iter = my/8
         remd = mod(my, 8)
         jrem = (remd .eq. 0)
         write (forydef, 8) remd
         do j = 1, iter
            js = 1 + 8*(j - 1)
            je = js + 7
            write (ch, 10) mlat(js:je)
         end do
         if (.not. jrem) write (ch, forydef) mlat(je + 1:my)

         write (ch, 13) 'xdef', nx, 'linear 0.0', dlon
         write (ch, 14) 'tdef', itau, 'linear', hh, 'Z', dd, mon(mn), yy, '1hr'
         write (ch, 12) 'zdef', lev, 'levels '
         iter = lev/8
         remd = mod(lev, 8)
         jrem = (remd .eq. 0)
         write (forzdef, 9) remd
         do j = 1, iter
            js = 1 + 8*(j - 1)
            je = js + 7
            write (ch, 11) prsl(js:je)
         end do
         if (.not. jrem) write (ch, forzdef) prsl(je + 1:my)
         write (ch, '(A4,1X,I2)') 'vars', nshum + 1
         do n = 1, nshum
            write (ch, 15) 'scale', n, lev1, '99', rpattern_shum_lenscale(n)/1000., 'km 2D Random Pattern'
         end do
         write (ch, 16) 'shum3d_dq  ', lev, '99', '3D delta q'
         write (ch, '(A7)') 'endvars'

         close (ch)
      end if

8     format("(", I4, "(2x,F11.7))")
9     format("(", I4, "(2x,F10.5))")
10    format(8(2x, F11.7))
11    format(8(2x, F10.5))
12    format(A4, 1X, I4, 1X, A6)
13    format(A4, 1X, I4, 1X, A10, 1X, F10.7)
14    format(A4, 1X, I4, 1X, A6, 1X, A2, A1, A2, A3, A4, 1X, A3)
15    format(A5, I1, 3X, I3, 1X, A2, 1X, F7.1, A20)
16    format(A6, 3X, I3, 1X, A2, 1X, A20)

      deallocate (mlat, prsl)
   end subroutine shumctl

!----
!----
   subroutine skebctl(tau)
!
      use const, only: idtg, sinl, sigma
      use rank, only: myrank

      integer :: nxj, j, k, itau, ch, iter, remd, js, je, kk, lev1, mn, n
      real :: pi, r2d, dlon, tau
      real, allocatable :: mlat(:), prsl(:)
      character(len=30) ::  forydef, forzdef

      character yy*4, dd*2, hh*2, mm*2, dtg*12
      character*3 mon(12)
      logical jrem

      data mon/'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep' &
         , 'Oct', 'Nov', 'Dec'/

      allocate (mlat(my), prsl(lev))

      itau = tau
      ch = 10
      lev1 = 1
      pi = 4.0*atan(1.0)
      r2d = 180./pi
      dlon = 360./float(nx)

      do j = 1, my
         mlat(j) = asin(sinl(j))*r2d
      end do

      do k = 1, lev
         kk = lev - k + 1
         prsl(kk) = sigma(k, 2) + sigma(k + 1, 2)
         prsl(kk) = prsl(kk) + (sigma(k, 1) + sigma(k + 1, 1))*1000.
         prsl(kk) = 0.5*prsl(kk)
      end do

      write (forydef, 8) my
      write (forzdef, 9) lev

      write (dtg, '(I12)') idtg
      read (dtg, '(A4,I2,A2,A2,A2)') yy, mn, dd, hh, mm

      if (myrank .eq. 0) then
         OPEN (UNIT=ch, FILE='skeb.ctl', STATUS='UNKNOWN' &
               , ACCESS='SEQUENTIAL')

         write (ch, '(A14)') 'dset ^skeb.dat'
         write (ch, '(A18)') 'options big_endian'
         write (ch, '(A12)') 'undef -999.0'
         write (ch, 12) 'ydef', my, 'levels'
         iter = my/8
         remd = mod(my, 8)
         jrem = (remd .eq. 0)
         write (forydef, 8) remd
         do j = 1, iter
            js = 1 + 8*(j - 1)
            je = js + 7
            write (ch, 10) mlat(js:je)
         end do
         if (.not. jrem) write (ch, forydef) mlat(je + 1:my)

         write (ch, 13) 'xdef', nx, 'linear 0.0', dlon
         write (ch, 14) 'tdef', itau, 'linear', hh, 'Z', dd, mon(mn), yy, '1hr'
         write (ch, 12) 'zdef', lev, 'levels '
         iter = lev/8
         remd = mod(lev, 8)
         jrem = (remd .eq. 0)
         write (forzdef, 9) remd
         do j = 1, iter
            js = 1 + 8*(j - 1)
            je = js + 7
            write (ch, 11) prsl(js:je)
         end do
         if (.not. jrem) write (ch, forzdef) prsl(je + 1:my)
         write (ch, '(A4,1X,I2)') 'vars', 6
         write (ch, 16) 'skebu  ', lev, '99', 'U-dir 3D Random Pattern'
         write (ch, 16) 'skebv  ', lev, '99', 'V-dir 3D Random Pattern'
         write (ch, 16) 'dissest', lev, '99', 'dissipation of KE      '
         write (ch, 16) 'keb    ', lev, '99', 'kE before SKEB         '
         write (ch, 16) 'kea    ', lev, '99', 'KE after SKEB          '
         write (ch, 16) 'dissdc ', lev, '99', 'dissipation from deep convenction'
         write (ch, '(A7)') 'endvars'

         close (ch)
      end if

8     format("(", I4, "(2x,F11.7))")
9     format("(", I4, "(2x,F10.5))")
10    format(8(2x, F11.7))
11    format(8(2x, F10.5))
12    format(A4, 1X, I4, 1X, A6)
13    format(A4, 1X, I4, 1X, A10, 1X, F10.7)
14    format(A4, 1X, I4, 1X, A6, 1X, A2, A1, A2, A3, A4, 1X, A3)
15    format(A5, I1, 3X, I3, 1X, A2, 1X, F7.1, A20)
16    format(A7, 3X, I3, 1X, A2, 1X, A23)

      deallocate (mlat, prsl)
   end subroutine skebctl
!
   subroutine spectrun_inp2d(jcap1, jtr, jtm, mlsort, ns, speci, speco, &
      rpattern_jtrun_max, async_id)
!
! use spectral truncation to change resoltuion
!
      implicit none
      integer lev, jtr, jcap1, jtm, ns, ml, rpattern_jtrun_max
      real(kind=RTYPE) speci(jcap1*(jcap1 + 1)/2, 2)
      real(kind=RTYPE) speco(jtr, jtm*ns*2)
      integer i, j, k, jj, jp, jr, j1, j2
      integer mlsort(rpattern_jtrun_max, rpattern_jtrun_max)
      integer :: async_id
!
      !$acc parallel loop collapse(2) async(async_id)
      do j = 1, jtm*ns*2
         do i = 1, jtr
            speco(i, j) = 0.0
         end do
      end do
      if (jcap1 .gt. jtr) then
         !$acc parallel loop gang private(jj, jp, jr, j1, j2) async(async_id)
         do j = 1, jcap1
            if (j .le. jtr) then
               jj = nlist(j)
               jp = (jj - 1)/jtm
               jr = mod(jj - 1, jtm) + 1
               j1 = jp*jtm*2 + jr
               j2 = j1 + jtm
            end if
            !$acc loop vector private(ml)
            do i = j, jcap1
               ml = mlsort(j, i)
               if (i .le. jtr) then
                  speco(i, j1) = speci(ml, 1)
                  speco(i, j2) = speci(ml, 2)
               end if
            end do
         end do
      else if (jcap1 .lt. jtr) then
         !$acc parallel loop gang private(jj, jp, jr, j1, j2) async(async_id)
         do j = 1, jtr
            jj = nlist(j)
            jp = (jj - 1)/jtm
            jr = mod(jj - 1, jtm) + 1
            j1 = jp*jtm*2 + jr
            j2 = j1 + jtm
            !$acc loop vector private(ml)
            do i = j, jtr
               if (i .le. jcap1) then
                  ml = mlsort(j, i)
                  speco(i, j1) = speci(ml, 1)
                  speco(i, j2) = speci(ml, 2)
               end if
            end do
         end do
      else      ! jcap1=jtr
         !$acc parallel loop gang private(jj, jp, jr, j1, j2) async(async_id)
         do j = 1, jtr
            jj = nlist(j)
            jp = (jj - 1)/jtm
            jr = mod(jj - 1, jtm) + 1
            j1 = jp*jtm*2 + jr
            j2 = j1 + jtm
            !$acc loop vector private(ml)
            do i = j, jtr
               ml = mlsort(j, i)
               speco(i, j1) = speci(ml, 1)
               speco(i, j2) = speci(ml, 2)
            end do
         end do
      end if

      return
   end subroutine spectrun_inp2d

 subroutine nccl_scatter_gpu(ainp, aout, len, nsize, nccl_col_comm, col_rank, async_id)
   use const, only: RTYPE, MPI_RTYPE
   use param, only: nx, my, lev, my_max
   use index, only: nxp
   use rank, only: nccl_comm_gfs, myrank
   use cudafor
   use openacc
   use nccl
   implicit none
   integer, intent(in):: len, nsize, col_rank
   real(kind=RTYPE) :: aout(len)
   real(kind=RTYPE) :: ainp(len, nsize)
   type(ncclComm), intent(in) :: nccl_col_comm
   integer i, pts
   integer(kind=cuda_stream_kind) :: stream
   integer :: async_id

   pts = len
   stream = acc_get_cuda_stream(async_id)
   if (col_rank == 0) then
      !$acc kernels async(async_id)
      aout(1:pts) = ainp(1:pts, 1)
      !$acc end kernels
   end if
   NCCLCHECK(ncclGroupStart())
   !$acc host_data use_device(ainp, aout)
   if (col_rank == 0) then
      do i = 1, nsize
         if ((i-1) /= col_rank) then
            NCCLCHECK(ncclSend(ainp(1,i), pts, ncclFloat64, i-1, nccl_col_comm, stream))
         end if
      end do
   end if

   if (col_rank /= 0) then
      NCCLCHECK(ncclRecv(aout, pts, ncclFloat64, 0, nccl_col_comm, stream))
   end if

   !$acc end host_data
   NCCLCHECK(ncclGroupEnd())
   return
 end subroutine nccl_scatter_gpu

end module mod_stochastic_physics_gpu
