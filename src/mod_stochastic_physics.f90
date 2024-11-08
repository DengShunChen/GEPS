module mod_stochastic_physics
  use mpe, only : mpe_bcast, mpe_double, mpe_integer
  use rank, only : myrank
  use index
  use param
  use const, only : aki, bki, first_call, dosppt, doshum, doskeb, dossst, &
                    poly, dpoly, wdfac, wcfac, onocos, radsq, weight,     &
                    RTYPE, rad, cosl, doskeb_dc 
  use mersenne_twister, only: random_setseed,random_gauss,random_stat
  implicit none
  private 

  type random_pattern
    real(kind=RTYPE), allocatable :: n2du(:,:,:)
    real(kind=RTYPE), allocatable :: n2dv(:,:,:)
    real(kind=RTYPE), allocatable :: n2d(:,:)
    real, allocatable :: kenorm(:,:)
    real(kind=RTYPE), allocatable :: spec(:,:)
    real(kind=RTYPE), allocatable :: varspec(:)
    real :: stdev ! stochastic physics tendency amplitude
    real :: decortau ! time scales
    real :: lenscale ! length scales
    real :: phi
    integer :: mlmax
    integer :: jtrun
    integer,allocatable :: mlsort(:,:)
    integer,allocatable :: msort(:)
    integer,allocatable :: lsort(:)
    type(random_stat),public :: rstate
    integer, public :: seed
  end type random_pattern

  integer,save :: recnsppt=1, recnskeb=1, recnshum=1
  real ::  dt
  logical, public :: ncep_seeds=.false.
  real,allocatable :: sl(:)

  type(random_pattern), public, save, allocatable, dimension(:) :: &
       rpattern_sppt, rpattern_shum, rpattern_skeb, rpattern_ssst

  ! SPPT
  integer :: nsppt
  real(kind=RTYPE), allocatable, save :: sppt3d(:,:,:)
  real :: sppt(5) = -999.             ! amplitude(0.~1.)
  real :: sppt_seed(5) = -999.        ! random seeds
  real :: sppt_decort(5) = -999.      ! time scales(seconds)
  real :: sppt_lscale(5) = -999.      ! length scales(meters)
  real, allocatable, dimension(:) :: vfact_sppt
  real, public :: sppt_sigtop1 = 0.1
  real, public :: sppt_sigtop2 = 0.025
  real, public :: sppt_sigbot1 = 0.975
  real, public :: sppt_sigbot2 = 0.9
  logical, public :: sppt_sfclimit=.false.
  logical, public :: sppt_logit=.true.

  ! SHUM
  integer :: nshum
  real(kind=RTYPE), allocatable, save :: shum3d(:,:,:)
  real :: shum(5) = -999.             ! amplitude(0.~1.)
  real :: shum_seed(5) = -999.         ! random seeds
  real :: shum_decort(5) = -999.      ! time scales(seconds)
  real :: shum_lscale(5) = -999.      ! length scales(meters)
  real, allocatable, dimension(:) :: vfact_shum
  real, public :: shum_sigefold = 0.2
  real(kind=RTYPE), allocatable, save :: shum3d_dq(:,:,:)

  ! SKEB
  integer :: nskeb,skeblevs
  real(kind=RTYPE), allocatable, save :: skeb3du(:,:,:),skeb3dv(:,:,:),diss_est(:,:,:)
  real(kind=RTYPE), allocatable, save :: diss_dc(:,:,:)    !dissipation from deep convection
  real(kind=RTYPE), allocatable, save :: kea(:,:,:),keb(:,:,:)
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
  real, public, allocatable,dimension(:,:) :: skeb_vwts
  integer, public, allocatable,dimension(:,:) :: skeb_vpts

  ! SSST
  integer :: nssst
  real(kind=RTYPE), allocatable, save :: ssst3d(:,:,:)
  real :: ssst(5) = -999.             ! amplitude(0.~1.)
  real :: ssst_seed(5) = -999.         ! random seeds
  real :: ssst_decort(5) = -999.      ! time scales(seconds)
  real :: ssst_lscale(5) = -999.      ! length scales(meters)
  real, allocatable, dimension(:) :: vfact_ssst

  public random_pattern

  public nsppt, sppt, sppt_seed, sppt_decort, sppt_lscale, sppt3d
  public nshum, shum, shum_seed, shum_decort, shum_lscale, shum3d, shum3d_dq
  public nskeb, skeb, skeb_seed, skeb_decort, skeb_lscale,        &
         skeb3du, skeb3dv, diss_est, diss_dc, &
         skeblevs, keb, kea
  public nssst, ssst, ssst_seed, ssst_decort, ssst_lscale, ssst3d

  public  init_stochastic_physics, &
           run_stochastic_physics, &
       destroy_stochastic_physics

  public spptout,shumout,skebout,skebest
  public avevar_sppt2d

contains

  subroutine init_stochastic_physics(dtau)
    implicit none
    integer :: n, k 
    real :: dtau

    allocate(sl(lev))
    ! calculation sigma values
    do k=1,lev
      sl(k)=0.5*(aki(k)/1013.0+bki(k)+aki(k+1)/1013.0+bki(k+1))
    enddo
 
    if (dosppt) then
      call init_sppt(dtau)
    endif

    if (doshum) then
      call init_shum(dtau)
    endif

    if (doskeb) then
      call init_skeb(dtau)
    endif

    if (dossst) then
      call init_ssst(dtau)
    endif

  end subroutine init_stochastic_physics

  subroutine run_stochastic_physics()
    implicit none
    integer :: k,n
!    real  :: glob(nx,my),temp(nxp,my_max)
!    real :: aves,vars,stds

    if (dosppt) then
      call get_random_pattern_run(rpattern_sppt,nsppt)
      call get_stochy_physics(rpattern_sppt,nsppt,lev,vfact_sppt,sppt3d)
      if (sppt_logit) sppt3d(:,:,:) = (2./(1.+exp(sppt3d(:,:,:))))-1.
    endif
    if (doshum) then
      call get_random_pattern_run(rpattern_shum,nshum)
      call get_stochy_physics(rpattern_shum,nshum,lev,vfact_shum,shum3d)
      if (sppt_logit) shum3d(:,:,:) = (2./(1.+exp(shum3d(:,:,:))))-1.
    endif
    if (doskeb) then
      if ( first_call ) then
        do k=skeblevs-1,1,-1
          call get_random_pattern_run_vect(rpattern_skeb,nskeb,k)
        enddo
      endif
      do n=1,nskeb
        do k=skeblevs,2,-1
          rpattern_skeb(n)%n2du(:,:,k)=rpattern_skeb(n)%n2du(:,:,k-1)
          rpattern_skeb(n)%n2dv(:,:,k)=rpattern_skeb(n)%n2dv(:,:,k-1)
        enddo
      enddo
      call get_random_pattern_run_vect(rpattern_skeb,nskeb,1)
      call get_stochy_physics_vect(rpattern_skeb,nskeb,lev,vfact_skeb,skeb3du,skeb3dv)
      first_call=.false.
    endif
    if (dossst) then
      call get_random_pattern_run(rpattern_ssst,nssst)
      call get_stochy_physics(rpattern_ssst,nssst,1  ,vfact_ssst,ssst3d)
      if (sppt_logit) ssst3d(:,1,:) = (2./(1.+exp(ssst3d(:,1,:))))-1.
    endif

  end subroutine run_stochastic_physics

  subroutine destroy_stochastic_physics()
    implicit none
    deallocate(sl)

    if (dosppt) then
      call get_random_pattern_destroy(rpattern_sppt,nsppt)
      deallocate(sppt3d)
      deallocate(rpattern_sppt)
      deallocate(vfact_sppt)
    endif
    if (doshum) then
      call get_random_pattern_destroy(rpattern_shum,nshum)
      deallocate(shum3d)
      deallocate(shum3d_dq)
      deallocate(rpattern_shum)
      deallocate(vfact_shum)
    endif
    if (doskeb) then
      call get_random_pattern_destroy(rpattern_skeb,nskeb)
      deallocate(skeb3du)
      deallocate(skeb3dv)
      deallocate(diss_est)
      deallocate(diss_dc)
      deallocate(rpattern_skeb)
      deallocate(vfact_skeb)
      deallocate(keb)
      deallocate(kea)
    endif
    if (dossst) then
      call get_random_pattern_destroy(rpattern_ssst,nssst)
      deallocate(ssst3d)
      deallocate(rpattern_ssst)
      deallocate(vfact_ssst)
    endif

  end subroutine destroy_stochastic_physics

  subroutine init_sppt(dtau)
    implicit none
    real :: dtau
    integer :: n, k 
                      
    do n=1,size(sppt)
      if (sppt(n) > 0) then
        nsppt=nsppt+1
      else
        exit
      endif
    enddo

    allocate(rpattern_sppt(nsppt))
    do n=1,nsppt
      rpattern_sppt(n)%stdev = sppt(n)
      rpattern_sppt(n)%decortau = sppt_decort(n)
      rpattern_sppt(n)%lenscale = sppt_lscale(n)
      rpattern_sppt(n)%seed = int(sppt_seed(n))
      if (myrank .eq. 0 ) then
        write(6,*)'mod_stochastic_physics : sppt : stdev  ',sppt(n)
        write(6,*)'mod_stochastic_physics : sppt : decort ',sppt_decort(n)
        write(6,*)'mod_stochastic_physics : sppt : lscale ',sppt_lscale(n)
        write(6,*)'mod_stochastic_physics : sppt : seed   ',sppt_seed(n)
      endif
    enddo

    allocate(sppt3d(nxp,lev,my_max))       
    call get_random_pattern_init(rpattern_sppt,nsppt,dtau,.false.)

    ! set up vfact_sppt
    allocate(vfact_sppt(lev))
    do k=1,lev
      if (sl(k) .lt. sppt_sigtop1 .and. sl(k) .gt. sppt_sigtop2) then
         vfact_sppt(k) = (sl(k)-sppt_sigtop2)/(sppt_sigtop1-sppt_sigtop2)
      else if (sl(k) .lt. sppt_sigtop2) then
          vfact_sppt(k) = 0.0
      else
          vfact_sppt(k) = 1.0
      endif
    enddo

    if (sppt_sfclimit) then
    ! vfact_sppt(lev-1)=vfact_sppt(lev-2)*0.5
    ! vfact_sppt(lev)=0.0
      do k=1,lev
        if (sl(k) .lt. sppt_sigbot1 .and. sl(k) .gt. sppt_sigbot2) then
           vfact_sppt(k) = 1. - (sl(k)-sppt_sigbot2)/(sppt_sigbot1-sppt_sigbot2)
        else if (sl(k) .gt. sppt_sigbot1) then
            vfact_sppt(k) = 0.0
        endif
      enddo
    endif

    do k=1,lev
      if (myrank == 0) print 301,'mod_stochastic_physics : k,sl,vfact_sppt',k,sl(k),vfact_sppt(k)
    enddo
    301 format ( A , I4 , 2F20.12 )

  end subroutine init_sppt

  subroutine init_shum(dtau)
    implicit none
    real :: dtau
    integer :: n, k 
      do n=1,size(shum)
        if (shum(n) > 0) then
          nshum=nshum + 1
        else
          exit
        endif
      enddo

      allocate(rpattern_shum(nshum))
      do n=1,nshum
        rpattern_shum(n)%stdev = shum(n)
        rpattern_shum(n)%decortau = shum_decort(n)
        rpattern_shum(n)%lenscale = shum_lscale(n)
        rpattern_shum(n)%seed = int(shum_seed(n))
        if (myrank .eq. 0 ) then
          write(6,*)'mod_stochastic_physics : shum : stdev  ',shum(n)
          write(6,*)'mod_stochastic_physics : shum : decort ',shum_decort(n)
          write(6,*)'mod_stochastic_physics : shum : lscale ',shum_lscale(n)
          write(6,*)'mod_stochastic_physics : shum : seed   ',shum_seed(n)
        endif
      enddo

      allocate(shum3d(nxp,lev,my_max))       
      allocate(shum3d_dq(nxp,lev,my_max))       
      shum3d_dq = 0.
      call get_random_pattern_init(rpattern_shum,nshum,dtau,.false.)

      allocate(vfact_shum(lev))
      do k=1,lev
         vfact_shum(k) = exp((sl(k)-1.)/shum_sigefold)
         if (sl(k).LT. 2*shum_sigefold) then
            vfact_shum(k)=0.0
         endif
        if (myrank == 0) print 301,'mod_stochastic_physics : k,sl,vfact_shum',k,sl(k),vfact_shum(k)
      enddo
    301 format ( A , I4 , 2F20.12 )
  end subroutine init_shum

  subroutine init_skeb(dtau)
    implicit none
    real :: dtau
    integer :: n, k, k2
    real, allocatable,dimension(:) :: skeb_vloc
                      
    do n=1,size(skeb)
      if (skeb(n) > 0) then
        nskeb=nskeb+1
      else
        exit
      endif
    enddo

    allocate(rpattern_skeb(nskeb))
    do n=1,nskeb
      if ( skebnorm .eq. 0 ) then ! stream function norm
        rpattern_skeb(n)%stdev = skeb(n)*1.111e3*sqrt(dtau)
      endif
      if ( skebnorm .eq. 1 ) then ! kinectic energy function norm
        rpattern_skeb(n)%stdev = skeb(n)*0.00222e3*sqrt(dtau)
      endif
      if ( skebnorm .eq. 2 ) then ! vorticity function norm
       rpattern_skeb(n)%stdev = skeb(n)*1.111e-9*sqrt(dtau)
      endif
      rpattern_skeb(n)%decortau = skeb_decort(n)
      rpattern_skeb(n)%lenscale = skeb_lscale(n)
      rpattern_skeb(n)%seed = int(skeb_seed(n))
      if (myrank .eq. 0 ) then
        write(6,*)'mod_stochastic_physics : skeb : stdev  ',skeb(n)
        write(6,*)'mod_stochastic_physics : skeb : decort ',skeb_decort(n)
        write(6,*)'mod_stochastic_physics : skeb : lscale ',skeb_lscale(n)
        write(6,*)'mod_stochastic_physics : skeb : seed   ',skeb_seed(n)
      endif
    enddo

    allocate(skeb3du(nxp,lev,my_max))
    allocate(skeb3dv(nxp,lev,my_max))
    allocate(diss_est(nxp,lev,my_max))
    allocate(diss_dc(nxp,lev,my_max))
    allocate(keb(nxp,lev,my_max))
    allocate(kea(nxp,lev,my_max))



    diss_est=0.
    diss_dc=0.
    kea=0.
    keb=0.

    ! for 3 time level scheme to keep skeblevs same and make sure can be
    ! reproduced when restart ( restart function not ready yet ).
    if ( first_call ) then
      skeblevs=nint(rpattern_skeb(1)%decortau/(2.*dtau)*skeb_vdof)
    else
      skeblevs=nint(rpattern_skeb(1)%decortau/dtau*skeb_vdof)
    endif

    call get_random_pattern_init(rpattern_skeb,nskeb,dtau,.true.)

    ! set up vfact_skeb
    allocate(vfact_skeb(lev))
    allocate(skeb_vloc(skeblevs))
    allocate(skeb_vwts(lev,2))
    allocate(skeb_vpts(lev,2))
    do k=1,lev
      if (sl(k) .lt. skeb_sigtop1 .and. sl(k) .gt. skeb_sigtop2) then
         vfact_skeb(k) = (sl(k)-skeb_sigtop2)/(skeb_sigtop1-skeb_sigtop2)
      else if (sl(k) .lt. skeb_sigtop2) then
          vfact_skeb(k) = 0.0
      else
          vfact_skeb(k) = 1.0
      endif
    enddo

    do k=1,lev
      if (myrank == 0) print 301,'mod_stochastic_physics : k,sl,vfact_skeb',k,sl(k),vfact_skeb(k)
    enddo
    301 format ( A , I4 , 2F20.12 )
    ! calculate vertical interpolation weights
    do k=1,skeblevs
      skeb_vloc(k)=sl(lev)-real(skeblevs-k)/real(skeblevs-1.0)*(sl(lev)-sl(1))
    enddo
    ! surface
    skeb_vwts(lev,2)=0.
    skeb_vpts(lev,1)=skeblevs-2
    ! top
    skeb_vwts(1,2)=1.
    skeb_vpts(1,1)=1
    ! internal
    do k=2,lev-1
      do k2=1,skeblevs-1
        if (sl(k) .LE. skeb_vloc(k2+1) .AND. sl(k) .GT. skeb_vloc(k2)) then
          skeb_vpts(k,1)=k2
          skeb_vwts(k,2)=(sl(k)-skeb_vloc(k2))/(skeb_vloc(k2+1)-skeb_vloc(k2))
        endif
      enddo
    enddo
    deallocate(skeb_vloc)
    skeb_vwts(:,1)=1.0-skeb_vwts(:,2)
    skeb_vpts(:,2)=skeb_vpts(:,1)+1
    if (myrank .eq. 0) then
      do k=1,lev
        print*,'skeb vpts ',skeb_vpts(k,1),skeb_vwts(k,2)
      enddo
    endif

  end subroutine init_skeb

  subroutine init_ssst(dtau)
    implicit none
    real :: dtau
    integer :: n, k 
                      
    do n=1,size(ssst)
      if (ssst(n) > 0) then
        nssst=nssst+1
      else
        exit
      endif
    enddo

    allocate(rpattern_ssst(nssst))
    do n=1,nssst
      rpattern_ssst(n)%stdev    = ssst(n)
      rpattern_ssst(n)%decortau = ssst_decort(n)
      rpattern_ssst(n)%lenscale = ssst_lscale(n)
      rpattern_ssst(n)%seed = int(ssst_seed(n))
      if (myrank .eq. 0 ) then
        write(6,*)'mod_stochastic_physics : ssst : stdev  ',ssst(n)
        write(6,*)'mod_stochastic_physics : ssst : decort ',ssst_decort(n)
        write(6,*)'mod_stochastic_physics : ssst : lscale ',ssst_lscale(n)
        write(6,*)'mod_stochastic_physics : ssst : seed   ',ssst_seed(n)
      endif
    enddo

    allocate(ssst3d(nxp,1,my_max))
    call get_random_pattern_init(rpattern_ssst,nssst,dtau,.false.)

    allocate(vfact_ssst(1))
    vfact_ssst = 1.

  end subroutine init_ssst

  subroutine get_random_pattern_init(rpattern,nscale,dt,skebrun)
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
    integer :: timearray(3),iseed
    integer :: n, k, nscale, ncx, ml, ms, ns, i, j
    real :: rerth, pi, var, correLsq, rkT, rnn1
    type(random_pattern), intent(inout) :: rpattern(nscale)
    integer :: irand
    real :: dt
    real(kind=RTYPE), allocatable :: noise(:,:)
    integer(8) count, count_rate, count_max, count_trunc
    integer(8) :: iscale = 10000000000_8
    integer :: count4 
    logical :: skebrun
 
    rerth = 6.3712e+6      ! radius of earth (m)
    pi = 4.*atan(1.)
!    radsq = rerth*rerth

    do n=1,nscale
      ncx = 2.*pi*rerth/rpattern(n)%lenscale
   !  rpattern(n)%jtrun = 2*((1+(ncx-1)/3)/2) 
      rpattern(n)%jtrun = 2*((1+(4*ncx-1)/4)/2)
      rpattern(n)%mlmax = rpattern(n)%jtrun*(rpattern(n)%jtrun+1)/2

      allocate(rpattern(n)%n2d(nxp,my_max))
      if ( doskeb ) then
        allocate(rpattern(n)%n2du(nxp,my_max,skeblevs))
        allocate(rpattern(n)%n2dv(nxp,my_max,skeblevs))
        allocate(rpattern(n)%kenorm(rpattern(n)%mlmax,2))
        rpattern(n)%n2du(:,:,:)=0.0
        rpattern(n)%n2dv(:,:,:)=0.0
      endif
      allocate(rpattern(n)%spec(rpattern(n)%mlmax,2))
      allocate(rpattern(n)%varspec(rpattern(n)%mlmax))
      allocate(rpattern(n)%msort(rpattern(n)%mlmax))
      allocate(rpattern(n)%lsort(rpattern(n)%mlmax))
      allocate(rpattern(n)%mlsort(rpattern(n)%jtrun,rpattern(n)%jtrun))
      allocate(noise(rpattern(n)%mlmax,2))

      ! Real random seeds
      if (myrank.eq.0) then     
        if(.not. ncep_seeds) then 
          call itime(timearray)
          iseed=irand(0)
          count4=irand(timearray(1)+(i*11467-iseed)*timearray(2) &
              +(iseed-i*23)*timearray(3))
        else
          call system_clock(count, count_rate, count_max)
          count_trunc = iscale*(count/iscale)
          count4 = count - count_trunc
        endif
      endif
!     call mpe_bcast(count4,1,0,mpe_double) 
      call mpe_bcast(count4,1,0,mpe_integer)
      if (rpattern(n)%seed == -999 ) then
        rpattern(n)%seed = count4
      endif
      if (myrank.eq.0) write(6,*)'scale',n,' : using seed :',rpattern(n)%seed
      call random_setseed(rpattern(n)%seed,rpattern(n)%rstate)

      ! horizontal decorrelation 
      correLsq = rpattern(n)%lenscale*rpattern(n)%lenscale
      rkT = 0.25*correLsq/radsq
      if (myrank.eq.0)  write(6,*)'mod_stochastic_physics : n,rkT = ',n,rkT
 
      ! time decorrelation 
      rpattern(n)%phi = exp(-dt/rpattern(n)%decortau)  
 
      call sortml(rpattern(n)%jtrun,rpattern(n)%mlmax, &
                  rpattern(n)%msort,rpattern(n)%lsort,rpattern(n)%mlsort)

      noise = 0.
      do ml=1,rpattern(n)%mlmax
        noise(ml,1) = 1./sqrt(float(2*(rpattern(n)%lsort(ml)-1))+1)
        noise(ml,2) = 1./sqrt(float(2*(rpattern(n)%lsort(ml)-1))+1)
        if (rpattern(n)%msort(ml) .eq. 1) then
          noise(ml,1) = sqrt(2.)/sqrt(float(2*(rpattern(n)%lsort(ml)-1))+1)
          noise(ml,2) = 0.
        endif
      enddo

      noise(1,1) = 0. 
      noise(1,2) = 0.
      noise = noise*sqrt(1./float(rpattern(n)%jtrun))

      ! set up the amplitude of noise
      do ml=1,rpattern(n)%mlmax
        rpattern(n)%varspec(ml) = sqrt(float(rpattern(n)%jtrun) & 
            * exp(-rkT*float(rpattern(n)%lsort(ml))*(float(rpattern(n)%lsort(ml)-1))))
      enddo

      noise(:,1) = noise(:,1)*rpattern(n)%varspec
      noise(:,2) = noise(:,2)*rpattern(n)%varspec

      ! get specral variance
      var=0.
      do ml=1,rpattern(n)%mlmax
        if (rpattern(n)%msort(ml) .ne. 1) then
          var = var + (noise(ml,1)**2 + noise(ml,2)**2)
        else
          var = var + 0.5*(noise(ml,1)**2 + noise(ml,2)**2)
        endif
      enddo
      rpattern(n)%varspec = rpattern(n)%varspec / sqrt(var)

#ifdef VERBOSE
      if (myrank.eq.0) write(6,*)'total variance =',var
      if (myrank.eq.0) write(6,*)'rpattern(n)%varspec =',rpattern(n)%varspec
#endif
      ! initialize spectrum coefficient 
      noise = 0. 
      call get_noise(rpattern(n),noise) 
      do ml=1,rpattern(n)%mlmax
        rpattern(n)%spec(ml,1) = rpattern(n)%stdev*rpattern(n)%varspec(ml)*noise(ml,1)
        rpattern(n)%spec(ml,2) = rpattern(n)%stdev*rpattern(n)%varspec(ml)*noise(ml,2)
      enddo
      rpattern(n)%spec(1,1) = 0.
      rpattern(n)%spec(1,2) = 0.

      deallocate(noise)
!
      if ( skebrun ) then
        rpattern(n)%kenorm(:,:) = 1.
        if ( skebnorm .eq. 0 ) then
          do j=1,rpattern(n)%jtrun
            do i=j,rpattern(n)%jtrun
              ml=rpattern(n)%mlsort(j,i)
              rnn1=float(i*(i+1))
              rpattern(n)%kenorm(ml,1) = rnn1/radsq
              rpattern(n)%kenorm(ml,2) = rnn1/radsq
            enddo
          enddo
          if (myrank .eq. 0 ) print *,'using streamfunction ',  &
             maxval(rpattern(n)%kenorm(:,1)),minval(rpattern(n)%kenorm(:,1))
        endif
        if ( skebnorm .eq. 1 ) then
          do j=1,rpattern(n)%jtrun
            do i=j,rpattern(n)%jtrun
              ml=rpattern(n)%mlsort(j,i)
              rnn1=float(i*(i+1))
              rpattern(n)%kenorm(ml,1) = sqrt(rnn1)/rerth
              rpattern(n)%kenorm(ml,2) = sqrt(rnn1)/rerth
            enddo
          enddo
          if (myrank .eq. 0 ) print *,'using kenorm ',  &
             maxval(rpattern(n)%kenorm(:,1)),minval(rpattern(n)%kenorm(:,1))
        endif
      endif
    enddo

    if (myrank.eq.0) then
      write(6,*)'mod_stochastic_physics : dt    = ',dt
      write(6,*)'mod_stochastic_physics : stdev = ',rpattern%stdev
      write(6,*)'mod_stochastic_physics : seed  = ',rpattern%seed
      write(6,*)'mod_stochastic_physics : jtrun = ',rpattern%jtrun
      write(6,*)'mod_stochastic_physics : mlmax = ',rpattern%mlmax
      write(6,*)'mod_stochastic_physics : tau   = ',rpattern%decortau
      write(6,*)'mod_stochastic_physics : phi   = ',rpattern%phi
    endif

  end subroutine get_random_pattern_init

  subroutine get_random_pattern_destroy(rpattern,nscale)
    implicit none
    integer :: n, nscale
    type(random_pattern), intent(inout) :: rpattern(nscale)

    do n=1,nscale
      deallocate(rpattern(n)%n2d)
      if ( doskeb ) then
        deallocate(rpattern(n)%n2du)
        deallocate(rpattern(n)%n2dv)
        deallocate(rpattern(n)%kenorm)
      endif
      deallocate(rpattern(n)%spec)
      deallocate(rpattern(n)%varspec)
      deallocate(rpattern(n)%msort)
      deallocate(rpattern(n)%lsort)
      deallocate(rpattern(n)%mlsort)
    enddo

  end subroutine get_random_pattern_destroy

  subroutine get_random_pattern_run(rpattern,nscale)
    implicit none
    integer :: nscale
    integer :: n, ii, i, jj, j, k, nxj
    type(random_pattern), intent(inout) :: rpattern(nscale)

    do n=1,nscale
      call gen_random_pattern_2d(rpattern(n))
    enddo

  end subroutine get_random_pattern_run

  subroutine get_random_pattern_run_vect(rpattern,nscale,k)
    implicit none
    integer :: nscale
    integer :: n, ii, i, jj, j, k, nxj, k2
    type(random_pattern), intent(inout) :: rpattern(nscale)

    do n=1,nscale
      call gen_random_pattern_2d_vect(rpattern(n),k)
    enddo

  end subroutine get_random_pattern_run_vect

  subroutine get_stochy_physics(rpattern,nscale,nlev,vfact,n3d)
!------------------------------------------------------------------------! 
!  purpose: To generate 3D SPPT strucutre
!  output: n3d
!------------------------------------------------------------------------! 
    implicit none
    integer :: n, ii, i, jj, j, k, nxj, nlev
    integer, intent(in) :: nscale
    real, intent(in) :: vfact(nlev) 
    type(random_pattern), intent(inout) :: rpattern(nscale)
    real(kind=RTYPE), intent(  out) :: n3d(nxp,nlev,my_max) 
 
    n3d = 0.
    do n=1,nscale
      do jj=1,jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do k=1,nlev
          do i=1,nxj
            n3d(i,k,jj)=n3d(i,k,jj)+rpattern(n)%n2d(i,jj)*vfact(k)
          enddo
        enddo
      enddo
    enddo

  end subroutine get_stochy_physics

  subroutine get_stochy_physics_vect(rpattern,nscale,nlev,vfact,n3du,n3dv)
!------------------------------------------------------------------------! 
!  purpose: To generate 3D SKEB strucutre
!  output: n3du n3dv
!------------------------------------------------------------------------! 
    implicit none
    integer :: n, ii, i, jj, j, k, nxj, nlev
    integer, intent(in) :: nscale
    real, intent(in) :: vfact(nlev) 
    type(random_pattern), intent(inout) :: rpattern(nscale)
    real(kind=RTYPE), intent(  out) :: n3du(nxp,nlev,my_max),n3dv(nxp,nlev,my_max)

 
    n3du = 0.
    n3dv = 0.
    do n=1,nscale
      do jj=1,jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do k=1,nlev
          do i=1,nxj
            n3du(i,k,jj)=n3du(i,k,jj)+vfact(k)                         &
               *(skeb_vwts(k,1)*rpattern(n)%n2du(i,jj,skeb_vpts(k,1))  &
               +skeb_vwts(k,2)*rpattern(n)%n2du(i,jj,skeb_vpts(k,2)))
            n3dv(i,k,jj)=n3dv(i,k,jj)+vfact(k)                         &
               *(skeb_vwts(k,1)*rpattern(n)%n2dv(i,jj,skeb_vpts(k,1))  &
               +skeb_vwts(k,2)*rpattern(n)%n2dv(i,jj,skeb_vpts(k,2)))
          enddo
        enddo
      enddo
    enddo

  end subroutine get_stochy_physics_vect
  
  subroutine get_noise(rpattern,noise)
    implicit none
    integer :: ml, ns, ms
    type(random_pattern), intent(inout) :: rpattern
    real(kind=RTYPE), intent(out) :: noise(rpattern%mlmax,2)
    real :: noise_gauss(2*rpattern%mlmax)
    real :: ave, var, std

    ! get white noise  with a gaussian (normal) distribution
    call random_gauss(noise_gauss,rpattern%rstate)
    noise_gauss(1) = 0.; noise_gauss(rpattern%mlmax+1) = 0.
    noise_gauss = noise_gauss*sqrt(1./float(rpattern%jtrun-1))

    noise=0.
    do ml=1,rpattern%mlmax
      ! set up to red noise 
      noise(ml,1) = noise_gauss(ml)/sqrt(float(2*(rpattern%lsort(ml)-1)+1))
      noise(ml,2) = noise_gauss(rpattern%mlmax+ml)/sqrt(float(2*(rpattern%lsort(ml)-1)+1))
      ! zero out the imagenary part when m(zonal wavenumber) equal to 1
      if (rpattern%msort(ml) .eq. 1) then
        noise(ml,1) = sqrt(2.)*noise(ml,1)
        noise(ml,2) = 0.
      endif
    enddo
  end subroutine get_noise

  subroutine gen_random_pattern_2d(rpattern)

    implicit none
    type(random_pattern), intent(inout) :: rpattern
    !real(kind=RTYPE), intent(out) :: sppt2d(nxp,my_max)
    integer :: ml, ns, ms
    real(kind=RTYPE), allocatable :: noise(:,:)
    real(kind=RTYPE), allocatable :: specp(:,:,:),bufr2d(:,:,:)

    allocate(bufr2d(jtrun,jtmax*nsizey,2)) 
    allocate(specp(jtrun,jtmax,2)) 

    if ( col_rank .eq. 0 ) then
      ! get noise
      allocate(noise(rpattern%mlmax,2)) 
      call get_noise(rpattern,noise) 
 
      !  radom pattern advance with first order AR
      rpattern%spec(:,1) = rpattern%phi*rpattern%spec(:,1) + & 
            sqrt(1.-rpattern%phi**2.)*rpattern%stdev*rpattern%varspec*noise(:,1)
      rpattern%spec(:,2) = rpattern%phi*rpattern%spec(:,2) + & 
            sqrt(1.-rpattern%phi**2.)*rpattern%stdev*rpattern%varspec*noise(:,2)

      ! ready for mpi_scatter random pattern
      call spectrun_inp2d(rpattern%jtrun,jtrun,jtmax         &
                         ,rpattern%mlsort,nsizey             &
                         ,rpattern%spec,bufr2d)

      deallocate(noise)
    endif

    !  mpi_scatter random pattern from root
    call mpe_scatter_sppt(bufr2d,specp,2*jtrun*jtmax,nsizey)

    ! transform spectral to physical space 
    call transr1(jtrun,jtmax,nx,my,my_max,poly,specp,rpattern%n2d,nsizey)

    deallocate(bufr2d)
    deallocate(specp)

  end subroutine gen_random_pattern_2d

  subroutine gen_random_pattern_2d_vect(rpattern,k)

    implicit none
    type(random_pattern), intent(inout) :: rpattern
    integer :: ml, ns, ms, k, i, j, jj, nxj
    real    :: xx
    real(kind=RTYPE), allocatable :: specpv(:,:,:),specpd(:,:,:)
    real(kind=RTYPE), allocatable :: bufr2d(:,:,:) ,noise(:,:) ,specf(:,:)

    allocate(bufr2d(jtrun,jtmax*nsizey,2)) 
    allocate(specpd(jtrun,jtmax,2)) !divergence
    allocate(specpv(jtrun,jtmax,2)) !vorticity

    specpd = 0.0
    specpv = 0.0

    if ( col_rank .eq. 0 ) then
      ! get noise
      allocate(noise(rpattern%mlmax,2)) 
      allocate(specf(rpattern%mlmax,2)) 
      call get_noise(rpattern,noise) 
 
      !  radom pattern advance with first order AR
      rpattern%spec(:,1) = rpattern%phi*rpattern%spec(:,1) + & 
            sqrt(1.-rpattern%phi**2.)*rpattern%stdev*rpattern%varspec*noise(:,1)
      rpattern%spec(:,2) = rpattern%phi*rpattern%spec(:,2) + & 
            sqrt(1.-rpattern%phi**2.)*rpattern%stdev*rpattern%varspec*noise(:,2)

      specf(:,1) = rpattern%kenorm(:,1) * rpattern%spec(:,1)
      specf(:,2) = rpattern%kenorm(:,2) * rpattern%spec(:,2)

      ! ready for mpi_scatter random pattern
      call spectrun_inp2d(rpattern%jtrun,jtrun,jtmax         &
                         ,rpattern%mlsort,nsizey             &
                         ,specf,bufr2d)

      deallocate(noise)
      deallocate(specf)
    endif

    !  mpi_scatter random pattern from root
    call mpe_scatter_sppt(bufr2d,specpv,2*jtrun*jtmax,nsizey)

    !  spectral transform for velocity components
    call tranuv1 (jtrun,jtmax,nx,my,my_max,1,onocos,wcfac,wdfac &
                , poly,dpoly,specpv,specpd,rpattern%n2du(1,1,k),rpattern%n2dv(1,1,k),nsizey)
    do jj =1,jlistnum
      j=jlist1(jj)
      nxj=nxdef_2d(j)
      xx=rad/cosl(j)
      do i =1,nxj
        rpattern%n2du(i,jj,k) = rpattern%n2du(i,jj,k) * xx
        rpattern%n2dv(i,jj,k) = rpattern%n2dv(i,jj,k) * xx
      enddo
    enddo

    deallocate(bufr2d)
    deallocate(specpd)
    deallocate(specpv)

  end subroutine gen_random_pattern_2d_vect

  SUBROUTINE avevar_sppt(data,n,ave,var,std)
      INTEGER n
      REAL ave,var,data(n)
      INTEGER j
      REAL s,ep,std
      ave=0.0
      do j=1,n
        ave=ave+data(j)
      enddo
      ave=ave/n
      var=0.0
      ep=0.0
      do j=1,n
        s=data(j)-ave
        ep=ep+s
        var=var+s*s
      enddo
      var=(var-ep**2/n)/n
      std=sqrt(var)
      
  end SUBROUTINE avevar_sppt
!
!
  SUBROUTINE avevar_sppt2d(data2d,n,m,ave,var,std)
    INTEGER :: n,m,nmdim
    REAL :: ave,var,data(n*m)
    REAL(kind=RTYPE) :: data2d(n,m)
    INTEGER :: i,j
    REAL :: s,ep,std

    nmdim=0
    do j=1,m
      do i=1,n
      nmdim=nmdim+1
      data(nmdim)=data2d(i,j)
      enddo
    enddo
    
    ave=0.0
    do j=1,nmdim
      ave=ave+data(j)
    enddo
    ave=ave/nmdim
    var=0.0
    ep=0.0
    do j=1,nmdim
      s=data(j)-ave
      ep=ep+s
      var=var+s*s
    enddo
!   var=(var-ep**2/n)/(n-1)    ! sample numners < 30
    var=(var-ep**2/nmdim)/nmdim
    std=sqrt(var) 
  end SUBROUTINE avevar_sppt2d
!
  subroutine spptout(tau)
    implicit none
    integer      :: i, j, k, jj, nxj, ihead, n
    integer      :: nxmy4
    real         :: tau
    real*4       :: glob4(nx,my)
    real(kind=RTYPE) :: glob(nx,my),temp(nxp,my_max)

    ihead=15
    nxmy4=nx*my*4
    if ( myrank .eq. 0 ) then
      open(ihead,file='sppt.dat',access='direct',form='unformatted'   &
                ,recl=nxmy4,status='unknown',convert='big_endian')
    endif

    do n=1,nsppt
      call unify_reduceintp(nx,my,my_max,rpattern_sppt(n)%n2d,glob)   
      if ( myrank .eq. 0 ) then
        glob4=glob
        write(ihead,rec=recnsppt) glob4
        recnsppt=recnsppt+1
      endif
    enddo
    do k=lev,1,-1
      do jj=1,jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i=1,nxj
        temp(i,jj)=sppt3d(i,k,jj)
        enddo
      enddo
      call unify_reduceintp(nx,my,my_max,temp,glob)   
      if ( myrank.eq.0 ) then
        glob4=glob
        write(ihead,rec=recnsppt) glob4
        recnsppt=recnsppt+1
      endif
    enddo

    if ( myrank.eq.0 ) close(ihead)

    !creating ctl file
    call spptctl(tau)
 
  end subroutine spptout
!----
  subroutine shumout(tau)
    implicit none
    integer      :: i, j, k, jj, nxj, ihead, n
    integer      :: nxmy4
    real         :: tau
    real*4       :: glob4(nx,my)
    real(kind=RTYPE) :: glob(nx,my),temp(nxp,my_max)

    ihead=15
    nxmy4=nx*my*4
    if ( myrank .eq. 0 ) then
      open(ihead,file='shum.dat',access='direct',form='unformatted'   &
                ,recl=nxmy4,status='unknown',convert='big_endian')
    endif

    do n=1,nshum
      call unify_reduceintp(nx,my,my_max,rpattern_shum(n)%n2d,glob)
      if ( myrank .eq. 0 ) then
        glob4=glob
        write(ihead,rec=recnshum) glob4
        recnshum=recnshum+1
      endif
    enddo
    do k=lev,1,-1
      do jj=1,jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i=1,nxj
        temp(i,jj)=shum3d_dq(i,k,jj)
        enddo
      enddo
      call unify_reduceintp(nx,my,my_max,temp,glob)
      if ( myrank.eq.0 ) then
        glob4=glob
        write(ihead,rec=recnshum) glob4
        recnshum=recnshum+1
      endif
    enddo

    if ( myrank.eq.0 ) close(ihead)

    !creating ctl file
    call shumctl(tau)

  end subroutine shumout
!----
  subroutine skebest(um,vm)
    use grid,  only : ut,vt
    implicit none

    integer :: n, ii, i, jj, j, k, nxj
    real    :: xx, axx
    real(kind=RTYPE) :: spectmp(levp,2,jtrun,jtmax),                 &
                        cc(nx+2,levp,1,my_max),                      &
                        um(nxp,lev,my_max),vm(nxp,lev,my_max),       &
                        dummy,temp

    ! estimate the dissipation of kinectic energy for SKEB
        do jj =1,jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          xx=radsq*onocos(j)
          do k=1,lev
            do i=1,nxj
              diss_est(i,k,jj)=((um(i,k,jj)*ut(i,k,jj)               &
                               +vm(i,k,jj)*vt(i,k,jj))               &
                               +0.5*(um(i,k,jj)**2.+vm(i,k,jj)**2.)) &
                               *xx
            enddo
          enddo
        enddo
        call joinrs(cc,diss_est,dummy,dummy,dummy,nx,my_max,lev      &
                   ,jlistnum,1,1)
        call tranrs(jtrun,jtmax,nx,my,my_max,levp,poly,weight,cc     &
                   ,spectmp,1,nsizey)
        ! apply spectral filter of the dissipation of kinetic energy
        call filter_skeb(jtrun,jtmax,levp,spectmp,skebfilt)
        call transr(jtrun,jtmax,nx,my,my_max,levp,poly,spectmp,cc,1,nsizey)
        call ujoinsr(cc,diss_est,dummy,dummy,dummy,nx,my_max,lev,jlistnum,1,1)
!
!        if ( myrank .eq. 0 ) print *,'intgrt: diss_est(1,72,1)=',diss_est(1,72,1)
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          xx=rad/cosl(j)
          axx=cosl(j)/rad
          do k = 1, lev
            do i = 1, nxj
              temp=diss_est(i,k,jj)
              if(doskeb_dc) then
                temp=diss_est(i,k,jj)+diss_dc(i,k,jj)
              endif
              !change virtual wind to real wind
              ut(i,k,jj)=ut(i,k,jj)*xx
              vt(i,k,jj)=vt(i,k,jj)*xx
              keb(i,k,jj)=0.5*(ut(i,k,jj)**2.+vt(i,k,jj)**2.)
!              ut(i,k,jj)=ut(i,k,jj)+skeb3du(i,k,jj)*diss_est(i,k,jj)
!              vt(i,k,jj)=vt(i,k,jj)+skeb3dv(i,k,jj)*diss_est(i,k,jj)
              ut(i,k,jj)=ut(i,k,jj)+skeb3du(i,k,jj)*temp
              vt(i,k,jj)=vt(i,k,jj)+skeb3dv(i,k,jj)*temp
              kea(i,k,jj)=0.5*(ut(i,k,jj)**2.+vt(i,k,jj)**2.)
              !change back to virtual wind
              ut(i,k,jj)=ut(i,k,jj)*axx
              vt(i,k,jj)=vt(i,k,jj)*axx
            enddo
          enddo
        enddo
!        if ( myrank .eq. 0 ) print *,'intgrt: keb(1,72,1)=',keb(1,72,1)
!        if ( myrank .eq. 0 ) print *,'intgrt: kea(1,72,1)=',kea(1,72,1)

  end subroutine skebest
!----
  subroutine skebout(tau)
    implicit none
    integer      :: i, j, k, jj, nxj, ihead, n
    integer      :: nxmy4
    real         :: tau,xx
    real(kind=RTYPE):: glob(nx,my),temp(nxp,my_max)
    real(kind=4) :: glob4(nx,my)


    ihead=15
    nxmy4=nx*my*4
    if ( myrank .eq. 0 ) then
      open(ihead,file='skeb.dat',access='direct',form='unformatted'   &
                ,recl=nxmy4,status='unknown',convert='big_endian')
    endif

    do k=lev,1,-1
      do jj=1,jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i=1,nxj
        temp(i,jj)=skeb3du(i,k,jj)
        enddo
      enddo
      call unify_reduceintp(nx,my,my_max,temp,glob)   
      if ( myrank.eq.0 ) then
        glob4=glob
        write(ihead,rec=recnskeb) glob4
        recnskeb=recnskeb+1
      endif
    enddo

    do k=lev,1,-1
      do jj=1,jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i=1,nxj
        temp(i,jj)=skeb3dv(i,k,jj)
        enddo
      enddo
      call unify_reduceintp(nx,my,my_max,temp,glob)
      if ( myrank.eq.0 ) then
        glob4=glob
        write(ihead,rec=recnskeb) glob4
        recnskeb=recnskeb+1
      endif
    enddo

    do k=lev,1,-1
      do jj=1,jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i=1,nxj
          temp(i,jj)=diss_est(i,k,jj)
        enddo
      enddo
      call unify_reduceintp(nx,my,my_max,temp,glob)
      if ( myrank.eq.0 ) then
        glob4=glob
        write(ihead,rec=recnskeb) glob4
        recnskeb=recnskeb+1
      endif
    enddo

    do k=lev,1,-1
      do jj=1,jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i=1,nxj
          temp(i,jj)=keb(i,k,jj)
        enddo
      enddo
      call unify_reduceintp(nx,my,my_max,temp,glob)
      if ( myrank.eq.0 ) then
        glob4=glob
        write(ihead,rec=recnskeb) glob4
        recnskeb=recnskeb+1
      endif
    enddo

    do k=lev,1,-1
      do jj=1,jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i=1,nxj
          temp(i,jj)=kea(i,k,jj)
        enddo
      enddo
      call unify_reduceintp(nx,my,my_max,temp,glob)
      if ( myrank.eq.0 ) then
        glob4=glob
        write(ihead,rec=recnskeb) glob4
        recnskeb=recnskeb+1
      endif
    enddo

    do k=lev,1,-1
      do jj=1,jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i=1,nxj
          temp(i,jj)=diss_dc(i,k,jj)
        enddo
      enddo
      call unify_reduceintp(nx,my,my_max,temp,glob)
      if ( myrank.eq.0 ) then
        glob4=glob
        write(ihead,rec=recnskeb) glob4
        recnskeb=recnskeb+1
      endif
    enddo
    if ( myrank.eq.0 ) close(ihead)

    !creating ctl file
    call skebctl(tau)
 
  end subroutine skebout
!
!----
  subroutine spptctl(tau)
!
    use const, only : idtg,sinl,sigma
    use rank, only : myrank

    integer :: nxj,j,k,itau,ch,iter,remd,js,je,kk,lev1, mn, n
    real :: pi, r2d, dlon, tau
    real, allocatable :: mlat(:),prsl(:)
    character(len=30) ::  forydef,forzdef
    
    character yy*4,dd*2,hh*2,mm*2,dtg*12
    character*3 mon(12)
    logical jrem

    data mon/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep' &
             ,'Oct','Nov','Dec'/

    allocate(mlat(my),prsl(lev))
   
    itau=tau
    ch=10
    lev1=1
    pi=4.0*atan(1.0)
    r2d=180./pi
    dlon=360./float(nx)

    do j=1,my
      mlat(j)=asin(sinl(j))*r2d
    enddo

    do k=1,lev
      kk=lev-k+1
      prsl(kk)=sigma(k,2)+sigma(k+1,2)
      prsl(kk)=prsl(kk)+(sigma(k,1)+sigma(k+1,1))*1000.
      prsl(kk)=0.5*prsl(kk)
    enddo

      write(forydef,8) my
      write(forzdef,9) lev

      write(dtg,'(I12)') idtg
      read(dtg,'(A4,I2,A2,A2,A2)') yy,mn,dd,hh,mm

      if ( myrank .eq. 0 ) then
      OPEN(UNIT=ch, FILE='sppt.ctl', STATUS='UNKNOWN'             &
         , ACCESS='SEQUENTIAL')

      write(ch,'(A14)') 'dset ^sppt.dat'
      write(ch,'(A18)') 'options big_endian'
      write(ch,'(A12)') 'undef -999.0'
      write(ch,12) 'ydef' ,my, 'levels'
      iter=my/8
      remd=mod(my,8)
      jrem=(remd .eq. 0)
      write(forydef,8) remd
      do j=1,iter
       js=1+8*(j-1)
       je=js+7
       write(ch,10) mlat(js:je)
      enddo
      if ( .not. jrem ) write(10,forydef) mlat(je+1:my)
       
      write(ch,13) 'xdef'  ,nx, 'linear 0.0',dlon
      write(ch,14) 'tdef',itau, 'linear',hh,'Z',dd,mon(mn),yy,'1hr'
      write(ch,12) 'zdef' ,lev, 'levels '
      iter=lev/8
      remd=mod(lev,8)
      jrem=(remd .eq. 0)
      write(forzdef,9) remd
      do j=1,iter
       js=1+8*(j-1)
       je=js+7
       write(ch,11) prsl(js:je)
      enddo
      if ( .not. jrem ) write(ch,forzdef) prsl(je+1:my)
      write(ch,'(A4,1X,I2)') 'vars',nsppt+1
      do n=1,nsppt
        write(ch,15) 'scale',n ,lev1,'99',rpattern_sppt(n)%lenscale/1000.,'km 2D Random Pattern'
      enddo
      write(ch,16) 'sppt  '   , lev,'99','3D Random Pattern'
      write(ch,'(A7)') 'endvars'

      close(ch)
      endif

8     format("(",I4,"(2x,F11.7))")
9     format("(",I4,"(2x,F10.5))")
10    format(8(2x,F11.7))
11    format(8(2x,F10.5))
12    format(A4,1X,I4,1X,A6)
13    format(A4,1X,I4,1X,A10,1X,F10.7)
14    format(A4,1X,I4,1X,A6,1X,A2,A1,A2,A3,A4,1X,A3)
15    format(A5,I1,3X,I3,1X,A2,1X,F7.1,A20)
16    format(A6,3X,I3,1X,A2,1X,A20)

    deallocate(mlat,prsl)
  end subroutine spptctl
!
!----
!----
  subroutine shumctl(tau)
!
    use const, only : idtg,sinl,sigma
    use rank, only : myrank

    integer :: nxj,j,k,itau,ch,iter,remd,js,je,kk,lev1, mn, n
    real :: pi, r2d, dlon, tau
    real, allocatable :: mlat(:),prsl(:)
    character(len=30) ::  forydef,forzdef

    character yy*4,dd*2,hh*2,mm*2,dtg*12
    character*3 mon(12)
    logical jrem

    data mon/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep' &
             ,'Oct','Nov','Dec'/

    allocate(mlat(my),prsl(lev))

    itau=tau
    ch=10
    lev1=1
    pi=4.0*atan(1.0)
    r2d=180./pi
    dlon=360./float(nx)

    do j=1,my
      mlat(j)=asin(sinl(j))*r2d
    enddo

    do k=1,lev
      kk=lev-k+1
      prsl(kk)=sigma(k,2)+sigma(k+1,2)
      prsl(kk)=prsl(kk)+(sigma(k,1)+sigma(k+1,1))*1000.
      prsl(kk)=0.5*prsl(kk)
    enddo

      write(forydef,8) my
      write(forzdef,9) lev

      write(dtg,'(I12)') idtg
      read(dtg,'(A4,I2,A2,A2,A2)') yy,mn,dd,hh,mm

      if ( myrank .eq. 0 ) then
      OPEN(UNIT=ch, FILE='shum.ctl', STATUS='UNKNOWN'             &
         , ACCESS='SEQUENTIAL')

      write(ch,'(A14)') 'dset ^shum.dat'
      write(ch,'(A18)') 'options big_endian'
      write(ch,'(A12)') 'undef -999.0'
      write(ch,12) 'ydef' ,my, 'levels'
      iter=my/8
      remd=mod(my,8)
      jrem=(remd .eq. 0)
      write(forydef,8) remd
      do j=1,iter
       js=1+8*(j-1)
       je=js+7
       write(ch,10) mlat(js:je)
      enddo
      if ( .not. jrem ) write(ch,forydef) mlat(je+1:my)

      write(ch,13) 'xdef'  ,nx, 'linear 0.0',dlon
      write(ch,14) 'tdef',itau, 'linear',hh,'Z',dd,mon(mn),yy,'1hr'
      write(ch,12) 'zdef' ,lev, 'levels '
      iter=lev/8
      remd=mod(lev,8)
      jrem=(remd .eq. 0)
      write(forzdef,9) remd
      do j=1,iter
       js=1+8*(j-1)
       je=js+7
       write(ch,11) prsl(js:je)
      enddo
      if ( .not. jrem ) write(ch,forzdef) prsl(je+1:my)
      write(ch,'(A4,1X,I2)') 'vars',nshum+1
      do n=1,nshum
        write(ch,15) 'scale',n ,lev1,'99',rpattern_shum(n)%lenscale/1000.,'km 2D Random Pattern'
      enddo
      write(ch,16) 'shum3d_dq  '   , lev,'99','3D delta q'
      write(ch,'(A7)') 'endvars'

      close(ch)
      endif

8     format("(",I4,"(2x,F11.7))")
9     format("(",I4,"(2x,F10.5))")
10    format(8(2x,F11.7))
11    format(8(2x,F10.5))
12    format(A4,1X,I4,1X,A6)
13    format(A4,1X,I4,1X,A10,1X,F10.7)
14    format(A4,1X,I4,1X,A6,1X,A2,A1,A2,A3,A4,1X,A3)
15    format(A5,I1,3X,I3,1X,A2,1X,F7.1,A20)
16    format(A6,3X,I3,1X,A2,1X,A20)

    deallocate(mlat,prsl)
  end subroutine shumctl

!----
!----
  subroutine skebctl(tau)
!
    use const, only : idtg,sinl,sigma
    use rank, only : myrank

    integer :: nxj,j,k,itau,ch,iter,remd,js,je,kk,lev1, mn, n
    real :: pi, r2d, dlon, tau
    real, allocatable :: mlat(:),prsl(:)
    character(len=30) ::  forydef,forzdef
    
    character yy*4,dd*2,hh*2,mm*2,dtg*12
    character*3 mon(12)
    logical jrem

    data mon/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep' &
             ,'Oct','Nov','Dec'/

    allocate(mlat(my),prsl(lev))
   
    itau=tau
    ch=10
    lev1=1
    pi=4.0*atan(1.0)
    r2d=180./pi
    dlon=360./float(nx)

    do j=1,my
      mlat(j)=asin(sinl(j))*r2d
    enddo

    do k=1,lev
      kk=lev-k+1
      prsl(kk)=sigma(k,2)+sigma(k+1,2)
      prsl(kk)=prsl(kk)+(sigma(k,1)+sigma(k+1,1))*1000.
      prsl(kk)=0.5*prsl(kk)
    enddo

      write(forydef,8) my
      write(forzdef,9) lev

      write(dtg,'(I12)') idtg
      read(dtg,'(A4,I2,A2,A2,A2)') yy,mn,dd,hh,mm

      if ( myrank .eq. 0 ) then
      OPEN(UNIT=ch, FILE='skeb.ctl', STATUS='UNKNOWN'             &
         , ACCESS='SEQUENTIAL')

      write(ch,'(A14)') 'dset ^skeb.dat'
      write(ch,'(A18)') 'options big_endian'
      write(ch,'(A12)') 'undef -999.0'
      write(ch,12) 'ydef' ,my, 'levels'
      iter=my/8
      remd=mod(my,8)
      jrem=(remd .eq. 0)
      write(forydef,8) remd
      do j=1,iter
       js=1+8*(j-1)
       je=js+7
       write(ch,10) mlat(js:je)
      enddo
      if ( .not. jrem ) write(ch,forydef) mlat(je+1:my)
       
      write(ch,13) 'xdef'  ,nx, 'linear 0.0',dlon
      write(ch,14) 'tdef',itau, 'linear',hh,'Z',dd,mon(mn),yy,'1hr'
      write(ch,12) 'zdef' ,lev, 'levels '
      iter=lev/8
      remd=mod(lev,8)
      jrem=(remd .eq. 0)
      write(forzdef,9) remd
      do j=1,iter
       js=1+8*(j-1)
       je=js+7
       write(ch,11) prsl(js:je)
      enddo
      if ( .not. jrem ) write(ch,forzdef) prsl(je+1:my)
      write(ch,'(A4,1X,I2)') 'vars',6
      write(ch,16) 'skebu  '   , lev,'99','U-dir 3D Random Pattern'
      write(ch,16) 'skebv  '   , lev,'99','V-dir 3D Random Pattern'
      write(ch,16) 'dissest'   , lev,'99','dissipation of KE      '
      write(ch,16) 'keb    '   , lev,'99','kE before SKEB         '
      write(ch,16) 'kea    '   , lev,'99','KE after SKEB          '
      write(ch,16) 'dissdc '   , lev,'99','dissipation from deep convenction'
      write(ch,'(A7)') 'endvars'

      close(ch)
      endif

8     format("(",I4,"(2x,F11.7))")
9     format("(",I4,"(2x,F10.5))")
10    format(8(2x,F11.7))
11    format(8(2x,F10.5))
12    format(A4,1X,I4,1X,A6)
13    format(A4,1X,I4,1X,A10,1X,F10.7)
14    format(A4,1X,I4,1X,A6,1X,A2,A1,A2,A3,A4,1X,A3)
15    format(A5,I1,3X,I3,1X,A2,1X,F7.1,A20)
16    format(A7,3X,I3,1X,A2,1X,A23)

    deallocate(mlat,prsl)
  end subroutine skebctl
!
  subroutine spectrun_inp2d(jcap1,jtr,jtm,mlsort,ns,speci,speco)
!
! use spectral truncation to change resoltuion
!
      implicit none
      integer lev,jtr,jcap1,jtm,ns,ml
      real(kind=RTYPE) speci(jcap1*(jcap1+1)/2,2)
      real(kind=RTYPE) speco(jtr,jtm*ns*2)
      integer i,j,k,jj,jp,jr,j1,j2
      integer mlsort(jcap1,jcap1)
!
      speco(:,:) = 0.0
      if( jcap1.gt.jtr ) then
          do j=1,jcap1
            if( j.le.jtr ) then
              jj=nlist(j)
              jp=(jj-1)/jtm
              jr=mod(jj-1,jtm)+1
              j1=jp*jtm*2+jr
              j2=j1+jtm
            endif
            do i=j,jcap1
              ml=mlsort(j,i)
              if( i.le.jtr ) then
                speco(i,j1) = speci(ml,1)
                speco(i,j2) = speci(ml,2)
              endif
            enddo
          enddo
      else if( jcap1.lt.jtr ) then
          do j=1,jtr
            jj=nlist(j)
            jp=(jj-1)/jtm
            jr=mod(jj-1,jtm)+1
            j1=jp*jtm*2+jr
            j2=j1+jtm
            do i=j,jtr
              if( i.le.jcap1 ) then
                ml=mlsort(j,i)
                speco(i,j1) = speci(ml,1)
                speco(i,j2) = speci(ml,2)
              endif
            enddo
          enddo
      else      ! jcap1=jtr
          do j=1,jtr
            jj=nlist(j)
            jp=(jj-1)/jtm
            jr=mod(jj-1,jtm)+1
            j1=jp*jtm*2+jr
            j2=j1+jtm
            do i=j,jtr
              ml=mlsort(j,i)
              speco(i,j1) = speci(ml,1)
              speco(i,j2) = speci(ml,2)
            enddo
          enddo
      endif

      return
  end subroutine spectrun_inp2d

end module mod_stochastic_physics
