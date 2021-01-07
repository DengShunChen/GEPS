module mod_stochastic_physics
  use mpe
  use rank, only : myrank
  use index
  use param
  use const, only : aki, bki, dosppt, doshum
  use mersenne_twister, only: random_setseed,random_gauss,random_stat
  implicit none
  private 

  type random_pattern
    real, allocatable :: n2d(:,:)
    real, allocatable :: spec(:,:)
    real, allocatable :: poly(:,:)
    real, allocatable :: varspec(:)
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

  integer :: recn=1
  real ::  dt
  logical, public :: ncep_seeds=.false.

  type(random_pattern), public, save, allocatable, dimension(:) :: &
       rpattern_sppt, rpattern_shum

  ! SPPT
  integer :: nsppt
  real, allocatable, save :: sppt3d(:,:,:)
  real :: sppt(5) = -999.             ! amplitude(0.~1.)
  real :: sppt_seed(5) = -999.        ! random seeds
  real :: sppt_decort(5) = -999.      ! time scales(seconds)
  real :: sppt_lscale(5) = -999.      ! length scales(meters)
  real, allocatable, dimension(:) :: vfact_sppt
  real, public :: sppt_sigtop1 = 0.1
  real, public :: sppt_sigtop2 = 0.025
  real, public :: sppt_sigbot1 = 0.975
  real, public :: sppt_sigbot2 = 0.9
  logical, public :: sppt_sfclimit=.true.
  logical, public :: sppt_logit=.false.

  ! SHUM
  integer :: nshum
  real, allocatable, save :: shum3d(:,:,:)
  real :: shum(5) = -999.             ! amplitude(0.~1.)
  real :: shum_seed(5) = -999.         ! random seeds
  real :: shum_decort(5) = -999.      ! time scales(seconds)
  real :: shum_lscale(5) = -999.      ! length scales(meters)
  real, allocatable, dimension(:) :: vfact_shum
  real, public :: shum_sigefold = 0.2

  public random_pattern

  public nsppt, sppt, sppt_seed, sppt_decort, sppt_lscale, sppt3d
  public nshum, shum, shum_seed, shum_decort, shum_lscale, shum3d

  public  init_stochastic_physics, &
           run_stochastic_physics, &
       destroy_stochastic_physics

  public spptout
  public avevar_sppt2d
 
contains

  subroutine init_stochastic_physics(dtau)
    implicit none
    integer :: n, k 
    real :: sl(lev)
    real :: dtau

    ! calculation sigma values
    do k=1,lev
      sl(k)=0.5*(aki(k)/1013.0+bki(k)+aki(k+1)/1013.0+bki(k+1))
    enddo
 
    if (dosppt) then
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
      call get_random_pattern_init(rpattern_sppt,nsppt,dtau)

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
        if (myrank == 0) print *,'mod_stochastic_physics : k,sl,vfact_sppt',k,sl(k),vfact_sppt(k)
      enddo
    endif

    if (doshum) then
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
      shum3d = 0.
      call get_random_pattern_init(rpattern_shum,nshum,dtau)

      allocate(vfact_shum(lev))
      do k=1,lev
         vfact_shum(k) = exp((sl(k)-1.)/shum_sigefold)
         if (sl(k).LT. 2*shum_sigefold) then
            vfact_shum(k)=0.0
         endif
        if (myrank == 0) print *,'mod_stochastic_physics : k,sl,vfact_shum',k,sl(k),vfact_shum(k)
      enddo
    endif

  end subroutine init_stochastic_physics

  subroutine run_stochastic_physics()
    implicit none
    integer :: k 
    real  :: glob(nx,my),temp(nxp,my_max)
    real :: aves,vars,stds

    if (dosppt) then
      call get_random_pattern_run(rpattern_sppt,nsppt)
      call get_stochy_physics(rpattern_sppt,nsppt,vfact_sppt,s  ppt3d)
      if (sppt_logit) sppt3d(:,:,:) = (2./(1.+exp(sppt3d(:,:,:))))-1.
    endif
    if (doshum) then
      call get_random_pattern_run(rpattern_shum,nshum)
      call get_stochy_physics(rpattern_shum,nshum,vfact_shum,shum3d)
    endif

  end subroutine run_stochastic_physics

  subroutine destroy_stochastic_physics()
    implicit none

    if (dosppt) then
      call get_random_pattern_destroy(rpattern_sppt,nsppt)
      deallocate(sppt3d)
      deallocate(rpattern_sppt)
      deallocate(vfact_sppt)
    endif
    if (doshum) then
      call get_random_pattern_destroy(rpattern_shum,nshum)
      deallocate(shum3d)
      deallocate(rpattern_shum)
      deallocate(vfact_shum)
    endif

  end subroutine destroy_stochastic_physics

  subroutine get_random_pattern_init(rpattern,nscale,dt)
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
    integer :: n, k, nscale, ncx, ml, ms, ns
    real :: rerth, pi, var, radsq, correLsq, rkT
    type(random_pattern), intent(inout) :: rpattern(nscale)
    integer :: irand, i
    real :: dt
    real, allocatable :: noise(:,:)
    integer(8) count, count_rate, count_max, count_trunc
    integer(8) :: iscale = 10000000000
    integer :: count4 
 
    rerth = 6.3712e+6      ! radius of earth (m)
    pi = 4.*atan(1.)
    radsq = rerth*rerth

    do n=1,nscale
      ncx = 2.*pi*rerth/rpattern(n)%lenscale
   !  rpattern(n)%jtrun = 2*((1+(ncx-1)/3)/2) 
      rpattern(n)%jtrun = 2*((1+(4*ncx-1)/4)/2)
      rpattern(n)%mlmax = rpattern(n)%jtrun*(rpattern(n)%jtrun+1)/2

      allocate(rpattern(n)%n2d(nxp,my_max))
      allocate(rpattern(n)%spec(rpattern(n)%mlmax,2))
      allocate(rpattern(n)%poly(rpattern(n)%mlmax,my/2))
      allocate(rpattern(n)%varspec(rpattern(n)%mlmax))
      allocate(rpattern(n)%msort(rpattern(n)%mlmax))
      allocate(rpattern(n)%lsort(rpattern(n)%mlmax))
      allocate(rpattern(n)%mlsort(rpattern(n)%jtrun,rpattern(n)%jtrun))
      allocate(noise(rpattern(n)%mlmax,2))

      ! get legendre polynomials 
      call get_legendre_poly(rpattern(n)%mlmax,rpattern(n)%jtrun,rpattern(n)%poly)

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
      call mpe_bcast(count4,1,0,mpe_double) 
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
      deallocate(rpattern(n)%spec)
      deallocate(rpattern(n)%poly)
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
    real :: rpattern2d(nx,my) 

    do n=1,nscale
      if (myrank .eq. 0 ) then
        call gen_random_pattern_2d(rpattern2d,rpattern(n))
      endif
      call mpe_bcast(rpattern2d,nx*my,0,mpe_double)
      do jj=1,jlistnum
        j=jlist1(jj)
        if( lreduce.eq.1 ) then
          call reducepick (rpattern2d(1,j),nxdef(j),nx,1)
        endif
        ii=nxjstart(j)
        nxj=nxdef_2d(j)
        do i=1,nxj
          rpattern(n)%n2d(i,jj)=rpattern2d(ii,j)
          ii=ii+1
        enddo
      enddo
    enddo
  end subroutine get_random_pattern_run

  subroutine get_stochy_physics(rpattern,nscale,vfact,n3d)
!------------------------------------------------------------------------! 
!  purpose: To generate 3D SPPT strucutre
!  output: n3d
!------------------------------------------------------------------------! 
    implicit none
    integer :: n, ii, i, jj, j, k, nxj
    integer, intent(in) :: nscale
    real, intent(in) :: vfact(lev) 
    type(random_pattern), intent(inout) :: rpattern(nscale)
    real, intent(  out) :: n3d(nxp,lev,my_max) 
 
    n3d = 0.
    do n=1,nscale
      do jj=1,jlistnum
        j=jlist1(jj)
        do k=1,lev
          nxj=nxdef_2d(j)
          do i=1,nxj
            n3d(i,k,jj)=n3d(i,k,jj)+rpattern(n)%n2d(i,jj)*vfact(k)
          enddo
        enddo
      enddo
    enddo

  end subroutine get_stochy_physics
  
  subroutine get_noise(rpattern,noise)
    implicit none
    integer :: ml, ns, ms
    type(random_pattern), intent(inout) :: rpattern
    real, intent(out) :: noise(rpattern%mlmax,2)
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

  subroutine gen_random_pattern_2d(sppt2d,rpattern)
    implicit none
    type(random_pattern), intent(inout) :: rpattern
    real, intent(out) :: sppt2d(nx,my)
    integer :: ml, ns, ms
    real, allocatable :: noise(:,:)

    ! get noise
    allocate(noise(rpattern%mlmax,2)) 
    call get_noise(rpattern,noise) 

    !  radom pattern advance with first order AR
    rpattern%spec(:,1) = rpattern%phi*rpattern%spec(:,1) + & 
          sqrt(1.-rpattern%phi**2.)*rpattern%stdev*rpattern%varspec*noise(:,1)
    rpattern%spec(:,2) = rpattern%phi*rpattern%spec(:,2) + & 
          sqrt(1.-rpattern%phi**2.)*rpattern%stdev*rpattern%varspec*noise(:,2)

    ! transform spectral to physical space 
    call transr_sppt(rpattern%jtrun,rpattern%mlmax,nx,my,1,rpattern%poly,rpattern%spec,sppt2d)

    deallocate(noise)

  end subroutine gen_random_pattern_2d

  subroutine get_legendre_poly(mlmax,jtrun,poly)
    implicit none
    integer, intent(in) :: mlmax,jtrun
    real,intent(out)   ::  poly(mlmax,my/2) 

    !   Spheric Harmonic Constants
    integer ::  msort(mlmax),lsort(mlmax),mlsort(mlmax,jtrun)
    real :: dpoly(mlmax,my/2),eps4(mlmax),cim(mlmax)
    real :: cosl(my),onocos(my)
    real :: weight(my),sinl(my)
    real :: cp,capa,rgas,pi,radsq,rad,one,onem,irad
    integer :: rl,rm,rlm,ml
    integer :: j,my2

    !   prepare spheric harmonic coefficients
    cp=1004.24
    capa= 1.0/3.5
    rgas= capa*cp
    pi  = 4.0*atan(1.0)
    rad = 6.371e6
    radsq= rad*rad
! 
!   build pointer arrays for locating zonal and total wavenumber
!   values in the one-dimensional spherical harmonic arrays.
    call sortml (jtrun,mlmax,msort,lsort,mlsort)
! 
    do ml = 1, mlmax
      rl = lsort(ml)
      rm = msort(ml)-1
      rlm= rl-1.0
      if (msort(ml).eq.1)  rm= 0.0
      if (lsort(ml).eq.1)  rlm= 0.0
      eps4(ml)= rl*rlm/radsq
      cim(ml)= rm
    enddo
 
!   gaussian quadrature weights and latitudes
    one = 1.0
    onem= -one
    call gausl3 (my,onem,one,weight,sinl)
! 
    my2= my/2
    do j = 1, my2
      sinl(my+1-j)  = -sinl(j)
      weight(my+1-j)= weight(j)
      onocos(j)     = 1.0/(1.0-sinl(j)*sinl(j))
      onocos(my+1-j)= onocos(j)
      cosl(j)       = 1.0/sqrt(onocos(j))
      cosl(my+1-j)  = cosl(j)
    enddo
! 
!   define associated legendre polynomials and their derivatives 
    call lgndr_sppt(my2,jtrun,mlmax,mlsort,sinl,poly,dpoly)

  end subroutine get_legendre_poly  

! FUNCTION fun_gasdev(idum)result(gasdev)
!     INTEGER :: idum
!     REAL :: gasdev
!     INTEGER,save :: iset=0
!     REAL :: fac,rsq,v1,v2,ran1
!     real, save :: gset

!     if (iset.eq.0) then
!  11   v1=2.*fun_ran1(idum)-1.
!       v2=2.*fun_ran1(idum)-1.
!       rsq=v1**2+v2**2
!       if(rsq.ge.1..or.rsq.eq.0.)goto 11
!       fac=sqrt(-2.*log(rsq)/rsq)
!       gset=v1*fac
!       gasdev=v2*fac
!       iset=1
!     else
!       gasdev=gset
!       iset=0
!     endif
! END FUNCTION fun_gasdev

! FUNCTION fun_ran1(idum)result(ran1)
!     INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
!     REAL ran1,AM,EPS,RNMX
!     PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
!     NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
!     INTEGER j,k,iv(NTAB),iy
!     SAVE iv,iy
!     DATA iv /NTAB*0/, iy /0/
!     if (idum.le.0.or.iy.eq.0) then
!       idum=max(-idum,1)
!       do j=NTAB+8,1,-1
!         k=idum/IQ
!         idum=IA*(idum-k*IQ)-IR*k
!         if (idum.lt.0) idum=idum+IM
!         if (j.le.NTAB) iv(j)=idum
!       enddo
!       iy=iv(1)
!     endif
!     k=idum/IQ
!     idum=IA*(idum-k*IQ)-IR*k
!     if (idum.lt.0) idum=idum+IM
!     j=1+iy/NDIV
!     iy=iv(j)
!     iv(j)=idum
!     ran1=min(AM*iy,RNMX)
!     return
! END FUNCTION fun_ran1

! FUNCTION fun_ran2(idum)result(ran2)
!     INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
!     REAL ran2,AM,EPS,RNMX
!     PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
!     IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,    &
!     NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
!     INTEGER idum2,j,k,iv(NTAB),iy
!     SAVE iv,iy,idum2
!     DATA idum2/123456789/, iv/NTAB*0/, iy/0/
!     if (idum.le.0) then
!       idum=max(-idum,1)
!       idum2=idum
!       do j=NTAB+8,1,-1
!         k=idum/IQ1
!         idum=IA1*(idum-k*IQ1)-k*IR1
!         if (idum.lt.0) idum=idum+IM1
!         if (j.le.NTAB) iv(j)=idum
!       enddo
!       iy=iv(1)
!     endif
!     k=idum/IQ1
!     idum=IA1*(idum-k*IQ1)-k*IR1
!     if (idum.lt.0) idum=idum+IM1
!     k=idum2/IQ2
!     idum2=IA2*(idum2-k*IQ2)-k*IR2
!     if (idum2.lt.0) idum2=idum2+IM2
!     j=1+iy/NDIV
!     iy=iv(j)-idum2
!     iv(j)=idum
!     if(iy.lt.1)iy=iy+IMM1
!     ran2=min(AM*iy,RNMX)
! END FUNCTION fun_ran2


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
    REAL :: ave,var,data2d(n,m),data(n*m)
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

  subroutine removegt2std(scaleval,nx,my,aves,stds)
      integer i,j,nx,my
      real scaleval(nx,my),aves,stds
      real vcheck2p,vcheck2n
!c
      vcheck2p=2.0*stds
      vcheck2n=-2.0*stds
      do j=1,my
      do i=1,nx
        if ((scaleval(i,j)-aves).gt.vcheck2p) scaleval(i,j)=vcheck2p
        if ((scaleval(i,j)-aves).lt.vcheck2n) scaleval(i,j)=vcheck2n 
      enddo
      enddo
  end subroutine removegt2std
  subroutine removegt3std(scaleval,nx,my,aves,stds)
      integer i,j,nx,my
      real scaleval(nx,my),aves,stds
      real vcheck3p,vcheck3n
!c
      vcheck3p=3.0*stds
      vcheck3n=-3.0*stds

      do j=1,my
      do i=1,nx
        if ((scaleval(i,j)-aves).gt.vcheck3p) scaleval(i,j)=vcheck3p
        if ((scaleval(i,j)-aves).lt.vcheck3n) scaleval(i,j)=vcheck3n 
      enddo
      enddo
!    
  end subroutine removegt3std
  subroutine lgndr_sppt(my2,jtrun,mlmax,mlsort,sinl,poly,dpoly)
!c
!c  generate legendre polynomials and their derivatives on the
!c  gaussian latitudes
!c
!c ***input***
!c
!c  my2:  number of gaussian latitudes from south pole and equator
!c  jtrun:  zonal wavenumber truncation limit
!c  mlmax: total number of triangular truncation spherical harmonics
!c  mlsort: pointer array of 1-d indexs at functions of zonal and
!c          total wavenumbers
!c  sinl: sin of gaussian latitudes
!c
!c  ***output***
!c
!c  poly: associated legendre coefficients
!c  dpoly: d(poly)/d(sinl)
!c
!c ******************************************************************
!c
!c ref= belousov, s. l., 1962= tables of normalized associated
!c        legendre polynomials. pergamon press, new york
!c
      integer :: j, n, np, kp, k, mp, m, nps, l, ml, m1, mk
      integer :: my2, jtrun, mlmax, jtrunp
      real :: poly(mlmax,my2),dpoly(mlmax,my2),sinl(my2)
      integer :: mlsort(jtrun,jtrun)
!c
!c      parameter (jtrunx= 100)
      real :: pnm(jtrun+1,jtrun+1),dpnm(jtrun+1,jtrun+1)
      real :: xx, sn, sn2i, rt2, c1, fn, fn2, fn2s, c3, s1, s2, c4, c5, c6, cf
      real :: a, b, fk, fm, fm1, fm2, fm3, c7, c8, c, d, e, fms, fnp, fnp2
      real :: theta, ang
!c
!c sinl is sin(latitude) = cos(colatitude)
!c pnm(np,mp) is legendre polynomial p(n,m) with np=n+1, mp=m+1
!c pnm(mp,np+1) is x derivative of p(n,m) with np=n+1, mp=m+1
!c
      jtrunp= jtrun+1
      do 1001 j=1,my2
      xx= sinl(j)
      sn= sqrt(1.0-xx*xx)
	sn2i = 1.0/(1.0 - xx*xx)
      rt2= sqrt(2.0)
	c1 = rt2
!c
	pnm(1,1) = 1.0/rt2
      theta=-atan(xx/sqrt(1.0-xx*xx))+2.0*atan(1.0)
!c
      do 20 n=1,jtrun
	np = n + 1
      fn=n
	fn2 = fn + fn
	fn2s = fn2*fn2
!c eq 22
      c1= c1*sqrt(1.0-1.0/fn2s)
      c3= c1/sqrt(fn*(fn+1.0))
	ang = fn*theta
	s1 = 0.0
	s2 = 0.0
	c4 = 1.0
	c5 = fn
	a = -1.0
	b = 0.0
!c
      do 27 kp=1,np,2
	k = kp - 1
      s2= s2+c5*sin(ang)*c4
      if (k.eq.n) c4 = 0.5*c4
      s1= s1+c4*cos(ang)
	a = a + 2.0
	b = b + 1.0
      fk=k
	ang = theta*(fn - fk - 2.0)
	c4 = (a*(fn - b + 1.0)/(b*(fn2 - a)))*c4
	c5 = c5 - 2.0
   27 continue
!c eq 19
	pnm(np,1) = s1*c1
!c eq 21
	pnm(np,2) = s2*c3
   20 continue
!c
      do 4 mp=3,jtrunp
	m = mp - 1
      fm= m
	fm1 = fm - 1.0
	fm2 = fm - 2.0
	fm3 = fm - 3.0
      c6= sqrt(1.0+1.0/(fm+fm))
!c eq 23
	pnm(mp,mp) = c6*sn*pnm(m,m)
      if (mp - jtrunp) 3,4,4
    3 continue
	nps = mp + 1
!c
      do 41 np=nps,jtrunp
	n = np - 1
      fn= n
	fn2 = fn + fn
	c7 = (fn2 + 1.0)/(fn2 - 1.0)
	c8 = (fm1 + fn)/((fm + fn)*(fm2 + fn))
      c= sqrt((fn2+1.0)*c8*(fm3+fn)/(fn2-3.0))
      d= -sqrt(c7*c8*(fn-fm1))
      e= sqrt(c7*(fn-fm)/(fn+fm))
!c eq 17
	pnm(np,mp) = c*pnm(np-2,mp-2) &
                 + xx*(d*pnm(np-1,mp-2) + e*pnm(np - 1,mp))
   41 continue
    4 continue
!c
      do 50 mp=1,jtrun
      fm= mp-1.0
	fms = fm*fm
      do 50 np=mp,jtrun
      fnp= np
	fnp2 = fnp + fnp
	cf = (fnp*fnp - fms)*(fnp2 - 1.0)/(fnp2 + 1.0)
      cf= sqrt(cf)
!c der
      dpnm(np,mp)   = -sn2i*(cf*pnm(np+1,mp) - fnp*xx*pnm(np,mp))
   50 continue
!c
      do 71 m=1,jtrun
      do 71 l=m,jtrun
      ml= mlsort(m,l)
      poly(ml,j)= pnm(l,m)
      dpoly(ml,j)=dpnm(l,m)
   71 continue
      dpoly(1,j)= 0.0
 1001 continue
      
  end subroutine lgndr_sppt

  subroutine transr_sppt (jtrun,mlmax,nx,my,ll,poly,s,r)
!   subroutine to transform a spectral coefficient field to
!   grid point form
! 
!  *** const ***
! 
!   jtrun: zonal wavenumber resolution limit
!   mlmax: number of spectral coefficients (horizontal field)
!   nx: e-w dimension no.
!   my: n-s dimension no.
!   ll: number of levels to transform
!   poly: legendre polynomials
! 
!  *** input variable ***
! 
!   s: spectral coefficient array to transform
! 
!  *** output variable ***
! 
!   r: 3-d output grid point fields
! 
!   **************************************

      use fftcom  
      real :: poly(mlmax,my/2),s(mlmax,2,ll),r(nx,ll,my)
      real :: cc(nx+2,my),work(nx*my,2)
      integer :: mlmax,my,nx,ll,jtrun
      integer :: mlx,i,k,m,j,jj,ml,mm,mp
      integer :: m1, l, mk, n

      mlx= (jtrun/2)*((jtrun+1)/2)
      do k=1,ll
        do m=1,(nx+2)*my/2
          cc(m,1)= 0.0
          cc(m,my/2+1)= 0.0
        enddo
        do j=1,my/2
          jj = my+1-j
          ml = 2*mlx
          do m=2,jtrun,2
            ml = ml+1
            mm = 2*m-1
            mp = mm+1
            cc(mm,j ) = poly(ml,j)*s(ml,1,k)
            cc(mp,j ) = poly(ml,j)*s(ml,2,k)
            cc(mm,jj) = cc(mm,j)
            cc(mp,jj) = cc(mp,j)
          enddo

          m1= 0
          do l=jtrun-1,1,-2
            do m=1,l
              mm = 2*m-1
              mp = mm+1
              ml = m+m1
              mk = ml+mlx
              cc(mm,j ) = cc(mm,j ) + poly(ml,j)*s(ml,1,k) + poly(mk,j)*s(mk,1,k)
              cc(mm,jj) = cc(mm,jj) + poly(ml,j)*s(ml,1,k) - poly(mk,j)*s(mk,1,k)
              cc(mp,j ) = cc(mp,j ) + poly(ml,j)*s(ml,2,k) + poly(mk,j)*s(mk,2,k)
              cc(mp,jj) = cc(mp,jj) + poly(ml,j)*s(ml,2,k) - poly(mk,j)*s(mk,2,k)
            enddo
            m1= m1+l
          enddo
        enddo
!       call fft991(cc,work,trigs,ifax,1,nx+3,nx,my,1)
        call rfftmlt(cc,work,trigs,ifax,1,nx+2,nx,my,1)
        do j=1,my
          do i=1,nx
            r(i,k,j)= cc(i,j)
          enddo
        enddo
      enddo
      
  end subroutine transr_sppt

  subroutine spptout(tau)
    implicit none
    integer      :: i, j, k, jj, nxj, ihead, n
    integer      :: nxmy4
    real         :: tau
    real         :: glob(nx,my),temp(nxp,my_max)
    real(kind=4) :: glob4(nx,my)

    ihead=15
    nxmy4=nx*my*4
    if ( myrank .eq. 0 ) then
      open(ihead,file='sppt.dat',access='direct',form='unformatted',recl=nxmy4,status='unknown')
    endif

    do n=1,nsppt
      call unify_reduceintp(nx,my,my_max,rpattern_sppt(n)%n2d,glob)   
      if ( myrank .eq. 0 ) then
        glob4=glob
        write(ihead,rec=recn) glob4
        recn=recn+1
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
        write(ihead,rec=recn) glob4
        recn=recn+1
      endif
    enddo

    if ( myrank.eq.0 ) close(ihead)

    !creating ctl file
    call spptctl(tau)
 
  end subroutine spptout
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
      prsl(kk)=0.5*prsl(kk)/1000.
    enddo

      write(forydef,8) my
      write(forzdef,9) lev

      write(dtg,'(I12)') idtg
      read(dtg,'(A4,I2,A2,A2,A2)') yy,mn,dd,hh,mm

      if ( myrank .eq. 0 ) then
      OPEN(UNIT=ch, FILE='sppt.ctl', STATUS='UNKNOWN'             &
         , ACCESS='SEQUENTIAL')

      write(ch,'(A23)') 'dset ^sppt.dat'
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
       write(10,10) mlat(js:je)
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
       write(10,11) prsl(js:je)
      enddo
      if ( .not. jrem ) write(10,forzdef) prsl(je+1:my)
      write(ch,'(A4,1X,I2)') 'vars',nsppt+1
      do n=1,nsppt
        write(ch,15) 'scale',n ,lev1,'99',rpattern_sppt(n)%lenscale/1000.,'km 2D Random Pattern'
      enddo
      write(ch,16) 'sppt  '   , lev,'99','3D Random Pattern'
      write(ch,'(A7)') 'endvars'

      close(ch)
      endif

8     format("(",I4,"(2x,F11.7))")
9     format("(",I4,"(2x,F7.5))")
10    format(8(2x,F11.7))
11    format(8(2x,F7.5))
12    format(A4,1X,I4,1X,A6)
13    format(A4,1X,I4,1X,A10,1X,F10.7)
14    format(A4,1X,I4,1X,A6,1X,A2,A1,A2,A3,A4,1X,A3)
15    format(A5,I1,3X,I3,1X,A2,1X,F7.1,A20)
16    format(A6,3X,I3,1X,A2,1X,A20)

    deallocate(mlat,prsl)
  end subroutine spptctl

end module mod_stochastic_physics
