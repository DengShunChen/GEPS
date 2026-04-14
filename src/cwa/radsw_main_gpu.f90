!!!!!  ==============================================================  !!!!!
!!!!!              sw-rrtm3 radiation package description              !!!!!
!!!!!  ==============================================================  !!!!!
!                                                                          !
!   this package includes ncep's modifications of the rrtm-sw radiation    !
!   code from aer inc.                                                     !
!                                                                          !
!   the sw-rrtm3 package includes these parts:                             !
!                                                                          !
!      'radsw_rrtm3_param.f'                                               !
!      'radsw_rrtm3_datatb.f'                                              !
!      'radsw_rrtm3_main.f'                                                !
!                                                                          !
!   the 'radsw_rrtm3_param.f' contains:                                    !
!                                                                          !
!      'module_radsw_parameters'  -- band parameters set up                !
!                                                                          !
!   the 'radsw_rrtm3_datatb.f' contains:                                   !
!                                                                          !
!      'module_radsw_ref'         -- reference temperature and pressure    !
!      'module_radsw_cldprtb'     -- cloud property coefficients table     !
!      'module_radsw_sflux'       -- spectral distribution of solar flux   !
!      'module_radsw_kgbnn'       -- absorption coeffients for 14          !
!                                    bands, where nn = 16-29               !
!                                                                          !
!   the 'radsw_rrtm3_main.f' contains:                                     !
!                                                                          !
!      'module_radsw_main'        -- main sw radiation transfer            !
!                                                                          !
!   in the main module 'module_radsw_main' there are only two              !
!   externally callable subroutines:                                       !
!                                                                          !
!      'swrad'      -- main sw radiation routine                           !
!         inputs:                                                          !
!           (plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                         !
!            clouds,icseed,aerosols,sfcalb,                                !
!            cosz,solcon,nday,idxday,                                      !
!            npts, nlay, nlp1, lprnt,                                      !
!         outputs:                                                         !
!            hswc,topflx,sfcflx,                                           !
!!        optional outputs:                                                !
!            hsw0,hswb,flxprf,fdncmp)                                      !
!           )                                                              !
!                                                                          !
!      'rswinit'    -- initialization routine                              !
!         inputs:                                                          !
!           ( me )                                                         !
!         outputs:                                                         !
!           (none)                                                         !
!                                                                          !
!   all the sw radiation subprograms become contained subprograms          !
!   in module 'module_radsw_main' and many of them are not directly        !
!   accessable from places outside the module.                             !
!                                                                          !
!    derived data type constructs used:                                    !
!                                                                          !
!     1. radiation flux at toa: (from module 'module_radsw_parameters')    !
!          topfsw_type   -  derived data type for toa rad fluxes           !
!            upfxc              total sky upward flux at toa               !
!            dnfxc              total sky downward flux at toa             !
!            upfx0              clear sky upward flux at toa               !
!                                                                          !
!     2. radiation flux at sfc: (from module 'module_radsw_parameters')    !
!          sfcfsw_type   -  derived data type for sfc rad fluxes           !
!            upfxc              total sky upward flux at sfc               !
!            dnfxc              total sky downward flux at sfc             !
!            upfx0              clear sky upward flux at sfc               !
!            dnfx0              clear sky downward flux at sfc             !
!                                                                          !
!     3. radiation flux profiles(from module 'module_radsw_parameters')    !
!          profsw_type    -  derived data type for rad vertical prof       !
!            upfxc              level upward flux for total sky            !
!            dnfxc              level downward flux for total sky          !
!            upfx0              level upward flux for clear sky            !
!            dnfx0              level downward flux for clear sky          !
!                                                                          !
!     4. surface component fluxes(from module 'module_radsw_parameters'    !
!          cmpfsw_type    -  derived data type for component sfc flux      !
!            uvbfc              total sky downward uv-b flux at sfc        !
!            uvbf0              clear sky downward uv-b flux at sfc        !
!            nirbm              surface downward nir direct beam flux      !
!            nirdf              surface downward nir diffused flux         !
!            visbm              surface downward uv+vis direct beam flx    !
!            visdf              surface downward uv+vis diffused flux      !
!                                                                          !
!   external modules referenced:                                           !
!                                                                          !
!       'module physpara'                                                  !
!       'module physcons'                                                  !
!       'mersenne_twister'                                                 !
!                                                                          !
!   compilation sequence is:                                               !
!                                                                          !
!      'radsw_rrtm3_param.f'                                               !
!      'radsw_rrtm3_datatb.f'                                              !
!      'radsw_rrtm3_main.f'                                                !
!                                                                          !
!   and all should be put in front of routines that use sw modules         !
!                                                                          !
!==========================================================================!
!                                                                          !
!   the original program declarations:                                     !
!                                                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                          !
!  copyright 2002-2007, atmospheric & environmental research, inc. (aer).  !
!  this software may be used, copied, or redistributed as long as it is    !
!  not sold and this copyright notice is reproduced on each copy made.     !
!  this model is provided as is without any express or implied warranties. !
!                       (http://www.rtweb.aer.com/)                        !
!                                                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                          !
! ************************************************************************ !
!                                                                          !
!                              rrtmg_sw                                    !
!                                                                          !
!                                                                          !
!                   a rapid radiative transfer model                       !
!                    for the solar spectral region                         !
!            atmospheric and environmental research, inc.                  !
!                        131 hartwell avenue                               !
!                        lexington, ma 02421                               !
!                                                                          !
!                           eli j. mlawer                                  !
!                        jennifer s. delamere                              !
!                         michael j. iacono                                !
!                         shepard a. clough                                !
!                                                                          !
!                                                                          !
!                       email:  miacono@aer.com                            !
!                       email:  emlawer@aer.com                            !
!                       email:  jdelamer@aer.com                           !
!                                                                          !
!        the authors wish to acknowledge the contributions of the          !
!        following people:  steven j. taubman, patrick d. brown,           !
!        ronald e. farren, luke chen, robert bergstrom.                    !
!                                                                          !
! ************************************************************************ !
!                                                                          !
!    references:                                                           !
!    (rrtm_sw/rrtmg_sw):                                                   !
!      clough, s.a., m.w. shephard, e.j. mlawer, j.s. delamere,            !
!      m.j. iacono, k. cady-pereira, s. boukabara, and p.d. brown:         !
!      atmospheric radiative transfer modeling: a summary of the aer       !
!      codes, j. quant. spectrosc. radiat. transfer, 91, 233-244, 2005.    !
!                                                                          !
!    (mcica):                                                              !
!      pincus, r., h. w. barker, and j.-j. morcrette: a fast, flexible,    !
!      approximation technique for computing radiative transfer in         !
!      inhomogeneous cloud fields, j. geophys. res., 108(d13), 4376,       !
!      doi:10.1029/2002jd003322, 2003.                                     !
!                                                                          !
! ************************************************************************ !
!                                                                          !
!    aer's revision history:                                               !
!     this version of rrtmg_sw has been modified from rrtm_sw to use a     !
!     reduced set of g-point intervals and a two-stream model for          !
!     application to gcms.                                                 !
!                                                                          !
! --  original version (derived from rrtm_sw)                              !
!        2002: aer. inc.                                                   !
! --  conversion to f90 formatting; addition of 2-stream radiative transfer!
!        feb 2003: j.-j. morcrette, ecmwf                                  !
! --  additional modifications for gcm application                         !
!        aug 2003: m. j. iacono, aer inc.                                  !
! --  total number of g-points reduced from 224 to 112.  original          !
!     set of 224 can be restored by exchanging code in module parrrsw.f90  !
!     and in file rrtmg_sw_init.f90.                                       !
!        apr 2004: m. j. iacono, aer, inc.                                 !
! --  modifications to include output for direct and diffuse               !
!     downward fluxes.  there are output as "true" fluxes without          !
!     any delta scaling applied.  code can be commented to exclude         !
!     this calculation in source file rrtmg_sw_spcvrt.f90.                 !
!        jan 2005: e. j. mlawer, m. j. iacono, aer, inc.                   !
! --  revised to add mcica capability.                                     !
!        nov 2005: m. j. iacono, aer, inc.                                 !
! --  reformatted for consistency with rrtmg_lw.                           !
!        feb 2007: m. j. iacono, aer, inc.                                 !
! --  modifications to formatting to use assumed-shape arrays.             !
!        aug 2007: m. j. iacono, aer, inc.                                 !
!                                                                          !
! ************************************************************************ !
!                                                                          !
!   ncep modifications history log:                                        !
!                                                                          !
!       sep 2003,  yu-tai hou        -- received aer's rrtm-sw gcm version !
!                    code (v224)                                           !
!       nov 2003,  yu-tai hou        -- corrected errors in direct/diffuse !
!                    surface alabedo components.                           !
!       jan 2004,  yu-tai hou        -- modified code into standard modular!
!                    f9x code for ncep models. the original three cloud    !
!                    control flags are simplified into two: iflagliq and   !
!                    iflagice. combined the org subr sw_224 and setcoef    !
!                    into radsw (the main program); put all kgb##together  !
!                    and reformat into a separated data module; combine    !
!                    reftra and vrtqdr as swflux; optimized taumol and all !
!                    taubgs to form a contained subroutines.               !
!       jun 2004,  yu-tai hou        -- modified code based on aer's faster!
!                    version rrtmg_sw (v2.0) with 112 g-points.            !
!       mar 2005,  yu-tai hou        -- modified to aer v2.3, correct cloud!
!                    scaling error, total sky properties are delta scaled  !
!                    after combining clear and cloudy parts. the testing   !
!                    criterion of ssa is saved before scaling. added cloud !
!                    layer rain and snow contributions. all cloud water    !
!                    partical contents are treated the same way as other   !
!                    atmos particles.                                      !
!       apr 2005,  yu-tai hou        -- modified on module structures (this!
!                    version of code was given back to aer in jun 2006)    !
!       nov 2006,  yu-tai hou        -- modified code to include the       !
!                    generallized aerosol optical property scheme for gcms.!
!       apr 2007,  yu-tai hou        -- added spectral band heating as an  !
!                    optional output to support the 500km model's upper    !
!                    stratospheric radiation calculations. restructure     !
!                    optional outputs for easy access by different models. !
!       oct 2008,  yu-tai hou        -- modified to include new features   !
!                    from aer's newer release v3.5-v3.61, including mcica  !
!                    sub-grid cloud option and true direct/diffuse fluxes  !
!                    without delta scaling. added rain/snow opt properties !
!                    support to cloudy sky calculations. simplified and    !
!                    unified sw and lw sub-column cloud subroutines into   !
!                    one module by using optional parameters.              !
!       mar 2009,  yu-tai hou        -- replaced the original random number!
!                    generator coming with the original code with ncep w3  !
!                    library to simplify the program and moved sub-column  !
!                    cloud subroutines inside the main module. added       !
!                    option of user provided permutation seeds that could  !
!                    be randomly generated from forecast time stamp.       !
!       mar 2009,  yu-tai hou        -- replaced random number generator   !
!                    programs coming from the original code with the ncep  !
!                    w3 library to simplify the program and moved sub-col  !
!                    cloud subroutines inside the main module. added       !
!                    option of user provided permutation seeds that could  !
!                    be randomly generated from forecast time stamp.       !
!       nov 2009,  yu-tai hou        -- updated to aer v3.7-v3.8 version.  !
!                    notice the input cloud ice/liquid are assumed as      !
!                    in-cloud quantities, not grid average quantities.     !
!       aug 2010,  yu-tai hou        -- uptimized code to improve efficiency
!                    splited subroutine spcvrt into two subs, spcvrc and   !
!                    spcvrm, to handling non-mcica and mcica type of calls.!
!       apr 2012,  b. ferrier and y. hou -- added conversion factor to fu's!
!                    cloud-snow optical property scheme.                   !
!       jul 2012,  s. moorthi and y. hou  -- eliminated the pointer array  !
!                     in subr 'spcvrt' for multi-threading issue running   !
!                     under intel's fortran compiler.                      !
!       nov 2012,  yu-tai hou        -- modified control parameters thru   !
!                     module 'physpara'.                                   !
!       jun 2013,  yu-tai hou        -- moving band 9 surface treatment    !
!                     back as in the rrtm2 version, spliting surface flux  !
!                     into two spectral regions (vis & nir), instead of    !
!                     designated it in nir region only.                    !
!                                                                          !
!!!!!  ==============================================================  !!!!!
!!!!!                         end descriptions                         !!!!!
!!!!!  ==============================================================  !!!!!


!========================================!
      module module_radsw_main_gpu           !
!........................................!
!
      use physpara,         only : iswrate, iswrgas, iswcliq, iswcice,  &
     &                             isubcsw, icldflg, iovrsw,  ivflip,   &
     &                             iswmode, kind_phys
      use physcons,         only : con_g, con_cp, con_avgd, con_amd,    &
     &                             con_amw, con_amo3
     use module_radiation_aerosols_gpu,only : setaer_sw_gpu

      use module_radsw_parameters
      use mersenne_twister, only : random_setseed, random_number,       &
     &                             random_stat
      use module_radsw_ref, only : preflog, tref
      use module_radsw_sflux
      use param,             only : my
      use index,             only : jlistnum, nxptot, nxjp_acc
      use rank,              only : myrank
      !use nvtx
!
      implicit none
!
      private
!
!  ---  version tag and last revision date
      character(40), parameter ::                                       &
     &   vtagsw='ncep sw v5.1  nov 2012 -rrtmg-sw v3.8   '
!    &   vtagsw='ncep sw v5.0  aug 2012 -rrtmg-sw v3.8   '
!    &   vtagsw='rrtmg-sw v3.8   nov 2009'
!    &   vtagsw='rrtmg-sw v3.7   nov 2009'
!    &   vtagsw='rrtmg-sw v3.61  oct 2008'
!    &   vtagsw='rrtmg-sw v3.5   oct 2008'
!    &   vtagsw='rrtm-sw 112v2.3 apr 2007'
!    &   vtagsw='rrtm-sw 112v2.3 mar 2005'
!    &   vtagsw='rrtm-sw 112v2.0 jul 2004'

!  ---  constant values
      real (kind=kind_phys), parameter :: eps     = 1.0e-6
      real (kind=kind_phys), parameter :: oneminus= 1.0 - eps
      real (kind=kind_phys), parameter :: bpade   = 1.0/0.278  ! pade approx constant
      real (kind=kind_phys), parameter :: stpfac  = 296.0/1013.0
      real (kind=kind_phys), parameter :: ftiny   = 1.0e-12
      real (kind=kind_phys), parameter :: s0      = 1368.22    ! internal solar const
                                                               ! adj through input
      real (kind=kind_phys), parameter :: f_zero  = 0.0
      real (kind=kind_phys), parameter :: f_one   = 1.0

!  ---  atomic weights for conversion from mass to volume mixing ratios
      real (kind=kind_phys), parameter :: amdw    = con_amd/con_amw
      real (kind=kind_phys), parameter :: amdo3   = con_amd/con_amo3

!  ---  band indices
      integer, dimension(nblow:nbhgh) :: nspa, nspb, idxebc, idxsfc

      data nspa(:) /  9, 9, 9, 9, 1, 9, 9, 1, 9, 1, 0, 1, 9, 1 /
      data nspb(:) /  1, 5, 1, 1, 1, 5, 1, 0, 1, 0, 0, 1, 5, 1 /

!     data idxsfc(:) / 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 1 /  ! band index for sfc flux
      data idxsfc(:) / 1, 1, 1, 1, 1, 1, 1, 1, 0, 2, 2, 2, 2, 1 /  ! band index for sfc flux
      data idxebc(:) / 5, 5, 4, 4, 3, 3, 2, 2, 1, 1, 1, 1, 1, 5 /  ! band index for cld prop

!  ---  band wavenumber intervals
!     real (kind=kind_phys), dimension(nblow:nbhgh):: wavenum1,wavenum2
!     data wavenum1(:)  /                                               &
!    &         2600.0, 3250.0, 4000.0, 4650.0, 5150.0, 6150.0, 7700.0,  &
!    &         8050.0,12850.0,16000.0,22650.0,29000.0,38000.0,  820.0 /
!     data wavenum2(:)  /                                               &
!              3250.0, 4000.0, 4650.0, 5150.0, 6150.0, 7700.0, 8050.0,  &
!    &        12850.0,16000.0,22650.0,29000.0,38000.0,50000.0, 2600.0 /
!     real (kind=kind_phys), dimension(nblow:nbhgh) :: delwave
!     data delwave(:)   /                                               &
!    &          650.0,  750.0,  650.0,  500.0, 1000.0, 1550.0,  350.0,  &
!    &         4800.0, 3150.0, 6650.0, 6350.0, 9000.0,12000.0, 1780.0 /

      integer, parameter :: nuvb = 27            !uv-b band index

!! ---  logical flags for optional output fields

      logical :: lhswb  = .false.
      logical :: lhsw0  = .false.
      logical :: lflxprf = .false.
      logical :: lflxprf_upfxc = .false.
      logical :: lflxprf_dnfxc = .false.
      logical :: lflxprf_upfx0 = .false.
      logical :: lflxprf_dnfx0 = .false.
      logical :: lfdncmp = .false.
      logical :: lfdncmp_uvbfc = .false.
      logical :: lfdncmp_uvbf0 = .false.
      logical :: lfdncmp_nirbm = .false.
      logical :: lfdncmp_nirdf = .false.
      logical :: lfdncmp_visbm = .false.
      logical :: lfdncmp_visdf = .false.
      

!  ---  those data will be set up only once by "rswinit"

      real (kind=kind_phys) :: exp_tbl(0:ntbmx)

!  ...  heatfac is the factor for heating rates
!       (in k/day, or k/sec set by subroutine 'rswinit')

      real (kind=kind_phys) :: heatfac

!  ---  the following variables are used for sub-column cloud scheme

      integer, parameter :: ipsdsw0 = 1          ! initial permutation seed

!  ---  for ngptsw packing recipe

!  ---  public accessable subprograms
      integer, parameter :: npacks = 4
      integer :: max_packs_nbands
      integer, dimension(:), allocatable :: pack_size, packs_nbands_list
      integer, dimension(:,:), allocatable :: pack_idx, recipe, reverse_recipe
      integer :: small_ngpts

      public swrad_gpu, rswinit_gpu, cldprop, setcoef, taumol, spcvrtc, &
         copyin_radsw_main_gpu


! =================
      contains
! =================

      subroutine ngptsw_packing_recipe_init(npacks, small_ngptsw)
      implicit none
      
      integer, intent(in) :: npacks
      integer, intent(out) :: small_ngptsw

      integer, dimension(4:11) :: max_pack_size
      integer :: np, nib
      data max_pack_size(4:11) /28, 24, 20, 18, 16, 0, 0, 12/

      small_ngptsw = max_pack_size(npacks)
      allocate(pack_size(npacks))
      allocate(packs_nbands_list(npacks))
      allocate(pack_idx(small_ngptsw, npacks))
      allocate(recipe(nbdsw, npacks))
      allocate(reverse_recipe(small_ngptsw, npacks))
      recipe = 0
      pack_size = 0
      pack_idx = 0

      select case (npacks)
      case (4)
         recipe(1:3, 1) = (/17, 18, 19/)
         recipe(1:3, 2) = (/24, 27, 29/)
         recipe(1:4, 3) = (/16, 20, 21, 22/)
         recipe(1:4, 4) = (/23, 25, 26, 28/)
         max_packs_nbands = 4
         packs_nbands_list = (/3, 3, 4, 4/)

      case (5)
         recipe(1:3, 1) = (/17, 18, 22/)
         recipe(1:2, 2) = (/20, 29/)
         recipe(1:3, 3) = (/16, 19, 21/)
         recipe(1:3, 4) = (/23, 24, 25/)
         recipe(1:3, 5) = (/26, 27, 28/)
         max_packs_nbands = 3
         packs_nbands_list = (/3, 2, 3, 3, 3/)

      case (6)
         recipe(1:2, 1) = (/17, 18/)
         recipe(1:2, 2) = (/19, 29/)
         recipe(1:2, 3) = (/20, 21/)
         recipe(1:3, 4) = (/22, 23, 24/)
         recipe(1:3, 5) = (/16, 25, 27/)
         recipe(1:2, 6) = (/26, 28/)
         max_packs_nbands = 3
         packs_nbands_list = (/2, 2, 2, 3, 3, 2/)
         
      case (7)
         recipe(1:2, 1) = (/16, 17/)
         recipe(1:2, 2) = (/25, 29/)
         recipe(1:2, 3) = (/18, 20/)
         recipe(1:2, 4) = (/19, 21/)
         recipe(1:3, 5) = (/22, 23, 26/)
         recipe(1:2, 6) = (/24, 27/)
         recipe(1:1, 7) = (/28/)
         max_packs_nbands = 3
         packs_nbands_list = (/2, 2, 2, 2, 3, 2, 1/)
         
      case (8)
         recipe(1:1, 1) = (/17/)
         recipe(1:1, 2) = (/29/)
         recipe(1:2, 3) = (/16, 20/)
         recipe(1:2, 4) = (/21, 25/)
         recipe(1:2, 5) = (/23, 26/)
         recipe(1:2, 6) = (/18, 19/)
         recipe(1:2, 7) = (/24, 28/)
         recipe(1:2, 8) = (/22, 27/)
         max_packs_nbands = 2
         packs_nbands_list = (/1, 1, 2, 2, 2, 2, 2, 2/)
      
      case (11)
         recipe(1:1, 1) = (/17/)
         recipe(1:1, 2) = (/29/)
         recipe(1:1, 3) = (/20/)
         recipe(1:1, 4) = (/21/)
         recipe(1:1, 5) = (/23/)
         recipe(1:2, 6) = (/18, 22/)
         recipe(1:1, 7) = (/19/)
         recipe(1:1, 8) = (/24/)
         recipe(1:1, 9) = (/27/)
         recipe(1:2, 10) = (/16, 25/)
         recipe(1:2, 11) = (/26, 28/)
         max_packs_nbands = 2
         packs_nbands_list = (/1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 2/)

      case default
         if (myrank .eq. 0) write(*,*) "[GPU:radsw] Not support for npacks = ", npacks
         stop
         
      end select

      do np = 1, npacks
         do nib = 1, nbdsw
            if (recipe(nib, np) .eq. 0) exit
            call add_ib_to_group(recipe(nib, np), np, nib)
         end do
      end do

      !if (myrank .eq. 0) then
      !   do np = 1, npacks
      !      write(*,*) pack_size(np), sum(pack_size)
      !      write(*,*) reverse_recipe(:,np)
      !   end do
      !end if

      return

      contains
      !-------------------------
      subroutine add_ib_to_group(ib, np, n)
      implicit none
      integer, intent(in) :: ib, np, n
      integer :: i
      
      do i = 1, ng(ib)
         pack_size(np)  = pack_size(np) + 1
         pack_idx(pack_size(np), np) = ib
         reverse_recipe(pack_size(np), np) = n
      end do


      end subroutine
      !-------------------------
      end subroutine

      subroutine copyin_radsw_main_gpu(async_id)
      ! must be called after executing rswinit_gpu
      implicit none
      
      integer, intent(in) :: async_id

      !$acc enter data copyin(nspa, nspb, idxebc, idxsfc, exp_tbl, ng, ngs, ngb, &
      !$acc&      wvnum1, wvnum2, pack_size, pack_idx, recipe, reverse_recipe, &
      !$acc&      packs_nbands_list) async(async_id)

      return
      end subroutine


!-----------------------------------
      subroutine swrad_gpu                                                  &
!...................................

!  ---  inputs for setaer_sw_gpu:
           ( prslk1, tvly, rhly, slmsk, tracer1, xlon, xlat, lsswr, lslwr, &
             me, map_jj, map_i, my_max, ntrac, &
!  ---  inputs:
     &       plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr_co2, gasvmr_other,                      &
     &       clouds,cfrac, cliqp, reliq, cicep, reice, cdat1, cdat2, cdat3, cdat4, &
             icseed,sfcalb,                             &
     &       cosz,solcon,nday,idxday, map_nday_ipt, map_nday_jj, offset_nday,                                   &
     &       myim, nlay, nlp1, lprnt, myrank, ix, kd,                            &
     &       nf_clds, nf_vgas, nf_albd, nf_aesw, max_nday_length, async_id, fulljj, blocks, smalljj, lhtrswb, &
!  ---  outputs:
     &       hswc,topflx_upfxc, topflx_dnfxc, topflx_upfx0, sfcflx_upfxc, &
         sfcflx_dnfxc, sfcflx_upfx0, sfcflx_dnfx0,                              &
!! ---  optional:
     &      hsw0,hswb,flxprf_upfxc, flxprf_dnfxc, flxprf_upfx0, flxprf_dnfx0, &
            fdncmp_uvbfc, fdncmp_uvbf0, fdncmp_nirbm, fdncmp_nirdf, fdncmp_visbm, fdncmp_visdf &
     &     )

!  ====================  defination of variables  ====================  !
!                                                                       !
!  input variables:                                                     !
!   plyr (npts,nlay) : model layer mean pressure in mb                  !
!   plvl (npts,nlp1) : model level pressure in mb                       !
!   tlyr (npts,nlay) : model layer mean temperature in k                !
!   tlvl (npts,nlp1) : model level temperature in k    (not in use)     !
!   qlyr (npts,nlay) : layer specific humidity in gm/gm   *see inside   !
!   olyr (npts,nlay) : layer ozone concentration in gm/gm               !
!   gasvmr(npts,nlay,:): atmospheric constent gases:                    !
!                      (check module_radiation_gases for definition)    !
!      gasvmr(:,:,1)  - co2 volume mixing ratio                         !
!      gasvmr(:,:,2)  - n2o volume mixing ratio                         !
!      gasvmr(:,:,3)  - ch4 volume mixing ratio                         !
!      gasvmr(:,:,4)  - o2  volume mixing ratio                         !
!      gasvmr(:,:,5)  - co  volume mixing ratio        (not used)       !
!      gasvmr(:,:,6)  - cfc11 volume mixing ratio      (not used)       !
!      gasvmr(:,:,7)  - cfc12 volume mixing ratio      (not used)       !
!      gasvmr(:,:,8)  - cfc22 volume mixing ratio      (not used)       !
!      gasvmr(:,:,9)  - ccl4  volume mixing ratio      (not used)       !
!   clouds(npts,nlay,:): cloud profile                                  !
!                      (check module_radiation_clouds for definition)   !
!                ---  for  iswcliq > 0  ---                             !
!       clouds(:,:,1)  -   layer total cloud fraction                   !
!       clouds(:,:,2)  -   layer in-cloud liq water path   (g/m**2)     !
!       clouds(:,:,3)  -   mean eff radius for liq cloud   (micron)     !
!       clouds(:,:,4)  -   layer in-cloud ice water path   (g/m**2)     !
!       clouds(:,:,5)  -   mean eff radius for ice cloud   (micron)     !
!       clouds(:,:,6)  -   layer rain drop water path      (g/m**2)     !
!       clouds(:,:,7)  -   mean eff radius for rain drop   (micron)     !
!       clouds(:,:,8)  -   layer snow flake water path     (g/m**2)     !
!       clouds(:,:,9)  -   mean eff radius for snow flake  (micron)     !
!                ---  for  iswcliq = 0  ---                             !
!       clouds(:,:,1)  -   layer total cloud fraction                   !
!       clouds(:,:,2)  -   layer cloud optical depth                    !
!       clouds(:,:,3)  -   layer cloud single scattering albedo         !
!       clouds(:,:,4)  -   layer cloud asymmetry factor                 !
!     icseed(npts)   : auxiliary special cloud related array            !
!                      when module variable isubcsw=2, it provides      !
!                      permutation seed for each column profile that    !
!                      are used for generating random numbers.          !
!                      when isubcsw /=2, it will not be used.           !
!   aerosols(npts,nlay,nbdsw,:) : aerosol optical properties            !
!                      (check module_radiation_aerosols for definition) !
!         (:,:,:,1)   - optical depth                                   !
!         (:,:,:,2)   - single scattering albedo                        !
!         (:,:,:,3)   - asymmetry parameter                             !
!   sfcalb(npts, : ) : surface albedo in fraction                       !
!                      (check module_radiation_surface for definition)  !
!         ( :, 1 )    - near ir direct beam albedo                      !
!         ( :, 2 )    - near ir diffused albedo                         !
!         ( :, 3 )    - uv+vis direct beam albedo                       !
!         ( :, 4 )    - uv+vis diffused albedo                          !
!   cosz  (npts)     : cosine of solar zenith angle                     !
!   solcon           : solar constant                      (w/m**2)     !
!   nday             : num of daytime points                            !
!   idxday(npts)     : index array for daytime points                   !
!   npts             : number of horizontal points                      !
!   nlay,nlp1        : vertical layer/lavel numbers                     !
!   lprnt            : logical check print flag                         !
!                                                                       !
!  output variables:                                                    !
!   hswc  (npts,nlay): total sky heating rates (k/sec or k/day)         !
!   topflx(npts)     : radiation fluxes at toa (w/m**2), components:    !
!                      (check module_radsw_parameters for definition)   !
!     upfxc            - total sky upward flux at toa                   !
!     dnflx            - total sky downward flux at toa                 !
!     upfx0            - clear sky upward flux at toa                   !
!   sfcflx(npts)     : radiation fluxes at sfc (w/m**2), components:    !
!                      (check module_radsw_parameters for definition)   !
!     upfxc            - total sky upward flux at sfc                   !
!     dnfxc            - total sky downward flux at sfc                 !
!     upfx0            - clear sky upward flux at sfc                   !
!     dnfx0            - clear sky downward flux at sfc                 !
!                                                                       !
!!optional outputs variables:                                           !
!   hswb(npts,nlay,nbdsw): spectral band total sky heating rates        !
!   hsw0  (npts,nlay): clear sky heating rates (k/sec or k/day)         !
!   flxprf(npts,nlp1): level radiation fluxes (w/m**2), components:     !
!                      (check module_radsw_parameters for definition)   !
!     dnfxc            - total sky downward flux at interface           !
!     upfxc            - total sky upward flux at interface             !
!     dnfx0            - clear sky downward flux at interface           !
!     upfx0            - clear sky upward flux at interface             !
!   fdncmp(npts)     : component surface downward fluxes (w/m**2):      !
!                      (check module_radsw_parameters for definition)   !
!     uvbfc            - total sky downward uv-b flux at sfc            !
!     uvbf0            - clear sky downward uv-b flux at sfc            !
!     nirbm            - downward surface nir direct beam flux          !
!     nirdf            - downward surface nir diffused flux             !
!     visbm            - downward surface uv+vis direct beam flux       !
!     visdf            - downward surface uv+vis diffused flux          !
!                                                                       !
!  external module variables:  (in physpara)                            !
!   iswrgas - control flag for rare gases (ch4,n2o,o2, etc.)            !
!           =0: do not include rare gases                               !
!           >0: include all rare gases                                  !
!   iswcliq - control flag for liq-cloud optical properties             !
!           =0: input cloud optical depth, fixed ssa, asy               !
!           =1: use hu and stamnes(1993) method for liq cld             !
!           =2: not used                                                !
!   iswcice - control flag for ice-cloud optical properties             !
!           *** if iswcliq==0, iswcice is ignored                       !
!           =1: use ebert and curry (1992) scheme for ice clouds        !
!           =2: use streamer v3.0 (2001) method for ice clouds          !
!           =3: use fu's method (1996) for ice clouds                   !
!   iswmode - control flag for 2-stream transfer scheme                 !
!           =1; delta-eddington    (joseph et al., 1976)                !
!           =2: pifm               (zdunkowski et al., 1980)            !
!           =3: discrete ordinates (liou, 1973)                         !
!   isubcsw - sub-column cloud approximation control flag               !
!           =0: no sub-col cld treatment, use grid-mean cld quantities  !
!           =1: mcica sub-col, prescribed seeds to get random numbers   !
!           =2: mcica sub-col, providing array icseed for random numbers!
!   iovrsw  - cloud overlapping control flag                            !
!           =0: random overlapping clouds                               !
!           =1: maximum/random overlapping clouds                       !
!           =2: maximum overlap cloud                                   !
!   ivflip  - control flg for direction of vertical index               !
!           =0: index from toa to surface                               !
!           =1: index from surface to toa                               !
!                                                                       !
!  module parameters, control variables:                                !
!     nblow,nbhgh      - lower and upper limits of spectral bands       !
!     maxgas           - maximum number of absorbing gaseous            !
!     ngptsw           - total number of g-point subintervals           !
!     ng##             - number of g-points in band (##=16-29)          !
!     ngb(ngptsw)      - band indices for each g-point                  !
!     bpade            - pade approximation constant (1/0.278)          !
!     nspa,nspb(nblow:nbhgh)                                            !
!                      - number of lower/upper ref atm's per band       !
!     ipsdsw0          - permutation seed for mcica sub-col clds        !
!                                                                       !
!  major local variables:                                               !
!     pavel  (nlay)         - layer pressures (mb)                      !
!     delp   (nlay)         - layer pressure thickness (mb)             !
!     tavel  (nlay)         - layer temperatures (k)                    !
!     coldry (nlay)         - dry air column amount                     !
!                                   (1.e-20*molecules/cm**2)            !
!     cldfrc (nlay)         - layer cloud fraction (norm by tot cld)    !
!     cldfmc (nlay,ngptsw)  - layer cloud fraction for g-point          !
!     taucw  (nlay,nbdsw)   - cloud optical depth                       !
!     ssacw  (nlay,nbdsw)   - cloud single scattering albedo (weighted) !
!     asycw  (nlay,nbdsw)   - cloud asymmetry factor         (weighted) !
!     tauaer (nlay,nbdsw)   - aerosol optical depths                    !
!     ssaaer (nlay,nbdsw)   - aerosol single scattering albedo          !
!     asyaer (nlay,nbdsw)   - aerosol asymmetry factor                  !
!     colamt (nlay,maxgas)  - column amounts of absorbing gases         !
!                             1 to maxgas are for h2o, co2, o3, n2o,    !
!                             ch4, o2, co, respectively (mol/cm**2)     !
!     facij  (nlay)         - indicator of interpolation factors        !
!                             =0/1: indicate lower/higher temp & height !
!     selffac(nlay)         - scale factor for self-continuum, equals   !
!                          (w.v. density)/(atm density at 296k,1013 mb) !
!     selffrac(nlay)        - factor for temp interpolation of ref      !
!                             self-continuum data                       !
!     indself(nlay)         - index of the lower two appropriate ref    !
!                             temp for the self-continuum interpolation !
!     forfac (nlay)         - scale factor for w.v. foreign-continuum   !
!     forfrac(nlay)         - factor for temp interpolation of ref      !
!                             w.v. foreign-continuum data               !
!     indfor (nlay)         - index of the lower two appropriate ref    !
!                             temp for the foreign-continuum interp     !
!     laytrop               - layer at which switch is made from one    !
!                             combination of key species to another     !
!     jp(nlay),jt(nlay),jt1(nlay)                                       !
!                           - lookup table indexes                      !
!     flxucb(nlp1,nbdsw)    - spectral bnd total-sky upward flx (w/m2)  !
!     flxdcb(nlp1,nbdsw)    - spectral bnd total-sky downward flx (w/m2)!
!     flxu0b(nlp1,nbdsw)    - spectral bnd clear-sky upward flx (w/m2)  !
!     flxd0b(nlp1,nbdsw)    - spectral b d clear-sky downward flx (w/m2)!
!                                                                       !
!                                                                       !
!  =====================    end of definitions    ====================  !

!  ---  inputs:
      integer,  intent(in) :: me, my_max, ntrac, max_nday_length
      logical,  intent(in) :: lsswr, lslwr
      real (kind=kind_phys), dimension(ix, my_max), intent(in) ::  slmsk,  &
               xlon, xlat
      real (kind=kind_phys), dimension(nxptot,nlay)  :: rhly, prslk1, tvly
      real (kind=kind_phys), dimension(nxptot,nlay,ntrac)   :: tracer1
      integer, dimension(nxptot) :: map_jj, map_i


      integer, intent(in) :: myim(fulljj), nlay, nlp1, nday(fulljj), myrank, ix, &
         nf_clds, nf_vgas, nf_albd, nf_aesw, fulljj, blocks, smalljj, kd
      integer, dimension(fulljj) :: offset_nday
      integer, dimension(ix*fulljj) :: map_nday_ipt, map_nday_jj

      integer, dimension(nxptot), intent(in) :: idxday
      integer, dimension(ix, fulljj), intent(in) :: icseed

      logical, intent(in) :: lprnt, lhtrswb

      real (kind=kind_phys), dimension(nxptot,nlp1), intent(in) ::        &
     &       plvl, tlvl
      real (kind=kind_phys), dimension(nxptot,nlay), intent(in) ::        &
     &       plyr, tlyr, qlyr, olyr
      real (kind=kind_phys), dimension(nxptot,nf_albd),    intent(in) :: sfcalb

      real (kind=kind_phys), dimension(nxptot,nlay),intent(in):: gasvmr_co2
      real (kind=kind_phys), dimension(nf_vgas),intent(in):: gasvmr_other
      real (kind=kind_phys), dimension(nxptot,nlay,nf_clds),intent(in):: clouds
      real (kind=kind_phys), dimension(nxptot,nlay),intent(in):: &
         cliqp, reliq, cicep, reice, cdat4
      real (kind=kind_phys), dimension(nxptot,nlay),intent(inout):: cfrac
      real (kind=kind_phys), dimension(nxptot,nlay),intent(inout):: cdat1, cdat2, cdat3

!      real (kind=kind_phys), dimension(ix,nlay,nbdsw,nf_aesw, fulljj),intent(in)::  &
!     &       aerosols


      real (kind=kind_phys), intent(in) :: cosz(ix, fulljj), solcon

!  ---  outputs:
      real (kind=kind_phys), dimension(ix,nlay, fulljj), intent(out) :: hswc

      !real (kind=kind_phys),    dimension(ix, fulljj, 3), intent(out) :: topflx
      real (kind=kind_phys),    dimension(ix, fulljj), intent(out) :: topflx_upfxc, &
         topflx_upfx0, topflx_dnfxc
      real (kind=kind_phys),    dimension(ix, fulljj), intent(out) :: sfcflx_upfxc, &
         sfcflx_dnfxc, sfcflx_upfx0, sfcflx_dnfx0

!! ---  optional outputs:
      real (kind=kind_phys), dimension(ix,nlay,nbdsw, fulljj), optional,      &
     &       intent(out) :: hswb

      real (kind=kind_phys), dimension(ix,nlay, fulljj),       optional,      &
     &       intent(out) :: hsw0
      real (kind=kind_phys),    dimension(ix,nlp1, fulljj),       optional,      &
     &       intent(out) :: flxprf_upfxc, flxprf_dnfxc, flxprf_upfx0, flxprf_dnfx0
      real (kind=kind_phys),    dimension(ix, fulljj),            optional,      &
     &       intent(out) :: fdncmp_uvbfc, fdncmp_uvbf0, fdncmp_nirbm, fdncmp_nirdf, &
        fdncmp_visbm, fdncmp_visdf

!  ---  locals:
      real (kind=kind_phys), dimension(:,:,:,:), allocatable :: cldfmc
      real (kind=kind_phys), dimension(max_nday_length, nlp1, nbdsw) :: fxupc, fxdnc,      &
     &       fxup0, fxdn0


      real (kind=kind_phys), dimension(max_nday_length, nlay)   :: delp,     &
     &       pavel, tavel, coldry, colmol, h2ovmr, o3vmr, temcol, rfdelp

      real (kind=kind_phys), dimension(max_nday_length, nlp1) :: fnet, flxdc, flxuc,     &
     &       flxd0, flxu0
      real (kind=kind_phys), dimension(max_nday_length, nlp1, nbdsw) :: fnetm

      real (kind=kind_phys), dimension(max_nday_length, 2) :: albbm, albdf, sfbmc,       &
     &       sfbm0, sfdfc, sfdf0

      real (kind=kind_phys) :: cosz1, tem0, tem1, tem2, s0fac
      real (kind=kind_phys), dimension(max_nday_length) :: ssolar, zcf0, zcf1, ftoau0, &
         ftoauc, ftoadc, fsfcu0, fsfcuc, fsfcd0, fsfcdc, suvbfc, suvbf0, sntz1

!  ---  column amount of absorbing gases:
!       (:,m) m = 1-h2o, 2-co2, 3-o3, 4-n2o, 5-ch4, 6-o2, 7-co
      real (kind=kind_phys), dimension(max_nday_length, nlay, maxgas) ::  colamt

      integer, dimension(ix, smalljj) :: ipseed

      integer :: i, ib, ipt, j1, k, kk, mb, jj
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real (kind=kind_phys), dimension(nlay)   :: cldfrc_im,     delp_im,     &
     &       pavel_im, tavel_im, coldry_im, colmol_im, h2ovmr_im, o3vmr_im, temcol_im,       &
     &       cliqp_im, reliq_im, cicep_im, reice_im, cdat1_im, cdat2_im, cdat3_im, cdat4_im,    &
     &       cfrac_im, fac00_im, fac01_im, fac10_im, fac11_im, forfac_im, forfrac_im,        &
     &       selffac_im, selffrac_im, rfdelp_im
      real (kind=kind_phys), dimension(2) :: albbm_im, albdf_im, sfbmc_im,       &
     &       sfbm0_im, sfdfc_im, sfdf0_im
      real (kind=kind_phys), dimension(nlay,nbdsw)  ::                  &
     &       tauae_im, ssaae_im, asyae_im, taucw_im, ssacw_im, asycw_im
      real (kind=kind_phys), dimension(nlay,ngptsw) :: cldfmc_im,          &
     &       taug_im, taur_im
      real (kind=kind_phys), dimension(nlp1,nbdsw):: fxupc_im, fxdnc_im,      &
     &       fxup0_im, fxdn0_im
      integer, dimension(nlay) :: indfor_im, indself_im, jp_im, jt_im, jt1_im
      real (kind=kind_phys), dimension(ngptsw) :: sfluxzen_im
      integer, dimension(ix) :: ipseed_im
      real (kind=kind_phys) ::  colamt_im(nlay,maxgas)
      integer :: async_id, jf, jb,  jjoffset, jbs, jbe, blocks_local, blockjj, &
         nday_length, jbe_nday, jbs_nday, n, m, r


      ! GPU: variable name changed: CPU - topflx%upfxc, GPU - topflx_upfxc
      ! GPU: variable name changed: CPU - topflx%dnfxc, GPU - topflx_dnfxc
      ! GPU: variable name changed: CPU - topflx%upfx0, GPU - topflx_upfx0
      ! GPU: variable name changed: CPU - sfcflx%upfxc, GPU - sfcflx_upfxc
      ! GPU: variable name changed: CPU - sfcflx%dnfxc, GPU - sfcflx_dnfxc
      ! GPU: variable name changed: CPU - sfcflx%upfx0, GPU - sfcflx_upfx0
      ! GPU: variable name changed: CPU - sfcflx%dnfx0, GPU - sfcflx_dnfx0
      ! GPU: variable name changed: CPU - flxprf%upfxc, GPU - flxprf_upfxc
      ! GPU: variable name changed: CPU - flxprf%dnfxc, GPU - flxprf_dnfxc
      ! GPU: variable name changed: CPU - flxprf%upfx0, GPU - flxprf_upfx0
      ! GPU: variable name changed: CPU - flxprf%dnfx0, GPU - flxprf_dnfx0
      ! GPU: variable name changed: CPU - fdncmp%uvbfc, GPU - fdncmp_uvbfc
      ! GPU: variable name changed: CPU - fdncmp%uvbf0, GPU - fdncmp_uvbf0
      ! GPU: variable name changed: CPU - fdncmp%nirbm, GPU - fdncmp_nirbm
      ! GPU: variable name changed: CPU - fdncmp%nirdf, GPU - fdncmp_nirdf
      ! GPU: variable name changed: CPU - fdncmp%visbm, GPU - fdncmp_visbm
      ! GPU: variable name changed: CPU - fdncmp%visdf, GPU - fdncmp_visdf
      ! GPU: variable name changed: CPU - aerosols(:,:,:,1), GPU - tauae1
      ! GPU: variable name changed: CPU - aerosols(:,:,:,2), GPU - ssaae1
      ! GPU: variable name changed: CPU - aerosols(:,:,:,3), GPU - asyae1
!
!===> ... begin here
!

      if (isubcsw > 0) allocate(cldfmc(ix, nlay,ngptsw, smalljj))

      lhswb  = present ( hswb )
      lhsw0  = present ( hsw0 )
      lflxprf_upfxc= present ( flxprf_upfxc )
      lflxprf_dnfxc= present ( flxprf_dnfxc )
      lflxprf_upfx0= present ( flxprf_upfx0 )
      lflxprf_dnfx0= present ( flxprf_upfxc )
      lflxprf = ((lflxprf_upfxc .and. lflxprf_dnfxc) .and. &
                 (lflxprf_upfx0 .and. lflxprf_dnfx0))
      
      lfdncmp_uvbfc= present ( fdncmp_uvbfc )
      lfdncmp_uvbf0= present ( fdncmp_uvbf0 )
      lfdncmp_nirbm= present ( fdncmp_nirbm )
      lfdncmp_nirdf= present ( fdncmp_nirdf )
      lfdncmp_visbm= present ( fdncmp_visbm )
      lfdncmp_visdf= present ( fdncmp_visdf )
      lfdncmp= (((lfdncmp_uvbfc .and. lfdncmp_uvbf0) .and. &
                 (lfdncmp_nirbm .and. lfdncmp_nirdf)) .and. &
                (lfdncmp_visbm .and. lfdncmp_visdf))
 
!      if (myrank .eq. 0) print *,'$$$$ swrad start $$$$'
!      if (myrank .eq. 0) print *,'$$$$ lhswb   =',lhswb
!      if (myrank .eq. 0) print *,'$$$$ lhsw0   =',lhsw0
!      if (myrank .eq. 0) print *,'$$$$ lflxprf =',lflxprf
!      if (myrank .eq. 0) print *,'$$$$ lfdncmp =',lfdncmp
!      if (myrank .eq. 0) print *,'$$$$ nlp1    =',nlp1

!  --- ...  compute solar constant adjustment factor according to solcon.
!      ***  s0, the solar constant at toa in w/m**2, is hard-coded with
!           each spectra band, the total flux is about 1368.22 w/m**2.

      s0fac = solcon / s0

         !  --- ...  initial output arrays
      !$acc data create(fxupc, fxdnc, fxup0, fxdn0, delp, pavel, &
      !$acc&     tavel, coldry, colmol, h2ovmr, o3vmr, temcol, rfdelp, fnet, &
      !$acc&     flxdc, flxuc, flxd0, flxu0, fnetm, albbm, albdf, sfbmc, sfbm0, &
      !$acc&     sfdfc, sfdf0, sntz1, ssolar, zcf0, zcf1, ftoau0, ftoauc, ftoadc, &
      !$acc&     fsfcu0, fsfcuc, fsfcd0, fsfcdc, suvbfc, suvbf0, colamt, ipseed) async(async_id)
      do jb = 1, blocks
      jjoffset = (fulljj-1)*(jb-1)/blocks
      jbs = jjoffset+1
      jbe = (fulljj-1)*jb/blocks
      blockjj = jbe - jbs + 1
      jbs_nday = offset_nday(jbs)+1
      jbe_nday = offset_nday(jbe+1)
      nday_length = jbe_nday - jbs_nday + 1
      !$acc parallel loop gang collapse(2) private(jf) async(async_id)
      do jj = 1, blockjj
         do k = 1, nlay
            jf = jjoffset+jj
            !$acc loop vector
            do i = 1, myim(jf)
               hswc(i, k, jf) = f_zero
            end do
         end do
      end do
      !$acc parallel loop collapse(2) private(jf) async(async_id)
      do jj = 1, blockjj
         do i = 1, ix
            jf = jjoffset+jj
            if (i .le. myim(jf)) then
               topflx_upfxc(i, jf) = f_zero
               topflx_upfx0(i, jf) = f_zero
               topflx_dnfxc(i, jf) = f_zero
               sfcflx_upfxc(i, jf) = f_zero
               sfcflx_dnfxc(i, jf) = f_zero
               sfcflx_upfx0(i, jf) = f_zero
               sfcflx_dnfx0(i, jf) = f_zero
            end if
         end do
      end do

         !! --- ...  initial optional outputs
      if ( lflxprf ) then
         !$acc parallel loop gang collapse(2) private(jf) async(async_id)
         do jj = 1, blockjj
            do k = 1, nlp1
               jf = jjoffset+jj
               !$acc loop vector
               do i = 1, myim(jf)
                  flxprf_upfxc(i, k, jf) = f_zero
                  flxprf_dnfxc(i, k, jf) = f_zero
                  flxprf_upfx0(i, k, jf) = f_zero
                  flxprf_dnfx0(i, k, jf) = f_zero
               end do
            end do
         end do
      endif

      if ( lfdncmp ) then
         !$acc parallel loop collapse(2) private(jf) async(async_id)
         do jj = 1, blockjj
            do i = 1, ix
               jf = jjoffset+jj
               if (i .le. myim(jf)) then
                  fdncmp_uvbfc(i, jf) = f_zero
                  fdncmp_uvbf0(i, jf) = f_zero
                  fdncmp_nirbm(i, jf) = f_zero
                  fdncmp_nirdf(i, jf) = f_zero
                  fdncmp_visbm(i, jf) = f_zero
                  fdncmp_visdf(i, jf) = f_zero
               end if
            end do
         end do
      endif

      if ( lhsw0 ) then
         !$acc parallel loop gang collapse(2) private(jf) async(async_id)
         do jj = 1, blockjj
            do k = 1, nlay
               jf = jjoffset+jj
               !$acc loop vector
               do i = 1, myim(jf)
                  hsw0(i,k, jf) = f_zero
               end do
            end do
         end do
      endif

      ! GPU: hswb has already initialized just before calling swrad
      !if ( lhswb ) then
      !   !$acc parallel loop gang collapse(3) private(jf) async(async_id)
      !   do jj = 1, blockjj
      !      do ib = 1, nbdsw
      !         do k = 1, nlay
      !            jf = jjoffset+jj
      !            !$acc loop vector
      !            do i = 1, myim(jf)
      !               hswb(i, k, ib, jf) = f_zero
      !            end do
      !         end do
      !      end do
      !   end do
      !endif

      !  --- ...  change random number seed value for each radiation invocation

      if     ( isubcsw == 1 ) then     ! advance prescribed permutation seed
         !$acc parallel loop collapse(2) private(jf) async(async_id)
         do jj = 1, blockjj
            do i = 1, ix
               jf = jjoffset+jj
               if (i .le. myim(jf)) then
                  ipseed(i, jj) = ipsdsw0 + i
               end if
            enddo
         end do
      elseif ( isubcsw == 2 ) then     ! use input array of permutaion seeds
         !$acc parallel loop collapse(2) private(jf) async(async_id)
         do jj = 1, blockjj
            do i = 1, ix
               jf = jjoffset+jj
               if (i .le. myim(jf)) then
                  ipseed(i, jj) = icseed(i, jf)
               end if
            enddo
         end do
      endif

         !     if ( lprnt ) then
         !       print *,'  in radsw, isubcsw, ipsdsw0,ipseed =',                &
         !    &           isubcsw, ipsdsw0, ipseed
         !     endif

         !  --- ...  loop over each daytime grid point
      !$acc parallel loop private(j1, jf, jj, m, r) async(async_id)
      do n = jbs_nday, nday_length + jbs_nday - 1
         jj = map_nday_jj(n) - jjoffset
         jf = jjoffset+jj
         j1 = idxday(n)
         m = n - jbs_nday + 1
         r = j1 + nxjp_acc(jf) - 1
         sntz1(m)  = f_one / cosz(j1, jf)
         ssolar(m) = s0fac * cosz(j1, jf)

         !  --- ...  surface albedo: bm,df - dir,dif;  1,2 - nir,uvv
         albbm(m, 1) = sfcalb(r,1)
         albdf(m, 1) = sfcalb(r,2)
         albbm(m, 2) = sfcalb(r,3)
         albdf(m, 2) = sfcalb(r,4)
      end do

            !  --- ...  prepare atmospheric profile for use in rrtm
            !           the vertical index of internal array is from surface to top

      if (ivflip == 0) then       ! input from toa to sfc

         tem1 = 100.0 * con_g
         tem2 = 1.0e-20 * 1.0e3 * con_avgd
         !$acc parallel loop collapse(2) private(jf, j1, kk, tem0, jj, m, r) async(async_id)
         do k = 1, nlay
            do n = jbs_nday, nday_length + jbs_nday - 1
               jj = map_nday_jj(n) - jjoffset
               jf = jjoffset+jj
               j1 = idxday(n)
               m = n - jbs_nday + 1
               r = j1 + nxjp_acc(jf) - 1
               kk = nlp1 - k
               pavel(m, k) = plyr(r,kk)
               tavel(m, k) = tlyr(r,kk)
               delp (m, k) = plvl(r,kk+1) - plvl(r,kk)

               !  --- ...  set absorber amount
               !test use
               !           h2ovmr(k)= max(f_zero,qlyr(j1,kk)*amdw)                     ! input mass mixing ratio
               !           h2ovmr(k)= max(f_zero,qlyr(j1,kk))                          ! input vol mixing ratio
               !           o3vmr (k)= max(f_zero,olyr(j1,kk))                          ! input vol mixing ratio
               !ncep model use
               h2ovmr(m, k)= max(f_zero,qlyr(r,kk)*amdw/(f_one-qlyr(r,kk))) ! input specific humidity
               o3vmr (m, k)= max(f_zero,olyr(r,kk)*amdo3)                    ! input mass mixing ratio

               tem0 = (f_one - h2ovmr(m, k))*con_amd + h2ovmr(m, k)*con_amw
               coldry(m, k) = tem2 * delp(m, k) / &
               (tem1*tem0*(f_one + h2ovmr(m, k)))
               temcol(m, k) = 1.0e-12 * coldry(m, k)

               colamt(m, k,1) = max(f_zero,    coldry(m, k)*h2ovmr(m, k))         ! h2o
               colamt(m, k,2) = max(temcol(m, k), coldry(m, k)*gasvmr_co2(r,kk))   ! co2
               colamt(m, k,3) = max(f_zero,    coldry(m, k)*o3vmr(m, k))          ! o3
               colmol(m, k)   = coldry(m, k) + colamt(m, k,1)
            end do
         end do

               !  --- ...  set up gas column amount, convert from volume mixing ratio
               !           to molec/cm2 based on coldry (scaled to 1.0e-20)

         if (iswrgas > 0) then
            !$acc parallel loop collapse(2) private(jj, m) async(async_id)
            do k = 1, nlay
               do n = jbs_nday, nday_length + jbs_nday - 1
                  jj = map_nday_jj(n) - jjoffset
                  m = n - jbs_nday + 1
                  colamt(m, k,4) = max(temcol(m, k), coldry(m, k)*gasvmr_other(2))  ! n2o
                  colamt(m, k,5) = max(temcol(m, k), coldry(m, k)*gasvmr_other(3))  ! ch4
                  colamt(m, k,6) = max(temcol(m, k), coldry(m, k)*gasvmr_other(4))  ! o2
                  !             colamt(k,7) = max(temcol(k), coldry(k)*gasvmr(j1,kk,5))  ! co - notused
               end do
            end do
         else
            !$acc parallel loop collapse(2) private(jj, m) async(async_id)
            do k = 1, nlay
               do n = 1, nday_length
                  jj = map_nday_jj(n) - jjoffset
                  m = n - jbs_nday + 1
                  colamt(m, k,4) = temcol(m, k)                                  ! n2o
                  colamt(m, k,5) = temcol(m, k)                                  ! ch4
                  colamt(m, k,6) = temcol(m, k)                                  ! o2
                  !             colamt(k,7) = temcol(k)                                  ! co - notused
               end do
            end do
         endif

               !  --- ...  set aerosol optical properties
         if (iswcliq <= 0) then    ! use prognostic cloud method
            !$acc parallel loop collapse(2) private(kk, jf, jj) async(async_id)
            do k = 1, nlay
               do n = jbs_nday, nday_length + jbs_nday - 1
                  jj = map_nday_jj(n) - jjoffset
                  kk = nlp1 - k
                  jf = jjoffset+jj
                  j1 = idxday(n)
                  r = j1 + nxjp_acc(jf) - 1
                  cfrac(n, k) = clouds(r,kk,1)      ! cloud fraction
                  cdat1(n, k) = clouds(r,kk,2)      ! cloud optical depth
                  cdat2(n, k) = clouds(r,kk,3)      ! cloud single scattering albedo
                  cdat3(n, k) = clouds(r,kk,4)      ! cloud asymmetry factor
               end do
            end do
         endif                    ! end if_iswcliq
      end if
      if (ivflip .ne. 0) then                        ! input from sfc to toa

         tem1 = 100.0 * con_g
         tem2 = 1.0e-20 * 1.0e3 * con_avgd
         !$acc parallel loop collapse(2) private(jf, j1, tem0, jj, m, r) async(async_id)
         do k = 1, nlay
            do n = jbs_nday, nday_length + jbs_nday - 1
               jj = map_nday_jj(n) - jjoffset
               jf = jjoffset+jj
               j1 = idxday(n)
               m = n - jbs_nday + 1
               r = j1 + nxjp_acc(jf) - 1
               pavel(m, k) = plyr(r,k)
               tavel(m, k) = tlyr(r,k)
               delp (m, k) = plvl(r,k) - plvl(r,k+1)

               !  --- ...  set absorber amount
               !test use
               !           h2ovmr(k)= max(f_zero,qlyr(j1,k)*amdw)                    ! input mass mixing ratio
               !           h2ovmr(k)= max(f_zero,qlyr(j1,k))                         ! input vol mixing ratio
               !           o3vmr (k)= max(f_zero,olyr(j1,k))                         ! input vol mixing ratio
               !ncep model use
               h2ovmr(m, k)= max(f_zero,qlyr(r,k)*amdw/(f_one-qlyr(r,k))) ! input specific humidity
               o3vmr (m, k)= max(f_zero,olyr(r,k)*amdo3)                   ! input mass mixing ratio

               tem0 = (f_one - h2ovmr(m, k))*con_amd + h2ovmr(m, k)*con_amw
               coldry(m, k) = tem2 * delp(m, k) / (tem1*tem0*(f_one + h2ovmr(m, k)))
               temcol(m, k) = 1.0e-12 * coldry(m, k)

               colamt(m, k,1) = max(f_zero,    coldry(m, k)*h2ovmr(m, k))         ! h2o
               colamt(m, k,2) = max(temcol(m, k), coldry(m, k)*gasvmr_co2(r,k))    ! co2
               colamt(m, k,3) = max(f_zero,    coldry(m, k)*o3vmr(m, k))          ! o3
               colmol(m, k)   = coldry(m, k) + colamt(m, k,1)
            end do
         end do

         !  --- ...  set up gas column amount, convert from volume mixing ratio
         !           to molec/cm2 based on coldry (scaled to 1.0e-20)

         if (iswrgas > 0) then
            !$acc parallel loop collapse(2) private(jj, m) async(async_id)
            do k = 1, nlay
               do n = jbs_nday, nday_length + jbs_nday - 1
                  jj = map_nday_jj(n) - jjoffset
                  m = n - jbs_nday + 1
                  colamt(m, k,4) = max(temcol(m, k), coldry(m, k)*gasvmr_other(2))  ! n2o
                  colamt(m, k,5) = max(temcol(m, k), coldry(m, k)*gasvmr_other(3))  ! ch4
                  colamt(m, k,6) = max(temcol(m, k), coldry(m, k)*gasvmr_other(4))  ! o2
                  !             colamt(k,7) = max(temcol(k), coldry(k)*gasvmr(j1,k,5))  ! co - notused
               end do
            end do
         else
            !$acc parallel loop collapse(2) private(jj, m) async(async_id)
            do k = 1, nlay
               do n = jbs_nday, nday_length + jbs_nday - 1
                  jj = map_nday_jj(n) - jjoffset
                  m = n - jbs_nday + 1
                  colamt(m, k,4) = temcol(m, k)                                 ! n2o
                  colamt(m, k,5) = temcol(m, k)                                 ! ch4
                  colamt(m, k,6) = temcol(m, k)                                 ! o2
                  !             colamt(k,7) = temcol(k)                                 ! co - notused
               end do
            end do
         endif

               !  --- ...  set aerosol optical properties
         if (iswcliq <= 0) then    ! use prognostic cloud method
            !$acc parallel loop collapse(2) private(jf, jj) async(async_id)
            do k = 1, nlay
               do n = jbs_nday, nday_length + jbs_nday - 1
                  jj = map_nday_jj(n) - jjoffset
                  jf = jjoffset+jj
                  j1 = idxday(n)
                  r = j1 + nxjp_acc(jf) - 1
                  cfrac(n, k) = clouds(r,k,1)       ! cloud fraction
                  cdat1(n, k) = clouds(r,k,2)       ! cloud optical depth
                  cdat2(n, k) = clouds(r,k,3)       ! cloud single scattering albedo
                  cdat3(n, k) = clouds(r,k,4)       ! cloud asymmetry factor
               end do
            end do
         endif                    ! end if_iswcliq
      endif                       ! if_ivflip
      !!$acc exit data delete(tauae1, ssaae1, asyae1) async(async_id)

            !  --- ...  compute fractions of clear sky view
      !$acc parallel loop private(jj, m) async(async_id)
      do n = jbs_nday, nday_length + jbs_nday - 1
         jj = map_nday_jj(n) - jjoffset
         m = n - jbs_nday + 1
         zcf0(m)   = f_one
         zcf1(m)   = f_one
      end do
      if (iovrsw == 0) then                    ! random overlapping
         !$acc parallel loop private(jf, j1, jj, m, r) async(async_id)
         do n = jbs_nday, nday_length + jbs_nday - 1
            jj = map_nday_jj(n) - jjoffset
            jf = jjoffset+jj
            j1 = idxday(n)
            m = n - jbs_nday + 1
            r = j1 + nxjp_acc(jf) - 1
            !$acc loop seq
            do k = 1, nlay
               zcf0(m) = zcf0(m) * (f_one - cfrac(r, k))
            enddo
         end do
      else if (iovrsw == 1) then               ! max/ran overlapping
         !$acc parallel loop private(jf, j1, jj, m, r) async(async_id)
         do n = jbs_nday, nday_length + jbs_nday - 1
            jj = map_nday_jj(n) - jjoffset
            jf = jjoffset+jj
            j1 = idxday(n)
            m = n - jbs_nday + 1
            r = j1 + nxjp_acc(jf) - 1
            !$acc loop seq
            do k = 1, nlay
               if (cfrac(r, k) > ftiny) then                ! cloudy layer
                  zcf1(m) = min ( zcf1(m), f_one-cfrac(r, k) )
               elseif (zcf1(m) < f_one) then                ! clear layer
                  zcf0(m) = zcf0(m) * zcf1(m)
                  zcf1(m) = f_one
               endif
            enddo
         end do
         !$acc parallel loop private(m) async(async_id)
         do n = jbs_nday, nday_length + jbs_nday - 1
            m = n - jbs_nday + 1
            zcf0(m) = zcf0(m) * zcf1(m)
         end do
      else if (iovrsw == 2) then               ! maximum overlapping
         !$acc parallel loop private(jf, j1, jj, m, r) async(async_id)
         do n = jbs_nday, nday_length + jbs_nday - 1
            jj = map_nday_jj(n) - jjoffset
            jf = jjoffset+jj
            j1 = idxday(n)
            m = n - jbs_nday + 1
            r = j1 + nxjp_acc(jf) - 1
            !$acc loop seq
            do k = 1, nlay
               zcf0(m) = min ( zcf0(m), f_one-cfrac(r, k) )
            enddo
         end do
      endif
      !$acc parallel loop private(m) async(async_id)
      do n = jbs_nday, nday_length + jbs_nday - 1
         m = n - jbs_nday + 1
         if (zcf0(m) <= ftiny) zcf0(m) = f_zero
         if (zcf0(m) > oneminus) zcf0(m) = f_one
         zcf1(m) = f_one - zcf0(m)
      end do

            !  --- ...  compute cloud optical properties
      !call nvtxStartRange("sw_cldprop")
      !call cldprop                                                  &
      !   !  ---  inputs:
      !   &     ( cfrac,cliqp,reliq, &
      !           cicep,reice,cdat1, &
      !           cdat2,cdat3,cdat4,     &
      !   &       zcf1, nlay, ipseed, ix, nday(jbs:jbe), idxday(jbs_nday:jbe_nday), &
      !           map_nday_ipt(jbs_nday:jbe_nday), map_nday_jj(jbs_nday:jbe_nday), &
      !           nday_length, max_nday_length, jjoffset, jbs_nday, &
      !           async_id, smalljj, blockjj,                                    &
      !   !  ---  outputs:
      !   &       taucw, ssacw, asycw, cldfrc, cldfmc                        &
      !   &     )
      !call nvtxEndRange
      if (isubcsw > 0) then
         !$acc parallel loop gang collapse(2) private(jj, ipt, m) async(async_id)
         do k = 1, nlay
            do n = jbs_nday, nday_length + jbs_nday - 1
               jj = map_nday_jj(n) - jjoffset
               ipt = map_nday_ipt(n)
               m = n - jbs_nday + 1
               if (zcf1(m) <= 0.) then     ! clear sky column
                  !$acc loop seq
                  do ib = 1, ngptsw
                     cldfmc(ipt, k,ib, jj)= f_zero
                  end do
               endif   ! end if_zcf1_block
            end do
         end do
      end if
      !call nvtxStartRange("sw_setcoef")
      !call setcoef                                                    &
      !   !  ---  inputs:
      !   &     ( pavel,tavel,h2ovmr, nlay,nlp1, ix, nday(jbs:jbe), async_id, smalljj, blockjj,                              &
      !           map_nday_ipt(jbs_nday:jbe_nday), map_nday_jj(jbs_nday:jbe_nday), &
      !           nday_length, max_nday_length, jjoffset, &
      !   !  ---  outputs:
      !   &       laytrop,jp,jt,jt1,fac00,fac01,fac10,fac11,                 &
      !   &       selffac,selffrac,indself,forfac,forfrac,indfor             &
      !   &     )
      !call nvtxEndRange

         !  --- ...  calculate optical depths for gaseous absorption and rayleigh
         !           scattering
      !call nvtxStartRange("sw_taumol")
      !call taumol                                                     &
      !   !  ---  inputs:
      !   &     ( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
      !   &       forfac,forfrac,indfor,selffac,selffrac,indself, &
      !           nlay, ix, nday(jbs:jbe), async_id, smalljj, blockjj,     &
      !   !  ---  outputs:
      !   &       sfluxzen, taug, taur                                       &
      !   &     )
      !call nvtxEndRange
      !call nvtxStartRange("sw_setaer")
      !call setaer_sw_gpu                                                       &
      !!  ---  inputs:
      !&     ( plvl,plyr,prslk1,tvly,rhly,slmsk,tracer1, &
      !xlon,xlat,        &
      !&       myim,nlay,nlp1,lsswr,lslwr,me,myrank, ix, map_jj, map_i, nxptot, &
      !        idxday(jbs_nday:jbe_nday), map_nday_jj(jbs_nday:jbe_nday), jbs_nday, &
      !        nday_length, jjoffset, max_nday_length, async_id,                          &
      !!  ---  outputs:
      !&       tauae, ssaae, asyae                                              &
      !!    &       faersw,faerlw,aerodp                                       &
      !&     )
      !call nvtxEndRange

         !  --- ...  call the 2-stream radiation transfer model
      if ( isubcsw <= 0 ) then     ! use standard cloud scheme
         !call nvtxStartRange("sw_spcvrtc")
         call spcvrtc                                                  &
            !  ---  inputs:
                  ( plvl, plyr, prslk1, tvly, rhly, slmsk, tracer1, xlon, xlat, &
                    myim, lsswr, lslwr, me, map_jj, map_i, ntrac, my_max, &
                    cfrac, cliqp, reliq, cicep, reice, cdat1, cdat2, cdat3, cdat4, ipseed, &
                    colamt,colmol, pavel, tavel, h2ovmr, &
            &       ssolar,cosz(:,jbs:jbe),sntz1,albbm, &
                    albdf,            &
            &       zcf1,zcf0, &
            &       nlay, nlp1, ix, idxday(jbs_nday:jbe_nday), nday(jbs:jbe), &
                    map_nday_ipt(jbs_nday:jbe_nday), map_nday_jj(jbs_nday:jbe_nday), &
                    nday_length, max_nday_length, jjoffset, jbs_nday, &
                    async_id, smalljj, blockjj,                                    &
            !  ---  outputs:
            &       fxupc,fxdnc,fxup0,fxdn0,                                   &
            &       ftoauc,ftoau0,ftoadc, &
                     fsfcuc,fsfcu0,fsfcdc,fsfcd0,          &
            &       sfbmc,sfdfc,sfbm0,sfdf0,suvbfc,suvbf0                      &
            &     )
         !call nvtxEndRange

      else                         ! use mcica cloud scheme
         !!!! GPU: spcvrtm is not ported to GPU. the indices of arguments should 
         !!!! GPU: be revisited.
         !do jj = 1, blockjj
         !   jf = jjoffset+jj
         !   do ipt = 1, nday(jf) ! lab_do_ipt
         !      j1 = idxday(ipt, jf)
         !      do k = 1, nlay
         !         do i = 1, ngptsw
         !            cldfmc_im(k, i) = cldfmc(ipt, k, i, jj)
         !            taug_im(k, i) = taug(ipt, k, i, jj)
         !            taur_im(k, i) = taur(ipt, k, i, jj)
         !         end do
         !      end do
         !      do i = 1, ngptsw
         !         sfluxzen_im(i) = sfluxzen(ipt, i, jj)
         !      end do
!
         !      call spcvrtm                                                  &
         !         !  ---  inputs:
         !         &     ( ssolar(ipt, jj),cosz(j1, jf),sntz1(ipt, jj),albbm_im, &
         !                 albdf_im,sfluxzen_im,cldfmc_im,            &
         !         &       zcf1(ipt, jj),zcf0(ipt, jj),taug_im,taur_im,tauae_im, &
         !                 ssaae_im,asyae_im,taucw_im,ssacw_im,asycw_im,   &
         !         &       nlay, nlp1,                                                &
         !         !  ---  outputs:
         !         &       fxupc_im,fxdnc_im,fxup0_im,fxdn0_im,                                   &
         !         &       ftoauc(ipt, jj),ftoau0(ipt, jj),ftoadc(ipt, jj), &
         !                 fsfcuc(ipt, jj),fsfcu0(ipt, jj),fsfcdc(ipt, jj),fsfcd0(ipt, jj),          &
         !         &       sfbmc_im,sfdfc_im,sfbm0_im,sfdf0_im,suvbfc(ipt, jj),suvbf0(ipt, jj)                      &
         !         &     )
         !      do k = 1, nlay
         !         do i = 1, ngptsw
         !            cldfmc(ipt, k, i, jj) = cldfmc_im(k, i)
         !            taug(ipt, k, i, jj) = taug_im(k, i)
         !            taur(ipt, k, i, jj) = taur_im(k, i)
         !         end do
         !      end do
         !      do i = 1, ngptsw
         !         sfluxzen(ipt, i, jj) = sfluxzen_im(i)
         !      end do
         !   end do
         !end do
      endif

         !  --- ...  sum up total spectral fluxes for total-sky
      !$acc parallel loop collapse(2) private(m) async(async_id)
      do k = 1, nlp1
         do n = jbs_nday, nday_length + jbs_nday - 1
            m = n - jbs_nday + 1
            flxuc(m, k) = f_zero
            flxdc(m, k) = f_zero
            !$acc loop seq
            do ib = 1, nbdsw
               flxuc(m, k) = flxuc(m, k) + fxupc(m, k,ib)
               flxdc(m, k) = flxdc(m, k) + fxdnc(m, k,ib)
            enddo
         end do
      end do

         !! --- ...  optional clear sky fluxes

      if ( lhsw0 .or. lflxprf ) then
         !$acc parallel loop collapse(2) private(m) async(async_id)
         do k = 1, nlp1
            do n = jbs_nday, nday_length + jbs_nday - 1
               m = n - jbs_nday + 1
               flxu0(m, k) = f_zero
               flxd0(m, k) = f_zero
               !$acc loop seq
               do ib = 1, nbdsw
                  flxu0(m, k) = flxu0(m, k) + fxup0(m, k,ib)
                  flxd0(m, k) = flxd0(m, k) + fxdn0(m, k,ib)
               enddo
            end do
         end do
      endif

         !  --- ...  prepare for final outputs
      !$acc parallel loop collapse(2) private(m) async(async_id)
      do k = 1, nlay
         do n = jbs_nday, nday_length + jbs_nday - 1
            m = n - jbs_nday + 1
            rfdelp(m, k) = heatfac / delp(m, k)
         end do
      end do

      if ( lfdncmp ) then
         !$acc parallel loop private(j1, jf, jj) async(async_id)
         do n = jbs_nday, nday_length + jbs_nday - 1
            jj = map_nday_jj(n) - jjoffset
            jf = jjoffset+jj
            j1 = idxday(n)
            m = n - jbs_nday + 1
            !! --- ...  optional uv-b surface downward flux
            fdncmp_uvbf0(j1, jf) = suvbf0(m)
            fdncmp_uvbfc(j1, jf) = suvbfc(m)

            !! --- ...  optional beam and diffuse sfc fluxes
            fdncmp_nirbm(j1, jf) = sfbmc(m, 1)
            fdncmp_nirdf(j1, jf) = sfdfc(m, 1)
            fdncmp_visbm(j1, jf) = sfbmc(m, 2)
            fdncmp_visdf(j1, jf) = sfdfc(m, 2)
         end do
      endif    ! end if_lfdncmp

         !  --- ...  toa and sfc fluxes
      !$acc parallel loop private(j1, jf, jj, m) async(async_id)
      do n = jbs_nday, nday_length + jbs_nday - 1
         jj = map_nday_jj(n) - jjoffset
         jf = jjoffset+jj
         j1 = idxday(n)
         m = n - jbs_nday + 1
         topflx_upfxc(j1, jf) = ftoauc(m)
         topflx_dnfxc(j1, jf) = ftoadc(m)
         topflx_upfx0(j1, jf) = ftoau0(m)

         sfcflx_upfxc(j1, jf) = fsfcuc(m)
         sfcflx_dnfxc(j1, jf) = fsfcdc(m)
         sfcflx_upfx0(j1, jf) = fsfcu0(m)
         sfcflx_dnfx0(j1, jf) = fsfcd0(m)
      end do
      
      if (ivflip == 0) then       ! output from toa to sfc

               !  --- ...  compute heating rates
         !$acc parallel loop private(j1, kk, jf, jj, m) async(async_id)
         do n = jbs_nday, nday_length + jbs_nday - 1
            jj = map_nday_jj(n) - jjoffset
            jf = jjoffset+jj
            j1 = idxday(n)
            m = n - jbs_nday + 1
            fnet(m, 1) = flxdc(m, 1) - flxuc(m, 1)
            !$acc loop seq
            do k = 2, nlp1
               kk = nlp1 - k + 1
               fnet(m, k) = flxdc(m, k) - flxuc(m, k)
               hswc(j1, kk-kd, jf) = (fnet(m, k)-fnet(m, k-1)) * rfdelp(m, k-1)
            enddo
         end do

               !! --- ...  optional flux profiles

         if ( lflxprf ) then
            !$acc parallel loop collapse(2) private(kk, jf, jj, j1, m) async(async_id)
            do k = 1, nlp1
               do n = jbs_nday, nday_length + jbs_nday - 1
                  jj = map_nday_jj(n) - jjoffset
                  jf = jjoffset+jj
                  kk = nlp1 - k + 1
                  j1 = idxday(n)
                  m = n - jbs_nday + 1
                  flxprf_upfxc(j1,kk-kd, jf) = flxuc(m, k)
                  flxprf_dnfxc(j1,kk-kd, jf) = flxdc(m, k)
                  flxprf_upfx0(j1,kk-kd, jf) = flxu0(m, k)
                  flxprf_dnfx0(j1,kk-kd, jf) = flxd0(m, k)
               end do
            end do
         endif

               !! --- ...  optional clear sky heating rates

         if ( lhsw0 ) then
            !$acc parallel loop private(j1, kk, jf, jj, m) async(async_id)
            do n = jbs_nday, nday_length + jbs_nday - 1
               jj = map_nday_jj(n) - jjoffset
               jf = jjoffset+jj
               j1 = idxday(n)
               m = n - jbs_nday + 1
               fnet(m, 1) = flxd0(m, 1) - flxu0(m, 1)
               !$acc loop seq
               do k = 2, nlp1
                  kk = nlp1 - k + 1
                  fnet(m, k) = flxd0(m, k) - flxu0(m, k)
                  hsw0(j1,kk-kd, jf) = (fnet(m, k)-fnet(m, k-1)) * rfdelp(m, k-1)
               enddo
            end do
         endif

               !! --- ...  optional spectral band heating rates

         if ( lhswb ) then
            !$acc parallel loop collapse(2) private(jf, j1, kk, jj, m) async(async_id)
            do mb = 1, nbdsw
               do n = jbs_nday, nday_length + jbs_nday - 1
                  m = n - jbs_nday + 1
                  jj = map_nday_jj(n) - jjoffset
                  jf = jjoffset+jj
                  j1 = idxday(n)
                  fnetm(m, 1, mb) = fxdnc(m, 1,mb) - fxupc(m, 1,mb)
                  !$acc loop seq
                  do k = 2, nlp1
                     kk = nlp1 - k + 1
                     fnetm(m, k, mb) = fxdnc(m, k,mb) - fxupc(m, k,mb)
                     if (lhtrswb) hswb(j1,kk-kd,mb, jf) &
                        = (fnetm(m, k, mb) - fnetm(m, k-1, mb)) * rfdelp(m, k-1)
                  enddo
               end do
            end do
         endif

      else                        ! output from sfc to toa

               !  --- ...  compute heating rates
         !$acc parallel loop private(j1, jf, jj, m) async(async_id)
         do n = jbs_nday, nday_length + jbs_nday - 1
            jj = map_nday_jj(n) - jjoffset
            jf = jjoffset+jj
            j1 = idxday(n)
            m = n - jbs_nday + 1
            fnet(m, 1) = flxdc(m, 1) - flxuc(m, 1)
            !$acc loop seq
            do k = 2, nlp1
               fnet(m, k) = flxdc(m, k) - flxuc(m, k)
               hswc(j1,k-1-kd, jf) = (fnet(m, k)-fnet(m, k-1)) * rfdelp(m, k-1)
            enddo
         end do

               !! --- ...  optional flux profiles

         if ( lflxprf ) then
            !$acc parallel loop collapse(2) private(jf, j1, jj, m) async(async_id)
            do k = 1, nlp1
               do n = jbs_nday, nday_length + jbs_nday - 1
                  jj = map_nday_jj(n) - jjoffset
                  jf = jjoffset+jj
                  j1 = idxday(n)
                  m = n - jbs_nday + 1
                  flxprf_upfxc(j1,k-kd, jf) = flxuc(m, k)
                  flxprf_dnfxc(j1,k-kd, jf) = flxdc(m, k)
                  flxprf_upfx0(j1,k-kd, jf) = flxu0(m, k)
                  flxprf_dnfx0(j1,k-kd, jf) = flxd0(m, k)
               end do
            end do
         endif

               !! --- ...  optional clear sky heating rates

         if ( lhsw0 ) then
            !$acc parallel loop private(j1, jf, jj, m) async(async_id)
            do n = jbs_nday, nday_length + jbs_nday - 1
               jj = map_nday_jj(n) - jjoffset
               jf = jjoffset+jj
               m = n - jbs_nday + 1
               j1 = idxday(n)
               fnet(m, 1) = flxd0(m, 1) - flxu0(m, 1)
               !$acc loop seq
               do k = 2, nlp1
                  fnet(m, k) = flxd0(m, k) - flxu0(m, k)
                  hsw0(j1,k-1-kd, jf) = (fnet(m, k)-fnet(m, k-1)) * rfdelp(m, k-1)
               enddo
            end do
         endif

               !! --- ...  optional spectral band heating rates

         if ( lhswb ) then
            !$acc parallel loop collapse(2) private(jf, j1, jj, ipt, m) async(async_id)
            do mb = 1, nbdsw
               do n = jbs_nday, nday_length + jbs_nday - 1
                  jj = map_nday_jj(n) - jjoffset
                  ipt = map_nday_ipt(n)
                  jf = jjoffset+jj
                  j1 = idxday(n)
                  m = n - jbs_nday + 1
                  fnetm(m, 1, mb) = fxdnc(ipt, 1,mb) - fxupc(ipt, 1,mb)
                  !$acc loop seq
                  do k = 1, nlay
                     fnetm(m, k+1, mb) = fxdnc(ipt, k+1,mb) - fxupc(ipt, k+1,mb)
                     if (lhtrswb) hswb(j1,k-kd,mb, jf) &
                        = (fnetm(m, k+1, mb) - fnetm(m, k, mb)) * rfdelp(m, k)
                  enddo
               enddo
            end do
         endif

      endif                       ! if_ivflip
      end do
      !$acc end data


!      if (myrank .eq. 0) then
!          print *,'$$$ SW : flxprf(1,61)%upfxc = ', flxprf(1,61)%upfxc
!          print *,'$$$ SW : flxprf(1,61)%dnfxc = ', flxprf(1,61)%dnfxc
!          print *,'$$$ SW : flxprf(1,61)%upfx0 = ', flxprf(1,61)%upfx0
!          print *,'$$$ SW : flxprf(1,61)%dnfx0 = ', flxprf(1,61)%dnfx0
!          print *,'$$$ SW : flxprf(1,1)%upfxc = ', flxprf(1,1)%upfxc
!          print *,'$$$ SW : flxprf(1,1)%dnfxc = ', flxprf(1,1)%dnfxc
!          print *,'$$$ SW : flxprf(1,1)%upfx0 = ', flxprf(1,1)%upfx0
!          print *,'$$$ SW : flxprf(1,1)%dnfx0 = ', flxprf(1,1)%dnfx0
!      endif
      if (isubcsw > 0) deallocate(cldfmc)

      return
!...................................
      end subroutine swrad_gpu
!-----------------------------------


!-----------------------------------
      subroutine rswinit_gpu                                                &
!...................................

!  ---  inputs:
     &     ( me ,myrank)
!  ---  outputs: (none)

!  ===================  program usage description  ===================  !
!                                                                       !
! purpose:  initialize non-varying module variables, conversion factors,!
! and look-up tables.                                                   !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                              !
!    me       - print control for parallel process                      !
!                                                                       !
!  outputs: (none)                                                      !
!                                                                       !
!  external module variables:  (in physpara)                            !
!   iswrate - heating rate unit selections                              !
!           =1: output in k/day                                         !
!           =2: output in k/second                                      !
!   iswrgas - control flag for rare gases (ch4,n2o,o2, etc.)            !
!           =0: do not include rare gases                               !
!           >0: include all rare gases                                  !
!   iswcliq - liquid cloud optical properties contrl flag               !
!           =0: input cloud opt depth from diagnostic scheme            !
!           >0: input cwp,rew, and other cloud content parameters       !
!   isubcsw - sub-column cloud approximation control flag               !
!           =0: no sub-col cld treatment, use grid-mean cld quantities  !
!           =1: mcica sub-col, prescribed seeds to get random numbers   !
!           =2: mcica sub-col, providing array icseed for random numbers!
!   icldflg - cloud scheme control flag                                 !
!           =0: diagnostic scheme gives cloud tau, omiga, and g.        !
!           =1: prognostic scheme gives cloud liq/ice path, etc.        !
!   iovrsw  - clouds vertical overlapping control flag                  !
!           =0: random overlapping clouds                               !
!           =1: maximum/random overlapping clouds                       !
!           =2: maximum overlap cloud                                   !
!   iswmode - control flag for 2-stream transfer scheme                 !
!           =1; delta-eddington    (joseph et al., 1976)                !
!           =2: pifm               (zdunkowski et al., 1980)            !
!           =3: discrete ordinates (liou, 1973)                         !
!                                                                       !
!  *******************************************************************  !
!                                                                       !
! definitions:                                                          !
!     arrays for 10000-point look-up tables:                            !
!     tau_tbl  clear-sky optical depth                                  !
!     exp_tbl  exponential lookup table for transmittance               !
!                                                                       !
!  *******************************************************************  !
!                                                                       !
!  ======================  end of description block  =================  !

!  ---  inputs:
      integer, intent(in) :: me, myrank

!  ---  outputs: none

!  ---  locals:
      real (kind=kind_phys), parameter :: expeps = 1.e-20

      integer :: i

      real (kind=kind_phys) :: tfn, tau

!
!===> ... begin here
!
      if ( iovrsw<0 .or. iovrsw>2 ) then
        print *,'  *** error in specification of cloud overlap flag',   &
     &          ' iovrsw=',iovrsw,' in rswinit !!'
        stop
      endif

      if (me == 0 .and. myrank == 0 ) then
        print *,' - using aer shortwave radiation, version: ',vtagsw

        if (iswmode == 1) then
          print *,'   --- delta-eddington 2-stream transfer scheme'
        else if (iswmode == 2) then
          print *,'   --- pifm 2-stream transfer scheme'
        else if (iswmode == 3) then
          print *,'   --- discrete ordinates 2-stream transfer scheme'
        endif

        if (iswrgas <= 0) then
          print *,'   --- rare gases absorption is not included in sw'
        else
          print *,'   --- include rare gases n2o, ch4, o2, absorptions',&
     &            ' in sw'
        endif

        if ( isubcsw == 0 ) then
          print *,'   --- using standard grid average clouds, no ',     &
     &            'sub-column clouds approximation applied'
        elseif ( isubcsw == 1 ) then
          print *,'   --- using mcica sub-colum clouds approximation ', &
     &            'with a prescribed sequence of permutation seeds'
        elseif ( isubcsw == 2 ) then
          print *,'   --- using mcica sub-colum clouds approximation ', &
     &            'with provided input array of permutation seeds'
        else
          print *,'  *** error in specification of sub-column cloud ',  &
     &            ' control flag isubcsw =',isubcsw,' !!'
          stop
        endif
      endif

!  --- ...  check cloud flags for consistency

      if ((icldflg == 0 .and. iswcliq /= 0) .or.                        &
     &    (icldflg == 1 .and. iswcliq == 0)) then
        print *,'  *** model cloud scheme inconsistent with sw',        &
     &          ' radiation cloud radiative property setup !!'
        stop
      endif

!  --- ...  setup constant factors for heating rate
!           the 1.0e-2 is to convert pressure from mb to n/m**2

      if (iswrate == 1) then
!       heatfac = 8.4391
!       heatfac = con_g * 86400. * 1.0e-2 / con_cp  !   (in k/day)
        heatfac = con_g * 864.0 / con_cp            !   (in k/day)
      else
        heatfac = con_g * 1.0e-2 / con_cp           !   (in k/second)
      endif

!  --- ...  define exponential lookup tables for transmittance. tau is
!           computed as a function of the tau transition function, and
!           transmittance is calculated as a function of tau.  all tables
!           are computed at intervals of 0.0001.  the inverse of the
!           constant used in the pade approximation to the tau transition
!           function is set to bpade.

      exp_tbl(0) = 1.0
      exp_tbl(ntbmx) = expeps

      do i = 1, ntbmx-1
        tfn = float(i) / float(ntbmx-i)
        tau = bpade * tfn
        exp_tbl(i) = exp( -tau )
      enddo

      call ngptsw_packing_recipe_init(npacks, small_ngpts)

      return
!...................................
      end subroutine rswinit_gpu
!-----------------------------------


!-----------------------------------
      subroutine cldprop                                                &
!...................................
!  ---  inputs:
     &     ( cfrac,cliqp,reliq,cicep,reice,cdat1,cdat2,cdat3,cdat4,     &
     &       cf1, nlay, ipseed, ix, nday, idxday, map_nday_ipt, map_nday_jj, &
             nday_length, max_nday_length, jjoffset, jbs_nday, &
             jb, small_ngpts, packs_nbands, offset_ng, &
             async_id, fulljj, blockjj,                                           &
!  ---  output:
     &       taucw, ssacw, asycw, cldfrc, cldfmc                        &
     &     )

!  ===================  program usage description  ===================  !
!                                                                       !
! purpose: compute the cloud optical properties for each cloudy layer   !
! and g-point interval.                                                 !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                        size  !
!    cfrac - real, layer cloud fraction                            nlay !
!        .....  for  iswcliq > 0 (prognostic cloud sckeme)  - - -       !
!    cliqp - real, layer in-cloud liq water path (g/m**2)          nlay !
!    reliq - real, mean eff radius for liq cloud (micron)          nlay !
!    cicep - real, layer in-cloud ice water path (g/m**2)          nlay !
!    reice - real, mean eff radius for ice cloud (micron)          nlay !
!    cdat1 - real, layer rain drop water path (g/m**2)             nlay !
!    cdat2 - real, effective radius for rain drop (micron)         nlay !
!    cdat3 - real, layer snow flake water path(g/m**2)             nlay !
!    cdat4 - real, mean eff radius for snow flake(micron)          nlay !
!        .....  for iswcliq = 0  (diagnostic cloud sckeme)  - - -       !
!    cdat1 - real, layer cloud optical depth                       nlay !
!    cdat2 - real, layer cloud single scattering albedo            nlay !
!    cdat3 - real, layer cloud asymmetry factor                    nlay !
!    cdat4 - real, optional use                                    nlay !
!    cliqp - real, not used                                        nlay !
!    cicep - real, not used                                        nlay !
!    reliq - real, not used                                        nlay !
!    reice - real, not used                                        nlay !
!                                                                       !
!    cf1   - real, effective total cloud cover at surface           1   !
!    nlay  - integer, vertical layer number                         1   !
!    ipseed- permutation seed for generating random numbers (isubcsw>0) !
!                                                                       !
!  outputs:                                                             !
!    taucw  - real, cloud optical depth, w/o delta scaled    nlay*nbdsw !
!    ssacw  - real, weighted cloud single scattering albedo  nlay*nbdsw !
!                             (ssa = ssacw / taucw)                     !
!    asycw  - real, weighted cloud asymmetry factor          nlay*nbdsw !
!                             (asy = asycw / ssacw)                     !
!    cldfrc - real, cloud fraction of grid mean value              nlay !
!    cldfmc - real, cloud fraction for each sub-column       nlay*ngptsw!
!                                                                       !
!                                                                       !
!  explanation of the method for each value of iswcliq, and iswcice.    !
!  set up in module "physpara"                                          !
!                                                                       !
!     iswcliq=0  : input cloud optical property (tau, ssa, asy).        !
!                  (used for diagnostic cloud method)                   !
!     iswcliq>0  : input cloud liq/ice path and effective radius, also  !
!                  require the user of 'iswcice' to specify the method  !
!                  used to compute aborption due to water/ice parts.    !
!  ...................................................................  !
!                                                                       !
!     iswcliq=1  : liquid water cloud optical properties are computed   !
!                  as in hu and stamnes (1993), j. clim., 6, 728-742.   !
!                                                                       !
!     iswcice used only when iswcliq > 0                                !
!                  the cloud ice path (g/m2) and ice effective radius   !
!                  (microns) are inputs.                                !
!     iswcice=1  : ice cloud optical properties are computed as in      !
!                  ebert and curry (1992), jgr, 97, 3831-3836.          !
!     iswcice=2  : ice cloud optical properties are computed as in      !
!                  streamer v3.0 (2001), key, streamer user's guide,    !
!                  cooperative institude for meteorological studies,95pp!
!     iswcice=3  : ice cloud optical properties are computed as in      !
!                  fu (1996), j. clim., 9.                              !
!                                                                       !
!  other cloud control module variables:                                !
!     isubcsw =0: standard cloud scheme, no sub-col cloud approximation !
!             >0: mcica sub-col cloud scheme using ipseed as permutation!
!                 seed for generating rundom numbers                    !
!                                                                       !
!  ======================  end of description block  =================  !
!
      use module_radsw_cldprtb

!  ---  inputs:
      integer, intent(in) :: nlay, fulljj, ipseed(ix, fulljj), ix, jbs_nday, &
         nday(fulljj), idxday(nday_length), blockjj, nday_length, jjoffset, &
         max_nday_length, jb(nbdsw), small_ngpts, packs_nbands, offset_ng
      real (kind=kind_phys), intent(in) :: cf1(max_nday_length)
      integer, dimension(nday_length) :: map_nday_ipt, map_nday_jj

      real (kind=kind_phys), dimension(nxptot, nlay), intent(in) :: cfrac, &
         cliqp, reliq, cicep, reice, cdat1, cdat2, cdat3, cdat4

!  ---  outputs:
      real (kind=kind_phys), dimension(ix, nlay,ngptsw, fulljj), intent(out) ::     &
     &       cldfmc
      real (kind=kind_phys), dimension(nday_length, nlay, max_packs_nbands),  intent(out) ::     &
     &       taucw, ssacw, asycw
      real (kind=kind_phys), dimension(max_nday_length, nlay), intent(out) :: cldfrc

!  ---  locals:
      real (kind=kind_phys), dimension(nblow:nbhgh) :: tauliq, tauice,  &
     &       ssaliq, ssaice, ssaran, ssasnw, asyliq, asyice,            &
     &       asyran, asysnw
     
      real (kind=kind_phys) :: dgeice, factor, fint, tauran, tausnw,    &
     &       cldliq, refliq, cldice, refice, cldran, cldsnw, refsnw,    &
     &       extcoliq, ssacoliq, asycoliq, extcoice, ssacoice, asycoice,&
     &       dgesnw
     
      real (kind=kind_phys), dimension(:,:,:), allocatable       :: cldf
      logical, allocatable, dimension(:,:,:,:) :: lcloudy
      integer :: ia, ib, ig, k, index, ipt, jj, j1, async_id, n, rr, nn
      real (kind=kind_phys), dimension(:), allocatable       :: cldf_im
      logical, allocatable, dimension(:,:) :: lcloudy_im

      !
      !===> ...  begin here
      !
      !$acc parallel loop collapse(3) async(async_id) 
      do nn = 1, packs_nbands
         do k = 1, nlay
            do n = 1, nday_length
               taucw (n, k,nn) = f_zero
               asycw (n, k,nn) = f_zero
               if (cf1(n) > f_zero) then
                  ssacw (n, k,nn) = f_one
               else
                  ssacw (n, k,nn) = f_zero
               end if
            enddo
         enddo
      end do


               !  --- ...  compute cloud radiative properties for a cloudy column
      if (iswcliq > 0) then ! lab_if_iswcliq
         !$acc parallel loop collapse(3) private(cldran, cldsnw, refsnw, dgesnw, tauran, tausnw, &
         !$acc&     ssaran, ssasnw, asyran, asysnw, cldliq, cldice, refliq, refice, &
         !$acc&     tauliq, ssaliq, asyliq, factor, index, fint, extcoliq, ssacoliq, &
         !$acc&     asycoliq, tauice, ssaice, asyice, ia, extcoice, ssacoice, asycoice, &
         !$acc&     dgeice, j1, rr)async(async_id)
         do nn = 1, packs_nbands
            do k = 1, nlay ! lab_do_k
               do n = 1, nday_length
                  j1 = idxday(n)
                  ib = jb(nn)
                  rr = j1 + nxjp_acc(map_nday_jj(n)) - 1
                  if (cf1(n) > f_zero) then     ! cloudy sky column
                     if (cfrac(rr, k) > ftiny) then ! lab_if_cld

                        !  --- ...  optical properties for rain and snow
                        cldran = cdat1(rr, k)
                        !           refran = cdat2(k)
                        cldsnw = cdat3(rr, k)
                        refsnw = cdat4(rr, k)
                        dgesnw = 1.0315 * refsnw        ! for fu's snow formula

                        tauran = cldran * a0r
                        !  ---  if use fu's formula it needs to be normalized by snow/ice density
                        !       !not use snow density = 0.1 g/cm**3 = 0.1 g/(mu * m**2)
                        !       use ice density = 0.9167 g/cm**3 = 0.9167 g/(mu * m**2)
                        !       1/0.9167 = 1.09087
                        !       factor 1.5396=8/(3*sqrt(3)) converts reff to generalized ice particle size
                        !       use newer factor value 1.0315
                        if (cldsnw>f_zero .and. refsnw>10.0_kind_phys) then
                           !             tausnw = cldsnw * (a0s + a1s/refsnw)
                           tausnw = cldsnw*1.09087*(a0s + a1s/dgesnw)     ! fu's formula
                        else
                           tausnw = f_zero
                        endif
                        ssaran(ib) = tauran * (f_one - b0r(ib))
                        ssasnw(ib) = tausnw * (f_one - (b0s(ib)+b1s(ib)*dgesnw))
                        asyran(ib) = ssaran(ib) * c0r(ib)
                        asysnw(ib) = ssasnw(ib) * c0s(ib)

                        cldliq = cliqp(rr, k)
                        cldice = cicep(rr, k)
                        refliq = reliq(rr, k)
                        refice = reice(rr, k)

                        !  --- ...  calculation of absorption coefficients due to water clouds.

                        if ( cldliq <= f_zero ) then
                           tauliq(ib) = f_zero
                           ssaliq(ib) = f_zero
                           asyliq(ib) = f_zero
                        else
                           if ( iswcliq == 1 ) then
                              factor = refliq - 1.5
                              index  = max( 1, min( 57, int( factor ) ))
                              fint   = factor - float(index)
                              extcoliq = max(f_zero,            extliq1(index,ib)   &
                              &              + fint*(extliq1(index+1,ib)-extliq1(index,ib)) )
                              ssacoliq = max(f_zero, min(f_one, ssaliq1(index,ib)   &
                              &              + fint*(ssaliq1(index+1,ib)-ssaliq1(index,ib)) ))

                              asycoliq = max(f_zero, min(f_one, asyliq1(index,ib)   &
                              &              + fint*(asyliq1(index+1,ib)-asyliq1(index,ib)) ))
                              !                 forcoliq = asycoliq * asycoliq

                              tauliq(ib) = cldliq     * extcoliq
                              ssaliq(ib) = tauliq(ib) * ssacoliq
                              asyliq(ib) = ssaliq(ib) * asycoliq
                           endif   ! end if_iswcliq_block
                        endif   ! end if_cldliq_block

                        !  --- ...  calculation of absorption coefficients due to ice clouds.

                        if ( cldice <= f_zero ) then
                           tauice(ib) = f_zero
                           ssaice(ib) = f_zero
                           asyice(ib) = f_zero
                        else

                           !  --- ...  ebert and curry approach for all particle sizes though somewhat
                           !           unjustified for large ice particles

                           if ( iswcice == 1 ) then
                              refice = min(130.0_kind_phys,max(13.0_kind_phys,refice))
                                 ia = idxebc(ib)           ! eb_&_c band index for ice cloud coeff

                                 extcoice = max(f_zero, abari(ia)+bbari(ia)/refice )
                                 ssacoice = max(f_zero, min(f_one,                     &
                                 &                             f_one-cbari(ia)-dbari(ia)*refice ))
                                 asycoice = max(f_zero, min(f_one,                     &
                                 &                                   ebari(ia)+fbari(ia)*refice ))
                                 !                 forcoice = asycoice * asycoice

                                 tauice(ib) = cldice     * extcoice
                                 ssaice(ib) = tauice(ib) * ssacoice
                                 asyice(ib) = ssaice(ib) * asycoice

                              !  --- ...  streamer approach for ice effective radius between 5.0 and 131.0 microns

                           elseif ( iswcice == 2 ) then
                              refice = min(131.0_kind_phys,max(5.0_kind_phys,refice))

                              factor = (refice - 2.0) / 3.0
                              index  = max( 1, min( 42, int( factor ) ))
                              fint   = factor - float(index)
                              extcoice = max(f_zero,            extice2(index,ib)   &
                              &                + fint*(extice2(index+1,ib)-extice2(index,ib)) )
                              ssacoice = max(f_zero, min(f_one, ssaice2(index,ib)   &
                              &                + fint*(ssaice2(index+1,ib)-ssaice2(index,ib)) ))
                              asycoice = max(f_zero, min(f_one, asyice2(index,ib)   &
                              &                + fint*(asyice2(index+1,ib)-asyice2(index,ib)) ))
                              !                 forcoice = asycoice * asycoice

                              tauice(ib) = cldice     * extcoice
                              ssaice(ib) = tauice(ib) * ssacoice
                              asyice(ib) = ssaice(ib) * asycoice

                              !  --- ...  fu's approach for ice effective radius between 4.8 and 135 microns
                              !           (generalized effective size from 5 to 140 microns)

                           elseif ( iswcice == 3 ) then
                              dgeice = max( 5.0, min( 140.0, 1.0315*refice ))

                              factor = (dgeice - 2.0) / 3.0
                              index  = max( 1, min( 45, int( factor ) ))
                              fint   = factor - float(index)
                              extcoice = max(f_zero,            extice3(index,ib)   &
                              &                + fint*(extice3(index+1,ib)-extice3(index,ib)) )
                              ssacoice = max(f_zero, min(f_one, ssaice3(index,ib)   &
                              &                + fint*(ssaice3(index+1,ib)-ssaice3(index,ib)) ))
                              asycoice = max(f_zero, min(f_one, asyice3(index,ib)   &
                              &                + fint*(asyice3(index+1,ib)-asyice3(index,ib)) ))
                              !                 fdelta   = max(f_zero, min(f_one, fdlice3(index,ib)   &
                              !    &                + fint*(fdlice3(index+1,ib)-fdlice3(index,ib)) ))
                              !                 forcoice = min( asycoice, fdelta+0.5/ssacoice )           ! see fu 1996 p. 2067

                              tauice(ib) = cldice     * extcoice
                              ssaice(ib) = tauice(ib) * ssacoice
                              asyice(ib) = ssaice(ib) * asycoice

                           endif   ! end if_iswcice_block
                        endif   ! end if_cldice_block
                        taucw(n, k,nn) = tauliq(ib)+tauice(ib)+tauran+tausnw
                        ssacw(n, k,nn) = ssaliq(ib)+ssaice(ib)+ssaran(ib)+ssasnw(ib)
                        asycw(n, k,nn) = asyliq(ib)+asyice(ib)+asyran(ib)+asysnw(ib)

                     endif  ! lab_if_cld
                  end if
               enddo  ! lab_do_k
            end do
         end do

      else  ! lab_if_iswcliq
         !$acc parallel loop collapse(3) private(j1, rr) async(async_id)
         do nn = 1, packs_nbands
            do k = 1, nlay
               do n = 1, nday_length
                  j1 = idxday(n)
                  rr = j1 + nxjp_acc(map_nday_jj(n)) - 1
                  if (cf1(n) > f_zero) then     ! cloudy sky column
                     if (cfrac(rr, k) > ftiny) then
                        taucw(n, k,nn) = cdat1(rr, k)
                        ssacw(n, k,nn) = cdat1(rr, k)    * cdat2(rr, k)
                        asycw(n, k,nn) = ssacw(n, k,nn) * cdat3(rr, k)
                     endif
                  end if
               enddo
            end do
         end do

      endif  ! lab_if_iswcliq

               !  ---  distribute cloud properties to each g-point

      if ( isubcsw > 0 ) then      ! mcica sub-col clouds approx
         allocate(cldf(ix, nlay, fulljj))
         allocate(cldf_im(nlay))
         allocate(lcloudy_im(nlay,ngptsw))
         allocate(lcloudy(ix, nlay,ngptsw, fulljj))
         do n = 1, nday_length
            jj = map_nday_jj(n) - jjoffset
            ipt = map_nday_ipt(n)
            j1 = idxday(n)
            if (cf1(n) > f_zero) then     ! cloudy sky column
               cldf(ipt, :, jj) = cfrac(n, :)
               where (cldf(ipt, :, jj) < ftiny)
                  cldf(ipt, :, jj) = f_zero
               end where

               !  --- ...  call sub-column cloud generator
               do k = 1, nlay
                  cldf_im(k) = cldf(ipt, k, jj)
                  do ib = 1, ngptsw
                     lcloudy_im(k,ig) = lcloudy(ipt, k,ig, jj)
                  end do
               end do
               
               call mcica_subcol                                               &
               !  ---  inputs:
               &     ( cldf_im, nlay, ipseed(j1, jj),                                        &
               !  ---  outputs:
               &       lcloudy_im                                                    &
               &     )
               do k = 1, nlay
                  cldf(ipt, k, jj) = cldf_im(k)
                  do ib = 1, ngptsw
                     lcloudy(ipt, k,ig, jj) = lcloudy_im(k,ig)
                  end do
               end do

               do ig = 1, ngptsw
                  do k = 1, nlay
                     if ( lcloudy(ipt, k,ig, jj) ) then
                        cldfmc(ipt, k,ig, jj) = f_one
                     else
                        cldfmc(ipt, k,ig, jj) = f_zero
                     endif
                  enddo
               enddo
            end if
         end do
         deallocate(cldf)
         deallocate(cldf_im)
         deallocate(lcloudy_im)
         deallocate(lcloudy)

      else                         ! non-mcica, normalize cloud
         !$acc parallel loop collapse(2) private(j1, jj, ipt, rr) async(async_id)
         do k = 1, nlay
            do n = 1, nday_length
               jj = map_nday_jj(n) - jjoffset
               ipt = map_nday_ipt(n)
               j1 = idxday(n)
               rr = j1 + nxjp_acc(map_nday_jj(n)) - 1
               if (cf1(n) > f_zero) then     ! cloudy sky column
                  cldfrc(n, k) = cfrac(rr, k) / cf1(n)
               else
                  cldfrc(n, k) = f_zero
               end if
            enddo
         end do
      endif   ! end if_isubcsw_block

      return
!...................................
      end subroutine cldprop
!-----------------------------------


! ----------------------------------
      subroutine mcica_subcol                                           &
! ..................................
!  ---  inputs:
     &    ( cldf, nlay, ipseed,                                         &
!  ---  outputs:
     &      lcloudy                                                     &
     &    )

!  ====================  defination of variables  ====================  !
!                                                                       !
!  input variables:                                                size !
!   cldf    - real, layer cloud fraction                           nlay !
!   nlay    - integer, number of model vertical layers               1  !
!   ipseed  - integer, permute seed for random num generator         1  !
!    ** note : if the cloud generator is called multiple times, need    !
!              to permute the seed between each call; if between calls  !
!              for lw and sw, use values differ by the number of g-pts. !
!                                                                       !
!  output variables:                                                    !
!   lcloudy - logical, sub-colum cloud profile flag array    nlay*ngptsw!
!                                                                       !
!  other control flags from module variables:                           !
!     iovrsw    : control flag for cloud overlapping method             !
!                 =0:random; =1:maximum/random; =2:maximum              !
!                                                                       !
!                                                                       !
!  =====================    end of definitions    ====================  !

      implicit none

!  ---  inputs:
      integer, intent(in) :: nlay, ipseed

      real (kind=kind_phys), dimension(nlay), intent(in) :: cldf

!  ---  outputs:
      logical, dimension(nlay,ngptsw), intent(out):: lcloudy

!  ---  locals:
      real (kind=kind_phys) :: cdfunc(nlay,ngptsw), tem1,               &
     &       rand2d(nlay*ngptsw), rand1d(ngptsw)

      type (random_stat) :: stat          ! for thread safe random generator

      integer :: k, n, k1
!
!===> ...  begin here
!
!  --- ...  advance randum number generator by ipseed values

      call random_setseed                                               &
!  ---  inputs:
     &    ( ipseed,                                                     &
!  ---  outputs:
     &      stat                                                        &
     &    )

!  --- ...  sub-column set up according to overlapping assumption

      select case ( iovrsw )

        case( 0 )        ! random overlap, pick a random value at every level

          call random_number                                            &
!  ---  inputs: ( none )
!  ---  outputs:
     &     ( rand2d, stat )

          k1 = 0
          do n = 1, ngptsw
            do k = 1, nlay
              k1 = k1 + 1
              cdfunc(k,n) = rand2d(k1)
            enddo
          enddo

        case( 1 )        ! max-ran overlap

          call random_number                                            &
!  ---  inputs: ( none )
!  ---  outputs:
     &     ( rand2d, stat )

          k1 = 0
          do n = 1, ngptsw
            do k = 1, nlay
              k1 = k1 + 1
              cdfunc(k,n) = rand2d(k1)
            enddo
          enddo

!  ---  first pick a random number for bottom/top layer.
!       then walk up the column: (aer's code)
!       if layer below is cloudy, use the same rand num in the layer below
!       if layer below is clear,  use a new random number

!  ---  from bottom up
          do k = 2, nlay
            k1 = k - 1
            tem1 = f_one - cldf(k1)

            do n = 1, ngptsw
              if ( cdfunc(k1,n) > tem1 ) then
                cdfunc(k,n) = cdfunc(k1,n)
              else
                cdfunc(k,n) = cdfunc(k,n) * tem1
              endif
            enddo
          enddo

!  ---  then walk down the column: (if use original author's method)
!       if layer above is cloudy, use the same rand num in the layer above
!       if layer above is clear,  use a new random number

!  ---  from top down
!         do k = nlay-1, 1, -1
!           k1 = k + 1
!           tem1 = f_one - cldf(k1)

!           do n = 1, ngptsw
!             if ( cdfunc(k1,n) > tem1 ) then
!               cdfunc(k,n) = cdfunc(k1,n)
!             else
!               cdfunc(k,n) = cdfunc(k,n) * tem1
!             endif
!           enddo
!         enddo

        case( 2 )        ! maximum overlap, pick same random numebr at every level

          call random_number                                            &
!  ---  inputs: ( none )
!  ---  outputs:
     &     ( rand1d, stat )

          do n = 1, ngptsw
            tem1 = rand1d(n)

            do k = 1, nlay
              cdfunc(k,n) = tem1
            enddo
          enddo

      end select

!  --- ...  generate subcolumns for homogeneous clouds

      do k = 1, nlay
        tem1 = f_one - cldf(k)

        do n = 1, ngptsw
          lcloudy(k,n) = cdfunc(k,n) >= tem1
        enddo
      enddo

      return
! ..................................
      end subroutine mcica_subcol
! ----------------------------------


! ----------------------------------
      subroutine setcoef                                                &
! ..................................
!  ---  inputs:
     &     ( pavel,tavel,h2ovmr, nlay,nlp1, ix, nday, &
             async_id, fulljj, blockjj,                             &
             map_nday_ipt, map_nday_jj, nday_length, max_nday_length, jjoffset, &
!  ---  outputs:
     &       laytrop,jp,jt,jt1,fac00,fac01,fac10,fac11,                 &
     &       selffac,selffrac,indself,forfac,forfrac,indfor             &
     &     )


!  ---  inputs:
      integer, intent(in) :: nlay, nlp1, ix, nday(fulljj), fulljj, blockjj, &
         nday_length, jjoffset, max_nday_length
      integer, dimension(nday_length) :: map_nday_ipt, map_nday_jj

      real (kind=kind_phys), dimension(:,:), intent(in) :: pavel, tavel,  &
     &       h2ovmr

!  ---  outputs:
      integer, dimension(max_nday_length, nlay), intent(out) :: indself, indfor,         &
     &       jp, jt, jt1
      integer, intent(out) :: laytrop(nday_length)

      real (kind=kind_phys), dimension(:, :), intent(out) :: fac00,     &
     &       fac01, fac10, fac11, selffac, selffrac, forfac, forfrac

!  ---  locals:
      real (kind=kind_phys) :: plog, fp, fp1, ft, ft1, tem1, tem2

      integer :: i, k, jp1, jj, ipt, async_id, n
!
!===> ... begin here
!
      !$acc parallel loop async(async_id)
      do n = 1, nday_length
         laytrop(n) = nlay
         !$acc loop seq
         do k = 1, nlay
            if (log(pavel(n, k)) > 4.56) laytrop(n) =  k
         end do
      end do
      !$acc parallel loop collapse(2) private(plog, jp1, fp, tem1, tem2, &
      !$acc&         ft, ft1, fp1) async(async_id)
      do k = 1, nlay
         do n = 1, nday_length
            forfac(n, k) = pavel(n, k)*stpfac / (tavel(n, k)*(f_one + h2ovmr(n, k)))

            !  --- ...  find the two reference pressures on either side of the
            !           layer pressure.  store them in jp and jp1.  store in fp the
            !           fraction of the difference (in ln(pressure)) between these
            !           two values that the layer pressure lies.

            plog  = log(pavel(n, k))
            jp(n, k) = max(1, min(58, int(36.0 - 5.0*(plog+0.04)) ))
            jp1   = jp(n, k) + 1
            fp    = 5.0 * (preflog(jp(n, k)) - plog)

            !  --- ...  determine, for each reference pressure (jp and jp1), which
            !          reference temperature (these are different for each reference
            !          pressure) is nearest the layer temperature but does not exceed it.
            !          store these indices in jt and jt1, resp. store in ft (resp. ft1)
            !          the fraction of the way between jt (jt1) and the next highest
            !          reference temperature that the layer temperature falls.

            tem1 = (tavel(n, k) - tref(jp(n, k))) / 15.0
            tem2 = (tavel(n, k) - tref(jp1  )) / 15.0
            jt (n, k) = max(1, min(4, int(3.0 + tem1) ))
            jt1(n, k) = max(1, min(4, int(3.0 + tem2) ))
            ft  = tem1 - float(jt (n, k) - 3)
            ft1 = tem2 - float(jt1(n, k) - 3)

            !  --- ...  we have now isolated the layer ln pressure and temperature,
            !           between two reference pressures and two reference temperatures
            !           (for each reference pressure).  we multiply the pressure
            !           fraction fp with the appropriate temperature fractions to get
            !           the factors that will be needed for the interpolation that yields
            !           the optical depths (performed in routines taugbn for band n).

            fp1 = f_one - fp
            fac10(n, k) = fp1 * ft
            fac00(n, k) = fp1 * (f_one - ft)
            fac11(n, k) = fp  * ft1
            fac01(n, k) = fp  * (f_one - ft1)
         end do
      end do
      !$acc parallel loop collapse(2) private(tem1, tem2) async(async_id)
      do k = 1, nlay
         do n = 1, nday_length
            !  --- ...  if the pressure is less than ~100mb, perform a different
            !           set of species interpolations.

            if ( log(pavel(n, k)) > 4.56 ) then

               !  --- ...  set up factors needed to separately include the water vapor
               !           foreign-continuum in the calculation of absorption coefficient.

               tem1 = (332.0 - tavel(n, k)) / 36.0
               indfor (n, k) = min(2, max(1, int(tem1)))
               forfrac(n, k) = tem1 - float(indfor(n, k))

               !  --- ...  set up factors needed to separately include the water vapor
               !           self-continuum in the calculation of absorption coefficient.

               tem2 = (tavel(n, k) - 188.0) / 7.2
               indself (n, k) = min(9, max(1, int(tem2)-7))
               selffrac(n, k) = tem2 - float(indself(n, k) + 7)
               selffac (n, k) = h2ovmr(n, k) * forfac(n, k)

            else

               !  --- ...  set up factors needed to separately include the water vapor
               !           foreign-continuum in the calculation of absorption coefficient.

               tem1 = (tavel(n, k) - 188.0) / 36.0
               indfor (n, k) = 3
               forfrac(n, k) = tem1 - f_one

               indself (n, k) = 0
               selffrac(n, k) = f_zero
               selffac (n, k) = f_zero

            endif

         enddo    ! end_do_k_loop
      end do

      return
! ..................................
      end subroutine setcoef
! ----------------------------------


!-----------------------------------
      subroutine spcvrtc                                                &
!...................................
!  ---  inputs:
           ( plvl, plyr, prslk1, tvly, rhly, slmsk, tracer1, xlon, xlat, &
             myim, lsswr, lslwr, me, map_jj, map_i, ntrac, my_max, &
      &      cfrac, cliqp, reliq, cicep, reice, cdat1, cdat2, cdat3, cdat4, ipseed, &
             colamt,colmol,pavel, tavel, h2ovmr, &
     &       ssolar,cosz,sntz,albbm,albdf,              &
     &       cf1,cf0,      &
     &       nlay, nlp1, ix, idxday, nday, map_nday_ipt, map_nday_jj, &
             nday_length, max_nday_length, jjoffset, jbs_nday, &
             async_id, fulljj, blockjj,                                                &
!  ---  outputs:
     &       fxupc,fxdnc,fxup0,fxdn0,                                   &
     &       ftoauc,ftoau0,ftoadc,fsfcuc,fsfcu0,fsfcdc,fsfcd0,          &
     &       sfbmc,sfdfc,sfbm0,sfdf0,suvbfc,suvbf0                      &
     &     )

!  ===================  program usage description  ===================  !
!                                                                       !
!   purpose:  computes the shortwave radiative fluxes using two-stream  !
!             method                                                    !
!                                                                       !
!   subprograms called:  swflux                                         !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                        size  !
!    ssolar  - real, incoming solar flux at top                    1    !
!    cosz    - real, cosine solar zenith angle                     1    !
!    sntz    - real, secant solar zenith angle                     1    !
!    albbm   - real, surface albedo for direct beam radiation      2    !
!    albdf   - real, surface albedo for diffused radiation         2    !
!    sfluxzen- real, spectral distribution of incoming solar flux ngptsw!
!    cldfrc  - real, layer cloud fraction                         nlay  !
!    cf1     - real, >0: cloudy sky, otherwise: clear sky          1    !
!    cf0     - real, =1-cf1                                        1    !
!    taug    - real, spectral optical depth for gases        nlay*ngptsw!
!    taur    - real, optical depth for rayleigh scattering   nlay*ngptsw!
!    tauae   - real, aerosols optical depth                  nlay*nbdsw !
!    ssaae   - real, aerosols single scattering albedo       nlay*nbdsw !
!    asyae   - real, aerosols asymmetry factor               nlay*nbdsw !
!    taucw   - real, weighted cloud optical depth            nlay*nbdsw !
!    ssacw   - real, weighted cloud single scat albedo       nlay*nbdsw !
!    asycw   - real, weighted cloud asymmetry factor         nlay*nbdsw !
!    nlay,nlp1 - integer,  number of layers/levels                 1    !
!                                                                       !
!  output variables:                                                    !
!    fxupc   - real, tot sky upward flux                     nlp1*nbdsw !
!    fxdnc   - real, tot sky downward flux                   nlp1*nbdsw !
!    fxup0   - real, clr sky upward flux                     nlp1*nbdsw !
!    fxdn0   - real, clr sky downward flux                   nlp1*nbdsw !
!    ftoauc  - real, tot sky toa upwd flux                         1    !
!    ftoau0  - real, clr sky toa upwd flux                         1    !
!    ftoadc  - real, toa downward (incoming) solar flux            1    !
!    fsfcuc  - real, tot sky sfc upwd flux                         1    !
!    fsfcu0  - real, clr sky sfc upwd flux                         1    !
!    fsfcdc  - real, tot sky sfc dnwd flux                         1    !
!    fsfcd0  - real, clr sky sfc dnwd flux                         1    !
!    sfbmc   - real, tot sky sfc dnwd beam flux (nir/uv+vis)       2    !
!    sfdfc   - real, tot sky sfc dnwd diff flux (nir/uv+vis)       2    !
!    sfbm0   - real, clr sky sfc dnwd beam flux (nir/uv+vis)       2    !
!    sfdf0   - real, clr sky sfc dnwd diff flux (nir/uv+vis)       2    !
!    suvbfc  - real, tot sky sfc dnwd uv-b flux                    1    !
!    suvbf0  - real, clr sky sfc dnwd uv-b flux                    1    !
!                                                                       !
!  internal variables:                                                  !
!    zrefb   - real, direct beam reflectivity for clear/cloudy    nlp1  !
!    zrefd   - real, diffuse reflectivity for clear/cloudy        nlp1  !
!    ztrab   - real, direct beam transmissivity for clear/cloudy  nlp1  !
!    ztrad   - real, diffuse transmissivity for clear/cloudy      nlp1  !
!    zldbt   - real, layer beam transmittance for clear/cloudy    nlp1  !
!    ztdbt   - real, lev total beam transmittance for clr/cld     nlp1  !
!                                                                       !
!  control parameters in module "physpara"                              !
!    iswmode - control flag for 2-stream transfer schemes               !
!              = 1 delta-eddington    (joseph et al., 1976)             !
!              = 2 pifm               (zdunkowski et al., 1980)         !
!              = 3 discrete ordinates (liou, 1973)                      !
!                                                                       !
!  *******************************************************************  !
!  original code description                                            !
!                                                                       !
!  method:                                                              !
!  -------                                                              !
!     standard delta-eddington, p.i.f.m., or d.o.m. layer calculations. !
!     kmodts  = 1 eddington (joseph et al., 1976)                       !
!             = 2 pifm (zdunkowski et al., 1980)                        !
!             = 3 discrete ordinates (liou, 1973)                       !
!                                                                       !
!  modifications:                                                       !
!  --------------                                                       !
!   original: h. barker                                                 !
!   revision: merge with rrtmg_sw: j.-j.morcrette, ecmwf, feb 2003      !
!   revision: add adjustment for earth/sun distance:mjiacono,aer,oct2003!
!   revision: bug fix for use of palbp and palbd: mjiacono, aer, nov2003!
!   revision: bug fix to apply delta scaling to clear sky: aer, dec2004 !
!   revision: code modified so that delta scaling is not done in cloudy !
!             profiles if routine cldprop is used; delta scaling can be !
!             applied by swithcing code below if cldprop is not used to !
!             get cloud properties. aer, jan 2005                       !
!   revision: uniform formatting for rrtmg: mjiacono, aer, jul 2006     !
!   revision: use exponential lookup table for transmittance: mjiacono, !
!             aer, aug 2007                                             !
!                                                                       !
!  *******************************************************************  !
!  ======================  end of description block  =================  !

!  ---  constant parameters:
      real (kind=kind_phys), parameter :: zcrit = 0.9999995 ! thresold for conservative scattering
      real (kind=kind_phys), parameter :: zsr3  = sqrt(3.0)
      real (kind=kind_phys), parameter :: od_lo = 0.06
      real (kind=kind_phys), parameter :: eps1  = 1.0e-8
      integer, parameter :: small_factor = 7
      integer, parameter :: small_ngptsw = int(ngptsw/small_factor)

!  ---  inputs:
      real (kind=kind_phys), dimension(nxptot,nlp1), intent(in) ::        &
     &       plvl
      real (kind=kind_phys), dimension(nxptot,nlay), intent(in) ::        &
     &       plyr
      integer,  intent(in) :: me, my_max, ntrac
      logical,  intent(in) :: lsswr, lslwr
      real (kind=kind_phys), dimension(ix, my_max), intent(in) ::  slmsk,  &
               xlon, xlat
      real (kind=kind_phys), dimension(nxptot,nlay)  :: rhly, prslk1, tvly
      real (kind=kind_phys), dimension(nxptot,nlay,ntrac)   :: tracer1
      integer, dimension(nxptot) :: map_jj, map_i
      integer, intent(in) :: myim(fulljj)
      real (kind=kind_phys), dimension(nxptot,nlay),intent(in):: &
         cliqp, reliq, cicep, reice, cdat4
      real (kind=kind_phys), dimension(nxptot,nlay),intent(inout):: cfrac
      real (kind=kind_phys), dimension(nxptot,nlay),intent(inout):: cdat1, cdat2, cdat3
      integer, dimension(ix, fulljj) :: ipseed
      integer, intent(in) :: nday_length, jjoffset, max_nday_length, jbs_nday
      integer, dimension(nday_length) :: map_nday_ipt, map_nday_jj
      real (kind=kind_phys), dimension(max_nday_length, nlay), intent(in) :: colmol
      real (kind=kind_phys), dimension(max_nday_length, nlay), intent(in) :: pavel, tavel, h2ovmr


      real (kind=kind_phys), dimension(max_nday_length, nlay, maxgas),intent(in) :: colamt

      integer, intent(in) :: nlay, nlp1, ix, idxday(nday_length), nday(fulljj), fulljj, blockjj


      

      real (kind=kind_phys), dimension(max_nday_length, 2),  intent(in) :: albbm, albdf

      real (kind=kind_phys), dimension(ix, fulljj), intent(in) :: cosz
      real (kind=kind_phys), dimension(max_nday_length), intent(in) :: sntz, &
         cf1, cf0, ssolar

!  ---  outputs:
      real (kind=kind_phys), dimension(max_nday_length, nlp1,nbdsw), intent(out) ::      &
     &       fxupc, fxdnc, fxup0, fxdn0

      real (kind=kind_phys), dimension(max_nday_length, 2), intent(out) :: sfbmc, sfdfc, &
     &       sfbm0, sfdf0

      real (kind=kind_phys), dimension(max_nday_length), intent(out) :: suvbfc, suvbf0, ftoadc,     &
     &       ftoauc, ftoau0, fsfcuc, fsfcu0, fsfcdc, fsfcd0

!  ---  locals:
      real (kind=kind_phys), dimension(nday_length, nlay, max_packs_nbands) :: &
         tauae, ssaae, asyae
      real (kind=kind_phys), dimension(:,:,:,:), allocatable :: cldfmc
      real (kind=kind_phys), dimension(nday_length, nlay, max_packs_nbands) :: taucw, ssacw, asycw
      real (kind=kind_phys), dimension(max_nday_length, nlay) :: cldfrc     
      integer, dimension(max_nday_length, nlay) :: indfor, indself,          &
     &       jp, jt, jt1
      integer :: laytrop(max_nday_length)
      real (kind=kind_phys), dimension(max_nday_length, nlay) :: &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac

      real (kind=kind_phys) :: ztau1, zssa1, zasy1, ztau0, zssa0,       &
     &       zasy0, zasy3, zssaw, zasyw, zgam1, zgam2, zgam3, zgam4,    &
     &       zc0, zc1, za1, za2, &
     &       zrk, zrk2, zrp, zrp1, zrm1, zb1, zb2, &
     &       zrpp, zrkg1, zrkg3, zrkg4, zexp1, zexm1, zexp2, zexm2,     &
     &       zexp3, zexp4, zden1, ze1r45, ftind, zrefb1,        &
     &       zrefd1, ztrab1, ztrad1, zr1, zr2, zr3, zr4, zr5, ztdbt0r,    &
     &       zt1, zt2, zt3, zf1, zf2, zldbt0, zfu, zfd, zrupbr, zrupbr1, zrupdr, zrupdr1, &
             ztdbtr, zrefbr, taurr, tauaer, ssaaer, coszr, zrefdr, &
             ztdnr, zrdndr, ztradr, ztrabr, zldbtr

      !real (kind=kind_phys), dimension(ix, small_ngptsw, fulljj) :: ztdbt0, zsolar, zfd0

      integer :: ib, ibd, jb, jg, k, kp, itind, ipt, jj, j1, j2, jg2, i, j, n, nn
      !real (kind=kind_phys), dimension(ix, nlp1, small_ngptsw, fulljj) :: zrdnd, ztdn
      integer :: async_id

      integer :: ng00, offset_ng
      real (kind=kind_phys), dimension(nday_length, nlay, small_ngpts) :: small_taur, &
     &       small_taug
      real (kind=kind_phys), dimension(nday_length, small_ngpts) :: small_sfluxzen
      real (kind=kind_phys), dimension(nday_length, nlp1, small_ngpts) :: zrefb, zrefd, ztrab,    &
     &       ztrad, zldbt, zfua, zfda, zrdnd, ztdn, ztdbt
      real (kind=kind_phys), dimension(nday_length, small_ngpts) :: ztdbt0, zsolar, zfd0
      integer, dimension(small_ngptsw, small_factor) :: ngbidx
      ng00 = maxval(ng)

!
!===> ...  begin here
!     
!  --- ... initialization of output fluxes
      !$acc parallel loop collapse(3) async(async_id)
      do ib = 1, nbdsw
         do k = 1, nlp1
            do n = 1, nday_length
               fxdnc(n, k,ib) = f_zero
               fxupc(n, k,ib) = f_zero
               fxdn0(n, k,ib) = f_zero
               fxup0(n, k,ib) = f_zero
            enddo
         end do
      end do
      !$acc parallel loop async(async_id)
      do n = 1, nday_length
         ftoadc(n) = f_zero
         ftoauc(n) = f_zero
         ftoau0(n) = f_zero
         fsfcuc(n) = f_zero
         fsfcu0(n) = f_zero
         fsfcdc(n) = f_zero
         fsfcd0(n) = f_zero

         !! --- ...  uv-b surface downward fluxes
         suvbfc(n)  = f_zero
         suvbf0(n)  = f_zero

         !! --- ...  output surface flux components
         sfbmc(n, 1) = f_zero
         sfbmc(n, 2) = f_zero
         sfdfc(n, 1) = f_zero
         sfdfc(n, 2) = f_zero
         sfbm0(n, 1) = f_zero
         sfbm0(n, 2) = f_zero
         sfdf0(n, 1) = f_zero
         sfdf0(n, 2) = f_zero
      end do

   do j2 = 1, npacks
      !$acc enter data create(tauae, ssaae, asyae) async(async_id)

      !$acc enter data create(taucw, ssacw, asycw, cldfrc) async(async_id)
      !$acc enter data create(laytrop,jp,jt,jt1,fac00,fac01,fac10,fac11,   &
      !$acc&     selffac,selffrac,indself,forfac,forfrac,indfor) async(async_id)


      !$acc enter data create(small_sfluxzen, small_taug, small_taur) async(async_id)
      !  --- ...  loop over all g-points in each band
      !call nvtxStartRange("sw_setcoef")
      call setcoef                                                    &
         !  ---  inputs:
         &     ( pavel,tavel,h2ovmr, nlay,nlp1, ix, nday, async_id, fulljj, blockjj,                              &
               map_nday_ipt, map_nday_jj, &
               nday_length, max_nday_length, jjoffset, &
         !  ---  outputs:
         &       laytrop,jp,jt,jt1,fac00,fac01,fac10,fac11,                 &
         &       selffac,selffrac,indself,forfac,forfrac,indfor             &
         &     )
      !call nvtxEndRange
      !call nvtxStartRange("sw_taumol")
      do n = 1, nbdsw
         jb = recipe(n, j2) ! jb is serial numbers of nbands
         if (jb .eq. 0) exit

         do i = 1, small_ngpts
            if (jb .eq. pack_idx(i, j2)) then
               offset_ng = i-1
               exit
            end if
         end do
         call taumol                                                     &
            !  ---  inputs:
            &     ( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
            &       forfac,forfrac,indfor,selffac,selffrac,indself, map_nday_ipt, map_nday_jj, &
                  nlay, ix, nday, async_id, fulljj, blockjj, jb, small_ngpts, j2, &
                  offset_ng, nday_length, max_nday_length, jjoffset,     &
            !  ---  outputs:
            &       small_sfluxzen, small_taug, small_taur                                       &
            &     )
         end do
         !call nvtxEndRange
         call cldprop                                                  &
            !  ---  inputs:
            &     ( cfrac,cliqp,reliq, &
                  cicep,reice,cdat1, &
                  cdat2,cdat3,cdat4,     &
            &       cf1, nlay, ipseed, ix, nday, idxday, &
                  map_nday_ipt, map_nday_jj, &
                  nday_length, max_nday_length, jjoffset, jbs_nday, &
                  recipe(:, j2), small_ngpts, packs_nbands_list(j2), offset_ng, &
                  async_id, fulljj, blockjj,                                    &
            !  ---  outputs:
            &       taucw, ssacw, asycw, cldfrc, cldfmc                        &
            &     )
      !$acc exit data delete(laytrop,jp,jt,jt1,fac00,fac01,fac10,fac11,   &
      !$acc&     selffac,selffrac,indself,forfac,forfrac,indfor) async(async_id)
      !call nvtxStartRange("sw_setaer")
      call setaer_sw_gpu                                                       &
      !  ---  inputs:
      &     ( plvl,plyr,prslk1,tvly,rhly,slmsk,tracer1, &
      xlon,xlat,        &
      &       myim,nlay,nlp1,lsswr,lslwr,me,myrank, ix, map_jj, map_i, nxptot, &
            idxday, map_nday_jj, jbs_nday, &
            nday_length, jjoffset, max_nday_length, recipe(:, j2), small_ngpts, packs_nbands_list(j2), offset_ng, max_packs_nbands, async_id,        &
      !  ---  outputs:
      &       tauae, ssaae, asyae                                              &
      !    &       faersw,faerlw,aerodp                                       &
      &     )
      !call nvtxEndRange
      !$acc enter data create(zrefb, zrefd, ztrab, ztrad, ztdbt, zldbt, zfd0, &
      !$acc&     ztdbt0, zsolar, zrdnd, ztdn, zfda, zfua) async(async_id)

      !call nvtxStartRange("sw_spcvrtc")
      !$acc parallel loop collapse(2) private(jj, jb, ib, ibd, za1, za2, &
      !$acc&         ztau1, zssa1, zasy1, zasy3, zgam1, zgam2, zgam3, &
      !$acc&         zgam4, zb1, zb2, ftind, itind, zrk, zrk2, zrp, zrp1, zrm1, &
      !$acc&         zrpp, zrkg1, zrkg3, zrkg4, zr1, zr2, zr3, zr4, zr5, zt1, &
      !$acc&         zt2, zt3, zexm1, zexp1, zexm2, zexp2, ze1r45, zden1, zexp3, &
      !$acc&         zexp4, ztdbtr, zrefdr, taurr, tauaer, ssaaer, coszr, zrefbr) async(async_id)
      do jg = 1, pack_size(j2) ! lab_do_jg
         do n = 1, nday_length
            jj = map_nday_jj(n) - jjoffset
            jb = pack_idx(jg, j2)
            ib = jb + 1 - nblow
            nn = reverse_recipe(jg, j2)
            ibd = idxsfc(jb)
            zsolar(n, jg) = ssolar(n) * small_sfluxzen(n, jg)

            !  --- ...  set up toa direct beam and surface values (beam and diff)
            ztdbtr = f_one
            ztdbt(n, nlp1, jg) = ztdbtr

            zldbt(n, 1, jg) = f_zero
            if (ibd /= 0) then
               zrefb(n, 1, jg) = albbm(n, ibd)
               zrefd(n, 1, jg) = albdf(n, ibd)
            else
               zrefb(n, 1, jg) = 0.5 * (albbm(n, 1) + albbm(n, 2))
               zrefd(n, 1, jg) = 0.5 * (albdf(n, 1) + albdf(n, 2))
            endif
            ztrab(n, 1, jg) = f_zero
            ztrad(n, 1, jg) = f_zero

         !  --- ...  compute clear-sky optical parameters, layer reflectance and transmittance
            j1 = idxday(n)
            coszr = cosz(j1, jj)
            ztdbt0r = f_one
            !$acc loop seq
            do k = nlay, 1, -1
               kp = k + 1
               taurr = small_taur(n, k,jg)
               tauaer = tauae(n, k,nn)
               ssaaer = ssaae(n, k,nn)

               ztau0 = max( ftiny, taurr+small_taug(n, k,jg)+tauaer )
               zssa0 = taurr + tauaer*ssaaer
               zasy0 = asyae(n, k,nn)*ssaaer*tauaer
               zssaw = min( oneminus, zssa0 / ztau0 )
               zasyw = zasy0 / max( ftiny, zssa0 )


               !  --- ...  delta scaling for clear-sky condition
               za1 = zasyw * zasyw
               za2 = zssaw * za1

               ztau1 = (f_one - za2) * ztau0
               zssa1 = (zssaw - za2) / (f_one - za2)
               !org      zasy1 = (zasyw - za1) / (f_one - za1)   ! this line is replaced by the next
               zasy1 = zasyw / (f_one + zasyw)         ! to reduce truncation error
               zasy3 = 0.75 * zasy1

               !  --- ...  general two-stream expressions
               if ( iswmode == 1 ) then
                  zgam1 = 1.75 - zssa1 * (f_one + zasy3)
                  zgam2 =-0.25 + zssa1 * (f_one - zasy3)
                  zgam3 = 0.5  - zasy3 * coszr
               elseif ( iswmode == 2 ) then               ! pifm
                  zgam1 = 2.0 - zssa1 * (1.25 + zasy3)
                  zgam2 = 0.75* zssa1 * (f_one- zasy1)
                  zgam3 = 0.5 - zasy3 * coszr
               elseif ( iswmode == 3 ) then               ! discrete ordinates
                  zgam1 = zsr3 * (2.0 - zssa1 * (1.0 + zasy1)) * 0.5
                  zgam2 = zsr3 * zssa1 * (1.0 - zasy1) * 0.5
                  zgam3 = (1.0 - zsr3 * zasy1 * coszr) * 0.5
               endif
               zgam4 = f_one - zgam3

               !  --- ...  compute homogeneous reflectance and transmittance

               if ( zssaw >= zcrit ) then    ! for conservative scattering
                  za1 = zgam1 * coszr - zgam3
                  za2 = zgam1 * ztau1

                  !  --- ...  use exponential lookup table for transmittance, or expansion
                  !           of exponential for low optical depth

                  zb1 = min ( ztau1*sntz(n) , 500.0 )
                  if ( zb1 <= od_lo ) then
                     zb2 = f_one - zb1 + 0.5*zb1*zb1
                  else
                     ftind = zb1 / (bpade + zb1)
                     itind = ftind*ntbmx + 0.5
                     zb2 = exp_tbl(itind)
                  endif

                  !      ...  collimated beam
                  zrefbr = max(f_zero, min(f_one,                          &
                  &                  (za2 - za1*(f_one - zb2))/(f_one + za2) ))
                  ztrab(n, kp, jg) = max(f_zero, min(f_one, f_one-zrefbr ))
                  zrefb(n, kp, jg) = zrefbr

                  !      ...  isotropic incidence
                  zrefdr = max(f_zero, min(f_one, za2/(f_one + za2) ))
                  ztrad(n, kp, jg) = max(f_zero, min(f_one, f_one-zrefdr ))
                  zrefd(n, kp, jg) = zrefdr

               else                          ! for non-conservative scattering
                  za1 = zgam1*zgam4 + zgam2*zgam3
                  za2 = zgam1*zgam3 + zgam2*zgam4
                  zrk = sqrt ( (zgam1 - zgam2) * (zgam1 + zgam2) )
                  zrk2= 2.0 * zrk

                  zrp  = zrk * coszr
                  zrp1 = f_one + zrp
                  zrm1 = f_one - zrp
                  zrpp = f_one - zrp*zrp
                  zrkg1= zrk + zgam1
                  zrkg3= zrk * zgam3
                  zrkg4= zrk * zgam4

                  zr1  = zrm1 * (za2 + zrkg3)
                  zr2  = zrp1 * (za2 - zrkg3)
                  zr3  = zrk2 * (zgam3 - za2*coszr)
                  zr4  = zrpp * zrkg1
                  zr5  = zrpp * (zrk - zgam1)

                  zt1  = zrp1 * (za1 + zrkg4)
                  zt2  = zrm1 * (za1 - zrkg4)
                  zt3  = zrk2 * (zgam4 + za1*coszr)

                  !  --- ...  use exponential lookup table for transmittance, or expansion
                  !           of exponential for low optical depth

                  zb1 = min ( zrk*ztau1, 500.0 )
                  if ( zb1 <= od_lo ) then
                     zexm1 = f_one - zb1 + 0.5*zb1*zb1
                  else
                     ftind = zb1 / (bpade + zb1)
                     itind = ftind*ntbmx + 0.5
                     zexm1 = exp_tbl(itind)
                  endif
                  zexp1 = f_one / zexm1

                  zb2 = min ( sntz(n)*ztau1, 500.0 )
                  if ( zb2 <= od_lo ) then
                     zexm2 = f_one - zb2 + 0.5*zb2*zb2
                  else
                     ftind = zb2 / (bpade + zb2)
                     itind = ftind*ntbmx + 0.5
                     zexm2 = exp_tbl(itind)
                  endif
                  zexp2 = f_one / zexm2
                  ze1r45 = zr4*zexp1 + zr5*zexm1

                  !      ...  collimated beam
                  if (ze1r45>=-eps1 .and. ze1r45<=eps1) then
                     zrefb(n, kp, jg) = eps1
                     ztrab(n, kp, jg) = zexm2
                  else
                     zden1 = zssa1 / ze1r45
                     zrefb(n, kp, jg) = max(f_zero, min(f_one,                        &
                     &                    (zr1*zexp1 - zr2*zexm1 - zr3*zexm2)*zden1 ))
                     ztrab(n, kp, jg) = max(f_zero, min(f_one, zexm2*(f_one           &
                     &                  - (zt1*zexp1 - zt2*zexm1 - zt3*zexp2)*zden1) ))
                  endif

                  !      ...  diffuse beam
                  zden1 = zr4 / (ze1r45 * zrkg1)
                  zrefd(n, kp, jg) = max(f_zero, min(f_one,                          &
                  &                  zgam2*(zexp1 - zexm1)*zden1 ))
                  ztrad(n, kp, jg) = max(f_zero, min(f_one, zrk2*zden1 ))
               endif    ! end if_zssaw_block

               !  --- ...  direct beam transmittance. use exponential lookup table
               !           for transmittance, or expansion of exponential for low
               !           optical depth

               zr1 = ztau1 * sntz(n)
               if ( zr1 <= od_lo ) then
                  zexp3 = f_one - zr1 + 0.5*zr1*zr1
               else
                  ftind = zr1 / (bpade + zr1)
                  itind = max(0, min(ntbmx, int(0.5+ntbmx*ftind) ))
                  zexp3 = exp_tbl(itind)
               endif

               ztdbtr  = zexp3 * ztdbtr
               ztdbt(n, k, jg) = ztdbtr
               zldbt(n, kp, jg) = zexp3

               !  --- ...  pre-delta-scaling clear and cloudy direct beam transmittance
               !           (must use 'orig', unscaled cloud optical depth)

               zr1 = ztau0 * sntz(n)
               if ( zr1 <= od_lo ) then
                  zexp4 = f_one - zr1 + 0.5*zr1*zr1
               else
                  ftind = zr1 / (bpade + zr1)
                  itind = max(0, min(ntbmx, int(0.5+ntbmx*ftind) ))
                  zexp4 = exp_tbl(itind)
               endif

               ztdbt0r = zexp4 * ztdbt0r
            enddo    ! end do_k_loop
            ztdbt0(n, jg) = ztdbt0r
         end do
      end do
               
               ! call swflux
               !  --- ...  link lowest layer with surface
      !$acc parallel loop collapse(2) private(jb, ib, ibd, kp, zden1, ztdbtr, &
      !$acc&         zfu, zfd, zrupbr, zrupbr1, zrupdr, zrupdr1, ztdnr, zrdndr, ztradr, zrefdr) async(async_id)
      do jg = 1, pack_size(j2) ! lab_do_jg
         do n = 1, nday_length
            jb = pack_idx(jg, j2)
            ib = jb + 1 - nblow
            ibd = idxsfc(jb)
            !  --- ...  upper boundary conditions

            ztdn (n, nlp1, jg) = f_one
            zrdnd(n, nlp1, jg) = f_zero
            ztdnr = ztrab(n, nlp1, jg)
            ztdn (n, nlay, jg) = ztdnr
            zrdndr = zrefd(n, nlp1, jg)
            zrdnd(n, nlay, jg) = zrdndr
            !  --- ...  pass from top to bottom
            !$acc loop seq
            do k = nlay, 2, -1
               ztdbtr = ztdbt(n, k, jg)
               ztradr = ztrad(n, k, jg)
               zrefdr = zrefd(n, k, jg)
               zden1 = f_one / (f_one - zrefdr*zrdndr)
               ztdnr = ztdbtr*ztrab(n, k, jg) &
                                 + ( ztradr *                 &
               &                 ( (ztdnr - ztdbtr) + ztdbtr *              &
               &                 zrefb(n, k, jg)*zrdndr )) * zden1
               ztdn (n, k-1, jg) = ztdnr
               zrdndr = zrefdr + ztradr &
                              *ztradr*zrdndr*zden1
               zrdnd(n, k-1, jg) =zrdndr
            enddo
            
            zrupbr = zrefb(n, 1, jg)        ! direct beam
            zrupdr = zrefd(n, 1, jg)        ! diffused
            ztdbtr = ztdbt(n, 1, jg)
            !  --- ...  up and down-welling fluxes at levels
            zden1 = f_one / (f_one - zrdndr*zrupdr)
            zfu = ( ztdbtr*zrupbr +                                &
            &             (ztdnr - ztdbtr) &
                           *zrupdr ) * zden1
            zfd = ztdbtr + ( ztdnr &
                           - ztdbtr +                    &
            &             ztdbtr*zrupbr &
                           *zrdndr ) * zden1
            ! end call swflux
            zfua(n, 1, jg) = zfu
            zfda(n, 1, jg) = zfd
            zfd0(n, jg) = zfd
            !  --- ...  pass from bottom to top
            !$acc loop seq
            do k = 1, nlay
               kp = k + 1
               ztdnr = ztdn(n, kp, jg)
               zrdndr = zrdnd(n, kp, jg)
               ztdbtr = ztdbt(n, kp, jg)
               ztradr = ztrad(n, kp, jg)
               zrefdr = zrefd(n, kp, jg)
               zden1 = f_one / ( f_one - zrupdr*zrefdr )
               zrupbr1 = zrefb(n, kp, jg) + ( ztradr *                         &
               &                ( (ztrab(n, kp, jg) - zldbt(n, kp, jg)) &
                              *zrupdr +              &
               &                zldbt(n, kp, jg)*zrupbr) ) * zden1
               zrupdr1 = zrefdr + ztradr &
                              *ztradr*zrupdr*zden1
               
               zden1 = f_one / (f_one - zrdndr*zrupdr1)
               zfu = ( ztdbtr*zrupbr1 +                                &
               &             (ztdnr - ztdbtr) &
                           *zrupdr1 ) * zden1
               zfd = ztdbtr + ( ztdnr &
                           - ztdbtr +                    &
               &             ztdbtr*zrupbr1 &
                           *zrdndr ) * zden1
               zfua(n, kp, jg) = zfu
               zfda(n, kp, jg) = zfd
               ! end call swflux
               !  --- ...  compute upward and downward fluxes at levels
               zrupbr = zrupbr1
               zrupdr = zrupdr1
            end do
         end do
      end do

      !$acc parallel loop collapse(2) private(jb, ib) async(async_id)
      do k = 1, nlp1
         do n = 1, nday_length
            !$acc loop seq
            do jg = 1, pack_size(j2) ! lab_do_jg
               jb = pack_idx(jg, j2)
               ib = jb + 1 - nblow
               !  --- ...  compute upward and downward fluxes at levels
               fxup0(n, k,ib) = fxup0(n, k,ib) + zsolar(n, jg)*zfua(n, k, jg)
               fxdn0(n, k,ib) = fxdn0(n, k,ib) + zsolar(n, jg)*zfda(n, k, jg)
            end do
         end do
      end do

               !! --- ...  surface downward beam/diffused flux components
      !$acc parallel loop private(jb, ib, ibd, zf1, zf2, zb1, zb2) async(async_id)
      do n = 1, nday_length
         !$acc loop seq 
         do jg = 1, pack_size(j2) ! lab_do_jg
            jb = pack_idx(jg, j2)
            ib = jb + 1 - nblow
            ibd = idxsfc(jb)
            zb1 = zsolar(n, jg)*ztdbt0(n, jg)
            zb2 = zsolar(n, jg)*(zfd0(n, jg) - ztdbt0(n, jg))
            if (ibd /= 0) then
               sfbm0(n, ibd) = sfbm0(n, ibd) + zb1
               sfdf0(n, ibd) = sfdf0(n, ibd) + zb2
            else
               zf1 = 0.5 * zb1
               zf2 = 0.5 * zb2
               sfbm0(n, 1) = sfbm0(n, 1) + zf1
               sfdf0(n, 1) = sfdf0(n, 1) + zf2
               sfbm0(n, 2) = sfbm0(n, 2) + zf1
               sfdf0(n, 2) = sfdf0(n, 2) + zf2
            endif
            !       sfbm0(ibd) = sfbm0(ibd) + zsolar*ztdbt0
            !       sfdf0(ibd) = sfdf0(ibd) + zsolar*(zfd(1) - ztdbt0)
         end do
      end do

               !  --- ...  compute total sky optical parameters, layer reflectance and transmittance
      !$acc parallel loop collapse(2) private(jj, zb1, zb2, zf1, zf2, kp, zc0, zc1, ztau0, &
      !$acc&         zssa0, zasy0, zldbt0, ftind, itind, zssaw, zasyw, za1, za2, &
      !$acc&         ztau1, zssa1, zasy1, zasy3, zgam1, zgam2, zgam3, zgam4, &
      !$acc&         zrefb1, zrefd1, ztrab1, ztrad1, zrk, zrk2, zrp, zrp1, zrm1, &
      !$acc&         zrpp, zrkg1, zrkg3, zrkg4, zr1, zr2, zr3, zr4, zr5, &
      !$acc&         zt1, zt2, zt3, zexm1, zexp1, zexm2, zexp2, ze1r45, zden1, &
      !$acc&         zexp3, zexp4, j1, jb, ib, ibd, ztdbt0r, zrefbr, zrefdr, ztrabr, &
      !$acc&         ztradr, zldbtr, taurr, tauaer, ssaaer) async(async_id)
      do jg = 1, pack_size(j2) ! lab_do_jg
         do n = 1, nday_length
            jj = map_nday_jj(n) - jjoffset
            jb = pack_idx(jg, j2)
            ib = jb + 1 - nblow
            nn = reverse_recipe(jg, j2)
            ibd = idxsfc(jb)
            j1 = idxday(n)
            if ( cf1(n) > eps ) then

               !  --- ...  set up toa direct beam and surface values (beam and diff)
               ztdbt0r = f_one
               zldbt(n, 1, jg) = f_zero
               ztdbtr = ztdbt(n, nlp1, jg)
               !$acc loop seq
               do k = nlay, 1, -1
                  kp = k + 1
                  zc0 = f_one - cldfrc(n, k)
                  zc1 = cldfrc(n, k)
                  taurr = small_taur(n, k,jg)
                  tauaer = tauae(n, k,nn)
                  ssaaer = ssaae(n, k,nn)

                  !  --- ...  saving clear-sky quantities for later total-sky usage
                  ztau0 = max( ftiny, taurr+small_taug(n, k,jg)+tauaer )
                  zssa0 = taurr + tauaer*ssaaer
                  zasy0 = asyae(n, k,nn)*ssaaer*tauaer
                  zr1 = ztau0 * sntz(n)
                  if ( zr1 <= od_lo ) then
                     zldbt0 = f_one - zr1 + 0.5*zr1*zr1
                  else
                     ftind = zr1 / (bpade + zr1)
                     itind = max(0, min(ntbmx, int(0.5+ntbmx*ftind) ))
                     zldbt0 = exp_tbl(itind)
                  endif

                  zldbtr = zldbt(n, kp, jg)
                  if ( zc1 > ftiny ) then          ! it is a cloudy-layer
                     ztau0 = ztau0 + taucw(n, k,nn)
                     zssa0 = zssa0 + ssacw(n, k,nn)
                     zasy0 = zasy0 + asycw(n, k,nn)
                     zssaw = min(oneminus, zssa0 / ztau0)
                     zasyw = zasy0 / max(ftiny, zssa0)

                     !  --- ...  delta scaling for total-sky condition
                     za1 = zasyw * zasyw
                     za2 = zssaw * za1

                     ztau1 = (f_one - za2) * ztau0
                     zssa1 = (zssaw - za2) / (f_one - za2)
                     !org          zasy1 = (zasyw - za1) / (f_one - za1)
                     zasy1 = zasyw / (f_one + zasyw)
                     zasy3 = 0.75 * zasy1

                     !  --- ...  general two-stream expressions
                     if ( iswmode == 1 ) then
                        zgam1 = 1.75 - zssa1 * (f_one + zasy3)
                        zgam2 =-0.25 + zssa1 * (f_one - zasy3)
                        zgam3 = 0.5  - zasy3 * cosz(j1, jj)
                     elseif ( iswmode == 2 ) then               ! pifm
                        zgam1 = 2.0 - zssa1 * (1.25 + zasy3)
                        zgam2 = 0.75* zssa1 * (f_one- zasy1)
                        zgam3 = 0.5 - zasy3 * cosz(j1, jj)
                     elseif ( iswmode == 3 ) then               ! discrete ordinates
                        zgam1 = zsr3 * (2.0 - zssa1 * (1.0 + zasy1)) * 0.5
                        zgam2 = zsr3 * zssa1 * (1.0 - zasy1) * 0.5
                        zgam3 = (1.0 - zsr3 * zasy1 * cosz(j1, jj)) * 0.5
                     endif
                     zgam4 = f_one - zgam3

                     zrefb1 = zrefb(n, kp, jg)
                     zrefd1 = zrefd(n, kp, jg)
                     ztrab1 = ztrab(n, kp, jg)
                     ztrad1 = ztrad(n, kp, jg)

                     !  --- ...  compute homogeneous reflectance and transmittance

                     if ( zssaw >= zcrit ) then    ! for conservative scattering
                        za1 = zgam1 * cosz(j1, jj) - zgam3
                        za2 = zgam1 * ztau1

                        !  --- ...  use exponential lookup table for transmittance, or expansion
                        !           of exponential for low optical depth

                        zb1 = min ( ztau1*sntz(n) , 500.0 )
                        if ( zb1 <= od_lo ) then
                           zb2 = f_one - zb1 + 0.5*zb1*zb1
                        else
                           ftind = zb1 / (bpade + zb1)
                           itind = ftind*ntbmx + 0.5
                           zb2 = exp_tbl(itind)
                        endif

                        !      ...  collimated beam
                        zrefbr = max(f_zero, min(f_one,                      &
                        &                      (za2 - za1*(f_one - zb2))/(f_one + za2) ))
                        ztrabr = max(f_zero, min(f_one, f_one-zrefbr))

                        !      ...  isotropic incidence
                        zrefdr = max(f_zero, min(f_one, za2 / (f_one+za2) ))
                        ztradr = max(f_zero, min(f_one, f_one - zrefdr ))

                     else                          ! for non-conservative scattering
                        za1 = zgam1*zgam4 + zgam2*zgam3
                        za2 = zgam1*zgam3 + zgam2*zgam4
                        zrk = sqrt ( (zgam1 - zgam2) * (zgam1 + zgam2) )
                        zrk2= 2.0 * zrk

                        zrp  = zrk * cosz(j1, jj)
                        zrp1 = f_one + zrp
                        zrm1 = f_one - zrp
                        zrpp = f_one - zrp*zrp
                        zrkg1= zrk + zgam1
                        zrkg3= zrk * zgam3
                        zrkg4= zrk * zgam4

                        zr1  = zrm1 * (za2 + zrkg3)
                        zr2  = zrp1 * (za2 - zrkg3)
                        zr3  = zrk2 * (zgam3 - za2*cosz(j1, jj))
                        zr4  = zrpp * zrkg1
                        zr5  = zrpp * (zrk - zgam1)

                        zt1  = zrp1 * (za1 + zrkg4)
                        zt2  = zrm1 * (za1 - zrkg4)
                        zt3  = zrk2 * (zgam4 + za1*cosz(j1, jj))

                        !  --- ...  use exponential lookup table for transmittance, or expansion
                        !           of exponential for low optical depth

                        zb1 = min ( zrk*ztau1, 500.0 )
                        if ( zb1 <= od_lo ) then
                           zexm1 = f_one - zb1 + 0.5*zb1*zb1
                        else
                           ftind = zb1 / (bpade + zb1)
                           itind = ftind*ntbmx + 0.5
                           zexm1 = exp_tbl(itind)
                        endif
                        zexp1 = f_one / zexm1

                        zb2 = min ( ztau1*sntz(n), 500.0 )
                        if ( zb2 <= od_lo ) then
                           zexm2 = f_one - zb2 + 0.5*zb2*zb2
                        else
                           ftind = zb2 / (bpade + zb2)
                           itind = ftind*ntbmx + 0.5
                           zexm2 = exp_tbl(itind)
                        endif
                        zexp2 = f_one / zexm2
                        ze1r45 = zr4*zexp1 + zr5*zexm1

                        !      ...  collimated beam
                        if ( ze1r45>=-eps1 .and. ze1r45<=eps1 ) then
                           zrefbr = eps1
                           ztrabr = zexm2
                        else
                           zden1 = zssa1 / ze1r45
                           zrefbr = max(f_zero, min(f_one,                    &
                           &                        (zr1*zexp1-zr2*zexm1-zr3*zexm2)*zden1 ))
                           ztrabr = max(f_zero, min(f_one, zexm2*(f_one -     &
                           &                        (zt1*zexp1-zt2*zexm1-zt3*zexp2)*zden1) ))
                        endif

                        !      ...  diffuse beam
                        zden1 = zr4 / (ze1r45 * zrkg1)
                        zrefdr = max(f_zero, min(f_one,                      &
                        &                      zgam2*(zexp1 - zexm1)*zden1 ))
                        ztradr = max(f_zero, min(f_one, zrk2*zden1 ))
                     endif    ! end if_zssaw_block

                     !  --- ...  combine clear and cloudy contributions for total sky
                     !           and calculate direct beam transmittances

                     zrefb(n, kp, jg) = zc0*zrefb1 + zc1*zrefbr
                     zrefd(n, kp, jg) = zc0*zrefd1 + zc1*zrefdr
                     ztrab(n, kp, jg) = zc0*ztrab1 + zc1*ztrabr
                     ztrad(n, kp, jg) = zc0*ztrad1 + zc1*ztradr

                     !  --- ...  direct beam transmittance. use exponential lookup table
                     !           for transmittance, or expansion of exponential for low
                     !           optical depth

                     zr1 = ztau1 * sntz(n)
                     if ( zr1 <= od_lo ) then
                        zexp3 = f_one - zr1 + 0.5*zr1*zr1
                     else
                        ftind = zr1 / (bpade + zr1)
                        itind = max(0, min(ntbmx, int(0.5+ntbmx*ftind) ))
                        zexp3 = exp_tbl(itind)
                     endif
                     
                     zldbtr = zc0*zldbtr + zc1*zexp3
                     ztdbtr = zldbtr * ztdbtr
                     zldbt(n, kp, jg) = zldbtr

                     !  --- ...  pre-delta-scaling clear and cloudy direct beam transmittance
                     !           (must use 'orig', unscaled cloud optical depth)

                     zr1 = ztau0 * sntz(n)
                     if ( zr1 <= od_lo ) then
                        zexp4 = f_one - zr1 + 0.5*zr1*zr1
                     else
                        ftind = zr1 / (bpade + zr1)
                        itind = max(0, min(ntbmx, int(0.5+ntbmx*ftind) ))
                        zexp4 = exp_tbl(itind)
                     endif

                     ztdbt0r = (zc0*zldbt0 + zc1*zexp4) * ztdbt0r

                  else     ! if_zc1_block  ---  it is a clear layer

                     !  --- ...  direct beam transmittance
                     ztdbtr = zldbtr * ztdbtr

                     !  --- ...  pre-delta-scaling clear and cloudy direct beam transmittance
                     ztdbt0r = zldbt0 * ztdbt0r

                  endif    ! end if_zc1_block
                  ztdbt(n, k, jg) = ztdbtr
               enddo   ! end do_k_loop
               ztdbt0(n, jg) = ztdbt0r
            end if
         end do
      end do
      !$acc exit data delete(small_sfluxzen, small_taug, small_taur) async(async_id)

                  !  --- ...  perform vertical quadrature
      !$acc parallel loop collapse(2) private(j1, kp, zden1, zfu, zfd, &
      !$acc&         jb, ib, zrupbr, zrupbr1, zrupdr, zrupdr1, ztdnr, zrdndr, &
      !$acc&         ztdbtr, ztradr, zrefdr, zldbtr) async(async_id)
      do jg = 1, pack_size(j2) ! lab_do_jg
         do n = 1, nday_length
            jb = pack_idx(jg, j2)
            ib = jb + 1 - nblow
            if ( cf1(n) > eps ) then
               ! call swflux

               !  --- ...  upper boundary conditions

               ztdn (n, nlp1, jg) = f_one
               zrdnd(n, nlp1, jg) = f_zero
               ztdnr = ztrab(n, nlp1, jg)
               ztdn (n, nlay, jg) = ztdnr
               zrdndr = zrefd(n, nlp1, jg)
               zrdnd(n, nlay, jg) = zrdndr

               !  --- ...  pass from top to bottom
               !$acc loop seq
               do k = nlay, 2, -1
                  ztdbtr = ztdbt(n, k, jg)
                  ztradr = ztrad(n, k, jg)
                  zrefdr = zrefd(n, k, jg)
                  zden1 = f_one / (f_one - zrefdr*zrdndr)
                  ztdnr = ztdbtr*ztrab(n, k, jg) &
                                    + ( ztradr *                 &
                  &                 ( (ztdnr - ztdbtr) &
                                    + ztdbtr *              &
                  &                 zrefb(n, k, jg)*zrdndr )) * zden1
                  ztdn (n, k-1, jg) = ztdnr
                  zrdndr = zrefdr + ztradr &
                                    *ztradr*zrdndr*zden1
                  zrdnd(n, k-1, jg) = zrdndr
               enddo

               !  --- ...  link lowest layer with surface

               zrupbr = zrefb(n, 1, jg)        ! direct beam
               zrupdr = zrefd(n, 1, jg)        ! diffused
               ztdbtr = ztdbt(n, 1, jg)
               
               zden1 = f_one / (f_one - zrdndr*zrupdr)
               zfu = ( ztdbtr*zrupbr +                                &
               &             (ztdnr - ztdbtr) &
                              *zrupdr ) * zden1
               zfd = ztdbtr + ( ztdnr &
                              - ztdbtr +                    &
               &             ztdbtr*zrupbr &
                              *zrdndr ) * zden1
                              ! end call swflux
               zfua(n, 1, jg) = zfu
               zfda(n, 1, jg) = zfd

               zfd0(n, jg) = zfd
               !  --- ...  pass from bottom to top
               !$acc loop seq
               do k = 1, nlay
                  kp = k + 1
                  zrdndr = zrdnd(n, kp, jg)
                  ztdnr = ztdn(n, kp, jg)
                  ztdbtr = ztdbt(n, kp, jg)
                  ztradr = ztrad(n, kp, jg)
                  zrefdr = zrefd(n, kp, jg)
                  zldbtr = zldbt(n, kp, jg)
                  zden1 = f_one / ( f_one - zrupdr*zrefdr )
                  zrupbr1 = zrefb(n, kp, jg) + ( ztradr *                         &
                  &                ( (ztrab(n, kp, jg) - zldbtr)*zrupdr +              &
                  &                zldbtr*zrupbr) ) * zden1
                  zrupdr1 = zrefdr + ztradr &
                                 *ztradr*zrupdr*zden1

               !  --- ...  up and down-welling fluxes at levels
                  zden1 = f_one / (f_one - zrdndr*zrupdr1)
                  zfu = ( ztdbtr*zrupbr1 +                                &
                  &             (ztdnr - ztdbtr) &
                                 *zrupdr1 ) * zden1
                  zfd = ztdbtr + ( ztdnr &
                                 - ztdbtr +                    &
                  &             ztdbtr*zrupbr1 &
                                 *zrdndr ) * zden1
                                 ! end call swflux
                  zfua(n, kp, jg) = zfu
                  zfda(n, kp, jg) = zfd
                  zrupbr = zrupbr1
                  zrupdr = zrupdr1
               end do
            end if
         end do
      end do
      !$acc parallel loop collapse(2) private(jb, ib) async(async_id)
      do k = 1, nlp1
         do n = 1, nday_length
            !$acc loop seq
            do jg = 1, pack_size(j2) ! lab_do_jg
               jb = pack_idx(jg, j2)
               ib = jb + 1 - nblow
               !  --- ...  compute upward and downward fluxes at levels
               fxupc(n, k,ib) = fxupc(n, k,ib) + zsolar(n, jg)*zfua(n, k, jg)
               fxdnc(n, k,ib) = fxdnc(n, k,ib) + zsolar(n, jg)*zfda(n, k, jg)
            end do
         end do
      end do
                  !! --- ...  surface downward beam/diffused flux components
      !$acc parallel loop private(jb, ibd, zb1, zb2, zf1, zf2) async(async_id)
      do n = 1, nday_length
         !$acc loop seq
         do jg = 1, pack_size(j2) ! lab_do_jg
            jb = pack_idx(jg, j2)
            ibd = idxsfc(jb)
            if ( cf1(n) > eps ) then
               zb1 = zsolar(n, jg)*ztdbt0(n, jg)
               zb2 = zsolar(n, jg)*(zfd0(n, jg) - ztdbt0(n, jg))

               if (ibd /= 0) then
                  sfbmc(n, ibd) = sfbmc(n, ibd) + zb1
                  sfdfc(n, ibd) = sfdfc(n, ibd) + zb2
               else
                  zf1 = 0.5 * zb1
                  zf2 = 0.5 * zb2
                  sfbmc(n, 1) = sfbmc(n, 1) + zf1
                  sfdfc(n, 1) = sfdfc(n, 1) + zf2
                  sfbmc(n, 2) = sfbmc(n, 2) + zf1
                  sfdfc(n, 2) = sfdfc(n, 2) + zf2
               endif
               !         sfbmc(ibd) = sfbmc(ibd) + zsolar*ztdbt0
               !         sfdfc(ibd) = sfdfc(ibd) + zsolar*(zfd(1) - ztdbt0)

            endif      ! end if_cf1_block

         enddo  ! lab_do_jg
      end do
      !call nvtxEndRange
      !$acc exit data delete(tauae, ssaae, asyae) async(async_id)
      !$acc exit data delete(taucw, ssacw, asycw, cldfrc) async(async_id)
      !$acc exit data delete(zrefb, zrefd, ztrab, ztrad, ztdbt, zldbt, zfd0, &
      !$acc&     ztdbt0, zsolar, zrdnd, ztdn, zfda, zfua) async(async_id)
   end do


            !  --- ...  end of g-point loop
      !$acc parallel loop async(async_id)
      do n = 1, nday_length
         !$acc loop seq 
         do ib = 1, nbdsw
            ftoadc(n) = ftoadc(n) + fxdn0(n, nlp1,ib)
            ftoau0(n) = ftoau0(n) + fxup0(n, nlp1,ib)
            fsfcu0(n) = fsfcu0(n) + fxup0(n, 1,ib)
            fsfcd0(n) = fsfcd0(n) + fxdn0(n, 1,ib)
         enddo
      end do

      !! --- ...  uv-b surface downward flux
      !$acc parallel loop collapse(3) async(async_id)
      do ib = 1, nbdsw
         do k = 1, nlp1
            do n = 1, nday_length
               if ( cf1(n) <= eps ) then       ! clear column, set total-sky=clear-sky fluxes
                  fxupc(n, k,ib) = fxup0(n, k,ib)
                  fxdnc(n, k,ib) = fxdn0(n, k,ib)
               else                             ! cloudy column, compute total-sky fluxes
                  fxupc(n, k,ib) = cf1(n)*fxupc(n, k,ib) + cf0(n)*fxup0(n, k,ib)
                  fxdnc(n, k,ib) = cf1(n)*fxdnc(n, k,ib) + cf0(n)*fxdn0(n, k,ib)
               end if
            enddo
         end do
      end do
      !$acc parallel loop private(ibd) async(async_id)
      do n = 1, nday_length
         ibd = nuvb - nblow + 1
         suvbf0(n) = fxdn0(n, 1,ibd)
         if ( cf1(n) <= eps ) then       ! clear column, set total-sky=clear-sky fluxes
            ftoauc(n) = ftoau0(n)
            fsfcuc(n) = fsfcu0(n)
            fsfcdc(n) = fsfcd0(n)

            !! --- ...  surface downward beam/diffused flux components
            sfbmc(n, 1) = sfbm0(n, 1)
            sfdfc(n, 1) = sfdf0(n, 1)
            sfbmc(n, 2) = sfbm0(n, 2)
            sfdfc(n, 2) = sfdf0(n, 2)

            !! --- ...  uv-b surface downward flux
            suvbfc(n) = suvbf0(n)
         else                        ! cloudy column, compute total-sky fluxes
            !! --- ...  uv-b surface downward flux
            ibd = nuvb - nblow + 1
            suvbfc(n) = fxdnc(n, 1,ibd)

            !! --- ...  surface downward beam/diffused flux components
            sfbmc(n, 1) = cf1(n)*sfbmc(n, 1) + cf0(n)*sfbm0(n, 1)
            sfbmc(n, 2) = cf1(n)*sfbmc(n, 2) + cf0(n)*sfbm0(n, 2)
            sfdfc(n, 1) = cf1(n)*sfdfc(n, 1) + cf0(n)*sfdf0(n, 1)
            sfdfc(n, 2) = cf1(n)*sfdfc(n, 2) + cf0(n)*sfdf0(n, 2)
         endif    ! end if_cf1_block
      end do
      !$acc parallel loop  async(async_id)
      do n = 1, nday_length
         !$acc loop seq
         do ib = 1, nbdsw
            if ( cf1(n) > eps ) then                        ! cloudy column, compute total-sky fluxes
               ftoauc(n) = ftoauc(n) + fxupc(n, nlp1,ib)
               fsfcuc(n) = fsfcuc(n) + fxupc(n, 1,ib)
               fsfcdc(n) = fsfcdc(n) + fxdnc(n, 1,ib)
            end if
         enddo
      end do

      return
!...................................
      end subroutine spcvrtc
!-----------------------------------


!-----------------------------------
      subroutine spcvrtm                                                &
!...................................
!  ---  inputs:
     &     ( ssolar,cosz,sntz,albbm,albdf,sfluxzen,cldfmc,              &
     &       cf1,cf0,taug,taur,tauae,ssaae,asyae,taucw,ssacw,asycw,     &
     &       nlay, nlp1,                                                &
!  ---  outputs:
     &       fxupc,fxdnc,fxup0,fxdn0,                                   &
     &       ftoauc,ftoau0,ftoadc,fsfcuc,fsfcu0,fsfcdc,fsfcd0,          &
     &       sfbmc,sfdfc,sfbm0,sfdf0,suvbfc,suvbf0                      &
     &     )

!  ===================  program usage description  ===================  !
!                                                                       !
!   purpose:  computes the shortwave radiative fluxes using two-stream  !
!             method of h. barker and mcica, the monte-carlo independent!
!             column approximation, for the representation of sub-grid  !
!             cloud variability (i.e. cloud overlap).                   !
!                                                                       !
!   subprograms called:  swflux                                         !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                        size  !
!    ssolar  - real, incoming solar flux at top                    1    !
!    cosz    - real, cosine solar zenith angle                     1    !
!    sntz    - real, secant solar zenith angle                     1    !
!    albbm   - real, surface albedo for direct beam radiation      2    !
!    albdf   - real, surface albedo for diffused radiation         2    !
!    sfluxzen- real, spectral distribution of incoming solar flux ngptsw!
!    cldfmc  - real, layer cloud fraction for g-point        nlay*ngptsw!
!    cf1     - real, >0: cloudy sky, otherwise: clear sky          1    !
!    cf0     - real, =1-cf1                                        1    !
!    taug    - real, spectral optical depth for gases        nlay*ngptsw!
!    taur    - real, optical depth for rayleigh scattering   nlay*ngptsw!
!    tauae   - real, aerosols optical depth                  nlay*nbdsw !
!    ssaae   - real, aerosols single scattering albedo       nlay*nbdsw !
!    asyae   - real, aerosols asymmetry factor               nlay*nbdsw !
!    taucw   - real, weighted cloud optical depth            nlay*nbdsw !
!    ssacw   - real, weighted cloud single scat albedo       nlay*nbdsw !
!    asycw   - real, weighted cloud asymmetry factor         nlay*nbdsw !
!    nlay,nlp1 - integer,  number of layers/levels                 1    !
!                                                                       !
!  output variables:                                                    !
!    fxupc   - real, tot sky upward flux                     nlp1*nbdsw !
!    fxdnc   - real, tot sky downward flux                   nlp1*nbdsw !
!    fxup0   - real, clr sky upward flux                     nlp1*nbdsw !
!    fxdn0   - real, clr sky downward flux                   nlp1*nbdsw !
!    ftoauc  - real, tot sky toa upwd flux                         1    !
!    ftoau0  - real, clr sky toa upwd flux                         1    !
!    ftoadc  - real, toa downward (incoming) solar flux            1    !
!    fsfcuc  - real, tot sky sfc upwd flux                         1    !
!    fsfcu0  - real, clr sky sfc upwd flux                         1    !
!    fsfcdc  - real, tot sky sfc dnwd flux                         1    !
!    fsfcd0  - real, clr sky sfc dnwd flux                         1    !
!    sfbmc   - real, tot sky sfc dnwd beam flux (nir/uv+vis)       2    !
!    sfdfc   - real, tot sky sfc dnwd diff flux (nir/uv+vis)       2    !
!    sfbm0   - real, clr sky sfc dnwd beam flux (nir/uv+vis)       2    !
!    sfdf0   - real, clr sky sfc dnwd diff flux (nir/uv+vis)       2    !
!    suvbfc  - real, tot sky sfc dnwd uv-b flux                    1    !
!    suvbf0  - real, clr sky sfc dnwd uv-b flux                    1    !
!                                                                       !
!  internal variables:                                                  !
!    zrefb   - real, direct beam reflectivity for clear/cloudy    nlp1  !
!    zrefd   - real, diffuse reflectivity for clear/cloudy        nlp1  !
!    ztrab   - real, direct beam transmissivity for clear/cloudy  nlp1  !
!    ztrad   - real, diffuse transmissivity for clear/cloudy      nlp1  !
!    zldbt   - real, layer beam transmittance for clear/cloudy    nlp1  !
!    ztdbt   - real, lev total beam transmittance for clr/cld     nlp1  !
!                                                                       !
!  control parameters in module "physpara"                              !
!    iswmode - control flag for 2-stream transfer schemes               !
!              = 1 delta-eddington    (joseph et al., 1976)             !
!              = 2 pifm               (zdunkowski et al., 1980)         !
!              = 3 discrete ordinates (liou, 1973)                      !
!                                                                       !
!  *******************************************************************  !
!  original code description                                            !
!                                                                       !
!  method:                                                              !
!  -------                                                              !
!     standard delta-eddington, p.i.f.m., or d.o.m. layer calculations. !
!     kmodts  = 1 eddington (joseph et al., 1976)                       !
!             = 2 pifm (zdunkowski et al., 1980)                        !
!             = 3 discrete ordinates (liou, 1973)                       !
!                                                                       !
!  modifications:                                                       !
!  --------------                                                       !
!   original: h. barker                                                 !
!   revision: merge with rrtmg_sw: j.-j.morcrette, ecmwf, feb 2003      !
!   revision: add adjustment for earth/sun distance:mjiacono,aer,oct2003!
!   revision: bug fix for use of palbp and palbd: mjiacono, aer, nov2003!
!   revision: bug fix to apply delta scaling to clear sky: aer, dec2004 !
!   revision: code modified so that delta scaling is not done in cloudy !
!             profiles if routine cldprop is used; delta scaling can be !
!             applied by swithcing code below if cldprop is not used to !
!             get cloud properties. aer, jan 2005                       !
!   revision: uniform formatting for rrtmg: mjiacono, aer, jul 2006     !
!   revision: use exponential lookup table for transmittance: mjiacono, !
!             aer, aug 2007                                             !
!                                                                       !
!  *******************************************************************  !
!  ======================  end of description block  =================  !

!  ---  constant parameters:
      real (kind=kind_phys), parameter :: zcrit = 0.9999995 ! thresold for conservative scattering
      real (kind=kind_phys), parameter :: zsr3  = sqrt(3.0)
      real (kind=kind_phys), parameter :: od_lo = 0.06
      real (kind=kind_phys), parameter :: eps1  = 1.0e-8

!  ---  inputs:
      integer, intent(in) :: nlay, nlp1

      real (kind=kind_phys), dimension(nlay,ngptsw), intent(in) ::      &
     &       taug, taur, cldfmc
      real (kind=kind_phys), dimension(nlay,nbdsw),  intent(in) ::      &
     &       taucw, ssacw, asycw, tauae, ssaae, asyae

      real (kind=kind_phys), dimension(ngptsw), intent(in) :: sfluxzen

      real (kind=kind_phys), dimension(2),  intent(in) :: albbm, albdf

      real (kind=kind_phys), intent(in) :: cosz, sntz, cf1, cf0, ssolar

!  ---  outputs:
      real (kind=kind_phys), dimension(nlp1,nbdsw), intent(out) ::      &
     &       fxupc, fxdnc, fxup0, fxdn0

      real (kind=kind_phys), dimension(2), intent(out) :: sfbmc, sfdfc, &
     &       sfbm0, sfdf0

      real (kind=kind_phys), intent(out) :: suvbfc, suvbf0, ftoadc,     &
     &       ftoauc, ftoau0, fsfcuc, fsfcu0, fsfcdc, fsfcd0

!  ---  locals:
      real (kind=kind_phys), dimension(nlay) :: ztaus, zssas, zasys,    &
     &       zldbt0

      real (kind=kind_phys), dimension(nlp1) :: zrefb, zrefd, ztrab,    &
     &       ztrad, ztdbt, zldbt, zfu, zfd

      real (kind=kind_phys) :: ztau1, zssa1, zasy1, ztau0, zssa0,       &
     &       zasy0, zasy3, zssaw, zasyw, zgam1, zgam2, zgam3, zgam4,    &
     &       za1, za2, zb1, zb2, zrk, zrk2, zrp, zrp1, zrm1, zrpp,      &
     &       zrkg1, zrkg3, zrkg4, zexp1, zexm1, zexp2, zexm2, zden1,    &
     &       zexp3, zexp4, ze1r45, ftind, zsolar, ztdbt0, zr1, zr2,     &
     &       zr3, zr4, zr5, zt1, zt2, zt3, zf1, zf2

      integer :: ib, ibd, jb, jg, k, kp, itind
!
!===> ...  begin here
!
!  --- ... initialization of output fluxes

      do ib = 1, nbdsw
        do k = 1, nlp1
          fxdnc(k,ib) = f_zero
          fxupc(k,ib) = f_zero
          fxdn0(k,ib) = f_zero
          fxup0(k,ib) = f_zero
        enddo
      enddo

      ftoadc = f_zero
      ftoauc = f_zero
      ftoau0 = f_zero
      fsfcuc = f_zero
      fsfcu0 = f_zero
      fsfcdc = f_zero
      fsfcd0 = f_zero

!! --- ...  uv-b surface downward fluxes
      suvbfc  = f_zero
      suvbf0  = f_zero

!! --- ...  output surface flux components
      sfbmc(1) = f_zero
      sfbmc(2) = f_zero
      sfdfc(1) = f_zero
      sfdfc(2) = f_zero
      sfbm0(1) = f_zero
      sfbm0(2) = f_zero
      sfdf0(1) = f_zero
      sfdf0(2) = f_zero

!  --- ...  loop over all g-points in each band

      lab_do_jg : do jg = 1, ngptsw

        jb = ngb(jg)
        ib = jb + 1 - nblow
        ibd = idxsfc(jb)         ! spectral band index

        zsolar = ssolar * sfluxzen(jg)

!  --- ...  set up toa direct beam and surface values (beam and diff)

        ztdbt(nlp1) = f_one
        ztdbt0   = f_one

        zldbt(1) = f_zero
        if (ibd /= 0) then
          zrefb(1) = albbm(ibd)
          zrefd(1) = albdf(ibd)
        else
          zrefb(1) = 0.5 * (albbm(1) + albbm(2))
          zrefd(1) = 0.5 * (albdf(1) + albdf(2))
        endif
        ztrab(1) = f_zero
        ztrad(1) = f_zero

!  --- ...  compute clear-sky optical parameters, layer reflectance and transmittance

        do k = nlay, 1, -1
          kp = k + 1

          ztau0 = max( ftiny, taur(k,jg)+taug(k,jg)+tauae(k,ib) )
          zssa0 = taur(k,jg) + tauae(k,ib)*ssaae(k,ib)
          zasy0 = asyae(k,ib)*ssaae(k,ib)*tauae(k,ib)
          zssaw = min( oneminus, zssa0 / ztau0 )
          zasyw = zasy0 / max( ftiny, zssa0 )

!  --- ...  saving clear-sky quantities for later total-sky usage
          ztaus(k) = ztau0
          zssas(k) = zssa0
          zasys(k) = zasy0

!  --- ...  delta scaling for clear-sky condition
          za1 = zasyw * zasyw
          za2 = zssaw * za1

          ztau1 = (f_one - za2) * ztau0
          zssa1 = (zssaw - za2) / (f_one - za2)
!org      zasy1 = (zasyw - za1) / (f_one - za1)   ! this line is replaced by the next
          zasy1 = zasyw / (f_one + zasyw)         ! to reduce truncation error
          zasy3 = 0.75 * zasy1

!  --- ...  general two-stream expressions
          if ( iswmode == 1 ) then
            zgam1 = 1.75 - zssa1 * (f_one + zasy3)
            zgam2 =-0.25 + zssa1 * (f_one - zasy3)
            zgam3 = 0.5  - zasy3 * cosz
          elseif ( iswmode == 2 ) then               ! pifm
            zgam1 = 2.0 - zssa1 * (1.25 + zasy3)
            zgam2 = 0.75* zssa1 * (f_one- zasy1)
            zgam3 = 0.5 - zasy3 * cosz
          elseif ( iswmode == 3 ) then               ! discrete ordinates
            zgam1 = zsr3 * (2.0 - zssa1 * (1.0 + zasy1)) * 0.5
            zgam2 = zsr3 * zssa1 * (1.0 - zasy1) * 0.5
            zgam3 = (1.0 - zsr3 * zasy1 * cosz) * 0.5
          endif
          zgam4 = f_one - zgam3

!  --- ...  compute homogeneous reflectance and transmittance

          if ( zssaw >= zcrit ) then    ! for conservative scattering
            za1 = zgam1 * cosz - zgam3
            za2 = zgam1 * ztau1

!  --- ...  use exponential lookup table for transmittance, or expansion
!           of exponential for low optical depth

            zb1 = min ( ztau1*sntz , 500.0 )
            if ( zb1 <= od_lo ) then
              zb2 = f_one - zb1 + 0.5*zb1*zb1
            else
              ftind = zb1 / (bpade + zb1)
              itind = ftind*ntbmx + 0.5
              zb2 = exp_tbl(itind)
            endif

!      ...  collimated beam
            zrefb(kp) = max(f_zero, min(f_one,                          &
     &                  (za2 - za1*(f_one - zb2))/(f_one + za2) ))
            ztrab(kp) = max(f_zero, min(f_one, f_one-zrefb(kp) ))

!      ...  isotropic incidence
            zrefd(kp) = max(f_zero, min(f_one, za2/(f_one + za2) ))
            ztrad(kp) = max(f_zero, min(f_one, f_one-zrefd(kp) ))

          else                          ! for non-conservative scattering
            za1 = zgam1*zgam4 + zgam2*zgam3
            za2 = zgam1*zgam3 + zgam2*zgam4
            zrk = sqrt ( (zgam1 - zgam2) * (zgam1 + zgam2) )
            zrk2= 2.0 * zrk

            zrp  = zrk * cosz
            zrp1 = f_one + zrp
            zrm1 = f_one - zrp
            zrpp = f_one - zrp*zrp
            zrkg1= zrk + zgam1
            zrkg3= zrk * zgam3
            zrkg4= zrk * zgam4

            zr1  = zrm1 * (za2 + zrkg3)
            zr2  = zrp1 * (za2 - zrkg3)
            zr3  = zrk2 * (zgam3 - za2*cosz)
            zr4  = zrpp * zrkg1
            zr5  = zrpp * (zrk - zgam1)

            zt1  = zrp1 * (za1 + zrkg4)
            zt2  = zrm1 * (za1 - zrkg4)
            zt3  = zrk2 * (zgam4 + za1*cosz)

!  --- ...  use exponential lookup table for transmittance, or expansion
!           of exponential for low optical depth

            zb1 = min ( zrk*ztau1, 500.0 )
            if ( zb1 <= od_lo ) then
              zexm1 = f_one - zb1 + 0.5*zb1*zb1
            else
              ftind = zb1 / (bpade + zb1)
              itind = ftind*ntbmx + 0.5
              zexm1 = exp_tbl(itind)
            endif
            zexp1 = f_one / zexm1

            zb2 = min ( sntz*ztau1, 500.0 )
            if ( zb2 <= od_lo ) then
              zexm2 = f_one - zb2 + 0.5*zb2*zb2
            else
              ftind = zb2 / (bpade + zb2)
              itind = ftind*ntbmx + 0.5
              zexm2 = exp_tbl(itind)
            endif
            zexp2 = f_one / zexm2
            ze1r45 = zr4*zexp1 + zr5*zexm1

!      ...  collimated beam
            if (ze1r45>=-eps1 .and. ze1r45<=eps1) then
              zrefb(kp) = eps1
              ztrab(kp) = zexm2
            else
              zden1 = zssa1 / ze1r45
              zrefb(kp) = max(f_zero, min(f_one,                        &
     &                    (zr1*zexp1 - zr2*zexm1 - zr3*zexm2)*zden1 ))
              ztrab(kp) = max(f_zero, min(f_one, zexm2*(f_one           &
     &                  - (zt1*zexp1 - zt2*zexm1 - zt3*zexp2)*zden1) ))
            endif

!      ...  diffuse beam
            zden1 = zr4 / (ze1r45 * zrkg1)
            zrefd(kp) = max(f_zero, min(f_one,                          &
     &                  zgam2*(zexp1 - zexm1)*zden1 ))
            ztrad(kp) = max(f_zero, min(f_one, zrk2*zden1 ))
          endif    ! end if_zssaw_block

!  --- ...  direct beam transmittance. use exponential lookup table
!           for transmittance, or expansion of exponential for low
!           optical depth

          zr1 = ztau1 * sntz
          if ( zr1 <= od_lo ) then
            zexp3 = f_one - zr1 + 0.5*zr1*zr1
          else
            ftind = zr1 / (bpade + zr1)
            itind = max(0, min(ntbmx, int(0.5+ntbmx*ftind) ))
            zexp3 = exp_tbl(itind)
          endif

          ztdbt(k)  = zexp3 * ztdbt(kp)
          zldbt(kp) = zexp3

!  --- ...  pre-delta-scaling clear and cloudy direct beam transmittance
!           (must use 'orig', unscaled cloud optical depth)

          zr1 = ztau0 * sntz
          if ( zr1 <= od_lo ) then
            zexp4 = f_one - zr1 + 0.5*zr1*zr1
          else
            ftind = zr1 / (bpade + zr1)
            itind = max(0, min(ntbmx, int(0.5+ntbmx*ftind) ))
            zexp4 = exp_tbl(itind)
          endif

          zldbt0(k) = zexp4
          ztdbt0 = zexp4 * ztdbt0
        enddo    ! end do_k_loop

        call swflux                                                     &
!  ---  inputs:
     &     ( zrefb,zrefd,ztrab,ztrad,zldbt,ztdbt,                       &
     &       nlay, nlp1,                                                &
!  ---  outputs:
     &       zfu, zfd                                                   &
     &     )

!  --- ...  compute upward and downward fluxes at levels
        do k = 1, nlp1
          fxup0(k,ib) = fxup0(k,ib) + zsolar*zfu(k)
          fxdn0(k,ib) = fxdn0(k,ib) + zsolar*zfd(k)
        enddo

!! --- ...  surface downward beam/diffuse flux components
        zb1 = zsolar*ztdbt0
        zb2 = zsolar*(zfd(1) - ztdbt0)

        if (ibd /= 0) then
          sfbm0(ibd) = sfbm0(ibd) + zb1
          sfdf0(ibd) = sfdf0(ibd) + zb2
        else
          zf1 = 0.5 * zb1
          zf2 = 0.5 * zb2
          sfbm0(1) = sfbm0(1) + zf1
          sfdf0(1) = sfdf0(1) + zf2
          sfbm0(2) = sfbm0(2) + zf1
          sfdf0(2) = sfdf0(2) + zf2
        endif
!       sfbm0(ibd) = sfbm0(ibd) + zsolar*ztdbt0
!       sfdf0(ibd) = sfdf0(ibd) + zsolar*(zfd(1) - ztdbt0)

!  --- ...  compute total sky optical parameters, layer reflectance and transmittance

        if ( cf1 > eps ) then

!  --- ...  set up toa direct beam and surface values (beam and diff)
          ztdbt0 = f_one
          zldbt(1) = f_zero

          do k = nlay, 1, -1
            kp = k + 1
            if ( cldfmc(k,jg) > ftiny ) then      ! it is a cloudy-layer

              ztau0 = ztaus(k) + taucw(k,ib)
              zssa0 = zssas(k) + ssacw(k,ib)
              zasy0 = zasys(k) + asycw(k,ib)
              zssaw = min(oneminus, zssa0 / ztau0)
              zasyw = zasy0 / max(ftiny, zssa0)

!  --- ...  delta scaling for total-sky condition
              za1 = zasyw * zasyw
              za2 = zssaw * za1

              ztau1 = (f_one - za2) * ztau0
              zssa1 = (zssaw - za2) / (f_one - za2)
!org          zasy1 = (zasyw - za1) / (f_one - za1)
              zasy1 = zasyw / (f_one + zasyw)
              zasy3 = 0.75 * zasy1

!  --- ...  general two-stream expressions
              if ( iswmode == 1 ) then
                zgam1 = 1.75 - zssa1 * (f_one + zasy3)
                zgam2 =-0.25 + zssa1 * (f_one - zasy3)
                zgam3 = 0.5  - zasy3 * cosz
              elseif ( iswmode == 2 ) then               ! pifm
                zgam1 = 2.0 - zssa1 * (1.25 + zasy3)
                zgam2 = 0.75* zssa1 * (f_one- zasy1)
                zgam3 = 0.5 - zasy3 * cosz
              elseif ( iswmode == 3 ) then               ! discrete ordinates
                zgam1 = zsr3 * (2.0 - zssa1 * (1.0 + zasy1)) * 0.5
                zgam2 = zsr3 * zssa1 * (1.0 - zasy1) * 0.5
                zgam3 = (1.0 - zsr3 * zasy1 * cosz) * 0.5
              endif
              zgam4 = f_one - zgam3

!  --- ...  compute homogeneous reflectance and transmittance

              if ( zssaw >= zcrit ) then    ! for conservative scattering
                za1 = zgam1 * cosz - zgam3
                za2 = zgam1 * ztau1

!  --- ...  use exponential lookup table for transmittance, or expansion
!           of exponential for low optical depth

                zb1 = min ( ztau1*sntz , 500.0 )
                if ( zb1 <= od_lo ) then
                  zb2 = f_one - zb1 + 0.5*zb1*zb1
                else
                  ftind = zb1 / (bpade + zb1)
                  itind = ftind*ntbmx + 0.5
                  zb2 = exp_tbl(itind)
                endif

!      ...  collimated beam
                zrefb(kp) = max(f_zero, min(f_one,                      &
     &                      (za2 - za1*(f_one - zb2))/(f_one + za2) ))
                ztrab(kp) = max(f_zero, min(f_one, f_one-zrefb(kp)))

!      ...  isotropic incidence
                zrefd(kp) = max(f_zero, min(f_one, za2 / (f_one+za2) ))
                ztrad(kp) = max(f_zero, min(f_one, f_one - zrefd(kp) ))

              else                          ! for non-conservative scattering
                za1 = zgam1*zgam4 + zgam2*zgam3
                za2 = zgam1*zgam3 + zgam2*zgam4
                zrk = sqrt ( (zgam1 - zgam2) * (zgam1 + zgam2) )
                zrk2= 2.0 * zrk

                zrp  = zrk * cosz
                zrp1 = f_one + zrp
                zrm1 = f_one - zrp
                zrpp = f_one - zrp*zrp
                zrkg1= zrk + zgam1
                zrkg3= zrk * zgam3
                zrkg4= zrk * zgam4

                zr1  = zrm1 * (za2 + zrkg3)
                zr2  = zrp1 * (za2 - zrkg3)
                zr3  = zrk2 * (zgam3 - za2*cosz)
                zr4  = zrpp * zrkg1
                zr5  = zrpp * (zrk - zgam1)

                zt1  = zrp1 * (za1 + zrkg4)
                zt2  = zrm1 * (za1 - zrkg4)
                zt3  = zrk2 * (zgam4 + za1*cosz)

!  --- ...  use exponential lookup table for transmittance, or expansion
!           of exponential for low optical depth

                zb1 = min ( zrk*ztau1, 500.0 )
                if ( zb1 <= od_lo ) then
                  zexm1 = f_one - zb1 + 0.5*zb1*zb1
                else
                  ftind = zb1 / (bpade + zb1)
                  itind = ftind*ntbmx + 0.5
                  zexm1 = exp_tbl(itind)
                endif
                zexp1 = f_one / zexm1

                zb2 = min ( ztau1*sntz, 500.0 )
                if ( zb2 <= od_lo ) then
                  zexm2 = f_one - zb2 + 0.5*zb2*zb2
                else
                  ftind = zb2 / (bpade + zb2)
                  itind = ftind*ntbmx + 0.5
                  zexm2 = exp_tbl(itind)
                endif
                zexp2 = f_one / zexm2
                ze1r45 = zr4*zexp1 + zr5*zexm1

!      ...  collimated beam
                if ( ze1r45>=-eps1 .and. ze1r45<=eps1 ) then
                  zrefb(kp) = eps1
                  ztrab(kp) = zexm2
                else
                  zden1 = zssa1 / ze1r45
                  zrefb(kp) = max(f_zero, min(f_one,                    &
     &                        (zr1*zexp1-zr2*zexm1-zr3*zexm2)*zden1 ))
                  ztrab(kp) = max(f_zero, min(f_one, zexm2*(f_one -     &
     &                        (zt1*zexp1-zt2*zexm1-zt3*zexp2)*zden1) ))
                endif

!      ...  diffuse beam
                zden1 = zr4 / (ze1r45 * zrkg1)
                zrefd(kp) = max(f_zero, min(f_one,                      &
     &                      zgam2*(zexp1 - zexm1)*zden1 ))
                ztrad(kp) = max(f_zero, min(f_one, zrk2*zden1 ))
              endif    ! end if_zssaw_block

!  --- ...  direct beam transmittance. use exponential lookup table
!           for transmittance, or expansion of exponential for low
!           optical depth

              zr1 = ztau1 * sntz
              if ( zr1 <= od_lo ) then
                zexp3 = f_one - zr1 + 0.5*zr1*zr1
              else
                ftind = zr1 / (bpade + zr1)
                itind = max(0, min(ntbmx, int(0.5+ntbmx*ftind) ))
                zexp3 = exp_tbl(itind)
              endif

              zldbt(kp) = zexp3
              ztdbt(k)  = zexp3 * ztdbt(kp)

!  --- ...  pre-delta-scaling clear and cloudy direct beam transmittance
!           (must use 'orig', unscaled cloud optical depth)

              zr1 = ztau0 * sntz
              if ( zr1 <= od_lo ) then
                zexp4 = f_one - zr1 + 0.5*zr1*zr1
              else
                ftind = zr1 / (bpade + zr1)
                itind = max(0, min(ntbmx, int(0.5+ntbmx*ftind) ))
                zexp4 = exp_tbl(itind)
              endif

              ztdbt0 = zexp4 * ztdbt0

            else     ! if_cldfmc_block  ---  it is a clear layer

!  --- ...  direct beam transmittance
              ztdbt(k) = zldbt(kp) * ztdbt(kp)

!  --- ...  pre-delta-scaling clear and cloudy direct beam transmittance
              ztdbt0 = zldbt0(k) * ztdbt0

            endif    ! end if_cldfmc_block
          enddo   ! end do_k_loop

!  --- ...  perform vertical quadrature

          call swflux                                                   &
!  ---  inputs:
     &     ( zrefb,zrefd,ztrab,ztrad,zldbt,ztdbt,                       &
     &       nlay, nlp1,                                                &
!  ---  outputs:
     &       zfu, zfd                                                   &
     &     )

!  --- ...  compute upward and downward fluxes at levels
          do k = 1, nlp1
            fxupc(k,ib) = fxupc(k,ib) + zsolar*zfu(k)
            fxdnc(k,ib) = fxdnc(k,ib) + zsolar*zfd(k)
          enddo

!! --- ...  surface downward beam/diffused flux components
          zb1 = zsolar*ztdbt0
          zb2 = zsolar*(zfd(1) - ztdbt0)

          if (ibd /= 0) then
           sfbmc(ibd) = sfbmc(ibd) + zb1
           sfdfc(ibd) = sfdfc(ibd) + zb2
          else
            zf1 = 0.5 * zb1
            zf2 = 0.5 * zb2
            sfbmc(1) = sfbmc(1) + zf1
            sfdfc(1) = sfdfc(1) + zf2
            sfbmc(2) = sfbmc(2) + zf1
            sfdfc(2) = sfdfc(2) + zf2
          endif
!         sfbmc(ibd) = sfbmc(ibd) + zsolar*ztdbt0
!         sfdfc(ibd) = sfdfc(ibd) + zsolar*(zfd(1) - ztdbt0)

        endif      ! end if_cf1_block

      enddo  lab_do_jg

!  --- ...  end of g-point loop

      do ib = 1, nbdsw
        ftoadc = ftoadc + fxdn0(nlp1,ib)
        ftoau0 = ftoau0 + fxup0(nlp1,ib)
        fsfcu0 = fsfcu0 + fxup0(1,ib)
        fsfcd0 = fsfcd0 + fxdn0(1,ib)
      enddo

!! --- ...  uv-b surface downward flux
      ibd = nuvb - nblow + 1
      suvbf0 = fxdn0(1,ibd)

      if ( cf1 <= eps ) then       ! clear column, set total-sky=clear-sky fluxes
        do ib = 1, nbdsw
          do k = 1, nlp1
            fxupc(k,ib) = fxup0(k,ib)
            fxdnc(k,ib) = fxdn0(k,ib)
          enddo
        enddo

        ftoauc = ftoau0
        fsfcuc = fsfcu0
        fsfcdc = fsfcd0

!! --- ...  surface downward beam/diffused flux components
        sfbmc(1) = sfbm0(1)
        sfdfc(1) = sfdf0(1)
        sfbmc(2) = sfbm0(2)
        sfdfc(2) = sfdf0(2)

!! --- ...  uv-b surface downward flux
        suvbfc = suvbf0
      else                        ! cloudy column, compute total-sky fluxes
        do ib = 1, nbdsw
          ftoauc = ftoauc + fxupc(nlp1,ib)
          fsfcuc = fsfcuc + fxupc(1,ib)
          fsfcdc = fsfcdc + fxdnc(1,ib)
        enddo

!! --- ...  uv-b surface downward flux
        suvbfc = fxdnc(1,ibd)
      endif    ! end if_cf1_block

      return
!...................................
      end subroutine spcvrtm
!-----------------------------------


!-----------------------------------
      subroutine swflux                                                 &
!...................................
!  ---  inputs:
     &     ( zrefb,zrefd,ztrab,ztrad,zldbt,ztdbt,                       &
     &       nlay, nlp1,                                                &
!  ---  outputs:
     &       zfu, zfd                                                   &
     &     )

!  ===================  program usage description  ===================  !
!                                                                       !
!   purpose:  computes the upward and downward radiation fluxes         !
!                                                                       !
!   interface:  "swflux" is called by "spcvrc" and "spcvrm"             !
!                                                                       !
!   subroutines called : none                                           !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  input variables:                                                     !
!    zrefb(nlp1)     - layer direct beam reflectivity                   !
!    zrefd(nlp1)     - layer diffuse reflectivity                       !
!    ztrab(nlp1)     - layer direct beam transmissivity                 !
!    ztrad(nlp1)     - layer diffuse transmissivity                     !
!    zldbt(nlp1)     - layer mean beam transmittance                    !
!    ztdbt(nlp1)     - total beam transmittance at levels               !
!    nlay, nlp1      - number of layers/levels                          !
!                                                                       !
!  output variables:                                                    !
!    zfu  (nlp1)     - upward flux at layer interface                   !
!    zfd  (nlp1)     - downward flux at layer interface                 !
!                                                                       !
!  *******************************************************************  !
!  ======================  end of description block  =================  !

!  ---  inputs:
      integer, intent(in) :: nlay, nlp1

      real (kind=kind_phys), dimension(nlp1), intent(in) :: zrefb,      &
     &       zrefd, ztrab, ztrad, ztdbt, zldbt

!  ---  outputs:
      real (kind=kind_phys), dimension(nlp1), intent(out) :: zfu, zfd

!  ---  locals:
      real (kind=kind_phys), dimension(nlp1) :: zrupb,zrupd,zrdnd,ztdn

      real (kind=kind_phys) :: zden1

      integer :: k, kp
!
!===> ... begin here
!

!  --- ...  link lowest layer with surface

        zrupb(1) = zrefb(1)        ! direct beam
        zrupd(1) = zrefd(1)        ! diffused

!  --- ...  pass from bottom to top

        do k = 1, nlay
          kp = k + 1

          zden1 = f_one / ( f_one - zrupd(k)*zrefd(kp) )
          zrupb(kp) = zrefb(kp) + ( ztrad(kp) *                         &
     &                ( (ztrab(kp) - zldbt(kp))*zrupd(k) +              &
     &                zldbt(kp)*zrupb(k)) ) * zden1
          zrupd(kp) = zrefd(kp) + ztrad(kp)*ztrad(kp)*zrupd(k)*zden1
        enddo

!  --- ...  upper boundary conditions

        ztdn (nlp1) = f_one
        zrdnd(nlp1) = f_zero
        ztdn (nlay) = ztrab(nlp1)
        zrdnd(nlay) = zrefd(nlp1)

!  --- ...  pass from top to bottom

        do k = nlay, 2, -1
          zden1 = f_one / (f_one - zrefd(k)*zrdnd(k))
          ztdn (k-1) = ztdbt(k)*ztrab(k) + ( ztrad(k) *                 &
     &                 ( (ztdn(k) - ztdbt(k)) + ztdbt(k) *              &
     &                 zrefb(k)*zrdnd(k) )) * zden1
          zrdnd(k-1) = zrefd(k) + ztrad(k)*ztrad(k)*zrdnd(k)*zden1
        enddo

!  --- ...  up and down-welling fluxes at levels

        do k = 1, nlp1
          zden1 = f_one / (f_one - zrdnd(k)*zrupd(k))
          zfu(k) = ( ztdbt(k)*zrupb(k) +                                &
     &             (ztdn(k) - ztdbt(k))*zrupd(k) ) * zden1
          zfd(k) = ztdbt(k) + ( ztdn(k) - ztdbt(k) +                    &
     &             ztdbt(k)*zrupb(k)*zrdnd(k) ) * zden1
        enddo

      return
!...................................
      end subroutine swflux
!-----------------------------------


!-----------------------------------
      subroutine taumol                                                 &
!...................................
!  ---  inputs:
     &     ( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, map_nday_ipt, map_nday_jj, &
             nlay, ix, nday, async_id, fulljj, blockjj, jb, ng00, np, &
             offset_ng, nday_length, max_nday_length, jjoffset,    &
!  ---  outputs:
     &       small_sfluxzen, small_taug, small_taur                                       &
     &     )

!  ==================   program usage description   ==================  !
!                                                                       !
!  description:                                                         !
!    calculate optical depths for gaseous absorption and rayleigh       !
!    scattering.                                                        !
!                                                                       !
!  subroutines called: taugb## (## = 16 - 29)                           !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                         size !
!    colamt  - real, column amounts of absorbing gases the index        !
!                    are for h2o, co2, o3, n2o, ch4, and o2,            !
!                    respectively (molecules/cm**2)          nlay*maxgas!
!    colmol  - real, total column amount (dry air+water vapor)     nlay !
!    facij   - real, for each layer, these are factors that are         !
!                    needed to compute the interpolation factors        !
!                    that multiply the appropriate reference k-         !
!                    values.  a value of 0/1 for i,j indicates          !
!                    that the corresponding factor multiplies           !
!                    reference k-value for the lower/higher of the      !
!                    two appropriate temperatures, and altitudes,       !
!                    respectively.                                 naly !
!    jp      - real, the index of the lower (in altitude) of the        !
!                    two appropriate ref pressure levels needed         !
!                    for interpolation.                            nlay !
!    jt, jt1 - integer, the indices of the lower of the two approp      !
!                    ref temperatures needed for interpolation (for     !
!                    pressure levels jp and jp+1, respectively)    nlay !
!    laytrop - integer, tropopause layer index                       1  !
!    forfac  - real, scale factor needed to foreign-continuum.     nlay !
!    forfrac - real, factor needed for temperature interpolation   nlay !
!    indfor  - integer, index of the lower of the two appropriate       !
!                    reference temperatures needed for foreign-         !
!                    continuum interpolation                       nlay !
!    selffac - real, scale factor needed to h2o self-continuum.    nlay !
!    selffrac- real, factor needed for temperature interpolation        !
!                    of reference h2o self-continuum data          nlay !
!    indself - integer, index of the lower of the two appropriate       !
!                    reference temperatures needed for the self-        !
!                    continuum interpolation                       nlay !
!    nlay    - integer, number of vertical layers                    1  !
!                                                                       !
!  output:                                                              !
!    sfluxzen- real, spectral distribution of incoming solar flux ngptsw!
!    taug    - real, spectral optical depth for gases        nlay*ngptsw!
!    taur    - real, opt depth for rayleigh scattering       nlay*ngptsw!
!                                                                       !
!  ===================================================================  !
!  ************     original subprogram description    ***************  !
!                                                                       !
!                  optical depths developed for the                     !
!                                                                       !
!                rapid radiative transfer model (rrtm)                  !
!                                                                       !
!            atmospheric and environmental research, inc.               !
!                        131 hartwell avenue                            !
!                        lexington, ma 02421                            !
!                                                                       !
!                                                                       !
!                           eli j. mlawer                               !
!                         jennifer delamere                             !
!                         steven j. taubman                             !
!                         shepard a. clough                             !
!                                                                       !
!                                                                       !
!                                                                       !
!                       email:  mlawer@aer.com                          !
!                       email:  jdelamer@aer.com                        !
!                                                                       !
!        the authors wish to acknowledge the contributions of the       !
!        following people:  patrick d. brown, michael j. iacono,        !
!        ronald e. farren, luke chen, robert bergstrom.                 !
!                                                                       !
!  *******************************************************************  !
!                                                                       !
!  taumol                                                               !
!                                                                       !
!    this file contains the subroutines taugbn (where n goes from       !
!    16 to 29).  taugbn calculates the optical depths and planck        !
!    fractions per g-value and layer for band n.                        !
!                                                                       !
!  output:  optical depths (unitless)                                   !
!           fractions needed to compute planck functions at every layer !
!           and g-value                                                 !
!                                                                       !
!  modifications:                                                       !
!                                                                       !
! revised: adapted to f90 coding, j.-j.morcrette, ecmwf, feb 2003       !
! revised: modified for g-point reduction, mjiacono, aer, dec 2003      !
! revised: reformatted for consistency with rrtmg_lw, mjiacono, aer,    !
!          jul 2006                                                     !
!                                                                       !
!  *******************************************************************  !
!  ======================  end of description block  =================  !

!  ---  inputs:
      integer, intent(in) :: nlay, laytrop(max_nday_length), ix, nday(fulljj), fulljj, &
         blockjj, jb, ng00, np, offset_ng, nday_length, jjoffset, max_nday_length

      integer, dimension(nday_length) :: map_nday_ipt, map_nday_jj
      integer, dimension(max_nday_length, nlay), intent(in) :: indfor, indself,          &
     &       jp, jt, jt1

      real (kind=kind_phys), dimension(max_nday_length, nlay),  intent(in) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac

      real (kind=kind_phys), dimension(max_nday_length, nlay, maxgas),intent(in) :: colamt

!  ---  outputs:
      real (kind=kind_phys), dimension(nday_length, nlay, ng00) :: small_taur, small_taug
      real (kind=kind_phys), dimension(nday_length, ng00) :: small_sfluxzen

!  ---  locals:
      !integer, parameter :: small_factor = 10
      !integer, parameter :: small_ix = int(ix/small_factor)
      real (kind=kind_phys) :: fsa, speccomb, specmult, colm1, colm2

      integer, dimension(nday_length, nlay,nblow:nbhgh) :: id0, id1

      integer :: ibd, j, jsa, k, klow, khgh, klim, ks, njb, ns, ipt, jj, i2, i3, n
      integer :: async_id
!
!===> ... begin here
!
!  --- ...  loop over each spectral band
      !do i2 = small_factor
      !  --- ...  indices for layer optical depth
      !$acc data create(id0, id1) async(async_id)
      !$acc parallel loop collapse(2) async(async_id)
      do k = 1, nlay
         do n = 1, nday_length
            if (k .le. laytrop(n)) then
               id0(n, k,jb) = ((jp(n, k)-1)*5 + (jt (n, k)-1)) * nspa(jb)
               id1(n, k,jb) = ( jp(n, k)   *5 + (jt1(n, k)-1)) * nspa(jb)
            else
               id0(n, k,jb) = ((jp(n, k)-13)*5 + (jt (n, k)-1)) * nspb(jb)
               id1(n, k,jb) = ((jp(n, k)-12)*5 + (jt1(n, k)-1)) * nspb(jb)
            end if
         enddo
      end do

      !  --- ...  calculate spectral flux at toa

      !$acc parallel loop private(ks, colm1, colm2, speccomb, &
      !$acc&         specmult, jsa, fsa) async(async_id)
      do n = 1, nday_length
         ibd = ibx(jb)
         njb = ng (jb)
         ns  = ngs(jb)
         if ((jb .eq. 16) .or. (jb .eq. 20) .or. (jb .eq. 23) .or. &
               (jb .eq. 25) .or. (jb .eq. 26) .or. (jb .eq. 29)) then
            !$acc loop seq
            do j = 1, njb
               !sfluxzen(ipt, ngs(jb)+j, jj) = sfluxref01(j,1,ibx(jb))
               small_sfluxzen(n, j + offset_ng) = sfluxref01(j,1,ibx(jb))
            enddo

         elseif (jb .eq. 27) then
            !$acc loop seq
            do j = 1, njb
               !sfluxzen(ipt, ngs(jb)+j, jj) = scalekur * sfluxref01(j,1,ibx(jb))
               small_sfluxzen(n, j + offset_ng) = scalekur * sfluxref01(j,1,ibx(jb))
            enddo

         elseif ((jb .eq. 17) .or. (jb .eq. 28)) then
            ks = nlay
            !$acc loop seq
            do k = laytrop(n), nlay-1 ! lab_do_k1
               if (jp(n, k)<layreffr(jb) .and. jp(n, k+1)>=layreffr(jb)) then
                  ks = k + 1
                  exit ! lab_do_k1
               endif
            enddo  ! lab_do_k1

            colm1 = colamt(n, ks,ix1(jb))
            colm2 = colamt(n, ks,ix2(jb))
            speccomb = colm1 + strrat(jb)*colm2
            specmult = specwt(jb) * min( oneminus, colm1/speccomb )
            jsa = 1 + int( specmult )
            fsa = mod(specmult, f_one)
            !$acc loop seq
            do j = 1, njb
               !sfluxzen(ipt, ngs(jb)+j, jj) = sfluxref02(j,jsa,ibd)                   &
               !&           + fsa * (sfluxref02(j,jsa+1,ibd) - sfluxref02(j,jsa,ibd))
               small_sfluxzen(n, j + offset_ng) = sfluxref02(j,jsa,ibd)                   &
               &           + fsa * (sfluxref02(j,jsa+1,ibd) - sfluxref02(j,jsa,ibd))
            enddo

         else
            ks = laytrop(n)
            !$acc loop seq
            do k = 1, laytrop(n)-1 ! lab_do_k2
               if (jp(n, k)<layreffr(jb) .and. jp(n, k+1)>=layreffr(jb)) then
                  ks = k + 1
                  exit ! lab_do_k2
               endif
            enddo  ! lab_do_k2

            colm1 = colamt(n, ks,ix1(jb))
            colm2 = colamt(n, ks,ix2(jb))
            speccomb = colm1 + strrat(jb)*colm2
            specmult = specwt(jb) * min( oneminus, colm1/speccomb )
            jsa = 1 + int( specmult )
            fsa = mod(specmult, f_one)
            !$acc loop seq
            do j = 1, njb
               !sfluxzen(ipt, ngs(jb)+j, jj) = sfluxref03(j,jsa,ibd)                   &
               !&           + fsa * (sfluxref03(j,jsa+1,ibd) - sfluxref03(j,jsa,ibd))
               small_sfluxzen(n, j + offset_ng) = sfluxref03(j,jsa,ibd)                   &
               &           + fsa * (sfluxref03(j,jsa+1,ibd) - sfluxref03(j,jsa,ibd))
            enddo
         end if
      enddo

      select case (jb)

         !  --- ...  call taumol## to calculate layer optical depth
      case(16)
         call taumol16( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
      &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
      &       async_id, fulljj, id0, id1, blockjj, ng00, &
      &       map_nday_ipt, map_nday_jj, nday_length, max_nday_length, jjoffset, &
      &       small_sfluxzen, small_taug, small_taur, np, offset_ng)
     
      case(17)
         call taumol17( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
      &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
      &       async_id, fulljj, id0, id1, blockjj, ng00, &
      &       map_nday_ipt, map_nday_jj, nday_length, max_nday_length, jjoffset, &
      &       small_sfluxzen, small_taug, small_taur, np, offset_ng)
      
      case(18)
      call taumol18( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
      &       async_id, fulljj, id0, id1, blockjj, ng00, &
      &       map_nday_ipt, map_nday_jj, nday_length, max_nday_length, jjoffset, &
      &       small_sfluxzen, small_taug, small_taur, np, offset_ng)
     
      case(19)
      call taumol19( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
      &       async_id, fulljj, id0, id1, blockjj, ng00, &
      &       map_nday_ipt, map_nday_jj, nday_length, max_nday_length, jjoffset, &
      &       small_sfluxzen, small_taug, small_taur, np, offset_ng)
     
      case(20)
         call taumol20( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
      &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
      &       async_id, fulljj, id0, id1, blockjj, ng00, &
      &       map_nday_ipt, map_nday_jj, nday_length, max_nday_length, jjoffset, &
      &       small_sfluxzen, small_taug, small_taur, np, offset_ng)
      
      case(21)
         call taumol21( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
      &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
      &       async_id, fulljj, id0, id1, blockjj, ng00, &
      &       map_nday_ipt, map_nday_jj, nday_length, max_nday_length, jjoffset, &
      &       small_sfluxzen, small_taug, small_taur, np, offset_ng)
      
      case(22)
         call taumol22( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
      &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
      &       async_id, fulljj, id0, id1, blockjj, ng00, &
      &       map_nday_ipt, map_nday_jj, nday_length, max_nday_length, jjoffset, &
      &       small_sfluxzen, small_taug, small_taur, np, offset_ng)
      
      case(23)
         call taumol23( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
      &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
      &       async_id, fulljj, id0, id1, blockjj, ng00, &
      &       map_nday_ipt, map_nday_jj, nday_length, max_nday_length, jjoffset, &
      &       small_sfluxzen, small_taug, small_taur, np, offset_ng)
      
      case(24)
         call taumol24( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
      &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
      &       async_id, fulljj, id0, id1, blockjj, ng00, &
      &       map_nday_ipt, map_nday_jj, nday_length, max_nday_length, jjoffset, &
      &       small_sfluxzen, small_taug, small_taur, np, offset_ng)
      
      case(25)
         call taumol25( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
      &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
      &       async_id, fulljj, id0, id1, blockjj, ng00, &
      &       map_nday_ipt, map_nday_jj, nday_length, max_nday_length, jjoffset, &
      &       small_sfluxzen, small_taug, small_taur, np, offset_ng)
      
      case(26)
         call taumol26( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
      &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
      &       async_id, fulljj, id0, id1, blockjj, ng00, &
      &       map_nday_ipt, map_nday_jj, nday_length, max_nday_length, jjoffset, &
      &       small_sfluxzen, small_taug, small_taur, np, offset_ng)
      
      case(27)
         call taumol27( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
      &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
      &       async_id, fulljj, id0, id1, blockjj, ng00, &
      &       map_nday_ipt, map_nday_jj, nday_length, max_nday_length, jjoffset, &
      &       small_sfluxzen, small_taug, small_taur, np, offset_ng)
      
      case(28)
         call taumol28( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
      &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
      &       async_id, fulljj, id0, id1, blockjj, ng00, &
      &       map_nday_ipt, map_nday_jj, nday_length, max_nday_length, jjoffset, &
      &       small_sfluxzen, small_taug, small_taur, np, offset_ng)

      case(29)
      call taumol29( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
      &       async_id, fulljj, id0, id1, blockjj, ng00, &
      &       map_nday_ipt, map_nday_jj, nday_length, max_nday_length, jjoffset, &
      &       small_sfluxzen, small_taug, small_taur, np, offset_ng)

      end select
      !$acc end data

!...................................
      end subroutine taumol
!-----------------------------------

!-----------------------------------
      subroutine taumol16                                                 &
!...................................
     &     ( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
      &       async_id, fulljj, id0, id1, blockjj, ng00, &
      &       map_nday_ipt, map_nday_jj, nday_length, max_nday_length, jjoffset, &
      &       small_sfluxzen, small_taug, small_taur, np, offset_ng)

!  ------------------------------------------------------------------  !
!     band 16:  2600-3250 cm-1 (low - h2o,ch4; high - ch4)             !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb16
      
      integer, dimension(nday_length) :: map_nday_ipt, map_nday_jj
      integer :: ng00, np, offset_ng, nday_length, jjoffset, max_nday_length
      real (kind=kind_phys), dimension(nday_length, nlay, ng00) :: small_taur, small_taug
      real (kind=kind_phys), dimension(nday_length, ng00) :: small_sfluxzen
      integer :: nlay, laytrop(max_nday_length), ix, nday(fulljj), blockjj
      integer, dimension(max_nday_length, nlay) :: indfor, indself,          &
     &       jp, jt, jt1
      real (kind=kind_phys), dimension(max_nday_length, nlay) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac
      real (kind=kind_phys), dimension(max_nday_length, nlay, maxgas) :: colamt
      integer, dimension(nday_length, nlay,nblow:nbhgh) :: id0, id1
!  ---  locals:
      integer :: ipt, jj, async_id, fulljj

      real (kind=kind_phys) :: speccomb, specmult, tauray, fs, fs1,     &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111

      integer :: ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14
      integer :: inds, indf, indsp, indfp, j, js, k, n

!
!===> ... begin here
!

!  --- ... compute the optical depth by interpolating in ln(pressure),
!          temperature, and appropriate species.  below laytrop, the water
!          vapor self-continuum is interpolated (in temperature) separately.
      !$acc parallel loop collapse(2) private(speccomb, specmult, js, fs, fs1, fac000, &
      !$acc&     fac010, fac100, fac110, fac001, fac011, fac101, fac111, &
      !$acc&     ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14, &
      !$acc&     inds, indf, indsp, indfp, tauray) async(async_id)
      do k = 1, nlay
         do n = 1, nday_length
            tauray = colmol(n, k) * rayl
            if (k .le. laytrop(n)) then
               speccomb = colamt(n, k,1) + strrat(16)*colamt(n, k,5)
               specmult = 8.0 * min( oneminus, colamt(n, k,1)/speccomb )

               js = 1 + int( specmult )
               fs = mod( specmult, f_one )
               fs1= f_one - fs
               fac000 = fs1 * fac00(n, k)
               fac010 = fs1 * fac10(n, k)
               fac100 = fs  * fac00(n, k)
               fac110 = fs  * fac10(n, k)
               fac001 = fs1 * fac01(n, k)
               fac011 = fs1 * fac11(n, k)
               fac101 = fs  * fac01(n, k)
               fac111 = fs  * fac11(n, k)

               ind01 = id0(n, k,16) + js
               ind02 = ind01 + 1
               ind03 = ind01 + 9
               ind04 = ind01 + 10
               ind11 = id1(n, k,16) + js
               ind12 = ind11 + 1
               ind13 = ind11 + 9
               ind14 = ind11 + 10
               inds = indself(n, k)
               indf = indfor (n, k)
               indsp= inds + 1
               indfp= indf + 1
               
               !$acc loop seq
               do j = 1, ng16
                  !taug(ipt, k,ns16+j, jj) = speccomb                                     &
                  !&        *( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)        &
                  !&        +  fac010 * absa(ind03,j) + fac110 * absa(ind04,j)        &
                  !&        +  fac001 * absa(ind11,j) + fac101 * absa(ind12,j)        &
                  !&        +  fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )      &
                  !&        + colamt(ipt, k,1, jj) * (selffac(ipt, k, jj) * (selfref(inds,j)            &
                  !&        + selffrac(ipt, k, jj) * (selfref(indsp,j)-selfref(inds,j)))       &
                  !&        + forfac(ipt, k, jj) * (forref(indf,j) + forfrac(ipt, k, jj)                &
                  !&        * (forref(indfp,j) - forref(indf,j))))
                  !taur(ipt, k,ns16+j, jj) = tauray
                  small_taug(n, k,j + offset_ng) = speccomb                                     &
                  &        *( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)        &
                  &        +  fac010 * absa(ind03,j) + fac110 * absa(ind04,j)        &
                  &        +  fac001 * absa(ind11,j) + fac101 * absa(ind12,j)        &
                  &        +  fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )      &
                  &        + colamt(n, k,1) * (selffac(n, k) * (selfref(inds,j)            &
                  &        + selffrac(n, k) * (selfref(indsp,j)-selfref(inds,j)))       &
                  &        + forfac(n, k) * (forref(indf,j) + forfrac(n, k)                &
                  &        * (forref(indfp,j) - forref(indf,j))))
                  small_taur(n, k,j + offset_ng) = tauray
               enddo
            else
               ind01 = id0(n, k,16) + 1
               ind02 = ind01 + 1
               ind11 = id1(n, k,16) + 1
               ind12 = ind11 + 1
               
               !$acc loop seq
               do j = 1, ng16
                  !taug(ipt, k,ns16+j, jj) = colamt(ipt, k,5, jj)                                  &
                  !&      * ( fac00(ipt, k, jj)*absb(ind01,j) + fac10(ipt, k, jj)*absb(ind02,j)         &
                  !&      +   fac01(ipt, k, jj)*absb(ind11,j) + fac11(ipt, k, jj)*absb(ind12,j) )
                  !taur(ipt, k,ns16+j, jj) = tauray
                  small_taug(n, k,j + offset_ng) = colamt(n, k,5)                                  &
                  &      * ( fac00(n, k)*absb(ind01,j) + fac10(n, k)*absb(ind02,j)         &
                  &      +   fac01(n, k)*absb(ind11,j) + fac11(n, k)*absb(ind12,j) )
                  small_taur(n, k,j + offset_ng) = tauray
               enddo
            end if
         end do
      end do

      return
!...................................
      end subroutine taumol16
!-----------------------------------


!-----------------------------------
      subroutine taumol17                                                 &
!...................................
     &     ( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
      &       async_id, fulljj, id0, id1, blockjj, ng00, &
      &       map_nday_ipt, map_nday_jj, nday_length, max_nday_length, jjoffset, &
      &       small_sfluxzen, small_taug, small_taur, np, offset_ng)

!  ------------------------------------------------------------------  !
!     band 17:  3250-4000 cm-1 (low - h2o,co2; high - h2o,co2)         !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb17
      
      integer, dimension(nday_length) :: map_nday_ipt, map_nday_jj
      integer :: ng00, np, offset_ng, nday_length, jjoffset, max_nday_length
      real (kind=kind_phys), dimension(nday_length, nlay, ng00) :: small_taur, small_taug
      real (kind=kind_phys), dimension(nday_length, ng00) :: small_sfluxzen
      integer :: nlay, laytrop(max_nday_length), ix, nday(fulljj), blockjj
      integer, dimension(max_nday_length, nlay) :: indfor, indself,          &
     &       jp, jt, jt1
      real (kind=kind_phys), dimension(max_nday_length, nlay) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac
      real (kind=kind_phys), dimension(max_nday_length, nlay, maxgas) :: colamt
      integer, dimension(nday_length, nlay,nblow:nbhgh) :: id0, id1
!  ---  locals:
      integer :: ipt, jj, async_id, fulljj
      real (kind=kind_phys) :: speccomb, specmult, tauray, fs, fs1,     &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111

      integer :: ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14
      integer :: inds, indf, indsp, indfp, j, js, k, n

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.
      
      !$acc parallel loop collapse(2) private(speccomb, specmult, js, fs, fs1, fac000, &
      !$acc&     fac010, fac100, fac110, fac001, fac011, fac101, fac111, &
      !$acc&     ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14, &
      !$acc&     inds, indf, indsp, indfp, tauray) async(async_id)
      do k = 1, nlay
         do n = 1, nday_length
            tauray = colmol(n, k) * rayl
            if (k .le. laytrop(n)) then
               speccomb = colamt(n, k,1) + strrat(17)*colamt(n, k,2)
               specmult = 8.0 * min(oneminus, colamt(n, k,1) / speccomb)

               js = 1 + int(specmult)
               fs = mod(specmult, f_one)
               fs1= f_one - fs
               fac000 = fs1 * fac00(n, k)
               fac010 = fs1 * fac10(n, k)
               fac100 = fs  * fac00(n, k)
               fac110 = fs  * fac10(n, k)
               fac001 = fs1 * fac01(n, k)
               fac011 = fs1 * fac11(n, k)
               fac101 = fs  * fac01(n, k)
               fac111 = fs  * fac11(n, k)

               ind01 = id0(n, k,17) + js
               ind02 = ind01 + 1
               ind03 = ind01 + 9
               ind04 = ind01 + 10
               ind11 = id1(n, k,17) + js
               ind12 = ind11 + 1
               ind13 = ind11 + 9
               ind14 = ind11 + 10

               inds = indself(n, k)
               indf = indfor (n, k)
               indsp= inds + 1
               indfp= indf + 1
               
               !$acc loop seq
               do j = 1, ng17
                  !taug(ipt, k,ns17+j, jj) = speccomb                                     &
                  !&        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
                  !&        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
                  !&        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
                  !&        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )     &
                  !&        + colamt(ipt, k,1, jj) * (selffac(ipt, k, jj) * (selfref(inds,j)            &
                  !&        + selffrac(ipt, k, jj) * (selfref(indsp,j)-selfref(inds,j)))       &
                  !&        + forfac(ipt, k, jj) * (forref(indf,j) + forfrac(ipt, k, jj)                &
                  !&        * (forref(indfp,j) - forref(indf,j))))
                  !taur(ipt, k,ns17+j, jj) = tauray
                  small_taug(n, k,j + offset_ng) = speccomb                                     &
                  &        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
                  &        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
                  &        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
                  &        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )     &
                  &        + colamt(n, k,1) * (selffac(n, k) * (selfref(inds,j)            &
                  &        + selffrac(n, k) * (selfref(indsp,j)-selfref(inds,j)))       &
                  &        + forfac(n, k) * (forref(indf,j) + forfrac(n, k)                &
                  &        * (forref(indfp,j) - forref(indf,j))))
                  small_taur(n, k,j + offset_ng) = tauray
               enddo
            else
               speccomb = colamt(n, k,1) + strrat(17)*colamt(n, k,2)
               specmult = 4.0 * min(oneminus, colamt(n, k,1) / speccomb)

               js = 1 + int(specmult)
               fs = mod(specmult, f_one)
               fs1= f_one - fs
               fac000 = fs1 * fac00(n, k)
               fac010 = fs1 * fac10(n, k)
               fac100 = fs  * fac00(n, k)
               fac110 = fs  * fac10(n, k)
               fac001 = fs1 * fac01(n, k)
               fac011 = fs1 * fac11(n, k)
               fac101 = fs  * fac01(n, k)
               fac111 = fs  * fac11(n, k)

               ind01 = id0(n, k,17) + js
               ind02 = ind01 + 1
               ind03 = ind01 + 5
               ind04 = ind01 + 6
               ind11 = id1(n, k,17) + js
               ind12 = ind11 + 1
               ind13 = ind11 + 5
               ind14 = ind11 + 6

               indf = indfor(n, k)
               indfp= indf + 1
               
               !$acc loop seq
               do j = 1, ng17
                  !taug(ipt, k,ns17+j, jj) = speccomb                                     &
                  !&        * ( fac000 * absb(ind01,j) + fac100 * absb(ind02,j)       &
                  !&        +   fac010 * absb(ind03,j) + fac110 * absb(ind04,j)       &
                  !&        +   fac001 * absb(ind11,j) + fac101 * absb(ind12,j)       &
                  !&        +   fac011 * absb(ind13,j) + fac111 * absb(ind14,j) )     &
                  !&        + colamt(ipt, k,1, jj) * forfac(ipt, k, jj) * (forref(indf,j)               &
                  !&        + forfrac(ipt, k, jj) * (forref(indfp,j) - forref(indf,j)))
                  !taur(ipt, k,ns17+j, jj) = tauray
                  small_taug(n, k,j + offset_ng) = speccomb                                     &
                  &        * ( fac000 * absb(ind01,j) + fac100 * absb(ind02,j)       &
                  &        +   fac010 * absb(ind03,j) + fac110 * absb(ind04,j)       &
                  &        +   fac001 * absb(ind11,j) + fac101 * absb(ind12,j)       &
                  &        +   fac011 * absb(ind13,j) + fac111 * absb(ind14,j) )     &
                  &        + colamt(n, k,1) * forfac(n, k) * (forref(indf,j)               &
                  &        + forfrac(n, k) * (forref(indfp,j) - forref(indf,j)))
                  small_taur(n, k,j + offset_ng) = tauray
               enddo
            end if
         enddo
      end do

      return
!...................................
      end subroutine taumol17
!-----------------------------------


!-----------------------------------
      subroutine taumol18                                                 &
!...................................
     &     ( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
      &       async_id, fulljj, id0, id1, blockjj, ng00, &
      &       map_nday_ipt, map_nday_jj, nday_length, max_nday_length, jjoffset, &
      &       small_sfluxzen, small_taug, small_taur, np, offset_ng)
     
!  ------------------------------------------------------------------  !
!     band 18:  4000-4650 cm-1 (low - h2o,ch4; high - ch4)             !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb18
      
      integer, dimension(nday_length) :: map_nday_ipt, map_nday_jj
      integer :: ng00, np, offset_ng, nday_length, jjoffset, max_nday_length
      real (kind=kind_phys), dimension(nday_length, nlay, ng00) :: small_taur, small_taug
      real (kind=kind_phys), dimension(nday_length, ng00) :: small_sfluxzen
      integer :: nlay, laytrop(max_nday_length), ix, nday(fulljj), blockjj
      integer, dimension(max_nday_length, nlay) :: indfor, indself,          &
     &       jp, jt, jt1
      real (kind=kind_phys), dimension(max_nday_length, nlay) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac
      real (kind=kind_phys), dimension(max_nday_length, nlay, maxgas) :: colamt
      integer, dimension(nday_length, nlay,nblow:nbhgh) :: id0, id1
!  ---  locals:
      integer :: ipt, jj, async_id, fulljj
      real (kind=kind_phys) :: speccomb, specmult, tauray, fs, fs1,     &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111

      integer :: ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14
      integer :: inds, indf, indsp, indfp, j, js, k, n

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.
      !$acc parallel loop collapse(2) private(speccomb, specmult, js, fs, fs1, fac000, &
      !$acc&     fac010, fac100, fac110, fac001, fac011, fac101, fac111, &
      !$acc&     ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14, &
      !$acc&     inds, indf, indsp, indfp, tauray) async(async_id)
      do k = 1, nlay
         do n = 1, nday_length
            tauray = colmol(n, k) * rayl
            if (k .le. laytrop(n)) then
               speccomb = colamt(n, k,1) + strrat(18)*colamt(n, k,5)
               specmult = 8.0 * min(oneminus, colamt(n, k,1) / speccomb)

               js = 1 + int(specmult)
               fs = mod(specmult, f_one)
               fs1= f_one - fs
               fac000 = fs1 * fac00(n, k)
               fac010 = fs1 * fac10(n, k)
               fac100 = fs  * fac00(n, k)
               fac110 = fs  * fac10(n, k)
               fac001 = fs1 * fac01(n, k)
               fac011 = fs1 * fac11(n, k)
               fac101 = fs  * fac01(n, k)
               fac111 = fs  * fac11(n, k)

               ind01 = id0(n, k,18) + js
               ind02 = ind01 + 1
               ind03 = ind01 + 9
               ind04 = ind01 + 10
               ind11 = id1(n, k,18) + js
               ind12 = ind11 + 1
               ind13 = ind11 + 9
               ind14 = ind11 + 10

               inds = indself(n, k)
               indf = indfor (n, k)
               indsp= inds + 1
               indfp= indf + 1
               
               !$acc loop seq
               do j = 1, ng18
                  !taug(ipt, k,ns18+j, jj) = speccomb                                     &
                  !&        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
                  !&        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
                  !&        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
                  !&        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )     &
                  !&        + colamt(ipt, k,1, jj) * (selffac(ipt, k, jj) * (selfref(inds,j)            &
                  !&        + selffrac(ipt, k, jj) * (selfref(indsp,j)-selfref(inds,j)))       &
                  !&        + forfac(ipt, k, jj) * (forref(indf,j) + forfrac(ipt, k, jj)                &
                  !&        * (forref(indfp,j) - forref(indf,j))))
                  !taur(ipt, k,ns18+j, jj) = tauray
                  small_taug(n, k,j + offset_ng) = speccomb                                     &
                  &        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
                  &        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
                  &        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
                  &        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )     &
                  &        + colamt(n, k,1) * (selffac(n, k) * (selfref(inds,j)            &
                  &        + selffrac(n, k) * (selfref(indsp,j)-selfref(inds,j)))       &
                  &        + forfac(n, k) * (forref(indf,j) + forfrac(n, k)                &
                  &        * (forref(indfp,j) - forref(indf,j))))
                  small_taur(n, k,j + offset_ng) = tauray
               enddo
            else
               ind01 = id0(n, k,18) + 1
               ind02 = ind01 + 1
               ind11 = id1(n, k,18) + 1
               ind12 = ind11 + 1
               
               !$acc loop seq
               do j = 1, ng18
                  !taug(ipt, k,ns18+j, jj) = colamt(ipt, k,5, jj)                                  &
                  !&        * ( fac00(ipt, k, jj)*absb(ind01,j) + fac10(ipt, k, jj)*absb(ind02,j)       &
                  !&        +   fac01(ipt, k, jj)*absb(ind11,j) + fac11(ipt, k, jj)*absb(ind12,j) )
                  !taur(ipt, k,ns18+j, jj) = tauray
                  small_taug(n, k,j + offset_ng) = colamt(n, k,5)                                  &
                  &        * ( fac00(n, k)*absb(ind01,j) + fac10(n, k)*absb(ind02,j)       &
                  &        +   fac01(n, k)*absb(ind11,j) + fac11(n, k)*absb(ind12,j) )
                  small_taur(n, k,j + offset_ng) = tauray
               enddo
            end if
         end do
      end do

      return
!...................................
      end subroutine taumol18
!-----------------------------------


!-----------------------------------
      subroutine taumol19                                                 &
!...................................
     &     ( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
      &       async_id, fulljj, id0, id1, blockjj, ng00, &
      &       map_nday_ipt, map_nday_jj, nday_length, max_nday_length, jjoffset, &
      &       small_sfluxzen, small_taug, small_taur, np, offset_ng)
     
!  ------------------------------------------------------------------  !
!     band 19:  4650-5150 cm-1 (low - h2o,co2; high - co2)             !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb19
      
      integer, dimension(nday_length) :: map_nday_ipt, map_nday_jj
      integer :: ng00, np, offset_ng, nday_length, jjoffset, max_nday_length
      real (kind=kind_phys), dimension(nday_length, nlay, ng00) :: small_taur, small_taug
      real (kind=kind_phys), dimension(nday_length, ng00) :: small_sfluxzen
      integer :: nlay, laytrop(max_nday_length), ix, nday(fulljj), blockjj
      integer, dimension(max_nday_length, nlay) :: indfor, indself,          &
     &       jp, jt, jt1
      real (kind=kind_phys), dimension(max_nday_length, nlay) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac
      real (kind=kind_phys), dimension(max_nday_length, nlay, maxgas) :: colamt
      integer, dimension(nday_length, nlay,nblow:nbhgh) :: id0, id1
!  ---  locals:
      integer :: ipt, jj, async_id, fulljj
      real (kind=kind_phys) :: speccomb, specmult, tauray, fs, fs1,     &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111

      integer :: ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14
      integer :: inds, indf, indsp, indfp, j, js, k, n

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.
      
      !$acc parallel loop gang collapse(2) private(speccomb, specmult, js, fs, fs1, fac000, &
      !$acc&     fac010, fac100, fac110, fac001, fac011, fac101, fac111, &
      !$acc&     ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14, &
      !$acc&     indsp, indfp, tauray) async(async_id)
      do k = 1, nlay
         do n = 1, nday_length
            tauray = colmol(n, k) * rayl
            if (k .le. laytrop(n)) then
               speccomb = colamt(n, k,1) + strrat(19)*colamt(n, k,2)
               specmult = 8.0 * min(oneminus, colamt(n, k,1) / speccomb)

               js = 1 + int(specmult)
               fs = mod(specmult, f_one)
               fs1= f_one - fs
               fac000 = fs1 * fac00(n, k)
               fac010 = fs1 * fac10(n, k)
               fac100 = fs  * fac00(n, k)
               fac110 = fs  * fac10(n, k)
               fac001 = fs1 * fac01(n, k)
               fac011 = fs1 * fac11(n, k)
               fac101 = fs  * fac01(n, k)
               fac111 = fs  * fac11(n, k)

               ind01 = id0(n, k,19) + js
               ind02 = ind01 + 1
               ind03 = ind01 + 9
               ind04 = ind01 + 10
               ind11 = id1(n, k,19) + js
               ind12 = ind11 + 1
               ind13 = ind11 + 9
               ind14 = ind11 + 10

               inds = indself(n, k)
               indf = indfor (n, k)
               indsp= inds + 1
               indfp= indf + 1
               
               !$acc loop seq
               do j = 1, ng19
                  !taug(ipt, k,ns19+j, jj) = speccomb                                     &
                  !&        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
                  !&        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
                  !&        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
                  !&        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )     &
                  !&        + colamt(ipt, k,1, jj) * (selffac(ipt, k, jj) * (selfref(inds,j)            &
                  !&        + selffrac(ipt, k, jj) * (selfref(indsp,j)-selfref(inds,j)))       &
                  !&        + forfac(ipt, k, jj) * (forref(indf,j) + forfrac(ipt, k, jj)                &
                  !&        * (forref(indfp,j) - forref(indf,j))))
                  !taur(ipt, k,ns19+j, jj) = tauray
                  small_taug(n, k,j + offset_ng) = speccomb                                     &
                  &        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
                  &        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
                  &        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
                  &        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )     &
                  &        + colamt(n, k,1) * (selffac(n, k) * (selfref(inds,j)            &
                  &        + selffrac(n, k) * (selfref(indsp,j)-selfref(inds,j)))       &
                  &        + forfac(n, k) * (forref(indf,j) + forfrac(n, k)                &
                  &        * (forref(indfp,j) - forref(indf,j))))
                  small_taur(n, k,j + offset_ng) = tauray
               enddo
            else
               ind01 = id0(n, k,19) + 1
               ind02 = ind01 + 1
               ind11 = id1(n, k,19) + 1
               ind12 = ind11 + 1
               
               !$acc loop seq
               do j = 1, ng19
                  !taug(ipt, k,ns19+j, jj) = colamt(ipt, k,2, jj)                                  &
                  !&        * ( fac00(ipt, k, jj)*absb(ind01,j) + fac10(ipt, k, jj)*absb(ind02,j)       &
                  !&        +   fac01(ipt, k, jj)*absb(ind11,j) + fac11(ipt, k, jj)*absb(ind12,j) ) 
                  !taur(ipt, k,ns19+j, jj) = tauray
                  small_taug(n, k,j + offset_ng) = colamt(n, k,2)                                  &
                  &        * ( fac00(n, k)*absb(ind01,j) + fac10(n, k)*absb(ind02,j)       &
                  &        +   fac01(n, k)*absb(ind11,j) + fac11(n, k)*absb(ind12,j) ) 
                  small_taur(n, k,j + offset_ng) = tauray
               enddo
            end if
         end do
      end do

!...................................
      end subroutine taumol19
!-----------------------------------


!-----------------------------------
      subroutine taumol20                                                 &
!...................................
     &     ( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
      &       async_id, fulljj, id0, id1, blockjj, ng00, &
      &       map_nday_ipt, map_nday_jj, nday_length, max_nday_length, jjoffset, &
      &       small_sfluxzen, small_taug, small_taur, np, offset_ng)
     
!  ------------------------------------------------------------------  !
!     band 20:  5150-6150 cm-1 (low - h2o; high - h2o)                 !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb20
      
      integer, dimension(nday_length) :: map_nday_ipt, map_nday_jj
      integer :: ng00, np, offset_ng, nday_length, jjoffset, max_nday_length
      real (kind=kind_phys), dimension(nday_length, nlay, ng00) :: small_taur, small_taug
      real (kind=kind_phys), dimension(nday_length, ng00) :: small_sfluxzen
      integer :: nlay, laytrop(max_nday_length), ix, nday(fulljj), blockjj
      integer, dimension(max_nday_length, nlay) :: indfor, indself,          &
     &       jp, jt, jt1
      real (kind=kind_phys), dimension(max_nday_length, nlay) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac
      real (kind=kind_phys), dimension(max_nday_length, nlay, maxgas) :: colamt
      integer, dimension(nday_length, nlay,nblow:nbhgh) :: id0, id1
!  ---  locals:
      integer :: ipt, jj, async_id, fulljj
      real (kind=kind_phys) :: tauray

      integer :: ind01, ind02, ind11, ind12
      integer :: inds, indf, indsp, indfp, j, k, n

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.

      !$acc parallel loop gang collapse(2) private(ind01, ind02, ind11, ind12, inds, indf, &
      !$acc&     indsp, indfp, tauray) async(async_id)
      do k = 1, nlay
         do n = 1, nday_length
            tauray = colmol(n, k) * rayl
            if (k .le. laytrop(n)) then
               ind01 = id0(n, k,20) + 1
               ind02 = ind01 + 1
               ind11 = id1(n, k,20) + 1
               ind12 = ind11 + 1

               inds = indself(n, k)
               indf = indfor (n, k)
               indsp= inds + 1
               indfp= indf + 1
               
               !$acc loop seq
               do j = 1, ng20
                  !taug(ipt, k,ns20+j, jj) = colamt(ipt, k,1, jj)                                  &
                  !&        * ( (fac00(ipt, k, jj)*absa(ind01,j) + fac10(ipt, k, jj)*absa(ind02,j)      &
                  !&        +    fac01(ipt, k, jj)*absa(ind11,j) + fac11(ipt, k, jj)*absa(ind12,j))     &
                  !&        +   selffac(ipt, k, jj) * (selfref(inds,j) + selffrac(ipt, k, jj)           &
                  !&        *   (selfref(indsp,j) - selfref(inds,j)))                 &
                  !&        +   forfac(ipt, k, jj) * (forref(indf,j) + forfrac(ipt, k, jj)              &
                  !&        *   (forref(indfp,j) - forref(indf,j))) )                 &
                  !&        + colamt(ipt, k,5, jj) * absch4(j)
                  !taur(ipt, k,ns20+j, jj) = tauray
                  small_taug(n, k,j + offset_ng) = colamt(n, k,1)                                  &
                  &        * ( (fac00(n, k)*absa(ind01,j) + fac10(n, k)*absa(ind02,j)      &
                  &        +    fac01(n, k)*absa(ind11,j) + fac11(n, k)*absa(ind12,j))     &
                  &        +   selffac(n, k) * (selfref(inds,j) + selffrac(n, k)           &
                  &        *   (selfref(indsp,j) - selfref(inds,j)))                 &
                  &        +   forfac(n, k) * (forref(indf,j) + forfrac(n, k)              &
                  &        *   (forref(indfp,j) - forref(indf,j))) )                 &
                  &        + colamt(n, k,5) * absch4(j)
                  small_taur(n, k,j + offset_ng) = tauray
               enddo
            else
               ind01 = id0(n, k,20) + 1
               ind02 = ind01 + 1
               ind11 = id1(n, k,20) + 1
               ind12 = ind11 + 1

               indf = indfor(n, k)
               indfp= indf + 1

               !$acc loop seq
               do j = 1, ng20
                  !taug(ipt, k,ns20+j, jj) = colamt(ipt, k,1, jj)                                  &
                  !&        * ( fac00(ipt, k, jj)*absb(ind01,j) + fac10(ipt, k, jj)*absb(ind02,j)       &
                  !&        +   fac01(ipt, k, jj)*absb(ind11,j) + fac11(ipt, k, jj)*absb(ind12,j)       &
                  !&        +   forfac(ipt, k, jj) * (forref(indf,j) + forfrac(ipt, k, jj)              &
                  !&        *   (forref(indfp,j) - forref(indf,j))) )                 &
                  !&        + colamt(ipt, k,5, jj) * absch4(j)
                  !taur(ipt, k,ns20+j, jj) = tauray
                  small_taug(n, k,j + offset_ng) = colamt(n, k,1)                                  &
                  &        * ( fac00(n, k)*absb(ind01,j) + fac10(n, k)*absb(ind02,j)       &
                  &        +   fac01(n, k)*absb(ind11,j) + fac11(n, k)*absb(ind12,j)       &
                  &        +   forfac(n, k) * (forref(indf,j) + forfrac(n, k)              &
                  &        *   (forref(indfp,j) - forref(indf,j))) )                 &
                  &        + colamt(n, k,5) * absch4(j)
                  small_taur(n, k,j + offset_ng) = tauray
               enddo
            end if
         end do
      end do

      return
!...................................
      end subroutine taumol20
!-----------------------------------


!-----------------------------------
      subroutine taumol21                                                 &
!...................................
     &     ( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
      &       async_id, fulljj, id0, id1, blockjj, ng00, &
      &       map_nday_ipt, map_nday_jj, nday_length, max_nday_length, jjoffset, &
      &       small_sfluxzen, small_taug, small_taur, np, offset_ng)
     
!  ------------------------------------------------------------------  !
!     band 21:  6150-7700 cm-1 (low - h2o,co2; high - h2o,co2)         !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb21
      
      integer, dimension(nday_length) :: map_nday_ipt, map_nday_jj
      integer :: ng00, np, offset_ng, nday_length, jjoffset, max_nday_length
      real (kind=kind_phys), dimension(nday_length, nlay, ng00) :: small_taur, small_taug
      real (kind=kind_phys), dimension(nday_length, ng00) :: small_sfluxzen
      integer :: nlay, laytrop(max_nday_length), ix, nday(fulljj), blockjj
      integer, dimension(max_nday_length, nlay) :: indfor, indself,          &
     &       jp, jt, jt1
      real (kind=kind_phys), dimension(max_nday_length, nlay) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac
      real (kind=kind_phys), dimension(max_nday_length, nlay, maxgas) :: colamt
      integer, dimension(nday_length, nlay,nblow:nbhgh) :: id0, id1
!  ---  locals:
      integer :: ipt, jj, async_id, fulljj
      real (kind=kind_phys) :: speccomb, specmult, tauray, fs, fs1,     &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111

      integer :: ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14
      integer :: inds, indf, indsp, indfp, j, js, k, n

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.
      
      !$acc parallel loop gang collapse(2) private(speccomb, specmult, js, fs, fs1, fac000, &
      !$acc&     fac010, fac100, fac110, fac001, fac011, fac101, fac111, &
      !$acc&     ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14, &
      !$acc&     inds, indf, indsp, indfp, tauray) async(async_id)
      do k = 1, nlay
         do n = 1, nday_length
            tauray = colmol(n, k) * rayl
            if (k .le. laytrop(n)) then
               speccomb = colamt(n, k,1) + strrat(21)*colamt(n, k,2)
               specmult = 8.0 * min(oneminus, colamt(n, k,1) / speccomb)

               js = 1 + int(specmult)
               fs = mod(specmult, f_one)
               fs1= f_one - fs
               fac000 = fs1 * fac00(n, k)
               fac010 = fs1 * fac10(n, k)
               fac100 = fs  * fac00(n, k)
               fac110 = fs  * fac10(n, k)
               fac001 = fs1 * fac01(n, k)
               fac011 = fs1 * fac11(n, k)
               fac101 = fs  * fac01(n, k)
               fac111 = fs  * fac11(n, k)

               ind01 = id0(n, k,21) + js
               ind02 = ind01 + 1
               ind03 = ind01 + 9
               ind04 = ind01 + 10
               ind11 = id1(n, k,21) + js
               ind12 = ind11 + 1
               ind13 = ind11 + 9
               ind14 = ind11 + 10

               inds = indself(n, k)
               indf = indfor (n, k)
               indsp= inds + 1
               indfp= indf + 1
               
               !$acc loop seq
               do j = 1, ng21
                  !taug(ipt, k,ns21+j, jj) = speccomb                                     &
                  !&        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
                  !&        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
                  !&        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
                  !&        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )     &
                  !&        + colamt(ipt, k,1, jj) * (selffac(ipt, k, jj) * (selfref(inds,j)            &
                  !&        + selffrac(ipt, k, jj) * (selfref(indsp,j) - selfref(inds,j)))     &
                  !&        + forfac(ipt, k, jj) * (forref(indf,j) + forfrac(ipt, k, jj)                &
                  !&        * (forref(indfp,j) - forref(indf,j))))
                  !taur(ipt, k,ns21+j, jj) = tauray
                  small_taug(n, k,j + offset_ng) = speccomb                                     &
                  &        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
                  &        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
                  &        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
                  &        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )     &
                  &        + colamt(n, k,1) * (selffac(n, k) * (selfref(inds,j)            &
                  &        + selffrac(n, k) * (selfref(indsp,j) - selfref(inds,j)))     &
                  &        + forfac(n, k) * (forref(indf,j) + forfrac(n, k)                &
                  &        * (forref(indfp,j) - forref(indf,j))))
                  small_taur(n, k,j + offset_ng) = tauray
               enddo
            else
               speccomb = colamt(n, k,1) + strrat(21)*colamt(n, k,2)
               specmult = 4.0 * min(oneminus, colamt(n, k,1) / speccomb)

               js = 1 + int(specmult)
               fs = mod(specmult, f_one)
               fs1= f_one - fs
               fac000 = fs1 * fac00(n, k)
               fac010 = fs1 * fac10(n, k)
               fac100 = fs  * fac00(n, k)
               fac110 = fs  * fac10(n, k)
               fac001 = fs1 * fac01(n, k)
               fac011 = fs1 * fac11(n, k)
               fac101 = fs  * fac01(n, k)
               fac111 = fs  * fac11(n, k)

               ind01 = id0(n, k,21) + js
               ind02 = ind01 + 1
               ind03 = ind01 + 5
               ind04 = ind01 + 6
               ind11 = id1(n, k,21) + js
               ind12 = ind11 + 1
               ind13 = ind11 + 5
               ind14 = ind11 + 6

               indf = indfor(n, k)
               indfp= indf + 1
               
               !$acc loop seq
               do j = 1, ng21
                  !taug(ipt, k,ns21+j, jj) = speccomb                                     &
                  !&        * ( fac000 * absb(ind01,j) + fac100 * absb(ind02,j)       &
                  !&        +   fac010 * absb(ind03,j) + fac110 * absb(ind04,j)       &
                  !&        +   fac001 * absb(ind11,j) + fac101 * absb(ind12,j)       &
                  !&        +   fac011 * absb(ind13,j) + fac111 * absb(ind14,j) )     &
                  !&        + colamt(ipt, k,1, jj) * forfac(ipt, k, jj) * (forref(indf,j)               &
                  !&        + forfrac(ipt, k, jj) * (forref(indfp,j) - forref(indf,j)))
                  !taur(ipt, k,ns21+j, jj) = tauray
                  small_taug(n, k,j + offset_ng) = speccomb                                     &
                  &        * ( fac000 * absb(ind01,j) + fac100 * absb(ind02,j)       &
                  &        +   fac010 * absb(ind03,j) + fac110 * absb(ind04,j)       &
                  &        +   fac001 * absb(ind11,j) + fac101 * absb(ind12,j)       &
                  &        +   fac011 * absb(ind13,j) + fac111 * absb(ind14,j) )     &
                  &        + colamt(n, k,1) * forfac(n, k) * (forref(indf,j)               &
                  &        + forfrac(n, k) * (forref(indfp,j) - forref(indf,j)))
                  small_taur(n, k,j + offset_ng) = tauray
               enddo
            end if
         enddo
      end do

!...................................
      end subroutine taumol21
!-----------------------------------


!-----------------------------------
      subroutine taumol22                                                 &
!...................................
     &     ( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
      &       async_id, fulljj, id0, id1, blockjj, ng00, &
      &       map_nday_ipt, map_nday_jj, nday_length, max_nday_length, jjoffset, &
      &       small_sfluxzen, small_taug, small_taur, np, offset_ng)
     
!  ------------------------------------------------------------------  !
!     band 22:  7700-8050 cm-1 (low - h2o,o2; high - o2)               !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb22
      
      integer, dimension(nday_length) :: map_nday_ipt, map_nday_jj
      integer :: ng00, np, offset_ng, nday_length, jjoffset, max_nday_length
      real (kind=kind_phys), dimension(nday_length, nlay, ng00) :: small_taur, small_taug
      real (kind=kind_phys), dimension(nday_length, ng00) :: small_sfluxzen
      integer :: nlay, laytrop(max_nday_length), ix, nday(fulljj), blockjj
      integer, dimension(max_nday_length, nlay) :: indfor, indself,          &
     &       jp, jt, jt1
      real (kind=kind_phys), dimension(max_nday_length, nlay) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac
      real (kind=kind_phys), dimension(max_nday_length, nlay, maxgas) :: colamt
      integer, dimension(nday_length, nlay,nblow:nbhgh) :: id0, id1
!  ---  locals:
      integer :: ipt, jj, async_id, fulljj
      real (kind=kind_phys) :: speccomb, specmult, tauray, fs, fs1,     &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111,  &
     &       o2adj, o2cont, o2tem

      integer :: ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14
      integer :: inds, indf, indsp, indfp, j, js, k, n

!
!===> ... begin here
!
!  --- ...  the following factor is the ratio of total o2 band intensity (lines
!           and mate continuum) to o2 band intensity (line only). it is needed
!           to adjust the optical depths since the k's include only lines.

      o2adj = 1.6
      o2tem = 4.35e-4 / (350.0*2.0)
      
!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.
      
      !$acc parallel loop gang collapse(2) private(speccomb, specmult, js, fs, fs1, fac000, &
      !$acc&     fac010, fac100, fac110, fac001, fac011, fac101, fac111, &
      !$acc&     ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14, &
      !$acc&     inds, indf, indsp, indfp, o2cont, tauray) async(async_id)
      do k = 1, nlay
         do n = 1, nday_length
            tauray = colmol(n, k) * rayl
            if (k .le. laytrop(n)) then
               o2cont   = o2tem * colamt(n, k,6)
               speccomb = colamt(n, k,1) + strrat(22)*colamt(n, k,6)
               specmult = 8.0 * min(oneminus, colamt(n, k,1) / speccomb)

               js = 1 + int(specmult)
               fs = mod(specmult, f_one)
               fs1= f_one - fs
               fac000 = fs1 * fac00(n, k)
               fac010 = fs1 * fac10(n, k)
               fac100 = fs  * fac00(n, k)
               fac110 = fs  * fac10(n, k)
               fac001 = fs1 * fac01(n, k)
               fac011 = fs1 * fac11(n, k)
               fac101 = fs  * fac01(n, k)
               fac111 = fs  * fac11(n, k)

               ind01 = id0(n, k,22) + js
               ind02 = ind01 + 1
               ind03 = ind01 + 9
               ind04 = ind01 + 10
               ind11 = id1(n, k,22) + js
               ind12 = ind11 + 1
               ind13 = ind11 + 9
               ind14 = ind11 + 10

               inds = indself(n, k)
               indf = indfor (n, k)
               indsp= inds + 1
               indfp= indf + 1
               
               !$acc loop seq
               do j = 1, ng22
                  !taug(ipt, k,ns22+j, jj) = speccomb                                     &
                  !&        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
                  !&        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
                  !&        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
                  !&        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )     &
                  !&        + colamt(ipt, k,1, jj) * (selffac(ipt, k, jj) * (selfref(inds,j)            &
                  !&        + selffrac(ipt, k, jj) * (selfref(indsp,j)-selfref(inds,j)))       &
                  !&        + forfac(ipt, k, jj) * (forref(indf,j) + forfrac(ipt, k, jj)                &
                  !&        * (forref(indfp,j) - forref(indf,j)))) + o2cont
                  !taur(ipt, k,ns22+j, jj) = tauray
                  small_taug(n, k,j + offset_ng) = speccomb                                     &
                  &        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
                  &        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
                  &        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
                  &        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )     &
                  &        + colamt(n, k,1) * (selffac(n, k) * (selfref(inds,j)            &
                  &        + selffrac(n, k) * (selfref(indsp,j)-selfref(inds,j)))       &
                  &        + forfac(n, k) * (forref(indf,j) + forfrac(n, k)                &
                  &        * (forref(indfp,j) - forref(indf,j)))) + o2cont
                  small_taur(n, k,j + offset_ng) = tauray
               enddo
            else
               o2cont = o2tem * colamt(n, k,6)

               ind01 = id0(n, k,22) + 1
               ind02 = ind01 + 1
               ind11 = id1(n, k,22) + 1
               ind12 = ind11 + 1
               
               !$acc loop seq
               do j = 1, ng22
                  !taug(ipt, k,ns22+j, jj) = colamt(ipt, k,6, jj) * o2adj                          &
                  !&        * ( fac00(ipt, k, jj)*absb(ind01,j) + fac10(ipt, k, jj)*absb(ind02,j)       &
                  !&        +   fac01(ipt, k, jj)*absb(ind11,j) + fac11(ipt, k, jj)*absb(ind12,j) )     &
                  !&        + o2cont
                  !taur(ipt, k,ns22+j, jj) = tauray
                  small_taug(n, k,j + offset_ng) = colamt(n, k,6) * o2adj                          &
                  &        * ( fac00(n, k)*absb(ind01,j) + fac10(n, k)*absb(ind02,j)       &
                  &        +   fac01(n, k)*absb(ind11,j) + fac11(n, k)*absb(ind12,j) )     &
                  &        + o2cont
                  small_taur(n, k,j + offset_ng) = tauray
               enddo
            end if
         enddo
      end do

      return
!...................................
      end subroutine taumol22
!-----------------------------------


!-----------------------------------
      subroutine taumol23                                                 &
!...................................
     &     ( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
      &       async_id, fulljj, id0, id1, blockjj, ng00, &
      &       map_nday_ipt, map_nday_jj, nday_length, max_nday_length, jjoffset, &
      &       small_sfluxzen, small_taug, small_taur, np, offset_ng)
     
!  ------------------------------------------------------------------  !
!     band 23:  8050-12850 cm-1 (low - h2o; high - nothing)            !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb23
      
      integer, dimension(nday_length) :: map_nday_ipt, map_nday_jj
      integer :: ng00, np, offset_ng, nday_length, jjoffset, max_nday_length
      real (kind=kind_phys), dimension(nday_length, nlay, ng00) :: small_taur, small_taug
      real (kind=kind_phys), dimension(nday_length, ng00) :: small_sfluxzen
      integer :: nlay, laytrop(max_nday_length), ix, nday(fulljj), blockjj
      integer, dimension(max_nday_length, nlay) :: indfor, indself,          &
     &       jp, jt, jt1
      real (kind=kind_phys), dimension(max_nday_length, nlay) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac
      real (kind=kind_phys), dimension(max_nday_length, nlay, maxgas) :: colamt
      integer, dimension(nday_length, nlay,nblow:nbhgh) :: id0, id1
!  ---  locals:
      integer :: ipt, jj, async_id, fulljj
      integer :: ind01, ind02, ind11, ind12
      integer :: inds, indf, indsp, indfp, j, k, n

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.
      !$acc parallel loop gang collapse(2) private(ind01, ind02, ind11, ind12, inds, indf, &
      !$acc&     indsp, indfp) async(async_id)
      do k = 1, nlay
         do n = 1, nday_length
            if (k .le. laytrop(n)) then
               ind01 = id0(n, k,23) + 1
               ind02 = ind01 + 1
               ind11 = id1(n, k,23) + 1
               ind12 = ind11 + 1

               inds = indself(n, k)
               indf = indfor (n, k)
               indsp= inds + 1
               indfp= indf + 1
               
               !$acc loop seq
               do j = 1, ng23
                  !taug(ipt, k,ns23+j, jj) = colamt(ipt, k,1, jj) * (givfac                        &
                  !&        * ( fac00(ipt, k, jj)*absa(ind01,j) + fac10(ipt, k, jj)*absa(ind02,j)       &
                  !&        +   fac01(ipt, k, jj)*absa(ind11,j) + fac11(ipt, k, jj)*absa(ind12,j) )     &
                  !&        + selffac(ipt, k, jj) * (selfref(inds,j) + selffrac(ipt, k, jj)             &
                  !&        * (selfref(indsp,j) - selfref(inds,j)))                   &
                  !&        + forfac(ipt, k, jj) * (forref(indf,j) + forfrac(ipt, k, jj)                &
                  !&        * (forref(indfp,j) - forref(indf,j))))
                  !taur(ipt, k,ns23+j, jj) = colmol(ipt, k, jj) * rayl(j)
                  small_taug(n, k,j + offset_ng) = colamt(n, k,1) * (givfac                        &
                  &        * ( fac00(n, k)*absa(ind01,j) + fac10(n, k)*absa(ind02,j)       &
                  &        +   fac01(n, k)*absa(ind11,j) + fac11(n, k)*absa(ind12,j) )     &
                  &        + selffac(n, k) * (selfref(inds,j) + selffrac(n, k)             &
                  &        * (selfref(indsp,j) - selfref(inds,j)))                   &
                  &        + forfac(n, k) * (forref(indf,j) + forfrac(n, k)                &
                  &        * (forref(indfp,j) - forref(indf,j))))
                  small_taur(n, k,j + offset_ng) = colmol(n, k) * rayl(j)
               enddo
            else
               !$acc loop seq
               do j = 1, ng23
                  !taug(ipt, k,ns23+j, jj) = f_zero
                  !taur(ipt, k,ns23+j, jj) = colmol(ipt, k, jj) * rayl(j)
                  small_taug(n, k,j + offset_ng) = f_zero
                  small_taur(n, k,j + offset_ng) = colmol(n, k) * rayl(j)
               enddo
            end if
         end do
      end do

!...................................
      end subroutine taumol23
!-----------------------------------


!-----------------------------------
      subroutine taumol24                                                 &
!...................................
     &     ( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
      &       async_id, fulljj, id0, id1, blockjj, ng00, &
      &       map_nday_ipt, map_nday_jj, nday_length, max_nday_length, jjoffset, &
      &       small_sfluxzen, small_taug, small_taur, np, offset_ng)
     
!  ------------------------------------------------------------------  !
!     band 24:  12850-16000 cm-1 (low - h2o,o2; high - o2)             !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb24
      
      integer, dimension(nday_length) :: map_nday_ipt, map_nday_jj
      integer :: ng00, np, offset_ng, nday_length, jjoffset, max_nday_length
      real (kind=kind_phys), dimension(nday_length, nlay, ng00) :: small_taur, small_taug
      real (kind=kind_phys), dimension(nday_length, ng00) :: small_sfluxzen
      integer :: nlay, laytrop(max_nday_length), ix, nday(fulljj), blockjj
      integer, dimension(max_nday_length, nlay) :: indfor, indself,          &
     &       jp, jt, jt1
      real (kind=kind_phys), dimension(max_nday_length, nlay) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac
      real (kind=kind_phys), dimension(max_nday_length, nlay, maxgas) :: colamt
      integer, dimension(nday_length, nlay,nblow:nbhgh) :: id0, id1
!  ---  locals:
      integer :: ipt, jj, async_id, fulljj
      real (kind=kind_phys) :: speccomb, specmult, fs, fs1,             &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111

      integer :: ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14
      integer :: inds, indf, indsp, indfp, j, js, k, n

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.
      !$acc parallel loop gang collapse(2) private(speccomb, specmult, js, fs, fs1, fac000, &
      !$acc&     fac010, fac100, fac110, fac001, fac011, fac101, fac111, &
      !$acc&     ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14, &
      !$acc&     inds, indf, indsp, indfp) async(async_id)
      do k = 1, nlay
         do n = 1, nday_length
            if (k .le. laytrop(n)) then
               speccomb = colamt(n, k,1) + strrat(24)*colamt(n, k,6)
               specmult = 8.0 * min(oneminus, colamt(n, k,1) / speccomb)

               js = 1 + int(specmult)
               fs = mod(specmult, f_one)
               fs1= f_one - fs
               fac000 = fs1 * fac00(n, k)
               fac010 = fs1 * fac10(n, k)
               fac100 = fs  * fac00(n, k)
               fac110 = fs  * fac10(n, k)
               fac001 = fs1 * fac01(n, k)
               fac011 = fs1 * fac11(n, k)
               fac101 = fs  * fac01(n, k)
               fac111 = fs  * fac11(n, k)

               ind01 = id0(n, k,24) + js
               ind02 = ind01 + 1
               ind03 = ind01 + 9
               ind04 = ind01 + 10
               ind11 = id1(n, k,24) + js
               ind12 = ind11 + 1
               ind13 = ind11 + 9
               ind14 = ind11 + 10

               inds = indself(n, k)
               indf = indfor (n, k)
               indsp= inds + 1
               indfp= indf + 1
               
               !$acc loop seq
               do j = 1, ng24
                  !taug(ipt, k,ns24+j, jj) = speccomb                                     &
                  !&        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
                  !&        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
                  !&        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
                  !&        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )     &
                  !&        + colamt(ipt, k,3, jj) * abso3a(j) + colamt(ipt, k,1, jj)                   &
                  !&        * (selffac(ipt, k, jj) * (selfref(inds,j) + selffrac(ipt, k, jj)            &
                  !&        * (selfref(indsp,j) - selfref(inds,j)))                   &
                  !&        + forfac(ipt, k, jj) * (forref(indf,j) + forfrac(ipt, k, jj)                &
                  !&        * (forref(indfp,j) - forref(indf,j))))
!
                  !taur(ipt, k,ns24+j, jj) = colmol(ipt, k, jj)                                    &
                  !&           * (rayla(j,js) + fs*(rayla(j,js+1) - rayla(j,js)))
                  small_taug(n, k,j + offset_ng) = speccomb                                     &
                  &        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
                  &        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
                  &        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
                  &        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )     &
                  &        + colamt(n, k,3) * abso3a(j) + colamt(n, k,1)                   &
                  &        * (selffac(n, k) * (selfref(inds,j) + selffrac(n, k)            &
                  &        * (selfref(indsp,j) - selfref(inds,j)))                   &
                  &        + forfac(n, k) * (forref(indf,j) + forfrac(n, k)                &
                  &        * (forref(indfp,j) - forref(indf,j))))

                  small_taur(n, k,j + offset_ng) = colmol(n, k)                                    &
                  &           * (rayla(j,js) + fs*(rayla(j,js+1) - rayla(j,js)))
               enddo
            else
               ind01 = id0(n, k,24) + 1
               ind02 = ind01 + 1
               ind11 = id1(n, k,24) + 1
               ind12 = ind11 + 1
               
               !$acc loop seq
               do j = 1, ng24
                  !taug(ipt, k,ns24+j, jj) = colamt(ipt, k,6, jj)                                  &
                  !&        * ( fac00(ipt, k, jj)*absb(ind01,j) + fac10(ipt, k, jj)*absb(ind02,j)       &
                  !&        +   fac01(ipt, k, jj)*absb(ind11,j) + fac11(ipt, k, jj)*absb(ind12,j) )     &
                  !&        + colamt(ipt, k,3, jj) * abso3b(j)
!
                  !taur(ipt, k,ns24+j, jj) = colmol(ipt, k, jj) * raylb(j)
                  small_taug(n, k,j + offset_ng) = colamt(n, k,6)                                  &
                  &        * ( fac00(n, k)*absb(ind01,j) + fac10(n, k)*absb(ind02,j)       &
                  &        +   fac01(n, k)*absb(ind11,j) + fac11(n, k)*absb(ind12,j) )     &
                  &        + colamt(n, k,3) * abso3b(j)

                  small_taur(n, k,j + offset_ng) = colmol(n, k) * raylb(j)
               enddo
            end if
         end do
      end do

      return
!...................................
      end subroutine taumol24
!-----------------------------------


!-----------------------------------
      subroutine taumol25                                                 &
!...................................
     &     ( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
      &       async_id, fulljj, id0, id1, blockjj, ng00, &
      &       map_nday_ipt, map_nday_jj, nday_length, max_nday_length, jjoffset, &
      &       small_sfluxzen, small_taug, small_taur, np, offset_ng)
     
!  ------------------------------------------------------------------  !
!     band 25:  16000-22650 cm-1 (low - h2o; high - nothing)           !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb25
      
      integer, dimension(nday_length) :: map_nday_ipt, map_nday_jj
      integer :: ng00, np, offset_ng, nday_length, jjoffset, max_nday_length
      real (kind=kind_phys), dimension(nday_length, nlay, ng00) :: small_taur, small_taug
      real (kind=kind_phys), dimension(nday_length, ng00) :: small_sfluxzen
      integer :: nlay, laytrop(max_nday_length), ix, nday(fulljj), blockjj
      integer, dimension(max_nday_length, nlay) :: indfor, indself,          &
     &       jp, jt, jt1
      real (kind=kind_phys), dimension(max_nday_length, nlay) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac
      real (kind=kind_phys), dimension(max_nday_length, nlay, maxgas) :: colamt
      integer, dimension(nday_length, nlay,nblow:nbhgh) :: id0, id1
!  ---  locals:
      integer :: ipt, jj, async_id, fulljj
      integer :: ind01, ind02, ind11, ind12
      integer :: j, k, n

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.
      
      !$acc parallel loop gang collapse(2) private(ind01, ind02, ind11, ind12) async(async_id)
      do k = 1, nlay
         do n = 1, nday_length
            if (k .le. laytrop(n)) then
               ind01 = id0(n, k,25) + 1
               ind02 = ind01 + 1
               ind11 = id1(n, k,25) + 1
               ind12 = ind11 + 1
               
               !$acc loop seq
               do j = 1, ng25
                  !taug(ipt, k,ns25+j, jj) = colamt(ipt, k,1, jj)                                  &
                  !&        * ( fac00(ipt, k, jj)*absa(ind01,j) + fac10(ipt, k, jj)*absa(ind02,j)       &
                  !&        +   fac01(ipt, k, jj)*absa(ind11,j) + fac11(ipt, k, jj)*absa(ind12,j) )     &
                  !&        + colamt(ipt, k,3, jj) * abso3a(j) 
                  !taur(ipt, k,ns25+j, jj) = colmol(ipt, k, jj) * rayl(j)
                  small_taug(n, k,j + offset_ng) = colamt(n, k,1)                                  &
                  &        * ( fac00(n, k)*absa(ind01,j) + fac10(n, k)*absa(ind02,j)       &
                  &        +   fac01(n, k)*absa(ind11,j) + fac11(n, k)*absa(ind12,j) )     &
                  &        + colamt(n, k,3) * abso3a(j) 
                  small_taur(n, k,j + offset_ng) = colmol(n, k) * rayl(j)
               enddo
            else
               !$acc loop seq
               do j = 1, ng25
                  !taug(ipt, k,ns25+j, jj) = colamt(ipt, k,3, jj) * abso3b(j) 
                  !taur(ipt, k,ns25+j, jj) = colmol(ipt, k, jj) * rayl(j)
                  small_taug(n, k,j + offset_ng) = colamt(n, k,3) * abso3b(j) 
                  small_taur(n, k,j + offset_ng) = colmol(n, k) * rayl(j)
               enddo
            end if
         end do
      end do

      return
!...................................
      end subroutine taumol25
!-----------------------------------


!-----------------------------------
      subroutine taumol26                                                 &
!...................................
     &     ( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
      &       async_id, fulljj, id0, id1, blockjj, ng00, &
      &       map_nday_ipt, map_nday_jj, nday_length, max_nday_length, jjoffset, &
      &       small_sfluxzen, small_taug, small_taur, np, offset_ng)
     
!  ------------------------------------------------------------------  !
!     band 26:  22650-29000 cm-1 (low - nothing; high - nothing)       !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb26
      
      integer, dimension(nday_length) :: map_nday_ipt, map_nday_jj
      integer :: ng00, np, offset_ng, nday_length, jjoffset, max_nday_length
      real (kind=kind_phys), dimension(nday_length, nlay, ng00) :: small_taur, small_taug
      real (kind=kind_phys), dimension(nday_length, ng00) :: small_sfluxzen
      integer :: nlay, laytrop(max_nday_length), ix, nday(fulljj), blockjj
      integer, dimension(max_nday_length, nlay) :: indfor, indself,          &
     &       jp, jt, jt1
      real (kind=kind_phys), dimension(max_nday_length, nlay) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac
      real (kind=kind_phys), dimension(max_nday_length, nlay, maxgas) :: colamt
      integer, dimension(nday_length, nlay,nblow:nbhgh) :: id0, id1
!  ---  locals:
      integer :: ipt, jj, async_id, fulljj
      integer :: j, k, n

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.
      !$acc parallel loop collapse(3) async(async_id)
      do j = 1, ng26
         do k = 1, nlay
            do n = 1, nday_length
               !taug(ipt, k,ns26+j, jj) = f_zero
               !taur(ipt, k,ns26+j, jj) = colmol(ipt, k, jj) * rayl(j) 
               small_taug(n, k,j + offset_ng) = f_zero
               small_taur(n, k,j + offset_ng) = colmol(n, k) * rayl(j) 
            enddo
         end do
      end do

      return
!...................................
      end subroutine taumol26
!-----------------------------------


!-----------------------------------
      subroutine taumol27                                                 &
!...................................
     &     ( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
      &       async_id, fulljj, id0, id1, blockjj, ng00, &
      &       map_nday_ipt, map_nday_jj, nday_length, max_nday_length, jjoffset, &
      &       small_sfluxzen, small_taug, small_taur, np, offset_ng)
     
!  ------------------------------------------------------------------  !
!     band 27:  29000-38000 cm-1 (low - o3; high - o3)                 !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb27
      
      integer, dimension(nday_length) :: map_nday_ipt, map_nday_jj
      integer :: ng00, np, offset_ng, nday_length, jjoffset, max_nday_length
      real (kind=kind_phys), dimension(nday_length, nlay, ng00) :: small_taur, small_taug
      real (kind=kind_phys), dimension(nday_length, ng00) :: small_sfluxzen
      integer :: nlay, laytrop(max_nday_length), ix, nday(fulljj), blockjj
      integer, dimension(max_nday_length, nlay) :: indfor, indself,          &
     &       jp, jt, jt1
      real (kind=kind_phys), dimension(max_nday_length, nlay) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac
      real (kind=kind_phys), dimension(max_nday_length, nlay, maxgas) :: colamt
      integer, dimension(nday_length, nlay,nblow:nbhgh) :: id0, id1
!  ---  locals:
      integer :: ipt, jj, async_id, fulljj
      integer :: ind01, ind02, ind11, ind12
      integer :: j, k, n
      real (kind=kind_phys) :: abs01, abs02, abs11, abs12

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.
      
      !$acc parallel loop collapse(2) private(ind01, ind02, ind11, ind12, &
      !$acc&         abs01, abs02, abs11, abs12) async(async_id)
      do k = 1, nlay
         do n = 1, nday_length
            ind01 = id0(n, k,27) + 1
            ind02 = ind01 + 1
            ind11 = id1(n, k,27) + 1
            ind12 = ind11 + 1
               
            !$acc loop seq
            do j = 1, ng27
               if (k .le. laytrop(n)) then
                  abs01 = absa(ind01,j)
                  abs02 = absa(ind02,j)
                  abs11 = absa(ind11,j)
                  abs12 = absa(ind12,j)
               else
                  abs01 = absb(ind01,j)
                  abs02 = absb(ind02,j)
                  abs11 = absb(ind11,j)
                  abs12 = absb(ind12,j)
               end if
               !taug(ipt, k,ns27+j, jj) = colamt(ipt, k,3, jj)                                  &
               !&        * ( fac00(ipt, k, jj)*abs01 + fac10(ipt, k, jj)*abs02       &
               !&        +   fac01(ipt, k, jj)*abs11 + fac11(ipt, k, jj)*abs12 )
               !taur(ipt, k,ns27+j, jj) = colmol(ipt, k, jj) * rayl(j)
               small_taug(n, k,j + offset_ng) = colamt(n, k,3)                                  &
               &        * ( fac00(n, k)*abs01 + fac10(n, k)*abs02       &
               &        +   fac01(n, k)*abs11 + fac11(n, k)*abs12 )
               small_taur(n, k,j + offset_ng) = colmol(n, k) * rayl(j)
            enddo
         end do
      end do

      return
!...................................
      end subroutine taumol27
!-----------------------------------


!-----------------------------------
      subroutine taumol28                                                 &
!...................................
     &     ( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
      &       async_id, fulljj, id0, id1, blockjj, ng00, &
      &       map_nday_ipt, map_nday_jj, nday_length, max_nday_length, jjoffset, &
      &       small_sfluxzen, small_taug, small_taur, np, offset_ng)
     
!  ------------------------------------------------------------------  !
!     band 28:  38000-50000 cm-1 (low - o3,o2; high - o3,o2)           !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb28
      
      integer, dimension(nday_length) :: map_nday_ipt, map_nday_jj
      integer :: ng00, np, offset_ng, nday_length, jjoffset, max_nday_length
      real (kind=kind_phys), dimension(nday_length, nlay, ng00) :: small_taur, small_taug
      real (kind=kind_phys), dimension(nday_length, ng00) :: small_sfluxzen
      integer :: nlay, laytrop(max_nday_length), ix, nday(fulljj), blockjj
      integer, dimension(max_nday_length, nlay) :: indfor, indself,          &
     &       jp, jt, jt1
      real (kind=kind_phys), dimension(max_nday_length, nlay) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac
      real (kind=kind_phys), dimension(max_nday_length, nlay, maxgas) :: colamt
      integer, dimension(nday_length, nlay,nblow:nbhgh) :: id0, id1
!  ---  locals:
      integer :: ipt, jj, async_id, fulljj
      real (kind=kind_phys) :: speccomb, specmult, tauray, fs, fs1,     &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111

      integer :: ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14
      integer :: j, js, k, n

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.

      !$acc parallel loop gang collapse(2) private(speccomb, specmult, js, fs, fs1, fac000, &
      !$acc&     fac010, fac100, fac110, fac001, fac011, fac101, fac111, &
      !$acc&     ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14, tauray) async(async_id)
      do k = 1, nlay
         do n = 1, nday_length
            tauray = colmol(n, k) * rayl
            if (k .le. laytrop(n)) then
               speccomb = colamt(n, k,3) + strrat(28)*colamt(n, k,6)
               specmult = 8.0 * min(oneminus, colamt(n, k,3) / speccomb)

               js = 1 + int(specmult)
               fs = mod(specmult, f_one)
               fs1= f_one - fs
               fac000 = fs1 * fac00(n, k)
               fac010 = fs1 * fac10(n, k)
               fac100 = fs  * fac00(n, k)
               fac110 = fs  * fac10(n, k)
               fac001 = fs1 * fac01(n, k)
               fac011 = fs1 * fac11(n, k)
               fac101 = fs  * fac01(n, k)
               fac111 = fs  * fac11(n, k)

               ind01 = id0(n, k,28) + js
               ind02 = ind01 + 1
               ind03 = ind01 + 9
               ind04 = ind01 + 10
               ind11 = id1(n, k,28) + js
               ind12 = ind11 + 1
               ind13 = ind11 + 9
               ind14 = ind11 + 10
               
               !$acc loop seq
               do j = 1, ng28
                  !taug(ipt, k,ns28+j, jj) = speccomb                                     &
                  !&        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
                  !&        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
                  !&        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
                  !&        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )
                  !taur(ipt, k,ns28+j, jj) = tauray
                  small_taug(n, k,j + offset_ng) = speccomb                                     &
                  &        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
                  &        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
                  &        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
                  &        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )
                  small_taur(n, k,j + offset_ng) = tauray
               enddo
            else
               speccomb = colamt(n, k,3) + strrat(28)*colamt(n, k,6)
               specmult = 4.0 * min(oneminus, colamt(n, k,3) / speccomb)

               js = 1 + int(specmult)
               fs = mod(specmult, f_one)
               fs1= f_one - fs
               fac000 = fs1 * fac00(n, k)
               fac010 = fs1 * fac10(n, k)
               fac100 = fs  * fac00(n, k)
               fac110 = fs  * fac10(n, k)
               fac001 = fs1 * fac01(n, k)
               fac011 = fs1 * fac11(n, k)
               fac101 = fs  * fac01(n, k)
               fac111 = fs  * fac11(n, k)

               ind01 = id0(n, k,28) + js
               ind02 = ind01 + 1
               ind03 = ind01 + 5
               ind04 = ind01 + 6
               ind11 = id1(n, k,28) + js
               ind12 = ind11 + 1
               ind13 = ind11 + 5
               ind14 = ind11 + 6
               
               !$acc loop seq
               do j = 1, ng28
                  !taug(ipt, k,ns28+j, jj) = speccomb                                     &
                  !&        * ( fac000 * absb(ind01,j) + fac100 * absb(ind02,j)       &
                  !&        +   fac010 * absb(ind03,j) + fac110 * absb(ind04,j)       &
                  !&        +   fac001 * absb(ind11,j) + fac101 * absb(ind12,j)       &
                  !&        +   fac011 * absb(ind13,j) + fac111 * absb(ind14,j) )
                  !taur(ipt, k,ns28+j, jj) = tauray
                  small_taug(n, k,j + offset_ng) = speccomb                                     &
                  &        * ( fac000 * absb(ind01,j) + fac100 * absb(ind02,j)       &
                  &        +   fac010 * absb(ind03,j) + fac110 * absb(ind04,j)       &
                  &        +   fac001 * absb(ind11,j) + fac101 * absb(ind12,j)       &
                  &        +   fac011 * absb(ind13,j) + fac111 * absb(ind14,j) )
                  small_taur(n, k,j + offset_ng) = tauray
               enddo
            end if
         end do
      end do

      return
!...................................
      end subroutine taumol28
!-----------------------------------


!-----------------------------------
      subroutine taumol29                                                 &
!...................................
     &     ( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
      &       async_id, fulljj, id0, id1, blockjj, ng00, &
      &       map_nday_ipt, map_nday_jj, nday_length, max_nday_length, jjoffset, &
      &       small_sfluxzen, small_taug, small_taur, np, offset_ng)
     
!  ------------------------------------------------------------------  !
!     band 29:  820-2600 cm-1 (low - h2o; high - co2)                  !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb29

      integer, dimension(nday_length) :: map_nday_ipt, map_nday_jj
      integer :: ng00, np, offset_ng, nday_length, jjoffset, max_nday_length
      real (kind=kind_phys), dimension(nday_length, nlay, ng00) :: small_taur, small_taug
      real (kind=kind_phys), dimension(nday_length, ng00) :: small_sfluxzen
      integer :: nlay, laytrop(max_nday_length), ix, nday(fulljj), blockjj
      integer, dimension(max_nday_length, nlay) :: indfor, indself,          &
     &       jp, jt, jt1
      real (kind=kind_phys), dimension(max_nday_length, nlay) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac
      real (kind=kind_phys), dimension(max_nday_length, nlay, maxgas) :: colamt
      integer, dimension(nday_length, nlay,nblow:nbhgh) :: id0, id1
!  ---  locals:
      integer :: ipt, jj, async_id, fulljj
      real (kind=kind_phys) :: tauray

      integer :: ind01, ind02, ind11, ind12
      integer :: inds, indf, indsp, indfp, j, k, n
      

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.
      
      !$acc parallel loop collapse(2) private(ind01, ind02, ind11, ind12, inds, indf, &
      !$acc&     indsp, indfp, tauray) async(async_id)
      do k = 1, nlay
         do n = 1, nday_length
            tauray = colmol(n, k) * rayl
            if (k .le. laytrop(n)) then
               ind01 = id0(n, k,29) + 1
               ind02 = ind01 + 1
               ind11 = id1(n, k,29) + 1
               ind12 = ind11 + 1

               inds = indself(n, k)
               indf = indfor (n, k)
               indsp= inds + 1
               indfp= indf + 1
               
               !$acc loop seq
               do j = 1, ng29
                  !taug(ipt, k,ns29+j, jj) = colamt(ipt, k,1, jj)                                  &
                  !&        * ( (fac00(ipt, k, jj)*absa(ind01,j) + fac10(ipt, k, jj)*absa(ind02,j)      &
                  !&        +    fac01(ipt, k, jj)*absa(ind11,j) + fac11(ipt, k, jj)*absa(ind12,j) )    &
                  !&        +  selffac(ipt, k, jj) * (selfref(inds,j) + selffrac(ipt, k, jj)            &
                  !&        *  (selfref(indsp,j) - selfref(inds,j)))                  &
                  !&        +  forfac(ipt, k, jj) * (forref(indf,j) + forfrac(ipt, k, jj)               &
                  !&        *  (forref(indfp,j) - forref(indf,j))))                   &
                  !&        +  colamt(ipt, k,2, jj) * absco2(j)
                  !taur(ipt, k,ns29+j, jj) = tauray
                  small_taug(n, k,j + offset_ng) = colamt(n, k,1)                                  &
                  &        * ( (fac00(n, k)*absa(ind01,j) + fac10(n, k)*absa(ind02,j)      &
                  &        +    fac01(n, k)*absa(ind11,j) + fac11(n, k)*absa(ind12,j) )    &
                  &        +  selffac(n, k) * (selfref(inds,j) + selffrac(n, k)            &
                  &        *  (selfref(indsp,j) - selfref(inds,j)))                  &
                  &        +  forfac(n, k) * (forref(indf,j) + forfrac(n, k)               &
                  &        *  (forref(indfp,j) - forref(indf,j))))                   &
                  &        +  colamt(n, k,2) * absco2(j)
                  small_taur(n, k,j + offset_ng) = tauray
               enddo
            else
               ind01 = id0(n, k,29) + 1
               ind02 = ind01 + 1
               ind11 = id1(n, k,29) + 1
               ind12 = ind11 + 1
               
               !$acc loop seq
               do j = 1, ng29
                  !taug(ipt, k,ns29+j, jj) = colamt(ipt, k,2, jj)                                  &
                  !&        * ( fac00(ipt, k, jj)*absb(ind01,j) + fac10(ipt, k, jj)*absb(ind02,j)       &
                  !&        +   fac01(ipt, k, jj)*absb(ind11,j) + fac11(ipt, k, jj)*absb(ind12,j) )     &
                  !&        + colamt(ipt, k,1, jj) * absh2o(j) 
                  !taur(ipt, k,ns29+j, jj) = tauray
                  small_taug(n, k,j + offset_ng) = colamt(n, k,2)                                  &
                  &        * ( fac00(n, k)*absb(ind01,j) + fac10(n, k)*absb(ind02,j)       &
                  &        +   fac01(n, k)*absb(ind11,j) + fac11(n, k)*absb(ind12,j) )     &
                  &        + colamt(n, k,1) * absh2o(j) 
                  small_taur(n, k,j + offset_ng) = tauray
               enddo
            end if
         end do
      end do

      return
!...................................
      end subroutine taumol29
!-----------------------------------


!
!........................................!
      end module module_radsw_main_gpu       !
!========================================!

