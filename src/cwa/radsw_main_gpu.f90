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

      use module_radsw_parameters
      use mersenne_twister, only : random_setseed, random_number,       &
     &                             random_stat
      use module_radsw_ref, only : preflog, tref
      use module_radsw_sflux
      use param,             only : my
      use index,             only : jlistnum
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
      logical :: lflxprf= .false.
      logical :: lfdncmp= .false.

!  ---  those data will be set up only once by "rswinit"

      real (kind=kind_phys) :: exp_tbl(0:ntbmx)

!  ...  heatfac is the factor for heating rates
!       (in k/day, or k/sec set by subroutine 'rswinit')

      real (kind=kind_phys) :: heatfac

!  ---  the following variables are used for sub-column cloud scheme

      integer, parameter :: ipsdsw0 = 1          ! initial permutation seed

!  ---  public accessable subprograms

      public swrad_gpu, rswinit_gpu, cldprop, setcoef, taumol, spcvrtc


! =================
      contains
! =================


!-----------------------------------
      subroutine swrad_gpu                                                  &
!...................................

!  ---  inputs:
     &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                      &
     &       clouds,icseed,aerosols,sfcalb,                             &
     &       cosz,solcon,nday,idxday,                                   &
     &       myim, nlay, nlp1, lprnt, myrank, ix,                           &
     &       nf_clds, nf_vgas, nf_albd, nf_aesw, async_id, fulljj, blocks, smalljj, &
!  ---  outputs:
     &       hswc,topflx,sfcflx                                         &
!! ---  optional:
     &,      hsw0,hswb,flxprf,fdncmp                                    &
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
      integer, intent(in) :: myim(fulljj), nlay, nlp1, nday(fulljj), myrank, ix, &
         nf_clds, nf_vgas, nf_albd, nf_aesw, fulljj, blocks, smalljj

      integer, dimension(ix, fulljj), intent(in) :: idxday, icseed

      logical, intent(in) :: lprnt

      real (kind=kind_phys), dimension(ix,nlp1, fulljj), intent(in) ::        &
     &       plvl, tlvl
      real (kind=kind_phys), dimension(ix,nlay, fulljj), intent(in) ::        &
     &       plyr, tlyr, qlyr, olyr
      real (kind=kind_phys), dimension(ix,nf_albd, fulljj),    intent(in) :: sfcalb

      real (kind=kind_phys), dimension(ix,nlay,nf_vgas, fulljj),intent(in):: gasvmr
      real (kind=kind_phys), dimension(ix,nlay,nf_clds, fulljj),intent(in):: clouds
      real (kind=kind_phys), dimension(ix,nlay,nbdsw,nf_aesw, fulljj),intent(in)::  &
     &       aerosols

      real (kind=kind_phys), intent(in) :: cosz(ix, fulljj), solcon

!  ---  outputs:
      real (kind=kind_phys), dimension(ix,nlay, fulljj), intent(out) :: hswc

      real (kind=kind_phys),    dimension(ix, fulljj, 3), intent(out) :: topflx
      real (kind=kind_phys),    dimension(ix, fulljj, 4), intent(out) :: sfcflx

!! ---  optional outputs:
      real (kind=kind_phys), dimension(ix,nlay,nbdsw, fulljj), optional,      &
     &       intent(out) :: hswb

      real (kind=kind_phys), dimension(ix,nlay, fulljj),       optional,      &
     &       intent(out) :: hsw0
      real (kind=kind_phys),    dimension(ix,nlp1, fulljj, 4),       optional,      &
     &       intent(out) :: flxprf
      real (kind=kind_phys),    dimension(ix, fulljj, 6),            optional,      &
     &       intent(out) :: fdncmp

!  ---  locals:
      real (kind=kind_phys), dimension(ix, nlay,ngptsw, smalljj) :: taug, taur
      real (kind=kind_phys), dimension(:,:,:,:), allocatable :: cldfmc
      real (kind=kind_phys), dimension(ix, nlp1,nbdsw, smalljj):: fxupc, fxdnc,      &
     &       fxup0, fxdn0

      real (kind=kind_phys), dimension(ix, nlay,nbdsw, smalljj)  ::                  &
     &       tauae, ssaae, asyae, taucw, ssacw, asycw

      real (kind=kind_phys), dimension(ix, ngptsw, smalljj) :: sfluxzen

      real (kind=kind_phys), dimension(ix, nlay, smalljj)   :: cldfrc,     delp,     &
     &       pavel, tavel, coldry, colmol, h2ovmr, o3vmr, temcol,       &
     &       cliqp, reliq, cicep, reice, cdat1, cdat2, cdat3, cdat4,    &
     &       cfrac, fac00, fac01, fac10, fac11, forfac, forfrac,        &
     &       selffac, selffrac, rfdelp

      real (kind=kind_phys), dimension(ix, nlp1, smalljj) :: fnet, flxdc, flxuc,     &
     &       flxd0, flxu0
      real (kind=kind_phys), dimension(ix, nlp1, nbdsw, smalljj) :: fnetm

      real (kind=kind_phys), dimension(ix, 2, smalljj) :: albbm, albdf, sfbmc,       &
     &       sfbm0, sfdfc, sfdf0

      real (kind=kind_phys) :: cosz1, sntz1(ix, smalljj), tem0, tem1, tem2, s0fac
      real (kind=kind_phys), dimension(ix, smalljj) :: ssolar, zcf0, zcf1, ftoau0, &
         ftoauc, ftoadc, fsfcu0, fsfcuc, fsfcd0, fsfcdc, suvbfc, suvbf0

!  ---  column amount of absorbing gases:
!       (:,m) m = 1-h2o, 2-co2, 3-o3, 4-n2o, 5-ch4, 6-o2, 7-co
      real (kind=kind_phys) ::  colamt(ix, nlay,maxgas, smalljj)

      integer, dimension(ix, smalljj) :: ipseed
      integer, dimension(ix, nlay, smalljj) :: indfor, indself, jp, jt, jt1

      integer :: i, ib, ipt, j1, k, kk, laytrop(ix, smalljj), mb, jj
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
      integer :: async_id, jf, jb,  jjoffset, jbs, jbe
!
!===> ... begin here
!
      if (isubcsw > 0) allocate(cldfmc(ix, nlay,ngptsw, smalljj))
      lhswb  = present ( hswb )
      lhsw0  = present ( hsw0 )
      lflxprf= present ( flxprf )
      lfdncmp= present ( fdncmp )
 
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
      !$acc data create(taug, taur, fxupc, fxdnc, fxup0, fxdn0, tauae, ssaae, &
      !$acc&     asyae, taucw, ssacw, asycw, sfluxzen, cldfrc, delp, pavel, &
      !$acc&     tavel, coldry, colmol, h2ovmr, o3vmr, temcol, cliqp, reliq, &
      !$acc&     cicep, reice, cdat1, cdat2, cdat3, cdat4, cfrac, fac00, fac01, &
      !$acc&     fac10, fac11, forfac, forfrac, selffac, selffrac, rfdelp, fnet, &
      !$acc&     flxdc, flxuc, flxd0, flxu0, fnetm, albbm, albdf, sfbmc, sfbm0, &
      !$acc&     sfdfc, sfdf0, sntz1, ssolar, zcf0, zcf1, ftoau0, ftoauc, ftoadc, &
      !$acc&     fsfcu0, fsfcuc, fsfcd0, fsfcdc, suvbfc, suvbf0, colamt, ipseed, &
      !$acc&     indfor, indself, jp, jt, jt1, laytrop) async(async_id)
      do jb = 1, blocks
      jjoffset = (jb-1)*smalljj
      jbs = jjoffset+1
      jbe = jjoffset+smalljj
      !$acc parallel loop gang collapse(2) private(jf) async(async_id)
      do jj = 1, smalljj
         do k = 1, nlay
            jf = jjoffset+jj
            !$acc loop vector
            do i = 1, myim(jf)
               hswc(i, k, jf) = f_zero
            end do
         end do
      end do
      !$acc parallel loop collapse(2) private(jf) async(async_id)
      do jj = 1, smalljj
         do i = 1, ix
            jf = jjoffset+jj
            if (i .le. myim(jf)) then
               topflx(i, jf, 1) = f_zero
               topflx(i, jf, 2) = f_zero
               topflx(i, jf, 3) = f_zero
               sfcflx(i, jf, 1) = f_zero
               sfcflx(i, jf, 2) = f_zero
               sfcflx(i, jf, 3) = f_zero
               sfcflx(i, jf, 4) = f_zero
            end if
         end do
      end do

         !! --- ...  initial optional outputs
      if ( lflxprf ) then
         !$acc parallel loop gang collapse(2) private(jf) async(async_id)
         do jj = 1, smalljj
            do k = 1, nlp1
               jf = jjoffset+jj
               !$acc loop vector
               do i = 1, myim(jf)
                  flxprf(i, k, jf, 1) = f_zero
                  flxprf(i, k, jf, 2) = f_zero
                  flxprf(i, k, jf, 3) = f_zero
                  flxprf(i, k, jf, 4) = f_zero
               end do
            end do
         end do
      endif

      if ( lfdncmp ) then
         !$acc parallel loop collapse(2) private(jf) async(async_id)
         do jj = 1, smalljj
            do i = 1, ix
               jf = jjoffset+jj
               if (i .le. myim(jf)) then
                  fdncmp(i, jf, 1) = f_zero
                  fdncmp(i, jf, 2) = f_zero
                  fdncmp(i, jf, 3) = f_zero
                  fdncmp(i, jf, 4) = f_zero
                  fdncmp(i, jf, 5) = f_zero
                  fdncmp(i, jf, 6) = f_zero
               end if
            end do
         end do
      endif

      if ( lhsw0 ) then
         !$acc parallel loop gang collapse(2) private(jf) async(async_id)
         do jj = 1, smalljj
            do k = 1, nlay
               jf = jjoffset+jj
               !$acc loop vector
               do i = 1, myim(jf)
                  hsw0(i,k, jf) = f_zero
               end do
            end do
         end do
      endif

      if ( lhswb ) then
         !$acc parallel loop gang collapse(3) private(jf) async(async_id)
         do jj = 1, smalljj
            do ib = 1, nbdsw
               do k = 1, nlay
                  jf = jjoffset+jj
                  !$acc loop vector
                  do i = 1, myim(jf)
                     hswb(i, k, ib, jf) = f_zero
                  end do
               end do
            end do
         end do
      endif

      !  --- ...  change random number seed value for each radiation invocation

      if     ( isubcsw == 1 ) then     ! advance prescribed permutation seed
         !$acc parallel loop collapse(2) private(jf) async(async_id)
         do jj = 1, smalljj
            do i = 1, ix
               jf = jjoffset+jj
               if (i .le. myim(jf)) then
                  ipseed(i, jj) = ipsdsw0 + i
               end if
            enddo
         end do
      elseif ( isubcsw == 2 ) then     ! use input array of permutaion seeds
         !$acc parallel loop collapse(2) private(jf) async(async_id)
         do jj = 1, smalljj
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
      !$acc parallel loop collapse(2) private(j1, jf) async(async_id)
      do jj = 1, smalljj
         do ipt = 1, ix
            jf = jjoffset+jj
            if (ipt .le. nday(jf)) then ! lab_do_ipt
               j1 = idxday(ipt, jf)

               sntz1(ipt, jj)  = f_one / cosz(j1, jf)
               ssolar(ipt, jj) = s0fac * cosz(j1, jf)

               !  --- ...  surface albedo: bm,df - dir,dif;  1,2 - nir,uvv
               albbm(ipt, 1, jj) = sfcalb(j1,1, jf)
               albdf(ipt, 1, jj) = sfcalb(j1,2, jf)
               albbm(ipt, 2, jj) = sfcalb(j1,3, jf)
               albdf(ipt, 2, jj) = sfcalb(j1,4, jf)
            end if
         end do
      end do

            !  --- ...  prepare atmospheric profile for use in rrtm
            !           the vertical index of internal array is from surface to top

      if (ivflip == 0) then       ! input from toa to sfc

         tem1 = 100.0 * con_g
         tem2 = 1.0e-20 * 1.0e3 * con_avgd
         !$acc parallel loop gang collapse(2) private(jf) async(async_id)
         do jj = 1, smalljj
            do k = 1, nlay
               jf = jjoffset+jj
               !$acc loop vector private(j1, kk, tem0)
               do ipt = 1, nday(jf) ! lab_do_ipt
                  j1 = idxday(ipt, jf)
                  kk = nlp1 - k
                  pavel(ipt, k, jj) = plyr(j1,kk, jf)
                  tavel(ipt, k, jj) = tlyr(j1,kk, jf)
                  delp (ipt, k, jj) = plvl(j1,kk+1, jf) - plvl(j1,kk, jf)

                  !  --- ...  set absorber amount
                  !test use
                  !           h2ovmr(k)= max(f_zero,qlyr(j1,kk)*amdw)                     ! input mass mixing ratio
                  !           h2ovmr(k)= max(f_zero,qlyr(j1,kk))                          ! input vol mixing ratio
                  !           o3vmr (k)= max(f_zero,olyr(j1,kk))                          ! input vol mixing ratio
                  !ncep model use
                  h2ovmr(ipt, k, jj)= max(f_zero,qlyr(j1,kk, jf)*amdw/(f_one-qlyr(j1,kk, jf))) ! input specific humidity
                  o3vmr (ipt, k, jj)= max(f_zero,olyr(j1,kk, jf)*amdo3)                    ! input mass mixing ratio

                  tem0 = (f_one - h2ovmr(ipt, k, jj))*con_amd + h2ovmr(ipt, k, jj)*con_amw
                  coldry(ipt, k, jj) = tem2 * delp(ipt, k, jj) / &
                  (tem1*tem0*(f_one + h2ovmr(ipt, k, jj)))
                  temcol(ipt, k, jj) = 1.0e-12 * coldry(ipt, k, jj)

                  colamt(ipt, k,1, jj) = max(f_zero,    coldry(ipt, k, jj)*h2ovmr(ipt, k, jj))         ! h2o
                  colamt(ipt, k,2, jj) = max(temcol(ipt, k, jj), coldry(ipt, k, jj)*gasvmr(j1,kk,1, jf))   ! co2
                  colamt(ipt, k,3, jj) = max(f_zero,    coldry(ipt, k, jj)*o3vmr(ipt, k, jj))          ! o3
                  colmol(ipt, k, jj)   = coldry(ipt, k, jj) + colamt(ipt, k,1, jj)
               enddo
            end do
         end do

               !  --- ...  set up gas column amount, convert from volume mixing ratio
               !           to molec/cm2 based on coldry (scaled to 1.0e-20)

         if (iswrgas > 0) then
            !$acc parallel loop gang collapse(2) private(jf) async(async_id)
            do jj = 1, smalljj
               do k = 1, nlay
                  jf = jjoffset+jj
                  !$acc loop vector private(j1, kk)
                  do ipt = 1, nday(jf) ! lab_do_ipt
                     j1 = idxday(ipt, jf)
                     kk = nlp1 - k
                     colamt(ipt, k,4, jj) = max(temcol(ipt, k, jj), coldry(ipt, k, jj)*gasvmr(j1,kk,2, jf))  ! n2o
                     colamt(ipt, k,5, jj) = max(temcol(ipt, k, jj), coldry(ipt, k, jj)*gasvmr(j1,kk,3, jf))  ! ch4
                     colamt(ipt, k,6, jj) = max(temcol(ipt, k, jj), coldry(ipt, k, jj)*gasvmr(j1,kk,4, jf))  ! o2
                     !             colamt(k,7) = max(temcol(k), coldry(k)*gasvmr(j1,kk,5))  ! co - notused
                  enddo
               end do
            end do
         else
            !$acc parallel loop gang collapse(2) private(jf) async(async_id)
            do jj = 1, smalljj
               do k = 1, nlay
                  jf = jjoffset+jj
                  !$acc loop vector private(j1)
                  do ipt = 1, nday(jf) ! lab_do_ipt
                     j1 = idxday(ipt, jf)
                     colamt(ipt, k,4, jj) = temcol(ipt, k, jj)                                  ! n2o
                     colamt(ipt, k,5, jj) = temcol(ipt, k, jj)                                  ! ch4
                     colamt(ipt, k,6, jj) = temcol(ipt, k, jj)                                  ! o2
                     !             colamt(k,7) = temcol(k)                                  ! co - notused
                  enddo
               end do
            end do
         endif

               !  --- ...  set aerosol optical properties
         !$acc parallel loop gang collapse(3) private(kk, jf) async(async_id)
         do jj = 1, smalljj
            do ib = 1, nbdsw
               do k = 1, nlay
                  kk = nlp1 - k
                  jf = jjoffset+jj
                  !$acc loop vector private(j1)
                  do ipt = 1, nday(jf) ! lab_do_ipt
                     j1 = idxday(ipt, jf)
                     tauae(ipt, k,ib, jj) = aerosols(j1,kk,ib,1, jf)
                     ssaae(ipt, k,ib, jj) = aerosols(j1,kk,ib,2, jf)
                     asyae(ipt, k,ib, jj) = aerosols(j1,kk,ib,3, jf)
                  enddo
               enddo
            end do
         end do

         if (iswcliq > 0) then    ! use prognostic cloud method
            !$acc parallel loop gang collapse(2) private(kk, jf) async(async_id)
            do jj = 1, smalljj
               do k = 1, nlay
                  kk = nlp1 - k
                  jf = jjoffset+jj
                  !$acc loop vector private(j1)
                  do ipt = 1, nday(jf) ! lab_do_ipt
                     j1 = idxday(ipt, jf)
                     cfrac(ipt, k, jj) = clouds(j1,kk,1, jf)      ! cloud fraction
                     cliqp(ipt, k, jj) = clouds(j1,kk,2, jf)      ! cloud liq path
                     reliq(ipt, k, jj) = clouds(j1,kk,3, jf)      ! liq partical effctive radius
                     cicep(ipt, k, jj) = clouds(j1,kk,4, jf)      ! cloud ice path
                     reice(ipt, k, jj) = clouds(j1,kk,5, jf)      ! ice partical effctive radius
                     cdat1(ipt, k, jj) = clouds(j1,kk,6, jf)      ! cloud rain drop path
                     cdat2(ipt, k, jj) = clouds(j1,kk,7, jf)      ! rain partical effctive radius
                     cdat3(ipt, k, jj) = clouds(j1,kk,8, jf)      ! cloud snow path
                     cdat4(ipt, k, jj) = clouds(j1,kk,9, jf)      ! snow partical effctive radius
                  enddo
               end do
            end do
         else                     ! use diagnostic cloud method
            !$acc parallel loop gang collapse(2) private(kk, jf) async(async_id)
            do jj = 1, smalljj
               do k = 1, nlay
                  kk = nlp1 - k
                  jf = jjoffset+jj
                  !$acc loop vector private(j1)
                  do ipt = 1, nday(jf) ! lab_do_ipt
                     j1 = idxday(ipt, jf)
                     cfrac(ipt, k, jj) = clouds(j1,kk,1, jf)      ! cloud fraction
                     cdat1(ipt, k, jj) = clouds(j1,kk,2, jf)      ! cloud optical depth
                     cdat2(ipt, k, jj) = clouds(j1,kk,3, jf)      ! cloud single scattering albedo
                     cdat3(ipt, k, jj) = clouds(j1,kk,4, jf)      ! cloud asymmetry factor
                  enddo
               end do
            end do
         endif                    ! end if_iswcliq
      end if
      if (ivflip .ne. 0) then                        ! input from sfc to toa

         tem1 = 100.0 * con_g
         tem2 = 1.0e-20 * 1.0e3 * con_avgd
         !$acc parallel loop gang collapse(2) private(jf) async(async_id)
         do jj = 1, smalljj
            do k = 1, nlay
               jf = jjoffset+jj
               !$acc loop vector private(j1, tem0)
               do ipt = 1, nday(jf) ! lab_do_ipt
                  j1 = idxday(ipt, jf)
                  pavel(ipt, k, jj) = plyr(j1,k, jf)
                  tavel(ipt, k, jj) = tlyr(j1,k, jf)
                  delp (ipt, k, jj) = plvl(j1,k, jf) - plvl(j1,k+1, jf)

                  !  --- ...  set absorber amount
                  !test use
                  !           h2ovmr(k)= max(f_zero,qlyr(j1,k)*amdw)                    ! input mass mixing ratio
                  !           h2ovmr(k)= max(f_zero,qlyr(j1,k))                         ! input vol mixing ratio
                  !           o3vmr (k)= max(f_zero,olyr(j1,k))                         ! input vol mixing ratio
                  !ncep model use
                  h2ovmr(ipt, k, jj)= max(f_zero,qlyr(j1,k, jf)*amdw/(f_one-qlyr(j1,k, jf))) ! input specific humidity
                  o3vmr (ipt, k, jj)= max(f_zero,olyr(j1,k, jf)*amdo3)                   ! input mass mixing ratio

                  tem0 = (f_one - h2ovmr(ipt, k, jj))*con_amd + h2ovmr(ipt, k, jj)*con_amw
                  coldry(ipt, k, jj) = tem2 * delp(ipt, k, jj) / (tem1*tem0*(f_one + h2ovmr(ipt, k, jj)))
                  temcol(ipt, k, jj) = 1.0e-12 * coldry(ipt, k, jj)

                  colamt(ipt, k,1, jj) = max(f_zero,    coldry(ipt, k, jj)*h2ovmr(ipt, k, jj))         ! h2o
                  colamt(ipt, k,2, jj) = max(temcol(ipt, k, jj), coldry(ipt, k, jj)*gasvmr(j1,k,1, jf))    ! co2
                  colamt(ipt, k,3, jj) = max(f_zero,    coldry(ipt, k, jj)*o3vmr(ipt, k, jj))          ! o3
                  colmol(ipt, k, jj)   = coldry(ipt, k, jj) + colamt(ipt, k,1, jj)
               enddo
            end do
         end do

         !  --- ...  set up gas column amount, convert from volume mixing ratio
         !           to molec/cm2 based on coldry (scaled to 1.0e-20)

         if (iswrgas > 0) then
            !$acc parallel loop gang collapse(2) private(jf) async(async_id)
            do jj = 1, smalljj
               do k = 1, nlay
                  jf = jjoffset+jj
                  !$acc loop vector private(j1)
                  do ipt = 1, nday(jf) ! lab_do_ipt
                     j1 = idxday(ipt, jf)
                     colamt(ipt, k,4, jj) = max(temcol(ipt, k, jj), coldry(ipt, k, jj)*gasvmr(j1,k,2, jf))  ! n2o
                     colamt(ipt, k,5, jj) = max(temcol(ipt, k, jj), coldry(ipt, k, jj)*gasvmr(j1,k,3, jf))  ! ch4
                     colamt(ipt, k,6, jj) = max(temcol(ipt, k, jj), coldry(ipt, k, jj)*gasvmr(j1,k,4, jf))  ! o2
                     !             colamt(k,7) = max(temcol(k), coldry(k)*gasvmr(j1,k,5))  ! co - notused
                  enddo
               end do
            end do
         else
            !$acc parallel loop gang collapse(2) private(jf) async(async_id)
            do jj = 1, smalljj
               do k = 1, nlay
                  jf = jjoffset+jj
                  !$acc loop vector private(j1)
                  do ipt = 1, nday(jf) ! lab_do_ipt
                     j1 = idxday(ipt, jf)
                     colamt(ipt, k,4, jj) = temcol(ipt, k, jj)                                 ! n2o
                     colamt(ipt, k,5, jj) = temcol(ipt, k, jj)                                 ! ch4
                     colamt(ipt, k,6, jj) = temcol(ipt, k, jj)                                 ! o2
                     !             colamt(k,7) = temcol(k)                                 ! co - notused
                  enddo
               end do
            end do
         endif

               !  --- ...  set aerosol optical properties
         !$acc parallel loop gang collapse(3) private(jf) async(async_id)
         do jj = 1, smalljj
            do ib = 1, nbdsw
               do k = 1, nlay
                  jf = jjoffset+jj
                  !$acc loop vector private(j1)
                  do ipt = 1, nday(jf) ! lab_do_ipt
                     j1 = idxday(ipt, jf)
                     tauae(ipt, k,ib, jj) = aerosols(j1,k,ib,1, jf)
                     ssaae(ipt, k,ib, jj) = aerosols(j1,k,ib,2, jf)
                     asyae(ipt, k,ib, jj) = aerosols(j1,k,ib,3, jf)
                  enddo
               enddo
            end do
         end do

         if (iswcliq > 0) then    ! use prognostic cloud method
            !$acc parallel loop gang collapse(2) private(jf) async(async_id)
            do jj = 1, smalljj
               do k = 1, nlay
                  jf = jjoffset+jj
                  !$acc loop vector private(j1)
                  do ipt = 1, nday(jf) ! lab_do_ipt
                     j1 = idxday(ipt, jf)
                     cfrac(ipt, k, jj) = clouds(j1,k,1, jf)       ! cloud fraction
                     cliqp(ipt, k, jj) = clouds(j1,k,2, jf)       ! cloud liq path
                     reliq(ipt, k, jj) = clouds(j1,k,3, jf)       ! liq partical effctive radius
                     cicep(ipt, k, jj) = clouds(j1,k,4, jf)       ! cloud ice path
                     reice(ipt, k, jj) = clouds(j1,k,5, jf)       ! ice partical effctive radius
                     cdat1(ipt, k, jj) = clouds(j1,k,6, jf)       ! cloud rain drop path
                     cdat2(ipt, k, jj) = clouds(j1,k,7, jf)       ! rain partical effctive radius
                     cdat3(ipt, k, jj) = clouds(j1,k,8, jf)       ! cloud snow path
                     cdat4(ipt, k, jj) = clouds(j1,k,9, jf)       ! snow partical effctive radius
                  enddo
               end do
            end do
         else                     ! use diagnostic cloud method
            !$acc parallel loop gang collapse(2) private(jf) async(async_id)
            do jj = 1, smalljj
               do k = 1, nlay
                  jf = jjoffset+jj
                  !$acc loop vector private(j1)
                  do ipt = 1, nday(jf) ! lab_do_ipt
                     j1 = idxday(ipt, jf)
                     cfrac(ipt, k, jj) = clouds(j1,k,1, jf)       ! cloud fraction
                     cdat1(ipt, k, jj) = clouds(j1,k,2, jf)       ! cloud optical depth
                     cdat2(ipt, k, jj) = clouds(j1,k,3, jf)       ! cloud single scattering albedo
                     cdat3(ipt, k, jj) = clouds(j1,k,4, jf)       ! cloud asymmetry factor
                  enddo
               end do
            end do
         endif                    ! end if_iswcliq
      endif                       ! if_ivflip

            !  --- ...  compute fractions of clear sky view
      !$acc parallel loop collapse(2) private(jf) async(async_id)
      do jj = 1, smalljj
         do ipt = 1, ix
            jf = jjoffset+jj
            if (ipt .le. nday(jf)) then ! lab_do_ipt
               zcf0(ipt, jj)   = f_one
               zcf1(ipt, jj)   = f_one
            end if
         end do
      end do
      if (iovrsw == 0) then                    ! random overlapping
         !$acc parallel loop collapse(2) private(jf) async(async_id)
         do jj = 1, smalljj
            do ipt = 1, ix
               jf = jjoffset+jj
               if (ipt .le. nday(jf)) then ! lab_do_ipt
                  !$acc loop seq
                  do k = 1, nlay
                     zcf0(ipt, jj) = zcf0(ipt, jj) * (f_one - cfrac(ipt, k, jj))
                  enddo
               end if
            end do
         end do
      else if (iovrsw == 1) then               ! max/ran overlapping
         !$acc parallel loop collapse(2) private(jf) async(async_id)
         do jj = 1, smalljj
            do ipt = 1, ix
               jf = jjoffset+jj
               if (ipt .le. nday(jf)) then ! lab_do_ipt
                  !$acc loop seq
                  do k = 1, nlay
                     if (cfrac(ipt, k, jj) > ftiny) then                ! cloudy layer
                        zcf1(ipt, jj) = min ( zcf1(ipt, jj), f_one-cfrac(ipt, k, jj) )
                     elseif (zcf1(ipt, jj) < f_one) then                ! clear layer
                        zcf0(ipt, jj) = zcf0(ipt, jj) * zcf1(ipt, jj)
                        zcf1(ipt, jj) = f_one
                     endif
                  enddo
               end if
            end do
         end do
         !$acc parallel loop collapse(2) private(jf) async(async_id)
         do jj = 1, smalljj
            do ipt = 1, ix
               jf = jjoffset+jj
               if (ipt .le. nday(jf)) then ! lab_do_ipt
                  zcf0(ipt, jj) = zcf0(ipt, jj) * zcf1(ipt, jj)
               end if
            end do
         end do
      else if (iovrsw == 2) then               ! maximum overlapping
         !$acc parallel loop collapse(2) private(jf) async(async_id)
         do jj = 1, smalljj
            do ipt = 1, ix
               jf = jjoffset+jj
               if (ipt .le. nday(jf)) then ! lab_do_ipt
                  !$acc loop seq
                  do k = 1, nlay
                     zcf0(ipt, jj) = min ( zcf0(ipt, jj), f_one-cfrac(ipt, k, jj) )
                  enddo
               end if
            end do
         end do
      endif
      !$acc parallel loop collapse(2) private(jf) async(async_id)
      do jj = 1, smalljj
         do ipt = 1, ix
            jf = jjoffset+jj
            if (ipt .le. nday(jf)) then ! lab_do_ipt
               if (zcf0(ipt, jj) <= ftiny) zcf0(ipt, jj) = f_zero
               if (zcf0(ipt, jj) > oneminus) zcf0(ipt, jj) = f_one
               zcf1(ipt, jj) = f_one - zcf0(ipt, jj)
            end if
         end do
      end do

            !  --- ...  compute cloud optical properties
      !call nvtxStartRange("sw_cldprop")
      call cldprop                                                  &
         !  ---  inputs:
         &     ( cfrac,cliqp,reliq,cicep,reice,cdat1,cdat2,cdat3,cdat4,     &
         &       zcf1, nlay, ipseed, ix, nday(jbs:jbe), idxday(:,jbs:jbe), async_id, smalljj,                                    &
         !  ---  outputs:
         &       taucw, ssacw, asycw, cldfrc, cldfmc                        &
         &     )
      !call nvtxEndRange
      if (isubcsw > 0) then
         !$acc parallel loop gang collapse(2) private(jf) async(async_id)
         do jj = 1, smalljj
            do k = 1, nlay
               jf = jjoffset+jj
               !$acc loop vector private(j1)
               do ipt = 1, nday(jf) ! lab_do_ipt
                  j1 = idxday(ipt, jf)
                  if (zcf1(ipt, jj) <= 0.) then     ! clear sky column
                     !$acc loop seq
                     do ib = 1, ngptsw
                        cldfmc(ipt, k,ib, jj)= f_zero
                     end do
                  endif   ! end if_zcf1_block
               enddo
            end do
         end do
      end if
      !call nvtxStartRange("sw_setcoef")
      call setcoef                                                    &
         !  ---  inputs:
         &     ( pavel,tavel,h2ovmr, nlay,nlp1, ix, nday(jbs:jbe), async_id, smalljj,                              &
         !  ---  outputs:
         &       laytrop,jp,jt,jt1,fac00,fac01,fac10,fac11,                 &
         &       selffac,selffrac,indself,forfac,forfrac,indfor             &
         &     )
      !call nvtxEndRange

         !  --- ...  calculate optical depths for gaseous absorption and rayleigh
         !           scattering
      !call nvtxStartRange("sw_taumol")
      call taumol                                                     &
         !  ---  inputs:
         &     ( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
         &       forfac,forfrac,indfor,selffac,selffrac,indself, &
                 nlay, ix, nday(jbs:jbe), async_id, smalljj,     &
         !  ---  outputs:
         &       sfluxzen, taug, taur                                       &
         &     )
      !call nvtxEndRange

         !  --- ...  call the 2-stream radiation transfer model
      if ( isubcsw <= 0 ) then     ! use standard cloud scheme
         !call nvtxStartRange("sw_spcvrtc")
         call spcvrtc                                                  &
            !  ---  inputs:
            &     ( ssolar,cosz(:,jbs:jbe),sntz1,albbm, &
                     albdf,sfluxzen,cldfrc,            &
            &       zcf1,zcf0,taug,taur,tauae, &
                     ssaae,asyae,taucw,ssacw,asycw,   &
            &       nlay, nlp1, ix, idxday(:,jbs:jbe), nday(jbs:jbe), async_id, smalljj,                                    &
            !  ---  outputs:
            &       fxupc,fxdnc,fxup0,fxdn0,                                   &
            &       ftoauc,ftoau0,ftoadc, &
                     fsfcuc,fsfcu0,fsfcdc,fsfcd0,          &
            &       sfbmc,sfdfc,sfbm0,sfdf0,suvbfc,suvbf0                      &
            &     )
         !call nvtxEndRange

      else                         ! use mcica cloud scheme
         do jj = 1, smalljj
            jf = jjoffset+jj
            do ipt = 1, nday(jf) ! lab_do_ipt
               j1 = idxday(ipt, jf)
               do k = 1, nlay
                  do i = 1, nbdsw
                     tauae_im(k, i) = tauae(ipt, k, i, jj)
                     ssaae_im(k, i) = ssaae(ipt, k, i, jj)
                     asyae_im(k, i) = asyae(ipt, k, i, jj)
                     taucw_im(k, i) = taucw(ipt, k, i, jj)
                     ssacw_im(k, i) = ssacw(ipt, k, i, jj)
                     asycw_im(k, i) = asycw(ipt, k, i, jj)
                  end do
                  do i = 1, ngptsw
                     cldfmc_im(k, i) = cldfmc(ipt, k, i, jj)
                     taug_im(k, i) = taug(ipt, k, i, jj)
                     taur_im(k, i) = taur(ipt, k, i, jj)
                  end do
               end do
               do k = 1, 2
                  albbm_im(k) = albbm(ipt, k, jj)
                  albdf_im(k) = albdf(ipt, k, jj)
                  sfbmc_im(k) = sfbmc(ipt, k, jj)
                  sfdfc_im(k) = sfdfc(ipt, k, jj)
                  sfbm0_im(k) = sfbm0(ipt, k, jj)
                  sfdf0_im(k) = sfdf0(ipt, k, jj)
               end do
               do k = 1, nlp1
                  do i = 1, nbdsw
                     fxupc_im(k, i) = fxupc(ipt, k, i, jj)
                     fxdnc_im(k, i) = fxdnc(ipt, k, i, jj)
                     fxup0_im(k, i) = fxup0(ipt, k, i, jj)
                     fxdn0_im(k, i) = fxdn0(ipt, k, i, jj)
                  end do
               end do
               do i = 1, ngptsw
                  sfluxzen_im(i) = sfluxzen(ipt, i, jj)
               end do

               call spcvrtm                                                  &
                  !  ---  inputs:
                  &     ( ssolar(ipt, jj),cosz(j1, jf),sntz1(ipt, jj),albbm_im, &
                          albdf_im,sfluxzen_im,cldfmc_im,            &
                  &       zcf1(ipt, jj),zcf0(ipt, jj),taug_im,taur_im,tauae_im, &
                          ssaae_im,asyae_im,taucw_im,ssacw_im,asycw_im,   &
                  &       nlay, nlp1,                                                &
                  !  ---  outputs:
                  &       fxupc_im,fxdnc_im,fxup0_im,fxdn0_im,                                   &
                  &       ftoauc(ipt, jj),ftoau0(ipt, jj),ftoadc(ipt, jj), &
                          fsfcuc(ipt, jj),fsfcu0(ipt, jj),fsfcdc(ipt, jj),fsfcd0(ipt, jj),          &
                  &       sfbmc_im,sfdfc_im,sfbm0_im,sfdf0_im,suvbfc(ipt, jj),suvbf0(ipt, jj)                      &
                  &     )
               do k = 1, nlay
                  do i = 1, nbdsw
                     tauae(ipt, k, i, jj) = tauae_im(k, i)
                     ssaae(ipt, k, i, jj) = ssaae_im(k, i)
                     asyae(ipt, k, i, jj) = asyae_im(k, i)
                     taucw(ipt, k, i, jj) = taucw_im(k, i)
                     ssacw(ipt, k, i, jj) = ssacw_im(k, i)
                     asycw(ipt, k, i, jj) = asycw_im(k, i)
                  end do
                  do i = 1, ngptsw
                     cldfmc(ipt, k, i, jj) = cldfmc_im(k, i)
                     taug(ipt, k, i, jj) = taug_im(k, i)
                     taur(ipt, k, i, jj) = taur_im(k, i)
                  end do
               end do
               do k = 1, 2
                  albbm(ipt, k, jj) = albbm_im(k)
                  albdf(ipt, k, jj) = albdf_im(k)
                  sfbmc(ipt, k, jj) = sfbmc_im(k)
                  sfdfc(ipt, k, jj) = sfdfc_im(k)
                  sfbm0(ipt, k, jj) = sfbm0_im(k)
                  sfdf0(ipt, k, jj) = sfdf0_im(k)
               end do
               do k = 1, nlp1
                  do i = 1, nbdsw
                     fxupc(ipt, k, i, jj) = fxupc_im(k, i)
                     fxdnc(ipt, k, i, jj) = fxdnc_im(k, i)
                     fxup0(ipt, k, i, jj) = fxup0_im(k, i)
                     fxdn0(ipt, k, i, jj) = fxdn0_im(k, i)
                  end do
               end do
               do i = 1, ngptsw
                  sfluxzen(ipt, i, jj) = sfluxzen_im(i)
               end do
            end do
         end do
      endif

         !  --- ...  sum up total spectral fluxes for total-sky
      !$acc parallel loop collapse(2) private(jf) async(async_id)
      do jj = 1, smalljj
         do k = 1, nlp1
            jf = jjoffset+jj
            !$acc loop vector
            do ipt = 1, nday(jf) ! lab_do_ipt
               flxuc(ipt, k, jj) = f_zero
               flxdc(ipt, k, jj) = f_zero
               !$acc loop seq
               do ib = 1, nbdsw
                  flxuc(ipt, k, jj) = flxuc(ipt, k, jj) + fxupc(ipt, k,ib, jj)
                  flxdc(ipt, k, jj) = flxdc(ipt, k, jj) + fxdnc(ipt, k,ib, jj)
               enddo
            enddo
         end do
      end do

         !! --- ...  optional clear sky fluxes

      if ( lhsw0 .or. lflxprf ) then
         !$acc parallel loop collapse(2) private(jf) async(async_id)
         do jj = 1, smalljj
            do k = 1, nlp1
               jf = jjoffset+jj
               !$acc loop vector
               do ipt = 1, nday(jf) ! lab_do_ipt
                  flxu0(ipt, k, jj) = f_zero
                  flxd0(ipt, k, jj) = f_zero
                  !$acc loop seq
                  do ib = 1, nbdsw
                     flxu0(ipt, k, jj) = flxu0(ipt, k, jj) + fxup0(ipt, k,ib, jj)
                     flxd0(ipt, k, jj) = flxd0(ipt, k, jj) + fxdn0(ipt, k,ib, jj)
                  enddo
               enddo
            end do
         end do
      endif

         !  --- ...  prepare for final outputs
      !$acc parallel loop collapse(2) private(jf) async(async_id)
      do jj = 1, smalljj
         do k = 1, nlay
            jf = jjoffset+jj
            !$acc loop vector
            do ipt = 1, nday(jf) ! lab_do_ipt
               rfdelp(ipt, k, jj) = heatfac / delp(ipt, k, jj)
            enddo
         end do
      end do

      if ( lfdncmp ) then
         !$acc parallel loop collapse(2) private(j1, jf) async(async_id)
         do jj = 1, smalljj
            do ipt = 1, ix
               jf = jjoffset+jj
               if (ipt .le. nday(jf)) then ! lab_do_ipt
                  j1 = idxday(ipt, jf)
                  !! --- ...  optional uv-b surface downward flux
                  fdncmp(j1, jf, 2) = suvbf0(ipt, jj)
                  fdncmp(j1, jf, 1) = suvbfc(ipt, jj)

                  !! --- ...  optional beam and diffuse sfc fluxes
                  fdncmp(j1, jf, 3) = sfbmc(ipt, 1, jj)
                  fdncmp(j1, jf, 4) = sfdfc(ipt, 1, jj)
                  fdncmp(j1, jf, 5) = sfbmc(ipt, 2, jj)
                  fdncmp(j1, jf, 6) = sfdfc(ipt, 2, jj)
               end if
            end do
         end do
      endif    ! end if_lfdncmp

         !  --- ...  toa and sfc fluxes
      !$acc parallel loop collapse(2) private(j1, jf) async(async_id)
      do jj = 1, smalljj
         do ipt = 1, ix
            jf = jjoffset+jj
            if (ipt .le. nday(jf)) then ! lab_do_ipt
               j1 = idxday(ipt, jf)
               topflx(j1, jf, 1) = ftoauc(ipt, jj)
               topflx(j1, jf, 2) = ftoadc(ipt, jj)
               topflx(j1, jf, 3) = ftoau0(ipt, jj)

               sfcflx(j1, jf, 1) = fsfcuc(ipt, jj)
               sfcflx(j1, jf, 2) = fsfcdc(ipt, jj)
               sfcflx(j1, jf, 3) = fsfcu0(ipt, jj)
               sfcflx(j1, jf, 4) = fsfcd0(ipt, jj)
            end if
         end do
      end do
      
      if (ivflip == 0) then       ! output from toa to sfc

               !  --- ...  compute heating rates
         !$acc parallel loop collapse(2) private(j1, kk, jf) async(async_id)
         do jj = 1, smalljj
            do ipt = 1, ix
               jf = jjoffset+jj
               if (ipt .le. nday(jf)) then ! lab_do_ipt
                  j1 = idxday(ipt, jf)
                  fnet(ipt, 1, jj) = flxdc(ipt, 1, jj) - flxuc(ipt, 1, jj)
                  !$acc loop seq
                  do k = 2, nlp1
                     kk = nlp1 - k + 1
                     fnet(ipt, k, jj) = flxdc(ipt, k, jj) - flxuc(ipt, k, jj)
                     hswc(j1, kk, jf) = (fnet(ipt, k, jj)-fnet(ipt, k-1, jj)) * rfdelp(ipt, k-1, jj)
                  enddo
               end if
            end do
         end do

               !! --- ...  optional flux profiles

         if ( lflxprf ) then
            !$acc parallel loop collapse(2) private(kk, jf) async(async_id)
            do jj = 1, smalljj
               do k = 1, nlp1
                  jf = jjoffset+jj
                  kk = nlp1 - k + 1
                  !$acc loop vector private(j1)
                  do ipt = 1, nday(jf) ! lab_do_ipt
                     j1 = idxday(ipt, jf)
                     flxprf(j1,kk, jf, 1) = flxuc(ipt, k, jj)
                     flxprf(j1,kk, jf, 2) = flxdc(ipt, k, jj)
                     flxprf(j1,kk, jf, 3) = flxu0(ipt, k, jj)
                     flxprf(j1,kk, jf, 4) = flxd0(ipt, k, jj)
                  enddo
               end do
            end do
         endif

               !! --- ...  optional clear sky heating rates

         if ( lhsw0 ) then
            !$acc parallel loop collapse(2) private(j1, kk, jf) async(async_id)
            do jj = 1, smalljj
               do ipt = 1, ix 
                  jf = jjoffset+jj
                  if (ipt .le. nday(jf)) then ! lab_do_ipt
                     j1 = idxday(ipt, jf)
                     fnet(ipt, 1, jj) = flxd0(ipt, 1, jj) - flxu0(ipt, 1, jj)
                     !$acc loop seq
                     do k = 2, nlp1
                        kk = nlp1 - k + 1
                        fnet(ipt, k, jj) = flxd0(ipt, k, jj) - flxu0(ipt, k, jj)
                        hsw0(j1,kk, jf) = (fnet(ipt, k, jj)-fnet(ipt, k-1, jj)) * rfdelp(ipt, k-1, jj)
                     enddo
                  end if
               end do
            end do
         endif

               !! --- ...  optional spectral band heating rates

         if ( lhswb ) then
            !$acc parallel loop collapse(2) private(jf) async(async_id)
            do jj = 1, smalljj
               do mb = 1, nbdsw
                  jf = jjoffset+jj
                  !$acc loop vector private(j1, kk)
                  do ipt = 1, nday(jf) ! lab_do_ipt
                     j1 = idxday(ipt, jf)
                     fnetm(ipt, 1, mb, jj) = fxdnc(ipt, 1,mb, jj) - fxupc(ipt, 1,mb, jj)
                     !$acc loop seq
                     do k = 2, nlp1
                        kk = nlp1 - k + 1
                        fnetm(ipt, k, mb, jj) = fxdnc(ipt, k,mb, jj) - fxupc(ipt, k,mb, jj)
                        hswb(j1,kk,mb, jf) = (fnetm(ipt, k, mb, jj) - fnetm(ipt, k-1, mb, jj)) * rfdelp(ipt, k-1, jj)
                     enddo
                  enddo
               end do
            end do
         endif

      else                        ! output from sfc to toa

               !  --- ...  compute heating rates
         !$acc parallel loop collapse(2) private(j1, jf) async(async_id)
         do jj = 1, smalljj
            do ipt = 1, ix
               jf = jjoffset+jj
               if (ipt .le. nday(jf)) then ! lab_do_ipt
                  j1 = idxday(ipt, jf)
                  fnet(ipt, 1, jj) = flxdc(ipt, 1, jj) - flxuc(ipt, 1, jj)
                  !$acc loop seq
                  do k = 2, nlp1
                     fnet(ipt, k, jj) = flxdc(ipt, k, jj) - flxuc(ipt, k, jj)
                     hswc(j1,k-1, jf) = (fnet(ipt, k, jj)-fnet(ipt, k-1, jj)) * rfdelp(ipt, k-1, jj)
                  enddo
               end if
            end do
         end do

               !! --- ...  optional flux profiles

         if ( lflxprf ) then
            !$acc parallel loop gang collapse(2) private(jf) async(async_id)
            do jj = 1, smalljj
               do k = 1, nlp1
                  jf = jjoffset+jj
                  !$acc loop vector private(j1)
                  do ipt = 1, nday(jf) ! lab_do_ipt
                     j1 = idxday(ipt, jf)
                     flxprf(j1,k, jf, 1) = flxuc(ipt, k, jj)
                     flxprf(j1,k, jf, 2) = flxdc(ipt, k, jj)
                     flxprf(j1,k, jf, 3) = flxu0(ipt, k, jj)
                     flxprf(j1,k, jf, 4) = flxd0(ipt, k, jj)
                  enddo
               end do
            end do
         endif

               !! --- ...  optional clear sky heating rates

         if ( lhsw0 ) then
            !$acc parallel loop collapse(2) private(j1, jf) async(async_id)
            do jj = 1, smalljj
               do ipt = 1, ix 
                  jf = jjoffset+jj
                  if (ipt .le. nday(jf)) then ! lab_do_ipt
                     j1 = idxday(ipt, jf)
                     fnet(ipt, 1, jj) = flxd0(ipt, 1, jj) - flxu0(ipt, 1, jj)
                     !$acc loop seq
                     do k = 2, nlp1
                        fnet(ipt, k, jj) = flxd0(ipt, k, jj) - flxu0(ipt, k, jj)
                        hsw0(j1,k-1, jf) = (fnet(ipt, k, jj)-fnet(ipt, k-1, jj)) * rfdelp(ipt, k-1, jj)
                     enddo
                  end if
               end do
            end do
         endif

               !! --- ...  optional spectral band heating rates

         if ( lhswb ) then
            !$acc parallel loop collapse(2) private(jf) async(async_id)
            do jj = 1, smalljj
               do mb = 1, nbdsw
                  jf = jjoffset+jj
                  !$acc loop vector private(j1)
                  do ipt = 1, nday(jf) ! lab_do_ipt
                     j1 = idxday(ipt, jf)
                     fnetm(ipt, 1, mb, jj) = fxdnc(ipt, 1,mb, jj) - fxupc(ipt, 1,mb, jj)
                     !$acc loop seq
                     do k = 1, nlay
                        fnetm(ipt, k+1, mb, jj) = fxdnc(ipt, k+1,mb, jj) - fxupc(ipt, k+1,mb, jj)
                        hswb(j1,k,mb, jf) = (fnetm(ipt, k+1, mb, jj) - fnetm(ipt, k, mb, jj)) * rfdelp(ipt, k, jj)
                     enddo
                  enddo
               end do
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

      return
!...................................
      end subroutine rswinit_gpu
!-----------------------------------


!-----------------------------------
      subroutine cldprop                                                &
!...................................
!  ---  inputs:
     &     ( cfrac,cliqp,reliq,cicep,reice,cdat1,cdat2,cdat3,cdat4,     &
     &       cf1, nlay, ipseed, ix, nday, idxday, async_id, fulljj,                                           &
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
      integer, intent(in) :: nlay, fulljj, ipseed(ix, fulljj), ix, &
         nday(fulljj), idxday(ix, fulljj)
      real (kind=kind_phys), intent(in) :: cf1(ix, fulljj)

      real (kind=kind_phys), dimension(ix, nlay, fulljj), intent(in) :: cliqp,      &
     &       reliq, cicep, reice, cdat1, cdat2, cdat3, cdat4, cfrac

!  ---  outputs:
      real (kind=kind_phys), dimension(ix, nlay,ngptsw, fulljj), intent(out) ::     &
     &       cldfmc
      real (kind=kind_phys), dimension(ix, nlay,nbdsw, fulljj),  intent(out) ::     &
     &       taucw, ssacw, asycw
      real (kind=kind_phys), dimension(ix, nlay, fulljj), intent(out) :: cldfrc

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
      integer :: ia, ib, ig, jb, k, index, ipt, jj, j1, async_id
      real (kind=kind_phys), dimension(:), allocatable       :: cldf_im
      logical, allocatable, dimension(:,:) :: lcloudy_im

      !
      !===> ...  begin here
      !
      !$acc parallel loop gang collapse(3) async(async_id) 
      do jj = 1, fulljj
         do ib = 1, nbdsw
            do k = 1, nlay
               !$acc loop vector 
               do ipt = 1, nday(jj) ! lab_do_ipt
                  taucw (ipt, k,ib, jj) = f_zero
                  asycw (ipt, k,ib, jj) = f_zero
                  if (cf1(ipt, jj) > f_zero) then
                     ssacw (ipt, k,ib, jj) = f_one
                  else
                     ssacw (ipt, k,ib, jj) = f_zero
                  end if
               enddo
            enddo
         end do
      end do


               !  --- ...  compute cloud radiative properties for a cloudy column
      if (iswcliq > 0) then ! lab_if_iswcliq
         !$acc parallel loop gang collapse(2) async(async_id)
         do jj = 1, fulljj
            do k = 1, nlay ! lab_do_k
               !$acc loop vector private(cldran, cldsnw, refsnw, dgesnw, tauran, tausnw, &
               !$acc&     ssaran, ssasnw, asyran, asysnw, cldliq, cldice, refliq, refice, &
               !$acc&     tauliq, ssaliq, asyliq, factor, index, fint, extcoliq, ssacoliq, &
               !$acc&     asycoliq, tauice, ssaice, asyice, ia, extcoice, ssacoice, asycoice, &
               !$acc&     dgeice)
               do ipt = 1, nday(jj) ! lab_do_ipt
                  if (cf1(ipt, jj) > f_zero) then     ! cloudy sky column
                     if (cfrac(ipt, k, jj) > ftiny) then ! lab_if_cld

                        !  --- ...  optical properties for rain and snow
                        cldran = cdat1(ipt, k, jj)
                        !           refran = cdat2(k)
                        cldsnw = cdat3(ipt, k, jj)
                        refsnw = cdat4(ipt, k, jj)
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
                        !$acc loop seq
                        do ib = nblow, nbhgh
                           ssaran(ib) = tauran * (f_one - b0r(ib))
                           ssasnw(ib) = tausnw * (f_one - (b0s(ib)+b1s(ib)*dgesnw))
                           asyran(ib) = ssaran(ib) * c0r(ib)
                           asysnw(ib) = ssasnw(ib) * c0s(ib)
                        enddo

                        cldliq = cliqp(ipt, k, jj)
                        cldice = cicep(ipt, k, jj)
                        refliq = reliq(ipt, k, jj)
                        refice = reice(ipt, k, jj)

                        !  --- ...  calculation of absorption coefficients due to water clouds.

                        if ( cldliq <= f_zero ) then
                           !$acc loop seq
                           do ib = nblow, nbhgh
                              tauliq(ib) = f_zero
                              ssaliq(ib) = f_zero
                              asyliq(ib) = f_zero
                           enddo
                        else
                           if ( iswcliq == 1 ) then
                              factor = refliq - 1.5
                              index  = max( 1, min( 57, int( factor ) ))
                              fint   = factor - float(index)
                              !$acc loop seq
                              do ib = nblow, nbhgh
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
                              enddo
                           endif   ! end if_iswcliq_block
                        endif   ! end if_cldliq_block

                        !  --- ...  calculation of absorption coefficients due to ice clouds.

                        if ( cldice <= f_zero ) then
                           !$acc loop seq
                           do ib = nblow, nbhgh
                              tauice(ib) = f_zero
                              ssaice(ib) = f_zero
                              asyice(ib) = f_zero
                           enddo
                        else

                           !  --- ...  ebert and curry approach for all particle sizes though somewhat
                           !           unjustified for large ice particles

                           if ( iswcice == 1 ) then
                              refice = min(130.0_kind_phys,max(13.0_kind_phys,refice))
                              !$acc loop seq
                              do ib = nblow, nbhgh
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
                              enddo

                              !  --- ...  streamer approach for ice effective radius between 5.0 and 131.0 microns

                           elseif ( iswcice == 2 ) then
                              refice = min(131.0_kind_phys,max(5.0_kind_phys,refice))

                              factor = (refice - 2.0) / 3.0
                              index  = max( 1, min( 42, int( factor ) ))
                              fint   = factor - float(index)
                              !$acc loop seq
                              do ib = nblow, nbhgh
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
                              enddo

                              !  --- ...  fu's approach for ice effective radius between 4.8 and 135 microns
                              !           (generalized effective size from 5 to 140 microns)

                           elseif ( iswcice == 3 ) then
                              dgeice = max( 5.0, min( 140.0, 1.0315*refice ))

                              factor = (dgeice - 2.0) / 3.0
                              index  = max( 1, min( 45, int( factor ) ))
                              fint   = factor - float(index)
                              !$acc loop seq
                              do ib = nblow, nbhgh
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
                              enddo

                           endif   ! end if_iswcice_block
                        endif   ! end if_cldice_block
                        !$acc loop seq
                        do ib = 1, nbdsw
                           jb = nblow + ib - 1
                           taucw(ipt, k,ib, jj) = tauliq(jb)+tauice(jb)+tauran+tausnw
                           ssacw(ipt, k,ib, jj) = ssaliq(jb)+ssaice(jb)+ssaran(jb)+ssasnw(jb)
                           asycw(ipt, k,ib, jj) = asyliq(jb)+asyice(jb)+asyran(jb)+asysnw(jb)
                        enddo

                     endif  ! lab_if_cld
                  end if
               enddo  ! lab_do_k
            end do
         end do

      else  ! lab_if_iswcliq
         !$acc parallel loop gang collapse(3) async(async_id)
         do jj = 1, fulljj
            do ib = 1, nbdsw
               do k = 1, nlay
                  !$acc loop vector
                  do ipt = 1, nday(jj) ! lab_do_ipt
                     if (cf1(ipt, jj) > f_zero) then     ! cloudy sky column
                        if (cfrac(ipt, k, jj) > ftiny) then
                           taucw(ipt, k,ib, jj) = cdat1(ipt, k, jj)
                           ssacw(ipt, k,ib, jj) = cdat1(ipt, k, jj)    * cdat2(ipt, k, jj)
                           asycw(ipt, k,ib, jj) = ssacw(ipt, k,ib, jj) * cdat3(ipt, k, jj)
                        endif
                     end if
                  enddo
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
         do jj = 1, fulljj
            do ipt = 1, nday(jj) ! lab_do_ipt
               if (cf1(ipt, jj) > f_zero) then     ! cloudy sky column
                  cldf(ipt, :, jj) = cfrac(ipt, :, jj)
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
                  j1 = idxday(ipt, jj)
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
         end do
         deallocate(cldf)
         deallocate(cldf_im)
         deallocate(lcloudy_im)
         deallocate(lcloudy)

      else                         ! non-mcica, normalize cloud
         !$acc parallel loop gang collapse(2) async(async_id)
         do jj = 1, fulljj
            do k = 1, nlay
               !$acc loop vector
               do ipt = 1, nday(jj) ! lab_do_ipt
                  if (cf1(ipt, jj) > f_zero) then     ! cloudy sky column
                     cldfrc(ipt, k, jj) = cfrac(ipt, k, jj) / cf1(ipt, jj)
                  else
                     cldfrc(ipt, k, jj) = f_zero
                  end if
               enddo
            end do
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
     &     ( pavel,tavel,h2ovmr, nlay,nlp1, ix, nday, async_id, fulljj,                             &
!  ---  outputs:
     &       laytrop,jp,jt,jt1,fac00,fac01,fac10,fac11,                 &
     &       selffac,selffrac,indself,forfac,forfrac,indfor             &
     &     )


!  ---  inputs:
      integer, intent(in) :: nlay, nlp1, ix, nday(fulljj), fulljj

      real (kind=kind_phys), dimension(:,:,:), intent(in) :: pavel, tavel,  &
     &       h2ovmr

!  ---  outputs:
      integer, dimension(ix, nlay, fulljj), intent(out) :: indself, indfor,         &
     &       jp, jt, jt1
      integer, intent(out) :: laytrop(ix, fulljj)

      real (kind=kind_phys), dimension(ix, nlay, fulljj), intent(out) :: fac00,     &
     &       fac01, fac10, fac11, selffac, selffrac, forfac, forfrac

!  ---  locals:
      real (kind=kind_phys) :: plog, fp, fp1, ft, ft1, tem1, tem2

      integer :: i, k, jp1, jj, ipt, async_id
!
!===> ... begin here
!
      !$acc parallel loop collapse(2) async(async_id)
      do jj = 1, fulljj
         do ipt = 1, ix
            if (ipt .le. nday(jj)) then
               laytrop(ipt, jj) = nlay
               !$acc loop seq
               do k = 1, nlay
                  if (log(pavel(ipt, k, jj)) > 4.56) laytrop(ipt, jj) =  k
               end do
            end if
         end do
      end do
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !$acc loop vector private(plog, jp1, fp, tem1, tem2, ft, ft1, fp1)
            do ipt = 1, nday(jj)

               forfac(ipt, k, jj) = pavel(ipt, k, jj)*stpfac / (tavel(ipt, k, jj)*(f_one + h2ovmr(ipt, k, jj)))

               !  --- ...  find the two reference pressures on either side of the
               !           layer pressure.  store them in jp and jp1.  store in fp the
               !           fraction of the difference (in ln(pressure)) between these
               !           two values that the layer pressure lies.

               plog  = log(pavel(ipt, k, jj))
               jp(ipt, k, jj) = max(1, min(58, int(36.0 - 5.0*(plog+0.04)) ))
               jp1   = jp(ipt, k, jj) + 1
               fp    = 5.0 * (preflog(jp(ipt, k, jj)) - plog)

               !  --- ...  determine, for each reference pressure (jp and jp1), which
               !          reference temperature (these are different for each reference
               !          pressure) is nearest the layer temperature but does not exceed it.
               !          store these indices in jt and jt1, resp. store in ft (resp. ft1)
               !          the fraction of the way between jt (jt1) and the next highest
               !          reference temperature that the layer temperature falls.

               tem1 = (tavel(ipt, k, jj) - tref(jp(ipt, k, jj))) / 15.0
               tem2 = (tavel(ipt, k, jj) - tref(jp1  )) / 15.0
               jt (ipt, k, jj) = max(1, min(4, int(3.0 + tem1) ))
               jt1(ipt, k, jj) = max(1, min(4, int(3.0 + tem2) ))
               ft  = tem1 - float(jt (ipt, k, jj) - 3)
               ft1 = tem2 - float(jt1(ipt, k, jj) - 3)

               !  --- ...  we have now isolated the layer ln pressure and temperature,
               !           between two reference pressures and two reference temperatures
               !           (for each reference pressure).  we multiply the pressure
               !           fraction fp with the appropriate temperature fractions to get
               !           the factors that will be needed for the interpolation that yields
               !           the optical depths (performed in routines taugbn for band n).

               fp1 = f_one - fp
               fac10(ipt, k, jj) = fp1 * ft
               fac00(ipt, k, jj) = fp1 * (f_one - ft)
               fac11(ipt, k, jj) = fp  * ft1
               fac01(ipt, k, jj) = fp  * (f_one - ft1)
            end do
         end do
      end do
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !$acc loop vector private(tem1, tem2)
            do ipt = 1, nday(jj)
               !  --- ...  if the pressure is less than ~100mb, perform a different
               !           set of species interpolations.

               if ( log(pavel(ipt, k, jj)) > 4.56 ) then

                  !  --- ...  set up factors needed to separately include the water vapor
                  !           foreign-continuum in the calculation of absorption coefficient.

                  tem1 = (332.0 - tavel(ipt, k, jj)) / 36.0
                  indfor (ipt, k, jj) = min(2, max(1, int(tem1)))
                  forfrac(ipt, k, jj) = tem1 - float(indfor(ipt, k, jj))

                  !  --- ...  set up factors needed to separately include the water vapor
                  !           self-continuum in the calculation of absorption coefficient.

                  tem2 = (tavel(ipt, k, jj) - 188.0) / 7.2
                  indself (ipt, k, jj) = min(9, max(1, int(tem2)-7))
                  selffrac(ipt, k, jj) = tem2 - float(indself(ipt, k, jj) + 7)
                  selffac (ipt, k, jj) = h2ovmr(ipt, k, jj) * forfac(ipt, k, jj)

               else

                  !  --- ...  set up factors needed to separately include the water vapor
                  !           foreign-continuum in the calculation of absorption coefficient.

                  tem1 = (tavel(ipt, k, jj) - 188.0) / 36.0
                  indfor (ipt, k, jj) = 3
                  forfrac(ipt, k, jj) = tem1 - f_one

                  indself (ipt, k, jj) = 0
                  selffrac(ipt, k, jj) = f_zero
                  selffac (ipt, k, jj) = f_zero

               endif

            enddo    ! end_do_k_loop
         end do
      end do

      return
! ..................................
      end subroutine setcoef
! ----------------------------------


!-----------------------------------
      subroutine spcvrtc                                                &
!...................................
!  ---  inputs:
     &     ( ssolar,cosz,sntz,albbm,albdf,sfluxzen,cldfrc,              &
     &       cf1,cf0,taug,taur,tauae,ssaae,asyae,taucw,ssacw,asycw,     &
     &       nlay, nlp1, ix, idxday, nday, async_id, fulljj,                                                &
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
      integer, parameter :: small_factor = 2
      integer, parameter :: small_ngptsw = int(ngptsw/small_factor)

!  ---  inputs:
      integer, intent(in) :: nlay, nlp1, ix, idxday(ix, fulljj), nday(fulljj), fulljj

      real (kind=kind_phys), dimension(ix, nlay,ngptsw, fulljj), intent(in) ::      &
     &       taug, taur
      real (kind=kind_phys), dimension(ix, nlay,nbdsw, fulljj),  intent(in) ::      &
     &       taucw, ssacw, asycw, tauae, ssaae, asyae

      real (kind=kind_phys), dimension(ix, ngptsw, fulljj), intent(in) :: sfluxzen
      real (kind=kind_phys), dimension(ix, nlay, fulljj),   intent(in) :: cldfrc

      real (kind=kind_phys), dimension(ix, 2, fulljj),  intent(in) :: albbm, albdf

      real (kind=kind_phys), dimension(ix, fulljj), intent(in) :: cosz, sntz, &
         cf1, cf0, ssolar

!  ---  outputs:
      real (kind=kind_phys), dimension(ix, nlp1,nbdsw, fulljj), intent(out) ::      &
     &       fxupc, fxdnc, fxup0, fxdn0

      real (kind=kind_phys), dimension(ix, 2, fulljj), intent(out) :: sfbmc, sfdfc, &
     &       sfbm0, sfdf0

      real (kind=kind_phys), dimension(ix, fulljj), intent(out) :: suvbfc, suvbf0, ftoadc,     &
     &       ftoauc, ftoau0, fsfcuc, fsfcu0, fsfcdc, fsfcd0

!  ---  locals:
      real (kind=kind_phys), dimension(ix, nlp1, small_ngptsw, fulljj) :: zrefb, zrefd, ztrab,    &
     &       ztrad, ztdbt, zldbt, zfua, zfda

      real (kind=kind_phys) :: ztau1, zssa1, zasy1, ztau0, zssa0,       &
     &       zasy0, zasy3, zssaw, zasyw, zgam1, zgam2, zgam3, zgam4,    &
     &       zc0, zc1, za1, za2, &
     &       zrk, zrk2, zrp, zrp1, zrm1, zb1, zb2, &
     &       zrpp, zrkg1, zrkg3, zrkg4, zexp1, zexm1, zexp2, zexm2,     &
     &       zexp3, zexp4, zden1, ze1r45, ftind, zrefb1,        &
     &       zrefd1, ztrab1, ztrad1, zr1, zr2, zr3, zr4, zr5, ztdbt0r,    &
     &       zt1, zt2, zt3, zf1, zf2, zldbt0, zfu, zfd, zrupbr, zrupbr1, zrupdr, zrupdr1

      real (kind=kind_phys), dimension(ix, small_ngptsw, fulljj) :: ztdbt0, zsolar, zfd0

      integer :: ib, ibd, jb, jg, k, kp, itind, ipt, jj, j1, j2, jg2
      real (kind=kind_phys), dimension(ix, nlp1, small_ngptsw, fulljj) :: zrdnd, ztdn
      integer :: async_id

!
!===> ...  begin here
!

      !$acc data create(zrefb, zrefd, ztrab, ztrad, ztdbt, zldbt, zfd0, &
      !$acc&     ztdbt0, zsolar, zrdnd, ztdn, zfda, zfua) async(async_id)
!  --- ... initialization of output fluxes
      !$acc parallel loop gang collapse(3) async(async_id)
      do jj = 1, fulljj
         do ib = 1, nbdsw
            do k = 1, nlp1
               !$acc loop vector
               do ipt = 1, nday(jj) ! lab_do_ipt
                  fxdnc(ipt, k,ib, jj) = f_zero
                  fxupc(ipt, k,ib, jj) = f_zero
                  fxdn0(ipt, k,ib, jj) = f_zero
                  fxup0(ipt, k,ib, jj) = f_zero
               enddo
            enddo
         end do
      end do
      !$acc parallel loop collapse(2) async(async_id)
      do jj = 1, fulljj
         do ipt = 1, ix
            if (ipt .le. nday(jj)) then ! lab_do_ipt
               ftoadc(ipt, jj) = f_zero
               ftoauc(ipt, jj) = f_zero
               ftoau0(ipt, jj) = f_zero
               fsfcuc(ipt, jj) = f_zero
               fsfcu0(ipt, jj) = f_zero
               fsfcdc(ipt, jj) = f_zero
               fsfcd0(ipt, jj) = f_zero

               !! --- ...  uv-b surface downward fluxes
               suvbfc(ipt, jj)  = f_zero
               suvbf0(ipt, jj)  = f_zero

               !! --- ...  output surface flux components
               sfbmc(ipt, 1, jj) = f_zero
               sfbmc(ipt, 2, jj) = f_zero
               sfdfc(ipt, 1, jj) = f_zero
               sfdfc(ipt, 2, jj) = f_zero
               sfbm0(ipt, 1, jj) = f_zero
               sfbm0(ipt, 2, jj) = f_zero
               sfdf0(ipt, 1, jj) = f_zero
               sfdf0(ipt, 2, jj) = f_zero
            end if
         end do
      end do

            !  --- ...  loop over all g-points in each band
   do j2 = 1, small_factor
      !$acc parallel loop gang collapse(2) private(jg2, jb, ib, ibd) async(async_id)
      do jj = 1, fulljj
         do jg = 1, small_ngptsw ! lab_do_jg
            jg2 = jg + (j2-1)*small_ngptsw
            jb = ngb(jg2)
            ib = jb + 1 - nblow
            ibd = idxsfc(jb)
            !$acc loop vector private(j1, kp, ztau0, zssa0, zasy0, zssaw, zasyw, &
            !$acc&     za1, za2, ztau1, zssa1, zasy1, zasy3, zgam1, zgam2, zgam3, &
            !$acc&     zgam4, zb1, zb2, ftind, itind, zrk, zrk2, zrp, zrp1, zrm1, &
            !$acc&     zrpp, zrkg1, zrkg3, zrkg4, zr1, zr2, zr3, zr4, zr5, zt1, &
            !$acc&     zt2, zt3, zexm1, zexp1, zexm2, zexp2, ze1r45, zden1, zexp3, zexp4)
            do ipt = 1, nday(jj) ! lab_do_ipt

               zsolar(ipt, jg, jj) = ssolar(ipt, jj) * sfluxzen(ipt, jg2, jj)

               !  --- ...  set up toa direct beam and surface values (beam and diff)

               ztdbt(ipt, nlp1, jg, jj) = f_one

               zldbt(ipt, 1, jg, jj) = f_zero
               if (ibd /= 0) then
                  zrefb(ipt, 1, jg, jj) = albbm(ipt, ibd, jj)
                  zrefd(ipt, 1, jg, jj) = albdf(ipt, ibd, jj)
               else
                  zrefb(ipt, 1, jg, jj) = 0.5 * (albbm(ipt, 1, jj) + albbm(ipt, 2, jj))
                  zrefd(ipt, 1, jg, jj) = 0.5 * (albdf(ipt, 1, jj) + albdf(ipt, 2, jj))
               endif
               ztrab(ipt, 1, jg, jj) = f_zero
               ztrad(ipt, 1, jg, jj) = f_zero

            !  --- ...  compute clear-sky optical parameters, layer reflectance and transmittance
               j1 = idxday(ipt, jj)
               ztdbt0r = f_one
               !$acc loop seq
               do k = nlay, 1, -1
                  kp = k + 1

                  ztau0 = max( ftiny, taur(ipt, k,jg2, jj)+taug(ipt, k,jg2, jj)+tauae(ipt, k,ib, jj) )
                  zssa0 = taur(ipt, k,jg2, jj) + tauae(ipt, k,ib, jj)*ssaae(ipt, k,ib, jj)
                  zasy0 = asyae(ipt, k,ib, jj)*ssaae(ipt, k,ib, jj)*tauae(ipt, k,ib, jj)
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

                  !  --- ...  compute homogeneous reflectance and transmittance

                  if ( zssaw >= zcrit ) then    ! for conservative scattering
                     za1 = zgam1 * cosz(j1, jj) - zgam3
                     za2 = zgam1 * ztau1

                     !  --- ...  use exponential lookup table for transmittance, or expansion
                     !           of exponential for low optical depth

                     zb1 = min ( ztau1*sntz(ipt, jj) , 500.0 )
                     if ( zb1 <= od_lo ) then
                        zb2 = f_one - zb1 + 0.5*zb1*zb1
                     else
                        ftind = zb1 / (bpade + zb1)
                        itind = ftind*ntbmx + 0.5
                        zb2 = exp_tbl(itind)
                     endif

                     !      ...  collimated beam
                     zrefb(ipt, kp, jg, jj) = max(f_zero, min(f_one,                          &
                     &                  (za2 - za1*(f_one - zb2))/(f_one + za2) ))
                     ztrab(ipt, kp, jg, jj) = max(f_zero, min(f_one, f_one-zrefb(ipt, kp, jg, jj) ))

                     !      ...  isotropic incidence
                     zrefd(ipt, kp, jg, jj) = max(f_zero, min(f_one, za2/(f_one + za2) ))
                     ztrad(ipt, kp, jg, jj) = max(f_zero, min(f_one, f_one-zrefd(ipt, kp, jg, jj) ))

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

                     zb2 = min ( sntz(ipt, jj)*ztau1, 500.0 )
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
                        zrefb(ipt, kp, jg, jj) = eps1
                        ztrab(ipt, kp, jg, jj) = zexm2
                     else
                        zden1 = zssa1 / ze1r45
                        zrefb(ipt, kp, jg, jj) = max(f_zero, min(f_one,                        &
                        &                    (zr1*zexp1 - zr2*zexm1 - zr3*zexm2)*zden1 ))
                        ztrab(ipt, kp, jg, jj) = max(f_zero, min(f_one, zexm2*(f_one           &
                        &                  - (zt1*zexp1 - zt2*zexm1 - zt3*zexp2)*zden1) ))
                     endif

                     !      ...  diffuse beam
                     zden1 = zr4 / (ze1r45 * zrkg1)
                     zrefd(ipt, kp, jg, jj) = max(f_zero, min(f_one,                          &
                     &                  zgam2*(zexp1 - zexm1)*zden1 ))
                     ztrad(ipt, kp, jg, jj) = max(f_zero, min(f_one, zrk2*zden1 ))
                  endif    ! end if_zssaw_block

                  !  --- ...  direct beam transmittance. use exponential lookup table
                  !           for transmittance, or expansion of exponential for low
                  !           optical depth

                  zr1 = ztau1 * sntz(ipt, jj)
                  if ( zr1 <= od_lo ) then
                     zexp3 = f_one - zr1 + 0.5*zr1*zr1
                  else
                     ftind = zr1 / (bpade + zr1)
                     itind = max(0, min(ntbmx, int(0.5+ntbmx*ftind) ))
                     zexp3 = exp_tbl(itind)
                  endif

                  ztdbt(ipt, k, jg, jj)  = zexp3 * ztdbt(ipt, kp, jg, jj)
                  zldbt(ipt, kp, jg, jj) = zexp3

                  !  --- ...  pre-delta-scaling clear and cloudy direct beam transmittance
                  !           (must use 'orig', unscaled cloud optical depth)

                  zr1 = ztau0 * sntz(ipt, jj)
                  if ( zr1 <= od_lo ) then
                     zexp4 = f_one - zr1 + 0.5*zr1*zr1
                  else
                     ftind = zr1 / (bpade + zr1)
                     itind = max(0, min(ntbmx, int(0.5+ntbmx*ftind) ))
                     zexp4 = exp_tbl(itind)
                  endif

                  ztdbt0r = zexp4 * ztdbt0r
               enddo    ! end do_k_loop
               ztdbt0(ipt, jg, jj) = ztdbt0r
            end do
         end do
      end do
               
               ! call swflux
               !  --- ...  link lowest layer with surface
      !$acc parallel loop gang collapse(2) private(jg2, jb, ib, ibd) async(async_id)
      do jj = 1, fulljj
         do jg = 1, small_ngptsw ! lab_do_jg
            jg2 = jg + (j2-1)*small_ngptsw
            jb = ngb(jg2)
            ib = jb + 1 - nblow
            ibd = idxsfc(jb)
            !$acc loop vector private(kp, zden1, zfu, zfd, zrupbr, zrupbr1, zrupdr, zrupdr1)
            do ipt = 1, nday(jj) ! lab_do_ipt
               !  --- ...  upper boundary conditions

               ztdn (ipt, nlp1, jg, jj) = f_one
               zrdnd(ipt, nlp1, jg, jj) = f_zero
               ztdn (ipt, nlay, jg, jj) = ztrab(ipt, nlp1, jg, jj)
               zrdnd(ipt, nlay, jg, jj) = zrefd(ipt, nlp1, jg, jj)

               !  --- ...  pass from top to bottom
               !$acc loop seq
               do k = nlay, 2, -1
                  zden1 = f_one / (f_one - zrefd(ipt, k, jg, jj)*zrdnd(ipt, k, jg, jj))
                  ztdn (ipt, k-1, jg, jj) = ztdbt(ipt, k, jg, jj)*ztrab(ipt, k, jg, jj) &
                                    + ( ztrad(ipt, k, jg, jj) *                 &
                  &                 ( (ztdn(ipt, k, jg, jj) - ztdbt(ipt, k, jg, jj)) + ztdbt(ipt, k, jg, jj) *              &
                  &                 zrefb(ipt, k, jg, jj)*zrdnd(ipt, k, jg, jj) )) * zden1
                  zrdnd(ipt, k-1, jg, jj) = zrefd(ipt, k, jg, jj) + ztrad(ipt, k, jg, jj) &
                                   *ztrad(ipt, k, jg, jj)*zrdnd(ipt, k, jg, jj)*zden1
               enddo
               
               zrupbr = zrefb(ipt, 1, jg, jj)        ! direct beam
               zrupdr = zrefd(ipt, 1, jg, jj)        ! diffused
               !  --- ...  up and down-welling fluxes at levels
               zden1 = f_one / (f_one - zrdnd(ipt, 1, jg, jj)*zrupdr)
               zfu = ( ztdbt(ipt, 1, jg, jj)*zrupbr +                                &
               &             (ztdn(ipt, 1, jg, jj) - ztdbt(ipt, 1, jg, jj)) &
                              *zrupdr ) * zden1
               zfd = ztdbt(ipt, 1, jg, jj) + ( ztdn(ipt, 1, jg, jj) &
                              - ztdbt(ipt, 1, jg, jj) +                    &
               &             ztdbt(ipt, 1, jg, jj)*zrupbr &
                              *zrdnd(ipt, 1, jg, jj) ) * zden1
               ! end call swflux
               zfua(ipt, 1, jg, jj) = zfu
               zfda(ipt, 1, jg, jj) = zfd
               zfd0(ipt, jg, jj) = zfd
               !  --- ...  pass from bottom to top
               !$acc loop seq
               do k = 1, nlay
                  kp = k + 1

                  zden1 = f_one / ( f_one - zrupdr*zrefd(ipt, kp, jg, jj) )
                  zrupbr1 = zrefb(ipt, kp, jg, jj) + ( ztrad(ipt, kp, jg, jj) *                         &
                  &                ( (ztrab(ipt, kp, jg, jj) - zldbt(ipt, kp, jg, jj)) &
                                   *zrupdr +              &
                  &                zldbt(ipt, kp, jg, jj)*zrupbr) ) * zden1
                  zrupdr1 = zrefd(ipt, kp, jg, jj) + ztrad(ipt, kp, jg, jj) &
                                   *ztrad(ipt, kp, jg, jj)*zrupdr*zden1
                  
                  zden1 = f_one / (f_one - zrdnd(ipt, kp, jg, jj)*zrupdr1)
                  zfu = ( ztdbt(ipt, kp, jg, jj)*zrupbr1 +                                &
                  &             (ztdn(ipt, kp, jg, jj) - ztdbt(ipt, kp, jg, jj)) &
                                *zrupdr1 ) * zden1
                  zfd = ztdbt(ipt, kp, jg, jj) + ( ztdn(ipt, kp, jg, jj) &
                                - ztdbt(ipt, kp, jg, jj) +                    &
                  &             ztdbt(ipt, kp, jg, jj)*zrupbr1 &
                                *zrdnd(ipt, kp, jg, jj) ) * zden1
                  zfua(ipt, kp, jg, jj) = zfu
                  zfda(ipt, kp, jg, jj) = zfd
                  ! end call swflux
                  !  --- ...  compute upward and downward fluxes at levels
                  zrupbr = zrupbr1
                  zrupdr = zrupdr1
               end do
            end do
         end do
      end do

      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlp1
            !$acc loop vector private(jg2, jb, ib)
            do ipt = 1, nday(jj) ! lab_do_ipt
               !$acc loop seq
               do jg = 1, small_ngptsw ! lab_do_jg
                  jg2 = jg + (j2-1)*small_ngptsw
                  jb = ngb(jg2)
                  ib = jb + 1 - nblow
                  !  --- ...  compute upward and downward fluxes at levels
                  fxup0(ipt, k,ib, jj) = fxup0(ipt, k,ib, jj) + zsolar(ipt, jg, jj)*zfua(ipt, k, jg, jj)
                  fxdn0(ipt, k,ib, jj) = fxdn0(ipt, k,ib, jj) + zsolar(ipt, jg, jj)*zfda(ipt, k, jg, jj)
               end do
            end do
         end do
      end do

               !! --- ...  surface downward beam/diffused flux components
      !$acc parallel loop collapse(2) private(jg2, jb, ib, ibd, zf1, zf2, zb1, zb2) async(async_id)
      do jj = 1, fulljj
         do ipt = 1, ix
            if (ipt .le. nday(jj)) then ! lab_do_ipt
               !$acc loop seq 
               do jg = 1, small_ngptsw ! lab_do_jg
                  jg2 = jg + (j2-1)*small_ngptsw
                  jb = ngb(jg2)
                  ib = jb + 1 - nblow
                  ibd = idxsfc(jb)
                  zb1 = zsolar(ipt, jg, jj)*ztdbt0(ipt, jg, jj)
                  zb2 = zsolar(ipt, jg, jj)*(zfd0(ipt, jg, jj) - ztdbt0(ipt, jg, jj))
                  if (ibd /= 0) then
                     sfbm0(ipt, ibd, jj) = sfbm0(ipt, ibd, jj) + zb1
                     sfdf0(ipt, ibd, jj) = sfdf0(ipt, ibd, jj) + zb2
                  else
                     zf1 = 0.5 * zb1
                     zf2 = 0.5 * zb2
                     sfbm0(ipt, 1, jj) = sfbm0(ipt, 1, jj) + zf1
                     sfdf0(ipt, 1, jj) = sfdf0(ipt, 1, jj) + zf2
                     sfbm0(ipt, 2, jj) = sfbm0(ipt, 2, jj) + zf1
                     sfdf0(ipt, 2, jj) = sfdf0(ipt, 2, jj) + zf2
                  endif
                  !       sfbm0(ibd) = sfbm0(ibd) + zsolar*ztdbt0
                  !       sfdf0(ibd) = sfdf0(ibd) + zsolar*(zfd(1) - ztdbt0)
               end do
            end if
         end do
      end do

               !  --- ...  compute total sky optical parameters, layer reflectance and transmittance
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do jg = 1, small_ngptsw ! lab_do_jg
            !$acc loop vector private(zb1, zb2, zf1, zf2, kp, zc0, zc1, ztau0, &
            !$acc&     zssa0, zasy0, zldbt0, ftind, itind, zssaw, zasyw, za1, za2, &
            !$acc&     ztau1, zssa1, zasy1, zasy3, zgam1, zgam2, zgam3, zgam4, &
            !$acc&     zrefb1, zrefd1, ztrab1, ztrad1, zrk, zrk2, zrp, zrp1, zrm1, &
            !$acc&     zrpp, zrkg1, zrkg3, zrkg4, zr1, zr2, zr3, zr4, zr5, &
            !$acc&     zt1, zt2, zt3, zexm1, zexp1, zexm2, zexp2, ze1r45, zden1, &
            !$acc&     zexp3, zexp4, j1, jg2, jb, ib, ibd, ztdbt0r)
            do ipt = 1, nday(jj) ! lab_do_ipt
               jg2 = jg + (j2-1)*small_ngptsw
               jb = ngb(jg2)
               ib = jb + 1 - nblow
               ibd = idxsfc(jb)
               j1 = idxday(ipt, jj)
               if ( cf1(ipt, jj) > eps ) then

                  !  --- ...  set up toa direct beam and surface values (beam and diff)
                  ztdbt0r = f_one
                  zldbt(ipt, 1, jg, jj) = f_zero
                  !$acc loop seq
                  do k = nlay, 1, -1
                     kp = k + 1
                     zc0 = f_one - cldfrc(ipt, k, jj)
                     zc1 = cldfrc(ipt, k, jj)

                     !  --- ...  saving clear-sky quantities for later total-sky usage
                     ztau0 = max( ftiny, taur(ipt, k,jg2, jj)+taug(ipt, k,jg2, jj)+tauae(ipt, k,ib, jj) )
                     zssa0 = taur(ipt, k,jg2, jj) + tauae(ipt, k,ib, jj)*ssaae(ipt, k,ib, jj)
                     zasy0 = asyae(ipt, k,ib, jj)*ssaae(ipt, k,ib, jj)*tauae(ipt, k,ib, jj)
                     zr1 = ztau0 * sntz(ipt, jj)
                     if ( zr1 <= od_lo ) then
                        zldbt0 = f_one - zr1 + 0.5*zr1*zr1
                     else
                        ftind = zr1 / (bpade + zr1)
                        itind = max(0, min(ntbmx, int(0.5+ntbmx*ftind) ))
                        zldbt0 = exp_tbl(itind)
                     endif
                     
                     if ( zc1 > ftiny ) then          ! it is a cloudy-layer
                        ztau0 = ztau0 + taucw(ipt, k,ib, jj)
                        zssa0 = zssa0 + ssacw(ipt, k,ib, jj)
                        zasy0 = zasy0 + asycw(ipt, k,ib, jj)
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

                        zrefb1 = zrefb(ipt, kp, jg, jj)
                        zrefd1 = zrefd(ipt, kp, jg, jj)
                        ztrab1 = ztrab(ipt, kp, jg, jj)
                        ztrad1 = ztrad(ipt, kp, jg, jj)

                        !  --- ...  compute homogeneous reflectance and transmittance

                        if ( zssaw >= zcrit ) then    ! for conservative scattering
                           za1 = zgam1 * cosz(j1, jj) - zgam3
                           za2 = zgam1 * ztau1

                           !  --- ...  use exponential lookup table for transmittance, or expansion
                           !           of exponential for low optical depth

                           zb1 = min ( ztau1*sntz(ipt, jj) , 500.0 )
                           if ( zb1 <= od_lo ) then
                              zb2 = f_one - zb1 + 0.5*zb1*zb1
                           else
                              ftind = zb1 / (bpade + zb1)
                              itind = ftind*ntbmx + 0.5
                              zb2 = exp_tbl(itind)
                           endif

                           !      ...  collimated beam
                           zrefb(ipt, kp, jg, jj) = max(f_zero, min(f_one,                      &
                           &                      (za2 - za1*(f_one - zb2))/(f_one + za2) ))
                           ztrab(ipt, kp, jg, jj) = max(f_zero, min(f_one, f_one-zrefb(ipt, kp, jg, jj)))

                           !      ...  isotropic incidence
                           zrefd(ipt, kp, jg, jj) = max(f_zero, min(f_one, za2 / (f_one+za2) ))
                           ztrad(ipt, kp, jg, jj) = max(f_zero, min(f_one, f_one - zrefd(ipt, kp, jg, jj) ))

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

                           zb2 = min ( ztau1*sntz(ipt, jj), 500.0 )
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
                              zrefb(ipt, kp, jg, jj) = eps1
                              ztrab(ipt, kp, jg, jj) = zexm2
                           else
                              zden1 = zssa1 / ze1r45
                              zrefb(ipt, kp, jg, jj) = max(f_zero, min(f_one,                    &
                              &                        (zr1*zexp1-zr2*zexm1-zr3*zexm2)*zden1 ))
                              ztrab(ipt, kp, jg, jj) = max(f_zero, min(f_one, zexm2*(f_one -     &
                              &                        (zt1*zexp1-zt2*zexm1-zt3*zexp2)*zden1) ))
                           endif

                           !      ...  diffuse beam
                           zden1 = zr4 / (ze1r45 * zrkg1)
                           zrefd(ipt, kp, jg, jj) = max(f_zero, min(f_one,                      &
                           &                      zgam2*(zexp1 - zexm1)*zden1 ))
                           ztrad(ipt, kp, jg, jj) = max(f_zero, min(f_one, zrk2*zden1 ))
                        endif    ! end if_zssaw_block

                        !  --- ...  combine clear and cloudy contributions for total sky
                        !           and calculate direct beam transmittances

                        zrefb(ipt, kp, jg, jj) = zc0*zrefb1 + zc1*zrefb(ipt, kp, jg, jj)
                        zrefd(ipt, kp, jg, jj) = zc0*zrefd1 + zc1*zrefd(ipt, kp, jg, jj)
                        ztrab(ipt, kp, jg, jj) = zc0*ztrab1 + zc1*ztrab(ipt, kp, jg, jj)
                        ztrad(ipt, kp, jg, jj) = zc0*ztrad1 + zc1*ztrad(ipt, kp, jg, jj)

                        !  --- ...  direct beam transmittance. use exponential lookup table
                        !           for transmittance, or expansion of exponential for low
                        !           optical depth

                        zr1 = ztau1 * sntz(ipt, jj)
                        if ( zr1 <= od_lo ) then
                           zexp3 = f_one - zr1 + 0.5*zr1*zr1
                        else
                           ftind = zr1 / (bpade + zr1)
                           itind = max(0, min(ntbmx, int(0.5+ntbmx*ftind) ))
                           zexp3 = exp_tbl(itind)
                        endif

                        zldbt(ipt, kp, jg, jj) = zc0*zldbt(ipt, kp, jg, jj) + zc1*zexp3
                        ztdbt(ipt, k, jg, jj) = zldbt(ipt, kp, jg, jj) * ztdbt(ipt, kp, jg, jj)

                        !  --- ...  pre-delta-scaling clear and cloudy direct beam transmittance
                        !           (must use 'orig', unscaled cloud optical depth)

                        zr1 = ztau0 * sntz(ipt, jj)
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
                        ztdbt(ipt, k, jg, jj) = zldbt(ipt, kp, jg, jj) * ztdbt(ipt, kp, jg, jj)

                        !  --- ...  pre-delta-scaling clear and cloudy direct beam transmittance
                        ztdbt0r = zldbt0 * ztdbt0r

                     endif    ! end if_zc1_block
                  enddo   ! end do_k_loop
                  ztdbt0(ipt, jg, jj) = ztdbt0r
               end if
            end do
         end do
      end do
                  !  --- ...  perform vertical quadrature
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do jg = 1, small_ngptsw ! lab_do_jg
            !$acc loop vector private(j1, kp, zden1, zfu, zfd, jg2, jb, ib, zrupbr, zrupbr1, zrupdr, zrupdr1)
            do ipt = 1, nday(jj) ! lab_do_ipt
               jg2 = jg + (j2-1)*small_ngptsw
               jb = ngb(jg2)
               ib = jb + 1 - nblow
               if ( cf1(ipt, jj) > eps ) then
                  ! call swflux

                  !  --- ...  upper boundary conditions

                  ztdn (ipt, nlp1, jg, jj) = f_one
                  zrdnd(ipt, nlp1, jg, jj) = f_zero
                  ztdn (ipt, nlay, jg, jj) = ztrab(ipt, nlp1, jg, jj)
                  zrdnd(ipt, nlay, jg, jj) = zrefd(ipt, nlp1, jg, jj)

                  !  --- ...  pass from top to bottom
                  !$acc loop seq
                  do k = nlay, 2, -1
                     zden1 = f_one / (f_one - zrefd(ipt, k, jg, jj)*zrdnd(ipt, k, jg, jj))
                     ztdn (ipt, k-1, jg, jj) = ztdbt(ipt, k, jg, jj)*ztrab(ipt, k, jg, jj) &
                                       + ( ztrad(ipt, k, jg, jj) *                 &
                     &                 ( (ztdn(ipt, k, jg, jj) - ztdbt(ipt, k, jg, jj)) &
                                       + ztdbt(ipt, k, jg, jj) *              &
                     &                 zrefb(ipt, k, jg, jj)*zrdnd(ipt, k, jg, jj) )) * zden1
                     zrdnd(ipt, k-1, jg, jj) = zrefd(ipt, k, jg, jj) + ztrad(ipt, k, jg, jj) &
                                       *ztrad(ipt, k, jg, jj)*zrdnd(ipt, k, jg, jj)*zden1
                  enddo

                  !  --- ...  link lowest layer with surface

                  zrupbr = zrefb(ipt, 1, jg, jj)        ! direct beam
                  zrupdr = zrefd(ipt, 1, jg, jj)        ! diffused
                  
                  zden1 = f_one / (f_one - zrdnd(ipt, 1, jg, jj)*zrupdr)
                  zfu = ( ztdbt(ipt, 1, jg, jj)*zrupbr +                                &
                  &             (ztdn(ipt, 1, jg, jj) - ztdbt(ipt, 1, jg, jj)) &
                                 *zrupdr ) * zden1
                  zfd = ztdbt(ipt, 1, jg, jj) + ( ztdn(ipt, 1, jg, jj) &
                                 - ztdbt(ipt, 1, jg, jj) +                    &
                  &             ztdbt(ipt, 1, jg, jj)*zrupbr &
                                 *zrdnd(ipt, 1, jg, jj) ) * zden1
                                 ! end call swflux
                  zfua(ipt, 1, jg, jj) = zfu
                  zfda(ipt, 1, jg, jj) = zfd

                  zfd0(ipt, jg, jj) = zfd
                  !  --- ...  pass from bottom to top
                  !$acc loop seq
                  do k = 1, nlay
                     kp = k + 1

                     zden1 = f_one / ( f_one - zrupdr*zrefd(ipt, kp, jg, jj) )
                     zrupbr1 = zrefb(ipt, kp, jg, jj) + ( ztrad(ipt, kp, jg, jj) *                         &
                     &                ( (ztrab(ipt, kp, jg, jj) - zldbt(ipt, kp, jg, jj))*zrupdr +              &
                     &                zldbt(ipt, kp, jg, jj)*zrupbr) ) * zden1
                     zrupdr1 = zrefd(ipt, kp, jg, jj) + ztrad(ipt, kp, jg, jj) &
                                      *ztrad(ipt, kp, jg, jj)*zrupdr*zden1

                  !  --- ...  up and down-welling fluxes at levels
                     zden1 = f_one / (f_one - zrdnd(ipt, kp, jg, jj)*zrupdr1)
                     zfu = ( ztdbt(ipt, kp, jg, jj)*zrupbr1 +                                &
                     &             (ztdn(ipt, kp, jg, jj) - ztdbt(ipt, kp, jg, jj)) &
                                    *zrupdr1 ) * zden1
                     zfd = ztdbt(ipt, kp, jg, jj) + ( ztdn(ipt, kp, jg, jj) &
                                    - ztdbt(ipt, kp, jg, jj) +                    &
                     &             ztdbt(ipt, kp, jg, jj)*zrupbr1 &
                                    *zrdnd(ipt, kp, jg, jj) ) * zden1
                                    ! end call swflux
                     zfua(ipt, kp, jg, jj) = zfu
                     zfda(ipt, kp, jg, jj) = zfd
                     zrupbr = zrupbr1
                     zrupdr = zrupdr1
                  end do
               end if
            end do
         end do
      end do
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlp1
            !$acc loop vector private(jg2, jb, ib)
            do ipt = 1, nday(jj) ! lab_do_ipt
               !$acc loop seq
               do jg = 1, small_ngptsw ! lab_do_jg
                  jg2 = jg + (j2-1)*small_ngptsw
                  jb = ngb(jg2)
                  ib = jb + 1 - nblow
                  !  --- ...  compute upward and downward fluxes at levels
                  fxupc(ipt, k,ib, jj) = fxupc(ipt, k,ib, jj) + zsolar(ipt, jg, jj)*zfua(ipt, k, jg, jj)
                  fxdnc(ipt, k,ib, jj) = fxdnc(ipt, k,ib, jj) + zsolar(ipt, jg, jj)*zfda(ipt, k, jg, jj)
               end do
            end do
         end do
      end do
                  !! --- ...  surface downward beam/diffused flux components
      !$acc parallel loop gang collapse(2) private(jg2, jb, ibd, zb1, zb2, zf1, zf2) async(async_id)
      do jj = 1, fulljj
         do ipt = 1, ix
            if (ipt .le. nday(jj)) then ! lab_do_ipt
               !$acc loop seq
               do jg = 1, small_ngptsw ! lab_do_jg
                  jg2 = jg + (j2-1)*small_ngptsw
                  jb = ngb(jg2)
                  ibd = idxsfc(jb)
                  if ( cf1(ipt, jj) > eps ) then
                     zb1 = zsolar(ipt, jg, jj)*ztdbt0(ipt, jg, jj)
                     zb2 = zsolar(ipt, jg, jj)*(zfd0(ipt, jg, jj) - ztdbt0(ipt, jg, jj))

                     if (ibd /= 0) then
                        sfbmc(ipt, ibd, jj) = sfbmc(ipt, ibd, jj) + zb1
                        sfdfc(ipt, ibd, jj) = sfdfc(ipt, ibd, jj) + zb2
                     else
                        zf1 = 0.5 * zb1
                        zf2 = 0.5 * zb2
                        sfbmc(ipt, 1, jj) = sfbmc(ipt, 1, jj) + zf1
                        sfdfc(ipt, 1, jj) = sfdfc(ipt, 1, jj) + zf2
                        sfbmc(ipt, 2, jj) = sfbmc(ipt, 2, jj) + zf1
                        sfdfc(ipt, 2, jj) = sfdfc(ipt, 2, jj) + zf2
                     endif
                     !         sfbmc(ibd) = sfbmc(ibd) + zsolar*ztdbt0
                     !         sfdfc(ibd) = sfdfc(ibd) + zsolar*(zfd(1) - ztdbt0)

                  endif      ! end if_cf1_block

               enddo  ! lab_do_jg
            end if
         end do
      end do
   end do


            !  --- ...  end of g-point loop
      !$acc parallel loop collapse(2) async(async_id)
      do jj = 1, fulljj
         do ipt = 1, ix
            if (ipt .le. nday(jj)) then ! lab_do_ipt
               !$acc loop seq 
               do ib = 1, nbdsw
                  ftoadc(ipt, jj) = ftoadc(ipt, jj) + fxdn0(ipt, nlp1,ib, jj)
                  ftoau0(ipt, jj) = ftoau0(ipt, jj) + fxup0(ipt, nlp1,ib, jj)
                  fsfcu0(ipt, jj) = fsfcu0(ipt, jj) + fxup0(ipt, 1,ib, jj)
                  fsfcd0(ipt, jj) = fsfcd0(ipt, jj) + fxdn0(ipt, 1,ib, jj)
               enddo
            end if
         end do
      end do

      !! --- ...  uv-b surface downward flux
      !$acc parallel loop gang collapse(3) async(async_id)
      do jj = 1, fulljj
         do ib = 1, nbdsw
            do k = 1, nlp1
               !$acc loop vector
               do ipt = 1, nday(jj) ! lab_do_ipt
                  if ( cf1(ipt, jj) <= eps ) then       ! clear column, set total-sky=clear-sky fluxes
                     fxupc(ipt, k,ib, jj) = fxup0(ipt, k,ib, jj)
                     fxdnc(ipt, k,ib, jj) = fxdn0(ipt, k,ib, jj)
                  end if
               enddo
            enddo
         end do
      end do
      !$acc parallel loop collapse(2) private(ibd) async(async_id)
      do jj = 1, fulljj
         do ipt = 1, ix
            if (ipt .le. nday(jj)) then ! lab_do_ipt
               ibd = nuvb - nblow + 1
               suvbf0(ipt, jj) = fxdn0(ipt, 1,ibd, jj)
               if ( cf1(ipt, jj) <= eps ) then       ! clear column, set total-sky=clear-sky fluxes
                  ftoauc(ipt, jj) = ftoau0(ipt, jj)
                  fsfcuc(ipt, jj) = fsfcu0(ipt, jj)
                  fsfcdc(ipt, jj) = fsfcd0(ipt, jj)

                  !! --- ...  surface downward beam/diffused flux components
                  sfbmc(ipt, 1, jj) = sfbm0(ipt, 1, jj)
                  sfdfc(ipt, 1, jj) = sfdf0(ipt, 1, jj)
                  sfbmc(ipt, 2, jj) = sfbm0(ipt, 2, jj)
                  sfdfc(ipt, 2, jj) = sfdf0(ipt, 2, jj)

                  !! --- ...  uv-b surface downward flux
                  suvbfc(ipt, jj) = suvbf0(ipt, jj)
               end if
            end if
         end do
      end do
      !$acc parallel loop gang collapse(3) async(async_id)
      do jj = 1, fulljj
         do ib = 1, nbdsw
            do k = 1, nlp1
               !$acc loop vector 
               do ipt = 1, nday(jj) ! lab_do_ipt
                  if ( cf1(ipt, jj) > eps ) then                        ! cloudy column, compute total-sky fluxes
                     fxupc(ipt, k,ib, jj) = cf1(ipt, jj)*fxupc(ipt, k,ib, jj) + cf0(ipt, jj)*fxup0(ipt, k,ib, jj)
                     fxdnc(ipt, k,ib, jj) = cf1(ipt, jj)*fxdnc(ipt, k,ib, jj) + cf0(ipt, jj)*fxdn0(ipt, k,ib, jj)
                  end if
               enddo
            enddo
         end do
      end do
      !$acc parallel loop collapse(2) async(async_id)
      do jj = 1, fulljj
         do ipt = 1, ix
            if (ipt .le. nday(jj)) then ! lab_do_ipt
               !$acc loop seq
               do ib = 1, nbdsw
                  if ( cf1(ipt, jj) > eps ) then                        ! cloudy column, compute total-sky fluxes
                     ftoauc(ipt, jj) = ftoauc(ipt, jj) + fxupc(ipt, nlp1,ib, jj)
                     fsfcuc(ipt, jj) = fsfcuc(ipt, jj) + fxupc(ipt, 1,ib, jj)
                     fsfcdc(ipt, jj) = fsfcdc(ipt, jj) + fxdnc(ipt, 1,ib, jj)
                  end if
               enddo
            end if
         end do
      end do
      !$acc parallel loop collapse(2) private(ibd) async(async_id)
      do jj = 1, fulljj
         do ipt = 1, ix
            if (ipt .le. nday(jj)) then ! lab_do_ipt
               if ( cf1(ipt, jj) > eps ) then                        ! cloudy column, compute total-sky fluxes
                  !! --- ...  uv-b surface downward flux
                  ibd = nuvb - nblow + 1
                  suvbfc(ipt, jj) = fxdnc(ipt, 1,ibd, jj)

                  !! --- ...  surface downward beam/diffused flux components
                  sfbmc(ipt, 1, jj) = cf1(ipt, jj)*sfbmc(ipt, 1, jj) + cf0(ipt, jj)*sfbm0(ipt, 1, jj)
                  sfbmc(ipt, 2, jj) = cf1(ipt, jj)*sfbmc(ipt, 2, jj) + cf0(ipt, jj)*sfbm0(ipt, 2, jj)
                  sfdfc(ipt, 1, jj) = cf1(ipt, jj)*sfdfc(ipt, 1, jj) + cf0(ipt, jj)*sfdf0(ipt, 1, jj)
                  sfdfc(ipt, 2, jj) = cf1(ipt, jj)*sfdfc(ipt, 2, jj) + cf0(ipt, jj)*sfdf0(ipt, 2, jj)
               endif    ! end if_cf1_block
            end if
         end do
      end do
      !$acc end data

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
     &       forfac,forfrac,indfor,selffac,selffrac,indself, &
             nlay, ix, nday, async_id, fulljj,     &
!  ---  outputs:
     &       sfluxzen, taug, taur                                       &
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
      integer, intent(in) :: nlay, laytrop(ix, fulljj), ix, nday(fulljj), fulljj

      integer, dimension(ix, nlay, fulljj), intent(in) :: indfor, indself,          &
     &       jp, jt, jt1

      real (kind=kind_phys), dimension(ix, nlay, fulljj),  intent(in) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac

      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj),intent(in) :: colamt

!  ---  outputs:
      real (kind=kind_phys), dimension(ix, ngptsw, fulljj), intent(out) :: sfluxzen

      real (kind=kind_phys), dimension(ix, nlay,ngptsw, fulljj), intent(out) ::     &
     &       taug, taur

!  ---  locals:
      !integer, parameter :: small_factor = 10
      !integer, parameter :: small_ix = int(ix/small_factor)
      real (kind=kind_phys) :: fsa, speccomb, specmult, colm1, colm2

      integer, dimension(ix, nlay,nblow:nbhgh, fulljj) :: id0, id1

      integer :: ibd, j, jb, jsa, k, klow, khgh, klim, ks, njb, ns, ipt, jj, i2, i3
      integer :: async_id
!
!===> ... begin here
!
!  --- ...  loop over each spectral band
      !do i2 = small_factor
      !  --- ...  indices for layer optical depth
      !$acc data create(id0, id1) async(async_id)
      !$acc parallel loop gang collapse(3) async(async_id)
      do jj = 1, fulljj
         do jb = nblow, nbhgh
            do k = 1, nlay
               !$acc loop vector
               do ipt = 1, nday(jj) ! lab_do_ipt
                  if (k .le. laytrop(ipt, jj)) then
                     id0(ipt, k,jb, jj) = ((jp(ipt, k, jj)-1)*5 + (jt (ipt, k, jj)-1)) * nspa(jb)
                     id1(ipt, k,jb, jj) = ( jp(ipt, k, jj)   *5 + (jt1(ipt, k, jj)-1)) * nspa(jb)
                  else
                     id0(ipt, k,jb, jj) = ((jp(ipt, k, jj)-13)*5 + (jt (ipt, k, jj)-1)) * nspb(jb)
                     id1(ipt, k,jb, jj) = ((jp(ipt, k, jj)-12)*5 + (jt1(ipt, k, jj)-1)) * nspb(jb)
                  end if
               enddo
            end do
         end do
      end do

      !  --- ...  calculate spectral flux at toa

      !$acc parallel loop collapse(3) private(ks, colm1, colm2, speccomb, &
      !$acc&         specmult, jsa, fsa) async(async_id)
      do jj = 1, fulljj
         do jb = nblow, nbhgh
            do ipt = 1, ix
               if (ipt .le. nday(jj)) then ! lab_do_ipt
                  ibd = ibx(jb)
                  njb = ng (jb)
                  ns  = ngs(jb)
                  if ((jb .eq. 16) .or. (jb .eq. 20) .or. (jb .eq. 23) .or. &
                        (jb .eq. 25) .or. (jb .eq. 26) .or. (jb .eq. 29)) then
                     !$acc loop seq
                     do j = 1, njb
                        sfluxzen(ipt, ngs(jb)+j, jj) = sfluxref01(j,1,ibx(jb))
                     enddo

                  elseif (jb .eq. 27) then
                     !$acc loop seq
                     do j = 1, njb
                        sfluxzen(ipt, ngs(jb)+j, jj) = scalekur * sfluxref01(j,1,ibx(jb))
                     enddo

                  elseif ((jb .eq. 17) .or. (jb .eq. 28)) then
                     ks = nlay
                     !$acc loop seq
                     do k = laytrop(ipt, jj), nlay-1 ! lab_do_k1
                        if (jp(ipt, k, jj)<layreffr(jb) .and. jp(ipt, k+1, jj)>=layreffr(jb)) then
                           ks = k + 1
                           exit ! lab_do_k1
                        endif
                     enddo  ! lab_do_k1

                     colm1 = colamt(ipt, ks,ix1(jb), jj)
                     colm2 = colamt(ipt, ks,ix2(jb), jj)
                     speccomb = colm1 + strrat(jb)*colm2
                     specmult = specwt(jb) * min( oneminus, colm1/speccomb )
                     jsa = 1 + int( specmult )
                     fsa = mod(specmult, f_one)
                     !$acc loop seq
                     do j = 1, njb
                        sfluxzen(ipt, ngs(jb)+j, jj) = sfluxref02(j,jsa,ibd)                   &
                        &           + fsa * (sfluxref02(j,jsa+1,ibd) - sfluxref02(j,jsa,ibd))
                     enddo

                  else
                     ks = laytrop(ipt, jj)
                     !$acc loop seq
                     do k = 1, laytrop(ipt, jj)-1 ! lab_do_k2
                        if (jp(ipt, k, jj)<layreffr(jb) .and. jp(ipt, k+1, jj)>=layreffr(jb)) then
                           ks = k + 1
                           exit ! lab_do_k2
                        endif
                     enddo  ! lab_do_k2

                     colm1 = colamt(ipt, ks,ix1(jb), jj)
                     colm2 = colamt(ipt, ks,ix2(jb), jj)
                     speccomb = colm1 + strrat(jb)*colm2
                     specmult = specwt(jb) * min( oneminus, colm1/speccomb )
                     jsa = 1 + int( specmult )
                     fsa = mod(specmult, f_one)
                     !$acc loop seq
                     do j = 1, njb
                        sfluxzen(ipt, ngs(jb)+j, jj) = sfluxref03(j,jsa,ibd)                   &
                        &           + fsa * (sfluxref03(j,jsa+1,ibd) - sfluxref03(j,jsa,ibd))
                     enddo
                  end if
               end if
            end do
         end do
      enddo


            !  --- ...  call taumol## to calculate layer optical depth
      call taumol16( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
     &       sfluxzen, taug, taur, async_id, fulljj, id0, id1)

      call taumol17( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
     &       sfluxzen, taug, taur, async_id, fulljj, id0, id1)

      call taumol18( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
     &       sfluxzen, taug, taur, async_id, fulljj, id0, id1)

      call taumol19( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
     &       sfluxzen, taug, taur, async_id, fulljj, id0, id1)

      call taumol20( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
     &       sfluxzen, taug, taur, async_id, fulljj, id0, id1)

      call taumol21( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
     &       sfluxzen, taug, taur, async_id, fulljj, id0, id1)

      call taumol22( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
     &       sfluxzen, taug, taur, async_id, fulljj, id0, id1)

      call taumol23( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
     &       sfluxzen, taug, taur, async_id, fulljj, id0, id1)

      call taumol24( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
     &       sfluxzen, taug, taur, async_id, fulljj, id0, id1)

      call taumol25( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
     &       sfluxzen, taug, taur, async_id, fulljj, id0, id1)

      call taumol26( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
     &       sfluxzen, taug, taur, async_id, fulljj, id0, id1)

      call taumol27( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
     &       sfluxzen, taug, taur, async_id, fulljj, id0, id1)

      call taumol28( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
     &       sfluxzen, taug, taur, async_id, fulljj, id0, id1)

      call taumol29( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
     &       sfluxzen, taug, taur, async_id, fulljj, id0, id1)
      !$acc end data


!...................................
      end subroutine taumol
!-----------------------------------

!-----------------------------------
      subroutine taumol16                                                 &
!...................................
     &     ( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
     &       sfluxzen, taug, taur, async_id, fulljj, id0, id1                                        &
     &     )

!  ------------------------------------------------------------------  !
!     band 16:  2600-3250 cm-1 (low - h2o,ch4; high - ch4)             !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb16

      integer :: nlay, laytrop(ix, fulljj), ix, nday(fulljj)
      integer, dimension(ix, nlay, fulljj) :: indfor, indself,          &
     &       jp, jt, jt1
      real (kind=kind_phys), dimension(ix, nlay, fulljj) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac
      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj) :: colamt
      real (kind=kind_phys), dimension(ix, ngptsw, fulljj) :: sfluxzen
      real (kind=kind_phys), dimension(ix, nlay,ngptsw, fulljj) ::     &
     &       taug, taur
      integer, dimension(ix, nlay,nblow:nbhgh, fulljj) :: id0, id1
!  ---  locals:
      integer :: ipt, jj, async_id, fulljj

      real (kind=kind_phys) :: speccomb, specmult, tauray, fs, fs1,     &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111

      integer :: ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14
      integer :: inds, indf, indsp, indfp, j, js, k

!
!===> ... begin here
!

!  --- ... compute the optical depth by interpolating in ln(pressure),
!          temperature, and appropriate species.  below laytrop, the water
!          vapor self-continuum is interpolated (in temperature) separately.
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !!$acc cache(absa, absb, selfref, forref)
            !$acc loop vector private(speccomb, specmult, js, fs, fs1, fac000, &
            !$acc&     fac010, fac100, fac110, fac001, fac011, fac101, fac111, &
            !$acc&     ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14, &
            !$acc&     inds, indf, indsp, indfp, tauray)
            do ipt = 1, nday(jj) ! lab_do_ipt
               tauray = colmol(ipt, k, jj) * rayl
               if (k .le. laytrop(ipt, jj)) then
                  speccomb = colamt(ipt, k,1, jj) + strrat(16)*colamt(ipt, k,5, jj)
                  specmult = 8.0 * min( oneminus, colamt(ipt, k,1, jj)/speccomb )

                  js = 1 + int( specmult )
                  fs = mod( specmult, f_one )
                  fs1= f_one - fs
                  fac000 = fs1 * fac00(ipt, k, jj)
                  fac010 = fs1 * fac10(ipt, k, jj)
                  fac100 = fs  * fac00(ipt, k, jj)
                  fac110 = fs  * fac10(ipt, k, jj)
                  fac001 = fs1 * fac01(ipt, k, jj)
                  fac011 = fs1 * fac11(ipt, k, jj)
                  fac101 = fs  * fac01(ipt, k, jj)
                  fac111 = fs  * fac11(ipt, k, jj)

                  ind01 = id0(ipt, k,16, jj) + js
                  ind02 = ind01 + 1
                  ind03 = ind01 + 9
                  ind04 = ind01 + 10
                  ind11 = id1(ipt, k,16, jj) + js
                  ind12 = ind11 + 1
                  ind13 = ind11 + 9
                  ind14 = ind11 + 10
                  inds = indself(ipt, k, jj)
                  indf = indfor (ipt, k, jj)
                  indsp= inds + 1
                  indfp= indf + 1
                  
                  !$acc loop seq
                  do j = 1, ng16
                     taug(ipt, k,ns16+j, jj) = speccomb                                     &
                     &        *( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)        &
                     &        +  fac010 * absa(ind03,j) + fac110 * absa(ind04,j)        &
                     &        +  fac001 * absa(ind11,j) + fac101 * absa(ind12,j)        &
                     &        +  fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )      &
                     &        + colamt(ipt, k,1, jj) * (selffac(ipt, k, jj) * (selfref(inds,j)            &
                     &        + selffrac(ipt, k, jj) * (selfref(indsp,j)-selfref(inds,j)))       &
                     &        + forfac(ipt, k, jj) * (forref(indf,j) + forfrac(ipt, k, jj)                &
                     &        * (forref(indfp,j) - forref(indf,j))))
                     taur(ipt, k,ns16+j, jj) = tauray
                  enddo
               else
                  ind01 = id0(ipt, k,16, jj) + 1
                  ind02 = ind01 + 1
                  ind11 = id1(ipt, k,16, jj) + 1
                  ind12 = ind11 + 1
                  
                  !$acc loop seq
                  do j = 1, ng16
                     taug(ipt, k,ns16+j, jj) = colamt(ipt, k,5, jj)                                  &
                     &      * ( fac00(ipt, k, jj)*absb(ind01,j) + fac10(ipt, k, jj)*absb(ind02,j)         &
                     &      +   fac01(ipt, k, jj)*absb(ind11,j) + fac11(ipt, k, jj)*absb(ind12,j) )
                     taur(ipt, k,ns16+j, jj) = tauray
                  enddo
               end if
            enddo
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
     &       sfluxzen, taug, taur, async_id, fulljj, id0, id1                                        &
     &     )

!  ------------------------------------------------------------------  !
!     band 17:  3250-4000 cm-1 (low - h2o,co2; high - h2o,co2)         !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb17

      integer :: nlay, laytrop(ix, fulljj), ix, nday(fulljj)
      integer, dimension(ix, nlay, fulljj) :: indfor, indself,          &
     &       jp, jt, jt1
      real (kind=kind_phys), dimension(ix, nlay, fulljj) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac
      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj) :: colamt
      real (kind=kind_phys), dimension(ix, ngptsw, fulljj) :: sfluxzen
      real (kind=kind_phys), dimension(ix, nlay,ngptsw, fulljj) ::     &
     &       taug, taur
      integer, dimension(ix, nlay,nblow:nbhgh, fulljj) :: id0, id1
!  ---  locals:
      integer :: ipt, jj, async_id, fulljj
      real (kind=kind_phys) :: speccomb, specmult, tauray, fs, fs1,     &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111

      integer :: ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14
      integer :: inds, indf, indsp, indfp, j, js, k

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.
      
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !!$acc cache(selfref, forref)
            !$acc loop vector private(speccomb, specmult, js, fs, fs1, fac000, &
            !$acc&     fac010, fac100, fac110, fac001, fac011, fac101, fac111, &
            !$acc&     ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14, &
            !$acc&     inds, indf, indsp, indfp, tauray)
            do ipt = 1, nday(jj) ! lab_do_ipt
               tauray = colmol(ipt, k, jj) * rayl
               if (k .le. laytrop(ipt, jj)) then
                  speccomb = colamt(ipt, k,1, jj) + strrat(17)*colamt(ipt, k,2, jj)
                  specmult = 8.0 * min(oneminus, colamt(ipt, k,1, jj) / speccomb)

                  js = 1 + int(specmult)
                  fs = mod(specmult, f_one)
                  fs1= f_one - fs
                  fac000 = fs1 * fac00(ipt, k, jj)
                  fac010 = fs1 * fac10(ipt, k, jj)
                  fac100 = fs  * fac00(ipt, k, jj)
                  fac110 = fs  * fac10(ipt, k, jj)
                  fac001 = fs1 * fac01(ipt, k, jj)
                  fac011 = fs1 * fac11(ipt, k, jj)
                  fac101 = fs  * fac01(ipt, k, jj)
                  fac111 = fs  * fac11(ipt, k, jj)

                  ind01 = id0(ipt, k,17, jj) + js
                  ind02 = ind01 + 1
                  ind03 = ind01 + 9
                  ind04 = ind01 + 10
                  ind11 = id1(ipt, k,17, jj) + js
                  ind12 = ind11 + 1
                  ind13 = ind11 + 9
                  ind14 = ind11 + 10

                  inds = indself(ipt, k, jj)
                  indf = indfor (ipt, k, jj)
                  indsp= inds + 1
                  indfp= indf + 1
                  
                  !$acc loop seq
                  do j = 1, ng17
                     taug(ipt, k,ns17+j, jj) = speccomb                                     &
                     &        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
                     &        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
                     &        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
                     &        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )     &
                     &        + colamt(ipt, k,1, jj) * (selffac(ipt, k, jj) * (selfref(inds,j)            &
                     &        + selffrac(ipt, k, jj) * (selfref(indsp,j)-selfref(inds,j)))       &
                     &        + forfac(ipt, k, jj) * (forref(indf,j) + forfrac(ipt, k, jj)                &
                     &        * (forref(indfp,j) - forref(indf,j))))
                     taur(ipt, k,ns17+j, jj) = tauray
                  enddo
               else
                  speccomb = colamt(ipt, k,1, jj) + strrat(17)*colamt(ipt, k,2, jj)
                  specmult = 4.0 * min(oneminus, colamt(ipt, k,1, jj) / speccomb)

                  js = 1 + int(specmult)
                  fs = mod(specmult, f_one)
                  fs1= f_one - fs
                  fac000 = fs1 * fac00(ipt, k, jj)
                  fac010 = fs1 * fac10(ipt, k, jj)
                  fac100 = fs  * fac00(ipt, k, jj)
                  fac110 = fs  * fac10(ipt, k, jj)
                  fac001 = fs1 * fac01(ipt, k, jj)
                  fac011 = fs1 * fac11(ipt, k, jj)
                  fac101 = fs  * fac01(ipt, k, jj)
                  fac111 = fs  * fac11(ipt, k, jj)

                  ind01 = id0(ipt, k,17, jj) + js
                  ind02 = ind01 + 1
                  ind03 = ind01 + 5
                  ind04 = ind01 + 6
                  ind11 = id1(ipt, k,17, jj) + js
                  ind12 = ind11 + 1
                  ind13 = ind11 + 5
                  ind14 = ind11 + 6

                  indf = indfor(ipt, k, jj)
                  indfp= indf + 1
                  
                  !$acc loop seq
                  do j = 1, ng17
                     taug(ipt, k,ns17+j, jj) = speccomb                                     &
                     &        * ( fac000 * absb(ind01,j) + fac100 * absb(ind02,j)       &
                     &        +   fac010 * absb(ind03,j) + fac110 * absb(ind04,j)       &
                     &        +   fac001 * absb(ind11,j) + fac101 * absb(ind12,j)       &
                     &        +   fac011 * absb(ind13,j) + fac111 * absb(ind14,j) )     &
                     &        + colamt(ipt, k,1, jj) * forfac(ipt, k, jj) * (forref(indf,j)               &
                     &        + forfrac(ipt, k, jj) * (forref(indfp,j) - forref(indf,j)))
                     taur(ipt, k,ns17+j, jj) = tauray
                  enddo
               end if
            enddo
         end do
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
     &       sfluxzen, taug, taur, async_id, fulljj, id0, id1                                        &
     &     )
     
!  ------------------------------------------------------------------  !
!     band 18:  4000-4650 cm-1 (low - h2o,ch4; high - ch4)             !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb18

      integer :: nlay, laytrop(ix, fulljj), ix, nday(fulljj)
      integer, dimension(ix, nlay, fulljj) :: indfor, indself,          &
     &       jp, jt, jt1
      real (kind=kind_phys), dimension(ix, nlay, fulljj) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac
      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj) :: colamt
      real (kind=kind_phys), dimension(ix, ngptsw, fulljj) :: sfluxzen
      real (kind=kind_phys), dimension(ix, nlay,ngptsw, fulljj) ::     &
     &       taug, taur
      integer, dimension(ix, nlay,nblow:nbhgh, fulljj) :: id0, id1
!  ---  locals:
      integer :: ipt, jj, async_id, fulljj
      real (kind=kind_phys) :: speccomb, specmult, tauray, fs, fs1,     &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111

      integer :: ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14
      integer :: inds, indf, indsp, indfp, j, js, k

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !!$acc cache(absa, absb, selfref, forref)
            !$acc loop vector private(speccomb, specmult, js, fs, fs1, fac000, &
            !$acc&     fac010, fac100, fac110, fac001, fac011, fac101, fac111, &
            !$acc&     ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14, &
            !$acc&     inds, indf, indsp, indfp, tauray)
            do ipt = 1, nday(jj) ! lab_do_ipt
               tauray = colmol(ipt, k, jj) * rayl
               if (k .le. laytrop(ipt, jj)) then
                  speccomb = colamt(ipt, k,1, jj) + strrat(18)*colamt(ipt, k,5, jj)
                  specmult = 8.0 * min(oneminus, colamt(ipt, k,1, jj) / speccomb)

                  js = 1 + int(specmult)
                  fs = mod(specmult, f_one)
                  fs1= f_one - fs
                  fac000 = fs1 * fac00(ipt, k, jj)
                  fac010 = fs1 * fac10(ipt, k, jj)
                  fac100 = fs  * fac00(ipt, k, jj)
                  fac110 = fs  * fac10(ipt, k, jj)
                  fac001 = fs1 * fac01(ipt, k, jj)
                  fac011 = fs1 * fac11(ipt, k, jj)
                  fac101 = fs  * fac01(ipt, k, jj)
                  fac111 = fs  * fac11(ipt, k, jj)

                  ind01 = id0(ipt, k,18, jj) + js
                  ind02 = ind01 + 1
                  ind03 = ind01 + 9
                  ind04 = ind01 + 10
                  ind11 = id1(ipt, k,18, jj) + js
                  ind12 = ind11 + 1
                  ind13 = ind11 + 9
                  ind14 = ind11 + 10

                  inds = indself(ipt, k, jj)
                  indf = indfor (ipt, k, jj)
                  indsp= inds + 1
                  indfp= indf + 1
                  
                  !$acc loop seq
                  do j = 1, ng18
                     taug(ipt, k,ns18+j, jj) = speccomb                                     &
                     &        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
                     &        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
                     &        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
                     &        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )     &
                     &        + colamt(ipt, k,1, jj) * (selffac(ipt, k, jj) * (selfref(inds,j)            &
                     &        + selffrac(ipt, k, jj) * (selfref(indsp,j)-selfref(inds,j)))       &
                     &        + forfac(ipt, k, jj) * (forref(indf,j) + forfrac(ipt, k, jj)                &
                     &        * (forref(indfp,j) - forref(indf,j))))
                     taur(ipt, k,ns18+j, jj) = tauray
                  enddo
               else
                  ind01 = id0(ipt, k,18, jj) + 1
                  ind02 = ind01 + 1
                  ind11 = id1(ipt, k,18, jj) + 1
                  ind12 = ind11 + 1
                  
                  !$acc loop seq
                  do j = 1, ng18
                     taug(ipt, k,ns18+j, jj) = colamt(ipt, k,5, jj)                                  &
                     &        * ( fac00(ipt, k, jj)*absb(ind01,j) + fac10(ipt, k, jj)*absb(ind02,j)       &
                     &        +   fac01(ipt, k, jj)*absb(ind11,j) + fac11(ipt, k, jj)*absb(ind12,j) )
                     taur(ipt, k,ns18+j, jj) = tauray
                  enddo
               end if
            enddo
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
     &       sfluxzen, taug, taur, async_id, fulljj, id0, id1                                        &
     &     )
     
!  ------------------------------------------------------------------  !
!     band 19:  4650-5150 cm-1 (low - h2o,co2; high - co2)             !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb19

      integer :: nlay, laytrop(ix, fulljj), ix, nday(fulljj)
      integer, dimension(ix, nlay, fulljj) :: indfor, indself,          &
     &       jp, jt, jt1
      real (kind=kind_phys), dimension(ix, nlay, fulljj) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac
      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj) :: colamt
      real (kind=kind_phys), dimension(ix, ngptsw, fulljj) :: sfluxzen
      real (kind=kind_phys), dimension(ix, nlay,ngptsw, fulljj) ::     &
     &       taug, taur
      integer, dimension(ix, nlay,nblow:nbhgh, fulljj) :: id0, id1
!  ---  locals:
      integer :: ipt, jj, async_id, fulljj
      real (kind=kind_phys) :: speccomb, specmult, tauray, fs, fs1,     &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111

      integer :: ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14
      integer :: inds, indf, indsp, indfp, j, js, k

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.
      
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !!$acc cache(absa, absb, selfref, forref)
            !$acc loop vector private(speccomb, specmult, js, fs, fs1, fac000, &
            !$acc&     fac010, fac100, fac110, fac001, fac011, fac101, fac111, &
            !$acc&     ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14, &
            !$acc&     indsp, indfp, tauray)
            do ipt = 1, nday(jj) ! lab_do_ipt
               tauray = colmol(ipt, k, jj) * rayl
               if (k .le. laytrop(ipt, jj)) then
                  speccomb = colamt(ipt, k,1, jj) + strrat(19)*colamt(ipt, k,2, jj)
                  specmult = 8.0 * min(oneminus, colamt(ipt, k,1, jj) / speccomb)

                  js = 1 + int(specmult)
                  fs = mod(specmult, f_one)
                  fs1= f_one - fs
                  fac000 = fs1 * fac00(ipt, k, jj)
                  fac010 = fs1 * fac10(ipt, k, jj)
                  fac100 = fs  * fac00(ipt, k, jj)
                  fac110 = fs  * fac10(ipt, k, jj)
                  fac001 = fs1 * fac01(ipt, k, jj)
                  fac011 = fs1 * fac11(ipt, k, jj)
                  fac101 = fs  * fac01(ipt, k, jj)
                  fac111 = fs  * fac11(ipt, k, jj)

                  ind01 = id0(ipt, k,19, jj) + js
                  ind02 = ind01 + 1
                  ind03 = ind01 + 9
                  ind04 = ind01 + 10
                  ind11 = id1(ipt, k,19, jj) + js
                  ind12 = ind11 + 1
                  ind13 = ind11 + 9
                  ind14 = ind11 + 10

                  inds = indself(ipt, k, jj)
                  indf = indfor (ipt, k, jj)
                  indsp= inds + 1
                  indfp= indf + 1
                  
                  !$acc loop seq
                  do j = 1, ng19
                     taug(ipt, k,ns19+j, jj) = speccomb                                     &
                     &        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
                     &        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
                     &        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
                     &        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )     &
                     &        + colamt(ipt, k,1, jj) * (selffac(ipt, k, jj) * (selfref(inds,j)            &
                     &        + selffrac(ipt, k, jj) * (selfref(indsp,j)-selfref(inds,j)))       &
                     &        + forfac(ipt, k, jj) * (forref(indf,j) + forfrac(ipt, k, jj)                &
                     &        * (forref(indfp,j) - forref(indf,j))))
                     taur(ipt, k,ns19+j, jj) = tauray
                  enddo
               else
                  ind01 = id0(ipt, k,19, jj) + 1
                  ind02 = ind01 + 1
                  ind11 = id1(ipt, k,19, jj) + 1
                  ind12 = ind11 + 1
                  
                  !$acc loop seq
                  do j = 1, ng19
                     taug(ipt, k,ns19+j, jj) = colamt(ipt, k,2, jj)                                  &
                     &        * ( fac00(ipt, k, jj)*absb(ind01,j) + fac10(ipt, k, jj)*absb(ind02,j)       &
                     &        +   fac01(ipt, k, jj)*absb(ind11,j) + fac11(ipt, k, jj)*absb(ind12,j) ) 
                     taur(ipt, k,ns19+j, jj) = tauray
                  enddo
               end if
            enddo
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
     &       sfluxzen, taug, taur, async_id, fulljj, id0, id1                                        &
     &     )
     
!  ------------------------------------------------------------------  !
!     band 20:  5150-6150 cm-1 (low - h2o; high - h2o)                 !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb20

      integer :: nlay, laytrop(ix, fulljj), ix, nday(fulljj)
      integer, dimension(ix, nlay, fulljj) :: indfor, indself,          &
     &       jp, jt, jt1
      real (kind=kind_phys), dimension(ix, nlay, fulljj) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac
      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj) :: colamt
      real (kind=kind_phys), dimension(ix, ngptsw, fulljj) :: sfluxzen
      real (kind=kind_phys), dimension(ix, nlay,ngptsw, fulljj) ::     &
     &       taug, taur
      integer, dimension(ix, nlay,nblow:nbhgh, fulljj) :: id0, id1
!  ---  locals:
      integer :: ipt, jj, async_id, fulljj
      real (kind=kind_phys) :: tauray

      integer :: ind01, ind02, ind11, ind12
      integer :: inds, indf, indsp, indfp, j, k

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.

      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !!$acc cache(absa, selfref, forref, absb, absch4)
            !$acc loop vector private(ind01, ind02, ind11, ind12, inds, indf, &
            !$acc&     indsp, indfp, tauray)
            do ipt = 1, nday(jj)
               tauray = colmol(ipt, k, jj) * rayl
               if (k .le. laytrop(ipt, jj)) then
                  ind01 = id0(ipt, k,20, jj) + 1
                  ind02 = ind01 + 1
                  ind11 = id1(ipt, k,20, jj) + 1
                  ind12 = ind11 + 1

                  inds = indself(ipt, k, jj)
                  indf = indfor (ipt, k, jj)
                  indsp= inds + 1
                  indfp= indf + 1
                  
                  !$acc loop seq
                  do j = 1, ng20
                     taug(ipt, k,ns20+j, jj) = colamt(ipt, k,1, jj)                                  &
                     &        * ( (fac00(ipt, k, jj)*absa(ind01,j) + fac10(ipt, k, jj)*absa(ind02,j)      &
                     &        +    fac01(ipt, k, jj)*absa(ind11,j) + fac11(ipt, k, jj)*absa(ind12,j))     &
                     &        +   selffac(ipt, k, jj) * (selfref(inds,j) + selffrac(ipt, k, jj)           &
                     &        *   (selfref(indsp,j) - selfref(inds,j)))                 &
                     &        +   forfac(ipt, k, jj) * (forref(indf,j) + forfrac(ipt, k, jj)              &
                     &        *   (forref(indfp,j) - forref(indf,j))) )                 &
                     &        + colamt(ipt, k,5, jj) * absch4(j)
                     taur(ipt, k,ns20+j, jj) = tauray
                  enddo
               else
                  ind01 = id0(ipt, k,20, jj) + 1
                  ind02 = ind01 + 1
                  ind11 = id1(ipt, k,20, jj) + 1
                  ind12 = ind11 + 1

                  indf = indfor(ipt, k, jj)
                  indfp= indf + 1

                  !$acc loop seq
                  do j = 1, ng20
                     taug(ipt, k,ns20+j, jj) = colamt(ipt, k,1, jj)                                  &
                     &        * ( fac00(ipt, k, jj)*absb(ind01,j) + fac10(ipt, k, jj)*absb(ind02,j)       &
                     &        +   fac01(ipt, k, jj)*absb(ind11,j) + fac11(ipt, k, jj)*absb(ind12,j)       &
                     &        +   forfac(ipt, k, jj) * (forref(indf,j) + forfrac(ipt, k, jj)              &
                     &        *   (forref(indfp,j) - forref(indf,j))) )                 &
                     &        + colamt(ipt, k,5, jj) * absch4(j)
                     taur(ipt, k,ns20+j, jj) = tauray
                  enddo
               end if
            enddo
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
     &       sfluxzen, taug, taur, async_id, fulljj, id0, id1                                        &
     &     )
     
!  ------------------------------------------------------------------  !
!     band 21:  6150-7700 cm-1 (low - h2o,co2; high - h2o,co2)         !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb21

      integer :: nlay, laytrop(ix, fulljj), ix, nday(fulljj)
      integer, dimension(ix, nlay, fulljj) :: indfor, indself,          &
     &       jp, jt, jt1
      real (kind=kind_phys), dimension(ix, nlay, fulljj) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac
      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj) :: colamt
      real (kind=kind_phys), dimension(ix, ngptsw, fulljj) :: sfluxzen
      real (kind=kind_phys), dimension(ix, nlay,ngptsw, fulljj) ::     &
     &       taug, taur
      integer, dimension(ix, nlay,nblow:nbhgh, fulljj) :: id0, id1
!  ---  locals:
      integer :: ipt, jj, async_id, fulljj
      real (kind=kind_phys) :: speccomb, specmult, tauray, fs, fs1,     &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111

      integer :: ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14
      integer :: inds, indf, indsp, indfp, j, js, k

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.
      
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !!$acc cache(absa, selfref, forref)
            !$acc loop vector private(speccomb, specmult, js, fs, fs1, fac000, &
            !$acc&     fac010, fac100, fac110, fac001, fac011, fac101, fac111, &
            !$acc&     ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14, &
            !$acc&     inds, indf, indsp, indfp, tauray)
            do ipt = 1, nday(jj) ! lab_do_ipt
               tauray = colmol(ipt, k, jj) * rayl
               if (k .le. laytrop(ipt, jj)) then
                  speccomb = colamt(ipt, k,1, jj) + strrat(21)*colamt(ipt, k,2, jj)
                  specmult = 8.0 * min(oneminus, colamt(ipt, k,1, jj) / speccomb)

                  js = 1 + int(specmult)
                  fs = mod(specmult, f_one)
                  fs1= f_one - fs
                  fac000 = fs1 * fac00(ipt, k, jj)
                  fac010 = fs1 * fac10(ipt, k, jj)
                  fac100 = fs  * fac00(ipt, k, jj)
                  fac110 = fs  * fac10(ipt, k, jj)
                  fac001 = fs1 * fac01(ipt, k, jj)
                  fac011 = fs1 * fac11(ipt, k, jj)
                  fac101 = fs  * fac01(ipt, k, jj)
                  fac111 = fs  * fac11(ipt, k, jj)

                  ind01 = id0(ipt, k,21, jj) + js
                  ind02 = ind01 + 1
                  ind03 = ind01 + 9
                  ind04 = ind01 + 10
                  ind11 = id1(ipt, k,21, jj) + js
                  ind12 = ind11 + 1
                  ind13 = ind11 + 9
                  ind14 = ind11 + 10

                  inds = indself(ipt, k, jj)
                  indf = indfor (ipt, k, jj)
                  indsp= inds + 1
                  indfp= indf + 1
                  
                  !$acc loop seq
                  do j = 1, ng21
                     taug(ipt, k,ns21+j, jj) = speccomb                                     &
                     &        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
                     &        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
                     &        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
                     &        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )     &
                     &        + colamt(ipt, k,1, jj) * (selffac(ipt, k, jj) * (selfref(inds,j)            &
                     &        + selffrac(ipt, k, jj) * (selfref(indsp,j) - selfref(inds,j)))     &
                     &        + forfac(ipt, k, jj) * (forref(indf,j) + forfrac(ipt, k, jj)                &
                     &        * (forref(indfp,j) - forref(indf,j))))
                     taur(ipt, k,ns21+j, jj) = tauray
                  enddo
               else
                  speccomb = colamt(ipt, k,1, jj) + strrat(21)*colamt(ipt, k,2, jj)
                  specmult = 4.0 * min(oneminus, colamt(ipt, k,1, jj) / speccomb)

                  js = 1 + int(specmult)
                  fs = mod(specmult, f_one)
                  fs1= f_one - fs
                  fac000 = fs1 * fac00(ipt, k, jj)
                  fac010 = fs1 * fac10(ipt, k, jj)
                  fac100 = fs  * fac00(ipt, k, jj)
                  fac110 = fs  * fac10(ipt, k, jj)
                  fac001 = fs1 * fac01(ipt, k, jj)
                  fac011 = fs1 * fac11(ipt, k, jj)
                  fac101 = fs  * fac01(ipt, k, jj)
                  fac111 = fs  * fac11(ipt, k, jj)

                  ind01 = id0(ipt, k,21, jj) + js
                  ind02 = ind01 + 1
                  ind03 = ind01 + 5
                  ind04 = ind01 + 6
                  ind11 = id1(ipt, k,21, jj) + js
                  ind12 = ind11 + 1
                  ind13 = ind11 + 5
                  ind14 = ind11 + 6

                  indf = indfor(ipt, k, jj)
                  indfp= indf + 1
                  
                  !$acc loop seq
                  do j = 1, ng21
                     taug(ipt, k,ns21+j, jj) = speccomb                                     &
                     &        * ( fac000 * absb(ind01,j) + fac100 * absb(ind02,j)       &
                     &        +   fac010 * absb(ind03,j) + fac110 * absb(ind04,j)       &
                     &        +   fac001 * absb(ind11,j) + fac101 * absb(ind12,j)       &
                     &        +   fac011 * absb(ind13,j) + fac111 * absb(ind14,j) )     &
                     &        + colamt(ipt, k,1, jj) * forfac(ipt, k, jj) * (forref(indf,j)               &
                     &        + forfrac(ipt, k, jj) * (forref(indfp,j) - forref(indf,j)))
                     taur(ipt, k,ns21+j, jj) = tauray
                  enddo
               end if
            enddo
         end do
      end do

!...................................
      end subroutine taumol21
!-----------------------------------


!-----------------------------------
      subroutine taumol22                                                 &
!...................................
     &     ( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay, ix, nday,     &
     &       sfluxzen, taug, taur, async_id, fulljj, id0, id1                                        &
     &     )
     
!  ------------------------------------------------------------------  !
!     band 22:  7700-8050 cm-1 (low - h2o,o2; high - o2)               !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb22

      integer :: nlay, laytrop(ix, fulljj), ix, nday(fulljj)
      integer, dimension(ix, nlay, fulljj) :: indfor, indself,          &
     &       jp, jt, jt1
      real (kind=kind_phys), dimension(ix, nlay, fulljj) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac
      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj) :: colamt
      real (kind=kind_phys), dimension(ix, ngptsw, fulljj) :: sfluxzen
      real (kind=kind_phys), dimension(ix, nlay,ngptsw, fulljj) ::     &
     &       taug, taur
      integer, dimension(ix, nlay,nblow:nbhgh, fulljj) :: id0, id1
!  ---  locals:
      integer :: ipt, jj, async_id, fulljj
      real (kind=kind_phys) :: speccomb, specmult, tauray, fs, fs1,     &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111,  &
     &       o2adj, o2cont, o2tem

      integer :: ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14
      integer :: inds, indf, indsp, indfp, j, js, k

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
      
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !!$acc cache(absa, selfref, forref, absb)
            !$acc loop vector private(speccomb, specmult, js, fs, fs1, fac000, &
            !$acc&     fac010, fac100, fac110, fac001, fac011, fac101, fac111, &
            !$acc&     ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14, &
            !$acc&     inds, indf, indsp, indfp, o2cont, tauray)
            do ipt = 1, nday(jj) ! lab_do_ipt
               tauray = colmol(ipt, k, jj) * rayl
               if (k .le. laytrop(ipt, jj)) then
                  o2cont   = o2tem * colamt(ipt, k,6, jj)
                  speccomb = colamt(ipt, k,1, jj) + strrat(22)*colamt(ipt, k,6, jj)
                  specmult = 8.0 * min(oneminus, colamt(ipt, k,1, jj) / speccomb)

                  js = 1 + int(specmult)
                  fs = mod(specmult, f_one)
                  fs1= f_one - fs
                  fac000 = fs1 * fac00(ipt, k, jj)
                  fac010 = fs1 * fac10(ipt, k, jj)
                  fac100 = fs  * fac00(ipt, k, jj)
                  fac110 = fs  * fac10(ipt, k, jj)
                  fac001 = fs1 * fac01(ipt, k, jj)
                  fac011 = fs1 * fac11(ipt, k, jj)
                  fac101 = fs  * fac01(ipt, k, jj)
                  fac111 = fs  * fac11(ipt, k, jj)

                  ind01 = id0(ipt, k,22, jj) + js
                  ind02 = ind01 + 1
                  ind03 = ind01 + 9
                  ind04 = ind01 + 10
                  ind11 = id1(ipt, k,22, jj) + js
                  ind12 = ind11 + 1
                  ind13 = ind11 + 9
                  ind14 = ind11 + 10

                  inds = indself(ipt, k, jj)
                  indf = indfor (ipt, k, jj)
                  indsp= inds + 1
                  indfp= indf + 1
                  
                  !$acc loop seq
                  do j = 1, ng22
                     taug(ipt, k,ns22+j, jj) = speccomb                                     &
                     &        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
                     &        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
                     &        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
                     &        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )     &
                     &        + colamt(ipt, k,1, jj) * (selffac(ipt, k, jj) * (selfref(inds,j)            &
                     &        + selffrac(ipt, k, jj) * (selfref(indsp,j)-selfref(inds,j)))       &
                     &        + forfac(ipt, k, jj) * (forref(indf,j) + forfrac(ipt, k, jj)                &
                     &        * (forref(indfp,j) - forref(indf,j)))) + o2cont
                     taur(ipt, k,ns22+j, jj) = tauray
                  enddo
               else
                  o2cont = o2tem * colamt(ipt, k,6, jj)

                  ind01 = id0(ipt, k,22, jj) + 1
                  ind02 = ind01 + 1
                  ind11 = id1(ipt, k,22, jj) + 1
                  ind12 = ind11 + 1
                  
                  !$acc loop seq
                  do j = 1, ng22
                     taug(ipt, k,ns22+j, jj) = colamt(ipt, k,6, jj) * o2adj                          &
                     &        * ( fac00(ipt, k, jj)*absb(ind01,j) + fac10(ipt, k, jj)*absb(ind02,j)       &
                     &        +   fac01(ipt, k, jj)*absb(ind11,j) + fac11(ipt, k, jj)*absb(ind12,j) )     &
                     &        + o2cont
                     taur(ipt, k,ns22+j, jj) = tauray
                  enddo
               end if
            enddo
         end do
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
     &       sfluxzen, taug, taur, async_id, fulljj, id0, id1                                        &
     &     )
     
!  ------------------------------------------------------------------  !
!     band 23:  8050-12850 cm-1 (low - h2o; high - nothing)            !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb23

      integer :: nlay, laytrop(ix, fulljj), ix, nday(fulljj)
      integer, dimension(ix, nlay, fulljj) :: indfor, indself,          &
     &       jp, jt, jt1
      real (kind=kind_phys), dimension(ix, nlay, fulljj) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac
      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj) :: colamt
      real (kind=kind_phys), dimension(ix, ngptsw, fulljj) :: sfluxzen
      real (kind=kind_phys), dimension(ix, nlay,ngptsw, fulljj) ::     &
     &       taug, taur
      integer, dimension(ix, nlay,nblow:nbhgh, fulljj) :: id0, id1
!  ---  locals:
      integer :: ipt, jj, async_id, fulljj
      integer :: ind01, ind02, ind11, ind12
      integer :: inds, indf, indsp, indfp, j, k

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !!$acc cache(absa, selfref, forref, rayl)
            !$acc loop vector private(ind01, ind02, ind11, ind12, inds, indf, &
            !$acc&     indsp, indfp)
            do ipt = 1, nday(jj) ! lab_do_ipt
               if (k .le. laytrop(ipt, jj)) then
                  ind01 = id0(ipt, k,23, jj) + 1
                  ind02 = ind01 + 1
                  ind11 = id1(ipt, k,23, jj) + 1
                  ind12 = ind11 + 1

                  inds = indself(ipt, k, jj)
                  indf = indfor (ipt, k, jj)
                  indsp= inds + 1
                  indfp= indf + 1
                  
                  !$acc loop seq
                  do j = 1, ng23
                     taug(ipt, k,ns23+j, jj) = colamt(ipt, k,1, jj) * (givfac                        &
                     &        * ( fac00(ipt, k, jj)*absa(ind01,j) + fac10(ipt, k, jj)*absa(ind02,j)       &
                     &        +   fac01(ipt, k, jj)*absa(ind11,j) + fac11(ipt, k, jj)*absa(ind12,j) )     &
                     &        + selffac(ipt, k, jj) * (selfref(inds,j) + selffrac(ipt, k, jj)             &
                     &        * (selfref(indsp,j) - selfref(inds,j)))                   &
                     &        + forfac(ipt, k, jj) * (forref(indf,j) + forfrac(ipt, k, jj)                &
                     &        * (forref(indfp,j) - forref(indf,j))))
                     taur(ipt, k,ns23+j, jj) = colmol(ipt, k, jj) * rayl(j)
                  enddo
               else
                  !$acc loop seq
                  do j = 1, ng23
                     taug(ipt, k,ns23+j, jj) = f_zero
                     taur(ipt, k,ns23+j, jj) = colmol(ipt, k, jj) * rayl(j)
                  enddo
               end if
            enddo
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
     &       sfluxzen, taug, taur, async_id, fulljj, id0, id1                                        &
     &     )
     
!  ------------------------------------------------------------------  !
!     band 24:  12850-16000 cm-1 (low - h2o,o2; high - o2)             !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb24

      integer :: nlay, laytrop(ix, fulljj), ix, nday(fulljj)
      integer, dimension(ix, nlay, fulljj) :: indfor, indself,          &
     &       jp, jt, jt1
      real (kind=kind_phys), dimension(ix, nlay, fulljj) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac
      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj) :: colamt
      real (kind=kind_phys), dimension(ix, ngptsw, fulljj) :: sfluxzen
      real (kind=kind_phys), dimension(ix, nlay,ngptsw, fulljj) ::     &
     &       taug, taur
      integer, dimension(ix, nlay,nblow:nbhgh, fulljj) :: id0, id1
!  ---  locals:
      integer :: ipt, jj, async_id, fulljj
      real (kind=kind_phys) :: speccomb, specmult, fs, fs1,             &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111

      integer :: ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14
      integer :: inds, indf, indsp, indfp, j, js, k

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !!$acc cache(rayla, absa, abso3a, selfref, forref, absb, raylb, abso3b)
            !$acc loop vector private(speccomb, specmult, js, fs, fs1, fac000, &
            !$acc&     fac010, fac100, fac110, fac001, fac011, fac101, fac111, &
            !$acc&     ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14, &
            !$acc&     inds, indf, indsp, indfp)
            do ipt = 1, nday(jj) ! lab_do_ipt
               if (k .le. laytrop(ipt, jj)) then
                  speccomb = colamt(ipt, k,1, jj) + strrat(24)*colamt(ipt, k,6, jj)
                  specmult = 8.0 * min(oneminus, colamt(ipt, k,1, jj) / speccomb)

                  js = 1 + int(specmult)
                  fs = mod(specmult, f_one)
                  fs1= f_one - fs
                  fac000 = fs1 * fac00(ipt, k, jj)
                  fac010 = fs1 * fac10(ipt, k, jj)
                  fac100 = fs  * fac00(ipt, k, jj)
                  fac110 = fs  * fac10(ipt, k, jj)
                  fac001 = fs1 * fac01(ipt, k, jj)
                  fac011 = fs1 * fac11(ipt, k, jj)
                  fac101 = fs  * fac01(ipt, k, jj)
                  fac111 = fs  * fac11(ipt, k, jj)

                  ind01 = id0(ipt, k,24, jj) + js
                  ind02 = ind01 + 1
                  ind03 = ind01 + 9
                  ind04 = ind01 + 10
                  ind11 = id1(ipt, k,24, jj) + js
                  ind12 = ind11 + 1
                  ind13 = ind11 + 9
                  ind14 = ind11 + 10

                  inds = indself(ipt, k, jj)
                  indf = indfor (ipt, k, jj)
                  indsp= inds + 1
                  indfp= indf + 1
                  
                  !$acc loop seq
                  do j = 1, ng24
                     taug(ipt, k,ns24+j, jj) = speccomb                                     &
                     &        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
                     &        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
                     &        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
                     &        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )     &
                     &        + colamt(ipt, k,3, jj) * abso3a(j) + colamt(ipt, k,1, jj)                   &
                     &        * (selffac(ipt, k, jj) * (selfref(inds,j) + selffrac(ipt, k, jj)            &
                     &        * (selfref(indsp,j) - selfref(inds,j)))                   &
                     &        + forfac(ipt, k, jj) * (forref(indf,j) + forfrac(ipt, k, jj)                &
                     &        * (forref(indfp,j) - forref(indf,j))))

                     taur(ipt, k,ns24+j, jj) = colmol(ipt, k, jj)                                    &
                     &           * (rayla(j,js) + fs*(rayla(j,js+1) - rayla(j,js)))
                  enddo
               else
                  ind01 = id0(ipt, k,24, jj) + 1
                  ind02 = ind01 + 1
                  ind11 = id1(ipt, k,24, jj) + 1
                  ind12 = ind11 + 1
                  
                  !$acc loop seq
                  do j = 1, ng24
                     taug(ipt, k,ns24+j, jj) = colamt(ipt, k,6, jj)                                  &
                     &        * ( fac00(ipt, k, jj)*absb(ind01,j) + fac10(ipt, k, jj)*absb(ind02,j)       &
                     &        +   fac01(ipt, k, jj)*absb(ind11,j) + fac11(ipt, k, jj)*absb(ind12,j) )     &
                     &        + colamt(ipt, k,3, jj) * abso3b(j)

                     taur(ipt, k,ns24+j, jj) = colmol(ipt, k, jj) * raylb(j)
                  enddo
               end if
            enddo
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
     &       sfluxzen, taug, taur, async_id, fulljj, id0, id1                                        &
     &     )
     
!  ------------------------------------------------------------------  !
!     band 25:  16000-22650 cm-1 (low - h2o; high - nothing)           !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb25

      integer :: nlay, laytrop(ix, fulljj), ix, nday(fulljj)
      integer, dimension(ix, nlay, fulljj) :: indfor, indself,          &
     &       jp, jt, jt1
      real (kind=kind_phys), dimension(ix, nlay, fulljj) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac
      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj) :: colamt
      real (kind=kind_phys), dimension(ix, ngptsw, fulljj) :: sfluxzen
      real (kind=kind_phys), dimension(ix, nlay,ngptsw, fulljj) ::     &
     &       taug, taur
      integer, dimension(ix, nlay,nblow:nbhgh, fulljj) :: id0, id1
!  ---  locals:
      integer :: ipt, jj, async_id, fulljj
      integer :: ind01, ind02, ind11, ind12
      integer :: j, k

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.
      
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !!$acc cache(absa, abso3a, rayl, abso3b)
            !$acc loop vector private(ind01, ind02, ind11, ind12)
            do ipt = 1, nday(jj) ! lab_do_ipt
               if (k .le. laytrop(ipt, jj)) then
                  ind01 = id0(ipt, k,25, jj) + 1
                  ind02 = ind01 + 1
                  ind11 = id1(ipt, k,25, jj) + 1
                  ind12 = ind11 + 1
                  
                  !$acc loop seq
                  do j = 1, ng25
                     taug(ipt, k,ns25+j, jj) = colamt(ipt, k,1, jj)                                  &
                     &        * ( fac00(ipt, k, jj)*absa(ind01,j) + fac10(ipt, k, jj)*absa(ind02,j)       &
                     &        +   fac01(ipt, k, jj)*absa(ind11,j) + fac11(ipt, k, jj)*absa(ind12,j) )     &
                     &        + colamt(ipt, k,3, jj) * abso3a(j) 
                     taur(ipt, k,ns25+j, jj) = colmol(ipt, k, jj) * rayl(j)
                  enddo
               else
                  !$acc loop seq
                  do j = 1, ng25
                     taug(ipt, k,ns25+j, jj) = colamt(ipt, k,3, jj) * abso3b(j) 
                     taur(ipt, k,ns25+j, jj) = colmol(ipt, k, jj) * rayl(j)
                  enddo
               end if
            enddo
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
     &       sfluxzen, taug, taur, async_id, fulljj, id0, id1                                        &
     &     )
     
!  ------------------------------------------------------------------  !
!     band 26:  22650-29000 cm-1 (low - nothing; high - nothing)       !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb26

      integer :: nlay, laytrop(ix, fulljj), ix, nday(fulljj)
      integer, dimension(ix, nlay, fulljj) :: indfor, indself,          &
     &       jp, jt, jt1
      real (kind=kind_phys), dimension(ix, nlay, fulljj) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac
      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj) :: colamt
      real (kind=kind_phys), dimension(ix, ngptsw, fulljj) :: sfluxzen
      real (kind=kind_phys), dimension(ix, nlay,ngptsw, fulljj) ::     &
     &       taug, taur
      integer, dimension(ix, nlay,nblow:nbhgh, fulljj) :: id0, id1
!  ---  locals:
      integer :: ipt, jj, async_id, fulljj
      integer :: j, k

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.
      !$acc parallel loop gang collapse(3) async(async_id)
      do jj = 1, fulljj
         do j = 1, ng26
            do k = 1, nlay
               !!$acc cache(rayl)
               !$acc loop vector
               do ipt = 1, nday(jj) ! lab_do_ipt
                  taug(ipt, k,ns26+j, jj) = f_zero
                  taur(ipt, k,ns26+j, jj) = colmol(ipt, k, jj) * rayl(j) 
               enddo
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
     &       sfluxzen, taug, taur, async_id, fulljj, id0, id1                                        &
     &     )
     
!  ------------------------------------------------------------------  !
!     band 27:  29000-38000 cm-1 (low - o3; high - o3)                 !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb27
!
      integer :: nlay, laytrop(ix, fulljj), ix, nday(fulljj)
      integer, dimension(ix, nlay, fulljj) :: indfor, indself,          &
     &       jp, jt, jt1
      real (kind=kind_phys), dimension(ix, nlay, fulljj) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac
      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj) :: colamt
      real (kind=kind_phys), dimension(ix, ngptsw, fulljj) :: sfluxzen
      real (kind=kind_phys), dimension(ix, nlay,ngptsw, fulljj) ::     &
     &       taug, taur
      integer, dimension(ix, nlay,nblow:nbhgh, fulljj) :: id0, id1
!  ---  locals:
      integer :: ipt, jj, async_id, fulljj
      integer :: ind01, ind02, ind11, ind12
      integer :: j, k
      real (kind=kind_phys) :: abs01, abs02, abs11, abs12

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.
      
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !!$acc cache(rayl, absa, absb)
            !$acc loop vector private(ind01, ind02, ind11, ind12, abs01, abs02, abs11, abs12)
            do ipt = 1, nday(jj) ! lab_do_ipt
               ind01 = id0(ipt, k,27, jj) + 1
               ind02 = ind01 + 1
               ind11 = id1(ipt, k,27, jj) + 1
               ind12 = ind11 + 1
                  
               !$acc loop seq
               do j = 1, ng27
                  if (k .le. laytrop(ipt, jj)) then
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
                  taug(ipt, k,ns27+j, jj) = colamt(ipt, k,3, jj)                                  &
                  &        * ( fac00(ipt, k, jj)*abs01 + fac10(ipt, k, jj)*abs02       &
                  &        +   fac01(ipt, k, jj)*abs11 + fac11(ipt, k, jj)*abs12 )
                  taur(ipt, k,ns27+j, jj) = colmol(ipt, k, jj) * rayl(j)
               enddo
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
     &       sfluxzen, taug, taur, async_id, fulljj, id0, id1                                        &
     &     )
     
!  ------------------------------------------------------------------  !
!     band 28:  38000-50000 cm-1 (low - o3,o2; high - o3,o2)           !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb28

      integer :: nlay, laytrop(ix, fulljj), ix, nday(fulljj)
      integer, dimension(ix, nlay, fulljj) :: indfor, indself,          &
     &       jp, jt, jt1
      real (kind=kind_phys), dimension(ix, nlay, fulljj) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac
      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj) :: colamt
      real (kind=kind_phys), dimension(ix, ngptsw, fulljj) :: sfluxzen
      real (kind=kind_phys), dimension(ix, nlay,ngptsw, fulljj) ::     &
     &       taug, taur
      integer, dimension(ix, nlay,nblow:nbhgh, fulljj) :: id0, id1
!  ---  locals:
      integer :: ipt, jj, async_id, fulljj
      real (kind=kind_phys) :: speccomb, specmult, tauray, fs, fs1,     &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111

      integer :: ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14
      integer :: j, js, k

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.

      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !!$acc cache(absa)
            !$acc loop vector private(speccomb, specmult, js, fs, fs1, fac000, &
            !$acc&     fac010, fac100, fac110, fac001, fac011, fac101, fac111, &
            !$acc&     ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14, tauray)
            do ipt = 1, nday(jj) ! lab_do_ipt
               tauray = colmol(ipt, k, jj) * rayl
               if (k .le. laytrop(ipt, jj)) then
                  speccomb = colamt(ipt, k,3, jj) + strrat(28)*colamt(ipt, k,6, jj)
                  specmult = 8.0 * min(oneminus, colamt(ipt, k,3, jj) / speccomb)

                  js = 1 + int(specmult)
                  fs = mod(specmult, f_one)
                  fs1= f_one - fs
                  fac000 = fs1 * fac00(ipt, k, jj)
                  fac010 = fs1 * fac10(ipt, k, jj)
                  fac100 = fs  * fac00(ipt, k, jj)
                  fac110 = fs  * fac10(ipt, k, jj)
                  fac001 = fs1 * fac01(ipt, k, jj)
                  fac011 = fs1 * fac11(ipt, k, jj)
                  fac101 = fs  * fac01(ipt, k, jj)
                  fac111 = fs  * fac11(ipt, k, jj)

                  ind01 = id0(ipt, k,28, jj) + js
                  ind02 = ind01 + 1
                  ind03 = ind01 + 9
                  ind04 = ind01 + 10
                  ind11 = id1(ipt, k,28, jj) + js
                  ind12 = ind11 + 1
                  ind13 = ind11 + 9
                  ind14 = ind11 + 10
                  
                  !$acc loop seq
                  do j = 1, ng28
                     taug(ipt, k,ns28+j, jj) = speccomb                                     &
                     &        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
                     &        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
                     &        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
                     &        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )
                     taur(ipt, k,ns28+j, jj) = tauray
                  enddo
               else
                  speccomb = colamt(ipt, k,3, jj) + strrat(28)*colamt(ipt, k,6, jj)
                  specmult = 4.0 * min(oneminus, colamt(ipt, k,3, jj) / speccomb)

                  js = 1 + int(specmult)
                  fs = mod(specmult, f_one)
                  fs1= f_one - fs
                  fac000 = fs1 * fac00(ipt, k, jj)
                  fac010 = fs1 * fac10(ipt, k, jj)
                  fac100 = fs  * fac00(ipt, k, jj)
                  fac110 = fs  * fac10(ipt, k, jj)
                  fac001 = fs1 * fac01(ipt, k, jj)
                  fac011 = fs1 * fac11(ipt, k, jj)
                  fac101 = fs  * fac01(ipt, k, jj)
                  fac111 = fs  * fac11(ipt, k, jj)

                  ind01 = id0(ipt, k,28, jj) + js
                  ind02 = ind01 + 1
                  ind03 = ind01 + 5
                  ind04 = ind01 + 6
                  ind11 = id1(ipt, k,28, jj) + js
                  ind12 = ind11 + 1
                  ind13 = ind11 + 5
                  ind14 = ind11 + 6
                  
                  !$acc loop seq
                  do j = 1, ng28
                     taug(ipt, k,ns28+j, jj) = speccomb                                     &
                     &        * ( fac000 * absb(ind01,j) + fac100 * absb(ind02,j)       &
                     &        +   fac010 * absb(ind03,j) + fac110 * absb(ind04,j)       &
                     &        +   fac001 * absb(ind11,j) + fac101 * absb(ind12,j)       &
                     &        +   fac011 * absb(ind13,j) + fac111 * absb(ind14,j) )
                     taur(ipt, k,ns28+j, jj) = tauray
                  enddo
               end if
            enddo
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
     &       sfluxzen, taug, taur, async_id, fulljj, id0, id1                                        &
     &     )
     
!  ------------------------------------------------------------------  !
!     band 29:  820-2600 cm-1 (low - h2o; high - co2)                  !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb29

      integer :: nlay, laytrop(ix, fulljj), ix, nday(fulljj)
      integer, dimension(ix, nlay, fulljj) :: indfor, indself,          &
     &       jp, jt, jt1
      real (kind=kind_phys), dimension(ix, nlay, fulljj) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac
      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj) :: colamt
      real (kind=kind_phys), dimension(ix, ngptsw, fulljj) :: sfluxzen
      real (kind=kind_phys), dimension(ix, nlay,ngptsw, fulljj) ::     &
     &       taug, taur
      integer, dimension(ix, nlay,nblow:nbhgh, fulljj) :: id0, id1
!  ---  locals:
      integer :: ipt, jj, async_id, fulljj
      real (kind=kind_phys) :: tauray

      integer :: ind01, ind02, ind11, ind12
      integer :: inds, indf, indsp, indfp, j, k

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.
      
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !$acc loop vector private(ind01, ind02, ind11, ind12, inds, indf, &
            !$acc&     indsp, indfp, tauray)
            do ipt = 1, nday(jj) ! lab_do_ipt
               !!$acc cache(absa, selfref, forref, absco2, absh2o, absb)
               tauray = colmol(ipt, k, jj) * rayl
               if (k .le. laytrop(ipt, jj)) then
                  ind01 = id0(ipt, k,29, jj) + 1
                  ind02 = ind01 + 1
                  ind11 = id1(ipt, k,29, jj) + 1
                  ind12 = ind11 + 1

                  inds = indself(ipt, k, jj)
                  indf = indfor (ipt, k, jj)
                  indsp= inds + 1
                  indfp= indf + 1
                  
                  !$acc loop seq
                  do j = 1, ng29
                     taug(ipt, k,ns29+j, jj) = colamt(ipt, k,1, jj)                                  &
                     &        * ( (fac00(ipt, k, jj)*absa(ind01,j) + fac10(ipt, k, jj)*absa(ind02,j)      &
                     &        +    fac01(ipt, k, jj)*absa(ind11,j) + fac11(ipt, k, jj)*absa(ind12,j) )    &
                     &        +  selffac(ipt, k, jj) * (selfref(inds,j) + selffrac(ipt, k, jj)            &
                     &        *  (selfref(indsp,j) - selfref(inds,j)))                  &
                     &        +  forfac(ipt, k, jj) * (forref(indf,j) + forfrac(ipt, k, jj)               &
                     &        *  (forref(indfp,j) - forref(indf,j))))                   &
                     &        +  colamt(ipt, k,2, jj) * absco2(j)
                     taur(ipt, k,ns29+j, jj) = tauray
                  enddo
               else
                  ind01 = id0(ipt, k,29, jj) + 1
                  ind02 = ind01 + 1
                  ind11 = id1(ipt, k,29, jj) + 1
                  ind12 = ind11 + 1
                  
                  !$acc loop seq
                  do j = 1, ng29
                     taug(ipt, k,ns29+j, jj) = colamt(ipt, k,2, jj)                                  &
                     &        * ( fac00(ipt, k, jj)*absb(ind01,j) + fac10(ipt, k, jj)*absb(ind02,j)       &
                     &        +   fac01(ipt, k, jj)*absb(ind11,j) + fac11(ipt, k, jj)*absb(ind12,j) )     &
                     &        + colamt(ipt, k,1, jj) * absh2o(j) 
                     taur(ipt, k,ns29+j, jj) = tauray
                  enddo
               end if
            enddo
         end do
      end do

      return
!...................................
      end subroutine taumol29
!-----------------------------------
!-----------------------------------
      subroutine spcvrtc_atomic                                                &
!...................................
!  ---  inputs:
     &     ( ssolar,cosz,sntz,albbm,albdf,sfluxzen,cldfrc,              &
     &       cf1,cf0,taug,taur,tauae,ssaae,asyae,taucw,ssacw,asycw,     &
     &       nlay, nlp1, ix, idxday, nday, async_id, fulljj,                                                &
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
      integer, parameter :: small_factor = 112
      integer, parameter :: small_ngptsw = int(ngptsw/small_factor)

!  ---  inputs:
      integer, intent(in) :: nlay, nlp1, ix, idxday(ix, fulljj), nday(fulljj), fulljj

      real (kind=kind_phys), dimension(ix, nlay,ngptsw, fulljj), intent(in) ::      &
     &       taug, taur
      real (kind=kind_phys), dimension(ix, nlay,nbdsw, fulljj),  intent(in) ::      &
     &       taucw, ssacw, asycw, tauae, ssaae, asyae

      real (kind=kind_phys), dimension(ix, ngptsw, fulljj), intent(in) :: sfluxzen
      real (kind=kind_phys), dimension(ix, nlay, fulljj),   intent(in) :: cldfrc

      real (kind=kind_phys), dimension(ix, 2, fulljj),  intent(in) :: albbm, albdf

      real (kind=kind_phys), dimension(ix, fulljj), intent(in) :: cosz, sntz, &
         cf1, cf0, ssolar

!  ---  outputs:
      real (kind=kind_phys), dimension(ix, nlp1,nbdsw, fulljj), intent(out) ::      &
     &       fxupc, fxdnc, fxup0, fxdn0

      real (kind=kind_phys), dimension(ix, 2, fulljj), intent(out) :: sfbmc, sfdfc, &
     &       sfbm0, sfdf0

      real (kind=kind_phys), dimension(ix, fulljj), intent(out) :: suvbfc, suvbf0, ftoadc,     &
     &       ftoauc, ftoau0, fsfcuc, fsfcu0, fsfcdc, fsfcd0

!  ---  locals:
      real (kind=kind_phys), dimension(ix, nlp1, small_ngptsw, fulljj) :: zrefb, zrefd, ztrab,    &
     &       ztrad, ztdbt, zldbt

      real (kind=kind_phys) :: ztau1, zssa1, zasy1, ztau0, zssa0,       &
     &       zasy0, zasy3, zssaw, zasyw, zgam1, zgam2, zgam3, zgam4,    &
     &       zc0, zc1, za1, za2, &
     &       zrk, zrk2, zrp, zrp1, zrm1, zb1, zb2, &
     &       zrpp, zrkg1, zrkg3, zrkg4, zexp1, zexm1, zexp2, zexm2,     &
     &       zexp3, zexp4, zden1, ze1r45, ftind, zrefb1,        &
     &       zrefd1, ztrab1, ztrad1, zr1, zr2, zr3, zr4, zr5, ztdbt0r,    &
     &       zt1, zt2, zt3, zf1, zf2, zldbt0, zfu, zfd, zrupbr, zrupbr1, zrupdr, zrupdr1

      real (kind=kind_phys), dimension(ix, small_ngptsw, fulljj) :: ztdbt0, zsolar, zfd0

      integer :: ib, ibd, jb, jg, k, kp, itind, ipt, jj, j1, j2, jg2
      real (kind=kind_phys), dimension(ix, nlp1, small_ngptsw, fulljj) :: zrdnd, ztdn
      integer :: async_id

!
!===> ...  begin here
!

      !$acc data create(zrefb, zrefd, ztrab, ztrad, ztdbt, zldbt, zfd0, &
      !$acc&     ztdbt0, zsolar, zrdnd, ztdn) async(async_id)
!  --- ... initialization of output fluxes
      !$acc parallel loop gang collapse(3) async(async_id)
      do jj = 1, fulljj
         do ib = 1, nbdsw
            do k = 1, nlp1
               !$acc loop vector
               do ipt = 1, nday(jj) ! lab_do_ipt
                  fxdnc(ipt, k,ib, jj) = f_zero
                  fxupc(ipt, k,ib, jj) = f_zero
                  fxdn0(ipt, k,ib, jj) = f_zero
                  fxup0(ipt, k,ib, jj) = f_zero
               enddo
            enddo
         end do
      end do
      !$acc parallel loop collapse(2) async(async_id)
      do jj = 1, fulljj
         do ipt = 1, ix
            if (ipt .le. nday(jj)) then ! lab_do_ipt
               ftoadc(ipt, jj) = f_zero
               ftoauc(ipt, jj) = f_zero
               ftoau0(ipt, jj) = f_zero
               fsfcuc(ipt, jj) = f_zero
               fsfcu0(ipt, jj) = f_zero
               fsfcdc(ipt, jj) = f_zero
               fsfcd0(ipt, jj) = f_zero

               !! --- ...  uv-b surface downward fluxes
               suvbfc(ipt, jj)  = f_zero
               suvbf0(ipt, jj)  = f_zero

               !! --- ...  output surface flux components
               sfbmc(ipt, 1, jj) = f_zero
               sfbmc(ipt, 2, jj) = f_zero
               sfdfc(ipt, 1, jj) = f_zero
               sfdfc(ipt, 2, jj) = f_zero
               sfbm0(ipt, 1, jj) = f_zero
               sfbm0(ipt, 2, jj) = f_zero
               sfdf0(ipt, 1, jj) = f_zero
               sfdf0(ipt, 2, jj) = f_zero
            end if
         end do
      end do

            !  --- ...  loop over all g-points in each band
   do j2 = 1, small_factor
      !$acc parallel loop gang collapse(2) private(jg2, jb, ib, ibd) async(async_id)
      do jj = 1, fulljj
         do jg = 1, small_ngptsw ! lab_do_jg
            jg2 = jg + (j2-1)*small_ngptsw
            jb = ngb(jg2)
            ib = jb + 1 - nblow
            ibd = idxsfc(jb)
            !$acc loop vector private(j1, kp, ztau0, zssa0, zasy0, zssaw, zasyw, &
            !$acc&     za1, za2, ztau1, zssa1, zasy1, zasy3, zgam1, zgam2, zgam3, &
            !$acc&     zgam4, zb1, zb2, ftind, itind, zrk, zrk2, zrp, zrp1, zrm1, &
            !$acc&     zrpp, zrkg1, zrkg3, zrkg4, zr1, zr2, zr3, zr4, zr5, zt1, &
            !$acc&     zt2, zt3, zexm1, zexp1, zexm2, zexp2, ze1r45, zden1, zexp3, zexp4)
            do ipt = 1, nday(jj) ! lab_do_ipt

               zsolar(ipt, jg, jj) = ssolar(ipt, jj) * sfluxzen(ipt, jg2, jj)

               !  --- ...  set up toa direct beam and surface values (beam and diff)

               ztdbt(ipt, nlp1, jg, jj) = f_one

               zldbt(ipt, 1, jg, jj) = f_zero
               if (ibd /= 0) then
                  zrefb(ipt, 1, jg, jj) = albbm(ipt, ibd, jj)
                  zrefd(ipt, 1, jg, jj) = albdf(ipt, ibd, jj)
               else
                  zrefb(ipt, 1, jg, jj) = 0.5 * (albbm(ipt, 1, jj) + albbm(ipt, 2, jj))
                  zrefd(ipt, 1, jg, jj) = 0.5 * (albdf(ipt, 1, jj) + albdf(ipt, 2, jj))
               endif
               ztrab(ipt, 1, jg, jj) = f_zero
               ztrad(ipt, 1, jg, jj) = f_zero

            !  --- ...  compute clear-sky optical parameters, layer reflectance and transmittance
               j1 = idxday(ipt, jj)
               ztdbt0r = f_one
               !$acc loop seq
               do k = nlay, 1, -1
                  kp = k + 1

                  ztau0 = max( ftiny, taur(ipt, k,jg2, jj)+taug(ipt, k,jg2, jj)+tauae(ipt, k,ib, jj) )
                  zssa0 = taur(ipt, k,jg2, jj) + tauae(ipt, k,ib, jj)*ssaae(ipt, k,ib, jj)
                  zasy0 = asyae(ipt, k,ib, jj)*ssaae(ipt, k,ib, jj)*tauae(ipt, k,ib, jj)
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

                  !  --- ...  compute homogeneous reflectance and transmittance

                  if ( zssaw >= zcrit ) then    ! for conservative scattering
                     za1 = zgam1 * cosz(j1, jj) - zgam3
                     za2 = zgam1 * ztau1

                     !  --- ...  use exponential lookup table for transmittance, or expansion
                     !           of exponential for low optical depth

                     zb1 = min ( ztau1*sntz(ipt, jj) , 500.0 )
                     if ( zb1 <= od_lo ) then
                        zb2 = f_one - zb1 + 0.5*zb1*zb1
                     else
                        ftind = zb1 / (bpade + zb1)
                        itind = ftind*ntbmx + 0.5
                        zb2 = exp_tbl(itind)
                     endif

                     !      ...  collimated beam
                     zrefb(ipt, kp, jg, jj) = max(f_zero, min(f_one,                          &
                     &                  (za2 - za1*(f_one - zb2))/(f_one + za2) ))
                     ztrab(ipt, kp, jg, jj) = max(f_zero, min(f_one, f_one-zrefb(ipt, kp, jg, jj) ))

                     !      ...  isotropic incidence
                     zrefd(ipt, kp, jg, jj) = max(f_zero, min(f_one, za2/(f_one + za2) ))
                     ztrad(ipt, kp, jg, jj) = max(f_zero, min(f_one, f_one-zrefd(ipt, kp, jg, jj) ))

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

                     zb2 = min ( sntz(ipt, jj)*ztau1, 500.0 )
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
                        zrefb(ipt, kp, jg, jj) = eps1
                        ztrab(ipt, kp, jg, jj) = zexm2
                     else
                        zden1 = zssa1 / ze1r45
                        zrefb(ipt, kp, jg, jj) = max(f_zero, min(f_one,                        &
                        &                    (zr1*zexp1 - zr2*zexm1 - zr3*zexm2)*zden1 ))
                        ztrab(ipt, kp, jg, jj) = max(f_zero, min(f_one, zexm2*(f_one           &
                        &                  - (zt1*zexp1 - zt2*zexm1 - zt3*zexp2)*zden1) ))
                     endif

                     !      ...  diffuse beam
                     zden1 = zr4 / (ze1r45 * zrkg1)
                     zrefd(ipt, kp, jg, jj) = max(f_zero, min(f_one,                          &
                     &                  zgam2*(zexp1 - zexm1)*zden1 ))
                     ztrad(ipt, kp, jg, jj) = max(f_zero, min(f_one, zrk2*zden1 ))
                  endif    ! end if_zssaw_block

                  !  --- ...  direct beam transmittance. use exponential lookup table
                  !           for transmittance, or expansion of exponential for low
                  !           optical depth

                  zr1 = ztau1 * sntz(ipt, jj)
                  if ( zr1 <= od_lo ) then
                     zexp3 = f_one - zr1 + 0.5*zr1*zr1
                  else
                     ftind = zr1 / (bpade + zr1)
                     itind = max(0, min(ntbmx, int(0.5+ntbmx*ftind) ))
                     zexp3 = exp_tbl(itind)
                  endif

                  ztdbt(ipt, k, jg, jj)  = zexp3 * ztdbt(ipt, kp, jg, jj)
                  zldbt(ipt, kp, jg, jj) = zexp3

                  !  --- ...  pre-delta-scaling clear and cloudy direct beam transmittance
                  !           (must use 'orig', unscaled cloud optical depth)

                  zr1 = ztau0 * sntz(ipt, jj)
                  if ( zr1 <= od_lo ) then
                     zexp4 = f_one - zr1 + 0.5*zr1*zr1
                  else
                     ftind = zr1 / (bpade + zr1)
                     itind = max(0, min(ntbmx, int(0.5+ntbmx*ftind) ))
                     zexp4 = exp_tbl(itind)
                  endif

                  ztdbt0r = zexp4 * ztdbt0r
               enddo    ! end do_k_loop
               ztdbt0(ipt, jg, jj) = ztdbt0r
            end do
         end do
      end do
               
               ! call swflux
               !  --- ...  link lowest layer with surface
      !$acc parallel loop gang collapse(2) private(jg2, jb, ib, ibd) async(async_id)
      do jj = 1, fulljj
         do jg = 1, small_ngptsw ! lab_do_jg
            jg2 = jg + (j2-1)*small_ngptsw
            jb = ngb(jg2)
            ib = jb + 1 - nblow
            ibd = idxsfc(jb)
            !$acc loop vector private(kp, zden1, zfu, zfd, zrupbr, zrupbr1, zrupdr, zrupdr1)
            do ipt = 1, nday(jj) ! lab_do_ipt
               !  --- ...  upper boundary conditions

               ztdn (ipt, nlp1, jg, jj) = f_one
               zrdnd(ipt, nlp1, jg, jj) = f_zero
               ztdn (ipt, nlay, jg, jj) = ztrab(ipt, nlp1, jg, jj)
               zrdnd(ipt, nlay, jg, jj) = zrefd(ipt, nlp1, jg, jj)

               !  --- ...  pass from top to bottom
               !$acc loop seq
               do k = nlay, 2, -1
                  zden1 = f_one / (f_one - zrefd(ipt, k, jg, jj)*zrdnd(ipt, k, jg, jj))
                  ztdn (ipt, k-1, jg, jj) = ztdbt(ipt, k, jg, jj)*ztrab(ipt, k, jg, jj) &
                                    + ( ztrad(ipt, k, jg, jj) *                 &
                  &                 ( (ztdn(ipt, k, jg, jj) - ztdbt(ipt, k, jg, jj)) + ztdbt(ipt, k, jg, jj) *              &
                  &                 zrefb(ipt, k, jg, jj)*zrdnd(ipt, k, jg, jj) )) * zden1
                  zrdnd(ipt, k-1, jg, jj) = zrefd(ipt, k, jg, jj) + ztrad(ipt, k, jg, jj) &
                                   *ztrad(ipt, k, jg, jj)*zrdnd(ipt, k, jg, jj)*zden1
               enddo
               
               zrupbr = zrefb(ipt, 1, jg, jj)        ! direct beam
               zrupdr = zrefd(ipt, 1, jg, jj)        ! diffused
               !  --- ...  up and down-welling fluxes at levels
               zden1 = f_one / (f_one - zrdnd(ipt, 1, jg, jj)*zrupdr)
               zfu = ( ztdbt(ipt, 1, jg, jj)*zrupbr +                                &
               &             (ztdn(ipt, 1, jg, jj) - ztdbt(ipt, 1, jg, jj)) &
                              *zrupdr ) * zden1
               zfd = ztdbt(ipt, 1, jg, jj) + ( ztdn(ipt, 1, jg, jj) &
                              - ztdbt(ipt, 1, jg, jj) +                    &
               &             ztdbt(ipt, 1, jg, jj)*zrupbr &
                              *zrdnd(ipt, 1, jg, jj) ) * zden1
               ! end call swflux
               !  --- ...  compute upward and downward fluxes at levels
               !$acc atomic
               fxup0(ipt, 1,ib, jj) = fxup0(ipt, 1,ib, jj) + zsolar(ipt, jg, jj)*zfu
               !$acc atomic
               fxdn0(ipt, 1,ib, jj) = fxdn0(ipt, 1,ib, jj) + zsolar(ipt, jg, jj)*zfd
               zfd0(ipt, jg, jj) = zfd
               !  --- ...  pass from bottom to top
               !$acc loop seq
               do k = 1, nlay
                  kp = k + 1

                  zden1 = f_one / ( f_one - zrupdr*zrefd(ipt, kp, jg, jj) )
                  zrupbr1 = zrefb(ipt, kp, jg, jj) + ( ztrad(ipt, kp, jg, jj) *                         &
                  &                ( (ztrab(ipt, kp, jg, jj) - zldbt(ipt, kp, jg, jj)) &
                                   *zrupdr +              &
                  &                zldbt(ipt, kp, jg, jj)*zrupbr) ) * zden1
                  zrupdr1 = zrefd(ipt, kp, jg, jj) + ztrad(ipt, kp, jg, jj) &
                                   *ztrad(ipt, kp, jg, jj)*zrupdr*zden1
                  
                  zden1 = f_one / (f_one - zrdnd(ipt, kp, jg, jj)*zrupdr1)
                  zfu = ( ztdbt(ipt, kp, jg, jj)*zrupbr1 +                                &
                  &             (ztdn(ipt, kp, jg, jj) - ztdbt(ipt, kp, jg, jj)) &
                                *zrupdr1 ) * zden1
                  zfd = ztdbt(ipt, kp, jg, jj) + ( ztdn(ipt, kp, jg, jj) &
                                - ztdbt(ipt, kp, jg, jj) +                    &
                  &             ztdbt(ipt, kp, jg, jj)*zrupbr1 &
                                *zrdnd(ipt, kp, jg, jj) ) * zden1
                  ! end call swflux
                  !  --- ...  compute upward and downward fluxes at levels
                  !$acc atomic
                  fxup0(ipt, kp,ib, jj) = fxup0(ipt, kp,ib, jj) + zsolar(ipt, jg, jj)*zfu
                  !$acc atomic
                  fxdn0(ipt, kp,ib, jj) = fxdn0(ipt, kp,ib, jj) + zsolar(ipt, jg, jj)*zfd
                  zrupbr = zrupbr1
                  zrupdr = zrupdr1
               end do
            end do
         end do
      end do

               !! --- ...  surface downward beam/diffused flux components
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do jg = 1, small_ngptsw ! lab_do_jg
            !$acc loop vector private(jg2, jb, ib, ibd, zf1, zf2, zb1, zb2)
            do ipt = 1, nday(jj) ! lab_do_ipt
               jg2 = jg + (j2-1)*small_ngptsw
               jb = ngb(jg2)
               ib = jb + 1 - nblow
               ibd = idxsfc(jb)
               zb1 = zsolar(ipt, jg, jj)*ztdbt0(ipt, jg, jj)
               zb2 = zsolar(ipt, jg, jj)*(zfd0(ipt, jg, jj) - ztdbt0(ipt, jg, jj))
               if (ibd /= 0) then
                  !$acc atomic
                  sfbm0(ipt, ibd, jj) = sfbm0(ipt, ibd, jj) + zb1
                  !$acc atomic
                  sfdf0(ipt, ibd, jj) = sfdf0(ipt, ibd, jj) + zb2
               else
                  zf1 = 0.5 * zb1
                  zf2 = 0.5 * zb2
                  !$acc atomic
                  sfbm0(ipt, 1, jj) = sfbm0(ipt, 1, jj) + zf1
                  !$acc atomic
                  sfdf0(ipt, 1, jj) = sfdf0(ipt, 1, jj) + zf2
                  !$acc atomic
                  sfbm0(ipt, 2, jj) = sfbm0(ipt, 2, jj) + zf1
                  !$acc atomic
                  sfdf0(ipt, 2, jj) = sfdf0(ipt, 2, jj) + zf2
               endif
               !       sfbm0(ibd) = sfbm0(ibd) + zsolar*ztdbt0
               !       sfdf0(ibd) = sfdf0(ibd) + zsolar*(zfd(1) - ztdbt0)
            end do
         end do
      end do

               !  --- ...  compute total sky optical parameters, layer reflectance and transmittance
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do jg = 1, small_ngptsw ! lab_do_jg
            !$acc loop vector private(zb1, zb2, zf1, zf2, kp, zc0, zc1, ztau0, &
            !$acc&     zssa0, zasy0, zldbt0, ftind, itind, zssaw, zasyw, za1, za2, &
            !$acc&     ztau1, zssa1, zasy1, zasy3, zgam1, zgam2, zgam3, zgam4, &
            !$acc&     zrefb1, zrefd1, ztrab1, ztrad1, zrk, zrk2, zrp, zrp1, zrm1, &
            !$acc&     zrpp, zrkg1, zrkg3, zrkg4, zr1, zr2, zr3, zr4, zr5, &
            !$acc&     zt1, zt2, zt3, zexm1, zexp1, zexm2, zexp2, ze1r45, zden1, &
            !$acc&     zexp3, zexp4, j1, jg2, jb, ib, ibd, ztdbt0r)
            do ipt = 1, nday(jj) ! lab_do_ipt
               jg2 = jg + (j2-1)*small_ngptsw
               jb = ngb(jg2)
               ib = jb + 1 - nblow
               ibd = idxsfc(jb)
               j1 = idxday(ipt, jj)
               if ( cf1(ipt, jj) > eps ) then

                  !  --- ...  set up toa direct beam and surface values (beam and diff)
                  ztdbt0r = f_one
                  zldbt(ipt, 1, jg, jj) = f_zero
                  !$acc loop seq
                  do k = nlay, 1, -1
                     kp = k + 1
                     zc0 = f_one - cldfrc(ipt, k, jj)
                     zc1 = cldfrc(ipt, k, jj)

                     !  --- ...  saving clear-sky quantities for later total-sky usage
                     ztau0 = max( ftiny, taur(ipt, k,jg2, jj)+taug(ipt, k,jg2, jj)+tauae(ipt, k,ib, jj) )
                     zssa0 = taur(ipt, k,jg2, jj) + tauae(ipt, k,ib, jj)*ssaae(ipt, k,ib, jj)
                     zasy0 = asyae(ipt, k,ib, jj)*ssaae(ipt, k,ib, jj)*tauae(ipt, k,ib, jj)
                     zr1 = ztau0 * sntz(ipt, jj)
                     if ( zr1 <= od_lo ) then
                        zldbt0 = f_one - zr1 + 0.5*zr1*zr1
                     else
                        ftind = zr1 / (bpade + zr1)
                        itind = max(0, min(ntbmx, int(0.5+ntbmx*ftind) ))
                        zldbt0 = exp_tbl(itind)
                     endif
                     
                     if ( zc1 > ftiny ) then          ! it is a cloudy-layer
                        ztau0 = ztau0 + taucw(ipt, k,ib, jj)
                        zssa0 = zssa0 + ssacw(ipt, k,ib, jj)
                        zasy0 = zasy0 + asycw(ipt, k,ib, jj)
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

                        zrefb1 = zrefb(ipt, kp, jg, jj)
                        zrefd1 = zrefd(ipt, kp, jg, jj)
                        ztrab1 = ztrab(ipt, kp, jg, jj)
                        ztrad1 = ztrad(ipt, kp, jg, jj)

                        !  --- ...  compute homogeneous reflectance and transmittance

                        if ( zssaw >= zcrit ) then    ! for conservative scattering
                           za1 = zgam1 * cosz(j1, jj) - zgam3
                           za2 = zgam1 * ztau1

                           !  --- ...  use exponential lookup table for transmittance, or expansion
                           !           of exponential for low optical depth

                           zb1 = min ( ztau1*sntz(ipt, jj) , 500.0 )
                           if ( zb1 <= od_lo ) then
                              zb2 = f_one - zb1 + 0.5*zb1*zb1
                           else
                              ftind = zb1 / (bpade + zb1)
                              itind = ftind*ntbmx + 0.5
                              zb2 = exp_tbl(itind)
                           endif

                           !      ...  collimated beam
                           zrefb(ipt, kp, jg, jj) = max(f_zero, min(f_one,                      &
                           &                      (za2 - za1*(f_one - zb2))/(f_one + za2) ))
                           ztrab(ipt, kp, jg, jj) = max(f_zero, min(f_one, f_one-zrefb(ipt, kp, jg, jj)))

                           !      ...  isotropic incidence
                           zrefd(ipt, kp, jg, jj) = max(f_zero, min(f_one, za2 / (f_one+za2) ))
                           ztrad(ipt, kp, jg, jj) = max(f_zero, min(f_one, f_one - zrefd(ipt, kp, jg, jj) ))

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

                           zb2 = min ( ztau1*sntz(ipt, jj), 500.0 )
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
                              zrefb(ipt, kp, jg, jj) = eps1
                              ztrab(ipt, kp, jg, jj) = zexm2
                           else
                              zden1 = zssa1 / ze1r45
                              zrefb(ipt, kp, jg, jj) = max(f_zero, min(f_one,                    &
                              &                        (zr1*zexp1-zr2*zexm1-zr3*zexm2)*zden1 ))
                              ztrab(ipt, kp, jg, jj) = max(f_zero, min(f_one, zexm2*(f_one -     &
                              &                        (zt1*zexp1-zt2*zexm1-zt3*zexp2)*zden1) ))
                           endif

                           !      ...  diffuse beam
                           zden1 = zr4 / (ze1r45 * zrkg1)
                           zrefd(ipt, kp, jg, jj) = max(f_zero, min(f_one,                      &
                           &                      zgam2*(zexp1 - zexm1)*zden1 ))
                           ztrad(ipt, kp, jg, jj) = max(f_zero, min(f_one, zrk2*zden1 ))
                        endif    ! end if_zssaw_block

                        !  --- ...  combine clear and cloudy contributions for total sky
                        !           and calculate direct beam transmittances

                        zrefb(ipt, kp, jg, jj) = zc0*zrefb1 + zc1*zrefb(ipt, kp, jg, jj)
                        zrefd(ipt, kp, jg, jj) = zc0*zrefd1 + zc1*zrefd(ipt, kp, jg, jj)
                        ztrab(ipt, kp, jg, jj) = zc0*ztrab1 + zc1*ztrab(ipt, kp, jg, jj)
                        ztrad(ipt, kp, jg, jj) = zc0*ztrad1 + zc1*ztrad(ipt, kp, jg, jj)

                        !  --- ...  direct beam transmittance. use exponential lookup table
                        !           for transmittance, or expansion of exponential for low
                        !           optical depth

                        zr1 = ztau1 * sntz(ipt, jj)
                        if ( zr1 <= od_lo ) then
                           zexp3 = f_one - zr1 + 0.5*zr1*zr1
                        else
                           ftind = zr1 / (bpade + zr1)
                           itind = max(0, min(ntbmx, int(0.5+ntbmx*ftind) ))
                           zexp3 = exp_tbl(itind)
                        endif

                        zldbt(ipt, kp, jg, jj) = zc0*zldbt(ipt, kp, jg, jj) + zc1*zexp3
                        ztdbt(ipt, k, jg, jj) = zldbt(ipt, kp, jg, jj) * ztdbt(ipt, kp, jg, jj)

                        !  --- ...  pre-delta-scaling clear and cloudy direct beam transmittance
                        !           (must use 'orig', unscaled cloud optical depth)

                        zr1 = ztau0 * sntz(ipt, jj)
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
                        ztdbt(ipt, k, jg, jj) = zldbt(ipt, kp, jg, jj) * ztdbt(ipt, kp, jg, jj)

                        !  --- ...  pre-delta-scaling clear and cloudy direct beam transmittance
                        ztdbt0r = zldbt0 * ztdbt0r

                     endif    ! end if_zc1_block
                  enddo   ! end do_k_loop
                  ztdbt0(ipt, jg, jj) = ztdbt0r
               end if
            end do
         end do
      end do
                  !  --- ...  perform vertical quadrature
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do jg = 1, small_ngptsw ! lab_do_jg
            !$acc loop vector private(j1, kp, zden1, zfu, zfd, jg2, jb, ib, zrupbr, zrupbr1, zrupdr, zrupdr1)
            do ipt = 1, nday(jj) ! lab_do_ipt
               jg2 = jg + (j2-1)*small_ngptsw
               jb = ngb(jg2)
               ib = jb + 1 - nblow
               if ( cf1(ipt, jj) > eps ) then
                  ! call swflux

                  !  --- ...  upper boundary conditions

                  ztdn (ipt, nlp1, jg, jj) = f_one
                  zrdnd(ipt, nlp1, jg, jj) = f_zero
                  ztdn (ipt, nlay, jg, jj) = ztrab(ipt, nlp1, jg, jj)
                  zrdnd(ipt, nlay, jg, jj) = zrefd(ipt, nlp1, jg, jj)

                  !  --- ...  pass from top to bottom
                  !$acc loop seq
                  do k = nlay, 2, -1
                     zden1 = f_one / (f_one - zrefd(ipt, k, jg, jj)*zrdnd(ipt, k, jg, jj))
                     ztdn (ipt, k-1, jg, jj) = ztdbt(ipt, k, jg, jj)*ztrab(ipt, k, jg, jj) &
                                       + ( ztrad(ipt, k, jg, jj) *                 &
                     &                 ( (ztdn(ipt, k, jg, jj) - ztdbt(ipt, k, jg, jj)) &
                                       + ztdbt(ipt, k, jg, jj) *              &
                     &                 zrefb(ipt, k, jg, jj)*zrdnd(ipt, k, jg, jj) )) * zden1
                     zrdnd(ipt, k-1, jg, jj) = zrefd(ipt, k, jg, jj) + ztrad(ipt, k, jg, jj) &
                                       *ztrad(ipt, k, jg, jj)*zrdnd(ipt, k, jg, jj)*zden1
                  enddo

                  !  --- ...  link lowest layer with surface

                  zrupbr = zrefb(ipt, 1, jg, jj)        ! direct beam
                  zrupdr = zrefd(ipt, 1, jg, jj)        ! diffused
                  
                  zden1 = f_one / (f_one - zrdnd(ipt, 1, jg, jj)*zrupdr)
                  zfu = ( ztdbt(ipt, 1, jg, jj)*zrupbr +                                &
                  &             (ztdn(ipt, 1, jg, jj) - ztdbt(ipt, 1, jg, jj)) &
                                 *zrupdr ) * zden1
                  zfd = ztdbt(ipt, 1, jg, jj) + ( ztdn(ipt, 1, jg, jj) &
                                 - ztdbt(ipt, 1, jg, jj) +                    &
                  &             ztdbt(ipt, 1, jg, jj)*zrupbr &
                                 *zrdnd(ipt, 1, jg, jj) ) * zden1
                                 ! end call swflux

                  !  --- ...  compute upward and downward fluxes at levels
                  !$acc atomic
                  fxupc(ipt, 1,ib, jj) = fxupc(ipt, 1,ib, jj) + zsolar(ipt, jg, jj)*zfu
                  !$acc atomic
                  fxdnc(ipt, 1,ib, jj) = fxdnc(ipt, 1,ib, jj) + zsolar(ipt, jg, jj)*zfd
                  zfd0(ipt, jg, jj) = zfd
                  !  --- ...  pass from bottom to top
                  !$acc loop seq
                  do k = 1, nlay
                     kp = k + 1

                     zden1 = f_one / ( f_one - zrupdr*zrefd(ipt, kp, jg, jj) )
                     zrupbr1 = zrefb(ipt, kp, jg, jj) + ( ztrad(ipt, kp, jg, jj) *                         &
                     &                ( (ztrab(ipt, kp, jg, jj) - zldbt(ipt, kp, jg, jj))*zrupdr +              &
                     &                zldbt(ipt, kp, jg, jj)*zrupbr) ) * zden1
                     zrupdr1 = zrefd(ipt, kp, jg, jj) + ztrad(ipt, kp, jg, jj) &
                                      *ztrad(ipt, kp, jg, jj)*zrupdr*zden1

                  !  --- ...  up and down-welling fluxes at levels
                     zden1 = f_one / (f_one - zrdnd(ipt, kp, jg, jj)*zrupdr1)
                     zfu = ( ztdbt(ipt, kp, jg, jj)*zrupbr1 +                                &
                     &             (ztdn(ipt, kp, jg, jj) - ztdbt(ipt, kp, jg, jj)) &
                                    *zrupdr1 ) * zden1
                     zfd = ztdbt(ipt, kp, jg, jj) + ( ztdn(ipt, kp, jg, jj) &
                                    - ztdbt(ipt, kp, jg, jj) +                    &
                     &             ztdbt(ipt, kp, jg, jj)*zrupbr1 &
                                    *zrdnd(ipt, kp, jg, jj) ) * zden1
                                    ! end call swflux

                     !  --- ...  compute upward and downward fluxes at levels
                     !$acc atomic
                     fxupc(ipt, kp,ib, jj) = fxupc(ipt, kp,ib, jj) + zsolar(ipt, jg, jj)*zfu
                     !$acc atomic
                     fxdnc(ipt, kp,ib, jj) = fxdnc(ipt, kp,ib, jj) + zsolar(ipt, jg, jj)*zfd
                     zrupbr = zrupbr1
                     zrupdr = zrupdr1
                  end do
               end if
            end do
         end do
      end do

                  !! --- ...  surface downward beam/diffused flux components
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do jg = 1, small_ngptsw ! lab_do_jg
            !$acc loop vector private(jg2, jb, ibd, zb1, zb2, zf1, zf2)
            do ipt = 1, nday(jj) ! lab_do_ipt
               jg2 = jg + (j2-1)*small_ngptsw
               jb = ngb(jg2)
               ibd = idxsfc(jb)
               if ( cf1(ipt, jj) > eps ) then
                  zb1 = zsolar(ipt, jg, jj)*ztdbt0(ipt, jg, jj)
                  zb2 = zsolar(ipt, jg, jj)*(zfd0(ipt, jg, jj) - ztdbt0(ipt, jg, jj))

                  if (ibd /= 0) then
                     !$acc atomic
                     sfbmc(ipt, ibd, jj) = sfbmc(ipt, ibd, jj) + zb1
                     !$acc atomic
                     sfdfc(ipt, ibd, jj) = sfdfc(ipt, ibd, jj) + zb2
                  else
                     zf1 = 0.5 * zb1
                     zf2 = 0.5 * zb2
                     !$acc atomic
                     sfbmc(ipt, 1, jj) = sfbmc(ipt, 1, jj) + zf1
                     !$acc atomic
                     sfdfc(ipt, 1, jj) = sfdfc(ipt, 1, jj) + zf2
                     !$acc atomic
                     sfbmc(ipt, 2, jj) = sfbmc(ipt, 2, jj) + zf1
                     !$acc atomic
                     sfdfc(ipt, 2, jj) = sfdfc(ipt, 2, jj) + zf2
                  endif
                  !         sfbmc(ibd) = sfbmc(ibd) + zsolar*ztdbt0
                  !         sfdfc(ibd) = sfdfc(ibd) + zsolar*(zfd(1) - ztdbt0)

               endif      ! end if_cf1_block

            enddo  ! lab_do_jg
         end do
      end do
   end do


            !  --- ...  end of g-point loop
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do ib = 1, nbdsw
            !$acc loop vector 
            do ipt = 1, nday(jj) ! lab_do_ipt
               !$acc atomic
               ftoadc(ipt, jj) = ftoadc(ipt, jj) + fxdn0(ipt, nlp1,ib, jj)
               !$acc atomic
               ftoau0(ipt, jj) = ftoau0(ipt, jj) + fxup0(ipt, nlp1,ib, jj)
               !$acc atomic
               fsfcu0(ipt, jj) = fsfcu0(ipt, jj) + fxup0(ipt, 1,ib, jj)
               !$acc atomic
               fsfcd0(ipt, jj) = fsfcd0(ipt, jj) + fxdn0(ipt, 1,ib, jj)
            enddo
         end do
      end do

      !! --- ...  uv-b surface downward flux
      !$acc parallel loop gang collapse(3) async(async_id)
      do jj = 1, fulljj
         do ib = 1, nbdsw
            do k = 1, nlp1
               !$acc loop vector
               do ipt = 1, nday(jj) ! lab_do_ipt
                  if ( cf1(ipt, jj) <= eps ) then       ! clear column, set total-sky=clear-sky fluxes
                     fxupc(ipt, k,ib, jj) = fxup0(ipt, k,ib, jj)
                     fxdnc(ipt, k,ib, jj) = fxdn0(ipt, k,ib, jj)
                  end if
               enddo
            enddo
         end do
      end do
      !$acc parallel loop collapse(2) private(ibd) async(async_id)
      do jj = 1, fulljj
         do ipt = 1, ix
            if (ipt .le. nday(jj)) then ! lab_do_ipt
               ibd = nuvb - nblow + 1
               suvbf0(ipt, jj) = fxdn0(ipt, 1,ibd, jj)
               if ( cf1(ipt, jj) <= eps ) then       ! clear column, set total-sky=clear-sky fluxes
                  ftoauc(ipt, jj) = ftoau0(ipt, jj)
                  fsfcuc(ipt, jj) = fsfcu0(ipt, jj)
                  fsfcdc(ipt, jj) = fsfcd0(ipt, jj)

                  !! --- ...  surface downward beam/diffused flux components
                  sfbmc(ipt, 1, jj) = sfbm0(ipt, 1, jj)
                  sfdfc(ipt, 1, jj) = sfdf0(ipt, 1, jj)
                  sfbmc(ipt, 2, jj) = sfbm0(ipt, 2, jj)
                  sfdfc(ipt, 2, jj) = sfdf0(ipt, 2, jj)

                  !! --- ...  uv-b surface downward flux
                  suvbfc(ipt, jj) = suvbf0(ipt, jj)
               end if
            end if
         end do
      end do
      !$acc parallel loop gang collapse(3) async(async_id)
      do jj = 1, fulljj
         do ib = 1, nbdsw
            do k = 1, nlp1
               !$acc loop vector 
               do ipt = 1, nday(jj) ! lab_do_ipt
                  if ( cf1(ipt, jj) > eps ) then                        ! cloudy column, compute total-sky fluxes
                     fxupc(ipt, k,ib, jj) = cf1(ipt, jj)*fxupc(ipt, k,ib, jj) + cf0(ipt, jj)*fxup0(ipt, k,ib, jj)
                     fxdnc(ipt, k,ib, jj) = cf1(ipt, jj)*fxdnc(ipt, k,ib, jj) + cf0(ipt, jj)*fxdn0(ipt, k,ib, jj)
                  end if
               enddo
            enddo
         end do
      end do
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do ib = 1, nbdsw
            !$acc loop vector
            do ipt = 1, nday(jj) ! lab_do_ipt
               if ( cf1(ipt, jj) > eps ) then                        ! cloudy column, compute total-sky fluxes
                  !$acc atomic
                  ftoauc(ipt, jj) = ftoauc(ipt, jj) + fxupc(ipt, nlp1,ib, jj)
                  !$acc atomic
                  fsfcuc(ipt, jj) = fsfcuc(ipt, jj) + fxupc(ipt, 1,ib, jj)
                  !$acc atomic
                  fsfcdc(ipt, jj) = fsfcdc(ipt, jj) + fxdnc(ipt, 1,ib, jj)
               end if
            enddo
         end do
      end do
      !$acc parallel loop collapse(2) private(ibd) async(async_id)
      do jj = 1, fulljj
         do ipt = 1, ix
            if (ipt .le. nday(jj)) then ! lab_do_ipt
               if ( cf1(ipt, jj) > eps ) then                        ! cloudy column, compute total-sky fluxes
                  !! --- ...  uv-b surface downward flux
                  ibd = nuvb - nblow + 1
                  suvbfc(ipt, jj) = fxdnc(ipt, 1,ibd, jj)

                  !! --- ...  surface downward beam/diffused flux components
                  sfbmc(ipt, 1, jj) = cf1(ipt, jj)*sfbmc(ipt, 1, jj) + cf0(ipt, jj)*sfbm0(ipt, 1, jj)
                  sfbmc(ipt, 2, jj) = cf1(ipt, jj)*sfbmc(ipt, 2, jj) + cf0(ipt, jj)*sfbm0(ipt, 2, jj)
                  sfdfc(ipt, 1, jj) = cf1(ipt, jj)*sfdfc(ipt, 1, jj) + cf0(ipt, jj)*sfdf0(ipt, 1, jj)
                  sfdfc(ipt, 2, jj) = cf1(ipt, jj)*sfdfc(ipt, 2, jj) + cf0(ipt, jj)*sfdf0(ipt, 2, jj)
               endif    ! end if_cf1_block
            end if
         end do
      end do
      !$acc end data

      return
!...................................
      end subroutine spcvrtc_atomic
!-----------------------------------


!
!........................................!
      end module module_radsw_main_gpu       !
!========================================!

