!!!!!  ==============================================================  !!!!!
!!!!!               lw-rrtm3 radiation package description             !!!!!
!!!!!  ==============================================================  !!!!!
!                                                                          !
!   this package includes ncep's modifications of the rrtm-lw radiation    !
!   code from aer inc.                                                     !
!                                                                          !
!    the lw-rrtm3 package includes these parts:                            !
!                                                                          !
!       'radlw_rrtm3_param.f'                                              !
!       'radlw_rrtm3_datatb.f'                                             !
!       'radlw_rrtm3_main.f'                                               !
!                                                                          !
!    the 'radlw_rrtm3_param.f' contains:                                   !
!                                                                          !
!       'module_radlw_parameters'  -- band parameters set up               !
!                                                                          !
!    the 'radlw_rrtm3_datatb.f' contains:                                  !
!                                                                          !
!       'module_radlw_avplank'     -- plank flux data                      !
!       'module_radlw_ref'         -- reference temperature and pressure   !
!       'module_radlw_cldprlw'     -- cloud property coefficients          !
!       'module_radlw_kgbnn'       -- absorption coeffients for 16         !
!                                     bands, where nn = 01-16              !
!                                                                          !
!    the 'radlw_rrtm3_main.f' contains:                                    !
!                                                                          !
!       'module_radlw_main'        -- main lw radiation transfer           !
!                                                                          !
!    in the main module 'module_radlw_main' there are only two             !
!    externally callable subroutines:                                      !
!                                                                          !
!                                                                          !
!       'lwrad'     -- main lw radiation routine                           !
!          inputs:                                                         !
!           (plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                         !
!            clouds,icseed,aerosols,sfemis,sfgtmp,                         !
!            npts, nlay, nlp1, lprnt,                                      !
!          outputs:                                                        !
!            hlwc,topflx,sfcflx,                                           !
!!         optional outputs:                                               !
!            hlw0,hlwb,flxprf)                                             !
!                                                                          !
!       'rlwinit'   -- initialization routine                              !
!          inputs:                                                         !
!           ( me )                                                         !
!          outputs:                                                        !
!           (none)                                                         !
!                                                                          !
!    all the lw radiation subprograms become contained subprograms         !
!    in module 'module_radlw_main' and many of them are not directly       !
!    accessable from places outside the module.                            !
!                                                                          !
!    derived data type constructs used:                                    !
!                                                                          !
!     1. radiation flux at toa: (from module 'module_radlw_parameters')    !
!          topflw_type   -  derived data type for toa rad fluxes           !
!            upfxc              total sky upward flux at toa               !
!            upfx0              clear sky upward flux at toa               !
!                                                                          !
!     2. radiation flux at sfc: (from module 'module_radlw_parameters')    !
!          sfcflw_type   -  derived data type for sfc rad fluxes           !
!            upfxc              total sky upward flux at sfc               !
!            upfx0              clear sky upward flux at sfc               !
!            dnfxc              total sky downward flux at sfc             !
!            dnfx0              clear sky downward flux at sfc             !
!                                                                          !
!     3. radiation flux profiles(from module 'module_radlw_parameters')    !
!          proflw_type    -  derived data type for rad vertical prof       !
!            upfxc              level upward flux for total sky            !
!            dnfxc              level downward flux for total sky          !
!            upfx0              level upward flux for clear sky            !
!            dnfx0              level downward flux for clear sky          !
!                                                                          !
!    external modules referenced:                                          !
!                                                                          !
!       'module physpara'                                                  !
!       'module physcons'                                                  !
!       'mersenne_twister'                                                 !
!                                                                          !
!    compilation sequence is:                                              !
!                                                                          !
!       'radlw_rrtm3_param.f'                                              !
!       'radlw_rrtm3_datatb.f'                                             !
!       'radlw_rrtm3_main.f'                                               !
!                                                                          !
!    and all should be put in front of routines that use lw modules        !
!                                                                          !
!==========================================================================!
!                                                                          !
!    the original aer's program declarations:                              !
!                                                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                          |
!  copyright 2002-2007, atmospheric & environmental research, inc. (aer).  |
!  this software may be used, copied, or redistributed as long as it is    |
!  not sold and this copyright notice is reproduced on each copy made.     |
!  this model is provided as is without any express or implied warranties. |
!                       (http://www.rtweb.aer.com/)                        |
!                                                                          |
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                          !
! ************************************************************************ !
!                                                                          !
!                              rrtmg_lw                                    !
!                                                                          !
!                                                                          !
!                   a rapid radiative transfer model                       !
!                       for the longwave region                            ! 
!             for application to general circulation models                !
!                                                                          !
!                                                                          !
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
!        following people:  steven j. taubman, karen cady-pereira,         !
!        patrick d. brown, ronald e. farren, luke chen, robert bergstrom.  !
!                                                                          !
! ************************************************************************ !
!                                                                          !
!    references:                                                           !
!    (rrtm_lw/rrtmg_lw):                                                   !
!      clough, s.a., m.w. shephard, e.j. mlawer, j.s. delamere,            !
!      m.j. iacono, k. cady-pereira, s. boukabara, and p.d. brown:         !
!      atmospheric radiative transfer modeling: a summary of the aer       !
!      codes, j. quant. spectrosc. radiat. transfer, 91, 233-244, 2005.    !
!                                                                          !
!      mlawer, e.j., s.j. taubman, p.d. brown, m.j. iacono, and s.a.       !
!      clough:  radiative transfer for inhomogeneous atmospheres: rrtm,    !
!      a validated correlated-k model for the longwave.  j. geophys. res., !
!      102, 16663-16682, 1997.                                             !
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
!     this version of rrtmg_lw has been modified from rrtm_lw to use a     !
!     reduced set of g-points for application to gcms.                     !
!                                                                          !
! --  original version (derived from rrtm_lw), reduction of g-points,      !
!     other revisions for use with gcms.                                   !
!        1999: m. j. iacono, aer, inc.                                     !
! --  adapted for use with ncar/cam3.                                      !
!        may 2004: m. j. iacono, aer, inc.                                 !
! --  revised to add mcica capability.                                     !
!        nov 2005: m. j. iacono, aer, inc.                                 !
! --  conversion to f90 formatting for consistency with rrtmg_sw.          !
!        feb 2007: m. j. iacono, aer, inc.                                 !
! --  modifications to formatting to use assumed-shape arrays.             !
!        aug 2007: m. j. iacono, aer, inc.                                 !
!                                                                          !
! ************************************************************************ !
!                                                                          !
!    ncep modifications history log:                                       !
!                                                                          !
!       nov 1999,  ken campana       -- received the original code from    !
!                    aer (1998 ncar ccm version), updated to link up with  !
!                    ncep mrf model                                        !
!       jun 2000,  ken campana       -- added option to switch random and  !
!                    maximum/random cloud overlap                          !
!           2001,  shrinivas moorthi -- further updates for mrf model      !
!       may 2001,  yu-tai hou        -- updated on trace gases and cloud   !
!                    property based on rrtm_v3.0 codes.                    !
!       dec 2001,  yu-tai hou        -- rewritten code into fortran 90 std !
!                    set ncep radiation structure standard that contains   !
!                    three plug-in compatable fortran program files:       !
!                    'radlw_param.f', 'radlw_datatb.f', 'radlw_main.f'     !
!                    fixed bugs in subprograms taugb14, taugb2, etc. added !
!                    out-of-bounds protections. (a detailed note of        !
!                    up_to_date modifications/corrections by ncep was sent !
!                    to aer in 2002)                                       !
!       jun 2004,  yu-tai hou        -- added mike iacono's apr 2004       !
!                    modification of variable diffusivity angles.          !
!       apr 2005,  yu-tai hou        -- minor modifications on module      !
!                    structures include rain/snow effect (this version of  !
!                    code was given back to aer in jun 2006)               !
!       mar 2007,  yu-tai hou        -- added aerosol effect for ncep      !
!                    models using the generallized aerosol optical property!
!                    scheme for gfs model.                                 !
!       apr 2007,  yu-tai hou        -- added spectral band heating as an  !
!                    optional output to support the 500 km gfs model's     !
!                    upper stratospheric radiation calculations. and       !
!                    restructure optional outputs for easy access by       !
!                    different models.                                     !
!       oct 2008,  yu-tai hou        -- modified to include new features   !
!                    from aer's newer release v4.4-v4.7, including the     !
!                    mcica sub-grid cloud option. add rain/snow optical    !
!                    properties support to cloudy sky calculations.        !
!                    correct errors in mcica cloud optical properties for  !
!                    ebert & curry scheme (ilwcice=1) that needs band      !
!                    index conversion. simplified and unified sw and lw    !
!                    sub-column cloud subroutines into one module by using !
!                    optional parameters.                                  !
!       mar 2009,  yu-tai hou        -- replaced the original random number!
!                    generator coming from the original code with ncep w3  !
!                    library to simplify the program and moved sub-column  !
!                    cloud subroutines inside the main module. added       !
!                    option of user provided permutation seeds that could  !
!                    be randomly generated from forecast time stamp.       !
!       oct 2009,  yu-tai hou        -- modified subrtines "cldprop" and   !
!                    "rlwinit" according updats from aer's rrtmg_lw v4.8.  !
!       nov 2009,  yu-tai hou        -- modified subrtine "taumol" according
!                    updats from aer's rrtmg_lw version 4.82. notice the   !
!                    cloud ice/liquid are assumed as in-cloud quantities,  !
!                    not as grid averaged quantities.                      !
!       jun 2010,  yu-tai hou        -- optimized code to improve efficiency
!       apr 2012,  b. ferrier and y. hou -- added conversion factor to fu's!
!                    cloud-snow optical property scheme.                   !
!       nov 2012,  yu-tai hou        -- modified control parameters thru   !
!                     module 'physpara'.                                   !  
!                                                                          !
!!!!!  ==============================================================  !!!!!
!!!!!                         end descriptions                         !!!!!
!!!!!  ==============================================================  !!!!!


!========================================!
      module module_radlw_main_gpu           !
!........................................!
!
      use physpara,         only : ilwrate, ilwrgas, ilwcliq, ilwcice,  &
     &                             isubclw, icldflg, iovrlw,  ivflip,   &
     &                             kind_phys
      use physcons,         only : con_g, con_cp, con_avgd, con_amd,    &
     &                             con_amw, con_amo3
      use mersenne_twister, only : random_setseed, random_number,       &
     &                             random_stat
     use module_radiation_aerosols_gpu,only : setaer_lw_gpu
      use module_radlw_parameters
!
      use module_radlw_avplank, only : totplnk
      use module_radlw_ref,     only : preflog, tref, chi_mls
      use param,             only : my
      use index,             only : jlistnum, nxptot, nxjp_acc
      !use nvtx
      use rank, only: myrank
!
      implicit none
!
      private
!
!  ...  version tag and last revision date
      character(40), parameter ::                                       &
     &   vtaglw='ncep lw v5.1  nov 2012 -rrtmg-lw v4.82  '
!    &   vtaglw='ncep lw v5.0  aug 2012 -rrtmg-lw v4.82  '
!    &   vtaglw='rrtmg-lw v4.82  nov 2009  '
!    &   vtaglw='rrtmg-lw v4.8   oct 2009  '
!    &   vtaglw='rrtmg-lw v4.71  mar 2009  '
!    &   vtaglw='rrtmg-lw v4.4   oct 2008  '
!    &   vtaglw='rrtm-lw v2.3g   mar 2007  '
!    &   vtaglw='rrtm-lw v2.3g   apr 2004  '

!  ---  constant values
      real (kind=kind_phys), parameter :: eps     = 1.0e-6
      real (kind=kind_phys), parameter :: oneminus= 1.0-eps
      real (kind=kind_phys), parameter :: cldmin  = 1.0e-80
      real (kind=kind_phys), parameter :: bpade   = 1.0/0.278  ! pade approx constant
      real (kind=kind_phys), parameter :: stpfac  = 296.0/1013.0
      real (kind=kind_phys), parameter :: wtdiff  = 0.5        ! weight for radiance to flux conversion
      real (kind=kind_phys), parameter :: tblint  = ntbl       ! lookup table conversion factor
      real (kind=kind_phys), parameter :: f_zero  = 0.0
      real (kind=kind_phys), parameter :: f_one   = 1.0

!  ...  atomic weights for conversion from mass to volume mixing ratios
      real (kind=kind_phys), parameter :: amdw    = con_amd/con_amw
      real (kind=kind_phys), parameter :: amdo3   = con_amd/con_amo3

!  ...  band indices
      integer, dimension(nbands) :: nspa, nspb

      data nspa / 1, 1, 9, 9, 9, 1, 9, 1, 9, 1, 1, 9, 9, 1, 9, 9 /
      data nspb / 1, 1, 5, 5, 5, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0 /

!  ...  band wavenumber intervals
!     real (kind=kind_phys) :: wavenum1(nbands), wavenum2(nbands)
!     data wavenum1/                                                    &
!    &         10.,  350.,  500.,  630.,  700.,  820.,  980., 1080.,    &
!err &       1180., 1390., 1480., 1800., 2080., 2250., 2390., 2600. /
!    &       1180., 1390., 1480., 1800., 2080., 2250., 2380., 2600. /
!     data wavenum2/                                                    &
!    &        350.,  500.,  630.,  700.,  820.,  980., 1080., 1180.,    &
!err &       1390., 1480., 1800., 2080., 2250., 2390., 2600., 3250. /
!    &       1390., 1480., 1800., 2080., 2250., 2380., 2600., 3250. /
!     real (kind=kind_phys) :: delwave(nbands)
!     data delwave / 340., 150., 130.,  70., 120., 160., 100., 100.,    &
!    &               210.,  90., 320., 280., 170., 130., 220., 650. /

!  ---  reset diffusivity angle for bands 2-3 and 5-9 to vary (between 1.50
!       and 1.80) as a function of total column water vapor.  the function
!       has been defined to minimize flux and cooling rate errors in these bands
!       over a wide range of precipitable water values.
      real (kind=kind_phys), dimension(nbands) :: a0, a1, a2

      data a0 / 1.66,  1.55,  1.58,  1.66,  1.54, 1.454,  1.89,  1.33,  &
     &         1.668,  1.66,  1.66,  1.66,  1.66,  1.66,  1.66,  1.66 /
      data a1 / 0.00,  0.25,  0.22,  0.00,  0.13, 0.446, -0.10,  0.40,  &
     &        -0.006,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00 /
      data a2 / 0.00, -12.0, -11.7,  0.00, -0.72,-0.243,  0.19,-0.062,  &
     &         0.414,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00 /

!! ---  logical flags for optional output fields

      logical :: lhlwb  = .false.
      logical :: lhlw0  = .false.
      logical :: lflxprf= .false.
      logical :: lflxprf_upfxc= .false.
      logical :: lflxprf_dnfxc= .false.
      logical :: lflxprf_upfx0= .false.
      logical :: lflxprf_dnfx0= .false.

!  ---  those data will be set up only once by "rlwinit"

!  ...  fluxfac, heatfac are factors for fluxes (in w/m**2) and heating
!       rates (in k/day, or k/sec set by subroutine 'rlwinit')
!       semiss0 are default surface emissivity for each bands

      real (kind=kind_phys) :: fluxfac, heatfac, semiss0(nbands)
      data semiss0(:) / nbands*1.0 /

      real (kind=kind_phys) :: tau_tbl(0:ntbl)  !clr-sky opt dep (for cldy transfer)
      real (kind=kind_phys) :: exp_tbl(0:ntbl)  !transmittance lookup table
      real (kind=kind_phys) :: tfn_tbl(0:ntbl)  !tau transition function; i.e. the
                                                !transition of planck func from mean lyr
                                                !temp to lyr boundary temp as a func of
                                                !opt dep. "linear in tau" method is used.

!  ---  the following variables are used for sub-column cloud scheme

      integer, parameter :: ipsdlw0 = ngptlw     ! initial permutation seed

!  ---  public accessable subprograms

      public lwrad_gpu, rlwinit_gpu, copyin_radlw_main_gpu, taumol, rtrnmr, &
         setcoef, cldprop


! ================
      contains
! ================
      subroutine copyin_radlw_main_gpu(async_id)
      ! must be called after executing rlwinit_gpu
      implicit none
      
      integer, intent(in) :: async_id

      !$acc enter data copyin(nspa, nspb, a0, a1, a2, semiss0, tau_tbl, exp_tbl, &
      !$acc&      tfn_tbl, ngb, wvnlw1, wvnlw2, delwave) async(async_id)

      return
      end subroutine
! --------------------------------
      subroutine lwrad_gpu                                                  &
! --------------------------------
      !  ---  inputs for setaer_lw_gpu:
      &     ( prslk1, tvly, rhly, slmsk, tracer1, xlon, xlat, lsswr, lslwr, me, ntrac, my_max,        &
!  ---  inputs:
     &       plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr_co2, gasvmr_other,                      &
     &       clouds, cldfrc0, clwp, relw, ciwp, reiw, cda1, cda2, cda3, cda4,icseed,sfemis,sfgtmp,                      &
     &       myim, map_jj, map_i, nlay, nlp1, lprnt, myrank, ix, kd,                       &
             nf_vgas, nf_clds, nf_aelw, max_nxjp_acc_length, async_id, &
             fulljj, blocks, smalljj, lhtrlwb, &
!  ---  outputs:
     &       hlwc,topflx_upfxc, topflx_upfx0, sfcflx_upfxc, &
             sfcflx_upfx0, sfcflx_dnfxc, sfcflx_dnfx0,                                          &
!! ---  optional:
     &       hlw0,hlwb, flxprf_upfxc, flxprf_dnfxc, flxprf_upfx0, flxprf_dnfx0    &
     &     )

!  ====================  defination of variables  ====================  !
!                                                                       !
!  input variables:                                                     !
!     plyr (npts,nlay) : layer mean pressures (mb)                      !
!     plvl (npts,nlp1) : interface pressures (mb)                       !
!     tlyr (npts,nlay) : layer mean temperature (k)                     !
!     tlvl (npts,nlp1) : interface temperatures (k)                     !
!     qlyr (npts,nlay) : layer specific humidity (gm/gm)   *see inside  !
!     olyr (npts,nlay) : layer ozone concentration (gm/gm) *see inside  !
!     gasvmr(npts,nlay,:): atmospheric gases amount:                    !
!                       (check module_radiation_gases for definition)   !
!       gasvmr(:,:,1)  -   co2 volume mixing ratio                      !
!       gasvmr(:,:,2)  -   n2o volume mixing ratio                      !
!       gasvmr(:,:,3)  -   ch4 volume mixing ratio                      !
!       gasvmr(:,:,4)  -   o2  volume mixing ratio                      !
!       gasvmr(:,:,5)  -   co  volume mixing ratio                      !
!       gasvmr(:,:,6)  -   cfc11 volume mixing ratio                    !
!       gasvmr(:,:,7)  -   cfc12 volume mixing ratio                    !
!       gasvmr(:,:,8)  -   cfc22 volume mixing ratio                    !
!       gasvmr(:,:,9)  -   ccl4  volume mixing ratio                    !
!     clouds(npts,nlay,:): layer cloud profiles:                        !
!                       (check module_radiation_clouds for definition)  !
!                ---  for  ilwcliq > 0  ---                             !
!       clouds(:,:,1)  -   layer total cloud fraction                   !
!       clouds(:,:,2)  -   layer in-cloud liq water path   (g/m**2)     !
!       clouds(:,:,3)  -   mean eff radius for liq cloud   (micron)     !
!       clouds(:,:,4)  -   layer in-cloud ice water path   (g/m**2)     !
!       clouds(:,:,5)  -   mean eff radius for ice cloud   (micron)     !
!       clouds(:,:,6)  -   layer rain drop water path      (g/m**2)     !
!       clouds(:,:,7)  -   mean eff radius for rain drop   (micron)     !
!       clouds(:,:,8)  -   layer snow flake water path     (g/m**2)     !
!       clouds(:,:,9)  -   mean eff radius for snow flake  (micron)     !
!                ---  for  ilwcliq = 0  ---                             !
!       clouds(:,:,1)  -   layer total cloud fraction                   !
!       clouds(:,:,2)  -   layer cloud optical depth                    !
!       clouds(:,:,3)  -   layer cloud single scattering albedo         !
!       clouds(:,:,4)  -   layer cloud asymmetry factor                 !
!     icseed(npts)   : auxiliary special cloud related array            !
!                      when module variable isubclw=2, it provides      !
!                      permutation seed for each column profile that    !
!                      are used for generating random numbers.          !
!                      when isubclw /=2, it will not be used.           !
!     aerosols(npts,nlay,nbands,:) : aerosol optical properties         !
!                       (check module_radiation_aerosols for definition)!
!        (:,:,:,1)     - optical depth                                  !
!        (:,:,:,2)     - single scattering albedo                       !
!        (:,:,:,3)     - asymmetry parameter                            !
!     sfemis (npts)  : surface emissivity                               !
!     sfgtmp (npts)  : surface ground temperature (k)                   !
!     npts           : total number of horizontal points                !
!     nlay, nlp1     : total number of vertical layers, levels          !
!     lprnt          : cntl flag for diagnostic print out               !
!                                                                       !
!  output variables:                                                    !
!     hlwc  (npts,nlay): total sky heating rate (k/day or k/sec)        !
!     topflx(npts)     : radiation fluxes at top, component:            !
!                        (check module_radlw_paramters for definition)  !
!        upfxc           - total sky upward flux at top (w/m2)          !
!        upfx0           - clear sky upward flux at top (w/m2)          !
!     sfcflx(npts)     : radiation fluxes at sfc, component:            !
!                        (check module_radlw_paramters for definition)  !
!        upfxc           - total sky upward flux at sfc (w/m2)          !
!        upfx0           - clear sky upward flux at sfc (w/m2)          !
!        dnfxc           - total sky downward flux at sfc (w/m2)        !
!        dnfx0           - clear sky downward flux at sfc (w/m2)        !
!                                                                       !
!! optional output variables:                                           !
!     hlwb(npts,nlay,nbands): spectral band total sky heating rates     !
!     hlw0  (npts,nlay): clear sky heating rate (k/day or k/sec)        !
!     flxprf(npts,nlp1): level radiative fluxes (w/m2), components:     !
!                        (check module_radlw_paramters for definition)  !
!        upfxc           - total sky upward flux                        !
!        dnfxc           - total sky dnward flux                        !
!        upfx0           - clear sky upward flux                        !
!        dnfx0           - clear sky dnward flux                        !
!                                                                       !
!  external module variables:  (in physpara)                            !
!   ilwrgas - control flag for rare gases (ch4,n2o,o2,cfcs, etc.)       !
!           =0: do not include rare gases                               !
!           >0: include all rare gases                                  !
!   ilwcliq - control flag for liq-cloud optical properties             !
!           =0: input cloud optical depth, ignor ilwcice                !
!           =1: input cld liqp & reliq, hu & stamnes (1993)             !
!           =2: not used                                                !
!   ilwcice - control flag for ice-cloud optical properties             !
!           *** if ilwcliq==0, ilwcice is ignored                       !
!           =1: input cld icep & reice, ebert & curry (1997)            !
!           =2: input cld icep & reice, streamer (1996)                 !
!           =3: input cld icep & reice, fu (1998)                       !
!   isubclw - sub-column cloud approximation control flag               !
!           =0: no sub-col cld treatment, use grid-mean cld quantities  !
!           =1: mcica sub-col, prescribed seeds to get random numbers   !
!           =2: mcica sub-col, providing array icseed for random numbers!
!   iovrlw  - cloud overlapping control flag                            !
!           =0: random overlapping clouds                               !
!           =1: maximum/random overlapping clouds                       !
!           =2: maximum overlap cloud (used for isubclw>0 only)         !
!   ivflip  - control flag for vertical index direction                 !
!           =0: vertical index from toa to surface                      !
!           =1: vertical index from surface to toa                      !
!                                                                       !
!  module parameters, control variables:                                !
!     nbands           - number of longwave spectral bands              !
!     maxgas           - maximum number of absorbing gaseous            !
!     maxxsec          - maximum number of cross-sections               !
!     ngptlw           - total number of g-point subintervals           !
!     ng##             - number of g-points in band (##=1-16)           !
!     ngb(ngptlw)      - band indices for each g-point                  !
!     bpade            - pade approximation constant (1/0.278)          !
!     nspa,nspb(nbands)- number of lower/upper ref atm's per band       !
!     delwave(nbands)  - longwave band width (wavenumbers)              !
!     ipsdlw0          - permutation seed for mcica sub-col clds        !
!                                                                       !
!  major local variables:                                               !
!     pavel  (nlay)         - layer pressures (mb)                      !
!     delp   (nlay)         - layer pressure thickness (mb)             !
!     tavel  (nlay)         - layer temperatures (k)                    !
!     tz     (0:nlay)       - level (interface) temperatures (k)        !
!     semiss (nbands)       - surface emissivity for each band          !
!     wx     (nlay,maxxsec) - cross-section molecules concentration     !
!     coldry (nlay)         - dry air column amount                     !
!                                   (1.e-20*molecules/cm**2)            !
!     cldfrc (0:nlp1)       - layer cloud fraction                      !
!     taucld (nbands,nlay)  - layer cloud optical depth for each band   !
!     cldfmc (ngptlw,nlay)  - layer cloud fraction for each g-point     !
!     tauaer (nbands,nlay)  - aerosol optical depths                    !
!     fracs  (ngptlw,nlay)  - planck fractions                          !
!     tautot (ngptlw,nlay)  - total optical depths (gaseous+aerosols)   !
!     colamt (nlay,maxgas)  - column amounts of absorbing gases         !
!                             1-maxgas are for watervapor, carbon       !
!                             dioxide, ozone, nitrous oxide, methane,   !
!                             oxigen, carbon monoxide, respectively     !
!                             (molecules/cm**2)                         !
!     pwvcm                 - column precipitable water vapor (cm)      !
!     secdiff(nbands)       - variable diffusivity angle defined as     !
!                             an exponential function of the column     !
!                             water amount in bands 2-3 and 5-9.        !
!                             this reduces the bias of several w/m2 in  !
!                             downward surface flux in high water       !
!                             profiles caused by using the constant     !
!                             diffusivity angle of 1.66.         (mji)  !
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
!     laytrop               - tropopause layer index at which switch is !
!                             made from one conbination kew species to  !
!                             another.                                  !
!     jp(nlay),jt(nlay),jt1(nlay)                                       !
!                           - lookup table indexes                      !
!     totuflux(0:nlay)      - total-sky upward longwave flux (w/m2)     !
!     totdflux(0:nlay)      - total-sky downward longwave flux (w/m2)   !
!     htr(nlay)             - total-sky heating rate (k/day or k/sec)   !
!     totuclfl(0:nlay)      - clear-sky upward longwave flux (w/m2)     !
!     totdclfl(0:nlay)      - clear-sky downward longwave flux (w/m2)   !
!     htrcl(nlay)           - clear-sky heating rate (k/day or k/sec)   !
!     fnet    (0:nlay)      - net longwave flux (w/m2)                  !
!     fnetc   (0:nlay)      - clear-sky net longwave flux (w/m2)        !
!                                                                       !
!                                                                       !
!  ======================    end of definitions    ===================  !

!  ---  inputs for setaer_lw_gpu
      integer, intent(in) :: me, ntrac, my_max
      logical, intent(in) :: lsswr, lslwr
      real (kind=kind_phys), dimension(ix, my_max),  intent(in) ::  slmsk,  &
         xlon, xlat
      real (kind=kind_phys), dimension(nxptot,nlay), intent(in)  :: rhly, prslk1, tvly
      real (kind=kind_phys), dimension(nxptot,nlay,ntrac), intent(in)   :: tracer1


!  ---  inputs:
      integer, intent(in) :: myim(fulljj), nlay, nlp1, myrank, fulljj, &
         nf_vgas, nf_clds, nf_aelw, ix, blocks, smalljj, kd, max_nxjp_acc_length
      integer, dimension(nxptot), intent(in) :: map_jj, map_i
      integer, intent(in) :: icseed(ix, fulljj)

      logical,  intent(in) :: lprnt, lhtrlwb

      real (kind=kind_phys), dimension(nxptot, nlp1), intent(in) :: plvl,  &
     &       tlvl
      real (kind=kind_phys), dimension(nxptot, nlay), intent(in) :: plyr,  &
     &       tlyr, qlyr, olyr

      real (kind=kind_phys), dimension(nxptot,nlay),intent(in):: gasvmr_co2
      real (kind=kind_phys), dimension(nf_vgas),intent(in):: gasvmr_other
      real (kind=kind_phys), dimension(nxptot,nlay,nf_clds),intent(in):: clouds
      real (kind=kind_phys), dimension(nxptot,nlay),intent(in):: cldfrc0
      real (kind=kind_phys), dimension(nxptot,nlay),intent(inout):: cda1
      real (kind=kind_phys), dimension(nxptot,nlay),intent(in):: &
         clwp, relw, ciwp, reiw, cda2, cda3, cda4

      real (kind=kind_phys), dimension(nxptot), intent(in) :: sfemis,     &
     &       sfgtmp

     ! real (kind=kind_phys), dimension(ix,nlay,nbands,nf_aelw, fulljj),intent(in):: &
     !&       aerosols

!  ---  outputs:
      real (kind=kind_phys), dimension(ix, nlay, fulljj), intent(out) :: hlwc

      real (kind=kind_phys),    dimension(ix, fulljj), intent(out) :: topflx_upfxc, topflx_upfx0
      real (kind=kind_phys),    dimension(ix, fulljj), intent(out) :: sfcflx_upfxc, sfcflx_upfx0, &
         sfcflx_dnfxc, sfcflx_dnfx0

!! ---  optional outputs:
      real (kind=kind_phys), dimension(ix,nlay,nbands, fulljj),optional,      &
     &       intent(out) :: hlwb
      real (kind=kind_phys), dimension(ix, nlay, fulljj),       optional,      &
     &       intent(out) :: hlw0
      real (kind=kind_phys),    dimension(ix, nlp1, fulljj),       optional,      &
     &       intent(out) :: flxprf_upfxc, flxprf_dnfxc, flxprf_upfx0, flxprf_dnfx0 

!  ---  locals:
      real (kind=kind_phys), dimension(max_nxjp_acc_length, 0:nlay) :: totuflux, totdflux !, tz
      real (kind=kind_phys), dimension(max_nxjp_acc_length, 0:nlay) :: totuclfl, totdclfl

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay) :: htr, htrcl
      real (kind=kind_phys), dimension(max_nxjp_acc_length, 0:nlp1) :: cldfrc
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay) :: tavel, delp,   &
     &       h2ovmr, o3vmr
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay)   :: pavel, &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, &
     &       selffac, selffrac, forfac, forfrac, minorfrac, scaleminor, &
     &       scaleminorn2

      real (kind=kind_phys), dimension(max_nxjp_acc_length, 0:nlay, nbands) :: pklev, pklay

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay, nbands) :: htrb
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay, nbands) :: taucld

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nbands) :: semiss, secdiff

!  ---  column amount of absorbing gases:
!       (:,m) m = 1-h2o, 2-co2, 3-o3, 4-n2o, 5-ch4, 6-o2, 7-co
      real (kind=kind_phys) :: colamtr

!  ---  column cfc cross-section amounts:
!       (:,m) m = 1-ccl4, 2-cfc11, 3-cfc12, 4-cfc22

!  ---  reference ratios of binary species parameter in lower atmosphere:
!       (:,m,:) m = 1-h2o/co2, 2-h2o/o3, 3-h2o/n2o, 4-h2o/ch4, 5-n2o/co2, 6-o3/co2

      real (kind=kind_phys) :: tem0, tem1, tem2, summol
      real (kind=kind_phys), dimension(max_nxjp_acc_length) :: pwvcm, stemp
      integer, dimension(max_nxjp_acc_length) :: ipseed, laytrop
      logical, dimension(max_nxjp_acc_length) :: lcf1

      integer, dimension(max_nxjp_acc_length, nlay) :: jp, jt, jt1, indself, indfor, indminor
      integer                  :: iplon, i, j, k, k1, jj
      !!!!!!!!!!!!!!!!!!!!!!!!!
      real (kind=kind_phys), dimension(nlay)   :: htr_im, htrcl_im

      real (kind=kind_phys), dimension(nlay)   :: pavel_im, tavel_im, delp_im,   &
     &       clwp_im, ciwp_im, relw_im, reiw_im, cda1_im, cda2_im, cda3_im, cda4_im,            &
     &       coldry_im, colbrd_im, h2ovmr_im, o3vmr_im, fac00_im, fac01_im, fac10_im, fac11_im, &
     &       selffac_im, selffrac_im, forfac_im, forfrac_im, minorfrac_im, scaleminor_im, &
     &       scaleminorn2_im, temcol_im
      real (kind=kind_phys), dimension(0:nlp1) :: cldfrc_im

      real (kind=kind_phys), dimension(0:nlay) :: totuflux_im, totdflux_im,   &
     &       totuclfl_im, totdclfl_im, tz_im
      real (kind=kind_phys), dimension(0:nlay, nbands) :: pklev_im, pklay_im

      real (kind=kind_phys), dimension(nlay,nbands) :: htrb_im
      real (kind=kind_phys), dimension(nlay,nbands) :: taucld_im, tauaer_im
      real (kind=kind_phys), dimension(nlay,ngptlw) :: fracs_im, tautot_im,   &
     &       cldfmc_im
      real (kind=kind_phys), dimension(nbands) :: semiss_im, secdiff_im
      real (kind=kind_phys) :: colamt_im(nlay,maxgas), wx_im(nlay,maxxsec), rfrate_im(nlay,nrates,2)
      integer, dimension(nlay) :: jp_im, jt_im, jt1_im, indself_im, indfor_im, indminor_im
      character(len=4) :: myrank_str
      integer :: async_id, ig, ib, n, m
      real (kind=kind_phys), allocatable, dimension(:,:,:)   :: cldf
      logical,allocatable,dimension(:,:,:,:) :: lcloudy
      real (kind=kind_phys), allocatable, dimension(:)   :: cldf_im
      logical,allocatable,dimension(:,:) :: lcloudy_im

      real (kind=kind_phys), dimension(:,:,:,:), allocatable :: cldfmc
      integer :: jb, jf, jjoffset, jbs, jbe, blocks_local, blockjj, jbs_nxjp_acc, &
         jbe_nxjp_acc, nxjp_acc_length, rr
      ! for GPU register
      real (kind=kind_phys) :: temcol

      real (kind=kind_phys), dimension(:,:,:), allocatable:: tauaer_input
      real (kind=kind_phys), dimension(:,:), allocatable:: tlvl_input, clwp_input, &
         relw_input, ciwp_input, reiw_input, cda1_input, cda2_input, cda3_input, cda4_input
      integer, dimension(nbands) :: ng_array, ns_array
      integer :: ng00
      data ng_array(:) /ng01, ng02, ng03, ng04, ng05, ng06, ng07, ng08, &
                        ng09, ng10, ng11, ng12, ng13, ng14, ng15, ng16/ 
      data ns_array(:) /ns01, ns02, ns03, ns04, ns05, ns06, ns07, ns08, &
                        ns09, ns10, ns11, ns12, ns13, ns14, ns15, ns16/ 
      ng00 = maxval(ng_array)
   
      
      
      ! GPU: variable name changed: CPU - topflx%upfxc, GPU - topflx_upfxc
      ! GPU: variable name changed: CPU - topflx%upfx0, GPU - topflx_upfx0
      ! GPU: variable name changed: CPU - sfcflx%upfxc, GPU - sfcflx_upfxc
      ! GPU: variable name changed: CPU - sfcflx%upfx0, GPU - sfcflx_upfx0
      ! GPU: variable name changed: CPU - sfcflx%dnfxc, GPU - sfcflx_dnfxc
      ! GPU: variable name changed: CPU - sfcflx%dnfx0, GPU - sfcflx_dnfx0
      ! GPU: variable name changed: CPU - flxprf%upfxc, GPU - flxprf_upfxc
      ! GPU: variable name changed: CPU - flxprf%dnfxc, GPU - flxprf_dnfxc
      ! GPU: variable name changed: CPU - flxprf%upfx0, GPU - flxprf_upfx0
      ! GPU: variable name changed: CPU - flxprf%dnfx0, GPU - flxprf_dnfx0
      ! GPU: variable name changed: CPU - rfrate,       GPU - <dismiss>
      ! GPU: variable name changed: CPU - htrb,         GPU - hlwb
      ! GPU: variable name changed: CPU - wx,           GPU - <dismiss>
      ! GPU: variable name changed: CPU - colamt,       GPU - colamtr
      ! GPU: variable name changed: CPU - htr,          GPU - hlwc
      ! GPU: variable name changed: CPU - htrcl,        GPU - hlw0
      
      ! GPU: GPU version dismiss the array 'aerosols'. Originally, 'tauaer' is 
      ! GPU: computed by tauaer = aerosols(1) * (1. - aerosols(2)) in this 
      ! GPU: subroutine (in the current CPU version). In this GPU version, 'tauaer'
      ! GPU: is computed at module_radiation_aerosols/setaer_gpu, and then passed
      ! GPU: into this subroutine.
      

      if (isubclw > 0) allocate(cldfmc(ix, nlay,ngptlw, fulljj))
      
      !
      !===> ... begin here
      !

      !  --- ...  initialization

      lhlwb  = present ( hlwb )
      lhlw0  = present ( hlw0 )
      lflxprf_upfxc= present ( flxprf_upfxc )
      lflxprf_dnfxc= present ( flxprf_dnfxc )
      lflxprf_upfx0= present ( flxprf_upfx0 )
      lflxprf_dnfx0= present ( flxprf_dnfx0 )
      lflxprf= ((lflxprf_upfxc .and. lflxprf_dnfxc) .and. (lflxprf_upfx0 .and. lflxprf_dnfx0))

      !      if (myrank .eq. 0) print *,'$$$$ lwrad start $$$$$'
      !      if (myrank .eq. 0) print *,'$$$$ lhlwb   = ', lhlwb
      !      if (myrank .eq. 0) print *,'$$$$ lhlw0   = ', lhlw0
      !      if (myrank .eq. 0) print *,'$$$$ lflxprf = ', lflxprf
      
      !$acc data create(lcf1) async(async_id)
      do jb = 1, blocks
         jjoffset = (fulljj-1)*(jb-1)/blocks
         jbs = jjoffset+1
         jbe = (fulljj-1)*jb/blocks
         blockjj = jbe - jbs + 1
         jbs_nxjp_acc = nxjp_acc(jbs)
         jbe_nxjp_acc = nxjp_acc(jbe+1)
         nxjp_acc_length = jbe_nxjp_acc - jbs_nxjp_acc

         !GPU: Be careful when calling rtrn and rtrnmc !!
         !GPU: Before calling rtrn and rtrnmc, taumol and setaer_lw_gpu must be 
         !GPU: called. Variabel tauaer, which is not on GPU in this situation, 
         !GPU: must be created/deleted at suitable location (be created before 
         !GPU: setaer_lw and be deleted after taumol).
         !!$acc enter data create(tauaer) async(async_id)
         !call nvtxStartRange("lw_setaer")
         !call setaer_lw_gpu                                                       &
         !!  ---  inputs:
         !&     ( plvl,plyr,prslk1,tvly,rhly,slmsk,tracer1, xlon,xlat,        &
         !&       myim,nlay,nlp1,lsswr,lslwr,me,myrank, ix, map_jj, map_i, nxptot, async_id,                          &
         !!  ---  outputs:
         !&       tauaer                                              &
         !!    &       faersw,faerlw,aerodp                                       &
         !&     )
         !!$acc enter data create(tauaer) async(async_id)
         !call nvtxEndRange

         !, jj  --- ...  change random number seed value for each radiation invocation
         !$acc enter data create(ipseed) async(async_id)
         if     ( isubclw == 1 ) then     ! advance prescribed permutation seed
            !$acc parallel loop private(i, m) async(async_id)
            do n = jbs_nxjp_acc, jbs_nxjp_acc + nxjp_acc_length - 1
               i = map_i(n)
               m = n - jbs_nxjp_acc + 1
               ipseed(m) = ipsdlw0 + i
            end do
         elseif ( isubclw == 2 ) then     ! use input array of permutaion seeds
            !$acc parallel loop private(jf, iplon, m) async(async_id)
            do n = jbs_nxjp_acc, jbs_nxjp_acc + nxjp_acc_length - 1
               jf = map_jj(n)
               i = map_i(n)
               m = n - jbs_nxjp_acc + 1
               ipseed(m) = icseed(i, jf)
            end do
         endif

         !     if ( lprnt ) then
         !       print *,'  in radlw, isubclw, ipsdlw0,ipseed =',                &
         !    &          isubclw, ipsdlw0, ipseed
         !     endif

         !  --- ...  loop over horizontal npts profiles
         !$acc enter data create(semiss, secdiff) async(async_id)

         !$acc parallel loop collapse(2) private(iplon) async(async_id)
         do j = 1, nbands
            do n = jbs_nxjp_acc, jbs_nxjp_acc + nxjp_acc_length - 1
               m = n - jbs_nxjp_acc + 1
               if (sfemis(n) > eps .and. sfemis(n) <= 1.0) then  ! input surface emissivity
                  semiss(m, j) = sfemis(n)
               else                                                      ! use default values
                  semiss(m, j) = semiss0(j)
               endif
            enddo
         end do

               !stemp = sfgtmp(iplon, jj)          ! surface ground temp

               !  --- ...  prepare atmospheric profile for use in rrtm
               !           the vertical index of internal array is from surface to top

               !  --- ...  molecular amounts are input or converted to volume mixing ratio
               !           and later then converted to molecular amount (molec/cm2) by the
               !           dry air column coldry (in molec/cm2) which is calculated from the
               !           layer pressure thickness (in mb), based on the hydrostatic equation
               !  --- ...  and includes a correction to account for h2o in the layer.
         !$acc enter data create(pavel, coldry, pwvcm, cldfrc, tavel, delp, &
         !$acc&      h2ovmr, o3vmr) async(async_id)




         if (ivflip == 0) then
            tem1 = 100.0 * con_g
            tem2 = 1.0e-20 * 1.0e3 * con_avgd
            !$acc parallel loop collapse(2) private(tem0, temcol, m) async(async_id)
            do k = 1, nlay
               do n = jbs_nxjp_acc, jbs_nxjp_acc + nxjp_acc_length - 1
                  m = n - jbs_nxjp_acc + 1
                  k1 = nlp1 - k
                  pavel(m, k)= plyr(n,k1)
                  delp(m, k) = plvl(n,k1+1) - plvl(n,k1)
                  tavel(m, k)= tlyr(n,k1)

                  !  --- ...  set absorber amount
                  !test use
                  !           h2ovmr(k)= max(f_zero,qlyr(iplon,k1)*amdw)                  ! input mass mixing ratio
                  !           h2ovmr(k)= max(f_zero,qlyr(iplon,k1))                       ! input vol mixing ratio
                  !           o3vmr (k)= max(f_zero,olyr(iplon,k1))                       ! input vol mixing ratio
                  !ncep model use
                  h2ovmr(m, k)= max(f_zero,qlyr(n,k1)                        &
                  &                           *amdw/(f_one-qlyr(n,k1)))          ! input specific humidity
                  o3vmr (m, k)= max(f_zero,olyr(n,k1)*amdo3)                 ! input mass mixing ratio

                  !  --- ...  tem0 is the molecular weight of moist air
                  tem0 = (f_one - h2ovmr(m, k))*con_amd + h2ovmr(m, k)*con_amw
                  coldry(m, k) = tem2*delp(m, k) / (tem1*tem0*(f_one+h2ovmr(m, k)))
                  temcol = 1.0e-12 * coldry(m, k)


               !  --- ...  set up col amount for rare gases, convert from volume mixing ratio
               !           to molec/cm2 based on coldry (scaled to 1.0e-20)
               end do
            end do

            !$acc parallel loop private(m) async(async_id)
            do n = jbs_nxjp_acc, jbs_nxjp_acc + nxjp_acc_length - 1
               m = n - jbs_nxjp_acc + 1
               cldfrc(m, 0)    = f_one       ! padding value only
               cldfrc(m, nlp1) = f_zero      ! padding value only
            end do
            !$acc parallel loop collapse(2) private(tem0, m) async(async_id)
            do k = 1, nlay
               do n = jbs_nxjp_acc, jbs_nxjp_acc + nxjp_acc_length - 1
                  m = n - jbs_nxjp_acc + 1
                  cldfrc(m, k) = cldfrc0(n, k)
               end do
            end do

            !  --- ...  compute precipitable water vapor for diffusivity angle adjustments
            !$acc parallel loop private(tem1, tem2, tem0, colamtr, m) async(async_id)
            do n = jbs_nxjp_acc, jbs_nxjp_acc + nxjp_acc_length - 1
               m = n - jbs_nxjp_acc + 1
               tem1 = f_zero
               tem2 = f_zero
               !$acc loop seq
               do k = 1, nlay
                  colamtr = max(f_zero,    coldry(m, k)*h2ovmr(m, k))          ! h2o
                  tem1 = tem1 + coldry(m, k) + colamtr
                  tem2 = tem2 + colamtr
               enddo

               tem0 = 10.0 * tem2 / (amdw * tem1 * con_g)
               pwvcm(m) = tem0 * plvl(n,nlp1)
            end do
            if (ilwcliq <= 0) then    ! use prognostic cloud method
               !$acc parallel loop collapse(2) private(k1, m) async(async_id)
               do k = 1, nlay
                  do n = jbs_nxjp_acc, jbs_nxjp_acc + nxjp_acc_length - 1
                     m = n - jbs_nxjp_acc + 1
                     k1 = nlp1 - k
                     cldfrc(m, k)  = clouds(n,k1,1)
                     cda1(n, k)  = clouds(n,k1,2)
                  end do
               end do
            endif                      ! end if_ilwcliq
         end if

         if (ivflip .ne. 0) then                        ! input from sfc to toa
            tem1 = 100.0 * con_g
            tem2 = 1.0e-20 * 1.0e3 * con_avgd
            !$acc parallel loop collapse(2) private(tem0, temcol, m) async(async_id)
            do k = 1, nlay
               do n = jbs_nxjp_acc, jbs_nxjp_acc + nxjp_acc_length - 1
                  m = n - jbs_nxjp_acc + 1
                  pavel(m, k)= plyr(n,k)
                  delp(m, k) = plvl(n,k) - plvl(n,k+1)
                  tavel(m, k)= tlyr(n,k)

                  !  --- ...  set absorber amount
                  !test use
                  !           h2ovmr(k)= max(f_zero,qlyr(iplon,k)*amdw)                   ! input mass mixing ratio
                  !           h2ovmr(k)= max(f_zero,qlyr(iplon,k))                        ! input vol mixing ratio
                  !           o3vmr (k)= max(f_zero,olyr(iplon,k))                        ! input vol mixing ratio
                  !ncep model use
                  h2ovmr(m, k)= max(f_zero,qlyr(n,k)                         &
                  &                           *amdw/(f_one-qlyr(n,k)))           ! input specific humidity
                  o3vmr (m, k)= max(f_zero,olyr(n,k)*amdo3)                  ! input mass mixing ratio

                  !  --- ...  tem0 is the molecular weight of moist air
                  tem0 = (f_one - h2ovmr(m, k))*con_amd + h2ovmr(m, k)*con_amw
                  coldry(m, k) = tem2*delp(m, k) / (tem1*tem0*(f_one+h2ovmr(m, k)))
                  temcol = 1.0e-12 * coldry(m, k)


               !  --- ...  set up col amount for rare gases, convert from volume mixing ratio
               !           to molec/cm2 based on coldry (scaled to 1.0e-20)
               end do
            enddo



                  !  --- ...  set aerosol optical properties
         
            !$acc parallel loop collapse(2) private(tem0, m) async(async_id)
            do k = 1, nlay
               do n = jbs_nxjp_acc, jbs_nxjp_acc + nxjp_acc_length - 1
                  m = n - jbs_nxjp_acc + 1
                  cldfrc(m, k) = cldfrc0(n, k)
               end do
            end do

            !  --- ...  compute precipitable water vapor for diffusivity angle adjustments
            !$acc parallel loop private(tem1, tem2, tem0, colamtr, m) async(async_id)
            do n = jbs_nxjp_acc, jbs_nxjp_acc + nxjp_acc_length - 1
               m = n - jbs_nxjp_acc + 1
               cldfrc(m, 0)    = f_one       ! padding value only
               cldfrc(m, nlp1) = f_zero      ! padding value only
               tem1 = f_zero
               tem2 = f_zero
               !$acc loop seq
               do k = 1, nlay
                  colamtr = max(f_zero,    coldry(m, k)*h2ovmr(m, k))          ! h2o
                  tem1 = tem1 + coldry(m, k) + colamtr
                  tem2 = tem2 + colamtr
               enddo

               tem0 = 10.0 * tem2 / (amdw * tem1 * con_g)
               pwvcm(m) = tem0 * plvl(n,1)
            end do
            if (ilwcliq <= 0) then    ! use prognostic cloud method
               !$acc parallel loop collapse(2) private(m) async(async_id)
               do k = 1, nlay
                  do n = jbs_nxjp_acc, jbs_nxjp_acc + nxjp_acc_length - 1
                     m = n - jbs_nxjp_acc + 1
                     cldfrc(m, k)  = clouds(n,k,1)
                     cda1(n, k)  = clouds(n,k,2)
                  end do
               end do
            endif                      ! end if_ilwcliq
         endif                       ! if_ivflip

               !  --- ...  compute column amount for broadening gases
         !$acc enter data create(colbrd) async(async_id)
         !$acc parallel loop collapse(2) private(k1, m) async(async_id)
         do k = 1, nlay
            do n = jbs_nxjp_acc, jbs_nxjp_acc + nxjp_acc_length - 1
               m = n - jbs_nxjp_acc + 1
               if (ivflip == 0) then
                  k1 = nlay + 1 - k
               else
                  k1 = k
               end if
               temcol = 1.0e-12 * coldry(m, k)
               summol = f_zero
               summol = summol + max(temcol, coldry(m, k)*gasvmr_co2(n,k1)) ! co2
               summol = summol + max(temcol, coldry(m, k)*o3vmr(m, k))           ! o3
               if (ilwrgas > 0) then
                  summol = summol + max(temcol, coldry(m, k)*gasvmr_other(2))  ! n2o
                  summol = summol + max(temcol, coldry(m, k)*gasvmr_other(3))  ! ch4
                  summol = summol + max(f_zero,    coldry(m, k)*gasvmr_other(4))  ! o2
                  summol = summol + max(f_zero,    coldry(m, k)*gasvmr_other(5))  ! co
               else
                  summol = summol + f_zero     ! n2o
                  summol = summol + f_zero     ! ch4
                  summol = summol + f_zero     ! o2
                  summol = summol + f_zero     ! co
               endif
               colbrd(m, k) = coldry(m, k) - summol
            enddo
         end do

         !  --- ...  compute diffusivity angle adjustments

         tem1 = 1.80
         tem2 = 1.50
         !$acc parallel loop collapse(2) private(m) async(async_id)
         do j = 1, nbands
            do n = jbs_nxjp_acc, jbs_nxjp_acc + nxjp_acc_length - 1
               m = n - jbs_nxjp_acc + 1
               if (j==1 .or. j==4 .or. j==10) then
                  secdiff(m, j) = 1.66
               else
                  secdiff(m, j) = min( tem1, max( tem2,                          &
                  &                   a0(j)+a1(j)*exp(a2(j)*pwvcm(m)) ))
               endif
            end do
         end do
         !$acc exit data delete(pwvcm) async(async_id)
         

         !     if (lprnt) then
         !      print *,'  coldry',coldry
         !      print *,' wx(*,1) ',(wx(k,1),k=1,nlay)
         !      print *,' wx(*,2) ',(wx(k,2),k=1,nlay)
         !      print *,' wx(*,3) ',(wx(k,3),k=1,nlay)
         !      print *,' wx(*,4) ',(wx(k,4),k=1,nlay)
         !      print *,' iplon ',iplon
         !      print *,'  pavel ',pavel
         !      print *,'  delp ',delp
         !      print *,'  tavel ',tavel
         !      print *,'  tz ',tz
         !      print *,' h2ovmr ',h2ovmr
         !      print *,' o3vmr ',o3vmr
         !     endif

         !  --- ...  for cloudy atmosphere, use cldprop to set cloud optical properties
         !$acc enter data create(taucld) async(async_id)
         !$acc parallel loop private(m) async(async_id)
         do n = jbs_nxjp_acc, jbs_nxjp_acc + nxjp_acc_length - 1
            m = n - jbs_nxjp_acc + 1
            lcf1(m) = .false.
            !$acc loop seq
            do k = 1, nlay ! lab_do_k0
               if ( cldfrc(m, k) > eps ) then
                  lcf1(m) = .true.
                  exit ! lab_do_k0
               endif
            enddo  ! lab_do_k0
         end do
         !call nvtxStartRange("lw_cldprop")
         call cldprop                                                  &
            !  ---  inputs:
            &     ( cldfrc,clwp,relw,ciwp, &
                  reiw,cda1,cda2, &
                  cda3,cda4,            &
            &       nlay, nlp1, ipseed, ix, myim(jbs:jbe), lcf1, map_jj(jbs_nxjp_acc:jbe_nxjp_acc-1), &
                  map_i(jbs_nxjp_acc:jbe_nxjp_acc-1), nxjp_acc_length, jjoffset, &
                  max_nxjp_acc_length, jbs_nxjp_acc, async_id, smalljj, blockjj,                                  &
            !  ---  outputs:
            &       taucld                                             &
            &     )
         if ( isubclw > 0 ) then      ! mcica sub-col clouds approx
            allocate(lcloudy(ix, ngptlw,nlay, fulljj))
            allocate(cldf_im(nlay))
            allocate(lcloudy_im(ngptlw,nlay))
            allocate(cldf(ix, nlay, fulljj))
            !$acc parallel loop collapse(3) private(jj, m) async(async_id)
            do ig = 1, ngptlw
               do k = 1, nlay
                  do n = jbs_nxjp_acc, jbs_nxjp_acc + nxjp_acc_length - 1
                     jj = map_jj(n) - jjoffset
                     m = n - jbs_nxjp_acc + 1
                     cldfmc(iplon, k,ig, jj) = f_zero
                  enddo
               end do
            end do
            do jj = 1, blockjj
               jf = jjoffset+jj
               do iplon = 1, myim(jf) ! lab_do_iplon
                  n = i + nxjp_acc(jj) - 1
                  if ( lcf1(m) ) then
                  !  ---  distribute cloud properties to each g-point

                     do k = 1, nlay
                        if ( cldfrc(n, k) < cldmin ) then
                           cldf(iplon, k, jj) = f_zero
                        else
                           cldf(iplon, k, jj) = cldfrc(n, k)
                        endif
                     enddo

                     !  --- ...  call sub-column cloud generator
                     do k = 1, nlay
                        cldf_im(k) = cldf(iplon, k, jj)
                        do ib = 1, ngptlw
                           lcloudy_im(ib, k) = lcloudy(iplon, ib, k, jj)
                        end do
                     end do
                     call mcica_subcol                                               &
                        !  ---  inputs:
                        &     ( cldf_im, nlay, ipseed(m),                                        &
                        !  ---  output:
                        &       lcloudy_im                                                    &
                        &     )
                     do k = 1, nlay
                        cldf(iplon, k, jj) = cldf_im(k)
                        do ib = 1, ngptlw
                           lcloudy(iplon, ib, k, jj) = lcloudy_im(ib, k)
                        end do
                     end do

                     do k = 1, nlay
                        do ig = 1, ngptlw
                           if ( lcloudy(iplon, ig,k, jj) ) then
                              cldfmc(iplon, k,ig, jj) = f_one
                           else
                              cldfmc(iplon, k,ig, jj) = f_zero
                           endif
                        enddo
                     enddo
                  end if
               end do
            end do
            deallocate(lcloudy)
            deallocate(cldf_im)
            deallocate(lcloudy_im)
            deallocate(cldf)
         end if
         !$acc exit data delete(ipseed) async(async_id)
         


         !call nvtxEndRange

               !     if (lprnt) then
               !      print *,' after cldprop'
               !      print *,' clwp',clwp
               !      print *,' ciwp',ciwp
               !      print *,' relw',relw
               !      print *,' reiw',reiw
               !      print *,' taucl',cda1
               !      print *,' cldfrac',cldfrc
               !     endif
         !call nvtxStartRange("lw_setcoef")
         !$acc enter data create(fac00, fac01, fac10, fac11, selffac, selffrac, &
         !$acc&      forfac, forfrac, minorfrac, scaleminor, scaleminorn2, laytrop, &
         !$acc&      jp, jt, jt1, indself, indfor, indminor) async(async_id)


         call setcoef                                                    &
            !  ---  inputs:
            &     ( pavel,tavel,sfgtmp,nf_vgas, h2ovmr, gasvmr_co2, &
                  o3vmr, gasvmr_other,coldry,colbrd,          &
            &       nlay, nlp1, ix, myim(jbs:jbe), map_jj(jbs_nxjp_acc:jbe_nxjp_acc-1), &
                  map_i(jbs_nxjp_acc:jbe_nxjp_acc-1), nxjp_acc_length, jjoffset, &
                  max_nxjp_acc_length, jbs_nxjp_acc, async_id, smalljj, blockjj,      &
            !  ---  outputs:
            &       laytrop,pklay,pklev,jp,jt,jt1,                             &
            &       fac00,fac01,fac10,fac11,                            &
            &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
            &       minorfrac,scaleminor,scaleminorn2,indminor                 &
            &     )
         !call nvtxEndRange

               !     if (lprnt) then
               !      print *,'laytrop',laytrop
               !      print *,'colh2o',(colamt(k,1),k=1,nlay)
               !      print *,'colco2',(colamt(k,2),k=1,nlay)
               !      print *,'colo3', (colamt(k,3),k=1,nlay)
               !      print *,'coln2o',(colamt(k,4),k=1,nlay)
               !      print *,'colch4',(colamt(k,5),k=1,nlay)
               !      print *,'fac00',fac00
               !      print *,'fac01',fac01
               !      print *,'fac10',fac10
               !      print *,'fac11',fac11
               !      print *,'jp',jp
               !      print *,'jt',jt
               !      print *,'jt1',jt1
               !      print *,'selffac',selffac
               !      print *,'selffrac',selffrac
               !      print *,'indself',indself
               !      print *,'forfac',forfac
               !      print *,'forfrac',forfrac
               !      print *,'indfor',indfor
               !     endif

               !  --- ...  calculate the gaseous optical depths and planck fractions for
               !           each longwave spectral band.
         !call taumol                                                     &
         !   !  ---  inputs:
         !   &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2(:,:,jbs:jbe), &
         !           o3vmr, gasvmr_other, colbrd,tauaer(:,:,:,jbs:jbe),              &
         !   &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
         !   &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
         !   &       minorfrac,scaleminor,scaleminorn2,indminor,                &
         !   &       nlay, ix, myim(jbs:jbe), async_id, smalljj, blockjj,                                                   &
         !   !  ---  outputs:
         !   &       fracs, tautot                                              &
         !   &     )

               !     if (lprnt) then
               !     print *,' after taumol'
               !     do k = 1, nlay
               !       write(6,121) k
               !121    format(' k =',i3,5x,'fracs')
               !       write(6,122) (fracs(j,k),j=1,ngptlw)
               !122    format(10e14.7)
               !       write(6,123) k
               !123    format(' k =',i3,5x,'tautot')
               !       write(6,122) (tautot(j,k),j=1,ngptlw)
               !     enddo
               !     endif

               !  --- ... call the radiative transfer routine based on cloud scheme
               !          selection. clear sky calculation is done at the same time.
         !$acc enter data create(totuclfl, totdclfl, stemp, totuflux, totdflux) async(async_id)


         if (isubclw <= 0) then
            if (iovrlw <= 0) then
               do jj = 1, blockjj
                  jf = jjoffset+jj
                  do iplon = 1, myim(jf) ! lab_do_iplon
                     m = i + nxjp_acc(jf) - jbs_nxjp_acc
                     do k = 1, nlay
                        delp_im(k) = delp(m, k)
                        htr_im(k) = htr(m, k)
                        htrcl_im(k) = htrcl(m, k)
                     end do
                     do k = 0, nlay
                        totuflux_im(k) = totuflux(m, k)
                        totdflux_im(k) = totdflux(m, k)
                        totuclfl_im(k) = totuclfl(m, k)
                        totdclfl_im(k) = totdclfl(m, k)
                     end do
                     do k = 0, nlp1
                        cldfrc_im(k) = cldfrc(m, k)
                     end do
                     do j = 1, nbands
                        semiss_im(j) = semiss(m, j)
                        secdiff_im(j) = secdiff(m, j)
                     end do
                     call rtrn                                                   &
                        !  ---  inputs:
                        &     ( semiss_im,delp_im,cldfrc_im,taucld_im,tautot_im,pklay_im,pklev_im,              &
                        &       fracs_im,secdiff_im,nlay,nlp1,                                   &
                        !  ---  outputs:
                        &       totuflux_im,totdflux_im,htr_im, totuclfl_im,totdclfl_im,htrcl_im, htrb_im       &
                        &     )
                     do k = 1, nlay
                        delp(m, k) = delp_im(k)
                        htr(m, k) = htr_im(k)
                        htrcl(m, k) = htrcl_im(k)
                     end do
                     do k = 0, nlay
                        totuflux(m, k) = totuflux_im(k)
                        totdflux(m, k) = totdflux_im(k)
                        totuclfl(m, k) = totuclfl_im(k)
                        totdclfl(m, k) = totdclfl_im(k)
                     end do
                     do k = 0, nlp1
                        cldfrc(m, k) = cldfrc_im(k)
                     end do
                     do j = 1, nbands
                        semiss(m, j) = semiss_im(j)
                        secdiff(m, j) = secdiff_im(j)
                     end do
                  end do
               end do

            else
               call rtrnmr                                                 &
                  !  ---  inputs:
                  &     (  plvl, plyr, prslk1, tvly, rhly, slmsk, tracer1, xlon, xlat, &
                           lsswr, lslwr, me, myrank, nxptot, my_max, ntrac, &
                        laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
                        o3vmr, gasvmr_other, colbrd,              &
                  &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
                  &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
                  &       minorfrac,scaleminor,scaleminorn2,indminor,                &
                  &       semiss,delp,cldfrc,taucld, stemp, tavel, &
                        tlvl, secdiff,nlay,nlp1, ix, &
                        myim(jbs:jbe), map_jj(jbs_nxjp_acc:jbe_nxjp_acc-1), &
                        map_i(jbs_nxjp_acc:jbe_nxjp_acc-1), nxjp_acc_length, jjoffset, &
                        max_nxjp_acc_length, jbs_nxjp_acc, async_id, ng_array, ns_array, ng00, &
                        smalljj, blockjj, lhtrlwb, kd,  &
                  !  ---  outputs:
                  &       totuflux,totdflux,hlwc(:,:,jbs:jbe), totuclfl, &
                        totdclfl,hlw0(:,:,jbs:jbe), hlwb       &
                  &     )

            endif   ! end if_iovrlw_block

         else
            do jj = 1, blockjj
               jf = jjoffset+jj
               do iplon = 1, myim(jf) ! lab_do_iplon
                  m = i + nxjp_acc(jf) - jbs_nxjp_acc
                  do k = 1, nlay
                     delp_im(k) = delp(m, k)
                     htr_im(k) = htr(m, k)
                     htrcl_im(k) = htrcl(m, k)
                     do j = 1, nbands
                        taucld_im(k, j) = taucld(n, k, j)
                        htrb_im(k, j) = htrb(n, k, j)
                     end do
                     do j = 1, ngptlw
                        !tautot_im(k, j) = tautot(iplon, k, j, jj)
                        !fracs_im(k, j) = fracs(iplon, k, j, jj)
                        cldfmc_im(k, j) = cldfmc(iplon, k, j, jj)
                     end do
                  end do
                  do k = 0, nlay
                     totuflux_im(k) = totuflux(m, k)
                     totdflux_im(k) = totdflux(m, k)
                     totuclfl_im(k) = totuclfl(m, k)
                     totdclfl_im(k) = totdclfl(m, k)
                     do j = 1, nbands
                        pklay(n, k, j) = pklay_im(k, j)
                        pklev(n, k, j) = pklev_im(k, j)
                     end do
                  end do
                  do j = 1, nbands
                     semiss_im(j) = semiss(m, j)
                     secdiff_im(j) = secdiff(m, j)
                  end do
                  call rtrnmc                                                   &
                     !  ---  inputs:
                     &     ( semiss_im,delp_im,cldfmc_im,taucld_im,tautot_im,pklay_im,pklev_im,              &
                     &       fracs_im,secdiff_im,nlay,nlp1,                                   &
                     !  ---  outputs:
                     &       totuflux_im,totdflux_im,htr_im, totuclfl_im,totdclfl_im,htrcl_im, htrb_im       &
                     &     )
                  do k = 1, nlay
                     delp(m, k) = delp_im(k)
                     htr(m, k) = htr_im(k)
                     htrcl(m, k) = htrcl_im(k)
                     do j = 1, nbands
                        taucld(n, k, j) = taucld_im(k, j)
                        htrb(n, k, j) = htrb_im(k, j)
                     end do
                     do j = 1, ngptlw
                        !tautot(iplon, k, j, jj) = tautot_im(k, j)
                        !fracs(iplon, k, j, jj) = fracs_im(k, j)
                        cldfmc(iplon, k, j, jj) = cldfmc_im(k, j)
                     end do
                  end do
                  do k = 0, nlay
                     totuflux(m, k) = totuflux_im(k)
                     totdflux(m, k) = totdflux_im(k)
                     totuclfl(m, k) = totuclfl_im(k)
                     totdclfl(m, k) = totdclfl_im(k)
                     do j = 1, nbands
                        pklay(n, k, j) = pklay_im(k, j)
                        pklev(n, k, j) = pklev_im(k, j)
                     end do
                  end do
                  do j = 1, nbands
                     semiss(m, j) = semiss_im(j)
                     secdiff(m, j) = secdiff_im(j)
                  end do
               end do
            end do
         endif   ! end if_isubclw_block
         !$acc exit data delete(pavel, coldry, cldfrc, tavel, delp, h2ovmr, o3vmr) async(async_id)
         !$acc exit data delete(stemp, taucld, colbrd, semiss, secdiff) async(async_id)
         !$acc exit data delete(fac00, fac01, fac10, fac11, selffac, selffrac, &
         !$acc&      forfac, forfrac, minorfrac, scaleminor, scaleminorn2, laytrop, &
         !$acc&      jp, jt, jt1, indself, indfor, indminor) async(async_id)


               !  --- ...  output total-sky and clear-sky fluxes and heating rates
         !$acc parallel loop private(jf, iplon, m) async(async_id)
         do n = jbs_nxjp_acc, jbs_nxjp_acc + nxjp_acc_length - 1
            jf = map_jj(n)
            iplon = map_i(n)
            m = n - jbs_nxjp_acc + 1
            topflx_upfxc(iplon, jf) = totuflux(m, nlay)
            topflx_upfx0(iplon, jf) = totuclfl(m, nlay)

            sfcflx_upfxc(iplon, jf) = totuflux(m, 0)
            sfcflx_upfx0(iplon, jf) = totuclfl(m, 0)
            sfcflx_dnfxc(iplon, jf) = totdflux(m, 0)
            sfcflx_dnfx0(iplon, jf) = totdclfl(m, 0)
         end do

         if (ivflip == 0) then       ! output from toa to sfc
                  !! --- ...  optional fluxes
            if ( lflxprf ) then
               !$acc parallel loop collapse(2) private(jf, iplon, m, k1) async(async_id)
               do k = 0, nlay
                  do n = jbs_nxjp_acc, jbs_nxjp_acc + nxjp_acc_length - 1
                     jf = map_jj(n)
                     iplon = map_i(n)
                     m = n - jbs_nxjp_acc + 1
                     k1 = nlp1 - k
                     flxprf_upfxc(iplon,k1, jf) = totuflux(m, k)
                     flxprf_dnfxc(iplon,k1, jf) = totdflux(m, k)
                     flxprf_upfx0(iplon,k1, jf) = totuclfl(m, k)
                     flxprf_dnfx0(iplon,k1, jf) = totdclfl(m, k)
                  end do
               end do
            endif

                  !! --- ...  optional clear sky heating rate

                  !! --- ...  optional spectral band heating rate

         else                        ! output from sfc to toa
               !! --- ...  optional fluxes
            if ( lflxprf ) then
               !$acc parallel loop collapse(2) private(jf, iplon, m) async(async_id)
               do k = 0, nlay
                  do n = jbs_nxjp_acc, jbs_nxjp_acc + nxjp_acc_length - 1
                     jf = map_jj(n)
                     iplon = map_i(n)
                     m = n - jbs_nxjp_acc + 1
                     flxprf_upfxc(iplon,k+1, jf) = totuflux(m, k)
                     flxprf_dnfxc(iplon,k+1, jf) = totdflux(m, k)
                     flxprf_upfx0(iplon,k+1, jf) = totuclfl(m, k)
                     flxprf_dnfx0(iplon,k+1, jf) = totdclfl(m, k)
                  end do
               end do
            endif

                  !! --- ...  optional clear sky heating rate

               !! --- ...  optional spectral band heating rate

         endif                       ! if_ivflip
         !$acc exit data delete(totuclfl, totdclfl, totuflux, totdflux) async(async_id)
         

      end do
      if (isubclw > 0) deallocate(cldfmc)
      !$acc end data 
      !!$acc wait(async_id)

!...................................
      end subroutine lwrad_gpu
!-----------------------------------



!-----------------------------------
      subroutine rlwinit_gpu                                                &
!...................................
!  ---  inputs:
     &     ( me , myrank)
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
!   ilwrate - heating rate unit selections                              !
!           =1: output in k/day                                         !
!           =2: output in k/second                                      !
!   ilwrgas - control flag for rare gases (ch4,n2o,o2,cfcs, etc.)       !
!           =0: do not include rare gases                               !
!           >0: include all rare gases                                  !
!   ilwcliq - liquid cloud optical properties contrl flag               !
!           =0: input cloud opt depth from diagnostic scheme            !
!           >0: input cwp,rew, and other cloud content parameters       !
!   isubclw - sub-column cloud approximation control flag               !
!           =0: no sub-col cld treatment, use grid-mean cld quantities  !
!           =1: mcica sub-col, prescribed seeds to get random numbers   !
!           =2: mcica sub-col, providing array icseed for random numbers!
!   icldflg - cloud scheme control flag                                 !
!           =0: diagnostic scheme gives cloud tau, omiga, and g.        !
!           =1: prognostic scheme gives cloud liq/ice path, etc.        !
!   iovrlw  - clouds vertical overlapping control flag                  !
!           =0: random overlapping clouds                               !
!           =1: maximum/random overlapping clouds                       !
!           =2: maximum overlap cloud (isubcol>0 only)                  !
!                                                                       !
!  *******************************************************************  !
!  original code description                                            !
!                                                                       !
!  original version:       michael j. iacono; july, 1998                !
!  first revision for ncar ccm:               september, 1998           !
!  second revision for rrtm_v3.0:             september, 2002           !
!                                                                       !
!  this subroutine performs calculations necessary for the initialization
!  of the longwave model.  lookup tables are computed for use in the lw !
!  radiative transfer, and input absorption coefficient data for each   !
!  spectral band are reduced from 256 g-point intervals to 140.         !
!                                                                       !
!  *******************************************************************  !
!                                                                       !
! definitions:                                                          !
!   arrays for 10000-point look-up tables:                              !
!   tau_tbl - clear-sky optical depth (used in cloudy radiative transfer!
!   exp_tbl - exponential lookup table for tansmittance                 !
!   tfn_tbl - tau transition function; i.e. the transition of the planck!
!             function from that for the mean layer temperature to that !
!             for the layer boundary temperature as a function of optical
!             depth. the "linear in tau" method is used to make the table
!                                                                       !
!  *******************************************************************  !
!                                                                       !
!  ======================  end of description block  =================  !

!  ---  inputs:
      integer, intent(in) :: me, myrank

!  ---  outputs: none

!  ---  locals:
      real (kind=kind_phys), parameter :: expeps = 1.e-20

      real (kind=kind_phys) :: tfn, pival, explimit

      integer               :: i

!
!===> ... begin here
!
      if ( iovrlw<0 .or. iovrlw>2 ) then
        print *,'  *** error in specification of cloud overlap flag',   &
     &          ' iovrlw=',iovrlw,' in rlwinit !!'
        stop
      elseif ( iovrlw==2 .and. isubclw==0 ) then
        if (me == 0 .and. myrank == 0 ) then
          print *,'  *** iovrlw=2 - maximum cloud overlap, is not yet', &
     &          ' available for isubclw=0 setting!!'
          print *,'      the program uses maximum/random overlap',      &
     &          ' instead.'
        endif

        iovrlw = 1
      endif

      if (me == 0 .and. myrank == 0 ) then
        print *,' - using aer longwave radiation, version: ', vtaglw

        if (ilwrgas > 0) then
          print *,'   --- include rare gases n2o, ch4, o2, cfcs ',      &
     &            'absorptions in lw'
        else
          print *,'   --- rare gases effect is not included in lw'
        endif

        if ( isubclw == 0 ) then
          print *,'   --- using standard grid average clouds, no ',     &
     &            'sub-column clouds approximation applied'
        elseif ( isubclw == 1 ) then
          print *,'   --- using mcica sub-colum clouds approximation ', &
     &            'with a prescribed sequence of permutaion seeds'
        elseif ( isubclw == 2 ) then
          print *,'   --- using mcica sub-colum clouds approximation ', &
     &            'with provided input array of permutation seeds'
        else
          print *,'  *** error in specification of sub-column cloud ',  &
     &            ' control flag isubclw =',isubclw,' !!'
          stop
        endif
      endif

!  --- ...  check cloud flags for consistency

      if ((icldflg == 0 .and. ilwcliq /= 0) .or.                        &
     &    (icldflg == 1 .and. ilwcliq == 0)) then
        print *,'  *** model cloud scheme inconsistent with lw',        &
     &          ' radiation cloud radiative property setup !!'
        stop
      endif

!  --- ...  setup default surface emissivity for each band here

      semiss0(:) = f_one

!  --- ...  setup constant factors for flux and heating rate
!           the 1.0e-2 is to convert pressure from mb to n/m**2

      pival = 2.0 * asin(f_one)
      fluxfac = pival * 2.0d4
!     fluxfac = 62831.85307179586                   ! = 2 * pi * 1.0e4

      if (ilwrate == 1) then
!       heatfac = 8.4391
!       heatfac = con_g * 86400. * 1.0e-2 / con_cp  !   (in k/day)
        heatfac = con_g * 864.0 / con_cp            !   (in k/day)
      else
        heatfac = con_g * 1.0e-2 / con_cp           !   (in k/second)
      endif

!  --- ...  compute lookup tables for transmittance, tau transition
!           function, and clear sky tau (for the cloudy sky radiative
!           transfer).  tau is computed as a function of the tau
!           transition function, transmittance is calculated as a
!           function of tau, and the tau transition function is
!           calculated using the linear in tau formulation at values of
!           tau above 0.01.  tf is approximated as tau/6 for tau < 0.01.
!           all tables are computed at intervals of 0.001.  the inverse
!           of the constant used in the pade approximation to the tau
!           transition function is set to b.

      tau_tbl(0) = f_zero
      exp_tbl(0) = f_one
      tfn_tbl(0) = f_zero

      tau_tbl(ntbl) = 1.e10
      exp_tbl(ntbl) = expeps
      tfn_tbl(ntbl) = f_one

      explimit = aint( -log(tiny(exp_tbl(0))) )

      do i = 1, ntbl-1
!org    tfn = float(i) / float(ntbl)
!org    tau_tbl(i) = bpade * tfn / (f_one - tfn)
        tfn = real(i, kind_phys) / real(ntbl-i, kind_phys)
        tau_tbl(i) = bpade * tfn
        if (tau_tbl(i) >= explimit) then
          exp_tbl(i) = expeps
        else
          exp_tbl(i) = exp( -tau_tbl(i) )
        endif

        if (tau_tbl(i) < 0.06) then
          tfn_tbl(i) = tau_tbl(i) / 6.0
        else
          tfn_tbl(i) = f_one - 2.0*( (f_one / tau_tbl(i))               &
     &               - ( exp_tbl(i) / (f_one - exp_tbl(i)) ) )
        endif
      enddo
      !$acc enter data copyin(nspa, nspb, ngb)

!...................................
      end subroutine rlwinit_gpu
!-----------------------------------


! ----------------------------
      subroutine cldprop                                                &
! ............................
!  ---  inputs:
     &     ( cfrac,cliqp,reliq,cicep,reice,cdat1,cdat2,cdat3,cdat4,     &
     &       nlay, nlp1, ipseed, ix, myim, lcf1, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                                          &
!  ---  outputs:
     &       taucld                                             &
     &     )

!  ===================  program usage description  ===================  !
!                                                                       !
! purpose:  compute the cloud optical depth(s) for each cloudy layer    !
! and g-point interval.                                                 !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       -size- !
!    cfrac - real, layer cloud fraction                          0:nlp1 !
!        .....  for ilwcliq > 0  (prognostic cloud sckeme)  - - -       !
!    cliqp - real, layer in-cloud liq water path (g/m**2)          nlay !
!    reliq - real, mean eff radius for liq cloud (micron)          nlay !
!    cicep - real, layer in-cloud ice water path (g/m**2)          nlay !
!    reice - real, mean eff radius for ice cloud (micron)          nlay !
!    cdat1 - real, layer rain drop water path  (g/m**2)            nlay !
!    cdat2 - real, effective radius for rain drop (microm)         nlay !
!    cdat3 - real, layer snow flake water path (g/m**2)            nlay !
!    cdat4 - real, effective radius for snow flakes (micron)       nlay !
!        .....  for ilwcliq = 0  (diagnostic cloud sckeme)  - - -       !
!    cdat1 - real, input cloud optical depth                       nlay !
!    cdat2 - real, layer cloud single scattering albedo            nlay !
!    cdat3 - real, layer cloud asymmetry factor                    nlay !
!    cdat4 - real, optional use                                    nlay !
!    cliqp - not used                                              nlay !
!    reliq - not used                                              nlay !
!    cicep - not used                                              nlay !
!    reice - not used                                              nlay !
!                                                                       !
!    nlay  - integer, number of vertical layers                      1  !
!    nlp1  - integer, number of vertical levels                      1  !
!    ipseed- permutation seed for generating random numbers (isubclw>0) !
!                                                                       !
!  outputs:                                                             !
!    cldfmc - real, cloud fraction for each sub-column       ngptlw*nlay!
!    taucld - real, cld opt depth for bands (non-mcica)      nbands*nlay!
!                                                                       !
!  explanation of the method for each value of ilwcliq, and ilwcice.    !
!    set up in module "module_radlw_cntr_para"                          !
!                                                                       !
!     ilwcliq=0  : input cloud optical property (tau, ssa, asy).        !
!                  (used for diagnostic cloud method)                   !
!     ilwcliq>0  : input cloud liq/ice path and effective radius, also  !
!                  require the user of 'ilwcice' to specify the method  !
!                  used to compute aborption due to water/ice parts.    !
!  ...................................................................  !
!                                                                       !
!     ilwcliq=1:   the water droplet effective radius (microns) is input!
!                  and the opt depths due to water clouds are computed  !
!                  as in hu and stamnes, j., clim., 6, 728-742, (1993). !
!                  the values for absorption coefficients appropriate for
!                  the spectral bands in rrtm have been obtained for a  !
!                  range of effective radii by an averaging procedure   !
!                  based on the work of j. pinto (private communication).
!                  linear interpolation is used to get the absorption   !
!                  coefficients for the input effective radius.         !
!                                                                       !
!     ilwcice=1:   the cloud ice path (g/m2) and ice effective radius   !
!                  (microns) are input and the optical depths due to ice!
!                  clouds are computed as in ebert and curry, jgr, 97,  !
!                  3831-3836 (1992).  the spectral regions in this work !
!                  have been matched with the spectral bands in rrtm to !
!                  as great an extent as possible:                      !
!                     e&c 1      ib = 5      rrtm bands 9-16            !
!                     e&c 2      ib = 4      rrtm bands 6-8             !
!                     e&c 3      ib = 3      rrtm bands 3-5             !
!                     e&c 4      ib = 2      rrtm band 2                !
!                     e&c 5      ib = 1      rrtm band 1                !
!     ilwcice=2:   the cloud ice path (g/m2) and ice effective radius   !
!                  (microns) are input and the optical depths due to ice!
!                  clouds are computed as in rt code, streamer v3.0     !
!                  (ref: key j., streamer user's guide, cooperative     !
!                  institute for meteorological satellite studies, 2001,!
!                  96 pp.) valid range of values for re are between 5.0 !
!                  and 131.0 micron.                                    !
!     ilwcice=3:   the ice generalized effective size (dge) is input and!
!                  the optical properties, are calculated as in q. fu,  !
!                  j. climate, (1998). q. fu provided high resolution   !
!                  tales which were appropriately averaged for the bands!
!                  in rrtm_lw. linear interpolation is used to get the  !
!                  coeff from the stored tables. valid range of values  !
!                  for deg are between 5.0 and 140.0 micron.            !
!                                                                       !
!  other cloud control module variables:                                !
!     isubclw =0: standard cloud scheme, no sub-col cloud approximation !
!             >0: mcica sub-col cloud scheme using ipseed as permutation!
!                 seed for generating rundom numbers                    !
!                                                                       !
!  ======================  end of description block  =================  !
!
      use module_radlw_cldprlw

!  ---  inputs:
      integer, intent(in) :: nlay, nlp1, ipseed(max_nxjp_acc_length), ix, myim(fulljj), &
         fulljj, blockjj, nxjp_acc_length, jjoffset, max_nxjp_acc_length, jbs_nxjp_acc
      integer, dimension(nxjp_acc_length), intent(in) :: map_jj, map_i

      real (kind=kind_phys), dimension(max_nxjp_acc_length, 0:nlp1), intent(in) :: cfrac
      real (kind=kind_phys), dimension(nxptot, nlay),   intent(in) :: cliqp,    &
     &       reliq, cicep, reice, cdat1, cdat2, cdat3, cdat4
     logical, dimension(max_nxjp_acc_length), intent(in) :: lcf1

!  ---  outputs:
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay, nbands),intent(out):: taucld

!  ---  locals:
      real (kind=kind_phys), dimension(nbands) :: tauliq, tauice
      real (kind=kind_phys), allocatable, dimension(:,:,:)   :: cldf

      real (kind=kind_phys) :: dgeice, factor, fint, tauran, tausnw,    &
     &       cldliq, refliq, cldice, refice

      logical,allocatable,dimension(:,:,:,:) :: lcloudy
      integer :: ia, ib, ig, k, index, iplon, jj

      real (kind=kind_phys), allocatable, dimension(:)   :: cldf_im
      logical,allocatable,dimension(:,:) :: lcloudy_im

      integer :: async_id, n, rr
!
!===> ...  begin here
!
      !$acc parallel loop collapse(2) async(async_id)
      do ib = 1, nbands
         do k = 1, nlay
            do n = 1, nxjp_acc_length
               taucld(n, k,ib) = f_zero
            enddo
         end do
      end do

               !  --- ...  compute cloud radiative properties for a cloudy column
      if (ilwcliq > 0) then ! lab_if_ilwcliq


         !  --- ...  calculation of absorption coefficients due to ice clouds.
         !$acc parallel loop collapse(2) private(cldice, refice, dgeice, factor, index, &
         !$acc&         fint, tausnw, tauran, tauliq, tauice, rr) async(async_id)
         do k = 1, nlay ! lab_do_k
            do n = 1, nxjp_acc_length
               rr = n + jbs_nxjp_acc - 1
               if ( lcf1(n) ) then
                  if (cfrac(n, k) > cldmin) then ! lab_if_cld
                     cldliq = cliqp(rr, k)
                     !           refliq = max(2.5e0, min(60.0e0, reliq(k) ))
                     !           refice = max(5.0e0, reice(k) )
                     refliq = reliq(rr, k)

                     !  --- ...  calculation of absorption coefficients due to water clouds.

                     if ( cldliq <= f_zero ) then
                        !$acc loop seq
                        do ib = 1, nbands
                           tauliq(ib) = f_zero
                        enddo
                     else
                        if ( ilwcliq == 1 ) then
                           factor = refliq - 1.5
                           index  = max( 1, min( 57, int( factor ) ))
                           fint   = factor - float(index)
                           !$acc loop seq
                           do ib = 1, nbands
                              tauliq(ib) = max(f_zero, cldliq*(absliq1(index,ib)    &
                              &              + fint*(absliq1(index+1,ib)-absliq1(index,ib)) ))
                           enddo
                        end if
                     endif   ! end if_ilwcliq_block


                     cldice = cicep(rr, k)
                     refice = reice(rr, k)
                     if ( cldice <= f_zero ) then
                        !$acc loop seq
                        do ib = 1, nbands
                           tauice(ib) = f_zero
                        enddo
                     else
                        if ( ilwcice == 1 ) then
                           !  --- ...  ebert and curry approach for all particle sizes though somewhat
                           !           unjustified for large ice particles

                           refice = min(130.0, max(13.0, real(refice) ))
                           !$acc loop seq
                           do ib = 1, nbands
                              ia = ipat(ib)             ! eb_&_c band index for ice cloud coeff
                              tauice(ib) = max(f_zero, cldice*(absice1(1,ia)        &
                              &                         + absice1(2,ia)/refice) )
                           enddo
                        elseif ( ilwcice == 2 ) then
                           !  --- ...  streamer approach for ice effective radius between 5.0 and 131.0 microns
                           !           and ebert and curry approach for ice eff radius greater than 131.0 microns.
                           !           no smoothing between the transition of the two methods.
                           factor = (refice - 2.0) / 3.0
                           index  = max( 1, min( 42, int( factor ) ))
                           fint   = factor - float(index)
                           !$acc loop seq
                           do ib = 1, nbands
                              tauice(ib) = max(f_zero, cldice*(absice2(index,ib)    &
                              &              + fint*(absice2(index+1,ib) - absice2(index,ib)) ))
                           enddo
                        elseif ( ilwcice == 3 ) then
                           !  --- ...  fu's approach for ice effective radius between 4.8 and 135 microns
                           !           (generalized effective size from 5 to 140 microns)


                           !               dgeice = max(5.0, 1.5396*refice)              ! v4.4 value
                           dgeice = max(5.0, 1.0315*refice)              ! v4.71 value
                           factor = (dgeice - 2.0) / 3.0
                           index  = max( 1, min( 45, int( factor ) ))
                           fint   = factor - float(index)
                           !$acc loop seq
                           do ib = 1, nbands
                              tauice(ib) = max(f_zero, cldice*(absice3(index,ib)    &
                              &              + fint*(absice3(index+1,ib) - absice3(index,ib)) ))
                           enddo
                        end if

                     endif   ! end if_cldice_block
                     if (cdat3(rr, k)>f_zero .and. cdat4(rr, k)>10.0_kind_phys) then
                        tausnw = abssnow0*1.05756*cdat3(rr, k)/cdat4(rr, k)      ! fu's formula
                     else
                        tausnw = f_zero
                     endif
                     tauran = absrain * cdat1(rr, k)                      ! ncar formula
                     !$acc loop seq
                     do ib = 1, nbands
                        taucld(n, k,ib) = tauice(ib) + tauliq(ib) + tauran + tausnw
                     enddo
                  end if
               end if
            end do
         end do

      else  ! lab_if_ilwcliq
         !$acc parallel loop collapse(2) async(async_id)
         do k = 1, nlay
            do n = 1, nxjp_acc_length
               if ( lcf1(n) ) then
                  if (cfrac(n, k) > cldmin) then
                     do ib = 1, nbands
                        taucld(n, k,ib) = cdat1(rr, k)
                     enddo
                  endif
               end if
            end do
         end do
      endif  ! lab_if_ilwcliq


      return
! ..................................
      end subroutine cldprop
! ----------------------------------


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
!   lcloudy - logical, sub-colum cloud profile flag array    ngptlw*nlay!
!                                                                       !
!  other control flags from module variables:                           !
!     iovrlw    : control flag for cloud overlapping method             !
!                 =0:random; =1:maximum/random: =2:maximum              !
!                                                                       !
!  =====================    end of definitions    ====================  !

      implicit none

!  ---  inputs:
      integer, intent(in) :: nlay, ipseed

      real (kind=kind_phys), dimension(nlay), intent(in) :: cldf

!  ---  outputs:
      logical, dimension(ngptlw,nlay), intent(out) :: lcloudy

!  ---  locals:
      real (kind=kind_phys) :: cdfunc(ngptlw,nlay), rand1d(ngptlw),     &
     &       rand2d(nlay*ngptlw), tem1

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

      select case ( iovrlw )

        case( 0 )        ! random overlap, pick a random value at every level

          call random_number                                            &
!  ---  inputs: ( none )
!  ---  outputs:
     &     ( rand2d, stat )

          k1 = 0
          do n = 1, ngptlw
            do k = 1, nlay
              k1 = k1 + 1
              cdfunc(n,k) = rand2d(k1)
            enddo
          enddo

        case( 1 )        ! max-ran overlap

          call random_number                                            &
!  ---  inputs: ( none )
!  ---  outputs:
     &     ( rand2d, stat )

          k1 = 0
          do n = 1, ngptlw
            do k = 1, nlay
              k1 = k1 + 1
              cdfunc(n,k) = rand2d(k1)
            enddo
          enddo

!  ---  first pick a random number for bottom (or top) layer.
!       then walk up the column: (aer's code)
!       if layer below is cloudy, use the same rand num in the layer below
!       if layer below is clear,  use a new random number

!  ---  from bottom up
          do k = 2, nlay
            k1 = k - 1
            tem1 = f_one - cldf(k1)

            do n = 1, ngptlw
              if ( cdfunc(n,k1) > tem1 ) then
                cdfunc(n,k) = cdfunc(n,k1)
              else
                cdfunc(n,k) = cdfunc(n,k) * tem1
              endif
            enddo
          enddo

!  ---  or walk down the column: (if use original author's method)
!       if layer above is cloudy, use the same rand num in the layer above
!       if layer above is clear,  use a new random number

!  ---  from top down
!         do k = nlay-1, 1, -1
!           k1 = k + 1
!           tem1 = f_one - cldf(k1)

!           do n = 1, ngptlw
!             if ( cdfunc(n,k1) > tem1 ) then
!               cdfunc(n,k) = cdfunc(n,k1)
!             else
!               cdfunc(n,k) = cdfunc(n,k) * tem1
!             endif
!           enddo
!         enddo

        case( 2 )        ! maximum overlap, pick same random numebr at every level

          call random_number                                            &
!  ---  inputs: ( none )
!  ---  outputs:
     &     ( rand1d, stat )

          do n = 1, ngptlw
            tem1 = rand1d(n)

            do k = 1, nlay
              cdfunc(n,k) = tem1
            enddo
          enddo

      end select

!  --- ...  generate subcolumns for homogeneous clouds

      do k = 1, nlay
        tem1 = f_one - cldf(k)

        do n = 1, ngptlw
          lcloudy(n,k) = cdfunc(n,k) >= tem1
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
     &     ( pavel,tavel,stemp,nf_vgas, h2ovmr, gasvmr_co2, o3vmr, gasvmr_other,coldry,colbrd,          &
     &       nlay, nlp1, ix, myim, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,   &
!  ---  outputs:
     &       laytrop,pklay,pklev,jp,jt,jt1,                             &
     &       fac00,fac01,fac10,fac11,                            &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor                 &
     &     )

!  ===================  program usage description  ===================  !
!                                                                       !
! purpose:  compute various coefficients needed in radiative transfer   !
!    calculations.                                                      !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       -size- !
!   pavel     - real, layer pressures (mb)                         nlay !
!   tavel     - real, layer temperatures (k)                       nlay !
!   tz        - real, level (interface) temperatures (k)         0:nlay !
!   stemp     - real, surface ground temperature (k)                1   !
!   h2ovmr    - real, layer w.v. volum mixing ratio (kg/kg)        nlay !
!   colamt    - real, column amounts of absorbing gases      nlay*maxgas!
!                 2nd indices range: 1-maxgas, for watervapor,          !
!                 carbon dioxide, ozone, nitrous oxide, methane,        !
!                 oxigen, carbon monoxide,etc. (molecules/cm**2)        !
!   coldry    - real, dry air column amount                        nlay !
!   colbrd    - real, column amount of broadening gases            nlay !
!   nlay/nlp1 - integer, total number of vertical layers, levels    1   !
!                                                                       !
!  outputs:                                                             !
!   laytrop   - integer, tropopause layer index (unitless)          1   !
!   pklay     - real, integrated planck func at lay temp   nbands*0:nlay!
!   pklev     - real, integrated planck func at lev temp   nbands*0:nlay!
!   jp        - real, indices of lower reference pressure          nlay !
!   jt, jt1   - real, indices of lower reference temperatures      nlay !
!   rfrate    - real, ref ratios of binary species param   nlay*nrates*2!
!     (:,m,:)m=1-h2o/co2,2-h2o/o3,3-h2o/n2o,4-h2o/ch4,5-n2o/co2,6-o3/co2!
!     (:,:,n)n=1,2: the rates of ref press at the 2 sides of the layer  !
!   facij     - real, factors multiply the reference ks,           nlay !
!                 i,j=0/1 for lower/higher of the 2 appropriate         !
!                 temperatures and altitudes.                           !
!   selffac   - real, scale factor for w. v. self-continuum        nlay !
!                 equals (w. v. density)/(atmospheric density           !
!                 at 296k and 1013 mb)                                  !
!   selffrac  - real, factor for temperature interpolation of      nlay !
!                 reference w. v. self-continuum data                   !
!   indself   - integer, index of lower ref temp for selffac       nlay !
!   forfac    - real, scale factor for w. v. foreign-continuum     nlay !
!   forfrac   - real, factor for temperature interpolation of      nlay !
!                 reference w.v. foreign-continuum data                 !
!   indfor    - integer, index of lower ref temp for forfac        nlay !
!   minorfrac - real, factor for minor gases                       nlay !
!   scaleminor,scaleminorn2                                             !
!             - real, scale factors for minor gases                nlay !
!   indminor  - integer, index of lower ref temp for minor gases   nlay !
!                                                                       !
!  ======================    end of definitions    ===================  !

!  ---  inputs:
      integer, intent(in) :: nlay, nlp1, ix,  myim(fulljj), fulljj, blockjj, &
         nxjp_acc_length, jjoffset, max_nxjp_acc_length, jbs_nxjp_acc
      integer, dimension(nxjp_acc_length), intent(in) :: map_jj, map_i
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: h2ovmr, o3vmr
      real (kind=kind_phys), dimension(nxptot,nlay), intent(in) :: gasvmr_co2
      real (kind=kind_phys), dimension(nf_vgas), intent(in) :: gasvmr_other
      integer, intent(in) :: nf_vgas

      real (kind=kind_phys) :: colamt1

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: tavel
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: pavel, coldry, colbrd

      real (kind=kind_phys), dimension(nxptot), intent(in) :: stemp

!  ---  outputs:
      integer, dimension(max_nxjp_acc_length, nlay), intent(out) :: jp, jt, jt1, indself,    &
     &       indfor, indminor

      integer, dimension(max_nxjp_acc_length), intent(out) :: laytrop

      real (kind=kind_phys), dimension(max_nxjp_acc_length,0:nlay, nbands), intent(out) ::   &
     &       pklev, pklay

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay),          intent(out) ::   &
     &       fac00, fac01, fac10, fac11, selffac, selffrac, forfac,     &
     &       forfrac, minorfrac, scaleminor, scaleminorn2

!  ---  locals:
      real (kind=kind_phys) :: tlvlfr, tlyrfr, plog, fp, ft, ft1,       &
     &       tem1, tem2, tavelr, tzr, pavelr, forfacr

      integer :: i, k, jp1, indlev, indlay, iplon, jj, async_id, jpr, jtr, jt1r, indminorr, tzi, taveli, n
      ! GPU: variable name changed: CPU - colamt,       GPU - colamt1

!
!===> ... begin here
!
!  --- ...  calculate information needed by the radiative transfer routine
!           that is specific to this atmosphere, especially some of the
!           coefficients and indices needed to compute the optical depths
!           by interpolating data from stored reference atmospheres.

            !  --- ...  begin layer loop
            !           calculate the integrated planck functions for each band at the
            !           surface, level, and layer temperatures.
      !$acc parallel loop private(jpr) async(async_id)
      do n = 1, nxjp_acc_length
         jpr = 0
         !$acc loop seq
         do k = 1, nlay
            if (log(pavel(n, k)) > 4.56) jpr = jpr + 1
         end do
         laytrop(n) = jpr
      end do


      
      !$acc parallel loop collapse(2) private(plog, &
      !$acc&         jp1, fp, tem1, tem2, ft, ft1, tavelr, pavelr, jpr, jtr, jt1r, &
      !$acc&         forfacr, indminorr) async(async_id)
      do k = 1, nlay
         do n = 1, nxjp_acc_length
            tavelr = tavel(n, k)
            pavelr = pavel(n, k)



            !  --- ...  find the two reference pressures on either side of the
            !           layer pressure. store them in jp and jp1. store in fp the
            !           fraction of the difference (in ln(pressure)) between these
            !           two values that the layer pressure lies.

            plog = log(pavelr)
            jpr = max(1, min(58, int(36.0 - 5.0*(plog+0.04)) ))
            jp1  = jpr + 1
            !  --- ...  limit pressure extrapolation at the top
            fp   = max(f_zero, min(f_one, 5.0*(preflog(jpr)-plog) ))
            !org    fp   = 5.0 * (preflog(jp(k)) - plog)

            !  --- ...  determine, for each reference pressure (jp and jp1), which
            !           reference temperature (these are different for each
            !           reference pressure) is nearest the layer temperature but does
            !           not exceed it. store these indices in jt and jt1, resp.
            !           store in ft (resp. ft1) the fraction of the way between jt
            !           (jt1) and the next highest reference temperature that the
            !           layer temperature falls.

            tem1 = (tavelr-tref(jpr)) / 15.0
            tem2 = (tavelr-tref(jp1  )) / 15.0
            jtr = max(1, min(4, int(3.0 + tem1) ))
            jt1r = max(1, min(4, int(3.0 + tem2) ))
            !  --- ...  restrict extrapolation ranges by limiting abs(det t) < 37.5 deg
            ft  = max(-0.5, min(1.5, tem1 - float(jtr - 3) ))
            ft1 = max(-0.5, min(1.5, tem2 - float(jt1r - 3) ))
            !org    ft  = tem1 - float(jt (k) - 3)
            !org    ft1 = tem2 - float(jt1(k) - 3)

            !  --- ...  we have now isolated the layer ln pressure and temperature,
            !           between two reference pressures and two reference temperatures
            !           (for each reference pressure).  we multiply the pressure
            !           fraction fp with the appropriate temperature fractions to get
            !           the factors that will be needed for the interpolation that yields
            !           the optical depths (performed in routines taugbn for band n)

            tem1 = f_one - fp
            fac10(n, k) = tem1 * ft
            fac00(n, k) = tem1 * (f_one - ft)
            fac11(n, k) = fp * ft1
            fac01(n, k) = fp * (f_one - ft1)

            forfacr = pavelr*stpfac / (tavelr*(1.0 + h2ovmr(n, k)))
            selffac(n, k) = h2ovmr(n, k) * forfacr

            !  --- ...  set up factors needed to separately include the minor gases
            !           in the calculation of absorption coefficient
            colamt1 = max(f_zero,    coldry(n, k)*h2ovmr(n, k))          ! h2o
            scaleminor(n, k) = pavelr / tavelr
            scaleminorn2(n, k) = (pavelr / tavelr)                         &
            &                  * (colbrd(n, k)/(coldry(n, k) + colamt1))
            tem1 = (tavelr - 180.8) / 7.2
            indminorr = min(18, max(1, int(tem1)))
            minorfrac(n, k) = tem1 - float(indminorr)

            !  --- ...  if the pressure is less than ~100mb, perform a different
            !           set of species interpolations.
            if (plog > 4.56) then

               tem1 = (332.0 - tavel(n, k)) / 36.0
               indfor(n, k) = min(2, max(1, int(tem1)))
               forfrac(n, k) = tem1 - float(indfor(n, k))

               !  --- ...  set up factors needed to separately include the water vapor
               !           self-continuum in the calculation of absorption coefficient.

               tem1 = (tavel(n, k) - 188.0) / 7.2
               indself(n, k) = min(9, max(1, int(tem1)-7))
               selffrac(n, k) = tem1 - float(indself(n, k) + 7)

            else

               tem1 = (tavel(n, k) - 188.0) / 36.0
               indfor(n, k) = 3
               forfrac(n, k) = tem1 - f_one

               indself(n, k) = 0
               selffrac(n, k) = f_zero

            endif

            !  --- ...  rescale selffac and forfac for use in taumol

            selffac(n, k) = colamt1 * selffac(n, k)
            forfacr = colamt1 * forfacr
            jp(n, k) = jpr
            jt(n, k) = jtr
            jt1(n, k) = jt1r
            forfac(n, k) = forfacr
            indminor(n, k) = indminorr
         end do
      end do

      return
! ..................................
      end subroutine setcoef
! ----------------------------------



! ----------------------------------
      subroutine rtrn                                                   &
! ..................................
!  ---  inputs:
     &     ( semiss,delp,cldfrc,taucld,tautot,pklay,pklev,              &
     &       fracs,secdif, nlay,nlp1,                                   &
!  ---  outputs:
     &       totuflux,totdflux,htr, totuclfl,totdclfl,htrcl, htrb       &
     &     )

!  ===================  program usage description  ===================  !
!                                                                       !
! purpose:  compute the upward/downward radiative fluxes, and heating   !
! rates for both clear or cloudy atmosphere.  clouds are assumed as     !
! randomly overlaping in a vertical colum.                              !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                     -size-   !
!   semiss  - real, lw surface emissivity                         nbands!
!   delp    - real, layer pressure thickness (mb)                  nlay !
!   cldfrc  - real, layer cloud fraction                         0:nlp1 !
!   taucld  - real, layer cloud opt depth                    nbands,nlay!
!   tautot  - real, total optical depth (gas+aerosols)       ngptlw,nlay!
!   pklay   - real, integrated planck func at lay temp     nbands*0:nlay!
!   pklev   - real, integrated planck func at lev temp     nbands*0:nlay!
!   fracs   - real, planck fractions                         ngptlw,nlay!
!   secdif  - real, secant of diffusivity angle                   nbands!
!   nlay    - integer, number of vertical layers                    1   !
!   nlp1    - integer, number of vertical levels (interfaces)       1   !
!                                                                       !
!  outputs:                                                             !
!   totuflux- real, total sky upward flux (w/m2)                 0:nlay !
!   totdflux- real, total sky downward flux (w/m2)               0:nlay !
!   htr     - real, total sky heating rate (k/sec or k/day)        nlay !
!   totuclfl- real, clear sky upward flux (w/m2)                 0:nlay !
!   totdclfl- real, clear sky downward flux (w/m2)               0:nlay !
!   htrcl   - real, clear sky heating rate (k/sec or k/day)        nlay !
!   htrb    - real, spectral band lw heating rate (k/day)    nlay*nbands!
!                                                                       !
!  module veriables:                                                    !
!   ngb     - integer, band index for each g-value                ngptlw!
!   fluxfac - real, conversion factor for fluxes (pi*2.e4)           1  !
!   heatfac - real, conversion factor for heating rates (g/cp*1e-2)  1  !
!   tblint  - real, conversion factor for look-up tbl (float(ntbl)   1  !
!   bpade   - real, pade approx constant (1/0.278)                   1  !
!   wtdiff  - real, weight for radiance to flux conversion           1  !
!   ntbl    - integer, dimension of look-up tables                   1  !
!   tau_tbl - real, clr-sky opt dep lookup table                 0:ntbl !
!   exp_tbl - real, transmittance lookup table                   0:ntbl !
!   tfn_tbl - real, tau transition function                      0:ntbl !
!                                                                       !
!  local variables:                                                     !
!    itgas  - integer, index for gases contribution look-up table    1  !
!    ittot  - integer, index for gases plus clouds  look-up table    1  !
!    reflct - real, surface reflectance                              1  !
!    atrgas - real, gaseous absorptivity                             1  !
!    atrtot - real, gaseous and cloud absorptivity                   1  !
!    odcld  - real, cloud optical depth                              1  !
!    efclrfr- real, effective clear sky fraction (1-efcldfr)       nlay !
!    odepth - real, optical depth of gaseous only                    1  !
!    odtot  - real, optical depth of gas and cloud                   1  !
!    gasfac - real, gas-only pade factor, used for planck fn         1  !
!    totfac - real, gas+cld pade factor, used for planck fn          1  !
!    bbdgas - real, gas-only planck function for downward rt         1  !
!    bbugas - real, gas-only planck function for upward rt           1  !
!    bbdtot - real, gas and cloud planck function for downward rt    1  !
!    bbutot - real, gas and cloud planck function for upward rt      1  !
!    gassrcu- real, upwd source radiance due to gas only            nlay!
!    totsrcu- real, upwd source radiance due to gas+cld             nlay!
!    gassrcd- real, dnwd source radiance due to gas only             1  !
!    totsrcd- real, dnwd source radiance due to gas+cld              1  !
!    radtotu- real, spectrally summed total sky upwd radiance        1  !
!    radclru- real, spectrally summed clear sky upwd radiance        1  !
!    radtotd- real, spectrally summed total sky dnwd radiance        1  !
!    radclrd- real, spectrally summed clear sky dnwd radiance        1  !
!    toturad- real, total sky upward radiance by layer     0:nlay*nbands!
!    clrurad- real, clear sky upward radiance by layer     0:nlay*nbands!
!    totdrad- real, total sky downward radiance by layer   0:nlay*nbands!
!    clrdrad- real, clear sky downward radiance by layer   0:nlay*nbands!
!    fnet   - real, net longwave flux (w/m2)                     0:nlay !
!    fnetc  - real, clear sky net longwave flux (w/m2)           0:nlay !
!                                                                       !
!                                                                       !
!  *******************************************************************  !
!  original code description                                            !
!                                                                       !
!  original version:   e. j. mlawer, et al. rrtm_v3.0                   !
!  revision for gcms:  michael j. iacono; october, 2002                 !
!  revision for f90:   michael j. iacono; june, 2006                    !
!                                                                       !
!  this program calculates the upward fluxes, downward fluxes, and      !
!  heating rates for an arbitrary clear or cloudy atmosphere. the input !
!  to this program is the atmospheric profile, all planck function      !
!  information, and the cloud fraction by layer.  a variable diffusivity!
!  angle (secdif) is used for the angle integration. bands 2-3 and 5-9  !
!  use a value for secdif that varies from 1.50 to 1.80 as a function   !
!  of the column water vapor, and other bands use a value of 1.66.  the !
!  gaussian weight appropriate to this angle (wtdiff=0.5) is applied    !
!  here.  note that use of the emissivity angle for the flux integration!
!  can cause errors of 1 to 4 w/m2 within cloudy layers.                !
!  clouds are treated with a random cloud overlap method.               !
!                                                                       !
!  *******************************************************************  !
!  ======================  end of description block  =================  !

!  ---  inputs:
      integer, intent(in) :: nlay, nlp1

      real (kind=kind_phys), dimension(0:nlp1), intent(in) :: cldfrc
      real (kind=kind_phys), dimension(nbands), intent(in) :: semiss,   &
     &       secdif
      real (kind=kind_phys), dimension(nlay),   intent(in) :: delp

      real (kind=kind_phys), dimension(nlay, nbands),intent(in):: taucld
      real (kind=kind_phys), dimension(nlay, ngptlw),intent(in):: fracs, &
     &       tautot

      real (kind=kind_phys), dimension(0:nlay, nbands), intent(in) ::    &
     &       pklev, pklay 

!  ---  outputs:
      real (kind=kind_phys), dimension(nlay), intent(out) :: htr, htrcl

      real (kind=kind_phys), dimension(nlay,nbands),intent(out) :: htrb

      real (kind=kind_phys), dimension(0:nlay), intent(out) ::          &
     &       totuflux, totdflux, totuclfl, totdclfl

!  ---  locals:
      real (kind=kind_phys), parameter :: rec_6 = 0.166667

      real (kind=kind_phys), dimension(0:nlay,nbands) :: clrurad,       &
     &       clrdrad, toturad, totdrad

      real (kind=kind_phys), dimension(nlay)   :: gassrcu, totsrcu,     &
     &       trngas, efclrfr, rfdelp
      real (kind=kind_phys), dimension(0:nlay) :: fnet, fnetc

      real (kind=kind_phys) :: totsrcd, gassrcd, tblind, odepth, odtot, &
     &       odcld, atrtot, atrgas, reflct, totfac, gasfac, flxfac,     &
     &       plfrac, blay, bbdgas, bbdtot, bbugas, bbutot, dplnku,      &
     &       dplnkd, radtotu, radclru, radtotd, radclrd, rad0,          &
     &       clfr, trng, gasu

      integer :: ittot, itgas, ib, ig, k
!
!===> ...  begin here
!
      do ib = 1, nbands
        do k = 0, nlay
          toturad(k,ib) = f_zero
          totdrad(k,ib) = f_zero
          clrurad(k,ib) = f_zero
          clrdrad(k,ib) = f_zero
        enddo
      enddo

      do k = 0, nlay
        totuflux(k) = f_zero
        totdflux(k) = f_zero
        totuclfl(k) = f_zero
        totdclfl(k) = f_zero
      enddo

!  --- ...  loop over all g-points

      do ig = 1, ngptlw
        ib = ngb(ig)

        radtotd = f_zero
        radclrd = f_zero

!  --- ...  downward radiative transfer loop.

        do k = nlay, 1, -1

!  --- ...  clear sky, gases contribution

          odepth = max( f_zero, secdif(ib)*tautot(k,ig) )
          if (odepth <= 0.06) then
            atrgas = odepth - 0.5*odepth*odepth
            trng   = f_one - atrgas
            gasfac = rec_6 * odepth
          else
            tblind = odepth / (bpade + odepth)
            itgas = tblint*tblind + 0.5
            trng  = exp_tbl(itgas)
            atrgas = f_one - trng
            gasfac = tfn_tbl(itgas)
            odepth = tau_tbl(itgas)
          endif

          plfrac = fracs(k,ig)
          blay = pklay(k,ib)

          dplnku = pklev(k,ib  ) - blay
          dplnkd = pklev(k-1,ib) - blay
          bbdgas = plfrac * (blay + dplnkd*gasfac)
          bbugas = plfrac * (blay + dplnku*gasfac)
          gassrcd= bbdgas * atrgas
          gassrcu(k)= bbugas * atrgas
          trngas(k) = trng

!  --- ...  total sky, gases+clouds contribution

          clfr = cldfrc(k)
          if (clfr >= eps) then
!  --- ...  cloudy layer

            odcld = secdif(ib) * taucld(k,ib)
            efclrfr(k) = f_one-(f_one - exp(-odcld))*clfr
            odtot = odepth + odcld
            if (odtot < 0.06) then
              totfac = rec_6 * odtot
              atrtot = odtot - 0.5*odtot*odtot
            else
              tblind = odtot / (bpade + odtot)
              ittot  = tblint*tblind + 0.5
              totfac = tfn_tbl(ittot)
              atrtot = f_one - exp_tbl(ittot)
            endif

            bbdtot = plfrac * (blay + dplnkd*totfac)
            bbutot = plfrac * (blay + dplnku*totfac)
            totsrcd= bbdtot * atrtot
            totsrcu(k)= bbutot * atrtot

!  --- ...  total sky radiance
            radtotd = radtotd*trng*efclrfr(k) + gassrcd                 &
     &              + clfr*(totsrcd - gassrcd)
            totdrad(k-1,ib) = totdrad(k-1,ib) + radtotd

!  --- ...  clear sky radiance
            radclrd = radclrd*trng + gassrcd
            clrdrad(k-1,ib) = clrdrad(k-1,ib) + radclrd

          else
!  --- ...  clear layer

!  --- ...  total sky radiance
            radtotd = radtotd*trng + gassrcd
            totdrad(k-1,ib) = totdrad(k-1,ib) + radtotd

!  --- ...  clear sky radiance
            radclrd = radclrd*trng + gassrcd
            clrdrad(k-1,ib) = clrdrad(k-1,ib) + radclrd

          endif   ! end if_clfr_block

        enddo   ! end do_k_loop

!  --- ...  spectral emissivity & reflectance
!           include the contribution of spectrally varying longwave emissivity
!           and reflection from the surface to the upward radiative transfer.
!     note: spectral and lambertian reflection are identical for the
!           diffusivity angle flux integration used here.

        reflct = f_one - semiss(ib)
        rad0 = semiss(ib) * fracs(1,ig) * pklay(0,ib)

!  --- ...  total sky radiance
        radtotu = rad0 + reflct*radtotd
        toturad(0,ib) = toturad(0,ib) + radtotu

!  --- ...  clear sky radiance
        radclru = rad0 + reflct*radclrd
        clrurad(0,ib) = clrurad(0,ib) + radclru

!  --- ...  upward radiative transfer

        do k = 1, nlay
          clfr = cldfrc(k)
          trng = trngas(k)
          gasu = gassrcu(k)

          if (clfr >= eps) then
!  --- ...  cloudy layer

!  --- ... total sky radiance
            radtotu = radtotu*trng*efclrfr(k) + gasu                    &
     &            + clfr*(totsrcu(k) - gasu)
            toturad(k,ib) = toturad(k,ib) + radtotu

!  --- ... clear sky radiance
            radclru = radclru*trng + gasu
            clrurad(k,ib) = clrurad(k,ib) + radclru

          else
!  --- ...  clear layer

!  --- ... total sky radiance
            radtotu = radtotu*trng + gasu
            toturad(k,ib) = toturad(k,ib) + radtotu

!  --- ... clear sky radiance
            radclru = radclru*trng + gasu
            clrurad(k,ib) = clrurad(k,ib) + radclru

          endif   ! end if_clfr_block

        enddo   ! end do_k_loop

      enddo   ! end do_ig_loop

!  --- ...  process longwave output from band for total and clear streams.
!           calculate upward, downward, and net flux.

      flxfac = wtdiff * fluxfac

      do k = 0, nlay
        do ib = 1, nbands
          totuflux(k) = totuflux(k) + toturad(k,ib)
          totdflux(k) = totdflux(k) + totdrad(k,ib)
          totuclfl(k) = totuclfl(k) + clrurad(k,ib)
          totdclfl(k) = totdclfl(k) + clrdrad(k,ib)
        enddo

        totuflux(k) = totuflux(k) * flxfac
        totdflux(k) = totdflux(k) * flxfac
        totuclfl(k) = totuclfl(k) * flxfac
        totdclfl(k) = totdclfl(k) * flxfac
      enddo

!  --- ...  calculate net fluxes and heating rates
      fnet(0) = totuflux(0) - totdflux(0)

      do k = 1, nlay
        rfdelp(k) = heatfac / delp(k)
        fnet(k) = totuflux(k) - totdflux(k)
        htr (k) = (fnet(k-1) - fnet(k)) * rfdelp(k)
      enddo

!! --- ...  optional clear sky heating rates
      if ( lhlw0 ) then
        fnetc(0) = totuclfl(0) - totdclfl(0)

        do k = 1, nlay
          fnetc(k) = totuclfl(k) - totdclfl(k)
          htrcl(k) = (fnetc(k-1) - fnetc(k)) * rfdelp(k)
        enddo
      endif

!! --- ...  optional spectral band heating rates
      if ( lhlwb ) then
        do ib = 1, nbands
          fnet(0) = (toturad(0,ib) - totdrad(0,ib)) * flxfac

          do k = 1, nlay
            fnet(k) = (toturad(k,ib) - totdrad(k,ib)) * flxfac
            htrb(k,ib) = (fnet(k-1) - fnet(k)) * rfdelp(k)
          enddo
        enddo
      endif

! ..................................
      end subroutine rtrn
! ----------------------------------


! ----------------------------------
      subroutine rtrnmr                                                 &
! ..................................
!  ---  inputs:
      &    ( plvl, plyr, prslk1, tvly, rhly, slmsk, tracer1, xlon, xlat, &
             lsswr, lslwr, me, myrank, nxptot, my_max, ntrac, &
             laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,              &
      &      fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
      &      selffac,selffrac,indself,forfac,forfrac,indfor,            &
      &      minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       semiss,delp,cldfrc,taucld,stemp, tavel, tlvl,               &
     &       secdif, nlay,nlp1, ix, myim, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, ng_array, ns_array, ng00, &
             fulljj, blockjj, lhtrlwb, kd,                                   &
!  ---  outputs:
     &       totuflux,totdflux,hlwc, totuclfl,totdclfl,hlw0, hlwb       &
     &     )

!  ===================  program usage description  ===================  !
!                                                                       !
! purpose:  compute the upward/downward radiative fluxes, and heating   !
! rates for both clear or cloudy atmosphere.  clouds are assumed as in  !
! maximum-randomly overlaping in a vertical colum.                      !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                     -size-   !
!   semiss  - real, lw surface emissivity                         nbands!
!   delp    - real, layer pressure thickness (mb)                  nlay !
!   cldfrc  - real, layer cloud fraction                         0:nlp1 !
!   taucld  - real, layer cloud opt depth                    nbands,nlay!
!   tautot  - real, total optical depth (gas+aerosols)       ngptlw,nlay!
!   pklay   - real, integrated planck func at lay temp     nbands*0:nlay!
!   pklev   - real, integrated planck func at lev temp     nbands*0:nlay!
!   fracs   - real, planck fractions                         ngptlw,nlay!
!   secdif  - real, secant of diffusivity angle                   nbands!
!   nlay    - integer, number of vertical layers                    1   !
!   nlp1    - integer, number of vertical levels (interfaces)       1   !
!                                                                       !
!  outputs:                                                             !
!   totuflux- real, total sky upward flux (w/m2)                 0:nlay !
!   totdflux- real, total sky downward flux (w/m2)               0:nlay !
!   htr     - real, total sky heating rate (k/sec or k/day)        nlay !
!   totuclfl- real, clear sky upward flux (w/m2)                 0:nlay !
!   totdclfl- real, clear sky downward flux (w/m2)               0:nlay !
!   htrcl   - real, clear sky heating rate (k/sec or k/day)        nlay !
!   htrb    - real, spectral band lw heating rate (k/day)    nlay*nbands!
!                                                                       !
!  module veriables:                                                    !
!   ngb     - integer, band index for each g-value                ngptlw!
!   fluxfac - real, conversion factor for fluxes (pi*2.e4)           1  !
!   heatfac - real, conversion factor for heating rates (g/cp*1e-2)  1  !
!   tblint  - real, conversion factor for look-up tbl (float(ntbl)   1  !
!   bpade   - real, pade approx constant (1/0.278)                   1  !
!   wtdiff  - real, weight for radiance to flux conversion           1  !
!   ntbl    - integer, dimension of look-up tables                   1  !
!   tau_tbl - real, clr-sky opt dep lookup table                 0:ntbl !
!   exp_tbl - real, transmittance lookup table                   0:ntbl !
!   tfn_tbl - real, tau transition function                      0:ntbl !
!                                                                       !
!  local variables:                                                     !
!    itgas  - integer, index for gases contribution look-up table    1  !
!    ittot  - integer, index for gases plus clouds  look-up table    1  !
!    reflct - real, surface reflectance                              1  !
!    atrgas - real, gaseous absorptivity                             1  !
!    atrtot - real, gaseous and cloud absorptivity                   1  !
!    odcld  - real, cloud optical depth                              1  !
!    odepth - real, optical depth of gaseous only                    1  !
!    odtot  - real, optical depth of gas and cloud                   1  !
!    gasfac - real, gas-only pade factor, used for planck fn         1  !
!    totfac - real, gas+cld pade factor, used for planck fn          1  !
!    bbdgas - real, gas-only planck function for downward rt         1  !
!    bbugas - real, gas-only planck function for upward rt           1  !
!    bbdtot - real, gas and cloud planck function for downward rt    1  !
!    bbutot - real, gas and cloud planck function for upward rt      1  !
!    gassrcu- real, upwd source radiance due to gas only            nlay!
!    totsrcu- real, upwd source radiance due to gas + cld           nlay!
!    gassrcd- real, dnwd source radiance due to gas only             1  !
!    totsrcd- real, dnwd source radiance due to gas + cld            1  !
!    radtotu- real, spectrally summed total sky upwd radiance        1  !
!    radclru- real, spectrally summed clear sky upwd radiance        1  !
!    radtotd- real, spectrally summed total sky dnwd radiance        1  !
!    radclrd- real, spectrally summed clear sky dnwd radiance        1  !
!    toturad- real, total sky upward radiance by layer     0:nlay*nbands!
!    clrurad- real, clear sky upward radiance by layer     0:nlay*nbands!
!    totdrad- real, total sky downward radiance by layer   0:nlay*nbands!
!    clrdrad- real, clear sky downward radiance by layer   0:nlay*nbands!
!    fnet   - real, net longwave flux (w/m2)                     0:nlay !
!    fnetc  - real, clear sky net longwave flux (w/m2)           0:nlay !
!                                                                       !
!                                                                       !
!  *******************************************************************  !
!  original code description                                            !
!                                                                       !
!  original version:   e. j. mlawer, et al. rrtm_v3.0                   !
!  revision for gcms:  michael j. iacono; october, 2002                 !
!  revision for f90:   michael j. iacono; june, 2006                    !
!                                                                       !
!  this program calculates the upward fluxes, downward fluxes, and      !
!  heating rates for an arbitrary clear or cloudy atmosphere. the input !
!  to this program is the atmospheric profile, all planck function      !
!  information, and the cloud fraction by layer.  a variable diffusivity!
!  angle (secdif) is used for the angle integration. bands 2-3 and 5-9  !
!  use a value for secdif that varies from 1.50 to 1.80 as a function   !
!  of the column water vapor, and other bands use a value of 1.66.  the !
!  gaussian weight appropriate to this angle (wtdiff=0.5) is applied    !
!  here.  note that use of the emissivity angle for the flux integration!
!  can cause errors of 1 to 4 w/m2 within cloudy layers.                !
!  clouds are treated with a maximum-random cloud overlap method.       !
!                                                                       !
!  *******************************************************************  !
!  ======================  end of description block  =================  !

!  ---  inputs for setaer_lw_gpu:
     integer, intent(in) :: me, myrank, nxptot, my_max, ntrac
     logical, intent(in) :: lsswr, lslwr
     real (kind=kind_phys), dimension(ix, my_max),  intent(in) ::  slmsk,  &
     xlon, xlat
     real (kind=kind_phys), dimension(nxptot,nlay), intent(in)  :: rhly, prslk1, tvly
     real (kind=kind_phys), dimension(nxptot,nlay,ntrac), intent(in)   :: tracer1
     real (kind=kind_phys), dimension(nxptot, nlp1), intent(in) :: plvl
     real (kind=kind_phys), dimension(nxptot, nlay), intent(in) :: plyr
!  ---  inputs:
     integer, intent(in) :: nlay, nlp1, ix, myim(fulljj), fulljj, blockjj, kd, &
     nxjp_acc_length, jjoffset, max_nxjp_acc_length, jbs_nxjp_acc, ng00
      integer, intent(in) :: laytrop(max_nxjp_acc_length)
      integer, dimension(nxjp_acc_length), intent(in) :: map_jj, map_i

      integer, dimension(max_nxjp_acc_length, nlay), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: h2ovmr, o3vmr
      real (kind=kind_phys), dimension(nxptot,nlay), intent(in) :: gasvmr_co2
      real (kind=kind_phys), dimension(nf_vgas), intent(in) :: gasvmr_other
      integer, intent(in) :: nf_vgas
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay, nbands) :: tauaer
     
      real (kind=kind_phys), dimension(nxptot, nlp1), intent(in) :: tlvl
      logical, intent(in) :: lhtrlwb
      real (kind=kind_phys), dimension(max_nxjp_acc_length), intent(in) :: stemp
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: tavel
      real (kind=kind_phys), dimension(max_nxjp_acc_length, 0:nlp1), intent(in) :: cldfrc
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nbands), intent(in) :: semiss,   &
     &       secdif
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay),   intent(in) :: delp

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay, nbands),intent(in):: taucld
      real (kind=kind_phys), dimension(nxjp_acc_length, 0:nlay, ng00) :: radtotd_2, radtotu_2


!  ---  outputs:
      real (kind=kind_phys), dimension(ix, nlay, fulljj), intent(out) :: hlwc, hlw0

      real (kind=kind_phys), dimension(ix, nlay,nbands, fulljj),intent(out) :: hlwb

      real (kind=kind_phys), dimension(max_nxjp_acc_length, 0:nlay), intent(out) ::          &
     &       totuflux, totdflux
      real (kind=kind_phys), dimension(max_nxjp_acc_length, 0:nlay), intent(out) :: totuclfl, totdclfl

!  ---  locals:
      real (kind=kind_phys), parameter :: rec_6 = 0.166667

      integer :: small_factor = 1
      integer :: small_fulljj

      real (kind=kind_phys) :: rfdelp

      real (kind=kind_phys) :: totsrcd, gassrcd, tblind, odepth, odtot, &
     &       odcld, atrtot, atrgas, reflct, totfac, gasfac, flxfac,     &
     &       plfrac, blay, bbdgas, bbdtot, bbugas, bbutot, dplnku,      &
     &       dplnkd, radtotu, radclru, rad0, rad,     &
     &       totradd, clrradd, totradu, clrradu, fmax, fmin, rat1, rat2,&
     &       radmod, clfr, trng, trnt, gasu, totu

      integer :: ittot, itgas, ib, ig, k, iplon, jj, iplon2, i2, jj2, k1, n

!  dimensions for cloud overlap adjustment
      real (kind=kind_phys), dimension(nxjp_acc_length, nlp1) :: faccld1u, faccld2u,     &
     &        facclr1u, facclr2u, faccmb1u, faccmb2u
      real (kind=kind_phys), dimension(nxjp_acc_length, 0:nlay) :: faccld1d, faccld2d,   &
     &        facclr1d, facclr2d, faccmb1d, faccmb2d
     real (kind=kind_phys), dimension(nxjp_acc_length, 0:nlay, nbands) :: toturad, totdrad

      integer :: async_id
      real (kind=kind_phys), dimension(nxjp_acc_length, nlay, ng00) :: gassrcu, totsrcu, &
         trngas, trntot, small_fracs, small_tautot
      
      ! GPU register
      logical :: lstcldr
      real (kind=kind_phys) :: fnet, fnet1, totufluxr, totdfluxr, fnet3, fnet4, clfrm, clfrp
      real (kind=kind_phys) :: fnet2r, fnet2r1

      integer, dimension(nbands), intent(in) :: ng_array, ns_array

      ! GPU: variable name changed: CPU - htrb,         GPU - hlwb
      ! GPU: variable name changed: CPU - htr,          GPU - hlwc
      ! GPU: variable name changed: CPU - htrcl,        GPU - hlw0

!
!===> ...  begin here
!
      !small_fulljj = ceiling(float(jlistnum)/float(small_factor))
      
      !$acc enter data create(tauaer) async(async_id)
      !call nvtxStartRange("lw_setaer")
      call setaer_lw_gpu                                                       &
      !  ---  inputs:
      &     ( plvl,plyr,prslk1,tvly,rhly,slmsk,tracer1, xlon,xlat,        &
      &       nlay,nlp1,lsswr,lslwr,me,myrank, ix, &
               map_jj, map_i, &
               nxptot, nxjp_acc_length, jjoffset, max_nxjp_acc_length, jbs_nxjp_acc, async_id,                          &
      !  ---  outputs:
      &       tauaer                                              &
      !    &       faersw,faerlw,aerodp                                       &
      &     )
      !call nvtxEndRange
      !$acc data create(faccld1u, faccld2u, facclr1u, facclr2u, faccmb1u, &
      !$acc&     faccmb2u, faccld1d, faccld2d, facclr1d, facclr2d, faccmb1d, &
      !$acc&     faccmb2d, toturad, totdrad, &
      !$acc&     gassrcu, totsrcu, small_tautot, small_fracs, &
      !$acc&     trngas, trntot, radtotd_2, radtotu_2) async(async_id)
      
      !$acc parallel loop collapse(2) async(async_id)
      do k = 1, nlp1
         do n = 1, nxjp_acc_length
            faccld1u(n, k) = f_zero
            faccld2u(n, k) = f_zero
            facclr1u(n, k) = f_zero
            facclr2u(n, k) = f_zero
            faccmb1u(n, k) = f_zero
            faccmb2u(n, k) = f_zero
         end do
      end do
      
      !$acc parallel loop collapse(2) async(async_id)
      do k = 0, nlay
         do n = 1, nxjp_acc_length
            faccld1d(n, k) = f_zero
            faccld2d(n, k) = f_zero
            facclr1d(n, k) = f_zero
            facclr2d(n, k) = f_zero
            faccmb1d(n, k) = f_zero
            faccmb2d(n, k) = f_zero
         end do
      end do

      !$acc parallel loop private(rat1, rat2, fmax, fmin, &
      !$acc&         lstcldr, clfr, clfrp, clfrm) async(async_id)
      do n = 1, nxjp_acc_length
         clfr = cldfrc(n, 1)
         clfrm = cldfrc(n, 0)
         lstcldr = clfr > eps
         rat1 = f_zero
         rat2 = f_zero
         !$acc loop seq
         do k = 1, nlay-1
            clfrp = cldfrc(n, k+1)
            if (clfr > eps) then
            !  --- ...  maximum/random cloud overlap

               if (clfrp >= clfr) then
                  if (lstcldr) then
                     if (clfr < f_one) then
                        facclr2u(n, k+1) = (clfrp - clfr)               &
                        &                        / (f_one - clfr)
                     endif
                     facclr2u(n, k) = f_zero
                     faccld2u(n, k) = f_zero
                  else
                     fmax = max(clfr, clfrm)
                     if (clfrp > fmax) then
                        facclr1u(n, k+1) = rat2
                        facclr2u(n, k+1) = (clfrp - fmax)/(f_one - fmax)
                     elseif (clfrp < fmax) then
                        facclr1u(n, k+1) = (clfrp - clfr)               &
                        &                        / (clfrm - clfr)
                     else
                        facclr1u(n, k+1) = rat2
                     endif
                  endif

                  if (facclr1u(n, k+1)>f_zero .or. facclr2u(n, k+1)>f_zero) then
                     rat1 = f_one
                     rat2 = f_zero
                  else
                     rat1 = f_zero
                     rat2 = f_zero
                  endif
               else
                  if (lstcldr) then
                     faccld2u(n, k+1) = (clfr - &
                                                   clfrp) / clfr
                     facclr2u(n, k) = f_zero
                     faccld2u(n, k) = f_zero
                  else
                     fmin = min(clfr, clfrm)
                     if (clfrp <= fmin) then
                        faccld1u(n, k+1) = rat1
                        faccld2u(n, k+1) = (fmin - clfrp) / fmin
                     else
                        faccld1u(n, k+1) = (clfr - clfrp)               &
                        &                        / (clfr - fmin)
                     endif
                  endif

                  if (faccld1u(n, k+1)>f_zero .or. faccld2u(n, k+1)>f_zero) then
                     rat1 = f_zero
                     rat2 = f_one
                  else
                     rat1 = f_zero
                     rat2 = f_zero
                  endif
               endif

               faccmb1u(n, k+1) = facclr1u(n, k+1) * &
                                          faccld2u(n, k) * clfrm
               faccmb2u(n, k+1) = faccld1u(n, k+1) * facclr2u(n, k)                   &
               &                  * (f_one - clfrm)
            endif
            lstcldr = clfrp>eps .and. clfr<=eps
            clfrm = clfr
            clfr = clfrp
         enddo
         
         clfr = cldfrc(n, nlay)
         clfrp = cldfrc(n, nlay+1)
         lstcldr = clfr > eps
         rat1 = f_zero
         rat2 = f_zero
         !$acc loop seq
         do k = nlay, 2, -1
            clfrm = cldfrc(n, k-1)
            if (clfr > eps) then

               if (clfrm >= clfr) then
                  if (lstcldr) then
                     if (clfr < f_one) then
                        facclr2d(n, k-1) = (clfrm - clfr)               &
                        &                        / (f_one - clfr)
                     endif

                     facclr2d(n, k) = f_zero
                     faccld2d(n, k) = f_zero
                  else
                     fmax = max(clfr, clfrp)

                     if (clfrm > fmax) then
                        facclr1d(n, k-1) = rat2
                        facclr2d(n, k-1) = (clfrm - fmax) / (f_one - fmax)
                     elseif (clfrm < fmax) then
                        facclr1d(n, k-1) = (clfrm - clfr)               &
                        &                        / (clfrp - clfr)
                     else
                        facclr1d(n, k-1) = rat2
                     endif
                  endif

                  if (facclr1d(n, k-1)>f_zero .or. facclr2d(n, k-1)>f_zero) then
                     rat1 = f_one
                     rat2 = f_zero
                  else
                     rat1 = f_zero
                     rat2 = f_zero
                  endif
               else
                  if (lstcldr) then
                     faccld2d(n, k-1) = (clfr - &
                                                clfrm) / clfr
                     facclr2d(n, k) = f_zero
                     faccld2d(n, k) = f_zero
                  else
                     fmin = min(clfr, clfrp)

                     if (clfrm <= fmin) then
                        faccld1d(n, k-1) = rat1
                        faccld2d(n, k-1) = (fmin - clfrm) / fmin
                     else
                        faccld1d(n, k-1) = (clfr - clfrm)               &
                        &                        / (clfr - fmin)
                     endif
                  endif

                  if (faccld1d(n, k-1)>f_zero .or. faccld2d(n, k-1)>f_zero) then
                     rat1 = f_zero
                     rat2 = f_one
                  else
                     rat1 = f_zero
                     rat2 = f_zero
                  endif
               endif

               faccmb1d(n, k-1) = facclr1d(n, k-1) * &
                                          faccld2d(n, k) * clfrp
               faccmb2d(n, k-1) = faccld1d(n, k-1) * facclr2d(n, k)                   &
               &                  * (f_one - clfrp)
            endif
            lstcldr = clfrm > eps .and. clfr<=eps
            clfrp = clfr
            clfr = clfrm
         enddo
      end do
      
      
      !  --- ...  initialize for radiative transfer.
      !$acc parallel loop collapse(2) async(async_id)
      do ib = 1, nbands
         do n = 1, nxjp_acc_length
            totdrad(n, nlay,ib) = f_zero
         end do
      end do

      !$acc parallel loop collapse(2) async(async_id)
      do k = 0, nlay
         do n = 1, nxjp_acc_length
            totuclfl(n, k) = f_zero
            totdclfl(n, k) = f_zero
         end do
      end do
      

      do ib = 1, nbands
         !call nvtxStartRange("lw_taumol")
         call taumol                                                     &
            !  ---  inputs:
            &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
                    o3vmr, gasvmr_other, colbrd,tauaer,              &
            &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
            &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
            &       minorfrac,scaleminor,scaleminorn2,indminor,                &
            &       nlay, ix, myim, map_jj, map_i, nxjp_acc_length, jjoffset, &
                    max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj, ib, ng00,                                                   &
            !  ---  outputs:
            &       small_fracs, small_tautot                                              &
            &     )
         !call nvtxEndRange
         !call nvtxStartRange("lw_rtrnmr")
         call rtrnmr_ngptlw(secdif, small_tautot, small_fracs, &
            totdrad, totdclfl, facclr1d, faccld1d, faccmb1d, &
            faccmb2d, facclr2d, faccld2d, semiss, toturad, totuclfl, cldfrc, &
            facclr1u, faccld1u, faccmb1u, faccmb2u, facclr2u, faccld2u, taucld, &
            nlay, myim, ix, ngptlw, nlp1, nbands, &
            ns_array(ib), ng_array(ib), ng00, map_jj, map_i, nxjp_acc_length, jjoffset, &
            max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj, &
            gassrcu, totsrcu, &
            trngas, trntot, radtotd_2, radtotu_2, stemp, tavel, tlvl)
         !call nvtxEndRange
      end do
      !$acc exit data delete(tauaer) async(async_id)

      flxfac = wtdiff * fluxfac
      !  --- ...  process longwave output from band for total and clear streams.
      !           calculate upward, downward, and net flux.
      !$acc parallel loop collapse(2) private(totufluxr, totdfluxr) async(async_id) 
      do k = 0, nlay
         do n = 1, nxjp_acc_length
            totufluxr = f_zero
            totdfluxr = f_zero
            !$acc loop seq
            do ib = 1, nbands
               totufluxr = totufluxr + toturad(n, k,ib)
               totdfluxr = totdfluxr + totdrad(n, k,ib)
            enddo
            totuflux(n, k) = totufluxr * flxfac
            totdflux(n, k) = totdfluxr * flxfac
         end do
      end do
      
      !$acc parallel loop private(fnet, fnet1, rfdelp) async(async_id)
      do n = 1, nxjp_acc_length
         !  --- ...  calculate net fluxes and heating rates
         fnet1 = totuflux(n, 0) - totdflux(n, 0)
         !$acc loop seq
         do k = 1, nlay
            rfdelp = heatfac / delp(n, k)
         enddo
      end do
      
      !$acc parallel loop private(fnet, fnet1, rfdelp, fnet4, fnet3, iplon, jj) async(async_id)
      do n = 1, nxjp_acc_length
         iplon = map_i(n)
         jj = map_jj(n) - jjoffset
         !! --- ...  optional clear sky heating rates
         if ( lhlw0 )fnet1 = totuclfl(n, 0) - totdclfl(n, 0)
         fnet4 = totuflux(n, 0) - totdflux(n, 0)
         !$acc loop seq
         do k = 1, nlay
            rfdelp = heatfac / delp(n, k)
            fnet = totuclfl(n, k) - totdclfl(n, k)
               
            fnet3 = totuflux(n, k) - totdflux(n, k)
            if (ivflip == 0) then
               k1 = nlp1 - k
               hlwc(iplon, k1-kd, jj) = (fnet4 - fnet3) * rfdelp
               if ( lhlw0 ) hlw0(iplon,k1-kd, jj) = (fnet1 - fnet) * rfdelp
            else
               hlwc(iplon, k-kd, jj) = (fnet4 - fnet3) * rfdelp
               if ( lhlw0 ) hlw0(iplon,k-kd, jj) = (fnet1 - fnet) * rfdelp
            end if
            fnet4 = fnet3
            fnet1 = fnet
         enddo
      end do
         
      !! --- ...  optional spectral band heating rates
      if ( lhlwb ) then
         !$acc parallel loop collapse(2) private(rfdelp, k1, iplon, jj, fnet2r, fnet2r1) async(async_id)
         do ib = 1, nbands
            do n = 1, nxjp_acc_length
               iplon = map_i(n)
               jj = map_jj(n) - jjoffset
               fnet2r = (toturad(n, 0,ib) - totdrad(n, 0,ib)) * flxfac
               !$acc loop seq
               do k = 1, nlay
                  rfdelp = heatfac / delp(n, k)
                  fnet2r1 = (toturad(n, k,ib) - totdrad(n, k,ib)) * flxfac
                  k1 = nlp1 - k
                  if (ivflip == 0) then
                     if (lhtrlwb) hlwb(iplon, k1-kd,ib, jj) = (fnet2r - fnet2r1) * rfdelp
                  else
                     if (lhtrlwb) hlwb(iplon, k,ib, jj) = (fnet2r - fnet2r1) * rfdelp
                  end if
                  fnet2r = fnet2r1
               enddo
            end do
         end do
      end if

      !$acc end data

! .................................
      end subroutine rtrnmr
! ---------------------------------

      subroutine rtrnmr_ngptlw(secdif, small_tautot, small_fracs, &
         totdrad, totdclfl, facclr1d, faccld1d, faccmb1d, &
         faccmb2d, facclr2d, faccld2d, semiss, toturad, totuclfl, cldfrc, &
         facclr1u, faccld1u, faccmb1u, faccmb2u, facclr2u, faccld2u, taucld, &
         nlay, myim, ix, ngptlw, nlp1, nbands, nslw, nglw, ng00, map_jj, map_i, &
         nxjp_acc_length, jjoffset, max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj, &
         gassrcu, totsrcu, &
         trngas, trntot, radtotd_2, radtotu_2, stemp, tavel, tlvl)

      implicit none
      integer :: nlay, nlp1, ix, myim(fulljj), ngptlw, nbands, nslw, nglw, fulljj, &
         blockjj, nxjp_acc_length, jjoffset, max_nxjp_acc_length, jbs_nxjp_acc
      integer, dimension(nxjp_acc_length), intent(in) :: map_jj, map_i
      real (kind=kind_phys), dimension(nxjp_acc_length, nlay, ng00) :: small_fracs, &
     &       small_tautot
      real (kind=kind_phys), dimension(max_nxjp_acc_length), intent(in) :: stemp
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: tavel
      real (kind=kind_phys), dimension(nxptot, nlp1), intent(in) :: tlvl

      real (kind=kind_phys), dimension(max_nxjp_acc_length, 0:nlp1) :: cldfrc
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nbands) :: semiss,   &
     &       secdif
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay, nbands) :: taucld
      real (kind=kind_phys), parameter :: rec_6 = 0.166667

!  ---  locals:


      real (kind=kind_phys),dimension(nxjp_acc_length, nlay, ng00) :: &
         gassrcu, totsrcu, trngas, trntot
      real (kind=kind_phys),dimension(nxjp_acc_length, 0:nlay, ng00) :: &
         radtotd_2, radtotu_2

      real (kind=kind_phys) :: totsrcd, gassrcd, tblind, odepth, odtot, &
     &       odcld, atrtot, atrgas, reflct, totfac, gasfac,     &
     &       plfrac, blay, bbdgas, bbdtot, bbugas, bbutot, dplnku,      &
     &       dplnkd, radtotu, radclru, rad0, rad,     &
     &       totradd, clrradd, totradu, clrradu, &
     &       radmod, clfr, trng, trnt, gasu, totu

      integer :: ittot, itgas, ib, ig, k, iplon, jj, iplon2, i2, jj2, ir, ng00, k1

!  dimensions for cloud overlap adjustment
      real (kind=kind_phys), dimension(nxjp_acc_length, nlp1) :: faccld1u, faccld2u,     &
     &        facclr1u, facclr2u, faccmb1u, faccmb2u
      real (kind=kind_phys), dimension(nxjp_acc_length, 0:nlay) :: faccld1d, faccld2d,   &
     &        facclr1d, facclr2d, faccmb1d, faccmb2d
      real (kind=kind_phys), dimension(max_nxjp_acc_length, 0:nlay) :: totdclfl, totuclfl
      real (kind=kind_phys), dimension(nxjp_acc_length, 0:nlay, nbands) :: toturad, totdrad

      real(kind=kind_phys) :: radtotd, radclrd, pklevr1, pklevr, clfr1, semissr, &
         secdifr, radtotdr, radtotur, flxfac, blayr
      logical :: lstcldr
      integer :: async_id, indlay, indlev, n, rr
      ! GPU: variable name changed: CPU - tz,          GPU - tlvl
      
      ib = ngb(nslw + 1)
      flxfac = wtdiff * fluxfac
      !$acc parallel loop gang collapse(2) private(odepth, atrgas, trng, &
      !$acc&         gasfac, tblind, itgas, plfrac, blay, radtotd, radclrd, &
      !$acc&         dplnku, dplnkd, bbdgas, bbugas, gassrcd, clfr, totradd, &
      !$acc&         clrradd, rad, odcld, odtot, totfac, atrtot, trnt, ittot, &
      !$acc&         bbdtot, bbutot, totsrcd, radmod, reflct, rad0, radtotu, &
      !$acc&         radclru, gasu, totradu, clrradu, totu, ir, pklevr1, pklevr, &
      !$acc&         clfr1, lstcldr, secdifr, semissr, indlay, indlev, blayr, k1, rr) async(async_id)
      do ig = nslw + 1, nslw + nglw
         do n = 1, nxjp_acc_length
            !$acc cache(totplnk)
            rr = n + jbs_nxjp_acc - 1
            k1 = 0
            if (ivflip == 0) k1 = nlay
            ir = ig - nslw
            !  --- ...  loop over all g-points
            radtotd = f_zero
            radclrd = f_zero
            !  --- ...  downward radiative transfer loop.
            
            dplnku = tavel(n, nlay)
            dplnkd = tlvl(rr,nlp1 - k1)
            indlay = min(180, max(1, int(dplnku-159.0) ))
            bbdgas = dplnku - int(dplnku)
            indlev = min(180, max(1, int(dplnkd-159.0) ))
            bbugas = dplnkd - int(dplnkd)
            blay = delwave(ib) * (totplnk(indlay,ib) + bbdgas         &
            &               * (totplnk(indlay+1,ib) - totplnk(indlay,ib)) )
            pklevr = delwave(ib) * (totplnk(indlev,ib) + bbugas         &
            &               * (totplnk(indlev+1,ib) - totplnk(indlev,ib)) )




            clfr = cldfrc(n, nlay)
            lstcldr = clfr > eps
            secdifr = secdif(n, ib)
            !$acc loop seq
            do k = nlay, 1, -1
               k1 = k
               if (ivflip == 0) k1 = nlp1 - k - 1
               if (k .eq. 1) then
                  dplnku = stemp(n)
               else
                  dplnku = tavel(n, k-1)
               endif
               dplnkd = tlvl(rr,k1)
               indlay = min(180, max(1, int(dplnku-159.0) ))
               bbdgas = dplnku - int(dplnku)
               indlev = min(180, max(1, int(dplnkd-159.0) ))
               bbugas = dplnkd - int(dplnkd)
               blayr = delwave(ib) * (totplnk(indlay,ib) + bbdgas         &
               &               * (totplnk(indlay+1,ib) - totplnk(indlay,ib)) )
               pklevr1 = delwave(ib) * (totplnk(indlev,ib) + bbugas         &
               &               * (totplnk(indlev+1,ib) - totplnk(indlev,ib)) )

               !  --- ...  clear sky, gases contribution

               odepth = max( f_zero, secdifr*small_tautot(n, k,ir) )
               if (odepth <= 0.06) then
                  atrgas = odepth - 0.5*odepth*odepth
                  trng   = f_one - atrgas
                  gasfac = rec_6 * odepth
               else
                  tblind = odepth / (bpade + odepth)
                  itgas = tblint*tblind + 0.5
                  trng  = exp_tbl(itgas)
                  atrgas = f_one - trng
                  gasfac = tfn_tbl(itgas)
                  odepth = tau_tbl(itgas)
               endif

               plfrac = small_fracs(n, k,ir)

               dplnku = pklevr - blay
               dplnkd = pklevr1 - blay
               bbdgas = plfrac * (blay + dplnkd*gasfac)
               bbugas = plfrac * (blay + dplnku*gasfac)
               gassrcd   = bbdgas * atrgas
               gassrcu(n, k, ir)= bbugas * atrgas
               trngas(n, k, ir) = trng
               pklevr = pklevr1

               !  --- ...  total sky, gases+clouds contribution

               if (lstcldr) then
                  totradd = clfr * radtotd
                  clrradd = radtotd - totradd
                  rad = f_zero
               endif
               clfr1 = cldfrc(n, k-1)
               lstcldr = clfr1 > eps .and. clfr<=eps


               if (clfr >= eps) then
                  !  --- ...  cloudy layer

                  odcld = secdifr * taucld(n, k,ib)
                  odtot = odepth + odcld
                  if (odtot < 0.06) then
                     totfac = rec_6 * odtot
                     atrtot = odtot - 0.5*odtot*odtot
                     trnt   = f_one - atrtot
                  else
                     tblind = odtot / (bpade + odtot)
                     ittot  = tblint*tblind + 0.5
                     totfac = tfn_tbl(ittot)
                     trnt   = exp_tbl(ittot)
                     atrtot = f_one - trnt
                  endif

                  bbdtot = plfrac * (blay + dplnkd*totfac)
                  bbutot = plfrac * (blay + dplnku*totfac)
                  totsrcd   = bbdtot * atrtot
                  totsrcu(n, k, ir)= bbutot * atrtot
                  trntot(n, k, ir) = trnt

                  totradd = totradd*trnt + clfr*totsrcd
                  clrradd = clrradd*trng + (f_one - clfr)*gassrcd

                  !  --- ...  total sky radiance
                  radtotd = totradd + clrradd
                  radtotd_2(n, k, ir) = radtotd

                  !  --- ...  clear sky radiance
                  radclrd = radclrd*trng + gassrcd
                  !$acc atomic
                  totdclfl(n, k-1) = totdclfl(n, k-1) + radclrd * flxfac

                  radmod = rad*(facclr1d(n, k-1)*trng + faccld1d(n, k-1)*trnt)      &
                  &             - faccmb1d(n, k-1)*gassrcd + faccmb2d(n, k-1)*totsrcd

                  rad = -radmod + facclr2d(n, k-1)*(clrradd + radmod)            &
                  &                    - faccld2d(n, k-1)*(totradd - radmod)
                  totradd = totradd + rad
                  clrradd = clrradd - rad

               else
                  !  --- ...  clear layer

                  !  --- ...  total sky radiance
                  radtotd = radtotd*trng + gassrcd
                  radtotd_2(n, k, ir) = radtotd

                  !  --- ...  clear sky radiance
                  radclrd = radclrd*trng + gassrcd
                  !$acc atomic
                  totdclfl(n, k-1) = totdclfl(n, k-1) + radclrd * flxfac

               endif   ! end if_clfr_block
               clfr = clfr1
               blay = blayr

            enddo   ! end do_k_loop

            !  --- ...  spectral emissivity & reflectance
            !           include the contribution of spectrally varying longwave emissivity
            !           and reflection from the surface to the upward radiative transfer.
            !     note: spectral and lambertian reflection are identical for the
            !           diffusivity angle flux integration used here.
            semissr = semiss(n, ib)
            reflct = f_one - semissr
            rad0 = semissr * small_fracs(n, 1,ir) * pklevr1

            !  --- ...  total sky radiance
            radtotu = rad0 + reflct*radtotd
            radtotu_2(n, 0, ir) = radtotu

            !  --- ...  clear sky radiance
            radclru = rad0 + reflct*radclrd
            !$acc atomic
            totuclfl(n, 0) = totuclfl(n, 0) + radclru * flxfac

            !  --- ...  upward radiative transfer loop.
            clfr = cldfrc(n, 1)
            lstcldr = clfr > eps
            !$acc loop seq
            do k = 1, nlay

               trng = trngas(n, k, ir)
               gasu = gassrcu(n, k, ir)

               if (lstcldr) then
                  totradu = clfr * radtotu
                  clrradu = radtotu - totradu
                  rad = f_zero
               endif
               clfr1 = cldfrc(n, k+1)
               lstcldr = clfr1>eps .and. clfr<=eps


               if (clfr >= eps) then
                  !  --- ...  cloudy layer

                  trnt = trntot(n, k, ir)
                  totu = totsrcu(n, k, ir)
                  totradu = totradu*trnt + clfr*totu
                  clrradu = clrradu*trng + (f_one - clfr)*gasu

                  !  --- ...  total sky radiance
                  radtotu = totradu + clrradu
                  radtotu_2(n, k, ir) = radtotu

                  !  --- ...  clear sky radiance
                  radclru = radclru*trng + gasu
                  !$acc atomic
                  totuclfl(n, k) = totuclfl(n, k) + radclru * flxfac

                  radmod = rad*(facclr1u(n, k+1)*trng + faccld1u(n, k+1)*trnt)      &
                  &             - faccmb1u(n, k+1)*gasu + faccmb2u(n, k+1)*totu
                  rad = -radmod + facclr2u(n, k+1)*(clrradu + radmod)            &
                  &                    - faccld2u(n, k+1)*(totradu - radmod)
                  totradu = totradu + rad
                  clrradu = clrradu - rad

               else
                  !  --- ...  clear layer

                  !  --- ...  total sky radiance
                  radtotu = radtotu*trng + gasu
                  radtotu_2(n, k, ir) = radtotu

                  !  --- ...  clear sky radiance
                  radclru = radclru*trng + gasu
                  !$acc atomic
                  totuclfl(n, k) = totuclfl(n, k) + radclru * flxfac

               endif   ! end if_clfr_block
               clfr = clfr1
            enddo   ! end do_k_loop
         end do
      end do

      !$acc parallel loop collapse(2) private(ir, radtotur, radtotdr) async(async_id)
      do k = 0, nlay
         do n = 1, nxjp_acc_length
            radtotdr = 0.
            radtotur = 0.
            !$acc loop seq
            do ig = nslw + 1, nslw + nglw
               ir = ig - nslw
               if (k .ge. 1) radtotdr = radtotdr + radtotd_2(n, k, ir)
               radtotur = radtotur + radtotu_2(n, k, ir)
            end do
            toturad(n, k,ib) =  radtotur
            if (k .ge. 1) totdrad(n, k-1, ib) = radtotdr
         end do
      end do




      return
      end subroutine


! ---------------------------------
      subroutine rtrnmc                                                 &
! .................................
!  ---  inputs:
     &     ( semiss,delp,cldfmc,taucld,tautot,pklay,pklev,              &
     &       fracs,secdif, nlay,nlp1,                                   &
!  ---  outputs:
     &       totuflux,totdflux,htr, totuclfl,totdclfl,htrcl, htrb       &
     &     )

!  ===================  program usage description  ===================  !
!                                                                       !
! purpose:  compute the upward/downward radiative fluxes, and heating   !
! rates for both clear or cloudy atmosphere.  clouds are treated with   !
! the mcica stochastic approach.                                        !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                     -size-   !
!   semiss  - real, lw surface emissivity                         nbands!
!   delp    - real, layer pressure thickness (mb)                  nlay !
!   cldfmc  - real, layer cloud fraction (sub-column)        ngptlw*nlay!
!   taucld  - real, layer cloud opt depth                    nbands*nlay!
!   tautot  - real, total optical depth (gas+aerosols)       ngptlw*nlay!
!   pklay   - real, integrated planck func at lay temp     nbands*0:nlay!
!   pklev   - real, integrated planck func at lev temp     nbands*0:nlay!
!   fracs   - real, planck fractions                         ngptlw*nlay!
!   secdif  - real, secant of diffusivity angle                   nbands!
!   nlay    - integer, number of vertical layers                    1   !
!   nlp1    - integer, number of vertical levels (interfaces)       1   !
!                                                                       !
!  outputs:                                                             !
!   totuflux- real, total sky upward flux (w/m2)                 0:nlay !
!   totdflux- real, total sky downward flux (w/m2)               0:nlay !
!   htr     - real, total sky heating rate (k/sec or k/day)        nlay !
!   totuclfl- real, clear sky upward flux (w/m2)                 0:nlay !
!   totdclfl- real, clear sky downward flux (w/m2)               0:nlay !
!   htrcl   - real, clear sky heating rate (k/sec or k/day)        nlay !
!   htrb    - real, spectral band lw heating rate (k/day)    nlay*nbands!
!                                                                       !
!  module veriables:                                                    !
!   ngb     - integer, band index for each g-value                ngptlw!
!   fluxfac - real, conversion factor for fluxes (pi*2.e4)           1  !
!   heatfac - real, conversion factor for heating rates (g/cp*1e-2)  1  !
!   tblint  - real, conversion factor for look-up tbl (float(ntbl)   1  !
!   bpade   - real, pade approx constant (1/0.278)                   1  !
!   wtdiff  - real, weight for radiance to flux conversion           1  !
!   ntbl    - integer, dimension of look-up tables                   1  !
!   tau_tbl - real, clr-sky opt dep lookup table                 0:ntbl !
!   exp_tbl - real, transmittance lookup table                   0:ntbl !
!   tfn_tbl - real, tau transition function                      0:ntbl !
!                                                                       !
!  local variables:                                                     !
!    itgas  - integer, index for gases contribution look-up table    1  !
!    ittot  - integer, index for gases plus clouds  look-up table    1  !
!    reflct - real, surface reflectance                              1  !
!    atrgas - real, gaseous absorptivity                             1  !
!    atrtot - real, gaseous and cloud absorptivity                   1  !
!    odcld  - real, cloud optical depth                              1  !
!    efclrfr- real, effective clear sky fraction (1-efcldfr)        nlay!
!    odepth - real, optical depth of gaseous only                    1  !
!    odtot  - real, optical depth of gas and cloud                   1  !
!    gasfac - real, gas-only pade factor, used for planck function   1  !
!    totfac - real, gas and cloud pade factor, used for planck fn    1  !
!    bbdgas - real, gas-only planck function for downward rt         1  !
!    bbugas - real, gas-only planck function for upward rt           1  !
!    bbdtot - real, gas and cloud planck function for downward rt    1  !
!    bbutot - real, gas and cloud planck function for upward rt      1  !
!    gassrcu- real, upwd source radiance due to gas                 nlay!
!    totsrcu- real, upwd source radiance due to gas+cld             nlay!
!    gassrcd- real, dnwd source radiance due to gas                  1  !
!    totsrcd- real, dnwd source radiance due to gas+cld              1  !
!    radtotu- real, spectrally summed total sky upwd radiance        1  !
!    radclru- real, spectrally summed clear sky upwd radiance        1  !
!    radtotd- real, spectrally summed total sky dnwd radiance        1  !
!    radclrd- real, spectrally summed clear sky dnwd radiance        1  !
!    toturad- real, total sky upward radiance by layer     0:nlay*nbands!
!    clrurad- real, clear sky upward radiance by layer     0:nlay*nbands!
!    totdrad- real, total sky downward radiance by layer   0:nlay*nbands!
!    clrdrad- real, clear sky downward radiance by layer   0:nlay*nbands!
!    fnet   - real, net longwave flux (w/m2)                     0:nlay !
!    fnetc  - real, clear sky net longwave flux (w/m2)           0:nlay !
!                                                                       !
!                                                                       !
!  *******************************************************************  !
!  original code description                                            !
!                                                                       !
!  original version:   e. j. mlawer, et al. rrtm_v3.0                   !
!  revision for gcms:  michael j. iacono; october, 2002                 !
!  revision for f90:   michael j. iacono; june, 2006                    !
!                                                                       !
!  this program calculates the upward fluxes, downward fluxes, and      !
!  heating rates for an arbitrary clear or cloudy atmosphere. the input !
!  to this program is the atmospheric profile, all planck function      !
!  information, and the cloud fraction by layer.  a variable diffusivity!
!  angle (secdif) is used for the angle integration. bands 2-3 and 5-9  !
!  use a value for secdif that varies from 1.50 to 1.80 as a function   !
!  of the column water vapor, and other bands use a value of 1.66.  the !
!  gaussian weight appropriate to this angle (wtdiff=0.5) is applied    !
!  here.  note that use of the emissivity angle for the flux integration!
!  can cause errors of 1 to 4 w/m2 within cloudy layers.                !
!  clouds are treated with the mcica stochastic approach and            !
!  maximum-random cloud overlap.                                        !
!                                                                       !
!  *******************************************************************  !
!  ======================  end of description block  =================  !

!  ---  inputs:
      integer, intent(in) :: nlay, nlp1

      real (kind=kind_phys), dimension(nbands), intent(in) :: semiss,   &
     &       secdif
      real (kind=kind_phys), dimension(nlay),   intent(in) :: delp

      real (kind=kind_phys), dimension(nlay, nbands),intent(in):: taucld
      real (kind=kind_phys), dimension(nlay, ngptlw),intent(in):: fracs, &
     &       tautot, cldfmc

      real (kind=kind_phys), dimension(0:nlay, nbands), intent(in) ::    &
     &       pklev, pklay

!  ---  outputs:
      real (kind=kind_phys), dimension(nlay), intent(out) :: htr, htrcl

      real (kind=kind_phys), dimension(nlay,nbands),intent(out) :: htrb

      real (kind=kind_phys), dimension(0:nlay), intent(out) ::          &
     &       totuflux, totdflux, totuclfl, totdclfl

!  ---  locals:
      real (kind=kind_phys), parameter :: rec_6 = 0.166667

      real (kind=kind_phys), dimension(0:nlay,nbands) :: clrurad,       &
     &       clrdrad, toturad, totdrad

      real (kind=kind_phys), dimension(nlay)   :: gassrcu, totsrcu,     &
     &       trngas, efclrfr, rfdelp
      real (kind=kind_phys), dimension(0:nlay) :: fnet, fnetc

      real (kind=kind_phys) :: totsrcd, gassrcd, tblind, odepth, odtot, &
     &       odcld, atrtot, atrgas, reflct, totfac, gasfac, flxfac,     &
     &       plfrac, blay, bbdgas, bbdtot, bbugas, bbutot, dplnku,      &
     &       dplnkd, radtotu, radclru, radtotd, radclrd, rad0,          &
     &       clfm, trng, gasu

      integer :: ittot, itgas, ib, ig, k
!
!===> ...  begin here
!
      do ib = 1, nbands
        do k = 0, nlay
          toturad(k,ib) = f_zero
          totdrad(k,ib) = f_zero
          clrurad(k,ib) = f_zero
          clrdrad(k,ib) = f_zero
        enddo
      enddo

      do k = 0, nlay
        totuflux(k) = f_zero
        totdflux(k) = f_zero
        totuclfl(k) = f_zero
        totdclfl(k) = f_zero
      enddo

!  --- ...  loop over all g-points

      do ig = 1, ngptlw
        ib = ngb(ig)

        radtotd = f_zero
        radclrd = f_zero

!  --- ...  downward radiative transfer loop.

        do k = nlay, 1, -1

!  --- ...  clear sky, gases contribution

          odepth = max( f_zero, secdif(ib)*tautot(k, ig) )
          if (odepth <= 0.06) then
            atrgas = odepth - 0.5*odepth*odepth
            trng   = f_one - atrgas
            gasfac = rec_6 * odepth
          else
            tblind = odepth / (bpade + odepth)
            itgas = tblint*tblind + 0.5
            trng  = exp_tbl(itgas)
            atrgas = f_one - trng
            gasfac = tfn_tbl(itgas)
            odepth = tau_tbl(itgas)
          endif

          plfrac = fracs(k, ig)
          blay = pklay(k, ib)

          dplnku = pklev(k, ib  ) - blay
          dplnkd = pklev(k-1, ib) - blay
          bbdgas = plfrac * (blay + dplnkd*gasfac)
          bbugas = plfrac * (blay + dplnku*gasfac)
          gassrcd= bbdgas * atrgas
          gassrcu(k)= bbugas * atrgas
          trngas(k) = trng

!  --- ...  total sky, gases+clouds contribution

          clfm = cldfmc(k, ig)
          if (clfm >= eps) then
!  --- ...  cloudy layer

            odcld = secdif(ib) * taucld(k, ib)
            efclrfr(k) = f_one - (f_one - exp(-odcld))*clfm
            odtot = odepth + odcld
            if (odtot < 0.06) then
              totfac = rec_6 * odtot
              atrtot = odtot - 0.5*odtot*odtot
            else
              tblind = odtot / (bpade + odtot)
              ittot  = tblint*tblind + 0.5
              totfac = tfn_tbl(ittot)
              atrtot = f_one - exp_tbl(ittot)
            endif

            bbdtot = plfrac * (blay + dplnkd*totfac)
            bbutot = plfrac * (blay + dplnku*totfac)
            totsrcd= bbdtot * atrtot
            totsrcu(k)= bbutot * atrtot

!  --- ...  total sky radiance
            radtotd = radtotd*trng*efclrfr(k) + gassrcd                 &
     &              + clfm*(totsrcd - gassrcd)
            totdrad(k-1,ib) = totdrad(k-1,ib) + radtotd

!  --- ...  clear sky radiance
            radclrd = radclrd*trng + gassrcd
            clrdrad(k-1,ib) = clrdrad(k-1,ib) + radclrd

          else
!  --- ...  clear layer

!  --- ...  total sky radiance
            radtotd = radtotd*trng + gassrcd
            totdrad(k-1,ib) = totdrad(k-1,ib) + radtotd

!  --- ...  clear sky radiance
            radclrd = radclrd*trng + gassrcd
            clrdrad(k-1,ib) = clrdrad(k-1,ib) + radclrd

          endif   ! end if_clfm_block

        enddo   ! end do_k_loop

!  --- ...  spectral emissivity & reflectance
!           include the contribution of spectrally varying longwave emissivity
!           and reflection from the surface to the upward radiative transfer.
!     note: spectral and lambertian reflection are identical for the
!           diffusivity angle flux integration used here.

        reflct = f_one - semiss(ib)
        rad0 = semiss(ib) * fracs(1,ig) * pklay(0,ib)

!  --- ... total sky radiance
        radtotu = rad0 + reflct*radtotd
        toturad(0,ib) = toturad(0,ib) + radtotu

!  --- ... clear sky radiance
        radclru = rad0 + reflct*radclrd
        clrurad(0,ib) = clrurad(0,ib) + radclru

!  --- ...  upward radiative transfer
!          toturad holds summed radiance for total sky stream
!          clrurad holds summed radiance for clear sky stream

        do k = 1, nlay
          clfm = cldfmc(k,ig)
          trng = trngas(k)
          gasu = gassrcu(k)

          if (clfm > eps) then
!  --- ...  cloudy layer

!  --- ... total sky radiance
            radtotu = radtotu*trng*efclrfr(k) + gasu                    &
     &              + clfm*(totsrcu(k) - gasu)
            toturad(k,ib) = toturad(k,ib) + radtotu

!  --- ... clear sky radiance
            radclru = radclru*trng + gasu
            clrurad(k,ib) = clrurad(k,ib) + radclru

          else
!  --- ...  clear layer

!  --- ... total sky radiance
            radtotu = radtotu*trng + gasu
            toturad(k,ib) = toturad(k,ib) + radtotu

!  --- ... clear sky radiance
            radclru = radclru*trng + gasu
            clrurad(k,ib) = clrurad(k,ib) + radclru

          endif   ! end if_clfm_block

        enddo   ! end do_k_loop

      enddo   ! end do_ig_loop

!  --- ...  process longwave output from band for total and clear streams.
!           calculate upward, downward, and net flux.

      flxfac = wtdiff * fluxfac

      do k = 0, nlay
        do ib = 1, nbands
          totuflux(k) = totuflux(k) + toturad(k,ib)
          totdflux(k) = totdflux(k) + totdrad(k,ib)
          totuclfl(k) = totuclfl(k) + clrurad(k,ib)
          totdclfl(k) = totdclfl(k) + clrdrad(k,ib)
        enddo

        totuflux(k) = totuflux(k) * flxfac
        totdflux(k) = totdflux(k) * flxfac
        totuclfl(k) = totuclfl(k) * flxfac
        totdclfl(k) = totdclfl(k) * flxfac
      enddo

!  --- ...  calculate net fluxes and heating rates
      fnet(0) = totuflux(0) - totdflux(0)

      do k = 1, nlay
        rfdelp(k) = heatfac / delp(k)
        fnet(k) = totuflux(k) - totdflux(k)
        htr (k) = (fnet(k-1) - fnet(k)) * rfdelp(k)
      enddo

!! --- ...  optional clear sky heating rates
      if ( lhlw0 ) then
        fnetc(0) = totuclfl(0) - totdclfl(0)

        do k = 1, nlay
          fnetc(k) = totuclfl(k) - totdclfl(k)
          htrcl(k) = (fnetc(k-1) - fnetc(k)) * rfdelp(k)
        enddo
      endif

!! --- ...  optional spectral band heating rates
      if ( lhlwb ) then
        do ib = 1, nbands
          fnet(0) = (toturad(0,ib) - totdrad(0,ib)) * flxfac

          do k = 1, nlay
            fnet(k) = (toturad(k,ib) - totdrad(k,ib)) * flxfac
            htrb(k,ib) = (fnet(k-1) - fnet(k)) * rfdelp(k)
          enddo
        enddo
      endif

! ..................................
      end subroutine rtrnmc
! ----------------------------------


! ----------------------------------
      subroutine taumol                                                 &
! ..................................
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj, ib, ng00,                                                         &
!  ---  outputs:
     &       small_fracs, small_tautot                                              &
     &     )

!  ************    original subprogram description    ***************   !
!                                                                       !
!                  optical depths developed for the                     !
!                                                                       !
!                rapid radiative transfer model (rrtm)                  !
!                                                                       !
!            atmospheric and environmental research, inc.               !
!                        131 hartwell avenue                            !
!                        lexington, ma 02421                            !
!                                                                       !
!                           eli j. mlawer                               !
!                         jennifer delamere                             !
!                         steven j. taubman                             !
!                         shepard a. clough                             !
!                                                                       !
!                       email:  mlawer@aer.com                          !
!                       email:  jdelamer@aer.com                        !
!                                                                       !
!        the authors wish to acknowledge the contributions of the       !
!        following people:  karen cady-pereira, patrick d. brown,       !
!        michael j. iacono, ronald e. farren, luke chen,                !
!        robert bergstrom.                                              !
!                                                                       !
!  revision for g-point reduction: michael j. iacono; aer, inc.         !
!                                                                       !
!     taumol                                                            !
!                                                                       !
!     this file contains the subroutines taugbn (where n goes from      !
!     1 to 16).  taugbn calculates the optical depths and planck        !
!     fractions per g-value and layer for band n.                       !
!                                                                       !
!  *******************************************************************  !
!  ==================   program usage description   ==================  !
!                                                                       !
!    call  taumol                                                       !
!       inputs:                                                         !
!          ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              !
!            rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  !
!            selffac,selffrac,indself,forfac,forfrac,indfor,            !
!            minorfrac,scaleminor,scaleminorn2,indminor,                !
!            nlay,                                                      !
!       outputs:                                                        !
!            fracs, tautot )                                            !
!                                                                       !
!  subprograms called:  taugb## (## = 01 -16)                           !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                        size  !
!     laytrop   - integer, tropopause layer index (unitless)        1   !
!                   layer at which switch is made for key species       !
!     pavel     - real, layer pressures (mb)                       nlay !
!     coldry    - real, column amount for dry air (mol/cm2)        nlay !
!     colamt    - real, column amounts of h2o, co2, o3, n2o, ch4,       !
!                   o2, co (mol/cm**2)                       nlay*maxgas!
!     colbrd    - real, column amount of broadening gases          nlay !
!     wx        - real, cross-section amounts(mol/cm2)      nlay*maxxsec!
!     tauaer    - real, aerosol optical depth               nbands*nlay !
!     rfrate    - real, reference ratios of binary species parameter    !
!     (:,m,:)m=1-h2o/co2,2-h2o/o3,3-h2o/n2o,4-h2o/ch4,5-n2o/co2,6-o3/co2!
!     (:,:,n)n=1,2: the rates of ref press at the 2 sides of the layer  !
!                                                          nlay*nrates*2!
!     facij     - real, factors multiply the reference ks, i,j of 0/1   !
!                   for lower/higher of the 2 appropriate temperatures  !
!                   and altitudes                                  nlay !
!     jp        - real, index of lower reference pressure          nlay !
!     jt, jt1   - real, indices of lower reference temperatures    nlay !
!                   for pressure levels jp and jp+1, respectively       !
!     selffac   - real, scale factor for water vapor self-continuum     !
!                   equals (water vapor density)/(atmospheric density   !
!                   at 296k and 1013 mb)                           nlay !
!     selffrac  - real, factor for temperature interpolation of         !
!                   reference water vapor self-continuum data      nlay !
!     indself   - integer, index of lower reference temperature for     !
!                   the self-continuum interpolation               nlay !
!     forfac    - real, scale factor for w. v. foreign-continuum   nlay !
!     forfrac   - real, factor for temperature interpolation of         !
!                   reference w.v. foreign-continuum data          nlay !
!     indfor    - integer, index of lower reference temperature for     !
!                   the foreign-continuum interpolation            nlay !
!     minorfrac - real, factor for minor gases                     nlay !
!     scaleminor,scaleminorn2                                           !
!               - real, scale factors for minor gases              nlay !
!     indminor  - integer, index of lower reference temperature for     !
!                   minor gases                                    nlay !
!     nlay      - integer, total number of layers                   1   !
!                                                                       !
!  outputs:                                                             !
!     fracs     - real, planck fractions                     ngptlw,nlay!
!     tautot    - real, total optical depth (gas+aerosols)   ngptlw,nlay!
!                                                                       !
!  internal variables:                                                  !
!     ng##      - integer, number of g-values in band ## (##=01-16) 1   !
!     nspa      - integer, for lower atmosphere, the number of ref      !
!                   atmos, each has different relative amounts of the   !
!                   key species for the band                      nbands!
!     nspb      - integer, same but for upper atmosphere          nbands!
!     absa      - real, k-values for lower ref atmospheres (no w.v.     !
!                   self-continuum) (cm**2/molecule)  nspa(##)*5*13*ng##!
!     absb      - real, k-values for high ref atmospheres (all sources) !
!                   (cm**2/molecule)               nspb(##)*5*13:59*ng##!
!     ka_m'mgas'- real, k-values for low ref atmospheres minor species  !
!                   (cm**2/molecule)                          mmn##*ng##!
!     kb_m'mgas'- real, k-values for high ref atmospheres minor species !
!                   (cm**2/molecule)                          mmn##*ng##!
!     selfref   - real, k-values for w.v. self-continuum for ref atmos  !
!                   used below laytrop (cm**2/mol)               10*ng##!
!     forref    - real, k-values for w.v. foreign-continuum for ref atmos
!                   used below/above laytrop (cm**2/mol)          4*ng##!
!                                                                       !
!  ******************************************************************   !

!  ---  inputs:
      integer, intent(in) :: nlay, laytrop(max_nxjp_acc_length), ix, myim(fulljj), &
         async_id, fulljj, blockjj, ib, ng00, nxjp_acc_length, jjoffset, &
         max_nxjp_acc_length, jbs_nxjp_acc
      integer, dimension(nxjp_acc_length), intent(in) :: map_jj, map_i

      integer, dimension(max_nxjp_acc_length, nlay), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: h2ovmr, o3vmr
      real (kind=kind_phys), dimension(nxptot,nlay), intent(in) :: gasvmr_co2
      real (kind=kind_phys), dimension(nf_vgas), intent(in) :: gasvmr_other
      integer, intent(in) :: nf_vgas


      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay, nbands), intent(in):: tauaer


!  ---  outputs:
      real (kind=kind_phys), dimension(nxjp_acc_length, nlay, ng00) :: small_fracs, &
     &       small_tautot

     ! GPU: In order to reduce the GPU memory usage, both of variables fracs and tautot 
     ! GPU: are replaced by variables small_fracs and small_tautot, respectively.
     ! GPU: The length of third dimension of those arrays is reduce from ngptlw = 140
     ! GPU: to ng00 = 16, which is the maximum length of ng01 ~ ng16.

     ! GPU: In this GPU version, taumol need to be called in the ib-loop of subroutine  
     ! GPU: rtrnmr, instead of being called in subroutine lwrad_gpu originally.
     ! GPU: The fracs and tautot are still remained in this subroutine API. If 
     ! GPU: subroutine ntrnmc or ntrn is needed in the future use and variables 
     ! GPU: fracs and tautot needed to compute at once, do the following steps:
     ! GPU: 1. uncomment taumol caller in lwrad_gpu.
     ! GPU: 2. uncomment frac and tautot in subroutines taugb*
     ! GPU: 3. comment small_frac and small_tautot in subroutines taugb*

      ! GPU: variable name changed: CPU - frac,          GPU - small_frac
      ! GPU: variable name changed: CPU - tautot,        GPU - small_tautot

!  ---  locals
      integer :: small_factor = 1
      integer :: small_ix
      real (kind=kind_phys) :: taug 
      !integer, dimension(:,:,:,:), allocatable :: mask

      integer :: ig, k, jj, i2


      small_ix = ceiling(float(ix)/float(small_factor))
      !allocate(mask(small_ix, nlay, ngptlw, fulljj))
      !mask = 0
!
!===> ...  begin here
!
      select case (ib)
      case (1)
         call taugb01 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )

     case (2)
         call taugb02 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )

     case (3)
         call taugb03 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )

     case (4)
         call taugb04 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )

     case (5)
         call taugb05 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )

     case (6)
         call taugb06 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )

     case (7)
         call taugb07 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )

     case (8)
         call taugb08 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )

     case (9)
         call taugb09 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )

     case (10)
         call taugb10 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )

     case (11)
         call taugb11 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )

     case (12)
         call taugb12 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )

     case (13)
         call taugb13 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )

     case (14)
         call taugb14 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )

     case (15)
         call taugb15 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )

     case (16)
         call taugb16 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )
     end select


! ..................................
      end subroutine taumol
!-----------------------------------

! ----------------------------------
      subroutine taugb01 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!  written by eli j. mlawer, atmospheric & environmental research.     !
!  revised by michael j. iacono, atmospheric & environmental research. !
!                                                                      !
!     band 1:  10-350 cm-1 (low key - h2o; low minor - n2)             !
!                          (high key - h2o; high minor - n2)           !
!                                                                      !
!  compute the optical depth by interpolating in ln(pressure) and      !
!  temperature.  below laytrop, the water vapor self-continuum and     !
!  foreign continuum is interpolated (in temperature) separately.      !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb01
      integer, intent(in) :: nlay, laytrop(nxjp_acc_length), ix, myim(fulljj), &
         async_id, ng00, i2, blockjj, fulljj, nxjp_acc_length, jjoffset, &
         max_nxjp_acc_length, jbs_nxjp_acc
      integer, dimension(nxjp_acc_length), intent(in) :: map_jj, map_i

      integer, dimension(max_nxjp_acc_length, nlay), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys) :: temcol, colamt1, colamt2, colamt3, &
                               colamt4, colamt5, colamt6, colamt7
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: h2ovmr, o3vmr
      real (kind=kind_phys), dimension(nxptot,nlay), intent(in) :: gasvmr_co2
      real (kind=kind_phys), dimension(nf_vgas), intent(in) :: gasvmr_other
      integer, intent(in) :: nf_vgas
      integer :: k1, n, rr

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay, nbands), intent(in):: tauaer

      real (kind=kind_phys) :: rfrate
!  ---  outputs:
      real (kind=kind_phys), dimension(nxjp_acc_length, nlay, ng00) :: small_fracs, &
     &       small_tautot


!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &           indm, indmp, ig, ib

      real (kind=kind_phys) :: pp, corradj, scalen2, tauself, taufor,   &
     &       taun2, taug
      integer :: jj, iplon, iplon2
      ! GPU: variable name changed: CPU - colamt(:,:,1),       GPU - colamt1
!
!===> ...  begin here
!
!  ---  minor gas mapping levels:
!     lower - n2, p = 142.5490 mbar, t = 215.70 k
!     upper - n2, p = 142.5490 mbar, t = 215.70 k

      !$acc parallel loop gang collapse(2) private(ind0, ind1, inds, indf, indm, &
      !$acc&         ind0p, ind1p, indsp, indfp, indmp, pp, scalen2, corradj, &
      !$acc&         tauself, taufor, taun2, taug, colamt1, rr) async(async_id)
      do k = 1, nlay
         do n = 1, nxjp_acc_length
            rr = n + jbs_nxjp_acc - 1
            colamt1 = max(f_zero,    coldry(n, k)*h2ovmr(n, k))          ! h2o
            if (k .le. laytrop(n)) then
               !  --- ...  lower atmosphere loop
               ind0 = ((jp(n, k)-1)*5 + (jt (n, k)-1)) * nspa(1) + 1
               ind1 = ( jp(n, k)   *5 + (jt1(n, k)-1)) * nspa(1) + 1
               inds = indself(n, k)
               indf = indfor(n, k)
               indm = indminor(n, k)

               ind0p = ind0 + 1
               ind1p = ind1 + 1
               indsp = inds + 1
               indfp = indf + 1
               indmp = indm + 1

               pp = pavel(n, k)
               scalen2 = colbrd(n, k) * scaleminorn2(n, k)
               if (pp < 250.0) then
                  corradj = f_one - 0.15 * (250.0-pp) / 154.4
               else
                  corradj = f_one
               endif
               !$acc loop seq
               do ig = 1, ng01
                  ib = ngb(ig)
                  tauself = selffac(n, k) * (selfref(ig,inds) + selffrac(n, k)        &
                  &            * (selfref(ig,indsp) - selfref(ig,inds)))
                  taufor  = forfac(n, k) * (forref(ig,indf) + forfrac(n, k)           &
                  &            * (forref(ig,indfp) -  forref(ig,indf))) 
                  taun2   = scalen2 * (ka_mn2(ig,indm) + minorfrac(n, k)           &
                  &            * (ka_mn2(ig,indmp) - ka_mn2(ig,indm)))

                  taug = corradj * (colamt1                           &
                  &            * (fac00(n, k)*absa(ig,ind0) + fac10(n, k)*absa(ig,ind0p)   &
                  &            +  fac01(n, k)*absa(ig,ind1) + fac11(n, k)*absa(ig,ind1p))  &
                  &            + tauself + taufor + taun2)

                  !fracs(iplon, k, ns01+ig, jj) = fracrefa(ig)
                  !tautot(iplon, k, ns01+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) = fracrefa(ig)
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)
               enddo
            else
               !  --- ...  upper atmosphere loop
               ind0 = ((jp(n, k)-13)*5 + (jt (n, k)-1)) * nspb(1) + 1
               ind1 = ((jp(n, k)-12)*5 + (jt1(n, k)-1)) * nspb(1) + 1
               indf = indfor(n, k)
               indm = indminor(n, k)

               ind0p = ind0 + 1
               ind1p = ind1 + 1
               indfp = indf + 1
               indmp = indm + 1

               scalen2 = colbrd(n, k) * scaleminorn2(n, k)
               corradj = f_one - 0.15 * (pavel(n, k) / 95.6)
               !$acc loop seq
               do ig = 1, ng01
                  ib = ngb(ig)
                  taufor = forfac(n, k) * (forref(ig,indf) + forfrac(n, k)            &
                  &           * (forref(ig,indfp) - forref(ig,indf))) 
                  taun2  = scalen2 * (kb_mn2(ig,indm) + minorfrac(n, k)            &
                  &           * (kb_mn2(ig,indmp) - kb_mn2(ig,indm)))

                  taug = corradj * (colamt1                           &
                  &           * (fac00(n, k)*absb(ig,ind0) + fac10(n, k)*absb(ig,ind0p)    &
                  &           +  fac01(n, k)*absb(ig,ind1) + fac11(n, k)*absb(ig,ind1p))   &
                  &           + taufor + taun2)

                  !fracs(iplon, k, ns01+ig, jj) = fracrefb(ig)
                  !tautot(iplon, k, ns01+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) = fracrefb(ig)
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)
               end do
            end if
         end do
      end do


! ..................................
      end subroutine taugb01
! ----------------------------------

! ----------------------------------
      subroutine taugb02 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!     band 2:  350-500 cm-1 (low key - h2o; high key - h2o)            !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb02
      integer, intent(in) :: nlay, laytrop(nxjp_acc_length), ix, myim(fulljj), &
         async_id, ng00, i2, blockjj, fulljj, nxjp_acc_length, jjoffset, &
         max_nxjp_acc_length, jbs_nxjp_acc
      integer, dimension(nxjp_acc_length), intent(in) :: map_jj, map_i

      integer, dimension(max_nxjp_acc_length, nlay), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys) :: temcol, colamt1, colamt2, colamt3, &
                               colamt4, colamt5, colamt6, colamt7
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: h2ovmr, o3vmr
      real (kind=kind_phys), dimension(nxptot,nlay), intent(in) :: gasvmr_co2
      real (kind=kind_phys), dimension(nf_vgas), intent(in) :: gasvmr_other
      integer, intent(in) :: nf_vgas
      integer :: k1, rr

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay, nbands), intent(in):: tauaer

      real (kind=kind_phys) :: rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(nxjp_acc_length, nlay, ng00) :: small_fracs, &
     &       small_tautot

!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &           ig, ib

      real (kind=kind_phys) :: corradj, tauself, taufor, taug
      integer :: jj, iplon, iplon2, n
      ! GPU: variable name changed: CPU - colamt,       GPU - colamt1
!
!===> ...  begin here
!
      !$acc parallel loop collapse(2) private(ind0, ind1, inds, indf, ind0p, &
      !$acc&         ind1p, indsp, indfp, corradj, tauself, taufor, taug, colamt1, &
      !$acc&         rr) async(async_id)
      do k = 1, nlay
         do n = 1, nxjp_acc_length
            rr = n + jbs_nxjp_acc - 1
            colamt1 = max(f_zero,    coldry(n, k)*h2ovmr(n, k))          ! h2o
            if (k .le. laytrop(n)) then
               !  --- ...  lower atmosphere loop
               ind0 = ((jp(n, k)-1)*5 + (jt (n, k)-1)) * nspa(2) + 1
               ind1 = ( jp(n, k)   *5 + (jt1(n, k)-1)) * nspa(2) + 1
               inds = indself(n, k)
               indf = indfor(n, k)

               ind0p = ind0 + 1
               ind1p = ind1 + 1
               indsp = inds + 1
               indfp = indf + 1

               corradj = f_one - 0.05 * (pavel(n, k) - 100.0) / 900.0
               !$acc loop seq
               do ig = 1, ng02
                  ib = ngb(ns02+ig)
                  tauself = selffac(n, k) * (selfref(ig,inds) + selffrac(n, k)        &
                  &            * (selfref(ig,indsp) - selfref(ig,inds)))
                  taufor  = forfac(n, k) * (forref(ig,indf) + forfrac(n, k)           &
                  &            * (forref(ig,indfp) - forref(ig,indf))) 

                  taug = corradj * (colamt1                      &
                  &            * (fac00(n, k)*absa(ig,ind0) + fac10(n, k)*absa(ig,ind0p)   &
                  &            +  fac01(n, k)*absa(ig,ind1) + fac11(n, k)*absa(ig,ind1p))  &
                  &            + tauself + taufor)

                  !fracs(iplon, k, ns02+ig, jj) = fracrefa(ig)
                  !tautot(iplon, k, ns02+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) = fracrefa(ig)
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)
               enddo
            else
               ind0 = ((jp(n, k)-13)*5 + (jt (n, k)-1)) * nspb(2) + 1
               ind1 = ((jp(n, k)-12)*5 + (jt1(n, k)-1)) * nspb(2) + 1
               indf = indfor(n, k)

               ind0p = ind0 + 1
               ind1p = ind1 + 1
               indfp = indf + 1
               !$acc loop seq
               do ig = 1, ng02
                  ib = ngb(ns02+ig)
                  taufor = forfac(n, k) * (forref(ig,indf) + forfrac(n, k)            &
                  &           * (forref(ig,indfp) - forref(ig,indf))) 

                  taug = colamt1                                 &
                  &           * (fac00(n, k)*absb(ig,ind0) + fac10(n, k)*absb(ig,ind0p)    &
                  &           +  fac01(n, k)*absb(ig,ind1) + fac11(n, k)*absb(ig,ind1p))   &
                  &           + taufor

                  !fracs(iplon, k, ns02+ig, jj) = fracrefb(ig)
                  !tautot(iplon, k, ns02+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) = fracrefb(ig)
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)
               enddo
            end if
         end do
      end do

! ..................................
      end subroutine taugb02
! ----------------------------------

! ----------------------------------
      subroutine taugb03 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!     band 3:  500-630 cm-1 (low key - h2o,co2; low minor - n2o)       !
!                           (high key - h2o,co2; high minor - n2o)     !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb03
      integer, intent(in) :: nlay, laytrop(nxjp_acc_length), ix, myim(fulljj), &
         async_id, ng00, i2, blockjj, fulljj, nxjp_acc_length, jjoffset, &
         max_nxjp_acc_length, jbs_nxjp_acc
      integer, dimension(nxjp_acc_length), intent(in) :: map_jj, map_i

      integer, dimension(max_nxjp_acc_length, nlay), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys) :: temcol, colamt1, colamt2, colamt3, &
                               colamt4, colamt5, colamt6, colamt7
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: h2ovmr, o3vmr
      real (kind=kind_phys), dimension(nxptot,nlay), intent(in) :: gasvmr_co2
      real (kind=kind_phys), dimension(nf_vgas), intent(in) :: gasvmr_other
      integer, intent(in) :: nf_vgas
      integer :: k1, rr

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay, nbands), intent(in):: tauaer

      real (kind=kind_phys) :: rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(nxjp_acc_length, nlay, ng00) :: small_fracs, &
     &       small_tautot

!  ---  locals:
      integer :: k, ind0, ind1, inds, indsp, indf, indfp, indm, indmp,  &
     &       id000, id010, id100, id110, id200, id210, jmn2o, jmn2op,   &
     &       id001, id011, id101, id111, id201, id211, jpl, jplp,       &
     &       ig, js, js1, jpr

      real (kind=kind_phys) ::  absn2o, ratn2o, adjfac, adjcoln2o,      &
     &      speccomb,       specparm,       specmult,       fs,         &
     &      speccomb1,      specparm1,      specmult1,      fs1,        &
     &      speccomb_mn2o,  specparm_mn2o,  specmult_mn2o,  fmn2o,      &
     &      speccomb_planck,specparm_planck,specmult_planck,fpl,        &
     &      refrat_planck_a, refrat_planck_b, refrat_m_a, refrat_m_b,   &
     &      fac000, fac100, fac200, fac010, fac110, fac210,             &
     &      fac001, fac101, fac201, fac011, fac111, fac211,             &
     &      tau_major, tau_major1, tauself, taufor, n2om1, n2om2,       &
     &      p, p4, fk0, fk1, fk2, taug, colamtr1, colamtr2, colamtr4, k_mn2or0, k_mn2or1
      integer :: jj, iplon, iplon2, ib, n
      ! GPU: variable name changed: CPU - colamt(:,:,1),       GPU - colamt1
      ! GPU: variable name changed: CPU - colamt(:,:,2),       GPU - colamt2
      ! GPU: variable name changed: CPU - colamt(:,:,4),       GPU - colamt4
!
!===> ...  begin here
!
!  --- ...  minor gas mapping levels:
!     lower - n2o, p = 706.272 mbar, t = 278.94 k
!     upper - n2o, p = 95.58 mbar, t = 215.7 k
      refrat_planck_a = chi_mls(1,9)/chi_mls(2,9)    ! p = 212.725 mb
      refrat_planck_b = chi_mls(1,13)/chi_mls(2,13)  ! p = 95.58   mb
      refrat_m_a      = chi_mls(1,3)/chi_mls(2,3)    ! p = 706.270 mb
      refrat_m_b      = chi_mls(1,13)/chi_mls(2,13)  ! p = 95.58   mb

     
      !$acc parallel loop collapse(2) private(speccomb, specparm, specmult, js, &
      !$acc&         fs, ind0, speccomb1, specparm1, specmult1, js1, fs1, ind1, &
      !$acc&         speccomb_mn2o, specparm_mn2o, specmult_mn2o, jmn2o, fmn2o, &
      !$acc&         speccomb_planck, specparm_planck, specmult_planck, jpl, &
      !$acc&         fpl, inds, indf, indm, indsp, indfp, indmp, jmn2op, jplp, p, &
      !$acc&         ratn2o, adjfac, adjcoln2o, p4, fk0, fk1, fk2, id000, id010, &
      !$acc&         id100, id110, id200, id210, fac000, fac100, fac200, fac010, &
      !$acc&         fac110, fac210, id001, id011, id101, id111, id201, id211, &
      !$acc&         fac001, fac101, fac201, fac011, fac111, fac211, tauself, &
      !$acc&         taufor, n2om1, n2om2, absn2o, tau_major, tau_major1, taug, &
      !$acc&         colamt1, colamt2, colamt4, k_mn2or0, k_mn2or1, jpr, rfrate, k1, temcol, rr) async(async_id)
      do k = 1, nlay
         do n = 1, nxjp_acc_length
            rr = n + jbs_nxjp_acc - 1
            if (ivflip == 0) then
               k1 = nlay + 1 - k
            else
               k1 = k
            end if
            temcol = 1.0e-12 * coldry(n, k)
            colamt1 = max(f_zero,    coldry(n, k)*h2ovmr(n, k))          ! h2o
            colamt2 = max(temcol, coldry(n, k)*gasvmr_co2(rr,k1)) ! co2
            if (ilwrgas > 0) then
               colamt4 = max(temcol, coldry(n, k)*gasvmr_other(2))  ! n2o
            else
               colamt4 = f_zero     ! n2o
            endif
            jpr = jp(n, k)
            if (k .le. laytrop(n)) then
            !  --- ...  lower atmosphere loop
               rfrate = chi_mls(1,jpr) / chi_mls(2,jpr)
               speccomb = colamt1 + rfrate*colamt2
               specparm = colamt1 / speccomb
               specmult = 8.0 * min(specparm, oneminus)
               js = 1 + int(specmult)
               fs = mod(specmult, f_one)        
               ind0 = ((jp(n, k)-1)*5 + (jt(n, k)-1)) * nspa(3) + js

               rfrate = chi_mls(1,jpr+1) / chi_mls(2,jpr+1)
               speccomb1 = colamt1 + rfrate*colamt2
               specparm1 = colamt1 / speccomb1
               specmult1 = 8.0 * min(specparm1, oneminus)
               js1 = 1 + int(specmult1)
               fs1 = mod(specmult1, f_one)
               ind1 = (jp(n, k)*5 + (jt1(n, k)-1)) * nspa(3) + js1

               speccomb_mn2o = colamt1 + refrat_m_a*colamt2
               specparm_mn2o = colamt1 / speccomb_mn2o
               specmult_mn2o = 8.0 * min(specparm_mn2o, oneminus)
               jmn2o = 1 + int(specmult_mn2o)
               fmn2o = mod(specmult_mn2o, f_one)

               speccomb_planck = colamt1 + refrat_planck_a*colamt2
               specparm_planck = colamt1 / speccomb_planck
               specmult_planck = 8.0 * min(specparm_planck, oneminus)
               jpl = 1 + int(specmult_planck)
               fpl = mod(specmult_planck, f_one)

               inds = indself(n, k)
               indf = indfor(n, k)
               indm = indminor(n, k)
               indsp = inds + 1
               indfp = indf + 1
               indmp = indm + 1
               jmn2op= jmn2o+ 1
               jplp  = jpl  + 1

               !  --- ...  in atmospheres where the amount of n2o is too great to be considered
               !           a minor species, adjust the column amount of n2o by an empirical factor
               !           to obtain the proper contribution.

               p = coldry(n, k) * chi_mls(4,jp(n, k)+1)
               ratn2o = colamt4 / p
               if (ratn2o > 1.5) then
                  adjfac = 0.5 + (ratn2o - 0.5)**0.65
                  adjcoln2o = adjfac * p
               else
                  adjcoln2o = colamt4
               endif

               if (specparm < 0.125) then
                  p = fs - f_one
                  p4 = p**4
                  fk0 = p4
                  fk1 = f_one - p - 2.0*p4
                  fk2 = p + p4
                  id000 = ind0
                  id010 = ind0 + 9
                  id100 = ind0 + 1
                  id110 = ind0 +10
                  id200 = ind0 + 2
                  id210 = ind0 +11
               else if (specparm > 0.875) then
                  p = -fs
                  p4 = p**4
                  fk0 = p4
                  fk1 = f_one - p - 2.0*p4
                  fk2 = p + p4
                  id000 = ind0 + 1
                  id010 = ind0 +10
                  id100 = ind0
                  id110 = ind0 + 9
                  id200 = ind0 - 1
                  id210 = ind0 + 8
               else
                  fk0 = f_one - fs
                  fk1 = fs
                  fk2 = f_zero
                  id000 = ind0
                  id010 = ind0 + 9
                  id100 = ind0 + 1
                  id110 = ind0 +10
                  id200 = ind0
                  id210 = ind0
               endif

               fac000 = fk0*fac00(n, k)
               fac100 = fk1*fac00(n, k)
               fac200 = fk2*fac00(n, k)
               fac010 = fk0*fac10(n, k)
               fac110 = fk1*fac10(n, k)
               fac210 = fk2*fac10(n, k)

               if (specparm1 < 0.125) then
                  p = fs1 - f_one
                  p4 = p**4
                  fk0 = p4
                  fk1 = f_one - p - 2.0*p4
                  fk2 = p + p4
                  id001 = ind1
                  id011 = ind1 + 9
                  id101 = ind1 + 1
                  id111 = ind1 +10
                  id201 = ind1 + 2
                  id211 = ind1 +11
               elseif (specparm1 > 0.875) then
                  p = -fs1
                  p4 = p**4
                  fk0 = p4
                  fk1 = f_one - p - 2.0*p4
                  fk2 = p + p4
                  id001 = ind1 + 1
                  id011 = ind1 +10
                  id101 = ind1
                  id111 = ind1 + 9
                  id201 = ind1 - 1
                  id211 = ind1 + 8
               else
                  fk0 = f_one - fs1
                  fk1 = fs1
                  fk2 = f_zero
                  id001 = ind1
                  id011 = ind1 + 9
                  id101 = ind1 + 1
                  id111 = ind1 +10
                  id201 = ind1
                  id211 = ind1
               endif

               fac001 = fk0*fac01(n, k)
               fac101 = fk1*fac01(n, k)
               fac201 = fk2*fac01(n, k)
               fac011 = fk0*fac11(n, k)
               fac111 = fk1*fac11(n, k)
               fac211 = fk2*fac11(n, k)
               !$acc loop seq
               do ig = 1, ng03
                  ib = ngb(ns03+ig)
                  k_mn2or0 = ka_mn2o(ig,jmn2o,indm)
                  k_mn2or1 = ka_mn2o(ig,jmn2o,indmp)
                  tauself = selffac(n, k)* (selfref(ig,inds) + selffrac(n, k)         &
                  &            * (selfref(ig,indsp) - selfref(ig,inds)))
                  taufor  = forfac(n, k) * (forref(ig,indf) + forfrac(n, k)           &
                  &            * (forref(ig,indfp) - forref(ig,indf)))
                  n2om1   = k_mn2or0 + fmn2o                      &
                  &            * (ka_mn2o(ig,jmn2op,indm) - k_mn2or0)
                  n2om2   = k_mn2or1 + fmn2o                     &
                  &            * (ka_mn2o(ig,jmn2op,indmp) - k_mn2or1)
                  absn2o  = n2om1 + minorfrac(n, k) * (n2om2 - n2om1)

                  tau_major = speccomb                                          &
                  &              * (fac000*absa(ig,id000) + fac010*absa(ig,id010)    &
                  &              +  fac100*absa(ig,id100) + fac110*absa(ig,id110)    &
                  &              +  fac200*absa(ig,id200) + fac210*absa(ig,id210))

                  tau_major1 = speccomb1                                        &
                  &              * (fac001*absa(ig,id001) + fac011*absa(ig,id011)    &
                  &              +  fac101*absa(ig,id101) + fac111*absa(ig,id111)    &
                  &              +  fac201*absa(ig,id201) + fac211*absa(ig,id211))

                  taug = tau_major + tau_major1                      &
                  &                    + tauself + taufor + adjcoln2o*absn2o

                  !fracs(iplon, k, ns03+ig, jj) = fracrefa(ig,jpl) + fpl                     &
                  !&                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
                  !tautot(iplon, k, ns03+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) = fracrefa(ig,jpl) + fpl                     &
                  &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)
               enddo     ! end do_k_loop
            else
               rfrate = chi_mls(1,jpr) / chi_mls(2,jpr)
               speccomb = colamt1 + rfrate*colamt2
               specparm = colamt1 / speccomb
               specmult = 4.0 * min(specparm, oneminus)
               js = 1 + int(specmult)
               fs = mod(specmult, f_one)
               ind0 = ((jp(n, k)-13)*5 + (jt(n, k)-1)) * nspb(3) + js

               rfrate = chi_mls(1,jpr+1) / chi_mls(2,jpr+1)
               speccomb1 = colamt1 + rfrate*colamt2
               specparm1 = colamt1 / speccomb1
               specmult1 = 4.0 * min(specparm1, oneminus)
               js1 = 1 + int(specmult1)
               fs1 = mod(specmult1, f_one)
               ind1 = ((jp(n, k)-12)*5 + (jt1(n, k)-1)) * nspb(3) + js1

               speccomb_mn2o = colamt1 + refrat_m_b*colamt2
               specparm_mn2o = colamt1 / speccomb_mn2o
               specmult_mn2o = 4.0 * min(specparm_mn2o, oneminus)
               jmn2o = 1 + int(specmult_mn2o)
               fmn2o = mod(specmult_mn2o, f_one)

               speccomb_planck = colamt1 + refrat_planck_b*colamt2
               specparm_planck = colamt1 / speccomb_planck
               specmult_planck = 4.0 * min(specparm_planck, oneminus)
               jpl = 1 + int(specmult_planck)
               fpl = mod(specmult_planck, f_one)

               indf = indfor(n, k)
               indm = indminor(n, k)
               indfp = indf + 1
               indmp = indm + 1
               jmn2op= jmn2o+ 1
               jplp  = jpl  + 1

               id000 = ind0
               id010 = ind0 + 5
               id100 = ind0 + 1
               id110 = ind0 + 6
               id001 = ind1
               id011 = ind1 + 5
               id101 = ind1 + 1
               id111 = ind1 + 6

               !  --- ...  in atmospheres where the amount of n2o is too great to be considered
               !           a minor species, adjust the column amount of n2o by an empirical factor
               !           to obtain the proper contribution.

               p = coldry(n, k) * chi_mls(4,jp(n, k)+1)
               ratn2o = colamt4 / p
               if (ratn2o > 1.5) then
                  adjfac = 0.5 + (ratn2o - 0.5)**0.65
                  adjcoln2o = adjfac * p
               else
                  adjcoln2o = colamt4
               endif

               fk0 = f_one - fs
               fk1 = fs
               fac000 = fk0*fac00(n, k)
               fac010 = fk0*fac10(n, k)
               fac100 = fk1*fac00(n, k)
               fac110 = fk1*fac10(n, k)

               fk0 = f_one - fs1
               fk1 = fs1
               fac001 = fk0*fac01(n, k)
               fac011 = fk0*fac11(n, k)
               fac101 = fk1*fac01(n, k)
               fac111 = fk1*fac11(n, k)
               !$acc loop seq
               do ig = 1, ng03
                  ib = ngb(ns03+ig)
                  k_mn2or0 = kb_mn2o(ig,jmn2o,indm)
                  k_mn2or1 = kb_mn2o(ig,jmn2o,indmp)
                  taufor = forfac(n, k) * (forref(ig,indf) + forfrac(n, k)            &
                  &           * (forref(ig,indfp) - forref(ig,indf))) 
                  n2om1  = k_mn2or0 + fmn2o                       &
                  &           * (kb_mn2o(ig,jmn2op,indm) - k_mn2or0)
                  n2om2  = k_mn2or1 + fmn2o                      &
                  &           * (kb_mn2o(ig,jmn2op,indmp) - k_mn2or1)
                  absn2o = n2om1 + minorfrac(n, k) * (n2om2 - n2om1)

                  tau_major = speccomb                                          &
                  &              * (fac000*absb(ig,id000) + fac010*absb(ig,id010)    &
                  &              +  fac100*absb(ig,id100) + fac110*absb(ig,id110))

                  tau_major1 = speccomb1                                        &
                  &              * (fac001*absb(ig,id001) + fac011*absb(ig,id011)    &
                  &              +  fac101*absb(ig,id101) + fac111*absb(ig,id111))

                  taug = tau_major + tau_major1                      &
                  &                    + taufor + adjcoln2o*absn2o            

                  !fracs(iplon, k, ns03+ig, jj) = fracrefb(ig,jpl) + fpl                     &
                  !&                     * (fracrefb(ig,jplp) - fracrefb(ig,jpl))
                  !tautot(iplon, k, ns03+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) = fracrefb(ig,jpl) + fpl                     &
                  &                     * (fracrefb(ig,jplp) - fracrefb(ig,jpl))
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)
               enddo
            end if
         end do
      end do

! ..................................
      end subroutine taugb03
! ----------------------------------

! ----------------------------------
      subroutine taugb04 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!     band 4:  630-700 cm-1 (low key - h2o,co2; high key - o3,co2)     !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb04
      integer, intent(in) :: nlay, laytrop(nxjp_acc_length), ix, myim(fulljj), &
         async_id, ng00, i2, blockjj, fulljj, nxjp_acc_length, jjoffset, &
         max_nxjp_acc_length, jbs_nxjp_acc
      integer, dimension(nxjp_acc_length), intent(in) :: map_jj, map_i

      integer, dimension(max_nxjp_acc_length, nlay), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys) :: temcol, colamt1, colamt2, colamt3, &
                               colamt4, colamt5, colamt6, colamt7
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: h2ovmr, o3vmr
      real (kind=kind_phys), dimension(nxptot,nlay), intent(in) :: gasvmr_co2
      real (kind=kind_phys), dimension(nf_vgas), intent(in) :: gasvmr_other
      integer, intent(in) :: nf_vgas
      integer :: k1, rr

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay, nbands), intent(in):: tauaer

      real (kind=kind_phys) :: rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(nxjp_acc_length, nlay, ng00) :: small_fracs, &
     &       small_tautot

!  ---  locals:
      integer :: k, ind0, ind1, inds, indsp, indf, indfp, jpl, jplp,    &
     &       id000, id010, id100, id110, id200, id210, ig, js, js1,     &
     &       id001, id011, id101, id111, id201, id211, jpr

      real (kind=kind_phys) :: tauself, taufor, p, p4, fk0, fk1, fk2,   &
     &      speccomb,       specparm,       specmult,       fs,         &
     &      speccomb1,      specparm1,      specmult1,      fs1,        &
     &      speccomb_planck,specparm_planck,specmult_planck,fpl,        &
     &      fac000, fac100, fac200, fac010, fac110, fac210,             &
     &      fac001, fac101, fac201, fac011, fac111, fac211,             &
     &      refrat_planck_a, refrat_planck_b, tau_major, tau_major1, taug
      integer :: jj, iplon, iplon2, ib, n
      ! GPU: variable name changed: CPU - colamt(:,:,1),       GPU - colamt1
      ! GPU: variable name changed: CPU - colamt(:,:,2),       GPU - colamt2
      ! GPU: variable name changed: CPU - colamt(:,:,3),       GPU - colamt3
!
!===> ...  begin here
!
      refrat_planck_a = chi_mls(1,11)/chi_mls(2,11)     ! p = 142.5940 mb
      refrat_planck_b = chi_mls(3,13)/chi_mls(2,13)     ! p = 95.58350 mb
      !$acc parallel loop collapse(2) private(speccomb, specparm, specmult, js, &
      !$acc&         fs, ind0, speccomb1, specparm1, specmult1, js1, fs1, ind1, &
      !$acc&         speccomb_planck, specparm_planck, specmult_planck, jpl, &
      !$acc&         fpl, inds, indf, indsp, indfp, jplp, p, k1, temcol, &
      !$acc&         p4, fk0, fk1, fk2, id000, id010, colamt1, colamt2, colamt3, &
      !$acc&         id100, id110, id200, id210, fac000, fac100, fac200, fac010, &
      !$acc&         fac110, fac210, id001, id011, id101, id111, id201, id211, &
      !$acc&         fac001, fac101, fac201, fac011, fac111, fac211, tauself, &
      !$acc&         taufor, tau_major, tau_major1, taug, jpr, rfrate, rr) async(async_id)
      do k = 1, nlay
         do n = 1, nxjp_acc_length
            jpr = jp(n, k)
            rr = n + jbs_nxjp_acc - 1
            if (ivflip == 0) then
               k1 = nlay + 1 - k
            else
               k1 = k
            end if
            temcol = 1.0e-12 * coldry(n, k)
            colamt1 = max(f_zero,    coldry(n, k)*h2ovmr(n, k))          ! h2o
            colamt2 = max(temcol, coldry(n, k)*gasvmr_co2(rr,k1)) ! co2
            colamt3 = max(temcol, coldry(n, k)*o3vmr(n, k))           ! o3
            if (k .le. laytrop(n)) then
               !  --- ...  lower atmosphere loop
               rfrate = chi_mls(1,jpr) / chi_mls(2,jpr)
               speccomb = colamt1 + rfrate*colamt2
               specparm = colamt1 / speccomb
               specmult = 8.0 * min(specparm, oneminus)
               js = 1 + int(specmult)
               fs = mod(specmult, f_one)
               ind0 = ((jp(n, k)-1)*5 + (jt(n, k)-1)) * nspa(4) + js

               rfrate = chi_mls(1,jpr+1) / chi_mls(2,jpr+1)
               speccomb1 = colamt1 + rfrate*colamt2
               specparm1 = colamt1 / speccomb1
               specmult1 = 8.0 * min(specparm1, oneminus)
               js1 = 1 + int(specmult1)
               fs1 = mod(specmult1, f_one)
               ind1 = ( jp(n, k)*5 + (jt1(n, k)-1)) * nspa(4) + js1

               speccomb_planck = colamt1 + refrat_planck_a*colamt2
               specparm_planck = colamt1 / speccomb_planck
               specmult_planck = 8.0 * min(specparm_planck, oneminus)
               jpl = 1 + int(specmult_planck)
               fpl = mod(specmult_planck, 1.0)

               inds = indself(n, k)
               indf = indfor(n, k)
               indsp = inds + 1
               indfp = indf + 1
               jplp  = jpl  + 1

               if (specparm < 0.125) then
                  p = fs - f_one
                  p4 = p**4
                  fk0 = p4
                  fk1 = f_one - p - 2.0*p4
                  fk2 = p + p4
                  id000 = ind0
                  id010 = ind0 + 9
                  id100 = ind0 + 1
                  id110 = ind0 +10
                  id200 = ind0 + 2
                  id210 = ind0 +11
               elseif (specparm > 0.875) then
                  p = -fs
                  p4 = p**4
                  fk0 = p4
                  fk1 = f_one - p - 2.0*p4
                  fk2 = p + p4
                  id000 = ind0 + 1
                  id010 = ind0 +10
                  id100 = ind0
                  id110 = ind0 + 9
                  id200 = ind0 - 1
                  id210 = ind0 + 8
               else
                  fk0 = f_one - fs
                  fk1 = fs
                  fk2 = f_zero
                  id000 = ind0
                  id010 = ind0 + 9
                  id100 = ind0 + 1
                  id110 = ind0 +10
                  id200 = ind0
                  id210 = ind0
               endif

               fac000 = fk0*fac00(n, k)
               fac100 = fk1*fac00(n, k)
               fac200 = fk2*fac00(n, k)
               fac010 = fk0*fac10(n, k)
               fac110 = fk1*fac10(n, k)
               fac210 = fk2*fac10(n, k)

               if (specparm1 < 0.125) then
                  p = fs1 - f_one
                  p4 = p**4
                  fk0 = p4
                  fk1 = f_one - p - 2.0*p4
                  fk2 = p + p4
                  id001 = ind1
                  id011 = ind1 + 9
                  id101 = ind1 + 1
                  id111 = ind1 +10
                  id201 = ind1 + 2
                  id211 = ind1 +11
               elseif (specparm1 > 0.875) then
                  p = -fs1
                  p4 = p**4
                  fk0 = p4
                  fk1 = f_one - p - 2.0*p4
                  fk2 = p + p4
                  id001 = ind1 + 1
                  id011 = ind1 +10
                  id101 = ind1
                  id111 = ind1 + 9
                  id201 = ind1 - 1
                  id211 = ind1 + 8
               else
                  fk0 = f_one - fs1
                  fk1 = fs1
                  fk2 = f_zero
                  id001 = ind1
                  id011 = ind1 + 9
                  id101 = ind1 + 1
                  id111 = ind1 +10
                  id201 = ind1
                  id211 = ind1
               endif

               fac001 = fk0*fac01(n, k)
               fac101 = fk1*fac01(n, k)
               fac201 = fk2*fac01(n, k)
               fac011 = fk0*fac11(n, k)
               fac111 = fk1*fac11(n, k)
               fac211 = fk2*fac11(n, k)
               !$acc loop seq
               do ig = 1, ng04
                  ib = ngb(ns04+ig)
                  tauself = selffac(n, k)* (selfref(ig,inds) + selffrac(n, k)         &
                  &            * (selfref(ig,indsp) - selfref(ig,inds)))
                  taufor  = forfac(n, k) * (forref(ig,indf) + forfrac(n, k)           &
                  &            * (forref(ig,indfp) - forref(ig,indf))) 

                  tau_major = speccomb                                          &
                  &              * (fac000*absa(ig,id000) + fac010*absa(ig,id010)    &
                  &              +  fac100*absa(ig,id100) + fac110*absa(ig,id110)    &
                  &              +  fac200*absa(ig,id200) + fac210*absa(ig,id210))

                  tau_major1 = speccomb1                                        &
                  &              * (fac001*absa(ig,id001) + fac011*absa(ig,id011)    &
                  &              +  fac101*absa(ig,id101) + fac111*absa(ig,id111)    &
                  &              +  fac201*absa(ig,id201) + fac211*absa(ig,id211))

                  taug = tau_major + tau_major1 + tauself + taufor

                  !fracs(iplon, k, ns04+ig, jj) = fracrefa(ig,jpl) + fpl                     &
                  !&                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
                  !tautot(iplon, k, ns04+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) = fracrefa(ig,jpl) + fpl                     &
                  &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)
               enddo     ! end do_k_loop
            else
               rfrate = chi_mls(3,jpr) / chi_mls(2,jpr)
               speccomb = colamt3 + rfrate*colamt2
               specparm = colamt3 / speccomb
               specmult = 4.0 * min(specparm, oneminus)
               js = 1 + int(specmult)
               fs = mod(specmult, f_one)
               ind0 = ((jp(n, k)-13)*5 + (jt(n, k)-1)) * nspb(4) + js

               rfrate = chi_mls(3,jpr+1) / chi_mls(2,jpr+1)
               speccomb1 = colamt3 + rfrate*colamt2
               specparm1 = colamt3 / speccomb1
               specmult1 = 4.0 * min(specparm1, oneminus)
               js1 = 1 + int(specmult1)
               fs1 = mod(specmult1, f_one)
               ind1 = ((jp(n, k)-12)*5 + (jt1(n, k)-1)) * nspb(4) + js1

               speccomb_planck = colamt3 + refrat_planck_b*colamt2
               specparm_planck = colamt3 / speccomb_planck
               specmult_planck = 4.0 * min(specparm_planck, oneminus)
               jpl = 1 + int(specmult_planck)
               fpl = mod(specmult_planck, f_one)
               jplp = jpl + 1

               id000 = ind0
               id010 = ind0 + 5
               id100 = ind0 + 1
               id110 = ind0 + 6
               id001 = ind1
               id011 = ind1 + 5
               id101 = ind1 + 1
               id111 = ind1 + 6

               fk0 = f_one - fs
               fk1 = fs
               fac000 = fk0*fac00(n, k)
               fac010 = fk0*fac10(n, k)
               fac100 = fk1*fac00(n, k)
               fac110 = fk1*fac10(n, k)

               fk0 = f_one - fs1
               fk1 = fs1
               fac001 = fk0*fac01(n, k)
               fac011 = fk0*fac11(n, k)
               fac101 = fk1*fac01(n, k)
               fac111 = fk1*fac11(n, k)
               !$acc loop seq
               do ig = 1, ng04
                  ib = ngb(ns04+ig)
                  tau_major =  speccomb                                         &
                  &              * (fac000*absb(ig,id000) + fac010*absb(ig,id010)    &
                  &              +  fac100*absb(ig,id100) + fac110*absb(ig,id110))
                  tau_major1 = speccomb1                                        &
                  &              * (fac001*absb(ig,id001) + fac011*absb(ig,id011)    &
                  &              +  fac101*absb(ig,id101) + fac111*absb(ig,id111))

                  taug =  tau_major + tau_major1

                  !fracs(iplon, k, ns04+ig, jj) = fracrefb(ig,jpl) + fpl                     &
                  !&                     * (fracrefb(ig,jplp) - fracrefb(ig,jpl))
                  small_fracs(n, k, ig) = fracrefb(ig,jpl) + fpl                     &
                  &                     * (fracrefb(ig,jplp) - fracrefb(ig,jpl))
               
                  !  --- ...  empirical modification to code to improve stratospheric cooling rates
               !           for co2. revised to apply weighting for g-point reduction in this band.
                  if (ig .eq. 8) taug = taug * 0.92
                  if (ig .eq. 9) taug = taug * 0.88
                  if (ig .eq. 10) taug = taug * 1.07
                  if (ig .eq. 11) taug = taug * 1.1
                  if (ig .eq. 12) taug = taug * 0.99
                  if (ig .eq. 13) taug = taug * 0.88
                  if (ig .eq. 14) taug = taug * 0.943
                  !tautot(iplon, k, ns04+ig, jj) = taug + tauaer(n, k, ib)
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)
               enddo
            end if
         end do
      end do

! ..................................
      end subroutine taugb04
! ----------------------------------

! ----------------------------------
      subroutine taugb05 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!     band 5:  700-820 cm-1 (low key - h2o,co2; low minor - o3, ccl4)  !
!                           (high key - o3,co2)                        !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb05
      integer, intent(in) :: nlay, laytrop(nxjp_acc_length), ix, myim(fulljj), &
         async_id, ng00, i2, blockjj, fulljj, nxjp_acc_length, jjoffset, &
         max_nxjp_acc_length, jbs_nxjp_acc
      integer, dimension(nxjp_acc_length), intent(in) :: map_jj, map_i

      integer, dimension(max_nxjp_acc_length, nlay), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys) :: temcol, colamt1, colamt2, colamt3, &
                               colamt4, colamt5, colamt6, colamt7, wx1
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: h2ovmr, o3vmr
      real (kind=kind_phys), dimension(nxptot,nlay), intent(in) :: gasvmr_co2
      real (kind=kind_phys), dimension(nf_vgas), intent(in) :: gasvmr_other
      integer, intent(in) :: nf_vgas
      integer :: k1, rr

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay, nbands), intent(in):: tauaer

      real (kind=kind_phys) :: rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(nxjp_acc_length, nlay, ng00) :: small_fracs, &
     &       small_tautot

!  ---  locals: 
      integer :: k, ind0, ind1, inds, indsp, indf, indfp, indm, indmp,  &
     &       id000, id010, id100, id110, id200, id210, jmo3, jmo3p,     &
     &       id001, id011, id101, id111, id201, id211, jpl, jplp,       &
     &       ig, js, js1, jpr

      real (kind=kind_phys)  :: tauself, taufor, o3m1, o3m2, abso3,     &
     &      speccomb,       specparm,       specmult,       fs,         &
     &      speccomb1,      specparm1,      specmult1,      fs1,        &
     &      speccomb_mo3,   specparm_mo3,   specmult_mo3,   fmo3,       &
     &      speccomb_planck,specparm_planck,specmult_planck,fpl,        &
     &      refrat_planck_a, refrat_planck_b, refrat_m_a,               &
     &      fac000, fac100, fac200, fac010, fac110, fac210,             &
     &      fac001, fac101, fac201, fac011, fac111, fac211,             &
     &      p0, p40, fk00, fk10, fk20, p1, p41, fk01, fk11, fk21, taug
      integer :: jj, iplon, iplon2, ib, n
      ! GPU: variable name changed: CPU - wx(:,:,1),           GPU - wx1
      ! GPU: variable name changed: CPU - colamt(:,:,1),       GPU - colamt1
      ! GPU: variable name changed: CPU - colamt(:,:,2),       GPU - colamt2
      ! GPU: variable name changed: CPU - colamt(:,:,3),       GPU - colamt3
!
!===> ...  begin here
!
!  --- ...  minor gas mapping level :
!     lower - o3, p = 317.34 mbar, t = 240.77 k
!     lower - ccl4

!  --- ...  calculate reference ratio to be used in calculation of planck
!           fraction in lower/upper atmosphere.
      refrat_planck_a = chi_mls(1,5)/chi_mls(2,5)      ! p = 473.420 mb
      refrat_planck_b = chi_mls(3,43)/chi_mls(2,43)    ! p = 0.2369  mb
      refrat_m_a = chi_mls(1,7)/chi_mls(2,7)           ! p = 317.348 mb

      !$acc parallel loop collapse(2) private(speccomb, specparm, specmult, js, &
      !$acc&         fs, ind0, speccomb1, specparm1, specmult1, js1, fs1, ind1, &
      !$acc&         speccomb_planck, specparm_planck, specmult_planck, jpl, &
      !$acc&         fpl, inds, indf, indm, indsp, indfp, indmp, jplp, &
      !$acc&         id000, id010, speccomb_mo3, specparm_mo3, specmult_mo3, &
      !$acc&         id100, id110, id200, id210, fac000, fac100, fac200, fac010, &
      !$acc&         fac110, fac210, id001, id011, id101, id111, id201, id211, &
      !$acc&         fac001, fac101, fac201, fac011, fac111, fac211, tauself, &
      !$acc&         taufor, jmo3, fmo3, jmo3p, p0, p40, fk00, fk10, fk20, p1, &
      !$acc&         p41, fk01, fk11, fk21, o3m1, o3m2, abso3, taug, jpr, rfrate, &
      !$acc&         temcol, k1, colamt1, colamt2, colamt3, wx1, rr) async(async_id)
      do k = 1, nlay
         do n = 1, nxjp_acc_length
            jpr = jp(n, k)
            rr = n + jbs_nxjp_acc - 1
            if (ivflip == 0) then
               k1 = nlay + 1 - k
            else
               k1 = k
            end if
            temcol = 1.0e-12 * coldry(n, k)
            colamt1 = max(f_zero,    coldry(n, k)*h2ovmr(n, k))          ! h2o
            colamt2 = max(temcol, coldry(n, k)*gasvmr_co2(rr,k1)) ! co2
            colamt3 = max(temcol, coldry(n, k)*o3vmr(n, k))           ! o3
            if (ilwrgas > 0) then
               wx1 = max( f_zero, coldry(n, k)*gasvmr_other(9) )   ! ccl4
            else
               wx1 = f_zero
            endif
            if (k .le. laytrop(n)) then
               !  --- ...  lower atmosphere loop
               rfrate = chi_mls(1,jpr) / chi_mls(2,jpr)
               speccomb = colamt1 + rfrate*colamt2
               specparm = colamt1 / speccomb
               specmult = 8.0 * min(specparm, oneminus)
               js = 1 + int(specmult)
               fs = mod(specmult, f_one)
               ind0 = ((jp(n, k)-1)*5 + (jt(n, k)-1)) * nspa(5) + js

               rfrate = chi_mls(1,jpr+1) / chi_mls(2,jpr+1)
               speccomb1 = colamt1 + rfrate*colamt2
               specparm1 = colamt1 / speccomb1
               specmult1 = 8.0 * min(specparm1, oneminus)
               js1 = 1 + int(specmult1)
               fs1 = mod(specmult1, f_one)
               ind1 = (jp(n, k)*5 + (jt1(n, k)-1)) * nspa(5) + js1

               speccomb_mo3 = colamt1 + refrat_m_a*colamt2
               specparm_mo3 = colamt1 / speccomb_mo3
               specmult_mo3 = 8.0 * min(specparm_mo3, oneminus)
               jmo3 = 1 + int(specmult_mo3)
               fmo3 = mod(specmult_mo3, f_one)

               speccomb_planck = colamt1 + refrat_planck_a*colamt2
               specparm_planck = colamt1 / speccomb_planck
               specmult_planck = 8.0 * min(specparm_planck, oneminus)
               jpl = 1 + int(specmult_planck)
               fpl = mod(specmult_planck, f_one)

               inds = indself(n, k)
               indf = indfor(n, k)
               indm = indminor(n, k)
               indsp = inds + 1
               indfp = indf + 1
               indmp = indm + 1
               jplp  = jpl  + 1
               jmo3p = jmo3 + 1

               if (specparm < 0.125 .and. specparm1 < 0.125) then
                  p0 = fs - f_one
                  p40 = p0**4
                  fk00 = p40
                  fk10 = f_one - p0 - 2.0*p40
                  fk20 = p0 + p40

                  p1 = fs1 - f_one
                  p41 = p1**4
                  fk01 = p41
                  fk11 = f_one - p1 - 2.0*p41
                  fk21 = p1 + p41

                  id000 = ind0
                  id010 = ind0 + 9
                  id100 = ind0 + 1
                  id110 = ind0 +10
                  id200 = ind0 + 2
                  id210 = ind0 +11

                  id001 = ind1
                  id011 = ind1 + 9
                  id101 = ind1 + 1
                  id111 = ind1 +10
                  id201 = ind1 + 2
                  id211 = ind1 +11
               elseif (specparm > 0.875 .and. specparm1 > 0.875) then
                  p0 = -fs
                  p40 = p0**4
                  fk00 = p40
                  fk10 = f_one - p0 - 2.0*p40
                  fk20 = p0 + p40

                  p1 = -fs1
                  p41 = p1**4
                  fk01 = p41
                  fk11 = f_one - p1 - 2.0*p41
                  fk21 = p1 + p41

                  id000 = ind0 + 1
                  id010 = ind0 +10
                  id100 = ind0
                  id110 = ind0 + 9
                  id200 = ind0 - 1
                  id210 = ind0 + 8

                  id001 = ind1 + 1
                  id011 = ind1 +10
                  id101 = ind1
                  id111 = ind1 + 9
                  id201 = ind1 - 1
                  id211 = ind1 + 8
               else
                  fk00 = f_one - fs
                  fk10 = fs
                  fk20 = f_zero

                  fk01 = f_one - fs1
                  fk11 = fs1
                  fk21 = f_zero

                  id000 = ind0
                  id010 = ind0 + 9
                  id100 = ind0 + 1
                  id110 = ind0 +10
                  id200 = ind0
                  id210 = ind0

                  id001 = ind1
                  id011 = ind1 + 9
                  id101 = ind1 + 1
                  id111 = ind1 +10
                  id201 = ind1
                  id211 = ind1
               endif

               fac000 = fk00 * fac00(n, k)
               fac100 = fk10 * fac00(n, k)
               fac200 = fk20 * fac00(n, k)
               fac010 = fk00 * fac10(n, k)
               fac110 = fk10 * fac10(n, k)
               fac210 = fk20 * fac10(n, k)

               fac001 = fk01 * fac01(n, k)
               fac101 = fk11 * fac01(n, k)
               fac201 = fk21 * fac01(n, k)
               fac011 = fk01 * fac11(n, k)
               fac111 = fk11 * fac11(n, k)
               fac211 = fk21 * fac11(n, k)
               !$acc loop seq
               do ig = 1, ng05
                  ib = ngb(ns05+ig)
                  tauself = selffac(n, k) * (selfref(ig,inds) + selffrac(n, k)        &
                  &            * (selfref(ig,indsp) - selfref(ig,inds)))
                  taufor  = forfac(n, k) * (forref(ig,indf) + forfrac(n, k)           &
                  &            * (forref(ig,indfp) - forref(ig,indf)))
                  o3m1    = ka_mo3(ig,jmo3,indm) + fmo3                         &
                  &            * (ka_mo3(ig,jmo3p,indm) -  ka_mo3(ig,jmo3,indm))
                  o3m2    = ka_mo3(ig,jmo3,indmp) + fmo3                        &
                  &            * (ka_mo3(ig,jmo3p,indmp) - ka_mo3(ig,jmo3,indmp))
                  abso3   = o3m1 + minorfrac(n, k)*(o3m2 - o3m1)

                  taug = speccomb                                    &
                  &            * (fac000*absa(ig,id000) + fac010*absa(ig,id010)      &
                  &            +  fac100*absa(ig,id100) + fac110*absa(ig,id110)      &
                  &            +  fac200*absa(ig,id200) + fac210*absa(ig,id210))     &
                  &            +     speccomb1                                       &
                  &            * (fac001*absa(ig,id001) + fac011*absa(ig,id011)      &
                  &            +  fac101*absa(ig,id101) + fac111*absa(ig,id111)      &
                  &            +  fac201*absa(ig,id201) + fac211*absa(ig,id211))     &
                  &            + tauself + taufor+abso3*colamt3+wx1*ccl4(ig)

                  !fracs(iplon, k, ns05+ig, jj) = fracrefa(ig,jpl) + fpl                     &
                  !&                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
                  !tautot(iplon, k, ns05+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) = fracrefa(ig,jpl) + fpl                     &
                  &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)
               enddo
            else
               rfrate = chi_mls(3,jpr) / chi_mls(2,jpr)
               speccomb = colamt3 + rfrate*colamt2
               specparm = colamt3 / speccomb
               specmult = 4.0 * min(specparm, oneminus)
               js = 1 + int(specmult)
               fs = mod(specmult, f_one)
               ind0 = ((jp(n, k)-13)*5 + (jt(n, k)-1)) * nspb(5) + js

               rfrate = chi_mls(3,jpr+1) / chi_mls(2,jpr+1)
               speccomb1 = colamt3 + rfrate*colamt2
               specparm1 = colamt3 / speccomb1
               specmult1 = 4.0 * min(specparm1, oneminus)
               js1 = 1 + int(specmult1)
               fs1 = mod(specmult1, f_one)
               ind1 = ((jp(n, k)-12)*5 + (jt1(n, k)-1)) * nspb(5) + js1

               speccomb_planck = colamt3 + refrat_planck_b*colamt2
               specparm_planck = colamt3 / speccomb_planck
               specmult_planck = 4.0 * min(specparm_planck, oneminus)
               jpl = 1 + int(specmult_planck)
               fpl = mod(specmult_planck, f_one)
               jplp= jpl + 1

               id000 = ind0
               id010 = ind0 + 5
               id100 = ind0 + 1
               id110 = ind0 + 6
               id001 = ind1
               id011 = ind1 + 5
               id101 = ind1 + 1
               id111 = ind1 + 6

               fk00 = f_one - fs
               fk10 = fs

               fk01 = f_one - fs1
               fk11 = fs1

               fac000 = fk00 * fac00(n, k)
               fac010 = fk00 * fac10(n, k)
               fac100 = fk10 * fac00(n, k)
               fac110 = fk10 * fac10(n, k)

               fac001 = fk01 * fac01(n, k)
               fac011 = fk01 * fac11(n, k)
               fac101 = fk11 * fac01(n, k)
               fac111 = fk11 * fac11(n, k)
               !$acc loop seq
               do ig = 1, ng05
                  ib = ngb(ns05+ig)
                  taug = speccomb                                    &
                  &                * (fac000*absb(ig,id000) + fac010*absb(ig,id010)  &
                  &                +  fac100*absb(ig,id100) + fac110*absb(ig,id110)) &
                  &                +     speccomb1                                   &
                  &                * (fac001*absb(ig,id001) + fac011*absb(ig,id011)  &
                  &                +  fac101*absb(ig,id101) + fac111*absb(ig,id111)) &
                  &                + wx1 * ccl4(ig)

                  !fracs(iplon, k, ns05+ig, jj) = fracrefb(ig,jpl) + fpl                     &
                  !&                     * (fracrefb(ig,jplp) - fracrefb(ig,jpl))
                  !tautot(iplon, k, ns05+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) = fracrefb(ig,jpl) + fpl                     &
                  &                     * (fracrefb(ig,jplp) - fracrefb(ig,jpl))
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)
               enddo
            end if
         end do
      end do

! ..................................
      end subroutine taugb05
! ----------------------------------

! ----------------------------------
      subroutine taugb06 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!     band 6:  820-980 cm-1 (low key - h2o; low minor - co2)           !
!                           (high key - none; high minor - cfc11, cfc12)
!  ------------------------------------------------------------------  !

      use module_radlw_kgb06
      integer, intent(in) :: nlay, laytrop(nxjp_acc_length), ix, myim(fulljj), &
         async_id, ng00, i2, blockjj, fulljj, nxjp_acc_length, jjoffset, &
         max_nxjp_acc_length, jbs_nxjp_acc
      integer, dimension(nxjp_acc_length), intent(in) :: map_jj, map_i

      integer, dimension(max_nxjp_acc_length, nlay), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys) :: temcol, colamt1, colamt2, colamt3, &
                               colamt4, colamt5, colamt6, colamt7, wx2, wx3
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: h2ovmr, o3vmr
      real (kind=kind_phys), dimension(nxptot,nlay), intent(in) :: gasvmr_co2
      real (kind=kind_phys), dimension(nf_vgas), intent(in) :: gasvmr_other
      integer, intent(in) :: nf_vgas
      integer :: k1, rr

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay, nbands), intent(in):: tauaer

      real (kind=kind_phys) :: rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(nxjp_acc_length, nlay, ng00) :: small_fracs, &
     &       small_tautot

!  ---  locals: 
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       indm, indmp, ig

      real (kind=kind_phys) :: ratco2, adjfac, adjcolco2, tauself,      &
     &      taufor, absco2, temp, taug
      integer :: jj, iplon, iplon2, ib, n
      ! GPU: variable name changed: CPU - wx(:,:,2),           GPU - wx2
      ! GPU: variable name changed: CPU - wx(:,:,3),           GPU - wx3
      ! GPU: variable name changed: CPU - colamt(:,:,1),       GPU - colamt1
      ! GPU: variable name changed: CPU - colamt(:,:,2),       GPU - colamt2
!
!===> ...  begin here
!
!  --- ...  minor gas mapping level:
!     lower - co2, p = 706.2720 mb, t = 294.2 k
!     upper - cfc11, cfc12

      !$acc parallel loop collapse(2) private(ind0, ind1, inds, indf, k1, temcol, &
      !$acc&         indm, indsp, indfp, indmp, ind0p, ind1p, temp,  ratco2, adjfac, &
      !$acc&         adjcolco2, tauself, taufor, absco2, taug, colamt1, colamt2, wx2, wx3, rr) async(async_id)
      do k = 1, nlay
         do n = 1, nxjp_acc_length
            rr = n + jbs_nxjp_acc - 1
            if (ivflip == 0) then
               k1 = nlay + 1 - k
            else
               k1 = k
            end if
            temcol = 1.0e-12 * coldry(n, k)
            colamt1 = max(f_zero,    coldry(n, k)*h2ovmr(n, k))          ! h2o
            colamt2 = max(temcol, coldry(n, k)*gasvmr_co2(rr,k1)) ! co2
            if (ilwrgas > 0) then
               wx2 = max( f_zero, coldry(n, k)*gasvmr_other(6) )   ! cf11
               wx3 = max( f_zero, coldry(n, k)*gasvmr_other(7) )   ! cf12
            else
               wx2 = f_zero
               wx3 = f_zero
            endif
            if (k .le. laytrop(n)) then
               !  --- ...  lower atmosphere loop
               ind0 = ((jp(n, k)-1)*5 + (jt (n, k)-1)) * nspa(6) + 1
               ind1 = ( jp(n, k)   *5 + (jt1(n, k)-1)) * nspa(6) + 1

               inds = indself(n, k)
               indf = indfor(n, k)
               indm = indminor(n, k)
               indsp = inds + 1
               indfp = indf + 1
               indmp = indm + 1
               ind0p = ind0 + 1
               ind1p = ind1 + 1

               !  --- ...  in atmospheres where the amount of co2 is too great to be considered
               !           a minor species, adjust the column amount of co2 by an empirical factor
               !           to obtain the proper contribution.

               temp   = coldry(n, k) * chi_mls(2,jp(n, k)+1)
               ratco2 = colamt2 / temp
               if (ratco2 > 3.0) then
                  adjfac = 2.0 + (ratco2-2.0)**0.77
                  adjcolco2 = adjfac * temp
               else
                  adjcolco2 = colamt2
               endif
               !$acc loop seq
               do ig = 1, ng06
                  ib = ngb(ns06+ig)
                  tauself = selffac(n, k) * (selfref(ig,inds) + selffrac(n, k)        &
                  &            * (selfref(ig,indsp) - selfref(ig,inds)))
                  taufor  = forfac(n, k) * (forref(ig,indf) + forfrac(n, k)           &
                  &            * (forref(ig,indfp) - forref(ig,indf)))
                  absco2  = ka_mco2(ig,indm) + minorfrac(n, k)                     &
                  &            * (ka_mco2(ig,indmp) - ka_mco2(ig,indm))

                  taug = colamt1                                 &
                  &            * (fac00(n, k)*absa(ig,ind0) + fac10(n, k)*absa(ig,ind0p)   &
                  &            +  fac01(n, k)*absa(ig,ind1) + fac11(n, k)*absa(ig,ind1p))  &
                  &            +  tauself + taufor + adjcolco2*absco2                &
                  &            +  wx2*cfc11adj(ig) + wx3*cfc12(ig)

                  !fracs(iplon, k, ns06+ig, jj) = fracrefa(ig)
                  !tautot(iplon, k, ns06+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) = fracrefa(ig)
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)
               enddo
            else
               !$acc loop seq
               do ig = 1, ng06
                  ib = ngb(ns06+ig)
                  taug = wx2*cfc11adj(ig) + wx3*cfc12(ig)

                  !fracs(iplon, k, ns06+ig, jj) = fracrefa(ig)
                  !tautot(iplon, k, ns06+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) = fracrefa(ig)
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)
               enddo
            end if
         end do
      end do

! ..................................
      end subroutine taugb06
! ----------------------------------

! ----------------------------------
      subroutine taugb07 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!     band 7:  980-1080 cm-1 (low key - h2o,o3; low minor - co2)       !
!                            (high key - o3; high minor - co2)         !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb07
      integer, intent(in) :: nlay, laytrop(nxjp_acc_length), ix, myim(fulljj), &
         async_id, ng00, i2, blockjj, fulljj, nxjp_acc_length, jjoffset, &
         max_nxjp_acc_length, jbs_nxjp_acc
      integer, dimension(nxjp_acc_length), intent(in) :: map_jj, map_i

      integer, dimension(max_nxjp_acc_length, nlay), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys) :: temcol, colamt1, colamt2, colamt3, &
                               colamt4, colamt5, colamt6, colamt7
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: h2ovmr, o3vmr
      real (kind=kind_phys), dimension(nxptot,nlay), intent(in) :: gasvmr_co2
      real (kind=kind_phys), dimension(nf_vgas), intent(in) :: gasvmr_other
      integer, intent(in) :: nf_vgas
      integer :: k1, rr

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay, nbands), intent(in):: tauaer

      real (kind=kind_phys) :: rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(nxjp_acc_length, nlay, ng00) :: small_fracs, &
     &       small_tautot

!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       id000, id010, id100, id110, id200, id210, indm, indmp,     &
     &       id001, id011, id101, id111, id201, id211, jmco2, jmco2p,   &
     &       jpl, jplp, ig, js, js1, jpr

      real (kind=kind_phys) :: tauself, taufor, co2m1, co2m2, absco2,   &
     &      speccomb,       specparm,       specmult,       fs,         &
     &      speccomb1,      specparm1,      specmult1,      fs1,        &
     &      speccomb_mco2,  specparm_mco2,  specmult_mco2,  fmco2,      &
     &      speccomb_planck,specparm_planck,specmult_planck,fpl,        &
     &      refrat_planck_a, refrat_m_a, ratco2, adjfac, adjcolco2,     &
     &      fac000, fac100, fac200, fac010, fac110, fac210,             &
     &      fac001, fac101, fac201, fac011, fac111, fac211,             &
     &      p0, p40, fk00, fk10, fk20, p1, p41, fk01, fk11, fk21, temp, taug
      integer :: jj, iplon, iplon2, ib, n
      ! GPU: variable name changed: CPU - colamt(:,:,1),       GPU - colamt1
      ! GPU: variable name changed: CPU - colamt(:,:,2),       GPU - colamt2
      ! GPU: variable name changed: CPU - colamt(:,:,3),       GPU - colamt3
!
!===> ...  begin here
!
!  --- ...  minor gas mapping level :
!     lower - co2, p = 706.2620 mbar, t= 278.94 k
!     upper - co2, p = 12.9350 mbar, t = 234.01 k

!  --- ...  calculate reference ratio to be used in calculation of planck
!           fraction in lower atmosphere.
      refrat_planck_a = chi_mls(1,3)/chi_mls(3,3)     ! p = 706.2620 mb
      refrat_m_a = chi_mls(1,3)/chi_mls(3,3)          ! p = 706.2720 mb

      !$acc parallel loop collapse(2) private(speccomb, specparm, specmult, js, &
      !$acc&         fs, ind0, speccomb1, specparm1, specmult1, js1, fs1, ind1, &
      !$acc&         speccomb_planck, specparm_planck, specmult_planck, jpl, &
      !$acc&         fpl, inds, indf, indm, indsp, indfp, indmp, jplp, &
      !$acc&         id000, id010, speccomb_mco2, specparm_mco2, specmult_mco2, &
      !$acc&         id100, id110, id200, id210, fac000, fac100, fac200, fac010, &
      !$acc&         fac110, fac210, id001, id011, id101, id111, id201, id211, &
      !$acc&         fac001, fac101, fac201, fac011, fac111, fac211, tauself, &
      !$acc&         taufor, p0, p40, fk00, fk10, fk20, p1, jmco2, fmco2, jmco2p, &
      !$acc&         p41, fk01, fk11, fk21, ind0p, ind1p, temp, ratco2, adjfac, &
      !$acc&         adjcolco2, co2m1, co2m2, absco2, taug, jpr, rfrate, &
      !$acc&         k1, temcol, colamt1, colamt2, colamt3, rr) async(async_id)
      do k = 1, nlay
         do n = 1, nxjp_acc_length
            jpr = jp(n, k)
            rr = n + jbs_nxjp_acc - 1
            if (ivflip == 0) then
               k1 = nlay + 1 - k
            else
               k1 = k
            end if
            temcol = 1.0e-12 * coldry(n, k)
            colamt1 = max(f_zero,    coldry(n, k)*h2ovmr(n, k))          ! h2o
            colamt2 = max(temcol, coldry(n, k)*gasvmr_co2(rr,k1)) ! co2
            colamt3 = max(temcol, coldry(n, k)*o3vmr(n, k))           ! o3
            if (k .le. laytrop(n)) then
               !  --- ...  lower atmosphere loop
               rfrate = chi_mls(1,jpr) / chi_mls(3,jpr)
               speccomb = colamt1 + rfrate*colamt3
               specparm = colamt1 / speccomb
               specmult = 8.0 * min(specparm, oneminus)
               js = 1 + int(specmult)
               fs = mod(specmult, f_one)
               ind0 = ((jp(n, k)-1)*5 + (jt(n, k)-1)) * nspa(7) + js

               rfrate = chi_mls(1,jpr+1) / chi_mls(3,jpr+1)
               speccomb1 = colamt1 + rfrate*colamt3
               specparm1 = colamt1 / speccomb1
               specmult1 = 8.0 * min(specparm1, oneminus)
               js1 = 1 + int(specmult1)
               fs1 = mod(specmult1, f_one)
               ind1 = (jp(n, k)*5 + (jt1(n, k)-1)) * nspa(7) + js1

               speccomb_mco2 = colamt1 + refrat_m_a*colamt3
               specparm_mco2 = colamt1 / speccomb_mco2
               specmult_mco2 = 8.0 * min(specparm_mco2, oneminus)
               jmco2 = 1 + int(specmult_mco2)
               fmco2 = mod(specmult_mco2, f_one)

               speccomb_planck = colamt1 + refrat_planck_a*colamt3
               specparm_planck = colamt1 / speccomb_planck
               specmult_planck = 8.0 * min(specparm_planck, oneminus)
               jpl = 1 + int(specmult_planck)
               fpl = mod(specmult_planck, f_one)

               inds = indself(n, k)
               indf = indfor(n, k)
               indm = indminor(n, k)
               indsp = inds + 1
               indfp = indf + 1
               indmp = indm + 1
               jplp  = jpl  + 1
               jmco2p= jmco2+ 1
               ind0p = ind0 + 1
               ind1p = ind1 + 1

               !  --- ...  in atmospheres where the amount of co2 is too great to be considered
               !           a minor species, adjust the column amount of co2 by an empirical factor
               !           to obtain the proper contribution.

               temp   = coldry(n, k) * chi_mls(2,jp(n, k)+1)
               ratco2 = colamt2 / temp
               if (ratco2 > 3.0) then
                  adjfac = 3.0 + (ratco2-3.0)**0.79
                  adjcolco2 = adjfac * temp
               else
                  adjcolco2 = colamt2
               endif

               if (specparm < 0.125 .and. specparm1 < 0.125) then
                  p0 = fs - f_one
                  p40 = p0**4
                  fk00 = p40
                  fk10 = f_one - p0 - 2.0*p40
                  fk20 = p0 + p40

                  p1 = fs1 - f_one
                  p41 = p1**4
                  fk01 = p41
                  fk11 = f_one - p1 - 2.0*p41
                  fk21 = p1 + p41

                  id000 = ind0
                  id010 = ind0 + 9
                  id100 = ind0 + 1
                  id110 = ind0 +10
                  id200 = ind0 + 2
                  id210 = ind0 +11

                  id001 = ind1
                  id011 = ind1 + 9
                  id101 = ind1 + 1
                  id111 = ind1 +10
                  id201 = ind1 + 2
                  id211 = ind1 +11
               elseif (specparm > 0.875 .and. specparm1 > 0.875) then
                  p0 = -fs
                  p40 = p0**4
                  fk00 = p40
                  fk10 = f_one - p0 - 2.0*p40
                  fk20 = p0 + p40

                  p1 = -fs1
                  p41 = p1**4
                  fk01 = p41
                  fk11 = f_one - p1 - 2.0*p41
                  fk21 = p1 + p41

                  id000 = ind0 + 1
                  id010 = ind0 +10
                  id100 = ind0
                  id110 = ind0 + 9
                  id200 = ind0 - 1
                  id210 = ind0 + 8

                  id001 = ind1 + 1
                  id011 = ind1 +10
                  id101 = ind1
                  id111 = ind1 + 9
                  id201 = ind1 - 1
                  id211 = ind1 + 8
               else
                  fk00 = f_one - fs
                  fk10 = fs
                  fk20 = f_zero

                  fk01 = f_one - fs1
                  fk11 = fs1
                  fk21 = f_zero

                  id000 = ind0
                  id010 = ind0 + 9
                  id100 = ind0 + 1
                  id110 = ind0 +10
                  id200 = ind0
                  id210 = ind0

                  id001 = ind1
                  id011 = ind1 + 9
                  id101 = ind1 + 1
                  id111 = ind1 +10
                  id201 = ind1
                  id211 = ind1
               endif

               fac000 = fk00 * fac00(n, k)
               fac100 = fk10 * fac00(n, k)
               fac200 = fk20 * fac00(n, k)
               fac010 = fk00 * fac10(n, k)
               fac110 = fk10 * fac10(n, k)
               fac210 = fk20 * fac10(n, k)

               fac001 = fk01 * fac01(n, k)
               fac101 = fk11 * fac01(n, k)
               fac201 = fk21 * fac01(n, k)
               fac011 = fk01 * fac11(n, k)
               fac111 = fk11 * fac11(n, k)
               fac211 = fk21 * fac11(n, k)
               !$acc loop seq
               do ig = 1, ng07
                  ib = ngb(ns07+ig)
                  tauself = selffac(n, k)* (selfref(ig,inds) + selffrac(n, k)         &
                  &            * (selfref(ig,indsp) - selfref(ig,inds)))
                  taufor  = forfac(n, k) * (forref(ig,indf) + forfrac(n, k)           &
                  &            * (forref(ig,indfp) - forref(ig,indf))) 
                  co2m1   = ka_mco2(ig,jmco2,indm) + fmco2                      &
                  &            * (ka_mco2(ig,jmco2p,indm) - ka_mco2(ig,jmco2,indm))
                  co2m2   = ka_mco2(ig,jmco2,indmp) + fmco2                     &
                  &            * (ka_mco2(ig,jmco2p,indmp) - ka_mco2(ig,jmco2,indmp))
                  absco2  = co2m1 + minorfrac(n, k) * (co2m2 - co2m1)

                  taug = speccomb                                    &
                  &                * (fac000*absa(ig,id000) + fac010*absa(ig,id010)  &
                  &                +  fac100*absa(ig,id100) + fac110*absa(ig,id110)  &
                  &                +  fac200*absa(ig,id200) + fac210*absa(ig,id210)) &
                  &                +     speccomb1                                   &
                  &                * (fac001*absa(ig,id001) + fac011*absa(ig,id011)  &
                  &                +  fac101*absa(ig,id101) + fac111*absa(ig,id111)  &
                  &                +  fac201*absa(ig,id201) + fac211*absa(ig,id211)) &
                  &                + tauself + taufor + adjcolco2*absco2

                  !fracs(iplon, k, ns07+ig, jj) = fracrefa(ig,jpl) + fpl                     &
                  !&                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
                  !tautot(iplon, k, ns07+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) = fracrefa(ig,jpl) + fpl                     &
                  &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)
               enddo
            else
               temp   = coldry(n, k) * chi_mls(2,jp(n, k)+1)
               ratco2 = colamt2 / temp
               if (ratco2 > 3.0) then
                  adjfac = 2.0 + (ratco2-2.0)**0.79
                  adjcolco2 = adjfac * temp
               else
                  adjcolco2 = colamt2
               endif

               ind0 = ((jp(n, k)-13)*5 + (jt (n, k)-1)) * nspb(7) + 1
               ind1 = ((jp(n, k)-12)*5 + (jt1(n, k)-1)) * nspb(7) + 1

               indm = indminor(n, k)
               indmp = indm + 1
               ind0p = ind0 + 1
               ind1p = ind1 + 1
               !$acc loop seq
               do ig = 1, ng07
                  ib = ngb(ns07+ig)
                  absco2 = kb_mco2(ig,indm) + minorfrac(n, k)                      &
                  &           * (kb_mco2(ig,indmp) - kb_mco2(ig,indm))

                  taug= colamt3                                 &
                  &            * (fac00(n, k)*absb(ig,ind0) + fac10(n, k)*absb(ig,ind0p)   &
                  &            +  fac01(n, k)*absb(ig,ind1) + fac11(n, k)*absb(ig,ind1p))  &
                  &            + adjcolco2 * absco2
               !  --- ...  empirical modification to code to improve stratospheric cooling rates
               !           for o3.  revised to apply weighting for g-point reduction in this band.
                  if (ig .eq. 6) taug = taug * 0.92
                  if (ig .eq. 7) taug = taug * 0.88
                  if (ig .eq. 8) taug = taug * 1.07
                  if (ig .eq. 9) taug = taug * 1.1
                  if (ig .eq. 10) taug = taug * 0.99
                  if (ig .eq. 11) taug = taug * 0.855
                  !fracs(iplon, k, ns07+ig, jj) = fracrefb(ig)
                  !tautot(iplon, k, ns07+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) = fracrefb(ig)
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)
               enddo
            end if
         end do
      end do

! ..................................
      end subroutine taugb07
! ----------------------------------

! ----------------------------------
      subroutine taugb08 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!     band 8:  1080-1180 cm-1 (low key - h2o; low minor - co2,o3,n2o)  !
!                             (high key - o3; high minor - co2, n2o)   !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb08
      integer, intent(in) :: nlay, laytrop(nxjp_acc_length), ix, myim(fulljj), &
         async_id, ng00, i2, blockjj, fulljj, nxjp_acc_length, jjoffset, &
         max_nxjp_acc_length, jbs_nxjp_acc
      integer, dimension(nxjp_acc_length), intent(in) :: map_jj, map_i

      integer, dimension(max_nxjp_acc_length, nlay), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys) :: temcol, colamt1, colamt2, colamt3, &
                               colamt4, colamt5, colamt6, colamt7, wx3, wx4
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: h2ovmr, o3vmr
      real (kind=kind_phys), dimension(nxptot,nlay), intent(in) :: gasvmr_co2
      real (kind=kind_phys), dimension(nf_vgas), intent(in) :: gasvmr_other
      integer, intent(in) :: nf_vgas
      integer :: k1, rr

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay, nbands), intent(in):: tauaer

      real (kind=kind_phys) :: rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(nxjp_acc_length, nlay, ng00) :: small_fracs, &
     &       small_tautot

!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       indm, indmp, ig

      real (kind=kind_phys) :: tauself, taufor, absco2, abso3, absn2o,  &
     &      ratco2, adjfac, adjcolco2, temp, taug
      integer :: jj, iplon, iplon2, ib, n
      ! GPU: variable name changed: CPU - wx(:,:,3),           GPU - wx3
      ! GPU: variable name changed: CPU - wx(:,:,4),           GPU - wx4
      ! GPU: variable name changed: CPU - colamt(:,:,1),       GPU - colamt1
      ! GPU: variable name changed: CPU - colamt(:,:,2),       GPU - colamt2
      ! GPU: variable name changed: CPU - colamt(:,:,3),       GPU - colamt3
      ! GPU: variable name changed: CPU - colamt(:,:,4),       GPU - colamt4
!
!===> ...  begin here
!
!  --- ...  minor gas mapping level:
!     lower - co2, p = 1053.63 mb, t = 294.2 k
!     lower - o3,  p = 317.348 mb, t = 240.77 k
!     lower - n2o, p = 706.2720 mb, t= 278.94 k
!     lower - cfc12,cfc11
!     upper - co2, p = 35.1632 mb, t = 223.28 k
!     upper - n2o, p = 8.716e-2 mb, t = 226.03 k

      !$acc parallel loop collapse(2) private(ind0, ind1, inds, indf, indm, &
      !$acc&         ind0p, ind1p, indsp, indfp, indmp, temp, ratco2, adjfac, adjcolco2, &
      !$acc&         tauself, taufor, absco2, abso3, absn2o, taug, &
      !$acc&         k1, temcol, colamt1, colamt2, colamt3, colamt4, wx3, wx4, rr) async(async_id)
      do k = 1, nlay
         do n = 1, nxjp_acc_length
            rr = n + jbs_nxjp_acc - 1
            if (ivflip == 0) then
               k1 = nlay + 1 - k
            else
               k1 = k
            end if
            temcol = 1.0e-12 * coldry(n, k)
            colamt1 = max(f_zero,    coldry(n, k)*h2ovmr(n, k))          ! h2o
            colamt2 = max(temcol, coldry(n, k)*gasvmr_co2(rr,k1)) ! co2
            colamt3 = max(temcol, coldry(n, k)*o3vmr(n, k))           ! o3
            if (ilwrgas > 0) then
               colamt4 = max(temcol, coldry(n, k)*gasvmr_other(2))  ! n2o
               wx3 = max( f_zero, coldry(n, k)*gasvmr_other(7) )   ! cf12
               wx4 = max( f_zero, coldry(n, k)*gasvmr_other(8) )   ! cf22
            else
               colamt4 = f_zero     ! n2o
               wx3 = f_zero
               wx4 = f_zero
            endif

            if (k .le. laytrop(n)) then
               !  --- ...  lower atmosphere loop
               ind0 = ((jp(n, k)-1)*5 + (jt (n, k)-1)) * nspa(8) + 1
               ind1 = ( jp(n, k)   *5 + (jt1(n, k)-1)) * nspa(8) + 1

               inds = indself(n, k)
               indf = indfor(n, k)
               indm = indminor(n, k)
               ind0p = ind0 + 1
               ind1p = ind1 + 1
               indsp = inds + 1
               indfp = indf + 1
               indmp = indm + 1

               !  --- ...  in atmospheres where the amount of co2 is too great to be considered
               !           a minor species, adjust the column amount of co2 by an empirical factor
               !           to obtain the proper contribution.

               temp   = coldry(n, k) * chi_mls(2,jp(n, k)+1)
               ratco2 = colamt2 / temp
               if (ratco2 > 3.0) then
                  adjfac = 2.0 + (ratco2-2.0)**0.65
                  adjcolco2 = adjfac * temp
               else
                  adjcolco2 = colamt2
               endif
               !$acc loop seq
               do ig = 1, ng08
                  ib = ngb(ns08+ig)
                  tauself = selffac(n, k) * (selfref(ig,inds) + selffrac(n, k)        &
                  &            * (selfref(ig,indsp) - selfref(ig,inds)))
                  taufor  = forfac(n, k) * (forref(ig,indf) + forfrac(n, k)           &
                  &            * (forref(ig,indfp) - forref(ig,indf)))
                  absco2  = (ka_mco2(ig,indm) + minorfrac(n, k)                    &
                  &            * (ka_mco2(ig,indmp) - ka_mco2(ig,indm)))
                  abso3   = (ka_mo3(ig,indm) + minorfrac(n, k)                     &
                  &            * (ka_mo3(ig,indmp) - ka_mo3(ig,indm)))
                  absn2o  = (ka_mn2o(ig,indm) + minorfrac(n, k)                    &
                  &            * (ka_mn2o(ig,indmp) - ka_mn2o(ig,indm)))

                  taug = colamt1                                 &
                  &            * (fac00(n, k)*absa(ig,ind0) + fac10(n, k)*absa(ig,ind0p)   &
                  &            +  fac01(n, k)*absa(ig,ind1) + fac11(n, k)*absa(ig,ind1p))  &
                  &            + tauself+taufor + adjcolco2*absco2                   &
                  &            + colamt3*abso3 + colamt4*absn2o              &
                  &            + wx3*cfc12(ig) + wx4*cfc22adj(ig)

                  !fracs(iplon, k, ns08+ig, jj) = fracrefa(ig)
                  !tautot(iplon, k, ns08+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) = fracrefa(ig)
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)
               enddo
            else
               ind0 = ((jp(n, k)-13)*5 + (jt (n, k)-1)) * nspb(8) + 1
               ind1 = ((jp(n, k)-12)*5 + (jt1(n, k)-1)) * nspb(8) + 1

               indm = indminor(n, k)
               ind0p = ind0 + 1
               ind1p = ind1 + 1
               indmp = indm + 1

               !  --- ...  in atmospheres where the amount of co2 is too great to be considered
               !           a minor species, adjust the column amount of co2 by an empirical factor
               !           to obtain the proper contribution.

               temp   = coldry(n, k) * chi_mls(2,jp(n, k)+1)
               ratco2 = colamt2 / temp
               if (ratco2 > 3.0) then
                  adjfac = 2.0 + (ratco2-2.0)**0.65
                  adjcolco2 = adjfac * temp
               else
                  adjcolco2 = colamt2
               endif
               !$acc loop seq
               do ig = 1, ng08
                  ib = ngb(ns08+ig)
                  absco2 = (kb_mco2(ig,indm) + minorfrac(n, k)                     &
                  &           * (kb_mco2(ig,indmp) - kb_mco2(ig,indm)))
                  absn2o = (kb_mn2o(ig,indm) + minorfrac(n, k)                     &
                  &           * (kb_mn2o(ig,indmp) - kb_mn2o(ig,indm)))

                  taug = colamt3                                 &
                  &           * (fac00(n, k)*absb(ig,ind0) + fac10(n, k)*absb(ig,ind0p)    &
                  &           +  fac01(n, k)*absb(ig,ind1) + fac11(n, k)*absb(ig,ind1p))   &
                  &           + adjcolco2*absco2 + colamt4*absn2o                &
                  &           + wx3*cfc12(ig) + wx4*cfc22adj(ig)

                  !fracs(iplon, k, ns08+ig, jj) = fracrefb(ig)
                  !tautot(iplon, k, ns08+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) = fracrefb(ig)
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)
               enddo
            end if
         end do
      end do

! ..................................
      end subroutine taugb08
! ----------------------------------

! ----------------------------------
      subroutine taugb09 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!     band 9:  1180-1390 cm-1 (low key - h2o,ch4; low minor - n2o)     !
!                             (high key - ch4; high minor - n2o)       !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb09
      integer, intent(in) :: nlay, laytrop(nxjp_acc_length), ix, myim(fulljj), &
         async_id, ng00, i2, blockjj, fulljj, nxjp_acc_length, jjoffset, &
         max_nxjp_acc_length, jbs_nxjp_acc
      integer, dimension(nxjp_acc_length), intent(in) :: map_jj, map_i

      integer, dimension(max_nxjp_acc_length, nlay), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys) :: temcol, colamt1, colamt2, colamt3, &
                               colamt4, colamt5, colamt6, colamt7
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: h2ovmr, o3vmr
      real (kind=kind_phys), dimension(nxptot,nlay), intent(in) :: gasvmr_co2
      real (kind=kind_phys), dimension(nf_vgas), intent(in) :: gasvmr_other
      integer, intent(in) :: nf_vgas
      integer :: k1, rr

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay, nbands), intent(in):: tauaer

      real (kind=kind_phys) :: rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(nxjp_acc_length, nlay, ng00) :: small_fracs, &
     &       small_tautot

!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       id000, id010, id100, id110, id200, id210, indm, indmp,     &
     &       id001, id011, id101, id111, id201, id211, jmn2o, jmn2op,   &
     &       jpl, jplp, ig, js, js1, jpr

      real (kind=kind_phys) :: tauself, taufor, n2om1, n2om2, absn2o,   &
     &       speccomb,       specparm,       specmult,       fs,        &
     &       speccomb1,      specparm1,      specmult1,      fs1,       &
     &       speccomb_mn2o,  specparm_mn2o,  specmult_mn2o,  fmn2o,     &
     &       speccomb_planck,specparm_planck,specmult_planck,fpl,       &
     &       refrat_planck_a, refrat_m_a, ratn2o, adjfac, adjcoln2o,    &
     &       fac000, fac100, fac200, fac010, fac110, fac210,            &
     &       fac001, fac101, fac201, fac011, fac111, fac211,            &
     &       p0, p40, fk00, fk10, fk20, p1, p41, fk01, fk11, fk21, temp, taug
      integer :: jj, iplon, iplon2, ib, n
      ! GPU: variable name changed: CPU - colamt(:,:,1),       GPU - colamt1
      ! GPU: variable name changed: CPU - colamt(:,:,4),       GPU - colamt4
      ! GPU: variable name changed: CPU - colamt(:,:,5),       GPU - colamt5
!
!===> ...  begin here
!
!  --- ...  minor gas mapping level :
!     lower - n2o, p = 706.272 mbar, t = 278.94 k
!     upper - n2o, p = 95.58 mbar, t = 215.7 k

!  --- ...  calculate reference ratio to be used in calculation of planck
!           fraction in lower/upper atmosphere.
      refrat_planck_a = chi_mls(1,9)/chi_mls(6,9)       ! p = 212 mb
      refrat_m_a = chi_mls(1,3)/chi_mls(6,3)            ! p = 706.272 mb
      !$acc parallel loop collapse(2) private(speccomb, specparm, specmult, js, &
      !$acc&         fs, ind0, speccomb1, specparm1, specmult1, js1, fs1, ind1, &
      !$acc&         speccomb_planck, specparm_planck, specmult_planck, jpl, &
      !$acc&         fpl, inds, indf, indm, indsp, indfp, indmp, jplp, &
      !$acc&         id000, id010, speccomb_mn2o, specparm_mn2o, specmult_mn2o, &
      !$acc&         id100, id110, id200, id210, fac000, fac100, fac200, fac010, &
      !$acc&         fac110, fac210, id001, id011, id101, id111, id201, id211, &
      !$acc&         fac001, fac101, fac201, fac011, fac111, fac211, tauself, &
      !$acc&         taufor, p0, p40, fk00, fk10, fk20, p1, jmn2o, fmn2o, &
      !$acc&         p41, fk01, fk11, fk21, temp, ratn2o, adjfac, &
      !$acc&         adjcoln2o, n2om1, n2om2, absn2o, taug, jpr, rfrate, &
      !$acc&         k1, temcol, colamt1, colamt4, colamt5, rr) async(async_id)
      do k = 1, nlay
         do n = 1, nxjp_acc_length
            rr = n + jbs_nxjp_acc - 1
            jpr = jp(n, k)
            temcol = 1.0e-12 * coldry(n, k)
            colamt1 = max(f_zero,    coldry(n, k)*h2ovmr(n, k))          ! h2o
            if (ilwrgas > 0) then
               colamt4 = max(temcol, coldry(n, k)*gasvmr_other(2))  ! n2o
               colamt5 = max(temcol, coldry(n, k)*gasvmr_other(3))  ! ch4
            else
               colamt4 = f_zero     ! n2o
               colamt5 = f_zero     ! ch4
            endif
            if (k .le. laytrop(n)) then
               !  --- ...  lower atmosphere loop
               rfrate = chi_mls(1,jpr) / chi_mls(6,jpr)
               speccomb = colamt1 + rfrate*colamt5
               specparm = colamt1 / speccomb
               specmult = 8.0 * min(specparm, oneminus)
               js = 1 + int(specmult)
               fs = mod(specmult, f_one)
               ind0 = ((jp(n, k)-1)*5 + (jt(n, k)-1)) * nspa(9) + js

               rfrate = chi_mls(1,jpr+1) / chi_mls(6,jpr+1)
               speccomb1 = colamt1 + rfrate*colamt5
               specparm1 = colamt1 / speccomb1
               specmult1 = 8.0 * min(specparm1, oneminus)
               js1 = 1 + int(specmult1)
               fs1 = mod(specmult1, f_one)
               ind1 = (jp(n, k)*5 + (jt1(n, k)-1)) * nspa(9) + js1

               speccomb_mn2o = colamt1 + refrat_m_a*colamt5
               specparm_mn2o = colamt1 / speccomb_mn2o
               specmult_mn2o = 8.0 * min(specparm_mn2o, oneminus)
               jmn2o = 1 + int(specmult_mn2o)
               fmn2o = mod(specmult_mn2o, f_one)

               speccomb_planck = colamt1 + refrat_planck_a*colamt5
               specparm_planck = colamt1 / speccomb_planck
               specmult_planck = 8.0 * min(specparm_planck, oneminus)
               jpl = 1 + int(specmult_planck)
               fpl = mod(specmult_planck, f_one)

               inds = indself(n, k)
               indf = indfor(n, k)
               indm = indminor(n, k)
               indsp = inds + 1
               indfp = indf + 1
               indmp = indm + 1
               jplp  = jpl  + 1
               jmn2op= jmn2o+ 1

               !  --- ...  in atmospheres where the amount of n2o is too great to be considered
               !           a minor species, adjust the column amount of n2o by an empirical factor
               !           to obtain the proper contribution.

               temp   = coldry(n, k) * chi_mls(4,jp(n, k)+1)
               ratn2o = colamt4 / temp
               if (ratn2o > 1.5) then
                  adjfac = 0.5 + (ratn2o-0.5)**0.65
                  adjcoln2o = adjfac * temp
               else
                  adjcoln2o = colamt4
               endif

               if (specparm < 0.125 .and. specparm1 < 0.125) then
                  p0 = fs - f_one
                  p40 = p0**4
                  fk00 = p40
                  fk10 = f_one - p0 - 2.0*p40
                  fk20 = p0 + p40

                  p1 = fs1 - f_one
                  p41 = p1**4
                  fk01 = p41
                  fk11 = f_one - p1 - 2.0*p41
                  fk21 = p1 + p41

                  id000 = ind0
                  id010 = ind0 + 9
                  id100 = ind0 + 1
                  id110 = ind0 +10
                  id200 = ind0 + 2
                  id210 = ind0 +11

                  id001 = ind1
                  id011 = ind1 + 9
                  id101 = ind1 + 1
                  id111 = ind1 +10
                  id201 = ind1 + 2
                  id211 = ind1 +11

               elseif (specparm > 0.875 .and. specparm1 > 0.875) then
                  p0 = -fs
                  p40 = p0**4
                  fk00 = p40
                  fk10 = f_one - p0 - 2.0*p40
                  fk20 = p0 + p40

                  p1 = -fs1
                  p41 = p1**4
                  fk01 = p41
                  fk11 = f_one - p1 - 2.0*p41
                  fk21 = p1 + p41

                  id000 = ind0 + 1
                  id010 = ind0 +10
                  id100 = ind0
                  id110 = ind0 + 9
                  id200 = ind0 - 1
                  id210 = ind0 + 8

                  id001 = ind1 + 1
                  id011 = ind1 +10
                  id101 = ind1
                  id111 = ind1 + 9
                  id201 = ind1 - 1
                  id211 = ind1 + 8
               else
                  fk00 = f_one - fs
                  fk10 = fs
                  fk20 = f_zero

                  fk01 = f_one - fs1
                  fk11 = fs1
                  fk21 = f_zero

                  id000 = ind0
                  id010 = ind0 + 9
                  id100 = ind0 + 1
                  id110 = ind0 +10
                  id200 = ind0
                  id210 = ind0

                  id001 = ind1
                  id011 = ind1 + 9
                  id101 = ind1 + 1
                  id111 = ind1 +10
                  id201 = ind1
                  id211 = ind1
               endif

               fac000 = fk00 * fac00(n, k)
               fac100 = fk10 * fac00(n, k)
               fac200 = fk20 * fac00(n, k)
               fac010 = fk00 * fac10(n, k)
               fac110 = fk10 * fac10(n, k)
               fac210 = fk20 * fac10(n, k)

               fac001 = fk01 * fac01(n, k)
               fac101 = fk11 * fac01(n, k)
               fac201 = fk21 * fac01(n, k)
               fac011 = fk01 * fac11(n, k)
               fac111 = fk11 * fac11(n, k)
               fac211 = fk21 * fac11(n, k)
               !$acc loop seq
               do ig = 1, ng09
                  ib = ngb(ns09+ig)
                  tauself = selffac(n, k)* (selfref(ig,inds) + selffrac(n, k)         &
                  &            * (selfref(ig,indsp) - selfref(ig,inds)))
                  taufor  = forfac(n, k) * (forref(ig,indf) + forfrac(n, k)           &
                  &            * (forref(ig,indfp) - forref(ig,indf))) 
                  n2om1   = ka_mn2o(ig,jmn2o,indm) + fmn2o                      &
                  &            * (ka_mn2o(ig,jmn2op,indm) - ka_mn2o(ig,jmn2o,indm))
                  n2om2   = ka_mn2o(ig,jmn2o,indmp) + fmn2o                     &
                  &            * (ka_mn2o(ig,jmn2op,indmp) - ka_mn2o(ig,jmn2o,indmp))
                  absn2o  = n2om1 + minorfrac(n, k) * (n2om2 - n2om1)

                  taug = speccomb                                    &
                  &                * (fac000*absa(ig,id000) + fac010*absa(ig,id010)  &
                  &                +  fac100*absa(ig,id100) + fac110*absa(ig,id110)  &
                  &                +  fac200*absa(ig,id200) + fac210*absa(ig,id210)) &
                  &                +     speccomb1                                   &
                  &                * (fac001*absa(ig,id001) + fac011*absa(ig,id011)  &
                  &                +  fac101*absa(ig,id101) + fac111*absa(ig,id111)  &
                  &                +  fac201*absa(ig,id201) + fac211*absa(ig,id211)) &
                  &                + tauself + taufor + adjcoln2o*absn2o            

                  !fracs(iplon, k, ns09+ig, jj) = fracrefa(ig,jpl) + fpl                     &
                  !&                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
                  !tautot(iplon, k, ns09+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) = fracrefa(ig,jpl) + fpl                     &
                  &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)
               enddo
            else
               ind0 = ((jp(n, k)-13)*5 + (jt (n, k)-1)) * nspb(9) + 1
               ind1 = ((jp(n, k)-12)*5 + (jt1(n, k)-1)) * nspb(9) + 1

               indm = indminor(n, k)
               ind0p = ind0 + 1
               ind1p = ind1 + 1
               indmp = indm + 1

               !  --- ...  in atmospheres where the amount of n2o is too great to be considered
               !           a minor species, adjust the column amount of n2o by an empirical factor
               !           to obtain the proper contribution.

               temp   = coldry(n, k) * chi_mls(4,jp(n, k)+1)
               ratn2o = colamt4 / temp
               if (ratn2o > 1.5) then
                  adjfac = 0.5 + (ratn2o - 0.5)**0.65
                  adjcoln2o = adjfac * temp
               else
                  adjcoln2o = colamt4
               endif
               !$acc loop seq
               do ig = 1, ng09
                  ib = ngb(ns09+ig)
                  absn2o = kb_mn2o(ig,indm) + minorfrac(n, k)                      &
                  &           * (kb_mn2o(ig,indmp) - kb_mn2o(ig,indm))

                  taug = colamt5                                 &
                  &           * (fac00(n, k)*absb(ig,ind0) + fac10(n, k)*absb(ig,ind0p)    &
                  &           +  fac01(n, k)*absb(ig,ind1) + fac11(n, k)*absb(ig,ind1p))   &
                  &           + adjcoln2o*absn2o

                  !fracs(iplon, k, ns09+ig, jj) = fracrefb(ig)
                  !tautot(iplon, k, ns09+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) = fracrefb(ig)
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)
               enddo
            end if
         end do
      end do

! ..................................
      end subroutine taugb09
! ----------------------------------

! ----------------------------------
      subroutine taugb10 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!     band 10:  1390-1480 cm-1 (low key - h2o; high key - h2o)         !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb10
      integer, intent(in) :: nlay, laytrop(nxjp_acc_length), ix, myim(fulljj), &
         async_id, ng00, i2, blockjj, fulljj, nxjp_acc_length, jjoffset, &
         max_nxjp_acc_length, jbs_nxjp_acc
      integer, dimension(nxjp_acc_length), intent(in) :: map_jj, map_i

      integer, dimension(max_nxjp_acc_length, nlay), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys) :: temcol, colamt1, colamt2, colamt3, &
                               colamt4, colamt5, colamt6, colamt7
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: h2ovmr, o3vmr
      real (kind=kind_phys), dimension(nxptot,nlay), intent(in) :: gasvmr_co2
      real (kind=kind_phys), dimension(nf_vgas), intent(in) :: gasvmr_other
      integer, intent(in) :: nf_vgas
      integer :: k1, rr

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay, nbands), intent(in):: tauaer

      real (kind=kind_phys) :: rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(nxjp_acc_length, nlay, ng00) :: small_fracs, &
     &       small_tautot

!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       ig

      real (kind=kind_phys) :: tauself, taufor, taug
      integer :: jj, iplon, iplon2, ib, n
      ! GPU: variable name changed: CPU - colamt(:,:,1),       GPU - colamt1
!
!===> ...  begin here
!
      !$acc parallel loop collapse(2) private(ind0, ind1, inds, indf, ind0p, &
      !$acc&         ind1p, indsp, indfp, tauself, taufor, taug, colamt1, rr) async(async_id)
      do k = 1, nlay
         do n = 1, nxjp_acc_length
            rr = n + jbs_nxjp_acc - 1
            colamt1 = max(f_zero,    coldry(n, k)*h2ovmr(n, k))          ! h2o
            if (k .le. laytrop(n)) then
               !  --- ...  lower atmosphere loop
               ind0 = ((jp(n, k)-1)*5 + (jt (n, k)-1)) * nspa(10) + 1
               ind1 = ( jp(n, k)   *5 + (jt1(n, k)-1)) * nspa(10) + 1

               inds = indself(n, k)
               indf = indfor(n, k)
               ind0p = ind0 + 1
               ind1p = ind1 + 1
               indsp = inds + 1
               indfp = indf + 1
               !$acc loop seq
               do ig = 1, ng10
                  ib = ngb(ns10+ig)
                  tauself = selffac(n, k) * (selfref(ig,inds) + selffrac(n, k)        &
                  &            * (selfref(ig,indsp) - selfref(ig,inds)))
                  taufor  = forfac(n, k) * (forref(ig,indf) + forfrac(n, k)           &
                  &            * (forref(ig,indfp) - forref(ig,indf))) 

                  taug = colamt1                                 &
                  &            * (fac00(n, k)*absa(ig,ind0) + fac10(n, k)*absa(ig,ind0p)   &
                  &            +  fac01(n, k)*absa(ig,ind1) + fac11(n, k)*absa(ig,ind1p))  &
                  &            + tauself + taufor

                  !fracs(iplon, k, ns10+ig, jj) = fracrefa(ig)
                  !tautot(iplon, k, ns10+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) = fracrefa(ig)
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)
               enddo
            else
               ind0 = ((jp(n, k)-13)*5 + (jt (n, k)-1)) * nspb(10) + 1
               ind1 = ((jp(n, k)-12)*5 + (jt1(n, k)-1)) * nspb(10) + 1

               indf = indfor(n, k)
               ind0p = ind0 + 1
               ind1p = ind1 + 1
               indfp = indf + 1
               !$acc loop seq
               do ig = 1, ng10
                  ib = ngb(ns10+ig)
                  taufor = forfac(n, k) * (forref(ig,indf) + forfrac(n, k)            &
                  &           * (forref(ig,indfp) - forref(ig,indf))) 

                  taug = colamt1                                 &
                  &            * (fac00(n, k)*absb(ig,ind0) + fac10(n, k)*absb(ig,ind0p)   &
                  &            +  fac01(n, k)*absb(ig,ind1) + fac11(n, k)*absb(ig,ind1p))  &
                  &            + taufor

                  !fracs(iplon, k, ns10+ig, jj) = fracrefb(ig)
                  !tautot(iplon, k, ns10+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) = fracrefb(ig)
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)
               enddo
            end if
         end do
      end do

! ..................................
      end subroutine taugb10
! ----------------------------------

! ----------------------------------
      subroutine taugb11 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!     band 11:  1480-1800 cm-1 (low - h2o; low minor - o2)             !
!                              (high key - h2o; high minor - o2)       !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb11
      integer, intent(in) :: nlay, laytrop(nxjp_acc_length), ix, myim(fulljj), &
         async_id, ng00, i2, blockjj, fulljj, nxjp_acc_length, jjoffset, &
         max_nxjp_acc_length, jbs_nxjp_acc
      integer, dimension(nxjp_acc_length), intent(in) :: map_jj, map_i

      integer, dimension(max_nxjp_acc_length, nlay), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys) :: temcol, colamt1, colamt2, colamt3, &
                               colamt4, colamt5, colamt6, colamt7
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: h2ovmr, o3vmr
      real (kind=kind_phys), dimension(nxptot,nlay), intent(in) :: gasvmr_co2
      real (kind=kind_phys), dimension(nf_vgas), intent(in) :: gasvmr_other
      integer, intent(in) :: nf_vgas
      integer :: k1, rr

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay, nbands), intent(in):: tauaer

      real (kind=kind_phys) :: rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(nxjp_acc_length, nlay, ng00) :: small_fracs, &
     &       small_tautot

!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       indm, indmp, ig

      real (kind=kind_phys) :: scaleo2, tauself, taufor, tauo2, taug
      integer :: jj, iplon, iplon2, ib, n
      ! GPU: variable name changed: CPU - colamt(:,:,1),       GPU - colamt1
      ! GPU: variable name changed: CPU - colamt(:,:,6),       GPU - colamt6
!
!===> ...  begin here
!
!  --- ...  minor gas mapping level :
!     lower - o2, p = 706.2720 mbar, t = 278.94 k
!     upper - o2, p = 4.758820 mbarm t = 250.85 k

      !$acc parallel loop collapse(2) private(ind0, ind1, inds, indf, ind0p, &
      !$acc&         ind1p, indsp, indfp, tauself, taufor, indm, indmp, scaleo2, &
      !$acc&         tauo2, taug, colamt1, colamt6, rr) async(async_id)
      do k = 1, nlay
         do n = 1, nxjp_acc_length
            rr = n + jbs_nxjp_acc - 1
            temcol = 1.0e-12 * coldry(n, k)
            colamt1 = max(f_zero,    coldry(n, k)*h2ovmr(n, k))          ! h2o
            if (ilwrgas > 0) then
               colamt6 = max(f_zero,    coldry(n, k)*gasvmr_other(4))  ! o2
            else
               colamt6 = f_zero     ! o2
            endif
            if (k .le. laytrop(n)) then
               !  --- ...  lower atmosphere loop
               ind0 = ((jp(n, k)-1)*5 + (jt (n, k)-1)) * nspa(11) + 1
               ind1 = ( jp(n, k)   *5 + (jt1(n, k)-1)) * nspa(11) + 1

               inds = indself(n, k)
               indf = indfor(n, k)
               indm = indminor(n, k)
               ind0p = ind0 + 1
               ind1p = ind1 + 1
               indsp = inds + 1
               indfp = indf + 1
               indmp = indm + 1

               scaleo2 = colamt6 * scaleminor(n, k)
               !$acc loop seq
               do ig = 1, ng11
                  ib = ngb(ns11+ig)
                  tauself = selffac(n, k) * (selfref(ig,inds) + selffrac(n, k)        &
                  &            * (selfref(ig,indsp) - selfref(ig,inds)))
                  taufor  = forfac(n, k) * (forref(ig,indf) + forfrac(n, k)           &
                  &            * (forref(ig,indfp) - forref(ig,indf)))
                  tauo2   = scaleo2 * (ka_mo2(ig,indm) + minorfrac(n, k)           &
                  &            * (ka_mo2(ig,indmp) - ka_mo2(ig,indm)))

                  taug = colamt1                                 &
                  &            * (fac00(n, k)*absa(ig,ind0) + fac10(n, k)*absa(ig,ind0p)   &
                  &            +  fac01(n, k)*absa(ig,ind1) + fac11(n, k)*absa(ig,ind1p))  &
                  &            + tauself + taufor + tauo2

                  !fracs(iplon, k, ns11+ig, jj) = fracrefa(ig)
                  !tautot(iplon, k, ns11+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) = fracrefa(ig)
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)
               enddo
            else
               ind0 = ((jp(n, k)-13)*5 + (jt (n, k)-1)) * nspb(11) + 1
               ind1 = ((jp(n, k)-12)*5 + (jt1(n, k)-1)) * nspb(11) + 1

               indf = indfor(n, k)
               indm = indminor(n, k)
               ind0p = ind0 + 1
               ind1p = ind1 + 1
               indfp = indf + 1
               indmp = indm + 1

               scaleo2 = colamt6 * scaleminor(n, k)
               !$acc loop seq
               do ig = 1, ng11
                  ib = ngb(ns11+ig)
                  taufor = forfac(n, k) * (forref(ig,indf) + forfrac(n, k)            &
                  &           * (forref(ig,indfp) - forref(ig,indf))) 
                  tauo2  = scaleo2 * (kb_mo2(ig,indm) + minorfrac(n, k)            &
                  &           * (kb_mo2(ig,indmp) - kb_mo2(ig,indm)))

                  taug = colamt1                                 &
                  &            * (fac00(n, k)*absb(ig,ind0) + fac10(n, k)*absb(ig,ind0p)   &
                  &            +  fac01(n, k)*absb(ig,ind1) + fac11(n, k)*absb(ig,ind1p))  &
                  &            + taufor + tauo2

                  !fracs(iplon, k, ns11+ig, jj) = fracrefb(ig)
                  !tautot(iplon, k, ns11+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) = fracrefb(ig)
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)
               enddo
            end if
         end do
      end do

! ..................................
      end subroutine taugb11
! ----------------------------------

! ----------------------------------
      subroutine taugb12 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!     band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)         !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb12
      integer, intent(in) :: nlay, laytrop(nxjp_acc_length), ix, myim(fulljj), &
         async_id, ng00, i2, blockjj, fulljj, nxjp_acc_length, jjoffset, &
         max_nxjp_acc_length, jbs_nxjp_acc
      integer, dimension(nxjp_acc_length), intent(in) :: map_jj, map_i

      integer, dimension(max_nxjp_acc_length, nlay), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys) :: temcol, colamt1, colamt2, colamt3, &
                               colamt4, colamt5, colamt6, colamt7
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: h2ovmr, o3vmr
      real (kind=kind_phys), dimension(nxptot,nlay), intent(in) :: gasvmr_co2
      real (kind=kind_phys), dimension(nf_vgas), intent(in) :: gasvmr_other
      integer, intent(in) :: nf_vgas
      integer :: k1, rr

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay, nbands), intent(in):: tauaer

      real (kind=kind_phys) :: rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(nxjp_acc_length, nlay, ng00) :: small_fracs, &
     &       small_tautot

!  ---  locals:
      integer :: k, ind0, ind1, inds, indsp, indf, indfp, jpl, jplp,    &
     &       id000, id010, id100, id110, id200, id210, ig, js, js1,     &
     &       id001, id011, id101, id111, id201, id211, jpr

      real (kind=kind_phys) :: tauself, taufor, refrat_planck_a,        &
     &       speccomb,       specparm,       specmult,       fs,        &
     &       speccomb1,      specparm1,      specmult1,      fs1,       &
     &       speccomb_planck,specparm_planck,specmult_planck,fpl,       &
     &       fac000, fac100, fac200, fac010, fac110, fac210,            &
     &       fac001, fac101, fac201, fac011, fac111, fac211,            &
     &       p0, p40, fk00, fk10, fk20, p1, p41, fk01, fk11, fk21, taug
      integer :: jj, iplon, iplon2, ib, n
      ! GPU: variable name changed: CPU - colamt(:,:,1),       GPU - colamt1
      ! GPU: variable name changed: CPU - colamt(:,:,2),       GPU - colamt2
!
!===> ...  begin here
!
!  --- ...  calculate reference ratio to be used in calculation of planck
!           fraction in lower/upper atmosphere.
      refrat_planck_a = chi_mls(1,10)/chi_mls(2,10)      ! p =   174.164 mb
      !$acc parallel loop collapse(2) private(speccomb, specparm, specmult, js, &
      !$acc&         fs, ind0, speccomb1, specparm1, specmult1, js1, fs1, ind1, &
      !$acc&         speccomb_planck, specparm_planck, specmult_planck, jpl, &
      !$acc&         fpl, inds, indf, indsp, indfp, jplp, &
      !$acc&         id000, id010, k1, temcol, colamt1, colamt2, &
      !$acc&         id100, id110, id200, id210, fac000, fac100, fac200, fac010, &
      !$acc&         fac110, fac210, id001, id011, id101, id111, id201, id211, &
      !$acc&         fac001, fac101, fac201, fac011, fac111, fac211, tauself, &
      !$acc&         taufor, p0, p40, fk00, fk10, fk20, p1, &
      !$acc&         p41, fk01, fk11, fk21, taug, jpr, rfrate, rr) async(async_id)
      do k = 1, nlay
         do n = 1, nxjp_acc_length
            jpr = jp(n, k)
            rr = n + jbs_nxjp_acc - 1
            if (ivflip == 0) then
               k1 = nlay + 1 - k
            else
               k1 = k
            end if
            temcol = 1.0e-12 * coldry(n, k)
            colamt1 = max(f_zero,    coldry(n, k)*h2ovmr(n, k))          ! h2o
            colamt2 = max(temcol, coldry(n, k)*gasvmr_co2(rr,k1)) ! co2
            if (k .le. laytrop(n)) then
               !  --- ...  lower atmosphere loop
               rfrate = chi_mls(1,jpr) / chi_mls(2,jpr)
               speccomb = colamt1 + rfrate*colamt2
               specparm = colamt1 / speccomb
               specmult = 8.0 * min(specparm, oneminus)
               js = 1 + int(specmult)
               fs = mod(specmult, f_one)
               ind0 = ((jp(n, k)-1)*5 + (jt(n, k)-1)) * nspa(12) + js

               rfrate = chi_mls(1,jpr+1) / chi_mls(2,jpr+1)
               speccomb1 = colamt1 + rfrate*colamt2
               specparm1 = colamt1 / speccomb1
               specmult1 = 8.0 * min(specparm1, oneminus)
               js1 = 1 + int(specmult1)
               fs1 = mod(specmult1, f_one)
               ind1 = (jp(n, k)*5 + (jt1(n, k)-1)) * nspa(12) + js1

               speccomb_planck = colamt1 + refrat_planck_a*colamt2
               specparm_planck = colamt1 / speccomb_planck
               if (specparm_planck >= oneminus) specparm_planck=oneminus
               specmult_planck = 8.0 * specparm_planck
               jpl = 1 + int(specmult_planck)
               fpl = mod(specmult_planck, f_one)

               inds = indself(n, k)
               indf = indfor(n, k)
               indsp = inds + 1
               indfp = indf + 1
               jplp  = jpl  + 1

               if (specparm < 0.125 .and. specparm1 < 0.125) then
                  p0 = fs - f_one
                  p40 = p0**4
                  fk00 = p40
                  fk10 = f_one - p0 - 2.0*p40
                  fk20 = p0 + p40

                  p1 = fs1 - f_one
                  p41 = p1**4
                  fk01 = p41
                  fk11 = f_one - p1 - 2.0*p41
                  fk21 = p1 + p41

                  id000 = ind0
                  id010 = ind0 + 9
                  id100 = ind0 + 1
                  id110 = ind0 +10
                  id200 = ind0 + 2
                  id210 = ind0 +11

                  id001 = ind1
                  id011 = ind1 + 9
                  id101 = ind1 + 1
                  id111 = ind1 +10
                  id201 = ind1 + 2
                  id211 = ind1 +11
               elseif (specparm > 0.875 .and. specparm1 > 0.875) then
                  p0 = -fs
                  p40 = p0**4
                  fk00 = p40
                  fk10 = f_one - p0 - 2.0*p40
                  fk20 = p0 + p40

                  p1 = -fs1
                  p41 = p1**4
                  fk01 = p41
                  fk11 = f_one - p1 - 2.0*p41
                  fk21 = p1 + p41

                  id000 = ind0 + 1
                  id010 = ind0 +10
                  id100 = ind0
                  id110 = ind0 + 9
                  id200 = ind0 - 1
                  id210 = ind0 + 8

                  id001 = ind1 + 1
                  id011 = ind1 +10
                  id101 = ind1
                  id111 = ind1 + 9
                  id201 = ind1 - 1
                  id211 = ind1 + 8
               else
                  fk00 = f_one - fs
                  fk10 = fs
                  fk20 = f_zero

                  fk01 = f_one - fs1
                  fk11 = fs1
                  fk21 = f_zero

                  id000 = ind0
                  id010 = ind0 + 9
                  id100 = ind0 + 1
                  id110 = ind0 +10
                  id200 = ind0
                  id210 = ind0

                  id001 = ind1
                  id011 = ind1 + 9
                  id101 = ind1 + 1
                  id111 = ind1 +10
                  id201 = ind1
                  id211 = ind1
               endif

               fac000 = fk00 * fac00(n, k)
               fac100 = fk10 * fac00(n, k)
               fac200 = fk20 * fac00(n, k)
               fac010 = fk00 * fac10(n, k)
               fac110 = fk10 * fac10(n, k)
               fac210 = fk20 * fac10(n, k)

               fac001 = fk01 * fac01(n, k)
               fac101 = fk11 * fac01(n, k)
               fac201 = fk21 * fac01(n, k)
               fac011 = fk01 * fac11(n, k)
               fac111 = fk11 * fac11(n, k)
               fac211 = fk21 * fac11(n, k)
               !$acc loop seq
               do ig = 1, ng12
                  ib = ngb(ns12+ig)
                  tauself = selffac(n, k)* (selfref(ig,inds) + selffrac(n, k)         &
                  &            * (selfref(ig,indsp) - selfref(ig,inds)))
                  taufor  = forfac(n, k) * (forref(ig,indf) + forfrac(n, k)           &
                  &            * (forref(ig,indfp) - forref(ig,indf))) 

                  taug = speccomb                                    &
                  &                * (fac000*absa(ig,id000) + fac010*absa(ig,id010)  &
                  &                +  fac100*absa(ig,id100) + fac110*absa(ig,id110)  &
                  &                +  fac200*absa(ig,id200) + fac210*absa(ig,id210)) &
                  &                +     speccomb1                                   &
                  &                * (fac001*absa(ig,id001) + fac011*absa(ig,id011)  &
                  &                +  fac101*absa(ig,id101) + fac111*absa(ig,id111)  &
                  &                +  fac201*absa(ig,id201) + fac211*absa(ig,id211)) &
                  &                + tauself + taufor

                  !fracs(iplon, k, ns12+ig, jj) = fracrefa(ig,jpl) + fpl                     &
                  !&                     *(fracrefa(ig,jplp) - fracrefa(ig,jpl))
                  !tautot(iplon, k, ns12+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) = fracrefa(ig,jpl) + fpl                     &
                  &                     *(fracrefa(ig,jplp) - fracrefa(ig,jpl))
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)
               enddo
            else
               !$acc loop seq
               do ig = 1, ng12
                  ib = ngb(ns12+ig)
                  taug = f_zero
                  !fracs(iplon, k, ns12+ig, jj) = f_zero
                  !tautot(iplon, k, ns12+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) = f_zero
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)
               enddo
            end if
         end do
      end do

! ..................................
      end subroutine taugb12
! ----------------------------------

! ----------------------------------
      subroutine taugb13 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!     band 13:  2080-2250 cm-1 (low key-h2o,n2o; high minor-o3 minor)  !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb13
      integer, intent(in) :: nlay, laytrop(nxjp_acc_length), ix, myim(fulljj), &
         async_id, ng00, i2, blockjj, fulljj, nxjp_acc_length, jjoffset, &
         max_nxjp_acc_length, jbs_nxjp_acc
      integer, dimension(nxjp_acc_length), intent(in) :: map_jj, map_i

      integer, dimension(max_nxjp_acc_length, nlay), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys) :: temcol, colamt1, colamt2, colamt3, &
                               colamt4, colamt5, colamt6, colamt7
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: h2ovmr, o3vmr
      real (kind=kind_phys), dimension(nxptot,nlay), intent(in) :: gasvmr_co2
      real (kind=kind_phys), dimension(nf_vgas), intent(in) :: gasvmr_other
      integer, intent(in) :: nf_vgas
      integer :: k1, rr

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay, nbands), intent(in):: tauaer

      real (kind=kind_phys) :: rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(nxjp_acc_length, nlay, ng00) :: small_fracs, &
     &       small_tautot

!  ---  locals:
      integer :: k, ind0, ind1, inds, indsp, indf, indfp, indm, indmp,  &
     &       id000, id010, id100, id110, id200, id210, jmco2, jpl,      &
     &       id001, id011, id101, id111, id201, id211, jmco2p, jplp,    &
     &       jmco, jmcop, ig, js, js1, jpr

      real (kind=kind_phys) :: tauself, taufor, co2m1, co2m2, absco2,   &
     &       speccomb,       specparm,       specmult,       fs,        &
     &       speccomb1,      specparm1,      specmult1,      fs1,       &
     &       speccomb_mco2,  specparm_mco2,  specmult_mco2,  fmco2,     &
     &       speccomb_mco,   specparm_mco,   specmult_mco,   fmco,      &
     &       speccomb_planck,specparm_planck,specmult_planck,fpl,       &
     &       refrat_planck_a, refrat_m_a, refrat_m_a3, ratco2,          &
     &       adjfac, adjcolco2, com1, com2, absco, abso3,               &
     &       fac000, fac100, fac200, fac010, fac110, fac210,            &
     &       fac001, fac101, fac201, fac011, fac111, fac211,            &
     &       p0, p40, fk00, fk10, fk20, p1, p41, fk01, fk11, fk21, temp, taug
      integer :: jj, iplon, iplon2, ib, n
      ! GPU: variable name changed: CPU - colamt(:,:,1),       GPU - colamt1
      ! GPU: variable name changed: CPU - colamt(:,:,2),       GPU - colamt2
      ! GPU: variable name changed: CPU - colamt(:,:,3),       GPU - colamt3
      ! GPU: variable name changed: CPU - colamt(:,:,4),       GPU - colamt4
      ! GPU: variable name changed: CPU - colamt(:,:,7),       GPU - colamt7
!
!===> ...  begin here
!
!  --- ...  minor gas mapping levels :
!     lower - co2, p = 1053.63 mb, t = 294.2 k
!     lower - co, p = 706 mb, t = 278.94 k
!     upper - o3, p = 95.5835 mb, t = 215.7 k

!  --- ...  calculate reference ratio to be used in calculation of planck
     !           fraction in lower/upper atmosphere.
      refrat_planck_a = chi_mls(1,5)/chi_mls(4,5)        ! p = 473.420 mb (level 5)
      refrat_m_a = chi_mls(1,1)/chi_mls(4,1)             ! p = 1053. (level 1)
      refrat_m_a3 = chi_mls(1,3)/chi_mls(4,3)            ! p = 706. (level 3)
      !$acc parallel loop collapse(2) private(speccomb, specparm, specmult, js, &
      !$acc&         fs, ind0, speccomb1, specparm1, specmult1, js1, fs1, ind1, &
      !$acc&         speccomb_planck, specparm_planck, specmult_planck, jpl, &
      !$acc&         fpl, inds, indf, indm, indsp, indfp, indmp, jplp, &
      !$acc&         id000, id010, speccomb_mco2, specparm_mco2, specmult_mco2, &
      !$acc&         id100, id110, id200, id210, fac000, fac100, fac200, fac010, &
      !$acc&         fac110, fac210, id001, id011, id101, id111, id201, id211, &
      !$acc&         fac001, fac101, fac201, fac011, fac111, fac211, tauself, &
      !$acc&         taufor, p0, p40, fk00, fk10, fk20, p1, jmco2, fmco2, &
      !$acc&         p41, fk01, fk11, fk21, temp, adjfac, &
      !$acc&         speccomb_mco, specparm_mco, specmult_mco, jmco, fmco, &
      !$acc&         jmco2p, jmcop, ratco2, adjcolco2, co2m1, co2m2, absco2, &
      !$acc&         com1, com2, absco, taug, jpr, rfrate, &
      !$acc&         k1, temcol, colamt1, colamt2, colamt3, colamt4, colamt7, rr) async(async_id)
      do k = 1, nlay
         do n = 1, nxjp_acc_length
            jpr = jp(n, k)
            rr = n + jbs_nxjp_acc - 1
            if (ivflip == 0) then
               k1 = nlay + 1 - k
            else
               k1 = k
            end if
            temcol = 1.0e-12 * coldry(n, k)
            colamt1 = max(f_zero,    coldry(n, k)*h2ovmr(n, k))          ! h2o
            colamt2 = max(temcol, coldry(n, k)*gasvmr_co2(rr,k1)) ! co2
            colamt3 = max(temcol, coldry(n, k)*o3vmr(n, k))           ! o3
            if (ilwrgas > 0) then
               colamt4 = max(temcol, coldry(n, k)*gasvmr_other(2))  ! n2o
               colamt7 = max(f_zero,    coldry(n, k)*gasvmr_other(5))  ! co
            else
               colamt4 = f_zero     ! n2o
               colamt7 = f_zero     ! co
            endif
            if (k .le. laytrop(n)) then
               !  --- ...  lower atmosphere loop
               rfrate = chi_mls(1,jpr) / chi_mls(4,jpr)
               speccomb = colamt1 + rfrate*colamt4
               specparm = colamt1 / speccomb
               specmult = 8.0 * min(specparm, oneminus)
               js = 1 + int(specmult)
               fs = mod(specmult, f_one)
               ind0 = ((jp(n, k)-1)*5 + (jt(n, k)-1)) * nspa(13) + js

               rfrate = chi_mls(1,jpr+1) / chi_mls(4,jpr+1)
               speccomb1 = colamt1 + rfrate*colamt4
               specparm1 = colamt1 / speccomb1
               specmult1 = 8.0 * min(specparm1, oneminus)
               js1 = 1 + int(specmult1)
               fs1 = mod(specmult1, f_one)
               ind1 = (jp(n, k)*5 + (jt1(n, k)-1)) * nspa(13) + js1

               speccomb_mco2 = colamt1 + refrat_m_a*colamt4
               specparm_mco2 = colamt1 / speccomb_mco2
               specmult_mco2 = 8.0 * min(specparm_mco2, oneminus)
               jmco2 = 1 + int(specmult_mco2)
               fmco2 = mod(specmult_mco2, f_one)

               !  --- ...  in atmospheres where the amount of co2 is too great to be considered
               !           a minor species, adjust the column amount of co2 by an empirical factor
               !           to obtain the proper contribution.

               speccomb_mco = colamt1 + refrat_m_a3*colamt4
               specparm_mco = colamt1 / speccomb_mco
               specmult_mco = 8.0 * min(specparm_mco, oneminus)
               jmco = 1 + int(specmult_mco)
               fmco = mod(specmult_mco, f_one)

               speccomb_planck = colamt1 + refrat_planck_a*colamt4
               specparm_planck = colamt1 / speccomb_planck
               specmult_planck = 8.0 * min(specparm_planck, oneminus)
               jpl = 1 + int(specmult_planck)
               fpl = mod(specmult_planck, f_one)

               inds = indself(n, k)
               indf = indfor(n, k)
               indm = indminor(n, k)
               indsp = inds + 1
               indfp = indf + 1
               indmp = indm + 1
               jplp  = jpl  + 1
               jmco2p= jmco2+ 1
               jmcop = jmco + 1

               !  --- ...  in atmospheres where the amount of co2 is too great to be considered
               !           a minor species, adjust the column amount of co2 by an empirical factor
               !           to obtain the proper contribution.

               temp   = coldry(n, k) * 3.55e-4
               ratco2 = colamt2 / temp
               if (ratco2 > 3.0) then
                  adjfac = 2.0 + (ratco2-2.0)**0.68
                  adjcolco2 = adjfac * temp
               else
                  adjcolco2 = colamt2
               endif

               if (specparm < 0.125 .and. specparm1 < 0.125) then
                  p0 = fs - f_one
                  p40 = p0**4
                  fk00 = p40
                  fk10 = f_one - p0 - 2.0*p40
                  fk20 = p0 + p40

                  p1 = fs1 - f_one
                  p41 = p1**4
                  fk01 = p41
                  fk11 = f_one - p1 - 2.0*p41
                  fk21 = p1 + p41

                  id000 = ind0
                  id010 = ind0 + 9
                  id100 = ind0 + 1
                  id110 = ind0 +10
                  id200 = ind0 + 2
                  id210 = ind0 +11

                  id001 = ind1
                  id011 = ind1 + 9
                  id101 = ind1 + 1
                  id111 = ind1 +10
                  id201 = ind1 + 2
                  id211 = ind1 +11
               elseif (specparm > 0.875 .and. specparm1 > 0.875) then
                  p0 = -fs
                  p40 = p0**4
                  fk00 = p40
                  fk10 = f_one - p0 - 2.0*p40
                  fk20 = p0 + p40

                  p1 = -fs1
                  p41 = p1**4
                  fk01 = p41
                  fk11 = f_one - p1 - 2.0*p41
                  fk21 = p1 + p41

                  id000 = ind0 + 1
                  id010 = ind0 +10
                  id100 = ind0
                  id110 = ind0 + 9
                  id200 = ind0 - 1
                  id210 = ind0 + 8

                  id001 = ind1 + 1
                  id011 = ind1 +10
                  id101 = ind1
                  id111 = ind1 + 9
                  id201 = ind1 - 1
                  id211 = ind1 + 8
               else
                  fk00 = f_one - fs
                  fk10 = fs
                  fk20 = f_zero

                  fk01 = f_one - fs1
                  fk11 = fs1
                  fk21 = f_zero

                  id000 = ind0
                  id010 = ind0 + 9
                  id100 = ind0 + 1
                  id110 = ind0 +10
                  id200 = ind0
                  id210 = ind0

                  id001 = ind1
                  id011 = ind1 + 9
                  id101 = ind1 + 1
                  id111 = ind1 +10
                  id201 = ind1
                  id211 = ind1
               endif

               fac000 = fk00 * fac00(n, k)
               fac100 = fk10 * fac00(n, k)
               fac200 = fk20 * fac00(n, k)
               fac010 = fk00 * fac10(n, k)
               fac110 = fk10 * fac10(n, k)
               fac210 = fk20 * fac10(n, k)

               fac001 = fk01 * fac01(n, k)
               fac101 = fk11 * fac01(n, k)
               fac201 = fk21 * fac01(n, k)
               fac011 = fk01 * fac11(n, k)
               fac111 = fk11 * fac11(n, k)
               fac211 = fk21 * fac11(n, k)
               !$acc loop seq
               do ig = 1, ng13
                  ib = ngb(ns13+ig)
                  tauself = selffac(n, k)* (selfref(ig,inds) + selffrac(n, k)         &
                  &            * (selfref(ig,indsp) - selfref(ig,inds)))
                  taufor  = forfac(n, k) * (forref(ig,indf) + forfrac(n, k)           &
                  &            * (forref(ig,indfp) - forref(ig,indf))) 
                  co2m1   = ka_mco2(ig,jmco2,indm) + fmco2                      &
                  &            * (ka_mco2(ig,jmco2p,indm) - ka_mco2(ig,jmco2,indm))
                  co2m2   = ka_mco2(ig,jmco2,indmp) + fmco2                     &
                  &            * (ka_mco2(ig,jmco2p,indmp) - ka_mco2(ig,jmco2,indmp))
                  absco2  = co2m1 + minorfrac(n, k) * (co2m2 - co2m1)
                  com1    = ka_mco(ig,jmco,indm) + fmco                         &
                  &            * (ka_mco(ig,jmcop,indm) - ka_mco(ig,jmco,indm))
                  com2    = ka_mco(ig,jmco,indmp) + fmco                        &
                  &            * (ka_mco(ig,jmcop,indmp) - ka_mco(ig,jmco,indmp))
                  absco   = com1 + minorfrac(n, k) * (com2 - com1)

                  taug = speccomb                                    &
                  &                * (fac000*absa(ig,id000) + fac010*absa(ig,id010)  &
                  &                +  fac100*absa(ig,id100) + fac110*absa(ig,id110)  &
                  &                +  fac200*absa(ig,id200) + fac210*absa(ig,id210)) &
                  &                +     speccomb1                                   &
                  &                * (fac001*absa(ig,id001) + fac011*absa(ig,id011)  &
                  &                +  fac101*absa(ig,id101) + fac111*absa(ig,id111)  &
                  &                +  fac201*absa(ig,id201) + fac211*absa(ig,id211)) &
                  &                + tauself + taufor + adjcolco2*absco2             &
                  &                + colamt7*absco

                  !fracs(iplon, k, ns13+ig, jj) = fracrefa(ig,jpl) + fpl                     &
                  !&                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
                  !tautot(iplon, k, ns13+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) = fracrefa(ig,jpl) + fpl                     &
                  &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)
               enddo
            else
               indm = indminor(n, k)
               indmp = indm + 1
               !$acc loop seq
               do ig = 1, ng13
                  ib = ngb(ns13+ig)
                  abso3 = kb_mo3(ig,indm) + minorfrac(n, k)                        &
                  &          * (kb_mo3(ig,indmp) - kb_mo3(ig,indm))

                  taug = colamt3*abso3

                  !fracs(iplon, k, ns13+ig, jj) =  fracrefb(ig)
                  !tautot(iplon, k, ns13+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) =  fracrefb(ig)
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)
               enddo
            end if
         end do
      end do

! ..................................
      end subroutine taugb13
! ----------------------------------

! ----------------------------------
      subroutine taugb14 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!     band 14:  2250-2380 cm-1 (low - co2; high - co2)                 !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb14
      integer, intent(in) :: nlay, laytrop(nxjp_acc_length), ix, myim(fulljj), &
         async_id, ng00, i2, blockjj, fulljj, nxjp_acc_length, jjoffset, &
         max_nxjp_acc_length, jbs_nxjp_acc
      integer, dimension(nxjp_acc_length), intent(in) :: map_jj, map_i

      integer, dimension(max_nxjp_acc_length, nlay), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys) :: temcol, colamt1, colamt2, colamt3, &
                               colamt4, colamt5, colamt6, colamt7
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: h2ovmr, o3vmr
      real (kind=kind_phys), dimension(nxptot,nlay), intent(in) :: gasvmr_co2
      real (kind=kind_phys), dimension(nf_vgas), intent(in) :: gasvmr_other
      integer, intent(in) :: nf_vgas
      integer :: k1, rr

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay, nbands), intent(in):: tauaer

      real (kind=kind_phys) :: rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(nxjp_acc_length, nlay, ng00) :: small_fracs, &
     &       small_tautot

!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       ig

      real (kind=kind_phys) :: tauself, taufor, taug
      integer :: jj, iplon, iplon2, ib, n
      ! GPU: variable name changed: CPU - colamt(:,:,1),       GPU - colamt1
      ! GPU: variable name changed: CPU - colamt(:,:,2),       GPU - colamt2
!
!===> ...  begin here
!
      !$acc parallel loop collapse(2) private(ind0, ind1, inds, indf, ind0p, rr, &
      !$acc&         ind1p, indsp, indfp, tauself, taufor, taug, k1, temcol, colamt2) async(async_id)
      do k = 1, nlay
         do n = 1, nxjp_acc_length
            rr = n + jbs_nxjp_acc - 1
            if (ivflip == 0) then
               k1 = nlay + 1 - k
            else
               k1 = k
            end if
            temcol = 1.0e-12 * coldry(n, k)
            colamt2 = max(temcol, coldry(n, k)*gasvmr_co2(rr,k1)) ! co2
            if (k .le. laytrop(n)) then
               !  --- ...  lower atmosphere loop
               ind0 = ((jp(n, k)-1)*5 + (jt (n, k)-1)) * nspa(14) + 1
               ind1 = ( jp(n, k)   *5 + (jt1(n, k)-1)) * nspa(14) + 1

               inds = indself(n, k)
               indf = indfor(n, k)
               ind0p = ind0 + 1
               ind1p = ind1 + 1
               indsp = inds + 1
               indfp = indf + 1
               !$acc loop seq
               do ig = 1, ng14
                  ib = ngb(ns14+ig)
                  tauself = selffac(n, k) * (selfref(ig,inds) + selffrac(n, k)        &
                  &            * (selfref(ig,indsp) - selfref(ig,inds)))
                  taufor  = forfac(n, k) * (forref(ig,indf) + forfrac(n, k)           &
                  &            * (forref(ig,indfp) - forref(ig,indf))) 

                  taug = colamt2                                 &
                  &            * (fac00(n, k)*absa(ig,ind0) + fac10(n, k)*absa(ig,ind0p)   &
                  &            +  fac01(n, k)*absa(ig,ind1) + fac11(n, k)*absa(ig,ind1p))  &
                  &            + tauself + taufor

                  !fracs(iplon, k, ns14+ig, jj) = fracrefa(ig)
                  !tautot(iplon, k, ns14+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) = fracrefa(ig)
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)

               enddo
            else
               ind0 = ((jp(n, k)-13)*5 + (jt (n, k)-1)) * nspb(14) + 1
               ind1 = ((jp(n, k)-12)*5 + (jt1(n, k)-1)) * nspb(14) + 1

               ind0p = ind0 + 1
               ind1p = ind1 + 1
               !$acc loop seq
               do ig = 1, ng14
                  ib = ngb(ns14+ig)
                  taug = colamt2                                 &
                  &             * (fac00(n, k)*absb(ig,ind0) + fac10(n, k)*absb(ig,ind0p)  &
                  &             +  fac01(n, k)*absb(ig,ind1) + fac11(n, k)*absb(ig,ind1p))

                  !fracs(iplon, k, ns14+ig, jj) = fracrefb(ig)
                  !tautot(iplon, k, ns14+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) = fracrefb(ig)
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)
               enddo
            end if
         end do
      end do

! ..................................
      end subroutine taugb14
! ----------------------------------

! ----------------------------------
      subroutine taugb15 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!     band 15:  2380-2600 cm-1 (low - n2o,co2; low minor - n2)         !
!                              (high - nothing)                        !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb15
      integer, intent(in) :: nlay, laytrop(nxjp_acc_length), ix, myim(fulljj), &
         async_id, ng00, i2, blockjj, fulljj, nxjp_acc_length, jjoffset, &
         max_nxjp_acc_length, jbs_nxjp_acc
      integer, dimension(nxjp_acc_length), intent(in) :: map_jj, map_i

      integer, dimension(max_nxjp_acc_length, nlay), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys) :: temcol, colamt1, colamt2, colamt3, &
                               colamt4, colamt5, colamt6, colamt7
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: h2ovmr, o3vmr
      real (kind=kind_phys), dimension(nxptot,nlay), intent(in) :: gasvmr_co2
      real (kind=kind_phys), dimension(nf_vgas), intent(in) :: gasvmr_other
      integer, intent(in) :: nf_vgas
      integer :: k1, rr

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay, nbands), intent(in):: tauaer

      real (kind=kind_phys) :: rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(nxjp_acc_length, nlay, ng00) :: small_fracs, &
     &       small_tautot

!  ---  locals:
      integer :: k, ind0, ind1, inds, indsp, indf, indfp, indm, indmp,  &
     &       id000, id010, id100, id110, id200, id210, jpl, jplp,       &
     &       id001, id011, id101, id111, id201, id211, jmn2, jmn2p,     &
     &       ig, js, js1, jpr

      real (kind=kind_phys) :: scalen2, tauself, taufor,                &
     &       speccomb,       specparm,       specmult,       fs,        &
     &       speccomb1,      specparm1,      specmult1,      fs1,       &
     &       speccomb_mn2,   specparm_mn2,   specmult_mn2,   fmn2,      &
     &       speccomb_planck,specparm_planck,specmult_planck,fpl,       &
     &       refrat_planck_a, refrat_m_a, n2m1, n2m2, taun2,            &
     &       fac000, fac100, fac200, fac010, fac110, fac210,            &
     &       fac001, fac101, fac201, fac011, fac111, fac211,            &
     &       p0, p40, fk00, fk10, fk20, p1, p41, fk01, fk11, fk21, taug
      integer :: jj, iplon, iplon2, ib, n
      ! GPU: variable name changed: CPU - colamt(:,:,2),       GPU - colamt2
      ! GPU: variable name changed: CPU - colamt(:,:,4),       GPU - colamt4
!
!===> ...  begin here
!
!  --- ...  minor gas mapping level :
!     lower - nitrogen continuum, p = 1053., t = 294.

!  --- ...  calculate reference ratio to be used in calculation of planck
!           fraction in lower atmosphere.
      refrat_planck_a = chi_mls(4,1)/chi_mls(2,1)      ! p = 1053. mb (level 1)
      refrat_m_a = chi_mls(4,1)/chi_mls(2,1)           ! p = 1053. mb
      !$acc parallel loop collapse(2) private(speccomb, specparm, specmult, js, &
      !$acc&         fs, ind0, speccomb1, specparm1, specmult1, js1, fs1, ind1, &
      !$acc&         speccomb_planck, specparm_planck, specmult_planck, jpl, &
      !$acc&         fpl, inds, indf, indm, indsp, indfp, indmp, jplp, &
      !$acc&         id000, id010, speccomb_mn2, specparm_mn2, specmult_mn2, &
      !$acc&         id100, id110, id200, id210, fac000, fac100, fac200, fac010, &
      !$acc&         fac110, fac210, id001, id011, id101, id111, id201, id211, &
      !$acc&         fac001, fac101, fac201, fac011, fac111, fac211, tauself, &
      !$acc&         taufor, p0, p40, fk00, fk10, fk20, p1, jmn2, fmn2, &
      !$acc&         p41, fk01, fk11, fk21, scalen2, jmn2p, n2m1, n2m2, taun2, taug, &
      !$acc&         jpr, rfrate, k1, temcol, colamt2, colamt4, rr) async(async_id)
      do k = 1, nlay
         do n = 1, nxjp_acc_length
            jpr = jp(n, k)
            rr = n + jbs_nxjp_acc - 1
            if (ivflip == 0) then
               k1 = nlay + 1 - k
            else
               k1 = k
            end if
            temcol = 1.0e-12 * coldry(n, k)
            colamt1 = max(f_zero,    coldry(n, k)*h2ovmr(n, k))          ! h2o
            colamt2 = max(temcol, coldry(n, k)*gasvmr_co2(rr,k1)) ! co2
            if (ilwrgas > 0) then
               colamt4 = max(temcol, coldry(n, k)*gasvmr_other(2))  ! n2o
            else
               colamt4 = f_zero     ! n2o
            endif
            if (k .le. laytrop(n)) then
               !  --- ...  lower atmosphere loop
               rfrate = chi_mls(4,jpr) / chi_mls(2,jpr)
               speccomb = colamt4 + rfrate*colamt2
               specparm = colamt4 / speccomb
               specmult = 8.0 * min(specparm, oneminus)
               js = 1 + int(specmult)
               fs = mod(specmult, f_one)
               ind0 = ((jp(n, k)-1)*5 + (jt(n, k)-1)) * nspa(15) + js

               rfrate = chi_mls(4,jpr+1) / chi_mls(2,jpr+1)
               speccomb1 = colamt4 + rfrate*colamt2
               specparm1 = colamt4 / speccomb1
               specmult1 = 8.0 * min(specparm1, oneminus)
               js1 = 1 + int(specmult1)
               fs1 = mod(specmult1, f_one)
               ind1 = (jp(n, k)*5 + (jt1(n, k)-1)) * nspa(15) + js1

               speccomb_mn2 = colamt4 + refrat_m_a*colamt2
               specparm_mn2 = colamt4 / speccomb_mn2
               specmult_mn2 = 8.0 * min(specparm_mn2, oneminus)
               jmn2 = 1 + int(specmult_mn2)
               fmn2 = mod(specmult_mn2, f_one)

               speccomb_planck = colamt4 + refrat_planck_a*colamt2
               specparm_planck = colamt4 / speccomb_planck
               specmult_planck = 8.0 * min(specparm_planck, oneminus)
               jpl = 1 + int(specmult_planck)
               fpl = mod(specmult_planck, f_one)

               scalen2 = colbrd(n, k) * scaleminor(n, k)

               inds = indself(n, k)
               indf = indfor(n, k)
               indm = indminor(n, k)
               indsp = inds + 1
               indfp = indf + 1
               indmp = indm + 1
               jplp  = jpl  + 1
               jmn2p = jmn2 + 1

               if (specparm < 0.125 .and. specparm1 < 0.125) then
                  p0 = fs - f_one
                  p40 = p0**4
                  fk00 = p40
                  fk10 = f_one - p0 - 2.0*p40
                  fk20 = p0 + p40

                  p1 = fs1 - f_one
                  p41 = p1**4
                  fk01 = p41
                  fk11 = f_one - p1 - 2.0*p41
                  fk21 = p1 + p41

                  id000 = ind0
                  id010 = ind0 + 9
                  id100 = ind0 + 1
                  id110 = ind0 +10
                  id200 = ind0 + 2
                  id210 = ind0 +11

                  id001 = ind1
                  id011 = ind1 + 9
                  id101 = ind1 + 1
                  id111 = ind1 +10
                  id201 = ind1 + 2
                  id211 = ind1 +11
               elseif (specparm > 0.875 .and. specparm1 > 0.875) then
                  p0 = -fs
                  p40 = p0**4
                  fk00 = p40
                  fk10 = f_one - p0 - 2.0*p40
                  fk20 = p0 + p40

                  p1 = -fs1
                  p41 = p1**4
                  fk01 = p41
                  fk11 = f_one - p1 - 2.0*p41
                  fk21 = p1 + p41

                  id000 = ind0 + 1
                  id010 = ind0 +10
                  id100 = ind0
                  id110 = ind0 + 9
                  id200 = ind0 - 1
                  id210 = ind0 + 8

                  id001 = ind1 + 1
                  id011 = ind1 +10
                  id101 = ind1
                  id111 = ind1 + 9
                  id201 = ind1 - 1
                  id211 = ind1 + 8
               else
                  fk00 = f_one - fs
                  fk10 = fs
                  fk20 = f_zero

                  fk01 = f_one - fs1
                  fk11 = fs1
                  fk21 = f_zero

                  id000 = ind0
                  id010 = ind0 + 9
                  id100 = ind0 + 1
                  id110 = ind0 +10
                  id200 = ind0
                  id210 = ind0

                  id001 = ind1
                  id011 = ind1 + 9
                  id101 = ind1 + 1
                  id111 = ind1 +10
                  id201 = ind1
                  id211 = ind1
               endif

               fac000 = fk00 * fac00(n, k)
               fac100 = fk10 * fac00(n, k)
               fac200 = fk20 * fac00(n, k)
               fac010 = fk00 * fac10(n, k)
               fac110 = fk10 * fac10(n, k)
               fac210 = fk20 * fac10(n, k)

               fac001 = fk01 * fac01(n, k)
               fac101 = fk11 * fac01(n, k)
               fac201 = fk21 * fac01(n, k)
               fac011 = fk01 * fac11(n, k)
               fac111 = fk11 * fac11(n, k)
               fac211 = fk21 * fac11(n, k)
               !$acc loop seq
               do ig = 1, ng15
                  ib = ngb(ns15+ig)
                  tauself = selffac(n, k)* (selfref(ig,inds) + selffrac(n, k)         &
                  &            * (selfref(ig,indsp) - selfref(ig,inds)))
                  taufor  = forfac(n, k) * (forref(ig,indf) + forfrac(n, k)           &
                  &            * (forref(ig,indfp) - forref(ig,indf))) 
                  n2m1    = ka_mn2(ig,jmn2,indm) + fmn2                         &
                  &            * (ka_mn2(ig,jmn2p,indm) - ka_mn2(ig,jmn2,indm))
                  n2m2    = ka_mn2(ig,jmn2,indmp) + fmn2                        &
                  &            * (ka_mn2(ig,jmn2p,indmp) - ka_mn2(ig,jmn2,indmp))
                  taun2   = scalen2 * (n2m1 + minorfrac(n, k) * (n2m2 - n2m1))

                  taug = speccomb                                    &
                  &                * (fac000*absa(ig,id000) + fac010*absa(ig,id010)  &
                  &                +  fac100*absa(ig,id100) + fac110*absa(ig,id110)  &
                  &                +  fac200*absa(ig,id200) + fac210*absa(ig,id210)) &
                  &                +     speccomb1                                   &
                  &                * (fac001*absa(ig,id001) + fac011*absa(ig,id011)  &
                  &                +  fac101*absa(ig,id101) + fac111*absa(ig,id111)  &
                  &                +  fac201*absa(ig,id201) + fac211*absa(ig,id211)) &
                  &                + tauself + taufor + taun2

                  !fracs(iplon, k, ns15+ig, jj) = fracrefa(ig,jpl) + fpl                     &
                  !&                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
                  !tautot(iplon, k, ns15+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) = fracrefa(ig,jpl) + fpl                     &
                  &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)

               enddo
            else
               !$acc loop seq
               do ig = 1, ng15
                  ib = ngb(ns15+ig)
                  taug = f_zero

                  !fracs(iplon, k, ns15+ig, jj) = f_zero
                  !tautot(iplon, k, ns15+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) = f_zero
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)

               enddo
            end if
         end do
      end do

! ..................................
      end subroutine taugb15
! ----------------------------------

! ----------------------------------
      subroutine taugb16 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,nf_vgas, h2ovmr, gasvmr_co2, &
             o3vmr, gasvmr_other, colbrd,tauaer,              &
     &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, ng00, i2, map_jj, map_i, nxjp_acc_length, jjoffset, &
             max_nxjp_acc_length, jbs_nxjp_acc, async_id, fulljj, blockjj,                    &
!  ---  outputs:
     &       small_fracs, small_tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!     band 16:  2600-3250 cm-1 (low key- h2o,ch4; high key - ch4)      !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb16
      integer, intent(in) :: nlay, laytrop(nxjp_acc_length), ix, myim(fulljj), &
         async_id, ng00, i2, blockjj, fulljj, nxjp_acc_length, jjoffset, &
         max_nxjp_acc_length, jbs_nxjp_acc
      integer, dimension(nxjp_acc_length), intent(in) :: map_jj, map_i

      integer, dimension(max_nxjp_acc_length, nlay), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys) :: temcol, colamt1, colamt2, colamt3, &
                               colamt4, colamt5, colamt6, colamt7
      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay), intent(in) :: h2ovmr, o3vmr
      real (kind=kind_phys), dimension(nxptot,nlay), intent(in) :: gasvmr_co2
      real (kind=kind_phys), dimension(nf_vgas), intent(in) :: gasvmr_other
      integer, intent(in) :: nf_vgas
      integer :: k1, rr

      real (kind=kind_phys), dimension(max_nxjp_acc_length, nlay, nbands), intent(in):: tauaer

      real (kind=kind_phys) :: rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(nxjp_acc_length, nlay, ng00) :: small_fracs, &
     &       small_tautot
!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       id000, id010, id100, id110, id200, id210, jpl, jplp,       &
     &       id001, id011, id101, id111, id201, id211, ig, js, js1, jpr

      real (kind=kind_phys) :: tauself, taufor, refrat_planck_a,        &
     &       speccomb,       specparm,       specmult,       fs,        &
     &       speccomb1,      specparm1,      specmult1,      fs1,       &
     &       speccomb_planck,specparm_planck,specmult_planck,fpl,       &
     &       fac000, fac100, fac200, fac010, fac110, fac210,            &
     &       fac001, fac101, fac201, fac011, fac111, fac211,            &
     &       p0, p40, fk00, fk10, fk20, p1, p41, fk01, fk11, fk21, taug
      integer :: jj, iplon, iplon2, ib, n
      ! GPU: variable name changed: CPU - colamt(:,:,1),       GPU - colamt1
      ! GPU: variable name changed: CPU - colamt(:,:,5),       GPU - colamt5
!
!===> ...  begin here
!
!  --- ...  calculate reference ratio to be used in calculation of planck
!           fraction in lower atmosphere.
      refrat_planck_a = chi_mls(1,6)/chi_mls(6,6)        ! p = 387. mb (level 6)
      !$acc parallel loop collapse(2) private(speccomb, specparm, specmult, js, &
      !$acc&         fs, ind0, speccomb1, specparm1, specmult1, js1, fs1, ind1, &
      !$acc&         speccomb_planck, specparm_planck, specmult_planck, jpl, &
      !$acc&         fpl, inds, indf, indsp, indfp, jplp, id000, id010, &
      !$acc&         id100, id110, id200, id210, fac000, fac100, fac200, fac010, &
      !$acc&         fac110, fac210, id001, id011, id101, id111, id201, id211, &
      !$acc&         fac001, fac101, fac201, fac011, fac111, fac211, tauself, &
      !$acc&         taufor, p0, p40, fk00, fk10, fk20, p1, &
      !$acc&         p41, fk01, fk11, fk21, taug, jpr, rfrate, k1, temcol, &
      !$acc&         colamt1, colamt5, rr)  async(async_id)
      do k = 1, nlay
         do n = 1, nxjp_acc_length
            jpr = jp(n, k)
            rr = n + jbs_nxjp_acc - 1
            temcol = 1.0e-12 * coldry(n, k)
            colamt1 = max(f_zero,    coldry(n, k)*h2ovmr(n, k))          ! h2o
            if (ilwrgas > 0) then
               colamt5 = max(temcol, coldry(n, k)*gasvmr_other(3))  ! ch4
            else
               colamt5 = f_zero     ! ch4
            endif
            if (k .le. laytrop(n)) then
               !  --- ...  lower atmosphere loop
               rfrate = chi_mls(1,jpr) / chi_mls(6,jpr)
               speccomb = colamt1 + rfrate*colamt5
               specparm = colamt1 / speccomb
               specmult = 8.0 * min(specparm, oneminus)
               js = 1 + int(specmult)
               fs = mod(specmult, f_one)
               ind0 = ((jp(n, k)-1)*5 + (jt(n, k)-1)) * nspa(16) + js

               rfrate = chi_mls(1,jpr+1) / chi_mls(6,jpr+1)
               speccomb1 = colamt1 + rfrate*colamt5
               specparm1 = colamt1 / speccomb1
               specmult1 = 8.0 * min(specparm1, oneminus)
               js1 = 1 + int(specmult1)
               fs1 = mod(specmult1, f_one)
               ind1 = (jp(n, k)*5 + (jt1(n, k)-1)) * nspa(16) + js1

               speccomb_planck = colamt1 + refrat_planck_a*colamt5
               specparm_planck = colamt1 / speccomb_planck
               specmult_planck = 8.0 * min(specparm_planck, oneminus)
               jpl = 1 + int(specmult_planck)
               fpl = mod(specmult_planck, f_one)

               inds = indself(n, k)
               indf = indfor(n, k)
               indsp = inds + 1
               indfp = indf + 1
               jplp  = jpl  + 1

               if (specparm < 0.125 .and. specparm1 < 0.125) then
                  p0 = fs - f_one
                  p40 = p0**4
                  fk00 = p40
                  fk10 = f_one - p0 - 2.0*p40
                  fk20 = p0 + p40

                  p1 = fs1 - f_one
                  p41 = p1**4
                  fk01 = p41
                  fk11 = f_one - p1 - 2.0*p41
                  fk21 = p1 + p41

                  id000 = ind0
                  id010 = ind0 + 9
                  id100 = ind0 + 1
                  id110 = ind0 +10
                  id200 = ind0 + 2
                  id210 = ind0 +11

                  id001 = ind1
                  id011 = ind1 + 9
                  id101 = ind1 + 1
                  id111 = ind1 +10
                  id201 = ind1 + 2
                  id211 = ind1 +11
               elseif (specparm > 0.875 .and. specparm1 > 0.875) then
                  p0 = -fs
                  p40 = p0**4
                  fk00 = p40
                  fk10 = f_one - p0 - 2.0*p40
                  fk20 = p0 + p40

                  p1 = -fs1
                  p41 = p1**4
                  fk01 = p41
                  fk11 = f_one - p1 - 2.0*p41
                  fk21 = p1 + p41

                  id000 = ind0 + 1
                  id010 = ind0 +10
                  id100 = ind0
                  id110 = ind0 + 9
                  id200 = ind0 - 1
                  id210 = ind0 + 8

                  id001 = ind1 + 1
                  id011 = ind1 +10
                  id101 = ind1
                  id111 = ind1 + 9
                  id201 = ind1 - 1
                  id211 = ind1 + 8
               else
                  fk00 = f_one - fs
                  fk10 = fs
                  fk20 = f_zero

                  fk01 = f_one - fs1
                  fk11 = fs1
                  fk21 = f_zero

                  id000 = ind0
                  id010 = ind0 + 9
                  id100 = ind0 + 1
                  id110 = ind0 +10
                  id200 = ind0
                  id210 = ind0

                  id001 = ind1
                  id011 = ind1 + 9
                  id101 = ind1 + 1
                  id111 = ind1 +10
                  id201 = ind1
                  id211 = ind1
               endif

               fac000 = fk00 * fac00(n, k)
               fac100 = fk10 * fac00(n, k)
               fac200 = fk20 * fac00(n, k)
               fac010 = fk00 * fac10(n, k)
               fac110 = fk10 * fac10(n, k)
               fac210 = fk20 * fac10(n, k)

               fac001 = fk01 * fac01(n, k)
               fac101 = fk11 * fac01(n, k)
               fac201 = fk21 * fac01(n, k)
               fac011 = fk01 * fac11(n, k)
               fac111 = fk11 * fac11(n, k)
               fac211 = fk21 * fac11(n, k)
               !$acc loop seq
               do ig = 1, ng16
                  ib = ngb(ns16+ig)
                  tauself = selffac(n, k)* (selfref(ig,inds) + selffrac(n, k)         &
                  &            * (selfref(ig,indsp) - selfref(ig,inds)))
                  taufor  = forfac(n, k) * (forref(ig,indf) + forfrac(n, k)           &
                  &            * (forref(ig,indfp) - forref(ig,indf))) 

                  taug = speccomb                                    &
                  &                * (fac000*absa(ig,id000) + fac010*absa(ig,id010)  &
                  &                +  fac100*absa(ig,id100) + fac110*absa(ig,id110)  &
                  &                +  fac200*absa(ig,id200) + fac210*absa(ig,id210)) &
                  &                +     speccomb1                                   &
                  &                * (fac001*absa(ig,id001) + fac011*absa(ig,id011)  &
                  &                +  fac101*absa(ig,id101) + fac111*absa(ig,id111)  &
                  &                +  fac201*absa(ig,id201) + fac211*absa(ig,id211)) &
                  &                + tauself + taufor

                  !fracs(iplon, k, ns16+ig, jj) = fracrefa(ig,jpl) + fpl                     &
                  !&                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
                  !tautot(iplon, k, ns16+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) = fracrefa(ig,jpl) + fpl                     &
                  &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)

               enddo
            else
               ind0 = ((jp(n, k)-13)*5 + (jt (n, k)-1)) * nspb(16) + 1
               ind1 = ((jp(n, k)-12)*5 + (jt1(n, k)-1)) * nspb(16) + 1

               ind0p = ind0 + 1
               ind1p = ind1 + 1
               !$acc loop seq
               do ig = 1, ng16
                  ib = ngb(ns16+ig)
                  taug = colamt5                                 &
                  &           * (fac00(n, k)*absb(ig,ind0) + fac10(n, k)*absb(ig,ind0p)    &
                  &           +  fac01(n, k)*absb(ig,ind1) + fac11(n, k)*absb(ig,ind1p))

                  !fracs(iplon, k, ns16+ig, jj) = fracrefb(ig)
                  !tautot(iplon, k, ns16+ig, jj) = taug + tauaer(n, k, ib)
                  small_fracs(n, k, ig) = fracrefb(ig)
                  small_tautot(n, k, ig) = taug + tauaer(n, k, ib)

               enddo
            end if
         end do
      end do

! ..................................
      end subroutine taugb16
! ----------------------------------



!
!........................................!
      end module module_radlw_main_gpu       !
!========================================!

