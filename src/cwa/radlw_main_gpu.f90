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

      use module_radlw_parameters
!
      use module_radlw_avplank, only : totplnk
      use module_radlw_ref,     only : preflog, tref, chi_mls
      use param,             only : my
      use index,             only : jlistnum
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

      public lwrad_gpu, rlwinit_gpu, taumol_interface_gpu, taumol, rtrnmr, &
         setcoef, cldprop


! ================
      contains
! ================
      subroutine taumol_interface_gpu                                                 &
! ..................................
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, async_id, fulljj,                                                         &
!  ---  outputs:
     &       fracs, tautot                                              &
     &     )
     implicit none
!  ---  inputs:
      integer, intent(in) :: nlay, laytrop(ix, fulljj), ix, myim(fulljj), async_id, fulljj

      integer, dimension(ix, nlay, fulljj), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(ix, nlay, fulljj), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2

      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj), intent(in):: colamt
      real (kind=kind_phys), dimension(ix, nlay,maxxsec, fulljj),intent(in):: wx

      real (kind=kind_phys), dimension(ix, nlay, nbands, fulljj), intent(in):: tauaer

      real (kind=kind_phys), dimension(ix, nlay,nrates,2, fulljj), intent(in) ::    &
     &       rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(ix, nlay, ngptlw, fulljj), intent(out) ::     &
     &       fracs, tautot

     call taumol( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, async_id, fulljj,                                                         &
!  ---  outputs:
     &       fracs, tautot                                              &
     &     )
     

      end subroutine
! --------------------------------
      subroutine lwrad_gpu                                                  &
! --------------------------------

!  ---  inputs:
     &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                      &
     &       clouds,icseed,aerosols,sfemis,sfgtmp,                      &
     &       myim, nlay, nlp1, lprnt, myrank, ix,                       &
             nf_vgas, nf_clds, nf_aelw, async_id, fulljj, blocks, smalljj, &
!  ---  outputs:
     &       hlwc,topflx,sfcflx                                         &
!! ---  optional:
     &,      hlw0,hlwb,flxprf                                           &
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

!  ---  inputs:
      integer, intent(in) :: myim(fulljj), nlay, nlp1, myrank, fulljj, &
         nf_vgas, nf_clds, nf_aelw, ix, blocks, smalljj
      integer, intent(in) :: icseed(ix, fulljj)

      logical,  intent(in) :: lprnt

      real (kind=kind_phys), dimension(ix, nlp1, fulljj), intent(in) :: plvl,  &
     &       tlvl
      real (kind=kind_phys), dimension(ix, nlay, fulljj), intent(in) :: plyr,  &
     &       tlyr, qlyr, olyr

      real (kind=kind_phys), dimension(ix,nlay,nf_vgas, fulljj),intent(in):: gasvmr
      real (kind=kind_phys), dimension(ix,nlay,nf_clds, fulljj),intent(in):: clouds

      real (kind=kind_phys), dimension(ix, fulljj), intent(in) :: sfemis,     &
     &       sfgtmp

      real (kind=kind_phys), dimension(ix,nlay,nbands,nf_aelw, fulljj),intent(in):: &
     &       aerosols

!  ---  outputs:
      real (kind=kind_phys), dimension(ix, nlay, fulljj), intent(out) :: hlwc

      type (topflw_type),    dimension(ix, fulljj), intent(out) :: topflx
      type (sfcflw_type),    dimension(ix, fulljj), intent(out) :: sfcflx

!! ---  optional outputs:
      real (kind=kind_phys), dimension(ix,nlay,nbands, fulljj),optional,      &
     &       intent(out) :: hlwb
      real (kind=kind_phys), dimension(ix, nlay, fulljj),       optional,      &
     &       intent(out) :: hlw0
      type (proflw_type),    dimension(ix, nlp1, fulljj),       optional,      &
     &       intent(out) :: flxprf

!  ---  locals:
      real (kind=kind_phys), dimension(ix, 0:nlp1, smalljj) :: cldfrc

      real (kind=kind_phys), dimension(ix, 0:nlay, smalljj) :: totuflux, totdflux,   &
     &       totuclfl, totdclfl, tz

      real (kind=kind_phys), dimension(ix, nlay, smalljj)   :: htr, htrcl

      real (kind=kind_phys), dimension(ix, nlay, smalljj)   :: pavel, tavel, delp,   &
     &       clwp, ciwp, relw, reiw, cda1, cda2, cda3, cda4,            &
     &       coldry, colbrd, h2ovmr, o3vmr, fac00, fac01, fac10, fac11, &
     &       selffac, selffrac, forfac, forfrac, minorfrac, scaleminor, &
     &       scaleminorn2, temcol

      real (kind=kind_phys), dimension(ix,0:nlay, nbands, smalljj) :: pklev, pklay

      real (kind=kind_phys), dimension(ix, nlay,nbands, smalljj) :: htrb
      real (kind=kind_phys), dimension(ix, nlay,nbands, smalljj) :: taucld, tauaer
      real (kind=kind_phys), dimension(ix, nlay,ngptlw, smalljj) :: fracs, tautot

      real (kind=kind_phys), dimension(ix, nbands, smalljj) :: semiss, secdiff

!  ---  column amount of absorbing gases:
!       (:,m) m = 1-h2o, 2-co2, 3-o3, 4-n2o, 5-ch4, 6-o2, 7-co
      real (kind=kind_phys) :: colamt(ix, nlay,maxgas, smalljj)

!  ---  column cfc cross-section amounts:
!       (:,m) m = 1-ccl4, 2-cfc11, 3-cfc12, 4-cfc22
      real (kind=kind_phys) :: wx(ix, nlay,maxxsec, smalljj)

!  ---  reference ratios of binary species parameter in lower atmosphere:
!       (:,m,:) m = 1-h2o/co2, 2-h2o/o3, 3-h2o/n2o, 4-h2o/ch4, 5-n2o/co2, 6-o3/co2
      real (kind=kind_phys) :: rfrate(ix, nlay,nrates,2, smalljj)

      real (kind=kind_phys) :: tem0, tem1, tem2, pwvcm(ix, smalljj), summol, stemp

      integer, dimension(ix, smalljj) :: ipseed
      integer, dimension(ix, nlay, smalljj) :: jp, jt, jt1, indself, indfor, indminor
      integer                  :: laytrop(ix, smalljj), iplon, i, j, k, k1, jj
      logical :: lcf1(ix, smalljj)
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
      integer :: async_id, ig, ib
      real (kind=kind_phys), allocatable, dimension(:,:,:)   :: cldf
      logical,allocatable,dimension(:,:,:,:) :: lcloudy
      real (kind=kind_phys), allocatable, dimension(:)   :: cldf_im
      logical,allocatable,dimension(:,:) :: lcloudy_im

      real (kind=kind_phys), dimension(:,:,:,:), allocatable :: cldfmc
      integer :: jb, jf, jjoffset, jbs, jbe
      
      if (isubclw > 0) allocate(cldfmc(ix, nlay,ngptlw, fulljj))


      !
      !===> ... begin here
      !

      !  --- ...  initialization

      lhlwb  = present ( hlwb )
      lhlw0  = present ( hlw0 )
      lflxprf= present ( flxprf )

      !      if (myrank .eq. 0) print *,'$$$$ lwrad start $$$$$'
      !      if (myrank .eq. 0) print *,'$$$$ lhlwb   = ', lhlwb
      !      if (myrank .eq. 0) print *,'$$$$ lhlw0   = ', lhlw0
      !      if (myrank .eq. 0) print *,'$$$$ lflxprf = ', lflxprf
      !$acc data create(cldfrc, totuflux, totdflux, totuclfl, totdclfl, tz, &
      !$acc&     htr, htrcl, pavel, tavel, delp, clwp, ciwp, relw, reiw, cda1, &
      !$acc&     cda2, cda3, cda4, coldry, colbrd, h2ovmr, o3vmr, fac00, fac01, &
      !$acc&     fac10, fac11, selffac, selffrac, forfac, forfrac, minorfrac, &
      !$acc&     scaleminor, scaleminorn2, temcol, pklev, pklay, htrb, taucld, tauaer, &
      !$acc&     fracs, tautot, semiss, secdiff, colamt, wx, rfrate, &
      !$acc&     pwvcm, ipseed, jp, jt, jt1, indself, indfor, indminor, laytrop, lcf1) &
      !$acc&     async(async_id)

      do jb = 1, blocks
         jjoffset = (jb-1)*smalljj
         jbs = jjoffset+1
         jbe = jjoffset+smalljj
      !$acc parallel loop gang collapse(3) async(async_id)
      do jj = 1, smalljj
         do j = 1, maxgas
            do k = 1, nlay
               !$acc loop vector private(jf)
               do i = 1, ix
                  jf = jjoffset+jj
                  colamt(i, k, j, jj) = f_zero
               end do
            end do
         end do 
      end do

      !, jj  --- ...  change random number seed value for each radiation invocation

      if     ( isubclw == 1 ) then     ! advance prescribed permutation seed
         !$acc parallel loop collapse(2) private(jf) async(async_id)
         do jj = 1, smalljj
            do i = 1, ix
               jf = jjoffset+jj
               if (i .le. myim(jf)) then
                  ipseed(i, jj) = ipsdlw0 + i
               end if
            enddo
         end do
      elseif ( isubclw == 2 ) then     ! use input array of permutaion seeds
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
      !       print *,'  in radlw, isubclw, ipsdlw0,ipseed =',                &
      !    &          isubclw, ipsdlw0, ipseed
      !     endif

      !  --- ...  loop over horizontal npts profiles
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, smalljj
         do j = 1, nbands
            jf = jjoffset+jj
            !$acc loop vector private(jf)
            do iplon = 1, myim(jf) ! lab_do_iplon
               if (sfemis(iplon, jf) > eps .and. sfemis(iplon, jf) <= 1.0) then  ! input surface emissivity
                  semiss(iplon, j, jj) = sfemis(iplon, jf)
               else                                                      ! use default values
                  semiss(iplon, j, jj) = semiss0(j)
               endif
            enddo
         end do
      end do

            !stemp = sfgtmp(iplon, jj)          ! surface ground temp

            !  --- ...  prepare atmospheric profile for use in rrtm
            !           the vertical index of internal array is from surface to top

            !  --- ...  molecular amounts are input or converted to volume mixing ratio
            !           and later then converted to molecular amount (molec/cm2) by the
            !           dry air column coldry (in molec/cm2) which is calculated from the
            !           layer pressure thickness (in mb), based on the hydrostatic equation
            !  --- ...  and includes a correction to account for h2o in the layer.
      if (ivflip == 0) then
         tem1 = 100.0 * con_g
         tem2 = 1.0e-20 * 1.0e3 * con_avgd
         !$acc parallel loop collapse(2) private(jf) async(async_id)
         do jj = 1, smalljj
            do iplon = 1, ix
               jf = jjoffset+jj
               if (iplon .le. myim(jf)) then ! lab_do_iplon
                  tz(iplon, 0, jj) = tlvl(iplon,nlp1, jf)
               end if
            end do
         end do

         !$acc parallel loop gang collapse(2) private(jf) async(async_id)
         do jj = 1, smalljj
            do k = 1, nlay
               jf = jjoffset+jj
               !$acc loop vector private(tem0)
               do iplon = 1, myim(jf) ! lab_do_iplon
                  k1 = nlp1 - k
                  pavel(iplon, k, jj)= plyr(iplon,k1, jf)
                  delp(iplon, k, jj) = plvl(iplon,k1+1, jf) - plvl(iplon,k1, jf)
                  tavel(iplon, k, jj)= tlyr(iplon,k1, jf)
                  tz(iplon, k, jj)   = tlvl(iplon,k1, jf)

                  !  --- ...  set absorber amount
                  !test use
                  !           h2ovmr(k)= max(f_zero,qlyr(iplon,k1)*amdw)                  ! input mass mixing ratio
                  !           h2ovmr(k)= max(f_zero,qlyr(iplon,k1))                       ! input vol mixing ratio
                  !           o3vmr (k)= max(f_zero,olyr(iplon,k1))                       ! input vol mixing ratio
                  !ncep model use
                  h2ovmr(iplon, k, jj)= max(f_zero,qlyr(iplon,k1, jf)                        &
                  &                           *amdw/(f_one-qlyr(iplon,k1, jf)))          ! input specific humidity
                  o3vmr (iplon, k, jj)= max(f_zero,olyr(iplon,k1, jf)*amdo3)                 ! input mass mixing ratio

                  !  --- ...  tem0 is the molecular weight of moist air
                  tem0 = (f_one - h2ovmr(iplon, k, jj))*con_amd + h2ovmr(iplon, k, jj)*con_amw
                  coldry(iplon, k, jj) = tem2*delp(iplon, k, jj) / (tem1*tem0*(f_one+h2ovmr(iplon, k, jj)))
                  temcol(iplon, k, jj) = 1.0e-12 * coldry(iplon, k, jj)

                  colamt(iplon, k,1, jj) = max(f_zero,    coldry(iplon, k, jj)*h2ovmr(iplon, k, jj))          ! h2o
                  colamt(iplon, k,2, jj) = max(temcol(iplon, k, jj), coldry(iplon, k, jj)*gasvmr(iplon,k1,1, jf)) ! co2
                  colamt(iplon, k,3, jj) = max(temcol(iplon, k, jj), coldry(iplon, k, jj)*o3vmr(iplon, k, jj))           ! o3

               !  --- ...  set up col amount for rare gases, convert from volume mixing ratio
               !           to molec/cm2 based on coldry (scaled to 1.0e-20)
                  if (ilwrgas > 0) then
                     k1 = nlp1 - k
                     colamt(iplon, k,4, jj)=max(temcol(iplon, k, jj), coldry(iplon, k, jj)*gasvmr(iplon,k1,2, jf))  ! n2o
                     colamt(iplon, k,5, jj)=max(temcol(iplon, k, jj), coldry(iplon, k, jj)*gasvmr(iplon,k1,3, jf))  ! ch4
                     colamt(iplon, k,6, jj)=max(f_zero,    coldry(iplon, k, jj)*gasvmr(iplon,k1,4, jf))  ! o2
                     colamt(iplon, k,7, jj)=max(f_zero,    coldry(iplon, k, jj)*gasvmr(iplon,k1,5, jf))  ! co

                     wx(iplon, k,1, jj) = max( f_zero, coldry(iplon, k, jj)*gasvmr(iplon,k1,9, jf) )   ! ccl4
                     wx(iplon, k,2, jj) = max( f_zero, coldry(iplon, k, jj)*gasvmr(iplon,k1,6, jf) )   ! cf11
                     wx(iplon, k,3, jj) = max( f_zero, coldry(iplon, k, jj)*gasvmr(iplon,k1,7, jf) )   ! cf12
                     wx(iplon, k,4, jj) = max( f_zero, coldry(iplon, k, jj)*gasvmr(iplon,k1,8, jf) )   ! cf22
                  else
                     colamt(iplon, k,4, jj) = f_zero     ! n2o
                     colamt(iplon, k,5, jj) = f_zero     ! ch4
                     colamt(iplon, k,6, jj) = f_zero     ! o2
                     colamt(iplon, k,7, jj) = f_zero     ! co

                     wx(iplon, k,1, jj) = f_zero
                     wx(iplon, k,2, jj) = f_zero
                     wx(iplon, k,3, jj) = f_zero
                     wx(iplon, k,4, jj) = f_zero
                  endif
               enddo
            end do
         end do

               !  --- ...  set aerosol optical properties
         !$acc parallel loop gang collapse(3) private(jf) async(async_id)
         do jj = 1, smalljj
            do j = 1, nbands
               do k = 1, nlay
                  jf = jjoffset+jj
                  !$acc loop vector private(k1)
                  do iplon = 1, myim(jf) ! lab_do_iplon
                     k1 = nlp1 - k
                     tauaer(iplon, k,j, jj) = aerosols(iplon,k1,j,1, jf)                      &
                     &                    * (f_one - aerosols(iplon,k1,j,2, jf))
                  enddo
               enddo
            end do
         end do
         if (ilwcliq > 0) then    ! use prognostic cloud method
            !$acc parallel loop gang collapse(2) private(jf) async(async_id)
            do jj = 1, smalljj
               do k = 1, nlay
                  jf = jjoffset+jj
                  !$acc loop vector private(k1)
                  do iplon = 1, myim(jf) ! lab_do_iplon
                     k1 = nlp1 - k
                     cldfrc(iplon, k, jj)= clouds(iplon,k1,1, jf)
                     clwp(iplon, k, jj)  = clouds(iplon,k1,2, jf)
                     relw(iplon, k, jj)  = clouds(iplon,k1,3, jf)
                     ciwp(iplon, k, jj)  = clouds(iplon,k1,4, jf)
                     reiw(iplon, k, jj)  = clouds(iplon,k1,5, jf)
                     cda1(iplon, k, jj)  = clouds(iplon,k1,6, jf)
                     cda2(iplon, k, jj)  = clouds(iplon,k1,7, jf)
                     cda3(iplon, k, jj)  = clouds(iplon,k1,8, jf)
                     cda4(iplon, k, jj)  = clouds(iplon,k1,9, jf)
                  enddo
               end do
            end do
         else                       ! use diagnostic cloud method
            !$acc parallel loop gang collapse(2) private(jf) async(async_id)
            do jj = 1, smalljj
               do k = 1, nlay
                  jf = jjoffset+jj
                  !$acc loop vector private(k1)
                  do iplon = 1, myim(jf) ! lab_do_iplon
                     k1 = nlp1 - k
                     cldfrc(iplon, k, jj)= clouds(iplon,k1,1, jf)
                     cda1(iplon, k, jj)  = clouds(iplon,k1,2, jf)
                  enddo
               end do
            end do
         endif                      ! end if_ilwcliq
         !$acc parallel loop collapse(2) private(jf) async(async_id)
         do jj = 1, smalljj
            do iplon = 1, ix
               jf = jjoffset+jj
               if (iplon .le. myim(jf)) then ! lab_do_iplon
                  cldfrc(iplon, 0, jj)    = f_one       ! padding value only
                  cldfrc(iplon, nlp1, jj) = f_zero      ! padding value only
               end if
            end do
         end do

         !  --- ...  compute precipitable water vapor for diffusivity angle adjustments
         !$acc parallel loop collapse(2) private(tem1, tem2, tem0, jf) async(async_id)
         do jj = 1, smalljj
            do iplon = 1, ix
               jf = jjoffset+jj
               if (iplon .le. myim(jf)) then ! lab_do_iplon
                  tem1 = f_zero
                  tem2 = f_zero
                  !$acc loop seq
                  do k = 1, nlay
                     tem1 = tem1 + coldry(iplon, k, jj) + colamt(iplon, k,1, jj)
                     tem2 = tem2 + colamt(iplon, k,1, jj)
                  enddo

                  tem0 = 10.0 * tem2 / (amdw * tem1 * con_g)
                  pwvcm(iplon, jj) = tem0 * plvl(iplon,nlp1, jf)
               end if
            end do
         end do
      end if

      if (ivflip .ne. 0) then                        ! input from sfc to toa
         tem1 = 100.0 * con_g
         tem2 = 1.0e-20 * 1.0e3 * con_avgd
         !$acc parallel loop collapse(2) private(jf) async(async_id)
         do jj = 1, smalljj
            do iplon = 1, ix
               jf = jjoffset+jj
               if (iplon .le. myim(jf)) then ! lab_do_iplon
                  tz(iplon, 0, jj) = tlvl(iplon,1, jf)
               end if
            end do
         end do
         !$acc parallel loop gang collapse(2) private(jf) async(async_id)
         do jj = 1, smalljj
            do k = 1, nlay
               jf = jjoffset+jj
               !$acc loop vector private(tem0)
               do iplon = 1, myim(jf) ! lab_do_iplon
                  pavel(iplon, k, jj)= plyr(iplon,k, jf)
                  delp(iplon, k, jj) = plvl(iplon,k, jf) - plvl(iplon,k+1, jf)
                  tavel(iplon, k, jj)= tlyr(iplon,k, jf)
                  tz(iplon, k, jj)   = tlvl(iplon,k+1, jf)

                  !  --- ...  set absorber amount
                  !test use
                  !           h2ovmr(k)= max(f_zero,qlyr(iplon,k)*amdw)                   ! input mass mixing ratio
                  !           h2ovmr(k)= max(f_zero,qlyr(iplon,k))                        ! input vol mixing ratio
                  !           o3vmr (k)= max(f_zero,olyr(iplon,k))                        ! input vol mixing ratio
                  !ncep model use
                  h2ovmr(iplon, k, jj)= max(f_zero,qlyr(iplon,k, jf)                         &
                  &                           *amdw/(f_one-qlyr(iplon,k, jf)))           ! input specific humidity
                  o3vmr (iplon, k, jj)= max(f_zero,olyr(iplon,k, jf)*amdo3)                  ! input mass mixing ratio

                  !  --- ...  tem0 is the molecular weight of moist air
                  tem0 = (f_one - h2ovmr(iplon, k, jj))*con_amd + h2ovmr(iplon, k, jj)*con_amw
                  coldry(iplon, k, jj) = tem2*delp(iplon, k, jj) / (tem1*tem0*(f_one+h2ovmr(iplon, k, jj)))
                  temcol(iplon, k, jj) = 1.0e-12 * coldry(iplon, k, jj)

                  colamt(iplon, k,1, jj) = max(f_zero,    coldry(iplon, k, jj)*h2ovmr(iplon, k, jj))          ! h2o
                  colamt(iplon, k,2, jj) = max(temcol(iplon, k, jj), coldry(iplon, k, jj)*gasvmr(iplon,k,1, jf))  ! co2
                  colamt(iplon, k,3, jj) = max(temcol(iplon, k, jj), coldry(iplon, k, jj)*o3vmr(iplon, k, jj))           ! o3

               !  --- ...  set up col amount for rare gases, convert from volume mixing ratio
               !           to molec/cm2 based on coldry (scaled to 1.0e-20)
                  if (ilwrgas > 0) then
                     colamt(iplon, k,4, jj)=max(temcol(iplon, k, jj), coldry(iplon, k, jj)*gasvmr(iplon,k,2, jf))  ! n2o
                     colamt(iplon, k,5, jj)=max(temcol(iplon, k, jj), coldry(iplon, k, jj)*gasvmr(iplon,k,3, jf))  ! ch4
                     colamt(iplon, k,6, jj)=max(f_zero,    coldry(iplon, k, jj)*gasvmr(iplon,k,4, jf))  ! o2
                     colamt(iplon, k,7, jj)=max(f_zero,    coldry(iplon, k, jj)*gasvmr(iplon,k,5, jf))  ! co

                     wx(iplon, k,1, jj) = max( f_zero, coldry(iplon, k, jj)*gasvmr(iplon,k,9, jf) )   ! ccl4
                     wx(iplon, k,2, jj) = max( f_zero, coldry(iplon, k, jj)*gasvmr(iplon,k,6, jf) )   ! cf11
                     wx(iplon, k,3, jj) = max( f_zero, coldry(iplon, k, jj)*gasvmr(iplon,k,7, jf) )   ! cf12
                     wx(iplon, k,4, jj) = max( f_zero, coldry(iplon, k, jj)*gasvmr(iplon,k,8, jf) )   ! cf22
                  else
                     colamt(iplon, k,4, jj) = f_zero     ! n2o
                     colamt(iplon, k,5, jj) = f_zero     ! ch4
                     colamt(iplon, k,6, jj) = f_zero     ! o2
                     colamt(iplon, k,7, jj) = f_zero     ! co

                     wx(iplon, k,1, jj) = f_zero
                     wx(iplon, k,2, jj) = f_zero
                     wx(iplon, k,3, jj) = f_zero
                     wx(iplon, k,4, jj) = f_zero
                  endif
               end do
            end do
         enddo



               !  --- ...  set aerosol optical properties
         !$acc parallel loop gang collapse(3) private(jf) async(async_id)
         do jj = 1, smalljj
            do j = 1, nbands
               do k = 1, nlay
                  jf = jjoffset+jj
                  !$acc loop vector
                  do iplon = 1, myim(jf) ! lab_do_iplon
                     tauaer(iplon, k,j, jj) = aerosols(iplon,k,j,1, jf)                       &
                     &                    * (f_one - aerosols(iplon,k,j,2, jf))
                  enddo
               enddo
            end do
         end do
         if (ilwcliq > 0) then    ! use prognostic cloud method
            !$acc parallel loop gang collapse(2) private(jf) async(async_id)
            do jj = 1, smalljj
               do k = 1, nlay
                  jf = jjoffset+jj
                  !$acc loop vector
                  do iplon = 1, myim(jf) ! lab_do_iplon
                     cldfrc(iplon, k, jj)= clouds(iplon,k,1, jf)
                     clwp(iplon, k, jj)  = clouds(iplon,k,2, jf)
                     relw(iplon, k, jj)  = clouds(iplon,k,3, jf)
                     ciwp(iplon, k, jj)  = clouds(iplon,k,4, jf)
                     reiw(iplon, k, jj)  = clouds(iplon,k,5, jf)
                     cda1(iplon, k, jj)  = clouds(iplon,k,6, jf)
                     cda2(iplon, k, jj)  = clouds(iplon,k,7, jf)
                     cda3(iplon, k, jj)  = clouds(iplon,k,8, jf)
                     cda4(iplon, k, jj)  = clouds(iplon,k,9, jf)
                  enddo
               end do
            end do
         else                       ! use diagnostic cloud method
            !$acc parallel loop gang collapse(2) private(jf) async(async_id)
            do jj = 1, smalljj
               do k = 1, nlay
                  jf = jjoffset+jj
                  !$acc loop vector
                  do iplon = 1, myim(jf) ! lab_do_iplon
                     cldfrc(iplon, k, jj)= clouds(iplon,k,1, jf)
                     cda1(iplon, k, jj)  = clouds(iplon,k,2, jf)
                  enddo
               end do
            end do
         endif                      ! end if_ilwcliq
      

         !  --- ...  compute precipitable water vapor for diffusivity angle adjustments
         !$acc parallel loop collapse(2) private(tem1, tem2, tem0, jf) async(async_id)
         do jj = 1, smalljj
            do iplon = 1, ix
               jf = jjoffset+jj
               if (iplon .le. myim(jf)) then ! lab_do_iplon
                  cldfrc(iplon, 0, jj)    = f_one       ! padding value only
                  cldfrc(iplon, nlp1, jj) = f_zero      ! padding value only
                  tem1 = f_zero
                  tem2 = f_zero
                  !$acc loop seq
                  do k = 1, nlay
                     tem1 = tem1 + coldry(iplon, k, jj) + colamt(iplon, k,1, jj)
                     tem2 = tem2 + colamt(iplon, k,1, jj)
                  enddo

                  tem0 = 10.0 * tem2 / (amdw * tem1 * con_g)
                  pwvcm(iplon, jj) = tem0 * plvl(iplon,1, jf)
               end if
            end do
         end do
      endif                       ! if_ivflip

            !  --- ...  compute column amount for broadening gases
      !$acc parallel loop gang collapse(2) private(jf) async(async_id)
      do jj = 1, smalljj
         do k = 1, nlay
            jf = jjoffset+jj
            !$acc loop vector private(summol)
            do iplon = 1, myim(jf)
               summol = f_zero
               !$acc loop seq
               do i = 2, maxgas
                  summol = summol + colamt(iplon, k,i, jj)
               enddo
               colbrd(iplon, k, jj) = coldry(iplon, k, jj) - summol
            enddo
         end do
      end do

      !  --- ...  compute diffusivity angle adjustments

      tem1 = 1.80
      tem2 = 1.50
      !$acc parallel loop gang collapse(2) private(jf) async(async_id)
      do jj = 1, smalljj
         do j = 1, nbands
            jf = jjoffset+jj
            !$acc loop vector
            do iplon = 1, myim(jf) ! lab_do_iplon
               if (j==1 .or. j==4 .or. j==10) then
                  secdiff(iplon, j, jj) = 1.66
               else
                  secdiff(iplon, j, jj) = min( tem1, max( tem2,                          &
                  &                   a0(j)+a1(j)*exp(a2(j)*pwvcm(iplon, jj)) ))
               endif
            enddo
         end do
      end do

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
      !$acc parallel loop collapse(2) private(jf) async(async_id)
      do jj = 1, smalljj
         do iplon = 1, ix
            jf = jjoffset+jj
            if (iplon .le. myim(jf)) then ! lab_do_iplon
               lcf1(iplon, jj) = .false.
               !$acc loop seq
               do k = 1, nlay ! lab_do_k0
                  if ( cldfrc(iplon, k, jj) > eps ) then
                     lcf1(iplon, jj) = .true.
                     exit ! lab_do_k0
                  endif
               enddo  ! lab_do_k0
            end if
         end do
      end do
      !call nvtxStartRange("lw_cldprop")
      call cldprop                                                  &
         !  ---  inputs:
         &     ( cldfrc,clwp,relw,ciwp,reiw,cda1,cda2,cda3,cda4,            &
         &       nlay, nlp1, ipseed, ix, myim(jbs:jbe), lcf1, async_id, smalljj,                                  &
         !  ---  outputs:
         &       taucld                                             &
         &     )
      if ( isubclw > 0 ) then      ! mcica sub-col clouds approx
         allocate(lcloudy(ix, ngptlw,nlay, fulljj))
         allocate(cldf_im(nlay))
         allocate(lcloudy_im(ngptlw,nlay))
         allocate(cldf(ix, nlay, fulljj))
         !$acc parallel loop gang collapse(3) private(jf) async(async_id)
         do jj = 1, smalljj
            do ig = 1, ngptlw
               do k = 1, nlay
                  jf = jjoffset+jj
                  !$acc loop vector
                  do iplon = 1, myim(jf) ! lab_do_iplon
                     cldfmc(iplon, k,ig, jj) = f_zero
                  enddo
               enddo
            end do
         end do
         do jj = 1, smalljj
            jf = jjoffset+jj
            do iplon = 1, myim(jf) ! lab_do_iplon
               if ( lcf1(iplon, jj) ) then
               !  ---  distribute cloud properties to each g-point

                  do k = 1, nlay
                     if ( cldfrc(iplon, k, jj) < cldmin ) then
                        cldf(iplon, k, jj) = f_zero
                     else
                        cldf(iplon, k, jj) = cldfrc(iplon, k, jj)
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
                     &     ( cldf_im, nlay, ipseed(iplon, jj),                                        &
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
      call setcoef                                                    &
         !  ---  inputs:
         &     ( pavel,tavel,tz,sfgtmp(:,jbs:jbe),h2ovmr,colamt,coldry,colbrd,          &
         &       nlay, nlp1, ix, myim(jbs:jbe), async_id, smalljj,                                               &
         !  ---  outputs:
         &       laytrop,pklay,pklev,jp,jt,jt1,                             &
         &       rfrate,fac00,fac01,fac10,fac11,                            &
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
      !      write(myrank_str,'(I3)') myrank
      !      open(unit=1000, file='lw_taumol_input.'//trim(adjustl(myrank_str)), form='unformatted', &
      !         access='stream', status='replace')
      !      write(1000) laytrop,pavel,coldry,colbrd,             &
      !   &       fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
      !   &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
      !   &       minorfrac,scaleminor,scaleminorn2,indminor,                &
      !   &       nlay, ix, myim, wx,tauaer, colamt, rfrate,          &
      !   !  ---  outputs:
      !   &       fracs, tautot  
      !      close(1000)
      !      write(*,*) nlay, ix, myim
      !call nvtxStartRange("taumol")
      call taumol                                                     &
         !  ---  inputs:
         &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
         &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
         &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
         &       minorfrac,scaleminor,scaleminorn2,indminor,                &
         &       nlay, ix, myim(jbs:jbe), async_id, smalljj,                                                   &
         !  ---  outputs:
         &       fracs, tautot                                              &
         &     )
      !call nvtxEndRange

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
      if (isubclw <= 0) then
         if (iovrlw <= 0) then
            do jj = 1, smalljj
               jf = jjoffset+jj
               do iplon = 1, myim(jf) ! lab_do_iplon
                  do k = 1, nlay
                     delp_im(k) = delp(iplon, k, jj)
                     htr_im(k) = htr(iplon, k, jj)
                     htrcl_im(k) = htrcl(iplon, k, jj)
                  end do
                  do k = 0, nlay
                     totuflux_im(k) = totuflux(iplon, k, jj)
                     totdflux_im(k) = totdflux(iplon, k, jj)
                     totuclfl_im(k) = totuclfl(iplon, k, jj)
                     totdclfl_im(k) = totdclfl(iplon, k, jj)
                  end do
                  do k = 0, nlp1
                     cldfrc_im(k) = cldfrc(iplon, k, jj)
                  end do
                  do j = 1, nbands
                     semiss_im(j) = semiss(iplon, j, jj)
                     secdiff_im(j) = secdiff(iplon, j, jj)
                  end do
                  call rtrn                                                   &
                     !  ---  inputs:
                     &     ( semiss,delp_im,cldfrc_im,taucld_im,tautot_im,pklay_im,pklev_im,              &
                     &       fracs_im,secdiff,nlay,nlp1,                                   &
                     !  ---  outputs:
                     &       totuflux_im,totdflux_im,htr_im, totuclfl_im,totdclfl_im,htrcl_im, htrb_im       &
                     &     )
                  do k = 1, nlay
                     delp(iplon, k, jj) = delp_im(k)
                     htr(iplon, k, jj) = htr_im(k)
                     htrcl(iplon, k, jj) = htrcl_im(k)
                  end do
                  do k = 0, nlay
                     totuflux(iplon, k, jj) = totuflux_im(k)
                     totdflux(iplon, k, jj) = totdflux_im(k)
                     totuclfl(iplon, k, jj) = totuclfl_im(k)
                     totdclfl(iplon, k, jj) = totdclfl_im(k)
                  end do
                  do k = 0, nlp1
                     cldfrc(iplon, k, jj) = cldfrc_im(k)
                  end do
                  do j = 1, nbands
                     semiss(iplon, j, jj) = semiss_im(j)
                     secdiff(iplon, j, jj) = secdiff_im(j)
                  end do
               end do
            end do

         else
            !call nvtxStartRange("rtrnmr")
            !write(myrank_str,'(I3)') myrank
            !open(unit=1000, file='lw_rtrnmr_input.'//trim(adjustl(myrank_str)), form='unformatted', &
            !   access='stream', status='replace')
            !write(1000) delp,cldfrc,              &
            !   &       nlay,nlp1, ix, myim, &
            !   &       totuflux,totdflux,htr, totuclfl,totdclfl,htrcl, &
            !   semiss, secdiff, taucld, pklev, pklay, htrb, fracs, tautot
            !close(1000)
            !stop
            call rtrnmr                                                 &
               !  ---  inputs:
               &     ( semiss,delp,cldfrc,taucld,tautot,pklay,pklev,              &
               &       fracs,secdiff,nlay,nlp1, ix, myim(jbs:jbe), async_id, smalljj,                                  &
               !  ---  outputs:
               &       totuflux,totdflux,htr, totuclfl,totdclfl,htrcl, htrb       &
               &     )
         endif   ! end if_iovrlw_block
            !call nvtxEndRange

      else
         do jj = 1, smalljj
            jf = jjoffset+jj
            do iplon = 1, myim(jf) ! lab_do_iplon
               do k = 1, nlay
                  delp_im(k) = delp(iplon, k, jj)
                  htr_im(k) = htr(iplon, k, jj)
                  htrcl_im(k) = htrcl(iplon, k, jj)
                  do j = 1, nbands
                     taucld_im(k, j) = taucld(iplon, k, j, jj)
                     htrb_im(k, j) = htrb(iplon, k, j, jj)
                  end do
                  do j = 1, ngptlw
                     tautot_im(k, j) = tautot(iplon, k, j, jj)
                     fracs_im(k, j) = fracs(iplon, k, j, jj)
                     cldfmc_im(k, j) = cldfmc(iplon, k, j, jj)
                  end do
               end do
               do k = 0, nlay
                  totuflux_im(k) = totuflux(iplon, k, jj)
                  totdflux_im(k) = totdflux(iplon, k, jj)
                  totuclfl_im(k) = totuclfl(iplon, k, jj)
                  totdclfl_im(k) = totdclfl(iplon, k, jj)
                  do j = 1, nbands
                     pklay(iplon, k, j, jj) = pklay_im(k, j)
                     pklev(iplon, k, j, jj) = pklev_im(k, j)
                  end do
               end do
               do j = 1, nbands
                  semiss_im(j) = semiss(iplon, j, jj)
                  secdiff_im(j) = secdiff(iplon, j, jj)
               end do
               call rtrnmc                                                   &
                  !  ---  inputs:
                  &     ( semiss_im,delp_im,cldfmc_im,taucld_im,tautot_im,pklay_im,pklev_im,              &
                  &       fracs_im,secdiff_im,nlay,nlp1,                                   &
                  !  ---  outputs:
                  &       totuflux_im,totdflux_im,htr_im, totuclfl_im,totdclfl_im,htrcl_im, htrb_im       &
                  &     )
               do k = 1, nlay
                  delp(iplon, k, jj) = delp_im(k)
                  htr(iplon, k, jj) = htr_im(k)
                  htrcl(iplon, k, jj) = htrcl_im(k)
                  do j = 1, nbands
                     taucld(iplon, k, j, jj) = taucld_im(k, j)
                     htrb(iplon, k, j, jj) = htrb_im(k, j)
                  end do
                  do j = 1, ngptlw
                     tautot(iplon, k, j, jj) = tautot_im(k, j)
                     fracs(iplon, k, j, jj) = fracs_im(k, j)
                     cldfmc(iplon, k, j, jj) = cldfmc_im(k, j)
                  end do
               end do
               do k = 0, nlay
                  totuflux(iplon, k, jj) = totuflux_im(k)
                  totdflux(iplon, k, jj) = totdflux_im(k)
                  totuclfl(iplon, k, jj) = totuclfl_im(k)
                  totdclfl(iplon, k, jj) = totdclfl_im(k)
                  do j = 1, nbands
                     pklay(iplon, k, j, jj) = pklay_im(k, j)
                     pklev(iplon, k, j, jj) = pklev_im(k, j)
                  end do
               end do
               do j = 1, nbands
                  semiss(iplon, j, jj) = semiss_im(j)
                  secdiff(iplon, j, jj) = secdiff_im(j)
               end do
            end do
         end do
      endif   ! end if_isubclw_block

            !  --- ...  output total-sky and clear-sky fluxes and heating rates
      !$acc parallel loop collapse(2) private(jf) async(async_id)
      do jj = 1, smalljj
         do iplon = 1, ix
            jf = jjoffset+jj
            if (iplon .le. myim(jf)) then ! lab_do_iplon
               topflx(iplon, jf)%upfxc = totuflux(iplon, nlay, jj)
               topflx(iplon, jf)%upfx0 = totuclfl(iplon, nlay, jj)

               sfcflx(iplon, jf)%upfxc = totuflux(iplon, 0, jj)
               sfcflx(iplon, jf)%upfx0 = totuclfl(iplon, 0, jj)
               sfcflx(iplon, jf)%dnfxc = totdflux(iplon, 0, jj)
               sfcflx(iplon, jf)%dnfx0 = totdclfl(iplon, 0, jj)
            end if
         end do
      end do

      if (ivflip == 0) then       ! output from toa to sfc
               !! --- ...  optional fluxes
         if ( lflxprf ) then
            !$acc parallel loop gang collapse(2) private(jf) async(async_id)
            do jj = 1, smalljj
               do k = 0, nlay
                  jf = jjoffset+jj
                  !$acc loop vector private(k1)
                  do iplon = 1, myim(jf) ! lab_do_iplon
                     k1 = nlp1 - k
                     flxprf(iplon,k1, jf)%upfxc = totuflux(iplon, k, jj)
                     flxprf(iplon,k1, jf)%dnfxc = totdflux(iplon, k, jj)
                     flxprf(iplon,k1, jf)%upfx0 = totuclfl(iplon, k, jj)
                     flxprf(iplon,k1, jf)%dnfx0 = totdclfl(iplon, k, jj)
                  enddo
               end do
            end do
         endif
         !$acc parallel loop gang collapse(2) private(jf) async(async_id)
         do jj = 1, smalljj
            do k = 1, nlay
               jf = jjoffset+jj
               !$acc loop vector private(k1)
               do iplon = 1, myim(jf) ! lab_do_iplon
                  k1 = nlp1 - k
                  hlwc(iplon,k1, jf) = htr(iplon, k, jj)
               enddo
            end do
         end do

               !! --- ...  optional clear sky heating rate
         if ( lhlw0 ) then
            !$acc parallel loop gang collapse(2) private(jf) async(async_id)
            do jj = 1, smalljj
               do k = 1, nlay
                  jf = jjoffset+jj
                  !$acc loop vector private(k1)
                  do iplon = 1, myim(jf) ! lab_do_iplon
                     k1 = nlp1 - k
                     hlw0(iplon,k1, jf) = htrcl(iplon, k, jj)
                  enddo
               end do
            end do
               
         endif

               !! --- ...  optional spectral band heating rate
         if ( lhlwb ) then
            !$acc parallel loop gang collapse(3) private(jf) async(async_id)
            do jj = 1, smalljj
               do j = 1, nbands
                  do k = 1, nlay
                     jf = jjoffset+jj
                     !$acc loop vector private(k1)
                     do iplon = 1, myim(jf) ! lab_do_iplon
                        k1 = nlp1 - k
                        hlwb(iplon,k1,j, jf) = htrb(iplon, k,j, jj)
                     enddo
                  enddo
               end do
            end do
         endif

      else                        ! output from sfc to toa
            !! --- ...  optional fluxes
         if ( lflxprf ) then
            !$acc parallel loop gang collapse(2) private(jf) async(async_id)
            do jj = 1, smalljj
               do k = 0, nlay
                  jf = jjoffset+jj
                  !$acc loop vector
                  do iplon = 1, myim(jf) ! lab_do_iplon
                     flxprf(iplon,k+1, jf)%upfxc = totuflux(iplon, k, jj)
                     flxprf(iplon,k+1, jf)%dnfxc = totdflux(iplon, k, jj)
                     flxprf(iplon,k+1, jf)%upfx0 = totuclfl(iplon, k, jj)
                     flxprf(iplon,k+1, jf)%dnfx0 = totdclfl(iplon, k, jj)
                  enddo
               end do
            end do
         endif
         !$acc parallel loop gang collapse(2) private(jf) async(async_id)
         do jj = 1, smalljj
            do k = 1, nlay
               jf = jjoffset+jj
               !$acc loop vector
               do iplon = 1, myim(jf) ! lab_do_iplon
                  hlwc(iplon,k, jf) = htr(iplon, k, jj)
               enddo
            end do
         end do

               !! --- ...  optional clear sky heating rate
         if ( lhlw0 ) then
            !$acc parallel loop gang collapse(2) private(jf) async(async_id)
            do jj = 1, smalljj
               do k = 1, nlay
                  jf = jjoffset+jj
                  !$acc loop vector
                  do iplon = 1, myim(jf) ! lab_do_iplon
                     hlw0(iplon,k, jf) = htrcl(iplon, k, jj)
                  enddo
               end do
            end do
         endif

            !! --- ...  optional spectral band heating rate
         if ( lhlwb ) then
            !$acc parallel loop gang collapse(3) private(jf) async(async_id)
            do jj = 1, smalljj
               do j = 1, nbands
                  do k = 1, nlay
                     jf = jjoffset+jj
                     !$acc loop vector
                     do iplon = 1, myim(jf) ! lab_do_iplon
                        hlwb(iplon,k,j, jf) = htrb(iplon, k,j, jj)
                     enddo
                  enddo
               end do
            end do
         endif

      endif                       ! if_ivflip
      end do
      !$acc end data
      if (isubclw > 0) deallocate(cldfmc)
      !$acc wait(async_id)

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
     &       nlay, nlp1, ipseed, ix, myim, lcf1, async_id, fulljj,                                          &
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
      integer, intent(in) :: nlay, nlp1, ipseed(ix, fulljj), ix, myim(fulljj), fulljj

      real (kind=kind_phys), dimension(ix, 0:nlp1, fulljj), intent(in) :: cfrac
      real (kind=kind_phys), dimension(ix, nlay, fulljj),   intent(in) :: cliqp,    &
     &       reliq, cicep, reice, cdat1, cdat2, cdat3, cdat4
     logical, dimension(ix, fulljj), intent(in) :: lcf1

!  ---  outputs:
      real (kind=kind_phys), dimension(ix, nlay, nbands, fulljj),intent(out):: taucld

!  ---  locals:
      real (kind=kind_phys), dimension(nbands) :: tauliq, tauice
      real (kind=kind_phys), allocatable, dimension(:,:,:)   :: cldf

      real (kind=kind_phys) :: dgeice, factor, fint, tauran, tausnw,    &
     &       cldliq, refliq, cldice, refice

      logical,allocatable,dimension(:,:,:,:) :: lcloudy
      integer :: ia, ib, ig, k, index, iplon, jj

      real (kind=kind_phys), allocatable, dimension(:)   :: cldf_im
      logical,allocatable,dimension(:,:) :: lcloudy_im

      integer :: async_id
!
!===> ...  begin here
!
      !$acc parallel loop gang collapse(3) async(async_id)
      do jj = 1, fulljj
         do ib = 1, nbands
            do k = 1, nlay
               do iplon = 1, myim(jj) ! lab_do_iplon
                  taucld(iplon, k,ib, jj) = f_zero
               enddo
            enddo
         end do
      end do

               !  --- ...  compute cloud radiative properties for a cloudy column
      if (ilwcliq > 0) then ! lab_if_ilwcliq


         !  --- ...  calculation of absorption coefficients due to ice clouds.
         !$acc parallel loop gang collapse(2) async(async_id)
         do jj = 1, fulljj
            do k = 1, nlay ! lab_do_k
               !$acc loop vector private(cldice, refice, dgeice, factor, index, &
               !$acc&     fint, tausnw, tauran, tauliq, tauice)
               do iplon = 1, myim(jj) ! lab_do_iplon
                  if ( lcf1(iplon, jj) ) then
                     if (cfrac(iplon, k, jj) > cldmin) then ! lab_if_cld
                        cldliq = cliqp(iplon, k, jj)
                        !           refliq = max(2.5e0, min(60.0e0, reliq(k) ))
                        !           refice = max(5.0e0, reice(k) )
                        refliq = reliq(iplon, k, jj)

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


                        cldice = cicep(iplon, k, jj)
                        refice = reice(iplon, k, jj)
                        if ( cldice <= f_zero ) then
                           do ib = 1, nbands
                              tauice(ib) = f_zero
                           enddo
                        else
                           if ( ilwcice == 1 ) then
                              !  --- ...  ebert and curry approach for all particle sizes though somewhat
                              !           unjustified for large ice particles

                              refice = min(130.0, max(13.0, real(refice) ))

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

                              do ib = 1, nbands
                                 tauice(ib) = max(f_zero, cldice*(absice3(index,ib)    &
                                 &              + fint*(absice3(index+1,ib) - absice3(index,ib)) ))
                              enddo
                           end if

                        endif   ! end if_cldice_block
                        if (cdat3(iplon, k, jj)>f_zero .and. cdat4(iplon, k, jj)>10.0_kind_phys) then
                           tausnw = abssnow0*1.05756*cdat3(iplon, k, jj)/cdat4(iplon, k, jj)      ! fu's formula
                        else
                           tausnw = f_zero
                        endif
                        tauran = absrain * cdat1(iplon, k, jj)                      ! ncar formula
                        do ib = 1, nbands
                           taucld(iplon, k,ib, jj) = tauice(ib) + tauliq(ib) + tauran + tausnw
                        enddo
                     end if
                  end if
               end do
            end do
         end do

      else  ! lab_if_ilwcliq
         !$acc parallel loop gang collapse(2) async(async_id)
         do jj = 1, fulljj
            do k = 1, nlay
               !$acc loop vector 
               do iplon = 1, myim(jj) ! lab_do_iplon
                  if ( lcf1(iplon, jj) ) then
                     if (cfrac(iplon, k, jj) > cldmin) then
                        do ib = 1, nbands
                           taucld(iplon, k,ib, jj) = cdat1(iplon, k, jj)
                        enddo
                     endif
                  end if
               enddo
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
     &     ( pavel,tavel,tz,stemp,h2ovmr,colamt,coldry,colbrd,          &
     &       nlay, nlp1, ix, myim, async_id, fulljj,                                                 &
!  ---  outputs:
     &       laytrop,pklay,pklev,jp,jt,jt1,                             &
     &       rfrate,fac00,fac01,fac10,fac11,                            &
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
      integer, intent(in) :: nlay, nlp1, ix,  myim(fulljj), fulljj

      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj),intent(in) :: colamt
      real (kind=kind_phys), dimension(ix, 0:nlay, fulljj),     intent(in):: tz

      real (kind=kind_phys), dimension(ix, nlay, fulljj), intent(in) :: pavel,      &
     &       tavel, h2ovmr, coldry, colbrd

      real (kind=kind_phys), dimension(ix, fulljj), intent(in) :: stemp

!  ---  outputs:
      integer, dimension(ix, nlay, fulljj), intent(out) :: jp, jt, jt1, indself,    &
     &       indfor, indminor

      integer, dimension(ix, fulljj), intent(out) :: laytrop

      real (kind=kind_phys), dimension(ix, nlay,nrates,2, fulljj), intent(out) ::   &
     &       rfrate
      real (kind=kind_phys), dimension(ix,0:nlay, nbands, fulljj), intent(out) ::   &
     &       pklev, pklay

      real (kind=kind_phys), dimension(ix, nlay, fulljj),          intent(out) ::   &
     &       fac00, fac01, fac10, fac11, selffac, selffrac, forfac,     &
     &       forfrac, minorfrac, scaleminor, scaleminorn2

!  ---  locals:
      real (kind=kind_phys) :: tlvlfr, tlyrfr, plog, fp, ft, ft1,       &
     &       tem1, tem2, tavelr, tzr, pavelr, forfacr

      integer :: i, k, jp1, indlev, indlay, iplon, jj, async_id, jpr, jtr, jt1r, indminorr, tzi, taveli
!
!===> ... begin here
!
!  --- ...  calculate information needed by the radiative transfer routine
!           that is specific to this atmosphere, especially some of the
!           coefficients and indices needed to compute the optical depths
!           by interpolating data from stored reference atmospheres.

      !$acc parallel loop collapse(2) private(indlay, indlev, tlyrfr, tlvlfr, &
      !$acc&         tem1, tem2) async(async_id)
      do jj = 1, fulljj
         do iplon = 1, ix ! lab_do_iplon
            if (iplon .le. myim(jj)) then
               indlay = min(180, max(1, int(stemp(iplon, jj)-159.0) ))
               indlev = min(180, max(1, int(tz(iplon, 0, jj)-159.0) ))
               tlyrfr = stemp(iplon, jj) - int(stemp(iplon, jj))
               tlvlfr = tz(iplon, 0, jj) - int(tz(iplon, 0, jj))
               !$acc loop seq
               do i = 1, nbands
                  tem1 = totplnk(indlay+1,i) - totplnk(indlay,i)
                  tem2 = totplnk(indlev+1,i) - totplnk(indlev,i)
                  pklay(iplon, 0,i, jj) = delwave(i) * (totplnk(indlay,i) + tlyrfr*tem1)
                  pklev(iplon, 0,i, jj) = delwave(i) * (totplnk(indlev,i) + tlvlfr*tem2)
               enddo
            end if
         end do
      end do

            !  --- ...  begin layer loop
            !           calculate the integrated planck functions for each band at the
            !           surface, level, and layer temperatures.

      !$acc parallel loop collapse(2) private(jpr) async(async_id)
      do jj = 1, fulljj
         do iplon = 1, ix ! lab_do_iplon
            if (iplon .le. myim(jj)) then
               jpr = 0
               !$acc loop seq
               do k = 1, nlay
                  if (log(pavel(iplon, k, jj)) > 4.56) jpr = jpr + 1
               end do
               laytrop(iplon, jj) = jpr
            end if
         end do
      end do


      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !$acc loop vector private(tavelr, tzr, indlay, tlyrfr, &
            !$acc&         indlev, tlvlfr)
            do iplon = 1, myim(jj) ! lab_do_iplon
               !!$acc cache(delwave, totplnk)
               !  --- ...  begin spectral band loop
               tavelr = tavel(iplon, k, jj)
               tzr = tz(iplon, k, jj)
               indlay = min(180, max(1, int(tavelr-159.0) ))
               tlyrfr = tavelr - int(tavelr)
               
               indlev = min(180, max(1, int(tzr-159.0) ))
               tlvlfr = tzr - int(tzr)
               !$acc loop seq
               do i = 1, nbands
                  pklay(iplon, k,i, jj) = delwave(i) * (totplnk(indlay,i) + tlyrfr         &
                  &               * (totplnk(indlay+1,i) - totplnk(indlay,i)) )
                  pklev(iplon, k,i, jj) = delwave(i) * (totplnk(indlev,i) + tlvlfr         &
                  &               * (totplnk(indlev+1,i) - totplnk(indlev,i)) )
               enddo
            end do
         end do
      end do
      
      !$acc parallel loop collapse(3) private(plog, &
      !$acc&         jp1, fp, tem1, tem2, ft, ft1, tavelr, pavelr, jpr, jtr, jt1r, &
      !$acc&         forfacr, indminorr) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            do iplon = 1, ix
               !!$acc cache(chi_mls)
               if (iplon .le. myim(jj)) then ! lab_do_iplon
                  tavelr = tavel(iplon, k, jj)
                  pavelr = pavel(iplon, k, jj)



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
                  fac10(iplon, k, jj) = tem1 * ft
                  fac00(iplon, k, jj) = tem1 * (f_one - ft)
                  fac11(iplon, k, jj) = fp * ft1
                  fac01(iplon, k, jj) = fp * (f_one - ft1)

                  forfacr = pavelr*stpfac / (tavelr*(1.0 + h2ovmr(iplon, k, jj)))
                  selffac(iplon, k, jj) = h2ovmr(iplon, k, jj) * forfacr

                  !  --- ...  set up factors needed to separately include the minor gases
                  !           in the calculation of absorption coefficient

                  scaleminor(iplon, k, jj) = pavelr / tavelr
                  scaleminorn2(iplon, k, jj) = (pavelr / tavelr)                         &
                  &                  * (colbrd(iplon, k, jj)/(coldry(iplon, k, jj) + colamt(iplon, k,1, jj)))
                  tem1 = (tavelr - 180.8) / 7.2
                  indminorr = min(18, max(1, int(tem1)))
                  minorfrac(iplon, k, jj) = tem1 - float(indminorr)

                  !  --- ...  if the pressure is less than ~100mb, perform a different
                  !           set of species interpolations.
                  if (plog > 4.56) then

                     tem1 = (332.0 - tavel(iplon, k, jj)) / 36.0
                     indfor(iplon, k, jj) = min(2, max(1, int(tem1)))
                     forfrac(iplon, k, jj) = tem1 - float(indfor(iplon, k, jj))

                     !  --- ...  set up factors needed to separately include the water vapor
                     !           self-continuum in the calculation of absorption coefficient.

                     tem1 = (tavel(iplon, k, jj) - 188.0) / 7.2
                     indself(iplon, k, jj) = min(9, max(1, int(tem1)-7))
                     selffrac(iplon, k, jj) = tem1 - float(indself(iplon, k, jj) + 7)

                     !  --- ...  setup reference ratio to be used in calculation of binary
                     !           species parameter in lower atmosphere.

                     rfrate(iplon, k,1,1, jj) = chi_mls(1,jpr) / chi_mls(2,jpr)
                     rfrate(iplon, k,1,2, jj) = chi_mls(1,jpr+1) / chi_mls(2,jpr+1)

                     rfrate(iplon, k,2,1, jj) = chi_mls(1,jpr) / chi_mls(3,jpr)
                     rfrate(iplon, k,2,2, jj) = chi_mls(1,jpr+1) / chi_mls(3,jpr+1)

                     rfrate(iplon, k,3,1, jj) = chi_mls(1,jpr) / chi_mls(4,jpr)
                     rfrate(iplon, k,3,2, jj) = chi_mls(1,jpr+1) / chi_mls(4,jpr+1)

                     rfrate(iplon, k,4,1, jj) = chi_mls(1,jpr) / chi_mls(6,jpr)
                     rfrate(iplon, k,4,2, jj) = chi_mls(1,jpr+1) / chi_mls(6,jpr+1)

                     rfrate(iplon, k,5,1, jj) = chi_mls(4,jpr) / chi_mls(2,jpr)
                     rfrate(iplon, k,5,2, jj) = chi_mls(4,jpr+1) / chi_mls(2,jpr+1)

                  else

                     tem1 = (tavel(iplon, k, jj) - 188.0) / 36.0
                     indfor(iplon, k, jj) = 3
                     forfrac(iplon, k, jj) = tem1 - f_one

                     indself(iplon, k, jj) = 0
                     selffrac(iplon, k, jj) = f_zero

                     !  --- ...  setup reference ratio to be used in calculation of binary
                     !           species parameter in upper atmosphere.

                     rfrate(iplon, k,1,1, jj) = chi_mls(1,jpr) / chi_mls(2,jpr)
                     rfrate(iplon, k,1,2, jj) = chi_mls(1,jpr+1) / chi_mls(2,jpr+1)

                     rfrate(iplon, k,6,1, jj) = chi_mls(3,jpr) / chi_mls(2,jpr)
                     rfrate(iplon, k,6,2, jj) = chi_mls(3,jpr+1) / chi_mls(2,jpr+1)

                  endif

                  !  --- ...  rescale selffac and forfac for use in taumol

                  selffac(iplon, k, jj) = colamt(iplon, k,1, jj) * selffac(iplon, k, jj)
                  forfacr = colamt(iplon, k,1, jj) * forfacr
                  jp(iplon, k, jj) = jpr
                  jt(iplon, k, jj) = jtr
                  jt1(iplon, k, jj) = jt1r
                  forfac(iplon, k, jj) = forfacr
                  indminor(iplon, k, jj) = indminorr
               end if

            enddo   ! end do_k layer loop
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
     &     ( semiss,delp,cldfrc,taucld,tautot,pklay,pklev,              &
     &       fracs,secdif, nlay,nlp1, ix, myim, async_id, fulljj,                                   &
!  ---  outputs:
     &       totuflux,totdflux,htr, totuclfl,totdclfl,htrcl, htrb       &
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

!  ---  inputs:
      integer, intent(in) :: nlay, nlp1, ix, myim(fulljj), fulljj

      real (kind=kind_phys), dimension(ix, 0:nlp1, fulljj), intent(in) :: cldfrc
      real (kind=kind_phys), dimension(ix, nbands, fulljj), intent(in) :: semiss,   &
     &       secdif
      real (kind=kind_phys), dimension(ix, nlay, fulljj),   intent(in) :: delp

      real (kind=kind_phys), dimension(ix, nlay, nbands, fulljj),intent(in):: taucld
      real (kind=kind_phys), dimension(ix, nlay, ngptlw, fulljj),intent(in):: fracs, &
     &       tautot

      real (kind=kind_phys), dimension(ix,0:nlay, nbands, fulljj), intent(in) ::    &
     &       pklev, pklay

!  ---  outputs:
      real (kind=kind_phys), dimension(ix, nlay, fulljj), intent(out) :: htr, htrcl

      real (kind=kind_phys), dimension(ix, nlay,nbands, fulljj),intent(out) :: htrb

      real (kind=kind_phys), dimension(ix, 0:nlay, fulljj), intent(out) ::          &
     &       totuflux, totdflux, totuclfl, totdclfl

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

      integer :: ittot, itgas, ib, ig, k, iplon, jj, iplon2, i2, jj2, ng00

!  dimensions for cloud overlap adjustment
      real (kind=kind_phys), dimension(ix, nlp1, fulljj) :: faccld1u, faccld2u,     &
     &        facclr1u, facclr2u, faccmb1u, faccmb2u
      real (kind=kind_phys), dimension(ix, 0:nlay, fulljj) :: faccld1d, faccld2d,   &
     &        facclr1d, facclr2d, faccmb1d, faccmb2d
     real (kind=kind_phys), dimension(ix, 0:nlay, nbands, fulljj) :: toturad, totdrad

      integer :: async_id
      real (kind=kind_phys), dimension(:,:,:,:), allocatable :: gassrcu, totsrcu, &
         trngas, trntot, radtotd_2, radtotu_2
      
      ! GPU register
      logical :: lstcldr
      real (kind=kind_phys) :: fnet, fnet1, totufluxr, totdfluxr, fnet3, fnet4, clfrm, clfrp
      real (kind=kind_phys), dimension(ix, 0:nlay, nbands, fulljj) :: fnet2

      integer, dimension(nbands) :: ng_array, ns_array
      data ng_array(:) /ng01, ng02, ng03, ng04, ng05, ng06, ng07, ng08, &
                        ng09, ng10, ng11, ng12, ng13, ng14, ng15, ng16/ 
      data ns_array(:) /ns01, ns02, ns03, ns04, ns05, ns06, ns07, ns08, &
                        ns09, ns10, ns11, ns12, ns13, ns14, ns15, ns16/ 

!
!===> ...  begin here
!
      !small_fulljj = ceiling(float(jlistnum)/float(small_factor))
      ng00 = maxval(ng_array)
      allocate(gassrcu(ix, nlay, ng00, fulljj))
      allocate(totsrcu(ix, nlay, ng00, fulljj))
      allocate(trngas(ix, nlay, ng00, fulljj))
      allocate(trntot(ix, nlay, ng00, fulljj))
      allocate(radtotd_2(ix, 0:nlay, ng00, fulljj))
      allocate(radtotu_2(ix, 0:nlay, ng00, fulljj))
      
      !$acc data create(faccld1u, faccld2u, facclr1u, facclr2u, faccmb1u, &
      !$acc&     faccmb2u, faccld1d, faccld2d, facclr1d, facclr2d, faccmb1d, &
      !$acc&     faccmb2d, toturad, totdrad, fnet2, &
      !$acc&     gassrcu, totsrcu, &
      !$acc&     trngas, trntot, radtotd_2, radtotu_2) async(async_id)

      !$acc parallel loop gang collapse(2) private(iplon) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlp1
            !$acc loop vector
            do iplon = 1, ix ! lab_do_iplon
               if (iplon .le. myim(jj)) then
                  faccld1u(iplon, k, jj) = f_zero
                  faccld2u(iplon, k, jj) = f_zero
                  facclr1u(iplon, k, jj) = f_zero
                  facclr2u(iplon, k, jj) = f_zero
                  faccmb1u(iplon, k, jj) = f_zero
                  faccmb2u(iplon, k, jj) = f_zero
               end if
            enddo
         end do
      end do
      
      !$acc parallel loop gang collapse(2) private(iplon) async(async_id)
      do jj = 1, fulljj
         do k = 0, nlay
            !$acc loop vector
            do iplon = 1, ix ! lab_do_iplon
               if (iplon .le. myim(jj)) then
                  faccld1d(iplon, k, jj) = f_zero
                  faccld2d(iplon, k, jj) = f_zero
                  facclr1d(iplon, k, jj) = f_zero
                  facclr2d(iplon, k, jj) = f_zero
                  faccmb1d(iplon, k, jj) = f_zero
                  faccmb2d(iplon, k, jj) = f_zero
               end if
            enddo
         end do
      end do

      !$acc parallel loop collapse(2) private(iplon, rat1, rat2, fmax, fmin, &
      !$acc&         lstcldr, clfr, clfrp, clfrm) async(async_id)
      do jj = 1, fulljj
         do iplon = 1, ix ! lab_do_iplon
            if (iplon .le. myim(jj)) then
               clfr = cldfrc(iplon, 1, jj)
               clfrm = cldfrc(iplon, 0, jj)
               lstcldr = clfr > eps
               rat1 = f_zero
               rat2 = f_zero
               !$acc loop seq
               do k = 1, nlay-1
                  clfrp = cldfrc(iplon, k+1, jj)
                  if (clfr > eps) then
                  !  --- ...  maximum/random cloud overlap

                     if (clfrp >= clfr) then
                        if (lstcldr) then
                           if (clfr < f_one) then
                              facclr2u(iplon, k+1, jj) = (clfrp - clfr)               &
                              &                        / (f_one - clfr)
                           endif
                           facclr2u(iplon, k, jj) = f_zero
                           faccld2u(iplon, k, jj) = f_zero
                        else
                           fmax = max(clfr, clfrm)
                           if (clfrp > fmax) then
                              facclr1u(iplon, k+1, jj) = rat2
                              facclr2u(iplon, k+1, jj) = (clfrp - fmax)/(f_one - fmax)
                           elseif (clfrp < fmax) then
                              facclr1u(iplon, k+1, jj) = (clfrp - clfr)               &
                              &                        / (clfrm - clfr)
                           else
                              facclr1u(iplon, k+1, jj) = rat2
                           endif
                        endif

                        if (facclr1u(iplon, k+1, jj)>f_zero .or. facclr2u(iplon, k+1, jj)>f_zero) then
                           rat1 = f_one
                           rat2 = f_zero
                        else
                           rat1 = f_zero
                           rat2 = f_zero
                        endif
                     else
                        if (lstcldr) then
                           faccld2u(iplon, k+1, jj) = (clfr - &
                                                         clfrp) / clfr
                           facclr2u(iplon, k, jj) = f_zero
                           faccld2u(iplon, k, jj) = f_zero
                        else
                           fmin = min(clfr, clfrm)
                           if (clfrp <= fmin) then
                              faccld1u(iplon, k+1, jj) = rat1
                              faccld2u(iplon, k+1, jj) = (fmin - clfrp) / fmin
                           else
                              faccld1u(iplon, k+1, jj) = (clfr - clfrp)               &
                              &                        / (clfr - fmin)
                           endif
                        endif

                        if (faccld1u(iplon, k+1, jj)>f_zero .or. faccld2u(iplon, k+1, jj)>f_zero) then
                           rat1 = f_zero
                           rat2 = f_one
                        else
                           rat1 = f_zero
                           rat2 = f_zero
                        endif
                     endif

                     faccmb1u(iplon, k+1, jj) = facclr1u(iplon, k+1, jj) * &
                                                faccld2u(iplon, k, jj) * clfrm
                     faccmb2u(iplon, k+1, jj) = faccld1u(iplon, k+1, jj) * facclr2u(iplon, k, jj)                   &
                     &                  * (f_one - clfrm)
                  endif
                  lstcldr = clfrp>eps .and. clfr<=eps
                  clfrm = clfr
                  clfr = clfrp
               enddo
               
               clfr = cldfrc(iplon, nlay, jj)
               clfrp = cldfrc(iplon, nlay+1, jj)
               lstcldr = clfr > eps
               rat1 = f_zero
               rat2 = f_zero
               !$acc loop seq
               do k = nlay, 2, -1
                  clfrm = cldfrc(iplon, k-1, jj)
                  if (clfr > eps) then

                     if (clfrm >= clfr) then
                        if (lstcldr) then
                           if (clfr < f_one) then
                              facclr2d(iplon, k-1, jj) = (clfrm - clfr)               &
                              &                        / (f_one - clfr)
                           endif

                           facclr2d(iplon, k, jj) = f_zero
                           faccld2d(iplon, k, jj) = f_zero
                        else
                           fmax = max(clfr, clfrp)

                           if (clfrm > fmax) then
                              facclr1d(iplon, k-1, jj) = rat2
                              facclr2d(iplon, k-1, jj) = (clfrm - fmax) / (f_one - fmax)
                           elseif (clfrm < fmax) then
                              facclr1d(iplon, k-1, jj) = (clfrm - clfr)               &
                              &                        / (clfrp - clfr)
                           else
                              facclr1d(iplon, k-1, jj) = rat2
                           endif
                        endif

                        if (facclr1d(iplon, k-1, jj)>f_zero .or. facclr2d(iplon, k-1, jj)>f_zero) then
                           rat1 = f_one
                           rat2 = f_zero
                        else
                           rat1 = f_zero
                           rat2 = f_zero
                        endif
                     else
                        if (lstcldr) then
                           faccld2d(iplon, k-1, jj) = (clfr - &
                                                      clfrm) / clfr
                           facclr2d(iplon, k, jj) = f_zero
                           faccld2d(iplon, k, jj) = f_zero
                        else
                           fmin = min(clfr, clfrp)

                           if (clfrm <= fmin) then
                              faccld1d(iplon, k-1, jj) = rat1
                              faccld2d(iplon, k-1, jj) = (fmin - clfrm) / fmin
                           else
                              faccld1d(iplon, k-1, jj) = (clfr - clfrm)               &
                              &                        / (clfr - fmin)
                           endif
                        endif

                        if (faccld1d(iplon, k-1, jj)>f_zero .or. faccld2d(iplon, k-1, jj)>f_zero) then
                           rat1 = f_zero
                           rat2 = f_one
                        else
                           rat1 = f_zero
                           rat2 = f_zero
                        endif
                     endif

                     faccmb1d(iplon, k-1, jj) = facclr1d(iplon, k-1, jj) * &
                                                faccld2d(iplon, k, jj) * clfrp
                     faccmb2d(iplon, k-1, jj) = faccld1d(iplon, k-1, jj) * facclr2d(iplon, k, jj)                   &
                     &                  * (f_one - clfrp)
                  endif
                  lstcldr = clfrm > eps .and. clfr<=eps
                  clfrp = clfr
                  clfr = clfrm
               enddo
            end if
         end do
      end do
      
      
      !  --- ...  initialize for radiative transfer.
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do ib = 1, nbands
            !$acc loop vector
            do iplon = 1, myim(jj) ! lab_do_iplon
               totdrad(iplon, nlay,ib, jj) = f_zero
            enddo
         end do
      end do

      !$acc parallel loop gang collapse(2) private(iplon) async(async_id)
      do jj = 1, fulljj
         do k = 0, nlay
            !$acc loop vector
            do iplon = 1, ix ! lab_do_iplon
               if (iplon .le. myim(jj)) then
                  totuclfl(iplon, k, jj) = f_zero
                  totdclfl(iplon, k, jj) = f_zero
               end if
            enddo
         end do
      end do
      

      do ib = 1, nbands
         call rtrnmr_ngptlw(secdif, tautot, fracs, pklay, pklev, &
            totdrad, totdclfl, facclr1d, faccld1d, faccmb1d, &
            faccmb2d, facclr2d, faccld2d, semiss, toturad, totuclfl, cldfrc, &
            facclr1u, faccld1u, faccmb1u, faccmb2u, facclr2u, faccld2u, taucld, &
            nlay, myim, ix, ngptlw, nlp1, nbands, &
            ns_array(ib), ng_array(ib), ng00, async_id, fulljj, &
            gassrcu, totsrcu, &
            trngas, trntot, radtotd_2, radtotu_2)
      end do

      flxfac = wtdiff * fluxfac
      !  --- ...  process longwave output from band for total and clear streams.
      !           calculate upward, downward, and net flux.
      !$acc parallel loop collapse(2) async(async_id) 
      do jj = 1, fulljj
         do k = 0, nlay
            !$acc loop vector private(totufluxr, totdfluxr)
            do iplon = 1, myim(jj) ! lab_do_iplon
               totufluxr = f_zero
               totdfluxr = f_zero
               !$acc loop seq
               do ib = 1, nbands
                  totufluxr = totufluxr + toturad(iplon, k,ib, jj)
                  totdfluxr = totdfluxr + totdrad(iplon, k,ib, jj)
               enddo
               totuflux(iplon, k, jj) = totufluxr * flxfac
               totdflux(iplon, k, jj) = totdfluxr * flxfac
            end do
         end do
      end do
      
      !$acc parallel loop collapse(2) private(fnet, fnet1, rfdelp) async(async_id)
      do jj = 1, fulljj
         do iplon = 1, ix ! lab_do_iplon
            if (iplon .le. myim(jj)) then
               !  --- ...  calculate net fluxes and heating rates
               fnet1 = totuflux(iplon, 0, jj) - totdflux(iplon, 0, jj)
               !$acc loop seq
               do k = 1, nlay
                  rfdelp = heatfac / delp(iplon, k, jj)
               enddo
            end if
         end do
      end do
      
      !$acc parallel loop collapse(2) private(fnet, fnet1, rfdelp, fnet4, fnet3) async(async_id)
      do jj = 1, fulljj
         do iplon = 1, ix ! lab_do_iplon
            if (iplon .le. myim(jj)) then
               !! --- ...  optional clear sky heating rates
               if ( lhlw0 )fnet1 = totuclfl(iplon, 0, jj) - totdclfl(iplon, 0, jj)
               fnet4 = totuflux(iplon, 0, jj) - totdflux(iplon, 0, jj)
               !$acc loop seq
               do k = 1, nlay
                  rfdelp = heatfac / delp(iplon, k, jj)
                  if ( lhlw0 ) then
                     fnet = totuclfl(iplon, k, jj) - totdclfl(iplon, k, jj)
                     htrcl(iplon, k, jj) = (fnet1 - fnet) * rfdelp
                     fnet1 = fnet
                  end if
                  fnet3 = totuflux(iplon, k, jj) - totdflux(iplon, k, jj)
                  htr (iplon, k, jj) = (fnet4 - fnet3) * rfdelp
                  fnet4 = fnet3
               enddo
            endif
         end do
      end do
         
      !! --- ...  optional spectral band heating rates
      if ( lhlwb ) then
         !$acc parallel loop collapse(3) private(rfdelp) async(async_id)
         do jj = 1, fulljj
            do ib = 1, nbands
               do iplon = 1, ix ! lab_do_iplon
                  if (iplon .le. myim(jj)) then
                     fnet2(iplon, 0, ib, jj) = (toturad(iplon, 0,ib, jj) - totdrad(iplon, 0,ib, jj)) * flxfac
                     !$acc loop seq
                     do k = 1, nlay
                        rfdelp = heatfac / delp(iplon, k, jj)
                        fnet2(iplon, k, ib, jj) = (toturad(iplon, k,ib, jj) - totdrad(iplon, k,ib, jj)) * flxfac
                        htrb(iplon, k,ib, jj) = (fnet2(iplon, k-1, ib, jj) - fnet2(iplon, k, ib, jj)) * rfdelp
                     enddo
                  endif
               enddo
            end do
         end do
      end if

      !$acc end data
      deallocate(gassrcu, totsrcu, &
         trngas, trntot, radtotd_2, radtotu_2)

! .................................
      end subroutine rtrnmr
! ---------------------------------

      subroutine rtrnmr_ngptlw(secdif, tautot, fracs, pklay, pklev, &
         totdrad, totdclfl, facclr1d, faccld1d, faccmb1d, &
         faccmb2d, facclr2d, faccld2d, semiss, toturad, totuclfl, cldfrc, &
         facclr1u, faccld1u, faccmb1u, faccmb2u, facclr2u, faccld2u, taucld, &
         nlay, myim, ix, ngptlw, nlp1, nbands, nslw, nglw, ng00, async_id, fulljj, &
         gassrcu, totsrcu, &
         trngas, trntot, radtotd_2, radtotu_2)

      implicit none
      integer :: nlay, nlp1, ix, myim(fulljj), ngptlw, nbands, nslw, nglw, fulljj

      real (kind=kind_phys), dimension(ix, 0:nlp1, fulljj) :: cldfrc
      real (kind=kind_phys), dimension(ix, nbands, fulljj) :: semiss,   &
     &       secdif

      real (kind=kind_phys), dimension(ix, nlay, ngptlw, fulljj) :: fracs, &
     &       tautot

      real (kind=kind_phys), dimension(ix,0:nlay, nbands, fulljj) ::    &
     &       pklev, pklay
      real (kind=kind_phys), dimension(ix, nlay, nbands, fulljj) :: taucld
      real (kind=kind_phys), parameter :: rec_6 = 0.166667

!  ---  locals:


      real (kind=kind_phys),dimension(ix, nlay, ng00, fulljj) :: &
         gassrcu, totsrcu, trngas, trntot
      real (kind=kind_phys),dimension(ix, 0:nlay, ng00, fulljj) :: &
         radtotd_2, radtotu_2

      real (kind=kind_phys) :: totsrcd, gassrcd, tblind, odepth, odtot, &
     &       odcld, atrtot, atrgas, reflct, totfac, gasfac,     &
     &       plfrac, blay, bbdgas, bbdtot, bbugas, bbutot, dplnku,      &
     &       dplnkd, radtotu, radclru, rad0, rad,     &
     &       totradd, clrradd, totradu, clrradu, &
     &       radmod, clfr, trng, trnt, gasu, totu

      integer :: ittot, itgas, ib, ig, k, iplon, jj, iplon2, i2, jj2, ir, ng00

!  dimensions for cloud overlap adjustment
      real (kind=kind_phys), dimension(ix, nlp1, fulljj) :: faccld1u, faccld2u,     &
     &        facclr1u, facclr2u, faccmb1u, faccmb2u
      real (kind=kind_phys), dimension(ix, 0:nlay, fulljj) :: faccld1d, faccld2d,   &
     &        facclr1d, facclr2d, faccmb1d, faccmb2d, totdclfl, totuclfl
      real (kind=kind_phys), dimension(ix, 0:nlay, nbands, fulljj) :: toturad, totdrad

      real(kind=kind_phys) :: radtotd, radclrd, pklevr1, pklevr, clfr1, semissr, &
         secdifr, radtotdr, radtotur, flxfac
      logical :: lstcldr
      integer :: async_id
      
      ib = ngb(nslw + 1)
      flxfac = wtdiff * fluxfac
      !$acc parallel loop gang collapse(3) private(odepth, atrgas, trng, &
      !$acc&         gasfac, tblind, itgas, plfrac, blay, radtotd, radclrd, &
      !$acc&         dplnku, dplnkd, bbdgas, bbugas, gassrcd, clfr, totradd, &
      !$acc&         clrradd, rad, odcld, odtot, totfac, atrtot, trnt, ittot, &
      !$acc&         bbdtot, bbutot, totsrcd, radmod, reflct, rad0, radtotu, &
      !$acc&         radclru, gasu, totradu, clrradu, totu, ir, pklevr1, pklevr, &
      !$acc&         clfr1, lstcldr, secdifr, semissr) async(async_id)
      do jj = 1, fulljj
         do ig = nslw + 1, nslw + nglw
            do iplon = 1, ix
               if (iplon .le. myim(jj)) then ! lab_do_iplon
                  ir = ig - nslw
                  !  --- ...  loop over all g-points
                  radtotd = f_zero
                  radclrd = f_zero
                  !  --- ...  downward radiative transfer loop.
                  pklevr = pklev(iplon, nlay,ib, jj)

                  clfr = cldfrc(iplon, nlay, jj)
                  lstcldr = clfr > eps
                  secdifr = secdif(iplon, ib, jj)
                  !$acc loop seq
                  do k = nlay, 1, -1
                     pklevr1 = pklev(iplon, k-1,ib, jj)

                  !  --- ...  clear sky, gases contribution

                     odepth = max( f_zero, secdifr*tautot(iplon, k,ig, jj) )
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

                     plfrac = fracs(iplon, k,ig, jj)
                     blay = pklay(iplon, k,ib, jj)

                     dplnku = pklevr - blay
                     dplnkd = pklevr1 - blay
                     bbdgas = plfrac * (blay + dplnkd*gasfac)
                     bbugas = plfrac * (blay + dplnku*gasfac)
                     gassrcd   = bbdgas * atrgas
                     gassrcu(iplon, k, ir, jj)= bbugas * atrgas
                     trngas(iplon, k, ir, jj) = trng
                     pklevr = pklevr1

                     !  --- ...  total sky, gases+clouds contribution

                     if (lstcldr) then
                        totradd = clfr * radtotd
                        clrradd = radtotd - totradd
                        rad = f_zero
                     endif
                     clfr1 = cldfrc(iplon, k-1, jj)
                     lstcldr = clfr1 > eps .and. clfr<=eps


                     if (clfr >= eps) then
                        !  --- ...  cloudy layer

                        odcld = secdifr * taucld(iplon, k,ib, jj)
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
                        totsrcu(iplon, k, ir, jj)= bbutot * atrtot
                        trntot(iplon, k, ir, jj) = trnt

                        totradd = totradd*trnt + clfr*totsrcd
                        clrradd = clrradd*trng + (f_one - clfr)*gassrcd

                        !  --- ...  total sky radiance
                        radtotd = totradd + clrradd
                        radtotd_2(iplon, k, ir, jj) = radtotd

                        !  --- ...  clear sky radiance
                        radclrd = radclrd*trng + gassrcd
                        !$acc atomic
                        totdclfl(iplon, k-1, jj) = totdclfl(iplon, k-1, jj) + radclrd * flxfac

                        radmod = rad*(facclr1d(iplon, k-1, jj)*trng + faccld1d(iplon, k-1, jj)*trnt)      &
                        &             - faccmb1d(iplon, k-1, jj)*gassrcd + faccmb2d(iplon, k-1, jj)*totsrcd

                        rad = -radmod + facclr2d(iplon, k-1, jj)*(clrradd + radmod)            &
                        &                    - faccld2d(iplon, k-1, jj)*(totradd - radmod)
                        totradd = totradd + rad
                        clrradd = clrradd - rad

                     else
                        !  --- ...  clear layer

                        !  --- ...  total sky radiance
                        radtotd = radtotd*trng + gassrcd
                        radtotd_2(iplon, k, ir, jj) = radtotd

                        !  --- ...  clear sky radiance
                        radclrd = radclrd*trng + gassrcd
                        !$acc atomic
                        totdclfl(iplon, k-1, jj) = totdclfl(iplon, k-1, jj) + radclrd * flxfac

                     endif   ! end if_clfr_block
                     clfr = clfr1

                  enddo   ! end do_k_loop

                  !  --- ...  spectral emissivity & reflectance
                  !           include the contribution of spectrally varying longwave emissivity
                  !           and reflection from the surface to the upward radiative transfer.
                  !     note: spectral and lambertian reflection are identical for the
                  !           diffusivity angle flux integration used here.
                  semissr = semiss(iplon, ib, jj)
                  reflct = f_one - semissr
                  rad0 = semissr * fracs(iplon, 1,ig, jj) * pklevr1

                  !  --- ...  total sky radiance
                  radtotu = rad0 + reflct*radtotd
                  radtotu_2(iplon, 0, ir, jj) = radtotu

                  !  --- ...  clear sky radiance
                  radclru = rad0 + reflct*radclrd
                  !$acc atomic
                  totuclfl(iplon, 0, jj) = totuclfl(iplon, 0, jj) + radclru * flxfac

                  !  --- ...  upward radiative transfer loop.
                  clfr = cldfrc(iplon, 1, jj)
                  lstcldr = clfr > eps
                  !$acc loop seq
                  do k = 1, nlay

                     trng = trngas(iplon, k, ir, jj)
                     gasu = gassrcu(iplon, k, ir, jj)

                     if (lstcldr) then
                        totradu = clfr * radtotu
                        clrradu = radtotu - totradu
                        rad = f_zero
                     endif
                     clfr1 = cldfrc(iplon, k+1, jj)
                     lstcldr = clfr1>eps .and. clfr<=eps


                     if (clfr >= eps) then
                        !  --- ...  cloudy layer

                        trnt = trntot(iplon, k, ir, jj)
                        totu = totsrcu(iplon, k, ir, jj)
                        totradu = totradu*trnt + clfr*totu
                        clrradu = clrradu*trng + (f_one - clfr)*gasu

                        !  --- ...  total sky radiance
                        radtotu = totradu + clrradu
                        radtotu_2(iplon, k, ir, jj) = radtotu

                        !  --- ...  clear sky radiance
                        radclru = radclru*trng + gasu
                        !$acc atomic
                        totuclfl(iplon, k, jj) = totuclfl(iplon, k, jj) + radclru * flxfac

                        radmod = rad*(facclr1u(iplon, k+1, jj)*trng + faccld1u(iplon, k+1, jj)*trnt)      &
                        &             - faccmb1u(iplon, k+1, jj)*gasu + faccmb2u(iplon, k+1, jj)*totu
                        rad = -radmod + facclr2u(iplon, k+1, jj)*(clrradu + radmod)            &
                        &                    - faccld2u(iplon, k+1, jj)*(totradu - radmod)
                        totradu = totradu + rad
                        clrradu = clrradu - rad

                     else
                        !  --- ...  clear layer

                        !  --- ...  total sky radiance
                        radtotu = radtotu*trng + gasu
                        radtotu_2(iplon, k, ir, jj) = radtotu

                        !  --- ...  clear sky radiance
                        radclru = radclru*trng + gasu
                        !$acc atomic
                        totuclfl(iplon, k, jj) = totuclfl(iplon, k, jj) + radclru * flxfac

                     endif   ! end if_clfr_block
                     clfr = clfr1
                  enddo   ! end do_k_loop
               end if
            enddo   ! end do_ig_loop
         end do
      end do

      !$acc parallel loop collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 0, nlay
            !$acc loop vector private(ir, radtotur, radtotdr)
            do iplon = 1, ix ! lab_do_iplon
               if (iplon .le. myim(jj)) then
                  radtotdr = 0.
                  radtotur = 0.
                  !$acc loop seq
                  do ig = nslw + 1, nslw + nglw
                     ir = ig - nslw
                     if (k .ge. 1) radtotdr = radtotdr + radtotd_2(iplon, k, ir, jj)
                     radtotur = radtotur + radtotu_2(iplon, k, ir, jj)
                  end do
                  toturad(iplon, k,ib, jj) =  radtotur
                  if (k .ge. 1) totdrad(iplon, k-1, ib, jj) = radtotdr
               end if
            end do
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
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, async_id, fulljj,                                                         &
!  ---  outputs:
     &       fracs, tautot                                              &
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
      integer, intent(in) :: nlay, laytrop(ix, fulljj), ix, myim(fulljj), async_id, fulljj

      integer, dimension(ix, nlay, fulljj), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(ix, nlay, fulljj), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2

      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj), intent(in):: colamt
      real (kind=kind_phys), dimension(ix, nlay,maxxsec, fulljj),intent(in):: wx

      real (kind=kind_phys), dimension(ix, nlay, nbands, fulljj), intent(in):: tauaer

      real (kind=kind_phys), dimension(ix, nlay,nrates,2, fulljj), intent(in) ::    &
     &       rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(ix, nlay, ngptlw, fulljj), intent(out) ::     &
     &       fracs, tautot

!  ---  locals
      integer :: small_factor = 1
      integer :: small_ix
      real (kind=kind_phys) :: taug 
      !integer, dimension(:,:,:,:), allocatable :: mask

      integer :: ib, ig, k, jj, i2


      small_ix = ceiling(float(ix)/float(small_factor))
      !allocate(mask(small_ix, nlay, ngptlw, fulljj))
      !mask = 0
!
!===> ...  begin here
!
      do i2 = 1, small_factor
         call taugb01 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )

         call taugb02 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )

         call taugb03 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )

         call taugb04 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )

         call taugb05 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )

         call taugb06 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )

         call taugb07 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )

         call taugb08 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )

         call taugb09 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )

         call taugb10 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )

         call taugb11 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )

         call taugb12 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )

         call taugb13 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )

         call taugb14 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )

         call taugb15 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )

         call taugb16 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )

         !write(*,*) "ARA", myrank, sum(mask)
      end do

! ..................................
      end subroutine taumol
!-----------------------------------

! ----------------------------------
      subroutine taugb01 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
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
      integer, intent(in) :: nlay, laytrop(ix, fulljj), ix, myim(fulljj), &
         async_id, small_ix, i2, fulljj

      integer, dimension(ix, nlay, fulljj), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(ix, nlay, fulljj), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj), intent(in):: colamt
      real (kind=kind_phys), dimension(ix, nlay,maxxsec, fulljj),intent(in):: wx

      real (kind=kind_phys), dimension(ix, nlay, nbands, fulljj), intent(in):: tauaer

      real (kind=kind_phys), dimension(ix, nlay,nrates,2, fulljj), intent(in) ::    &
     &       rfrate
!  ---  outputs:
      real (kind=kind_phys), dimension(ix, nlay, ngptlw, fulljj), intent(out) ::     &
     &       fracs, tautot

!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &           indm, indmp, ig, ib

      real (kind=kind_phys) :: pp, corradj, scalen2, tauself, taufor,   &
     &       taun2, taug
      integer :: jj, iplon, iplon2
!
!===> ...  begin here
!
!  ---  minor gas mapping levels:
!     lower - n2, p = 142.5490 mbar, t = 215.70 k
!     upper - n2, p = 142.5490 mbar, t = 215.70 k

      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !$acc loop vector private(ind0, ind1, inds, indf, indm, &
            !$acc&         ind0p, ind1p, indsp, indfp, indmp, pp, scalen2, corradj, &
            !$acc&         tauself, taufor, taun2, iplon, taug)
            do iplon = 1, myim(jj)
               if (k .le. laytrop(iplon, jj)) then
                  !  --- ...  lower atmosphere loop
                  ind0 = ((jp(iplon, k, jj)-1)*5 + (jt (iplon, k, jj)-1)) * nspa(1) + 1
                  ind1 = ( jp(iplon, k, jj)   *5 + (jt1(iplon, k, jj)-1)) * nspa(1) + 1
                  inds = indself(iplon, k, jj)
                  indf = indfor(iplon, k, jj)
                  indm = indminor(iplon, k, jj)

                  ind0p = ind0 + 1
                  ind1p = ind1 + 1
                  indsp = inds + 1
                  indfp = indf + 1
                  indmp = indm + 1

                  pp = pavel(iplon, k, jj)
                  scalen2 = colbrd(iplon, k, jj) * scaleminorn2(iplon, k, jj)
                  if (pp < 250.0) then
                     corradj = f_one - 0.15 * (250.0-pp) / 154.4
                  else
                     corradj = f_one
                  endif
                  !$acc loop seq
                  do ig = 1, ng01
                     ib = ngb(ig)
                     tauself = selffac(iplon, k, jj) * (selfref(ig,inds) + selffrac(iplon, k, jj)        &
                     &            * (selfref(ig,indsp) - selfref(ig,inds)))
                     taufor  = forfac(iplon, k, jj) * (forref(ig,indf) + forfrac(iplon, k, jj)           &
                     &            * (forref(ig,indfp) -  forref(ig,indf))) 
                     taun2   = scalen2 * (ka_mn2(ig,indm) + minorfrac(iplon, k, jj)           &
                     &            * (ka_mn2(ig,indmp) - ka_mn2(ig,indm)))

                     taug = corradj * (colamt(iplon, k,1, jj)                           &
                     &            * (fac00(iplon, k, jj)*absa(ig,ind0) + fac10(iplon, k, jj)*absa(ig,ind0p)   &
                     &            +  fac01(iplon, k, jj)*absa(ig,ind1) + fac11(iplon, k, jj)*absa(ig,ind1p))  &
                     &            + tauself + taufor + taun2)

                     fracs(iplon, k, ig, jj) = fracrefa(ig)
                     tautot(iplon, k, ns01+ig, jj) = taug + tauaer(iplon, k, ib, jj)
                  enddo
               else
                  !  --- ...  upper atmosphere loop
                  ind0 = ((jp(iplon, k, jj)-13)*5 + (jt (iplon, k, jj)-1)) * nspb(1) + 1
                  ind1 = ((jp(iplon, k, jj)-12)*5 + (jt1(iplon, k, jj)-1)) * nspb(1) + 1
                  indf = indfor(iplon, k, jj)
                  indm = indminor(iplon, k, jj)

                  ind0p = ind0 + 1
                  ind1p = ind1 + 1
                  indfp = indf + 1
                  indmp = indm + 1

                  scalen2 = colbrd(iplon, k, jj) * scaleminorn2(iplon, k, jj)
                  corradj = f_one - 0.15 * (pavel(iplon, k, jj) / 95.6)
                  !$acc loop seq
                  do ig = 1, ng01
                     ib = ngb(ig)
                     taufor = forfac(iplon, k, jj) * (forref(ig,indf) + forfrac(iplon, k, jj)            &
                     &           * (forref(ig,indfp) - forref(ig,indf))) 
                     taun2  = scalen2 * (kb_mn2(ig,indm) + minorfrac(iplon, k, jj)            &
                     &           * (kb_mn2(ig,indmp) - kb_mn2(ig,indm)))

                     taug = corradj * (colamt(iplon, k,1, jj)                           &
                     &           * (fac00(iplon, k, jj)*absb(ig,ind0) + fac10(iplon, k, jj)*absb(ig,ind0p)    &
                     &           +  fac01(iplon, k, jj)*absb(ig,ind1) + fac11(iplon, k, jj)*absb(ig,ind1p))   &
                     &           + taufor + taun2)

                     fracs(iplon, k, ig, jj) = fracrefb(ig)
                     tautot(iplon, k, ns01+ig, jj) = taug + tauaer(iplon, k, ib, jj)
                  end do
               end if
            end do
         end do
      end do


! ..................................
      end subroutine taugb01
! ----------------------------------

! ----------------------------------
      subroutine taugb02 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!     band 2:  350-500 cm-1 (low key - h2o; high key - h2o)            !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb02
      integer, intent(in) :: nlay, laytrop(ix, fulljj), ix, myim(fulljj), &
         async_id, small_ix, i2, fulljj

      integer, dimension(ix, nlay, fulljj), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(ix, nlay, fulljj), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj), intent(in):: colamt
      real (kind=kind_phys), dimension(ix, nlay,maxxsec, fulljj),intent(in):: wx

      real (kind=kind_phys), dimension(ix, nlay, nbands, fulljj), intent(in):: tauaer

      real (kind=kind_phys), dimension(ix, nlay,nrates,2, fulljj), intent(in) ::    &
     &       rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(ix, nlay, ngptlw, fulljj), intent(out) ::     &
     &       fracs, tautot
!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &           ig, ib

      real (kind=kind_phys) :: corradj, tauself, taufor, taug
      integer :: jj, iplon, iplon2
!
!===> ...  begin here
!
     !$acc parallel loop gang collapse(2) async(async_id)
     do jj = 1, fulljj
         do k = 1, nlay
            !$acc loop vector private(iplon, ind0, ind1, inds, indf, ind0p, &
            !$acc&         ind1p, indsp, indfp, corradj, tauself, taufor, taug)
            do iplon = 1, myim(jj)
               if (k .le. laytrop(iplon, jj)) then
                  !  --- ...  lower atmosphere loop
                  ind0 = ((jp(iplon, k, jj)-1)*5 + (jt (iplon, k, jj)-1)) * nspa(2) + 1
                  ind1 = ( jp(iplon, k, jj)   *5 + (jt1(iplon, k, jj)-1)) * nspa(2) + 1
                  inds = indself(iplon, k, jj)
                  indf = indfor(iplon, k, jj)

                  ind0p = ind0 + 1
                  ind1p = ind1 + 1
                  indsp = inds + 1
                  indfp = indf + 1

                  corradj = f_one - 0.05 * (pavel(iplon, k, jj) - 100.0) / 900.0
                  !$acc loop seq
                  do ig = 1, ng02
                     ib = ngb(ns02+ig)
                     tauself = selffac(iplon, k, jj) * (selfref(ig,inds) + selffrac(iplon, k, jj)        &
                     &            * (selfref(ig,indsp) - selfref(ig,inds)))
                     taufor  = forfac(iplon, k, jj) * (forref(ig,indf) + forfrac(iplon, k, jj)           &
                     &            * (forref(ig,indfp) - forref(ig,indf))) 

                     taug = corradj * (colamt(iplon, k,1, jj)                      &
                     &            * (fac00(iplon, k, jj)*absa(ig,ind0) + fac10(iplon, k, jj)*absa(ig,ind0p)   &
                     &            +  fac01(iplon, k, jj)*absa(ig,ind1) + fac11(iplon, k, jj)*absa(ig,ind1p))  &
                     &            + tauself + taufor)

                     fracs(iplon, k, ns02+ig, jj) = fracrefa(ig)
                     tautot(iplon, k, ns02+ig, jj) = taug + tauaer(iplon, k, ib, jj)
                  enddo
               else
                  ind0 = ((jp(iplon, k, jj)-13)*5 + (jt (iplon, k, jj)-1)) * nspb(2) + 1
                  ind1 = ((jp(iplon, k, jj)-12)*5 + (jt1(iplon, k, jj)-1)) * nspb(2) + 1
                  indf = indfor(iplon, k, jj)

                  ind0p = ind0 + 1
                  ind1p = ind1 + 1
                  indfp = indf + 1
                  !$acc loop seq
                  do ig = 1, ng02
                     ib = ngb(ns02+ig)
                     taufor = forfac(iplon, k, jj) * (forref(ig,indf) + forfrac(iplon, k, jj)            &
                     &           * (forref(ig,indfp) - forref(ig,indf))) 

                     taug = colamt(iplon, k,1, jj)                                 &
                     &           * (fac00(iplon, k, jj)*absb(ig,ind0) + fac10(iplon, k, jj)*absb(ig,ind0p)    &
                     &           +  fac01(iplon, k, jj)*absb(ig,ind1) + fac11(iplon, k, jj)*absb(ig,ind1p))   &
                     &           + taufor

                     fracs(iplon, k, ns02+ig, jj) = fracrefb(ig)
                     tautot(iplon, k, ns02+ig, jj) = taug + tauaer(iplon, k, ib, jj)
                  enddo
               end if
            end do
         end do
      end do

! ..................................
      end subroutine taugb02
! ----------------------------------

! ----------------------------------
      subroutine taugb03 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!     band 3:  500-630 cm-1 (low key - h2o,co2; low minor - n2o)       !
!                           (high key - h2o,co2; high minor - n2o)     !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb03
      integer, intent(in) :: nlay, laytrop(ix, fulljj), ix, myim(fulljj), &
         async_id, small_ix, i2, fulljj

      integer, dimension(ix, nlay, fulljj), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(ix, nlay, fulljj), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj), intent(in):: colamt
      real (kind=kind_phys), dimension(ix, nlay,maxxsec, fulljj),intent(in):: wx

      real (kind=kind_phys), dimension(ix, nlay, nbands, fulljj), intent(in):: tauaer

      real (kind=kind_phys), dimension(ix, nlay,nrates,2, fulljj), intent(in) ::    &
     &       rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(ix, nlay, ngptlw, fulljj), intent(out) ::     &
     &       fracs, tautot
!  ---  locals:
      integer :: k, ind0, ind1, inds, indsp, indf, indfp, indm, indmp,  &
     &       id000, id010, id100, id110, id200, id210, jmn2o, jmn2op,   &
     &       id001, id011, id101, id111, id201, id211, jpl, jplp,       &
     &       ig, js, js1

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
      integer :: jj, iplon, iplon2, ib
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

     
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !$acc loop vector private(iplon, speccomb, specparm, specmult, js, &
            !$acc&         fs, ind0, speccomb1, specparm1, specmult1, js1, fs1, ind1, &
            !$acc&         speccomb_mn2o, specparm_mn2o, specmult_mn2o, jmn2o, fmn2o, &
            !$acc&         speccomb_planck, specparm_planck, specmult_planck, jpl, &
            !$acc&         fpl, inds, indf, indm, indsp, indfp, indmp, jmn2op, jplp, p, &
            !$acc&         ratn2o, adjfac, adjcoln2o, p4, fk0, fk1, fk2, id000, id010, &
            !$acc&         id100, id110, id200, id210, fac000, fac100, fac200, fac010, &
            !$acc&         fac110, fac210, id001, id011, id101, id111, id201, id211, &
            !$acc&         fac001, fac101, fac201, fac011, fac111, fac211, tauself, &
            !$acc&         taufor, n2om1, n2om2, absn2o, tau_major, tau_major1, taug, &
            !$acc&         colamtr1, colamtr2, colamtr4, k_mn2or0, k_mn2or1)
            do iplon = 1, myim(jj)
               colamtr1 = colamt(iplon, k, 1, jj)
               colamtr2 = colamt(iplon, k, 2, jj)
               colamtr4 = colamt(iplon, k, 4, jj)
               if (k .le. laytrop(iplon, jj)) then
               !  --- ...  lower atmosphere loop
                  speccomb = colamtr1 + rfrate(iplon, k,1,1, jj)*colamtr2
                  specparm = colamtr1 / speccomb
                  specmult = 8.0 * min(specparm, oneminus)
                  js = 1 + int(specmult)
                  fs = mod(specmult, f_one)        
                  ind0 = ((jp(iplon, k, jj)-1)*5 + (jt(iplon, k, jj)-1)) * nspa(3) + js

                  speccomb1 = colamtr1 + rfrate(iplon, k,1,2, jj)*colamtr2
                  specparm1 = colamtr1 / speccomb1
                  specmult1 = 8.0 * min(specparm1, oneminus)
                  js1 = 1 + int(specmult1)
                  fs1 = mod(specmult1, f_one)
                  ind1 = (jp(iplon, k, jj)*5 + (jt1(iplon, k, jj)-1)) * nspa(3) + js1

                  speccomb_mn2o = colamtr1 + refrat_m_a*colamtr2
                  specparm_mn2o = colamtr1 / speccomb_mn2o
                  specmult_mn2o = 8.0 * min(specparm_mn2o, oneminus)
                  jmn2o = 1 + int(specmult_mn2o)
                  fmn2o = mod(specmult_mn2o, f_one)

                  speccomb_planck = colamtr1 + refrat_planck_a*colamtr2
                  specparm_planck = colamtr1 / speccomb_planck
                  specmult_planck = 8.0 * min(specparm_planck, oneminus)
                  jpl = 1 + int(specmult_planck)
                  fpl = mod(specmult_planck, f_one)

                  inds = indself(iplon, k, jj)
                  indf = indfor(iplon, k, jj)
                  indm = indminor(iplon, k, jj)
                  indsp = inds + 1
                  indfp = indf + 1
                  indmp = indm + 1
                  jmn2op= jmn2o+ 1
                  jplp  = jpl  + 1

                  !  --- ...  in atmospheres where the amount of n2o is too great to be considered
                  !           a minor species, adjust the column amount of n2o by an empirical factor
                  !           to obtain the proper contribution.

                  p = coldry(iplon, k, jj) * chi_mls(4,jp(iplon, k, jj)+1)
                  ratn2o = colamtr4 / p
                  if (ratn2o > 1.5) then
                     adjfac = 0.5 + (ratn2o - 0.5)**0.65
                     adjcoln2o = adjfac * p
                  else
                     adjcoln2o = colamtr4
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

                  fac000 = fk0*fac00(iplon, k, jj)
                  fac100 = fk1*fac00(iplon, k, jj)
                  fac200 = fk2*fac00(iplon, k, jj)
                  fac010 = fk0*fac10(iplon, k, jj)
                  fac110 = fk1*fac10(iplon, k, jj)
                  fac210 = fk2*fac10(iplon, k, jj)

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

                  fac001 = fk0*fac01(iplon, k, jj)
                  fac101 = fk1*fac01(iplon, k, jj)
                  fac201 = fk2*fac01(iplon, k, jj)
                  fac011 = fk0*fac11(iplon, k, jj)
                  fac111 = fk1*fac11(iplon, k, jj)
                  fac211 = fk2*fac11(iplon, k, jj)
                  !$acc loop seq
                  do ig = 1, ng03
                     ib = ngb(ns03+ig)
                     k_mn2or0 = ka_mn2o(ig,jmn2o,indm)
                     k_mn2or1 = ka_mn2o(ig,jmn2o,indmp)
                     tauself = selffac(iplon, k, jj)* (selfref(ig,inds) + selffrac(iplon, k, jj)         &
                     &            * (selfref(ig,indsp) - selfref(ig,inds)))
                     taufor  = forfac(iplon, k, jj) * (forref(ig,indf) + forfrac(iplon, k, jj)           &
                     &            * (forref(ig,indfp) - forref(ig,indf)))
                     n2om1   = k_mn2or0 + fmn2o                      &
                     &            * (ka_mn2o(ig,jmn2op,indm) - k_mn2or0)
                     n2om2   = k_mn2or1 + fmn2o                     &
                     &            * (ka_mn2o(ig,jmn2op,indmp) - k_mn2or1)
                     absn2o  = n2om1 + minorfrac(iplon, k, jj) * (n2om2 - n2om1)

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

                     fracs(iplon, k, ns03+ig, jj) = fracrefa(ig,jpl) + fpl                     &
                     &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
                     tautot(iplon, k, ns03+ig, jj) = taug + tauaer(iplon, k, ib, jj)
                  enddo     ! end do_k_loop
               else
                  speccomb = colamtr1 + rfrate(iplon, k,1,1, jj)*colamtr2
                  specparm = colamtr1 / speccomb
                  specmult = 4.0 * min(specparm, oneminus)
                  js = 1 + int(specmult)
                  fs = mod(specmult, f_one)
                  ind0 = ((jp(iplon, k, jj)-13)*5 + (jt(iplon, k, jj)-1)) * nspb(3) + js

                  speccomb1 = colamtr1 + rfrate(iplon, k,1,2, jj)*colamtr2
                  specparm1 = colamtr1 / speccomb1
                  specmult1 = 4.0 * min(specparm1, oneminus)
                  js1 = 1 + int(specmult1)
                  fs1 = mod(specmult1, f_one)
                  ind1 = ((jp(iplon, k, jj)-12)*5 + (jt1(iplon, k, jj)-1)) * nspb(3) + js1

                  speccomb_mn2o = colamtr1 + refrat_m_b*colamtr2
                  specparm_mn2o = colamtr1 / speccomb_mn2o
                  specmult_mn2o = 4.0 * min(specparm_mn2o, oneminus)
                  jmn2o = 1 + int(specmult_mn2o)
                  fmn2o = mod(specmult_mn2o, f_one)

                  speccomb_planck = colamtr1 + refrat_planck_b*colamtr2
                  specparm_planck = colamtr1 / speccomb_planck
                  specmult_planck = 4.0 * min(specparm_planck, oneminus)
                  jpl = 1 + int(specmult_planck)
                  fpl = mod(specmult_planck, f_one)

                  indf = indfor(iplon, k, jj)
                  indm = indminor(iplon, k, jj)
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

                  p = coldry(iplon, k, jj) * chi_mls(4,jp(iplon, k, jj)+1)
                  ratn2o = colamtr4 / p
                  if (ratn2o > 1.5) then
                     adjfac = 0.5 + (ratn2o - 0.5)**0.65
                     adjcoln2o = adjfac * p
                  else
                     adjcoln2o = colamtr4
                  endif

                  fk0 = f_one - fs
                  fk1 = fs
                  fac000 = fk0*fac00(iplon, k, jj)
                  fac010 = fk0*fac10(iplon, k, jj)
                  fac100 = fk1*fac00(iplon, k, jj)
                  fac110 = fk1*fac10(iplon, k, jj)

                  fk0 = f_one - fs1
                  fk1 = fs1
                  fac001 = fk0*fac01(iplon, k, jj)
                  fac011 = fk0*fac11(iplon, k, jj)
                  fac101 = fk1*fac01(iplon, k, jj)
                  fac111 = fk1*fac11(iplon, k, jj)
                  !$acc loop seq
                  do ig = 1, ng03
                     ib = ngb(ns03+ig)
                     k_mn2or0 = kb_mn2o(ig,jmn2o,indm)
                     k_mn2or1 = kb_mn2o(ig,jmn2o,indmp)
                     taufor = forfac(iplon, k, jj) * (forref(ig,indf) + forfrac(iplon, k, jj)            &
                     &           * (forref(ig,indfp) - forref(ig,indf))) 
                     n2om1  = k_mn2or0 + fmn2o                       &
                     &           * (kb_mn2o(ig,jmn2op,indm) - k_mn2or0)
                     n2om2  = k_mn2or1 + fmn2o                      &
                     &           * (kb_mn2o(ig,jmn2op,indmp) - k_mn2or1)
                     absn2o = n2om1 + minorfrac(iplon, k, jj) * (n2om2 - n2om1)

                     tau_major = speccomb                                          &
                     &              * (fac000*absb(ig,id000) + fac010*absb(ig,id010)    &
                     &              +  fac100*absb(ig,id100) + fac110*absb(ig,id110))

                     tau_major1 = speccomb1                                        &
                     &              * (fac001*absb(ig,id001) + fac011*absb(ig,id011)    &
                     &              +  fac101*absb(ig,id101) + fac111*absb(ig,id111))

                     taug = tau_major + tau_major1                      &
                     &                    + taufor + adjcoln2o*absn2o            

                     fracs(iplon, k, ns03+ig, jj) = fracrefb(ig,jpl) + fpl                     &
                     &                     * (fracrefb(ig,jplp) - fracrefb(ig,jpl))
                     tautot(iplon, k, ns03+ig, jj) = taug + tauaer(iplon, k, ib, jj)
                  enddo
               end if
            end do
         end do
      end do

! ..................................
      end subroutine taugb03
! ----------------------------------

! ----------------------------------
      subroutine taugb04 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!     band 4:  630-700 cm-1 (low key - h2o,co2; high key - o3,co2)     !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb04
      integer, intent(in) :: nlay, laytrop(ix, fulljj), ix, myim(fulljj), &
         async_id, small_ix, i2, fulljj

      integer, dimension(ix, nlay, fulljj), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(ix, nlay, fulljj), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj), intent(in):: colamt
      real (kind=kind_phys), dimension(ix, nlay,maxxsec, fulljj),intent(in):: wx

      real (kind=kind_phys), dimension(ix, nlay, nbands, fulljj), intent(in):: tauaer

      real (kind=kind_phys), dimension(ix, nlay,nrates,2, fulljj), intent(in) ::    &
     &       rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(ix, nlay, ngptlw, fulljj), intent(out) ::     &
     &       fracs, tautot
!  ---  locals:
      integer :: k, ind0, ind1, inds, indsp, indf, indfp, jpl, jplp,    &
     &       id000, id010, id100, id110, id200, id210, ig, js, js1,     &
     &       id001, id011, id101, id111, id201, id211

      real (kind=kind_phys) :: tauself, taufor, p, p4, fk0, fk1, fk2,   &
     &      speccomb,       specparm,       specmult,       fs,         &
     &      speccomb1,      specparm1,      specmult1,      fs1,        &
     &      speccomb_planck,specparm_planck,specmult_planck,fpl,        &
     &      fac000, fac100, fac200, fac010, fac110, fac210,             &
     &      fac001, fac101, fac201, fac011, fac111, fac211,             &
     &      refrat_planck_a, refrat_planck_b, tau_major, tau_major1, taug
      integer :: jj, iplon, iplon2, ib
!
!===> ...  begin here
!
      refrat_planck_a = chi_mls(1,11)/chi_mls(2,11)     ! p = 142.5940 mb
      refrat_planck_b = chi_mls(3,13)/chi_mls(2,13)     ! p = 95.58350 mb
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !$acc loop vector private(iplon, speccomb, specparm, specmult, js, &
            !$acc&         fs, ind0, speccomb1, specparm1, specmult1, js1, fs1, ind1, &
            !$acc&         speccomb_planck, specparm_planck, specmult_planck, jpl, &
            !$acc&         fpl, inds, indf, indsp, indfp, jplp, p, &
            !$acc&         p4, fk0, fk1, fk2, id000, id010, &
            !$acc&         id100, id110, id200, id210, fac000, fac100, fac200, fac010, &
            !$acc&         fac110, fac210, id001, id011, id101, id111, id201, id211, &
            !$acc&         fac001, fac101, fac201, fac011, fac111, fac211, tauself, &
            !$acc&         taufor, tau_major, tau_major1, taug)
            do iplon = 1, myim(jj)
               if (k .le. laytrop(iplon, jj)) then
                  !  --- ...  lower atmosphere loop
                  speccomb = colamt(iplon, k,1, jj) + rfrate(iplon, k,1,1, jj)*colamt(iplon, k,2, jj)
                  specparm = colamt(iplon, k,1, jj) / speccomb
                  specmult = 8.0 * min(specparm, oneminus)
                  js = 1 + int(specmult)
                  fs = mod(specmult, f_one)
                  ind0 = ((jp(iplon, k, jj)-1)*5 + (jt(iplon, k, jj)-1)) * nspa(4) + js

                  speccomb1 = colamt(iplon, k,1, jj) + rfrate(iplon, k,1,2, jj)*colamt(iplon, k,2, jj)
                  specparm1 = colamt(iplon, k,1, jj) / speccomb1
                  specmult1 = 8.0 * min(specparm1, oneminus)
                  js1 = 1 + int(specmult1)
                  fs1 = mod(specmult1, f_one)
                  ind1 = ( jp(iplon, k, jj)*5 + (jt1(iplon, k, jj)-1)) * nspa(4) + js1

                  speccomb_planck = colamt(iplon, k,1, jj) + refrat_planck_a*colamt(iplon, k,2, jj)
                  specparm_planck = colamt(iplon, k,1, jj) / speccomb_planck
                  specmult_planck = 8.0 * min(specparm_planck, oneminus)
                  jpl = 1 + int(specmult_planck)
                  fpl = mod(specmult_planck, 1.0)

                  inds = indself(iplon, k, jj)
                  indf = indfor(iplon, k, jj)
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

                  fac000 = fk0*fac00(iplon, k, jj)
                  fac100 = fk1*fac00(iplon, k, jj)
                  fac200 = fk2*fac00(iplon, k, jj)
                  fac010 = fk0*fac10(iplon, k, jj)
                  fac110 = fk1*fac10(iplon, k, jj)
                  fac210 = fk2*fac10(iplon, k, jj)

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

                  fac001 = fk0*fac01(iplon, k, jj)
                  fac101 = fk1*fac01(iplon, k, jj)
                  fac201 = fk2*fac01(iplon, k, jj)
                  fac011 = fk0*fac11(iplon, k, jj)
                  fac111 = fk1*fac11(iplon, k, jj)
                  fac211 = fk2*fac11(iplon, k, jj)
                  !$acc loop seq
                  do ig = 1, ng04
                     ib = ngb(ns04+ig)
                     tauself = selffac(iplon, k, jj)* (selfref(ig,inds) + selffrac(iplon, k, jj)         &
                     &            * (selfref(ig,indsp) - selfref(ig,inds)))
                     taufor  = forfac(iplon, k, jj) * (forref(ig,indf) + forfrac(iplon, k, jj)           &
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

                     fracs(iplon, k, ns04+ig, jj) = fracrefa(ig,jpl) + fpl                     &
                     &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
                     tautot(iplon, k, ns04+ig, jj) = taug + tauaer(iplon, k, ib, jj)
                  enddo     ! end do_k_loop
               else
                  speccomb = colamt(iplon, k,3, jj) + rfrate(iplon, k,6,1, jj)*colamt(iplon, k,2, jj)
                  specparm = colamt(iplon, k,3, jj) / speccomb
                  specmult = 4.0 * min(specparm, oneminus)
                  js = 1 + int(specmult)
                  fs = mod(specmult, f_one)
                  ind0 = ((jp(iplon, k, jj)-13)*5 + (jt(iplon, k, jj)-1)) * nspb(4) + js

                  speccomb1 = colamt(iplon, k,3, jj) + rfrate(iplon, k,6,2, jj)*colamt(iplon, k,2, jj)
                  specparm1 = colamt(iplon, k,3, jj) / speccomb1
                  specmult1 = 4.0 * min(specparm1, oneminus)
                  js1 = 1 + int(specmult1)
                  fs1 = mod(specmult1, f_one)
                  ind1 = ((jp(iplon, k, jj)-12)*5 + (jt1(iplon, k, jj)-1)) * nspb(4) + js1

                  speccomb_planck = colamt(iplon, k,3, jj) + refrat_planck_b*colamt(iplon, k,2, jj)
                  specparm_planck = colamt(iplon, k,3, jj) / speccomb_planck
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
                  fac000 = fk0*fac00(iplon, k, jj)
                  fac010 = fk0*fac10(iplon, k, jj)
                  fac100 = fk1*fac00(iplon, k, jj)
                  fac110 = fk1*fac10(iplon, k, jj)

                  fk0 = f_one - fs1
                  fk1 = fs1
                  fac001 = fk0*fac01(iplon, k, jj)
                  fac011 = fk0*fac11(iplon, k, jj)
                  fac101 = fk1*fac01(iplon, k, jj)
                  fac111 = fk1*fac11(iplon, k, jj)
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

                     fracs(iplon, k, ns04+ig, jj) = fracrefb(ig,jpl) + fpl                     &
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
                     tautot(iplon, k, ns04+ig, jj) = taug + tauaer(iplon, k, ib, jj)
                  enddo
               end if
            end do
         end do
      end do

! ..................................
      end subroutine taugb04
! ----------------------------------

! ----------------------------------
      subroutine taugb05 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!     band 5:  700-820 cm-1 (low key - h2o,co2; low minor - o3, ccl4)  !
!                           (high key - o3,co2)                        !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb05
      integer, intent(in) :: nlay, laytrop(ix, fulljj), ix, myim(fulljj), &
         async_id, small_ix, i2, fulljj

      integer, dimension(ix, nlay, fulljj), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(ix, nlay, fulljj), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj), intent(in):: colamt
      real (kind=kind_phys), dimension(ix, nlay,maxxsec, fulljj),intent(in):: wx

      real (kind=kind_phys), dimension(ix, nlay, nbands, fulljj), intent(in):: tauaer

      real (kind=kind_phys), dimension(ix, nlay,nrates,2, fulljj), intent(in) ::    &
     &       rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(ix, nlay, ngptlw, fulljj), intent(out) ::     &
     &       fracs, tautot
!  ---  locals: 
      integer :: k, ind0, ind1, inds, indsp, indf, indfp, indm, indmp,  &
     &       id000, id010, id100, id110, id200, id210, jmo3, jmo3p,     &
     &       id001, id011, id101, id111, id201, id211, jpl, jplp,       &
     &       ig, js, js1

      real (kind=kind_phys)  :: tauself, taufor, o3m1, o3m2, abso3,     &
     &      speccomb,       specparm,       specmult,       fs,         &
     &      speccomb1,      specparm1,      specmult1,      fs1,        &
     &      speccomb_mo3,   specparm_mo3,   specmult_mo3,   fmo3,       &
     &      speccomb_planck,specparm_planck,specmult_planck,fpl,        &
     &      refrat_planck_a, refrat_planck_b, refrat_m_a,               &
     &      fac000, fac100, fac200, fac010, fac110, fac210,             &
     &      fac001, fac101, fac201, fac011, fac111, fac211,             &
     &      p0, p40, fk00, fk10, fk20, p1, p41, fk01, fk11, fk21, taug
      integer :: jj, iplon, iplon2, ib
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

      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !$acc loop vector private(iplon, speccomb, specparm, specmult, js, &
            !$acc&         fs, ind0, speccomb1, specparm1, specmult1, js1, fs1, ind1, &
            !$acc&         speccomb_planck, specparm_planck, specmult_planck, jpl, &
            !$acc&         fpl, inds, indf, indm, indsp, indfp, indmp, jplp, &
            !$acc&         id000, id010, speccomb_mo3, specparm_mo3, specmult_mo3, &
            !$acc&         id100, id110, id200, id210, fac000, fac100, fac200, fac010, &
            !$acc&         fac110, fac210, id001, id011, id101, id111, id201, id211, &
            !$acc&         fac001, fac101, fac201, fac011, fac111, fac211, tauself, &
            !$acc&         taufor, jmo3, fmo3, jmo3p, p0, p40, fk00, fk10, fk20, p1, &
            !$acc&         p41, fk01, fk11, fk21, o3m1, o3m2, abso3, taug)
            do iplon = 1, myim(jj)
               if (k .le. laytrop(iplon, jj)) then
                  !  --- ...  lower atmosphere loop
                  speccomb = colamt(iplon, k,1, jj) + rfrate(iplon, k,1,1, jj)*colamt(iplon, k,2, jj)
                  specparm = colamt(iplon, k,1, jj) / speccomb
                  specmult = 8.0 * min(specparm, oneminus)
                  js = 1 + int(specmult)
                  fs = mod(specmult, f_one)
                  ind0 = ((jp(iplon, k, jj)-1)*5 + (jt(iplon, k, jj)-1)) * nspa(5) + js

                  speccomb1 = colamt(iplon, k,1, jj) + rfrate(iplon, k,1,2, jj)*colamt(iplon, k,2, jj)
                  specparm1 = colamt(iplon, k,1, jj) / speccomb1
                  specmult1 = 8.0 * min(specparm1, oneminus)
                  js1 = 1 + int(specmult1)
                  fs1 = mod(specmult1, f_one)
                  ind1 = (jp(iplon, k, jj)*5 + (jt1(iplon, k, jj)-1)) * nspa(5) + js1

                  speccomb_mo3 = colamt(iplon, k,1, jj) + refrat_m_a*colamt(iplon, k,2, jj)
                  specparm_mo3 = colamt(iplon, k,1, jj) / speccomb_mo3
                  specmult_mo3 = 8.0 * min(specparm_mo3, oneminus)
                  jmo3 = 1 + int(specmult_mo3)
                  fmo3 = mod(specmult_mo3, f_one)

                  speccomb_planck = colamt(iplon, k,1, jj) + refrat_planck_a*colamt(iplon, k,2, jj)
                  specparm_planck = colamt(iplon, k,1, jj) / speccomb_planck
                  specmult_planck = 8.0 * min(specparm_planck, oneminus)
                  jpl = 1 + int(specmult_planck)
                  fpl = mod(specmult_planck, f_one)

                  inds = indself(iplon, k, jj)
                  indf = indfor(iplon, k, jj)
                  indm = indminor(iplon, k, jj)
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

                  fac000 = fk00 * fac00(iplon, k, jj)
                  fac100 = fk10 * fac00(iplon, k, jj)
                  fac200 = fk20 * fac00(iplon, k, jj)
                  fac010 = fk00 * fac10(iplon, k, jj)
                  fac110 = fk10 * fac10(iplon, k, jj)
                  fac210 = fk20 * fac10(iplon, k, jj)

                  fac001 = fk01 * fac01(iplon, k, jj)
                  fac101 = fk11 * fac01(iplon, k, jj)
                  fac201 = fk21 * fac01(iplon, k, jj)
                  fac011 = fk01 * fac11(iplon, k, jj)
                  fac111 = fk11 * fac11(iplon, k, jj)
                  fac211 = fk21 * fac11(iplon, k, jj)
                  !$acc loop seq
                  do ig = 1, ng05
                     ib = ngb(ns05+ig)
                     tauself = selffac(iplon, k, jj) * (selfref(ig,inds) + selffrac(iplon, k, jj)        &
                     &            * (selfref(ig,indsp) - selfref(ig,inds)))
                     taufor  = forfac(iplon, k, jj) * (forref(ig,indf) + forfrac(iplon, k, jj)           &
                     &            * (forref(ig,indfp) - forref(ig,indf)))
                     o3m1    = ka_mo3(ig,jmo3,indm) + fmo3                         &
                     &            * (ka_mo3(ig,jmo3p,indm) -  ka_mo3(ig,jmo3,indm))
                     o3m2    = ka_mo3(ig,jmo3,indmp) + fmo3                        &
                     &            * (ka_mo3(ig,jmo3p,indmp) - ka_mo3(ig,jmo3,indmp))
                     abso3   = o3m1 + minorfrac(iplon, k, jj)*(o3m2 - o3m1)

                     taug = speccomb                                    &
                     &            * (fac000*absa(ig,id000) + fac010*absa(ig,id010)      &
                     &            +  fac100*absa(ig,id100) + fac110*absa(ig,id110)      &
                     &            +  fac200*absa(ig,id200) + fac210*absa(ig,id210))     &
                     &            +     speccomb1                                       &
                     &            * (fac001*absa(ig,id001) + fac011*absa(ig,id011)      &
                     &            +  fac101*absa(ig,id101) + fac111*absa(ig,id111)      &
                     &            +  fac201*absa(ig,id201) + fac211*absa(ig,id211))     &
                     &            + tauself + taufor+abso3*colamt(iplon, k,3, jj)+wx(iplon, k,1, jj)*ccl4(ig)

                     fracs(iplon, k, ns05+ig, jj) = fracrefa(ig,jpl) + fpl                     &
                     &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
                     tautot(iplon, k, ns05+ig, jj) = taug + tauaer(iplon, k, ib, jj)
                  enddo
               else
                  speccomb = colamt(iplon, k,3, jj) + rfrate(iplon, k,6,1, jj)*colamt(iplon, k,2, jj)
                  specparm = colamt(iplon, k,3, jj) / speccomb
                  specmult = 4.0 * min(specparm, oneminus)
                  js = 1 + int(specmult)
                  fs = mod(specmult, f_one)
                  ind0 = ((jp(iplon, k, jj)-13)*5 + (jt(iplon, k, jj)-1)) * nspb(5) + js

                  speccomb1 = colamt(iplon, k,3, jj) + rfrate(iplon, k,6,2, jj)*colamt(iplon, k,2, jj)
                  specparm1 = colamt(iplon, k,3, jj) / speccomb1
                  specmult1 = 4.0 * min(specparm1, oneminus)
                  js1 = 1 + int(specmult1)
                  fs1 = mod(specmult1, f_one)
                  ind1 = ((jp(iplon, k, jj)-12)*5 + (jt1(iplon, k, jj)-1)) * nspb(5) + js1

                  speccomb_planck = colamt(iplon, k,3, jj) + refrat_planck_b*colamt(iplon, k,2, jj)
                  specparm_planck = colamt(iplon, k,3, jj) / speccomb_planck
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

                  fac000 = fk00 * fac00(iplon, k, jj)
                  fac010 = fk00 * fac10(iplon, k, jj)
                  fac100 = fk10 * fac00(iplon, k, jj)
                  fac110 = fk10 * fac10(iplon, k, jj)

                  fac001 = fk01 * fac01(iplon, k, jj)
                  fac011 = fk01 * fac11(iplon, k, jj)
                  fac101 = fk11 * fac01(iplon, k, jj)
                  fac111 = fk11 * fac11(iplon, k, jj)
                  !$acc loop seq
                  do ig = 1, ng05
                     ib = ngb(ns05+ig)
                     taug = speccomb                                    &
                     &                * (fac000*absb(ig,id000) + fac010*absb(ig,id010)  &
                     &                +  fac100*absb(ig,id100) + fac110*absb(ig,id110)) &
                     &                +     speccomb1                                   &
                     &                * (fac001*absb(ig,id001) + fac011*absb(ig,id011)  &
                     &                +  fac101*absb(ig,id101) + fac111*absb(ig,id111)) &
                     &                + wx(iplon, k,1, jj) * ccl4(ig)

                     fracs(iplon, k, ns05+ig, jj) = fracrefb(ig,jpl) + fpl                     &
                     &                     * (fracrefb(ig,jplp) - fracrefb(ig,jpl))
                     tautot(iplon, k, ns05+ig, jj) = taug + tauaer(iplon, k, ib, jj)
                  enddo
               end if
            end do
         end do
      end do

! ..................................
      end subroutine taugb05
! ----------------------------------

! ----------------------------------
      subroutine taugb06 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!     band 6:  820-980 cm-1 (low key - h2o; low minor - co2)           !
!                           (high key - none; high minor - cfc11, cfc12)
!  ------------------------------------------------------------------  !

      use module_radlw_kgb06
      integer, intent(in) :: nlay, laytrop(ix, fulljj), ix, myim(fulljj), &
         async_id, small_ix, i2, fulljj

      integer, dimension(ix, nlay, fulljj), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(ix, nlay, fulljj), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj), intent(in):: colamt
      real (kind=kind_phys), dimension(ix, nlay,maxxsec, fulljj),intent(in):: wx

      real (kind=kind_phys), dimension(ix, nlay, nbands, fulljj), intent(in):: tauaer

      real (kind=kind_phys), dimension(ix, nlay,nrates,2, fulljj), intent(in) ::    &
     &       rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(ix, nlay, ngptlw, fulljj), intent(out) ::     &
     &       fracs, tautot
!  ---  locals: 
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       indm, indmp, ig

      real (kind=kind_phys) :: ratco2, adjfac, adjcolco2, tauself,      &
     &      taufor, absco2, temp, taug
      integer :: jj, iplon, iplon2, ib
!
!===> ...  begin here
!
!  --- ...  minor gas mapping level:
!     lower - co2, p = 706.2720 mb, t = 294.2 k
!     upper - cfc11, cfc12

      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !$acc loop vector private(iplon, ind0, ind1, inds, indf, &
            !$acc&         indm, indsp, indfp, indmp, ind0p, ind1p, temp,  ratco2, adjfac, &
            !$acc&         adjcolco2, tauself, taufor, absco2, taug)
            do iplon = 1, myim(jj)
               if (k .le. laytrop(iplon, jj)) then
                  !  --- ...  lower atmosphere loop
                  ind0 = ((jp(iplon, k, jj)-1)*5 + (jt (iplon, k, jj)-1)) * nspa(6) + 1
                  ind1 = ( jp(iplon, k, jj)   *5 + (jt1(iplon, k, jj)-1)) * nspa(6) + 1

                  inds = indself(iplon, k, jj)
                  indf = indfor(iplon, k, jj)
                  indm = indminor(iplon, k, jj)
                  indsp = inds + 1
                  indfp = indf + 1
                  indmp = indm + 1
                  ind0p = ind0 + 1
                  ind1p = ind1 + 1

                  !  --- ...  in atmospheres where the amount of co2 is too great to be considered
                  !           a minor species, adjust the column amount of co2 by an empirical factor
                  !           to obtain the proper contribution.

                  temp   = coldry(iplon, k, jj) * chi_mls(2,jp(iplon, k, jj)+1)
                  ratco2 = colamt(iplon, k,2, jj) / temp
                  if (ratco2 > 3.0) then
                     adjfac = 2.0 + (ratco2-2.0)**0.77
                     adjcolco2 = adjfac * temp
                  else
                     adjcolco2 = colamt(iplon, k,2, jj)
                  endif
                  !$acc loop seq
                  do ig = 1, ng06
                     ib = ngb(ns06+ig)
                     tauself = selffac(iplon, k, jj) * (selfref(ig,inds) + selffrac(iplon, k, jj)        &
                     &            * (selfref(ig,indsp) - selfref(ig,inds)))
                     taufor  = forfac(iplon, k, jj) * (forref(ig,indf) + forfrac(iplon, k, jj)           &
                     &            * (forref(ig,indfp) - forref(ig,indf)))
                     absco2  = ka_mco2(ig,indm) + minorfrac(iplon, k, jj)                     &
                     &            * (ka_mco2(ig,indmp) - ka_mco2(ig,indm))

                     taug = colamt(iplon, k,1, jj)                                 &
                     &            * (fac00(iplon, k, jj)*absa(ig,ind0) + fac10(iplon, k, jj)*absa(ig,ind0p)   &
                     &            +  fac01(iplon, k, jj)*absa(ig,ind1) + fac11(iplon, k, jj)*absa(ig,ind1p))  &
                     &            +  tauself + taufor + adjcolco2*absco2                &
                     &            +  wx(iplon, k,2, jj)*cfc11adj(ig) + wx(iplon, k,3, jj)*cfc12(ig)

                     fracs(iplon, k, ns06+ig, jj) = fracrefa(ig)
                     tautot(iplon, k, ns06+ig, jj) = taug + tauaer(iplon, k, ib, jj)
                  enddo
               else
                  !$acc loop seq
                  do ig = 1, ng06
                     ib = ngb(ns06+ig)
                     taug = wx(iplon, k,2, jj)*cfc11adj(ig) + wx(iplon, k,3, jj)*cfc12(ig)

                     fracs(iplon, k, ns06+ig, jj) = fracrefa(ig)
                     tautot(iplon, k, ns06+ig, jj) = taug + tauaer(iplon, k, ib, jj)
                  enddo
               end if
            end do
         end do
      end do

! ..................................
      end subroutine taugb06
! ----------------------------------

! ----------------------------------
      subroutine taugb07 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!     band 7:  980-1080 cm-1 (low key - h2o,o3; low minor - co2)       !
!                            (high key - o3; high minor - co2)         !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb07
      integer, intent(in) :: nlay, laytrop(ix, fulljj), ix, myim(fulljj), &
         async_id, small_ix, i2, fulljj

      integer, dimension(ix, nlay, fulljj), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(ix, nlay, fulljj), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj), intent(in):: colamt
      real (kind=kind_phys), dimension(ix, nlay,maxxsec, fulljj),intent(in):: wx

      real (kind=kind_phys), dimension(ix, nlay, nbands, fulljj), intent(in):: tauaer

      real (kind=kind_phys), dimension(ix, nlay,nrates,2, fulljj), intent(in) ::    &
     &       rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(ix, nlay, ngptlw, fulljj), intent(out) ::     &
     &       fracs, tautot
!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       id000, id010, id100, id110, id200, id210, indm, indmp,     &
     &       id001, id011, id101, id111, id201, id211, jmco2, jmco2p,   &
     &       jpl, jplp, ig, js, js1

      real (kind=kind_phys) :: tauself, taufor, co2m1, co2m2, absco2,   &
     &      speccomb,       specparm,       specmult,       fs,         &
     &      speccomb1,      specparm1,      specmult1,      fs1,        &
     &      speccomb_mco2,  specparm_mco2,  specmult_mco2,  fmco2,      &
     &      speccomb_planck,specparm_planck,specmult_planck,fpl,        &
     &      refrat_planck_a, refrat_m_a, ratco2, adjfac, adjcolco2,     &
     &      fac000, fac100, fac200, fac010, fac110, fac210,             &
     &      fac001, fac101, fac201, fac011, fac111, fac211,             &
     &      p0, p40, fk00, fk10, fk20, p1, p41, fk01, fk11, fk21, temp, taug
      integer :: jj, iplon, iplon2, ib
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

      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !$acc loop vector private(iplon, speccomb, specparm, specmult, js, &
            !$acc&         fs, ind0, speccomb1, specparm1, specmult1, js1, fs1, ind1, &
            !$acc&         speccomb_planck, specparm_planck, specmult_planck, jpl, &
            !$acc&         fpl, inds, indf, indm, indsp, indfp, indmp, jplp, &
            !$acc&         id000, id010, speccomb_mco2, specparm_mco2, specmult_mco2, &
            !$acc&         id100, id110, id200, id210, fac000, fac100, fac200, fac010, &
            !$acc&         fac110, fac210, id001, id011, id101, id111, id201, id211, &
            !$acc&         fac001, fac101, fac201, fac011, fac111, fac211, tauself, &
            !$acc&         taufor, p0, p40, fk00, fk10, fk20, p1, jmco2, fmco2, jmco2p, &
            !$acc&         p41, fk01, fk11, fk21, ind0p, ind1p, temp, ratco2, adjfac, &
            !$acc&         adjcolco2, co2m1, co2m2, absco2, taug)
            do iplon = 1, myim(jj)
               if (k .le. laytrop(iplon, jj)) then
                  !  --- ...  lower atmosphere loop
                  speccomb = colamt(iplon, k,1, jj) + rfrate(iplon, k,2,1, jj)*colamt(iplon, k,3, jj)
                  specparm = colamt(iplon, k,1, jj) / speccomb
                  specmult = 8.0 * min(specparm, oneminus)
                  js = 1 + int(specmult)
                  fs = mod(specmult, f_one)
                  ind0 = ((jp(iplon, k, jj)-1)*5 + (jt(iplon, k, jj)-1)) * nspa(7) + js

                  speccomb1 = colamt(iplon, k,1, jj) + rfrate(iplon, k,2,2, jj)*colamt(iplon, k,3, jj)
                  specparm1 = colamt(iplon, k,1, jj) / speccomb1
                  specmult1 = 8.0 * min(specparm1, oneminus)
                  js1 = 1 + int(specmult1)
                  fs1 = mod(specmult1, f_one)
                  ind1 = (jp(iplon, k, jj)*5 + (jt1(iplon, k, jj)-1)) * nspa(7) + js1

                  speccomb_mco2 = colamt(iplon, k,1, jj) + refrat_m_a*colamt(iplon, k,3, jj)
                  specparm_mco2 = colamt(iplon, k,1, jj) / speccomb_mco2
                  specmult_mco2 = 8.0 * min(specparm_mco2, oneminus)
                  jmco2 = 1 + int(specmult_mco2)
                  fmco2 = mod(specmult_mco2, f_one)

                  speccomb_planck = colamt(iplon, k,1, jj) + refrat_planck_a*colamt(iplon, k,3, jj)
                  specparm_planck = colamt(iplon, k,1, jj) / speccomb_planck
                  specmult_planck = 8.0 * min(specparm_planck, oneminus)
                  jpl = 1 + int(specmult_planck)
                  fpl = mod(specmult_planck, f_one)

                  inds = indself(iplon, k, jj)
                  indf = indfor(iplon, k, jj)
                  indm = indminor(iplon, k, jj)
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

                  temp   = coldry(iplon, k, jj) * chi_mls(2,jp(iplon, k, jj)+1)
                  ratco2 = colamt(iplon, k,2, jj) / temp
                  if (ratco2 > 3.0) then
                     adjfac = 3.0 + (ratco2-3.0)**0.79
                     adjcolco2 = adjfac * temp
                  else
                     adjcolco2 = colamt(iplon, k,2, jj)
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

                  fac000 = fk00 * fac00(iplon, k, jj)
                  fac100 = fk10 * fac00(iplon, k, jj)
                  fac200 = fk20 * fac00(iplon, k, jj)
                  fac010 = fk00 * fac10(iplon, k, jj)
                  fac110 = fk10 * fac10(iplon, k, jj)
                  fac210 = fk20 * fac10(iplon, k, jj)

                  fac001 = fk01 * fac01(iplon, k, jj)
                  fac101 = fk11 * fac01(iplon, k, jj)
                  fac201 = fk21 * fac01(iplon, k, jj)
                  fac011 = fk01 * fac11(iplon, k, jj)
                  fac111 = fk11 * fac11(iplon, k, jj)
                  fac211 = fk21 * fac11(iplon, k, jj)
                  !$acc loop seq
                  do ig = 1, ng07
                     ib = ngb(ns07+ig)
                     tauself = selffac(iplon, k, jj)* (selfref(ig,inds) + selffrac(iplon, k, jj)         &
                     &            * (selfref(ig,indsp) - selfref(ig,inds)))
                     taufor  = forfac(iplon, k, jj) * (forref(ig,indf) + forfrac(iplon, k, jj)           &
                     &            * (forref(ig,indfp) - forref(ig,indf))) 
                     co2m1   = ka_mco2(ig,jmco2,indm) + fmco2                      &
                     &            * (ka_mco2(ig,jmco2p,indm) - ka_mco2(ig,jmco2,indm))
                     co2m2   = ka_mco2(ig,jmco2,indmp) + fmco2                     &
                     &            * (ka_mco2(ig,jmco2p,indmp) - ka_mco2(ig,jmco2,indmp))
                     absco2  = co2m1 + minorfrac(iplon, k, jj) * (co2m2 - co2m1)

                     taug = speccomb                                    &
                     &                * (fac000*absa(ig,id000) + fac010*absa(ig,id010)  &
                     &                +  fac100*absa(ig,id100) + fac110*absa(ig,id110)  &
                     &                +  fac200*absa(ig,id200) + fac210*absa(ig,id210)) &
                     &                +     speccomb1                                   &
                     &                * (fac001*absa(ig,id001) + fac011*absa(ig,id011)  &
                     &                +  fac101*absa(ig,id101) + fac111*absa(ig,id111)  &
                     &                +  fac201*absa(ig,id201) + fac211*absa(ig,id211)) &
                     &                + tauself + taufor + adjcolco2*absco2

                     fracs(iplon, k, ns07+ig, jj) = fracrefa(ig,jpl) + fpl                     &
                     &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
                     tautot(iplon, k, ns07+ig, jj) = taug + tauaer(iplon, k, ib, jj)
                  enddo
               else
                  temp   = coldry(iplon, k, jj) * chi_mls(2,jp(iplon, k, jj)+1)
                  ratco2 = colamt(iplon, k,2, jj) / temp
                  if (ratco2 > 3.0) then
                     adjfac = 2.0 + (ratco2-2.0)**0.79
                     adjcolco2 = adjfac * temp
                  else
                     adjcolco2 = colamt(iplon, k,2, jj)
                  endif

                  ind0 = ((jp(iplon, k, jj)-13)*5 + (jt (iplon, k, jj)-1)) * nspb(7) + 1
                  ind1 = ((jp(iplon, k, jj)-12)*5 + (jt1(iplon, k, jj)-1)) * nspb(7) + 1

                  indm = indminor(iplon, k, jj)
                  indmp = indm + 1
                  ind0p = ind0 + 1
                  ind1p = ind1 + 1
                  !$acc loop seq
                  do ig = 1, ng07
                     ib = ngb(ns07+ig)
                     absco2 = kb_mco2(ig,indm) + minorfrac(iplon, k, jj)                      &
                     &           * (kb_mco2(ig,indmp) - kb_mco2(ig,indm))

                     taug= colamt(iplon, k,3, jj)                                 &
                     &            * (fac00(iplon, k, jj)*absb(ig,ind0) + fac10(iplon, k, jj)*absb(ig,ind0p)   &
                     &            +  fac01(iplon, k, jj)*absb(ig,ind1) + fac11(iplon, k, jj)*absb(ig,ind1p))  &
                     &            + adjcolco2 * absco2
                  !  --- ...  empirical modification to code to improve stratospheric cooling rates
                  !           for o3.  revised to apply weighting for g-point reduction in this band.
                     if (ig .eq. 6) taug = taug * 0.92
                     if (ig .eq. 7) taug = taug * 0.88
                     if (ig .eq. 8) taug = taug * 1.07
                     if (ig .eq. 9) taug = taug * 1.1
                     if (ig .eq. 10) taug = taug * 0.99
                     if (ig .eq. 11) taug = taug * 0.855
                     fracs(iplon, k, ns07+ig, jj) = fracrefb(ig)
                     tautot(iplon, k, ns07+ig, jj) = taug + tauaer(iplon, k, ib, jj)
                  enddo
               end if
            end do
         end do
      end do

! ..................................
      end subroutine taugb07
! ----------------------------------

! ----------------------------------
      subroutine taugb08 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!     band 8:  1080-1180 cm-1 (low key - h2o; low minor - co2,o3,n2o)  !
!                             (high key - o3; high minor - co2, n2o)   !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb08
      integer, intent(in) :: nlay, laytrop(ix, fulljj), ix, myim(fulljj), &
         async_id, small_ix, i2, fulljj

      integer, dimension(ix, nlay, fulljj), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(ix, nlay, fulljj), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj), intent(in):: colamt
      real (kind=kind_phys), dimension(ix, nlay,maxxsec, fulljj),intent(in):: wx

      real (kind=kind_phys), dimension(ix, nlay, nbands, fulljj), intent(in):: tauaer

      real (kind=kind_phys), dimension(ix, nlay,nrates,2, fulljj), intent(in) ::    &
     &       rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(ix, nlay, ngptlw, fulljj), intent(out) ::     &
     &       fracs, tautot
!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       indm, indmp, ig

      real (kind=kind_phys) :: tauself, taufor, absco2, abso3, absn2o,  &
     &      ratco2, adjfac, adjcolco2, temp, taug
      integer :: jj, iplon, iplon2, ib
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

      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !$acc loop vector private(iplon, ind0, ind1, inds, indf, indm, &
            !$acc&         ind0p, ind1p, indsp, indfp, indmp, temp, ratco2, adjfac, adjcolco2, &
            !$acc&         tauself, taufor, absco2, abso3, absn2o, taug)
            do iplon = 1, myim(jj)
               if (k .le. laytrop(iplon, jj)) then
                  !  --- ...  lower atmosphere loop
                  ind0 = ((jp(iplon, k, jj)-1)*5 + (jt (iplon, k, jj)-1)) * nspa(8) + 1
                  ind1 = ( jp(iplon, k, jj)   *5 + (jt1(iplon, k, jj)-1)) * nspa(8) + 1

                  inds = indself(iplon, k, jj)
                  indf = indfor(iplon, k, jj)
                  indm = indminor(iplon, k, jj)
                  ind0p = ind0 + 1
                  ind1p = ind1 + 1
                  indsp = inds + 1
                  indfp = indf + 1
                  indmp = indm + 1

                  !  --- ...  in atmospheres where the amount of co2 is too great to be considered
                  !           a minor species, adjust the column amount of co2 by an empirical factor
                  !           to obtain the proper contribution.

                  temp   = coldry(iplon, k, jj) * chi_mls(2,jp(iplon, k, jj)+1)
                  ratco2 = colamt(iplon, k,2, jj) / temp
                  if (ratco2 > 3.0) then
                     adjfac = 2.0 + (ratco2-2.0)**0.65
                     adjcolco2 = adjfac * temp
                  else
                     adjcolco2 = colamt(iplon, k,2, jj)
                  endif
                  !$acc loop seq
                  do ig = 1, ng08
                     ib = ngb(ns08+ig)
                     tauself = selffac(iplon, k, jj) * (selfref(ig,inds) + selffrac(iplon, k, jj)        &
                     &            * (selfref(ig,indsp) - selfref(ig,inds)))
                     taufor  = forfac(iplon, k, jj) * (forref(ig,indf) + forfrac(iplon, k, jj)           &
                     &            * (forref(ig,indfp) - forref(ig,indf)))
                     absco2  = (ka_mco2(ig,indm) + minorfrac(iplon, k, jj)                    &
                     &            * (ka_mco2(ig,indmp) - ka_mco2(ig,indm)))
                     abso3   = (ka_mo3(ig,indm) + minorfrac(iplon, k, jj)                     &
                     &            * (ka_mo3(ig,indmp) - ka_mo3(ig,indm)))
                     absn2o  = (ka_mn2o(ig,indm) + minorfrac(iplon, k, jj)                    &
                     &            * (ka_mn2o(ig,indmp) - ka_mn2o(ig,indm)))

                     taug = colamt(iplon, k,1, jj)                                 &
                     &            * (fac00(iplon, k, jj)*absa(ig,ind0) + fac10(iplon, k, jj)*absa(ig,ind0p)   &
                     &            +  fac01(iplon, k, jj)*absa(ig,ind1) + fac11(iplon, k, jj)*absa(ig,ind1p))  &
                     &            + tauself+taufor + adjcolco2*absco2                   &
                     &            + colamt(iplon, k,3, jj)*abso3 + colamt(iplon, k,4, jj)*absn2o              &
                     &            + wx(iplon, k,3, jj)*cfc12(ig) + wx(iplon, k,4, jj)*cfc22adj(ig)

                     fracs(iplon, k, ns08+ig, jj) = fracrefa(ig)
                     tautot(iplon, k, ns08+ig, jj) = taug + tauaer(iplon, k, ib, jj)
                  enddo
               else
                  ind0 = ((jp(iplon, k, jj)-13)*5 + (jt (iplon, k, jj)-1)) * nspb(8) + 1
                  ind1 = ((jp(iplon, k, jj)-12)*5 + (jt1(iplon, k, jj)-1)) * nspb(8) + 1

                  indm = indminor(iplon, k, jj)
                  ind0p = ind0 + 1
                  ind1p = ind1 + 1
                  indmp = indm + 1

                  !  --- ...  in atmospheres where the amount of co2 is too great to be considered
                  !           a minor species, adjust the column amount of co2 by an empirical factor
                  !           to obtain the proper contribution.

                  temp   = coldry(iplon, k, jj) * chi_mls(2,jp(iplon, k, jj)+1)
                  ratco2 = colamt(iplon, k,2, jj) / temp
                  if (ratco2 > 3.0) then
                     adjfac = 2.0 + (ratco2-2.0)**0.65
                     adjcolco2 = adjfac * temp
                  else
                     adjcolco2 = colamt(iplon, k,2, jj)
                  endif
                  !$acc loop seq
                  do ig = 1, ng08
                     ib = ngb(ns08+ig)
                     absco2 = (kb_mco2(ig,indm) + minorfrac(iplon, k, jj)                     &
                     &           * (kb_mco2(ig,indmp) - kb_mco2(ig,indm)))
                     absn2o = (kb_mn2o(ig,indm) + minorfrac(iplon, k, jj)                     &
                     &           * (kb_mn2o(ig,indmp) - kb_mn2o(ig,indm)))

                     taug = colamt(iplon, k,3, jj)                                 &
                     &           * (fac00(iplon, k, jj)*absb(ig,ind0) + fac10(iplon, k, jj)*absb(ig,ind0p)    &
                     &           +  fac01(iplon, k, jj)*absb(ig,ind1) + fac11(iplon, k, jj)*absb(ig,ind1p))   &
                     &           + adjcolco2*absco2 + colamt(iplon, k,4, jj)*absn2o                &
                     &           + wx(iplon, k,3, jj)*cfc12(ig) + wx(iplon, k,4, jj)*cfc22adj(ig)

                     fracs(iplon, k, ns08+ig, jj) = fracrefb(ig)
                     tautot(iplon, k, ns08+ig, jj) = taug + tauaer(iplon, k, ib, jj)
                  enddo
               end if
            end do
         end do
      end do

! ..................................
      end subroutine taugb08
! ----------------------------------

! ----------------------------------
      subroutine taugb09 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!     band 9:  1180-1390 cm-1 (low key - h2o,ch4; low minor - n2o)     !
!                             (high key - ch4; high minor - n2o)       !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb09
      integer, intent(in) :: nlay, laytrop(ix, fulljj), ix, myim(fulljj), &
         async_id, small_ix, i2, fulljj

      integer, dimension(ix, nlay, fulljj), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(ix, nlay, fulljj), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj), intent(in):: colamt
      real (kind=kind_phys), dimension(ix, nlay,maxxsec, fulljj),intent(in):: wx

      real (kind=kind_phys), dimension(ix, nlay, nbands, fulljj), intent(in):: tauaer

      real (kind=kind_phys), dimension(ix, nlay,nrates,2, fulljj), intent(in) ::    &
     &       rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(ix, nlay, ngptlw, fulljj), intent(out) ::     &
     &       fracs, tautot
!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       id000, id010, id100, id110, id200, id210, indm, indmp,     &
     &       id001, id011, id101, id111, id201, id211, jmn2o, jmn2op,   &
     &       jpl, jplp, ig, js, js1

      real (kind=kind_phys) :: tauself, taufor, n2om1, n2om2, absn2o,   &
     &       speccomb,       specparm,       specmult,       fs,        &
     &       speccomb1,      specparm1,      specmult1,      fs1,       &
     &       speccomb_mn2o,  specparm_mn2o,  specmult_mn2o,  fmn2o,     &
     &       speccomb_planck,specparm_planck,specmult_planck,fpl,       &
     &       refrat_planck_a, refrat_m_a, ratn2o, adjfac, adjcoln2o,    &
     &       fac000, fac100, fac200, fac010, fac110, fac210,            &
     &       fac001, fac101, fac201, fac011, fac111, fac211,            &
     &       p0, p40, fk00, fk10, fk20, p1, p41, fk01, fk11, fk21, temp, taug
      integer :: jj, iplon, iplon2, ib
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
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !$acc loop vector private(iplon, speccomb, specparm, specmult, js, &
            !$acc&         fs, ind0, speccomb1, specparm1, specmult1, js1, fs1, ind1, &
            !$acc&         speccomb_planck, specparm_planck, specmult_planck, jpl, &
            !$acc&         fpl, inds, indf, indm, indsp, indfp, indmp, jplp, &
            !$acc&         id000, id010, speccomb_mn2o, specparm_mn2o, specmult_mn2o, &
            !$acc&         id100, id110, id200, id210, fac000, fac100, fac200, fac010, &
            !$acc&         fac110, fac210, id001, id011, id101, id111, id201, id211, &
            !$acc&         fac001, fac101, fac201, fac011, fac111, fac211, tauself, &
            !$acc&         taufor, p0, p40, fk00, fk10, fk20, p1, jmn2o, fmn2o, &
            !$acc&         p41, fk01, fk11, fk21, temp, ratn2o, adjfac, &
            !$acc&         adjcoln2o, n2om1, n2om2, absn2o, taug)
            do iplon = 1, myim(jj)
               if (k .le. laytrop(iplon, jj)) then
                  !  --- ...  lower atmosphere loop
                  speccomb = colamt(iplon, k,1, jj) + rfrate(iplon, k,4,1, jj)*colamt(iplon, k,5, jj)
                  specparm = colamt(iplon, k,1, jj) / speccomb
                  specmult = 8.0 * min(specparm, oneminus)
                  js = 1 + int(specmult)
                  fs = mod(specmult, f_one)
                  ind0 = ((jp(iplon, k, jj)-1)*5 + (jt(iplon, k, jj)-1)) * nspa(9) + js

                  speccomb1 = colamt(iplon, k,1, jj) + rfrate(iplon, k,4,2, jj)*colamt(iplon, k,5, jj)
                  specparm1 = colamt(iplon, k,1, jj) / speccomb1
                  specmult1 = 8.0 * min(specparm1, oneminus)
                  js1 = 1 + int(specmult1)
                  fs1 = mod(specmult1, f_one)
                  ind1 = (jp(iplon, k, jj)*5 + (jt1(iplon, k, jj)-1)) * nspa(9) + js1

                  speccomb_mn2o = colamt(iplon, k,1, jj) + refrat_m_a*colamt(iplon, k,5, jj)
                  specparm_mn2o = colamt(iplon, k,1, jj) / speccomb_mn2o
                  specmult_mn2o = 8.0 * min(specparm_mn2o, oneminus)
                  jmn2o = 1 + int(specmult_mn2o)
                  fmn2o = mod(specmult_mn2o, f_one)

                  speccomb_planck = colamt(iplon, k,1, jj) + refrat_planck_a*colamt(iplon, k,5, jj)
                  specparm_planck = colamt(iplon, k,1, jj) / speccomb_planck
                  specmult_planck = 8.0 * min(specparm_planck, oneminus)
                  jpl = 1 + int(specmult_planck)
                  fpl = mod(specmult_planck, f_one)

                  inds = indself(iplon, k, jj)
                  indf = indfor(iplon, k, jj)
                  indm = indminor(iplon, k, jj)
                  indsp = inds + 1
                  indfp = indf + 1
                  indmp = indm + 1
                  jplp  = jpl  + 1
                  jmn2op= jmn2o+ 1

                  !  --- ...  in atmospheres where the amount of n2o is too great to be considered
                  !           a minor species, adjust the column amount of n2o by an empirical factor
                  !           to obtain the proper contribution.

                  temp   = coldry(iplon, k, jj) * chi_mls(4,jp(iplon, k, jj)+1)
                  ratn2o = colamt(iplon, k,4, jj) / temp
                  if (ratn2o > 1.5) then
                     adjfac = 0.5 + (ratn2o-0.5)**0.65
                     adjcoln2o = adjfac * temp
                  else
                     adjcoln2o = colamt(iplon, k,4, jj)
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

                  fac000 = fk00 * fac00(iplon, k, jj)
                  fac100 = fk10 * fac00(iplon, k, jj)
                  fac200 = fk20 * fac00(iplon, k, jj)
                  fac010 = fk00 * fac10(iplon, k, jj)
                  fac110 = fk10 * fac10(iplon, k, jj)
                  fac210 = fk20 * fac10(iplon, k, jj)

                  fac001 = fk01 * fac01(iplon, k, jj)
                  fac101 = fk11 * fac01(iplon, k, jj)
                  fac201 = fk21 * fac01(iplon, k, jj)
                  fac011 = fk01 * fac11(iplon, k, jj)
                  fac111 = fk11 * fac11(iplon, k, jj)
                  fac211 = fk21 * fac11(iplon, k, jj)
                  !$acc loop seq
                  do ig = 1, ng09
                     ib = ngb(ns09+ig)
                     tauself = selffac(iplon, k, jj)* (selfref(ig,inds) + selffrac(iplon, k, jj)         &
                     &            * (selfref(ig,indsp) - selfref(ig,inds)))
                     taufor  = forfac(iplon, k, jj) * (forref(ig,indf) + forfrac(iplon, k, jj)           &
                     &            * (forref(ig,indfp) - forref(ig,indf))) 
                     n2om1   = ka_mn2o(ig,jmn2o,indm) + fmn2o                      &
                     &            * (ka_mn2o(ig,jmn2op,indm) - ka_mn2o(ig,jmn2o,indm))
                     n2om2   = ka_mn2o(ig,jmn2o,indmp) + fmn2o                     &
                     &            * (ka_mn2o(ig,jmn2op,indmp) - ka_mn2o(ig,jmn2o,indmp))
                     absn2o  = n2om1 + minorfrac(iplon, k, jj) * (n2om2 - n2om1)

                     taug = speccomb                                    &
                     &                * (fac000*absa(ig,id000) + fac010*absa(ig,id010)  &
                     &                +  fac100*absa(ig,id100) + fac110*absa(ig,id110)  &
                     &                +  fac200*absa(ig,id200) + fac210*absa(ig,id210)) &
                     &                +     speccomb1                                   &
                     &                * (fac001*absa(ig,id001) + fac011*absa(ig,id011)  &
                     &                +  fac101*absa(ig,id101) + fac111*absa(ig,id111)  &
                     &                +  fac201*absa(ig,id201) + fac211*absa(ig,id211)) &
                     &                + tauself + taufor + adjcoln2o*absn2o            

                     fracs(iplon, k, ns09+ig, jj) = fracrefa(ig,jpl) + fpl                     &
                     &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
                     tautot(iplon, k, ns09+ig, jj) = taug + tauaer(iplon, k, ib, jj)
                  enddo
               else
                  ind0 = ((jp(iplon, k, jj)-13)*5 + (jt (iplon, k, jj)-1)) * nspb(9) + 1
                  ind1 = ((jp(iplon, k, jj)-12)*5 + (jt1(iplon, k, jj)-1)) * nspb(9) + 1

                  indm = indminor(iplon, k, jj)
                  ind0p = ind0 + 1
                  ind1p = ind1 + 1
                  indmp = indm + 1

                  !  --- ...  in atmospheres where the amount of n2o is too great to be considered
                  !           a minor species, adjust the column amount of n2o by an empirical factor
                  !           to obtain the proper contribution.

                  temp   = coldry(iplon, k, jj) * chi_mls(4,jp(iplon, k, jj)+1)
                  ratn2o = colamt(iplon, k,4, jj) / temp
                  if (ratn2o > 1.5) then
                     adjfac = 0.5 + (ratn2o - 0.5)**0.65
                     adjcoln2o = adjfac * temp
                  else
                     adjcoln2o = colamt(iplon, k,4, jj)
                  endif
                  !$acc loop seq
                  do ig = 1, ng09
                     ib = ngb(ns09+ig)
                     absn2o = kb_mn2o(ig,indm) + minorfrac(iplon, k, jj)                      &
                     &           * (kb_mn2o(ig,indmp) - kb_mn2o(ig,indm))

                     taug = colamt(iplon, k,5, jj)                                 &
                     &           * (fac00(iplon, k, jj)*absb(ig,ind0) + fac10(iplon, k, jj)*absb(ig,ind0p)    &
                     &           +  fac01(iplon, k, jj)*absb(ig,ind1) + fac11(iplon, k, jj)*absb(ig,ind1p))   &
                     &           + adjcoln2o*absn2o

                     fracs(iplon, k, ns09+ig, jj) = fracrefb(ig)
                     tautot(iplon, k, ns09+ig, jj) = taug + tauaer(iplon, k, ib, jj)
                  enddo
               end if
            end do
         end do
      end do

! ..................................
      end subroutine taugb09
! ----------------------------------

! ----------------------------------
      subroutine taugb10 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!     band 10:  1390-1480 cm-1 (low key - h2o; high key - h2o)         !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb10
      integer, intent(in) :: nlay, laytrop(ix, fulljj), ix, myim(fulljj), &
         async_id, small_ix, i2, fulljj

      integer, dimension(ix, nlay, fulljj), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(ix, nlay, fulljj), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj), intent(in):: colamt
      real (kind=kind_phys), dimension(ix, nlay,maxxsec, fulljj),intent(in):: wx

      real (kind=kind_phys), dimension(ix, nlay, nbands, fulljj), intent(in):: tauaer

      real (kind=kind_phys), dimension(ix, nlay,nrates,2, fulljj), intent(in) ::    &
     &       rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(ix, nlay, ngptlw, fulljj), intent(out) ::     &
     &       fracs, tautot
!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       ig

      real (kind=kind_phys) :: tauself, taufor, taug
      integer :: jj, iplon, iplon2, ib
!
!===> ...  begin here
!
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !$acc loop vector private(iplon, ind0, ind1, inds, indf, ind0p, &
            !$acc&         ind1p, indsp, indfp, tauself, taufor, taug)
            do iplon = 1, myim(jj)
               if (k .le. laytrop(iplon, jj)) then
                  !  --- ...  lower atmosphere loop
                  ind0 = ((jp(iplon, k, jj)-1)*5 + (jt (iplon, k, jj)-1)) * nspa(10) + 1
                  ind1 = ( jp(iplon, k, jj)   *5 + (jt1(iplon, k, jj)-1)) * nspa(10) + 1

                  inds = indself(iplon, k, jj)
                  indf = indfor(iplon, k, jj)
                  ind0p = ind0 + 1
                  ind1p = ind1 + 1
                  indsp = inds + 1
                  indfp = indf + 1
                  !$acc loop seq
                  do ig = 1, ng10
                     ib = ngb(ns10+ig)
                     tauself = selffac(iplon, k, jj) * (selfref(ig,inds) + selffrac(iplon, k, jj)        &
                     &            * (selfref(ig,indsp) - selfref(ig,inds)))
                     taufor  = forfac(iplon, k, jj) * (forref(ig,indf) + forfrac(iplon, k, jj)           &
                     &            * (forref(ig,indfp) - forref(ig,indf))) 

                     taug = colamt(iplon, k,1, jj)                                 &
                     &            * (fac00(iplon, k, jj)*absa(ig,ind0) + fac10(iplon, k, jj)*absa(ig,ind0p)   &
                     &            +  fac01(iplon, k, jj)*absa(ig,ind1) + fac11(iplon, k, jj)*absa(ig,ind1p))  &
                     &            + tauself + taufor

                     fracs(iplon, k, ns10+ig, jj) = fracrefa(ig)
                     tautot(iplon, k, ns10+ig, jj) = taug + tauaer(iplon, k, ib, jj)
                  enddo
               else
                  ind0 = ((jp(iplon, k, jj)-13)*5 + (jt (iplon, k, jj)-1)) * nspb(10) + 1
                  ind1 = ((jp(iplon, k, jj)-12)*5 + (jt1(iplon, k, jj)-1)) * nspb(10) + 1

                  indf = indfor(iplon, k, jj)
                  ind0p = ind0 + 1
                  ind1p = ind1 + 1
                  indfp = indf + 1
                  !$acc loop seq
                  do ig = 1, ng10
                     ib = ngb(ns10+ig)
                     taufor = forfac(iplon, k, jj) * (forref(ig,indf) + forfrac(iplon, k, jj)            &
                     &           * (forref(ig,indfp) - forref(ig,indf))) 

                     taug = colamt(iplon, k,1, jj)                                 &
                     &            * (fac00(iplon, k, jj)*absb(ig,ind0) + fac10(iplon, k, jj)*absb(ig,ind0p)   &
                     &            +  fac01(iplon, k, jj)*absb(ig,ind1) + fac11(iplon, k, jj)*absb(ig,ind1p))  &
                     &            + taufor

                     fracs(iplon, k, ns10+ig, jj) = fracrefb(ig)
                     tautot(iplon, k, ns10+ig, jj) = taug + tauaer(iplon, k, ib, jj)
                  enddo
               end if
            end do
         end do
      end do

! ..................................
      end subroutine taugb10
! ----------------------------------

! ----------------------------------
      subroutine taugb11 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!     band 11:  1480-1800 cm-1 (low - h2o; low minor - o2)             !
!                              (high key - h2o; high minor - o2)       !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb11
      integer, intent(in) :: nlay, laytrop(ix, fulljj), ix, myim(fulljj), &
         async_id, small_ix, i2, fulljj

      integer, dimension(ix, nlay, fulljj), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(ix, nlay, fulljj), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj), intent(in):: colamt
      real (kind=kind_phys), dimension(ix, nlay,maxxsec, fulljj),intent(in):: wx

      real (kind=kind_phys), dimension(ix, nlay, nbands, fulljj), intent(in):: tauaer

      real (kind=kind_phys), dimension(ix, nlay,nrates,2, fulljj), intent(in) ::    &
     &       rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(ix, nlay, ngptlw, fulljj), intent(out) ::     &
     &       fracs, tautot
!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       indm, indmp, ig

      real (kind=kind_phys) :: scaleo2, tauself, taufor, tauo2, taug
      integer :: jj, iplon, iplon2, ib
!
!===> ...  begin here
!
!  --- ...  minor gas mapping level :
!     lower - o2, p = 706.2720 mbar, t = 278.94 k
!     upper - o2, p = 4.758820 mbarm t = 250.85 k

      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !$acc loop vector private(iplon, ind0, ind1, inds, indf, ind0p, &
            !$acc&         ind1p, indsp, indfp, tauself, taufor, indm, indmp, scaleo2, &
            !$acc&         tauo2, taug)
            do iplon = 1, myim(jj)
               if (k .le. laytrop(iplon, jj)) then
                  !  --- ...  lower atmosphere loop
                  ind0 = ((jp(iplon, k, jj)-1)*5 + (jt (iplon, k, jj)-1)) * nspa(11) + 1
                  ind1 = ( jp(iplon, k, jj)   *5 + (jt1(iplon, k, jj)-1)) * nspa(11) + 1

                  inds = indself(iplon, k, jj)
                  indf = indfor(iplon, k, jj)
                  indm = indminor(iplon, k, jj)
                  ind0p = ind0 + 1
                  ind1p = ind1 + 1
                  indsp = inds + 1
                  indfp = indf + 1
                  indmp = indm + 1

                  scaleo2 = colamt(iplon, k,6, jj) * scaleminor(iplon, k, jj)
                  !$acc loop seq
                  do ig = 1, ng11
                     ib = ngb(ns11+ig)
                     tauself = selffac(iplon, k, jj) * (selfref(ig,inds) + selffrac(iplon, k, jj)        &
                     &            * (selfref(ig,indsp) - selfref(ig,inds)))
                     taufor  = forfac(iplon, k, jj) * (forref(ig,indf) + forfrac(iplon, k, jj)           &
                     &            * (forref(ig,indfp) - forref(ig,indf)))
                     tauo2   = scaleo2 * (ka_mo2(ig,indm) + minorfrac(iplon, k, jj)           &
                     &            * (ka_mo2(ig,indmp) - ka_mo2(ig,indm)))

                     taug = colamt(iplon, k,1, jj)                                 &
                     &            * (fac00(iplon, k, jj)*absa(ig,ind0) + fac10(iplon, k, jj)*absa(ig,ind0p)   &
                     &            +  fac01(iplon, k, jj)*absa(ig,ind1) + fac11(iplon, k, jj)*absa(ig,ind1p))  &
                     &            + tauself + taufor + tauo2

                     fracs(iplon, k, ns11+ig, jj) = fracrefa(ig)
                     tautot(iplon, k, ns11+ig, jj) = taug + tauaer(iplon, k, ib, jj)
                  enddo
               else
                  ind0 = ((jp(iplon, k, jj)-13)*5 + (jt (iplon, k, jj)-1)) * nspb(11) + 1
                  ind1 = ((jp(iplon, k, jj)-12)*5 + (jt1(iplon, k, jj)-1)) * nspb(11) + 1

                  indf = indfor(iplon, k, jj)
                  indm = indminor(iplon, k, jj)
                  ind0p = ind0 + 1
                  ind1p = ind1 + 1
                  indfp = indf + 1
                  indmp = indm + 1

                  scaleo2 = colamt(iplon, k,6, jj) * scaleminor(iplon, k, jj)
                  !$acc loop seq
                  do ig = 1, ng11
                     ib = ngb(ns11+ig)
                     taufor = forfac(iplon, k, jj) * (forref(ig,indf) + forfrac(iplon, k, jj)            &
                     &           * (forref(ig,indfp) - forref(ig,indf))) 
                     tauo2  = scaleo2 * (kb_mo2(ig,indm) + minorfrac(iplon, k, jj)            &
                     &           * (kb_mo2(ig,indmp) - kb_mo2(ig,indm)))

                     taug = colamt(iplon, k,1, jj)                                 &
                     &            * (fac00(iplon, k, jj)*absb(ig,ind0) + fac10(iplon, k, jj)*absb(ig,ind0p)   &
                     &            +  fac01(iplon, k, jj)*absb(ig,ind1) + fac11(iplon, k, jj)*absb(ig,ind1p))  &
                     &            + taufor + tauo2

                     fracs(iplon, k, ns11+ig, jj) = fracrefb(ig)
                     tautot(iplon, k, ns11+ig, jj) = taug + tauaer(iplon, k, ib, jj)
                  enddo
               end if
            end do
         end do
      end do

! ..................................
      end subroutine taugb11
! ----------------------------------

! ----------------------------------
      subroutine taugb12 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!     band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)         !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb12
      integer, intent(in) :: nlay, laytrop(ix, fulljj), ix, myim(fulljj), &
         async_id, small_ix, i2, fulljj

      integer, dimension(ix, nlay, fulljj), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(ix, nlay, fulljj), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj), intent(in):: colamt
      real (kind=kind_phys), dimension(ix, nlay,maxxsec, fulljj),intent(in):: wx

      real (kind=kind_phys), dimension(ix, nlay, nbands, fulljj), intent(in):: tauaer

      real (kind=kind_phys), dimension(ix, nlay,nrates,2, fulljj), intent(in) ::    &
     &       rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(ix, nlay, ngptlw, fulljj), intent(out) ::     &
     &       fracs, tautot
!  ---  locals:
      integer :: k, ind0, ind1, inds, indsp, indf, indfp, jpl, jplp,    &
     &       id000, id010, id100, id110, id200, id210, ig, js, js1,     &
     &       id001, id011, id101, id111, id201, id211

      real (kind=kind_phys) :: tauself, taufor, refrat_planck_a,        &
     &       speccomb,       specparm,       specmult,       fs,        &
     &       speccomb1,      specparm1,      specmult1,      fs1,       &
     &       speccomb_planck,specparm_planck,specmult_planck,fpl,       &
     &       fac000, fac100, fac200, fac010, fac110, fac210,            &
     &       fac001, fac101, fac201, fac011, fac111, fac211,            &
     &       p0, p40, fk00, fk10, fk20, p1, p41, fk01, fk11, fk21, taug
      integer :: jj, iplon, iplon2, ib
!
!===> ...  begin here
!
!  --- ...  calculate reference ratio to be used in calculation of planck
!           fraction in lower/upper atmosphere.
      refrat_planck_a = chi_mls(1,10)/chi_mls(2,10)      ! p =   174.164 mb
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !$acc loop vector private(iplon, speccomb, specparm, specmult, js, &
            !$acc&         fs, ind0, speccomb1, specparm1, specmult1, js1, fs1, ind1, &
            !$acc&         speccomb_planck, specparm_planck, specmult_planck, jpl, &
            !$acc&         fpl, inds, indf, indsp, indfp, jplp, &
            !$acc&         id000, id010, &
            !$acc&         id100, id110, id200, id210, fac000, fac100, fac200, fac010, &
            !$acc&         fac110, fac210, id001, id011, id101, id111, id201, id211, &
            !$acc&         fac001, fac101, fac201, fac011, fac111, fac211, tauself, &
            !$acc&         taufor, p0, p40, fk00, fk10, fk20, p1, &
            !$acc&         p41, fk01, fk11, fk21, taug)
            do iplon = 1, myim(jj)
               if (k .le. laytrop(iplon, jj)) then
                  !  --- ...  lower atmosphere loop
                  speccomb = colamt(iplon, k,1, jj) + rfrate(iplon, k,1,1, jj)*colamt(iplon, k,2, jj)
                  specparm = colamt(iplon, k,1, jj) / speccomb
                  specmult = 8.0 * min(specparm, oneminus)
                  js = 1 + int(specmult)
                  fs = mod(specmult, f_one)
                  ind0 = ((jp(iplon, k, jj)-1)*5 + (jt(iplon, k, jj)-1)) * nspa(12) + js

                  speccomb1 = colamt(iplon, k,1, jj) + rfrate(iplon, k,1,2, jj)*colamt(iplon, k,2, jj)
                  specparm1 = colamt(iplon, k,1, jj) / speccomb1
                  specmult1 = 8.0 * min(specparm1, oneminus)
                  js1 = 1 + int(specmult1)
                  fs1 = mod(specmult1, f_one)
                  ind1 = (jp(iplon, k, jj)*5 + (jt1(iplon, k, jj)-1)) * nspa(12) + js1

                  speccomb_planck = colamt(iplon, k,1, jj) + refrat_planck_a*colamt(iplon, k,2, jj)
                  specparm_planck = colamt(iplon, k,1, jj) / speccomb_planck
                  if (specparm_planck >= oneminus) specparm_planck=oneminus
                  specmult_planck = 8.0 * specparm_planck
                  jpl = 1 + int(specmult_planck)
                  fpl = mod(specmult_planck, f_one)

                  inds = indself(iplon, k, jj)
                  indf = indfor(iplon, k, jj)
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

                  fac000 = fk00 * fac00(iplon, k, jj)
                  fac100 = fk10 * fac00(iplon, k, jj)
                  fac200 = fk20 * fac00(iplon, k, jj)
                  fac010 = fk00 * fac10(iplon, k, jj)
                  fac110 = fk10 * fac10(iplon, k, jj)
                  fac210 = fk20 * fac10(iplon, k, jj)

                  fac001 = fk01 * fac01(iplon, k, jj)
                  fac101 = fk11 * fac01(iplon, k, jj)
                  fac201 = fk21 * fac01(iplon, k, jj)
                  fac011 = fk01 * fac11(iplon, k, jj)
                  fac111 = fk11 * fac11(iplon, k, jj)
                  fac211 = fk21 * fac11(iplon, k, jj)
                  !$acc loop seq
                  do ig = 1, ng12
                     ib = ngb(ns12+ig)
                     tauself = selffac(iplon, k, jj)* (selfref(ig,inds) + selffrac(iplon, k, jj)         &
                     &            * (selfref(ig,indsp) - selfref(ig,inds)))
                     taufor  = forfac(iplon, k, jj) * (forref(ig,indf) + forfrac(iplon, k, jj)           &
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

                     fracs(iplon, k, ns12+ig, jj) = fracrefa(ig,jpl) + fpl                     &
                     &                     *(fracrefa(ig,jplp) - fracrefa(ig,jpl))
                     tautot(iplon, k, ns12+ig, jj) = taug + tauaer(iplon, k, ib, jj)
                  enddo
               else
                  !$acc loop seq
                  do ig = 1, ng12
                     ib = ngb(ns12+ig)
                     taug = f_zero
                     fracs(iplon, k, ns12+ig, jj) = f_zero
                     tautot(iplon, k, ns12+ig, jj) = taug + tauaer(iplon, k, ib, jj)
                  enddo
               end if
            end do
         end do
      end do

! ..................................
      end subroutine taugb12
! ----------------------------------

! ----------------------------------
      subroutine taugb13 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!     band 13:  2080-2250 cm-1 (low key-h2o,n2o; high minor-o3 minor)  !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb13
      integer, intent(in) :: nlay, laytrop(ix, fulljj), ix, myim(fulljj), &
         async_id, small_ix, i2, fulljj

      integer, dimension(ix, nlay, fulljj), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(ix, nlay, fulljj), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj), intent(in):: colamt
      real (kind=kind_phys), dimension(ix, nlay,maxxsec, fulljj),intent(in):: wx

      real (kind=kind_phys), dimension(ix, nlay, nbands, fulljj), intent(in):: tauaer

      real (kind=kind_phys), dimension(ix, nlay,nrates,2, fulljj), intent(in) ::    &
     &       rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(ix, nlay, ngptlw, fulljj), intent(out) ::     &
     &       fracs, tautot
!  ---  locals:
      integer :: k, ind0, ind1, inds, indsp, indf, indfp, indm, indmp,  &
     &       id000, id010, id100, id110, id200, id210, jmco2, jpl,      &
     &       id001, id011, id101, id111, id201, id211, jmco2p, jplp,    &
     &       jmco, jmcop, ig, js, js1

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
      integer :: jj, iplon, iplon2, ib
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
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !$acc loop vector private(iplon, speccomb, specparm, specmult, js, &
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
            !$acc&         com1, com2, absco, taug)
            do iplon = 1, myim(jj)
               if (k .le. laytrop(iplon, jj)) then
                  !  --- ...  lower atmosphere loop
                  speccomb = colamt(iplon, k,1, jj) + rfrate(iplon, k,3,1, jj)*colamt(iplon, k,4, jj)
                  specparm = colamt(iplon, k,1, jj) / speccomb
                  specmult = 8.0 * min(specparm, oneminus)
                  js = 1 + int(specmult)
                  fs = mod(specmult, f_one)
                  ind0 = ((jp(iplon, k, jj)-1)*5 + (jt(iplon, k, jj)-1)) * nspa(13) + js

                  speccomb1 = colamt(iplon, k,1, jj) + rfrate(iplon, k,3,2, jj)*colamt(iplon, k,4, jj)
                  specparm1 = colamt(iplon, k,1, jj) / speccomb1
                  specmult1 = 8.0 * min(specparm1, oneminus)
                  js1 = 1 + int(specmult1)
                  fs1 = mod(specmult1, f_one)
                  ind1 = (jp(iplon, k, jj)*5 + (jt1(iplon, k, jj)-1)) * nspa(13) + js1

                  speccomb_mco2 = colamt(iplon, k,1, jj) + refrat_m_a*colamt(iplon, k,4, jj)
                  specparm_mco2 = colamt(iplon, k,1, jj) / speccomb_mco2
                  specmult_mco2 = 8.0 * min(specparm_mco2, oneminus)
                  jmco2 = 1 + int(specmult_mco2)
                  fmco2 = mod(specmult_mco2, f_one)

                  !  --- ...  in atmospheres where the amount of co2 is too great to be considered
                  !           a minor species, adjust the column amount of co2 by an empirical factor
                  !           to obtain the proper contribution.

                  speccomb_mco = colamt(iplon, k,1, jj) + refrat_m_a3*colamt(iplon, k,4, jj)
                  specparm_mco = colamt(iplon, k,1, jj) / speccomb_mco
                  specmult_mco = 8.0 * min(specparm_mco, oneminus)
                  jmco = 1 + int(specmult_mco)
                  fmco = mod(specmult_mco, f_one)

                  speccomb_planck = colamt(iplon, k,1, jj) + refrat_planck_a*colamt(iplon, k,4, jj)
                  specparm_planck = colamt(iplon, k,1, jj) / speccomb_planck
                  specmult_planck = 8.0 * min(specparm_planck, oneminus)
                  jpl = 1 + int(specmult_planck)
                  fpl = mod(specmult_planck, f_one)

                  inds = indself(iplon, k, jj)
                  indf = indfor(iplon, k, jj)
                  indm = indminor(iplon, k, jj)
                  indsp = inds + 1
                  indfp = indf + 1
                  indmp = indm + 1
                  jplp  = jpl  + 1
                  jmco2p= jmco2+ 1
                  jmcop = jmco + 1

                  !  --- ...  in atmospheres where the amount of co2 is too great to be considered
                  !           a minor species, adjust the column amount of co2 by an empirical factor
                  !           to obtain the proper contribution.

                  temp   = coldry(iplon, k, jj) * 3.55e-4
                  ratco2 = colamt(iplon, k,2, jj) / temp
                  if (ratco2 > 3.0) then
                     adjfac = 2.0 + (ratco2-2.0)**0.68
                     adjcolco2 = adjfac * temp
                  else
                     adjcolco2 = colamt(iplon, k,2, jj)
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

                  fac000 = fk00 * fac00(iplon, k, jj)
                  fac100 = fk10 * fac00(iplon, k, jj)
                  fac200 = fk20 * fac00(iplon, k, jj)
                  fac010 = fk00 * fac10(iplon, k, jj)
                  fac110 = fk10 * fac10(iplon, k, jj)
                  fac210 = fk20 * fac10(iplon, k, jj)

                  fac001 = fk01 * fac01(iplon, k, jj)
                  fac101 = fk11 * fac01(iplon, k, jj)
                  fac201 = fk21 * fac01(iplon, k, jj)
                  fac011 = fk01 * fac11(iplon, k, jj)
                  fac111 = fk11 * fac11(iplon, k, jj)
                  fac211 = fk21 * fac11(iplon, k, jj)
                  !$acc loop seq
                  do ig = 1, ng13
                     ib = ngb(ns13+ig)
                     tauself = selffac(iplon, k, jj)* (selfref(ig,inds) + selffrac(iplon, k, jj)         &
                     &            * (selfref(ig,indsp) - selfref(ig,inds)))
                     taufor  = forfac(iplon, k, jj) * (forref(ig,indf) + forfrac(iplon, k, jj)           &
                     &            * (forref(ig,indfp) - forref(ig,indf))) 
                     co2m1   = ka_mco2(ig,jmco2,indm) + fmco2                      &
                     &            * (ka_mco2(ig,jmco2p,indm) - ka_mco2(ig,jmco2,indm))
                     co2m2   = ka_mco2(ig,jmco2,indmp) + fmco2                     &
                     &            * (ka_mco2(ig,jmco2p,indmp) - ka_mco2(ig,jmco2,indmp))
                     absco2  = co2m1 + minorfrac(iplon, k, jj) * (co2m2 - co2m1)
                     com1    = ka_mco(ig,jmco,indm) + fmco                         &
                     &            * (ka_mco(ig,jmcop,indm) - ka_mco(ig,jmco,indm))
                     com2    = ka_mco(ig,jmco,indmp) + fmco                        &
                     &            * (ka_mco(ig,jmcop,indmp) - ka_mco(ig,jmco,indmp))
                     absco   = com1 + minorfrac(iplon, k, jj) * (com2 - com1)

                     taug = speccomb                                    &
                     &                * (fac000*absa(ig,id000) + fac010*absa(ig,id010)  &
                     &                +  fac100*absa(ig,id100) + fac110*absa(ig,id110)  &
                     &                +  fac200*absa(ig,id200) + fac210*absa(ig,id210)) &
                     &                +     speccomb1                                   &
                     &                * (fac001*absa(ig,id001) + fac011*absa(ig,id011)  &
                     &                +  fac101*absa(ig,id101) + fac111*absa(ig,id111)  &
                     &                +  fac201*absa(ig,id201) + fac211*absa(ig,id211)) &
                     &                + tauself + taufor + adjcolco2*absco2             &
                     &                + colamt(iplon, k,7, jj)*absco

                     fracs(iplon, k, ns13+ig, jj) = fracrefa(ig,jpl) + fpl                     &
                     &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
                     tautot(iplon, k, ns13+ig, jj) = taug + tauaer(iplon, k, ib, jj)
                  enddo
               else
                  indm = indminor(iplon, k, jj)
                  indmp = indm + 1
                  !$acc loop seq
                  do ig = 1, ng13
                     ib = ngb(ns13+ig)
                     abso3 = kb_mo3(ig,indm) + minorfrac(iplon, k, jj)                        &
                     &          * (kb_mo3(ig,indmp) - kb_mo3(ig,indm))

                     taug = colamt(iplon, k,3, jj)*abso3

                     fracs(iplon, k, ns13+ig, jj) =  fracrefb(ig)
                     tautot(iplon, k, ns13+ig, jj) = taug + tauaer(iplon, k, ib, jj)
                  enddo
               end if
            end do
         end do
      end do

! ..................................
      end subroutine taugb13
! ----------------------------------

! ----------------------------------
      subroutine taugb14 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!     band 14:  2250-2380 cm-1 (low - co2; high - co2)                 !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb14
      integer, intent(in) :: nlay, laytrop(ix, fulljj), ix, myim(fulljj), &
         async_id, small_ix, i2, fulljj

      integer, dimension(ix, nlay, fulljj), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(ix, nlay, fulljj), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj), intent(in):: colamt
      real (kind=kind_phys), dimension(ix, nlay,maxxsec, fulljj),intent(in):: wx

      real (kind=kind_phys), dimension(ix, nlay, nbands, fulljj), intent(in):: tauaer

      real (kind=kind_phys), dimension(ix, nlay,nrates,2, fulljj), intent(in) ::    &
     &       rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(ix, nlay, ngptlw, fulljj), intent(out) ::     &
     &       fracs, tautot
!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       ig

      real (kind=kind_phys) :: tauself, taufor, taug
      integer :: jj, iplon, iplon2, ib
!
!===> ...  begin here
!
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !$acc loop vector private(iplon, ind0, ind1, inds, indf, ind0p, &
            !$acc&         ind1p, indsp, indfp, tauself, taufor, taug)
            do iplon = 1, myim(jj)
               if (k .le. laytrop(iplon, jj)) then
                  !  --- ...  lower atmosphere loop
                  ind0 = ((jp(iplon, k, jj)-1)*5 + (jt (iplon, k, jj)-1)) * nspa(14) + 1
                  ind1 = ( jp(iplon, k, jj)   *5 + (jt1(iplon, k, jj)-1)) * nspa(14) + 1

                  inds = indself(iplon, k, jj)
                  indf = indfor(iplon, k, jj)
                  ind0p = ind0 + 1
                  ind1p = ind1 + 1
                  indsp = inds + 1
                  indfp = indf + 1
                  !$acc loop seq
                  do ig = 1, ng14
                     ib = ngb(ns14+ig)
                     tauself = selffac(iplon, k, jj) * (selfref(ig,inds) + selffrac(iplon, k, jj)        &
                     &            * (selfref(ig,indsp) - selfref(ig,inds)))
                     taufor  = forfac(iplon, k, jj) * (forref(ig,indf) + forfrac(iplon, k, jj)           &
                     &            * (forref(ig,indfp) - forref(ig,indf))) 

                     taug = colamt(iplon, k,2, jj)                                 &
                     &            * (fac00(iplon, k, jj)*absa(ig,ind0) + fac10(iplon, k, jj)*absa(ig,ind0p)   &
                     &            +  fac01(iplon, k, jj)*absa(ig,ind1) + fac11(iplon, k, jj)*absa(ig,ind1p))  &
                     &            + tauself + taufor

                     fracs(iplon, k, ns14+ig, jj) = fracrefa(ig)
                     tautot(iplon, k, ns14+ig, jj) = taug + tauaer(iplon, k, ib, jj)

                  enddo
               else
                  ind0 = ((jp(iplon, k, jj)-13)*5 + (jt (iplon, k, jj)-1)) * nspb(14) + 1
                  ind1 = ((jp(iplon, k, jj)-12)*5 + (jt1(iplon, k, jj)-1)) * nspb(14) + 1

                  ind0p = ind0 + 1
                  ind1p = ind1 + 1
                  !$acc loop seq
                  do ig = 1, ng14
                     ib = ngb(ns14+ig)
                     taug = colamt(iplon, k,2, jj)                                 &
                     &             * (fac00(iplon, k, jj)*absb(ig,ind0) + fac10(iplon, k, jj)*absb(ig,ind0p)  &
                     &             +  fac01(iplon, k, jj)*absb(ig,ind1) + fac11(iplon, k, jj)*absb(ig,ind1p))

                     fracs(iplon, k, ns14+ig, jj) = fracrefb(ig)
                     tautot(iplon, k, ns14+ig, jj) = taug + tauaer(iplon, k, ib, jj)
                  enddo
               end if
            end do
         end do
      end do

! ..................................
      end subroutine taugb14
! ----------------------------------

! ----------------------------------
      subroutine taugb15 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!     band 15:  2380-2600 cm-1 (low - n2o,co2; low minor - n2)         !
!                              (high - nothing)                        !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb15
      integer, intent(in) :: nlay, laytrop(ix, fulljj), ix, myim(fulljj), &
         async_id, small_ix, i2, fulljj

      integer, dimension(ix, nlay, fulljj), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(ix, nlay, fulljj), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj), intent(in):: colamt
      real (kind=kind_phys), dimension(ix, nlay,maxxsec, fulljj),intent(in):: wx

      real (kind=kind_phys), dimension(ix, nlay, nbands, fulljj), intent(in):: tauaer

      real (kind=kind_phys), dimension(ix, nlay,nrates,2, fulljj), intent(in) ::    &
     &       rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(ix, nlay, ngptlw, fulljj), intent(out) ::     &
     &       fracs, tautot
!  ---  locals:
      integer :: k, ind0, ind1, inds, indsp, indf, indfp, indm, indmp,  &
     &       id000, id010, id100, id110, id200, id210, jpl, jplp,       &
     &       id001, id011, id101, id111, id201, id211, jmn2, jmn2p,     &
     &       ig, js, js1

      real (kind=kind_phys) :: scalen2, tauself, taufor,                &
     &       speccomb,       specparm,       specmult,       fs,        &
     &       speccomb1,      specparm1,      specmult1,      fs1,       &
     &       speccomb_mn2,   specparm_mn2,   specmult_mn2,   fmn2,      &
     &       speccomb_planck,specparm_planck,specmult_planck,fpl,       &
     &       refrat_planck_a, refrat_m_a, n2m1, n2m2, taun2,            &
     &       fac000, fac100, fac200, fac010, fac110, fac210,            &
     &       fac001, fac101, fac201, fac011, fac111, fac211,            &
     &       p0, p40, fk00, fk10, fk20, p1, p41, fk01, fk11, fk21, taug
      integer :: jj, iplon, iplon2, ib
!
!===> ...  begin here
!
!  --- ...  minor gas mapping level :
!     lower - nitrogen continuum, p = 1053., t = 294.

!  --- ...  calculate reference ratio to be used in calculation of planck
!           fraction in lower atmosphere.
      refrat_planck_a = chi_mls(4,1)/chi_mls(2,1)      ! p = 1053. mb (level 1)
      refrat_m_a = chi_mls(4,1)/chi_mls(2,1)           ! p = 1053. mb
      !$acc parallel loop gang collapse(2) async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !$acc loop vector private(iplon, speccomb, specparm, specmult, js, &
            !$acc&         fs, ind0, speccomb1, specparm1, specmult1, js1, fs1, ind1, &
            !$acc&         speccomb_planck, specparm_planck, specmult_planck, jpl, &
            !$acc&         fpl, inds, indf, indm, indsp, indfp, indmp, jplp, &
            !$acc&         id000, id010, speccomb_mn2, specparm_mn2, specmult_mn2, &
            !$acc&         id100, id110, id200, id210, fac000, fac100, fac200, fac010, &
            !$acc&         fac110, fac210, id001, id011, id101, id111, id201, id211, &
            !$acc&         fac001, fac101, fac201, fac011, fac111, fac211, tauself, &
            !$acc&         taufor, p0, p40, fk00, fk10, fk20, p1, jmn2, fmn2, &
            !$acc&         p41, fk01, fk11, fk21, scalen2, jmn2p, n2m1, n2m2, taun2, taug)
            do iplon = 1, myim(jj)
               if (k .le. laytrop(iplon, jj)) then
                  !  --- ...  lower atmosphere loop
                  speccomb = colamt(iplon, k,4, jj) + rfrate(iplon, k,5,1, jj)*colamt(iplon, k,2, jj)
                  specparm = colamt(iplon, k,4, jj) / speccomb
                  specmult = 8.0 * min(specparm, oneminus)
                  js = 1 + int(specmult)
                  fs = mod(specmult, f_one)
                  ind0 = ((jp(iplon, k, jj)-1)*5 + (jt(iplon, k, jj)-1)) * nspa(15) + js

                  speccomb1 = colamt(iplon, k,4, jj) + rfrate(iplon, k,5,2, jj)*colamt(iplon, k,2, jj)
                  specparm1 = colamt(iplon, k,4, jj) / speccomb1
                  specmult1 = 8.0 * min(specparm1, oneminus)
                  js1 = 1 + int(specmult1)
                  fs1 = mod(specmult1, f_one)
                  ind1 = (jp(iplon, k, jj)*5 + (jt1(iplon, k, jj)-1)) * nspa(15) + js1

                  speccomb_mn2 = colamt(iplon, k,4, jj) + refrat_m_a*colamt(iplon, k,2, jj)
                  specparm_mn2 = colamt(iplon, k,4, jj) / speccomb_mn2
                  specmult_mn2 = 8.0 * min(specparm_mn2, oneminus)
                  jmn2 = 1 + int(specmult_mn2)
                  fmn2 = mod(specmult_mn2, f_one)

                  speccomb_planck = colamt(iplon, k,4, jj) + refrat_planck_a*colamt(iplon, k,2, jj)
                  specparm_planck = colamt(iplon, k,4, jj) / speccomb_planck
                  specmult_planck = 8.0 * min(specparm_planck, oneminus)
                  jpl = 1 + int(specmult_planck)
                  fpl = mod(specmult_planck, f_one)

                  scalen2 = colbrd(iplon, k, jj) * scaleminor(iplon, k, jj)

                  inds = indself(iplon, k, jj)
                  indf = indfor(iplon, k, jj)
                  indm = indminor(iplon, k, jj)
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

                  fac000 = fk00 * fac00(iplon, k, jj)
                  fac100 = fk10 * fac00(iplon, k, jj)
                  fac200 = fk20 * fac00(iplon, k, jj)
                  fac010 = fk00 * fac10(iplon, k, jj)
                  fac110 = fk10 * fac10(iplon, k, jj)
                  fac210 = fk20 * fac10(iplon, k, jj)

                  fac001 = fk01 * fac01(iplon, k, jj)
                  fac101 = fk11 * fac01(iplon, k, jj)
                  fac201 = fk21 * fac01(iplon, k, jj)
                  fac011 = fk01 * fac11(iplon, k, jj)
                  fac111 = fk11 * fac11(iplon, k, jj)
                  fac211 = fk21 * fac11(iplon, k, jj)
                  !$acc loop seq
                  do ig = 1, ng15
                     ib = ngb(ns15+ig)
                     tauself = selffac(iplon, k, jj)* (selfref(ig,inds) + selffrac(iplon, k, jj)         &
                     &            * (selfref(ig,indsp) - selfref(ig,inds)))
                     taufor  = forfac(iplon, k, jj) * (forref(ig,indf) + forfrac(iplon, k, jj)           &
                     &            * (forref(ig,indfp) - forref(ig,indf))) 
                     n2m1    = ka_mn2(ig,jmn2,indm) + fmn2                         &
                     &            * (ka_mn2(ig,jmn2p,indm) - ka_mn2(ig,jmn2,indm))
                     n2m2    = ka_mn2(ig,jmn2,indmp) + fmn2                        &
                     &            * (ka_mn2(ig,jmn2p,indmp) - ka_mn2(ig,jmn2,indmp))
                     taun2   = scalen2 * (n2m1 + minorfrac(iplon, k, jj) * (n2m2 - n2m1))

                     taug = speccomb                                    &
                     &                * (fac000*absa(ig,id000) + fac010*absa(ig,id010)  &
                     &                +  fac100*absa(ig,id100) + fac110*absa(ig,id110)  &
                     &                +  fac200*absa(ig,id200) + fac210*absa(ig,id210)) &
                     &                +     speccomb1                                   &
                     &                * (fac001*absa(ig,id001) + fac011*absa(ig,id011)  &
                     &                +  fac101*absa(ig,id101) + fac111*absa(ig,id111)  &
                     &                +  fac201*absa(ig,id201) + fac211*absa(ig,id211)) &
                     &                + tauself + taufor + taun2

                     fracs(iplon, k, ns15+ig, jj) = fracrefa(ig,jpl) + fpl                     &
                     &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
                     tautot(iplon, k, ns15+ig, jj) = taug + tauaer(iplon, k, ib, jj)

                  enddo
               else
                  !$acc loop seq
                  do ig = 1, ng15
                     ib = ngb(ns15+ig)
                     taug = f_zero

                     fracs(iplon, k, ns15+ig, jj) = f_zero
                     tautot(iplon, k, ns15+ig, jj) = taug + tauaer(iplon, k, ib, jj)

                  enddo
               end if
            end do
         end do
      end do

! ..................................
      end subroutine taugb15
! ----------------------------------

! ----------------------------------
      subroutine taugb16 &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay, ix, myim, small_ix, i2, async_id, fulljj,                    &
!  ---  outputs:
     &       fracs, tautot                                        &
     &     )
! ..................................

!  ------------------------------------------------------------------  !
!     band 16:  2600-3250 cm-1 (low key- h2o,ch4; high key - ch4)      !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb16
      integer, intent(in) :: nlay, laytrop(ix, fulljj), ix, myim(fulljj), &
         async_id, small_ix, i2, fulljj

      integer, dimension(ix, nlay, fulljj), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(ix, nlay, fulljj), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2
      real (kind=kind_phys), dimension(ix, nlay,maxgas, fulljj), intent(in):: colamt
      real (kind=kind_phys), dimension(ix, nlay,maxxsec, fulljj),intent(in):: wx

      real (kind=kind_phys), dimension(ix, nlay, nbands, fulljj), intent(in):: tauaer

      real (kind=kind_phys), dimension(ix, nlay,nrates,2, fulljj), intent(in) ::    &
     &       rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(ix, nlay, ngptlw, fulljj), intent(out) ::     &
     &       fracs, tautot
!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       id000, id010, id100, id110, id200, id210, jpl, jplp,       &
     &       id001, id011, id101, id111, id201, id211, ig, js, js1

      real (kind=kind_phys) :: tauself, taufor, refrat_planck_a,        &
     &       speccomb,       specparm,       specmult,       fs,        &
     &       speccomb1,      specparm1,      specmult1,      fs1,       &
     &       speccomb_planck,specparm_planck,specmult_planck,fpl,       &
     &       fac000, fac100, fac200, fac010, fac110, fac210,            &
     &       fac001, fac101, fac201, fac011, fac111, fac211,            &
     &       p0, p40, fk00, fk10, fk20, p1, p41, fk01, fk11, fk21, taug
      integer :: jj, iplon, iplon2, ib
!
!===> ...  begin here
!
!  --- ...  calculate reference ratio to be used in calculation of planck
!           fraction in lower atmosphere.
      refrat_planck_a = chi_mls(1,6)/chi_mls(6,6)        ! p = 387. mb (level 6)
      !$acc parallel loop gang collapse(2)async(async_id)
      do jj = 1, fulljj
         do k = 1, nlay
            !$acc loop vector private(iplon, speccomb, specparm, specmult, js, &
            !$acc&         fs, ind0, speccomb1, specparm1, specmult1, js1, fs1, ind1, &
            !$acc&         speccomb_planck, specparm_planck, specmult_planck, jpl, &
            !$acc&         fpl, inds, indf, indsp, indfp, jplp, id000, id010, &
            !$acc&         id100, id110, id200, id210, fac000, fac100, fac200, fac010, &
            !$acc&         fac110, fac210, id001, id011, id101, id111, id201, id211, &
            !$acc&         fac001, fac101, fac201, fac011, fac111, fac211, tauself, &
            !$acc&         taufor, p0, p40, fk00, fk10, fk20, p1, &
            !$acc&         p41, fk01, fk11, fk21, taug) 
            do iplon = 1, myim(jj)
               if (k .le. laytrop(iplon, jj)) then
                  !  --- ...  lower atmosphere loop
                  speccomb = colamt(iplon, k,1, jj) + rfrate(iplon, k,4,1, jj)*colamt(iplon, k,5, jj)
                  specparm = colamt(iplon, k,1, jj) / speccomb
                  specmult = 8.0 * min(specparm, oneminus)
                  js = 1 + int(specmult)
                  fs = mod(specmult, f_one)
                  ind0 = ((jp(iplon, k, jj)-1)*5 + (jt(iplon, k, jj)-1)) * nspa(16) + js

                  speccomb1 = colamt(iplon, k,1, jj) + rfrate(iplon, k,4,2, jj)*colamt(iplon, k,5, jj)
                  specparm1 = colamt(iplon, k,1, jj) / speccomb1
                  specmult1 = 8.0 * min(specparm1, oneminus)
                  js1 = 1 + int(specmult1)
                  fs1 = mod(specmult1, f_one)
                  ind1 = (jp(iplon, k, jj)*5 + (jt1(iplon, k, jj)-1)) * nspa(16) + js1

                  speccomb_planck = colamt(iplon, k,1, jj) + refrat_planck_a*colamt(iplon, k,5, jj)
                  specparm_planck = colamt(iplon, k,1, jj) / speccomb_planck
                  specmult_planck = 8.0 * min(specparm_planck, oneminus)
                  jpl = 1 + int(specmult_planck)
                  fpl = mod(specmult_planck, f_one)

                  inds = indself(iplon, k, jj)
                  indf = indfor(iplon, k, jj)
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

                  fac000 = fk00 * fac00(iplon, k, jj)
                  fac100 = fk10 * fac00(iplon, k, jj)
                  fac200 = fk20 * fac00(iplon, k, jj)
                  fac010 = fk00 * fac10(iplon, k, jj)
                  fac110 = fk10 * fac10(iplon, k, jj)
                  fac210 = fk20 * fac10(iplon, k, jj)

                  fac001 = fk01 * fac01(iplon, k, jj)
                  fac101 = fk11 * fac01(iplon, k, jj)
                  fac201 = fk21 * fac01(iplon, k, jj)
                  fac011 = fk01 * fac11(iplon, k, jj)
                  fac111 = fk11 * fac11(iplon, k, jj)
                  fac211 = fk21 * fac11(iplon, k, jj)
                  !$acc loop seq
                  do ig = 1, ng16
                     ib = ngb(ns16+ig)
                     tauself = selffac(iplon, k, jj)* (selfref(ig,inds) + selffrac(iplon, k, jj)         &
                     &            * (selfref(ig,indsp) - selfref(ig,inds)))
                     taufor  = forfac(iplon, k, jj) * (forref(ig,indf) + forfrac(iplon, k, jj)           &
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

                     fracs(iplon, k, ns16+ig, jj) = fracrefa(ig,jpl) + fpl                     &
                     &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
                     tautot(iplon, k, ns16+ig, jj) = taug + tauaer(iplon, k, ib, jj)

                  enddo
               else
                  ind0 = ((jp(iplon, k, jj)-13)*5 + (jt (iplon, k, jj)-1)) * nspb(16) + 1
                  ind1 = ((jp(iplon, k, jj)-12)*5 + (jt1(iplon, k, jj)-1)) * nspb(16) + 1

                  ind0p = ind0 + 1
                  ind1p = ind1 + 1
                  !$acc loop seq
                  do ig = 1, ng16
                     ib = ngb(ns16+ig)
                     taug = colamt(iplon, k,5, jj)                                 &
                     &           * (fac00(iplon, k, jj)*absb(ig,ind0) + fac10(iplon, k, jj)*absb(ig,ind0p)    &
                     &           +  fac01(iplon, k, jj)*absb(ig,ind1) + fac11(iplon, k, jj)*absb(ig,ind1p))

                     fracs(iplon, k, ns16+ig, jj) = fracrefb(ig)
                     tautot(iplon, k, ns16+ig, jj) = taug + tauaer(iplon, k, ib, jj)

                  enddo
               end if
            end do
         end do
      end do

! ..................................
      end subroutine taugb16
! ----------------------------------



!
!........................................!
      end module module_radlw_main_gpu       !
!========================================!

