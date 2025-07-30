!CWB2016
!ocl nosimd
!!!!!              module_radiation_clouds description             !!!!!
!!!!!  ==========================================================  !!!!!
!                                                                      !
!    the 'radiation_clouds.f' contains:                                !
!                                                                      !
!       'module_radiation_clouds' ---  compute cloud related quantities!
!                for radiation computations                            !
!                                                                      !
!    the following are the externally accessable subroutines:          !
!                                                                      !
!       'cld_init'           --- initialization routine                !
!          inputs:                                                     !
!           ( si, nlay, me, myrank )                                   !
!          outputs:                                                    !
!           ( none )                                                   !
!                                                                      !
!       'progcld1'           --- zhao/moorthi prognostic cloud scheme  !
!          inputs:                                                     !
!           (plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,                   !
!            xlat,xlon,slmsk,                                          !
!            ix, nlay, nlp1,                                           !
!          outputs:                                                    !
!            clouds,clds,mtop,mbot)                                    !
!                                                                      !
!       'progcld2'           --- ferrier prognostic cloud microphysics !
!          inputs:                                                     !
!           (plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,                   !
!            xlat,xlon,slmsk, f_ice,f_rain,r_rime,flgmin,              !
!            ix, nlay, nlp1,                                           !
!          outputs:                                                    !
!            clouds,clds,mtop,mbot)                                    !
!       'progcld3'           --- zhao/moorthi prognostic cloud + pdfcld!
!          inputs:                                                     !
!           (plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw, cnvw,cnvc,        !
!            xlat,xlon,slmsk,                                          !
!            ix, nlay, nlp1,                                           !
!            deltaq,sup,kdt,me,                                        !
!          outputs:                                                    !
!            clouds,clds,mtop,mbot)                                    !
!                                                                      !
!       'diagcld1'           --- diagnostic cloud calc routine         !
!          inputs:                                                     !
!           (plyr,plvl,tlyr,rhly,vvel,cv,cvt,cvb,                      !
!            xlat,xlon,slmsk,                                          !
!            ix, nlay, nlp1,                                           !
!          outputs:                                                    !
!            clouds,clds,mtop,mbot)                                    !
!                                                                      !
!    internal accessable only subroutines:                             !
!       'gethml'             --- get diagnostic hi, mid, low clouds    !
!                                                                      !
!       'rhtable'            --- rh lookup table for diag cloud scheme !
!                                                                      !
!                                                                      !
!    cloud array description:                                          !
!                ---  for prognostic cloud: icldflg=1  ---             !
!          clouds(:,:,1)  -  layer total cloud fraction                !
!          clouds(:,:,2)  -  layer cloud liq water path                !
!          clouds(:,:,3)  -  mean effective radius for liquid cloud    !
!          clouds(:,:,4)  -  layer cloud ice water path                !
!          clouds(:,:,5)  -  mean effective radius for ice cloud       !
!          clouds(:,:,6)  -  layer rain drop water path                !
!          clouds(:,:,7)  -  mean effective radius for rain drop       !
!   **     clouds(:,:,8)  -  layer snow flake water path               !
!          clouds(:,:,9)  -  mean effective radius for snow flake      !
!   ** fu's scheme need to be normalized by snow density (g/m**3/1.0e6)!
!                ---  for diagnostic cloud: icldflg=0  ---             !
!          clouds(:,:,1)  -  layer total cloud fraction                !
!          clouds(:,:,2)  -  layer cloud optical depth                 !
!          clouds(:,:,3)  -  layer cloud single scattering albedo      !
!          clouds(:,:,4)  -  layer cloud asymmetry factor              !
!                                                                      !
!    external modules referenced:                                      !
!                                                                      !
!       'module physpara'            in 'physpara.f'                   !
!       'module physcons'            in 'physcons.f'                   !
!       'module module_microphysics' in 'module_bfmicrophysics.f'      !
!                                                                      !
!                                                                      !
! program history log:                                                 !
!      nov 1992,   y.h., k.a.c, a.k. - cloud parameterization          !
!        'cldjms' patterned after slingo and slingo's work (jgr,       !
!        1992), stratiform clouds are allowed in any layer except      !
!        the surface and upper stratosphere. the relative humidity     !
!        criterion may cery in different model layers.                 !
!      mar 1993,   kenneth campana   - created original crhtab tables  !
!        for critical rh look up references.
!      feb 1994,   kenneth campana   - modified to use only one table  !
!        for all forecast hours.                                       !
!      oct 1995,   kenneth campana   - tuned cloud rh curves           !
!        rh-cld relation from tables created using mitchell-hahn       !
!        tuning technique on airforce rtneph observations.             !
!      nov 1995,   kenneth campana   - the bl relationships used       !
!        below llyr, except in marine stratus regions.                 !
!      apr 1996,   kenneth campana   - save bl cld amt in cld(,5)      !
!      aug 1997,   kenneth campana   - smooth out last bunch of bins   !
!        of the rh look up tables                                      !
!      dec 1998,   s. moorthi        - added prognostic cloud method   !
!      apr 2003,   yu-tai hou        - rewritten in frotran 90         !
!        modulized form 'module_rad_clouds' from combining the original!
!        subroutines 'cldjms', 'cldprp', and 'gcljms'. and seperated   !
!        prognostic and diagnostic methods into two packages.          !
!      --- 2003,   s. moorthi        - adapted b. ferrier's prognostic !
!        cloud scheme to ncep gfs model as additional option.          !
!      apr 2004,   yu-tai hou        - separated calculation of the    !
!        averaged h,m,l,bl cloud amounts from each of the cld schemes  !
!        to become an shared individule subprogram 'gethml'.           !
!      may 2004,   yu-tai hou        - rewritten ferrier's scheme as a !
!        separated program 'progcld2' in the cloud module.             !
!      apr 2005,   yu-tai hou        - modified cloud array and module !
!        structures.                                                   !
!      dec 2008,   yu-tai hou        - changed low-cld calculation,    !
!        now cantains clds from sfc layer and upward to the low/mid    !
!        boundary (include bl-cld). h,m,l clds domain boundaries are   !
!        adjusted for better agreement with observations.              !
!      jan 2011,   yu-tai hou        - changed virtual temperature     !
!        as input variable instead of originally computed inside the   !
!        two prognostic cld schemes 'progcld1' and 'progcld2'.         !
!      aug 2012,   yu-tai hou        - modified subroutine cld_init    !
!        to pass all fixed control variables at the start. and set     !
!        their correponding internal module variables to be used by    !
!        module subroutines.                                           !
!      nov 2012,   yu-tai hou        - modified control parameters     !
!        thru module 'physpara'.                                       !
!      apr 2013,   h.lin/y.hou       - corrected error in calculating  !
!        llyr for different vertical indexing directions.              !
!      jul 2013, r. sun/h. pan       - modified to use pdf cloud and   !
! 	 convective cloud cover and water for radiation                !
!                                                                      !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!


!========================================!
      module module_radiation_clouds     !
!........................................!
!
      use physpara,            only : icldflg, icmphys, iovrsw, iovrlw, &
     &                                lcrick, lcnorm, lnoprec,          &
     &                                ivflip, kind_phys, kind_io4
      use physcons,            only : con_fvirt, con_ttp, con_rocp,     &
     &                                con_t0c, con_pi, con_g, con_rd,   & 
     &                                con_thgni
      use module_microphysics, only : rsipath2
      use module_iounitdef,    only : nicltun
!
      implicit   none
!
      private

!  ---  version tag and last revision date
      character(40), parameter ::                                       &
     &   vtagcld='ncep-radiation_clouds    v5.1  nov 2012 '
!    &   vtagcld='ncep-radiation_clouds    v5.0  aug 2012 '

!  ---  set constant parameters
      real (kind=kind_phys), parameter :: gfac=1.0e5/con_g              &
     &,                                   gord=con_g/con_rd
!      integer, parameter, public :: nf_clds = 11  ! number of fields in cloud array
      integer, parameter, public :: nf_clds = 9   ! number of fields in cloud array
!#endif
      integer, parameter, public :: nk_clds = 3   ! number of cloud vertical domains

!  ---  pressure limits of cloud domain interfaces (low,mid,high) in mb (0.1kpa)
      real (kind=kind_phys), save :: ptopc(nk_clds+1,2)

!org  data ptopc / 1050., 642., 350., 0.0,  1050., 750., 500., 0.0 /
      data ptopc / 1050., 650., 400., 0.0,  1050., 750., 500., 0.0 /

!     real (kind=kind_phys), parameter :: climit = 0.01
      real (kind=kind_phys), parameter :: climit = 0.001, climit2=0.05
      real (kind=kind_phys), parameter :: ovcst  = 1.0 - 1.0e-8

!  ---  set default quantities as parameters (for prognostic cloud)
      real (kind=kind_phys), parameter :: reliq_def = 10.0    ! default liq radius to 10 micron
      real (kind=kind_phys), parameter :: reice_def = 50.0    ! default ice radius to 50 micron
      real (kind=kind_phys), parameter :: rrain_def = 1000.0  ! default rain radius to 1000 micron
      real (kind=kind_phys), parameter :: rsnow_def = 250.0   ! default snow radius to 250 micron

!  ---  set look-up table dimensions and other parameters (for diagnostic cloud)
      integer, parameter :: nbin=100     ! rh in one percent interval
      integer, parameter :: nlon=2       ! =1,2 for eastern and western hemispheres
      integer, parameter :: nlat=4       ! =1,4 for 60n-30n,30n-equ,equ-30s,30s-60s
      integer, parameter :: mcld=4       ! =1,4 for bl,low,mid,hi cld type
      integer, parameter :: nseal=2      ! =1,2 for land,sea

      real (kind=kind_phys), parameter :: cldssa_def = 0.99   ! default cld single scat albedo
      real (kind=kind_phys), parameter :: cldasy_def = 0.84   ! default cld asymmetry factor

!  ---  xlabdy: lat bndry between tuning regions, +/- xlim for transition
!       xlobdy: lon bndry between tuning regions
      real (kind=kind_phys), parameter :: xlim=5.0
      real (kind=kind_phys)            :: xlabdy(3), xlobdy(3)

      data xlabdy / 30.0,  0.0, -30.0 /,  xlobdy / 0.0, 180., 360. /

!  ---  low cloud vertical velocity adjustment boundaries in mb/sec
      real (kind=kind_phys), parameter :: vvcld1= 0.0003e0
      real (kind=kind_phys), parameter :: vvcld2=-0.0005e0

!  ---  those data will be set up by "cld_init"
!     rhcl : tuned rh relation table for diagnostic cloud scheme

      real (kind=kind_phys) :: rhcl(nbin,nlon,nlat,mcld,nseal)

      integer  :: llyr   = 2           ! upper limit of boundary layer clouds
      integer  :: iovr   = 1           ! cloud over lapping method for diagnostic 3-domain
                                       ! output calc (see iovrsw/iovrlw description)

      ! Default ice crystal sizes vs. temperature following Kristjansson and Mitchell
      real (kind=kind_phys), dimension(95), parameter :: retab =(/      &
     &                   5.92779, 6.26422, 6.61973, 6.99539, 7.39234,   &
     &           7.81177, 8.25496, 8.72323, 9.21800, 9.74075, 10.2930,  &
     &           10.8765, 11.4929, 12.1440, 12.8317, 13.5581, 14.2319,  &
     &           15.0351, 15.8799, 16.7674, 17.6986, 18.6744, 19.6955,  &
     &           20.7623, 21.8757, 23.0364, 24.2452, 25.5034, 26.8125,  &
     &           27.7895, 28.6450, 29.4167, 30.1088, 30.7306, 31.2943,  &
     &           31.8151, 32.3077, 32.7870, 33.2657, 33.7540, 34.2601,  &
     &           34.7892, 35.3442, 35.9255, 36.5316, 37.1602, 37.8078,  &
     &           38.4720, 39.1508, 39.8442, 40.5552, 41.2912, 42.0635,  &
     &           42.8876, 43.7863, 44.7853, 45.9170, 47.2165, 48.7221,  &
     &           50.4710, 52.4980, 54.8315, 57.4898, 60.4785, 63.7898,  &
     &           65.5604, 71.2885, 75.4113, 79.7368, 84.2351, 88.8833,  &
     &           93.6658, 98.5739, 103.603, 108.752, 114.025, 119.424,  &
     &           124.954, 130.630, 136.457, 142.446, 148.608, 154.956,  &
     &           161.503, 168.262, 175.248, 182.473, 189.952, 197.699,  &
     &           205.728, 214.055, 222.694, 231.661, 240.971, 250.639/)

      public progcld1, progcld2, progcld3, progcld4, diagcld1, cld_init,&
             progcld5, progcld5o, progclduni, progcld6,                 &
             progcld_thompson, progcld_gce


! =================
      contains
! =================


!-----------------------------------
      subroutine cld_init                                               &
!...................................

!  ---  inputs:
     &     ( si, nlay, me, myrank )
!  ---  outputs:
!          ( none )

!  ===================================================================  !
!                                                                       !
! abstract: cld_init is an initialization program for cloud-radiation   !
!   calculations. it sets up boundary layer cloud top.                  !
!                                                                       !
!                                                                       !
! inputs:                                                               !
!   si (l+1)        : model vertical sigma layer interface              !
!   nlay            : vertical layer number                             !
!   me              : print control flag                                !
!                                                                       !
!  outputs: (none)                                                      !
!           to module variables                                         !
!                                                                       !
!  external module variables: (in physpara)                             !
!   icldflg         : cloud optical property scheme control flag        !
!                     =0: model use diagnostic cloud method             !
!                     =1: model use prognostic cloud method             !
!   icmphys         : cloud microphysics scheme control flag            !
!                     =1: zhao/carr/sundqvist microphysics cloud        !
!                     =2: ferrier microphysics cloud scheme             !
!		      =3: zhao/carr/sundqvist microphysics cloud +pdfcld!
!   iovrsw/iovrlw   : sw/lw control flag for cloud overlapping scheme   !
!                     =0: random overlapping clouds                     !
!                     =1: max/ran overlapping clouds                    !
!   ivflip          : control flag for direction of vertical index      !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!  usage:       call cld_init                                           !
!                                                                       !
!  subroutines called:    rhtable                                       !
!                                                                       !
!  ===================================================================  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: nlay, me, myrank

      real (kind=kind_phys), intent(in) :: si(:)

!  ---  outputs: (none)

!  ---  locals:
      integer :: k, kl, ier

!
!===> ...  begin here
!
!  ---  set up module variables

      iovr    = max( iovrsw, iovrlw )    !cld ovlp used for diag hml cld output

      if (me == 0 .and. myrank ==0) print *, vtagcld      !print out version tag

      if ( icldflg == 0 ) then
        if (me == 0 .and. myrank ==0) print *,' - using diagnostic cloud method'

!  ---  set up tuned rh table

        call rhtable( me, ier )

        if (ier < 0) then
          write(6,99) ier
  99      format(3x,' *** error in finding tuned rh table ***'          &
     &,         /3x,'     stop at calling subroutine rhtable !!'/)
          stop 99
        endif
      else
        if (myrank ==0) then
          print *,' - using prognostic cloud method'
          if (icmphys == 1) then
            print *,'   --- zhao/carr/sundqvist microphysics'
          elseif (icmphys == 2) then
            print *,'   --- ferrier cloud microphysics'
          elseif (icmphys == 3) then
            print *,'   --- zhao/carr/sundqvist + pdf cloud'
          elseif (icmphys == 6) then
            print *,'   --- WSM6 microphysics'
          elseif (icmphys == 8) then
            print *,'   --- Thompson microphysics'
          elseif (icmphys == 18) then
            print *,'   --- 2M Thompson microphysics'
          elseif (icmphys == 11) then
            print *,'   --- GFDL microphysics version 1'
          elseif (icmphys == 12) then
            print *,'   --- GFDL microphysics version 2'
          elseif (icmphys == 13) then
            print *,'   --- GFDL microphysics version 3'
          elseif (icmphys == 15) then
            print *,'   --- Goddard (GCE) 3ICE microphysics'
          elseif (icmphys == 16) then
            print *,'   --- Goddard (GCE) 4ICE microphysics'
          else
            print *,'  !!! error in cloud microphysc specification!!!', &
     &              '  icmphys (np3d) =',icmphys
            stop
          endif
        endif
      endif

!  ---  compute llyr - the top of bl cld and is topmost non cld(low) layer
!       for stratiform (at or above lowest 0.1 of the atmosphere)

      if ( ivflip == 0 ) then    ! data from toa to sfc
        lab_do_k0 : do k = nlay, 2, -1
          kl = k
          if (si(k) < 0.9e0) exit lab_do_k0
        enddo  lab_do_k0

        llyr = kl
      else                      ! data from sfc to top
        lab_do_k1 : do k = 2, nlay
          kl = k
          if (si(k) < 0.9e0) exit lab_do_k1
        enddo  lab_do_k1

        llyr = kl - 1
      endif                     ! end_if_ivflip

!
      return
!...................................
      end subroutine cld_init
!-----------------------------------


!-----------------------------------
      subroutine progcld1                                               &
!...................................

!  ---  inputs:
     &     ( plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,                    &
     &       xlat,xlon,slmsk,                                           &
     &       ix, nlay, nlp1, myrank,                                    &
!  ---  outputs:
     &       clouds,clds,mtop,mbot                                      &
     &      )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:    progcld1    computes cloud related quantities using    !
!   zhao/moorthi's prognostic cloud microphysics scheme.                !
!                                                                       !
! abstract:  this program computes cloud fractions from cloud           !
!   condensates, calculates liquid/ice cloud droplet effective radius,  !
!   and computes the low, mid, high, total and boundary layer cloud     !
!   fractions and the vertical indices of low, mid, and high cloud      !
!   top and base.  the three vertical cloud domains are set up in the   !
!   initial subroutine "cld_init".                                      !
!                                                                       !
! usage:         call progcld1                                          !
!                                                                       !
! subprograms called:   gethml                                          !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (ix,nlay) : model layer mean pressure in mb (100pa)           !
!   plvl  (ix,nlp1) : model level pressure in mb (100pa)                !
!   tlyr  (ix,nlay) : model layer mean temperature in k                 !
!   tvly  (ix,nlay) : model layer virtual temperature in k              !
!   qlyr  (ix,nlay) : layer specific humidity in gm/gm                  !
!   qstl  (ix,nlay) : layer saturate humidity in gm/gm                  !
!   rhly  (ix,nlay) : layer relative humidity (=qlyr/qstl)              !
!   clw   (ix,nlay) : layer cloud condensate amount                     !
!   xlat  (ix)      : grid latitude in radians, default to pi/2 -> -pi/2!
!                     range, otherwise see in-line comment              !
!   xlon  (ix)      : grid longitude in radians  (not used)             !
!   slmsk (ix)      : sea/land mask array (sea:0,land:1,sea-ice:2)      !
!   ix              : horizontal dimention                              !
!   nlay,nlp1       : vertical layer/level dimensions                   !
!                                                                       !
! output variables:                                                     !
!   clouds(ix,nlay,nf_clds) : cloud profiles                            !
!      clouds(:,:,1) - layer total cloud fraction                       !
!      clouds(:,:,2) - layer cloud liq water path         (g/m**2)      !
!      clouds(:,:,3) - mean eff radius for liq cloud      (micron)      !
!      clouds(:,:,4) - layer cloud ice water path         (g/m**2)      !
!      clouds(:,:,5) - mean eff radius for ice cloud      (micron)      !
!      clouds(:,:,6) - layer rain drop water path         not assigned  !
!      clouds(:,:,7) - mean eff radius for rain drop      (micron)      !
!  *** clouds(:,:,8) - layer snow flake water path        not assigned  !
!      clouds(:,:,9) - mean eff radius for snow flake     (micron)      !
!  *** fu's scheme need to be normalized by snow density (g/m**3/1.0e6) !
!   clds  (ix,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (ix,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (ix,3)    : vertical indices for low, mid, hi cloud bases     !
!                                                                       !
! module variables:                                                     !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   lcrick          : control flag for eliminating crick                !
!                     =t: apply layer smoothing to eliminate crick      !
!                     =f: do not apply layer smoothing                  !
!   lcnorm          : control flag for in-cld condensate                !
!                     =t: normalize cloud condensate                    !
!                     =f: not normalize cloud condensate                !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs
      integer,  intent(in) :: ix, nlay, nlp1, myrank

      real (kind=kind_phys), dimension(:,:), intent(in) :: plvl, plyr,  &
     &       tlyr, tvly, qlyr, qstl, rhly, clw

      real (kind=kind_phys), dimension(:),   intent(in) :: xlat, xlon,  &
     &       slmsk

!  ---  outputs
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: clouds

      real (kind=kind_phys), dimension(:,:),   intent(out) :: clds

      integer,               dimension(:,:),   intent(out) :: mtop,mbot

!  ---  local variables:
      real (kind=kind_phys), dimension(ix,nlay) :: cldtot, cldcnv,      &
     &       cwp, cip, crp, csp, rew, rei, res, rer, delp, tem2d, clwf

      real (kind=kind_phys) :: ptop1(ix,nk_clds+1)

      real (kind=kind_phys) :: clwmin, clwm, clwt, onemrh, value,       &
     &       tem1, tem2, tem3

      integer :: i, k, id, nf
!
!===> ... begin here
!
      do nf=1,nf_clds
        do k=1,nlay
          do i=1,ix
            clouds(i,k,nf) = 0.0
          enddo
        enddo
      enddo
!     clouds(:,:,:) = 0.0

      do k = 1, nlay
        do i = 1, ix
          cldtot(i,k) = 0.0
          cldcnv(i,k) = 0.0
          cwp   (i,k) = 0.0
          cip   (i,k) = 0.0
          crp   (i,k) = 0.0
          csp   (i,k) = 0.0
          rew   (i,k) = reliq_def            ! default liq radius to 10 micron
          rei   (i,k) = reice_def            ! default ice radius to 50 micron
          rer   (i,k) = rrain_def            ! default rain radius to 1000 micron
          res   (i,k) = rsnow_def            ! default snow radius to 250 micron
          tem2d (i,k) = min( 1.0, max( 0.0, (con_ttp-tlyr(i,k))*0.05 ) )
          clwf(i,k)   = 0.0
        enddo
      enddo
!
      if ( lcrick ) then
        do i = 1, ix
          clwf(i,1)    = 0.75*clw(i,1)    + 0.25*clw(i,2)
          clwf(i,nlay) = 0.75*clw(i,nlay) + 0.25*clw(i,nlay-1)
        enddo
        do k = 2, nlay-1
          do i = 1, ix
            clwf(i,k) = 0.25*clw(i,k-1) + 0.5*clw(i,k) + 0.25*clw(i,k+1)
          enddo
        enddo
      else
        do k = 1, nlay
          do i = 1, ix
            clwf(i,k) = clw(i,k)
          enddo
        enddo
      endif

!  ---  find top pressure for each cloud domain for given latitude
!       ptopc(k,i): top presure of each cld domain (k=1-4 are sfc,l,m,h;
!  ---  i=1,2 are low-lat (<45 degree) and pole regions)

      do id = 1, 4
        tem1 = ptopc(id,2) - ptopc(id,1)

        do i =1, ix
          tem2 = xlat(i) / con_pi        ! if xlat in pi/2 -> -pi/2 range
!         tem2 = 0.5 - xlat(i)/con_pi    ! if xlat in 0 -> pi range

          ptop1(i,id) = ptopc(id,1) + tem1*max( 0.0, 4.0*abs(tem2)-1.0 )
        enddo
      enddo

!  ---  compute liquid/ice condensate path in g/m**2

      if ( ivflip == 0 ) then          ! input data from toa to sfc
        do k = 1, nlay
          do i = 1, ix
            delp(i,k) = plvl(i,k+1) - plvl(i,k)
            clwt     = max(0.0, clwf(i,k)) * gfac * delp(i,k)
            cip(i,k) = clwt * tem2d(i,k)
            cwp(i,k) = clwt - cip(i,k)
          enddo
        enddo
      else                             ! input data from sfc to toa
        do k = 1, nlay
          do i = 1, ix
            delp(i,k) = plvl(i,k) - plvl(i,k+1)
            clwt     = max(0.0, clwf(i,k)) * gfac * delp(i,k)
            cip(i,k) = clwt * tem2d(i,k)
            cwp(i,k) = clwt - cip(i,k)
          enddo
        enddo
      endif                            ! end_if_ivflip

!  ---  effective liquid cloud droplet radius over land

      do i = 1, ix
        if (nint(slmsk(i)) == 1) then
          do k = 1, nlay
            rew(i,k) = 5.0 + 5.0 * tem2d(i,k)
          enddo
        endif
      enddo

!  ---  layer cloud fraction

      if ( ivflip == 0 ) then              ! input data from toa to sfc

        clwmin = 0.0
        do k = nlay, 1, -1
        do i = 1, ix
          clwt = 1.0e-6 * (plyr(i,k)*0.001)
!         clwt = 2.0e-6 * (plyr(i,k)*0.001)

          if (clwf(i,k) > clwt) then
            onemrh= max( 1.e-10, 1.0-rhly(i,k) )
            clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )

!           tem1  = min(max(sqrt(sqrt(onemrh*qstl(i,k))),0.0001),1.0)
!           tem1  = 2000.0 / tem1

            tem1  = min(max((onemrh*qstl(i,k))**0.49,0.0001),1.0)  !jhan
            tem1  = 100.0 / tem1
!
!           tem1  = 2000.0 / tem1
!           tem1  = 1000.0 / tem1
!

            value = max( min( tem1*(clwf(i,k)-clwm), 50.0 ), 0.0 )
            tem2  = sqrt( sqrt(rhly(i,k)) )
            cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
          endif
         enddo
         enddo

      else                                 ! input data from sfc to toa

        clwmin = 0.0
        do k = 1, nlay
        do i = 1, ix
          clwt = 1.0e-6 * (plyr(i,k)*0.001)
!         clwt = 2.0e-6 * (plyr(i,k)*0.001)

          if (clwf(i,k) > clwt) then
            onemrh= max( 1.e-10, 1.0-rhly(i,k) )
            clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )

!           tem1  = min(max(sqrt(sqrt(onemrh*qstl(i,k))),0.0001),1.0)
!           tem1  = 2000.0 / tem1

            tem1  = min(max((onemrh*qstl(i,k))**0.49,0.0001),1.0)  !jhan
            tem1  = 100.0 / tem1
!
!           tem1  = 2000.0 / tem1
!           tem1  = 1000.0 / tem1
!

            value = max( min( tem1*(clwf(i,k)-clwm), 50.0 ), 0.0 )
            tem2  = sqrt( sqrt(rhly(i,k)) )

            cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
          endif
        enddo
        enddo

      endif                                ! end_if_flip

      do k = 1, nlay
        do i = 1, ix
          if (cldtot(i,k) < climit) then
            cldtot(i,k) = 0.0
            cwp(i,k)    = 0.0
            cip(i,k)    = 0.0
            crp(i,k)    = 0.0
            csp(i,k)    = 0.0
          endif
        enddo
      enddo

      if ( lcnorm ) then
        do k = 1, nlay
          do i = 1, ix
            if (cldtot(i,k) >= climit) then
              tem1 = 1.0 / max(climit2, cldtot(i,k))
              cwp(i,k) = cwp(i,k) * tem1
              cip(i,k) = cip(i,k) * tem1
              crp(i,k) = crp(i,k) * tem1
              csp(i,k) = csp(i,k) * tem1
            endif
          enddo
        enddo
      endif

!  ---  effective ice cloud droplet radius

      do k = 1, nlay
        do i = 1, ix
          tem2 = tlyr(i,k) - con_ttp

          if (cip(i,k) > 0.0) then
!           tem3 = gord * cip(i,k) * (plyr(i,k)/delp(i,k)) / tvly(i,k)
            tem3 = gord * cip(i,k) * plyr(i,k) / (delp(i,k)*tvly(i,k))

            if (tem2 < -50.0) then
              rei(i,k) = (1250.0/9.917) * tem3 ** 0.109
            elseif (tem2 < -40.0) then
              rei(i,k) = (1250.0/9.337) * tem3 ** 0.08
            elseif (tem2 < -30.0) then
              rei(i,k) = (1250.0/9.208) * tem3 ** 0.055
            else
              rei(i,k) = (1250.0/9.387) * tem3 ** 0.031
            endif
!           rei(i,k)   = max(20.0, min(rei(i,k), 300.0))
!           rei(i,k)   = max(10.0, min(rei(i,k), 100.0))
            rei(i,k)   = max(10.0, min(rei(i,k), 150.0))
!           rei(i,k)   = max(5.0, min(rei(i,k), 130.0))
          endif
        enddo
      enddo

!
      do k = 1, nlay
        do i = 1, ix
          clouds(i,k,1) = cldtot(i,k)
          clouds(i,k,2) = cwp(i,k)
          clouds(i,k,3) = rew(i,k)
          clouds(i,k,4) = cip(i,k)
          clouds(i,k,5) = rei(i,k)
!         clouds(i,k,6) = 0.0
          clouds(i,k,7) = rer(i,k)
!         clouds(i,k,8) = 0.0
!          clouds(i,k,9) = rei(i,k)
          clouds(i,k,9) = res(i,k)
        enddo
      enddo


!  ---  compute low, mid, high, total, and boundary layer cloud fractions
!       and clouds top/bottom layer indices for low, mid, and high clouds.
!       the three cloud domain boundaries are defined by ptopc.  the cloud
!       overlapping method is defined by control flag 'iovr', which may
!       be different for lw and sw radiation programs.

      call gethml                                                       &
!  ---  inputs:
     &     ( plyr, ptop1, cldtot, cldcnv,                               &
     &       ix,nlay,                                                   &
!  ---  outputs:
     &       clds, mtop, mbot                                           &
     &     )


!
      return
!...................................
      end subroutine progcld1
!-----------------------------------

!-----------------------------------
      subroutine progcld3                                               &
!...................................

!  ---  inputs:
     &     ( plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,cnvw,cnvc,          &
     &       xlat,xlon,slmsk,                                           &
     &       ix, nlay, nlp1,                                            &
     &       deltaq,sup,kdt,me,                                         &
!  ---  outputs:
     &       clouds,clds,mtop,mbot                                      &
     &      )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:    progcld3    computes cloud related quantities using    !
!   zhao/moorthi's prognostic cloud microphysics scheme.                !
!                                                                       !
! abstract:  this program computes cloud fractions from cloud           !
!   condensates, calculates liquid/ice cloud droplet effective radius,  !
!   and computes the low, mid, high, total and boundary layer cloud     !
!   fractions and the vertical indices of low, mid, and high cloud      !
!   top and base.  the three vertical cloud domains are set up in the   !
!   initial subroutine "cld_init".                                      !
!                                                                       !
! usage:         call progcld3                                          !
!                                                                       !
! subprograms called:   gethml                                          !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (ix,nlay) : model layer mean pressure in mb (100pa)           !
!   plvl  (ix,nlp1) : model level pressure in mb (100pa)                !
!   tlyr  (ix,nlay) : model layer mean temperature in k                 !
!   tvly  (ix,nlay) : model layer virtual temperature in k              !
!   qlyr  (ix,nlay) : layer specific humidity in gm/gm                  !
!   qstl  (ix,nlay) : layer saturate humidity in gm/gm                  !
!   rhly  (ix,nlay) : layer relative humidity (=qlyr/qstl)              !
!   clw   (ix,nlay) : layer cloud condensate amount                     !
!   xlat  (ix)      : grid latitude in radians, default to pi/2 -> -pi/2!
!                     range, otherwise see in-line comment              !
!   xlon  (ix)      : grid longitude in radians  (not used)             !
!   slmsk (ix)      : sea/land mask array (sea:0,land:1,sea-ice:2)      !
!   ix              : horizontal dimention                              !
!   nlay,nlp1       : vertical layer/level dimensions                   !
!   cnvw  (ix,nlay) : layer convective cloud condensate                 !
!   cnvc  (ix,nlay) : layer convective cloud cover                      !
!   deltaq(ix,nlay) : half total water distribution width               !
!   sup             : supersaturation                                   ! 

!                                                                       !
! output variables:                                                     !
!   clouds(ix,nlay,nf_clds) : cloud profiles                            !
!      clouds(:,:,1) - layer total cloud fraction                       !
!      clouds(:,:,2) - layer cloud liq water path         (g/m**2)      !
!      clouds(:,:,3) - mean eff radius for liq cloud      (micron)      !
!      clouds(:,:,4) - layer cloud ice water path         (g/m**2)      !
!      clouds(:,:,5) - mean eff radius for ice cloud      (micron)      !
!      clouds(:,:,6) - layer rain drop water path         not assigned  !
!      clouds(:,:,7) - mean eff radius for rain drop      (micron)      !
!  *** clouds(:,:,8) - layer snow flake water path        not assigned  !
!      clouds(:,:,9) - mean eff radius for snow flake     (micron)      !
!  *** fu's scheme need to be normalized by snow density (g/m**3/1.0e6) !
!   clds  (ix,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (ix,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (ix,3)    : vertical indices for low, mid, hi cloud bases     !
!                                                                       !
! module variables:                                                     !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   lcrick          : control flag for eliminating crick                !
!                     =t: apply layer smoothing to eliminate crick      !
!                     =f: do not apply layer smoothing                  !
!   lcnorm          : control flag for in-cld condensate                !
!                     =t: normalize cloud condensate                    !
!                     =f: not normalize cloud condensate                !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs
      integer,  intent(in) :: ix, nlay, nlp1,kdt

      real (kind=kind_phys), dimension(:,:), intent(in) :: plvl, plyr,  &
     &       tlyr, tvly, qlyr, qstl, rhly, clw
!     &       tlyr, tvly, qlyr, qstl, rhly, clw, cnvw, cnvc
!      real (kind=kind_phys), dimension(:,:), intent(in) :: deltaq
      real (kind=kind_phys), dimension(:,:) :: deltaq, cnvw, cnvc
      real (kind=kind_phys) qtmp,qsc,rhs
      real (kind=kind_phys), intent(in) :: sup
      real (kind=kind_phys), parameter :: epsq = 1.0e-12

      real (kind=kind_phys), dimension(:),   intent(in) :: xlat, xlon,  &
     &       slmsk
      integer :: me  

!  ---  outputs
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: clouds

      real (kind=kind_phys), dimension(:,:),   intent(out) :: clds

      integer,               dimension(:,:),   intent(out) :: mtop,mbot

!  ---  local variables:
      real (kind=kind_phys), dimension(ix,nlay) :: cldtot, cldcnv,      &
     &       cwp, cip, crp, csp, rew, rei, res, rer, delp, tem2d, clwf

      real (kind=kind_phys) :: ptop1(ix,nk_clds+1)

      real (kind=kind_phys) :: clwmin, clwm, clwt, onemrh, value,       &
     &       tem1, tem2, tem3

      integer :: i, k, id, nf

!
!===> ... begin here
!
      do nf=1,nf_clds
        do k=1,nlay
          do i=1,ix
            clouds(i,k,nf) = 0.0
          enddo
        enddo
      enddo
!     clouds(:,:,:) = 0.0

      do k = 1, nlay
        do i = 1, ix
          cldtot(i,k) = 0.0
          cldcnv(i,k) = 0.0
          cwp   (i,k) = 0.0
          cip   (i,k) = 0.0
          crp   (i,k) = 0.0
          csp   (i,k) = 0.0
          rew   (i,k) = reliq_def            ! default liq radius to 10 micron
          rei   (i,k) = reice_def            ! default ice radius to 50 micron
          rer   (i,k) = rrain_def            ! default rain radius to 1000 micron
          res   (i,k) = rsnow_def            ! default snow radius to 250 micron
          tem2d (i,k) = min( 1.0, max( 0.0, (con_ttp-tlyr(i,k))*0.05 ) )
          clwf(i,k)   = 0.0
        enddo
      enddo
!
      if ( lcrick ) then
        do i = 1, ix
          clwf(i,1)    = 0.75*clw(i,1)    + 0.25*clw(i,2)
          clwf(i,nlay) = 0.75*clw(i,nlay) + 0.25*clw(i,nlay-1)
        enddo
        do k = 2, nlay-1
          do i = 1, ix
            clwf(i,k) = 0.25*clw(i,k-1) + 0.5*clw(i,k) + 0.25*clw(i,k+1)
          enddo
        enddo
      else
        do k = 1, nlay
          do i = 1, ix
            clwf(i,k) = clw(i,k)
          enddo
        enddo
      endif

      if(kdt==1) then
        do k = 1, nlay
          do i = 1, ix
            deltaq(i,k) = (1.-0.95)*qstl(i,k)
          enddo
        enddo
      endif

!  ---  find top pressure for each cloud domain for given latitude
!       ptopc(k,i): top presure of each cld domain (k=1-4 are sfc,l,m,h;
!  ---  i=1,2 are low-lat (<45 degree) and pole regions)

      do id = 1, 4
        tem1 = ptopc(id,2) - ptopc(id,1)

        do i =1, ix
          tem2 = xlat(i) / con_pi        ! if xlat in pi/2 -> -pi/2 range
!         tem2 = 0.5 - xlat(i)/con_pi    ! if xlat in 0 -> pi range

          ptop1(i,id) = ptopc(id,1) + tem1*max( 0.0, 4.0*abs(tem2)-1.0 )
        enddo
      enddo

!  ---  compute liquid/ice condensate path in g/m**2

      if ( ivflip == 0 ) then          ! input data from toa to sfc
        do k = 1, nlay
          do i = 1, ix
            delp(i,k) = plvl(i,k+1) - plvl(i,k)
            clwt  = max(0.0,(clwf(i,k)+cnvw(i,k))) * gfac * delp(i,k)
            cip(i,k) = clwt * tem2d(i,k)
            cwp(i,k) = clwt - cip(i,k)
          enddo
        enddo
      else                             ! input data from sfc to toa
        do k = 1, nlay
          do i = 1, ix
            delp(i,k) = plvl(i,k) - plvl(i,k+1)
            clwt = max(0.0,(clwf(i,k)+cnvw(i,k))) * gfac * delp(i,k)
            cip(i,k) = clwt * tem2d(i,k)
            cwp(i,k) = clwt - cip(i,k)
          enddo
        enddo
      endif                            ! end_if_ivflip

!  ---  effective liquid cloud droplet radius over land

      do i = 1, ix
        if (nint(slmsk(i)) == 1) then
          do k = 1, nlay
            rew(i,k) = 5.0 + 5.0 * tem2d(i,k)
          enddo
        endif
      enddo

!  ---  layer cloud fraction

      if ( ivflip == 0 ) then              ! input data from toa to sfc
          do k = nlay, 1, -1
          do i = 1, ix
            tem1 = tlyr(i,k) - 273.16
            if(tem1 < con_thgni) then  ! for pure ice, has to be consistent with gscond
              qsc = sup * qstl(i,k)
              rhs = sup
            else
              qsc = qstl(i,k)
              rhs = 1.0
            endif
           if(rhly(i,k) >= rhs) then
             cldtot(i,k) = 1.0
           else
             qtmp = qlyr(i,k) + clwf(i,k) - qsc
             if(deltaq(i,k) > epsq) then
               if(qtmp < -deltaq(i,k) .or. clwf(i,k) < epsq) then
!               if(qtmp < -deltaq(i,k)) then
                 cldtot(i,k) = 0.0
               elseif(qtmp >= deltaq(i,k)) then
                 cldtot(i,k) = 1.0
               else
                 cldtot(i,k) = 0.5*qtmp/deltaq(i,k) + 0.5
                 cldtot(i,k) = max(cldtot(i,k),0.0)
                 cldtot(i,k) = min(cldtot(i,k),1.0)
               endif
             else
               if(qtmp.gt.0) then
                 cldtot(i,k) = 1.0
               else
                 cldtot(i,k) = 0.0
               endif
             endif
           endif
           cldtot(i,k) = cnvc(i,k)+(1-cnvc(i,k))*cldtot(i,k)
           cldtot(i,k) = max(cldtot(i,k),0.)
           cldtot(i,k) = min(cldtot(i,k),1.)
          enddo
          enddo
      else                                 ! input data from sfc to toa
          do k = 1, nlay
          do i = 1, ix
            tem1 = tlyr(i,k) - 273.16
            if(tem1 < con_thgni) then  ! for pure ice, has to be consistent with gscond
              qsc = sup * qstl(i,k)
              rhs = sup
            else
              qsc = qstl(i,k)
              rhs = 1.0
            endif
           if(rhly(i,k) >= rhs) then
             cldtot(i,k) = 1.0
           else
             qtmp = qlyr(i,k) + clwf(i,k) - qsc
             if(deltaq(i,k) > epsq) then
!              if(qtmp <= -deltaq(i,k) .or. cwmik < epsq) then
               if(qtmp <= -deltaq(i,k)) then
                 cldtot(i,k) = 0.0
               elseif(qtmp >= deltaq(i,k)) then
                 cldtot(i,k) = 1.0
               else
                 cldtot(i,k) = 0.5*qtmp/deltaq(i,k) + 0.5
                 cldtot(i,k) = max(cldtot(i,k),0.0)
                 cldtot(i,k) = min(cldtot(i,k),1.0)
               endif
             else
               if(qtmp > 0.) then
                 cldtot(i,k) = 1.0
               else
                 cldtot(i,k) = 0.0
               endif
             endif
           endif
           cldtot(i,k) = cnvc(i,k) + (1-cnvc(i,k))*cldtot(i,k)
           cldtot(i,k) = max(cldtot(i,k),0.)
           cldtot(i,k) = min(cldtot(i,k),1.)

          enddo
          enddo
      endif                                ! end_if_flip

      do k = 1, nlay
        do i = 1, ix
          if (cldtot(i,k) < climit) then
            cldtot(i,k) = 0.0
            cwp(i,k)    = 0.0
            cip(i,k)    = 0.0
            crp(i,k)    = 0.0
            csp(i,k)    = 0.0
          endif
        enddo
      enddo

      if ( lcnorm ) then
        do k = 1, nlay
          do i = 1, ix
            if (cldtot(i,k) >= climit) then
              tem1 = 1.0 / max(climit2, cldtot(i,k))
              cwp(i,k) = cwp(i,k) * tem1
              cip(i,k) = cip(i,k) * tem1
              crp(i,k) = crp(i,k) * tem1
              csp(i,k) = csp(i,k) * tem1
            endif
          enddo
        enddo
      endif

!  ---  effective ice cloud droplet radius

      do k = 1, nlay
        do i = 1, ix
          tem2 = tlyr(i,k) - con_ttp

          if (cip(i,k) > 0.0) then
!           tem3 = gord * cip(i,k) * (plyr(i,k)/delp(i,k)) / tvly(i,k)
            tem3 = gord * cip(i,k) * plyr(i,k) / (delp(i,k)*tvly(i,k))

            if (tem2 < -50.0) then
              rei(i,k) = (1250.0/9.917) * tem3 ** 0.109
            elseif (tem2 < -40.0) then
              rei(i,k) = (1250.0/9.337) * tem3 ** 0.08
            elseif (tem2 < -30.0) then
              rei(i,k) = (1250.0/9.208) * tem3 ** 0.055
            else
              rei(i,k) = (1250.0/9.387) * tem3 ** 0.031
            endif
!           rei(i,k)   = max(20.0, min(rei(i,k), 300.0))
!           rei(i,k)   = max(10.0, min(rei(i,k), 100.0))
            rei(i,k)   = max(10.0, min(rei(i,k), 150.0))
!           rei(i,k)   = max(5.0, min(rei(i,k), 130.0))
          endif
        enddo
      enddo

!
      do k = 1, nlay
        do i = 1, ix
          clouds(i,k,1) = cldtot(i,k)
          clouds(i,k,2) = cwp(i,k)
          clouds(i,k,3) = rew(i,k)
          clouds(i,k,4) = cip(i,k)
          clouds(i,k,5) = rei(i,k)
!         clouds(i,k,6) = 0.0
          clouds(i,k,7) = rer(i,k)
!         clouds(i,k,8) = 0.0
          clouds(i,k,9) = rei(i,k)
        enddo
      enddo


!  ---  compute low, mid, high, total, and boundary layer cloud fractions
!       and clouds top/bottom layer indices for low, mid, and high clouds.
!       the three cloud domain boundaries are defined by ptopc.  the cloud
!       overlapping method is defined by control flag 'iovr', which may
!       be different for lw and sw radiation programs.

      call gethml                                                       &
!  ---  inputs:
     &     ( plyr, ptop1, cldtot, cldcnv,                               &
     &       ix,nlay,                                                   &
!  ---  outputs:
     &       clds, mtop, mbot                                           &
     &     )


!
      return
!...................................
      end subroutine progcld3
!-----------------------------------

!-----------------------------------
      subroutine progcld2                                               &
!...................................

!  ---  inputs:
     &     ( plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,                    &
     &       xlat,xlon,slmsk, f_ice,f_rain,r_rime,flgmin,               &
     &       ix, nlay, nlp1,                                            &
!  ---  outputs:
     &       clouds,clds,mtop,mbot                                      &
     &      )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:    progcld2    computes cloud related quantities using    !
!   ferrier's prognostic cloud microphysics scheme.                     !
!                                                                       !
! abstract:  this program computes cloud fractions from cloud           !
!   condensates, calculates liquid/ice cloud droplet effective radius,  !
!   and computes the low, mid, high, total and boundary layer cloud     !
!   fractions and the vertical indices of low, mid, and high cloud      !
!   top and base.  the three vertical cloud domains are set up in the   !
!   initial subroutine "cld_init".                                      !
!                                                                       !
! usage:         call progcld2                                          !
!                                                                       !
! subprograms called:   gethml                                          !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (ix,nlay) : model layer mean pressure in mb (100pa)           !
!   plvl  (ix,nlp1) : model level pressure in mb (100pa)                !
!   tlyr  (ix,nlay) : model layer mean temperature in k                 !
!   tvly  (ix,nlay) : model layer virtual temperature in k              !
!   qlyr  (ix,nlay) : layer specific humidity in gm/gm                  !
!   qstl  (ix,nlay) : layer saturate humidity in gm/gm                  !
!   rhly  (ix,nlay) : layer relative humidity (=qlyr/qstl)              !
!   clw   (ix,nlay) : layer cloud condensate amount                     !
!   f_ice (ix,nlay) : fraction of layer cloud ice  (ferrier micro-phys) !
!   f_rain(ix,nlay) : fraction of layer rain water (ferrier micro-phys) !
!   r_rime(ix,nlay) : mass ratio of total ice to unrimed ice (>=1)      !
!   flgmin(ix)      : minimim large ice fraction                        !
!   xlat  (ix)      : grid latitude in radians, default to pi/2 -> -pi/2!
!                     range, otherwise see in-line comment              !
!   xlon  (ix)      : grid longitude in radians  (not used)             !
!   slmsk (ix)      : sea/land mask array (sea:0,land:1,sea-ice:2)      !
!   ix              : horizontal dimention                              !
!   nlay,nlp1       : vertical layer/level dimensions                   !
!                                                                       !
! output variables:                                                     !
!   clouds(ix,nlay,nf_clds) : cloud profiles                            !
!      clouds(:,:,1) - layer total cloud fraction                       !
!      clouds(:,:,2) - layer cloud liq water path         (g/m**2)      !
!      clouds(:,:,3) - mean eff radius for liq cloud      (micron)      !
!      clouds(:,:,4) - layer cloud ice water path         (g/m**2)      !
!      clouds(:,:,5) - mean eff radius for ice cloud      (micron)      !
!      clouds(:,:,6) - layer rain drop water path         (g/m**2)      !
!      clouds(:,:,7) - mean eff radius for rain drop      (micron)      !
!  *** clouds(:,:,8) - layer snow flake water path        (g/m**2)      !
!      clouds(:,:,9) - mean eff radius for snow flake     (micron)      !
!  *** fu's scheme need to be normalized by snow density (g/m**3/1.0e6) !
!   clds  (ix,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (ix,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (ix,3)    : vertical indices for low, mid, hi cloud bases     !
!                                                                       !
! external module variables:                                            !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   lcrick          : control flag for eliminating crick                !
!                     =t: apply layer smoothing to eliminate crick      !
!                     =f: do not apply layer smoothing                  !
!   lcnorm          : control flag for in-cld condensate                !
!                     =t: normalize cloud condensate                    !
!                     =f: not normalize cloud condensate                !
!   lnoprec         : precip effect in radiation flag (ferrier scheme)  !
!                     =t: snow/rain has no impact on radiation          !
!                     =f: snow/rain has impact on radiation             !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  constants

!  ---  inputs
      integer,  intent(in) :: ix, nlay, nlp1

      real (kind=kind_phys), dimension(:,:), intent(in) :: plvl, plyr,  &
     &       tlyr, tvly, qlyr, qstl, rhly, clw, f_ice, f_rain, r_rime

      real (kind=kind_phys), dimension(:),   intent(in) :: xlat, xlon,  &
     &       slmsk
      real (kind=kind_phys), dimension(:), intent(in) :: flgmin

!  ---  outputs
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: clouds

      real (kind=kind_phys), dimension(:,:),   intent(out) :: clds

      integer,               dimension(:,:),   intent(out) :: mtop,mbot

!  ---  local variables:
      real (kind=kind_phys), dimension(ix,nlay) :: cldtot, cldcnv,      &
     &       cwp, cip, crp, csp, rew, rei, res, rer, tem2d, clw2,       &
     &       qcwat, qcice, qrain, fcice, frain, rrime, rsden, clwf

      real (kind=kind_phys) :: ptop1(ix,nk_clds+1)

      real (kind=kind_phys) :: clwmin, clwm, clwt, onemrh, value,       &
     &       tem1, tem2, tem3

      integer :: i, k, id

!
!===> ... begin here
!
!     clouds(:,:,:) = 0.0

      do k = 1, nlay
        do i = 1, ix
          cldtot(i,k) = 0.0
          cldcnv(i,k) = 0.0
          cwp   (i,k) = 0.0
          cip   (i,k) = 0.0
          crp   (i,k) = 0.0
          csp   (i,k) = 0.0
          rew   (i,k) = reliq_def            ! default liq radius to 10 micron
          rei   (i,k) = reice_def            ! default ice radius to 50 micron
          rer   (i,k) = rrain_def            ! default rain radius to 1000 micron
          res   (i,k) = rsnow_def            ! default snow radius to 250 micron
          fcice (i,k) = max(0.0, min(1.0, f_ice(i,k)))
          frain (i,k) = max(0.0, min(1.0, f_rain(i,k)))
          rrime (i,k) = max(1.0, r_rime(i,k))
          tem2d (i,k) = tlyr(i,k) - con_t0c
        enddo
      enddo
!
      if ( lcrick ) then
        do i = 1, ix
          clwf(i,1)    = 0.75*clw(i,1)    + 0.25*clw(i,2)
          clwf(i,nlay) = 0.75*clw(i,nlay) + 0.25*clw(i,nlay-1)
        enddo
        do k = 2, nlay-1
          do i = 1, ix
            clwf(i,k) = 0.25*clw(i,k-1) + 0.5*clw(i,k) + 0.25*clw(i,k+1)
          enddo
        enddo
      else
        do k = 1, nlay
          do i = 1, ix
            clwf(i,k) = clw(i,k)
          enddo
        enddo
      endif

!  ---  find top pressure for each cloud domain for given latitude
!       ptopc(k,i): top presure of each cld domain (k=1-4 are sfc,l,m,h;
!  ---  i=1,2 are low-lat (<45 degree) and pole regions)

      do id = 1, 4
        tem1 = ptopc(id,2) - ptopc(id,1)

        do i =1, ix
          tem2 = xlat(i) / con_pi        ! if xlat in pi/2 -> -pi/2 range
!         tem2 = 0.5 - xlat(i)/con_pi    ! if xlat in 0 -> pi range

          ptop1(i,id) = ptopc(id,1) + tem1*max( 0.0, 4.0*abs(tem2)-1.0 )
        enddo
      enddo

!  ---  separate cloud condensate into liquid, ice, and rain types, and
!       save the liquid+ice condensate in array clw2 for later
!       calculation of cloud fraction

      do k = 1, nlay
        do i = 1, ix
          if (tem2d(i,k) > -40.0) then
            qcice(i,k) = clwf(i,k) * fcice(i,k)
            tem1       = clwf(i,k) - qcice(i,k)
            qrain(i,k) = tem1 * frain(i,k)
            qcwat(i,k) = tem1 - qrain(i,k)
            clw2 (i,k) = qcwat(i,k) + qcice(i,k)
          else
            qcice(i,k) = clwf(i,k)
            qrain(i,k) = 0.0
            qcwat(i,k) = 0.0
            clw2 (i,k) = clwf(i,k)
          endif
        enddo
      enddo

      call  rsipath2                                                    &
!  ---  inputs:
     &     ( plyr, plvl, tlyr, qlyr, qcwat, qcice, qrain, rrime,        &
     &       ix, nlay, ivflip, flgmin,                                  &
!  ---  outputs:
     &       cwp, cip, crp, csp, rew, rer, res, rsden                   &
     &     )


      if ( ivflip == 0 ) then          ! input data from toa to sfc
        do k = 1, nlay
          do i = 1, ix
            tem2d(i,k) = (con_g * plyr(i,k))                            &
     &                 / (con_rd* (plvl(i,k+1) - plvl(i,k)))
          enddo
        enddo
      else                             ! input data from sfc to toa
        do k = 1, nlay
          do i = 1, ix
            tem2d(i,k) = (con_g * plyr(i,k))                            &
     &                 / (con_rd* (plvl(i,k) - plvl(i,k+1)))
          enddo
        enddo
      endif                            ! end_if_ivflip

!  ---  layer cloud fraction

      if ( ivflip == 0 ) then              ! input data from toa to sfc

        clwmin = 0.0
        do k = nlay, 1, -1
        do i = 1, ix
!         clwt = 1.0e-6 * (plyr(i,k)*0.001)
          clwt = 2.0e-6 * (plyr(i,k)*0.001)

          if (clw2(i,k) > clwt) then
            onemrh= max( 1.e-10, 1.0-rhly(i,k) )
            clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )

            tem1  = min(max((onemrh*qstl(i,k))**0.49,0.0001),1.0)    !jhan
            tem1  = 100.0 / tem1

!           tem1  = min(max(sqrt(sqrt(onemrh*qstl(i,k))),0.0001),1.0)
!           tem1  = 2000.0 / tem1
!
!           tem1  = min(max(sqrt(sqrt(onemrh*qstl(i,k))),0.0001),1.0)
!           tem1  = 2200.0 / tem1
!           tem1  = 2400.0 / tem1
!           tem1  = 2500.0 / tem1
!           tem1  = min(max(sqrt(onemrh*qstl(i,k)),0.0001),1.0)
!           tem1  = 2000.0 / tem1
!           tem1  = 1000.0 / tem1
!           tem1  = 100.0 / tem1

            value = max( min( tem1*(clw2(i,k)-clwm), 50.0 ), 0.0 )
            tem2  = sqrt( sqrt(rhly(i,k)) )

            cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
          endif
        enddo
        enddo

      else                                 ! input data from sfc to toa

        clwmin = 0.0e-6
        do k = 1, nlay
        do i = 1, ix
!         clwt = 1.0e-6 * (plyr(i,k)*0.001)
          clwt = 2.0e-6 * (plyr(i,k)*0.001)

          if (clw2(i,k) > clwt) then
            onemrh= max( 1.e-10, 1.0-rhly(i,k) )
            clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )

            tem1  = min(max((onemrh*qstl(i,k))**0.49,0.0001),1.0)   !jhan
            tem1  = 100.0 / tem1

!           tem1  = min(max(sqrt(sqrt(onemrh*qstl(i,k))),0.0001),1.0)
!           tem1  = 2000.0 / tem1
!
!           tem1  = min(max(sqrt(sqrt(onemrh*qstl(i,k))),0.0001),1.0)
!           tem1  = 2200.0 / tem1
!           tem1  = 2400.0 / tem1
!           tem1  = 2500.0 / tem1
!           tem1  = min(max(sqrt(onemrh*qstl(i,k)),0.0001),1.0)
!           tem1  = 2000.0 / tem1
!           tem1  = 1000.0 / tem1
!           tem1  = 100.0 / tem1

            value = max( min( tem1*(clw2(i,k)-clwm), 50.0 ), 0.0 )
            tem2  = sqrt( sqrt(rhly(i,k)) )

            cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
          endif
        enddo
        enddo

      endif                                ! end_if_flip

      do k = 1, nlay
        do i = 1, ix
          if (cldtot(i,k) < climit) then
            cldtot(i,k) = 0.0
            cwp(i,k)    = 0.0
            cip(i,k)    = 0.0
            crp(i,k)    = 0.0
            csp(i,k)    = 0.0
          endif
        enddo
      enddo

!     when lnoprec = .true. snow/rain has no impact on radiation
      if ( lnoprec ) then
        do k = 1, nlay
          do i = 1, ix
            crp(i,k) = 0.0
            csp(i,k) = 0.0
          enddo
        enddo
      endif
!
      if ( lcnorm ) then
        do k = 1, nlay
          do i = 1, ix
            if (cldtot(i,k) >= climit) then
              tem1 = 1.0 / max(climit2, cldtot(i,k))
              cwp(i,k) = cwp(i,k) * tem1
              cip(i,k) = cip(i,k) * tem1
              crp(i,k) = crp(i,k) * tem1
              csp(i,k) = csp(i,k) * tem1
            endif
          enddo
        enddo
      endif

!  ---  effective ice cloud droplet radius

      do k = 1, nlay
        do i = 1, ix
          tem1 = tlyr(i,k) - con_ttp
          tem2 = cip(i,k)

          if (tem2 > 0.0) then
            tem3 = tem2d(i,k) * tem2 / tvly(i,k)

            if (tem1 < -50.0) then
              rei(i,k) = (1250.0/9.917) * tem3 ** 0.109
            elseif (tem1 < -40.0) then
              rei(i,k) = (1250.0/9.337) * tem3 ** 0.08
            elseif (tem1 < -30.0) then
              rei(i,k) = (1250.0/9.208) * tem3 ** 0.055
            else
              rei(i,k) = (1250.0/9.387) * tem3 ** 0.031
            endif

!           if (lprnt .and. k == l) print *,' reil=',rei(i,k),' icec=', &
!    &        icec,' cip=',cip(i,k),' tem=',tem,' delt=',delt

            rei(i,k)   = max(10.0, min(rei(i,k), 300.0))
!           rei(i,k)   = max(20.0, min(rei(i,k), 300.0))
!!!!        rei(i,k)   = max(30.0, min(rei(i,k), 300.0))
!           rei(i,k)   = max(50.0, min(rei(i,k), 300.0))
!           rei(i,k)   = max(100.0, min(rei(i,k), 300.0))
          endif
        enddo
      enddo
!
      do k = 1, nlay
        do i = 1, ix
          clouds(i,k,1) = cldtot(i,k)
          clouds(i,k,2) = cwp(i,k)
          clouds(i,k,3) = rew(i,k)
          clouds(i,k,4) = cip(i,k)
          clouds(i,k,5) = rei(i,k)
          clouds(i,k,6) = crp(i,k)
          clouds(i,k,7) = rer(i,k)
!         clouds(i,k,8) = csp(i,k)               !ncar scheme
          clouds(i,k,8) = csp(i,k) * rsden(i,k)  !fu's scheme
          clouds(i,k,9) = rei(i,k)
        enddo
      enddo


!  ---  compute low, mid, high, total, and boundary layer cloud fractions
!       and clouds top/bottom layer indices for low, mid, and high clouds.
!       the three cloud domain boundaries are defined by ptopc.  the cloud
!       overlapping method is defined by control flag 'iovr', which may
!       be different for lw and sw radiation programs.

      call gethml                                                       &
!  ---  inputs:
     &     ( plyr, ptop1, cldtot, cldcnv,                               &
     &       ix,nlay,                                                   &
!  ---  outputs:
     &       clds, mtop, mbot                                           &
     &     )


!
      return
!...................................
      end subroutine progcld2
!-----------------------------------


!-----------------------------------
      subroutine progcld4                                               &
     &     ( plyr,plvl,tlyr,qlyr,qstl,rhly,clw,                         &    !  ---  inputs:
     &       xlat,xlon,slmsk,                                           & 
     &       ntrac,ntcw,ntiw,ntrw,ntsw,ntgl,                            &            
     &       ix, nlay, nlp1,                                            &
     &       uni_cld, lmfshal, lmfdeep2, cldcov,                        &    
     &       re_cloud,re_ice,re_snow,                                   & 
     &       clouds,clds,mtop,mbot                                      &    !  ---  outputs:
     &      )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:    progcld4    computes cloud related quantities using    !
!   Thompson/WSM6 cloud microphysics scheme.                !
!                                                                       !
! abstract:  this program computes cloud fractions from cloud           !
!   condensates,                                                        !
!   and computes the low, mid, high, total and boundary layer cloud     !
!   fractions and the vertical indices of low, mid, and high cloud      !
!   top and base.  the three vertical cloud domains are set up in the   !
!   initial subroutine "cld_init".                                      !
!                                                                       !
! usage:         call progcld4                                          !
!                                                                       !
! subprograms called:   gethml                                          !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  definition of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (IX,NLAY) : model layer mean pressure in mb (100Pa)           !
!   plvl  (IX,NLP1) : model level pressure in mb (100Pa)                !
!   tlyr  (IX,NLAY) : model layer mean temperature in k                 !
!   tvly  (IX,NLAY) : model layer virtual temperature in k              !
!   qlyr  (IX,NLAY) : layer specific humidity in gm/gm                  !
!   qstl  (IX,NLAY) : layer saturate humidity in gm/gm                  !
!   rhly  (IX,NLAY) : layer relative humidity (=qlyr/qstl)              !
!   clw   (IX,NLAY,ntrac) : layer cloud condensate amount               !
!   xlat  (IX)      : grid latitude in radians, default to pi/2 -> -pi/2!
!                     range, otherwise see in-line comment              !
!   xlon  (IX)      : grid longitude in radians  (not used)             !
!   slmsk (IX)      : sea/land mask array (sea:0,land:1,sea-ice:2)      !
!   IX              : horizontal dimention                              !
!   NLAY,NLP1       : vertical layer/level dimensions                   !
!   uni_cld         : logical - true for cloud fraction from shoc       !
!   lmfshal         : logical - true for mass flux shallow convection   !
!   lmfdeep2        : logical - true for mass flux deep convection      !
!   cldcov          : layer cloud fraction (used when uni_cld=.true.    !
!                                                                       !
! output variables:                                                     !
!   clouds(IX,NLAY,NF_CLDS) : cloud profiles                            !
!      clouds(:,:,1) - layer total cloud fraction                       !
!      clouds(:,:,2) - layer cloud liq water path         (g/m**2)      !
!      clouds(:,:,3) - mean eff radius for liq cloud      (micron)      !
!      clouds(:,:,4) - layer cloud ice water path         (g/m**2)      !
!      clouds(:,:,5) - mean eff radius for ice cloud      (micron)      !
!      clouds(:,:,6) - layer rain drop water path         not assigned  !
!      clouds(:,:,7) - mean eff radius for rain drop      (micron)      !
!  *** clouds(:,:,8) - layer snow flake water path        not assigned  !
!      clouds(:,:,9) - mean eff radius for snow flake     (micron)      !
!  *** fu's scheme need to be normalized by snow density (g/m**3/1.0e6) !
!   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !
!                                                                       !
! module variables:                                                     !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   lmfshal         : mass-flux shallow conv scheme flag                !
!   lmfdeep2        : scale-aware mass-flux deep conv scheme flag       !
!   lcrick          : control flag for eliminating CRICK                !
!                     =t: apply layer smoothing to eliminate CRICK      !
!                     =f: do not apply layer smoothing                  !
!   lcnorm          : control flag for in-cld condensate                !
!                     =t: normalize cloud condensate                    !
!                     =f: not normalize cloud condensate                !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs
      integer,  intent(in) :: IX, NLAY, NLP1
      integer,  intent(in) :: ntrac, ntcw, ntiw, ntrw, ntsw, ntgl

      logical, intent(in)  :: uni_cld, lmfshal, lmfdeep2

      real (kind=kind_phys), dimension(:,:), intent(in) :: plvl, plyr,  &
     &       tlyr, qlyr, qstl, rhly, cldcov,                            &
     &       re_cloud, re_ice, re_snow 

      real (kind=kind_phys), dimension(:,:,:), intent(in) :: clw

      real (kind=kind_phys), dimension(:),   intent(in) :: xlat, xlon,  &
     &       slmsk

!  ---  outputs
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: clouds

      real (kind=kind_phys), dimension(:,:),   intent(out) :: clds

      integer,               dimension(:,:),   intent(out) :: mtop,mbot

!  ---  local variables:
      real (kind=kind_phys), dimension(IX,NLAY) :: cldtot, cldcnv,      &
     &       cwp, cip, crp, csp, rew, rei, res, rer, delp, tem2d, clwf

      real (kind=kind_phys) :: ptop1(IX,NK_CLDS+1)

      real (kind=kind_phys) :: clwmin, clwm, clwt, onemrh, value,       &
     &       tem1, tem2, tem3

      integer :: i, k, id, nf

!  ---  constant values
!     real (kind=kind_phys), parameter :: xrc3 = 200.
      real (kind=kind_phys), parameter :: xrc3 = 100.

!
!===> ... begin here
!
      do nf=1,nf_clds
        do k=1,nlay
          do i=1,ix
            clouds(i,k,nf) = 0.0
          enddo
        enddo
      enddo
!     clouds(:,:,:) = 0.0

      do k = 1, NLAY
        do i = 1, IX
          cldtot(i,k) = 0.0
          cldcnv(i,k) = 0.0
          cwp   (i,k) = 0.0
          cip   (i,k) = 0.0
          crp   (i,k) = 0.0
          csp   (i,k) = 0.0
          rew   (i,k) = re_cloud(i,k)
          rei   (i,k) = re_ice(i,k) 
          rer   (i,k) = rrain_def            ! default rain radius to 1000 micron
          res   (i,k) = re_snow(i,K) 
!         tem2d (i,k) = min( 1.0, max( 0.0, (con_ttp-tlyr(i,k))*0.05 ) )
          clwf(i,k)   = 0.0
        enddo
      enddo
!
!      if ( lcrick ) then
!        do i = 1, IX
!          clwf(i,1)    = 0.75*clw(i,1)    + 0.25*clw(i,2)
!          clwf(i,nlay) = 0.75*clw(i,nlay) + 0.25*clw(i,nlay-1)
!        enddo
!        do k = 2, NLAY-1
!          do i = 1, IX
!            clwf(i,K) = 0.25*clw(i,k-1) + 0.5*clw(i,k) + 0.25*clw(i,k+1)
!          enddo
!        enddo
!      else
!        do k = 1, NLAY
!          do i = 1, IX
!            clwf(i,k) = clw(i,k)
!          enddo
!        enddo
!      endif

        do k = 1, NLAY
          do i = 1, IX
            clwf(i,k) = clw(i,k,ntcw) +  clw(i,k,ntiw) + clw(i,k,ntsw)
          enddo
        enddo
!> -# Find top pressure for each cloud domain for given latitude.
!     ptopc(k,i): top presure of each cld domain (k=1-4 are sfc,L,m,h;
!  ---  i=1,2 are low-lat (<45 degree) and pole regions)

      do id = 1, 4
        tem1 = ptopc(id,2) - ptopc(id,1)

        do i =1, IX
          tem2 = xlat(i) / con_pi        ! if xlat in pi/2 -> -pi/2 range
!         tem2 = 0.5 - xlat(i)/con_pi    ! if xlat in 0 -> pi range

          ptop1(i,id) = ptopc(id,1) + tem1*max( 0.0, 4.0*abs(tem2)-1.0 )
        enddo
      enddo

!> -# Compute cloud liquid/ice condensate path in \f$ g/m^2 \f$ .

      if ( ivflip == 0 ) then          ! input data from toa to sfc
        do k = 1, NLAY
          do i = 1, IX
            delp(i,k) = plvl(i,k+1) - plvl(i,k)
            cwp(i,k) = max(0.0, clw(i,k,ntcw) * gfac * delp(i,k))
            cip(i,k) = max(0.0, clw(i,k,ntiw) * gfac * delp(i,k))
            crp(i,k) = max(0.0, clw(i,k,ntrw) * gfac * delp(i,k))
            csp(i,k) = max(0.0, (clw(i,k,ntsw)+clw(i,k,ntgl)) *         &
     &                  gfac * delp(i,k))
          enddo
        enddo
      else                             ! input data from sfc to toa
        do k = 1, NLAY
          do i = 1, IX
            delp(i,k) = plvl(i,k) - plvl(i,k+1)
            cwp(i,k) = max(0.0, clw(i,k,ntcw) * gfac * delp(i,k))
            cip(i,k) = max(0.0, clw(i,k,ntiw) * gfac * delp(i,k))
            crp(i,k) = max(0.0, clw(i,k,ntrw) * gfac * delp(i,k))
            csp(i,k) = max(0.0, (clw(i,k,ntsw)+clw(i,k,ntgl)) *         &
     &                  gfac * delp(i,k))
          enddo
        enddo
      endif                            ! end_if_ivflip

      if (uni_cld) then     ! use unified sgs clouds generated outside
        do k = 1, NLAY
          do i = 1, IX
            cldtot(i,k) = cldcov(i,k)
          enddo
        enddo

      else

!> -# Calculate layer cloud fraction.

      if ( ivflip == 0 ) then              ! input data from toa to sfc

        clwmin = 0.0
        if (.not. lmfshal) then
          do k = NLAY, 1, -1
          do i = 1, IX
            clwt = 1.0e-6 * (plyr(i,k)*0.001)
!           clwt = 2.0e-6 * (plyr(i,k)*0.001)

            if (clwf(i,k) > clwt) then

              onemrh= max( 1.e-10, 1.0-rhly(i,k) )
              clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )

              tem1  = min(max(sqrt(sqrt(onemrh*qstl(i,k))),0.0001),1.0)
              tem1  = 2000.0 / tem1
!             tem1  = 1000.0 / tem1

              value = max( min( tem1*(clwf(i,k)-clwm), 50.0 ), 0.0 )
              tem2  = sqrt( sqrt(rhly(i,k)) )

              cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
            endif
          enddo
          enddo
        else
          do k = NLAY, 1, -1
          do i = 1, IX
            clwt = 1.0e-6 * (plyr(i,k)*0.001)
!           clwt = 2.0e-6 * (plyr(i,k)*0.001)

            if (clwf(i,k) > clwt) then
              onemrh= max( 1.e-10, 1.0-rhly(i,k) )
              clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )
!
              tem1  = min(max((onemrh*qstl(i,k))**0.49,0.0001),1.0)  !jhan
              if (lmfdeep2) then
                tem1  = xrc3 / tem1
              else
                tem1  = 100.0 / tem1
              endif
!
              value = max( min( tem1*(clwf(i,k)-clwm), 50.0 ), 0.0 )
              tem2  = sqrt( sqrt(rhly(i,k)) )
              cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
            endif
          enddo
          enddo
        endif

      else                                 ! input data from sfc to toa

        clwmin = 0.0
        if (.not. lmfshal) then
          do k = 1, NLAY
          do i = 1, IX
            clwt = 1.0e-6 * (plyr(i,k)*0.001)
!           clwt = 2.0e-6 * (plyr(i,k)*0.001)

            if (clwf(i,k) > clwt) then

              onemrh= max( 1.e-10, 1.0-rhly(i,k) )
              clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )

              tem1  = min(max(sqrt(sqrt(onemrh*qstl(i,k))),0.0001),1.0)
              tem1  = 2000.0 / tem1

!             tem1  = 1000.0 / tem1

              value = max( min( tem1*(clwf(i,k)-clwm), 50.0 ), 0.0 )
              tem2  = sqrt( sqrt(rhly(i,k)) )

              cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
            endif
          enddo
          enddo
        else
          do k = 1, NLAY
          do i = 1, IX
            clwt = 1.0e-6 * (plyr(i,k)*0.001)
!           clwt = 2.0e-6 * (plyr(i,k)*0.001)

            if (clwf(i,k) > clwt) then
              onemrh= max( 1.e-10, 1.0-rhly(i,k) )
              clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )
!
              tem1  = min(max((onemrh*qstl(i,k))**0.49,0.0001),1.0)  !jhan
              if (lmfdeep2) then
                tem1  = xrc3 / tem1
              else
                tem1  = 100.0 / tem1
              endif
!
              value = max( min( tem1*(clwf(i,k)-clwm), 50.0 ), 0.0 )
              tem2  = sqrt( sqrt(rhly(i,k)) )

              cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
            endif
          enddo
          enddo
        endif

      endif                                ! end_if_flip
      endif                                ! if (uni_cld) then

      do k = 1, NLAY
        do i = 1, IX
          if (cldtot(i,k) < climit) then
            cldtot(i,k) = 0.0
            cwp(i,k)    = 0.0
            cip(i,k)    = 0.0
            crp(i,k)    = 0.0
            csp(i,k)    = 0.0
          endif
        enddo
      enddo

      if ( lcnorm ) then
        do k = 1, NLAY
          do i = 1, IX
            if (cldtot(i,k) >= climit) then
              tem1 = 1.0 / max(climit2, cldtot(i,k))
              cwp(i,k) = cwp(i,k) * tem1
              cip(i,k) = cip(i,k) * tem1
              crp(i,k) = crp(i,k) * tem1
              csp(i,k) = csp(i,k) * tem1
            endif
          enddo
        enddo
      endif

!
      do k = 1, NLAY
        do i = 1, IX
          clouds(i,k,1) = cldtot(i,k)
          clouds(i,k,2) = cwp(i,k)
          clouds(i,k,3) = rew(i,k)
          clouds(i,k,4) = cip(i,k)
          clouds(i,k,5) = rei(i,k)
          clouds(i,k,6) = crp(i,k)  ! added for Thompson 
          clouds(i,k,7) = rer(i,k)
          clouds(i,k,8) = csp(i,k)  ! added for Thompson 
          clouds(i,k,9) = res(i,k)
        enddo
      enddo


!> -# Call gethml() to compute low,mid,high,total, and boundary layer
!!    cloud fractions and clouds top/bottom layer indices for low, mid,
!!    and high clouds.
!  ---  compute low, mid, high, total, and boundary layer cloud fractions
!       and clouds top/bottom layer indices for low, mid, and high clouds.
!       The three cloud domain boundaries are defined by ptopc.  The cloud
!       overlapping method is defined by control flag 'iovr', which may
!       be different for lw and sw radiation programs.

      call gethml                                                       &
!  ---  inputs:
     &     ( plyr, ptop1, cldtot, cldcnv,                               &
     &       ix,nlay,                                                   &
!  ---  outputs:
     &       clds, mtop, mbot                                           &
     &     )


!
      return
!...................................
      end subroutine progcld4
!-----------------------------------


!-----------------------------------
      subroutine progcld5                                               &
!...................................
!
!  ---  inputs:
           ( plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,cnvw,cnvc,          &
             xlat,xlon,slmsk,IX, NLAY, NLP1,                            &
             cldcov,                                                    &
!  ---  outputs:
     &       clouds,clds,mtop,mbot                                      &
     &      )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:    progcld5    computes cloud related quantities using    !
!   GFDL Lin MP prognostic cloud microphysics scheme.                   !
!                                                                       !
! abstract:  this program computes cloud fractions from cloud           !
!   condensates, calculates liquid/ice cloud droplet effective radius,  !
!   and computes the low, mid, high, total and boundary layer cloud     !
!   fractions and the vertical indices of low, mid, and high cloud      !
!   top and base.  the three vertical cloud domains are set up in the   !
!   initial subroutine "cld_init".                                      !
!                                                                       !
! usage:         call progcld5                                          !
!                                                                       !
! subprograms called:   gethml                                          !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  definition of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (IX,NLAY) : model layer mean pressure in mb (100Pa)           !
!   plvl  (IX,NLP1) : model level pressure in mb (100Pa)                !
!   tlyr  (IX,NLAY) : model layer mean temperature in k                 !
!   tvly  (IX,NLAY) : model layer virtual temperature in k              !
!   qlyr  (IX,NLAY) : layer specific humidity in gm/gm                  !
!   qstl  (IX,NLAY) : layer saturate humidity in gm/gm                  !
!   rhly  (IX,NLAY) : layer relative humidity (=qlyr/qstl)              !
!   clw   (IX,NLAY) : layer cloud condensate amount                     !
!   cnvw  (IX,NLAY) : layer convective cloud condensate                 !
!   cnvc  (IX,NLAY) : layer convective cloud cover                      !
!   xlat  (IX)      : grid latitude in radians, default to pi/2 -> -pi/2!
!                     range, otherwise see in-line comment              !
!   xlon  (IX)      : grid longitude in radians  (not used)             !
!   slmsk (IX)      : sea/land mask array (sea:0,land:1,sea-ice:2)      !
!   IX              : horizontal dimention                              !
!   NLAY,NLP1       : vertical layer/level dimensions                   !
!   cldcov(IX,NLAY) : total cloud fraction                              !
!                                                                       !
! output variables:                                                     !
!   clouds(IX,NLAY,NF_CLDS) : cloud profiles                            !
!      clouds(:,:,1) - layer total cloud fraction, =cldcov(IX,NLAY)     !
!      clouds(:,:,2) - layer cloud liq water path         (g/m**2)      !
!      clouds(:,:,3) - mean eff radius for liq cloud      (micron)      !
!      clouds(:,:,4) - layer cloud ice water path         (g/m**2)      !
!      clouds(:,:,5) - mean eff radius for ice cloud      (micron)      !
!      clouds(:,:,6) - layer rain drop water path         not assigned  !
!      clouds(:,:,7) - mean eff radius for rain drop      (micron)      !
!  *** clouds(:,:,8) - layer snow flake water path        not assigned  !
!      clouds(:,:,9) - mean eff radius for snow flake     (micron)      !
!  *** fu's scheme need to be normalized by snow density (g/m**3/1.0e6) !
!   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !
!                                                                       !
! module variables:                                                     !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   lsashal         : control flag for shallow convection               !
!   lcrick          : control flag for eliminating CRICK                !
!                     =t: apply layer smoothing to eliminate CRICK      !
!                     =f: do not apply layer smoothing                  !
!   lcnorm          : control flag for in-cld condensate                !
!                     =t: normalize cloud condensate                    !
!                     =f: not normalize cloud condensate                !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs
      integer,  intent(in) :: IX, NLAY, NLP1

      real (kind=kind_phys), dimension(:,:), intent(in) :: plvl, plyr,  &
     &       tlyr, tvly, qlyr, qstl, rhly, clw, cnvw, cnvc, cldcov!,     &
!             re_cloud, re_ice, re_snow

      real (kind=kind_phys), dimension(:),   intent(in) :: xlat, xlon,  &
     &       slmsk

!  ---  outputs
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: clouds

      real (kind=kind_phys), dimension(:,:),   intent(out) :: clds

      integer,               dimension(:,:),   intent(out) :: mtop,mbot

!  ---  local variables:
      real (kind=kind_phys), dimension(IX,NLAY) :: cldcnv,              &
     &       cwp, cip, crp, csp, rew, rei, res, rer, delp, tem2d, clwf, &
     &       cldtot

      real (kind=kind_phys) :: ptop1(IX,NK_CLDS+1)

      real (kind=kind_phys) :: clwmin, clwm, clwt, onemrh, value,       &
     &       tem1, tem2, tem3

      integer :: i, k, id, nf

!
!===> ... begin here
!
      cldtot = cldcov
      do nf=1,nf_clds
        do k=1,nlay
          do i=1,ix
            clouds(i,k,nf) = 0.0
          enddo
        enddo
      enddo
!     clouds(:,:,:) = 0.0

      do k = 1, NLAY
        do i = 1, IX
          cldcnv(i,k) = 0.0
          cwp   (i,k) = 0.0
          cip   (i,k) = 0.0
          crp   (i,k) = 0.0
          csp   (i,k) = 0.0
          rew   (i,k) = reliq_def            ! default liq radius to 10 micron
          rei   (i,k) = reice_def            ! default ice radius to 50 micron
          rer   (i,k) = rrain_def            ! default rain radius to 1000 micron
          res   (i,k) = rsnow_def            ! default snow radius to 250 micron
          tem2d (i,k) = min( 1.0, max( 0.0, (con_ttp-tlyr(i,k))*0.05 ) )
          clwf(i,k)   = 0.0
        enddo
      enddo
!
      if ( lcrick ) then
        do i = 1, IX
          clwf(i,1)    = 0.75*clw(i,1)    + 0.25*clw(i,2)
          clwf(i,nlay) = 0.75*clw(i,nlay) + 0.25*clw(i,nlay-1)
        enddo
        do k = 2, NLAY-1
          do i = 1, IX
            clwf(i,K) = 0.25*clw(i,k-1) + 0.5*clw(i,k) + 0.25*clw(i,k+1)
          enddo
        enddo
      else
        do k = 1, NLAY
          do i = 1, IX
            clwf(i,k) = clw(i,k)
          enddo
        enddo
      endif

!  ---  find top pressure for each cloud domain for given latitude
!       ptopc(k,i): top presure of each cld domain (k=1-4 are sfc,L,m,h;
!  ---  i=1,2 are low-lat (<45 degree) and pole regions)

      do id = 1, 4
        tem1 = ptopc(id,2) - ptopc(id,1)

        do i =1, IX
          tem2 = xlat(i) / con_pi        ! if xlat in pi/2 -> -pi/2 range
!         tem2 = 0.5 - xlat(i)/con_pi    ! if xlat in 0 -> pi range

          ptop1(i,id) = ptopc(id,1) + tem1*max( 0.0, 4.0*abs(tem2)-1.0 )
        enddo
      enddo

!  ---  compute liquid/ice condensate path in g/m**2

      if ( ivflip == 0 ) then          ! input data from toa to sfc
        do k = 1, NLAY
          do i = 1, IX
            delp(i,k) = plvl(i,k+1) - plvl(i,k)
            clwt     = max(0.0,(clwf(i,k)+cnvw(i,k))) * gfac * delp(i,k)
            cip(i,k) = clwt * tem2d(i,k)
            cwp(i,k) = clwt - cip(i,k)
          enddo
        enddo
      else                             ! input data from sfc to toa
        do k = 1, NLAY
          do i = 1, IX
            delp(i,k) = plvl(i,k) - plvl(i,k+1)
            clwt     = max(0.0,(clwf(i,k)+cnvw(i,k))) * gfac * delp(i,k)
            cip(i,k) = clwt * tem2d(i,k)
            cwp(i,k) = clwt - cip(i,k)
          enddo
        enddo
      endif                            ! end_if_ivflip

!  ---  effective liquid cloud droplet radius over land

      do i = 1, IX
        if (nint(slmsk(i)) == 1) then
          do k = 1, NLAY
            rew(i,k) = 5.0 + 5.0 * tem2d(i,k)
          enddo
        endif
      enddo

      do k = 1, NLAY
        do i = 1, IX
          if (cldtot(i,k) < climit) then
            cwp(i,k)    = 0.0
            cip(i,k)    = 0.0
            crp(i,k)    = 0.0
            csp(i,k)    = 0.0
          endif
        enddo
      enddo

      if ( lcnorm ) then
        do k = 1, NLAY
          do i = 1, IX
            if (cldtot(i,k) >= climit) then
              tem1 = 1.0 / max(climit2, cldtot(i,k))
              cwp(i,k) = cwp(i,k) * tem1
              cip(i,k) = cip(i,k) * tem1
              crp(i,k) = crp(i,k) * tem1
              csp(i,k) = csp(i,k) * tem1
            endif
          enddo
        enddo
      endif

!  ---  effective ice cloud droplet radius

      do k = 1, NLAY
        do i = 1, IX
          tem2 = tlyr(i,k) - con_ttp

          if (cip(i,k) > 0.0) then
            tem3 = gord * cip(i,k) * plyr(i,k) / (delp(i,k)*tvly(i,k))

            if (tem2 < -50.0) then
              rei(i,k) = (1250.0/9.917) * tem3 ** 0.109
            elseif (tem2 < -40.0) then
              rei(i,k) = (1250.0/9.337) * tem3 ** 0.08
            elseif (tem2 < -30.0) then
              rei(i,k) = (1250.0/9.208) * tem3 ** 0.055
            else
              rei(i,k) = (1250.0/9.387) * tem3 ** 0.031
            endif
!           rei(i,k)   = max(20.0, min(rei(i,k), 300.0))
!           rei(i,k)   = max(10.0, min(rei(i,k), 100.0))
            rei(i,k)   = max(10.0, min(rei(i,k), 150.0))
!           rei(i,k)   = max(5.0,  min(rei(i,k), 130.0))
          endif
        enddo
      enddo

!
      do k = 1, NLAY
        do i = 1, IX
          clouds(i,k,1) = cldtot(i,k)
          clouds(i,k,2) = cwp(i,k)
          clouds(i,k,3) = rew(i,k)
          clouds(i,k,4) = cip(i,k)
          clouds(i,k,5) = rei(i,k)
!         clouds(i,k,6) = 0.0
          clouds(i,k,7) = rer(i,k)
!         clouds(i,k,8) = 0.0
          clouds(i,k,9) = res(i,k)
        enddo
      enddo


!  ---  compute low, mid, high, total, and boundary layer cloud fractions
!       and clouds top/bottom layer indices for low, mid, and high clouds.
!       The three cloud domain boundaries are defined by ptopc.  The cloud
!       overlapping method is defined by control flag 'iovr', which may
!       be different for lw and sw radiation programs.

      call gethml                                                       &
!  ---  inputs:
     &     ( plyr, ptop1, cldtot, cldcnv,                               &
     &       IX,NLAY,                                                   &
!  ---  outputs:
     &       clds, mtop, mbot                                           &
     &     )


!
      return
!...................................
      end subroutine progcld5
!-----------------------------------



!-----------------------------------
      subroutine progcld5o                                              &
!...................................
!
!  ---  inputs:
     &     ( plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,                    &
     &       xlat,xlon,slmsk,                                           &
     &       ntrac,ntcw,ntiw,ntrw,ntsw,ntgl,cldcov,                     &
     &       effr_cw,effr_iw,effr_sw,effr_rw,effr_in,                   &
     &       IX, NLAY, NLP1,                                            &
!  ---  outputs:
     &       clouds,clds,mtop,mbot                                      &
     &      )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:    progcld5o   computes cloud related quantities using    !
!   GFDL Lin MP prognostic cloud microphysics scheme. Moist species     !
!   from MP are fed into the corresponding arrays for calcuation of     !
!                                                                       !
! abstract:  this program computes cloud fractions from cloud           !
!   condensates, calculates liquid/ice cloud droplet effective radius,  !
!   and computes the low, mid, high, total and boundary layer cloud     !
!   fractions and the vertical indices of low, mid, and high cloud      !
!   top and base.  the three vertical cloud domains are set up in the   !
!   initial subroutine "cld_init".                                      !
!                                                                       !
! usage:         call progcld5o                                         !
!                                                                       !
! subprograms called:   gethml                                          !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  definition of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (IX,NLAY) : model layer mean pressure in mb (100Pa)           !
!   plvl  (IX,NLP1) : model level pressure in mb (100Pa)                !
!   tlyr  (IX,NLAY) : model layer mean temperature in k                 !
!   tvly  (IX,NLAY) : model layer virtual temperature in k              !
!   qlyr  (IX,NLAY) : layer specific humidity in gm/gm                  !
!   qstl  (IX,NLAY) : layer saturate humidity in gm/gm                  !
!   rhly  (IX,NLAY) : layer relative humidity (=qlyr/qstl)              !
!   clw   (IX,NLAY,NTRAC) : layer cloud condensate amount               !
!   xlat  (IX)      : grid latitude in radians, default to pi/2 -> -pi/2!
!                     range, otherwise see in-line comment              !
!   xlon  (IX)      : grid longitude in radians  (not used)             !
!   slmsk (IX)      : sea/land mask array (sea:0,land:1,sea-ice:2)      !
!   ntrac           : number of tracers
!   ntcw            : index for liquid water
!   ntiw            : index for ice water
!   ntrw            : index for rain water
!   ntsw            : index for snow water
!   ntgl            : index for graupel
!   cldcov(IX,NLAY) : total cloud fraction
!   effr_cw(IX,NLAY): effective radius of liquid water (mm)             !
!   effr_iw(IX,NLAY): effective radius of ice water (mm)                !
!   effr_rw(IX,NLAY): effective radius of rain water (mm)               !
!   effr_sw(IX,NLAY): effective radius of snow water (mm)               !
!   IX              : horizontal dimention                              !
!   NLAY,NLP1       : vertical layer/level dimensions                   !
!                                                                       !
! output variables:                                                     !
!   clouds(IX,NLAY,NF_CLDS) : cloud profiles                            !
!      clouds(:,:,1) - layer total cloud fraction                       !
!      clouds(:,:,2) - layer cloud liq water path         (g/m**2)      !
!      clouds(:,:,3) - mean eff radius for liq cloud      (micron)      !
!      clouds(:,:,4) - layer cloud ice water path         (g/m**2)      !
!      clouds(:,:,5) - mean eff radius for ice cloud      (micron)      !
!      clouds(:,:,6) - layer rain drop water path         not assigned  !
!      clouds(:,:,7) - mean eff radius for rain drop      (micron)      !
!  *** clouds(:,:,8) - layer snow flake water path        not assigned  !
!      clouds(:,:,9) - mean eff radius for snow flake     (micron)      !
!  *** fu's scheme need to be normalized by snow density (g/m**3/1.0e6) !
!   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !
!                                                                       !
! module variables:                                                     !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   lsashal         : control flag for shallow convection               !
!   lcrick          : control flag for eliminating CRICK                !
!                     =t: apply layer smoothing to eliminate CRICK      !
!                     =f: do not apply layer smoothing                  !
!   lcnorm          : control flag for in-cld condensate                !
!                     =t: normalize cloud condensate                    !
!                     =f: not normalize cloud condensate                !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs
      integer,  intent(in) :: IX, NLAY, NLP1
      integer,  intent(in) :: ntrac, ntcw, ntiw, ntrw, ntsw, ntgl

      real (kind=kind_phys), dimension(:,:), intent(in) :: plvl, plyr,  &
     &       tlyr, tvly, qlyr, qstl, rhly, cldcov


      real (kind=kind_phys), dimension(:,:,:), intent(in) :: clw
      real (kind=kind_phys), dimension(:),   intent(in) :: xlat, xlon,  &
     &       slmsk

      real (kind=kind_phys), dimension(:,:), intent(in) :: effr_cw,     &
     &       effr_iw, effr_rw, effr_sw

      logical, intent(in) :: effr_in

!  ---  outputs
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: clouds

      real (kind=kind_phys), dimension(:,:),   intent(out) :: clds

      integer,               dimension(:,:),   intent(out) :: mtop,mbot

!  ---  local variables:
      real (kind=kind_phys), dimension(IX,NLAY) :: cldcnv,              &
     &       cwp, cip, crp, csp, rew, rei, res, rer, delp, tem2d

      real (kind=kind_phys) :: ptop1(IX,NK_CLDS+1)

      real (kind=kind_phys) :: clwmin, clwm, clwt, onemrh, value,       &
     &       tem1, tem2, tem3
      real (kind=kind_phys), dimension(IX,NLAY) :: cldtot

      integer :: i, k, id, nf

!
!===> ... begin here
!
      do nf=1,nf_clds
        do k=1,nlay
          do i=1,ix
            clouds(i,k,nf) = 0.0
          enddo
        enddo
      enddo
!     clouds(:,:,:) = 0.0

      do k = 1, NLAY
        do i = 1, IX
          cldcnv(i,k) = 0.0
          cwp   (i,k) = 0.0
          cip   (i,k) = 0.0
          crp   (i,k) = 0.0
          csp   (i,k) = 0.0
          if ( effr_in ) then
            rew (i,k) = effr_cw (i,k)
            rei (i,k) = effr_iw (i,k)
            rer (i,k) = effr_rw (i,k)
            res (i,k) = effr_sw (i,k)
          else
            rew (i,k) = reliq_def            ! default liq radius to 10 micron
            rei (i,k) = reice_def            ! default ice radius to 50 micron
            rer (i,k) = rrain_def            ! default rain radius to 1000 micron
            res (i,k) = rsnow_def            ! default snow radius to 250 micron
          endif
          tem2d (i,k) = min( 1.0, max( 0.0, (con_ttp-tlyr(i,k))*0.05 ) )
          cldtot(i,k) = cldcov(i,k)
        enddo
      enddo

!  ---  find top pressure for each cloud domain for given latitude
!       ptopc(k,i): top presure of each cld domain (k=1-4 are sfc,L,m,h;
!  ---  i=1,2 are low-lat (<45 degree) and pole regions)

      do id = 1, 4
        tem1 = ptopc(id,2) - ptopc(id,1)

        do i =1, IX
          tem2 = xlat(i) / con_pi        ! if xlat in pi/2 -> -pi/2 range
!         tem2 = 0.5 - xlat(i)/con_pi    ! if xlat in 0 -> pi range

          ptop1(i,id) = ptopc(id,1) + tem1*max( 0.0, 4.0*abs(tem2)-1.0 )
        enddo
      enddo

!  ---  compute liquid/ice condensate path in g/m**2

      if ( ivflip == 0 ) then          ! input data from toa to sfc
        do k = 1, NLAY
          do i = 1, IX
            delp(i,k) = plvl(i,k+1) - plvl(i,k)
            cwp(i,k) = max(0.0, clw(i,k,ntcw) * gfac * delp(i,k))
            cip(i,k) = max(0.0, clw(i,k,ntiw) * gfac * delp(i,k))
            crp(i,k) = max(0.0, clw(i,k,ntrw) * gfac * delp(i,k))
            csp(i,k) = max(0.0, (clw(i,k,ntsw)+clw(i,k,ntgl)) *         &
     &                  gfac * delp(i,k))
          enddo
        enddo
      else                             ! input data from sfc to toa
        do k = 1, NLAY
          do i = 1, IX
            delp(i,k) = plvl(i,k) - plvl(i,k+1)
            cwp(i,k) = max(0.0, clw(i,k,ntcw) * gfac * delp(i,k))
            cip(i,k) = max(0.0, clw(i,k,ntiw) * gfac * delp(i,k))
            crp(i,k) = max(0.0, clw(i,k,ntrw) * gfac * delp(i,k))
            csp(i,k) = max(0.0, (clw(i,k,ntsw)+clw(i,k,ntgl)) *         &
     &                  gfac * delp(i,k))
          enddo
        enddo
      endif                            ! end_if_ivflip

!  ---  effective liquid cloud droplet radius over land

      if ( .not. effr_in ) then
      do i = 1, IX
        if (nint(slmsk(i)) == 1) then
          do k = 1, NLAY
            rew(i,k) = 5.0 + 5.0 * tem2d(i,k)
          enddo
        endif
      enddo
      endif

      do k = 1, NLAY
        do i = 1, IX
          if (cldtot(i,k) < climit) then
            cldtot(i,k) = 0.0 
            cwp(i,k)    = 0.0
            cip(i,k)    = 0.0
            crp(i,k)    = 0.0
            csp(i,k)    = 0.0
          endif
        enddo
      enddo

      if ( lcnorm ) then
        do k = 1, NLAY
          do i = 1, IX
            if (cldtot(i,k) >= climit) then
              tem1 = 1.0 / max(climit2, cldtot(i,k))
              cwp(i,k) = cwp(i,k) * tem1
              cip(i,k) = cip(i,k) * tem1
              crp(i,k) = crp(i,k) * tem1
              csp(i,k) = csp(i,k) * tem1
            endif
          enddo
        enddo
      endif

!  ---  effective ice cloud droplet radius

      if ( .not. effr_in ) then
      do k = 1, NLAY
        do i = 1, IX
          tem2 = tlyr(i,k) - con_ttp

          if (cip(i,k) > 0.0) then
            tem3 = gord * cip(i,k) * plyr(i,k) / (delp(i,k)*tvly(i,k))

            if (tem2 < -50.0) then
              rei(i,k) = (1250.0/9.917) * tem3 ** 0.109
            elseif (tem2 < -40.0) then
              rei(i,k) = (1250.0/9.337) * tem3 ** 0.08
            elseif (tem2 < -30.0) then
              rei(i,k) = (1250.0/9.208) * tem3 ** 0.055
            else
              rei(i,k) = (1250.0/9.387) * tem3 ** 0.031
            endif
!           rei(i,k)   = max(20.0, min(rei(i,k), 300.0))
!           rei(i,k)   = max(10.0, min(rei(i,k), 100.0))
            rei(i,k)   = max(10.0, min(rei(i,k), 150.0))
!           rei(i,k)   = max(5.0,  min(rei(i,k), 130.0))
          endif
        enddo
      enddo
      endif

!
      do k = 1, NLAY
        do i = 1, IX
          clouds(i,k,1) = cldtot(i,k)
          clouds(i,k,2) = cwp(i,k)
          clouds(i,k,3) = rew(i,k)
          clouds(i,k,4) = cip(i,k)
          clouds(i,k,5) = rei(i,k)
          clouds(i,k,6) = crp(i,k) 
          clouds(i,k,7) = rer(i,k)
          clouds(i,k,8) = csp(i,k)
          clouds(i,k,9) = res(i,k)
        enddo
      enddo


!  ---  compute low, mid, high, total, and boundary layer cloud fractions
!       and clouds top/bottom layer indices for low, mid, and high clouds.
!       The three cloud domain boundaries are defined by ptopc.  The cloud
!       overlapping method is defined by control flag 'iovr', which may
!       be different for lw and sw radiation programs.

      call gethml                                                       &
!  ---  inputs:
     &     ( plyr, ptop1, cldtot, cldcnv,                               &
     &       IX,NLAY,                                                   &
!  ---  outputs:
     &       clds, mtop, mbot                                           &
     &     )


!
      return
!...................................
      end subroutine progcld5o
!-----------------------------------

!-----------------------------------
      subroutine progcld6                                               &
!...................................

!  ---  inputs:
     &     ( plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,cnvw,cnvc,              &
     &       qw,qr,qi,qs,qg,                                            &
     &       cldcov,slmsk,                                              &
     &       effr_cw,effr_iw,effr_sw,effr_rw,effr_in,                   &
     &       xlat,xlon,IX,NLAY,NLP1,                                    &
!  ---  outputs:
     &       clouds,clds,mtop,mbot                                      &
     &      )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:    progcld1    computes cloud related quantities using    !
!   zhao/moorthi's prognostic cloud microphysics scheme.                !
!                                                                       !
! abstract:  this program computes cloud fractions from cloud           !
!   condensates, calculates liquid/ice cloud droplet effective radius,  !
!   and computes the low, mid, high, total and boundary layer cloud     !
!   fractions and the vertical indices of low, mid, and high cloud      !
!   top and base.  the three vertical cloud domains are set up in the   !
!   initial subroutine "cld_init".                                      !
!                                                                       !
! usage:         call progcld1                                          !
!                                                                       !
! subprograms called:   gethml                                          !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  definition of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (IX,NLAY) : model layer mean pressure in mb (100Pa)           !
!   plvl  (IX,NLP1) : model level pressure in mb (100Pa)                !
!   tlyr  (IX,NLAY) : model layer mean temperature in k                 !
!   tvly  (IX,NLAY) : model layer virtual temperature in k              !
!   qlyr  (IX,NLAY) : layer specific humidity in gm/gm                  !
!   qstl  (IX,NLAY) : layer saturate humidity in gm/gm                  !
!   rhly  (IX,NLAY) : layer relative humidity (=qlyr/qstl)              !
!   cnvw  (ix,nlay) : layer convective cloud condensate                 !
!   cnvc  (ix,nlay) : layer convective cloud cover                      !
!   qw    (ix,nlay) : cloud water mixing ratio (kg/kg)                  !
!   qr    (ix,nlay) : rain water mixing ratio (kg/kg)                   !
!   qi    (ix,nlay) : cloud ice mixing ratio (kg/kg)                    !
!   qs    (ix,nlay) : snow water mixing ratio (kg/kg)                   !
!   qg    (ix,nlay) : graupel mixing ratio (kg/kg)                      !
!   qa    (ix,nlay) : aerosol mixing ratio (kg/kg)                      !
!   cldtot(ix,nlay) : layer cloud fraction                              !
!   slmsk (IX)      : sea/land mask array (sea:0,land:1,sea-ice:2)      !
!   snowd (ix)      : snow depth                                        !
!   xlat  (IX)      : grid latitude in radians, default to pi/2 -> -pi/2!
!                     range, otherwise see in-line comment              !
!   xlon  (IX)      : grid longitude in radians  (not used)             !
!   IX              : horizontal dimention                              !
!   NLAY,NLP1       : vertical layer/level dimensions                   !
!                                                                       !
! output variables:                                                     !
!   clouds(IX,NLAY,NF_CLDS) : cloud profiles                            !
!      clouds(:,:,1) - layer total cloud fraction                       !
!      clouds(:,:,2) - layer cloud liq water path         (g/m**2)      !
!      clouds(:,:,3) - mean eff radius for liq cloud      (micron)      !
!      clouds(:,:,4) - layer cloud ice water path         (g/m**2)      !
!      clouds(:,:,5) - mean eff radius for ice cloud      (micron)      !
!      clouds(:,:,6) - layer rain drop water path         not assigned  !
!      clouds(:,:,7) - mean eff radius for rain drop      (micron)      !
!  *** clouds(:,:,8) - layer snow flake water path        not assigned  !
!      clouds(:,:,9) - mean eff radius for snow flake     (micron)      !
!  *** fu's scheme need to be normalized by snow density (g/m**3/1.0e6) !
!   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !
!                                                                       !
! module variables:                                                     !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   lsashal         : control flag for shallow convection               !
!   lcrick          : control flag for eliminating CRICK                !
!                     =t: apply layer smoothing to eliminate CRICK      !
!                     =f: do not apply layer smoothing                  !
!   lcnorm          : control flag for in-cld condensate                !
!                     =t: normalize cloud condensate                    !
!                     =f: not normalize cloud condensate                !
!                                                                       !
!  ====================    end of description    =====================  !

      use physpara,           only : icmphys
!      use cld_eff_rad_v2,     only : cld_eff_rad
!      use module_mp_gfdl_v3,  only : cld_eff_rad_v3
!
      implicit none

!  ---  inputs
      integer,  intent(in) :: IX, NLAY, NLP1

      real (kind=kind_phys), dimension(:,:), intent(in) :: plvl, plyr,  &
     &       tlyr, tvly, qlyr, qstl, rhly,  cnvw, cnvc

      real (kind=kind_phys), dimension(:,:), intent(inout) ::           &
     &       qw, qr, qi, qs, qg
      real (kind=kind_phys), dimension(:,:), intent(inout) :: cldcov

      real (kind=kind_phys), dimension(:),   intent(in) :: xlat, xlon,  &
     &       slmsk

      real (kind=kind_phys), dimension(:,:), intent(in) :: effr_cw,     &
     &       effr_iw, effr_rw, effr_sw

      logical, intent(in) :: effr_in

!  ---  outputs
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: clouds

      real (kind=kind_phys), dimension(:,:),   intent(out) :: clds

      integer,               dimension(:,:),   intent(out) :: mtop,mbot

!  ---  local variables:
      real (kind=kind_phys), dimension(IX,NLAY) :: cldcnv, cldtot,      &
     &       cwp, cip, crp, csp, cgp, rew, rei, res, rer, reg, delp,    &
     &       tem2d, clwf

      real (kind=kind_phys) :: ptop1(IX,NK_CLDS+1)

      real (kind=kind_phys) :: clwmin, clwm, clwt, onemrh, value,       &
     &       tem1, tem2, tem3

      real, parameter :: qmin = 1.0e-12  !minimum mass mixing ratio (kg / kg)

      integer :: i, k, id, nf

!
!===> ... begin here
!
      do nf=1,nf_clds
        do k=1,nlay
          do i=1,ix
            clouds(i,k,nf) = 0.0
          enddo
        enddo
      enddo

      do k = 1, NLAY
        do i = 1, IX
          cldcnv(i,k) = 0.0
          cwp   (i,k) = 0.0
          cip   (i,k) = 0.0
          crp   (i,k) = 0.0
          csp   (i,k) = 0.0
          if ( effr_in ) then
            rew (i,k) = effr_cw (i,k)
            rei (i,k) = effr_iw (i,k)
            rer (i,k) = effr_rw (i,k)
            res (i,k) = effr_sw (i,k)
          else
            rew (i,k) = reliq_def            ! default liq radius to 10 micron
            rei (i,k) = reice_def            ! default ice radius to 50 micron
            rer (i,k) = rrain_def            ! default rain radius to 1000 micron
            res (i,k) = rsnow_def            ! default snow radius to 250 micron
          endif
          tem2d (i,k) = min( 1.0, max( 0.0, (con_ttp-tlyr(i,k))*0.05 ) )
          cldtot(i,k) = cldcov(i,k)
        enddo
      enddo

!  ---  compute liquid/ice condensate path in g/m**2

      do k = 1, NLAY
        do i = 1, IX
          if ( ivflip == 0 ) then          ! input data from toa to sfc
            delp(i,k) = plvl(i,k+1) - plvl(i,k)
          else                             ! input data from sfc to toa
            delp(i,k) = plvl(i,k) - plvl(i,k+1)
          endif

          if ( qw(i,k) .gt. qmin ) cwp(i,k) = qw(i,k) * gfac * delp(i,k)
          if ( qi(i,k) .gt. qmin ) cip(i,k) = qi(i,k) * gfac * delp(i,k)
          if ( qr(i,k) .gt. qmin ) crp(i,k) = qr(i,k) * gfac * delp(i,k)
          if ( qs(i,k)+qg(i,k) .gt. qmin ) &
                         csp(i,k) = (qs(i,k)+qg(i,k)) * gfac * delp(i,k)
        enddo
      enddo

!      if ( icmphys .eq. 12 ) then
!      call cld_eff_rad (1, IX, 1, NLAY, slmsk, plyr*100,                &
!     &                  abs(plvl(:,1:NLAY)-plvl(:,2:NLAY+1))*100,       &
!     &                  tlyr, qw, qi, qr, qs, qg, qa, cwp, cip, crp,    &
!     &                  csp, cgp, rew, rei, rer, res, reg, cldtot,      &
!     &                  snowd, cnvw=cnvw, cnvc=cnvc)
!      cldcnv = 0.0
!      elseif ( icmphys .eq. 13 ) then
!      call cld_eff_rad_v3 (1, IX, 1, NLAY, slmsk, plyr*100,             &
!     &                  abs(plvl(:,1:NLAY)-plvl(:,2:NLAY+1))*100,       &
!     &                  tlyr, qw, qi, qr, qs, qg, qa, cwp, cip, crp,    &
!     &                  csp, cgp, rew, rei, rer, res, reg, cldtot,      &
!     &                  snowd, cnvw=cnvw, cnvc=cnvc)
!      cldcnv = 0.0
!      endif


!  ---  find top pressure for each cloud domain for given latitude
!       ptopc(k,i): top presure of each cld domain (k=1-4 are sfc,L,m,h;
!  ---  i=1,2 are low-lat (<45 degree) and pole regions)

      do id = 1, 4
        tem1 = ptopc(id,2) - ptopc(id,1)

        do i =1, IX
          tem2 = xlat(i) / con_pi        ! if xlat in pi/2 -> -pi/2 range
!         tem2 = 0.5 - xlat(i)/con_pi    ! if xlat in 0 -> pi range

          ptop1(i,id) = ptopc(id,1) + tem1*max( 0.0, 4.0*abs(tem2)-1.0 )
        enddo
      enddo

      do k = 1, NLAY
        do i = 1, IX
          clouds(i,k,1) = cldtot(i,k)
          clouds(i,k,2) = cwp(i,k)
          clouds(i,k,3) = rew(i,k)
          clouds(i,k,4) = cip(i,k)
          clouds(i,k,5) = rei(i,k)
          clouds(i,k,6) = crp(i,k)
          clouds(i,k,7) = rer(i,k)
          clouds(i,k,8) = csp(i,k)
          clouds(i,k,9) = res(i,k)
!          clouds(i,k,10) = cgp(i,k)
!          clouds(i,k,11) = reg(i,k)
        enddo
      enddo


!  ---  compute low, mid, high, total, and boundary layer cloud fractions
!       and clouds top/bottom layer indices for low, mid, and high clouds.
!       The three cloud domain boundaries are defined by ptopc.  The cloud
!       overlapping method is defined by control flag 'iovr', which may
!       be different for lw and sw radiation programs.

      call gethml                                                       &
!  ---  inputs:
     &     ( plyr, ptop1, cldtot, cldcnv,                               &
     &       IX,NLAY,                                                   &
!  ---  outputs:
     &       clds, mtop, mbot                                           &
     &     )


!
      return
!...................................
      end subroutine progcld6
!-----------------------------------

!............................................
!> This subroutine added by G. Thompson specifically to account for
!! explicit (microphysics-produced) cloud liquid water, cloud ice, and
!! snow with 100% cloud fraction.  Also, a parameterization for cloud
!! fraction less than 1.0 but greater than 0.0 follows Mocko and Cotton
!! (1996) from Sundqvist et al. (1989) with cloud fraction increasing
!! as RH increases above a critical value.  In locations with non-zero
!! (but less than 1.0) cloud fraction, there MUST be a value assigned
!! to cloud liquid water and ice or else there is zero impact in the
!! RRTMG radiation scheme.
      subroutine progcld_thompson                                       &
!  --- inputs
           ( plyr,plvl,tlyr,qlyr,qstl,rhly,clw,                         &
             xlat,xlon,slmsk,                                           &
             ntrac,ntcw,ntiw,ntrw,ntsw,ntgl,                            &
             IX, NLAY, NLP1,                                            &
             uni_cld, lmfshal, lmfdeep2, cldcov,                        &
             re_cloud,re_ice,re_snow,                                   &
!             lwp_ex, iwp_ex, lwp_fc, iwp_fc, dzlay,                    &
!             gridkm,                                                    &
!             cldtot, cldcnv,                                            &
!   --- outputs:
!            cld_frac, cld_lwp, cld_reliq, cld_iwp,                     &
!            cld_reice,cld_rwp, cld_rerain,cld_swp, cld_resnow          &
             clouds,clds,mtop,mbot                                      &
           )
! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:    progcld_thompson    computes cloud related quantities  !
!   using Thompson cloud microphysics scheme.                           !
!                                                                       !
! abstract:  this program computes cloud fractions from cloud           !
!   condensates,                                                        !
!   and computes the low, mid, high, total and boundary layer cloud     !
!   fractions and the vertical indices of low, mid, and high cloud      !
!   top and base.  the three vertical cloud domains are set up in the   !
!   initial subroutine "cld_init".                                      !
!                                                                       !
! usage:         call progcld_thompson                                  !
!                                                                       !
! subprograms called:   gethml                                          !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  definition of variables  ====================  !
!                                                                       !
!                                                                       !
! input variables:                                                      !
!   plyr  (IX,NLAY) : model layer mean pressure in mb (100Pa)           !
!   plvl  (IX,NLP1) : model level pressure in mb (100Pa)                !
!   tlyr  (IX,NLAY) : model layer mean temperature in k                 !
!   tvly  (IX,NLAY) : model layer virtual temperature in k              !
!   qlyr  (IX,NLAY) : layer specific humidity in gm/gm                  !
!   qstl  (IX,NLAY) : layer saturate humidity in gm/gm                  !
!   rhly  (IX,NLAY) : layer relative humidity (=qlyr/qstl)              !
!   clw   (IX,NLAY,ntrac) : layer cloud condensate amount               !
!   xlat  (IX)      : grid latitude in radians, default to pi/2 -> -pi/2!
!                     range, otherwise see in-line comment              !
!   xlon  (IX)      : grid longitude in radians  (not used)             !
!   slmsk (IX)      : sea/land mask array (sea:0,land:1,sea-ice:2)      !
!   dz    (ix,nlay) : layer thickness (km)                              !
!   delp  (ix,nlay) : model layer pressure thickness in mb (100Pa)      !
!   gridkm (IX)     : grid length in km                                 !
!   IX              : horizontal dimention                              !
!   NLAY,NLP1       : vertical layer/level dimensions                   !
!   uni_cld         : logical - true for cloud fraction from shoc       !
!   lmfshal         : logical - true for mass flux shallow convection   !
!   lmfdeep2        : logical - true for mass flux deep convection      !
!   cldcov          : layer cloud fraction (used when uni_cld=.true.    !
!                                                                       !
! output variables:                                                     !
!   cloud profiles:                                                     !
!      cld_frac  (:,:) - layer total cloud fraction                     !
!      cld_lwp   (:,:) - layer cloud liq water path       (g/m**2)      !
!      cld_reliq (:,:) - mean eff radius for liq cloud    (micron)      !
!      cld_iwp   (:,:) - layer cloud ice water path       (g/m**2)      !
!      cld_reice (:,:) - mean eff radius for ice cloud    (micron)      !
!      cld_rwp   (:,:) - layer rain drop water path       not assigned  !
!      cld_rerain(:,:) - mean eff radius for rain drop    (micron)      !
!  *** cld_swp   (:,:) - layer snow flake water path      not assigned  !
!      cld_resnow(:,:) - mean eff radius for snow flake   (micron)      !
!                                                                       !
! module variables:                                                     !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   lmfshal         : mass-flux shallow conv scheme flag                !
!   lmfdeep2        : scale-aware mass-flux deep conv scheme flag       !
!   lcrick          : control flag for eliminating CRICK                !
!                     =t: apply layer smoothing to eliminate CRICK      !
!                     =f: do not apply layer smoothing                  !
!   lcnorm          : control flag for in-cld condensate                !
!                     =t: normalize cloud condensate                    !
!                     =f: not normalize cloud condensate                !
!                                                                       !
!  ====================    end of description    =====================  !
!
      use module_mp_thompson_new, only : cfflag_thom

      implicit none
!  ---  inputs
      integer,  intent(in) :: IX, NLAY, NLP1
      integer,  intent(in) :: ntrac, ntcw, ntiw, ntrw, ntsw, ntgl
      logical, intent(in)  :: uni_cld, lmfshal, lmfdeep2
      real (kind=kind_phys), dimension(:,:), intent(in) :: plvl, plyr,  &
     &       tlyr, qlyr, qstl, rhly, cldcov, re_cloud, re_ice, re_snow
!      real (kind=kind_phys), dimension(:,:), intent(in) :: delp, dz
!      real (kind=kind_phys), dimension(:,:), intent(in) :: dzlay
!      real (kind=kind_phys), dimension(:), intent(inout) ::             &
!     &       lwp_ex, iwp_ex, lwp_fc, iwp_fc
      real (kind=kind_phys), dimension(:,:,:), intent(in) :: clw
      real (kind=kind_phys), dimension(:),   intent(in) :: xlat, xlon,  &
     &       slmsk
!      real(kind=kind_phys), dimension(:), intent(in) :: gridkm
!  --- inputs/outputs
!      real (kind=kind_phys), dimension(:,:), intent(inout) ::            &
!     &   cld_frac, cld_lwp, cld_reliq, cld_iwp, cld_reice,               &
!     &   cld_rwp, cld_rerain, cld_swp, cld_resnow
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: clouds
      real (kind=kind_phys), dimension(:,:),   intent(out) :: clds
      integer,               dimension(:,:),   intent(out) :: mtop,mbot
!  ---  local variables:
      real (kind=kind_phys), dimension(IX,NLAY) :: cldtot, cldcnv,      &
     &       cwp, cip, crp, csp, rew, rei, res, rer, clwf
      real (kind=kind_phys), dimension(IX,NLAY) :: dz, delp
      real (kind=kind_phys), dimension(NLAY) :: cldfra1d, qv1d,         &
     &                                 qc1d, qi1d, qs1d, dz1d, p1d, t1d,&
     &                                 rh1d, qst1d
      real(kind=kind_phys) :: gridkm
      real (kind=kind_phys), dimension(IX,NK_CLDS+1) :: ptop1
      real (kind=kind_phys) :: clwmin, tem1, tem2
      real (kind=kind_phys) :: corr, xland, snow_mass_factor
      real (kind=kind_phys), parameter :: max_relh = 1.5
      real (kind=kind_phys), parameter :: snow_max_radius = 130.0
      integer :: i, k, k2, id, nf, idx_rei
!
!===> ... begin here
!
      clwmin = 1.0E-9
      do k = 1, NLAY
        do i = 1, IX
          cldtot(i,k) = 0.0
          cldcnv(i,k) = 0.0
          cwp   (i,k) = 0.0
          cip   (i,k) = 0.0
          crp   (i,k) = 0.0
          csp   (i,k) = 0.0
          rew   (i,k) = re_cloud(i,k)
          rei   (i,k) = re_ice(i,k)
          rer   (i,k) = rrain_def            ! default rain radius to 1000 micron
          res   (i,k) = re_snow(i,K)
          clwf  (i,k) = clw(i,k,ntcw) + clw(i,k,ntiw) + clw(i,k,ntsw)
        enddo
      enddo
!> - Compute cloud liquid/ice condensate path in \f$ g/m^2 \f$ .
!> - Since using Thompson MP, assume 1 percent of snow is actually in
!!   ice sizes.
      if ( ivflip == 0 ) then          ! input data from toa to sfc
          do k = 1, NLAY
            do i = 1, IX
              delp(i,k) = plvl(i,k+1) - plvl(i,k)
            enddo
          enddo
      else                             ! input data from sfc to toa
          do k = 1, NLAY
            do i = 1, IX
              delp(i,k) = plvl(i,k) - plvl(i,k+1)
            enddo
          enddo
      endif

      do k = 1, NLAY-1
        do i = 1, IX
          cwp(i,k) = max(0.0, clw(i,k,ntcw) * gfac * delp(i,k))
          crp(i,k) = 0.0
          snow_mass_factor = 0.90
          cip(i,k) = max(0.0, (clw(i,k,ntiw)                            &
     &            +(1.0-snow_mass_factor)*clw(i,k,ntsw))*gfac*delp(i,k))
          if (re_snow(i,k) .gt. snow_max_radius)then
             snow_mass_factor = min(snow_mass_factor,                   &
     &                              (snow_max_radius/re_snow(i,k))      &
     &                             *(snow_max_radius/re_snow(i,k)))
             res(i,k) = snow_max_radius
          endif
          csp(i,k) = max(0.,snow_mass_factor*clw(i,k,ntsw)*gfac*delp(i,k))
        enddo
      enddo

!> - Sum the liquid water and ice paths that come from explicit micro
!      do i = 1, IX
!         lwp_ex(i) = 0.0
!         iwp_ex(i) = 0.0
!         do k = 1, NLAY-1
!            lwp_ex(i) = lwp_ex(i) + cwp(i,k)
!            iwp_ex(i) = iwp_ex(i) + cip(i,k) + csp(i,k)
!         enddo
!      enddo

!> - Now determine the cloud fraction.  Here, we will use the scheme of
!!   G. Thompson that implements a variannt of Mocko and Cotton (1995)
!!   based on work within HWRF and WRF.  Where the bulk microphysics
!!   scheme already has explicit clouds, assume cloud fraction of one,
!!   but, otherwise, use a Sundqvist et al (1989) scheme and RH-critical
!!   to account for sub-grid-scale clouds, include those in the water
!!   and ice paths _seen_ by the radiation scheme, but do not actually
!!   include these fake clouds into anything other than radiation.
      if ( cfflag_thom .eq. 1 ) then
        call cloud_fraction_XuRandall                                   &
               ( IX, NLAY, plyr, clwf, rhly, qstl,                      &
                 lmfshal, lmfdeep2, cldtot )
      elseif ( cfflag_thom .eq. 2 ) then
!        do i = 1, IX
!          if (slmsk(i)-0.5 .gt. 0.5 .and. slmsk(i)+0.5 .lt. 1.5) then
!             xland = 1.0
!          else
!             xland = 2.0
!          endif
!          if (ivflip .eq. 1) then
!            do k = 1, NLAY
!               qv1d(k) = qlyr(i,k)
!               qc1d(k) = max(0.0, clw(i,k,ntcw))
!               qi1d(k) = max(0.0, clw(i,k,ntiw))
!               qs1d(k) = max(0.0, clw(i,k,ntsw))
!               rh1d(k) = rhly(i,k)
!               qst1d(k)= qstl(i,k)
!               p1d(k) = plyr(i,k)*100.0
!               t1d(k) = tlyr(i,k)
!            enddo
!          else
!            do k = NLAY, 1, -1
!               k2 = NLAY - k + 1
!               qv1d(k2) = qlyr(i,k)
!               qc1d(k2) = max(0.0, clw(i,k,ntcw))
!               qi1d(k2) = max(0.0, clw(i,k,ntiw))
!               qs1d(k2) = max(0.0, clw(i,k,ntsw))
!               rh1d(k2) = rhly(i,k)
!               qst1d(k2)= qstl(i,k)
!               p1d(k2) = plyr(i,k)*100.0
!               t1d(k2) = tlyr(i,k)
!            enddo
!          endif
!          call cal_cldfra3(cldfra1d, qv1d, qc1d, qi1d, qs1d, dz1d,      &
!     &                     p1d, t1d, xland, gridkm(i),                  &
!     &                     .false., max_relh, 1, nlay, .false.)
!          do k = 1, NLAY
!            cldtot(i,k) = cldfra1d(k)
!            if (qc1d(k).gt.clwmin .and. cldfra1d(k).lt.ovcst) then
!               cwp(i,k) = qc1d(k) * gfac * delp(i,k)
!               if ((xland-1.5).GT.0.) then
!                  rew(i,k) = 9.5   ! over ocean
!               else
!                  rew(i,k) = 5.5   ! over land
!               endif
!            endif
!            if (qi1d(k).gt.clwmin .and. cldfra1d(k).lt.ovcst) then
!               cip(i,k) = qi1d(k) * gfac * delp(i,k)
!               idx_rei = int(t1d(k)-179.)
!               idx_rei = min(max(idx_rei,1),75)
!               corr = t1d(k) - int(t1d(k))
!               rei(i,K) = max(5.0, retab(idx_rei)*(1.-corr) +           &
!     &                                 retab(idx_rei+1)*corr)
!            endif
!          enddo
!        enddo
        do i = 1, IX
          do k = 1, NLAY
            cldtot(i,k) = cldcov(i,k)
          enddo
        enddo
      endif

!  ---  find top pressure for each cloud domain for given latitude
!       ptopc(k,i): top presure of each cld domain (k=1-4 are sfc,L,m,h;
!  ---  i=1,2 are low-lat (<45 degree) and pole regions)

      do id = 1, 4
        tem1 = ptopc(id,2) - ptopc(id,1)

        do i =1, IX
          tem2 = xlat(i) / con_pi        ! if xlat in pi/2 -> -pi/2 range
!         tem2 = 0.5 - xlat(i)/con_pi    ! if xlat in 0 -> pi range

          ptop1(i,id) = ptopc(id,1) + tem1*max( 0.0, 4.0*abs(tem2)-1.0 )
        enddo
      enddo

      do k = 1, NLAY
        do i = 1, IX
!          cld_frac(i,k)   = cldtot(i,k)
!          cld_lwp(i,k)    = cwp(i,k)
!          cld_reliq(i,k)  = rew(i,k)
!          cld_iwp(i,k)    = cip(i,k)
!          cld_reice(i,k)  = rei(i,k)
!          cld_rwp(i,k)    = crp(i,k)  ! added for Thompson
!          cld_rerain(i,k) = rer(i,k)
!          cld_swp(i,k)    = csp(i,k)  ! added for Thompson
!          cld_resnow(i,k) = res(i,k)
          clouds(i,k,1) = cldtot(i,k)
          clouds(i,k,2) = cwp(i,k)
          clouds(i,k,3) = rew(i,k)
          clouds(i,k,4) = cip(i,k)
          clouds(i,k,5) = rei(i,k)
          clouds(i,k,6) = crp(i,k)
          clouds(i,k,7) = rer(i,k)
          clouds(i,k,8) = csp(i,k)
          clouds(i,k,9) = res(i,k)
        enddo
      enddo

!  ---  compute low, mid, high, total, and boundary layer cloud fractions
!       and clouds top/bottom layer indices for low, mid, and high clouds.
!       The three cloud domain boundaries are defined by ptopc.  The cloud
!       overlapping method is defined by control flag 'iovr', which may
!       be different for lw and sw radiation programs.

      call gethml                                                       &
!  ---  inputs:
     &     ( plyr, ptop1, cldtot, cldcnv,                               &
     &       IX,NLAY,                                                   &
!  ---  outputs:
     &       clds, mtop, mbot                                           &
     &     )

!> - Sum the liquid water and ice paths that come from fractional clouds
!      do i = 1, IX
!         lwp_fc(i) = 0.0
!         iwp_fc(i) = 0.0
!         do k = 1, NLAY
!            lwp_fc(i) = lwp_fc(i) + cwp(i,k)
!            iwp_fc(i) = iwp_fc(i) + cip(i,k) + csp(i,k)
!         enddo
!         lwp_fc(i) = MAX(0.0, lwp_fc(i) - lwp_ex(i))
!         iwp_fc(i) = MAX(0.0, iwp_fc(i) - iwp_ex(i))
!         lwp_fc(i) = lwp_fc(i)*1.E-3
!         iwp_fc(i) = iwp_fc(i)*1.E-3
!         lwp_ex(i) = lwp_ex(i)*1.E-3
!         iwp_ex(i) = iwp_ex(i)*1.E-3
!      enddo
!
      return
!............................................
      end subroutine progcld_thompson
!............................................

!-----------------------------------
      subroutine progcld_gce                                            &
!...................................
!
!  ---  inputs:
     &     ( plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,                    &
     &       xlat,xlon,slmsk,ntrac,                                     &
     &       effr_cw,effr_iw,effr_sw,effr_rw,effr_in,                   &
     &       IX, NLAY, NLP1, lmfshal, lmfdeep2,                         &
!  ---  outputs:
     &       clouds,clds,mtop,mbot,cldcov                               &
     &      )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:    progcld_gce   computes cloud related quantities using  !
!   Goddard (GCE) 4ICE MP prognostic cloud microphysics scheme. Moist   !
!   species from MP are fed into the corresponding arrays for           !
!   calcuation of     !
!                                                                       !
! abstract:  this program computes cloud fractions from cloud           !
!   condensates, calculates liquid/ice cloud droplet effective radius,  !
!   and computes the low, mid, high, total and boundary layer cloud     !
!   fractions and the vertical indices of low, mid, and high cloud      !
!   top and base.  the three vertical cloud domains are set up in the   !
!   initial subroutine "cld_init".                                      !
!                                                                       !
! usage:         call progcld_gce                                       !
!                                                                       !
! subprograms called:   gethml                                          !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  definition of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (IX,NLAY) : model layer mean pressure in mb (100Pa)           !
!   plvl  (IX,NLP1) : model level pressure in mb (100Pa)                !
!   tlyr  (IX,NLAY) : model layer mean temperature in k                 !
!   tvly  (IX,NLAY) : model layer virtual temperature in k              !
!   qlyr  (IX,NLAY) : layer specific humidity in gm/gm                  !
!   qstl  (IX,NLAY) : layer saturate humidity in gm/gm                  !
!   rhly  (IX,NLAY) : layer relative humidity (=qlyr/qstl)              !
!   clw   (IX,NLAY,NTRAC) : layer cloud condensate amount               !
!   xlat  (IX)      : grid latitude in radians, default to pi/2 -> -pi/2!
!                     range, otherwise see in-line comment              !
!   xlon  (IX)      : grid longitude in radians  (not used)             !
!   slmsk (IX)      : sea/land mask array (sea:0,land:1,sea-ice:2)      !
!   ntrac           : number of tracers
!   ntcw            : index for liquid water
!   ntiw            : index for ice water
!   ntrw            : index for rain water
!   ntsw            : index for snow water
!   ntgl            : index for graupel
!   nthl            : index for hail
!   cldcov(IX,NLAY) : total cloud fraction
!   effr_cw(IX,NLAY): effective radius of liquid water (mm)             !
!   effr_iw(IX,NLAY): effective radius of ice water (mm)                !
!   effr_rw(IX,NLAY): effective radius of rain water (mm)               !
!   effr_sw(IX,NLAY): effective radius of snow water (mm)               !
!   IX              : horizontal dimention                              !
!   NLAY,NLP1       : vertical layer/level dimensions                   !
!                                                                       !
! output variables:                                                     !
!   clouds(IX,NLAY,NF_CLDS) : cloud profiles                            !
!      clouds(:,:,1) - layer total cloud fraction                       !
!      clouds(:,:,2) - layer cloud liq water path         (g/m**2)      !
!      clouds(:,:,3) - mean eff radius for liq cloud      (micron)      !
!      clouds(:,:,4) - layer cloud ice water path         (g/m**2)      !
!      clouds(:,:,5) - mean eff radius for ice cloud      (micron)      !
!      clouds(:,:,6) - layer rain drop water path         not assigned  !
!      clouds(:,:,7) - mean eff radius for rain drop      (micron)      !
!  *** clouds(:,:,8) - layer snow flake water path        not assigned  !
!      clouds(:,:,9) - mean eff radius for snow flake     (micron)      !
!  *** fu's scheme need to be normalized by snow density (g/m**3/1.0e6) !
!   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !
!                                                                       !
! module variables:                                                     !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   lsashal         : control flag for shallow convection               !
!   lcrick          : control flag for eliminating CRICK                !
!                     =t: apply layer smoothing to eliminate CRICK      !
!                     =f: do not apply layer smoothing                  !
!   lcnorm          : control flag for in-cld condensate                !
!                     =t: normalize cloud condensate                    !
!                     =f: not normalize cloud condensate                !
!                                                                       !
!  ====================    end of description    =====================  !
!
      use radn, only : ntcw, ntiw, ntrw, ntsw, ntgl, nthl
      implicit none

!  ---  inputs
      integer,  intent(in) :: IX, NLAY, NLP1
      integer,  intent(in) :: ntrac

      real (kind=kind_phys), dimension(:,:), intent(in) :: plvl, plyr,  &
     &       tlyr, tvly, qlyr, qstl, rhly


      real (kind=kind_phys), dimension(:,:,:), intent(in) :: clw
      real (kind=kind_phys), dimension(:),   intent(in) :: xlat, xlon,  &
     &       slmsk

      real (kind=kind_phys), dimension(:,:), intent(in) :: effr_cw,     &
     &       effr_iw, effr_rw, effr_sw

      logical, intent(in) :: effr_in, lmfshal, lmfdeep2

!  ---  outputs
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: clouds

      real (kind=kind_phys), dimension(:,:),   intent(out) :: clds
      real (kind=kind_phys), dimension(:,:),   intent(out) :: cldcov

      integer,               dimension(:,:),   intent(out) :: mtop,mbot

!  ---  local variables:
      real (kind=kind_phys), dimension(IX,NLAY) :: cldcnv,              &
     &       cwp, cip, crp, csp, rew, rei, res, rer, delp, tem2d, clwf

      real (kind=kind_phys) :: ptop1(IX,NK_CLDS+1)

      real (kind=kind_phys) :: clwmin, clwm, clwt, onemrh, value,       &
     &       tem1, tem2, tem3

      integer :: i, k, id, nf

!
!===> ... begin here
!
      do nf=1,nf_clds
        do k=1,nlay
          do i=1,ix
            clouds(i,k,nf) = 0.0
          enddo
        enddo
      enddo
!     clouds(:,:,:) = 0.0

      do k = 1, NLAY
        do i = 1, IX
          cldcnv(i,k) = 0.0
          cwp   (i,k) = 0.0
          cip   (i,k) = 0.0
          crp   (i,k) = 0.0
          csp   (i,k) = 0.0
          if ( effr_in ) then
            rew (i,k) = effr_cw (i,k)
            rei (i,k) = effr_iw (i,k)
            rer (i,k) = effr_rw (i,k)
            res (i,k) = effr_sw (i,k)
          else
            rew (i,k) = reliq_def            ! default liq radius to 10 micron
            rei (i,k) = reice_def            ! default ice radius to 50 micron
            rer (i,k) = rrain_def            ! default rain radius to 1000 micron
            res (i,k) = rsnow_def            ! default snow radius to 250 micron
          endif
          tem2d (i,k) = min( 1.0, max( 0.0, (con_ttp-tlyr(i,k))*0.05 ) )
          cldcov(i,k) = 0.0
!          clwf  (i,k) = clw(i,k,ntcw)+clw(i,k,ntiw)+clw(i,k,ntrw)+      &
!                        clw(i,k,ntsw)+clw(i,k,ntgl)+clw(i,k,nthl)
!          clwf  (i,k) = clw(i,k,ntcw)+clw(i,k,ntiw)+clw(i,k,ntsw)
          clwf  (i,k) = clw(i,k,ntcw)+clw(i,k,ntiw)+clw(i,k,ntrw)+      &
                        clw(i,k,ntsw)+clw(i,k,ntgl)
        enddo
      enddo

!  ---  find top pressure for each cloud domain for given latitude
!       ptopc(k,i): top presure of each cld domain (k=1-4 are sfc,L,m,h;
!  ---  i=1,2 are low-lat (<45 degree) and pole regions)

      do id = 1, 4
        tem1 = ptopc(id,2) - ptopc(id,1)

        do i =1, IX
          tem2 = xlat(i) / con_pi        ! if xlat in pi/2 -> -pi/2 range
!         tem2 = 0.5 - xlat(i)/con_pi    ! if xlat in 0 -> pi range

          ptop1(i,id) = ptopc(id,1) + tem1*max( 0.0, 4.0*abs(tem2)-1.0 )
        enddo
      enddo

!  ---  compute liquid/ice condensate path in g/m**2

      if ( ivflip == 0 ) then          ! input data from toa to sfc
        do k = 1, NLAY
          do i = 1, IX
            delp(i,k) = plvl(i,k+1) - plvl(i,k)
            cwp(i,k) = max(0.0, clw(i,k,ntcw) * gfac * delp(i,k))
            cip(i,k) = max(0.0, clw(i,k,ntiw) * gfac * delp(i,k))
            crp(i,k) = max(0.0, clw(i,k,ntrw) * gfac * delp(i,k))
!            csp(i,k) = max(0.0, ( clw(i,k,ntsw)+clw(i,k,ntgl)+          &
!     &                 clw(i,k,nthl) ) * gfac * delp(i,k))
            csp(i,k) = max(0.0, ( clw(i,k,ntsw)+clw(i,k,ntgl) )         &
     &                  * gfac * delp(i,k))
          enddo
        enddo
      else                             ! input data from sfc to toa
        do k = 1, NLAY
          do i = 1, IX
            delp(i,k) = plvl(i,k) - plvl(i,k+1)
            cwp(i,k) = max(0.0, clw(i,k,ntcw) * gfac * delp(i,k))
            cip(i,k) = max(0.0, clw(i,k,ntiw) * gfac * delp(i,k))
            crp(i,k) = max(0.0, clw(i,k,ntrw) * gfac * delp(i,k))
            csp(i,k) = max(0.0, ( clw(i,k,ntsw)+clw(i,k,ntgl) )         &
     &                  * gfac * delp(i,k))
          enddo
        enddo
      endif                            ! end_if_ivflip

!  --- calculate cloud fraction :
      call cloud_fraction_XuRandall                                     &
               ( IX, NLAY, plyr, clwf, rhly, qstl,                      &
                 lmfshal, lmfdeep2, cldcov )

      do k = 1, NLAY
        do i = 1, IX
          ! impose minimum cloudiness if substantial liquid exists
          tem1 = clw(i,k,ntcw) + clw(i,k,ntrw)
          if ( tem1 .ge. 1.e-8 ) then
            cldcov(i,k) = max(0.05,cldcov(i,k))
          endif
          ! saturation cloudiness
          if ( rhly(i,k).gt.0.99 .and. tlyr(i,k).gt.288.16 ) then
            cldcov(i,k) = 1.0
          endif
        enddo
      enddo

!  ---  find top pressure for each cloud domain for given latitude
!       ptopc(k,i): top presure of each cld domain (k=1-4 are sfc,L,m,h;
!  ---  i=1,2 are low-lat (<45 degree) and pole regions)

      do id = 1, 4
        tem1 = ptopc(id,2) - ptopc(id,1)

        do i =1, IX
          tem2 = xlat(i) / con_pi        ! if xlat in pi/2 -> -pi/2 range
!         tem2 = 0.5 - xlat(i)/con_pi    ! if xlat in 0 -> pi range

          ptop1(i,id) = ptopc(id,1) + tem1*max( 0.0, 4.0*abs(tem2)-1.0 )
        enddo
      enddo

!  ---  effective liquid cloud droplet radius over land

      if ( .not. effr_in ) then
      do i = 1, IX
        if (nint(slmsk(i)) == 1) then
          do k = 1, NLAY
            rew(i,k) = 5.0 + 5.0 * tem2d(i,k)
          enddo
        endif
      enddo
      endif

      do k = 1, NLAY
        do i = 1, IX
          if (cldcov(i,k) < climit) then
            cldcov(i,k) = 0.0 
            cwp(i,k)    = 0.0
            cip(i,k)    = 0.0
            crp(i,k)    = 0.0
            csp(i,k)    = 0.0
          endif
        enddo
      enddo

      if ( lcnorm ) then
        do k = 1, NLAY
          do i = 1, IX
            if (cldcov(i,k) >= climit) then
              tem1 = 1.0 / max(climit2, cldcov(i,k))
              cwp(i,k) = cwp(i,k) * tem1
              cip(i,k) = cip(i,k) * tem1
              crp(i,k) = crp(i,k) * tem1
              csp(i,k) = csp(i,k) * tem1
            endif
          enddo
        enddo
      endif

!  ---  effective ice cloud droplet radius

      if ( .not. effr_in ) then
      do k = 1, NLAY
        do i = 1, IX
          tem2 = tlyr(i,k) - con_ttp

          if (cip(i,k) > 0.0) then
            tem3 = gord * cip(i,k) * plyr(i,k) / (delp(i,k)*tvly(i,k))

            if (tem2 < -50.0) then
              rei(i,k) = (1250.0/9.917) * tem3 ** 0.109
            elseif (tem2 < -40.0) then
              rei(i,k) = (1250.0/9.337) * tem3 ** 0.08
            elseif (tem2 < -30.0) then
              rei(i,k) = (1250.0/9.208) * tem3 ** 0.055
            else
              rei(i,k) = (1250.0/9.387) * tem3 ** 0.031
            endif
!           rei(i,k)   = max(20.0, min(rei(i,k), 300.0))
!           rei(i,k)   = max(10.0, min(rei(i,k), 100.0))
            rei(i,k)   = max(10.0, min(rei(i,k), 150.0))
!           rei(i,k)   = max(5.0,  min(rei(i,k), 130.0))
          endif
        enddo
      enddo
      endif

!
      do k = 1, NLAY
        do i = 1, IX
          clouds(i,k,1) = cldcov(i,k)
          clouds(i,k,2) = cwp(i,k)
          clouds(i,k,3) = rew(i,k)
          clouds(i,k,4) = cip(i,k)
          clouds(i,k,5) = rei(i,k)
          clouds(i,k,6) = crp(i,k) 
          clouds(i,k,7) = rer(i,k)
          clouds(i,k,8) = csp(i,k)
          clouds(i,k,9) = res(i,k)
        enddo
      enddo


!  ---  compute low, mid, high, total, and boundary layer cloud fractions
!       and clouds top/bottom layer indices for low, mid, and high clouds.
!       The three cloud domain boundaries are defined by ptopc.  The cloud
!       overlapping method is defined by control flag 'iovr', which may
!       be different for lw and sw radiation programs.

      call gethml                                                       &
!  ---  inputs:
     &     ( plyr, ptop1, cldcov, cldcnv,                               &
     &       IX,NLAY,                                                   &
!  ---  outputs:
     &       clds, mtop, mbot                                           &
     &     )


!
      return
!...................................
      end subroutine progcld_gce
!-----------------------------------

!-----------------------------------
      subroutine progclduni                                             &
!  ---  inputs:
     &     ( plyr,plvl,tlyr,tvly,ccnd,ncnd,                             &
     &       xlat,xlon,slmsk, IX, NLAY, NLP1, cldcov,                   &
     &       effrl,effri,effrs,effrr,effr_in,                           &
!  ---  outputs:
     &       clouds,clds,mtop,mbot                                      &
     &      )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:    progclduni    computes cloud related quantities using    !
!   for unified cloud microphysics scheme.                !
!                                                                       !
! abstract:  this program computes cloud fractions from cloud           !
!   condensates, calculates liquid/ice cloud droplet effective radius,  !
!   and computes the low, mid, high, total and boundary layer cloud     !
!   fractions and the vertical indices of low, mid, and high cloud      !
!   top and base.  the three vertical cloud domains are set up in the   !
!   initial subroutine "cld_init".                                      !
!                                                                       !
! usage:         call progclduni                                          !
!                                                                       !
! subprograms called:   gethml                                          !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  definition of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (IX,NLAY) : model layer mean pressure in mb (100Pa)           !
!   plvl  (IX,NLP1) : model level pressure in mb (100Pa)                !
!   tlyr  (IX,NLAY) : model layer mean temperature in k                 !
!   tvly  (IX,NLAY) : model layer virtual temperature in k              !
!   ccnd  (IX,NLAY,ncnd) : layer cloud condensate amount                     !
!   ncnd            : number of layer cloud condensate types            !
!   xlat  (IX)      : grid latitude in radians, default to pi/2 -> -pi/2!
!                     range, otherwise see in-line comment              !
!   xlon  (IX)      : grid longitude in radians  (not used)             !
!   slmsk (IX)      : sea/land mask array (sea:0,land:1,sea-ice:2)      !
!   IX              : horizontal dimention                              !
!   NLAY,NLP1       : vertical layer/level dimensions                   !
!   cldcov          : unified cloud fracrion from moist physics         !
!   effrl (ix,nlay) : effective radius for liquid water                 !
!   effri (ix,nlay) : effective radius for ice water                    !
!   effrr (ix,nlay) : effective radius for rain water                   !
!   effrs (ix,nlay) : effective radius for snow water                   !
!   effr_in         : logical - if .true. use input effective radii     !
!                                                                       !
! output variables:                                                     !
!   clouds(IX,NLAY,NF_CLDS) : cloud profiles                            !
!      clouds(:,:,1) - layer total cloud fraction                       !
!      clouds(:,:,2) - layer cloud liq water path         (g/m**2)      !
!      clouds(:,:,3) - mean eff radius for liq cloud      (micron)      !
!      clouds(:,:,4) - layer cloud ice water path         (g/m**2)      !
!      clouds(:,:,5) - mean eff radius for ice cloud      (micron)      !
!      clouds(:,:,6) - layer rain drop water path         not assigned  !
!      clouds(:,:,7) - mean eff radius for rain drop      (micron)      !
!  *** clouds(:,:,8) - layer snow flake water path        not assigned  !
!      clouds(:,:,9) - mean eff radius for snow flake     (micron)      !
!  *** fu's scheme need to be normalized by snow density (g/m**3/1.0e6) !
!   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !
!                                                                       !
! module variables:                                                     !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   lmfshal         : mass-flux shallow conv scheme flag                !
!   lmfdeep2        : scale-aware mass-flux deep conv scheme flag       !
!   lcrick          : control flag for eliminating CRICK                !
!                     =t: apply layer smoothing to eliminate CRICK      !
!                     =f: do not apply layer smoothing                  !
!   lcnorm          : control flag for in-cld condensate                !
!                     =t: normalize cloud condensate                    !
!                     =f: not normalize cloud condensate                !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs
      integer,  intent(in) :: IX, NLAY, NLP1, ncnd
      logical,  intent(in) :: effr_in

      real (kind=kind_phys), dimension(:,:,:), intent(in) :: ccnd
      real (kind=kind_phys), dimension(:,:),   intent(in) :: plvl, plyr,&
     &       tlyr, tvly, cldcov, effrl, effri, effrr, effrs

      real (kind=kind_phys), dimension(:),   intent(in) :: xlat, xlon,  &
     &       slmsk

!  ---  outputs
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: clouds

      real (kind=kind_phys), dimension(:,:),   intent(out) :: clds

      integer,               dimension(:,:),   intent(out) :: mtop,mbot

!  ---  local variables:
      real (kind=kind_phys), dimension(IX,NLAY) :: cldcnv, cwp, cip,    &
     &       crp, csp, rew, rei, res, rer, delp, tem2d, cldtot
      real (kind=kind_phys), dimension(IX,NLAY,ncnd) :: cndf

      real (kind=kind_phys) :: ptop1(IX,NK_CLDS+1)

      real (kind=kind_phys) :: tem1, tem2, tem3

      integer :: i, k, id, nf, n

!
!===> ... begin here
!
      do nf=1,nf_clds
        do k=1,nlay
          do i=1,ix
            clouds(i,k,nf) = 0.0
          enddo
        enddo
      enddo
!     clouds(:,:,:) = 0.0

      if (effr_in) then
        do k = 1, NLAY
          do i = 1, IX
            cldtot(i,k) = cldcov(i,k)
            cldcnv(i,k) = 0.0
            cwp   (i,k) = 0.0
            cip   (i,k) = 0.0
            crp   (i,k) = 0.0
            csp   (i,k) = 0.0
            rew   (i,k) = effrl (i,k)
            rei   (i,k) = effri (i,k)
            rer   (i,k) = effrr (i,k)
            res   (i,k) = effrs (i,k)
            tem2d (i,k) = min( 1.0, max( 0.0,(con_ttp-tlyr(i,k))*0.05))
          enddo
        enddo
      else
        do k = 1, NLAY
          do i = 1, IX
            cldcnv(i,k) = 0.0
            cwp   (i,k) = 0.0
            cip   (i,k) = 0.0
            crp   (i,k) = 0.0
            csp   (i,k) = 0.0
            rew   (i,k) = reliq_def            ! default liq radius to 10 micron
            rei   (i,k) = reice_def            ! default ice radius to 50 micron
            rer   (i,k) = rrain_def            ! default rain radius to 1000 micron
            res   (i,k) = rsnow_def            ! default snow radius to 250 micron
            tem2d (i,k) = min(1.0, max(0.0, (con_ttp-tlyr(i,k))*0.05))
          enddo
        enddo
      endif
!
      do n=1,ncnd
        do k = 1, NLAY
          do i = 1, IX
            cndf(i,k,n) = ccnd(i,k,n)
          enddo
        enddo
      enddo
      if ( lcrick ) then
        do n=1,ncnd
          do i = 1, IX
            cndf(i,1,n)    = 0.75*ccnd(i,1,n)    + 0.25*ccnd(i,2,n)
            cndf(i,nlay,n) = 0.75*ccnd(i,nlay,n) + 0.25*ccnd(i,nlay-1,n)
          enddo
          do k = 2, NLAY-1
            do i = 1, IX
              cndf(i,K,n) = 0.25 * (ccnd(i,k-1,n) + ccnd(i,k+1,n))      &
     &                    + 0.5  *  ccnd(i,k,n)
            enddo
          enddo
        enddo
      endif

!> -# Find top pressure for each cloud domain for given latitude.
!     ptopc(k,i): top presure of each cld domain (k=1-4 are sfc,L,m,h;
!  ---  i=1,2 are low-lat (<45 degree) and pole regions)

      do id = 1, 4
        tem1 = ptopc(id,2) - ptopc(id,1)

        do i =1, IX
          tem2 = xlat(i) / con_pi        ! if xlat in pi/2 -> -pi/2 range
!         tem2 = 0.5 - xlat(i)/con_pi    ! if xlat in 0 -> pi range

          ptop1(i,id) = ptopc(id,1) + tem1*max( 0.0, 4.0*abs(tem2)-1.0 )
        enddo
      enddo

!> -# Compute cloud liquid/ice condensate path in \f$ g/m^2 \f$ .

      if ( ivflip == 0 ) then          ! input data from toa to sfc
        if (ncnd == 2) then
          do k = 1, NLAY
            do i = 1, IX
              delp(i,k) = plvl(i,k+1) - plvl(i,k)
              tem1      = gfac * delp(i,k)
              cwp(i,k)  = cndf(i,k,1) * tem1
              cip(i,k)  = cndf(i,k,2) * tem1
            enddo
          enddo
        elseif (ncnd == 4 .or. ncnd == 5) then
          do k = 1, NLAY
            do i = 1, IX
              delp(i,k) = plvl(i,k+1) - plvl(i,k)
              tem1      = gfac * delp(i,k)
              cwp(i,k)  = cndf(i,k,1) * tem1
              cip(i,k)  = cndf(i,k,2) * tem1
              crp(i,k)  = cndf(i,k,3) * tem1
              csp(i,k)  = cndf(i,k,4) * tem1
            enddo
          enddo
        endif
      else                             ! input data from sfc to toa
        if (ncnd == 2) then
          do k = 1, NLAY
            do i = 1, IX
              delp(i,k) = plvl(i,k) - plvl(i,k+1)
              tem1      = gfac * delp(i,k)
              cwp(i,k)  = cndf(i,k,1) * tem1
              cip(i,k)  = cndf(i,k,2) * tem1
            enddo
          enddo
        elseif (ncnd == 4 .or. ncnd == 5) then
          do k = 1, NLAY
            do i = 1, IX
              delp(i,k) = plvl(i,k) - plvl(i,k+1)
              tem1      = gfac * delp(i,k)
              cwp(i,k)  = cndf(i,k,1) * tem1
              cip(i,k)  = cndf(i,k,2) * tem1
              crp(i,k)  = cndf(i,k,3) * tem1
              csp(i,k)  = cndf(i,k,4) * tem1
            enddo
          enddo
        endif

      endif                            ! end_if_ivflip

!> -# Compute effective liquid cloud droplet radius over land.

      if(.not. effr_in) then
        do i = 1, IX
          if (nint(slmsk(i)) == 1) then
            do k = 1, NLAY
              rew(i,k) = 5.0 + 5.0 * tem2d(i,k)
            enddo
          endif
        enddo
      endif

      do k = 1, NLAY
        do i = 1, IX
          if (cldtot(i,k) < climit) then
            cwp(i,k)    = 0.0
            cip(i,k)    = 0.0
            crp(i,k)    = 0.0
            csp(i,k)    = 0.0
          endif
        enddo
      enddo

      if ( lcnorm ) then
        do k = 1, NLAY
          do i = 1, IX
            if (cldtot(i,k) >= climit) then
              tem1 = 1.0 / max(climit2, cldtot(i,k))
              cwp(i,k) = cwp(i,k) * tem1
              cip(i,k) = cip(i,k) * tem1
              crp(i,k) = crp(i,k) * tem1
              csp(i,k) = csp(i,k) * tem1
            endif
          enddo
        enddo
      endif

!> -# Compute effective ice cloud droplet radius following Heymsfield 
!!    and McFarquhar (1996) \cite heymsfield_and_mcfarquhar_1996.

      if(.not. effr_in) then
        do k = 1, NLAY
          do i = 1, IX
            tem2 = tlyr(i,k) - con_ttp

            if (cip(i,k) > 0.0) then
              tem3 = gord * cip(i,k) * plyr(i,k) / (delp(i,k)*tvly(i,k))

              if (tem2 < -50.0) then
                rei(i,k) = (1250.0/9.917) * tem3 ** 0.109
              elseif (tem2 < -40.0) then
                rei(i,k) = (1250.0/9.337) * tem3 ** 0.08
              elseif (tem2 < -30.0) then
                rei(i,k) = (1250.0/9.208) * tem3 ** 0.055
              else
                rei(i,k) = (1250.0/9.387) * tem3 ** 0.031
              endif
!             rei(i,k)   = max(20.0, min(rei(i,k), 300.0))
!             rei(i,k)   = max(10.0, min(rei(i,k), 100.0))
              rei(i,k)   = max(10.0, min(rei(i,k), 150.0))
!             rei(i,k)   = max(5.0,  min(rei(i,k), 130.0))
            endif
          enddo
        enddo
      endif

!
      do k = 1, NLAY
        do i = 1, IX
          clouds(i,k,1) = cldtot(i,k)
          clouds(i,k,2) = cwp(i,k)
          clouds(i,k,3) = rew(i,k)
          clouds(i,k,4) = cip(i,k)
          clouds(i,k,5) = rei(i,k)
          clouds(i,k,6) = crp(i,k)
          clouds(i,k,7) = rer(i,k)
          clouds(i,k,8) = csp(i,k)
          clouds(i,k,9) = res(i,k)
        enddo
      enddo


!> -# Call gethml() to compute low,mid,high,total, and boundary layer
!!    cloud fractions and clouds top/bottom layer indices for low, mid,
!!    and high clouds.
!  ---  compute low, mid, high, total, and boundary layer cloud fractions
!       and clouds top/bottom layer indices for low, mid, and high clouds.
!       The three cloud domain boundaries are defined by ptopc.  The cloud
!       overlapping method is defined by control flag 'iovr', which may
!       be different for lw and sw radiation programs.

      call gethml                                                       &
!  ---  inputs:
     &     ( plyr, ptop1, cldtot, cldcnv,                               &
     &       IX,NLAY,                                                   &
!  ---  outputs:
     &       clds, mtop, mbot                                           &
     &     )


!
      return
!...................................
      end subroutine progclduni
!-----------------------------------


!-----------------------------------
      subroutine diagcld1                                               &
!...................................

!  ---  inputs:
     &     ( plyr,plvl,tlyr,rhly,vvel,cv,cvt,cvb,                       &
     &       xlat,xlon,slmsk,                                           &
     &       ix, nlay, nlp1,                                            &
!  ---  outputs:
     &       clouds,clds,mtop,mbot                                      &
     &      )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:    diagcld1    computes cloud fractions for radiation     !
!   calculations.                                                       !
!                                                                       !
! abstract:  clouds are diagnosed from layer relative humidity, and     !
!   estimate cloud optical depth from temperature and layer thickness.  !
!   then computes the low, mid, high, total and boundary layer cloud    !
!   fractions and the vertical indices of low, mid, and high cloud top  !
!   and base.  the three vertical cloud domains are set up in the       !
!   initial subroutine "cld_init".                                      !
!                                                                       !
! usage:         call diagcld1                                          !
!                                                                       !
! subprograms called:   gethml                                          !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (ix,nlay) : model layer mean pressure in mb (100pa)           !
!   plvl  (ix,nlp1) : model level pressure in mb (100pa)                !
!   tlyr  (ix,nlay) : model layer mean temperature in k                 !
!   rhly  (ix,nlay) : layer relative humidity                           !
!   vvel  (ix,nlay) : layer mean vertical velocity in mb/sec            !
!   clw   (ix,nlay) : layer cloud condensate amount         (not used)  !
!   xlat  (ix)      : grid latitude in radians, default to pi/2 -> -pi/2!
!                     range, otherwise see in-line comment              !
!   xlon  (ix)      : grid longitude in radians, ok for both 0->2pi or  !
!                     -pi -> +pi ranges                                 !
!   slmsk (ix)      : sea/land mask array (sea:0,land:1,sea-ice:2)      !
!   cv    (ix)      : fraction of convective cloud                      !
!   cvt, cvb (ix)   : conv cloud top/bottom pressure in mb              !
!   ix              : horizontal dimention                              !
!   nlay,nlp1       : vertical layer/level dimensions                   !
!                                                                       !
! output variables:                                                     !
!   clouds(ix,nlay,nf_clds) : cloud profiles                            !
!      clouds(:,:,1) - layer total cloud fraction                       !
!      clouds(:,:,2) - layer cloud optical depth                        !
!      clouds(:,:,3) - layer cloud single scattering albedo             !
!      clouds(:,:,4) - layer cloud asymmetry factor                     !
!   clds  (ix,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (ix,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (ix,3)    : vertical indices for low, mid, hi cloud bases     !
!                                                                       !
! external module variables:                                            !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs
      integer,  intent(in) :: ix, nlay, nlp1

      real (kind=kind_phys), dimension(:,:), intent(in) :: plvl, plyr,  &
     &       tlyr, rhly, vvel

      real (kind=kind_phys), dimension(:),   intent(in) :: xlat, xlon,  &
     &       slmsk, cv, cvt, cvb

!  ---  outputs
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: clouds

      real (kind=kind_phys), dimension(:,:),   intent(out) :: clds

      integer,               dimension(:,:),   intent(out) :: mtop,mbot

!  ---  local variables:
      real (kind=kind_phys), dimension(ix,nlay) :: cldtot, cldcnv,      &
     &       cldtau, taufac, dthdp, tem2d

      real (kind=kind_phys) :: ptop1(ix,nk_clds+1)

      real (kind=kind_phys) :: cr1, cr2, crk, pval, cval, omeg, value,  &
     &       tem1, tem2

      integer, dimension(ix):: idom, kcut

!  ---  for rh-cl calculation
      real (kind=kind_phys) :: xlatdg, xlondg, xlnn, xlss, xrgt, xlft,  &
     &       rhcla(nbin,nlon,mcld,nseal),  rhcld(ix,nbin,mcld)

      integer :: ireg, ib, ic, id, id1, il, is, nhalf

      integer :: i, j, k, klowt
!     integer :: klowb

      logical :: notstop

!
!===> ... begin here
!
      clouds(:,:,:) = 0.0

      tem1 = 180.0 / con_pi

      lab_do_i_ix : do i = 1, ix

        xlatdg = xlat(i) * tem1                    ! if xlat in pi/2 -> -pi/2 range
!       xlatdg = 90.0 - xlat(i)*tem1               ! if xlat in 0 -> pi range

        xlondg = xlon(i) * tem1
        if (xlondg < 0.0) xlondg = xlondg + 360.0  ! if in -180->180, chg to 0->360 range

        ireg = 4

!  ---  get rh-cld relation for this lat

        lab_do_j : do j = 1, 3
          if (xlatdg > xlabdy(j)) then
            ireg = j
            exit lab_do_j
          endif
        enddo  lab_do_j

        do is = 1, nseal
          do ic = 1, mcld
            do il = 1, nlon
              do ib = 1, nbin
                rhcla(ib,il,ic,is) = rhcl(ib,il,ireg,ic,is)
              enddo
            enddo
          enddo
        enddo

!  ---  linear transition between latitudinal regions...
        do j = 1, 3
          xlnn = xlabdy(j) + xlim
          xlss = xlabdy(j) - xlim

          if (xlatdg < xlnn .and. xlatdg > xlss) then
            do is = 1, nseal
              do ic = 1, mcld
                do il = 1, nlon
                  do ib = 1, nbin
                    rhcla(ib,il,ic,is) = rhcl(ib,il,j+1,ic,is)          &
     &                + (rhcl(ib,il,j,ic,is)-rhcl(ib,il,j+1,ic,is))     &
     &                * (xlatdg-xlss) / (xlnn-xlss)
                  enddo
                enddo
              enddo
            enddo
          endif

        enddo        ! end_j_loop

!  ---  get rh-cld relationship for each grid point, interpolating
!       longitudinally between regions if necessary..

        if (slmsk(i) < 1.0) then
          is = 2
        else
          is = 1
        endif

!  ---  which hemisphere (e,w)

        if (xlondg > 180.e0) then
          il = 2
        else
          il = 1
        endif

        do ic = 1, mcld
          do ib = 1, nbin
            rhcld(i,ib,ic) = rhcla(ib,il,ic,is)
          enddo

          lab_do_k : do k = 1, 3
            tem2 = abs(xlondg - xlobdy(k))

            if (tem2 < xlim) then
              id = il
              id1= id + 1
              if (id1 > nlon) id1 = 1

              xlft = xlobdy(k) - xlim
              xrgt = xlobdy(k) + xlim

              do ib = 1, nbin
                rhcld(i,ib,ic) = rhcla(ib,id1,ic,is)                    &
     &            + (rhcla(ib,id,ic,is) - rhcla(ib,id1,ic,is))          &
     &            * (xlondg-xrgt)/(xlft-xrgt)
              enddo
              exit lab_do_k
            endif

          enddo  lab_do_k

        enddo   ! end_do_ic_loop
      enddo  lab_do_i_ix

!  ---  find top pressure for each cloud domain

      do j = 1, 4
        tem1 = ptopc(j,2) - ptopc(j,1)

        do i = 1, ix
          tem2 = xlat(i) / con_pi        ! if xlat in pi/2 -> -pi/2 range
!         tem2 = 0.5 - xlat(i)/con_pi    ! if xlat in 0 -> pi range

          ptop1(i,j) = ptopc(j,1) + tem1*max( 0.0, 4.0*abs(tem2)-1.0 )
        enddo
      enddo

!  ---  stratiform cloud optical depth

      do k = 1, nlay
        do i = 1, ix
          tem1 = tlyr(i,k) - con_ttp
          if (tem1 <= -10.0) then
            cldtau(i,k) = max( 0.1e-3, 2.0e-6*(tem1+82.5)**2 )
          else
            cldtau(i,k) = min( 0.08, 6.949e-3*tem1+0.08 )
          endif
        enddo
      enddo

!  ---  potential temperature and its lapse rate

      do k = 1, nlay
        do i = 1, ix
          cldtot(i,k) = 0.0
          cldcnv(i,k) = 0.0
          tem1        = (plyr(i,k)*0.001) ** (-con_rocp)
          tem2d(i,k)  = tem1 * tlyr(i,k)
        enddo
      enddo

      do k = 1, nlay-1
        do i = 1, ix
          dthdp(i,k) = (tem2d(i,k+1)-tem2d(i,k))/(plyr(i,k+1)-plyr(i,k))
        enddo
      enddo
!
!===> ... diagnostic method to find cloud amount cldtot, cldcnv
!

      if ( ivflip == 0 ) then                 ! input data from toa to sfc

!  ---  find the lowest low cloud top sigma level, computed for each lat cause
!       domain definition changes with latitude...

!       klowb = 1
        klowt = 1
        do k = 1, nlay
          do i = 1, ix
!           if (plvl(i,k) < ptop1(i,2))  klowb = k
            if (plvl(i,k) < ptop1(i,2))  klowt = max(klowt,k)
            taufac(i,k) = plvl(i,k+1) - plvl(i,k)
          enddo
        enddo

        do i = 1, ix

!  ---  find the stratosphere cut off layer for high cloud (about 250mb).
!       it is assumed to be above the layer with dthdp less than -0.25 in
!       the high cloud domain

          kcut(i) = 2
          lab_do_kcut0 : do k = klowt-1, 2, -1
            if (plyr(i,k) <= ptop1(i,3) .and.                           &
     &          dthdp(i,k) < -0.25e0) then
              kcut(i) = k
              exit lab_do_kcut0
            endif
          enddo  lab_do_kcut0

!  ---  put convective cloud into 'cldcnv', no merge at this point..

          if (cv(i) >= climit .and. cvt(i) < cvb(i)) then
            id  = nlay
            id1 = nlay

            lab_do_k_cvt0 : do k = 2, nlay
              if (cvt(i) <= plyr(i,k)) then
                id = k - 1
                exit lab_do_k_cvt0
              endif
            enddo  lab_do_k_cvt0

            lab_do_k_cvb0 : do k = nlay-1, 1, -1
              if (cvb(i) >= plyr(i,k)) then
                id1 = k + 1
                exit lab_do_k_cvb0
              endif
            enddo  lab_do_k_cvb0

            tem1 = plyr(i,id1) - plyr(i,id)
            do k = id, id1
              cldcnv(i,k) = cv(i)
              taufac(i,k) = taufac(i,k) * max( 0.25, 1.0-0.125*tem1 )
              cldtau(i,k) = 0.06
            enddo
          endif

        enddo                ! end_do_i_loop

!  ---  calculate stratiform cloud and put into array 'cldtot' using
!       the cloud-rel.humidity relationship from table look-up..where
!       tables obtained using k.mitchell frequency distribution tuning
!bl       (observations are daily means from us af rtneph).....k.a.c.
!bl       tables created without lowest 10 percent of atmos.....k.a.c.
!      (observations are synoptic using -6,+3 window from rtneph)
!       tables are created with lowest 10-percent-of-atmos, and are
!  ---  now used..  25 october 1995 ... kac.

        do k = nlay-1, 2, -1

          if (k < llyr) then
            do i = 1, ix
              idom(i) = 0
            enddo

            do i = 1, ix
              lab_do_ic0 : do ic = 2, 4
                if(plyr(i,k) >= ptop1(i,ic)) then
                  idom(i) = ic
                  exit lab_do_ic0
                endif
              enddo  lab_do_ic0
            enddo
          else
            do i = 1, ix
              idom(i) = 1
            enddo
          endif

          do i = 1, ix
            id = idom(i)
            nhalf = (nbin + 1) / 2

            if (id <= 0 .or. k < kcut(i)) then
              cldtot(i,k) = 0.0
            elseif (rhly(i,k) <= rhcld(i,1,id)) then
              cldtot(i,k) = 0.0
            elseif (rhly(i,k) >= rhcld(i,nbin,id)) then
              cldtot(i,k) = 1.0
            else
              ib = nhalf
              crk = rhly(i,k)

              notstop = .true.
              do while ( notstop )
                nhalf = (nhalf + 1) / 2
                cr1 = rhcld(i,ib,  id)
                cr2 = rhcld(i,ib+1,id)

                if (crk <= cr1) then
                  ib = max( ib-nhalf, 1 )
                elseif (crk > cr2) then
                  ib = min( ib+nhalf, nbin-1 )
                else
                  cldtot(i,k) = 0.01 * (ib + (crk - cr1)/(cr2 - cr1))
                  notstop = .false.
                endif
              enddo      ! end_do_while
            endif
          enddo          ! end_do_i_loop

        enddo            ! end_do_k_loop

! --- vertical velocity adjustment on low clouds

        value = vvcld1 - vvcld2
        do k = klowt, llyr+1
          do i = 1, ix

            omeg = vvel(i,k)
            cval = cldtot(i,k)
            pval = plyr(i,k)

! --- vertical velocity adjustment on low clouds

            if (cval >= climit .and. pval >= ptop1(i,2)) then
              if (omeg >= vvcld1) then
                cldtot(i,k) = 0.0
              elseif (omeg > vvcld2) then
                tem1 = (vvcld1 - omeg) / value
                cldtot(i,k) = cldtot(i,k) * sqrt(tem1)
              endif
            endif

          enddo     ! end_do_i_loop
        enddo       ! end_do_k_loop

      else                                    ! input data from sfc to toa

!  ---  find the lowest low cloud top sigma level, computed for each lat cause
!       domain definition changes with latitude...

!       klowb = nlay
        klowt = nlay
        do k = nlay, 1, -1
          do i = 1, ix
!           if (plvl(i,k) < ptop1(i,2))  klowb = k
            if (plvl(i,k) < ptop1(i,2))  klowt = min(klowt,k)
            taufac(i,k) = plvl(i,k) - plvl(i,k+1)       ! dp for later cal cldtau use
          enddo
        enddo

        do i = 1, ix

!  ---  find the stratosphere cut off layer for high cloud (about 250mb).
!       it is assumed to be above the layer with dthdp less than -0.25 in
!       the high cloud domain

          kcut(i) = nlay - 1
          lab_do_kcut1 : do k = klowt+1, nlay-1
            if (plyr(i,k) <= ptop1(i,3) .and.                           &
     &          dthdp(i,k) < -0.25e0) then
              kcut(i) = k
              exit lab_do_kcut1
            endif
          enddo  lab_do_kcut1

!  ---  put convective cloud into 'cldcnv', no merge at this point..

          if (cv(i) >= climit .and. cvt(i) < cvb(i)) then
            id  = 1
            id1 = 1

            lab_do_k_cvt : do k = nlay-1, 1, -1
              if (cvt(i) <= plyr(i,k)) then
                id = k + 1
                exit lab_do_k_cvt
              endif
            enddo  lab_do_k_cvt

            lab_do_k_cvb : do k = 2, nlay
              if (cvb(i) >= plyr(i,k)) then
                id1 = k - 1
                exit lab_do_k_cvb
              endif
            enddo  lab_do_k_cvb

            tem1 = plyr(i,id1) - plyr(i,id)
            do k = id1, id
              cldcnv(i,k) = cv(i)
              taufac(i,k) = taufac(i,k) * max( 0.25, 1.0-0.125*tem1 )
              cldtau(i,k) = 0.06
            enddo
          endif

        enddo     ! end_do_i_loop

!  ---  calculate stratiform cloud and put into array 'cldtot' using
!       the cloud-rel.humidity relationship from table look-up..where
!       tables obtained using k.mitchell frequency distribution tuning
!bl       (observations are daily means from us af rtneph).....k.a.c.
!bl       tables created without lowest 10 percent of atmos.....k.a.c.
!      (observations are synoptic using -6,+3 window from rtneph)
!       tables are created with lowest 10-percent-of-atmos, and are
!  ---  now used..  25 october 1995 ... kac.

        do k = 2, nlay-1

          if (k > llyr) then
            do i = 1, ix
              idom(i) = 0
            enddo

            do i = 1, ix
              lab_do_ic1 : do ic = 2, 4
                if(plyr(i,k) >= ptop1(i,ic)) then
                  idom(i) = ic
                  exit lab_do_ic1
                endif
              enddo  lab_do_ic1
            enddo
          else
            do i = 1, ix
              idom(i) = 1
            enddo
          endif

          do i = 1, ix
            id = idom(i)
            nhalf = (nbin + 1) / 2

            if (id <= 0 .or. k > kcut(i)) then
              cldtot(i,k) = 0.0
            elseif (rhly(i,k) <= rhcld(i,1,id)) then
              cldtot(i,k) = 0.0
            elseif (rhly(i,k) >= rhcld(i,nbin,id)) then
              cldtot(i,k) = 1.0
            else
              ib = nhalf
              crk = rhly(i,k)

              notstop = .true.
              do while ( notstop )
                nhalf = (nhalf + 1) / 2
                cr1 = rhcld(i,ib,  id)
                cr2 = rhcld(i,ib+1,id)

                if (crk <= cr1) then
                  ib = max( ib-nhalf, 1 )
                elseif (crk > cr2) then
                  ib = min( ib+nhalf, nbin-1 )
                else
                  cldtot(i,k) = 0.01 * (ib + (crk - cr1)/(cr2 - cr1))
                  notstop = .false.
                endif
              enddo      ! end_do_while
            endif
          enddo          ! end_do_i_loop

        enddo            ! end_do_k_loop

! --- vertical velocity adjustment on low clouds

        value = vvcld1 - vvcld2
        do k = llyr-1, klowt
          do i = 1, ix

            omeg = vvel(i,k)
            cval = cldtot(i,k)
            pval = plyr(i,k)

! --- vertical velocity adjustment on low clouds

            if (cval >= climit .and. pval >= ptop1(i,2)) then
              if (omeg >= vvcld1) then
                cldtot(i,k) = 0.0
              elseif (omeg > vvcld2) then
                tem1 = (vvcld1 - omeg) / value
                cldtot(i,k) = cldtot(i,k) * sqrt(tem1)
              endif
            endif

          enddo     ! end_do_i_loop
        enddo       ! end_do_k_loop

      endif                                   ! end_if_ivflip

!  ---  diagnostic cloud optical depth
!     cldtau = cldtau * taufac

      where (cldtot < climit)
        cldtot = 0.0
      endwhere
      where (cldcnv < climit)
        cldcnv = 0.0
      endwhere

      where (cldtot < climit .and. cldcnv < climit)
        cldtau = 0.0
      endwhere

      do k = 1, nlay
        do i = 1, ix
          clouds(i,k,1) = max(cldtot(i,k), cldcnv(i,k))
          clouds(i,k,2) = cldtau(i,k) * taufac(i,k)
          clouds(i,k,3) = cldssa_def
          clouds(i,k,4) = cldasy_def
        enddo
      enddo

!
!===> ... compute low, mid, high, total, and boundary layer cloud fractions
!         and clouds top/bottom layer indices for low, mid, and high clouds.
!         the three cloud domain boundaries are defined by ptopc.  the cloud
!         overlapping method is defined by control flag 'iovr', which may
!         be different for lw and sw radiation programs.

      call gethml                                                       &
!  ---  inputs:
     &     ( plyr, ptop1, cldtot, cldcnv,                               &
     &       ix, nlay,                                                  &
!  ---  outputs:
     &       clds, mtop, mbot                                           &
     &     )

!
      return
!...................................
      end subroutine diagcld1
!-----------------------------------


!-----------------------------------                                    !
      subroutine gethml                                                 &
!...................................                                    !

!  ---  inputs:
     &     ( plyr, ptop1, cldtot, cldcnv,                               &
     &       ix, nlay,                                                  &
!  ---  outputs:
     &       clds, mtop, mbot                                           &
     &     )

!  ===================================================================  !
!                                                                       !
! abstract: compute high, mid, low, total, and boundary cloud fractions !
!   and cloud top/bottom layer indices for model diagnostic output.     !
!   the three cloud domain boundaries are defined by ptopc.  the cloud  !
!   overlapping method is defined by control flag 'iovr', which is also !
!   used by lw and sw radiation programs.                               !
!                                                                       !
! usage:         call gethml                                            !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (ix,nlay) : model layer mean pressure in mb (100pa)           !
!   ptop1 (ix,4)    : pressure limits of cloud domain interfaces        !
!                     (sfc,low,mid,high) in mb (100pa)                  !
!   cldtot(ix,nlay) : total or straiform cloud profile in fraction      !
!   cldcnv(ix,nlay) : convective cloud (for diagnostic scheme only)     !
!   ix              : horizontal dimention                              !
!   nlay            : vertical layer dimensions                         !
!                                                                       !
! output variables:                                                     !
!   clds  (ix,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (ix,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (ix,3)    : vertical indices for low, mid, hi cloud bases     !
!                                                                       !
! external module variables:  (in physpara)                             !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!                                                                       !
! internal module variables:                                            !
!   iovr            : control flag for cloud overlap                    !
!                     =0 random overlapping clouds                      !
!                     =1 max/ran overlapping clouds                     !
!                                                                       !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none!

!  ---  inputs:
      integer, intent(in) :: ix, nlay

      real (kind=kind_phys), dimension(:,:), intent(in) :: plyr, ptop1, &
     &       cldtot, cldcnv

!  ---  outputs
      real (kind=kind_phys), dimension(:,:), intent(out) :: clds

      integer,               dimension(:,:), intent(out) :: mtop, mbot

!  ---  local variables:
      real (kind=kind_phys) :: cl1(ix), cl2(ix)
      real (kind=kind_phys) :: pcur, pnxt, ccur, cnxt

      integer, dimension(ix):: idom, kbt1, kth1, kbt2, kth2
      integer :: i, k, id, id1, kstr, kend, kinc

!
!===> ... begin here
!
      clds(:,:) = 0.0

      do i = 1, ix
        cl1(i) = 1.0
        cl2(i) = 1.0
      enddo

!  ---  total and bl clouds, where cl1, cl2 are fractions of clear-sky view
!       layer processed from surface and up

      if ( ivflip == 0 ) then                   ! input data from toa to sfc
        kstr = nlay
        kend = 1
        kinc = -1
      else                                      ! input data from sfc to toa
        kstr = 1
        kend = nlay
        kinc = 1
      endif                                     ! end_if_ivflip

      if ( iovr == 0 ) then                     ! random overlap

        do k = kstr, kend, kinc
          do i = 1, ix
            ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))
            if (ccur >= climit) cl1(i) = cl1(i) * (1.0 - ccur)
          enddo

          if (k == llyr) then
            do i = 1, ix
              clds(i,5) = 1.0 - cl1(i)          ! save bl cloud
            enddo
          endif
        enddo

        do i = 1, ix
          clds(i,4) = 1.0 - cl1(i)              ! save total cloud
        enddo

      else                                      ! max/ran overlap

        do k = kstr, kend, kinc
          do i = 1, ix
            ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))
            if (ccur >= climit) then             ! cloudy layer
              cl2(i) = min( cl2(i), (1.0 - ccur) )
            else                                ! clear layer
              cl1(i) = cl1(i) * cl2(i)
              cl2(i) = 1.0
            endif
          enddo

          if (k == llyr) then
            do i = 1, ix
              clds(i,5) = 1.0 - cl1(i) * cl2(i) ! save bl cloud
            enddo
          endif
        enddo

        do i = 1, ix
          clds(i,4) = 1.0 - cl1(i) * cl2(i)     ! save total cloud
        enddo

      endif                                     ! end_if_iovr

!  ---  high, mid, low clouds, where cl1, cl2 are cloud fractions
!       layer processed from one layer below llyr and up
!  ---  change! layer processed from surface to top, so low clouds will
!       contains both bl and low clouds.

      if ( ivflip == 0 ) then                   ! input data from toa to sfc

        do i = 1, ix
          cl1 (i) = 0.0
          cl2 (i) = 0.0
          kbt1(i) = nlay
          kbt2(i) = nlay
          kth1(i) = 0
          kth2(i) = 0
          idom(i) = 1
          mbot(i,1) = nlay
          mtop(i,1) = nlay
          mbot(i,2) = nlay - 1
          mtop(i,2) = nlay - 1
          mbot(i,3) = nlay - 1
          mtop(i,3) = nlay - 1
        enddo

!org    do k = llyr-1, 1, -1
        do k = nlay, 1, -1
          do i = 1, ix
            id = idom(i)
            id1= id + 1

            pcur = plyr(i,k)
            ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))

            if (k > 1) then
              pnxt = plyr(i,k-1)
              cnxt = min( ovcst, max( cldtot(i,k-1), cldcnv(i,k-1) ))
            else
              pnxt = -1.0
              cnxt = 0.0
            endif

            if (pcur < ptop1(i,id1)) then
              id = id + 1
              id1= id1 + 1
              idom(i) = id
            endif

            if (ccur >= climit) then
              if (kth2(i) == 0) kbt2(i) = k
              kth2(i) = kth2(i) + 1

              if ( iovr == 0 ) then
                cl2(i) = cl2(i) + ccur - cl2(i)*ccur
              else
                cl2(i) = max( cl2(i), ccur )
              endif

              if (cnxt < climit .or. pnxt < ptop1(i,id1)) then
                kbt1(i) = nint( (cl1(i)*kbt1(i) + cl2(i)*kbt2(i) )      &
     &                  / (cl1(i) + cl2(i)) )
                kth1(i) = nint( (cl1(i)*kth1(i) + cl2(i)*kth2(i) )      &
     &                  / (cl1(i) + cl2(i)) )
                cl1 (i) = cl1(i) + cl2(i) - cl1(i)*cl2(i)

                kbt2(i) = k - 1
                kth2(i) = 0
                cl2 (i) = 0.0
              endif   ! end_if_cnxt_or_pnxt
            endif     ! end_if_ccur

            if (pnxt < ptop1(i,id1)) then
              clds(i,id) = cl1(i)
              mtop(i,id) = min( kbt1(i), kbt1(i)-kth1(i)+1 )
              mbot(i,id) = kbt1(i)

              cl1 (i) = 0.0
              kbt1(i) = k - 1
              kth1(i) = 0

              if (id1 <= nk_clds) then
                mbot(i,id1) = kbt1(i)
                mtop(i,id1) = kbt1(i)
              endif
            endif     ! end_if_pnxt

          enddo       ! end_do_i_loop
        enddo         ! end_do_k_loop

      else                                      ! input data from sfc to toa

        do i = 1, ix
          cl1 (i) = 0.0
          cl2 (i) = 0.0
          kbt1(i) = 1
          kbt2(i) = 1
          kth1(i) = 0
          kth2(i) = 0
          idom(i) = 1
          mbot(i,1) = 1
          mtop(i,1) = 1
          mbot(i,2) = 2
          mtop(i,2) = 2
          mbot(i,3) = 2
          mtop(i,3) = 2
        enddo

!org    do k = llyr+1, nlay
        do k = 1, nlay
          do i = 1, ix
            id = idom(i)
            id1= id + 1

            pcur = plyr(i,k)
            ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))

            if (k < nlay) then
              pnxt = plyr(i,k+1)
              cnxt = min( ovcst, max( cldtot(i,k+1), cldcnv(i,k+1) ))
            else
              pnxt = -1.0
              cnxt = 0.0
            endif

            if (pcur < ptop1(i,id1)) then
              id = id + 1
              id1= id1 + 1
              idom(i) = id
            endif

            if (ccur >= climit) then
              if (kth2(i) == 0) kbt2(i) = k
              kth2(i) = kth2(i) + 1

              if ( iovr == 0 ) then
                cl2(i) = cl2(i) + ccur - cl2(i)*ccur
              else
                cl2(i) = max( cl2(i), ccur )
              endif

              if (cnxt < climit .or. pnxt < ptop1(i,id1)) then
                kbt1(i) = nint( (cl1(i)*kbt1(i) + cl2(i)*kbt2(i))       &
     &                  / (cl1(i) + cl2(i)) )
                kth1(i) = nint( (cl1(i)*kth1(i) + cl2(i)*kth2(i))       &
     &                  / (cl1(i) + cl2(i)) )
                cl1 (i) = cl1(i) + cl2(i) - cl1(i)*cl2(i)

                kbt2(i) = k + 1
                kth2(i) = 0
                cl2 (i) = 0.0
              endif     ! end_if_cnxt_or_pnxt
            endif       ! end_if_ccur

            if (pnxt < ptop1(i,id1)) then
              clds(i,id) = cl1(i)
              mtop(i,id) = max( kbt1(i), kbt1(i)+kth1(i)-1 )
              mbot(i,id) = kbt1(i)

              cl1 (i) = 0.0
              kbt1(i) = k + 1
              kth1(i) = 0

              if (id1 <= nk_clds) then
                mbot(i,id1) = kbt1(i)
                mtop(i,id1) = kbt1(i)
              endif
            endif     ! end_if_pnxt

          enddo       ! end_do_i_loop
        enddo         ! end_do_k_loop

      endif                                     ! end_if_ivflip

!
      return
!...................................
      end subroutine gethml
!-----------------------------------


!-----------------------------------
!> Cloud fraction scheme by G. Thompson (NCAR-RAL), not intended for
!! combining with any cumulus or shallow cumulus parameterization
!! scheme cloud fractions.  This is intended as a stand-alone for
!! cloud fraction and is relatively good at getting widespread stratus
!! and stratoCu without caring whether any deep/shallow Cu param schemes
!! is making sub-grid-spacing clouds/precip.  Under the hood, this
!! scheme follows Mocko and Cotton (1995) in application of the
!! Sundqvist et al (1989) scheme but using a grid-scale dependent
!! RH threshold, one each for land v. ocean points based on
!! experiences with HWRF testing.
      SUBROUTINE cal_cldfra3(CLDFRA, qv, qc, qi, qs, dz,                &
     &                 p, t, XLAND, gridkm,                             &
     &                 modify_qvapor, max_relh,                         &
     &                 kts,kte, debug_flag)
!
      USE module_mp_thompson_new   , ONLY : rsif, rslf
      IMPLICIT NONE
!
      INTEGER, INTENT(IN):: kts, kte
      LOGICAL, INTENT(IN):: modify_qvapor
      REAL, DIMENSION(kts:kte), INTENT(INOUT):: qv, qc, qi, cldfra
      REAL, DIMENSION(kts:kte), INTENT(IN):: p, t, dz, qs
      REAL, INTENT(IN):: gridkm, XLAND, max_relh
      LOGICAL, INTENT(IN):: debug_flag

!..Local vars.
      REAL:: RH_00L, RH_00O, RH_00
      REAL:: entrmnt=0.5
      INTEGER:: k
      REAL:: TC, qvsi, qvsw, RHUM, delz, var_temp
      REAL, DIMENSION(kts:kte):: qvs, rh, rhoa
      integer:: ndebug = 0

      character*512 dbg_msg

!+---+

!..Initialize cloud fraction, compute RH, and rho-air.

         DO k = kts,kte
            CLDFRA(K) = 0.0

            qvsw = rslf(P(k), t(k))
            qvsi = rsif(P(k), t(k))

            tc = t(k) - 273.15
            if (tc .ge. -12.0) then
               qvs(k) = qvsw
            elseif (tc .lt. -35.0) then
               qvs(k) = qvsi
            else
               qvs(k) = qvsw - (qvsw-qvsi)*(-12.0-tc)/(-12.0+35.)
            endif

            if (modify_qvapor) then
               if (qc(k).gt.1.E-8) then
                  qv(k) = MAX(qv(k), qvsw)
                  qvs(k) = qvsw
               endif
               if (qc(k).le.1.E-8 .and. qi(k).ge.1.E-9) then
                  qv(k) = MAX(qv(k), qvsi*1.005)    !..To ensure a tiny bit ice supersaturated
                  qvs(k) = qvsi
               endif
            endif

            rh(k) = MAX(0.01, qv(k)/qvs(k))
            rhoa(k) = p(k)/(287.0*t(k))

         ENDDO


!..First cut scale-aware. Higher resolution should require closer to
!.. saturated grid box for higher cloud fraction.  Simple functions
!.. chosen based on Mocko and Cotton (1995) starting point and desire
!.. to get near 100% RH as grid spacing moves toward 1.0km, but higher
!.. RH over ocean required as compared to over land.

      DO k = kts,kte

         delz = MAX(100., dz(k))
         RH_00L = 0.77+MIN(0.22,SQRT(1./(50.0+gridkm*gridkm*delz*0.01)))
         RH_00O = 0.85+MIN(0.14,SQRT(1./(50.0+gridkm*gridkm*delz*0.01)))
         RHUM = rh(k)

         if (qc(k).ge.1.E-6 .or. qi(k).ge.1.E-7                         &
     &                    .or. (qs(k).gt.1.E-6 .and. t(k).lt.273.)) then
            CLDFRA(K) = 1.0
         elseif (((qc(k)+qi(k)).gt.1.E-10) .and.                        &
     &                                    ((qc(k)+qi(k)).lt.1.E-6)) then
            var_temp = MIN(0.99, 0.1*(11.0 + log10(qc(k)+qi(k))))
            CLDFRA(K) = var_temp*var_temp
         else

            IF ((XLAND-1.5).GT.0.) THEN                                  !--- Ocean
               RH_00 = RH_00O
            ELSE                                                         !--- Land
               RH_00 = RH_00L
            ENDIF

            tc = MAX(-80.0, t(k) - 273.15)
            if (tc .lt. -24.0) RH_00 = RH_00L

            if (tc .gt. 20.0) then
               CLDFRA(K) = 0.0
            elseif (tc .ge. -24.0) then
               RHUM = MIN(rh(k), 1.0)
               var_temp = SQRT(SQRT((1.001-RHUM)/(1.001-RH_00)))
               CLDFRA(K) = MAX(0., 1.0-var_temp)
            else
               if (max_relh.gt.1.12 .or. (.NOT.(modify_qvapor)) ) then
!..For HRRR model, the following look OK.
                  RHUM = MIN(rh(k), 1.45)
                  RH_00 = RH_00 + (1.45-RH_00)*(-24.0-tc)/(-24.0+85.)
                  var_temp = SQRT(SQRT((1.46-RHUM)/(1.46-RH_00)))
                  CLDFRA(K) = MAX(0., 1.0-var_temp)
               else
!..but for the GFS model, RH is way lower.
                  RHUM = MIN(rh(k), 1.05)
                  RH_00 = RH_00 + (1.05-RH_00)*(-24.0-tc)/(-24.0+85.)
                  var_temp = SQRT(SQRT((1.06-RHUM)/(1.06-RH_00)))
                  CLDFRA(K) = MAX(0., 1.0-var_temp)
               endif
            endif
            if (CLDFRA(K).gt.0.) CLDFRA(K)=MAX(0.01,MIN(CLDFRA(K),0.99))

         endif
         if (cldfra(k).gt.0.0 .and. p(k).lt.7000.0) CLDFRA(K) = 0.0
      ENDDO

      call find_cloudLayers(qvs, cldfra, T, P, Dz, entrmnt,             &
     &                      debug_flag, qc, qi, qs, kts,kte)

!..Do a final total column adjustment since we may have added more than 1 mm
!.. LWP/IWP for multiple cloud decks.

      call adjust_cloudFinal(cldfra, qc, qi, rhoa, dz, kts,kte)

!..Intended for cold start model runs, we use modify_qvapor to ensure that cloudy
!.. areas are actually saturated such that the inserted clouds do not evaporate a
!.. timestep later.

      if (modify_qvapor) then
         DO k = kts,kte
            if (cldfra(k).gt.0.2 .and. cldfra(k).lt.1.0) then
               qv(k) = MAX(qv(k),qvs(k))
            endif
         ENDDO
      endif

!...................................                                    !
      END SUBROUTINE cal_cldfra3
!-----------------------------------                                    !


!-----------------------------------                                    !
!>From cloud fraction array, find clouds of multi-level depth and compute
!! a reasonable value of LWP or IWP that might be contained in that depth,
!! unless existing LWC/IWC is already there.
      SUBROUTINE find_cloudLayers(qvs1d, cfr1d, T1d, P1d, Dz1d, entrmnt,&
     &                            debugfl, qc1d, qi1d, qs1d, kts,kte)
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN):: kts, kte
      LOGICAL, INTENT(IN):: debugfl
      REAL, INTENT(IN):: entrmnt
      REAL, DIMENSION(kts:kte), INTENT(IN):: qs1d,qvs1d,T1d,P1d,Dz1d
      REAL, DIMENSION(kts:kte), INTENT(INOUT):: cfr1d, qc1d, qi1d

!..Local vars.
      REAL, DIMENSION(kts:kte):: theta
      REAL:: theta1, theta2, delz
      INTEGER:: k, k2, k_tropo, k_m12C, k_cldb, k_cldt, kbot, k_p200
      LOGICAL:: in_cloud
      character*512 dbg_msg

!+---+

      k_m12C = 0
      k_p200 = 0
      DO k = kte, kts, -1
         theta(k) = T1d(k)*((100000.0/P1d(k))**(287.05/1004.))
         if (T1d(k)-273.16 .gt. -12.0 .and. P1d(k).gt.10100.0)          &
     &                   k_m12C = MAX(k_m12C, k)
         if (P1d(k).gt.19999.0 .and. k_p200.eq.0) k_p200 = k
      ENDDO
      if (k_m12C .le. kts) k_m12C = kts

!..Find tropopause height, best surrogate, because we would not really
!.. wish to put fake clouds into the stratosphere.  The 10/1500 ratio
!.. d(Theta)/d(Z) approximates a vertical line on typical SkewT chart
!.. near typical (mid-latitude) tropopause height.  Since messy data
!.. could give us a false signal of such a transition, do the check over 
!.. three K-level change, not just a level-to-level check.  This method
!.. has potential failure in arctic-like conditions with extremely low
!.. tropopause height, as would any other diagnostic, so ensure resulting
!.. k_tropo level is above 700hPa.

      if ( (kte-k_p200) .lt. 3) k_p200 = kte-3
      DO k = k_p200-2, kts, -1
         theta1 = theta(k)
         theta2 = theta(k+2)
         delz = 0.5*dz1d(k) + dz1d(k+1) + 0.5*dz1d(k+2)
         if ( (((theta2-theta1)/delz).lt.10./1500.) .OR.                &
     &                       P1d(k).gt.70000.) EXIT
      ENDDO
      k_tropo = MAX(kts+2, MIN(k+2, kte-1))

      if (k_tropo .gt. k_p200) then
         DO k = kte-3, k_p200-2, -1
            theta1 = theta(k)
            theta2 = theta(k+2)
            delz = 0.5*dz1d(k) + dz1d(k+1) + 0.5*dz1d(k+2)
            if ( (((theta2-theta1)/delz).lt.10./1500.) .AND.            &
     &                                    P1d(k).gt.9000.) EXIT
         ENDDO
         k_tropo = MAX(k_p200-1, MIN(k+2, kte-1))
      endif

!..Eliminate possible fractional clouds above supposed tropopause.
      DO k = k_tropo+1, kte
         if (cfr1d(k).gt.0.0 .and. cfr1d(k).lt.1.0) then
            cfr1d(k) = 0.
         endif
      ENDDO

!..Be a bit more conservative with lower cloud fraction in scenario with
!.. well-mixed convective boundary layer below LCL.

      kbot = kts+1
      DO k = kbot, k_m12C
         if ( (theta(k)-theta(k-1)) .gt. 0.010E-3*Dz1d(k)) EXIT
      ENDDO
      kbot = MAX(kts+1, k-2)
      DO k = kts, kbot
         if (cfr1d(k).gt.0.0 .and. cfr1d(k).lt.1.0)                     &
     &              cfr1d(k) = MAX(0.01,0.33*cfr1d(k))
      ENDDO
      DO k = kts,k_tropo
         if (cfr1d(k).gt.0.0) kbot = MIN(k,kbot)
      ENDDO

!..Starting below tropo height, if cloud fraction greater than 1 percent,
!.. compute an approximate total layer depth of cloud, determine a total 
!.. liquid water/ice path (LWP/IWP), then reduce that amount with tuning 
!.. parameter to represent entrainment factor, then divide up LWP/IWP
!.. into delta-Z weighted amounts for individual levels per cloud layer. 

      k_cldb = k_tropo
      in_cloud = .false.
      k = k_tropo
      DO WHILE (.not. in_cloud .AND. k.gt.k_m12C+1)
         k_cldt = 0
         if (cfr1d(k).ge.0.01) then
            in_cloud = .true.
            k_cldt = MAX(k_cldt, k)
         endif
         if (in_cloud) then
            DO k2 = k_cldt-1, k_m12C, -1
               if (cfr1d(k2).lt.0.01 .or. k2.eq.k_m12C) then
                  k_cldb = k2+1
                  goto 87
               endif
            ENDDO
 87         continue
            in_cloud = .false.
         endif
         if ((k_cldt - k_cldb + 1) .ge. 2) then
            call adjust_cloudIce(cfr1d, qi1d, qs1d, qvs1d, T1d, Dz1d,   &
     &                           entrmnt, k_cldb,k_cldt,kts,kte)
            k = k_cldb
         elseif ((k_cldt - k_cldb + 1) .eq. 1) then
            if (cfr1d(k_cldb).gt.0.and.cfr1d(k_cldb).lt.1.)             &
     &        qi1d(k_cldb)=qi1d(k_cldb)+0.05*qvs1d(k_cldb)*cfr1d(k_cldb)
            k = k_cldb
         endif
         k = k - 1
      ENDDO

      k_cldb = k_m12C + 3
      in_cloud = .false.
      k = k_m12C + 2
      DO WHILE (.not. in_cloud .AND. k.gt.kbot)
         k_cldt = 0
         if (cfr1d(k).ge.0.01) then
            in_cloud = .true.
            k_cldt = MAX(k_cldt, k)
         endif
         if (in_cloud) then
            DO k2 = k_cldt-1, kbot, -1
               if (cfr1d(k2).lt.0.01 .or. k2.eq.kbot) then
                  k_cldb = k2+1
                  goto 88
               endif
            ENDDO
 88         continue
            in_cloud = .false.
         endif
         if ((k_cldt - k_cldb + 1) .ge. 2) then
            call adjust_cloudH2O(cfr1d, qc1d, qvs1d, T1d, Dz1d,         &
     &                           entrmnt, k_cldb,k_cldt,kts,kte)
            k = k_cldb
         elseif ((k_cldt - k_cldb + 1) .eq. 1) then
            if (cfr1d(k_cldb).gt.0.and.cfr1d(k_cldb).lt.1.)             &
     &        qc1d(k_cldb)=qc1d(k_cldb)+0.05*qvs1d(k_cldb)*cfr1d(k_cldb)
            k = k_cldb
         endif
         k = k - 1
      ENDDO


      END SUBROUTINE find_cloudLayers


!+---+-----------------------------------------------------------------+

!>
      SUBROUTINE adjust_cloudIce(cfr,qi,qs,qvs,T,dz,entr, k1,k2,kts,kte)
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN):: k1,k2, kts,kte
      REAL, INTENT(IN):: entr
      REAL, DIMENSION(kts:kte), INTENT(IN):: cfr, qs, qvs, T, dz
      REAL, DIMENSION(kts:kte), INTENT(INOUT):: qi
      REAL:: iwc, max_iwc, tdz, this_iwc, this_dz
      INTEGER:: k

      tdz = 0.
      do k = k1, k2
         tdz = tdz + dz(k)
      enddo
!     max_iwc = ABS(qvs(k2)-qvs(k1))
      max_iwc = MAX(0.0, qvs(k1)-qvs(k2))

      do k = k1, k2
         max_iwc = MAX(1.E-6, max_iwc - (qi(k)+qs(k)))
      enddo
      max_iwc = MIN(1.E-4, max_iwc)

      this_dz = 0.0
      do k = k1, k2
         if (k.eq.k1) then
            this_dz = this_dz + 0.5*dz(k)
         else
            this_dz = this_dz + dz(k)
         endif
         this_iwc = max_iwc*this_dz/tdz
         iwc = MAX(1.E-6, this_iwc*(1.-entr))
         if (cfr(k).gt.0.0.and.cfr(k).lt.1.0.and.T(k).ge.203.16) then
            qi(k) = qi(k) + cfr(k)*cfr(k)*iwc
         endif
      enddo

      END SUBROUTINE adjust_cloudIce

!+---+-----------------------------------------------------------------+

!>
      SUBROUTINE adjust_cloudH2O(cfr, qc, qvs,T,dz,entr, k1,k2,kts,kte)
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN):: k1,k2, kts,kte
      REAL, INTENT(IN):: entr
      REAL, DIMENSION(kts:kte), INTENT(IN):: cfr, qvs, T, dz
      REAL, DIMENSION(kts:kte), INTENT(INOUT):: qc
      REAL:: lwc, max_lwc, tdz, this_lwc, this_dz
      INTEGER:: k

      tdz = 0.
      do k = k1, k2
         tdz = tdz + dz(k)
      enddo
!     max_lwc = ABS(qvs(k2)-qvs(k1))
      max_lwc = MAX(0.0, qvs(k1)-qvs(k2))
!     print*, ' max_lwc = ', max_lwc, ' over DZ=',tdz

      do k = k1, k2
         max_lwc = MAX(1.E-6, max_lwc - qc(k))
      enddo
      max_lwc = MIN(1.E-4, max_lwc)

      this_dz = 0.0
      do k = k1, k2
         if (k.eq.k1) then
            this_dz = this_dz + 0.5*dz(k)
         else
            this_dz = this_dz + dz(k)
         endif
         this_lwc = max_lwc*this_dz/tdz
         lwc = MAX(1.E-6, this_lwc*(1.-entr))
         if (cfr(k).gt.0.0.and.cfr(k).lt.1.0.and.T(k).ge.258.16) then
            qc(k) = qc(k) + cfr(k)*cfr(k)*lwc
         endif
      enddo

      END SUBROUTINE adjust_cloudH2O

!+---+-----------------------------------------------------------------+

!> Do not alter any grid-explicitly resolved hydrometeors, rather only
!! the supposed amounts due to the cloud fraction scheme.
      SUBROUTINE adjust_cloudFinal(cfr, qc, qi, Rho,dz, kts,kte)
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN):: kts,kte
      REAL, DIMENSION(kts:kte), INTENT(IN):: cfr, Rho, dz
      REAL, DIMENSION(kts:kte), INTENT(INOUT):: qc, qi
      REAL:: lwp, iwp, xfac
      INTEGER:: k

      lwp = 0.
      iwp = 0.
      do k = kts, kte
         if (cfr(k).gt.0.0 .and. cfr(k).lt.1.0) then
            lwp = lwp + qc(k)*Rho(k)*dz(k)
            iwp = iwp + qi(k)*Rho(k)*dz(k)
         endif
      enddo

      if (lwp .gt. 1.0) then
         xfac = 1.0/lwp
         do k = kts, kte
            if (cfr(k).gt.0.0 .and. cfr(k).lt.1.0) then
               qc(k) = qc(k)*xfac
            endif
         enddo
      endif

      if (iwp .gt. 1.0) then
         xfac = 1.0/iwp
         do k = kts, kte
            if (cfr(k).gt.0.0 .and. cfr(k).lt.1.0) then
               qi(k) = qi(k)*xfac
            endif
         enddo
      endif

      END SUBROUTINE adjust_cloudFinal

!+---+-----------------------------------------------------------------+

!> This subroutine computes the Xu-Randall cloud fraction scheme.
      subroutine cloud_fraction_XuRandall                               &
!  ---  inputs:
     &     ( IX, NLAY, plyr, clwf, rhly, qstl, lmfshal, lmfdeep2,       &
!  ---  outputs:
     &       cldtot ) 

!  ---  inputs:
      integer, intent(in) :: IX, NLAY
      real (kind=kind_phys), dimension(:,:), intent(in) :: plyr, clwf,  &
     &                                                     rhly, qstl
      logical, intent(in) :: lmfshal, lmfdeep2

!  ---  outputs
      real (kind=kind_phys), dimension(:,:), intent(inout) :: cldtot

!  ---  local variables:

       real (kind=kind_phys) :: clwmin, clwm, clwt, onemrh, value,      &
     &       tem1, tem2
       real (kind=kind_phys), parameter :: xrc3 = 100.
       integer :: i, k

!> - Compute layer cloud fraction.

       if ( ivflip == 0 ) then              ! input data from toa to sfc

         clwmin = 0.0
         if (.not. lmfshal) then
           do k = NLAY, 1, -1
           do i = 1, IX
             clwt = 1.0e-6 * (plyr(i,k)*0.001)

             if (clwf(i,k) > clwt) then

               onemrh= max( 1.e-10, 1.0-rhly(i,k) )
               clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )

               tem1  = min(max(sqrt(sqrt(onemrh*qstl(i,k))),0.0001),1.0)
               tem1  = 2000.0 / tem1

               value = max( min( tem1*(clwf(i,k)-clwm), 50.0 ), 0.0 )
               tem2  = sqrt( sqrt(rhly(i,k)) )

               cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
             endif
          enddo
          enddo
        else
          do k = NLAY, 1, -1
          do i = 1, IX
            clwt = 1.0e-6 * (plyr(i,k)*0.001)

            if (clwf(i,k) > clwt) then
              onemrh= max( 1.e-10, 1.0-rhly(i,k) )
              clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )

              tem1  = min(max((onemrh*qstl(i,k))**0.49,0.0001),1.0)  !jhan
              if (lmfdeep2) then
                tem1  = xrc3 / tem1
              else
                tem1  = 100.0 / tem1
              endif

              value = max( min( tem1*(clwf(i,k)-clwm), 50.0 ), 0.0 )
              tem2  = sqrt( sqrt(rhly(i,k)) )
              cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
            endif
          enddo
          enddo
        endif

      else                                 ! input data from sfc to toa

        clwmin = 0.0
        if (.not. lmfshal) then
          do k = 1, NLAY
          do i = 1, IX
            clwt = 1.0e-6 * (plyr(i,k)*0.001)

            if (clwf(i,k) > clwt) then

              onemrh= max( 1.e-10, 1.0-rhly(i,k) )
              clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )

              tem1  = min(max(sqrt(sqrt(onemrh*qstl(i,k))),0.0001),1.0)
              tem1  = 2000.0 / tem1

              value = max( min( tem1*(clwf(i,k)-clwm), 50.0 ), 0.0 )
              tem2  = sqrt( sqrt(rhly(i,k)) )

              cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
            endif
          enddo
          enddo
        else
          do k = 1, NLAY
          do i = 1, IX
            clwt = 1.0e-6 * (plyr(i,k)*0.001)

            if (clwf(i,k) > clwt) then
              onemrh= max( 1.e-10, 1.0-rhly(i,k) )
              clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )
!
              tem1  = min(max((onemrh*qstl(i,k))**0.49,0.0001),1.0)  !jhan
              if (lmfdeep2) then
                tem1  = xrc3 / tem1
              else
                tem1  = 100.0 / tem1
              endif

              value = max( min( tem1*(clwf(i,k)-clwm), 50.0 ), 0.0 )
              tem2  = sqrt( sqrt(rhly(i,k)) )

              cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
            endif
          enddo
          enddo
        endif
      endif

      end subroutine cloud_fraction_XuRandall

      subroutine cloud_fraction1d_XuRandall                             &
!  ---  inputs:
     &     ( NLAY, plyr, clwf, rhly, qstl, lmfshal, lmfdeep2,           &
!  ---  outputs:
     &       cldtot ) 

!  ---  inputs:
      integer, intent(in) :: NLAY
      real (kind=kind_phys), dimension(:), intent(in)   :: plyr, clwf,  &
     &                                                     rhly, qstl
      logical, intent(in) :: lmfshal, lmfdeep2

!  ---  outputs
      real (kind=kind_phys), dimension(:), intent(inout) :: cldtot

!  ---  local variables:

       real (kind=kind_phys) :: clwmin, clwm, clwt, onemrh, value,      &
     &       tem1, tem2
       real (kind=kind_phys), parameter :: xrc3 = 100.
       integer :: k

!> - Compute layer cloud fraction.

      clwmin = 0.0
      if (.not. lmfshal) then
        do k = 1, NLAY
          clwt = 1.0e-6 * (plyr(k)*0.001)

          if (clwf(k) > clwt) then

            onemrh= max( 1.e-10, 1.0-rhly(k) )
            clwm  = clwmin / max( 0.01, plyr(k)*0.001 )

            tem1  = min(max(sqrt(sqrt(onemrh*qstl(k))),0.0001),1.0)
            tem1  = 2000.0 / tem1

            value = max( min( tem1*(clwf(k)-clwm), 50.0 ), 0.0 )
            tem2  = sqrt( sqrt(rhly(k)) )

            cldtot(k) = max( tem2*(1.0-exp(-value)), 0.0 )
          endif
        enddo
      else
        do k = 1, NLAY
          clwt = 1.0e-6 * (plyr(k)*0.001)

          if (clwf(k) > clwt) then
            onemrh= max( 1.e-10, 1.0-rhly(k) )
            clwm  = clwmin / max( 0.01, plyr(k)*0.001 )

            tem1  = min(max((onemrh*qstl(k))**0.49,0.0001),1.0)  !jhan
            if (lmfdeep2) then
              tem1  = xrc3 / tem1
            else
              tem1  = 100.0 / tem1
            endif

            value = max( min( tem1*(clwf(k)-clwm), 50.0 ), 0.0 )
            tem2  = sqrt( sqrt(rhly(k)) )

            cldtot(k) = max( tem2*(1.0-exp(-value)), 0.0 )
          endif
        enddo
      endif

      end subroutine cloud_fraction1d_XuRandall
!+---+-----------------------------------------------------------------+


!-----------------------------------                                    !
      subroutine rhtable                                                &
!...................................                                    !

!  ---  inputs:
     &     ( me                                                         &
!  ---  outputs:
     &,      ier )

!  ===================================================================  !
!                                                                       !
! abstract: cld-rh relations obtained from mitchell-hahn procedure,     !
!   here read cld/rh tuning tables for day 0,1,...,5 and merge into 1   !
!   file.                                                               !
!                                                                       !
! inputs:                                                               !
!   me              : check print control flag                          !
!                                                                       !
! outputs:                                                              !
!   ier             : error flag                                        !
!                                                                       !
!  ===================================================================  !
!
      implicit none!

!  ---  inputs:
      integer, intent(in) :: me

!  ---  output:
      integer, intent(out) :: ier

!  ---  locals:
      real (kind=kind_phys), dimension(nbin,nlon,nlat,mcld,nseal) ::    &
     &      rhfd, rtnfd, rhcf, rtncf, rhcla

      real (kind=kind_io4), dimension(nbin,nlon,nlat,mcld,nseal) ::     &
     &      rhfd4, rtnfd4

      real(kind=kind_io4)  :: fhour

      real(kind=kind_phys) :: binscl, cfrac, clsat, rhsat, cstem

      integer, dimension(nlon,nlat,mcld,nseal) :: kpts, kkpts

      integer :: icdays(15), idate(4), nbdayi, isat

      integer :: i, i1, j, k, l, m, id, im, iy

!
!===> ...  begin here
!

      ier = 1

      rewind nicltun

      binscl = 1.0 / nbin

!  ---  array initializations

      do m=1,nseal
       do l=1,mcld
        do k=1,nlat
         do j=1,nlon
          do i=1,nbin
            rhcf (i,j,k,l,m) = 0.0
            rtncf(i,j,k,l,m) = 0.0
            rhcla(i,j,k,l,m) = -0.1
          enddo
         enddo
        enddo
       enddo
      enddo

      kkpts = 0

!  ---  read the data off the rotating file

      read (nicltun,err=998,end=999) nbdayi, icdays

      if (me == 0) print 11, nbdayi
  11  format('   from rhtable days on file =',i5)

      do i = 1, nbdayi
       id = icdays(i) / 10000
       im = (icdays(i)-id*10000) / 100
       iy = icdays(i)-id*10000-im*100
       if (me == 0) print 51, id,im,iy
  51   format('   from rhtable archv data from da,mo,yr=',3i4)
      enddo

      read (nicltun,err=998,end=999) fhour,idate

      do i1 = 1, nbdayi
        read (nicltun) rhfd4
        rhfd = rhfd4

        read (nicltun) rtnfd4
        rtnfd = rtnfd4

        read (nicltun) kpts

        do m=1,nseal
         do l=1,mcld
          do k=1,nlat
           do j=1,nlon
            do i=1,nbin
              rhcf (i,j,k,l,m) = rhcf (i,j,k,l,m) + rhfd (i,j,k,l,m)
              rtncf(i,j,k,l,m) = rtncf(i,j,k,l,m) + rtnfd(i,j,k,l,m)
            enddo
           enddo
          enddo
         enddo
        enddo

        kkpts = kkpts + kpts

      enddo     ! end_do_i1_loop

      do m = 1, nseal
       do l = 1, mcld
        do k = 1, nlat
         do j = 1, nlon

!  ---  compute the cumulative frequency distribution

           do i = 2, nbin
             rhcf (i,j,k,l,m) = rhcf (i-1,j,k,l,m) + rhcf (i,j,k,l,m)
             rtncf(i,j,k,l,m) = rtncf(i-1,j,k,l,m) + rtncf(i,j,k,l,m)
           enddo   ! end_do_i_loop

           if (kkpts(j,k,l,m) > 0) then
             do i = 1, nbin
               rhcf (i,j,k,l,m)= rhcf (i,j,k,l,m)/kkpts(j,k,l,m)
               rtncf(i,j,k,l,m)=min(1., rtncf(i,j,k,l,m)/kkpts(j,k,l,m))
             enddo

!  ---  cause we mix calculations of rh retune with cray and ibm words
!       the last value of rhcf is close to but ne 1.0,
!  ---  so we reset it in order that the 360 loop gives complete tabl

             rhcf(nbin,j,k,l,m) = 1.0

             do i = 1, nbin
               lab_do_i1 : do i1 = 1, nbin
                 if (rhcf(i1,j,k,l,m) >= rtncf(i,j,k,l,m)) then
                   rhcla(i,j,k,l,m) = i1 * binscl
                   exit  lab_do_i1
                 endif
               enddo  lab_do_i1
             enddo

           else                   ! if_kkpts
!  ---  no critical rh

             do i = 1, nbin
               rhcf (i,j,k,l,m) = -0.1
               rtncf(i,j,k,l,m) = -0.1
             enddo

             if (me == 0) then
               print 210, k,j,m
 210           format('  no crit rh for lat=',i3,' and lon band=',i3,   &
     &                ' land(=1) sea=',i3//'  model rh ',' obs rtcld')
               do i = 1, nbin
                 print 203, rhcf(i,j,k,l,m), rtncf(i,j,k,l,m)
 203             format(2f10.2)
               enddo
             endif

           endif               ! if_kkpts

         enddo    ! end_do_j_loop
        enddo     ! end_do_k_loop
       enddo      ! end_do_l_loop
      enddo       ! end_do_m_loop

      do m = 1, nseal
       do l = 1, mcld
        do k = 1, nlat
         do j = 1, nlon

           isat = 0
           do i = 1, nbin-1
             cfrac = binscl * (i - 1)

             if (rhcla(i,j,k,l,m) < 0.0) then
               print 1941, i,m,l,k,j
 1941          format('  neg rhcla for it,nsl,nc,lat,lon=',5i4          &
     &,               '...stoppp..')
               stop
             endif

             if (rtncf(i,j,k,l,m) >= 1.0) then
               if (isat <= 0) then
                 isat  = i
                 rhsat = rhcla(i,j,k,l,m)
                 clsat = cfrac
               endif

               rhcla(i,j,k,l,m) = rhsat + (1.0 - rhsat)                 &
     &                         * (cfrac - clsat) / (1.0 - clsat)
             endif
           enddo

           rhcla(nbin,j,k,l,m) = 1.0

         enddo    ! end_do_j_loop
        enddo     ! end_do_k_loop
       enddo      ! end_do_l_loop
      enddo       ! end_do_m_loop

!  ---  smooth out the table as it reaches rh=1.0, via linear interpolation
!       between location of rh ge .98 and the nbin bin (where rh=1.0)
!       previously rh=1.0 occurred for many of the latter bins in the
!  ---  table, thereby giving a cloud value of less then 1.0 for rh=1.0

      rhcl = rhcla

      do m = 1, nseal
       do l = 1, mcld
        do k = 1, nlat
         do j = 1, nlon

           lab_do_i : do i = 1, nbin - 2
             cfrac = binscl * i

             if (rhcla(i,j,k,l,m) >= 0.98) then
               do i1 = i, nbin
                 cstem = binscl * i1

                 rhcl(i1,j,k,l,m) = rhcla(i,j,k,l,m)                    &
     &                    + (rhcla(nbin,j,k,l,m) - rhcla(i,j,k,l,m))    &
     &                    * (cstem - cfrac) / (1.0 - cfrac)
               enddo
               exit  lab_do_i
             endif
           enddo  lab_do_i

         enddo    ! end_do_j_loop
        enddo     ! end_do_k_loop
       enddo      ! end_do_l_loop
      enddo       ! end_do_m_loop

      if (me == 0) then
        print *,'completed rhtable for cloud tuninig tables'
      endif
      return

 998  print 988
 988  format(' from rhtable error reading tables')
      ier = -1
      return

 999  print 989
 989  format(' from rhtable e.o.f reading tables')
      ier = -1
      return

!...................................
      end subroutine rhtable
!-----------------------------------


!
!........................................!
      end module module_radiation_clouds !
!========================================!
