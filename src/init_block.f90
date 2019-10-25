      subroutine init_block
!
! modify to f90 by C-H Lee and sort by River Chen in 2015
!
      use param
      use const
      use noah
      use radn
      use physcons, only :con_cp  ,con_rerth,con_omega,con_g   , &
                          con_sbc ,con_solr ,con_hvap ,con_hfus, &
                          con_tice

      implicit none

      cp   = con_cp
      rad  = con_rerth
      omega= con_omega
      grav = con_g
      stbo= con_sbc
      s0  = con_solr
      hltm= con_hvap
      tice= con_tice
      hice= con_hfus!
!!      cp=1004.24
!!      rad=6.371e6
!!      omega=7.292e-5
!!      grav=9.80616

      dt=900.0
      tfilt=0.02

      frad=1.0
      ptop=1.0
      ptmean=1000.0
      ksgeo=3
      ktcup=16

      ktpbl=16
      ktshl=16
      njump=2
      ldiag=1
      idg=16
      jdg=16

!!      stbo=5.669e-8
!!      s0=1368.3
!!      hltm=2.52e6
!!      tice=273.15

!!      hice=3.336e5
      evaprh=0.98
      hfilt=1.5e15

      nnmiit=3
      nnmivm=3
      cutfreq=1.0
      itypbl=0

      taup=6.
      ktrop=8
      taureg=6.
!
      lsimpl=.true.
      lzadv=.true.
      yesdia=.true.

      dopbl=.true.
      docup=.true.
      dorad=.true.

      dolsp=.true.
      dograv=.true.
      docgrav=.true.
      doshl=.true.

      donnmi=.true.
      dodry=.false.
      ozon=.true.

      hdiff=.true.
      cstar=.false.

      update=.true.
      doincr=.true.
!
      doo3l=.true.
!
      domfc=384.
      otgreen=6.
      out_green=.false.
      out_hp=.false.
! pdf cloud
      pdfcloud=.false.
! sppt
      dosppt=.false.
      dospptout=.false.
!---------------------------------------------------------------------------
! sppt parameters, de_corretime unit is hour
!
      de_corretime_500=6.
      de_corretime_1000=72.    ! 3d=3.*24.
      de_corretime_2000=720.   ! 30d=30.*24.
!
! perturbation amplitude factor, basically do not great than 0.17
!
      facsppt500=0.15
      facsppt1000=0.05
      facsppt2000=0.015
!
! specify the default option for cup and pbl
!
      nmcup=3
      nmpbl=2
      nmland=2
      nmshl=1

      cgw=1.0e-4
!
!---for using forecast daily sst, sea ice fraction, snow depth
      ldailyFCTsst=.false.
      ldailyFCTicesndpt=.false.
      lFCTweight=.false.
      lopgsst=.false.
      dailyClm_option=-99
!--for SIT
      do_sit=.false.
      fsit=-99.
!
! specify the default option for orographic and convective gwd
!
      nmgwor=2
      nmgwcv=2
      mtnvar=14
!
! specify the default option for reduced grids
! numreduce : -99 for full grids, 1 to 4 proper for reduced grids
!
      numreduce=-99
      if ( lev .eq. 60 ) then
!
! L60 hybrid coordinate(op7)
         tmean=(/272.1, 265.1, 253.9, 245.4, 238.7, 234.5, 230.7, 227.6, &
                 226.2, 225.0, 223.9, 222.8, 221.7, 220.7, 219.8, 218.8, &
                 217.9, 217.0, 216.7, 216.7, 216.7, 216.7, 216.7, 216.7, &
                 216.7, 216.7, 216.7, 216.7, 216.7, 217.0, 220.6, 224.9, &
                 229.1, 233.2, 237.2, 241.0, 244.8, 248.4, 251.8, 255.1, &
                 258.2, 261.1, 263.9, 266.4, 268.8, 271.0, 273.0, 274.9, &
                 276.6, 278.1, 279.5, 280.8, 281.9, 283.0, 283.9, 284.7, &
                 285.4, 286.1, 286.7, 287.2/)
!
! L60: hybrid coordinate(op7)
          aki=(/  .00000,   .85400,  1.82570,  2.93106,  4.18811,  5.61725 &
            ,    7.24145,  9.08662, 11.18186, 13.55986, 16.25721, 19.31478 &
            ,   22.77809, 26.69766, 31.12936, 36.13468, 41.78099, 48.14158 &
            ,   55.29569, 63.32831, 72.32859, 82.26131, 92.84074,103.75080 &
            ,  114.67568,125.30212,135.32210,144.43608,152.35664,158.81248 &
            ,  163.55255,166.35379,167.12546,165.92496,162.84922,158.02912 &
            ,  151.62865,143.84211,134.88958,125.01073,114.45724,103.48445 &
            ,   92.34269, 81.26901, 70.47979, 60.16483, 50.48301, 41.55990 &
            ,   33.48706, 26.32304, 20.09558, 14.80485, 10.42730,  6.91977 &
            ,    4.22363,  2.26873,   .97701,   .26557,   .02193,   .00000 &
            ,    0.00000/)
!
          bki=(/.00000000,.00000000,.00000000,.00000000,.00000000,.00000000 &
            ,  .00000000,.00000000,.00000000,.00000000,.00000000,.00000000  &
            ,  .00000000,.00000000,.00000000,.00000000,.00000000,.00000000  &
            ,  .00000000,.00000000,.00000104,.00013304,.00077978,.00235711  &
            ,  .00528030,.00995912,.01679191,.02615862,.03841285,.05387320  &
            ,  .07281429,.09545413,.12184549,.15185637,.18527624,.22182071  &
            ,  .26113359,.30279276,.34632012,.39119502,.43687069,.48279246  &
            ,  .52841650,.57322786,.61675643,.65858995,.69838324,.73586353  &
            ,  .77083188,.80316112,.83279100,.85972113,.88400268,.90572942  &
            ,  .92502870,.94205290,.95697164,.96996497,.98124505,.99058760  &
            ,  1.00000000/)
      else if ( lev .eq. 72 ) then
!
! L72 hybrid coordinate
   tmean=(/ 272.07001,270.56509,260.00107,252.01338,245.39818,239.86288, &
            235.97810,232.73827,229.76340,227.40314,226.29744,225.25049, &
            224.25146,223.29295,222.36887,221.47464,220.60663,219.76190, &
            218.93806,218.13326,217.34613,216.64999,216.64999,216.64999, &
            216.64999,216.64999,216.64999,216.64999,216.64999,216.64999, &
            216.64999,216.64999,216.64999,216.64999,216.95261,220.14218, &
            223.91943,227.64014,231.29068,234.86462,238.34796,241.73264, &
            245.00887,248.16821,251.20274,254.10707,256.87320,259.49893, &
            261.98141,264.31860,266.51178,268.55948,270.46774,272.23737, &
            273.87387,275.38428,276.76971,278.04099,279.20291,280.26285, &
            281.22516,282.10135,282.89340,283.61090,284.25690,284.84167, &
            285.36823,285.83975,286.26508,286.64282,286.98969,287.29077 /)
!
! L72: hybrid coordinate
    aki=(/   0.00000,  0.63614,  1.35033,  2.15204,  3.05182,  4.06147,  &
             5.19414,  6.46451,  7.88892,  9.48552, 11.27451, 13.27825,  &
            15.52153, 18.03171, 20.83899, 23.97656, 27.48082, 31.39157,  &
            35.75213, 40.60953, 46.01449, 52.02149, 58.68863, 66.07746,  &
            74.24635, 83.11719, 92.47566,102.10102,111.77279,121.27200,  &
           130.38263,138.89340,146.59965,153.30553,158.82621,162.99021,  &
           165.64246,166.70406,166.20951,164.22232,160.82692,156.12838,  &
           150.25131,143.33774,135.54422,127.03812,117.99331,108.58553,  &
            98.98758, 89.36473, 79.87055, 70.64337, 61.80350, 53.45145,  &
            45.66689, 38.50863, 32.01524, 26.20638, 21.08455, 16.63717,  &
            12.83884,  9.65364,  7.03738,  4.93965,  3.30576,  2.07840,  &
             1.19905,  0.60919,  0.25125,  0.06934,  0.00691,  0.00000,  &
             0.00000 /)

    bki=(/ .00000000,.00000000,.00000000,.00000000,.00000000,.00000000,  &
           .00000000,.00000000,.00000000,.00000000,.00000000,.00000000,  &
           .00000000,.00000000,.00000000,.00000000,.00000000,.00000000,  &
           .00000000,.00000000,.00000000,.00000000,.00000000,.00000000,  &
           .00000624,.00016402,.00075667,.00207486,.00440868,.00804493,  &
           .01326388,.02033548,.02951509,.04103882,.05511859,.07193703,  &
           .09164181,.11428347,.13977611,.16798052,.19871153,.23173872,  &
           .26678885,.30355047,.34168025,.38081121,.42056203,.46054736,  &
           .50038827,.53972232,.57821275,.61555624,.65148892,.68579044,  &
           .71828604,.74884679,.77738809,.80386683,.82827754,.85064780,  &
           .87103331,.88951283,.90618329,.92115521,.93454851,.94648889,  &
           .95710466,.96652416,.97487364,.98227564,.98885067,.99472349,  &
           1.0000000 /)
      endif
!
      tmeans=300.

!-- for hybrid coordinates, ptmeans reset for numerical stability
!      ptmeans=800.
      ptmeans=600.
!
      ifilin ='ifilin'
      ifilout='ifilout'
      cwbout ='cwbout'
      bckfile='bckfile'
      phyout ='phyout'
      namlsts='namlsts'
      crdate ='crdate'
      ocards ='ocards'
      cntrl  ='gfsctl'
!
!dms
!t512l60
      ggdef='gh0g'
      gmdef='ghmg'
      gsdef='gh0s'
!soil
      tsat=(/.421,.464,.468,.434,.406,.465,.404,.439,.421/)
      ref=(/.283,.387,.412,.312,.338,.382,.315,.329,.283/)
      wlt=(/.029,.119,.139,.047,.010,.103,.069,.066,.029/)
!---------------------------------------------------------------------------
! Default values for some radiation controls
!---------------------------------------------------------------------------
! radiation NCEP-RRTMG scheme
!     levr     = 60 ! vertical layers for radiation scheme
!     ictm     = 1  ! ictm=0 => use data at initial cond time, if not
!                   !           available, use latest, no extrapolation.
!                   ! ictm=1 => use data at the forecast time, if
!                   !           not available, use latest and extrapolation.
!                   ! ictm=yyyy0 => use yyyy data for the forecast time,
!                   !               no further data extrapolation.
!                   ! ictm=yyyy1 => use yyyy data for the fcst. if needed,
!                   !               do extrapolation to match the fcst time.
!                   ! ictm=-1 => use user provided external data.
!                   !            for the fcst time, no extrapolation.
!                   ! ictm=-2 => same as ictm=0, but add seasonal cycle
!                   !            from climatology. no extrapolation.
!                   !
!     isol     = 1  ! use noaa old yearly solar constant table with 11-year
!                   ! cycle (range : 1944-2006)
!     ico2     = 2  ! use observed co2 monthly 2-D data table in 15 degree
!                   ! horizontal resolution
!     iaer     =111 ! opac-climatology aerosol scheme, include recorded
!                   ! stratospheric volcanic aerosol effect and
!                   ! toposperic sw/lw areosol effects.
!     ialb     = 0  ! ialb=0, surface vegetation type based climatology scheme,
!                   !         monthly data in 1 degree horizontal resolution.
!                   ! ialb=1, use modis based alb
!     iems     = 1  ! surface type based on climatology in 1 degree horizontal
!                   ! resolution.
!     ntcw     = 1  ! ntcw=0, no cloud condensate calculated
!                   ! ntcw>0, include microphysics cloud scheme
!     num_p3d  = 4  ! num_p3d=4, Zhao Microphysics cloud scheme (default)
!              = 3  ! num_p3d=3, Brad Ferrier's Microphysics cloud scheme
!              = 5  ! num_p3d=5, WSM6 Microphysics cloud scheme
!     ntoz     = 0  ! use climatological ozone profile
!              > 0  ! use interactive ozone profile
!     iovr_sw  = 1  ! sw: maximum-random overlapping vertical cloud layer
!     iovr_lw  = 1  ! lw: maximum-random overlapping vertical cloud layer
!     isubc_sw = 2  ! sw with mcica sub-col approximation provided random seed
!     isubc_lw = 2  ! lw with mcica sub-col approximation provided random seed
!                   ! =1 => sub-grid cloud with prescribed seeds
!                   ! =2 => sub-grid cloud with randomly generated
!     icice_sw = 3  ! sw cloud optical property for cloud ice
!     icice_lw = 3  ! lw cloud optical property for cloud ice
!                   ! =1 ,ebert & curry (1992) method
!                   ! =2 ,streamer v3 (2001) method
!                   ! =3 ,fu (1996) method
!     icliq_sw = 1  ! sw cloud optical property for cloud water
!     icliq_lw = 1  ! lw cloud optical property for cloud water
!                   ! =0 ,diagnostic cld opt depth
!                   ! =1 ,hu & stamnes (1993) method
!     idate(8)      ! ncep absolute date and time of initial conditions
!     iflip    = 0  ! = 0 ; index from toa to surface
!              = 1  ! = 1 ; index from surface to toa
!     me       = 0  ! = 0 ; print out control flag 'on'
!              = 1  ! = 1 ; print out control flag 'off'
!     ioutsigr = 0  ! = 0 ; output radiation sigma parameters
!              = 1  ! = 1 ; not output
!     sashal       = .true
!     crick_proof  = .false.
!     ccnorm       = .false.
!     norad_precip = .false.   ! This is effective only for
!     Ferrier/Moorthi
!
!---------------------------------------------------------------------------
      levr=lev
      ictm=1
      isol=1
      ico2=2
      iaer=111
      ialb=0
      iems=1
      ntcw=2
      num_p3d=4
      me=1
!---------------------------------------------------------------------------
      irad=2
      iflip=1
      ntiw=3
      ntrw=4
      ntsw=5
      ntgl=6
      ntoz=3
      nclds=1
      ioutsigr=0
!---------------------------------------------------------------------------
      idate(1)=0
      idate(2)=0
      idate(3)=0
      idate(4)=0
      idate(5)=0
      idate(6)=0
      idate(7)=0
      idate(8)=0
!---------------------------------------------------------------------------
      sashal=.true.
      crick_proof=.false.
      ccnorm=.false.
      norad_precip=.false.
!---------------------------------------------------------------------------
      return
      end
