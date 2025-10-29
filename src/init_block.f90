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
      use leapyr

      implicit none

      cp   = con_cp
      rad  = con_rerth
      omega= con_omega
      grav = con_g
      stbo= con_sbc
      s0  = con_solr
      hltm= con_hvap
      tice= con_tice
      hice= con_hfus
!!      cp=1004.24
!!      rad=6.371e6
!!      omega=7.292e-5
!!      grav=9.80616

      dt=900.0

!  order of horozontal diffusion
      hord=4
!
      frad=1.0
      ptmean=1000.0
      ksgeo=3
      ktcup=16

      ktpbl=16
      ktshl=16
      njump=2
      ldiag=1
      idg=16
      jdg=16
      qmin=1.0e-20

!!      stbo=5.669e-8
!!      s0=1368.3
!!      hltm=2.52e6
!!      tice=273.15

!!      hice=3.336e5
      evaprh=0.98
      hfilt=1.
      factop=60.

      nnmiit=3
      nnmivm=3
      cutfreq=1.0
      itypbl=0

      taup=6.
      ktrop=8
      taureg=6.
!
      lsimpl=.true.
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
! about of input and output set
      domfc    =    384.
      otgreen  =      6.
      out_green= .false.
      out_hp   = .false.
      outgrb2  =      0  !output grib2 format
      outdms   =      1  !output dmskey 
      outfv3   = .false.
! pdf cloud
      pdfcloud=.false.
! stochastic physics
      dosppt=.false.
      dospptout=.false.
      doshum=.false.
      doshumout=.false.
! update low boundary condition
      doclx=.false.
! Semi-Lagrangian Averaging of Physical Parametrizations
      doslavepp=.false.
! dy-core two loop sequence
      two_loop=.false.
! dy-core two time level
      ttl=.true.
      itter=2
! dry air mass correction
      mass_dp=.false.
      dpprt  =.false.
! output data for RSM (Also, RSM compiling flag is necessary)
      outrsm=.false.
      rsmoutinv=6
      rlon1=100.
      rlon2=150.
      rlat1=5.
      rlat2=40.
      rgrdsz=0.25
      rsmsfcmgrhr=24
!---------------------------------------------------------------------------
!
! specify the default option for cup and pbl
!
      nmcup=6
      nmpbl=4
      nmland=2
      nmshl=3

      cgw=1.0e-4
!
!---for using forecast daily sst, sea ice fraction, snow depth
      ldailyFCTsst=.false.
      ldailyFCTicesndpt=.false.
      lFCTweight=.false.
      mom4ice=.false.
      lopgsst=.false.
      dailyClm_option=-99
!--for SIT
      do_sit=.false.
      fsit=-99.
      dSITdt_intv=-99.
      weightSIT=1.
      updatetg=24.
!
! specify the default option for orographic and convective gwd
!
      nmgwor=2
      tofd=.false.
      nmgwcv=2
      mtnvar=14
      cmbk = 1.0
      cgwd = 1.2
!-for Cloud Micro Physics
      nmmiph=12
      ntinc=7   ! tracer index for ice number concentration
      ntrnc=8   ! tracer index for rain number concentration
!     ntlnc=9   ! tracer index for liquid number concentration
! for using aerosol climatology
      naero = 1
! for GCE 3ice
      SL_sedi = .false.
      sat_predict = .true.
      new_saturation = .true.
      use_cpm = .false.
      use_declination = .false.


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
      ptop=0.1
!  sponge layer 
      spl1=10.
      spl2=100.
      vd=0.01
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
      ptop=0.1
!  sponge layer 
      spl1=5.
      spl2=50.
      vd=0.01
      else if ( lev .eq. 128 ) then
!
! L128 hybrid coordinate
   tmean=(/ 272.07000,272.07000,272.07000,272.07000,272.07000,272.07000, &
            272.07000,272.07000,272.07000,272.07000,272.07000,272.07000, &
            270.87453,266.82119,262.94369,259.24990,255.53608,251.95015, &
            248.51807,245.20174,242.11459,239.25858,237.08906,235.17648, &
            233.34808,231.61195,229.94067,228.32980,227.31186,226.67073, &
            226.04708,225.43891,224.84470,224.26301,223.69289,223.13332, &
            222.58351,222.04276,221.51039,220.98583,220.46863,219.95843, &
            219.45476,218.95740,218.46592,217.98024,217.50008,217.02529, &
            216.65000,216.65000,216.65000,216.65000,216.65000,216.65000, &
            216.65000,216.65000,216.65000,216.65000,216.65000,216.65000, &
            216.65000,216.65000,216.65000,216.65000,216.65000,216.65000, &
            216.65000,216.65000,216.65000,216.98519,218.79468,221.16790, &
            223.52246,225.85513,228.16483,230.44784,232.70132,234.92467, &
            237.11320,239.26497,241.37789,243.45012,245.47867,247.46152, &
            249.39684,251.28284,253.11786,254.90058,256.62974,258.30419, &
            259.92348,261.48677,262.99365,264.44381,265.83740,267.17456, &
            268.45597,269.68204,270.85337,271.97115,273.03641,274.05008, &
            275.01374,275.92862,276.79607,277.61809,278.39579,279.13099, &
            279.82524,280.48027,281.09786,281.67957,282.22703,282.74192, &
            283.22576,283.68006,284.10634,284.50642,284.88135,285.23245, &
            285.56147,285.86927,286.15730,286.42657,286.67835,286.91348, &
            287.13331,287.33843 /)
!
! L128: hybrid coordinate
    aki=(/   0.00000,  0.02197,  0.03378,  0.05082,  0.07487,  0.10809,  &
             0.15305,  0.21268,  0.29024,  0.38929,  0.51353,  0.66674,  &
             0.85264,  1.07473,  1.33621,  1.63986,  1.98794,  2.38224,  &
             2.82398,  3.31398,  3.85267,  4.44029,  5.07708,  5.76344,  &
             6.50025,  7.28907,  8.13245,  9.03427, 10.00000, 11.03644,  &
            12.14921, 13.34374, 14.62584, 16.00168, 17.47782, 19.06127,  &
            20.75943, 22.58019, 24.53192, 26.62347, 28.86419, 31.26399,  &
            33.83330, 36.58312, 39.52499, 42.67106, 46.03402, 49.62714,  &
            53.46428, 57.55980, 61.92864, 66.58620, 71.54836, 76.82692,  &
            82.39076, 88.18442, 94.15203,100.23763,106.38533,112.53935,  &
           118.64421,124.64488,130.48688,136.11654,141.48113,146.52906,  &
           151.21015,155.47579,159.27919,162.57562,165.32260,167.48016,  &
           169.01158,169.90105,170.15713,169.79275,168.82412,167.27074,  &
           165.15541,162.50417,159.34622,155.71373,151.64165,147.16742,  &
           142.33063,137.17270,131.73642,126.06555,120.20433,114.19705,  &
           108.08753,101.91870, 95.73212, 89.56761, 83.46284, 77.45304,  &
            71.57067, 65.84527, 60.30320, 54.96762, 49.85837, 44.99199,  &
            40.38177, 36.03782, 31.96720, 28.17408, 24.65992, 21.42367,  &
            18.46199, 15.76946, 13.33883, 11.16123,  9.22638,  7.52284,  &
             6.03817,  4.75916,  3.67197,  2.76235,  2.01571,  1.41732,  &
             0.95190,  0.60249,  0.35193,  0.18327,  0.07986,  0.02538,  &
             0.00393,  0.00001,  0.00000 /)

    bki=(/ .00000000,.00000000,.00000000,.00000000,.00000000,.00000000,  &
           .00000000,.00000000,.00000000,.00000000,.00000000,.00000000,  &
           .00000000,.00000000,.00000000,.00000000,.00000000,.00000000,  &
           .00000000,.00000000,.00000000,.00000000,.00000000,.00000000,  &
           .00000000,.00000000,.00000000,.00000000,.00000000,.00000000,  &
           .00000000,.00000000,.00000000,.00000000,.00000000,.00000000,  &
           .00000000,.00000000,.00000000,.00000000,.00000000,.00000000,  &
           .00000000,.00000000,.00000000,.00000000,.00000000,.00000000,  &
           .00000000,.00000000,.00000000,.00000000,.00000000,.00000447,  &
           .00006115,.00024241,.00062120,.00127063,.00226373,.00367320,  &
           .00557114,.00802871,.01111581,.01490065,.01944939,.02482570,  &
           .03109023,.03830025,.04650905,.05576552,.06611363,.07759200,  &
           .09023283,.10404410,.11900588,.13509022,.15226209,.17047936,  &
           .18969276,.20984599,.23087597,.25271315,.27528202,.29850170,  &
           .32228657,.34654717,.37119100,.39612350,.42124904,.44647195,  &
           .47169749,.49683292,.52178840,.54647788,.57081995,.59473849,  &
           .61816332,.64103066,.66328347,.68487175,.70575261,.72589030,  &
           .74525610,.76382818,.78159126,.79853634,.81466026,.82996528,  &
           .84445864,.85815203,.87106113,.88320512,.89460619,.90528910,  &
           .91528070,.92460955,.93330554,.94139951,.94892295,.95590771,  &
           .96238622,.96839253,.97396074,.97912445,.98391663,.98836950,  &
           .99251445,.99638192,1.0000000 /) 
      ptop=0.01
!  sponge layer 
      spl1=1.
      spl2=50.
      vd=0.005
      endif
!
      tmeans=300.
!      tmeans=350.

!-- for hybrid coordinates, ptmeans reset for numerical stability
!      ptmeans=800.
      ptmeans=600.
!
! for forward weighting Semi-Implicit
!
      alpha=0.7
!
! for Robert time filter in three time level
!
      tfilt=0.04
!
! for two time level 
!
!    coefficient of horizontal difussion for mid-point wind
!
      mwhd=1.

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
!-- for sit
      ifilin_ncep   = 'ifilin_ncep'
      ifilin_sst    = 'ifilin_sst'
      ifilin_nc     = '.' ! 'ifilin_nc'  change to path
      ifilin_ClmANA = 'ifilin_ClmANA'
      ifilin_ClmFCT = 'ifilin_ClmFCT'
!-- for grib2 output path
      ifilout_grb  = '.'
!-- for aerosol climatology
      ifilin_aero   = 'ifilin_aero'
!
!dms
!t512l60
      ggdef='gh0g'
      gmdef='ghmg'
      gsdef='gh0s'
!soil
!      tsat=(/.421,.464,.468,.434,.406,.465,.404,.439,.421/)
!      ref=(/.283,.387,.412,.312,.338,.382,.315,.329,.283/)
!      wlt=(/.029,.119,.139,.047,.100,.103,.069,.066,.029/)
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
!              = 5  ! num_p3d=5, GFDL Microphysics cloud scheme with effective radii
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
      me=1
!---------------------------------------------------------------------------
      irad=2
      iflip=1
      iovr_sw  = 1  ! sw: maximum-random overlapping vertical cloud layer
      iovr_lw  = 1  ! lw: maximum-random overlapping vertical cloud layer
      isubc_sw = 0  ! sw with mcica sub-col approximation provided random seed
      isubc_lw = 0  ! lw with mcica sub-col approximation provided random seed
      icice_sw = 3  ! sw cloud optical property for cloud ice
      icice_lw = 3  ! lw cloud optical property for cloud ice
      icliq_sw = 1  ! sw cloud optical property for cloud water
      icliq_lw = 1  ! lw cloud optical property for cloud water
      ntrw=3
      ntiw=4
      ntsw=5
      ntgl=6
      nthl=7  ! hail
      ntoz=3
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
! for leapyear
      leap=.false.
      leapm1=.false.
      yrd=365
!---------------------------------------------------------------------------
      sashal=.true.
      crick_proof=.false.
      ccnorm=.false.
      norad_precip=.false.
!---------------------------------------------------------------------------
! for aeroclx index
      naso4 = 1  ! Sulphate
      nadu1 = 2  ! Dust bin 001
      nadu2 = 3  ! Dust bin 002
      nadu3 = 4  ! Dust bin 003
      nadu4 = 5  ! Dust bin 004
      nadu5 = 6  ! Dust bin 005
      nass1 = 7  ! Sea Salt bin 001
      nass2 = 8  ! Sea Salt bin 002
      nass3 = 9  ! Sea Salt bin 003
      nass4 = 10 ! Sea Salt bin 004
      nass5 = 11 ! Sea Salt bin 005
      nablc = 12 ! Hydrophilic Black Carbon
      nabbc = 13 ! Hydrophobic Black Carbon
      naolc = 14 ! Hydrophilic Organic Carbon
      naobc = 15 ! Hydrophobic Organic Carbon
      namsa = 16 ! Methanesulphonic acid
      nadms = 17 ! Dimethylsulphide
      naso2 = 18 ! Sulphur dioxide
!---------------------------------------------------------------------------
      return
      end
