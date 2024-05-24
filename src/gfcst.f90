      program gfcst
!
! main program of CWBGFS
! modify to f90 bt C-H Lee and sort by River Chen in 2015
!
      use param
      use const
!
      implicit none

      integer  no
!
!  logical io units:
!
!  input namelist file = 'namlsts'
!  date-time-group file = 'crdate'
!  output directives file = 'ocards'
!  model history file = 'cwbout'
!  output model history for diabatic variables='phyout'
!  input file of path/file names='filist'
!

      call mpe_init
!
!     get model constants
!
      call cons
!
!  read in initial data and prepare for initialization/forecast
!
      call getrdy
!
!  initialization would be done here, taking initial spectral fields
!  out of 'getrdy' and preparing them for 'intgrt'.  for identification
!  purposes only, initialization history file data is assigned
!  "tau"=1.0.
!
      if (donnmi.and.taui.lt.1.0) then
         no=2*((jtrun+1)/2)+(jtrun/2)+10
         call initial(no,jtrun,jtmax,lev,nx,my,my_max,mlmax)
      endif
!
!  recompute eigen value for semi-implicit scheme
!
      call matrix_hybrid_cwb ( cp,sigma,dsigma,ptop,ptmeans,tmeans,spalm   &
                  , eigval,evecin,evectr,arrhyd,arsddt,pmcor,tmcor )
!
!  time integration
!
      if ( ttl ) then
#ifdef USE_CUDA
        call intgrt_gpu
#else
        call intgrt
#endif
      else
        call intgrt_3tl
      endif

!
      call mpe_finalize
      call dmsexit(0)
!
      stop
      end
