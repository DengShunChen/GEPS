  module const
! modify to f90 by C-H Lee and sort by River Chen in 2015
    use param
    use mpi,   only : MPI_REAL4, MPI_REAL8

    implicit none
 
    public

!CWB2021 for single precison test
#ifdef SP
      integer, parameter ::     RTYPE=4
      integer, parameter :: MPI_RTYPE=MPI_REAL4
      character(len=1), parameter  :: kflag='R'
#else
      integer, parameter ::     RTYPE=8
      integer, parameter :: MPI_RTYPE=MPI_REAL8
      character(len=1), parameter  :: kflag='H'
#endif
 
    real(kind=RTYPE), dimension(:)  , allocatable, save  :: aki,bki
    real(kind=RTYPE), dimension(:,:), allocatable, save  :: sigma,dsigma
 
    integer, allocatable, save ::  mlsort(:,:)
    integer, allocatable, save ::  msort(:),lsort(:)
 
    integer :: numout,ipadding,jm2,ksgeo,                     &
            ktpbl,ktshl,ktcup,julian,ldiag,idg,jdg,njump,  &
            nnmiit,nnmivm,itypbl,                          &
            nmgwor,nmgwcv,mtnvar,                          &
            ktrop,ncpu,nmcup,nmpbl,nmland,numreduce,nmshl, &
            nmmiph
 
    common/constI/                                         &
            numout,ipadding,jm2,ksgeo,                     &
            ktpbl,ktshl,ktcup,julian,ldiag,idg,jdg,njump,  &
            nnmiit,nnmivm,itypbl,                          &
            nmgwor,nmgwcv,mtnvar,                          &
            ktrop,ncpu,nmcup,nmpbl,nmland,numreduce,nmshl, &
            nmmiph
 
!    real, dimension(:), allocatable, save  ::              &
!         weight,sinl,cosl,cor,onocos
!         tmean,spalm,eigval,pmcor,tmeans
 
!    real, dimension(:,:), allocatable, save  :: evecin,    &
!         evectr,arrhyd,arsddt,tmcor

    real(kind=RTYPE), dimension(:), allocatable, save ::   &
         weight,sinl,cosl,cor,onocos,                      &
         tmean,spalm,eigval,pmcor,tmeans
    real(kind=RTYPE), dimension(:,:), allocatable, save :: &
         evecin,evectr,arrhyd,arsddt,tmcor
!
    real ::                                                &
         capa,cp,rad,radsq,grav,omega,rgas,stbo,s0,hltm,   &
         ptop,dt,tau,taui,taue,tauo,                       &
         hours,frad,evaprh,qgini,                          &
         tice,hice,cutfreq,taup,hfilt,                     &
         taureg,cgw,domfc,otgreen,cgwd,cmbk,spl1,spl2
    real(kind=RTYPE) :: ptmean,ptmeans,qmin
    !sit
    real :: fsit         !fsit>0., turn on sit_vdiff when mod(tau/fsit)<0.001
                         !default fsit<=0., turn on sit_vdiff every tau
    real dSITdt_intv    !interval of returning ave. dSITdt to dSST/dt
    real weightSIT      !the weighting of SIT tendency
    real updatetg       !the time interval to update tg 
          
 
    common/constR/                                         &
         capa,cp,rad,radsq,grav,omega,rgas,stbo,s0,hltm,   &
         ptop,ptmean,dt,tau,taui,taue,tauo,                &
         hours,frad,evaprh,qgini,                          &
         tice,hice,cutfreq,taup,hfilt,ptmeans,             &
         taureg,cgw,fsit,domfc,otgreen,spl1,spl2,          &
         dSITdt_intv,weightSIT,updatetg
    logical :: lsimpl,lzadv, yesdia,dopbl, docup, dorad,   &
            dolsp, dograv,doshl, dodry, donnmi,ozon,       &
            restrt,hdiff, cstar, update,doincr,hybrid,     &
            doo3l, docgrav, doclx, tofd, doslavepp,        &
            two_loop

    ! for stochastic physics
    logical :: dosppt       =.false.
    logical :: dospptout    =.false.
    logical :: doskebout    =.false.
    logical :: doshum       =.false.
    logical :: doskeb       =.false.
    logical :: dossst       =.false.
    logical :: use_zmtnblck =.false.
         
    logical :: out_green,out_hp

    !for Semi-Lagrangain
    logical :: ndsladvh2

! output data for RSM (Also, RSM compiling flag is necessary)
    !for RSM output
    logical :: outrsm
    integer :: rsmoutinv
    real    :: rlon1, rlon2, rlat1, rlat2, rgrdsz

    !for horizontal diffusion
    integer :: hdk1,hdk2(3),hord
    real    :: vd

    !for pdf cloud
    logical :: pdfcloud

    !for 2dMPI
    logical :: idg_jdg_owner
    integer :: idg_listnum,jdg_listnum

    ! daily forecast sst, sea ice fraction, water equivlent snow depth, time weighting
    logical :: ldailyFCTsst,ldailyFCTicesndpt,lFCTweight
    integer :: dailyClm_option
    logical :: lopgsst

    ! sit
    logical :: do_sit

    ! SKEB
    logical :: first_call

    !for output 
    integer :: outgrb2    !output grib2 format
    integer :: outdms     !output dmskey

    common/constL/lsimpl,lzadv,yesdia,dopbl,docup,dorad,   &
            dolsp, dograv,doshl, dodry, donnmi,ozon,       &
            restrt,hdiff, cstar, update,doincr,hybrid,     &
            doo3l,ndsladvh2,docgrav,out_green,out_hp,      &
            ldailyFCTsst,ldailyFCTicesndpt,lFCTweight,     &
            dailyClm_option,lopgsst,do_sit,tofd


    character(len=80) ifilin,cwbout,bckfile,namlsts, &
            ifilout,crdate,ocards,phyout,cntrl, &
            ifilin_ncep,ifilin_sst,ifilin_nc,   &
            ifilin_ClmANA,ifilin_ClmFCT,ifilout_grb

    common/files/ifilin,cwbout,bckfile,namlsts, &
            ifilout,crdate,ocards,phyout,cntrl, &
            ifilin_ncep,ifilin_sst,ifilin_nc,   &
            ifilin_ClmANA,ifilin_ClmFCT

    character(len=16), dimension(:), allocatable, save  :: outdir

    !dms34
    integer(kind=8) ::  idtg,idtg2
    common/constI8/idtg,idtg2

    !dms34
    !ggdef, gmdef, gsdef : grid system definition for gg, gm, gs
    character(len=4)  ::  ggdef,gmdef,gsdef
    common/dmskey34/ggdef,gmdef,gsdef

    real(kind=RTYPE), dimension(:,:,:), allocatable, save  :: poly,dpoly
    real(kind=RTYPE), dimension(:,:)  , allocatable, save  :: eps4,wdfac,wcfac
    real(kind=RTYPE), dimension(:)    , allocatable, save  :: cim,eps4L   ! for 2dMPI

    contains 

      subroutine allocate_const_array
        implicit none
        integer  ierr

        allocate (poly(jtrun,my/2,jtmax),dpoly(jtrun,my/2,jtmax),          &
                  eps4(jtrun,jtmax),wdfac(jtrun,jtmax),wcfac(jtrun,jtmax), &
                  cim(jtmax), stat=ierr)
        if (ierr/= 0) then
            write(6,*) 'mod_const : allocate fail 1 '
            stop
        end if

        allocate (aki(lev+1),bki(lev+1),                         &
                  sigma(lev+1,2),dsigma(lev,2), stat= ierr)
        if (ierr/= 0) then
            write(6,*) 'mod_const : allocate fail 2 '
            stop
        end if

        allocate (mlsort(jtrun,jtrun),msort(mlmax),lsort(mlmax), &
                  stat= ierr)
        if (ierr/= 0) then
            write(6,*) 'mod_const : allocate fail 3 '
            stop
        end if

        allocate (weight(my),sinl(my),cosl(my),           &
        cor(my),onocos(my),                               &
        tmean(lev),spalm(lev),eigval(lev),evecin(lev,lev),&
        evectr(lev,lev),arrhyd(lev,lev),arsddt(lev,lev),  &
        pmcor(lev),tmcor(lev,lev),tmeans(lev),            &
                  stat= ierr)
        if (ierr/= 0) then
            write(6,*) 'mod_const : allocate fail 4 '
            stop
        end if


        allocate (outdir(nout),stat= ierr)
        if (ierr/= 0) then
            write(6,*) 'mod_const : allocate fail 5 '
            stop
        end if
      end subroutine

      subroutine deallocate_const_array
        implicit none
     
        deallocate(poly,dpoly,eps4,wdfac,wcfac,cim)
        deallocate(aki,bki,sigma,dsigma)
        deallocate(mlsort,msort,lsort)
        deallocate(weight,sinl,cosl,cor,onocos,          &
                   tmean,spalm,eigval,pmcor,tmeans)
        deallocate(outdir)
      end subroutine
  end module const
