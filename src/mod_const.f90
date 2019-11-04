      module const
!
! modify to f90 by C-H Lee and sort by River Chen in 2015
!
      use param

      implicit none

      public

      real, dimension(:)  , allocatable, save  :: aki,bki
      real, dimension(:,:), allocatable, save  :: sigma,dsigma
 
      integer, allocatable, save ::  mlsort(:,:)
      integer, allocatable, save ::  msort(:),lsort(:)

      integer numout,ipadding,jm2,ksgeo,                     &
              ktpbl,ktshl,ktcup,julian,ldiag,idg,jdg,njump,  &
              nnmiit,nnmivm,itypbl,                          &
              nmgwor,nmgwcv,mtnvar,                          &
              ktrop,ncpu,nmcup,nmpbl,nmland,numreduce,nmshl

      common/constI/                                         &
              numout,ipadding,jm2,ksgeo,                     &
              ktpbl,ktshl,ktcup,julian,ldiag,idg,jdg,njump,  &
              nnmiit,nnmivm,itypbl,                          &
              nmgwor,nmgwcv,mtnvar,                          &
              ktrop,ncpu,nmcup,nmpbl,nmland,numreduce,nmshl

      real, dimension(:), allocatable, save  ::              &
           weight,sinl,cosl,cor,onocos,sig,dsig,             &
           tmean,spalm,eigval,pmcor,tmeans

      real, dimension(:,:), allocatable, save  :: evecin,    &
           evectr,arrhyd,arsddt,tmcor

      real                                                   &
           capa,cp,rad,radsq,grav,omega,rgas,stbo,s0,hltm,   &
           ptop,ptmean,tfilt,dt,tau,taui,taue,tauo,          &
           hours,frad,evaprh,qgini,                          &
           tice,hice,cutfreq,taup,hfilt,ptmeans,             &
           taureg,cgw,domfc,otgreen
!sit
      real fsit           !fsit>0., turn on sit_vdiff when mod(tau/fsit)<0.001
                          !default fsit<=0., turn on sit_vdiff every tau
           
 
      common/constR/                                         &
           capa,cp,rad,radsq,grav,omega,rgas,stbo,s0,hltm,   &
           ptop,ptmean,tfilt,dt,tau,taui,taue,tauo,          &
           hours,frad,evaprh,qgini,                          &
           tice,hice,cutfreq,taup,hfilt,ptmeans,             &
           taureg,cgw,fsit,domfc,otgreen
! sppt parameters
      real                                                       &
           de_corretime_500,de_corretime_1000,de_corretime_2000, &
           facsppt500,facsppt1000,facsppt2000
           

 
      logical lsimpl,lzadv, yesdia,dopbl, docup, dorad,      &
              dolsp, dograv,doshl, dodry, donnmi,ozon,       &
              restrt,hdiff, cstar, update,doincr,hybrid,     &
              doo3l, dosppt, dospptout,   docgrav
      logical out_green,out_hp
!for Semi-Lagrangain
      logical ndsladvh2

!for horizontal diffusion
      integer hdktop,hdk1,hdk2,hdk3
      real    factop,coefu
!for pdf cloud
      logical pdfcloud

!for 2dMPI
      logical idg_jdg_owner
      integer idg_listnum,jdg_listnum
! daily forecast sst, sea ice fraction, water equivlent snow depth, time weighting
      logical ldailyFCTsst,ldailyFCTicesndpt,lFCTweight,lday_chtg
      integer dailyClm_option
      logical lopgsst
! sit
      logical do_sit
      logical ltgtest

      common/constL/lsimpl,lzadv,yesdia,dopbl,docup,dorad,   &
              dolsp, dograv,doshl, dodry, donnmi,ozon,       &
              restrt,hdiff, cstar, update,doincr,hybrid,     &
              doo3l,ndsladvh2,docgrav,out_green,out_hp,      &
              ldailyFCTsst,ldailyFCTicesndpt,lFCTweight,     &
              dailyClm_option,lopgsst,do_sit,ltgtest,lday_chtg

 
      character*80 ifilin,cwbout,bckfile,namlsts, &
              ifilout,crdate,ocards,phyout,cntrl, &
              ifilin_ncep,ifilin_sst,ifilin_nc,   &
              ifilin_ClmANA,ifilin_ClmFCT

      common/files/ifilin,cwbout,bckfile,namlsts, &
              ifilout,crdate,ocards,phyout,cntrl, &
              ifilin_ncep,ifilin_sst,ifilin_nc,   &
              ifilin_ClmANA,ifilin_ClmFCT
 
      character(len=16), dimension(:), allocatable, save  :: outdir

!dms34
      integer*8      idtg,idtg2
      common/constI8/idtg,idtg2

!dms34
!ggdef, gmdef, gsdef : grid system definition for gg, gm, gs
      character*4     ggdef,gmdef,gsdef
      common/dmskey34/ggdef,gmdef,gsdef

      real, dimension(:,:,:), allocatable, save  :: poly,dpoly
      real, dimension(:,:)  , allocatable, save  :: eps4,wdfac,wcfac
      real, dimension(:)    , allocatable, save  :: cim
      real, dimension(:)    , allocatable, save  :: eps4L   ! for 2dMPI

      contains 

         subroutine allocate_const_array

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
           cor(my),onocos(my),sig(lev+1),dsig(lev),          &
           tmean(lev),spalm(lev),eigval(lev),evecin(lev,lev),&
           evectr(lev,lev),arrhyd(lev,lev),arsddt(lev,lev),  &
           pmcor(lev),tmcor(lev,lev),tmeans(lev),            &
                     stat= ierr)

           if (ierr/= 0) then
               write(6,*) 'mod_const : allocate fail 4 '
               stop
           end if

!CWB2015
           sig=0.

           allocate (outdir(nout),stat= ierr)

           if (ierr/= 0) then
               write(6,*) 'mod_const : allocate fail 5 '
               stop
           end if

           return

         end subroutine

         subroutine deallocate_const_array

           deallocate(poly,dpoly,eps4,wdfac,wcfac,cim)
           deallocate(aki,bki,sigma,dsigma)
           deallocate(mlsort,msort,lsort)
           deallocate(weight,sinl,cosl,cor,onocos,sig,dsig,  &
                      tmean,spalm,eigval,pmcor,tmeans)
           deallocate(outdir)

           return

         end subroutine

      end module const
