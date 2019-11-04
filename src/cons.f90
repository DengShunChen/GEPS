      subroutine cons
!
!***********************************************************************
!  this subroutine defines several important constants and arrays
!  used by the spectral forecast model.  they include physical
!  parameters, model vertical structure, matrix operators used in
!  the semi-implicit algorthm, and polynomial arrays used in the
!  spherical harmonic transforms.
!
!
! **** output: ***
!
!  restrt: logical variable set to .true. if initial tau is greater
!          than zero.
!
!  modify to f90 in 2015 by C-H Lee and sort by River Chen in 2015
!***********************************************************************
!
      use param
      use mpe
      use rank
      use index
      use const
      use spec
      use fftcom
      use mod_typhoon
      use mod_sit_control,       ONLY:sit_nml
!-----------------------------------------------------------------------
      use radn
      use noah
!-----------------------------------------------------------------------
!

! for ECHAM4 Tiedtke cumulus scheme
      USE mo_cumulus_flux, only : cuparam
      USE mo_constants,    only : inicon
      USE mo_convect_tables, only : set_lookup_tables
! for WSM6
      use module_mp_wsm6, only : wsm6init

      implicit  none

      integer i,j,k,n,jj,nn,istat,istat1,istat2,istat3,irstat,ii,Wntyph,nc,Wltyph
      integer ix,jy,ip,istat_r,istat_w,ierror,ltyph,io
      integer ifromtau,itotau,itau,itaui,l,nxj,m,mf,my2
      integer tflag,Wflag

      real    pi,rm,rl,rlm,one,onem,r2d, pnm_max,pnmcut,sumreduce
      real    reducefactor,d2r,cew,clon,cns,clat

      integer*8 idtg8

      real      pnm(jtrun+1,jtrun+1)

!
      namelist /modlst/ ksgeo,ptop,ptmean,tfilt,dt,taui,taue            &
                      , tauo,frad,ktpbl,ktshl,ktcup,njump,evaprh,lsimpl &
                      , lzadv,yesdia,dopbl,docup,dorad,dolsp,dograv     &
                      , doshl,dodry,donnmi,idg,jdg,ldiag,nnmiit,nnmivm  &
                      , cutfreq,hdiff,itypbl,cstar,taup,hfilt           &
                      , ptmeans,update,taureg,doincr,numreduce          &
                      , nmcup,nmpbl,nmland,nmshl,cgw,ggdef,gmdef        &
                      , nmgwor,nmgwcv,mtnvar,docgrav                    &
                      , ictm,isol,ico2,iaer,ialb,irad,iems,ntcw         &
                      , num_p3d,ntoz,iovr_sw,iovr_lw,isubc_sw,isubc_lw  &
                      , sashal,crick_proof,ccnorm,norad_precip,me,doo3l &
                      , ioutsigr,domfc,out_green,isot,ivegsrc           &
                      , otgreen,out_hp,dosppt,dospptout                 &
                      , de_corretime_500,de_corretime_1000              &
                      , de_corretime_2000                               &
                      , facsppt500,facsppt1000,facsppt2000,ndsladvh2    &
                      , ldailyFCTsst,ldailyFCTicesndpt,lFCTweight       &
                      , dailyClm_option,lopgsst,do_sit,fsit,pdfcloud    &
                      , ltgtest,lday_chtg
!
      real    si(lev+1)
      logical flag
      character*10 fulldtg,Wfulldtg
      character*80 filist
      character cdtg*12
      character*80 pathname,logicname,truefile
      character*64 type_r,type_w,argument
      namelist /filst/ ifilin,cwbout,bckfile,namlsts &
                     , ifilout,crdate,ocards,phyout,cntrl &
                     , ifilin_ncep, ifilin_sst, ifilin_nc &
                     , ifilin_ClmANA,ifilin_ClmFCT

      namelist /typ/ write_tau, write_mem, trk_intv, min_trk_pres
      integer istat4,istat5,istat6,istat7

      data pathname/'NWPETCGLB'/
      data logicname/'filist'/
!

! for ECHAM4 Tiedtke cumulus scheme
      call cuparam
      call inicon
      call set_lookup_tables

      do 100 k = 2, lev
      dsig(k-1)= sig(k) - sig(k-1)
  100 continue
!
      dsig(lev)= 1.0 - sig(lev)
!
      capa= 1.0/3.5
      rgas= capa*cp
      pi  = 4.0*atan(1.0)
      radsq= rad*rad
!
      call gpvs
!----------------------------------------------------------------!
!
!  read namlist of path/file name(operation)
!
      call getfname(pathname,logicname,truefile,istat)
      if(istat.ne.0)then
        print *,'getfname : error','RANK=',myrank
        call mpe_finalize
        call dmsexit(-1)
      else
        if(myrank .eq. 0) print *,truefile
      endif
!
      open (unit=12,file=trim(truefile),form='formatted')
!!      open (unit=77,file='output.dat',form='formatted')
!
      read (12,filst,end=110)
!
  110 continue
      close(12)
!
      if(myrank .eq. 0) print filst
!
!  read namelist for model parameters
!
      open (unit=1,file=trim(namlsts),form='formatted')
!
      read (1,modlst,end=120)
  120 continue
      read (1,typ,end=121)
  121 continue
      close(1)
!
      close(1)
      open (unit=1,file=trim(namlsts),form='formatted')

      if(do_sit) then
        read (1,sit_nml,end=130)
        if(myrank .eq. 0) then
          print *, 'chlee debug...'
          print sit_nml
        endif
  130 continue
      endif
!
      open (unit=2,file=trim(crdate),form='formatted')
!
! read in idtg*12
      read(2,'(i8.8)')idtg8
      close(2)
! transfer idtg8 to idtg*12
      if(idtg8.gt.60000000)then
        idtg = 200000000000 + idtg8*100
      else
        idtg = 200000000000 + idtg8*100
      endif
!
      write(cdtg,900)idtg
 900  format(i12.12)
!
      if(myrank .eq. 0) print*,' dtg=',cdtg
!-----------------------------------------------------------------------
      read(cdtg,'(i4,i2,i2,i2,i2)') idate(1),idate(2),idate(3),idate(5)&
                                   ,ii
!-----------------------------------------------------------------------
!
      call days (cdtg,julian,hours)
!
      flag =.false.
      if(myrank .eq. 0)then
        call recmsg('gfs',ifromtau,itotau,istat)
        flag =.true.
      endif
!ch   call mpe_broadcast(istat,1,flag,mpe_integer)
      call mpe_bcast(istat,1,0,mpe_integer)

      if(istat.eq.-1)then
        if(myrank.eq.0) print *,'RECMSG ERROR','RANK=',myrank
        call mpe_finalize
        call dmsexit(-1)
      else if(istat.eq.1)then
        call mpe_finalize
        call dmsexit(0)
      endif
!
!ch   call mpe_broadcast(ifromtau,1,flag,mpe_integer)
!ch   call mpe_broadcast(itotau,1,flag,mpe_integer)
      call mpe_bcast(ifromtau,1,0,mpe_integer)
      call mpe_bcast(itotau,1,0,mpe_integer)
!
      taui=float(ifromtau)
      taue=float(itotau)
      tauo=float(itotau)+0.000001
!
      tau   = 0.0
      restrt=.false.
      if (taui.gt.0.0)  then
        itaui=ifromtau
        itau= int(taui) + 0.001
        itau= min(itau,itaui)
        taui= float(itau)
        restrt=.true.
        if(myrank .eq. 0) print*,' restarting at tau=',itau
        hours = hours+taui
        julian= julian+hours/24.0+0.001
        hours = mod(hours,24.)
      endif
!
      if(myrank .eq. 0) print modlst
      if(myrank .eq. 0) print typ
!
      if (taui .ge. taue)  then
        if(myrank .eq. 0)  &
         print *,'******** aborting run, taui.ge.taue ******* '
        stop
      endif

!
!  build pointer arrays for locating zonal and total wavenumber
!  values in the one-dimensional spherical harmonic arrays.
!
!  original allocation system only for spectral space
!
      call sortml (jtrun,mlmax,msort,lsort,mlsort)
!
!  new allocation system in parallization for both spectral 
!  and grid spaces
!

      call make_list

      do 150 m =1,mlistnum
       mf=mlist(m)
!
       rm=mf-1
       if (mf.eq.1)  rm= 0.0
       cim(m) = rm
!
      do 150 l=mf,jtrun
       rl = l
       rlm= rl-1.0
       if (mf.eq.1)  rm= 0.0
       eps4(l,m)= rl*rlm/radsq
       if ( l.eq.1 ) then
        wdfac(1,m) = 0. ; wcfac(1,m) = 0.
       else
        wdfac(l,m) = 1.0/(radsq*eps4(l,m))
        wcfac(l,m) = cim(m)*wdfac(l,m)
       endif
  150 continue
!
      do k=1,lev+1
        sigma(k,1) = bki(k)
        sigma(k,2) = aki(k)
      enddo
      do k=1,lev
        dsigma(k,1) = sigma(k+1,1) - sigma(k,1)
        dsigma(k,2) = sigma(k+1,2) - sigma(k,2)
      enddo
!--------------------------------------------------------------------------
!   for rrtmg
!--------------------------------------------------------------------------
      do k=1,lev+1
         si(lev+2-k)=sigma(k,1)+sigma(k,2)/1000.
      enddo
!--------------------------------------------------------------------------
!     build matrices for semi-implicit and normal mode initialization
!
       call matrix_hybrid_cwb ( cp,sigma,dsigma,ptop,ptmean,tmean,spalm &
                  , eigval,evecin,evectr,arrhyd,arsddt,pmcor,tmcor )
!
!  gaussian quadrature weights and latitudes
!
      one = 1.0
      onem= -one
      call gausl3 (my,onem,one,weight,sinl)
!
      my2= my/2
!
      do 180 j = 1, my2
      sinl(my+1-j)  = -sinl(j)
      weight(my+1-j)= weight(j)
      onocos(j)     = 1.0/(1.0-sinl(j)*sinl(j))
      onocos(my+1-j)= onocos(j)
      cosl(j)       = 1.0/sqrt(onocos(j))
      cosl(my+1-j)  = cosl(j)
  180 continue
!  define ndslfv
!
      call ndslfv_init(nx,my,ncld,cosl,weight(1))
!
!  horizontal diffusion settings
      if ( lev .eq. 80 ) then
        hdktop = 46
        hdk1   = 20
        hdk2   = 36
        hdk3   = 41
      else if ( lev .eq. 72 ) then
        hdktop =  1
        hdk1   = 20
        hdk2   = 39
        hdk3   = 43
      else if ( lev .eq. 60 ) then
        hdktop =  1
        hdk1   = 17
        hdk2   = 34
        hdk3   = 37
      endif
      factop = 1.5
!!      factop = 2.5
      
      coefu=factop/float(hdk3-hdk2)
!
!
!  define coriolis parameter for each latitude
!
      do 190 j=1,my
      cor(j)= 2.0*omega*sinl(j)
  190 continue
!
      r2d=180./pi
      if( numreduce.gt.0 ) then
        lreduce=1
        pnm_max=0.0
        do j=1,my2
          call pnmy (jtrun,sinl(j),pnm)
          do m=1,jtrun
            do n=m,jtrun
              pnm_max = max ( pnm_max, abs(pnm(n,m)) )
            enddo
          enddo
        enddo
        pnmcut = pnm_max / (10.**numreduce)
        if(myrank.eq.0)print *,' pnm_max pnmcut ',pnm_max,pnmcut
        do j=1,my2
          call pnmy (jtrun,sinl(j),pnm)
          call reducegrid(pnm,jtrun,pnmcut,j,mtrundef(j),nxdef(j),  &
                          octahedral)
          if(myrank.eq.0)print *,'j=',j,' mtrundef,nxdef=',mtrundef(j), &
                             nxdef(j),asin(sinl(j))*r2d
        enddo
      else
        lreduce=0
        do j=1,my2
          mtrundef(j)=jtrun
          nxdef(j)=nx
        enddo
      endif
!
      do j=1,my2
        jj = my + 1 - j
        mtrundef(jj)=mtrundef(j)
        nxdef(jj)   =nxdef(j)
      enddo

!for 2dMPI
      call make_list_nx  ! making nx index for reduce/non_reduce

!
      sumreduce=0.0
      do j=1,my
        sumreduce=sumreduce+nxdef(j)
      enddo
      reducefactor=sumreduce/float(nx*my)
      if(myrank.eq.0)print *,' numreduce lreduce reducefactor' &
                    ,numreduce,lreduce,reducefactor
!
!  initialize ifax and trigs for rfftmlt routine
!
      call fftfax (nx,ifax,trigs)
!
      do j=1,my
        nxj = nxdef(j)
        call fftfax (nxj,ifaxj(1,j),trigsj(1,j))
      enddo
!
      lessl_fft=.false.
      if (ibm_fft.eq.1) then
        nn = nx
        if (iand(nn,1).eq.1) go to 111   ! not an even number
        if (mod(nn,9).eq.0) then
          if (mod(nn/9,3).eq.0) go to 111   ! radix of 3**i, i>2
        endif
        if (mod(nn,5).eq.0) then
          if (mod(nn/5,5).eq.0) go to 111   ! radix of 5**i, i>1
        endif
        if (mod(nn,7).eq.0) then
          if (mod(nn/7,7).eq.0) go to 111   ! radix of 7**i, i>1
        endif
        if (mod(nn,11).eq.0) then
          if (mod(nn/11,11).eq.0) go to 111 ! radix of 11**i, i>1
        endif
        if (mod(nn,13).eq.0) go to 111   ! radix of 13
        if (mod(nn,17).eq.0) go to 111
        if (mod(nn,19).eq.0) go to 111
        if (mod(nn,23).eq.0) go to 111
        if (mod(nn,29).eq.0) go to 111
        lessl_fft=.true.
  111   continue
      endif
!
!  define associated legendre polynomials and their derivatives
!
      call lgndr (my2,jtrun,jtmax,sinl,poly,dpoly)
!
!  specify the dms read-in and write-out only for 34 keys
!
      if( myrank .eq. 0 ) then
       type_r="RORDER"//char(0)
       type_w="WORDER"//char(0)
       argument="34"//char(0)
       call dmscfg(type_r,argument,istat_r)
       call dmscfg(type_w,argument,istat_w)
       istat = abs(istat_r) + abs(istat_w)
      endif
!ch   call mpe_broadcast(istat,1,flag,mpe_integer)
      call mpe_bcast(istat,1,0,mpe_integer)
      if(istat.ne.0)then
        if(myrank.eq.0) print *,'dmscfg error'
        call mpe_finalize
        call dmsexit(-1)
      endif
!
!  open back ground dmsfile
!
      if(myrank.eq.0) then
      call dmsmsg("ALL",istat)
      call dmsopn(bckfile,"r",istat1)
!
!  open the input file.  this too will be replaced by the appropriate
!  dbms operation when available
!
      call dmsopn(ifilin,"w",istat2)

      call dmsopn(ifilout,"w",istat3)

      istat = abs(istat1) + abs(istat2) + abs(istat3)
!
! ldailyFCTsst=true, restore sst, snow depth, sea ice fraction from ncep
! data
! open ncep data dms
!
       if(ldailyFCTsst) then
          istat4=0
          call dmsopn(ifilin_sst,"r",istat4)
          istat = istat + abs(istat4)
        endif
        if(ldailyFCTicesndpt) then
          istat5=0
          call dmsopn(ifilin_ncep,"r",istat5)
          istat = istat + abs(istat5)
        endif
        if(dailyClm_option .ge. 1) then
          istat6=0
          istat7=0
          call dmsopn(ifilin_ClmANA,"r",istat6)
          if(dailyClm_option .eq. 2) then
            call dmsopn(ifilin_ClmFCT,"r",istat7)
          endif
          istat = istat + abs(istat6)+abs(istat7)
        endif

      end if
!ch   call mpe_broadcast(istat,1,flag,mpe_integer)
      call mpe_bcast(istat,1,0,mpe_integer)
      if(istat.ne.0)then
        if(myrank.eq.0) then
         print *,'dmsopn error'
         if(istat1.ne.0) print*,' BCKFILE=',bckfile,' dms open failed !'
         if(istat2.ne.0) print*,' IFILEIN=',ifilin,' dms open failed!'
         if(istat3.ne.0) print*,' IFILEOUT=',ifilout,' dms open failed!'
         if(istat4.ne.0) print*,' IFILE_SST=',ifilin_sst,' dms open failed!'
         if(istat5.ne.0) print*,' IFILE_NCEP=',ifilin_ncep,' dms open failed!'
         if(istat6.ne.0) print*,' IFILE_ClmANA=',ifilin_ClmANA,' dmsopen failed!'
         if(istat7.ne.0) print*,' IFILE_ClmFCT=',ifilin_ClmFCT,' dmsopen failed!'
        endif
        call mpe_finalize
        call dmsexit(-1)
      endif
!
!
! read file of output directives specifying desired output
! fields.
!
      open (unit=4,file=trim(ocards),form='formatted')
!
      do 80 k=1,nout
      numout= k
      read (4,800,iostat=io)  outdir(numout)
  800 format (a16)
      if (io.ne.0)  go to 85
      if (outdir(numout).eq.'nomodata')  go to 85
   80 continue
   85 numout= numout-1
      close(4)
!-----------------------------------------------------------------------
!  for WSM6
!-----------------------------------------------------------------------
      if (dolsp .and. ncld .eq. 7) then
        call wsm6init()
        ntoz=ncld
        ntcw=2
        num_p3d=5
        nclds=2
      endif

!-----------------------------------------------------------------------
!  for rrtmg scheme : rad_initialize
!-----------------------------------------------------------------------
      if (irad .eq. 2) then
       call rad_initialize (si,lev,ictm, isol, ico2, iaer, ialb,       &
       iems, ntcw, num_p3d, ntoz, iovr_sw, iovr_lw, isubc_sw, isubc_lw, &
       icliq_sw, icice_sw, icliq_lw, icice_lw, sashal, crick_proof,     &
       ccnorm, norad_precip, idate, iflip, me, myrank)

      if(myrank .eq. 0) print *,'after rad_initialize ..'
      if(myrank .eq. 0) print *,'ntoz=',ntoz,' iflip=',iflip
      if(myrank .eq. 0) print *,'me=',me,' lev=',lev,' ictm=',ictm,   &
         ' isol=', isol,' ico2=',ico2, ' iaer=',iaer, ' ialb=',ialb,    &
         ' iems=', iems,' ntcw=',ntcw,' num_p3d=', num_p3d,             &
         ' iovr_sw=',iovr_sw,' iovr_lw=', iovr_lw,                      &
         ' isubc_sw=',isubc_sw,' isubc_lw=', isubc_lw,                  &
         ' icliq_sw=',icliq_sw,' icice_sw=', icice_sw,                  &
         ' icliq_lw=',icliq_lw,' icice_lw=', icice_lw
      endif
!-----------------------------------------------------------------------
!
! check if typhoon exit
!
      tflag=0
      Wflag=0
      call getenv('GLB_TYPHINI',typhpath)
      call getenv('GLB_WTYPHINI',Wtyphpath)
!
      write(otyphfile,'(a7,i8.8,a4)')'typhoon',idtg8,'.dat'    ! for old a8 tyname
      write(oWtyphfile,'(a8,i8.8,a4)')'Wtyphoon',idtg8,'.dat'  ! for old a8 tyname
      write(ntyphfile,'(a7,i8.8,a4)')'typhoon',idtg8,'.txt'    ! for new a15 tyname
      write(nWtyphfile,'(a8,i8.8,a4)')'Wtyphoon',idtg8,'.txt'  ! for new a15 tyname
!
      call chlen(typhpath,64,ltyph)
      call chlen(Wtyphpath,64,Wltyph)
!
      typhoon=.false.
      inquire (file=typhpath(1:ltyph)//'/'//otyphfile,exist=olexist)
      inquire (file=Wtyphpath(1:Wltyph)//'/'//oWtyphfile,exist=oWlexist)
      inquire (file=typhpath(1:ltyph)//'/'//ntyphfile,exist=nlexist)
      inquire (file=Wtyphpath(1:Wltyph)//'/'//nWtyphfile,exist=nWlexist)

      if(nlexist)then
        tflag=1
        if(myrank.eq.0) print*,'typhfile= ',typhpath(1:ltyph)//'/'//ntyphfile,' exist'
        open( unit=14, file=typhpath(1:ltyph)//'/'//ntyphfile, &
             status='old', iostat=ierror )
        if( ierror .ne. 0 ) then
          if(myrank.eq.0)print*,' open typhoon file fail, model go on'
          if(myrank.eq.0)print *,typhpath(1:ltyph),ntyphfile
          go to 666
        end if
      else if(olexist)then
        tflag=2
        if(myrank.eq.0) print*,'typhfile=',typhpath(1:ltyph)//'/'//otyphfile,' exist'
        open( unit=14, file=typhpath(1:ltyph)//'/'//otyphfile, &
             status='old', iostat=ierror )
        if( ierror .ne. 0 ) then
          if(myrank.eq.0)print*,' open typhoon file fail, model go on'
          if(myrank.eq.0)print *,typhpath(1:ltyph),otyphfile
          go to 666
        end if
      endif

      if(nWlexist)then
        Wflag=1
        if(myrank.eq.0) print*,'Wtyphfile= ',Wtyphpath(1:Wltyph)//'/'//nWtyphfile,' exist'
        open( unit=15, file=Wtyphpath(1:Wltyph)//'/'//nWtyphfile, &
             status='old', iostat=ierror )
        if( ierror .ne. 0 ) then
          if(myrank.eq.0)print*,' open typhoon file fail, model go on'
          if(myrank.eq.0)print *,Wtyphpath(1:Wltyph),nWtyphfile
          go to 666
        end if
      else if(oWlexist)then
        Wflag=2
        if(myrank.eq.0) print*,'Wtyphfile=',Wtyphpath(1:Wltyph)//'/'//oWtyphfile,' exist'
        open( unit=15, file=Wtyphpath(1:Wltyph)//'/'//oWtyphfile, &
             status='old', iostat=ierror )
        if( ierror .ne. 0 ) then
          if(myrank.eq.0)print*,' open typhoon file fail, model go on'
          if(myrank.eq.0)print *,Wtyphpath(1:Wltyph),oWtyphfile
          go to 666
        end if
      endif
!
        pi=4.0*atan(1.0)
        d2r=pi/180.
        r2d=1./d2r
        do j=1,my
          tlon(1,j)=0.
!          nxj=nxdef(j)
          nxj=nx        ! findtrack do in full grid
          do i=2,nxj
            tlon(i,j)=tlon(1,j)+float(i-1)*360./nxj
          enddo
          tlat(j)=asin(sinl(j))*r2d
        enddo
!
      ntyph=0
      if(tflag .eq.1)then        ! for a15 tyname
        read(14,'(a10)')fulldtg
        read(14,'(i2)')ntyph
        if(myrank.eq.0)print*,' number of typhoons = ', ntyph
      else if(tflag .eq.2)then   ! for a8 tyname
        read(14,'(i2)')ntyph
        if(myrank.eq.0)print*,' number of typhoons = ', ntyph
      endif

      if(tflag .eq.1 .or. tflag .eq.2)then
        do n = 1, ntyph
          if(tflag .eq.1)then        ! for a15 tyname
          read(14,'(a15,1x,f4.1,1x,a1,1x,f5.1,1x,a1,1x)',iostat=irstat) &
               typhnam,clat,cns,clon,cew
          else if(tflag .eq.2)then   ! for a8 tyname
          read(14,'(a8,1x,f4.1,1x,a1,1x,f5.1,1x,a1,1x)',iostat=irstat)  &
               typhnam,clat,cns,clon,cew
          endif

          if(irstat.ne.0)then
           if(myrank.eq.0)print*,'read data error, set no typhoon exist'
           go to 666
          endif
!
          typname(n)=typhnam
          nrec(n)=0
          clattyp(n)=clat
          clontyp(n)=clon
          typhoon=.true.
          if(myrank.eq.0)print*,'typh-name= ',typhnam
          if(myrank.eq.0)print*,' at clon,clat= ',clontyp(n),clattyp(n)
          jytyp(1,n)=(clat-tlat(1))/(180./my)+1.005
!          nxj=nxdef(jytyp(1,n))
          nxj=nx        ! findtrack do in full grid
          ixtyp(1,n)=clon/(360./nxj)+1.005
          ix=ixtyp(1,n)
          jy=jytyp(1,n)
          tflat(0,1,n)=clattyp(n)
          tflon(0,1,n)=clontyp(n)
!          xshift(n)=tlon(ix,jy)-clon
!          yshift(n)=tlat(jy)-clat
          do ip=2,5
            ixtyp(ip,n)=ixtyp(1,n)
            jytyp(ip,n)=jytyp(1,n)
            tflat(0,ip,n)=clattyp(n)
            tflon(0,ip,n)=clontyp(n)
          enddo
          if(myrank.eq.0)print*,' tflat = ',tflat(0,1:5,n)
          if(myrank.eq.0)print*,' tflon = ',tflon(0,1:5,n)
          if(myrank.eq.0)print*,' position at ix,jy= ',ix,jy
          if(myrank.eq.0)print*,' position at tlat,tlon= ',tlat(jy) &
                               ,  tlon(ix,jy)
!
        enddo  ! end of do n
        close(14)
!
      endif     ! end if(tflag=1 or tflag=2)
!
      Wntyph=0
      if(Wflag .eq.1)then        ! for a15 tyname
        read(15,'(a10)')Wfulldtg
        read(15,'(i2)')Wntyph
        if(myrank.eq.0)print*,' number of Wtyphoons = ', Wntyph
      else if(Wflag .eq.2)then   ! for a8 tyname
        read(15,'(i2)')Wntyph
        if(myrank.eq.0)print*,' number of Wtyphoons = ', Wntyph
      endif

      if(Wflag .eq.1 .or. Wflag .eq.2)then
        do n = 1, Wntyph
          if(Wflag .eq.1)then        ! for a15 tyname
          read(15,'(a15,1x,f4.1,1x,a1,1x,f5.1,1x,a1,1x)',iostat=irstat) &
               typhnam,clat,cns,clon,cew
          else if(Wflag .eq.2)then   ! for a8 tyname
          read(15,'(a8,1x,f4.1,1x,a1,1x,f5.1,1x,a1,1x)',iostat=irstat)  &
               typhnam,clat,cns,clon,cew
          endif

          if(irstat.ne.0)then
           if(myrank.eq.0)print*,'read data error, set no typhoon exist'
           go to 666
          endif
!
          nc=n+ntyph
          typname(nc)=typhnam
          nrec(nc)=0
          clattyp(nc)=clat
          clontyp(nc)=360.-clon
          typhoon=.true.
          if(myrank.eq.0)print*,'typh-name= ',typhnam
          if(myrank.eq.0)print*,' at Atlantic clon,clat= ',clontyp(nc),clattyp(nc)
          jytyp(1,nc)=(clat-tlat(1))/(180./my)+1.005
!          nxj=nxdef(jytyp(1,nc))
          nxj=nx        ! findtrack do in full grid
          ixtyp(1,nc)=clontyp(nc)/(360./nxj)+1.005
          ix=ixtyp(1,nc)
          jy=jytyp(1,nc)
          tflat(0,1,nc)=clattyp(nc)
          tflon(0,1,nc)=clontyp(nc)
          do ip=2,5
            ixtyp(ip,nc)=ixtyp(1,nc)
            jytyp(ip,nc)=jytyp(1,nc)
            tflat(0,ip,nc)=clattyp(nc)
            tflon(0,ip,nc)=clontyp(nc)
          enddo
          if(myrank.eq.0)print*,' tflat = ',tflat(0,1:5,nc)
          if(myrank.eq.0)print*,' tflon = ',tflon(0,1:5,nc)
          if(myrank.eq.0)print*,' position at ix,jy= ',ix,jy
          if(myrank.eq.0)print*,' position at tlat,tlon= ',tlat(jy) &
                               ,  tlon(ix,jy)
!
        enddo  ! end of do n
        ntyph=nc
        close(15)
!
      endif     ! end if(Wflag=1 or Wflag=2)
!
 666  continue

!for 2dMPI >>

      allocate (eps4L(jtp),                 &
                plnowL(jtp,2),              &
                ploldL(jtp,2),              &
                pltenL(jtp,2), stat=ierror) 

      if (ierror/= 0) then
          write(6,*) 'cons : allocate fail 1 '
          stop
      end if

      call mpe2d_reshape_eps4(eps4, eps4L)

!   checking owner of idg and jdg, so no need for rcup/rlsp/tg to call mpe_unify
!   in diabat to print out vertical profile at selected point (idg,jdg)

     idg_jdg_owner=.false.
     do jj=1,jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
     do ii=1,nxj
        i=map2to1(ii,j)
        if((i .eq. idg) .and. (j .eq. jdg))then
          idg_listnum=ii
          jdg_listnum=jj
          idg_jdg_owner=.true.
          goto 990
        endif
     enddo
     enddo
990  continue

!for 2dMPI <<

      return
      end
