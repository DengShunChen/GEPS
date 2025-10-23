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
! #ifdef USE_CUDA
!       use mod_ndslfv_monoadv_gpu, only: allocate_ndslfv_array_gpu
! #endif
!-----------------------------------------------------------------------
      use radn
      use noah
!-----------------------------------------------------------------------
      ! for ECHAM4 Tiedtke cumulus scheme
      USE mo_cumulus_flux, only : cuparam
      USE mo_constants,    only : inicon
      USE mo_convect_tables, only : set_lookup_tables
      ! for stochastic_physics
      use mod_stochastic_physics, only : ncep_seeds, &
                   sppt, sppt_seed, sppt_decort, sppt_lscale, &
                   sppt_sigtop1, sppt_sigtop2, sppt_sigbot1, sppt_sigbot2, &
                   sppt_sfclimit, sppt_logit, &
                   shum, shum_seed, shum_decort, shum_lscale, &
                   shum_sigefold, &
                   skeb, skeb_seed, skeb_decort, skeb_lscale, &
                   skeb_sigtop1, skeb_sigtop2, skeb_sigbot1, skeb_sigbot2, &
                   skeb_vdof,skebnorm, skebfilt, &
                   ssst, ssst_seed, ssst_decort, ssst_lscale, &
                   init_stochastic_physics
      use mod_grb2_param, only:grbmem,grbnumm

      implicit  none

      integer i,j,k,n,jj,nn,istat,istat1,istat2,istat3,irstat,ii,Wntyph,nc,Wltyph
      integer ix,jy,ip,istat_r,istat_w,ierror,ltyph,io
      integer ifromtau,itotau,itau,itaui,l,nxj,m,mf,my2,mf1,mf2,m2,m1len,m2len
      integer jlistnum_fj, llistnum_fj, j_str, m_str, ind
      integer tflag,Wflag,ntrac_req

      real    pi,rm,rl,rlm,one,onem,r2d, pnm_max,pnmcut,sumreduce
      real    reducefactor,d2r,cew,clon,cns,clat,prslp

      integer*8 idtg8

      real(kind=RTYPE) pnm(jtrun+1,jtrun+1)

!
      namelist /modlst/ ksgeo,ptmean,dt,taui,taue                       &
                      , tauo,frad,ktpbl,ktshl,ktcup,njump,evaprh,lsimpl &
                      , yesdia,dopbl,docup,dorad,dolsp,dograv           &
                      , doshl,dodry,donnmi,idg,jdg,ldiag,nnmiit,nnmivm  &
                      , cutfreq,hdiff,itypbl,cstar,taup,hfilt           &
                      , ptmeans,update,taureg,doincr,numreduce          &
                      , nmcup,nmpbl,nmland,nmshl,cgw,ggdef,gmdef        &
                      , nmgwor,tofd,nmgwcv,mtnvar,docgrav               &
                      , ictm,isol,ico2,iaer,ialb,irad,iems,ntcw         &
                      , ntoz,iovr_sw,iovr_lw,isubc_sw,isubc_lw          &
                      , sashal,crick_proof,ccnorm,norad_precip,me,doo3l &
                      , ioutsigr,domfc,out_green,isot,ivegsrc           &
                      , otgreen,out_hp,dosppt,dospptout,doshum,doshumout &
                      , dossst &
                      , doskeb,doskeb_dc,doskebout,ndsladvh2,hord       &
                      , ldailyFCTsst,ldailyFCTicesndpt,lFCTweight       &
                      , dailyClm_option,lopgsst,do_sit,fsit,pdfcloud,updatetg       &
! output data for RSM (Also, RSM compiling flag is necessary)
                      , outrsm,rsmoutinv,rlon1,rlon2,rlat1,rlat2,rgrdsz,rsmsfcmgrhr &
!
                      , cmbk,cgwd,nmmiph,spl1,spl2                      &
                      , weightSIT,dSITdt_intv,mwhd,doclx,doslavepp      &
                      , outdms,outgrb2,alpha,two_loop,ttl,tfilt,factop  &
                      , mass_dp,dpprt,outfv3,itter,vd,dorst,naero,mom4ice
!
      real    si(lev+1)
      logical flag
      character*10 fulldtg,Wfulldtg
      character*255 filist
      character cdtg*12
      character*255 pathname,logicname,truefile
      character*64 type_r,type_w,argument
      integer istat4,istat5,istat6,istat7,istat8
      integer :: ios
      data pathname/'NWPETCGLB'/
      data logicname/'filist'/


      namelist /filst/ ifilin,cwbout,bckfile,namlsts &
                     , ifilout,crdate,ocards,phyout,cntrl &
                     , ifilin_ncep, ifilin_sst, ifilin_nc &
                     , ifilin_ClmANA,ifilin_ClmFCT,ifilout_grb &
                     , ifilin_aero

      namelist /typ/ write_tau, write_mem, trk_intv, min_trk_pres

      namelist /stochy_physics/ ncep_seeds, &
                   use_zmtnblck, &
                   sppt, sppt_seed, sppt_decort, sppt_lscale, &
                   sppt_sigtop1, sppt_sigtop2, sppt_sigbot1, sppt_sigbot2, &
                   sppt_sfclimit, sppt_logit, &
                   shum, shum_seed, shum_decort, shum_lscale, &
                   shum_sigefold, &
                   skeb, skeb_seed, skeb_decort, skeb_lscale, &
                   skeb_sigtop1, skeb_sigtop2, skeb_sigbot1, skeb_sigbot2, &
                   skeb_vdof, skebnorm, skebfilt, &
                   ssst, ssst_seed, ssst_decort, ssst_lscale
      namelist /grb_conf/ grbmem,grbnumm
      namelist /gce_3ice/ SL_sedi, sat_predict, new_saturation, &
                          use_cpm, use_declination

! for ECHAM4 Tiedtke cumulus scheme
      call cuparam
      call inicon
      call set_lookup_tables

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
      ! read stochastic_physics
      read (1,stochy_physics,end=122)
  122 continue
      read (1,gce_3ice,end=124)
  124 continue
      read (1,grb_conf,end=123)
  123 continue
      close(1)

!
      open (unit=1,file=trim(namlsts),form='formatted')

      if(do_sit) then
        read (1,sit_nml,end=130)
        if(myrank .eq. 0) then
          print *, 'chlee debug...'
          print sit_nml
        endif
        if(.NOT. ldailyFCTsst)then
          print *,'do_sit=.true., auto set: ldailyFCTsst=.true.' &
              ,',dailyClm_option=1, check ifilin_sst & ifilin_ClmANA'&
              ,'is in filelist.'
          ldailyFCTsst=.true.
          dailyClm_option=1
        endif
  130 continue
      endif
      close(1)

!
      open (unit=2,file=trim(crdate),form='formatted')
!
! read in idtg*12
      read(2,'(i8.8)')idtg8
      close(2)
! transfer idtg8 to idtg*12
      if(idtg8.gt.60000000)then
        idtg = 190000000000_8 + idtg8*100
      else
        idtg = 200000000000_8 + idtg8*100
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
      call leapyear(idate(1))
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
        hours = hours+taui-(dt/3600.)  ! for restart
        julian= julian+hours/24.0+0.001
        hours = mod(hours,24.)
      endif
!
      if(myrank .eq. 0) print modlst
      if(myrank .eq. 0) print typ
      if(myrank .eq. 0) print stochy_physics
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
      coslr(j)      = 1./onocos(j)
      coslr(my+1-j) = coslr(j)
      cosl(j)       = 1.0/sqrt(onocos(j))
      cosl(my+1-j)  = cosl(j)
  180 continue
!  define ndslfv
!
      call ndslfv_init(nx,my,ncld,cosl,weight(1))
!
!  horizontal diffusion settings for sponge layer
      do k = 1, lev
        prslp=sigma(k,2)+sigma(k,1)*1000.+ptop
        if ( prslp .le. spl1  ) hdk1=k
        if ( prslp .le. spl2  ) hdk2(1)=k
        if ( prslp .le.  50.  ) hdk2(2)=k
        if ( prslp .le. 200.  ) hdk2(3)=k
      enddo
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
#ifdef VERBOSE
          if(myrank.eq.0)print *,'j=',j,' mtrundef,nxdef=',mtrundef(j), &
                             nxdef(j),asin(sinl(j))*r2d
#endif
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
! #ifdef USE_CUDA
!       call allocate_ndslfv_array_gpu
! #endif
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
#ifdef USE_CUDA
      tcolt_jlist = 0
      do m = 1, mlistnum
         mf = mlist(m)
         do j = 1, my/2
            if (mf .le. mtrundef(j)) then
               tcolt_jlist(1, m) = tcolt_jlist(1, m) + 1
            end if
         end do
         do j = 1, my/2
            if (mf .le. mtrundef(j)) then
               tcolt_jlist(2, m) = j
               exit
            end if
         end do
      end do

      poly_mlist(1) = 1
      do m = 1, mlistnum-1
         mf = mlist(m)
         llistnum_fj = jtrun - mf + 1
         poly_mlist(m+1) = llistnum_fj*my + poly_mlist(m)
      end do

      polyf = 0.
      dpolyf = 0.
      do m = 1, mlistnum
         m_str = poly_mlist(m)
         jlistnum_fj = tcolt_jlist(1, m)
         j_str = tcolt_jlist(2, m)

         mf = mlist(m)
         llistnum_fj = jtrun - mf + 1
         do j = 1, jlistnum_fj
            do l = mf, jtrun
               ind = (l - mf + 1) + (j - 1)*llistnum_fj + m_str - 1
               jj = j_str + j -1
               polyf(ind) = poly(l, jj, m)
               dpolyf(ind) = dpoly(l, jj, m)

               ind = ind + jlistnum_fj*llistnum_fj
               jj = my2 - j + 1
               polyf(ind) = (-1)**(l - mf)*poly(l, jj, m)
               dpolyf(ind) = (-1)**(l - mf + 1)*dpoly(l, jj, m)
            end do

         end do
      end do
#endif
!
!  specify the dms read-in and write-out only for 34 keys
!
      if( myrank .eq. 0 ) then
       type_r="RORDER"//char(0)
       type_w="WORDER"//char(0)
#ifdef I38K
       argument="38"//char(0)
#else
       argument="34"//char(0)
#endif
       call dmscfg(type_r,argument,istat_r)
#ifdef O38K
       argument="38"//char(0)
#else
       argument="34"//char(0)
#endif
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
      endif
!
!  open the input file.  this too will be replaced by the appropriate
!  dbms operation when available
!
!      if(col_rank .eq. 0) call dmsopn(ifilin,"w",istat2)
      if(col_rank .eq. 0) call dmsopn(ifilin,"r",istat2)
!
      if(myrank .lt. lev) call dmsopn(ifilout,"w",istat3)
!
      if(myrank.eq.0) then
      istat = abs(istat1) + abs(istat2) + abs(istat3)
!
! ldailyFCTsst=true, restore sst, snow depth, sea ice fraction from ncep
! data
! open ncep data dms
!
       istat4=0; istat5=0; istat6=0; istat7=0; istat8=0
        if(ldailyFCTsst) then
          call dmsopn(ifilin_sst,"r",istat4)
          istat = istat + abs(istat4)
        endif
        if(ldailyFCTicesndpt) then
          call dmsopn(ifilin_ncep,"r",istat5)
          istat = istat + abs(istat5)
        endif
        if(dailyClm_option .ge. 1) then
          call dmsopn(ifilin_ClmANA,"r",istat6)
          if(dailyClm_option .eq. 2) then
            call dmsopn(ifilin_ClmFCT,"r",istat7)
          endif
          istat = istat + abs(istat6)+abs(istat7)
        endif
#ifdef Readaeroclx
        call dmsopn(ifilin_aero,"r",istat8)
        istat = istat + abs(istat8)
#endif

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
#ifdef Readaeroclx
         if(istat8.ne.0) print*,' IFILE_AERO=',ifilin_aero,' dmsopen failed!'
#endif
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
!  for cloud microphysics initialization
!-----------------------------------------------------------------------
      if ( nmmiph .eq. 11 .or. nmmiph .eq. 12 .or. nmmiph .eq. 13 ) then
        ntrac_req = 6   ! only six species of hydrometeors for GFDL MP
      elseif ( nmmiph .eq. 18 ) then
        ntrac_req = 6   ! only six species of hydrometeors for 2M Thompson MP
      elseif ( nmmiph .eq. 8 ) then
        ntrac_req = 6   ! only six species of hydrometeors for Thompson MP
      elseif ( nmmiph .eq. 15 ) then
        ntrac_req = 6   ! only six species of hydrometeors for Goddard 3ICE MP
      elseif ( nmmiph .eq. 16 ) then
        ntrac_req = 6   ! only six species of hydrometeors for Goddard 4ICE MP
      else
        ntrac_req = nmmiph
      endif
      if ( ntoz .gt. 0 ) then
        ntrac_req = ntrac_req + 1
        ntoz = ncld
      endif
      if ( dolsp ) then
        if ( ncld .lt. ntrac_req ) then
           if ( myrank .ge. 0 ) print *,'not enough number of tracers'
           call mpe_finalize
           call dmsexit(-1)
        endif
!
        if ( nmmiph.eq.6 .or.                                           &
             nmmiph.eq.8 .or. nmmiph.eq.18 .or.                         &
             nmmiph.eq.11 .or. nmmiph.eq.12 .or. nmmiph.eq.13 .or.      &
             nmmiph.eq.15 .or. nmmiph.eq.16 )                           &
          call mp_init(nmmiph,myrank)
!
      endif

!-----------------------------------------------------------------------
!  for rrtmg scheme : rad_initialize
!-----------------------------------------------------------------------
      if (irad .eq. 2) then
       call rad_initialize (si,lev,ictm, isol, ico2, iaer, ialb,        &
       iems, ntcw, nmmiph, ntoz, iovr_sw, iovr_lw, isubc_sw, isubc_lw,  &
       icliq_sw, icice_sw, icliq_lw, icice_lw, sashal, crick_proof,     &
       ccnorm, norad_precip, idate, iflip, me, myrank)

      if(myrank .eq. 0) print *,'after rad_initialize ..'
      if(myrank .eq. 0) print *,'ntoz=',ntoz,' iflip=',iflip
      if(myrank .eq. 0) print *,'me=',me,' lev=',lev,' ictm=',ictm,   &
         ' isol=', isol,' ico2=',ico2, ' iaer=',iaer, ' ialb=',ialb,    &
         ' iems=', iems,' ntcw=',ntcw,' nmmiph=', nmmiph,               &
         ' iovr_sw=',iovr_sw,' iovr_lw=', iovr_lw,                      &
         ' isubc_sw=',isubc_sw,' isubc_lw=', isubc_lw,                  &
         ' icliq_sw=',icliq_sw,' icice_sw=', icice_sw,                  &
         ' icliq_lw=',icliq_lw,' icice_lw=', icice_lw
      endif
!-----------------------------------------------------------------------
!  for stochastic_physics initialization
!-----------------------------------------------------------------------
      call init_stochastic_physics(dt)
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

        tlon(1)=0.
!          nxj=nxdef(j)
!        nxj=nx        ! findtrack do in full grid
        do i=2,nx
          tlon(i)=tlon(1)+float(i-1)*360./nx
        enddo
        do j=1,my
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
!          xshift(n)=tlon(ix)-clon
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
                               ,  tlon(ix)
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
                               ,  tlon(ix)
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
