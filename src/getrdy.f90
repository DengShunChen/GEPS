      subroutine getrdy
!
!***********************************************************************
!  subroutine to prepare initial conditins for model start at tau=0
!  reads in grid point fields and tranform to spectral coefficients.
!
! **** input ****
!
!  restrt:  logical variable. when true = initial tau .gt. zero
!  restart flag.
!
!  modify to f90 by C-H Lee and sort by River Chen
!
!***********************************************************************
!
      use param
      use mpe
      use rank
      use index
      use const
      use grid
      use spec
      use fftcom
      use phygrid
      use mod_typhoon
      use noah
      use namelist_soilveg
!-----------------------------------------------------------------------
      use ozne_def
      use radn
      use albn
!-----------------------------------------------------------------------
      use mod_sitgrid
      use mod_sit_vdiff,     only:sit_vdiff_init,sit_vdiff,SICEDFN &
                               ,maskid,ctfreez,cal_ratioBlending
      use mod_sit_control,   only:xmissing,sit_nml,lgodas,ldailysst &
                               ,locaf,locaf0,lwoa0,lsitstart,sit_domain_w &
                               ,sit_domain_e,sit_domain_s,sit_domain_n &
                               ,sit_domain_extgrd,lpre6hr_sit,bathydepth&
                               ,lsftobswt
      use mod_eos_ocean,     only: tmelts,api
      use mod_sst,           only:read_woa0,read_godas,read_dailygodas &
                                ,read_dailyFCT                         &
                                ,ifilin_ocaf,read_ocaf,read_ocaf0      &
                                ,wtfn12,wsfn12,time_weights,mask1st
      USE mo_netcdf,         ONLY:lkvl,set_ocndepth
! for Thompson MP
      use physcons,          only: con_rd,con_eps
      use module_mp_thompson_make_number_concentrations,                &
                             only: make_IceNumber, make_RainNumber
!-----------------------------------------------------------------------
      use mod_grb2_param , only : grbid ,grbfile ,opn_grb2 ,cls_grb2,grbnxmy

      implicit   none

!  local working array
!
      real      sst(nxp,my_max),ww1(nx,my),ww2(nxp,my_max),     &
                wk1(nxp,lev,my_max),pklev(nxp,my_max)
      real(kind=RTYPE) cc(nx+2,levp,1,my_max),dummy,ww3(nx,my_max),    &
                       ww4(nx,my)
!byl                wss3(levp,2,3,jtrun,jtmax),cc3(nx+2,levp,3,my_max)

      character lrec*26,rfile*55,ctau*6,topostd*4,topohgt*4,key*34
!helio>
      character f71*50
!helio<
#ifdef RSM
      character*12 dtgrsm
      integer idtgrsm
#endif
!
! restart  : read(7) work array
!
!!      real, dimension(:), allocatable :: work_io
!
!
! restart  : read(10) work array
!
      real, dimension(:,:,:), allocatable :: tm1,tm2,tm3,tm4
!
! for ncld>=2
      real, dimension(:,:,:), allocatable :: tmc1,tmc2
      real, dimension(:,:,:), allocatable :: tmc3,tmc4
      real, dimension(:,:,:), allocatable :: tmc5,tmc6
! add for sfcuvt
      real    tx(nxp), qx(nxp),ux(nxp),vx(nxp),  &
              qs(nxp),tsx(nxp),hs(nxp),ps(nxp),  &
              qsfc(nxp),tgp(nxp)
      logical flg, snow
!
! add for soil
      real, dimension(:,:,:), allocatable :: temp1,temp2,temp3
!     integer   ls(nxp,my_max),icex(nx,my)
      integer   ls(nxp,my_max)
      logical   ncepsnow,ncepice
!soil
      logical wrestrt
      data wrestrt/.false./
!
      integer lmax,nxmy,mlmax2,i,j,jj,k,m,mf,n,nxj,ios,lcwb,lphy,itaui, &
              lncrec,istat,itaup,isnow,njump1,njump2,njump3,lvlw,lvlw1, &
              lvlw2,lvlw3,ii,icwarn
      real    fact,xxaa,taux,q1,dsigp,pi,xx,wet
!xb110>
!      real    flash(nxp,my_max)
!xb110<

      real  t_surf
      real  sitlat(nxp)
      real  sitlon(nxp,my_max)
      real, parameter:: specified_ice_thickness  = 2.0
      real lontest(nxp,my_max)
      integer nxjpart      
! for io quilting
      character:: keydoit*34
! for Thompson MP
      real  tem,rho,ttr,ttv

      lmax=26
!
      nxmy   = nx*my
      mlmax2 = mlmax * 2
! noah
      call set_soilveg(isot,ivegsrc)
! osu
      call landpack(maxsmc,dfkt,xktk,dfk)

! ------------------------------------------------------------
!   read ozone prognostic parameters
! ------------------------------------------------------------
      if (doo3l)then
      if (ntoz .gt. 0) then
      call read_ozplin(myrank)
      if (myrank.eq.0) print *, 'after read_ozplin'
      if (myrank.eq.0) then
            print *,' pl_coeff=',pl_coeff
            print *,' latsozp=',latsozp,' levozp=',levozp,     &
                    ' timeoz=',timeoz
!           print * ,'ozplin(1,1,1,1)=',ozplin(1,1,1,1)
      endif
      endif
      endif
!------------------------------------------------------------
!helio>
       if ( nco .gt. 999 ) then
        write(f71,105) nco,nx,my
       else
        if ( nx .gt. 999 .and. my .gt. 999 ) write(f71,106) nco,nx,my
        if ( nx .gt. 999 .and. my .le. 999 ) write(f71,107) nco,nx,my
        if ( nx .le. 999 .and. my .le. 999 ) write(f71,108) nco,nx,my
       endif
 105  format('global_idw.t',i4.4,'.',i4.4,'.',i4.4,'.dat')
 106  format('global_idw.t',i3.3,'.',i4.4,'.',i4.4,'.dat')
 107  format('global_idw.t',i3.3,'.',i4.4,'.',i3.3,'.dat')
 108  format('global_idw.t',i3.3,'.',i3.3,'.',i3.3,'.dat')

      open(71,file=f71,form='unformatted',access='direct',recl=8*nx*my)

      do k=1,8
        read(71,rec=k) ww1
        do jj = 1, jlistnum
          j=jlist1(jj)
          outp(:,jj,k) = ww1(:,j)
        enddo
      enddo

      close(71)

!helio<
!------------------------------------------------------------

      if (restrt) then
!
! restart gfcst: read restart file and pdiff file
!
      call chlen (cwbout,48,lcwb)
      call chlen (phyout,48,lphy)
      itaui = taui + 0.001
      write (ctau,800) itaui
  800 format('tau',i4.4)
!
      rfile = cwbout(1:lcwb)//ctau
!
      if(myrank .eq. 0) then
      print*, ' restart at tau, lcwb, rfile = ', taui, lcwb, rfile
      endif

      flg=.false.
      if(myrank .eq. 0) then
       flg=.true.
       open (unit=7,file=rfile,form='unformatted',status='old' &
            ,iostat=ios)
      endif
!
!ch   call mpe_broadcast(ios,1,flg,mpe_integer)
      call mpe_bcast(ios,1,0,mpe_integer)
      if(ios .ne. 0) goto 820
!
!!      if(myrank .eq. 0) then
!!      read (7)  work_io
!!      close (7)
!!      endif
!
!!      call scatter_spec(work_io,vornow,divnow,temnow,qnow,plnow,   &
!!         vorold,divold,temold,qold,plold,dsqgeo,spgeo,trefs,       &
!!         lev,ncld,jtrun,jtmax,my,nsize)
!
      rfile = phyout(1:lphy)//ctau
!
      if(myrank .eq. 0) then
      print*, ' restart at tau, lphy, rfile = ', taui, lphy, rfile
      endif
!
      open (unit=10,file=rfile,form='unformatted',status='old',err=820)

      allocate (tm1(nx,lev,my))
      allocate (tm2(nx,lev,my))
      allocate (tm3(nx,lev,my))
      allocate (tm4(nx,lev,my))
      allocate (tmc1(nx,lev,my))
      allocate (tmc2(nx,lev,my))
      allocate (tmc3(nx,lev,my))
      allocate (tmc4(nx,lev,my))
      allocate (tmc5(nx,lev,my))
      allocate (tmc6(nx,lev,my))
      allocate (temp1(nx,lev,my))
      allocate (temp2(nx,lev,my))
!
!      read(10) snr,gwr,tg,tm1,tm2,ss,rs,tm3,tm4,ustar,tstar,qstar   &
!              , hflux,qflux,raincu,rainlp,totalp,curate,plcl,cumtop &
!--
!reduceg : base on sgeo to induce spgeo, different from normal version
!--
!              , tgclim,gwet,z0,alb,land,ice,ocean,gwclim,acld,sgeo  &
!              , temp1,temp2,canopy,sigmaf,istyp,ivegtyp,rld         &
!              , tmc1,tmc2,tmc3,tmc4,fpsp,fpsp1
!      close (10)
!
      do jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
       do k = 1, lev
        do i = 1, nxj
!        e(i,k,jj)     = tm1(i,k,j)
!        eps(i,k,jj)   = tm2(i,k,j)
!        o3l(i,k,jj)   = tm3(i,k,j)
!        dtrad(i,k,jj) = tm4(i,k,j)
!        ftp(i,k,jj)   = tmc1(i,k,j)
!        fqp(i,k,jj)   = tmc2(i,k,j)
!        ftp1(i,k,jj)   = tmc3(i,k,j)
!        fqp1(i,k,jj)   = tmc4(i,k,j)
        enddo
       enddo
      enddo
!soil
      do jj=1,jlistnum
       j=jlist1(jj)
       nxj=nxdef_2d(j)
       do k = 1, km_soil
        do i = 1, nxj
!        smc(i,k,jj)=temp1(i,k,j)
!        stc(i,k,jj)=temp2(i,k,j)
        enddo
       enddo
      enddo
!
      deallocate (tm1)
      deallocate (tm2)
      deallocate (tm3)
      deallocate (tm4)
      deallocate (tmc1)
      deallocate (tmc2)
      deallocate (tmc3)
      deallocate (tmc4)
      deallocate (tmc5)
      deallocate (tmc6)
      deallocate (temp1)
      deallocate (temp2)

      lncrec=nx*my
      call syslbl('x00dif',idtg,0,ggdef,lrec)
      call dmsread(nx,my,lrec,lncrec,'H',ifilout,ww1,istat)
      call unify_reducepick(nx,my,my_max,ww1,pdiff)

      call syslbl('h00100',idtg,0,ggdef,lrec)
      call dmsread(nx,my,lrec,lncrec,'H',ifilin,ww1,istat)
      call unify_reducepick(nx,my,my_max,ww1,t1000)

      call syslbl('x00tsv',idtg,0,ggdef,lrec)
      call dmsread(nx,my,lrec,lncrec,'H',ifilout,ww1,istat)
      call unify_reducepick(nx,my,my_max,ww1,tsave)
!
!  read dmsdata for standard deviation of terrain field
!
!dms    istdno=99
!dms    istdno=0
!c      topostd='gbkf'   ! responding to istdno=99
        topostd='gbk0'   ! responding to istdno=0
        write(lrec,'("s00062",a4,a4,12x)')topostd,ggdef
        call dmsread(nx,my,lrec,nxmy,'H',bckfile,ww1,istat)
!
        do jj=1,jlistnum
          j=jlist1(jj)
          ii=nxjstart(j)
          nxj=nxdef_2d(j)
        if( lreduce.eq.1 ) call reducepickr (ww1(1,j),nxdef(j),nx,1)
          do i=1,nxj
            std(i,jj)=ww1(ii,j)*ww1(ii,j)
            if(ww1(ii,j).le.0. .or. ocean(i,jj)) std(i,jj)=0.
            ii=ii+1
          enddo
        enddo
!
      go to 900
  820 continue
      if ( rfile(1:lcwb) .eq. cwbout(1:lcwb) )  close (7)
      if ( rfile(1:lphy) .eq. phyout(1:lphy) )  close (10)
!
      if(myrank .eq. 0) then
        print*, ' restart failed at tau, rfile = ', taui, rfile
        print*, ' the gfcst run will start from tau = 0 again'
      endif
!
      taui = 0.0
      restrt = .false.
      donnmi = .true.
  900 continue
!
      endif     ! end of (restrt=true)

      if ( .not. restrt )  then
!
! zero out array for ncld>=2
!
        fpsp=0.
        fpsp1=0.
        ftp=0.
        fqp=0.
        ftp1=0.
        fqp1=0.
        itaui=0
!
! new start gfcst: read climate data, initialize parameters
!
        call readclx( nx,my,my_max,julian,land,ocean,ice,tgclim,gwclim  &
                   ,z0,alb,sst,bckfile,sigmaf,istyp,ivegtyp,ls   &
                   ,shdmax,shdmin,slopetyp,snoalb,ggdef,isot,ivegsrc )
!
!  read sst analysis data
!
        call syslbl('w00100',idtg,0,ggdef,lrec)
        call dmsread(nx,my,lrec,nxmy,'H',ifilin,ww1,istat)
        call unify_reducepick(nx,my,my_max,ww1,sst)
! ------------------------------------------------------------
!   read new albedo
!-------------------------------------------------------------
      if (irad .eq. 2) then
        call readalb(bckfile,nx,my,my_max,julian,ggdef,         &
                   alvsf,alvwf,alnsf,alnwf,facsf,facwf)
        if (myrank.eq.0) print *, 'irad=2, readalb ok!!'
        do jj=1,jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do i=1,nxj
            alvsf(i,jj)=alvsf(i,jj)*0.01
            alvwf(i,jj)=alvwf(i,jj)*0.01
            alnsf(i,jj)=alnsf(i,jj)*0.01
            alnwf(i,jj)=alnwf(i,jj)*0.01
            facsf(i,jj)=facsf(i,jj)*0.01
            facwf(i,jj)=facwf(i,jj)*0.01
          enddo
        enddo
      endif
! ------------------------------------------------------------
      if(.not.cstar)then
!
        if(taup .gt. 72.)then
          if(myrank .eq. 0) then
          print *,'*** warning update taup greater than 72. ***'
          endif
          call mpe_finalize
          call dmsexit(-1)
        endif
!
        itaup=int(taup+0.0001)
!
!  dtgfix12 needs idtg of long integer for dms 34 keys
!
        call dtgfix12(idtg,idtg2,-itaup)
        call rdpbl(nx,my,my_max,snr,gwr,tg,zice,ifilin,idtg2,itaup,ggdef)
!
!  if true, read ncep snow depth and sea-ice analysis data              
!                                                         
        ncepsnow =.true.            
        ncepice  =.true.
!                     
!  update snow depth with ncep's snow analysis at 00Z
!  otherwise, keep model's cycle
!
        isnow = idtg - (idtg/10000)*10000
!                           
        if( ncepsnow  .and. isnow.eq.0 )then
          call syslbl('b00650',idtg,0,ggdef,lrec)
          call dmsread(nx,my,lrec,nxmy,'H',ifilin,ww1,istat)
          call unify_reducepick(nx,my,my_max,ww1,snr)
          if( myrank .eq. 0 ) print*, &
            "update snow depth with ncep's snow analysis, at dtg=",idtg
        endif
!         
! should prepare a data set
! temporary setting for test run
        do jj=1,jlistnum
          j=jlist1(jj)
          ii=nxjstart(j)
          nxj=nxdef_2d(j)
          do i=1,nxj
! slopetyp
! read from background
!           slopetyp(i,jj) = 1
!
! sndepth(mm), snow real depth
! set to be 3 times of snr(water equivlent snow depth(mm))
!           sndepth(i,jj)=snr(i,jj)*3.    !noah1
            sndepth(i,jj)=snr(i,jj)*8.    !noah1x
!           sndepth(i,jj)=snr(i,jj)*10.
!
! sncover,snow fraction of a grid
! just set an inital value and will be recalculated in sflx.f
            sncover(i,jj)=min(1., snr(i,jj)/400.)
          enddo
        enddo
!
        if( ncepice )then
!
!          call syslbl('w00090',idtg,0,ggdef,lrec)
!          call dmsreadi(nx,my,lrec,nxmy,'I',ifilin,icex,istat)
!          if( lreduce.eq.1 ) call reducepicki (icex,nxdef,nx,my)
!
          call syslbl('w00091',idtg,0,ggdef,lrec)
          call dmsread(nx,my,lrec,nxmy,'H',ifilin,ww1,istat)
          call unify_reducepick(nx,my,my_max,ww1,cice)
!
          if( myrank .eq. 0 ) then
             print*,"get ncep's sea ice analysis, at dtg=",idtg
             call syslbl('w00092',idtg,0,ggdef,lrec)
             write(key,'(a26,a1,i7.7)') lrec,'H',nxmy
             call dmschkr (ifilin,key//char(0),istat)
          endif
          call mpe_bcast(istat,1,0,mpe_integer)
!
          if ( istat .eq. 0 ) then
            call dmsread(nx,my,lrec,nxmy,'H',ifilin,ww1,istat)
            call unify_reducepick(nx,my,my_max,ww1,zice)
            if( myrank .eq. 0 ) &
               print*,"get sea ice thickness from ncep analysis,",     &
               " at dtg=",idtg
          else
            if( myrank .eq. 0 ) &
               print*,"get sea ice thickness from model 6hr forecast,",&
               " initial at dtg=",idtg2
          endif
!
! reset albedo and tgclim at seaice grids                 
!
          icwarn=0
          do jj=1,jlistnum
            j=jlist1(jj)
            nxj=nxdef_2d(j)
            do i = 1, nxj
              ice(i,jj) = .false.               
              if( ls(i,jj).eq.0 )then
                alb(i,jj)   = 0.09
                ocean(i,jj) = .true.
              endif
              if( ls(i,jj).eq.0 .and. cice(i,jj).ge.0.5 )then
                ice(i,jj)    = .true.
                ocean(i,jj)  = .false.                     
                alb(i,jj)    = 0.55
                tgclim(i,jj) = 271.2
! for noah
                z0(i,jj)=0.00001

! initial sea ice temperature is set as tg(grid average temperature from first guess)
!
                cice(i,jj)    = max(0.5,cice(i,jj))
                xtice(i,jj)   = tg(i,jj)
! zice read from first guess or set to be 1 m over south hemesphere.
!                                     and 3 m over north hemesphere
                if( zice(i,jj) .lt. 1.*cice(i,jj) ) then 
                  zice(i,jj) = max(zice(i,jj),1.*cice(i,jj))
                  icwarn = icwarn + 1 
                endif
              endif            
            enddo             
          enddo
!
          call mpe_global_sum(icwarn,1,mpe_integer)
          if ( myrank .eq. 0 .and. icwarn .gt. 0 ) then
          print*,achar(27)//"[1;31m==================   Warnig!!!   ==================="//achar(27)//'[1;m'
          print*,achar(27)//"[1;31m=  ice thickness not consistent with sea ice mask  ="//achar(27)//'[1;m'
          print*,achar(27)//"[1;31m=  set the thickness to 1 meter for first guess    ="//achar(27)//'[1;m'
          print*,achar(27)//"[1;31m===================================================="//achar(27)//'[1;m'
          endif 
! 
        endif     !end of (ncepice)
!soil             
        call rdsoil(nx,my,my_max,km_soil,smc,stc,slc,canopy &
                  ,ifilin,idtg2,itaup,ggdef,gmdef)

!
!put upper and lower bound to soil moisture
!slc(liquid soilmoisture) should be smaller than smc, then its upper
!limit is smc
! 201704 by phon
!
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do i = 1, nxj
            slc(i,1,jj)=min( smc(i,1,jj),slc(i,1,jj) )
            slc(i,2,jj)=min( smc(i,2,jj),slc(i,2,jj) )
            slc(i,3,jj)=min( smc(i,3,jj),slc(i,3,jj) )
            slc(i,4,jj)=min( smc(i,4,jj),slc(i,4,jj) )
          enddo
        enddo
!                                                      
! blending smc(i,2,jj) with climate value
!                                       
        fact=exp(-1./60.*itaup/24.)      
        do jj = 1, jlistnum      
          j=jlist1(jj)            
          nxj=nxdef_2d(j)
          do i = 1, nxj
            xxaa=slc(i,4,jj)/smc(i,4,jj)
            smc(i,4,jj)=smc(i,4,jj)*fact + gwclim(i,jj)*(1.-fact)
            slc(i,4,jj)=smc(i,4,jj)*xxaa
          enddo                           
        enddo                          
!soil
!
        do jj = 1, jlistnum      
          j=jlist1(jj)            
          nxj=nxdef_2d(j)
          do i = 1,nxj
            totalp(i,jj)=0.
!            flash(i,jj)=0.   !xb110, flash density
!            ustar(i,jj)=0.1
            ustar(i,jj)=sqrt(0.14) !make sure z0 will be 0.0002 over ocean
            tstar(i,jj)=0.025
            qstar(i,jj)=0.0
            hflux(i,jj)=-ustar(i,jj)*tstar(i,jj)
            qflux(i,jj)=-ustar(i,jj)*qstar(i,jj)
          enddo
        enddo
!
        do 48 jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
        do 48 k=1,lev
        do 48 i=1,nxj
          e(i,k,jj)=1.0e-4
          eps(i,k,jj)=1.0e-7
          dtrad(i,k,jj)=0.
          asl(i,k,jj)=0.
          atl(i,k,jj)=0.
 48     continue
!
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do i=1,nxj
            if(ocean(i,jj))tg(i,jj)=sst(i,jj)
          enddo
        enddo
!
      else
!
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do i=1,nxj
            snr(i,jj)=0.
            totalp(i,jj)=0.
!            flash(i,jj)=0.   !xb110, flash density
            ustar(i,jj)=0.1
            tstar(i,jj)=0.025
            qstar(i,jj)=0.0
            hflux(i,jj)=-ustar(i,jj)*tstar(i,jj)
            qflux(i,jj)=-ustar(i,jj)*qstar(i,jj)
            curate(i,jj)=0.
            cumtop(i,jj)=1.0
            plcl(i,jj)=1.0
            ss(i,jj)=0.
            rs(i,jj)=0.
!soil        gwr(i,jj)= gwet(i,jj) * 20.0
          enddo
        enddo
!
        do 45 jj = 1, jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
        do 45 k=1,lev
        do 45 i=1,nxj
         e(i,k,jj)=1.0e-4
         eps(i,k,jj)=1.0e-7
         dtrad(i,k,jj)=0.
         asl(i,k,jj)=0.
         atl(i,k,jj)=0.
 45     continue
!
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do i=1,nxj
            if(ocean(i,jj))then
              tg(i,jj)=sst(i,jj)
              gwr(i,jj)=20.      ! soil
            else
              tg(i,jj)=tgclim(i,jj)
              gwr(i,jj)=min(1., gwclim(i,jj)/0.37)*20.  ! soil
            endif
            canopy(i,jj)=0.1   ! soil
          enddo
        enddo
!
        do jj = 1, jlistnum
          j=jlist1(jj)    
          nxj=nxdef_2d(j)
          do k = 1, km_soil 
            do i = 1, nxj
              stc(i,k,jj)=tgclim(i,jj)
              smc(i,k,jj)=gwclim(i,jj)
              slc(i,k,jj)=gwclim(i,jj)
            enddo                 
          enddo                 
        enddo
!          
      endif     ! end of  (.not.cstar)
!
!     read terrain geopotential from data base
!
      if (ksgeo.lt.0)  then
!byl        call zilch (spgeo,jtrun*jtmax*2)
        do jj = 1, jlistnum
          j=jlist1(jj)    
          nxj=nxdef_2d(j)
          do i = 1, nxj
            std(i,jj) = 0.
            sgeo(i,jj)= 0.
          enddo 
        enddo
      else
        if(ksgeo.eq.99) then
         topohgt='gbkf'
        else
         write(topohgt,'(a3,i1.1)')'gbk',ksgeo
        end if
        write(lrec,'("s00060",a4,a4,12x)')topohgt,ggdef
        call dmsread(nx,my,lrec,nxmy,'H',bckfile,ww1,istat)
        if(istat.ne.0)then
          call mpe_finalize
          call dmsexit(-1)
        endif
        do jj = 1, jlistnum
          j=jlist1(jj)
          ii=nxjstart(j)
!byl          nxj=nxdef(j)
          nxj=nxdef_2d(j)
          if( lreduce.eq.1 ) call reducepick (ww1(1,j),nxdef(j),nx,1)
          do i = 1, nxj
!byl            ww3(i,jj) = ww1(i,j)*grav
            sgeo(i,jj) = ww1(ii,j)*grav
            ii=ii+1
          enddo
        enddo


!ch     call tranrs1(jtrun,jtmax,nx,my,my_max,poly,weight,sgeo,spgeo,nsize)
!ch     call transr1(jtrun,jtmax,nx,my,my_max,poly,spgeo,sgeo,nsize)
!ch     call mpe_unify_1(ww1,sgeo,nx,my,2,mpe_double)

        call mpe2d_unify_nx(ww3,sgeo)
        call tranrs1(jtrun,jtmax,nx,my,my_max,poly,weight,ww3,spgeo,nsizey)
        call transr1(jtrun,jtmax,nx,my,my_max,poly,spgeo,sgeo,nsizey)
        call unify_reduceintp(nx,my,my_max,sgeo,ww4)
        call qmaxn3 (ww4,'sgeo',' ',1,1,1,nx,my,1)
!dms    istdno=99
!dms    istdno=0
!c      topostd='gbkf'   ! responding to istdno=99
        topostd='gbk0'   ! responding to istdno=0
        write(lrec,'("s00062",a4,a4,12x)')topostd,ggdef
        call dmsread(nx,my,lrec,nxmy,'H',bckfile,ww1,istat)
        if(istat.ne.0)then
          call mpe_finalize
          call dmsexit(-1)
        endif
        do jj = 1, jlistnum
          j=jlist1(jj)    
          ii=nxjstart(j)
          nxj=nxdef_2d(j)
          if( lreduce.eq.1 ) call reducepickr (ww1(1,j),nxdef(j),nx,1)
          do i=1,nxj
            std(i,jj)=ww1(ii,j)*ww1(ii,j)
!byl            if(std(i,jj).le.0. .or. ocean(i,jj)) std(i,jj)=0.
            if(ww1(ii,j).le.0. .or. ocean(i,jj)) std(i,jj)=0.
            ii=ii+1
          enddo
        enddo
      endif      ! end of (ksgeo.lt.0)
!
!     laplacian of terrain geopotential for divergence equation
!
!byl      do 200 m=1,mlistnum
!byl       mf=mlist(m)
!byl      do 200 n=mf,jtrun
!byl       dsqgeo(n,m,1)= spgeo(n,m,1)*eps4(n,m)
!byl       dsqgeo(n,m,2)= spgeo(n,m,2)*eps4(n,m)
!byl  200 continue
!
      if(doincr)then
       call incrini
      endif
!
      if(myrank .eq. 0) print *,'idtg=',idtg
      taux=0.
      call  sigful( nx,my,my_max,lev,ncld,lmax,jtrun,jtmax,ifilin    &
             , ifilout,cstar,ktrop,idtg,ptop,taux,capa,grav,rgas,rad &
             , cp,weight,poly,sigma,cosl,phi,tt,ut,vt,qt,o3l,pt,sgeo &
             , pdiff,tsave,t1000,plt,pk,pk2,taup                     &
             , ggdef,gmdef)
!
! reset update cycle tau,if it is abnormal
!
!      if(taup .gt. taureg)taup=6.
!
!  grid to spectral transforms for:
!
!  temperature
!  specific humidity
!  terrain pressure
!
!byl      call joinrs(cc,tt,qt,dummy,dummy,nx,my_max,lev,jlistnum,2,ncld)
!      call tranrs(jtrun,jtmax,nx,my,my_max,levp,poly,weight,cc  &
!                 ,wss,1+ncld,nsizey)
!      call ujoinrs(wss,temnow,qnow,dummy,dummy,jtrun,jtmax,levp &
!byl                 ,mlistnum,2,ncld)
      call joinrs(cc,tt,dummy,dummy,dummy,nx,my_max,lev,jlistnum,1,1)
      call tranrs(jtrun,jtmax,nx,my,my_max,levp,poly,weight,cc  &
                 ,temnow,1,nsizey)
      call mpe2d_unify_nx(ww3,pt) !2dMPI
!ch   call tranrs1(jtrun,jtmax,nx,my,my_max,poly,weight,pt     &
      call tranrs1(jtrun,jtmax,nx,my,my_max,poly,weight,ww3    &
                 ,plnow,nsizey)
!
!  compute vorticity and divergence from u and v
!
      call trandv ( jtrun,jtmax,nx,my,my_max,lev,ut,vt,weight,cim &
                   ,onocos,poly,dpoly,vornow,divnow,nsizey)
!
!  set both time levels equal at tau=0
!
      if(doincr)then
        do m=1,mlistnum
          mf=mlist(m)
          if (mf.eq.1) then
            do k = 1, levp
              divnow(k,1,1,m)= 0.0
              divnow(k,2,1,m)= 0.0
              vornow(k,1,1,m)= 0.0
              vornow(k,2,1,m)= 0.0
            enddo
          endif
        enddo
      else
        do m=1,mlistnum
          mf=mlist(m)
          if (mf.eq.1) then
            do k = 1, levp
              divnow(k,1,1,m)= 0.0
              divnow(k,2,1,m)= 0.0
              vornow(k,1,1,m)= 0.0
              vornow(k,2,1,m)= 0.0
            enddo
          endif
        enddo
        do m=1,mlistnum
          mf=mlist(m)
          do n=mf,jtrun
            do i=1,2
            do k=1,levp
              divold(k,i,n,m)= divnow(k,i,n,m)
              vorold(k,i,n,m)= vornow(k,i,n,m)
              temold(k,i,n,m)= temnow(k,i,n,m)
            enddo
            enddo
          enddo
        enddo
!byl        do m=1,mlistnum
!          mf=mlist(m)
!          do n=mf,jtrun
!            do k=1,levp*ncld*2
!              qold(k,1,n,m)  = qnow(k,1,n,m)
!            enddo
!          enddo
!byl        enddo
        do m=1,mlistnum
          mf=mlist(m)
          do n=mf,jtrun
            plold(n,m,1)= plnow(n,m,1)
            plold(n,m,2)= plnow(n,m,2)
          enddo
        enddo
      endif  !end of (doincr)
!
!     use wk1 and wk2 as work arrays to store up and vp
!
!!      call tranuv ( jtrun,jtmax,nx,my,my_max,levp,onocos,wcfac,wdfac &
!!                   ,poly,dpoly,vornow,divnow,wk1,wk2,nsizey)
!
!
!  write initial spectral coefficients to history file
!
      if(wrestrt)then
        call chlen (cwbout,48,lcwb)
        call chlen (phyout,48,lphy)
        rfile = cwbout(1:lcwb)//'tau000'
!!        allocate (work_io((( (7+2*ncld)*2*lev+8)*jtrun*jtmax+my*lev)*nsize))
!!        call gather_spec(work_io,vornow,divnow,temnow,qnow,plnow,   &
!!           vorold,divold,temold,qold,plold,dsqgeo,spgeo,trefs,      &
!!           lev,ncld,jtrun,jtmax,my,nsize)
!         if(myrank .eq. 0) then
!           open (unit=7,file=rfile,form='unformatted')
!cc         write (7) work_io
!cc         call flush (7)
!           close (7)
!         endif
!!        deallocate (work_io)
!
        rfile = phyout(1:lphy)//'tau000'
        allocate (tm1(nx,lev,my))
        allocate (tm2(nx,lev,my))
        allocate (tm3(nx,lev,my))
        allocate (tm4(nx,lev,my))
        allocate (tmc1(nx,lev,my))
        allocate (tmc2(nx,lev,my))
        allocate (tmc3(nx,lev,my))
        allocate (tmc4(nx,lev,my))
        allocate (tmc5(nx,lev,my))
        allocate (tmc6(nx,lev,my))
        allocate (temp1(nx,lev,my))
        allocate (temp2(nx,lev,my))
        allocate (temp3(nx,lev,my))
!        call unify_grid                                                   &
!            ( snr,gwr,tg,tm1,tm2,ss,rs,tm3,tm4,ustar,tstar,qstar          &
!            , hflux,qflux,raincu,rainlp,totalp,curate,plcl,cumtop         &
!            , tgclim,gwet,z0,alb,land,ice,ocean,gwclim,acld               &
!            , tmc1,tmc2,tmc3,tmc4,tmc5,tmc6,fpsp,ftp,fqp,fpsp1,ftp1,fqp1  &
!            , e,eps,o3l,dtrad,pt,ptend,ww1,ww2,nx,my,my_max,lev           &
!            , smc,slc,stc,canopy,sigmaf,istyp,ivegtyp,km_soil             &
!            , temp1,temp2,temp3,rld,zice,asl,atl)
!        if(myrank .eq. 0) then
!          open (unit=10,file=rfile,form='unformatted')
!cc        write(10) snr,gwr,tg,tm1,tm2,ss,rs,tm3,tm4,ustar,tstar,qstar
!cc  1          , hflux,qflux,raincu,rainlp,totalp,curate,plcl,cumtop
!cc  2          , tgclim,gwet,z0,alb,land,ice,ocean,gwclim,acld,sgeo
!cc  3          , temp1,temp2,canopy,sigmaf,istyp,ivegtyp,rld
!hmhj3          , tmc1,tmc2,fpsp
!cc        call flush (10)
!          close (10)
!        endif
        deallocate (tm1)
        deallocate (tm2)
        deallocate (tm3)
        deallocate (tm4)
        deallocate (tmc1)
        deallocate (tmc2)
        deallocate (tmc3)
        deallocate (tmc4)
        deallocate (tmc5)
        deallocate (tmc6)
        deallocate (temp1)
        deallocate (temp2)
        deallocate (temp3)
      endif    ! end of (wrestrt) for write initial spectral coefficients
!
      if(myrank .eq. 0) print*,' history file written at tau= 0'
!
      endif   ! end of ( .not. restrt )
!
!  spectral to real transforms of spherical harmonic coefficients of
!  model's dependent variables
!
!  vorticity
!  divergence
!  temperature
!  specific humidity ( on grid point space in NDSL version)
!  terrain pressure
!
!!      call joinsr(wss3,vornow,divnow,temnow,dummy,jtrun,jtmax,levp &
!!                 ,mlistnum,3,1)
!!      call transr(jtrun,jtmax,nx,my,my_max,levp,poly,wss3,cc3     &
!!                 ,3,nsizey)
!!      call ujoinsr(cc3,rvor,rdiv,tt,dummy,nx,my_max,lev,jlistnum,3,1)
      call transr(jtrun,jtmax,nx,my,my_max,levp,poly,vornow,cc,1,nsizey)
      call ujoinsr(cc,rvor,dummy,dummy,dummy,nx,my_max,lev,jlistnum,1,1)
      call transr(jtrun,jtmax,nx,my,my_max,levp,poly,divnow,cc,1,nsizey)
      call ujoinsr(cc,rdiv,dummy,dummy,dummy,nx,my_max,lev,jlistnum,1,1)
      call transr(jtrun,jtmax,nx,my,my_max,levp,poly,temnow,cc,1,nsizey)
      call ujoinsr(cc,tt,dummy,dummy,dummy,nx,my_max,lev,jlistnum,1,1)
      call transr1(jtrun,jtmax,nx,my,my_max,poly,plnow,pt,nsizey)
!
      do 160 jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        call prexp_hybrid_cwb ( nxjp(j),nxp,lev,ptop,sigma,pt(1,jj) &
                          ,pk(1,1,jj),pk2(1,1,jj),plt(1,1,jj) )
 160  continue
!
!  compute globel moisture budget and p-coordinate variables
!
      qgini = 0.
!
      do 420 jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        q1 = 0.0
        do 400 k = 1, lev
        do 400 i = 1, nxj
          dsigp = dsigma(k,1)*pt(i,jj)+dsigma(k,2)
          q1        = q1 + qt(i,k,jj)*dsigp
 400    continue
        qgini = qgini + q1*cosl(j)*100./grav
 420  continue
!
      call mpe_global_sum(qgini,1,mpe_double)
      if(myrank .eq. 0) print *,'in getrdy qgini=',qgini
!
!  compute the spectrally truncated velocity coefficients
!
      call tranuv ( jtrun,jtmax,nx,my,my_max,levp,onocos,wcfac,wdfac &
                   ,poly,dpoly,vornow,divnow,ut,vt,nsizey)
!
!  zonal and meridional gradients of terrain pressure
!
      call trngra ( jtrun,jtmax,nx,my,my_max,cim,poly,dpoly,plnow, &
                    dlpl,dtpl,nsizey)
!
!  compute spline coefficients for long wave radiation calculation
!  compute lontitude and latitude of grid points for zenth angle
!  calculation in by radiation parameterization. those array are
!  passed to subroutine diabat by common block /radcon/
!
      njump1 = njump + 1
      if ( mod( nxp,njump1) .ne. 0 )  njump1 = njump
      njump2 = njump + 2
      if ( mod( nxp,njump2) .ne. 0 )  njump2 = njump1
      njump3 = njump + 3
      if ( mod( nxp,njump3) .ne. 0 )  njump3 = njump2
!
      if ( lreduce .eq. 1 ) then
        njump1=njump
        njump2=njump
        njump3=njump
      endif
!
      lvlw = nxp/njump
      lvlw1= nxp/njump1
      lvlw2= nxp/njump2
      lvlw3= nxp/njump3
      call splinc (lvlw ,nxp,il(1,1),ib(1,1),cof(1,1))
      call splinc (lvlw1,nxp,il(1,2),ib(1,2),cof(1,2))
      call splinc (lvlw2,nxp,il(1,3),ib(1,3),cof(1,3))
      call splinc (lvlw3,nxp,il(1,4),ib(1,4),cof(1,4))
!
      if(myrank .eq. 0) then
      print *, "  from getrdy: njump= ",njump,njump1,njump2,njump3
      endif
!
      pi = 4.0*atan(1.0)
      do 520 j = 1, my
        xlat(j) = asin(sinl(j))*180./pi
 520  continue

      do jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef(j)
        xlon(1,jj)=0.
        do i=2,nxj
          xlon(i,jj)=xlon(1,jj)+float(i-1)*360./nxj
          if(xlon(i,jj).gt.180. .and. xlon(i,jj).lt.360.)then
            xlon(i,jj)=-180.0+abs(xlon(i,jj)-180.)
          else if(xlon(i,jj).ge.360.)then
            xlon(i,jj)=xlon(i,jj)-360.
          endif
        enddo
      enddo


!---------------------------------
! read forecast sst
!---------------------------------
      IF(ldailyFCTsst .OR. ldailyFCTicesndpt .OR. dailyClm_option .ge.1 ) THEN
        CALL read_dailyFCT(idtg,taui,dt,tg,cice,sndepth,xlon,xlat,ocean)
      ENDIF

!---------------------------------
!0.0 initial_sit
!---------------------------------
      if(do_sit) then
        CALL set_ocndepth()
        if(myrank .eq. 0) print *,'end set_ocndepth'
        CALL allocate_sitgrid_array(nxp,my_max)
        if(myrank .eq. 0) print *,'end allocate_sitgrid_array'
        if(lwoa0) CALL read_woa0
        if(myrank .eq. 0) print *,'end read_woa0'
        IF (lgodas) then
          IF(ldailysst)then
            CALL read_dailygodas(idtg,taui,dt)
            if(myrank .eq. 0) print *,'end read_dailygodas'
          ELSE
            CALL read_godas(idtg)
            if(myrank .eq. 0) print *,'end read_godas'
          ENDIF
        ENDIF
        if(locaf0) then
          if(myrank.eq.0) then
            print *,'ready in call ifilin_ocaf'
          endif
          call read_ocaf0(nx,my,lkvl,ggdef,idtg)
        endif

        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do ii=1,nxj
            i=nxjstart(j)+ii-1
            sitmask(ii,jj)=0.
!    !  1.0 set geological data
            sitcor(ii,jj)   = 2*(7.292e-5)*sin(xlat(j)*api/180.)
            sitlat(ii)      = xlat(j)
            IF(xlon(i,jj) .LT. 0.) then
              sitlon(ii,jj)=xlon(i,jj)+360.
            ELSE
              sitlon(ii,jj)=xlon(i,jj)
            ENDIF
            call cal_ratioBlending(myrank,ii,jj,sitlat(ii),sitlon(ii,jj),ratioSIT(ii,jj))


!    !  2.0 set sst, sss and sic
            obswtb(ii,jj)   = tg(ii,jj)
            obsseaice(ii,jj)= cice(ii,jj)
            seaice(ii,jj)   = cice(ii,jj)
            t_surf         = obswtb(ii,jj)
            if ( obsseaice(ii,jj) > 0.5 ) then
              thickness(ii,jj) = specified_ice_thickness
!              ice_mask (ii,jj) = .true.
              t_surf          = MIN( t_surf, ctfreez )
            else
              thickness(ii,jj) = 0.0
!              ice_mask (ii,jj) = .false.
              t_surf          = MAX( t_surf, ctfreez )
            endif
            tsi(ii,jj)= MIN( t_surf, ctfreez )

!    !  4. set additional SIT input variables
            if(land(ii,jj) ) then
              sitlclass(ii,jj)= 1.                   ! land landclass
              bathy(ii,jj)    = 0.                   ! 0 m for the first guess, need to read terrain data later
              wlvl(ii,jj)     = bathy(ii,jj)-1.       ! set water level at 1 m below the bathy for land grids
              ocnmask(ii,jj)  = xmissing             ! not coupled to 3-D ocn model
              obox_mask(ii,jj)= xmissing             ! not coupled to 3-D ocn model
              sni(ii,jj)      = 0.                   ! assuming initially no snow over seaice (m swe)
              slm(ii,jj)      = 1.                   ! land fraction
              obswsb(ii,jj)   = 0.                   ! set observed SSS at 0 PSU over land water
              tsl(ii,jj)      = tg(ii,jj)
              tslm(ii,jj)     = tg(ii,jj)
              tslm1(ii,jj)    = tg(ii,jj)
            endif
            if(ice(ii,jj) ) then                     ! sea ice
              sitlclass(ii,jj)= 2.                   ! water landclass
              bathy(ii,jj)    = bathydepth                ! 0 m for the first guess, need to read terrain data later
              wlvl(ii,jj)     = 0.                   ! set water level at 1 m below the bathy for land grids
              ocnmask(ii,jj)  = xmissing             ! not coupled to 3-D ocn model
              obox_mask(ii,jj)= xmissing             ! not coupled to 3-D ocn model
              sni(ii,jj)      = 0.                   ! assuming initially no snow over seaice (m swe)
              slm(ii,jj)      = 0.                   ! land fraction
              obswsb(ii,jj)   = 36.3                 ! set observed SSS at 0 PSU over land water
              tsl(ii,jj)      = (1-cice(ii,jj))*tg(ii,jj)+cice(ii,jj)*tsi(ii,jj)
              tslm(ii,jj)     = tsl(ii,jj)
              tslm1(ii,jj)    = tsl(ii,jj)
            endif
            if(ocean(ii,jj) )then
              sitlclass(ii,jj)= 2.                   ! water landclass
              bathy(ii,jj)    = bathydepth                ! 200 m depth for the first guess, need to read terrain data later)
              wlvl(ii,jj)     = 0.                   ! set water level at 0 m
              ocnmask(ii,jj)  = xmissing             ! not coupled to 3-D ocn model
              obox_mask(ii,jj)= xmissing             ! not coupled to 3-D ocn model
              sni(ii,jj)      = 0.                   ! assuming no snow over seaice (m swe). It can be read from NCEP data.
              slm(ii,jj)      = 0.                   ! land fraction
              obswsb(ii,jj)   = 36.3                 ! set observed SSS at 36.3 PSU
              tsl(ii,jj)      = tg(ii,jj)
              tslm(ii,jj)     = tg(ii,jj)
              tslm1(ii,jj)    = tg(ii,jj)
!!Ocea n within 40N-40S
             if( (sitlat(ii).GE.(sit_domain_s-sit_domain_extgrd)) &
               .AND. (sitlat(ii).LE.(sit_domain_n+sit_domain_extgrd)) ) then
               if ((sitlon(ii,jj).GE.(sit_domain_w-sit_domain_extgrd))&
               .AND. (sitlon(ii,jj).LE.(sit_domain_e+sit_domain_extgrd)))then

                 sitmask(ii,jj)=1.
                 if(locaf0 .and. mask1st(ii,jj) .eq. 0.) sitmask(ii,jj)=0.
               endif
             endif
            endif

            ctfreez2(ii,jj) = tmelts(obswsb(ii,jj))
            tsw(ii,jj)      = tg(ii,jj)

          enddo



!      ! 5.0 set additional SIT ocn profle t,s,u,v and tke
         call time_weights(idtg,0.)
         lsitstart=.true.
          call sit_vdiff_init ( nxjp(j), nxp, jj, j,                   &
!       ! 0-INPUT only, original ATM/SIT variabels
             sitlat, sitlon(:,jj),                                     &
             sitmask(:,jj), bathy(:,jj), wlvl(:,jj),                   &
             ocnmask(:,jj), obox_mask(:,jj),                           &
             sni(:,jj), thickness(:,jj), tsi(:,jj),                    &
             obsseaice(:,jj), obswtb(:,jj), obswsb(:,jj),              &
             ctfreez2(:,jj),                                           &
!       ! 2-d SIT vars
             sitwtb(:,jj), sitwub(:,jj), sitwvb(:,jj),                 &
             sitwsb(:,jj),                                             &
             subfluxw(:,jj), wsubsal(:,jj),                            &
             sitcc(:,jj), sithc(:,jj), engwac(:,jj),                   &
             sc(:,jj), saltwac(:,jj),                                  &
             wtfns(:,jj), wsfns(:,jj),                                 &
!       ! 3-d SIT vars: snow/ice
             zsi(:,jj,0:1), silw(:,jj,0:1), tsnic(:,jj,0:3),           &
!       ! 3-d SIT vars: water column
             obswt(:,jj,0:lkvl+1), obsws(:,jj,0:lkvl+1), obswu(:,jj,0:lkvl+1), &
             obswv(:,jj,0:lkvl+1),                                     &
             sitwt(:,jj,0:lkvl+1), sitwu(:,jj,0:lkvl+1), sitwv(:,jj,0:lkvl+1), &
             sitww(:,jj,0:lkvl+1), sitws(:,jj,0:lkvl+1),               &
             sitwtke(:,jj,0:lkvl+1), wlmx(:,jj,0:lkvl+1),   &
             wldisp(:,jj,0:lkvl+1), wkm(:,jj,0:lkvl+1), wkh(:,jj,0:lkvl+1), &
             wrho1000(:,jj,0:lkvl+1),                                  &
             wtfn(:,jj,0:lkvl+1), wsfn(:,jj,0:lkvl+1),                 &
             wtfn0(:,jj,0:lkvl+1), wsfn0(:,jj,0:lkvl+1),               &
             awufl(:,jj,0:lkvl+1), awvfl(:,jj,0:lkvl+1), awtfl(:,jj,0:lkvl+1), &
             awsfl(:,jj,0:lkvl+1), awtfl0(:,jj,0:lkvl+1),awsfl0(:,jj,0:lkvl+1),&
             awtkefl(:,jj,0:lkvl+1),             &
!       ! 4- OUTPUT only, original ATM variabels
             seaice(:,jj),                                             &
             grndcapc(:,jj), grndhflx(:,jj), grndflux(:,jj),           &
             sftobswt(:,jj,0:lkvl+1) )

        enddo


        if(locaf) then
          if(myrank.eq.0) then
            print *,'ready in call ifilin_ocaf'
          endif
          call read_ocaf(nx,my,lkvl,ggdef)
        endif


        if(restrt) then
          call readrerun_sitgrid1(itaui)
          call readrerun_sitgrid2(itaui)
          call readrerun_sitgrid3(itaui)
          restrt= .false.
        endif

        if(lpre6hr_sit) then
          call readpre6hr_sit(nx,my,itaup,ifilin,idtg,ggdef)
        endif


      endif   !end of(do_sit)

!------------------------------
! end do_sit
!------------------------------
!----
! for u10 v10 t2 being output at tau=0 (6/20/2003)
!
      if( .not. restrt ) then
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
!
!         (ut, vt) : scaled velocity
!         tt : virtual potential temperature
!         tg : real temperature
!         (ux, vx) : real velocity at the lowest sigma level
!         tx : potential temperature at the lowest sigma level
!         tsx: surface skin real temperature
!         tgx: potential temperature at the ground
!
          xx = rad/cosl(j)
          do i = 1, nxj
            ps(i) = pt(i,jj) + ptop
            hs(i) = cp*tt(i,lev,jj)*(pk2(i,lev,jj)-pk(i,lev,jj))/grav
            ux(i) = ut(i,lev,jj) * xx  ! back to real velocity
            vx(i) = vt(i,lev,jj) * xx  ! back to real velocity
            tx(i) = tt(i,lev,jj) / (1.0+0.608*qt(i,lev,jj))
            qx(i) = qt(i,lev,jj)
            tsx(i) = tt(i,lev,jj) * pk2(i,lev,jj) / (1.0+0.608*qt(i,lev,jj))
            tgp(i) = tg(i,jj)/pk2(i,lev,jj)   ! tg in term of potential temp
          enddo
          call qsatq ( nxjp(j), tg(1,jj), ps, qs )
          do i = 1, nxj
            if(ocean(i,jj) .or. istyp(i,jj).eq.0)then
              qsfc(i) = qs(i)
            else
              wet = ( smc(i,1,jj)-wltsmc(istyp(i,jj)) ) / &
                    ( refsmc(istyp(i,jj))-wltsmc(istyp(i,jj)) )
              qsfc(i) = wet*qs(i)+(1.0-wet)*qx(i)
              qsfc(i) = min(qs(i),qsfc(i) )
            endif
          enddo
          call sfcuvt( nxjp(j),nxp,dt,grav,rgas,cp,hltm                 &
                    ,tgp,z0(1,jj),ocean(1,jj),ps,tsx                    &
                    ,hs,ux,vx,tx,qx,ustar(1,jj),tstar(1,jj),qstar(1,jj) &
                    ,hflux(1,jj),qflux(1,jj),qsfc,t2(1,jj),q2(1,jj)     &
                    ,rh2(1,jj),rh10(1,jj),u10(1,jj),v10(1,jj) )
        enddo
!        call mpe_unify(ustar,nx,my,2,mpe_double)
!        call mpe_unify(tstar,nx,my,2,mpe_double)
!        call mpe_unify(qstar,nx,my,2,mpe_double)
!        call mpe_unify(t2,nx,my,2,mpe_double) 
!        call mpe_unify(rh2,nx,my,2,mpe_double)
!        call mpe_unify(u10,nx,my,2,mpe_double) 
!        call mpe_unify(v10,nx,my,2,mpe_double) 
!        call mpe_unify(hflux,nx,my,2,mpe_double)
!        call mpe_unify(qflux,nx,my,2,mpe_double)
      endif    ! end of ( .not. restrt ) for u10 v10 t2 being output at tau=0

! for Thompson : 1st guess number concentration where mass non-zero
      if ( .not.restrt .and. nmmiph.eq.18 ) then
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do k = 1, lev
            do i = 1, nxj
              ! virtual temperature :
              ttv = tt(i,k,jj)*pk(i,k,jj)
              ! real temperature :
              ttr = ttv/(1.0+0.608*qt(i,k,jj))
              ! air density :
              rho = plt(i,k,jj)*100./(con_rd*ttv)

              tem = qt(i,(ntiw-1)*lev+k,jj)
              if ( tem .gt. 0. ) then
                 qt(i,(ntinc-1)*lev+k,jj) =                             &
                      make_IceNumber(tem*rho,ttr)/rho
              endif

              tem = qt(i,(ntrw-1)*lev+k,jj)
              if ( tem .gt. 0. ) then
                 qt(i,(ntrnc-1)*lev+k,jj) =                             &
                      make_RainNumber(tem*rho,ttr)/rho
              endif
            enddo
          enddo
        enddo
      endif
!
!  output initial fileds
!
      if (.not.restrt)then
        wk1 = 0.
        ww1 = 0.
        ww2 = 0.
        raintot=0.
        raincu=0.
        rainlp=0.
        raincu6=0.
        rainlp6=0.
        gfx=0.
        sld=0.
        rld=0.
!        flash=0.   !xb110, flash density

!!       open grib2 file
        if( outgrb2 == 1 .and. myrank == 0 )then
          if(io_quilting)then 
              grbnxmy=nx*my
              write( keydoit,'(A14,I12.12,A8)') &
              "OPEN..0000....",idtg,"H...DOIT"
              ntag=ntag+1
              call mpe_send_key(keydoit,ntag,istat)
          else
            grbid=233  ! 231 outflds  232 out24  233 mfc
 133                    format( A  ,A ,I10.10 ,A       )
            write(grbfile,133 )trim(ifilout_grb),'/GFS_',idtg/100 ,'_0000.grb2'
            if(myrank==0) print*,'OutFileName= ',trim(grbfile)
            call opn_grb2(nx,my,idtg, 0 ,istat)
          endif
        endif

        call outflds ( 0,nx,my,my_max,lev,ncld,lmax,numout,idtg,ifilout &
             , outdir,ktrop,ptop,capa,cp,rgas,grav,sigma,sgeo           &
             , ptend,pt,plt,pk,pk2,phi,ut,vt,vvel                       &
             , tt,qt,rdiv,rvor,tg,gwr,z0,hflux,qflux,snr                &
             , raintot,raincu,rainlp,plcl,cumtop,ss,rs,alb,gwclim       &
             , acld,cosl,wk1,ww2,ww2,t2,q2,rh2,rh10,u10,v10,gfx,rld,sld &
!byl             , km_soil,smc,slc,stc,canopy,ggdef,slp,v850,v700,h850,h500 &
             , km_soil,smc,slc,stc,canopy,ggdef,typtrk                  &
!             , ctot,chig,cmid,clow,hpbl,.true.,flash,do_sit)
             , ctot,chig,cmid,clow,hpbl,.true.,do_sit)

! add 40m 100m output for green energy plan

      if(out_green)then

        
        call  outflds_green(0,nx,my,my_max,lev,ncld                     &
              , idtg,ifilout,cp,rgas,grav,t2,u10,v10,ss,pk           &
              , sgeo,pt,plt,ptop,ut,vt,tt,qt,cosl,raincu6,rainlp6       &
              , ggdef)
      endif
!
!#ifdef RSM_sigp
        if(outgrb2==1.and.myrank==0)then
            if(io_quilting)then 
              keydoit(1:4)='CLSE'
              ntag=ntag+1
              call mpe_send_key(keydoit,ntag,istat)
            else
              call cls_grb2(istat)
            endif
        endif

#ifdef RSM
       if(outrsm) then
        if(myrank.eq.0)print*,' output: rsm date',idtg
        write(dtgrsm,'(I12)') idtg
        read(dtgrsm,'(I10,I2)')idtgrsm,ii   ! ii is dummy integer
#ifdef CWB_MPMD
        call send_idate(idtgrsm)
#else
        call wrte_idate(idtgrsm)
#endif
        call rsmout_sigp( itaui,nx,my,my_max,lev,ncld        &
                     , idtg,ptop,rad,grav,cosl           &
                     , pt,sgeo,snr,gwr,tg,pk             &
                     , ut,vt,tt,qt,km_soil,smc,stc       &
                     , ice,land,ocean,xlon,xlat)
       endif
#endif
!#ifdef RSM
!      if (outrsm) then
!        if(myrank.eq.0)print*,' output: rsm date',idtg
!        write(dtgrsm,'(I12)') idtg
!        read(dtgrsm,'(I10,I2)')idtgrsm,ii   ! ii is dummy integer
!#ifdef CWB_MPMD
!        call send_idate(idtgrsm)
!#else
!        call wrte_idate(idtgrsm)
!#endif
!        call rsmout(idtg,0,nx,my,my_max,lev,ncld      &
!                , ptop,cp,rgas,grav,sgeo,pdiff        &
!                , t1000,pt,plt,pk,pk2,phi,ut,vt       &
!                , tt,qt,tg,snr,cosl                   &
!                , km_soil,smc,stc                     &
!                , ice,land,ocean)
!      endif
!#endif
!
!
        if(typhoon)then
          do n=1,ntyph
            i=ixtyp(1,n)
            j=jytyp(1,n)
          do m=1,5 !(1:slp 2:v850 3:v700 4:h850 5:h500)
            call unify_reduceintp(nx,my,my_max,typtrk(1,1,m),ww4)
            tensity(0,m,n)=( ww4(i,j+1)+ww4(i+1,j+1)    &
                           + ww4(i,j  )+ww4(i+1,j  ) )/4.
          enddo
!byl            tensity(0,2,n)=( v850(i,j+1)+v850(i+1,j+1)  &
!byl                           + v850(i,j  )+v850(i+1,j  ) )/4.
!byl            tensity(0,3,n)=( v700(i,j+1)+v700(i+1,j+1)  &
!byl                           + v700(i,j  )+v700(i+1,j  ) )/4.
!byl            tensity(0,4,n)=( h850(i,j+1)+h850(i+1,j+1)  &
!byl                           + h850(i,j  )+h850(i+1,j  ) )/4.
!byl            tensity(0,5,n)=( h500(i,j+1)+h500(i+1,j+1)  &
!byl                           + h500(i,j  )+h500(i+1,j  ) )/4.
          enddo
        endif

      endif   ! end of (.not.restrt) for output initial fileds
!
      return
      end
