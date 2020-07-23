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

      implicit   none

!  local working array
!
      real      sst(nxp,my_max),ww1(nx,my),ww2(nxp,my_max),     &
                rh2100(nxp,my_max),rh10100(nxp,my_max),         &
                wk1(nxp,lev,my_max),pklev(nxp,my_max),          &
                cc(nx+2,levp,1,my_max),ww3(nx,my_max)
!byl                wss3(levp,2,3,jtrun,jtmax),cc3(nx+2,levp,3,my_max)

      character lrec*26,rfile*55,ctau*6,topostd*4,topohgt*4
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
      real    fact,xxaa,taux,dummy,q1,dsigp,pi,xx,wet
!xb110>
!      real    flash(nxp,my_max)
!xb110<

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
!!         uzm,lev,ncld,jtrun,jtmax,my,nsize)
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
!byl      if( lreduce.eq.1 ) call reducepick (ww1,nxdef,nx,my)
      do jj=1,jlistnum
        j=jlist1(jj)
        ii=nxjstart(j)
        nxj=nxdef_2d(j)
        if( lreduce.eq.1 ) call reducepick (ww1(1,j),nxdef(j),nx,1)
        do i=1,nxj
          pdiff(i,jj)=ww1(ii,j)
          ii=ii+1
        enddo
      enddo

      call syslbl('h00100',idtg,0,ggdef,lrec)
      call dmsread(nx,my,lrec,lncrec,'H',ifilin,ww1,istat)
!byl      if( lreduce.eq.1 ) call reducepick (ww1,nxdef,nx,my)
      do jj=1,jlistnum
        j=jlist1(jj)
        ii=nxjstart(j)
        nxj=nxdef_2d(j)
        if( lreduce.eq.1 ) call reducepick (ww1(1,j),nxdef(j),nx,1)
        do i=1,nxj
          t1000(i,jj)=ww1(ii,j)
          ii=ii+1
        enddo
      enddo

      call syslbl('x00tsv',idtg,0,ggdef,lrec)
      call dmsread(nx,my,lrec,lncrec,'H',ifilout,ww1,istat)
!byl      if( lreduce.eq.1 ) call reducepick (ww1,nxdef,nx,my)
      do jj=1,jlistnum
        j=jlist1(jj)
        ii=nxjstart(j)
        nxj=nxdef_2d(j)
        if( lreduce.eq.1 ) call reducepick (ww1(1,j),nxdef(j),nx,1)
        do i=1,nxj
          tsave(i,jj)=ww1(ii,j)
          ii=ii+1
        enddo
      enddo
!
!  read dmsdata for standard deviation of terrain field
!
!dms    istdno=99
!dms    istdno=0
!c      topostd='gbkf'   ! responding to istdno=99
        topostd='gbk0'   ! responding to istdno=0
        write(lrec,'("s00062",a4,a4,12x)')topostd,ggdef
        call dmsread(nx,my,lrec,nxmy,'H',bckfile,ww1,istat)
!byl        if( lreduce.eq.1 ) call reducepick (ww1,nxdef,nx,my)
!
        do jj=1,jlistnum
          j=jlist1(jj)
          ii=nxjstart(j)
          nxj=nxdef_2d(j)
        if( lreduce.eq.1 ) call reducepick (ww1(1,j),nxdef(j),nx,1)
          do i=1,nxj
            std(i,jj)=ww1(ii,j)*ww1(ii,j)
!byl            if(std(i,jj).le.0. .or. ocean(i,jj)) std(i,jj)=0.
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
!byl        if( lreduce.eq.1 ) call reducepick (ww1,nxdef,nx,my)
        do jj=1,jlistnum
          j=jlist1(jj)
          ii=nxjstart(j)
          nxj=nxdef_2d(j)
        if( lreduce.eq.1 ) call reducepick (ww1(1,j),nxdef(j),nx,1)
          do i=1,nxj
            sst(i,jj)=ww1(ii,j)
            ii=ii+1
          enddo
        enddo
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
!byl          if( lreduce.eq.1 ) call reducepick (ww1,nxdef,nx,my)
          do jj=1,jlistnum
            j=jlist1(jj)
            ii=nxjstart(j)
            nxj=nxdef_2d(j)
            if( lreduce.eq.1 ) call reducepick (ww1(1,j),nxdef(j),nx,1)
            do i=1,nxj
              snr(i,jj)=ww1(ii,j)
              ii=ii+1
            enddo
          enddo
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
!byl          if( lreduce.eq.1 ) call reducepick (ww1,nxdef,nx,my)
          do jj=1,jlistnum
            j=jlist1(jj)
            ii=nxjstart(j)
            nxj=nxdef_2d(j)
            if( lreduce.eq.1 ) call reducepick (ww1(1,j),nxdef(j),nx,1)
            do i=1,nxj
              cice(i,jj)=ww1(ii,j)
              ii=ii+1
            enddo
          enddo
          if( myrank .eq. 0 ) &
             print*,"get ncep's sea ice analysis, at dtg=",idtg
          if ( ncepicthk ) then
            call syslbl('w00092',idtg,0,ggdef,lrec)
            call dmsread(nx,my,lrec,nxmy,'H',ifilin,ww1,istat)
!byl          if( lreduce.eq.1 ) call reducepick (ww1,nxdef,nx,my)
            do jj=1,jlistnum
              j=jlist1(jj)
              ii=nxjstart(j)
              nxj=nxdef_2d(j)
              if( lreduce.eq.1 ) call reducepick(ww1(1,j),nxdef(j),nx,1)
              do i=1,nxj
                zice(i,jj)=ww1(ii,j)
                ii=ii+1
              enddo
            enddo
            if( myrank .eq. 0 ) &
               print*,"get ncep's sea ice thickness analysis, at dtg=",idtg
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
          print*,"==================   Warnig!!!   ==================="
          print*,"=  ice thickness not consistent with sea ice mask  ="
          print*,"=  set the thickness to 1 meter for first guess    ="
          print*,"================================= =================="
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
!byl            ustar(i,jj)=0.1
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
!byl        if( lreduce.eq.1 ) call reducepick (ww1,nxdef,nx,my)
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
        call mpe2d_unify(ww1,sgeo)

        call qmaxn3 (ww1,'sgeo',' ',1,1,1,nx,my,1)
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
!byl        if( lreduce.eq.1 ) call reducepick (ww1,nxdef,nx,my)
        do jj = 1, jlistnum
          j=jlist1(jj)    
          ii=nxjstart(j)
          nxj=nxdef_2d(j)
          if( lreduce.eq.1 ) call reducepick (ww1(1,j),nxdef(j),nx,1)
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
             , pdiff,tsave,t1000,plt,pk,pk2,trefs,taup               &
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
!!      call uzmean ( nx,my,my_max,lev,wk1,uzm )
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
!!           uzm,lev,ncld,jtrun,jtmax,my,nsize)
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
!
      pi = 4.0*atan(1.0)
      do 520 j = 1, my
        xlat(j) = asin(sinl(j))*180./pi
 520  continue
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
       rh2100=rh2*100.
       rh10100=rh10*100.
!        flash=0.   !xb110, flash density
        call outflds ( 0,nx,my,my_max,lev,ncld,lmax,numout,idtg,ifilout &
             , outdir,ktrop,ptop,capa,cp,rgas,grav,sigma,sgeo           &
             , ptend,pt,plt,pk,pk2,phi,ut,vt,sd                         &
             , tt,qt,rdiv,rvor,tg,gwr,z0,hflux,qflux,snr                &
             , raintot,raincu,rainlp,plcl,cumtop,ss,rs,alb,gwclim       &
             , acld,cosl,wk1,ww2,ww2,t2,q2,rh2100,rh10100,u10,v10,gfx,rld,sld &
!byl             , km_soil,smc,slc,stc,canopy,ggdef,slp,v850,v700,h850,h500 &
             , km_soil,smc,slc,stc,canopy,ggdef,typtrk                  &
!             , ctot,chig,cmid,clow,hpbl,.true.,flash,do_sit)
             , ctot,chig,cmid,clow,hpbl,.true.,do_sit)

! add 40m 100m output for green energy plan
      if(out_green)then
          do jj = 1, jlistnum
            j=jlist1(jj)
            nxj=nxdef_2d(j)
            do i = 1,nxj
              pklev(i,jj) = pk(i,lev,jj)
            enddo
          enddo
        call  outflds_green(0,nx,my,my_max,lev,ncld                     &
              , idtg,ifilout,cp,rgas,grav,t2,u10,v10,ss,pklev           &
              , sgeo,pt,plt,ptop,ut,vt,tt,qt,cosl,raincu6,rainlp6       &
              , ggdef)
      endif
!
!
        if(typhoon)then
          do n=1,ntyph
            i=ixtyp(1,n)
            j=jytyp(1,n)
          do m=1,5 !(1:slp 2:v850 3:v700 4:h850 5:h500)
            call unify_reduceintp(nx,my,my_max,typtrk(1,1,m),ww1)
            tensity(0,m,n)=( ww1(i,j+1)+ww1(i+1,j+1)    &
                           + ww1(i,j  )+ww1(i+1,j  ) )/4.
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
