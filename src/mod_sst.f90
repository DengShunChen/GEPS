#define LCWBGFS T
#if defined(LCWBGFS)
#define nlon nx
#define ngl my
#define nlat my
#define mpp_root_pe() 0
#define p_parallel_io (myrank .eq. 0)
#define p_pe myrank
#define myrank_check 55 
#define jj_check 4
#define ii_check 1
#endif

!    read wtfn12, wsfn12 data
      MODULE mod_sst

        use param
        use mpe
        use index
        use const,             ONLY:ggdef,ifilin_nc,ifilin_sst,ifilin_ncep  &
                                   ,ldailyFCTsst,ldailyFCTicesndpt,ifilin   &
                                   ,dailyClm_option,ifilin_ClmANA,ifilin_ClmFCT
        use mod_sit_control,   ONLY: lwarning_msg,xmissing,lwoa0,lsitstart &
                                    ,lamip,lmixedlayer,ngodas,nwoa0     &
                                    ,lgodas,ldailysst,locaf0 
        USE mo_netcdf,         ONLY: IO_info_print,FILE_INFO,NETCDF     &
                                    ,io_open,io_open_unit,io_close      &
                                    ,lkvl,sit_zdepth,message_text       
        USE mod_eos_ocean,      ONLY:API,IDAYLEN
    
        IMPLICIT NONE

        PUBLIC :: ifilin_ocaf
        PUBLIC :: read_godas, read_woa0, read_ocaf,read_dailygodas,read_dailyFCT
        PUBLIC :: read_ocaf0
        PUBLIC :: nodepth, odepths, ot12, os12, ou12, ov12,mixedlayer12
        PUBLIC :: nodepth0, odepth0, ot0, os0, ou0, ov0, mixedlayer0  
        PUBLIC :: lou, lov
        PUBLIC :: nwdepth, wdepths, wtfn12, wsfn12,wtfn1st,wsfn1st
        PUBLIC :: albice,albsn,albw,csn,cice,rhosn,rhoice,xkice,xksn,xkw,         &
                  omegas,wcri
        PUBLIC :: deallocate_ocaf_array,deallocate_woa0_array,deallocate_godas_array
!ps				  ,tol,wlvlref,dpthmx,init_sit_ocean

        PUBLIC :: opgsst
        PUBLIC :: nmw1,nmw2,wgt1,wgt2
        PUBLIC :: now1,now2,wgto1,wgto2             !! GODAS MONTHLY/PENTAD Data
        PUBLIC :: obswtbnmw1,obswtbnmw2,obswtbwgt1,obswtbwgt2     !obswtb Data
        PUBLIC :: dailyFCTsst,dailyFCTcice,dailyFCTsndepth
        PUBLIC :: ANAsstT0,dailyClmANAsst,dailyClmFCTsst
        PUBLIC :: deallocate_dailyFCT_array

        INCLUDE 'netcdf.inc'


        character*60, save:: ifilin_ocaf ='OCAFDMS'

        !! memory pointer for GODAS+Ishii WORLD OCEAN data (lgodas) (GODAS+Ishii)
        INTEGER               :: nodepth        ! number of depths of the godas data (=24)
        REAL, ALLOCATABLE :: odepths(:)     ! depths of the godas data (m)  
!!!  INTEGER               :: gpdepth        ! number of depths of daily fodas pentad data (=40)
!!!  REAL(dp), ALLOCATABLE :: odepths(:)    ! depths of the godas pentad data (m)  
        REAL, ALLOCATABLE :: ot12(:,:,:,:) ! (nlon,nodepth,ngl,0:13) in global coordinates,
                                           ! observed water tempeature profile (K): "ot"       
        REAL, ALLOCATABLE :: os12(:,:,:,:) ! (nlon,nodepth,ngl,0:13) in global coordinates,
                                           ! observed salinity (0/00): "os"
        REAL, ALLOCATABLE :: ou12(:,:,:,:) ! (nlon,nodepth,ngl,0:13) in global coordinates,
                                           ! observed u current (m/s): "ou"
        REAL, ALLOCATABLE :: ov12(:,:,:,:) ! (nlon,nodepth,ngl,0:13) in global coordinates,
                                           ! observed v current (m/s): "ou"
        REAL, ALLOCATABLE :: mixedlayer12(:,:,:) ! (nlon,nodepth,ngl,0:13) in global coordinates,
                                           ! observed v current (m/s): "ou"
                                          
        !! memory pointer for Initial WORLD OCEAN ATLAS 2005 data (lwoa0)(http://www.nodc.noaa.gov/OC5/WOA05/pr_woa05.html)
        INTEGER               :: nodepth0         ! number of depths of the woa data (=24)
        REAL, ALLOCATABLE :: odepth0(:)       ! depths of the woa data (m)  
        REAL, ALLOCATABLE :: ot0(:,:,:)    ! (nlon,nodepth0,ngl,0:13) in global coordinates,
                                           ! observed water tempeature profile (K): "ot"       
        REAL, ALLOCATABLE :: os0(:,:,:)    ! (nlon,nodepth0,ngl,0:13) in global coordinates,
                                           ! observed salinity (0/00): "os"
        REAL, ALLOCATABLE :: ou0(:,:,:)    ! (nlon,nodepth0,ngl,0:13) in global coordinates,
                                           ! observed u-componet current (m/s): "ou"
        REAL, ALLOCATABLE :: ov0(:,:,:)    ! (nlon,nodepth0,ngl,0:13) in global coordinates,
                                           ! observed v-componet current (m/s): "ov"                                        
        REAL, ALLOCATABLE :: mixedlayer0(:,:) ! (nlon,ngl) in global coordinates,
                                           ! observed mixedlayer (m): "ou"
        
        LOGICAL :: lou=.FALSE.                    ! u- current available ?
        LOGICAL :: lov=.FALSE.                    ! v- current available ?

  !! memory pointer for sit flux correction term
        INTEGER           :: nwdepth         ! number of depths of the woa data
        REAL, ALLOCATABLE :: wdepths(:)      ! depths of the woa data (m)  
        REAL, ALLOCATABLE :: wtfn12(:,:,:,:) ! (nlon,nwdepth,ngl,0:13) in global coordinates,
                                             ! observed water tempeature profile (K): "ot"       
        REAL, ALLOCATABLE :: wsfn12(:,:,:,:) ! (nlon,nwdepth,ngl,0:13) in global coordinates,
        REAL, ALLOCATABLE :: wtfn1st(:,:,:) ! (nlon,nwdepth,ngl) in global coordinates,
                                             ! observed water tempeature profile (K): "ot"       
        REAL, ALLOCATABLE :: wsfn1st(:,:,:) ! (nlon,nwdepth,ngl) in global coordinates,

  !! memory pointer for pentad GODAS OCEAN data (lgodas & ldailysst)
        REAL      :: timevals_godas(3) = 0.  ! absoulte time (e.g., 19971003.25) GODAS PENTAD Data
        INTEGER   :: files_godas(3)          ! fileid for pentad data at day-1, day+0 and day+1, respectively
        INTEGER   :: tsID_godas(3)           ! record id for pentad data at day-1, day+0 and day+1, in its repective file, respectively
        INTEGER   :: nts_godas(3)            ! # of timestamps for pentad data at day-1, day+0 and day+1, in its repective file, respectively

  !! memory pointer for daily FCTsst data (ldailyFCTsst)
        REAL      :: timevals_dailyFCT(3)= 0.   ! absoulte time (e.g.,19971003.25) daily forecast sst Data
        REAL, ALLOCATABLE :: dailyFCTsst(:,:,:)    ! (nlon,ngl,2) at ydate, ydate+1 day in global coordinates,
                                                   ! forecast daily water tempeature (K)
        REAL, ALLOCATABLE :: dailyFCTcice(:,:,:)    ! (nlon,ngl,2) at ydate, ydate+1 day in global coordinates,
                                                   ! forecast sea ice fration
        REAL, ALLOCATABLE :: dailyFCTsndepth(:,:,:)    ! (nlon,ngl,2) at ydate, ydate+1 day in global coordinates,
                                                   ! forecast daily snow depth (mm)
        REAL, ALLOCATABLE :: ANAsstT0(:,:)       ! (nlon,ngl), analysis SST at tau=0 
        REAL, ALLOCATABLE :: dailyClmANAsst(:,:,:)    ! (nlon,ngl,2) at tau=0, ydate, ydate+1 day in global coordinates,
                                                   ! reanalysis daily climatology water tempeature (K)
        REAL, ALLOCATABLE :: dailyClmFCTsst(:,:,:)    ! (nlon,ngl,2) at ydate, ydate+1 day in global coordinates,
                                                   ! forecast daily climatology water tempeature (K)
       

  !*    1.0 COEFFICIENTS IN sit_ocean MODEL

        REAL :: albice = 0.2,             & ! albedo of glacier ice (0.2 - 0.4) (Pielke, 1984). Lake ice is the smallest.
                albsn  = 0.7,             & ! albedo of snow decreases with age (0.4 - 0.95) (Pielke, 1984)
                albw   = 0.06,            & !	  (Brutsaert, 1982; Tsuang, 1990; Gaspar et al., 1990)
                csn    = 2116.,           &
                cice   = 2116.,           & ! Cw = 4217.7-2.55*(Tw-tmelt) (4178.4 is used)
                                      ! csn = cice=104.369+7.369*TSN (Marks, 1988)
                                      ! csn = 2116 j kg k-1 at 0 ! is used for simplicity.
                rhosn  = 300.,            & !	rhosn = dry snow density. Although snow density changes with age,
                                      !	it is not sensitive to snowmelt runoff. A constant value is assumed.
                rhoice = 917.,            & !	rhoice = density of ice (917 kg/m3)
                omegas = 2.*API/IDAYLEN,  & ! ANGULAR VELOCITY OF EARTH in respect to sun
                xkice = 1.2E-6,           & ! Dickinson et al. (1986)
                xksn = 2.0E-7,            & !	  Effect heat diffusivity of snow
                                      !       (conduction heat diffusivity + vapor transfer)
                                      ! xksn=kcon+De*Lv/cp*dqsat/dT
                                      !	kcon (conduction) = (3.2238E-8)*rhosn/csn (Yeh, 1965)
                                      ! This value increases with age of snow, and is an
                                      ! important tunning parameter for snow melting rate.
                                      !   =(1E-7 ~ 4E-7) (Dickinson et al., 1986)
                xkw = 1.50E-4,            & !  (Kondo, 1979, Tsuang, 1990)
                wcri = 0.1                    !     minimum thickness of a water layer. Thickness less than
                                      !     this thickness is treated as thin layer. That is T,s,U,V,
                                      !     TKE is the same as the layer underneath.
  
      INTEGER, PARAMETER :: nerr = 6     ! error output stream


!for opgsst
      REAL, dimension(:,:,:),allocatable,save :: opgsst
!for time_interpolation
      REAL wgt1,wgt2,obswtbwgt1,obswtbwgt2,wgto1,wgto2
      INTEGER nmw1,nmw2,obswtbnmw1,obswtbnmw2,now1,now2


        
      CONTAINS

!---------------------------------------------------------
        subroutine allocate_opgsst_array

          integer  ierr
          allocate ( opgsst(nxp,my_max,0:13), stat=ierr)
             if (ierr/= 0) then
                 write(6,*) 'mod_opgsst : allocate fail 1 '
                 stop
             end if

             return

        end subroutine


        subroutine deallocate_opgsst_array

           deallocate (opgsst)
           return

        end subroutine

        SUBROUTINE read_opgsst(idtg1,ggdef,ocean,ice)

          use index
          use rank

          character*4 ggdef
          integer*8 idtg1
          character*12 cdtg
          integer   iyyyy,mm,ddhhmn,iyy,imm,istat
          integer lncrec
          character lrec*26
          integer i,j,ii,jj,nxj,kmm
          real,dimension(:,:), allocatable:: temp1
          logical ocean(nxp,my_max),ice(nxp,my_max)
 

          allocate(temp1(nx,my))
!
  11     format('W00100','0000',a4,i4.4,i2.2,'010000')  ! W10

!           if(myrank.eq.0) call dmsopn(ifilin_sst,"r",istat)
!           call mpe_broadcast(istat,1,flag,mpe_integer)
!           if(istat.ne.0)then
!             if(myrank.eq.0) print *,'dmsopn sst error'
!               call mpe_finalize
!               call dmsexit(-1)
!           endif


            write(cdtg,'(i12)') idtg1
            read(cdtg,'(i4,i2,i8)')iyyyy,mm,ddhhmn

!              iy=idtg_sst/1000000

           lncrec=nx*my

           do kmm=0,13
             if( kmm.eq.0 ) then
               imm=12
               iyy=iyyyy-1
             elseif (kmm.gt.12) then
               imm=1
               iyy=iyyyy+1
             else
               imm=i
               iyy=iyyyy
             endif

             write(lrec,11) ggdef,iyy,imm
             call dmsread(nx,my,lrec,lncrec,'H',ifilin_sst,temp1(:,:),istat)
!             if( lreduce.eq.1 ) call reducepick(temp1(1,imm,1),nxdef,nx,my)
             do jj = 1, jlistnum
               j=jlist1(jj)
               i=nxjstart(j)
               nxj=nxdef_2d(j)
               if( lreduce.eq.1 )call reducepick (temp1(1,j),nxdef(j),nx,1)
               do ii = 1, nxj
                 i=nxjstart(j)+ii-1
                 if(ocean(ii,jj) .or. ice(ii,jj))then
                  opgsst(ii,jj,kmm)=temp1(i,j)+273.16
                 endif
               end do
             end do

           enddo

           deallocate(temp1)

        END SUBROUTINE read_opgsst

!----------------------------------------------------------

        subroutine deallocate_dailyFCT_array

           if(ldailyFCTsst) deallocate (dailyFCTsst)
           if(ldailyFCTicesndpt) then
             deallocate (dailyFCTcice) 
             deallocate (dailyFCTsndepth)
           endif
           if(dailyClm_option .ge. 1) then
             deallocate (ANAsstT0)
             deallocate (dailyClmANAsst)
             if(dailyClm_option .eq. 2) then
               deallocate (dailyClmFCTsst)
             endif
           endif 
           return

        end subroutine deallocate_dailyFCT_array


        SUBROUTINE read_dailyFCT(idtg1,tau,dtx,tg1,cice1,sndepth1)

          USE index
          USE rank

          INTEGER*8 :: idtg1,idtg_temp,idtg_FCT
          REAL      :: tau,dtx,tauleft
          REAL      :: tautemp
          INTEGER   :: icurrenttau
          REAL      :: ydate,ydate2
          INTEGER       :: yr, mo, dy, hr, mn
          character*12 cdtg
          real tg1(nxp,my_max),cice1(nxp,my_max),sndepth1(nxp,my_max)
           

          icurrenttau=int(tau)
          tauleft=float(int((tau-int(tau)+0.001)*3600./dtx))*dtx   !(sec)
          if(tauleft .eq. 3600.) then
          icurrenttau=icurrenttau+1
          tauleft=0.
          endif
          call dtgfix12(idtg1,idtg_temp,icurrenttau)
          write(cdtg,'(i12)')idtg_temp
          read(cdtg,'(i4,i2,i2,i2,i2)')yr,mo,dy,hr,mn
          ydate=float(yr*10000+mo*100+dy)+float(hr)/24.+tauleft/3600./24.
          if(myrank .eq. 0) print *,"read_dailyFCT,ydate=",ydate


          tautemp=tau+24.     !tau +24 hr
          icurrenttau=int(tautemp)
          tauleft=float(int((tautemp-int(tautemp)+0.001)*3600./dtx))*dtx !(sec)
          if(tauleft .eq. 3600.) then
          icurrenttau=icurrenttau+1
          tauleft=0.
          endif
          call dtgfix12(idtg1,idtg_FCT,icurrenttau)
          write(cdtg,'(i12)')idtg_FCT
          read(cdtg,'(i4,i2,i2,i2,i2)')yr,mo,dy,hr,mn
!          ydate2=float(yr*10000+mo*100+dy)+float(hr/24)
          ydate2=float(yr*10000+mo*100+dy)
          write(cdtg,'(i4,i2.2,i2.2,i2.2,i2.2)')yr,mo,dy,00,00
          read(cdtg,'(i12)')idtg_FCT
          if(myrank .eq. 0) print *,"read_dailyFCT,ydate2=",ydate2

          IF(ldailyFCTsst) THEN
            IF (.NOT. ALLOCATED(dailyFCTsst)) ALLOCATE (dailyFCTsst(nxp,my_max,2))
          ENDIF
          IF(ldailyFCTicesndpt) THEN
            IF (.NOT. ALLOCATED(dailyFCTcice)) ALLOCATE (dailyFCTcice(nxp,my_max,2))
            IF (.NOT. ALLOCATED(dailyFCTsndepth)) ALLOCATE (dailyFCTsndepth(nxp,my_max,2))
          ENDIF
          IF(dailyClm_option .ge. 1) THEN
            IF (.NOT. ALLOCATED(ANAsstT0)) ALLOCATE (ANAsstT0(nxp,my_max))
            IF (.NOT. ALLOCATED(dailyClmANAsst)) ALLOCATE (dailyClmANAsst(nxp,my_max,0:2))
            if(dailyClm_option .eq. 2) then
              IF (.NOT. ALLOCATED(dailyClmFCTsst)) ALLOCATE (dailyClmFCTsst(nxp,my_max,0:2))
            endif
          ENDIF


          IF ( (tau .eq. 0.) .OR. lsitstart) THEN
     !!! warm/cold start
            timevals_dailyFCT(0)=ydate
            timevals_dailyFCT(1)=ydate
            if(ldailyFCTsst) dailyFCTsst(:,:,1)=tg1(:,:)
            if(dailyClm_option .ge. 1) ANAsstT0(:,:)=tg1(:,:)
            if(ldailyFCTicesndpt) then
              dailyFCTcice(:,:,1)=cice1(:,:)
              dailyFCTsndepth(:,:,1)=sndepth1(:,:)
            endif
            if(ldailyFCTsst .or. ldailyFCTicesndpt) CALL read_dailyFCT_dayp1(idtg_FCT)
            if(dailyClm_option .ge. 1)then
              CALL read_dailyClm_2days(idtg1,idtg_FCT,icurrenttau,dailyClm_option)
            endif

            timevals_dailyFCT(2)=ydate2
            if( myrank .eq. myrank_check) then
              print*,"read_dailyFCT: tg1=",tg1(ii_check,jj_check)
              if(ldailyFCTsst) then
                print*,"1. dailyFCTsst(ii_check,jj_check,1)=" &
                      , dailyFCTsst(ii_check,jj_check,1)      &
                      ,",dailyFCTsst(ii_check,jj_check,2)="   &
                      , dailyFCTsst(ii_check,jj_check,2)
              endif
              if(dailyClm_option .ge. 1)then
                print*,"dailyClmANAssst(",ii_check,",",jj_check,",0)="  &
                      , dailyClmANAsst(ii_check,jj_check,0)     &
                      ,",dailyClmANAsst(",ii_check,",",jj_check,",1)="  &
                      , dailyClmANAsst(ii_check,jj_check,1)     &
                      ,",dailyClmANAsst(",ii_check,",",jj_check,",2)="  &
                      , dailyClmANAsst(ii_check,jj_check,2)
                if(dailyClm_option .eq. 2)then
                  print*,",dailyClmFCTsst(",ii_check,",",jj_check,",0)="&
                       , dailyClmFCTsst(ii_check,jj_check,0)            &
                       ,",dailyClmFCTsst(",ii_check,",",jj_check,",1)=" &
                       , dailyClmFCTsst(ii_check,jj_check,1)             &
                       ,",dailyClmFCTsst(",ii_check,",",jj_check,",2)=" &
                       , dailyClmFCTsst(ii_check,jj_check,2)
                endif
              endif
            endif
          ELSEIF (ydate.LT.timevals_dailyFCT(2)) THEN
     ! data were read. Note that initial value of  timevals_godas=0.
            RETURN
          ELSE
     ! note that initial value of  timevals_godas=0.
     ! Shift left
            if(ldailyFCTsst) dailyFCTsst(:,:,1)=dailyFCTsst(:,:,2)
            if(ldailyFCTicesndpt) then
              dailyFCTcice(:,:,1)=dailyFCTcice(:,:,2)
              dailyFCTsndepth(:,:,1)=dailyFCTsndepth(:,:,2)
            endif
            timevals_dailyFCT(1)=timevals_dailyFCT(2)
            if(ldailyFCTsst .or. ldailyFCTicesndpt) then
              CALL read_dailyFCT_dayp1(idtg_FCT)    ! read next FCST data
            endif
            if(dailyClm_option .ge. 1)then
              dailyClmANAsst(:,:,1)=dailyClmANAsst(:,:,2)
              if(dailyClm_option .eq. 2) then
                dailyClmFCTsst(:,:,1)=dailyClmFCTsst(:,:,2)
              endif
              CALL read_dailyClm_dayp1(idtg1,idtg_FCT,icurrenttau,dailyClm_option)
            endif

            timevals_dailyFCT(2)=ydate2
            if( myrank .eq. myrank_check) then
              print*,"read_dailyFCT: tg1=",tg1(ii_check,jj_check)
              if(ldailyFCTsst) then
                print*,"1. dailyFCTsst(ii_check,jj_check,1)=" &
                      , dailyFCTsst(ii_check,jj_check,1)      &
                      ,",dailyFCTsst(ii_check,jj_check,2)="   &
                      , dailyFCTsst(ii_check,jj_check,2)
              endif
              if(dailyClm_option .ge. 1)then
                print*,"dailyClmANAssst(",ii_check,",",jj_check,",0)="  &
                      , dailyClmANAsst(ii_check,jj_check,0)     &
                      ,",dailyClmANAsst(",ii_check,",",jj_check,",1)="  &
                      , dailyClmANAsst(ii_check,jj_check,1)     &
                      ,",dailyClmANAsst(",ii_check,",",jj_check,",2)="  &
                      , dailyClmANAsst(ii_check,jj_check,2)
                if(dailyClm_option .eq. 2)then
                  print*,",dailyClmFCTsst(",ii_check,",",jj_check,",0)="&
                       , dailyClmFCTsst(ii_check,jj_check,0)            &
                       ,",dailyClmFCTsst(",ii_check,",",jj_check,",1)=" &
                       , dailyClmFCTsst(ii_check,jj_check,1)            &
                       ,",dailyClmFCTsst(",ii_check,",",jj_check,",2)=" &
                       , dailyClmFCTsst(ii_check,jj_check,2)
                endif
              endif
            endif

          ENDIF

      CONTAINS
     !-----------------------------------
        SUBROUTINE read_dailyFCT_dayp1(idtg1)

          INTEGER*8 idtg1
          INTEGER iyyyy,imm,idd,ihh,imn
          INTEGER lncrec
          character lrec*26
          REAL ssttemp(nx,my),cicetemp(nx,my),sndpttemp(nx,my)
          INTEGER ii,jj,nxj,i,j,istat,nxjtot
          INTEGER itemp,jtemp
          INTEGER iitemp,jjtemp
          REAL sst_suntemp,sst_counttemp 
          REAL cice_suntemp,cice_counttemp,sndpt_suntemp,sndpt_counttemp

          write(cdtg,'(i12)') idtg1
          read(cdtg,'(i4.4,i2.2,i2.2,i2.2,i2.2)') iyyyy,imm,idd,ihh,imn

          lncrec=nx*my
          
          ssttemp=0.
          cicetemp=0.
          sndpttemp=0.

   11     format('W00100','0000',a4,i4.4,i2.2,i2.2,i2.2,'00')  ! sea surface temperature
          write(lrec,11) ggdef,iyyyy,imm,idd,ihh
          call dmsread(nx,my,lrec,lncrec,'H',ifilin_sst,ssttemp(:,:),istat)
          if(ldailyFCTicesndpt)then
   12       format('W00091','0000',a4,i4.4,i2.2,i2.2,i2.2,'00')  ! sea ice fration
            write(lrec,12) ggdef,iyyyy,imm,idd,ihh
            call dmsread(nx,my,lrec,lncrec,'H',ifilin_ncep,cicetemp(:,:),istat)
   13       format('B00650','0000',a4,i4.4,i2.2,i2.2,i2.2,'00')  ! water equivlent snow depth
            write(lrec,13) ggdef,iyyyy,imm,idd,ihh
            call dmsread(nx,my,lrec,lncrec,'H',ifilin_ncep,sndpttemp(:,:),istat)
          endif

!          if( lreduce.eq.1 ) then
!            call reducepick(ssttemp(1,1),nxdef,nx,my)
!            if(ldailyFCTicesndpt)then
!              call reducepick(cicetemp(1,1),nxdef,nx,my)
!              call reducepick(sndpttemp(1,1),nxdef,nx,my)
!            endif
!          endif 

          DO jj = 1, jlistnum
            j=jlist1(jj)
            i=nxjstart(j)
            nxj=nxdef_2d(j)
            if( lreduce.eq.1 ) then
              call reducepick (ssttemp(1,j),nxdef(j),nx,1)
              if(ldailyFCTicesndpt)then
                call reducepick (cicetemp(1,j),nxdef(j),nx,1)
                call reducepick (sndpttemp(1,j),nxdef(j),nx,1)
              endif
            endif

            DO ii=1,nxj
              i=nxjstart(j)+ii-1
              if( myrank .eq. myrank_check .AND. i .eq. ii_check .AND. jj .eq. jj_check) then
                 print*,"read_dailyFCT_dayp1: ssttemp(",ii_check,"," &
                       ,j,")=",ssttemp(ii_check,j)
              endif

              dailyFCTsst(ii,jj,2)=MERGE(ssttemp(i,j),xmissing,  &
                   (ssttemp(i,j).GE.271. .AND. ssttemp(i,j) .LE. 400.))
              if(ldailyFCTicesndpt)then
                dailyFCTcice(ii,jj,2)=MERGE(cicetemp(i,j),xmissing,  &
                                         (cicetemp(i,j) .GE. 0))
                dailyFCTsndepth(ii,jj,2)=MERGE(sndpttemp(i,j),xmissing, &
                                         (sndpttemp(i,j) .GE. 0))
              endif

              if( (ssttemp(i,j).EQ.xmissing) .OR. (ldailyFCTicesndpt .AND. &
                ((cicetemp(i,j).EQ.xmissing).OR.(sndpttemp(i,j).EQ.xmissing))) ) then
                sst_suntemp=0.
                sst_counttemp=0.
                cice_suntemp=0.
                cice_counttemp=0.
                sndpt_suntemp=0.
                sndpt_counttemp=0.
                nxjtot=nxdef(j)
                do itemp=i-1, i+1
                  do jtemp=j-1, j+1
                    if (itemp .eq. 0) then
                      iitemp=nxjtot
                    elseif (itemp .gt. nxjtot) then
                      iitemp= 1
                    else
                      iitemp=itemp
                    endif
                    if (jtemp .eq. 0) then
                      jjtemp=my
                    elseif (jtemp .gt. my) then
                      jjtemp=1
                    else
                      jjtemp=jtemp
                    endif
                    if (ssttemp(iitemp,jjtemp).GE.271. .AND. ssttemp(iitemp,jjtemp) .LE. 400.)then
                      sst_suntemp=sst_suntemp+ssttemp(iitemp,jjtemp)
                      sst_counttemp=sst_counttemp+1.
                    endif
                    if(ldailyFCTicesndpt)then
                      if (cicetemp(iitemp,jjtemp).GE.0.)then
                        cice_suntemp=cice_suntemp+cicetemp(iitemp,jjtemp)
                        cice_counttemp=cice_counttemp+1.
                      endif
                      if (sndpttemp(iitemp,jjtemp).GE.0.)then
                        sndpt_suntemp=sndpt_suntemp+sndpttemp(iitemp,jjtemp)
                        sndpt_counttemp=sndpt_counttemp+1.
                      endif
                    endif
                  enddo
                enddo
                if(ssttemp(i,j).LT.271. .OR. ssttemp(i,j) .GT. 400.)then
                  ssttemp(i,j)=sst_suntemp/sst_counttemp
                  if (ssttemp(i,j).GE.271. .AND. ssttemp(i,j) .LE. 400.)then
                    dailyFCTsst(ii,jj,2)=ssttemp(i,j)
                  else
                    dailyFCTsst(ii,jj,2)=xmissing
                  endif
                endif
                if(ldailyFCTicesndpt)then
                  if(cicetemp(i,j).LT. 0.)then
                    cicetemp(i,j)=cice_suntemp/cice_counttemp
                    if(cicetemp(i,j).GE. 0.)then
                      dailyFCTcice(ii,jj,2)=cicetemp(i,j)
                    endif
                  endif
                  if(sndpttemp(i,j).LT. 0.)then
                    sndpttemp(i,j)=sndpt_suntemp/sndpt_counttemp
                    if(sndpttemp(i,j).GE. 0.)then
                      dailyFCTsndepth(ii,jj,2)=sndpttemp(i,j)  
                    endif
                  endif
                endif
              endif

              if (myrank .eq. myrank_check .AND. ii .eq. ii_check .AND. jj.eq. jj_check) then
                print*,"myrank=",myrank,",i=",i,",j=",j   &
                      ,",ii=",ii, ",jj=",jj               &
                      ,",ssttemp(i,j)=",ssttemp(i,j)      &
                      ,",dailyFCTsst(ii,jj,2)=",dailyFCTsst(ii,jj,2)
              endif

            ENDDO  !end do ii
          ENDDO    !end do jj


          if( myrank .eq. myrank_check) then
            print*,"read_dailyFCT_dayp1: dailyFCTsst(",ii_check,",",jj_check,",2)=" &
                   ,dailyFCTsst(ii_check,jj_check,2)
          endif
        
        END SUBROUTINE read_dailyFCT_dayp1


        SUBROUTINE read_dailyClm_2days(idtg1,idtg_2,itau,ioption)

          INTEGER*8 idtg1,idtg_2
          INTEGER itau
          INTEGER iyyyy,imm,idd,ihh,imn
          INTEGER iyyyy2,imm2,idd2,ihh2,imn2
          INTEGER lncrec
          character lrec*26
          REAL sstANA0(nx,my),sstANA1(nx,my)
          REAL sstFCT0(nx,my),sstFCT1(nx,my)
          INTEGER i,j,ii,jj,nxj,istat
          INTEGER ioption

          write(cdtg,'(i12)') idtg1
          read(cdtg,'(i4.4,i2.2,i2.2,i2.2,i2.2)') iyyyy,imm,idd,ihh,imn
          write(cdtg,'(i12)') idtg_2
          read(cdtg,'(i4.4,i2.2,i2.2,i2.2,i2.2)') iyyyy2,imm2,idd2,ihh2,imn2

          lncrec=nx*my
          sstANA0=0.
          sstANA1=0.
          sstFCT0=0.
          sstFCT1=0.

   11     format('W00100',4x,a4,4x,i2.2,i2.2,4x)  ! sea surface temperature
          write(lrec,11) ggdef,imm,idd
          call dmsread(nx,my,lrec,lncrec,'H',ifilin_ClmANA,sstANA0(:,:),istat)

   12     format('W00100',4x,a4,4x,i2.2,i2.2,4x)  ! sea surface temperature
          write(lrec,12) ggdef,imm2,idd2
          call dmsread(nx,my,lrec,lncrec,'H',ifilin_ClmANA,sstANA1(:,:),istat)

!          if( lreduce.eq.1 ) then
!            call reducepick(sstANA0(1,1),nxdef,nx,my)
!            call reducepick(sstANA1(1,1),nxdef,nx,my)
!          endif


          sstANA0=MERGE(sstANA0,xmissing,(sstANA0.GE.271. .AND. sstANA0.LE.400.))
          sstANA1=MERGE(sstANA1,xmissing,(sstANA1.GE.271. .AND. sstANA1.LE.400.))

          CALL fill_missing2(sstANA0(:,:),nx,my,1,.FALSE.)
          CALL fill_missing2(sstANA1(:,:),nx,my,1,.FALSE.)


          if(ioption .eq. 2) then
   13       format('W00100',i4.4,a4,4x,i2.2,i2.2,4x)  ! sea surface temperature
            write(lrec,13) itau,ggdef,imm,idd
            call dmsread(nx,my,lrec,lncrec,'H',ifilin_ClmFCT,sstFCT0(:,:),istat)

   14       format('W00100',i4.4,a4,4x,i2.2,i2.2,4x)  ! sea surface temperature
            write(lrec,14) itau,ggdef,imm2,idd2
            call dmsread(nx,my,lrec,lncrec,'H',ifilin_ClmFCT,sstFCT1(:,:),istat)


!            if( lreduce.eq.1 ) then
!              call reducepick(sstFCT0(1,1),nxdef,nx,my)
!              call reducepick(sstFCT1(1,1),nxdef,nx,my)
!            endif 

            sstFCT0=MERGE(sstFCT0,xmissing,(sstFCT0.GE.271. .AND. sstFCT0.LE.400.))
            sstFCT1=MERGE(sstFCT1,xmissing,(sstFCT1.GE.271. .AND. sstFCT1.LE.400.))

            CALL fill_missing2(sstFCT0(:,:),nx,my,1,.FALSE.)
            CALL fill_missing2(sstFCT1(:,:),nx,my,1,.FALSE.)
          endif


          DO jj = 1, jlistnum
            j=jlist1(jj)
            nxj=nxdef_2d(j)
            i=nxjstart(j)
            if( lreduce.eq.1 ) then
              call reducepick(sstANA0(1,j),nxdef(j),nx,1)
              call reducepick(sstANA1(1,j),nxdef(j),nx,1)
              call reducepick(sstFCT0(1,j),nxdef(j),nx,1)
              call reducepick(sstFCT1(1,j),nxdef(j),nx,1)
            end if
            DO ii=1,nxj
              i=nxjstart(j)+ii-1
              dailyClmANAsst(ii,jj,0)=sstANA0(i,j)
              dailyClmANAsst(ii,jj,1)=sstANA0(i,j)
              dailyClmANAsst(ii,jj,2)=sstANA1(i,j)
              if(ioption .eq. 2) then
                dailyClmFCTsst(ii,jj,0)=sstFCT0(i,j)
                dailyClmFCTsst(ii,jj,1)=sstFCT0(i,j)
                dailyClmFCTsst(ii,jj,2)=sstFCT1(i,j)
              endif

              if (myrank.eq.myrank_check .AND. i.eq.ii_check .AND. jj.eq.jj_check) then
                print*,"myrank=",myrank,",i=",i,",j=",j                &
                  ,",ii=",ii,",jj=",jj                                 &
                  ,",dailyClmANAsst(ii,jj,0)=",dailyClmANAsst(ii,jj,0)   &
                  ,",dailyClmANAsst(ii,jj,1)=",dailyClmANAsst(ii,jj,1)   &
                  ,",dailyClmANAsst(ii,jj,2)=",dailyClmANAsst(ii,jj,2)
                if(ioption .eq. 2) then
                  print*,",dailyClmFCTsst(ii,jj,0)=",dailyClmFCTsst(ii,jj,0)   &
                    ,",dailyClmFCTsst(ii,jj,1)=",dailyClmFCTsst(ii,jj,1)   &
                    ,",dailyClmFCTsst(ii,jj,2)=",dailyClmFCTsst(ii,jj,2)
                endif
              endif
            ENDDO  !end do ii
          ENDDO    !end do jj
        
        END SUBROUTINE read_dailyClm_2days

               
        SUBROUTINE read_dailyClm_dayp1(idtg1,idtg_2,itau,ioption)

          INTEGER*8 idtg1,idtg_2
          INTEGER itau
          INTEGER iyyyy,imm,idd,ihh,imn
          INTEGER iyyyy2,imm2,idd2,ihh2,imn2
          INTEGER lncrec
          character lrec*26
          REAL sstANA(nx,my),sstFCT(nx,my)
          INTEGER i,j,ii,jj,nxj,istat
          INTEGER ioption

          write(cdtg,'(i12)') idtg1
          read(cdtg,'(i4.4,i2.2,i2.2,i2.2,i2.2)') iyyyy,imm,idd,ihh,imn
          write(cdtg,'(i12)') idtg_2
          read(cdtg,'(i4.4,i2.2,i2.2,i2.2,i2.2)') iyyyy2,imm2,idd2,ihh2,imn2

          lncrec=nx*my

          sstANA=0.          
          sstFCT=0.          


   11     format('W00100',4x,a4,4x,i2.2,i2.2,4x)  ! sea surface temperature
          write(lrec,11) ggdef,imm2,idd2
          call dmsread(nx,my,lrec,lncrec,'H',ifilin_ClmANA,sstANA(:,:),istat)

!          if( lreduce.eq.1 ) then
!            call reducepick(sstANA(1,1),nxdef,nx,my)
!          endif

          sstANA=MERGE(sstANA,xmissing,(sstANA.GE.271. .AND. sstANA.LE.400.))
          CALL fill_missing2(sstANA(:,:),nx,my,1,.FALSE.)



          if(ioption .eq. 2) then
   13       format('W00100',i4.4,a4,4x,i2.2,i2.2,4x)  ! sea surface temperature
            write(lrec,13) itau,ggdef,imm,idd
            call dmsread(nx,my,lrec,lncrec,'H',ifilin_ClmFCT,sstFCT(:,:),istat)

!            if( lreduce.eq.1 ) then
!              call reducepick(sstFCT(1,1),nxdef,nx,my)
!            endif

            sstFCT=MERGE(sstFCT,xmissing,(sstFCT.GE.271. .AND. sstFCT.LE.400.))
            CALL fill_missing2(sstFCT(:,:),nx,my,1,.FALSE.)
          endif


          DO jj = 1, jlistnum
            j=jlist1(jj)
            nxj=nxdef_2d(j)
            if( lreduce.eq.1 ) then 
              call reducepick(sstANA(1,j),nxdef(j),nx,1)
              call reducepick(sstFCT(1,j),nxdef(j),nx,1)
            end if
            DO ii=1,nxj
              i=nxjstart(j)+ii-1
              dailyClmANAsst(ii,jj,2)=sstANA(i,j)
              if(ioption .eq. 2) then
                dailyClmFCTsst(ii,jj,2)=sstFCT(i,j)
              endif

              if (myrank.eq.myrank_check .AND. ii.eq.ii_check .AND. jj.eq.jj_check) then
                print*,"myrank=",myrank,",i=",i,",j=",j              &
                  ,",ii=",ii,",jj=",jj                                   &
                  ,",dailyClmANAsst(ii,jj,2)=",dailyClmANAsst(ii,jj,2)
                if(ioption .eq. 2)then
                  print*,",dailyClmFCTsst(ii,jj,2)=",dailyClmFCTsst(ii,jj,2)
                endif
              endif
            ENDDO
          ENDDO

        
        END SUBROUTINE read_dailyClm_dayp1
           

      END SUBROUTINE read_dailyFCT

!---------------------------------------------------------------
        subroutine deallocate_ocaf_array  

           deallocate (wdepths,wtfn12,wsfn12)
           if(locaf0) deallocate (wtfn1st,wsfn1st)
           return
        
        end subroutine


!------------------------------------------------------------------------------
        subroutine deallocate_woa0_array

           deallocate (odepth0)
           deallocate (ot0,os0,ou0,ov0,mixedlayer0)
           return

        end subroutine

!------------------------------------------------------------------------------
        subroutine deallocate_godas_array  
    
           deallocate (odepths)
           deallocate (ot12,os12,ou12,ov12)
           if(lmixedlayer) deallocate(mixedlayer12)
           return
        
        end subroutine

!------------------------------------------------------------------------------		
        SUBROUTINE read_ocaf(nx,my,lkvl,ggdef)
        
          use index
          use rank
        
        
          integer nx,my,lkvl
          character*4 ggdef
!          integer*8 idtg1
!          character*12 cdtg
!          integer   iyyyy,mmddhhmn,iyy,imm,istat
          integer istat,k,mm,lncrec
          integer j,jj,nxj,i,ii
          character lrec*26
          REAL, ALLOCATABLE, TARGET :: temp1(:,:,:),temp2(:,:,:)
          logical flag
!       
          if(.not. ALLOCATED(temp1)) allocate (temp1(nx,12,my))
          if(.not. ALLOCATED(temp2)) allocate (temp2(nx,12,my))

          if(.not. ALLOCATED(wdepths)) ALLOCATE(wdepths(1:lkvl+2)) 
          if(.not. ALLOCATED(wtfn12)) ALLOCATE(wtfn12(nxp,1:lkvl+2,my_max,0:13)) 
          if(.not. ALLOCATED(wsfn12)) ALLOCATE(wsfn12(nxp,1:lkvl+2,my_max,0:13)) 

          nwdepth=lkvl+2
          wdepths(1:lkvl+2)=sit_zdepth(0:lkvl+1)
          wtfn12=0.
          wsfn12=0.

         
  11     format(i3.3,'TFM','gbck',a4,4x,i2.2,6x)  ! ???TFM
  12     format(i3.3,'SFM','gbck',a4,4x,i2.2,6x)  ! ???SFM

          flag=.false.
          if(myrank.eq.0) then
            call dmsmsg("ALL",istat)
            print *,'ready to open ifilin_ocaf'
            call dmsopn(ifilin_ocaf,"r",istat)
            flag=.true.
          endif
!          call mpe_broadcast(istat,1,flag,mpe_integer)
          call mpe_bcast(istat,1,0,mpe_integer)
          if(istat.ne.0)then
            if(myrank.eq.0) print *,'dmsopn ocaf error'
            stop
            call mpe_finalize
            call dmsexit(-1)
          endif


          do k=1,lkvl+2
        
!           if(myrank.eq.0) call dmsopn(ifilin_sst,"r",istat)
!           call mpe_broadcast(istat,1,flag,mpe_integer)
!           if(istat.ne.0)then
!             if(myrank.eq.0) print *,'dmsopn sst error'
!               call mpe_finalize
!               call dmsexit(-1)
!           endif	  
        

!            write(cdtg,'(i12)') idtg1
!            read(cdtg,'(i4,i8)')iyyyy,mmddhhmn
        
!	       iy=idtg_sst/1000000
        
           lncrec=nx*my
           do mm=1,12
             write(lrec,11) k-1,ggdef,mm
             if(myrank .eq. 0) print *, 'lrec11=',lrec
             call dmsread(nx,my,lrec,lncrec,'H',ifilin_ocaf,temp1(:,mm,:),istat)
             write(lrec,12) k-1,ggdef,mm
             if(myrank .eq. 0) print *, 'lrec12=',lrec
             call dmsread(nx,my,lrec,lncrec,'H',ifilin_ocaf,temp2(:,mm,:),istat)
             if( myrank .eq. 72) then
              print *,"ocaf 1: wtfn(914,",mm,",265)=",temp1(914,mm,265)
              print *,"ocaf 2: wsfn(914,",mm,",265)=",temp2(914,mm,265)
             endif

!             if( lreduce.eq.1 ) call reducepick (temp1(1,mm,1),nxdef,nx,my) 
!             if( lreduce.eq.1 ) call reducepick (temp2(1,mm,1),nxdef,nx,my)
           enddo


           do jj=1,jlistnum
              j=jlist1(jj)
              nxj=nxdef_2d(j)
              do mm = 1, 12
                if( lreduce.eq.1 )then
                  call reducepick (temp1(1,mm,j),nxdef(j),nx,1)
                  call reducepick (temp2(1,mm,j),nxdef(j),nx,1)
                endif
                do ii = 1, nxj
                  i=nxjstart(j)+ii-1
                  wtfn12(ii,k,jj,mm)=temp1(i,mm,j)
                  wsfn12(ii,k,jj,mm)=temp2(i,mm,j)

                  if(mm .eq. 1) then
                    wtfn12(ii,k,jj,13)=temp1(i,mm,j)
                    wsfn12(ii,k,jj,13)=temp2(i,mm,j)
                  endif
                  if(mm .eq. 12) then
                    wtfn12(ii,k,jj,0)=temp1(i,mm,j)
                    wsfn12(ii,k,jj,0)=temp2(i,mm,j)
                  endif


                enddo  !end do ii

              enddo  !end do mm
            enddo    !end do jj

         enddo !end of k

         deallocate(temp1,temp2)
         call dmscls(ifilin_ocaf,istat)
   
        END SUBROUTINE read_ocaf

!------------------------------------------------------------------------------
        SUBROUTINE read_ocaf0(nx,my,lkvl,ggdef,idtg1)

          use index
          use rank


          integer nx,my,lkvl
          character*4 ggdef
          integer*8 idtg1
          character*12 cdtg
!          integer   iyyyy,mm,ddhhmn,iyy,imm,istat
          integer istat,k,lncrec
          integer j,jj,nxj,i,ii
          character lrec*26
          REAL, ALLOCATABLE, TARGET :: temp1(:,:),temp2(:,:)
          logical flag
!
          if(.not. ALLOCATED(temp1)) allocate (temp1(nx,my))
          if(.not. ALLOCATED(temp2)) allocate (temp2(nx,my))

          if(.not. ALLOCATED(wdepths)) ALLOCATE(wdepths(1:lkvl+2))
          if(.not. ALLOCATED(wtfn1st)) ALLOCATE(wtfn1st(nxp,1:lkvl+2,my_max))
!          if(.not. ALLOCATED(wsfn0)) ALLOCATE(wsfn0(nxp,1:lkvl+2,my_max))

          nwdepth=lkvl+2
          wdepths(1:lkvl+2)=sit_zdepth(0:lkvl+1)
          wtfn1st=0.
!          wsfn0=0.

          write(cdtg,'(i12)') idtg1
!          read(cdtg,'(i4,i2,i8)')iyyyy,mm,ddhhmn


  11     format(i3.3,'TFN','0000',a4,a12)  ! ???TFM
  12     format(i3.3,'SFN','0000',a4,a12)  ! ???SFM

!          flag=.false.
!          if(myrank.eq.0) then
!            call dmsmsg("ALL",istat)
!            print *,'ready to open ifilin_ocaf'
!            call dmsopn(ifilin_ocaf,"r",istat)
!            flag=.true.
!          endif
!          call mpe_broadcast(istat,1,flag,mpe_integer)
!          if(istat.ne.0)then
!            if(myrank.eq.0) print *,'dmsopn ocaf error'
!            stop
!            call mpe_finalize
!            call dmsexit(-1)
!          endif


!          do k=1,lkvl+2
          do k=1,1
!           if(myrank.eq.0) call dmsopn(ifilin_sst,"r",istat)
!           call mpe_broadcast(istat,1,flag,mpe_integer)
!           if(istat.ne.0)then
!             if(myrank.eq.0) print *,'dmsopn sst error'
!               call mpe_finalize
!               call dmsexit(-1)
!           endif


!            write(cdtg,'(i12)') idtg1
!            read(cdtg,'(i4,i8)')iyyyy,mmddhhmn

!              iy=idtg_sst/1000000

           lncrec=nx*my
           write(lrec,11) k-1,ggdef,cdtg
           if(myrank .eq. 0) print *, 'lrec11=',lrec
           call dmsread(nx,my,lrec,lncrec,'H',ifilin,temp1(:,:),istat)
!           write(lrec,12) k-1,ggdef,cdtg
!           if(myrank .eq. 0) print *, 'lrec12=',lrec
!           call dmsread(nx,my,lrec,lncrec,'H',ifilin,temp2(:,:),istat)
!           if( myrank .eq. 72) then
!             print *,"ocaf 1: wtfn(914,265)=",temp1(914,265)
!             print *,"ocaf 2: wsfn(914,265)=",temp2(914,265)
!           endif

!           if( lreduce.eq.1 ) call reducepick(temp1(1,1),nxdef,nx,my)
!           if( lreduce.eq.1 ) call reducepick(temp2(1,1),nxdef,nx,my)


           do jj=1,jlistnum
              j=jlist1(jj)
              nxj=nxdef_2d(j)
              if( lreduce.eq.1 )then
                call reducepick (temp1(1,j),nxdef(j),nx,1)
!                call reducepick (temp2(1,j),nxdef(j),nx,1)
              endif
              do ii = 1, nxj
!                  if(sitmask(i,jj).eq. 1)then
                 i=nxjstart(j)+ii-1
                 wtfn1st(i,k,jj)=temp1(i,j)
!                 wsfn0(i,k,jj)=temp2(i,j)
!                  endif
!                  if((myrank .EQ.72 ) .AND. (jj .EQ. 3) .AND. (i .EQ.
!                  913))then
!                    print
!                    *,"wtfn12(913,",k,",3,",mm,")=",wtfn12(i,k,jj,mm) &
!                           ,"wsfn12(913,",k,",3,",mm,")=",wsfn12(i,k,jj,mm)
!                  endif
!                  if((myrank .EQ.72 ) .AND. (jj .EQ. 3) .AND. (i .EQ.
!                  914))then
!                    print
!                    *,"wtfn12(914,",k,",3,",mm,")=",wtfn12(i,k,jj,mm) &
!                           ,"wsfn12(914,",k,",3,",mm,")=",wsfn12(i,k,jj,mm)
!                  endif
              enddo   !end of ii
            enddo    !end of jj
         
          enddo !end of k

          deallocate(temp1,temp2)

!         call dmscls(ifilin_ocaf,istat)

        END SUBROUTINE read_ocaf0


!--------------------------------------------------

        SUBROUTINE read_godas(idtg1)
! Ben-Jei Tsuang, NCHU, June 2009, Read GODAS data
! (http://http://cfs.ncep.noaa.gov/cfs/godas/)
! An initial ocean dataset is generated by combining
! script:  irish3::/tcrg/u40bjt00/get/get_godas.sh
!
! history:
!   2008/8/19: modified the code from mo_so4.f90
!   2009/10/23: modified the code from read_woa0
  
!          USE mo_control,       ONLY: ngl, nlon, ngodas, lamip, lmlo
!          USE mo_mpi,           ONLY: p_parallel_io, p_bcast, p_io, p_pe 
!          USE mo_exception,     ONLY: finish, message, message_text
          
!          USE mo_io,            ONLY: io_open_unit, io_close, io_read, &
!                                      io_var_id, io_file_id, io_open, &
!                                      woanc0, woanc1, woanc2
          USE mo_netCDF,        ONLY: io_inq_dimid, io_inq_dimlen,     &
                                      io_inq_varid, io_get_var_double, &
                                      io_get_vara_double, io_get_att_double
!          USE mo_decomposition, ONLY: lc => local_decomposition, &
!                                      gl_dc => global_decomposition
!          USE mo_transpose,     ONLY: scatter_gp
!          USE mo_filename,      ONLY: NETCDF
          
          
          !  Local scalars: 
          
          ! number of codes read from godas file
          INTEGER, PARAMETER :: nrec = 5
          CHARACTER (12) :: cname, cwoa(nrec)
          
          REAL, ALLOCATABLE, TARGET :: zin(:,:,:,:)
!          REAL, ALLOCATABLE, TARGET :: zintemp(:,:)
!          REAL, POINTER :: gl_woa(:,:,:,:)
          REAL      :: missing_value
          
          INTEGER               :: io_nlon  ! number of longitudes in NetCDF file
          INTEGER               :: io_ngl   ! number of latitudes in NetCDF file
          INTEGER               :: io_ndepth  ! number of odepths in NetCDF file
          INTEGER               :: io_ntime  ! number of timesteps in NetCDF file
          INTEGER, DIMENSION(4) :: io_start ! start index for NetCDF-read
          INTEGER, DIMENSION(4) :: io_count ! number of iterations for NetCDF-read
          
          INTEGER               :: jk ,i     ! loop index
          INTEGER               :: irec      ! variable index
          ! Read GODAS data file
          ! ===============
          INTEGER :: IO_file_id0, IO_file_id1, IO_file_id2
!ps          CHARACTER (12) :: fn0, fn1, fn2
          CHARACTER (80) :: fn0, fn1, fn2
          INTEGER*4       :: ihy0, ihy1, ihy2
          LOGICAL       :: lex0, lex1, lex2
          INTEGER       :: status
!ps          
          INTEGER*8:: idtg1
          INTEGER, PARAMETER :: nerr = 6     ! error output stream
          TYPE (FILE_INFO), save :: woanc0, woanc1, woanc2
          INTEGER, PARAMETER:: IO_READ=1
          LOGICAL :: lmlo=.FALSE.
          INTEGER ::  io_var_id     ! IO_file_id, IO_dim_id
          LOGICAL flag
          INTEGER :: j,jj,nxj,im,ii
          INTEGER :: lnc
!ps
          ! Read world ocean atlas data file
          ! ===============

          flag=.false.
          IF (p_parallel_io) THEN
            WRITE(nerr,'(/,A,I2)') ' Read GODAS 1.0 '
            IF(lamip .AND. .NOT. lmlo) THEN
              WRITE(nerr,*) 'This is an AMIP run with lgodas enable (lamip = .true. & lgodas = .true. ).'
              WRITE(nerr,*) 'Read data from NCEP Global Ocean Data Assimilation System (GODAS)'
              WRITE(nerr,*) '(http://www.cpc.ncep.noaa.gov/products/GODAS/)'
              WRITE(nerr,*) 'or from Ishii dataset'
              WRITE(nerr,*) '(http://dss.ucar.edu/datasets/ds285.3/docs/) '
              
!ps              CALL set_years(ihy0, ihy1, ihy2)
              ihy1=idtg1/100000000
              ihy0=ihy1-1
              ihy2=ihy1+1              
              call chlen (ifilin_nc,80,lnc)
              WRITE (fn0, '(A,A6,i4)') ifilin_nc(1:lnc),'/godas',ihy0
              WRITE (fn1, '(A,A6,i4)') ifilin_nc(1:lnc),'/godas',ihy1
              WRITE (fn2, '(A,A6,i4)') ifilin_nc(1:lnc),'/godas',ihy2
 
!ps              WRITE (fn0, '("godas",i4)') ihy0
!ps              WRITE (fn1, '("godas",i4)') ihy1
!ps              WRITE (fn2, '("godas",i4)') ihy2
              
              WRITE (nerr, '(/)')
            
              WRITE(message_text,*) 'fn0: ', TRIM(fn0),' fn1: ',TRIM(fn1), ' fn2: ',TRIM(fn2)
              WRITE(nerr,*) 'read_godas',message_text
              
              INQUIRE (file=trim(fn0), exist=lex0)
              INQUIRE (file=trim(fn1), exist=lex1)
              INQUIRE (file=trim(fn2), exist=lex2)

              IF ( .NOT. lex0 ) THEN
                WRITE (message_text,*) 'Could not open file <',fn0,'>'
                WRITE(nerr,*) '', message_text
                ! CALL finish ('read_godas', 'run terminated.')
              ELSE
                CALL IO_open (fn0, woanc0, IO_READ)
              END IF
              
              IF ( .NOT. lex1 ) THEN
                WRITE (message_text,*) 'Could not open file <',fn1,'>'
                WRITE(nerr,*) '',message_text
                ! CALL finish ('read_godas', 'run terminated.')
              ELSE
                CALL IO_open (fn1, woanc1, IO_READ)      
              END IF
              
              IF ( .NOT. lex2 ) THEN
                WRITE (message_text,*) 'Could not open file <',fn2,'>'
                WRITE(nerr,*) '',message_text
                ! CALL finish ('read_godas', 'run terminated.')
              ELSE
                CALL IO_open (fn2, woanc2, IO_READ)
              END IF
               
            ELSE
              WRITE(nerr,'(/,A,I2)') ' Read GODAS data from unit ', ngodas
              WRITE(nerr,'(/,A,I2)') ' (http://www.cpc.ncep.noaa.gov/products/GODAS/) '
          
              WRITE(nerr,*) 'This is no AMIP run (lamip = .false.)'
              INQUIRE (ngodas, exist=lex1)
              ! unit ngodas=98      
              WRITE(message_text,*) 'lex1: ', lex1
              print *,'read_godas',message_text
              IF (lex1) THEN
                woanc1%format = NETCDF
                CALL IO_open_unit (ngodas, woanc1, IO_READ)
                ! has to be fixed...
                !          CALL IO_read_header(sstnc1)
                !          CALL IO_info_print(sstnc1)
                WRITE(nerr,'(/,A,I2)') ' Read GODAS 2.0: open successfully ',woanc1%file_id
              ELSE
                WRITE (message_text,*) 'Could not open unit <',ngodas,'>'
                !WRITE(nerr,*) '',message_text
                !WRITE(nerr,*) 'read_godas', 'Could not open godas file'
              ENDIF
          
             
            END IF  !end IF(lamip .AND. .NOT. lmlo)
            
            !WRITE(nerr,'(/,A,I2)') ' Read GODAS 3.0 '
            IF (lex1) THEN
              ! Check resolution
              CALL io_inq_dimid  (woanc1%file_id, 'lat', io_var_id)
              CALL io_inq_dimlen (woanc1%file_id, io_var_id, io_ngl)
              CALL io_inq_dimid  (woanc1%file_id, 'lon', io_var_id)
              CALL io_inq_dimlen (woanc1%file_id, io_var_id, io_nlon)
              CALL io_inq_dimid  (woanc1%file_id, 'time', io_var_id)
              CALL io_inq_dimlen (woanc1%file_id, io_var_id, io_ntime)
              
              !WRITE(nerr,'(/,A,I2)') ' Read GODAS 4.0 '
              WRITE(nerr,*) 'number of world ocean atlas data = ' &
			    ,'io_ngl=',io_ngl,',io_nlon=',io_nlon &
			    ,'io_ntime=',io_ntime

              IF (io_ngl/=ngl) THEN
                 WRITE(nerr,*) 'read_godas: unexpected resolution ',io_nlon,io_ngl
                 WRITE(nerr,*) 'expected number of latitudes = ',ngl
                 WRITE(nerr,*) 'number of latitudes of world ocean atlas data = ',io_ngl
                 WRITE(nerr,*)  'read_godas','unexpected resolution'
                 stop
              END IF
              CALL io_inq_dimid  (woanc1%file_id, 'depth', io_var_id)
              CALL io_inq_dimlen (woanc1%file_id, io_var_id, nodepth)
              !WRITE(nerr,'(/,A,I2)') ' Read GODAS 4.1 '
              WRITE(nerr,'(/,A,I2)') 'number of odepths = ',nodepth
              flag=.true.
            ELSE
              nodepth=nodepth0
            ENDIF
            IF(lamip .AND. .NOT. lmlo) THEN
              IF (lex0) THEN  
                CALL io_inq_dimid  (woanc0%file_id, 'depth', io_var_id)
                CALL io_inq_dimlen (woanc0%file_id, io_var_id, io_ndepth)
                IF (io_ndepth/=nodepth) THEN
                   WRITE(nerr,*) 'read_godas: inconsistent number of odepths between godas files',io_ndepth
                   WRITE(nerr,*) 'expected number of odepths = ',nodepth
                   WRITE (message_text,*) 'Read nodepth error in file <',fn0,'>'
                   WRITE(nerr,*) '',message_text
                   WRITE(nerr,*) 'read_godas','unexpected resolution'
                   stop
                END IF
              ENDIF
              IF (lex2) THEN          
                CALL io_inq_dimid  (woanc2%file_id, 'depth', io_var_id)
                CALL io_inq_dimlen (woanc2%file_id, io_var_id, io_ndepth)
                IF (io_ndepth/=nodepth) THEN
                   WRITE(nerr,*) 'read_godas: inconsistent number of ', &
                                  'odepths between godas files',io_ndepth
                   WRITE(nerr,*) 'expected number of odepths = ',nodepth
                   WRITE (message_text,*) 'Read nodepth error in file <',fn2,'>'
                   WRITE(nerr,*) '',message_text
                   WRITE(nerr,*) 'read_godas','unexpected resolution'
                   stop
                ENDIF
              ENDIF
            ENDIF
          END IF  !end IF (p_parallel_io) 
          !! RETURN
!ps          CALL p_bcast (nodepth, p_io)
!          call mpe_broadcast(nodepth,1,flag,mpe_integer)
          call mpe_bcast(nodepth,1,0,mpe_integer)
          
!          WRITE(nerr,'(/,A,I2)') ' Read GODAS 5.0, nodepth= ',nodepth
          
          !     Allocate memory for ot12 and os12 per PE
          
          IF (.NOT. ALLOCATED(odepths)) ALLOCATE (odepths(nodepth))
          IF (.NOT. ALLOCATED(ot12)) ALLOCATE (ot12(nxp, nodepth, my_max, 0:13))
          IF (.NOT. ALLOCATED(os12)) ALLOCATE (os12(nxp, nodepth, my_max, 0:13))
          IF (.NOT. ALLOCATED(ou12)) ALLOCATE (ou12(nxp, nodepth, my_max, 0:13))
          IF (.NOT. ALLOCATED(ov12)) ALLOCATE (ov12(nxp, nodepth, my_max, 0:13))
          IF(lmixedlayer) THEN
            IF(.NOT. ALLOCATED(mixedlayer12)) ALLOCATE (mixedlayer12(nxp, my_max,0:13))
          ENDIF
          
          ! Read odepths
          flag=.false.
          IF (p_parallel_io) THEN
            IF (lex1) THEN
              !WRITE(nerr,'(/,A,I2)') ' Read GODAS 6.0 '
              CALL io_inq_dimid  (woanc1%file_id, 'depth', io_var_id)
              CALL io_get_var_double (woanc1%file_id, io_var_id, odepths)
              WRITE(nerr,*) 'number of odepths = ',odepths
              flag=.true.
            else
              odepths=odepth0
            endif
          ENDIF
!ps          CALL p_bcast (odepths(1:nodepth), p_io)
!          call mpe_broadcast(odepths,nodepth,flag,mpe_double)
          call mpe_bcast(odepths,nodepth,0,mpe_double)
          IF ( lwarning_msg.GE.3 ) THEN       
              WRITE (nerr,*) 'read_godas:: pe=',p_pe,',  odepths=',odepths(1:nodepth)
!!          print *,'pe=',p_pe,',  odepths=',odepths(1:nodepth)
          ENDIF
          
          ! Codes read from GODAS
          
          cwoa(1) = 'ot'    ! ocean temperature profile (K)
          cwoa(2) = 'os'    ! ocean salinity profile (0/00)
          cwoa(3) = 'ou'    ! ocean u-component current (m/s)
          cwoa(4) = 'ov'    ! ocean v-component current (m/s)
          cwoa(5) = 'mixedlayer'    ! ocean mixed layer (m)
          DO irec = 1, nrec
            IF (.NOT. ALLOCATED(zin)) ALLOCATE (zin(nlon,nodepth,ngl,0:13))
            IF (p_parallel_io) THEN
            !WRITE(nerr,'(/,A,I2)') ' Read GODAS 7.0 '
            !     Allocate memory for godas global fields
!            IF (.NOT. ALLOCATED(zin)) ALLOCATE (zin(nlon,nodepth,ngl,0:13))
              zin=xmissing
              cname = cwoa(irec)
              ! read world ocean atlas data
              IF (lex1) THEN
                status = NF_INQ_VARID (woanc1%file_id, cname, io_var_id)
                IF (status /= NF_NOERR) THEN
                  WRITE(nerr,*) 'IO_INQ_VARID :', woanc1%file_id, cname
                  WRITE(nerr,*) 'IO_INQ_VARID', NF_STRERROR(status)
                ELSE
                  IF (irec == 3) lou=.TRUE.
                  IF (irec == 4) lov=.TRUE.
                  CALL io_get_att_double (woanc1%file_id, io_var_id, &
                                            '_FillValue', missing_value)
!!          WRITE(nerr,*) 'pe=',p_pe,', read GODAS 7.1 _FillValue= ',missing_value
                  IF(irec .eq. 5) then
                      io_start(:) = (/       1,   1,  1,        1 /)
                      io_count(:) = (/ io_nlon, ngl,  1, io_ntime /)
                    ! for depth jk: read io_nlon longitudes, ngl latitudes and 12 months
                      CALL io_get_vara_double(woanc1%file_id,io_var_id,io_start,io_count, &
                                         zin(1:io_nlon,jk,:,1:io_ntime))
                      DO im=1,12
                        print *,'GODAS: irec=',irec,',jk=',jk, &
                                   ',zin(212,',jk,',235,',im,')=',zin(212,jk,ngl+1-235,im)
                      ENDDO
                  ELSE
                    DO jk = 1, nodepth
                      io_start(:) = (/       1,   1, jk,        1 /)
                      io_count(:) = (/ io_nlon, ngl,  1, io_ntime /)
                    ! for depth jk: read io_nlon longitudes, ngl latitudes and 12 months
                      CALL io_get_vara_double(woanc1%file_id,io_var_id,io_start,io_count, &
                                         zin(1:io_nlon,jk,:,1:io_ntime))
                      DO im=1,12
                        print *,'GODAS: irec=',irec,',jk=',jk, &
                                   ',zin(212,',jk,',235,',im,')=',zin(212,jk,ngl+1-235,im)
                      ENDDO
                    END DO
                    !!          WRITE(nerr,'(/,A,I2)') ' Read GODAS 7.2 '
                  ENDIF
                ENDIF
              ENDIF  !end lex1
  
  
              IF(lamip .AND. .NOT. lmlo) THEN
                IF (lex0) THEN
                  status = NF_INQ_VARID (woanc0%file_id, cname, io_var_id)
                  IF (status /= NF_NOERR) THEN
                    WRITE(nerr,*) 'IO_INQ_VARID :', woanc0%file_id, cname
                    WRITE(nerr,*) 'IO_INQ_VARID', NF_STRERROR(status)
                  ELSE
                    ! read world ocean atlas data december of last year
                    CALL io_inq_varid (woanc0%file_id, cname, io_var_id)
                    IF(irec .eq. 5) then
                      io_start(:) = (/       1,   1,  1, 12 /)
                      io_count(:) = (/ io_nlon, ngl,  1, 1 /)
                      ! for depth jk: read io_nlon longitudes, ngl latitudes and 1 months
                      CALL io_get_vara_double(woanc0%file_id,io_var_id,io_start,io_count, &
                                             zin(1:io_nlon,jk,:,0))
                      print *,'GODAS: irec=',irec,',jk=',jk, &
                                      ',zin(212,',jk,',235,0)=',zin(212,jk,ngl+1-235,0)
                    ELSE
                      DO jk = 1, nodepth
                        io_start(:) = (/       1,   1, jk, 12 /)
                        io_count(:) = (/ io_nlon, ngl,  1,  1 /)
                      ! for depth jk: read io_nlon longitudes, ngl latitudes and 1 months
                        CALL io_get_vara_double(woanc0%file_id,io_var_id,io_start,io_count, &
                                             zin(1:io_nlon,jk,:,0))
                        print *,'GODAS: irec=',irec,',jk=',jk, &
                                      ',zin(212,',jk,',235,0)=',zin(212,jk,ngl+1-235,0)
                      ENDDO
                    ENDIF
                  ENDIF
                ENDIF  !end lex0

                IF (lex2) THEN
                  status = NF_INQ_VARID (woanc2%file_id, cname, io_var_id)
                  IF (status /= NF_NOERR) THEN
                    WRITE(nerr,*) 'IO_INQ_VARID :', woanc2%file_id, cname
                    WRITE(nerr,*) 'IO_INQ_VARID', NF_STRERROR(status)
                  ELSE
                    ! read world ocean atlas data january of next year
                    CALL io_inq_varid (woanc2%file_id, cname, io_var_id)
                    IF(irec .eq. 5) then
                      io_start(:) = (/       1,   1,  1, 1 /)
                      io_count(:) = (/ io_nlon, ngl,  1, 1 /)
                      ! for depth jk: read io_nlon longitudes, ngl latitudes and 1 months
                        CALL io_get_vara_double(woanc2%file_id,io_var_id,io_start,io_count, &
                                             zin(1:io_nlon,jk,:,13))
                        print *,'GODAS: irec=',irec,',jk=',jk, &
                                     ',zin(212,',jk,',235,13)=',zin(212,jk,ngl+1-235,13)
                    ELSE
                      DO jk = 1, nodepth
                        io_start(:) = (/       1,   1, jk,  1 /)
                        io_count(:) = (/ io_nlon, ngl,  1,  1 /)
                      ! for depth jk: read io_nlon longitudes, ngl latitudes and 1 months
                        CALL io_get_vara_double(woanc2%file_id,io_var_id,io_start,io_count, &
                                             zin(1:io_nlon,jk,:,13))
                        print *,'GODAS: irec=',irec,',jk=',jk, &
                                     ',zin(212,',jk,',235,13)=',zin(212,jk,ngl+1-235,13)
                      ENDDO
                    ENDIF
                  ENDIF
                ENDIF  !end lex2
              ELSE            
                ! copy December to month 0
                zin(:,:,:,0)  = zin(:,:,:,12)
                ! copy January to month 13
                zin(:,:,:,13)  = zin(:,:,:,1)
              ENDIF  !end lamip&not lmlo

!ps                zin=MERGE(zin,xmissing,zin.NE.missing_value)
              zin=MERGE(zin,xmissing,zin.GT.-8.E+33)
              IF (lsitstart .AND. (.NOT.lwoa0)) THEN
                WRITE(nerr,'(/,A,5I5)') &
                     ' Read GODAS 7.3: CALL fill_missing: ',io_nlon,ngl,nodepth
!ps                CALL print_cpu_time()
                DO im = 0, 13
                  CALL fill_missing2(zin(1:io_nlon,:,:,im),io_nlon,ngl,1,.FALSE.)
                  DO jk = 1, nodepth
                  print *,'GODAS after fillmissing: irec=',irec,',jk=',jk, &
                      ',zin(212,',jk,',235,',im,')=',zin(212,jk,ngl+1-235,im)
                  ENDDO
                ENDDO
                WRITE(nerr,'(/,A,I2)') &
                      ' Read GODAS 7.4: After CALL fill_missing '
!ps                CALL print_cpu_time()
              ENDIF
            ENDIF  !end p_parallel_io
              !WRITE(nerr,'(/,A,I2)') ' Read GODAS 8.0 '

          
!ps            NULLIFY (gl_woa)
!ps            DO i = 0, 13
!ps              IF (p_pe == p_io) gl_woa => zin(:,:,:,i:i)
!ps              IF (irec == 1) THEN
!ps                CALL scatter_gp(gl_woa, ot12(:,:,:,i:i), gl_dc)
!ps              ELSE IF (irec == 2) THEN
!ps                CALL scatter_gp(gl_woa, os12(:,:,:,i:i), gl_dc)
!ps              ELSE IF (irec == 3) THEN
!ps                CALL scatter_gp(gl_woa, ou12(:,:,:,i:i), gl_dc)
!ps              ELSE IF (irec == 4) THEN
!ps                CALL scatter_gp(gl_woa, ov12(:,:,:,i:i), gl_dc)
!ps              END IF
!ps            END DO
              
            DO jk = 1, nodepth
              DO im=0, 13
                flag=.false.
                if(myrank .eq. 0) flag=.true.
!                call mpe_broadcast(zin(:,jk,:,im),nx*my,flag,mpe_double)
                call mpe_bcast(zin(:,jk,:,im),nx*my,0,mpe_double)
                if (myrank .eq. 44) then
                  print *,'GODAS after broadcast: jk=',jk,',im=',im, & 
                      ',zin(212,',jk,',235,',im,')=',zin(212,jk,ngl+1-235,im)
                endif
!                if( lreduce.eq.1 ) call reducepick(zin(1,jk,1,im),nxdef,nx,my)
    
                DO jj = 1, jlistnum
                  j=jlist1(jj)
                  nxj=nxdef_2d(j)
                  if( lreduce.eq.1 )call reducepick(zin(1,jk,j,im),nxdef(j),nx,1)
                  DO ii=1,nxj
                    i=nxjstart(j)+ii-1
                    IF(irec .eq. 1) THEN
                      ot12(ii,jk,jj,im) = zin(i,jk,ngl-j+1,im)
!                      if((myrank.eq.44).AND.(i.eq.212).AND.(jj.eq.3).AND.(j.eq.235)) then
!                        print *,'ot12(',i,',',jk,',',jj,',',im,')=',ot12(i,jk,jj,im) 
!                      endif
                    ELSE IF(irec .eq. 2) THEN
                      os12(ii,jk,jj,im) = zin(i,jk,ngl-j+1,im)
                    ELSE IF(irec .eq. 3) THEN
                      ou12(ii,jk,jj,im) = zin(i,jk,ngl-j+1,im)
                    ELSE IF(irec .eq. 4) THEN
                      ov12(ii,jk,jj,im) = zin(i,jk,ngl-j+1,im)
                    ELSE IF(irec .eq. 5) THEN
                      mixedlayer12(ii,jj,im) = zin(i,1,ngl-j+1,im)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
            ENDDO

          END DO  !end do irec=1, nrec
  
!ps          CALL p_bcast (lou, p_io)    
!ps          CALL p_bcast (lov, p_io)    
          IF (lwarning_msg.GE.3) THEN    
            WRITE(nerr,*) 'pe=',p_pe,   &
             ', read GODAS 9.1, os12(212,0,3,0:13)= ',os12(212,0,3,0:13)
          END IF
          
          IF (p_parallel_io) THEN
            IF (lex1) CALL io_close (woanc1)
            IF(lamip .AND. .NOT. lmlo) THEN
              IF (lex0) CALL io_close (woanc0)
              IF (lex2) CALL io_close (woanc2)       
            END IF
!            DEALLOCATE (zin)
            WRITE(nerr,*) 'read GODAS'
          END IF
          
!          DEALLOCATE (zintemp)
          DEALLOCATE (zin)
        END SUBROUTINE read_godas


!------------------------------------------------------------------------------		
!  woa0
! ----------------------------------------------------------------------
!
! Ben-Jei Tsuang, NCHU, AUG 2008, Read WORLD OCEAN ATLAS 2005 data
! (woa05)
! (http://www.nodc.noaa.gov/OC5/WOA05/pr_woa05.html)
! Change woa ocean temperature variable from "t0112an1" to "ot".
! Change woa ocean salinity variable from "s0112an1" to "os"
! Change the unit of ocean temperature of woa05 from degree C to
! Kelvin
! ######### script >>
! cdo -f nc chname,t0112an1,ot t0112an1.nc xx
! cdo -f nc addc,273.15 xx ot.nc
! cdo -f nc chname,s0112an1,os s0112an1.nc os.nc
! cdo -f nc merge ot.nc os.nc woa05_clim.nc
! cdo -f nc interpolate,t21grid woa05_clim.nc T21_woa05_clim.nc
! cdo -f nc interpolate,t31grid woa05_clim.nc T31_woa05_clim.nc
! cdo -f nc interpolate,t42grid woa05_clim.nc T42_woa05_clim.nc
! cdo -f nc interpolate,t63grid woa05_clim.nc T63_woa05_clim.nc
! cdo -f nc interpolate,t85grid woa05_clim.nc T85_woa05_clim.nc
! cdo -f nc interpolate,t106grid woa05_clim.nc T106_woa05_clim.nc
!
! cp -pr T21_woa05_clim.nc /u1/u40bjt00/MPI/T21/amip2/
! cp -pr T31_woa05_clim.nc /u1/u40bjt00/MPI/T31/amip2/
! cp -pr T42_woa05_clim.nc /u1/u40bjt00/MPI/T42/amip2/
! cp -pr T63_woa05_clim.nc /u1/u40bjt00/MPI/T63/amip2/
! cp -pr T85_woa05_clim.nc /u1/u40bjt00/MPI/T85/amip2/
! cp -pr T106_woa05_clim.nc /u1/u40bjt00/MPI/T106/amip2/
! ######### << script
!
      SUBROUTINE read_woa0
! ----------------------------------------------------------------------
! Ben-Jei Tsuang, NCHU, June 2009, Read background initial WORLD OCEAN
! ATLAS 2005 data (woa05)
! (http://www.nodc.noaa.gov/OC5/WOA05/pr_woa05.html)
! An initial ocean dataset is generated by combining
!   ishii (from surface to 700 m depth), woa monthy data (form > 700
!   to
!   1500 m depth), and woa annual data (form > 1500 to 5500 m depth).
! Ishii data: http://dss.ucar.edu/datasets/ds285.3/docs/
! script:  irish3::/tcrg/u40bjt00/data/woa05/interp_echam.sh
!
! history:
!   2009/6/10: modified the code from read_woa0
!
! ----------------------------------------------------------------------
!    USE mo_control,       ONLY: ngl, nlon, nwoa0, lamip, lmlo
!    USE mo_mpi,           ONLY: p_parallel_io, p_bcast, p_io, p_pe,p_parallel,p_all_comm,p_barrier
!    USE mo_exception,     ONLY: finish, message, message_text

!    USE mo_io,            ONLY: io_open_unit, io_close, io_read, &
!                                io_var_id, io_file_id, io_open, &
!                                woa0nc1
     USE mo_netCDF,        ONLY: io_inq_dimid, io_inq_dimlen,     &
                                io_inq_varid, io_get_var_double, &
                                io_get_vara_double, io_get_att_double
!    USE mo_decomposition, ONLY: lc => local_decomposition, &
!                                gl_dc => global_decomposition
!    USE mo_transpose,     ONLY: scatter_gp
!    USE mo_filename,      ONLY: NETCDF


!  Local scalars:

! number of codes read from WOA file
      INTEGER, PARAMETER :: nrec = 5
      CHARACTER (12) :: cname, cwoa0(nrec)

      REAL, ALLOCATABLE, TARGET :: zin(:,:,:)
      REAL, POINTER :: gl_woa0(:,:,:)
      REAL      :: missing_value

      INTEGER               :: io_nlon  ! number of longitudes in NetCDF file
      INTEGER               :: io_ngl   ! number of latitudes in NetCDF file
      INTEGER               :: io_ndepth  ! number of odepth0 in NetCDF file
!!!    INTEGER, DIMENSION(3) :: io_start ! start index for NetCDF-read
!!!    INTEGER, DIMENSION(3) :: io_count ! number of iterations for
!NetCDF-read
!!!    INTEGER, DIMENSION(4) :: io_start ! start index for NetCDF-read
!!!    INTEGER, DIMENSION(4) :: io_count ! number of iterations for
!NetCDF-read
      INTEGER, ALLOCATABLE, TARGET :: io_start(:) ! start index for NetCDF-read
      INTEGER, ALLOCATABLE, TARGET :: io_count(:) ! number of iterations for NetCDF-read

      INTEGER               :: jk ,i     ! loop index
      INTEGER               :: irec      ! variable index
     ! Read world ocean atlas data file
     ! ===============
      INTEGER :: IO_file_id0, IO_file_id1, IO_file_id2
!      CHARACTER (12) :: fn0, fn1, fn2
      CHARACTER (80) :: fn0, fn1, fn2
      LOGICAL       :: lex1
      INTEGER       :: status, ncid, ndims, nvars, ngatts, unlimdimid

!ps
      INTEGER, PARAMETER :: nerr = 6     ! error output stream
      TYPE (FILE_INFO), save :: woa0nc1 
      INTEGER, PARAMETER:: IO_READ=1
      INTEGER :: io_var_id         !,IO_file_id,  IO_dim_id
      LOGICAL flag
      INTEGER :: j,jj,nxj,im,ii
      INTEGER :: lnc
!ps
     ! Read world ocean atlas data file
     ! ===============
!!    WRITE(nerr,'(/,A,I2)') ' Read WOA0 1.0 '
      flag=.false.
      IF (p_parallel_io) THEN
        WRITE(nerr,'(/,A,I2)') ' Read WOA0 from unit ', nwoa0
        WRITE(nerr,'(A,I2)') ' Read initial ocean temp. and sal. profiles from'
        WRITE(nerr,'(A,I2)') ' WORLD OCEAN ATLAS 2005 data (http://www.nodc.noaa.gov/OC5/WOA05/pr_woa05.html)'

        INQUIRE (nwoa0, exist=lex1)
        ! unit nwoa0=97
        IF (lex1) THEN
          woa0nc1%format = NETCDF
          CALL IO_open_unit (nwoa0, woa0nc1, IO_READ)
        ! has to be fixed...
        !          CALL IO_read_header(sstnc1)
        !          CALL IO_info_print(sstnc1)
          WRITE(nerr,'(/,A,I2)') ' Read WOA0 2.0: open successfully!',woa0nc1%file_id
        ELSE
!       CALL finish ('read_WOA0', 'Could not open WOA0 file')
          WRITE (message_text,*) 'Could not open unit<',nwoa0,'>'
          WRITE(nerr,*) '',message_text
          WRITE(nerr,*) 'read_woa0', 'Could not open woa0 file'
          stop
        ENDIF
        WRITE(nerr,'(/,A,I2)') ' Read WOA0 3.0 '

      ! Check resolution
        CALL io_inq_dimid  (woa0nc1%file_id, 'lat', io_var_id)
        CALL io_inq_dimlen (woa0nc1%file_id, io_var_id, io_ngl)
        CALL io_inq_dimid  (woa0nc1%file_id, 'lon', io_var_id)
        CALL io_inq_dimlen (woa0nc1%file_id, io_var_id, io_nlon)
!!      WRITE(nerr,'(/,A,I2)') ' Read WOA0 4.0 '
!!      WRITE(nerr,'(/,A,I2)') 'number of latitudes of world ocean atlas
!data = ',io_ngl
        IF (io_ngl/=ngl) THEN
          WRITE(nerr,*) 'read_WOA0: unexpected resolution',io_nlon,io_ngl
          WRITE(nerr,*) 'expected number of latitudes = ',ngl
         WRITE(nerr,*) 'number of latitudes of world ocean atlas data =',io_ngl
!         CALL finish ('read_WOA0','unexpected resolution')
          WRITE(nerr,*) 'read_WOA0','unexpected resolution'
          stop
        END IF
        CALL io_inq_dimid  (woa0nc1%file_id, 'depth', io_var_id)
        CALL io_inq_dimlen (woa0nc1%file_id, io_var_id, nodepth0)
!!      WRITE(nerr,'(/,A,I2)') ' Read WOA0 4.1 '
!!      WRITE(nerr,'(/,A,I2)') 'number of odepth0 = ',nodepth0
        flag=.true.
      END IF
     !! RETURN
!ps     CALL p_bcast (nodepth0, p_io)
!      CALL mpe_broadcast(nodepth0,1,flag,mpe_integer)
      CALL mpe_bcast(nodepth0,1,0,mpe_integer)
!!    WRITE(nerr,'(/,A,I2)') ' Read WOA0 5.0, nodepth0= ',nodepth0

     !     Allocate memory for ot0, os0, ou0, ov0 per PE

      IF (.NOT. ALLOCATED(odepth0)) ALLOCATE (odepth0(nodepth0))
      IF (.NOT. ALLOCATED(ot0)) ALLOCATE (ot0(nxp, nodepth0, my_max))
      IF (.NOT. ALLOCATED(os0)) ALLOCATE (os0(nxp, nodepth0, my_max))
      IF (.NOT. ALLOCATED(ou0)) ALLOCATE (ou0(nxp, nodepth0, my_max))
      IF (.NOT. ALLOCATED(ov0)) ALLOCATE (ov0(nxp, nodepth0, my_max))
      IF (.NOT. ALLOCATED(mixedlayer0)) ALLOCATE (mixedlayer0(nxp, my_max))

     ! Read odepth0

      IF (p_parallel_io) THEN
        CALL io_inq_dimid  (woa0nc1%file_id, 'depth', io_var_id)
        CALL io_get_var_double (woa0nc1%file_id, io_var_id, odepth0)
!!      WRITE(nerr,'(/,A,I2)') ' Read WOA0 6.0 '
      ENDIF
!ps      CALL p_bcast (odepth0(1:nodepth0), p_io)
!      CALL mpe_broadcast(odepth0,nodepth0,flag,mpe_double)
      CALL mpe_bcast(odepth0,nodepth0,0,mpe_double)

      IF ( lwarning_msg.GE.3 ) THEN
        WRITE (nerr,*) 'read_WOA0:: pe=',p_pe,', odepth0=',odepth0(1:nodepth0)
!!    print *,'pe=',p_pe,',  odepth0=',odepth0(1:nodepth0)
      ENDIF

     ! Codes read from WORLD OCEAN ATLAS 2005

      cwoa0( 1) = 'ot'    ! ocean temperature profile (K)
      cwoa0( 2) = 'os'    ! ocean salinity profile (0/00)
      cwoa0( 3) = 'ou'    ! ocean u-component current (m/s)
      cwoa0( 4) = 'ov'    ! ocean v-component current (m/s)
      cwoa0( 5) = 'mixedlayer'    ! ocean mixed layer (m)
      DO irec = 1, nrec
        IF (.NOT. ALLOCATED(zin)) ALLOCATE (zin(nlon,nodepth0,ngl))
      IF (p_parallel_io) THEN
!!      WRITE(nerr,'(/,A,I2)') ' Read WOA0 7.0 '
      !     Allocate memory for WOA0 global fields
!ps        IF (.NOT. ALLOCATED(zin)) ALLOCATE (zin(nlon,nodepth0,ngl))
        cname = cwoa0(irec)
        ! read world ocean atlas data
!!!
!!!        CALL io_inq_varid (woa0nc1%file_id, cname, io_var_id)
!!!
        status = NF_INQ_VARID (woa0nc1%file_id, cname, io_var_id)
        IF (status /= NF_NOERR) THEN
          WRITE(nerr,*) 'IO_INQ_VARID :', woa0nc1%file_id, cname
!ps          CALL message ('IO_INQ_VARID', NF_STRERROR(status))
          WRITE(nerr,*) 'IO_INQ_VARID', NF_STRERROR(status)
          zin=xmissing
        ELSE
          IF (irec == 3) lou=.TRUE.
          IF (irec == 4) lov=.TRUE.
          IF (irec == 5) lmixedlayer=.TRUE.
          CALL io_get_att_double (woa0nc1%file_id, io_var_id,'_FillValue', missing_value)
!!          WRITE(nerr,*) 'pe=',p_pe,', read WOA0 7.1 _FillValue=
!',missing_value
          status = NF_INQ_VARNDIMS (woa0nc1%file_id, io_var_id, ndims)
          IF (status /= NF_NOERR) THEN
            WRITE(nerr,*) 'NF_INQ_VARNDIMS :', woa0nc1%file_id, nwoa0,cname
!            CALL message ('NF_INQ_VARNDIMS', NF_STRERROR(status))
            WRITE(nerr,*) 'NF_INQ_VARNDIMS', NF_STRERROR(status)
!            CALL finish  ('NF_INQ_VARNDIMS', 'Run terminated.')
            stop
          ELSE
            IF ((ndims.GE.3).AND.(ndims.LE.4)) THEN
              IF (.NOT. ALLOCATED(io_start)) ALLOCATE (io_start(ndims))
              IF (.NOT. ALLOCATED(io_count)) ALLOCATE (io_count(ndims))
            ELSE
              WRITE(nerr,*) 'NF_INQ_VARNDIMS :', woa0nc1%file_id, nwoa0,cname
              WRITE(nerr,*) 'ndims=', ndims
              WRITE(nerr,*) 'ndims should be within 3 -4'
!              CALL finish  ('NF_INQ_VARNDIMS', 'Run terminated.')
              stop
            ENDIF
          ENDIF

          DO jk = 1, nodepth0
            IF (ndims  == 3) THEN
              io_start(:) = (/       1,   1, jk /)
              io_count(:) = (/ io_nlon, ngl,  1 /)
            ELSEIF (ndims  == 4) THEN
              io_start(:) = (/       1,  1, jk, 1 /)
              io_count(:) = (/ io_nlon, ngl ,  1, 1 /)
            ELSE
            ENDIF
            ! for depth jk: read io_nlon longitudes, ngl latitudes and
            ! 12 months
            IF(irec .EQ. 5) then
             IF(jk .EQ. 1)then
              CALL io_get_vara_double(woa0nc1%file_id,io_var_id,io_start,io_count,&
                 &                  zin(1:io_nlon,jk,:))
             ENDIF
            ELSE
              CALL io_get_vara_double(woa0nc1%file_id,io_var_id,io_start,io_count,&
                 &                  zin(1:io_nlon,jk,:))
            ENDIF
            print *,"woa0: irec=",irec,"zin:(768,",jk,",277)=",zin(768,jk,ngl+1-277)
          END DO
!!          WRITE(nerr,'(/,A,I2)') ' Read WOA0 7.2 '
          zin(1:io_nlon,:,:)=MERGE(zin(1:io_nlon,:,:),xmissing,zin(1:io_nlon,:,:).NE.missing_value)

          IF (lsitstart) THEN
            IF (lwarning_msg.GE.2) THEN
               WRITE(nerr,'(/,A,5I5)') &
                 ' Read WOA0 7.3: CALL fill_missing: ',io_nlon,ngl,nodepth0
            ENDIF
            !!! CALL print_cpu_time()
            CALL fill_missing2(zin(1:io_nlon,1:nodepth0,:),io_nlon,ngl,nodepth0,.TRUE.)
            IF (lwarning_msg.GE.2) THEN
              WRITE(nerr,'(/,A,5I5)') &
                 ' Read WOA0 7.4: After CALL fill_missing '
            ENDIF
            !!! CALL print_cpu_time()
          ENDIF
          DEALLOCATE (io_start)
          DEALLOCATE (io_count)
        END IF
      ENDIF
!!      WRITE(nerr,'(/,A,I2)') ' Read WOA0 8.0 '

!ps      NULLIFY (gl_woa0)
!ps      IF (p_pe == p_io) gl_woa0 => zin(:,:,:)
!ps      IF (irec == 1) THEN
!ps        CALL scatter_gp(gl_woa0, ot0(:,:,:), gl_dc)
!ps      ELSE IF (irec == 2) THEN
!ps        CALL scatter_gp(gl_woa0, os0(:,:,:), gl_dc)
!ps      ELSE IF (irec == 3) THEN
!ps        CALL scatter_gp(gl_woa0, ou0(:,:,:), gl_dc)
!ps      ELSE IF (irec == 4) THEN
!ps        CALL scatter_gp(gl_woa0, ov0(:,:,:), gl_dc)
!ps      END IF
        DO jk = 1, nodepth0
          flag=.false.
          if(myrank .eq. 0) flag=.true.
!          call mpe_broadcast(zin(:,jk,:),nx*my,flag,mpe_double)
          call mpe_bcast(zin(:,jk,:),nx*my,0,mpe_double)
!          if( lreduce.eq.1 ) call reducepick (zin(1,jk,1),nxdef,nx,my)
!          IF(myrank.eq. myrank_check) THEN
!            print *,"after broadcast,woa0: irec=",irec,  &
!                         ",zin:(212,",jk,",235)=",zin(212,jk,ngl+1-235)
!          ENDIF

          DO jj = 1, jlistnum
            j=jlist1(jj)
            nxj=nxdef_2d(j)
            if( lreduce.eq.1 )call reducepick (zin(1,jk,j),nxdef(j),nx,1)
            DO ii=1,nxj
              i=nxjstart(j)+ii-1
              IF (irec .eq. 1) THEN
                ot0(ii,jk,jj)=zin(i,jk,ngl-j+1)
!                IF((myrank.eq.44) .and. (jj.eq.3) .and. (ii.eq.1))THEN
                IF((i.eq.768) .and. (j.eq.492) )THEN
                  print *,"myrank=",myrank
                  print *,"WOA0:ot0(",ii,",",jk,",",jj,")=",ot0(ii,jk,jj)
                ENDIF
              ELSE IF(irec .eq. 2) THEN
                os0(ii,jk,jj)=zin(i,jk,ngl-j+1)
!                IF((myrank.eq.44) .and. (jj.eq.3) .and. (i.eq.1))THEN
                IF((i.eq.768) .and. (j.eq.492) )THEN
                  print *,"myrank=",myrank
                  print *,"WOA0:os0(",ii,",",jk,",",jj,")=",os0(ii,jk,jj)
                ENDIF
              ELSE IF(irec .eq. 3) THEN
                ou0(ii,jk,jj)=zin(i,jk,ngl-j+1)
!                IF((myrank.eq.44) .and. (jj.eq.3) .and. (i.eq.1))THEN
!                  print *,"WOA0:ou0(",ii,",",jk,",",jj,")=",ou0(ii,jk,jj)
!                ENDIF
              ELSE IF(irec .eq. 4) THEN
                ov0(ii,jk,jj)=zin(i,jk,ngl-j+1)
!                IF((myrank.eq.44) .and. (jj.eq.3) .and. (i.eq.1))THEN
!                  print *,"WOA0:ov0(",ii,",",jk,",",jj,")=",ov0(ii,jk,jj)
!                ENDIF
              ELSE IF(irec .eq. 5) THEN
                mixedlayer0(ii,jj)=zin(i,jk,ngl-j+1)
!                IF((myrank.eq.44) .and. (jj.eq.3) .and. (i.eq.1))THEN
                IF((i.eq.768) .and. (j.eq.492) )THEN
                  print *,"WOA0:mixedlayer0(",ii,",",jj,")=",mixedlayer0(ii,jj)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO             

               

      END DO  !end do irec=1, nrec
      IF (lwarning_msg.GE.3) THEN
!ps        CALL p_barrier(p_all_comm)
!ps      IF (p_parallel_io) WRITE(nerr,*) 'pe=',p_pe,', read WOA0 9.1, os0(1,:,1)= ',os0(1,:,1)
        IF (p_parallel_io) THEN

        END IF
      END IF
!ps      CALL p_bcast (lou, p_io)
!ps      CALL p_bcast (lov, p_io)
      IF (p_parallel_io) THEN
        CALL io_close (woa0nc1)
        DEALLOCATE (zin)
        WRITE(nerr,*) 'read WOA0'
      END IF

      END SUBROUTINE read_woa0



!----------------------------------------------------------------------
! ----------------------------------------------------------------------
      SUBROUTINE read_dailygodas(idtg1,tau,dtx)

    ! U. Schulzweida, MPI, March 2007
    ! Ben-Jei Tsuang, NCHU, Sep 2015

!ps       USE mo_doctor,        ONLY: nout
!ps       USE mo_control,       ONLY: nist
!ps       USE mo_exception,     ONLY: finish, message, message_text
!ps       USE mo_io,           ONLY: io_open_unit, io_close, io_open, &
!ps                                io_var_id, io_file_id, io_read, &
!ps                                gpnc0, gpnc1, gpnc2
!ps       USE mo_mpi,           ONLY: p_parallel_io, p_bcast, p_io, p_pe
!ps       USE mo_transpose,     ONLY: scatter_gp
       USE mo_netcdf,        ONLY: io_open_unit, io_close, io_open, &
                                io_inq_dimid, io_inq_dimlen,     &
                                io_inq_varid, io_get_var_double, &
                                io_get_vara_double, io_get_att_double

   
       REAL      :: ydate
       INTEGER, PARAMETER:: LAST_RECORD=-999
       INTEGER, PARAMETER:: DAY_MINUS1=1
       INTEGER, PARAMETER:: THE_DATE=2
       INTEGER, PARAMETER:: DAY_PLUS1=3

       CHARACTER (90) :: fn(3)
       INTEGER       :: ihy0, ihy1, ihy2
       INTEGER       :: yrori, modyhrmn
       INTEGER       :: yr, mo, dy, hr, mn
!ps
       INTEGER*8 :: idtg1,idtg_temp
       REAL      :: tau,dtx,tauleft,hrleft
       INTEGER   :: icurrenttau
       INTEGER :: lnc
       TYPE (FILE_INFO), save :: gpnc0, gpnc1, gpnc2
       INTEGER :: io_var_id         !,IO_file_id,  IO_dim_id
       INTEGER, PARAMETER:: IO_READ=1
       LOGICAL :: flag=.false.
       character*12 cdtg
!ps

!ps       CALL get_date_components(next_date, yr, mo, dy, hr, mn, se)
!ps       ydate = yr*10000.+mo*100.+dy+(hr+mn/60.+se/3600.)/24.
       icurrenttau=int(tau)
       tauleft=float(int((tau-int(tau)+0.001)*3600./dtx))*dtx   !(sec)
       if(tauleft .eq. 3600.) then
        icurrenttau=icurrenttau+1
        tauleft=0.
       endif
       call dtgfix12(idtg1,idtg_temp,icurrenttau)
       write(cdtg,'(i12)')idtg1
       read(cdtg,'(i4,i8)')yrori,modyhrmn
       write(cdtg,'(i12)')idtg_temp
       read(cdtg,'(i4,i2,i2,i2,i2)')yr,mo,dy,hr,mn
       ydate=float(yr*10000+mo*100+dy)+float(hr)/24.+tauleft/3600./24.
       if(myrank .eq. 0) print *,"read_dailygodas,ydate=",ydate 
       
       IF (p_parallel_io) THEN
!ps         CALL set_years(ihy0, ihy1, ihy2)
         ihy1=idtg_temp/100000000
         ihy0=ihy1-1
         ihy2=ihy1+1
         CALL chlen (ifilin_nc,80,lnc)
         WRITE (fn(1), '(A,A11,i4)') ifilin_nc(1:lnc),'/dailygodas',ihy0
         WRITE (fn(2), '(A,A11,i4)') ifilin_nc(1:lnc),'/dailygodas',ihy1
         WRITE (fn(3), '(A,A11,i4)') ifilin_nc(1:lnc),'/dailygodas',ihy2
!ps         WRITE (fn(1), '("dailygodas",i4)') ihy0
!ps         WRITE (fn(2), '("dailygodas",i4)') ihy1
!ps         WRITE (fn(3), '("dailygodas",i4)') ihy2
!ps         WRITE (nerr, '(/)')
         WRITE(nerr,*) 'dailygodas ydate: ', ydate
         WRITE(nerr,*) 'dailygodas ydate: fn(1): ' &
           , TRIM(fn(1)),' fn(2): ', TRIM(fn(2)),' fn(3): ', TRIM(fn(3))
       ENDIF

!       IF (lresume .OR. lstart) THEN
       IF ((tau .eq. 0.) .OR. lsitstart) THEN
    !!! warm/cold start
         CALL read_godas_3days
       ELSEIF (ydate.LE.timevals_godas(2)) THEN
    ! data were read. Note that initial value of  timevals_godas=0.
         RETURN 
       ELSEIF (ydate.LE.timevals_godas(3)) THEN
    ! note that initial value of  timevals_godas=0.
    ! Shift left
        ot12(:,:,:,1)=ot12(:,:,:,2)
        os12(:,:,:,1)=os12(:,:,:,2)
        ou12(:,:,:,1)=ou12(:,:,:,2)
        ov12(:,:,:,1)=ov12(:,:,:,2)
        mixedlayer12(:,:,1)=mixedlayer12(:,:,2)
        timevals_godas(1)=timevals_godas(2)
        files_godas(1)=files_godas(2)
        nts_godas(1)=nts_godas(2)
        tsID_godas(1)=tsID_godas(2)
      
        ot12(:,:,:,2)=ot12(:,:,:,3)
        os12(:,:,:,2)=os12(:,:,:,3)
        ou12(:,:,:,2)=ou12(:,:,:,3)
        ov12(:,:,:,2)=ov12(:,:,:,3)
        mixedlayer12(:,:,2)=mixedlayer12(:,:,3)
        timevals_godas(2)=timevals_godas(3)
        files_godas(2)=files_godas(3)
        nts_godas(2)=nts_godas(3)
        tsID_godas(2)=tsID_godas(3)
        if( (files_godas(3).eq.3) .AND. (yrori.lt.yr)) then
           files_godas(2)=2
        endif 


        CALL read_godas_dayp1(DAY_PLUS1)    ! read day+1 data
      ELSE
        CALL read_godas_3days
      ENDIF

      CONTAINS
    !-----------------------------------
       SUBROUTINE read_godas_3days
     
       INTEGER               :: io_nlon  ! number of longitudes in NetCDF file
       INTEGER               :: io_ngl   ! number of latitudes in NetCDF file
       INTEGER               :: io_ndepth  ! number of odepths in NetCDF file
       REAL, ALLOCATABLE :: timevals1(:), timevals2(:)
     
       INTEGER       :: i, jk, k, j, m
       LOGICAL       :: lex1, lex2
       INTEGER       :: istat
       INTEGER       :: ndimid, nts1, tsID
       INTEGER       :: ndimid2, nts2
     
       REAL :: timevals_godas(3) = 0.  ! absoulte time (e.g., 19971003.25) GODAS PENTAD Data 

       istat=0
       flag=.false.
       IF (p_parallel_io) THEN
         INQUIRE (file=fn(2), exist=lex1)
         IF ( .NOT. lex1 ) THEN
           WRITE (message_text,*) 'Could not open file <',fn(2),'>'
!ps           CALL message('',message_text)
           WRITE (nerr,*) message_text
!ps           CALL finish ('read_dailygodas', 'run terminated.')
           istat=-1
           nodepth=40               ! number of depths of daily fodas pentad data (=40)
         ELSE         
           CALL IO_open (fn(2), gpnc1, IO_READ)
         !! WRITE(nerr,'(/,A,I2)') ' Read Pentad GODAS 3.0 '
           ! Check resolution
           CALL io_inq_dimid  (gpnc1%file_id, 'lat', io_var_id)
           CALL io_inq_dimlen (gpnc1%file_id, io_var_id, io_ngl)
           CALL io_inq_dimid  (gpnc1%file_id, 'lon', io_var_id)
           CALL io_inq_dimlen (gpnc1%file_id, io_var_id, io_nlon)
           CALL io_inq_dimid  (gpnc1%file_id, 'time', io_var_id)
           !!! CALL io_inq_dimlen (gpnc1%file_id, io_var_id, io_ntime)
           !!      WRITE(nerr,'(/,A,I2)') ' Read Pentad GODAS 4.0 '
           !!      WRITE(nerr,'(/,A,I2)') 'number of latitudes of world ocean atlas data = ',io_ngl
       
!ps           IF ( (io_ngl/=lc%nlat).OR.(io_nlon/=lc%nlon) ) THEN
           IF ( (io_ngl/=nlat).OR.(io_nlon/=nlon) ) THEN
              WRITE(nerr,*) 'read_godas_3days: unexpected resolution ',io_nlon,io_ngl
              WRITE(nerr,*) 'expected number of latitudes = ',nlat
              WRITE(nerr,*) 'number of latitudes of pentad GODAS data data = ',io_ngl
              WRITE(nerr,*) 'expected number of longitude = ',nlon
              WRITE(nerr,*) 'number of latitudes of pentad GODAS data = ',io_nlon
!ps              CALL finish ('read_dailygodas','unexpected resolution')
              WRITE(nerr,*) 'read_godas_3days','unexpected resolution'
              istat=-1
           ELSE
             CALL io_inq_dimid  (gpnc1%file_id, 'depth', io_var_id)
             CALL io_inq_dimlen (gpnc1%file_id, io_var_id, nodepth)
           !!      WRITE(nerr,'(/,A,I2)') ' Read Pentad GODAS 4.1 '
           !!      WRITE(nerr,'(/,A,I2)') 'number of odepths = ',nodepth
           ENDIF
         ENDIF
         flag=.true.
       ENDIF  !end if (p_parallel_io)

!       CALL mpe_broadcast (istat, 1, flag, mpe_integer)
       CALL mpe_bcast (istat, 1, 0, mpe_integer)
       if (myrank.eq.3) print *,'myrank3 istat=',istat
       if(istat .ne. 0)then
        if(myrank.eq.0) print *,'read_godas_3days fn2',' terminated.'
        stop
        call mpe_finalize
        call dmsexit(-1)
       endif
        
!       CALL p_bcast (nodepth, p_io)
!       CALL mpe_broadcast (nodepth, 1, flag, mpe_integer)
       CALL mpe_bcast (nodepth, 1, 0, mpe_integer)
       if(myrank .eq. 0) print *,"broadcast nodepth"     
   
       IF (.NOT. ALLOCATED(odepths)) ALLOCATE (odepths(nodepth))
       IF (.NOT. ALLOCATED(ot12)) ALLOCATE (ot12(nxp, nodepth, my_max,3))
       IF (.NOT. ALLOCATED(os12)) ALLOCATE (os12(nxp, nodepth, my_max,3))
       IF (.NOT. ALLOCATED(ou12)) ALLOCATE (ou12(nxp, nodepth, my_max,3))
       IF (.NOT. ALLOCATED(ov12)) ALLOCATE (ov12(nxp, nodepth, my_max,3))
       IF (.NOT. ALLOCATED(mixedlayer12)) ALLOCATE (mixedlayer12(nxp, my_max,3))
       
       ! Read odepths
       flag=.false.
       IF (p_parallel_io) THEN
         IF (lex1) THEN
           !!      WRITE(nerr,'(/,A,I2)') ' Read Pentad GODAS 6.0 '
           CALL io_inq_dimid  (gpnc1%file_id, 'depth', io_var_id)
           CALL io_get_var_double (gpnc1%file_id, io_var_id, odepths)
         else
           odepths=odepth0
         endif
         flag=.true.
       ENDIF
!       CALL p_bcast (odepths(1:nodepth), p_io)
!       CALL mpe_broadcast (odepths,nodepth, flag, mpe_double)
       CALL mpe_bcast (odepths,nodepth, 0, mpe_double)
       if(myrank .eq. 0) print *,"broadcast odepths",odepths(:)     
       IF ( lwarning_msg.GE.3 ) THEN       
           WRITE (nerr,*) 'read_godas:: pe=',p_pe,',  odepths=',odepths(1:nodepth)
!!       print *,'pe=',p_pe,',  odepths=',odepths(1:nodepth)
       ENDIF
       
       !     Read dailygodas-file
       istat=0
       flag=.false.
       IF (p_parallel_io) THEN
         CALL IO_INQ_DIMID (gpnc1%file_id, 'time', ndimid)
         CALL IO_INQ_DIMLEN (gpnc1%file_id, ndimid, nts1)
         CALL IO_INQ_VARID (gpnc1%file_id, 'time', io_var_id)
       
         IF ( nts1 .lt. 1 ) then
!ps           CALL finish('read_dailygodas', 'To few time steps < 1')
           WRITE(nerr,*) 'read_godas_3days', 'To few time steps < 1'
           istat=-1
         ELSE
           ALLOCATE (timevals1(nts1))
           CALL IO_GET_VAR_DOUBLE(gpnc1%file_id, io_var_id, timevals1)
           tsID=1
           DO WHILE ( (ydate.GT.INT(timevals1(tsID))).AND.(tsID.LE.nts1) )
             tsID=tsID+1
           ENDDO
           if(myrank .eq. 0) print*,"tsID=",tsID,",nts1=",nts1
           IF ( tsID .gt. nts1 ) THEN
             INQUIRE (file=fn(3), exist=lex2)
             IF ( .NOT. lex2 ) THEN
               WRITE (message_text,*) 'Could not open file <',fn(3),'>'
!ps             CALL message('',message_text)
!ps             CALL finish ('read_dailygodas', 'run terminated.')
               WRITE(nerr,*) message_text
               WRITE(nerr,*) 'read_godas_3days fn3 ', 'run terminated.'
               istat=-1
             ELSE
               CALL IO_open (fn(3), gpnc2, IO_READ)
               CALL IO_INQ_DIMID (gpnc2%file_id, 'time', ndimid2)
               CALL IO_INQ_DIMLEN (gpnc2%file_id, ndimid2, nts2)
               IF ( nts2 .lt. 1 ) THEN
                 WRITE (message_text,*) 'File <',fn(3),'>'
!ps             CALL message('',message_text)
!ps             CALL finish('read_dailygodas:', 'To few time steps < 1')
                 WRITE (nerr,*) message_text
                 WRITE (nerr,*) 'read_dailygodas:','To few time steps<1'
                 istat=-1
               ELSE
                 ALLOCATE (timevals2(nts2))
                 CALL IO_INQ_VARID (gpnc2%file_id, 'time', io_var_id)
                 CALL IO_GET_VAR_DOUBLE(gpnc2%file_id, io_var_id, timevals2)
                 tsID=1
                 DO WHILE ( (ydate.GT.INT(timevals2(tsID))).AND.(tsID.LE.nts2) )
                   tsID=tsID+1
                 ENDDO
                 IF ( tsID .gt. nts2 ) THEN
!ps             CALL finish('read_dailygodas', 'Date not found')
                   WRITE(nerr,*) 'read_dailygodas', 'Date not found'
                   istat=-1
                 ELSE
                   IF ( lwarning_msg.GE.1 ) THEN
                     WRITE (nerr,*) 'read_godas_3days: timevals2_date=',timevals2(tsID)
                   ENDIF  
                   tsID_godas(2)=tsID
                   files_godas(2)=3
                   nts_godas(2)=nts2
                 ENDIF
               ENDIF
               CALL IO_close(gpnc2)
               IF (ALLOCATED(timevals2)) DEALLOCATE(timevals2)
             ENDIF
           ELSE
             tsID_godas(2)=tsID
             files_godas(2)=2
             nts_godas(2)=nts1
             IF ( lwarning_msg.GE.1 ) THEN
               WRITE (nerr,*) 'read_godas_3days: timevals1_date=',timevals1(tsID)
             ENDIF
           ENDIF
         ENDIF
         CALL IO_close(gpnc1)
         IF (ALLOCATED(timevals1)) DEALLOCATE(timevals1)
         flag=.true.
       ENDIF

!       CALL mpe_broadcast (istat, 1, flag, mpe_integer)
       CALL mpe_bcast (istat, 1, 0, mpe_integer)
       if(istat .ne. 0)then
        if(myrank.eq.0) print *,'read_godas_3days.2', 'run terminated.'
        stop
        call mpe_finalize
        call dmsexit(-1)
       endif

       if(myrank .eq. 0) print*,'read to read_goads_dayp1' 
       !!! read day-1, day and day+1 data
       CALL read_godas_dayp1(THE_DATE)     ! day+0
       CALL read_godas_dayp1(DAY_MINUS1)   ! day-1
       CALL read_godas_dayp1(DAY_PLUS1)    ! day+1
       END SUBROUTINE read_godas_3days

    !-----------------------------------

       SUBROUTINE read_godas_dayp1(dayID)
       INTEGER, INTENT(in):: dayID               ! dayID = DAY_MINUS1, THE_DAY_MINUS1, or DAY_PLUS1
       REAL, ALLOCATABLE, TARGET :: zot(:,:,:)
       REAL, ALLOCATABLE, TARGET :: zos(:,:,:)
       REAL, ALLOCATABLE, TARGET :: zou(:,:,:)
       REAL, ALLOCATABLE, TARGET :: zov(:,:,:)
       REAL, ALLOCATABLE, TARGET :: zmixedlayer(:,:)
!!!       REAL(dp), ALLOCATABLE, TARGET :: zow(:,:,:)
!ps       REAL, POINTER :: gl_ot(:,:,:)
!ps       REAL, POINTER :: gl_os(:,:,:)
!ps       REAL, POINTER :: gl_ou(:,:,:)
!ps       REAL, POINTER :: gl_ov(:,:,:)
!!!       REAL(dp), POINTER :: gl_ow(:,:,:)
!ps
       INTEGER jk,jj,j,nxj,i,ii
       INTEGER istat


       IF(.NOT. ALLOCATED(zot)) ALLOCATE (zot(nx,nodepth,my))
       IF(.NOT. ALLOCATED(zos)) ALLOCATE (zos(nx,nodepth,my))
       IF(.NOT. ALLOCATED(zou)) ALLOCATE (zou(nx,nodepth,my))
       IF(.NOT. ALLOCATED(zov)) ALLOCATE (zov(nx,nodepth,my))
       IF(.NOT. ALLOCATED(zmixedlayer)) ALLOCATE (zmixedlayer(nx,my))

       flag=.false.     
       IF (p_parallel_io) THEN
!ps         IF(.NOT. ALLOCATED(zot)) ALLOCATE (zot(nx,nodepth,my))
!ps         IF(.NOT. ALLOCATED(zos)) ALLOCATE (zos(nx,nodepth,my))
!ps         IF(.NOT. ALLOCATED(zou)) ALLOCATE (zou(nx,nodepth,my))
!ps         IF(.NOT. ALLOCATED(zov)) ALLOCATE (zov(nx,nodepth,my))
!ps         IF(.NOT. ALLOCATED(zmixedlayer)) ALLOCATE (zmixedlayer(nx,my))
         IF (dayID .EQ. DAY_PLUS1) THEN
         ! read Day+1 data
           IF (tsID_godas(2).EQ.nts_godas(2)) THEN
             files_godas(3)=files_godas(2)+1
             tsID_godas(3)=1                   !first record of file 3
           ELSE
             files_godas(3)=files_godas(2)
             tsID_godas(3)=tsID_godas(2)+1
           ENDIF
         ELSEIF (dayID .EQ. DAY_MINUS1) THEN
         ! read Day-1 data
           IF (tsID_godas(2).EQ.1) THEN
             files_godas(1)=files_godas(2)-1
             tsID_godas(1)=LAST_RECORD                !last record of file 1 !!!????? NEED to CODE
           ELSE
             files_godas(1)=files_godas(2)
             tsID_godas(1)=tsID_godas(2)-1
           ENDIF
!!!         ELSEIF (dayID .EQ. THE_DATE) THEN
!!!         ! read Day+0 data
!!!           IF (tsID_godas(2).EQ.1) THEN
!!!             files_godas(2)=files_godas(2)-1
!!!             tsID_godas(2)=LAST_RECORD                !last record of file 1 !!!????? NEED to CODE
!!!           ELSE
!!!             files_godas(2)=files_godas(2)
!!!             tsID_godas(2)=tsID_godas(2)-1
!!!           ENDIF
         ENDIF
         CALL read_godas_1record (fn(files_godas(dayID)),nx,nodepth,my,tsID_godas(dayID), &
           nts_godas(dayID),timevals_godas(dayID),zot(:,:,:),zos(:,:,:),zou(:,:,:),zov(:,:,:),zmixedlayer(:,:),istat)
         print*,'end_read_godas_1record,dayID=',dayID,",",timevals_godas(dayID)
         flag=.true.
       ENDIF   !end if (p_parallel_io)

!       CALL mpe_broadcast (istat, 1, flag, mpe_integer)
       CALL mpe_bcast (istat, 1, 0, mpe_integer)
       if(istat .ne. 0)then
        if(myrank.eq.0) print *,'read_dailygodas', 'run terminated.'
        stop
        call mpe_finalize
        call dmsexit(-1)
       endif

!       CALL mpe_broadcast(timevals_godas,3,flag,mpe_double)
!       CALL mpe_broadcast(tsID_godas,3,flag,mpe_integer)
!       CALL mpe_broadcast(files_godas,3,flag,mpe_integer)
!       CALL mpe_broadcast(nts_godas,3,flag,mpe_integer)
       CALL mpe_bcast(timevals_godas,3,0,mpe_double)
       CALL mpe_bcast(tsID_godas,3,0,mpe_integer)
       CALL mpe_bcast(files_godas,3,0,mpe_integer)
       CALL mpe_bcast(nts_godas,3,0,mpe_integer)
       
!ps       NULLIFY (gl_ot)
!ps       NULLIFY (gl_os)
!ps       NULLIFY (gl_ou)
!ps       NULLIFY (gl_ov)
!!!       NULLIFY (gl_ow)
       !WRITE(nerr,*) 'read_dailygodas before scatter_gp nodepth=',nodepth,' p_pe=',p_pe,' p_io=',p_io
!ps       IF (p_pe == p_io) gl_ot => zot(:,:,:)
!ps       CALL scatter_gp (gl_ot, ot12(:,:,:,dayID), global_decomposition)
!ps       IF (p_pe == p_io) gl_os => zos(:,:,:)
!ps       CALL scatter_gp (gl_os, os12(:,:,:,dayID), global_decomposition)
!ps       IF (p_pe == p_io) gl_ou => zou(:,:,:)
!ps       CALL scatter_gp (gl_ou, ou12(:,:,:,dayID), global_decomposition)
!ps       IF (p_pe == p_io) gl_ov => zov(:,:,:)
!ps       CALL scatter_gp (gl_ov, ov12(:,:,:,dayID), global_decomposition)
!!!       IF (p_pe == p_io) gl_ow => zow(:,:,:)
!!!       CALL scatter_gp (gl_ow, ow12(:,:,:,dayID), global_decomposition)
       !WRITE(nerr,*) 'ot12(15,1:40,15,1)=',ot12(15,1:40,15,1)

       DO jk = 1, nodepth
         flag=.false.
         if(myrank .eq. 0) flag=.true.
!         call mpe_broadcast(zot(:,jk,:),nx*my,flag,mpe_double)
!         call mpe_broadcast(zos(:,jk,:),nx*my,flag,mpe_double)
!         call mpe_broadcast(zou(:,jk,:),nx*my,flag,mpe_double)
!         call mpe_broadcast(zov(:,jk,:),nx*my,flag,mpe_double)
         call mpe_bcast(zot(:,jk,:),nx*my,0,mpe_double)
         call mpe_bcast(zos(:,jk,:),nx*my,0,mpe_double)
         call mpe_bcast(zou(:,jk,:),nx*my,0,mpe_double)
         call mpe_bcast(zov(:,jk,:),nx*my,0,mpe_double)
         if(jk .eq. 1) then
!           call mpe_broadcast(zmixedlayer(:,:),nx*my,flag,mpe_double)
           call mpe_bcast(zmixedlayer(:,:),nx*my,0,mpe_double)
         endif
!         if( lreduce.eq.1 ) then
!           call reducepick(zot(1,jk,1),nxdef,nx,my)
!           call reducepick(zos(1,jk,1),nxdef,nx,my)
!           call reducepick(zou(1,jk,1),nxdef,nx,my)
!           call reducepick(zov(1,jk,1),nxdef,nx,my)
!           if(jk .eq. 1) then
!             call reducepick(zmixedlayer(1,1),nxdef,nx,my)
!           endif
!         endif
         DO jj = 1, jlistnum
           j=jlist1(jj)
           nxj=nxdef_2d(j)
           if( lreduce.eq.1 ) then
             call reducepick(zot(1,jk,j),nxdef(j),nx,1)
             call reducepick(zos(1,jk,j),nxdef(j),nx,1)
             call reducepick(zou(1,jk,j),nxdef(j),nx,1)
             call reducepick(zov(1,jk,j),nxdef(j),nx,1)
             if(jk .eq. 1) then
               call reducepick(zmixedlayer(1,j),nxdef(j),nx,1)
             endif
           endif
           DO ii=1,nxj
             i=nxjstart(j)+ii-1
             ot12(ii,jk,jj,dayID) = zot(i,jk,ngl-j+1)
!             if((myrank.eq.44).AND.(i.eq.212).AND.(jj.eq.3).AND.(j.eq.235)) then
!                print *,'ot12(',ii,',',jk,',',jj,',',dayID,')=',ot12(ii,jk,jj,dayID)
!             endif
             os12(ii,jk,jj,dayID) = zos(i,jk,ngl-j+1)
             ou12(ii,jk,jj,dayID) = zou(i,jk,ngl-j+1)
             ov12(ii,jk,jj,dayID) = zov(i,jk,ngl-j+1)
             if(jk .eq. 1) then
               mixedlayer12(ii,jj,dayID) = zmixedlayer(i,ngl-j+1)
             endif
           ENDDO
         ENDDO
       ENDDO


!ps       IF (p_parallel_io) THEN
         DEALLOCATE (zot)
         DEALLOCATE (zos)
         DEALLOCATE (zou)
         DEALLOCATE (zov)
         DEALLOCATE (zmixedlayer)
!!!         DEALLOCATE (zow)
!ps       ENDIF
       END SUBROUTINE read_godas_dayp1
    !-----------------------------------

       SUBROUTINE read_godas_1record(fname,nlon,nodepth,nlat,tsID,nts,recdate,zot,zos,zou,zov,zmixedlayer,istat)
     
!ps       CHARACTER (14), INTENT(in):: fname
       CHARACTER (*), INTENT(in):: fname
       INTEGER, INTENT(in):: nlon, nodepth, nlat
       INTEGER, INTENT(in out):: tsID
       INTEGER, INTENT(out):: nts
       INTEGER, INTENT(out):: istat
       REAL, INTENT(out):: recdate    ! record data in absolute time foremat, e.g., 20130911.1350
       REAL, DIMENSION(nlon,nodepth,nlat), INTENT(out):: zot, zos, zou, zov
       REAL, DIMENSION(nlon,nlat), INTENT(out):: zmixedlayer
     
       REAL      :: missing_value
       INTEGER               :: io_ntime  ! number of timesteps in NetCDF file
       INTEGER, DIMENSION(4) :: io_start ! start index for NetCDF-read
       INTEGER, DIMENSION(4) :: io_count ! number of iterations for NetCDF-read
       INTEGER, DIMENSION(4) :: imixed_start ! start index for NetCDF-read
       INTEGER, DIMENSION(4) :: imixed_count ! number of iterations for NetCDF-read
     
       INTEGER       :: i, jk, k, j, m
       LOGICAL       :: lex2
       !!! INTEGER       :: start(4), COUNT(4), nvarid, ndimid, nts, tsID
       INTEGER       :: otid2,osid2,ouid2,ovid2,ndimid2,mixedid2
       REAL, ALLOCATABLE :: timevals2(:)
     
     
       ! read one-record godas data
       istat=0
       IF (.NOT.p_parallel_io) RETURN    !!! Only for p_parallel_io, else return
       INQUIRE (file=fname, exist=lex2)
       IF ( .NOT. lex2 ) THEN
         WRITE (message_text,*) 'Could not open file <',fname,'>'
!ps         CALL message('',message_text)
!ps         CALL finish ('read_dailygodas', 'run terminated.')
         WRITE (nerr,*) message_text
         WRITE (nerr,*) 'read_dailygodas', 'run terminated.'
         istat=-1
       ELSE
         CALL IO_open (fname, gpnc2, IO_READ)
         CALL IO_INQ_DIMID (gpnc2%file_id, 'time', ndimid2)
         CALL IO_INQ_DIMLEN (gpnc2%file_id, ndimid2, nts)
         WRITE (nerr,*) 'read_godas_1record:1. nts=',nts
         IF ( nts .lt. 1 ) THEN
           WRITE (message_text,*) 'File <',fname,'>'
!ps         CALL message('',message_text)
!ps         CALL finish('read_godas_1record:', 'To few time steps < 1')
           WRITE(nerr,*) message_text
           WRITE(nerr,*) 'read_godas_1record:', 'To few time steps < 1'
           istat=-1
         ELSE
           ALLOCATE (timevals2(nts))
           CALL IO_INQ_VARID (gpnc2%file_id, 'time', io_var_id)
           CALL IO_GET_VAR_DOUBLE(gpnc2%file_id, io_var_id, timevals2)
           recdate=timevals2(tsID)
           
           WRITE (nerr,*) 'read_godas_1record: nts=',nts &
                        ,',timevals(tsID)=',timevals2(tsID)
           IF(tsID .lt. nts) THEN
             print *, 'read_godas_1record: timevals(tsID+1)=',timevals2(tsID+1)
           ENDIF
           DEALLOCATE (timevals2)
           CALL IO_INQ_VARID (gpnc2%file_id, 'ot', otid2)
           CALL IO_INQ_VARID (gpnc2%file_id, 'os', osid2)
           CALL IO_INQ_VARID (gpnc2%file_id, 'ou', ouid2)
           CALL IO_INQ_VARID (gpnc2%file_id, 'ov', ovid2)
           CALL IO_INQ_VARID (gpnc2%file_id, 'mixedlayer', mixedid2)
!!!       CALL IO_INQ_VARID (gpnc2%file_id, 'ow', owid2)
           CALL io_get_att_double (gpnc2%file_id, otid2, '_FillValue', missing_value)
           IF (tsID.EQ.LAST_RECORD) tsID=nts     !!! modify tsID for LAST_RECORD
           IF (tsID .gt. nts) THEN
         !!! data out range, set to be missing value
             zot(:,:,:)=xmissing
             zos(:,:,:)=xmissing
             zou(:,:,:)=xmissing
             zov(:,:,:)=xmissing
             zmixedlayer(:,:)=xmissing
             WRITE (nerr,*) 'read_dailygodas, date=',recdate &
                           ,'DATA out of range!!! Set to MISSING DATA!'
           ELSE
             DO jk=1, nodepth
               io_start(:) = (/ 1, 1, jk, tsID /)
               io_count(:) = (/ nlon, nlat, 1, 1 /)
               CALL IO_GET_VARA_DOUBLE (gpnc2%file_id,otid2,io_start,io_count, zot(1:nlon,jk,1:nlat))
               CALL IO_GET_VARA_DOUBLE (gpnc2%file_id,osid2,io_start,io_count, zos(1:nlon,jk,1:nlat))
               CALL IO_GET_VARA_DOUBLE (gpnc2%file_id,ouid2,io_start,io_count, zou(1:nlon,jk,1:nlat))
               CALL IO_GET_VARA_DOUBLE (gpnc2%file_id,ovid2,io_start,io_count, zov(1:nlon,jk,1:nlat))
               imixed_start(:) = (/ 1, 1, 1, tsID /)
               imixed_count(:) = (/ nlon, nlat, 1, 1 /)
               if(jk .eq. 1) then
                 CALL IO_GET_VARA_DOUBLE (gpnc2%file_id,mixedid2,imixed_start,imixed_count, zmixedlayer(1:nlon,1:nlat))
               endif
             ENDDO
             zot(:,:,:)=MERGE(zot(:,:,:),xmissing,zot(:,:,:).NE.missing_value)
             zos(:,:,:)=MERGE(zos(:,:,:),xmissing,zos(:,:,:).NE.missing_value)
             zou(:,:,:)=MERGE(zou(:,:,:),xmissing,zou(:,:,:).NE.missing_value)
             zov(:,:,:)=MERGE(zov(:,:,:),xmissing,zov(:,:,:).NE.missing_value)
             zmixedlayer(:,:)=MERGE(zmixedlayer(:,:),xmissing,zmixedlayer(:,:).NE.missing_value)
!!!         zow(:,:,:,3)=MERGE(zow(:,:,:,3),xmissing,zow(:,:,:).NE.missing_value)
           ENDIF
         ENDIF
         CALL IO_close(gpnc2)
         IF ( lwarning_msg.GE.2 ) THEN
            WRITE (nerr,*) 'read_dailygodas, date=',recdate
         ENDIF
       ENDIF
       END SUBROUTINE read_godas_1record
      END SUBROUTINE read_dailygodas
  ! ----------------------------------------------------------------------

!------------------------------------------------------------------------------		
      SUBROUTINE fill_missing2(finout,io_nlon,io_ngl,ndepth,ldeep)
!
!   fill missing value using nearest data points
!   Method: Inverse Distance Weighting Interpolation
!   finout: input/output field
!   
      IMPLICIT NONE
      INTEGER, PARAMETER :: nerr = 6     ! error output stream
      INTEGER,INTENT(in):: io_nlon  ! number of longitudes in NetCDF file
      INTEGER,INTENT(in):: io_ngl   ! number of latitudes in NetCDF file  
      INTEGER,INTENT(in):: ndepth   ! number of depths in NetCDF file
      LOGICAL,INTENT(in):: ldeep    ! fill missig for deeper levels
      REAL,INTENT(in out)::finout(io_nlon,ndepth,io_ngl)
      INTEGER::i,j,k,irad(io_nlon,io_ngl),jrad,nmissing,nfilled,nmax_lon_radius,nmax_lat_radius
      REAL::sumw,aspect,dlon,dlat
      REAL, DIMENSION(-io_nlon/2:(3*io_nlon)/2+1,-io_ngl/2:(3*io_ngl)/2+1):: finm
      INTEGER,  DIMENSION(-io_nlon/2:(3*io_nlon)/2+1,-io_ngl/2:(3*io_ngl)/2+1):: inm
      REAL, PARAMETER:: MAX_LAT_RADIUS=2.   !! maximum radius in lat dir. 2 deg
      REAL, PARAMETER:: MAX_LON_RADIUS=10.  !! maximum radius in lon dir. 10 deg
!
      !!! WRITE(nerr,'(/,A,I2)') ' fill_missing 0.0'
      finm=xmissing
      dlon=360./DBLE(io_nlon)
      dlat=180./DBLE(io_ngl)
      nmax_lon_radius=MAX_LON_RADIUS/dlon
      nmax_lat_radius=MAX_LAT_RADIUS/dlat
      aspect=(MAX_LAT_RADIUS/dlat)/(MAX_LON_RADIUS/dlon)
     !!! WRITE(nerr,'(/,A,I7,A,I7,A,F9.2)') ' fill_missing 1.0: nmax_lon_radius=',nmax_lon_radius, &
     !!!  ' nmax_lat_radius=',nmax_lat_radius,' aspect=',aspect 
!
      irad=0
      DO k = 1,ndepth
      nmissing=0
      nfilled=0
      !!! WRITE(nerr,'(/,A,I2)') ' fill_missing 1.1: k=', k
      finm(1:io_nlon,1:io_ngl)=finout(:,k,:)
      !!! cylinic for longitude (-180 - 0)
      !!! WRITE(nerr,'(/,A,I2)') ' fill_missing 1.2: k=', k
      finm(-io_nlon/2:0,1:io_ngl)=finout(io_nlon/2:io_nlon,k,:)
      !!! cylinic for longitude (360 - 720)
      !!! WRITE(nerr,'(/,A,I2)') ' fill_missing 1.3: k=', k
      finm(io_nlon+1:(3*io_nlon)/2+1,1:io_ngl)=finout(1:io_nlon/2+1,k,:)
      inm=MERGE(1,0,finm.NE.xmissing) ! unit mask
      DO j = 1,io_ngl
        DO i = 1,io_nlon
          IF (finout(i,k,j).EQ.xmissing) THEN
            nmissing=nmissing+1
            IF (k.EQ.1) THEN
            ! Determine Radius
              sumw=0.
!!!              WRITE(nerr,'(/,A,6I7)') ' fill_missing 5.1:',i,j,k
!!!              DO WHILE ( (sumw.LE.0._dp).AND.(irad(i,j).LT.io_nlon/2-1) )
              DO WHILE ( (sumw.LE.0.).AND.(irad(i,j).LT.nmax_lon_radius) )
                irad(i,j)=irad(i,j)+1
                jrad=irad(i,j)*aspect
                sumw=DBLE(SUM(inm(i-irad(i,j):i+irad(i,j),j-jrad:j+jrad)))
              ENDDO
            ELSE
              ! Using the same radius as surface grid, it assumes the depth of the missing value grid is the same as the maximun depth within the radius
              jrad=irad(i,j)*aspect
              sumw=DBLE(SUM(inm(i-irad(i,j):i+irad(i,j),j-jrad:j+jrad)))
              IF (ldeep) THEN
                DO WHILE ( (sumw.LE.0.).AND.(irad(i,j).LT.nmax_lon_radius) )
                  irad(i,j)=irad(i,j)+1
                  jrad=irad(i,j)*aspect
                  sumw=DBLE(SUM(inm(i-irad(i,j):i+irad(i,j),j-jrad:j+jrad)))
                ENDDO
              ENDIF
            ENDIF
            IF (sumw.GT.0.) THEN
              nfilled=nfilled+1
!!!              WRITE(nerr,'(/,A,6I7)') ' fill_missing 5.2:',i,j,k,irad(i,j),jrad,INT(sumw)
              finout(i,k,j)=SUM(finm(i-irad(i,j):i+irad(i,j),j-jrad:j+jrad),mask=finm(i-irad(i,j):i+irad(i,j),j-jrad:j+jrad).NE.xmissing)/sumw
            ELSE
              ! value of deeper layers are missing (not fill) as well.
              ! EXIT
            ENDIF
!!!            WRITE(nerr,'(/,A,6I7)') ' fill_missing 5.3:',i,j,k,irad(i,j),jrad,INT(sumw)          
          ENDIF
        ENDDO
      ENDDO
      IF ( lwarning_msg.GE.2 ) then 
        WRITE(nerr,'(/,A,I5,A,I9,A,I9,A)') ' fill_missing: k=',k,&
             ',   ',nmissing,' missing value grids, where ',nfilled,&
             ' grids filled.'
      ENDIF
    ENDDO
  END SUBROUTINE fill_missing2
!
!--------------------------------------------------------------

!-----
!-----------------------------------------------------------------------------
        SUBROUTINE time_weights(idtg1,tauleft)
!*****  *********************************************

        ! calculates weighting factores for clsst2 and ozone
        !-

!            USE mo_interpo, ONLY : wgt1, wgt2, nmw1, nmw2, nmw1cl,
!            nmw2cl, &
!                               wgtasd1, wgtd2, ndw1, ndw2
!         IMPLICIT NONE

!        TYPE (time_native) :: date_monm1, date_mon, date_monp1
!          INTEGER*8 :: idtg1
!          REAL      :: tauleft
!          REAL      :: wgt1,wgt2
!          INTEGER   :: nmw1,nmw2
           use rank

          integer*8 idtg1
          real tauleft
          character*12 cdtg
          INTEGER   yyyy, mm, dd, hh, mn
          INTEGER   yr,mo,dy,hr
          INTEGER   seconds, isec
          INTEGER   imp1, imm1, imp1cl, imm1cl, imlenm1, imlen, imlenp1
          REAL      zdayl, zsec
          REAL      zmohlf, zmohlfp1, zmohlfm1
!          REAL     zdh, zdhp1, zdhm1
          INTEGER, PARAMETER :: NDAYLEN=86400
          REAL      ydate,ydate1,ydate2,ydate3
          REAL      ydate_jd,ydate1_jd,ydate2_jd,ydate3_jd

          ! ***  set calendar related parameters
          ! -----------------------------------

          write(cdtg,'(i12)')idtg1
          read(cdtg,'(i4,i2,i2,i2,i2)')yyyy,mm,dd,hh,mn
          yr=int(yyyy)
          mo=int(mm)
          dy=int(dd)
          hr=int(hh)
          seconds=int(tauleft)

!        CALL TC_get (date_mon, yr, mo, dy, hr, mn, se)

          ! month index for AMIP data  (0..13)
          imp1 = mo+1
          imm1 = mo-1

          ! month index for cyclic climatological data (1..12)
          imp1cl = mo+1
          imm1cl = mo-1
          IF (imp1cl > 12) imp1cl= 1
          IF (imm1cl <  1) imm1cl=12

          ! *** determine length of months and position within current
          ! month -------

!        CALL TC_set(yr, imm1cl, 1, 0, 0, 0, date_monm1)
!        CALL TC_set(yr, imp1cl, 1, 0, 0, 0, date_monp1)
          imlenm1 = Get_JulianMonLen(yr,imm1cl)
          imlen   = Get_JulianMonLen(yr,mo)
          imlenp1 = Get_JulianMonLen(yr,imp1cl)

          zdayl    = REAL(NDAYLEN)
          zmohlfm1 = REAL(imlenm1*zdayl*0.5)
          zmohlf   = REAL(imlen  *zdayl*0.5)
          zmohlfp1 = REAL(imlenp1*zdayl*0.5)


          ! *** weighting factors for first/second half of month
          ! -------------------

          nmw1   = mo
!          nmw1cl = mo

          ! seconds in the present month
!        CALL TC_get (next_date, days, seconds)
          isec = (dy-1)*NDAYLEN + hr*3600+ seconds
          zsec = REAL(isec)

          IF(zsec <= zmohlf) THEN                     ! first part of month
            wgt1   = (zmohlfm1+zsec)/(zmohlfm1+zmohlf)
            wgt2   = 1.-wgt1
            nmw2   = imm1
!            nmw2cl = imm1cl
          ELSE                                        ! second part of month
            wgt2   = (zsec-zmohlf)/(zmohlf+zmohlfp1)
            wgt1   = 1.-wgt2
            nmw2   = imp1
!            nmw2cl = imp1cl
          ENDIF

!          ! *** weighting factors for first/second half of day
!          -------------------
!
!          ndw1   = 2
!
!          zsec = REAL(seconds,dp)
!          zdh   = 12._dp*3600._dp
!          zdhm1 = zdh
!          zdhp1 = zdh
!          IF( zsec <= zdh ) THEN                     ! first part of day
!            wgtd1  = (zdhm1+zsec)/(zdhm1+zdh)
!            wgtd2  = 1._dp-wgtd1
!            ndw2   = 1
!          ELSE                                       ! second part of day
!            wgtd2  = (zsec-zdh)/(zdh+zdhp1)
!            wgtd1  = 1._dp-wgtd2
!            ndw2   = 3
!          ENDIF


      ! *** weighting factors for GODAS PENTAD data -------------------
         IF (lgodas) THEN
           IF (ldailysst) THEN
!ps          CALL get_date_components(next_date, yr, mo, dy, hr, mn, se)
!ps          ydate =yr*10000._dp+mo*100._dp+dy+(hr+mn/60._dp+se/3600._dp)/24._dp
             ydate=real(yr*10000.+mo*100.+dy)+real(hr/24.)+tauleft/3600./24.
             ydate1=timevals_godas(1)
             ydate2=timevals_godas(2)
             ydate3=timevals_godas(3)
             call date2JulianDay(ydate,ydate_jd) 
             call date2JulianDay(ydate1,ydate1_jd)
             call date2JulianDay(ydate2,ydate2_jd)
             call date2JulianDay(ydate3,ydate3_jd)

             IF (ydate_jd .LE. ydate2_jd) THEN
               wgto2=(ydate_jd-ydate1_jd)/(ydate2_jd-ydate1_jd)
               wgto1=1.-wgto2
               now1=1
               now2=2
             ELSEIF (ydate_jd .LE. ydate3_jd) THEN
               wgto2=(ydate_jd-ydate2_jd)/(ydate3_jd-ydate2_jd)
               wgto1=1.-wgto2
               now1=2
               now2=3
             ELSE
!ps               CALL finish('time_weights', 'GODAS PENTAD Date not found')
               print *,'ydate=',ydate,'timevals_godas(2)=',timevals_godas(2) &
                      ,'timevals_godas(3)=',timevals_godas(3) 
               print *,'time_weights', 'GODAS PENTAD Date not found'
               stop
               CALL mpe_finalize
               CALL dmsexit(-1)
             ENDIF
          !!! IF(wgtd(1).GT.1._dp .OR. wgtd(2).GT.1._dp )THEN
          !!!   WRITE(nerr,*) 'get_5dwgtd yr, mo, dy, hr, mn, se=',yr,
          !mo, dy, hr, mn, se,' ydate=',ydate,' ID=',ID,'
          !timevals_godas=',INT(timevals_godas(0:2)),&
          !!!               ' JD(0:2)=',JD(0:2),' wgtd(1:2)=',wgtd(1:2)
          !!! ENDIF
!             if(myrank .eq. 49) then
!               print *,'time weight dailysst'    &
!                   ,',timevals_godas(1)=',timevals_godas(1)   &
!                   ,',timevals_godas(2)=',timevals_godas(2)   &
!                   ,',timevals_godas(3)=',timevals_godas(3)   &
!                   ,',now1=',now1,',now2=',now2        &
!                   ,',wgto1=',wgto1,',wgto2=',wgto2
!             endif

           ELSE
          ! monthly godas data
             wgto2=wgt2
             wgto1=wgt1
             now1=nmw1
             now2=nmw2
           ENDIF
         ENDIF

         obswtbwgt1=wgto1
         obswtbwgt2=wgto2
         obswtbnmw1=now1
         obswtbnmw2=now2

      ! *** weighting factors for first/second half of day
         IF (ldailyFCTsst .OR. ldailyFCTicesndpt .OR. (dailyClm_option.ge. 1)) THEN
          ydate=real(yr*10000.+mo*100.+dy)+real(hr/24.)+tauleft/3600./24.
          ydate1=timevals_dailyFCT(1)
          ydate2=timevals_dailyFCT(2)
          call date2JulianDay(ydate,ydate_jd) 
          call date2JulianDay(ydate1,ydate1_jd)
          call date2JulianDay(ydate2,ydate2_jd)

          obswtbwgt1 = (ydate2_jd-ydate_jd)/(ydate2_jd-ydate1_jd)
          obswtbwgt2 = 1.-obswtbwgt1
          obswtbnmw1 = 1
          obswtbnmw2 = 2
        
          if(myrank .eq. 49) then
            print *,'time weight, idtg1=',idtg1,',yr=',yr,',mo=',mo &
                   ,',dy=',dy,',hr=',hr,',tauleft=',tauleft  & 
                   ,',ydate=',ydate,',ydate1=',ydate1        &
                   ,',ydate2=',ydate2,',ydatejd=',ydate_jd   &
                  ,',ydate1jd=',ydate1_jd,',ydate2jd=',ydate2_jd &
                  ,',obswtbwgt1=',obswtbwgt1   &
                  ,',obswtbwgt2=',obswtbwgt2
          endif

         ENDIF



!-----------------------------------------------------
        CONTAINS

          SUBROUTINE date2JulianDay(zdate,jd)
            real, INTENT(IN):: zdate
            real, INTENT(OUT):: jd
            integer:: yr,mo,dy
            real:: zsec

            call get_date_component(zdate, yr, mo, dy, zsec)
            jd=Set_JulianDay(yr, mo, dy, zsec) 

          END SUBROUTINE date2JulianDay

          SUBROUTINE get_date_component(zdate, kyr, kmo, kdy, zsec)

          real, INTENT(IN):: zdate
          integer:: kyr,kmo,kdy
          integer:: tmp,yyyy,mm,dd
          character*8 :: cdtg
          real:: zsec        !second [real] input (seconds of the day)


          tmp=int(zdate)
          write(cdtg,'(i8)')tmp
          read(cdtg,'(i4,i2,i2)')yyyy,mm,dd
          kyr=int(yyyy)
          kmo=int(mm)
          kdy=int(dd)
          zsec=zdate-real(tmp)

          END SUBROUTINE get_date_component
 


          FUNCTION Get_JulianMonLen(ky, km) RESULT(idmax)
        !+
        !
        ! Get_JulianMonLen [function, integer]
        !    get the length of a months in a Julian year
        !    (
        !    year  [integer] input (Calendar year)
        !    month [integer] input (month of the year)
        !    )
        !
        !-
           INTEGER, INTENT(in) :: km, ky
           INTEGER :: idmax

           SELECT CASE(km)
           CASE(1,3,5,7,8,10,12);  idmax = 31
           CASE(4,6,9,11);         idmax = 30
           CASE(2)
             IF ( (MOD(ky,4)==0 .AND. MOD(ky,100)/=0) .OR. MOD(ky,400)==0 ) THEN
               ! leap year found
               idmax = 29
             ELSE
               idmax = 28
             END IF

           CASE default
             print*,'mo_time_weight:Get_JulianMonLen, month invalid'

           END SELECT
!           Get_JulianMonLen = idmax

          END FUNCTION Get_JulianMonLen


          FUNCTION Set_JulianDay(ky, km, kd, zsec) RESULT(zd)
     !+
     !
     ! Set_JulianDay  [subroutine]
     !    convert year, month, day, seconds into Julian calendar day
     !    (
     !    year   [integer] input (calendar year)
     !    month  [integer] input (month of the year)
     !    day    [integer] input (day of the month)
     !    second [integer] input (seconds of the day)
     !    date   [julian_date] output (Julian day)
     !    )
     !
     !-
     !
          INTEGER, INTENT(IN) :: ky
          INTEGER, INTENT(IN) :: km
          INTEGER, INTENT(IN) :: kd
          REAL,    INTENT(IN) :: zsec
     !
     ! for reference: 1. January 1998 00 UTC === Julian Day 2450814.5
     !
          INTEGER :: ib, iy, im, idmax
          REAL :: zd

          IF ( zsec > 1. ) THEN
            print*,'Set_JulianDay: invalid number of seconds'
          ENDIF
          

          IF (km <= 2) THEN
            iy = ky-1
            im = km+12
          ELSE
            iy = ky
            im = km
          ENDIF

         ib = INT(iy/400)-INT(iy/100)
        ! check the length of the month
         idmax = Get_JulianMonLen (ky, km)

         IF (kd < 1 .OR. idmax < kd) &
           print *,'Set_JulianDay: day in months invalid'

         zd = real(365.25*iy)+INT(30.6001*(im+1)) &
             +REAL(ib)+1720996.5+REAL(kd)+zsec

         END FUNCTION Set_JulianDay


        END SUBROUTINE time_weights







      END MODULE mod_sst
