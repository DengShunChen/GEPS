#define LCWBGFS T
#if defined(LCWBGFS)
#define nlon nx
#define ngl my
#define nlat my
#define mpp_root_pe() 0
#define p_pe myrank
#define sstmin 270
#endif

!    read wtfn12, wsfn12 data
      MODULE mod_sst

        use param
        use mpe
        use index
        use const,             ONLY:ggdef,ifilin_nc,ifilin_sst,ifilin_ncep  &
                                   ,ldailyFCTsst,ldailyFCTicesndpt,ifilin   &
                                   ,dailyClm_option,ifilin_ClmANA,ifilin_ClmFCT &
                                   ,do_sit
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
        PUBLIC :: nwdepth, wdepths, wtfn12, wsfn12,wtfn1st,mask1st
        PUBLIC :: albice,albsn,albw,csn,cice,rhosn,rhoice,xkice,xksn,xkw,         &
                  omegas,wcri
        PUBLIC :: deallocate_ocaf_array,deallocate_woa0_array,deallocate_godas_array
!ps				  ,tol,wlvlref,dpthmx,init_sit_ocean

        PUBLIC :: opgsst
        PUBLIC :: nmw1,nmw2,wgt1,wgt2
        PUBLIC :: now1,now2,wgto1,wgto2             !! GODAS MONTHLY/PENTAD Data
        PUBLIC :: obswtbnmw1,obswtbnmw2,obswtbwgt1,obswtbwgt2     !obswtb Data
        PUBLIC :: dailyFCTsst,dailyFCTcice,dailyFCTsndepth
        PUBLIC :: dFCTsstdt,dFCTcicedt,dFCTsndepthdt
        PUBLIC :: ANAsstT0,dailyClmANAsst,dailyClmFCTsst
        PUBLIC :: deallocate_dailyFCT_array
        PUBLIC :: outtseadiffFCT24
        PUBLIC :: obswtbold,obswtbnow,obswtbnew
        PUBLIC :: tseaold,tseanow,tseanew,dtseadt


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
        REAL, ALLOCATABLE :: mask1st(:,:) ! (nlon,nwdepth,ngl) in global coordinates,

  !! memory pointer for pentad GODAS OCEAN data (lgodas & ldailysst)
        REAL      :: timevals_godas(3) = 0.  ! absoulte time (e.g., 19971003.25) GODAS PENTAD Data
        INTEGER   :: files_godas(3)          ! fileid for pentad data at day-1, day+0 and day+1, respectively
        INTEGER   :: tsID_godas(3)           ! record id for pentad data at day-1, day+0 and day+1, in its repective file, respectively
        INTEGER   :: nts_godas(3)            ! # of timestamps for pentad data at day-1, day+0 and day+1, in its repective file, respectively

  !! memory pointer for daily FCTsst data (ldailyFCTsst)
        REAL      :: timevals_dailyFCT(2)= 0.       ! absoulte time (e.g.,19971003.25) daily forecast sst Data
        REAL, ALLOCATABLE :: dailyFCTsst(:,:,:)     ! (nlon,ngl,2) at ydate, ydate+1 day in global coordinates,
                                                    ! forecast daily water tempeature (K)
        REAL, ALLOCATABLE :: dailyFCTcice(:,:,:)    ! (nlon,ngl,2) at ydate, ydate+1 day in global coordinates,
                                                    ! forecast sea ice fration
        REAL, ALLOCATABLE :: dailyFCTsndepth(:,:,:) ! (nlon,ngl,2) at ydate, ydate+1 day in global coordinates,
                                                    ! forecast daily snow depth (mm)
        REAL, ALLOCATABLE :: ANAsstT0(:,:)          ! (nlon,ngl), analysis SST at tau=0 
        REAL, ALLOCATABLE :: dailyClmANAsst(:,:,:)  ! (nlon,ngl,2) at tau=0, ydate, ydate+1 day in global coordinates,
        REAL, ALLOCATABLE :: dailyClmFCTsst(:,:,:)  ! (nlon,ngl,2) at ydate, ydate+1 day in global coordinates,
                                                    ! forecast climatology daily water tempeature (K)
        REAL, ALLOCATABLE :: dFCTsstdt(:,:)         ! d(dailyFCTsst)/dt (K/s)
        REAL, ALLOCATABLE :: dFCTcicedt(:,:)        ! d(dailyFCTcice)/dt (K/s)
        REAL, ALLOCATABLE :: dFCTsndepthdt(:,:)     ! d(dailyFCTsndepth)/dt (K/s)
        REAL, ALLOCATABLE :: obswtbold(:,:)           ! calculated obs. sst per (n-1) timestep
        REAL, ALLOCATABLE :: obswtbnow(:,:)           ! calculated obs. sst now (n) timestep 
        REAL, ALLOCATABLE :: obswtbnew(:,:)           ! calculated obs. sst next (n+1) timestep 
        REAL, ALLOCATABLE :: tseadiffFCT(:,:)       ! the change of tg from dta*dFCTsstdt 
        REAL, ALLOCATABLE :: tseadiffFCT24(:,:)     ! the average change of tg from dta*dFCTsstdt 
        REAL, ALLOCATABLE :: dtseadt(:,:)           ! the change rate of tg from FCTsst and SIT
        REAL, ALLOCATABLE :: tseaold(:,:)             ! sst at pre. timestep (n-1)
        REAL, ALLOCATABLE :: tseanow(:,:)             ! sst at now timestep (n) 
        REAL, ALLOCATABLE :: tseanew(:,:)             ! sst at next timestep (n+1)


!for opgsst
      REAL, dimension(:,:,:),allocatable,save :: opgsst
!for time_interpolation
      REAL wgt1,wgt2,obswtbwgt1,obswtbwgt2,wgto1,wgto2
      INTEGER nmw1,nmw2,obswtbnmw1,obswtbnmw2,now1,now2


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



        
      CONTAINS


!---------------------------------------------------------
        subroutine allocate_opgsst_array

          integer  ierr
          allocate ( opgsst(nxp,my_max,0:13), stat=ierr)
             if (ierr/= 0) then
                 write(6,*) 'mod_opgsst : allocate fail 1 '
                 call dmsexit(-1)
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
          use const, only: ihdgi

          character*4 ggdef
          integer*8 idtg1
          character*12 cdtg
          integer   iyyyy,mm,ddhhmn,iyy,imm,istat
          integer lncrec
          integer i,j,ii,jj,nxj,kmm
          real,dimension(:,:), allocatable:: temp1
          logical ocean(nxp,my_max),ice(nxp,my_max)
 

          allocate(temp1(nx,my))
!
#ifdef I38K
  11     format('W00100','000000',a4,i4.4,i2.2,'010000')  ! W10
#else
  11     format('W00100','0000',a4,i4.4,i2.2,'010000')  ! W10
#endif

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

             write(ihdgi,11) ggdef,iyy,imm
             call dmsread(nx,my,lncrec,'H',ifilin_sst,temp1(:,:),istat)
!             if( lreduce.eq.1 ) call reducepickr(temp1(1,imm,1),nxdef,nx,my)
             do jj = 1, jlistnum
               j=jlist1(jj)
               i=nxjstart(j)
               nxj=nxdef_2d(j)
               if( lreduce.eq.1 )call reducepickr (temp1(1,j),nxdef(j),nx,1)
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

           deallocate (obswtbold,obswtbnow,obswtbnew)
           deallocate (tseaold,tseanow,tseanew,dtseadt)
           deallocate (tseadiffFCT,tseadiffFCT24)
           deallocate (dailyFCTsst,dFCTsstdt)
           if(ldailyFCTicesndpt) then
             deallocate (dailyFCTcice,dailyFCTsndepth) 
             deallocate (dFCTcicedt,dFCTsndepthdt)
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


        SUBROUTINE read_dailyFCT(idtg1,tau,dtx,tg1,cice1,sndepth1  &
                                ,plon,plat,locean)

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
          real plon(nx,my_max),plat(my)
          logical  locean(nxp,my_max)

          icurrenttau=int(tau)
          tauleft=float(int((tau-int(tau)+0.001)*3600./dtx))*dtx   !(sec)
          if(tauleft > 3599.) then
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
          if(tauleft > 3599. ) then
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


          IF (.NOT. ALLOCATED(dailyFCTsst)) ALLOCATE (dailyFCTsst(nxp,my_max,2))
          IF (.NOT. ALLOCATED(dFCTsstdt)) ALLOCATE (dFCTsstdt(nxp,my_max))
          IF (.NOT. ALLOCATED(obswtbold)) ALLOCATE (obswtbold(nxp,my_max))
          IF (.NOT. ALLOCATED(obswtbnow)) ALLOCATE (obswtbnow(nxp,my_max))
          IF (.NOT. ALLOCATED(obswtbnew)) ALLOCATE (obswtbnew(nxp,my_max))
          IF (.NOT. ALLOCATED(tseaold)) ALLOCATE (tseaold(nxp,my_max))
          IF (.NOT. ALLOCATED(tseanow)) ALLOCATE (tseanow(nxp,my_max))
          IF (.NOT. ALLOCATED(tseanew)) ALLOCATE (tseanew(nxp,my_max))
          IF (.NOT. ALLOCATED(dtseadt)) then 
            ALLOCATE (dtseadt(nxp,my_max))
            dtseadt=0.
          ENDIF
          IF (.NOT. ALLOCATED(tseadiffFCT)) then
            ALLOCATE (tseadiffFCT(nxp,my_max))
            tseadiffFCT=0.
          ENDIF
          IF (.NOT. ALLOCATED(tseadiffFCT24)) then
            ALLOCATE (tseadiffFCT24(nxp,my_max))
            tseadiffFCT24=0.
          ENDIF

          IF(ldailyFCTicesndpt) THEN
            IF (.NOT. ALLOCATED(dailyFCTcice)) ALLOCATE (dailyFCTcice(nxp,my_max,2))
            IF (.NOT. ALLOCATED(dailyFCTsndepth)) ALLOCATE (dailyFCTsndepth(nxp,my_max,2))
            IF (.NOT. ALLOCATED(dFCTcicedt)) ALLOCATE (dFCTcicedt(nxp,my_max))
            IF (.NOT. ALLOCATED(dFCTsndepthdt)) ALLOCATE (dFCTsndepthdt(nxp,my_max))
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
            timevals_dailyFCT(1)=ydate
            dailyFCTsst(:,:,1)=tg1(:,:)
            obswtbold(:,:)=tg1(:,:)
            obswtbnow(:,:)=tg1(:,:)
            tseaold(:,:)=tg1(:,:)
            tseanow(:,:)=tg1(:,:)
            if(dailyClm_option .ge. 1) ANAsstT0(:,:)=tg1(:,:)
            if(ldailyFCTicesndpt) then
              dailyFCTcice(:,:,1)=cice1(:,:)
              dailyFCTsndepth(:,:,1)=sndepth1(:,:)
            endif
            if(ldailyFCTsst .or. ldailyFCTicesndpt)then
              CALL read_dailyFCT_dayp1(idtg_FCT,locean)
            endif
            if(dailyClm_option .ge. 1)then
              CALL read_dailyClm_2days(idtg1,idtg_FCT,icurrenttau,plon,plat,locean)
            endif

            timevals_dailyFCT(2)=ydate2
          ELSEIF (ydate.LT.timevals_dailyFCT(2)) THEN
     ! data were read. Note that initial value of  timevals_godas=0.
            RETURN
          ELSE
     ! note that initial value of  timevals_godas=0.
     ! Shift left
            dailyFCTsst(:,:,1)=dailyFCTsst(:,:,2)
            if(ldailyFCTicesndpt) then
              dailyFCTcice(:,:,1)=dailyFCTcice(:,:,2)
              dailyFCTsndepth(:,:,1)=dailyFCTsndepth(:,:,2)
            endif
            timevals_dailyFCT(1)=timevals_dailyFCT(2)
            if(ldailyFCTsst .or. ldailyFCTicesndpt) then
              CALL read_dailyFCT_dayp1(idtg_FCT,locean)    ! read next FCST data
            endif
            if(dailyClm_option .ge. 1)then
              dailyClmANAsst(:,:,1)=dailyClmANAsst(:,:,2)
              if(dailyClm_option .eq. 2) then
                dailyClmFCTsst(:,:,1)=dailyClmFCTsst(:,:,2)
              endif
              CALL read_dailyClm_dayp1(idtg1,idtg_FCT,icurrenttau,plon,plat,locean)
            endif

            timevals_dailyFCT(2)=ydate2

          ENDIF
          



      CONTAINS
     !-----------------------------------
        SUBROUTINE read_dailyFCT_dayp1(idtg1,locean)

          use const, only: ihdgi
          INTEGER*8 idtg1
          INTEGER iyyyy,imm,idd,ihh,imn
          INTEGER lncrec
          REAL ssttemp(nx,my),cicetemp(nx,my),sndpttemp(nx,my)
          INTEGER ii,jj,nxj,i,j,istat,nxjtot
          INTEGER itemp,jtemp
          INTEGER iitemp,jjtemp
          REAL sst_suntemp,sst_counttemp 
          REAL cice_suntemp,cice_counttemp,sndpt_suntemp,sndpt_counttemp
          LOGICAL locean(nxp,my_max)

          write(cdtg,'(i12)') idtg1
          read(cdtg,'(i4.4,i2.2,i2.2,i2.2,i2.2)') iyyyy,imm,idd,ihh,imn

          lncrec=nx*my
          
          ssttemp=0.
          cicetemp=0.
          sndpttemp=0.

#ifdef I38K
   11     format('W00100','000000',a4,i4.4,i2.2,i2.2,i2.2,'00')  ! sea surface temperature
   12     format('W00091','000000',a4,i4.4,i2.2,i2.2,i2.2,'00')  ! sea ice fration
   13     format('B00650','000000',a4,i4.4,i2.2,i2.2,i2.2,'00')  ! water equivlent snow depth
#else
   11     format('W00100','0000',a4,i4.4,i2.2,i2.2,i2.2,'00')  ! sea surface temperature
   12     format('W00091','0000',a4,i4.4,i2.2,i2.2,i2.2,'00')  ! sea ice fration
   13     format('B00650','0000',a4,i4.4,i2.2,i2.2,i2.2,'00')  ! water equivlent snow depth
#endif
          write(ihdgi,11) ggdef,iyyyy,imm,idd,ihh
          call dmsread(nx,my,lncrec,'H',ifilin_sst,ssttemp(:,:),istat)
          if(ldailyFCTicesndpt)then
            write(ihdgi,12) ggdef,iyyyy,imm,idd,ihh
            call dmsread(nx,my,lncrec,'H',ifilin_ncep,cicetemp(:,:),istat)
            write(ihdgi,13) ggdef,iyyyy,imm,idd,ihh
            call dmsread(nx,my,lncrec,'H',ifilin_ncep,sndpttemp(:,:),istat)
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
              call reducepickr (ssttemp(1,j),nxdef(j),nx,1)
              if(ldailyFCTicesndpt)then
                call reducepickr (cicetemp(1,j),nxdef(j),nx,1)
                call reducepickr (sndpttemp(1,j),nxdef(j),nx,1)
              endif
            endif

            DO ii=1,nxj
              i=nxjstart(j)+ii-1
              dailyFCTsst(ii,jj,2)=MERGE(ssttemp(i,j),xmissing,  &
                   (ssttemp(i,j).GE.sstmin .AND. ssttemp(i,j) .LE. 400.))
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
                      jjtemp=1 
                    elseif (jtemp .gt. my) then
                      jjtemp=my
                    else
                      jjtemp=jtemp
                    endif
                    if (ssttemp(iitemp,jjtemp).GE.sstmin .AND. ssttemp(iitemp,jjtemp) .LE. 400.)then
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
                if(ssttemp(i,j).LT.sstmin .OR. ssttemp(i,j) .GT. 400.)then
                  ssttemp(i,j)=sst_suntemp/sst_counttemp
                  if (ssttemp(i,j).GE.sstmin .AND. ssttemp(i,j) .LE. 400.)then
                    dailyFCTsst(ii,jj,2)=ssttemp(i,j)
                  else
                    dailyFCTsst(ii,jj,2)=dailyFCTsst(ii,jj,1)
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
              if(locean(ii,jj) .AND. dailyFCTsst(ii,jj,1).ge.sstmin &
                 .AND. dailyFCTsst(ii,jj,2).ge.sstmin)then 
                dFCTsstdt(ii,jj)=(dailyFCTsst(ii,jj,2)-dailyFCTsst(ii,jj,1))/(24.*3600.)
              else
                dFCTsstdt(ii,jj)=0.
              endif
              if(ldailyFCTicesndpt)then
                dFCTcicedt(ii,jj)=(dailyFCTcice(ii,jj,2)-dailyFCTcice(ii,jj,1))/(24.*3600.)
                dFCTsndepthdt(ii,jj)=(dailyFCTsndepth(ii,jj,2)-dailyFCTsndepth(ii,jj,1))/(24.*3600.)
              endif

            ENDDO  !end do ii
          ENDDO    !end do jj

        
        END SUBROUTINE read_dailyFCT_dayp1


        SUBROUTINE read_dailyClm_2days(idtg1,idtg_2,itau,plon,plat,locean)

          use const, only: ihdgi
          INTEGER*8 idtg1,idtg_2
          INTEGER itau
          INTEGER iyyyy,imm,idd,ihh,imn
          INTEGER iyyyy2,imm2,idd2,ihh2,imn2
          INTEGER lncrec
          REAL sstANA0(nx,my),sstANA1(nx,my)
          REAL sstFCT0(nx,my),sstFCT1(nx,my)
          INTEGER i,j,ii,jj,nxj,istat
          REAL wweight
          REAL plon(nx,my_max),plat(my)
          LOGICAL locean(nxp,my_max)

          write(cdtg,'(i12)') idtg1
          read(cdtg,'(i4.4,i2.2,i2.2,i2.2,i2.2)') iyyyy,imm,idd,ihh,imn
          write(cdtg,'(i12)') idtg_2
          read(cdtg,'(i4.4,i2.2,i2.2,i2.2,i2.2)') iyyyy2,imm2,idd2,ihh2,imn2

          lncrec=nx*my
          sstANA0=0.
          sstANA1=0.
          sstFCT0=0.
          sstFCT1=0.

#ifdef I38K
   11     format('W00100',6x,a4,4x,i2.2,i2.2,4x)  ! sea surface temperature
   13     format('W00100',i6.6,a4,4x,i2.2,i2.2,4x)  ! sea surface temperature
#else
   11     format('W00100',4x,a4,4x,i2.2,i2.2,4x)  ! sea surface temperature
   13     format('W00100',i4.4,a4,4x,i2.2,i2.2,4x)  ! sea surface temperature
#endif
          !dailyClm_option>=1, read climatology ana. sst
          write(ihdgi,11) ggdef,imm,idd
          call dmsread(nx,my,lncrec,'H',ifilin_ClmANA,sstANA0(:,:),istat)

          write(ihdgi,11) ggdef,imm2,idd2
          call dmsread(nx,my,lncrec,'H',ifilin_ClmANA,sstANA1(:,:),istat)

!          if( lreduce.eq.1 ) then
!            call reducepickr(sstANA0(1,1),nxdef,nx,my)
!            call reducepickr(sstANA1(1,1),nxdef,nx,my)
!          endif


          sstANA0=MERGE(sstANA0,xmissing,(sstANA0.GE.sstmin .AND. sstANA0.LE.400.))
          sstANA1=MERGE(sstANA1,xmissing,(sstANA1.GE.sstmin .AND. sstANA1.LE.400.))

!ps          CALL fill_missing2(sstANA0(:,:),nx,my,1,.FALSE.)
!ps          CALL fill_missing2(sstANA1(:,:),nx,my,1,.FALSE.)


          if(dailyClm_option .eq. 2) then     !dailyClm_option=2, read forcast climatology sst
            write(ihdgi,13) itau,ggdef,imm,idd
            call dmsread(nx,my,lncrec,'H',ifilin_ClmFCT,sstFCT0(:,:),istat)

            write(ihdgi,13) itau+24,ggdef,imm,idd
            call dmsread(nx,my,lncrec,'H',ifilin_ClmFCT,sstFCT1(:,:),istat)


!            if( lreduce.eq.1 ) then
!              call reducepickr(sstFCT0(1,1),nxdef,nx,my)
!              call reducepickr(sstFCT1(1,1),nxdef,nx,my)
!            endif 

            sstFCT0=MERGE(sstFCT0,xmissing,(sstFCT0.GE.sstmin .AND. sstFCT0.LE.400.))
            sstFCT1=MERGE(sstFCT1,xmissing,(sstFCT1.GE.sstmin .AND. sstFCT1.LE.400.))

!ps            CALL fill_missing2(sstFCT0(:,:),nx,my,1,.FALSE.)
!ps            CALL fill_missing2(sstFCT1(:,:),nx,my,1,.FALSE.)
          endif


          DO jj = 1, jlistnum
            j=jlist1(jj)
            nxj=nxdef_2d(j)
            i=nxjstart(j)
            if( lreduce.eq.1 ) then
              call reducepickr(sstANA0(1,j),nxdef(j),nx,1)
              call reducepickr(sstANA1(1,j),nxdef(j),nx,1)
              call reducepickr(sstFCT0(1,j),nxdef(j),nx,1)
              call reducepickr(sstFCT1(1,j),nxdef(j),nx,1)
            end if
            DO ii=1,nxj
              wweight=0.
              i=nxjstart(j)+ii-1
              dailyClmANAsst(ii,jj,0)=sstANA0(i,j)
              dailyClmANAsst(ii,jj,1)=sstANA0(i,j)
              dailyClmANAsst(ii,jj,2)=sstANA1(i,j)
              if(locean(ii,jj) .AND. dailyClm_option .eq. 1) then   !dailyClm_option=1
                if(ldailyFCTsst)then
                !similar to (Yuejian Zhu, operational)             
                !  wweight=min(float(itau)/24./35.,1.)
                !do weigting when (tau>=14.*24. & abs(plat)>=40.)
                  if((itau .gt. 14*24) .AND. (abs(plat(j)) .gt. 40.)) then
                    wweight=min(max((float(itau)-14.*24.)/(24.*30.),0.),1.)
                  endif
                  dailyFCTsst(ii,jj,2)=(1.-wweight)*dailyFCTsst(ii,jj,2)    &
                                      +wweight*dailyClmANAsst(ii,jj,2)
                else
                !persistent anomaly sst
                !SSTf_t=[SSTa_t0-SSTc_t0]*exp(-(t-t0)/90)+SSTc_t
                !(Yuejian Zhu, operational)
                  wweight=min(max(exp(-float(itau)/(90.*24.)),0.),1.)
                  dailyFCTsst(ii,jj,2)=wweight*(ANAsstT0(ii,jj)-dailyClmANAsst(ii,jj,0))  &
                                      +dailyClmANAsst(ii,jj,2)
                endif
              endif 
              if(ldailyFCTsst .AND. locean(ii,jj) .AND. dailyClm_option .eq. 2) then   !dailyClm_option=2
                dailyClmFCTsst(ii,jj,0)=sstFCT0(i,j)
                dailyClmFCTsst(ii,jj,1)=sstFCT0(i,j)
                dailyClmFCTsst(ii,jj,2)=sstFCT1(i,j)
                !idea from Yuejian Zhu(2018 JGR)
                 wweight=min(max(float(itau)/24./35.,1.),0.)
                 dailyFCTsst(ii,jj,2)=(1.-wweight)*(ANAsstT0(ii,jj)-dailyClmANAsst(ii,jj,0)    &
                     +dailyClmANAsst(ii,jj,2) )+wweight*(dailyFCTsst(ii,jj,2)    &
                     -(dailyClmFCTsst(ii,jj,2)-dailyClmANAsst(ii,jj,2)))
              endif
              if(locean(ii,jj) .AND. dailyFCTsst(ii,jj,1).ge.sstmin &
                 .AND.dailyFCTsst(ii,jj,2).ge.sstmin)then
                dFCTsstdt(ii,jj)=(dailyFCTsst(ii,jj,2)-dailyFCTsst(ii,jj,1))/(24.*3600.)
              else
                dFCTsstdt(ii,jj)=0.
              endif

            ENDDO  !end do ii
          ENDDO    !end do jj
        
        END SUBROUTINE read_dailyClm_2days

               
        SUBROUTINE read_dailyClm_dayp1(idtg1,idtg_2,itau,plon,plat,locean)

          use const, only: ihdgi
          INTEGER*8 idtg1,idtg_2
          INTEGER itau
          INTEGER iyyyy,imm,idd,ihh,imn
          INTEGER iyyyy2,imm2,idd2,ihh2,imn2
          INTEGER lncrec
          REAL sstANA(nx,my),sstFCT(nx,my)
          INTEGER i,j,ii,jj,nxj,istat
          REAL wweight
          REAL plon(nx,my_max),plat(my)
          LOGICAL locean(nxp,my_max)
          real::ssttend

          write(cdtg,'(i12)') idtg1
          read(cdtg,'(i4.4,i2.2,i2.2,i2.2,i2.2)') iyyyy,imm,idd,ihh,imn
          write(cdtg,'(i12)') idtg_2
          read(cdtg,'(i4.4,i2.2,i2.2,i2.2,i2.2)') iyyyy2,imm2,idd2,ihh2,imn2

          lncrec=nx*my
          sstANA=0.          
          sstFCT=0.          

#ifdef I38K
   11     format('W00100',6x,a4,4x,i2.2,i2.2,4x)  ! sea surface temperature
   13     format('W00100',i6.6,a4,4x,i2.2,i2.2,4x)  ! sea surface temperature
#else
   11     format('W00100',4x,a4,4x,i2.2,i2.2,4x)  ! sea surface temperature
   13     format('W00100',i4.4,a4,4x,i2.2,i2.2,4x)  ! sea surface temperature
#endif
          write(ihdgi,11) ggdef,imm2,idd2
          call dmsread(nx,my,lncrec,'H',ifilin_ClmANA,sstANA(:,:),istat)

!          if( lreduce.eq.1 ) then
!            call reducepickr(sstANA(1,1),nxdef,nx,my)
!          endif

          sstANA=MERGE(sstANA,xmissing,(sstANA.GE.sstmin .AND. sstANA.LE.400.))
!ps          CALL fill_missing2(sstANA(:,:),nx,my,1,.FALSE.)



          if(dailyClm_option .eq. 2) then
            write(ihdgi,13) itau,ggdef,imm,idd
            call dmsread(nx,my,lncrec,'H',ifilin_ClmFCT,sstFCT(:,:),istat)

!            if( lreduce.eq.1 ) then
!              call reducepickr(sstFCT(1,1),nxdef,nx,my)
!            endif

            sstFCT=MERGE(sstFCT,xmissing,(sstFCT.GE.sstmin .AND. sstFCT.LE.400.))
!ps            CALL fill_missing2(sstFCT(:,:),nx,my,1,.FALSE.)
          endif


          DO jj = 1, jlistnum
            j=jlist1(jj)
            nxj=nxdef_2d(j)
            if( lreduce.eq.1 ) then 
              call reducepickr(sstANA(1,j),nxdef(j),nx,1)
              call reducepickr(sstFCT(1,j),nxdef(j),nx,1)
            end if
            DO ii=1,nxj
             wweight=0.
             i=nxjstart(j)+ii-1
             dailyClmANAsst(ii,jj,2)=sstANA(i,j)
             if(locean(ii,jj) )then
              if( dailyClm_option .eq. 1) then   !dailyClm_option=1
                if(ldailyFCTsst)then
                !similar to (Yuejian Zhu, operational)             
                !  wweight=min(float(itau)/24./35.,1.)
                !do weigting when (tau>=14.*24. & abs(plat)>=40.)
                  if((itau .gt. 14*24) .AND. (abs(plat(j)) .gt. 40.)) then
                    wweight=min(max((float(itau)-14.*24.)/(24.*30.),0.),1.)
                  endif
                  dailyFCTsst(ii,jj,2)=(1.-wweight)*dailyFCTsst(ii,jj,2) &
                                      +wweight*dailyClmANAsst(ii,jj,2)
                else
                !persistent anomaly sst
                !SSTf_t=[SSTa_t0-SSTc_t0]*exp(-(t-t0)/90)+SSTc_t
                !(Yuejian Zhu, operational)
                  wweight=min(max(exp(-float(itau)/(90.*24.)),0.),1.)
                  dailyFCTsst(ii,jj,2)=wweight*(ANAsstT0(ii,jj)-dailyClmANAsst(ii,jj,0))  &
                                       +dailyClmANAsst(ii,jj,2)
                endif
              elseif( ldailyFCTsst .AND. dailyClm_option .eq. 2) then   !dailyClm_option=2
                dailyClmFCTsst(ii,jj,2)=sstFCT(i,j)
                !idea from Yuejian Zhu(2018 JGR)
                 wweight=min(max(float(itau)/24./35.,0.),1.)
                 dailyFCTsst(ii,jj,2)=(1.-wweight)*(ANAsstT0(ii,jj)-dailyClmANAsst(ii,jj,0)   &
                     +dailyClmANAsst(ii,jj,2))+wweight*(dailyFCTsst(ii,jj,2)    &
                     -(dailyClmFCTsst(ii,jj,2)-dailyClmANAsst(ii,jj,2)))
              endif
              if( dailyFCTsst(ii,jj,1).ge.sstmin .AND. &
                  dailyFCTsst(ii,jj,2).ge.sstmin )then
                !remove large tendency 20231207
                ssttend =  dailyFCTsst(ii,jj,2)-dailyFCTsst(ii,jj,1)
                ssttend = min( 10. , max( -10. , ssttend ))
                dFCTsstdt(ii,jj)=ssttend/(24.*3600.)
              else
                dFCTsstdt(ii,jj)=0.
              endif
             endif !if locean

            ENDDO
          ENDDO

        
        END SUBROUTINE read_dailyClm_dayp1
           

      END SUBROUTINE read_dailyFCT

!---------------------------------------------------------------
      SUBROUTINE outtseadiffFCT24(nx,my,my_max,dt24,itau,idtg,ggdef)

      use index 
      use mpe
      use const, only: kflag,RTYPE,outdms,outgrb2,ihdgo,ihdgo2
      use mod_grb2_param,only: wrt_grb2_v2

      implicit none

      integer   nx,my,my_max,itau
      real      dt24
      real(kind=RTYPE) glob(nx,my),wrk(nxp,my_max)
      integer*8 idtg
      character*4  ggdef
      integer   imax,jmax,lenc,j,nxj,i,istat,jj

      imax=nx
      jmax=my
      lenc= imax*jmax
!
      do jj = 1,jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i=1,nxj
         wrk(i,jj)=tseadiffFCT24(i,jj)/dt24
        enddo
      enddo
      call unify_reduceintp(nx,my,my_max,wrk,glob)
      call syslbl_w ('w0001f',idtg,itau,ggdef)
      if(outdms.gt.0) &
      call dmswrit(imax,jmax,lenc,kflag,glob,istat)
      if( outgrb2==1.and.myrank==0 )then
      ihdgo2 = ihdgo
      call wrt_grb2_v2(itau,10,3,199,6,168,0,2,glob)
      endif
      tseadiffFCT24=0.

      END SUBROUTINE 
!---------------------------------------------------------------
        subroutine deallocate_ocaf_array  

           deallocate (wdepths,wtfn12,wsfn12)
           if(locaf0) deallocate (wtfn1st,mask1st)
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
          use const, only: ihdgi 
        
          integer nx,my,lkvl
          character*4 ggdef
!          integer*8 idtg1
!          character*12 cdtg
!          integer   iyyyy,mmddhhmn,iyy,imm,istat
          integer istat,k,mm,lncrec
          integer j,jj,nxj,i,ii
          REAL, ALLOCATABLE, TARGET :: temp1(:,:),temp2(:,:)
          logical flag
!       
          if(.not. ALLOCATED(temp1)) allocate (temp1(nx,my))
          if(.not. ALLOCATED(temp2)) allocate (temp2(nx,my))

          if(.not. ALLOCATED(wdepths)) ALLOCATE(wdepths(1:lkvl+2)) 
          if(.not. ALLOCATED(wtfn12)) ALLOCATE(wtfn12(nxp,1:lkvl+2,my_max,0:13)) 
          if(.not. ALLOCATED(wsfn12)) ALLOCATE(wsfn12(nxp,1:lkvl+2,my_max,0:13)) 

          nwdepth=lkvl+2
          wdepths(1:lkvl+2)=sit_zdepth(0:lkvl+1)
          wtfn12=0.
          wsfn12=0.

#ifdef I38K         
  11     format(i3.3,'TFM','  gbck',a4,4x,i2.2,6x)  ! ???TFM
  12     format(i3.3,'SFM','  gbck',a4,4x,i2.2,6x)  ! ???SFM
#else
  11     format(i3.3,'TFM','gbck',a4,4x,i2.2,6x)  ! ???TFM
  12     format(i3.3,'SFM','gbck',a4,4x,i2.2,6x)  ! ???SFM
#endif

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
            call mpe_finalize
            call dmsexit(-1)
            stop
          endif


          lncrec=nx*my
          do k=1,lkvl+2
            do mm=1,12
              write(ihdgi,11) k-1,ggdef,mm
              if(myrank .eq. 0) print *, 'ihdg11=',ihdgi
              call dmsread(nx,my,lncrec,'H',ifilin_ocaf,temp1(:,:),istat)
              write(ihdgi,12) k-1,ggdef,mm
              if(myrank .eq. 0) print *, 'ihdg12=',ihdgi
              call dmsread(nx,my,lncrec,'H',ifilin_ocaf,temp2(:,:),istat)
              if( myrank .eq. 72) then
                print *,"ocaf: mm=",mm,",wtfn(914,265)=",temp1(914,265)
                print *,"ocaf: mm=",mm,",wsfn(914,265)=",temp2(914,265)
              endif

              do jj=1,jlistnum
                j=jlist1(jj)
                nxj=nxdef_2d(j)
                if( lreduce.eq.1 )then
                  call reducepickr (temp1(1,j),nxdef(j),nx,1)
                  call reducepickr (temp2(1,j),nxdef(j),nx,1)
                endif
                do ii = 1, nxj
                  i=nxjstart(j)+ii-1
                  wtfn12(ii,k,jj,mm)=temp1(i,j)
                  wsfn12(ii,k,jj,mm)=temp2(i,j)

                  if(mm .eq. 1) then
                    wtfn12(ii,k,jj,13)=temp1(i,j)
                    wsfn12(ii,k,jj,13)=temp2(i,j)
                  endif
                  if(mm .eq. 12) then
                    wtfn12(ii,k,jj,0)=temp1(i,j)
                    wsfn12(ii,k,jj,0)=temp2(i,j)
                  endif
                enddo  !end do ii
              enddo    !end do jj
            enddo    !end do mm
          enddo !end of k

         deallocate(temp1,temp2)
         call dmscls(ifilin_ocaf,istat)
   
        END SUBROUTINE read_ocaf

!------------------------------------------------------------------------------
        SUBROUTINE read_ocaf0(nx,my,lkvl,ggdef,idtg1)

          use index
          use rank
          use const, only: ihdgi 

          integer nx,my,lkvl
          character*4 ggdef
          integer*8 idtg1
          character*12 cdtg
!          integer   iyyyy,mm,ddhhmn,iyy,imm,istat
          integer istat,k,lncrec
          integer j,jj,nxj,i,ii
          REAL, ALLOCATABLE, TARGET :: temp1(:,:),temp2(:,:)
          logical flag
!
          if(.not. ALLOCATED(temp1)) allocate (temp1(nx,my))
          if(.not. ALLOCATED(temp2)) allocate (temp2(nx,my))

          if(.not. ALLOCATED(wdepths)) ALLOCATE(wdepths(1:lkvl+2))
          if(.not. ALLOCATED(wtfn1st)) ALLOCATE(wtfn1st(nxp,1:lkvl+2,my_max))
          if(.not. ALLOCATED(mask1st)) ALLOCATE(mask1st(nxp,my_max))

          nwdepth=lkvl+2
          wdepths(1:lkvl+2)=sit_zdepth(0:lkvl+1)
          wtfn1st=0.
          mask1st=0.

          write(cdtg,'(i12)') idtg1
!          read(cdtg,'(i4,i2,i8)')iyyyy,mm,ddhhmn

#ifdef I38K
  11     format(i3.3,'TFN','000000',a4,a12)  ! ???TFM
  12     format(i3.3,'MSK','000000',a4,a12)  ! ???SFM
#else
  11     format(i3.3,'TFN','0000',a4,a12)  ! ???TFM
  12     format(i3.3,'MSK','0000',a4,a12)  ! ???SFM
#endif

!          do k=1,lkvl+2
          do k=1,1

           lncrec=nx*my
           write(ihdgi,11) k-1,ggdef,cdtg
           if(myrank .eq. 0) print *, 'ihdg11=',ihdgi
           call dmsread(nx,my,lncrec,'H',ifilin,temp1(:,:),istat)
           write(ihdgi,12) k-1,ggdef,cdtg
           if(myrank .eq. 0) print *, 'ihdg12=',ihdgi
           call dmsread(nx,my,lncrec,'H',ifilin,temp2(:,:),istat)

!           if( lreduce.eq.1 ) call reducepickr(temp1(1,1),nxdef,nx,my)
!           if( lreduce.eq.1 ) call reducepickr(temp2(1,1),nxdef,nx,my)


           do jj=1,jlistnum
              j=jlist1(jj)
              nxj=nxdef_2d(j)
              if( lreduce.eq.1 )then
                call reducepickr (temp1(1,j),nxdef(j),nx,1)
                call reducepickr (temp2(1,j),nxdef(j),nx,1)
              endif
              do ii = 1, nxj
                 i=nxjstart(j)+ii-1
                 wtfn1st(ii,k,jj)=temp1(i,j)
                 mask1st(ii,jj)=temp2(i,j)
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
          
          REAL, ALLOCATABLE, TARGET :: zin(:,:)
          REAL, ALLOCATABLE, TARGET :: zintemp(:,:)
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
          INTEGER :: lnc,zindepth
!ps
          ! Read world ocean atlas data file
          ! ===============

          flag=.false.
          IF (myrank .eq. 0) THEN
            WRITE(nerr,'(/,A,I2)') 'Read GODAS 1.0 '
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
              WRITE(nerr,'(/,A,I2)') 'Read GODAS data from unit ', ngodas
!              WRITE(nerr,'(/,A,I2)') '(http://www.cpc.ncep.noaa.gov/products/GODAS/) '
          
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
                WRITE(nerr,'(/,A,I2)') 'Read GODAS 2.0: open successfully ',woanc1%file_id
              ELSE
                WRITE (message_text,*) 'Could not open unit <',ngodas,'>'
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
!              WRITE(nerr,*) 'number of world ocean atlas data = ' &
!                            ,'io_ngl=',io_ngl,',io_nlon=',io_nlon &
!                            ,'io_ntime=',io_ntime

              IF (io_ngl/=ngl) THEN
                 WRITE(nerr,*) 'read_godas: unexpected resolution ',io_nlon,io_ngl
                 WRITE(nerr,*) 'expected number of latitudes = ',ngl
                 WRITE(nerr,*) 'number of latitudes of ocean data = ',io_ngl
                 WRITE(nerr,*)  'read_godas','unexpected resolution'
                 call dmsexit(-1)
                 stop
              END IF
              CALL io_inq_dimid  (woanc1%file_id, 'depth', io_var_id)
              CALL io_inq_dimlen (woanc1%file_id, io_var_id, nodepth)
              !WRITE(nerr,'(/,A,I2)') ' Read GODAS 4.1 '
              !WRITE(nerr,'(/,A,I2)') 'number of odepths = ',nodepth
              flag=.true.
            ELSE
              nodepth=nodepth0
            ENDIF

            IF(lamip .AND. .NOT. lmlo) THEN
              IF (lex0) THEN  
                CALL io_inq_dimid  (woanc0%file_id, 'depth', io_var_id)
                CALL io_inq_dimlen (woanc0%file_id, io_var_id, io_ndepth)
                IF (io_ndepth/=nodepth) THEN
                   WRITE(nerr,*) 'read_godas: inconsistent number of ' &
                                 ,'odepths between godas files',io_ndepth
                   WRITE(nerr,*) 'expected number of odepths = ',nodepth
                   WRITE (message_text,*) 'Read nodepth error in file <',fn0,'>'
                   WRITE(nerr,*) '',message_text
                   WRITE(nerr,*) 'read_godas','unexpected resolution'
                   call dmsexit(-1)
                   stop
                END IF
              ENDIF
              IF (lex2) THEN          
                CALL io_inq_dimid  (woanc2%file_id, 'depth', io_var_id)
                CALL io_inq_dimlen (woanc2%file_id, io_var_id, io_ndepth)
                IF (io_ndepth/=nodepth) THEN
                   WRITE(nerr,*) 'read_godas: inconsistent number of ' &
                                 ,'odepths between godas files',io_ndepth
                   WRITE(nerr,*) 'expected number of odepths = ',nodepth
                   WRITE (message_text,*) 'Read nodepth error in file <',fn2,'>'
                   WRITE(nerr,*) '',message_text
                   WRITE(nerr,*) 'read_godas','unexpected resolution'
                   call dmsexit(-1)
                   stop
                ENDIF
              ENDIF
            ENDIF
          END IF  !end IF (myrank .eq. 0) 
          !! RETURN
          call mpe_bcast(nodepth,1,0,mpe_integer)
          
!          WRITE(nerr,'(/,A,I2)') ' Read GODAS 5.0, nodepth= ',nodepth
          
          !     Allocate memory for ot12 and os12 per PE
          IF (.NOT. ALLOCATED(odepths)) ALLOCATE (odepths(nodepth))
          IF (.NOT. ALLOCATED(ot12)) ALLOCATE (ot12(nxp, nodepth, my_max, 0:13))
          IF (.NOT. ALLOCATED(os12)) ALLOCATE (os12(nxp, nodepth, my_max, 0:13))
          IF (.NOT. ALLOCATED(ou12)) ALLOCATE (ou12(nxp, nodepth, my_max, 0:13))
          IF (.NOT. ALLOCATED(ov12)) ALLOCATE (ov12(nxp, nodepth, my_max, 0:13))
          IF (.NOT. ALLOCATED(mixedlayer12)) ALLOCATE (mixedlayer12(nxp, my_max,0:13))
          
          ! Read odepths
          flag=.false.
          IF (myrank .eq. 0) THEN
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
          call mpe_bcast(odepths,nodepth,0,mpe_double)
          IF ( lwarning_msg.GE.3 ) THEN       
              WRITE (nerr,*) 'read_godas:: pe=',p_pe,',  odepths=',odepths(1:nodepth)
          ENDIF
          
          ! Codes read from GODAS
          cwoa(1) = 'ot'    ! ocean temperature profile (K)
          cwoa(2) = 'os'    ! ocean salinity profile (0/00)
          cwoa(3) = 'ou'    ! ocean u-component current (m/s)
          cwoa(4) = 'ov'    ! ocean v-component current (m/s)
          cwoa(5) = 'mixedlayer'    ! ocean mixed layer (m)
          IF (.NOT. ALLOCATED(zin)) ALLOCATE (zin(nlon,ngl))
          DO irec = 1, nrec
           !WRITE(nerr,'(/,A,I2)') ' Read GODAS 7.0 '
            zindepth=nodepth
            if(irec .eq. 5) zindepth=1
            cname = cwoa(irec)
            DO jk = 1, zindepth
              DO im=1, io_ntime
                flag=.false.
                IF (myrank .eq. 0) THEN
                  IF (.NOT. ALLOCATED(zintemp)) ALLOCATE (zintemp(nlon,ngl))
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
                      io_start(:) = (/       1,   1, jk, im /)
                      io_count(:) = (/ io_nlon, ngl,  1,  1 /)
                    ! for depth jk: read io_nlon longitudes, ngl latitudes and 12 months
                      CALL io_get_vara_double(woanc1%file_id,io_var_id,io_start,io_count, &
                                        zintemp(1:io_nlon,1:ngl))
                      zintemp=MERGE(zintemp,xmissing,zintemp.NE.missing_value)
                      Do j=1,ngl
                        zin(:,j)=zintemp(:,ngl-j+1)
                      ENDDO
                    ENDIF
                  ENDIF !end if(lex1)
                  flag=.true.
                ENDIF !end if(myrank.eq.0)
                call mpe_bcast(zin(:,:),nx*my,0,mpe_double)

                DO jj = 1, jlistnum
                  j=jlist1(jj)
                  nxj=nxdef_2d(j)
                  if( lreduce.eq.1 )call reducepickr(zin(1,j),nxdef(j),nx,1)
                  DO ii=1,nxj
                    i=nxjstart(j)+ii-1
                    IF(irec .eq. 1) THEN
                      ot12(ii,jk,jj,im) = zin(i,j)
                    ELSE IF(irec .eq. 2) THEN
                      os12(ii,jk,jj,im) = zin(i,j)
                    ELSE IF(irec .eq. 3) THEN
                      ou12(ii,jk,jj,im) = zin(i,j)
                    ELSE IF(irec .eq. 4) THEN
                      ov12(ii,jk,jj,im) = zin(i,j)
                    ELSE IF(irec .eq. 5) THEN
                      mixedlayer12(ii,jj,im) = zin(i,j)
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO  !end im (1-12)

              ! copy December to month 0
              ! copy January to month 13
              IF(irec .eq. 1) THEN
                ot12(:,jk,:,0) = ot12(:,jk,:,12)
                ot12(:,jk,:,13) = ot12(:,jk,:,1)
              ELSE IF(irec .eq. 2) THEN
                os12(:,jk,:,0) = os12(:,jk,:,12)
                os12(:,jk,:,13) = os12(:,jk,:,1)
              ELSE IF(irec .eq. 3) THEN
                ou12(:,jk,:,0) = ou12(:,jk,:,12)
                ou12(:,jk,:,13) = ou12(:,jk,:,1)
              ELSE IF(irec .eq. 4) THEN
                ov12(:,jk,:,0) = ov12(:,jk,:,12)
                ov12(:,jk,:,13) = ov12(:,jk,:,1)
              ELSE IF(irec .eq. 5) THEN
                mixedlayer12(:,:,0) = mixedlayer12(:,:,12)
                mixedlayer12(:,:,13)= mixedlayer12(:,:,1)
              ENDIF

  
              IF(lamip .AND. .NOT. lmlo) THEN
                flag=.false.
                IF(myrank .eq. 0) THEN
                  IF (lex0) THEN
                    status = NF_INQ_VARID (woanc0%file_id, cname, io_var_id)
                    IF (status /= NF_NOERR) THEN
                      WRITE(nerr,*) 'IO_INQ_VARID :', woanc0%file_id, cname
                      WRITE(nerr,*) 'IO_INQ_VARID', NF_STRERROR(status)
                    ELSE
                    ! read world ocean atlas data december of last year
                      CALL io_inq_varid (woanc0%file_id, cname, io_var_id)
                      io_start(:) = (/       1,   1, jk, 12 /)
                      io_count(:) = (/ io_nlon, ngl,  1,  1 /)
                      ! for depth jk: read io_nlon longitudes, ngl latitudes and 1 months
                      CALL io_get_vara_double(woanc0%file_id,io_var_id,io_start,io_count, &
                                              zintemp(1:io_nlon,1:io_ngl))
                      Do j=1,ngl
                        zin(:,j)=zintemp(:,ngl-j+1)
                      ENDDO
                    ENDIF
                  ENDIF  !end lex0
                  flag=.true.
                ENDIF !end (myrank .eq. 0)
                call mpe_bcast(zin(:,:),nx*my,0,mpe_double)

                DO jj = 1, jlistnum
                  j=jlist1(jj)
                  nxj=nxdef_2d(j)
                  if( lreduce.eq.1 )call reducepickr(zin(1,j),nxdef(j),nx,1)
                  DO ii=1,nxj
                    i=nxjstart(j)+ii-1
                    IF(irec .eq. 1) THEN
                      ot12(ii,jk,jj,0) = zin(i,j)
                    ELSE IF(irec .eq. 2) THEN
                      os12(ii,jk,jj,0) = zin(i,j)
                    ELSE IF(irec .eq. 3) THEN
                      ou12(ii,jk,jj,0) = zin(i,j)
                    ELSE IF(irec .eq. 4) THEN
                      ov12(ii,jk,jj,0) = zin(i,j)
                    ELSE IF(irec .eq. 5) THEN
                      mixedlayer12(ii,jj,0) = zin(i,j)
                    ENDIF
                  ENDDO
                ENDDO


                flag=.false.
                IF (myrank .eq.0 ) THEN
                  IF (lex2) THEN
                    status = NF_INQ_VARID (woanc2%file_id, cname, io_var_id)
                    IF (status /= NF_NOERR) THEN
                      WRITE(nerr,*) 'IO_INQ_VARID :', woanc2%file_id, cname
                      WRITE(nerr,*) 'IO_INQ_VARID', NF_STRERROR(status)
                    ELSE
                      ! read world ocean atlas data january of next year
                      CALL io_inq_varid (woanc2%file_id, cname, io_var_id)
                      io_start(:) = (/       1,   1, jk,  1 /)
                      io_count(:) = (/ io_nlon, ngl,  1,  1 /)
                      ! for depth jk: read io_nlon longitudes, ngl latitudes and 1 months
                      CALL io_get_vara_double(woanc2%file_id,io_var_id,io_start,io_count, &
                                           zintemp(1:io_nlon,1:io_ngl))
                      Do j=1,ngl
                        zin(:,j)=zintemp(:,ngl-j+1)
                      ENDDO
                    ENDIF
                  ENDIF  !end lex2
                  flag=.true.
                ENDIF
                call mpe_bcast(zin(:,:),nx*my,0,mpe_double)

                DO jj = 1, jlistnum
                  j=jlist1(jj)
                  nxj=nxdef_2d(j)
                  if( lreduce.eq.1 ) call reducepickr(zin(1,j),nxdef(j),nx,1)
                  DO ii=1,nxj
                    i=nxjstart(j)+ii-1
                    IF(irec .eq. 1) THEN
                      ot12(ii,jk,jj,13) = zin(i,j)
                    ELSE IF(irec .eq. 2) THEN
                      os12(ii,jk,jj,13) = zin(i,j)
                    ELSE IF(irec .eq. 3) THEN
                      ou12(ii,jk,jj,13) = zin(i,j)
                    ELSE IF(irec .eq. 4) THEN
                      ov12(ii,jk,jj,13) = zin(i,j)
                    ELSE IF(irec .eq. 5) THEN
                      mixedlayer12(ii,jj,13) = zin(i,j)
                    ENDIF
                  ENDDO
                ENDDO

              ENDIF  !end if(lamip .AND. .NOT. lmlo)
              !WRITE(nerr,'(/,A,I2)') ' Read GODAS 8.0 '
            END DO !end do jk=1,zindepth
          END DO  !end do irec=1, nrec
  
          IF (lwarning_msg.GE.3) THEN    
            WRITE(nerr,*) 'pe=',p_pe,   &
             ', read GODAS 9.1, os12(212,0,3,0:13)= ',os12(212,0,3,0:13)
          END IF
          
          IF (myrank .eq. 0) THEN
            IF (lex1) CALL io_close (woanc1)
            IF(lamip .AND. .NOT. lmlo) THEN
              IF (lex0) CALL io_close (woanc0)
              IF (lex2) CALL io_close (woanc2)       
            END IF
            DEALLOCATE (zintemp)
            WRITE(nerr,*) 'read GODAS'
          END IF
          
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

      REAL, ALLOCATABLE, TARGET :: zin(:,:)
      REAL, ALLOCATABLE, TARGET :: zintemp(:,:)
      REAL      :: missing_value

      INTEGER               :: io_nlon  ! number of longitudes in NetCDF file
      INTEGER               :: io_ngl   ! number of latitudes in NetCDF file
      INTEGER               :: io_ndepth  ! number of odepth0 in NetCDF file
      INTEGER, DIMENSION(4) :: io_start ! start index for NetCDF-read
      INTEGER, DIMENSION(4) :: io_count ! number of iterations for NetCDF-read

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
      INTEGER :: lnc,istat
      INTEGER :: zindepth
!ps
     ! Read world ocean atlas data file
     ! ===============
      flag=.false.
      IF (myrank .eq. 0) THEN
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
          WRITE (message_text,*) 'Could not open unit<',nwoa0,'>'
          WRITE(nerr,*) '',message_text
          WRITE(nerr,*) 'read_woa0', 'Could not open woa0 file'
          call dmsexit(-1)
          stop
        ENDIF
        WRITE(nerr,'(/,A,I2)') ' Read WOA0 3.0 '

      ! Check resolution
        CALL io_inq_dimid  (woa0nc1%file_id, 'lat', io_var_id)
        CALL io_inq_dimlen (woa0nc1%file_id, io_var_id, io_ngl)
        CALL io_inq_dimid  (woa0nc1%file_id, 'lon', io_var_id)
        CALL io_inq_dimlen (woa0nc1%file_id, io_var_id, io_nlon)
!        print *,"read W0A0 4.0"
        IF (io_ngl/=ngl) THEN
          WRITE(nerr,*) 'read_WOA0: unexpected resolution',io_nlon,io_ngl
          WRITE(nerr,*) 'expected number of latitudes = ',ngl
         WRITE(nerr,*) 'number of latitudes of world ocean atlas data =',io_ngl
          WRITE(nerr,*) 'read_WOA0','unexpected resolution'
          call dmsexit(-1)
          stop
        END IF
        CALL io_inq_dimid  (woa0nc1%file_id, 'depth', io_var_id)
        CALL io_inq_dimlen (woa0nc1%file_id, io_var_id, nodepth0)

        print *,"read W0A0 4.0"
        IF ( lwarning_msg.GE.3 ) THEN
          WRITE (nerr,*) 'read_WOA0:: pe=',p_pe,',odepth0=',odepth0(1:nodepth0)
        ENDIF
        flag=.true.
      END IF  !end if(myrank .eq. 0)
     ! RETURN
      CALL mpe_bcast(nodepth0,1,0,mpe_integer)
     ! Allocate memory for ot0, os0, ou0, ov0 per PE
      IF (.NOT. ALLOCATED(odepth0)) ALLOCATE (odepth0(nodepth0))
      IF (.NOT. ALLOCATED(ot0)) ALLOCATE (ot0(nxp, nodepth0, my_max))
      IF (.NOT. ALLOCATED(os0)) ALLOCATE (os0(nxp, nodepth0, my_max))
      IF (.NOT. ALLOCATED(ou0)) ALLOCATE (ou0(nxp, nodepth0, my_max))
      IF (.NOT. ALLOCATED(ov0)) ALLOCATE (ov0(nxp, nodepth0, my_max))
      IF (.NOT. ALLOCATED(mixedlayer0)) ALLOCATE (mixedlayer0(nxp, my_max))
      if(myrank.eq.0) print *,"W0A0, end allocate" 

      flag=.false.
      IF (myrank .eq. 0) THEN
        CALL io_inq_dimid  (woa0nc1%file_id, 'depth', io_var_id)
        CALL io_get_var_double (woa0nc1%file_id, io_var_id, odepth0)
!!      WRITE(nerr,'(/,A,I2)') ' Read WOA0 6.0 '
        flag=.true.
      ENDIF
      CALL mpe_bcast(odepth0,nodepth0,0,mpe_double)


     ! Codes read from WORLD OCEAN ATLAS 2005
      cwoa0( 1) = 'ot'    ! ocean temperature profile (K)
      cwoa0( 2) = 'os'    ! ocean salinity profile (0/00)
      cwoa0( 3) = 'ou'    ! ocean u-component current (m/s)
      cwoa0( 4) = 'ov'    ! ocean v-component current (m/s)
      cwoa0( 5) = 'mixedlayer'    ! ocean mixed layer (m)
      IF (.NOT. ALLOCATED(zin)) ALLOCATE (zin(nlon,ngl))

      DO irec = 1, nrec
        cname = cwoa0(irec)
        zindepth=nodepth0
        if(irec .eq. 5) zindepth=1
        DO jk = 1, zindepth
          flag=.false.
          IF (myrank .eq. 0) THEN
            IF (.NOT. ALLOCATED(zintemp)) ALLOCATE (zintemp(nlon,ngl))
            ! read world ocean atlas data
            status = NF_INQ_VARID (woa0nc1%file_id, cname, io_var_id)
            IF (status /= NF_NOERR) THEN
              WRITE(nerr,*) 'IO_INQ_VARID :', woa0nc1%file_id, cname
              WRITE(nerr,*) 'IO_INQ_VARID', NF_STRERROR(status)
            ELSE
              CALL io_get_att_double (woa0nc1%file_id, io_var_id, &
                                  '_FillValue', missing_value)

              io_start(:) = (/       1,  1, jk, 1 /)
              io_count(:) = (/ io_nlon, io_ngl ,  1, 1 /)
              ! for depth jk: read io_nlon longitudes, ngl latitudes and
              ! 12 months
              CALL io_get_vara_double(woa0nc1%file_id,io_var_id,io_start,io_count,&
                                      zintemp(1:io_nlon,1:io_ngl))
              print *,"woa0: irec=",irec,",jk=",jk,"zintemp:(768,277)=",zintemp(768,ngl+1-277)
            ENDIF
            zintemp(1:io_nlon,:)=MERGE(zintemp(1:io_nlon,:),xmissing,zintemp(1:io_nlon,:).NE.missing_value)
            Do j=1,ngl
              zin(:,j)=zintemp(:,ngl-j+1)
            ENDDO
            flag=.true.
          ENDIF

          call mpe_bcast(zin(:,:),nx*my,0,mpe_double)
          DO jj = 1, jlistnum
            j=jlist1(jj)
            nxj=nxdef_2d(j)
            if( lreduce.eq.1 )call reducepickr (zin(1,j),nxdef(j),nx,1)
            DO ii=1,nxj
              i=nxjstart(j)+ii-1
              IF (irec .eq. 1) THEN
                ot0(ii,jk,jj)=zin(i,j)
              ELSE IF(irec .eq. 2) THEN
                os0(ii,jk,jj)=zin(i,j)
              ELSE IF(irec .eq. 3) THEN
                ou0(ii,jk,jj)=zin(i,j)
              ELSE IF(irec .eq. 4) THEN
                ov0(ii,jk,jj)=zin(i,j)
              ELSE IF(irec .eq. 5) THEN
                mixedlayer0(ii,jj)=zin(i,j)
              ENDIF
            ENDDO
          ENDDO

        ENDDO  !end do jk
      END DO  !end do irec=1, nrec

      IF (myrank .eq. 0) THEN
        CALL io_close (woa0nc1)
        WRITE(nerr,*) 'end read WOA0'
        DEALLOCATE (zintemp)
      END IF
      DEALLOCATE (zin)

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
       if(tauleft > 3599.) then
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
       
       IF (myrank .eq. 0) THEN
!ps         CALL set_years(ihy0, ihy1, ihy2)
         ihy1=idtg_temp/100000000
         ihy0=ihy1-1
         ihy2=ihy1+1
         CALL chlen (ifilin_nc,80,lnc)
         WRITE (fn(1), '(A,A11,i4)') ifilin_nc(1:lnc),'/dailygodas',ihy0
         WRITE (fn(2), '(A,A11,i4)') ifilin_nc(1:lnc),'/dailygodas',yrori !ihy1
         WRITE (fn(3), '(A,A11,i4)') ifilin_nc(1:lnc),'/dailygodas',ihy2
!ps         WRITE (fn(1), '("dailygodas",i4)') ihy0
!ps         WRITE (fn(2), '("dailygodas",i4)') ihy1
!ps         WRITE (fn(3), '("dailygodas",i4)') ihy2
!ps         WRITE (nerr, '(/)')
!ps         WRITE(nerr,*) 'dailygodas ydate: ', ydate
!ps         WRITE(nerr,*) 'dailygodas ydate: fn(1): ' &
!ps           , TRIM(fn(1)),' fn(2): ', TRIM(fn(2)),' fn(3): ', TRIM(fn(3))
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
       IF (myrank .eq. 0) THEN
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
       ENDIF  !end if (myrank .eq. 0)

!       CALL mpe_broadcast (istat, 1, flag, mpe_integer)
       CALL mpe_bcast (istat, 1, 0, mpe_integer)
       if (myrank.eq.3) print *,'myrank3 istat=',istat
       if(istat .ne. 0)then
        if(myrank.eq.0) print *,'read_godas_3days fn2',' terminated.'
        call mpe_finalize
        call dmsexit(-1)
        stop
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
       IF (myrank .eq. 0) THEN
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
       IF (myrank .eq. 0) THEN
         CALL IO_INQ_DIMID (gpnc1%file_id, 'time', ndimid)
         CALL IO_INQ_DIMLEN (gpnc1%file_id, ndimid, nts1)
         CALL IO_INQ_VARID (gpnc1%file_id, 'time', io_var_id)
       
         IF ( nts1 .lt. 1 ) then
           WRITE(nerr,*) 'read_godas_3days', 'To few time steps < 1'
           istat=-1
         ELSE
           ALLOCATE (timevals1(nts1))
           CALL IO_GET_VAR_DOUBLE(gpnc1%file_id, io_var_id, timevals1)
           tsID=1
           DO WHILE ( (ydate.GT.INT(timevals1(tsID))).AND.(tsID.LT.nts1) )
             tsID=tsID+1
           ENDDO
           if(myrank .eq. 0) print*,"tsID=",tsID,",nts1=",nts1
!           IF ( tsID .gt. nts1 ) THEN
!             INQUIRE (file=fn(3), exist=lex2)
!             IF ( .NOT. lex2 ) THEN
!               WRITE (message_text,*) 'Could not open file <',fn(3),'>'
!               WRITE(nerr,*) message_text
!               WRITE(nerr,*) 'read_godas_3days fn3 ', 'run terminated.'
!               istat=-1
!             ELSE
!               CALL IO_open (fn(3), gpnc2, IO_READ)
!               CALL IO_INQ_DIMID (gpnc2%file_id, 'time', ndimid2)
!               CALL IO_INQ_DIMLEN (gpnc2%file_id, ndimid2, nts2)
!               IF ( nts2 .lt. 1 ) THEN
!                 WRITE (message_text,*) 'File <',fn(3),'>'
!                 WRITE (nerr,*) message_text
!                 WRITE (nerr,*) 'read_dailygodas:','To few time steps<1'
!                 istat=-1
!               ELSE
!                 ALLOCATE (timevals2(nts2))
!                 CALL IO_INQ_VARID (gpnc2%file_id, 'time', io_var_id)
!                 CALL IO_GET_VAR_DOUBLE(gpnc2%file_id, io_var_id, timevals2)
!                 tsID=1
!                 DO WHILE ( (ydate.GT.INT(timevals2(tsID))).AND.(tsID.LE.nts2) )
!                   tsID=tsID+1
!                 ENDDO
!                 IF ( tsID .gt. nts2 ) THEN
!                   WRITE(nerr,*) 'read_dailygodas', 'Date not found'
!                   istat=-1
!                 ELSE
!                   IF ( lwarning_msg.GE.1 ) THEN
!                     WRITE (nerr,*) 'read_godas_3days: timevals2_date=',timevals2(tsID)
!                   ENDIF  
!                   tsID_godas(2)=tsID
!                   files_godas(2)=3
!                   nts_godas(2)=nts2
!                 ENDIF
!               ENDIF
!               CALL IO_close(gpnc2)
!               IF (ALLOCATED(timevals2)) DEALLOCATE(timevals2)
!             ENDIF
!           ELSE
             tsID_godas(2)=tsID
             files_godas(2)=2
             nts_godas(2)=nts1
             IF ( lwarning_msg.GE.1 ) THEN
               WRITE (nerr,*) 'read_godas_3days: timevals1_date=',timevals1(tsID)
             ENDIF
!           ENDIF
         ENDIF
         CALL IO_close(gpnc1)
         IF (ALLOCATED(timevals1)) DEALLOCATE(timevals1)
         flag=.true.
       ENDIF

!       CALL mpe_broadcast (istat, 1, flag, mpe_integer)
       CALL mpe_bcast (istat, 1, 0, mpe_integer)
       if(istat .ne. 0)then
        if(myrank.eq.0) print *,'read_godas_3days.2', 'run terminated.'
        call mpe_finalize
        call dmsexit(-1)
        stop
       endif

       if(myrank .eq. 0) print*,'read to read_goads_dayp1' 
       !!! read day-1, day and day+1 data
       CALL read_godas_dayp1(THE_DATE)     ! day+0
       CALL read_godas_dayp1(DAY_MINUS1)   ! day-1
       CALL read_godas_dayp1(DAY_PLUS1)    ! day+1
       END SUBROUTINE read_godas_3days

    !-----------------------------------

       SUBROUTINE read_godas_dayp1(dayID)

       ! number of codes read from godas file
       INTEGER, PARAMETER :: nrec = 5

       INTEGER, INTENT(in):: dayID               ! dayID = DAY_MINUS1, THE_DAY_MINUS1, or DAY_PLUS1
       REAL, ALLOCATABLE, TARGET :: zin(:,:)
       REAL, ALLOCATABLE, TARGET :: zintemp(:,:)

       INTEGER jk,jj,j,nxj,i,ii
       INTEGER istat

       LOGICAL lex2
       INTEGER irec
       CHARACTER(12) :: cname, cgodas(nrec)
       INTEGER ndimid2,varid2
       REAL, ALLOCATABLE:: timevals2(:)
       REAL:: missing_value
       INTEGER zindepth
       INTEGER, DIMENSION(4) :: io_start ! start index for NetCDF-read
       INTEGER, DIMENSION(4) :: io_count ! number of iterations for NetCDF-read


       IF (dayID .EQ. DAY_PLUS1) THEN
       ! read Day+1 data
         IF (tsID_godas(2).EQ.nts_godas(2)) THEN
           files_godas(3)=files_godas(2)!+1
           tsID_godas(3)=nts_godas(2)         !first record of file 3
         ELSE
           files_godas(3)=files_godas(2)
           tsID_godas(3)=tsID_godas(2)+1
         ENDIF
       ELSEIF (dayID .EQ. DAY_MINUS1) THEN
       ! read Day-1 data
         IF (tsID_godas(2).EQ.1) THEN
           files_godas(1)=files_godas(2)!-1
           tsID_godas(1)=1 !LAST_RECORD     !last record of file 
         ELSE
           files_godas(1)=files_godas(2)
           tsID_godas(1)=tsID_godas(2)-1
         ENDIF
       ENDIF

       flag=.false.
       istat=-1
       IF (myrank .eq. 0) THEN
       INQUIRE (file=fn(files_godas(dayID)), exist=lex2)
       IF ( .NOT. lex2 ) THEN
         WRITE (message_text,*) 'Could not open file<',files_godas(dayID),'>'
         WRITE (nerr,*) message_text
         WRITE (nerr,*) 'read_dailygodas', 'run terminated.'
         istat=-1
       ELSE
         CALL IO_open (fn(files_godas(dayID)), gpnc2, IO_READ)
         CALL IO_INQ_DIMID (gpnc2%file_id, 'time', ndimid2)
         CALL IO_INQ_DIMLEN (gpnc2%file_id, ndimid2, nts_godas(dayID))
         IF ( nts_godas(dayID) .lt. 1 ) THEN
           WRITE (message_text,*) 'File <',fn(files_godas(dayID)),'>'
           WRITE(nerr,*) message_text
           WRITE(nerr,*) 'read_godas_dayp1:', 'To few time steps < 1'
           istat=-1
         ELSE
           IF(.NOT. ALLOCATED(timevals2)) ALLOCATE(timevals2(nts_godas(dayID)))
           CALL IO_INQ_VARID (gpnc2%file_id, 'time', io_var_id)
           CALL IO_GET_VAR_DOUBLE(gpnc2%file_id, io_var_id, timevals2)
           !!! modify tsID for LAST_RECORD
           IF (tsID_godas(dayID).EQ.LAST_RECORD) tsID_godas(dayID)=nts_godas(dayID)
           timevals_godas(dayID)=timevals2(tsID_godas(dayID))
           !!! data out range, set to be missing value
           IF (tsID_godas(dayID) .gt. nts_godas(dayID)) THEN
             WRITE (nerr,*) 'read_dailygodas,date=',timevals_godas(dayID) &
                           ,'DATA out of range!!! Set to MISSING DATA!'
             istat=-1
           ELSE
             istat=0
           ENDIF
           DEALLOCATE (timevals2)
         ENDIF
       ENDIF
       flag=.true.
       ENDIF   !!endif (myrank.eq.0)

!       CALL mpe_broadcast (istat, 1, flag, mpe_integer)
       CALL mpe_bcast (istat, 1, 0, mpe_integer)
       if(istat .ne. 0)then
         if(myrank.eq.0) print *,'read_dailygodas', 'run terminated.'
         call mpe_finalize
         call dmsexit(-1)
         stop
       endif

       CALL mpe_bcast(timevals_godas,3,0,mpe_double)
       CALL mpe_bcast(tsID_godas,3,0,mpe_integer)
       CALL mpe_bcast(files_godas,3,0,mpe_integer)
       CALL mpe_bcast(nts_godas,3,0,mpe_integer)

       IF(.NOT. ALLOCATED(zin)) ALLOCATE (zin(nx,my))
       ! Codes read from GODAS
       cgodas(1) = 'ot'    ! ocean temperature profile (K)
       cgodas(2) = 'os'    ! ocean salinity profile (0/00)
       cgodas(3) = 'ou'    ! ocean u-component current (m/s)
       cgodas(4) = 'ov'    ! ocean v-component current (m/s)
       cgodas(5) = 'mixedlayer'    ! ocean mixed layer (m)
       DO irec= 1, nrec
         cname=cgodas(irec)
         zindepth=nodepth
         if(irec.eq.5) zindepth=1
         DO jk=1, zindepth
           flag=.false.
           IF(myrank .eq. 0) then
             IF(.NOT. ALLOCATED(zintemp)) ALLOCATE (zintemp(nx,my))
             CALL IO_INQ_VARID (gpnc2%file_id, cname, varid2)
             CALL io_get_att_double (gpnc2%file_id,varid2,'_FillValue', missing_value)
             io_start(:) = (/ 1, 1, jk, tsID_godas(dayID) /)
             io_count(:) = (/ nlon, nlat, 1, 1 /)
             CALL IO_GET_VARA_DOUBLE(gpnc2%file_id,varid2,io_start,io_count,zintemp(1:nlon,1:nlat))
             zintemp(1:nlon,:)=MERGE(zintemp(1:nlon,:),xmissing,zintemp(1:nlon,:).NE.missing_value)
             Do j=1,nlat
               zin(:,j)=zintemp(:,nlat-j+1)
             ENDDO
             flag=.true.
           ENDIF   !end if (myrank .eq. 0)
           call mpe_bcast(zin(:,:),nx*my,0,mpe_double)

           DO jj = 1, jlistnum
             j=jlist1(jj)
             nxj=nxdef_2d(j)
             if( lreduce.eq.1 ) call reducepickr(zin(1,j),nxdef(j),nx,1)
             DO ii=1,nxj
               i=nxjstart(j)+ii-1
               if(irec .eq. 1) then
                 ot12(ii,jk,jj,dayID) = zin(i,j)
               else if(irec .eq. 2) then
                 os12(ii,jk,jj,dayID) = zin(i,j)
               else if(irec .eq. 3) then
                 ou12(ii,jk,jj,dayID) = zin(i,j)
               else if(irec .eq. 4) then
                 ov12(ii,jk,jj,dayID) = zin(i,j)
               else if(irec .eq. 5) then
                 mixedlayer12(ii,jj,dayID) = zin(i,j)
               endif
             ENDDO
           ENDDO  !!end jj
         ENDDO  !!end jk
       ENDDO   !!end irec

       IF(myrank .eq. 0) THEN
         CALL IO_close(gpnc2)
         DEALLOCATE (zintemp)
       ENDIF

       DEALLOCATE (zin)

       END SUBROUTINE read_godas_dayp1

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
!      IF ( lwarning_msg.GE.2 ) then 
!        WRITE(nerr,'(/,A,I5,A,I9,A,I9,A)') ' fill_missing: k=',k,&
!             ',   ',nmissing,' missing value grids, where ',nfilled,&
!             ' grids filled.'
!      ENDIF
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
               IF (ydate_jd .EQ. ydate2_jd) THEN
                 wgto2=1
               else
                 wgto2=(ydate_jd-ydate1_jd)/(ydate2_jd-ydate1_jd)
               endif
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
               CALL mpe_finalize
               CALL dmsexit(-1)
               stop
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
        
!          if(myrank .eq. 49) then
!            print *,'time weight, idtg1=',idtg1,',yr=',yr,',mo=',mo &
!                   ,',dy=',dy,',hr=',hr,',tauleft=',tauleft  & 
!                   ,',ydate=',ydate,',ydate1=',ydate1        &
!                   ,',ydate2=',ydate2,',ydatejd=',ydate_jd   &
!                  ,',ydate1jd=',ydate1_jd,',ydate2jd=',ydate2_jd &
!                  ,',obswtbwgt1=',obswtbwgt1   &
!                  ,',obswtbwgt2=',obswtbwgt2
!          endif

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
             if(myrank==0) then
             print*,'mo_time_weight:Get_JulianMonLen, month invalid' &
                   , 'ky=',ky,'km=',km
             endif

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

         if(myrank==0)then
         IF (kd < 1 .OR. idmax < kd) &
           print *,'Set_JulianDay: day in months invalid' &
                  ,' kd=', kd,'idmax= ',idmax
         endif

         zd = real(365.25*iy)+INT(30.6001*(im+1)) &
             +REAL(ib)+1720996.5+REAL(kd)+zsec

         END FUNCTION Set_JulianDay


        END SUBROUTINE time_weights







      END MODULE mod_sst
