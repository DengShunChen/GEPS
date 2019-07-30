      module mod_sitgrid
!
! modify to f90 by C-H Lee and sort by River Chen in 2015
!
!      USE param
      USE mod_sit_control,     ONLY: xmissing
      USE mo_netcdf,         ONLY: lkvl
  
    
      implicit none

      public
!      LOGICAL, SAVE :: do_sit     = .TRUE. ! .true. for calculation of upper ocean temperature profile using sit model

!    ! 2-d from ATM vars
       real, dimension(:,:),allocatable,save ::          &
              sitcor,       slm, sitlclass,              &
             sitmask,  sitmask2, bathy,  wlvl,           &
             ocnmask, obox_mask,                         &
                 sni,       tsi,       tsw,              &
                 tsl,      tslm,     tslm1,              &
                 ocu,       ocv,  ctfreez2

!    ! 2-d SIT vars
       real, dimension(:,:),allocatable,save ::          &
              sitwtb,    sitwub,    sitwvb,              &
              sitwsb,    fluxiw,      pme2,              &
            subfluxw,   wsubsal,                         & 
               sitcc,     sithc,    engwac,              &
                  sc,   saltwac,                         &
               wtfns,     wsfns

!    ! 4- OUTPUT only, original ATM variabels
       real, dimension(:,:),allocatable,save ::          &
              seaice,  grndcapc,                         &
            grndhflx,  grndflux

!    ! other
       real, dimension(:,:),allocatable,save ::          &
               fluxw,    dfluxs,     soflw,              &
               fluxi,    sofli,                          &
           thickness, obsseaice,    obswtb,              &
              obswsb,       rsf,       ssf,              &
               disch,     ustrw,     vstrw,              &
               evapw,   wind10w

   
!     & ! 3-d SIT vars: snow/ice
      real, dimension(:,:,:),allocatable,save ::         &
                 zsi,  silw,    tsnic
!    ! 3-d SIT vars: water column
      real, dimension(:,:,:),allocatable,save ::         &
               obswt,   obsws,    obswu,  obswv,         &
               sitwt,   sitwu,    sitwv,  sitww,         &   
               sitws, sitwtke,     wlmx, wldisp,         &
                 wkm,   wkh, wrho1000, sftobswt
      real, dimension(:,:,:),allocatable,save ::         &
                wtfn,   wsfn,  wtfn0,  wsfn0,             & 
               awufl,  awvfl,  awtfl, awsfl,              &
              awtfl0, awsfl0, awtkefl
      real:: dtsitmon

! for SIT var. every tau mean
      real:: dtsittau
      real, dimension(:,:,:),allocatable,save ::         &
            sitwttau, sitwstau, sitwutau, sitwvtau

! for SIT var. daily mean
      real:: dtsit24
      real, dimension(:,:,:),allocatable,save ::         &
             sitwt24, sitws24, sitwu24, sitwv24
  
  
! for store vars. which sit_vdiff need during "fsit" period 
      real dtfsit          ! accumulate dtx when do varfsit=varfsit+var*dtx
      real,dimension(:,:), allocatable, save::           &
              ssfsit,  rsfsit, hfluxfsit, qfluxfsit,     &
             u10fsit, v10fsit,  rlspfsit,  rcupfsit,     &
           ustarfsit,  t2fsit,   rh2fsit,   pstfsit,     &
            cicefsit, snrfsit,  zicefsit, xticefsit,     &
          obswtbfsit,  tgfsit
            
! for store latest 2 steps sitwt
      real,dimension(:,:,:,:), allocatable, save::       &
            oldsitwt, oldsitwu, oldsitwv, oldsitww,      &
            oldsitws, oldsitwtke
     
      real,dimension(:,:), allocatable, save::           &
              tgold, dtswdt

 
      contains 

        subroutine allocate_sitgrid_array(nx,my_max)

           integer  ierr
           integer  nx,my_max

!    ! 2-d from ATM vars
           allocate (                                                  &
         sitcor(nx,my_max),       slm(nx,my_max), sitlclass(nx,my_max),&
        sitmask(nx,my_max),  sitmask2(nx,my_max),     bathy(nx,my_max),&
           wlvl(nx,my_max),   ocnmask(nx,my_max), obox_mask(nx,my_max),&
            sni(nx,my_max),       tsi(nx,my_max),      tsw(nx,my_max), &
            tsl(nx,my_max),      tslm(nx,my_max),    tslm1(nx,my_max), &
            ocu(nx,my_max),       ocv(nx,my_max), ctfreez2(nx,my_max), &
                stat=ierr)
           if (ierr/= 0) then
             write(6,*) 'mod_sitgrid_2d : allocate fail 1 '
             stop
           end if

!    ! 2-d SIT vars
           allocate (                                                  &
          sitwtb(nx,my_max),  sitwub(nx,my_max), sitwvb(nx,my_max),    &
          sitwsb(nx,my_max),  fluxiw(nx,my_max),   pme2(nx,my_max),    &
        subfluxw(nx,my_max), wsubsal(nx,my_max),                       &
           sitcc(nx,my_max),   sithc(nx,my_max), engwac(nx,my_max),    &
              sc(nx,my_max), saltwac(nx,my_max),                       &
           wtfns(nx,my_max),   wsfns(nx,my_max), stat=ierr)
          if (ierr/= 0) then
            write(6,*) 'mod_sitgrid_2d : allocate fail 2 '
            stop
          end if

!    ! 4- OUTPUT only, original ATM variabels
          allocate (                                                   &
           seaice(nx,my_max), grndcapc(nx,my_max),                     &
         grndhflx(nx,my_max), grndflux(nx,my_max), stat=ierr)
          if (ierr/= 0) then
            write(6,*) 'mod_sitgrid_2d : allocate fail 3 '
            stop
          end if

!    ! other
          allocate(                                                    &
            fluxw(nx,my_max),    dfluxs(nx,my_max),  soflw(nx,my_max), &
            fluxi(nx,my_max),     sofli(nx,my_max),                    &
        thickness(nx,my_max), obsseaice(nx,my_max), obswtb(nx,my_max), &
           obswsb(nx,my_max),       rsf(nx,my_max),    ssf(nx,my_max), &
            disch(nx,my_max),     ustrw(nx,my_max),  vstrw(nx,my_max), &
            evapw(nx,my_max),   wind10w(nx,my_max), stat=ierr)

           if (ierr/= 0) then
               write(6,*) 'mod_sitgrid_2d : allocate fail 4 '
               stop
           end if
   
   
!     & ! 3-d SIT vars: snow/ice
           allocate (                                                  &
          zsi(nx,my_max,0:1), silw(nx,my_max,0:1),                     &
        tsnic(nx,my_max,0:3), stat=ierr)
           if (ierr/= 0) then
               write(6,*) 'mod_sitgrid_3d : allocate fail 1 '
               stop
           end if

!     & ! 3-d SIT vars: water column (part1)
           allocate (                                                  &
               obswt(nx,my_max,0:lkvl+1),   obsws(nx,my_max,0:lkvl+1), &
               obswu(nx,my_max,0:lkvl+1),   obswv(nx,my_max,0:lkvl+1), &
               sitwt(nx,my_max,0:lkvl+1),   sitwu(nx,my_max,0:lkvl+1), &
               sitwv(nx,my_max,0:lkvl+1),   sitww(nx,my_max,0:lkvl+1), &
               sitws(nx,my_max,0:lkvl+1), sitwtke(nx,my_max,0:lkvl+1), &
                wlmx(nx,my_max,0:lkvl+1),  wldisp(nx,my_max,0:lkvl+1), &
                 wkm(nx,my_max,0:lkvl+1),     wkh(nx,my_max,0:lkvl+1), &
            wrho1000(nx,my_max,0:lkvl+1),sftobswt(nx,my_max,0:lkvl+1), &
             stat=ierr)
           if (ierr/= 0) then
               write(6,*) 'mod_sitgrid_3d : allocate fail 2 '
               stop
           end if

!     & ! 3-d SIT vars: water column (part2)
           allocate (                                                  &
                 wtfn(nx,my_max,0:lkvl+1),  wsfn(nx,my_max,0:lkvl+1),  &
                wtfn0(nx,my_max,0:lkvl+1), wsfn0(nx,my_max,0:lkvl+1),  &
                awufl(nx,my_max,0:lkvl+1), awvfl(nx,my_max,0:lkvl+1),  &
                awtfl(nx,my_max,0:lkvl+1), awsfl(nx,my_max,0:lkvl+1),  &
               awtfl0(nx,my_max,0:lkvl+1),awsfl0(nx,my_max,0:lkvl+1),  &
              awtkefl(nx,my_max,0:lkvl+1), stat=ierr)

           if (ierr/= 0) then
               write(6,*) 'mod_sitgrid_3d : allocate fail 3 '
               stop
           end if

! for SIT var. every tau mean
           allocate (         &
             sitwttau(nx,my_max,0:lkvl+1), sitwstau(nx,my_max,0:lkvl+1), &
             sitwutau(nx,my_max,0:lkvl+1), sitwvtau(nx,my_max,0:lkvl+1), &
             stat=ierr)

! for SIT var. daily mean
           allocate (         &
             sitwt24(nx,my_max,0:lkvl+1), sitws24(nx,my_max,0:lkvl+1), &
             sitwu24(nx,my_max,0:lkvl+1), sitwv24(nx,my_max,0:lkvl+1), &
             stat=ierr)

           if (ierr/= 0) then
               write(6,*) 'mod_sitgrid_sit24 : allocate fail 1 '
               stop
           end if


! for store vars. which sit_vdiff need during "fsit" period
           allocate(                                                   &
          ssfsit(nx,my_max),    rsfsit(nx,my_max), hfluxfsit(nx,my_max), &
       qfluxfsit(nx,my_max),   u10fsit(nx,my_max),   v10fsit(nx,my_max), &
        rlspfsit(nx,my_max),  rcupfsit(nx,my_max), ustarfsit(nx,my_max), &
          t2fsit(nx,my_max),   rh2fsit(nx,my_max),   pstfsit(nx,my_max), &
        cicefsit(nx,my_max),   snrfsit(nx,my_max),  zicefsit(nx,my_max), &
       xticefsit(nx,my_max),obswtbfsit(nx,my_max),    tgfsit(nx,my_max), &
             stat=ierr)

          if (ierr/= 0) then
               write(6,*) 'mod_sitgrid_varfsit : allocate fail 1 '
               stop
           end if

! for store latest 2 steps tsw
!          allocate( oldtsw(nx,my_max,0:1), stat=ierr)
          allocate(                                                    &
     oldsitwt(nx,my_max,0:lkvl+1,0:1),oldsitwu(nx,my_max,0:lkvl+1,0:1),&
     oldsitwv(nx,my_max,0:lkvl+1,0:1),oldsitww(nx,my_max,0:lkvl+1,0:1),&
     oldsitws(nx,my_max,0:lkvl+1,0:1),oldsitwtke(nx,my_max,0:lkvl+1,0:1), &
          stat=ierr)

          if (ierr/= 0) then
               write(6,*) 'mod_sitgrid_oldsit : allocate fail 1 '
               stop
           end if

           allocate ( tgold(nx,my_max), dtswdt(nx,my_max), stat=ierr)

          if (ierr/= 0) then
               write(6,*) 'mod_sitgrid_dtswdt : allocate fail 1 '
               stop
           end if


           sitcor=xmissing
           sitlclass=xmissing
           sitmask=xmissing
           sitmask2=xmissing
           bathy=xmissing
           wlvl=xmissing
           ocnmask=xmissing
           obox_mask=xmissing
           sni=xmissing
           tsi=xmissing
           tsw=xmissing
           tsl=xmissing
           tslm=xmissing
           tslm1=xmissing
           ocu=xmissing
           ocv=xmissing
           slm=xmissing                       
           ctfreez2=xmissing
!     !     2-d SIT vars
           sitwtb=xmissing
           sitwub=xmissing
           sitwvb=xmissing
           sitwsb=xmissing
           fluxiw=xmissing
           pme2=xmissing
           subfluxw=xmissing
           wsubsal=xmissing
           sitcc=xmissing
           sithc=xmissing
           engwac=xmissing
           sc=xmissing
           saltwac=xmissing
           wtfns=xmissing
           wsfns=xmissing
!     !     3-d SIT vars: snow/ice
           zsi=0.
           silw=5.
           tsnic=xmissing
!     !     3-d SIT vars: water column
           obswt=xmissing
           obsws=xmissing
           obswu=xmissing
           obswv=xmissing
           sitwt=288.
           sitwu=0.
           sitwv=0.
           sitww=0.
           sitws=0.
           sitwtke=0.
           wlmx=xmissing
           wldisp=xmissing
           wkm=0.
           wkh=0.
           wrho1000=xmissing
           sftobswt=0.
           dtfsit=0.
           dtsitmon=0. 
           wtfn=0.
           wsfn=0.                   
           wtfn0=0.
           wsfn0=0.                   
           awufl=0.
           awvfl=0.
           awtfl=0.
           awsfl=0.
           awtfl0=0.
           awsfl0=0.
           awtkefl=0.
!     !     4- OUTPUT only, original ATM variabels
           seaice=xmissing                                             
           grndcapc=xmissing
           grndhflx=xmissing
           grndflux=xmissing
!     !    other
           fluxw=0.
           dfluxs=0.
           fluxi=0.
           soflw=xmissing
           sofli=xmissing
           thickness=xmissing
           obsseaice=xmissing
           obswtb=xmissing
           obswsb=xmissing
           rsf=0.
           ssf=0.
           disch=0.
           ustrw=xmissing
           vstrw=xmissing
           evapw=xmissing
           wind10w=xmissing
! for SIT var. daily mean
           dtsittau=0.
           sitwttau=0.
           sitwstau=0.
           sitwutau=0.
           sitwvtau=0.

! for SIT var. daily mean
           dtsit24=0.
           sitwt24=0.
           sitws24=0.
           sitwu24=0.
           sitwv24=0.

! for store vars. which sit_vdiff need during "fsit" period
           ssfsit=0.
           rsfsit=0.
           hfluxfsit=0.
           qfluxfsit=0.
           u10fsit=0.
           v10fsit=0.
           rlspfsit=0.
           rcupfsit=0.
           ustarfsit=0.
           t2fsit=0.
           rh2fsit=0.
           pstfsit=0.
           cicefsit=0.
           snrfsit=0.
           zicefsit=0.
           xticefsit=0.
           obswtbfsit=0.
           tgfsit=0.

! for store latest 2 steps tsw
!           oldtsw=xmissing
           oldsitwt=xmissing
           oldsitwu=xmissing
           oldsitwv=xmissing
           oldsitww=xmissing
           oldsitws=xmissing
           oldsitwtke=xmissing
           tgold=0.
           dtswdt=0.

           return

        end subroutine allocate_sitgrid_array

 
 
        subroutine deallocate_sitgrid_array

!    ! 1-input only, original ATM/SIT variabels
           deallocate (                                  &
              sitcor,       slm, sitlclass,              &
             sitmask,  sitmask2,     bathy,      wlvl,   &
             ocnmask, obox_mask,                         &
                 sni,       tsi,       tsw,              &
                 tsl,      tslm,     tslm1,              &
                 ocu,       ocv,  ctfreez2)
!    ! 2-d SIT vars
           deallocate (                                  &
              sitwtb,    sitwub,    sitwvb,              &
              sitwsb,    fluxiw,      pme2,              &
            subfluxw,   wsubsal,                         & 
               sitcc,     sithc,    engwac,              &
                  sc,   saltwac,                         &
               wtfns,     wsfns)

!    ! 4- OUTPUT only, original ATM variabels
           deallocate (                                  &
              seaice,  grndcapc,                         &
            grndhflx,  grndflux)

!    ! other
           deallocate (                                  &
               fluxw,    dfluxs,     soflw,              &
               fluxi,    sofli,                          &
           thickness, obsseaice,    obswtb,              &
              obswsb,       rsf,       ssf,              &
               disch,     ustrw,     vstrw,              &
               evapw,   wind10w)

!    ! 3-d SIT vars: water column
           deallocate (                                  &
                 zsi,  silw,    tsnic)
           deallocate (                                  &
               obswt,   obsws,  obswu,  obswv,           &
               sitwt,   sitwu,  sitwv,  sitww,           &   
               sitws, sitwtke,   wlmx, wldisp,           &
                 wkm,   wkh, wrho1000, sftobswt)
           deallocate (                                  &
                wtfn,   wsfn, wtfn0, wsfn0,              &
               awufl,  awvfl, awtfl, awsfl,              &
              awtfl0, awsfl0, awtkefl )

! for SIT var. every tau mean
           deallocate (                                  &
             sitwttau, sitwstau, sitwutau, sitwvtau) 

! for SIT var. daily mean
           deallocate (                                  &
             sitwt24, sitws24, sitwu24, sitwv24)

! for store vars. which sit_vdiff need during "fsit" period
           deallocate (                                  &
              ssfsit,  rsfsit, hfluxfsit, qfluxfsit,     &
             u10fsit, v10fsit,  rlspfsit,  rcupfsit,     &
           ustarfsit,  t2fsit,   rh2fsit,   pstfsit,     &
            cicefsit, snrfsit,  zicefsit, xticefsit,     &
          obswtbfsit,  tgfsit)

! for store latest 2 steps tsw
           deallocate (                                  &
            oldsitwt, oldsitwu, oldsitwv,oldsitww,       &
            oldsitws,oldsitwtke)

           deallocate ( tgold, dtswdt)

           return

        end subroutine deallocate_sitgrid_array

      end module mod_sitgrid
