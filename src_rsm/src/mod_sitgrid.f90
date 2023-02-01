      module mod_sitgrid
!
! modify to f90 by C-H Lee and sort by River Chen in 2015
!
!      USE param
!      USE index
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
              dtswdt, tseadiffSIT, tseadiffSIT24,        &
              sumdSITdt, countdSITdt, ratioSIT

 
      contains 

        subroutine allocate_sitgrid_array(nxp,my_max)

           integer  ierr
           integer  nxp,my_max

!    ! 2-d from ATM vars
           allocate (                                                  &
         sitcor(nxp,my_max),       slm(nxp,my_max), sitlclass(nxp,my_max),&
        sitmask(nxp,my_max),  sitmask2(nxp,my_max),     bathy(nxp,my_max),&
           wlvl(nxp,my_max),   ocnmask(nxp,my_max), obox_mask(nxp,my_max),&
            sni(nxp,my_max),       tsi(nxp,my_max),      tsw(nxp,my_max), &
            tsl(nxp,my_max),      tslm(nxp,my_max),    tslm1(nxp,my_max), &
            ocu(nxp,my_max),       ocv(nxp,my_max), ctfreez2(nxp,my_max), &
                stat=ierr)
           if (ierr/= 0) then
             write(6,*) 'mod_sitgrid_2d : allocate fail 1 '
             stop
           end if

!    ! 2-d SIT vars
           allocate (                                                  &
          sitwtb(nxp,my_max),  sitwub(nxp,my_max), sitwvb(nxp,my_max),    &
          sitwsb(nxp,my_max),  fluxiw(nxp,my_max),   pme2(nxp,my_max),    &
        subfluxw(nxp,my_max), wsubsal(nxp,my_max),                       &
           sitcc(nxp,my_max),   sithc(nxp,my_max), engwac(nxp,my_max),    &
              sc(nxp,my_max), saltwac(nxp,my_max),                       &
           wtfns(nxp,my_max),   wsfns(nxp,my_max), stat=ierr)
          if (ierr/= 0) then
            write(6,*) 'mod_sitgrid_2d : allocate fail 2 '
            stop
          end if

!    ! 4- OUTPUT only, original ATM variabels
          allocate (                                                   &
           seaice(nxp,my_max), grndcapc(nxp,my_max),                     &
         grndhflx(nxp,my_max), grndflux(nxp,my_max), stat=ierr)
          if (ierr/= 0) then
            write(6,*) 'mod_sitgrid_2d : allocate fail 3 '
            stop
          end if

!    ! other
          allocate(                                                    &
            fluxw(nxp,my_max),    dfluxs(nxp,my_max),  soflw(nxp,my_max), &
            fluxi(nxp,my_max),     sofli(nxp,my_max),                    &
        thickness(nxp,my_max), obsseaice(nxp,my_max), obswtb(nxp,my_max), &
           obswsb(nxp,my_max),       rsf(nxp,my_max),    ssf(nxp,my_max), &
            disch(nxp,my_max),     ustrw(nxp,my_max),  vstrw(nxp,my_max), &
            evapw(nxp,my_max),   wind10w(nxp,my_max), stat=ierr)

           if (ierr/= 0) then
               write(6,*) 'mod_sitgrid_2d : allocate fail 4 '
               stop
           end if
   
   
!     & ! 3-d SIT vars: snow/ice
           allocate (                                                  &
          zsi(nxp,my_max,0:1), silw(nxp,my_max,0:1),                     &
        tsnic(nxp,my_max,0:3), stat=ierr)
           if (ierr/= 0) then
               write(6,*) 'mod_sitgrid_3d : allocate fail 1 '
               stop
           end if

!     & ! 3-d SIT vars: water column (part1)
           allocate (                                                  &
               obswt(nxp,my_max,0:lkvl+1),   obsws(nxp,my_max,0:lkvl+1), &
               obswu(nxp,my_max,0:lkvl+1),   obswv(nxp,my_max,0:lkvl+1), &
               sitwt(nxp,my_max,0:lkvl+1),   sitwu(nxp,my_max,0:lkvl+1), &
               sitwv(nxp,my_max,0:lkvl+1),   sitww(nxp,my_max,0:lkvl+1), &
               sitws(nxp,my_max,0:lkvl+1), sitwtke(nxp,my_max,0:lkvl+1), &
                wlmx(nxp,my_max,0:lkvl+1),  wldisp(nxp,my_max,0:lkvl+1), &
                 wkm(nxp,my_max,0:lkvl+1),     wkh(nxp,my_max,0:lkvl+1), &
            wrho1000(nxp,my_max,0:lkvl+1),sftobswt(nxp,my_max,0:lkvl+1), &
             stat=ierr)
           if (ierr/= 0) then
               write(6,*) 'mod_sitgrid_3d : allocate fail 2 '
               stop
           end if

!     & ! 3-d SIT vars: water column (part2)
           allocate (                                                  &
                 wtfn(nxp,my_max,0:lkvl+1),  wsfn(nxp,my_max,0:lkvl+1),  &
                wtfn0(nxp,my_max,0:lkvl+1), wsfn0(nxp,my_max,0:lkvl+1),  &
                awufl(nxp,my_max,0:lkvl+1), awvfl(nxp,my_max,0:lkvl+1),  &
                awtfl(nxp,my_max,0:lkvl+1), awsfl(nxp,my_max,0:lkvl+1),  &
               awtfl0(nxp,my_max,0:lkvl+1),awsfl0(nxp,my_max,0:lkvl+1),  &
              awtkefl(nxp,my_max,0:lkvl+1), stat=ierr)

           if (ierr/= 0) then
               write(6,*) 'mod_sitgrid_3d : allocate fail 3 '
               stop
           end if

! for SIT var. every tau mean
           allocate (         &
             sitwttau(nxp,my_max,0:lkvl+1), sitwstau(nxp,my_max,0:lkvl+1), &
             sitwutau(nxp,my_max,0:lkvl+1), sitwvtau(nxp,my_max,0:lkvl+1), &
             stat=ierr)

! for SIT var. daily mean
           allocate (         &
             sitwt24(nxp,my_max,0:lkvl+1), sitws24(nxp,my_max,0:lkvl+1), &
             sitwu24(nxp,my_max,0:lkvl+1), sitwv24(nxp,my_max,0:lkvl+1), &
             stat=ierr)

           if (ierr/= 0) then
               write(6,*) 'mod_sitgrid_sit24 : allocate fail 1 '
               stop
           end if


! for store vars. which sit_vdiff need during "fsit" period
           allocate(                                                   &
          ssfsit(nxp,my_max),    rsfsit(nxp,my_max), hfluxfsit(nxp,my_max), &
       qfluxfsit(nxp,my_max),   u10fsit(nxp,my_max),   v10fsit(nxp,my_max), &
        rlspfsit(nxp,my_max),  rcupfsit(nxp,my_max), ustarfsit(nxp,my_max), &
          t2fsit(nxp,my_max),   rh2fsit(nxp,my_max),   pstfsit(nxp,my_max), &
        cicefsit(nxp,my_max),   snrfsit(nxp,my_max),  zicefsit(nxp,my_max), &
       xticefsit(nxp,my_max),obswtbfsit(nxp,my_max),    tgfsit(nxp,my_max), &
             stat=ierr)

          if (ierr/= 0) then
               write(6,*) 'mod_sitgrid_varfsit : allocate fail 1 '
               stop
           end if

! for store latest 2 steps tsw
!          allocate( oldtsw(nx,my_max,0:1), stat=ierr)
          allocate(                                                    &
     oldsitwt(nxp,my_max,0:lkvl+1,0:1),oldsitwu(nxp,my_max,0:lkvl+1,0:1),&
     oldsitwv(nxp,my_max,0:lkvl+1,0:1),oldsitww(nxp,my_max,0:lkvl+1,0:1),&
     oldsitws(nxp,my_max,0:lkvl+1,0:1),oldsitwtke(nxp,my_max,0:lkvl+1,0:1), &
          stat=ierr)

          if (ierr/= 0) then
               write(6,*) 'mod_sitgrid_oldsit : allocate fail 1 '
               stop
           end if

           allocate ( dtswdt(nxp,my_max), tseadiffSIT(nxp,my_max),     &
                     tseadiffSIT24(nxp,my_max), sumdSITdt(nxp,my_max), &
                     countdSITdt(nxp,my_max), ratioSIT(nxp,my_max), stat=ierr)

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
           dtswdt=0.
           tseadiffSIT=0.
           tseadiffSIT24=0.
           sumdSITdt=0.
           countdSITdt=0.
           ratioSIT=0.
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

           deallocate ( dtswdt, tseadiffSIT, tseadiffSIT24, &
                        sumdSITdt, countdSITdt, ratioSIT)

           return

        end subroutine deallocate_sitgrid_array

      end module mod_sitgrid
