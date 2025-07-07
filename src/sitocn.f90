subroutine sitocn(j,jj,nxj, dt ,dta ,tau ,idtg ,xlon ,xlat,   &
           cp,grav,hltm,lrun_sitvdiff ,ic_sit,                &
           ss_adj,rs_adj ,hflux,qflux,ustar,u10,v10,t2,rh2,   &
           rcup,rlsp,pst,cice,snr,zice,xtice,ocean,itimestep, &
           tg)

use param,only:nx,my,my_max
use mpe
use rank
use index             !nxp ,nxjstart ,nxjp ,jlistnum
use const,                 ONLY:do_sit,ldailyFCTsst,dailyClm_option,      &
                                pdfcloud,cmbk,cgwd, fsit, dosppt, doshum, dossst, &
                                use_zmtnblck,ldailyFCTicesndpt,dSITdt_intv, &
                                weightSIT,bckfile,ggdef,doclx,doslavepp,    &
                                RTYPE,qmin,updatetg
use mod_sitgrid
USE mod_sit_vdiff,         ONLY:sit_vdiff,ctfreez
USE mod_sit_control,       ONLY:ftrigsit,ltrigsit,lsitstart,lsftobswt &
                               ,timebl_option,timebl_start
USE mod_eos_ocean,         ONLY:AirVaporPressure,CalcSm
USE mod_sst,               ONLY:time_weights,now1,now2,wgto1,wgto2 &
                               ,obswtbnmw1,obswtbnmw2,obswtbwgt1   &
                               ,obswtbwgt2,obswtbold,obswtbnow     &
                               ,obswtbnew,dFCTsstdt                &
                               ,tseadiffFCT,tseadiffFCT24,dtseadt  &
                               ,tseanow,tseaold,tseanew
USE mo_netcdf,             ONLY:lkvl
use mod_stochastic_physics, only : ssst3d
implicit  none
!input
integer::j,jj,nxj
integer::itimestep
real::dt ,dta ,tau ,cp,grav ,hltm
real::xlat(my), xlon(nx,my_max)
integer(kind=8)  :: idtg
logical lrun_sitvdiff
integer ic_sit
real::ss_adj(nxp),rs_adj(nxp) , &
      hflux(nxp,my_max),qflux(nxp,my_max) ,ustar(nxp,my_max), &
      u10(nxp,my_max),v10(nxp,my_max) ,   &
      t2(nxp,my_max),rh2(nxp,my_max),     &
      rcup(nxp,my_max),rlsp(nxp,my_max),  &
      cice(nxp,my_max),snr(nxp,my_max),zice(nxp,my_max),xtice(nxp,my_max)
logical:: ocean(nxp,my_max)
real(kind=RTYPE)::pst(nxp,my_max)
! in/out
real::tg(nxp,my_max)
!local variable
      real liw, t_surf
      real sitlat(nxp)
      real sitlon(nxp,my_max)
      integer istep,icurrenttau
      real rnuw,z0w,d0w,Pairvapor,rhoa,ps1w,ustra,vstra
real::dtsit ,tauleft 
integer*8:: idtg_sitvdiff
character:: cdtg*12
integer yr, mo, dy, hr, mn
integer::i,ii
real tauhr ,dtaup,dtx_tau
real, dimension(nxp) :: sstm, rstm,hfluxtm,qfluxtm, &
                                    u10tm,v10tm, rlsptm, rcuptm, &
                                  ustartm, t2tm,  rh2tm,  psttm, &
                                   cicetm,snrtm, zicetm,xticetm, &
                                 obswtbtm, tgtm

INTEGER, PARAMETER :: nerr = 6
!======start
!if(myrank==0)then
!  print*, j,jj,nxj, dt ,dta ,tau ,idtg 
!  print*,cp,grav,hltm,lrun_sitvdiff ,ic_sit 
!endif
dtx_tau = dt/3600.

  dtsit=dt
  istep= int(tau/(dt/3600.)+0.01)
  tauleft=float(int((tau-int(tau)+0.001)*3600./dt))*dt
  icurrenttau=int(tau)
  if(tauleft > 3590.0   )then
     icurrenttau=icurrenttau+1
     tauleft=0.0
  endif
  call dtgfix12(idtg,idtg_sitvdiff,icurrenttau)
  call time_weights(idtg_sitvdiff,tauleft)
  write(cdtg,'(i12)')idtg_sitvdiff
  read(cdtg,'(i4,i2,i2,i2,i2)')yr,mo,dy,hr,mn
!  ssttau = mod(tau+0.001, 24.)
  tauhr=float(hr)+tauleft/3600.

!
!
  if (jj .eq. 1) then
    dtfsit=dtfsit+dtsit
  endif

  do ii = 1, nxj
    i=nxjstart(j)+ii-1
    sitlat(ii)=xlat(j)
    IF(xlon(i,jj).LT. 0.) THEN
      sitlon(ii,jj)=xlon(i,jj)+360.
    ELSE
      sitlon(ii,jj)=xlon(i,jj)
    ENDIF

    if (.NOT. lrun_sitvdiff .or. ic_sit .eq. 1) then
      ssfsit(ii,jj)    =ssfsit(ii,jj)    +ss_adj(ii)      *dtsit
      rsfsit(ii,jj)    =rsfsit(ii,jj)    +rs_adj(ii)      *dtsit
      hfluxfsit(ii,jj) =hfluxfsit(ii,jj) +hflux(ii,jj)    *dtsit
      qfluxfsit(ii,jj) =qfluxfsit(ii,jj) +qflux(ii,jj)    *dtsit
      u10fsit(ii,jj)   =u10fsit(ii,jj)   +u10(ii,jj)      *dtsit
      v10fsit(ii,jj)   =v10fsit(ii,jj)   +v10(ii,jj)      *dtsit
      rlspfsit(ii,jj)  =rlspfsit(ii,jj)  +rlsp(ii,jj)     *dtsit
      rcupfsit(ii,jj)  =rcupfsit(ii,jj)  +rcup(ii,jj)     *dtsit
      ustarfsit(ii,jj) =ustarfsit(ii,jj) +ustar(ii,jj)    *dtsit
      t2fsit(ii,jj)    =t2fsit(ii,jj)    +t2(ii,jj)       *dtsit
      rh2fsit(ii,jj)   =rh2fsit(ii,jj)   +rh2(ii,jj)      *dtsit
      pstfsit(ii,jj)   =pstfsit(ii,jj)   +pst(ii,jj)      *dtsit
      cicefsit(ii,jj)  =cicefsit(ii,jj)  +cice(ii,jj)     *dtsit
      snrfsit(ii,jj)   =snrfsit(ii,jj)   +snr(ii,jj)      *dtsit
      zicefsit(ii,jj)  =zicefsit(ii,jj)  +zice(ii,jj)     *dtsit
      xticefsit(ii,jj) =xticefsit(ii,jj) +xtice(ii,jj)    *dtsit
      obswtbfsit(ii,jj)=obswtbfsit(ii,jj)+obswtbnow(ii,jj)*dtsit
      tgfsit(ii,jj)    =tgfsit(ii,jj)    +tg(ii,jj)       *dtsit
    endif

    if(lrun_sitvdiff)then
      if(ic_sit .eq. 1) then
        sstm(ii)    =ssfsit(ii,jj)    /dtfsit
        rstm(ii)    =rsfsit(ii,jj)    /dtfsit
        hfluxtm(ii) =hfluxfsit(ii,jj) /dtfsit
        qfluxtm(ii) =qfluxfsit(ii,jj) /dtfsit
        u10tm(ii)   =u10fsit(ii,jj)   /dtfsit
        v10tm(ii)   =v10fsit(ii,jj)   /dtfsit
        rlsptm(ii)  =rlspfsit(ii,jj)  /dtfsit
        rcuptm(ii)  =rcupfsit(ii,jj)  /dtfsit
        ustartm(ii) =ustarfsit(ii,jj) /dtfsit
        t2tm(ii)    =t2fsit(ii,jj)    /dtfsit
        rh2tm(ii)   =rh2fsit(ii,jj)   /dtfsit
        psttm(ii)   =pstfsit(ii,jj)   /dtfsit
        cicetm(ii)  =cicefsit(ii,jj)  /dtfsit
        snrtm(ii)   =snrfsit(ii,jj)   /dtfsit
        zicetm(ii)  =zicefsit(ii,jj)  /dtfsit
        xticetm(ii) =xticefsit(ii,jj) /dtfsit
        obswtbtm(ii)=obswtbfsit(ii,jj)/dtfsit
        tgtm(ii)    =tgfsit(ii,jj)    /dtfsit

        ssfsit(ii,jj)    =0.
        rsfsit(ii,jj)    =0.
        hfluxfsit(ii,jj) =0.
        qfluxfsit(ii,jj) =0.
        u10fsit(ii,jj)   =0.
        v10fsit(ii,jj)   =0.
        rlspfsit(ii,jj)  =0.
        rcupfsit(ii,jj)  =0.
        ustarfsit(ii,jj) =0.
        t2fsit(ii,jj)    =0.
        rh2fsit(ii,jj)   =0.
        pstfsit(ii,jj)   =0.
        cicefsit(ii,jj)  =0.
        snrfsit(ii,jj)   =0.
        zicefsit(ii,jj)  =0.
        xticefsit(ii,jj) =0.
        obswtbfsit(ii,jj)=0.
        tgfsit(ii,jj)    =0.
      else
        sstm(ii)=ss_adj(ii)
        rstm(ii)=rs_adj(ii)
        hfluxtm(ii)=hflux(ii,jj)
        qfluxtm(ii)=qflux(ii,jj)
        u10tm(ii)=u10(ii,jj)
        v10tm(ii)=v10(ii,jj)
        rlsptm(ii)=rlsp(ii,jj)
        rcuptm(ii)=rcup(ii,jj)
        ustartm(ii)=ustar(ii,jj)
        t2tm(ii)=t2(ii,jj)
        rh2tm(ii)=rh2(ii,jj)
        psttm(ii)=pst(ii,jj)
        cicetm(ii)=cice(ii,jj)
        snrtm(ii)=snr(ii,jj)
        zicetm(ii)=zice(ii,jj)
        xticetm(ii)=xtice(ii,jj)
        obswtbtm(ii)=obswtbnow(ii,jj)
        tgtm(ii)=tg(ii,jj)
        dtfsit=dtsit
      endif

      fluxw(ii,jj)=-(1.-cicetm(ii))*(sstm(ii)-rstm(ii)-hfluxtm(ii)-qfluxtm(ii))
      fluxi(ii,jj)=-cicetm(ii)*(sstm(ii)-rstm(ii)-hfluxtm(ii)-qfluxtm(ii))
      soflw(ii,jj)=(1.-cicetm(ii))*sstm(ii)
      sofli(ii,jj)=cicetm(ii)*sstm(ii)
      wind10w(ii,jj)=sqrt(u10tm(ii)**2+v10tm(ii)**2)
      obsseaice(ii,jj)=cicetm(ii)
      evapw(ii,jj)=-qfluxtm(ii)/hltm
      rsf(ii,jj)=(rlsptm(ii)+rcuptm(ii))/dta

      !ustrw, vstrw
      rnuw = 1.14E-6
      z0w = 0.11*rnuw/ustartm(ii)+0.1*SQRT(rnuw*ustartm(ii)/grav) &
             +0.015*ustartm(ii)**2/grav  ! Kraus & Businger(1994, page 146)
      d0w = 0.

      Pairvapor= AirVaporPressure(t2tm(ii),rh2tm(ii))
      rhoa = 100.*psttm(ii)/(287.04*t2tm(ii))*(1.0-0.378*Pairvapor/psttm(ii))
           ! 287.04 [k.g-1.K-1]  Gas constant of dry air
      liw = -(0.4*grav*((hfluxtm(ii)/(t2tm(ii)*cp)+0.61*(qfluxtm(ii)/hltm)))/(ustartm(ii)**3.)*rhoa)
      ps1w = log((10.0-d0w)/z0w)-CalcSm((10.0-d0w)*liw,z0w,liw)+CalcSm(z0w*liw,z0w,liw)
      ustra = u10tm(ii)*0.4/ps1w
      vstra = v10tm(ii)*0.4/ps1w
      ustrw(ii,jj) = rhoa*ustra**2
      vstrw(ii,jj) = rhoa*vstra**2

      !in&out
      seaice(ii,jj)=cicetm(ii)
      thickness(ii,jj)=zicetm(ii)
      sni(ii,jj)=snrtm(ii)/1000.   !mm->m
      t_surf   =obswtbtm(ii)
      t_surf = MAX( t_surf, ctfreez )
      tsi(ii,jj)= MIN( t_surf, ctfreez )
      obswtb(ii,jj)=obswtbtm(ii)
      tsw(ii,jj)=tg(ii,jj)
      dtswdt(ii,jj)=0.

    endif    !end (lrun_sitvdiff)
  enddo    !end i=1,nxj

        if(lrun_sitvdiff) then
         call sit_vdiff ( nxjp(j), nxp, jj, istep, dtfsit,             &
              sitlat, sitlon(:,jj), tau, tauhr,                        &
              sitcor(:,jj), slm(:,jj), sitlclass(:,jj),                &
              sitmask(:,jj), bathy(:,jj), wlvl(:,jj),   &
              ocnmask(:,jj), obox_mask(:,jj),                          &
!            ! - same as lake and ml_ocean
              fluxw(:,jj), dfluxs(:,jj), soflw(:,jj),                  &
              fluxi(:,jj), sofli(:,jj),                                &
!            ! - 1D from mo_memory_g3b (wind stress)
              ustrw(:,jj), vstrw(:,jj),                                &
!            ! -    water mass variables(rain, snow, evap and runoff):
              rsf(:,jj), ssf(:,jj), evapw(:,jj), disch(:,jj),          &
!            ! - 1D from mo_memory_g3b
!              t2(:,j), wind10w(:,jj),                                 &
              t2tm(:), wind10w(:,jj),                               &
!            ! - 1D from mo_memory_g3b (sit variables)
!              obsseaice(:,jj), obswtbtm(:,jj), obswsb(:,jj),           &
              obsseaice(:,jj), obswtb(:,jj), obswsb(:,jj),           &
              sitwtb(:,jj), sitwub(:,jj), sitwvb(:,jj),                &
              sitwsb(:,jj), fluxiw(:,jj), pme2(:,jj),                  &
              subfluxw(:,jj), wsubsal(:,jj),                           &
              sitcc(:,jj), sithc(:,jj), engwac(:,jj),                  &
              sc(:,jj), saltwac(:,jj), wtfns(:,jj), wsfns(:,jj),       &
!            ! 3-d SIT vars: snow/ice
              zsi(:,jj,0:1), silw(:,jj,0:1), tsnic(:,jj,0:3),          &
!            ! 3-d SIT vars: water column
              obswt(:,jj,0:lkvl+1),obsws(:,jj,0:lkvl+1),               &
              obswu(:,jj,0:lkvl+1),obswv(:,jj,0:lkvl+1),               &
              sitwt(:,jj,0:lkvl+1),sitwu(:,jj,0:lkvl+1),               &
              sitwv(:,jj,0:lkvl+1),sitww(:,jj,0:lkvl+1),               &
              sitws(:,jj,0:lkvl+1),sitwtke(:,jj,0:lkvl+1),             &
              wlmx(:,jj,0:lkvl+1),wldisp(:,jj,0:lkvl+1),               &
              wkm(:,jj,0:lkvl+1),wkh(:,jj,0:lkvl+1),                   &
              wrho1000(:,jj,0:lkvl+1),wtfn(:,jj,0:lkvl+1),             &
              wsfn(:,jj,0:lkvl+1), wtfn0(:,jj,0:lkvl+1),               &
              wsfn0(:,jj,0:lkvl+1),awufl(:,jj,0:lkvl+1),               &
              awvfl(:,jj,0:lkvl+1),awtfl(:,jj,0:lkvl+1),               &
              awsfl(:,jj,0:lkvl+1),awtfl0(:,jj,0:lkvl+1),              &
              awsfl0(:,jj,0:lkvl+1),awtkefl(:,jj,0:lkvl+1),            &
!             ! final output only
!              cice(:,j), snr(:,j), zice(:,j), xtice(:,j), tg(:,j),    &
              seaice(:,jj),sni(:,jj),thickness(:,jj),xticetm(:),    &
              tsw(:,jj), tsl(:,jj), tslm(:,jj), tslm1(:,jj),           &
              ocu(:,jj), ocv(:,jj),  ctfreez2(:,jj),                   &
!            ! implicit with vdiff
              grndcapc(:,jj), grndhflx(:,jj), grndflux(:,jj),          &
              oldsitwt(:,jj,0:lkvl+1,0:1),oldsitwu(:,jj,0:lkvl+1,0:1), &
              oldsitwv(:,jj,0:lkvl+1,0:1),oldsitww(:,jj,0:lkvl+1,0:1), &
             oldsitws(:,jj,0:lkvl+1,0:1),oldsitwtke(:,jj,0:lkvl+1,0:1),&
              dtswdt(:,jj), sftobswt(:,jj,0:lkvl+1) )
        
          if(jj .eq. 1) then
            dtsittau=dtsittau+dtfsit
            dtsitmon=dtsitmon+dtfsit
            dtsit24 =dtsit24 +dtfsit
          endif
          call storesittau(nxjp(j),jj,nxp,my_max,lkvl,sitwt,sitws,sitwu,sitwv,dtfsit)
          call  storesit24(nxjp(j),jj,nxp,my_max,lkvl,sitwt,sitws,sitwu,sitwv,dtfsit)

        endif  !end lrun_sitvdiff
!
        if(jj .eq. jlistnum) then
          if (lrun_sitvdiff) dtfsit=0.
        endif
!=======================================================================
! Cal. SST tendency (dtsea/dt)
!=======================================================================
       if(ldailyFCTsst .OR. ldailyFCTicesndpt .OR. dailyClm_option.ge.1) then
        if(itimestep .le. 1) then
          do ii = 1,nxj
            i=nxjstart(j)+ii-1
            dtseadt(ii,jj)=0.
            dFCTsstdt(ii,jj)=0. 
            obswtbnow(ii,jj)=dta*dFCTsstdt(ii,jj)+obswtbold(ii,jj)
            if(ocean(ii,jj))then
              tseadiffFCT(ii,jj)=dta*dFCTsstdt(ii,jj)
              tseanow(ii,jj)=dta*dFCTsstdt(ii,jj)+ tseaold(ii,jj)
            endif
          end do
        endif
        if (tau .ge. 24.) then
           do ii = 1, nxj
              i=nxjstart(j)+ii-1
              dtseadt(ii,jj)    =0.
              tseadiffFCT(ii,jj)=0.
              if (ocean(ii, jj)) then
                obswtbnew(ii,jj)=dta*dFCTsstdt(ii,jj)+obswtbold(ii,jj)
                obswtbold(ii,jj)=obswtbnew(ii,jj)
                obswtbnow(ii,jj)=obswtbnew(ii,jj)
                ! sea surface temperature tendency 
                dtseadt(ii,jj)=dFCTsstdt(ii,jj)
                tseadiffFCT(ii,jj)=dta*dFCTsstdt(ii,jj)
                if(do_sit) then
                  if(sitmask(ii,jj) .EQ. 1.) then
                    if(lrun_sitvdiff .AND. ltrigsit )then
                      tseadiffSIT(ii,jj)=0.
                      sumdSITdt(ii,jj)=sumdSITdt(ii,jj)+dtswdt(ii,jj)
                      countdSITdt(ii,jj)=countdSITdt(ii,jj)+1.
                    endif
                   
                    dtaup = mod(tau+0.001, dSITdt_intv)
                    if( (dSITdt_intv .lt. 0.) .OR. (dtaup .lt. dtx_tau) )then
                      if(countdSITdt(ii,jj) .ge. 1.) then
                        tseadiffFCT(ii,jj)=dta * (1.-weightSIT*ratioSIT(ii,jj)) * dFCTsstdt(ii,jj)
                        tseadiffSIT(ii,jj)=dta * weightSIT*ratioSIT(ii,jj) * & 
                                                (sumdSITdt(ii,jj)/countdSITdt(ii,jj))
                        tseadiffSIT24(ii,jj)=tseadiffSIT24(ii,jj)+tseadiffSIT(ii,jj)/dta*dt
                        ! sea surface temperature tendency
                        dtseadt(ii,jj)=(1.-weightSIT*ratioSIT(ii,jj)) * dFCTsstdt(ii,jj) &
                                    + weightSIT * ratioSIT(ii,jj) * &
                                      (sumdSITdt(ii,jj)/countdSITdt(ii,jj))
                        sumdSITdt(ii,jj)=0.
                        countdSITdt(ii,jj)=0.
                      endif
                    endif
                  endif ! end if sitmask(ii,jj) .EQ. 1
                endif ! end if do_sit
                if (dossst) then 
                  dtseadt(ii,jj) = dtseadt(ii,jj) * (ssst3d(ii,1,jj) + 1.)
                endif 
                tseadiffFCT24(ii,jj)=tseadiffFCT24(ii,jj)+tseadiffFCT(ii,jj)/dta*dt
              end if !end if(ocean)
           end do
        end if
       end if

       if(ldailyFCTsst .OR. ldailyFCTicesndpt .OR. dailyClm_option.ge.1) then
         if (tau .ge. 24.) then
           do ii = 1, nxj
              i=nxjstart(j)+ii-1
              if (ocean(ii, jj)) then
                !
                ! update tg, dtsea/dt (W00100)
                !
                tseanew(ii, jj) = dta*dtseadt(ii, jj) + tseaold(ii, jj)
                !tseaold(ii,jj)=tseanow(ii,jj) + tfilt*(tseaold(ii,jj)     &
                !          -2.0*tseanow(ii,jj)+tseanew(ii,jj) )
                tseaold(ii, jj) = tseanew(ii, jj)
                tseanow(ii, jj) = tseanew(ii, jj)

                dtaup = mod(tau + 0.001, updatetg)
                if (dtaup .lt. dtx_tau) then
                   tg(ii, jj) = tseanow(ii, jj)
                end if
              endif
           enddo
         endif
       endif

end subroutine
