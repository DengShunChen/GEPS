#define myrank_check 55
#define ii_check 1
#define jj_check 4

       subroutine intgrt_3tl
!
!***********************************************************************
!  this subroutine is the basic time stepping driver.  it does the
!  gaussian quadrature integrations to compute the adiabatic spectral
!  tendencies for the model variables. these tendencies are adjusted
!  with a semi-implicit method to make the model stable for long
!  time steps.  time steps are made with a time filter to remove
!  any computational modes.
!
!  modify to f90 by C-H Lee and sort by River Chen in 2015
!
!***********************************************************************
!
      use param
      use mpe
      use rank
      use index
      use const
      use spec
      use grid
      use fftcom
      use phygrid
      use mod_typhoon
      use noah
      use namelist_soilveg
!-----------------------------------------------------------------------
      USE mod_sitgrid
      USE mod_sit_vdiff,       ONLY:sit_vdiff_end
      USE mod_sit_control,     ONLY:xmissing,lgodas,locaf,lwoa0   &
                                    ,loutsit24,outsitmean,lsitstart   &
                                    ,ldailysst
      USE mod_sst,             ONLY:deallocate_ocaf_array,deallocate_woa0_array &
                                 ,deallocate_godas_array,read_dailygodas        &
                                 ,time_weights,nmw1,nmw2,wgt1,wgt2              &
                                 ,read_dailyFCT,dailyFCTsst,dailyFCTcice        &
                                 ,dailyFCTsndepth,deallocate_dailyFCT_array     &
                                 ,ANAsstT0,dailyClmANAsst,dailyClmFCTsst        &
                                 ,allocate_opgsst_array,read_opgsst             &
                                 ,opgsst,deallocate_opgsst_array                &
                                 ,obswtbnmw1,obswtbnmw2,obswtbwgt1,obswtbwgt2   &
                                 ,sit_zdepth
      USE mo_netcdf,           ONLY:lkvl,cleanup_netcdf
!-----------------------------------------------------------------------
      use raddiag
      use radn
      use albn
!-----------------------------------------------------------------------
      use mod_stochastic_physics, only : spptout,skebout, &
                  run_stochastic_physics,                 &
                  destroy_stochastic_physics,             &
                  skeb3du,skeb3dv,diss_est,skebfilt,      &
                  keb,kea
!
!-----------------------------------------------------------------------

      implicit  none
!
!  local work array
!
      integer   nfxr
!  for Semi-Lagrangian
      real(kind=RTYPE) ndsldta,ndsldtah,                           &
                deldm(nxp,my_max),ddtemp(nxp,lev,my_max),          &
                pten(nxp,lev,my_max),pdot(nxp,lev+1,latpart),      &
                qvadv(nxp,lev*ncld,my_max),diveng(nxp,lev,my_max), &
                vdmerd(nxp,lev,my_max),vdzonl(nxp,lev,my_max),     &
                vdmerdr(nxp,lev,my_max),vdzonlr(nxp,lev,my_max),   &
                ttm_sl(nx,levp,my_max),                            &
                pten_sl(nx,levp,my_max),                           &
                qm_sl(nx,levp*ncld,my_max),                        &
                vvm_sl(nx,levp,my_max),                            &
                uum_sl(nx,levp,my_max),ddtempr(nxp,lev,my_max)
!
      integer   ierr,itt,year,yrd
!
      real(kind=RTYPE) glob(nx,my)
      real      hf24(nxp,my_max),qf24(nxp,my_max),ss24(nxp,my_max),rs24(nxp,my_max), &
                asol24(nxp,my_max),olr24(nxp,my_max),rain24(nxp,my_max),             &
                drag(nxp,lev,my_max),ugws(nxp,my_max),vgws(nxp,my_max),              &
                sdpbl(nxp,my_max),rain1(nxp,my_max),                                 &
                pklev(nxp,my_max),totallp(nxp,my_max),rainlp24(nxp,my_max)

       integer  kn

      real      tmin(nxp,my_max),tmax(nxp,my_max),td(nxp,my_max),temp
!

      real(kind=RTYPE) pltemp(jtrun,jtmax,2),cc(nx+2,levp,1,my_max)
      real      ww1(nx,my_max)
!byl      real      cc3(nx+2,levp,3,my_max),wss3(levp,2,3,jtrun,jtmax)
!
      character rfile*55, ctau*6
!!      real      tbar(lev),qbar(lev*ncld),qbrrow(ncld,my)
      integer,  parameter :: ktop=4
      real      fac(ktop), wkj(my,4), wkmf(jtrun), windmax3, sumtot  &
               ,sumtoti  ,sumtott(nxp,my_max)  ,qtot(nxp,my_max)     &
               ,sumwati  ,sumwat   ,sumwatt(nxp,my_max)
      data      windmax3/130./
!
      logical   forward, histim, tchange, flag
!
      logical   wrestrt
      data      wrestrt/.false./
!
! for topographic gravity wave drag
!
      real hprime_b(nxp,mtnvar,my_max)

! for due point temperature
      real tda,tdb
      parameter (tda=17.27,tdb=237.7)

!
!   restart  : write(7) work array
!
!!      real, dimension(:), allocatable :: work_io
!
!  restart  : write(10) work array
!
      real, dimension(:,:,:), allocatable :: tm1,tm2,tm3,tm4
!
      real, dimension(:,:,:), allocatable :: tmc1,tmc2
      real, dimension(:,:,:), allocatable :: tmc3,tmc4
      real, dimension(:,:,:), allocatable :: tmc5,tmc6
      real, dimension(:,:,:), allocatable :: temp1,temp2,temp3
!-------------------------------------------------------------------
      real, dimension(:,:,:), allocatable :: tmr1,tmr2
      real, dimension(:,:,:), allocatable :: tmr3,tmr4
      real, dimension(:,:,:), allocatable :: tmr5,tmr6
      real, dimension(:,:,:), allocatable :: tmr7,tmr8
      real, dimension(:,:,:), allocatable :: tms1,tms2,tms3,tms4
!-------------------------------------------------------------------
!
      integer i,j,k,m,n,jj,kk,mf,kw,nxj,nml,lmax,leng,nxmy, &
              jlim,mlst,mlmax2,itaui,itaue,itauo,itaup,     &
              ntau,itau,lcwb,lphy,ifromtau,itotau,istat,    &
              istst,ii,nxjf !,n_stable,n_unstable,nc_stable

      real    www,dtx,dta,thdai,tkei,tpei,dsigp,            &
              cosw,tengi,dt24,tg2,dtx_tau,sqhaf,            &
              dummy,dt1,sptend,wmax,xx,facw,dtaup!!,          &
!!              sptendmax2,sptendmax1,dt_chg
! sppt variables
!            by John Tseng 2017/12/13 
!            modified by PangYen Liu for 2D-MPI 2019/02/20
!      real sppt2d(nx,my)
!      real rold500(mlmax_c,2),rold1000(mlmax_c,2),rold2000(mlmax_c,2)
!      real rold500(jtrun_c,jtmax_c,2),rold1000(jtrun_c,jtmax_c,2),  &
!           rold2000(jtrun_c,jtmax_c,2)
      integer recn

! for io quilting
      character*34 keydoit,keydone
      data keydoit/"DOIT..........................DOIT"/
      data keydone/"DONE..........................DONE"/
!
!for sst_restore_tau>0., update sst(W00100), seaice(W00091), snowdepth(B00650)
      integer*8 idtg_sst,idtg1_sst,idtg_temp
      integer icurrenttau,yyyymmdd,hhii
      logical lsstrestore,iceold(nxp,my),oceanold(nxp,my)
      character*12 cdtg
      real    ssttemp,cicetemp,snrtemp
      real    sst(nx,my),ssttau,tautemp
      integer yr, mo, dy, hr, mn

!for opgsst sst
      integer*8 idtg1_opgsst
      integer icurrentyear,inextyear,inexttau
      logical lnewyear
      real tauleft

!for SIT
      integer icurrentyymmdd,inextyymmdd
      logical lnewday
      integer icurrentyymm, ibeforeyymm
      logical lnewyymm
      integer ic_sit,nc_sit
      logical turn_sit,lrun_sitvdiff
      real wweight
      integer kkk
!
!xb110>
!byl      real rmr(nxp,lev,my_max),smr(nxp,lev,my_max)
!for lightning scheme from ECMWF
      real flash(nxp,my_max),flash24(nxp,my_max)
!xb110<

#ifdef TIMING
! for timing
      real*8 tm_1,tm_2,tm_use,mpi_wtime
!CWB2020
      tm_1=mpi_wtime()
      tm_2=mpi_wtime()
#endif
!      fsit=-99.             !fsit>0., turn on sit_vdiff when
!      mod(tau/fsit)<0.001
                             !default fsit<=0., turn on sit_vdiff every tau
      ic_sit=-99             !if fsit>0., store now tau is the ic_sit times when sit_vidff is turn on
      nc_sit=1               !if fsit>0., when mod(tau/fsit)<0.001, turn on sit_vdiff for "nc_sit" timesteps
      turn_sit=.false.       !turn_sit=.true., will run sit_vdiff in some tau
      lrun_sitvdiff=.false.  !lrun_sitvdiff=.true., run sit_vdiff in this tau
!
      lmax  = 16
      nfxr  = 33
!
      year = idate(1)
      yrd  = 365
      if ( mod(year,4) .eq. 0 ) yrd = 366
!
! output initialization field
!
      if(wrestrt)then
!
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
!-------------------------------------------------------------------
        allocate (tmr1(nx,lev+1,my))
        allocate (tmr2(nx,lev+1,my))
        allocate (tmr3(nx,lev+1,my))
        allocate (tmr4(nx,lev+1,my))
        allocate (tmr5(nx,lev+1,my))
        allocate (tmr6(nx,lev+1,my))
        allocate (tmr7(nx,lev+1,my))
        allocate (tmr8(nx,lev+1,my))
        allocate (tms1(nx,lev,my))
        allocate (tms2(nx,lev,my))
        allocate (tms3(nx,lev,my))
        allocate (tms4(nx,lev,my))
!-------------------------------------------------------------------
!        call unify_gridr                                                  &
!          ( tmr1,tmr2,tmr3,tmr4,tmr5,tmr6,tmr7,tmr8,tms1,tms2,tms3,tms4,  &
!            fusl,fdsl,fuir,fdir,fuslr,fdslr,fuirr,fdirr,clds,sd,          &
!            asl_clr,atl_clr,asol,olr,sld,                                 &
!            asol_clr,olr_clr,rld_clr,sld_clr,ss_clr,rs_clr,               &
!            nx,my,my_max,lev)
!-------------------------------------------------------------------
!        call unify_grid                                                   &
!            ( snr,gwr,tg,tm1,tm2,ss,rs,tm3,tm4,ustar,tstar,qstar          &
!            , hflux,qflux,raincu,rainlp,totalp,curate,plcl,cumtop         &
!            , tgclim,gwet,z0,alb,land,ice,ocean,gwclim,acld               &
!            , tmc1,tmc2,tmc3,tmc4,tmc5,tmc6,fpsp,ftp,fqp,fpsp1,ftp1,fqp1  &
!            , e,eps,o3l,dtrad,pt,ptend,ugws,vgws,nx,my,my_max,lev         &
!            , smc,slc,stc,canopy,sigmaf,istyp,ivegtyp,km_soil,temp1       &
!            , temp2,temp3,rld,zice,asl,atl)
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
        deallocate (temp3)
!-------------------------------------------------------------------
        deallocate (tmr1)
        deallocate (tmr2)
        deallocate (tmr3)
        deallocate (tmr4)
        deallocate (tmr5)
        deallocate (tmr6)
        deallocate (tmr7)
        deallocate (tmr8)
        deallocate (tms1)
        deallocate (tms2)
        deallocate (tms3)
        deallocate (tms4)
!-------------------------------------------------------------------
! output gwr or gwet
! in new soil model, gwr did not exist
! but for output required,
! define gwet as near surfce 0.5m soil layer wetness
! define gwr as near surfce 0.5m soil layer moist contain
! gwr=gwet*gwrcc
! gwet(ground wetness) get from smc1*0.2+smc2*0.8
! saturate gwrcc set to be 20mm as original land mode setting
!
!        do jj = 1, jlistnum
!          j=jlist1(jj)
!          nxj=nxdef(j)
!          do i = 1, nxj
!            if( istyp(i,jj).ne.0 )then
!              www = smc(i,1,jj)*0.2+smc(i,2,jj)*0.8
!              gwet(i,jj) = ( www-wlt(istyp(i,jj))) /  &
!                          ( ref(istyp(i,jj))-wlt(istyp(i,jj)) )
!            else
!              gwet(i,jj) = 1.
!            endif
!            gwr(i,jj) = max(0.,min(1.,gwet(i,jj)))*20.
!          enddo
!        enddo
!!        call mpe_unify(gwet,nx,my,2,mpe_double)
!!        call mpe_unify(gwr,nx,my,2,mpe_double)
!
!       rh2100=rh2*100. 
!       rh10100=rh10*100. 
!        call  outflds( 1,nx,my,my_max,lev,ncld                                     &
!                     , lmax,numout,idtg,ifilout,outdir                             &
!                     , ktrop,ptop,capa,cp,rgas,grav,sigma,sgeo                     &
!                     , ptend,pt,plt,pk,pk2,phi,ut,vt,sd                            &
!                     , tt,qt,rdiv,rvor,tg,gwr,z0,hflux,qflux,snr                   &
!                     , raintot,raincu,rainlp,asol,olr,ss,rs,alb,gwclim             &
!                     , acld,cosl,drag,ugws,vgws,t2,rh2100,rh10100,u10,v10,gfx,rld,sld &
!                     , km_soil,smc,slc,stc,canopy,ggdef,slp,v850,v700,h850,h500    &
!                     , ctot,chig,cmid,clow,hpbl,.true.,flash,do_sit)
!        call  outsigs ( 1,nx,my,my_max,lev,ncld                                    &
!                    , idtg,ifilout,ptop,rad,grav                                   &
!                    , cp,cosl,pt,sgeo,snr,gwr,tg,pk,pk2                            &
!                    , ut,vt,tt,qt,phi,rdiv                                         &
!                    , km_soil,smc,slc,stc,canopy,zice,ggdef,gmdef )
!
      endif     ! end of (wrestrt)
!
      do jj = 1, jlistnum
        j=jlist1(jj)
!ch     nxj=nxdef(j)
        call prexp_hybrid_cwb ( nxjp(j),nxp,lev,ptop,sigma,pt(1,jj), &
                          pk(1,1,jj),pk2(1,1,jj),plt(1,1,jj) )
      enddo
!      call  outsigs ( 0,nx,my,my_max,lev,ncld                                    &
!                    , idtg,ifilout,ptop,rad,grav                                 &
!                    , cp,cosl,pt,sgeo,snr,gwr,tg,pk,pk2                          &
!                    , ut,vt,tt,qt,phi,rdiv                                       &
!                    , km_soil,smc,slc,stc,canopy,zice,ggdef,gmdef )

      raintot=0.
      raincu=0.
      rainlp=0.
      raincu6=0.
      rainlp6=0.
      raincu3=0.
      rainlp3=0.
      raincu1=0.
      rainlp1=0.
      gfx=0.
      rld=0.
      sld=0.
      recn=1
!
! read mountant variables for topographic gravity wave drag
!
      if(yesdia .and. dograv .and. nmgwor .eq. 2) then
         call read_mtnvar(nx,my,mtnvar,hprime_b)
!
         if( myrank .eq. 0 ) &
           print*,'read mtnvar=14 hprime_b=',(hprime_b(1,i,1),i=1,mtnvar)
      endif
!

!#ifndef NO_OUT
!       rh2100=rh2*100. 
!       rh10100=rh10*100. 
!      call  outflds( 0,nx,my,my_max,lev,ncld                                       &
!                     , lmax,numout,idtg,ifilout,outdir                             &
!                     , ktrop,ptop,capa,cp,rgas,grav,sigma,sgeo                     &
!                     , ptend,pt,plt,pk,pk2,phi,ut,vt,sd                            &
!                     , tt,qt,rdiv,rvor,tg,gwr,z0,hflux,qflux,snr                   &
!                     , raintot,raincu,rainlp,asol,olr,ss,rs,alb,gwclim             &
!                     , acld,cosl,drag,ugws,vgws,t2,rh2100,rh10100,u10,v10,gfx,rld,sld &
!                     , km_soil,smc,slc,stc,canopy,ggdef,slp,v850,v700,h850,h500    &
!                     , ctot,chig,cmid,clow,hpbl,.true.,flash,do_sit)
!#endif
!      if(typhoon)then
!        do n=1,ntyph
!            i=ixtyp(1,n)
!            j=jytyp(1,n)
!            tensity(0,1,n)=( slp(i,j+1)+slp(i+1,j+1)    &
!                           + slp(i,j  )+slp(i+1,j  ) )/4.
!            tensity(0,2,n)=( v850(i,j+1)+v850(i+1,j+1)  &
!                           + v850(i,j  )+v850(i+1,j  ) )/4.
!            tensity(0,3,n)=( v700(i,j+1)+v700(i+1,j+1)  &
!                           + v700(i,j  )+v700(i+1,j  ) )/4.
!            tensity(0,4,n)=( h850(i,j+1)+h850(i+1,j+1)  &
!                           + h850(i,j  )+h850(i+1,j  ) )/4.
!            tensity(0,5,n)=( h500(i,j+1)+h500(i+1,j+1)  &
!                           + h500(i,j  )+h500(i+1,j  ) )/4.
!        enddo
!      endif
!
! add 40m 100m output for green energy plan
!      if(out_green)then
          do jj = 1, jlistnum
            j=jlist1(jj)
            nxj=nxdef_2d(j)
            do i = 1,nxj
              pklev(i,jj) = pk(i,lev,jj)
            enddo
          enddo
!#ifndef NO_OUT
!        call  outflds_green(0,nx,my,my_max,lev,ncld                &
!              , idtg,ifilout,cp,rgas,grav,t2,u10,v10,ss,pklev      &
!              , sgeo,pt,plt,ptop,ut,vt,tt,qt,cosl,raincu6,rainlp6  &
!              , ggdef,.true.)
!#endif
!      endif
!
      tchange=.false.
      nxmy  = nx*my
      nml   = nx*my*lev
      mlmax2= mlmax*2
      leng  = mlmax*2*lev
!
      itaui=taui+0.1
      itaue=taue+0.1
      itauo=tauo+0.1
      itaup=taup+0.1
      tau=taui
      dtx=dt
!
      forward = itaui .eq. 0
      if (forward)  then
        dta = dtx
      else
        dta = 2*dtx
      endif
!
!  compute initial moisture and potential temperature
!
      qgini   = 0.0
      thdai   = 0.0
      tkei    = 0.0
      tpei    = 0.0
      do i = 1, 4
      do j = 1, my
        wkj(j,i) = 0.
      enddo
      enddo
!
      do jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do k = 1, lev
          do i = 1,nxj
            dsigp = dsigma(k,1)*pt(i,jj)+dsigma(k,2)
            wkj(j,1) = wkj(j,1) + qt(i,k,jj)*dsigp*cosl(j)
            wkj(j,2) = wkj(j,2) + tt(i,k,jj)*dsigp*cosl(j) &
                      /(1.+0.608*qt(i,k,jj))
            wkj(j,3) = wkj(j,3) + (ut(i,k,jj)**2+vt(i,k,jj)**2)*dsigp &
                      *radsq/(2.*cosl(j))
            wkj(j,4) = wkj(j,4) + cp*tt(i,k,jj)*pk(i,k,jj)*dsigp &
                      *cosl(j)
          enddo
        enddo
      enddo
      call mpe_unify(wkj,my,4,1,mpe_double)
!
      cosw=0.
      do j = 1, my
        nxj=nxdef(j)
        cosw = cosw + cosl(j)*nxj
        qgini = qgini + wkj(j,1)
        thdai = thdai + wkj(j,2)
        tkei  = tkei  + wkj(j,3)
        tpei  = tpei  + wkj(j,4)
      enddo
      tengi = tkei+tpei
      qgini = qgini * 100./grav/cosw
      cosw  = cosw*lev
      thdai = thdai/cosw
      tengi = tengi/cosw
      !
      sumtoti=0.
      sumwati=0.
      sumtott=0.
      sumwatt=0.
      do jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        nxjf=nxdef(j)
        do i = 1,nxj
          do k = 1, lev
            dsigp = dsigma(k,1)*pt(i,jj)+dsigma(k,2)
            qtot=0.
            do n = 1, ncld-1
               kk=k+(n-1)*lev
               qtot(i,jj)=qtot(i,jj)+qt(i,kk,jj)
            enddo
            sumtott(i,jj) = sumtott(i,jj)+dsigp
            sumwatt(i,jj) = sumwatt(i,jj)+dsigp*qtot(i,jj)
          enddo
          sumtoti = sumtoti + sumtott(i,jj)*cosl(j)*nx/nxjf
          sumwati = sumwati + sumwatt(i,jj)*cosl(j)*nx/nxjf
        enddo
      enddo
      call mpe_global_sum(sumtoti,1,mpe_double)
      call mpe_global_sum(sumwati,1,mpe_double)
      if( myrank .eq. 0 ) then
        open(35,file='pdry.txt',form='formatted',status='unknown', &
           position='append')
        write(35,*)sumtoti-sumwati,sumwati,sumtoti
        close(35)
        print*,'dry air mass at initial = ',sumtoti-sumwati,' hPa'
        print*,'water  mass at initial = ',sumwati,' hPa'
        print*,'total air mass at initial = ',sumtoti,' hPa'
      endif
!
      if( myrank .eq. 0 ) &
        print*,'qgini, thdai, tengi= ',qgini, thdai, tengi
!
      if(myrank .eq. 0) print*,' beginning integration '
!
!  zero out precip arrays
!
      dt24 = 0.
      do jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i = 1,nxj
          ss24(i,jj)   = 0.
          rs24(i,jj)   = 0.
          hf24(i,jj)   = 0.
          qf24(i,jj)   = 0.
          asol24(i,jj) = 0.
          olr24(i,jj)  = 0.
          rain24(i,jj) = 0.
          rainlp24(i,jj) = 0.
          totallp(i,jj)= 0.
          raincu(i,jj) = 0.
          rainlp(i,jj) = 0.
          raincu6(i,jj)= 0.
          rainlp6(i,jj)= 0.
          raincu3(i,jj)= 0.
          rainlp3(i,jj)= 0.
          raincu1(i,jj)= 0.
          rainlp1(i,jj)= 0.
          raintot(i,jj)= 0.
          runoff(i,jj) = 0.  ! soil
          tmax(i,jj)   = 0.
          tmin(i,jj)   = 0.
          flash24(i,jj)= 0.  !xb110, flash density
        enddo
      enddo
!
      if( itaui .eq. 0 ) then   ! when restart, don't zero out
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do  i = 1,nxj
            ss(i,jj) = 0.
            rs(i,jj) = 0.
            sld(i,jj) = 0.
            tg2 = tg(i,jj)*tg(i,jj)     !soil
            rld(i,jj)= stbo*(tg2*tg2)   !soil
          enddo
        enddo 
      endif
!xb110>
!byl      rmr = 0.
!byl      smr = 0.
!xb110<
!
!***********************************************************************
!     start time integration iterations
!***********************************************************************
!
! sppt
      itimestep=1
!
!!      n_stable=0
!!      n_unstable=0
!!      if(nx.eq.960)then
!!        dt_chg=225.
!!        nc_stable=15
!!        sptendmax2=0.461
!!        sptendmax1=0.4375
!!      else if(nx.eq.1536)then
!!        dt_chg=120.
!!        nc_stable=20
!!        sptendmax2=0.430
!!        sptendmax1=0.4075
!!      else
!!        dt_chg=dt
!!        nc_stable=20
!!        sptendmax2=0.430
!!        sptendmax1=0.4075
!!      endif
!
!      if(typhoon)then
!        dt_trk=6.
!      else
!        dt_trk=tauo
!      endif
! store first idtg to idtg_sst
      icurrenttau=int(tau)
      call dtgfix12(idtg,idtg_sst,icurrenttau)
      if(myrank .eq. 0) print*,'idtg_sst=',idtg_sst

! read opgsst data
      if(lopgsst) then
        call allocate_opgsst_array

        if(myrank .eq. 0) print *,'myrank=',myrank,'idtg_sst=',idtg_sst
          call read_opgsst(idtg_sst,ggdef,ocean,ice)
          icurrentyear=idtg_sst/100000000

      endif   !end lopgsst
!
 10   continue
!
      dtx_tau=dtx/3600.
!
      if(myrank .eq. 0) then 
         print *,'forcast begin tau=',itaui,' to tau=',itaue

! for io quilting
      if(io_quilting)then
         ntag=ntag+1
         call mpe_send_key(keydoit,ntag,istat)
      endif

      endif

!CWB2020 fix the bug for negative timing info
!#ifdef TIMING
!      tm_1=mpi_wtime()
!#endif
 100  continue
!
      tau=tau+dtx/3600.

#ifdef TIMING
      tm_use=tm_2-tm_1
      if(myrank .eq. 0) print 265, tau,tm_use
  265 format(' tau= ',f8.3,',     Timing=',f8.3,' elapse seconds')
      tm_1=mpi_wtime()
#else
      if(myrank .eq. 0) print 265, tau
  265 format(' tau= ',f8.3)
#endif
!
      if(do_sit) then
        if(fsit .le. 0) then
          lrun_sitvdiff=.true.
          ic_sit=-99
          turn_sit=.true.
        else
          if(.NOT. forward) then
            turn_sit= mod(tau+0.001, fsit) .lt. dtx_tau
            if( turn_sit ) ic_sit=0
            if( turn_sit .or. (ic_sit .ge. 0 .AND. ic_sit .lt. nc_sit)) then
              lrun_sitvdiff= .true.
              ic_sit=ic_sit+1
            else
              lrun_sitvdiff= .false.
            endif
          endif
        endif
        if(myrank .eq. 0) then
          print *,'final turn_sit=',turn_sit,',run_sitvdiff=',lrun_sitvdiff
        endif
      endif

!
!!      if(tau.lt.12.)then
!!         hfiltx=hfilt*4.
!!      else if(tau.ge.12. .and. tau.le.24.)then
!!         hfiltx=hfilt*3.
!!      else if(tau.gt.24. .and. tau.le.36.)then
!!         hfiltx=hfilt*2.
!!      else
!!       hfiltx=hfilt
!!      endif
!
!  global mean tempertures (tbar) and specific humid (qbar)
!
!!      sqhaf = sqrt(0.5)
!!      flag  =.false.
!!      do m = 1, mlistnum
!!        mf=mlist(m)
!!        if ( mf.eq.1 ) then
!!          flag=.true.
!!          do k = 1, lev
!!            tbar(k) = sqhaf*temnow(k,1,1,m)
!!          enddo
!!          do k = 1, lev*ncld
!!            qbar(k) = sqhaf*qnow(k,1,1,m)
!!          enddo
!!        endif
!!      enddo
!!      call mpe_broadcast(tbar,lev,flag,mpe_double)
!!      call mpe_broadcast(qbar,lev*ncld,flag,mpe_double)
!
      do m=1,mlistnum
        mf=mlist(m)
        do n=mf,jtrun
          do i = 1, 2
          plten(n,m,i) = 0.0
          pltemp(n,m,i) = 0.0
          do k = 1, levp
            divten(k,i,n,m) = 0.0
            vorten(k,i,n,m) = 0.0
            temten(k,i,n,m) = 0.0
            hldten(k,i,n,m) = 0.0
          enddo
          enddo
        enddo
      enddo
!!      do m=1,mlistnum
!!        mf=mlist(m)
!!        do n=mf,jtrun
!!          do k = 1, levp*ncld*2
!!            qten(k,1,n,m)   = 0.0
!!          enddo
!!        enddo
!!      enddo
 
! Transfer Spectral to Gridpoint for u,v,t,q,ps at n-1
        call transr(jtrun,jtmax,nx,my,my_max,levp,poly,temold,cc,1,nsizey)
        call ujoinsr(cc,ttp,dummy,dummy,dummy,nx,my_max,lev,jlistnum,1,1)
        call tranuv(jtrun,jtmax,nx,my,my_max,levp,onocos,wcfac,wdfac    &
                   ,poly,dpoly,vorold,divold,up,vp,nsizey)
!!        call transr1(jtrun,jtmax,nx,my,my_max,poly,plold,ptp,nsizey)

! for Semi-Lagrangian advection
!
       pdot=0.
       vdmerd=0.
       vdzonl=0.
       vdmerdr=0.
       vdzonlr=0.
       qvadv=0.
       ddtemp=0.
       ddtempr=0.
       pten=0.
!
       ndsldta = 0.5*dta
       ndsldtah= ndsldta/float(itter)

!ch>
! transpose partial to full: ut -> ut_sl, vt -> vt_sl, up -> uum_sl, vp -> vvm_sl, ttp -> ttm_sl, qm -> qm_sl

#ifdef MULTIPLE
      call mpe2d_transpose_ndsl_p2f_multi(ut   ,vt   ,up   ,vp   ,ttp   ,qm   , &
                                    ut_sl,vt_sl,uum_sl,vvm_sl,ttm_sl,qm_sl, &
                                    nxp,nx,levf,levp,ncld,myf,my_max,jlistnum,jlen,nsizex,row_comm,6)
#else
      call mpe2d_transpose_ndsl_p2f(ut,ut_sl,    &
                                    nxp,nx,levf,levp,1,   myf,my_max,jlistnum,jlen,nsizex,row_comm)
      call mpe2d_transpose_ndsl_p2f(vt,vt_sl,    &
                                    nxp,nx,levf,levp,1,   myf,my_max,jlistnum,jlen,nsizex,row_comm)
      call mpe2d_transpose_ndsl_p2f(up,uum_sl,    &
                                    nxp,nx,levf,levp,1,   myf,my_max,jlistnum,jlen,nsizex,row_comm)
      call mpe2d_transpose_ndsl_p2f(vp,vvm_sl,    &
                                    nxp,nx,levf,levp,1,   myf,my_max,jlistnum,jlen,nsizex,row_comm)
      call mpe2d_transpose_ndsl_p2f(ttp,ttm_sl,  &
                                    nxp,nx,levf,levp,1,   myf,my_max,jlistnum,jlen,nsizex,row_comm)
      call mpe2d_transpose_ndsl_p2f(qm,qm_sl,    &
                                    nxp,nx,levf,levp,ncld,myf,my_max,jlistnum,jlen,nsizex,row_comm)
#endif


!!    call mpe2d_unify_nx(pt_sl,pt)
!!    call mpe2d_unify_nx(ptp_sl,ptp)
!ch<

!     Semi-Lagrangian
!       Horizontal Advection
        do itt = 1,itter
        call ndslfv_monoadvh(ttm_sl,qm_sl,pten_sl,uum_sl,vvm_sl  &
                             ,nxdef,ndsldtah,xy,levp)
        enddo


!ch>
! transpose full to partial: ttm_sl -> ddtemp,  pten_sl -> pten, uum_sl -> vdzonl
!                            vvm_sl -> vdmerd,  qm_sl -> qvadv

#ifdef MULTIPLE
      call mpe2d_transpose_ndsl_f2p_multi(ttm_sl,pten_sl,uum_sl,vvm_sl,qm_sl, &
                                    ddtemp   ,pten   ,vdzonl   ,vdmerd   ,qvadv   , &
                                    nxp,nx,levf,levp,ncld,myf,my_max,jlistnum,jlen,nsizex,row_comm)
#else
      call mpe2d_transpose_ndsl_f2p(ttm_sl,ddtemp, &
                                    nxp,nx,levf,levp,1,   myf,my_max,jlistnum,jlen,nsizex,row_comm)
      call mpe2d_transpose_ndsl_f2p(pten_sl,pten,     &
                                    nxp,nx,levf,levp,1,   myf,my_max,jlistnum,jlen,nsizex,row_comm)
      call mpe2d_transpose_ndsl_f2p(uum_sl,vdzonl, &
                                    nxp,nx,levf,levp,1,   myf,my_max,jlistnum,jlen,nsizex,row_comm)
      call mpe2d_transpose_ndsl_f2p(vvm_sl,vdmerd, &
                                    nxp,nx,levf,levp,1,   myf,my_max,jlistnum,jlen,nsizex,row_comm)
      call mpe2d_transpose_ndsl_f2p(qm_sl,qvadv,   &
                                    nxp,nx,levf,levp,ncld,myf,my_max,jlistnum,jlen,nsizex,row_comm)
#endif

!ch<

!
!  the gaussian quadrature loop for spectral tendencies.  subroutine
!  'gridnl' is called for companion mirror image gaussian latitudes
!  these non-linear contributions are then combined in 'rstran' using
!  the symmetry properties of the spherical harmonics
!
      do jj = 1, jlistnum

        j=jlist1(jj)
        nxj=nxdef_2d(j)
!
!   new p**capa quantities were computed in previous diabat call
!
        if (.not.yesdia) call prexp_hybrid_cwb ( nxjp(j),nxp,lev,ptop,sigma,pt(1,jj), &
                               pk(1,1,jj),pk2(1,1,jj),plt(1,1,jj) )
!
!
!       Calculate Vertical velocity & Stream Functions
!
        call gridnl_hybrid_ndsl (nxjp(j),nxp,lev,ncld                  &
        , cp,radsq,ut(1,1,jj),vt(1,1,jj),rdiv(1,1,jj),tt(1,1,jj)       &
        , qt(1,1,jj),phi(1,1,jj),pt(1,jj),dtpl(1,jj),dlpl(1,jj),sinl(j)&
        , pk(1,1,jj),pk2(1,1,jj),dsigma,sigma,onocos(j),cor(j)         &
        , diveng(1,1,jj),vdmerdr(1,1,jj),vdzonlr(1,1,jj),pten(1,1,jj)  &
        , deldm(1,jj),sdpbl(1,jj),sd(1,1,jj),pdot(1,1,jj),vvel(1,1,jj) &
        , sgeo(1,jj) )
!
      enddo !jj = 1,jlistnum

        call joinrs(cc,diveng,dummy,dummy,dummy,nx,my_max,lev,jlistnum,1,1)
        call tranrs(jtrun,jtmax,nx,my,my_max,levp,poly,weight,cc       &
                   ,hldten,1,nsizey)
        call trngra3(jtrun,jtmax,nx,levp,my,my_max,cim,poly,dpoly      &
                   ,hldten,dlphi,dtphi,nsizey)
!
!
!       update all horizontal informations
!
        call ndslfv_update_3tl(nxjp,vdzonl,vdmerd,vdzonlr,vdmerdr,ndsldta)
!
!       Vertical Advection
!
        call ndslfv_monoadvv(ddtemp,qvadv,vdzonl,vdmerd,pdot,pt,nxjp,ndsldta,.false.)


        call mpe2d_unify_nx(ww1,deldm) 
        call tranrs1(jtrun,jtmax,nx,my,my_max,poly,weight,ww1         &
                 ,plten,nsizey)
!
!  before a time step,save previous time u,v,t,q as
!  previous time u,v,t,q which used by physical parameterization
!
      do jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do k = 1, lev
          do i = 1,nxj
            vdzonlr(i,k,jj) = up(i,k,jj)
            vdmerdr(i,k,jj) = vp(i,k,jj)
            ddtempr(i,k,jj) = ttp(i,k,jj)
            up(i,k,jj) = ut(i,k,jj)
            vp(i,k,jj) = vt(i,k,jj)
            ttp(i,k,jj)= tt(i,k,jj)
            ut(i,k,jj) = vdzonl(i,k,jj)
            vt(i,k,jj) = vdmerd(i,k,jj)
            tt(i,k,jj) = ddtemp(i,k,jj)
!!            uum(i,k,jj)= ut(i,k,jj)
!!            vvm(i,k,jj)= vt(i,k,jj)
!!            ttm(i,k,jj)= tt(i,k,jj)
          enddo
        enddo
        do k = 1, lev*ncld
          do i = 1,nxj
            qp(i,k,jj) = qt(i,k,jj)
            qt(i,k,jj) = qvadv(i,k,jj)
!            qmt(i,k,jj) = qt(i,k,jj)
          enddo
        enddo
        do i = 1,nxj
          ptp(i,jj)= pt(i,jj)
        enddo
      enddo

!
      do m = 1, mlistnum
        mf=mlist(m)
        do n = mf, jtrun
          pltemp(n,m,1)= dta*plten(n,m,1)+plold(n,m,1)
          pltemp(n,m,2)= dta*plten(n,m,2)+plold(n,m,2)
        enddo
      enddo
!
!
!  for physical parameterization,output spectrum u,v,t,q to grid point
!
        call transr1(jtrun,jtmax,nx,my,my_max,poly,pltemp,pt,nsizey)
!
            sumtot=0.
      sumwat=0.
      sumtott=0.
      sumwatt=0.
      do jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        nxjf=nxdef(j)
        do i = 1,nxj
          do k = 1, lev
            dsigp = dsigma(k,1)*pt(i,jj)+dsigma(k,2)
            qtot=0.
            do n = 1, ncld-1
               kk=k+(n-1)*lev
               qtot(i,jj)=qtot(i,jj)+qt(i,kk,jj)
            enddo
            sumtott(i,jj) = sumtott(i,jj)+dsigp
            sumwatt(i,jj) = sumwatt(i,jj)+dsigp*qtot(i,jj)
          enddo
          sumtot = sumtot + sumtott(i,jj)*cosl(j)*nx/nxjf
          sumwat = sumwat + sumwatt(i,jj)*cosl(j)*nx/nxjf
        enddo
      enddo
      call mpe_global_sum(sumtot,1,mpe_double)
      call mpe_global_sum(sumwat,1,mpe_double)
      if( myrank .eq. 0 ) then
        open(35,file='pdry.txt',form='formatted',status='unknown', &
           position='append')
        write(35,*)sumtot-sumwat,sumwat,sumtot
        close(35)
        print*,'dry air mass = ',sumtot-sumwat,' hPa'
        print*,'water mass = ',sumwat,' hPa'
        print*,'total air mass = ',sumtot,' hPa'
      endif
!
        !  stochastic_physics
        call run_stochastic_physics()       

        if (yesdia)  then
!
          call diabat ( forward,docup,dodry,dolsp,dopbl,dorad,doshl,dograv,tofd &
                      , nx,my,my_max,lev,ncld,nmcup,nmpbl,nmland,nmshl,cgw      &
                      , idg,jdg,ldiag,dtx,tau,hours,julian,year,yrd             &
                      , frad,ozon,njump,itypbl,ktcup,ktpbl,ktshl,grav           &
                      , rgas,cp,stbo,s0,evaprh,hltm,ptop,sigma,dsigma,il,ib     &
                      , cof,xlat,xlon,sgeo,z0,alb,land,ocean,ice                &
                      , snr,tg,tgclim,curate,plcl,cumtop,totalp,raintot,raincu  &
                      , rainlp,raincu6,rainlp6,raincu3,rainlp3,raincu1,rainlp1  &
                      , hflux,qflux,ustar,tstar,qstar,e                         &
                      , eps,o3l,dtrad,ss,rs,plt,pk,pk2                          &
                      , ptp,    up,    vp,   ttp,qp                             &
                      , pt ,    ut,    vt,    tt,qt                             &
                      , gwclim,tice,hice,qgini,thdai,tengi                      &
                      , acld,std,asol,olr,drag,ugws,vgws                        &
                      , sdpbl,t2,q2,rh2,rh10,u10,v10,gfx                        &
                      , fm,fh,fm10,fh2,srflag                                   &
                      , rld,km_soil,smc,stc,canopy,runoff                       &
                      , sigmaf,istyp,ivegtyp,wltsmc,refsmc,maxsmc,dfkt,xktk,dfk &
                      , ftp,fqp,fpsp,ftp1,fqp1,fpsp1,deltaq,cnvwr,cnvcr,pdot    &
                      , shdmax,shdmin,snoalb                                    &
                      , slopetyp,sld,slc,zice,cice,xtice,sncover,sndepth        &
                      , ctot,chig,cmid,clow,hpbl,asl,atl,cosz                   &
                      , nmgwor,nmgwcv,hprime_b,mtnvar,docgrav,nmmiph            &
!--------------------------------------------------------------------------------
                      , fusl,fdsl,fuir,fdir                                     &
                      , fuslr,fdslr,fuirr,fdirr                                 &
                      , asl_clr,atl_clr,clds                                    &
                      , ss_clr,rs_clr,asol_clr,olr_clr,sld_clr,rld_clr          &
                      , alvsf,alvwf,alnsf,alnwf,facsf,facwf                     &
                      , idtg,doo3l,nfxr,sfalb,sfemis,isot,ivegsrc               &
                      , itimestep,lrun_sitvdiff,ic_sit                          &
!xb110>
!byl                      , rmr,smr,flash)
                      , flash,tsflw,vvel,totallp)      
!--------------------------------------------------------------------------------
!
! add reynolds stress
!
          call rayleifr(nx,my,my_max,lev,rad,cosl,dt,ut,vt)
        endif    ! end of (yesdia)
        itimestep=itimestep+1   ! for sppt time evolution)
!
!  after phyical parameterization,transform grid point u,v,t,q to
!  spectrum
!
          do jj = 1, jlistnum
            j=jlist1(jj)
            nxj=nxdef_2d(j)
            do k = 1, lev
              do i = 1,nxj
                vdzonl(i,k,jj) = ( ut(i,k,jj) - vdzonlr(i,k,jj) ) / dta
                vdmerd(i,k,jj) = ( vdmerdr(i,k,jj) - vt(i,k,jj) ) / dta
                ddtemp(i,k,jj) = ( tt(i,k,jj) - ddtempr(i,k,jj) ) / dta
              enddo
            enddo
          enddo
!
          call joinrs(cc,ddtemp,dummy,dummy,dummy,nx,my_max,lev        &
                   ,jlistnum,1,1)
          call tranrs(jtrun,jtmax,nx,my,my_max,levp,poly,weight,cc     &
                   ,temten,1,nsizey)
          call rstrandz (jtrun,jtmax,nx,my,my_max,levp,vdmerd,vdzonl   &
                   ,weight,cim,onocos,poly,dpoly,divten,vorten,nsizey)
!
      if (lsimpl)  then
!
!   apply semi-implicit adjustemts to above explicitly computed
!   tendencies to stablize integration for long time steps
!
      call siimpl ( jtrun,jtmax,lev,dta,ptmeans,dsigma,spalm,eps4,eigval &
                  , evecin,evectr,arrhyd,arsddt,temold,divold,plold      &
                  , temnow,divnow,plnow,temten,divten,plten,alpha)
!
      endif
!
!  zero out global mean tendencies for divergence, vorticity, and
!  terrain pressure to ensure consistency with gauss's theorem.
!
      mlst=ilist(1)
      if(mlst .ne. 0) then
        plten(1,mlst,1) = 0.0
        plten(1,mlst,2) = 0.0
      endif
      do m = 1, mlistnum
         mf=mlist(m)
         if ( mf.eq.1 ) then
         do i = 1, 2
           do k = 1, levp
             divten(k,i,1,m)= 0.0
             vorten(k,i,1,m)= 0.0
           enddo
         enddo
         endif
      enddo
!
      call transr1( jtrun,jtmax,nx,my,my_max,poly,plten,ptend,nsizey)
      do mf = 1, jtrun
        wkmf(mf) = 0.
      enddo
      sptend= 0.0
      do m = 1,mlistnum
        mf=mlist(m)
        do n = mf, jtrun
          if(n.ne.1)then
            wkmf(mf)= wkmf(mf) + plten(n,m,1)**2 + plten(n,m,2)**2
          endif
        enddo
      enddo
      call mpe_unify(wkmf,1,jtrun,3,mpe_double) 
      do mf = 1, jtrun
        sptend= sptend + wkmf(mf)
      enddo
!
      sptend= sqrt(0.5*sptend)*3600.0
      if(myrank .eq. 0)  &
         print *, '  surf pres tend rms =',sptend,' mb/hrs'
!
!  take a time step
!
      if (forward)  then
!
        do m = 1, mlistnum
          mf=mlist(m)
          do n = mf, jtrun
            do i = 1, 2
            do k = 1, levp
              vornow(k,i,n,m)= dta*vorten(k,i,n,m) + vorold(k,i,n,m)
              divnow(k,i,n,m)= dta*divten(k,i,n,m) + divold(k,i,n,m)
              temnow(k,i,n,m)= dta*temten(k,i,n,m) + temold(k,i,n,m)
            enddo
            enddo
          enddo
        enddo
!
        do m = 1, mlistnum
          mf=mlist(m)
          do n = mf, jtrun
            plnow(n,m,1)= dta*plten(n,m,1)+plold(n,m,1)
            plnow(n,m,2)= dta*plten(n,m,2)+plold(n,m,2)
          enddo
        enddo

!
      if (hdiff) call hdiffu_3tl ( dta,my,my_max,nx,jtrun,jtmax,lev,ncld   &
                             , hfiltx,rad,cosl,up,vp,vornow,divnow,temnow  &
                             , eps4,trefs)
!
        forward=.false.
        dta= 2.0*dtx
        lsitstart=.false.
!
      else ! if ( forward )
!
!  leap-frog/time filtered steps
!
        do m =1,mlistnum
          mf=mlist(m)
          do n  = mf,jtrun
            do i  = 1, 2
            do k  = 1, levp
              vorten(k,i,n,m)= dta*vorten(k,i,n,m) + vorold(k,i,n,m)
              divten(k,i,n,m)= dta*divten(k,i,n,m) + divold(k,i,n,m)
              temten(k,i,n,m)= dta*temten(k,i,n,m) + temold(k,i,n,m)
            enddo
            enddo
          enddo
        enddo
!
        do m =1,mlistnum
          mf=mlist(m)
          do n  = mf,jtrun
            plten(n,m,1)= dta*plten(n,m,1)+plold(n,m,1)
            plten(n,m,2)= dta*plten(n,m,2)+plold(n,m,2)
          enddo
        enddo
!
        if (hdiff) call hdiffu_3tl ( dta,my,my_max,nx,jtrun,jtmax,lev,ncld   &
                               , hfiltx,rad,cosl,up,vp,vorten,divten,temten  &
                               , eps4,trefs)
!
! accumulate some flux every time step to output point (24 hour)
! 1994 11 11
!
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do i = 1,nxj
            hf24(i,jj)  = hf24(i,jj)+hflux(i,jj)*dtx
            qf24(i,jj)  = qf24(i,jj)+qflux(i,jj)*dtx
            ss24(i,jj)  = ss24(i,jj)+ss(i,jj)*dtx
            rs24(i,jj)  = rs24(i,jj)+rs(i,jj)*dtx
            asol24(i,jj)= asol24(i,jj)+asol(i,jj)*dtx
            olr24(i,jj) = olr24(i,jj)+olr(i,jj)*dtx
            rain24(i,jj)= rain24(i,jj)+totalp(i,jj)
            rainlp24(i,jj)= rainlp24(i,jj)+totallp(i,jj)
            flash24(i,jj)=flash24(i,jj)+flash(i,jj)*dtx
          enddo
        enddo
        dt24 = dt24 + dtx
!
!  robert time filter
!
        do m = 1, mlistnum
          mf=mlist(m)
          do n = mf, jtrun
            do i = 1, 2
            do k = 1, levp
              vorold(k,i,n,m)= vornow(k,i,n,m) + tfilt*( vorold(k,i,n,m) &
                             - 2.0*vornow(k,i,n,m)+vorten(k,i,n,m) )
              vornow(k,i,n,m)= vorten(k,i,n,m)
              divold(k,i,n,m)= divnow(k,i,n,m) + tfilt*( divold(k,i,n,m) &
                             - 2.0*divnow(k,i,n,m)+divten(k,i,n,m) )
              divnow(k,i,n,m)= divten(k,i,n,m)
              temold(k,i,n,m)= temnow(k,i,n,m) + tfilt*( temold(k,i,n,m) &
                             - 2.0*temnow(k,i,n,m)+temten(k,i,n,m) )
              temnow(k,i,n,m)= temten(k,i,n,m)
            enddo
            enddo
          enddo
        enddo
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do k = 1, lev*ncld
            do i = 1, nxj
              qm(i,k,jj) = qp(i,k,jj) + tfilt*( qm(i,k,jj)           &
                         - 2.0*qp(i,k,jj) + qt(i,k,jj) )
            enddo
          enddo
        enddo
        do m =1,mlistnum
          mf=mlist(m)
          do n  = mf,jtrun
            plold(n,m,1)  = plnow(n,m,1) + tfilt*( plold(n,m,1) &
                          - 2.0*plnow(n,m,1)+plten(n,m,1) )
            plnow(n,m,1)  = plten(n,m,1)
            plold(n,m,2)  = plnow(n,m,2) + tfilt*( plold(n,m,2) &
                          - 2.0*plnow(n,m,2)+plten(n,m,2) )
            plnow(n,m,2)  = plten(n,m,2)
          enddo
        enddo
!
        if ( mod(itimestep,2) .eq. 0 ) xy = -1 * xy
      endif     ! end of (forward)
!
!  from updated spectral variables compute corresponding grid
!  point fields
!
!!      call joinsr(wss3,vornow,divnow,temnow,dummy,jtrun,jtmax,levp &
!!                 ,mlistnum,3,1)
!!      call transr(jtrun,jtmax,nx,my,my_max,levp,poly,wss3,cc3,3,nsizey)
!!      call ujoinsr(cc3,rvor,rdiv,tt,dummy,nx,my_max,lev,jlistnum,3,1)

!!      call transr(jtrun,jtmax,nx,my,my_max,levp,poly,vornow,cc,1,nsizey)
!!      call ujoinsr(cc,rvor,dummy,dummy,dummy,nx,my_max,lev,jlistnum,1,1)
      call transr(jtrun,jtmax,nx,my,my_max,levp,poly,divnow,cc,1,nsizey)
      call ujoinsr(cc,rdiv,dummy,dummy,dummy,nx,my_max,lev,jlistnum,1,1)
      call transr(jtrun,jtmax,nx,my,my_max,levp,poly,temnow,cc,1,nsizey)
      call ujoinsr(cc,tt,dummy,dummy,dummy,nx,my_max,lev,jlistnum,1,1)
      call transr1(jtrun,jtmax,nx,my,my_max,poly,plnow,pt,nsizey)
!
!  zonal and meridional gradients of terrain pressure
!
      call trngra (jtrun,jtmax,nx,my,my_max,cim,poly,dpoly,plnow &
                   ,dlpl,dtpl,nsizey)
!
!  velocity components
!
      call tranuv (jtrun,jtmax,nx,my,my_max,levp,onocos,wcfac,wdfac &
                  , poly,dpoly,vornow,divnow,ut,vt,nsizey)
!
!   computing new p**capa quantities
!
      do jj = 1, jlistnum
        j=jlist1(jj)
!ch     nxj=nxdef(j)
        call prexp_hybrid_cwb ( nxjp(j),nxp,lev,ptop,sigma,pt(1,jj), &
                          pk(1,1,jj),pk2(1,1,jj),plt(1,1,jj) )
      enddo
!c    make sure moisture field is positive
!
!      do jj =1, jlistnum
!        j=jlist1(jj)
!        nxj=nxdef(j)
!        do n=1,ncld
!          qbrrow(n,j)=0.
!        enddo
!        call postq (nxj,nx,lev,ncld,dsigma,pt(1,jj),qt(1,1,jj),qbrrow(1,j),cosl(j))
!      enddo
!
!  detact instability occure or not
!
!!      if(tau.gt.24.) then
!!!        sptendmax2=0.439
!!!        sptendmax1=0.409
!!        if(sptend.le.sptendmax1)n_stable=n_stable+1
!!        if(sptend.gt.sptendmax2)n_unstable=n_unstable+1
!!        if(myrank .eq. 0) print *,'n_stable=',n_stable,' n_unstable=',n_unstable
!!        if( mod(tau+0.001, 1.) .lt. dtx_tau)then
!!          if(n_stable .gt. nc_stable)then
!!            dtx=dt_chg
!!            dta=dtx*2
!!            dtx_tau=dtx/3600.
!!            if(myrank .eq. 0)print *,'** stable change dt=',dtx,' sec**'
!!          else if(n_unstable .gt. nc_stable)then
!!            dtx=dt
!!            dta=dtx*2
!!            dtx_tau=dtx/3600.
!!            hfiltx=hfilt*.4
!!            if(myrank .eq. 0)print *,'** unstable change dt=',dtx,' sec**'
!!          else
!!            if(myrank .eq. 0)print *,'** keep dt=',dtx,' sec**'
!!          endif
!!          n_unstable=0
!!          n_stable=0
!!        endif
!!      endif
!
! output sit var. at "outsitmean" interval
!
      if(do_sit .and. (outsitmean .gt. 0.)) then
        dtaup = mod(tau+0.001, outsitmean)
        if( dtaup .lt. dtx_tau ) then
          ntau=tau+0.001
          if(myrank .eq. 0) print *,'outsitmean at tau=',tau
          call writesitmean(nx,my,lkvl,ntau,idtg,ggdef)
        endif
      endif
!
! output flux at 24 hour interval
!
      itau  = NINT( tau )
      dtaup = mod(tau+0.001, 24.)

      if( dtaup .lt. dtx_tau ) then
        ntau=tau+0.001
        if(myrank .eq. 0) print *,'out24 at tau=',tau
!        call mpe_unify(hf24,nx,my,2,mpe_double)
!        call mpe_unify(qf24,nx,my,2,mpe_double)
!        call mpe_unify(ss24,nx,my,2,mpe_double)
!        call mpe_unify(rs24,nx,my,2,mpe_double)
!        call mpe_unify(asol24,nx,my,2,mpe_double)
!        call mpe_unify(olr24,nx,my,2,mpe_double)
!        call mpe_unify(rain24,nx,my,2,mpe_double)
#ifndef NO_OUT
        call out24(nx,my,my_max,hf24,qf24,ss24,rs24,asol24,olr24,rain24,rainlp24,dt24 &
                  ,glob,itau,idtg,ggdef,flash24)
#endif
        if(do_sit)then
          if(loutsit24)then
            call outsit24(nx,my,lkvl,ntau,idtg,ggdef)
          endif
          call dtgfix12(idtg,idtg_temp,ntau-1)
          ibeforeyymm=idtg_temp/1000000
          call dtgfix12(idtg,idtg_temp,ntau)
          icurrentyymm=idtg_temp/1000000
          lnewyymm=icurrentyymm/=ibeforeyymm
          if(myrank .eq. 0) then
            print *,'ibeforeyymm=',ibeforeyymm,'icurrentyymm=',icurrentyymm &
                   ,'lnewyymm=',lnewyymm
          endif
          if( lnewyymm ) then
            call outsitmon(nx,my,lkvl,ntau,idtg,ggdef)
          endif
       endif
!
!  zero set arrays
!
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do  i = 1,nxj
            hf24(i,jj)   = 0.
            qf24(i,jj)   = 0.
            ss24(i,jj)   = 0.
            rs24(i,jj)   = 0.
            asol24(i,jj) = 0.
            olr24(i,jj)  = 0.
            rain24(i,jj) = 0.
            rainlp24(i,jj) = 0.
            flash24(i,jj)= 0.
          enddo
        enddo
        dt24 = 0.
      endif 
!
      dtaup= mod(tau+0.001, tauo)
      histim=(dtaup .lt. dtx_tau)
!       if(myrank.eq.0)print *,'chkhis dtaup,tauo,dtx_tau=',dtaup,tauo,dtx_tau

!  for tracker
      dt_trk=real(trk_intv)
      dtaup= mod(tau+0.001, dt_trk)
      ltrack=(dtaup .lt. dtx_tau)
!       if(myrank.eq.0)print *,'chkltr dtaup,dt_trk,dtx_tau=',dtaup,dt_trk,dtx_tau
!
!--- histim (start)
!
! output gwr or gwet
! in new soil model, gwr did not exist
! but for output required,
! define gwet as near surfce 0.5m soil layer wetness
! define gwr as near surfce 0.5m soil layer moist contain
! gwr=gwet*gwrcc
! gwet(ground wetness) get from smc1*0.2+smc2*0.8
! saturate gwrcc set to be 20mm as original land mode setting
!
      do jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i = 1,nxj
          if(istyp(i,jj).ne.0)then
            www = smc(i,1,jj)*0.2+smc(i,2,jj)*0.8
            gwet(i,jj) = ( www-wltsmc(istyp(i,jj)) ) / &
                        ( refsmc(istyp(i,jj))-wltsmc(istyp(i,jj)) )
          else
            gwet(i,jj) = 1.0
          endif
          gwr(i,jj) = max(0.,min(1.,gwet(i,jj)))*20.0
        enddo
      enddo
!      call mpe_unify(gwet,nx,my,2,mpe_double)
!      call mpe_unify(gwr,nx,my,2,mpe_double)
!
      if ( histim .or. ltrack )  then

        itau = tau + 0.1
        if(myrank .eq. 0) print *,' history file written at tau= ', itau
!
        if(wrestrt)then
          call chlen (cwbout,48,lcwb)
          call chlen (phyout,48,lphy)
          write (ctau,800) itau
  800     format('tau',i3.3)
          rfile = cwbout(1:lcwb)//ctau

!!          allocate (work_io((( (7+2*ncld)*2*lev+8)*jtrun*jtmax+my*lev)*nsize))
 
!!          call gather_spec(work_io,vornow,divnow,temnow,qnow,plnow,    &
!!             vorold,divold,temold,qold,plold,dsqgeo,spgeo,trefs,       &
!!             uzm,lev,ncld,jtrun,jtmax,my,nsize)
!
!  add if check to let restart output performed every 24h
!
!          if( mod(float(itau),24.) .lt. 0.01 ) then
!            if(myrank .eq. 0) then
!              open (unit=7,file=rfile,form='unformatted')
!cc            write (7) work_io
!cc            call flush(7)
!              close (7)
!            endif
!          endif
!
!!          deallocate (work_io)
          rfile = phyout(1:lphy)//ctau
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
!
!          call unify_grid                                                  &
!              ( snr,gwr,tg,tm1,tm2,ss,rs,tm3,tm4,ustar,tstar,qstar         &
!              , hflux,qflux,raincu,rainlp,totalp,curate,plcl,cumtop        &
!              , tgclim,gwet,z0,alb,land,ice,ocean,gwclim,acld              &
!              , tmc1,tmc2,tmc3,tmc4,tmc5,tmc6,fpsp,ftp,fqp,fpsp1,ftp1,fqp1 &
!              , e,eps,o3l,dtrad,pt,ptend,ugws,vgws,nx,my,my_max,lev        &
!              , smc,slc,stc,canopy,sigmaf,istyp,ivegtyp,km_soil,temp1      &
!              , temp2,temp3,rld,zice,asl,atl)
!
!-------------------------------------------------------------------
         allocate (tmr1(nx,lev+1,my))
         allocate (tmr2(nx,lev+1,my))
         allocate (tmr3(nx,lev+1,my))
         allocate (tmr4(nx,lev+1,my))
         allocate (tmr5(nx,lev+1,my))
         allocate (tmr6(nx,lev+1,my))
         allocate (tmr7(nx,lev+1,my))
         allocate (tmr8(nx,lev+1,my))
         allocate (tms1(nx,lev,my))
         allocate (tms2(nx,lev,my))
         allocate (tms3(nx,lev,my))
         allocate (tms4(nx,lev,my))
!-------------------------------------------------------------------
!         call unify_gridr                                                 &
!          ( tmr1,tmr2,tmr3,tmr4,tmr5,tmr6,tmr7,tmr8,tms1,tms2,tms3,tms4   &
!          , fusl,fdsl,fuir,fdir,fuslr,fdslr,fuirr,fdirr,clds,sd           &
!          , asl_clr,atl_clr,asol,olr,sld                                  &
!          , asol_clr,olr_clr,rld_clr,sld_clr,ss_clr,rs_clr                &
!          , nx,my,my_max,lev)
!-------------------------------------------------------------------
!
!          if( mod(float(itau),24.) .lt. 0.01 ) then
!            if(myrank .eq. 0) then
!              open (unit=10,file=rfile,form='unformatted')
!cc           write(10) snr,gwr,tg,tm1,tm2,ss,rs,tm3,tm4,ustar,tstar,qstar
!cc  1                , hflux,qflux,raincu,rainlp,totalp,curate,plcl,cumtop
!cc  2                , tgclim,gwet,z0,alb,land,ice,ocean,gwclim,acld
!cc  3                , temp1,temp2,canopy,sigmaf,istyp,ivegtyp,rld
!hmhj3                , tmc1,tmc2,fpsp
!cc           call flush(10)
!              close (10)
!            endif
!          endif
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
          deallocate (temp3)
!---------------------------------------------------------------------
          deallocate (tmr1)
          deallocate (tmr2)
          deallocate (tmr3)
          deallocate (tmr4)
          deallocate (tmr5)
          deallocate (tmr6)
          deallocate (tmr7)
          deallocate (tmr8)
          deallocate (tms1)
          deallocate (tms2)
          deallocate (tms3)
          deallocate (tms4)
!---------------------------------------------------------------------
!
          if(do_sit) then
            if(myrank .eq. 0) print *, 'ready rerun_sitgrid1'
              call rerun_sitgrid1(itau)
            if(myrank .eq. 0) print *, 'ready rerun_sitgrid2'
              call rerun_sitgrid2(itau)
            if(myrank .eq. 0) print *, 'ready rerun_sitgrid3'
              call rerun_sitgrid3(itau)
          endif
        endif   ! end of (wrestrt)
!
!        dtaup = mod( tauo+0.001, taup )
!ds        dtaup = mod( tau+0.001, taup )
!ds       if(tau.le.(taureg+0.001) .and. dtaup.lt.0.01)then
       if(tau.le.(taureg+0.001) .and. histim )then
!jh        if( histim ) then 
#ifndef NO_OUT
        call outsigs ( itau,nx,my,my_max,lev,ncld        &
                     , idtg,ptop,rad,grav                &
                     , cp,cosl,pt,sgeo,snr,gwr,tg,pk,pk2 &
                     , ut,vt,tt,qt,phi,rdiv,km_soil,smc  &
                     , slc,stc,canopy,zice,ggdef,gmdef )
#endif
        endif
!-------------------------------------------------------------------------------

!      if (myrank .eq. 0) print *,'ioutsigr =',ioutsigr
      if (ioutsigr .eq. 0 .and. histim ) then
      if (myrank .eq. 0) print *,'outsigr start !!!'
#ifndef NO_OUT
      call outsigr ( itau,nx,my,my_max,lev                                     &
                   , idtg                                                      &
                   , fusl,fdsl,fuir,fdir                                       &
                   , fuslr,fdslr,fuirr,fdirr                                   &
                   , asl,atl,asl_clr,atl_clr                                   &
                   , dtrad,clds,sd                                             &
                   , ss,rs,olr,asol,sld,rld                                    &
                   , ss_clr,rs_clr,olr_clr                                     &
                   , asol_clr,sld_clr,rld_clr                                  &
                   , cice,xtice,snr,sncover,snoalb                             &
                   , ctot,chig,cmid,clow                                       &
                   , ggdef,gmdef)
#endif
      if (myrank .eq. 0) print *,'outsigr ok !!!'
      endif
!-------------------------------------------------------------------------------
!  write operational fields (p-levels and surface)
!
#ifndef NO_OUT
!
       call transr(jtrun,jtmax,nx,my,my_max,levp,poly,vornow,cc,1,nsizey)
       call ujoinsr(cc,rvor,dummy,dummy,dummy,nx,my_max,lev,jlistnum,1,1)
       call  outflds( itau,nx,my,my_max,lev,ncld                               &
                    , lmax,numout,idtg,outdir                                  &
                    , ktrop,ptop,capa,cp,rgas,grav,sigma,sgeo                  &
                    , ptend,pt,plt,pk,pk2,phi,ut,vt,vvel                       &
                    , tt,qt,rdiv,rvor,tg,gwr,z0,hflux,qflux,snr                &
                    , raintot,raincu,rainlp,asol,olr,ss,rs,alb,gwclim          &
                    , acld,cosl,drag,ugws,vgws,t2,q2,rh2,rh10,u10,v10,gfx,rld,sld &
!byl                    , km_soil,smc,slc,stc,canopy,ggdef,slp,v850,v700,h850,h500 &
                    , km_soil,smc,slc,stc,canopy,ggdef,typtrk                  &
!xb110                    , ctot,chig,cmid,clow,hpbl,histim,flash,do_sit)
                    , ctot,chig,cmid,clow,hpbl,histim,do_sit)
#endif
!
!        if(typhoon .and. ltrack)then
         if(typhoon .and. ltrack .and. itau .le. 384 )then
          if(myrank .eq. 0)print *,' calling tracking,  tau= ',tau
!          if(myrank .eq. 0)print *,'dt_trk=',dt_trk
          call tracking(tau,dt_trk,dt,nx,my,                            &
                  ntyph,typname,ixtyp,jytyp,tlon,tlat,tflon,tflat,idtg, &
                  nrec,typhoon,tensity)
        endif
!
!  zero out precip arrays
!  
!  add 12-hour check to let prep. amount be of 12-hour accumulation for
!  every 12-hour output, but prep. amount still 24-hour accumulation for
!  every 24-hour output without 12-hour output points
!  (in order to be consistent to the rule followed by nfs, 11/30/2000)
!
        if( mod(float(itau)+0.00001, 12.) .lt. 0.01 ) then
!        if( histim .and. (mod(float(itau)+0.00001, 12.) .lt. 0.01) ) then
          raincu=0.
          rainlp=0.
          runoff=0.
!          call zilch (raincu,nxmy)
!          call zilch (rainlp,nxmy)
!          call zilch (runoff,nxmy)
        endif
!
        ifromtau = taui
        itotau   = taue
        flag =.false.
        if(myrank .eq. 0)then
          flag =.true.
!CWB2016 
          if(.not. io_quilting)then
!CWB2017           call sendmsg ('gfs',ifromtau,itotau,istat)
            if(itau.eq.itotau) call sendmsg ('gfs',ifromtau,itotau,istat)
          else
            istat=0
          endif
        endif
!ch     call mpe_broadcast(istat,1,flag,mpe_integer)
        call mpe_bcast(istat,1,0,mpe_integer)

        if (istat.eq.-1) then
          print *,' SENDMSG ERROR ', 'RANK=',myrank
          call dmsexit(-1)
        endif
!
        if(myrank .eq. 0) print *,' history written at tau=',itau

      endif     ! end of (histim)

! out green energy plan
!      if(out_green)then
!        if (mod(tau+0.00001, otgreen) .lt. 0.01) then
      if( out_green .and. mod( itau , nint(otgreen) ) == 0 )then
#ifndef NO_OUT
        call  outflds_green(nint(tau),nx,my,my_max,lev,ncld                &
                          , idtg,cp,rgas,grav,t2,u10,v10,ss,pk             &
                          , sgeo,pt,plt,ptop,ut,vt,tt,qt,cosl,raincu6,rainlp6  &
                          , ggdef)
#endif
!          if (mod(tau+0.00001, 6. ) .lt. 0.01) then
          if ( mod( itau , 6 ) == 0 ) then
           raincu6=0.
           rainlp6=0.
          endif
!        endif
      endif
! output f006 data for FV3
      if (abs(tau+0.00001-6.) .lt. 0.01) then
#ifndef NO_OUT
        if (myrank .eq. 0) print *,'output FV3 data !!!'
        call outflds_fv3(nint(tau),nx,my,my_max,idtg,ggdef                     &
                          ,q2,fm,fh,fm10,fh2,srflag,ustar)
#endif
      endif
! output 3hr accumulated precipitation for CWB home page
      if(out_hp)then
        if (mod(tau+0.00001, 3.) .lt. 0.01) then
#ifndef NO_OUT
        call  outflds_hp(nint(tau),nx,my,my_max,idtg                           &
                        ,raincu3,rainlp3,ggdef)
#endif
          raincu3=0.
          rainlp3=0.
        endif
      endif

! calculate Tmax Tmin @ 2m from T2
          do jj = 1, jlistnum
            j=jlist1(jj)
            nxj=nxdef_2d(j)
            do i=1,nxj
             if (tmax(i,jj).eq.0.0)tmax(i,jj) = t2(i,jj)
             if (tmin(i,jj).eq.0.0)tmin(i,jj) = t2(i,jj)
             if (t2(i,jj).gt. tmax(i,jj)) then
              tmax(i,jj)= t2(i,jj)
             elseif (t2(i,jj).lt. tmin(i,jj)) then
              tmin(i,jj)= t2(i,jj)
             endif
            enddo
          enddo
!
!kc
!  output t2,raintot,u10,v10,ctot at 1 hour interval within 192hr. 
!kc
!      if(domfc)then
!==xb118
! change domfc type from logical to real
!==xb118

        dtaup = mod(tau+0.001, 1.)
!        if(myrank .eq. 0) print*,'domfc at tau,dtaup=',tau,dtaup
!        if(tau.le.192. .and. dtaup.lt.0.01)then
        if(tau.le.(domfc+0.001) .and. dtaup.lt.0.01)then
          ntau=tau+0.001
          if(myrank .eq. 0) print *,'out1 at tau=',tau
          do jj = 1, jlistnum
            j=jlist1(jj)
            nxj=nxdef_2d(j)
            do i=1,nxj
!!              slpty(i,jj)= pt(i,jj)+pdiff(i,jj)
              rain1(i,jj)= raincu1(i,jj)+rainlp1(i,jj)
              temp= tda*(t2(i,jj)-273.15)/(tdb+(t2(i,jj)-273.15))+     &
                     log(rh2(i,jj))
              td(i,jj)= (tdb*temp/(tda-temp))+273.15
            enddo
          enddo
#ifndef NO_OUT
          call out2d_mfc(nx,lev,my,my_max,itau,idtg                 &
                      ,raincu1,rainlp1,raintot,glob,t2,q2,rh2,rh10  &
                      ,u10,v10,tmax,tmin,rld,sld,ctot,pt,ggdef)
#endif
          raincu1=0.0
          rainlp1=0.0
          tmax=0.0
          tmin=0.0
        endif
!      endif
!
        if (dosppt .and.  dospptout .and. mod( itau , 1 ) == 0 ) then
           call spptout(tau)
        endif
        if (doskeb .and.  doskebout .and. mod( itau , 1 ) == 0 ) then
           call skebout(tau)
        endif
!
!for update sst(W00100), seaice(W00091), snowdepth(B00650)
!
      if(ldailyFCTsst .OR. ldailyFCTicesndpt .OR. (dailyClm_option.ge.1)) then
!        CALL read_dailyFCT(idtg,tau,dtx,tg,cice,sndepth)
!        icurrenttau=int(tau)
!        tauleft=float(int((tau-int(tau)+0.001)*3600./dtx))*dtx
        tautemp=tau+dtx/3600.
        CALL read_dailyFCT(idtg,tautemp,dtx,tg,cice,sndepth,xlon,xlat,ocean)
        icurrenttau=int(tautemp)
        tauleft=float(int((tautemp-int(tautemp)+0.001)*3600./dtx))*dtx

        if(myrank .eq. 0) then
          print *,'icurrenttau=',icurrenttau, &
                  ',tauleft=',tauleft
        endif
        if(tauleft > 3599.) then
          icurrenttau=icurrenttau+1
          tauleft=0.
        endif
        call dtgfix12(idtg_sst,idtg_temp,icurrenttau)
        call time_weights(idtg_temp,tauleft)
        write(cdtg,'(i12)')idtg_temp
        read(cdtg,'(i4,i2,i2,i2,i2)')yr,mo,dy,hr,mn
!        ssttau = mod(tau+0.001, 24.)
        tautemp=float(hr)+tauleft/3600.
        ssttau = mod(tautemp+0.001, 24.)

        if( lopgsst .AND. (ssttau .lt. dtx_tau) ) then
!  read opgsst data
          if(myrank .eq. 0) then
            print *,'ready opgsst_weights,icurrenttau=',icurrenttau, &
                    ',tauleft=',tauleft
          endif
          do jj=1,jlistnum
            j=jlist1(jj)
            nxj=nxdef_2d(j)
            do ii = 1, nxj
              sst(ii,jj)=wgt1*opgsst(ii,jj,nmw1)+wgt2*opgsst(ii,jj,nmw2)
              if((myrank .eq.2).AND.(jj.EQ.4).AND.(i.EQ.670)) then
                print *,'myrank=',myrank,'tauleft=',tauleft &
                       ,'icurrenttau=',icurrenttau &
                       ,'tautemp=',tautemp,'ssttau=',ssttau &
                       ,'idtg1_sst=',idtg1_sst,'wgt1=',wgt1,'wgt2=',wgt2 &
                       ,'nmw1=',nmw1,'nmw2=',nmw2,'tg=',tg(ii,jj) &
                       ,'opgsst(ii,jj,nmw1)=',opgsst(ii,jj,nmw1) &
                       ,'opgsst(ii,jj,nmw2)=',opgsst(ii,jj,nmw2)
              endif
            enddo
          enddo
!          call mpe_unify(sst,nx,my,2,mpe_double)

        endif

        if(lFCTweight .OR. ssttau .lt. dtx_tau) then
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do ii=1,nxj
            if(dailyClm_option .eq. 1 )then       !persistent anomaly sst
             ssttemp=ANAsstT0(ii,jj)-dailyClmANAsst(ii,jj,0) &
                     +obswtbwgt1*dailyClmANAsst(ii,jj,obswtbnmw1) &
                     +obswtbwgt2*dailyClmANAsst(ii,jj,obswtbnmw2)
            elseif(dailyClm_option .eq. 2 )then   !idea from Yuejian Zhu(2018 JGR)
              wweight=min(tautemp/24./35.,1.)
              ssttemp=(1.-wweight)*(ANAsstT0(ii,jj)-dailyClmANAsst(ii,jj,0) &
                     +dailyClmANAsst(ii,jj,1) ) +wweight*(dailyFCTsst(ii,jj,1) &
                     -(dailyClmFCTsst(ii,jj,1)-dailyClmANAsst(ii,jj,1)))
            else
              ssttemp=obswtbwgt1*dailyFCTsst(ii,jj,obswtbnmw1) &
                      +obswtbwgt2*dailyFCTsst(ii,jj,obswtbnmw2)
            endif
            sst(ii,jj)=merge(ssttemp,tg(ii,jj),(ssttemp.GE.271. .AND. ssttemp.LT.400.))
            if(ldailyFCTicesndpt)then
                snrtemp=obswtbwgt1*dailyFCTsndepth(ii,jj,1)+obswtbwgt2*dailyFCTsndepth(ii,jj,2)
                cicetemp=obswtbwgt1*dailyFCTcice(ii,jj,1)+obswtbwgt2*dailyFCTcice(ii,jj,2)
                snr(ii,jj)=max(0.,snrtemp)
                sndepth(ii,jj)=max(0.,snrtemp)*8.   !same as getrdy.f90
                sncover(ii,jj)=min(1.,max(0.,snrtemp)/400.)   !same as getrdy.f90
                cice(ii,jj)=max(0.,cicetemp)
                oceanold(ii,jj)=ocean(ii,jj)
                iceold(ii,jj)=ice(ii,jj)
                if( ocean(ii,jj) .or. ice(ii,jj) )then
!                  alb(i,j)=0.09
                  ocean(ii,jj)=.true.
                  ice(ii,jj)= .false.
                  if( cice(ii,jj) .ge. 0.5) then
                    ice(ii,jj)=.true.
                    ocean(ii,jj)=.false.
                    if(.not. iceold(ii,jj))then
                      alb(ii,jj)=0.55
                      tgclim(ii,jj)=271.2
                      xtice(ii,jj)=max(271.2,sst(ii,j))
                      cice(ii,jj)=max(0.5,cice(ii,jj))
                      if( myrank .eq. 0 ) print *, "ocean to ice(",ii,",",jj,"),sst=",sst(ii,jj)
                    else
                      if(.not. oceanold(i,jj))then
                        alb(ii,jj)=0.09
                      endif
                    endif
                  endif
                endif
            endif
            if(ocean(ii,jj)) then
              if(do_sit .and. (sitmask(ii,jj) .eq. 1.)) then
                obswtb(ii,jj)=max(271.,sst(ii,jj))
!                  tg(i,jj)=max(271.,sst(i,j))
!                  tgold(i,jj)=max(271.,sst(i,j))
!                  tsw(i,jj)=max(271.,sst(i,j))
              else
                  tg(ii,jj)=max(271.,sst(ii,jj))
              endif
            endif
            if(myrank .eq. myrank_check .AND. jj .eq. jj_check .AND. ii .eq. ii_check) then
              print *,"ssttau=",ssttau,",dtx_tau=",dtx_tau   &
               ,",obswtbwgt1=",obswtbwgt1,",obswtbwgt2=",obswtbwgt2 &
               ,",tautemp=",tautemp,",ssttau=",ssttau &
               ,",tg=",tg(ii,jj)
              if(ldailyFCTsst) then
                print *,",dailyFCTsst1=",dailyFCTsst(ii,jj,1)   &
                       ,",dailyFCTsst2=",dailyFCTsst(ii,jj,2)
              endif
              if(dailyClm_option .ge. 1) then
                print *,",ANAsstT0=",ANAsstT0(ii,jj)   &
                 ,",dailyClmANAsst0=",dailyClmANAsst(ii,jj,0) &
                 ,",dailyClmANAsst1=",dailyClmANAsst(ii,jj,1)
                if (dailyClm_option .eq. 2) then
                  print *,",dailyClmFCTsst0=",dailyClmFCTsst(ii,jj,0) &
                    ,",dailyClmFCTsst1=",dailyClmFCTsst(ii,jj,1)
                endif
              endif
              if(do_sit) print *,"obswtb=",obswtb(ii,jj)
            endif
          end do
        end do
        endif !end if(lFCTweight .OR. ssttau .lt. dtx_tau) then

      endif   !end ldailyFCT

!      inexttau = int(tau+dtx/3600.+0.001)
!      call dtgfix12(idtg,idtg_temp,inexttau)
!      inextyymmdd=idtg_temp/10000
!      lnewday=icurrentyymmdd/=inextyymmdd
!      if(lnewday)then
        if(do_sit .AND. lgodas .AND. ldailysst) then
          CALL read_dailygodas(idtg,tau,dtx)
        endif
!      endif
!
!
! new year, read obs sst
!
      if(lopgsst) then
        inexttau = int(tau+dtx/3600.+0.001)
        call dtgfix12(idtg,idtg1_sst,inexttau)
        inextyear=idtg1_sst/100000000
        lnewyear=icurrentyear/=inextyear
        if( lnewyear ) then
          call read_opgsst(idtg1_sst,ggdef,ocean,ice)
          icurrentyear=idtg1_sst/100000000
        endif
      endif   !end lopgsst
!
      itau=tau+0.001
#ifdef TIMING
      tm_2=mpi_wtime()
#endif
      if(itau.lt.itaue) go to 100
!
      flag =.false.
      if(myrank .eq. 0)then
        call recmsg ('gfs',ifromtau,itotau,istat)
        flag =.true.
      endif
!ch   call mpe_broadcast(istat,1,flag,mpe_integer)
      call mpe_bcast(istat,1,0,mpe_integer)

      if (istat.eq.-1)  then
        print *,' RECMSG ERROR ', 'RANK=',myrank
        call dmsexit(-1)
      else if (istat.eq.1)  then
        if(myrank .eq. 0) then
           print *,' finished integration '

! for io quilting
      if(io_quilting)then
         ntag=ntag+1
         call mpe_send_key(keydoit,ntag,istat)
         ntag=ntag+1
         call mpe_send_key(keydone,ntag,istat)
      endif

           call dmscls(bckfile,istat)
           call dmscls(ifilin,istst)
           call dmscls(ifilout,istst)
           if(ldailyFCTsst)      call dmscls(ifilin_sst,istat)
           if(ldailyFCTicesndpt) call dmscls(ifilin_ncep,istat)
           if(dailyClm_option .ge. 1) then
             call dmscls(ifilin_ClmANA,istat)
             if(dailyClm_option .eq. 2) then
               call dmscls(ifilin_ClmFCT,istat)
             endif
           endif
           close(77)
        end if
        if(ldailyFCTsst .OR. ldailyFCTicesndpt .OR. (dailyClm_option.ge.1) ) then
          call deallocate_dailyFCT_array
        endif
        if(lopgsst) then
          call deallocate_opgsst_array
        endif
        if(do_sit) then
           print *,'myrank=',myrank,'ready sit_vdiff_end'
           if(locaf) call deallocate_ocaf_array
           if(lwoa0) call deallocate_woa0_array
           if(lgodas) call deallocate_godas_array
           call deallocate_sitgrid_array
           call cleanup_netcdf
           call sit_vdiff_end
        endif
        return
      endif
!
!ch   call mpe_broadcast(ifromtau,1,flag,mpe_integer)
!ch   call mpe_broadcast(itotau,1,flag,mpe_integer)
      call mpe_bcast(ifromtau,1,0,mpe_integer)
      call mpe_bcast(itotau,1,0,mpe_integer)
!
      taui = float(ifromtau)
      taue = float(itotau)
      tauo = float(itotau)
      itaui=taui+0.1
      itaue=taue+0.1
      itauo=tauo+0.1
!
      go to 10
!
      return
      end
