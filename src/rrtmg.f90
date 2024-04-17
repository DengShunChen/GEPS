      subroutine rrtmg                                              &
! -------------------------------------------------------------------
!    -  inputs: 
           ( sigma,pst,plt,std,tt,qt,o3l,sd,tg,                     &
             slimsk,cice,xtice,snr,sncover,snoalb,z0,               &
             alvsg,alnsg,alvwg,alnwg,facsg,facwg,                   &
             curate,icsdswg,icsdlwg,                                &
             sinlj,coslj,xlatj,xlonr,jdat,d2r,xkapa,                &
             ptop,dtlw,dtsw,lsswr,lslwr,lssav,                      &
             nfxr,j,                                                &
             nx,nxj,lev,ncld,lprnt,ipt,kdt,                         &
             uni_cloud,lmfshal,lmfdeep2,                            &
             deltaq,sup,cnvw,cnvc,                                  &
             ftp,ftp1,fqp,fqp1,nmmiph,                              &
!    -  outputs:
             asol,olr,ss,rs,sld,rld,tsflwr,                         &
             ctot,chig,cmid,clow,                                   &
             cldcov,htrsw,htrlw,                                    &
             fusl,fdsl,fuir,fdir,                                   &
! -------------------------------------------------------------------
             fuslr,fdslr,fuirr,fdirr,                               &
             htrsw0,htrlw0,cosz,                                    &
             asol_clr,olr_clr,ss_clr,rs_clr,                        &
             sld_clr,rld_clr,sfalb_g,semis_g )
! -------------------------------------------------------------------
! --- for RRTMG scheme :
!
      use physpara
      use module_radiation_driver, only : grrad
! -------------------------------------------------------------------
      use mpe
      use rank
      use index
      use radn
      use const, only: RTYPE
!      use noah, only:ioutsigr
! -------------------------------------------------------------------
! --- for rrtmg input :
!
      integer i,k,kc,n
      integer ntrac,nfxr,nx,nxj,lev,ipt,ncld,kdt,nmmiph,nclds
! --- 3d parameters
!
      real    plt(nx,lev),std(nx),tg(nx),                         &
              tt(nx,lev),sd(nx,lev)
      real(kind=RTYPE)    qt(nx,lev*ncld),o3l(nx,lev),pst(nx),    &
                          sigma(lev+1,2)
!
! --- 2d parameters
!
      real    slimsk(nx),cice(nx),xtice(nx),snr(nx),sncover(nx),  &
              snoalb(nx),z0(nx)
      real    alvsg(nx),alvwg(nx),alnsg(nx),alnwg(nx),facsg(nx),  &
              facwg(nx),curate(nx),xlonr(nx),tsflwr(nx),cosz(nx)
      integer icsdlwg(nx),icsdswg(nx),jdat(8),j
      real    sinlj,coslj,xlatj,ptop,dtlw,dtsw,d2r,xkapa      
      logical lsswr,lslwr,lssav,lprnt
      logical uni_cloud,lmfshal,lmfdeep2
      real    www,cmax,cmin,imax,imin,tem1,tem2


! --- for grrad input/output (local) :
!
! --- 3d
!
      real    prsi(nx,lev+1)                      !for even levels
      real    prslk(nx,lev),prsl(nx,lev)          !for odd levels
      real    qgrs(nx,lev),tgrs(nx,lev)           !for odd levels
!      real    tracer(nx,lev,ncld),vvl(nx,lev)
      real    vvl(nx,lev)
      real, dimension(:,:,:), allocatable ::  tracer

      real    slmsk(nxj),xlon(nxj),xlat(nxj),tsfc(nxj),            &
              snowd(nxj),sncovr(nxj),zorl(nxj),hprim(nxj),         &
              cv(nxj),cvt(nxj),cvb(nxj),snalb(nxj)

      real    alvsf(nxj),alvwf(nxj),alnsf(nxj),alnwf(nxj),         &
              facsf(nxj),facwf(nxj),                               &
              fice(nxj),tisfc(nxj),sinlat(nxj),coslat(nxj)

      real    sfalb(nxj),coszen(nxj),coszdg(nxj)
      real    tsflw(nxj),semis(nxj) 

      integer icsdlw(nxj),icsdsw(nxj)
! --- for pdf cloud
      real    sup
      real    deltaq(nx,lev),cnvw(nx,lev),cnvc(nx,lev)
! --- for WSM6 & Thompson & GFDL MP
      real    ftp(nx,lev),ftp1(nx,lev),fqp(nx,lev),fqp1(nx,lev)
      real    phy3d(nxj,lev,5)

! -------------------------------------------------------------------
! --- for rrtmg output:
      real    asol(nx),olr(nx),ss(nx),rs(nx),sld(nx),rld(nx)
      real    sfalb_g(nx),semis_g(nx)
!
! --- 3d
!
      real    htrsw(nx,lev),htrlw(nx,lev)
      real(kind=RTYPE)    fusl(nx,lev+1),fdsl(nx,lev+1)
      real(kind=RTYPE)    fuir(nx,lev+1),fdir(nx,lev+1)
     
      real    dummy1(nx,lev),dummy2(nx,lev)
      real    work1(nx,lev+1),work2(nx,lev+1)
      real    work3(nx,lev+1),work4(nx,lev+1)
!
! --- 2d
!
!      real    dtrad(nx,lev)
      real    ctot(nx),chig(nx),cmid(nx),clow(nx),csbl(nx)

! for WSM6 & Thompson & GFDL MP
      real    cldcov(nx,lev)   ! input/output layer cloud fraction
      real    dummy3(nxj,lev)
      

!
! --- for clear sky
!
      real    asol_clr(nx),olr_clr(nx),ss_clr(nx),rs_clr(nx)
      real    sld_clr(nx),rld_clr(nx)
!
! --- 3d 
!
      real    htrsw0(nx,lev),htrlw0(nx,lev)
      real(kind=RTYPE)    fuslr(nx,lev+1),fdslr(nx,lev+1)
      real(kind=RTYPE)    fuirr(nx,lev+1),fdirr(nx,lev+1)
!
      real    dummy4(nx,lev),dummy5(nx,lev)
      real    work5(nx,lev+1),work6(nx,lev+1)
      real    work7(nx,lev+1),work8(nx,lev+1)
!
! --- 2d
!
      real    fluxr(nx,nfxr)
!
!     if (myrank .eq. 0) print *,'### in rrtmg.f ###'
!     if (myrank .eq. 0) print *,'### j=',j
! -------------------------------------------------------------------
!    to set variables  for grrad input
! -------------------------------------------------------------------

      if(ntoz.eq.0)then
        ntrac=ncld+1
      else
        ntrac=ncld
      endif
      allocate(tracer(nx,lev,ntrac))

      do i=1,nxj
        prsi(i,lev+1)=(sigma(1,1)*pst(i)+sigma(1,2)+ptop)*0.1
      enddo

      do k=1,lev
         kc=lev-k+1
      do i=1,nxj
         prsi(i,kc)=(sigma(k+1,1)*pst(i)+sigma(k+1,2)+ptop)*0.1
         prsl(i,kc)=plt(i,k)*0.1 ! change to cb
         prslk(i,kc)=(plt(i,k)/1000.)**xkapa
         tgrs(i,kc)=tt(i,k)
         qgrs(i,kc)=qt(i,k)
         vvl(i,kc)=sd(i,k)*0.1                !cb/sec
         dummy3(i,kc)=cldcov(i,k)
      enddo
      enddo

      do n = 1, ntrac
      do k = 1, lev
         kc=lev-k+1
      do i = 1, nxj
         tracer(i,kc,n) = qt(i,k+(n-1)*lev)
      enddo
      enddo
      enddo
!     if (myrank .eq. 0) print *,'### j=',j
!     if (myrank .eq. 0) print *,'tracer(1,60,1)=',tracer(1,60,1) 
!     if (myrank .eq. 0) print *,'tracer(1,60,2)=',tracer(1,60,2) 

!      do k = 1, lev
!         kc=lev-k+1
!      do i = 1, nxj
!         tracer(i,kc,ntoz) = o3l(i,k)
!      enddo
!      enddo
!
!      no need the reduction of O3 concentration over model top.
!      02/23/2024 proposed by Jen-Her Chen
!      
!      do k=1,8
!         fac_o3=k*0.1
!         if(fac_o3 .le. 0.3) fac_o3=0.3
!         kc=lev-k+1
!      do i = 1, nxj
!!       tracer(i,kc,ntoz) = o3l(i,k)*fac_o3
!       tracer(i,kc,ntoz) = tracer(i,kc,ntoz)*fac_o3
!      end do
!      end do
      if ( nmmiph.eq.6 .or. nmmiph.eq.8 .or. nmmiph.eq.18 ) then
        nclds=3
! for MP WSM6 & Thompson effective radius
        do k = 1, lev
          do i = 1, nxj
            phy3d(i,k,1) = ftp(i,k)
            phy3d(i,k,2) = ftp1(i,k)
            phy3d(i,k,3) = fqp(i,k)
          enddo
        enddo
      endif

      if ( nmmiph.eq.11 .or. nmmiph.eq.12 .or. nmmiph.eq.13 ) then
! for MP GFDL effective radius
        nclds = 5  ! number of effective cloud condensates used in radiation processes
        do k = 1, lev
          do i = 1, nxj
            phy3d(i,k,1) = ftp(i,k)    ! effective radius for liquid water (micron)
            phy3d(i,k,2) = ftp1(i,k)   ! effective radius for ice water    (micron)
            phy3d(i,k,3) = fqp(i,k)    ! effective radius for snow water   (micron)
            phy3d(i,k,4) = fqp1(i,k)   ! effective radius for rain water   (micron)
          enddo
        enddo
      endif

      if ( nmmiph.eq.15 .or. nmmiph.eq.16) then
! for MP Goddard (GCE) effective radius
        nclds = 6
        do k = 1, lev
          do i = 1, nxj
            phy3d(i,k,1) = ftp(i,k)    ! effective radius for liquid water (micron)
            phy3d(i,k,2) = ftp1(i,k)   ! effective radius for ice water    (micron)
            phy3d(i,k,3) = fqp(i,k)    ! effective radius for snow water   (micron)
            phy3d(i,k,4) = fqp1(i,k)   ! effective radius for rain water   (micron)
          enddo
        enddo
      endif

      if ( nmmiph.eq.2 ) then
        nclds=1
        phy3d=0.
      endif
!
!     if (myrank .eq. 0) print *,'### j=',j
!     if (myrank .eq. 0) print *,'tracer(1,60,3)=',tracer(1,60,3) 
 
!--------------------------------------------------------------------
!    to set xlat(nx), sinlat(nx), coslat(nx) for grrad : input
!    to set qgrs(nx,lev), tgrs(nx,lev) for grrad : input
!--------------------------------------------------------------------

      do i = 1, nxj
         slmsk(i)=slimsk(i)
         xlon(i)=xlonr(i)
         xlat(i)=xlatj*d2r
!        tsfc(i)=tt(i,lev)          !surface temp in k
         tsfc(i)=tg(i)
         snowd(i)=snr(i)
         sncovr(i)=sncover(i)
         snalb(i)=snoalb(i)
         zorl(i)=z0(i)*100.         !surface roughness in cm
         hprim(i)=std(i)
         alvsf(i)=alvsg(i)
         alnsf(i)=alnsg(i)
         alvwf(i)=alvwg(i)
         alnwf(i)=alnwg(i)
         facsf(i)=facsg(i)
         facwf(i)=facwg(i)
         fice(i)=cice(i)
         tisfc(i)=xtice(i)
         sinlat(i)=sinlj
         coslat(i)=coslj
         cvt(i)=0.   
         cvb(i)=0. 
         icsdsw(i)=icsdswg(i)
         icsdlw(i)=icsdlwg(i)
      enddo

!     if (myrank .eq. 0) print *,'### j=',j
!     if (myrank .eq. 0) print *,'xlatj=',xlatj
!     if (myrank .eq. 0) print *,'xlatj=',xlat(1)

!     if (myrank .eq. 0) print *,'### j=',j
!     if (myrank .eq. 0) print *,'sinlj=',sinlj
!     if (myrank .eq. 0) print *,'coslj=',coslj
!     if (myrank .eq. 0) print *,'sinlat(1)=',sinlat(1)
!     if (myrank .eq. 0) print *,'coslat(1)=',coslat(1)

!     if (myrank .eq. 0) print *,'### j=',j
!     if (myrank .eq. 0) print *,'qgrs(1,1)=',qgrs(1,1)
!     if (myrank .eq. 0) print *,'qgrs(1,lev)=',qgrs(1,lev)
!     if (myrank .eq. 0) print *,'qgrs(ipt,1)=',qgrs(ipt,1)
!     if (myrank .eq. 0) print *,'qgrs(ipt,lev)=',qgrs(ipt,lev)
!     if (myrank .eq. 0) print *,'tgrs(1,1)=',tgrs(1,1)
!     if (myrank .eq. 0) print *,'tgrs(1,lev)=',tgrs(1,lev)
!     if (myrank .eq. 0) print *,'tgrs(ipt,1)=',tgrs(ipt,1)
!     if (myrank .eq. 0) print *,'tgrs(ipt,lev)=',tgrs(ipt,lev)
!------------------------------------------------------------------------

      cmax=-1.0e+25
      cmin=1.0e+25
      imin=1
      imax=1
      do i = 1, nxj
      www = curate(i) * 0.0416667 * 0.1
      www = max ( www, 5.531e-04 )
      cv(i) = 0.93 + 0.124*log(www)
      if(cv(i) .gt. cmax) then
         cmax = cv(i)
         imax = i
      endif
      if(cv(i) .lt. cmin) then
         cmin = cv(i)
         imin = i
      endif
      enddo

!     if (myrank .eq. 0) then
!     if (myrank .eq. 0) print *,'### j=',j
!         print *,'### iter=',iter
!         print *,'### cloud fraction : cv-max =',cmax,' imax=',imax
!         print *,'### cloud fraction : cv-min =',cmin,' imin=',imin
!     endif

      tem1= 0.1
      tem2= 0.8
      do i = 1, nxj
      if ( cv(i) .lt. tem1 )  cv(i) = 0.0
      if ( cv(i) .gt. tem2 )  cv(i) = tem2 + 0.01
      enddo
!
!--------------------------------------------------------------------
!      if (myrank .eq. 0) then
!     print *,'ptop=', ptop
!     pstmax=-1.0e10
!     pstmin=1.0e10
!     imax=1
!     imin=1
!     do i=1,nxj
!        if (pst(i) .gt. pstmax) then
!            pstmax=pst(i)
!            imax=i
!        endif
!        if (pst(i) .lt. pstmin ) then
!            pstmin=pst(i)
!            imin=i
!        endif
!     enddo
!     if (myrank .eq. 0) print *,'### j=',j
!     print *,'pst-max=',pstmax,' imax=',imax
!     print *,'pst-min=',pstmin,' imin=',imin
!     endif
!
!     if (myrank .eq. 0) then
!     pmax=-1.0e10
!     pmin=1.0e10
!     imax=1
!     imin=1
!     kmax=1
!     kmin=1
!     do k=1,lev+1
!     do i=1,nxj
!        if (prsi(i,k) .gt. pmax) then
!            pmax=prsi(i,k)
!            imax=i
!            kmax=k
!        endif
!        if (prsi(i,k) .lt. pmin) then
!            pmin=prsi(i,k)
!            imin=i
!            kmin=k
!        endif
!     enddo
!     enddo
!     if (myrank .eq. 0) print *,'### j=',j
!     print *,'### prsi : pressure at even level(ptop/pbot)###'
!     print *,'prsi-max=',pmax,' imax=',imax,' kmax=',kmax
!     print *,'prsi-min=',pmin,' imin=',imin,' kmin=',kmin
!     print *,'prsi(ipt,lev+1)=',prsi(ipt,lev+1)
!     endif
!--------------------------------------------------------------------
!     if (myrank .eq. 0) then
!     pmax=-1.0e10
!     pmin=1.0e10
!     imax=1
!     imin=1
!     kmax=1
!     kmin=1
!     do k=1,lev
!     do i=1,nxj
!        if (prsl(i,k) .gt. pmax) then
!            pmax=prsl(i,k)
!            imax=i
!            kmax=k
!        endif
!        if (prsl(i,k) .lt. pmin) then
!            pmin=prsl(i,k)
!            imin=i
!            kmin=k
!        endif
!     enddo
!     enddo
!     if (myrank .eq. 0) print *,'### j=',j
!     print *,'### prsl : pressure at odd level(mid-level)##'
!     print *,'prsl-max=',pmax,' imax=',imax,' kmax=',kmax
!     print *,'prsl-min=',pmin,' imin=',imin,' kmin=',kmin
!     print *,'prsl(ipt,lev)=',prsl(ipt,lev)
!     endif
!--------------------------------------------------------------------
!     if (myrank .eq. 0) then
!     pmax=-1.0e10
!     pmin=1.0e10
!     imax=1
!     imin=1
!     kmax=1
!     kmin=1
!     do k=1,lev
!     do i=1,nxj
!        if (prslk(i,k) .gt. pmax) then
!            pmax=prslk(i,k)
!            imax=i
!            kmax=k
!        endif
!        if (prslk(i,k) .lt. pmin) then
!            pmin=prslk(i,k)
!            imin=i
!            kmin=k
!        endif
!     enddo
!     enddo
!     if (myrank .eq. 0) print *,'### j=',j
!     print *,'### prslk : (p_odd/1000.)**xkapa ###'
!     print *,'prslk-max=',pmax,' imax=',imax,' kmax=',kmax
!     print *,'prslk-min=',pmin,' imin=',imin,' kmin=',kmin
!     print *,'prslk(ipt,lev)=',prslk(ipt,lev)
!     endif
!--------------------------------------------------------------------
!     if (myrank .eq. 0) then
!     pmax=-1.0e10
!     pmin=1.0e10
!     imax=1
!     imin=1
!     kmax=1
!     kmin=1
!     do k=1,lev
!     do i=1,nxj
!        if (plt(i,k) .gt. pmax) then
!            pmax=plt(i,k)
!            imax=i
!            kmax=k
!        endif
!        if (plt(i,k) .lt. pmin) then
!            pmin=plt(i,k)
!            imin=i
!            kmin=k
!        endif
!     enddo
!     enddo
!     if (myrank .eq. 0) print *,'### j=',j
!     print *,'plt-max=',pmax,' imax=',imax,' kmax=',kmax
!     print *,'plt-min=',pmin,' imin=',imin,' kmin=',kmin
!     print *,'plt(ipt,lev)=',plt(ipt,lev)
!     endif
!--------------------------------------------------------------------
          call grrad(prsi,prsl,prslk,tgrs,qgrs,tracer,vvl,slmsk,     &
             xlon,xlat,tsfc,snowd,sncovr,snalb,zorl,hprim,           &
             alvsf,alnsf,alvwf,alnwf,facsf,facwf,fice,tisfc,         &
             sinlat,coslat,solhr,jdat,solcon,                        &
             cv,cvt,cvb,                                             &
             icsdsw,icsdlw,ntcw,nclds,ntoz,ntrac,nfxr,               &
             dtlw,dtsw,lsswr,lslwr,lssav,                            &
             nx,nxj,lev,me,lprnt,ipt,kdt,myrank,                     &
             ntiw,ntrw,ntsw,ntgl,uni_cloud,lmfshal,lmfdeep2,         &
             deltaq,sup,cnvw,cnvc,phy3d,                             &
!  ---  outputs:
             dummy1,sfalb,coszen,coszdg,                             &
             dummy2,tsflw,semis,                                     &
!  ---  input/output:
             dummy3,fluxr,                                           &
!  ---  optional outputs:
             dummy4,dummy5,                                          &
             work1,work2,work3,work4,                                &
             work5,work6,work7,work8                                 &
             )

!       if (myrank .eq. 0) print *,'grrad ok!'

       do i=1,nxj
!
! total sky
!
          ss(i)   = fluxr(i,4)-fluxr(i,5)
          rs(i)   = fluxr(i,7)-fluxr(i,6)
          asol(i) = fluxr(i,1)-fluxr(i,2)
          olr(i)  = fluxr(i,3)
          sld(i)  = fluxr(i,4)
          rld(i)  = fluxr(i,6)
! clear sky
          ss_clr(i)   = fluxr(i,24)-fluxr(i,25)
          rs_clr(i)   = fluxr(i,27)-fluxr(i,26)
          asol_clr(i) = fluxr(i,1 )-fluxr(i,22)
          olr_clr(i)  = fluxr(i,23)
          sld_clr(i)  = fluxr(i,24)
          rld_clr(i)  = fluxr(i,26)
!
! cloud fraction
!
          chig(i)=fluxr(i,8)   ! high cloud fraction
          cmid(i)=fluxr(i,9)   ! middle cloud fraction
          clow(i)=fluxr(i,10)  ! low cloud fraction
          ctot(i)=fluxr(i,20)  ! total cloud fraction
!
          sfalb_g(i)=sfalb(i)
          semis_g(i)=semis(i)  ! surface emissivity
       enddo

       do k = 1, lev
          kc=lev-k+1
       do i = 1, nxj
          htrsw(i,kc)=dummy1(i,k)*86400.
          htrlw(i,kc)=dummy2(i,k)*86400.
          cldcov(i,kc)=dummy3(i,k)
          htrsw0(i,kc)=dummy4(i,k)*86400.
          htrlw0(i,kc)=dummy5(i,k)*86400.
       enddo
       enddo

!       do k = 1, lev
       do i = 1, nxj
!          dtrad(i,k)=htrsw(i,k)+htrlw(i,k)
          tsflwr(i)=tsflw(i)
          cosz(i)  =coszen(i)
       enddo
!       enddo

!       if(ioutsigr == 0)then
       do k = 1, lev+1
          kc=lev-k+2
       do i = 1, nxj
          fusl(i,kc)=work1(i,k)
          fdsl(i,kc)=work2(i,k)
          fuir(i,kc)=work3(i,k)
          fdir(i,kc)=work4(i,k)
          fuslr(i,kc)=work5(i,k)
          fdslr(i,kc)=work6(i,k)
          fuirr(i,kc)=work7(i,k)
          fdirr(i,kc)=work8(i,k)
       enddo
       enddo
!       endif ! ioutsigr .eq. 0
       deallocate(tracer)

       return
       end
