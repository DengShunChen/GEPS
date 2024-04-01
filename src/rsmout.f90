    subroutine rsmout(idtg,itau,nx,my,my_max,lev,ncld            &
             , ptop,cp,rgas,grav,sgeo,pdiff                        &
             , t1000,pt,plt,pk,pk2,phi,ut,vt                       &
             , tt,qt,tg,snr,cosl                                   &
             , km,smc,stc                                          &
             , ice,land,ocean)

#ifdef RSM
!
!  output driver subroutine to process sigma level data to standard
!  pressure surfaces and standard grids
!
      use mpe
      use rank
      use index
      use const, only: rlon1, rlon2, rlat1, rlat2, rgrdsz, RTYPE
      use noah,  only: cice
      implicit  none

      integer(kind=8)  :: idtg
      integer   itau,nx,my,my_max,lev,ncld,km
      real      ptop,cp,rgas,grav

      real      pdiff(nxp,my_max)                                     &
              , t1000(nxp,my_max),plt(nxp,lev,my_max)                 &
              , tg(nxp,my_max),snr(nxp,my_max)                        &
!soil
              , smc(nxp,km,my_max),stc(nxp,km,my_max)
      real(kind=RTYPE) ut(nxp,lev,my_max),vt(nxp,lev,my_max)          &
              ,        tt(nxp,lev,my_max),qt(nxp,lev*ncld,my_max)     &
              ,        phi(nxp,lev,my_max),sgeo(nxp,my_max)           &
              ,        pt(nxp,my_max)                                 &
              ,        pk(nxp,lev,my_max),pk2(nxp,lev,my_max)         &
              ,        cosl(my)

      logical   land(nxp,my_max),ocean(nxp,my_max),ice(nxp,my_max)

!
! local work arrays
!
      real      slp(nxp,my_max),tmp(nxp,lev,my_max),plog(nxp,lev,my_max) &
              , pllp(nxp,my_max),slmsk(nxp,my_max)
      real(kind=RTYPE) glob(nx,my)
      real      tens(lev+1)             
!
      integer,  parameter :: lpout = 47 
      real      wrk1(nxp,lev),pout(lpout),pkout(lpout),phistd(lpout) &
              , bt1(nxp,my_max),bt2(nxp,my_max)
      real(kind=RTYPE) pres3d(nxp,my_max,lpout),hld1(nxp,my_max)     &
              ,        hld2(nxp,my_max)
!
      real(kind=RTYPE) soil_xy(nxp,my_max,8)   ! the last dim is changable
!
      integer   nxmy,jj,j,nxj,k,i,n,nk,ngq,ntt,kk,ntrac,ii
      integer   llts
      real      rad,ograv,alaps,rdg,ttb,ttp,ttt,ttt1,ttt2,anlslp
      real      apha,pl1000,splog,ax,bx,cx,dx,tmid,tsf,tadia,xxx
!
!yj rsm dimension 
      ! global: lon 0~359.75, lat -90~90
      integer  nx2,my2,ierr
      ! default: lon 60~179.75, lat -20~70
      integer nxs, mys
      integer x1, x2, y1, y2
      real,allocatable ::  rsmoutp(:,:),                          &
     &     temp_gfs(:,:,:),spfh_gfs(:,:,:),clwr_gfs(:,:,:),           &
     &     rain_gfs(:,:,:),qice_gfs(:,:,:),snow_gfs(:,:,:),grpl_gfs(:,:,:), &
     &     ozon_gfs(:,:,:),geop_gfs(:,:),u_gfs(:,:,:),v_gfs(:,:,:),     &   
     &     tg_gfs(:,:),smc_gfs(:,:,:),snr_gfs(:,:),stc_gfs(:,:,:),      &
     &     ice_gfs(:,:),                                               &
!yj2019
     &     terr_gfs(:,:),slmsk_gfs(:,:)

      real,allocatable :: xr(:),yr(:)
!
      character cdum1*14
      character cdum2*12
      data cdum1/'-------------'/
      data cdum2/'--rsmout---'/
!
      data pout/ 1.0, 2.0, 3.0, 5.0, 7.0                          &
               ,10.0,20.0,30.0,50.0,70.0,100.0,125.0,150.0,175.0  &
               ,200.0,225.0,250.0,275.0,300.0,325.0,350.0,375.0   &
               ,400.0,425.0,450.0,475.0,500.0,525.0,550.0,575.0   &
               ,600.0,625.0,650.0,675.0,700.0,725.0,750.0,775.0   &
               ,800.0,825.0,850.0,875.0,900.0,925.0,950.0,975.0   &
               ,1000.0/
!
      data rad/6.371e6/
!
      nx2=nint(360./rgrdsz)
      my2=nint(180./rgrdsz+1.)
      nxs=nint((rlon2-rlon1)/rgrdsz+1.)
      mys=nint((rlat2-rlat1)/rgrdsz+1.)
!
      if (myrank.eq.0) print*,'allocate rsmoutp'
      allocate (rsmoutp(nx2,my2))
      allocate (xr(nx2),yr(my2))
      if (myrank.eq.0) print*,'allocate all gfs variables for RSM'
      if (myrank.eq.0) then
        x1=nint(  (rlon1-0.)     /rgrdsz+1.)
        x2=nint(  (rlon2-0.)     /rgrdsz+1.)
        y1=nint(  (rlat1-(-90.)) /rgrdsz+1.)
        y2=nint(  (rlat2-(-90.)) /rgrdsz+1.)
        print*,'nxs, mys=',nxs,mys
        print*,'rlon1, rlon2, rlat1, rlat2, rgrdsz',                &
                  rlon1, rlon2, rlat1, rlat2, rgrdsz       
        if ( (x2-x1+1).ne.nxs  .or. (y2-y1+1).ne.mys) then
          print*,'wrong dimension for output in rsmout.f90'
          print*,'check: mod(360.,rgrdsz) should be zero'
!          call mpe_finalize
          stop
        endif
        allocate (                                                    &
     &     temp_gfs(nxs,mys,lpout),spfh_gfs(nxs,mys,lpout),           &
     &     clwr_gfs(nxs,mys,lpout),ozon_gfs(nxs,mys,lpout),           &
     &     rain_gfs(nxs,mys,lpout),qice_gfs(nxs,mys,lpout),           &
     &     snow_gfs(nxs,mys,lpout),grpl_gfs(nxs,mys,lpout),           &
     &     geop_gfs(nxs,mys),u_gfs(nxs,mys,lpout),v_gfs(nxs,mys,lpout), &
     &     tg_gfs(nxs,mys),smc_gfs(nxs,mys,km),snr_gfs(nxs,mys),         &
     &     stc_gfs(nxs,mys,km),ice_gfs(nxs,mys),                      &
!yj2019     
     &     terr_gfs(nxs,mys),slmsk_gfs(nxs,mys),                      &
     &     stat=ierr)
        if (ierr/=0) stop "rsmout: allocate fail"
        if (ierr/=0) call mpe_finalize               
          temp_gfs=0.0
          spfh_gfs=0.0
          clwr_gfs=0.0
          rain_gfs=0.0
          qice_gfs=0.0
          snow_gfs=0.0
          grpl_gfs=0.0
          ozon_gfs=0.0
          geop_gfs=0.0
          u_gfs=0.0
          v_gfs=0.0
          tg_gfs=0.0
          smc_gfs=0.0
          snr_gfs=0.0
          stc_gfs=0.0
          ice_gfs=0.0
          terr_gfs=0.0
          slmsk_gfs=0.0
      endif
      if (myrank.eq.0) print*,'after allocate all gfs variables for RSM'
!
      nxmy = nx*my
! initializing 
      pllp=0.
      bt1=0.
      bt2=0.
      glob=0.
      pres3d=0.
      soil_xy=0.
!
      if(myrank .eq. 0) print*,'   in rsmout for tau= ',itau
!
!
      if(ncld.ge.2)then
        ntrac=2
      else
        ntrac=ncld
      endif
      do n=1,ntrac
        nk=(n-1)*lev
        ngq = 0
        ntt = 0
        do jj =1,jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do k = 1, lev
            kk=nk+k
            do i = 1,nxj
              if ( qt(i,kk,jj) .lt. 0.0 )  then
                ngq = ngq + 1
                qt(i,kk,jj) = 0.0
              endif
            enddo
          enddo
!         ntt = ntt + nxj*lev
          ntt = ntt + nxjp(j)*lev
        enddo
        call mpe_global_sum(ngq,1,mpe_integer)
        call mpe_global_sum(ntt,1,mpe_integer)
        if( n.eq.1 .and. myrank.eq.0 .and. ngq.ne.0 )  print 901, ngq
        if( n.eq.2 .and. myrank.eq.0 .and. ngq.ne.0 )  print 902, ntt-ngq
      enddo

  901 format ( 1x, " *** warning: in outflds, there are grid points"  &
             , " with q < 0. , total number = ",i8,// )
  902 format ( 1x, " *** checking: in outflds, there are grid points" &
             , " with qc > 0. , total number = ",i8,// )
!
!                 compute interpolation coeffs
!
      ograv= 1.0/grav
!
!  generate structure variables
!
!  first, the log(p) quantities
!
      do k=1,lpout
        pkout(k)= log(pout(k))
      enddo
!
!  generate std height at p levels
!
      call geostd(lpout,pout,phistd)
!
      alaps = 0.0065
      rdg = rgas/grav
!
!  hydrostatic equation
!
      do jj =1,jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i=1,nxj
          phi(i,lev,jj)= cp*tt(i,lev,jj)*(pk2(i,lev,jj)-pk(i,lev,jj)) &
                       + sgeo(i,jj)
        enddo
        do k=lev-1,1,-1
          do i=1,nxj
            phi(i,k,jj)= phi(i,k+1,jj)+cp*(tt(i,k,jj)*(pk2(i,k,jj)-pk(i,k,jj)) &
                       + tt(i,k+1,jj)*(pk(i,k+1,jj)-pk2(i,k,jj)))
          enddo
        enddo
        do k=1,lev
          do i=1,nxj
            plog(i,k,jj)= log(plt(i,k,jj))
          enddo
        enddo
      enddo
!
!  compute sea level pressure and update pdiff
!  The method is based on one used by ecmwf, reseach manual 2 (1988)
!
!  llts layer's temperature is used to derive an alternative
!  surface skin temperature
!
!      if( itau .gt. 0 )then
!
      llts = lev-5
!
      do jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i=1,nxj
          ttb  = tt(i,lev,jj)*pk(i,lev,jj)/(1.0+0.608*qt(i,lev,jj))
          ttp  = tt(i,llts,jj)*pk(i,llts,jj)/(1.0+0.608*qt(i,llts,jj))
          ttt1 = ttb + alaps*rdg*ttb*   &
              ((pt(i,jj)+ptop)/plt(i,lev,jj)-1.0)
          ttt2 = ttp + alaps*(phi(i,llts,jj)-sgeo(i,jj))/grav
          hld1(i,jj) = 0.25*ttt1 + 0.75*ttt2
          hld2(i,jj) = hld1(i,jj) + alaps*sgeo(i,jj)/grav
          if( sgeo(i,jj) .lt. 0.1 ) then
            anlslp = pt(i,jj) + ptop
          else if( hld1(i,jj) .le. 290.5 .and. hld2(i,jj) .gt. 290.5 ) then
            apha = rgas*(290.5-hld1(i,jj))/sgeo(i,jj)
            ttt = sgeo(i,jj)/(rgas*hld1(i,jj))
            anlslp = (pt(i,jj)+ptop)*exp( ttt*(1.0-0.5*apha*ttt+0.333333*  &
                     apha*ttt*apha*ttt) )
          else if( hld1(i,jj) .gt. 290.5 .and. hld2(i,jj) .gt. 290.5 ) then
            hld1(i,jj) = (hld1(i,jj)+290.5)*0.5
            anlslp = (pt(i,jj)+ptop)*exp( sgeo(i,jj)/(rgas*hld1(i,jj)) )
          else if( hld1(i,jj) .lt. 255.0 .and. hld2(i,jj) .lt. 255.0 ) then
            hld1(i,jj) = (hld1(i,jj)+255.0)*0.5
            anlslp = (pt(i,jj)+ptop)*exp( sgeo(i,jj)/(rgas*hld1(i,jj)) )
          else
            apha = alaps * rdg
            ttt = sgeo(i,jj)/(rgas*hld1(i,jj))
            anlslp = (pt(i,jj)+ptop)*exp( ttt*(1.0-0.5*apha*ttt+0.333333*  &
                     apha*ttt*apha*ttt) )
          endif
          pllp(i,jj) = anlslp - pt(i,jj)
          pdiff(i,jj)=pllp(i,jj)
        enddo
      enddo
!
!      else
!
!  pllp equal to pdiff at the initial time, tau=0
!
!      do jj = 1, jlistnum
!        j=jlist1(jj)
!        nxj=nxdef(j)
!        do i = 1, nxj
!          pllp(i,j) = pdiff(i,jj)
!        enddo
!      enddo
!
!      endif     ! end of ( itau .gt. 0 )
!
!      call mpe_unify(pllp,nx,my,2,mpe_double)
!      call mpe_unify(pt,nx,my,2,mpe_double)
!
!
!  obtain the the bottom pressure for the following interpolations
!  bt2 will be used in geoptential interpolations
!
      pl1000= log(1000.1)
!
      do jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i = 1, nxj
          slp(i,jj)= pt(i,jj)+pdiff(i,jj)
          if( slp(i,jj) .le. 1000.1) then
            pllp(i,jj) = 1000.1
          else
            pllp(i,jj) = slp(i,jj)
          endif
        enddo
!
        call geostd (nxjp(j),pllp(1,jj),bt2(1,jj))
!
        do i = 1, nxj
          pllp(i,jj) = log(pllp(i,jj))
        enddo
      enddo
!
!      call mpe_unify(pllp,nx,my,2,mpe_double)
!      call mpe_unify(bt2,nx,my,2,mpe_double)
!
!  temperature output
!
      do jj =1,jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
!
!  bottom boundary condition for temperature:
!  limits for 1000 mb boundary condition are (1) standard lapse rate
!  from bottom model level for terrain height and (2) layer mean temp
!  of slp and 1000 mb.  Boundary condition is linear combination
!  of these two cases weighted with subterrainean thickness.
!
        do k = 1, lev
          do i = 1,nxj
            tmp(i,k,jj) = tt(i,k,jj)*pk(i,k,jj)/(1.0+0.608*qt(i,k,jj))
          enddo
        enddo
!
        do i=1,nxj
          slp(i,jj) = max( slp(i,jj), pt(i,jj) )
          splog= log(slp(i,jj))
          dx= log(slp(i,jj)/(pt(i,jj)+ptop))
!          bx= tmp(i,lev,jj)*(1.0/sig(lev+1))**0.19023 + sgeo(i,jj)*0.00065
          bx= tmp(i,lev,jj) + sgeo(i,jj)*0.00065
          if ( sgeo(i,jj) .ge. 1000.0 )  then
            tmid = sgeo(i,jj)/(dx*rgas)
            ax = tmp(i,lev,jj)+ (tmid-tmp(i,lev,jj))*(splog-plog(i,lev,jj)) &
                / (log((slp(i,jj)+pt(i,jj)+ptop)*0.5)-plog(i,lev,jj))
          else
            ax = bx
          endif
          cx = (sgeo(i,jj)-1000.0)/(10000.0-1000.0)
          cx = max (0.0, min (1.0, cx))
          dx = cx*cx
          tsf = dx*ax + (1.0-dx)*bx
          tadia= tmp(i,lev,jj)*(slp(i,jj)/plt(i,lev,jj))**0.25
          tsf  = max( tmp(i,lev,jj), min( tsf,tadia ) )
          if ( slp(i,jj) .gt. 1000.1 )  then
            bt1(i,jj)= tsf
            t1000(i,jj) = tsf
          else
            bt1(i,jj)= tsf + (tsf-tmp(i,lev,jj))*(pl1000-splog)  &
                     / (splog - plog(i,lev,jj))
            t1000(i,jj) = bt1(i,jj)
          endif
        enddo
      enddo
!
!!      call mpe_unify(bt1,nx,my,2,mpe_double)
!      if(myrank.eq.0) then
!        print*,'check bottom temperature'
!        call qmaxn3(bt1,cdum1,cdum2,1,1,1,nx,my,1)
!      endif

!
      if(myrank.eq.0)print*,' rsmout : start tempout'
! yj replace tempout        
      do k = 1, lev+1
        tens(k) = 1.0
      end do
      tens(lev+1)= 0.0
      tens(lev)= 0.0
      call voterp(nx,my,my_max,lev,lpout,plog,pllp,tmp,bt1 &
                  ,pkout,pres3d,tens)

      do k=1,lpout
        call unify_reduceintp(nx,my,my_max,pres3d(1,1,k),glob)
        if(myrank.eq.0) then
!          call qmax2d(glob,1,1,nx,my)
          call xyintpo('gg',nx,my,'ga',nx2,my2,glob,rsmoutp,0 &
                        ,xr,yr,.true.)
!          print*,'x1,x2,nxs,mys,k=',x1,x2,nxs,mys,k
!          call qmax2d(rsmoutp,1,1,nx2,my2)
          jj=0
          do j=y1,y2
            jj=jj+1
            ii=0
            do i=x1,x2
              ii=ii+1
              temp_gfs(ii,jj,lpout-k+1)=rsmoutp(i,j)
            enddo
          enddo
!          print*,'k=,temp_gfs(2,2,k)=',k,temp_gfs(2,2,lpout-k+1)
!          print*,'====check temp_gfs after partail for temp===='
          call qmax2d(temp_gfs(1,1,lpout-k+1),1,1,nxs,mys)
        endif
      enddo    

!
!  moisture output as mixing ratio
!
        if(ncld.lt.3) then
          print*,'Err: ntrac must include spehum,O3,and clwr for RSM'
          stop
        endif
!
      do ntrac=1,ncld
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do k = 1, lev
            kk = (ntrac-1)*lev+k
            do i = 1,nxj
              tmp(i,k,jj)=qt(i,kk,jj)
            enddo
          enddo
          do i = 1,nxj
            bt1(i,jj)=tmp(i,lev,jj)
          enddo
        enddo
!!        call mpe_unify(bt1,nx,my,2,mpe_double)
        if(myrank.eq.0)print*,' rsmout : start shumout2'
        if(myrank.eq.0)print*,' now is ntrac = ',ntrac, &
!                               '(1:spfh,2:cw,3:o3)'
                    '(1:spfh,2:cw,3:rain,4:ice,5:snow,6:graupel,7:o3)'
! yj replace shumuot2 out
        call voterp(nx,my,my_max,lev,lpout,plog,pllp,tmp,bt1 &
             ,pkout,pres3d,tens)
      glob=0.
        do k=1,lpout
! relative humidity, clwr, and O3 must be smaller or equal 1.0
          do jj = 1, jlistnum
            j=jlist1(jj)
            nxj=nxdef_2d(j)
            do i=1,nxj
              pres3d(i,jj,k)= max(pres3d(i,jj,k),0.0)
            enddo !end --loop i
          enddo !end --loop jj
          call unify_reduceintp(nx,my,my_max,pres3d(1,1,k),glob)

         if(myrank.eq.0) then
          call xyintpo('gg',nx,my,'ga',nx2,my2,glob &
                       ,rsmoutp,0,xr,yr,.true.)
!          jj=0
!          do j=y1,y2
!            jj=jj+1
!            ii=0
!            do i=x1,x2
!              ii=ii+1
           if (ntrac.eq.1) spfh_gfs(:,:,lpout-k+1)=max(rsmoutp(x1:x2,y1:y2),0.0)
           if (ntrac.eq.2) clwr_gfs(:,:,lpout-k+1)=max(rsmoutp(x1:x2,y1:y2),0.0)
           if (ntrac.eq.3) rain_gfs(:,:,lpout-k+1)=max(rsmoutp(x1:x2,y1:y2),0.0)
           if (ntrac.eq.4) qice_gfs(:,:,lpout-k+1)=max(rsmoutp(x1:x2,y1:y2),0.0)
           if (ntrac.eq.5) snow_gfs(:,:,lpout-k+1)=max(rsmoutp(x1:x2,y1:y2),0.0)
           if (ntrac.eq.6) grpl_gfs(:,:,lpout-k+1)=max(rsmoutp(x1:x2,y1:y2),0.0)
           if (ntrac.eq.ncld) ozon_gfs(:,:,lpout-k+1)=max(rsmoutp(x1:x2,y1:y2),0.0)
!            enddo
!          enddo
          if (ntrac.eq.1)call qmax2d(spfh_gfs(1,1,lpout-k+1),1,1,nxs,mys)
          if (ntrac.eq.2)call qmax2d(clwr_gfs(1,1,lpout-k+1),1,1,nxs,mys)
!          if (ntrac.eq.3)call qmax2d(ozon_gfs(1,1,lpout-k+1),1,1,nxs,mys)
          if (ntrac.eq.3)call qmax2d(rain_gfs(1,1,lpout-k+1),1,1,nxs,mys)
          if (ntrac.eq.4)call qmax2d(qice_gfs(1,1,lpout-k+1),1,1,nxs,mys)
          if (ntrac.eq.5)call qmax2d(snow_gfs(1,1,lpout-k+1),1,1,nxs,mys)
          if (ntrac.eq.6)call qmax2d(grpl_gfs(1,1,lpout-k+1),1,1,nxs,mys)
          if (ntrac.eq.ncld)call qmax2d(ozon_gfs(1,1,lpout-k+1),1,1,nxs,mys)
         endif ! end -- myrank
        enddo !end --loop k
      enddo    !end --loop ntrac

!
!  geopotential height output
!
        do jj =1,jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do i=1,nxj
            if(slp(i,jj).le.1000.1) then
              bt1(i,jj)= t1000(i,jj)*rgas*log(slp(i,jj)*0.001)
            else
              bt1(i,jj) = min( sgeo(i,jj), 0.0 )
            endif
          enddo
!
!  compute standard geopotentials, subtract them from sigma
!  level values. vertical interpolation will be done on
!  these deviation from standard values
!
          do k=1,lev
            call geostd (nxjp(j),plt(1,k,jj),wrk1(1,k))
            do i=1,nxj
              tmp(i,k,jj) = ograv*phi(i,k,jj)-wrk1(i,k)
            enddo
          enddo
!
!  bottom boundary conditions
!
          do i=1,nxj
            bt1(i,jj) = ograv*bt1(i,jj) - bt2(i,jj)
          enddo
        enddo
!!        call mpe_unify(bt1,nx,my,2,mpe_double)
!
        if(myrank.eq.0)print*,' rsmout : start geopout'
! yj replace geopout # output unit [gpm]
        call voterp(nx,my,my_max,lev,lpout,plog,pllp,tmp,bt1 &
            ,pkout,pres3d,tens)
      glob=0.
      do k=lpout,lpout
!yj      do k=1,lpout
        do  jj=1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do  i=1,nxj
         pres3d(i,jj,k)= pres3d(i,jj,k)+phistd(k)
        enddo
        enddo
        call unify_reduceintp(nx,my,my_max,pres3d(1,1,k),glob)
!       call smth9(nx,my,glob,glob2,1)
!       glob=glob2
!       call smth9(nx,my,glob,glob2,2)
!       glob=glob2

        if(myrank.eq.0) then
          call xyintpo('gg',nx,my,'ga',nx2,my2,glob &
                       ,rsmoutp,0,xr,yr,.true.)
          geop_gfs(:,:)=rsmoutp(x1:x2,y1:y2)
          call qmax2d(geop_gfs(1,1),1,1,nxs,mys)
        endif
      enddo    

! wind --- u and v
        do jj =1,jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          do i=1,nxj
            bt1(i,jj)= ut(i,lev,jj)
            bt2(i,jj)= vt(i,lev,jj)
          enddo
        enddo
!        call mpe_unify(bt1,nx,my,2,mpe_double)
!        call mpe_unify(bt2,nx,my,2,mpe_double)
        if(myrank.eq.0)print*,' rsmout : start windout'
!  yj replace windout
! first: do the u components
      tmp=ut
      call voterp(nx,my,my_max,lev,lpout,plog,pllp,tmp,bt1 &
                 ,pkout,pres3d,tens)
!      if( lreduce.eq.1 ) call reduceintp (bt1,nxdef,nx,my)
!      do i =1,nx
!        do j =1,my
!          pklzl(i,j)=pllp(i,j)
!        enddo
!      enddo
!      if( lreduce.eq.1 ) call reduceintp(pklzl,nxdef,nx,my)

      do k=1,lpout
!        
!  below ground level extrapolate surface wind downward
!  deweight wind with cos latitude, earth radius
!
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          xxx= rad/cosl(j)
          do i=1,nxj
           if(pkout(k).gt.pllp(i,jj)) pres3d(i,jj,k)= bt1(i,jj)
           pres3d(i,jj,k)=pres3d(i,jj,k)*xxx
          enddo
        enddo
        call unify_reduceintp(nx,my,my_max,pres3d(1,1,k),glob)
!        do i=1,nxmy
!          if(pkout(k).gt.pklzl(i,1)) pres3d(i,1,k)= bt1(i,1)
!        enddo
!
!  deweight wind with cos latitude, earth radius
!
!        glob=0.
!        do j=1,my
!          xx= rad/cosl(j)
!          do i=1,nx
!            glob(i,j)= pres3d(i,j,k)*xx
!          enddo
!        enddo 
! wind flag is "1" in xyintpo        
        if(myrank.eq.0) then
        call xyintpo('gg',nx,my,'ga',nx2,my2,glob &
                     ,rsmoutp,1,xr,yr,.true.)
          u_gfs(:,:,lpout-k+1)=rsmoutp(x1:x2,y1:y2)
          call qmax2d(u_gfs(1,1,lpout-k+1),1,1,nxs,mys)
        endif
      enddo  
!  now the v components
      tmp=vt
      call voterp(nx,my,my_max,lev,lpout,plog,pllp,tmp,bt2 &
                 ,pkout,pres3d,tens)
!      if( lreduce.eq.1 ) call reduceintp (bt2,nxdef,nx,my)
!
      do k=1,lpout
!        
!  below ground level extrapolate surface wind downward
!  deweight wind with cos latitude, earth radius
!
        do jj = 1, jlistnum
          j=jlist1(jj)
          nxj=nxdef_2d(j)
          xxx= rad/cosl(j)
          do i=1,nxj
           if(pkout(k).gt.pllp(i,jj)) pres3d(i,jj,k)= bt2(i,jj)
           pres3d(i,jj,k)=pres3d(i,jj,k)*xxx
          enddo
        enddo
        call unify_reduceintp(nx,my,my_max,pres3d(1,1,k),glob)
!        call qmaxn3(glob,cdum1,cdum2,1,1,1,nx,my,1)

!        do i=1,nxmy
!          if(pkout(k).gt.pklzl(i,1)) pres3d(i,1,k)=bt2(i,1)
!        enddo
!        glob=0.
!        do j=1,my
!          xx= rad/cosl(j)
!          do i=1,nx
!            glob(i,j)= pres3d(i,j,k)*xx
!          enddo
!        enddo
! wind flag is "1" in xyintpo        
        if(myrank.eq.0) then
        call xyintpo('gg',nx,my,'ga',nx2,my2,glob &
                     ,rsmoutp,1,xr,yr,.true.)
!          jj=0
!          do j=y1,y2
!            jj=jj+1
!            ii=0
!            do i=x1,x2
!              ii=ii+1
!              v_gfs(ii,jj,lpout-k+1)=rsmoutp(i,j)
!            enddo
!          enddo
          v_gfs(:,:,lpout-k+1)=rsmoutp(x1:x2,y1:y2)
          call qmax2d(v_gfs(1,1,lpout-k+1),1,1,nxs,mys)
        endif
      enddo  
!
!output the terr geopotential height,
! turn off by yjchen2019July, terr is not used and some index is wrong below
       do jj=1, jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
         do i =1,nxj
           hld1(i,jj)=sgeo(i,jj)*ograv
         enddo
       enddo  
!byl       call mpe2d_unify(glob,hld1)
!!      do j=1,my
!!       nxj=nxdef(j)
!!      do i=1,nxj
!!        glob(i,j)=glob(i,j)*ograv
!!      enddo
!!      enddo
!byl       if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)      
       call unify_reduceintp(nx,my,my_max,hld1,glob)
       if(myrank.eq.0) then
       call xyintpo('gg',nx,my,'ga',nx2,my2,glob &
                   ,rsmoutp,0,xr,yr,.true.)
!!yj2019         
!         print*,'rsmout output terrain gfs'
!          jj=0
!          do j=y1,y2
!            jj=jj+1
!            ii=0
!            do i=x1,x2
!              ii=ii+1
!              terr_gfs(ii,jj)=max(rsmoutp(i,j),0.0)
!            enddo
!          enddo
!          call qmaxn3(terr_gfs(1,1),cdum1,cdum2,1,1,1,nxs,mys,1)
           terr_gfs(:,:)=max(rsmoutp(x1:x2,y1:y2),0.0)
           call qmax2d(terr_gfs(1,1),1,1,nxs,mys)
        endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  !!!!!!!!!!!!!!!---below output the surface file for rsm (nsig2)---!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
      do jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i = 1, nxj
          soil_xy(i,jj,1) = smc(i,1,jj)
          soil_xy(i,jj,2) = smc(i,2,jj)
          soil_xy(i,jj,3) = smc(i,3,jj)
          soil_xy(i,jj,4) = smc(i,4,jj)
          soil_xy(i,jj,5) = stc(i,1,jj)
          soil_xy(i,jj,6) = stc(i,2,jj)
          soil_xy(i,jj,7) = stc(i,3,jj)
          soil_xy(i,jj,8) = stc(i,4,jj)
        enddo
      enddo
!
! output some 2-dimension veriable to dmsfile
!
      if(myrank .eq. 0) print *,' rsmout : output surface file'
!replace out2d
! need mpe unify
! ***tg***
!byl      call mpe2d_unify(glob,tg)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)      
      hld1=tg
      call unify_reduceintp(nx,my,my_max,hld1,glob)
      if(myrank.eq.0) then
      call xyintpo('gg',nx,my,'ga',nx2,my2,glob &
                  ,rsmoutp,0,xr,yr,.true.)
!          jj=0
!          do j=y1,y2
!            jj=jj+1
!            ii=0
!            do i=x1,x2
!              ii=ii+1
!              tg_gfs(ii,jj)=rsmoutp(i,j)
!            enddo
!          enddo
          tg_gfs(:,:)=rsmoutp(x1:x2,y1:y2)
          call qmax2d(tg_gfs(1,1),1,1,nxs,mys)
      endif
! ***0-0.1m soil moisture content [fraction]***
!byl      call mpe2d_unify(glob,soil_xy(1,1,1))
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)      
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,1),glob)
      if(myrank.eq.0) then
      call xyintpo('gg',nx,my,'ga',nx2,my2,glob &
                  ,rsmoutp,0,xr,yr,.true.)
!          jj=0
!          do j=y1,y2
!            jj=jj+1
!            ii=0
!            do i=x1,x2
!              ii=ii+1
!              smc_gfs(ii,jj,1)=rsmoutp(i,j)
!            enddo
!          enddo
          smc_gfs(:,:,1)=rsmoutp(x1:x2,y1:y2)
          call qmax2d(smc_gfs(1,1,1),1,1,nxs,mys)
      endif
! ***0.1-0.4m soil moisture content [fraction]***
!byl      call mpe2d_unify(glob,soil_xy(1,1,2))
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)      
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,2),glob)
      if(myrank.eq.0) then
      call xyintpo('gg',nx,my,'ga',nx2,my2,glob &
                  ,rsmoutp,0,xr,yr,.true.)
          smc_gfs(:,:,2)=rsmoutp(x1:x2,y1:y2)
          call qmax2d(smc_gfs(1,1,2),1,1,nxs,mys)
      endif
! ***0.4-1m soil moisture content [fraction]***
!byl      call mpe2d_unify(glob,soil_xy(1,1,3))
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)      
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,3),glob)
      if(myrank.eq.0) then
      call xyintpo('gg',nx,my,'ga',nx2,my2,glob &
                  ,rsmoutp,0,xr,yr,.true.)
!          jj=0
!          do j=y1,y2
!            jj=jj+1
!            ii=0
!            do i=x1,x2
!              ii=ii+1
!              smc_gfs(ii,jj,3)=rsmoutp(i,j)
!            enddo
!          enddo
          smc_gfs(:,:,3)=rsmoutp(x1:x2,y1:y2)
          call qmax2d(smc_gfs(1,1,3),1,1,nxs,mys)
      endif
! ***below 1m soil moisture content [fraction]***
!byl      call mpe2d_unify(glob,soil_xy(1,1,4))
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)      
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,4),glob)
      if(myrank.eq.0) then
      call xyintpo('gg',nx,my,'ga',nx2,my2,glob &
                  ,rsmoutp,0,xr,yr,.true.)
          smc_gfs(:,:,4)=rsmoutp(x1:x2,y1:y2)
          call qmax2d(smc_gfs(1,1,4),1,1,nxs,mys)
      endif
! ***snr***
!byl      call mpe2d_unify(glob,snr)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)      
      hld1=snr
      call unify_reduceintp(nx,my,my_max,hld1,glob)
      if(myrank.eq.0) then
      call xyintpo('gg',nx,my,'ga',nx2,my2,glob &
                  ,rsmoutp,0,xr,yr,.true.)
!          jj=0
!          do j=y1,y2
!            jj=jj+1
!            ii=0
!            do i=x1,x2
!              ii=ii+1
!              snr_gfs(ii,jj)=rsmoutp(i,j)
!            enddo
!          enddo
          snr_gfs(:,:)=rsmoutp(x1:x2,y1:y2)
          call qmax2d(snr_gfs(1,1),1,1,nxs,mys)
      endif
! ***0-0.1m soil temperature***
!byl      call mpe2d_unify(glob,soil_xy(1,1,5))
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)      
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,5),glob)
      if(myrank.eq.0) then
      call xyintpo('gg',nx,my,'ga',nx2,my2,glob &
                  ,rsmoutp,0,xr,yr,.true.)
!          jj=0
!          do j=y1,y2
!            jj=jj+1
!            ii=0
!            do i=x1,x2
!              ii=ii+1
!              stc_gfs(ii,jj,1)=rsmoutp(i,j)
!            enddo
!          enddo
          stc_gfs(:,:,1)=rsmoutp(x1:x2,y1:y2)
          call qmax2d(stc_gfs(1,1,1),1,1,nxs,mys)
      endif
! ***0.1-0.4m soil temperature***
!byl      call mpe2d_unify(glob,soil_xy(1,1,6))
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)      
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,6),glob)
      if(myrank.eq.0) then
      call xyintpo('gg',nx,my,'ga',nx2,my2,glob &
                  ,rsmoutp,0,xr,yr,.true.)
          stc_gfs(:,:,2)=rsmoutp(x1:x2,y1:y2)
          call qmax2d(stc_gfs(1,1,2),1,1,nxs,mys)
      endif
! ***0.4-1m soil temperature***
!byl      call mpe2d_unify(glob,soil_xy(1,1,7))
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)      
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,7),glob)
      if(myrank.eq.0) then
      call xyintpo('gg',nx,my,'ga',nx2,my2,glob &
                  ,rsmoutp,0,xr,yr,.true.)
!          jj=0
!          do j=y1,y2
!            jj=jj+1
!            ii=0
!            do i=x1,x2
!              ii=ii+1
!              stc_gfs(ii,jj,3)=rsmoutp(i,j)
!            enddo
!          enddo
          stc_gfs(:,:,3)=rsmoutp(x1:x2,y1:y2)
          call qmax2d(stc_gfs(1,1,3),1,1,nxs,mys)
      endif
! ***below 1m soil temperature***
!byl      call mpe2d_unify(glob,soil_xy(1,1,8))
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)      
      call unify_reduceintp(nx,my,my_max,soil_xy(1,1,8),glob)
      if(myrank.eq.0) then
      call xyintpo('gg',nx,my,'ga',nx2,my2,glob &
                  ,rsmoutp,0,xr,yr,.true.)
          stc_gfs(:,:,4)=rsmoutp(x1:x2,y1:y2)
          call qmax2d(stc_gfs(1,1,4),1,1,nxs,mys)
      endif
! ***land and sea mask***(land=1,sea=0)***
       glob=0.
       slmsk=0.
       do jj=1,jlistnum
         j=jlist1(jj)
         nxj=nxdef_2d(j)
       do i=1,nxj
         if(land(i,jj))slmsk(i,jj)=1.0
         if(ocean(i,jj))slmsk(i,jj)=0.0
       enddo
       enddo
!byl       call mpe2d_unify(glob,slmsk)
!byl       if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)
      hld1=slmsk
      call unify_reduceintp(nx,my,my_max,hld1,glob)
!yj2019
       if(myrank.eq.0) then
       call xyintpo('gg',nx,my,'ga',nx2,my2,glob &
                   ,rsmoutp,0,xr,yr,.true.)
            jj=0
            do j=y1,y2
              jj=jj+1
              ii=0
              do i=x1,x2
                ii=ii+1
                slmsk_gfs(ii,jj)=max(min(rsmoutp(i,j),1.0),0.0)
              enddo
            enddo
!           call qmaxn3(slmsk_gfs(1,1),cdum1,'--slmsk-----',1,1,1,nxs,mys,1)
!           slmsk_gfs(:,:)=rsmoutp(x1:x2,y1:y2)
           call qmax2d(slmsk_gfs(1,1),1,1,nxs,mys)
       endif
!! ***ice***(simk in RSM-csfcfld(:,13))
      glob=0.
!byl      call mpe2d_unify(glob,slmsk)
!byl      if( lreduce.eq.1 ) call reduceintp (glob,nxdef,nx,my)      
      hld1=cice
      call unify_reduceintp(nx,my,my_max,hld1,glob)
      if(myrank.eq.0) then
      call xyintpo('gg',nx,my,'ga',nx2,my2,glob &
                  ,rsmoutp,0,xr,yr,.true.)
      ice_gfs(:,:)=rsmoutp(x1:x2,y1:y2)
      call qmax2d(ice_gfs(1,1),1,1,nxs,mys)
      endif
! ***
#if defined(CWB_MPMD) || defined(CWBSUM)
      if(myrank.eq.0) then
      print*,'before into send_data in rsmout.f90'
      call send_data(float(itau),nxs,mys,km,lpout,ncld,temp_gfs,spfh_gfs &
     &              ,clwr_gfs,rain_gfs,qice_gfs,snow_gfs,grpl_gfs   &
     &              ,ozon_gfs,geop_gfs,u_gfs,v_gfs      &   
     &              ,tg_gfs,smc_gfs,snr_gfs,stc_gfs,ice_gfs     &
!yj2019     
     &              ,terr_gfs,slmsk_gfs)
      endif
#else
      if(myrank.eq.0) then
      print*,'before into wrte_data in rsmout.f90'
      call wrte_data(idtg,float(itau),nxs,mys,km,lpout,ncld,temp_gfs,spfh_gfs &
     &              ,clwr_gfs,rain_gfs,qice_gfs,snow_gfs,grpl_gfs   &
     &              ,ozon_gfs,geop_gfs,u_gfs,v_gfs      &   
     &              ,tg_gfs,smc_gfs,snr_gfs,stc_gfs,ice_gfs     &
!yj2019     
     &              ,terr_gfs,slmsk_gfs) 
      endif
#endif      

!      
      if (myrank.eq.0) then
      deallocate ( slmsk_gfs,terr_gfs,ice_gfs,stc_gfs, & 
                   snr_gfs,smc_gfs,tg_gfs,v_gfs,u_gfs,  &
                   geop_gfs,ozon_gfs,grpl_gfs,snow_gfs, &
                   qice_gfs,rain_gfs,clwr_gfs,spfh_gfs, &
                   temp_gfs, stat=ierr)
      if (ierr/=0) stop "rsmout: deallocate fail"
      endif
      deallocate (xr, yr, stat=ierr)                                
      deallocate (rsmoutp, stat=ierr)                                
      if (ierr/=0) stop "rsmout: deallocate fail rsmoutp"
!      
#endif
    end subroutine rsmout
