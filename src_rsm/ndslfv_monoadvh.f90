      subroutine ndslfv_monoadvh(ddtemp,qvadv,pten,vdzonl,vdmerd     &
                , lonsperlat,deltim,xy,levs)
      use param
      use const, only : RTYPE
      implicit none
      real(kind=RTYPE) pten(nx,levs,my_max)
      real(kind=RTYPE) ddtemp(nx,levs,my_max),qvadv(nx,levs*ncld,my_max)
      real(kind=RTYPE) vdmerd(nx,levs,my_max),vdzonl(nx,levs,my_max)
      integer,intent(in):: lonsperlat(my)
      real(kind=RTYPE),   intent(in):: deltim
      integer xy,levs
      
      if(xy .eq.   0) call ndslfv_monoadvh2(ddtemp,qvadv,pten,vdzonl    &
                             ,vdmerd,lonsperlat,deltim,levs)
      if(xy .gt. 0.5) call ndslfv_monoadvh2_xy(ddtemp,qvadv,pten,vdzonl &
                             ,vdmerd,lonsperlat,deltim,levs)
      if(xy .lt.-0.5) call ndslfv_monoadvh2_yx(ddtemp,qvadv,pten,vdzonl &
                             ,vdmerd,lonsperlat,deltim,levs)
!      xy = -1 * xy
      return
      end
!
      subroutine ndslfv_monoadvh_fgnl(vdzonl,vdmerd,ddtemp            &
                , lonsperlat,deltim,xy,levs,nvars,forward)
      use param
      use const, only : RTYPE
      implicit none
      real(kind=RTYPE) vdmerd(nx,levs,my_max),vdzonl(nx,levs,my_max)
      real(kind=RTYPE) ddtemp(nx,levs,my_max)
      integer,intent(in):: lonsperlat(my)
      real(kind=RTYPE),   intent(in):: deltim
      integer xy,levs,nvars
      logical forward

      if(xy .eq.   0) call ndslfv_monoadvh2_fgnl(vdzonl,vdmerd,ddtemp    &
                             ,lonsperlat,deltim,levs,nvars,forward)
      if(xy .gt. 0.5) call ndslfv_monoadvh2_fgnl_xy(vdzonl,vdmerd,ddtemp &
                             ,lonsperlat,deltim,levs,nvars,forward)
      if(xy .lt.-0.5) call ndslfv_monoadvh2_fgnl_yx(vdzonl,vdmerd,ddtemp &
                             ,lonsperlat,deltim,levs,nvars,forward)
!      xy = -1 * xy
      return
      end
!
      subroutine ndslfv_monoadvh2(ddtemp,qvadv,pten,vdzonl,vdmerd    &
                , lonsperlat,deltim,levs)
!
! a routine to do non-iteration semi-Lagrangain advection
! considering advection  with monotonicity in interpolation
! contact: hann-ming henry juang
! program log
! 2011 02 20 : henry juang, created for ndsl advection
! 2013 06 20 : Henry Juang correct wind direction for north-south advection
!
      use param
      use grid
      use index
      use const
      use rank

      implicit none

      real(kind=RTYPE) pten(nx,levs,my_max)
      real(kind=RTYPE) ddtemp(nx,levs,my_max),qvadv(nx,levs,ncld,my_max)
      real(kind=RTYPE) vdmerd(nx,levs,my_max),vdzonl(nx,levs,my_max)
!      real(kind=RTYPE) plev(lonfull,lev+1)
!      integer,intent(in):: global_lats_a(my)
      integer,intent(in):: lonsperlat(my)
      real(kind=RTYPE),   intent(in):: deltim

      real(kind=RTYPE)      uulon(lonfull,levs,latpart)
      real(kind=RTYPE)      vvlon(lonfull,levs,latpart)
      real(kind=RTYPE)      qqlon(lonfull,levs*ndslhvar,latpart)
      real(kind=RTYPE)      rrlon(lonfull,levs*ndslhvar,latpart)

      real(kind=RTYPE)      vvlat(latfull,levs,lonpart)
      real(kind=RTYPE)      qqlat(latfull,levs*ndslhvar,lonpart)
      real(kind=RTYPE)      rrlat(latfull,levs*ndslhvar,lonpart)
      real(kind=RTYPE)      xr    (lonfull,levs)
      real(kind=RTYPE)      xcp   (lonfull,levs)
      real(kind=RTYPE)      sumrq (lonfull,levs)
      real(kind=RTYPE)      xkappa(lonfull,levs)
      real(kind=RTYPE)      rdt2, rkt, pi, cons0, cons1, rma, rm2a

!      logical   lprint

      integer mono,mass,levs
      integer nlevs,nvars!,levh
      integer i,j,n,k,lon,lan,lat,lons_lat,irc,kk,KL
      integer kuu, kvv, kqq, ktt, kup, nqq
      integer ku , kv , kq , kt , kp
!
!      lprint = .false.

!      if( lprint ) print *,' enter ndslfv_monoadvh '
!
      mono  = 1
      mass  = 0
      cons0 = 0.0
      cons1 = 1.0
      qqlon = 0.
      rrlon = 0.
      uulon = 0.
      vvlon = 0.
!
!      levh = ncld * lev
!
!      kuu = 1
!      kvv = kuu + lev
!      kuu = kvv + lev
      kuu = 1
      kvv = kuu + levs
      ktt = kvv + levs
      kup = ktt + levs
!!      kqq = kup + levs
      kqq = ktt + levs

      nvars = ndslhvar
      nlevs = nvars * levs
!
      rdt2 = 0.5 / deltim
!
! =================================================================
!   prepare wind and variable in flux form with gaussina weight
! =================================================================
!

!$omp parallel do                                                      &
!$omp private(lan,lat,lons_lat,rma,rm2a,i,n,kk,k,kt,kv,ku,kp,kq) &
!$omp schedule(dynamic)
      do lan=1,jlistnum

!        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lat = jlist1(lan)
        lons_lat = lonsperlat(lat)
!        rma  = 1. / cosglat(lat) / con_rerth
!        rm2a = rma / cosglat(lat)
        rma  = 1. / cosl(lat) 
        rm2a = rma / cosl(lat)
!
! wind at time step n
        do k=1,levs
          kk=levs-k+1
          do i=1,lons_lat
!ch         uulon(i,k,lan) = ut(i,kk,lan) * rm2a
            uulon(i,k,lan) = ut_sl(i,k ,lan) * rm2a
!hmhj use real(kind=RTYPE) wind
!            vvlon(i,k,lan) = grid_gr(ilan,kvg) * rma
!hmhj use virtual wind
!            vvlon(i,k,lan) = grid_gr(ilan,kvg) / con_rerth
!for CWB GFS
!ch         vvlon(i,k,lan) = vt(i,kk,lan) * rma
            vvlon(i,k,lan) = vt_sl(i,k ,lan) * rma
          enddo
        enddo
!        if( lprint ) then
!          call mymaxmin(uulon(1,1,lan),lons_lat,lonfull,1,' uu1 in deg')
!          call mymaxmin(uulon(1,5,lan),lons_lat,lonfull,1,' uu5 in deg')
!          call mymaxmin(vvlon(1,1,lan),lons_lat,lonfull,1,' vv in deg')
!        endif
!
! u v t at n-1
!
        do k=1,levs
          ku=kuu+k-1
          kv=kvv+k-1
          kt=ktt+k-1
          kp=kup+k-1
          kk=levs-k+1
!ch>
!       KL=lev-Llist(k)+1
!ch<
          do i=1,lons_lat
!ch         qqlon(i,ku,lan) = uum(i,kk,lan)
!ch         qqlon(i,kv,lan) = vvm(i,kk,lan)
!ch         qqlon(i,kt,lan) = ttm(i,kk,lan)
!!            qqlon(i,kp,lan) = dsigma(kk,1)*ptp(i,lan)+dsigma(kk,2)

            qqlon(i,ku,lan) = vdzonl(i,k,lan)
            qqlon(i,kv,lan) = vdmerd(i,k,lan)
            qqlon(i,kt,lan) = ddtemp(i,k,lan)
!!            qqlon(i,kp,lan) = pten(i,k,lan)
!!            qqlon(i,kp,lan) = dsigma(KL,1)*ptp_sl(i,lan)+dsigma(KL,2)
          enddo
        enddo
! rq at n-1
        do n=1,ncld
        do k=1,levs
          kq=kqq+k-1+(n-1)*levs
          kk=levs-k+1+(n-1)*levs
          do i=1,lons_lat
!ch         qqlon(i,kq,lan) = qm(i,kk,lan)
            qqlon(i,kq,lan) = qvadv(i,k,n,lan)
          enddo
        enddo
        enddo
! add surface pressure perturbation for removing resonance
!       rkt = con_g / ( con_rd * 300.0 )
!       do i=1,lons_lat
!         ilan=i+jlonf
!         qqlon(i,kuu,lan) = log(plev(i,1))+grid_gr(ilan,g_gz)*rkt
!       enddo
!
! save qqlon into n+1 for later as tendency
! z d t at n+1
!        do k=1,lev
!          kk=lev-k+1
!          ku=kuu+k-1
!          kv=kvv+k-1
!          kt=ktt+k-1
!          kug=g_d+k-1
!          kvg=g_z+k-1
!          ktg=g_t+k-1
!          do i=1,lons_lat
!            ilan=i+jlonf
!            grid_gr(ilan,kug) = qqlon(i,ku,lan)
!            grid_gr(ilan,kvg) = qqlon(i,kv,lan)
!            grid_gr(ilan,ktg) = qqlon(i,kt,lan)
!          enddo
!        enddo
! change h to theta
!
! ----- prepare xr, xcp, xkapa
!
!       xr    = cons0
!       xcp   = cons0
!       sumrq = cons0
!
!       do n=1,ntrac
!         nqq = kqq + (n-1)*levs
!         if( ri(n) .ne. cons0 .and. cpi(n) .ne. cons0 ) then
!           do k=1,lev
!             kq=nqq+k-1
!             do i=1,lons_lat
!               xr   (i,k) = xr   (i,k) + qqlon(i,kq,lan)*ri(n)
!               xcp  (i,k) = xcp  (i,k) + qqlon(i,kq,lan)*cpi(n)
!               sumrq(i,k) = sumrq(i,k) + qqlon(i,kq,lan)
!             enddo
!           enddo
!         endif
!       enddo
!       do k=1,lev
!         do i=1,lons_lat
!           xr (i,k)   = ( cons1 - sumrq(i,k) )*ri(0)  + xr (i,k)
!           xcp(i,k)   = ( cons1 - sumrq(i,k) )*cpi(0) + xcp(i,k)
!           xkappa(i,k) = xr(i,k) / xcp(i,k)
!         enddo
!       enddo
!
!       do k=1,lev
!         kt=ktt+k-1
!         ku=kuu+k-1
!         do i=1,lons_lat
!           pi = (qqlon(i,ku,lan)*0.005)**xkappa(i,k)
!           qqlon(i,kt,lan) = qqlon(i,kt,lan)/pi
!         enddo
!       enddo
!
! rq  no need for tendency
!
!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        ilan=1+jlonf
!        call mymaxmin(grid_gr(ilan,g_ttm),lons_lat,lonfull,1,' n-1 t ')
!        call mymaxmin(grid_gr(ilan,g_rrm ),lons_lat,lonfull,1,' n-1 q ')
!        call mymaxmin(grid_gr(ilan,g_zem),lons_lat,lonfull,1,' n-1 z ')
!        call mymaxmin(grid_gr(ilan,g_dim),lons_lat,lonfull,1,' n-1 d ')
!        print *,' ------------------------------------------- '
!        call mymaxmin(qqlon(1,ktt,lan),lons_lat,lonfull,1,' red t ')
!        call mymaxmin(qqlon(1,kqq,lan),lons_lat,lonfull,1,' red q ')
!        call mymaxmin(qqlon(1,kvv,lan),lons_lat,lonfull,1,' red z ')
!        call mymaxmin(qqlon(1,kuu,lan),lons_lat,lonfull,1,' red d ')
!        endif

!       call cyclic_cell_intpx(levs,lons_lat,lonf,uulon(1,1,lan))

!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        call mymaxmin(qqlon(1,ktt,lan),lonfull,lonfull,1,' full t ')
!        call mymaxmin(qqlon(1,kqq,lan),lonfull,lonfull,1,' full q ')
!        call mymaxmin(qqlon(1,kvv,lan),lonfull,lonfull,1,' full z ')
!        call mymaxmin(qqlon(1,kuu,lan),lonfull,lonfull,1,' full d ')
!        endif

        rrlon(:,:,lan) = qqlon(:,:,lan)
!
! first set positive advection in east-west direction
!
        call cyclic_cell_massadvx(lons_lat,lonfull,levs,nvars,deltim,   &
                         uulon(1,1,lan),rrlon(1,1,lan),mass)
!       call cyclic_mono_advectx (lonfull,levs,nvars,deltim,               &
!    &                   uulon(1,1,lan),rrlon(1,1,lan),mono)

        call cyclic_cell_intpx(nlevs,lons_lat,lonfull,rrlon(1,1,lan))
        call cyclic_cell_intpx(nlevs,lons_lat,lonfull,qqlon(1,1,lan))
        call cyclic_cell_intpx(levs,lons_lat,lonfull,vvlon(1,1,lan))

!        if( lprint ) then
!        print *,' done cyclic_massadvx with mass= ',mass
!        print *,' ------------------------------------------- '
!        call mymaxmin(rrlon(1,ktt,lan),lonfull,lonfull,1,' advx t ')
!        call mymaxmin(rrlon(1,kqq,lan),lonfull,lonfull,1,' advx q ')
!        call mymaxmin(rrlon(1,kvv,lan),lonfull,lonfull,1,' advx z ')
!        call mymaxmin(rrlon(1,kuu,lan),lonfull,lonfull,1,' advx d ')
!        print *,' done the first x adv for lan=',lan
!        print *,' =========================================== '
!        endif

      enddo
!$omp end parallel do

! ---------------------------------------------------------------------
! mpi para from east-west full grid to north-south full grid
! ---------------------------------------------------------------------
!
! para vvlon, qqlon, and rrlon to vvlat, qqlat, rrlat
!       if( lprint ) print *,' ndslfv_advect transport from we to ns '

       call para_we2ns(vvlon,vvlat,levs,my)
       call para_we2ns(rrlon,rrlat,nlevs,my)
       call para_we2ns(qqlon,qqlat,nlevs,my)

!       if( lprint ) then
!       print *,' ------------ after we2ns ---------------------- '
!       do lon=1,mylonlen
!        print *,'  lon=',lon
!        call mymaxmin(vvlat(1,1  ,lon),latfull,latfull,1,' we2ns v')
!        call mymaxmin(rrlat(1,kqq,lon),latfull,latfull,1,' we2ns r')
!        call mymaxmin(qqlat(1,kqq,lon),latfull,latfull,1,' we2ns q')
!       enddo
!       endif
!
! ---------------------------------------------------------------------
! -------------- in north-soutn great circle -------------------
! ---------------------------------------------------------------------

!      if ( myrank .eq. 0 ) lprint = .true.
!       if( lprint ) then
!       print *,' ndslfv_advect adv loop in y '
!       print *,' mylonlen=',mylonlen
!       endif

!$omp parallel do private(lon,k,j,ku,kv) &
!$omp schedule(dynamic)
       do lon=1,mylonlen
!
! convert wind before advy
        do k=1,levs
          ku=kuu+k-1
          kv=kvv+k-1
          do j=1,lathalf
            rrlat(j,ku,lon) = -rrlat(j,ku,lon)
            rrlat(j,kv,lon) = -rrlat(j,kv,lon)
            qqlat(j,ku,lon) = -qqlat(j,ku,lon)
            qqlat(j,kv,lon) = -qqlat(j,kv,lon)
          enddo
        enddo
!
!        if( lprint ) print *,' lon=',lon

! first set advection in north-south direction in great circle through two poles
!
!        call fixend_cell_massadvy(latfull,lathalf,levs,nvars,deltim, &
!                         vvlat(1,1,lon),rrlat(1,1,lon),mass)
        call cyclic_cell_massadvy(latfull,levs,nvars,deltim,         &
                         vvlat(1,1,lon),rrlat(1,1,lon),mass)
!       call cyclic_mono_advecty (latfull,levs,nvars,deltim,
!    &                   vvlat(1,1,lon),rrlat(1,1,lon),mono)

!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        call mymaxmin(rrlat(1,kuu,lon),latfull,latfull,1,' advry u ')
!        call mymaxmin(rrlat(1,kvv,lon),latfull,latfull,1,' advry v ')
!        call mymaxmin(rrlat(1,ktt,lon),latfull,latfull,1,' advry t ')
!        call mymaxmin(rrlat(1,kqq,lon),latfull,latfull,1,' advry q ')
!        call mymaxmin(rrlat(1,kvv,lon),latfull,latfull,1,' advry z ')
!        call mymaxmin(rrlat(1,kuu,lon),latfull,latfull,1,' advry d ')
!        endif
!
! second set advection in north-south direction in great circle through two poles
!
!        call fixend_cell_massadvy(latfull,lathalf,levs,nvars,deltim, &
!                         vvlat(1,1,lon),qqlat(1,1,lon),mass)
        call cyclic_cell_massadvy(latfull,levs,nvars,deltim,            &
                         vvlat(1,1,lon),qqlat(1,1,lon),mass)
!       call cyclic_mono_advecty (my,levs,nvars,deltim,
!    &                   vvlat(1,1,lon),qqlat(1,1,lon),mono)

!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        call mymaxmin(qqlat(1,kuu,lon),latfull,latfull,1,' advqy u ')
!        call mymaxmin(qqlat(1,kvv,lon),latfull,latfull,1,' advqy v ')
!        call mymaxmin(qqlat(1,ktt,lon),latfull,latfull,1,' advqy t ')
!        call mymaxmin(qqlat(1,kqq,lon),latfull,latfull,1,' advqy q ')
!        call mymaxmin(qqlat(1,kvv,lon),latfull,latfull,1,' advqy z ')
!        call mymaxmin(qqlat(1,kuu,lon),latfull,latfull,1,' advqy d ')
!        print *,' done with y at lon=',lon
!        endif
! convert wind back after advy
        do k=1,levs
          ku=kuu+k-1
          kv=kvv+k-1
          do j=1,lathalf
            rrlat(j,ku,lon) = -rrlat(j,ku,lon)
            rrlat(j,kv,lon) = -rrlat(j,kv,lon)
            qqlat(j,ku,lon) = -qqlat(j,ku,lon)
            qqlat(j,kv,lon) = -qqlat(j,kv,lon)
          enddo
        enddo

       enddo
!$omp end parallel do
!
! ----------------------------------------------------------------------
! mpi para from north-south direction to east-west direeectory
! ----------------------------------------------------------------------
!
! para qqlat and rrlat to qqlon and rrlon

!       if( lprint ) print *,' ndslfv_advect transport from ns to we '

       call para_ns2we(rrlat,rrlon,nlevs,my)
       call para_ns2we(qqlat,qqlon,nlevs,my)
!       if( lprint ) then
!       print *,' ------------ after ns2we ---------------------- '
!       do lan=1,jlistnum
!        print *,'  lan=',lan
!        call mymaxmin(rrlon(1,kqq,lan),lonfull,lonfull,1,' ns2we r')
!        call mymaxmin(qqlon(1,kqq,lan),lonfull,lonfull,1,' ns2we q')
!       enddo
!       endif

! ---------------------------------------------------------------
! ---------------- back to east-west direction ------------------
! ---------------------------------------------------------------
!      print *,' ndslfv_advect adv loop in x for last '

!$omp parallel do                                                &
!$omp private(lan,lat,lons_lat,i,k,kk,n,kt,kv,ku,kp,kq)    &
!$omp schedule(dynamic)

      do lan=1,jlistnum

!        lat = global_lats_a(ipt_jlistnum-1+lan)
        lat = jlist1(lan)
        lons_lat = lonsperlat(lat)
!
! mass conserving interpolation from full grid to reduced grid
!
        call cyclic_cell_intpx(nlevs,lonfull,lons_lat,qqlon(1,1,lan))
        call cyclic_cell_intpx(nlevs,lonfull,lons_lat,rrlon(1,1,lan))

!
! second set advection in x for the second of the pair
!
        call cyclic_cell_massadvx(lons_lat,lonfull,levs,nvars,deltim,   &
                         uulon(1,1,lan),qqlon(1,1,lan),mass)
!       call cyclic_mono_advectx (lonfull,levs,nvars,deltim,               &
!    &                   uulon(1,1,lan),qqlon(1,1,lan),mono)
!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        call mymaxmin(qqlon(1,kqq,lan),lonfull,lonfull,1,' adv x q ')
!        endif

        do k=1,nlevs
          do i=1,lons_lat
            qqlon(i,k,lan) = 0.5 * ( qqlon(i,k,lan) + rrlon(i,k,lan) )
          enddo
        enddo

!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        call mymaxmin(qqlon(1,kuu,lan),lons_lat,lonfull,1,' do redu u ')
!        call mymaxmin(qqlon(1,kvv,lan),lons_lat,lonfull,1,' do redu v ')
!        call mymaxmin(qqlon(1,ktt,lan),lons_lat,lonfull,1,' do redu t ')
!        call mymaxmin(qqlon(1,kqq,lan),lons_lat,lonfull,1,' do redu q ')
!        call mymaxmin(qqlon(1,kvv,lan),lons_lat,lonfull,1,' do redu z')
!        call mymaxmin(qqlon(1,kuu,lan),lons_lat,lonfull,1,' do redu d')
!        print *,' finish horizonatal advection at lan=',lan
!        endif

! change theta to h
!
! ----- prepare xr, xcp, xkapa
!
!       xr    = cons0
!       xcp   = cons0
!       sumrq = cons0
!
!       do n=1,ntrac
!         nqq = kqq + (n-1)*levs
!         if( ri(n) .ne. cons0 .and. cpi(n) .ne. cons0 ) then
!           do k=1,lev
!             kq=nqq+k-1
!             do i=1,lons_lat
!               xr   (i,k) = xr   (i,k) + qqlon(i,kq,lan)*ri(n)
!               xcp  (i,k) = xcp  (i,k) + qqlon(i,kq,lan)*cpi(n)
!               sumrq(i,k) = sumrq(i,k) + qqlon(i,kq,lan)
!             enddo
!           enddo
!         endif
!       enddo
!       do k=1,lev
!         do i=1,lons_lat
!           xr (i,k)   = ( cons1 - sumrq(i,k) )*ri(0)  + xr (i,k)
!           xcp(i,k)   = ( cons1 - sumrq(i,k) )*cpi(0) + xcp(i,k)
!           xkappa(i,k) = xr(i,k) / xcp(i,k)
!         enddo
!       enddo
!
!       do k=1,lev
!         kt=ktt+k-1
!         ku=kuu+k-1
!         do i=1,lons_lat
!           pi = (qqlon(i,ku,lan)*0.005)**xkappa(i,k)
!           qqlon(i,kt,lan) = qqlon(i,kt,lan)*pi
!         enddo
!       enddo
! u v t update at n+1
        do k=1,levs
          kk=levs-k+1
          ku=kuu+k-1
          kv=kvv+k-1
          kt=ktt+k-1
          kp=kup+k-1
!ch>
!       KL=lev-Llist(k)+1:w
!ch<
          do i=1,lons_lat
!ch         vdzonl(i,kk,lan) = (qqlon(i,ku,lan)-uum(i,kk,lan))*rdt2
!ch         vdmerd(i,kk,lan) = (qqlon(i,kv,lan)-vvm(i,kk,lan))*rdt2
!ch         ddtemp(i,kk,lan) = qqlon(i,kt,lan)
!ch         pten(i,kk,lan)    =(qqlon(i,kp,lan)-(dsigma(kk,1)          &
!ch                            *ptp(i,lan)+dsigma(kk,2)))*rdt2

            vdzonl(i,k ,lan) = qqlon(i,ku,lan)
            vdmerd(i,k ,lan) = qqlon(i,kv,lan)
            ddtemp(i,k ,lan) = qqlon(i,kt,lan)
!!              pten(i,k ,lan) = (qqlon(i,kp,lan)-pten(i,k,lan))*rdt2
!            pten(i,k ,lan)    =(qqlon(i,kp,lan)-(dsigma(KL,1)          &
!                               *ptp_sl(i,lan)+dsigma(KL,2)))*rdt2
          enddo
        enddo
! rq update
        do n=1,ncld
        do k=1,levs
          kk=levs-k+1+(n-1)*levs
          kq=kqq+k-1+(n-1)*levs
          do i=1,lons_lat
!ch         qvadv(i,kk,lan) = qqlon(i,kq,lan)
            qvadv(i,k,n,lan) = qqlon(i,kq,lan)
          enddo
        enddo
        enddo
!        if (myrank .eq. 0 ) then
!        call mymaxmin(qqlon(1,ktt,lan),lons_lat,lonfull,1,' tend temp')
!        endif
!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        ilan=1+jlonf
!        call mymaxmin(grid_gr(ilan,g_u),lons_lat,lonfull,1,' tend u ')
!        call mymaxmin(grid_gr(ilan,g_v),lons_lat,lonfull,1,' tend v ')
!        call mymaxmin(grid_gr(ilan,g_tt),lons_lat,lonfull,1,' tend t ')
!        call mymaxmin(grid_gr(ilan,g_rt),lons_lat,lonfull,1,' tend q ')
!        call mymaxmin(grid_gr(ilan,g_z),lons_lat,lonfull,1,'tend z')
!        do k=1,lev,10
!        kug=g_d+k-1
!        print *,' k=',k
!        call mymaxmin(grid_gr(ilan,kug),lons_lat,lonfull,1,' tend d ')
!        enddo
!        print *,' finish horizonatal advection at lan=',lan
!        endif
!
      enddo
!$omp end parallel do

!
! ===============================
!
      return
      end
!
!
      subroutine ndslfv_monoadvh2_xy(ddtemp,qvadv,pten,vdzonl,vdmerd  &
                , lonsperlat,deltim,levs) 
!
! a routine to do non-iteration semi-Lagrangain advection
! considering advection  with monotonicity in interpolation
! contact: hann-ming henry juang
! program log
! 2011 02 20 : henry juang, created for ndsl advection
! 2013 06 20 : Henry Juang correct wind direction for north-south advection
!
      use param
      use grid
      use index
      use const
      use rank

      implicit none

      real(kind=RTYPE) pten(nx,levs,my_max)
      real(kind=RTYPE) ddtemp(nx,levs,my_max),qvadv(nx,levs,ncld,my_max)
      real(kind=RTYPE) vdmerd(nx,levs,my_max),vdzonl(nx,levs,my_max)
!      real(kind=RTYPE) plev(lonfull,lev+1)
!      integer,intent(in):: global_lats_a(my)
      integer,intent(in):: lonsperlat(my)
      real(kind=RTYPE),   intent(in):: deltim

      real(kind=RTYPE)      uulon(lonfull,levs,latpart)
      real(kind=RTYPE)      vvlon(lonfull,levs,latpart)
      real(kind=RTYPE)      qqlon(lonfull,levs*ndslhvar,latpart)
!      real(kind=RTYPE)      rrlon(lonfull,lev*ndslhvar,latpart)

      real(kind=RTYPE)      vvlat(latfull,levs,lonpart)
      real(kind=RTYPE)      qqlat(latfull,levs*ndslhvar,lonpart)
!      real(kind=RTYPE)      rrlat(latfull,levs*ndslhvar,lonpart)
      real(kind=RTYPE)      xr    (lonfull,levs)
      real(kind=RTYPE)      xcp   (lonfull,levs)
      real(kind=RTYPE)      sumrq (lonfull,levs)
      real(kind=RTYPE)      xkappa(lonfull,levs)
      real(kind=RTYPE)      rdt2, rkt, pi, cons0, cons1, rma, rm2a

!      logical   lprint

      integer mono,mass,levs
      integer nlevs,nvars!,levh
      integer i,j,n,k,lon,lan,lat,lons_lat,irc,kk,KL
      integer kuu, kvv, kqq, ktt, kup, nqq
      integer ku , kv , kq , kt , kp
!
!      lprint = .false.

!      if( lprint ) print *,' enter ndslfv_monoadvh '
!
      mono  = 1
      mass  = 0
      cons0 = 0.0
      cons1 = 1.0
      qqlon = 0.
      uulon = 0.
      vvlon = 0.
!
!      levh = ncld * lev
!
!      kuu = 1
!      kvv = kuu + lev
!      kuu = kvv + lev
      kuu = 1
      kvv = kuu + levs
      ktt = kvv + levs
      kup = ktt + levs
!!      kqq = kup + levs
      kqq = ktt + levs

      nvars = ndslhvar
      nlevs = nvars * levs
!
      rdt2 = 0.5 / deltim
!
! =================================================================
!   prepare wind and variable in flux form with gaussina weight
! =================================================================
!

!$omp parallel do                                                      &
!$omp private(lan,lat,lons_lat,rma,rm2a,i,n,kk,k,kt,kv,ku,kp,kq) &
!$omp schedule(dynamic)
      do lan=1,jlistnum

!        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lat = jlist1(lan)
        lons_lat = lonsperlat(lat)
!        rma  = 1. / cosglat(lat) / con_rerth
!        rm2a = rma / cosglat(lat)
        rma  = 1. / cosl(lat) 
        rm2a = rma / cosl(lat)
!
! wind at time step n
        do k=1,levs
          kk=levs-k+1
          do i=1,lons_lat
!ch         uulon(i,k,lan) = ut(i,kk,lan) * rm2a
            uulon(i,k,lan) = ut_sl(i,k ,lan) * rm2a
!hmhj use real(kind=RTYPE) wind
!            vvlon(i,k,lan) = grid_gr(ilan,kvg) * rma
!hmhj use virtual wind
!            vvlon(i,k,lan) = grid_gr(ilan,kvg) / con_rerth
!for CWB GFS
!ch         vvlon(i,k,lan) = vt(i,kk,lan) * rma
            vvlon(i,k,lan) = vt_sl(i,k ,lan) * rma
          enddo
        enddo
!        if( lprint ) then
!          call mymaxmin(uulon(1,1,lan),lons_lat,lonfull,1,' uu1 in deg')
!          call mymaxmin(uulon(1,5,lan),lons_lat,lonfull,1,' uu5 in deg')
!          call mymaxmin(vvlon(1,1,lan),lons_lat,lonfull,1,' vv in deg')
!        endif
!
! u v t at n-1
!
        do k=1,levs
          ku=kuu+k-1
          kv=kvv+k-1
          kt=ktt+k-1
          kp=kup+k-1
          kk=levs-k+1
!ch>
!       KL=lev-Llist(k)+1
!ch<
          do i=1,lons_lat
!ch         qqlon(i,ku,lan) = uum(i,kk,lan)
!ch         qqlon(i,kv,lan) = vvm(i,kk,lan)
!ch         qqlon(i,kt,lan) = ttm(i,kk,lan)
!!          qqlon(i,kp,lan) = dsigma(kk,1)*ptp(i,lan)+dsigma(kk,2)

            qqlon(i,ku,lan) = vdzonl(i,k,lan)
            qqlon(i,kv,lan) = vdmerd(i,k,lan)
            qqlon(i,kt,lan) = ddtemp(i,k,lan)
!!            qqlon(i,kp,lan) = pten(i,k,lan)
!!            qqlon(i,kp,lan) = dsigma(KL,1)*ptp_sl(i,lan)+dsigma(KL,2)
          enddo
        enddo
! rq at n-1
        do n=1,ncld
        do k=1,levs
          kq=kqq+k-1+(n-1)*levs
          kk=levs-k+1+(n-1)*levs
          do i=1,lons_lat
!ch         qqlon(i,kq,lan) = qm(i,kk,lan)
            qqlon(i,kq,lan) = qvadv(i,k,n,lan)
          enddo
        enddo
        enddo
! add surface pressure perturbation for removing resonance
!       rkt = con_g / ( con_rd * 300.0 )
!       do i=1,lons_lat
!         ilan=i+jlonf
!         qqlon(i,kuu,lan) = log(plev(i,1))+grid_gr(ilan,g_gz)*rkt
!       enddo
!
! save qqlon into n+1 for later as tendency
! z d t at n+1
!        do k=1,lev
!          kk=lev-k+1
!          ku=kuu+k-1
!          kv=kvv+k-1
!          kt=ktt+k-1
!          kug=g_d+k-1
!          kvg=g_z+k-1
!          ktg=g_t+k-1
!          do i=1,lons_lat
!            ilan=i+jlonf
!            grid_gr(ilan,kug) = qqlon(i,ku,lan)
!            grid_gr(ilan,kvg) = qqlon(i,kv,lan)
!            grid_gr(ilan,ktg) = qqlon(i,kt,lan)
!          enddo
!        enddo
! change h to theta
!
! ----- prepare xr, xcp, xkapa
!
!       xr    = cons0
!       xcp   = cons0
!       sumrq = cons0
!
!       do n=1,ntrac
!         nqq = kqq + (n-1)*levs
!         if( ri(n) .ne. cons0 .and. cpi(n) .ne. cons0 ) then
!           do k=1,lev
!             kq=nqq+k-1
!             do i=1,lons_lat
!               xr   (i,k) = xr   (i,k) + qqlon(i,kq,lan)*ri(n)
!               xcp  (i,k) = xcp  (i,k) + qqlon(i,kq,lan)*cpi(n)
!               sumrq(i,k) = sumrq(i,k) + qqlon(i,kq,lan)
!             enddo
!           enddo
!         endif
!       enddo
!       do k=1,lev
!         do i=1,lons_lat
!           xr (i,k)   = ( cons1 - sumrq(i,k) )*ri(0)  + xr (i,k)
!           xcp(i,k)   = ( cons1 - sumrq(i,k) )*cpi(0) + xcp(i,k)
!           xkappa(i,k) = xr(i,k) / xcp(i,k)
!         enddo
!       enddo
!
!       do k=1,lev
!         kt=ktt+k-1
!         ku=kuu+k-1
!         do i=1,lons_lat
!           pi = (qqlon(i,ku,lan)*0.005)**xkappa(i,k)
!           qqlon(i,kt,lan) = qqlon(i,kt,lan)/pi
!         enddo
!       enddo
!
! rq  no need for tendency
!
!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        ilan=1+jlonf
!        call mymaxmin(grid_gr(ilan,g_ttm),lons_lat,lonfull,1,' n-1 t ')
!        call mymaxmin(grid_gr(ilan,g_rrm ),lons_lat,lonfull,1,' n-1 q ')
!        call mymaxmin(grid_gr(ilan,g_zem),lons_lat,lonfull,1,' n-1 z ')
!        call mymaxmin(grid_gr(ilan,g_dim),lons_lat,lonfull,1,' n-1 d ')
!        print *,' ------------------------------------------- '
!        call mymaxmin(qqlon(1,ktt,lan),lons_lat,lonfull,1,' red t ')
!        call mymaxmin(qqlon(1,kqq,lan),lons_lat,lonfull,1,' red q ')
!        call mymaxmin(qqlon(1,kvv,lan),lons_lat,lonfull,1,' red z ')
!        call mymaxmin(qqlon(1,kuu,lan),lons_lat,lonfull,1,' red d ')
!        endif

!       call cyclic_cell_intpx(levs,lons_lat,lonf,uulon(1,1,lan))

!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        call mymaxmin(qqlon(1,ktt,lan),lonfull,lonfull,1,' full t ')
!        call mymaxmin(qqlon(1,kqq,lan),lonfull,lonfull,1,' full q ')
!        call mymaxmin(qqlon(1,kvv,lan),lonfull,lonfull,1,' full z ')
!        call mymaxmin(qqlon(1,kuu,lan),lonfull,lonfull,1,' full d ')
!        endif

!!        rrlon(:,:,lan) = qqlon(:,:,lan)
!
! first set positive advection in east-west direction
!
        call cyclic_cell_massadvx(lons_lat,lonfull,levs,nvars,deltim,   &
                         uulon(1,1,lan),qqlon(1,1,lan),mass)
!       call cyclic_mono_advectx (lonfull,levs,nvars,deltim,               &
!    &                   uulon(1,1,lan),rrlon(1,1,lan),mono)

!!        call cyclic_cell_intpx(nlevs,lons_lat,lonfull,rrlon(1,1,lan))
        call cyclic_cell_intpx(nlevs,lons_lat,lonfull,qqlon(1,1,lan))
        call cyclic_cell_intpx(levs,lons_lat,lonfull,vvlon(1,1,lan))

!        if( lprint ) then
!        print *,' done cyclic_massadvx with mass= ',mass
!        print *,' ------------------------------------------- '
!        call mymaxmin(rrlon(1,ktt,lan),lonfull,lonfull,1,' advx t ')
!        call mymaxmin(rrlon(1,kqq,lan),lonfull,lonfull,1,' advx q ')
!        call mymaxmin(rrlon(1,kvv,lan),lonfull,lonfull,1,' advx z ')
!        call mymaxmin(rrlon(1,kuu,lan),lonfull,lonfull,1,' advx d ')
!        print *,' done the first x adv for lan=',lan
!        print *,' =========================================== '
!        endif

      enddo
!$omp end parallel do

! ---------------------------------------------------------------------
! mpi para from east-west full grid to north-south full grid
! ---------------------------------------------------------------------
!
! para vvlon, qqlon, and rrlon to vvlat, qqlat, rrlat
!       if( lprint ) print *,' ndslfv_advect transport from we to ns '

       call para_we2ns(vvlon,vvlat,levs,my)
!!       call para_we2ns(rrlon,rrlat,nlevs,my)
       call para_we2ns(qqlon,qqlat,nlevs,my)

!       if( lprint ) then
!       print *,' ------------ after we2ns ---------------------- '
!       do lon=1,mylonlen
!        print *,'  lon=',lon
!        call mymaxmin(vvlat(1,1  ,lon),latfull,latfull,1,' we2ns v')
!        call mymaxmin(rrlat(1,kqq,lon),latfull,latfull,1,' we2ns r')
!        call mymaxmin(qqlat(1,kqq,lon),latfull,latfull,1,' we2ns q')
!       enddo
!       endif
!
! ---------------------------------------------------------------------
! -------------- in north-soutn great circle -------------------
! ---------------------------------------------------------------------

!      if ( myrank .eq. 0 ) lprint = .true.
!       if( lprint ) then
!       print *,' ndslfv_advect adv loop in y '
!       print *,' mylonlen=',mylonlen
!       endif

!$omp parallel do private(lon,k,j,ku,kv) &
!$omp schedule(dynamic)
       do lon=1,mylonlen
!
!        if( lprint ) print *,' lon=',lon
! convert wind before advy
        do k=1,levs
          ku=kuu+k-1
          kv=kvv+k-1
          do j=1,lathalf
            qqlat(j,ku,lon) = -qqlat(j,ku,lon)
            qqlat(j,kv,lon) = -qqlat(j,kv,lon)
          enddo
        enddo
! first set advection in north-south direction in great circle through two poles
!
!!        call fixend_cell_massadvy(latfull,lathalf,lev,nvars,deltim, &
!!                         vvlat(1,1,lon),rrlat(1,1,lon),mass)
!       call cyclic_cell_massadvy(latfull,levs,nvars,deltim,
!    &                   vvlat(1,1,lon),rrlat(1,1,lon),mass)
!       call cyclic_mono_advecty (latfull,levs,nvars,deltim,
!    &                   vvlat(1,1,lon),rrlat(1,1,lon),mono)

!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        call mymaxmin(rrlat(1,kuu,lon),latfull,latfull,1,' advry u ')
!        call mymaxmin(rrlat(1,kvv,lon),latfull,latfull,1,' advry v ')
!        call mymaxmin(rrlat(1,ktt,lon),latfull,latfull,1,' advry t ')
!        call mymaxmin(rrlat(1,kqq,lon),latfull,latfull,1,' advry q ')
!        call mymaxmin(rrlat(1,kvv,lon),latfull,latfull,1,' advry z ')
!        call mymaxmin(rrlat(1,kuu,lon),latfull,latfull,1,' advry d ')
!        endif
!
! second set advection in north-south direction in great circle through two poles
!
!        call fixend_cell_massadvy(latfull,lathalf,levs,nvars,deltim, &
!                         vvlat(1,1,lon),qqlat(1,1,lon),mass)
        call cyclic_cell_massadvy(latfull,levs,nvars,deltim,          &
                        vvlat(1,1,lon),qqlat(1,1,lon),mass)
!       call cyclic_mono_advecty (my,levs,nvars,deltim,
!    &                   vvlat(1,1,lon),qqlat(1,1,lon),mono)

!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        call mymaxmin(qqlat(1,kuu,lon),latfull,latfull,1,' advqy u ')
!        call mymaxmin(qqlat(1,kvv,lon),latfull,latfull,1,' advqy v ')
!        call mymaxmin(qqlat(1,ktt,lon),latfull,latfull,1,' advqy t ')
!        call mymaxmin(qqlat(1,kqq,lon),latfull,latfull,1,' advqy q ')
!        call mymaxmin(qqlat(1,kvv,lon),latfull,latfull,1,' advqy z ')
!        call mymaxmin(qqlat(1,kuu,lon),latfull,latfull,1,' advqy d ')
!        print *,' done with y at lon=',lon
!        endif

! convert wind back after advy
        do k=1,levs
          ku=kuu+k-1
          kv=kvv+k-1
          do j=1,lathalf
            qqlat(j,ku,lon) = -qqlat(j,ku,lon)
            qqlat(j,kv,lon) = -qqlat(j,kv,lon)
          enddo
        enddo


       enddo
!$omp end parallel do
!
! ----------------------------------------------------------------------
! mpi para from north-south direction to east-west direeectory
! ----------------------------------------------------------------------
!
! para qqlat and rrlat to qqlon and rrlon

!       if( lprint ) print *,' ndslfv_advect transport from ns to we '

!!       call para_ns2we(rrlat,rrlon,nlevs,my)
       call para_ns2we(qqlat,qqlon,nlevs,my)
!       if( lprint ) then
!       print *,' ------------ after ns2we ---------------------- '
!       do lan=1,jlistnum
!        print *,'  lan=',lan
!        call mymaxmin(rrlon(1,kqq,lan),lonfull,lonfull,1,' ns2we r')
!        call mymaxmin(qqlon(1,kqq,lan),lonfull,lonfull,1,' ns2we q')
!       enddo
!       endif

! ---------------------------------------------------------------
! ---------------- back to east-west direction ------------------
! ---------------------------------------------------------------
!      print *,' ndslfv_advect adv loop in x for last '

!$omp parallel do                                                &
!$omp private(lan,lat,lons_lat,i,k,kk,n,kt,kv,ku,kp,kq)    &
!$omp schedule(dynamic)

      do lan=1,jlistnum

!        lat = global_lats_a(ipt_jlistnum-1+lan)
        lat = jlist1(lan)
        lons_lat = lonsperlat(lat)
!
! mass conserving interpolation from full grid to reduced grid
!
        call cyclic_cell_intpx(nlevs,lonfull,lons_lat,qqlon(1,1,lan))
!!        call cyclic_cell_intpx(nlevs,lonfull,lons_lat,rrlon(1,1,lan))

!
! second set advection in x for the second of the pair
!
!!        call cyclic_cell_massadvx(lons_lat,lonfull,levs,nvars,deltim,   &
!!                         uulon(1,1,lan),qqlon(1,1,lan),mass)
!       call cyclic_mono_advectx (lonfull,levs,nvars,deltim,               &
!    &                   uulon(1,1,lan),qqlon(1,1,lan),mono)
!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        call mymaxmin(qqlon(1,kqq,lan),lonfull,lonfull,1,' adv x q ')
!        endif

!!        do k=1,nlevs
!!          do i=1,lons_lat
!!            qqlon(i,k,lan) = 0.5 * ( qqlon(i,k,lan) + rrlon(i,k,lan) )
!!          enddo
!!        enddo

!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        call mymaxmin(qqlon(1,kuu,lan),lons_lat,lonfull,1,' do redu u ')
!        call mymaxmin(qqlon(1,kvv,lan),lons_lat,lonfull,1,' do redu v ')
!        call mymaxmin(qqlon(1,ktt,lan),lons_lat,lonfull,1,' do redu t ')
!        call mymaxmin(qqlon(1,kqq,lan),lons_lat,lonfull,1,' do redu q ')
!        call mymaxmin(qqlon(1,kvv,lan),lons_lat,lonfull,1,' do redu z')
!        call mymaxmin(qqlon(1,kuu,lan),lons_lat,lonfull,1,' do redu d')
!        print *,' finish horizonatal advection at lan=',lan
!        endif

! change theta to h
!
! ----- prepare xr, xcp, xkapa
!
!       xr    = cons0
!       xcp   = cons0
!       sumrq = cons0
!
!       do n=1,ntrac
!         nqq = kqq + (n-1)*levs
!         if( ri(n) .ne. cons0 .and. cpi(n) .ne. cons0 ) then
!           do k=1,lev
!             kq=nqq+k-1
!             do i=1,lons_lat
!               xr   (i,k) = xr   (i,k) + qqlon(i,kq,lan)*ri(n)
!               xcp  (i,k) = xcp  (i,k) + qqlon(i,kq,lan)*cpi(n)
!               sumrq(i,k) = sumrq(i,k) + qqlon(i,kq,lan)
!             enddo
!           enddo
!         endif
!       enddo
!       do k=1,lev
!         do i=1,lons_lat
!           xr (i,k)   = ( cons1 - sumrq(i,k) )*ri(0)  + xr (i,k)
!           xcp(i,k)   = ( cons1 - sumrq(i,k) )*cpi(0) + xcp(i,k)
!           xkappa(i,k) = xr(i,k) / xcp(i,k)
!         enddo
!       enddo
!
!       do k=1,lev
!         kt=ktt+k-1
!         ku=kuu+k-1
!         do i=1,lons_lat
!           pi = (qqlon(i,ku,lan)*0.005)**xkappa(i,k)
!           qqlon(i,kt,lan) = qqlon(i,kt,lan)*pi
!         enddo
!       enddo
! u v t update at n+1
        do k=1,levs
          kk=levs-k+1
          ku=kuu+k-1
          kv=kvv+k-1
          kt=ktt+k-1
          kp=kup+k-1
!ch>
!        KL=lev-Llist(k)+1
!ch<
          do i=1,lons_lat
!ch         vdzonl(i,kk,lan) = (qqlon(i,ku,lan)-uum(i,kk,lan))*rdt2
!ch         vdmerd(i,kk,lan) = (qqlon(i,kv,lan)-vvm(i,kk,lan))*rdt2
!ch         ddtemp(i,kk,lan) = qqlon(i,kt,lan)
!!          pten(i,kk,lan)    =(qqlon(i,kp,lan)-(dsigma(kk,1)          &
!!                             *ptp(i,lan)+dsigma(kk,2)))*rdt2

            vdzonl(i,k ,lan) = qqlon(i,ku,lan)
            vdmerd(i,k ,lan) = qqlon(i,kv,lan)
            ddtemp(i,k ,lan) = qqlon(i,kt,lan)
!!              pten(i,k ,lan) = (qqlon(i,kp,lan)-pten(i,k,lan))*rdt2
!!            pten(i,k ,lan)    =(qqlon(i,kp,lan)-(dsigma(KL,1)          &
!!                               *ptp_sl(i,lan)+dsigma(KL,2)))*rdt2
          enddo
        enddo
! rq update
        do n=1,ncld
        do k=1,levs
          kk=levs-k+1+(n-1)*levs
          kq=kqq+k-1+(n-1)*levs
          do i=1,lons_lat
!ch         qvadv(i,kk,lan) = qqlon(i,kq,lan)
            qvadv(i,k,n,lan) = qqlon(i,kq,lan)
          enddo
        enddo
        enddo
!        if (myrank .eq. 0 ) then
!        call mymaxmin(qqlon(1,ktt,lan),lons_lat,lonfull,1,' tend temp')
!        endif
!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        ilan=1+jlonf
!        call mymaxmin(grid_gr(ilan,g_u),lons_lat,lonfull,1,' tend u ')
!        call mymaxmin(grid_gr(ilan,g_v),lons_lat,lonfull,1,' tend v ')
!        call mymaxmin(grid_gr(ilan,g_tt),lons_lat,lonfull,1,' tend t ')
!        call mymaxmin(grid_gr(ilan,g_rt),lons_lat,lonfull,1,' tend q ')
!        call mymaxmin(grid_gr(ilan,g_z),lons_lat,lonfull,1,'tend z')
!        do k=1,lev,10
!        kug=g_d+k-1
!        print *,' k=',k
!        call mymaxmin(grid_gr(ilan,kug),lons_lat,lonfull,1,' tend d ')
!        enddo
!        print *,' finish horizonatal advection at lan=',lan
!        endif
!
      enddo
!$omp end parallel do

!
! ===============================
!
      return
      end
!
!
      subroutine ndslfv_monoadvh2_yx(ddtemp,qvadv,pten,vdzonl,vdmerd  &
                , lonsperlat,deltim,levs)
!
! a routine to do non-iteration semi-Lagrangain advection
! considering advection  with monotonicity in interpolation
! contact: hann-ming henry juang
! program log
! 2011 02 20 : henry juang, created for ndsl advection
! 2013 06 20 : Henry Juang correct wind direction for north-south advection
!
      use param
      use grid
      use index
      use const
      use rank

      implicit none

      real(kind=RTYPE) pten(nx,levs,my_max)
      real(kind=RTYPE) ddtemp(nx,levs,my_max),qvadv(nx,levs,ncld,my_max)
      real(kind=RTYPE) vdmerd(nx,levs,my_max),vdzonl(nx,levs,my_max)
!      real(kind=RTYPE) plev(lonfull,levs+1)
!      integer,intent(in):: global_lats_a(my)
      integer,intent(in):: lonsperlat(my)
      real(kind=RTYPE),   intent(in):: deltim

      real(kind=RTYPE)      uulon(lonfull,levs,latpart)
      real(kind=RTYPE)      vvlon(lonfull,levs,latpart)
      real(kind=RTYPE)      qqlon(lonfull,levs*ndslhvar,latpart)
!      real(kind=RTYPE)      rrlon(lonfull,levs*ndslhvar,latpart)

      real(kind=RTYPE)      vvlat(latfull,levs,lonpart)
      real(kind=RTYPE)      qqlat(latfull,levs*ndslhvar,lonpart)
!      real(kind=RTYPE)      rrlat(latfull,levs*ndslhvar,lonpart)
      real(kind=RTYPE)      xr    (lonfull,levs)
      real(kind=RTYPE)      xcp   (lonfull,levs)
      real(kind=RTYPE)      sumrq (lonfull,levs)
      real(kind=RTYPE)      xkappa(lonfull,levs)
      real(kind=RTYPE)      rdt2, rkt, pi, cons0, cons1, rma, rm2a

!      logical   lprint

      integer mono,mass,levs
      integer nlevs,nvars!,levh
      integer i,j,n,k,lon,lan,lat,lons_lat,irc,kk,KL
      integer kuu, kvv, kqq, ktt, kup, nqq
      integer ku , kv , kq , kt , kp
!
!      lprint = .false.

!      if( lprint ) print *,' enter ndslfv_monoadvh '
!
      mono  = 1
      mass  = 0
      cons0 = 0.0
      cons1 = 1.0
      qqlon = 0.
      uulon = 0.
      vvlon = 0.
!
!      levh = ncld * levs
!
!      kuu = 1
!      kvv = kuu + levs
!      kuu = kvv + levs
      kuu = 1
      kvv = kuu + levs
      ktt = kvv + levs
      kup = ktt + levs
!!      kqq = kup + levs
      kqq = ktt + levs

      nvars = ndslhvar
      nlevs = nvars * levs
!
      rdt2 = 0.5 / deltim
!
! =================================================================
!   prepare wind and variable in flux form with gaussina weight
! =================================================================
!

!$omp parallel do                                                      &
!$omp private(lan,lat,lons_lat,rma,rm2a,i,n,kk,k,kt,kv,ku,kp,kq) &
!$omp schedule(dynamic)
      do lan=1,jlistnum

!        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lat = jlist1(lan)
        lons_lat = lonsperlat(lat)
!        rma  = 1. / cosglat(lat) / con_rerth
!        rm2a = rma / cosglat(lat)
        rma  = 1. / cosl(lat) 
        rm2a = rma / cosl(lat)
!
! wind at time step n
        do k=1,levs
          kk=levs-k+1
          do i=1,lons_lat
!ch         uulon(i,k,lan) = ut(i,kk,lan) * rm2a
            uulon(i,k,lan) = ut_sl(i,k ,lan) * rm2a
!hmhj use real(kind=RTYPE) wind
!            vvlon(i,k,lan) = grid_gr(ilan,kvg) * rma
!hmhj use virtual wind
!            vvlon(i,k,lan) = grid_gr(ilan,kvg) / con_rerth
!for CWB GFS
!ch         vvlon(i,k,lan) = vt(i,kk,lan) * rma
            vvlon(i,k,lan) = vt_sl(i,k ,lan) * rma
          enddo
        enddo
!        if( lprint ) then
!          call mymaxmin(uulon(1,1,lan),lons_lat,lonfull,1,' uu1 in deg')
!          call mymaxmin(uulon(1,5,lan),lons_lat,lonfull,1,' uu5 in deg')
!          call mymaxmin(vvlon(1,1,lan),lons_lat,lonfull,1,' vv in deg')
!        endif
!
! u v t at n-1
!
        do k=1,levs
          ku=kuu+k-1
          kv=kvv+k-1
          kt=ktt+k-1
          kp=kup+k-1
          kk=levs-k+1
!ch>
!       KL=lev-Llist(k)+1
!ch<
          do i=1,lons_lat
!ch         qqlon(i,ku,lan) = uum(i,kk,lan)
!ch         qqlon(i,kv,lan) = vvm(i,kk,lan)
!ch         qqlon(i,kt,lan) = ttm(i,kk,lan)
!!          qqlon(i,kp,lan) = dsigma(kk,1)*ptp(i,lan)+dsigma(kk,2)

            qqlon(i,ku,lan) = vdzonl(i,k,lan)
            qqlon(i,kv,lan) = vdmerd(i,k,lan)
            qqlon(i,kt,lan) = ddtemp(i,k,lan)
!!            qqlon(i,kp,lan) = pten(i,k,lan)
!!            qqlon(i,kp,lan) = dsigma(KL,1)*ptp_sl(i,lan)+dsigma(KL,2)
          enddo
        enddo
! rq at n-1
        do n=1,ncld
        do k=1,levs
          kq=kqq+k-1+(n-1)*levs
          kk=levs-k+1+(n-1)*levs
          do i=1,lons_lat
!ch         qqlon(i,kq,lan) = qm(i,kk,lan)
            qqlon(i,kq,lan) = qvadv(i,k,n,lan)
          enddo
        enddo
        enddo
! add surface pressure perturbation for removing resonance
!       rkt = con_g / ( con_rd * 300.0 )
!       do i=1,lons_lat
!         ilan=i+jlonf
!         qqlon(i,kuu,lan) = log(plev(i,1))+grid_gr(ilan,g_gz)*rkt
!       enddo
!
! save qqlon into n+1 for later as tendency
! z d t at n+1
!        do k=1,levs
!          kk=levs-k+1
!          ku=kuu+k-1
!          kv=kvv+k-1
!          kt=ktt+k-1
!          kug=g_d+k-1
!          kvg=g_z+k-1
!          ktg=g_t+k-1
!          do i=1,lons_lat
!            ilan=i+jlonf
!            grid_gr(ilan,kug) = qqlon(i,ku,lan)
!            grid_gr(ilan,kvg) = qqlon(i,kv,lan)
!            grid_gr(ilan,ktg) = qqlon(i,kt,lan)
!          enddo
!        enddo
! change h to theta
!
! ----- prepare xr, xcp, xkapa
!
!       xr    = cons0
!       xcp   = cons0
!       sumrq = cons0
!
!       do n=1,ntrac
!         nqq = kqq + (n-1)*levs
!         if( ri(n) .ne. cons0 .and. cpi(n) .ne. cons0 ) then
!           do k=1,levs
!             kq=nqq+k-1
!             do i=1,lons_lat
!               xr   (i,k) = xr   (i,k) + qqlon(i,kq,lan)*ri(n)
!               xcp  (i,k) = xcp  (i,k) + qqlon(i,kq,lan)*cpi(n)
!               sumrq(i,k) = sumrq(i,k) + qqlon(i,kq,lan)
!             enddo
!           enddo
!         endif
!       enddo
!       do k=1,levs
!         do i=1,lons_lat
!           xr (i,k)   = ( cons1 - sumrq(i,k) )*ri(0)  + xr (i,k)
!           xcp(i,k)   = ( cons1 - sumrq(i,k) )*cpi(0) + xcp(i,k)
!           xkappa(i,k) = xr(i,k) / xcp(i,k)
!         enddo
!       enddo
!
!       do k=1,levs
!         kt=ktt+k-1
!         ku=kuu+k-1
!         do i=1,lons_lat
!           pi = (qqlon(i,ku,lan)*0.005)**xkappa(i,k)
!           qqlon(i,kt,lan) = qqlon(i,kt,lan)/pi
!         enddo
!       enddo
!
! rq  no need for tendency
!
!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        ilan=1+jlonf
!        call mymaxmin(grid_gr(ilan,g_ttm),lons_lat,lonfull,1,' n-1 t ')
!        call mymaxmin(grid_gr(ilan,g_rrm ),lons_lat,lonfull,1,' n-1 q ')
!        call mymaxmin(grid_gr(ilan,g_zem),lons_lat,lonfull,1,' n-1 z ')
!        call mymaxmin(grid_gr(ilan,g_dim),lons_lat,lonfull,1,' n-1 d ')
!        print *,' ------------------------------------------- '
!        call mymaxmin(qqlon(1,ktt,lan),lons_lat,lonfull,1,' red t ')
!        call mymaxmin(qqlon(1,kqq,lan),lons_lat,lonfull,1,' red q ')
!        call mymaxmin(qqlon(1,kvv,lan),lons_lat,lonfull,1,' red z ')
!        call mymaxmin(qqlon(1,kuu,lan),lons_lat,lonfull,1,' red d ')
!        endif

!       call cyclic_cell_intpx(levs,lons_lat,lonf,uulon(1,1,lan))

!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        call mymaxmin(qqlon(1,ktt,lan),lonfull,lonfull,1,' full t ')
!        call mymaxmin(qqlon(1,kqq,lan),lonfull,lonfull,1,' full q ')
!        call mymaxmin(qqlon(1,kvv,lan),lonfull,lonfull,1,' full z ')
!        call mymaxmin(qqlon(1,kuu,lan),lonfull,lonfull,1,' full d ')
!        endif

!!        rrlon(:,:,lan) = qqlon(:,:,lan)
!
! first set positive advection in east-west direction
!
!!        call cyclic_cell_massadvx(lons_lat,lonfull,lev,nvars,deltim,   &
!!                         uulon(1,1,lan),qqlon(1,1,lan),mass)
!       call cyclic_mono_advectx (lonfull,lev,nvars,deltim,               &
!    &                   uulon(1,1,lan),rrlon(1,1,lan),mono)

!!        call cyclic_cell_intpx(nlevs,lons_lat,lonfull,rrlon(1,1,lan))
        call cyclic_cell_intpx(nlevs,lons_lat,lonfull,qqlon(1,1,lan))
        call cyclic_cell_intpx(levs,lons_lat,lonfull,vvlon(1,1,lan))

!        if( lprint ) then
!        print *,' done cyclic_massadvx with mass= ',mass
!        print *,' ------------------------------------------- '
!        call mymaxmin(rrlon(1,ktt,lan),lonfull,lonfull,1,' advx t ')
!        call mymaxmin(rrlon(1,kqq,lan),lonfull,lonfull,1,' advx q ')
!        call mymaxmin(rrlon(1,kvv,lan),lonfull,lonfull,1,' advx z ')
!        call mymaxmin(rrlon(1,kuu,lan),lonfull,lonfull,1,' advx d ')
!        print *,' done the first x adv for lan=',lan
!        print *,' =========================================== '
!        endif

      enddo
!$omp end parallel do

! ---------------------------------------------------------------------
! mpi para from east-west full grid to north-south full grid
! ---------------------------------------------------------------------
!
! para vvlon, qqlon, and rrlon to vvlat, qqlat, rrlat
!       if( lprint ) print *,' ndslfv_advect transport from we to ns '

       call para_we2ns(vvlon,vvlat,levs,my)
!!       call para_we2ns(rrlon,rrlat,nlevs,my)
       call para_we2ns(qqlon,qqlat,nlevs,my)

!       if( lprint ) then
!       print *,' ------------ after we2ns ---------------------- '
!       do lon=1,mylonlen
!        print *,'  lon=',lon
!        call mymaxmin(vvlat(1,1  ,lon),latfull,latfull,1,' we2ns v')
!        call mymaxmin(rrlat(1,kqq,lon),latfull,latfull,1,' we2ns r')
!        call mymaxmin(qqlat(1,kqq,lon),latfull,latfull,1,' we2ns q')
!       enddo
!       endif
!
! ---------------------------------------------------------------------
! -------------- in north-soutn great circle -------------------
! ---------------------------------------------------------------------

!      if ( myrank .eq. 0 ) lprint = .true.
!       if( lprint ) then
!       print *,' ndslfv_advect adv loop in y '
!       print *,' mylonlen=',mylonlen
!       endif

!$omp parallel do private(lon,k,j,ku,kv) &
!$omp schedule(dynamic)
       do lon=1,mylonlen
!
!        if( lprint ) print *,' lon=',lon
! convert wind before advy
        do k=1,levs
          ku=kuu+k-1
          kv=kvv+k-1
          do j=1,lathalf
            qqlat(j,ku,lon) = -qqlat(j,ku,lon)
            qqlat(j,kv,lon) = -qqlat(j,kv,lon)
          enddo
        enddo

! first set advection in north-south direction in great circle through two poles
!
!!        call fixend_cell_massadvy(latfull,lathalf,lev,nvars,deltim, &
!!                         vvlat(1,1,lon),rrlat(1,1,lon),mass)
!       call cyclic_cell_massadvy(latfull,lev,nvars,deltim,
!    &                   vvlat(1,1,lon),rrlat(1,1,lon),mass)
!       call cyclic_mono_advecty (latfull,lev,nvars,deltim,
!    &                   vvlat(1,1,lon),rrlat(1,1,lon),mono)

!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        call mymaxmin(rrlat(1,kuu,lon),latfull,latfull,1,' advry u ')
!        call mymaxmin(rrlat(1,kvv,lon),latfull,latfull,1,' advry v ')
!        call mymaxmin(rrlat(1,ktt,lon),latfull,latfull,1,' advry t ')
!        call mymaxmin(rrlat(1,kqq,lon),latfull,latfull,1,' advry q ')
!        call mymaxmin(rrlat(1,kvv,lon),latfull,latfull,1,' advry z ')
!        call mymaxmin(rrlat(1,kuu,lon),latfull,latfull,1,' advry d ')
!        endif
!
! second set advection in north-south direction in great circle through two poles
!
!        call fixend_cell_massadvy(latfull,lathalf,levs,nvars,deltim, &
!                         vvlat(1,1,lon),qqlat(1,1,lon),mass)
        call cyclic_cell_massadvy(latfull,levs,nvars,deltim,            &
                         vvlat(1,1,lon),qqlat(1,1,lon),mass)
!       call cyclic_mono_advecty (my,lev,nvars,deltim,
!    &                   vvlat(1,1,lon),qqlat(1,1,lon),mono)

!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        call mymaxmin(qqlat(1,kuu,lon),latfull,latfull,1,' advqy u ')
!        call mymaxmin(qqlat(1,kvv,lon),latfull,latfull,1,' advqy v ')
!        call mymaxmin(qqlat(1,ktt,lon),latfull,latfull,1,' advqy t ')
!        call mymaxmin(qqlat(1,kqq,lon),latfull,latfull,1,' advqy q ')
!        call mymaxmin(qqlat(1,kvv,lon),latfull,latfull,1,' advqy z ')
!        call mymaxmin(qqlat(1,kuu,lon),latfull,latfull,1,' advqy d ')
!        print *,' done with y at lon=',lon
!        endif
! convert wind back after advy
        do k=1,levs
          ku=kuu+k-1
          kv=kvv+k-1
          do j=1,lathalf
            qqlat(j,ku,lon) = -qqlat(j,ku,lon)
            qqlat(j,kv,lon) = -qqlat(j,kv,lon)
          enddo
        enddo


       enddo
!$omp end parallel do
!
! ----------------------------------------------------------------------
! mpi para from north-south direction to east-west direeectory
! ----------------------------------------------------------------------
!
! para qqlat and rrlat to qqlon and rrlon

!       if( lprint ) print *,' ndslfv_advect transport from ns to we '

!!       call para_ns2we(rrlat,rrlon,nlevs,my)
       call para_ns2we(qqlat,qqlon,nlevs,my)
!       if( lprint ) then
!       print *,' ------------ after ns2we ---------------------- '
!       do lan=1,jlistnum
!        print *,'  lan=',lan
!        call mymaxmin(rrlon(1,kqq,lan),lonfull,lonfull,1,' ns2we r')
!        call mymaxmin(qqlon(1,kqq,lan),lonfull,lonfull,1,' ns2we q')
!       enddo
!       endif

! ---------------------------------------------------------------
! ---------------- back to east-west direction ------------------
! ---------------------------------------------------------------
!      print *,' ndslfv_advect adv loop in x for last '

!$omp parallel do                                                &
!$omp private(lan,lat,lons_lat,i,k,kk,n,kt,kv,ku,kq)    &
!$omp schedule(dynamic)

      do lan=1,jlistnum

!        lat = global_lats_a(ipt_jlistnum-1+lan)
        lat = jlist1(lan)
        lons_lat = lonsperlat(lat)
!
! mass conserving interpolation from full grid to reduced grid
!
        call cyclic_cell_intpx(nlevs,lonfull,lons_lat,qqlon(1,1,lan))
!!        call cyclic_cell_intpx(nlevs,lonfull,lons_lat,rrlon(1,1,lan))

!
! second set advection in x for the second of the pair
!
        call cyclic_cell_massadvx(lons_lat,lonfull,levs,nvars,deltim,   &
                         uulon(1,1,lan),qqlon(1,1,lan),mass)
!       call cyclic_mono_advectx (lonfull,lev,nvars,deltim,               &
!    &                   uulon(1,1,lan),qqlon(1,1,lan),mono)
!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        call mymaxmin(qqlon(1,kqq,lan),lonfull,lonfull,1,' adv x q ')
!        endif

!!        do k=1,nlevs
!!          do i=1,lons_lat
!!            qqlon(i,k,lan) = 0.5 * ( qqlon(i,k,lan) + rrlon(i,k,lan) )
!!          enddo
!!        enddo

!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        call mymaxmin(qqlon(1,kuu,lan),lons_lat,lonfull,1,' do redu u ')
!        call mymaxmin(qqlon(1,kvv,lan),lons_lat,lonfull,1,' do redu v ')
!        call mymaxmin(qqlon(1,ktt,lan),lons_lat,lonfull,1,' do redu t ')
!        call mymaxmin(qqlon(1,kqq,lan),lons_lat,lonfull,1,' do redu q ')
!        call mymaxmin(qqlon(1,kvv,lan),lons_lat,lonfull,1,' do redu z')
!        call mymaxmin(qqlon(1,kuu,lan),lons_lat,lonfull,1,' do redu d')
!        print *,' finish horizonatal advection at lan=',lan
!        endif

! change theta to h
!
! ----- prepare xr, xcp, xkapa
!
!       xr    = cons0
!       xcp   = cons0
!       sumrq = cons0
!
!       do n=1,ntrac
!         nqq = kqq + (n-1)*levs
!         if( ri(n) .ne. cons0 .and. cpi(n) .ne. cons0 ) then
!           do k=1,lev
!             kq=nqq+k-1
!             do i=1,lons_lat
!               xr   (i,k) = xr   (i,k) + qqlon(i,kq,lan)*ri(n)
!               xcp  (i,k) = xcp  (i,k) + qqlon(i,kq,lan)*cpi(n)
!               sumrq(i,k) = sumrq(i,k) + qqlon(i,kq,lan)
!             enddo
!           enddo
!         endif
!       enddo
!       do k=1,lev
!         do i=1,lons_lat
!           xr (i,k)   = ( cons1 - sumrq(i,k) )*ri(0)  + xr (i,k)
!           xcp(i,k)   = ( cons1 - sumrq(i,k) )*cpi(0) + xcp(i,k)
!           xkappa(i,k) = xr(i,k) / xcp(i,k)
!         enddo
!       enddo
!
!       do k=1,lev
!         kt=ktt+k-1
!         ku=kuu+k-1
!         do i=1,lons_lat
!           pi = (qqlon(i,ku,lan)*0.005)**xkappa(i,k)
!           qqlon(i,kt,lan) = qqlon(i,kt,lan)*pi
!         enddo
!       enddo
! u v t update at n+1
        do k=1,levs
          kk=levs-k+1
          ku=kuu+k-1
          kv=kvv+k-1
          kt=ktt+k-1
          kp=kup+k-1
!ch>
!       KL=lev-Llist(k)+1
!ch<
          do i=1,lons_lat
!ch         vdzonl(i,kk,lan) = (qqlon(i,ku,lan)-uum(i,kk,lan))*rdt2
!ch         vdmerd(i,kk,lan) = (qqlon(i,kv,lan)-vvm(i,kk,lan))*rdt2
!ch         ddtemp(i,kk,lan) = qqlon(i,kt,lan)
!!          pten(i,kk,lan)    =(qqlon(i,kp,lan)-(dsigma(kk,1)          &
!!                             *ptp(i,lan)+dsigma(kk,2)))*rdt2

            vdzonl(i,k ,lan) = qqlon(i,ku,lan)
            vdmerd(i,k ,lan) = qqlon(i,kv,lan)
            ddtemp(i,k ,lan) = qqlon(i,kt,lan)
!!              pten(i,k ,lan) = (qqlon(i,kp,lan)-pten(i,k,lan))*rdt2
!!          pten(i,k ,lan)    =(qqlon(i,kp,lan)-(dsigma(KL,1)          &
!!                             *ptp_sl(i,lan)+dsigma(KL,2)))*rdt2
          enddo
        enddo
! rq update
        do n=1,ncld
        do k=1,levs
          kk=levs-k+1+(n-1)*levs
          kq=kqq+k-1+(n-1)*levs
          do i=1,lons_lat
!ch         qvadv(i,kk,lan) = qqlon(i,kq,lan)
            qvadv(i,k,n,lan) = qqlon(i,kq,lan)
          enddo
        enddo
        enddo
!        if (myrank .eq. 0 ) then
!        call mymaxmin(qqlon(1,ktt,lan),lons_lat,lonfull,1,' tend temp')
!        endif
!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        ilan=1+jlonf
!        call mymaxmin(grid_gr(ilan,g_u),lons_lat,lonfull,1,' tend u ')
!        call mymaxmin(grid_gr(ilan,g_v),lons_lat,lonfull,1,' tend v ')
!        call mymaxmin(grid_gr(ilan,g_tt),lons_lat,lonfull,1,' tend t ')
!        call mymaxmin(grid_gr(ilan,g_rt),lons_lat,lonfull,1,' tend q ')
!        call mymaxmin(grid_gr(ilan,g_z),lons_lat,lonfull,1,'tend z')
!        do k=1,lev,10
!        kug=g_d+k-1
!        print *,' k=',k
!        call mymaxmin(grid_gr(ilan,kug),lons_lat,lonfull,1,' tend d ')
!        enddo
!        print *,' finish horizonatal advection at lan=',lan
!        endif
!
      enddo
!$omp end parallel do

!
! ===============================
!
      return
      end
! ------------------------------
      subroutine ndslfv_monoadvh2_fgnl(vdzonl,vdmerd,ddtemp    &
                , lonsperlat,deltim,levs,nvars,forward)
!
! a routine to do non-iteration semi-Lagrangain advection
! considering advection  with monotonicity in interpolation
! contact: hann-ming henry juang
! program log
! 2011 02 20 : henry juang, created for ndsl advection
! 2013 06 20 : Henry Juang correct wind direction for north-south advection
!
      use param
      use grid
      use index
      use const
      use rank

      implicit none

      integer nvars
      real(kind=RTYPE) ddtemp(nx,levs,my_max)
      real(kind=RTYPE) vdmerd(nx,levs,my_max),vdzonl(nx,levs,my_max)
      integer,intent(in):: lonsperlat(my)
      real(kind=RTYPE),   intent(in):: deltim

      real(kind=RTYPE)      uulon(lonfull,levs,latpart)
      real(kind=RTYPE)      vvlon(lonfull,levs,latpart)
      real(kind=RTYPE)      qqlon(lonfull,levs*nvars,latpart)
      real(kind=RTYPE)      rrlon(lonfull,levs*nvars,latpart)

      real(kind=RTYPE)      vvlat(latfull,levs,lonpart)
      real(kind=RTYPE)      qqlat(latfull,levs*nvars,lonpart)
      real(kind=RTYPE)      rrlat(latfull,levs*nvars,lonpart)
      real(kind=RTYPE)      xr    (lonfull,levs)
      real(kind=RTYPE)      xcp   (lonfull,levs)
      real(kind=RTYPE)      sumrq (lonfull,levs)
      real(kind=RTYPE)      xkappa(lonfull,levs)
      real(kind=RTYPE)      rdt2, rkt, pi, cons0, cons1, rma, rm2a

!      logical   lprint

      integer mono,mass,levs
      integer nlevs!,levh
      integer i,j,n,k,lon,lan,lat,lons_lat,irc,kk,KL
      integer kuu, kvv, ktt, kup, nqq
      integer ku , kv , kt,  kp
      logical forward
!
!      lprint = .false.

!      if( lprint ) print *,' enter ndslfv_monoadvh '
!
      mono  = 1
      mass  = 0
      cons0 = 0.0
      cons1 = 1.0
!
!      levh = ncld * lev
!
!      kuu = 1
!      kvv = kuu + lev
!      kuu = kvv + lev
      kuu = 1
      kvv = kuu + levs
      ktt = kvv + levs
      kup = ktt + levs

!      nvars = 3
      nlevs = nvars * levs
!
      rdt2 = 0.5 / deltim
!
! =================================================================
!   prepare wind and variable in flux form with gaussina weight
! =================================================================
!

!$omp parallel do                                                      &
!$omp private(lan,lat,lons_lat,rma,rm2a,i,n,kk,k,kv,ku,kt) &
!$omp schedule(dynamic)
      do lan=1,jlistnum

!        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lat = jlist1(lan)
        lons_lat = lonsperlat(lat)
!        rma  = 1. / cosglat(lat) / con_rerth
!        rm2a = rma / cosglat(lat)
        rma  = 1. / cosl(lat) 
        rm2a = rma / cosl(lat)
!
! wind at time step n
        do k=1,levs
          kk=levs-k+1
          do i=1,lons_lat
!ch         uulon(i,k,lan) = ut(i,kk,lan) * rm2a
            uulon(i,k,lan) = ut_sl(i,k ,lan) * rm2a
!hmhj use real(kind=RTYPE) wind
!            vvlon(i,k,lan) = grid_gr(ilan,kvg) * rma
!hmhj use virtual wind
!            vvlon(i,k,lan) = grid_gr(ilan,kvg) / con_rerth
!for CWB GFS
!ch         vvlon(i,k,lan) = vt(i,kk,lan) * rma
            vvlon(i,k,lan) = vt_sl(i,k ,lan) * rma
          enddo
        enddo
!        if( lprint ) then
!          call mymaxmin(uulon(1,1,lan),lons_lat,lonfull,1,' uu1 in deg')
!          call mymaxmin(uulon(1,5,lan),lons_lat,lonfull,1,' uu5 in deg')
!          call mymaxmin(vvlon(1,1,lan),lons_lat,lonfull,1,' vv in deg')
!        endif
!
! u v t at n-1
!
        do k=1,levs
          ku=kuu+k-1
          kv=kvv+k-1
          kt=ktt+k-1
          kk=levs-k+1
!ch>
!       KL=lev-Llist(k)+1
!ch<
          do i=1,lons_lat
!ch         qqlon(i,ku,lan) = uum(i,kk,lan)
!ch         qqlon(i,kv,lan) = vvm(i,kk,lan)
!ch         qqlon(i,kt,lan) = ttm(i,kk,lan)
!!            qqlon(i,kp,lan) = dsigma(kk,1)*ptp(i,lan)+dsigma(kk,2)

            qqlon(i,ku,lan) = vdzonl(i,k,lan)
            qqlon(i,kv,lan) = vdmerd(i,k,lan)
            if ( nvars .ge. 3 ) qqlon(i,kt,lan) = ddtemp(i,k,lan)
!!            qqlon(i,kp,lan) = dsigma(KL,1)*ptp_sl(i,lan)+dsigma(KL,2)
          enddo
        enddo
! add surface pressure perturbation for removing resonance
!       rkt = con_g / ( con_rd * 300.0 )
!       do i=1,lons_lat
!         ilan=i+jlonf
!         qqlon(i,kuu,lan) = log(plev(i,1))+grid_gr(ilan,g_gz)*rkt
!       enddo
!
! save qqlon into n+1 for later as tendency
! z d t at n+1
!        do k=1,lev
!          kk=lev-k+1
!          ku=kuu+k-1
!          kv=kvv+k-1
!          kt=ktt+k-1
!          kug=g_d+k-1
!          kvg=g_z+k-1
!          ktg=g_t+k-1
!          do i=1,lons_lat
!            ilan=i+jlonf
!            grid_gr(ilan,kug) = qqlon(i,ku,lan)
!            grid_gr(ilan,kvg) = qqlon(i,kv,lan)
!            grid_gr(ilan,ktg) = qqlon(i,kt,lan)
!          enddo
!        enddo
! change h to theta
!
! ----- prepare xr, xcp, xkapa
!
!       xr    = cons0
!       xcp   = cons0
!       sumrq = cons0
!
!       do n=1,ntrac
!         nqq = kqq + (n-1)*levs
!         if( ri(n) .ne. cons0 .and. cpi(n) .ne. cons0 ) then
!           do k=1,lev
!             kq=nqq+k-1
!             do i=1,lons_lat
!               xr   (i,k) = xr   (i,k) + qqlon(i,kq,lan)*ri(n)
!               xcp  (i,k) = xcp  (i,k) + qqlon(i,kq,lan)*cpi(n)
!               sumrq(i,k) = sumrq(i,k) + qqlon(i,kq,lan)
!             enddo
!           enddo
!         endif
!       enddo
!       do k=1,lev
!         do i=1,lons_lat
!           xr (i,k)   = ( cons1 - sumrq(i,k) )*ri(0)  + xr (i,k)
!           xcp(i,k)   = ( cons1 - sumrq(i,k) )*cpi(0) + xcp(i,k)
!           xkappa(i,k) = xr(i,k) / xcp(i,k)
!         enddo
!       enddo
!
!       do k=1,lev
!         kt=ktt+k-1
!         ku=kuu+k-1
!         do i=1,lons_lat
!           pi = (qqlon(i,ku,lan)*0.005)**xkappa(i,k)
!           qqlon(i,kt,lan) = qqlon(i,kt,lan)/pi
!         enddo
!       enddo
!
! rq  no need for tendency
!
!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        ilan=1+jlonf
!        call mymaxmin(grid_gr(ilan,g_ttm),lons_lat,lonfull,1,' n-1 t ')
!        call mymaxmin(grid_gr(ilan,g_rrm ),lons_lat,lonfull,1,' n-1 q ')
!        call mymaxmin(grid_gr(ilan,g_zem),lons_lat,lonfull,1,' n-1 z ')
!        call mymaxmin(grid_gr(ilan,g_dim),lons_lat,lonfull,1,' n-1 d ')
!        print *,' ------------------------------------------- '
!        call mymaxmin(qqlon(1,ktt,lan),lons_lat,lonfull,1,' red t ')
!        call mymaxmin(qqlon(1,kqq,lan),lons_lat,lonfull,1,' red q ')
!        call mymaxmin(qqlon(1,kvv,lan),lons_lat,lonfull,1,' red z ')
!        call mymaxmin(qqlon(1,kuu,lan),lons_lat,lonfull,1,' red d ')
!        endif

!       call cyclic_cell_intpx(levs,lons_lat,lonf,uulon(1,1,lan))

!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        call mymaxmin(qqlon(1,ktt,lan),lonfull,lonfull,1,' full t ')
!        call mymaxmin(qqlon(1,kqq,lan),lonfull,lonfull,1,' full q ')
!        call mymaxmin(qqlon(1,kvv,lan),lonfull,lonfull,1,' full z ')
!        call mymaxmin(qqlon(1,kuu,lan),lonfull,lonfull,1,' full d ')
!        endif

        rrlon(:,:,lan) = qqlon(:,:,lan)
!
! first set positive advection in east-west direction
!
        call cyclic_cell_massadvxl(lons_lat,lonfull,levs,nvars,deltim, &
                         uulon(1,1,lan),rrlon(1,1,lan),mass,forward)
!       call cyclic_mono_advectx (lonfull,levs,nvars,deltim,               &
!    &                   uulon(1,1,lan),rrlon(1,1,lan),mono)

        call cyclic_cell_intpx(nlevs,lons_lat,lonfull,rrlon(1,1,lan))
        call cyclic_cell_intpx(nlevs,lons_lat,lonfull,qqlon(1,1,lan))
        call cyclic_cell_intpx(levs,lons_lat,lonfull,vvlon(1,1,lan))

!        if( lprint ) then
!        print *,' done cyclic_massadvx with mass= ',mass
!        print *,' ------------------------------------------- '
!        call mymaxmin(rrlon(1,ktt,lan),lonfull,lonfull,1,' advx t ')
!        call mymaxmin(rrlon(1,kqq,lan),lonfull,lonfull,1,' advx q ')
!        call mymaxmin(rrlon(1,kvv,lan),lonfull,lonfull,1,' advx z ')
!        call mymaxmin(rrlon(1,kuu,lan),lonfull,lonfull,1,' advx d ')
!        print *,' done the first x adv for lan=',lan
!        print *,' =========================================== '
!        endif

      enddo
!$omp end parallel do

! ---------------------------------------------------------------------
! mpi para from east-west full grid to north-south full grid
! ---------------------------------------------------------------------
!
! para vvlon, qqlon, and rrlon to vvlat, qqlat, rrlat
!       if( lprint ) print *,' ndslfv_advect transport from we to ns '

       call para_we2ns(vvlon,vvlat,levs,my)
       call para_we2ns(rrlon,rrlat,nlevs,my)
       call para_we2ns(qqlon,qqlat,nlevs,my)

!       if( lprint ) then
!       print *,' ------------ after we2ns ---------------------- '
!       do lon=1,mylonlen
!        print *,'  lon=',lon
!        call mymaxmin(vvlat(1,1  ,lon),latfull,latfull,1,' we2ns v')
!        call mymaxmin(rrlat(1,kqq,lon),latfull,latfull,1,' we2ns r')
!        call mymaxmin(qqlat(1,kqq,lon),latfull,latfull,1,' we2ns q')
!       enddo
!       endif
!
! ---------------------------------------------------------------------
! -------------- in north-soutn great circle -------------------
! ---------------------------------------------------------------------

!      if ( myrank .eq. 0 ) lprint = .true.
!       if( lprint ) then
!       print *,' ndslfv_advect adv loop in y '
!       print *,' mylonlen=',mylonlen
!       endif

!$omp parallel do private(lon,k,j,ku,kv) &
!$omp schedule(dynamic)
       do lon=1,mylonlen
!
! convert wind before advy
        do k=1,levs
          ku=kuu+k-1
          kv=kvv+k-1
          do j=1,lathalf
            rrlat(j,ku,lon) = -rrlat(j,ku,lon)
            rrlat(j,kv,lon) = -rrlat(j,kv,lon)
            qqlat(j,ku,lon) = -qqlat(j,ku,lon)
            qqlat(j,kv,lon) = -qqlat(j,kv,lon)
          enddo
        enddo
!
!        if( lprint ) print *,' lon=',lon

! first set advection in north-south direction in great circle through two poles
!
!        call fixend_cell_massadvy(latfull,lathalf,levs,nvars,deltim, &
!                         vvlat(1,1,lon),rrlat(1,1,lon),mass)
        call cyclic_cell_massadvyl(latfull,levs,nvars,deltim,         &
                         vvlat(1,1,lon),rrlat(1,1,lon),mass,forward)
!       call cyclic_mono_advecty (latfull,levs,nvars,deltim,
!    &                   vvlat(1,1,lon),rrlat(1,1,lon),mono)

!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        call mymaxmin(rrlat(1,kuu,lon),latfull,latfull,1,' advry u ')
!        call mymaxmin(rrlat(1,kvv,lon),latfull,latfull,1,' advry v ')
!        call mymaxmin(rrlat(1,ktt,lon),latfull,latfull,1,' advry t ')
!        call mymaxmin(rrlat(1,kqq,lon),latfull,latfull,1,' advry q ')
!        call mymaxmin(rrlat(1,kvv,lon),latfull,latfull,1,' advry z ')
!        call mymaxmin(rrlat(1,kuu,lon),latfull,latfull,1,' advry d ')
!        endif
!
! second set advection in north-south direction in great circle through two poles
!
!        call fixend_cell_massadvy(latfull,lathalf,levs,nvars,deltim, &
!                         vvlat(1,1,lon),qqlat(1,1,lon),mass)
        call cyclic_cell_massadvyl(latfull,levs,nvars,deltim,         &
                         vvlat(1,1,lon),qqlat(1,1,lon),mass,forward)
!       call cyclic_mono_advecty (my,levs,nvars,deltim,
!    &                   vvlat(1,1,lon),qqlat(1,1,lon),mono)

!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        call mymaxmin(qqlat(1,kuu,lon),latfull,latfull,1,' advqy u ')
!        call mymaxmin(qqlat(1,kvv,lon),latfull,latfull,1,' advqy v ')
!        call mymaxmin(qqlat(1,ktt,lon),latfull,latfull,1,' advqy t ')
!        call mymaxmin(qqlat(1,kqq,lon),latfull,latfull,1,' advqy q ')
!        call mymaxmin(qqlat(1,kvv,lon),latfull,latfull,1,' advqy z ')
!        call mymaxmin(qqlat(1,kuu,lon),latfull,latfull,1,' advqy d ')
!        print *,' done with y at lon=',lon
!        endif
! convert wind back after advy
        do k=1,levs
          ku=kuu+k-1
          kv=kvv+k-1
          do j=1,lathalf
            rrlat(j,ku,lon) = -rrlat(j,ku,lon)
            rrlat(j,kv,lon) = -rrlat(j,kv,lon)
            qqlat(j,ku,lon) = -qqlat(j,ku,lon)
            qqlat(j,kv,lon) = -qqlat(j,kv,lon)
          enddo
        enddo

       enddo
!$omp end parallel do
!
! ----------------------------------------------------------------------
! mpi para from north-south direction to east-west direeectory
! ----------------------------------------------------------------------
!
! para qqlat and rrlat to qqlon and rrlon

!       if( lprint ) print *,' ndslfv_advect transport from ns to we '

       call para_ns2we(rrlat,rrlon,nlevs,my)
       call para_ns2we(qqlat,qqlon,nlevs,my)
!       if( lprint ) then
!       print *,' ------------ after ns2we ---------------------- '
!       do lan=1,jlistnum
!        print *,'  lan=',lan
!        call mymaxmin(rrlon(1,kqq,lan),lonfull,lonfull,1,' ns2we r')
!        call mymaxmin(qqlon(1,kqq,lan),lonfull,lonfull,1,' ns2we q')
!       enddo
!       endif

! ---------------------------------------------------------------
! ---------------- back to east-west direction ------------------
! ---------------------------------------------------------------
!      print *,' ndslfv_advect adv loop in x for last '

!$omp parallel do                                                &
!$omp private(lan,lat,lons_lat,i,k,kk,n,kv,ku,kt)    &
!$omp schedule(dynamic)

      do lan=1,jlistnum

!        lat = global_lats_a(ipt_jlistnum-1+lan)
        lat = jlist1(lan)
        lons_lat = lonsperlat(lat)
!
! mass conserving interpolation from full grid to reduced grid
!
        call cyclic_cell_intpx(nlevs,lonfull,lons_lat,qqlon(1,1,lan))
        call cyclic_cell_intpx(nlevs,lonfull,lons_lat,rrlon(1,1,lan))

!
! second set advection in x for the second of the pair
!
        call cyclic_cell_massadvxl(lons_lat,lonfull,levs,nvars,deltim, &
                         uulon(1,1,lan),qqlon(1,1,lan),mass,forward)
!       call cyclic_mono_advectx (lonfull,levs,nvars,deltim,               &
!    &                   uulon(1,1,lan),qqlon(1,1,lan),mono)
!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        call mymaxmin(qqlon(1,kqq,lan),lonfull,lonfull,1,' adv x q ')
!        endif

        do k=1,nlevs
          do i=1,lons_lat
            qqlon(i,k,lan) = 0.5 * ( qqlon(i,k,lan) + rrlon(i,k,lan) )
          enddo
        enddo

!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        call mymaxmin(qqlon(1,kuu,lan),lons_lat,lonfull,1,' do redu u ')
!        call mymaxmin(qqlon(1,kvv,lan),lons_lat,lonfull,1,' do redu v ')
!        call mymaxmin(qqlon(1,ktt,lan),lons_lat,lonfull,1,' do redu t ')
!        call mymaxmin(qqlon(1,kqq,lan),lons_lat,lonfull,1,' do redu q ')
!        call mymaxmin(qqlon(1,kvv,lan),lons_lat,lonfull,1,' do redu z')
!        call mymaxmin(qqlon(1,kuu,lan),lons_lat,lonfull,1,' do redu d')
!        print *,' finish horizonatal advection at lan=',lan
!        endif

! change theta to h
!
! ----- prepare xr, xcp, xkapa
!
!       xr    = cons0
!       xcp   = cons0
!       sumrq = cons0
!
!       do n=1,ntrac
!         nqq = kqq + (n-1)*levs
!         if( ri(n) .ne. cons0 .and. cpi(n) .ne. cons0 ) then
!           do k=1,lev
!             kq=nqq+k-1
!             do i=1,lons_lat
!               xr   (i,k) = xr   (i,k) + qqlon(i,kq,lan)*ri(n)
!               xcp  (i,k) = xcp  (i,k) + qqlon(i,kq,lan)*cpi(n)
!               sumrq(i,k) = sumrq(i,k) + qqlon(i,kq,lan)
!             enddo
!           enddo
!         endif
!       enddo
!       do k=1,lev
!         do i=1,lons_lat
!           xr (i,k)   = ( cons1 - sumrq(i,k) )*ri(0)  + xr (i,k)
!           xcp(i,k)   = ( cons1 - sumrq(i,k) )*cpi(0) + xcp(i,k)
!           xkappa(i,k) = xr(i,k) / xcp(i,k)
!         enddo
!       enddo
!
!       do k=1,lev
!         kt=ktt+k-1
!         ku=kuu+k-1
!         do i=1,lons_lat
!           pi = (qqlon(i,ku,lan)*0.005)**xkappa(i,k)
!           qqlon(i,kt,lan) = qqlon(i,kt,lan)*pi
!         enddo
!       enddo
! u v t update at n+1
        do k=1,levs
          kk=levs-k+1
          ku=kuu+k-1
          kv=kvv+k-1
          kt=ktt+k-1
!!          kp=kup+k-1
!ch>
!       KL=lev-Llist(k)+1:w
!ch<
          do i=1,lons_lat
!ch         vdzonl(i,kk,lan) = (qqlon(i,ku,lan)-uum(i,kk,lan))*rdt2
!ch         vdmerd(i,kk,lan) = (qqlon(i,kv,lan)-vvm(i,kk,lan))*rdt2
!ch         ddtemp(i,kk,lan) = qqlon(i,kt,lan)
!ch         pten(i,kk,lan)    =(qqlon(i,kp,lan)-(dsigma(kk,1)          &
!ch                            *ptp(i,lan)+dsigma(kk,2)))*rdt2

            vdzonl(i,k ,lan) = qqlon(i,ku,lan)
            vdmerd(i,k ,lan) = qqlon(i,kv,lan)
            if ( nvars .ge. 3 ) ddtemp(i,k ,lan) = qqlon(i,kt,lan)
!            pten(i,k ,lan)    =(qqlon(i,kp,lan)-(dsigma(KL,1)          &
!                               *ptp_sl(i,lan)+dsigma(KL,2)))*rdt2
          enddo
        enddo
!        if (myrank .eq. 0 ) then
!        call mymaxmin(qqlon(1,ktt,lan),lons_lat,lonfull,1,' tend temp')
!        endif
!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        ilan=1+jlonf
!        call mymaxmin(grid_gr(ilan,g_u),lons_lat,lonfull,1,' tend u ')
!        call mymaxmin(grid_gr(ilan,g_v),lons_lat,lonfull,1,' tend v ')
!        call mymaxmin(grid_gr(ilan,g_tt),lons_lat,lonfull,1,' tend t ')
!        call mymaxmin(grid_gr(ilan,g_rt),lons_lat,lonfull,1,' tend q ')
!        call mymaxmin(grid_gr(ilan,g_z),lons_lat,lonfull,1,'tend z')
!        do k=1,lev,10
!        kug=g_d+k-1
!        print *,' k=',k
!        call mymaxmin(grid_gr(ilan,kug),lons_lat,lonfull,1,' tend d ')
!        enddo
!        print *,' finish horizonatal advection at lan=',lan
!        endif
!
      enddo
!$omp end parallel do

!
! ===============================
!
      return
      end
! ------------------------------
      subroutine ndslfv_monoadvh2_fgnl_xy(vdzonl,vdmerd,ddtemp    &
                , lonsperlat,deltim,levs,nvars,forward)
!
! a routine to do non-iteration semi-Lagrangain advection
! considering advection  with monotonicity in interpolation
! contact: hann-ming henry juang
! program log
! 2011 02 20 : henry juang, created for ndsl advection
! 2013 06 20 : Henry Juang correct wind direction for north-south advection
!
      use param
      use grid
      use index
      use const
      use rank

      implicit none

      integer nvars
      real(kind=RTYPE) ddtemp(nx,levs,my_max)
      real(kind=RTYPE) vdmerd(nx,levs,my_max),vdzonl(nx,levs,my_max)
      integer,intent(in):: lonsperlat(my)
      real(kind=RTYPE),   intent(in):: deltim

      real(kind=RTYPE)      uulon(lonfull,levs,latpart)
      real(kind=RTYPE)      vvlon(lonfull,levs,latpart)
      real(kind=RTYPE)      qqlon(lonfull,levs*nvars,latpart)
!      real(kind=RTYPE)      rrlon(lonfull,levs*nvars,latpart)

      real(kind=RTYPE)      vvlat(latfull,levs,lonpart)
      real(kind=RTYPE)      qqlat(latfull,levs*nvars,lonpart)
!      real(kind=RTYPE)      rrlat(latfull,levs*nvars,lonpart)
      real(kind=RTYPE)      xr    (lonfull,levs)
      real(kind=RTYPE)      xcp   (lonfull,levs)
      real(kind=RTYPE)      sumrq (lonfull,levs)
      real(kind=RTYPE)      xkappa(lonfull,levs)
      real(kind=RTYPE)      rdt2, rkt, pi, cons0, cons1, rma, rm2a

!      logical   lprint

      integer mono,mass,levs
      integer nlevs!,levh
      integer i,j,n,k,lon,lan,lat,lons_lat,irc,kk,KL
      integer kuu, kvv, ktt, kup, nqq
      integer ku , kv , kt,  kp
      logical forward
!
!      lprint = .false.

!      if( lprint ) print *,' enter ndslfv_monoadvh '
!
      mono  = 1
      mass  = 0
      cons0 = 0.0
      cons1 = 1.0
!
!      levh = ncld * lev
!
!      kuu = 1
!      kvv = kuu + lev
!      kuu = kvv + lev
      kuu = 1
      kvv = kuu + levs
      ktt = kvv + levs
      kup = ktt + levs

!      nvars = 3
      nlevs = nvars * levs
!
      rdt2 = 0.5 / deltim
!
! =================================================================
!   prepare wind and variable in flux form with gaussina weight
! =================================================================
!

!$omp parallel do                                                      &
!$omp private(lan,lat,lons_lat,rma,rm2a,i,n,kk,k,kv,ku,kt) &
!$omp schedule(dynamic)
      do lan=1,jlistnum

!        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lat = jlist1(lan)
        lons_lat = lonsperlat(lat)
!        rma  = 1. / cosglat(lat) / con_rerth
!        rm2a = rma / cosglat(lat)
        rma  = 1. / cosl(lat) 
        rm2a = rma / cosl(lat)
!
! wind at time step n
        do k=1,levs
          kk=levs-k+1
          do i=1,lons_lat
!ch         uulon(i,k,lan) = ut(i,kk,lan) * rm2a
            uulon(i,k,lan) = ut_sl(i,k ,lan) * rm2a
!hmhj use real(kind=RTYPE) wind
!            vvlon(i,k,lan) = grid_gr(ilan,kvg) * rma
!hmhj use virtual wind
!            vvlon(i,k,lan) = grid_gr(ilan,kvg) / con_rerth
!for CWB GFS
!ch         vvlon(i,k,lan) = vt(i,kk,lan) * rma
            vvlon(i,k,lan) = vt_sl(i,k ,lan) * rma
          enddo
        enddo
!        if( lprint ) then
!          call mymaxmin(uulon(1,1,lan),lons_lat,lonfull,1,' uu1 in deg')
!          call mymaxmin(uulon(1,5,lan),lons_lat,lonfull,1,' uu5 in deg')
!          call mymaxmin(vvlon(1,1,lan),lons_lat,lonfull,1,' vv in deg')
!        endif
!
! u v t at n-1
!
        do k=1,levs
          ku=kuu+k-1
          kv=kvv+k-1
          kt=ktt+k-1
          kk=levs-k+1
!ch>
!       KL=lev-Llist(k)+1
!ch<
          do i=1,lons_lat
!ch         qqlon(i,ku,lan) = uum(i,kk,lan)
!ch         qqlon(i,kv,lan) = vvm(i,kk,lan)
!ch         qqlon(i,kt,lan) = ttm(i,kk,lan)
!!            qqlon(i,kp,lan) = dsigma(kk,1)*ptp(i,lan)+dsigma(kk,2)

            qqlon(i,ku,lan) = vdzonl(i,k,lan)
            qqlon(i,kv,lan) = vdmerd(i,k,lan)
            if ( nvars .ge. 3 ) qqlon(i,kt,lan) = ddtemp(i,k,lan)
!!            qqlon(i,kp,lan) = dsigma(KL,1)*ptp_sl(i,lan)+dsigma(KL,2)
          enddo
        enddo
! add surface pressure perturbation for removing resonance
!       rkt = con_g / ( con_rd * 300.0 )
!       do i=1,lons_lat
!         ilan=i+jlonf
!         qqlon(i,kuu,lan) = log(plev(i,1))+grid_gr(ilan,g_gz)*rkt
!       enddo
!
! save qqlon into n+1 for later as tendency
! z d t at n+1
!        do k=1,lev
!          kk=lev-k+1
!          ku=kuu+k-1
!          kv=kvv+k-1
!          kt=ktt+k-1
!          kug=g_d+k-1
!          kvg=g_z+k-1
!          ktg=g_t+k-1
!          do i=1,lons_lat
!            ilan=i+jlonf
!            grid_gr(ilan,kug) = qqlon(i,ku,lan)
!            grid_gr(ilan,kvg) = qqlon(i,kv,lan)
!            grid_gr(ilan,ktg) = qqlon(i,kt,lan)
!          enddo
!        enddo
! change h to theta
!
! ----- prepare xr, xcp, xkapa
!
!       xr    = cons0
!       xcp   = cons0
!       sumrq = cons0
!
!       do n=1,ntrac
!         nqq = kqq + (n-1)*levs
!         if( ri(n) .ne. cons0 .and. cpi(n) .ne. cons0 ) then
!           do k=1,lev
!             kq=nqq+k-1
!             do i=1,lons_lat
!               xr   (i,k) = xr   (i,k) + qqlon(i,kq,lan)*ri(n)
!               xcp  (i,k) = xcp  (i,k) + qqlon(i,kq,lan)*cpi(n)
!               sumrq(i,k) = sumrq(i,k) + qqlon(i,kq,lan)
!             enddo
!           enddo
!         endif
!       enddo
!       do k=1,lev
!         do i=1,lons_lat
!           xr (i,k)   = ( cons1 - sumrq(i,k) )*ri(0)  + xr (i,k)
!           xcp(i,k)   = ( cons1 - sumrq(i,k) )*cpi(0) + xcp(i,k)
!           xkappa(i,k) = xr(i,k) / xcp(i,k)
!         enddo
!       enddo
!
!       do k=1,lev
!         kt=ktt+k-1
!         ku=kuu+k-1
!         do i=1,lons_lat
!           pi = (qqlon(i,ku,lan)*0.005)**xkappa(i,k)
!           qqlon(i,kt,lan) = qqlon(i,kt,lan)/pi
!         enddo
!       enddo
!
! rq  no need for tendency
!
!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        ilan=1+jlonf
!        call mymaxmin(grid_gr(ilan,g_ttm),lons_lat,lonfull,1,' n-1 t ')
!        call mymaxmin(grid_gr(ilan,g_rrm ),lons_lat,lonfull,1,' n-1 q ')
!        call mymaxmin(grid_gr(ilan,g_zem),lons_lat,lonfull,1,' n-1 z ')
!        call mymaxmin(grid_gr(ilan,g_dim),lons_lat,lonfull,1,' n-1 d ')
!        print *,' ------------------------------------------- '
!        call mymaxmin(qqlon(1,ktt,lan),lons_lat,lonfull,1,' red t ')
!        call mymaxmin(qqlon(1,kqq,lan),lons_lat,lonfull,1,' red q ')
!        call mymaxmin(qqlon(1,kvv,lan),lons_lat,lonfull,1,' red z ')
!        call mymaxmin(qqlon(1,kuu,lan),lons_lat,lonfull,1,' red d ')
!        endif

!       call cyclic_cell_intpx(levs,lons_lat,lonf,uulon(1,1,lan))

!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        call mymaxmin(qqlon(1,ktt,lan),lonfull,lonfull,1,' full t ')
!        call mymaxmin(qqlon(1,kqq,lan),lonfull,lonfull,1,' full q ')
!        call mymaxmin(qqlon(1,kvv,lan),lonfull,lonfull,1,' full z ')
!        call mymaxmin(qqlon(1,kuu,lan),lonfull,lonfull,1,' full d ')
!        endif

!        rrlon(:,:,lan) = qqlon(:,:,lan)
!
! first set positive advection in east-west direction
!
        call cyclic_cell_massadvxl(lons_lat,lonfull,levs,nvars,deltim, &
                         uulon(1,1,lan),qqlon(1,1,lan),mass,forward)
!       call cyclic_mono_advectx (lonfull,levs,nvars,deltim,               &
!    &                   uulon(1,1,lan),rrlon(1,1,lan),mono)

!        call cyclic_cell_intpx(nlevs,lons_lat,lonfull,rrlon(1,1,lan))
        call cyclic_cell_intpx(nlevs,lons_lat,lonfull,qqlon(1,1,lan))
        call cyclic_cell_intpx(levs,lons_lat,lonfull,vvlon(1,1,lan))

!        if( lprint ) then
!        print *,' done cyclic_massadvx with mass= ',mass
!        print *,' ------------------------------------------- '
!        call mymaxmin(rrlon(1,ktt,lan),lonfull,lonfull,1,' advx t ')
!        call mymaxmin(rrlon(1,kqq,lan),lonfull,lonfull,1,' advx q ')
!        call mymaxmin(rrlon(1,kvv,lan),lonfull,lonfull,1,' advx z ')
!        call mymaxmin(rrlon(1,kuu,lan),lonfull,lonfull,1,' advx d ')
!        print *,' done the first x adv for lan=',lan
!        print *,' =========================================== '
!        endif

      enddo
!$omp end parallel do

! ---------------------------------------------------------------------
! mpi para from east-west full grid to north-south full grid
! ---------------------------------------------------------------------
!
! para vvlon, qqlon, and rrlon to vvlat, qqlat, rrlat
!       if( lprint ) print *,' ndslfv_advect transport from we to ns '

       call para_we2ns(vvlon,vvlat,levs,my)
!       call para_we2ns(rrlon,rrlat,nlevs,my)
       call para_we2ns(qqlon,qqlat,nlevs,my)

!       if( lprint ) then
!       print *,' ------------ after we2ns ---------------------- '
!       do lon=1,mylonlen
!        print *,'  lon=',lon
!        call mymaxmin(vvlat(1,1  ,lon),latfull,latfull,1,' we2ns v')
!        call mymaxmin(rrlat(1,kqq,lon),latfull,latfull,1,' we2ns r')
!        call mymaxmin(qqlat(1,kqq,lon),latfull,latfull,1,' we2ns q')
!       enddo
!       endif
!
! ---------------------------------------------------------------------
! -------------- in north-soutn great circle -------------------
! ---------------------------------------------------------------------

!      if ( myrank .eq. 0 ) lprint = .true.
!       if( lprint ) then
!       print *,' ndslfv_advect adv loop in y '
!       print *,' mylonlen=',mylonlen
!       endif

!$omp parallel do private(lon,k,j,ku,kv) &
!$omp schedule(dynamic)
       do lon=1,mylonlen
!
! convert wind before advy
        do k=1,levs
          ku=kuu+k-1
          kv=kvv+k-1
          do j=1,lathalf
!            rrlat(j,ku,lon) = -rrlat(j,ku,lon)
!            rrlat(j,kv,lon) = -rrlat(j,kv,lon)
            qqlat(j,ku,lon) = -qqlat(j,ku,lon)
            qqlat(j,kv,lon) = -qqlat(j,kv,lon)
          enddo
        enddo
!
!        if( lprint ) print *,' lon=',lon

! first set advection in north-south direction in great circle through two poles
!
!        call fixend_cell_massadvy(latfull,lathalf,levs,nvars,deltim, &
!                         vvlat(1,1,lon),rrlat(1,1,lon),mass)
!!        call cyclic_cell_massadvyl(latfull,levs,nvars,deltim,         &
!!                         vvlat(1,1,lon),rrlat(1,1,lon),mass,forward)
!       call cyclic_mono_advecty (latfull,levs,nvars,deltim,
!    &                   vvlat(1,1,lon),rrlat(1,1,lon),mono)

!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        call mymaxmin(rrlat(1,kuu,lon),latfull,latfull,1,' advry u ')
!        call mymaxmin(rrlat(1,kvv,lon),latfull,latfull,1,' advry v ')
!        call mymaxmin(rrlat(1,ktt,lon),latfull,latfull,1,' advry t ')
!        call mymaxmin(rrlat(1,kqq,lon),latfull,latfull,1,' advry q ')
!        call mymaxmin(rrlat(1,kvv,lon),latfull,latfull,1,' advry z ')
!        call mymaxmin(rrlat(1,kuu,lon),latfull,latfull,1,' advry d ')
!        endif
!
! second set advection in north-south direction in great circle through two poles
!
!        call fixend_cell_massadvy(latfull,lathalf,levs,nvars,deltim, &
!                         vvlat(1,1,lon),qqlat(1,1,lon),mass)
        call cyclic_cell_massadvyl(latfull,levs,nvars,deltim,         &
                         vvlat(1,1,lon),qqlat(1,1,lon),mass,forward)
!       call cyclic_mono_advecty (my,levs,nvars,deltim,
!    &                   vvlat(1,1,lon),qqlat(1,1,lon),mono)

!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        call mymaxmin(qqlat(1,kuu,lon),latfull,latfull,1,' advqy u ')
!        call mymaxmin(qqlat(1,kvv,lon),latfull,latfull,1,' advqy v ')
!        call mymaxmin(qqlat(1,ktt,lon),latfull,latfull,1,' advqy t ')
!        call mymaxmin(qqlat(1,kqq,lon),latfull,latfull,1,' advqy q ')
!        call mymaxmin(qqlat(1,kvv,lon),latfull,latfull,1,' advqy z ')
!        call mymaxmin(qqlat(1,kuu,lon),latfull,latfull,1,' advqy d ')
!        print *,' done with y at lon=',lon
!        endif
! convert wind back after advy
        do k=1,levs
          ku=kuu+k-1
          kv=kvv+k-1
          do j=1,lathalf
!!            rrlat(j,ku,lon) = -rrlat(j,ku,lon)
!!            rrlat(j,kv,lon) = -rrlat(j,kv,lon)
            qqlat(j,ku,lon) = -qqlat(j,ku,lon)
            qqlat(j,kv,lon) = -qqlat(j,kv,lon)
          enddo
        enddo

       enddo
!$omp end parallel do
!
! ----------------------------------------------------------------------
! mpi para from north-south direction to east-west direeectory
! ----------------------------------------------------------------------
!
! para qqlat and rrlat to qqlon and rrlon

!       if( lprint ) print *,' ndslfv_advect transport from ns to we '

!       call para_ns2we(rrlat,rrlon,nlevs,my)
       call para_ns2we(qqlat,qqlon,nlevs,my)
!       if( lprint ) then
!       print *,' ------------ after ns2we ---------------------- '
!       do lan=1,jlistnum
!        print *,'  lan=',lan
!        call mymaxmin(rrlon(1,kqq,lan),lonfull,lonfull,1,' ns2we r')
!        call mymaxmin(qqlon(1,kqq,lan),lonfull,lonfull,1,' ns2we q')
!       enddo
!       endif

! ---------------------------------------------------------------
! ---------------- back to east-west direction ------------------
! ---------------------------------------------------------------
!      print *,' ndslfv_advect adv loop in x for last '

!$omp parallel do                                                &
!$omp private(lan,lat,lons_lat,i,k,kk,n,kv,ku,kt)    &
!$omp schedule(dynamic)

      do lan=1,jlistnum

!        lat = global_lats_a(ipt_jlistnum-1+lan)
        lat = jlist1(lan)
        lons_lat = lonsperlat(lat)
!
! mass conserving interpolation from full grid to reduced grid
!
        call cyclic_cell_intpx(nlevs,lonfull,lons_lat,qqlon(1,1,lan))
!!        call cyclic_cell_intpx(nlevs,lonfull,lons_lat,rrlon(1,1,lan))

!
! second set advection in x for the second of the pair
!
!!        call cyclic_cell_massadvxl(lons_lat,lonfull,levs,nvars,deltim, &
!!                         uulon(1,1,lan),qqlon(1,1,lan),mass,forward)
!       call cyclic_mono_advectx (lonfull,levs,nvars,deltim,               &
!    &                   uulon(1,1,lan),qqlon(1,1,lan),mono)
!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        call mymaxmin(qqlon(1,kqq,lan),lonfull,lonfull,1,' adv x q ')
!        endif

!!        do k=1,nlevs
!!          do i=1,lons_lat
!!            qqlon(i,k,lan) = 0.5 * ( qqlon(i,k,lan) + rrlon(i,k,lan) )
!!          enddo
!!        enddo

!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        call mymaxmin(qqlon(1,kuu,lan),lons_lat,lonfull,1,' do redu u ')
!        call mymaxmin(qqlon(1,kvv,lan),lons_lat,lonfull,1,' do redu v ')
!        call mymaxmin(qqlon(1,ktt,lan),lons_lat,lonfull,1,' do redu t ')
!        call mymaxmin(qqlon(1,kqq,lan),lons_lat,lonfull,1,' do redu q ')
!        call mymaxmin(qqlon(1,kvv,lan),lons_lat,lonfull,1,' do redu z')
!        call mymaxmin(qqlon(1,kuu,lan),lons_lat,lonfull,1,' do redu d')
!        print *,' finish horizonatal advection at lan=',lan
!        endif

! change theta to h
!
! ----- prepare xr, xcp, xkapa
!
!       xr    = cons0
!       xcp   = cons0
!       sumrq = cons0
!
!       do n=1,ntrac
!         nqq = kqq + (n-1)*levs
!         if( ri(n) .ne. cons0 .and. cpi(n) .ne. cons0 ) then
!           do k=1,lev
!             kq=nqq+k-1
!             do i=1,lons_lat
!               xr   (i,k) = xr   (i,k) + qqlon(i,kq,lan)*ri(n)
!               xcp  (i,k) = xcp  (i,k) + qqlon(i,kq,lan)*cpi(n)
!               sumrq(i,k) = sumrq(i,k) + qqlon(i,kq,lan)
!             enddo
!           enddo
!         endif
!       enddo
!       do k=1,lev
!         do i=1,lons_lat
!           xr (i,k)   = ( cons1 - sumrq(i,k) )*ri(0)  + xr (i,k)
!           xcp(i,k)   = ( cons1 - sumrq(i,k) )*cpi(0) + xcp(i,k)
!           xkappa(i,k) = xr(i,k) / xcp(i,k)
!         enddo
!       enddo
!
!       do k=1,lev
!         kt=ktt+k-1
!         ku=kuu+k-1
!         do i=1,lons_lat
!           pi = (qqlon(i,ku,lan)*0.005)**xkappa(i,k)
!           qqlon(i,kt,lan) = qqlon(i,kt,lan)*pi
!         enddo
!       enddo
! u v t update at n+1
        do k=1,levs
          kk=levs-k+1
          ku=kuu+k-1
          kv=kvv+k-1
          kt=ktt+k-1
!!          kp=kup+k-1
!ch>
!       KL=lev-Llist(k)+1:w
!ch<
          do i=1,lons_lat
!ch         vdzonl(i,kk,lan) = (qqlon(i,ku,lan)-uum(i,kk,lan))*rdt2
!ch         vdmerd(i,kk,lan) = (qqlon(i,kv,lan)-vvm(i,kk,lan))*rdt2
!ch         ddtemp(i,kk,lan) = qqlon(i,kt,lan)
!ch         pten(i,kk,lan)    =(qqlon(i,kp,lan)-(dsigma(kk,1)          &
!ch                            *ptp(i,lan)+dsigma(kk,2)))*rdt2

            vdzonl(i,k ,lan) = qqlon(i,ku,lan)
            vdmerd(i,k ,lan) = qqlon(i,kv,lan)
            if ( nvars .ge. 3 ) ddtemp(i,k ,lan) = qqlon(i,kt,lan)
!            pten(i,k ,lan)    =(qqlon(i,kp,lan)-(dsigma(KL,1)          &
!                               *ptp_sl(i,lan)+dsigma(KL,2)))*rdt2
          enddo
        enddo
!        if (myrank .eq. 0 ) then
!        call mymaxmin(qqlon(1,ktt,lan),lons_lat,lonfull,1,' tend temp')
!        endif
!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        ilan=1+jlonf
!        call mymaxmin(grid_gr(ilan,g_u),lons_lat,lonfull,1,' tend u ')
!        call mymaxmin(grid_gr(ilan,g_v),lons_lat,lonfull,1,' tend v ')
!        call mymaxmin(grid_gr(ilan,g_tt),lons_lat,lonfull,1,' tend t ')
!        call mymaxmin(grid_gr(ilan,g_rt),lons_lat,lonfull,1,' tend q ')
!        call mymaxmin(grid_gr(ilan,g_z),lons_lat,lonfull,1,'tend z')
!        do k=1,lev,10
!        kug=g_d+k-1
!        print *,' k=',k
!        call mymaxmin(grid_gr(ilan,kug),lons_lat,lonfull,1,' tend d ')
!        enddo
!        print *,' finish horizonatal advection at lan=',lan
!        endif
!
      enddo
!$omp end parallel do

!
! ===============================
!
      return
      end
! ------------------------------
      subroutine ndslfv_monoadvh2_fgnl_yx(vdzonl,vdmerd,ddtemp    &
                , lonsperlat,deltim,levs,nvars,forward)
!
! a routine to do non-iteration semi-Lagrangain advection
! considering advection  with monotonicity in interpolation
! contact: hann-ming henry juang
! program log
! 2011 02 20 : henry juang, created for ndsl advection
! 2013 06 20 : Henry Juang correct wind direction for north-south advection
!
      use param
      use grid
      use index
      use const
      use rank

      implicit none

      integer nvars
      real(kind=RTYPE) ddtemp(nx,levs,my_max)
      real(kind=RTYPE) vdmerd(nx,levs,my_max),vdzonl(nx,levs,my_max)
      integer,intent(in):: lonsperlat(my)
      real(kind=RTYPE),   intent(in):: deltim

      real(kind=RTYPE)      uulon(lonfull,levs,latpart)
      real(kind=RTYPE)      vvlon(lonfull,levs,latpart)
      real(kind=RTYPE)      qqlon(lonfull,levs*nvars,latpart)
!      real(kind=RTYPE)      rrlon(lonfull,levs*nvars,latpart)

      real(kind=RTYPE)      vvlat(latfull,levs,lonpart)
      real(kind=RTYPE)      qqlat(latfull,levs*nvars,lonpart)
!      real(kind=RTYPE)      rrlat(latfull,levs*nvars,lonpart)
      real(kind=RTYPE)      xr    (lonfull,levs)
      real(kind=RTYPE)      xcp   (lonfull,levs)
      real(kind=RTYPE)      sumrq (lonfull,levs)
      real(kind=RTYPE)      xkappa(lonfull,levs)
      real(kind=RTYPE)      rdt2, rkt, pi, cons0, cons1, rma, rm2a

!      logical   lprint

      integer mono,mass,levs
      integer nlevs!,levh
      integer i,j,n,k,lon,lan,lat,lons_lat,irc,kk,KL
      integer kuu, kvv, ktt, kup, nqq
      integer ku , kv , kt,  kp
      logical forward
!
!      lprint = .false.

!      if( lprint ) print *,' enter ndslfv_monoadvh '
!
      mono  = 1
      mass  = 0
      cons0 = 0.0
      cons1 = 1.0
!
!      levh = ncld * lev
!
!      kuu = 1
!      kvv = kuu + lev
!      kuu = kvv + lev
      kuu = 1
      kvv = kuu + levs
      ktt = kvv + levs
      kup = ktt + levs

!      nvars = 3
      nlevs = nvars * levs
!
      rdt2 = 0.5 / deltim
!
! =================================================================
!   prepare wind and variable in flux form with gaussina weight
! =================================================================
!

!$omp parallel do                                                      &
!$omp private(lan,lat,lons_lat,rma,rm2a,i,n,kk,k,kv,ku,kt) &
!$omp schedule(dynamic)
      do lan=1,jlistnum

!        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lat = jlist1(lan)
        lons_lat = lonsperlat(lat)
!        rma  = 1. / cosglat(lat) / con_rerth
!        rm2a = rma / cosglat(lat)
        rma  = 1. / cosl(lat) 
        rm2a = rma / cosl(lat)
!
! wind at time step n
        do k=1,levs
          kk=levs-k+1
          do i=1,lons_lat
!ch         uulon(i,k,lan) = ut(i,kk,lan) * rm2a
            uulon(i,k,lan) = ut_sl(i,k ,lan) * rm2a
!hmhj use real(kind=RTYPE) wind
!            vvlon(i,k,lan) = grid_gr(ilan,kvg) * rma
!hmhj use virtual wind
!            vvlon(i,k,lan) = grid_gr(ilan,kvg) / con_rerth
!for CWB GFS
!ch         vvlon(i,k,lan) = vt(i,kk,lan) * rma
            vvlon(i,k,lan) = vt_sl(i,k ,lan) * rma
          enddo
        enddo
!        if( lprint ) then
!          call mymaxmin(uulon(1,1,lan),lons_lat,lonfull,1,' uu1 in deg')
!          call mymaxmin(uulon(1,5,lan),lons_lat,lonfull,1,' uu5 in deg')
!          call mymaxmin(vvlon(1,1,lan),lons_lat,lonfull,1,' vv in deg')
!        endif
!
! u v t at n-1
!
        do k=1,levs
          ku=kuu+k-1
          kv=kvv+k-1
          kt=ktt+k-1
          kk=levs-k+1
!ch>
!       KL=lev-Llist(k)+1
!ch<
          do i=1,lons_lat
!ch         qqlon(i,ku,lan) = uum(i,kk,lan)
!ch         qqlon(i,kv,lan) = vvm(i,kk,lan)
!ch         qqlon(i,kt,lan) = ttm(i,kk,lan)
!!            qqlon(i,kp,lan) = dsigma(kk,1)*ptp(i,lan)+dsigma(kk,2)

            qqlon(i,ku,lan) = vdzonl(i,k,lan)
            qqlon(i,kv,lan) = vdmerd(i,k,lan)
            if ( nvars .ge. 3 ) qqlon(i,kt,lan) = ddtemp(i,k,lan)
!!            qqlon(i,kp,lan) = dsigma(KL,1)*ptp_sl(i,lan)+dsigma(KL,2)
          enddo
        enddo
! add surface pressure perturbation for removing resonance
!       rkt = con_g / ( con_rd * 300.0 )
!       do i=1,lons_lat
!         ilan=i+jlonf
!         qqlon(i,kuu,lan) = log(plev(i,1))+grid_gr(ilan,g_gz)*rkt
!       enddo
!
! save qqlon into n+1 for later as tendency
! z d t at n+1
!        do k=1,lev
!          kk=lev-k+1
!          ku=kuu+k-1
!          kv=kvv+k-1
!          kt=ktt+k-1
!          kug=g_d+k-1
!          kvg=g_z+k-1
!          ktg=g_t+k-1
!          do i=1,lons_lat
!            ilan=i+jlonf
!            grid_gr(ilan,kug) = qqlon(i,ku,lan)
!            grid_gr(ilan,kvg) = qqlon(i,kv,lan)
!            grid_gr(ilan,ktg) = qqlon(i,kt,lan)
!          enddo
!        enddo
! change h to theta
!
! ----- prepare xr, xcp, xkapa
!
!       xr    = cons0
!       xcp   = cons0
!       sumrq = cons0
!
!       do n=1,ntrac
!         nqq = kqq + (n-1)*levs
!         if( ri(n) .ne. cons0 .and. cpi(n) .ne. cons0 ) then
!           do k=1,lev
!             kq=nqq+k-1
!             do i=1,lons_lat
!               xr   (i,k) = xr   (i,k) + qqlon(i,kq,lan)*ri(n)
!               xcp  (i,k) = xcp  (i,k) + qqlon(i,kq,lan)*cpi(n)
!               sumrq(i,k) = sumrq(i,k) + qqlon(i,kq,lan)
!             enddo
!           enddo
!         endif
!       enddo
!       do k=1,lev
!         do i=1,lons_lat
!           xr (i,k)   = ( cons1 - sumrq(i,k) )*ri(0)  + xr (i,k)
!           xcp(i,k)   = ( cons1 - sumrq(i,k) )*cpi(0) + xcp(i,k)
!           xkappa(i,k) = xr(i,k) / xcp(i,k)
!         enddo
!       enddo
!
!       do k=1,lev
!         kt=ktt+k-1
!         ku=kuu+k-1
!         do i=1,lons_lat
!           pi = (qqlon(i,ku,lan)*0.005)**xkappa(i,k)
!           qqlon(i,kt,lan) = qqlon(i,kt,lan)/pi
!         enddo
!       enddo
!
! rq  no need for tendency
!
!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        ilan=1+jlonf
!        call mymaxmin(grid_gr(ilan,g_ttm),lons_lat,lonfull,1,' n-1 t ')
!        call mymaxmin(grid_gr(ilan,g_rrm ),lons_lat,lonfull,1,' n-1 q ')
!        call mymaxmin(grid_gr(ilan,g_zem),lons_lat,lonfull,1,' n-1 z ')
!        call mymaxmin(grid_gr(ilan,g_dim),lons_lat,lonfull,1,' n-1 d ')
!        print *,' ------------------------------------------- '
!        call mymaxmin(qqlon(1,ktt,lan),lons_lat,lonfull,1,' red t ')
!        call mymaxmin(qqlon(1,kqq,lan),lons_lat,lonfull,1,' red q ')
!        call mymaxmin(qqlon(1,kvv,lan),lons_lat,lonfull,1,' red z ')
!        call mymaxmin(qqlon(1,kuu,lan),lons_lat,lonfull,1,' red d ')
!        endif

!       call cyclic_cell_intpx(levs,lons_lat,lonf,uulon(1,1,lan))

!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        call mymaxmin(qqlon(1,ktt,lan),lonfull,lonfull,1,' full t ')
!        call mymaxmin(qqlon(1,kqq,lan),lonfull,lonfull,1,' full q ')
!        call mymaxmin(qqlon(1,kvv,lan),lonfull,lonfull,1,' full z ')
!        call mymaxmin(qqlon(1,kuu,lan),lonfull,lonfull,1,' full d ')
!        endif

!        rrlon(:,:,lan) = qqlon(:,:,lan)
!
! first set positive advection in east-west direction
!
!!        call cyclic_cell_massadvxl(lons_lat,lonfull,levs,nvars,deltim, &
!!                         uulon(1,1,lan),rrlon(1,1,lan),mass,forward)
!       call cyclic_mono_advectx (lonfull,levs,nvars,deltim,               &
!    &                   uulon(1,1,lan),rrlon(1,1,lan),mono)

!!        call cyclic_cell_intpx(nlevs,lons_lat,lonfull,rrlon(1,1,lan))
        call cyclic_cell_intpx(nlevs,lons_lat,lonfull,qqlon(1,1,lan))
        call cyclic_cell_intpx(levs,lons_lat,lonfull,vvlon(1,1,lan))

!        if( lprint ) then
!        print *,' done cyclic_massadvx with mass= ',mass
!        print *,' ------------------------------------------- '
!        call mymaxmin(rrlon(1,ktt,lan),lonfull,lonfull,1,' advx t ')
!        call mymaxmin(rrlon(1,kqq,lan),lonfull,lonfull,1,' advx q ')
!        call mymaxmin(rrlon(1,kvv,lan),lonfull,lonfull,1,' advx z ')
!        call mymaxmin(rrlon(1,kuu,lan),lonfull,lonfull,1,' advx d ')
!        print *,' done the first x adv for lan=',lan
!        print *,' =========================================== '
!        endif

      enddo
!$omp end parallel do

! ---------------------------------------------------------------------
! mpi para from east-west full grid to north-south full grid
! ---------------------------------------------------------------------
!
! para vvlon, qqlon, and rrlon to vvlat, qqlat, rrlat
!       if( lprint ) print *,' ndslfv_advect transport from we to ns '

       call para_we2ns(vvlon,vvlat,levs,my)
!!       call para_we2ns(rrlon,rrlat,nlevs,my)
       call para_we2ns(qqlon,qqlat,nlevs,my)

!       if( lprint ) then
!       print *,' ------------ after we2ns ---------------------- '
!       do lon=1,mylonlen
!        print *,'  lon=',lon
!        call mymaxmin(vvlat(1,1  ,lon),latfull,latfull,1,' we2ns v')
!        call mymaxmin(rrlat(1,kqq,lon),latfull,latfull,1,' we2ns r')
!        call mymaxmin(qqlat(1,kqq,lon),latfull,latfull,1,' we2ns q')
!       enddo
!       endif
!
! ---------------------------------------------------------------------
! -------------- in north-soutn great circle -------------------
! ---------------------------------------------------------------------

!      if ( myrank .eq. 0 ) lprint = .true.
!       if( lprint ) then
!       print *,' ndslfv_advect adv loop in y '
!       print *,' mylonlen=',mylonlen
!       endif

!$omp parallel do private(lon,k,j,ku,kv) &
!$omp schedule(dynamic)
       do lon=1,mylonlen
!
! convert wind before advy
        do k=1,levs
          ku=kuu+k-1
          kv=kvv+k-1
          do j=1,lathalf
!!            rrlat(j,ku,lon) = -rrlat(j,ku,lon)
!!            rrlat(j,kv,lon) = -rrlat(j,kv,lon)
            qqlat(j,ku,lon) = -qqlat(j,ku,lon)
            qqlat(j,kv,lon) = -qqlat(j,kv,lon)
          enddo
        enddo
!
!        if( lprint ) print *,' lon=',lon

! first set advection in north-south direction in great circle through two poles
!
!        call fixend_cell_massadvy(latfull,lathalf,levs,nvars,deltim, &
!                         vvlat(1,1,lon),rrlat(1,1,lon),mass)
!!        call cyclic_cell_massadvyl(latfull,levs,nvars,deltim,         &
!!                         vvlat(1,1,lon),rrlat(1,1,lon),mass,forward)
!       call cyclic_mono_advecty (latfull,levs,nvars,deltim,
!    &                   vvlat(1,1,lon),rrlat(1,1,lon),mono)

!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        call mymaxmin(rrlat(1,kuu,lon),latfull,latfull,1,' advry u ')
!        call mymaxmin(rrlat(1,kvv,lon),latfull,latfull,1,' advry v ')
!        call mymaxmin(rrlat(1,ktt,lon),latfull,latfull,1,' advry t ')
!        call mymaxmin(rrlat(1,kqq,lon),latfull,latfull,1,' advry q ')
!        call mymaxmin(rrlat(1,kvv,lon),latfull,latfull,1,' advry z ')
!        call mymaxmin(rrlat(1,kuu,lon),latfull,latfull,1,' advry d ')
!        endif
!
! second set advection in north-south direction in great circle through two poles
!
!        call fixend_cell_massadvy(latfull,lathalf,levs,nvars,deltim, &
!                         vvlat(1,1,lon),qqlat(1,1,lon),mass)
        call cyclic_cell_massadvyl(latfull,levs,nvars,deltim,         &
                         vvlat(1,1,lon),qqlat(1,1,lon),mass,forward)
!       call cyclic_mono_advecty (my,levs,nvars,deltim,
!    &                   vvlat(1,1,lon),qqlat(1,1,lon),mono)

!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        call mymaxmin(qqlat(1,kuu,lon),latfull,latfull,1,' advqy u ')
!        call mymaxmin(qqlat(1,kvv,lon),latfull,latfull,1,' advqy v ')
!        call mymaxmin(qqlat(1,ktt,lon),latfull,latfull,1,' advqy t ')
!        call mymaxmin(qqlat(1,kqq,lon),latfull,latfull,1,' advqy q ')
!        call mymaxmin(qqlat(1,kvv,lon),latfull,latfull,1,' advqy z ')
!        call mymaxmin(qqlat(1,kuu,lon),latfull,latfull,1,' advqy d ')
!        print *,' done with y at lon=',lon
!        endif
! convert wind back after advy
        do k=1,levs
          ku=kuu+k-1
          kv=kvv+k-1
          do j=1,lathalf
!!            rrlat(j,ku,lon) = -rrlat(j,ku,lon)
!!            rrlat(j,kv,lon) = -rrlat(j,kv,lon)
            qqlat(j,ku,lon) = -qqlat(j,ku,lon)
            qqlat(j,kv,lon) = -qqlat(j,kv,lon)
          enddo
        enddo

       enddo
!$omp end parallel do
!
! ----------------------------------------------------------------------
! mpi para from north-south direction to east-west direeectory
! ----------------------------------------------------------------------
!
! para qqlat and rrlat to qqlon and rrlon

!       if( lprint ) print *,' ndslfv_advect transport from ns to we '

!!       call para_ns2we(rrlat,rrlon,nlevs,my)
       call para_ns2we(qqlat,qqlon,nlevs,my)
!       if( lprint ) then
!       print *,' ------------ after ns2we ---------------------- '
!       do lan=1,jlistnum
!        print *,'  lan=',lan
!        call mymaxmin(rrlon(1,kqq,lan),lonfull,lonfull,1,' ns2we r')
!        call mymaxmin(qqlon(1,kqq,lan),lonfull,lonfull,1,' ns2we q')
!       enddo
!       endif

! ---------------------------------------------------------------
! ---------------- back to east-west direction ------------------
! ---------------------------------------------------------------
!      print *,' ndslfv_advect adv loop in x for last '

!$omp parallel do                                                &
!$omp private(lan,lat,lons_lat,i,k,kk,n,kv,ku,kt)    &
!$omp schedule(dynamic)

      do lan=1,jlistnum

!        lat = global_lats_a(ipt_jlistnum-1+lan)
        lat = jlist1(lan)
        lons_lat = lonsperlat(lat)
!
! mass conserving interpolation from full grid to reduced grid
!
        call cyclic_cell_intpx(nlevs,lonfull,lons_lat,qqlon(1,1,lan))
!!        call cyclic_cell_intpx(nlevs,lonfull,lons_lat,rrlon(1,1,lan))

!
! second set advection in x for the second of the pair
!
        call cyclic_cell_massadvxl(lons_lat,lonfull,levs,nvars,deltim, &
                         uulon(1,1,lan),qqlon(1,1,lan),mass,forward)
!       call cyclic_mono_advectx (lonfull,levs,nvars,deltim,               &
!    &                   uulon(1,1,lan),qqlon(1,1,lan),mono)
!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        call mymaxmin(qqlon(1,kqq,lan),lonfull,lonfull,1,' adv x q ')
!        endif

!!        do k=1,nlevs
!!          do i=1,lons_lat
!!            qqlon(i,k,lan) = 0.5 * ( qqlon(i,k,lan) + rrlon(i,k,lan) )
!!          enddo
!!        enddo

!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        call mymaxmin(qqlon(1,kuu,lan),lons_lat,lonfull,1,' do redu u ')
!        call mymaxmin(qqlon(1,kvv,lan),lons_lat,lonfull,1,' do redu v ')
!        call mymaxmin(qqlon(1,ktt,lan),lons_lat,lonfull,1,' do redu t ')
!        call mymaxmin(qqlon(1,kqq,lan),lons_lat,lonfull,1,' do redu q ')
!        call mymaxmin(qqlon(1,kvv,lan),lons_lat,lonfull,1,' do redu z')
!        call mymaxmin(qqlon(1,kuu,lan),lons_lat,lonfull,1,' do redu d')
!        print *,' finish horizonatal advection at lan=',lan
!        endif

! change theta to h
!
! ----- prepare xr, xcp, xkapa
!
!       xr    = cons0
!       xcp   = cons0
!       sumrq = cons0
!
!       do n=1,ntrac
!         nqq = kqq + (n-1)*levs
!         if( ri(n) .ne. cons0 .and. cpi(n) .ne. cons0 ) then
!           do k=1,lev
!             kq=nqq+k-1
!             do i=1,lons_lat
!               xr   (i,k) = xr   (i,k) + qqlon(i,kq,lan)*ri(n)
!               xcp  (i,k) = xcp  (i,k) + qqlon(i,kq,lan)*cpi(n)
!               sumrq(i,k) = sumrq(i,k) + qqlon(i,kq,lan)
!             enddo
!           enddo
!         endif
!       enddo
!       do k=1,lev
!         do i=1,lons_lat
!           xr (i,k)   = ( cons1 - sumrq(i,k) )*ri(0)  + xr (i,k)
!           xcp(i,k)   = ( cons1 - sumrq(i,k) )*cpi(0) + xcp(i,k)
!           xkappa(i,k) = xr(i,k) / xcp(i,k)
!         enddo
!       enddo
!
!       do k=1,lev
!         kt=ktt+k-1
!         ku=kuu+k-1
!         do i=1,lons_lat
!           pi = (qqlon(i,ku,lan)*0.005)**xkappa(i,k)
!           qqlon(i,kt,lan) = qqlon(i,kt,lan)*pi
!         enddo
!       enddo
! u v t update at n+1
        do k=1,levs
          kk=levs-k+1
          ku=kuu+k-1
          kv=kvv+k-1
          kt=ktt+k-1
!!          kp=kup+k-1
!ch>
!       KL=lev-Llist(k)+1:w
!ch<
          do i=1,lons_lat
!ch         vdzonl(i,kk,lan) = (qqlon(i,ku,lan)-uum(i,kk,lan))*rdt2
!ch         vdmerd(i,kk,lan) = (qqlon(i,kv,lan)-vvm(i,kk,lan))*rdt2
!ch         ddtemp(i,kk,lan) = qqlon(i,kt,lan)
!ch         pten(i,kk,lan)    =(qqlon(i,kp,lan)-(dsigma(kk,1)          &
!ch                            *ptp(i,lan)+dsigma(kk,2)))*rdt2

            vdzonl(i,k ,lan) = qqlon(i,ku,lan)
            vdmerd(i,k ,lan) = qqlon(i,kv,lan)
            if ( nvars .ge. 3 ) ddtemp(i,k ,lan) = qqlon(i,kt,lan)
!            pten(i,k ,lan)    =(qqlon(i,kp,lan)-(dsigma(KL,1)          &
!                               *ptp_sl(i,lan)+dsigma(KL,2)))*rdt2
          enddo
        enddo
!        if (myrank .eq. 0 ) then
!        call mymaxmin(qqlon(1,ktt,lan),lons_lat,lonfull,1,' tend temp')
!        endif
!        if( lprint ) then
!        print *,' ------------------------------------------- '
!        ilan=1+jlonf
!        call mymaxmin(grid_gr(ilan,g_u),lons_lat,lonfull,1,' tend u ')
!        call mymaxmin(grid_gr(ilan,g_v),lons_lat,lonfull,1,' tend v ')
!        call mymaxmin(grid_gr(ilan,g_tt),lons_lat,lonfull,1,' tend t ')
!        call mymaxmin(grid_gr(ilan,g_rt),lons_lat,lonfull,1,' tend q ')
!        call mymaxmin(grid_gr(ilan,g_z),lons_lat,lonfull,1,'tend z')
!        do k=1,lev,10
!        kug=g_d+k-1
!        print *,' k=',k
!        call mymaxmin(grid_gr(ilan,kug),lons_lat,lonfull,1,' tend d ')
!        enddo
!        print *,' finish horizonatal advection at lan=',lan
!        endif
!
      enddo
!$omp end parallel do

!
! ===============================
!
      return
      end
