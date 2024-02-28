      subroutine ndslfv_monoadvv (ddtemp,qvadv,vdzonl,vdmerd,pdot      &
                                , pt,lonsperlat,deltim,forward)
!
! a routine to do non-iteration semi-Lagrangain advection
! considering advection  with monotonicity in interpolation
! contact: hann-ming henry juang
! program log:
! 2011 02 20 : henry juang, initial implemented into nems as NDSL with mass_dp
! 2013 09 30 : henry juang, add option of theta advection, (used later)
!
!
      use param
      use grid, only : latpart,ndslvvar
      use index
      use rank
      use const

      implicit none

!ch   real(kind=RTYPE) pdot(lonfull,lev+1,latpart)
!ch   real(kind=RTYPE) plev(lonfull,lev+1)
      real(kind=RTYPE) pdot(nxp,    lev+1,latpart)
      real(kind=RTYPE) plev(nxp,    lev+1)
      real(kind=RTYPE) pt(nxp,latpart)
!      integer,intent(in):: global_lats_a(my)
      integer,intent(in):: lonsperlat(my)
      real(kind=RTYPE),   intent(in):: deltim

!ch   real(kind=RTYPE)      qqlon(lonfull,lev*ndslvvar,latpart)
      real(kind=RTYPE)      qqlon(nxp,    lev*ndslvvar,latpart)
      real(kind=RTYPE)      ddtemp(nxp,lev,my_max),qvadv(nxp,lev*ncld,my_max),    &
                vdmerd(nxp,lev,my_max),vdzonl(nxp,lev,my_max)
      real(kind=RTYPE)      rdt2,dt2

!     real(kind=RTYPE)      xksav(lonfull,levs      ,latpart)
!     real(kind=RTYPE)      stsav(lonfull,levs      ,latpart)
!     real(kind=RTYPE)      ttsav(lonfull,levs      ,latpart)

!     real(kind=RTYPE)      xr    (lonfull,levs)
!     real(kind=RTYPE)      xcp   (lonfull,levs)
!     real(kind=RTYPE)      sumrq (lonfull,levs)
!     real(kind=RTYPE)      xkappa(lonfull,levs)
!     real(kind=RTYPE)      kappa, pi, ply, hh
!      real(kind=RTYPE)      cons0, cons1

!      logical 	lprint

      integer mono,mass,nvars
      integer ii,i,n,k,kk,lon,lan,lat,lons_lat
      integer kqq, ktt, kuu, kvv
      integer kq , kt , ku , kv
      logical forward
!
!     lprint = .false.

!     if( lprint ) then
!       print *,' enter ndslfv_advect  with monotonicity '
!     endif
!
      mono  = 1
      mass  = 0
!      cons0 = 0.0
!      cons1 = 1.0
!
!      levh = ncld * lev
!
      kuu = 1
      kvv = kuu + lev
      ktt = kvv + lev
      kqq = ktt + lev

      nvars = ndslvvar
!
      rdt2 = 0.5 / deltim
      dt2 =  2. * deltim
!
!$omp parallel do                                                &
!$omp private(lan,lat,lons_lat,plev,i,k,n,kk,mass,ku,kv,kt,kq)   &
!$omp schedule(dynamic)

!      do lan=1,lats_node_a
      do lan=1,jlistnum

!        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lat = jlist1(lan)
        lons_lat = lonsperlat(lat)

        plev(:,lev+1) = 0.0
        do k=lev,1,-1
          kk=lev-k+1
          do i=1,lons_lat
            plev(i,k)=plev(i,k+1)+dsigma(kk,1)*pt(i,lan)+dsigma(kk,2)
          enddo
        enddo

!       if( lprint ) then
!       do k=1,lev
!       print *,' k= ',k
!       call mymaxmin(pdot(1,k,lan),lons_lat,lonfull,1,' pdot ')
!       call mymaxmin(plev(1,k),lons_lat,lonfull,1,' plev ')
!       enddo
!       endif
!
! d z t at n+1*
!       kappa = con_rd / con_cp
        do k=1,lev
          kk=lev-k+1
          ku=kuu+k-1
          kv=kvv+k-1
          kt=ktt+k-1
          do i=1,lons_lat
            qqlon(i,ku,lan) = vdzonl(i,kk,lan)
            qqlon(i,kv,lan) = vdmerd(i,kk,lan)
            qqlon(i,kt,lan) = ddtemp(i,kk,lan)
          enddo
        enddo
! rq at n+1*
        do n=1,ncld
        do k=1,lev
          kk=lev-k+1+(n-1)*lev
          kq=kqq+k-1+(n-1)*lev
          do i=1,lons_lat
            qqlon(i,kq,lan) = qvadv(i,kk,lan)
          enddo
        enddo
        enddo
!
! ----- prepare xr, xcp, xkapa
!
!       xr    = cons0
!       xcp   = cons0
!       sumrq = cons0
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
!           xksav(i,k,lan) = xkappa(i,k)
!         enddo
!       enddo
!
! change h to theta
!
!       do k=1,levs
!         kt=ktt+k-1
!         do i=1,lons_lat
!           pi = ((plev(i,k)+plev(i,k+1))*0.005)**xkappa(i,k)
!           ttsav(i,k ,lan) = qqlon(i,kt,lan)
!           qqlon(i,kt,lan) = qqlon(i,kt,lan)/pi
!           stsav(i,k ,lan) = qqlon(i,kt,lan)
!         enddo
!       enddo
!
!
        mass=0
!ch     call vertical_cell_advect (lons_lat,lonfull,lev,nvars, &
        call vertical_cell_advect (lons_lat,nxp,    lev,nvars, &
                  deltim,plev,pdot(1,1,lan),qqlon(1,1,lan),mass,forward)
!
! dp with mass conserving
!       do k=1,levs
!         kp =g_dp +k-1
!         kpg=g_dpn+k-1
!         do i=1,lon_dim
!           ilan=i+jlonf
!           rrlon(i,k,lan) = grid_gr(ilan,kpg)/grid_gr(ilan,kp )
!           rrlon(i,k,lan) = 1.0
!         enddo
!       enddo
!       mass=1
!       call vertical_cell_advect (lons_lat,lonfull,levs,1,
!    &            deltim,plev,pdot(1,1,lan),rrlon(1,1,lan),mass)
!       do k=1,levs
!         kp =g_dp +k-1
!         kpg=g_dpn+k-1
!         do i=1,lon_dim
!           ilan=i+jlonf
!           grid_gr(ilan,kpg) = rrlon(i,k ,lan)*grid_gr(ilan,kp )
!         enddo
!       enddo
!
! ----- prepare xr, xcp, xkapa
!
!       xr    = cons0
!       xcp   = cons0
!       sumrq = cons0
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
! add change of p, theta and kappa to h
!
!       do k=1,levs
!         kt=ktt+k-1
!         do i=1,lons_lat
!           ply = (plev(i,k)+plev(i,k+1))*0.5
!           pi = ((plev(i,k)+plev(i,k+1))*0.005)**xkappa(i,k)
!           hh = ttsav(i,k,lan)
!           ttsav(i,k,lan) = ttsav(i,k,lan)
!    &                      + xksav(i,k,lan)*hh/ply*
!    &                       (pdot(i,k,lan)+pdot(i,k+1,lan))*deltim
!           ttsav(i,k,lan) = ttsav(i,k,lan)
!    &                      + pi*(qqlon(i,kt,lan)-stsav(i,k,lan))
!    &                      + hh*log(ply/100.)*
!    &                        (xkappa(i,k)-xksav(i,k,lan))
!           qqlon(i,kt,lan)= ttsav(i,k,lan)
!         enddo
!       enddo
!
! u v t tendency at n
        do k=1,lev
          kk=lev-k+1
          ku=kuu+k-1
          kv=kvv+k-1
          kt=ktt+k-1
          do i=1,lons_lat
!!            vdzonl(i,kk,lan) = (qqlon(i,ku,lan)-up(i,kk,lan))*rdt2
!!            vdmerd(i,kk,lan) = (vp(i,kk,lan)-qqlon(i,kv,lan))*rdt2
!!            ddtemp(i,kk,lan) = (qqlon(i,kt,lan)-ttp(i,kk,lan))*rdt2
            vdzonl(i,kk,lan) = qqlon(i,ku,lan)
            vdmerd(i,kk,lan) = qqlon(i,kv,lan)
            ddtemp(i,kk,lan) = qqlon(i,kt,lan)
          enddo
        enddo
! rq tendency at n
        do n=1,ncld
        do k=1,lev
          kk=lev-k+1+(n-1)*lev
          kq=kqq+k-1+(n-1)*lev
          do i=1,lons_lat
!byl no need tendency for Tracers
!!            qvadv(i,kk,lan) = (qqlon(i,kq,lan)-qm(i,kk,lan))*rdt2
            qvadv(i,kk,lan) = qqlon(i,kq,lan)
          enddo
        enddo
        enddo
!
!       if( lprint ) then
!       call mymaxmin(qqlon(1,kqq,lan),lons_lat,lonfull,1,' q vertadv')
!       print *,' ------------------------------------------- '
!       print *,' finish updating n+1* in grid_gr at lan=',lan
!       endif
!
      enddo
!$omp end parallel do

!
! ===============================
!
      return
      end
!
! -------------------------------
      subroutine ndslfv_monoadvv_fgnl(vdzonl,vdmerd,ddtemp,pdot      &
                                , pt,lonsperlat,deltim,nvars,forward)
!
! a routine to do non-iteration semi-Lagrangain advection
! considering advection  with monotonicity in interpolation
! contact: hann-ming henry juang
! program log:
! 2011 02 20 : henry juang, initial implemented into nems as NDSL with mass_dp
! 2013 09 30 : henry juang, add option of theta advection, (used later)
!
!
      use param
      use grid, only : latpart
      use index
      use rank
      use const

      implicit none

      real(kind=RTYPE) pdot(nxp,    lev+1,latpart)
      real(kind=RTYPE) plev(nxp,    lev+1)
      real(kind=RTYPE) pt(nxp,latpart)
      integer,intent(in):: lonsperlat(my)
      real(kind=RTYPE),   intent(in):: deltim

      real(kind=RTYPE)      qqlon(nxp,    lev*nvars,latpart)
      real(kind=RTYPE)      vdmerd(nxp,lev,my_max),vdzonl(nxp,lev,my_max)
      real(kind=RTYPE)      ddtemp(nxp,lev,my_max)
      real(kind=RTYPE)      rdt2,dt2

      integer mono,mass,nvars
      integer ii,i,n,k,kk,lon,lan,lat,lons_lat
      integer kuu, kvv, ktt
      integer ku , kv, kt
      logical forward
!
!     lprint = .false.

!     if( lprint ) then
!       print *,' enter ndslfv_advect  with monotonicity '
!     endif
!
      mono  = 1
      mass  = 0
!      cons0 = 0.0
!      cons1 = 1.0
!
!      levh = ncld * lev
!
      kuu = 1
      kvv = kuu + lev
      ktt = kvv + lev

!      nvars = 2
!
      rdt2 = 0.5 / deltim
      dt2 =  2. * deltim
!
!$omp parallel do                                                &
!$omp private(lan,lat,lons_lat,plev,i,k,n,kk,mass,ku,kv,kt)      &
!$omp schedule(dynamic)

!      do lan=1,lats_node_a
      do lan=1,jlistnum

!        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lat = jlist1(lan)
        lons_lat = lonsperlat(lat)

        plev(:,lev+1) = 0.0
        do k=lev,1,-1
          kk=lev-k+1
          do i=1,lons_lat
            plev(i,k)=plev(i,k+1)+dsigma(kk,1)*pt(i,lan)+dsigma(kk,2)
          enddo
        enddo

!       if( lprint ) then
!       do k=1,lev
!       print *,' k= ',k
!       call mymaxmin(pdot(1,k,lan),lons_lat,lonfull,1,' pdot ')
!       call mymaxmin(plev(1,k),lons_lat,lonfull,1,' plev ')
!       enddo
!       endif
!
! d z t at n+1*
!       kappa = con_rd / con_cp
        do k=1,lev
          kk=lev-k+1
          ku=kuu+k-1
          kv=kvv+k-1
          kt=ktt+k-1
          do i=1,lons_lat
            qqlon(i,ku,lan) = vdzonl(i,kk,lan)
            qqlon(i,kv,lan) = vdmerd(i,kk,lan)
            if ( nvars .ge. 3 ) qqlon(i,kt,lan) = ddtemp(i,kk,lan)
          enddo
        enddo
!
!
        mass=0
!
        call vertical_cell_advect (lons_lat,nxp,    lev,nvars, &
                  deltim,plev,pdot(1,1,lan),qqlon(1,1,lan),mass,forward)
!
!
! u v t tendency at n
        do k=1,lev
          kk=lev-k+1
          ku=kuu+k-1
          kv=kvv+k-1
          kt=ktt+k-1
          do i=1,lons_lat
            vdzonl(i,kk,lan) = qqlon(i,ku,lan)
            vdmerd(i,kk,lan) = qqlon(i,kv,lan)
            if ( nvars .ge. 3 ) ddtemp(i,kk,lan) = qqlon(i,kt,lan)
          enddo
        enddo
!
!       if( lprint ) then
!       call mymaxmin(qqlon(1,kqq,lan),lons_lat,lonfull,1,' q vertadv')
!       print *,' ------------------------------------------- '
!       print *,' finish updating n+1* in grid_gr at lan=',lan
!       endif
!
      enddo
!$omp end parallel do

!
! ===============================
!
      return
      end
