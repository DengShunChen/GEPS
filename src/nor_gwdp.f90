      subroutine nor_gwdp(j,nxj,nx,lev,u,v,t,q,plt,pk,pk2,hi,dt,grav,rgas,sinl &
      ,                   cosl,avgdrag_u,avgdrag_v,cp,ptop)
!
! reference : Scinocca(2002,2003) JAS and Andrew Orr(2010) J. of Climate
!
!---- input:
!---- nxj    : number of point in lontitude of reduced gauslian grid.
!---- nx     : dimension of lontitude.
!---- lev    : dimension of vertical.
!---- dt     : 2.*dt dt: time step of model integration
!---- grav   : gravity of earth
!     rgas   : dry gas constant (=287.)
!---- sinl   : sin of latitude
!     t     : temperature  (k)
!     q     : specific humidity (kg/kg)
!     plt   : pressure at odd-level (mb)
!     pk    : p**kcapa at odd-level ((=p/1000.)**0.286)
!     pk2   : p**kcapa at even-level ((=p/1000.)**0.286)
!     hi    : geopotential height (m2/s2)
!---- u      : velocity of x-component
!---- v      : velocity of y-component
!---- output:
!---  u,v
!-----------------------------------------
!---- den    : density of atmosphere
!---- sn2    : Brunt-Vasisla frequency
!---- p2     : pressure at even-level (mb)
!---- nc        : number of spectral elements of the gravity wave phase speed.
!---- nphi      : number of azimuths of the gravity waves.
!---- cmin      : minimun phase speed of the gravity waves(m/s).
!---- cmax      : the maximun phase speed of the gravity waves(m/s).
!---- taul      : constant=0.6 s/m
!---- cstar     : constant=1.
!---- pw        : constant=1.5 exponent power of gravity wave frequency
!---- s         : constant=1. low-m exponent
!---- lunch_lev     : lunch level of gravity wave
!---- wavelenth : characteric vertical wave-length
!---- fluxtoatl : total lunch E-P flux in each nc direction
!
!  modify to f90 by C-H Lee and sort by River Chen in 2015
!

!ch>
      use index, only : nxdef
      use const, only : RTYPE
!ch<

      implicit  none

      integer   j,nxj,nx,lev
      real      dt,grav,rgas,sinl,cosl,cp

      real      plt(nx,lev),                                  &
                avgdrag_u(lev),avgdrag_v(lev)
      real(kind=RTYPE) u(nx,lev),v(nx,lev),t(nx,lev),         &
                       q(nx,lev),hi(nx,lev),pk(nx,lev),       &
                       pk2(nx,lev)

! local array
      integer,  parameter :: nc=20,nphi=4
      real      p2(nx,lev),den(nx,lev),sn2(nx,lev),pt(nx,lev)
      real      fluxsat(nc,nx,lev,nphi),flux(nc,nx,lev,nphi),dflux(nx,lev,nphi),&
                fluxsat_nc(nx,lev,nphi),flux_nc(nx,lev,nphi),           &
                dflux_u(nx,lev),dflux_v(nx,lev),                        &
                a(nx),cn(nc),uv_lunch(nx,nphi),uv(nx,lev,nphi),         &
                fxin(nc,nphi),fxout(nc),drag_u(nx,lev),drag_v(nx,lev)
      real      mstar
!
      integer   ilon,jlat,kk,lunch_lev,ltop,i,k,nc1,n,nn,ndir
      real      cmin,cmax,taul,cstar,pw,s,wavelenth,fluxtotal,ptop
      real      pi,omega,acor,cor,wtcosl,dz1,s1,s2,p,corp,xmin,xmax
      real      b1,b2,dc,xbar,x,s3,uhat,fxtotal,chat,flux1,flux2,dflux1
      real      uhat_z0,uhat_z1,dpk,wsold,wsnew

!      data cmin/0.25/, cmax/2000./ taul/0.6/
!      data cstar/1./, pw/1./, s/1./, wavelenth/2000./
!      data fluxtotal/3.75e-4/, lunch_lev/44/, ltop/1/
      data cmin/0.25/, cmax/1000./ taul/0.25/, lunch_lev/44/
      data cstar/1./, ltop/1/, wavelenth/2000./
      data fluxtotal/3.75e-4/, pw/1.5/, s/0./    !CG1
!
      ilon=1142
      jlat=nx/4
      kk=31
!
      pi=4.*atan(1.)
      omega=7.292e-5
      cor=2.*omega*sinl
      acor=abs(cor)
!
!ch   wtcosl=nxj*cosl
      wtcosl=nxdef(j)*cosl
!
! ... pressure at even-level (p2)
! ... potential temperature (pt)
! ... density (den)
!
      p2=0.
      pt=0.
      den=0.
      sn2=0.

      do k=1,lunch_lev+1
      do i=1,nxj
        p2(i,k)=1000.0*pk2(i,k)*pk2(i,k)*pk2(i,k)*sqrt(pk2(i,k))
        pt(i,k)=t(i,k)*(1.0+0.608*q(i,k))/pk(i,k)
        den(i,k)=plt(i,k)*100./(rgas*t(i,k))
      enddo
      enddo
!
! ... brunt-vaisala frequency (sn2)
!
      do k=1,lunch_lev
      do i=1,nxj
        dz1=(hi(i,k)-hi(i,k+1))/grav
        s1=pt(i,k)-pt(i,k+1)
        s2=(pt(i,k)+pt(i,k+1))/2.
        sn2(i,k)=grav*s1/(s2*dz1)
        if(sn2(i,k) .le. acor)sn2(i,k)=acor
      enddo
      enddo
!      do i=1,nxj
!        sn2(i,lev)=sn2(i,lev-1)
!      enddo
!      if(j.eq.jlat)then
!          call findmaxf(nxj,nx,lev,sn2,vmax,iv,kv,vmin,ivm,kvm)
!          print *,' sn2=',vmax,iv,kv,' vmin=',vmin,ivm,kvm
!          call findmaxf(nxj,nx,lev,den,vmax,iv,kv,vmin,ivm,kvm)
!          print *,' den=',vmax,iv,kv,' vmin=',vmin,ivm,kvm
!      endif
!
!  coefficient of A(indepedent of height) in eq-(4) (Orr,2010)
!
      mstar=2.*pi/wavelenth
      a=0.
      p=2.-pw
      corp=acor**p
!      if(j.eq.jlat)print *,'mstar,p,cor,corp=',mstar,p,cor,corp
!
      do i=1,nxj
        a(i)=cstar*mstar**3*(sn2(i,lunch_lev)**p-corp)/p
        if(a(i).le.0.)a(i)=0.
      enddo
!      if(j.eq.jlat)then
!          call findmaxf(nxj,nx,1,a,vmax,iv,kv,vmin,ivm,kvm)
!          print *,' a=',vmax,iv,kv,' vmin=',vmin,ivm,kvm
!      endif
!
      uv_lunch=0.
      uv=0.
      k=lunch_lev
      do i=1,nxj
        uv_lunch(i,1)=u(i,k)
        uv_lunch(i,2)=u(i,k)
        uv_lunch(i,3)=v(i,k)
        uv_lunch(i,4)=v(i,k)
        uv(i,k,1)=uv_lunch(i,1)
        uv(i,k,2)=uv_lunch(i,2)
        uv(i,k,3)=uv_lunch(i,3)
        uv(i,k,4)=uv_lunch(i,4)
      enddo
!      if(j.eq.jlat)then
!        do nn=1,nphi
!          print *,' lunch nn,u,v=',nn,uv_lunch(ilon,nn)
!        enddo
!      endif
!
! do coordinate stretch, eq-(28),(29),(30) (Scinocca,2003)
! cn: transformed wave phase speed
!
      xmin=1./cmin
      xmax=1./cmax
      b1=(xmax-xmin)/(exp((xmax-xmin)/taul)-1)
      b2=xmin-b1
!
      dc=xmax-xmin
      nc1=nc-1
      do n=1,nc
        xbar=xmin+dc*(n-1)/nc1
        x=b1*exp((xbar-xmin)/taul)+b2
        cn(n)=1./x
      enddo
!
!      if(j.eq.jlat)then
!          print *,' b1,b2=',b1,b2
!          print *,' cn=',(cn(n),n=1,nc)
!      endif
!
      s3=s+3
      fluxsat=0.
      flux=0.
      dflux=0.
      fluxsat_nc=0.
      flux_nc=0.
!
!  for lunch level, assign fluxtotal value to spectral flux density and
!  saturate flux density and assign u,v at lunch  level to variable 'uv'
!  for pw=1,p=2-pw=1 and uhat=0, the formula of flux density is simply rewriten.
!
      k=lunch_lev
      uhat=0.
      do nn=1,nphi
        ndir=1.
        if(nn.eq.2 .or. nn.eq.4)ndir=-1.
        fxtotal=fluxtotal*ndir
        do i=1,nxj
          do n=1,nc
            chat=cn(n)*ndir
            flux1=cstar*den(i,k)*a(i)*(chat-uhat)/sn2(i,k) &
!    &           *abs((chat-uhat)/chat)**p
                 *abs((chat-uhat)/chat)
            flux2=1./(1.+(mstar*(chat-uhat)/sn2(i,k))**s3)
            fluxsat(n,i,k,nn)=flux1
            fxin(n,nn)=flux1*flux2
          enddo
!
!          if(j.eq.jlat .and. i.eq.ilon)then
!          if(nn.eq.1)then
!            print *,'bef norm fxin at lunch=',(fxin(n,nn),n=1,nc)
!          endif
!          endif
!
          call normalize(fxtotal,nc,fxin(1,nn),flux(1,i,k,nn))
!
!          if(j.eq.jlat .and. i.eq.ilon)then
!          if(nn.eq.1)then
!            print *,'aft norm fxin at lunch=',(flux(n,i,k,nn),n=1,nc)
!          endif
!          endif
!
          do n=1,nc
            fluxsat_nc(i,k,nn)=fluxsat_nc(i,k,nn)+fluxsat(n,i,k,nn)
            flux_nc(i,k,nn)=flux_nc(i,k,nn)+flux(n,i,k,nn)
          enddo
        enddo
      enddo
!
!CWB2015 turn off for reproducible
!!$omp  parallel do default(none) &
!!$omp  private(nn,ndir,k,i,uv,uhat_z1,uhat_z0,n,chat,flux1,flux2,   &
!!$omp  fxout,dflux1,fxtotal) &
!!$omp  shared(ltop,lunch_lev,nxj,p,fxin, &
!!$omp  u,v,uv_lunch,cn,cstar,den,a,sn2,s3,mstar,fluxsat,flux,dflux, &
!!$omp  fluxsat_nc,flux_nc)
      do nn=1,nphi
        ndir=1.
        if(nn.eq.2 .or. nn.eq.4)ndir=-1.
        do k=lunch_lev-1,ltop,-1
          do i=1,nxj
            uv(i,k,1)=u(i,k)
            uv(i,k,2)=u(i,k)
            uv(i,k,3)=v(i,k)
            uv(i,k,4)=v(i,k)
!
! compute spectral flux(eq-4) and saturated spectral flux(eq-6) Orr(2010)
! uhat_z1: upper level, uhat_z0: below level.
! for pw=1,p=2-pw=1, the formula of flux density is simply rewriten.
!
            uhat_z1=uv(i,k,nn)-uv_lunch(i,nn)
            uhat_z0=uv(i,k+1,nn)-uv_lunch(i,nn)
            do n=1,nc
              chat=cn(n)*ndir
              flux1=cstar*den(i,k)*a(i)*(chat-uhat_z1)/sn2(i,k) &
!    &             *abs((chat-uhat_z1)/chat)**p
                   *abs((chat-uhat_z1)/chat)
              flux2=1./(1.+(mstar*(chat-uhat_z1)/sn2(i,k))**s3)
              fluxsat(n,i,k,nn)=flux1
              fxin(n,nn)=flux1*flux2
            enddo
!
            fxtotal=flux_nc(i,k+1,nn)
            call normalize(fxtotal,nc,fxin(1,nn),flux(1,i,k,nn))
!
            do n=1,nc
!
! encounter critical level and do critical level filtering
!
              chat=cn(n)*ndir
              if(chat.le.uhat_z1 .and. chat.ge.uhat_z0)then
                dflux1=flux(n,i,k,nn)-flux(n,i,k+1,nn)
                dflux(i,k,nn)=dflux(i,k,nn)+dflux1
                fluxsat(n,i,k,nn)=fluxsat(n,i,k,nn)-dflux1
                flux(n,i,k,nn)=flux(n,i,k,nn)-dflux1
              else
                flux(n,i,k,nn)=flux(n,i,k+1,nn)
              endif
            enddo
          enddo
        enddo
        if(nn.eq.1 .or. nn.eq.3)then
        do k=lunch_lev-1,ltop,-1
          do i=1,nxj
            do n=1,nc
!
! nonlinear dispersion mechanism (saturation theorey)
!
              if( flux(n,i,k,nn) .gt. fluxsat(n,i,k,nn) )then
                dflux1=flux(n,i,k,nn)-fluxsat(n,i,k,nn)
                dflux(i,k,nn)=dflux(i,k,nn)+dflux1
                flux(n,i,k,nn)=fluxsat(n,i,k,nn)
              endif
            enddo
          enddo
        enddo
        else
        do k=lunch_lev-1,ltop,-1
          do i=1,nxj
            do n=1,nc
!
! nonlinear dispersion mechanism (saturation theorey)
!
              if( flux(n,i,k,nn) .lt. fluxsat(n,i,k,nn) )then
                dflux1=flux(n,i,k,nn)-fluxsat(n,i,k,nn)
                dflux(i,k,nn)=dflux(i,k,nn)+dflux1
                flux(n,i,k,nn)=fluxsat(n,i,k,nn)
              endif
            enddo
          enddo
        enddo
        endif
!
        do k=lunch_lev-1,ltop,-1
          do i=1,nxj
            do n=1,nc
              fluxsat_nc(i,k,nn)=fluxsat_nc(i,k,nn)+fluxsat(n,i,k,nn)
              flux_nc(i,k,nn)=flux_nc(i,k,nn)+flux(n,i,k,nn)
            enddo
          enddo
        enddo
      enddo
!!$omp end parallel do
!
      dflux_u=0.
      dflux_v=0.
      do k=lunch_lev-1,ltop,-1
        do i=1,nxj
          dflux_u(i,k)=dflux(i,k,1)+dflux(i,k,2)
          dflux_v(i,k)=dflux(i,k,3)+dflux(i,k,4)
        enddo
      enddo
!      if(j.eq.jlat)then
!          call findmaxf(nxj,nx,lev,dflux_u,vmax,iv,kv,vmin,ivm,kvm)
!          print *,' dflux_u=',vmax,iv,kv,' vmin=',vmin,ivm,kvm
!          call findmaxf(nxj,nx,lev,dflux_v,vmax,iv,kv,vmin,ivm,kvm)
!          print *,' dflux_v=',vmax,iv,kv,' vmin=',vmin,ivm,kvm
!      endif
!
!      if(j.eq.jlat)then
!        print *,' orig u,v'
!        call findmaxf(nxj,nx,lev,u,umax,iu,ku,umin,ium,kum)
!        call findmaxf(nxj,nx,lev,v,vmax,iv,kv,vmin,ivm,kvm)
!        print *,' max_u=',umax,iu,ku,' umin=',umin,ium,kum
!        print *,' max_v=',vmax,iv,kv,' vmin=',vmin,ivm,kvm
!      endif
!
      avgdrag_u=0.
      avgdrag_v=0.
      drag_u=0.
      drag_v=0.
      do k=ltop+1,lunch_lev-1
      do i=1,nxj
        dpk=(p2(i,k-1)-p2(i,k))*100.
        drag_u(i,k)=-grav*(dflux_u(i,k-1)-dflux_u(i,k))/dpk
        drag_v(i,k)=-grav*(dflux_v(i,k-1)-dflux_v(i,k))/dpk
        wsold=0.5*(u(i,k)*u(i,k)+v(i,k)*v(i,k))
!        wsold=u(i,k)*u(i,k)+v(i,k)*v(i,k)
        u(i,k)=u(i,k)+dt*drag_u(i,k)
        v(i,k)=v(i,k)+dt*drag_v(i,k)
!        wsnew=u(i,k)*u(i,k)+v(i,k)*v(i,k)
!        t(i,k)=t(i,k)+dt*(wsold-wsnew)/(4.*cp*dt)
        wsnew=0.5*(u(i,k)*u(i,k)+v(i,k)*v(i,k))
        t(i,k)=t(i,k)+(wsold-wsnew)/cp
        avgdrag_u(k)=avgdrag_u(k)+drag_u(i,k)*cosl
        avgdrag_v(k)=avgdrag_v(k)+drag_v(i,k)*cosl
      enddo
        avgdrag_u(k)=avgdrag_u(k)/wtcosl
        avgdrag_v(k)=avgdrag_v(k)/wtcosl
      enddo
      k=1
      do i=1,nxj
        wsold=0.5*(u(i,k)*u(i,k)+v(i,k)*v(i,k))
!        wsold=u(i,k)*u(i,k)+v(i,k)*v(i,k)
        dpk=(ptop-p2(i,k))*100.
        drag_u(i,k)=grav*dflux_u(i,k)/dpk
        drag_v(i,k)=grav*dflux_v(i,k)/dpk
        u(i,k)=u(i,k)+dt*drag_u(i,k)
        v(i,k)=v(i,k)+dt*drag_v(i,k)
!        u(i,k)=u(i,k)+dt*drag_u(i,k+1)*2.
!        v(i,k)=v(i,k)+dt*drag_v(i,k+1)*2.
!        wsnew=u(i,k)*u(i,k)+v(i,k)*v(i,k)
!        t(i,k)=t(i,k)+dt*(wsold-wsnew)/(4.*cp*dt)
        wsnew=0.5*(u(i,k)*u(i,k)+v(i,k)*v(i,k))
        t(i,k)=t(i,k)+(wsold-wsnew)/cp
        avgdrag_u(k)=avgdrag_u(k)+drag_u(i,k+1)*2.*cosl
        avgdrag_v(k)=avgdrag_v(k)+drag_v(i,k+1)*2.*cosl
      enddo
      avgdrag_u(k)=avgdrag_u(k)/wtcosl
      avgdrag_v(k)=avgdrag_v(k)/wtcosl
!
!      if(j.eq.jlat)then
!          print *,' drag u,v'
!          call findmaxf(nxj,nx,lev,u,umax,iu,ku,umin,ium,kum)
!          call findmaxf(nxj,nx,lev,v,vmax,iv,kv,vmin,ivm,kvm)
!          print *,' max_u=',umax,iu,ku,' umin=',umin,ium,kum
!          print *,' max_v=',vmax,iv,kv,' vmin=',vmin,ivm,kvm
!      endif
!
      return
      end
!
      subroutine findmaxf(nxj,nx,lev,f,fmax,i1,k1,fmin,i1m,k1m)
      implicit  none
      integer   nxj,nx,lev,i1,k1,i1m,k1m
      real      f(nx,lev),fmax,fmin
      integer   k,i
      
      fmax=-99999.
      do k=1,lev
        do i=1,nxj
          if(f(i,k).ge.fmax)then
            fmax=f(i,k)
            i1=i
            k1=k
          endif
        enddo
      enddo
      fmin=99999.
      do k=1,lev
        do i=1,nxj
          if(f(i,k).le.fmin)then
            fmin=f(i,k)
            i1m=i
            k1m=k
          endif
        enddo
      enddo
      return
      end
