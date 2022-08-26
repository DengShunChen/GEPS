!      subroutine cumastr_driv(nx,lev,dt,g,r,cp,hltm,ptop,land,topo     &
      subroutine cumastr_driv(nxj,nx,lev,dt,g,r,cp,hltm,ptop,land,topo     &
                           , phi,u,v,t,q,ut,vt,tt,qt,rcup,pk,pk2,sd   &
                           , qflux,kcbot,kctop,ncld,sigma,plt,pt,j    &
!xb110>
                           ,kcnv)
!xb110<
!c
!c#######################################################################
!c                     subroutine description
!c
!c 1. parameter specification
!c
!c     nx : dimension of horizontal-direction
!c    lev : dimension of z-direction
!c     dt : time step for ut, vt, tt, qt updates                (s)
!c      g : gravity                                             (m/s2)
!c      r : air gas constant                                    (j/kg)
!c     cp : air specific heat constant                          (j/kg)
!c   hltm : water vapor latent heat constant                    (j/kg)
!c   ptop : model top pressure                                  (mb)
!c   land : logic for bare soil land    (nx)  (=true for bare land)
!c   topo : terrain geopotential   (nx)                         (m2/s2)
!c    phi : odd level geopotential (nx,lev)                     (m2/s2)
!c      u : x-direction velocity   (nx,lev) at current time lvl (m/s)
!c      v : y-direction velocity   (nx,lev) at current time lvl (m/s)
!c      t : potent temp on odd lvl (nx,lev) at current time lvl (k)
!c      q : moisture    on odd lvl (nx,lev) at current time lvl (kg/kg)
!c     ut : x-direction velocity   (nx,lev) at past/future  lvl (m/s)
!c     vt : y-direction velocity   (nx,lev) at past/future  lvl (m/s)
!c     tt : real temp.  on odd lvl (nx,lev) at past/future  lvl (k)
!c     qt : moisture    on odd lvl (nx,lev) at past/future  lvl (kg/kg)
!c     pk : p**capa on odd levels  (exner func)
!c    pk2 : p**capa on even levels (exner func)
!c     sd : vertical velocity (nx,lev) mb/s
!c   qflux: upward surface moisture flux     (nx)               (w/m2)
!c   ocean: logic for open water      (nx)  (=true for open water)
!c#####################################################################
!c
      use mpe
      use rank
      use index
      use mo_constants, only:alv
      use const, only:RTYPE
!
      implicit none
!c input & output variable
!c
!      integer nx,lev,ncld,j,jj
      integer nx,nxj,lev,ncld,j,jj
      integer kcbot(nx),kctop(nx)
      real  u(nx,lev),v(nx,lev),t(nx,lev)                            &
          , qflux(nx),sd(nx,lev)                                     &
          , plt(nx,lev),ttmp(nx)
      real(kind=RTYPE) ut(nx,lev),vt(nx,lev),tt(nx,lev)              &
          ,            q(nx,lev*ncld),qt(nx,lev*ncld),phi(nx,lev)    &
          ,            topo(nx),pt(nx),sigma(lev+1,2)                &
          ,            pk(nx,lev),pk2(nx,lev)
!c
!c
!c  local work arrays
      real*8 pkxmb(nx,lev)
      real*8 pk2x(nx,lev),pkx(nx,lev)
!c
!ccumastr
      real*8 papp1(nx,lev),paphp1(nx,lev+1),pgeo(nx,lev)              &
           ,zqsat(nx,lev),pverv(nx,lev)                              &
           ,ztp1(nx,lev),zqp1(nx,lev),zup1(nx,lev),zvp1(nx,lev)      &
           ,ptte(nx,lev),pqte(nx,lev),pvom(nx,lev),pvol(nx,lev)      &
           ,zxp1(nx,lev),pxtec(nx,lev),ptop,dth,hltm,g,dt,cp,r
      real*8 pqhfl(nx),prsfc(nx),pssfc(nx),rcup(nx)
!
      real*8 rhoh2o
!c
      integer klevp1,klevm1,k,i,kc,ncldq
      logical land(nx),ldland(nx)
!xb110>
      integer kcnv(nx)
!xb110<
!c------------------------------------------------------------
      ncldq=2
      rhoh2o=1000.
!
      do k=1,lev
      do i=1,nxj
        papp1(i,k)=plt(i,k)*100.  ! from mb to pa
        paphp1(i,k)= (sigma(k,1)*pt(i)+sigma(k,2)+ptop)*100. ! from mb to pa
        pgeo(i,k) = phi(i,k) 
        pverv(i,k)=sd(i,k)*100.  ! from mb to pa
      enddo
        ttmp(:)=tt(:,k)
        call qsatq(nxj,ttmp,plt(1,k),zqsat(1,k))
      enddo
      do i=1,nxj
        paphp1(i,lev+1)= (sigma(lev+1,1)*pt(i)+sigma(lev+1,2)+ptop)*100.
      enddo

!      do k=1,lev
!!      call vlog(pk2x(1,k),pk2(1,k),nx)
!!      call vlog(pkx(1,k), pk(1,k), nx)
!      call vlog(pk2x(1,k),pk2(1,k),nxj)
!      call vlog(pkx(1,k), pk(1,k), nxj)
!!      do i=1,nx
!      do i=1,nxj
!        pk2x(i,k)=pk2x(i,k)*(cp/r)
!        pkx(i,k) = pkx(i,k)*(cp/r)
!      enddo
!!      call vexp(pk2x(1,k),pk2x(1,k),nx)
!!      call vexp(pkx(1,k), pkx(1,k), nx)
!      call vexp(pk2x(1,k),pk2x(1,k),nxj)
!      call vexp(pkx(1,k), pkx(1,k), nxj)
!      enddo
!c
!!      do  i=1,nx
!      do  i=1,nxj
!      paphp1(i,1)=ptop*100.
!!c     paphp1(i,lev+1)=ptop*100.
!      enddo
!c
!      do  k=1,lev
!!c     kc=lev-k+1
!      kc=k
!!      do  i=1,nx
!      do  i=1,nxj
!      pkxmb(i,k)=pkx(i,k)*1000. ! change into mb
!      papp1(i,kc)=pkx(i,k) * 1000.*100.    !hpa to pa
!      paphp1(i,k+1)=pk2x(i,k) * 1000.*100.     !hpa to pa
!!c     paphp1(i,kc)=pk2x(i,k) * 1000.*100.     !hpa to pa
!      enddo
!      enddo
!c
!      do  k = 1, lev
!!c     kc=lev-k+1
!      kc=k
!!      do  i = 1, nx
!      do  i = 1, nxj
!      pgeo(i,kc) = phi(i,k) 
!      pverv(i,kc)=sd(i,k)*100.  ! from mb/s to pa/s
!      enddo
!      enddo
!c
!      do  k=1,lev
!!c     kc=lev-k+1
!      kc=k
!!      call qsatq(nx,tt(1,k),pkxmb(1,k),zqsat(1,kc))
!      call qsatq(nxj,tt(1,k),pkxmb(1,k),zqsat(1,kc))
!      enddo
!c
!      do  i=1,nx
      do  i=1,nxj
      pqhfl(i)=qflux(i)/alv     !xb110 20190715, convert to upward surface moisture flux (kgm^-2s^-1)
      ldland(i)=land(i)
      enddo
!c                                                                      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c transfer t from potential temp to real temp
      do  k = 1, lev
!      do  i = 1, nx
      do  i = 1, nxj
      t(i,k) = t(i,k)*pk(i,k)
      enddo
      enddo
!c
      dth = dt
!c
      do k=1,lev
!c     kc=lev-k+1
      kc=k
!      do i=1,nx
      do i=1,nxj
        ztp1(i,kc)=tt(i,k)
        zqp1(i,kc)=qt(i,k)
        zxp1(i,kc)=qt(i,k+(ncldq-1)*lev)
        zup1(i,kc)=ut(i,k)
        zvp1(i,kc)=vt(i,k)
!c
!c       ptte(i,kc)=(tt(i,k)-t(i,k))/dth
        pqte(i,kc)=(qt(i,k)-q(i,k))/dth
!c       pvom(i,kc)=(ut(i,k)-u(i,k))/dth
!c       pvol(i,kc)=(vt(i,k)-v(i,k))/dth
        pxtec(i,kc)=(zxp1(i,k)-q(i,k+(ncldq-1)*lev))/dth
      enddo
      enddo
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!co
      klevp1=lev+1
      klevm1=lev-1
!c
!      call cumastr(nx,lev,klevp1,klevm1,ztp1,zqp1,zxp1,zup1,zvp1,ldland,   &
      call cumastr(nxj,nx,lev,klevp1,klevm1,ztp1,zqp1,zxp1,zup1,zvp1,ldland,   &
                  pverv,zqsat,pqhfl,papp1,paphp1,pgeo,ptte,pqte,pxtec,    &
                  pvom,pvol,prsfc,pssfc,kcbot,kctop,dth,j                ,&
!xb110>
                  kcnv)
!xb110<
!c
      do k=1,lev
!c     kc=lev-k+1
      kc=k
!      do i=1,nx
      do i=1,nxj
        tt(i,k)=tt(i,k)+ptte(i,kc)*dt
        qt(i,k)=qt(i,k)+pqte(i,kc)*dt
        ut(i,k)=ut(i,k)+pvom(i,kc)*dt
        vt(i,k)=vt(i,k)+pvol(i,kc)*dt
        qt(i,k+(ncldq-1)*lev)=qt(i,k+(ncldq-1)*lev)+pxtec(i,kc)*dt
      enddo
      enddo
!c
!      do i=1,nx
      do i=1,nxj
!        rcup(i)=(prsfc(i)+pssfc(i))*dt
        rcup(i)=(prsfc(i)+pssfc(i))*dt/rhoh2o
      enddo
!c
      return
      end
