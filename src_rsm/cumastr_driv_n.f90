      subroutine cumastr_driv_n  &
                   (nxj  ,nx   ,lev  ,dt   ,g    ,&
                    r    ,cp   ,hltm ,ptop ,land ,&
                    topo ,phi  ,u    ,v    ,t    ,&
                    q    ,ut   ,vt   ,tt   ,qt   ,&
                    rcup ,pk   ,pk2  ,sd   ,qflux,&
                    kcbot,kctop,ncld ,sigma,      &
                    plt  ,pt   ,j    ,lndj ,hfx  ,&
                    garea ,kcnv ,flash)
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
      use mo_constants, only:vtmpc1,alv
      use const, only:RTYPE
      implicit none
!c input & output variable
      integer nx,nxj,lev,ncld,j,jj
      integer kcbot(nx),kctop(nx)
      real  u(nx,lev),v(nx,lev),t(nx,lev)                            &
          , qflux(nx),sd(nx,lev)                                     &

          , plt(nx,lev)
      real(kind=RTYPE) ut(nx,lev),vt(nx,lev),tt(nx,lev)              &
          ,            q(nx,lev*ncld),qt(nx,lev*ncld),phi(nx,lev)    &
          ,            topo(nx),pt(nx),sigma(lev+1,2)                &
          ,            pk(nx,lev),pk2(nx,lev)
!c  local work arrays
      real pkxmb(nx,lev)
      real pk2x(nx,lev),pkx(nx,lev)
!ccumastr
      real papp1(nx,lev),paphp1(nx,lev+1),pgeo(nx,lev)              &
           ,zqsat(nx,lev),pverv(nx,lev)                              &
           ,ztp1(nx,lev),zqp1(nx,lev),zup1(nx,lev),zvp1(nx,lev)      &
           ,ptte(nx,lev),pqte(nx,lev),pvom(nx,lev),pvol(nx,lev)      &
           ,zxp1(nx,lev),pxtec(nx,lev),ptop,dth,hltm,g,dt,cp,r       &
           ,ztu(nx,lev),zqu(nx,lev),zmfu(nx,lev),zmfd(nx,lev)        &
           ,pcte(nx,lev)
      real pqhfl(nx),prsfc(nx),pssfc(nx),rcup(nx),zrain(nx)        &
            ,phhfl(nx),hfx(nx)
!
      real rhoh2o,d2r,xlatj,xlat,tt1
!c
      integer klevp1,klevm1,k,i,kc,ncldq
      integer lndj(nx)
      logical land(nx),ldland(nx)
!xb110>
      real zew,zqs,zcor,foeewm,mdlon
      logical locum(nx)
      integer kcnv(nx)
      real flash(nx),garea(nx),dx(nx)  
!xb110<
      ncldq=2
      rhoh2o=1000.
!
      do i = 1,nxj
       dx(i) = sqrt(garea(i))
      end do

      do k=1,lev
      do i=1,nxj
        papp1(i,k)=plt(i,k)*100.  ! from mb to pa
        paphp1(i,k)= (sigma(k,1)*pt(i)+sigma(k,2)+ptop)*100. ! from mb to pa
        pgeo(i,k) = phi(i,k) 
        pverv(i,k)=sd(i,k)*100.  ! from mb to pa
      enddo
      enddo
      do i=1,nxj
        paphp1(i,lev+1)= (sigma(lev+1,1)*pt(i)+sigma(lev+1,2)+ptop)*100.
      enddo

      do  i=1,nxj
      zrain(i) = 0.0
      locum(i) =.false.
      phhfl(i) = hfx(i)
      pqhfl(i) = qflux(i)/alv   !xb110 20190715, convert to upward surface moisture flux (kgm^-2s^-1)
      ldland(i)= land(i)
      enddo

      dth = dt

      do k=1,lev
      kc=k
      do i=1,nxj
        ztp1(i,kc)=tt(i,k)
        zqp1(i,kc)=qt(i,k)
        zxp1(i,kc)=qt(i,k+(ncldq-1)*lev)
        zup1(i,kc)=ut(i,k)
        zvp1(i,kc)=vt(i,k)
        pcte(i,kc)=0.0
        pvom(i,kc)=0.0
        pvol(i,kc)=0.0
        tt1=tt(i,k)
          zew  = foeewm(tt1)
          zqs  = zew/papp1(i,kc)
          zqs  = min(0.5,zqs)
          zcor = 1./(1.-vtmpc1*zqs)
          zqsat(i,kc)=zqs*zcor

        ptte(i,kc)=(tt(i,k)-t(i,k))/dth
        pqte(i,kc)=(qt(i,k)-q(i,k))/dth
        pvom(i,kc)=(ut(i,k)-u(i,k))/dth
        pvol(i,kc)=(vt(i,k)-v(i,k))/dth
        pxtec(i,kc)=(zxp1(i,k)-q(i,k+(ncldq-1)*lev))/dth
      enddo
      enddo
!============================================================
      klevp1=lev+1
      klevm1=lev-1

      call cumastr_n  &
             (nxj,   nx,   lev,  klevp1,klevm1, &
              ztp1,  zqp1, zxp1, zup1,  zvp1,   &
              ldland,pverv,zqsat,pqhfl, papp1,  &
              paphp1,pgeo, ptte, pqte,  pxtec,  &
              pvom,  pvol, prsfc,pssfc, kcbot,  &
              kctop, dth,  j,    ztu,   zqu,    &
              zmfu,  zmfd, zrain,pcte,  phhfl,  &
              lndj,  locum,dx   , kcnv,   &
              flash)

      do k=1,lev
      kc=k
      do i=1,nxj
        tt(i,k)=t(i,k)+ptte(i,kc)*dt
        qt(i,k)=max(q(i,k)+pqte(i,kc)*dt,0.)
        ut(i,k)=u(i,k)+pvom(i,kc)*dt
        vt(i,k)=v(i,k)+pvol(i,kc)*dt
        qt(i,k+(ncldq-1)*lev)=max(q(i,k+(ncldq-1)*lev)+pxtec(i,kc)*dt,0.)
      enddo
      enddo

      do i=1,nxj
        rcup(i)=(prsfc(i)+pssfc(i))*dt/rhoh2o
      enddo

      return
      end
