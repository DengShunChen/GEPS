      subroutine gridnl_hybrid_ndsl ( nxj,nx,lev,ncld                &
                        , cp,radsq,ut,vt,rdiv,tt,qt,phi,pt           &
                        , dtpl,dlpl,sinl,pk,pk2,dsigma,sigma,onocos  &
                        , cor,diveng,vdmerd,vdzonl,pten,deldm,sdpbl  &
                        , sd,pdot,vvel,sgeo)
!
!  real to spectral transformation, compute non-linear contributions
!  to spectral tendencies
!
!  **** input variables ****
!
!  nx:  e-w dimension
!  lev: no. of vertical levels
!  cp: specific heat of air
!  radsq: radius of earth squared
!  ut: gridpt e-w velocity (scaled)
!  vt: gridpt n-s velocity (scaled)
!  tt: gridpt virtual potential temp
!  qt: gridpt specific humidity
!  rdiv: gridpt divergence
!  pt: gridpt terrain pressure
!  dtpl: d(pt)/dy
!  dlpl: d(pt)/dx
!  pk: exner function on full(odd) levels
!  pk2: exner function on half(even) levels
!  dsigma: thickness of sigma layers
!  sigma: sigma level values
!  onocos: 1.0/(cos(lat)**2)
!  cor: coriolis on each gaussian latitude
!  pten: horizontal adv of terrain pressure
!
! *** output variables ***
!
!  deldm: terrain pressure tendency
!  sd: vertical velocity
!  pdot: vertical velocity
!  vvel: vertical velocity at mean layer(Pa/s)
!  diveng: energy term of divergence equation
!  vdmerd: meridional advection term of divergence/vorticity tends
!  vdzonl: zonal advection term of divergence/vorticity tends
!
!  modify to f90 by C-H Lee and sort by River Chen in 2015
!
! ******************************************************************
!
      use const, only : RTYPE
!
      implicit  none

      integer   nxj,nx,lev,ncld
      real      cp,radsq
      real(kind=RTYPE) onocos,cor,sinl

!
      real(kind=RTYPE) diveng(nx,lev),vdmerd(nx,lev),vdzonl(nx,lev), &
                pdot(nx,lev+1),pten(nx,lev),dlpl(nx),dtpl(nx),       &
                rdiv(nx,lev),ut(nx,lev),vt(nx,lev),tt(nx,lev),       &
                qt(nx,lev*ncld),phi(nx,lev),pt(nx),sgeo(nx),         &
                deldm(nx),spal(nx,lev),sd(nx,lev),sdpbl(nx),         &
                dsigma(lev,2),sigma(lev+1,2),pk(nx,lev),pk2(nx,lev), &
                cg(nx,lev),vvel(nx,lev)
!
      logical   flag(nx)
!
      integer   k,i,kbgn,kk,kkp1
      real      px,px_pbl


!CWB2014 fixed undefined value problem in diabat line 665
      sd=0.
      deldm=0.
!
!  surface pressure tendency
!
      do 2 k=1,lev-1
      do 2 i=1,nxj
        cg(i,k) = ut(i,k)*dlpl(i)*onocos+vt(i,k)*dtpl(i)
        deldm(i)= deldm(i)-dsigma(k,1)*cg(i,k)                        &
                 -rdiv(i,k)*(dsigma(k,2)+dsigma(k,1)*pt(i))
!!        deldm(i)= deldm(i)+pten(i,k)-rdiv(i,k)*(dsigma(k,2)+dsigma(k,1)*pt(i))
!        deldm(i)= deldm(i)-dsig(k)*(ut(i,k)*dlpl(i)*onocos
!     *           +vt(i,k)*dtpl(i)+rdiv(i,k)*pt(i))
        sd(i,k+1)= deldm(i)
    2 continue
!
      k= lev
      do 24 i=1,nxj
        cg(i,k) = ut(i,k)*dlpl(i)*onocos+vt(i,k)*dtpl(i)
        deldm(i)= deldm(i)-dsigma(k,1)*cg(i,k)                        &
                 -rdiv(i,k)*(dsigma(k,2)+dsigma(k,1)*pt(i))
!!        deldm(i)= deldm(i)+pten(i,k)-rdiv(i,k)*(dsigma(k,2)+dsigma(k,1)*pt(i))
!        deldm(i)= deldm(i)-dsig(k)*(ut(i,k)*dlpl(i)*onocos
!     *           +vt(i,k)*dtpl(i)+rdiv(i,k)*pt(i))
   24 continue
!
!  vertical velocity
!
      sd(:,    1) = 0.0
      pdot(:,lev+1) = 0.0
      pdot(:,    1) = 0.0
      do 3 k=2,lev
      do 3 i=1,nxj
        kk=lev-k+2
        sd(i,k)= sd(i,k)-sigma(k,1)*deldm(i)
!        sd(i,k)= sd(i,k)-sig(k)*deldm(i)
        pdot(i,kk)=sd(i,k)
    3 continue
      do k=1,lev
        kk  =lev-k+1
        kkp1=lev-k+2
        do i=1,nxj
          vvel(i,k)=0.5*((sigma(k,1)+sigma(k+1,1))*(cg(i,k)+deldm(i))  &
                  + (pdot(i,kk)+pdot(i,kkp1)))
        enddo
      enddo
!
!  obtain vertical velocity within low layers
!
      do i = 1, nxj
        flag(i) = .true.
        sdpbl(i) = sd(i,lev-1)
      enddo
!
      px_pbl = 850.
      kbgn = (lev*2)/3 - 1
      do 4 k = kbgn, lev-1
      do 4 i = 1, nxj
        px = sigma(k,1) * pt(i) + sigma(k,2)
!        px = sig(k) * pt(i)
        if( px .ge. px_pbl .and. flag(i) ) then
          sdpbl(i) = sd(i,k)
          flag(i) = .false.
        endif
    4 continue
!
!
      call vstruc_hybrid_cwb(nxj,nx,lev,cp,radsq,sigma,dsigma,pt,tt,qt &
                   ,pk,pk2,spal,phi,ncld)
!
      do 13 k=1,lev
      do 13 i=1,nxj
        vdmerd(i,k)= -spal(i,k)*dtpl(i)/onocos-ut(i,k)*cor             &
                     -(ut(i,k)*ut(i,k)+vt(i,k)*vt(i,k))*onocos*sinl
        vdzonl(i,k)= -spal(i,k)*dlpl(i)+vt(i,k)*cor

        diveng(i,k)= sgeo(i)+phi(i,k)
  13  continue
      return
      end
!
!----------------------------------------------------------------------
      subroutine gridnl_hybrid_ndsl_2tl ( nxj,nx,lev,ncld            &
                        , cp,radsq,ut,vt,rdiv,tt,qt,phi,pt           &
                        , dtpl,dlpl,sinl,pk,pk2,dsigma,sigma,onocos  &
                        , cor,diveng,vdmerd,vdzonl,pten,deldm,sdpbl  &
                        , sd,pdot,sgeo,step)
!
!  real to spectral transformation, compute non-linear contributions
!  to spectral tendencies
!
!  **** input variables ****
!
!  nx:  e-w dimension
!  lev: no. of vertical levels
!  cp: specific heat of air
!  radsq: radius of earth squared
!  ut: gridpt e-w velocity (scaled)
!  vt: gridpt n-s velocity (scaled)
!  tt: gridpt virtual potential temp
!  qt: gridpt specific humidity
!  rdiv: gridpt divergence
!  pt: gridpt terrain pressure
!  dtpl: d(pt)/dy
!  dlpl: d(pt)/dx
!  pk: exner function on full(odd) levels
!  pk2: exner function on half(even) levels
!  dsigma: thickness of sigma layers
!  sigma: sigma level values
!  onocos: 1.0/(cos(lat)**2)
!  cor: coriolis on each gaussian latitude
!  pten: horizontal adv of terrain pressure
!
! *** output variables ***
!
!  deldm: terrain pressure tendency
!  sd: vertical velocity
!  pdot: vertical velocity
!  diveng: energy term of divergence equation
!  vdmerd: meridional advection term of divergence/vorticity tends
!  vdzonl: zonal advection term of divergence/vorticity tends
!
!  modify to f90 by C-H Lee and sort by River Chen in 2015
!
! ******************************************************************
!
      use const, only : RTYPE
      implicit  none

      integer   nxj,nx,lev,ncld
      real      cp,radsq
      real(kind=RTYPE) onocos,cor,sinl

!
      real(kind=RTYPE) diveng(nx,lev),vdmerd(nx,lev),vdzonl(nx,lev),   &
                       pdot(nx,lev+1),pten(nx,lev),dlpl(nx),dtpl(nx),  &
                       rdiv(nx,lev),ut(nx,lev),vt(nx,lev),tt(nx,lev),  &
                       qt(nx,lev*ncld),phi(nx,lev),pt(nx),sgeo(nx),    &
                       deldm(nx),spal(nx,lev),sd(nx,lev),sdpbl(nx),    &
                       dsigma(lev,2),sigma(lev+1,2),pk(nx,lev),        &
                       pk2(nx,lev)
!
      logical   flag(nx)
!
      integer   k,i,kbgn,kk,step
      real      px,px_pbl

      if ( step .eq. 2 ) then
!CWB2014 fixed undefined value problem in diabat line 665
      sd=0.
!
!  surface pressure tendency
!
      k= 1
      do 22 i=1,nxj
        deldm(i)= -dsigma(k,1)*(ut(i,k)*dlpl(i)*onocos+vt(i,k)*dtpl(i))   &
                  -rdiv(i,k)*(dsigma(k,2)+dsigma(k,1)*pt(i))
!        deldm(i)= pten(i,k)-rdiv(i,k)*(dsigma(k,2)+dsigma(k,1)*pt(i))
!        deldm(i)= -dsig(k)*(ut(i,k)*dlpl(i)*onocos+vt(i,k)*dtpl(i)
!     *            +rdiv(i,k)*pt(i))
        sd(i,k+1)= deldm(i)
   22 continue
!
      do 2 k=2,lev-1
      do 2 i=1,nxj
        deldm(i)= deldm(i)-dsigma(k,1)*(ut(i,k)*dlpl(i)*onocos            &
                 +vt(i,k)*dtpl(i))-rdiv(i,k)*(dsigma(k,2)+dsigma(k,1)*pt(i))
!        deldm(i)= deldm(i)+pten(i,k)-rdiv(i,k)*(dsigma(k,2)+dsigma(k,1)*pt(i))
!        deldm(i)= deldm(i)-dsig(k)*(ut(i,k)*dlpl(i)*onocos
!     *           +vt(i,k)*dtpl(i)+rdiv(i,k)*pt(i))
        sd(i,k+1)= deldm(i)
    2 continue
!
      k= lev
      do 24 i=1,nxj
        deldm(i)= deldm(i)-dsigma(k,1)*(ut(i,k)*dlpl(i)*onocos            &
                 +vt(i,k)*dtpl(i))-rdiv(i,k)*(dsigma(k,2)+dsigma(k,1)*pt(i))
!        deldm(i)= deldm(i)+pten(i,k)-rdiv(i,k)*(dsigma(k,2)+dsigma(k,1)*pt(i))
!        deldm(i)= deldm(i)-dsig(k)*(ut(i,k)*dlpl(i)*onocos
!     *           +vt(i,k)*dtpl(i)+rdiv(i,k)*pt(i))
   24 continue
!
!  vertical velocity
!
      sd(:,    1) = 0.0
      pdot(:,lev+1) = 0.0
      pdot(:,    1) = 0.0
      do 3 k=2,lev
      do 3 i=1,nxj
        kk=lev-k+2
        sd(i,k)= sd(i,k)-sigma(k,1)*deldm(i)
!        sd(i,k)= sd(i,k)-sig(k)*deldm(i)
        pdot(i,kk)=sd(i,k)
    3 continue
!
!  obtain vertical velocity within low layers
!
      do i = 1, nxj
        flag(i) = .true.
        sdpbl(i) = sd(i,lev-1)
      enddo
!
      px_pbl = 850.
      kbgn = (lev*2)/3 - 1
      do 4 k = kbgn, lev-1
      do 4 i = 1, nxj
        px = sigma(k,1) * pt(i) + sigma(k,2)
!        px = sig(k) * pt(i)
        if( px .ge. px_pbl .and. flag(i) ) then
          sdpbl(i) = sd(i,k)
          flag(i) = .false.
        endif
    4 continue
!
      else
!
      call vstruc_hybrid_cwb(nxj,nx,lev,cp,radsq,sigma,dsigma,pt,tt,qt &
                   ,pk,pk2,spal,phi,ncld)
!
      do 13 k=1,lev
      do 13 i=1,nxj
        vdmerd(i,k)= -spal(i,k)*dtpl(i)/onocos-ut(i,k)*cor             &
                     -(ut(i,k)*ut(i,k)+vt(i,k)*vt(i,k))*onocos*sinl
        vdzonl(i,k)= -spal(i,k)*dlpl(i)+vt(i,k)*cor

        diveng(i,k)= sgeo(i)+phi(i,k)
  13  continue
      endif
      return
      end
