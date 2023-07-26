      subroutine gridnl_hybrid ( nxj,nx,lev,ncld                        &
                        , cp,radsq,ut,vt,rdiv,rvor,tt,qt,phi            &
                        , pt,dtpl,dlpl,pk,pk2,dsigma,sigma,tbar,qbar    &
                        , onocos,cor,deldm,ddtemp,tadvu,qvadv           &
                        , qadvu,diveng,tadvv,qadvv,vdmerd,vdzonl        &
                        , sdpbl,sd,ptop )
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
!  rdiv: gridpt divergence
!  rvor: gridpt vorticity
!  tt: gridpt virtual potential temp
!  qt: gridpt specific humidity
!  phi: geopotential
!  pt: gridpt terrain pressure
!  dtpl: d(pt)/dy
!  dlpl: d(pt)/dx
!  pk: exner function on full(odd) levels
!  pk2: exner function on half(even) levels
!  dsigma: thickness of sigma layers
!  sigma: sigma level values
!  tbar: global mean temperature on each odd model level
!  qbar: global mean spec humid on each odd model level
!  onocos: 1.0/(cos(lat)**2)
!  cor: coriolis on each gaussian latitude
!
! *** output variables ***
!
!  deldm: terrain pressure tendency
!  ddtemp: vertical advection of temp tendency
!  tadvu: zonal advection part of temp tendency
!  qvadv: vertical advection of moisture tendency
!  qadvu: zonal advection part of moisture tendency
!  diveng: energy term of divergence equation
!  tadvv: meridianal advection part of temp tendency
!  qadvv: meridional advection part of moisture tendency
!  vdmerd: meridional advection term of divergence/vorticity tends
!  vdzonl: zonal advection term of divergence/vorticity tends
!
!  modify to f90 by C-H Lee and sort by River Chen in 2015
!
! ******************************************************************
!
      implicit  none

      integer   nxj,nx,lev,ncld
      real      cp,radsq,onocos,cor,ptop

      real      ut(nx,lev),vt(nx,lev),rdiv(nx,lev),rvor(nx,lev)         &
      , tt(nx,lev),qt(nx,lev*ncld),pt(nx),dlpl(nx),dtpl(nx),pk(nx,lev)  &
      , phi(nx,lev),pk2(nx,lev),dsigma(lev,2),sigma(lev+1,2),tbar(lev)  &
      , qbar(lev*ncld)
!
      real      deldm(nx),ddtemp(nx  ,lev),tadvu(nx  ,lev)              &
      , qvadv(nx  ,lev*ncld),qadvu(nx  ,lev*ncld),diveng(nx  ,lev)      &
      , tadvv(nx  ,lev),qadvv(nx  ,lev*ncld),vdmerd(nx  ,lev)           &
      , vdzonl(nx  ,lev)
!
      real      sdpbl(nx)
      logical   flag(nx)
!
      real      sd(nx,lev),that(nx,lev),odpsig(nx,lev),                 &
                qhat(nx,lev*ncld),spal(nx,lev)

      integer   k,i,kbgn,n,nk,kk
      real      px,px_pbl

!
!CWB2014 fixed undefined value problem in rstrantq and rstrandz loop 23
      tadvv=0.
      tadvu=0.
      qadvv=0.
      qadvu=0.
      vdmerd=0.
      vdzonl=0.

!CWB2014 fixed undefined value problem in diabat line 665
      sd=0.
!
!  surface pressure tendency
!
      k= 1
      do 22 i=1,nxj
        deldm(i)= -dsigma(k,1)*(ut(i,k)*dlpl(i)*onocos+vt(i,k)*dtpl(i))   &
                  -rdiv(i,k)*(dsigma(k,2)+dsigma(k,1)*pt(i))
!        deldm(i)= -dsig(k)*(ut(i,k)*dlpl(i)*onocos+vt(i,k)*dtpl(i)
!     *            +rdiv(i,k)*pt(i))
        sd(i,k+1)= deldm(i)
   22 continue
!
      do 2 k=2,lev-1
      do 2 i=1,nxj
        deldm(i)= deldm(i)-dsigma(k,1)*(ut(i,k)*dlpl(i)*onocos            &
                 +vt(i,k)*dtpl(i))-rdiv(i,k)*(dsigma(k,2)+dsigma(k,1)*pt(i))
!        deldm(i)= deldm(i)-dsig(k)*(ut(i,k)*dlpl(i)*onocos
!     *           +vt(i,k)*dtpl(i)+rdiv(i,k)*pt(i))
        sd(i,k+1)= deldm(i)
    2 continue
!
      k= lev
      do 24 i=1,nxj
        deldm(i)= deldm(i)-dsigma(k,1)*(ut(i,k)*dlpl(i)*onocos            &
                 +vt(i,k)*dtpl(i))-rdiv(i,k)*(dsigma(k,2)+dsigma(k,1)*pt(i))
!        deldm(i)= deldm(i)-dsig(k)*(ut(i,k)*dlpl(i)*onocos
!     *           +vt(i,k)*dtpl(i)+rdiv(i,k)*pt(i))
   24 continue
!
!  vertical velocity
!
      sd(:,    1) = 0.0
      do 3 k=2,lev
      do 3 i=1,nxj
        sd(i,k)= sd(i,k)-sigma(k,1)*deldm(i)
!        sd(i,k)= sd(i,k)-sig(k)*deldm(i)
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
      call vstruc_hybrid_cwb(nxj,nx,lev,cp,radsq,sigma,dsigma,pt,tt,qt  &
                   ,ptop,pk,pk2,spal,odpsig,that,qhat,phi,ncld)
!
!  temperature tendency
!
      k=1
      do 6 i=1,nxj
        ddtemp(i,k)= rdiv(i,k)*(tt(i,k)-tbar(k))+sd(i,k+1)                &
                    *(tt(i,k)-that(i,k+1))*odpsig(i,k)
    6 continue
!
      do 61 k=2,lev-1
      do 61 i=1,nxj
        ddtemp(i,k)= rdiv(i,k)*(tt(i,k)-tbar(k))+(sd(i,k+1)*(tt(i,k)      &
                    -that(i,k+1))+sd(i,k)*(that(i,k)-tt(i,k)))*odpsig(i,k)
   61 continue
!
      k=lev
      do 62 i=1,nxj
        ddtemp(i,k)= rdiv(i,k)*(tt(i,k)-tbar(k))+sd(i,k)                  &
                    *(that(i,k)-tt(i,k))*odpsig(i,k)
   62 continue
!
!  horizontal advection of temperature
!
      do 8 k=1,lev
      do 8 i=1,nxj
        tadvv(i,k)= (tt(i,k)-tbar(k))*vt(i,k)
        tadvu(i,k)= (tt(i,k)-tbar(k))*ut(i,k)
    8 continue
!
!  'vdmerd': meridional contribution to vorticity and divergence tendenc
!  'vdzonl': zonal contribution to vorticity and divergence tendencies
!
      k=1
      do 12 i=1,nxj
        vdmerd(i,k)= spal(i,k)*dtpl(i)/onocos+0.5*sd(i,k+1)               &
                  *(vt(i,k+1)-vt(i,k))*odpsig(i,k)+ut(i,k)*(rvor(i,k)+cor)
        vdzonl(i,k)= vt(i,k)*(rvor(i,k)+cor)-spal(i,k)*dlpl(i)            &
                  -0.5*sd(i,k+1)*(ut(i,k+1)-ut(i,k))*odpsig(i,k)
   12 continue
!
      do 120 k=2,lev-1
      do 120 i=1,nxj
        vdmerd(i,k)= spal(i,k)*dtpl(i)/onocos+0.5*(sd(i,k+1)              &
                    *(vt(i,k+1)-vt(i,k))+sd(i,k)*(vt(i,k)-vt(i,k-1)))     &
                    *odpsig(i,k)+ut(i,k)*(rvor(i,k)+cor)
        vdzonl(i,k)= vt(i,k)*(rvor(i,k)+cor)-spal(i,k)*dlpl(i)            &
        -0.5*(sd(i,k+1)*(ut(i,k+1)-ut(i,k))+sd(i,k)*(ut(i,k)-ut(i,k-1)))  &
                   *odpsig(i,k)
  120 continue
!
      k=lev
      do 121 i=1,nxj
        vdmerd(i,k)= spal(i,k)*dtpl(i)/onocos                             &
        +0.5*sd(i,k)*(vt(i,k)-vt(i,k-1))*odpsig(i,k)                      &
        +ut(i,k)*(rvor(i,k)+cor)
        vdzonl(i,k)= vt(i,k)*(rvor(i,k)+cor)-spal(i,k)*dlpl(i)            &
        -0.5*sd(i,k)*(ut(i,k)-ut(i,k-1))*odpsig(i,k)
  121 continue
!
!  vertical advection of moisture
!
      do n=1,ncld
        nk=(n-1)*lev
        k=1
        kk=nk+k
        do i=1,nxj
          qvadv(i,kk)= rdiv(i,k)*(qt(i,kk)-qbar(kk))+sd(i,k+1)              &
                      *(qt(i,kk)-qhat(i,kk+1))*odpsig(i,k)
        enddo
        do k=2,lev-1
          kk=nk+k
          do i=1,nxj
            qvadv(i,kk)= rdiv(i,k)*(qt(i,kk)-qbar(kk))+(sd(i,k+1)*(qt(i,kk)   &
                    -qhat(i,kk+1))+sd(i,k)*(qhat(i,kk)-qt(i,kk)))*odpsig(i,k)
          enddo
        enddo
        k= lev
        kk=nk+k
        do i=1,nxj
          qvadv(i,kk)= rdiv(i,k)*(qt(i,kk)-qbar(kk))+sd(i,k)                &
                      *(qhat(i,kk)-qt(i,kk))*odpsig(i,k)
        enddo
      enddo
!
!  horizontal moisture advection
!
      do n=1,ncld
        nk=(n-1)*lev
        do k=1,lev
          kk=nk+k
          do i=1,nxj
            qadvv(i,kk)= vt(i,k)*(qt(i,kk)-qbar(kk))
            qadvu(i,kk)= ut(i,k)*(qt(i,kk)-qbar(kk))
          enddo
        enddo
      enddo
!
!  laplacian of energy term in divergence equation
!
      do 85 k=1,lev
      do 85 i=1,nxj
        diveng(i,k)= 0.5*radsq*onocos*(ut(i,k)*ut(i,k)                    &
                    +vt(i,k)*vt(i,k))+phi(i,k)
   85 continue
!
      return
      end
