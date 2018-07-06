      subroutine vstruc_hybrid ( nxj,nx,lev,cp,radsq,sigma,dsigma,pt,tt,qt &
                        ,ptop,pk,pk2,spal,odpsig,that,qhat,phi,ncld )
!
!  subroutine to compute several intermediate pressure dependent
!  variables and parameters
!
! *** input ****
!
!  nx: e-w dimension no.
!  lev: number of vertical levels
!  pk: odd (full) level p**capa
!  pk2: even( half) level p**capa
!  cp: specific heat of air
!  radsq: (rad of earth)**2
!  dsig: sigma layer thicknesses
!  sig: sigma levels
!  pt: terrain pressure
!  tt: virtual potential temperature
!  qt: specific humidity
!
! *** output ***
!
!  phi: geopotential
!  spal: energy conversion term (sigma*ps*alpha)
!  odpsig: reciprical of layer mass
!  that: even level potential temp (thickness temp)
!  qhat: even level specific humidity
!
! **************************************************
!
      implicit none

      integer  nxj,nx,lev,ncld

      real     pk(nx,lev),pk2(nx,lev),spal(nx,lev),odpsig(nx,lev)       &
      , tt(nx,lev),qt(nx,lev*ncld),that(nx,lev),qhat(nx,lev*ncld)       &
      , phi(nx,lev),pt(nx),dsigma(lev,2),sigma(lev+1,2)
      real     dpk(nx,lev),dphi(nx,lev)

      real     capa,pk2top,top,cpr2,cp,radsq,ptop
      integer  i,n,k,nk,kk
!
!  compute time dependent pressure variables
!
      capa=1.0/3.5
      pk2top=(ptop/1000.)**capa

      cpr2= cp/radsq
      do 115 i=1,nxj
      odpsig(i,1)= 1.0/(dsigma(1,2)+dsigma(1,1)*pt(i))
      dpk   (i,1)= pk2(i,1)-pk2top
      spal  (i,1)= cpr2*tt(i,1)*(sigma(2,1)+sigma(1,1))/2.              &
                     *dpk(i,1)*odpsig(i,1)
  115 continue
!
!
!  half-level specific humidity, interpolate in p**capa
!
      do n=1,ncld
      nk=(n-1)*lev
      do k=1,lev-1
      kk=nk+k
      do i=1,nxj
      qhat(i,kk+1)= qt(i,kk+1)+(qt(i,kk)-qt(i,kk+1))                    &
       *(pk(i,k+1)-pk2(i,k))/(pk(i,k+1)-pk(i,k))
      enddo
      enddo
      enddo
!
      do 160 k=1,lev-1
      do 160 i=1,nxj
!
!  recipricol of layer pressure depth
!
      odpsig(i,k+1)= 1.0/(dsigma(k+1,2)+dsigma(k+1,1)*pt(i))
      dpk   (i,k+1)= pk2(i,k+1) - pk2(i,k)
!
!  half level potential temperature, defined as weighted combination
!  of full level thicknesses, not an interpolation
!
!  half level potential temperature, interpolate in p**capa
!
!  in hybrid coord. proposed by Hnery Juang, we can choose
!  a better way to define that
!
      that(i,k+1)= tt(i,k+1)+(tt(i,k)-tt(i,k+1))                        &
       *(pk(i,k+1)-pk2(i,k))/(pk(i,k+1)-pk(i,k))
!
!  energy conversion term for terrain pressure contribution to
!  horizontal pressure gradient
!
      spal(i,k+1)= cpr2*tt(i,k+1)*(sigma(k+1,1)+sigma(k+2,1))/2.        &
                 *dpk(i,k+1)*odpsig(i,k+1)
!
  160 continue
!
!  hydrostatic equation: note that terrain geopotential is excluded.
!  it it constant forcing term that is included when laplacian of
!  geopotential is computed in divergence equation.
!
      do k=1,lev
        do i=1,nxj
          dphi(i,k)= cp*tt(i,k)*dpk(i,k)
        enddo
      enddo
!
      do i=1,nxj
        phi(i,lev)= dphi(i,lev)/2.
      enddo

      do k=lev-1,1,-1
        do i=1,nxj
          phi(i,k)= phi(i,k+1)+(dphi(i,k)+dphi(i,k+1))/2.
        enddo
      enddo
!
      return
      end
