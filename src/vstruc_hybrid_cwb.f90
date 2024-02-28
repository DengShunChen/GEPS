      subroutine vstruc_hybrid_cwb (nxj,nx,lev,cp,radsq,sigma,dsigma   &
                  ,pt,tt,qt,pk,pk2,spal,phi,ncld)
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
!
! **************************************************
!
      use const, only : RTYPE
!
      implicit  none

      integer   nxj,nx,lev,ncld

      real(kind=RTYPE) tt(nx,lev),qt(nx,lev*ncld),phi(nx,lev),pt(nx)   &
      , odpsig(nx,lev),spal(nx,lev),that(nx,lev),dsigma(lev,2)         &
      , sigma(lev+1,2),pk(nx,lev),pk2(nx,lev)

      real      cpr2,cp,radsq
      integer   i,n,nk,k,kk
!
!  compute time dependent pressure variables
!
      cpr2= cp/radsq
      do 115 i=1,nxj
      odpsig(i,1)= 1.0/(dsigma(1,2)+dsigma(1,1)*pt(i))
      spal(i,1)= cpr2*tt(i,1)*sigma(2,1)*(pk2(i,1)-pk(i,1))*odpsig(i,1)
!      odpsig(i,1)= 1.0/(dsig(1)*pt(i))
!      spal(i,1)= cpr2*tt(i,1)*sig(2)*(pk2(i,1)-pk(i,1))*odpsig(i,1)
  115 continue
!
!
!  half-level specific humidity, interpolate in p**capa
!
!      do n=1,ncld
!      nk=(n-1)*lev
!      do k=1,lev-1
!      kk=nk+k
!      do i=1,nxj
!      qhat(i,kk+1)= qt(i,kk+1)+(qt(i,kk)-qt(i,kk+1))     &
!       *(pk(i,k+1)-pk2(i,k))/(pk(i,k+1)-pk(i,k))
!      enddo
!      enddo
!      enddo
!
      do 160 k=1,lev-1
      do 160 i=1,nxj
!
!  recipricol of layer pressure depth
!
      odpsig(i,k+1)= 1.0/(dsigma(k+1,2)+dsigma(k+1,1)*pt(i))
!      odpsig(i,k+1)= 1.0/(dsig(k+1)*pt(i))
!
!  half-level specific humidity, interpolate in p**capa
!
!x    qhat(i,k+1)= qt(i,k+1)+(qt(i,k)-qt(i,k+1))
!x   * *(pk(i,k+1)-pk2(i,k))/(pk(i,k+1)-pk(i,k))
!
!  half level potential temperature, defined as weighted combination
!  of full level thicknesses, not an interpolation
!
      that(i,k)= tt(i,k)-(tt(i,k)-tt(i,k+1))   &
       *(pk(i,k+1)-pk2(i,k))/(pk(i,k+1)-pk(i,k))
!
!  geopotential thicknesses
!
      phi(i,k)= that(i,k)*(pk(i,k+1)-pk(i,k))
!
!  energy conversion term for terrain pressure contribution to
!  horizontal pressure gradient
!
      spal(i,k+1)= cpr2*tt(i,k+1)*(sigma(k+1,1)*(pk(i,k+1)-pk2(i,k)) &
       +sigma(k+2,1)*(pk2(i,k+1)-pk(i,k+1)))*odpsig(i,k+1)
!      spal(i,k+1)= cpr2*tt(i,k+1)*(sig(k+1)*(pk(i,k+1)-pk2(i,k))
!     * +sig(k+2)*(pk2(i,k+1)-pk(i,k+1)))*odpsig(i,k+1)
!
  160 continue
!
!  hydrostatic equation: note that terrain geopotential is excluded.
!  it it constant forcing term that is included when laplacian of
!  geopotential is computed in divergence equation.
!
      do 100 i=1,nxj
      phi(i,lev)= cp*tt(i,lev)*(pk2(i,lev)-pk(i,lev))
  100 continue
!
      do 105 k=lev-1,1,-1
      do 105 i=1,nxj
      phi(i,k)= phi(i,k+1)+cp*phi(i,k)
  105 continue
!
      return
      end
