      subroutine radsl99( ns,nx,lev,fluxcl,s0,u0,as            &
                       ,pu,pa,ta,wa,oa,cld,clwc,cre,ciwc,cde   &
                       ,cufrac,fusl,fdsl,fuslr,fdslr )
!
!#####################################################################
!
!  (1)Compute solar radiation heating rates and net fluxes at the top
!  and surface using the delta-eddington approximation.
!
!  (2)The treatment of spetral bands(11 bands), gas absorption and cloud
!  optical property totally follows the context of Chou and Suarez (1999).
!  But, for simplicity we ingor the minor gas absorption by co2 and o2,
!  aerosol and rain drop effects. They can be introduced into the code easily.
!
!  (3)The treatment of the cloud overlapping, however, is based on the
!  method proposed by Geleyn and Hollingsworth (1979, contrib. atmos. phys.).
!
!
!  ++ input/output variables :
!   (1) input :
!       ns    : number of grid points to calculate sw fluxes
!       nx    : horizontal dimension of input array
!       lev   : number of vertical layers
!       fluxcl : logical const. to calculate clear-sky sw fluxes for
!                performing cloud radiative forcing purpose
!       s0    : solar const.
!       u0(nx) : cos(zenith angle)
!       as(nx,ibsl) : ground surface albedo
!       pu(nx*(lev+1)) : presure at level levels
!       pa(nx*lev) : presure at odd levels ( prediction levels )
!       ta(nx*lev) : temperature at odd levels
!       wa(nx*lev) : moisture (mixing ratio) specific humility ( kg/kg )
!       oa(nx*lev) : ozone mixing ratio from climatology ( kg/kg )
!       cld(nx*lev) : cloud fraction at each layer
!       clwc(nx*lev): liquid water content at each layer ( g/m3 )
!       cre(nx*lev) : effective size of water clouds ( um )
!       ciwc(nx*lev): ice water content at each layer ( g/m3 )
!       cde(nx*lev) : effective size of ice clouds ( um )
!       cufrac(nx) : cumulus-related cloud fraction for chooseing
!                    cloud overlapping approximation
!
!  (2) output : (nx*(lev+1))
!     fusl  : all-sky solar upward flux at each level ( w/m2 )
!     fdsl  : all-sky solar downward flux at each level ( w/m2 )
!     fuslr : clear-sky solar upward flux at each level ( w/m2 )
!     fdslr : clear-sky solar downward flux at each level ( w/m2 )
!
!#########################################################################
!
      implicit none
      integer, parameter :: nk=10, ibo3=8, ibsl=11, igsl=38 
!
      integer  ns,nx,lev
      logical fluxcl
      real    s0 
      real    u0(nx),as(nx,ibsl),pu(nx*(lev+1)),pa(nx*lev),           &
              ta(nx*lev),wa(nx*lev),oa(nx*lev),cld(nx*lev),           &
              clwc(nx*lev),cre(nx*lev),ciwc(nx*lev),cde(nx*lev),      &
              cufrac(nx),fusl(nx*(lev+1)),fdsl(nx*(lev+1)),           &
              fuslr(nx*(lev+1)),fdslr(nx*(lev+1))
!
!  local working arrays
!
!cc   include 'param.h'
!
      integer igb(ibsl)

      real    hks0(igsl),rayk(ibsl),co3k(ibo3),cwk(ibo3),             &
              ch2ok(nk)
      real    zw1(nx*lev,2),zw2(nx*lev,2),zw(nx*lev,2),ztau(nx*lev,2),&
              zttau(nx*(lev+1),2),zffu(nx*(lev+1)),zffd(nx*(lev+1))
      real    tauray(nx*lev),tauic(nx*lev),tauwc(nx*lev),             &
              taugas(nx*lev),taurn(nx*lev)
      real    gic(nx*lev),gwc(nx*lev),grn(nx*lev),                    &
              wic(nx*lev),wwc(nx*lev),wrn(nx*lev)
      real    cldb(nx*lev,4),waa(nx*lev),oaa(nx*lev),dp(nx*lev),      &
              dz(nx*lev),clwp(nx*lev),ciwp(nx*lev)
      real    aib(ibsl,2),awb(ibsl,2),aia(ibsl,3),awa(ibsl,3),        &
              aig(ibsl,3),awg(ibsl,3),ara(ibsl),arg(ibsl)

      integer nsn,nsn1,i,j,jk,jj,m,igt,ib,ig
      real    h476,cld_maxran,vmin,dzmax,dzz,scal
      real    reffi,reffw,tics,twcs,ttts,ggc
!
!  igb : the number of k intervals contained in the band 1 to band 11
!
      data igb /8*1,10,10,10/
!
!  hks0 : the frational extra-terrestrial solar flux in each k interval
!
      data hks0/                                                        &
! band 1 - band 8
        .00057, .00367, .00083, .00417, .00600, .00556, .05913, .39081  &
! band 9
      , .20673, .03497, .03011, .02260, .01336, .00696, .00441, .00115  &
      , .00026, .00000                                                  &
! band 10
      , .08236, .01157, .01133, .01143, .01240, .01258, .01381, .00650  &
      , .00244, .00094                                                  &
! band 11
      , .01074, .00360, .00411, .00421, .00389, .00326, .00499, .00465  &
      , .00245, .00145/
!
!  rayk : the extinction coefficient for Rayleigh scattering.
!         unit: 1/mb.

      data rayk /.00604, .00170, .00222, .00132, .00107, .00091         &
               , .00055, .00012, .000016, .0000018, .000000/

!
!  co3k : the ozone absorption coefficient. unit=1/(cm-atm)stp
!
!     The effect of the absorption due to ozone in the region
!     with wavelength grester than 0.7 micron is taken into
!     consideration by adding the ozone absorption coefficient in
!     the 0.4-0.7 micron band by a value of 0.0033.  The effect is
!     to increase the atmospheric absorption by 0.5 W/m^2 for the MLS
!     and a solar zenith angle of 60 degrees.
!
      data co3k /30.47, 187.2,  301.9,   42.83                          &
               , 7.09,  1.25,   0.0345,  0.0572/
!
!  cwk : the water vapor absorption coefficient. unit: cm**2/g
!
      data cwk /7*0.0, 0.00070/
!
!  ch2ok: water vapor absorption coefficient for 10 k-intervals
!         for ir band.  unit: cm^2/gm
!
      data ch2ok /0.0010, 0.0133, 0.0422, 0.1334, 0.4217                &
               , 1.334,  5.623,  31.62,  177.8,  1000.0/
!
!  coefficients for computing the extinction coefficient of
!  ice clouds from b=aib(*,1)+aib(*,2)/reff
!
      data aib /11*0.000333, 11*2.52/
!
!  coefficients for computing the extinction coefficient of
!  water clouds from b=awb(*,1)+awb(*,2)/reff

      data awb /8*-0.00659, -0.0101, -0.0166, -0.0339                   &
              , 8*    1.65,    1.72,    1.85,    2.16/
!
!  coefficients for computing the single scattering albedo of
!  ice clouds from ssa=1-(aia(*,1)+aia(*,2)*reff+aia(*,3)*reff**2)

      data aia/8*0., -.00000260, .00215346, .08938331                   &
              ,8*0.,  .00000746, .00073709, .00299387                   &
              ,8*0.,  .00000000,-.00000134,-.00001038/
!
!  coefficients for computing the single scattering albedo of
!  liquid clouds from ssa=1-(awa(*,1)+awa(*,2)*reff+awa(*,3)*reff**2)

      data awa/8*0., .00000007,-.00019934, .01209318                    &
              ,8*0., .00000845, .00088757, .01784739                    &
              ,8*0.,-.00000004,-.00000650,-.00036910/
!
!  single-scattering coalbedo for rain drops fixed at a certain
!  drop distribution.  Data is provided by Q. Fu.

      data ara/8*0., .029, .342, .466/
!
!  coefficients for computing the asymmetry factor of ice clouds
!  from asycl=aig(*,1)+aig(*,2)*reff+aig(*,3)*reff**2

      data aig/8* .74625000, .74935228, .76098937, .84090400            &
              ,8* .00105410, .00119715, .00141864, .00126222            &
              ,8*-.00000264,-.00000367,-.00000396,-.00000385/

!  coefficients for computing the asymmetry factor of liquid clouds
!  from asycl=awg(*,1)+awg(*,2)*reff+awg(*,3)*reff**2

      data awg/8* .82562000, .79375035, .74513197, .83530748            &
              ,8* .00529000, .00832441, .01370071, .00257181            &
              ,8*-.00014866,-.00023263,-.00038203, .00005519/
!
!  asymmetry factor for rain drops fixed at a certain
!  drop distribution.  Data is provided by Q. Fu.

      data arg/8* 0.883, 0.891, 0.948, 0.971/
!
!  476.16 = 2.24e4 / M * 10.0 / 9.8, where M = 48 for o3
!
      data h476/476.16/
!
!  if "cufrac(i)" .gt. "cld_maxran", assume maximum-random cloud overlapping
!  if "cufrac(i)" .le. "cld_maxran", assume random cloud overlapping
!  so, set cld_maxran equal to less than 0. if assume the maximum-random cloud
!  overlapping approximation always
!
      data cld_maxran / 0. /  ! zero means maximum-random overlap always
!
      nsn  = ns*lev
      nsn1 = ns*(lev+1)
      vmin = 1.0e-8
!
      do 10 i = 1, nsn1
      fusl(i) = 0.
      fdsl(i) = 0.
      fuslr(i) = 0.
      fdslr(i) = 0.
   10 continue
!
!  pu : pressure at even level
!  dp : layer thickness in pressure (hPa)
!  dz : layer thickness (m)
!
      dzmax = 1000.    ! avoid too thick cloud depth
!
      do 25 j = 1, lev
      jk = j*ns
      jj = jk-ns
      do 25 i = 1, ns
      dp(i+jj) = pu(i+jk) - pu(i+jj)
      dz(i+jj) = 29.2674 * ta(i+jj) * log( pu(i+jk)/pu(i+jj) )
      dzz = min( dzmax, dz(i+jj) )
      clwp(i+jj) = dzz * clwc(i+jj)
      ciwp(i+jj) = dzz * ciwc(i+jj)
   25 continue
!
!  compute scaled water vapor amount, unit is g/cm**2
!  compute ozone amount,unit is (cm-atm)stp
!  the constant h476(476.16) is a conversion factor from g/cm**2
!  to (cm-atm)stp
!
      do 30 i = 1, nsn
      scal = exp( 0.8*log(pa(i)/300.) )
      waa(i) = 1.02*wa(i)*scal*(1.+0.00135*(ta(i)-240.))*dp(i)
      oaa(i) = h476*oa(i)*dp(i)
   30 continue
!
!  compute cloud overlapping variables depending on the overlapping
!  approximation
!
      do 35 m = 1, 4
      do 35 i = 1, nsn
      cldb(i,m) = 1.0
   35 continue
!
      jj = nsn - ns
      do 40 i = 1, ns
      cldb(i,1) = 1.0 - cld(i)
      cldb(i,3) = cld(i)
      cldb(i+jj,2) = 1.0 - cld(i+jj)
      cldb(i+jj,4) = cld(i+jj)
   40 continue
!
      do 45 j = 1, lev-1
      jk = j*ns
      jj = jk-ns
      do 45 i = 1, ns
!
      if ( cufrac(i) .ge. cld_maxran ) then
!
!  maximum-random cloud overlapping approximation
!
      cldb(i+jk,1) = max( vmin, 1.0-max( cld(i+jk),cld(i+jj) ) ) / &
                     max( vmin, 1.0-cld(i+jj) )
      cldb(i+jk,3) = max( vmin, min( cld(i+jk),cld(i+jj) ) ) /     &
                     max( vmin, cld(i+jj) )
      cldb(i+jj,2) = max( vmin, 1.0-max( cld(i+jj),cld(i+jk) ) ) / &
                     max( vmin, 1.0-cld(i+jk) )
      cldb(i+jj,4) = max( vmin, min( cld(i+jj),cld(i+jk) ) ) /     &
                     max( vmin, cld(i+jk) )
      else
!
!  random cloud overlapping approximation
!
      cldb(i+jk,1) = 1.0 - cld(i+jk)
      cldb(i+jk,3) = cld(i+jk)
      cldb(i+jj,2) = 1.0 - cld(i+jj)
      cldb(i+jj,4) = cld(i+jj)
      end if
!
   45 continue

      igt = 0

      do 500 ib = 1, ibsl

!
!  compute rayleigh optical thickness
!
       do 50 i = 1, nsn
       tauray(i) = rayk(ib) * dp(i)
   50  continue
!
!
!  compute cloud optical properties
!
      do 60 i = 1, nsn
      tauwc(i) = 0.
      tauic(i) = 0.
      taurn(i) = 0.
      wwc(i)   = 0.
      wic(i)   = 0.
      wrn(i)   = 0.
      gwc(i)   = 0.
      gic(i)   = 0.
      grn(i)   = 0.
   60 continue

      do 70 i = 1, nsn
      if( cld(i) .gt. 0.01 .and. ciwp(i) .gt. 1.e-5 ) then
       reffi = min( cde(i), 130. )
!cc    tauic(i) = ciwc(i)*dz(i)*( aib(ib,1)+aib(ib,2)/reffi )
       tauic(i) = ciwp(i)*( aib(ib,1)+aib(ib,2)/reffi )
       wic(i) = 1. - ( aia(ib,1)+(aia(ib,2)+aia(ib,3)*reffi)*reffi )
       gic(i) = aig(ib,1)+(aig(ib,2)+aig(ib,3)*reffi)*reffi
      end if
!
      if( cld(i) .gt. 0.01 .and. clwp(i) .gt. 1.e-5 ) then
       reffw = min( cre(i), 20. )
!cc    tauwc(i) = clwc(i)*dz(i)*( awb(ib,1)+awb(ib,2)/reffw )
       tauwc(i) = clwp(i)*( awb(ib,1)+awb(ib,2)/reffw )
       wwc(i) = 1. - ( awa(ib,1)+(awa(ib,2)+awa(ib,3)*reffw)*reffw )
       gwc(i) = awg(ib,1)+(awg(ib,2)+awg(ib,3)*reffw)*reffw
      end if
!rain taurn(i) = crnc(i)*dz(i)*0.00307
!rain wrn(i) = 1. - ara(ib)
!rain grn(i) = arg(ib)
   70 continue
!
      do 400 ig = 1, igb(ib)
      igt = igt + 1

      if( ib .le. ibo3 ) then
!
!  compute ozone and water vapor optical thickness
!
       do 80 i = 1, nsn
       taugas(i) = co3k(ib)*oaa(i) + cwk(ib)*waa(i)
   80  continue
!
      else
!
       do 90 i = 1, nsn
       taugas(i) = ch2ok(ig) * waa(i)
   90  continue
!
      end if
!
!  combine single-scattering properties due to ice and water clouds
!  and rayleigh molecules along with nongray gasous absorption
!  m=1 : clear
!  m=2 : overcast
!
      do 215 m = 1, 2
      do 215 i = 1, ns
      zttau(i,m) = 0.
  215 continue
!
!  ztau : layer effective optical thickness
!  zw  : layer mean single scattering albedo
!  zw1 : the first expansion coefficient of the phase function
!  zw2 : the second expansion coefficient of the phase function
!  the Henyey-Greenstein approximation is applied
!
      do 220 i = 1, nsn
      ztau(i,1) = taugas(i) + tauray(i) + vmin
      ztau(i,2) = taugas(i) + tauray(i) + tauwc(i) + tauic(i) + vmin
!rain ztau(i,2) = taugas(i)+tauray(i)+tauwc(i)+tauic(i)+taurn(i)
      zw(i,1) = tauray(i) / ztau(i,1)
      zw(i,2) = ( tauic(i)*wic(i)+tauwc(i)*wwc(i)+tauray(i) )/ztau(i,2)
!rain zw(i,2) = ( tauic(i)*wic(i)+tauwc(i)*wwc(i)+tauray(i)
!rain            +taurn(i)*wrn(i) )/ztau(i,2)
  220 continue
!
      do 230 i = 1, nsn
      zw1(i,1) = 0.
      zw2(i,1) = 0.
  230 continue
!
      do 235 i = 1, nsn
      tics = tauic(i) * wic(i)
      twcs = tauwc(i) * wwc(i)
!rain trns = taurn(i) * wrn(i)
      ttts = tics + twcs + vmin
!rain ttts = tics + twcs + trns + vmin
      ggc = ( tics*gic(i)+twcs*gwc(i) )/ttts
!rain ggc = ( tics*gic(i)+twcs*gwc(i)+trns*grn(i) )/ttts
      zw1(i,2) = 3. * ggc
      zw2(i,2) = 5. * ggc * ggc
  235 continue
!
      do 250 m = 1, 2
      do 250 j = 1, lev
      jk = j*ns
      jj = jk-ns
      do 250 i = 1, ns
      zttau(i+jk,m) = zttau(i+jj,m) + ztau(i+jj,m)
  250 continue
!
!  compute normalized solar fluxes
!
!  partly cloud fraction case
!
      call rcmpsl2c( ns,nx,lev,u0,as(1,ib),zw1,zw2,zw,zttau,cld,cldb &
                    ,zffu,zffd )
!
!
!  multiply the solar energy contained in each k-interval
!
      do 300 i = 1, nsn1
      fusl(i) = fusl(i) + zffu(i)*hks0(igt)*s0
      fdsl(i) = fdsl(i) + zffd(i)*hks0(igt)*s0
  300 continue
!
      if( fluxcl ) then
!
!  clear-sky case
!
      call rcmpsl2 ( ns,nx,lev,u0,as(1,ib),zw1(1,1),zw2(1,1),zw(1,1) &
                    ,zttau(1,1),zffu,zffd )
      do 320 i = 1, nsn1
      fuslr(i) = fuslr(i) + zffu(i)*hks0(igt)*s0
      fdslr(i) = fdslr(i) + zffd(i)*hks0(igt)*s0
  320 continue
      end if
!
  400 continue
  500 continue
!cc
      return
      end
