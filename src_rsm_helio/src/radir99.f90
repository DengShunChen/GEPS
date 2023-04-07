      subroutine radir99( nr,nx,lev,fluxcl,ee,pu                      &
                       ,pa,ta,tg,wa,oa,umco2,umch4,umn2o,cld,clwc     &
                       ,cre,ciwc,cde,cufrac,fuir,fdir,fuirr,fdirr )
!
!#####################################################################
!
!  (1)compute infrared radiation fluxes at each model level using 
!  correlated k-distribution method to resolve gas absorption 
!  coefficients (Fu and Liou, 1992) and delta two and four-stream 
!  approximation to compute fluxes (Fu et al., 1998; Fong, 1998).
!
!  (2)the delta two and four-stream approximation includes the
!  partial cloud treatment consistent with that used in the solar
!  radiation model.
!
!  ++ input/output variables :
!     (1) input :
!      nr    : number of grid points to calculate ir flux ( <= nx )
!      nx    : horizontal dimension of input array
!      lev   : number of vertical layers
!      fluxcl : logical const. to perform cloud forcing calculation
!      ee(nx,ibir) : ir surface emissivity at each band
!      pu(nx*(lev+1)) : presure at level levels
!      pa(nx*lev)     : presure at odd levels ( prediction levels )
!      ta(nx*lev)     : temperature at odd levels
!      tg(nx)         : surface ground temperature ( k )
!      wa(nx*lev)     : moisture (mixing ratio) specific humility ( kg/kg )
!      oa(nx*lev)     : ozone mixing ratio from climatology ( kg/kg )
!      umco2 : constant concentration of co2 ( ppmv )
!      umch4 : constant concentration of ch4 ( ppmv )
!      umn2o : constant concentration of n2o ( ppmv )
!      cld(nx*lev)    : cloud fraction at each layer
!      clwc(nx*lev)   : liquid water content at each layer ( g/m3 )
!      cre(nx*lev)    : effective size of water clouds ( um )
!      ciwc(nx*lev)   : ice water content at each layer ( g/m3 )
!      cde(nx*lev)    : effective size of ice clouds ( um )
!      cufrac(nx) : cumulus-related cloud fraction for chooseing
!                    cloud overlapping approximation
!
!     (2) output : (nx*(lev+1))
!      fuir   : all-sky ir upward flux at each level ( w/m2 )
!      fdir   : all-sky ir downward flux at each level ( w/m2 )
!      fuirr  : clear-sky ir upward flux at each level ( w/m2 )
!      fdirr  : clear-sky ir downward flux at each level ( w/m2 )
!
!#########################################################################
!
!  ibir : total band number in ir region
!  igir : total k-interval (g-value) number in ir region
! 
      implicit none
      integer, parameter :: ibir = 12, igir = 67 

      integer  nr,nx,lev
      logical  fluxcl
      real     ee(nx,ibir),pu(nx*(lev+1)),pa(nx*lev),ta(nx*lev),     &
               tg(nx),wa(nx*lev),oa(nx*lev),cld(nx*lev),clwc(nx*lev),&
               cre(nx*lev),ciwc(nx*lev),cde(nx*lev),                 &
               cufrac(nx),fuir(nx*(lev+1)),fdir(nx*(lev+1)),         &
               fuirr(nx*(lev+1)),fdirr(nx*(lev+1))
!
!  common block, data comes from "rkdisir.f" 
!
      real     hk1(2), c1h2o(3,19,2)
      real     hk2(3), c2h2o(3,19,3)
      real     hk3(4), c3h2o(3,19,4)
      real     hk4(4),c4h2o(3,19,4),c4ch4(3,19),c4n2o(3,19)
      real     hk5(3),c5h2o(3,19,3),c5ch4(3,19),c5n2o(3,19)
      real     hk6(5),c6o3(3,19,5),c6h2o(3,19)
      real     hk7(2), c7h2o(3,19,2)
      real     hk8(10), c8hca(3,19,10), c8hcb(3,19,10)
      real     hk9(12), c9hca(3,19,12), c9hcb(3,19,12)
      real     hk10(7), c10h2o(3,19,7)
      real     hk11(7), c11h2o(3,19,7)
      real     hk12(8), c12h2o(3,19,8)

      common /ir1/ hk1, c1h2o
      common /ir2/ hk2, c2h2o
      common /ir3/ hk3, c3h2o
      common /ir4/ hk4,c4h2o,c4ch4,c4n2o
      common /ir5/ hk5,c5h2o,c5ch4,c5n2o
      common /ir6/ hk6,c6o3,c6h2o
      common /ir7/ hk7, c7h2o
      common /ir8/ hk8, c8hca, c8hcb
      common /ir9/ hk9, c9hca, c9hcb
      common /ir10/ hk10, c10h2o
      common /ir11/ hk11, c11h2o
      common /ir12/ hk12, c12h2o
!
!  local working arrays
!
      logical ozon
      real    hk(igir),vv(12),pc(5,12)
      integer igb(12)
      real    dp(nx*lev),tu(nx*(lev+2))
      real    zw1(nx*lev,2),zw2(nx*lev,2),                      &
              zw(nx*lev,2),ztau(nx*lev,2),zttau(nx*(lev+1),2),  &
              zbf(nx*(lev+2)),zffu(nx*(lev+1)),zffd(nx*(lev+1))
      real    taugas(nx*lev,igir),tauic(nx*lev),tauwc(nx*lev),  &
              taucon(nx*lev)
      real    wic1(nx*lev),wic2(nx*lev),wwc1(nx*lev),wwc2(nx*lev),&
              wic(nx*lev),wwc(nx*lev)
      real    fkg(nx*lev,igir),fuq(nx),wk1(nx*lev),wk2(nx*(lev+2),3),&
              dz(nx*lev)
      real    cldb(nx*lev,4)

      integer nrn,nrn1,i,j,jk,jj,jje,m,igt,ig,ib
      real    pi,vmin,umn2o,umch4,umco2,h1269,h476,cld_maxran,dzmax_km,dlogpa
      real    dlogp,x1,x2,rch4,rn2o,rco2,x,y,z,r,trg,s,pe,w,fk,tics,twcs,ttts

!
!  igb : the max. number of g-values contained in the band 1-12 
!
      data igb /2,3,4,4,3,5,2,10,12,7,7,8/
!
! "taucon" are the optical depthes due to water vapor continuum absorp-
! tion in lev layers for a given band ib. We include continuum absorp-
! tion in the 280 to 1250 cm**-1 region. vv(11)-vv(17) are the central
! wavenumbers of each band in this region.
!
      data vv / 4*0.0,1175.0,1040.0,890.0,735.0,605.0,470.0,340.0,0.0 /
!
!  coefficients for computing the blackbody intensity function
!
       data pc /                                                        &
       9.598466e-9, -7.287285e-6, 2.090916e-3,-2.677473e-1, 1.286976e+1,& 
       7.211260e-9, -4.792760e-6, 1.202248e-3,-1.344489e-1, 5.639765e+0,&
       7.572943e-9, -2.626821e-6, 2.216726e-5, 7.445105e-2,-6.670347e+0,&
       4.730481e-11, 2.706418e-6,-1.351777e-3, 2.282604e-1,-1.300338e+1,&
      -3.297827e-9,  5.815023e-6,-2.261827e-3, 3.382716e-1,-1.779915e+1,&
      -4.643682e-9,  6.200410e-6,-2.132982e-3, 2.915934e-1,-1.426057e+1,&
      -7.868328e-9,  9.028216e-6,-2.644667e-3, 3.072150e-1,-1.274817e+1,&
      -4.047419e-9,  3.969190e-6,-7.282698e-4, 2.771434e-2, 1.444413e+0,&
      -1.476049e-9,  9.066730e-7, 4.060300e-4,-1.256673e-1, 8.467850e+0,&
       5.730060e-10,-1.242197e-6, 1.034899e-3,-1.776529e-1, 9.197895e+0,&
       8.067974e-10,-1.134579e-6, 6.529138e-4,-8.009666e-2, 2.954545e+0,&
       3.569705e-10,-4.619423e-7, 2.387227e-4,-1.647910e-3,-7.843924e-1/
!
!  476.16 = 2.24e4 / M * 10.0 / 9.8, where M = 48 for o3
!  1269.841 = 2.24e4 / M * 10.0 / 9.8, where M = 18 for h2o
!
      data h476/476.16/, h1269/1269.841/
!
!  if "cufrac(i)" .gt. "cld_maxran", assume maximum-random cloud overlapping
!  if "cufrac(i)" .le. "cld_maxran", assume random cloud overlapping
!  so, set cld_maxran equal to less than 0. if assume the maximum-random cloud
!  overlapping approximation always
!
      data cld_maxran / 0. /
!
      nrn  = nr*lev
      nrn1 = nr*(lev+1)
      pi = 4.*atan(1.)
      vmin = 1.0e-8
!
!
!  obtain the temperature (tu) at even levels 
!
!  lev+1 : the surface skin temperature
!  lev+2 : the ground temperature
!  dp : layer thickness in pressure (hPa)
!  dz : layer thickness (km)
!
!  assume the temperature at the top even level equal to that at
!  the top odd level
!
      do 20 i = 1, nr
      tu(i) = ta(i)
   20 continue
!
! note : dz in ir radiation uses unit of km instead of m in solar
!
      dzmax_km = 1.   ! avoid too thick cloud depth
!
      do 25 j = 1, lev
      jk = j*nr
      jj = jk-nr
      do 25 i = 1, nr
      dp(i+jj) = pu(i+jk) - pu(i+jj)
      dz(i+jj) = 0.0292674 * ta(i+jj) * log( pu(i+jk) / pu(i+jj) )
!
      dz(i+jj) = min( dzmax_km, dz(i+jj) )
!
   25 continue
!
      do 30 i = 1, nr*(lev-1)
      dlogpa = log ( pa(i+nr) ) - log ( pa(i) )
      dlogp  = log ( pu(i+nr) ) - log ( pa(i) )
      tu(i+nr) = ta(i) + ( ta(i+nr) - ta(i) ) * dlogp / dlogpa
   30 continue
!
!  assume the lowest model layer is well mixed (i.e. theta = const. )
!
      jje = (lev-1)*nr 
      do 35 i = 1, nr
      x1 = log ( pa(i+jje) )
      x2 = log ( pu(i+nrn) )
      tu(i+nrn)  = ta(i+jje) * exp ( 0.286*(x2-x1) )
      tu(i+nrn) = (2.0*tg(i)+tu(i+nrn))/3.0
      tu(i+nrn1) = tg(i)
   35 continue
!
!  compute cloud optical properties
!
      do 60 m = 1, 4
      do 60 i = 1, nrn
      cldb(i,m) = 1.0
   60 continue
!
!  compute cloud overlapping variables assuming mixed maximum
!  and random overlap based on the method proposed by Geleyn and 
!  Hollingsworth (1979, contrib. atmos. phys.)
!
      jj = nrn - nr
      do 68 i = 1, nr
      cldb(i,1) = 1.0 - cld(i)
      cldb(i,3) = cld(i)
      cldb(i+jj,2) = 1.0 - cld(i+jj)
      cldb(i+jj,4) = cld(i+jj)
   68 continue
      do 69 j = 1, lev-1
      jk = j*nr
      jj = jk-nr
      do 69 i = 1, nr
      if ( cufrac(i) .ge. cld_maxran ) then
      cldb(i+jk,1) = max( vmin, 1.0-max( cld(i+jk),cld(i+jj) ) ) /  &
                     max( vmin, 1.0-cld(i+jj) )
      cldb(i+jk,3) = max( vmin, min( cld(i+jk),cld(i+jj) ) ) /      &
                     max( vmin, cld(i+jj) )
      cldb(i+jj,2) = max( vmin, 1.0-max( cld(i+jj),cld(i+jk) ) ) /  &
                     max( vmin, 1.0-cld(i+jk) )
      cldb(i+jj,4) = max( vmin, min( cld(i+jj),cld(i+jk) ) ) /      &
                     max( vmin, cld(i+jk) )
      else
      cldb(i+jk,1) = 1.0 - cld(i+jk)
      cldb(i+jk,3) = cld(i+jk)
      cldb(i+jj,2) = 1.0 - cld(i+jj)
      cldb(i+jj,4) = cld(i+jj)
      end if
   69 continue
!
      igt = 0
      ozon = .false.
!****************************************************************
!
!  compute the gas optical depth for a given g-value in the band 1-12
!
!  (1) band 1 ( 0 - 280 cm^-1 )
!
      call rqkir ( nr,nx,lev,igb(1),ozon,c1h2o,pa,ta,fkg )   
!    
      do 210 ig = 1, igb(1)
      igt = igt + 1
      hk(igt) = hk1(ig)
      do 210  i = 1, nrn
      taugas(i,igt) = fkg(i,ig) * wa(i) * dp(i) * h1269 
  210 continue
!
!  (2) band 2 ( 280 - 400 cm^-1 )
!
      call rqkir ( nr,nx,lev,igb(2),ozon,c2h2o,pa,ta,fkg )   
!
      do 220 ig = 1, igb(2)
      igt = igt + 1
      hk(igt) = hk2(ig)
      do 220  i = 1, nrn
      taugas(i,igt) = fkg(i,ig) * wa(i) * dp(i) * h1269
  220 continue
!
!  (3) band 3 ( 400 - 540 cm^-1 )
!
      call rqkir ( nr,nx,lev,igb(3),ozon,c3h2o,pa,ta,fkg )
!
      do 230 ig = 1, igb(3)
      igt = igt + 1
      hk(igt) = hk3(ig)
      do 230  i = 1, nrn
      taugas(i,igt) = fkg(i,ig) * wa(i) * dp(i) * h1269
  230 continue
!
!  (4) band 4 ( 540 - 670 cm^-1 )
!
      call rqkir ( nr,nx,lev,igb(4),ozon,c4h2o,pa,ta,fkg(1,1) )
      call rqkir ( nr,nx,lev,1,ozon,c4ch4,pa,ta,fkg(1,11) )
      call rqkir ( nr,nx,lev,1,ozon,c4n2o,pa,ta,fkg(1,12) )
!
!  1.26238e-3 = 2.24e4 / M * 10.0 / 9.8 * 1.6e-6 * M / 28.97, where
!               M = 16 for CH4.
!  2.20918e-4 = 2.24e4 / M * 10.0 / 9.8 * 0.28e-6 * M / 28.97, where
!                M = 44 for N2O.
!
      rch4 = umch4 / 1.6
      rn2o = umn2o / 0.28
      do 240 i = 1, nrn
!err  wk1(i) = ( fkg(i,11)*1.26238e-3 + fkg(i,12)*2.20918e-4 ) * dp(i)
      wk1(i) = (fkg(i,11)*1.26238e-3*rch4 + fkg(i,12)*2.20918e-4*rn2o) * dp(i)
  240 continue
      do 241 ig = 1, igb(4)
      igt = igt + 1
      hk(igt) = hk4(ig)
      do 241  i = 1, nrn
      taugas(i,igt) = fkg(i,ig) * wa(i) * dp(i) * h1269 + wk1(i)
  241 continue
!
!  (5) band 5 ( 670 - 800 cm^-1 )
!
      call rqkir ( nr,nx,lev,igb(5),ozon,c5h2o,pa,ta,fkg(1,1) )
      call rqkir ( nr,nx,lev,1,ozon,c5ch4,pa,ta,fkg(1,11) )
      call rqkir ( nr,nx,lev,1,ozon,c5n2o,pa,ta,fkg(1,12) )
!
      do 250 i = 1, nrn
!err  wk1(i) = ( fkg(i,11)*1.26238e-3 + fkg(i,12)*2.20918e-4 ) * dp(i)
      wk1(i) = (fkg(i,11)*1.26238e-3*rch4 + fkg(i,12)*2.20918e-4*rn2o) * dp(i)
  250 continue
      do 251 ig = 1, igb(5)
      igt = igt + 1
      hk(igt) = hk5(ig)
      do 251  i = 1, nrn
      taugas(i,igt) = fkg(i,ig) * wa(i) * dp(i) * h1269 + wk1(i)
  251 continue
!
!  (6) band 6 ( 800 - 980 cm^-1 )
!
      ozon = .true.
      call rqkir ( nr,nx,lev,igb(6),ozon,c6o3,pa,ta,fkg(1,1) )
      ozon = .false.
      call rqkir ( nr,nx,lev,1,ozon,c6h2o,pa,ta,fkg(1,11) )
!
      do 260 i = 1, nrn
      wk1(i) = fkg(i,11) * wa(i) * dp(i) * h1269
  260 continue
      do 261 ig = 1, igb(6)
      igt = igt + 1
      hk(igt) = hk6(ig)
      do 261  i = 1, nrn
      taugas(i,igt) = fkg(i,ig) * oa(i) * dp(i) * h476 + wk1(i)
  261 continue
!
!  (7) band 7 ( 980 - 1100 cm^-1 )
! 
      call rqkir ( nr,nx,lev,igb(7),ozon,c7h2o,pa,ta,fkg )
!
      do 270 ig = 1, igb(7)
      igt = igt + 1
      hk(igt) = hk7(ig)
      do 270  i = 1, nrn
      taugas(i,igt) = fkg(i,ig) * wa(i) * dp(i) * h1269
  270 continue
!
!  (8) band 8 ( 1100 - 1250 cm^-1 )
!
      call rqkir ( nr,nx,lev,igb(8),ozon,c8hca,pa,ta,fkg(1,1) )
      call rqkir ( nr,nx,lev,igb(8),ozon,c8hcb,pa,ta,fkg(1,igb(8)+1) )
!
      rco2 = umco2 / 330.
      do 280 ig = 1, igb(8)
      igt = igt + 1
      hk(igt) = hk8(ig)
      do 280  i = 1, nrn
      taugas(i,igt) = (fkg(i,ig)*rco2 + wa(i)*fkg(i,ig+igb(8))) * dp(i)
  280 continue
!
!  (9) band 9 ( 1250 - 1400 cm^-1 )
!
      call rqkir ( nr,nx,lev,igb(9),ozon,c9hca,pa,ta,fkg(1,1) )
      call rqkir ( nr,nx,lev,igb(9),ozon,c9hcb,pa,ta,fkg(1,igb(9)+1) )
!
      do 290 ig = 1, igb(9)
      igt = igt + 1
      hk(igt) = hk9(ig)
      do 290  i = 1, nrn
      taugas(i,igt) = (fkg(i,ig)*rco2 + wa(i)*fkg(i,ig+igb(9))) * dp(i)
  290 continue
!
!  (10) band 10 ( 1400 - 1700 cm^-1 )
!
      call rqkir ( nr,nx,lev,igb(10),ozon,c10h2o,pa,ta,fkg )
!
      do 300 ig = 1, igb(10)
      igt = igt + 1
      hk(igt) = hk10(ig)
      do 300  i = 1, nrn
      taugas(i,igt) = fkg(i,ig) * wa(i) * dp(i) * h1269
  300 continue
!
!  (11) band 11 ( 1700 - 1900 cm^-1 )
!
      call rqkir ( nr,nx,lev,igb(11),ozon,c11h2o,pa,ta,fkg )
!
      do 310 ig = 1, igb(11)
      igt = igt + 1
      hk(igt) = hk11(ig)
      do 310  i = 1, nrn
      taugas(i,igt) = fkg(i,ig) * wa(i) * dp(i) * h1269
  310 continue
!
!  (12) band 12 ( 1900 - 2200 cm^-1 )
!
      call rqkir ( nr,nx,lev,igb(12),ozon,c12h2o,pa,ta,fkg )
!
      do 320 ig = 1, igb(12)
      igt = igt + 1
      hk(igt) = hk12(ig)
      do 320  i = 1, nrn
      taugas(i,igt) = fkg(i,ig) * wa(i) * dp(i) * h1269
  320 continue
!
!*******************************************************************
!
      do 380 i = 1, nrn1
      fuir(i) = 0.
      fdir(i) = 0.
      fuirr(i) = 0.
      fdirr(i) = 0.
  380 continue
!
      igt = 0

      do 600 ib = 1, ibir

!
!  compute the blackbody intensity function integrated over a given
!  band at temperature t, the units are w/m2/sr
!
      do 400 i = 1, nr*(lev+2)
      wk2(i,1) = tu(i) * tu(i)
      wk2(i,2) = wk2(i,1) * tu(i)
      wk2(i,3) = wk2(i,1) * wk2(i,1)
  400 continue
      do 405 i = 1, nr*(lev+2)
      zbf(i) = pc(1,ib)*wk2(i,3) + pc(2,ib)*wk2(i,2) +       &
               pc(3,ib)*wk2(i,1) + pc(4,ib)*tu(i) + pc(5,ib)
      zbf(i) = max ( 1.e-3, zbf(i) )
  405 continue
!
!  compute the optical depth due to water vapor continuum absorption
!  the units of fk are cm2/g, see eq. (a.19) of fu (1991)
!
      x = 4.18
      y = 5577.8
      z = 0.00787
      r = 0.002
      trg = 10. / 9.80616
      if( ib .gt. 4 .and. ib .lt. 12 ) then
        s = ( x + y * exp ( -z * vv(ib) ) ) / 1013.25
        do 410  i = 1, nrn
        pe = pa(i) * wa(i) / ( 0.622 + 0.378 * wa(i) )
        w = exp ( 1800.0 / ta(i) - 6.08108 )
        fk = s * ( pe + r * pa(i) ) * w
        taucon(i) = ( fk * wa(i) )* dp(i) * trg
  410   continue
      else
        do 412 i = 1, nrn
        taucon(i) = 0.
  412   continue
      end if
!
!  obtain cloud optical properties
!
      call roptir( nr,nx,lev,ib,dz,cld,clwc,ciwc,cre,cde,tauwc,tauic  &
                  ,wwc,wic,wwc1,wwc2,wic1,wic2 )
!

      do 500 ig = 1, igb(ib)
      igt = igt + 1
!
!  combine single-scattering properties due to ice and water clouds
!  along with nongray gasous absorption
!  m=1 : clear part
!  m=2 : cloudy part
!
      do 415 m = 1, 2
      do 415 i = 1, nr
      zttau(i,m) = 0.
  415 continue
!
      do 420 i = 1, nrn
	ztau(i,1) = taugas(i,igt) + taucon(i)
	ztau(i,2) = taugas(i,igt) + taucon(i) + tauwc(i) + tauic(i) 
	zw(i,1) = 0.
	zw(i,2) = ( tauic(i)*wic(i) + tauwc(i)*wwc(i) ) / ztau(i,2)
  420 continue
!
      do 430 m = 1, 2
      do 430 i = 1, nrn
      if ( zw(i,m) .lt. 1.0e-10 ) then
	zw1(i,m) = 0.
	zw2(i,m) = 0.
      else
	tics = tauic(i) * wic(i)
	twcs = tauwc(i) * wwc(i)
	ttts = tics + twcs 
	zw1(i,m) = ( tics*wic1(i)+twcs*wwc1(i) ) / ttts
	zw2(i,m) = ( tics*wic2(i)+twcs*wwc2(i) ) / ttts
      end if
  430 continue
!     
      do 440 m = 1, 2
      do 440 j = 1, lev
      jk = j*nr
      jj = jk-nr
      do 440 i = 1, nr
      zttau(i+jk,m) = zttau(i+jj,m) + ztau(i+jj,m)
  440 continue
!
!  compute ir fluxes 
!
      call rcmpir24c(nr,nx,lev,zbf,ee(1,ib),zw1,zw2,zw,zttau,cld,cldb  &
                   ,zffu,zffd )
!
!
!  multiply the ir energy fraction contained in a g-value
!
      do 450 i = 1, nrn1
      fuir(i) = fuir(i) + zffu(i)*hk(igt)
      fdir(i) = fdir(i) + zffd(i)*hk(igt)
  450 continue
!
      if( fluxcl ) then
!
!  clear-sky case
!
      call rcmpir24(nr,nx,lev,zbf,ee(1,ib),zw1(1,1),zw2(1,1),zw(1,1)  &
                   ,zttau(1,1),zffu,zffd )
      do 480 i = 1, nrn1
      fuirr(i) = fuirr(i) + zffu(i)*hk(igt)
      fdirr(i) = fdirr(i) + zffd(i)*hk(igt)
  480 continue
!
      end if
!
  500 continue
  600 continue
!
!  fuq is the surface emitted flux in the band 0 - 280 cm**-1 with a
!  hk of 0.03, which isn't absorbed by atmosphere. 
!  ( note : the summation of all hk in band 12 only equalto 0.97 )
!
      do 750 i = 1, nr 
      fuq(i) = zbf(i+nrn1) * 0.03 * pi * ee(i,12)
  750 continue
      do 780 j = 1, lev+1
      jj = (j-1)*nr
      do 780 i = 1, nr
      fuir(i+jj) = fuir(i+jj) + fuq(i)
      fuirr(i+jj) = fuirr(i+jj) + fuq(i)
  780 continue
!cc
      return
      end
