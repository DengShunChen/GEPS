      subroutine radtn99 ( fluxcl,ozon,nxj,nx,lev,ncld,lvlw,julian    &
                       , stbo,s0,g                                    &
                       , cp,ptop,dsigma,sinl,cosz,alb,tg,curate,pst   &
                       , pl,tt,qtx,o3lx,plcl,cumtop,ss,rs             &
                       , asol_clr,olr_clr,ss_clr,rs_clr               &
                       , dtrad,asr,alr,asr_clr,alr_clr                &
                       , xsr,xlr,acld,aflxd,aflxu,il,ib,cof           &
                       , sdpbl,ctot,rld,sld,chig,cmid,clow,asl,atl    &
!---------------------------------------------------------------------
                       , afusl,afdsl,afuir,afdir                      &
                       , afuslr,afdslr,afuirr,afdirr                  &
                       , asl_clr,atl_clr,clds,rld_clr,sld_clr)
!######################################################################c
!                                                                      c
! a. function:                                                         c
!  main radiation driver, does both long and short wave radiation      c
!  calculations.
!  solar radiation scheme is a two-scheme method and parameterizations
!  of gas absorption and cloud optical property are based on fu and liou
!  , jas 1993.
!  longwave radiation scheme is delta two and four-stream based on
!  work of Fu et.al. 1998.
!                                                                      c
! b. common block specification:                                       c
!  none                                                                c
!                                                                      c
! c. parameter specification:                                          c
!  none                                                                c
!                                                                      c
! d. input/output variables                                            c
!  (1) input variables:                                                c
!
!  fluxcl: logical constant to perform cloud forcing calculation       c
!          if .T., perform calculation for both clear and cloudy skys
!          if .F., perform calculation only for cloudy sky
!  ozon  : logical constant to include ozon in sw/lw rad calculation   c
!  nx    : horizontal dimension of the model                           c
!  lev   : vertical layers of the model                                c
!  lvlw  : horizontal dimension for coarse resolution long wave calcu. c
!  julian: julian day ( julian = 1 for january 1 )                     c
!  stbo  : stefan-boltzman constant (5.669e-8)                         c
!  s0    : solar constant (1368.3)                                     c
!  g     : gravity constant (9.806)                                    c
!  cp    : air heat constant (1004.)                                   c
!  ptop  : model top pressure (mb) (must be > 0.0)                     c
!  dsigma(lev): sigma layer thickness                                    c
!  sinl  : sin of latitude of current ring                             c
!  cosz(nx): zenith angle at grid points                               c
!  alb(nx) : surface albedo                                            c
!  tg(nx)  : ground real temperature (k)                               c
!  curate(nx): cup precipitation rate (mm/day) -- cwb std units --     c
!  pst(nx): = (terrain pressure - ptop) (mb)                           c
!  pl(nx,lev) : pressure at prediction levels (mb)                     c
!  tt(nx,lev) : real temperature at prediction levels (k)              c
!  qt(nx,lev) : specific humidity at prediction levels (kg/kg)         c
!  o3l(nx,lev): ozone mixing ratio from climatology    (kg/kg)         c
!  plcl(nx)  : pressure at cloud bottom (mb)                           c
!  cumtop(nx): pressure at cloud top (mb)                              c
!  il(nx)    : address index for spline interp. (expand to fine grids) c
!  ib(lvlw)  : address index for compress to coarse grids              c
!  cof(lvlw,3): constant coefficients for spline interpolation         c
!
!  (2) output variables:                                               c
!  plcl(nx)  : total column solar radiation absorption (input changed) c
!  cumtop(nx): outgoing long wave radiation at model top (input changed)
!  ss(nx): net short wave flux down to ground surface (wat/m**2)       c
!  rs(nx): net long  wave flux up from ground surface (wat/m**2)       c
!  dtrad(nx,lev): total temp change due to radiation (k/day)           c
!  asr(lev): average temp change due to short wave radtn (k/day)       c
!  alr(lev): average temp change due to long  wave radtn (k/day)       c
!  xsr(lev): maximum heating     due to short wave radtn (k/day)       c
!  xlr(lev): maximum cooling     due to long  wave radtn (k/day)       c
!  acld(lev): average total cloud fraction assuming random overlap     c
!  aflxd(lev+2): average downward long wave flux (wat/m**2)            c
!  aflxu(lev+2): average upward   long wave flux (wat/m**2)            c
!
!  ** for clear-sky
!  asol_clr(nx), olr_clr(nx), ss_clr(nx), rs_clr(nx)
!  asr_clr(lev), alr_clr(lev)
!                                                                      c
!  modify to f90 by C-H Lee and sort by River Chen in 2015
!                                                                      c
!######################################################################c
!
      use const, only: RTYPE
!
      implicit none
      integer, parameter :: ibsl = 11, ibir = 12 
!
      logical fluxcl,ozon
      integer nxj,nx,lev,ncld,lvlw,julian
      real    stbo,s0,cp,ptop

      real    cosz(nx),alb(nx),tg(nx),curate(nx),                      &
              pl(nx,lev),qt(nx,lev),o3l(nx,lev),plcl(nx),              &
              cumtop(nx),ss(nx),rs(nx),dtrad(nx,lev),asr(lev),alr(lev),&
              xsr(lev),xlr(lev),acld(lev),aflxd(lev+2),aflxu(lev+2),   &
              cof(lvlw,3),ctot(nx),chig(nx),cmid(nx),clow(nx),         &
              ttmp(nx,lev)

      integer il(nx),ib(lvlw)

      real(kind=RTYPE) qtx(nx,lev*ncld),o3lx(nx,lev),pst(nx),sdpbl(nx) &
              ,        dsigma(lev,2),sigma(lev+1,2),tt(nx,lev)
!
!  clear part
!
      real    ss_clr(nx),rs_clr(nx),asol_clr(nx),olr_clr(nx), &
              asr_clr(lev),alr_clr(lev)
!
!     local work arrays
!
      real    asl(nx,lev),atl(nx,lev),clds(nx,lev),fluxdd(nx,lev+3), &
              fluxuu(nx,lev+2)
      real    wk1(nx),wk2(nx)
      real    cvclds(nx,lev),cufrac(nx)
!
!  working arrays for new radiation scheme
!
      real    pa(nx*lev),ta(nx*lev),wa(nx*lev),oa(nx*lev),tga(nx),  &
              clda(nx*lev),cufa(nx),clwca(nx*lev),ciwca(nx*lev),    &
              crea(nx*lev),cdea(nx*lev)
      real    pu(nx*(lev+1)),u0(nx),as(nx,ibsl),ee(nx,ibir)
      real    fusl(nx*(lev+1)),fdsl(nx*(lev+1)),fuslr(nx*(lev+1)),  &
              fdslr(nx*(lev+1)),fuir(nx*(lev+1)),fdir(nx*(lev+1)),  &
              fuirr(nx*(lev+1)),fdirr(nx*(lev+1))
      real    asl_clr(nx,lev),atl_clr(nx,lev)
      integer ibs(nx)

!cmy----
      real(kind=RTYPE) afusl(nx,lev+1) ,afdsl(nx,lev+1),afuslr(nx,lev+1), &
                       afdslr(nx,lev+1),afuir(nx,lev+1),afdir(nx,lev+1) , &
                       afuirr(nx,lev+1),afdirr(nx,lev+1)
!cmy----

      real    clwc(nx,lev),ciwc(nx,lev),cre(nx,lev),cde(nx,lev),    &
              cldb(nx,lev)
      real    qtr(nx,lev)
!
!  for soil model with the need to obtain surface downward ir flux
!
      real    rld(nx),sld(nx)
!----------------------------------------------------------------------------
!  for clear part
!----------------------------------------------------------------------------
      real    rld_clr(nx),sld_clr(nx)
!----------------------------------------------------------------------------
!
      logical lwfull

      integer i,j,k,nxlev,nqt,ls,lhig,lmid,nsol,jk,jj,nb,nsn,nxjl,nrn
      real    sinl,g,umn2o,umch4,umco2,fac,qtz,fac_o3,www,ptopc
      real    cld_maxran,vmin,xxx,xim
      real    cldmaxv
!
      real    d2r,r2d,deg_ju,deg_lat,cosl2,arg

!
! --  ls: number of top levels in which no clouds are allowed
!
      data ls /1/  ! don't change
!
!  for co2, ch4 and n2o, uniform mixing is assumed through the
!  atmosphere with concentrations of 330, 1.6 and  0.28 ppmv,
!  respectively.
!
!ccc  data umco2 / 330.0 /, umch4 / 1.6 /, umn2o / 0.28 /
!
!  increase the co2 concentrations for 330 to 360 based on IPCC report
!  ( ~375 ppm at the year of 2005 )
!  2008/10/7

      data umco2 / 360.0 /, umch4 / 1.6 /, umn2o / 0.28 /
!
      nxlev= nx*lev
      lwfull= lvlw.eq.nxj
      fac= 86400.0*g/(cp*100.0)
!
!  compute sigma values at even levels
!
      sigma(1,1) = 0.0
      sigma(1,2) = 0.0
      do 130 k = 1, lev
      sigma(k+1,1) = sigma(k,1) + dsigma(k,1)
      sigma(k+1,2) = sigma(k,2) + dsigma(k,2)
  130 continue
!
! --- ensure qt > 1.0e-12 (i.e. positive specific humidity)
!
      do 135 k = 1, lev
      do 135 i = 1, nxj
      qt(i,k) = max ( 1.0e-12, qtx(i,k) )
      o3l(i,k)= max ( 1.0e-12, o3lx(i,k) )
      qtr(i,k)= 0.
  135 continue
!
      if( ncld .ge. 2 )then
      do k = 1, lev
      do i = 1, nxj
       qtr(i,k)= qtx(i,lev+k)
      enddo
      enddo
      endif
!
!  avoid too much water vapor at upper layers
!  2002/10/31
!
!  modify this water vapor limit with an opposite direction
!  of setting a min amount for water vapor of upper layers
!  because that too much warm bias was found in summer
!  2008/10/7
!
      nqt=0
!
      if( nqt .eq. 2 )then
        ttmp(:,:)=tt(:,:)
        call qsatq_2d ( nxj,nx,lev,ttmp,pl,asl )
        do k = 1, lev/2
        do i = 1, nxj
!jh          if( pl(i,k) .lt. 150. )then
!jh          if( pl(i,k) .lt. 70. )then
          if( pl(i,k) .lt. 50. )then
            qtz = min( 5.e-6, asl(i,k)*0.2 )
            qt(i,k) = max( qtz,  qt(i,k) )
          endif
!jh         if( pl(i,k) .lt. 50. ) qt(i,k) = min( 5.e-5, qt(i,k) )
         if( pl(i,k) .lt. 30. ) qt(i,k) = min( 5.e-5, qt(i,k) )
        enddo
        enddo
      endif
!
!  reduce o3 amount at the top layer
!  2002/10/31
!
!  modify the reducing factor for 0.5 to 0.8
!  2008/10/7
!
      if(ncld.eq.2)then
        do k=1,8
          fac_o3=k*0.1
          if(fac_o3.le.0.3)fac_o3=0.3
        do i = 1, nxj
!         o3l(i,k) = o3l(i,1)*0.8
         o3l(i,k) = o3l(i,k)*fac_o3
        end do
        end do
      endif
!
!  initialize cloud optical parameters
!
      cre = 5.0
      cde = 20.
      clwc = 0.
      ciwc = 0.
!
!  zero out arrays for clear parts
!
      ss_clr = 0.
      rs_clr = 0.
      asol_clr = 0.
      olr_clr = 0.
!--------------------------------------------------------------------------
      rld_clr = 0.
      sld_clr = 0.
!--------------------------------------------------------------------------
      asl_clr = 0.
      atl_clr = 0.
      fuslr = 0.
      fdslr = 0.
      afuslr = 0.
      afdslr = 0.
!
!
!     diagnose total cloud fraction and cloud optical property
!     (hold relative humidity in atl, qsat in asl)
!
      ttmp(:,:)=tt(:,:)
      call qsatq_2d ( nxj,nx,lev,ttmp,pl,asl )
!
      do 140 k = 1, lev
      do 140 i = 1, nxj
      atl(i,k) = qt(i,k) / asl(i,k)
      if ( atl(i,k) .gt. 1.0 )  atl(i,k) = 1.0
  140 continue
!
!     convert curate from 'mm/day'  to  'cm/hr' for clouds diagnose
!
      do 146 i = 1, nxj
      www = curate(i) * 0.0416667 * 0.1
      www = max ( www, 5.531e-04 )
      wk1(i) = 0.93 + 0.124*log(www)
  146 continue
!
      d2r   = 3.141592654 / 180.0
      r2d   = 180. / 3.141592654
      deg_lat=asin(sinl)*r2d
      deg_ju=23.45*sin(d2r*(360./365.)*(julian+284.))
      arg=deg_lat-deg_ju
      if(arg.gt.90.)then
        arg=89.999
      else if(arg.lt.-90.)then
        arg=-89.999
      endif
      cosl2=cos(d2r*arg)**2
      ptopc= 200. - 130.*cosl2
      cldmaxv= 0.90 - 0.20*(1.-cosl2)
!      ptopc = 250. - 165.*(1.-sinl*sinl)
!
      do k=1,lev
        do i=1,nxj
!          clds(i,k)=0.999-0.08*cosl2    !a3 (in rcloud: name rhc)
!          if(clds(i,k).ge.0.98)clds(i,k)=0.98 ! (Zhao and Carr 1997)
          clds(i,k)=asl(i,k)   ! staturation water vapor (Xu and Randall 1996)
        enddo
      enddo
!
      call rcloud( nxj,nx,lev,ls,ptop,ptopc,pst,pl,tt,atl,wk1,plcl &
                 , cumtop,clds,cvclds,fluxuu,sdpbl,ncld,qtr        &
                 , clwc,ciwc,cre,cde,cldmaxv )
!
!     compute composite total cloud fraction
!     (implicitly assumes probability of stable and convective
!      clouds are independent to each other)
!
      do 150 k=ls+1,lev
      do 150 i=1,nxj
      clds(i,k)= clds(i,k)+cvclds(i,k)-clds(i,k)*cvclds(i,k)
  150 continue
      do i = 1, nxj
       cufrac(i) = wk1(i)
      end do
!
!  compute total cloudiness
!
!  if "cufrac(i)" .gt. "cld_maxran", assumes maximum-random overlapping
!  if "cufrac(i)" .le. "cld_maxran", assumes random overlapping
!
!  cld_maxran=0. means using max-ran overlapping all the time
!
      cld_maxran = 0.  ! must match with the value used by radsl99
!                         and radir99
      vmin = 1.e-8
      lhig = 37
      lmid = 49
!
      do i = 1, nxj
      cldb(i,1) = 1.0 - clds(i,1)
      ctot(i) = 1.0 - clds(i,1)
      end do
!
      do k = 1, lev-1
      do i = 1, nxj
      if ( cufrac(i) .ge. cld_maxran ) then
      cldb(i,k+1) = max( vmin, 1.0-max( clds(i,k+1), clds(i,k) ) ) /  &
                    max( vmin, 1.0-clds(i,k) )
      else
      cldb(i,k+1) = 1.0 - clds(i,k+1)
      end if
      end do
      end do
!
      do k = 2, lev
      do i = 1, nxj
      ctot(i) = ctot(i) * cldb(i,k)
      end do
!
!  calculate high, middle and low cloudiness
!
! high cloudiness
!
      if( k .eq. lhig )then   ! high cloud above 450hPa
      do i = 1, nxj
       chig(i) = 1.0 - ctot(i)
      end do
      endif
!
      end do
!
! total cloudiness
!
      do i = 1, nxj
      ctot(i) = 1.0 - ctot(i)
      end do
!
! mid cloudiness
!
      do i = 1, nxj
      cldb(i,lhig+1) = 1.0 - clds(i,lhig+1)
      cmid(i) = 1.0 - clds(i,lhig+1)
      end do
!
      do k = lhig+1, lmid-1
      do i = 1, nxj
      if ( cufrac(i) .ge. cld_maxran ) then
      cldb(i,k+1) = max( vmin, 1.0-max( clds(i,k+1), clds(i,k) ) ) /  &
                    max( vmin, 1.0-clds(i,k) )
      else
      cldb(i,k+1) = 1.0 - clds(i,k+1)
      end if
      end do
      end do
!
      do k = lhig+2, lmid
      do i = 1, nxj
      cmid(i) = cmid(i) * cldb(i,k)
      end do
      end do
!
      do i = 1, nxj
      cmid(i) = 1.0 - cmid(i)
      end do
!
! low cloudiness
!
      do i = 1, nxj
      cldb(i,lmid+1) = 1.0 - clds(i,lmid+1)
      clow(i) = 1.0 - clds(i,lmid+1)
      end do
!
      do k = lmid+1, lev-1
      do i = 1, nxj
      if ( cufrac(i) .ge. cld_maxran ) then
      cldb(i,k+1) = max( vmin, 1.0-max( clds(i,k+1), clds(i,k) ) ) /  &
                    max( vmin, 1.0-clds(i,k) )
      else
      cldb(i,k+1) = 1.0 - clds(i,k+1)
      end if
      end do
      end do
!
      do k = lmid+2, lev
      do i = 1, nxj
      clow(i) = clow(i) * cldb(i,k)
      end do
      end do
!
      do i = 1, nxj
      clow(i) = 1.0 - clow(i)
      end do
!
!ccccccccccccccccccccccccccccccccccccccc
!                                      c
!     short wave (solar) radiation     c
!                                      c
!ccccccccccccccccccccccccccccccccccccccc
!
!     (note:  input field 'plcl' is used to store 'dswtop')
!
      do i = 1, nxj
       ss(i) = 0.
       plcl(i) = 0.
       sld(i) = 0.
      end do
      do k = 1, lev
      do i = 1, nxj
       asl(i,k) = 0.
      end do
      end do
      do i = 1, nx*(lev+1)
       fusl(i)  = 0.
       fdsl(i)  = 0.
!cmy-------------------------------------------------------------------
       afusl(i,1)  = 0.
       afdsl(i,1)  = 0.
!cmy-------------------------------------------------------------------
      end do
!
      nsol = 0
      do i = 1, nxj
      if ( cosz(i) .ge. 0.01 ) then
           nsol = nsol + 1
           ibs(nsol) = i
      end if
      end do
!
      if ( nsol .eq. 0 ) go to 190
!
      do i = 1, nsol
      u0(i) = cosz(ibs(i))
      cufa(i) = cufrac(ibs(i))
      pu(i) = ptop
      end do
!
      do j = 1, lev
      jk = j*nsol
      jj = jk-nsol
      do i = 1, nsol
      pu(i+jk)    = pu(i+jj) + dsigma(j,1) * pst(ibs(i)) + dsigma(j,2)
      pa(i+jj)    = pl(ibs(i),j)
      ta(i+jj)    = tt(ibs(i),j)
      wa(i+jj)    = qt(ibs(i),j)
      oa(i+jj)    = o3l(ibs(i),j)
      clda(i+jj)  = clds(ibs(i),j)
      clwca(i+jj) = clwc(ibs(i),j)
      crea(i+jj)  = cre(ibs(i),j)
      ciwca(i+jj) = ciwc(ibs(i),j)
      cdea(i+jj)  = cde(ibs(i),j)
      end do
      end do
!
      do nb = 1, ibsl
      do  i = 1, nsol
      as(i,nb) = alb(ibs(i))
      end do
      end do
!
      call radsl99( nsol,nx,lev,fluxcl,s0,u0,as                 &
                  , pu,pa,ta,wa,oa,clda,clwca,crea,ciwca,cdea   &
                  , cufa,fusl,fdsl,fuslr,fdslr )
!
!  compute solar radiation heating rate at each level
!
      nsn = nsol*lev
      do i = 1, nsol
      ss(ibs(i))       = fdsl(i+nsn) - fusl(i+nsn)
      plcl(ibs(i))     = fdsl(i) - fusl(i)
! noah
      sld(ibs(i))       = fdsl(i+nsn)
!
      end do
!
      do j = 1, lev
      jk = j*nsol
      jj = jk-nsol
      do i = 1, nsol
      xxx = fac / ( dsigma(j,1)*pst(ibs(i))+dsigma(j,2) )
      asl(ibs(i),j) = xxx * ( (fdsl(jj+i)-fusl(jj+i)) -      &
                         (fdsl(jk+i)-fusl(jk+i)) )
      end do
      end do
!cmy-------------------------------------------------------------------
      do j = 1, lev+1
      jk = j*nsol
      jj = jk-nsol
      do i = 1, nsol
      afusl(ibs(i),j) = fusl(jj+i)
      afdsl(ibs(i),j) = fdsl(jj+i)
      end do
      end do
!cmy-------------------------------------------------------------------
!
      if( fluxcl ) then
       do i = 1, nsol
        ss_clr(ibs(i))   = fdslr(i+nsn) - fuslr(i+nsn)
        asol_clr(ibs(i)) = fdslr(i) - fuslr(i)
        sld_clr(ibs(i))  = fdslr(i+nsn)
       end do
!
       do j = 1, lev
       jk = j*nsol
       jj = jk-nsol
       do i = 1, nsol
        xxx = fac / ( dsigma(j,1)*pst(ibs(i))+dsigma(j,2) )
        asl_clr(ibs(i),j) = xxx * ( (fdslr(jj+i)-fuslr(jj+i)) -  &
                            (fdslr(jk+i)-fuslr(jk+i)) )
       end do
       end do
!cmy-------------------------------------------------------------------
      do j = 1, lev+1
      jk = j*nsol
      jj = jk-nsol
      do i = 1, nsol
      afuslr(ibs(i),j) = fuslr(jj+i)
      afdslr(ibs(i),j) = fdslr(jj+i)
      end do
      end do
!cmy-------------------------------------------------------------------

      end if  ! fluxcl=.true.
!
  190 continue
!
!ccccccccccccccccccccccccccccccccccccccc
!                                      c
!     long wave radiation              c
!                                      c
!ccccccccccccccccccccccccccccccccccccccc
!
!
!  assume the ground subject to blackbody, independent of infrared
!  spectrum bands
!
      do nb = 1, ibir
      do  i = 1, lvlw
       ee(i,nb) = 1.0
      end do
      end do
!
      if ( lwfull ) then
!
      do i = 1, nxj
       tga(i) = tg(i)
       pu(i)  = ptop
       cufa(i)= cufrac(i)
      end do
!
      do j = 1, lev
       jk = j*nxj
       jj = jk-nxj
      do i = 1, nxj
       pu(i+jk) = pu(i+jj) + dsigma(j,1) * pst(i) + dsigma(j,2)
       pa(i+jj) = pl(i,j)
       ta(i+jj) = tt(i,j)
       wa(i+jj) = qt(i,j)
       oa(i+jj) = o3l(i,j)
       clda(i+jj)  = clds(i,j)
       clwca(i+jj) = clwc(i,j)
       crea(i+jj)  = cre(i,j)
       ciwca(i+jj) = ciwc(i,j)
       cdea(i+jj)  = cde(i,j)
      end do
      end do
!
      call radir99 ( nxj,nx,lev,fluxcl,ee,pu                      &
                  ,pa,ta,tga,wa,oa,umco2,umch4,umn2o,clda,clwca   &
                  ,crea,ciwca,cdea,cufa,fuir,fdir,fuirr,fdirr )
!
      else
!
      do i = 1, lvlw
       tga(i) = tg(ib(i))
       pu(i)  = ptop
       cufa(i) = cufrac(ib(i))
      end do
!
      do j = 1, lev
       jk = j*lvlw
       jj = jk-lvlw
      do i = 1, lvlw
       pu(i+jk) = pu(i+jj) + dsigma(j,1) * pst(ib(i)) + dsigma(j,2)
       pa(i+jj) = pl(ib(i),j)
       ta(i+jj) = tt(ib(i),j)
       wa(i+jj) = qt(ib(i),j)
       oa(i+jj) = o3l(ib(i),j)
       clda(i+jj)  = clds(ib(i),j)
       clwca(i+jj) = clwc(ib(i),j)
       crea(i+jj)  = cre(ib(i),j)
       ciwca(i+jj) = ciwc(ib(i),j)
       cdea(i+jj)  = cde(ib(i),j)
      end do
      end do
!
      call radir99 ( lvlw,nx,lev,fluxcl,ee,pu                    &
                  ,pa,ta,tga,wa,oa,umco2,umch4,umn2o,clda,clwca  &
                  ,crea,ciwca,cdea,cufa,fuir,fdir,fuirr,fdirr )
!
      end if

!
!  compute infrared radiation heating rate at each level
!

!----
      if (lwfull) then
!----
      nxjl= nxj*lev
      do i = 1, nxj
       cumtop(i) = fuir(i) - fdir(i)
       rs(i)  = fuir(i+nxjl) - fdir(i+nxjl)
       rld(i) = fdir(i+nxjl)
      end do
!
      do j = 1, lev
       jk = j*nxj
       jj = jk-nxj
      do i = 1, nxj
       xxx = fac / ( dsigma(j,1)*pst(i)+dsigma(j,2) )
       atl(i,j) = xxx * ( (fdir(i+jj)-fuir(i+jj)) -       &
                 (fdir(i+jk)-fuir(i+jk)) )
      end do
      end do
!
!cmy--------------------------------------------------------------------
      do j = 1, lev+1
       jk = j*nxj
       jj = jk-nxj
      do i = 1, nxj
       afuir(i,j) = fuir(i+jj)
       afdir(i,j) = fdir(i+jj)
      end do
      end do
!cmy--------------------------------------------------------------------

      if ( fluxcl) then
!
      do i = 1, nxj
       olr_clr(i) = fuirr(i) - fdirr(i)
       rs_clr(i)  = fuirr(i+nxlev) - fdirr(i+nxlev)
       rld_clr(i) = fdirr(i+nxjl)
      end do
!
      do j = 1, lev
       jk = j*nxj
       jj = jk-nxj
      do i = 1, nxj
       xxx = fac / ( dsigma(j,1)*pst(i)+dsigma(j,2) )
       atl_clr(i,j) = xxx * ( (fdirr(i+jj)-fuirr(i+jj)) -   &
                     (fdirr(i+jk)-fuirr(i+jk)) )
      end do
      end do
!
!cmy-------------------------------------------------------------------
!
      do j = 1, lev+1
       jk = j*nxj
       jj = jk-nxj
      do i = 1, nxj
       afuirr(i,j) = fuirr(i+jj)
       afdirr(i,j) = fdirr(i+jj)
      end do
      end do
!
!cmy-------------------------------------------------------------------

      end if

!----
      else
!----

!
!     expend long wave radiation calculation to full grids
!

      nrn = lvlw * lev
      do i = 1, lvlw
       wk1(i) = fuir(i) - fdir(i)
       wk2(i) = fuir(i+nrn) - fdir(i+nrn)
      end do
      call splin1 ( wk1, cumtop, lvlw, nxj, il, cof, 1.0 )
      call splin1 ( wk2, rs, lvlw, nxj, il, cof, 1.0 )
!
      do i = 1, lvlw
       wk1(i) = fdir(i+nrn)
      end do
      call splin1 ( wk1, rld, lvlw, nxj, il, cof, 1.0 )
!
      do j = 1, lev
       jk = j*lvlw
       jj = jk-lvlw
      do i = 1, lvlw
       xxx = fac / ( dsigma(j,1)*pst(ib(i))+dsigma(j,2) )
       wk1(i) = xxx * ( (fdir(i+jj)-fuir(i+jj)) -        &
               (fdir(i+jk)-fuir(i+jk)) )
      end do
      call splin1 ( wk1, atl(1,j), lvlw, nxj, il, cof, 1.0 )
      end do
!
!cmy-------------------------------------------------------------------
      do j = 1, lev+1
       jk = j*lvlw
       jj = jk-lvlw
      do i = 1, lvlw
       wk1(i) = fuir(i+jj)
       wk2(i) = fdir(i+jj)
      end do
      call splin1_sp ( wk1, afuir(1,j), lvlw, nxj, il, cof, 1.0 )
      call splin1_sp ( wk2, afdir(1,j), lvlw, nxj, il, cof, 1.0 )
      end do
!cmy-------------------------------------------------------------------

      if ( fluxcl ) then
!
      do i = 1, lvlw
       wk1(i) = fuirr(i) - fdirr(i)
       wk2(i) = fuirr(i+nrn) - fdirr(i+nrn)
      end do
      call splin1 ( wk1, olr_clr, lvlw, nxj, il, cof, 1.0 )
      call splin1 ( wk2, rs_clr, lvlw, nxj, il, cof, 1.0 )
!
      do j = 1, lev
       jk = j*lvlw
       jj = jk-lvlw
      do i = 1, lvlw
       xxx = fac / ( dsigma(j,1)*pst(ib(i))+dsigma(j,2) )
       wk1(i) = xxx * ( (fdirr(i+jj)-fuirr(i+jj)) -       &
               (fdirr(i+jk)-fuirr(i+jk)) )
      end do
      call splin1 ( wk1, atl_clr(1,j), lvlw, nxj, il, cof, 1.0 )
      end do
!
!cmy-------------------------------------------------------------------
      do j = 1, lev+1
       jk = j*lvlw
       jj = jk-lvlw
      do i = 1, lvlw
       wk1(i) = fuirr(i+jj)
       wk2(i) = fdirr(i+jj)
      end do
      call splin1_sp ( wk1, afuirr(1,j), lvlw, nxj, il, cof, 1.0 )
      call splin1_sp ( wk2, afdirr(1,j), lvlw, nxj, il, cof, 1.0 )
      end do
!cmy-------------------------------------------------------------------

      end if

!----
      endif
!----

!
!     compute total heating rate: dtrad
!     and average/max values for diagnostics
!
      do k = 1, lev
      do i = 1, nxj
       dtrad(i,k)     = asl(i,k) + atl(i,k)
      end do
      end do
!
      xim = 1.0 / float(nxj)
!
      do k = 1, lev
       xsr(k) = asl(1,k)
       xlr(k) = atl(1,k)
       asr(k) = asl(1,k)
       alr(k) = atl(1,k)
       asr_clr(k) = asl_clr(1,k)
       alr_clr(k) = atl_clr(1,k)
       acld(k)= clds(1,k)
      end do
      do k = 1, lev
      do i = 2, nxj
       xsr(k) = max ( xsr(k),asl(i,k) )
       xlr(k) = min ( xlr(k),atl(i,k) )
       asr(k) = asr(k) + asl(i,k)
       alr(k) = alr(k) + atl(i,k)
       asr_clr(k) = asr_clr(k) + asl_clr(i,k)
       alr_clr(k) = alr_clr(k) + atl_clr(i,k)
       acld(k)= acld(k)+ clds(i,k)
      end do
      end do
      do k = 1, lev
       asr(k) = asr(k) * xim
       alr(k) = alr(k) * xim
       asr_clr(k) = asr_clr(k) * xim
       alr_clr(k) = alr_clr(k) * xim
       acld(k)= acld(k)* xim
      end do
!
      return
      end
