      subroutine diagrd (ozon,nx,lev,my,my_max,njump,julian,stbo,s0,cosl &
              ,tt,qt,plcl,cumtop,ss,rs,asr,alr,xsr,xlr,acld,aflxd,aflxu  &
              ,ncld)
      use mpe
      use rank
      use index
!######################################################################c
!                                                                      c
! a. function:                                                         c
!  print radiation diagnoses when raddia=.true.                        c
!                                                                      c
! b. common block specification:                                       c
!  none                                                                c
!                                                                      c
! c. parameter specification:                                          c
!  none                                                                c
!                                                                      c
! d. input/output variables                                            c
!  (1) input variables:                                                c
!  ozon  : logical constant to include ozon in sw/lw rad calculation   c
!  nx    : i-direction horizontal dimension of the model               c
!  lev    : vertical layers of the model                               c
!  my    : j-direction horizontal dimension of the model               c
!  njump : grid interval for long wave radiative calculation           c
!  julian: julian day ( julian = 1 for january 1 )                     c
!  stbo  : stefan-boltzman constant (5.669e-8)                         c
!  s0    : solar constant (1368.3)                                     c
!  cosl  : cos(lat)                                                    c
!  tt(nx,lev,my) : real temperature at prediction levels (k)           c
!  qt(nx,lev,my) : specific humidity at prediction levels (kg/kg)      c
!  plcl(nx,my)  : total column solar radiation absorption              c
!  cumtop(nx,my): outgoing long wave radiation at model top            c
!  ss(nx,my): net short wave flux down to ground surface (wat/m**2)    c
!  rs(nx,my): net long  wave flux up from ground surface (wat/m**2)    c
!  asr(lev,my): lat-mean temp change due to short wave radtn (k/day)   c
!  alr(lev,my): lat-maen temp change due to long  wave radtn (k/day)   c
!  xsr(lev,my): maxnxum short wave heating for each lat. ring (k/day)  c
!  xlr(lev,my): maxnxum long  wave colling for each lat. ring (k/day)  c
!  acld(lev,my): lat-mean total cloud fraction assuming random overlap c
!  aflxd(lev+2,my): lat-mean downward long wave flux (wat/m**2)        c
!  aflxu(lev+2,my): lat-mean upward   long wave flux (wat/m**2)        c
!                                                                      c
! e. calling modules:                                                  c
!    diabat                                                            c
!                                                                      c
! f. usage:                                                            c
!     call diagrd ( ozon,nx,lev,my,njump,julian,stbo,s0,cosl,tt,qt     c
!    1       , plcl,cumtop,ss,rs,asr,alr,xsr,xlr,acld,aflxd,aflxu )    c
!                                                                      c
! g. modules called:                                                   c
!  none                                                                c
!                                                                      c
! h. lnxitation:                                                       c
!  no                                                                  c
!                                                                      c
! i. date:                                                             c
!  created: oct 5 , 1992                                               c
!                                                                      c
! j. author:                                                           c
!  c.-s. liou, c. fong                                                 c
!  modify to f90 by C-H Lee and sort by River Chen in 2015
!                                                                      c
!######################################################################c
!
!      use index, only : nxp

      implicit  none

      integer   nx,lev,my,my_max,njump,julian,ncld

      real      stbo,s0

      logical ozon
!
      real      cosl(my),tt(nxp,lev,my_max),qt(nxp,lev*ncld,my_max),&
                plcl(nxp,my_max),cumtop(nxp,my_max),ss(nxp,my_max),rs(nxp,my_max),    &
                asr(lev,my),alr(lev,my),xsr(lev,my),xlr(lev,my),  &
                acld(lev,my),aflxd(lev+2,my),aflxu(lev+2,my)

      logical flag
!

      integer   i,j,k,jj,nxj,levp1
      real      xmy,xstop,xltop,xssfc,xlsfc,zfluxu,zfluxd,zasr,zalr, &
                zclds,ztt,zqt,xmsr,xnlr,zttj,zqtj,ztadd
!
      if(myrank .eq. 0) print 900, ozon,nx,lev,my,njump,julian,stbo,s0
  900 format (1h0,15x,'-------layer-mean radtn diagnosis : ---------'  &
             , //,2x, 'ozon =',l2,';  nx,lev,my =',3i4                 &
             , '; njump =', i3,';  julian =', i3                       &
             , ';  stbo,s0 =', 2e14.6 )
!
      levp1 = lev + 1
      xmy = 1./float(my)
!
      call glbmean_2d ( cosl,plcl,  xstop )
      call glbmean_2d ( cosl,cumtop,xltop )
      call glbmean_2d ( cosl,ss, xssfc )
      call glbmean_2d ( cosl,rs, xlsfc )
!
      if(myrank .eq. 0) print 910, xstop, xltop, xssfc, xlsfc
  910 format( 2x, 'top:  down-sw, upward-lw flux =', 2f7.2            &
            , ';   sfc:  down-sw, upward-lw flux =', 2f7.2, //        &
            , 4x, '   asl     atl     tadd    clds     tt      qt  '  &
            ,     '  mxasl   mnatl   fluxu   fluxd ')
!
      xmy = 0.0
      do 100 j = 1, my
      xmy = xmy + cosl(j)
  100 continue
      xmy = 1.0/xmy
!     ximy= xmy/float(nx)
!
      zfluxu = 0.0
      zfluxd = 0.0
      do 810 jj  = 1, jlistnum
      j=jlist1(jj)
      zfluxu = zfluxu + aflxu(1,j)*cosl(j)
      zfluxd = zfluxd + aflxd(1,j)*cosl(j)
  810 continue
!
      call mpe_global_sum(zfluxu,1,mpe_double)
      call mpe_global_sum(zfluxd,1,mpe_double)
!
      zfluxu = zfluxu * xmy
      zfluxd = zfluxd * xmy
      if(myrank .eq. 0) print 915, zfluxu,zfluxd
  915 format( 67x, 2f8.2 )
      do 890 k = 1, lev
      zasr = 0.0
      zalr = 0.0
      zclds = 0.0
      ztt = 0.0
      zqt = 0.0
      zfluxu = 0.0
      zfluxd = 0.0
!
!     xmsr = xsr(k,1)
!     xnlr = xlr(k,1)
      xmsr =-1.797693134862316d+308
      xnlr = 1.797693134862316d+308
!
      do 840 jj  = 1, jlistnum
      j=jlist1(jj)
      zasr  = zasr  + asr(k,j)*cosl(j)
      zalr  = zalr  + alr(k,j)*cosl(j)
      zclds = zclds + acld(k,j)*cosl(j)
      zfluxu = zfluxu + aflxu(k+1,j)*cosl(j)
      zfluxd = zfluxd + aflxd(k+1,j)*cosl(j)
      xmsr = max( xmsr, xsr(k,j) )
      xnlr = min( xnlr, xlr(k,j) )
      zttj = 0.0
      zqtj = 0.0
      nxj=nxdef(j)
!     do 845 i = 1,   nxj
      do 845 i = 1, nxdef_2d(j)
      zttj  = zttj  + tt(i,k,jj)*cosl(j)
      zqtj  = zqtj  + qt(i,k,jj)*cosl(j)
  845 continue
      ztt = ztt + zttj/float(nxj)
      zqt = zqt + zqtj/float(nxj)
  840 continue
!
      call mpe_global_sum(zasr,1,mpe_double)
      call mpe_global_sum(zalr,1,mpe_double)
      call mpe_global_sum(zclds,1,mpe_double)
      call mpe_global_sum(zfluxu,1,mpe_double)
      call mpe_global_sum(zfluxd,1,mpe_double)
      call mpe_global_sum(ztt,1,mpe_double)
      call mpe_global_sum(zqt,1,mpe_double)
      call mpe_global_max(xmsr,1,mpe_double)
      call mpe_global_min(xnlr,1,mpe_double)
!
      zasr = zasr * xmy
      zalr = zalr * xmy
      ztadd = zasr + zalr
      zclds = zclds * 100.0 * xmy
!x    ztt = ztt * ximy
!x    zqt = zqt * ximy * 1000.0
      ztt = ztt * xmy
      zqt = zqt * xmy * 1000.0
      zfluxu = zfluxu * xmy
      zfluxd = zfluxd * xmy
      if(myrank .eq. 0) print 920,k,zasr,zalr,ztadd,zclds,ztt,zqt,xmsr, &
                        xnlr,zfluxu,zfluxd
  920 format ( 1x, i2, 10f8.2 )
  890 continue
!
      zfluxu = 0.0
      zfluxd = 0.0
      do 895 jj  = 1, jlistnum
      j=jlist1(jj)
      zfluxu = zfluxu + aflxu(lev+2,j)*cosl(j)
      zfluxd = zfluxd + aflxd(lev+2,j)*cosl(j)
  895 continue
!
      call mpe_global_sum(zfluxu,1,mpe_double)
      call mpe_global_sum(zfluxd,1,mpe_double)
!
      zfluxu = zfluxu * xmy
      zfluxd = zfluxd * xmy
      if(myrank .eq. 0) print 915, zfluxu,zfluxd
!
      return
      end
