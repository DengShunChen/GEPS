      subroutine cupcwb ( jcup,nxj,nx,my,lev,ktcup,dt,g,r,cp,hltm,etop &
                        , prevap,topo,ps,pl,pk,pk2,tl,ql,phi,plcl      &
                        , cumtop,rcup,ncup,ptop,dsigma,nxx,dtcupx      &
                        , dtcupz,dqcupz,idg,dtcupd,dqcupd              &
                        , nlcl,nnegl,nosat,nwork,ntcup,nflx )
!
!***********************************************************************
!
!  1.   function description
!         perform cumulus parameterization
!
!  2.   common block specification
!          none
!
!  3.   parameter specification
!   input :
!     jcup      :   latitude index
!     nx        :   no of points in logitude direction
!     my        :   no of latitude rings
!     lev       :   vertical levels
!     ktcup     :   highest sigma level allowed for the top of
!                   highest cloud type
!     dt        :   time step used in this cup call
!     g         :   gravity (9.806)
!     r         :   gas constant (287.0)
!     cp        :   heat specific constant of the air (=1004.)
!     hltm      :   latent heat constant of water vapor (2.52e+6)
!     etop      :   coeff for cloud top enhanced entrain
!     prevap    :   fraction of falling precip evaporated
!     ptop      :   pressure at model top           (mb)
!     dsigma    :   sigma layer thickness (between even levels)
!     topo(nx)  :   terrain geopotential (m2/s2)
!     ps(nx)    :   terrain pressure - ptop (mb)
!     pl(nx,lev) :   pressure on model level
!     pk(nx,lev) :   p**capa on model levels
!     pk2(nx,lev):   p**capa on even levels
!     tl(nx,lev) :   temperature on model level
!     ql(nx,lev) :   water vapor mixing ratio on model level
!     phi(nx,lev):   odd level geopotential
!     idg       :   i-index address of selected point to print profile
!   output:
!     plcl(nx)  :   pressure value on cloud base (mb)
!     cumtop(nx):   pressure value on cloud top  (mb)
!     rcup(nx)  :   precipitation  due to cumulus convection (mm/call)
!     ncup      :   no of convective points in this call
!     nxx       :   2-d address index of the max heating point
!     dtcupx    :   max cup heating rate occurred in this call (k/day)
!     dtcupz(lev):   horiz. average of cup heating rate         (k/day)
!     dqcupz(lev):   horiz. average of cup moistening rate      (1/day)
!     dtcupd(lev):   cup heating    rate at selected point idg  (k/day)
!     dqcupd(lev):   cup moistening rate at selected point idg  (1/day)
!     nlcl(lev)  : no of clouds at level k with tops above lcl
!     nnegl(lev) : no of clouds at level k with positive lamda
!     nosat(lev) : no of clouds at level k with liquid water at top
!     nwork(lev) : no of clouds at level k with work func above ref
!     ntcup(lev) : no of clouds at level k with non-zero cbase masflx
!     nflx(lev)  : no of clouds exceeding mass flux limit
!
!   local work arrays
!
!     qls(im,lm):   saturation mixing ratio on model level
!     pl2(im,lm):   pressure on even level
!     pmassl(im,lm): mass in between model even sigma levels  (kg)
!                  = ps(i)*dsig(l)*100.0/g
!     dtcup(im,lm): temp. change due to cumulus convection    (k/call)
!     dqcup(im,lm): moisture change due to cumulus convection (1/call)
!
!  4.  calling modules
!          diabat
!
!  5.  usage
!
!        call cupcwb (jcup,nx,lev,ktcup,dt,g,r,cp,hltm,etop,prevap,ptop
!       1            ,dsigma,topo,ps,pl,pk,pk2,tl,ql,plcl,cumtop,rcup
!       2            ,ncup,imx,dtcupx,dtcupz,dqcupz,idg,dtcupd,dqcupd
!       3            ,nlcl,nnegl,nosat,nwork,ntcup,nflx )
!
!  6.  modules called
!          creeve , cup92
!
!  7.  date
!        created    by   c-s chen    ( cwb )    1992
!        modified   by   c-s liou and t rosmond 1993
!        modify to f90 by C-H Lee and sort by River Chen in 2015
!
!***********************************************************************
!
      use paramt
      use const, only : RTYPE

      implicit  none

      integer   jcup,nxj,nx,my,lev,ktcup,idg,nxx,ncup

      real      pl(nx,lev),tltmp(nx)                         &
        ,plcl(nx),cumtop(nx),rcup(nx)                        &
        ,dtcupz(lev),dqcupz(lev),dtcupd(lev),dqcupd(lev)
      real(kind=RTYPE) tl(nx,lev),ql(nx,lev),phi(nx,lev)     &
        ,       topo(nx),ps(nx),dsigma(lev,2),pk(nx,lev)     &
        ,       pk2(nx,lev)

      integer   nlcl(lev),nnegl(lev),nosat(lev),nwork(lev),ntcup(lev),nflx(lev)
!
!     local work arrays
!
      logical cexist(im)
      real      pl2(im,lm),tl2(im,lm),ql2(im,lm),pkc(im*lm)     &
       , pkc2(im*lm),qls(im,lm),zl(im,lm),ssl(im,lm),hhl(im,lm) &
       , qls2(im,lm),zl2(im,lm),phil2(im,lm),ssl2(im,lm)        &
       , hhl2(im,lm),hhls(im,lm),hhls2(im,lm),pmassl(im,lm)     &
       , gamwe(im,lm),plc(im*lm),tlc(im*lm),qlc(im,lm)          &
       , hc(im,lm),qc(im,lm),cu(im,lm),e(im,lm),flxmas(im,lm)   &
       , flxmasc(im*lm),dtcupc(im*lm),dqcupc(im*lm),rcupc(im)
!
      integer   kcbase(im),kctop(im)

      real      dtcup(im,lm),dqcup(im,lm),hbase(im),qbase(im),cubase(im)

      integer   nxlev,i,ii,l,k,lt,lb,m,n,nxx0,ilsum,ktmin
      real      ocp,xkapa,xnx,fac,qmin,const1
      real      dtcupx,ptop,prevap,etop,hltm,cp,r,g
      real      dt,gamfac,sige,sigd,dtcupx0
!
      e(1:im,1:lm)=0.0
!
      nxlev = nx*lev
      ocp  = 1.0/cp
      xkapa = 1.0/3.5
      xnx   = 1.0/float(nxj)
      fac   = 86400.0 / dt
      qmin  = 1.0e-10
      const1= cp*0.608/hltm
      gamfac= hltm*5417.9827*ocp

!CWB2015
      pl2=0.

!
!     define even level variables and mass between even levels (pmassl)
!
      do k = 1, lev
      tltmp(:)=tl(:,k)
      call qsatq (nxj,tltmp,pl(1,k),qls(1,k))
      enddo
!
      sige = 0.0
      sigd = 0.0
      do 20 l = 1, lev
        sige = sige + dsigma(l,1)
        sigd = sigd + dsigma(l,2)
      do 20 i = 1, nxj
        pl2(i,l) = ps(i)*sige +sigd + ptop
        pmassl(i,l) = (ps(i)*dsigma(l,1)+dsigma(l,2))*100.0/g
        qlc(i,l)= max(qmin, min(ql(i,l),qls(i,l)*0.999))
  20  continue
!
      call creeve ( nxj,nx,lev,g,r,cp,hltm,pl,pl2,tl,qlc,topo,pk,pk2 &
                  , qls,ssl,hhl,hhls,zl,phi,tl2,ql2,qls2,ssl2,hhl2   &
                  , hhls2,zl2,phil2 )
!
      do 12 k = 1,lev
      do 12 i = 1,nxj
        gamwe(i,k) = gamfac*qls2(i,k)/(tl2(i,k)*tl2(i,k))
        qc(i,k)= 0.0
        hc(i,k)= 0.0
        cu(i,k)= 0.0
   12 continue
!
!     call dcbase to determine cloud base level and cloud base condition
!
      call dcbase ( nxj,nx,lev,hltm,kcbase,zl2,qlc,qls2,hhl,hhls2,gamwe &
                  , e,qc,hc,cu,hbase,qbase,cubase )
!
! --- find the deepest cloud type index
!

!     do 14 i = 1, nxj
!     kctop(i)= 0
!  14 continue

!CWB2015
      kctop= 1
!
      ktmin = lev
      do 1012 l = lev-1, ktcup, -1
      do 1012 i = 1,nxj
!CWB2015
!     if ((hbase(i).gt.hhls2(i,l)).and.(hbase(i).le.hhls2(i,l-1)) &
!        .and.( kctop(i) .eq. 0 )) then
      if ((hbase(i).gt.hhls2(i,l)).and.(hbase(i).le.hhls2(i,l-1)) &
         .and.( kctop(i) .eq. 1 )) then
      kctop(i) = l
      ktmin= min(ktmin,l)
      endif
 1012 continue
!
      do 1014 i= 1, nxj
!CWB2015
!     cexist(i)= kctop(i).ne.0.and.kcbase(i).ne.0
      cexist(i)= kctop(i).ne.1.and.kcbase(i).ne.1
      cexist(i)= cexist(i).and.(kctop(i).le.kcbase(i))
 1014 continue
!
      do i = 1, nxj
       plcl(i) = 0.
       cumtop(i) = 0.
       rcup(i) = 0.
      enddo
      do k = 1, lev
      do i = 1, nxj
       dtcup(i,k) = 0.
       dqcup(i,k) = 0.
      enddo 
      enddo
!
      n= ilsum( nxj,cexist,1 )
      ncup= n
      if (n.eq.0) go to 250
!
!  compress thermodynamic arrays by eliminating those grid
!  points where neither a cloud top or cloud base were found
!
!      do 500 k = ktmin-1,lev
      do 500 k = 1,lev
      m = nx*(k-1)
!dir$ ivdep
!ocl novrec
      do 500 i = 1, nxj
      if (cexist(i)) then
      m = m+1
      plc(m)     = pl(i,k)
      pkc(m)     = pk(i,k)
      pkc2(m)    = pk2(i,k)
      tlc(m)     = tl(i,k)
      qlc(m,1)   = qlc(i,k)
      zl(m,1)    = zl(i,k)
      qls(m,1)   = qls(i,k)
      ssl(m,1)   = ssl(i,k)
      hhls(m,1)  = hhls(i,k)
      pl2(m,1)   = pl2(i,k)
      zl2(m,1)   = zl2(i,k)
      tl2(m,1)   = tl2(i,k)
      ql2(m,1)   = ql2(i,k)
      qls2(m,1)  = qls2(i,k)
      ssl2(m,1)  = ssl2(i,k)
      hhl2(m,1)  = hhl2(i,k)
      hhls2(m,1) = hhls2(i,k)
      pmassl(m,1)= pmassl(i,k)
      hc(m,1)    = hc(i,k)
      qc(m,1)    = qc(i,k)
      cu(m,1)    = cu(i,k)
      e(m,1)     = e(i,k)
      endif
  500 continue
!
      m = 0
!dir$ ivdep
!ocl novrec
      do 510 i = 1, nxj
      if (cexist(i)) then
      m = m+1
      kctop(m) = kctop(i)
      kcbase(m)= kcbase(i)
      hbase(m) = hbase(i)
      qbase(m) = qbase(i)
      cubase(m)= cubase(i)
      endif
  510 continue
!
!     compute cumulus heating and moistening rate using a-s param.
!
      call cup92 ( n,nx,lev,ktmin,dt,g,r,cp,hltm,etop,prevap,plc &
                 , pkc,zl,tlc,qlc,qls,ssl,hhl,hhls,pl2,pkc2      &
                 , zl2,tl2,ql2,qls2,ssl2,hhl2,hhls2,pmassl       &
                 , hbase,qbase,cubase,e,hc,qc,cu,kcbase,kctop    &
                 , dtcupc,dqcupc,flxmasc,rcupc,ncup,nlcl         &
                 , nnegl,nosat,nwork,ntcup,nflx )
!
!     add heating and moistening rate to tl and ql
!
!  expand compressed cumulus change arrays back into full model
!  arrays
!
      do 140 k = ktmin, lev
      m = nx*(k-1)
      do 140 i = 1, nxj
      if (cexist(i)) then
      m = m+1
      dtcup(i,k) = dtcupc(m)
      dqcup(i,k) = dqcupc(m)
      flxmas(i,k)= flxmasc(m)
      endif
  140 continue
!
!     determine cloud base (plcl) and cloud top (cumtop) pressure
!
      m = 0
      do 200 i = 1,nxj
      if ( cexist(i) ) then
      m = m+1
      rcup(i) = rcupc(m)
      lb = kcbase(m)
      lt = kctop(m)
        plcl(i)  = pl2(m,lb)
        cumtop(i)= pl(i,lt)
      endif
  200 continue
!
!  don't allow cumlus to create negative moisture, conserve
!  moist static energy
!
      ii = nx*(ktmin-1) + 1
      do 145 k = ktmin, lev
      do 145 i = 1, nxj
      dqcup(i,k)= max(dqcup(i,k),-(ql(i,k)-qmin))
      dtcup(i,k)= ocp*(dtcup(i,k)-hltm*dqcup(i,k))
  145 continue
!
      do 150 k = ktmin, lev
      do 150 i = 1, nxj
      tl(i,k) = tl(i,k) + dtcup(i,k)
      ql(i,k) = ql(i,k) + dqcup(i,k)
  150 continue
!
!     compute max heating rate of this call  (k/day)
!
  250 continue
!x    call maxp ( nxlev,dtcup,dtcupx,nxx )
      dtcupx=dtcup(1,1)
      nxx=1
      do k=1,lev
      call maxp ( nxj,dtcup(1,k),dtcupx0,nxx0 )
      if( dtcupx0.gt.dtcupx ) then
        dtcupx = dtcupx0
        nxx = (k-1)*nx+nxx0
      endif
      enddo
      dtcupx = dtcupx * fac
!
!     compute layer maen heating (k/day) and moistening rates
!
      do 350 l = 1, lev
        dtcupz(l) = 0.0
        dqcupz(l) = 0.0
        do 300 i = 1, nxj
          dtcupz(l) = dtcupz(l) + dtcup(i,l)
          dqcupz(l) = dqcupz(l) + dqcup(i,l)
  300   continue
        dtcupz(l) = dtcupz(l)*xnx*fac
        dqcupz(l) = dqcupz(l)*xnx*fac
  350 continue
!
!     save dtcup, dqcup profile to dtcupd, dqcupd; if idg is not 0
!
      if ( idg .ne. 0 )  then
        do 400 l = 1, lev
        dtcupd(l) = dtcup(idg,l)*fac
        dqcupd(l) = dqcup(idg,l)*fac
  400   continue
      endif
!
      return
      end
