      subroutine creeve ( nxj,nx,lev,g,r,cp,hltm,pl,pl2,tl,ql,topo  &
                        , pk,pk2,qls,sl,hl,hls,zl,phil,tl2,ql2      &
                        , qls2,sl2,hl2,hls2,zl2,phil2 )
!
!***********************************************************************
!
!    1.   function description
!           interpolate variable from model level to even level
!           we define the local large-scale environment in
!           terms of temperature, water vapor mixing ratio,
!           saturation water vapor mixing ratio, dry static
!           energy, moist static energy, saturation moist
!           static energy and heights of pressure levels.
!
!    2.   common block specification
!            none
!
!    3.   parameter specification
!     input :
!       nx          :  dimension of horizontal direction
!       lev         :  position index in the latitudinal direction
!       g           :  gravity (9.806)
!       r           :  gas constant (287.0)
!       cp          :  heat specific constant of the air (=1004.)
!       hltm        :  latent heat constant of water vapor (2.52e+6)
!       pl(nx,lev)  :  pressure     on odd level  ( mb  )
!       pl2(nx,lev) :  pressure     on even level ( mb  )
!       tl(nx,lev)  :  temperature  on odd level  ( k   )
!       ql(nx,lev)  :  mixing ration on odd level (kg/kg)
!       topo(nx)    :  geopotential height for toporaphy
!     output:
!       qls(nx,lev)  : saturation mixing ration on odd level (kg/kg)
!       sl(nx,lev)   : dry static energy on model level
!       hl(nx,lev)   : moist static energy on model level
!       hls(nx,lev)  : saturation moist static energy on model level
!       zl(nx,lev)   : geopotential hight on model level
!       phil(nx,lev) : geopotential on model level
!       tl2(nx,lev)  : temperature  on even level  ( k   )
!       ql2(nx,lev)  : mixing ration on even level (kg/kg)
!       qls2(nx,lev) : saturation mixing ration on even level (kg/kg)
!       sl2(nx,lev)  : dry static energy on even level
!       hl2(nx,lev)  : moist static energy on even level
!       hls2(nx,lev) : saturation moist static energy on even level
!       zl2(nx,lev)  : geopotential hight on even level
!       phil2(nx,lev): geopotential on even level
!
!    4.  calling modules
!            cupcwb
!
!    5.  usage
!        call creeve ( nx,lev,g,r,cp,hltm,pl,pl2,tl,ql,topo,pk,pk2
!       1            , qls,sl,hl,hls,zl,phil,tl2,ql2,qls2,sl2,hl2,hls2
!       2            , zl2,phil2 )
!
!    6.  modules called
!            none
!
!    7.  date
!          created    1992           by      c-s chen  ( cwb  )
!          modified   1993           by      c-s liou
!          modify to f90 by C-H Lee and sort by River Chen in 2015
!
!***********************************************************************
!
      use const, only: RTYPE
!
      implicit  none
      integer   nxj,nx,lev,levm1,nxlev,i,l,k,lk,lkp

      real      pl(nx,lev),ql(nx,lev),qls(nx,lev),sl(nx,lev)            &
              , hl(nx,lev),hls(nx,lev),zl(nx,lev)                       &
              , pl2(nx,lev)                                             &
              , tl2(nx,lev),ql2(nx,lev),qls2(nx,lev),sl2(nx,lev)        &
              , hl2(nx,lev),hls2(nx,lev),zl2(nx,lev),phil2(nx,lev)
      real(kind=RTYPE) tl(nx,lev),phil(nx,lev),topo(nx),pk(nx,lev)      &
              , pk2(nx,lev)
      real      tmp1(nx,lev)
!
      real      ograv,g,hltm,cp,r,ratio,amean,bmean

      ograv = 1.0/g
      levm1 = lev-1
      nxlev = nx*lev
!
!     (2) even levels (phil2)
!
      do 100 i = 1, nxj
      phil2(i,lev) = topo(i)
  100 continue
      do 200 l = 1, levm1
      k = lev - l
      do 200 i = 1, nxj
      phil2(i,k) = phil2(i,k+1) + cp*(pk2(i,k+1)-pk2(i,k)) &
                          *tl(i,k+1)*(1.0+0.608*ql(i,k+1))/pk(i,k+1)
  200 continue
      do 220 k = 1, lev
      do 220 i = 1, nxj
      zl(i,k)  = (phil(i,k)-topo(i)) *ograv
      zl2(i,k) = (phil2(i,k)-topo(i))*ograv
  220 continue
!
!     linear-z interpolation temp and log(q) to even levels
!
      do l = 1, levm1
      do i = 1, nxj
        ql2(i,l)=ql(i,l+1)/ql(i,l) 
      enddo
       call vlog(ql2(1,l),ql2(1,l),nxj)
      enddo
      do l = 1, levm1
      do i = 1, nxj
       tmp1(i,l)=(phil2(i,l)-phil(i,l))/(phil(i,l+1)-phil(i,l))
       ql2(i,l)=tmp1(i,l)*ql2(i,l)
      enddo
       call vexp(ql2(1,l),ql2(1,l),nxj)
      enddo

      do 300 l = 1, levm1
      do 300 i = 1, nxj
      ratio = tmp1(i,l)
      tl2(i,l) = tl(i,l) + ratio*(tl(i,l+1)-tl(i,l))
      ql2(i,l) = ql(i,l) * ql2(i,l)
  300 continue
!
      do i=1,nxj
        ql2(i,lev) = ql(i,lev)/ql2(i,levm1)
      enddo
      call vlog(ql2(1,lev),ql2(1,lev),nxj)
      do i=1,nxj
       tmp1(i,lev)=(phil2(i,lev)-phil(i,lev))/(phil(i,lev)-phil2(i,levm1))
       ql2(i,lev)=tmp1(i,lev)*ql2(i,lev)
      enddo
      call vexp(ql2(1,lev),ql2(1,lev),nxj)

      do 320 i = 1, nxj
      ratio =tmp1(i,lev)
      tl2(i,lev) = tl(i,lev) + ratio*(tl(i,lev)-tl2(i,levm1))
      ql2(i,lev) = ql(i,lev)*ql2(i,lev) 
  320 continue
!
!     compute qls2 and ensure ql2 is not supersatuated
!
      do 340 k = 1, lev
      call qsatq ( nxj,tl2(1,k),pl2(1,k),qls2(1,k) )
      do 340 i = 1, nxj
      ql2(i,k) = min(ql2(i,k),qls2(i,k))
  340 continue
!
!     compute sl, sl2 and ensure no static instability
!
      do 400 k = 1, lev
      do 400 i = 1, nxj
      sl(i,k) = cp*tl(i,k)  + phil(i,k)
      sl2(i,k)= cp*tl2(i,k) + phil2(i,k)
  400 continue
      do 410 l = 2, lev
      lk = lev - l + 1
      lkp= lk + 1
      do 415 i = 1, nxj
      if ( sl(i,lk) .lt. sl(i,lkp) )  then
        amean = 0.5*(sl(i,lk)+sl(i,lkp))
        sl(i,lk) = amean
        sl(i,lkp)= amean
      endif
      if ( sl2(i,lk) .lt. sl2(i,lkp) )  then
        bmean = 0.5*(sl2(i,lk)+sl2(i,lkp))
        sl2(i,lk) = bmean
        sl2(i,lkp)= bmean
      endif
  415 continue
  410 continue
!
!     compute hl, hls, hl2 and hls2 at odd and even levels
!
      do 420 k = 1, lev
      do 420 i = 1, nxj
      hls(i,k) = sl(i,k)  + hltm*qls(i,k)
      hl(i,k)  = sl(i,k)  + hltm*ql(i,k)
      hl2(i,k) = sl2(i,k) + hltm*ql2(i,k)
      hls2(i,k)= sl2(i,k) + hltm*qls2(i,k)
  420 continue
!
      return
      end
