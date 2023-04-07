      subroutine lsp ( t,q,pl,pst,dsigma,g,nxj,nx,lev,evaprh,rlsp,cp,hltm &
                     , lsppt,ipass )
!
!       large scale precipitation calculation     c.s. liou
!
!     input variables:  t, q, pl, pms, nx, lev
!       t: real temperature (nx,lev),                            (k)
!       q: specific humidity (nx,lev),                           (kg/kg)
!       pl: pressure on odd sigma levels                        (mb)
!      pst: (terrain pressure -ptop)                            (mb)
!     dsigma: sigma layer thickness
!        g: gravity
!       nx: horizontal dimension,   lev: vertical dimension,
!   evaprh: relative humidity that stops evaporation,
!           (set to zero if no evaporation is allowed)
!       cp: specific heat constant for air
!     hltm: latent heat relase constant for water vapor
!
!     output variables:  t, q, rlsp
!        t: adjusted real temperature after lsp process,        (k)
!        q: adjusted specific humidity after lsp process,       (kg/kg)
!     rlsp: lsp precipitation rate per call.                  (mm/call)
!    lsppt: no. of saturated points at each layer
!    ipass: no. of iterations at each layer for the converge
!
!     by c. s. liou     6/92
!     modify to f90 by C-H Lee and sort by River Chen in 2015
!
!     local work array:
!       pms: mass between even sigma levels= dsigma*pst*100.0/g  (kg)
!       ws,wq,wt,wp,qs,rain,cld,lsppt,ipass
!     external function called:  qsatq
!
      use paramt
      use const, only: RTYPE

      implicit  none

      integer   nxj,nx,lev,inx,ic
      integer   lsppt(lev),ipass(lev)

      real      g,evaprh,cp,hltm 
      real      pl(nx,lev),rlsp(nx),ttmp(nx)
      real(kind=RTYPE) t(nx,lev),q(nx,lev),pst(nx),dsigma(lev,2)
!
!     local work arrays
!
      real      pms(im,lm)
      real      ws(im),wq(im),wt(im),wp(im),qs(im)
      logical   cld(im)
!
      integer   icmax,l,i,nsat,nrain,lp
      real      gamfac,ww,wmax,wev,hlcp

      gamfac = hltm * 5417.983 / cp
      hlcp = hltm / cp
      icmax = 5
!
!     compute mass in between sigma layer  (kg, = pst*dsig*100/g)
!
      do 100 l = 1, lev
      do 100 i = 1, nxj
      pms(i,l) = (pst(i)*dsigma(l,1)+dsigma(l,2))*100.0/g
 100  continue
      do 150 i = 1, nxj
      rlsp(i) = 0.0
  150 continue
!
!     check saturated points for condensation level by level
!
      do 800 l = 1, lev
!
!     set up cloud index array 'cld'
!
      ttmp(:)=t(:,l)
      call qsatq (nxj,ttmp, pl(1,l), qs )
      do 200 i = 1, nxj
      cld(i) = q(i,l) .gt. qs(i)
  200 continue
      nsat = 0
      do 220 i = 1, nxj
      if ( cld(i) )  nsat = nsat + 1
  220 continue
      lsppt(l) = nsat
      ipass(l) = 0
      if ( nsat .eq. 0 )  go to 550
!
!     pack saturated points into work arrays ws, wq, wt, wp
!
      inx = 0
      do 240 i = 1, nxj
      if ( cld(i) )  then
         inx = inx + 1
         ws(inx) = qs(i)
         wq(inx) = q(i,l)
         wt(inx) = t(i,l)
         wp(inx) = pl(i,l)
      endif
  240 continue
!
!     iterate for new saturation
!
      do 390 ic = 1, icmax
      do 300 i = 1, nsat
      ww = (wq(i) - ws(i)) / (1.0 + gamfac*ws(i)/(wt(i)*wt(i)))
      wq(i) = wq(i) - ww
      wt(i) = wt(i) + ww*hlcp
  300 continue
      call qsatq (nsat, wt, wp, ws )
      wmax = 0.0
      do 320 i = 1, nsat
      wmax = max( wmax, abs( wq(i)/ws(i) - 1.0 ) )
  320 continue
      if ( wmax .lt. 0.002 )  go to 400
  390 continue
      print 992, l, wmax
  992 format (1x,'***warning: lsp not converge in 5 iter, (l,wmax) = ' &
             , i3,1x,e13.6 )
  400 continue
      ipass(l) = ic
      do 410 i = 1, nxj
      ws(i) = q(i,l)
  410 continue
!
!     expand adjusted q and t to full grids
!
      inx = 0
      do 500 i = 1, nxj
      if ( cld(i) )  then
         inx = inx + 1
         q(i,l) = wq(inx)
         t(i,l) = wt(inx)
      endif
  500 continue
!
!     update rain fall amounts generated in this call only
!
      do 520 i = 1, nxj
      wp(i) =  ws(i) - q(i,l)
      rlsp(i) = rlsp(i) + pms(i,l)*wp(i)
  520 continue
!
  550 continue
!
      nrain = 0
      do 560 i = 1, nxj
      if ( rlsp(i) .gt. 0.0 )  nrain = nrain + 1
  560 continue
!
!     evaporation into the next lower layer
!
      if ((nrain.eq.0).or.(evaprh.le.0.0).or.(l.ge.(lev-1))) go to 800
      lp = l + 1
      ttmp(:)=t(:,lp)
      call qsatq (nxj,ttmp, pl(1,lp), ws )
      do 600 i = 1, nxj
      wq(i) = pms(i,lp) * q(i,lp)
      wt(i) = wq(i) + rlsp(i)
      wp(i) = evaprh*ws(i)*pms(i,lp)
  600 continue
!
      do 700 i = 1, nxj
      if ( wq(i) .lt. wp(i) )  then
         if ( wt(i) .le. wp(i) )  then
!
!     case (1) all falling rain evapolates in the next lower layer
!
            ws(i) = rlsp(i)/pms(i,lp)
            q(i,lp) = q(i,lp) + ws(i)
            t(i,lp) = t(i,lp) - ws(i)*hlcp
            rlsp(i) = 0.0
         else
!
!     case (2) part of falling rain evaporates in the next lower layer
!     so that the relative humidity is evaprh, rest rain falls through
!
            wev = wp(i) - wq(i)
            rlsp(i) = rlsp(i) - wev
            ws(i) = wev / pms(i,lp)
            q(i,lp) = q(i,lp) + ws(i)
            t(i,lp) = t(i,lp) - ws(i)*hlcp
         endif
      endif
  700 continue
!
  800 continue
!
!
      return
      end
