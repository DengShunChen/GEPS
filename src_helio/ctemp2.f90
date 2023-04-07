      subroutine ctemp2 ( mn,lm,g,r,pr,zz,tp,zc,mode )
!
!     convert geopotential height to temperature on lm*2-1 p-levels
!
!  input parameters:
!
!     mn       :  horizontal dimension
!     lm       :  vertical dimension
!     g        :  gravity (=9.806)
!     r        :  gas constant (=287.04)
!     pr(lm)   :  presure (mb)
!     zz(mn,lm):  geopotential height of the pr pressure levels (m)
!     zc       :  converge criterion of hydrostatic iterations (m)
!     mode     :  mode to compute temperature at the bottom (lm) level;
!                 = 0     , using std atmosphere lapse rate
!                 = others, extrapolating from (lm-1) and (lm-1/2)
!
!  output parameter:
!
!     tp(mn,lm*2-1): temperature at the pressure levels and levels in
!                    between them  ( tp(lm*2-1) is for pr(lm) level )
!
!  modify to f90 by C-H Lee and sort by River Chen in 2015
!
      use rank

      implicit  none

      integer   mn,lm,mode
      real      pr(lm),zz(mn,lm),tp(mn,lm*2-1)
!
!  work arrays
!
      real      zdif(mn),znew(mn,lm),tt(mn,lm),plog(lm*2-1),zhyd(lm*2-1)

      integer   lmd2,lmp1,lm2m1,lm2m2,lm2m3,lmtd2,kt,kb,i,l,l1,l2,l3,iter
      real      zc,r,g,pxxx,dlp,tpx,xl1,zmax
!
      lmd2 = lm / 2
      lmp1 = lm + 1
      lm2m1= lm*2 - 1
      lm2m2= lm*2 - 2
      lm2m3= lm*2 - 3
      lmtd2= lm2m1 / 2
!
!***********************************************************************
!          flip vertical index of pr and zz for ctemp simulation
!***********************************************************************
!
      do 100 kt = 1, lmd2
      kb = lmp1 - kt
      pxxx   = pr(kt)
      pr(kt) = pr(kb)
      pr(kb) = pxxx
      do 110 i = 1,mn
      zdif(i)  = zz(i,kt)
      zz(i,kt) = zz(i,kb)
      zz(i,kb) = zdif(i)
  110 continue
  100 continue
!
!***********************************************************************
!          compute hydrostatic temperature at middle p-levels
!                 ( use tt as the work array to store it )
!***********************************************************************
!
      do 200 l = 1, lm2m1, 2
      plog(l) = log(pr((l+1)/2))
  200 continue
      do 210 l = 2, lm2m2, 2
      plog(l) = 0.5*(plog(l-1)+plog(l+1))
  210 continue
      do 220 l = 1, lm2m2
      dlp = plog(l)-plog(l+1)
      zhyd(l) = dlp*r/g
  220 continue
      do 240 l = 2, lm
      l1 = 2*l- 3
      l2 = l1 + 1
      do 230 i = 1, mn
      tt(i,l-1) = (zz(i,l)-zz(i,l-1))/(zhyd(l1)+zhyd(l2))
  230 continue
  240 continue
!
      tpx = (pr(1)/exp(plog(2)))**0.19023
!
      do 390 iter = 1 , 10
!
      do 310 l = 2, lm2m2, 2
      do 310 i = 1, mn
      tp(i,l)  = tt(i,l/2)
  310 continue
      do 325 l = 3, lm2m2, 2
      l1 = (l-1)/2
      xl1 = (plog(l-1)-plog(l))/(plog(l-1)-plog(l+1))
      do 320 i = 1,mn
      tp(i,l) = tt(i,l1)+(tt(i,l1+1)-tt(i,l1))*xl1
  320 continue
  325 continue
      if ( mode .eq. 0 )  then
         do 330 i = 1, mn
         tp(i,1) = tt(i,1) * tpx
  330    continue
      else
         do 335 i = 1, mn
         tp(i,1) = tp(i,2)+(tp(i,2)-tp(i,3))*(plog(1)-plog(2))       &
                 / (plog(2)-plog(3))
  335    continue
      endif
      do 340 i = 1, mn
      tp(i,lm2m1) = tp(i,lm2m2) + 1.5e-3*( zz(i,lm) - (zz(i,lm-1)    &
                  + 0.5*(tp(i,lm2m2)+tp(i,lm2m3))*zhyd(lm2m3)) )
  340 continue
!
      do 350 i = 1, mn
      znew(i,1) = zz(i,1)
  350 continue
      zmax = 0.0
      do 370 l = 2, lm
      l1 = 2*l- 3
      l2 = l1 + 1
      l3 = l2 + 1
      do 360 i = 1, mn
      znew(i,l)  = znew(i,l-1) + 0.5*(tp(i,l1)+tp(i,l2))*zhyd(l1)    &
                               + 0.5*(tp(i,l2)+tp(i,l3))*zhyd(l2)
      zdif(i)    = znew(i,l) - znew(i,l-1) - zz(i,l) + zz(i,l-1)
      tt(i,l-1)  = tt(i,l-1) - zdif(i)/(zhyd(l1)+zhyd(l2))
      zdif(i)    = abs(zdif(i))
  360 continue
      do 365 i = 1, mn
      zmax = max( zmax, zdif(i) )
  365 continue
  370 continue
      if ( zmax .lt. zc ) go to 400
  390 continue
      if(myrank .eq. 0) print 800, zc, zmax
  800 format(1x,'** warning: ctemp2 does not converage in 10'        &
            ,' iterations !!  zcrit =',f6.2,', zdifmax = ',e13.6)
  400 continue
!
!***********************************************************************
!          flip vertical index of pr and zz back to normal index
!***********************************************************************
!
      do 420 kt = 1, lmd2
      kb = lmp1 - kt
      pxxx   = pr(kt)
      pr(kt) = pr(kb)
      pr(kb) = pxxx
      do 410 i = 1, mn
      zdif(i)  = zz(i,kt)
      zz(i,kt) = zz(i,kb)
      zz(i,kb) = zdif(i)
  410 continue
  420 continue
!
!***********************************************************************
!          flip vertical index of temperatures
!***********************************************************************
!
      do 440 kt = 1, lmtd2
      kb = lm2m1+1 - kt
      do 430 i = 1, mn
      zdif(i)  = tp(i,kt)
      tp(i,kt) = tp(i,kb)
      tp(i,kb) = zdif(i)
  430 continue
  440 continue
!
      return
      end

      subroutine ctemp3 ( mn,lm,g,r,pr,zz,tp,zc,mode )
!
!     convert geopotential height to temperature on lm*2-1 p-levels
!
!  input parameters:
!
!     mn       :  horizontal dimension
!     lm       :  vertical dimension
!     g        :  gravity (=9.806)
!     r        :  gas constant (=287.04)
!     pr(mn,lm)   :  presure (mb)
!     zz(mn,lm):  geopotential height of the pr pressure levels (m)
!     zc       :  converge criterion of hydrostatic iterations (m)
!     mode     :  mode to compute temperature at the bottom (lm) level;
!                 = 0     , using std atmosphere lapse rate
!                 = others, extrapolating from (lm-1) and (lm-1/2)
!
!  output parameter:
!
!     tp(mn,lm*2-1): temperature at the pressure levels and levels in
!                    between them  ( tp(lm*2-1) is for pr(lm) level )
!
!  modify to f90 by C-H Lee and sort by River Chen in 2015
!
      use rank

      dimension pr(mn,lm),zz(mn,lm),tp(mn,lm*2-1)
!
!  work arrays
!
      dimension zdif(mn),pdif(mn),znew(mn,lm),tt(mn,lm)            &
              , plog(mn,lm*2-1),zhyd(mn,lm*2-1)
!
      lmd2 = lm / 2
      lmp1 = lm + 1
      lm2m1= lm*2 - 1
      lm2m2= lm*2 - 2
      lm2m3= lm*2 - 3
      lmtd2= lm2m1 / 2
!
!***********************************************************************
!          flip vertical index of pr and zz for ctemp simulation
!***********************************************************************
!
      do 100 kt = 1, lmd2
      kb = lmp1 - kt
      do 110 i = 1,mn
      zdif(i)  = zz(i,kt)
      zz(i,kt) = zz(i,kb)
      zz(i,kb) = zdif(i)
      pdif(i)  = pr(i,kt)
      pr(i,kt) = pr(i,kb)
      pr(i,kb) = pdif(i)
  110 continue
  100 continue
!
!***********************************************************************
!          compute hydrostatic temperature at middle p-levels
!                 ( use tt as the work array to store it )
!***********************************************************************
!
      do 200 l = 1, lm2m1, 2
      do 200 i = 1, mn
      plog(i,l) = log(pr(i,(l+1)/2))
  200 continue
      do 210 l = 2, lm2m2, 2
      do 210 i = 1, mn
      plog(i,l) = 0.5*(plog(i,l-1)+plog(i,l+1))
  210 continue
      do 220 l = 1, lm2m2
      do 220 i = 1, mn
      dlp = plog(i,l)-plog(i,l+1)
      zhyd(i,l) = dlp*r/g
  220 continue
      do 240 l = 2, lm
      l1 = 2*l- 3
      l2 = l1 + 1
      do 230 i = 1, mn
      tt(i,l-1) = (zz(i,l)-zz(i,l-1))/(zhyd(i,l1)+zhyd(i,l2))
  230 continue
  240 continue
!
      do i = 1, mn
!cc   pdif(i) = (pr(i,1)/exp(plog(i,2)))**0.19023
      fff = 0.19023 * ( log(pr(i,1)) - plog(i,2) )
      pdif(i) = exp( fff )
      end do
!
      do 390 iter = 1 , 10
!
      do 310 l = 2, lm2m2, 2
      do 310 i = 1, mn
      tp(i,l)  = tt(i,l/2)
  310 continue
      do 325 l = 3, lm2m2, 2
      l1 = (l-1)/2
      do 320 i = 1,mn
      xl1 = (plog(i,l-1)-plog(i,l))/(plog(i,l-1)-plog(i,l+1))
      tp(i,l) = tt(i,l1)+(tt(i,l1+1)-tt(i,l1))*xl1
  320 continue
  325 continue
      if ( mode .eq. 0 )  then
         do 330 i = 1, mn
         tp(i,1) = tt(i,1) * pdif(i)
  330    continue
      else
         do 335 i = 1, mn
         tp(i,1) = tp(i,2)+(tp(i,2)-tp(i,3))*(plog(i,1)-plog(i,2))    &
                 / (plog(i,2)-plog(i,3))
  335    continue
      endif
      do 340 i = 1, mn
      tp(i,lm2m1) = tp(i,lm2m2) + 1.5e-3*( zz(i,lm) - (zz(i,lm-1)     &
                  + 0.5*(tp(i,lm2m2)+tp(i,lm2m3))*zhyd(i,lm2m3)) )
  340 continue
!
      do 350 i = 1, mn
      znew(i,1) = zz(i,1)
  350 continue
      zmax = 0.0
      do 370 l = 2, lm
      l1 = 2*l- 3
      l2 = l1 + 1
      l3 = l2 + 1
      do 360 i = 1, mn
      znew(i,l)  = znew(i,l-1) + 0.5*(tp(i,l1)+tp(i,l2))*zhyd(i,l1)   &
                               + 0.5*(tp(i,l2)+tp(i,l3))*zhyd(i,l2)
      zdif(i)    = znew(i,l) - znew(i,l-1) - zz(i,l) + zz(i,l-1)
      tt(i,l-1)  = tt(i,l-1) - zdif(i)/(zhyd(i,l1)+zhyd(i,l2))
      zdif(i)    = abs(zdif(i))
  360 continue
      do 365 i = 1, mn
      zmax = max( zmax, zdif(i) )
  365 continue
  370 continue
      if ( zmax .lt. zc ) go to 400
  390 continue
      if(myrank .eq. 0) print 800, zc, zmax
  800 format(1x,'** warning: ctemp2 does not converage in 10'         &
            ,' iterations !!  zcrit =',f6.2,', zdifmax = ',e13.6)
  400 continue
!
!***********************************************************************
!          flip vertical index of pr and zz back to normal index
!***********************************************************************
!
      do 420 kt = 1, lmd2
      kb = lmp1 - kt
      do 410 i = 1, mn
      zdif(i)  = zz(i,kt)
      zz(i,kt) = zz(i,kb)
      zz(i,kb) = zdif(i)
      pdif(i)  = pr(i,kt)
      pr(i,kt) = pr(i,kb)
      pr(i,kb) = pdif(i)
  410 continue
  420 continue
!
!***********************************************************************
!          flip vertical index of temperatures
!***********************************************************************
!
      do 440 kt = 1, lmtd2
      kb = lm2m1+1 - kt
      do 430 i = 1, mn
      zdif(i)  = tp(i,kt)
      tp(i,kt) = tp(i,kb)
      tp(i,kb) = zdif(i)
  430 continue
  440 continue
!
      return
      end
