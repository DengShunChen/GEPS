      subroutine coszenpm ( nx,my,my_max,julian,time,alat,alon,prad,cosz)
!
!     compute cos(zenith angle) from day, zulu-time, alat and alon
!     for global model only ( alat(my),  alon(nx,my) )
!
!     process a mean cos within a period of time
!
!     prad : the time interval to call radiation ( unit: hr )
!
! modify to f90 in 2015 by C-H Lee and sort by River Chen in 2015
!
      use index

      implicit  none
      integer   nx,my,my_max,julian,jj,i,j,nxj,ii
      real      alat(my),alon(nx,my_max),cosz(nxp,my_max)
      real      time,pi,d2r,beta,xlat,sinxl,cosxl,hdif,dtzu1
      real      timex,prad,dtzu2,dtzu3,dtzu4,dlon,sinz0,cosz0
      real      timeg,plon,hourg,cosz1,cosz2,cosz3,cosz4
!
      pi = 3.141592654
      d2r = pi/180.0
      beta = 2.*pi*julian/365.
      xlat = 0.006918 - 0.399912*cos(beta) + 0.070257*sin(beta)         &
           - 0.006758*cos(2.*beta) + 0.000907*sin(2.*beta)              &
           - 0.002697*cos(3.*beta) + 0.001480*sin(3.*beta)
      sinxl = sin(xlat)
      cosxl = cos(xlat)
!
      hdif = 0.000075 + 0.001868*cos(beta) - 0.032077*sin(beta)         &
           - 0.014615*cos(2.*beta) - 0.040849*sin(2.*beta)
!
      dtzu1 = time - 12.0
      if(dtzu1 .lt. 0.) dtzu1=dtzu1+24.
!
      timex = time+prad/3.
      timex = mod ( timex, 24. )
      dtzu2 = timex - 12.0
      if(dtzu2 .lt. 0.) dtzu2=dtzu2+24.
!
      timex = time+prad*2./3.
      timex = mod ( timex, 24. )
      dtzu3 = timex - 12.0
      if(dtzu3 .lt. 0.) dtzu3=dtzu3+24.
!
      timex = time+prad
      timex = mod ( timex, 24. )
      dtzu4 = timex - 12.0
      if(dtzu4 .lt. 0.) dtzu4=dtzu4+24.
!
      do 200 jj = 1, jlistnum
       j=jlist1(jj)
       nxj=nxdef(j)
       dlon = 360./float(nxj)
       sinz0  = sinxl*sin(alat(j)*d2r)
       cosz0  = cosxl*cos(alat(j)*d2r)
!ch      do 100 i = 1, nxj
         do 100 i = 1, nxdef_2d(j)
!ch      plon = (i-1)*dlon
         ii=map2to1(i,j)
         plon = (ii-1)*dlon
         timeg = (dtzu1*15.0 + plon) * d2r
         hourg = timeg + hdif
         cosz1 = cosz0*cos(hourg) + sinz0
!
         timeg = (dtzu2*15.0 + plon) * d2r
         hourg = timeg + hdif
         cosz2 = cosz0*cos(hourg) + sinz0
!
         timeg = (dtzu3*15.0 + plon) * d2r
         hourg = timeg + hdif
         cosz3 = cosz0*cos(hourg) + sinz0
!
         timeg = (dtzu4*15.0 + plon) * d2r
         hourg = timeg + hdif
         cosz4 = cosz0*cos(hourg) + sinz0
!
         cosz(i,jj) = (cosz1+cosz2+cosz3+cosz4)/4.
!
  100    continue
  200 continue
!
      return
      end
