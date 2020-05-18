      subroutine ccoszen ( nx,my,my_max,julian,time,alat,alon,cosz)
!
!     compute cos(zenith angle) from day, zulu-time, alat and alon
!     for global model only ( alat(my),  alon(nx,my) )
!
!  modify to f90 in 2015 by C-H Lee and sort by River Chen in 2015
!
      use index

      implicit  none

      integer   nx,my,my_max,julian
      real      time
      real      alat(my),alon(nx,my_max),cosz(nxp,my_max)
!
      integer   i,jj,j,nxj,ii
      real      pi,d2r,beta,xlat,sinxl,cosxl,hdif,dtzu
      real      timeg,hourg,plon,sinz0,cosz0,dlon
!
      pi = 4.0*atan(1.0)
      d2r = pi/180.0
      beta = 2.*pi*julian/365.
      xlat = 0.006918 - 0.399912*cos(beta) + 0.070257*sin(beta) &
           - 0.006758*cos(2.*beta) + 0.000907*sin(2.*beta)      &
           - 0.002697*cos(3.*beta) + 0.001480*sin(3.*beta)
      sinxl = sin(xlat)
      cosxl = cos(xlat)
!
      hdif = 0.000075 + 0.001868*cos(beta) - 0.032077*sin(beta) &
           - 0.014615*cos(2.*beta) - 0.040849*sin(2.*beta)
      dtzu = time - 12.0
      if(dtzu .lt. 0.) dtzu=dtzu+24.
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
         timeg = (dtzu*15.0 + plon) * d2r
         hourg = timeg + hdif
         cosz(i,jj) = cosz0*cos(hourg) + sinz0
  100    continue
  200 continue
!
      return
      end
