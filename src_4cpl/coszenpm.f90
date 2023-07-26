      subroutine coszenpm ( nx,my,my_max,julian,solhr,sinl,cosl,xlonr,prad, &
                            slag,sdec,cdec,cosz)
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
      integer   nx,my,julian,jj,i,j,nxj,it,my_max,ii
      real      sinl(my),cosl(my),xlonr(nxp,my_max),cosz(nxp,my_max)
      real      solhr,pi,d2r,prad,pid
!      real      timex,prad,dlon,sinz0,cosz0
!      real      timeg,plon,hourg
      real      slag,sdec,cdec,ss,cc,cns,solang,rstp
!
      real      coszen(nxp),coszdg(nxp),rst,coszn
      integer   istsun(nxp),nst

!
      pi = 4.*atan(1.0)
!      d2r = pi/180.0
!      beta = 2.*pi*julian/365.
      nst = 6
      pid = pi / 12.
      rst = prad / float(nst-1)

      solang = solhr - 12.         ! solar angle at present time
!      rstp = 1.0 / float(nst)



      do jj = 1, jlistnum
        j   = jlist1(jj)
        nxj = nxdef_2d(j)
        ii  = nxjstart(j)
        ss  = sinl(j) * sdec
        cc  = cosl(j) * cdec
        do i = 1, nxj
          cosz(i,jj)= 0.0
          coszen(i) = 0.0
          istsun(i) = 0
        enddo

        do it = 1, nst
          cns = pid * ( solang + float(it-1)*rst ) + slag
          do i = 1, nxj
            coszn = ss + cc * cos(cns + xlonr(i,jj))
            coszen(i) = coszen(i) + max(0.0, coszn)
            if (coszn > 0.0001) istsun(i) = istsun(i) + 1
          enddo
        enddo

!  --- ...  compute time averages

        do i = 1, nxj
!          coszdg(i) = coszen(i) * rstp
          if (istsun(i) > 0) cosz(i,jj) = coszen(i) / istsun(i)
        enddo
      enddo

!      xlat = 0.006918 - 0.399912*cos(beta) + 0.070257*sin(beta)         &
!           - 0.006758*cos(2.*beta) + 0.000907*sin(2.*beta)              &
!           - 0.002697*cos(3.*beta) + 0.001480*sin(3.*beta)
!      sinxl = sin(xlat)
!      cosxl = cos(xlat)

!      hdif = 0.000075 + 0.001868*cos(beta) - 0.032077*sin(beta)         &
!           - 0.014615*cos(2.*beta) - 0.040849*sin(2.*beta)
!      sinxl = sdec
!      cosxl = cdec
!      hdif  = slag
!
!      do 200 jj = 1, jlistnum

!       j=jlist1(jj)
!       nxj=nxdef_2d(j)
!       ii=nxjstart(j)
!      dlon = 360./float(nxj)
!       sinz0  = sinxl*sinl(j)
!       cosz0  = cosxl*cosl(j)

!       do i=1,nxj
!          istsun(i)=0
!          coszen(i)=0.0
!       enddo
!
!       do 100 it=1,nst
!          timex = time+rst*float(it-1)
!          timex = mod ( timex, 24. )
!          dtzu = timex - 12.0
!          if(dtzu .lt. 0.) dtzu=dtzu+24.
!
!        do 100 i = 1, nxj
!           plon =  alon(ii,jj)
!           if (plon .lt. 0.)plon=plon+360.
!           timeg = (dtzu*15.0 + plon) * d2r
!           hourg = timeg + hdif
!           coszn = cosz0*cos(hourg) + sinz0
!           if (coszn > 0.0001) istsun(i)=istsun(i)+1
!           coszen(i) = coszen(i)+ max(0.0,coszn)
!           ii=ii+1
!  100   continue

!        do 150 i = 1, nxj
!           coszdg(i)=coszen(i)/float(nst)
!           cosz(i,jj)=coszen(i)
!           if (istsun(i) > 0 ) cosz(i,jj)=coszen(i)/float(istsun(i))
!  150   continue
!  200 continue
!
      return
      end
