      subroutine reducegrid(qtt,jtrun,qttcut,jf,lcapd,lonfd,octah)
!
! subroutine greduceg		programmer: hann-ming juang
!                                  modify : jen_her Chen
!
! purpose: global reduce grid initial routine to compute lcapd and lonfd
!          the maximal value of qtt and the accuray of the digit (ndigit)
!          are used to determine the resolution lcapd and grid point
!          lonfd for each latitude. The lonfd is determined by the
!          checks of factors only 2, 3 and 5 with at least one 2, and
!          it is larger than lcapd*3+1.
!
! input:
!	qtt	coefficient for spectral transform at given lat
!	lnt2	dimension of qtt at given lat
!	jtrun	wave dimension
!	qttcut	minimax value of accuracy for qtt
!       jf      position of y-direction at given lat
!       octah   option of Octahedral Gaussian reduced grid
! output:
!	lcapd	wave resolution for reduced grid
!	lonfd	number of reduced grid point for latitude
!
      use const, only: RTYPE
!
      implicit  none

      integer, parameter :: nibm=40
!
      integer  jtrun,j,k,l,m,n,lonfd,lcapd,ind,mwave,jf
      integer  need,lonfi,lonff,lonfo,ii,lonf,jtime,ktime,ltime,        &
               mtime,ntime

      real(kind=RTYPE) qtt(jtrun+1,jtrun+1)
      real     qttcut
!
!!!      integer  ibmfft(nibm)
      logical octah
!!      data ibm/.true./
!!      data ibm/.false./
!!      data ibmfft/2560,2304,2112,2048,1920,1680,1536,1440,1320,1280     &
!!      ,           1152,1056, 960, 840, 768, 640, 576, 512, 480, 384     &
!!      ,            360, 330, 288, 256, 240, 192, 180, 168, 144, 132     &
!!      ,            120,  96,  72,  64,  48,  36,  24,  16,  12,   8/
!
!     Calculate the minimum grid length at given latitude for resolving spectrum
!
      ind=0
      mwave=0
!
      if( octah ) then
       do m=1,jtrun
         need=0
         do n=m,jtrun
           if(abs(qtt(n,m)).ge.qttcut) then
              need=1
           endif
         enddo
         mwave=mwave+need
       enddo
!byl       lcapd=mwave
!
!byl       lonfi=2*(lcapd-1)+1
!byl       lonff=lonfi+mod(lonfi,2)
!     Octahedral Gaussian reduced grid with Cubic truncation
       lonfo=16+4*jf
       lcapd=min((lonfo-1)/2+1,jtrun)
!byl       if ( lonfo .lt. lonff )                                          &
!byl           print *,'Warning!! not enough grid lengh at j=',jf
!byl         go to 200
      else
       do m=1,jtrun
         need=0
         do n=m,jtrun
           if(abs(qtt(n,m)).ge.qttcut) then
              need=1
           endif
         enddo
         mwave=mwave+need
       enddo
       lcapd=mwave
!
       lonfi=2*(lcapd-1)+1
       lonff=lonfi+mod(lonfi,2)
#ifdef USE_FFTW
        lonfo=lonff
#else
       do 100 ii=0,100,2
        lonf=lonff+ii
        lonfo=lonf
        jtime=nint(log(float(lonf))/log(2.))
        do j=1,jtime
          if( mod(lonf,2).eq.0 ) then
            lonf=lonf/2
            if( lonf.eq.1 ) go to 200
          endif
        enddo
        ktime=nint(log(float(lonf))/log(3.))
        do k=1,ktime
          if( mod(lonf,3).eq.0 ) then
            lonf=lonf/3
            if( lonf.eq.1 ) go to 200
          endif
        enddo
        ltime=nint(log(float(lonf))/log(5.))
        do l=1,ltime
          if( mod(lonf,5).eq.0 ) then
            lonf=lonf/5
            if( lonf.eq.1 ) go to 200
          endif
        enddo
        mtime=nint(log(float(lonf))/log(7.))
        do l=1,mtime
          if( mod(lonf,7).eq.0 ) then
            lonf=lonf/7
            if( lonf.eq.1 ) go to 200
          endif
        enddo
        ntime=nint(log(float(lonf))/log(11.))
        do l=1,ntime
          if( mod(lonf,11).eq.0 ) then
            lonf=lonf/11
            if( lonf.eq.1 ) go to 200
          endif
        enddo
 100   continue
!
#endif
      endif
!
 200  lonfd=lonfo
!
      return
      end
