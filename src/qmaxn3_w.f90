      subroutine qmaxn3_w (fld,i1,j1,k1,im,jm,lm)
!
!  subroutine to print max and min in 2-d layer (or a subarray)
!  within a 3-d field stored with n-s index slowest varying
!
!  ***input***
!
!  fld: input array
!  ihdglen1,ihdglen2: 2* 8 character caption
!  i1: starting index of first dimension (e-w)
!  j1: starting index of third dimension (n-s)
!  k1: starting index of second dimension (vertical)
!  im: first dimension
!  jm: third dimension
!  lm: second dimension
!
! *****************************************************************
!
      use rank
      use index
      use const, only: RTYPE,ihdgo,ihdgleno1,ihdgleno2

      implicit   none
      integer    i1,j1,k1,im,jm,lm
      real(kind=RTYPE) fld(im,lm,jm)
!ch   character*14 t1, t2

      real       xmin,xmax
      integer    imin,jmin,imax,jmax,j,i
!
#ifdef O38K
      ihdgleno1=ihdgo(1:16)
      ihdgleno2=ihdgo(17:28)
#else
      ihdgleno1=ihdgo(1:14)
      ihdgleno2=ihdgo(15:26)
#endif
      xmin= 1.0e25
      xmax= -1.0e25
      imin=1
      jmin=1
      imax=1
      jmax=1
!
      do 10 j=j1,jm
      do 10 i=i1,im
      if (fld(i,k1,j).le.xmin) then
      xmin= fld(i,k1,j)
      jmin= j
      imin= i
      endif
      if (fld(i,k1,j).gt.xmax) then
      xmax= fld(i,k1,j)
      jmax= j
      imax= i
      endif
   10 continue
!
      if(myrank .eq. 0) print 9000, ihdgleno1,ihdgleno2
      if(myrank .eq. 0) print 8995, imax,jmax,xmax,imin,jmin,xmin
!
!ch 9000 format (2a14)
 9000 format (a16,a12)
 8995 format(' imax=',i4,' jmax=',i4,' xlarg=',g20.12  &
      ,' imin=',i4,' jmin=',i4,' xsmal=',g20.12)
!
      return
      end
