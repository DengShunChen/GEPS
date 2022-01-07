      subroutine qmaxn3 (fld,t1,t2,i1,j1,k1,im,jm,lm)
!
!  subroutine to print max and min in 2-d layer (or a subarray)
!  within a 3-d field stored with n-s index slowest varying
!
!  ***input***
!
!  fld: input array
!  t1,t2: 2* 8 character caption
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

      implicit   none
      integer    i1,j1,k1,im,jm,lm
      real      fld(im,lm,jm)
!ch   character*14 t1, t2
      character t1*16, t2*12

      real       xmin,xmax
      integer    imin,jmin,imax,jmax,j,i
!
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
      if(myrank .eq. 0) print 9000, t1, t2
      if(myrank .eq. 0) print 8995, imax,jmax,xmax,imin,jmin,xmin
!
!ch 9000 format (2a14)
 9000 format (a16,a12)
 8995 format(' imax=',i4,' jmax=',i4,' xlarg=',g20.12  &
      ,' imin=',i4,' jmin=',i4,' xsmal=',g20.12)
!
      return
      end
