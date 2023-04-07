      subroutine qmaxn3p(fld,t1,t2,i1,j1,k1,im,jm,lm)
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
      use mpe
      use rank
      use index
!
      implicit   none
      integer    i1,j1,k1,im,jm,lm
      real       fld(im,lm,jm)
!
      character*8 t1, t2

      real       xmin,xmax
      integer    imin,jmin,imax,jmax,jj,j,i,idummy
!
      xmin= 1.0e25
      xmax= -1.0e25
      imin=1
      jmin=1
      imax=1
      jmax=1
!
      do 20 jj =1,jlistnum
      j=jlist1(jj)
      if(j .ge. j1) then
      do 10 i=i1,im
      if (fld(i,k1,jj).le.xmin) then
      xmin= fld(i,k1,jj)
      jmin= j
      imin= i
      endif
      if (fld(i,k1,jj).gt.xmax) then
      xmax= fld(i,k1,jj)
      jmax= j
      imax= i
      endif
   10 continue
      endif
   20 continue
!
      call mpe_global_maxloc(xmax,2,jmax,imax,idummy,idummy,mpe_double)
      call mpe_global_minloc(xmin,2,jmin,imin,idummy,idummy,mpe_double)
!
      if(myrank .eq. 0) print 9000, t1, t2
 9000 format (1h0,2a8)
      if(myrank .eq. 0) print 8995, imax,jmax,xmax,imin,jmin,xmin
 8995 format(' imax=',i4,' jmax=',i4,' xlarg=',g20.12   &
      ,' imin=',i4,' jmin=',i4,' xsmal=',g20.12)
!
      return
      end
