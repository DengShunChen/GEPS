      subroutine qmax2d (fld,i1,j1,im,jm)
!
!  subroutine to print max and min in 2-d layer (or a subarray)
!  within a 3-d field stored with n-s index slowest varying
!
!  ***input***
!
!  fld: input array
!  i1: starting index of first dimension (e-w)
!  j1: starting index of third dimension (n-s)
!  im: first dimension
!  jm: third dimension
!
! *****************************************************************
!
      use rank
      use index

      real fld(im,jm)
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
      if (fld(i,j).lt.xmin) then
      xmin= fld(i,j)
      jmin= j
      imin= i
      endif
      if (fld(i,j).gt.xmax) then
      xmax= fld(i,j)
      jmax= j
      imax= i
      endif
   10 continue
!
      if(myrank .eq. 0) print 8995, imax,jmax,xmax,imin,jmin,xmin
!
 8995 format(' imax=',i4,' jmax=',i4,' xmax=',g20.12                   &
      ,' imin=',i4,' jmin=',i4,' xmin=',g20.12)
!
      return
      end
