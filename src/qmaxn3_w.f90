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
      real(kind=RTYPE) tmp
!ch   character*14 t1, t2

      real       xmin,xmax
      integer    imin,jmin,imax,jmax,j,i
      integer,parameter:: async_id =1
!
#ifdef O38K
      ihdgleno1=ihdgo(1:16)
      ihdgleno2=ihdgo(17:28)
#else
      ihdgleno1=ihdgo(1:14)
      ihdgleno2=ihdgo(15:26)
#endif
      xmin=  1.0e25
      xmax= -1.0e25
      imin=-1
      jmin=-1
      imax=-1
      jmax=-1
!

#ifdef USE_CUDA
!$acc data copyin(xmin,xmax,imin,jmin,imax,jmax) async(async_id)

!$acc parallel loop collapse(2) present(xmin,xmax)  &
!$acc&  reduction(min:xmin) reduction(max:xmax)  async(async_id) 
      do j = j1, jm
        do i = i1, im
          xmin = min(xmin, fld(i,k1,j) )
          xmax = max(xmax, fld(i,k1,j) )
        end do
      end do

!!$acc parallel loop collapse(2) &
!!$acc& present(xmin,xmax,imin,jmin,imax,jmax) async(async_id)
!      do j = j1, jm
!        do i = i1, im
!          if( xmin == fld(i,k1,j) ) then
!            jmin = j
!            imin=  i
!          endif
!          if( xmax == fld(i,k1,j) ) then
!            jmax= j
!            imax= i
!          endif
!        end do
!      end do

!$acc wait(async_id)
!$acc update self(xmin,xmax,imin,jmin,imax,jmax) async(async_id)
!$acc wait(async_id)
#else

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

#endif
!
!

      if(myrank .eq. 0) print 9000, ihdgleno1,ihdgleno2
      if(myrank .eq. 0) print 8995, imax,jmax,xmax,imin,jmin,xmin
!
!ch 9000 format (2a14)
 9000 format (a16,a12)
 8995 format(' imax=',i4,' jmax=',i4,' xlarg=',g20.12  &
      ,' imin=',i4,' jmin=',i4,' xsmal=',g20.12)

!
!$acc end data
      return
      end
