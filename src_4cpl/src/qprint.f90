      subroutine qprint ( a,t1,t2,ic,jc,m,n )
!
!  qprint prints an array a(m,n) starting at address a(ic,jc).
!  values are automatically scaled to allow integer format printing.
!
!     a     = input of m x n array
!     t1,t2 = title 1 and 2  (1 word ascii, or 8 characters for each)

!     ic,jc = lower left corner coords to be printed
!             ( up to 43 x n points printed )
!
      implicit  none
      integer   ic,jc,m,n
      real      a(m,n)
      character*8 t1, t2

      integer   ix(43),ie,jl,i,j,kp,kpm,ii,jli
      real      xm,af,fk
!
!     determine grid limits
!
      ie = min( ic+42, m )
      jl = n
!
!     index backwards checking for max
!
      xm = 0.0
      do 100 j = jc, jl
      do 100 i = ic, ie
      af = a(i,j)
      xm = max( xm, abs(af) )
  100 continue
!
!     determine scaling factor limits
!
      if ( xm .lt. 1.e-200 ) xm=99.0
!     xm = alog10( 99./xm )
      xm = dlog10( 99./xm ) ! CWB 2015
      kp = xm
      if( xm .lt. 0.0 ) kp = kp-1
!
!     print scaling constants
!
      kpm = - kp
      print 900, t1, t2, kpm, (i,i=ic,ie,2)
  900 format(1h0,a8,a8,'   k=',i3,2x,'( true=print*10**k )', /1x,22i6)
      fk=10.**kp
!
!     quick print field
!
      ii = ie - ic + 1
      do 200 j = jc,jl
      jli = jl - j + jc
         do 400 i = ic, ie
         af = a(i,jli)
         ix(i-ic+1) = af*fk + sign(.5,af)
  400    continue
      print 950, jli, (ix(i),i=1,ii), jli
  950 format(i4,44i3)
  200 continue
      return
      end
