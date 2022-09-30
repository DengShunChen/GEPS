      subroutine mtxprt (a,m,n,lab,fmt)
!
!  matrix print routine
!
!  ***input***
!
!  a: matrix to be printed
!  m: first dimension
!  n: second dimension
!  lab: caption to be printed (8 character)
!  fmt: format used for printing
!
! *************************************************************
!

      use rank
      use const, only : RTYPE

      implicit  none

      integer   m,n,nx,kk,is,ie,i,j
      real(kind=RTYPE) a(m,n)
      character*8 lab,fmt
      character*16 cfmt

      nx= (n-1)/15+1
      do 5 kk=1,nx
      is= 1+(kk-1)*15
      ie= min(is+14,n)
!
      if(myrank .eq. 0) print 100, lab
!
  100 format(/'0matrix ',a8)
      write(cfmt,200) fmt
  200 format('(1x,',a8,')')
!
      if(myrank .eq. 0) then
      do 1 i=1,m
      print cfmt,(a(i,j),j=is,ie)
    1 continue
      endif
!
    5 continue
      return
      end
