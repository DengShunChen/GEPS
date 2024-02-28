      subroutine trdivv (m,n,a,y)
!
!  tri-diagonal gaussian elimination subroutine call by bicubv
!
      use const, only: RTYPE
!
      implicit   none
      integer    m,n,nm,i,j,k
      real(kind=RTYPE) a(m,n), y(m,n)
      real(kind=RTYPE) c(m,n-1)
!
      nm = n-1
      do 201 i=1,m
      c(i,1)= 0.5/(1.0+a(i,1))
      y(i,1) = y(i,1)*c(i,1)
  201 continue
!
! gaussian elimination
!
      do 101 j=2,nm
!ocl novrec
      do 101 i=1,m
      c(i,j)= 1.0/(2.0+a(i,j)*(2.0-c(i,j-1)))
      y(i,j) = (y(i,j)-a(i,j)*y(i,j-1))*c(i,j)
  101 continue
!
      do 202 i=1,m
      y(i,n)= (y(i,n)-a(i,n)*y(i,nm))/(2.0+a(i,n)*(2.0-c(i,nm)))
  202 continue
!
! backwards substitution
!
      do 104 k=nm,1,-1
!ocl novrec
      do 104 i=1,m
      y(i,k)= y(i,k)-c(i,k)*y(i,k+1)
  104 continue
      return
      end
