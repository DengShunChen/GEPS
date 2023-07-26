      subroutine trdih (m,n,y)
!
!  vectorized tri-diagonal gaussian elimintaion solver
!
      implicit  none
      integer   m,n,i,j,k
      real      y(m,n),c(2000)
!
! gaussian elimination
!
      c(1) = 0.25
!
      do 201 i=1,m
      y(i,1)= y(i,1)*c(1)
  201 continue
!
      do 103 j=2,n-1
      c(j)= 1.0/(4.0-c(j-1))
!dir$ ivdep
      do 103 i=1,m
      y(i,j)= (y(i,j)-y(i,j-1))*c(j)
  103 continue
!
!dir$ ivdep
      do 202 i=1,m
      y(i,n)=(y(i,n)-y(i,n-1))/(4.0-c(n-1))
  202 continue
!
! backwards substitution
!
      do 104 k=n-1,1,-1
!dir$ ivdep
      do 104 i=1,m
      y(i,k)= y(i,k)-c(k)*y(i,k+1)
  104 continue
      return
      end
