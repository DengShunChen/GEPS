      subroutine trdiv2 (m,lm,a,b,c,x)
!
!  tri-diagonal gaussian elimination subroutine called by varpht
!
      implicit  none
      integer   m,lm,i,j,k
      real      work
      real      a(m,lm)   ,b(m,lm)   ,c(m,lm),  x(m,lm)
!
      do 90 i=1,m
      work = 1.0/b(i,1)
      c(i,1)= c(i,1)*work
      x(i,1) = x(i,1)*work
   90 continue
!
! gaussian elimination
!
      do 101 j=2,lm-1
!ocl novrec
      do 103 i=1,m
	  work = 1.0/(b(i,j)-a(i,j)*c(i,j-1))
	  c(i,j) = c(i,j)*work
      x(i,j)= (x(i,j)-a(i,j)*x(i,j-1))*work
  103 continue
  101 continue
!
      do 105 i=1,m
      x(i,lm) = (x(i,lm)-a(i,lm)*x(i,lm-1))/(b(i,lm)-a(i,lm)*c(i,lm-1))
  105 continue
!
! backwards substitution
!
      do 104 k= lm-1,1,-1
!ocl novrec
      do 104 i=1,m
	  x(i,k) = x(i,k)-c(i,k)*x(i,k+1)
  104 continue
      return
      end
