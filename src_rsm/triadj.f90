      subroutine triadj ( y,b,c,im,lm,lm2,work,stable )
!
!     tri-diagonal gaussian elimination for dryadj subroutine
!
      implicit  none

      integer   im,lm,lm2,lmm1,i,l,k,nd

      real      b(lm),c(lm),y(im,lm),work(im,lm2),big,eps

      logical   stable(im,lm)

      data eps/-1.0e-20/, big/ -1.0e+20/
!
      lmm1 = lm-1
!
!   to ensure no adjustment for stable layers, put zero forcing
!   and large diagonal values for those layers
!
      do 100 i = 1, im
      work(i,1)= 1.0 / b(1)
      if( stable(i,1) ) work(i,1) = eps
      work(i,lm) = c(1) * work(i,1)
      y(i,1) = y(i,1) * work(i,1)
  100 continue
!
!     gaussian elimination
!
      do 200 l = 2, lmm1
      nd = l+lmm1
      do 220 i = 1, im
      work(i,l) = b(l)
      if ( stable(i,l) ) work(i,l) = big
      work(i,l) = 1.0 / (work(i,l) - work(i,nd-1))
      work(i,nd) = c(l) * work(i,l)
  220 continue
  200 continue
      do 300 l = 2, lmm1
      do 300 i = 1, im
      y(i,l) = (y(i,l) - y(i,l-1))*work(i,l)
  300 continue
      nd = lmm1+lmm1
      do 400 i = 1, im
      y(i,lm) = y(i,lm) - y(i,lmm1)
      work(i,1) = b(lm)
      if( stable(i,lm) )  work(i,1) = big
      y(i,lm) = y(i,lm) / (work(i,1) - work(i,nd))
  400 continue
!
!     backwards substitution
!
      do 500 l = 1, lmm1
      k = lm - l
      nd = k + lmm1
      do 520 i = 1, im
      y(i,k) = y(i,k) - work(i,nd)*y(i,k+1)
  520 continue
  500 continue
      return
      end
