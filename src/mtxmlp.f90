      subroutine mtxmlp (a,b,c,lm)
!
!  matrix multiply routine
!
!  ***input***
!
!  a:  first array in matrix product
!  b:  second array in matrix product
!  lm: order of matrices
!
! ***output***
!
!  c: matrix product
!
! **************************************************************
!
      use const, only : RTYPE
!
      implicit  none
      integer   lm,j,k,i
      real(kind=RTYPE) a(lm,lm),b(lm,lm),c(lm,lm)
!
      call zilch(c,lm*lm)
!
      do 1 j=1,lm
      do 1 k=1,lm
      do 1 i=1,lm
    1 c(i,j)= c(i,j)+a(i,k)*b(k,j)
!
      return
      end
