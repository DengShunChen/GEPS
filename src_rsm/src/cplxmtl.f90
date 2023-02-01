! 
      subroutine cplxmtl(a,b,nn,c) 
! 
!  purpose : multiplize two array which b is complex number (c=a*b) 
!------------------------------------------------------------------- 
!  **** input *****
!  a    : array of first 
!  b    : array of complex
!  nn   : dimension of array 
!  **** output ***** 
!  c    : array of aXb 
!--------------------------------------------------------------------  
      implicit  none
      integer   nn,i,j
      real      a(nn,nn),b(nn,2),c(nn,2)
      do 90 i=1,nn
        c(i,1)=0.
        c(i,2)=0.
 90   continue
      do 100 i=1,nn
      do 100 j=1,nn
        c(i,1)=c(i,1)+a(i,j)*b(j,1)
        c(i,2)=c(i,2)+a(i,j)*b(j,2)
 100  continue
!
      return 
      end 
