      subroutine nnmi(no,x,epos,eval,evec,nn,c,d,bal,cutfreq)
!
       use rank
!
!  purpose : to get delta x vector at each iteration
!-------------------------------------------------------------------
!  **** input ****
!  x    : vector matrix of variable array( tendency )
!  epos : inverse matrix array of eigenvector matrix
!  eval : eigenval array of coeffienent matrix
!  evec : eigenvector matrix array of coefficient matrix 
!  nn   : dimension of array
!  c    : work array
!  d    : work array
!  cutfreq : cut of frequency
!  **** output ****
!  x    : vector matrix of variable array (correction)
!  bal  : convergence indicator
!
!  modify to f90 by C-H Lee and sort by River Chen in 2005
!
!--------------------------------------------------------------------
      implicit  none
      integer   no,nn
      real      x(nn,2),epos(nn,nn),evec(nn,nn),eval(nn)
      real      c(nn,2),d(nn,2),bal,cutfreq

      integer   i,j
!
! set cut of frequency ( now , freq=1.0 means 24 hours ) 
!
!      data freq/1.0/
!
!   matrix multiplication to get alpha, (5.54)
!
       c=0.
!
       call cplxmtl(epos,x,nn,c)
!
      do 110 i=1,nn
        d(i,1)=0.
        d(i,2)=0.
 110  continue
!
!   compute BAL, by (5.64)
!
      do 120 j=1,nn
        if(abs(eval(j)).ge.cutfreq)then
!        if(abs(eval(j)).le.cutfreq)then
          bal=bal+c(j,1)**2+c(j,2)**2
          d(j,1)=-c(j,2)/eval(j)
          d(j,2)= c(j,1)/eval(j)
        endif
 120  continue
!
!    transform alpha to get x, (5.55)
!
       call cplxmtl(evec,d,nn,x)
!
      return
      end
