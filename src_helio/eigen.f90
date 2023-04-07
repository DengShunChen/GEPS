      subroutine eigen(mx,nn,eval,evec,epos,wk)
!
!  purpose: find the eigenvecotrs and eigenvalues of matrix 
!------------------------------------------------------------------
!  **** input *****
!  mx    : coefficient matrix  mx(nn,nn)
!  nn    : dimension of coefficient matrix
!  wk    : working aaray  wk(nn,nn)
!  **** output ****
!  eval  : eigen value,   eval(nn)
!  evec  : eigen vector,  evec(nn,nn)
!  epos  : transpose of eigen vector,  epos(nn,nn)
!------------------------------------------------------------------
      real mx(nn,nn),eval(nn),evec(nn,nn),epos(nn,nn),wk(nn,nn)
      dimension wr(nn),wi(nn),z(nn,nn),fv1(nn),iv1(nn),ipp(nn)
!
!  find the eigenvalue and eigenvector
!--------
      call eigrs(mx,nn,12,eval,evec,nn,wk,ier)
!---------
!      call rg (nn,nn,mx,wr,wi,1,z,iv1,fv1,ierr)
!fuji      call deig1(mx,nn,nn,0,wr,wi,z,fv1,ierr)
!      call indexx (nn,wr,ipp)
!sun-no: i.e. cray only  ( make eigval be decending order)
!      do 40 j = 1, nn
!      jj = nn -j + 1
!      eval(jj)   = wr(ipp(j))
!      do 40 i = 1, nn
!      evec(i,jj) = z(i,ipp(j))
!   40 continue
! ------
!  find the inverse matrix of eigenvector matrix
!
      do 100 i=1,nn
        do 100 j=1,nn
          epos(j,i)=evec(i,j)
 100  continue
!
      return
      end
