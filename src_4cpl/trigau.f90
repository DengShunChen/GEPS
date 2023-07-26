      subroutine trigau ( aa, mn, mt, kk, kt, mode )
!
!##############################################################
!
! 1. function specification
!
!     solve a tri-diagonal system using gaussian elimination
!
! 2. common block
!
! 3. parameter specification
!
!   aa(in) : matrix of coefficients and r.h.s. forcing
!            1st column = elements left to the diagonal,
!            2nd column = diagonal elements,
!            3rd column = elements right to the diagonal,
!            4th column = r.h.s. forcings
!  aa(out) : 1, 2, 3, columns = upper tri-diagonal matrix (reformed),
!            4th column = solutions of the tri-diagonal system
!     mn   : horizontal dimension of aa
!     mt   : number of matrix sets (horizontal points) to be updated
!     kk   : vertical dimension of aa
!     kt   : number of equations (vertical levels) to be solved
!   mode   : calling mode, = 0, initial call for a new matrix
!                          = 1, subsequent calls for forcing change only
!  ** note: for mode = 0 the matrix is reformed to an upper tri-diag one
!
! 4. local veriable
!
! 5. calling modules
!
!    mixpbl,eeps
!
! 6. usage
!
!    call trigau (aa,mn,kk,kt,mode)
!
! 7. modules called
!
! 8. limitation
!
! 9. date & author
!
!      c.-s. liou,  november 1991
!
!#######################################################################
!
      implicit none

      integer  mn,mt,kk,kt,mode,i,k,kr,kt1

      real     aa(mn,kk,4)
!
!     eliminate all elements below the diagonal
!
      if ( mode .eq. 0 )  then
!ibm ...
         do i= 1,mt
           aa(i,1,2)=1.0/aa(i,1,2)
         enddo
!ibm

         do 100 k = 2, kt
         do 100 i = 1, mt
!ibm        aa(i,k,1) = aa(i,k,1) / aa(i,k-1,2)
!ibm        aa(i,k,2) = aa(i,k,2) - aa(i,k,1)*aa(i,k-1,3)
         aa(i,k,1) = aa(i,k,1) * aa(i,k-1,2)
         aa(i,k,2) = 1.0/(aa(i,k,2) - aa(i,k,1)*aa(i,k-1,3))
  100    continue
      endif
      do 120 k = 2, kt
      do 120 i = 1, mt
      aa(i,k,4) = aa(i,k,4) - aa(i,k,1)*aa(i,k-1,4)
  120 continue
!
!     back substitution
!
      do 200 i = 1, mt
!ibm     aa(i,kt,4) = aa(i,kt,4) / aa(i,kt,2)
      aa(i,kt,4) = aa(i,kt,4) * aa(i,kt,2)
  200 continue
      kt1 = kt - 1
      do 220 k = 1, kt1
      kr = kt1 - k + 1
         do 240 i = 1, mt
!ibm        aa(i,kr,4) = (aa(i,kr,4)-aa(i,kr,3)*aa(i,kr+1,4)) /aa(i,kr,2)
         aa(i,kr,4) = (aa(i,kr,4)-aa(i,kr,3)*aa(i,kr+1,4)) *aa(i,kr,2)
  240    continue
  220 continue
!
      return
      end
