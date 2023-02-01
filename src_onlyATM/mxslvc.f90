      subroutine mxslvc ( aa, mn, mt, kk )
!
!###############################################################
!
!  1. function specification
!
!    solve a matrix of a band-structured type coming from
!    two-stream approximations with cloud overlapping treatment
!
!  2. input/output variables :
!
!    aa (in) : matrix of coefficients and r.h.s. forcing
!              1st-5th  columns = elements left to the diagonal
!                  6th  column  = diagonal elements
!             7th-11th  columns = elements right to the diagonal
!                 12th  column  = r.h.s. forcings
!                                 need change two lows each other
!    aa (out): 6th-11th columns = upper tri-diagonal matrix
!                  12th colimn  = solutions of the matrix
!     mn   : horizontal dimension of aa
!     mt   : number of matrix sets (horizontal points) to be updated
!     kk   : vertical dimension of aa, 
!            4 * lev ,  lev : number of model layers
!
!  3. date & author
!
!     c. t. fong   july, 1996
!     modify to f90 by C-H Lee and sort by River Chen in 2015
!
!#####################################################################
!

      use fj_pad

      implicit  none

      integer   mn, mt, kk
      real      aa(mn+npad,kk,12)
!
!  local working arrays
!
      integer   lev,k,kini,i,m,km,jm,j,ln,kl,jl,jj
      real      vmin,tmp,zzz

      lev = kk / 4
      vmin = 1.0e-5  
!
!    eliminate all elements below the diagonal
!
      if ( lev .eq. 1 ) then
      print*,' only one model layer not allowed, job stop '
      stop
      end if
!
      do 500 k = 2, lev-1
      kini = (k-1)*4 
!
!ocl norecurrence
       do i = 1, mt
          aa(i,kini+1,6)   = aa(i,kini+1,6) - aa(i,kini+1,4)*aa(i,kini-1,8) &
               - aa(i,kini+1,5)*aa(i,kini,7)
          aa(i,kini+1,7) = aa(i,kini+1,7) - aa(i,kini+1,4)*aa(i,kini-1,9)   &
               - aa(i,kini+1,5)*aa(i,kini,8)
          aa(i,kini+1,12)   = aa(i,kini+1,12) - aa(i,kini+1,4)*aa(i,kini-1,12) &
               - aa(i,kini+1,5)*aa(i,kini,12)
          aa(i,kini+1,5) = 0.
          aa(i,kini+1,4) = 0.
          aa(i,kini+2,5)   = aa(i,kini+2,5) - aa(i,kini+2,3)*aa(i,kini-1,8) &
               - aa(i,kini+2,4)*aa(i,kini,7)
          aa(i,kini+2,6) = aa(i,kini+2,6) - aa(i,kini+2,3)*aa(i,kini-1,9) &
               - aa(i,kini+2,4)*aa(i,kini,8)
          aa(i,kini+2,12)   = aa(i,kini+2,12) - aa(i,kini+2,3)*aa(i,kini-1,12) &
               - aa(i,kini+2,4)*aa(i,kini,12)
          aa(i,kini+2,4) = 0.
          aa(i,kini+2,3) = 0.
          aa(i,kini+3,4)   = aa(i,kini+3,4) - aa(i,kini+3,2)*aa(i,kini-1,8) &
               - aa(i,kini+3,3)*aa(i,kini,7)
          aa(i,kini+3,5) = aa(i,kini+3,5) - aa(i,kini+3,2)*aa(i,kini-1,9) &
               - aa(i,kini+3,3)*aa(i,kini,8)
          aa(i,kini+3,12)   = aa(i,kini+3,12) - aa(i,kini+3,2)*aa(i,kini-1,12) &
               - aa(i,kini+3,3)*aa(i,kini,12)
          aa(i,kini+3,3) = 0.
          aa(i,kini+3,2) = 0.
          aa(i,kini+4,3)   = aa(i,kini+4,3) - aa(i,kini+4,1)*aa(i,kini-1,8) &
               - aa(i,kini+4,2)*aa(i,kini,7)
          aa(i,kini+4,4) = aa(i,kini+4,4) - aa(i,kini+4,1)*aa(i,kini-1,9) &
               - aa(i,kini+4,2)*aa(i,kini,8)
          aa(i,kini+4,12)   = aa(i,kini+4,12) - aa(i,kini+4,1)*aa(i,kini-1,12) &
               - aa(i,kini+4,2)*aa(i,kini,12)
          aa(i,kini+4,2) = 0.
          aa(i,kini+4,1) = 0.

          tmp = 1.d0/aa(i,kini+2,6)

!ibm   zzz = aa(i,km,7) / aa(i,km+1,6)
          zzz = aa(i,kini+1,7)*tmp
          aa(i,kini+1,6)  = aa(i,kini+1,6)  - aa(i,kini+2,5)*zzz
          aa(i,kini+1,10) = aa(i,kini+1,10) - aa(i,kini+2,9)*zzz
          aa(i,kini+1,11) = aa(i,kini+1,11) - aa(i,kini+2,10)*zzz
          aa(i,kini+1,12) = aa(i,kini+1,12) - aa(i,kini+2,12)*zzz
          aa(i,kini+1,7)  = 0.
!
          zzz  = 1.d0/aa(i,kini+1,6)
          aa(i,kini+1,10) = aa(i,kini+1,10)*zzz
          aa(i,kini+1,11) = aa(i,kini+1,11)*zzz
          aa(i,kini+1,12) = aa(i,kini+1,12)*zzz
          aa(i,kini+1,6) = 1.0
!
          aa(i,kini+2,9) = aa(i,kini+2,9) - aa(i,kini+2,5)*aa(i,kini+1,10)
          aa(i,kini+2,10) = aa(i,kini+2,10) - aa(i,kini+2,5)*aa(i,kini+1,11)
          aa(i,kini+2,12) = aa(i,kini+2,12) - aa(i,kini+2,5)*aa(i,kini+1,12)
          aa(i,kini+2,5) = 0.

          aa(i,kini+3,8) = aa(i,kini+3,8) - aa(i,kini+3,4)*aa(i,kini+1,10)
          aa(i,kini+3,9) = aa(i,kini+3,9) - aa(i,kini+3,4)*aa(i,kini+1,11)
          aa(i,kini+3,12) = aa(i,kini+3,12) - aa(i,kini+3,4)*aa(i,kini+1,12)
          aa(i,kini+3,4) = 0.

          aa(i,kini+4,7) = aa(i,kini+4,7) - aa(i,kini+4,3)*aa(i,kini+1,10)
          aa(i,kini+4,8) = aa(i,kini+4,8) - aa(i,kini+4,3)*aa(i,kini+1,11)
          aa(i,kini+4,12) = aa(i,kini+4,12) - aa(i,kini+4,3)*aa(i,kini+1,12)
          aa(i,kini+4,3) = 0.
!
!ibm
          zzz = 1.d0/aa(i,kini+2,6)
          aa(i,kini+2,9) = aa(i,kini+2,9)*zzz
          aa(i,kini+2,10) = aa(i,kini+2,10)*zzz
          aa(i,kini+2,12) = aa(i,kini+2,12)*zzz
          aa(i,kini+2,6) = 1.0
          
          aa(i,kini+3,8) = aa(i,kini+3,8) - aa(i,kini+3,5)*aa(i,kini+2,9)
          aa(i,kini+3,9) = aa(i,kini+3,9) - aa(i,kini+3,5)*aa(i,kini+2,10)
          aa(i,kini+3,12) = aa(i,kini+3,12) - aa(i,kini+3,5)*aa(i,kini+2,12)
          aa(i,kini+3,5) = 0.

          aa(i,kini+4,7) = aa(i,kini+4,7) - aa(i,kini+4,4)*aa(i,kini+2,9)
          aa(i,kini+4,8) = aa(i,kini+4,8) - aa(i,kini+4,4)*aa(i,kini+2,9+1)
          aa(i,kini+4,12) = aa(i,kini+4,12) - aa(i,kini+4,4)*aa(i,kini+2,12)
          aa(i,kini+4,4) = 0.
          
       enddo
!
  500 continue

!
!  when k = lev , only 6x4 matrix left
!
       kini = (lev-1) * 4
!ocl norecurrence,prefetch
       do 510  m = 1, 4     !  m : 1 ~ 4 ( diagonal increasing )
       km = kini + m
       jm = 7 - m
       do 510 i = 1, mt
       aa(i,km,jm)   = aa(i,km,jm) - aa(i,km,jm-2)*aa(i,kini-1,8) &
                                   - aa(i,km,jm-1)*aa(i,kini,7)
       aa(i,km,jm+1) = aa(i,km,jm+1) - aa(i,km,jm-2)*aa(i,kini-1,9) &
                                     - aa(i,km,jm-1)*aa(i,kini,8)
       aa(i,km,12)   = aa(i,km,12) - aa(i,km,jm-2)*aa(i,kini-1,12) &
                                   - aa(i,km,jm-1)*aa(i,kini,12)
       aa(i,km,jm-1) = 0.
       aa(i,km,jm-2) = 0.
  510  continue
!
!  only 4x4 matrix left
!
       do 650 m = 1, 3
       km = kini + m
       jm = 10 - m
!ibm
!fj       call vrec(tmp1(1),aa(1,km,6),mt)
!ibm
!ocl norecurrence,prefetch
	 do 610 j = 7, jm
         do 610 i = 1, mt
!ibm     aa(i,km,j) = aa(i,km,j) / aa(i,km,6)
!fj         aa(i,km,j) = aa(i,km,j)*tmp1(i)
         aa(i,km,j) = aa(i,km,j)*(1.d0/aa(i,km,6))

  610    continue

	 do 615 i = 1, mt
!ibm     aa(i,km,12) = aa(i,km,12) / aa(i,km,6)
!fj         aa(i,km,12) = aa(i,km,12)*tmp1(i)
         aa(i,km,12) = aa(i,km,12)*(1.d0/aa(i,km,6))
	 aa(i,km,6) = 1.0
  615    continue
!
!ocl norecurrence,prefetch
         do 640 ln = m+1, 4
         kl = kini + ln
!
         jl = 6 - ( ln - m ) 
         jm = 10 - ln
	 do 620 j = jl+1, jm
	 jj = 6 + (j-jl)
         do 620 i = 1, mt
         aa(i,kl,j) = aa(i,kl,j) - aa(i,kl,jl)*aa(i,km,jj)
  620    continue
!
         do 630 i = 1, mt
         aa(i,kl,12) = aa(i,kl,12) - aa(i,kl,jl)*aa(i,km,12)
         aa(i,kl,jl) = 0.
  630    continue
  640    continue
  650  continue
!
!ibm
!fj       call vrec(tmp1(1),aa(1,kk,6),mt)
!ibm
!ocl norecurrence,prefetch
       do 680 i = 1, mt
!ibm   aa(i,kk,12) = aa(i,kk,12) / aa(i,kk,6)
!fj       aa(i,kk,12) = aa(i,kk,12)*tmp1(i)
       aa(i,kk,12) = aa(i,kk,12)*(1.d0/aa(i,kk,6))
  680  continue  
!
!  back substitution
!
!
!  k = lev 
!
      kini = lev*4 + 1
!
      do 720  m = 1, 3
      km = kini - m
!
!ocl norecurrence,prefetch
       do 710 ln = m+1, 4
       kl = kini - ln
       jl = 6 + ( ln - m )
       do 710 i = 1, mt
       aa(i,kl,12) = aa(i,kl,12) - aa(i,kl,jl)*aa(i,km,12)
  710  continue
!
  720 continue
!
!  k = lev-1 ~ 1
!
      do 810 k = lev-1, 1, -1
      kini = k*4 + 1 
!
!ocl norecurrence,prefetch
      do 800 m = 1, 4
      km = kini - m
      jm = 6 + m
      do 800 i = 1, mt
      aa(i,km,12) = aa(i,km,12) - aa(i,km,jm)*aa(i,kini,12) &
                                - aa(i,km,jm+1)*aa(i,kini+1,12)
  800 continue
  810 continue
!
      return
      end
