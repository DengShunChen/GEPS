subroutine coftrix_refactor(m, k, mx, nn, a, b, c, jtrun, lev, isym)
   !
   !  purpose: to create the symmetric and anti symmetric matrixes for NNMI
   !-----------------------------------------------------------------------
   !  **** input ****
   !  m     : horizontal wave number
   !  k     : vertical level
   !  nn    : dimension of cofficient matrix
   !  a     : coefficient array
   !  b     : coefficient array
   !  c     : coefficient array
   !  jtrun : number of horizontol mode (symmetric or asymmetric)
   !  lev   : total vertical levels
   !  isym  : = +1  : symmetric case
   !          =-1  : asymmetric case
   !  **** output ****
   !    mx    : coefficient matrix   (nn by nn) in size
   !
   ! modify to f90 in 2015 by C-H Lee
   ! sort in 2015 by River Chen
   !------------------------------------------------------------------------
   real mx(nn, nn)
   !      real a(jtrun,jtrun,lev),b(jtrun,jtrun,lev),c(jtrun,jtrun,lev)
   real a(jtrun, jtrun), b(jtrun, jtrun), c(jtrun, jtrun)
   do i = 1, nn
      do j = 1, nn
         mx(i, j) = 0.
      end do
   end do
   if (isym .eq. 1) then
      !
      ! symmterix case
      !
      if (mod(nn, 3) .eq. 0) then
         i = m
         jj = 1
         ni = int(nn/3)
         do j = 1, ni
            mx(jj, jj) = 0.
            !byl            mx(jj+1,jj+1)=c(m,i,k)
            !byl            mx(jj+2,jj+2)=c(m,i+1,k)
            mx(jj + 1, jj + 1) = c(m, i)
            mx(jj + 2, jj + 2) = c(m, i + 1)
            i = i + 2
            jj = jj + 3
         end do
         i = m
         jj = 1
         ni = int((nn - 1)/3)
         do j = 1, ni
            !byl            mx(jj,jj+1)=b(m,i,k)
            !byl            mx(jj+1,jj+2)=a(m,i+1,k)
            mx(jj, jj + 1) = b(m, i)
            mx(jj + 1, jj + 2) = a(m, i + 1)
            mx(jj + 2, jj + 3) = 0.
            i = i + 2
            jj = jj + 3
         end do
         !byl            mx(nn-2,nn-1)=b(m,i,k)
         !byl            mx(nn-1,nn)=a(m,i+1,k)
         mx(nn - 2, nn - 1) = b(m, i)
         mx(nn - 1, nn) = a(m, i + 1)
         i = m
         jj = 1
         ni = int((nn - 2)/3)
         do j = 1, ni
            mx(jj, jj + 2) = 0.
            mx(jj + 1, jj + 3) = 0.
            !byl            mx(jj+2,jj+4)=a(m,i+2,k)
            mx(jj + 2, jj + 4) = a(m, i + 2)
            i = i + 2
            jj = jj + 3
         end do
         mx(nn - 2, nn) = 0.
      else
         i = m
         jj = 1
         ni = int(nn/3)
         do j = 1, ni
            mx(jj, jj) = 0.
            !byl            mx(jj+1,jj+1)=c(m,i,k)
            !byl            mx(jj+2,jj+2)=c(m,i+1,k)
            mx(jj + 1, jj + 1) = c(m, i)
            mx(jj + 2, jj + 2) = c(m, i + 1)
            i = i + 2
            jj = jj + 3
         end do
         mx(nn - 1, nn - 1) = 0.
         !byl            mx(nn,nn)=c(m,i,k)
         mx(nn, nn) = c(m, i)
         i = m
         jj = 1
         ni = int((nn - 1)/3)
         do j = 1, ni
            !byl            mx(jj,jj+1)=b(m,i,k)
            !byl            mx(jj+1,jj+2)=a(m,i+1,k)
            mx(jj, jj + 1) = b(m, i)
            mx(jj + 1, jj + 2) = a(m, i + 1)
            mx(jj + 2, jj + 3) = 0.
            i = i + 2
            jj = jj + 3
         end do
         !byl            mx(nn-1,nn)=b(m,i,k)
         mx(nn - 1, nn) = b(m, i)
         i = m
         jj = 1
         ni = int((nn - 2)/3)
         do j = 1, ni
            mx(jj, jj + 2) = 0.
            mx(jj + 1, jj + 3) = 0.
            !byl            mx(jj+2,jj+4)=a(m,i+2,k)
            mx(jj + 2, jj + 4) = a(m, i + 2)
            i = i + 2
            jj = jj + 3
         end do
      end if
      !
      !  anti-symmetrix case
      !
   else if (isym .eq. -1) then
      if (mod(nn, 3) .eq. 0) then
         i = m
         jj = 1
         ni = int(nn/3)
         do j = 1, ni
            !byl            mx(jj,jj)=c(m,i,k)
            mx(jj, jj) = c(m, i)
            mx(jj + 1, jj + 1) = 0.
            !byl            mx(jj+2,jj+2)=c(m,i+1,k)
            mx(jj + 2, jj + 2) = c(m, i + 1)
            i = i + 2
            jj = jj + 3
         end do
         i = m
         jj = 1
         ni = int((nn - 1)/3)
         do j = 1, ni
            mx(jj, jj + 1) = 0.
            !byl            mx(jj+1,jj+2)=b(m,i+1,k)
            !byl            mx(jj+2,jj+3)=a(m,i+2,k)
            mx(jj + 1, jj + 2) = b(m, i + 1)
            mx(jj + 2, jj + 3) = a(m, i + 2)
            i = i + 2
            jj = jj + 3
         end do
         mx(nn - 2, nn - 1) = 0.
         !byl            mx(nn-1,nn)=b(m,i+1,k)
         mx(nn - 1, nn) = b(m, i + 1)
         i = m
         jj = 1
         ni = int((nn - 2)/3)
         do j = 1, ni
            !byl            mx(jj,jj+2)=a(m,i+1,k)
            mx(jj, jj + 2) = a(m, i + 1)
            mx(jj + 1, jj + 3) = 0.
            mx(jj + 2, jj + 4) = 0.
            i = i + 2
            jj = jj + 3
         end do
         !byl            mx(nn-2,nn)=a(m,i+1,k)
         mx(nn - 2, nn) = a(m, i + 1)
      else
         i = m
         jj = 1
         ni = int(nn/3)
         do j = 1, ni
            !byl            mx(jj,jj)=c(m,i,k)
            mx(jj, jj) = c(m, i)
            mx(jj + 1, jj + 1) = 0.
            !byl            mx(jj+2,jj+2)=c(m,i+1,k)
            mx(jj + 2, jj + 2) = c(m, i + 1)
            i = i + 2
            jj = jj + 3
         end do
         !byl            mx(nn,nn)=c(m,i,k)
         mx(nn, nn) = c(m, i)
         i = m
         jj = 1
         ni = int((nn - 1)/3)
         do j = 1, ni
            mx(jj, jj + 1) = 0.
            !byl            mx(jj+1,jj+2)=b(m,i+1,k)
            !byl            mx(jj+2,jj+3)=a(m,i+2,k)
            mx(jj + 1, jj + 2) = b(m, i + 1)
            mx(jj + 2, jj + 3) = a(m, i + 2)
            i = i + 2
            jj = jj + 3
         end do
         i = m
         jj = 1
         ni = int((nn - 2)/3)
         do j = 1, ni
            !byl            mx(jj,jj+2)=a(m,i+1,k)
            mx(jj, jj + 2) = a(m, i + 1)
            mx(jj + 1, jj + 3) = 0.
            mx(jj + 2, jj + 4) = 0.
            i = i + 2
            jj = jj + 3
         end do
         if (nn .ne. 1) then
            !byl            mx(nn-3,nn-1)=a(m,i+1,k)
            mx(nn - 3, nn - 1) = a(m, i + 1)
            mx(nn - 2, nn) = 0.
         end if
      end if
   end if
   !
   !   copy upper half to lower half
   !
   ! jj = 0
   ! do n = nn - 1, 1, -1
   !    jj = jj + 1
   !    do i = 1, n
   !       j = jj + i
   !       mx(j, i) = mx(i, j)
   !    end do
   ! end do
   !      call smetrix(mx,nn)
   !
   return
end subroutine coftrix_refactor
