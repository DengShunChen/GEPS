      subroutine normalize(fxtotal,nx,x,y)

      implicit none

      integer  nx
      real     x(nx)
      real     y(nx)
      real     fxtotal,newsum

!
! new sum should be added
!
      newsum=fxtotal
!
      call reassignsum(x,nx,newsum,y)
!
      return
      end
!
!
      subroutine reassignsum(x,n,newsum,y)

      implicit  none
      integer   n,i
      real      x(n),y(n),ave,var
      real      newsum

      call avevar(x,n,ave,var)
      do i=1,n
        y(i)=x(i)-ave
      enddo
      newsum=newsum/n
      do i=1,n
        y(i)=y(i)+newsum
      enddo
      return
      end
!
!
      SUBROUTINE avevar(data,n,ave,var)

      implicit  none
      INTEGER   n
      REAL      ave,var,data(n)
      INTEGER   j
      REAL      s,ep

      ave=0.0
      do 11 j=1,n
        ave=ave+data(j)
11    continue
      ave=ave/n
      var=0.0
      ep=0.0
      do 12 j=1,n
        s=data(j)-ave
        ep=ep+s
        var=var+s*s
12    continue
      var=(var-ep**2/n)/(n-1)
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software $!6)$3D#21)8.
