      subroutine glbmean (nx,my,my_max,cosl,x,gmean)
      use mpe
      use index
!
      implicit  none

      integer   nx,my,my_max
      real      cosl(my),x(nx,my_max),gmean

      integer   i,j,jj,nxj
      real      gsum(my),denom
!
      denom= 0.0
      do 100 j=1,my
        denom= denom + cosl(j)
  100 continue
      denom= 1.0/denom
!
      gmean= 0.0
      gsum=0.
      do 200 jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i = 1,nxj
          gsum(j)= gsum(j) + cosl(j)*x(i,jj)
        enddo
!ch     gsum(j) = gsum(j) / float(nxj)
        gsum(j) = gsum(j) / float(nxdef(j))
  200 continue
!

!ch   call mpe_unify(gsum,1,my,2,mpe_double)
      call mpe_unify(gsum,my,1,1,mpe_double)
!
      do j=1,my
        gmean= gmean+ gsum(j)
      enddo
      gmean= gmean*denom
!
      return
      end
!-----------------------------------------------
!ch for 2dMPI with array(nx_partial,my_partial)

      subroutine glbmean_2d (cosl,x,gmean)
      use mpe
      use param, only : my, my_max
      use index
!
      implicit  none

      real      cosl(my),x(nxp,my_max),gmean

      integer   i,j,jj,nxj
      real      gsum(my),denom
!
      denom= 0.0
      do 100 j=1,my
        denom= denom + cosl(j)
  100 continue
      denom= 1.0/denom
!
      gmean= 0.0
      gsum=0.
      do 200 jj = 1, jlistnum
        j=jlist1(jj)
        nxj=nxdef_2d(j)
        do i = 1,nxj
          gsum(j)= gsum(j) + cosl(j)*x(i,jj)
        enddo
        gsum(j) = gsum(j) / float(nxdef(j))
  200 continue
!
      call mpe_unify(gsum,my,1,1,mpe_double)
!
      do j=1,my
        gmean= gmean+ gsum(j)
      enddo
      gmean= gmean*denom
!
      return
      end
