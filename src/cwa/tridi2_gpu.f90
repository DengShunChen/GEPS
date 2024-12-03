      subroutine tridi2_gpu(ix,l,n,cl,cm,cu,r1,r2,au,a1,a2)
      !$acc routine vector
!sela %include dbtridi2;
!c
      use machine     , only : kind_phys
      use rank, only : myrank
      implicit none
      integer             k,n,l,i,ix
      real(kind=kind_phys) fk
!c
      real(kind=kind_phys) cl(ix,2:n),cm(ix,n),cu(ix,n-1),r1(ix,n),r2(ix,n), &
                au(ix,n-1),a1(ix,n),a2(ix,n)
!-----------------------------------------------------------------------      
      
      !$acc loop vector
      do i=1,l
        fk      = 1./cm(i,1)
        au(i,1) = fk*cu(i,1)
        a1(i,1) = fk*r1(i,1)
        a2(i,1) = fk*r2(i,1)
      enddo

      !$acc loop vector
      do i=1,l
        !$acc loop seq
        do k=2,n-1
          fk      = 1./(cm(i,k)-cl(i,k)*au(i,k-1))
          au(i,k) = fk*cu(i,k)
          a1(i,k) = fk*(r1(i,k)-cl(i,k)*a1(i,k-1))
          a2(i,k) = fk*(r2(i,k)-cl(i,k)*a2(i,k-1))
        enddo
      enddo
      !$acc loop vector
      do i=1,l
        fk      = 1./(cm(i,n)-cl(i,n)*au(i,n-1))
        a1(i,n) = fk*(r1(i,n)-cl(i,n)*a1(i,n-1))
        a2(i,n) = fk*(r2(i,n)-cl(i,n)*a2(i,n-1))
      enddo
      !$acc loop vector
      do i=1,l
        !$acc loop seq
        do k=n-1,1,-1
          a1(i,k) = a1(i,k)-au(i,k)*a1(i,k+1)
          a2(i,k) = a2(i,k)-au(i,k)*a2(i,k+1)
        enddo
      enddo
      
      

      
!-----------------------------------------------------------------------
      return
      end



