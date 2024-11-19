!fpp$ noconcur r
!-----------------------------------------------------------------------
      subroutine tridin_gpu(ix,l,n,nt,cl,cm,cu,r1,r2,au,a1,a2,fk,fkk)
      !$acc routine vector
!sela %include dbtridi2;
!c
!     use machine     , only : kind_phys
      implicit none
      integer             is,k,kk,n,nt,l,i,ix
      real fk(ix)
!c
      real cl(ix,2:n), cm(ix,n), cu(ix,n-1),                               &
                           r1(ix,n),   r2(ix,n*nt),                       &
                           au(ix,n-1), a1(ix,n), a2(ix,n*nt),              &
                           fkk(ix,2:n-1)
!-----------------------------------------------------------------------
      !$acc loop vector private(i)
      do i=1,l
        fk(i)   = 1./cm(i,1)
        au(i,1) = fk(i)*cu(i,1)
        a1(i,1) = fk(i)*r1(i,1)
      enddo
      !$acc loop vector collapse(2) private(k,i,is)
      do k = 1, nt
        do i = 1, l
          is = (k-1) * n
          a2(i,1+is) = fk(i) * r2(i,1+is)
        enddo
      enddo
      !$acc loop vector private(i,k)
      do i=1,l
        !$acc loop seq
        do k=2,n-1
          fkk(i,k) = 1./(cm(i,k)-cl(i,k)*au(i,k-1))
          au(i,k)  = fkk(i,k)*cu(i,k)
          a1(i,k)  = fkk(i,k)*(r1(i,k)-cl(i,k)*a1(i,k-1))
        enddo
      enddo
      !$acc loop vector collapse(2) private(kk,i,k,is)
      do i=1,l
        do kk = 1, nt
          is = (kk-1) * n
          !$acc loop seq
          do k=2,n-1
            a2(i,k+is) = fkk(i,k)*(r2(i,k+is)-cl(i,k)*a2(i,k+is-1))
          enddo
        enddo
      enddo
      !$acc loop vector private(i)
      do i=1,l
        fk(i)   = 1./(cm(i,n)-cl(i,n)*au(i,n-1))
        a1(i,n) = fk(i)*(r1(i,n)-cl(i,n)*a1(i,n-1))
      enddo
      !$acc loop vector collapse(2) private(i,k,is)
      do k = 1, nt
        do i = 1, l
          is = (k-1) * n
          a2(i,n+is) = fk(i)*(r2(i,n+is)-cl(i,n)*a2(i,n+is-1))
        enddo
      enddo
      !$acc loop vector private(i,k)
      do i=1,l
        !$acc loop seq
        do k=n-1,1,-1
          a1(i,k) = a1(i,k) - au(i,k)*a1(i,k+1)
        enddo
      enddo
      !$acc loop vector collapse(2) private(i,kk,k,is)
      do i=1,l
        do kk = 1, nt
          is = (kk-1) * n
          !$acc loop seq
          do k=n-1,1,-1
            a2(i,k+is) = a2(i,k+is) - au(i,k)*a2(i,k+is+1)
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
      return
      end
