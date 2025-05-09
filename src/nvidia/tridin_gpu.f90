      subroutine tridin_gpu(l, n, nt, cl, cm, cu, r1, au, a1, async_id)
!sela %include dbtridin;
!c
         use machine, only: kind_phys
         use rank, only: myrank
         implicit none
         integer k, n, l, i, ix, kk, nt, async_id
         real(kind=kind_phys) fk(l, n), clr, aur, au_t
!c
         real(kind=kind_phys) cl(l, n), cm(l, n), cu(l, n), r1(l, n, nt+1), &
            au(l, n), a1(l, n, nt+1)
!-----------------------------------------------------------------------

         !$acc enter data create(fk) async(async_id)
         !$acc parallel loop gang vector private(i, k, aur) async(async_id)
         do i = 1, l
            fk(i, 1) = 1./cm(i, 1)
            aur = fk(i, 1)*cu(i, 1)
            au(i, 1) = aur
            !$acc loop seq
            do k = 2, n - 1
               fk(i, k) = 1./(cm(i, k) - cl(i, k)*aur)
               aur = fk(i, k)*cu(i, k)
               au(i, k) = aur
            end do
            fk(i, n) = 1./(cm(i, n) - cl(i, n)*aur)
         end do


         
         !$acc parallel loop gang vector private(i,k,clr, aur) async(async_id)
         do i = 1, l
            clr = fk(i, 1)
            !$acc loop seq
            do kk = 1, nt+1
               a1(i, 1, kk) = clr*r1(i, 1, kk)
            end do
               !$acc loop seq
               do k = 2, n - 1
                  clr = cl(i, k)
                  aur = fk(i, k)
                  !$acc loop seq
                  do kk = 1, nt+1
                     a1(i, k, kk) = aur*(r1(i, k, kk) - clr*a1(i, k-1, kk))
                  end do
               end do
            clr = cl(i, n)
            aur = fk(i, n)
            !$acc loop seq
            do kk = 1, nt+1
               a1(i, n, kk) = aur*(r1(i, n, kk) - clr*a1(i, n-1, kk))
            end do
            
            !$acc loop seq
            do k = n - 1, 1, -1
               aur = au(i, k)
               !$acc loop seq
               do kk = 1, nt+1
                  a1(i, k, kk) = a1(i, k, kk) - aur*a1(i, k+1, kk)
               end do
            end do
            
         end do
         !$acc exit data delete(fk) async(async_id)
!-----------------------------------------------------------------------
         return
      end

