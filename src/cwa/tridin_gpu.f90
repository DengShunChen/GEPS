      subroutine tridin_gpu(l, n, nt, cl, cm, cu, r1, r2, au, a1, a2, async_id)
!sela %include dbtridin;
!c
         use machine, only: kind_phys
         use rank, only: myrank
         implicit none
         integer k, n, l, i, ix, kk, nt, async_id
         real(kind=kind_phys) fk, clr, aur
!c
         real(kind=kind_phys) cl(l, n), cm(l, n), cu(l, n), r1(l, n), r2(l, n, nt), &
            au(l, n), a1(l, n), a2(l, n, nt)
!-----------------------------------------------------------------------

         !$acc parallel loop gang vector private(i,fk,k) async(async_id)
         do i = 1, l
            fk = 1./cm(i, 1)
            au(i, 1) = fk*cu(i, 1)
            a1(i, 1) = fk*r1(i, 1)
            !$acc loop seq
            do kk = 1, nt
               a2(i, 1, kk) = fk*r2(i, 1, kk)
            end do

            !$acc loop seq
            do k = 2, n - 1
               clr = cl(i, k)
               fk = 1./(cm(i, k) - clr*au(i, k - 1))
               au(i, k) = fk*cu(i, k)
               a1(i, k) = fk*(r1(i, k) - clr*a1(i, k - 1))
               !$acc loop seq
               do kk = 1, nt
                  a2(i, k, kk) = fk*(r2(i, k, kk) - clr*a2(i, k - 1, kk))
               end do
            end do

            clr = cl(i, n)
            fk = 1./(cm(i, n) - clr*au(i, n - 1))
            a1(i, n) = fk*(r1(i, n) - clr*a1(i, n - 1))
            !$acc loop seq
            do kk = 1, nt
               a2(i, n, kk) = fk*(r2(i, n, kk) - clr*a2(i, n - 1, kk))
            end do

            !$acc loop seq
            do k = n - 1, 1, -1
               aur = au(i, k)
               a1(i, k) = a1(i, k) - aur*a1(i, k + 1)
               !$acc loop seq
               do kk = 1, nt
                  a2(i, k, kk) = a2(i, k, kk) - aur*a2(i, k + 1, kk)
               end do
            end do
         end do

!-----------------------------------------------------------------------
         return
      end

