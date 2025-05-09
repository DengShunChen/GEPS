      subroutine tridi2_gpu(l, n, cl, cm, cu, r1, au, a1, async_id)
!sela %include dbtridi2;
!c
         use machine, only: kind_phys
         use rank, only: myrank
         implicit none
         integer k, n, l, i, ix, async_id, kk
         real(kind=kind_phys) fk(l, n), clr, au_t
!c
         real(kind=kind_phys) cl(l, n), cm(l, n), cu(l, n), r1(l, n, 2), &
            au(l, n), a1(l, n, 2)
!-----------------------------------------------------------------------
         !$acc enter data create(fk) async(async_id)
         !$acc parallel loop gang vector private(i, k, au_t) async(async_id)
         do i = 1, l
            fk(i, 1) = 1./cm(i, 1)
            au_t = fk(i, 1)*cu(i, 1)
            au(i, 1) = au_t
            !$acc loop seq
            do k = 2, n - 1
               fk(i, k) = 1./(cm(i, k) - cl(i, k)*au_t)
               au_t = fk(i, k)*cu(i, k)
               au(i, k) = au_t
            end do
            fk(i, n) = 1./(cm(i, n) - cl(i, n)*au_t)
         end do
               
               
         !$acc parallel loop gang vector private(i,k,kk, clr, au_t) async(async_id)

         do i = 1, l
            clr = fk(i, 1)
            a1(i, 1, 1) = clr*r1(i, 1, 1)
            a1(i, 1, 2) = clr*r1(i, 1, 2)

            !$acc loop seq
            do k = 2, n - 1
               clr = fk(i, k)
               au_t = cl(i, k)
               a1(i, k, 1) = clr*(r1(i, k, 1) - au_t*a1(i, k - 1, 1))
               a1(i, k, 2) = clr*(r1(i, k, 2) - au_t*a1(i, k - 1, 2))
            end do
            clr = fk(i, n)
            au_t = cl(i, n)
            a1(i, n, 1) = clr*(r1(i, n, 1) - au_t*a1(i, n - 1, 1))
            a1(i, n, 2) = clr*(r1(i, n, 2) - au_t*a1(i, n - 1, 2))

            !$acc loop seq
            do k = n - 1, 1, -1
               clr = au(i, k)
               a1(i, k, 1) = a1(i, k, 1) - clr*a1(i, k + 1, 1)
               a1(i, k, 2) = a1(i, k, 2) - clr*a1(i, k + 1, 2)
            end do
         end do
         
         !$acc exit data delete(fk) async(async_id)
         
!-----------------------------------------------------------------------
         return
      end

