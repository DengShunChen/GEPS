subroutine get_phi_gpu(nxjp, nxp, lev, ptop, cp, r, g, sgeo, pk, pk2, tt, qt, &
                       phii, phi)
   !
   use const, only: RTYPE
   use param, only: my_max, ncld
   use index, only: jlistnum, jlist1
   !
   ! input:
   !    sgeo:
   !    pk:
   !    pk2:
   !    tt:
   !    qt:
   ! output:
   !    phii
   !    phi:

   integer, intent(in):: nxjp(my), nxp, lev
   real, intent(in):: ptop, cp, r, g
   real(kind=RTYPE), dimension(nxp, my_max), intent(in):: sgeo
   real(kind=RTYPE), dimension(nxp, lev, my_max), intent(in) :: tt, pk, pk2
   real(kind=RTYPE), dimension(nxp, lev*ncld, my_max), intent(in) :: qt

   real, dimension(nxp, lev + 1, my_max) :: phii
   real(kind=RTYPE), dimension(nxp, lev, my_max) :: phi

   real :: ppd, ppp, ppu, ttv, dhgt, pkd, pk2d
   real, dimension(nxp, lev, my_max):: pk2x, dhgtz, theda

   integer i, k, kc, j, jj, nxj
   integer async_id

   async_id = 1
   !$acc enter data create(pk2x, dhgtz, theda) async(async_id)
   !
   ! geopotential height at model interface
   !
   !$acc parallel loop collapse(3) private(j, nxj, pk2d) async(async_id)
   do jj = 1, jlistnum
      do k = 1, lev
         do i = 1, nxp
            j = jlist1(jj)
            nxj = nxjp(j)
            if (i .le. nxj) then
               pk2d = pk2(i, k, jj)
               pk2d = log(pk2d)
               pk2d = pk2d*(cp/r)
               pk2x(i, k, jj) = exp(pk2d)
            end if
         end do
      end do
   end do

   !$acc parallel loop collapse(3) private(j, nxj, ppd, ppu, dhgt, pkd, ppp, ttv) async(async_id)
   do jj = 1, jlistnum
      do k = 1, lev
         do i = 1, nxp
            j = jlist1(jj)
            nxj = nxjp(j)
            if (i .le. nxj) then
               ppd = pk2x(i, k, jj)*1000.
               if (k .eq. 1) then
                  ppu = ptop
               else
                  ppu = pk2x(i, k - 1, jj)*1000.
               end if
               dhgt = (ppd - ppu)*100./g

               pkd = pk(i, k, jj)
               ppp = log(pkd)*(cp/r)
               ppp = exp(ppp)*1000.
               ttv = tt(i, k, jj)*(1.+0.608*qt(i, k, jj))
               !      ttv = tt(i,k, jj)
               dhgtz(i, k, jj) = dhgt*r*ttv/(100.*ppp)
               ! <<<< ===================================
               ! theda(i, k, jj) = tt(i, k, jj)*(1.0 + 0.608*qt(i, k, jj))/pk(i, k, jj)
               ! ========================================
               theda(i, k, jj) = ttv/pkd
               ! =================================== >>>>
            end if
         end do
      end do
   end do

   !$acc parallel loop collapse(2) private(j, nxj, kc) async(async_id)
   do jj = 1, jlistnum
      do i = 1, nxp
         j = jlist1(jj)
         nxj = nxjp(j)
         if (i .le. nxj) then
            phii(i, 1, jj) = 0.
            phii(i, 2, jj) = phii(i, 1, jj) + dhgtz(i, lev, jj)*g

            phi(i, lev, jj) = sgeo(i, jj) + cp*theda(i, lev, jj)*(pk2(i, lev, jj) - pk(i, lev, jj))

            !$acc loop seq
            do k = lev - 1, 1, -1
               kc = lev - k + 2
               phii(i, kc, jj) = phii(i, kc - 1, jj) + dhgtz(i, k, jj)*g

               !
               ! geopotential height at model layer
               !
               phi(i, k, jj) = phi(i, k + 1, jj) + cp*( &
                               theda(i, k, jj)*(pk2(i, k, jj) - pk(i, k, jj)) &
                               + theda(i, k + 1, jj)*(pk(i, k + 1, jj) - pk2(i, k, jj)))
            end do
         end if
      end do
   end do
   !$acc exit data delete(pk2x, dhgtz, theda) async(async_id)
   !
   return
end subroutine get_phi_gpu
