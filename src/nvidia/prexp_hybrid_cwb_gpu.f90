subroutine prexp_hybrid_cwb_gpu(nxj, nx, lev, ptop, sigma, pt, pk, pk2, plt)
   !$acc routine vector
!
!  subroutine to compute p to the kapa on odd and even levels
!
! *** input ****
!
!  nx: e-w dimension no.
!  lev: number of vertical levels
!  sig: sigma levels
!  pt: terrain pressure
!
! *** output ***
!
!  pk: odd (full) level p**capa
!  pk2: even( half) level p**capa
!
! **************************************************
!
   use const, only: RTYPE

   implicit none
   integer nxj, nx, lev
   real ptop
   real plt(nx, lev)
   real(kind=RTYPE) pt(nx), sigma(lev + 1, 2), pk2(nx, lev), pk(nx, lev)
   integer k, i
   real capa, capap1, opok, ptopk
   real rt, rbot, rtop

   capa = 1.0/3.5
   capap1 = 1.0 + capa
   opok = 1.0/1000.0**capa
   ptopk = ptop*opok*ptop**capa

   !$acc loop vector collapse(2) private(rt)
   do k = 1, lev
   do i = 1, nxj
      rt = sigma(k + 1, 1)*pt(i) + sigma(k + 1, 2) + ptop
      rt = log(rt)
      rt = capa*rt
      rt = exp(rt)
      pk2(i, k) = opok*rt
   end do
   end do

   !$acc loop vector private(rbot, rtop)
   do i = 1, nxj
      rbot = sigma(2, 1)*pt(i) + sigma(2, 2) + ptop
      pk(i, 1) = (rbot*pk2(i, 1) - ptopk)/(capap1*(rbot - ptop))
      plt(i, 1) = 1000.0*pk(i, 1)*pk(i, 1)*pk(i, 1)*sqrt(pk(i, 1))
      !$acc loop seq
      do k = 2, lev
         rtop = rbot
         rbot = sigma(k + 1, 1)*pt(i) + sigma(k + 1, 2) + ptop
         pk(i, k) = (rbot*pk2(i, k) - rtop*pk2(i, k - 1))/(capap1*(rbot - rtop))
         plt(i, k) = 1000.0*pk(i, k)*pk(i, k)*pk(i, k)*sqrt(pk(i, k))
      end do
   end do

   return
end
