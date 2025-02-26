subroutine prexp_hybrid_cwb_gpu_refactor(nxdef_2d, nxp, lev, ptop, sigma, pt, pk, pk2, plt)
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
   !  plt:
   !
   ! *** present ***
   ! jlist1, nxdef_2d
   !
   ! **************************************************
   !
   use const, only: RTYPE
   use param, only: my, my_max
   use index, only: jlistnum, jlist1
   !
   implicit none
   real, intent(out):: plt(nxp, lev, my_max)
   real(kind=RTYPE), intent(out):: pk2(nxp, lev, my_max), &
                                   pk(nxp, lev, my_max)
   real, intent(in):: ptop
   real(kind=RTYPE), intent(in):: pt(nxp, my_max), sigma(lev + 1, 2)
   integer, intent(in):: nxdef_2d(my), nxp, lev
   real pk2_top, pk2_bot
   !
   real pkbot, pktop
   !
   !  compute  pressure variables
   !
   integer k, i, ktop, kbot, j, jj, nxj
   real capa, capap1, opok, ptopk
   real rt

   integer, parameter:: async_id = 1

   capa = 1.0/3.5
   capap1 = 1.0 + capa
   opok = 1.0/1000.0**capa
   !
   ptopk = ptop*opok*ptop**capa

   !$acc parallel loop collapse(2) async(async_id) &
   !$acc& private(j, nxj)
   do jj = 1, jlistnum
      do k = 1, lev
         j = jlist1(jj)
         nxj = nxdef_2d(j)
         if (k .gt. 1) then
            !$acc loop vector private(rt, pktop, pkbot, pk2_top, pk2_bot)
            do i = 1, nxj
               pktop = sigma(k, 1)*pt(i, jj) + sigma(k, 2) + ptop
               pk2_top = capa*log(pktop)
               pk2_top = opok*exp(pk2_top)
               pkbot = sigma(k + 1, 1)*pt(i, jj) + sigma(k + 1, 2) + ptop
               pk2_bot = capa*log(pkbot)
               pk2_bot = opok*exp(pk2_bot)

               rt = (pkbot*pk2_bot - pktop*pk2_top) &
                    /(capap1*(pkbot - pktop))

               pk2(i, k, jj) = pk2_bot
               pk(i, k, jj) = rt
               plt(i, k, jj) = 1000.0*rt**3.*sqrt(rt)
            end do
         else
            !$acc loop vector private(rt, pktop, pkbot, pk2_top, pk2_bot)
            do i = 1, nxj
               pkbot = sigma(k + 1, 1)*pt(i, jj) + sigma(k + 1, 2) + ptop
               pk2_bot = capa*log(pkbot)
               pk2_bot = opok*exp(pk2_bot)

               rt = (pkbot*pk2_bot - ptopk) &
                    /(capap1*(pkbot - ptop))

               pk2(i, k, jj) = pk2_bot
               pk(i, k, jj) = rt
               plt(i, k, jj) = 1000.0*rt**3.*sqrt(rt)
            end do
         end if
      end do
   end do

   !
   return
end subroutine prexp_hybrid_cwb_gpu_refactor
