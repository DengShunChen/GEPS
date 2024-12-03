subroutine vstruc_hybrid_cwb_gpu_refactor(nxjp, nxp, lev, cp, radsq, sigma, dsigma, &
                                          pt, tt, qt, pk, pk2, spal, phi, ncld)
   !
   !  subroutine to compute several intermediate pressure dependent
   !  variables and parameters
   !
   ! *** input ****
   !
   !  nx: e-w dimension no.
   !  lev: number of vertical levels
   !  pk: odd (full) level p**capa
   !  pk2: even( half) level p**capa
   !  cp: specific heat of air
   !  radsq: (rad of earth)**2
   !  dsig: sigma layer thicknesses
   !  sig: sigma levels
   !  pt: terrain pressure
   !  tt: virtual potential temperature
   !  qt: specific humidity
   !
   ! *** output ***
   !
   !  phi: geopotential
   !  spal: energy conversion term (sigma*ps*alpha)
   !  odpsig: reciprical of layer mass
   !
   ! **************************************************
   !
   use const, only: RTYPE
   use param, only: my, my_max
   use index, only: jlistnum, jlist1
   use openacc
   use cudafor
   !
   implicit none

   integer nxjp(my), nxp, lev, ncld

   real(kind=RTYPE), intent(in):: tt(nxp, lev, my_max), &
                                  qt(nxp, lev*ncld, my_max), pt(nxp, my_max), &
                                  dsigma(lev, 2), sigma(lev + 1, 2), &
                                  pk(nxp, lev, my_max), pk2(nxp, lev, my_max)
   real(kind=RTYPE), intent(out):: phi(nxp, lev, my_max), spal(nxp, lev, my_max)
   real(kind=RTYPE) odpsig(nxp, lev, my_max)
   real(kind=RTYPE) rt, pk_t, pk_nxt

   real cpr2, cp, radsq
   integer i, n, nk, k, kk, nxj, j, jj
   integer async_id
   async_id = 1
   !$acc enter data create(odpsig) async(async_id)
   !
   !  compute time dependent pressure variables
   !
   cpr2 = cp/radsq
   !
   !
   !  half-level specific humidity, interpolate in p**capa
   !
   !      do n=1,ncld
   !      nk=(n-1)*lev
   !      do k=1,lev-1
   !      kk=nk+k
   !      do i=1,nxj
   !      qhat(i,kk+1)= qt(i,kk+1)+(qt(i,kk)-qt(i,kk+1))     &
   !       *(pk(i,k+1)-pk2(i,k))/(pk(i,k+1)-pk(i,k))
   !      enddo
   !      enddo
   !      enddo
   !
   !$acc parallel loop collapse(3) private(j, nxj, rt, pk_t, pk_nxt) async(async_id)
   do jj = 1, jlistnum
      do k = 1, lev - 1
         do i = 1, nxp
            j = jlist1(jj)
            nxj = nxjp(j)
            if (i .le. nxj) then
               pk_t = pk(i, k, jj)
               pk_nxt = pk(i, k + 1, jj)
               if (k .eq. 1) then
                  odpsig(i, 1, jj) = 1.0/(dsigma(1, 2) + dsigma(1, 1)*pt(i, jj))
                  spal(i, 1, jj) = cpr2*tt(i, 1, jj)*sigma(2, 1)*(pk2(i, 1, jj) - pk(i, 1, jj))*odpsig(i, 1, jj)
               end if
               !
               !  recipricol of layer pressure depth
               !
               odpsig(i, k + 1, jj) = 1.0/(dsigma(k + 1, 2) + dsigma(k + 1, 1)*pt(i, jj))
               !
               !  half-level specific humidity, interpolate in p**capa
               !
               !
               !  half level potential temperature, defined as weighted combination
               !  of full level thicknesses, not an interpolation
               !
               rt = tt(i, k, jj) - (tt(i, k, jj) - tt(i, k + 1, jj)) &
                    *(pk_nxt - pk2(i, k, jj))/(pk_nxt - pk_t)
               !
               !  geopotential thicknesses
               !
               phi(i, k, jj) = rt*(pk_nxt - pk_t)
               !
               !  energy conversion term for terrain pressure contribution to
               !  horizontal pressure gradient
               !
               spal(i, k + 1, jj) = cpr2*tt(i, k + 1, jj)*(sigma(k + 1, 1)*(pk_nxt - pk2(i, k, jj)) &
                                                           + sigma(k + 2, 1)*(pk2(i, k + 1, jj) - pk_nxt)) &
                                    *odpsig(i, k + 1, jj)
            end if
         end do
      end do
   end do
   !
   !  hydrostatic equation: note that terrain geopotential is excluded.
   !  it it constant forcing term that is included when laplacian of
   !  geopotential is computed in divergence equation.
   !
   !$acc parallel loop collapse(2) private(j, nxj) async(async_id)
   do jj = 1, jlistnum
      do i = 1, nxp
         j = jlist1(jj)
         nxj = nxjp(j)
         if (i .le. nxj) then
            phi(i, lev, jj) = cp*tt(i, lev, jj)*(pk2(i, lev, jj) - pk(i, lev, jj))
            !$acc loop seq
            do k = lev - 1, 1, -1
               phi(i, k, jj) = phi(i, k + 1, jj) + cp*phi(i, k, jj)
            end do
         end if
      end do
   end do
   !$acc exit data delete(odpsig) async(async_id)
   !
   return
end subroutine vstruc_hybrid_cwb_gpu_refactor
