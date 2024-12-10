subroutine gridnl_hybrid_ndsl_gpu_refactor(nxjp, nxp, lev, ncld, cp, radsq, &
                                           ut, vt, rdiv, tt, qt, phi, pt, dtpl, dlpl, &
                                           sinl, pk, pk2, dsigma, sigma, onocos, &
                                           cor, diveng, vdmerd, vdzonl, pten, deldm, &
                                           sdpbl, sd, pdot, vvel, sgeo)
   !
   !  real to spectral transformation, compute non-linear contributions
   !  to spectral tendencies
   !
   !  **** input variables ****
   !
   !  nx:  e-w dimension
   !  lev: no. of vertical levels
   !  cp: specific heat of air
   !  radsq: radius of earth squared
   !  ut: gridpt e-w velocity (scaled)
   !  vt: gridpt n-s velocity (scaled)
   !  tt: gridpt virtual potential temp
   !  qt: gridpt specific humidity
   !  rdiv: gridpt divergence
   !  pt: gridpt terrain pressure
   !  dtpl: d(pt)/dy
   !  dlpl: d(pt)/dx
   !  pk: exner function on full(odd) levels
   !  pk2: exner function on half(even) levels
   !  dsigma: thickness of sigma layers
   !  sigma: sigma level values
   !  onocos: 1.0/(cos(lat)**2)
   !  cor: coriolis on each gaussian latitude
   !  pten: horizontal adv of terrain pressure (not used)
   !
   !  phi:
   !  sgeo:
   ! *** output variables ***
   !
   !  deldm: terrain pressure tendency
   !  sd: vertical velocity
   !  pdot: vertical velocity
   !  vvel: vertical velocity at mean layer(Pa/s)
   !  diveng: energy term of divergence equation
   !  vdmerd: meridional advection term of divergence/vorticity tends
   !  vdzonl: zonal advection term of divergence/vorticity tends
   !
   !  modify to f90 by C-H Lee and sort by River Chen in 2015
   !
   ! ******************************************************************
   !
   use const, only: RTYPE
   use param, only: my, my_max
   use index, only: jlistnum, jlist1
   use grid, only: latpart
   use openacc
   use cudafor
   !
   implicit none

   integer nxjp(my), nxp, lev, ncld
   real cp, radsq
   ! << input >>
   real(kind=RTYPE), dimension(my), intent(in):: onocos, cor, sinl
   real(kind=RTYPE), dimension(nxp, my_max), intent(in):: &
      pt, dtpl, dlpl, sgeo, pten
   real(kind=RTYPE), dimension(nxp, lev, my_max), intent(in):: &
      ut, vt, rdiv, tt, phi, pk, pk2
   real(kind=RTYPE), intent(in):: dsigma(lev, 2), &
                                  sigma(lev + 1, 2), &
                                  qt(nxp, lev*ncld, my_max)
   ! << output >>
   real(kind=RTYPE), dimension(nxp, my_max), intent(out):: &
      deldm, sdpbl
   real(kind=RTYPE), dimension(nxp, lev, my_max), intent(out):: &
      diveng, vdmerd, vdzonl, sd, vvel
   real(kind=RTYPE), dimension(nxp, lev + 1, latpart), intent(out):: &
      pdot

   real(kind=RTYPE), dimension(nxp, lev, my_max):: spal, cg
   !
   real(kind=RTYPE) deldm_t, ut_t, vt_t
   logical flag(nxp, my_max)
   !
   integer k, i, kbgn, kk, kkp1, jj, j, nxj
   real px, px_pbl
   !
   integer async_id, istat
   integer(kind=cuda_stream_kind) stream
   async_id = 1
   stream = acc_get_cuda_stream(async_id)
   !$acc enter data create(spal, cg, flag) async(async_id)

   !CWB2014 fixed undefined value problem in diabat line 665
   !$acc host_data use_device(sd, deldm, pdot)
   istat = cudaMemsetAsync(sd, real(0.0, RTYPE), size(sd), stream)
   istat = cudaMemsetAsync(deldm, real(0.0, RTYPE), size(deldm), stream)
   istat = cudaMemsetAsync(pdot, real(0.0, RTYPE), size(pdot), stream)
   !$acc end host_data
   !
   ! << vstruc_hybrid_cwb >>
   call vstruc_hybrid_cwb_gpu_refactor(nxjp, nxp, lev, cp, radsq, &
                                       sigma, dsigma, pt, tt, qt, pk, pk2, &
                                       spal, phi, ncld)
   !$acc parallel loop collapse(3) private(j, nxj, ut_t, vt_t) async(async_id)
   do jj = 1, jlistnum
      do k = 1, lev
         do i = 1, nxp
            j = jlist1(jj)
            nxj = nxjp(j)
            if (i .le. nxj) then
               ut_t = ut(i, k, jj)
               vt_t = vt(i, k, jj)

               cg(i, k, jj) = ut_t*dlpl(i, jj)*onocos(j) &
                              + vt_t*dtpl(i, jj)

               vdmerd(i, k, jj) = -spal(i, k, jj)*dtpl(i, jj)/onocos(j) &
                                  - ut_t*cor(j) &
                                  - (ut_t*ut_t + vt_t*vt_t)*onocos(j)*sinl(j)
               vdzonl(i, k, jj) = -spal(i, k, jj)*dlpl(i, jj) + vt_t*cor(j)

               diveng(i, k, jj) = sgeo(i, jj) + phi(i, k, jj)
            end if
         end do
      end do
   end do
   !
   !  surface pressure tendency
   !
   !$acc parallel loop collapse(2) private(j, nxj) async(async_id)
   do jj = 1, jlistnum
      do i = 1, nxp
         j = jlist1(jj)
         nxj = nxjp(j)
         if (i .le. nxj) then
            !$acc loop seq
            do k = 1, lev
               deldm(i, jj) = deldm(i, jj) - dsigma(k, 1)*cg(i, k, jj) &
                              - rdiv(i, k, jj)*(dsigma(k, 2) + dsigma(k, 1)*pt(i, jj))
               if (k .le. lev - 1) then
                  sd(i, k + 1, jj) = deldm(i, jj)
               end if
            end do
         end if
      end do
   end do
   !
   !  vertical velocity
   !
   !$acc parallel loop collapse(3) private(j, nxj, kk) async(async_id)
   do jj = 1, jlistnum
      do k = 2, lev
         do i = 1, nxp
            j = jlist1(jj)
            nxj = nxjp(j)
            if (i .le. nxj) then
               kk = lev - k + 2
               sd(i, k, jj) = sd(i, k, jj) - sigma(k, 1)*deldm(i, jj)
               pdot(i, kk, jj) = sd(i, k, jj)
            end if
         end do
      end do
   end do
   !$acc parallel loop collapse(3) private(j, nxj, kk, kkp1) async(async_id)
   do jj = 1, jlistnum
      do k = 1, lev
         do i = 1, nxp
            j = jlist1(jj)
            nxj = nxjp(j)
            if (i .le. nxj) then
               kk = lev - k + 1
               kkp1 = lev - k + 2

               vvel(i, k, jj) = 0.5*((sigma(k, 1) + sigma(k + 1, 1)) &
                                     *(cg(i, k, jj) + deldm(i, jj)) &
                                     + (pdot(i, kk, jj) + pdot(i, kkp1, jj)))
            end if
         end do
      end do
   end do
   !
   !  obtain vertical velocity within low layers
   !
   !$acc parallel loop collapse(2) private(j, nxj) async(async_id)
   do jj = 1, jlistnum
      do i = 1, nxp
         j = jlist1(jj)
         nxj = nxjp(j)
         if (i .le. nxj) then
            flag(i, jj) = .true.
            sdpbl(i, jj) = sd(i, lev - 1, jj)
         end if
      end do
   end do
   !
   px_pbl = 850.
   kbgn = (lev*2)/3 - 1
   !$acc parallel loop collapse(2) private(j, nxj) async(async_id)
   do jj = 1, jlistnum
      do i = 1, nxp
         j = jlist1(jj)
         nxj = nxjp(j)
         if (i .le. nxj) then
            !$acc loop seq
            do k = kbgn, lev - 1
               px = sigma(k, 1)*pt(i, jj) + sigma(k, 2)
               if (px .ge. px_pbl .and. flag(i, jj)) then
                  sdpbl(i, jj) = sd(i, k, jj)
                  flag(i, jj) = .false.
               end if
            end do
         end if
      end do
   end do
   !
   !
   !$acc exit data delete(spal, cg, flag) async(async_id)
   return
   !
end subroutine gridnl_hybrid_ndsl_gpu_refactor
