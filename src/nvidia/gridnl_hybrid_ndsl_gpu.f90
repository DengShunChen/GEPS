subroutine gridnl_hybrid_ndsl_gpu(nxj, nx, lev, ncld &
                                  , cp, radsq, ut, vt, rdiv, tt, qt, phi, pt &
                                  , dtpl, dlpl, sinl, pk, pk2, dsigma, sigma, onocos &
                                  , cor, diveng, vdmerd, vdzonl, pten, deldm, sdpbl &
                                  , sd, pdot, vvel, sgeo)
   !$acc routine vector
   !$acc routine(vstruc_hybrid_cwb_gpu) vector
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
!  pten: horizontal adv of terrain pressure
!
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
   use openacc
   use cudafor

   implicit none

   integer nxj, nx, lev, ncld
   real cp, radsq, cpr2
   real(kind=RTYPE) onocos, cor, sinl

   real(kind=RTYPE) diveng(nx, lev), vdmerd(nx, lev), vdzonl(nx, lev), &
      pdot(nx, lev + 1), pten(nx, lev), dlpl(nx), dtpl(nx), &
      rdiv(nx, lev), ut(nx, lev), vt(nx, lev), tt(nx, lev), &
      qt(nx, lev*ncld), phi(nx, lev), pt(nx), sgeo(nx), &
      deldm(nx), spal(nx, lev), sd(nx, lev), sdpbl(nx), &
      dsigma(lev, 2), sigma(lev + 1, 2), pk(nx, lev), pk2(nx, lev), &
      cg(nx, lev), vvel(nx, lev)
   real(kind=RTYPE) odpsigt, thatt, spalt

   logical flag

   integer k, i, kbgn, kk, kkp1
   real px, px_pbl

   px_pbl = 850.
   kbgn = (lev*2)/3 - 1
   cpr2 = cp/radsq

   !$acc loop vector private(flag, px, odpsigt, thatt, spalt)
   do i = 1, nxj
      sd(i, 1) = 0.
      deldm(i) = 0.
      !$acc loop seq
      do k = 1, lev
         vvel(i, k) = ut(i, k)*dlpl(i)*onocos + vt(i, k)*dtpl(i)
         deldm(i) = deldm(i) - dsigma(k, 1)*vvel(i, k) - rdiv(i, k)*(dsigma(k, 2) + dsigma(k, 1)*pt(i))
         if (k .lt. lev) then
            sd(i, k + 1) = deldm(i)
         end if
      end do
      pdot(i, 1) = 0.0
      !$acc loop seq
      do k = 2, lev
         kk = lev - k + 2
         sd(i, k) = sd(i, k) - sigma(k, 1)*deldm(i)
         pdot(i, kk) = sd(i, k)
      end do
      pdot(i, lev + 1) = 0.0
      !$acc loop seq
      do k = 1, lev
         kk = lev - k + 1
         kkp1 = lev - k + 2
         vvel(i, k) = 0.5*((sigma(k, 1) + sigma(k + 1, 1))*(vvel(i, k) + deldm(i)) &
                           + (pdot(i, kk) + pdot(i, kkp1)))
      end do
      sdpbl(i) = sd(i, lev - 1)
      flag = .true.
      !$acc loop seq
      do k = kbgn, lev - 1
         px = sigma(k, 1)*pt(i) + sigma(k, 2)
         if (px .ge. px_pbl .and. flag) then
            sdpbl(i) = sd(i, k)
            flag = .false.
         end if
      end do
      odpsigt = 1.0/(dsigma(1, 2) + dsigma(1, 1)*pt(i))
      spalt = cpr2*tt(i, 1)*sigma(2, 1)*(pk2(i, 1) - pk(i, 1))*odpsigt
      vdmerd(i, 1) = -spalt*dtpl(i)/onocos - ut(i, 1)*cor &
                     - (ut(i, 1)*ut(i, 1) + vt(i, 1)*vt(i, 1))*onocos*sinl
      vdzonl(i, 1) = -spalt*dlpl(i) + vt(i, 1)*cor
      !$acc loop seq
      do k = 1, lev - 1
         odpsigt = 1.0/(dsigma(k + 1, 2) + dsigma(k + 1, 1)*pt(i))
         thatt = tt(i, k) - (tt(i, k) - tt(i, k + 1)) &
                 *(pk(i, k + 1) - pk2(i, k))/(pk(i, k + 1) - pk(i, k))
         phi(i, k) = thatt*(pk(i, k + 1) - pk(i, k))
         spalt = cpr2*tt(i, k + 1)*(sigma(k + 1, 1)*(pk(i, k + 1) - pk2(i, k)) &
                                    + sigma(k + 2, 1)*(pk2(i, k + 1) - pk(i, k + 1)))*odpsigt
         vdmerd(i, k + 1) = -spalt*dtpl(i)/onocos - ut(i, k + 1)*cor &
                            - (ut(i, k + 1)*ut(i, k + 1) + vt(i, k + 1)*vt(i, k + 1))*onocos*sinl
         vdzonl(i, k + 1) = -spalt*dlpl(i) + vt(i, k + 1)*cor
      end do
      phi(i, lev) = cp*tt(i, lev)*(pk2(i, lev) - pk(i, lev))
      diveng(i, lev) = sgeo(i) + phi(i, lev)
      !$acc loop seq
      do k = lev - 1, 1, -1
         phi(i, k) = phi(i, k + 1) + cp*phi(i, k)
         diveng(i, k) = sgeo(i) + phi(i, k)
      end do
   end do
   return
end
