      subroutine lightning_ec_gpu &
         (nxjp, klon, klev, ptu, pqu, ztenh, zqenh &
          , plu, rho, kcbot, kctop, pgeo, papn, pf &
          , ft, lndj, kuo)

!reference:
!Lopez P., 2016: A lightning parameterization for the ECMWF Integrated Forecasting System. Mon. Wea. Rev., 144, 3057-3075.

!input-------------------------------------------
!     nxj:
!     klon:
!     klev:
!     ptu: cloud temperature (K)
!     pqu: cloud specific humidity (kg/kg)
!     ztenh: environment temperature (K)
!     zqenh: environment specific humidity (kg/kg)
!     plu:convective cloud water (kg/kg)
!     kcbot: cloud base level
!     kctop: cloud top level
!     pgeo: geopotential height (m)
!     papn: pressure field (hPa)
!     pf: convective snow flux (kgm^-2s^-1))
!     rho: air parcel density (kgm-3)
!output------------------------------------------
!     ft: flash density (km^2day^-1)

!      use mpe
!      use rank
         use physcons, only: vtmpc1 => con_fvirt, g => con_g
         use const, only: RTYPE
         use param, only: my_max, my
         use index, only: jlistnum, jlist1

         implicit none
         integer jl, jk, klev, klon, j, jj, nxjp(my), myim(my_max)
         integer kcbot(klon, my_max), kctop(klon, my_max)
         integer lndj(klon, my_max), kuo(klon, my_max)
         real ptu(klon, klev, my_max), pqu(klon, klev, my_max) &
            , ztenh(klon, klev, my_max), zqenh(klon, klev, my_max) &
            , papn(klon, klev, my_max), pf(klon, klev, my_max) &
            , rho(klon, klev, my_max), plu(klon, klev, my_max), pap(klon, klev, my_max)
         real(kind=RTYPE) pgeo(klon, klev, my_max)
         real qgraup(klon, klev, my_max), qsnow(klon, klev, my_max)
         real, parameter::  vgraup = 3.0, &     ! typical fall speed for graupel (ms^-1)
                           vsnow = 0.5, &     ! typical fall speed for snow (ms^-1)
                           alpha = 32.4
         real beta, zdz, zdp
         real charg(klon, my_max), &     !charging rate
            ft(klon, my_max), &     !lightning flash density (flashes km^-2 day^-1)
            cape(klon, my_max)
         integer :: async_id = 1
!-------------------------------------------------------------------
         !$acc enter data create(pap, qgraup, qsnow, charg, cape, myim) async(async_id)
         !$acc parallel loop async(async_id) private(j)
         do jj = 1, jlistnum
            j = jlist1(jj)
            myim(jj) = nxjp(j)
         end do
         
         !$acc parallel loop gang collapse(2) async(async_id) private(beta)
         do jj = 1, jlistnum
            do jk = 1, klev
               !$acc loop vector
               do jl = 1, myim(jj)
                  if (lndj(jl, jj) .eq. 1) then
                     beta = 0.7    !over land
                  else
                     beta = 0.45   !over ocean or ice
                  end if
                  pap(jl, jk, jj) = papn(jl, jk, jj)*100. !mb to Pa
                  qgraup(jl, jk, jj) = beta*pf(jl, jk, jj)/(rho(jl, jk, jj)*vgraup)
                  qsnow(jl, jk, jj) = (1 - beta)*pf(jl, jk, jj)/(rho(jl, jk, jj)*vsnow)
               end do
            end do
         end do

         !$acc parallel loop collapse(2) async(async_id)
         do jj = 1, jlistnum
            do jl = 1, klon
               if (jl .le. myim(jj)) then
                  charg(jl, jj) = 0.
                  cape(jl, jj) = 0.
                  ft(jl, jj) = 0.
               end if
            end do
         end do
         
         !$acc parallel loop gang collapse(2) private(zdz) async(async_id)
         do jj = 1, jlistnum
            do jl = 1, klon
               if (jl .le. myim(jj)) then
                  if (kuo(jl, jj) .eq. 1) then
                     !$acc loop vector
                     do jk = 1, klev - 1

                        if (ztenh(jl, jk, jj) .ge. 248.15 .and. ztenh(jl, jk, jj) .le. 273.15) then
                           zdz = (pgeo(jl, jk, jj) - pgeo(jl, jk + 1, jj))/g
                           !$acc atomic
                           charg(jl, jj) = charg(jl, jj) + (qgraup(jl, jk, jj)* &
                                                    (plu(jl, jk, jj) + qsnow(jl, jk, jj)))*rho(jl, jk, jj)*zdz
                        end if

                        if (jk .le. kcbot(jl, jj) .and. jk .gt. kctop(jl, jj)) then
                           zdp = pap(jl, jk, jj) - pap(jl, jk - 1, jj)
                           !$acc atomic
                           cape(jl, jj) = cape(jl, jj) + &
                                      ((ptu(jl, jk, jj) - ztenh(jl, jk, jj))/ztenh(jl, jk, jj) + &
                                       vtmpc1*(pqu(jl, jk, jj) - zqenh(jl, jk, jj)) - plu(jl, jk, jj))*zdp
                        end if
                     end do
                  cape(jl, jj) = max(cape(jl, jj), 0.)
                  ft(jl, jj) = alpha*charg(jl, jj)*sqrt(cape(jl, jj))* &
                           (min(pgeo(jl, kcbot(jl, jj), jj)/g*1.0e-3, 1.8))**2
                  end if
               end if
            end do
         end do
         !$acc exit data delete(pap, qgraup, qsnow, charg, cape, myim) async(async_id)

!       if (myrank .eq. 0) print*,"maxft=",maxval(ft),   &
!                                 "minft=",minval(ft)
!       if (myrank .eq. 0) print*,"maxcape=",maxval(cape),   &
!                                 "mincape=",minval(cape)
         return
      end subroutine lightning_ec_gpu
