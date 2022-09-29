      subroutine lightning_ec                                           &
                 (nxj,klon,klev,ptu,pqu,ztenh,zqenh                     &
                 ,plu,rho,kcbot,kctop,pgeo,papn,pf                  &
                 ,ft,lndj,kuo)

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
      use physcons, only:vtmpc1 => con_fvirt, g => con_g
      use const,    only:RTYPE

      implicit none
      integer jl,jk,klev,klon,nxj
      integer kcbot(klon),kctop(klon)
      integer lndj(klon),kuo(klon)
      real    ptu(klon,klev)  ,pqu(klon,klev)                   &
             ,ztenh(klon,klev),zqenh(klon,klev)                 &
             ,papn(klon,klev)  ,pf(klon,klev)                   &
             ,rho(klon,klev),plu(klon,klev), pap(klon,klev)
      real(kind=RTYPE) pgeo(klon,klev)
      real    qgraup(klon,klev),qsnow(klon,klev)
      real,parameter::  vgraup= 3.0, &     ! typical fall speed for graupel (ms^-1)
                        vsnow = 0.5, &     ! typical fall speed for snow (ms^-1)
                        alpha = 32.4
      real    beta,zdz,zdp
      real    charg(klon),           &     !charging rate
              ft(klon),              &     !lightning flash density (flashes km^-2 day^-1)
              cape(klon)
!-------------------------------------------------------------------
      do jk = 1,klev
       do jl = 1,nxj
         if (lndj(jl) .eq. 1) then
             beta = 0.7    !over land
          else
             beta = 0.45   !over ocean or ice
         end if
         pap(jl,jk) = papn(jl,jk)*100. !mb to Pa
         qgraup(jl,jk) = beta*pf(jl,jk)/(rho(jl,jk)*vgraup)
         qsnow(jl,jk) = (1-beta)*pf(jl,jk)/(rho(jl,jk)*vsnow)
        end do
       end do

       charg = 0.
       cape = 0.
       ft = 0.
       do jl = 1,nxj
        if ( kuo(jl) .eq. 1 )then
        do jk = 1,klev-1      

          if (ztenh(jl,jk) .ge. 248.15 .and. ztenh(jl,jk) .le. 273.15) then
            zdz = (pgeo(jl,jk)-pgeo(jl,jk+1))/g             
            charg(jl) = charg(jl) + (qgraup(jl,jk)*                     &
                        (plu(jl,jk) + qsnow(jl,jk)))*rho(jl,jk)*zdz
          end if

          if (jk .le. kcbot(jl) .and. jk .gt. kctop(jl)) then
            zdp = pap(jl,jk) - pap(jl,jk-1)
            cape(jl) = cape(jl) +                                       &
                 ((ptu(jl,jk)-ztenh(jl,jk))/ztenh(jl,jk) +              &
                  vtmpc1*(pqu(jl,jk)-zqenh(jl,jk))-plu(jl,jk))*zdp
          end if
        end do
        cape(jl) = max(cape(jl),0.)
        ft(jl) = alpha*charg(jl)*sqrt(cape(jl))*                        & 
                 (min(pgeo(jl,kcbot(jl))/g*1.0e-3,1.8))**2
         end if
       end do

!       if (myrank .eq. 0) print*,"maxft=",maxval(ft),   &
!                                 "minft=",minval(ft)
!       if (myrank .eq. 0) print*,"maxcape=",maxval(cape),   &
!                                 "mincape=",minval(cape)
      return
      end subroutine lightning_ec
