    subroutine cududv_n(nxj,klon,klev,ktopm2,ktype,kcbot,kctop,ldcum,  &
                    ztmst,paph,puen,pven,pmfu,pmfd,puu,pud,pvu,pvd,ptenu, &
                    ptenv,jin)
!---------------------------------------------------------------------------
  USE shr_kind_mod, only: r8 => shr_kind_r8
  USE mo_constants, ONLY: g ! gravity acceleration
!---------------------------------------------------------------------------
    implicit none
    integer  klon,klev,ktopm2,nxj,jin
    integer  ktype(klon), kcbot(klon), kctop(klon)
    logical  ldcum(klon)
    real     ztmst
    real     paph(klon,klev+1)
    real     puen(klon,klev),  pven(klon,klev),&
             pmfu(klon,klev),  pmfd(klon,klev),&
             puu(klon,klev),   pud(klon,klev),&
             pvu(klon,klev),   pvd(klon,klev)
    real     ptenu(klon,klev), ptenv(klon,klev)

!local variables
    real     zuen(klon,klev) , zven(klon,klev) , zmfuu(klon,klev), &
             zmfdu(klon,klev), zmfuv(klon,klev), zmfdv(klon,klev)

    integer  ik , ikb , jk , jl
    real     zzp, zdtdt
   
    real     zdudt(klon,klev), zdvdt(klon,klev), zdp(klon,klev)
!
!xb110>
    zuen = 0.
    zven = 0.
    zdp  = 0.
    zmfuu = 0.
    zmfuv = 0.
    zmfdu = 0.
    zmfdv = 0.
    zdudt = 0.
    zdvdt = 0.
!xb110<
    do jk = 1 , klev
       do jl = 1, nxj         
        if ( ldcum(jl) ) then
          zuen(jl,jk) = puen(jl,jk)
          zven(jl,jk) = pven(jl,jk)
          zdp(jl,jk) = g/(paph(jl,jk+1)-paph(jl,jk))
        end if
      end do
    end do
!----------------------------------------------------------------------
!*    1.0          CALCULATE FLUXES AND UPDATE U AND V TENDENCIES
! ----------------------------------------------
    do jk = ktopm2 , klev
      ik = jk - 1
       do jl = 1, nxj         
        if ( ldcum(jl) ) then
          zmfuu(jl,jk) = pmfu(jl,jk)*(puu(jl,jk)-zuen(jl,ik))
          zmfuv(jl,jk) = pmfu(jl,jk)*(pvu(jl,jk)-zven(jl,ik))
          zmfdu(jl,jk) = pmfd(jl,jk)*(pud(jl,jk)-zuen(jl,ik))
          zmfdv(jl,jk) = pmfd(jl,jk)*(pvd(jl,jk)-zven(jl,ik))
        end if
      end do
    end do
    ! linear fluxes below cloud
      do jk = ktopm2 , klev
       do jl = 1, nxj         
          if ( ldcum(jl) .and. jk > kcbot(jl) ) then
            ikb = kcbot(jl)
            zzp = ((paph(jl,klev+1)-paph(jl,jk))/(paph(jl,klev+1)-paph(jl,ikb)))
            if ( ktype(jl) == 3 ) zzp = zzp*zzp
            zmfuu(jl,jk) = zmfuu(jl,ikb)*zzp
            zmfuv(jl,jk) = zmfuv(jl,ikb)*zzp
            zmfdu(jl,jk) = zmfdu(jl,ikb)*zzp
            zmfdv(jl,jk) = zmfdv(jl,ikb)*zzp
          end if
        end do
      end do
!----------------------------------------------------------------------
!*    2.0          COMPUTE TENDENCIES
! ------------------
    do jk = ktopm2 , klev
      if ( jk < klev ) then
        ik = jk + 1
       do jl = 1, nxj         
          if ( ldcum(jl) ) then
            zdudt(jl,jk) = zdp(jl,jk) * &
        &                  (zmfuu(jl,ik)-zmfuu(jl,jk)+zmfdu(jl,ik)-zmfdu(jl,jk))
            zdvdt(jl,jk) = zdp(jl,jk) * &
        &                  (zmfuv(jl,ik)-zmfuv(jl,jk)+zmfdv(jl,ik)-zmfdv(jl,jk))
          end if
        end do
      else
       do jl = 1, nxj         
          if ( ldcum(jl) ) then
            zdudt(jl,jk) = -zdp(jl,jk)*(zmfuu(jl,jk)+zmfdu(jl,jk))
            zdvdt(jl,jk) = -zdp(jl,jk)*(zmfuv(jl,jk)+zmfdv(jl,jk))
          end if
        end do
      end if
    end do
!---------------------------------------------------------------------
!*  3.0        UPDATE TENDENCIES
!   -----------------
    do jk = ktopm2 , klev
       do jl = 1, nxj         
        if ( ldcum(jl) ) then
          ptenu(jl,jk) = ptenu(jl,jk) + zdudt(jl,jk)
          ptenv(jl,jk) = ptenv(jl,jk) + zdvdt(jl,jk)
        end if
      end do
    end do
!----------------------------------------------------------------------
    return
    end subroutine cududv_n
