      subroutine cuentr_n(nxj,klon,klev,kk,kcbot,ldcum,ldwork, &
                    pgeoh,pmfu,pdmfen,pdmfde)
!--------------------------------------------------------------------
  USE mo_constants, ONLY : g,  & ! gravity acceleration
                           rd, & ! gas constant for dry air
                           zrg        ! 1.0/g
!--------------------------------------------------------------------
       implicit none
       integer  klon,klev,kk,nxj
       integer  kcbot(klon)
       logical  ldcum(klon)
       logical  ldwork
       real  pgeoh(klon,klev+1)
       real  pmfu(klon,klev) 
       real  pdmfen(klon) 
       real  pdmfde(klon) 
       logical  llo1
       integer  jl
       real  zdz , zmf
       real  zentr(klon)
    !
    !* 1. CALCULATE ENTRAINMENT AND DETRAINMENT RATES
    ! -------------------------------------------
    if ( ldwork ) then
       do jl = 1, nxj         
        pdmfen(jl) = 0.
        pdmfde(jl) = 0.
        zentr(jl) = 0.
      end do
      !
      !*  1.1 SPECIFY ENTRAINMENT RATES
      !   -------------------------
       do jl = 1, nxj         
        if ( ldcum(jl) ) then
          zdz = (pgeoh(jl,kk)-pgeoh(jl,kk+1))*zrg
          zmf = pmfu(jl,kk+1)*zdz
          llo1 = kk < kcbot(jl)
          if ( llo1 ) then
            pdmfen(jl) = zentr(jl)*zmf
            pdmfde(jl) = 0.75e-4*zmf
          end if
        end if
      end do
    end if
    end subroutine cuentr_n
