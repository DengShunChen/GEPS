     subroutine cudtdq_n(nxj,klon,klev,ktopm2,kctop,kdtop,ldcum, &
                     lddraf,ztmst,paph,pgeoh,pgeo,pten,ptenh,pqen,  &
                     pqenh,pqsen,plglac,plude,pmfu,pmfd,pmfus,pmfds, &
                     pmfuq,pmfdq,pmful,pdmfup,pdmfdp,pdpmel,ptent,ptenq,pcte,pxtec)
!-----------------------------------------------------------------------------
  USE shr_kind_mod, only: r8 => shr_kind_r8
  USE mo_constants, ONLY: rhoh2o, &! density of liquid water
                          tmelt,  &! temperature of fusion of ice
                          alv,    &! latent heat for vaporisation
                          g,      &! gravity acceleration
                          alf,    &! latent heat for fusion
                          als,    &! latent heat for sublimation
                          rcpd     ! rcpd=1./cpd
!-----------------------------------------------------------------------------
    implicit none
    integer  klon,klev,ktopm2,nxj
    integer  kctop(klon),  kdtop(klon)
    logical  ldcum(klon),  lddraf(klon)
    real     ztmst,foelhm
    real     paph(klon,klev+1), pgeoh(klon,klev+1)
    real     pgeo(klon,klev),   pten(klon,klev), &
             pqen(klon,klev),   ptenh(klon,klev),&
             pqenh(klon,klev),  pqsen(klon,klev),&
             plglac(klon,klev), plude(klon,klev)
    real     pmfu(klon,klev),   pmfd(klon,klev),&
             pmfus(klon,klev),  pmfds(klon,klev),&
             pmfuq(klon,klev),  pmfdq(klon,klev),&
             pmful(klon,klev),  pdmfup(klon,klev),&
             pdpmel(klon,klev), pdmfdp(klon,klev)
    real     ptent(klon,klev),  ptenq(klon,klev)
    real     pcte(klon,klev),   pxtec(klon,klev)

! local variables
    integer  jk , ik , jl
    real     zalv , zzp
    real     zmfus(klon,klev) , zmfuq(klon,klev) 
    real     zmfds(klon,klev) , zmfdq(klon,klev)
    real     zdtdt(klon,klev) , zdqdt(klon,klev) , zdp(klon,klev)
    !*    1.0          SETUP AND INITIALIZATIONS
    ! -------------------------
    pcte  = 0.
    pxtec = 0.
!xb110>
    zdp   = 0.
    zmfus = 0.
    zmfds = 0.
    zmfuq = 0.
    zmfdq = 0.
    zdtdt = 0.
    zdqdt = 0.
!xb110<

    do jk = 1 , klev
       do jl = 1, nxj         
        if ( ldcum(jl) ) then
          zdp(jl,jk) = g/(paph(jl,jk+1)-paph(jl,jk))
          zmfus(jl,jk) = pmfus(jl,jk)
          zmfds(jl,jk) = pmfds(jl,jk)
          zmfuq(jl,jk) = pmfuq(jl,jk)
          zmfdq(jl,jk) = pmfdq(jl,jk)
        end if
      end do
    end do
    !-----------------------------------------------------------------------
    !*    2.0          COMPUTE TENDENCIES
    ! ------------------
    do jk = ktopm2 , klev
      if ( jk < klev ) then
       do jl = 1, nxj         
          if ( ldcum(jl) ) then
            zalv = foelhm(pten(jl,jk))
         zdtdt(jl,jk) = zdp(jl,jk)*rcpd * &
        (zmfus(jl,jk+1)-zmfus(jl,jk)+zmfds(jl,jk+1) - &
         zmfds(jl,jk)+alf*plglac(jl,jk)-alf*pdpmel(jl,jk) - &
         zalv*(pmful(jl,jk+1)-pmful(jl,jk)-plude(jl,jk)-pdmfup(jl,jk)-pdmfdp(jl,jk)))
            zdqdt(jl,jk) = zdp(jl,jk)*(zmfuq(jl,jk+1) - &
      &        zmfuq(jl,jk)+zmfdq(jl,jk+1)-zmfdq(jl,jk)+pmful(jl,jk+1) - &
      &        pmful(jl,jk)-plude(jl,jk)-pdmfup(jl,jk)-pdmfdp(jl,jk))
          end if
        end do
      else
       do jl = 1, nxj         
          if ( ldcum(jl) ) then
            zalv = foelhm(pten(jl,jk))
            zdtdt(jl,jk) = -zdp(jl,jk)*rcpd * &
     &         (zmfus(jl,jk)+zmfds(jl,jk)+alf*pdpmel(jl,jk) - &
     &          zalv*(pmful(jl,jk)+pdmfup(jl,jk)+pdmfdp(jl,jk)))
            zdqdt(jl,jk) = -zdp(jl,jk)*(zmfuq(jl,jk) + &
     &         zmfdq(jl,jk)+(pmful(jl,jk)+pdmfup(jl,jk)+pdmfdp(jl,jk)))
          end if
        end do
      end if
    end do
  !---------------------------------------------------------------
  !*  3.0          UPDATE TENDENCIES
  !   -----------------
    do jk = ktopm2 , klev
       do jl = 1, nxj         
       if ( ldcum(jl) ) then
         ptent(jl,jk) = ptent(jl,jk) + zdtdt(jl,jk)
         ptenq(jl,jk) = ptenq(jl,jk) + zdqdt(jl,jk)
         pcte(jl,jk)  = zdp(jl,jk)*plude(jl,jk)
         pxtec(jl,jk) = (g/(paph(jl,jk+1)-paph(jl,jk)))*plude(jl,jk)
       end if
     end do
    end do

    return
  end subroutine cudtdq_n
