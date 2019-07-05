      subroutine cubasmc_n &
     &    (nxj,      klon,     klev,     klevm1,  kk,     pten,&
     &     pqen,     pqsen,    puen,    pven,   pverv,&
     &     pgeo,     pgeoh,    ldcum,   ktype,  klab,  plrain,&
     &     pmfu,     pmfub,    kcbot,   ptu,&
     &     pqu,      plu,      puu,     pvu,    pmfus,&
     &     pmfuq,    pmful,    pdmfup)
!      m.tiedtke         e.c.m.w.f.     12/89
!      c.zhang           iprc           05/2012
!***purpose.
!   --------
!          this routine calculates cloud base values
!          for midlevel convection
!***interface
!   ---------
!          this routine is called from *cuasc*.
!          input are environmental values t,q etc
!          it returns cloudbase values for midlevel convection
!***method.
!   -------
!          s. tiedtke (1989)
!***externals
!   ---------
!          none
! ----------------------------------------------------------------
  USE mo_constants,    ONLY: cpd,     &! specific heat at constantpressure
                             rcpd,    &! rcpd=1./cpd
                             g,       &! gravity acceleration
                             zrg       ! 1.0/g
  USE mo_cumulus_flux, ONLY: cmfcmin, &! minimum massflux value (for safety)
                             cmfcmax, &! maximum massflux value allowed for
                             entrmid, &! entrainment rate for midlevel convect.
                             lmfdudv, &! true if cumulus friction is switched on
                             lmfmid

!-----------------------------------------------------------------
      implicit none
      real     pten(klon,klev),        pqen(klon,klev),&
     &         puen(klon,klev),        pven(klon,klev),&
     &         pqsen(klon,klev),       pverv(klon,klev),&
     &         pgeo(klon,klev),        pgeoh(klon,klev+1)
      real     ptu(klon,klev),         pqu(klon,klev),&
     &         puu(klon,klev),         pvu(klon,klev),&
     &         plu(klon,klev),         pmfu(klon,klev),&
     &         pmfub(klon),   &         
     &         pmfus(klon,klev),       pmfuq(klon,klev),&
     &         pmful(klon,klev),       pdmfup(klon,klev),&
     &         plrain(klon,klev)
      integer  ktype(klon),            kcbot(klon),&
     &         klab(klon,klev)
      logical  ldcum(klon)
! local variabels
      integer  jl,kk,klev,klon,klevp1,klevm1,nxj
      real     zzzmb
!--------------------------------------------------------
!*    1.      calculate entrainment and detrainment rates
! -------------------------------------------------------
       do jl = 1, nxj         
          if(.not.ldcum(jl) .and. klab(jl,kk+1).eq.0) then
            if(lmfmid .and. pqen(jl,kk) .gt. 0.80*pqsen(jl,kk).and. &
              pgeo(jl,kk)*zrg .gt. 5.0e2  .and.  &
     &        pgeo(jl,kk)*zrg .lt. 1.0e4 )  then
            ptu(jl,kk+1)=(cpd*pten(jl,kk)+pgeo(jl,kk)-pgeoh(jl,kk+1))&
     &                          *rcpd
            pqu(jl,kk+1)=pqen(jl,kk)
            plu(jl,kk+1)=0.
            zzzmb=max(cmfcmin,-pverv(jl,kk)*zrg)
            zzzmb=min(zzzmb,cmfcmax)
            pmfub(jl)=zzzmb
            pmfu(jl,kk+1)=pmfub(jl)
            pmfus(jl,kk+1)=pmfub(jl)*(cpd*ptu(jl,kk+1)+pgeoh(jl,kk+1))
            pmfuq(jl,kk+1)=pmfub(jl)*pqu(jl,kk+1)
            pmful(jl,kk+1)=0.
            pdmfup(jl,kk+1)=0.
            kcbot(jl)=kk
            klab(jl,kk+1)=1
            plrain(jl,kk+1)=0.0
            ktype(jl)=3
          end if
         end if
        end do
      return
      end subroutine cubasmc_n
