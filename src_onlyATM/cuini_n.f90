      subroutine cuini_n &
     &    (nxj,     klon,     klev,     klevp1,   klevm1,   pten,&
     &     pqen,     pqsen,   pxen,     puen,     pven,     pverv,&
     &     pgeo,     paph,     pgeoh,    ptenh,    pqenh,&
     &     pqsenh,   pxenh,    klwmin,   ptu,      pqu,      ptd,&
     &     pqd,      puu,      pvu,      pud,      pvd,&
     &     pmfu,     pmfd,     pmfus,    pmfds,    pmfuq,&
     &     pmfdq,    pdmfup,   pdmfdp,   pdpmel,   plu,&
     &     plude,    klab,     jin)
!------------------------------------------------------------------
  USE shr_kind_mod, only: r8 => shr_kind_r8
  USE mo_constants, only: cpd, &! specific heat at constant pressure
                          rcpd  ! rcpd=1./cpd
!------------------------------------------------------------------
      implicit none
!      m.tiedtke         e.c.m.w.f.     12/89
!***purpose
!   -------
!          this routine interpolates large-scale fields of t,q etc.
!          to half levels (i.e. grid for massflux scheme),
!          and initializes values for updrafts and downdrafts
!***interface
!   ---------
!          this routine is called from *cumastr*.
!***method.
!  --------
!          for extrapolation to half levels see tiedtke(1989)
!***externals
!   ---------
!          *cuadjtq* to specify qs at half levels
! ----------------------------------------------------------------
      integer  klon,klev,klevp1,klevm1,nxj
      real     pten(klon,klev),        pqen(klon,klev),&
               pxen(klon,klev),                        &
     &         puen(klon,klev),        pven(klon,klev),&
     &         pqsen(klon,klev),       pverv(klon,klev),&
     &         pgeo(klon,klev),        pgeoh(klon,klevp1),&
     &         paph(klon,klevp1),      ptenh(klon,klev),&
     &         pqenh(klon,klev),       pqsenh(klon,klev), &
               pxenh(klon,klev)
      real     ptu(klon,klev),         pqu(klon,klev),&
     &         ptd(klon,klev),         pqd(klon,klev),&
     &         puu(klon,klev),         pud(klon,klev),&
     &         pvu(klon,klev),         pvd(klon,klev),&
     &         pmfu(klon,klev),        pmfd(klon,klev),&
     &         pmfus(klon,klev),       pmfds(klon,klev),&
     &         pmfuq(klon,klev),       pmfdq(klon,klev),&
     &         pdmfup(klon,klev),      pdmfdp(klon,klev),&
     &         plu(klon,klev),         plude(klon,klev)
      real     zwmax(klon),            zph(klon), &
     &         pdpmel(klon,klev)
      integer  klab(klon,klev),        klwmin(klon)
      logical  loflag(klon)
!  local variables
      integer  jl,jk,jin
      integer  icall,ik
      real zzs, zdp          
!------------------------------------------------------------
!*    1.       specify large scale parameters at half levels
!*             adjust temperature fields if staticly unstable
!*             find level of maximum vertical velocity
! -----------------------------------------------------------
!xb110>
      pgeoh = 0.
      ptenh = 0.
      pqenh = 0.
      pqsenh = 0.
      pxenh = 0.
      zwmax = 0.
      klwmin = 0
!xb110<
      zdp = 0.5         
      do jk=2,klev
       do jl = 1, nxj        
        pgeoh(jl,jk) = pgeo(jl,jk) + (pgeo(jl,jk-1)-pgeo(jl,jk))*zdp 
        ptenh(jl,jk)=(max(cpd*pten(jl,jk-1)+pgeo(jl,jk-1), &
     &             cpd*pten(jl,jk)+pgeo(jl,jk))-pgeoh(jl,jk))*rcpd
        pqenh(jl,jk) = pqen(jl,jk-1)
        pqsenh(jl,jk)= pqsen(jl,jk-1)
        zph(jl)=paph(jl,jk)
        loflag(jl)=.true.
      end do

      if ( jk >= klev-1 .or. jk < 2 ) cycle
      ik=jk
      icall=0
      call cuadjtq_n(nxj,klon,klev,ik,zph,ptenh,pqsenh,loflag,icall)   
      do jl = 1, nxj          
        pxenh(jl,jk) = (pxen(jl,jk)+pxen(jl,jk-1))*zdp              
        pqenh(jl,jk)=min(pqen(jl,jk-1),pqsen(jl,jk-1)) &
     &            +(pqsenh(jl,jk)-pqsen(jl,jk-1))
        pqenh(jl,jk)=max(pqenh(jl,jk),0.)
      end do
      end do

      do jl = 1, nxj             
        ptenh(jl,klev) = (cpd*pten(jl,klev)+pgeo(jl,klev)-pgeoh(jl,klev))*rcpd
        pxenh(jl,klev) = pxen(jl,klev)
        pqenh(jl,klev) = pqen(jl,klev)
        ptenh(jl,1) = pten(jl,1)
        pxenh(jl,1) = pxen(jl,1)             
        pqenh(jl,1) = pqen(jl,1)
        pgeoh(jl,1) = pgeo(jl,1)          
        klwmin(jl) = klev
        zwmax(jl) = 0.
      end do

      do jk=klevm1,2,-1
       do jl = 1, nxj                
        zzs=max(cpd*ptenh(jl,jk)+pgeoh(jl,jk), &
     &        cpd*ptenh(jl,jk+1)+pgeoh(jl,jk+1))
        ptenh(jl,jk)=(zzs-pgeoh(jl,jk))*rcpd
      end do
      end do

      do jk=klev,3,-1
       do jl = 1,nxj                
        if(pverv(jl,jk).lt.zwmax(jl)) then
           zwmax(jl)=pverv(jl,jk)
           klwmin(jl)=jk
        end if
      end do
      end do
!-----------------------------------------------------------
!*    2.0      initialize values for updrafts and downdrafts
!-----------------------------------------------------------
      do jk=1,klev
      ik=jk-1
      if(jk.eq.1) ik=1
        do jl = 1, nxj                
      ptu(jl,jk)=ptenh(jl,jk)
      ptd(jl,jk)=ptenh(jl,jk)
      pqu(jl,jk)=pqenh(jl,jk)
      pqd(jl,jk)=pqenh(jl,jk)
      plu(jl,jk)=0.
      puu(jl,jk)=puen(jl,ik)
      pud(jl,jk)=puen(jl,ik)
      pvu(jl,jk)=pven(jl,ik)
      pvd(jl,jk)=pven(jl,ik)
      klab(jl,jk)=0
      end do
      end do
      return
      end subroutine cuini_n
